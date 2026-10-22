#include "VoiceChain.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include "rnnoise.h"

#ifdef GV_WITH_DFN
extern "C" { // API C de libDF (DeepFilterNet, MIT/Apache)
DFState* df_create(const char* path, float atten_lim, const char* log_level);
float df_process_frame(DFState* st, float* input, float* output);
void df_free(DFState* st);
}
#else // sans libDF (WASM, web) : mode Max indisponible, maxAvailable() = false
static DFState* df_create(const char*, float, const char*) { return nullptr; }
static float df_process_frame(DFState*, float*, float*) { return 0; }
static void df_free(DFState*) {}
#endif

namespace {
constexpr double kPi = 3.14159265358979323846;

// ---- Reglages fixes ; le reste derive des macros (applyMacros) ----
// EQ : {type, Hz, Q}. Index 0-1 = soustractif (avant compresseur), 2-4 = additif (apres de-esser)
constexpr struct { Biquad::Type t; double f, q; } kEq[VoiceChain::kNumEq] = {
    {Biquad::Peak,      200.0,  1.0}, // boue
    {Biquad::Peak,      350.0,  1.4}, // carton
    {Biquad::LowShelf,  120.0,  0.7}, // chaleur
    {Biquad::Peak,     4000.0,  1.0}, // presence
    {Biquad::HighShelf, 11000.0, 0.7}, // air
};
constexpr int kNumSubEq = 2;
constexpr float kCompAttMs = 15.0f, kCompRelMs = 100.0f;
constexpr double kDsSplitHz = 5000.0, kDsDetectHz = 6500.0;
constexpr float kDsMaxCutDb = 8.0f;
constexpr float kLvMaxBoostDb = 15.0f, kLvMaxCutDb = 12.0f, kLvDbPerSec = 2.0f, kLvGateLufs = -50.0f;
constexpr float kLimCeilDb = -1.0f; // ponytail: plafond sur crete echantillon, pas true-peak ; ajouter suréchantillonnage x4 si l'ecart gene
constexpr float kLimLookaheadMs = 1.5f, kLimRelMs = 80.0f;

float dbToLin(float db) { return std::pow(10.0f, db / 20.0f); }
float linToDb(float x) { return 20.0f * std::log10(x + 1e-9f); }
float coefMs(float ms, double sr) { return std::exp(-1.0f / (0.001f * ms * (float)sr)); }

// Forme directe transposee II
inline float run(const Biquad& b, float* z, float x)
{
    const float y = b.b0 * x + z[0];
    z[0] = b.b1 * x - b.a1 * y + z[1];
    z[1] = b.b2 * x - b.a2 * y;
    return y;
}
} // namespace

// Formules RBJ (Audio EQ Cookbook)
void Biquad::design(Type t, double sr, double f, double q, double gainDb)
{
    const double A = std::pow(10.0, gainDb / 40.0);
    const double w = 2.0 * kPi * f / sr, cw = std::cos(w), alpha = std::sin(w) / (2.0 * q);
    const double sA = 2.0 * std::sqrt(A) * alpha;
    double B0, B1, B2, A0, A1, A2;
    switch (t) {
    case HighPass:
        B0 = (1 + cw) / 2; B1 = -(1 + cw); B2 = (1 + cw) / 2;
        A0 = 1 + alpha; A1 = -2 * cw; A2 = 1 - alpha;
        break;
    case Peak:
        B0 = 1 + alpha * A; B1 = -2 * cw; B2 = 1 - alpha * A;
        A0 = 1 + alpha / A; A1 = -2 * cw; A2 = 1 - alpha / A;
        break;
    case LowShelf:
        B0 = A * ((A + 1) - (A - 1) * cw + sA); B1 = 2 * A * ((A - 1) - (A + 1) * cw);
        B2 = A * ((A + 1) - (A - 1) * cw - sA);
        A0 = (A + 1) + (A - 1) * cw + sA; A1 = -2 * ((A - 1) + (A + 1) * cw);
        A2 = (A + 1) + (A - 1) * cw - sA;
        break;
    case BandPass: // gain 0 dB au centre
        B0 = alpha; B1 = 0; B2 = -alpha;
        A0 = 1 + alpha; A1 = -2 * cw; A2 = 1 - alpha;
        break;
    default: // HighShelf
        B0 = A * ((A + 1) + (A - 1) * cw + sA); B1 = -2 * A * ((A - 1) + (A + 1) * cw);
        B2 = A * ((A + 1) + (A - 1) * cw - sA);
        A0 = (A + 1) - (A - 1) * cw + sA; A1 = 2 * ((A - 1) - (A + 1) * cw);
        A2 = (A + 1) - (A - 1) * cw - sA;
        break;
    }
    b0 = (float)(B0 / A0); b1 = (float)(B1 / A0); b2 = (float)(B2 / A0);
    a1 = (float)(A1 / A0); a2 = (float)(A2 / A0);
}

void VoiceChain::prepare(double sampleRate, int numChannels, const char* dfnModel)
{
    sampleRate_ = sampleRate > 0 ? sampleRate : 48000.0;
    numChannels_ = std::clamp(numChannels, 1, kMaxChannels);
    const double sr = sampleRate_;

    smoothCoef_ = coefMs(20.0f, sr);
    currentGain_ = targetGain_;

    dsSplit_.design(Biquad::HighPass, sr, kDsSplitHz, 0.7071, 0);
    // passe-bande via pic etroit : on n'utilise que la difference (y - x) = bande isolee
    dsDetect_.design(Biquad::Peak, sr, kDsDetectHz, 1.0, 20.0);

    compAtt_ = coefMs(kCompAttMs, sr);
    compRel_ = coefMs(kCompRelMs, sr);
    dsAtt_ = coefMs(2.0f, sr);
    dsRel_ = coefMs(60.0f, sr);
    limRel_ = coefMs(kLimRelMs, sr);
    abC_ = coefMs(3000.0f, sr);
    abSmooth_ = coefMs(50.0f, sr);
    dry_.assign((size_t)kMaxChannels * kDryLen, 0.0f);
    dryPos_ = 0;
    abWet_ = abDry_ = 0.0f;
    abGain_ = 1.0f;

    // filtre K (BS.1770) approche en RBJ : plateau +4 dB @1681 Hz, passe-haut 38 Hz
    kw1_.design(Biquad::HighShelf, sr, 1681.0, 0.7071, 4.0);
    kw2_.design(Biquad::HighPass, sr, 38.0, 0.5, 0);
    lvMomC_ = coefMs(400.0f, sr);
    lvShortC_ = coefMs(3000.0f, sr);
    lvStep_ = kLvDbPerSec / (float)sr;
    std::memset(kwZ_, 0, sizeof(kwZ_));
    lvGainDb_ = 0.0f;
    applyMacros(pending_);
    lvMom_ = lvShort_ = std::pow(10.0f, (lvTarget_ + 0.691f) / 10.0f); // demarre "a la cible" : gain 0 dB

    lookahead_ = std::max(1, (int)std::lround(kLimLookaheadMs * 0.001 * sr));
    delay_.assign((size_t)kMaxChannels * lookahead_, 0.0f);
    reqGain_.assign((size_t)lookahead_ + 1, 1.0f); // couvre l'echantillon sortant + la fenetre

    std::fill(&hpfZ_[0][0], &hpfZ_[0][0] + sizeof(hpfZ_) / sizeof(float), 0.0f);
    std::fill(&eqZ_[0][0][0], &eqZ_[0][0][0] + sizeof(eqZ_) / sizeof(float), 0.0f);
    std::memset(corrZ_, 0, sizeof(corrZ_));
    std::fill(&dsSplitZ_[0][0], &dsSplitZ_[0][0] + sizeof(dsSplitZ_) / sizeof(float), 0.0f);
    std::fill(&dsDetectZ_[0][0], &dsDetectZ_[0][0] + sizeof(dsDetectZ_) / sizeof(float), 0.0f);
    for (int c = 0; c < kMaxChannels; ++c)
        dn_[c].reset(c < numChannels_ ? rnnoise_create(nullptr) : nullptr); // ponytail: 1 RNNoise par canal ; mono si le CPU gene
    std::memset(dnIn_, 0, sizeof(dnIn_));
    std::memset(dnOut_, 0, sizeof(dnOut_));
    std::memset(dnHist_, 0, sizeof(dnHist_));
    dnPos_ = dnHistIdx_ = 0;
    for (auto& d : dr_) d.prepare();

    // DeepFilterNet : libDF fait abort() si le modele manque -> on verifie avant.
    // ponytail: RNNoise et DFN supposent 48 kHz ; a d'autres frequences ils tournent quand meme, mal
    const bool dfnOk = dfnModel && std::ifstream(dfnModel).good();
    for (int c = 0; c < kMaxChannels; ++c) // ponytail: 1 DFN par canal (~2x CPU en stereo) ; downmix mono si ca rame
        dfn_[c].reset(dfnOk && c < numChannels_ ? df_create(dfnModel, 100.0f, nullptr) : nullptr);
    compEnvDb_ = -120.0f;
    dsEnv_ = 0.0f;
    pos_ = 0;
    limGain_ = 1.0f;
}

void VoiceChain::DnFree::operator()(DenoiseState* s) const { rnnoise_destroy(s); }
void VoiceChain::DfFree::operator()(DFState* s) const { df_free(s); }

void VoiceChain::denoiseFrame(int nc)
{
    const int h = dnHistIdx_;
    dnHistIdx_ = (dnHistIdx_ + 1) % kHist;
    for (int c = 0; c < nc; ++c) {
        std::memcpy(dnHist_[c][h], dnIn_[c], sizeof(dnIn_[c]));
        if (denoise_ && engine_ == Max && dfn_[c]) {
            // DFN : echelle -1..1, retard interne 3 trames (mesure : 1440 echantillons)
            // ponytail: libDF alloue dans le thread audio (tract) ; inevitable sans reecrire l'inference
            df_process_frame(dfn_[c].get(), dnIn_[c], dnOut_[c]);
        } else if (denoise_ && engine_ == Fast) {
            // RNNoise : echelle int16, retard interne 1 trame
            float tmp[kDnFrame];
            for (int k = 0; k < kDnFrame; ++k) tmp[k] = dnIn_[c][k] * 32768.0f;
            rnnoise_process_frame(dn_[c].get(), dnOut_[c], tmp);
            for (float& v : dnOut_[c]) v *= 1.0f / 32768.0f;
        } else {
            // hors service (ou modele DFN absent) : meme retard que le moteur choisi
            std::memcpy(dnOut_[c], dnHist_[c][(h - engineDelayFrames() + kHist) % kHist], sizeof(dnOut_[c]));
            continue;
        }
        if (dnWet_ < 1.0f) { // signal sec aligne = meme retard que le moteur
            const float* dry = dnHist_[c][(h - engineDelayFrames() + kHist) % kHist];
            for (int k = 0; k < kDnFrame; ++k)
                dnOut_[c][k] = dnWet_ * dnOut_[c][k] + (1.0f - dnWet_) * dry[k];
        }
    }
    for (int c = 0; drActive_ && c < nc; ++c) dr_[c].process(dnOut_[c]);
}

// Macros -> parametres. 0.5 partout = reference (HPF 80 Hz, boue -2.5, carton -1.5, chaleur +1.5,
// presence +2, air +2, comp 3:1 @-20 dB +4 dB, de-esser -28 dB, debruitage plein).
void VoiceChain::applyMacros(const Macros& m)
{
    const double sr = sampleRate_;
    const float up = std::max(0.0f, m.clean - 0.5f) * 2.0f; // 0..1 sur la moitie haute
    // Nettoyage : bas = moins de debruitage (remelange sec), haut = HPF plus haut + EQ soustractif plus marque
    dnWet_ = std::min(1.0f, 0.6f + 0.8f * m.clean);
    hpf_.design(Biquad::HighPass, sr, 80.0 + 40.0 * up, 0.7071, 0);
    eq_[0].design(kEq[0].t, sr, kEq[0].f, kEq[0].q, -2.5 - 2.5 * up);
    eq_[1].design(kEq[1].t, sr, kEq[1].f, kEq[1].q, -1.5 - 1.5 * up);
    // Chaleur : plateau grave -1.5 .. +4.5 dB
    eq_[2].design(kEq[2].t, sr, kEq[2].f, kEq[2].q, -1.5 + 6.0 * m.warmth);
    // Clarte : presence et air 0 .. +4 dB ; le de-esser se resserre quand on ajoute des aigus
    eq_[3].design(kEq[3].t, sr, kEq[3].f, kEq[3].q, 4.0 * m.clarity);
    eq_[4].design(kEq[4].t, sr, kEq[4].f, kEq[4].q, 4.0 * m.clarity);
    dsThreshDb_ = -28.0f - 6.0f * (m.clarity - 0.5f) * 2.0f;
    // Punch : ratio 1.5..4.5, seuil -14..-26 dB, compensation automatique (4 dB a la reference)
    compRatio_ = 1.5f + 3.0f * m.punch;
    compThreshDb_ = -14.0f - 12.0f * m.punch;
    compMakeup_ = dbToLin(0.3f * -compThreshDb_ * (1.0f - 1.0f / compRatio_));
    lvTarget_ = targetLufs(m.dest);
    // Correction du micro, avant le compresseur et le de-esser (ils travaillent sur un micro deja egalise)
    corr_[0].design(Biquad::LowShelf, sr, 150.0, 0.7, m.corrLow);
    corr_[1].design(Biquad::Peak, sr, 4000.0, 0.7, m.corrPresence);
    corr_[2].design(Biquad::HighShelf, sr, 7000.0, 0.7, m.corrAir);
    const bool dr = stages_.dereverb && m.rt60 > 0.05f;
    for (auto& d : dr_) {
        if (dr && !drActive_) d.prepare(); // repart d'un etat propre ; pas d'allocation (tailles inchangees)
        d.setRt60(dr ? m.rt60 : 0.0f);
    }
    // ponytail: la latence saute de 10 ms sans fondu (petit trou ou doublon audible) ; reglage rare, pas en direct
    drActive_ = dr;
}

void VoiceChain::setGainDb(float db)
{
    targetGain_ = dbToLin(db);
}

void VoiceChain::process(float* const* channels, int numChannels, size_t numFrames)
{
    if (numChannels <= 0 || delay_.empty())
        return;
    const int nc = std::min(numChannels, numChannels_);
    const float ceil = dbToLin(kLimCeilDb);
    if (dirty_.exchange(false))
        applyMacros(pending_);
    const Stages on = stages_; // ponytail: copie sans verrou comme pending_
    const int L = lookahead_, R = L + 1;
    float x[kMaxChannels] = {};

    for (size_t i = 0; i < numFrames; ++i) {
        currentGain_ = targetGain_ + smoothCoef_ * (currentGain_ - targetGain_);

        // 1-2. Trim + passe-haut, 3. debruitage (FIFO de trames) ; brut garde pour Avant/Apres
        for (int c = 0; c < nc; ++c) {
            const float raw = channels[c] ? channels[c][i] : 0.0f;
            dry_[(size_t)c * kDryLen + dryPos_] = raw;
            const float s = raw * currentGain_;
            x[c] = dnOut_[c][dnPos_];
            dnIn_[c][dnPos_] = on.hpf ? run(hpf_, hpfZ_[c], s) : s;
        }
        if (++dnPos_ == kDnFrame) {
            denoiseFrame(nc);
            dnPos_ = 0;
        }

        // 5. EQ soustractif
        float peak = 0.0f;
        for (int c = 0; c < nc; ++c) {
            float s = x[c];
            for (int k = 0; on.eq && k < 3; ++k)
                s = run(corr_[k], corrZ_[k][c], s);
            for (int k = 0; on.eq && k < kNumSubEq; ++k)
                s = run(eq_[k], eqZ_[k][c], s);
            x[c] = s;
            peak = std::max(peak, std::fabs(s));
        }

        // 6. Compresseur (detecteur crete en dB, canaux lies)
        const float lvl = linToDb(peak);
        compEnvDb_ = lvl + (lvl > compEnvDb_ ? compAtt_ : compRel_) * (compEnvDb_ - lvl);
        const float over = compEnvDb_ - compThreshDb_;
        const float compGain = !on.compressor ? 1.0f
                             : compMakeup_ * (over > 0 ? dbToLin(-over * (1.0f - 1.0f / compRatio_)) : 1.0f);

        // 7. De-esser : detection 6.5 kHz, attenuation de la bande > 5 kHz seulement
        float hi[kMaxChannels] = {}, sib = 0.0f;
        for (int c = 0; c < nc; ++c) {
            x[c] *= compGain;
            hi[c] = run(dsSplit_, dsSplitZ_[c], x[c]);
            sib = std::max(sib, std::fabs(run(dsDetect_, dsDetectZ_[c], x[c]) - x[c]) * 0.1f); // /10 = -20 dB du pic
        }
        dsEnv_ = sib + (sib > dsEnv_ ? dsAtt_ : dsRel_) * (dsEnv_ - sib);
        const float dsOver = linToDb(dsEnv_) - dsThreshDb_;
        const float dsGain = on.deesser && dsOver > 0 ? dbToLin(-std::min(kDsMaxCutDb, dsOver * 0.75f)) : 1.0f;

        // 8. EQ additif + mesure de sonie (filtre K)
        float ms = 0.0f;
        for (int c = 0; c < nc; ++c) {
            float s = x[c] + (dsGain - 1.0f) * hi[c]; // bas + dsGain*haut, reconstruction exacte si dsGain=1
            for (int k = kNumSubEq; on.eq && k < kNumEq; ++k)
                s = run(eq_[k], eqZ_[k][c], s);
            x[c] = s;
            const float k = run(kw2_, kwZ_[c][1], run(kw1_, kwZ_[c][0], s));
            ms += k * k;
        }
        if (nc > 2) ms *= 2.0f / nc; // ponytail: canaux > 2 moyennes ; BS.1770 ponderera les surround si un jour utile

        // 9. Leveler : court terme 3 s mis a jour seulement quand on parle, gain qui glisse a 2 dB/s
        lvMom_ = ms + lvMomC_ * (lvMom_ - ms);
        const float mom = -0.691f + 10.0f * std::log10(lvMom_ + 1e-12f);
        if (mom > kLvGateLufs) { // ponytail: porte absolue seule (le debruitage rend les pauses tres basses) ; porte relative BS.1770 si besoin
            lvShort_ = ms + lvShortC_ * (lvShort_ - ms);
            const float st = -0.691f + 10.0f * std::log10(lvShort_ + 1e-12f);
            const float want = std::clamp(lvTarget_ - st, -kLvMaxCutDb, kLvMaxBoostDb);
            lvGainDb_ += std::clamp(want - lvGainDb_, -lvStep_, lvStep_);
        }
        const float lvGain = on.leveler ? dbToLin(lvGainDb_) : 1.0f;
        peak = 0.0f;
        for (int c = 0; c < nc; ++c) {
            x[c] *= lvGain;
            peak = std::max(peak, std::fabs(x[c]));
        }

        // 10. Limiteur a anticipation : gain = min du gain requis sur la fenetre, relache lente
        reqGain_[pos_ % R] = on.limiter && peak > ceil ? ceil / peak : 1.0f;
        const float target = *std::min_element(reqGain_.begin(), reqGain_.end()); // ponytail: O(L) par echantillon, L=72 ; deque monotone si L grandit
        limGain_ = target < limGain_ ? target : target + limRel_ * (limGain_ - target);

        // Avant/Apres : brut aligne sur la latence, volume egalise sur 3 s quand on parle
        const int dryRead = (dryPos_ - latencyFor(engine_, drActive_) + kDryLen) % kDryLen;
        float wetE = 0.0f, dryE = 0.0f;
        const int dpos = pos_ % L;
        for (int c = 0; c < nc; ++c) {
            float* d = &delay_[(size_t)c * L + dpos];
            float out = *d * limGain_;
            *d = x[c];
            const float dry = dry_[(size_t)c * kDryLen + dryRead];
            wetE += out * out;
            dryE += dry * dry;
            if (bypass_)
                out = std::clamp(dry * abGain_, -1.0f, 1.0f);
            if (channels[c])
                channels[c][i] = out;
        }
        if (wetE > 1e-7f * nc) { // ~-70 dBFS : on n'egalise que sur la voix
            abWet_ = wetE + abC_ * (abWet_ - wetE);
            abDry_ = dryE + abC_ * (abDry_ - dryE);
        }
        const float want = std::clamp(std::sqrt(abWet_ / (abDry_ + 1e-12f)), 0.03f, 30.0f); // +-30 dB
        abGain_ = want + abSmooth_ * (abGain_ - want);
        dryPos_ = (dryPos_ + 1) % kDryLen;
        pos_ = (pos_ + 1) % (L * R);
    }
}
