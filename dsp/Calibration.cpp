#include "Calibration.h"
#include "VoiceChain.h" // Biquad
#include <algorithm>
#include <cmath>

void Capture::push(const float* const* ch, int nc, size_t n)
{
    if (!active_.load(std::memory_order_relaxed) || nc <= 0)
        return;
    size_t pos = size_.load(std::memory_order_relaxed);
    const size_t end = std::min(want_, pos + n);
    for (size_t i = 0; pos < end; ++i, ++pos) {
        float s = 0.0f;
        for (int c = 0; c < nc; ++c)
            if (ch[c]) s += ch[c][i];
        buf_[pos] = s / (float)nc;
    }
    size_.store(pos, std::memory_order_release);
}

namespace {
// Profil de reference d'une voix "neutre", bandes en dB relatives a la bande 1 kHz (filtres passe-bande ci-dessous).
// Depart LTASS (Byrne 1994 : +3 / -9 / -14), recale le 2026-10-10 sur les calibrations reelles : milieu entre
// le Vocaster (juge trop brillant) et le micro sans fil (juge trop sombre) -> -0.5 / -6.5 / -11.75,
// puis decale de 1 dB vers le brillant (micro sans fil encore juge un peu etouffe apres correction).
// ponytail: un seul locuteur ; re-centrer avec d'autres voix si les corrections derivent
constexpr float kRefLow = -0.5f, kRefPresence = -5.5f, kRefSibilance = -10.5f;
constexpr float kCorrStrength = 0.75f, kCorrMaxDb = 6.0f; // on corrige 75 % de l'ecart, +-6 dB max
// Niveau d'entree vise apres Gain. -20 poussait la voix 12 dB dans le compresseur (trop compresse a l'ecoute) ;
// -30 = niveau auquel le son de reference a ete valide (gain 0, micros a -30..-39 LUFS) ; le leveler fait le reste.
constexpr float kTargetInputLufs = -30.0f;
// De-reverb active au-dela de 0.3 s estime (l'estimation sous-evalue : piece a 0.8 s mesuree ~0.45 s) ; 0.6 s = tres reverberant
constexpr float kRt60On = 0.3f, kRt60Red = 0.6f;

std::vector<float> filt(const std::vector<float>& x, const Biquad& b1, const Biquad* b2 = nullptr)
{
    std::vector<float> y(x.size());
    float z1[2] = {}, z2[2] = {};
    auto run = [](const Biquad& b, float* z, float v) {
        const float o = b.b0 * v + z[0];
        z[0] = b.b1 * v - b.a1 * o + z[1];
        z[1] = b.b2 * v - b.a2 * o;
        return o;
    };
    for (size_t i = 0; i < x.size(); ++i) {
        float v = run(b1, z1, x[i]);
        y[i] = b2 ? run(*b2, z2, v) : v;
    }
    return y;
}

// Energie moyenne par bloc (taille 'len', pas 'hop')
std::vector<double> blockMs(const std::vector<float>& x, size_t len, size_t hop)
{
    std::vector<double> out;
    for (size_t s = 0; s + len <= x.size(); s += hop) {
        double e = 0;
        for (size_t i = s; i < s + len; ++i) e += (double)x[i] * x[i];
        out.push_back(e / len);
    }
    return out;
}

float lufs(double ms) { return (float)(-0.691 + 10.0 * std::log10(ms + 1e-20)); }
float db(double pow) { return (float)(10.0 * std::log10(pow + 1e-20)); }

float percentile(std::vector<float> v, float p)
{
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    return v[(size_t)std::clamp(p * (v.size() - 1), 0.0f, (float)(v.size() - 1))];
}

int macro(float v) { return (int)std::lround(std::clamp(v, 0.0f, 100.0f)); }

// RT60 aveugle : une piece ne s'eteint pas plus vite que son RT60, la parole si. Les chutes d'energie les plus
// rapides (fenetres de 80 ms) donnent donc la decroissance de la piece ; sur voix seche elles sont bien plus raides.
// ponytail: heuristique (pas le modele de Lollmann) ; percentile et fenetres a recaler sur de vraies pieces
float estimateRt60(const std::vector<float>& x, double sr, float floorDb)
{
    constexpr size_t W = 8; // 8 blocs de 10 ms
    const size_t b10 = (size_t)(0.01 * sr);
    const auto ms = blockMs(x, b10, b10);
    std::vector<float> e;
    for (size_t i = 2; i < ms.size(); ++i) e.push_back(db((ms[i] + ms[i - 1] + ms[i - 2]) / 3)); // lisse sur 30 ms
    std::vector<float> rts;
    for (size_t i = 0; i + W <= e.size(); ++i) {
        if (e[i + W - 1] < floorDb + 10.0f || e[i] - e[i + W - 1] < 3.0f) continue; // au-dessus du bruit, qui chute
        double num = 0, den = 0, me = 0;
        for (size_t k = 0; k < W; ++k) me += e[i + k] / W;
        for (size_t k = 0; k < W; ++k) { num += (k - (W - 1) / 2.0) * (e[i + k] - me); den += (k - (W - 1) / 2.0) * (k - (W - 1) / 2.0); }
        rts.push_back((float)(-60.0 / (num / den / 0.01)));
    }
    return rts.size() < 20 ? 0.0f : percentile(rts, 0.1f);
}
} // namespace

CalibrationReport analyzeCalibration(const std::vector<float>& noise, const std::vector<float>& speech,
                                     double sr)
{
    CalibrationReport r;
    Biquad k1, k2;
    k1.design(Biquad::HighShelf, sr, 1681.0, 0.7071, 4.0);
    k2.design(Biquad::HighPass, sr, 38.0, 0.5, 0);
    const size_t b400 = (size_t)(0.4 * sr), b100 = (size_t)(0.1 * sr), b20 = (size_t)(0.02 * sr);

    // Bruit de fond : sonie moyenne, sans la premiere demi-seconde (clic du debut)
    if (noise.size() > b400 + (size_t)(0.5 * sr)) {
        std::vector<float> n(noise.begin() + (size_t)(0.5 * sr), noise.end());
        double e = 0;
        for (double m : blockMs(filt(n, k1, &k2), b400, b100)) e += m;
        r.noiseLufs = lufs(e / std::max<size_t>(1, blockMs(n, b400, b100).size()));
    }

    // Parole : sonie integree BS.1770 (blocs 400 ms, pas 100 ms, portes -70 absolue et -10 relative)
    const auto kms = blockMs(filt(speech, k1, &k2), b400, b100);
    std::vector<size_t> active;
    double sum = 0;
    size_t cnt = 0;
    for (double m : kms) if (lufs(m) > -70) { sum += m; ++cnt; }
    if (cnt == 0) return r;
    const float rel = lufs(sum / cnt) - 10.0f;
    sum = 0;
    for (size_t i = 0; i < kms.size(); ++i)
        if (lufs(kms[i]) > -70 && lufs(kms[i]) > rel) { active.push_back(i); sum += kms[i]; }
    if (active.size() < 10) return r; // moins d'une seconde de parole utile
    r.speechLufs = lufs(sum / active.size());
    r.snrDb = r.speechLufs - r.noiseLufs;

    for (float s : speech) r.peakDb = std::max(r.peakDb, 20.0f * std::log10(std::fabs(s) + 1e-9f));
    size_t clips = 0;
    for (float s : speech) if (std::fabs(s) > 0.99f) ++clips;
    r.clipPct = 100.0f * clips / std::max<size_t>(1, speech.size());

    // Dynamique : ecart p95 - p10 de la sonie des blocs de parole
    std::vector<float> st;
    for (size_t i : active) st.push_back(lufs(kms[i]));
    r.spreadDb = percentile(st, 0.95f) - percentile(st, 0.10f);

    // Spectre sur les blocs de parole : bandes 150 Hz / 1 kHz / 4 kHz / 7 kHz
    auto bandDb = [&](double hz, double q, std::vector<float>* sig = nullptr) {
        Biquad b; b.design(Biquad::BandPass, sr, hz, q, 0);
        auto y = filt(speech, b);
        const auto ms = blockMs(y, b400, b100);
        double e = 0;
        for (size_t i : active) e += ms[i];
        if (sig) *sig = std::move(y);
        return db(e / active.size());
    };
    std::vector<float> sib;
    const float mid = bandDb(1000, 0.7);
    r.lowDb = bandDb(150, 0.7) - mid - kRefLow;
    r.presenceDb = bandDb(4000, 0.7) - mid - kRefPresence;
    r.sibilanceDb = bandDb(7000, 1.0, &sib) - mid - kRefSibilance;
    // Pics de sifflantes : blocs de 20 ms, p98 au-dessus de la moyenne de la bande
    {
        const auto ms = blockMs(sib, b20, b20);
        std::vector<float> d;
        double e = 0;
        for (double m : ms) e += m;
        const float avg = db(e / std::max<size_t>(1, ms.size()));
        for (double m : ms) d.push_back(db(m) - avg);
        r.sibPeakDb = percentile(d, 0.98f);
    }

    // Reverberation : estimee sur le signal passe-haut (le grondement grave n'est pas de la reverb)
    {
        Biquad hp; hp.design(Biquad::HighPass, sr, 200.0, 0.7071, 0);
        r.rt60 = estimateRt60(filt(speech, hp), sr, r.noiseLufs + 0.691f);
    }

    // Voyants
    r.noise = r.snrDb >= 35 ? Light::Green : r.snrDb >= 25 ? Light::Orange : Light::Red;
    r.level = (r.speechLufs >= -35 && r.speechLufs <= -12) ? Light::Green
            : (r.speechLufs >= -45 && r.speechLufs <= -8) ? Light::Orange : Light::Red;
    r.clipping = (r.clipPct > 0.05f || r.peakDb > -0.1f) ? Light::Red : r.peakDb > -3.0f ? Light::Orange : Light::Green;
    const float sibExcess = std::max(r.sibilanceDb, r.sibPeakDb - 14.0f);
    r.sibilance = sibExcess > 8 ? Light::Red : sibExcess > 4 ? Light::Orange : Light::Green;
    r.reverb = r.rt60 > kRt60Red ? Light::Red : r.rt60 > kRt60On ? Light::Orange : Light::Green;

    // Reglages : on reste pres de 50 (reference validee), ecarts bornes
    r.gainDb = std::clamp(kTargetInputLufs - r.speechLufs, -24.0f, 24.0f);
    // Retour d'ecoute 2026-10-10 : le son a 50 doit rester la signature ; corrections douces et bornees,
    // et jamais moins de debruitage que la reference (Nettoyage < 50 juge trop faible).
    // Le timbre du micro est egalise par la correction cachee ; les macros timbre restent a 50 (signature).
    r.clean = macro(std::clamp(50.0f + (35.0f - r.snrDb) * 1.5f, 50.0f, 85.0f));
    r.warmth = 50;
    r.clarity = 50;
    r.punch = macro(std::clamp(50.0f + 2.0f * (r.spreadDb - 8.0f), 40.0f, 60.0f));
    auto corr = [](float dev) { return std::clamp(-kCorrStrength * dev, -kCorrMaxDb, kCorrMaxDb); };
    r.corrLow = corr(r.lowDb);
    r.corrPresence = corr(r.presenceDb);
    r.corrAir = corr(r.sibilanceDb);
    r.dereverbRt60 = r.rt60 > kRt60On ? std::min(r.rt60, 1.5f) : 0.0f;
    r.valid = true;
    return r;
}
