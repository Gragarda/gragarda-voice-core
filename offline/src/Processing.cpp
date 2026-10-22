#include "gv/Processing.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include "gv/Dsp.h"
#include "Dereverb.h"
#include "rnnoise.h"

namespace gv {

namespace {
// ponytail: lecteur JSON minimal (objet plat cle -> nombre/chaine) ; une vraie lib si les presets s'imbriquent
bool findValue(const std::string& j, const char* key, size_t& at)
{
    const std::string k = std::string("\"") + key + "\"";
    size_t p = j.find(k);
    if (p == std::string::npos) return false;
    p = j.find(':', p + k.size());
    if (p == std::string::npos) return false;
    at = j.find_first_not_of(" \t\r\n", p + 1);
    return at != std::string::npos;
}
void num(const std::string& j, const char* key, float& v)
{
    size_t at;
    if (findValue(j, key, at)) v = std::strtof(j.c_str() + at, nullptr);
}

float pct(std::vector<float> v, float p)
{
    if (v.empty()) return -120.0f;
    const size_t k = std::min(v.size() - 1, (size_t)(p * (v.size() - 1) + 0.5f));
    std::nth_element(v.begin(), v.begin() + k, v.end());
    return v[k];
}

constexpr int kFrame = 480; // RNNoise : 10 ms a 48 kHz, retard interne 1 trame

void denoise(AudioBuffer& b, float wet)
{
    const size_t n = b.frames();
    DenoiseState* st = rnnoise_create(nullptr);
    float in[kFrame], out[kFrame];
    for (auto& c : b.ch) {
        rnnoise_destroy(st); // etat neuf par canal
        st = rnnoise_create(nullptr);
        std::vector<float> y(n, 0.0f);
        for (size_t f = 0; f * kFrame < n + kFrame; ++f) {
            for (int k = 0; k < kFrame; ++k) {
                const size_t i = f * kFrame + k;
                in[k] = i < n ? c[i] * 32768.0f : 0.0f;
            }
            rnnoise_process_frame(st, out, in);
            if (f == 0) continue;
            for (int k = 0; k < kFrame; ++k) {
                const size_t i = (f - 1) * kFrame + k;
                if (i < n) y[i] = out[k] / 32768.0f;
            }
        }
        for (size_t i = 0; i < n; ++i) c[i] = wet * y[i] + (1.0f - wet) * c[i];
    }
    rnnoise_destroy(st);
}

// Gain lisse applique a tous les canaux
void applyGain(AudioBuffer& b, const std::vector<float>& g)
{
    for (auto& c : b.ch)
        for (size_t i = 0; i < c.size(); ++i) c[i] *= g[i];
}

std::vector<float> maxAbs(const AudioBuffer& b)
{
    std::vector<float> m(b.frames(), 0.0f);
    for (const auto& c : b.ch)
        for (size_t i = 0; i < m.size(); ++i) m[i] = std::max(m[i], std::fabs(c[i]));
    return m;
}

// De-reverb (dsp/Dereverb) sur tout le fichier, retard d'une trame compense
void dereverb(AudioBuffer& b, float rt60)
{
    const size_t n = b.frames(), F = Dereverb::kFrame;
    std::vector<float> fr(F);
    for (auto& c : b.ch) {
        Dereverb d;
        d.prepare();
        d.setRt60(rt60);
        std::vector<float> y(n);
        for (size_t f = 0; f * F < n + F; ++f) {
            for (size_t k = 0; k < F; ++k) fr[k] = f * F + k < n ? c[f * F + k] : 0.0f;
            d.process(fr.data());
            for (size_t k = 0; f > 0 && k < F; ++k)
                if ((f - 1) * F + k < n) y[(f - 1) * F + k] = fr[k];
        }
        c.swap(y);
    }
}

void deess(AudioBuffer& b, float thrDb)
{
    const double sr = b.sampleRate;
    if (sr < 20000) return;
    Biquad hp;
    hp.design(Biquad::HighPass, sr, 5000, 0.7071);
    std::vector<std::vector<float>> hi;
    for (const auto& c : b.ch) { hi.push_back(c); hp.process(hi.back()); }
    const float att = coefMs(1, sr), rel = coefMs(50, sr);
    float eh = 0, ef = 0;
    for (size_t i = 0; i < b.frames(); ++i) {
        float ph = 0, pf = 0;
        for (size_t c = 0; c < hi.size(); ++c) { ph += hi[c][i] * hi[c][i]; pf += b.ch[c][i] * b.ch[c][i]; }
        eh = ph + (ph > eh ? att : rel) * (eh - ph);
        ef = pf + (pf > ef ? att : rel) * (ef - pf);
        const float over = 10.0f * std::log10((eh + 1e-14f) / (ef + 1e-14f)) - thrDb;
        if (over <= 0) continue;
        const float g = dbToLin(-std::min(8.0f, 0.75f * over));
        for (size_t c = 0; c < hi.size(); ++c) b.ch[c][i] += (g - 1.0f) * hi[c][i]; // bas + g * haut
    }
}

// Enveloppe crete en dB (attaque 10 ms, relachement 120 ms), canaux lies
std::vector<float> peakEnvDb(const AudioBuffer& b)
{
    const float att = coefMs(10, b.sampleRate), rel = coefMs(120, b.sampleRate);
    auto e = maxAbs(b);
    float env = -120.0f;
    for (float& v : e) {
        const float l = linToDb(v);
        env = l + (l > env ? att : rel) * (env - l);
        v = env;
    }
    return e;
}

void compress(AudioBuffer& b, const std::vector<float>& env, float thrDb, float ratio)
{
    std::vector<float> g(env.size());
    for (size_t i = 0; i < g.size(); ++i)
        g[i] = env[i] > thrDb ? dbToLin(-(env[i] - thrDb) * (1.0f - 1.0f / ratio)) : 1.0f;
    applyGain(b, g);
}
} // namespace

Preset parsePreset(const std::string& j)
{
    Preset p;
    size_t at;
    if (findValue(j, "name", at) && j[at] == '"') p.name = j.substr(at + 1, j.find('"', at + 1) - at - 1);
    num(j, "hpfHz", p.hpfHz); num(j, "denoiseMin", p.denoiseMin); num(j, "levelMaxDb", p.levelMaxDb);
    num(j, "lowShelfHz", p.lowShelfHz); num(j, "lowShelfDb", p.lowShelfDb);
    num(j, "peak1Hz", p.peak1Hz); num(j, "peak1Db", p.peak1Db); num(j, "peak1Q", p.peak1Q);
    num(j, "peak2Hz", p.peak2Hz); num(j, "peak2Db", p.peak2Db); num(j, "peak2Q", p.peak2Q);
    num(j, "highShelfHz", p.highShelfHz); num(j, "highShelfDb", p.highShelfDb);
    num(j, "deessOffsetDb", p.deessOffsetDb); num(j, "targetRangeDb", p.targetRangeDb);
    num(j, "targetLufs", p.targetLufs); num(j, "ceilingDbtp", p.ceilingDbtp); num(j, "dereverbRt60", p.dereverbRt60);
    auto flag = [&](const char* k, bool& v) { float f = v ? 1.0f : 0.0f; num(j, k, f); v = f != 0.0f; };
    flag("pops", p.use.pops); flag("hpf", p.use.hpf); flag("notch", p.use.notch); flag("denoise", p.use.denoise);
    flag("dereverb", p.use.dereverb); flag("eq", p.use.eq); flag("deesser", p.use.deesser); flag("leveler", p.use.leveler);
    flag("compressor", p.use.compressor); flag("loudness", p.use.loudness); flag("limiter", p.use.limiter);
    return p;
}

void limitTruePeak(AudioBuffer& b, float ceiling, float lookaheadMs, float releaseMs)
{
    const size_t n = b.frames(), W = (size_t)std::lround(lookaheadMs * 0.001 * b.sampleRate);
    if (!n) return;
    auto req = truePeakEnvelope(b);
    for (float& v : req) v = v > ceiling ? ceiling / v : 1.0f;
    // h[i] = min(req[i..i+W]) ; puis moyenne de h sur [i-W+1, i] : chaque terme <= req[i], donc le gain aussi
    std::vector<float> h(n);
    std::deque<size_t> dq;
    for (size_t k = n; k-- > 0;) {
        while (!dq.empty() && req[dq.back()] >= req[k]) dq.pop_back();
        dq.push_back(k);
        while (dq.front() > k + W) dq.pop_front();
        h[k] = req[dq.front()];
    }
    const float rel = coefMs(releaseMs, b.sampleRate);
    double acc = 0;
    float g = 1.0f;
    for (size_t i = 0; i < n; ++i) {
        acc += h[i];
        if (i >= W) acc -= h[i - W];
        const float avg = (float)(acc / std::min(i + 1, W));
        g = std::min(avg, 1.0f - (1.0f - g) * rel); // remonte lentement, jamais au-dessus du requis
        req[i] = g;
    }
    applyGain(b, req);
}

// Nivellement lent multi-personnes : chaque passage de parole (fenetre 1 s autour, parole seulement) est ramene
// vers le niveau median de la video, +-maxDb. Une personne loin du micro et une proche finissent au meme volume.
// Les pauses gardent le gain de la parole voisine (pas de remontee du bruit dans les silences).
// speechDb : seuil au-dessus duquel un bloc de 100 ms compte comme de la parole.
static float levelSpeakers(AudioBuffer& b, float speechDb, float maxDb)
{
    const size_t B = (size_t)(0.1 * b.sampleRate), n = b.frames() / B; // blocs de 100 ms
    if (n < 20) return 0;
    std::vector<float> lvl(n);
    for (size_t k = 0; k < n; ++k) {
        double e = 0;
        for (const auto& c : b.ch)
            for (size_t i = k * B; i < (k + 1) * B; ++i) e += (double)c[i] * c[i];
        lvl[k] = (float)(10 * std::log10(e / (B * b.ch.size()) + 1e-12));
    }
    std::vector<float> sp;
    for (float l : lvl)
        if (l > speechDb) sp.push_back(l);
    if (sp.size() < 10) return 0;
    const float target = pct(sp, 0.5f);
    // niveau local de la parole : moyenne en puissance des blocs de parole a +-5 blocs
    std::vector<float> gDb(n, NAN);
    for (size_t k = 0; k < n; ++k) {
        if (lvl[k] <= speechDb) continue;
        double e = 0;
        int cnt = 0;
        for (size_t i = k >= 5 ? k - 5 : 0; i < std::min(n, k + 6); ++i)
            if (lvl[i] > speechDb) { e += std::pow(10.0, lvl[i] / 10); ++cnt; }
        gDb[k] = std::clamp(target - (float)(10 * std::log10(e / cnt)), -maxDb, maxDb);
    }
    // pauses : interpolation entre les gains de parole voisins
    size_t last = n;
    for (size_t k = 0; k < n; ++k) {
        if (std::isnan(gDb[k])) continue;
        if (last == n) std::fill(gDb.begin(), gDb.begin() + k, gDb[k]);
        else for (size_t i = last + 1; i < k; ++i) gDb[i] = gDb[last] + (gDb[k] - gDb[last]) * (i - last) / (k - last);
        last = k;
    }
    if (last == n) return 0;
    std::fill(gDb.begin() + last + 1, gDb.end(), gDb[last]);
    // lissage (moyenne glissante 500 ms) puis interpolation lineaire par echantillon
    std::vector<float> s(n);
    float maxAbs = 0;
    for (size_t k = 0; k < n; ++k) {
        float acc = 0;
        int cnt = 0;
        for (size_t i = k >= 2 ? k - 2 : 0; i < std::min(n, k + 3); ++i) { acc += gDb[i]; ++cnt; }
        s[k] = acc / cnt;
        maxAbs = std::max(maxAbs, std::fabs(s[k]));
    }
    std::vector<float> g(b.frames());
    for (size_t i = 0; i < g.size(); ++i) {
        const float pos = (float)i / B - 0.5f; // centre des blocs
        const size_t k = (size_t)std::clamp(pos, 0.0f, (float)(n - 1));
        const float f = std::clamp(pos - k, 0.0f, 1.0f);
        g[i] = dbToLin(k + 1 < n ? s[k] + (s[k + 1] - s[k]) * f : s[n - 1]);
    }
    applyGain(b, g);
    return maxAbs;
}

ProcessInfo processChain(AudioBuffer& b, const Preset& p, const AnalysisResult& a, const Denoiser& externalDenoise)
{
    ProcessInfo info;
    const double sr = b.sampleRate;
    if (!b.frames()) return info;

    // 1. Passe-haut 18 dB/oct : Butterworth ordre 3 = 1er ordre + biquad Q=1
    Biquad hp;
    hp.design(Biquad::HighPass, sr, p.hpfHz, 1.0);
    for (auto& c : b.ch) if (p.use.hpf) { highPass1(c, sr, p.hpfHz); hp.process(c); }

    // 1b. Raies tonales stables : coupe-bande etroit (Q 30, ~+-12 Hz a 720 Hz), 2 passes = ~ -40 dB au centre
    for (float f : p.use.notch ? a.toneHz : std::vector<float>{}) {
        Biquad nf;
        nf.design(Biquad::Notch, sr, f, 30.0);
        for (auto& c : b.ch) { nf.process(c); nf.process(c); }
    }

    // 2. Debruitage. Externe (web : DeepFilterNet 3) : toujours, il separe la parole de tout le reste, y compris
    // une respiration pres du micro pendant que d'autres parlent. Sinon RNNoise, dose selon le bruit :
    // rien sous -70 dBFS (evite la "voix sous l'eau"), plein a partir de -45.
    if (!p.use.denoise) {
        info.denoiser = "off";
    } else if (externalDenoise) {
        for (auto& c : b.ch) externalDenoise(c);
        info.denoiser = "deepfilternet";
        info.denoiseWet = 1;
    } else {
        info.denoiseWet = std::max(p.denoiseMin, std::clamp((a.noiseDb + 70.0f) / 25.0f, 0.0f, 1.0f));
        if (sr != 48000) info.denoiseWet = 0; // ponytail: RNNoise entraine a 48 kHz ; le web extrait toujours en 48 kHz
        if (info.denoiseWet > 0.01f) denoise(b, info.denoiseWet);
    }

    // 3. De-reverb : RT60 mesure sur la voix debruitee (le bruit masquerait les queues), mono passe-haut 200 Hz
    // (le grondement grave n'est pas de la reverb). Fond = celui de la prise brute : prudent, il est plus haut.
    {
        std::vector<float> m(b.frames(), 0.0f);
        for (const auto& c : b.ch)
            for (size_t i = 0; i < m.size(); ++i) m[i] += c[i] / b.channels();
        Biquad hp200;
        hp200.design(Biquad::HighPass, sr, 200.0, 0.7071);
        hp200.process(m);
        info.rt60 = estimateRt60(m, sr, a.noiseDb);
        info.dereverb = info.rt60 > kRt60On || p.dereverbRt60 > 0;
        info.dereverb = info.dereverb && p.use.dereverb;
        if (info.dereverb) dereverb(b, p.dereverbRt60 > 0 ? p.dereverbRt60 : std::min(info.rt60, 1.5f));
    }

    // Passe 3 : la voix est propre, on la remesure. Tout ce qui suit (de-esser, nivellement, compresseur,
    // loudness) se regle sur le vrai fond apres debruitage, pas sur le bruit de la prise brute.
    info.clean = analyze(b);
    const AnalysisResult& cl = info.clean;
    // parole = au-dessus du fond propre, et a moins de 25 dB de la voix mediane (apres DeepFilterNet le fond
    // tombe tres bas : un residu a -90 dBFS n'est pas une personne lointaine)
    const float speechDb = std::max(cl.noiseDb + 12.0f, cl.voiceDb - 25.0f);

    // 4. EQ du preset
    const struct { Biquad::Type t; float f, db, q; } eq[] = {
        {Biquad::LowShelf, p.lowShelfHz, p.lowShelfDb, 0.7071f},
        {Biquad::Peak, p.peak1Hz, p.peak1Db, p.peak1Q},
        {Biquad::Peak, p.peak2Hz, p.peak2Db, p.peak2Q},
        {Biquad::HighShelf, p.highShelfHz, p.highShelfDb, 0.7071f},
    };
    for (const auto& e : eq) {
        if (!p.use.eq || e.db == 0 || e.f >= sr / 2) continue;
        Biquad f;
        f.design(e.t, sr, e.f, e.q, e.db);
        for (auto& c : b.ch) f.process(c);
    }

    // 5. De-esser : seuil relatif a la sibilance moyenne mesuree
    info.deessThreshDb = cl.sibilanceDb + p.deessOffsetDb;
    if (p.use.deesser) deess(b, info.deessThreshDb);

    // 5b. Nivellement multi-personnes (lent), puis plage dynamique remesuree pour doser le compresseur
    if (p.use.leveler) info.levelDb = levelSpeakers(b, speechDb, p.levelMaxDb);
    float p10, p95;
    loudnessRange(b, p10, p95);

    // 6. Compresseur 3:1 : seuil place pour ramener la plage restante vers celle du preset
    const float excess = std::clamp(p95 - p10 - p.targetRangeDb, 0.0f, 12.0f);
    if (p.use.compressor && excess > 0) {
        const auto env = peakEnvDb(b);
        std::vector<float> voiced;
        for (size_t i = 0; i < env.size(); i += 48) // ponytail: 1 echantillon sur 48 suffit pour un percentile
            if (env[i] > speechDb) voiced.push_back(env[i]);
        // a 3:1, un niveau X dB au-dessus du seuil perd 2X/3 dB : on veut perdre `excess` au p95
        info.compThreshDb = pct(voiced, 0.95f) - 1.5f * excess;
        info.compressed = !voiced.empty();
        if (info.compressed) compress(b, env, info.compThreshDb, 3.0f);
    }

    // 7-8. Loudness -14 LUFS puis limiteur true peak ; le limiteur fait perdre un peu de loudness -> on recale
    const double l0 = integratedLufs(b);
    if (!std::isfinite(l0)) return info; // silence
    const float ceiling = dbToLin(p.ceilingDbtp - 0.1f); // marge pour l'interpolation entre gains
    const auto pre = b.ch;
    constexpr float kMaxGainDb = 36; // prise lointaine debruitee : la voix propre peut etre a -36 LUFS
    float gainDb = p.use.loudness ? std::clamp(p.targetLufs - (float)l0, -20.0f, kMaxGainDb) : 0.0f;
    // Secante : sous limitation, +1 dB de gain ne rapporte qu'une fraction de dB de loudness
    float prevGain = 0, prevLufs = (float)l0;
    for (int it = 0; it < 8; ++it) { // ponytail: 8 passes de limiteur max ; converge en 2-3 sur de la voix
        b.ch = pre;
        for (auto& c : b.ch)
            for (float& v : c) v *= dbToLin(gainDb);
        if (p.use.limiter) limitTruePeak(b, ceiling);
        info.outLufs = (float)integratedLufs(b);
        const float err = p.targetLufs - info.outLufs;
        if (!p.use.loudness || std::fabs(err) < 0.1f || gainDb >= kMaxGainDb) break;
        const float slope = std::clamp((info.outLufs - prevLufs) / (gainDb - prevGain + 1e-6f), 0.1f, 1.0f);
        prevGain = gainDb;
        prevLufs = info.outLufs;
        gainDb = std::min(kMaxGainDb, gainDb + err / slope);
    }
    info.gainDb = gainDb;
    info.outTruePeakDb = linToDb(truePeak(b));
    return info;
}

std::string ProcessInfo::toJson() const
{
    char s[512];
    std::snprintf(s, sizeof(s),
                  "{\"denoiser\":\"%s\",\"denoiseWet\":%.2f,\"levelDb\":%.1f,\"rt60\":%.2f,\"dereverb\":%s,\"deessThreshDb\":%.1f,\"compressed\":%s,"
                  "\"compThreshDb\":%.1f,\"gainDb\":%.2f,\"outLufs\":%.2f,\"outTruePeakDb\":%.2f}",
                  denoiser, denoiseWet, levelDb, rt60, dereverb ? "true" : "false", deessThreshDb, compressed ? "true" : "false", compThreshDb, gainDb,
                  outLufs, outTruePeakDb);
    return std::string(s, s + std::strlen(s) - 1) + ",\"clean\":" + clean.toJson() + "}";
}

} // namespace gv
