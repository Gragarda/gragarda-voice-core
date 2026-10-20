#include "gv/Analysis.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include "gv/Dsp.h"

namespace gv {

namespace {
float pct(std::vector<float> v, float p)
{
    if (v.empty()) return -120.0f;
    const size_t k = std::min(v.size() - 1, (size_t)(p * (v.size() - 1) + 0.5f));
    std::nth_element(v.begin(), v.begin() + k, v.end());
    return v[k];
}

std::vector<float> mono(const AudioBuffer& b)
{
    std::vector<float> m(b.frames(), 0.0f);
    for (const auto& c : b.ch)
        for (size_t i = 0; i < m.size(); ++i) m[i] += c[i] / b.channels();
    return m;
}

// Puissance moyenne par fenetre de w echantillons
std::vector<double> power(const std::vector<float>& x, size_t w)
{
    std::vector<double> p(x.size() / w, 0.0);
    for (size_t i = 0; i < p.size() * w; ++i) p[i / w] += (double)x[i] * x[i];
    for (double& v : p) v /= w;
    return p;
}

float db(double pw) { return (float)(10.0 * std::log10(pw + 1e-14)); }
float finite(double v) { return std::isfinite(v) ? (float)v : -120.0f; }
} // namespace

AnalysisResult analyze(const AudioBuffer& b)
{
    AnalysisResult r;
    if (b.frames() < (size_t)b.sampleRate) return r; // < 1 s : rien de fiable
    const double sr = b.sampleRate;
    auto m = mono(b);
    highPass1(m, sr, 20.0); // un offset DC n'est pas du bruit

    // Bruit de fond et niveau de voix (fenetres de 50 ms)
    // ponytail: niveau seul, pas de spectre du bruit (RNNoise n'en a pas besoin) ; l'ajouter pour un debruiteur spectral
    const auto p50 = power(m, (size_t)(sr * 0.05));
    std::vector<float> lv;
    for (double v : p50) lv.push_back(db(v));
    {
        std::vector<float> s = lv;
        std::sort(s.begin(), s.end());
        const size_t q = std::max<size_t>(1, s.size() / 5);
        double acc = 0;
        for (size_t i = 0; i < q; ++i) acc += std::pow(10.0, s[i] / 10.0);
        r.noiseDb = db(acc / q);
    }
    std::vector<float> voice;
    std::vector<bool> isVoice(lv.size());
    for (size_t i = 0; i < lv.size(); ++i)
        if ((isVoice[i] = lv[i] > r.noiseDb + 10.0f && lv[i] > -70.0f)) voice.push_back(lv[i]);
    r.voiceDb = pct(voice, 0.5f);

    r.lufs = finite(integratedLufs(b));
    r.truePeakDb = linToDb(truePeak(b));

    // Plage dynamique : percentiles 10/95 de la loudness court terme, portes absolue -70 et relative -20
    loudnessRange(b, r.stP10, r.stP95);
    r.dynamicRangeDb = r.stP95 - r.stP10;

    // Sibilance : bande 5-9 kHz sur les fenetres de voix
    if (sr > 20000) {
        auto band = m;
        Biquad f;
        f.design(Biquad::HighPass, sr, 5000, 0.7071); f.process(band); f.process(band);
        f.design(Biquad::LowPass, sr, 9000, 0.7071); f.process(band); f.process(band);
        const auto pb = power(band, (size_t)(sr * 0.05));
        double sb = 0, sf = 0;
        for (size_t i = 0; i < pb.size() && i < p50.size(); ++i)
            if (isVoice[i]) { sb += pb[i]; sf += p50[i]; }
        if (sf > 0) r.sibilanceDb = db(sb / sf);
    }

    // Raies tonales stables : spectre moyen (resolution ~6 Hz) des 20 % de fenetres de 170 ms les plus calmes,
    // pic > 15 dB au-dessus de la mediane locale (+-30 Hz), entre 60 Hz et 8 kHz. 5 raies max, les plus fortes.
    {
        const size_t N = 8192;
        const auto pw = power(m, N);
        std::vector<size_t> order(pw.size());
        for (size_t i = 0; i < order.size(); ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&](size_t x, size_t y) { return pw[x] < pw[y]; });
        order.resize(std::max<size_t>(1, order.size() / 5));
        const auto S = welch(m, N, order);
        const double hz = sr / N;
        std::vector<std::pair<float, float>> peaks; // (proeminence dB, Hz)
        const size_t w = (size_t)(30 / hz);
        for (size_t k = (size_t)(60 / hz); k < std::min(S.size() - w, (size_t)(8000 / hz)); ++k) {
            if (S[k] < S[k - 1] || S[k] < S[k + 1]) continue;
            std::vector<float> nb;
            for (size_t i = k - w; i <= k + w; ++i)
                if (i + 2 < k || i > k + 2) nb.push_back((float)S[i]);
            const float prom = db(S[k] / (pct(nb, 0.5f) + 1e-20));
            if (prom > 15) peaks.push_back({prom, (float)(k * hz)});
        }
        std::sort(peaks.rbegin(), peaks.rend());
        for (size_t i = 0; i < peaks.size() && r.toneHz.size() < 5; ++i) r.toneHz.push_back(peaks[i].second);
    }

    // Resonance : temps pour perdre 20 dB apres une fin de mot (enveloppe 10 ms)
    const auto p10 = power(m, (size_t)(sr * 0.01));
    std::vector<float> env;
    for (double v : p10) env.push_back(db(v));
    const float hi = r.voiceDb - 6.0f;
    std::vector<float> decays;
    for (size_t i = 1; i < env.size(); ++i) {
        if (!(env[i - 1] >= hi && env[i] < hi)) continue;
        for (size_t j = i; j < env.size() && j < i + 50; ++j) // 500 ms max
            if (env[j] < env[i - 1] - 20.0f) { decays.push_back((j - i + 1) * 10.0f); i = j; break; }
    }
    if (!decays.empty()) r.decayMs = pct(decays, 0.5f);
    return r;
}

std::string AnalysisResult::toJson() const
{
    std::string tones = "[";
    for (size_t i = 0; i < toneHz.size(); ++i) tones += (i ? "," : "") + std::to_string((int)std::lround(toneHz[i]));
    tones += "]";
    char s[512];
    std::snprintf(s, sizeof(s),
                  "{\"noiseDb\":%.1f,\"voiceDb\":%.1f,\"lufs\":%.2f,\"truePeakDb\":%.2f,\"stP10\":%.1f,"
                  "\"stP95\":%.1f,\"dynamicRangeDb\":%.1f,\"sibilanceDb\":%.1f,\"decayMs\":%.0f,\"toneHz\":%s}",
                  noiseDb, voiceDb, lufs, truePeakDb, stP10, stP95, dynamicRangeDb, sibilanceDb, decayMs, tones.c_str());
    return s;
}

} // namespace gv
