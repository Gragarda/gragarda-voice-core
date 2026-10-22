#include "gv/Dsp.h"
#include <algorithm>
#include <limits>

namespace gv {

// Formules RBJ (Audio EQ Cookbook), reprises de gragarda-voice-core
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
    case LowPass:
        B0 = (1 - cw) / 2; B1 = 1 - cw; B2 = (1 - cw) / 2;
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
    case Notch:
        B0 = 1; B1 = -2 * cw; B2 = 1;
        A0 = 1 + alpha; A1 = -2 * cw; A2 = 1 - alpha;
        break;
    default: // HighShelf
        B0 = A * ((A + 1) + (A - 1) * cw + sA); B1 = -2 * A * ((A - 1) + (A + 1) * cw);
        B2 = A * ((A + 1) + (A - 1) * cw - sA);
        A0 = (A + 1) - (A - 1) * cw + sA; A1 = 2 * ((A - 1) - (A + 1) * cw);
        A2 = (A + 1) - (A - 1) * cw - sA;
        break;
    }
    b0 = B0 / A0; b1 = B1 / A0; b2 = B2 / A0; a1 = A1 / A0; a2 = A2 / A0;
}

void Biquad::process(std::vector<float>& x) const
{
    double z0 = 0, z1 = 0; // forme directe transposee II
    for (float& v : x) {
        const double in = v, y = b0 * in + z0;
        z0 = b1 * in - a1 * y + z1;
        z1 = b2 * in - a2 * y;
        v = (float)y;
    }
}

void highPass1(std::vector<float>& x, double sr, double f)
{
    const double a = std::exp(-2.0 * kPi * f / sr);
    double xp = 0, yp = 0;
    for (float& v : x) {
        const double y = 0.5 * (1 + a) * (v - xp) + a * yp;
        xp = v;
        yp = y;
        v = (float)y;
    }
}

namespace {
// Energie ponderee K (somme des canaux) par segment de 100 ms. Coefficients exacts BS.1770 pour tout sr
// (meme calcul que libebur128).
std::vector<double> kSegments(const AudioBuffer& b, size_t& segLen)
{
    const double sr = b.sampleRate;
    segLen = (size_t)std::lround(sr * 0.1);
    double K = std::tan(kPi * 1681.974450955533 / sr), Q = 0.7071752369554196;
    const double Vh = std::pow(10.0, 3.999843853973347 / 20.0), Vb = std::pow(Vh, 0.4996667741545416);
    double a0 = 1 + K / Q + K * K;
    Biquad s1;
    s1.b0 = (Vh + Vb * K / Q + K * K) / a0; s1.b1 = 2 * (K * K - Vh) / a0; s1.b2 = (Vh - Vb * K / Q + K * K) / a0;
    s1.a1 = 2 * (K * K - 1) / a0; s1.a2 = (1 - K / Q + K * K) / a0;
    K = std::tan(kPi * 38.13547087602444 / sr); Q = 0.5003270373238773;
    a0 = 1 + K / Q + K * K;
    Biquad s2;
    s2.b0 = 1; s2.b1 = -2; s2.b2 = 1; s2.a1 = 2 * (K * K - 1) / a0; s2.a2 = (1 - K / Q + K * K) / a0;

    const size_t n = b.frames(), nseg = n / segLen;
    std::vector<double> seg(nseg, 0.0);
    for (const auto& c : b.ch) {
        double z[4] = {};
        for (size_t i = 0; i < nseg * segLen; ++i) {
            const double x = c[i];
            const double y1 = s1.b0 * x + z[0];
            z[0] = s1.b1 * x - s1.a1 * y1 + z[1];
            z[1] = s1.b2 * x - s1.a2 * y1;
            const double y2 = s2.b0 * y1 + z[2];
            z[2] = s2.b1 * y1 - s2.a1 * y2 + z[3];
            z[3] = s2.b2 * y1 - s2.a2 * y2;
            seg[i / segLen] += y2 * y2;
        }
    }
    return seg; // ponytail: canaux surround non ponderes (G=1 partout) ; ok pour mono/stereo telephone
}

double toLufs(double ms) { return -0.691 + 10.0 * std::log10(ms + 1e-30); }
} // namespace

double integratedLufs(const AudioBuffer& b)
{
    size_t segLen = 0;
    const auto seg = kSegments(b, segLen);
    std::vector<double> blocks; // 400 ms, recouvrement 75 %
    for (size_t j = 0; j + 4 <= seg.size(); ++j)
        blocks.push_back((seg[j] + seg[j + 1] + seg[j + 2] + seg[j + 3]) / (4.0 * segLen));
    double sum = 0;
    int cnt = 0;
    for (double e : blocks)
        if (toLufs(e) > -70.0) { sum += e; ++cnt; }
    if (!cnt)
        return -std::numeric_limits<double>::infinity();
    const double rel = toLufs(sum / cnt) - 10.0;
    sum = 0;
    cnt = 0;
    for (double e : blocks)
        if (toLufs(e) > -70.0 && toLufs(e) > rel) { sum += e; ++cnt; }
    return cnt ? toLufs(sum / cnt) : -std::numeric_limits<double>::infinity();
}

std::vector<float> shortTermLufs(const AudioBuffer& b)
{
    size_t segLen = 0;
    const auto seg = kSegments(b, segLen);
    std::vector<float> out;
    double acc = 0;
    for (size_t j = 0; j < seg.size(); ++j) {
        acc += seg[j];
        if (j >= 30) acc -= seg[j - 30];
        if (j >= 29) out.push_back((float)toLufs(acc / (30.0 * segLen)));
    }
    return out;
}

std::vector<double> welch(const std::vector<float>& x, size_t n, const std::vector<size_t>& frames)
{
    std::vector<double> out(n / 2 + 1, 0.0);
    const size_t nf = x.size() / n;
    std::vector<size_t> idx = frames;
    if (idx.empty())
        for (size_t f = 0; f < nf; ++f) idx.push_back(f);
    std::vector<double> re(n), im(n);
    size_t used = 0;
    for (size_t f : idx) {
        if (f >= nf) continue;
        for (size_t i = 0; i < n; ++i) {
            re[i] = x[f * n + i] * (0.5 - 0.5 * std::cos(2 * kPi * i / n));
            im[i] = 0;
        }
        // FFT radix 2 en place
        for (size_t i = 1, j = 0; i < n; ++i) {
            size_t bit = n >> 1;
            for (; j & bit; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) { std::swap(re[i], re[j]); std::swap(im[i], im[j]); }
        }
        for (size_t len = 2; len <= n; len <<= 1) {
            const double a = -2 * kPi / len;
            for (size_t i = 0; i < n; i += len)
                for (size_t k = 0; k < len / 2; ++k) {
                    const double c = std::cos(a * k), s = std::sin(a * k);
                    const size_t u = i + k, v = i + k + len / 2;
                    const double tr = re[v] * c - im[v] * s, ti = re[v] * s + im[v] * c;
                    re[v] = re[u] - tr; im[v] = im[u] - ti;
                    re[u] += tr; im[u] += ti;
                }
        }
        for (size_t k = 0; k <= n / 2; ++k) out[k] += re[k] * re[k] + im[k] * im[k];
        ++used;
    }
    for (double& v : out) v /= std::max<size_t>(1, used);
    return out;
}

void loudnessRange(const AudioBuffer& b, float& p10, float& p95)
{
    p10 = p95 = -120.0f;
    const auto st = shortTermLufs(b);
    double e = 0;
    int cnt = 0;
    for (float v : st)
        if (v > -70.0f) { e += std::pow(10.0, v / 10.0); ++cnt; }
    if (!cnt) return;
    const double gate = 10.0 * std::log10(e / cnt) - 20.0;
    std::vector<float> g;
    for (float v : st)
        if (v > -70.0f && v > gate) g.push_back(v);
    std::sort(g.begin(), g.end());
    p10 = g[(size_t)(0.10 * (g.size() - 1) + 0.5)];
    p95 = g[(size_t)(0.95 * (g.size() - 1) + 0.5)];
}

std::vector<float> truePeakEnvelope(const AudioBuffer& b)
{
    // Interpolation sinc fenetree (Hann, 16 points) aux positions i + 1/4, 2/4, 3/4
    constexpr int kTaps = 8;
    float coef[3][2 * kTaps];
    for (int p = 1; p <= 3; ++p)
        for (int k = -kTaps + 1; k <= kTaps; ++k) {
            const double t = p / 4.0 - k;
            const double s = std::sin(kPi * t) / (kPi * t);
            coef[p - 1][k + kTaps - 1] = (float)(s * 0.5 * (1 + std::cos(kPi * t / kTaps)));
        }
    const size_t n = b.frames();
    std::vector<float> env(n, 0.0f);
    for (const auto& c : b.ch)
        for (size_t i = 0; i < n; ++i) {
            float m = std::fabs(c[i]);
            if (i + 1 >= kTaps && i + kTaps < n) { // ponytail: bords (8 echantillons) sans interpolation
                const float* x = &c[i - kTaps + 1];
                for (int p = 0; p < 3; ++p) {
                    float v = 0;
                    for (int k = 0; k < 2 * kTaps; ++k) v += x[k] * coef[p][k];
                    m = std::max(m, std::fabs(v));
                }
            }
            env[i] = std::max(env[i], m);
        }
    return env;
}

float truePeak(const AudioBuffer& b)
{
    const auto env = truePeakEnvelope(b);
    return env.empty() ? 0.0f : *std::max_element(env.begin(), env.end());
}

} // namespace gv
