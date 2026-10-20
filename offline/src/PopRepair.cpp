#include "gv/PopRepair.h"
#include <algorithm>
#include <cmath>
#include "gv/Dsp.h"

namespace gv {

namespace {
size_t ms(const AudioBuffer& b, double v) { return (size_t)std::lround(v * 0.001 * b.sampleRate); }

// Micro fondu sortant puis entrant de 5 ms : inaudible en plein silence
void fadeGap(std::vector<float>& x, size_t p, size_t w)
{
    for (size_t k = 0; k < w; ++k) {
        const float g = (float)k / w;
        x[p - 1 - k] *= g;
        x[p + k] *= g;
    }
}

// Rampe : l'offset d'apres la coupe part du niveau d'avant et rejoint le sien en w echantillons
void dcRamp(std::vector<float>& x, size_t p, size_t w40, size_t w)
{
    const float d = dcDelta(x, p, w40);
    for (size_t k = 0; k < w; ++k) x[p + k] -= d * (1.0f - (float)k / w);
}

// Pont cubique d'Hermite sur ]a, b[ : raccorde valeurs et pentes des deux bords
void hermite(std::vector<float>& x, size_t a, size_t b)
{
    const float L = (float)(b - a), y0 = x[a], y1 = x[b];
    const float m0 = (x[a] - x[a - 1]) * L, m1 = (x[b + 1] - x[b]) * L;
    for (size_t i = a + 1; i < b; ++i) {
        const float t = (i - a) / L, t2 = t * t, t3 = t2 * t;
        x[i] = (2 * t3 - 3 * t2 + 1) * y0 + (t3 - 2 * t2 + t) * m0 + (-2 * t3 + 3 * t2) * y1 + (t3 - t2) * m1;
    }
}

// Prediction lineaire (autocorrelation + Levinson-Durbin) : x[t] ~ sum a[k] x[t-1-k]
std::vector<double> lpc(const float* x, size_t n, int order)
{
    std::vector<double> r(order + 1, 0.0), a(order, 0.0), tmp(order);
    for (int k = 0; k <= order; ++k)
        for (size_t i = k; i < n; ++i) r[k] += (double)x[i] * x[i - k];
    r[0] = r[0] * 1.0001 + 1e-12; // leger bruit blanc : stabilite
    double err = r[0];
    for (int i = 0; i < order; ++i) {
        double acc = r[i + 1];
        for (int k = 0; k < i; ++k) acc -= a[k] * r[i - k];
        const double kk = acc / err;
        tmp = a;
        a[i] = kk;
        for (int k = 0; k < i; ++k) a[k] = tmp[k] - kk * tmp[i - 1 - k];
        err *= 1 - kk * kk;
        if (err <= 0) break;
    }
    return a;
}

// Prolonge h (historique, ordre chronologique) de m echantillons
std::vector<float> extrapolate(std::vector<float> h, size_t m, int order)
{
    const auto a = lpc(h.data(), h.size(), order);
    const size_t n0 = h.size();
    for (size_t t = 0; t < m; ++t) {
        double v = 0;
        for (int k = 0; k < order; ++k) v += a[k] * h[n0 + t - 1 - k];
        h.push_back((float)v);
    }
    return std::vector<float>(h.begin() + n0, h.end());
}

// Fondu croise a puissance constante centre sur p, entre l'avant prolonge par LPC vers l'avant
// et l'apres prolonge par LPC vers l'arriere : chaque cote reste lisse, pas de coude
void crossfade(std::vector<float>& x, size_t p, size_t half, size_t hist)
{
    constexpr int kOrder = 32;
    const auto fwd = extrapolate(std::vector<float>(x.begin() + (p - hist), x.begin() + p), half, kOrder);
    std::vector<float> after(x.begin() + p, x.begin() + (p + hist));
    std::reverse(after.begin(), after.end());
    const auto bwd = extrapolate(after, half, kOrder); // bwd[k] ~ x[p-1-k]
    for (size_t k = 0; k < 2 * half; ++k) {
        const size_t i = p - half + k;
        const float t = (k + 0.5f) / (2 * half);
        const float a = i < p ? x[i] : fwd[i - p], c = i >= p ? x[i] : bwd[p - 1 - i];
        x[i] = std::cos(t * (float)kPi / 2) * a + std::sin(t * (float)kPi / 2) * c;
    }
}

// Plus grand saut "naturel" autour de p (hors +-2 ms)
float neighborPeak(const AudioBuffer& b, size_t p)
{
    float j = 0;
    for (size_t i = p - ms(b, 12); i < p + ms(b, 12); ++i)
        if (i + ms(b, 2) < p || i > p + ms(b, 2)) j = std::max(j, jumpAt(b, i));
    return j;
}

float residual(const AudioBuffer& b, size_t p)
{
    float j = 0;
    for (size_t i = p - ms(b, 2); i <= p + ms(b, 2); ++i) j = std::max(j, jumpAt(b, i));
    return j;
}
} // namespace

void repairPops(AudioBuffer& b, const std::vector<Pop>& pops)
{
    const size_t n = b.frames(), margin = ms(b, 45);
    const size_t w40 = ms(b, 40), dcLen = ms(b, 15), half = ms(b, 0.5), xfHalf = ms(b, 4), fade = ms(b, 5), hist = ms(b, 20);
    for (const Pop& pop : pops) {
        const size_t p = pop.pos;
        if (p < margin || p + margin > n) continue; // le detecteur ignore deja les 50 ms des bords
        const float ref = localJumpRef(b, p);
        if (pop.type == PopType::Silence) {
            for (auto& c : b.ch) fadeGap(c, p, fade);
            continue;
        }
        if (pop.type == PopType::DcStep || pop.type == PopType::Both)
            for (auto& c : b.ch) dcRamp(c, p, w40, dcLen);
        if (pop.type == PopType::DcStep)
            continue;
        if (pop.type == PopType::Click) {
            std::vector<std::vector<float>> saved;
            for (auto& c : b.ch) saved.emplace_back(c.begin() + (p - margin), c.begin() + (p + margin));
            // ponytail: Hermite <= 1 ms ; interpolation AR/LPC ordre 32 si des clics longs restent audibles
            for (auto& c : b.ch) hermite(c, p - half - 1, p + half);
            // la voix elle-meme a des pics de j > 3x sa mediane : on tolere aussi le pic naturel voisin
            if (residual(b, p) <= std::max(3 * ref, 1.2f * neighborPeak(b, p)))
                continue;
            for (size_t c = 0; c < saved.size(); ++c) // garde-fou : l'interpolation ne suffit pas -> fondu croise
                std::copy(saved[c].begin(), saved[c].end(), b.ch[c].begin() + (p - margin));
        }
        for (auto& c : b.ch) crossfade(c, p, xfHalf, hist);
    }
}

} // namespace gv
