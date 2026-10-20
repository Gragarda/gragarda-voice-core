#include "gv/PopDetector.h"
#include <algorithm>
#include <cmath>
#include "gv/Dsp.h"

namespace gv {

namespace {
constexpr float kCutFactor = 8.0f, kScanFactor = 20.0f;
constexpr float kDcStep = 0.02f;      // ~ -34 dBFS
constexpr float kSilenceRms = 0.003f; // ~ -50 dBFS
constexpr float kRefFloor = 1e-4f;    // evite x8 sur du bruit numerique en plein silence

size_t ms(const AudioBuffer& b, double v) { return (size_t)std::lround(v * 0.001 * b.sampleRate); }

float median(std::vector<float>& v)
{
    if (v.empty()) return 0.0f;
    std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
    return v[v.size() / 2];
}

// Niveau moyen avant/apres p : moyenne ponderee Hann sur 40 ms de chaque cote. La spec dit moyenne sur 5 ms,
// mais la voix (forme d'onde asymetrique, f0 ~100-250 Hz) donne alors des "sauts DC" partout ; la fenetre
// de Hann rejette les sinusoides au-dela de ~50 Hz, seul l'offset passe.
float dcStepAt(const AudioBuffer& b, size_t p)
{
    float step = 0;
    for (const auto& c : b.ch) step = std::max(step, std::fabs(dcDelta(c, p, ms(b, 40))));
    return step;
}
} // namespace

float dcDelta(const std::vector<float>& c, size_t p, size_t w)
{
    if (p < w || p + w > c.size()) return 0.0f;
    double before = 0, after = 0, norm = 0;
    for (size_t k = 0; k < w; ++k) {
        const double h = 0.5 - 0.5 * std::cos(2 * kPi * (k + 0.5) / w);
        before += h * c[p - w + k];
        after += h * c[p + k];
        norm += h;
    }
    return (float)((after - before) / norm);
}

namespace {
// Ecart-type autour de p (chaque cote centre sur sa moyenne : un saut DC en plein silence reste du silence)
float rmsAround(const AudioBuffer& b, size_t p)
{
    const size_t w = ms(b, 20), n = b.frames();
    if (p < w || p + w > n) return 1.0f;
    double s = 0;
    for (const auto& c : b.ch)
        for (size_t a : {p - w, p}) {
            double m = 0, q = 0;
            for (size_t i = a; i < a + w; ++i) { m += c[i]; q += (double)c[i] * c[i]; }
            m /= w;
            s += q / w - m * m;
        }
    return (float)std::sqrt(s / (2.0 * b.ch.size()));
}

// Un saut entre x[c-1] et x[c] rend j grand en c ET c+1 : on ramene p sur c (premier echantillon du plan suivant)
size_t refine(const AudioBuffer& b, size_t p)
{
    return p > 2 && jumpAt(b, p - 1) > 0.5f * jumpAt(b, p) ? p - 1 : p;
}

// Type de plop a p, sachant si le saut d'echantillon depasse le seuil
bool classify(const AudioBuffer& b, size_t p, bool click, bool allowDcOnly, Pop& out)
{
    const float dc = dcStepAt(b, p);
    const bool isDc = dc > kDcStep;
    if (!click && !(isDc && allowDcOnly))
        return false;
    out.pos = p;
    if (click && rmsAround(b, p) < kSilenceRms) out.type = PopType::Silence;
    else if (click && isDc) out.type = PopType::Both;
    else if (click) out.type = PopType::Click;
    else { out.type = PopType::DcStep; out.severity = linToDb(dc); }
    return true;
}
} // namespace

float jumpAt(const AudioBuffer& b, size_t i)
{
    if (i < 2 || i >= b.frames()) return 0.0f;
    float j = 0;
    for (const auto& c : b.ch)
        j = std::max(j, std::fabs(c[i] - (2.0f * c[i - 1] - c[i - 2])));
    return j;
}

float localJumpRef(const AudioBuffer& b, size_t p)
{
    const size_t near = ms(b, 2), far = ms(b, 12), n = b.frames();
    std::vector<float> v;
    v.reserve(2 * (far - near));
    for (size_t i = p > far ? p - far : 2; i + near < p; ++i) v.push_back(jumpAt(b, i));
    for (size_t i = p + near + 1; i <= p + far && i < n; ++i) v.push_back(jumpAt(b, i));
    return std::max(kRefFloor, median(v));
}

std::vector<Pop> detectPops(const AudioBuffer& b, const std::vector<double>& cutsSec)
{
    std::vector<Pop> pops;
    const size_t n = b.frames(), dedupe = ms(b, 5);
    if (n < 4) return pops;
    auto taken = [&](size_t p) {
        for (const Pop& q : pops)
            if ((p > q.pos ? p - q.pos : q.pos - p) < dedupe) return true;
        return false;
    };

    // 1. Zones de cut : +-60 ms (decalage audio/video jusqu'a ~1 image), seuil x8 ou saut DC
    const size_t win = ms(b, 60);
    for (double t : cutsSec) {
        const size_t c = (size_t)std::max(0.0, t * b.sampleRate);
        if (c >= n) continue;
        size_t p = 0;
        float best = -1;
        for (size_t i = c > win ? c - win : 2; i < std::min(n, c + win); ++i) {
            const float j = jumpAt(b, i);
            if (j > best) { best = j; p = i; }
        }
        p = refine(b, p);
        const float ref = localJumpRef(b, p);
        Pop pop;
        pop.origin = PopOrigin::Cut;
        pop.severity = best / ref;
        if (!taken(p) && classify(b, p, best > kCutFactor * ref, true, pop))
            pops.push_back(pop);
    }

    // 2. Scan du fichier par blocs de 10 ms, seuil x20 de la plus grande mediane (bloc, voisins) :
    // une consonne explosive releve la mediane de son propre bloc, un saut d'un echantillon non.
    const size_t blk = ms(b, 10), nb = n / blk;
    std::vector<float> med(nb), peak(nb), tmp(blk);
    std::vector<size_t> arg(nb);
    for (size_t k = 0; k < nb; ++k) {
        peak[k] = -1;
        for (size_t i = 0; i < blk; ++i) {
            tmp[i] = jumpAt(b, k * blk + i);
            if (tmp[i] > peak[k]) { peak[k] = tmp[i]; arg[k] = k * blk + i; }
        }
        med[k] = median(tmp);
    }
    for (size_t k = 0; k < nb; ++k) {
        const float ref = std::max({kRefFloor, med[k], k ? med[k - 1] : 0.0f, k + 1 < nb ? med[k + 1] : 0.0f});
        if (peak[k] <= kScanFactor * ref || taken(arg[k])) continue;
        Pop pop;
        pop.origin = PopOrigin::Scan;
        pop.severity = peak[k] / ref;
        if (classify(b, refine(b, arg[k]), true, false, pop))
            pops.push_back(pop);
    }
    std::sort(pops.begin(), pops.end(), [](const Pop& a, const Pop& c) { return a.pos < c.pos; });
    return pops;
}

const char* toString(PopType t)
{
    switch (t) {
    case PopType::Click: return "click";
    case PopType::DcStep: return "dc_step";
    case PopType::Both: return "both";
    default: return "silence";
    }
}

} // namespace gv
