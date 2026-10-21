#include "Dereverb.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include "kiss_fft.h"

void Dereverb::FftFree::operator()(kiss_fft_state* s) const { opus_fft_free(s, 0); }

void Dereverb::prepare()
{
    if (!fft_) fft_.reset(opus_fft_alloc(2 * kFrame, nullptr, nullptr, 0));
    win_.resize(2 * kFrame);
    for (int k = 0; k < 2 * kFrame; ++k) // racine de Hann periodique : analyse x synthese = Hann, somme 1 a 50 %
        win_[k] = (float)std::sin(3.14159265358979323846 * k / (2 * kFrame));
    prev_.assign(kFrame, 0.0f);
    ola_.assign(kFrame, 0.0f);
    pow_.assign((size_t)kHist * kBins, 0.0f);
    idx_ = 0;
}

// Queue en -60 dB / RT60 ; apres 50 ms elle vaut 10^(-6 * 0.05 / RT60) en puissance
// ponytail: bascule on/off sans fondu (clic possible) ; reglage pose par la calibration, pas en direct
void Dereverb::setRt60(float s) { decay_ = s > 0.05f ? std::pow(10.0f, -6.0f * 0.01f * kHist / s) : 0.0f; }

void Dereverb::process(float* frame)
{
    constexpr int N = 2 * kFrame;
    constexpr float kFloor = 0.0625f, kSmooth = 0.5f; // gain mini -12 dB (evite les artefacts "musicaux")
    if (decay_ <= 0.0f) { // hors service : simple retard d'une trame
        float tmp[kFrame];
        std::memcpy(tmp, frame, sizeof(tmp));
        std::memcpy(frame, prev_.data(), sizeof(tmp));
        std::memcpy(prev_.data(), tmp, sizeof(tmp));
        return;
    }
    kiss_fft_cpx in[N], X[N];
    for (int k = 0; k < N; ++k) {
        in[k].r = win_[k] * (k < kFrame ? prev_[k] : frame[k - kFrame]);
        in[k].i = 0.0f;
    }
    std::memcpy(prev_.data(), frame, kFrame * sizeof(float));
    opus_fft(fft_.get(), in, X, 0);
    float* pOld = &pow_[(size_t)idx_ * kBins]; // trame t - kHist, remplacee par t
    const float* pPrev = &pow_[(size_t)((idx_ + kHist - 1) % kHist) * kBins];
    idx_ = (idx_ + 1) % kHist;
    for (int b = 0; b < kBins; ++b) {
        const float p = X[b].r * X[b].r + X[b].i * X[b].i;
        const float g = std::sqrt(std::max(kFloor, 1.0f - decay_ * pOld[b] / (p + 1e-12f)));
        pOld[b] = p + kSmooth * (pPrev[b] - p);
        X[b].r *= g; X[b].i *= g;
        if (b > 0 && b < kFrame) X[N - b] = {X[b].r, -X[b].i}; // symetrie hermitienne
    }
    opus_ifft(fft_.get(), X, in, 0);
    for (int k = 0; k < kFrame; ++k) {
        frame[k] = ola_[k] + win_[k] * in[k].r;
        ola_[k] = win_[k + kFrame] * in[k + kFrame].r;
    }
}

// ponytail: heuristique (pas le modele de Lollmann) ; percentile et fenetres a recaler sur de vraies pieces
float estimateRt60(const std::vector<float>& x, double sr, float floorDb)
{
    constexpr size_t W = 8; // 8 blocs de 10 ms
    const size_t b10 = (size_t)(0.01 * sr);
    std::vector<double> ms;
    for (size_t s = 0; s + b10 <= x.size(); s += b10) {
        double e = 0;
        for (size_t i = s; i < s + b10; ++i) e += (double)x[i] * x[i];
        ms.push_back(e / b10);
    }
    std::vector<float> e;
    for (size_t i = 2; i < ms.size(); ++i) // lisse sur 30 ms
        e.push_back((float)(10.0 * std::log10((ms[i] + ms[i - 1] + ms[i - 2]) / 3 + 1e-20)));
    std::vector<float> rts;
    for (size_t i = 0; i + W <= e.size(); ++i) {
        if (e[i + W - 1] < floorDb + 10.0f || e[i] - e[i + W - 1] < 3.0f) continue; // au-dessus du bruit, qui chute
        double num = 0, den = 0, me = 0;
        for (size_t k = 0; k < W; ++k) me += e[i + k] / W;
        for (size_t k = 0; k < W; ++k) { num += (k - (W - 1) / 2.0) * (e[i + k] - me); den += (k - (W - 1) / 2.0) * (k - (W - 1) / 2.0); }
        rts.push_back((float)(-60.0 / (num / den / 0.01)));
    }
    if (rts.size() < 20) return 0.0f;
    const size_t k = rts.size() / 10; // percentile 10 %
    std::nth_element(rts.begin(), rts.begin() + k, rts.end());
    return rts[k];
}
