#pragma once
#include <cmath>
#include <vector>
#include "gv/AudioBuffer.h"

// Briques communes : filtres, mesures BS.1770 (LUFS, true peak).
namespace gv {

constexpr double kPi = 3.14159265358979323846;
inline float dbToLin(float db) { return std::pow(10.0f, db / 20.0f); }
inline float linToDb(float x) { return 20.0f * std::log10(x + 1e-12f); }
inline float coefMs(float ms, double sr) { return std::exp(-1.0f / (0.001f * ms * (float)sr)); }

struct Biquad {
    double b0 = 1, b1 = 0, b2 = 0, a1 = 0, a2 = 0;
    enum Type { HighPass, LowPass, Peak, LowShelf, HighShelf, Notch };
    void design(Type t, double sr, double f, double q, double gainDb = 0);
    void process(std::vector<float>& x) const; // etat neuf, en place
};

void highPass1(std::vector<float>& x, double sr, double f); // 1er ordre, en place

// Loudness integree BS.1770-4 (portes -70 LUFS absolue et -10 LU relative). -inf si silence.
double integratedLufs(const AudioBuffer& b);
// Loudness court terme (fenetre 3 s, pas 100 ms)
std::vector<float> shortTermLufs(const AudioBuffer& b);
// Densite spectrale de puissance (Welch, Hann) de x, fenetres de n (puissance de 2) ; out[k] pour k*sr/n Hz.
// Seules les fenetres dont l'index est dans `frames` sont moyennees (vide = toutes).
std::vector<double> welch(const std::vector<float>& x, size_t n, const std::vector<size_t>& frames = {});

// Percentiles 10 / 95 de la loudness court terme, portes absolue -70 et relative -20 (comme LRA). -120 si silence.
void loudnessRange(const AudioBuffer& b, float& p10, float& p95);
// Crete par echantillon (max canaux, suréchantillonnage x4) et true peak global, lineaire
std::vector<float> truePeakEnvelope(const AudioBuffer& b);
float truePeak(const AudioBuffer& b);

} // namespace gv
