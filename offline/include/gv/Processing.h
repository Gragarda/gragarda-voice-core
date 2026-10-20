#pragma once
#include <string>
#include "gv/Analysis.h"
#include "gv/AudioBuffer.h"

namespace gv {

// Reglages d'un preset (core/presets/*.json, objet JSON plat). Valeurs par defaut = voix neutre.
struct Preset {
    std::string name = "default";
    float hpfHz = 80;           // passe-haut 18 dB/oct
    float denoiseMin = 0;       // dosage minimal du debruitage 0..1 (sinon selon le bruit mesure)
    float expanderDb = 0;       // attenuation max des queues de reverberation (0 = off)
    float lowShelfHz = 150, lowShelfDb = 0;
    float peak1Hz = 400, peak1Db = 0, peak1Q = 1;
    float peak2Hz = 3000, peak2Db = 0, peak2Q = 1;
    float highShelfHz = 10000, highShelfDb = 0;
    float deessOffsetDb = 6;    // le de-esser agit au-dessus de la sibilance moyenne mesuree + offset
    float targetRangeDb = 8;    // plage dynamique visee par le compresseur
    float targetLufs = -14, ceilingDbtp = -1;
};
Preset parsePreset(const std::string& json);

struct ProcessInfo {
    float denoiseWet = 0, expanderDb = 0, deessThreshDb = 0, compThreshDb = 0, gainDb = 0;
    bool compressed = false;
    float outLufs = -120, outTruePeakDb = -120;
    std::string toJson() const;
};

// Passe 2 : HPF -> debruitage -> expandeur -> EQ -> de-esser -> compresseur -> loudness -> limiteur true peak
ProcessInfo processChain(AudioBuffer& b, const Preset& p, const AnalysisResult& a);

// Limiteur hors temps reel : anticipation, crete suréchantillonnee x4 <= ceiling (lineaire)
void limitTruePeak(AudioBuffer& b, float ceiling, float lookaheadMs = 5, float releaseMs = 80);

} // namespace gv
