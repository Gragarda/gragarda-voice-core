#pragma once
#include <functional>
#include <string>
#include <vector>
#include "gv/Analysis.h"
#include "gv/AudioBuffer.h"

namespace gv {

// Reglages d'un preset (core/presets/*.json, objet JSON plat). Valeurs par defaut = voix neutre.
struct Preset {
    std::string name = "default";
    float hpfHz = 80;           // passe-haut 18 dB/oct
    float denoiseMin = 0;       // dosage minimal du debruitage 0..1 (sinon selon le bruit mesure)
    float levelMaxDb = 10;      // nivellement multi-personnes : correction max
    float lowShelfHz = 150, lowShelfDb = 0;
    float peak1Hz = 400, peak1Db = 0, peak1Q = 1;
    float peak2Hz = 3000, peak2Db = 0, peak2Q = 1;
    float highShelfHz = 10000, highShelfDb = 0;
    float deessOffsetDb = 6;    // le de-esser agit au-dessus de la sibilance moyenne mesuree + offset
    float targetRangeDb = 8;    // plage dynamique visee par le compresseur
    float targetLufs = -14, ceilingDbtp = -1;
    // Options avancees : cles JSON "pops", "hpf", ... a 0 pour couper un etage (le web les ajoute en tete du preset)
    struct Stages {
        bool pops = true, hpf = true, notch = true, denoise = true, dereverb = true, eq = true, deesser = true,
             leveler = true, compressor = true, loudness = true, limiter = true;
    } use;
};
Preset parsePreset(const std::string& json);

struct ProcessInfo {
    float denoiseWet = 0, deessThreshDb = 0, compThreshDb = 0, gainDb = 0;
    bool compressed = false;
    float rt60 = 0;        // reverberation estimee sur la voix debruitee, s (0 = non mesurable)
    bool dereverb = false; // de-reverb applique (rt60 > 0.3 s)
    const char* denoiser = "rnnoise";
    float levelDb = 0; // correction max du nivellement multi-personnes
    float outLufs = -120, outTruePeakDb = -120;
    AnalysisResult clean; // passe 3 : mesures sur la voix nettoyee (apres debruitage)
    std::string toJson() const;
};

// Nettoyage (regle par l'analyse de la prise brute `a`) : HPF -> coupe-bandes (raies) -> debruitage -> de-reverb,
// puis passe 3 = nouvelle analyse sur la voix propre, qui regle : EQ -> de-esser -> nivellement -> compresseur -> loudness -> limiteur true peak
// Debruiteur externe optionnel (mono 48 kHz, en place, latence deja compensee) a la place de RNNoise
using Denoiser = std::function<void(std::vector<float>&)>;
ProcessInfo processChain(AudioBuffer& b, const Preset& p, const AnalysisResult& a, const Denoiser& externalDenoise = {});

// Limiteur hors temps reel : anticipation, crete suréchantillonnee x4 <= ceiling (lineaire)
void limitTruePeak(AudioBuffer& b, float ceiling, float lookaheadMs = 5, float releaseMs = 80);

} // namespace gv
