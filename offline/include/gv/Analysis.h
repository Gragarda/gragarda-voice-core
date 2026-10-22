#pragma once
#include <string>
#include <vector>
#include "gv/AudioBuffer.h"

namespace gv {

// Mesures d'une piste. Appelee deux fois : sur la prise brute de-plopee (regle le nettoyage) puis sur la voix
// nettoyee (regle la dynamique), cf. processChain.
struct AnalysisResult {
    float noiseDb = -120;     // RMS dBFS des 20 % de fenetres de 50 ms les plus calmes
    float voiceDb = -120;     // RMS median des fenetres de voix (> bruit + 10 dB)
    float lufs = -120;        // loudness integree BS.1770
    float truePeakDb = -120;
    float stP10 = -120, stP95 = -120, dynamicRangeDb = 0; // loudness court terme, percentiles (portes LRA)
    float sibilanceDb = -60;  // energie 5-9 kHz / energie totale sur la voix
    std::vector<float> toneHz;  // raies stables du bruit de fond (sifflement, ventilation, secteur...)
    float decayMs = 0;        // chute de 20 dB apres les fins de mots (mediane) : resonance de la piece
    std::string toJson() const;
};

AnalysisResult analyze(const AudioBuffer& b);

} // namespace gv
