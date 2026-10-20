#pragma once
#include <string>
#include <vector>
#include "gv/Analysis.h"
#include "gv/PopDetector.h"
#include "gv/Processing.h"

namespace gv {

struct Report {
    int sampleRate = 48000;
    std::string preset;
    std::vector<Pop> pops;
    AnalysisResult analysis;
    ProcessInfo process;
    std::string toJson() const;
};

// Point d'entree unique, hors temps reel, en place :
// plops (guides par les cuts video, secondes) -> analyse de la prise brute -> nettoyage -> analyse de la voix
// propre -> dynamique (EQ, de-esser, nivellement, compresseur, loudness, limiteur).
// progress(etape) optionnel : "pops", "analysis", "processing".
Report runPipeline(AudioBuffer& b, const std::vector<double>& cutsSec, const std::string& presetJson,
                   void (*progress)(const char* stage, void* user) = nullptr, void* user = nullptr,
                   const Denoiser& externalDenoise = {});

} // namespace gv
