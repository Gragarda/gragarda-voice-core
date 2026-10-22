#include "gv/Pipeline.h"
#include <cstdio>
#include "gv/PopRepair.h"

namespace gv {

Report runPipeline(AudioBuffer& b, const std::vector<double>& cutsSec, const std::string& presetJson,
                   void (*progress)(const char*, void*), void* user, const Denoiser& externalDenoise)
{
    auto step = [&](const char* s) { if (progress) progress(s, user); };
    Report r;
    r.sampleRate = b.sampleRate;
    const Preset p = parsePreset(presetJson);
    r.preset = p.name;
    step("pops");
    r.pops = detectPops(b, cutsSec);
    repairPops(b, r.pops); // AVANT toute dynamique
    step("analysis");
    r.analysis = analyze(b);
    step("processing");
    r.process = processChain(b, p, r.analysis, externalDenoise);
    return r;
}

std::string Report::toJson() const
{
    std::string s = "{\"preset\":\"" + preset + "\",\"popsRemoved\":" + std::to_string(pops.size()) + ",\"pops\":[";
    char buf[160];
    for (size_t i = 0; i < pops.size(); ++i) {
        const Pop& p = pops[i];
        std::snprintf(buf, sizeof(buf), "%s{\"time\":%.4f,\"type\":\"%s\",\"origin\":\"%s\",\"severity\":%.1f}",
                      i ? "," : "", (double)p.pos / sampleRate, toString(p.type),
                      p.origin == PopOrigin::Cut ? "cut" : "scan", p.severity);
        s += buf;
    }
    return s + "],\"analysis\":" + analysis.toJson() + ",\"process\":" + process.toJson() + "}";
}

} // namespace gv
