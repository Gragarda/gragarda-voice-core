#pragma once
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include "gv/AudioBuffer.h"

// Lecture / ecriture WAV (PCM 16/24/32 bits, float 32) pour gv-cli et les tests. Pas utilise par le web.
namespace gv {

inline bool readWav(const std::string& path, AudioBuffer& b)
{
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return false;
    char id[4];
    uint32_t sz;
    if (std::fread(id, 1, 4, f) != 4 || std::memcmp(id, "RIFF", 4) || std::fread(&sz, 4, 1, f) != 1 ||
        std::fread(id, 1, 4, f) != 4 || std::memcmp(id, "WAVE", 4)) { std::fclose(f); return false; }
    uint16_t fmt = 0, nch = 0, bits = 0;
    uint32_t sr = 0;
    bool ok = false;
    while (std::fread(id, 1, 4, f) == 4 && std::fread(&sz, 4, 1, f) == 1) {
        if (!std::memcmp(id, "fmt ", 4)) {
            std::vector<uint8_t> c(sz);
            if (std::fread(c.data(), 1, sz, f) != sz || sz < 16) break;
            std::memcpy(&fmt, &c[0], 2); std::memcpy(&nch, &c[2], 2);
            std::memcpy(&sr, &c[4], 4); std::memcpy(&bits, &c[14], 2);
            if (fmt == 0xFFFE && sz >= 26) std::memcpy(&fmt, &c[24], 2); // WAVE_FORMAT_EXTENSIBLE
        } else if (!std::memcmp(id, "data", 4) && nch) {
            const size_t bps = bits / 8, n = sz / (bps * nch);
            std::vector<uint8_t> d(n * bps * nch);
            if (std::fread(d.data(), 1, d.size(), f) != d.size()) break;
            b.sampleRate = (int)sr;
            b.ch.assign(nch, std::vector<float>(n));
            for (size_t i = 0; i < n; ++i)
                for (int c = 0; c < nch; ++c) {
                    const uint8_t* p = &d[(i * nch + c) * bps];
                    float v = 0;
                    if (fmt == 3 && bits == 32) std::memcpy(&v, p, 4);
                    else if (bits == 16) v = (int16_t)(p[0] | p[1] << 8) / 32768.0f;
                    else if (bits == 24) v = (int32_t)((uint32_t)p[0] << 8 | (uint32_t)p[1] << 16 | (uint32_t)p[2] << 24) / 2147483648.0f;
                    else if (bits == 32) { int32_t s; std::memcpy(&s, p, 4); v = s / 2147483648.0f; }
                    b.ch[c][i] = v;
                }
            ok = true;
            break;
        } else {
            std::fseek(f, (long)(sz + (sz & 1)), SEEK_CUR);
        }
    }
    std::fclose(f);
    return ok;
}

inline bool writeWav(const std::string& path, const AudioBuffer& b) // float 32
{
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    const uint16_t nch = (uint16_t)b.channels(), fmt = 3, bits = 32, align = nch * 4;
    const uint32_t sr = b.sampleRate, rate = sr * align, data = (uint32_t)(b.frames() * align), riff = 36 + data, fsz = 16;
    std::fwrite("RIFF", 1, 4, f); std::fwrite(&riff, 4, 1, f); std::fwrite("WAVEfmt ", 1, 8, f);
    std::fwrite(&fsz, 4, 1, f); std::fwrite(&fmt, 2, 1, f); std::fwrite(&nch, 2, 1, f);
    std::fwrite(&sr, 4, 1, f); std::fwrite(&rate, 4, 1, f); std::fwrite(&align, 2, 1, f);
    std::fwrite(&bits, 2, 1, f); std::fwrite("data", 1, 4, f); std::fwrite(&data, 4, 1, f);
    for (size_t i = 0; i < b.frames(); ++i)
        for (const auto& c : b.ch) std::fwrite(&c[i], 4, 1, f);
    std::fclose(f);
    return true;
}

} // namespace gv
