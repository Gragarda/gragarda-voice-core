#pragma once
#include <cstddef>
#include <vector>

namespace gv {

// Audio planaire float -1..1 : ch[c][i]. Tout le fichier en memoire (traitement hors temps reel).
struct AudioBuffer {
    int sampleRate = 48000;
    std::vector<std::vector<float>> ch;

    size_t frames() const { return ch.empty() ? 0 : ch[0].size(); }
    int channels() const { return (int)ch.size(); }
};

} // namespace gv
