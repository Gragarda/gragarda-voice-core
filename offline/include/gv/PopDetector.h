#pragma once
#include <vector>
#include "gv/AudioBuffer.h"

namespace gv {

enum class PopType { Click, DcStep, Both, Silence };
enum class PopOrigin { Cut, Scan };

struct Pop {
    size_t pos = 0;      // echantillon du point de coupe
    PopType type = PopType::Click;
    float severity = 0;  // saut / reference locale (ou saut DC en dB pour DcStep seul)
    PopOrigin origin = PopOrigin::Scan;
};

// Cuts video (secondes) = zones de recherche sensibles (x8) ; ailleurs scan strict (x20).
// cutsSec vide = mode "scan seul".
std::vector<Pop> detectPops(const AudioBuffer& b, const std::vector<double>& cutsSec);

// Erreur de prediction lineaire j(i) = |x[i] - (2x[i-1] - x[i-2])|, max sur les canaux
float jumpAt(const AudioBuffer& b, size_t i);
// Mediane de j sur 10 ms de chaque cote de p, en excluant +-2 ms autour de p
float localJumpRef(const AudioBuffer& b, size_t p);

// Saut de niveau moyen a p (apres - avant), moyennes ponderees Hann sur w echantillons de chaque cote
float dcDelta(const std::vector<float>& x, size_t p, size_t w);

const char* toString(PopType t);

} // namespace gv
