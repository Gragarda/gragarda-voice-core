#pragma once
#include <vector>
#include "gv/PopDetector.h"

namespace gv {

// Repare en place. Doit passer AVANT tout traitement dynamique (un compresseur grossit les plops).
// Garde-fou : si l'interpolation laisse un saut > 3x la reference, on passe au fondu croise
// entre prolongements LPC (ordre 32).
void repairPops(AudioBuffer& b, const std::vector<Pop>& pops);

} // namespace gv
