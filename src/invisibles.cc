#include "invisibles.h"
#include "input.h"
#include "momentum.h"
#include "variables.h"

namespace yam2 {
Invisibles mkInvisibles(const InputKinematics &inp, const Variables &var) {
    const auto k1 = FourMomentum(inp.minv(), var.k1x(), var.k1y(), var.k1z());

    const double k2x = inp.ptmiss().px() - var.k1x();
    const double k2y = inp.ptmiss().py() - var.k1y();
    const auto k2 = FourMomentum(inp.minv(), k2x, k2y, var.k2z());

    return {k1, k2};
}
}  // namespace yam2
