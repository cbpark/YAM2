/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "invisibles.h"

#include "input.h"      // InputKinematics
#include "momentum.h"   // FourMomentum
#include "variables.h"  // Variables

namespace yam2 {
Invisibles mkInvisibles(const InputKinematics &inp, const Variables &var,
                        double scale) {
    const double k1x = var.k1x();
    const double k1y = var.k1y();
    const double k1z = var.k1z();
    const double k2x = inp.ptmiss().px() - k1x;
    const double k2y = inp.ptmiss().py() - k1y;
    double k2z;
    if (inp.ptot_z() && var.dimension() == Variables::Dim::THREE) {
        k2z = inp.ptot_z().value_or(0.0) - inp.p1().pz() - inp.p2().pz() - k1z;
    } else {
        k2z = var.k2z();
    }

    const auto m{inp.minv()};

    FourMomentum k1{m, k1x, k1y, k1z};
    k1 *= scale;

    FourMomentum k2{m, k2x, k2y, k2z};
    k2 *= scale;

    return {k1, k2};
}
}  // namespace yam2
