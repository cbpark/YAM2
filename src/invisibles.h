/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#ifndef YAM2_SRC_INVISIBLES_H_
#define YAM2_SRC_INVISIBLES_H_

#include "input.h"      // InputKinematics
#include "momentum.h"   // FourMomentum
#include "variables.h"  // Variables

namespace yam2 {
class Invisibles {
private:
    FourMomentum k1_, k2_;

    Invisibles(const FourMomentum &k1, const FourMomentum &k2)
        : k1_(k1), k2_(k2) {}

public:
    Invisibles() = delete;

    FourMomentum k1() const { return k1_; }
    FourMomentum k2() const { return k2_; }

    friend Invisibles mkInvisibles(const InputKinematics &inp,
                                   const Variables &var, double scale);
};

inline Invisibles mkInvisibles(const InputKinematics &inp, const Variables &var,
                               double scale = 1.0) {
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

#endif  // YAM2_SRC_INVISIBLES_H_
