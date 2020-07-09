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
    FourMomentum k1{inp.minv(), var.k1x(), var.k1y(), var.k1z()};
    k1 *= scale;

    const double k2x = inp.ptmiss().px() - var.k1x();
    const double k2y = inp.ptmiss().py() - var.k1y();
    FourMomentum k2{inp.minv(), k2x, k2y, var.k2z()};
    k2 *= scale;

    return {k1, k2};
}
}  // namespace yam2

#endif  // YAM2_SRC_INVISIBLES_H_
