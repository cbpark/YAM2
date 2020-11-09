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

Invisibles mkInvisibles(const InputKinematics &inp, const Variables &var,
                        double scale = 1.0);
}  // namespace yam2

#endif  // YAM2_SRC_INVISIBLES_H_
