/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#ifndef YAM2_SRC_INPUT_H_
#define YAM2_SRC_INPUT_H_

#include <optional>
#include <ostream>
#include <vector>

#include "momentum.h"

namespace yam2 {
class InputKinematics {
private:
    /** p1 = a1 + b1 */
    FourMomentum p1_;
    /** p2 = a2 + b2 */
    FourMomentum p2_;
    /** q1 = b1 */
    FourMomentum q1_;
    /** q2 = b2 */
    FourMomentum q2_;
    TransverseMomentum ptmiss_;
    /** the input mass for the invisible particle */
    Mass minv_;
    /** the input mass for the relative particle */
    Mass mrel_;
    double scale_;

    InputKinematics(const FourMomentum &p1, const FourMomentum &p2,
                    const FourMomentum &q1, const FourMomentum &q2,
                    const TransverseMomentum &ptmiss, const Mass &minv,
                    const Mass &mrel, double scale)
        : p1_(p1),
          p2_(p2),
          q1_(q1),
          q2_(q2),
          ptmiss_(ptmiss),
          minv_(minv),
          mrel_(mrel),
          scale_(scale) {}

    InputKinematics(const FourMomentum &p1, const FourMomentum &p2,
                    const FourMomentum &q1, const FourMomentum &q2,
                    const TransverseMomentum &ptmiss, const Mass &minv,
                    double scale)
        : InputKinematics(p1, p2, q1, q2, ptmiss, minv, Mass{0}, scale) {}

public:
    InputKinematics() = delete;

    FourMomentum p1() const { return p1_; }
    FourMomentum p2() const { return p2_; }
    FourMomentum q1() const { return q1_; }
    FourMomentum q2() const { return q2_; }

    TransverseMomentum ptmiss() const { return ptmiss_; }

    Mass minv() const { return minv_; }
    double minv_square() const { return minv_.square(); }

    Mass mrel() const { return mrel_; }
    double mrel_square() const { return mrel_.square(); }

    double scale() const { return scale_; }

    friend std::optional<InputKinematics> mkInput(
        const std::vector<FourMomentum> &as,
        const std::vector<FourMomentum> &bs, const TransverseMomentum &ptmiss,
        const Mass &minv, const std::optional<Mass> &mrel);

    friend std::ostream &operator<<(std::ostream &os, const InputKinematics &p);
};

std::optional<InputKinematics> mkInput(const std::vector<FourMomentum> &as,
                                       const std::vector<FourMomentum> &bs,
                                       const TransverseMomentum &ptmiss,
                                       const Mass &minv,
                                       const std::optional<Mass> &mrel);
}  // namespace yam2

#endif  // YAM2_SRC_INPUT_H_
