#ifndef YAM2_SRC_INPUT_H_
#define YAM2_SRC_INPUT_H_

#include "momentum.h"

#include <optional>
#include <vector>

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
    /** squared mass of the invisible particle */
    double minvsq_;
    double scale_;

    InputKinematics(const FourMomentum &p1, const FourMomentum &p2,
                    const FourMomentum &q1, const FourMomentum &q2,
                    const TransverseMomentum &ptmiss, double minvsq,
                    double scale)
        : p1_(p1),
          p2_(p2),
          q1_(q1),
          q2_(q2),
          ptmiss_(ptmiss),
          minvsq_(minvsq),
          scale_(scale) {}

public:
    InputKinematics() = delete;

    FourMomentum p1() const { return p1_; }
    FourMomentum p2() const { return p2_; }
    FourMomentum q1() const { return q1_; }
    FourMomentum q2() const { return q2_; }
    TransverseMomentum ptmiss() const { return ptmiss_; }
    double minvsq() const { return minvsq_; }
    double scale() const { return scale_; }

    friend std::optional<InputKinematics> mkInput(
        const std::vector<FourMomentum> &as,
        const std::vector<FourMomentum> &bs, const TransverseMomentum &ptmiss,
        double minv);
};

std::optional<InputKinematics> mkInput(const std::vector<FourMomentum> &as,
                                       const std::vector<FourMomentum> &bs,
                                       const TransverseMomentum &ptmiss,
                                       double minv);
}  // namespace yam2

#endif  // YAM2_SRC_INPUT_H_
