#ifndef YAM2_SRC_YAM2_H_
#define YAM2_SRC_YAM2_H_

#include "invisibles.h"
#include "momentum.h"

namespace yam2 {
class M2Solution {
private:
    double m2_;
    Invisibles ksol_;

public:
    M2Solution() = delete;

    double m2() const { return m2_; }
    FourMomentum k1() const { return ksol_.k1(); }
    FourMomentum k2() const { return ksol_.k2(); }
};
}  // namespace yam2

#endif  // namespace yam2
