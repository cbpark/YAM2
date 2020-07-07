#include "momentum.h"

#include <ostream>

namespace yam2 {
std::ostream &operator<<(std::ostream &os, const FourMomentum &p) {
    os << "e = " << p.t_ << ", px = " << p.x_ << ", py = " << p.y_
       << ", pz = " << p.z_;
    return os;
}
}  // namespace yam2
