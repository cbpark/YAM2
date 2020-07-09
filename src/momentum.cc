#include "momentum.h"

#include <ostream>

namespace yam2 {
std::ostream& operator<<(std::ostream& os, const FourMomentum& p) {
    os << "e = " << p.t_ << ", px = " << p.x_ << ", py = " << p.y_
       << ", pz = " << p.z_;
    return os;
}

std::ostream& operator<<(std::ostream& os, const TransverseMomentum& p) {
    os << "px = " << p.x_ << ", py = " << p.y_;
    return os;
}

double invariantMass(const FourMomentum& p, const FourMomentum& k) {
    const auto psum = p + k;
    return psum.m();
}
}  // namespace yam2
