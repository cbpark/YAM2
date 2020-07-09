/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "input.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <optional>
#include <vector>

#include "momentum.h"  // FourMomentum, TransverseMomentum, Mass

using std::vector;

namespace yam2 {
double getScale(const FourMomentum &p) { return p.pt(); }
double getScale(const TransverseMomentum &p) { return p.pt(); }

std::optional<InputKinematics> mkInput(const vector<FourMomentum> &as,
                                       const vector<FourMomentum> &bs,
                                       const TransverseMomentum &ptmiss,
                                       const Mass &minv) {
    if (as.size() != 2 || bs.size() != 2) { return {}; }

    vector<FourMomentum> ps;
    std::transform(as.cbegin(), as.cend(), bs.cbegin(), std::back_inserter(ps),
                   [](const auto &a, const auto &b) { return a + b; });

    const auto p1 = ps.front(), p2 = ps.back();
    const auto q1 = bs.front(), q2 = bs.back();
    const double m1 = p1.m(), m2 = p2.m();
    const double scalesq = getScale(p1) + getScale(p2) + getScale(ptmiss) +
                           m1 * m1 + m2 * m2 + 2 * minv.square();

    if (scalesq <= 0) { return {}; }

    const double scale = std::sqrt(scalesq);
    const double s = 1.0 / scale;

    return {{p1 * s, p2 * s, q1 * s, q2 * s, ptmiss * s, minv * s, scale}};
}

std::ostream &operator<<(std::ostream &os, const InputKinematics &p) {
    os << "p1: " << p.p1_ << '\n'
       << "p2: " << p.p2_ << '\n'
       << "q1: " << p.q1_ << '\n'
       << "q2: " << p.q2_ << '\n'
       << "ptmiss: " << p.ptmiss_ << '\n'
       << "m(invisible): " << p.minv_.value;
    return os;
}
}  // namespace yam2
