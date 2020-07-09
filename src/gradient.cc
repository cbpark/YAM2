#include "gradient.h"

#include "input.h"
#include "invisibles.h"
#include "momentum.h"
#include "variables.h"

#include <cmath>
#include <tuple>
#include <utility>

namespace yam2 {
double safeDivisor(double x) {
    const double eps0 = 1.0e-8;
    if (x >= 0.0) { return std::fmax(eps0, x); }
    return std::min(-eps0, x);
}

std::tuple<Gradients, double, double> m2Grad(const InputKinematics &inp,
                                             const Invisibles &ks,
                                             const FourMomentum &p1,
                                             const FourMomentum &p2,
                                             const Variables &var) {
    const auto k1 = ks.k1();
    const double m1 = invariantMass(p1, k1);
    const double m1inverse = 1.0 / safeDivisor(m1);
    const double r1 = p1.e() / safeDivisor(k1.e());
    auto d1 = Gradient(r1 * var.k1x() - p1.px(), r1 * var.k1y() - p1.py(),
                       r1 * var.k1z() - p1.pz(), 0);
    d1 *= m1inverse;

    const auto k2 = ks.k2();
    const double m2 = invariantMass(p2, k2);
    const double m2inverse = 1.0 / safeDivisor(m2);
    const double r2 = p2.e() / safeDivisor(k2.e());
    auto d2 = Gradient(r2 * (var.k1x() - inp.ptmiss().px()) + p2.px(),
                       r2 * (var.k1y() - inp.ptmiss().py()) + p2.py(), 0,
                       r2 * var.k2z() - p2.pz());
    d2 *= m2inverse;

    return {std::make_pair(d1, d2), m1, m2};
}
}  // namespace yam2
