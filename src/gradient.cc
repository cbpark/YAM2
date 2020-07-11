/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "gradient.h"

#include <cmath>
#include <tuple>

#include "input.h"       // InputKinematics
#include "invisibles.h"  // Invisibles
#include "momentum.h"    // FourMomentum
#include "variables.h"   // Variables

namespace yam2 {
double safeDivisor(double x) {
    const double eps0 = 1.0e-8;
    if (x >= 0.0) { return std::fmax(eps0, x); }
    return std::fmin(-eps0, x);
}

std::tuple<Gradients, double, double> m2Grad(const InputKinematics &inp,
                                             const FourMomentum &p1,
                                             const FourMomentum &p2,
                                             const Invisibles &ks,
                                             const Variables &var) {
    const auto k1 = ks.k1();
    const double m1 = invariantMass(p1, k1);
    const double m1inverse = 1.0 / safeDivisor(m1);
    const double r1 = p1.e() / safeDivisor(k1.e());
    Gradient d1{r1 * var.k1x() - p1.px(), r1 * var.k1y() - p1.py(),
                r1 * var.k1z() - p1.pz(), 0};
    d1 *= m1inverse;

    const auto k2 = ks.k2();
    const double m2 = invariantMass(p2, k2);
    const double m2inverse = 1.0 / safeDivisor(m2);
    const double r2 = p2.e() / safeDivisor(k2.e());
    Gradient d2{r2 * (var.k1x() - inp.ptmiss().px()) + p2.px(),
                r2 * (var.k1y() - inp.ptmiss().py()) + p2.py(), 0,
                r2 * var.k2z() - p2.pz()};
    d2 *= m2inverse;

    return {std::make_pair(d1, d2), m1, m2};
}

Gradient mtotGrad(const InputKinematics &inp, const FourMomentum &p1,
                  const FourMomentum &p2, const Invisibles &ks,
                  const Variables &var, double m) {
    const auto p = p1 + p2;
    const double eV = p.e(), pz = p.pz();
    const double e1 = ks.k1().e(), e2 = ks.k2().e();
    const double e1e2m = safeDivisor(e1 * e2 * m), e1m = safeDivisor(e1 * m),
                 e2m = safeDivisor(e2 * m);
    const double fac = -(eV + e1 + e2) / e1e2m;
    const double k1z = var.k1z(), k2z = var.k2z();

    const double dk1x = fac * (e1 * inp.ptmiss().px() - (e1 + e2) * var.k1x());
    const double dk1y = fac * (e1 * inp.ptmiss().py() - (e1 + e2) * var.k1y());
    const double dk1z = ((eV + e2) * k1z - (pz + k2z) * e1) / e1m;
    const double dk2z = ((eV + e1) * k2z - (pz + k1z) * e2) / e2m;

    return {dk1x, dk1y, dk1z, dk2z};
}
}  // namespace yam2
