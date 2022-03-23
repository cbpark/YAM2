/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "gradient.h"

#include <cmath>
#include <ostream>
#include <tuple>
#include "input.h"       // InputKinematics
#include "invisibles.h"  // Invisibles
#include "momentum.h"    // FourMomentum
#include "variables.h"   // Variables

using std::pair;

constexpr double EPS0 = 1.0e-20;

namespace yam2 {
double safeDivisor(double x) {
    if (x >= 0.0) { return std::fmax(EPS0, x); }
    return std::fmin(-EPS0, x);
}

pair<Gradient, double> m2Func1(const InputKinematics &, const FourMomentum &p1,
                               const Invisibles &ks) {
    const auto k1 = ks.k1();
    const double r1 = p1.e() / safeDivisor(k1.e());

    const double dk1x = r1 * k1.px() - p1.px();
    const double dk1y = r1 * k1.py() - p1.py();
    const double dk1z = r1 * k1.pz() - p1.pz();
    const double dk2z = 0.0;

    Gradient d1{dk1x, dk1y, dk1z, dk2z};

    const double m1 = invariantMass(p1, k1);
    d1 *= 1.0 / safeDivisor(m1);

    return {d1, m1};
}

pair<Gradient, double> m2Func2(const InputKinematics &inp,
                               const FourMomentum &p2, const Invisibles &ks) {
    const auto k1 = ks.k1(), k2 = ks.k2();
    const double r2 = p2.e() / safeDivisor(k2.e());

    const double dk1x = r2 * (k1.px() - inp.ptmiss().px()) + p2.px();
    const double dk1y = r2 * (k1.py() - inp.ptmiss().py()) + p2.py();
    double dk1z, dk2z;
    if (!inp.ptot_z()) {
        dk1z = 0.0;
        dk2z = r2 * k2.pz() - p2.pz();
    } else {
        dk1z = -r2 * k2.pz() + p2.pz();
        dk2z = 0.0;
    }

    Gradient d2{dk1x, dk1y, dk1z, dk2z};

    const double m2 = invariantMass(p2, k2);
    d2 *= 1.0 / safeDivisor(m2);

    return {d2, m2};
}

std::tuple<Gradients, double, double> m2Func(const InputKinematics &inp,
                                             const FourMomentum &p1,
                                             const FourMomentum &p2,
                                             const Invisibles &ks) {
    const auto &[d1, m1] = m2Func1(inp, p1, ks);
    const auto &[d2, m2] = m2Func2(inp, p2, ks);
    return {{d1, d2}, m1, m2};
}

Gradient mtotGrad(const InputKinematics &inp, const FourMomentum &p1,
                  const FourMomentum &p2, const Invisibles &ks, double m) {
    const auto p = p1 + p2;
    const double eV = p.e(), pz = p.pz();

    const auto k1 = ks.k1(), k2 = ks.k2();
    const double e1 = k1.e(), e2 = k2.e();
    const double k1z = k1.pz(), k2z = k2.pz();
    const double e1e2m = safeDivisor(e1 * e2 * m), e1m = safeDivisor(e1 * m),
                 e2m = safeDivisor(e2 * m);
    const double fac = -(eV + e1 + e2) / e1e2m;

    const double dk1x = fac * (e1 * inp.ptmiss().px() - (e1 + e2) * k1.px());
    const double dk1y = fac * (e1 * inp.ptmiss().py() - (e1 + e2) * k1.py());
    double dk1z, dk2z;
    if (!inp.ptot_z()) {
        dk1z = ((eV + e2) * k1z - (pz + k2z) * e1) / e1m;
        dk2z = ((eV + e1) * k2z - (pz + k1z) * e2) / e2m;
    } else {
        dk1z = fac * (e1 * (k1z + k2z) - (e1 + e2) * k1z);
        dk2z = 0.0;
    }

    return {dk1x, dk1y, dk1z, dk2z};
}

std::ostream &operator<<(std::ostream &os, const Gradient &g) {
    os << "dk1x = " << g.dk1x_ << ", dk1y = " << g.dk1y_
       << ", dk1z = " << g.dk1z_ << ", dk2z = " << g.dk2z_;
    return os;
}
}  // namespace yam2
