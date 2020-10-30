/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "input.h"

#include <nlopt.hpp>  // for the NLopt API

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <optional>
#include <vector>

#include "gradient.h"    // mTotGrad
#include "invisibles.h"  // mkInvisibles
#include "momentum.h"    // FourMomentum, TransverseMomentum, Mass
#include "variables.h"   // NLoptVar, mkVariables

using std::vector;

namespace yam2 {
NLoptVar InputKinematics::initial_guess(double eps, unsigned int neval) {
    // nlopt::opt algorithm{nlopt::LN_NELDERMEAD, 4};
    // nlopt::opt algorithm{nlopt::LD_SLSQP, 4};
    nlopt::opt algorithm{nlopt::LD_TNEWTON, 4};
    algorithm.set_min_objective(deltaSqrtS, this);
    const double epsf = eps * 0.1;
    // algorithm.set_ftol_rel(epsf);
    algorithm.set_ftol_abs(epsf);
    algorithm.set_maxeval(neval);

    const NLoptVar x0{0.5 * ptmiss_.px(), 0.5 * ptmiss_.py(), 0.0, 0.0};
    auto x{x0};
    double minf;
    auto result = algorithm.optimize(x, minf);
    if (result < 0) { return x0; }
    return x;
}

std::optional<InputKinematics> mkInput(const vector<FourMomentum> &as,
                                       const vector<FourMomentum> &bs,
                                       const TransverseMomentum &ptmiss,
                                       const Mass &minv, const double sqrt_s,
                                       const std::optional<Mass> &mrel) {
    if (as.size() != 2 || bs.size() != 2) {
        std::cerr << "mkInput: Invalid number of visible particles.\n";
        return {};
    }

    vector<FourMomentum> ps;
    std::transform(as.cbegin(), as.cend(), bs.cbegin(), std::back_inserter(ps),
                   [](const auto &a, const auto &b) { return a + b; });

    const auto p1 = ps.front(), p2 = ps.back();
    const auto q1 = bs.front(), q2 = bs.back();
    const double e1 = p1.e(), e2 = p2.e();
    const double etmiss_sq = ptmiss.ptsq() + 2.0 * minv.square();
    const double scalesq = e1 * e1 + e2 * e2 + etmiss_sq;

    if (scalesq <= 0.0) {
        std::cerr << "mkInput: found zero or negative energy scale.\n";
        return {};
    }

    const double scale = 8.0 * std::sqrt(scalesq);  // scale > 0

    if (!mrel) {
        return {{p1 / scale, p2 / scale, q1 / scale, q2 / scale, ptmiss / scale,
                 minv / scale, sqrt_s / scale, scale}};
    }
    return {{p1 / scale, p2 / scale, q1 / scale, q2 / scale, ptmiss / scale,
             minv / scale, mrel.value() / scale, sqrt_s / scale, scale}};
}

std::ostream &operator<<(std::ostream &os, const InputKinematics &p) {
    os << "p1: " << p.p1_ << '\n'
       << "p2: " << p.p2_ << '\n'
       << "q1: " << p.q1_ << '\n'
       << "q2: " << p.q2_ << '\n'
       << "ptmiss: " << p.ptmiss_ << '\n'
       << "m(invisible): " << p.minv_.value << '\n'
       << "sqrt(s): " << p.sqrt_s_;
    return os;
}

double deltaSqrtS(const NLoptVar &x, NLoptVar &grad, void *input) {
    const auto var = mkVariables(x);  // this will not be empty.
    const auto var_val = var.value();
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    const auto ks = mkInvisibles(*inp, var_val);

    const auto p1 = inp->p1(), p2 = inp->p2();
    const auto ptot = p1 + p2 + ks.k1() + ks.k2();
    const double m_tot = ptot.m();
    if (!grad.empty()) {
        const auto grad_mtot = mtotGrad(*inp, p1, p2, ks, var_val, m_tot);
        grad = grad_mtot.gradient();
    }
    return m_tot - inp->sqrt_s();
}
}  // namespace yam2
