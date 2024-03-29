/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "YAM2/input.h"

#include <nlopt.hpp>  // for the NLopt API

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <optional>
#include <vector>

#include "YAM2/invisibles.h"  // mkInvisibles
#include "YAM2/momentum.h"    // FourMomentum, TransverseMomentum, Mass
#include "YAM2/variables.h"   // NLoptVar, mkVariables
#include "gradient.h"         // mTotGrad

using std::vector;

namespace yam2 {
NLoptVar InputKinematics::initial_guess(double eps, unsigned int neval) {
    // nlopt::opt algorithm{nlopt::LN_NELDERMEAD, 4};
    nlopt::opt algorithm{nlopt::LD_SLSQP, 4};
    // nlopt::opt algorithm{nlopt::LD_TNEWTON, 4};
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

void InputKinematics::show(std::ostream &os) const {
    os << "p1: " << p1_ << '\n'
       << "p2: " << p2_ << '\n'
       << "q1: " << q1_ << '\n'
       << "q2: " << q2_ << '\n'
       << "ptmiss: " << ptmiss_ << '\n'
       << "M(invisible1): " << minv1_.value << '\n'
       << "M(invisible2): " << minv2_.value << '\n';

    if (mparent1_) { os << "M(parent1): " << mparent1_.value().value << '\n'; }
    if (mparent2_) { os << "M(parent2): " << mparent2_.value().value << '\n'; }
    if (mrel1_) { os << "M(relative1): " << mrel1_.value().value << '\n'; }
    if (mrel2_) { os << "M(relative2): " << mrel2_.value().value << '\n'; }
    if (sqrt_s_ > 0.0) { os << "sqrt(s): " << sqrt_s_ << '\n'; }
    if (ptot_z_) { os << "Pz: " << ptot_z_.value() << '\n'; }

    os << "scale: " << scale_;
}

std::optional<InputKinematics> mkInput(
    const vector<FourMomentum> &as, const vector<FourMomentum> &bs,
    const TransverseMomentum &ptmiss, const Mass &minv1, const Mass &minv2,
    const std::optional<Mass> &mparent1, const std::optional<Mass> &mparent2,
    const std::optional<Mass> &mrel1, const std::optional<Mass> &mrel2,
    double sqrt_s, const std::optional<double> ptot_z) {
    if (as.size() != 2 || bs.size() != 2) {
        std::cerr << "mkInput: Invalid number of visible particles.\n";
        return {};
    }

    if (as[0].msq() < -1.0e-10 || as[1].msq() < -1.0e-10 ||
        bs[0].msq() < -1.0e-10 || bs[1].msq() < -1.0e-10) {
        std::cerr << "mkInput: Unphysical input.\n";
        return {};
    }

    vector<FourMomentum> ps;
    std::transform(as.cbegin(), as.cend(), bs.cbegin(), std::back_inserter(ps),
                   [](const auto &a, const auto &b) { return a + b; });

    const auto p1 = ps.front(), p2 = ps.back();
    const auto q1 = bs.front(), q2 = bs.back();

    double sval;
    if (mparent1 && mparent2) {
        sval = 0.5 * (mparent1.value().value + mparent2.value().value);
    } else {
        const double e1 = p1.e(), e2 = p2.e();
        const double etmiss_sq =
            ptmiss.ptsq() + minv1.square() + minv2.square();
        const double scalesq = e1 * e1 + e2 * e2 + etmiss_sq;

        if (scalesq <= 0.0) {
            std::cerr << "mkInput: unphysical energy scale.\n";
            return {};
        }
        sval = 8.0 * std::sqrt(scalesq);  // scale > 0
    }
    if (sval <= 0.0) {
        std::cerr << "mkInput: found zero or negative energy scale.\n";
        return {};
    }

    return {{p1 / sval, p2 / sval, q1 / sval, q2 / sval, ptmiss / sval,
             minv1 / sval, minv2 / sval, scaleIfExists(mparent1, sval),
             scaleIfExists(mparent2, sval), scaleIfExists(mrel1, sval),
             scaleIfExists(mrel2, sval), sqrt_s / sval,
             scaleIfExists(ptot_z, sval), sval}};
}

std::ostream &operator<<(std::ostream &os, const InputKinematics &p) {
    p.show(os);
    return os;
}

void InputKinematicsWithVertex::show(std::ostream &os) const {
    InputKinematics::show(os);
    os << "\nvertex1: " << vertex1_ << '\n';
    os << "vertex2: " << vertex2_ << '\n';
    os << "delta_theta: " << delta_theta_;
}

std::optional<InputKinematicsWithVertex> mkInputWithVertex(
    const std::optional<InputKinematics> &input_kinematics,
    const SpatialMomentum &vertex1, const SpatialMomentum &vertex2,
    double delta_theta) {
    if (!input_kinematics) {
        std::cerr << "mkInputWithVertex: invalid input kinematics.\n";
        return {};
    }

    if (delta_theta < 0.0) {
        std::cerr << "mkInputWithVertex: delta_theta must be positive.\n";
        return {};
    }

    const auto v1 = vertex1.normalize();
    const auto v2 = vertex2.normalize();

    return {{input_kinematics.value().p1(),
             input_kinematics.value().p2(),
             input_kinematics.value().q1(),
             input_kinematics.value().q2(),
             input_kinematics.value().ptmiss(),
             input_kinematics.value().minv1(),
             input_kinematics.value().minv2(),
             v1,
             v2,
             delta_theta,
             {input_kinematics.value().mparent1()},
             {input_kinematics.value().mparent2()},
             {input_kinematics.value().mrel1()},
             {input_kinematics.value().mrel2()},
             input_kinematics.value().sqrt_s(),
             input_kinematics.value().ptot_z(),
             input_kinematics.value().scale()}};
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
        const auto grad_mtot = mtotGrad(*inp, p1, p2, ks, m_tot);
        grad_mtot.set_gradient(grad);
    }

    return m_tot - inp->sqrt_s();
}
}  // namespace yam2
