/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "YAM2/constraint.h"

#include <cmath>  // std::pow, std::tan, std::acos
#include <optional>
#include <utility>            // std::pair
#include "YAM2/input.h"       // InputKinematics, deltaSqrtS
#include "YAM2/invisibles.h"  // mkInvisibles
#include "YAM2/momentum.h"    // SpatialMomentum
#include "YAM2/variables.h"   // mkVariables, NLoptVar
#include "gradient.h"         // m2Func, m2Func1, m2Func2

namespace yam2 {
enum class DecayChainOf { Parent1, Parent2 };

double constraintEq(const InputKinematics &inp, const FourMomentum &p1,
                    const FourMomentum &p2, const NLoptVar &x, NLoptVar &grad) {
    const auto var = mkVariables(x);
    const auto var_val = var.value();
    const auto ks = mkInvisibles(inp, var_val);

    const auto &[grads, m1, m2] = m2Func(inp, p1, p2, ks);
    const auto &[grad1, grad2] = grads;

    if (!grad.empty()) {
        const auto dgrad = grad1 - grad2;
        dgrad.set_gradient(grad);
    }
    return m1 - m2;
}

double constraintA(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraintEq(*inp, inp->p1(), inp->p2(), x, grad);
}

double constraintB(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraintEq(*inp, inp->q1(), inp->q2(), x, grad);
}

double constraintAP(const InputKinematics &inp, const FourMomentum &p,
                    const std::optional<Mass> &mparent, GradFunc fmGrad,
                    const NLoptVar &x, NLoptVar &grad) {
    const auto var = mkVariables(x);
    const auto var_val = var.value();
    const auto ks = mkInvisibles(inp, var_val);

    const auto &[grad_, m_] = fmGrad(inp, p, ks);
    if (!grad.empty()) { grad_.set_gradient(grad); }
    return m_ - mparent.value_or(Mass{m_}).value;
}

double constraintAPUpper(const InputKinematics &inp, const FourMomentum &p,
                         const std::optional<Mass> &mparent, GradFunc fmGrad,
                         const NLoptVar &x, NLoptVar &grad) {
    const auto var = mkVariables(x);
    const auto var_val = var.value();
    const auto ks = mkInvisibles(inp, var_val);

    const auto &[grad_, m_] = fmGrad(inp, p, ks);
    if (!grad.empty()) { grad_.set_gradient(grad); }
    return m_ - mparent.value_or(Mass{m_}).value * (1.0 + inp.eps_constraint());
}

double constraintAPLower(const InputKinematics &inp, const FourMomentum &p,
                         const std::optional<Mass> &mparent, GradFunc fmGrad,
                         const NLoptVar &x, NLoptVar &grad) {
    const auto var = mkVariables(x);
    const auto var_val = var.value();
    const auto ks = mkInvisibles(inp, var_val);

    auto [grad_, m_] = fmGrad(inp, p, ks);
    grad_ = -grad_;
    if (!grad.empty()) { grad_.set_gradient(grad); }
    return -m_ +
           mparent.value_or(Mass{m_}).value * (1.0 - inp.eps_constraint());
}

double constraintA1(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraintAP(*inp, inp->p1(), inp->mparent1(), m2Func1, x, grad);
}

double constraintA1Upper(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraintAPUpper(*inp, inp->p1(), inp->mparent1(), m2Func1, x,
                             grad);
}

double constraintA1Lower(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraintAPLower(*inp, inp->p1(), inp->mparent1(), m2Func1, x,
                             grad);
}

double constraintA2(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraintAP(*inp, inp->p2(), inp->mparent2(), m2Func2, x, grad);
}

double constraintA2Upper(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraintAPUpper(*inp, inp->p2(), inp->mparent2(), m2Func2, x,
                             grad);
}

double constraintA2Lower(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraintAPLower(*inp, inp->p2(), inp->mparent2(), m2Func2, x,
                             grad);
}

double constraintR(const InputKinematics &inp, const FourMomentum &q,
                   const std::optional<Mass> &mrel, GradFunc fmGrad,
                   const NLoptVar &x, NLoptVar &grad) {
    const auto var = mkVariables(x);
    const auto var_val = var.value();
    const auto ks = mkInvisibles(inp, var_val);
    const auto &[grad_, m_] = fmGrad(inp, q, ks);

    if (!grad.empty()) { grad_.set_gradient(grad); }
    return m_ - mrel.value_or(Mass{m_}).value;
}

double constraintR1(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraintR(*inp, inp->q1(), inp->mrel1(), m2Func1, x, grad);
}

double constraintR2(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraintR(*inp, inp->q2(), inp->mrel2(), m2Func2, x, grad);
}

double constraintSqrtS(const NLoptVar &x, NLoptVar &grad, void *input) {
    return deltaSqrtS(x, grad, input);
}

SpatialMomentum getParent(const NLoptVar &x, const InputKinematics &inp,
                          const DecayChainOf &decay_chain) {
    const auto var = mkVariables(x);
    const auto var_val = var.value();
    const auto ks = mkInvisibles(inp, var_val);
    if (decay_chain == DecayChainOf::Parent1) {
        return (inp.p1() + ks.k1()).three_momentum();
    } else {
        return (inp.p2() + ks.k2()).three_momentum();
    }
}

double constraintVertex1X(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematicsWithVertex *>(input);
    auto parent1 = getParent(x, *inp, DecayChainOf::Parent1);
    auto parent1_norm_vec = parent1.normalize();

    if (!grad.empty()) {
        double parent1_norm = parent1.norm();

        double dk1x =
            (1 - parent1_norm_vec.px() * parent1_norm_vec.px()) / parent1_norm;
        double dk1y =
            -parent1_norm_vec.px() * parent1_norm_vec.py() / parent1_norm;
        double dk1z =
            -parent1_norm_vec.px() * parent1_norm_vec.pz() / parent1_norm;
        double dk2z = 0.0;

        auto grad_ = Gradient(dk1x, dk1y, dk1z, dk2z);
        grad_.set_gradient(grad);
    }

    return parent1_norm_vec.px() - inp->vertex1().px();
}

double constraintVertex1Y(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematicsWithVertex *>(input);
    auto parent1 = getParent(x, *inp, DecayChainOf::Parent1);
    auto parent1_norm_vec = parent1.normalize();

    if (!grad.empty()) {
        double parent1_norm = parent1.norm();

        double dk1x =
            -parent1_norm_vec.px() * parent1_norm_vec.py() / parent1_norm;
        double dk1y =
            (1 - parent1_norm_vec.py() * parent1_norm_vec.py()) / parent1_norm;
        double dk1z =
            -parent1_norm_vec.py() * parent1_norm_vec.pz() / parent1_norm;
        double dk2z = 0.0;

        auto grad_ = Gradient(dk1x, dk1y, dk1z, dk2z);
        grad_.set_gradient(grad);
    }

    return parent1_norm_vec.py() - inp->vertex1().py();
}

double constraintVertex1Theta(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematicsWithVertex *>(input);
    auto parent1 = getParent(x, *inp, DecayChainOf::Parent1);
    double theta = parent1.theta();

    if (!grad.empty()) {
        double pt = parent1.pt();
        double p2 = parent1.p2();

        double dk1x = parent1.px() * parent1.pz() / pt / p2;
        double dk1y = parent1.py() * parent1.pz() / pt / p2;
        double dk1z = -pt / p2;
        double dk2z = 0.0;

        // double r = std::pow(std::tan(theta), 2);
        // double denom = (1.0 + r) * parent1.pz();

        // double dk1x = parent1.px() / (pt * denom);
        // double dk1y = parent1.py() / (pt * denom);
        // double dk1z = -pt / (denom * parent1.pz());
        // double dk2z = 0.0;

        auto grad_ = Gradient(dk1x, dk1y, dk1z, dk2z);
        grad_.set_gradient(grad);
    }

    return theta - inp->vertex1().theta();
}

double constraintVertex1Phi(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematicsWithVertex *>(input);
    auto parent1 = getParent(x, *inp, DecayChainOf::Parent1);

    if (!grad.empty()) {
        double pt2 = parent1.pt2();
        double dk1x = -parent1.py() / pt2;
        double dk1y = parent1.px() / pt2;
        double dk1z = 0.0;
        double dk2z = 0.0;

        auto grad_ = Gradient(dk1x, dk1y, dk1z, dk2z);
        grad_.set_gradient(grad);
    }

    return parent1.phi() - inp->vertex1().phi();
}

double acos(double x) {
    if (x < -1.0) { return M_PI; }
    if (x > 1.0) { return 0.0; }
    return std::acos(x);
}

std::pair<double, Gradient> deltaTheta(const NLoptVar &x,
                                       const InputKinematicsWithVertex &inp,
                                       const DecayChainOf &decay_chain) {
    auto parent = getParent(x, inp, decay_chain);
    auto parent_norm_vec = parent.normalize();

    double vdot;
    if (decay_chain == DecayChainOf::Parent1) {
        vdot = parent_norm_vec.dot(inp.vertex1());
    } else {
        vdot = parent_norm_vec.dot(inp.vertex2());
    }

    double d = 1.0 - vdot * vdot;
    if (d < 1.0e-10) { d = 1.0e-10; }
    double dacos = -1.0 / std::sqrt(d);

    double dk1x, dk1y, dk1z, dk2z;
    if (decay_chain == DecayChainOf::Parent1) {
        dk1x = inp.vertex1().px() - parent_norm_vec.px() * vdot;
        dk1y = inp.vertex1().py() - parent_norm_vec.py() * vdot;
        dk1z = inp.vertex1().pz() - parent_norm_vec.pz() * vdot;
        dk2z = 0.0;
    } else {
        dk1x = -inp.vertex2().px() + parent_norm_vec.px() * vdot;
        dk1y = -inp.vertex2().py() + parent_norm_vec.py() * vdot;
        dk1z = 0.0;
        dk2z = inp.vertex2().pz() - parent_norm_vec.pz() * vdot;
    }

    Gradient grad{dk1x, dk1y, dk1z, dk2z};
    grad *= dacos / parent.norm();

    return {acos(vdot), grad};
}

double constraintVertex(const NLoptVar &x, NLoptVar &grad, void *input,
                        const DecayChainOf &decay_chain) {
    const auto *inp = reinterpret_cast<InputKinematicsWithVertex *>(input);
    const auto &[delta_theta, d] = deltaTheta(x, *inp, decay_chain);

    if (!grad.empty()) { d.set_gradient(grad); }

    return delta_theta;
}

double constraintVertex1(const NLoptVar &x, NLoptVar &grad, void *input) {
    return constraintVertex(x, grad, input, DecayChainOf::Parent1);
}

double constraintVertex2(const NLoptVar &x, NLoptVar &grad, void *input) {
    return constraintVertex(x, grad, input, DecayChainOf::Parent2);
}

double constraintVertex1Upper(const NLoptVar &x, NLoptVar &grad, void *input) {
    const auto *inp = reinterpret_cast<InputKinematicsWithVertex *>(input);
    const auto &[delta_theta, d] = deltaTheta(x, *inp, DecayChainOf::Parent1);

    if (!grad.empty()) { d.set_gradient(grad); }

    return delta_theta - inp->delta_theta();
}

double constraintVertex1Lower(const NLoptVar &x, NLoptVar &grad, void *input) {
    const auto *inp = reinterpret_cast<InputKinematicsWithVertex *>(input);
    auto [delta_theta, d] = deltaTheta(x, *inp, DecayChainOf::Parent1);

    if (!grad.empty()) {
        d = -d;
        d.set_gradient(grad);
    }

    return -delta_theta - inp->delta_theta();
}
}  // namespace yam2
