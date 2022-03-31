/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "YAM2/constraint.h"

#include <cmath>              // std::pow, std::tan, std::acos
#include <utility>            // std::pair
#include "YAM2/input.h"       // InputKinematics, deltaSqrtS
#include "YAM2/invisibles.h"  // mkInvisibles
#include "YAM2/momentum.h"    // SpatialMomentum
#include "YAM2/variables.h"   // mkVariables, NLoptVar
#include "gradient.h"         // m2Func, m2Func1, m2Func2

namespace yam2 {
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
                    GradFunc fmGrad, const NLoptVar &x, NLoptVar &grad) {
    const auto var = mkVariables(x);
    const auto var_val = var.value();
    const auto ks = mkInvisibles(inp, var_val);

    const auto &[grad_, m_] = fmGrad(inp, p, ks);
    if (!grad.empty()) { grad_.set_gradient(grad); }
    return m_ - inp.mparent().value_or(Mass{m_}).value;
}

double constraintA1(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraintAP(*inp, inp->p1(), m2Func1, x, grad);
}

double constraintA2(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraintAP(*inp, inp->p2(), m2Func2, x, grad);
}

double constraintR(const InputKinematics &inp, const FourMomentum &q,
                   GradFunc fmGrad, const NLoptVar &x, NLoptVar &grad) {
    const auto var = mkVariables(x);
    const auto var_val = var.value();
    const auto ks = mkInvisibles(inp, var_val);
    const auto &[grad_, m_] = fmGrad(inp, q, ks);

    if (!grad.empty()) { grad_.set_gradient(grad); }
    return m_ - inp.mrel().value_or(Mass{m_}).value;
}

double constraintR1(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraintR(*inp, inp->q1(), m2Func1, x, grad);
}

double constraintR2(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraintR(*inp, inp->q2(), m2Func2, x, grad);
}

double constraintSqrtS(const NLoptVar &x, NLoptVar &grad, void *input) {
    return deltaSqrtS(x, grad, input);
}

SpatialMomentum getParent1(const NLoptVar &x, const InputKinematics &inp) {
    const auto var = mkVariables(x);
    const auto var_val = var.value();
    const auto ks = mkInvisibles(inp, var_val);
    return (inp.p1() + ks.k1()).three_momentum();
}

double constraintVertex1Theta(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematicsWithVertex *>(input);
    auto parent1 = getParent1(x, *inp);
    double theta = parent1.theta();

    if (!grad.empty()) {
        double pt = parent1.pt();
        double r = std::pow(std::tan(theta), 2);
        double denom = (1.0 + r) * parent1.pz();

        double dk1x = parent1.px() / (pt * denom);
        double dk1y = parent1.py() / (pt * denom);
        double dk1z = -pt / (denom * parent1.pz());
        double dk2z = 0.0;

        auto grad_ = Gradient(dk1x, dk1y, dk1z, dk2z);
        grad_.set_gradient(grad);
    }

    return theta - inp->vertex1().theta();
}

double constraintVertex1Phi(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematicsWithVertex *>(input);
    auto parent1 = getParent1(x, *inp);

    if (!grad.empty()) {
        double px2 = std::pow(parent1.px(), 2);
        double py2px2 = 1.0 + std::pow(parent1.py(), 2) / px2;

        double dk1x = -parent1.py() / px2 / py2px2;
        double dk1y = 1.0 / parent1.px() / py2px2;
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

std::pair<double, Gradient> dotVertex1(const NLoptVar &x,
                                       const InputKinematicsWithVertex &inp) {
    auto parent1 = getParent1(x, inp);
    auto parent1_norm_vec = parent1.normalize();

    double vdot = parent1_norm_vec.dot(inp.vertex1());

    double d = 1.0 - vdot * vdot;
    if (d < 1.0e-10) { d = 1.0e-10; }
    double dacos = -1.0 / std::sqrt(d);

    double dk1x = inp.vertex1().px() - parent1_norm_vec.px() * vdot;
    double dk1y = inp.vertex1().py() - parent1_norm_vec.py() * vdot;
    double dk1z = inp.vertex1().pz() - parent1_norm_vec.pz() * vdot;
    double dk2z = 0.0;
    Gradient grad{dk1x, dk1y, dk1z, dk2z};
    grad *= dacos / parent1.norm();

    return {vdot, grad};
}

double constraintVertex1Upper(const NLoptVar &x, NLoptVar &grad, void *input) {
    const auto *inp = reinterpret_cast<InputKinematicsWithVertex *>(input);
    const auto &[vdot, d] = dotVertex1(x, *inp);

    if (!grad.empty()) { d.set_gradient(grad); }

    return acos(vdot) - inp->delta_theta_max();
}

double constraintVertex1Lower(const NLoptVar &x, NLoptVar &grad, void *input) {
    const auto *inp = reinterpret_cast<InputKinematicsWithVertex *>(input);
    auto [vdot, d] = dotVertex1(x, *inp);

    if (!grad.empty()) {
        d = -d;
        d.set_gradient(grad);
    }

    return -acos(vdot) - inp->delta_theta_max();
}
}  // namespace yam2
