/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "constraint.h"

#include "gradient.h"    // m2Func, m2Func1, m2Func2
#include "input.h"       // InputKinematics, deltaSqrtS
#include "invisibles.h"  // mkInvisibles
#include "variables.h"   // mkVariables, NLoptVar

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
        grad = dgrad.gradient();
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

    const auto &[gradi, mi] = fmGrad(inp, p, ks);
    if (!grad.empty()) { grad = gradi.gradient(); }
    return mi - inp.mparent().value_or(Mass{mi}).value;
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
    const auto &[gradi, mi] = fmGrad(inp, q, ks);

    if (!grad.empty()) { grad = gradi.gradient(); }
    return mi - inp.mrel().value_or(Mass{mi}).value;
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
}  // namespace yam2
