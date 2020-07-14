/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "constraint.h"

#include "gradient.h"    // m2Grad, m2Grad1, m2Grad2
#include "input.h"       // InputKinematics
#include "invisibles.h"  // mkInvisibles
#include "variables.h"   // mkVariables, NLoptVar

namespace yam2 {
double constraint(const InputKinematics &inp, const FourMomentum &p1,
                  const FourMomentum &p2, const NLoptVar &x, NLoptVar &grad) {
    const auto var = mkVariables(x);
    const auto var_val = var.value();
    const auto ks = mkInvisibles(inp, var_val);

    const auto &[grads, m1, m2] = m2Grad(inp, p1, p2, ks, var_val);
    const auto &[grad1, grad2] = grads;

    if (!grad.empty()) {
        const auto dgrad = grad1 - grad2;
        grad = dgrad.gradient();
    }
    return m1 - m2;
}

double constraintA(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraint(*inp, inp->p1(), inp->p2(), x, grad);
}

double constraintB(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraint(*inp, inp->q1(), inp->q2(), x, grad);
}

double constraintR(const InputKinematics &inp, const FourMomentum &q,
                   GradFunc fmGrad, const NLoptVar &x, NLoptVar &grad) {
    const auto var = mkVariables(x);
    const auto var_val = var.value();
    const auto ks = mkInvisibles(inp, var_val);
    const auto &[gradi, mi] = fmGrad(inp, q, ks, var_val);

    if (!grad.empty()) { grad = gradi.gradient(); }
    return mi - inp.mrel().value;
}

double constraintR1(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraintR(*inp, inp->q1(), m2Grad1, x, grad);
}

double constraintR2(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    return constraintR(*inp, inp->q2(), m2Grad2, x, grad);
}
}  // namespace yam2
