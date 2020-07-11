/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "constraint.h"

#include "gradient.h"    // m2Grad
#include "input.h"       // InputKinematics
#include "invisibles.h"  // mkInvisibles
#include "variables.h"   // mkVariables, NLoptVar

namespace yam2 {
double constraint(const FourMomentum &p1, const FourMomentum &p2,
                  const NLoptVar &x, NLoptVar &grad, void *input) {
    const auto var = mkVariables(x);
    const auto var_val = var.value();
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    const auto ks = mkInvisibles(*inp, var_val);

    const auto &[grads, m1, m2] = m2Grad(*inp, p1, p2, ks, var_val);
    const auto &[grad1, grad2] = grads;

    if (!grad.empty()) {
        const auto dgrad = grad1 - grad2;
        grad = dgrad.gradient();
    }
    return m1 - m2;
}

double constraintA(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    const auto p1 = inp->p1();
    const auto p2 = inp->p2();
    return constraint(p1, p2, x, grad, input);
}

double constraintB(const NLoptVar &x, NLoptVar &grad, void *input) {
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    const auto q1 = inp->q1();
    const auto q2 = inp->q2();
    return constraint(q1, q2, x, grad, input);
}
}  // namespace yam2
