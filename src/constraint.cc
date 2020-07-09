/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "constraint.h"

#include <vector>

#include "gradient.h"    // m2Grad
#include "input.h"       // InputKinematics
#include "invisibles.h"  // mkInvisibles
#include "variables.h"   // mkVariables

using std::vector;

namespace yam2 {
double constraintA(const vector<double> &x, vector<double> &grad, void *input) {
    const auto var = mkVariables(x);
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    const auto ks = mkInvisibles(*inp, var.value());

    const auto &[grads, m1, m2] =
        m2Grad(*inp, ks, inp->p1(), inp->p2(), var.value());
    const auto &[grad1, grad2] = grads;

    if (!grad.empty()) {
        const auto dgrad = grad1 - grad2;
        grad = dgrad.gradient();
    }
    return m1 - m2;
}

double constraintB(const vector<double> &x, vector<double> &grad, void *input) {
    const auto var = mkVariables(x);
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    const auto ks = mkInvisibles(*inp, var.value());

    const auto &[grads, m1, m2] =
        m2Grad(*inp, ks, inp->q1(), inp->q2(), var.value());
    const auto &[grad1, grad2] = grads;

    if (!grad.empty()) {
        const auto dgrad = grad1 - grad2;
        grad = dgrad.gradient();
    }
    return m1 - m2;
}
}  // namespace yam2
