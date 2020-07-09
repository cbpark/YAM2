/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "object.h"

#include "gradient.h"    // m2Grad
#include "input.h"       // InputKinematics
#include "invisibles.h"  // mkInvisibles
#include "variables.h"   // mkVariables

namespace yam2 {
double m2ObjF(const std::vector<double> &x, std::vector<double> &grad,
              void *input) {
    const auto var = mkVariables(x);
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    const auto ks = mkInvisibles(*inp, var.value());

    const auto &[grads, m1, m2] =
        m2Grad(*inp, ks, inp->p1(), inp->p2(), var.value());
    const auto &[grad1, grad2] = grads;

    if (!grad.empty()) {
        if (m1 < m2) {
            grad = grad2.gradient();
        } else {
            grad = grad1.gradient();
        }
    }
    return m1 < m2 ? m2 : m1;
}
}  // namespace yam2
