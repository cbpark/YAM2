#include "constraint.h"

#include "gradient.h"    // m2Grad
#include "input.h"       // InputKinematics
#include "invisibles.h"  // mkInvisibles
#include "momentum.h"    // FourMomentum
#include "object.h"      // NLFunc
#include "variables.h"   // mkVariables

#include <utility>  // std::as_const
#include <vector>

using std::as_const;
using std::vector;

namespace yam2 {
NLoptFunc constraint(const FourMomentum &p1, const FourMomentum &p2,
                     const InputKinematics &inp) {
    auto consF = [&inp = as_const(inp), &p1 = as_const(p1), &p2 = as_const(p2)](
                     const vector<double> &x, vector<double> &grad, void *) {
        const auto var = mkVariables(x);
        const auto ks = mkInvisibles(inp, var.value());

        const auto &[grads, m1, m2] = m2Grad(inp, ks, p1, p2, var.value());
        const auto &[grad1, grad2] = grads;

        if (!grad.empty()) {
            const auto dgrad = grad1 - grad2;
            grad = dgrad.gradient();
        }
        return m1 - m2;
    };

    return consF;
}

NLoptFunc constraintA(const InputKinematics &inp) {
    return constraint(inp.p1(), inp.p2(), inp);
}

NLoptFunc constraintB(const InputKinematics &inp) {
    return constraint(inp.q1(), inp.q2(), inp);
}
}  // namespace yam2
