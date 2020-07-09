#include "object.h"

#include "gradient.h"    // m2Grad
#include "input.h"       // InputKinematics
#include "invisibles.h"  // mkInvisibles
#include "variables.h"   // mkVariables

#include <utility>  // std::as_const

using std::vector;

namespace yam2 {
NLoptFunc m2ObjF(const InputKinematics &inp) {
    auto objf = [&inp = std::as_const(inp)](const vector<double> &x,
                                            vector<double> &grad, void *) {
        const auto var = mkVariables(x);
        const auto ks = mkInvisibles(inp, var.value());

        const auto &[grads, m1, m2] =
            m2Grad(inp, ks, inp.p1(), inp.p2(), var.value());
        const auto &[grad1, grad2] = grads;

        if (!grad.empty()) {
            if (m1 < m2) {
                grad = grad2.gradient();
            } else {
                grad = grad1.gradient();
            }
        }
        return m1 < m2 ? m2 : m1;
    };

    return objf;
}
}  // namespace yam2
