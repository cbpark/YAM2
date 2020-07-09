#include "object.h"

#include "gradient.h"
#include "input.h"
#include "invisibles.h"
#include "variables.h"

namespace yam2 {
FGType m2ObjF(const InputKinematics &inp, const Variables &var) {
    const auto ks = mkInvisibles(inp, var);
    const auto &[grads, m1, m2] = m2Grad(inp, ks, inp.p1(), inp.p2(), var);
    const auto &[grad1, grad2] = grads;

    if (m1 < m2) { return {m2, grad2}; }
    return {m1, grad1};
}
}  // namespace yam2
