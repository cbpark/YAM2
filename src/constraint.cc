#include "constraint.h"

#include "gradient.h"
#include "input.h"
#include "invisibles.h"
#include "momentum.h"
#include "variables.h"

namespace yam2 {
FGType constraintF(const FourMomentum &p1, const FourMomentum &p2,
                   const InputKinematics &inp, const Variables &var) {
    const auto ks = mkInvisibles(inp, var);
    const auto m2grad = m2Grad(inp, ks, p1, p2, var);
    const auto &[grads, m1, m2] = m2grad;
    const auto &[grad1, grad2] = grads;
    return {m1 - m2, grad1 - grad2};
}

FGType constrantA(const InputKinematics &inp, const Variables &var) {
    return constraintF(inp.p1(), inp.p2(), inp, var);
}

FGType constrantB(const InputKinematics &inp, const Variables &var) {
    return constraintF(inp.q1(), inp.q2(), inp, var);
}
}  // namespace yam2
