#ifndef YAM2_SRC_CONSTRAINT_H_
#define YAM2_SRC_CONSTRAINT_H_

#include "gradient.h"
#include "input.h"
#include "momentum.h"
#include "variables.h"

#include <functional>
#include <utility>
#include <vector>

namespace yam2 {
using Constraints = std::vector<
    std::function<FGType(const InputKinematics &, const Variables &)>>;

/** constraint: (p1 + k1)^2 = (p2 + k2)^2 */
FGType constrantA(const InputKinematics &inp, const Variables &var);

/** constraint: (q1 + k1)^2 = (q2 + k2)^2 */
FGType constrantB(const InputKinematics &inp, const Variables &var);
}  // namespace yam2

#endif  // YAM2_SRC_CONSTRAINT_H_
