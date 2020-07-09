#ifndef YAM2_SRC_CONSTRAINT_H_
#define YAM2_SRC_CONSTRAINT_H_

#include "momentum.h"
#include "input.h"
#include "variables.h"
#include "gradient.h"

#include <utility>

namespace yam2 {
/** constraint: (p1 + k1)^2 = (p2 + k2)^2 */
FGType constrantA(const InputKinematics& inp, const Variables& var);

/** constraint: (q1 + k1)^2 = (q2 + k2)^2 */
FGType constrantB(const InputKinematics& inp, const Variables& var);
}  // namespace yam2

#endif  // YAM2_SRC_CONSTRAINT_H_
