#ifndef YAM2_SRC_CONSTRAINT_H_
#define YAM2_SRC_CONSTRAINT_H_

#include "input.h"   // InputKinematics
#include "object.h"  // NLFunc

#include <vector>

namespace yam2 {
using Constraints = std::vector<NLoptFunc>;

/** constraint: (p1 + k1)^2 = (p2 + k2)^2 */
NLoptFunc constraintA(const InputKinematics &inp);

/** constraint: (q1 + k1)^2 = (q2 + k2)^2 */
NLoptFunc constraintB(const InputKinematics &inp);
}  // namespace yam2

#endif  // YAM2_SRC_CONSTRAINT_H_
