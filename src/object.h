#ifndef YAM2_SRC_OBJECT_H_
#define YAM2_SRC_OBJECT_H_

#include "gradient.h"
#include "input.h"
#include "variables.h"

namespace yam2 {
FGType m2ObjF(const InputKinematics& inp, const Variables& var);
}  // namespace yam2

#endif  // YAM2_SRC_OBJECT_H_
