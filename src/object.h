#ifndef YAM2_SRC_OBJECT_H_
#define YAM2_SRC_OBJECT_H_

#include "input.h"

#include <functional>
#include <vector>

namespace yam2 {
using NLoptFunc = std::function<double(const std::vector<double> &,
                                       std::vector<double> &, void *)>;

NLoptFunc m2ObjF(const InputKinematics &inp);
}  // namespace yam2

#endif  // YAM2_SRC_OBJECT_H_
