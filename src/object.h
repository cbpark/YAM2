/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#ifndef YAM2_SRC_OBJECT_H_
#define YAM2_SRC_OBJECT_H_

#include <vector>

namespace yam2 {
double m2ObjF(const std::vector<double> &x, std::vector<double> &grad,
              void *input);
}  // namespace yam2

#endif  // YAM2_SRC_OBJECT_H_
