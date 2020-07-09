/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#ifndef YAM2_SRC_CONSTRAINT_H_
#define YAM2_SRC_CONSTRAINT_H_

#include <vector>

namespace yam2 {
typedef double (*Constraint)(const std::vector<double> &x,
                             std::vector<double> &grad, void *input);

/** constraint: (p1 + k1)^2 = (p2 + k2)^2 */
double constraintA(const std::vector<double> &x, std::vector<double> &grad,
                   void *input);

/** constraint: (q1 + k1)^2 = (q2 + k2)^2 */
double constraintB(const std::vector<double> &x, std::vector<double> &grad,
                   void *input);
}  // namespace yam2

#endif  // YAM2_SRC_CONSTRAINT_H_
