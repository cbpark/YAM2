/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#ifndef YAM2_SRC_CONSTRAINT_H_
#define YAM2_SRC_CONSTRAINT_H_

#include <vector>
#include "variables.h"  // NLoptVar

namespace yam2 {
typedef double (*Constraint)(const NLoptVar &x, NLoptVar &grad, void *input);

using Constraints = std::vector<Constraint>;

/** constraint: (p1 + k1)^2 = (p2 + k2)^2 */
double constraintA(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: (q1 + k1)^2 = (q2 + k2)^2 */
double constraintB(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: (q1 + k1)^2 = M_{R}^2 */
double constraintR1(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: (q2 + k2)^2 = M_{R}^2 */
double constraintR2(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: (p1 + p2 + k1 + k2)^2 = s */
double constraintSqrtS(const NLoptVar &x, NLoptVar &grad, void *input);
}  // namespace yam2

#endif  // YAM2_SRC_CONSTRAINT_H_
