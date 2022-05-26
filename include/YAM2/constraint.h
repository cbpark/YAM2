/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#ifndef YAM2_SRC_CONSTRAINT_H_
#define YAM2_SRC_CONSTRAINT_H_

#include <vector>
#include "variables.h"  // NLoptVar

namespace yam2 {
using Constraint = double (*)(const NLoptVar &x, NLoptVar &grad, void *input);

using Constraints = std::vector<Constraint>;

/** constraint: (p1 + k1)^2 = (p2 + k2)^2 */
double constraintA(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: (q1 + k1)^2 = (q2 + k2)^2 */
double constraintB(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: (p1 + k1)^2 = M_{A1}^2 */
double constraintA1(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: (p1 + k1)^2 < M_{A1}^2 (1 + epsilon) */
double constraintA1Upper(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: (p1 + k1)^2 > M_{A1}^2 (1 - epsilon) */
double constraintA1Lower(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: (p2 + k2)^2 = M_{A2}^2 */
double constraintA2(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: (p2 + k2)^2 < M_{A2}^2 (1 + epsilon) */
double constraintA2Upper(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: (p2 + k2)^2 > M_{A2}^2 (1 - epsilon) */
double constraintA2Lower(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: (q1 + k1)^2 = M_{R}^2 */
double constraintR1(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: (q2 + k2)^2 = M_{R}^2 */
double constraintR2(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: (p1 + p2 + k1 + k2)^2 = s */
double constraintSqrtS(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: acos(parent1 . v) = 0 */
double constraintVertex1(const NLoptVar &x, NLoptVar &grad, void *input);

double constraintVertex1X(const NLoptVar &x, NLoptVar &grad, void *input);

double constraintVertex1Y(const NLoptVar &x, NLoptVar &grad, void *input);

double constraintVertex1Z(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: Delta theta(parent1, vertex1) = 0 */
double constraintVertex1Theta(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: Delta phi(parent1, vertex1) = 0 */
double constraintVertex1Phi(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: acos(parent2 . v) = 0 */
double constraintVertex2(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: acos(parent1 . v) < delta_theta */
double constraintVertex1Upper(const NLoptVar &x, NLoptVar &grad, void *input);

/** constraint: acos(parent1 . v) > -delta_theta */
double constraintVertex1Lower(const NLoptVar &x, NLoptVar &grad, void *input);
}  // namespace yam2

#endif  // YAM2_SRC_CONSTRAINT_H_
