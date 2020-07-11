/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#ifndef YAM2_SRC_YAM2_H_
#define YAM2_SRC_YAM2_H_

#include <ostream>

#include "input.h"

#include "invisibles.h"  // Invisibles, mkInvisibles
#include "momentum.h"    // FourMomentum
#include "variables.h"   // Variables

namespace yam2 {
class M2Solution {
private:
    double m2_;
    Invisibles ksol_;
    int neval_objf_;

public:
    M2Solution() = delete;
    M2Solution(const InputKinematics &inp, const Variables &sol, double m2,
               int neval_objf)
        : m2_(m2),
          ksol_(mkInvisibles(inp, sol, inp.scale())),
          neval_objf_(neval_objf) {}

    double m2() const { return m2_; }
    FourMomentum k1() const { return ksol_.k1(); }
    FourMomentum k2() const { return ksol_.k2(); }
    int neval_objf() const { return neval_objf_; }

    friend std::ostream &operator<<(std::ostream &os, const M2Solution &sol);
};

/** the default tolerance value */
constexpr double EPS = 1.0e-3;

/** the maximum number of interations */
constexpr int NEVAL = 5000;

/** M2XX variable using the SQP method */
std::optional<M2Solution> m2XXSQP(const std::optional<InputKinematics> &inp,
                                  double eps = EPS, int neval = NEVAL);
/** M2CX variable using the SQP method */
std::optional<M2Solution> m2CXSQP(const std::optional<InputKinematics> &inp,
                                  double eps = EPS, int neval = NEVAL);
/** M2XC variable using the SQP method */
std::optional<M2Solution> m2XCSQP(const std::optional<InputKinematics> &inp,
                                  double eps = EPS, int neval = NEVAL);
/** M2CC variable using the SQP method */
std::optional<M2Solution> m2CCSQP(const std::optional<InputKinematics> &inp,
                                  double eps = EPS, int neval = NEVAL);

/** M2XX variable using the augmented Lagrangian method with the BFGS update */
std::optional<M2Solution> m2XXAugLagBFGS(
    const std::optional<InputKinematics> &inp, double eps = EPS,
    int neval = NEVAL);
/** M2CX variable using the augmented Lagrangian method with the BFGS update */
std::optional<M2Solution> m2CXAugLagBFGS(
    const std::optional<InputKinematics> &inp, double eps = EPS,
    int neval = NEVAL);
/** M2XC variable using the augmented Lagrangian method with the BFGS update */
std::optional<M2Solution> m2XCAugLagBFGS(
    const std::optional<InputKinematics> &inp, double eps = EPS,
    int neval = NEVAL);
/** M2CC variable using the augmented Lagrangian method with the BFGS update */
std::optional<M2Solution> m2CCAugLagBFGS(
    const std::optional<InputKinematics> &inp, double eps = EPS,
    int neval = NEVAL);

/**
 * M2XX variable using the augmented Lagrangian method
 * with the Nelder-Mead simplex
 */
std::optional<M2Solution> m2XXAugLagNMSimplex(
    const std::optional<InputKinematics> &inp, double eps = EPS,
    int neval = NEVAL);
/**
 * M2CX variable using the augmented Lagrangian method
 * with the Nelder-Mead simplex
 */
std::optional<M2Solution> m2CXAugLagNMSimplex(
    const std::optional<InputKinematics> &inp, double eps = EPS,
    int neval = NEVAL);
/**
 * M2XC variable using the augmented Lagrangian method
 * with the Nelder-Mead simplex
 */
std::optional<M2Solution> m2XCAugLagNMSimplex(
    const std::optional<InputKinematics> &inp, double eps = EPS,
    int neval = NEVAL);
/**
 * M2CC variable using the augmented Lagrangian method
 * with the Nelder-Mead simplex
 */
std::optional<M2Solution> m2CCAugLagNMSimplex(
    const std::optional<InputKinematics> &inp, double eps = EPS,
    int neval = NEVAL);
}  // namespace yam2

#endif  // namespace yam2
