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
    unsigned int neval_objf_;

public:
    M2Solution() = delete;
    M2Solution(const InputKinematics &inp, const Variables &sol, double m2,
               unsigned int neval_objf)
        : m2_(m2),
          ksol_(mkInvisibles(inp, sol, inp.scale())),
          neval_objf_(neval_objf) {}

    double m2() const { return m2_; }
    FourMomentum k1() const { return ksol_.k1(); }
    FourMomentum k2() const { return ksol_.k2(); }
    unsigned int neval_objf() const { return neval_objf_; }
    void set_neval_objf(int n) { neval_objf_ = n; }

    friend bool operator<(const M2Solution &sol1, const M2Solution &sol2) {
        return sol1.m2_ < sol2.m2_;
    }

    friend bool operator>(const M2Solution &sol1, const M2Solution &sol2) {
        return sol2 < sol1;
    }

    friend std::ostream &operator<<(std::ostream &os, const M2Solution &sol);
};

/** the default tolerance value */
constexpr double EPS = 1.0e-3;

/** the default tolerance value for M2XX
 *
 *  it needs the smaller tolerance than others as the M2XX variables is
 *  actually obtained by unconstrained minimization.
 */
constexpr double EPSX = 1.0e-6;

/** the default tolerance value for M2Cons */
constexpr double EPSCONS = 1.0e-9;

/** the maximum number of iterations */
constexpr int NEVAL = 5000;

/** The M2XX variable.
 *
 *  It returns the smaller value between the results from SQP and
 *  augmented Lagrangian methods.
 */
std::optional<M2Solution> m2XX(const std::optional<InputKinematics> &inp,
                               double eps = EPSX, unsigned int neval = NEVAL);

/** The M2CX variable */
std::optional<M2Solution> m2CX(const std::optional<InputKinematics> &inp,
                               double eps = EPS, unsigned int neval = NEVAL);

/** The M2XC variable */
std::optional<M2Solution> m2XC(const std::optional<InputKinematics> &inp,
                               double eps = EPS, unsigned int neval = NEVAL);

/** The M2CC variable */
std::optional<M2Solution> m2CC(const std::optional<InputKinematics> &inp,
                               double eps = EPS, unsigned int neval = NEVAL);

/** The M2CR variable */
std::optional<M2Solution> m2CR(const std::optional<InputKinematics> &inp,
                               double eps = EPS, unsigned int neval = NEVAL);

/** The M2Cons variable */
std::optional<M2Solution> m2Cons(const std::optional<InputKinematics> &inp,
                                 double eps = EPSCONS,
                                 unsigned int neval = NEVAL);

/** The M2Cons variable */
std::optional<M2Solution> m2CCons(const std::optional<InputKinematics> &inp,
                                  double eps = EPSCONS,
                                  unsigned int neval = NEVAL);

/** The M2 variable with the vertex constraints */
std::optional<M2Solution> m2VertexEq(
    const std::optional<InputKinematicsWithVertex> &inp, double eps = EPSCONS,
    unsigned int neval = NEVAL);

/** The M2 variable with the sqrt(s) and vertex constraints */
std::optional<M2Solution> m2ConsVertexEq(
    const std::optional<InputKinematicsWithVertex> &inp, double eps = EPSCONS,
    unsigned int neval = NEVAL);

/** M2XX variable using the SQP method */
std::optional<M2Solution> m2XXSQP(const std::optional<InputKinematics> &inp,
                                  double eps = EPSX,
                                  unsigned int neval = NEVAL);

/** M2CX variable using the SQP method */
std::optional<M2Solution> m2CXSQP(const std::optional<InputKinematics> &inp,
                                  double eps = EPS, unsigned int neval = NEVAL);

/** M2XC variable using the SQP method */
std::optional<M2Solution> m2XCSQP(const std::optional<InputKinematics> &inp,
                                  double eps = EPS, unsigned int neval = NEVAL);

/** M2CC variable using the SQP method */
std::optional<M2Solution> m2CCSQP(const std::optional<InputKinematics> &inp,
                                  double eps = EPS, unsigned int neval = NEVAL);

/** M2CR variable using the SQP method */
std::optional<M2Solution> m2CRSQP(const std::optional<InputKinematics> &inp,
                                  double eps = EPS, unsigned int neval = NEVAL);

/** M2Cons variable using the SQP method */
std::optional<M2Solution> m2ConsSQP(const std::optional<InputKinematics> &inp,
                                    double eps = EPSCONS,
                                    unsigned int neval = NEVAL);

/** M2CCons variable using the SQP method */
std::optional<M2Solution> m2CConsSQP(const std::optional<InputKinematics> &inp,
                                     double eps = EPSCONS,
                                     unsigned int neval = NEVAL);

/** M2VertexEq variable using the SQP method */
std::optional<M2Solution> m2VertexEqSQP(
    const std::optional<InputKinematicsWithVertex> &inp, double eps = EPSCONS,
    unsigned int neval = NEVAL);

/** M2ConsVertexEq variable using the SQP method */
std::optional<M2Solution> m2ConsVertexEqSQP(
    const std::optional<InputKinematicsWithVertex> &inp, double eps = EPSCONS,
    unsigned int neval = NEVAL);

/** M2XX variable using the augmented Lagrangian method with the BFGS update */
std::optional<M2Solution> m2XXAugLagBFGS(
    const std::optional<InputKinematics> &inp, double eps = EPSX,
    unsigned int neval = NEVAL);

/** M2CX variable using the augmented Lagrangian method with the BFGS update */
std::optional<M2Solution> m2CXAugLagBFGS(
    const std::optional<InputKinematics> &inp, double eps = EPS,
    unsigned int neval = NEVAL);

/** M2XC variable using the augmented Lagrangian method with the BFGS update */
std::optional<M2Solution> m2XCAugLagBFGS(
    const std::optional<InputKinematics> &inp, double eps = EPS,
    unsigned int neval = NEVAL);

/** M2CC variable using the augmented Lagrangian method with the BFGS update */
std::optional<M2Solution> m2CCAugLagBFGS(
    const std::optional<InputKinematics> &inp, double eps = EPS,
    unsigned int neval = NEVAL);

/** M2CR variable using the augmented Lagrangian method with the BFGS update */
std::optional<M2Solution> m2CRAugLagBFGS(
    const std::optional<InputKinematics> &inp, double eps = EPS,
    unsigned int neval = NEVAL);

/**
 * M2Cons variable using the augmented Lagrangian method with the BFGS update
 */
std::optional<M2Solution> m2ConsAugLagBFGS(
    const std::optional<InputKinematics> &inp, double eps = EPSCONS,
    unsigned int neval = NEVAL);

/**
 * M2CCons variable using the augmented Lagrangian method with the BFGS update
 */
std::optional<M2Solution> m2CConsAugLagBFGS(
    const std::optional<InputKinematics> &inp, double eps = EPSCONS,
    unsigned int neval = NEVAL);

/**
 * M2VertexEq variable using the augmented Lagrangian method with the BFGS
 * update
 */
std::optional<M2Solution> m2VertexEqAugLagBFGS(
    const std::optional<InputKinematics> &inp, double eps = EPSCONS,
    unsigned int neval = NEVAL);

/**
 * M2ConsVertexEq variable using the augmented Lagrangian method with the BFGS
 * update
 */
std::optional<M2Solution> m2ConsVertexEqAugLagBFGS(
    const std::optional<InputKinematics> &inp, double eps = EPSCONS,
    unsigned int neval = NEVAL);

/**
 * M2XX variable using the augmented Lagrangian method
 * with the Nelder-Mead simplex
 */
std::optional<M2Solution> m2XXAugLagNMSimplex(
    const std::optional<InputKinematics> &inp, double eps = EPSX,
    unsigned int neval = NEVAL);

/**
 * M2CX variable using the augmented Lagrangian method
 * with the Nelder-Mead simplex
 */
std::optional<M2Solution> m2CXAugLagNMSimplex(
    const std::optional<InputKinematics> &inp, double eps = EPS,
    unsigned int neval = NEVAL);

/**
 * M2XC variable using the augmented Lagrangian method
 * with the Nelder-Mead simplex
 */
std::optional<M2Solution> m2XCAugLagNMSimplex(
    const std::optional<InputKinematics> &inp, double eps = EPS,
    unsigned int neval = NEVAL);

/**
 * M2CC variable using the augmented Lagrangian method
 * with the Nelder-Mead simplex
 */
std::optional<M2Solution> m2CCAugLagNMSimplex(
    const std::optional<InputKinematics> &inp, double eps = EPS,
    unsigned int neval = NEVAL);

/**
 * M2CR variable using the augmented Lagrangian method
 * with the Nelder-Mead simplex
 */
std::optional<M2Solution> m2CRAugLagNMSimplex(
    const std::optional<InputKinematics> &inp, double eps = EPS,
    unsigned int neval = NEVAL);

/**
 * M2Cons variable using the augmented Lagrangian method
 * with the Nelder-Mead simplex
 */
std::optional<M2Solution> m2ConsAugLagNMSimplex(
    const std::optional<InputKinematics> &inp, double eps = EPSCONS,
    unsigned int neval = NEVAL);

/**
 * M2CCons variable using the augmented Lagrangian method
 * with the Nelder-Mead simplex
 */
std::optional<M2Solution> m2CConsAugLagNMSimplex(
    const std::optional<InputKinematics> &inp, double eps = EPSCONS,
    unsigned int neval = NEVAL);

/**
 * M2VertexEq variable using the augmented Lagrangian method
 * with the Nelder-Mead simplex
 */
std::optional<M2Solution> m2VertexEqAugLagNMSimplex(
    const std::optional<InputKinematics> &inp, double eps = EPSCONS,
    unsigned int neval = NEVAL);

/**
 * M2ConsVertexEq variable using the augmented Lagrangian method
 * with the Nelder-Mead simplex
 */
std::optional<M2Solution> m2ConsVertexEqAugLagNMSimplex(
    const std::optional<InputKinematics> &inp, double eps = EPSCONS,
    unsigned int neval = NEVAL);
}  // namespace yam2

#endif  // namespace yam2
