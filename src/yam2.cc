/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "yam2.h"

#include <nlopt.hpp>  // for the NLopt API

#include <exception>  // std::exception
#include <functional>
#include <optional>
#include <ostream>
#include <tuple>  // std::tuple

#include "constraint.h"  // Constraints
#include "gradient.h"    // m2Func
#include "input.h"       // InputKinematics
#include "invisibles.h"  // mkInvisibles
#include "variables.h"   // Variables, mkVariables, NLoptVar

using std::optional;

/** to count the number of object function evaluations. */
unsigned int neval_objf = 0;

namespace yam2 {
using OptInp = optional<InputKinematics>;
using OptM2 = optional<M2Solution>;

/**
 * The objective function for the M2 variable
 */
double m2ObjF(const NLoptVar &x, NLoptVar &grad, void *input) {
    ++neval_objf;

    const auto var = mkVariables(x);
    const auto var_val = var.value();
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    const auto ks = mkInvisibles(*inp, var_val);

    const auto &[grads, m1, m2] =
        m2Func(*inp, inp->p1(), inp->p2(), ks, var_val);
    const auto &[grad1, grad2] = grads;

    if (!grad.empty()) {
        if (m1 < m2) {
            grad = grad2.gradient();
        } else {
            grad = grad1.gradient();
        }
    }
    return m1 < m2 ? m2 : m1;
}

/**
 *  If NLopt throws an exception (such as nlopt::roundoff_limited),
 *  rescale the tolerance (* 10) and re-run the optimizer.
 */
std::tuple<nlopt::result, double, NLoptVar> doOptimize(
    InputKinematics &inp, nlopt::opt &algorithm,
    optional<nlopt::opt> &subproblem, const NLoptVar &x0, double epsf) {
    if (!subproblem) {
        algorithm.set_ftol_rel(epsf);
        algorithm.set_ftol_abs(epsf);
    } else {
        subproblem.value().set_ftol_rel(epsf);
        subproblem.value().set_ftol_abs(epsf);
        // perhaps, the outer algorithm doesn't need to rescale the tolerance.
        // (check!)
        algorithm.set_local_optimizer(subproblem.value());
    }

    nlopt::result result;
    double minf;
    auto x{x0};
    try {
        neval_objf = 0;  // reset the number of evaluations
        result = algorithm.optimize(x, minf);
    } catch (std::exception &) {
        epsf *= 10.0;
        doOptimize(inp, algorithm, subproblem, x0, epsf);
    }
    return {result, minf, x};
}
/**
 *  'Constriant' type is the function pointer defined in 'constrain.h'.
 *  'InputKinematics' type is defined in 'input.h'
 */
OptM2 m2SQP(const Constraints &cfs, const OptInp &inp, double eps,
            unsigned int neval) {
    if (!inp) { return {}; }

    auto inpv = inp.value();
    nlopt::opt algorithm{nlopt::LD_SLSQP, 4};  // the subproblem is BFGS.
    algorithm.set_min_objective(m2ObjF, &inpv);
    algorithm.set_maxeval(neval);

    for (const auto &cf : cfs) {
        algorithm.add_equality_constraint(cf, &inpv, eps);
    }

    // x: variables = (k1x, k1y, k1z, k2z).
    auto x0 = inpv.initial_guess(eps, neval);
    optional<nlopt::opt> subproblem;
    const double epsf = eps * 1.0e-3;
    // minf = the minimum value of the objective function.
    const auto &[result, minf, x] =
        doOptimize(inpv, algorithm, subproblem, x0, epsf);

    // if M2 <= 0, it has been failed to find the minimum.
    if (result < 0 || minf <= 0) { return {}; }  // if failed, return empty.

    const auto sol_vars = mkVariables(x);
    // the solutions will be all back to the original scale.
    return M2Solution{inpv, sol_vars.value(), minf * inpv.scale(), neval_objf};
}

OptM2 m2XXSQP(const OptInp &inp, double eps, unsigned int neval) {
    return m2SQP(Constraints(), inp, eps, neval);
}

OptM2 m2CXSQP(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint{constraintA};
    return m2SQP(constraint, inp, eps, neval);
}

OptM2 m2XCSQP(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint{constraintB};
    return m2SQP(constraint, inp, eps, neval);
}

OptM2 m2CCSQP(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint{constraintA, constraintB};
    return m2SQP(constraint, inp, eps, neval);
}

OptM2 m2CRSQP(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint{constraintA, constraintR1, constraintR2};
    return m2SQP(constraint, inp, eps, neval);
}

OptM2 m2AugLag(const nlopt::algorithm &subopt, const Constraints &cfs,
               const OptInp &inp, double eps, unsigned int neval) {
    if (!inp) { return {}; }

    auto inpv = inp.value();
    auto subproblem = std::make_optional<nlopt::opt>(subopt, 4);
    subproblem.value().set_min_objective(m2ObjF, &inpv);
    subproblem.value().set_maxeval(neval * 10);

    nlopt::opt algorithm{nlopt::AUGLAG_EQ, 4};
    algorithm.set_min_objective(m2ObjF, &inpv);
    algorithm.set_local_optimizer(subproblem.value());
    algorithm.set_maxeval(neval);

    const double epsf = eps * 1.0e-3;
    algorithm.set_ftol_rel(epsf);
    algorithm.set_ftol_abs(epsf);

    for (const auto &cf : cfs) {
        algorithm.add_equality_constraint(cf, &inpv, eps);
    }

    auto x0 = inpv.initial_guess(eps, neval);
    const auto &[result, minf, x] =
        doOptimize(inpv, algorithm, subproblem, x0, epsf);

    if (result < 0 || minf <= 0) { return {}; }

    const auto sol_vars = mkVariables(x);
    return M2Solution{inpv, sol_vars.value(), minf * inpv.scale(), neval_objf};
}

OptM2 m2AugLagBFGS(const Constraints &cfs, const OptInp &inp, double eps,
                   unsigned int neval) {
    // Here, the SLSQP method is the BFGS in essence.
    return m2AugLag(nlopt::LD_SLSQP, cfs, inp, eps, neval);
}

OptM2 m2XXAugLagBFGS(const OptInp &inp, double eps, unsigned int neval) {
    return m2AugLagBFGS(Constraints(), inp, eps, neval);
}

OptM2 m2CXAugLagBFGS(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint{constraintA};
    return m2AugLagBFGS(constraint, inp, eps, neval);
}

OptM2 m2XCAugLagBFGS(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint{constraintB};
    return m2AugLagBFGS(constraint, inp, eps, neval);
}

OptM2 m2CCAugLagBFGS(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint{constraintA, constraintB};
    return m2AugLagBFGS(constraint, inp, eps, neval);
}

OptM2 m2CRAugLagBFGS(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint{constraintA, constraintR1, constraintR2};
    return m2AugLagBFGS(constraint, inp, eps, neval);
}

OptM2 m2AugLagNMSimplex(const Constraints &cfs, const OptInp &inp, double eps,
                        unsigned int neval) {
    return m2AugLag(nlopt::LN_NELDERMEAD, cfs, inp, eps, neval);
}

OptM2 m2XXAugLagNMSimplex(const OptInp &inp, double eps, unsigned int neval) {
    return m2AugLagNMSimplex(Constraints(), inp, eps, neval);
}

OptM2 m2CXAugLagNMSimplex(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint{constraintA};
    return m2AugLagNMSimplex(constraint, inp, eps, neval);
}

OptM2 m2XCAugLagNMSimplex(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint{constraintB};
    return m2AugLagNMSimplex(constraint, inp, eps, neval);
}

OptM2 m2CCAugLagNMSimplex(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint{constraintA, constraintB};
    return m2AugLagNMSimplex(constraint, inp, eps, neval);
}

OptM2 m2CRAugLagNMSimplex(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint{constraintA, constraintR1, constraintR2};
    return m2AugLagNMSimplex(constraint, inp, eps, neval);
}

using M2Func = std::function<OptM2(const OptInp &, double, int)>;

OptM2 m2(M2Func fSQP, M2Func fAugLagBFGS, M2Func fAugLagNMSimplex,
         const OptInp &inp, double eps, unsigned int neval) {
    auto m2_sqp = fSQP(inp, eps, neval);
    auto m2_auglag_bfgs = fAugLagBFGS(inp, eps, neval);

    // add the number of obj fn evals for both methods.
    unsigned int neval_objf_tot = 0;
    if (m2_sqp) { neval_objf_tot += m2_sqp.value().neval_objf(); }
    if (m2_auglag_bfgs) {
        neval_objf_tot += m2_auglag_bfgs.value().neval_objf();
    }

    // if both SQP and A-BFGS were successful
    if (m2_sqp && m2_auglag_bfgs) {
        m2_sqp.value().set_neval_objf(neval_objf_tot);
        m2_auglag_bfgs.value().set_neval_objf(neval_objf_tot);
        // return the one with smaller M2.
        return m2_sqp.value() < m2_auglag_bfgs.value() ? m2_sqp
                                                       : m2_auglag_bfgs;
    } else if (m2_sqp && !m2_auglag_bfgs) {
        return m2_sqp;
    } else if (!m2_sqp && m2_auglag_bfgs) {
        return m2_auglag_bfgs;
    }

    // if both SQP and A-BFGS failed, use A-Simplex.
    return fAugLagNMSimplex(inp, eps, neval);
}

OptM2 m2XX(const OptInp &inp, double eps, unsigned int neval) {
    return m2(m2XXSQP, m2XXAugLagBFGS, m2XXAugLagNMSimplex, inp, eps, neval);
}

OptM2 m2CX(const OptInp &inp, double eps, unsigned int neval) {
    return m2(m2CXSQP, m2CXAugLagBFGS, m2CXAugLagNMSimplex, inp, eps, neval);
}

OptM2 m2XC(const OptInp &inp, double eps, unsigned int neval) {
    return m2(m2XCSQP, m2XCAugLagBFGS, m2XCAugLagNMSimplex, inp, eps, neval);
}

OptM2 m2CC(const OptInp &inp, double eps, unsigned int neval) {
    return m2(m2CCSQP, m2CCAugLagBFGS, m2CCAugLagNMSimplex, inp, eps, neval);
}

OptM2 m2CR(const OptInp &inp, double eps, unsigned int neval) {
    return m2(m2CRSQP, m2CRAugLagBFGS, m2CRAugLagNMSimplex, inp, eps, neval);
}

std::ostream &operator<<(std::ostream &os, const M2Solution &sol) {
    os << "-- found minimum after " << sol.neval_objf() << " evaluations:\n";
    os << "M2: " << sol.m2() << '\n'
       << "k1: " << sol.k1() << '\n'
       << "k2: " << sol.k2();
    return os;
}
}  // namespace yam2
