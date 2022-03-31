/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "YAM2/yam2.h"

#include <nlopt.hpp>  // for the NLopt API

#include <algorithm>  // std::any_of
#include <cmath>      // std::is_nan, std::fabs
#include <exception>  // std::exception
#include <optional>
#include <ostream>
#include <tuple>  // std::tuple
#include <vector>

#include "YAM2/constraint.h"  // Constraints
#include "YAM2/input.h"       // InputKinematics
#include "YAM2/invisibles.h"  // mkInvisibles
#include "YAM2/variables.h"   // Variables, mkVariables, NLoptVar
#include "gradient.h"         // m2Func

using std::optional;

/** to count the number of object function evaluations. */
unsigned int neval_objf = 0;

namespace yam2 {
using OptInp = optional<InputKinematics>;
using OptInpWithVertex = optional<InputKinematicsWithVertex>;
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

    const auto &[grads, m1, m2] = m2Func(*inp, inp->p1(), inp->p2(), ks);
    const auto &[grad1, grad2] = grads;

    if (!grad.empty()) {
        if (m1 < m2) {
            grad2.set_gradient(grad);
        } else {
            grad1.set_gradient(grad);
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
    } catch (std::exception &e) {
        std::cerr << "doOptimize: exception from NLopt (" << e.what() << ")\n";
        epsf *= 10.0;
        doOptimize(inp, algorithm, subproblem, x0, epsf);
    }

    // the solution has NaN?
    bool nan_sol = std::any_of(x.cbegin(), x.cend(),
                               [](double xv) { return std::isnan(xv); });
    if (nan_sol) {
        std::cerr
            << "doOptimize: error (NaN solution)! we increase tolerance ...\n";
        doOptimize(inp, algorithm, subproblem, x0, epsf * 10.0);
    }

    // unphysical solution?
    bool invalid_sol = std::any_of(
        x.cbegin(), x.cend(), [](double xv) { return std::fabs(xv) > 1.0e10; });
    if (invalid_sol) {
        std::cerr << "doOptimize: error (invalid solution)! we increase "
                     "tolerance ...\n";
        doOptimize(inp, algorithm, subproblem, x0, epsf * 10.0);
    }

    return {result, minf, x};
}

/**
 *  'Constriant' type is the function pointer defined in 'constrain.h'.
 *  'InputKinematics' type is defined in 'input.h'
 */
template <typename Input>
OptM2 m2SQP(const Constraints &cfs_eq, const Constraints &cfs_ineq,
            const optional<Input> &inp, double eps, unsigned int neval) {
    if (!inp) { return {}; }

    auto inpv = inp.value();
    unsigned int dim;
    if (!inpv.ptot_z()) {
        dim = 4;
    } else {
        dim = 3;
    }

    nlopt::opt algorithm{nlopt::LD_SLSQP, dim};  // the subproblem is BFGS.
    algorithm.set_min_objective(m2ObjF, &inpv);
    algorithm.set_maxeval(neval);

    for (const auto &cf : cfs_eq) {
        algorithm.add_equality_constraint(cf, &inpv, eps);
    }
    for (const auto &cf : cfs_ineq) {
        algorithm.add_inequality_constraint(cf, &inpv, eps);
    }

    // x: variables = (k1x, k1y, k1z, k2z).
    auto x0 = inpv.initial_guess(eps, neval);
    x0.resize(dim);
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

template OptM2 m2SQP<InputKinematics>(const Constraints &cfs,
                                      const Constraints &cfs_ineq,
                                      const OptInp &inp, double eps,
                                      unsigned int neval);

template OptM2 m2SQP<InputKinematicsWithVertex>(const Constraints &cfs,
                                                const Constraints &cfs_ineq,
                                                const OptInpWithVertex &inp,
                                                double eps, unsigned int neval);

OptM2 m2XXSQP(const OptInp &inp, double eps, unsigned int neval) {
    return m2SQP(Constraints(), Constraints(), inp, eps, neval);
}

OptM2 m2CXSQP(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint_eq{constraintA};
    return m2SQP(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2XCSQP(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint_eq{constraintB};
    return m2SQP(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2CCSQP(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint_eq{constraintA, constraintB};
    return m2SQP(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2CRSQP(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint_eq{constraintA, constraintR1, constraintR2};
    return m2SQP(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2ConsSQP(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint_eq{constraintSqrtS};
    return m2SQP(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2CConsSQP(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint_eq{constraintSqrtS, constraintA1,
                                    constraintA2};
    return m2SQP(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2VertexEqSQP(const OptInpWithVertex &inp, double eps,
                    unsigned int neval) {
    const Constraints constraint_eq{constraintVertex1Theta,
                                    constraintVertex1Phi};
    return m2SQP(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2ConsVertexEqSQP(const OptInpWithVertex &inp, double eps,
                        unsigned int neval) {
    const Constraints constraint_eq{constraintSqrtS, constraintVertex1Theta,
                                    constraintVertex1Phi};
    return m2SQP(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2VertexIneqSQP(const OptInpWithVertex &inp, double eps,
                      unsigned int neval) {
    const Constraints constraint_ineq{constraintVertex1Upper,
                                      constraintVertex1Lower};
    return m2SQP(Constraints(), constraint_ineq, inp, eps, neval);
}

OptM2 m2ConsVertexIneqSQP(const OptInpWithVertex &inp, double eps,
                          unsigned int neval) {
    const Constraints constraint_eq{constraintSqrtS};
    const Constraints constraint_ineq{constraintVertex1Upper,
                                      constraintVertex1Lower};
    return m2SQP(constraint_eq, constraint_ineq, inp, eps, neval);
}

template <typename Input>
OptM2 m2AugLag(const nlopt::algorithm &subopt, const Constraints &cfs_eq,
               const Constraints &cfs_ineq, const optional<Input> &inp,
               double eps, unsigned int neval) {
    if (!inp) { return {}; }

    auto inpv = inp.value();
    unsigned int dim;
    if (!inpv.ptot_z()) {
        dim = 4;
    } else {
        dim = 3;
    }

    auto subproblem = std::make_optional<nlopt::opt>(subopt, dim);
    subproblem.value().set_min_objective(m2ObjF, &inpv);
    subproblem.value().set_maxeval(neval * 10);

    nlopt::opt algorithm{nlopt::AUGLAG, dim};
    algorithm.set_min_objective(m2ObjF, &inpv);
    algorithm.set_local_optimizer(subproblem.value());
    algorithm.set_maxeval(neval);

    const double epsf = eps * 1.0e-3;
    algorithm.set_ftol_rel(epsf);
    algorithm.set_ftol_abs(epsf);

    for (const auto &cf : cfs_eq) {
        algorithm.add_equality_constraint(cf, &inpv, eps);
    }
    for (const auto &cf : cfs_ineq) {
        algorithm.add_inequality_constraint(cf, &inpv, eps);
    }

    auto x0 = inpv.initial_guess(eps, neval);
    x0.resize(dim);
    const auto &[result, minf, x] =
        doOptimize(inpv, algorithm, subproblem, x0, epsf);

    if (result < 0 || minf <= 0) { return {}; }

    const auto sol_vars = mkVariables(x);
    return M2Solution{inpv, sol_vars.value(), minf * inpv.scale(), neval_objf};
}

template OptM2 m2AugLag<InputKinematics>(const nlopt::algorithm &subopt,
                                         const Constraints &cfs_eq,
                                         const Constraints &cfs_ineq,
                                         const OptInp &inp, double eps,
                                         unsigned int neval);

template OptM2 m2AugLag<InputKinematicsWithVertex>(
    const nlopt::algorithm &subopt, const Constraints &cfs_eq,
    const Constraints &cfs_ineq, const OptInpWithVertex &inp, double eps,
    unsigned int neval);

template <typename Input>
OptM2 m2AugLagBFGS(const Constraints &cfs_eq, const Constraints &cfs_ineq,
                   const optional<Input> &inp, double eps, unsigned int neval) {
    // Here, the SLSQP method is the BFGS in essence.
    return m2AugLag(nlopt::LD_SLSQP, cfs_eq, cfs_ineq, inp, eps, neval);
}

template OptM2 m2AugLagBFGS<InputKinematics>(const Constraints &cfs_eq,
                                             const Constraints &cfs_ineq,
                                             const OptInp &inp, double eps,
                                             unsigned int neval);

template OptM2 m2AugLagBFGS<InputKinematicsWithVertex>(
    const Constraints &cfs_eq, const Constraints &cfs_ineq,
    const OptInpWithVertex &inp, double eps, unsigned int neval);

OptM2 m2XXAugLagBFGS(const OptInp &inp, double eps, unsigned int neval) {
    return m2AugLagBFGS(Constraints(), Constraints(), inp, eps, neval);
}

OptM2 m2CXAugLagBFGS(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint_eq{constraintA};
    return m2AugLagBFGS(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2XCAugLagBFGS(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint_eq{constraintB};
    return m2AugLagBFGS(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2CCAugLagBFGS(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint_eq{constraintA, constraintB};
    return m2AugLagBFGS(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2CRAugLagBFGS(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint_eq{constraintA, constraintR1, constraintR2};
    return m2AugLagBFGS(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2ConsAugLagBFGS(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint_eq{constraintSqrtS};
    return m2AugLagBFGS(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2CConsAugLagBFGS(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint_eq{constraintSqrtS, constraintA1,
                                    constraintA2};
    return m2AugLagBFGS(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2VertexEqAugLagBFGS(const OptInpWithVertex &inp, double eps,
                           unsigned int neval) {
    const Constraints constraint_eq{constraintVertex1Theta,
                                    constraintVertex1Phi};
    return m2AugLagBFGS(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2ConsVertexEqAugLagBFGS(const OptInpWithVertex &inp, double eps,
                               unsigned int neval) {
    const Constraints constraint_eq{constraintSqrtS, constraintVertex1Theta,
                                    constraintVertex1Phi};
    return m2AugLagBFGS(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2VertexIneqAugLagBFGS(const OptInpWithVertex &inp, double eps,
                             unsigned int neval) {
    const Constraints constraint_ineq{constraintVertex1Upper,
                                      constraintVertex1Lower};
    return m2AugLagBFGS(Constraints(), constraint_ineq, inp, eps, neval);
}

OptM2 m2ConsVertexIneqAugLagBFGS(const OptInpWithVertex &inp, double eps,
                                 unsigned int neval) {
    const Constraints constraint_eq{constraintSqrtS};
    const Constraints constraint_ineq{constraintVertex1Upper,
                                      constraintVertex1Lower};
    return m2AugLagBFGS(constraint_eq, constraint_ineq, inp, eps, neval);
}

OptM2 m2AugLagNMSimplex(const Constraints &cfs_eq, const Constraints &cfs_ineq,
                        const OptInp &inp, double eps, unsigned int neval) {
    return m2AugLag(nlopt::LN_NELDERMEAD, cfs_eq, cfs_ineq, inp, eps, neval);
}

OptM2 m2XXAugLagNMSimplex(const OptInp &inp, double eps, unsigned int neval) {
    return m2AugLagNMSimplex(Constraints(), Constraints(), inp, eps, neval);
}

OptM2 m2CXAugLagNMSimplex(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint_eq{constraintA};
    return m2AugLagNMSimplex(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2XCAugLagNMSimplex(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint_eq{constraintB};
    return m2AugLagNMSimplex(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2CCAugLagNMSimplex(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint_eq{constraintA, constraintB};
    return m2AugLagNMSimplex(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2CRAugLagNMSimplex(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint_eq{constraintA, constraintR1, constraintR2};
    return m2AugLagNMSimplex(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2ConsAugLagNMSimplex(const OptInp &inp, double eps, unsigned int neval) {
    const Constraints constraint_eq{constraintSqrtS};
    return m2AugLagNMSimplex(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2CConsAugLagNMSimplex(const OptInp &inp, double eps,
                             unsigned int neval) {
    const Constraints constraint_eq{constraintSqrtS, constraintA1,
                                    constraintA2};
    return m2AugLagNMSimplex(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2VertexEqAugLagNMSimplex(const OptInpWithVertex &inp, double eps,
                                unsigned int neval) {
    const Constraints constraint_eq{constraintVertex1Theta,
                                    constraintVertex1Phi};
    return m2AugLagNMSimplex(constraint_eq, Constraints(), inp, eps, neval);
}

OptM2 m2ConsVertexEqAugLagNMSimplex(const OptInpWithVertex &inp, double eps,
                                    unsigned int neval) {
    const Constraints constraint_eq{constraintSqrtS, constraintVertex1Theta,
                                    constraintVertex1Phi};
    return m2AugLagNMSimplex(constraint_eq, Constraints(), inp, eps, neval);
}

template <typename Input>
using M2Func = OptM2 (*)(const optional<Input> &, double, unsigned int);

template <typename Input>
OptM2 m2MinStrategy1(M2Func<Input> fSQP, M2Func<Input> fAugLagBFGS,
                     M2Func<Input> fAugLagNMSimplex, const optional<Input> &inp,
                     double eps, unsigned int neval) {
    auto m2_sqp = fSQP(inp, eps, neval);
    auto m2_auglag_bfgs = fAugLagBFGS(inp, eps, neval);

    // add the number of obj fn evals for both methods.
    unsigned int neval_objf_tot = 0;
    if (m2_sqp) { neval_objf_tot += m2_sqp.value().neval_objf(); }
    if (m2_auglag_bfgs) {
        neval_objf_tot += m2_auglag_bfgs.value().neval_objf();
    }

    // if both SQP and AugLagBFGS were successful
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

    // if both SQP and AugLagBFGS failed, use AugLagNMSimplex.
    return fAugLagNMSimplex(inp, eps, neval);
}

template OptM2 m2MinStrategy1<InputKinematics>(
    M2Func<InputKinematics> fSQP, M2Func<InputKinematics> fAugLagBFGS,
    M2Func<InputKinematics> fAugLagNMSimplex, const OptInp &inp, double eps,
    unsigned int neval);

template OptM2 m2MinStrategy1<InputKinematicsWithVertex>(
    M2Func<InputKinematicsWithVertex> fSQP,
    M2Func<InputKinematicsWithVertex> fAugLagBFGS,
    M2Func<InputKinematicsWithVertex> fAugLagNMSimplex,
    const OptInpWithVertex &inp, double eps, unsigned int neval);

// Strategy 2: if the first algorithm works and it doesn't exceed the parent
// particle mass (input), return the result.
template <typename Input>
OptM2 m2MinStrategy2(const std::vector<M2Func<Input>> &f_algos,
                     const optional<Input> &inp, double eps,
                     unsigned int neval) {
    for (auto f_algo : f_algos) {
        auto m2sol = f_algo(inp, eps, neval);
        if (m2sol && inp) {
            auto inpv = inp.value();
            auto p_parent1 = inpv.p1() + m2sol.value().k1();
            auto p_parent2 = inpv.p2() + m2sol.value().k2();

            if (p_parent1.m() < inpv.mparent().value_or(Mass{1.0e+10}).value *
                                    inpv.scale() * 1.05 &&
                p_parent2.m() < inpv.mparent().value_or(Mass{1.0e+10}).value *
                                    inpv.scale() * 1.05) {
                return {m2sol};
            }
        }
    }

    return {};
}

template OptM2 m2MinStrategy2<InputKinematics>(
    const std::vector<M2Func<InputKinematics>> &f_algos, const OptInp &inp,
    double eps, unsigned int neval);

template OptM2 m2MinStrategy2<InputKinematicsWithVertex>(
    const std::vector<M2Func<InputKinematicsWithVertex>> &f_algos,
    const OptInpWithVertex &inp, double eps, unsigned int neval);

// Strategy 3: if the first algorithm in the list works, return the result.
template <typename Input>
OptM2 m2MinStrategy3(const std::vector<M2Func<Input>> &f_algos,
                     const optional<Input> &inp, double eps,
                     unsigned int neval) {
    for (auto f_algo : f_algos) {
        auto m2sol = f_algo(inp, eps, neval);
        if (m2sol) { return {m2sol}; }
    }

    return {};
}

template OptM2 m2MinStrategy3<InputKinematics>(
    const std::vector<M2Func<InputKinematics>> &f_algos, const OptInp &inp,
    double eps, unsigned int neval);

template OptM2 m2MinStrategy3<InputKinematicsWithVertex>(
    const std::vector<M2Func<InputKinematicsWithVertex>> &f_algos,
    const OptInpWithVertex &inp, double eps, unsigned int neval);

OptM2 m2XX(const OptInp &inp, double eps, unsigned int neval) {
    return m2MinStrategy1(m2XXSQP, m2XXAugLagBFGS, m2XXAugLagNMSimplex, inp,
                          eps, neval);
}

OptM2 m2CX(const OptInp &inp, double eps, unsigned int neval) {
    return m2MinStrategy1(m2CXSQP, m2CXAugLagBFGS, m2CXAugLagNMSimplex, inp,
                          eps, neval);
}

OptM2 m2XC(const OptInp &inp, double eps, unsigned int neval) {
    return m2MinStrategy1(m2XCSQP, m2XCAugLagBFGS, m2XCAugLagNMSimplex, inp,
                          eps, neval);
}

OptM2 m2CC(const OptInp &inp, double eps, unsigned int neval) {
    return m2MinStrategy1(m2CCSQP, m2CCAugLagBFGS, m2CCAugLagNMSimplex, inp,
                          eps, neval);
}

OptM2 m2CR(const OptInp &inp, double eps, unsigned int neval) {
    return m2MinStrategy1(m2CRSQP, m2CRAugLagBFGS, m2CRAugLagNMSimplex, inp,
                          eps, neval);
}

OptM2 m2Cons(const OptInp &inp, double eps, unsigned int neval) {
    return m2MinStrategy1(m2ConsSQP, m2ConsAugLagBFGS, m2ConsAugLagNMSimplex,
                          inp, eps, neval);
}

OptM2 m2CCons(const OptInp &inp, double eps, unsigned int neval) {
    std::vector<M2Func<InputKinematics>> f_algos{m2CConsSQP, m2CConsAugLagBFGS,
                                                 m2CConsAugLagNMSimplex};
    return m2MinStrategy2(f_algos, inp, eps, neval);
}

OptM2 m2VertexEq(const OptInpWithVertex &inp, double eps, unsigned int neval) {
    return m2MinStrategy1(m2VertexEqSQP, m2VertexEqAugLagBFGS,
                          m2VertexEqAugLagNMSimplex, inp, eps, neval);
}

OptM2 m2ConsVertexEq(const OptInpWithVertex &inp, double eps,
                     unsigned int neval) {
    std::vector<M2Func<InputKinematicsWithVertex>> f_algos{
        m2ConsVertexEqSQP, m2ConsVertexEqAugLagBFGS,
        m2ConsVertexEqAugLagNMSimplex};
    return m2MinStrategy2(f_algos, inp, eps, neval);
}

OptM2 m2VertexIneq(const OptInpWithVertex &inp, double eps,
                   unsigned int neval) {
    std::vector<M2Func<InputKinematicsWithVertex>> f_algos{
        m2VertexIneqAugLagBFGS, m2VertexIneqSQP};
    return m2MinStrategy3(f_algos, inp, eps, neval);
}

OptM2 m2ConsVertexIneq(const OptInpWithVertex &inp, double eps,
                       unsigned int neval) {
    std::vector<M2Func<InputKinematicsWithVertex>> f_algos{
        m2ConsVertexIneqAugLagBFGS, m2ConsVertexIneqSQP};
    return m2MinStrategy3(f_algos, inp, eps, neval);
}

std::ostream &operator<<(std::ostream &os, const M2Solution &sol) {
    os << "-- found minimum after " << sol.neval_objf() << " evaluations:\n";
    os << "M2: " << sol.m2() << '\n'
       << "k1: " << sol.k1() << '\n'
       << "k2: " << sol.k2();
    return os;
}
}  // namespace yam2
