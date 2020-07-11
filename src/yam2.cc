/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "yam2.h"

#include <nlopt.hpp>  // for the NLopt API

#include <exception>  // std::exception
#include <memory>
#include <optional>
#include <ostream>
#include <tuple>  // std::tuple

#include "constraint.h"  // Constraints
#include "gradient.h"    // m2Grad
#include "input.h"       // InputKinematics
#include "invisibles.h"  // mkInvisibles
#include "variables.h"   // Variables, mkVariables, NLoptVar

using std::optional;

/** to count the number of object function evaluations. */
int neval_objf = 0;

namespace yam2 {
double mTot(const NLoptVar &x, NLoptVar &grad, void *input) {
    const auto var = mkVariables(x);
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    const auto ks = mkInvisibles(*inp, var.value());

    const auto p1 = inp->p1(), p2 = inp->p2();
    const auto tot = p1 + p2 + ks.k1() + ks.k2();
    const double m_tot = tot.m();
    if (!grad.empty()) {
        const auto grad_mtot = mtotGrad(*inp, p1, p2, ks, var.value(), m_tot);
        grad = grad_mtot.gradient();
    }

    return m_tot;
}

/**
 *  the initial configuration that minimizes M_{tot}^2 = (p1 + p2 + k1 + k2)^2.
 */
NLoptVar initialGuessMtot(InputKinematics &inp, double eps, int neval) {
    // nlopt::opt algorithm{nlopt::LD_SLSQP, 4};
    // the simplex method seems to be fine enough.
    nlopt::opt algorithm{nlopt::LN_NELDERMEAD, 4};
    algorithm.set_min_objective(mTot, &inp);
    const double epsf = eps * 0.1;  // we don't have to be very accurate here.
    algorithm.set_ftol_rel(epsf);
    algorithm.set_ftol_abs(epsf);
    algorithm.set_maxeval(neval);

    const NLoptVar x0{0.5 * inp.ptmiss().px(), 0.5 * inp.ptmiss().py(), 0.0,
                      0.0};
    auto x{x0};
    double minf;
    auto result = algorithm.optimize(x, minf);
    if (result < 0) { return x0; }
    return x;
}

/**
 * The objective function for the M2 variable
 */
double m2ObjF(const NLoptVar &x, NLoptVar &grad, void *input) {
    ++neval_objf;

    const auto var = mkVariables(x);
    auto *const inp = reinterpret_cast<InputKinematics *>(input);
    const auto ks = mkInvisibles(*inp, var.value());

    const auto &[grads, m1, m2] =
        m2Grad(*inp, inp->p1(), inp->p2(), ks, var.value());
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
/*
 * 'Constriant' type is the function pointer defined in 'constrain.h'.
 * 'InputKinematics' type is defined in 'input.h'
 */
optional<M2Solution> m2SQP(const Constraints &cfs,
                           const optional<InputKinematics> &inp, double eps,
                           int neval) {
    if (!inp) { return {}; }

    auto inpv = inp.value();
    nlopt::opt algorithm{nlopt::LD_SLSQP, 4};  // the subproblem is BFGS.
    algorithm.set_min_objective(m2ObjF, &inpv);
    algorithm.set_maxeval(neval);

    const double epsf = eps * 1.0e-3;
    algorithm.set_ftol_rel(epsf);
    algorithm.set_ftol_abs(epsf);

    for (const auto &cf : cfs) {
        algorithm.add_equality_constraint(cf, &inpv, eps);
    }

    // x: variables = (k1x, k1y, k1z, k2z).
    auto x0 = initialGuessMtot(inpv, eps, neval);
    optional<nlopt::opt> subproblem;
    // minf = the minimum value of the objective function.
    const auto &[result, minf, x] =
        doOptimize(inpv, algorithm, subproblem, x0, epsf);

    if (result < 0) { return {}; }  // if failed, return empty.

    const auto sol_vars = mkVariables(x);
    // the solutions will be all back to the original scale.
    const M2Solution sol{inpv, sol_vars.value(), minf * inpv.scale(),
                         neval_objf};

    return sol;
}

optional<M2Solution> m2XXSQP(const optional<InputKinematics> &inp, double eps,
                             int neval) {
    return m2SQP(Constraints(), inp, eps, neval);
}

optional<M2Solution> m2CXSQP(const optional<InputKinematics> &inp, double eps,
                             int neval) {
    const Constraints constraint{constraintA};
    return m2SQP(constraint, inp, eps, neval);
}

optional<M2Solution> m2XCSQP(const optional<InputKinematics> &inp, double eps,
                             int neval) {
    const Constraints constraint{constraintB};
    return m2SQP(constraint, inp, eps, neval);
}

optional<M2Solution> m2CCSQP(const optional<InputKinematics> &inp, double eps,
                             int neval) {
    const Constraints constraint{constraintA, constraintB};
    return m2SQP(constraint, inp, eps, neval);
}

optional<M2Solution> m2AugLag(const nlopt::algorithm &subopt,
                              const Constraints &cfs,
                              const optional<InputKinematics> &inp, double eps,
                              int neval) {
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

    auto x0 = initialGuessMtot(inpv, eps, neval);
    const auto &[result, minf, x] =
        doOptimize(inpv, algorithm, subproblem, x0, epsf);

    if (result < 0) { return {}; }

    const auto sol_vars = mkVariables(x);
    const M2Solution sol{inpv, sol_vars.value(), minf * inpv.scale(),
                         neval_objf};
    return sol;
}

optional<M2Solution> m2AugLagBFGS(const Constraints &cfs,
                                  const optional<InputKinematics> &inp,
                                  double eps, int neval) {
    // Here, the SLSQP method is the BFGS in essence.
    return m2AugLag(nlopt::LD_SLSQP, cfs, inp, eps, neval);
}

optional<M2Solution> m2XXAugLagBFGS(const optional<InputKinematics> &inp,
                                    double eps, int neval) {
    return m2AugLagBFGS(Constraints(), inp, eps, neval);
}

optional<M2Solution> m2CXAugLagBFGS(const optional<InputKinematics> &inp,
                                    double eps, int neval) {
    const Constraints constraint{constraintA};
    return m2AugLagBFGS(constraint, inp, eps, neval);
}

optional<M2Solution> m2XCAugLagBFGS(const optional<InputKinematics> &inp,
                                    double eps, int neval) {
    const Constraints constraint{constraintB};
    return m2AugLagBFGS(constraint, inp, eps, neval);
}

optional<M2Solution> m2CCAugLagBFGS(const optional<InputKinematics> &inp,
                                    double eps, int neval) {
    const Constraints constraint{constraintA, constraintB};
    return m2AugLagBFGS(constraint, inp, eps, neval);
}

optional<M2Solution> m2AugLagNMSimplex(const Constraints &cfs,
                                       const optional<InputKinematics> &inp,
                                       double eps, int neval) {
    return m2AugLag(nlopt::LN_NELDERMEAD, cfs, inp, eps, neval);
}

optional<M2Solution> m2XXAugLagNMSimplex(const optional<InputKinematics> &inp,
                                         double eps, int neval) {
    return m2AugLagNMSimplex(Constraints(), inp, eps, neval);
}

optional<M2Solution> m2CXAugLagNMSimplex(const optional<InputKinematics> &inp,
                                         double eps, int neval) {
    const Constraints constraint{constraintA};
    return m2AugLagNMSimplex(constraint, inp, eps, neval);
}

optional<M2Solution> m2XCAugLagNMSimplex(const optional<InputKinematics> &inp,
                                         double eps, int neval) {
    const Constraints constraint{constraintB};
    return m2AugLagNMSimplex(constraint, inp, eps, neval);
}

optional<M2Solution> m2CCAugLagNMSimplex(const optional<InputKinematics> &inp,
                                         double eps, int neval) {
    const Constraints constraint{constraintA, constraintB};
    return m2AugLagNMSimplex(constraint, inp, eps, neval);
}

std::ostream &operator<<(std::ostream &os, const M2Solution &sol) {
    os << "-- found minimum after " << sol.neval_objf() << " evaluations:\n";
    os << "M2: " << sol.m2() << '\n'
       << "k1: " << sol.k1() << '\n'
       << "k2: " << sol.k2();
    return os;
}
}  // namespace yam2
