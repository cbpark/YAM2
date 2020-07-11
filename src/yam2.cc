/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "yam2.h"

#include <nlopt.hpp>  // for the NLopt API
#include <optional>
#include <ostream>

#include "constraint.h"  // Constraint
#include "gradient.h"    // m2Grad
#include "input.h"       // InputKinematics
#include "invisibles.h"  // mkInvisibles
#include "variables.h"   // Variables, mkVariables

using std::optional;
using std::vector;

int neval_objf = 0;  // to count the number of object function evaluations.

namespace yam2 {
double mTot(const vector<double> &x, vector<double> &grad, void *input) {
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
vector<double> initialGuessMtot(InputKinematics &inp, double eps, int neval) {
    // nlopt::opt algorithm{nlopt::LD_SLSQP, 4};
    // the simplex method seems to be fine enough.
    nlopt::opt algorithm{nlopt::LN_NELDERMEAD, 4};
    algorithm.set_min_objective(mTot, &inp);
    const double epsf = eps * 0.1;  // we don't have to be very accurate here.
    algorithm.set_ftol_rel(epsf);
    algorithm.set_ftol_abs(epsf);
    algorithm.set_maxeval(neval);

    const vector<double> x0{0.5 * inp.ptmiss().px(), 0.5 * inp.ptmiss().py(),
                            0.0, 0.0};
    auto x = x0;
    double minf;
    auto result = algorithm.optimize(x, minf);
    if (result < 0) { return x0; }
    return x;
}

/**
 * The objective function for the M2 variable
 */
double m2ObjF(const vector<double> &x, vector<double> &grad, void *input) {
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

/*
 * 'Constriant' type is the function pointer defined in 'constrain.h'.
 * 'InputKinematics' type is defined in 'input.h'
 */
optional<M2Solution> m2SQP(const vector<Constraint> &cfs,
                           const optional<InputKinematics> &inp, double eps,
                           int neval) {
    if (!inp) { return {}; }

    neval_objf = 0;

    auto inpv = inp.value();
    nlopt::opt algorithm{nlopt::LD_SLSQP, 4};  // the subproblem is BFGS.
    algorithm.set_min_objective(m2ObjF, &inpv);

    const double epsf = eps * 1.0e-3;
    algorithm.set_ftol_rel(epsf);
    algorithm.set_ftol_abs(epsf);
    algorithm.set_maxeval(neval);

    for (const auto &cf : cfs) {
        algorithm.add_equality_constraint(cf, &inpv, eps);
    }

    // x = variables = (k1x, k1y, k1z, k2z).
    // auto x = initialGuess(inpv);
    auto x = initialGuessMtot(inpv, eps, neval);
    double minf;  // minf = the minimum value of the objective function.
    auto result = algorithm.optimize(x, minf);
    if (result < 0) { return {}; }  // if failed, return nothing.

    minf *= inpv.scale();  // back to original scale.
    const auto sol_vars = mkVariables(x);
    // the solutions will be all back to the original scale.
    const M2Solution sol{inpv, sol_vars.value(), minf, neval_objf};

    return sol;
}

optional<M2Solution> m2XXSQP(const optional<InputKinematics> &inp, double eps,
                             int neval) {
    return m2SQP(vector<Constraint>(), inp, eps, neval);
}

optional<M2Solution> m2CXSQP(const optional<InputKinematics> &inp, double eps,
                             int neval) {
    const vector<Constraint> constraint{constraintA};
    return m2SQP(constraint, inp, eps, neval);
}

optional<M2Solution> m2XCSQP(const optional<InputKinematics> &inp, double eps,
                             int neval) {
    const vector<Constraint> constraint{constraintB};
    return m2SQP(constraint, inp, eps, neval);
}

optional<M2Solution> m2CCSQP(const optional<InputKinematics> &inp, double eps,
                             int neval) {
    const vector<Constraint> constraint{constraintA, constraintB};
    return m2SQP(constraint, inp, eps, neval);
}

optional<M2Solution> m2AugLag(const nlopt::algorithm &subopt,
                              const vector<Constraint> &cfs,
                              const optional<InputKinematics> &inp, double eps,
                              int neval) {
    if (!inp) { return {}; }

    neval_objf = 0;

    auto inpv = inp.value();
    nlopt::opt subproblem{subopt, 4};
    subproblem.set_min_objective(m2ObjF, &inpv);

    const double epsf = eps * 1.0e-3;
    subproblem.set_ftol_rel(epsf);
    subproblem.set_ftol_abs(epsf);
    subproblem.set_maxeval(neval * 5);

    nlopt::opt algorithm{nlopt::AUGLAG_EQ, 4};
    algorithm.set_min_objective(m2ObjF, &inpv);
    algorithm.set_local_optimizer(subproblem);

    algorithm.set_ftol_rel(epsf);
    algorithm.set_ftol_abs(epsf);
    algorithm.set_maxeval(neval * 5);

    for (const auto &cf : cfs) {
        algorithm.add_equality_constraint(cf, &inpv, eps);
    }

    auto x = initialGuessMtot(inpv, eps, neval);
    double minf;
    auto result = algorithm.optimize(x, minf);
    if (result < 0) { return {}; }

    minf *= inpv.scale();
    const auto sol_vars = mkVariables(x);
    const M2Solution sol{inpv, sol_vars.value(), minf, neval_objf};

    return sol;
}

optional<M2Solution> m2AugLagBFGS(const vector<Constraint> &cfs,
                                  const optional<InputKinematics> &inp,
                                  double eps, int neval) {
    // Here, the SLSQP method is the BFGS in essence.
    return m2AugLag(nlopt::LD_SLSQP, cfs, inp, eps, neval);
}

optional<M2Solution> m2XXAugLagBFGS(const optional<InputKinematics> &inp,
                                    double eps, int neval) {
    return m2AugLagBFGS(vector<Constraint>(), inp, eps, neval);
}

optional<M2Solution> m2CXAugLagBFGS(const optional<InputKinematics> &inp,
                                    double eps, int neval) {
    const vector<Constraint> constraint{constraintA};
    return m2AugLagBFGS(constraint, inp, eps, neval);
}

optional<M2Solution> m2XCAugLagBFGS(const optional<InputKinematics> &inp,
                                    double eps, int neval) {
    const vector<Constraint> constraint{constraintB};
    return m2AugLagBFGS(constraint, inp, eps, neval);
}

optional<M2Solution> m2CCAugLagBFGS(const optional<InputKinematics> &inp,
                                    double eps, int neval) {
    const vector<Constraint> constraint{constraintA, constraintB};
    return m2AugLagBFGS(constraint, inp, eps, neval);
}

optional<M2Solution> m2AugLagNMSimplex(const vector<Constraint> &cfs,
                                       const optional<InputKinematics> &inp,
                                       double eps, int neval) {
    return m2AugLag(nlopt::LN_NELDERMEAD, cfs, inp, eps, neval);
}

optional<M2Solution> m2XXAugLagNMSimplex(const optional<InputKinematics> &inp,
                                         double eps, int neval) {
    return m2AugLagNMSimplex(vector<Constraint>(), inp, eps, neval);
}

optional<M2Solution> m2CXAugLagNMSimplex(const optional<InputKinematics> &inp,
                                         double eps, int neval) {
    const vector<Constraint> constraint{constraintA};
    return m2AugLagNMSimplex(constraint, inp, eps, neval);
}

optional<M2Solution> m2XCAugLagNMSimplex(const optional<InputKinematics> &inp,
                                         double eps, int neval) {
    const vector<Constraint> constraint{constraintB};
    return m2AugLagNMSimplex(constraint, inp, eps, neval);
}

optional<M2Solution> m2CCAugLagNMSimplex(const optional<InputKinematics> &inp,
                                         double eps, int neval) {
    const vector<Constraint> constraint{constraintA, constraintB};
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
