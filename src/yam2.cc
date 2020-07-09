/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "yam2.h"

#include <nlopt.hpp>  // for the NLopt API
#include <optional>
#include <ostream>

#include "constraint.h"  // Constraint
#include "input.h"       // InputKinematics
#include "object.h"      // m2ObjF
#include "variables.h"   // Variables, initialGuess

using std::optional;
using std::vector;

namespace yam2 {
/*
 * 'Constriant' type is the function pointer defined in 'constrain.h'.
 * 'InputKinematics' type is defined in 'input.h'
 */
optional<M2Solution> m2SQP(const vector<Constraint> &cfs,
                           const optional<InputKinematics> &inp, double eps) {
    if (!inp) { return {}; }

    auto inpv = inp.value();
    nlopt::opt algorithm{nlopt::LD_SLSQP, 4};  // the subproblem is BFGS.
    algorithm.set_min_objective(m2ObjF, &inpv);

    const double epsf = eps * 1.0e-2;
    algorithm.set_ftol_rel(epsf);
    algorithm.set_ftol_abs(epsf);

    for (const auto &cf : cfs) {
        algorithm.add_equality_constraint(cf, &inpv, eps);
    }

    auto x = initialGuess(inpv);  // x = variables = (k1x, k1y, k1z, k2z).
    double minf;  // minf = the minimum value of the objective function.
    auto result = algorithm.optimize(x, minf);
    if (result < 0) { return {}; }  // if failed, return nothing.

    minf *= inpv.scale();  // back to original scale.
    const auto sol_vars = mkVariables(x);
    // the solutions will be all back to the original scale.
    const M2Solution sol{inpv, sol_vars.value(), minf};
    return sol;
}

optional<M2Solution> m2XXSQP(const optional<InputKinematics> &inp, double eps) {
    return m2SQP(vector<Constraint>(), inp, eps);
}

optional<M2Solution> m2CXSQP(const optional<InputKinematics> &inp, double eps) {
    const vector<Constraint> constraint{constraintA};
    return m2SQP(constraint, inp, eps);
}

optional<M2Solution> m2XCSQP(const optional<InputKinematics> &inp, double eps) {
    const vector<Constraint> constraint{constraintB};
    return m2SQP(constraint, inp, eps);
}

optional<M2Solution> m2CCSQP(const optional<InputKinematics> &inp, double eps) {
    const vector<Constraint> constraint{constraintA, constraintB};
    return m2SQP(constraint, inp, eps);
}

optional<M2Solution> m2AugLag(const vector<Constraint> &cfs,
                              const optional<InputKinematics> &inp,
                              double eps) {
    if (!inp) { return {}; }

    auto inpv = inp.value();
    nlopt::opt subproblem{nlopt::LN_NELDERMEAD, 4};
    subproblem.set_min_objective(m2ObjF, &inpv);

    const double epsf = eps * 1.0e-2;
    subproblem.set_ftol_rel(epsf);
    subproblem.set_ftol_abs(epsf);

    nlopt::opt algorithm{nlopt::LD_AUGLAG_EQ, 4};
    algorithm.set_min_objective(m2ObjF, &inpv);
    algorithm.set_local_optimizer(subproblem);

    algorithm.set_ftol_rel(epsf);
    algorithm.set_ftol_abs(epsf);

    for (const auto &cf : cfs) {
        algorithm.add_equality_constraint(cf, &inpv, eps);
    }

    auto x = initialGuess(inpv);
    double minf;
    auto result = algorithm.optimize(x, minf);
    if (result < 0) { return {}; }

    minf *= inpv.scale();
    const auto sol_vars = mkVariables(x);
    const M2Solution sol{inpv, sol_vars.value(), minf};
    return sol;
}

optional<M2Solution> m2XXAugLag(const optional<InputKinematics> &inp,
                                double eps) {
    return m2AugLag(vector<Constraint>(), inp, eps);
}

optional<M2Solution> m2CXAugLag(const optional<InputKinematics> &inp,
                                double eps) {
    const vector<Constraint> constraint{constraintA};
    return m2AugLag(constraint, inp, eps);
}

optional<M2Solution> m2XCAugLag(const optional<InputKinematics> &inp,
                                double eps) {
    const vector<Constraint> constraint{constraintB};
    return m2AugLag(constraint, inp, eps);
}

optional<M2Solution> m2CCAugLag(const optional<InputKinematics> &inp,
                                double eps) {
    const vector<Constraint> constraint{constraintA, constraintB};
    return m2AugLag(constraint, inp, eps);
}

std::ostream &operator<<(std::ostream &os, const M2Solution &sol) {
    os << "M2: " << sol.m2() << '\n'
       << "k1: " << sol.k1() << '\n'
       << "k2: " << sol.k2();
    return os;
}
}  // namespace yam2
