/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include "yam2.h"

#include <nlopt.hpp>  // for the NLopt API
#include <optional>
#include <ostream>

#include "constraint.h"  // Constraints
#include "input.h"       // InputKinematics
#include "object.h"      // m2ObjF
#include "variables.h"   // Variables, initialGuess

namespace yam2 {
std::optional<M2Solution> m2SQP(const std::vector<Constraint> &cfs,
                                const std::optional<InputKinematics> &inp,
                                double eps) {
    if (!inp) { return {}; }

    auto inpv = inp.value();

    auto opt = nlopt::opt{nlopt::LD_SLSQP, 4};
    opt.set_min_objective(m2ObjF, &inpv);

    const double epsf = eps * 1.0e-2;
    opt.set_ftol_rel(epsf);
    opt.set_ftol_abs(epsf);

    for (const auto &cf : cfs) { opt.add_equality_constraint(cf, &inpv, eps); }

    auto x = initialGuess(inpv);
    double minf;
    auto result = opt.optimize(x, minf);
    if (result < 0) { return {}; }

    minf *= inpv.scale();
    const auto sol_vars = mkVariables(x);
    const M2Solution sol{inpv, sol_vars.value(), minf};
    return sol;
}

std::optional<M2Solution> m2XXSQP(const std::optional<InputKinematics> &inp,
                                  double eps) {
    return m2SQP(std::vector<Constraint>(), inp, eps);
}

std::optional<M2Solution> m2CXSQP(const std::optional<InputKinematics> &inp,
                                  double eps) {
    const std::vector<Constraint> constraint{constraintA};
    return m2SQP(constraint, inp, eps);
}

std::optional<M2Solution> m2XCSQP(const std::optional<InputKinematics> &inp,
                                  double eps) {
    const std::vector<Constraint> constraint{constraintB};
    return m2SQP(constraint, inp, eps);
}

std::optional<M2Solution> m2CCSQP(const std::optional<InputKinematics> &inp,
                                  double eps) {
    const std::vector<Constraint> constraint{constraintA, constraintB};
    return m2SQP(constraint, inp, eps);
}

std::ostream &operator<<(std::ostream &os, const M2Solution &sol) {
    os << "M2: " << sol.m2() << '\n'
       << "k1: " << sol.k1() << '\n'
       << "k2: " << sol.k2();
    return os;
}
}  // namespace yam2
