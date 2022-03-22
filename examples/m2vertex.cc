/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include <iostream>
#include "yam2.h"

int main() {
    /**
     *  For B1 + B2 --> v1 X1 + v2 X2,
     *
     *  where X1 and X2 are invisible particles.
     */
    // FourMomentum{energy, px, py, pz}
    const yam2::FourMomentum v1{4.533286, -0.234916, 0.282849, 0.230556};
    const yam2::FourMomentum v2{
        3.472357,
        -1.195970,
        -0.0656993,
        -1.419437,
    };

    // TransverseMomentum{px, py}
    const yam2::TransverseMomentum ptmiss{1.430886, -0.217149};
    const yam2::Mass m_invis{0.0};  // M_{X1} = M_{X2} = 0.

    const yam2::Mass m_parent{5.279};

    const double sqrt_s = 10.583;

    const yam2::SpatialMomentum vertex1{0.317964, -0.0472086, -0.169791};
    const yam2::SpatialMomentum vertex2{-0.317964, 0.0472086, 0.169791};
    const double delta_theta_max = 0.15;

    // we have only one-step decay, so it's necessary to fake the second step.
    const auto zero = yam2::FourMomentum();

    const auto input =
        yam2::mkInput({v1, v2}, {zero, zero}, ptmiss, m_invis, vertex1, vertex2,
                      delta_theta_max, {m_parent}, {}, sqrt_s, {});
    if (!input) {
        std::cerr << "Invalid input.\n";
        return 1;
    }
    std::cout << "-- process information (scaled):\n" << input.value() << '\n';

    // // the default tolerance is 1.0e-9.
    // // const auto m2sol = yam2::m2ConsSQP(input, 1.0e-9);
    // const auto m2sol = yam2::m2ConsSQP(input);

    // // the other available methods are:
    // // - the augmented Lagrangian method with BFGS:
    // // const auto m2sol = yam2::m2ConsAugLagBFGS(input);
    // // - the augmented Lagrangian method with Nelder-Mead simplex:
    // // const auto m2sol = yam2::m2ConsAugLagNMSimplex(input);
    // // - the combination of all the methods (caution: slow but more
    // accurate):
    // // const auto m2sol = yam2::m2Cons(input);
    // if (!m2sol) {
    //     std::cerr << "Failed to find minimum.\n";
    //     return 1;
    // } else {
    //     // std::cout << m2sol.value() << '\n';

    //     /* k1 and k2 are the M2 solutions to the momenta of the invisible
    //      * particles B1 and B2.
    //      */
    //     std::cout << "M2Cons = " << m2sol.value().m2() << '\n'
    //               << "where \n"
    //               << "  k1: " << m2sol.value().k1() << '\n'
    //               << "  k2: " << m2sol.value().k2() << '\n'
    //               << "found after " << m2sol.value().neval_objf()
    //               << " evaluations.\n";
    // }
}
