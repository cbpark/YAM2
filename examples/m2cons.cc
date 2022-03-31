/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include <iostream>
#include "YAM2/yam2.h"

int main() {
    /**
     *  For Y --> A1 + A2 --> a1 B1 + a2 B2,
     *
     *  where Y is an on-shell particle, B1 and B2 are invisible particles.
     */
    // FourMomentum{energy, px, py, pz}
    const yam2::FourMomentum a1{51.6869, 16.5476, 44.5988, -20.1837};
    const yam2::FourMomentum a2{9.4859, -2.4573, -7.7910, -4.8206};

    // TransverseMomentum{px, py}
    const yam2::TransverseMomentum ptmiss{-14.0903, -36.8077};
    const yam2::Mass m_invis{0.0};  // M_{B1} = M_{B2} = 0.

    const double mY = 125.0;

    // we have only one-step decay, so it's necessary to fake the second step.
    const auto zero = yam2::FourMomentum();

    // `{}` corresponds to the relative particle mass. In the decay topology
    // of interest, we don't have a relative particle.
    const auto input =
        yam2::mkInput({a1, a2}, {zero, zero}, ptmiss, m_invis, {}, {}, mY);
    // If you want set the longitudinal momentum of the total system (say 10
    // GeV), the input should be given by
    // const auto input =
    //     yam2::mkInput({a1, a2}, {zero, zero}, ptmiss, m_invis, {}, {}, mY,
    //     {10.0});
    if (!input) {
        std::cerr << "Invalid input.\n";
        return 1;
    }
    std::cout << "-- process information (scaled):\n" << input.value() << '\n';

    // the default tolerance is 1.0e-9.
    // const auto m2sol = yam2::m2ConsSQP(input, 1.0e-9);
    const auto m2sol = yam2::m2ConsSQP(input);

    // the other available methods are:
    // - the augmented Lagrangian method with BFGS:
    // const auto m2sol = yam2::m2ConsAugLagBFGS(input);
    // - the augmented Lagrangian method with Nelder-Mead simplex:
    // const auto m2sol = yam2::m2ConsAugLagNMSimplex(input);
    // - the combination of all the methods (caution: slow but more accurate):
    // const auto m2sol = yam2::m2Cons(input);
    if (!m2sol) {
        std::cerr << "Failed to find minimum.\n";
        return 1;
    } else {
        // std::cout << m2sol.value() << '\n';

        /* k1 and k2 are the M2 solutions to the momenta of the invisible
         * particles B1 and B2.
         */
        std::cout << "M2Cons = " << m2sol.value().m2() << '\n'
                  << "where \n"
                  << "  k1: " << m2sol.value().k1() << '\n'
                  << "  k2: " << m2sol.value().k2() << '\n'
                  << "found after " << m2sol.value().neval_objf()
                  << " evaluations.\n";
    }
}
