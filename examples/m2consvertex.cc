/*
 *  Copyright (c) 2022 Chan Beom Park <cbpark@gmail.com>
 */

#include <iostream>
#include "YAM2/yam2.h"

int main() {
    /**
     *  For B1 + B2 --> v1 X1 + v2 X2,
     *
     *  where X1 and X2 are invisible particles.
     */
    // FourMomentum{energy, px, py, pz}
    const yam2::FourMomentum v1{4.53329, -0.23492, 0.28285, 0.23056};
    const yam2::FourMomentum v2{3.47236, -1.19597, -0.065699, -1.41944};

    // TransverseMomentum{px, py}
    const yam2::TransverseMomentum ptmiss{1.43089, -0.21715};

    const yam2::Mass m_invis{0.0};  // M_{X1} = M_{X2} = 0.

    const yam2::Mass m_parent{5.279};  // M_{B1} = M_{B2} = 5.279

    // sqrt(s) of B1 + B2.
    const double sqrt_s = 10.583;

    // the longitudinal momentum of B1 + B2.
    const double ptot_z = 0.0;  // CM frame

    const yam2::SpatialMomentum vertex1{0.31796, -0.047209, -0.16979};
    const yam2::SpatialMomentum vertex2{-0.31796, 0.047209, 0.16979};

    const auto input =
        yam2::mkInputWithVertex(v1, v2, ptmiss, m_invis, vertex1, vertex2, 0.0,
                                {m_parent}, sqrt_s, {ptot_z});
    if (!input) {
        std::cerr << "Invalid input.\n";
        return 1;
    }
    std::cout << "-- process information (scaled):\n" << input.value() << '\n';

    // const auto m2sol = yam2::m2ConsVertexEqSQP(input, 1.0e-6, 1000);
    const auto m2sol = yam2::m2ConsVertexEq(input, 1.0e-6, 1000);

    if (!m2sol) {
        std::cerr << "Failed to find minimum.\n";
        return 1;
    } else {
        // std::cout << m2sol.value() << '\n';

        std::cout << "M2ConsVertexEq = " << m2sol.value().m2() << '\n'
                  << "where \n"
                  << "  k1: " << m2sol.value().k1() << '\n'
                  << "  k2: " << m2sol.value().k2() << '\n'
                  << "found after " << m2sol.value().neval_objf()
                  << " evaluations.\n";
    }
}
