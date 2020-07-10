/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#include <iostream>
#include "yam2.h"

int main() {
    /**
     *  For A1 + A2 --> a1 B1 + a2 B2 --> a1 b1 C1 + a2 b2 C2
     */
    const yam2::FourMomentum a1{80.0160, 50.1283, 33.2622, 52.5386};
    // FourMomentum{energy, px, py, pz}
    const yam2::FourMomentum a2{75.2869, -73.8496, 13.7107, 1.8279};
    const yam2::FourMomentum b1{31.8737, 27.0471, 11.7883, -12.0587};
    const yam2::FourMomentum b2{69.5370, 10.7832, 34.5499, -59.3751};

    // TransverseMomentum{px, py}
    const yam2::TransverseMomentum ptmiss{-14.1090, -93.3111};

    // M_{C1} = M_{C2}
    const yam2::Mass m_invisible{0.0};

    const auto input = yam2::mkInput({a1, a2}, {b1, b2}, ptmiss, m_invisible);
    if (!input) {
        std::cerr << "Invalid input.\n";
        return 1;
    }

    // tolerance is set to be 1.0e-3.
    // -- if not set, it will use the default value (1.0e-2).
    auto m2sol = yam2::m2CCSQP(input.value(), 1.0e-3);
    if (!m2sol) {
        std::cerr << "Failed to find minimum.\n";
        return 1;
    } else {
        // std::cout << m2sol.value() << '\n';

        // k1 and k2 are the M2 solutions to the momenta of C1 and C2.
        std::cout << "M2 = " << m2sol.value().m2() << '\n'
                  << "where \n"
                  << "  k1: " << m2sol.value().k1() << '\n'
                  << "  k2: " << m2sol.value().k2() << '\n';
    }
}
