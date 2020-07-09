#ifndef YAM2_SRC_GRADIENT_H_
#define YAM2_SRC_GRADIENT_H_

#include "input.h"
#include "invisibles.h"
#include "momentum.h"
#include "variables.h"

#include <array>
#include <tuple>

namespace yam2 {
class Gradient {
private:
    double dk1x_, dk1y_, dk1z_, dk2z_;

public:
    Gradient() = delete;
    Gradient(double dk1x, double dk1y, double dk1z, double dk2z)
        : dk1x_(dk1x), dk1y_(dk1y), dk1z_(dk1z), dk2z_(dk2z) {}

    double dk1x() const { return dk1x_; }
    double dk1y() const { return dk1y_; }
    double dk1z() const { return dk1z_; }
    double dk2z() const { return dk2z_; }

    Gradient &operator*=(double a) {
        this->dk1x_ *= a;
        this->dk1y_ *= a;
        this->dk1z_ *= a;
        this->dk2z_ *= a;
        return *this;
    }

    Gradient operator*(double a) {
        return {a * dk1x_, a * dk1y_, a * dk1z_, a * dk2z_};
    }

    Gradient &operator-=(const Gradient &g) {
        this->dk1x_ -= g.dk1x();
        this->dk1y_ -= g.dk1y();
        this->dk1z_ -= g.dk1z();
        this->dk2z_ -= g.dk2z();
        return *this;
    }

    friend Gradient operator-(Gradient g1, const Gradient &g2) {
        g1 -= g2;
        return g1;
    }
};

using FGType = std::pair<double, Gradient>;

using Gradients = std::pair<Gradient, Gradient>;

std::tuple<Gradients, double, double> m2Grad(const InputKinematics &inp,
                                             const Invisibles &ks,
                                             const FourMomentum &p1,
                                             const FourMomentum &p2,
                                             const Variables &var);
}  // namespace yam2

#endif  // YAM2_SRC_GRADIENT_H_
