/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#ifndef YAM2_SRC_GRADIENT_H_
#define YAM2_SRC_GRADIENT_H_

#include <array>
#include <functional>
#include <ostream>
#include <tuple>
#include <utility>
#include "YAM2/input.h"       // InputKinematics
#include "YAM2/invisibles.h"  // Invisibles
#include "YAM2/momentum.h"    // FourMomentum
#include "YAM2/variables.h"   // Variables, NLoptVar

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

    // NLoptVar gradient() const { return {dk1x_, dk1y_, dk1z_, dk2z_}; }

    void set_gradient(NLoptVar &grad) const {
        grad[0] = dk1x_;
        grad[1] = dk1y_;
        grad[2] = dk1z_;
        if (grad.size() > 3) { grad[3] = dk2z_; }
    }

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

    Gradient operator-() const { return {-dk1x_, -dk1y_, -dk1z_, -dk2z_}; }

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

    friend std::ostream &operator<<(std::ostream &os, const Gradient &g);
};

using Gradients = std::pair<Gradient, Gradient>;

using GradFunc = std::function<std::pair<Gradient, double>(
    const InputKinematics &, const FourMomentum &, const Invisibles &)>;

std::pair<Gradient, double> m2Func1(const InputKinematics &inp,
                                    const FourMomentum &p1,
                                    const Invisibles &ks);

std::pair<Gradient, double> m2Func2(const InputKinematics &inp,
                                    const FourMomentum &p2,
                                    const Invisibles &ks);

std::tuple<Gradients, double, double> m2Func(const InputKinematics &inp,
                                             const FourMomentum &p1,
                                             const FourMomentum &p2,
                                             const Invisibles &ks);

Gradient mtotGrad(const InputKinematics &inp, const FourMomentum &p1,
                  const FourMomentum &p2, const Invisibles &ks, double m);
}  // namespace yam2

#endif  // YAM2_SRC_GRADIENT_H_
