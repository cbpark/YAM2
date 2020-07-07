#ifndef YAM2_SRC_MOMENTUM_H_
#define YAM2_SRC_MOMENTUM_H_

#include <cmath>
#include <ostream>

namespace yam2 {
struct Mass {
    double value;
    Mass() : value(0) {}
    explicit Mass(double v) : value(v) {}
    double square() const { return value * value; }
};

class FourMomentum {
private:
    double t_, x_, y_, z_;

public:
    FourMomentum() = delete;
    FourMomentum(double t, double x, double y, double z)
        : t_(t), x_(x), y_(y), z_(z) {}

    double e() const { return t_; }
    double px() const { return x_; }
    double py() const { return y_; }
    double pz() const { return z_; }

    double m2() const { return t_ * t_ - x_ * x_ - y_ * y_ - z_ * z_; }

    double m() const {
        const double mSq = m2();
        return mSq >= 0 ? std::sqrt(mSq) : 0;
    }

    double transverseEnergy() const {
        return std::sqrt(x_ * x_ + y_ * y_ + m2());
    }

    double dot(const FourMomentum &p) const {
        return e() * p.e() - px() * p.px() - py() * p.py() - pz() * p.pz();
    }

    friend std::ostream &operator<<(std::ostream &os, const FourMomentum &p);
};

std::ostream &operator<<(std::ostream &os, const FourMomentum &p);
}  // namespace yam2

#endif  // YAM2_SRC_MOMENTUM_H_
