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

class TransverseMomentum {
private:
    double x_, y_;

public:
    TransverseMomentum() = delete;
    TransverseMomentum(double x, double y) : x_(x), y_(y) {}

    double px() const { return x_; }
    double py() const { return y_; }
    double pt() const { return std::hypot(x_, y_); }

    TransverseMomentum operator*(double a) const { return {a * x_, a * y_}; }

    friend std::ostream &operator<<(std::ostream &os,
                                    const TransverseMomentum &p);
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

    double pt() const { return std::sqrt(x_ * x_ + y_ * y_); }

    double transverseEnergy() const { return std::sqrt(t_ * t_ - z_ * z_); }

    TransverseMomentum transverseVector() const { return {x_, y_}; }

    double dot(const FourMomentum &p) const {
        return t_ * p.e() - x_ * p.px() - y_ * p.py() - z_ * p.pz();
    }

    FourMomentum operator*(double a) const {
        return {a * t_, a * x_, a * y_, a * z_};
    }

    FourMomentum &operator+=(const FourMomentum &p) {
        this->t_ += p.e();
        this->x_ += p.px();
        this->y_ += p.py();
        this->z_ += p.pz();
        return *this;
    }

    friend FourMomentum operator+(FourMomentum p1, const FourMomentum &p2) {
        p1 += p2;
        return p1;
    }

    friend std::ostream &operator<<(std::ostream &os, const FourMomentum &p);
};

template <template <typename> class F>
FourMomentum sum(const F<FourMomentum> &ps) {
    FourMomentum psum{0, 0, 0, 0};
    for (const auto &p : ps) { psum += p; }
    return psum;
}
}  // namespace yam2

#endif  // YAM2_SRC_MOMENTUM_H_
