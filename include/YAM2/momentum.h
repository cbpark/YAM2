/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#ifndef YAM2_SRC_MOMENTUM_H_
#define YAM2_SRC_MOMENTUM_H_

#include <cmath>  // std::sqrt, std::hypot, std::atan2
#include <ostream>
#include <vector>
#ifdef HAS_ROOT
#include "Math/LorentzVector.h"
#include "TLorentzVector.h"
#endif

namespace yam2 {
struct Mass {
    double value;
    Mass() : value(0) {}
    explicit Mass(double v) : value(v) {}

    double square() const { return value * value; }
    Mass operator*(double a) const { return Mass{value * a}; }
    Mass operator/(double a) const { return Mass{value / a}; }
};

class TransverseMomentum {
private:
    double x_, y_;

public:
    TransverseMomentum() = delete;
    TransverseMomentum(double x, double y) : x_(x), y_(y) {}

    double px() const { return x_; }
    double py() const { return y_; }
    double ptsq() const { return x_ * x_ + y_ * y_; }

    TransverseMomentum operator*(double a) const { return {a * x_, a * y_}; }

    TransverseMomentum operator/(double a) const { return {x_ / a, y_ / a}; }

    friend std::ostream &operator<<(std::ostream &os,
                                    const TransverseMomentum &p);
};

class SpatialMomentum {
private:
    double x_, y_, z_;

public:
    SpatialMomentum() = delete;
    SpatialMomentum(double x, double y, double z) : x_(x), y_(y), z_(z) {}

    double px() const { return x_; }
    double py() const { return y_; }
    double pz() const { return z_; }

    double p2() const { return x_ * x_ + y_ * y_ + z_ * z_; }

    double norm() const { return std::hypot(x_, y_, z_); }
    double p() const { return this->norm(); }

    double pt() const { return std::hypot(x_, y_); }

    double pt2() const { return x_ * x_ + y_ * y_; }

    double theta() const { return std::atan2(pt(), z_); }

    double phi() const { return std::atan2(y_, x_); }

    SpatialMomentum normalize() const {
        double norm = this->norm();
        if (norm < 1.0e-10) { norm = 1.0e-10; }
        return {x_ / norm, y_ / norm, z_ / norm};
    }

    double dot(const SpatialMomentum &v) const {
        return x_ * v.px() + y_ * v.py() + z_ * v.pz();
    }

    friend std::ostream &operator<<(std::ostream &os, const SpatialMomentum &p);
};

class FourMomentum {
private:
    double t_ = 0.0;
    double x_ = 0.0;
    double y_ = 0.0;
    double z_ = 0.0;

public:
    FourMomentum() = default;

    FourMomentum(double t, double x, double y, double z)
        : t_(t), x_(x), y_(y), z_(z) {}

    FourMomentum(Mass m, double x, double y, double z) : x_(x), y_(y), z_(z) {
        t_ = std::sqrt(x * x + y * y + z * z + m.square());
    }

    double e() const { return t_; }
    double px() const { return x_; }
    double py() const { return y_; }
    double pz() const { return z_; }

    double msq() const { return t_ * t_ - x_ * x_ - y_ * y_ - z_ * z_; }

    double m() const {
        const double mSq = msq();
        return mSq >= 0.0 ? std::sqrt(mSq) : std::sqrt(-mSq);
    }

    SpatialMomentum three_momentum() const { return {x_, y_, z_}; }

    FourMomentum &operator*=(double a) {
        this->t_ *= a;
        this->x_ *= a;
        this->y_ *= a;
        this->z_ *= a;
        return *this;
    }

    FourMomentum operator*(double a) const {
        return {a * t_, a * x_, a * y_, a * z_};
    }

    FourMomentum operator/(double a) const {
        return {t_ / a, x_ / a, y_ / a, z_ / a};
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

inline double invariantMass(const FourMomentum &p, const FourMomentum &k) {
    const auto psum = p + k;
    return psum.m();
}

inline FourMomentum sum(const std::vector<FourMomentum> &ps) {
    FourMomentum psum{0.0, 0.0, 0.0, 0.0};
    for (const auto &p : ps) { psum += p; }
    return psum;
}

#ifdef HAS_ROOT
inline ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> toLorentzVector(
    const FourMomentum &p) {
    return {p.px(), p.py(), p.pz(), p.e()};
}

inline TLorentzVector toTLorentzVector(const FourMomentum &p) {
    return {p.px(), p.py(), p.pz(), p.e()};
}
#endif
}  // namespace yam2

#endif  // YAM2_SRC_MOMENTUM_H_
