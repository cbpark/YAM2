/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#ifndef YAM2_SRC_INPUT_H_
#define YAM2_SRC_INPUT_H_

#include <optional>
#include <ostream>
#include <vector>

#include "momentum.h"
#include "variables.h"

namespace yam2 {
class InputKinematics {
private:
    /** p1 = a1 + b1 */
    FourMomentum p1_;
    /** p2 = a2 + b2 */
    FourMomentum p2_;
    /** q1 = b1 */
    FourMomentum q1_;
    /** q2 = b2 */
    FourMomentum q2_;
    /** missing transverse momentum */
    TransverseMomentum ptmiss_;
    /** the input mass for the invisible particle */
    Mass minv_;
    /** the input mass for the parent particle */
    std::optional<Mass> mparent_;
    /** the input mass for the relative particle */
    std::optional<Mass> mrel_;
    /** collision energy (sqrt(s)) */
    double sqrt_s_;
    /** the longitudinal momentum of the total system */
    std::optional<double> ptot_z_;
    /** an energy scale of the process */
    double scale_;
    /** the precision of the equality constraint of the on-shell mass constraint
     * for the parent particle masses. It can be set by set_eps_constraint.
     */
    double eps_constraint_ = 1.0e-2;

protected:
    InputKinematics(const FourMomentum &p1, const FourMomentum &p2,
                    const FourMomentum &q1, const FourMomentum &q2,
                    const TransverseMomentum &ptmiss, const Mass &minv,
                    const std::optional<Mass> &mparent,
                    const std::optional<Mass> &mrel, double sqrt_s,
                    const std::optional<double> ptot_z, double scale)
        : p1_(p1),
          p2_(p2),
          q1_(q1),
          q2_(q2),
          ptmiss_(ptmiss),
          minv_(minv),
          mparent_(mparent),
          mrel_(mrel),
          sqrt_s_(sqrt_s),
          ptot_z_(ptot_z),
          scale_(scale) {}

public:
    InputKinematics() = delete;
    virtual ~InputKinematics() {}

    FourMomentum p1() const { return p1_; }
    FourMomentum p2() const { return p2_; }
    FourMomentum q1() const { return q1_; }
    FourMomentum q2() const { return q2_; }

    TransverseMomentum ptmiss() const { return ptmiss_; }

    Mass minv() const { return minv_; }
    double minv_square() const { return minv_.square(); }

    std::optional<Mass> mparent() const { return mparent_; }

    std::optional<Mass> mrel() const { return mrel_; }

    double sqrt_s() const { return sqrt_s_; }

    std::optional<double> ptot_z() const { return ptot_z_; }

    double scale() const { return scale_; }

    double eps_constraint() const { return eps_constraint_; }
    void set_eps_constraint(double eps) { eps_constraint_ = eps; }

    /** the initial guess configuration */
    NLoptVar initial_guess(double eps, unsigned int neval);

    virtual void show(std::ostream &os) const;

    friend std::optional<InputKinematics> mkInput(
        const std::vector<FourMomentum> &as,
        const std::vector<FourMomentum> &bs, const TransverseMomentum &ptmiss,
        const Mass &minv, const std::optional<Mass> &mparent,
        const std::optional<Mass> &mrel, double sqrt_s,
        const std::optional<double> ptot_z);

    friend std::ostream &operator<<(std::ostream &os, const InputKinematics &p);
};

std::optional<InputKinematics> mkInput(
    const std::vector<FourMomentum> &as, const std::vector<FourMomentum> &bs,
    const TransverseMomentum &ptmiss, const Mass &minv,
    const std::optional<Mass> &mparent = {Mass{0.0}},
    const std::optional<Mass> &mrel = {Mass{0.0}}, double sqrt_s = 0.0,
    const std::optional<double> ptot_z = {});

class InputKinematicsWithVertex : public InputKinematics {
private:
    SpatialMomentum vertex1_;
    SpatialMomentum vertex2_;
    double delta_theta_;

    InputKinematicsWithVertex(const FourMomentum &p1, const FourMomentum &p2,
                              const FourMomentum &q1, const FourMomentum &q2,
                              const TransverseMomentum &ptmiss,
                              const Mass &minv, const SpatialMomentum &vertex1,
                              const SpatialMomentum &vertex2,
                              double delta_theta,
                              const std::optional<Mass> &mparent,
                              const std::optional<Mass> &mrel, double sqrt_s,
                              const std::optional<double> ptot_z, double scale)
        : InputKinematics(p1, p2, q1, q2, ptmiss, minv, mparent, mrel, sqrt_s,
                          ptot_z, scale),
          vertex1_(vertex1),
          vertex2_(vertex2),
          delta_theta_(delta_theta) {}

public:
    InputKinematicsWithVertex() = delete;
    virtual ~InputKinematicsWithVertex() {}

    SpatialMomentum vertex1() const { return vertex1_; }
    SpatialMomentum vertex2() const { return vertex2_; }
    double delta_theta() const { return delta_theta_; }

    virtual void show(std::ostream &os) const;

    friend std::optional<InputKinematicsWithVertex> mkInputWithVertex(
        const std::vector<FourMomentum> &as,
        const std::vector<FourMomentum> &bs, const TransverseMomentum &ptmiss,
        const Mass &minv, const SpatialMomentum &vertex1,
        const SpatialMomentum &vertex2, double delta_theta,
        const std::optional<Mass> &mparent, const std::optional<Mass> &mrel,
        double sqrt_s, const std::optional<double> ptot_z);

    friend std::optional<InputKinematicsWithVertex> mkInputWithVertex(
        const std::optional<InputKinematics> &input_kinematics,
        const SpatialMomentum &vertex1, const SpatialMomentum &vertex2,
        double delta_theta);
};

std::optional<InputKinematicsWithVertex> mkInputWithVertex(
    const std::vector<FourMomentum> &as, const std::vector<FourMomentum> &bs,
    const TransverseMomentum &ptmiss, const Mass &minv,
    const SpatialMomentum &vertex1, const SpatialMomentum &vertex2,
    double delta_theta = 0.0, const std::optional<Mass> &mparent = {Mass{0.0}},
    const std::optional<Mass> &mrel = {Mass{0.0}}, double sqrt_s = 0.0,
    const std::optional<double> ptot_z = {});

std::optional<InputKinematicsWithVertex> mkInputWithVertex(
    const std::optional<InputKinematics> &input_kinematics,
    const SpatialMomentum &vertex1, const SpatialMomentum &vertex2,
    double delta_theta);

/**
 *  the difference between the total invariant mass and the collision energy
 */
double deltaSqrtS(const NLoptVar &x, NLoptVar &grad, void *input);

template <typename T>
std::optional<T> scaleIfExists(const std::optional<T> &x, const double scale) {
    if (x) {
        return {x.value() / scale};
    } else {
        return x;
    }
}
}  // namespace yam2

#endif  // YAM2_SRC_INPUT_H_
