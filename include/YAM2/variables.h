/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#ifndef YAM2_SRC_VARIABLES_H_
#define YAM2_SRC_VARIABLES_H_

#include <iostream>
#include <optional>
#include <vector>

namespace yam2 {
using NLoptVar = std::vector<double>;

class Variables {
public:
    enum class Dim { THREE, FOUR };

private:
    double k1x_, k1y_, k1z_, k2z_;
    Variables::Dim dim_;

    Variables(const double k1x, const double k1y, const double k1z,
              const double k2z, const Variables::Dim dim)
        : k1x_(k1x), k1y_(k1y), k1z_(k1z), k2z_(k2z), dim_(dim) {}

public:
    Variables() = delete;

    double k1x() const { return k1x_; }
    double k1y() const { return k1y_; }
    double k1z() const { return k1z_; }
    double k2z() const { return k2z_; }

    Variables::Dim dimension() const { return dim_; }

    friend std::optional<Variables> mkVariables(const NLoptVar &ks);
};

inline std::optional<Variables> mkVariables(const NLoptVar &ks) {
    switch (ks.size()) {
    case 3:
        return {{ks[0], ks[1], ks[2], 0.0, Variables::Dim::THREE}};

    case 4:
        return {{ks[0], ks[1], ks[2], ks[3], Variables::Dim::FOUR}};

    default:
        std::cerr << "mkVariables: invalid number of variables.\n";
        return {};
    }
}
}  // namespace yam2

#endif  // YAM2_SRC_VARIABLES_H_
