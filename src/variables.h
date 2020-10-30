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
private:
    double k1x_, k1y_, k1z_, k2z_;

    explicit Variables(const NLoptVar &ks)
        : k1x_(ks[0]), k1y_(ks[1]), k1z_(ks[2]), k2z_(ks[3]) {}

public:
    Variables() = delete;

    double k1x() const { return k1x_; }
    double k1y() const { return k1y_; }
    double k1z() const { return k1z_; }
    double k2z() const { return k2z_; }
    NLoptVar variables() const { return {k1x_, k1y_, k1z_, k2z_}; }

    friend std::optional<Variables> mkVariables(const NLoptVar &ks);
};

inline std::optional<Variables> mkVariables(const NLoptVar &ks) {
    if (ks.size() != 4) {
        std::cerr << "mkVariables: invalid number of variables.\n";
        return {};
    }
    return Variables{ks};
}
}  // namespace yam2

#endif  // YAM2_SRC_VARIABLES_H_
