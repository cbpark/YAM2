/*
 *  Copyright (c) 2020 Chan Beom Park <cbpark@gmail.com>
 */

#ifndef YAM2_SRC_VARIABLES_H_
#define YAM2_SRC_VARIABLES_H_

#include <optional>
#include <vector>

namespace yam2 {
class Variables {
private:
    double k1x_, k1y_, k1z_, k2z_;

    explicit Variables(const std::vector<double> &ks)
        : k1x_(ks[0]), k1y_(ks[1]), k1z_(ks[2]), k2z_(ks[3]) {}

public:
    Variables() = delete;

    double k1x() const { return k1x_; }
    double k1y() const { return k1y_; }
    double k1z() const { return k1z_; }
    double k2z() const { return k2z_; }
    std::vector<double> variables() const { return {k1x_, k1y_, k1z_, k2z_}; }

    friend std::optional<Variables> mkVariables(const std::vector<double> &ks);
};

inline std::optional<Variables> mkVariables(const std::vector<double> &ks) {
    if (ks.size() != 4) { return {}; }
    Variables var{ks};
    return var;
}
}  // namespace yam2

#endif  // YAM2_SRC_VARIABLES_H_
