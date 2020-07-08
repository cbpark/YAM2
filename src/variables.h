#ifndef YAM2_SRC_VARIABLES_H_
#define YAM2_SRC_VARIABLES_H_

namespace yam2 {
class Variables {
private:
    double k1x_, k1y_, k1z_, k2z_;

public:
    Variables() = delete;
    Variables(double k1x, double k1y, double k1z, double k2z)
        : k1x_(k1x), k1y_(k1y), k1z_(k1z), k2z_(k2z) {}

    double k1x() const { return k1x_; }
    double k1y() const { return k1y_; }
    double k1z() const { return k1z_; }
    double k2z() const { return k2z_; }
};
}  // namespace yam2

#endif  // YAM2_SRC_VARIABLES_H_
