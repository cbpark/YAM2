#include "yam2.h"

#include "constraint.h"
#include "input.h"

#include <functional>
#include <optional>

namespace yam2 {
std::optional<M2Solution> m2SQP(const Constraints &cfs,
                                const std::optional<InputKinematics> &inp) {
    if (!inp) { return {}; }


}
}  // namespace yam2
