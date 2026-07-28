#ifndef GRIP_ROUNDING_H
#define GRIP_ROUNDING_H

namespace grip {
namespace detail {

template <typename Integer>
inline Integer round_to_integer(double value) noexcept {
    return value > 0.0
        ? static_cast<Integer>(value + 0.5)
        : -static_cast<Integer>(0.5 - value);
}

} // namespace detail
} // namespace grip

#endif
