#ifndef GRIP_POINT_ND_H
#define GRIP_POINT_ND_H

#include <cmath>
#include <cstddef>
#include <vector>

namespace gripnd {

class PointND {
public:
    PointND() {}
    explicit PointND(std::size_t dim) : coord_(dim, 0.0) {}

    std::size_t dim() const { return coord_.size(); }

    double &operator[](std::size_t idx) { return coord_[idx]; }
    double operator[](std::size_t idx) const { return coord_[idx]; }

    void fill(double value)
    {
        for(double &x : coord_)
            x = value;
    }

    double norm2() const
    {
        double out = 0.0;
        for(double x : coord_)
            out += x * x;
        return out;
    }

    double norm() const
    {
        return std::sqrt(norm2());
    }

private:
    std::vector<double> coord_;
};

inline double squared_distance(const PointND &lhs, const PointND &rhs)
{
    double out = 0.0;
    for(std::size_t d = 0; d < lhs.dim(); d++){
        const double delta = lhs[d] - rhs[d];
        out += delta * delta;
    }
    return out;
}

} // namespace gripnd

#endif
