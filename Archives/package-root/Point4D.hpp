#ifndef POINT4D_HPP
#define POINT4D_HPP

#include <array>
#include <algorithm>  // for std::copy
#include <ranges>     // if using std::ranges::copy (C++20)
#include <iostream>
#include <cmath>

#include "Point.hpp"

#define ROUND(a) ((a)>0 ? (int)((a)+0.5) : -(int)(0.5-(a)))
#define ROUND_L(a) ((a)>0 ? (unsigned long)((a)+0.5) : -(unsigned long)(0.5-(a)))

template<typename T = double>
class Point4D {
public:
    // Constructors
    Point4D() : coord{0, 0, 0, 0} {}
    Point4D(T x, T y, T z, T w) : coord{x, y, z, w} {}

    Point4D(const T* _coord, size_t _dim);
    Point4D(const Point4D<T>& init);
    ~Point4D() {}

    // Add begin() and end() methods
    T* begin() { return coord.data(); }
    T* end() { return coord.data() + 4; }
    const T* begin() const { return coord.data(); }
    const T* end() const { return coord.data() + 4; }

    // Other methods
    unsigned long norm() const;
    float fnorm() const;
    float fnorm2() const;
    unsigned long norm2() const;
    T scal_pr(const Point4D<T>& u, const Point4D<T>& v, const Point4D<T>& w) const;

    // Accessors
    T getX() const { return coord[0]; }
    T getY() const { return coord[1]; }
    T getZ() const { return coord[2]; }
    T getW() const { return coord[3]; }
     size_t getDim() { return dim; }

    // Modifiers
    void setX(T x) { coord[0] = x; }
    void setY(T y) { coord[1] = y; }
    void setZ(T z) { coord[2] = z; }
    void setW(T w) { coord[3] = w; }
    void set_to_zero();

    // [] operator
    T& operator[](int i) { return coord[i]; }
    const T& operator[](int i) const { return coord[i]; }

    // Assignment operator
    Point4D<T>& operator=(const Point4D<T>& rhs);

    // Comparison
    bool operator<(const Point4D<T>& p) const;
    bool operator==(const Point4D<T>& right) const;

    //
    // this (op) constant  oparators
    //
    template<typename U>
    Point4D<T> operator*(U constant) const {
        return Point4D<T>(coord[0] * constant, coord[1] * constant, coord[2] * constant, coord[3] * constant);
    }

    template<typename U>
    Point4D<T> operator/(U constant) const {
        return Point4D<T>(coord[0] / constant, coord[1] / constant, coord[2] / constant, coord[3] / constant);
    }

    template<typename U>
    Point4D<T> operator+(U constant) const {
        return Point4D<T>(coord[0] + constant, coord[1] + constant, coord[2] + constant, coord[3] + constant);
    }

    template<typename U>
    Point4D<T> operator-(U constant) const {
        return Point4D<T>(coord[0] - constant, coord[1] - constant, coord[2] - constant, coord[3] - constant);
    }

    //
    // (op)= constant   operators
    //
    template<typename U>
    Point4D<T>& operator*=(U constant) {
        coord[0] *= constant;
        coord[1] *= constant;
        coord[2] *= constant;
        coord[3] *= constant;
        return *this;
    }

    template<typename U>
    Point4D<T>& operator/=(U constant) {
        coord[0] /= constant;
        coord[1] /= constant;
        coord[2] /= constant;
        coord[3] /= constant;
        return *this;
    }

    template<typename U>
    Point4D<T>& operator+=(U constant) {
        coord[0] += constant;
        coord[1] += constant;
        coord[2] += constant;
        coord[3] += constant;
        return *this;
    }

    template<typename U>
    Point4D<T>& operator-=(U constant) {
        coord[0] -= constant;
        coord[1] -= constant;
        coord[2] -= constant;
        coord[3] -= constant;
        return *this;
    }

    //
    // this (op) other_point   operators
    //
    Point4D<T> operator+(const Point4D<T>& p) const {
        Point4D<T> result(*this);
        for (size_t i = 0; i < dim; ++i) {
            result.coord[i] += p.coord[i];
        }
        return result;
    }

    Point4D<T> operator-(const Point4D<T>& p) const {
        Point4D<T> result(*this);
        for (size_t i = 0; i < dim; ++i) {
            result.coord[i] -= p.coord[i];
        }
        return result;
    }

    // scalar product operator
    T operator*(const Point4D<T>& p) const {
        T result = 0;
        for (size_t i = 0; i < dim; ++i) {
            result += coord[i] * p.coord[i];
        }
        return result;
    }

    // negation operator
    Point4D<T> operator-() const {
        return (*this) * static_cast<T>(-1);
    }

    // this (op)= other_point   operators
    Point4D<T>& operator+=(const Point4D<T>& p) {
        for (size_t i = 0; i < dim; ++i) {
            coord[i] += p.coord[i];
        }
        return *this;
    }

    Point4D<T>& operator-=(const Point4D<T>& p) {
        for (size_t i = 0; i < dim; ++i) {
            coord[i] -= p.coord[i];
        }
        return *this;
    }

    // Friend functions
    friend std::ostream& operator<<(std::ostream& output, const Point4D<T>& p);

    template<typename U>
    friend Point4D<T> operator*(U constant, const Point4D<T>& p) {
        return p * constant;
    }

    template<typename U>
    friend Point4D<T> operator+(U constant, const Point4D<T>& p) {
        return p + constant;
    }

private:
    size_t dim = 4;
    std::array<T, 4> coord;
};

template<class T>
T Point4D<T>::scal_pr(const Point4D<T>& u, const Point4D<T>& v, const Point4D<T>& w) const {
    T uv_x = v.getX() - u.getX();
    T uv_y = v.getY() - u.getY();
    T uv_z = v.getZ() - u.getZ();
    T uv_w = v.getW() - u.getW();

    T uw_x = w.getX() - u.getX();
    T uw_y = w.getY() - u.getY();
    T uw_z = w.getZ() - u.getZ();
    T uw_w = w.getW() - u.getW();

    return uv_x * uw_x + uv_y * uw_y + uv_z * uw_z + uv_w * uw_w;
}

// Implementation of member functions

template<class T>
Point4D<T>::Point4D(const T* _coord, size_t _dim) {
    size_t copy_size = std::min(_dim, dim);
    std::copy(_coord, _coord + copy_size, coord.begin());
    if (_dim < dim) {
        std::fill(coord.begin() + _dim, coord.end(), T());
    }
}

template<class T>
Point4D<T>::Point4D(const Point4D<T>& init) {
    std::copy(init.coord.begin(), init.coord.end(), coord.begin());
}

template<class T>
Point4D<T>& Point4D<T>::operator=(const Point4D<T>& rhs) {
    if (this != &rhs) {
        std::copy(rhs.coord.begin(), rhs.coord.end(), coord.begin());
    }
    return *this;
}

template<class T>
bool Point4D<T>::operator<(const Point4D<T>& p) const {
    return std::lexicographical_compare(coord, coord + dim, p.coord, p.coord + dim);
}

template<class T>
bool Point4D<T>::operator==(const Point4D<T>& right) const {
    return std::equal(coord, coord + dim, right.coord);
}

template<class T>
unsigned long Point4D<T>::norm() const {
    return ROUND_L(std::sqrt(static_cast<double>(norm2())));
}

template<class T>
float Point4D<T>::fnorm() const {
    return std::sqrt(fnorm2());
}

template<class T>
unsigned long Point4D<T>::norm2() const {
    unsigned long result = 0;
    for (size_t i = 0; i < dim; ++i) {
        result += static_cast<unsigned long>(coord[i]) * static_cast<unsigned long>(coord[i]);
    }
    return result;
}

template<class T>
float Point4D<T>::fnorm2() const {
    float result = 0;
    for (size_t i = 0; i < dim; ++i) {
        result += static_cast<float>(coord[i]) * static_cast<float>(coord[i]);
    }
    return result;
}

template<class T>
void Point4D<T>::set_to_zero() {
    std::fill(coord.begin(), coord.end(), T());
}

// Friend function implementations

template<typename T>
std::ostream& operator<<(std::ostream& output, const Point4D<T>& p) {
    output << '(';
    for (size_t i = 0; i < Point4D<T>::dim - 1; ++i) {
        output << p.coord[i] << ", ";
    }
    output << p.coord[Point4D<T>::dim - 1] << ')';
    return output;
}

#endif
