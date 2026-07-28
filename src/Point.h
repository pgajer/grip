#ifndef POINT_HPP
#define POINT_HPP

#include <cmath>
#include <cstddef>
#include <algorithm>

using coord_t = double;

template <class T = coord_t>
class Point {
public:
    // Constructors
    Point();
    Point(T x, T y);
    Point(const T* _coord, size_t _dim);
    Point(T x, T y, T z);
    Point(const Point<T>& init);
    ~Point() { delete[] coord; }

    // Other methods
    double norm() const;
    double fnorm() const;
    double fnorm2() const;
    double norm2() const;
    T scal_pr(const Point<T>& u, const Point<T>& v, const Point<T>& w) const;

    // Accessors
    T getX() const { return coord[0]; }
    T getY() const { return coord[1]; }
    T getZ() const { return dim > 2 ? coord[2] : T(); }
    size_t getDim() { return dim; }

    // Modifiers
    void setX(T x) { coord[0] = x; }
    void setY(T y) { coord[1] = y; }
    void setZ(T z) { if (dim > 2) coord[2] = z; }
    void set_to_zero();

    // Assignment operator
    Point<T>& operator=(const Point<T>& rhs);

    // Comparison
    bool operator<(const Point<T>& p) const;
    bool operator==(const Point<T>& right) const;

    // Constant operators.
    template<typename U>
    Point<T> operator*(U constant) const {
        Point<T> result(*this);
        for (size_t i = 0; i < dim; ++i) {
            result.coord[i] *= constant;
        }
        return result;
    }

    template<typename U>
    Point<T> operator/(U constant) const {
        Point<T> result(*this);
        for (size_t i = 0; i < dim; ++i) {
            result.coord[i] /= constant;
        }
        return result;
    }

    template<typename U>
    Point<T> operator+(U constant) const {
        Point<T> result(*this);
        for (size_t i = 0; i < dim; ++i) {
            result.coord[i] += constant;
        }
        return result;
    }

    template<typename U>
    Point<T> operator-(U constant) const {
        Point<T> result(*this);
        for (size_t i = 0; i < dim; ++i) {
            result.coord[i] -= constant;
        }
        return result;
    }

    //
    // (op)= operators
    //
    template<typename U>
    Point<T>& operator*=(U constant) {
        for (size_t i = 0; i < dim; ++i) {
            coord[i] *= constant;
        }
        return *this;
    }

    template<typename U>
    Point<T>& operator/=(U constant) {
        for (size_t i = 0; i < dim; ++i) {
            coord[i] /= constant;
        }
        return *this;
    }

    template<typename U>
    Point<T>& operator+=(U constant) {
        for (size_t i = 0; i < dim; ++i) {
            coord[i] += constant;
        }
        return *this;
    }

    template<typename U>
    Point<T>& operator-=(U constant) {
        for (size_t i = 0; i < dim; ++i) {
            coord[i] -= constant;
        }
        return *this;
    }

    //
    // this (op) other_point   operators
    //
    Point<T> operator+(const Point<T>& p) const {
        Point<T> result(*this);
        for (size_t i = 0; i < dim; ++i) {
            result.coord[i] += p.coord[i];
        }
        return result;
    }

    Point<T> operator-(const Point<T>& p) const {
        Point<T> result(*this);
        for (size_t i = 0; i < dim; ++i) {
            result.coord[i] -= p.coord[i];
        }
        return result;
    }

    // scalar product operator
    T operator*(const Point<T>& p) const {
        T result = 0;
        for (size_t i = 0; i < dim; ++i) {
            result += coord[i] * p.coord[i];
        }
        return result;
    }

    // negation operator
    Point<T> operator-() const {
        return (*this) * static_cast<T>(-1);
    }

    // this (op)= other_point   operators
    Point<T>& operator+=(const Point<T>& p) {
        for (size_t i = 0; i < dim; ++i) {
            coord[i] += p.coord[i];
        }
        return *this;
    }

    Point<T>& operator-=(const Point<T>& p) {
        for (size_t i = 0; i < dim; ++i) {
            coord[i] -= p.coord[i];
        }
        return *this;
    }


    // Friend functions
    template<typename U>
    friend Point<T> operator*(U constant, const Point<T>& p) {
        return p * constant;
    }

    template<typename U>
    friend Point<T> operator+(U constant, const Point<T>& p) {
        return p + constant;
    }

private:
    size_t dim;
    T* coord;
};

// Implementation of member functions

template<class T>
Point<T>::Point() : dim(3), coord(new T[3]()) {}

template<class T>
Point<T>::Point(T x, T y) : dim(2), coord(new T[2]{x, y}) {}

template<class T>
Point<T>::Point(const T* _coord, size_t _dim) : dim(_dim), coord(new T[_dim]) {
    std::copy(_coord, _coord + dim, coord);
}

template<class T>
Point<T>::Point(T x, T y, T z) : dim(3), coord(new T[3]{x, y, z}) {}

template<class T>
Point<T>::Point(const Point<T>& init) : dim(init.dim), coord(new T[init.dim]) {
    std::copy(init.coord, init.coord + dim, coord);
}

template<class T>
Point<T>& Point<T>::operator=(const Point<T>& rhs) {
    if (this != &rhs) {
        if (dim != rhs.dim) {
            delete[] coord;
            dim = rhs.dim;
            coord = new T[dim];
        }
        std::copy(rhs.coord, rhs.coord + dim, coord);
    }
    return *this;
}

template<class T>
bool Point<T>::operator<(const Point<T> &p) const {
    for(size_t ind = 0; ind < std::min(dim, p.dim); ind++) {
        if(coord[ind] < p.coord[ind]) return true;
        if(coord[ind] > p.coord[ind]) return false;
    }
    return dim < p.dim;
}

template<class T>
bool Point<T>::operator==(const Point<T> &right) const {
    if(dim != right.dim) return false;
    for(size_t ind = 0; ind < dim; ind++)
        if(coord[ind] != right.coord[ind]) return false;
    return true;
}

template<class T>
double Point<T>::norm() const {
    return std::sqrt(norm2());
}

template<class T>
double Point<T>::fnorm() const {
    return std::sqrt(norm2());
}

template<class T>
double Point<T>::norm2() const {
    double result = 0.0;
    for(size_t ind = 0; ind < dim; ind++){
        result += static_cast<double>(coord[ind]) * static_cast<double>(coord[ind]);
    }
    return result;
}

template<class T>
double Point<T>::fnorm2() const {
    return norm2();
}

template<class T>
T Point<T>::scal_pr(const Point<T> &u, const Point<T> &v, const Point<T> &w) const {
    T uv_x = v.getX() - u.getX();
    T uv_y = v.getY() - u.getY();
    T uv_z = v.getZ() - u.getZ();

    T uw_x = w.getX() - u.getX();
    T uw_y = w.getY() - u.getY();
    T uw_z = w.getZ() - u.getZ();

    return uv_x * uw_x + uv_y * uw_y + uv_z * uw_z;
}

template<class T>
void Point<T>::set_to_zero() {
    std::fill(coord, coord + dim, T());
}

template<typename T, typename U>
Point<T> operator*(U constant, const Point<T>& p) {
    return p * constant;
}

#endif
