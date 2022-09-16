#pragma once

#include <ostream>

template <typename T>
struct Vec2 {
    T d[2];

    Vec2(const T u, const T v) : d{u, v} {}

    Vec2(const T a[2]) : d{a[0], a[1]} {}

    Vec2() : d{0, 0} {}

    T& operator[](int i) { return d[i]; }

    T operator()(int i) const { return d[i]; }

    Vec2<T>& operator=(double s) {
        d[0] = s;
        d[1] = s;
        return (*this);
    }

    Vec2<T>& operator+=(Vec2<T> o) {
        d[0] += o[0];
        d[1] += o[1];
        return (*this);
    }

    Vec2<T>& operator-=(Vec2<T> o) {
        d[0] -= o[0];
        d[1] -= o[1];
        return (*this);
    }
};

// vec2-vec2 operations
template <typename T>  // addition of two vec3s
Vec2<T> operator+(const Vec2<T>& a, const Vec2<T>& b) {
    return Vec2<T>(a(0) + b(0), a(1) + b(1));
}
template <typename T>  // subtraction of two vec3s
Vec2<T> operator-(const Vec2<T>& a, const Vec2<T>& b) {
    return Vec2<T>(a(0) - b(0), a(1) - b(1));
}
template <typename T>  // element-wise multiplication of two vec3s
Vec2<T> operator*(const Vec2<T>& a, const Vec2<T>& b) {
    return Vec2<T>(a(0) * b(0), a(1) * b(1));
}
template <typename T>  // element wise division of two vec3s
Vec2<T> operator/(const Vec2<T>& a, const Vec2<T>& b) {
    return Vec2<T>(a(0) / b(0), a(1) / b(1));
}

// vec3 - scalar operations
template <typename T>  // scalar multiplication
Vec2<T> operator*(const Vec2<T>& a, T s) {
    return Vec2<T>(a(0) * s, a(1) * s);
}
template <typename T>  // scalar multiplication 2
Vec2<T> operator*(T s, const Vec2<T>& a) {
    return Vec2<T>(a(0) * s, a(1) * s);
}

// output
template <typename T>  // ostream output
std::ostream& operator<<(std::ostream& out, Vec2<T>& v) {
    out << v[0] << " " << v[1];
    return out;
}

using Vec2d = Vec2<double>;
using Vec2i = Vec2<int>;