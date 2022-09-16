#pragma once

#include <ostream>

template <typename T>
struct Vec3 {
    T d[3];

    Vec3(const T u, const T v, const T w) : d{u, v, w} {}

    Vec3(const T a[3]) : d{a[0], a[1], a[2]} {}

    Vec3() : d{0, 0, 0} {}

    T& operator[](int i) { return d[i]; }

    T operator()(int i) const { return d[i]; }

    Vec3<T>& operator=(double s) {
        d[0] = s;
        d[1] = s;
        d[2] = s;
        return (*this);
    }

    Vec3<T>& operator+=(Vec3<T> o) {
        d[0] += o[0];
        d[1] += o[1];
        d[2] += o[2];
        return (*this);
    }

    Vec3<T>& operator-=(Vec3<T> o) {
        d[0] -= o[0];
        d[1] -= o[1];
        d[2] -= o[2];
        return (*this);
    }
};

// vec3-vec3 operations
template <typename T>  // addition of two vec3s
Vec3<T> operator+(const Vec3<T>& a, const Vec3<T>& b) {
    return Vec3<T>(a(0) + b(0), a(1) + b(1), a(2) + b(2));
}
template <typename T>  // subtraction of two vec3s
Vec3<T> operator-(const Vec3<T>& a, const Vec3<T>& b) {
    return Vec3<T>(a(0) - b(0), a(1) - b(1), a(2) - b(2));
}
template <typename T>  // element-wise multiplication of two vec3s
Vec3<T> operator*(const Vec3<T>& a, const Vec3<T>& b) {
    return Vec3<T>(a(0) * b(0), a(1) * b(1), a(2) * b(2));
}
template <typename T>  // element wise division of two vec3s
Vec3<T> operator/(const Vec3<T>& a, const Vec3<T>& b) {
    return Vec3<T>(a(0) / b(0), a(1) / b(1), a(2) / b(2));
}

// vec3 - scalar operations
template <typename T>  // scalar multiplication
Vec3<T> operator*(const Vec3<T>& a, T s) {
    return Vec3<T>(a(0) * s, a(1) * s, a(2) * s);
}
template <typename T>  // scalar multiplication 2
Vec3<T> operator*(T s, const Vec3<T>& a) {
    return Vec3<T>(a(0) * s, a(1) * s, a(2) * s);
}

// output
template <typename T>  // ostream output
std::ostream& operator<<(std::ostream& out, Vec3<T>& v) {
    out << v[0] << " " << v[1] << " " << v[2];
    return out;
}

using Vec3d = Vec3<double>;
using Vec3i = Vec3<int>;
