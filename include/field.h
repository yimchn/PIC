#pragma once

#include <string>

#include "vector.h"

template <typename T>
struct Field {
    int ni, nj;  // allocated dimensions
    T **data;    /*data held by this field*/

    /*constructor*/
    Field(int ni, int nj) : ni{ni}, nj{nj} {
        // allocate memory for a 3D array
        data = new T *[ni];
        for (int i = 0; i < ni; i++) {
            data[i] = new T[nj];
        }

        clear();
    }

    // copy constructor
    Field(const Field &other) : Field{other.ni, other.nj} {
        for (int i = 0; i < ni; i++)
            for (int j = 0; j < nj; j++) {
                data[i][j] = other(i, j);
            }
    }

    // move constructor
    Field(Field &&other) : ni{other.ni}, nj{other.nj} {
        data = other.data;     // steal the data
        other.data = nullptr;  // invalidate
    }

    // move assignment operator
    Field &operator=(Field &&f) { return *this; }

    // destructor: release memory
    ~Field() {
        // don't do anything if data is not allocated (or was moved away)
        if (data == nullptr) {
            return;
        }

        for (int i = 0; i < ni; i++) {
            delete[] data[i];
        }

        delete[] data;
    }

    // overloaded operator [] to allow direct access to data
    T *operator[](int i) { return data[i]; }

    /*returns data[i][j][k] marked as const to signal no data change*/
    T operator()(int i, int j) const { return data[i][j]; }

    /*sets all values to some scalar*/
    void operator=(double s) {
        for (int i = 0; i < ni; i++)
            for (int j = 0; j < nj; j++) {
                data[i][j] = s;
            }
    }

    /*performs element by element division by another field*/
    void operator/=(const Field &other) {
        for (int i = 0; i < ni; i++)
            for (int j = 0; j < nj; j++) {
                if (other.data[i][j] != 0) {
                    data[i][j] /= other.data[i][j];
                } else {
                    data[i][j] = 0;
                }
            }
    }

    /*increments values by data from another field*/
    Field &operator+=(const Field &other) {
        for (int i = 0; i < ni; i++)
            for (int j = 0; j < nj; j++) {
                data[i][j] += other(i, j);
            }
        return (*this);
    }

    /*performs element by element division by another field*/
    Field &operator*=(double s) {
        for (int i = 0; i < ni; i++)
            for (int j = 0; j < nj; j++) {
                data[i][j] *= s;
            }
        return (*this);
    }

    // multiplication operator, returns f*s
    friend Field<T> operator*(double s, const Field<T> &f) {
        Field<T> r(f);
        return r *= s;
    }

    /*sets all data to zero*/
    void clear() { (*this) = 0; }

    /* scatters scalar value onto a field at logical coordinate lc*/
    void scatter(Vec2d lc, T value);

    /* gathers field value at logical coordinate lc*/
    T gather(Vec2d lc);

    template <typename S>
    friend std::ostream &operator<<(std::ostream &out, Field<S> &f);
};

template <typename T>
void Field<T>::scatter(Vec2d lc, T value) {
    int i = (int)lc[0];
    double di = lc[0] - i;

    int j = (int)lc[1];
    double dj = lc[1] - j;

    data[i][j] += (T)value * (1 - di) * (1 - dj);
    data[i + 1][j] += (T)value * (di) * (1 - dj);
    data[i][j + 1] += (T)value * (1 - di) * (dj);
    data[i + 1][j + 1] += (T)value * (di) * (dj);
}

template <typename T>
T Field<T>::gather(Vec2d lc) {
    int i = (int)lc[0];
    double di = lc[0] - i;

    int j = (int)lc[1];
    double dj = lc[1] - j;

    T val = data[i][j] * (1 - di) * (1 - dj) +
            data[i + 1][j] * (di) * (1 - dj) +
            data[i][j + 1] * (1 - di) * (dj) + data[i + 1][j + 1] * (di) * (dj);

    return val;
}

/*writes out data to a file stream*/
template <typename T>
std::ostream &operator<<(std::ostream &out, Field<T> &f) {
    // for (int k=0;k<f.nk;k++,out<<"\n")
    for (int j = 0; j < f.nj; j++, out << "\n") {
        for (int i = 0; i < f.ni; i++) out << f.data[i][j] << " ";
    }
    return out;
}
