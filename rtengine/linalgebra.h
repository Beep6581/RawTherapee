/* -*- C++ -*-
 *  This file is part of ART
 *
 *  Copyright (c) 2022 Alberto Griggio
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
#include <cmath>

namespace rtengine {

template <class T>
class Vec3 {
public:
    Vec3() { data_[0] = data_[1] = data_[2] = T(); }
    Vec3(T a, T b, T c) { data_[0] = a; data_[1] = b; data_[2] = c; }

    template <class T2>
    Vec3(T2 const a[3]) { data_[0] = a[0]; data_[1] = a[1]; data_[2] = a[2]; }

    Vec3 &operator=(const Vec3 &a) = default;

    template <class T2>
    Vec3 &operator=(T2 const a[3])
    {
        data_[0] = a[0]; data_[1] = a[1]; data_[2] = a[2];
        return *this;
    }

    T operator[](int i) const { return data_[i]; }
    T &operator[](int i) { return data_[i]; }
    operator const T *() const { return data_; }
    operator T *() { return data_; }

private:
    T data_[3];
};

typedef Vec3<float> Vec3f;


template <class T>
class Mat33 {
public:
    Mat33()
    {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                data_[i][j] = T();
            }
        }
    }

    Mat33(T a00, T a01, T a02,
          T a10, T a11, T a12,
          T a20, T a21, T a22)
    {
        data_[0][0] = a00; data_[0][1] = a01; data_[0][2] = a02;
        data_[1][0] = a10; data_[1][1] = a11; data_[1][2] = a12;
        data_[2][0] = a20; data_[2][1] = a21; data_[2][2] = a22;
    }

    template <class T2>
    Mat33(const T2 m[3][3])
    {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                data_[i][j] = m[i][j];
            }
        }
    }

    Mat33 &operator=(const Mat33 &m) = default;

    template <class T2>
    Mat33 &operator=(const T2 m[3][3])
    {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                data_[i][j] = m[i][j];
            }
        }
        return *this;
    }

    T const *operator[](int i) const { return data_[i]; }
    T *operator[](int i) { return data_[i]; }
    typedef const T(*Data)[3];
    operator Data() const { return data_; }

private:
    T data_[3][3];
};


typedef Mat33<float> Mat33f;


template <class T>
Mat33<T> identity()
{
    return Mat33<T>(1, 0, 0, 0, 1, 0, 0, 0, 1);
}


template <class T>
Mat33<T> diagonal(T a, T b, T c)
{
    return Mat33<T>(a, 0, 0, 0, b, 0, 0, 0, c);
}


template <class T>
Mat33<T> transpose(T const m[3][3])
{
    return Mat33<T>(m[0][0], m[1][0], m[2][0],
                    m[0][1], m[1][1], m[2][1],
                    m[0][2], m[1][2], m[2][2]);
}

template <class T>
Mat33<T> transpose(const Mat33<T> &m)
{
    return transpose(static_cast<typename Mat33<T>::Data>(m));
}


template <class T>
bool inverse(T const m[3][3], Mat33<T> &out)
{
    const T res00 = m[1][1] * m[2][2] - m[2][1] * m[1][2];
    const T res10 = m[2][0] * m[1][2] - m[1][0] * m[2][2];
    const T res20 = m[1][0] * m[2][1] - m[2][0] * m[1][1];

    const T det = m[0][0] * res00 + m[0][1] * res10 + m[0][2] * res20;

    if (std::abs(det) >= 1.0e-10) {
        out[0][0] = res00 / det;
        out[0][1] = (m[2][1] * m[0][2] - m[0][1] * m[2][2]) / det;
        out[0][2] = (m[0][1] * m[1][2] - m[1][1] * m[0][2]) / det;
        out[1][0] = res10 / det;
        out[1][1] = (m[0][0] * m[2][2] - m[2][0] * m[0][2]) / det;
        out[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) / det;
        out[2][0] = res20 / det;
        out[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) / det;
        out[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) / det;
        return true;
    } else {
        return false;
    }
}

template <class T>
Mat33<T> inverse(const Mat33<T> &m)
{
    Mat33<T> res;
    inverse(static_cast<typename Mat33<T>::Data>(m), res);
    return res;
}

template <class T>
Mat33<T> inverse(T const m[3][3])
{
    Mat33<T> res;
    inverse(m, res);
    return res;
}

template <class T>
bool inverse(const Mat33<T> &m, Mat33<T> &out)
{
    return inverse(static_cast<typename Mat33<T>::Data>(m), out);
}


template <class T>
Mat33<T> dot_product(T const a[3][3], T const b[3][3])
{
    Mat33<T> res;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            res[i][j] = 0;

            for (int k = 0; k < 3; ++k) {
                res[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    return res;
}

template <class T>
Mat33<T> dot_product(const Mat33<T> &a, T const b[3][3])
{
    return dot_product(static_cast<typename Mat33<T>::Data>(a), b);
}

template <class T>
Mat33<T> dot_product(T const a[3][3], const Mat33<T> &b)
{
    return dot_product(a, static_cast<typename Mat33<T>::Data>(b));
}

template <class T>
Mat33<T> dot_product(const Mat33<T> &a, const Mat33<T> &b)
{
    return dot_product(static_cast<typename Mat33<T>::Data>(a), static_cast<typename Mat33<T>::Data>(b));
}


template <class T>
Vec3<T> dot_product(T const a[3][3], T const b[3])
{
    Vec3<T> res;

    for (int i = 0; i < 3; ++i) {
        res[i] = 0;
        for (int k = 0; k < 3; ++k) {
            res[i] += a[i][k] * b[k];
        }
    }

    return res;
}


template <class T>
Vec3<T> dot_product(const Mat33<T> &a, T const b[3])
{
    return dot_product(static_cast<typename Mat33<T>::Data>(a), b);
}

template <class T>
Vec3<T> dot_product(T const a[3][3], const Vec3<T> &b)
{
    return dot_product(a, static_cast<T const *>(b));
}

template <class T>
Vec3<T> dot_product(const Mat33<T> &a, const Vec3<T> &b)
{
    return dot_product(static_cast<typename Mat33<T>::Data>(a), static_cast<T const *>(b));
}


template <class T>
Mat33<T> operator*(const Mat33<T> &m, T v)
{
    return Mat33<T>(m[0][0] * v, m[0][1] * v, m[0][2] * v,
                    m[1][0] * v, m[1][1] * v, m[1][2] * v,
                    m[2][0] * v, m[2][1] * v, m[2][2] * v);
}

template <class T>
Vec3<T> operator*(const Vec3<T> &a, T v)
{
    return Vec3<T>(a[0] * v, a[1] * v, a[2] * v);
}

} // namespace rtengine
