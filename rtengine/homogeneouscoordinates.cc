/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Lawrence Lee <billee@ucdavis.edu>
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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <cmath>

#include "homogeneouscoordinates.h"

namespace rtengine
{

template <typename T>
homogeneous::Vector<T> operator*(const homogeneous::Matrix<T>& a, const homogeneous::Vector<T>& b)
{
    homogeneous::Vector<T> prod;

    prod.fill(0);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            prod[i] += a[i][j] * b[j];
        }
    }

    return prod;
}

template <typename T>
homogeneous::Matrix<T> operator*(const homogeneous::Matrix<T>& a, const homogeneous::Matrix<T>& b)
{
    homogeneous::Matrix<T> prod;

    for (int i = 0; i < 4; i++) {
        prod[i].fill(0);

        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                prod[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    return prod;
}

namespace homogeneous
{

template <typename T>
Matrix<T> projectionMatrix(T location, Axis normal)
{
    Matrix<T> matrix;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            matrix[i][j] = 0;
        }
    }

    matrix[0][0] = location;
    matrix[1][1] = location;
    matrix[2][2] = location;
    matrix[3][3] = 0;

    switch (normal) {
        case X:
            matrix[3][0] = 1;
            break;

        case Y:
            matrix[3][1] = 1;
            break;

        case Z:
            matrix[3][2] = 1;
            break;
    }

    return matrix;
}

template <typename T>
Matrix<T> rotationMatrix(double radians, Axis axis)
{
    Matrix<T> matrix;

    switch (axis) {
        case X:
            matrix[0][0] = 1;
            matrix[0][1] = 0;
            matrix[0][2] = 0;
            matrix[1][0] = 0;
            matrix[1][1] = cos(radians);
            matrix[1][2] = -sin(radians);
            matrix[2][0] = 0;
            matrix[2][1] = sin(radians);
            matrix[2][2] = cos(radians);
            break;

        case Y:
            matrix[0][0] = cos(radians);
            matrix[0][1] = 0;
            matrix[0][2] = sin(radians);
            matrix[1][0] = 0;
            matrix[1][1] = 1;
            matrix[1][2] = 0;
            matrix[2][0] = -sin(radians);
            matrix[2][1] = 0;
            matrix[2][2] = cos(radians);
            break;

        case Z:
            matrix[0][0] = cos(radians);
            matrix[0][1] = -sin(radians);
            matrix[0][2] = 0;
            matrix[1][0] = sin(radians);
            matrix[1][1] = cos(radians);
            matrix[1][2] = 0;
            matrix[2][0] = 0;
            matrix[2][1] = 0;
            matrix[2][2] = 1;
            break;
    }

    matrix[0][3] = 0;
    matrix[1][3] = 0;
    matrix[2][3] = 0;
    matrix[3][0] = 0;
    matrix[3][1] = 0;
    matrix[3][2] = 0;
    matrix[3][3] = 1;

    return matrix;
}

template <typename T>
Matrix<T> scaleMatrix(T x, T y, T z)
{
    Matrix<T> matrix;

    for (int i = 0; i < 4; i++) {
        matrix[i].fill(0);
    }

    matrix[0][0] = x;
    matrix[1][1] = y;
    matrix[2][2] = z;
    matrix[3][3] = 1;

    return matrix;
}

template <typename T>
Matrix<T> translationMatrix(T x, T y, T z)
{
    Matrix<T> matrix;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++) {
            matrix[i][j] = 0;
        }
    }

    matrix[0][0] = 1;
    matrix[1][1] = 1;
    matrix[2][2] = 1;
    matrix[0][3] = x;
    matrix[1][3] = y;
    matrix[2][3] = z;
    matrix[3][3] = 1;

    return matrix;
}

}

}
