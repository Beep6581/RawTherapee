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
#pragma once

#include <array>

namespace rtengine
{

namespace homogeneous
{

enum Axis {X, Y, Z};

template <typename T>
using Matrix = std::array<std::array<T, 4>, 4>;

/**
 * 3 dimensional homogeneous vector.
 */
template <typename T>
using Vector = std::array<T, 4>;

/**
 * Creates a 3 dimensional transformation matrix for projection onto a plane.
 * @param location Distance from the origin to the plane.
 * @param normal Direction of the plane's normal.
 */
template <typename T>
Matrix<T> projectionMatrix(T location, Axis normal);

/**
 * Creates a 3 dimensional transformation matrix for rotation.
 * @param radians Rotation angle.
 * @param axis Axis of rotation.
 */
template <typename T>
Matrix<T> rotationMatrix(double radians, Axis axis);

/**
 * Creates a 3 dimensional transformation matrix for scaling.
 * @param x Scale in x-direction
 * @param y Scale in y-direction
 * @param z Scale in z-direction
 */
template <typename T>
Matrix<T> scaleMatrix(T x, T y, T z);

/**
 * Creates a 3 dimensional transformation matrix for translation.
 * @param x Translation in the the x-direction.
 * @param y Translation in the the y-direction.
 * @param z Translation in the the z-direction.
 */
template <typename T>
Matrix<T> translationMatrix(T x, T y, T z);

}

template <typename T>
homogeneous::Vector<T> operator*(const homogeneous::Matrix<T>& a, const homogeneous::Vector<T>& b);

template <typename T>
homogeneous::Matrix<T> operator*(const homogeneous::Matrix<T>& a, const homogeneous::Matrix<T>& b);

}

#include "homogeneouscoordinates.cc"
