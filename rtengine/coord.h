/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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

#ifndef __COORD__
#define __COORD__

namespace rtengine
{

struct Coord;
struct PolarCoord;

// Do not confuse with Coord2D, this one is used by the UI.
struct Coord
{
    int x = 0;
    int y = 0;

    Coord () = default;
    Coord (const int x, const int y);
    Coord (const Coord& other) = default;
    Coord (const PolarCoord& other);

    Coord& operator= (const Coord& other) = default;
    Coord& operator= (const PolarCoord& other);

    void get (int& x, int& y) const;
    void set (const int x, const int y);

    bool clip (const int width, const int height);

    Coord& operator+= (const Coord& other);
    Coord& operator-= (const Coord& other);
    Coord& operator*= (const double scale);

};

bool operator== (const Coord& lhs, const Coord& rhs);
bool operator!= (const Coord& lhs, const Coord& rhs);

const Coord operator+ (const Coord& lhs, const Coord& rhs);
const Coord operator- (const Coord& lhs, const Coord& rhs);
const Coord operator* (const Coord& lhs, const Coord& rhs);

struct PolarCoord
{
    double radius = 0.0;
    double angle = 0.0;

    PolarCoord () = default;
    PolarCoord (const double radius, const double angle);
    PolarCoord (const PolarCoord& other) = default;
    PolarCoord (const Coord& other);

    PolarCoord& operator= (const PolarCoord& other) = default;
    PolarCoord& operator= (const Coord& other);

    void get (double& radius, double& angle) const;
    void set (const double radius, const double angle);

    PolarCoord& operator+= (const PolarCoord& other);
    PolarCoord& operator-= (const PolarCoord& other);
    PolarCoord& operator*= (const double scale);

};

bool operator== (const PolarCoord& lhs, const PolarCoord& rhs);
bool operator!= (const PolarCoord& lhs, const PolarCoord& rhs);

const PolarCoord operator+ (const PolarCoord& lhs, const PolarCoord& rhs);
const PolarCoord operator- (const PolarCoord& lhs, const PolarCoord& rhs);
const PolarCoord operator* (const PolarCoord& lhs, const double rhs);
const PolarCoord operator* (const double lhs, const PolarCoord& rhs);

inline Coord::Coord (const int x, const int y) : x (x), y (y)
{
}

inline Coord::Coord (const PolarCoord& other)
{
    *this = other;
}

inline void Coord::get (int& x, int& y) const
{
    x = this->x;
    y = this->y;
}

inline void Coord::set (const int x, const int y)
{
    this->x = x;
    this->y = y;
}

inline Coord& Coord::operator+= (const Coord& other)
{
    x += other.x;
    y += other.y;
    return *this;
}

inline Coord& Coord::operator-= (const Coord& other)
{
    x -= other.x;
    y -= other.y;
    return *this;
}

inline Coord& Coord::operator*= (const double scale)
{
    x *= scale;
    y *= scale;
    return *this;
}

inline bool operator== (const Coord& lhs, const Coord& rhs)
{
    return lhs.x == rhs.x && lhs.y == rhs.y;
}

inline bool operator!= (const Coord& lhs, const Coord& rhs)
{
    return !(lhs == rhs);
}

inline const Coord operator+ (const Coord& lhs, const Coord& rhs)
{
    return Coord (lhs) += rhs;
}

inline const Coord operator- (const Coord& lhs, const Coord& rhs)
{
    return Coord (lhs) -= rhs;
}

inline const Coord operator* (const Coord& lhs, const double rhs)
{
    return Coord (lhs) *= rhs;
}

inline const Coord operator* (const double lhs, const Coord& rhs)
{
    return Coord (rhs) *= lhs;
}

inline PolarCoord::PolarCoord (const double radius, const double angle) : radius (radius), angle (angle)
{
}

inline PolarCoord::PolarCoord (const Coord& other)
{
    *this = other;
}

inline void PolarCoord::get (double& radius, double& angle) const
{
    radius = this->radius;
    angle = this->angle;
}

inline void PolarCoord::set (const double radius, const double angle)
{
    this->radius = radius;
    this->angle = angle;
}

inline PolarCoord& PolarCoord::operator+= (const PolarCoord& other)
{
    *this = Coord (*this) + Coord (other);
    return *this;
}

inline PolarCoord &PolarCoord::operator-= (const PolarCoord &other)
{
    *this = Coord (*this) - Coord (other);
    return *this;
}

inline PolarCoord &PolarCoord::operator*= (const double scale)
{
    radius *= scale;
    return *this;
}

inline bool operator== (const PolarCoord& lhs, const PolarCoord& rhs)
{
    return lhs.radius == rhs.radius && lhs.angle == rhs.angle;
}

inline bool operator!= (const PolarCoord& lhs, const PolarCoord& rhs)
{
    return !(lhs == rhs);
}

inline const PolarCoord operator+ (const PolarCoord& lhs, const PolarCoord& rhs)
{
    return PolarCoord (lhs) += rhs;
}

inline const PolarCoord operator- (const PolarCoord& lhs, const PolarCoord& rhs)
{
    return PolarCoord (lhs) -= rhs;
}

inline const PolarCoord operator* (const PolarCoord& lhs, const double rhs)
{
    return PolarCoord (lhs) *= rhs;
}

inline const PolarCoord operator* (const double lhs, const PolarCoord& rhs)
{
    return PolarCoord (rhs) *= lhs;
}

}

#endif
