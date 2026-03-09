/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2025 Daniel Gao <daniel.gao.work@gmail.com>
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

#include <cmath>
#include <type_traits>
#include <vector>

namespace rtengine {
namespace procparams {

template<typename T>
class Threshold final
{
public:
    Threshold(T _bottom, T _top, bool _start_at_one) :
        Threshold(_bottom, _top, 0, 0, _start_at_one, false)
    {
    }

    Threshold(T _bottom_left, T _top_left, T _bottom_right, T _top_right, bool _start_at_one) :
        Threshold(_bottom_left, _top_left, _bottom_right, _top_right, _start_at_one, true)
    {
    }

    template<typename U = T>
    typename std::enable_if<std::is_floating_point<U>::value, bool>::type operator ==(const Threshold<U>& rhs) const
    {
        if (is_double) {
            return
                std::fabs(bottom_left - rhs.bottom_left) < 1e-10
                && std::fabs(top_left - rhs.top_left) < 1e-10
                && std::fabs(bottom_right - rhs.bottom_right) < 1e-10
                && std::fabs(top_right - rhs.top_right) < 1e-10;
        } else {
            return
                std::fabs(bottom_left - rhs.bottom_left) < 1e-10
                && std::fabs(top_left - rhs.top_left) < 1e-10;
        }
    }

    template<typename U = T>
    typename std::enable_if<std::is_integral<U>::value, bool>::type operator ==(const Threshold<U>& rhs) const
    {
        if (is_double) {
            return
                bottom_left == rhs.bottom_left
                && top_left == rhs.top_left
                && bottom_right == rhs.bottom_right
                && top_right == rhs.top_right;
        } else {
            return
                bottom_left == rhs.bottom_left
                && top_left == rhs.top_left;
        }
    }

    template<typename U = T>
    typename std::enable_if<std::is_integral<U>::value, bool>::type operator !=(const Threshold<U>& rhs) const
    {
        return !(*this == rhs);
    }

    T getBottom() const
    {
        return bottom_left;
    }

    T getTop() const
    {
        return top_left;
    }

    T getBottomLeft() const
    {
        return bottom_left;
    }

    T getTopLeft() const
    {
        return top_left;
    }

    T getBottomRight() const
    {
        return bottom_right;
    }

    T getTopRight() const
    {
        return top_right;
    }

    void setValues(T bottom, T top)
    {
        bottom_left = bottom;
        top_left = top;
    }

    void setValues(T bottom_left, T top_left, T bottom_right, T top_right)
    {
        this->bottom_left = bottom_left;
        this->top_left = top_left;
        this->bottom_right = bottom_right;
        this->top_right = top_right;
    }

    bool isDouble() const
    {
        return is_double;
    }

    std::vector<T> toVector() const
    {
        if (is_double) {
            return {
                bottom_left,
                top_left,
                bottom_right,
                top_right
            };
        } else {
            return {
                bottom_left,
                top_left
            };
        }
    }

    // RT: Type of the returned value
    // RV: Type of the value on the X axis
    // RV2: Type of the maximum value on the Y axis
    template <typename RT, typename RV, typename RV2>
    RT multiply(RV x, RV2 y_max) const
    {
        const double val = x;

        if (init_eql) {
            if (is_double) {
                if (val == static_cast<double>(bottom_right) && static_cast<double>(bottom_right) == static_cast<double>(top_right)) {
                    // This handles the special case where the 2 right values are the same, then bottom one is sent back,
                    // useful if one wants to keep the bottom value even beyond the x max bound
                    return 0;
                }

                if (val >= static_cast<double>(top_right)) {
                    return y_max;
                }

                if (val > static_cast<double>(bottom_right)) {
                    return static_cast<double>(y_max * (val - static_cast<double>(bottom_right)) / (static_cast<double>(top_right) - static_cast<double>(bottom_right)));
                }
            }

            if (val >= static_cast<double>(bottom_left)) {
                return 0;
            }

            if (val > static_cast<double>(top_left)) {
                return static_cast<double>(y_max * (1. - (val - static_cast<double>(bottom_left)) / (static_cast<double>(top_left) - static_cast<double>(bottom_left))));
            }

            return y_max;
        } else {
            if (is_double) {
                if (val == static_cast<double>(bottom_right) && static_cast<double>(bottom_right) == static_cast<double>(top_right)) {
                    // This handles the special case where the 2 right values are the same, then top one is sent back,
                    // useful if one wants to keep the top value even beyond the x max bound
                    return y_max;
                }

                if (val >= static_cast<double>(bottom_right)) {
                    return 0;
                }

                if (val > static_cast<double>(top_right)) {
                    return static_cast<double>(y_max * (1.0 - (val - static_cast<double>(top_right)) / (static_cast<double>(bottom_right) - static_cast<double>(top_right))));
                }
            }

            if (val >= static_cast<double>(top_left)) {
                return y_max;
            }

            if (val > static_cast<double>(bottom_left)) {
                return static_cast<double>(y_max * (val - static_cast<double>(bottom_left)) / (static_cast<double>(top_left) - static_cast<double>(bottom_left)));
            }

            return 0;
        }
    }

private:
    Threshold(T _bottom_left, T _top_left, T _bottom_right, T _top_right, bool _start_at_one, bool _is_double) :
        bottom_left(_bottom_left),
        top_left(_top_left),
        bottom_right(_bottom_right),
        top_right(_top_right),
        init_eql(_start_at_one),
        is_double(_is_double)
    {
    }

    T bottom_left;
    T top_left;
    T bottom_right;
    T top_right;
    bool init_eql;
    bool is_double;
};

} // namespace procparams
} // namespace rtengine
