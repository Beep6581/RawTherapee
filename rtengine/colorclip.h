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
inline double tightestroot (double L, double a, double b, double r1, double r2, double r3);

#ifndef __COLORCLIP__
#define __COLORCLIP__

#include <cmath>
#include "median.h"

// gives back the tightest >0 amplification by which color clipping occures

inline double tightestroot (double L, double a, double b, double r1, double r2, double r3)
{

    double an = a / 500.0, bn = b / 200.0, p = (L + 16.0) / 116.0;

    double coeff3 = r1 * an * an * an - r3 * bn * bn * bn;
    double coeff2 = 3.0 * p * (r1 * an * an + r3 * bn * bn);
    double coeff1 = 3.0 * p * p * (r1 * an - r3 * bn);
    double coeff0 = p * p * p * (r1 + r2 + r3) - 1.0;

    double a1 = coeff2 / coeff3;
    double a2 = coeff1 / coeff3;
    double a3 = coeff0 / coeff3;

    double Q = (a1 * a1 - 3.0 * a2) / 9.0;
    double R = (2.0 * a1 * a1 * a1 - 9.0 * a1 * a2 + 27.0 * a3) / 54.0;
    double Qcubed = Q * Q * Q;
    double d = Qcubed - R * R;

//  printf ("input L=%g, a=%g, b=%g\n", L, a, b);
//  printf ("c1=%g, c2=%g, c3=%g, c4=%g\n", coeff3, coeff2, coeff1, coeff0);


    /* Three real roots */
    if (d >= 0) {
        double theta = acos(R / sqrt(Qcubed));
        double sqrtQ = sqrt(Q);
        double x0 = -2.0 * sqrtQ * cos( theta               / 3.0) - a1 / 3.0;
        double x1 = -2.0 * sqrtQ * cos((theta + 2.0 * M_PI) / 3.0) - a1 / 3.0;
        double x2 = -2.0 * sqrtQ * cos((theta + 4.0 * M_PI) / 3.0) - a1 / 3.0;

//    printf ("3 roots: %g, %g, %g\n", x0, x1, x2);

        SORT3 (x0, x1, x2, a1, a2, a3);

        if (a1 > 0) {
            return a1;
        }

        if (a2 > 0) {
            return a2;
        }

        if (a3 > 0) {
            return a3;
        }

        return -1;
    }

    /* One real root */
    else {
//    double e = pow(sqrt(-d) + fabs(R), 1.0 / 3.0);
        double e = exp (1.0 / 3.0 * log (sqrt(-d) + fabs(R)));

        if (R > 0) {
            e = -e;
        }

        double x0 = (e + Q / e) - a1 / 3.0;

//    printf ("1 root: %g\n", x0);

        if (x0 < 0) {
            return -1;
        } else {
            return x0;
        }
    }
}


/*******************************************************************************
 * FindCubicRoots
 * --------------
 *
 * Copyright (C) 1997-2001 Ken Turkowski. <turk_at_computer.org>
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * ------------------------------------------------------------------------
 *
 *  Solve:
 *      coeff[3] * x^3 + coeff[2] * x^2 + coeff[1] * x + coeff[0] = 0
 *
 *  returns:
 *      3 - 3 real roots
 *      1 - 1 real root (2 complex conjugate)
 *
 *******************************************************************************/

/*long
FindCubicRoots(const FLOAT coeff[4], FLOAT x[3])
{
    FLOAT a1 = coeff[2] / coeff[3];
    FLOAT a2 = coeff[1] / coeff[3];
    FLOAT a3 = coeff[0] / coeff[3];

    double_t Q = (a1 * a1 - 3 * a2) / 9;
    double_t R = (2 * a1 * a1 * a1 - 9 * a1 * a2 + 27 * a3) / 54;
    double_t Qcubed = Q * Q * Q;
    double_t d = Qcubed - R * R;

    // Three real roots
    if (d >= 0) {
        double_t theta = acos(R / sqrt(Qcubed));
        double_t sqrtQ = sqrt(Q);
        x[0] = -2 * sqrtQ * cos( theta           / 3) - a1 / 3;
        x[1] = -2 * sqrtQ * cos((theta + 2 * pi) / 3) - a1 / 3;
        x[2] = -2 * sqrtQ * cos((theta + 4 * pi) / 3) - a1 / 3;
        return (3);
    }

    // One real root
    else {
        double_t e = pow(sqrt(-d) + fabs(R), 1. / 3.);
        if (R > 0)
            e = -e;
        x[0] = (e + Q / e) - a1 / 3.;
        return (1);
    }
}
*/
#endif
