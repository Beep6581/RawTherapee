/* Copyright (C) 1997-2001 Ken Turkowski. <turk_at_computer.org>
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
 */

#include <math.h>


#define FLOAT float
#define double_t double

/*******************************************************************************
 * FindCubicRoots
 *
 *	Solve:
 *		coeff[3] * x^3 + coeff[2] * x^2 + coeff[1] * x + coeff[0] = 0
 *
 *	returns:
 *		3 - 3 real roots
 *		1 - 1 real root (2 complex conjugate)
 *******************************************************************************/

long
FindCubicRoots(const FLOAT coeff[4], FLOAT x[3])
{
	FLOAT a1 = coeff[2] / coeff[3];
	FLOAT a2 = coeff[1] / coeff[3];
	FLOAT a3 = coeff[0] / coeff[3];

	double_t Q = (a1 * a1 - 3 * a2) / 9;
	double_t R = (2 * a1 * a1 * a1 - 9 * a1 * a2 + 27 * a3) / 54;
	double_t Qcubed = Q * Q * Q;
	double_t d = Qcubed - R * R;

	/* Three real roots */
	if (d >= 0) {
		double_t theta = acos(R / sqrt(Qcubed));
		double_t sqrtQ = sqrt(Q);
		x[0] = -2 * sqrtQ * cos( theta           / 3) - a1 / 3;
		x[1] = -2 * sqrtQ * cos((theta + 2 * 3.14159265) / 3) - a1 / 3;
		x[2] = -2 * sqrtQ * cos((theta + 4 * 3.14159265) / 3) - a1 / 3;
		return (3);
	}

	/* One real root */
	else {
		double_t e = pow(sqrt(-d) + fabs(R), 1. / 3.);
		if (R > 0)
			e = -e;
		x[0] = (e + Q / e) - a1 / 3.;
		return (1);
	}
}
