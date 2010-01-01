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
/* Copyright (C) 1997-2001 Ken Turkowski. <turk_at_computer.org>
 *
 * All rights reserved.
 *
 * Warranty Information
 *  Even though I have reviewed this software, I make no warranty
 *  or representation, either express or implied, with respect to this
 *  software, its quality, accuracy, merchantability, or fitness for a
 *  particular purpose.  As a result, this software is provided "as is,"
 *  and you, its user, are assuming the entire risk as to its quality
 *  and accuracy.
 *
 * This code may be used and freely distributed as long as it includes
 * this copyright notice and the above warranty information.
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
