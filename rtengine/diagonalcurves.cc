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
#include <glib.h>
#include <glib/gstdio.h>
#include "curves.h"
#include <cmath>
#include <vector>
#include "mytime.h"
#include <cstring>

#define CLIPD(a) ((a)>0.0?((a)<1.0?(a):1.0):0.0)

namespace rtengine
{

DiagonalCurve::DiagonalCurve (const std::vector<double>& p, int poly_pn)
{

    ppn = poly_pn > 65500 ? 65500 : poly_pn;

    if (ppn < 500) {
        hashSize = 100;    // Arbitrary cut-off value, but mutliple of 10
    }

    if (ppn < 50) {
        hashSize = 10;    // Arbitrary cut-off value, but mutliple of 10
    }

    if (p.size() < 3) {
        kind = DCT_Empty;
    } else {
        bool identity = true;
        kind = (DiagonalCurveType)p[0];

        if (kind == DCT_Linear || kind == DCT_Spline || kind == DCT_NURBS) {
            N = (p.size() - 1) / 2;
            x = new double[N];
            y = new double[N];
            int ix = 1;

            for (int i = 0; i < N; i++) {
                x[i] = p[ix++];
                y[i] = p[ix++];

                if (std::fabs(x[i] - y[i]) >= 0.000009) {
                    identity = false;
                }
            }

            if (x[0] != 0.0f || x[N - 1] != 1.0f)
                // Special (and very rare) case where all points are on the identity line but
                // not reaching the limits
            {
                identity = false;
            }

            if(x[0] == 0.f && x[1] == 0.f)
                // Avoid crash when first two points are at x = 0 (git Issue 2888)
            {
                x[1] = 0.01f;
            }

            if(x[0] == 1.f && x[1] == 1.f)
                // Avoid crash when first two points are at x = 1 (100 in gui) (git Issue 2923)
            {
                x[0] = 0.99f;
            }

            if (!identity) {
                if (kind == DCT_Spline && N > 2) {
                    spline_cubic_set ();
                } else if (kind == DCT_NURBS && N > 2) {
                    NURBS_set ();
                    fillHash();
                } else {
                    kind = DCT_Linear;
                }
            }
        } else if (kind == DCT_Parametric) {
            if ((p.size() == 8 || p.size() == 9) && (p.at(4) != 0.0f || p.at(5) != 0.0f || p.at(6) != 0.0f || p.at(7) != 0.0f)) {
                identity = false;

                x = new double[9];
                x[0] = p[0];

                for (int i = 1; i < 4; i++) {
                    x[i] = min(max(p[i], 0.001), 0.99);
                }

                for (int i = 4; i < 8; i++) {
                    x[i] = (p[i] + 100.0) / 200.0;
                }

                if (p.size() < 9) {
                    x[8] = 1.0;
                } else {
                    x[8] = p[8] / 100.0;
                }

                mc = -xlog(2.0) / xlog(x[2]);
                double mbase = pfull (0.5, x[8], x[6], x[5]);
                mfc = mbase <= 1e-14 ? 0.0 : xexp(xlog(mbase) / mc);        // value of the curve at the center point
                msc = -xlog(2.0) / xlog(x[1] / x[2]);
                mhc = -xlog(2.0) / xlog((x[3] - x[2]) / (1 - x[2]));
            }
        }

        if (identity) {
            kind = DCT_Empty;
        }
    }
}

DiagonalCurve::~DiagonalCurve ()
{

    delete [] x;
    delete [] y;
    delete [] ypp;
    poly_x.clear();
    poly_y.clear();
}

void DiagonalCurve::spline_cubic_set ()
{

    double* u = new double[N - 1];
    delete [] ypp;      // TODO: why do we delete ypp here since it should not be allocated yet?
    ypp = new double [N];

    ypp[0] = u[0] = 0.0;    /* set lower boundary condition to "natural" */

    for (int i = 1; i < N - 1; ++i) {
        double sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        double p = sig * ypp[i - 1] + 2.0;
        ypp[i] = (sig - 1.0) / p;
        u[i] = ((y[i + 1] - y[i])
                / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]));
        u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }

    ypp[N - 1] = 0.0;

    for (int k = N - 2; k >= 0; --k) {
        ypp[k] = ypp[k] * ypp[k + 1] + u[k];
    }

    delete [] u;
}

void DiagonalCurve::NURBS_set ()
{

    int nbSubCurvesPoints = N + (N - 3) * 2;

    std::vector<double> sc_x(nbSubCurvesPoints);  // X sub-curve points (  XP0,XP1,XP2,  XP2,XP3,XP4,  ...)
    std::vector<double> sc_y(nbSubCurvesPoints);  // Y sub-curve points (  YP0,YP1,YP2,  YP2,YP3,YP4,  ...)
    std::vector<double> sc_length(N + 2);         // Length of the subcurves
    double total_length = 0.;

    // Create the list of Bezier sub-curves
    // NURBS_set is called if N > 2 and non identity only

    int j = 0;
    int k = 0;

    for (int i = 0; i < N - 1;) {
        double length;
        double dx;
        double dy;

        // first point (on the curve)
        if (!i) {
            sc_x[j] = x[i];
            sc_y[j++] = y[i++];
        } else {
            sc_x[j] = (x[i - 1] + x[i]) / 2.;
            sc_y[j++] = (y[i - 1] + y[i]) / 2.;
        }

        // second point (control point)
        sc_x[j] = x[i];
        sc_y[j] = y[i++];

        dx = sc_x[j] - sc_x[j - 1];
        dy = sc_y[j] - sc_y[j - 1];
        length = sqrt(dx * dx + dy * dy);
        j++;

        // third point (on the curve)
        if (i == N - 1) {
            sc_x[j] = x[i];
            sc_y[j] = y[i];
        } else {
            sc_x[j] =  (x[i - 1] + x[i]) / 2.;
            sc_y[j] =  (y[i - 1] + y[i]) / 2.;
        }

        dx = sc_x[j] - sc_x[j - 1];
        dy = sc_y[j] - sc_y[j - 1];
        length += sqrt(dx * dx + dy * dy);
        j++;

        // Storing the length of all sub-curves and the total length (to have a better distribution
        // of the points along the curve)
        sc_length[k++] = length;
        total_length += length;
    }

    poly_x.clear();
    poly_y.clear();
    unsigned int sc_xsize = j - 1;
    j = 0;

    // adding the initial horizontal segment, if any
    if (x[0] > 0.) {
        poly_x.push_back(0.);
        poly_y.push_back(y[0]);
    }

    // adding the initial horizontal segment, if any
    // create the polyline with the number of points adapted to the X range of the sub-curve
    for (unsigned int i = 0; i < sc_xsize /*sc_x.size()*/; i += 3) {
        // TODO: Speeding-up the interface by caching the polyline, instead of rebuilding it at each action on sliders !!!
        nbr_points = (int)(((double)(ppn + N - 2) * sc_length[i / 3] ) / total_length);

        if (nbr_points < 0) {
            for(unsigned int it = 0; it < sc_x.size(); it += 3) { // used unsigned int instead of size_t to avoid %zu in printf
                printf("sc_length[%u/3]=%f \n", it, sc_length[it / 3]);
            }

            printf("NURBS diagonal curve: error detected!\n i=%u nbr_points=%d ppn=%d N=%d sc_length[i/3]=%f total_length=%f", i, nbr_points, ppn, N, sc_length[i / 3], total_length);
            exit(0);
        }

        // increment along the curve, not along the X axis
        increment = 1.0 / (double)(nbr_points - 1);
        x1 = sc_x[i];
        y1 = sc_y[i];
        x2 = sc_x[i + 1];
        y2 = sc_y[i + 1];
        x3 = sc_x[i + 2];
        y3 = sc_y[i + 2];
        firstPointIncluded = !i;
        AddPolygons ();
    }

    // adding the final horizontal segment, always (see under)
    poly_x.push_back(3.0);      // 3.0 is a hack for optimization purpose of the getVal method (the last value has to be beyond the normal range)
    poly_y.push_back(y[N - 1]);

    fillDyByDx();
}

double DiagonalCurve::getVal (double t) const
{

    switch (kind) {

    case DCT_Parametric : {
        if (t <= 1e-14) {
            return 0.0;
        }

        double tv = xexp(mc * xlog(t));
        double base = pfull (tv, x[8], x[6], x[5]);
        double stretched = base <= 1e-14 ? 0.0 : xexp(xlog(base) / mc);

        if (t < x[2]) {
            // add shadows effect:
            double stv = xexp(msc * xlog(stretched / mfc));
            double sbase = pfull (stv, x[8], x[7], 0.5);
            return mfc * (sbase <= 1e-14 ? 0.0 : xexp(xlog(sbase) / msc));
        } else {
            // add highlights effect:
            double htv = xexp(mhc * xlog((stretched - mfc) / (1 - mfc)));
            double hbase = pfull (htv, x[8], 0.5, x[4]);
            return mfc + (1 - mfc) * (hbase <= 1e-14 ? 0.0 : xexp(xlog(hbase) / mhc));
        }

        break;
    }

    case DCT_Linear :
    case DCT_Spline : {
        // values under and over the first and last point
        if (t > x[N - 1]) {
            return y[N - 1];
        } else if (t < x[0]) {
            return y[0];
        }

        // do a binary search for the right interval:
        unsigned int k_lo = 0, k_hi = N - 1;

        while (k_hi > 1 + k_lo) {
            unsigned int k = (k_hi + k_lo) / 2;

            if (x[k] > t) {
                k_hi = k;
            } else {
                k_lo = k;
            }
        }

        double h = x[k_hi] - x[k_lo];

        // linear
        if (kind == DCT_Linear) {
            return y[k_lo] + (t - x[k_lo]) * ( y[k_hi] - y[k_lo] ) / h;
        }
        // spline curve
        else { // if (kind==Spline) {
            double a = (x[k_hi] - t) / h;
            double b = (t - x[k_lo]) / h;
            double r = a * y[k_lo] + b * y[k_hi] + ((a * a * a - a) * ypp[k_lo] + (b * b * b - b) * ypp[k_hi]) * (h * h) * 0.1666666666666666666666666666666;
            return CLIPD(r);
        }

        break;
    }

    case DCT_NURBS : {
        // get the hash table entry by rounding the value (previously multiplied by "hashSize")
        unsigned short int i = (unsigned short int)(t * hashSize);

        if (UNLIKELY(i > (hashSize + 1))) {
            //printf("\nOVERFLOW: hash #%d is used while seeking for value %.8f, corresponding polygon's point #%d (out of %d point) x value: %.8f\n\n", i, t, hash.at(i), poly_x.size(), poly_x[hash.at(i)]);
            printf("\nOVERFLOW: hash #%d is used while seeking for value %.8f\n\n", i, t);
            return t;
        }

        unsigned int k_lo;
        unsigned int k_hi;

        k_lo = hash.at(i).smallerValue;
        k_hi = hash.at(i).higherValue;

        // do a binary search for the right interval :
        while (k_hi > 1 + k_lo) {
            unsigned int k = (k_hi + k_lo) / 2;

            if (poly_x[k] > t) {
                k_hi = k;
            } else {
                k_lo = k;
            }
        }

        return poly_y[k_lo] + (t - poly_x[k_lo]) * dyByDx[k_lo];
    }

    case DCT_Empty :
    default:
        // all other (unknown) kind
        return t;
    }
}

void DiagonalCurve::getVal (const std::vector<double>& t, std::vector<double>& res) const
{

    res.resize (t.size());

    for (unsigned int i = 0; i < t.size(); i++) {
        res[i] = getVal(t[i]);
    }
}

}
