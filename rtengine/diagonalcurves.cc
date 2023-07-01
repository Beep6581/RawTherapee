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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
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
        hashSize = 100;    // Arbitrary cut-off value, but multiple of 10
    }

    if (ppn < 50) {
        hashSize = 10;    // Arbitrary cut-off value, but multiple of 10
    }

    if (p.size() < 3) {
        kind = DCT_Empty;
    } else {
        bool identity = true;
        kind = (DiagonalCurveType)p[0];

        if (kind == DCT_Linear || kind == DCT_Spline || kind == DCT_NURBS || kind == DCT_CatumullRom) {
            N = (p.size() - 1) / 2;
            x = new double[N];
            y = new double[N];
            int ix = 1;

            for (int i = 0; i < N; i++) {
                x[i] = p[ix++];
                y[i] = p[ix++];

                if (std::fabs(x[i] - y[i]) >= 0.000009) {
                    // the smallest possible difference between x and y curve point values is ~ 0.00001
                    // checking against >= 0.000009 is a bit saver than checking against >= 0.00001
                    identity = false;
                }
            }

            if (x[0] != 0.0 || x[N - 1] != 1.0)
                // Special (and very rare) case where all points are on the identity line but
                // not reaching the limits
            {
                identity = false;
            }

            if(x[0] == 0.0 && x[1] == 0.0)
                // Avoid crash when first two points are at x = 0 (git Issue 2888)
            {
                x[1] = 0.01;
            }

            if(x[0] == 1.0 && x[1] == 1.0)
                // Avoid crash when first two points are at x = 1 (100 in gui) (git Issue 2923)
            {
                x[0] = 0.99;
            }

            if (!identity) {
                if (kind == DCT_Spline && N > 2) {
                    spline_cubic_set ();
                } else if (kind == DCT_NURBS && N > 2) {
                    NURBS_set ();
                    fillHash();
                } else if (kind == DCT_CatumullRom && N > 2) {
                    catmull_rom_set();
                } else {
                    kind = DCT_Linear;
                }
            }
        } else if (kind == DCT_Parametric) {
            if ((p.size() == 8 || p.size() == 9) && (p.at(4) != 0.0 || p.at(5) != 0.0 || p.at(6) != 0.0 || p.at(7) != 0.0)) {
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
                double mbase = pfull_alt (0.5, x[6], x[5]);
                mfc = xexp(xlog(mbase) / mc);        // value of the curve at the center point
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


/*****************************************************************************
 * Catmull Rom Spline
 * (https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline)
 *****************************************************************************/

namespace {

inline double pow2(double x)
{
    return x*x;
}


inline double catmull_rom_tj(double ti,
                             double xi, double yi,
                             double xj, double yj)
{
    // see https://github.com/Beep6581/RawTherapee/pull/4701#issuecomment-414054187
    static constexpr double alpha = 0.375;
    return pow(sqrt(pow2(xj-xi) + pow2(yj-yi)), alpha) + ti;
}


inline void catmull_rom_spline(int n_points,
                               double p0_x, double p0_y,
                               double p1_x, double p1_y,
                               double p2_x, double p2_y,
                               double p3_x, double p3_y,
                               std::vector<double> &res_x,
                               std::vector<double> &res_y)
{
    res_x.reserve(n_points);
    res_y.reserve(n_points);

    double t0 = 0;
    double t1 = catmull_rom_tj(t0, p0_x, p0_y, p1_x, p1_y);
    double t2 = catmull_rom_tj(t1, p1_x, p1_y, p2_x, p2_y);
    double t3 = catmull_rom_tj(t2, p2_x, p2_y, p3_x, p3_y);

    double space = (t2-t1) / n_points;

    res_x.push_back(p1_x);
    res_y.push_back(p1_y);

    // special case, a segment at 0 or 1 is computed exactly
    if (p1_y == p2_y && (p1_y == 0 || p1_y == 1)) {
        for (int i = 1; i < n_points-1; ++i) {
            double t = p1_x + space * i;
            if (t >= p2_x) {
                break;
            }
            res_x.push_back(t);
            res_y.push_back(p1_y);
        }
    } else {
        for (int i = 1; i < n_points-1; ++i) {
            double t = t1 + space * i;
        
            double c = (t1 - t)/(t1 - t0);
            double d = (t - t0)/(t1 - t0);
            double A1_x = c * p0_x + d * p1_x;
            double A1_y = c * p0_y + d * p1_y;

            c = (t2 - t)/(t2 - t1);
            d = (t - t1)/(t2 - t1);
            double A2_x = c * p1_x + d * p2_x;
            double A2_y = c * p1_y + d * p2_y;

            c = (t3 - t)/(t3 - t2);
            d = (t - t2)/(t3 - t2);
            double A3_x = c * p2_x + d * p3_x;
            double A3_y = c * p2_y + d * p3_y;

            c = (t2 - t)/(t2 - t0);
            d = (t - t0)/(t2 - t0);
            double B1_x = c * A1_x + d * A2_x;
            double B1_y = c * A1_y + d * A2_y;

            c = (t3 - t)/(t3 - t1);
            d = (t - t1)/(t3 - t1);
            double B2_x = c * A2_x + d * A3_x;
            double B2_y = c * A2_y + d * A3_y;        

            c = (t2 - t)/(t2 - t1);
            d = (t - t1)/(t2 - t1);
            double C_x = c * B1_x + d * B2_x;
            double C_y = c * B1_y + d * B2_y;

            res_x.push_back(C_x);
            res_y.push_back(C_y);
        }
    }

    res_x.push_back(p2_x);
    res_y.push_back(p2_y);
}


inline void catmull_rom_reflect(double px, double py, double cx, double cy,
                                double &rx, double &ry)
{
#if 0
    double dx = px - cx;
    double dy = py - cy;
    rx = cx - dx;
    ry = cy - dy;
#else
    // see https://github.com/Beep6581/RawTherapee/pull/4701#issuecomment-414054187
    static constexpr double epsilon = 1e-5;
    double dx = px - cx;
    double dy = py - cy;
    rx = cx - dx * 0.01;
    ry = dx > epsilon ? (dy / dx) * (rx - cx) + cy : cy;
#endif    
}


void catmull_rom_chain(int n_points, int n_cp, double *x, double *y,
                       std::vector<double> &res_x, std::vector<double> &res_y)
{
    double x_first, y_first;
    double x_last, y_last;
    catmull_rom_reflect(x[1], y[1], x[0], y[0], x_first, y_first);
    catmull_rom_reflect(x[n_cp-2], y[n_cp-2], x[n_cp-1], y[n_cp-1], x_last, y_last);

    int segments = n_cp - 1;

    res_x.reserve(n_points);
    res_y.reserve(n_points);

    for (int i = 0; i < segments; ++i) {
        int n = max(int(n_points * (x[i+1] - x[i]) + 0.5), 2);
        catmull_rom_spline(
            n, i == 0 ? x_first : x[i-1], i == 0 ? y_first : y[i-1],
            x[i], y[i], x[i+1], y[i+1],
            i == segments-1 ? x_last : x[i+2],
            i == segments-1 ? y_last : y[i+2],
            res_x, res_y);
    }
}

} // namespace


void DiagonalCurve::catmull_rom_set()
{
    int n_points = max(ppn * 65, 65000);
    poly_x.clear();
    poly_y.clear();
    catmull_rom_chain(n_points, N, x, y, poly_x, poly_y);
}

/*****************************************************************************/

double DiagonalCurve::getVal (double t) const
{

    switch (kind) {

    case DCT_Parametric : {
        if (t <= 1e-14) {
            return 0.0;
        }

        double tv = xexp(max(mc * xlog(t),-236.0)); // prevents numerical issues when calling pfull, at the cost of minor artifacts
        double base = pfull_alt (tv, x[6], x[5]);
        double stretched = xexp(xlog(base) / mc);

        if (t < x[2]) {
            // add shadows effect:
            double stv = xexp(max(msc * xlog(stretched / mfc),-236.0));
            double sbase = pfull_alt (stv, x[7], 0.5);
            return mfc * xexp(xlog(sbase) / msc);
        } else {
            // add highlights effect:
            double htv = xexp(max(mhc * xlog((stretched - mfc) / (1.0 - mfc)),-236.0));
            if (htv < 1e-6) {
                return stretched; // this part of the curve isn't affected by highlight, return the base curve
            } else {
                double hbase = pfull_alt (htv, 0.5, x[4]);
                return mfc + (1.0 - mfc) * xexp(xlog(hbase) / mhc);
            }
        }

        break;
    }

    case DCT_Linear :
    case DCT_Spline :
    {
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

    case DCT_CatumullRom: {
        auto it = std::lower_bound(poly_x.begin(), poly_x.end(), t);
        if (it == poly_x.end()) {
            return poly_y.back();
        }
        auto d = it - poly_x.begin();
        if (it+1 < poly_x.end() && t - *it > *(it+1) - t) {
            ++d;
        }
        return LIM01(*(poly_y.begin() + d));
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
