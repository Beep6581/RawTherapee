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
#ifndef __CURVES_H__
#define __CURVES_H__

#include <glibmm.h>
#include <map>
#include <string>
#include <math.h>

namespace rtengine {

class Curve {

  protected:
    int N;
    double* x;
    double* y;
    double* ypp;
    Glib::ustring name;
    bool islinear;
    bool isempty;

  protected:
    void d3_np_fs (double a[], double b[]);
    void spline_cubic_set ();

  public:

    Curve (const char* iname, int iN, double ix[], double iy[]);
    Curve (const char* iname, const char* descr);
    Curve (const std::vector<double>& points);
   ~Curve ();

    double getVal (double x);
    Glib::ustring getName ();
    bool   isEmpty () { return isempty; }
};

class CurveFactory {

  protected:

    static std::map<std::string, Curve*> curves;
    static int gammatab[65536];
    static int igammatab_srgb[65536];
    static int gammatab_srgb[65536];
    static double solve_upper (double m, double c, double deriv);
    static double solve_lower (double m, double c, double deriv);
    static double dupper (const double b, const double m, const double c);
    static double dlower (const double b, const double m, const double c);
    
    static inline double basel (double x, double m1, double m2) {
        if (x==0.0)
            return 0.0;
        double k = sqrt ((m1-1.0)*(m1-m2)/2) / (1.0-m2);
        double l = (m1-m2) / (1.0-m2) + k;
        double lx = log(x);
        return m2*x + (1.0-m2)*(2.0 - exp(k*lx))*exp(l*lx);
    }
    static inline double baseu (double x, double m1, double m2) {
        return 1.0 - basel(1.0-x, m1, m2);
    }
    static inline double cupper (double x, double m, double hr) {
        if (hr>1.0)
            return baseu (x, m, 2.0*(hr-1.0)/m);
        double x1 = (1.0-hr)/m;
        double x2 = x1 + hr;
        if (x>=x2) return 1.0;
        if (x<x1) return x*m;
        return 1.0 - hr + hr*baseu((x-x1)/hr, m, 0);
    }
    static inline double clower (double x, double m, double sr) {
        return 1.0 - cupper(1.0-x, m, sr);
    }
    static inline double basecurve (double x, double a, double b, double D, double hr, double sr) { // a: slope, b: black, D: max. x value
        double m = b+0.5/a<D ? b+0.5/a : D;
        double y = (D-b)*a<0.5 ? (D-b)*a : 0.5;
        if (x<=m)
            return b==0 ? x*a : clower (x/m, a*m/y, sr) * y;
        else if (b+1.0/a<D)
            return y+(1.0-y)*cupper((x-m)/(D-m), a*(D-m)/(1.0-y), hr);
        else
            return y+(x-m)*a;
    }
    static inline double brightnessbase (double x, double amount) {
        if (x<0.5)
            return x + amount*cupper(2.0*x, 4.5, 0.0)/3.0;
        else
            return x + amount*cupper(2.0-2.0*x, 1.5, 0.0)/3.0;
    }
    static inline double brightness (double x, double amount) {
        if (amount==0)
            return x;
        else if (amount>0)
            return brightnessbase (x, amount);
        else 
            return 1.0 - brightnessbase (1.0-x, -amount);
    }
    

  public:

    static inline double softClip         (double x, double d1, double d2, double a, double b, double c, double d);
    static inline double contrast         (double x, double a);
    static inline double centercontrast   (double x, double b, double m);
    static inline double brightness       (double x, double a, double bd1, double bd2);
    static inline double gamma            (double x, double gamma, double start, double slope, double mul, double add){
                                            return (x <= start ? x*slope : exp(log(x)/gamma)*mul-add);
                                          }
    static inline double gamma2            (double x) {
                                            return x <= 0.00304 ? x*12.92 : 1.055*exp(log(x)/2.4)-0.055;
                                          }
    static inline int    gamma_srgb       (int x) { return gammatab_srgb[x]; }
    static inline int    gamma            (int x) { return gammatab[x]; }
    static inline int    igamma_srgb      (int x) { return igammatab_srgb[x]; }
    static inline double igamma2           (double x) {
                                            return x <= 0.03928 ? x/12.92 : exp(log((x+0.055)/1.055)*2.4);
                                          }
    static inline double levels           (double x, double b_lower, double b_upper, double m, double cmax);


  public:
    static void    loadCurves  (Glib::ustring fname);
    static void updateCurve2 (int* curve, int* ohistogram, const std::vector<double>& cpoints, double ecomp, double br, int black, double hlcompr, double shcompr, double contr, double gamma_, bool igamma, int skip=1);
    static void updateCurve3 (int* curve, int* ohistogram, const std::vector<double>& cpoints, double defmul, double ecomp, int black, double hlcompr, double shcompr, double br, double contr, double gamma_, bool igamma, int skip=1);
    static std::vector<Glib::ustring> curveNames  ();
};
};

#endif
