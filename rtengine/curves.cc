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
#include <curves.h>
#include <math.h>
#include <vector>
#include <mytime.h>
#include <string.h>

#undef CLIPD
#define CLIPD(a) ((a)>0.0?((a)<1.0?(a):1.0):0.0)


namespace rtengine {

Curve::Curve (const std::vector<double>& p) : x(NULL), y(NULL), ypp(NULL) {

    if (p.size()<3) {
        kind = 0;
    }
    else {
        kind = p[0];
        if (kind==-1 || kind==1) {
            N = (p.size()-1)/2;
            x = new double[N];
            y = new double[N];
            int ix = 1;
            for (int i=0; i<N; i++) {
                x[i] = p[ix++];
                y[i] = p[ix++];
            }
            if (kind==1)
                spline_cubic_set ();
        }
        if (kind==2) {
            if (p.size()!=8 && p.size()!=9)
                kind = 0;
            else {
                x = new double[9];
                for (int i=0; i<4; i++)
                    x[i] = p[i];
                for (int i=4; i<8; i++)
                    x[i] = (p[i]+100.0)/200.0;
                if (p.size()<9)
                    x[8] = 1.0;
                else
                    x[8] = p[8]/100.0;
            }
        }
    }
}  

Curve::~Curve () {

    delete [] x;
    delete [] y;
    delete [] ypp;
}

void Curve::spline_cubic_set () {

    double* u = new double[N-1];
    delete [] ypp;
    ypp = new double [N];

    ypp[0] = u[0] = 0.0;	/* set lower boundary condition to "natural" */

    for (int i = 1; i < N - 1; ++i) {
        double sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        double p = sig * ypp[i - 1] + 2.0;
        ypp[i] = (sig - 1.0) / p;
        u[i] = ((y[i + 1] - y[i])
          / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]));
        u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }

    ypp[N - 1] = 0.0;
    for (int k = N - 2; k >= 0; --k)
        ypp[k] = ypp[k] * ypp[k + 1] + u[k];

    delete [] u;
}

double Curve::getVal (double t) {

    if (!kind)
        return t;

    if (kind==2) {
        
        if (t<=1e-14)
            return 0.0;
        double c = -log(2.0)/log(x[2]);
        double tv = exp(c*log(t));
        double base = pfull (tv, x[8], x[6], x[5]);
        double stretched = base<=1e-14 ? 0.0 : exp(log(base)/c);
        
        base = pfull (0.5, x[8], x[6], x[5]);   
        double fc = base<=1e-14 ? 0.0 : exp(log(base)/c);   // value of the curve at the center point
        if (t<x[2]) {
            // add shadows effect:
            double sc = -log(2.0)/log(x[1]/x[2]);
            double stv = exp(sc*log(stretched/fc));
            double sbase = pfull (stv, x[8], x[7], 0.5);
            double sstretched = fc*(sbase<=1e-14 ? 0.0 : exp(log(sbase)/sc));
            return sstretched;
        }
        else {
            // add highlights effect:
            double hc = -log(2.0)/log((x[3]-x[2])/(1-x[2]));
            double htv = exp(hc*log((stretched-fc)/(1-fc)));
            double hbase = pfull (htv, x[8], 0.5, x[4]);
            double hstretched = fc + (1-fc)*(hbase<=1e-14 ? 0.0 : exp(log(hbase)/hc));
            return hstretched;
        }
    }
    else {
        if (t>x[N-1])
            return y[N-1];
        else if (t<x[0])
            return y[0];

        /* do a binary search for the right interval: */
        int k_lo = 0, k_hi = N - 1;
        while (k_hi - k_lo > 1){
            int k = (k_hi + k_lo) / 2;
            if (x[k] > t)
                k_hi = k;
            else
                k_lo = k;
        }

        double h = x[k_hi] - x[k_lo];
        if (kind==-1)
            return y[k_lo] + (t - x[k_lo]) * ( y[k_hi] - y[k_lo] ) / h;
        else if (kind==1) {
            double a = (x[k_hi] - t) / h;
            double b = (t - x[k_lo]) / h;
            double r = a*y[k_lo] + b*y[k_hi] + ((a*a*a - a)*ypp[k_lo] + (b*b*b - b)*ypp[k_hi]) * (h*h)/6.0;
   	        if (r < 0.0) return 0.0;
	        if (r > 1.0) return 1.0;
            return r;
        }
        else
            return t;
    }
}

void Curve::getVal (const std::vector<double>& t, std::vector<double>& res) {
    
// TODO!!!! can be made much faster!!! Binary search of getVal(double) at each point can be avoided

    res.resize (t.size());
    for (int i=0; i<t.size(); i++)
        res[i] = getVal(t[i]);
}

double CurveFactory::centercontrast (double x, double b, double m) {

  if (b==0)
    return x;
  if (b>0) {
    if (x>m) 
      return m + (1.0-m) * tanh (b*(x-m)/(1.0-m)) / tanh (b);
    else
      return m + m * tanh (b*(x-m)/m) / tanh (b);
  }
  else {
    if (x>m) 
      return 2.0*x - m - (1.0-m) * tanh (b*(x-m)/(1.0-m)) / tanh (b);
    else
      return 2.0*x - m - m * tanh (b*(x-m)/m) / tanh (b);
  }
}
/*
void CurveFactory::updateCurve3 (int* curve, int* ohistogram, const std::vector<double>& points, double defmul, double ecomp, int black, double hlcompr, double shcompr, double br, double contr, double gamma_, bool igamma, int skip) {

    double def_mul = pow (2.0, defmul);

    // compute parameters of the gamma curve
    double start = exp(gamma_*log( -0.099 / ((1.0/gamma_-1.0)*1.099 )));
    double slope = 1.099 * pow (start, 1.0/gamma_-1) - 0.099/start;
    double mul = 1.099;
    double add = 0.099;

    // theoretical maximum of the curve
    double D = gamma_>0 ? gamma (def_mul, gamma_, start, slope, mul, add) : def_mul;

    double a = pow (2.0, ecomp);
    double b = black / 65535.0;

    // curve without contrast
    double* dcurve = new double[65536];
    
    bool needcontrast = contr>0.00001 || contr<-0.00001;
    bool needigamma = !needcontrast && igamma && gamma_>0;

    // create a curve if needed
    Curve* tcurve = NULL;
    if (points.size()>0 && points[0]!=0)
        tcurve = new Curve (points);

    for (int i=0; i<=0xffff; i+= i<0xffff-skip ? skip : 1 ) {

        double val = (double)i / 65535.0;
        val *= def_mul;
        if (gamma_>0) 
            val = gamma (val, gamma_, start, slope, mul, add);
  
        val = basecurve (val, a, b, D, hlcompr/100.0, shcompr/100.0);
        val = brightness (val, br/100.0);
        
        if (tcurve) 
            val = tcurve->getVal (val);

	    if (needigamma)
            val = igamma2 (val);

        if (val>1.0)
            val = 1.0;
        else if (val<0.0)
            val = 0.0;
        dcurve[i] = val;
    }
    delete tcurve;
/*            
if (igamma) {
  FILE* f = fopen ("curve.txt","wt");
  for (int i=0; i<65536; i++)
//    fprintf (f, "%g\t%g\n", i/65535.0, basel(i/65535.0, 2, 0));
    fprintf (f, "%g\t%g\n", i/65535.0, clower(i/65535.0, 0.500015/0.5, 1.5));
//    fprintf (f, "%g\t%g\n", i/65535.0, basecurve(i/65535.0, 1.25701, 0, 1.47694, 1.0, 1.0));
//    fprintf (f, "%g\t%g\n", i/65535.0, dcurve[i]);
  fclose (f);
}
*/
/*
    int prev = 0;
    for (int i=1; i<=0xffff-skip; i++) {
        if (i%skip==0) {
            prev+=skip;
            continue;
        }
        dcurve[i] = ( dcurve[prev] * (skip - i%skip) + dcurve[prev+skip] * (i%skip) ) / skip;
    }

    if (needcontrast) {  
        // compute mean luminance of the image with the curve applied
        int sum = 0;
        double avg = 0;
        for (int i=0; i<=0xffff; i++) {
          avg += dcurve[i] * ohistogram[i];
          sum += ohistogram[i];
        }
        avg /= sum;

        // compute contrast parameter
        double contr_b = contr / 20;
        if (contr_b>=0 && contr_b < 0.00001)
          contr_b = 0.00001;
        else if (contr_b<0 && contr_b > -0.00001)
          contr_b = -0.00001;

        // apply contrast enhancement
        for (int i=0; i<=0xffff; i++) {
          double val = centercontrast (dcurve[i], contr_b, avg);
          if (igamma && gamma_>0)
            val = igamma2 (val);
          if (val>1.0) val = 1.0;
          if (val<0.0) val = 0.0;
          curve[i] = (int) (65535.0 * val);
        }
    }
    else 
        for (int i=0; i<=0xffff; i++) 
            curve[i] = (int) (65535.0 * dcurve[i]);
    delete [] dcurve;
}*/

void CurveFactory::complexCurve (double ecomp, double black, double hlcompr, double shcompr, double br, double contr, double defmul, double gamma_, bool igamma, const std::vector<double>& curvePoints, unsigned int* histogram, unsigned int* outCurve, unsigned int* outBeforeCCurveHistogram, int skip) {

    double def_mul = pow (2.0, defmul);

    // compute parameters of the gamma curve
    double start = exp(gamma_*log( -0.099 / ((1.0/gamma_-1.0)*1.099 )));
    double slope = 1.099 * pow (start, 1.0/gamma_-1) - 0.099/start;
    double mul = 1.099;
    double add = 0.099;

    // theoretical maximum of the curve
    double D = gamma_>0 ? gamma (def_mul, gamma_, start, slope, mul, add) : def_mul;

    // a: slope of the curve, black: starting point at the x axis
    double a = pow (2.0, ecomp);

    // curve without contrast
    double* dcurve = new double[65536];
    
    // check if contrast curve is needed
    bool needcontrast = contr>0.00001 || contr<-0.00001;
    
    // check if inverse gamma is needed at the end
    bool needigamma = !needcontrast && igamma && gamma_>0;

    // create a curve if needed
    Curve* tcurve = NULL;
    if (curvePoints.size()>0 && curvePoints[0]!=0)
        tcurve = new Curve (curvePoints);

    // clear array that stores histogram valid before applying the custom curve
    if (outBeforeCCurveHistogram)
        memset (outBeforeCCurveHistogram, 0, 256*sizeof(int));

    for (int i=0; i<=0xffff; i+= i<0xffff-skip ? skip : 1 ) {

        // change to [0,1] rage
        double val = (double)i / 65535.0;

        // apply default multiplier (that is >1 if highlight recovery is on)
        val *= def_mul;
        
        // gamma correction
        if (gamma_>0) 
            val = gamma (val, gamma_, start, slope, mul, add);
  
        // apply base curve, thus, exposure compensation and black point with shadow and highlight protection
        val = basecurve (val, a, black, D, hlcompr/100.0, shcompr/100.0);

        // apply brightness curve
        val = brightness (val, br/100.0);
        
        // apply custom/parametric curve, if any
        if (tcurve) {
            if (outBeforeCCurveHistogram) {
                double hval = val;
//                if (needigamma)
//                    hval = igamma2 (hval);
                int hi = (int)(255.0*CLIPD(hval));
                outBeforeCCurveHistogram[hi]+=histogram[i] ;
            }
            val = tcurve->getVal (val);
        }

	    // if inverse gamma is needed, do it (standard sRGB inverse gamma is applied)
        if (needigamma)
            val = igamma2 (val);

        // store result in a temporary array
        dcurve[i] = CLIPD(val);
    }
    delete tcurve;
    
    // if skip>1, let apply linear interpolation in the skipped points of the curve
    int prev = 0;
    for (int i=1; i<=0xffff-skip; i++) {
        if (i%skip==0) {
            prev+=skip;
            continue;
        }
        dcurve[i] = ( dcurve[prev] * (skip - i%skip) + dcurve[prev+skip] * (i%skip) ) / skip;
    }

    if (needcontrast) {  
        // compute mean luminance of the image with the curve applied
        int sum = 0;
        double avg = 0;
        for (int i=0; i<=0xffff; i++) {
          avg += dcurve[i] * histogram[i];
          sum += histogram[i];
        }
        avg /= sum;

        // compute contrast parameter
        double contr_b = contr / 20;
        if (contr_b>=0 && contr_b < 0.00001)
          contr_b = 0.00001;
        else if (contr_b<0 && contr_b > -0.00001)
          contr_b = -0.00001;

        // apply contrast enhancement
        for (int i=0; i<=0xffff; i++) {
          double val = centercontrast (dcurve[i], contr_b, avg);
          if (igamma && gamma_>0)
            val = igamma2 (val);
          outCurve[i] = (int) (65535.0 * CLIPD(val));
        }
    }
    else 
        for (int i=0; i<=0xffff; i++) 
            outCurve[i] = (int) (65535.0 * dcurve[i]);
    delete [] dcurve;
}


int CurveFactory::gammatab [65536];
int CurveFactory::igammatab_srgb [65536];
int CurveFactory::gammatab_srgb [65536];

void CurveFactory::init () {

  for (int i=0; i<65536; i++)
    gammatab_srgb[i] = (int)(65535 * gamma2 (i/65535.0));
  for (int i=0; i<65536; i++)
    igammatab_srgb[i] = (int)(65535 * igamma2 (i/65535.0));
  for (int i=0; i<65536; i++)
    gammatab[i] = (int)(65535 * pow (i/65535.0, 0.454545));
    
/*    FILE* f = fopen ("c.txt", "wt");
    for (int i=0; i<256; i++)
        fprintf (f, "%g %g\n", i/255.0, clower (i/255.0, 2.0, 1.0));
    fclose (f);*/
}

}

