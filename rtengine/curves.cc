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
#include "curves.h"
#include <math.h>
#include <string.h>
#include "mytime.h"
#include <iostream>
#include <stdio.h>

#undef CLIPD
#define CLIPD(a) ((a)>0.0?((a)<1.0?(a):1.0):0.0)


namespace rtengine {

Curve::Curve (const FloatVector& p) : x(NULL), y(NULL), ypp(NULL) {

    if (p.size()<3) {
        kind = 0;
    }
    else {
        kind = p[0];
        if (kind==-1 || kind==1) {
            N = (p.size()-1)/2;
            x = new float[N];
            y = new float[N];
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
                x = new float[9];
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

    float* u = new float[N-1];
    delete [] ypp;
    ypp = new float [N];

    ypp[0] = u[0] = 0.0;	/* set lower boundary condition to "natural" */

    for (int i = 1; i < N - 1; ++i) {
        float sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        float p = sig * ypp[i - 1] + 2.0;
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

float Curve::getVal (float t) {

    if (!kind)
        return t;

    if (kind==2) {
        
        if (t<=1e-14)
            return 0.0;
        float c = -log(2.0)/log(x[2]);
        float tv = exp(c*log(t));
        float base = pfull (tv, x[8], x[6], x[5]);
        float stretched = base<=1e-14 ? 0.0 : exp(log(base)/c);
        
        base = pfull (0.5, x[8], x[6], x[5]);   
        float fc = base<=1e-14 ? 0.0 : exp(log(base)/c);   // value of the curve at the center point
        if (t<x[2]) {
            // add shadows effect:
            float sc = -log(2.0)/log(x[1]/x[2]);
            float stv = exp(sc*log(stretched/fc));
            float sbase = pfull (stv, x[8], x[7], 0.5);
            float sstretched = fc*(sbase<=1e-14 ? 0.0 : exp(log(sbase)/sc));
            return sstretched;
        }
        else {
            // add highlights effect:
            float hc = -log(2.0)/log((x[3]-x[2])/(1-x[2]));
            float htv = exp(hc*log((stretched-fc)/(1-fc)));
            float hbase = pfull (htv, x[8], 0.5, x[4]);
            float hstretched = fc + (1-fc)*(hbase<=1e-14 ? 0.0 : exp(log(hbase)/hc));
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

        float h = x[k_hi] - x[k_lo];
        if (kind==-1)
            return y[k_lo] + (t - x[k_lo]) * ( y[k_hi] - y[k_lo] ) / h;
        else if (kind==1) {
            float a = (x[k_hi] - t) / h;
            float b = (t - x[k_lo]) / h;
            float r = a*y[k_lo] + b*y[k_hi] + ((a*a*a - a)*ypp[k_lo] + (b*b*b - b)*ypp[k_hi]) * (h*h)/6.0;
   	        if (r < 0.0) return 0.0;
	        if (r > 1.0) return 1.0;
            return r;
        }
        else
            return t;
    }
}

void Curve::getVal (const FloatVector& t, FloatVector& res) {
    
// TODO!!!! can be made much faster!!! Binary search of getVal(float) at each point can be avoided

    res.resize (t.size());
    for (int i=0; i<t.size(); i++)
        res[i] = getVal(t[i]);
}

float CurveFactory::centercontrast (float x, float b, float m) {

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

void CurveFactory::complexCurve (float ecomp, float black, float hlcompr, float shcompr, float br, float contr, float gamma_, bool igamma, const FloatVector& curvePoints, unsigned int* histogram, int curveSize, int curveScale, float* outCurve, unsigned int* outBeforeCCurveHistogram, int skip) {

    // compute parameters of the gamma curve
    float start = exp(gamma_*log( -0.099 / ((1.0/gamma_-1.0)*1.099 )));
    float slope = 1.099 * pow (start, 1.0/gamma_-1) - 0.099/start;
    float mul = 1.099;
    float add = 0.099;

    // theoretical maximum of the curve. A bit different than before, result is a bit different! (Brighter.)
    int maxVal = curveSize-1;
    while (maxVal>curveScale && histogram[maxVal]==0)
		maxVal--;
    double def_mul = (double) maxVal / curveScale;
    float D = gamma_>0 ? gamma (def_mul, gamma_, start, slope, mul, add) : def_mul;

    // a: slope of the curve, black: starting point at the x axis
    float a = pow (2.0, ecomp);

    // curve without contrast
    float* dcurve = new float[curveSize];
    memset (dcurve, 0, curveSize*sizeof(float));
    
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

    for (int i=0; i<curveSize; i+= i<curveSize-skip ? skip : 1 ) {

        // change to [0,1] rage
        float val = (float)i / curveScale;

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
                float hval = val;
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
    for (int i=1; i<=maxVal-skip; i++) {
        if (i%skip==0) {
            prev+=skip;
            continue;
        }
        dcurve[i] = ( dcurve[prev] * (skip - i%skip) + dcurve[prev+skip] * (i%skip) ) / skip;
    }

    memset (outCurve, 0, curveSize*sizeof(float));
    if (needcontrast) {  
        // compute mean luminance of the image with the curve applied
        int sum = 0;
        float avg = 0;
        for (int i=0; i<=maxVal; i++) {
          avg += dcurve[i] * histogram[i];
          sum += histogram[i];
        }
        avg /= sum;

        // compute contrast parameter
        float contr_b = contr / 20;
        if (contr_b>=0 && contr_b < 0.00001)
          contr_b = 0.00001;
        else if (contr_b<0 && contr_b > -0.00001)
          contr_b = -0.00001;

        // apply contrast enhancement
        for (int i=0; i<curveSize; i++) {
          float val = centercontrast (dcurve[i], contr_b, avg);
          if (igamma && gamma_>0)
            val = igamma2 (val);
          outCurve[i] = val;
        }
    }
    else 
        for (int i=0; i<curveSize; i++) 
            outCurve[i] = dcurve[i];
    delete [] dcurve;
/*
if (gamma_>0) {    
    FILE* f = fopen ("c2.txt", "wt");
    for (int i=0; i<maxVal; i++)
        fprintf (f, "%g %g\n", (float)i/curveScale, outCurve[i]);
    fclose (f);
}*/
}


int CurveFactory::igammatab_srgb [65536];
int CurveFactory::gammatab_srgb [65536];

void CurveFactory::init () {

  for (int i=0; i<65536; i++)
    gammatab_srgb[i] = (int)(65535 * gamma2 (i/65535.0));
  for (int i=0; i<65536; i++)
    igammatab_srgb[i] = (int)(65535 * igamma2 (i/65535.0));

/*    FILE* f = fopen ("c.txt", "wt");
    for (int i=0; i<256; i++)
        fprintf (f, "%g %g\n", i/255.0, clower (i/255.0, 2.0, 1.0));
    fclose (f);*/
}

}

