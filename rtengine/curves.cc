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

namespace rtengine {

Curve::Curve (const char* iname, const char* descr) : islinear(false), isempty(false) {

  ypp = NULL;
  name = iname;
  char* buffer = new char[strlen(descr)+1];
  strcpy (buffer, descr);
  char* token = strtok (buffer, ",; \t\n");
  std::vector<double> xv;
  std::vector<double> yv;
  while (token) {
    double xd = atof (token);
    token = strtok (NULL, ",; \t\n");
    if (token) {
      double yd = atof (token);
      xv.push_back (xd);
      yv.push_back (yd);
    }
    token = strtok (NULL, ",; \t\n");
  }
  N = xv.size ();
  x = new double[N];
  y = new double[N];
  for (int i=0; i<N; i++) {
    x[i] = xv[i];
    y[i] = yv[i];
  }
  delete [] buffer;
  spline_cubic_set ();
}

Curve::Curve (const char* iname, int iN, double ix[], double iy[]) : islinear(false), isempty(false) {

  ypp = NULL;
  N = iN;
  name = iname;
  x = new double[N];
  y = new double[N];
  for (int i=0; i<N; i++) {
    x[i] = ix[i];
    y[i] = iy[i];
  }
  spline_cubic_set ();
}

Curve::Curve (const std::vector<double>& p) {

  x = NULL;
  y = NULL;
  ypp = NULL;
  name = "custom";
  isempty = true;
  N = p.size()/2;
  if (N<2)
    return;
  int ix = 0;
  islinear = p[ix++]<0;
  x = new double[N];
  y = new double[N];
  for (int i=0; i<N; i++) {
    x[i] = p[ix++];
    y[i] = p[ix++];
  }
  if (N==2 && x[0]==0.0 && y[0]==0.0 && x[1]==1.0 && y[1]==1.0)
    isempty = true;
  else {
    isempty = false;
    spline_cubic_set ();
  }
}

Curve::~Curve () {

    if (x)
      delete [] x;
    if (y)
      delete [] y;
    if (ypp)
      delete [] ypp;
}
void Curve::d3_np_fs (double a[], double b[]) {

/*  ypp = new double [N];

  for (int i=0; i<N; i++)
    ypp[i] = b[i];

  for (int i=1; i<N; i++) {
    double xmult = a[2+(i-1)*3] / a[1+(i-1)*3];
    a[1+i*3] = a[1+i*3] - xmult * a[0+i*3];
    ypp[i] = ypp[i] - xmult * ypp[i-1];
  }

  ypp[N-1] = ypp[N-1] / a[1+(N-1)*3];

  for (int i=N-2; 0<=i; i--)
    ypp[i] = (ypp[i] - a[0+(i+1)*3] * ypp[i+1]) / a[1+i*3];*/
}

void Curve::spline_cubic_set () {

/*  double *a;
  double *b;
  int i;

  a = new double [3*N];
  b = new double [N];
//
//  Set up the first equation.
//
  b[0] = 0;
  a[1+0*3] = 1.0E+00;
//  a[0+1*3] = -1.0E+00;
  a[0+1*3] = 0.0E+00;

//
//  Set up the intermediate equations.
//
  for (int i=1; i<N-1; i++)
  {
    b[i] = ( y[i+1] - y[i] ) / ( x[i+1] - x[i] )
      - ( y[i] - y[i-1] ) / ( x[i] - x[i-1] );
    a[2+(i-1)*3] = (x[i] - x[i-1]) / 6.0E+00;
    a[1+ i   *3] = (x[i+1] - x[i-1]) / 3.0E+00;
    a[0+(i+1)*3] = (x[i+1] - x[i]) / 6.0E+00;
  }
//
//  Set up the last equation.
//
  b[N-1] = 0;
  a[2+(N-2)*3] = 0.0E+00;
//  a[2+(N-2)*3] = -1.0E+00;
  a[1+(N-1)*3] = 1.0E+00;

//
//  Solve the linear system.
//
  d3_np_fs (a, b);

  delete [] a;
  delete [] b;*/
  
  double* u = new double[N-1];
  if (ypp)
    delete [] ypp;
  ypp = new double [N];

  ypp[0] = u[0] = 0.0;	/* set lower boundary condition to "natural" */

  for (int i = 1; i < N - 1; ++i)
    {
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
//
//  Determine the interval [ T(I), T(I+1) ] that contains TVAL.
//  Values below T[0] or above T[N-1] use extrapolation.
//
  if (isempty)
    return t;

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
  if (islinear)
      return y[k_lo] + (t - x[k_lo]) * ( y[k_hi] - y[k_lo] ) / h;
  else {
      double a = (x[k_hi] - t) / h;
      double b = (t - x[k_lo]) / h;
      return a*y[k_lo] + b*y[k_hi] + ((a*a*a - a)*ypp[k_lo] + (b*b*b - b)*ypp[k_hi]) * (h*h)/6.0;
  }

/*
  if (t>x[N-1])
    return y[N-1];
  else if (t<x[0])
    return y[0];
    
  int ival = N - 2;

  for (int i=0; i<N-1; i++)
    if (t < x[i+1]) {
      ival = i;
      break;
    }
//
//  In the interval I, the polynomial is in terms of a normalized
//  coordinate between 0 and 1.
//
  double dt = t - x[ival];
  double h = x[ival+1] - x[ival];

  if (islinear) {
    return y[ival] + dt * ( y[ival+1] - y[ival] ) / h;
  }
  else
    return y[ival]
    + dt * ( ( y[ival+1] - y[ival] ) / h
	   - ( ypp[ival+1] / 6.0E+00 + ypp[ival] / 3.0E+00 ) * h
    + dt * ( 0.5E+00 * ypp[ival]
    + dt * ( ( ypp[ival+1] - ypp[ival] ) / ( 6.0E+00 * h ) ) ) );
*/
}

Glib::ustring Curve::getName () {

  return name;
}

std::map<std::string, Curve*> CurveFactory::curves;

/*double CurveFactory::centercontrast (double x, double b, double m) {

  if (b==0)
    return x;
  if (b>0) {
    if (x>m) 
      return m + (1.0-m) * tanh (b*(x-m)/(1.0-m)) / tanh (b);
    else
      return m - m * tanh (b*(m-x)/m) / tanh (b);
  }
  else {
    if (x>m) 
      return 2.0*x - m - (1.0-m) * tanh (b*(x-m)/(1.0-m)) / tanh (b);
    else
      return 2.0*x - m + m * tanh (b*(m-x)/m) / tanh (b);
  }
}
*/

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



double CurveFactory::contrast (double x, double a) {

  if (a==0)
    return x;
  else if (a>0) {
    double s = (1.0+exp(-0.5*a)) / (1.0+exp(-(x-0.5)*a)) * (exp(0.5*a)-exp(-(x-0.5)*a)) / (exp(0.5*a)-exp(-0.5*a));
    return s;
  }
  else {
    double s = (1.0+exp(-0.5*a)) / (1.0+exp(-(x-0.5)*a)) * (exp(0.5*a)-exp(-(x-0.5)*a)) / (exp(0.5*a)-exp(-0.5*a));
    return 2*x - s;
  }
}

double CurveFactory::brightness (double x, double a, double bd1, double bd2) {

  if (a==1)
    return x;
  else if (a<1)
    return a*x;
  else {
    if (x < 1.0/a-bd1) 
      return a*x;
    else if (x > 1.0/a+bd2)
      return 1;
    else {
      double d = bd1+bd2;
      double s = - (-a*d*(1.0+a*(bd2-x))*(1.0+a*(bd2-x))*(-1.0+a*(bd1+x))-(2.0+a*(bd1+3.0*bd2-2.0*x))*(-1.0+a*(bd1+x))*(-1.0+a*(bd1+x)) + (-1.0+a*bd1)*(1.0+a*(bd2-x))*(1.0+a*(bd2-x))*(-2.0+a*(3.0*bd1+bd2+2.0*x))) / (a*a*a*d*d*d);
      return s;
    }
  }
}

double CurveFactory::softClip (double x, double d1, double d2, double a, double b, double c, double d) {

  if (x<1.0-d1)
    return x;
  else if (x>1.0+d2)
    return 1.0;
  else
    return a*x*x*x + b*x*x + c*x + d;
}

double CurveFactory::dlower (const double b, const double m, const double c) {

  return b / (tanh(b) * 2.0 * m);
}

double CurveFactory::dupper (const double b, const double m, const double c) {

  return b / (tanh(b) * 2.0 * (c-m));
}

double CurveFactory::solve_lower (double m, double c, double deriv) {

  double b_u = 2.0*m*deriv;
  double b_l = 0.0;

  double b;
  while (b_u-b_l > 0.0000001) {
    b = (b_u+b_l) / 2.0;
    if (dlower(b,m,c)<deriv)
      b_l = b;
    else
      b_u = b;
  }
  return b;
}

double CurveFactory::solve_upper (double m, double c, double deriv) {

  double b_u = 2.0*(c-m)*deriv;
  double b_l = 0.0;

  double b;
  while (b_u-b_l > 0.0000001) {
    b = (b_u+b_l) / 2.0;
    if (dupper(b,m,c)<deriv)
      b_l = b;
    else
      b_u = b;
  }
  return b;
}


double CurveFactory::levels (double x, double b_lower, double b_upper, double m, double cmax) {

  if (x<=m) 
    return (1.0 + tanh (b_lower*(x-m)/m) / tanh (b_lower)) / 2.0;
  else
    return (1.0 + tanh (b_upper*(x-m)/(cmax-m)) / tanh (b_upper)) / 2.0;
}

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

    for (int i=0; i<=0xffff; i+= i<0xffff-skip ? skip : 1 ) {

        double val = (double)i / 65535.0;
        val *= def_mul;
        if (gamma_>0) 
            val = gamma (val, gamma_, start, slope, mul, add);
  
        val = basecurve (val, a, b, D, hlcompr/100.0, shcompr/100.0);
        val = brightness (val, br/100.0);
        
//        if (tcurve) 
//            val = tcurve->getVal (val);

	    if (needigamma)
            val = igamma2 (val);

        if (val>1.0)
            val = 1.0;
        else if (val<0.0)
            val = 0.0;
        dcurve[i] = val;
    }
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
}

void CurveFactory::updateCurve2 (int* curve, int* ohistogram, const std::vector<double>& points, double ecomp, double br, int black, double hlcompr, double shcompr, double contr, double gamma_, bool igamma, int skip) {

  double ec_mul = pow (2, ecomp);
  double bl = black / 65535.0;
  double hi = pow (2.0,-br) + bl;

  // compute parameters of the gamma curve
  double start = exp(gamma_*log( -0.099 / ((1.0/gamma_-1.0)*1.099 )));
  double slope = 1.099 * pow (start, 1.0/gamma_-1) - 0.099/start;
  double mul = 1.099;
  double add = 0.099;

  // compute parameters of the "levels" curve
  shcompr /= 100.0;
  hlcompr /= 100.0;
  double correction = hlcompr<0.85 ? 1.0 : hlcompr+0.15;
  correction = 1.0;
  double d = pow (2.0,br);
  double m = 0.5 / d + bl;
  double c = (gamma_>0 ? gamma (ec_mul, gamma_, start, slope, mul, add) : ec_mul) * correction;
//  double c = (gamma_>0 ? gamma (ec_mul, gamma_, start, slope, mul, add) : gamma2(ec_mul)) * correction;

  double b_upper = solve_upper (m, c, d);
  double b_lower = solve_lower (m, c, d);

  // generate curve without contrast (in double)

//  Curve* tcurve = curves[type];
  Curve* tcurve = new Curve (points);
  if (tcurve->isEmpty()) {
    delete tcurve;
    tcurve = NULL;
  }

  double* dcurve = new double[65536];
    
  double bltanh = tanh (b_lower);
  double butanh = tanh (b_upper);
  
  if (d * (c - bl) < 1)
    hlcompr = 0;  

  bool needcontrast = contr>0.00001 || contr<-0.00001;
  bool needigamma = !needcontrast && igamma && gamma_>0;

  for (int i=0; i<=0xffff; i+= i<0xffff-skip ? skip : 1 ) {
    double val = (double)i / 65535.0;
    val *= ec_mul;
    if (gamma_>0) 
      val = gamma (val, gamma_, start, slope, mul, add);
  
//    double sval = levels (val, b_lower, b_upper, m, c);
//  Acceleration:
    double sval;
    if (val <= m) {
        double ttag = 2.0 / (1.0 + exp (-2.0*b_lower*(val-m)/m)) - 1.0;
//        sval = (1.0 + tanh (b_lower*(val-m)/m) / bltanh) / 2.0;
        sval = (1.0 + ttag / bltanh) / 2.0;
    }
    else {
        double ttag = 2.0 / (1.0 + exp (-2.0*b_upper*(val-m)/(c-m))) - 1.0;
//        sval = (1.0 + tanh (b_upper*(val-m)/(c-m)) / butanh) / 2.0;
        sval = (1.0 + ttag / butanh) / 2.0;
    }
        
    if (val<bl)
      val = shcompr * sval;
    else if (val>hi)
      val = (1.0 - hlcompr) + hlcompr * sval;
    else if (val<m) 
      val = (1.0 - shcompr) * d * (val - bl) + shcompr * sval;
    else
      val = (1.0 - hlcompr) * d * (val - bl) + hlcompr * sval;
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
//if (igamma) {
//  FILE* f = fopen ("curve.txt","wt");
//  for (int i=0; i<65536; i++)
//    fprintf (f, "%d\t%d\n", i, curve[i]);
//  fclose (f);
//}

    delete [] dcurve;
}

int CurveFactory::gammatab [65536];
int CurveFactory::igammatab_srgb [65536];
int CurveFactory::gammatab_srgb [65536];

void CurveFactory::loadCurves (Glib::ustring fname) {

  for (int i=0; i<65536; i++)
    gammatab_srgb[i] = (int)(65535 * gamma2 (i/65535.0));
  for (int i=0; i<65536; i++)
    igammatab_srgb[i] = (int)(65535 * igamma2 (i/65535.0));
  for (int i=0; i<65536; i++)
    gammatab[i] = (int)(65535 * pow (i/65535.0, 0.454545));

  FILE* f = g_fopen (fname.c_str(), "rt");
  if (!f)
    return;

  setlocale (LC_ALL, "C");

  char* buffer = new char[1024];
  while (buffer = fgets(buffer, 1024, f)) {
    int es = 0;
    int llen = strlen(buffer);
    for (es = 0; es<llen && buffer[es]!='='; es++);
    if (es<llen) {
      buffer[es] = 0;
      Curve* c = new Curve (strtok(buffer," \t"), buffer+es+1);
      curves[c->getName()] = c;
    }
  }
  delete buffer;

  setlocale (LC_ALL, "");
}

std::vector<Glib::ustring> CurveFactory::curveNames () {

  std::vector<Glib::ustring> ret;

  int ix = 0;
  for (std::map<std::string, Curve*>::iterator i = curves.begin(); i!=curves.end(); i++)
    ret.push_back (i->second->getName());

  return ret;
}
}

