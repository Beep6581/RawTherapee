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

Curve::Curve (const std::vector<double>& p, int poly_pn) : x(NULL), y(NULL), ypp(NULL) {

	ppn = poly_pn;

    if (p.size()<3) {
        kind = Empty;
    }
    else {
        kind = (CurveType)p[0];
        if (kind==Linear || kind==Spline || kind==NURBS) {
            N = (p.size()-1)/2;
            x = new double[N];
            y = new double[N];
            int ix = 1;
            for (int i=0; i<N; i++) {
                x[i] = p[ix++];
                y[i] = p[ix++];
            }
            if (kind==Spline)
                spline_cubic_set ();
            else if (kind==NURBS && N > 2)
						NURBS_set ();
					else kind=Linear;
        }
        else if (kind==Parametric) {
            if (p.size()!=8 && p.size()!=9)
                kind = Empty;
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
    poly_x.clear();
    poly_y.clear();
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

void Curve::NURBS_set () {

	int nbSubCurvesPoints = N + (N-3)*2;

    std::vector<double> sc_x(nbSubCurvesPoints);  // X sub-curve points (  XP0,XP1,XP2,  XP2,XP3,XP4,  ...)
    std::vector<double> sc_y(nbSubCurvesPoints);  // Y sub-curve points (  YP0,YP1,YP2,  YP2,YP3,YP4,  ...)
    std::vector<double> sc_length(N+2);           // Length of the subcurves
    double total_length=0.;

    // Create the list of Bezier sub-curves
    // NURBS_set is called if N > 2 only

    int j = 0;
    int k = 0;
    for (int i = 0; i < N-1;) {
        double length;
        double dx;
        double dy;

    	// first point (on the curve)
    	if (!i) {
    		sc_x[j] = x[i];
    		sc_y[j++] = y[i++];
    	}
    	else {
    		sc_x[j] = (x[i-1] + x[i]) / 2.;
    		sc_y[j++] = (y[i-1] + y[i]) / 2.;
    	}

		// second point (control point)
		sc_x[j] = x[i];
		sc_y[j] = y[i++];
		length = sqrt(pow(sc_x[j] - sc_x[j-1],2) + pow(sc_y[j] - sc_y[j-1],2));
		j++;

    	// third point (on the curve)
		if (i==N-1) {
			sc_x[j] = x[i];
			sc_y[j] = y[i];
		}
		else {
			sc_x[j] =  (x[i-1] + x[i]) / 2.;
			sc_y[j] =  (y[i-1] + y[i]) / 2.;
		}
		dx = sc_x[j] - sc_x[j-1];
		dy = sc_y[j] - sc_y[j-1];
		length += sqrt(dx*dx + dy*dy);
		j++;

		// Storing the length of all sub-curves and the total length (to have a better distribution
		// of the points along the curve)
	    sc_length[k++] = length;
	    total_length += length;
    }

    poly_x.clear();
   	poly_y.clear();
   	unsigned int sc_xsize=j-1;
    j = 0;
    // create the polyline with the number of points adapted to the X range of the sub-curve
    for (unsigned int i=0; i < sc_xsize /*sc_x.size()*/; i+=3) {
    	// TODO: Speeding-up the interface by caching the polyline, instead of rebuilding it at each action on sliders !!!
    	int nbr_points = (int)(((double)(ppn+N-2) * sc_length[i/3] )/ total_length);
    	if (nbr_points<0){
    		for(int it=0;it < sc_x.size(); it+=3) printf("sc_length[%d/3]=%f \n",it,sc_length[it/3]);
    		printf("NURBS: error detected!\n i=%d nbr_points=%d ppn=%d N=%d sc_length[i/3]=%f total_length=%f",i,nbr_points,ppn,N,sc_length[i/3],total_length);
    		exit(0);
    	}
    	// increment along the curve, not along the X axis
    	double increment = 1.0 / (double)(nbr_points-1);
    	if (!i) {
    		poly_x.push_back( sc_x[i]);
    		poly_y.push_back(sc_y[i]);
    	}
    	for (k=1; k<(nbr_points-1); k++) {
    		double t = k*increment;
    		double t2 = t*t;
    		double tr = 1.-t;
    		double tr2 = tr*tr;
    		double tr2t = tr*2*t;

    		// adding a point to the polyline
    		poly_x.push_back( tr2*sc_x[i] + tr2t*sc_x[i+1] + t2*sc_x[i+2]);
    		poly_y.push_back( tr2*sc_y[i] + tr2t*sc_y[i+1] + t2*sc_y[i+2]);
    	}
    	// adding the last point of the sub-curve
    	poly_x.push_back( sc_x[i+2]);
    	poly_y.push_back(sc_y[i+2]);
    }
}

double Curve::getVal (double t) {

    switch (kind) {

    case Empty :
        return t;
        break;

    case Parametric : {
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
        break;
    }
    case Linear :
    case Spline : {
    	// values under and over the first and last point
        if (t>x[N-1])
            return y[N-1];
        else if (t<x[0])
            return y[0];

        // do a binary search for the right interval:
        int k_lo = 0, k_hi = N - 1;
        while (k_hi - k_lo > 1){
            int k = (k_hi + k_lo) / 2;
            if (x[k] > t)
                k_hi = k;
            else
                k_lo = k;
        }

        double h = x[k_hi] - x[k_lo];
        // linear
        if (kind==Linear)
            return y[k_lo] + (t - x[k_lo]) * ( y[k_hi] - y[k_lo] ) / h;
        // spline curve
        else { // if (kind==Spline) {
            double a = (x[k_hi] - t) / h;
            double b = (t - x[k_lo]) / h;
            double r = a*y[k_lo] + b*y[k_hi] + ((a*a*a - a)*ypp[k_lo] + (b*b*b - b)*ypp[k_hi]) * (h*h)/6.0;
            return CLIPD(r);
        }
        break;
    }
    case NURBS : {
    	// values under and over the first and last point
        if (t>x[N-1])
            return y[N-1];
        else if (t<x[0])
            return y[0];
        else if (N == 2)
            return y[0] + (t - x[0]) * ( y[1] - y[0] ) / (x[1] - x[0]);

        // do a binary search for the right interval:
        int k_lo = 0, k_hi = poly_x.size() - 1;
        while (k_hi - k_lo > 1){
            int k = (k_hi + k_lo) / 2;
            if (poly_x[k] > t)
                k_hi = k;
            else
                k_lo = k;
        }

        double h = poly_x[k_hi] - poly_x[k_lo];
        return poly_y[k_lo] + (t - poly_x[k_lo]) * ( poly_y[k_hi] - poly_y[k_lo] ) / h;
		break;
    }
    default:
    	// all other (unknown) kind
		return t;
    }
}

void Curve::getVal (const std::vector<double>& t, std::vector<double>& res) {

// TODO!!!! can be made much faster!!! Binary search of getVal(double) at each point can be avoided

    res.resize (t.size());
    for (unsigned int i=0; i<t.size(); i++)
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

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	void CurveFactory::complexsgnCurve (double satclip, double satcompr, double saturation, const std::vector<double>& curvePoints, float* outCurve, int skip) {
				
		// check if contrast curve is needed
		bool needsaturation = (saturation<-0.0001 || saturation>0.0001);
		
		// curve without contrast
		float* dcurve = new float[65536];
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		std::vector<double> satcurvePoints;
		satcurvePoints.push_back((double)((CurveType)NURBS));
		if (saturation>0) {
			satcurvePoints.push_back(0); //black point.  Value in [0 ; 1] range
			satcurvePoints.push_back(0); //black point.  Value in [0 ; 1] range
			
			satcurvePoints.push_back(0.25+saturation/500.0); //toe point
			satcurvePoints.push_back(0.25-saturation/500.0); //value at toe point
			
			satcurvePoints.push_back(0.75-saturation/500.0); //shoulder point
			satcurvePoints.push_back(0.75+saturation/500.0); //value at shoulder point
			
			satcurvePoints.push_back(1); // white point
			satcurvePoints.push_back(1); // value at white point
		} else {
			satcurvePoints.push_back(0); 
			satcurvePoints.push_back(-0.5*(saturation/100.0)); 
			
			satcurvePoints.push_back(1); 
			satcurvePoints.push_back(1+saturation/200.0); 
		}
		Curve* satcurve = NULL;
		satcurve = new Curve (satcurvePoints, CURVES_MIN_POLY_POINTS/skip); // Actually, CURVES_MIN_POLY_POINTS = 1000,
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		// create a curve if needed
		Curve* tcurve = NULL;
		if (curvePoints.size()>0 && curvePoints[0]!=0)
			tcurve = new Curve (curvePoints, CURVES_MIN_POLY_POINTS/skip);
		
		for (int i=0; i<=0xffff; i+= i<0xffff-skip ? skip : 1 ) {
			
			// change to [0,1] range
			float val = (float)i / 65535.0;
			
			// apply saturation curve
			if (needsaturation)
				val = satcurve->getVal (val);
			
			// apply custom/parametric/NURBS curve, if any
			if (tcurve) {
				val = tcurve->getVal (val);
			}
			
			// store result in a temporary array
			dcurve[i] = (val);
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
		 
		for (int i=0; i<=0xffff; i++) 
			outCurve[i] = (65535.0 * dcurve[i]);
		delete [] dcurve;
		delete satcurve;
	}

	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	void CurveFactory::complexCurve (double ecomp, double black, double hlcompr, double shcompr, \
									 double br, double contr, double defmul, double gamma_, bool igamma, \
									 const std::vector<double>& curvePoints, unsigned int* histogram, \
									 float* hlCurve, float* shCurve, float* outCurve, \
									 unsigned int* outBeforeCCurveHistogram, int skip) {
		
		//double def_mul = pow (2.0, defmul);
		
		/*printf ("def_mul= %f ecomp= %f black= %f  hlcompr= %f shcompr= %f br= %f contr= %f defmul= %f  \
				gamma= %f, skip= %d \n",def_mul,ecomp,black,hlcompr,shcompr,br,contr,defmul,gamma_,skip);*/
		
		// compute parameters of the gamma curve
		double start = exp(gamma_*log( -0.099 / ((1.0/gamma_-1.0)*1.099 )));
		double slope = 1.099 * pow (start, 1.0/gamma_-1) - 0.099/start;
		double mul = 1.099;
		double add = 0.099;
		
		// a: slope of the curve, black: starting point at the x axis
		double a = pow (2.0, ecomp);
		
		// curve without contrast
		float* dcurve = new float[2*65536];
		
		// check if contrast curve is needed
		bool needcontrast = contr>0.00001 || contr<-0.00001;
		
		// check if inverse gamma is needed at the end
		bool needigamma = igamma && gamma_>0;
		
		// create a curve if needed
		Curve* tcurve = NULL;
		if (curvePoints.size()>0 && curvePoints[0]!=0)
			tcurve = new Curve (curvePoints, CURVES_MIN_POLY_POINTS/skip);
		
		// clear array that stores histogram valid before applying the custom curve
		if (outBeforeCCurveHistogram)
			memset (outBeforeCCurveHistogram, 0, 256*sizeof(int));
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// tone curve base. a: slope (from exp.comp.), b: black, def_mul: max. x value (can be>1), hr,sr: highlight,shadow recovery
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		std::vector<double> brightcurvePoints;
		brightcurvePoints.push_back((double)((CurveType)NURBS));
		
		brightcurvePoints.push_back(0); //black point.  Value in [0 ; 1] range
		brightcurvePoints.push_back(0); //black point.  Value in [0 ; 1] range
		
		if(br>0) {
			brightcurvePoints.push_back(0.1); //toe point
			brightcurvePoints.push_back(0.1+br/150.0); //value at toe point
			
			brightcurvePoints.push_back(0.7); //shoulder point
			brightcurvePoints.push_back(MIN(1.0,0.7+br/300.0)); //value at shoulder point
		} else {
			brightcurvePoints.push_back(0.1-br/150.0); //toe point
			brightcurvePoints.push_back(0.1); //value at toe point
			
			brightcurvePoints.push_back(MIN(1.0,0.7-br/300.0)); //shoulder point
			brightcurvePoints.push_back(0.7); //value at shoulder point
		}
		brightcurvePoints.push_back(1); // white point
		brightcurvePoints.push_back(1); // value at white point
		
		Curve* brightcurve = NULL;
		brightcurve = new Curve (brightcurvePoints, CURVES_MIN_POLY_POINTS/skip); // Actually, CURVES_MIN_POLY_POINTS = 1000,
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		float exp_scale = a;
		float scale = 65536.0;
		float comp = (ecomp)*hlcompr/100.0;
		int shoulder = round((scale/exp_scale)*(shcompr/100.0));
		//printf ("exp_scale= %f comp= %f def_mul=%f a= %f \n",exp_scale,comp,def_mul,a);
		
		for (int i=0; i<0x10000; i++) {
			
			// change to [0,1] range
			float val = (float)i / 65535.0;
			
			// apply default multiplier (that is >1 if highlight recovery is on)
			// val *= def_mul;
			
			// apply base curve, thus, exposure compensation and black point with shadow and highlight protection
			//val = basecurve (val*def_mul, a, 0, def_mul, hlcompr/100.0, 0);
			
			//hlCurve[i] = (65535.0 * CLIPD(val));
			
			if ((hlcompr>0)&&(exp_scale>1.0))
			{
				if (i>shoulder) {
					float Y = (float)(i-shoulder)*exp_scale/(scale-shoulder);
					float R = (float)(i-shoulder)*comp/(scale-shoulder);
					hlCurve[i] = log(1+Y*comp)/R;
				} else {
					hlCurve[i]=exp_scale;
				}
			} else {
				hlCurve[i]=exp_scale;
			}
			
			//%%%%%%%%%%%%%%%%%%%%%%%%%%
			// change to [0,1] range
			val = (double)i / 65535.0;
			
			val = basecurve (val, 1, black, 1, 0, 1.5*shcompr/100.0);
			
			shCurve[i] = (65535.0 * CLIPD(val));
			
			//%%%%%%%%%%%%%%%%%%%%%%%%%%
			// change to [0,1] range
			val = (double)i / 65535.0;
			
			// gamma correction
			if (gamma_>0) 
				val = gamma (val, gamma_, start, slope, mul, add);
			
			// apply brightness curve
			//val = brightness (val, br/100.0);
			val = brightcurve->getVal (val);
			
			// store result in a temporary array
			dcurve[i] = CLIPD(val);
		}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		if (needcontrast) {  
			// compute mean luminance of the image with the curve applied
			int sum = 0;
			float avg = 0; 
			//double sqavg = 0;
			for (int i=0; i<=0xffff; i++) {
				avg += dcurve[(int)shCurve[(int)hlCurve[i]*i]] * histogram[i];
				//sqavg += dcurve[i]*dcurve[i] * histogram[i];
				sum += histogram[i];
			}
			avg /= sum;
			//sqavg /= sum;
			//double stddev = sqrt(sqavg-avg*avg);
			
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			std::vector<double> contrastcurvePoints;
			contrastcurvePoints.push_back((double)((CurveType)NURBS));
			
			contrastcurvePoints.push_back(0); //black point.  Value in [0 ; 1] range
			contrastcurvePoints.push_back(0); //black point.  Value in [0 ; 1] range
			
			contrastcurvePoints.push_back(avg-avg*(0.6-contr/250.0)); //toe point
			contrastcurvePoints.push_back(avg-avg*(0.6+contr/250.0)); //value at toe point
			
			contrastcurvePoints.push_back(avg+(1-avg)*(0.6-contr/250.0)); //shoulder point
			contrastcurvePoints.push_back(avg+(1-avg)*(0.6+contr/250.0)); //value at shoulder point
			
			contrastcurvePoints.push_back(1); // white point
			contrastcurvePoints.push_back(1); // value at white point
			
			Curve* contrastcurve = NULL;
			contrastcurve = new Curve (contrastcurvePoints, CURVES_MIN_POLY_POINTS/skip); // Actually, CURVES_MIN_POLY_POINTS = 1000,
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			// apply contrast enhancement
			for (int i=0; i<=0xffff; i++) {
				//double val = centercontrast (dcurve[i], contr_b, avg);
				dcurve[i]  = contrastcurve->getVal (dcurve[i]);
			}
			delete contrastcurve;
		}
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		for (int i=0; i<=0xffff; i++) {
			float val;
			
			// apply custom/parametric/NURBS curve, if any
			if (tcurve) {
				if (outBeforeCCurveHistogram) {
					float hval = dcurve[(int)shCurve[(int)(hlCurve[i])]];
					//if (needigamma)
					//	hval = igamma2 (hval);
					int hi = (int)(255.0*CLIPD(hval));
					outBeforeCCurveHistogram[hi]+=histogram[i] ;
				}
				val = tcurve->getVal (dcurve[i]);
			} else {
				val = (dcurve[i]);
			}
			
			// if inverse gamma is needed, do it (standard sRGB inverse gamma is applied)
			if (needigamma)
				val = igamma2 (val);
			
			outCurve[i] = (65535.0 * val);
		}
		
		
		delete [] dcurve;
		delete tcurve;
		delete brightcurve; 
		if (outBeforeCCurveHistogram) {
			//for (int i=0; i<256; i++) printf("i= %d bchist= %d \n",i,outBeforeCCurveHistogram[i]);
		}
		
	}

	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void CurveFactory::complexLCurve (double br, double contr, const std::vector<double>& curvePoints, \
									 unsigned int* histogram, float* outCurve, \
									 unsigned int* outBeforeCCurveHistogram, int skip) {
		
		// curve without contrast
		float* dcurve = new float[2*65536];
		
		// check if contrast curve is needed
		bool needcontrast = contr>0.00001 || contr<-0.00001;
		
		// create a curve if needed
		Curve* tcurve = NULL;
		if (curvePoints.size()>0 && curvePoints[0]!=0)
			tcurve = new Curve (curvePoints, CURVES_MIN_POLY_POINTS/skip);
		
		// clear array that stores histogram valid before applying the custom curve
		if (outBeforeCCurveHistogram)
			memset (outBeforeCCurveHistogram, 0, 256*sizeof(int));
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// tone curve base. a: slope (from exp.comp.), b: black, def_mul: max. x value (can be>1), hr,sr: highlight,shadow recovery
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		std::vector<double> brightcurvePoints;
		brightcurvePoints.push_back((double)((CurveType)NURBS));
		
		brightcurvePoints.push_back(0); //black point.  Value in [0 ; 1] range
		brightcurvePoints.push_back(0); //black point.  Value in [0 ; 1] range
		
		if(br>0) {
			brightcurvePoints.push_back(0.1); //toe point
			brightcurvePoints.push_back(0.1+br/150.0); //value at toe point
			
			brightcurvePoints.push_back(0.7); //shoulder point
			brightcurvePoints.push_back(MIN(1.0,0.7+br/300.0)); //value at shoulder point
		} else {
			brightcurvePoints.push_back(0.1-br/150.0); //toe point
			brightcurvePoints.push_back(0.1); //value at toe point
			
			brightcurvePoints.push_back(MIN(1.0,0.7-br/300.0)); //shoulder point
			brightcurvePoints.push_back(0.7); //value at shoulder point
		}
		brightcurvePoints.push_back(1); // white point
		brightcurvePoints.push_back(1); // value at white point
		
		Curve* brightcurve = NULL;
		brightcurve = new Curve (brightcurvePoints, CURVES_MIN_POLY_POINTS/skip); // Actually, CURVES_MIN_POLY_POINTS = 1000,
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		for (int i=0; i<0x10000; i++) {
			
			// change to [0,1] range
			float val = (float)i / 65535.0;
			
			// apply brightness curve
			val = brightcurve->getVal (val);
			
			// store result in a temporary array
			dcurve[i] = CLIPD(val);
		}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		if (needcontrast) {  
			// compute mean luminance of the image with the curve applied
			int sum = 0;
			float avg = 0; 
			//float sqavg = 0;
			for (int i=0; i<0x10000; i++) {
				avg += dcurve[i] * histogram[i];
				//sqavg += dcurve[i]*dcurve[i] * histogram[i];
				sum += histogram[i];
			}
			avg /= sum;
			//sqavg /= sum;
			//float stddev = sqrt(sqavg-avg*avg);
			
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			std::vector<double> contrastcurvePoints;
			contrastcurvePoints.push_back((double)((CurveType)NURBS));
			
			contrastcurvePoints.push_back(0); //black point.  Value in [0 ; 1] range
			contrastcurvePoints.push_back(0); //black point.  Value in [0 ; 1] range
			
			contrastcurvePoints.push_back(avg-avg*(0.6-contr/250.0)); //toe point
			contrastcurvePoints.push_back(avg-avg*(0.6+contr/250.0)); //value at toe point
			
			contrastcurvePoints.push_back(avg+(1-avg)*(0.6-contr/250.0)); //shoulder point
			contrastcurvePoints.push_back(avg+(1-avg)*(0.6+contr/250.0)); //value at shoulder point
			
			contrastcurvePoints.push_back(1); // white point
			contrastcurvePoints.push_back(1); // value at white point
			
			Curve* contrastcurve = NULL;
			contrastcurve = new Curve (contrastcurvePoints, CURVES_MIN_POLY_POINTS/skip); // Actually, CURVES_MIN_POLY_POINTS = 1000,
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			// apply contrast enhancement
			for (int i=0; i<0x10000; i++) {
				dcurve[i]  = contrastcurve->getVal (dcurve[i]);
			}
			delete contrastcurve;
		}
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		for (int i=0; i<0x10000; i++) {
			float val;
			
			// apply custom/parametric/NURBS curve, if any
			if (tcurve) {
				if (outBeforeCCurveHistogram) {
					float hval = dcurve[i];
					int hi = (int)(255.0*CLIPD(hval));
					outBeforeCCurveHistogram[hi]+=histogram[i] ;
				}
				val = tcurve->getVal (dcurve[i]);
			} else {
				val = (dcurve[i]);
			}
			
			outCurve[i] = (65535.0 * val);
		}
		
		
		delete [] dcurve;
		delete tcurve;
		delete brightcurve; 
		if (outBeforeCCurveHistogram) {
			//for (int i=0; i<256; i++) printf("i= %d bchist= %d \n",i,outBeforeCCurveHistogram[i]);
		}
		
	}
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	

int *CurveFactory::gammatab = 0;
int *CurveFactory::igammatab_srgb = 0;
int *CurveFactory::gammatab_srgb = 0;

void CurveFactory::init () {
	
	gammatab = new int[65536];
	igammatab_srgb = new int[65536];
	gammatab_srgb = new int[65536];

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

void CurveFactory::cleanup () {

  delete [] gammatab;
  delete [] igammatab_srgb;
  delete [] gammatab_srgb;
}

}
