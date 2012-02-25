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

#include "array2D.h"
#include "LUT.h"

#undef CLIPD
#define CLIPD(a) ((a)>0.0f?((a)<1.0f?(a):1.0f):0.0f)
#define CLIP(a) ((a)<65535 ? (a) : (65535))


namespace rtengine {

	Curve::Curve () {
		x = 0;
		y = 0;
		ypp = 0;
		hashSize = 1000;  // has to be initiallised to the maximum value
	}
	
	void Curve::AddPolygons ()
	{
		if (firstPointIncluded) {
			poly_x.push_back(x1);
			poly_y.push_back(y1);
		}
		for (int k=1; k<(nbr_points-1); k++) {
			double t = k*increment;
			double t2 = t*t;
			double tr = 1.-t;
			double tr2 = tr*tr;
			double tr2t = tr*2*t;
			
			// adding a point to the polyline
			poly_x.push_back( tr2*x1 + tr2t*x2 + t2*x3);
			poly_y.push_back( tr2*y1 + tr2t*y2 + t2*y3);
		}
		// adding the last point of the sub-curve
		poly_x.push_back(x3);
		poly_y.push_back(y3);
	}

	void Curve::fillHash() {
    	hash.resize(hashSize+2);

    	unsigned int polyIter = 0;
    	double const increment = 1./hashSize;
    	double milestone = 0.;

		for (unsigned short i=0; i<(hashSize+1);) {
    		while(poly_x[polyIter] <= milestone) ++polyIter;
    		hash.at(i).smallerValue = polyIter-1;
    		++i;
    		milestone = i*increment;
    	}
		milestone = 0.;
		polyIter = 0;
		for (unsigned int i=0; i<(hashSize+1);) {
    		while(poly_x[polyIter] < (milestone+increment)) ++polyIter;
    		hash.at(i).higherValue = polyIter;
    		++i;
    		milestone = i*increment;
    	}
    	hash.at(hashSize+1).smallerValue = poly_x.size()-1;
    	hash.at(hashSize+1).higherValue = poly_x.size();

    	/*
    	 * Uncoment the code below to dump the polygon points and the hash table in files
    	if (poly_x.size() > 500) {
    		printf("Files generated (%d points)\n", poly_x.size());
			FILE* f = fopen ("hash.txt", "wt");
			for (unsigned int i=0; i<hashSize;i++) {
				unsigned short s = hash.at(i).smallerValue;
				unsigned short h = hash.at(i).higherValue;
				fprintf (f, "%d: %d<%d (%.5f<%.5f)\n", i, s, h, poly_x[s], poly_x[h]);
			}
			fclose (f);
			f = fopen ("poly_x.txt", "wt");
			for (unsigned int i=0; i<poly_x.size();i++) {
				fprintf (f, "%d: %.5f, %.5f\n", i, poly_x[i], poly_y[i]);
			}
			fclose (f);
    	}
    	*/

	}

    // Wikipedia sRGB: Unlike most other RGB color spaces, the sRGB gamma cannot be expressed as a single numerical value.
    // The overall gamma is approximately 2.2, consisting of a linear (gamma 1.0) section near black, and a non-linear section elsewhere involving a 2.4 exponent 
    // and a gamma (slope of log output versus log input) changing from 1.0 through about 2.3.
    const double CurveFactory::sRGBGamma = 2.2;
    const double CurveFactory::sRGBGammaCurve = 2.4;

	void fillCurveArray(DiagonalCurve* diagCurve, LUTf &outCurve, int skip, bool needed) {
		if (needed) {
			LUTf lutCurve (65536);

			for (int i=0; i<=0xffff; i+= i<0xffff-skip ? skip : 1 ) {
				// change to [0,1] range
				double val = (double)i / 65535.0;
				// apply custom/parametric/NURBS curve, if any
				val = diagCurve->getVal (val);
				// store result in a temporary array
				lutCurve[i] = (val);
			}

			// if skip>1, let apply linear interpolation in the skipped points of the curve
			if (skip > 1) {
				int prev = 0;
				for (int i=1; i<=0xffff-skip; i++) {
					if (i%skip==0) {
						prev+=skip;
						continue;
					}
					lutCurve[i] = ( lutCurve[prev] * (skip - i%skip) + lutCurve[prev+skip] * (i%skip) ) / skip;
				}
			}

			for (int i=0; i<=0xffff; i++) {
				outCurve[i] = (65535.0 * lutCurve[i]);
			}
		}
		else {
			for (int i=0; i<=0xffff; i++) {
				outCurve[i] = (float)i;
			}
		}
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	void CurveFactory::complexsgnCurve (double saturation, bool satlimit, double satlimthresh,
										const std::vector<double>& acurvePoints, const std::vector<double>& bcurvePoints,
										LUTf & aoutCurve, LUTf & boutCurve, LUTf & satCurve, int skip) {
		
		//colormult = chroma_scale for Lab manipulations
		
		//-----------------------------------------------------

		bool needed;
		DiagonalCurve* dCurve = NULL;

		// check if contrast curve is needed
		needed = (saturation<-0.0001 || saturation>0.0001);

		// Filling the curve if needed
		if (needed) {

			//%%%%%%%%%%%%%%%%% Saturation curve's control points %%%%%%%%%%%%%%%%%
		    std::vector<double> satcurvePoints;
			satcurvePoints.push_back((double)DCT_NURBS);
			if (saturation>0) {
				double satslope = (0.5+2*saturation/500.0)/(0.5-2*saturation/500.0);
				double scale = (satlimthresh/100.1);
				if (!satlimit) scale=100/100.1;

				satcurvePoints.push_back(0); //black point.  Value in [0 ; 1] range
				satcurvePoints.push_back(0); //black point.  Value in [0 ; 1] range

				//if (satlimit) {
				satcurvePoints.push_back(0.5-0.5*scale); //toe point
				satcurvePoints.push_back(0.5-0.5*scale); //value at toe point

				satcurvePoints.push_back(0.5-(0.5/satslope)*scale); //toe point
				satcurvePoints.push_back(0.5-0.5*scale); //value at toe point

				satcurvePoints.push_back(0.5+(0.5/satslope)*scale); //shoulder point
				satcurvePoints.push_back(0.5+0.5*scale); //value at shoulder point

				satcurvePoints.push_back(0.5+0.5*scale); //shoulder point
				satcurvePoints.push_back(0.5+0.5*scale); //value at shoulder point
				/*} else {
				 satcurvePoints.push_back(0.25+saturation/500.0); //toe point
				 satcurvePoints.push_back(0.25-saturation/500.0); //value at toe point

				 satcurvePoints.push_back(0.75-saturation/500.0); //shoulder point
				 satcurvePoints.push_back(0.75+saturation/500.0); //value at shoulder point
				 }*/

				satcurvePoints.push_back(1); // white point
				satcurvePoints.push_back(1); // value at white point
			} else {
				satcurvePoints.push_back(0);
				satcurvePoints.push_back(-(saturation/200.0));

				satcurvePoints.push_back(1);
				satcurvePoints.push_back(1+saturation/200.0);
			}
			dCurve = new DiagonalCurve (satcurvePoints, CURVES_MIN_POLY_POINTS/skip);
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			fillCurveArray(dCurve, satCurve, skip, needed);

			delete dCurve;
			dCurve = NULL;
		}
		else {
			fillCurveArray(NULL, satCurve, skip, needed);
		}

		//-----------------------------------------------------

		needed = false;
		// create a curve if needed
		if (!acurvePoints.empty() && acurvePoints[0]!=0) {
			dCurve = new DiagonalCurve (acurvePoints, CURVES_MIN_POLY_POINTS/skip);
			if (dCurve && !dCurve->isIdentity())
				needed = true;
		}
		fillCurveArray(dCurve, aoutCurve, skip, needed);
		if (dCurve) {
			delete dCurve;
			dCurve = NULL;
		}

		//-----------------------------------------------------

		needed = false;
		if (!bcurvePoints.empty() && bcurvePoints[0]!=0) {
			dCurve = new DiagonalCurve (bcurvePoints, CURVES_MIN_POLY_POINTS/skip);
			if (dCurve && !dCurve->isIdentity())
				needed = true;
		}
		fillCurveArray(dCurve, boutCurve, skip, needed);
		if (dCurve) {
			delete dCurve;
			dCurve = NULL;
		}
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	void CurveFactory::complexCurve (double ecomp, double black, double hlcompr, double hlcomprthresh,
									 double shcompr, double br, double contr, double gamma_, bool igamma_,
									 const std::vector<double>& curvePoints, LUTu & histogram, LUTu & histogramCropped,
									 LUTf & hlCurve, LUTf & shCurve, LUTf & outCurve,
									 LUTu & outBeforeCCurveHistogram, int skip) {
		
		
		//double def_mul = pow (2.0, defmul);
		
		/*printf ("def_mul= %f ecomp= %f black= %f  hlcompr= %f shcompr= %f br= %f contr= %f defmul= %f
				gamma= %f, skip= %d \n",def_mul,ecomp,black,hlcompr,shcompr,br,contr,defmul,gamma_,skip);*/
		
		// compute parameters of the gamma curve
		/*double start = exp(gamma_*log( -0.099 / ((1.0/gamma_-1.0)*1.099 )));
		double slope = 1.099 * pow (start, 1.0/gamma_-1) - 0.099/start;
		double mul = 1.099;
		double add = 0.099;
		// gamma BT709*/
		
		//normalize gamma to sRGB
		double start = exp(gamma_*log( -0.055 / ((1.0/gamma_-1.0)*1.055 )));
		double slope = 1.055 * pow (start, 1.0/gamma_-1) - 0.055/start;
		double mul = 1.055;
		double add = 0.055;
		
		// a: slope of the curve, black: starting point at the x axis
		double a = pow (2.0, ecomp);
		
		// curve without contrast
		LUTf dcurve(0x10000);
		
		// check if inverse gamma is needed at the end
		bool needigamma = igamma_ && gamma_>0;
		
		// clear array that stores histogram valid before applying the custom curve
		outBeforeCCurveHistogram.clear();
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// tone curve base. a: slope (from exp.comp.), b: black, def_mul: max. x value (can be>1), hr,sr: highlight,shadow recovery
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		DiagonalCurve* brightcurve = NULL;

		// check if brightness curve is needed
		if (br>0.00001 || br<-0.00001) {

			std::vector<double> brightcurvePoints;
			brightcurvePoints.push_back((double)DCT_NURBS);

			brightcurvePoints.push_back(0.); //black point.  Value in [0 ; 1] range
			brightcurvePoints.push_back(0.); //black point.  Value in [0 ; 1] range
			
			if(br>0) {
				brightcurvePoints.push_back(0.1); //toe point
				brightcurvePoints.push_back(0.1+br/150.0); //value at toe point

				brightcurvePoints.push_back(0.7); //shoulder point
				brightcurvePoints.push_back(MIN(1.0,0.7+br/300.0)); //value at shoulder point
			} else {
				brightcurvePoints.push_back(MAX(0.0,0.1-br/150.0)); //toe point
				brightcurvePoints.push_back(0.1); //value at toe point

				brightcurvePoints.push_back(0.7-br/300.0); //shoulder point
				brightcurvePoints.push_back(0.7); //value at shoulder point
			}
			brightcurvePoints.push_back(1.); // white point
			brightcurvePoints.push_back(1.); // value at white point
			
			brightcurve = new DiagonalCurve (brightcurvePoints, CURVES_MIN_POLY_POINTS/skip);
		}
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		float exp_scale = a;
		float scale = 65536.0;
		float comp = (MAX(0,ecomp) + 1.0)*hlcompr/100.0;
		float shoulder = ((scale/MAX(1,exp_scale))*(hlcomprthresh/200.0))+0.1;
		//printf("shoulder = %e\n",shoulder);
		//printf ("exp_scale= %f comp= %f def_mul=%f a= %f \n",exp_scale,comp,def_mul,a);
		
		for (int i=0; i<0x10000; i++) {
			
			// change to [0,1] range
			float val = (float)i-shoulder;
			
			if (comp>0.0)
			{
				if (val>0.0) {
					float R = val*comp/(scale-shoulder);
					hlCurve[i] = log(1.0+R*exp_scale)/R;
				} else {
					hlCurve[i]=exp_scale;
				}
			} else {
				hlCurve[i]=exp_scale;
			}
			
			//%%%%%%%%%%%%%%%%%%%%%%%%%%
			// change to [0,1] range
			if (i!=0) {
				val = (float)i / 65535.0f;
			} else {
				val = 1.0/65535.0;
			}
			
			float	val2 = basecurve (val, 1.0, black, 1.0, 0.0, 1.5*shcompr/100.0);
			shCurve[i] = CLIPD(val2)/val;

			//%%%%%%%%%%%%%%%%%%%%%%%%%%
			// change to [0,1] range
			val = (double)i / 65535.0;
			
			// gamma correction
			if (gamma_>1)
				val = gamma (val, gamma_, start, slope, mul, add);
			
			// apply brightness curve
			if (brightcurve)
				val = brightcurve->getVal (val);  // TODO: getVal(double) is very slow! Optimize with a LUTf

			// store result in a temporary array
			dcurve[i] = CLIPD(val);
		}

		if (brightcurve)
			delete brightcurve;

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		// check if contrast curve is needed
		if (contr>0.00001 || contr<-0.00001) {

			// compute mean luminance of the image with the curve applied
			int sum = 0;
			float avg = 0; 
			//double sqavg = 0;
			for (int i=0; i<=0xffff; i++) {
				float fi=i;
				fi = hlCurve[fi]*fi;
				avg += dcurve[(int)(shCurve[fi]*fi)] * histogram[i];
				//sqavg += dcurve[i]*dcurve[i] * histogram[i];
				sum += histogram[i];
			}
			avg /= sum;
			//sqavg /= sum;
			//double stddev = sqrt(sqavg-avg*avg);
			
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			std::vector<double> contrastcurvePoints;
			contrastcurvePoints.push_back((double)DCT_NURBS);
			
			contrastcurvePoints.push_back(0); //black point.  Value in [0 ; 1] range
			contrastcurvePoints.push_back(0); //black point.  Value in [0 ; 1] range
			
			contrastcurvePoints.push_back(avg-avg*(0.6-contr/250.0)); //toe point
			contrastcurvePoints.push_back(avg-avg*(0.6+contr/250.0)); //value at toe point
			
			contrastcurvePoints.push_back(avg+(1-avg)*(0.6-contr/250.0)); //shoulder point
			contrastcurvePoints.push_back(avg+(1-avg)*(0.6+contr/250.0)); //value at shoulder point
			
			contrastcurvePoints.push_back(1); // white point
			contrastcurvePoints.push_back(1); // value at white point
			
			DiagonalCurve* contrastcurve = new DiagonalCurve (contrastcurvePoints, CURVES_MIN_POLY_POINTS/skip);
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			// apply contrast enhancement
			for (int i=0; i<=0xffff; i++) {
				dcurve[i]  = contrastcurve->getVal (dcurve[i]);
			}

			delete contrastcurve;
		}

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		// create a curve if needed
		bool histNeeded = false;
		DiagonalCurve* tcurve = NULL;
		if (!curvePoints.empty() && curvePoints[0]!=0) {
			tcurve = new DiagonalCurve (curvePoints, CURVES_MIN_POLY_POINTS/skip);
			if (outBeforeCCurveHistogram /*&& histogramCropped*/)
				histNeeded = true;
		}
		if (tcurve && tcurve->isIdentity()) {
			delete tcurve;
			tcurve = NULL;
		}

		for (int i=0; i<=0xffff; i++) {
			float val;

			if (histNeeded) {
				float fi=i;
				float hval = hlCurve[i]*fi;
				hval = dcurve[shCurve[hval]*hval];
				//if (needigamma)
				//	hval = igamma2 (hval);
				int hi = (int)(255.0*(hval));
				outBeforeCCurveHistogram[hi] += histogram/*Cropped*/[i] ;
			}

			// apply custom/parametric/NURBS curve, if any
			if (tcurve) {
				val = tcurve->getVal (dcurve[i]);  // TODO: getVal(double) is very slow! Optimize with a LUTf
			} else {
				val = (dcurve[i]);
			}

			// if inverse gamma is needed, do it (standard sRGB inverse gamma is applied)
			if (needigamma)
				val = igamma (val, gamma_, start, slope, mul, add);

			outCurve[i] = (65535.0 * val);
		}

		if (tcurve)
			delete tcurve;

		/*if (outBeforeCCurveHistogram) {
			for (int i=0; i<256; i++) printf("i= %d bchist= %d \n",i,outBeforeCCurveHistogram[i]);
		}*/
	}

	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void CurveFactory::complexLCurve (double br, double contr, const std::vector<double>& curvePoints,
									 LUTu & histogram, LUTu & histogramCropped, LUTf & outCurve,
									 LUTu & outBeforeCCurveHistogram, int skip) {
		
		// curve without contrast
		LUTf dcurve(65536,0);
		
		// clear array that stores histogram valid before applying the custom curve
		if (outBeforeCCurveHistogram)
			outBeforeCCurveHistogram.clear();
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// tone curve base. a: slope (from exp.comp.), b: black, def_mul: max. x value (can be>1), hr,sr: highlight,shadow recovery
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		// check if brightness curve is needed
		if (br>0.00001 || br<-0.00001) {

			std::vector<double> brightcurvePoints;
			brightcurvePoints.push_back((double)((CurveType)DCT_NURBS));

			brightcurvePoints.push_back(0.); // black point.  Value in [0 ; 1] range
			brightcurvePoints.push_back(0.); // black point.  Value in [0 ; 1] range
			
			if (br>0) {
				brightcurvePoints.push_back(0.1); // toe point
				brightcurvePoints.push_back(0.1+br/150.0); //value at toe point

				brightcurvePoints.push_back(0.7); // shoulder point
				brightcurvePoints.push_back(MIN(1.0,0.7+br/300.0)); //value at shoulder point
			} else {
				brightcurvePoints.push_back(0.1-br/150.0); // toe point
				brightcurvePoints.push_back(0.1); // value at toe point

				brightcurvePoints.push_back(MIN(1.0,0.7-br/300.0)); // shoulder point
				brightcurvePoints.push_back(0.7); // value at shoulder point
			}
			brightcurvePoints.push_back(1.); // white point
			brightcurvePoints.push_back(1.); // value at white point
			
			DiagonalCurve* brightcurve = new DiagonalCurve (brightcurvePoints, CURVES_MIN_POLY_POINTS/skip);
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			// Applying brightness curve
			for (int i=0; i<32768; i++) { // L values range up to 32767, higher values are for highlight overflow

				// change to [0,1] range
				float val = (float)i / 32767.0;

				// apply brightness curve
				val = brightcurve->getVal (val);

				// store result in a temporary array
				dcurve[i] = CLIPD(val);
			}
			delete brightcurve;
		}
		else {
			for (int i=0; i<32768; i++) { // L values range up to 32767, higher values are for highlight overflow
				// set the identity curve in the temporary array
				dcurve[i] = (float)i / 32767.0;
			}
		}
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		// check if contrast curve is needed
		if (contr>0.00001 || contr<-0.00001) {

			// compute mean luminance of the image with the curve applied
			int sum = 0;
			float avg = 0; 
			//float sqavg = 0;
			for (int i=0; i<32768; i++) {
				avg += dcurve[i] * histogram[i];
				//sqavg += dcurve[i]*dcurve[i] * histogram[i];
				sum += histogram[i];
			}
			avg /= sum;
			//sqavg /= sum;
			//float stddev = sqrt(sqavg-avg*avg);
			
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			std::vector<double> contrastcurvePoints;
			contrastcurvePoints.push_back((double)((CurveType)DCT_NURBS));
			
			contrastcurvePoints.push_back(0.); // black point.  Value in [0 ; 1] range
			contrastcurvePoints.push_back(0.); // black point.  Value in [0 ; 1] range
			
			contrastcurvePoints.push_back(avg-avg*(0.6-contr/250.0)); // toe point
			contrastcurvePoints.push_back(avg-avg*(0.6+contr/250.0)); // value at toe point
			
			contrastcurvePoints.push_back(avg+(1-avg)*(0.6-contr/250.0)); // shoulder point
			contrastcurvePoints.push_back(avg+(1-avg)*(0.6+contr/250.0)); // value at shoulder point
			
			contrastcurvePoints.push_back(1.); // white point
			contrastcurvePoints.push_back(1.); // value at white point
			
			DiagonalCurve* contrastcurve = new DiagonalCurve (contrastcurvePoints, CURVES_MIN_POLY_POINTS/skip);
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			// apply contrast enhancement
			for (int i=0; i<32768; i++) {
				dcurve[i]  = contrastcurve->getVal (dcurve[i]);
			}

			delete contrastcurve;
		}
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		// create a curve if needed
		DiagonalCurve* tcurve = NULL;
		bool histNeeded = false;
		if (!curvePoints.empty() && curvePoints[0]!=0) {
			tcurve = new DiagonalCurve (curvePoints, CURVES_MIN_POLY_POINTS/skip);
			if (outBeforeCCurveHistogram /*&& histogramCropped*/)
				histNeeded = true;
		}
		if (tcurve && tcurve->isIdentity()) {
			delete tcurve;
			tcurve = NULL;
		}

		if (tcurve) {
			// L values go up to 32767, last stop is for highlight overflow
			for (int i=0; i<32768; i++) {
				float val;

				if (histNeeded) {
					float hval = dcurve[i];
					int hi = (int)(255.0*CLIPD(hval));
					outBeforeCCurveHistogram[hi]+=histogram/*Cropped*/[i] ;
				}

				// apply custom/parametric/NURBS curve, if any
				val = tcurve->getVal (dcurve[i]);

				outCurve[i] = (32767.0 * val);
			}
		}
		else {
			// Skip the slow getval method if no curve is used (or an identity curve)
			// L values go up to 32767, last stop is for highlight overflow
			for (int i=0; i<32768; i++) {
				if (histNeeded) {
					float hval = dcurve[i];
					int hi = (int)(255.0*CLIPD(hval));
					outBeforeCCurveHistogram[hi]+=histogram/*Cropped*/[i] ;
				}

				outCurve[i] = 32767.0*dcurve[i];
			}
		}
		for (int i=32768; i<65536; i++) outCurve[i]=(float)i;

		if (tcurve)
			delete tcurve;

		/*if (outBeforeCCurveHistogram) {
			for (int i=0; i<256; i++) printf("i= %d bchist= %d \n",i,outBeforeCCurveHistogram[i]);
		}*/
	}


	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void CurveFactory::RGBCurve (const std::vector<double>& curvePoints, LUTf & outCurve, int skip) {
		
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						
		// create a curve if needed
		DiagonalCurve* tcurve = NULL;
		if (!curvePoints.empty() && curvePoints[0]!=0) {
			tcurve = new DiagonalCurve (curvePoints, CURVES_MIN_POLY_POINTS/skip);
		}
		if (tcurve && tcurve->isIdentity()) {
			delete tcurve;
			tcurve = NULL;
		}
		
		if (tcurve) {
			for (int i=0; i<65536; i++) {				
				// apply custom/parametric/NURBS curve, if any
				float val = tcurve->getVal ((float)i/65536.0f);
				outCurve[i] = (65536.0f * val);
			}
		}
		else {
			// Skip the slow getval method if no curve is used (or an identity curve)
			for (int i=0; i<65536; i++) {
				outCurve[i] = i;
			}
		}
		
		if (tcurve)
			delete tcurve;

	}
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	

LUTf CurveFactory::gammatab;
LUTf CurveFactory::igammatab_srgb;
LUTf CurveFactory::gammatab_srgb;

void CurveFactory::init () {
	
	gammatab(65536,0);
	igammatab_srgb(65536,0);
	gammatab_srgb(65536,0);

  for (int i=0; i<65536; i++)
    gammatab_srgb[i] = (65535.0 * gamma2 (i/65535.0));
  for (int i=0; i<65536; i++)
    igammatab_srgb[i] = (65535.0 * igamma2 (i/65535.0));
  for (int i=0; i<65536; i++)
    gammatab[i] = (65535.0 * pow (i/65535.0, 0.454545));
    
/*    FILE* f = fopen ("c.txt", "wt");
    for (int i=0; i<256; i++)
        fprintf (f, "%g %g\n", i/255.0, clower (i/255.0, 2.0, 1.0));
    fclose (f);*/
}

}
