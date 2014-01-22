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
#include <cmath>
#include <vector>
#include <cstring>
#include <algorithm>

#include "rt_math.h"

#include "mytime.h"
#include "array2D.h"
#include "LUT.h"
#include "curves.h"

#undef CLIPD
#define CLIPD(a) ((a)>0.0f?((a)<1.0f?(a):1.0f):0.0f)

using namespace std;

namespace rtengine {

	Curve::Curve () : N(0), x(NULL), y(NULL), ypp(NULL), hashSize(1000 /* has to be initialized to the maximum value */ ) {}
	
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
			for (size_t i=0; i<poly_x.size();i++) {
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

void CurveFactory::updatechroma (
		const std::vector<double>& cccurvePoints,
		LUTu & histogramC, LUTu & outBeforeCCurveHistogramC,//for chroma
		int skip)
{
	LUTf dCcurve(65536,0);
	float val;
	for (int i=0; i<48000; i++) {//32768*1.414  + ...
			val = (double)i / 47999.0;
			dCcurve[i] = CLIPD(val);
		}

	outBeforeCCurveHistogramC.clear();
	bool histNeededC = false;
	

	if (!cccurvePoints.empty() && cccurvePoints[0]!=0) {
		if (outBeforeCCurveHistogramC /*&& histogramCropped*/)
				histNeededC = true;	
	}
	for (int i=0; i<=48000; i++) {//32768*1.414  + ...
		float val;
			if (histNeededC) {
				float hval = dCcurve[i];
				int hi = (int)(255.0*CLIPD(hval)); //
				outBeforeCCurveHistogramC[hi] += histogramC[i] ;
			}
	}
}
	
	
void CurveFactory::curveLightBrightColor (
		procparams::ColorAppearanceParams::eTCModeId curveMode1, const std::vector<double>& curvePoints1,
		procparams::ColorAppearanceParams::eTCModeId curveMode2, const std::vector<double>& curvePoints2,
		procparams::ColorAppearanceParams::eCTCModeId curveMode3, const std::vector<double>& curvePoints3,
		LUTu & histogram, LUTu & histogramCropped, LUTu & outBeforeCCurveHistogram,//for Luminance  
		LUTu & histogramC, LUTu & outBeforeCCurveHistogramC,//for chroma
		ColorAppearance & customColCurve1,
		ColorAppearance & customColCurve2,
		ColorAppearance & customColCurve3,
		int skip)
{
	LUTf dcurve(65536,0);
	LUTf dCcurve(65536,0);
	
	float val;
	for (int i=0; i<32768; i++) {
			val = (double)i / 32767.0;
			dcurve[i] = CLIPD(val);
		}
	for (int i=0; i<48000; i++) {  //# 32768*1.414  approximation maxi for chroma
			val = (double)i / 47999.0;
			dCcurve[i] = CLIPD(val);
		}

	outBeforeCCurveHistogram.clear();
	outBeforeCCurveHistogramC.clear();
	bool histNeededC = false;
	
	bool histNeeded = false;
	DiagonalCurve* tcurve = NULL;
	customColCurve3.Reset();

	if (!curvePoints3.empty() && curvePoints3[0]>DCT_Linear && curvePoints3[0]<DCT_Unchanged) {
		tcurve = new DiagonalCurve (curvePoints3, CURVES_MIN_POLY_POINTS/skip);
		if (outBeforeCCurveHistogramC /*&& histogramCropped*/)
				histNeededC = true;
		
	}
	if (tcurve) {
		if (tcurve->isIdentity()) {
			delete tcurve;
			tcurve = NULL;
		}
		else
			customColCurve3.Set(tcurve);
		delete tcurve;
		tcurve = NULL;
	}

	customColCurve2.Reset();

	if (!curvePoints2.empty() && curvePoints2[0]>DCT_Linear && curvePoints2[0]<DCT_Unchanged) {
		tcurve = new DiagonalCurve (curvePoints2, CURVES_MIN_POLY_POINTS/skip);
		if (outBeforeCCurveHistogram /*&& histogramCropped*/)
				histNeeded = true;
		
	}
	if (tcurve) {
		if (tcurve->isIdentity()) {
			delete tcurve;
			tcurve = NULL;
		}
		else
			customColCurve2.Set(tcurve);
		delete tcurve;
		tcurve = NULL;
	}
	// create first curve if needed
	customColCurve1.Reset();

	if (!curvePoints1.empty() && curvePoints1[0]>DCT_Linear && curvePoints1[0]<DCT_Unchanged) {
		tcurve = new DiagonalCurve (curvePoints1, CURVES_MIN_POLY_POINTS/skip);
		if (outBeforeCCurveHistogram /*&& histogramCropped*/)
				histNeeded = true;
		
	}
	if (tcurve) {
		if (tcurve->isIdentity()) {
			delete tcurve;
			tcurve = NULL;
		}
		else  {
			customColCurve1.Set(tcurve);
			delete tcurve;
			tcurve = NULL;
		}
	}
	for (int i=0; i<=32768; i++) {
		float val;

			if (histNeeded) {
				float hval = dcurve[i];
				int hi = (int)(255.0*CLIPD(hval));
				outBeforeCCurveHistogram[hi] += histogram[i] ;
			}
	}
	for (int i=0; i<=48000; i++) {//32768*1.414  + ...
		float val;
			if (histNeededC) {
				float hval = dCcurve[i];
				int hi = (int)(255.0*CLIPD(hval)); //
				outBeforeCCurveHistogramC[hi] += histogramC[i] ;
			}
	}
	
	if (tcurve) delete tcurve;

}

void CurveFactory::curveBW (
		const std::vector<double>& curvePointsbw, const std::vector<double>& curvePointsbw2,
		LUTu & histogrambw, LUTu & outBeforeCCurveHistogrambw,//for Luminance
		ToneCurve & customToneCurvebw1, ToneCurve & customToneCurvebw2, int skip)
{
	LUTf dcurve(65536,0);
	
	float val;
	for (int i=0; i<32768; i++) {
			val = (double)i / 32767.0;
			dcurve[i] = CLIPD(val);
		}

	outBeforeCCurveHistogrambw.clear();
	bool histNeeded = false;
	
	DiagonalCurve* tcurve = NULL;
	customToneCurvebw2.Reset();

	if (!curvePointsbw2.empty() && curvePointsbw2[0]>DCT_Linear && curvePointsbw2[0]<DCT_Unchanged) {
		tcurve = new DiagonalCurve (curvePointsbw2, CURVES_MIN_POLY_POINTS/skip);
		if (outBeforeCCurveHistogrambw /*&& histogramCropped*/)
				histNeeded = true;
		
	}
	if (tcurve) {
		if (!tcurve->isIdentity())
			customToneCurvebw2.Set(tcurve);
		delete tcurve;
		tcurve = NULL;
	}

	customToneCurvebw1.Reset();

	if (!curvePointsbw.empty() && curvePointsbw[0]>DCT_Linear && curvePointsbw[0]<DCT_Unchanged) {
		tcurve = new DiagonalCurve (curvePointsbw, CURVES_MIN_POLY_POINTS/skip);
		if (outBeforeCCurveHistogrambw /*&& histogramCropped*/)
				histNeeded = true;
		
	}
	if (tcurve) {
		if (!tcurve->isIdentity())
			customToneCurvebw1.Set(tcurve);
		delete tcurve;
		tcurve = NULL;
	}
	// create first curve if needed

	for (int i=0; i<=32768; i++) {
		float val;

			if (histNeeded) {
				float hval = dcurve[i];
				int hi = (int)(255.0*CLIPD(hval));
				outBeforeCCurveHistogrambw[hi] += histogrambw[i] ;
			}
	}
	
	if (tcurve) delete tcurve;

}
// add curve Lab : C=f(L)
void CurveFactory::curveCL ( bool & clcutili,const std::vector<double>& clcurvePoints, LUTf & clCurve, LUTu & histogramcl, LUTu & outBeforeCLurveHistogram,int skip){
		bool needed;
		DiagonalCurve* dCurve = NULL;
		LUTf dCcurve(65536,0);
	
		float val;
		for (int i=0; i<50000; i++) {  //# 32768*1.414  approximation maxi for chroma
				dCcurve[i] = (float)i / 49999.0;
		}
		
		if (outBeforeCLurveHistogram)		
			outBeforeCLurveHistogram.clear();
		bool histNeededCL = false;

		needed = false;
		if (!clcurvePoints.empty() && clcurvePoints[0]!=0) {
			dCurve = new DiagonalCurve (clcurvePoints, CURVES_MIN_POLY_POINTS/skip);
			if (outBeforeCLurveHistogram)
				histNeededCL = true;
			
			if (dCurve && !dCurve->isIdentity())
				{needed = true;clcutili=true;}
		}
		for (int i=0; i<=50000; i++) {//32768*1.414  + ...
			float val;
			if (histNeededCL) {
				float hval = dCcurve[i];
				int hi = (int)(255.0*CLIPD(hval)); //
				outBeforeCLurveHistogram[hi] += histogramcl[i] ;
			}
		}
		
		fillCurveArray(dCurve, clCurve, skip, needed);
		if (dCurve) {
			delete dCurve;
			dCurve = NULL;
		}
}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	void CurveFactory::complexsgnCurve ( bool & autili,  bool & butili, bool & ccutili, bool & cclutili, double saturation, double rstprotection,
										const std::vector<double>& acurvePoints, const std::vector<double>& bcurvePoints,const std::vector<double>& cccurvePoints,
										const std::vector<double>& lccurvePoints, LUTf & aoutCurve, LUTf & boutCurve, LUTf & satCurve, LUTf & lhskCurve,
										LUTu & histogramC, LUTu & histogramLC, LUTu & outBeforeCCurveHistogram,LUTu & outBeforeLCurveHistogram, //for chroma
										int skip) {
		
		
		//-----------------------------------------------------

		bool needed;
		DiagonalCurve* dCurve = NULL;
		LUTf dCcurve(65536,0);
	
		for (int i=0; i<48000; i++) {  //# 32768*1.414  approximation maxi for chroma
				dCcurve[i] = (float)i / 47999.0;
		}
		if (outBeforeCCurveHistogram)		
			outBeforeCCurveHistogram.clear();
		bool histNeededC = false;
		
		if (outBeforeLCurveHistogram)		
			outBeforeLCurveHistogram.clear();
		bool histNeededLC = false;
		
		//-----------------------------------------------------

		needed = false;
		// create a curve if needed
		if (!acurvePoints.empty() && acurvePoints[0]!=0) {
			dCurve = new DiagonalCurve (acurvePoints, CURVES_MIN_POLY_POINTS/skip);
			if (dCurve && !dCurve->isIdentity()) {
				needed = true;
				autili=true;
			}
		}
		fillCurveArray(dCurve, aoutCurve, skip, needed);
		//if(autili) aoutCurve.dump("acurve");

		if (dCurve) {
			delete dCurve;
			dCurve = NULL;
		}

		//-----------------------------------------------------

		needed = false;
		if (!bcurvePoints.empty() && bcurvePoints[0]!=0) {
			dCurve = new DiagonalCurve (bcurvePoints, CURVES_MIN_POLY_POINTS/skip);
			if (dCurve && !dCurve->isIdentity()) {
				needed = true;
				butili=true;
			}
		}
		fillCurveArray(dCurve, boutCurve, skip, needed);
		if (dCurve) {
			delete dCurve;
			dCurve = NULL;
		}
		
		//-----------------------------------------------
		needed = false;
		if (!cccurvePoints.empty() && cccurvePoints[0]!=0) {
			dCurve = new DiagonalCurve (cccurvePoints, CURVES_MIN_POLY_POINTS/skip);
			if (outBeforeCCurveHistogram /*&& histogramCropped*/)
				histNeededC = true;
			
			if (dCurve && !dCurve->isIdentity())
				{needed = true;ccutili=true;}
		}
		for (int i=0; i<=48000; i++) {//32768*1.414  + ...
			float val;
			if (histNeededC) {
				float hval = dCcurve[i];
				int hi = (int)(255.0*CLIPD(hval)); //
				outBeforeCCurveHistogram[hi] += histogramC[i] ;
			}
		}
		
		fillCurveArray(dCurve, satCurve, skip, needed);
		if (dCurve) {
			delete dCurve;
			dCurve = NULL;
		}
		//----------------------------
		needed = false;
		if (!lccurvePoints.empty() && lccurvePoints[0]!=0) {
			dCurve = new DiagonalCurve (lccurvePoints, CURVES_MIN_POLY_POINTS/skip);
			if (outBeforeLCurveHistogram /*&& histogramCropped*/)
				histNeededLC = true;
			
			if (dCurve && !dCurve->isIdentity())
				{needed = true;cclutili=true;}
		}
		for (int i=0; i<=48000; i++) {//32768*1.414  + ...
			float val;
			if (histNeededLC) {
				float hval = dCcurve[i];
				int hi = (int)(255.0*CLIPD(hval)); //
				outBeforeLCurveHistogram[hi] += histogramLC[i] ;
			}
		}
		
		
		fillCurveArray(dCurve, lhskCurve, skip, needed);
		if (dCurve) {
			delete dCurve;
			dCurve = NULL;
		}
		
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	void CurveFactory::complexCurve (double ecomp, double black, double hlcompr, double hlcomprthresh,
									 double shcompr, double br, double contr, double gamma_, bool igamma_,
									 procparams::ToneCurveParams::eTCModeId curveMode, const std::vector<double>& curvePoints,
									 procparams::ToneCurveParams::eTCModeId curveMode2, const std::vector<double>& curvePoints2,
									 LUTu & histogram, LUTu & histogramCropped,
									 LUTf & hlCurve, LUTf & shCurve, LUTf & outCurve,
									 LUTu & outBeforeCCurveHistogram,
									 ToneCurve & customToneCurve1,
									 ToneCurve & customToneCurve2,
									 
									 int skip) {
		
		
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
		bool needigamma = igamma_ && gamma_>1.;
		
		// clear array that stores histogram valid before applying the custom curve
		outBeforeCCurveHistogram.clear();
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// tone curve base. a: slope (from exp.comp.), b: black, def_mul: max. x value (can be>1), hr,sr: highlight,shadow recovery
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		DiagonalCurve* brightcurve = NULL;

		// check if brightness curve is needed
		if (br>0.00001 || br<-0.00001) {

			std::vector<double> brightcurvePoints;
			brightcurvePoints.resize(9);
			brightcurvePoints.at(0) = double(DCT_NURBS);

			brightcurvePoints.at(1) = 0.; //black point.  Value in [0 ; 1] range
			brightcurvePoints.at(2) = 0.; //black point.  Value in [0 ; 1] range
			
			if(br>0) {
				brightcurvePoints.at(3) = 0.1; //toe point
				brightcurvePoints.at(4) = 0.1+br/150.0; //value at toe point

				brightcurvePoints.at(5) = 0.7; //shoulder point
				brightcurvePoints.at(6) = min(1.0,0.7+br/300.0); //value at shoulder point
			} else {
				brightcurvePoints.at(3) = max(0.0,0.1-br/150.0); //toe point
				brightcurvePoints.at(4) = 0.1; //value at toe point

				brightcurvePoints.at(5) = 0.7-br/300.0; //shoulder point
				brightcurvePoints.at(6) = 0.7; //value at shoulder point
			}
			brightcurvePoints.at(7) = 1.; // white point
			brightcurvePoints.at(8) = 1.; // value at white point
			
			brightcurve = new DiagonalCurve (brightcurvePoints, CURVES_MIN_POLY_POINTS/skip);
		}
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		float exp_scale = a;
		float scale = 65536.0;
		float comp = (max(0.0,ecomp) + 1.0)*hlcompr/100.0;
		float shoulder = ((scale/max(1.0f,exp_scale))*(hlcomprthresh/200.0))+0.1;
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
			if (gamma_>1.)
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
			contrastcurvePoints.resize(9);
			contrastcurvePoints.at(0) = double(DCT_NURBS);
			
			contrastcurvePoints.at(1) = 0; //black point.  Value in [0 ; 1] range
			contrastcurvePoints.at(2) = 0; //black point.  Value in [0 ; 1] range
			
			contrastcurvePoints.at(3) = avg-avg*(0.6-contr/250.0); //toe point
			contrastcurvePoints.at(4) = avg-avg*(0.6+contr/250.0); //value at toe point
			
			contrastcurvePoints.at(5) = avg+(1-avg)*(0.6-contr/250.0); //shoulder point
			contrastcurvePoints.at(6) = avg+(1-avg)*(0.6+contr/250.0); //value at shoulder point
			
			contrastcurvePoints.at(7) = 1.; // white point
			contrastcurvePoints.at(8) = 1.; // value at white point
			
			DiagonalCurve* contrastcurve = new DiagonalCurve (contrastcurvePoints, CURVES_MIN_POLY_POINTS/skip);
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			// apply contrast enhancement
			for (int i=0; i<=0xffff; i++) {
				dcurve[i]  = contrastcurve->getVal (dcurve[i]);
			}

			delete contrastcurve;
		}

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		// create second curve if needed
		bool histNeeded = false;
		DiagonalCurve* tcurve = NULL;
		customToneCurve2.Reset();

		if (!curvePoints2.empty() && curvePoints2[0]>DCT_Linear && curvePoints2[0]<DCT_Unchanged) {
			tcurve = new DiagonalCurve (curvePoints2, CURVES_MIN_POLY_POINTS/skip);
			if (outBeforeCCurveHistogram /*&& histogramCropped*/)
				histNeeded = true;
		}
		if (tcurve) {
			if (!tcurve->isIdentity())
				customToneCurve2.Set(tcurve);
			delete tcurve;
			tcurve = NULL;
		}

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		// create first curve if needed
		customToneCurve1.Reset();

		if (!curvePoints.empty() && curvePoints[0]>DCT_Linear && curvePoints[0]<DCT_Unchanged) {
			tcurve = new DiagonalCurve (curvePoints, CURVES_MIN_POLY_POINTS/skip);
			if (outBeforeCCurveHistogram /*&& histogramCropped*/)
				histNeeded = true;
		}
		if (tcurve) {
			if (!tcurve->isIdentity())
				customToneCurve1.Set(tcurve);
			delete tcurve;
			tcurve = NULL;
		}
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		// create curve bw
		// curve 2
/*		DiagonalCurve* tbwcurve = NULL;
		customToneCurvebw2.Reset();

		if (!curvePointsbw2.empty() && curvePointsbw2[0]>DCT_Linear && curvePointsbw2[0]<DCT_Unchanged) {
			tbwcurve = new DiagonalCurve (curvePointsbw2, CURVES_MIN_POLY_POINTS/skip);
			if (outBeforeCCurveHistogram )
				histNeeded = true;
		}
		if (tbwcurve) {
			if (tbwcurve->isIdentity()) {
				delete tbwcurve;
				tbwcurve = NULL;
			}
			else
				customToneCurvebw2.Set(tbwcurve);
			delete tbwcurve;
			tbwcurve = NULL;
		}
		
		customToneCurvebw1.Reset();

		if (!curvePointsbw.empty() && curvePointsbw[0]>DCT_Linear && curvePointsbw[0]<DCT_Unchanged) {
			tbwcurve = new DiagonalCurve (curvePointsbw, CURVES_MIN_POLY_POINTS/skip);
			if (outBeforeCCurveHistogram )
				histNeeded = true;
		}
		if (tbwcurve) {
			if (tbwcurve->isIdentity()) {
				delete tbwcurve;
				tbwcurve = NULL;
			}
			else if (curveModeb != procparams::BlackWhiteParams::TC_MODE_STD_BW) {
				customToneCurvebw1.Set(tbwcurve);
				delete tbwcurve;
				tbwcurve = NULL;
			}
		}
	*/	
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		for (int i=0; i<=0xffff; i++) {
			float val = dcurve[i];;

			if (histNeeded) {
				float fi=i;
				float hval = hlCurve[i]*fi;
				hval = dcurve[shCurve[hval]*hval];
				//if (needigamma)
				//	hval = igamma2 (hval);
				int hi = (int)(255.0*(hval));
				outBeforeCCurveHistogram[hi] += histogram/*Cropped*/[i] ;
			}

			// if inverse gamma is needed, do it (standard sRGB inverse gamma is applied)
			if (needigamma)
				val = igamma (val, gamma_, start, slope, mul, add);

			outCurve[i] = (65535.0 * val);
		}

		if (tcurve) delete tcurve;

		/*if (outBeforeCCurveHistogram) {
			for (int i=0; i<256; i++) printf("i= %d bchist= %d \n",i,outBeforeCCurveHistogram[i]);
		}*/
	}

	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	void CurveFactory::complexLCurve (double br, double contr, const std::vector<double>& curvePoints,
									 LUTu & histogram, LUTu & histogramCropped, LUTf & outCurve,
									 LUTu & outBeforeCCurveHistogram, int skip, bool & utili) {
		
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
			utili=true;

			std::vector<double> brightcurvePoints;
			brightcurvePoints.resize(9);
			brightcurvePoints.at(0) = double(DCT_NURBS);

			brightcurvePoints.at(1) = 0.; // black point.  Value in [0 ; 1] range
			brightcurvePoints.at(2) = 0.; // black point.  Value in [0 ; 1] range
			
			if (br>0) {
				brightcurvePoints.at(3) = 0.1; // toe point
				brightcurvePoints.at(4) = 0.1+br/150.0; //value at toe point

				brightcurvePoints.at(5) = 0.7; // shoulder point
				brightcurvePoints.at(6) = min(1.0,0.7+br/300.0); //value at shoulder point
			} else {
				brightcurvePoints.at(3) = 0.1-br/150.0; // toe point
				brightcurvePoints.at(4) = 0.1; // value at toe point

				brightcurvePoints.at(5) = min(1.0,0.7-br/300.0); // shoulder point
				brightcurvePoints.at(6) = 0.7; // value at shoulder point
			}
			brightcurvePoints.at(7) = 1.; // white point
			brightcurvePoints.at(8) = 1.; // value at white point
			
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
			utili=true;

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
			//		printf("avg=%f\n",avg);

			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			std::vector<double> contrastcurvePoints;
			contrastcurvePoints.resize(9);
			contrastcurvePoints.at(0) = double(DCT_NURBS);
			
			contrastcurvePoints.at(1) = 0.; // black point.  Value in [0 ; 1] range
			contrastcurvePoints.at(2) = 0.; // black point.  Value in [0 ; 1] range
			
			contrastcurvePoints.at(3) = avg-avg*(0.6-contr/250.0); // toe point
			contrastcurvePoints.at(4) = avg-avg*(0.6+contr/250.0); // value at toe point
			
			contrastcurvePoints.at(5) = avg+(1-avg)*(0.6-contr/250.0); // shoulder point
			contrastcurvePoints.at(6) = avg+(1-avg)*(0.6+contr/250.0); // value at shoulder point
			
			contrastcurvePoints.at(7) = 1.; // white point
			contrastcurvePoints.at(8) = 1.; // value at white point
			
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
			utili=true;//if active

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
			if (!outCurve)
				outCurve(65536, 0);
			for (int i=0; i<65536; i++) {
				// apply custom/parametric/NURBS curve, if any
				float val = tcurve->getVal ((float)i/65536.0f);
				outCurve[i] = (65536.0f * val);
			}
			delete tcurve;
		}
		// let the LUTf empty for identity curves
		else {
			outCurve.reset();
		}

	}
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void ColorAppearance::Reset() {
    lutColCurve.reset();
}

// Fill a LUT with X/Y, ranged 0xffff
void ColorAppearance::Set(Curve *pCurve) {
    lutColCurve(65536);
    for (int i=0; i<65536; i++) lutColCurve[i] = pCurve->getVal(double(i)/65535.) * 65535.;
}

void ToneCurve::Reset() {
    lutToneCurve.reset();
}

// Fill a LUT with X/Y, ranged 0xffff
void ToneCurve::Set(Curve *pCurve) {
    lutToneCurve(65536);
    for (int i=0; i<65536; i++) lutToneCurve[i] = pCurve->getVal(double(i)/65535.) * 65535.;
}

}
