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
#include "opthelper.h"
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

	/** @ brief Return the number of control points of the curve
	 * This method return the number of control points of a curve. Not suitable for parametric curves.
	 * @return number of control points of the curve. 0 will be sent back for Parametric curves
	 */
	int Curve::getSize () const {
		return N;
	}

	/** @ brief Return the a control point's value
	 * This method return a control points' value. Not suitable for parametric curves.
	 * @param cpNum id of the control points we're interested in
	 * @param x Y value of the control points, or -1 if invalid
	 * @param y Y value of the control points, or -1 if invalid
	 */
	void Curve::getControlPoint(int cpNum, double &x, double &y) const {
		if (this->x && cpNum < N) {
			x = this->x[cpNum];
			y = this->y[cpNum];
		}
		else {
			x = y = -1.;
		}
	}

    // Wikipedia sRGB: Unlike most other RGB color spaces, the sRGB gamma cannot be expressed as a single numerical value.
    // The overall gamma is approximately 2.2, consisting of a linear (gamma 1.0) section near black, and a non-linear section elsewhere involving a 2.4 exponent 
    // and a gamma (slope of log output versus log input) changing from 1.0 through about 2.3.
    const double CurveFactory::sRGBGamma = 2.2;
    const double CurveFactory::sRGBGammaCurve = 2.4;

	SSEFUNCTION void fillCurveArray(DiagonalCurve* diagCurve, LUTf &outCurve, int skip, bool needed) {

		if (needed) {
			LUTf lutCurve (65536);

			for (int i=0; i<=0xffff; i+= i<0xffff-skip ? skip : 1 ) {
				// change to [0,1] range
				float val = (float)i / 65535.f;
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
				outCurve[i] = (65535.f * lutCurve[i]);
			}
		}
		else {
#ifdef __SSE2__
			__m128 fourv = _mm_set1_ps(4.f);
			__m128 iv = _mm_set_ps(3.f,2.f,1.f,0.f);
			for (int i=0; i<=0xfffc; i+=4) {
				_mm_storeu_ps(&outCurve[i],iv);
				iv += fourv;
			}
#else
			for (int i=0; i<=0xffff; i++) {
				outCurve[i] = (float)i;
			}
#endif
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
	for (int i=0; i<48000; i++) {//32768*1.414  + ...
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
	if (histNeeded) {
		for (int i=0; i<32768; i++) {
			double hval = CLIPD((double)i / 32767.0);
			int hi = (int)(255.0*hval);
			outBeforeCCurveHistogram[hi] += histogram[i] ;
		}
	}
	if (histNeededC) {
		for (int i=0; i<48000; i++) {//32768*1.414  + ...
			double hval = CLIPD((double)i / 47999.0);
			int hi = (int)(255.0*hval); //
			outBeforeCCurveHistogramC[hi] += histogramC[i] ;
		}
	}
	
	if (tcurve) delete tcurve;

}
// add curve Denoise : C=f(C)
void CurveFactory::denoiseCC ( bool & ccdenoiseutili,const std::vector<double>& cccurvePoints, LUTf & NoiseCCcurve,int skip){
		bool needed;
		DiagonalCurve* dCurve = NULL;
		LUTf dCcurve(65536,0);
	
		float val;
		for (int i=0; i<48000; i++) {  
				dCcurve[i] = (float)i / 47999.0;
		}
		
		needed = false;
		if (!cccurvePoints.empty() && cccurvePoints[0]!=0) {
			dCurve = new DiagonalCurve (cccurvePoints, CURVES_MIN_POLY_POINTS/skip);
			
			if (dCurve && !dCurve->isIdentity())
				{needed = true;ccdenoiseutili=true;}
		}
		fillCurveArray(dCurve, NoiseCCcurve, skip, needed);
		//NoiseCCcurve.dump("Noise");
		
		if (dCurve) {
			delete dCurve;
			dCurve = NULL;
		}
}


void CurveFactory::curveBW (
		const std::vector<double>& curvePointsbw, const std::vector<double>& curvePointsbw2,
		LUTu & histogrambw, LUTu & outBeforeCCurveHistogrambw,//for Luminance
		ToneCurve & customToneCurvebw1, ToneCurve & customToneCurvebw2, int skip)
{

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
	if (histNeeded) {
		LUTf dcurve(65536,0);
		
		float val;
		for (int i=0; i<32768; i++) {
				val = (float)i / 32767.f;
				dcurve[i] = CLIPD(val);
			}

		for (int i=0; i<32768; i++) {
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
		if(histNeededCL)
			for (int i=0; i<50000; i++) {//32768*1.414  + ...
				int hi = (int)(255.0*CLIPD((float)i / 49999.0)); //
				outBeforeCLurveHistogram[hi] += histogramcl[i] ;
			}

		
		fillCurveArray(dCurve, clCurve, skip, needed);
		if (dCurve) {
			delete dCurve;
			dCurve = NULL;
		}
}

// add curve Colortoning : C=f(L)
void CurveFactory::curveToningCL ( bool & clctoningutili,const std::vector<double>& clcurvePoints, LUTf & clToningCurve,int skip){
		bool needed;
		DiagonalCurve* dCurve = NULL;
		
		needed = false;
		if (!clcurvePoints.empty() && clcurvePoints[0]!=0) {
			dCurve = new DiagonalCurve (clcurvePoints, CURVES_MIN_POLY_POINTS/skip);
			
			if (dCurve && !dCurve->isIdentity())
				{needed = true;clctoningutili=true;}
		}
		fillCurveArray(dCurve, clToningCurve, skip, needed);
	//	clToningCurve.dump("CLToning");
		if (dCurve) {
			delete dCurve;
			dCurve = NULL;
		}
}

// add curve Colortoning : CLf(L)
void CurveFactory::curveToningLL ( bool & llctoningutili,const std::vector<double>& llcurvePoints, LUTf & llToningCurve, int skip){
		bool needed;
		DiagonalCurve* dCurve = NULL;

		needed = false;
		if (!llcurvePoints.empty() && llcurvePoints[0]!=0) {
			dCurve = new DiagonalCurve (llcurvePoints, CURVES_MIN_POLY_POINTS/skip);
			
			if (dCurve && !dCurve->isIdentity())
				{needed = true;llctoningutili=true;}
		}
		fillCurveArray(dCurve, llToningCurve, skip, needed);
//		llToningCurve.dump("LLToning");

		if (dCurve) {
			delete dCurve;
			dCurve = NULL;
		}
}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	void CurveFactory::complexsgnCurve (float adjustr, bool & autili,  bool & butili, bool & ccutili, bool & cclutili, double saturation, double rstprotection,
										const std::vector<double>& acurvePoints, const std::vector<double>& bcurvePoints,const std::vector<double>& cccurvePoints,
										const std::vector<double>& lccurvePoints, LUTf & aoutCurve, LUTf & boutCurve, LUTf & satCurve, LUTf & lhskCurve,
										LUTu & histogramC, LUTu & histogramLC, LUTu & outBeforeCCurveHistogram,LUTu & outBeforeLCurveHistogram, //for chroma
										int skip) {
		
		
		//-----------------------------------------------------

		bool needed;
		DiagonalCurve* dCurve = NULL;
		LUTf dCcurve(65536,0);
		int k=48000;//32768*1.41
		if(outBeforeCCurveHistogram || outBeforeLCurveHistogram) {
			for (int i=0; i<k*adjustr; i++) {  //# 32768*1.414  approximation maxi for chroma
					dCcurve[i] = (float)i / (k*adjustr-1);
			}
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
		if (histNeededC) {
			for (int i=0; i<k*adjustr; i++) {//32768*1.414  + ...
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
		if (histNeededLC) {
			for (int i=0; i<k*adjustr; i++) {//32768*1.414  + ...
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

	SSEFUNCTION void CurveFactory::complexCurve (double ecomp, double black, double hlcompr, double hlcomprthresh,
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
		
		if (comp<=0.0f)
			for (int i=0; i<0x10000; i++)
				hlCurve[i]=exp_scale;
		else {
			for (int i=0; i<=shoulder; i++)
				hlCurve[i]=exp_scale;
			float scalemshoulder = scale - shoulder;

			// SSE makes more sense than omp here
#ifdef __SSE2__
			int i;
			__m128 exp_scalev = _mm_set1_ps(exp_scale);
			__m128 scalemshoulderv = _mm_set1_ps(scalemshoulder);
			__m128 compv = _mm_set1_ps(comp);
			__m128 valv = _mm_set_ps(4.f,3.f,2.f,1.f);
			__m128 onev = _mm_set1_ps(1.f);
			__m128 fourv = _mm_set1_ps(4.f);
			for (i=shoulder+1; i<0xFFFD; i+=4) {
				// change to [0,1] range
				__m128 Rv = valv*compv/(scalemshoulderv);
				_mm_storeu_ps(&hlCurve[i],xlogf(onev+Rv*exp_scalev)/Rv);
				valv += fourv;
			}
			for (; i<0x10000; i++) {
				// change to [0,1] range
				float val = (float)i-shoulder;
				float R = val*comp/(scalemshoulder);
				hlCurve[i] = xlogf(1.f+R*exp_scale)/R;
			}

#else
			for (int i=shoulder+1; i<0x10000; i++) {
				// change to [0,1] range
				float val = (float)i-shoulder;
				float R = val*comp/(scalemshoulder);
				hlCurve[i] = xlogf(1.f+R*exp_scale)/R;
			}
#endif

		}


		// curve without contrast
		LUTf dcurve(0x10000);

		//%%%%%%%%%%%%%%%%%%%%%%%%%%
		// change to [0,1] range

		float val = 1.f/65535.f;
		float	val2 = simplebasecurve (val, black, 0.015*shcompr);
		shCurve[0] = CLIPD(val2)/val;
		val = 0.0;
				// gamma correction
		if (gamma_>1.)
			val = gamma (val, gamma_, start, slope, mul, add);
		
		// apply brightness curve
		if (brightcurve)
			val = brightcurve->getVal (val);  // TODO: getVal(double) is very slow! Optimize with a LUTf

		// store result in a temporary array
		dcurve[0] = CLIPD(val);

#pragma omp parallel for
		for (int i=1; i<0x10000; i++) {
			float val;
				val = (float)i / 65535.0f;
			
			float	val2 = simplebasecurve (val, black, 0.015*shcompr);
			shCurve[i] = CLIPD(val2)/val;

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
			unsigned int sum = 0;
			float avg = 0; 
			//double sqavg = 0;
			for (int i=0; i<=0xffff; i++) {
				float fi=i;
				fi *= hlCurve[i];
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
			float val = dcurve[i];

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

			outCurve[i] = (65535.f * val);
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

			DiagonalCurve* contrastcurve = NULL;

			// compute mean luminance of the image with the curve applied
			int sum = 0;
			float avg = 0; 
			//float sqavg = 0;
			for (int i=0; i<32768; i++) {
				avg += dcurve[i] * histogram[i];
				//sqavg += dcurve[i]*dcurve[i] * histogram[i];
				sum += histogram[i];
			}
			if(sum) {
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
				
				contrastcurve = new DiagonalCurve (contrastcurvePoints, CURVES_MIN_POLY_POINTS/skip);
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			} else {
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				// sum has an invalid value (next to 0, producing a division by zero, so we create a fake contrast curve, producing a white image
				std::vector<double> contrastcurvePoints;
				contrastcurvePoints.resize(5);
				contrastcurvePoints.at(0) = double(DCT_NURBS);

				contrastcurvePoints.at(1) = 0.; // black point.  Value in [0 ; 1] range
				contrastcurvePoints.at(2) = 1.; // black point.  Value in [0 ; 1] range

				contrastcurvePoints.at(3) = 1.; // white point
				contrastcurvePoints.at(4) = 1.; // value at white point

				contrastcurve = new DiagonalCurve (contrastcurvePoints, CURVES_MIN_POLY_POINTS/skip);
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			}
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

void OpacityCurve::Reset() {
	lutOpacityCurve.reset();
}

void OpacityCurve::Set(const Curve *pCurve) {
	if (pCurve->isIdentity()) {
		lutOpacityCurve.reset(); // raise this value if the quality suffers from this number of samples
		return;
	}
	lutOpacityCurve(501); // raise this value if the quality suffers from this number of samples

	for (int i=0; i<501; i++) lutOpacityCurve[i] = pCurve->getVal(double(i)/500.);
	//lutOpacityCurve.dump("opacity");
}

void OpacityCurve::Set(const std::vector<double> &curvePoints, bool &opautili) {
	FlatCurve* tcurve = NULL;

	if (!curvePoints.empty() && curvePoints[0]>FCT_Linear && curvePoints[0]<FCT_Unchanged) {
		tcurve = new FlatCurve (curvePoints, false, CURVES_MIN_POLY_POINTS/2);
		tcurve->setIdentityValue(0.);
	}
	if (tcurve) {
		Set(tcurve);opautili=true;
		delete tcurve;
		tcurve = NULL;
	}
}

NoiseCurve::NoiseCurve() : sum(0.f) {};

void NoiseCurve::Reset() {
	lutNoiseCurve.reset();
	sum = 0.f;
}

void NoiseCurve::Set(const Curve &pCurve) {
	if (pCurve.isIdentity()) {
		Reset(); // raise this value if the quality suffers from this number of samples
		return;
	}
	lutNoiseCurve(501); // raise this value if the quality suffers from this number of samples
	sum=0.f;
	for (int i=0; i<501; i++) {
		lutNoiseCurve[i] = pCurve.getVal(double(i)/500.); 
		if(lutNoiseCurve[i] < 0.01f)
			lutNoiseCurve[i] = 0.01f;//avoid 0.f for wavelet : under 0.01f quasi no action for each value
		sum += lutNoiseCurve[i]; //minima for Wavelet about 6.f or 7.f quasi no action
	}
	//lutNoisCurve.dump("Nois");
}

void NoiseCurve::Set(const std::vector<double> &curvePoints) {

	if (!curvePoints.empty() && curvePoints[0]>FCT_Linear && curvePoints[0]<FCT_Unchanged) {
		FlatCurve tcurve(curvePoints, false, CURVES_MIN_POLY_POINTS/2);
		tcurve.setIdentityValue(0.);
		Set(tcurve);
	} else {
		Reset();
	}
}


void ColorGradientCurve::Reset() {
	lut1.reset();
	lut2.reset();
	lut3.reset();
}

void ColorGradientCurve::SetXYZ(const Curve *pCurve, const double xyz_rgb[3][3], const double rgb_xyz[3][3], float satur, float lumin) {
	if (pCurve->isIdentity()) {
		lut1.reset();
		lut2.reset();
		lut3.reset();
		return;
	}
	if (!lut1) {
		lut1(501);
		lut2(501);
		lut3(501);
	}

	float r, g, b, xx, yy, zz;
	float lr1,lr2;
	int upperBound = lut1.getUpperBound();
	if (pCurve->isIdentity()) {
		Color::hsv2rgb(0.5f, satur, lumin, r, g, b);
		Color::rgbxyz(r, g, b, xx, yy, zz, xyz_rgb);
		
		for (int i=0; i<=500; ++i) {
			// WARNING: set the identity value according to what is set in the GUI
			lut1[i] = xx;
			lut2[i] = yy;
			lut3[i] = zz;
		}
		return;
	}

	int nPoints = pCurve->getSize();
	int ptNum = 0;
	double nextX, nextY;
	pCurve->getControlPoint(ptNum, nextX, nextY);
	double prevY = nextY;
	double dY = 0.;
	low=nextX;
	lr1=(0.5f+low)/2.f;//optimize use of gamut in low light..one can optimize more using directly low ?
	//lr1=low;
	for (int i=0; i<=upperBound; ++i) {
		double x = double(i)/double(upperBound);

		if (x > nextX) {
			++ptNum;
			if (ptNum < nPoints) {
				prevY = nextY;
				pCurve->getControlPoint(ptNum, nextX, nextY);
				dY = nextY - prevY;
				high=nextX;
				lr2=(0.5f + high)/2.f;//optimize use of gamut in high light..one can optimize more using directly high ?
				//lr2=high;
			}
		}

		if (!ptNum) {
			Color::hsv2rgb(float(prevY), satur, lr1, r, g, b);
			Color::rgbxyz(r, g, b, xx, yy, zz, xyz_rgb);
			lut1[i] = xx;
			lut2[i] = yy;
			lut3[i] = zz;
		}
		else if (ptNum >= nPoints) {
			Color::hsv2rgb(float(nextY), satur, lr2, r, g, b);
			Color::rgbxyz(r, g, b, xx, yy, zz, xyz_rgb);
			lut1[i] = xx;
			lut2[i] = yy;
			lut3[i] = zz;
		}
		else {
			double currY = pCurve->getVal(x) - prevY;
			if (dY > 0.000001 || dY < -0.000001) {
				float r1, g1, b1, r2, g2, b2, ro, go, bo;
				Color::hsv2rgb(float(prevY), satur, lr1, r1, g1, b1);
				Color::hsv2rgb(float(nextY), satur, lr2, r2, g2, b2);
				bool chr = false;
				bool lum = true;
				LUTf dum;
				float X1,X2,Y1,Y2,Z1,Z2,L1,a_1,b_1,c1,h1;
				Color::rgbxyz(r2, g2, b2, X2, Y2, Z2, xyz_rgb);
				Color::rgbxyz(r1, g1, b1, X1, Y1, Z1, xyz_rgb);
				//I use XYZ to mix color 1 and 2 rather than rgb (gamut) and rather than Lab artifacts
				X1 = X1 + (X2-X1)*currY/dY; if(X1<0.f) X1=0.f;//negative value not good
				Y1 = Y1 + (Y2-Y1)*currY/dY; if(Y1<0.f) Y1=0.f;
				Z1 = Z1 + (Z2-Z1)*currY/dY; if(Z1<0.f) Z1=0.f;
				Color::XYZ2Lab(X1, Y1, Z1, L1, a_1, b_1);//prepare to gamut control
				Color::Lab2Lch(a_1, b_1, c1, h1);
				float Lr=L1/327.68f;
				float RR,GG,BB;
				#ifndef NDEBUG
					bool neg=false;
					bool more_rgb=false;
					//gamut control : Lab values are in gamut
					Color::gamutLchonly(h1,Lr,c1, RR, GG, BB, xyz_rgb, false, 0.15f, 0.96f, neg, more_rgb);
				#else
					Color::gamutLchonly(h1,Lr,c1, RR, GG, BB, xyz_rgb, false, 0.15f, 0.96f);
				#endif
				L1=Lr*327.68f;
				float a,b,X,Y,Z;
				// converting back to rgb
				Color::Lch2Lab(c1, h1, a, b);
				Color::Lab2XYZ(L1, a, b, X, Y, Z);
				lut1[i] = X;
				lut2[i] = Y;
				lut3[i] = Z;
			}
			else {
				Color::hsv2rgb(float(nextY), satur, lumin, r, g, b);
				Color::rgbxyz(r, g, b, xx, yy, zz, xyz_rgb);
				lut1[i] = xx;
				lut2[i] = yy;
				lut3[i] = zz;
			}
		}
	}
	/*
	#ifndef NDEBUG
	lutRed.dump("red");
	lutGreen.dump("green");
	lutBlue.dump("blue");
	#endif
	*/
}

void ColorGradientCurve::SetXYZ(const std::vector<double> &curvePoints, const double xyz_rgb[3][3], const double rgb_xyz[3][3], float satur, float lumin) {
	FlatCurve* tcurve = NULL;

	if (!curvePoints.empty() && curvePoints[0]>FCT_Linear && curvePoints[0]<FCT_Unchanged) {
		tcurve = new FlatCurve (curvePoints, false, CURVES_MIN_POLY_POINTS/2);
	}
	if (tcurve) {
		SetXYZ(tcurve, xyz_rgb, rgb_xyz, satur, lumin);
		delete tcurve;
		tcurve = NULL;
	}
}

void ColorGradientCurve::SetRGB(const Curve *pCurve, const double xyz_rgb[3][3], const double rgb_xyz[3][3]) {
	if (pCurve->isIdentity()) {
		lut1.reset();
		lut2.reset();
		lut3.reset();
		return;
	}
	if (!lut1) {
		lut1(501);
		lut2(501);
		lut3(501);
	}

	float r, g, b;

	int upperBound = lut1.getUpperBound();

	int nPoints = pCurve->getSize();
	int ptNum = 0;
	double nextX, nextY;
	pCurve->getControlPoint(ptNum, nextX, nextY);
	double prevY = nextY;
	double dY = 0.;
	Color::eInterpolationDirection dir = Color::ID_DOWN;
	for (int i=0; i<=upperBound; ++i) {
		double x = double(i)/double(upperBound);

		if (x > nextX) {
			++ptNum;
			if (ptNum < nPoints) {
				prevY = nextY;
				pCurve->getControlPoint(ptNum, nextX, nextY);
				dY = nextY - prevY;
				dir = Color::getHueInterpolationDirection(prevY, nextY, Color::IP_SHORTEST);
			}
		}

		if (!ptNum) {
			Color::hsv2rgb(float(prevY), 1.f, 1.f, r, g, b);
			lut1[i] = r;
			lut2[i] = g;
			lut3[i] = b;
		}
		else if (ptNum >= nPoints) {
			Color::hsv2rgb(float(nextY), 1.f, 1.f, r, g, b);
			lut1[i] = r;
			lut2[i] = g;
			lut3[i] = b;
		}
		else {
			double currY = pCurve->getVal(x) - prevY;
			if (dY > 0.0000001 || dY < -0.0000001) {
				#if 1
				float ro, go, bo;
				double h2 = Color::interpolateHueHSV(prevY, nextY, currY/dY, dir);
				Color::hsv2rgb(h2, 1.f, 1.f, ro, go, bo);
				#else
				float r1, g1, b1, r2, g2, b2, ro, go, bo;
				Color::hsv2rgb(float(prevY), 1., 1., r1, g1, b1);
				Color::hsv2rgb(float(nextY), 1., 1., r2, g2, b2);
				Color::interpolateRGBColor(currY/dY, r1, g1, b1, r2, g2, b2, Color::CHANNEL_LIGHTNESS|Color::CHANNEL_CHROMATICITY|Color::CHANNEL_HUE, xyz_rgb, rgb_xyz, ro, go, bo);
				#endif
				lut1[i] = ro;
				lut2[i] = go;
				lut3[i] = bo;
			}
			else {
				Color::hsv2rgb(float(nextY), 1.f, 1.f, r, g, b);
				lut1[i] = r;
				lut2[i] = g;
				lut3[i] = b;
			}
		}
	}
	/*
	#ifndef NDEBUG
	lut1.dump("red");
	lut2.dump("green");
	lut3.dump("blue");
	#endif
	*/
}

void ColorGradientCurve::SetRGB(const std::vector<double> &curvePoints, const double xyz_rgb[3][3], const double rgb_xyz[3][3]) {
	FlatCurve* tcurve = NULL;

	if (!curvePoints.empty() && curvePoints[0]>FCT_Linear && curvePoints[0]<FCT_Unchanged) {
		tcurve = new FlatCurve (curvePoints, false, CURVES_MIN_POLY_POINTS/2);
	}
	if (tcurve) {
		SetRGB(tcurve, xyz_rgb, rgb_xyz);
		delete tcurve;
		tcurve = NULL;
	}
}

void ColorGradientCurve::getVal(float index, float &r, float &g, float &b) const {
	r = lut1[index*500.f];
	g = lut2[index*500.f];
	b = lut3[index*500.f];
}


}
