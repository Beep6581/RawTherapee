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
#ifdef _OPENMP
#include <omp.h>
#endif

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
	const float gamma_ = Color::sRGBGammaCurve;
	const float start = expf(gamma_*logf( -0.055 / ((1.0/gamma_-1.0)*1.055 )));
	const float slope = 1.055 * powf (start, 1.0/gamma_-1) - 0.055/start;
	const float mul = 1.055;
	const float add = 0.055;

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
			customToneCurvebw2.Set(tcurve, gamma_, start, slope, mul, add);
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
			customToneCurvebw1.Set(tcurve, gamma_, start, slope, mul, add);
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
// add curve Lab wavelet : Cont=f(L)
void CurveFactory::curveWavContL ( bool & wavcontlutili,const std::vector<double>& wavclcurvePoints, LUTf & wavclCurve, /*LUTu & histogramwavcl, LUTu & outBeforeWavCLurveHistogram,*/int skip){
		bool needed;
		DiagonalCurve* dCurve = NULL;
		
	//	if (outBeforeWavCLurveHistogram)		
	//		outBeforeWavCLurveHistogram.clear();
		bool histNeededCL = false;

		needed = false;
		if (!wavclcurvePoints.empty() && wavclcurvePoints[0]!=0) {
			dCurve = new DiagonalCurve (wavclcurvePoints, CURVES_MIN_POLY_POINTS/skip);
		//	if (outBeforeWavCLurveHistogram)
		//		histNeededCL = true;
			
			if (dCurve && !dCurve->isIdentity())
				{needed = true;wavcontlutili=true;}
		}
	//	if(histNeededCL)
			for (int i=0; i<32768; i++) {//32768*1.414  + ...
				int hi = (int)(255.0*CLIPD((float)i / 32767.0)); //
			//	outBeforeWavCLurveHistogram[hi] += histogramwavcl[i] ;
			}

		
		fillCurveArray(dCurve, wavclCurve, skip, needed);
	//	wavclCurve.dump("WL");

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
									 double shcompr, double br, double contr, 
									 procparams::ToneCurveParams::eTCModeId curveMode, const std::vector<double>& curvePoints,
									 procparams::ToneCurveParams::eTCModeId curveMode2, const std::vector<double>& curvePoints2,
									 LUTu & histogram, LUTu & histogramCropped,
									 LUTf & hlCurve, LUTf & shCurve, LUTf & outCurve,
									 LUTu & outBeforeCCurveHistogram,
									 ToneCurve & customToneCurve1,
									 ToneCurve & customToneCurve2,
									 int skip) {

		// the curve shapes are defined in sRGB gamma, but the output curves will operate on linear floating point data,
		// hence we do both forward and inverse gamma conversions here.
		const float gamma_ = Color::sRGBGammaCurve;
		const float start = expf(gamma_*logf( -0.055 / ((1.0/gamma_-1.0)*1.055 )));
		const float slope = 1.055 * powf (start, 1.0/gamma_-1) - 0.055/start;
		const float mul = 1.055;
		const float add = 0.055;
		
		// a: slope of the curve, black: starting point at the x axis
		const float a = powf (2.0, ecomp);
		
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

		hlCurve.setClip(LUT_CLIP_BELOW); // used LUT_CLIP_BELOW, because we want to have a baseline of 2^expcomp in this curve. If we don't clip the lut we get wrong values, see Issue 2621 #14 for details
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

			for (int i=shoulder+1; i<0x10000; i++) {
				// change to [0,1] range
				float val = (float)i-shoulder;
				float R = val*comp/(scalemshoulder);
				hlCurve[i] = xlog(1.0+R*exp_scale)/R; // don't use xlogf or 1.f here. Leads to errors caused by too low precision
			}
		}


		// curve without contrast
		LUTf dcurve(0x10000);

		//%%%%%%%%%%%%%%%%%%%%%%%%%%
		// change to [0,1] range
		shCurve.setClip(LUT_CLIP_ABOVE); // used LUT_CLIP_ABOVE, because the curve converges to 1.0 at the upper end and we don't want to exceed this value.
		float val = 1.f/65535.f;
		float	val2 = simplebasecurve (val, black, 0.015*shcompr);
		shCurve[0] = CLIPD(val2)/val;
		val = 0.0;
		// gamma correction
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
				customToneCurve2.Set(tcurve, gamma_, start, slope, mul, add);
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
				customToneCurve1.Set(tcurve, gamma_, start, slope, mul, add);
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

				// RGB curves are defined with sRGB gamma, but operate on linear data
				float val = float(i)/65535.f;
				val = CurveFactory::gamma2 (val);
				val = tcurve->getVal(val);
				val = CurveFactory::igamma2 (val);
				//float val = tcurve->getVal ((float)i/65536.0f);
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
void ToneCurve::Set(Curve *pCurve, float gamma, float start, float slope, float mul, float add) {
	lutToneCurve(65536);
	if (gamma <= 0.0 || gamma == 1.) {
		for (int i=0; i<65536; i++) lutToneCurve[i] = (float)pCurve->getVal(float(i)/65535.f) * 65535.f;
	} else {
		// apply gamma, that is 'pCurve' is defined with the given gamma and here we convert it to a curve in linear space
		for (int i=0; i<65536; i++) {
			float val = float(i)/65535.f;
			val = CurveFactory::gamma (val, gamma, start, slope, mul, add);
			val = pCurve->getVal(val);
			val = CurveFactory::igamma (val, gamma, start, slope, mul, add);
			lutToneCurve[i] = val * 65535.f;
		}
	}
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


WavCurve::WavCurve() : sum(0.f) {};

void WavCurve::Reset() {
	lutWavCurve.reset();
	sum = 0.f;
}

void WavCurve::Set(const Curve &pCurve) {
	if (pCurve.isIdentity()) {
		Reset(); // raise this value if the quality suffers from this number of samples
		return;
	}
	lutWavCurve(501); // raise this value if the quality suffers from this number of samples
	sum=0.f;
	for (int i=0; i<501; i++) {
		lutWavCurve[i] = pCurve.getVal(double(i)/500.);
		if(lutWavCurve[i] < 0.02f)
			lutWavCurve[i] = 0.02f;//avoid 0.f for wavelet : under 0.01f quasi no action for each value
		sum += lutWavCurve[i];
	}
	//lutWavCurve.dump("wav");
}
void WavCurve::Set(const std::vector<double> &curvePoints) {
	
	if (!curvePoints.empty() && curvePoints[0]>FCT_Linear && curvePoints[0]<FCT_Unchanged) {
		FlatCurve tcurve(curvePoints, false, CURVES_MIN_POLY_POINTS/2);
		tcurve.setIdentityValue(0.);
		Set(tcurve);
	} else {
		Reset();
	}
}


WavOpacityCurveRG::WavOpacityCurveRG(){};

void WavOpacityCurveRG::Reset() {
	lutOpacityCurveRG.reset();
}

void WavOpacityCurveRG::Set(const Curve &pCurve) {
	if (pCurve.isIdentity()) {
	    Reset(); // raise this value if the quality suffers from this number of samples
		return;
	}
	lutOpacityCurveRG(501); // raise this value if the quality suffers from this number of samples

	for (int i=0; i<501; i++) lutOpacityCurveRG[i] = pCurve.getVal(double(i)/500.);
}

void WavOpacityCurveRG::Set(const std::vector<double> &curvePoints) {
	if (!curvePoints.empty() && curvePoints[0]>FCT_Linear && curvePoints[0]<FCT_Unchanged) {
		FlatCurve tcurve(curvePoints, false, CURVES_MIN_POLY_POINTS/2);
		tcurve.setIdentityValue(0.);
		Set(tcurve);
	} else {
		Reset();
	}

}

WavOpacityCurveBY::WavOpacityCurveBY(){};

void WavOpacityCurveBY::Reset() {
	lutOpacityCurveBY.reset();
}

void WavOpacityCurveBY::Set(const Curve &pCurve) {
	if (pCurve.isIdentity()) {
		lutOpacityCurveBY.reset(); // raise this value if the quality suffers from this number of samples
		return;
	}
	lutOpacityCurveBY(501); // raise this value if the quality suffers from this number of samples

	for (int i=0; i<501; i++) lutOpacityCurveBY[i] = pCurve.getVal(double(i)/500.);
}

void WavOpacityCurveBY::Set(const std::vector<double> &curvePoints) {
	if (!curvePoints.empty() && curvePoints[0]>FCT_Linear && curvePoints[0]<FCT_Unchanged) {
		FlatCurve tcurve(curvePoints, false, CURVES_MIN_POLY_POINTS/2);
		tcurve.setIdentityValue(0.);
		Set(tcurve);
	} else {
		Reset();
	}
}

WavOpacityCurveW::WavOpacityCurveW(){};

void WavOpacityCurveW::Reset() {
	lutOpacityCurveW.reset();
}

void WavOpacityCurveW::Set(const Curve &pCurve) {
	if (pCurve.isIdentity()) {
		lutOpacityCurveW.reset(); // raise this value if the quality suffers from this number of samples
		return;
	}
	lutOpacityCurveW(501); // raise this value if the quality suffers from this number of samples

	for (int i=0; i<501; i++) lutOpacityCurveW[i] = pCurve.getVal(double(i)/500.);
}

void WavOpacityCurveW::Set(const std::vector<double> &curvePoints) {
	if (!curvePoints.empty() && curvePoints[0]>FCT_Linear && curvePoints[0]<FCT_Unchanged) {
		FlatCurve tcurve(curvePoints, false, CURVES_MIN_POLY_POINTS/2);
		tcurve.setIdentityValue(0.);
		Set(tcurve);
	} else {
		Reset();
	}
}

WavOpacityCurveWL::WavOpacityCurveWL(){};

void WavOpacityCurveWL::Reset() {
	lutOpacityCurveWL.reset();
}

void WavOpacityCurveWL::Set(const Curve &pCurve) {
	if (pCurve.isIdentity()) {
		lutOpacityCurveWL.reset(); // raise this value if the quality suffers from this number of samples
		return;
	}
	lutOpacityCurveWL(501); // raise this value if the quality suffers from this number of samples

	for (int i=0; i<501; i++) lutOpacityCurveWL[i] = pCurve.getVal(double(i)/500.);
}

void WavOpacityCurveWL::Set(const std::vector<double> &curvePoints) {
	if (!curvePoints.empty() && curvePoints[0]>FCT_Linear && curvePoints[0]<FCT_Unchanged) {
		FlatCurve tcurve(curvePoints, false, CURVES_MIN_POLY_POINTS/2);
		tcurve.setIdentityValue(0.);
		Set(tcurve);
	} else {
		Reset();
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

// this is a generic cubic spline implementation, to clean up we could probably use something already existing elsewhere
void PerceptualToneCurve::cubic_spline(const float x[], const float y[], const int len, const float out_x[], float out_y[], const int out_len) {
    int i, j;

    float **A = (float **)malloc(2 * len * sizeof(*A));
    float *As = (float *)calloc(1, 2 * len * 2 * len * sizeof(*As));
    float *b = (float *)calloc(1, 2*len*sizeof(*b));
    float *c = (float *)calloc(1, 2*len*sizeof(*c));
    float *d = (float *)calloc(1, 2*len*sizeof(*d));

    for (i = 0; i < 2*len; i++) {
        A[i] = &As[2*len*i];
    }

    for (i = len-1; i > 0; i--) {
        b[i] = (y[i] - y[i-1]) / (x[i] - x[i-1]);
        d[i-1] = x[i] - x[i-1];
    }

    for (i = 1; i < len-1; i++) {
        A[i][i] = 2 * (d[i-1] + d[i]);
        if (i > 1) {
            A[i][i-1] = d[i-1];
            A[i-1][i] = d[i-1];
        }
        A[i][len-1] = 6 * (b[i+1] - b[i]);
    }

    for(i = 1; i < len-2; i++) {
        float v = A[i+1][i] / A[i][i];
        for(j = 1; j <= len-1; j++) {
            A[i+1][j] -= v * A[i][j];
        }
    }

    for(i = len-2; i > 0; i--) {
        float acc = 0;
        for(j = i; j <= len-2; j++) {
            acc += A[i][j]*c[j];
        }
        c[i] = (A[i][len-1] - acc) / A[i][i];
    }

    for (i = 0; i < out_len; i++) {
        float x_out = out_x[i];
        float y_out = 0;
        for (j = 0; j < len-1; j++) {
            if (x[j] <= x_out && x_out <= x[j+1]) {
                float v = x_out - x[j];
                y_out = y[j] +
                    ((y[j+1] - y[j]) / d[j] - (2 * d[j] * c[j] + c[j+1] * d[j]) / 6) * v +
                    (c[j] * 0.5) * v*v +
                    ((c[j+1] - c[j]) / (6 * d[j])) * v*v*v;
            }
        }
        out_y[i] = y_out;
    }
    free(A);
    free(As);
    free(b);
    free(c);
    free(d);
}

// generic function for finding minimum of f(x) in the a-b range using the interval halving method
float PerceptualToneCurve::find_minimum_interval_halving(float (*func)(float x, void *arg), void *arg, float a, float b, float tol, int nmax) {
	float L = b - a;
	float x = (a + b) * 0.5;
	for (int i = 0; i < nmax; i++) {
		float f_x = func(x, arg);
		if ((b - a) * 0.5 < tol) {
			return x;
		}
		float x1 = a + L/4;
		float f_x1 = func(x1, arg);
		if (f_x1 < f_x) {
			b = x;
			x = x1;
		} else {
			float x2 = b - L/4;
			float f_x2 = func(x2, arg);
			if (f_x2 < f_x) {
				a = x;
				x = x2;
			} else {
				a = x1;
				b = x2;
			}
		}
		L = b - a;
	}
	return x;
}

struct find_tc_slope_fun_arg {
	const ToneCurve * tc;
};

float PerceptualToneCurve::find_tc_slope_fun(float k, void *arg)
{
	struct find_tc_slope_fun_arg *a = (struct find_tc_slope_fun_arg *)arg;
	float areasum = 0;
	const int steps = 10;
	for (int i = 0; i < steps; i++) {
		float x = 0.1 + ((float)i / (steps-1)) * 0.5; // testing (sRGB) range [0.1 - 0.6], ie ignore highligths and dark shadows
		float y = CurveFactory::gamma2(a->tc->lutToneCurve[CurveFactory::igamma2(x) * 65535] / 65535.0);
		float y1 = k * x;
		if (y1 > 1) y1 = 1;
		areasum += (y - y1) * (y - y1); // square is a rough approx of (twice) the area, but it's fine for our purposes
	}
	return areasum;
}

float PerceptualToneCurve::get_curve_val(float x, float range[2], float lut[], size_t lut_size) {
	float xm = (x - range[0]) / (range[1] - range[0]) * (lut_size - 1);
	if (xm <= 0) {
		return lut[0];
	}
	int idx = (int)xm;
	if (idx >= lut_size-1) {
		return lut[lut_size-1];
	}
	float d = xm - (float)idx; // [0 .. 1]
	return (1.0 - d) * lut[idx] + d * lut[idx+1];
}

// calculate a single value that represents the contrast of the tone curve
float PerceptualToneCurve::calculateToneCurveContrastValue(void) const {

	// find linear y = k*x the best approximates the curve, which is the linear scaling/exposure component that does not contribute any contrast

	// Note: the analysis is made on the gamma encoded curve, as the LUT is linear we make backwards gamma to
	struct find_tc_slope_fun_arg arg = { this };
	float k = find_minimum_interval_halving(find_tc_slope_fun, &arg, 0.1, 5.0, 0.01, 20); // normally found in 8 iterations
	//fprintf(stderr, "average slope: %f\n", k);

        float maxslope = 0;
        {
		// look at midtone slope
		const float xd = 0.07;
		const float tx[] = { 0.30, 0.35, 0.40, 0.45 }; // we only look in the midtone range
		for (int i = 0; i < sizeof(tx)/sizeof(tx[0]); i++) {
			float x0 = tx[i] - xd;
			float y0 = CurveFactory::gamma2(lutToneCurve[CurveFactory::igamma2(x0) * 65535.f] / 65535.f) - k * x0;
			float x1 = tx[i] + xd;
			float y1 = CurveFactory::gamma2(lutToneCurve[CurveFactory::igamma2(x1) * 65535.f] / 65535.f) - k * x1;
			float slope = 1.0 + (y1 - y0) / (x1 - x0);
			if (slope > maxslope) {
				maxslope = slope;
			}
		}

		// look at slope at (light) shadows and (dark) highlights
		float e_maxslope = 0;
		{
			const float tx[] = { 0.20, 0.25, 0.50, 0.55 }; // we only look in the midtone range
			for (int i = 0; i < sizeof(tx)/sizeof(tx[0]); i++) {
				float x0 = tx[i] - xd;
				float y0 = CurveFactory::gamma2(lutToneCurve[CurveFactory::igamma2(x0) * 65535.f] / 65535.f) - k * x0;
				float x1 = tx[i] + xd;
				float y1 = CurveFactory::gamma2(lutToneCurve[CurveFactory::igamma2(x1) * 65535.f] / 65535.f) - k * x1;
				float slope = 1.0 + (y1 - y0) / (x1 - x0);
				if (slope > e_maxslope) {
					e_maxslope = slope;
				}
			}
		}
		//fprintf(stderr, "%.3f %.3f\n", maxslope, e_maxslope);
		// midtone slope is more important for contrast, but weigh in some slope from brights and darks too.
		maxslope = maxslope * 0.7 + e_maxslope * 0.3;
        }
	return maxslope;
}

void PerceptualToneCurve::Apply(float &r, float &g, float &b, PerceptualToneCurveState & state) const {
	float x,y,z;
	cmsCIEXYZ XYZ;
	cmsJCh JCh;

	int thread_idx = 0;
#ifdef _OPENMP
	thread_idx = omp_get_thread_num();
#endif
	if (!state.isProphoto) {
		// convert to prophoto space to make sure the same result is had regardless of working color space
		float newr = state.Working2Prophoto[0][0]*r + state.Working2Prophoto[0][1]*g + state.Working2Prophoto[0][2]*b;
		float newg = state.Working2Prophoto[1][0]*r + state.Working2Prophoto[1][1]*g + state.Working2Prophoto[1][2]*b;
		float newb = state.Working2Prophoto[2][0]*r + state.Working2Prophoto[2][1]*g + state.Working2Prophoto[2][2]*b;
		r = newr;
		g = newg;
		b = newb;
	}

	const AdobeToneCurve& adobeTC = static_cast<const AdobeToneCurve&>((const ToneCurve&)*this);
	float ar = r;
	float ag = g;
	float ab = b;
	adobeTC.Apply(ar, ag, ab);

	if (ar >= 65535.f && ag >= 65535.f && ab >= 65535.f) {
		// clip fast path, will also avoid strange colors of clipped highlights
		r = g = b = 65535.f;
		return;
	}
	if (ar <= 0.f && ag <= 0.f && ab <= 0.f) {
		r = g = b = 0;
		return;
	}

	// ProPhoto constants for luminance, that is xyz_prophoto[1][]
	const float Yr = 0.2880402f;
	const float Yg = 0.7118741f;
	const float Yb = 0.0000857f;

	// we use the Adobe (RGB-HSV hue-stabilized) curve to decide luminance, which generally leads to a less contrasty result
	// compared to a pure luminance curve. We do this to be more compatible with the most popular curves.
        float oldLuminance = r*Yr + g*Yg + b*Yb;
        float newLuminance = ar*Yr + ag*Yg + ab*Yb;
        float Lcoef = newLuminance/oldLuminance;
        r = LIM<float>(r*Lcoef, 0.f, 65535.f);
        g = LIM<float>(g*Lcoef, 0.f, 65535.f);
        b = LIM<float>(b*Lcoef, 0.f, 65535.f);

	// move to JCh so we can modulate chroma based on the global contrast-related chroma scaling factor
	Color::Prophotoxyz(r,g,b,x,y,z);
	XYZ = (cmsCIEXYZ){ .X = x * 100.0f/65535, .Y = y * 100.0f/65535, .Z = z * 100.0f/65535 };
	cmsCIECAM02Forward(h02[thread_idx], &XYZ, &JCh);
	if (!isfinite(JCh.J) || !isfinite(JCh.C) || !isfinite(JCh.h)) {
		// this can happen for dark noise colors or colors outside human gamut. Then we just return the curve's result.
		if (!state.isProphoto) {
			float newr = state.Prophoto2Working[0][0]*r + state.Prophoto2Working[0][1]*g + state.Prophoto2Working[0][2]*b;
			float newg = state.Prophoto2Working[1][0]*r + state.Prophoto2Working[1][1]*g + state.Prophoto2Working[1][2]*b;
			float newb = state.Prophoto2Working[2][0]*r + state.Prophoto2Working[2][1]*g + state.Prophoto2Working[2][2]*b;
			r = newr;
			g = newg;
			b = newb;
		}
		return;
	}

        float cmul = state.cmul_contrast; // chroma scaling factor

	// depending on color, the chroma scaling factor can be fine-tuned below

	{ // decrease chroma scaling sligthly of extremely saturated colors
		float saturated_scale_factor = 0.95;
		const float lolim = 35; // lower limit, below this chroma all colors will keep original chroma scaling factor
		const float hilim = 60; // high limit, above this chroma the chroma scaling factor is multiplied with the saturated scale factor value above
		if (JCh.C < lolim) {
			// chroma is low enough, don't scale
			saturated_scale_factor = 1.0;
		} else if (JCh.C < hilim) {
			// S-curve transition between low and high limit
			float x = (JCh.C - lolim) / (hilim - lolim); // x = [0..1], 0 at lolim, 1 at hilim
			if (x < 0.5) {
				x = 0.5 * powf(2*x, 2);
			} else {
				x = 0.5 + 0.5 * (1-powf(1-2*(x-0.5), 2));
			}
			saturated_scale_factor = 1.0*(1.0-x) + saturated_scale_factor*x;
		} else {
			// do nothing, high saturation color, keep scale factor
		}
		cmul *= saturated_scale_factor;
	}

	{ // increase chroma scaling slightly of shadows
		float nL = CurveFactory::gamma2(newLuminance / 65535); // apply gamma so we make comparison and transition with a more perceptual lightness scale
		float dark_scale_factor = 1.20;
		//float dark_scale_factor = 1.0 + state.debug.p2 / 100.0f;
		const float lolim = 0.15;
		const float hilim = 0.50;
		if (nL < lolim) {
			// do nothing, keep scale factor
		} else if (nL < hilim) {
			// S-curve transition
			float x = (nL - lolim) / (hilim - lolim); // x = [0..1], 0 at lolim, 1 at hilim
			if (x < 0.5) {
				x = 0.5 * powf(2*x, 2);
			} else {
				x = 0.5 + 0.5 * (1-powf(1-2*(x-0.5), 2));
			}
			dark_scale_factor = dark_scale_factor*(1.0-x) + 1.0*x;
		} else {
			dark_scale_factor = 1.0;
		}
		cmul *= dark_scale_factor;
	}

	{ // to avoid strange CIECAM02 chroma errors on close-to-shadow-clipping colors we reduce chroma scaling towards 1.0 for black colors
		float dark_scale_factor = 1.0 / cmul;
		const float lolim = 4;
		const float hilim = 7;
		if (JCh.J < lolim) {
			// do nothing, keep scale factor
		} else if (JCh.J < hilim) {
			// S-curve transition
			float x = (JCh.J - lolim) / (hilim - lolim);
			if (x < 0.5) {
				x = 0.5 * powf(2*x, 2);
			} else {
				x = 0.5 + 0.5 * (1-powf(1-2*(x-0.5), 2));
			}
			dark_scale_factor = dark_scale_factor*(1.0-x) + 1.0*x;
		} else {
			dark_scale_factor = 1.0;
		}
		cmul *= dark_scale_factor;
	}

	JCh.C *= cmul;
	cmsCIECAM02Reverse(h02[thread_idx], &JCh, &XYZ);
	if (!isfinite(XYZ.X) || !isfinite(XYZ.Y) || !isfinite(XYZ.Z)) {
		// can happen for colors on the rim of being outside gamut, that worked without chroma scaling but not with. Then we return only the curve's result.
		if (!state.isProphoto) {
			float newr = state.Prophoto2Working[0][0]*r + state.Prophoto2Working[0][1]*g + state.Prophoto2Working[0][2]*b;
			float newg = state.Prophoto2Working[1][0]*r + state.Prophoto2Working[1][1]*g + state.Prophoto2Working[1][2]*b;
			float newb = state.Prophoto2Working[2][0]*r + state.Prophoto2Working[2][1]*g + state.Prophoto2Working[2][2]*b;
			r = newr;
			g = newg;
			b = newb;
		}
		return;
	}
	Color::xyz2Prophoto(XYZ.X,XYZ.Y,XYZ.Z,r,g,b);
	r *= 655.35;
	g *= 655.35;
	b *= 655.35;
        r = LIM<float>(r, 0.f, 65535.f);
        g = LIM<float>(g, 0.f, 65535.f);
        b = LIM<float>(b, 0.f, 65535.f);

	{ // limit saturation increase in rgb space to avoid severe clipping and flattening in extreme highlights

		// we use the RGB-HSV hue-stable "Adobe" curve as reference. For S-curve contrast it increases
		// saturation greatly, but desaturates extreme highlights and thus provide a smooth transition to
		// the white point. However the desaturation effect is quite strong so we make a weighting
		float ah,as,av,h,s,v;
		Color::rgb2hsv(ar,ag,ab, ah,as,av);
		Color::rgb2hsv(r,g,b, h,s,v);

		float sat_scale = as <= 0.0 ? 1.0 : s / as; // saturation scale compared to Adobe curve
		float keep = 0.2;
		const float lolim = 1.00; // only mix in the Adobe curve if we have increased saturation compared to it
		const float hilim = 1.20;
		if (sat_scale < lolim) {
			// saturation is low enough, don't desaturate
			keep = 1.0;
		} else if (sat_scale < hilim) {
			// S-curve transition
			float x = (sat_scale - lolim) / (hilim - lolim); // x = [0..1], 0 at lolim, 1 at hilim
			if (x < 0.5) {
				x = 0.5 * powf(2*x, 2);
			} else {
				x = 0.5 + 0.5 * (1-powf(1-2*(x-0.5), 2));
			}
			keep = 1.0*(1.0-x) + keep*x;
		} else {
			// do nothing, very high increase, keep minimum amount
		}
		if (keep < 1.0) {
			// mix in some of the Adobe curve result
			r = r * keep + (1.0 - keep) * ar;
			g = g * keep + (1.0 - keep) * ag;
			b = b * keep + (1.0 - keep) * ab;
		}
	}

	if (!state.isProphoto) {
		float newr = state.Prophoto2Working[0][0]*r + state.Prophoto2Working[0][1]*g + state.Prophoto2Working[0][2]*b;
		float newg = state.Prophoto2Working[1][0]*r + state.Prophoto2Working[1][1]*g + state.Prophoto2Working[1][2]*b;
		float newb = state.Prophoto2Working[2][0]*r + state.Prophoto2Working[2][1]*g + state.Prophoto2Working[2][2]*b;
		r = newr;
		g = newg;
		b = newb;
	}
}

cmsContext * PerceptualToneCurve::c02;
cmsHANDLE * PerceptualToneCurve::h02;
float PerceptualToneCurve::cf_range[2];
float PerceptualToneCurve::cf[1000];

void PerceptualToneCurve::init() {

	{ // init ciecam02 state, used for chroma scalings
		cmsViewingConditions vc;
		vc.whitePoint = *cmsD50_XYZ();
		vc.whitePoint.X *= 100;
		vc.whitePoint.Y *= 100;
		vc.whitePoint.Z *= 100;
		vc.Yb = 20;
		vc.La = 20;
		vc.surround = AVG_SURROUND;
		vc.D_value = 1.0;

		int thread_count = 1;
#ifdef _OPENMP
		thread_count = omp_get_max_threads();
#endif
		h02 = (cmsHANDLE *)malloc(sizeof(h02[0]) * (thread_count + 1));
		c02 = (cmsContext *)malloc(sizeof(c02[0]) * (thread_count + 1));
		h02[thread_count] = NULL;
		c02[thread_count] = NULL;
		// little cms requires one state per thread, for thread safety
		for (int i = 0; i < thread_count; i++) {
			c02[i] = cmsCreateContext(NULL, NULL);
			h02[i] = cmsCIECAM02Init(c02[i], &vc);
		}
	}

	{ // init contrast-value-to-chroma-scaling conversion curve

		// contrast value in the left column, chroma scaling in the right. Handles for a spline.
		// Put the columns in a file (without commas) and you can plot the spline with gnuplot: "plot 'curve.txt' smooth csplines"
		// A spline can easily get overshoot issues so if you fine-tune the values here make sure that the resulting spline is smooth afterwards, by
		// plotting it for example.
		const float p[] = {
			0.60, 0.70, // lowest contrast
			0.70, 0.80,
			0.90, 0.94,
			0.99, 1.00,
			1.00, 1.00, // 1.0 (linear curve) to 1.0, no scaling
			1.07, 1.00,
			1.08, 1.00,
			1.11, 1.02,
			1.20, 1.08,
			1.30, 1.12,
			1.80, 1.20,
			2.00, 1.22  // highest contrast
		};

		const size_t in_len = sizeof(p)/sizeof(p[0])/2;
		float in_x[in_len];
		float in_y[in_len];
		for (size_t i = 0; i < in_len; i++) {
			in_x[i]= p[2*i+0];
			in_y[i]= p[2*i+1];
		}
		const size_t out_len = sizeof(cf)/sizeof(cf[0]);
		float out_x[out_len];
		for (size_t i = 0; i < out_len; i++) {
			out_x[i] = in_x[0] + (in_x[in_len-1] - in_x[0]) * (float)i / (out_len-1);
		}
		cubic_spline(in_x, in_y, in_len, out_x, cf, out_len);
		cf_range[0] = in_x[0];
		cf_range[1] = in_x[in_len-1];
	}
}

void PerceptualToneCurve::cleanup() {
	for (int i = 0; h02[i] != NULL; i++) {
		cmsCIECAM02Done(h02[i]);
		cmsDeleteContext(c02[i]);
	}
	free(h02);
	free(c02);
}

void PerceptualToneCurve::initApplyState(PerceptualToneCurveState & state, Glib::ustring workingSpace) const {

	// Get the curve's contrast value, and convert to a chroma scaling
        const float contrast_value = calculateToneCurveContrastValue();
        state.cmul_contrast = get_curve_val(contrast_value, cf_range, cf, sizeof(cf)/sizeof(cf[0]));
	//fprintf(stderr, "contrast value: %f => chroma scaling %f\n", contrast_value, state.cmul_contrast);

	// Create state for converting to/from prophoto (if necessary)
	if (workingSpace == "ProPhoto") {
		state.isProphoto = true;
	} else {
		state.isProphoto = false;
		TMatrix Work = iccStore->workingSpaceMatrix(workingSpace);
		memset(state.Working2Prophoto, 0, sizeof(state.Working2Prophoto));
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
				for (int k=0; k<3; k++)
					state.Working2Prophoto[i][j] += prophoto_xyz[i][k] * Work[k][j];
		Work = iccStore->workingSpaceInverseMatrix (workingSpace);
		memset(state.Prophoto2Working, 0, sizeof(state.Prophoto2Working));
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
				for (int k=0; k<3; k++)
					state.Prophoto2Working[i][j] += Work[i][k] * xyz_prophoto[k][j];
	}
}

}
