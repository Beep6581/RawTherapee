/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Copyright (c) 2011 Jacques Desmis  <jdesmis@gmail.com>
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

#include "rt_math.h"
//#include <algorithm>

#include "rtengine.h"
#include "improcfun.h"
#include "iccstore.h"
#include "mytime.h"
#include "../rtgui/thresholdselector.h"
#include "curves.h"
#include "color.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace rtengine {

using namespace procparams;

#define SAT(a,b,c) ((float)max(a,b,c)-(float)min(a,b,c))/(float)max(a,b,c)

extern const Settings* settings;

void fillCurveArrayVib(DiagonalCurve* diagCurve, LUTf &outCurve) {

	if (diagCurve) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i=0; i<=0xffff; i++ ) {
			// change to [0,1] range
			// apply custom/parametric/NURBS curve, if any
			// and store result in a temporary array
			outCurve[i] = 65535.f*diagCurve->getVal( double(i)/65535.0 );
		}
	}
	else {
		for (int i=0; i<=0xffff; i++) {
			outCurve[i] = float(i);
		}
	}
}


/*
 * Vibrance correction
 * copyright (c)2011  Jacques Desmis <jdesmis@gmail.com> and Jean-Christophe Frisch <natureh@free.fr>
 *
 */
void ImProcFunctions::vibrance (LabImage* lab) {
	if (!params->vibrance.enabled)
		return;

//	int skip=1; //scale==1 ? 1 : 16;
	bool skinCurveIsSet=false;
	DiagonalCurve* dcurve = NULL;
	dcurve = new DiagonalCurve (params->vibrance.skintonescurve, CURVES_MIN_POLY_POINTS);
	if (dcurve) {
		if (!dcurve->isIdentity()) {
			skinCurveIsSet = true;
		}
		else {
			delete dcurve;
			dcurve = NULL;
		}
	}

	if (!skinCurveIsSet && !params->vibrance.pastels && !params->vibrance.saturated) {
		if (dcurve) {
			delete dcurve;
			dcurve = NULL;
		}
		return;
	}

	const int width = lab->W;
	const int height = lab->H;

#ifdef _DEBUG
	MyTime t1e,t2e;
	t1e.set();
	int negat=0, moreRGB=0, negsat=0, moresat=0;
#endif

	// skin hue curve
	// I use diagonal because I think it's better
	LUTf skin_curve (65536,0);
	if(skinCurveIsSet)
		fillCurveArrayVib(dcurve, skin_curve);
	if (dcurve) {
		delete dcurve;
		dcurve = NULL;
	}


// skin_curve.dump("skin_curve");

	const float chromaPastel = float(params->vibrance.pastels)   / 100.0f;
	const float chromaSatur  = float(params->vibrance.saturated) / 100.0f;
	const float p00=0.07f;
	const float limitpastelsatur =    (static_cast<float>(params->vibrance.psthreshold.value[ThresholdSelector::TS_TOPLEFT])    / 100.0f)*(1.0f-p00) + p00;
	const float maxdp=(limitpastelsatur-p00)/4.0f;
	const float maxds=(1.0-limitpastelsatur)/4.0f;
	const float p0 = p00+maxdp;
	const float p1 = p00+2.0f*maxdp;
	const float p2 = p00+3.0f*maxdp;
	const float s0 = limitpastelsatur + maxds;
	const float s1 = limitpastelsatur + 2.0f*maxds;
	const float s2 = limitpastelsatur + 3.0f*maxds;
	const float transitionweighting = static_cast<float>(params->vibrance.psthreshold.value[ThresholdSelector::TS_BOTTOMLEFT]) / 100.0f;
	float chromamean=0.0f;
	if(chromaPastel != chromaSatur){
		//if sliders pastels and saturated are different: transition with a double linear interpolation: between p2 and limitpastelsatur, and between limitpastelsatur and s0
		//modify the "mean" point in function of double threshold  => differential transition
		chromamean = maxdp * (chromaSatur-chromaPastel) / (s0-p2) + chromaPastel;
		// move chromaMean up or down depending on transitionCtrl
		if (transitionweighting > 0.0f) {
			chromamean = (chromaSatur-chromamean) * transitionweighting + chromamean;
		}
		else if (transitionweighting < 0.0f) {
			chromamean = (chromamean-chromaPastel)  * transitionweighting + chromamean;
		}
	}
	const float chromaPastel_a = (chromaPastel-chromamean)/(p2-limitpastelsatur);
	const float chromaPastel_b = chromaPastel-chromaPastel_a*p2;

	const float chromaSatur_a=(chromaSatur-chromamean)/(s0-limitpastelsatur);
	const float chromaSatur_b=chromaSatur-chromaSatur_a*s0;

	const float dhue=0.15f;//hue transition
	const float dchr=20.0f;//chroma transition
	const float skbeg=-0.05f;//begin hue skin
	const float skend=1.60f;//end hue skin
	const float xx=0.5f;//soft : between 0.3 and 1.0
	const float ask=65535.0f/(skend-skbeg);
	const float bsk=-skbeg*ask;


	const bool highlight = params->toneCurve.hrenabled;//Get the value if "highlight reconstruction" is activated
	const bool protectskins = params->vibrance.protectskins;
	const bool avoidcolorshift = params->vibrance.avoidcolorshift;

	TMatrix wiprof = iccStore->workingSpaceInverseMatrix (params->icm.working);
	//inverse matrix user select
	const double wip[3][3] = {
		{wiprof[0][0],wiprof[0][1],wiprof[0][2]},
		{wiprof[1][0],wiprof[1][1],wiprof[1][2]},
		{wiprof[2][0],wiprof[2][1],wiprof[2][2]}
	};


#ifdef _DEBUG
    MunsellDebugInfo* MunsDebugInfo = NULL;
    if (avoidcolorshift)
    	MunsDebugInfo = new MunsellDebugInfo();
#pragma omp parallel default(shared) firstprivate(lab, width, height, chromaPastel, chromaSatur, highlight, limitpastelsatur, transitionweighting, protectskins, avoidcolorshift, MunsDebugInfo) reduction(+: negat, moreRGB, negsat, moresat) if (multiThread)
#else
#pragma omp parallel default(shared) if (multiThread)
#endif
{

	float sathue[5],sathue2[4];// adjust sat in function of hue
	
/*
	// Fitting limitpastelsatur into the real 0.07->1.0 range
//	limitpastelsatur = limitpastelsatur*(1.0f-p00) + p00;
	float p0,p1,p2;//adapt limit of pyramid to psThreshold
	float s0,s1,s2;
*/

#ifdef _OPENMP
	if (settings->verbose && omp_get_thread_num()==0) {
#else
	if (settings->verbose) {
#endif
		printf("vibrance:  p0=%1.2f  p1=%1.2f  p2=%1.2f  s0=%1.2f s1=%1.2f s2=%1.2f\n", p0,p1,p2,s0,s1,s2);
		printf("           pastel=%f   satur=%f   limit= %1.2f   chromamean=%0.5f\n",1.0f+chromaPastel,1.0f+chromaSatur, limitpastelsatur, chromamean);
	}

#pragma omp for schedule(dynamic, 16)
	for (int i=0; i<height; i++)
		for (int j=0; j<width; j++) {
			float LL=lab->L[i][j]/327.68f;
			float CC=sqrt(SQR(lab->a[i][j])+ SQR(lab->b[i][j]))/327.68f;
			float HH=xatan2f(lab->b[i][j],lab->a[i][j]);
			
			float satredu=1.0f; //reduct sat in function of skin
			if(protectskins) {
				Color::SkinSat (LL, HH, CC, satredu);// for skin colors
			}
			// here we work on Chromaticity and Hue
			// variation of Chromaticity  ==> saturation via RGB
			// Munsell correction, then conversion to Lab
			float Lprov=LL;
			float Chprov=CC;
			float R, G, B;
			float2 sincosval;
			if(CC==0.0f) {
				sincosval.y = 1.f;
				sincosval.x = 0.0f;
			} else {
				sincosval.y = lab->a[i][j]/(CC*327.68f);
				sincosval.x = lab->b[i][j]/(CC*327.68f);
			}

#ifdef _DEBUG
			bool neg=false;
			bool more_rgb=false;
			//gamut control : Lab values are in gamut
			Color::gamutLchonly(HH, sincosval, Lprov, Chprov, R, G, B, wip, highlight, 0.15f, 0.98f, neg, more_rgb);
			if(neg) negat++;
			if(more_rgb) moreRGB++;
#else
			//gamut control : Lab values are in gamut
			Color::gamutLchonly(HH, sincosval, Lprov, Chprov, R, G, B, wip, highlight, 0.15f, 0.98f);
#endif
			if(Chprov > 6.0f) {
				const float saturation=SAT(R,G,B);
				if(saturation>0.0f) {
					if(satredu!=1.0f) {
						// for skin, no differentiation
						sathue [0]=sathue [1]=sathue [2]=sathue [3]=sathue[4]=1.0f;
						sathue2[0]=sathue2[1]=sathue2[2]=sathue2[3]          =1.0f;
					} else {
						//double pyramid: LL and HH
						//I try to take into account: Munsell response (human vision) and Gamut..(less response for red): preferably using Prophoto or WideGamut
						//blue: -1.80 -3.14  green = 2.1 3.14   green-yellow=1.4 2.1  red:0 1.4  blue-purple:-0.7  -1.4   purple: 0 -0.7
						//these values allow a better and differential response
						if(LL < 20.0f) {//more for blue-purple, blue and red modulate
							if     (/*HH> -3.1415f &&*/ HH< -1.5f   ) {sathue[0]=1.3f;sathue[1]=1.2f;sathue[2]=1.1f;sathue[3]=1.05f;sathue[4]=0.4f;sathue2[0]=1.05f;sathue2[1]=1.1f ;sathue2[2]=1.05f;sathue2[3]=1.0f;}//blue
							else if(/*HH>=-1.5f    &&*/ HH< -0.7f   ) {sathue[0]=1.6f;sathue[1]=1.4f;sathue[2]=1.3f;sathue[3]=1.2f ;sathue[4]=0.4f;sathue2[0]=1.2f ;sathue2[1]=1.15f;sathue2[2]=1.1f ;sathue2[3]=1.0f;}//blue purple  1.2 1.1
							else if(/*HH>=-0.7f    &&*/ HH<  0.0f   ) {sathue[0]=1.2f;sathue[1]=1.0f;sathue[2]=1.0f;sathue[3]=1.0f ;sathue[4]=0.4f;sathue2[0]=1.0f ;sathue2[1]=1.0f ;sathue2[2]=1.0f ;sathue2[3]=1.0f;}//purple
				//			else if(  HH>= 0.0f    &&   HH<= 1.4f   ) {sathue[0]=1.1f;sathue[1]=1.1f;sathue[2]=1.1f;sathue[3]=1.0f ;sathue[4]=0.4f;sathue2[0]=1.0f ;sathue2[1]=1.0f ;sathue2[2]=1.0f ;sathue2[3]=1.0f;}//red   0.8 0.7
							else if(/*HH>= 0.0f    &&*/ HH<= 1.4f   ) {sathue[0]=1.3f;sathue[1]=1.2f;sathue[2]=1.1f;sathue[3]=1.0f ;sathue[4]=0.4f;sathue2[0]=1.0f ;sathue2[1]=1.0f ;sathue2[2]=1.0f ;sathue2[3]=1.0f;}//red   0.8 0.7
							else if(/*HH>  1.4f    &&*/ HH<= 2.1f   ) {sathue[0]=1.0f;sathue[1]=1.0f;sathue[2]=1.0f;sathue[3]=1.0f ;sathue[4]=0.4f;sathue2[0]=1.0f ;sathue2[1]=1.0f ;sathue2[2]=1.0f ;sathue2[3]=1.0f;}//green yellow 1.2 1.1
							else /*if(HH>  2.1f    && HH<= 3.1415f)*/ {sathue[0]=1.4f;sathue[1]=1.3f;sathue[2]=1.2f;sathue[3]=1.15f;sathue[4]=0.4f;sathue2[0]=1.15f;sathue2[1]=1.1f ;sathue2[2]=1.05f;sathue2[3]=1.0f;}//green
						}
						else if (LL< 50.0f) {//more for blue and green, less for red and green-yellow
							if     (/*HH> -3.1415f &&*/ HH< -1.5f   ) {sathue[0]=1.5f;sathue[1]=1.4f;sathue[2]=1.3f;sathue[3]=1.2f ;sathue[4]=0.4f;sathue2[0]=1.2f ;sathue2[1]=1.1f ;sathue2[2]=1.05f;sathue2[3]=1.0f;}//blue
							else if(/*HH>=-1.5f    &&*/ HH< -0.7f   ) {sathue[0]=1.3f;sathue[1]=1.2f;sathue[2]=1.1f;sathue[3]=1.05f;sathue[4]=0.4f;sathue2[0]=1.05f;sathue2[1]=1.05f;sathue2[2]=1.0f ;sathue2[3]=1.0f;}//blue purple  1.2 1.1
							else if(/*HH>=-0.7f    &&*/ HH<  0.0f   ) {sathue[0]=1.2f;sathue[1]=1.0f;sathue[2]=1.0f;sathue[3]=1.0f ;sathue[4]=0.4f;sathue2[0]=1.0f ;sathue2[1]=1.0f ;sathue2[2]=1.0f ;sathue2[3]=1.0f;}//purple
				//			else if(  HH>= 0.0f    &&   HH<= 1.4f   ) {sathue[0]=0.8f;sathue[1]=0.8f;sathue[2]=0.8f;sathue[3]=0.8f ;sathue[4]=0.4f;sathue2[0]=0.8f ;sathue2[1]=0.8f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//red   0.8 0.7
							else if(/*HH>= 0.0f    &&*/ HH<= 1.4f   ) {sathue[0]=1.1f;sathue[1]=1.0f;sathue[2]=0.9f;sathue[3]=0.8f ;sathue[4]=0.4f;sathue2[0]=0.8f ;sathue2[1]=0.8f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//red   0.8 0.7
							else if(/*HH>  1.4f    &&*/ HH<= 2.1f   ) {sathue[0]=1.1f;sathue[1]=1.1f;sathue[2]=1.1f;sathue[3]=1.05f;sathue[4]=0.4f;sathue2[0]=0.9f ;sathue2[1]=0.8f ;sathue2[2]=0.7f ;sathue2[3]=0.6f;}//green yellow 1.2 1.1
							else /*if(HH>  2.1f    && HH<= 3.1415f)*/ {sathue[0]=1.5f;sathue[1]=1.4f;sathue[2]=1.3f;sathue[3]=1.2f ;sathue[4]=0.4f;sathue2[0]=1.2f ;sathue2[1]=1.1f ;sathue2[2]=1.05f;sathue2[3]=1.0f;}//green

						}
						else if (LL< 80.0f) {//more for green, less for red and green-yellow
							if     (/*HH> -3.1415f &&*/ HH< -1.5f   ) {sathue[0]=1.3f;sathue[1]=1.2f;sathue[2]=1.15f;sathue[3]=1.1f ;sathue[4]=0.3f;sathue2[0]=1.1f ;sathue2[1]=1.1f ;sathue2[2]=1.05f;sathue2[3]=1.0f;}//blue
							else if(/*HH>=-1.5f    &&*/ HH< -0.7f   ) {sathue[0]=1.3f;sathue[1]=1.2f;sathue[2]=1.15f;sathue[3]=1.1f ;sathue[4]=0.3f;sathue2[0]=1.1f ;sathue2[1]=1.05f;sathue2[2]=1.0f ;sathue2[3]=1.0f;}//blue purple  1.2 1.1
							else if(/*HH>=-0.7f    &&*/ HH<  0.0f   ) {sathue[0]=1.2f;sathue[1]=1.0f;sathue[2]=1.0f ;sathue[3]=1.0f ;sathue[4]=0.3f;sathue2[0]=1.0f ;sathue2[1]=1.0f ;sathue2[2]=1.0f ;sathue2[3]=1.0f;}//purple
				//			else if(  HH>= 0.0f    &&   HH<= 1.4f   ) {sathue[0]=0.8f;sathue[1]=0.8f;sathue[2]=0.8f ;sathue[3]=0.8f ;sathue[4]=0.3f;sathue2[0]=0.8f ;sathue2[1]=0.8f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//red   0.8 0.7
							else if(/*HH>= 0.0f    &&*/ HH<= 1.4f   ) {sathue[0]=1.1f;sathue[1]=1.0f;sathue[2]=0.9f ;sathue[3]=0.8f ;sathue[4]=0.3f;sathue2[0]=0.8f ;sathue2[1]=0.8f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//red   0.8 0.7
							else if(/*HH>  1.4f    &&*/ HH<= 2.1f   ) {sathue[0]=1.3f;sathue[1]=1.2f;sathue[2]=1.1f ;sathue[3]=1.05f;sathue[4]=0.3f;sathue2[0]=1.0f ;sathue2[1]=0.9f ;sathue2[2]=0.8f ;sathue2[3]=0.7f;}//green yellow 1.2 1.1
							else /*if(HH>  2.1f    && HH<= 3.1415f)*/ {sathue[0]=1.6f;sathue[1]=1.4f;sathue[2]=1.3f ;sathue[3]=1.25f;sathue[4]=0.3f;sathue2[0]=1.25f;sathue2[1]=1.2f ;sathue2[2]=1.15f;sathue2[3]=1.05f;}//green - even with Prophoto green are too "little"  1.5 1.3
						}
						else /*if (LL>=80.0f)*/ {//more for green-yellow, less for red and purple
							if     (/*HH> -3.1415f &&*/ HH< -1.5f   ) {sathue[0]=1.0f;sathue[1]=1.0f;sathue[2]=0.9f;sathue[3]=0.8f;sathue[4]=0.2f;sathue2[0]=0.8f;sathue2[1]=0.8f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//blue
							else if(/*HH>=-1.5f    &&*/ HH< -0.7f   ) {sathue[0]=1.0f;sathue[1]=1.0f;sathue[2]=0.9f;sathue[3]=0.8f;sathue[4]=0.2f;sathue2[0]=0.8f;sathue2[1]=0.8f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//blue purple  1.2 1.1
							else if(/*HH>=-0.7f    &&*/ HH<  0.0f   ) {sathue[0]=1.2f;sathue[1]=1.0f;sathue[2]=1.0f;sathue[3]=0.9f;sathue[4]=0.2f;sathue2[0]=0.9f;sathue2[1]=0.9f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//purple
				//			else if(  HH>= 0.0f    &&   HH<= 1.4f   ) {sathue[0]=0.8f;sathue[1]=0.8f;sathue[2]=0.8f;sathue[3]=0.8f;sathue[4]=0.2f;sathue2[0]=0.8f;sathue2[1]=0.8f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//red   0.8 0.7
							else if(/*HH>= 0.0f    &&*/ HH<= 1.4f   ) {sathue[0]=1.1f;sathue[1]=1.0f;sathue[2]=0.9f;sathue[3]=0.8f;sathue[4]=0.2f;sathue2[0]=0.8f;sathue2[1]=0.8f ;sathue2[2]=0.8f ;sathue2[3]=0.8f;}//red   0.8 0.7
							else if(/*HH>  1.4f    &&*/ HH<= 2.1f   ) {sathue[0]=1.6f;sathue[1]=1.5f;sathue[2]=1.4f;sathue[3]=1.2f;sathue[4]=0.2f;sathue2[0]=1.1f;sathue2[1]=1.05f;sathue2[2]=1.0f ;sathue2[3]=1.0f;}//green yellow 1.2 1.1
							else /*if(HH>  2.1f    && HH<= 3.1415f)*/ {sathue[0]=1.4f;sathue[1]=1.3f;sathue[2]=1.2f;sathue[3]=1.1f;sathue[4]=0.2f;sathue2[0]=1.1f;sathue2[1]=1.05f;sathue2[2]=1.05f;sathue2[3]=1.0f;}//green
						}
					}
					float chmodpastel,chmodsat;
					// variables to improve transitions
					float pa, pb;// transition = pa*saturation + pb
					float chl00 = chromaPastel*satredu*sathue[4];
					float chl0  = chromaPastel*satredu*sathue[0];
					float chl1  = chromaPastel*satredu*sathue[1];
					float chl2  = chromaPastel*satredu*sathue[2];
					float chl3  = chromaPastel*satredu*sathue[3];
					float chs0  = chromaSatur*satredu*sathue2[0];
					float chs1  = chromaSatur*satredu*sathue2[1];
					float chs2  = chromaSatur*satredu*sathue2[2];
					float chs3  = chromaSatur*satredu*sathue2[3];
					float s3    = 1.0f;
					// We handle only positive values here ;  improve transitions
					if      (saturation < p00)               chmodpastel = chl00 ;   //neutral tones
					else if (saturation < p0 )               { pa=(chl00-chl0)/(p00-p0);              pb=chl00-pa*p00;             chmodpastel = pa*saturation + pb; }
					else if (saturation < p1)                { pa=(chl0-chl1)/(p0-p1);                pb=chl0-pa*p0;               chmodpastel = pa*saturation + pb; }
					else if (saturation < p2)                { pa=(chl1-chl2)/(p1-p2);                pb=chl1-pa*p1;               chmodpastel = pa*saturation + pb; }
					else if (saturation < limitpastelsatur)  { pa=(chl2- chl3)/(p2-limitpastelsatur); pb=chl2-pa*p2;               chmodpastel = pa*saturation + pb; }
					else if (saturation < s0)                { pa=(chl3-chs0)/(limitpastelsatur-s0) ; pb=chl3-pa*limitpastelsatur; chmodsat    = pa*saturation + pb; }
					else if (saturation < s1)                { pa=(chs0-chs1)/(s0-s1);                pb=chs0-pa*s0;               chmodsat    = pa*saturation + pb; }
					else if (saturation < s2)                { pa=(chs1-chs2)/(s1-s2);                pb=chs1-pa*s1;               chmodsat    = pa*saturation + pb; }
					else                                     { pa=(chs2-chs3)/(s2-s3);                pb=chs2-pa*s2;               chmodsat    = pa*saturation + pb; }

					if(chromaPastel != chromaSatur){

						// Pastels
						if(saturation > p2 && saturation < limitpastelsatur) {
							float newchromaPastel = chromaPastel_a*saturation + chromaPastel_b;
							chmodpastel = newchromaPastel*satredu*sathue[3];
						}

						// Saturated
						if(saturation < s0 && saturation >=limitpastelsatur) {
							float newchromaSatur=chromaSatur_a*saturation + chromaSatur_b;
							chmodsat = newchromaSatur*satredu*sathue2[0];
						}
					}// end transition

					if (saturation <= limitpastelsatur) {
						if (chmodpastel >  2.0f )
							chmodpastel = 2.0f;   //avoid too big values
						else if(chmodpastel < -0.93f)
							chmodpastel =-0.93f;  //avoid negative values

						Chprov *=(1.0f+chmodpastel);
						if(Chprov<6.0f)
							Chprov=6.0f;
					}
					else { //if (saturation > limitpastelsatur)
						if (chmodsat >  1.8f )
							chmodsat = 1.8f;        //saturated
						else if(chmodsat < -0.93f)
							chmodsat =-0.93f;

						Chprov *= 1.0f+chmodsat;
						if(Chprov < 6.0f)
							Chprov=6.0f;
					}					
				}
			}

			bool hhModified = false;
			// Vibrance's Skin curve
			if(skinCurveIsSet) {
				if (HH>skbeg && HH<skend) {
					if(Chprov < 60.0f) {//skin hue  : todo ==> transition
						float HHsk=ask*HH+bsk;
						float Hn=(skin_curve[HHsk]-bsk)/ask;
						float Hc=(Hn*xx+HH*(1.0f-xx));
						HH=Hc;
						hhModified = true;
					}
					else if(Chprov < (60.0f+dchr)) {//transition chroma
						float HHsk=ask*HH+bsk;
						float Hn=(skin_curve[HHsk]-bsk)/ask;
						float Hc=(Hn*xx+HH*(1.0f-xx));
						float aa= (HH-Hc)/dchr ; float bb= HH-(60.0f+dchr)*aa;
						HH=aa*Chprov+bb;
						hhModified = true;
					}
				}
				//transition hue
				else if(HH>(skbeg-dhue) && HH<=skbeg && Chprov < (60.0f+dchr*0.5f)) {
					float HHsk=ask*skbeg+bsk;
					float Hn=(skin_curve[HHsk]-bsk)/ask;
					float Hcc=(Hn*xx+skbeg*(1.0f-xx));
					float adh=(Hcc-(skbeg-dhue))/(dhue);
					float bdh=Hcc-adh*skbeg;
					HH=adh*HH+bdh;
					hhModified = true;
				}
				else if(HH>=skend && HH<(skend+dhue) && Chprov < (60.0f+dchr*0.5f)) {
					float HHsk=ask*skend+bsk;
					float Hn=(skin_curve[HHsk]-bsk)/ask;
					float Hcc=(Hn*xx+skend*(1.0f-xx));
					float adh=(skend+dhue-Hcc)/(dhue);
					float bdh=Hcc-adh*skend;
					HH=adh*HH+bdh;
					hhModified = true;
				}
			} // end skin hue

			//Munsell correction
//			float2 sincosval;
			if(!avoidcolorshift && hhModified)
				 sincosval = xsincosf(HH);
			float aprovn,bprovn;
			bool inGamut;
			do {
				inGamut=true;
				if(avoidcolorshift) {
					float correctionHue=0.0f;
					float correctlum=0.0f;

#ifdef _DEBUG
					Color::AllMunsellLch(/*lumaMuns*/false, Lprov,Lprov,HH,Chprov,CC,correctionHue,correctlum, MunsDebugInfo);
#else
					Color::AllMunsellLch(/*lumaMuns*/false, Lprov,Lprov,HH,Chprov,CC,correctionHue,correctlum);
#endif
					if(correctionHue != 0.f || hhModified) {
						sincosval = xsincosf(HH+correctionHue);
						hhModified = false;
					}
				}
				aprovn=Chprov*sincosval.y;
				bprovn=Chprov*sincosval.x;

				float fyy = (0.00862069f *Lprov )+ 0.137932f;
				float fxx = (0.002f * aprovn) + fyy;
				float fzz = fyy - (0.005f * bprovn);
				float xx_ = 65535.f * Color::f2xyz(fxx)*Color::D50x;
			//	float yy_ = 65535.0f * Color::f2xyz(fyy);
				float zz_ = 65535.f * Color::f2xyz(fzz)*Color::D50z;
				float yy_ = 65535.f * ((Lprov>Color::epskap) ? fyy*fyy*fyy : Lprov/Color::kappa);

				Color::xyz2rgb(xx_,yy_,zz_,R,G,B,wip);

				if(R<0.0f || G<0.0f || B<0.0f) {
#ifdef _DEBUG
					negsat++;
#endif
					Chprov *= 0.98f;
					inGamut = false;
				}

				// if "highlight reconstruction" enabled don't control Gamut for highlights
				if((!highlight) && (R>65535.0f || G>65535.0f || B>65535.0f)) {
#ifdef _DEBUG
					moresat++;
#endif
					Chprov *= 0.98f;
					inGamut = false;
				}
			} while (!inGamut);
			
			//put new values in Lab
			lab->L[i][j]=Lprov*327.68f;
			lab->a[i][j]=aprovn*327.68f;
			lab->b[i][j]=bprovn*327.68f;
	}

} // end of parallelization

#ifdef _DEBUG
	t2e.set();
	if (settings->verbose) {
		printf("Vibrance (performed in %d usec):\n", t2e.etime(t1e));
		printf("   Gamut: G1negat=%iiter G165535=%iiter G2negsat=%iiter G265535=%iiter\n",negat,moreRGB,negsat,moresat);
		if (MunsDebugInfo)
        printf("   Munsell chrominance: MaxBP=%1.2frad  MaxRY=%1.2frad  MaxGY=%1.2frad  MaxRP=%1.2frad  depass=%i\n", MunsDebugInfo->maxdhue[0], MunsDebugInfo->maxdhue[1], MunsDebugInfo->maxdhue[2], MunsDebugInfo->maxdhue[3], MunsDebugInfo->depass);
	}
	if (MunsDebugInfo)
		delete MunsDebugInfo;
#endif

}


}
