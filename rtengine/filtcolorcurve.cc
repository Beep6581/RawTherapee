/*
 * filtcolorcurve.cc
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#include "filtcolorcurve.h"
#include "rtengine.h"
#include "macros.h"
#include "colorclip.h"
#include "util.h"

#define CURVE_LUTSIZE 50000
#define CURVE_LUTSCALE 100.0

namespace rtengine {

ColorCurveFilterDescriptor colorCurveFilterDescriptor;

ColorCurveFilterDescriptor::ColorCurveFilterDescriptor ()
	: FilterDescriptor ("ColorCurve", MultiImage::Lab, MultiImage::Lab) {

	addTriggerEvent (EvCBBoost);
	addTriggerEvent (EvCBAvoidClip);
    addTriggerEvent (EvCBSatLimiter);
    addTriggerEvent (EvCBSatLimit);
    addTriggerEvent (EvCShiftA);
    addTriggerEvent (EvCShiftB);
}

void ColorCurveFilterDescriptor::getDefaultParameters (ProcParams& defProcParams) const {

	defProcParams.setInteger ("ColorBoostAmount", 0);
	defProcParams.setBoolean ("ColorBoostAvoidClip", false);
	defProcParams.setBoolean ("ColorBoostEnableSatLimiter", false);
	defProcParams.setDouble  ("ColorBoostSaturationLimit",  50);
	defProcParams.setDouble  ("ColorShiftCieaChannel",  0.0);
	defProcParams.setDouble  ("ColorShiftCiebChannel",  0.0);
}

void ColorCurveFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new ColorCurveFilter ());
}

ColorCurveFilter::ColorCurveFilter ()
	: Filter (&colorCurveFilterDescriptor), curve (NULL) {
}

ColorCurveFilter::~ColorCurveFilter () {

    delete [] curve;
}

void ColorCurveFilter::generateCurve (float boost, float limit) {

    if (!curve)
        curve = new float [CURVE_LUTSIZE];

    // generate color multiplier lookup table

    float c = boost;
    float d = limit;
    float alpha = 0.5;
    float threshold1 = alpha * d;
    float threshold2 = c*d*(alpha+1.0) - d;
    #pragma omp parallel for if (multiThread)
    for (int i=0; i<CURVE_LUTSIZE; i++) { // lookup table stores multipliers with a 0.25 chrominance resolution
        double chrominance = (double)i/CURVE_LUTSCALE;
        if (chrominance < threshold1)
            curve[i] = c;
        else if (chrominance < d)
            curve[i] = (c / (2.0*d*(alpha-1.0)) * (chrominance-d)*(chrominance-d) + c*d/2.0 * (alpha+1.0) ) / chrominance;
        else if (chrominance < threshold2)
            curve[i] = (1.0 / (2.0*d*(c*(alpha+1.0)-2.0)) * (chrominance-d)*(chrominance-d) + c*d/2.0 * (alpha+1.0) ) / chrominance;
        else
            curve[i] = 1.0;
    }
}

void ColorCurveFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer) {

	float boost =  (procParams->getInteger ("ColorBoostAmount") + 100.0) / 100.0;
	float shift_a = procParams->getDouble  ("ColorShiftCieaChannel");
	float shift_b = procParams->getDouble  ("ColorShiftCiebChannel");
	float satlimit = procParams->getDouble  ("ColorBoostSaturationLimit");
	bool  enalimiter = procParams->getBoolean ("ColorBoostEnableSatLimiter");
	bool  avoidclip = procParams->getBoolean ("ColorBoostAvoidClip");

	float* myCurve;

    // curve is only generated once: in the root filter chain
    if (enalimiter && boost > 1) {
        Filter* p = getParentFilter ();
        if (!p) {
            generateCurve (boost, satlimit);
            myCurve = curve;
        }
        else {
            Filter* root = p;
            while (root->getParentFilter())
                root = root->getParentFilter();
            myCurve = ((ColorCurveFilter*)root)->curve;
        }
    }

    #pragma omp parallel for if (multiThread)
    for (int i=0; i<sourceImage->height; i++) {
		for (int j=0; j<sourceImage->width; j++) {

            float oL = sourceImage->cieL[i][j];
		    float oa = sourceImage->ciea[i][j] + shift_a;
            float ob = sourceImage->cieb[i][j] + shift_b;

            float allowed_mul = boost;
            if (enalimiter && boost > 1) {
                float chroma = CURVE_LUTSCALE * sqrt(oa*oa + ob*ob);
                allowed_mul =  lutInterp<float,CURVE_LUTSIZE>(curve, chroma);
            }

            if (allowed_mul >= 1.0 && procParams->colorBoost.avoidclip) {
                float cclip = FLT_MAX;
                float cr = tightestroot (oL, oa, ob, 3.079935, -1.5371515, -0.54278342);
                float cg = tightestroot (oL, oa, ob, -0.92123418, 1.87599, 0.04524418);
                float cb = tightestroot (oL, oa, ob, 0.052889682, -0.20404134, 1.15115166);
                if (cr>1.0 && cr<cclip) cclip = cr;
                if (cg>1.0 && cg<cclip) cclip = cg;
                if (cb>1.0 && cb<cclip) cclip = cb;
                if (cclip<FLT_MAX) {
                    allowed_mul = -cclip + 2.0*cclip / (1.0+exp(-2.0*allowed_mul/cclip));
                    if (allowed_mul<1.0)
                        allowed_mul = 1.0;
                }
            }

            targetImage->ciea[i][j] = oa * allowed_mul;
            targetImage->cieb[i][j] = ob * allowed_mul;
            targetImage->cieL[i][j] = sourceImage->cieL[i][j];
        }
	}
}

}
