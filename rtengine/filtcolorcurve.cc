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

void ColorCurveFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new ColorCurveFilter ());
}

ColorCurveFilter::ColorCurveFilter ()
	: Filter (&colorCurveFilterDescriptor), curve (NULL) {
}

ColorCurveFilter::~ColorCurveFilter () {

    delete [] curve;
}

void ColorCurveFilter::generateCurve () {

    if (!curve)
        curve = new double [107701];

    // generate color multiplier lookup table
    double boost_a = (procParams->colorBoost.amount + 100.0) / 100.0;
    double boost_b = (procParams->colorBoost.amount + 100.0) / 100.0;
    double c = std::max (boost_a, boost_b);

    double d = 1077 * procParams->colorBoost.saturationlimit / 100.0;
    double alpha = 0.5;
    double threshold1 = alpha * d;
    double threshold2 = c*d*(alpha+1.0) - d;
    #pragma omp parallel for if (multiThread)
    for (int i=0; i<107700; i++) { // lookup table stores multipliers with a 0.25 chrominance resolution
        double chrominance = (double)i/1000;
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

void ColorCurveFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

    double* myCurve;

    // apply color curve on the image
    double boost_a = (procParams->colorBoost.amount + 100.0) / 100.0;
    double boost_b = (procParams->colorBoost.amount + 100.0) / 100.0;
    double cmul, amul = 1.0, bmul = 1.0;
    if (boost_a > boost_b) {
        cmul = boost_a;
        if (boost_a > 0)
            bmul = boost_b / boost_a;
    }
    else {
        cmul = boost_b;
        if (boost_b > 0)
            amul = boost_a / boost_b;
    }

    // curve is only generated once: in the root filter chain
    if (procParams->colorBoost.enable_saturationlimiter && cmul > 1) {
        Filter* p = getParentFilter ();
        if (!p) {
            generateCurve ();
            myCurve = curve;
        }
        else {
            Filter* root = p;
            while (root->getParentFilter())
                root = root->getParentFilter();
            myCurve = ((ColorCurveFilter*)root)->curve;
        }
    }

    double shift_a = procParams->colorShift.a * 16384 / 500;
    double shift_b = procParams->colorShift.b * 16384 / 200;

    #pragma omp parallel for if (multiThread)
    for (int i=0; i<sourceImage->height; i++) {
		for (int j=0; j<sourceImage->width; j++) {

            int oL = sourceImage->cieL[i][j];
		    int oa = sourceImage->ciea[i][j];
            int ob = sourceImage->cieb[i][j];

            double allowed_cmul = cmul;
            if (procParams->colorBoost.enable_saturationlimiter && cmul > 1) {
                double sa = (oa + shift_a)*500/16384;
                double sb = (ob + shift_b)*200/16384;
                int chroma = (int)(1000.0 * sqrt(sa*sa + sb*sb));
                allowed_cmul = curve [MIN(chroma,107700)];
            }

            if (allowed_cmul >= 1.0 && procParams->colorBoost.avoidclip) {
                double cclip = 100000;
                double cr = tightestroot ((double)oL/655.35, (double)(oa + shift_a)*500/16384*amul, (double)(ob + shift_b)*200/16384*bmul, 3.079935, -1.5371515, -0.54278342);
                double cg = tightestroot ((double)oL/655.35, (double)(oa + shift_a)*500/16384*amul, (double)(ob + shift_b)*200/16384*bmul, -0.92123418, 1.87599, 0.04524418);
                double cb = tightestroot ((double)oL/655.35, (double)(oa + shift_a)*500/16384*amul, (double)(ob + shift_b)*200/16384*bmul, 0.052889682, -0.20404134, 1.15115166);
                if (cr>1.0 && cr<cclip) cclip = cr;
                if (cg>1.0 && cg<cclip) cclip = cg;
                if (cb>1.0 && cb<cclip) cclip = cb;
                if (cclip<100000) {
                    allowed_cmul = -cclip + 2.0*cclip / (1.0+exp(-2.0*allowed_cmul/cclip));
                    if (allowed_cmul<1.0)
                        allowed_cmul = 1.0;
                }
            }

            int nna = (int)((oa + shift_a) * allowed_cmul * amul);
            int nnb = (int)((ob + shift_b) * allowed_cmul * bmul);

            targetImage->ciea[i][j] = CLIPTO(nna,-32768,32767);
            targetImage->cieb[i][j] = CLIPTO(nnb,-32768,32767);
            targetImage->cieL[i][j] = sourceImage->cieL[i][j];
        }
	}
}

}
