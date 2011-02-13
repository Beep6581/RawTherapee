/*
 * filtlumacurve
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#include "filtlumacurve.h"
#include "rtengine.h"
#include "macros.h"
#include <string.h>
#include "curves.h"
#include "util.h"

#define LC_LUTSIZE 	2*65536
#define LC_LUTSCALE 65536

namespace rtengine {

LumaCurveFilterDescriptor lumaCurveFilterDescriptor;

LumaCurveFilterDescriptor::LumaCurveFilterDescriptor ()
	: FilterDescriptor ("LuminanceCurve", MultiImage::Lab, MultiImage::Lab) {

    addTriggerEvent (EvLBrightness);
	addTriggerEvent (EvLContrast);
    addTriggerEvent (EvLBlack);
    addTriggerEvent (EvLHLCompr);
    addTriggerEvent (EvLSHCompr);
    addTriggerEvent (EvLCurve);
}

void LumaCurveFilterDescriptor::getDefaultParameters (ProcParams& defProcParams) const {

	defProcParams.setFloat   	("LumaCurve", "Brightness", 0);
	defProcParams.setFloat   	("LumaCurve", "Contrast", 0);
	FloatVector lcurve;
	defProcParams.setFloatVector ("LumaCurve", "CustomCurve",  lcurve);
}

void LumaCurveFilterDescriptor::createAndAddToList (Filter* tail) const {

    tail->addNext (new LumaCurveFilter ());
}

LumaCurveFilter::LumaCurveFilter ()
	: Filter (&lumaCurveFilterDescriptor), curve (NULL), histogram (NULL) {
}

LumaCurveFilter::~LumaCurveFilter () {

    delete [] curve;
    delete [] histogram;
}

void LumaCurveFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer) {

    Filter* p = getParentFilter ();

    float* myCurve;

    // curve and histogram is only generated once: in the root filter chain
    if (!p) {
    	if (!histogram)
            histogram = new unsigned int [LC_LUTSIZE];
        memset (histogram, 0, LC_LUTSIZE*sizeof(unsigned int));
        for (int i=0; i<sourceImage->height; i++)
            for (int j=0; j<sourceImage->width; j++)
                histogram[CLIPTO((int)(LC_LUTSCALE*sourceImage->cieL[i][j]/100.0),0,LC_LUTSIZE-1)]++;

    	float brightness    = procParams->getFloat  ("LumaCurve", "Brightness");
    	float contrast      = procParams->getFloat  ("LumaCurve", "Contrast");
    	FloatVector& ccurve = procParams->getFloatVector ("LumaCurve", "CustomCurve");

    	if (!curve) {
            curve = new float [LC_LUTSIZE];
            CurveFactory::complexCurve (0.0, 0.0, 0.0, 0.0, brightness, contrast, 0.0, false, ccurve, histogram, LC_LUTSIZE, LC_LUTSCALE, curve, NULL, getScale ()==1 ? 1 : 16);
        }
        else if (isTriggerEvent (events))
            CurveFactory::complexCurve (0.0, 0.0, 0.0, 0.0, brightness, contrast, 0.0, false, ccurve, histogram, LC_LUTSIZE, LC_LUTSCALE, curve, NULL, getScale ()==1 ? 1 : 16);
        myCurve = curve;
    }
    else {
        Filter* root = p;
        while (root->getParentFilter())
            root = root->getParentFilter();
        myCurve = ((LumaCurveFilter*)root)->curve;
    }

    // apply curve
    #pragma omp parallel for if (multiThread)
	for (int i=0; i<sourceImage->height; i++) {
		for (int j=0; j<sourceImage->width; j++) {
			targetImage->cieL[i][j] = 100.0 * lutInterp<float,LC_LUTSIZE> (myCurve, LC_LUTSCALE*sourceImage->cieL[i][j]/100.0);
			targetImage->ciea[i][j] = sourceImage->ciea[i][j];
			targetImage->cieb[i][j] = sourceImage->cieb[i][j];
		}
	}
}

}
