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

	defProcParams.setFloat   ("LumaCurveBrightness", 0);
	defProcParams.setFloat   ("LumaCurveContrast", 0);
	FloatList lcurve;
	defProcParams.setFloatList ("LumaCurveCustomCurve",  lcurve);
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
            histogram = new unsigned int [65536];
        memset (histogram, 0, 65536*sizeof(unsigned int));
        for (int i=0; i<sourceImage->height; i++)
            for (int j=0; j<sourceImage->width; j++)
                histogram[CLIP((int)(655.35*sourceImage->cieL[i][j]))]++;

    	float brightness  = procParams->getFloat  ("LumaCurveBrightness");
    	float contrast    = procParams->getFloat  ("LumaCurveContrast");
    	FloatList& ccurve = procParams->getFloatList ("LumaCurveCustomCurve");

    	if (!curve) {
            curve = new float [65536];
            CurveFactory::complexCurve (0.0, 0.0, 0.0, 0.0, brightness, contrast, 0.0, false, ccurve, histogram, curve, NULL, getScale ()==1 ? 1 : 16);
        }
        else if (isTriggerEvent (events))
            CurveFactory::complexCurve (0.0, 0.0, 0.0, 0.0, brightness, contrast, 0.0, false, ccurve, histogram, curve, NULL, getScale ()==1 ? 1 : 16);
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
			targetImage->cieL[i][j] = lutInterp<float,65536> (myCurve, 655.35*sourceImage->cieL[i][j]);
			targetImage->ciea[i][j] = sourceImage->ciea[i][j];
			targetImage->cieb[i][j] = sourceImage->cieb[i][j];
		}
	}
}

}
