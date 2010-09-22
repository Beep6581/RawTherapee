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

void LumaCurveFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

    Filter* p = getParentFilter ();

    unsigned int* myCurve;

    // curve and histogram is only generated once: in the root filter chain
    if (!p) {

        if (!histogram)
            histogram = new unsigned int [65536];
        memset (histogram, 0, 65536*sizeof(unsigned int));
        for (int i=0; i<sourceImage->height; i++)
            for (int j=0; j<sourceImage->width; j++)
                histogram[sourceImage->cieL[i][j]]++;

        if (!curve) {
            curve = new unsigned int [65536];
            CurveFactory::complexCurve (0.0, 0.0, 0.0, 0.0, procParams->lumaCurve.brightness, procParams->lumaCurve.contrast, 0.0, 0.0, false, procParams->lumaCurve.curve, histogram, curve, NULL, getScale ()==1 ? 1 : 16);
        }
        else if (isTriggerEvent (events))
            CurveFactory::complexCurve (0.0, 0.0, 0.0, 0.0, procParams->lumaCurve.brightness, procParams->lumaCurve.contrast, 0.0, 0.0, false, procParams->lumaCurve.curve, histogram, curve, NULL, getScale ()==1 ? 1 : 16);
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
			targetImage->cieL[i][j] = myCurve[sourceImage->cieL[i][j]];
			targetImage->ciea[i][j] = sourceImage->ciea[i][j];
			targetImage->cieb[i][j] = sourceImage->cieb[i][j];
		}
	}
}

}
