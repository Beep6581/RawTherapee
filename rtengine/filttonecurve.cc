/*
 * filttonecurve.cc
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#include "filttonecurve.h"
#include "rtengine.h"
#include "macros.h"
#include "iccstore.h"
#include "improcfuns.h"
#include <string.h>
#include "filterchain.h"
#include "curves.h"

namespace rtengine {

class PreToneCurveFilterDescriptor : public FilterDescriptor {
    public:
        PreToneCurveFilterDescriptor ()
            : FilterDescriptor ("PreToneCurve", MultiImage::RGB, MultiImage::RGB) {
            addTriggerEvent (EvAutoExp);
            addTriggerEvent (EvClip);
        }
        void createAndAddToList (Filter* tail) const {}
};

ToneCurveFilterDescriptor toneCurveFilterDescriptor;
PreToneCurveFilterDescriptor preToneCurveFilterDescriptor;

ToneCurveFilterDescriptor::ToneCurveFilterDescriptor ()
	: FilterDescriptor ("ToneCurve", MultiImage::RGB, MultiImage::RGB) {

    addTriggerEvent (EvBrightness);
	addTriggerEvent (EvContrast);
    addTriggerEvent (EvExpComp);
    addTriggerEvent (EvBlack);
    addTriggerEvent (EvHLCompr);
    addTriggerEvent (EvSHCompr);
    addTriggerEvent (EvToneCurve);
}

void ToneCurveFilterDescriptor::createAndAddToList (Filter* tail) const {

	PreToneCurveFilter* ptcf = new PreToneCurveFilter ();
    tail->addNext (ptcf);
    tail->addNext (new ToneCurveFilter (ptcf));
}

PreToneCurveFilter::PreToneCurveFilter ()
    : Filter (&preToneCurveFilterDescriptor), histogram (NULL) {
}

PreToneCurveFilter::~PreToneCurveFilter () {

    delete [] histogram;
}

void PreToneCurveFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

    Filter* p = getParentFilter ();
    if (!p) {

        // update histogram histogram
        if (!histogram)
            histogram = new unsigned int [65536];

        TMatrix wprof = iccStore.workingSpaceMatrix (procParams->icm.working);
        int mulr = round(32768.0 * wprof[0][1]);
        int mulg = round(32768.0 * wprof[1][1]);
        int mulb = round(32768.0 * wprof[2][1]);

        memset (histogram, 0, 65536*sizeof(int));
        for (int i=0; i<sourceImage->height; i++)
            for (int j=0; j<sourceImage->width; j++) {
                int y = (mulr * sourceImage->r[i][j] + mulg * sourceImage->g[i][j] + mulb * sourceImage->b[i][j]) >> 15;
                histogram[CLIP(y)]++;
            }

        // calculate auto exposure parameters
        if (procParams->toneCurve.autoexp) {
            unsigned int aehist[65536]; int aehistcompr;
            ImageSource* imgsrc = getFilterChain ()->getImageSource ();
            imgsrc->getAEHistogram (aehist, aehistcompr);
            ImProcFunctions::calcAutoExp (aehist, aehistcompr, imgsrc->getDefGain(), procParams->toneCurve.clip, procParams->toneCurve.expcomp, procParams->toneCurve.black);
        }
    }

    // we have to copy image data if input and output are not the same
    if (sourceImage != targetImage)
        targetImage->copyFrom (sourceImage);
}

unsigned int* PreToneCurveFilter::getHistogram () {

    return histogram;
}

ToneCurveFilter::ToneCurveFilter (PreToneCurveFilter* ptcf)
	: Filter (&toneCurveFilterDescriptor), curve (NULL), ptcFilter (ptcf), bchistogram (NULL) {
}

ToneCurveFilter::~ToneCurveFilter () {

    delete [] curve;
    delete [] bchistogram;
}

void ToneCurveFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

    Filter* p = getParentFilter ();
    ImageSource* imgsrc = getFilterChain ()->getImageSource ();

    unsigned int* myCurve;

    // curve is only generated once: in the root filter chain
    if (!p) {
        if (!curve) {
            curve = new unsigned int [65536];
            CurveFactory::complexCurve (procParams->toneCurve.expcomp, procParams->toneCurve.black/65535.0, procParams->toneCurve.hlcompr, procParams->toneCurve.shcompr, procParams->toneCurve.brightness, procParams->toneCurve.contrast, imgsrc->getDefGain(), imgsrc->isRaw() ? 2.2 : 0.0 , true, procParams->toneCurve.curve, ptcFilter->getHistogram(), curve, NULL, getScale ()==1 ? 1 : 16);
        }
        else if (isTriggerEvent (events) || ptcFilter->isTriggerEvent (events))
            CurveFactory::complexCurve (procParams->toneCurve.expcomp, procParams->toneCurve.black/65535.0, procParams->toneCurve.hlcompr, procParams->toneCurve.shcompr, procParams->toneCurve.brightness, procParams->toneCurve.contrast, imgsrc->getDefGain(), imgsrc->isRaw() ? 2.2 : 0.0 , true, procParams->toneCurve.curve, ptcFilter->getHistogram(), curve, NULL, getScale ()==1 ? 1 : 16);
        myCurve = curve;
    }
    else {
        Filter* root = p;
        while (root->getParentFilter())
            root = root->getParentFilter();
        myCurve = ((ToneCurveFilter*)root)->curve;
    }

    // apply curve
    #pragma omp parallel for if (multiThread)
	for (int i=0; i<sourceImage->height; i++) {
		for (int j=0; j<sourceImage->width; j++) {
			targetImage->r[i][j] = myCurve[sourceImage->r[i][j]];
			targetImage->g[i][j] = myCurve[sourceImage->g[i][j]];
			targetImage->b[i][j] = myCurve[sourceImage->b[i][j]];
		}
	}
}

}
