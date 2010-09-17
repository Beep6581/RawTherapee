/*
 * filtcolordenoise.cc
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#include "filtcolordenoise.h"
#include "rtengine.h"
#include "macros.h"
#include "bilateral2.h"

namespace rtengine {

ColorDenoiseFilterDescriptor colorDenoiseFilterDescriptor;

ColorDenoiseFilterDescriptor::ColorDenoiseFilterDescriptor ()
	: FilterDescriptor ("ColorDenoiser", MultiImage::Lab, MultiImage::Lab) {

    addTriggerEvent (EvCDNEnabled);
    addTriggerEvent (EvCDNRadius);
    addTriggerEvent (EvCDNEdgeTolerance);
    addTriggerEvent (EvCDNEdgeSensitive);
}

void ColorDenoiseFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new ColorDenoiseFilter ());
}

ColorDenoiseFilter::ColorDenoiseFilter ()
	: Filter (&colorDenoiseFilterDescriptor) {
}

void ColorDenoiseFilter::getReqiredBufferSize (int& w, int& h) {

    if (params->colorDenoise.enabled && sourceImage->width>=8 && sourceImage->height>=8) {
        int sw = getSourceImageView().getPixelWidth ();
        int sh = getSourceImageView().getPixelHeight ();
        if (sh > sw) {
            w = 2;  // since we need double buffer and rt gives int buffer
            h = sh*omp_get_max_threads();
        }
        else {
            w = sw*omp_get_max_threads();
            h = 2;
        }
    }
    else {
        w = sw;
        h = sh;
    }
}

void ColorDenoiseFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

    if (params->colorDenoise.enabled && sourceImage->width>=8 && sourceImage->height>=8) {

        double scale = getScale ();

        gaussHorizontal<short> (&sourceImage->getBufferView(sourceImage->ciea), &targetImage->getBufferView(targetImage->ciea), (double*)buffer->data, procParams->colorDenoise.amount / 10.0 / scale, multiThread);
        gaussHorizontal<short> (&sourceImage->getBufferView(sourceImage->cieb), &targetImage->getBufferView(targetImage->cieb), (double*)buffer->data, procParams->colorDenoise.amount / 10.0 / scale, multiThread);
        gaussVertical<short>   (&targetImage->getBufferView(targetImage->ciea), &targetImage->getBufferView(targetImage->ciea), (double*)buffer->data, procParams->colorDenoise.amount / 10.0 / scale, multiThread);
        gaussVertical<short>   (&targetImage->getBufferView(targetImage->cieb), &targetImage->getBufferView(targetImage->cieb), (double*)buffer->data, procParams->colorDenoise.amount / 10.0 / scale, multiThread);
    }
    else if (sourceImage != targetImage)
        targetImage->copyFrom (sourceImage);
}

}
