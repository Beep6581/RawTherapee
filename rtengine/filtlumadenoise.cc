/*
 * filtlumadenoise.cc
 *
 *  Created on: Aug 27, 2010
 *      Author: gabor
 */

#include "filtlumadenoise.h"
#include "rtengine.h"
#include "macros.h"
#include "bilateral2.h"

namespace rtengine {

LumaDenoiseFilterDescriptor lumaDenoiseFilterDescriptor;

LumaDenoiseFilterDescriptor::LumaDenoiseFilterDescriptor ()
	: FilterDescriptor ("LuminanceDenoiser", MultiImage::Lab, MultiImage::Lab) {

    addTriggerEvent (EvLDNEnabled);
    addTriggerEvent (EvLDNRadius);
    addTriggerEvent (EvLDNEdgeTolerance);
}

void LumaDenoiseFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new LumaDenoiseFilter ());
}

LumaDenoiseFilter::LumaDenoiseFilter ()
	: Filter (&lumaDenoiseFilterDescriptor) {
}

void LumaDenoiseFilter::getReqiredBufferSize (int& w, int& h) {

    w = getSourceImageView().getPixelWidth();
    h = getSourceImageView().getPixelHeight();
}

void LumaDenoiseFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

    if (params->lumaDenoise.enabled && sourceImage->width>=8 && sourceImage->height>=8) {
        bilateral<unsigned short, unsigned int> (sourceImage->cieL, targetImage->cieL, (unsigned short**)(buffer->rows), sourceImage->width, sourceImage->height, procParams->lumaDenoise.radius / getScale(), procParams->lumaDenoise.edgetolerance, multiThread);
    }
    else if (sourceImage != targetImage)
        targetImage->copyFrom (sourceImage);
}

}
