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

    applyOnThumbnail = false;
}

void LumaDenoiseFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new LumaDenoiseFilter ());
}

LumaDenoiseFilter::LumaDenoiseFilter ()
	: Filter (&lumaDenoiseFilterDescriptor) {
}

Dim LumaDenoiseFilter::getReqiredBufferSize () {

    return getScaledTargetImageView().getSize();
}

void LumaDenoiseFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

    if (getTargetImageView().skip==1 && procParams->lumaDenoise.enabled && sourceImage->width>=8 && sourceImage->height>=8) {
        bilateral<unsigned short, unsigned int> (sourceImage->cieL, targetImage->cieL, (unsigned short**)(buffer->rows), sourceImage->width, sourceImage->height, procParams->lumaDenoise.radius / getScale(), procParams->lumaDenoise.edgetolerance, multiThread);
        if (sourceImage != targetImage)
        	for (int i=0; i<targetImage->height; i++)
        		for (int j=0; j<targetImage->width; j++) {
        			targetImage->ciea[i][j] = sourceImage->ciea[i][j];
        			targetImage->cieb[i][j] = sourceImage->cieb[i][j];
        		}
    }
    else if (sourceImage != targetImage)
        targetImage->copyFrom (sourceImage);
}

}
