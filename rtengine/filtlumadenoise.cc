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

void LumaDenoiseFilterDescriptor::getDefaultParameters (ProcParams& defProcParams) const {

	defProcParams.setFloat   ("LumaDenoise", "Radius", 1.9);
	defProcParams.setFloat   ("LumaDenoise", "EdgeTolerance", 3.0);
	defProcParams.setBoolean ("LumaDenoise", "Enabled", false);
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

void LumaDenoiseFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer) {

	bool enabled = procParams->getBoolean ("LumaDenoise", "Enabled");
    if (getTargetImageView().skip==1 && enabled && sourceImage->width>=8 && sourceImage->height>=8) {
    	float radius = procParams->getFloat ("LumaDenoise", "Radius");
    	float etol = procParams->getFloat   ("LumaDenoise", "EdgeTolerance");

    	bilateral (sourceImage->cieL, targetImage->cieL, buffer->rows, sourceImage->width, sourceImage->height, radius / getScale(), etol, multiThread);

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
