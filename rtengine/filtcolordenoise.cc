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

    applyOnThumbnail = false;
}

void ColorDenoiseFilterDescriptor::getDefaultParameters (ProcParams& defProcParams) const {

	defProcParams.setBoolean ("ColorDenoiseEnabled", false);
	defProcParams.setDouble  ("ColorDenoiseAmount",  50);
}

void ColorDenoiseFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new ColorDenoiseFilter ());
}

ColorDenoiseFilter::ColorDenoiseFilter ()
	: Filter (&colorDenoiseFilterDescriptor) {
}

Dim ColorDenoiseFilter::getReqiredBufferSize () {

    Dim sdim = getScaledTargetImageView().getSize();

    if (procParams->colorDenoise.enabled && sdim.width >= 8 && sdim.height >= 8) {
        if (sdim.height > sdim.width)
            return Dim (2, sdim.height*omp_get_max_threads()); // since we need double buffer and rt gives float buffer
        else
            return Dim (sdim.width*omp_get_max_threads(), 2);
    }
    else
        return sdim;
}

void ColorDenoiseFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer) {

    Dim sdim = getScaledTargetImageView().getSize();
    if (getTargetImageView().skip==1 && procParams->getBoolean ("ColorDenoiseEnabled") && sdim.width >= 8 && sdim.height >= 8) {

        double scale = getScale ();
    	double radius = procParams->getDouble  ("ColorDenoiseAmount") / 10.0 / scale;

        Buffer<float> ta = targetImage->getBufferView(targetImage->ciea);
        Buffer<float> tb = targetImage->getBufferView(targetImage->cieb);
        Buffer<float> sa = sourceImage->getBufferView(targetImage->ciea);
        Buffer<float> sb = sourceImage->getBufferView(targetImage->cieb);

        gaussHorizontal<float> (&sa, &ta, sdim, (double*)buffer->data, radius, multiThread);
        gaussHorizontal<float> (&sb, &tb, sdim, (double*)buffer->data, radius, multiThread);
        gaussVertical<float>   (&ta, &ta, sdim, (double*)buffer->data, radius, multiThread);
        gaussVertical<float>   (&tb, &tb, sdim, (double*)buffer->data, radius, multiThread);

        if (sourceImage != targetImage)
        	for (int i=0; i<targetImage->height; i++)
        		for (int j=0; j<targetImage->width; j++)
        			targetImage->cieL[i][j] = sourceImage->cieL[i][j];
    }
    else if (sourceImage != targetImage)
        targetImage->copyFrom (sourceImage);
}

}
