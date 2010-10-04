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
            return Dim (2, sdim.height*omp_get_max_threads()); // since we need double buffer and rt gives int buffer
        else
            return Dim (sdim.width*omp_get_max_threads(), 2);
    }
    else
        return sdim;
}

void ColorDenoiseFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

    Dim sdim = getScaledTargetImageView().getSize();
    if (getTargetImageView().skip==1 && procParams->colorDenoise.enabled && sdim.width >= 8 && sdim.height >= 8) {

        double scale = getScale ();

        Buffer<short> ta = targetImage->getBufferView(targetImage->ciea);
        Buffer<short> tb = targetImage->getBufferView(targetImage->cieb);
        Buffer<short> sa = sourceImage->getBufferView(targetImage->ciea);
        Buffer<short> sb = sourceImage->getBufferView(targetImage->cieb);

        gaussHorizontal<short> (&sa, &ta, (double*)buffer->data, procParams->colorDenoise.amount / 10.0 / scale, multiThread);
        gaussHorizontal<short> (&sb, &tb, (double*)buffer->data, procParams->colorDenoise.amount / 10.0 / scale, multiThread);
        gaussVertical<short>   (&ta, &ta, (double*)buffer->data, procParams->colorDenoise.amount / 10.0 / scale, multiThread);
        gaussVertical<short>   (&tb, &tb, (double*)buffer->data, procParams->colorDenoise.amount / 10.0 / scale, multiThread);
    }
    else if (sourceImage != targetImage)
        targetImage->copyFrom (sourceImage);
}

}
