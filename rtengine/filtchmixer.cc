/*
 * filtchmixer.cc
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#include "filtchmixer.h"
#include "rtengine.h"
#include "macros.h"

namespace rtengine {

ColorMixerFilterDescriptor colorMixerFilterDescriptor;

ColorMixerFilterDescriptor::ColorMixerFilterDescriptor ()
	: FilterDescriptor ("ColorMixer", MultiImage::RGB, MultiImage::RGB) {

	addTriggerEvent (EvChMixer);
}

void ColorMixerFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new ColorMixerFilter ());
}

ColorMixerFilter::ColorMixerFilter ()
	: Filter (&colorMixerFilterDescriptor) {
}

void ColorMixerFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer) {

    bool mixchannels = procParams->chmixer.red[0]!=100 || procParams->chmixer.red[1]!=0     || procParams->chmixer.red[2]!=0
                    || procParams->chmixer.green[0]!=0 || procParams->chmixer.green[1]!=100 || procParams->chmixer.green[2]!=0
                    || procParams->chmixer.blue[0]!=0  || procParams->chmixer.blue[1]!=0    || procParams->chmixer.blue[2]!=100;

    if (mixchannels)
        #pragma omp parallel for if (multiThread)
        for (int i=0; i<sourceImage->height; i++)
            for (int j=0; j<sourceImage->width; j++) {
                float r = sourceImage->r[i][j];
                float g = sourceImage->g[i][j];
                float b = sourceImage->b[i][j];
                targetImage->r[i][j] = (r*procParams->chmixer.red[0]   + g*procParams->chmixer.red[1]   + b*procParams->chmixer.red[2]) / 100.0;
                targetImage->g[i][j] = (r*procParams->chmixer.green[0] + g*procParams->chmixer.green[1] + b*procParams->chmixer.green[2]) / 100.0;
                targetImage->b[i][j] = (r*procParams->chmixer.blue[0]  + g*procParams->chmixer.blue[1]  + b*procParams->chmixer.blue[2]) / 100.0;
            }
    else if (sourceImage != targetImage)
        targetImage->copyFrom (sourceImage);
}

}
