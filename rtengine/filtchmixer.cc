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

void ColorMixerFilterDescriptor::getDefaultParameters (ProcParams& defProcParams) const {

	IntList r(3), g(3), b(3);
	r[0] = 100; r[1] = 0;   r[2] = 0;
	g[0] = 0;   g[1] = 100; g[2] = 0;
	b[0] = 0;   b[1] = 0;   b[2] = 100;

	defProcParams.setIntegerList ("ChMixerRed",   r);
	defProcParams.setIntegerList ("ChMixerGreen", g);
	defProcParams.setIntegerList ("ChMixerBlue",  b);
}

void ColorMixerFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new ColorMixerFilter ());
}

ColorMixerFilter::ColorMixerFilter ()
	: Filter (&colorMixerFilterDescriptor) {
}

void ColorMixerFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer) {

	IntList& red   = procParams->getIntegerList ("ChMixerRed");
	IntList& green = procParams->getIntegerList ("ChMixerGreen");
	IntList& blue  = procParams->getIntegerList ("ChMixerBlue");

    bool mixchannels = red[0]!=100 || red[1]!=0     || red[2]!=0
                    || green[0]!=0 || green[1]!=100 || green[2]!=0
                    || blue[0]!=0  || blue[1]!=0    || blue[2]!=100;

    if (mixchannels)
        #pragma omp parallel for if (multiThread)
        for (int i=0; i<sourceImage->height; i++)
            for (int j=0; j<sourceImage->width; j++) {
                float r = sourceImage->r[i][j];
                float g = sourceImage->g[i][j];
                float b = sourceImage->b[i][j];
                targetImage->r[i][j] = (r*red[0]   + g*red[1]   + b*red[2])   / 100.0;
                targetImage->g[i][j] = (r*green[0] + g*green[1] + b*green[2]) / 100.0;
                targetImage->b[i][j] = (r*blue[0]  + g*blue[1]  + b*blue[2])  / 100.0;
            }
    else if (sourceImage != targetImage)
        targetImage->copyFrom (sourceImage);
}

}
