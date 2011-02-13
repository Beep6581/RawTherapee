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

	IntVector r, g, b;

	r.push_back (100); r.push_back (0);  r.push_back (0);
	g.push_back (0); g.push_back (100);  g.push_back (0);
	b.push_back (0); b.push_back (0);  b.push_back (100);

	defProcParams.setIntegerVector ("ChMixer", "Red",   r);
	defProcParams.setIntegerVector ("ChMixer", "Green", g);
	defProcParams.setIntegerVector ("ChMixer", "Blue",  b);
}

void ColorMixerFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new ColorMixerFilter ());
}

ColorMixerFilter::ColorMixerFilter ()
	: Filter (&colorMixerFilterDescriptor) {
}

void ColorMixerFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer) {

	IntVector& red   = procParams->getIntegerVector ("ChMixer", "Red");
	IntVector& green = procParams->getIntegerVector ("ChMixer", "Green");
	IntVector& blue  = procParams->getIntegerVector ("ChMixer", "Blue");

    bool mixchannels = red.size()==3 && green.size()==3 && blue.size()==3 &&
    				  (red[0]!=100 || red[1]!=0     || red[2]!=0
                    || green[0]!=0 || green[1]!=100 || green[2]!=0
                    || blue[0]!=0  || blue[1]!=0    || blue[2]!=100);

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
