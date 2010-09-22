/*
 * filtwb.cc
 *
 *  Created on: Aug 27, 2010
 *      Author: gabor
 */

#include "filtwb.h"
#include "rtengine.h"
#include "macros.h"
#include "filterchain.h"

namespace rtengine {

WhiteBalanceFilterDescriptor whiteBalanceFilterDescriptor;

WhiteBalanceFilterDescriptor::WhiteBalanceFilterDescriptor ()
	: FilterDescriptor ("WhiteBalance", MultiImage::RGB, MultiImage::RGB) {

	addTriggerEvent (EvWBMethod);
	addTriggerEvent (EvWBTemp);
	addTriggerEvent (EvWBGreen);
}

void WhiteBalanceFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new WhiteBalanceFilter ());
}

WhiteBalanceFilter::WhiteBalanceFilter ()
	: Filter (&whiteBalanceFilterDescriptor) {
}

void WhiteBalanceFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

    ImageSource* imgsrc = getFilterChain ()->getImageSource ();

	// determine which color temperature to use
	ColorTemp wb;
    if (procParams->wb.method=="Camera")
        wb = imgsrc->getCamWB ();
    else if (procParams->wb.method=="Auto")
        wb = imgsrc->getAutoWB ();
    else
    	wb = ColorTemp (procParams->wb.temperature, procParams->wb.green);

    // write back actual values so gui can read it out
    procParams->wb.temperature = wb.getTemp ();
    procParams->wb.green = wb.getGreen ();

    // compute channel multipliers
    double r, g, b, rreq, greq, breq, rcam, gcam, bcam;
    // wb multipliers are available now in rgb space. transform it back to camera space
    wb.getMultipliers (r, g, b);
    imgsrc->getRGBToCamMatrix().transform (r, g, b, rreq, greq, breq);

    // get camera wb and find out how much we have to multiply the channels relative to that
    imgsrc->getCamWB().getMultipliers (r, g, b);
    imgsrc->getRGBToCamMatrix().transform (r, g, b, rreq, greq, breq);
    r = rcam / rreq;
    g = gcam / greq;
    b = bcam / breq;

    // ensure that the luminance is not changed (approximately) when applying the channel multipliers
    double mul_lum = 0.299*r + 0.587*g + 0.114*b;
    r /= mul_lum;
    g /= mul_lum;
    b /= mul_lum;

    // multiply each channel with its multiplier
    #pragma omp parallel for if (multiThread)
	for (int i=0; i<sourceImage->height; i++) {
		for (int j=0; j<sourceImage->width; j++) {
			targetImage->r[i][j] = CLIP (r * sourceImage->r[i][j]);
			targetImage->g[i][j] = CLIP (g * sourceImage->g[i][j]);
			targetImage->b[i][j] = CLIP (b * sourceImage->b[i][j]);
		}
	}
}

}
