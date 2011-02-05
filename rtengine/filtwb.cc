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

void WhiteBalanceFilterDescriptor::getDefaultParameters (ProcParams& defProcParams) const {

	defProcParams.setString ("WhiteBalance", "Method", "Camera");
	defProcParams.setFloat  ("WhiteBalance", "Temperature", 6504);
	defProcParams.setFloat  ("WhiteBalance", "Green", 1.00102);
}

void WhiteBalanceFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new WhiteBalanceFilter ());
}

WhiteBalanceFilter::WhiteBalanceFilter ()
	: Filter (&whiteBalanceFilterDescriptor) {
}

void WhiteBalanceFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer) {

	String method = procParams->getString ("WhiteBalance", "Method");
	float temp    = procParams->getFloat ("WhiteBalance", "Temperature");
	float green   = procParams->getFloat ("WhiteBalance", "Green");

    ImageSource* imgsrc = getFilterChain ()->getImageSource ();

	// determine which color temperature to use
	ColorTemp wb;
    if (method=="Camera")
        wb = imgsrc->getCamWB ();
    else if (method=="Auto")
        wb = imgsrc->getAutoWB ();
    else
    	wb = ColorTemp (temp, green);

    // write back actual values so gui can read it out
    procParams->setFloat ("WhiteBalance", "Temperature", wb.getTemp ());
    procParams->setFloat ("WhiteBalance", "Green", wb.getGreen ());

    // compute channel multipliers
    float r, g, b, rreq, greq, breq, rcam, gcam, bcam;
    // wb multipliers are available now in rgb space. transform it back to camera space
    wb.getMultipliers (r, g, b);
    imgsrc->getRGBToCamMatrix().transform (r, g, b, rreq, greq, breq);

    // get camera wb and find out how much we have to multiply the channels relative to that
    imgsrc->getCamWB().getMultipliers (r, g, b);
    imgsrc->getRGBToCamMatrix().transform (r, g, b, rcam, gcam, bcam);
    r = rcam / rreq;
    g = gcam / greq;
    b = bcam / breq;

    // ensure that the luminance is not changed (approximately) when applying the channel multipliers
    float mul_lum = 0.299*r + 0.587*g + 0.114*b;
    r /= mul_lum;
    g /= mul_lum;
    b /= mul_lum;

    // multiply each channel with its multiplier
    #pragma omp parallel for if (multiThread)
	for (int i=0; i<sourceImage->height; i++) {
		for (int j=0; j<sourceImage->width; j++) {
			targetImage->r[i][j] = r * sourceImage->r[i][j];
			targetImage->g[i][j] = g * sourceImage->g[i][j];
			targetImage->b[i][j] = b * sourceImage->b[i][j];
		}
	}
}

}
