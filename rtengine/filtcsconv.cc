/*
 * filtcsconv.cc
 *
 *  Created on: Sep 23, 2010
 *      Author: gabor
 */

#include "filtcsconv.h"
#include "rtengine.h"
#include "macros.h"
#include "filterchain.h"
#include "iccstore.h"
#include "iccmatrices.h"
#include "settings.h"
#include "curves.h"

namespace rtengine {

ColorSpaceConvFilterDescriptor colorSpaceConvFilterDescriptor;

ColorSpaceConvFilterDescriptor::ColorSpaceConvFilterDescriptor ()
	: FilterDescriptor ("ColorSpaceConversion", MultiImage::RGB, MultiImage::RGB) {

	addTriggerEvent (EvIProfile);
	addTriggerEvent (EvWProfile);
}

void ColorSpaceConvFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new ColorSpaceConvFilter ());
}

ColorSpaceConvFilter::ColorSpaceConvFilter ()
	: Filter (&colorSpaceConvFilterDescriptor), hTransform (NULL) {
}

ColorSpaceConvFilter::~ColorSpaceConvFilter () {

    if (hTransform)
        cmsDeleteTransform (hTransform);
}

void ColorSpaceConvFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

    ImageSource* imgsrc = getFilterChain ()->getImageSource ();

    cmsHPROFILE in  = NULL;
    cmsHPROFILE out = iccStore.workingSpace (procParams->icm.working);;

    bool done = false;

    MultiImage* inImg = sourceImage;

    if (procParams->icm.input == "(embedded)")
        in = imgsrc->getEmbeddedProfile ();
    else if (procParams->icm.input != "" && procParams->icm.input != "(camera)" && procParams->icm.input != "(none)") {
        in = iccStore.getProfile (procParams->icm.input);
        if (in && procParams->icm.gammaOnInput) {
            for (int i=0; i<sourceImage->height; i++)
                for (int j=0; j<sourceImage->width; j++) {
                    targetImage->r[i][j] = CurveFactory::gamma (sourceImage->r[i][j]);
                    targetImage->g[i][j] = CurveFactory::gamma (sourceImage->g[i][j]);
                    targetImage->b[i][j] = CurveFactory::gamma (sourceImage->b[i][j]);
                }
            inImg = targetImage;
        }
    }
    if (!in && !imgsrc->isRaw())
        in = imgsrc->getEmbeddedProfile ();

    if (!in) {
        if (imgsrc->isRaw()) {
            // do the color transform "by hand" to avoid calling slow lcms2
            TMatrix working = iccStore.workingSpaceInverseMatrix (procParams->icm.working);
            Matrix33 mat = imgsrc->getCamToRGBMatrix();
            mat.multiply (sRGB_d50);
            mat.multiply (working);
            #pragma omp parallel for if (multiThread)
            for (int i=0; i<sourceImage->height; i++)
                for (int j=0; j<sourceImage->width; j++) {
                    double newr, newg, newb;
                    mat.transform (sourceImage->r[i][j], sourceImage->g[i][j], sourceImage->b[i][j], newr, newg, newb);
                    targetImage->r[i][j] = CLIP(newr);
                    targetImage->g[i][j] = CLIP(newg);
                    targetImage->b[i][j] = CLIP(newb);
                }
            return;
        }
        else if (sourceImage!=targetImage) {
            targetImage->copyFrom (sourceImage);
            return;
        }
    }

    if (hTransform && (trIn!=in || trOut!=out)) {
        cmsDeleteTransform (hTransform);
        hTransform = NULL;
    }

    if (!hTransform)
        hTransform = cmsCreateTransform (in, TYPE_RGB_16_PLANAR, out, TYPE_RGB_16_PLANAR, Settings::settings->colorimetricIntent, 0);

    cmsDoTransform (hTransform, inImg->getData(), targetImage->getData(), inImg->getAllocWidth() * inImg->getAllocHeight() / 2);
}

}
