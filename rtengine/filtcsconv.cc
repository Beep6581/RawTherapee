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

void ColorSpaceConvFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer) {

    ImageSource* imgsrc = getFilterChain ()->getImageSource ();

	String workingProfile = procParams->getString ("ColorManagementWorkingProfile");
	String inputProfile   = procParams->getString ("ColorManagementInputProfile");
	bool gammaOnInput     = procParams->getBoolean ("ColorManagementGammaOnInput");

    cmsHPROFILE in  = NULL;
    cmsHPROFILE out = iccStore->workingSpace (workingProfile);

    bool done = false;

    MultiImage* inImg = sourceImage;
    if (inputProfile == "(embedded)")
        in = imgsrc->getEmbeddedProfile ();
    else if (inputProfile != "" && inputProfile != "(camera)" && inputProfile != "(none)") {
        in = iccStore->getProfile (inputProfile);
        if (in && gammaOnInput) {
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
            Matrix33 working = iccStore->workingSpaceInverseMatrix (workingProfile);
            Matrix33 mat = imgsrc->getCamToRGBMatrix();
            mat.multiply (iccStore->workingSpaceMatrix("sRGB"));
            mat.multiply (working);
            #pragma omp parallel for if (multiThread)
            for (int i=0; i<sourceImage->height; i++)
                for (int j=0; j<sourceImage->width; j++)
                    mat.transform (sourceImage->r[i][j], sourceImage->g[i][j], sourceImage->b[i][j], targetImage->r[i][j], targetImage->g[i][j], targetImage->b[i][j]);
            return;
        }
        else if (sourceImage!=targetImage) {
            targetImage->copyFrom (sourceImage);
            return;
        }
        else
			return;
    }

    if (hTransform && (trIn!=in || trOut!=out)) {
        cmsDeleteTransform (hTransform);
        hTransform = NULL;
    }

    if (!hTransform)
        hTransform = cmsCreateTransform (in, (FLOAT_SH(1)|COLORSPACE_SH(PT_RGB)|CHANNELS_SH(3)|BYTES_SH(4)|PLANAR_SH(1)), out, (FLOAT_SH(1)|COLORSPACE_SH(PT_RGB)|CHANNELS_SH(3)|BYTES_SH(4)|PLANAR_SH(1)), Settings::settings->colorimetricIntent, cmsFLAGS_NOCACHE);

    cmsDoTransform (hTransform, inImg->getData(), targetImage->getData(), inImg->getAllocWidth() * inImg->getAllocHeight() / 2);
}

}
