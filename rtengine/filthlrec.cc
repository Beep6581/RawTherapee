/*
 * filthlrec.cc
 *
 *  Created on: Sep 23, 2010
 *      Author: gabor
 */

#include "filthlrec.h"
#include "rtengine.h"
#include "macros.h"
#include "iccstore.h"
#include "filterchain.h"
#include "util.h"

namespace rtengine {

HighlightRecoveryFilterDescriptor highlightRecoveryFilterDescriptor;

HighlightRecoveryFilterDescriptor::HighlightRecoveryFilterDescriptor ()
	: FilterDescriptor ("HighlightRecovery", MultiImage::RGB, MultiImage::RGB) {

	addTriggerEvent (EvHREnabled);
	addTriggerEvent (EvHRAmount);
	addTriggerEvent (EvHRMethod);

    applyOnThumbnail = false;
    applyOnStdImage  = false;
}

void HighlightRecoveryFilterDescriptor::getDefaultParameters (ProcParams& defProcParams) const {

	defProcParams.setString  ("HighlightRecoveryMethod", "Luminance");
	defProcParams.setBoolean ("HighlightRecoveryEnabled", true);
}

void HighlightRecoveryFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new HighlightRecoveryFilter ());
}

HighlightRecoveryFilter::HighlightRecoveryFilter ()
	: Filter (&highlightRecoveryFilterDescriptor) {
}

void HighlightRecoveryFilter::luminance (MultiImage* sourceImage, MultiImage* targetImage) {

    #pragma omp parallel for if (multiThread)
    for (int i=0; i<sourceImage->height; i++)
        for (int j=0; j<sourceImage->width; j++) {
            float r = sourceImage->r[i][j], g = sourceImage->g[i][j], b = sourceImage->b[i][j];
            if (r>1.0 || g>1.0 || b>1.0) {
            	float ro = MIN (r, 1.0);
            	float go = MIN (g, 1.0);
            	float bo = MIN (b, 1.0);
            	float L = r + g + b;
            	float C = 1.732050808 * (r - g);
            	float H = 2 * b - r - g;
            	float Co = 1.732050808 * (ro - go);
            	float Ho = 2 * bo - ro - go;
                if (r!=g && g!=b) {
                	float ratio = sqrt ((Co*Co+Ho*Ho) / (C*C+H*H));
                    C *= ratio;
                    H *= ratio;
                }
                float rr = L / 3.0 - H / 6.0 + C / 3.464101615;
                float gr = L / 3.0 - H / 6.0 - C / 3.464101615;
                float br = L / 3.0 + H / 3.0;
                targetImage->r[i][j] = rr;
                targetImage->g[i][j] = gr;
                targetImage->b[i][j] = br;
            }
            else {
                targetImage->r[i][j] = r;
                targetImage->g[i][j] = g;
                targetImage->b[i][j] = b;
            }
        }
}

void HighlightRecoveryFilter::cieblend (MultiImage* sourceImage, MultiImage* targetImage, Matrix33 cam) {

    static bool crTableReady = false;
    static float fv[0x20000];
    if (!crTableReady) {
        for (int ix=0; ix < 0x20000; ix++) {
        	float rx = ix / 65535.0;
            fv[ix] = rx > 0.008856 ? exp(1.0/3 * log(rx)) : 7.787*rx + 16/116.0;
        }
        crTableReady = true;
    }

    cam.multiply (iccStore.workingSpaceMatrix("sRGB"));
    Matrix33 icam = cam.inverse ();

    #pragma omp parallel for if (multiThread)
    for (int i=0; i<sourceImage->height; i++)
        for (int j=0; j<sourceImage->width; j++) {
            float r = sourceImage->r[i][j], g = sourceImage->g[i][j], b = sourceImage->b[i][j];
            if (r>1.0 || g>1.0 || b>1.0) {
            	float ro = MIN (r, 1.0);
            	float go = MIN (g, 1.0);
            	float bo = MIN (b, 1.0);
            	float xx, yy, zz;
                cam.transform (r, g, b, xx, yy, zz);
                float fy = lutInterp<float,0x20000>(fv, 65535.0*yy);
                // compute LCH decompostion of the clipped pixel (only color information, thus C and H will be used)
                float x, y, z;
                cam.transform (ro, go, bo, x, y, z);
                x = lutInterp<float,0x20000>(fv, 65535.0*x);
                y = lutInterp<float,0x20000>(fv, 65535.0*y);
                z = lutInterp<float,0x20000>(fv, 65535.0*z);
                // convert back to rgb
                float fz = fy - y + z;
                float fx = fy + x - y;
                float zr = (fz<=0.206893) ? ((116.0*fz-16.0)/903.3) : (fz * fz * fz);
                float xr = (fx<=0.206893) ? ((116.0*fx-16.0)/903.3) : (fx * fx * fx);
                x = xr;
                y = yy;
                z = zr;
                float rr, gr, br;
                icam.transform (x, y, z, rr, gr, br);
                targetImage->r[i][j] = rr;
                targetImage->g[i][j] = gr;
                targetImage->b[i][j] = br;
            }
            else {
                targetImage->r[i][j] = r;
                targetImage->g[i][j] = g;
                targetImage->b[i][j] = b;
            }
        }
}

void HighlightRecoveryFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer) {

	String method = procParams->getString ("HighlightRecoveryMethod");
	bool enabled  = procParams->getBoolean ("HighlightRecoveryEnabled");

    ImageSource* imgsrc = getFilterChain ()->getImageSource ();

    if (imgsrc->isRaw() && enabled && method == "Luminance")
        luminance (sourceImage, targetImage);
    else if (imgsrc->isRaw() && enabled && method == "CIELab blending")
        cieblend (sourceImage, targetImage, imgsrc->getCamToRGBMatrix());
    else if (sourceImage!=targetImage)
        targetImage->copyFrom (sourceImage);
}

}
