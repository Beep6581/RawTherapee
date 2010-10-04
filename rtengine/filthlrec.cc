/*
 * filthlrec.cc
 *
 *  Created on: Sep 23, 2010
 *      Author: gabor
 */

#include "filthlrec.h"
#include "rtengine.h"
#include "macros.h"
#include "iccmatrices.h"
#include "filterchain.h"

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

void HighlightRecoveryFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new HighlightRecoveryFilter ());
}

HighlightRecoveryFilter::HighlightRecoveryFilter ()
	: Filter (&highlightRecoveryFilterDescriptor) {
}

void HighlightRecoveryFilter::luminance (MultiImage* sourceImage, MultiImage* targetImage, int maxval) {

    #pragma omp parallel for if (multiThread)
    for (int i=0; i<sourceImage->height; i++)
        for (int j=0; j<sourceImage->width; j++) {
            int r = sourceImage->r[i][j], g = sourceImage->g[i][j], b = sourceImage->b[i][j];
            if (r>maxval || g>maxval || b>maxval) {
                int ro = MIN (r, maxval);
                int go = MIN (g, maxval);
                int bo = MIN (b, maxval);
                double L = r + g + b;
                double C = 1.732050808 * (r - g);
                double H = 2 * b - r - g;
                double Co = 1.732050808 * (ro - go);
                double Ho = 2 * bo - ro - go;
                if (r!=g && g!=b) {
                    double ratio = sqrt ((Co*Co+Ho*Ho) / (C*C+H*H));
                    C *= ratio;
                    H *= ratio;
                }
                int rr = L / 3.0 - H / 6.0 + C / 3.464101615;
                int gr = L / 3.0 - H / 6.0 - C / 3.464101615;
                int br = L / 3.0 + H / 3.0;
                targetImage->r[i][j] = CLIP(rr);
                targetImage->g[i][j] = CLIP(gr);
                targetImage->b[i][j] = CLIP(br);
            }
            else {
                targetImage->r[i][j] = r;
                targetImage->g[i][j] = g;
                targetImage->b[i][j] = b;
            }
        }
}

void HighlightRecoveryFilter::cieblend (MultiImage* sourceImage, MultiImage* targetImage, int maxval, Matrix33 cam) {

    static bool crTableReady = false;
    static double fv[0x10000];
    if (!crTableReady) {
        for (int ix=0; ix < 0x10000; ix++) {
            double rx = ix / 65535.0;
            fv[ix] = rx > 0.008856 ? exp(1.0/3 * log(rx)) : 7.787*rx + 16/116.0;
        }
        crTableReady = true;
    }

    cam.multiply (sRGB_d50);
    Matrix33 icam = cam.inverse ();

    #pragma omp parallel for if (multiThread)
    for (int i=0; i<sourceImage->height; i++)
        for (int j=0; j<sourceImage->width; j++) {
            int r = sourceImage->r[i][j], g = sourceImage->g[i][j], b = sourceImage->b[i][j];
            if (r>maxval || g>maxval || b>maxval) {
                int ro = MIN (r, maxval);
                int go = MIN (g, maxval);
                int bo = MIN (b, maxval);
                double xx, yy, zz;
                cam.transform (r, g, b, xx, yy, zz);
                double fy = fv[CLIP((int)yy)];
                // compute LCH decompostion of the clipped pixel (only color information, thus C and H will be used)
                double x, y, z;
                cam.transform (ro, go, bo, x, y, z);
                x = fv[CLIP((int)x)];
                y = fv[CLIP((int)y)];
                z = fv[CLIP((int)z)];
                // convert back to rgb
                double fz = fy - y + z;
                double fx = fy + x - y;
                double zr = (fz<=0.206893) ? ((116.0*fz-16.0)/903.3) : (fz * fz * fz);
                double xr = (fx<=0.206893) ? ((116.0*fx-16.0)/903.3) : (fx * fx * fx);
                x = xr*65535.0 - 0.5;
                y = yy;
                z = zr*65535.0 - 0.5;
                double rr, gr, br;
                icam.transform (x, y, z, rr, gr, br);
                targetImage->r[i][j] = CLIP(rr);
                targetImage->g[i][j] = CLIP(gr);
                targetImage->b[i][j] = CLIP(br);
            }
            else {
                targetImage->r[i][j] = r;
                targetImage->g[i][j] = g;
                targetImage->b[i][j] = b;
            }
        }
}

void HighlightRecoveryFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

    ImageSource* imgsrc = getFilterChain ()->getImageSource ();

    if (imgsrc->isRaw() && procParams->hlrecovery.enabled && procParams->hlrecovery.method == "Luminance")
        luminance (sourceImage, targetImage, 65535.0 / imgsrc->getDefGain());
    else if (imgsrc->isRaw() && procParams->hlrecovery.enabled && procParams->hlrecovery.method == "CIELab blending")
        cieblend (sourceImage, targetImage, 65535.0 / imgsrc->getDefGain(), imgsrc->getCamToRGBMatrix());
    else if (sourceImage!=targetImage)
        targetImage->copyFrom (sourceImage);
}

}
