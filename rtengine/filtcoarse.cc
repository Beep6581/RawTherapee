/*
 * filtcoarse.cc
 *
 *  Created on: Sep 17, 2010
 *      Author: gabor
 */

#include "filtcoarse.h"
#include "rtengine.h"
#include "macros.h"

namespace rtengine {

CoarseTransformFilterDescriptor coarseTransformFilterDescriptor;

CoarseTransformFilterDescriptor::CoarseTransformFilterDescriptor ()
	: FilterDescriptor ("CoarseTransform", MultiImage::RGB, MultiImage::RGB) {

	addTriggerEvent (EvCTRotate);
    addTriggerEvent (EvCTHFlip);
    addTriggerEvent (EvCTVFlip);
}

void CoarseTransformFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new CoarseTransformFilter ());
}

CoarseTransformFilter::CoarseTransformFilter ()
	: Filter (&coarseTransformFilterDescriptor) {
}

ImageView CoarseTransformFilter::calculateTargetImageView (const ImageView& requestedImView) {

    return requestedImView;
}

ImageView CoarseTransformFilter::calculateSourceImageView (const ImageView& requestedImView) {

    int x1, y1, x2, y2;
    reverseTransPoint (requestedImView.x, requestedImView.y, x1, y1);
    reverseTransPoint (requestedImView.x + requestedImView.w - 1, requestedImView.y + requestedImView.h - 1, x2, y2);

    return ImageView (std::min(x1,x2), std::min(y1,y2), ABS(x2-x1)+1, ABS(y2-y1)+1, 1);
}

Dim CoarseTransformFilter::getFullImageSize () {

    int ow, oh;
    Dim prevd = getPreviousFilter()->getFullImageSize ();
    if (procParams->coarse.rotate==90 || procParams->coarse.rotate==270)
        return Dim (prevd.height, prevd.width);
    else
        return prevd;
}

Dim CoarseTransformFilter::getReqiredBufferSize () {

    return getScaledTargetImageView().getSize();
}

void CoarseTransformFilter::reverseTransPoint (int x, int y, int& xv, int& yv) {

    Dim pfs = getPreviousFilter()->getFullImageSize ();
    int ow = pfs.width, oh = pfs.height;

    int nx = x, ny = y;
    if (procParams->coarse.vflip)
        ny = oh - 1 - ny;
    if (procParams->coarse.hflip)
        nx = ow - 1 - nx;
    if (procParams->coarse.rotate == 90) {
        xv = ny;
        yv = oh - nx;
    }
    else if (procParams->coarse.rotate == 180) {
        xv = ow - nx;
        yv = oh - ny;
    }
    else if (procParams->coarse.rotate == 270) {
        xv = ow - ny;
        yv = nx;
    }
    else {
        xv = nx;
        yv = ny;
    }
}

void CoarseTransformFilter::hflip (MultiImage* image) {

    int width  = image->width;
    int height = image->height;

    #pragma omp parallel for if (multiThread)
    for (int i=0; i<height; i++) {
        unsigned short v;
        for (int j=0; j<width; j++) {
            v = image->r[i][width-1-j];
            image->r[i][width-1-j] = image->r[i][j];
            image->r[i][j] = v;
            v = image->g[i][width-1-j];
            image->g[i][width-1-j] = image->g[i][j];
            image->g[i][j] = v;
            v = image->b[i][width-1-j];
            image->b[i][width-1-j] = image->b[i][j];
            image->b[i][j] = v;
      }
    }
}

void CoarseTransformFilter::vflip (MultiImage* image) {

    int width  = image->width;
    int height = image->height;

    #pragma omp parallel for if (multiThread)
    for (int i=0; i<height/2; i++) {
        unsigned short v;
        for (int j=0; j<width; j++) {
            v = image->r[i][j];
            image->r[i][j] = image->r[height-1-i][j];
            image->r[height-1-i][j] = v;
            v = image->g[i][j];
            image->g[i][j] = image->g[height-1-i][j];
            image->g[height-1-i][j] = v;
            v = image->b[i][j];
            image->b[i][j] = image->b[height-1-i][j];
            image->b[height-1-i][j] = v;
        }
    }
}

void CoarseTransformFilter::rotate90 (unsigned short** si, unsigned short** ti, int sW, int sH, Buffer<int>* buffer) {

    #pragma omp parallel for if (multiThread)
    for (int i=0; i<sW; i++)
        for (int j=0; j<sH; j++)
            buffer->rows[i][j] = si[sH-1-j][i];

    #pragma omp parallel for if (multiThread)
    for (int i=0; i<sW; i++)
        for (int j=0; j<sH; j++)
            ti[i][j] = buffer->rows[i][j];
}

void CoarseTransformFilter::rotate180 (unsigned short** si, unsigned short** ti, int sW, int sH, Buffer<int>* buffer) {

    #pragma omp parallel for if (multiThread)
    for (int i=0; i<sH; i++)
        for (int j=0; j<sW; j++)
            buffer->rows[i][j] = si[sH-1-i][sW-1-j];

    #pragma omp parallel for if (multiThread)
    for (int i=0; i<sH; i++)
        for (int j=0; j<sW; j++)
            ti[i][j] = buffer->rows[i][j];
}

void CoarseTransformFilter::rotate270 (unsigned short** si, unsigned short** ti, int sW, int sH, Buffer<int>* buffer) {

    #pragma omp parallel for if (multiThread)
    for (int i=0; i<sW; i++)
        for (int j=0; j<sH; j++)
            buffer->rows[i][j] = si[j][sW-1-i];

    #pragma omp parallel for if (multiThread)
    for (int i=0; i<sW; i++)
        for (int j=0; j<sH; j++)
            ti[i][j] = buffer->rows[i][j];
}


void CoarseTransformFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

    if (procParams->coarse.rotate==90) {
        rotate90 (sourceImage->r, targetImage->r, sourceImage->width, sourceImage->height, buffer);
        rotate90 (sourceImage->g, targetImage->g, sourceImage->width, sourceImage->height, buffer);
        rotate90 (sourceImage->b, targetImage->b, sourceImage->width, sourceImage->height, buffer);
    }
    else if (procParams->coarse.rotate==180) {
        rotate180 (sourceImage->r, targetImage->r, sourceImage->width, sourceImage->height, buffer);
        rotate180 (sourceImage->g, targetImage->g, sourceImage->width, sourceImage->height, buffer);
        rotate180 (sourceImage->b, targetImage->b, sourceImage->width, sourceImage->height, buffer);
    }
    else if (procParams->coarse.rotate==270) {
        rotate270 (sourceImage->r, targetImage->r, sourceImage->width, sourceImage->height, buffer);
        rotate270 (sourceImage->g, targetImage->g, sourceImage->width, sourceImage->height, buffer);
        rotate270 (sourceImage->b, targetImage->b, sourceImage->width, sourceImage->height, buffer);
    }
    else if (targetImage!=sourceImage)
        targetImage->copyFrom (sourceImage);

    if (procParams->coarse.hflip)
        hflip (targetImage);
    if (procParams->coarse.vflip)
        vflip (targetImage);
}

}
