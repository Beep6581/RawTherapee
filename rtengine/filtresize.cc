/*
 * filtcoarse.cc
 *
 *  Created on: Sep 17, 2010
 *      Author: gabor
 */

#include "filtresize.h"
#include "rtengine.h"
#include "macros.h"

namespace rtengine {

ResizeFilterDescriptor resizeFilterDescriptor;

ResizeFilterDescriptor::ResizeFilterDescriptor ()
	: FilterDescriptor ("Resize", MultiImage::RGB, MultiImage::RGB, true) {

	addTriggerEvent (EvResizeScale);
    addTriggerEvent (EvResizeMethod);
    addTriggerEvent (EvResizeSpec);
    addTriggerEvent (EvResizeWidth);
    addTriggerEvent (EvResizeHeight);
    addTriggerEvent (EvResizeEnabled);
}

void ResizeFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new ResizeFilter ());
}

ResizeFilter::ResizeFilter ()
	: Filter (&resizeFilterDescriptor) {
}

ImageView ResizeFilter::calculateTargetImageView (const ImageView& requestedImView) {

    return requestedImView;
}

ImageView ResizeFilter::calculateSourceImageView (const ImageView& requestedImView) {

    int x1, y1, x2, y2;
    reverseTransPoint (requestedImView.x, requestedImView.y, x1, y1);
    reverseTransPoint (requestedImView.x + requestedImView.w - 1, requestedImView.y + requestedImView.h - 1, x2, y2);

    return ImageView (std::min(x1,x2), std::min(y1,y2), ABS(x2-x1)+1, ABS(y2-y1)+1, 1);
}

double ResizeFilter::getScale () {

    int s = getTargetImageView().skip;
    if (s > 0)
        return parent->getScale () / getResizeScale ();
    else
        return parent->getScale ();
}


double ResizeFilter::getResizeScale () {

    if (procParams->resize.enabled) {
        int ow, oh;
        getPreviousFilter()->getFullImageSize (ow, oh);
        if (procParams->resize.dataspec==1)
            return procParams->resize.width / ow;
        else if (procParams->resize.dataspec==2)
            return procParams->resize.height / oh;
        else if (procParams->resize.dataspec==0)
            return procParams->resize.scale;
        else
           return 1.0;
    }
    else
        return 1.0;
}

void ResizeFilter::getFullImageSize (int& w, int& h) {

    int ow, oh;
    getPreviousFilter()->getFullImageSize (ow, oh);
    if (procParams->resize.enabled) {
        if (procParams->resize.dataspec==1) {
            w = procParams->resize.width;
            h = oh * w / ow;
        }
        else if (procParams->resize.dataspec==2) {
            h = procParams->resize.height;
            w = ow * h / oh;
        }
        else if (procParams->resize.dataspec==0) {
            h = oh * procParams->resize.scale;
            w = ow * procParams->resize.scale;
        }
        else {
            w = oh;
            h = ow;
        }
    }
    else {
        w = ow;
        h = oh;
    }
}

void ResizeFilter::reverseTransPoint (int x, int y, int& xv, int& yv) {

    xv = x / getResizeScale ();
    yv = y / getResizeScale ();
}

void ResizeFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

}

}
