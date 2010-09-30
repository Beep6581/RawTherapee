/*
 * filtcoarse.cc
 *
 *  Created on: Sep 17, 2010
 *      Author: gabor
 */

#include "filtresize.h"
#include "rtengine.h"
#include "macros.h"
#include "filterchain.h"

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

    double s = getParentFilter()->getScale();

    // if we are processing a thumbnail, do not apply resize filter, just update the "Scale", that is the
    // ratio of the image obtained compared to the image requested
    if (getFilterChain()->getImageSource()->isThumbnail())
        return s*getResizeScale ();
    else
        return s;
}

double ResizeFilter::getResizeScale () {

    if (procParams->resize.enabled) {
        Dim pdim = getPreviousFilter()->getFullImageSize ();
        if (procParams->resize.dataspec==1)
            return procParams->resize.width / pdim.width;
        else if (procParams->resize.dataspec==2)
            return procParams->resize.height / pdim.height;
        else if (procParams->resize.dataspec==0)
            return procParams->resize.scale;
        else
           return 1.0;
    }
    else
        return 1.0;
}

Dim ResizeFilter::getFullImageSize () {

    Dim pdim = getPreviousFilter()->getFullImageSize ();
    if (procParams->resize.enabled) {
        if (procParams->resize.dataspec==1)
            return Dim (procParams->resize.width, pdim.height * procParams->resize.width / pdim.width);
        else if (procParams->resize.dataspec==2)
            return Dim (pdim.width * procParams->resize.height / pdim.height);
        else if (procParams->resize.dataspec==0)
            return Dim (pdim.width * procParams->resize.scale, pdim.height * procParams->resize.scale);
        else
            return pdim;
    }
    else
        return pdim;
}

void ResizeFilter::reverseTransPoint (int x, int y, int& xv, int& yv) {

    xv = x / getResizeScale ();
    yv = y / getResizeScale ();
}

void ResizeFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

}

}
