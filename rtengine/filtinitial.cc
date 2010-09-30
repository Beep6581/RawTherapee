/*
 * filtinitial.cc
 *
 *  Created on: Sep 22, 2010
 *      Author: gabor
 */

#include "filtinitial.h"

namespace rtengine {

RawInitialFilterDescriptor rawInitialFilterDescriptor;
StdInitialFilterDescriptor stdInitialFilterDescriptor;

RawInitialFilterDescriptor::RawInitialFilterDescriptor ()
    : FilterDescriptor ("RawInitial", MultiImage::Invalid, MultiImage::Raw) {
}

StdInitialFilterDescriptor::StdInitialFilterDescriptor ()
    : FilterDescriptor ("StdInitial", MultiImage::Invalid, MultiImage::RGB) {
}

void RawInitialFilterDescriptor::createAndAddToList (Filter* tail) const {}
void StdInitialFilterDescriptor::createAndAddToList (Filter* tail) const {}

InitialFilter::InitialFilter (ImageSource* imgs)
    : Filter (imgs->isRaw() ? (FilterDescriptor*)&rawInitialFilterDescriptor : (FilterDescriptor*)&stdInitialFilterDescriptor), imgsrc (imgs) {}

Dim InitialFilter::getFullImageSize () {

    return imgsrc->getFullImageSize ();
}
void InitialFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

    imgsrc->getImage (getTargetImageView (), targetImage);
}

}
