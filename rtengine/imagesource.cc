/*
 * imagesource.cc
 *
 *  Created on: Aug 19, 2010
 *      Author: gabor
 */

#include "imagesource.h"
#include "multiimage.h"

namespace rtengine {

ImageSource::ImageSource (FilterDescriptor* descr)
	: references (1), Filter (descr, NULL) {
}

void ImageSource::increaseRef () {

	references++;
}

void ImageSource::decreaseRef () {

	references--;
	if (!references)
		delete this;
}
