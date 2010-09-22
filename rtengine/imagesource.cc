/*
 * imagesource.cc
 *
 *  Created on: Aug 19, 2010
 *      Author: gabor
 */

#include "imagesource.h"

namespace rtengine {

ImageSource::ImageSource ()
	: references (1) {
}

void ImageSource::increaseRef () {

	references++;
}

void ImageSource::decreaseRef () {

	references--;
	if (!references)
		delete this;
}
}
