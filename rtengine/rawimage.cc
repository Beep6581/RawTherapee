/*
 * rawimage.cc
 *
 *  Created on: Aug 22, 2010
 *      Author: gabor
 */

#include "rawimage.h"

RawImage::RawImage ()
	: width(-1), height(-1), filter(0), allocation(NULL), data(NULL),
	  profileLength(0), profileData(NULL) {
}

RawImage::~RawImage() {

	delete [] allocation;
	delete [] data;
	delete [] profileData;
}

