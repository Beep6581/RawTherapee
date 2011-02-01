/*
 * stdimagesource.cc
 *
 *  Created on: Aug 19, 2010
 *      Author: gabor
 */

#include "stdimagesource.h"
#include "multiimage.h"
#include <algorithm>
#include "macros.h"
#include "matrix33.h"
#include "curves.h"

namespace rtengine {

StdImageSource::StdImageSource ()
	: ImageSource (), img (NULL), idata (NULL), autoWBComputed (false), embProfile (NULL) {
}

StdImageSource::~StdImageSource () {

	delete img;
	delete idata;
    if (embProfile)
        cmsCloseProfile(embProfile);
}

int StdImageSource::load (const String& fileName, ProgressListener* plistener) {

    this->fileName = fileName;

    delete img;
    delete idata;
    if (embProfile)
        cmsCloseProfile(embProfile);
    img = NULL;
    idata = NULL;
    embProfile = NULL;

    img = Image::load (fileName);
    if (!img)
        return 1;

    idata = new ImageData (fileName);

    autoWBComputed = false;

	unsigned char* profileData = NULL;
	int profileLength = 0;
	img->getEmbeddedICCProfile (profileLength, profileData);
	if (profileData)
		embProfile = cmsOpenProfileFromMem (profileData, profileLength);

    return 0;
}

Matrix33 StdImageSource::getCamToRGBMatrix () {

	return Matrix33 ();
}

Matrix33 StdImageSource::getRGBToCamMatrix () {

	return Matrix33 ();
}

cmsHPROFILE StdImageSource::getEmbeddedProfile () {

	return embProfile;
}

ColorTemp StdImageSource::getCamWB () {

	return ColorTemp (1.0, 1.0, 1.0);
}

ColorTemp StdImageSource::getAutoWB () {

	if (autoWBComputed)
		return autoWB;

	float avg_r = 0;
    float avg_g = 0;
    float avg_b = 0;
    int n = 0;
    int p = 6;

	unsigned char* idata = img->getData ();
	int pitch = img->getScanLineSize ();
	int height = img->getHeight (), width = img->getWidth ();

    // detect tonal range
    int minv = 0xffff;
    int maxv = 0;
    for (int i=1; i<height-1; i++) {
		FIRGB16* pixel = (FIRGB16*)(idata + (height-i-1)*pitch);
        for (int j=1; j<width-1; j++) {
        	maxv = MAX(maxv, pixel[j].red);
        	maxv = MAX(maxv, pixel[j].green);
        	maxv = MAX(maxv, pixel[j].blue);
        	minv = MAX(minv, pixel[j].red);
        	minv = MAX(minv, pixel[j].green);
        	minv = MAX(minv, pixel[j].blue);
        }
	}
    // adjust to to 2%
    int upper = minv + (maxv+minv) * 0.98;
    int lower = minv + (maxv+minv) * 0.02;

    int v;
    for (int i=1; i<height-1; i++) {
		FIRGB16* pixel = (FIRGB16*)(idata + (height-i-1)*pitch);
        for (int j=1; j<width-1; j++) {
            if (pixel[j].red>upper || pixel[j].green>upper || pixel[j].blue>upper
            		|| pixel[j].red<lower || pixel[j].green<lower || pixel[j].blue<lower)
                continue;
            v = 1.0; for (int k=0; k<p; k++) v *= pixel[j].red;
            avg_r += p;
            v = 1.0; for (int k=0; k<p; k++) v *= pixel[j].green;
            avg_g += v;
            v = 1.0; for (int k=0; k<p; k++) v *= pixel[j].blue;
            avg_b += v;
            n++;
        }
	}

    autoWB = ColorTemp (pow(avg_r/n, 1.0/p), pow(avg_g/n, 1.0/p), pow(avg_b/n, 1.0/p));
    autoWBComputed = true;

    return autoWB;
}

ColorTemp StdImageSource::getSpotWB (std::vector<Coord2D> red, std::vector<Coord2D> green, std::vector<Coord2D>& blue) {

    int x; int y;
    float reds = 0, greens = 0, blues = 0;
    int rn = 0, gn = 0, bn = 0;
    
	unsigned char* idata = img->getData ();
	int pitch = img->getScanLineSize ();
	int height = img->getHeight (), width = img->getWidth ();

    for (int i=0; i<red.size(); i++) {
        if (red[i].x >= 0 && red[i].y >= 0 && round(red[i].x) < width && round(red[i].y) < height) {
			x = (int)round(red[i].x);
			y = (int)round(red[i].y);
			reds += ((FIRGB16*)(idata + (height-y-1)*pitch))[x].red;
            rn++;
        }
        if (green[i].x >= 0 && green[i].y >= 0 && round(green[i].x) < width &&round(green[i].y) < height) {
			x = (int)round(green[i].x);
			y = (int)round(green[i].y);
			greens += ((FIRGB16*)(idata + (height-y-1)*pitch))[x].green;
            gn++;
        }
        if (blue[i].x >= 0 && blue[i].y >= 0 && round(blue[i].x) < width && round(blue[i].y) < height) {
			x = (int)round(blue[i].x);
			y = (int)round(blue[i].y);
			blues += ((FIRGB16*)(idata + (height-y-1)*pitch))[x].blue;
            bn++;
        }
    }

    return ColorTemp (reds/rn, greens/gn, blues/bn);
}

void StdImageSource::getAEHistogram (unsigned int* histogram, int& histcompr) {

    histcompr = 3;

    memset (histogram, 0, (65536>>histcompr)*sizeof(int));

	unsigned char* idata = img->getData ();
	int pitch = img->getScanLineSize ();
	int height = img->getHeight (), width = img->getWidth ();
    
    for (int i=0; i<height; i++) {
		FIRGB16* pixel = (FIRGB16*)(idata + (height-i-1)*pitch);
        for (int j=0; j<width; j++) {
            histogram[CurveFactory::igamma_srgb (pixel[j].red)>>histcompr]++;
            histogram[CurveFactory::igamma_srgb (pixel[j].green)>>histcompr]++;
            histogram[CurveFactory::igamma_srgb (pixel[j].blue)>>histcompr]++;
        }
	}
}

Dim StdImageSource::getFullImageSize () {

    return Dim (img->getWidth(), img->getHeight());
}

void StdImageSource::getImage (const ImageView& view, MultiImage* targetImage) {

	int x = 0, y = 0;

	unsigned char* idata = img->getData ();
	int pitch = img->getScanLineSize ();
	int height = img->getHeight (), width = img->getWidth ();
    
	for (int i=view.y; i<view.y+view.h; i+=view.skip) {
		x = 0;
		FIRGB16* pixel = (FIRGB16*)(idata + (height-i-1)*pitch);
		for (int j=view.x; j<view.x+view.w; j+=view.skip) {
			targetImage->r[y][x] = pixel[j].red / 65535.0;
			targetImage->g[y][x] = pixel[j].green / 65535.0;
			targetImage->b[y][x] = pixel[j].blue / 65535.0;
			x++;
		}
		y++;
	}
}

}
