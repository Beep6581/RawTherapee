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

namespace rtengine {

FilterDescriptor stdImageSourceFilterDescriptor ("ImageSource", MultiImage::Invalid, MultiImage::RGB);

StdImageSource::StdImageSource ()
	: ImageSource (&stdImageSourceFilterDescriptor), img (NULL), idata (NULL), autoWBComputed (false) {
}

virtual StdImageSource::~StdImageSource () {

	delete img;
	delete idata;
}

int StdImageSource::load (const Glib::ustring& fileName, ProgressListener* plistener = NULL) {

    fileName = fname;

    delete img;
    delete idata;
    img = NULL;
    idata = NULL;

    img = new Image16 ();
    if (plistener) {
        plistener->setProgressStr ("Loading...");
        plistener->setProgress (0.0);
        img->setProgressListener (plistener);
    }

    int error = img->load (fname);
    if (error) {
        delete img;
        img = NULL;
        return error;
    }

    idata = new ImageData (fname);

    if (plistener) {
        plistener->setProgressStr ("Ready.");
        plistener->setProgress (1.0);
    }

    autoWBComputed = false;

    return 0;
}

Matrix33 StdImageSource::getCamToRGBMatrix () {

	double mat[3][3] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
	return Matrix (mat);
}

Matrix33 StdImageSource::getRGBToCamMatrix () {

	double mat[3][3] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
	return Matrix (mat);
}

cmsHPROFILE StdImageSource::getEmbeddedProfile () {

	return img->getEmbeddedProfile ();
}

ColorTemp StdImageSource::getCamWB () {

	return ColorTemp (1.0, 1.0, 1.0);
}

ColorTemp StdImageSource::getAutoWB () {

	if (autoWBComputed)
		return autoWB;

	double avg_r = 0;
    double avg_g = 0;
    double avg_b = 0;
    int n = 0;
    int p = 6;

    // detect tonal range
    int minv = 0xffff;
    int maxv = 0;
    for (int i=1; i<img->height-1; i++)
        for (int j=1; j<img->width-1; j++) {
        	maxv = MAX(maxv, img->r[i][j]);
        	maxv = MAX(maxv, img->g[i][j]);
        	maxv = MAX(maxv, img->b[i][j]);
        	minv = MAX(minv, img->r[i][j]);
        	minv = MAX(minv, img->g[i][j]);
        	minv = MAX(minv, img->b[i][j]);
        }
    // adjust to to 2%
    int upper = minv + (maxv+minv) * 0.98;
    int lower = minv + (maxv+minv) * 0.02;

    for (int i=1; i<img->height-1; i++)
        for (int j=1; j<img->width-1; j++) {
            if (img->r[i][j]>upper || img->g[i][j]>upper || img->b[i][j]>upper
            		|| img->r[i][j]<lower || img->g[i][j]<lower || img->b[i][j]<lower)
                continue;
            avg_r += intpow((double)img->r[i][j], p);
            avg_g += intpow((double)img->g[i][j], p);
            avg_b += intpow((double)img->b[i][j], p);
            n++;
        }

    autoWB = ColorTemp (pow(avg_r/n, 1.0/p), pow(avg_g/n, 1.0/p), pow(avg_b/n, 1.0/p));
    autoWBComputed = true;

    return autoWB;
}

ColorTemp StdImageSource::getSpotWB (std::vector<Coord2D> red, std::vector<Coord2D> green, std::vector<Coord2D>& blue) {

    int x; int y;
    double reds = 0, greens = 0, blues = 0;
    int rn = 0, gn = 0, bn = 0;
    for (int i=0; i<red.size(); i++) {
        if (red[i].x >= 0 && red[i].y >= 0 && red[i].x < img->width && red[i].y < img->height) {
            reds += img->r[red[i].y][red[i].x];
            rn++;
        }
        if (green[i].x >= 0 && green[i].y >= 0 && green[i].x < img->width && green[i].y < img->height) {
            greens += img->g[green[i].y][green[i].x];
            gn++;
        }
        transformPixel (blue[i].x, blue[i].y, tran, x, y);
        if (blue[i].x >= 0 && blue[i].y >= 0 && blue[i].x < img->width && blue[i].y < img->height) {
            blues += img->b[blue[i].y][blue[i].x];
            bn++;
        }
    }

    return ColorTemp (reds/rn, greens/gn, blues/bn);
}

void StdImageSource::getAEHistogram (unsigned int* histogram, int& histcompr) {

    histcompr = 3;

    memset (histogram, 0, (65536>>histcompr)*sizeof(int));

    for (int i=0; i<img->height; i++)
        for (int j=0; j<img->width; j++) {
            histogram[CurveFactory::igamma_srgb (img->r[i][j])>>histcompr]++;
            histogram[CurveFactory::igamma_srgb (img->g[i][j])>>histcompr]++;
            histogram[CurveFactory::igamma_srgb (img->b[i][j])>>histcompr]++;
        }
}

void StdImageSource::getFullImageSize (int& w, int& h) {

	w = img->width;
	h = img->height;
}

void StdImageSource::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

	ImageView& view = getTargetImageView ();
	int x = 0, y = 0;
	for (int i=view.y; i<view.y+view.h; i+=view.skip) {
		for (int j=view.x; j<view.x+view.w; j+=view.skip) {
			targetImage->r[y][x] = img->r[i][j];
			targetImage->g[y][x] = img->g[i][j];
			targetImage->b[y][x] = img->b[i][j];
			x++;
		}
		y++;
	}
}

}
