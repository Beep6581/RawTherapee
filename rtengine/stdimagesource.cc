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
	: ImageSource (), img (NULL), idata (NULL), autoWBComputed (false) {
}

StdImageSource::~StdImageSource () {

	delete img;
	delete idata;
}

int StdImageSource::load (const Glib::ustring& fileName, ProgressListener* plistener) {

    this->fileName = fileName;

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

    int error = img->load (fileName);
    if (error) {
        delete img;
        img = NULL;
        return error;
    }

    idata = new ImageData (fileName);

    if (plistener) {
        plistener->setProgressStr ("Ready.");
        plistener->setProgress (1.0);
    }

    autoWBComputed = false;

    return 0;
}

Matrix33 StdImageSource::getCamToRGBMatrix () {

	return Matrix33 ();
}

Matrix33 StdImageSource::getRGBToCamMatrix () {

	return Matrix33 ();
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

	float avg_r = 0;
    float avg_g = 0;
    float avg_b = 0;
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

    int v;
    for (int i=1; i<img->height-1; i++)
        for (int j=1; j<img->width-1; j++) {
            if (img->r[i][j]>upper || img->g[i][j]>upper || img->b[i][j]>upper
            		|| img->r[i][j]<lower || img->g[i][j]<lower || img->b[i][j]<lower)
                continue;
            v = 1.0; for (int k=0; k<p; k++) v *= img->r[i][j];
            avg_r += p;
            v = 1.0; for (int k=0; k<p; k++) v *= img->g[i][j];
            avg_g += v;
            v = 1.0; for (int k=0; k<p; k++) v *= img->b[i][j];
            avg_b += v;
            n++;
        }

    autoWB = ColorTemp (pow(avg_r/n, 1.0/p), pow(avg_g/n, 1.0/p), pow(avg_b/n, 1.0/p));
    autoWBComputed = true;

    return autoWB;
}

ColorTemp StdImageSource::getSpotWB (std::vector<Coord2D> red, std::vector<Coord2D> green, std::vector<Coord2D>& blue) {

    int x; int y;
    float reds = 0, greens = 0, blues = 0;
    int rn = 0, gn = 0, bn = 0;
    for (int i=0; i<red.size(); i++) {
        if (red[i].x >= 0 && red[i].y >= 0 && round(red[i].x) < img->width && round(red[i].y) < img->height) {
            reds += img->r[(int)round(red[i].y)][(int)round(red[i].x)];
            rn++;
        }
        if (green[i].x >= 0 && green[i].y >= 0 && round(green[i].x) < img->width &&round(green[i].y) < img->height) {
            greens += img->g[(int)round(green[i].y)][(int)round(green[i].x)];
            gn++;
        }
        if (blue[i].x >= 0 && blue[i].y >= 0 && round(blue[i].x) < img->width && round(blue[i].y) < img->height) {
            blues += img->b[(int)round(blue[i].y)][(int)round(blue[i].x)];
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

Dim StdImageSource::getFullImageSize () {

    return Dim (img->width, img->height);
}

void StdImageSource::getImage (const ImageView& view, MultiImage* targetImage) {

	int x = 0, y = 0;
	for (int i=view.y; i<view.y+view.h; i+=view.skip) {
		x = 0;
		for (int j=view.x; j<view.x+view.w; j+=view.skip) {
			targetImage->r[y][x] = img->r[i][j] / 65535.0;
			targetImage->g[y][x] = img->g[i][j] / 65535.0;
			targetImage->b[y][x] = img->b[i][j] / 65535.0;
			x++;
		}
		y++;
	}
}

}
