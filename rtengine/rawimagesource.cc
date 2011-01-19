/*
 * rawimagesource.cc
 *
 *  Created on: Aug 22, 2010
 *      Author: gabor
 */

#include "rawimagesource.h"
#include <algorithm>
#include <lcms2.h>
#include "macros.h"
#include "matrix33.h"

namespace rtengine {

RawImageSource::RawImageSource ()
	: ImageSource (), img (NULL), idata (NULL), autoWBComputed (false),
	  embProfile (NULL), border (0) {
}

RawImageSource::~RawImageSource () {

	delete img;
	delete idata;
	if (embProfile)
		cmsCloseProfile (embProfile);
}

int RawImageSource::load (const Glib::ustring& fileName, ProgressListener* plistener) {

    this->fileName = fileName;

    delete img;
    delete idata;
    img = NULL;
    idata = NULL;
    embProfile = NULL;

    img = new RawImage ();
    if (plistener) {
        plistener->setProgressStr ("Loading...");
        plistener->setProgress (0.0);
    }

    int error = img->load (fileName);
    if (error) {
        delete img;
        img = NULL;
        return error;
    }

    idata = new ImageData (fileName, &img->rml);

    if (img->profileData)
        embProfile = cmsOpenProfileFromMem (img->profileData, img->profileLength);

    if (plistener) {
        plistener->setProgressStr ("Ready.");
        plistener->setProgress (1.0);
    }

    autoWBComputed = false;

    return 0;
}

Matrix33 RawImageSource::getCamToRGBMatrix () {

	if (img)
		return img->cam_srgb;
	else
		return Matrix33 ();
}
Matrix33 RawImageSource::getRGBToCamMatrix () {

	if (img)
		return img->srgb_cam;
	else
		return Matrix33 ();
}

cmsHPROFILE RawImageSource::getEmbeddedProfile () {

	return embProfile;
}

ColorTemp RawImageSource::getCamWB () {

	if (img)
		return img->rgbSpaceTemp;
	else
		return ColorTemp (1.0, 1.0, 1.0);
}

ColorTemp RawImageSource::getAutoWB () {

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
        	maxv = MAX(maxv, img->defgain*img->data[i][j]);
        	minv = MIN(minv, img->defgain*img->data[i][j]);
        }
    maxv = CLIP(maxv);

    // adjust to to 2%
    int upper = minv + (maxv+minv) * 0.98;
    int lower = minv + (maxv+minv) * 0.02;

    int rn = 0, gn = 0, bn = 0;

    for (int i=32; i<img->height-32; i++)
        for (int j=32; j<img->width-32; j++) {
            if (!img->filter) {
            	float d = CLIP(img->defgain*img->data[i][3*j]);
                if (d>upper || d < lower)
                    continue;
                avg_r += d*d*d*d*d*d; rn++;
                d = CLIP(img->defgain*img->data[i][3*j+1]);
                if (d>upper || d < lower)
                    continue;
                avg_g += d*d*d*d*d*d; gn++;
                d = CLIP(img->defgain*img->data[i][3*j+2]);
                if (d>upper || d < lower)
                    continue;
                avg_b += d*d*d*d*d*d; bn++;
            }
            else {
            	float d = CLIP(img->defgain*img->data[i][j]);
                if (d>upper || d < lower)
                    continue;
                float dp = d*d*d*d*d*d;
                if (img->isRed (i,j)) {
                    avg_r += dp;
                    rn++;
                }
                else if (img->isGreen (i,j)) {
                    avg_g += dp;
                    gn++;
                }
                else if (img->isBlue (i,j)) {
                    avg_b += dp;
                    bn++;
                }
            }
        }

    //TODO: CAN BE SIMPLIFIED: matrix multiplication can be avoided by using rgbSpaceTemp instead of camSpaceTemp
    float camwb_red, camwb_green, camwb_blue;
    img->camSpaceTemp.getMultipliers (camwb_red, camwb_green, camwb_blue);

    float reds   = pow (avg_r/rn, 1.0/6.0) * camwb_red;
    float greens = pow (avg_g/gn, 1.0/6.0) * camwb_green;
    float blues  = pow (avg_b/bn, 1.0/6.0) * camwb_blue;

    img->cam_srgb.transform (reds, greens, blues);

    autoWB = ColorTemp (pow(avg_r/n, 1.0/p), pow(avg_g/n, 1.0/p), pow(avg_b/n, 1.0/p));
    autoWBComputed = true;

    return autoWB;
}

ColorTemp RawImageSource::getSpotWB (std::vector<Coord2D> red, std::vector<Coord2D> green, std::vector<Coord2D>& blue) {

    int x; int y;
    int d[9][2] = {0,0, -1,-1, -1,0, -1,1, 0,-1, 0,1, 1,-1, 1,0, 1,1};
    float reds = 0, greens = 0, blues = 0;
    int rn = 0, gn = 0, bn = 0;

    if (!img->filter) {
        for (int i=0; i<red.size(); i++) {
            x = red[i].x + border;
            y = red[i].y + border;
            if (x>=0 && y>=0 && x<img->width && y<img->height) {
                reds += img->data[y][3*x];
                rn++;
            }
            x = green[i].x + border;
            y = green[i].y + border;
            if (x>=0 && y>=0 && x<img->width && y<img->height) {
                greens += img->data[y][3*x+1];
                gn++;
            }
            x = blue[i].x + border;
            y = blue[i].y + border;
            if (x>=0 && y>=0 && x<img->width && y<img->height) {
                blues += img->data[y][3*x+2];
                bn++;
            }
        }
    }
    else {
        for (int i=0; i<red.size(); i++) {
            x = red[i].x + border;
            y = red[i].y + border;
            for (int k=0; k<9; k++) {
                int xv = x + d[k][0];
                int yv = y + d[k][1];
                if (img->isRed(yv,xv) && xv>=0 && yv>=0 && xv<img->width && yv<img->height) {
                    reds += img->data[yv][xv];
                    rn++;
                    break;
                }
            }
            x = green[i].x + border;
            y = green[i].y + border;
            for (int k=0; k<9; k++) {
                int xv = x + d[k][0];
                int yv = y + d[k][1];
                if (img->isGreen(yv,xv) && xv>=0 && yv>=0 && xv<img->width && yv<img->height) {
                    greens += img->data[yv][xv];
                    gn++;
                    break;
                }
            }
            x = blue[i].x + border;
            y = blue[i].y + border;
            for (int k=0; k<9; k++) {
                int xv = x + d[k][0];
                int yv = y + d[k][1];
                if (img->isBlue(yv,xv) && xv>=0 && yv>=0 && xv<img->width && yv<img->height) {
                    blues += img->data[yv][xv];
                    bn++;
                    break;
                }
            }
        }
    }

    //TODO: CAN BE SIMPLIFIED: matrix multiplication can be avoided by using rgbSpaceTemp instead of camSpaceTemp
    float camwb_red, camwb_green, camwb_blue;
    img->camSpaceTemp.getMultipliers (camwb_red, camwb_green, camwb_blue);

    reds = reds/rn * camwb_red;
    greens = greens/gn * camwb_green;
    blues = blues/bn * camwb_blue;

    img->cam_srgb.transform (reds, greens, blues);

    return ColorTemp (reds, greens, blues);
}

void RawImageSource::getAEHistogram (unsigned int* histogram, int& histcompr) {

    histcompr = 3;

    memset (histogram, 0, (65536>>histcompr)*sizeof(int));

    for (int i=border; i<img->height-border; i++) {
        int start, end;
        if (img->fujiWidth) {
            int fw = img->fujiWidth;
            start = ABS(fw-i) + border;
            end = MIN(img->height+ img->width-fw-i, fw+i) - border;
        }
        else {
            start = border;
            end = img->width-border;
        }
        if (img->filter)
            for (int j=start; j<end; j++)
                if (img->isGreen(i,j))
                    histogram[img->data[i][j]>>histcompr]+=2;
                else
                    histogram[img->data[i][j]>>histcompr]+=4;
        else
            for (int j=start; j<3*end; j++) {
                    histogram[img->data[i][j+0]>>histcompr]++;
                    histogram[img->data[i][j+1]>>histcompr]++;
                    histogram[img->data[i][j+2]>>histcompr]++;
            }
    }
}

Dim RawImageSource::getFullImageSize () {

    return Dim (img->width - 2*border, img->height - 2*border);
}

// TODO: Fuji/D1X and non-bayer raw files are broken!!!

void RawImageSource::getImage (const ImageView& view, MultiImage* targetImage) {

	int x = 0, y = 0;
	for (int i=view.y; i<view.y+view.h; i+=view.skip) {
	    x = 0;
		for (int j=view.x; j<view.x+view.w; j+=view.skip)
			targetImage->raw[y][x++] = img->defgain * img->data[i+border][j+border] / 65535.0;
		y++;
	}
	targetImage->rawFilter = img->filter;
}

}
