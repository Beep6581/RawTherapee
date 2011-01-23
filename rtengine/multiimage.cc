#include "multiimage.h"
#include <string.h>
#include "iccstore.h"
#include <math.h>
#include "macros.h"
#include "image16.h"
#include "util.h"

#define XYZ2LAB_LUTSIZE 2*65536-1
#define XYZ2LAB_LUTSCALE 65535

#define epsilon 0.00885645 //216/24389
#define kappa 903.2963 //24389/27
#define kappainv 0.00110706 //inverse of kappa
#define kapeps 8 // kappa*epsilon
#define Lab2xyz(f) (( (g=(f)*(f)*(f)) > epsilon) ? (g) : (116*(f)-16)*kappainv)

namespace rtengine {

float* MultiImage::xyz2labCache;
bool MultiImage::labConversionCacheInitialized = false;

MultiImage::MultiImage (int w, int h, ColorSpace cs)
	: data(NULL), width(w), height(h), colorSpace(cs), allocWidth(w), allocHeight(h) {

    ch[0] = ch[1] = ch[2] = NULL;

	initArrays ();
	if (!labConversionCacheInitialized)
		initLabConversionCache ();
}

MultiImage::MultiImage (const MultiImage& other)
	: data(NULL), width(other.width), height(other.height), colorSpace(other.colorSpace), allocWidth(other.allocWidth), allocHeight(other.allocHeight) {

    ch[0] = ch[1] = ch[2] = NULL;

	initArrays ();
	if (colorSpace==Raw)
		memcpy (data, other.data, allocWidth*allocHeight*sizeof(float));
	else if (colorSpace==RGB || colorSpace==Lab)
		memcpy (data, other.data, 3*allocWidth*allocHeight*sizeof(float));
}

void MultiImage::initArrays () {

    delete [] data;
    for (int i=0; i<3; i++)
        delete [] ch[i];

    data = new float[allocWidth*allocHeight*sizeof(float)*3];
    for (int i=0; i<3; i++) {
        ch[i] = new float* [allocHeight];
        for (int j=0; j<allocHeight; j++)
            ch[i][j] = data + i*allocWidth*allocHeight + j*allocWidth;
    }

    r = x = cieL = raw = ch[0];
    g = y = ciea = ch[1];
    b = z = cieb = ch[2];
}

MultiImage::~MultiImage () {

	delete [] data;

	for (int i=0; i<3; i++)
        delete [] ch[i];
}

bool MultiImage::setDimensions (int w, int h) {

	if (w > allocWidth || h > allocHeight)
		return false;
	else {
		width = w;
		height = h;
	}
}

bool MultiImage::copyFrom (MultiImage* other) {

    if (width!=other->width || height!=other->height || colorSpace!=other->colorSpace)
        return false;
    else {
        int channels = colorSpace==Raw ? 1 : 3;

        if (allocWidth==other->allocWidth && allocHeight==other->allocHeight)
            memcpy (data, other->data, channels*allocWidth*allocHeight*sizeof(float));
        else
            for (int k=0; k<channels; k++)
                for (int i=0; i<height; i++)
                    memcpy (ch[k][i], other->ch[k][i], width * sizeof(float));
        return true;
    }
}

bool MultiImage::copyFrom (MultiImage* other, int ofsx, int ofsy, int skip) {

	if (ofsx<0 || ofsy<0 || ofsx + (width-1)*skip >= other->width || ofsy + (height-1)*skip >= other->height || skip<1)
		return false;

	if (colorSpace!=other->colorSpace) {
        colorSpace = other->colorSpace;
	    initArrays ();
	}
    int channels = colorSpace==Raw ? 1 : 3;
    for (int k=0; k<channels; k++)
        for (int i=0; i<height; i++)
            for (int j=0; j<width; j++)
                ch[k][i][j] = other->raw[ofsy + i*skip][ofsx + j*skip];
	return true;
}


Buffer<float> MultiImage::getBufferView (float** channel) {

    return Buffer<float> (width, height, data, channel);
}

void MultiImage::convertTo (ColorSpace cs, bool multiThread, std::string workingColorSpace) {

	if (colorSpace==cs || colorSpace==Invalid || cs==Invalid)
		return;

	if (colorSpace==Raw) {
		// convert to RGB by simple bilinear demosaicing
		float* tmpData = new float[allocWidth*allocHeight*sizeof(float)*3];
		float** tmpch[3];
	    for (int i=0; i<3; i++) {
	        tmpch[i] = new float* [allocHeight];
	        for (int j=0; j<allocHeight; j++)
	            tmpch[i][j] = tmpData + i*allocWidth*allocHeight + j*allocWidth;
	    }
		// bilinear demosaicing of the center of the image
		#pragma omp parallel for if (multiThread)
		for (int i=1; i<allocHeight-1; i++)
			for (int j=1; j<allocWidth-1; j++)
				if (raw_isRed(i,j)) {
				    tmpch[0][i][j] = raw[i][j];
				    tmpch[1][i][j] = (raw[i-1][j] + raw[i+1][j] + raw[i][j-1] + raw[i][j+1]) / 4.0;
				    tmpch[2][i][j] = (raw[i-1][j-1] + raw[i+1][j-1] + raw[i-1][j+1] + raw[i+1][j+1]) / 4.0;
				}
				else if (raw_isBlue(i,j)) {
				    tmpch[0][i][j] = (raw[i-1][j-1] + raw[i+1][j-1] + raw[i-1][j+1] + raw[i+1][j+1]) / 4.0;
				    tmpch[1][i][j] = (raw[i-1][j] + raw[i+1][j] + raw[i][j-1] + raw[i][j+1]) / 4.0;
				    tmpch[2][i][j] = raw[i][j];
				}
		#pragma omp parallel for if (multiThread)
		for (int i=1; i<allocHeight-1; i++)
			for (int j=1; j<allocWidth-1; j++)
				if (raw_isGreen(i,j)) {
				    tmpch[0][i][j] = (tmpch[0][i-1][j] + tmpch[0][i+1][j] + tmpch[0][i][j-1] + tmpch[0][i][j+1]) / 4.0;
					tmpch[1][i][j] = raw[i][j];
					tmpch[2][i][j] = (tmpch[2][i-1][j] + tmpch[2][i+1][j] + tmpch[2][i][j-1] + tmpch[2][i][j+1]) / 4.0;
				}
		// demosaicing borders less efficiently
		#pragma omp parallel for if (multiThread)
		for (int i=0; i<allocHeight; i++)
			for (int j=0; j<allocWidth; j++)
				if (i==0 || j==0 || i==allocHeight-1 || j==allocWidth-1) {
					float r_ = 0, g_ = 0, b_ = 0; int rn = 0, gn = 0, bn = 0;
					for (int x=-1; x<=1; x++)
						for (int y=-1; y<=1; y++)
							if (i+x>=0 && j+y>=0 && i+x<allocHeight && j+y<allocWidth) {
								if (raw_isRed(i+x,j+y))
									r_ += raw[i+x][j+y], rn++;
								else if (raw_isGreen(i+x,j+y))
									g_ += raw[i+x][j+y], gn++;
								else if (raw_isBlue(i+x,j+y))
									b_ += raw[i+x][j+y], bn++;
							}
                    tmpch[0][i][j] = r_ / rn;
                    tmpch[1][i][j] = g_ / rn;
                    tmpch[2][i][j] = b_ / rn;
				}
		delete [] data;
        data = tmpData;
        for (int i=0; i<3; i++) {
            delete [] ch[i];
            ch[i] = tmpch[i];
        }
        r = x = cieL = raw = ch[0];
        g = y = ciea = ch[1];
        b = z = cieb = ch[2];
		// convert to the desired color space
		convertTo (cs, multiThread, workingColorSpace);
	}
	else if (cs==Raw) {
		// convert to rgb first
		convertTo (RGB, multiThread, workingColorSpace);
		// Do mosaicing
		#pragma omp parallel for if (multiThread)
		for (int i=0; i<allocHeight; i++)
			for (int j=0; j<allocWidth; j++)
				if (raw_isRed(i,j))
				    raw[i][j] = r[i][j];
				else if (raw_isGreen(i,j))
				    raw[i][j] = g[i][j];
				else if (raw_isBlue(i,j))
				    raw[i][j] = b[i][j];
	}
	else if (colorSpace==RGB && cs==Lab) {
		convertTo (XYZ, multiThread, workingColorSpace);
		convertTo (Lab, multiThread, workingColorSpace);
	}
	else if (colorSpace==Lab && cs==RGB) {
		convertTo (XYZ, multiThread, workingColorSpace);
		convertTo (RGB, multiThread, workingColorSpace);
	}
    else if (colorSpace==RGB && cs==XYZ) {
	    Matrix33 wprof = iccStore->workingSpaceMatrix (workingColorSpace);
		#pragma omp parallel for if (multiThread)
	    for (int i=0; i<height; i++)
        	for (int j=0; j<width; j++) {
        		float r_ = r[i][j], g_ = g[i][j], b_ = b[i][j];
        		wprof.transform (r_, g_, b_, x[i][j], y[i][j], z[i][j]);
            }
    }
    else if (colorSpace==XYZ && cs==RGB) {
    	Matrix33 iwprof = iccStore->workingSpaceInverseMatrix (workingColorSpace);
        #pragma omp parallel for if (multiThread)
        for (int i=0; i<height; i++) {
            for (int j=0; j<width; j++) {
            	float x_ = x[i][j], y_ = y[i][j], z_ = z[i][j];
        		iwprof.transform (x_, y_, z_, r[i][j], g[i][j], b[i][j]);
            }
        }
    }
    else if (colorSpace==Lab && cs==XYZ) {
	    // calculate white point tristimulus
    	Matrix33 wprof = iccStore->workingSpaceMatrix (workingColorSpace);
	    float xn = wprof.rowsum (0);
	    float yn = wprof.rowsum (1);
	    float zn = wprof.rowsum (2);

	    #pragma omp parallel for if (multiThread)
        for (int i=0; i<height; i++) {
        	float g;
            for (int j=0; j<width; j++) {
            	float fy = (cieL[i][j] + 16.0) / 116.0; // (L+16)/116
            	float fx = fy + ciea[i][j] / 500.0;
            	float fz = fy - cieb[i][j] / 200.0;
				y[i][j] = Lab2xyz(fy) * yn;
            	x[i][j] = Lab2xyz(fx) * xn;
				z[i][j] = Lab2xyz(fz) * zn;
            }
        }
    }
    else if (colorSpace==XYZ && cs==Lab) {
	    // calculate white point tristimulus
    	Matrix33 wprof = iccStore->workingSpaceMatrix (workingColorSpace);
	    float xn = wprof.rowsum (0);
	    float yn = wprof.rowsum (1);
	    float zn = wprof.rowsum (2);

		#pragma omp parallel for if (multiThread)
		for (int i=0; i<height; i++)
			for (int j=0; j<width; j++) {
				float x_ = x[i][j] * XYZ2LAB_LUTSCALE / xn, y_ = y[i][j] * XYZ2LAB_LUTSCALE / yn, z_ = z[i][j] * XYZ2LAB_LUTSCALE / zn;
				float cy = lutInterp<float,XYZ2LAB_LUTSIZE>(xyz2labCache, y_);
				cieL[i][j] = 116.0 * cy - 16.0;
				ciea[i][j] = 500.0 * (lutInterp<float,XYZ2LAB_LUTSIZE>(xyz2labCache, x_) - cy);
				cieb[i][j] = 200.0 * (cy - lutInterp<float,XYZ2LAB_LUTSIZE>(xyz2labCache, z_));
			}
    }

	colorSpace = cs;
}

void MultiImage::switchTo  (ColorSpace cs) {

	if (colorSpace==cs || colorSpace==Invalid || cs==Invalid)
		return;

	colorSpace = cs;
}

void MultiImage::initLabConversionCache () {

    xyz2labCache = new float[XYZ2LAB_LUTSIZE];

    const int threshold = (int)(epsilon*XYZ2LAB_LUTSCALE);
    for (int i=0; i<XYZ2LAB_LUTSIZE; i++)
        if (i>threshold)
        	xyz2labCache[i] = exp(1.0/3.0 * log((double)i/XYZ2LAB_LUTSCALE));
        else
        	xyz2labCache[i] = (kappa * (double)i/XYZ2LAB_LUTSCALE + 16.0) / 116.0;

	labConversionCacheInitialized = true;
}

Image16* MultiImage::createImage () {

    if (colorSpace == RGB || colorSpace == XYZ) {
        Image16* img = new Image16 (width, height);
        for (int i=0; i<height; i++) {
        	for (int j=0; j<width; j++) {
        		img->r[i][j] = CLIP((int)(r[i][j]*65535.0));
    			img->g[i][j] = CLIP((int)(g[i][j]*65535.0));
    			img->b[i][j] = CLIP((int)(b[i][j]*65535.0));
        	}
        }
        return img;
    }
    else if (colorSpace == Lab) {
        Image16* img = new Image16 (width, height);
        for (int i=0; i<height; i++) {
        	for (int j=0; j<width; j++) {
        		img->r[i][j] = CLIP((int)(r[i][j]/100*65535));
    			img->g[i][j] = CLIP((int)(g[i][j]/500*65535.0));
    			img->b[i][j] = CLIP((int)(b[i][j]/200*65535.0));
        	}
        }
        return img;
    }
    else
        return NULL;
}

}
