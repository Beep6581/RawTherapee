#include <multiimage.h>
#include <string.h>
#include <iccstore.h>

#undef CLIPTO
#undef XYZ_MAXVAL

#define XYZ_MAXVAL 2*65536-1
#define CLIPTO(a,b,c) ((a)>(b)?((a)<(c)?(a):(c)):(b))

namespace rtengine {

int* MultiImage::cacheL;
int* MultiImage::cacheab;
int* MultiImage::ycache;
int* MultiImage::xzcache;
bool MultiImage::labConversionCacheInitialized = false;

MultiImage::MultiImage (int w, int h, ColorSpace cs = RGB)
	: data(NULL), width(w), height(h), colorSpace(cs), allocWidth(w), allocHeight(h) {

	initArrays ();
	if (!labConversionCacheInitialized)
		initLabConversionCache ();
}

MultiImage::MultiImage (const MultiImage& other)
	: data(NULL), width(other.width), height(other.height), colorSpace(other.colorSpace), allocWidth(w), allocHeight(h) {

	initArrays ();
	if (colorSpace==Raw)
		memcpy (data, other.data, allocWidth*allocHeight*sizeof(unsigned short));
	else if (colorSpace==RGB || colorSpace==Lab)
		memcpy (data, other.data, 3*allocWidth*allocHeight*sizeof(unsigned short));
}

void MultiImage::initArrays () {

	if (colorSpace!=RGB)
		r = g = b = NULL;
	if (colorSpace!=Lab)
		cieL = ciea = cieb = NULL;
	if (colorSpace!=Raw)
		raw = NULL;

	if (colorSpace==Raw) {
		data = new unsigned short*[allocWidth*allocHeight*sizeof(unsigned short)];
		raw = new unsigned short* [allocHeight];
		for (int i=0; i<allocHeight; i++)
			raw[i] = data + i*allocWidth;
	}
	else if (colorSpace==RGB) {
		data = new unsigned short*[allocWidth*allocHeight*sizeof(unsigned short)*3];
		r = new unsigned short* [allocHeight];
		g = new unsigned short* [allocHeight];
		b = new unsigned short* [allocHeight];
		for (int i=0; i<allocHeight; i++) {
			r[i] = data + i*allocWidth;
			g[i] = data + allocWidth*allocHeight + i*allocWidth;
			b[i] = data + 2*allocWidth*allocHeight + i*allocWidth;
		}
	}
	else if (colorSpace==Lab) {
		data = new unsigned short*[allocWidth*allocHeight*sizeof(unsigned short)*3];
		cieL = new unsigned short* [allocHeight];
		ciea = new short* [allocHeight];
		cieb = new short* [allocHeight];
		for (int i=0; i<allocHeight; i++) {
			cieL[i] = data + i*allocWidth;
			ciea[i] = (short*)(data + allocWidth*allocHeight + i*allocWidth);
			cieb[i] = (short*)(data + 2*allocWidth*allocHeight + i*allocWidth);
		}
	}
}

MultiImage::~MultiImage () {

	delete [] data;
	delete [] r;
	delete [] g;
	delete [] b;
	delete [] cieL;
	delete [] ciea;
	delete [] cieb;
	delete [] raw;
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
    else if (allocWidth==other->allocWidth && allocHeight==other->allocHeight) {
        if (colorSpace==Raw)
            memcpy (data, other->data, allocWidth*allocHeight*sizeof(unsigned short));
        else if (colorSpace==RGB || colorSpace==Lab)
            memcpy (data, other->data, 3*allocWidth*allocHeight*sizeof(unsigned short));
        return true;
    }
    else {
        if (colorSpace==Raw)
            for (int i=0; i<height; i++)
                memcpy (raw[i], other->raw[i], width * sizeof(unsigned short));
        else if (colorSpace==RGB)
            for (int i=0; i<height; i++) {
                memcpy (r[i], other->r[i], width * sizeof(unsigned short));
                memcpy (g[i], other->g[i], width * sizeof(unsigned short));
                memcpy (b[i], other->b[i], width * sizeof(unsigned short));
            }
        else if (colorSpace==Lab)
            for (int i=0; i<height; i++) {
                memcpy (cieL[i], other->cieL[i], width * sizeof(unsigned short));
                memcpy (ciea[i], other->ciea[i], width * sizeof(unsigned short));
                memcpy (cieb[i], other->cieb[i], width * sizeof(unsigned short));
            }
        return true;
    }
}

bool MultiImage::copyFrom (MultiImage* other, int ofsx, int ofsy, int skip) {

	if (ofsx<0 || ofsy<0 || ofsx + (width-1)*skip >= other->width || ofsy + (height-1)*skip >= other->height || skip<1)
		return false;

	if (colorSpace!=other->colorSpace) {
		delete [] r;
		delete [] g;
		delete [] b;
		delete [] cieL;
		delete [] ciea;
		delete [] cieb;
		delete [] raw;
		r = g = b = NULL;
		cieL = ciea = cieb = NULL;
		raw = NULL;
		colorSpace = other->colorSpace;
		if (colorSpace==Raw) {
			raw = new unsigned short* [allocHeight];
			for (int i=0; i<allocHeight; i++)
				raw[i] = data + i*allocWidth;
		}
		else if (colorSpace==RGB) {
			r = new unsigned short* [allocHeight];
			g = new unsigned short* [allocHeight];
			b = new unsigned short* [allocHeight];
			for (int i=0; i<allocHeight; i++) {
				r[i] = data + i*allocWidth;
				g[i] = data + allocWidth*allocHeight + i*allocWidth;
				b[i] = data + 2*allocWidth*allocHeight + i*allocWidth;
			}
		}
		else if (colorSpace==Lab) {
			cieL = new unsigned short* [allocHeight];
			ciea = new short* [allocHeight];
			cieb = new short* [allocHeight];
			for (int i=0; i<allocHeight; i++) {
				cieL[i] = data + i*allocWidth;
				ciea[i] = (short*)(data + allocWidth*allocHeight + i*allocWidth);
				cieb[i] = (short*)(data + 2*allocWidth*allocHeight + i*allocWidth);
			}
		}
	}

	if (colorSpace==Raw) {
		for (int i=0; i<height; i++)
			for (int j=0; j<width; j++)
				raw[i][j] = other->raw[ofsy + i*skip][ofsx + j*skip];
		// TODO!!! UPDATE FILTER
	}
	else if (colorSpace==RGB)
		for (int i=0; i<height; i++)
			for (int j=0; j<width; j++) {
				r[i][j] = other->r[ofsy + i*skip][ofsx + j*skip];
				g[i][j] = other->g[ofsy + i*skip][ofsx + j*skip];
				b[i][j] = other->b[ofsy + i*skip][ofsx + j*skip];
		}
	else if (colorSpace==Lab)
		for (int i=0; i<height; i++)
			for (int j=0; j<width; j++) {
				cieL[i][j] = other->cieL[ofsy + i*skip][ofsx + j*skip];
				ciea[i][j] = other->ciea[ofsy + i*skip][ofsx + j*skip];
				cieb[i][j] = other->cieb[ofsy + i*skip][ofsx + j*skip];
		}
	return true;
}


Buffer<unsigned short> getBufferView (unsigned short** channel) {

    return Buffer<unsigned short> (width, height, data, channel);
}

void MultiImage::convertTo (ColorSpace cs, bool multiThread, std::string workingColorSpace) {

	if (colorSpace==cs || colorSpace==Invalid || cs==Invalid)
		return;

	if (colorSpace==Raw) {
		// convert to RGB by simple bilinear demosaicing
		unsigned short* tmpData = new unsigned short [3*allocWidth*allocHeight*sizeof(unsigned short)];
		r = new unsigned short* [allocHeight];
		g = new unsigned short* [allocHeight];
		b = new unsigned short* [allocHeight];
		for (int i=0; i<allocHeight; i++) {
			r[i] = data + i*allocWidth;
			g[i] = data + allocWidth*allocHeight + i*allocWidth;
			b[i] = data + 2*allocWidth*allocHeight + i*allocWidth;
		}
		// bilinear demosaicing of the center of the image
		#pragma omp parallel for if (multiThread)
		for (int i=1; i<allocHeight-1; i++)
			for (int j=1; j<allocWidth-1; j++)
				if (raw_isRed(i,j)) {
					r[i][j] = tmpData[i][j];
					g[i][j] = (tmpData[i-1][j] + tmpData[i+1][j] + tmpData[i][j-1] + tmpData[i][j+1]) >> 2;
					b[i][j] = (tmpData[i-1][j-1] + tmpData[i+1][j-1] + tmpData[i-1][j+1] + tmpData[i+1][j+1]) >> 2;
				}
				else if (raw_isBlue(i,j)) {
					r[i][j] = (tmpData[i-1][j-1] + tmpData[i+1][j-1] + tmpData[i-1][j+1] + tmpData[i+1][j+1]) >> 2;
					g[i][j] = (tmpData[i-1][j] + tmpData[i+1][j] + tmpData[i][j-1] + tmpData[i][j+1]) >> 2;
					b[i][j] = tmpData[i][j];
				}
		#pragma omp parallel for if (multiThread)
		for (int i=1; i<allocHeight-1; i++)
			for (int j=1; j<allocWidth-1; j++)
				if (raw_isGreen(i,j)) {
					r[i][j] = (r[i-1][j] + r[i+1][j] + r[i][j-1] + r[i][j+1]) >> 2;
					g[i][j] = tmpData[i][j];
					b[i][j] = (b[i-1][j] + b[i+1][j] + b[i][j-1] + b[i][j+1]) >> 2;
				}
		// demosaicing borders less efficiently
		#pragma omp parallel for if (multiThread)
		for (int i=0; i<allocHeight; i++)
			for (int j=0; j<allocWidth; j++)
				if (i==0 || j==0 || i==allocHeight-1 || j==allocWidth-1) {
					int r_ = 0, g_ = 0, b_ = 0, rn = 0, gn = 0, bn = 0;
					for (int x=-1; x<=1; x++)
						for (int y=-1; y<=1; y++)
							if (i+x>=0 && j+y>=0 && i+x<allocHeight && j+y<allocWidth) {
								if (raw_isRed(i+x,j+y))
									r_ += tmpData[i+x][j+y], rn++;
								else if (raw_isGreen(i+x,j+y))
									g_ += tmpData[i+x][j+y], gn++;
								else if (raw_isBlue(i+x,j+y))
									b_ += tmpData[i+x][j+y], bn++;
							}
				}
		delete [] data;
		data = tmpData;
		raw = NULL;
		// convert to the desired color space
		convertTo (cs, workingColorSpace);
	}
	else if (cs==Raw) {
		// convert to rgb first
		convertTo (RGB, workingColorSpace);
		// Do mosaicing
		unsigned short* tmpData = new unsigned short [allocWidth*allocHeight*sizeof(unsigned short)];
		#pragma omp parallel for if (multiThread)
		for (int i=0; i<allocHeight; i++)
			for (int j=0; j<allocWidth; j++)
				if (raw_isRed(i,j))
					tmpData[i][j] = r[i][j];
				else if (raw_isGreen(i,j))
					tmpData[i][j] = g[i][j];
				else if (raw_isBlue(i,j))
					tmpData[i][j] = b[i][j];
		delete [] data;
		r = g = b = NULL;
		raw = data = tmpData;
	}
	else if (colorSpace==RGB && cs==Lab) {
		cieL = r, ciea = (short**)g, cieb = (short**)b;
	    TMatrix wprof = iccStore.workingSpaceMatrix (workingColorSpace);
	    // calculate white point tristimulus
	    double xn = wprof[0][0] + wprof[1][0] + wprof[2][0];
	    double yn = wprof[0][1] + wprof[1][1] + wprof[2][1];
	    double zn = wprof[0][2] + wprof[1][2] + wprof[2][2];
	    int toxyz[3][3] = {								// white point: D50 (tristimulus: 0.96422, 1.0, 0.82521)
	        floor(32768.0 * wprof[0][0] / xn),
	        floor(32768.0 * wprof[0][1] / yn),
	        floor(32768.0 * wprof[0][2] / zn),
	        floor(32768.0 * wprof[1][0] / xn),
	        floor(32768.0 * wprof[1][1] / yn),
	        floor(32768.0 * wprof[1][2] / zn),
	        floor(32768.0 * wprof[2][0] / xn),
	        floor(32768.0 * wprof[2][1] / yn),
	        floor(32768.0 * wprof[2][2] / zn)};

		#pragma omp parallel for if (multiThread)
	    for (int i=0; i<height; i++)
        	for (int j=0; j<width; j++) {
        		int r_ = r[i][j], g_ = g[i][j], b_ = b[i][j];
				int x = (toxyz[0][0] * r_ + toxyz[1][0] * g_ + toxyz[2][0] * b_) >> 15;
				int y = (toxyz[0][1] * r_ + toxyz[1][1] * g_ + toxyz[2][1] * b_) >> 15;
				int z = (toxyz[0][2] * r_ + toxyz[1][2] * g_ + toxyz[2][2] * b_) >> 15;

				x = CLIPTO(x,0,2*65536-1);
				y = CLIPTO(y,0,2*65536-1);
				z = CLIPTO(z,0,2*65536-1);

				cieL[i][j] = cacheL[y];
				ciea[i][j] = cacheab[x] - cacheab[y];
				cieb[i][j] = cacheab[y] - cacheab[z];
        	}
        r = g = b = NULL;
	}
	else if (colorSpace==Lab && cs==RGB) {
		r = cieL, g = (unsigned short**)ciea, b = (unsigned short**)cieb;
	    TMatrix wprof = iccStore.workingSpaceMatrix (workingColorSpace);
	    // calculate the white point tristimulus
	    double xn = wprof[0][0] + wprof[1][0] + wprof[2][0];
	    double yn = wprof[0][1] + wprof[1][1] + wprof[2][1];
	    double zn = wprof[0][2] + wprof[1][2] + wprof[2][2];
	    TMatrix iwprof = iccStore.workingSpaceInverseMatrix (workingColorSpace);
	    int torgb[3][3] = {								// white point: D50 (tristimulus: 0.96422, 1.0, 0.82521)
	        floor(32768.0 * iwprof[0][0] * xn),
	        floor(32768.0 * iwprof[0][1] * xn),
	        floor(32768.0 * iwprof[0][2] * xn),
	        floor(32768.0 * iwprof[1][0] * yn),
	        floor(32768.0 * iwprof[1][1] * yn),
	        floor(32768.0 * iwprof[1][2] * yn),
	        floor(32768.0 * iwprof[2][0] * zn),
	        floor(32768.0 * iwprof[2][1] * zn),
	        floor(32768.0 * iwprof[2][2] * zn)};

		#pragma omp parallel for if (multiThread)
		for (int i=0; i<height; i++) {
			for (int j=0; j<width; j++) {
				int L_ = cieL[i][j], y_;

				int x_ = 141558 + (L_ + 10486 + 464 * ciea[i][j] / 100);
				int z_ = 141558 + (L_ + 10486 - 464 * cieb[i][j] / 100);

				x_ = xzcache[x_];
				y_ = ycache[L_];
				z_ = xzcache[z_];

				/* XYZ to RGB */
				int r_ = (torgb[0][0] * x_ + torgb[1][0] * y_ + torgb[2][0] * z_) >> 15;
				int g_ = (torgb[0][1] * x_ + torgb[1][1] * y_ + torgb[2][1] * z_) >> 15;
				int b_ = (torgb[0][2] * x_ + torgb[1][2] * y_ + torgb[2][2] * z_) >> 15;

				r[i][j] = CLIPTO(r_,0,65535);
				g[i][j] = CLIPTO(g_,0,65535);
				b[i][j] = CLIPTO(b_,0,65535);
			}
		}
		cieL = ciea = cieb = NULL;
	}
	colorSpace = cs;
}

void MultiImage::switchTo  (ColorSpace cs) {

	if (colorSpace==cs || colorSpace==Invalid || cs==Invalid || colorSpace==Raw || cs==Raw)
		return;

	if (cs==RGB) {
		r = cieL;
		g = (unsigned short**)ciea;
		b = (unsigned short**)cieb;
		cieL = NULL;
		ciea = NULL;
		cieb = NULL;

	}
	else if (cs==Lab) {
		cieL = r;
		ciea = (short**)g;
		cieb = (short**)b;
		r = NULL;
		g = NULL;
		b = NULL;
	}
	colorSpace = cs;
}

void MultiImage::initLabConversionCache () {

    cacheL = new int[XYZ_MAXVAL+1];
    cacheab = new int[XYZ_MAXVAL+1];

    const int threshold = (int)(0.008856*65535);
    for (int i=0; i<XYZ_MAXVAL; i++)
        if (i>threshold) {
            cacheL[i] = (int)round(655.35 * (116.0 * exp(1.0/3.0 * log(i/65535.0)) - 16.0));
            cacheab[i] = (int)round(16384.0 * exp(1.0/3.0 * log(i/65535.0)));
        }
        else {
            cacheL[i] = (int)round(655.35 * 903.3 * i/65535.0);
            cacheab[i] = (int)round(16384.0 * (7.787*i/65535.0+16.0/116.0));
        }

    double fY;
    ycache = new int[65536];
    for (int i=0; i<65536; i++)
        ycache[i] = (int)round(65535.0 * ((fY=(i/655.35+16.0)/116.0) > 0.206893034 ? fY*fY*fY : 0.001107019*i/655.35));

    double fX;
    xzcache = new int[369623];
    for (int i=-141558; i<228066; i++)
        xzcache[i+141558] = (int)round(65535.0 * ((fX=i/76020.6) > 0.206893034 ? fX*fX*fX : 0.128414183*fX - 0.017712301));

	labConversionCacheInitialized = true;
}

}
