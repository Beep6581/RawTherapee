#ifndef _MULTIIMAGE_H_
#define _MULTIIMAGE_H_

#include <string>
#include "buffer.h"
#include "image.h"

namespace rtengine {

class MultiImage {

	static float* xyz2labCache;
	static bool labConversionCacheInitialized;

	float* data;
	int allocWidth, allocHeight;

	void initLabConversionCache ();
	void initArrays ();

	float **ch[3];

public:
	unsigned rawFilter;
	int width, height;
	enum ColorSpace { Invalid, RGB, Lab, XYZ, Raw } colorSpace;	// RGB: linear, in "working" color space!

	float** r;
	float** g;
	float** b;
	float** cieL;		/// stores cieL*655.35
	float** ciea;		/// stores ciea*16384/500
	float** cieb;		/// stores cieb*16384/200
	float** x;
	float** y;
	float** z;
	float** raw;

	MultiImage (int w, int h, ColorSpace cs = RGB);
	MultiImage (const MultiImage& other);
	~MultiImage ();

	void convertTo (ColorSpace cs, bool multiThread=true, std::string workingColorSpace="sRGB");	// "workingColorSpace": the name of the working color space
	void switchTo  (ColorSpace cs);	// just switches the "colorSpace" without doing any conversion
	bool setDimensions (int w, int h); // sets dimensions without reallocating everything
	bool copyFrom (MultiImage* other);
	bool copyFrom (MultiImage* other, int ofsx, int ofsy, int skip);
    Buffer<float> getBufferView (float** channel);

	inline bool raw_isRed (int row, int col) {
		return (rawFilter >> (((row << 1 & 14) + (col & 1)) << 1) & 3)==0;
	}
	inline bool raw_isGreen (int row, int col) {
		return (rawFilter >> (((row << 1 & 14) + (col & 1)) << 1) & 3)==1;
	}
	inline bool raw_isBlue (int row, int col) {
		return (rawFilter >> (((row << 1 & 14) + (col & 1)) << 1) & 3)==2;
	}

	int getAllocWidth ()	{ return allocWidth; }
	int getAllocHeight () 	{ return allocHeight; }
	float* getData () 		{ return data; }

	Image* createImage ();
};
}
#endif
