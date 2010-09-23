#ifndef _MULTIIMAGE_H_
#define _MULTIIMAGE_H_

#include <string>
#include "buffer.h"

namespace rtengine {

class MultiImage {

	static int *cacheL, *cacheab;
    static int *ycache, *xzcache;
	static bool labConversionCacheInitialized;

	unsigned short* data;
	int allocWidth, allocHeight;

	void initLabConversionCache ();
	void initArrays ();

public:
	unsigned rawFilter;
	int width, height;
	enum ColorSpace { Invalid, RGB, Lab, Raw } colorSpace;	// RGB: linear, in "working" color space!

	unsigned short** r;
	unsigned short** g;
	unsigned short** b;
	unsigned short** cieL;		/// stores cieL*655.35
	short** ciea;				/// stores ciea*16384/500
	short** cieb;				/// stores cieb*16384/200
	unsigned short** raw;

	MultiImage (int w, int h, ColorSpace cs = RGB);
	MultiImage (const MultiImage& other);
	~MultiImage ();

	void convertTo (ColorSpace cs, bool multiThread=true, std::string workingColorSpace="sRGB");	// "workingColorSpace": the name of the working color space
	void switchTo  (ColorSpace cs);	// just switches the "colorSpace" without doing any conversion
	bool setDimensions (int w, int h); // sets dimensions without reallocating everything
	bool copyFrom (MultiImage* other);
	bool copyFrom (MultiImage* other, int ofsx, int ofsy, int skip);
    Buffer<unsigned short> getBufferView (unsigned short** channel);
    Buffer<short> getBufferView (short** channel);

	inline bool raw_isRed (int row, int col) {
		return (rawFilter >> (((row << 1 & 14) + (col & 1)) << 1) & 3)==0;
	}
	inline bool raw_isGreen (int row, int col) {
		return (rawFilter >> (((row << 1 & 14) + (col & 1)) << 1) & 3)==1;
	}
	inline bool raw_isBlue (int row, int col) {
		return (rawFilter >> (((row << 1 & 14) + (col & 1)) << 1) & 3)==2;
	}

	int getAllocWidth () { return allocWidth; }
	int getAllocHeight () { return allocHeight; }
	unsigned short* getData () { return data; }
};
}
#endif
