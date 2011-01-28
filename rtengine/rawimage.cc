/*
 * rawimage.cc
 *
 *  Created on: Aug 22, 2010
 *      Author: gabor
 */

#include "rawimage.h"
#include <libraw/libraw.h>
#include "macros.h"

#define FC(row,col) \
	(filter >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)

namespace rtengine {

RawImage::RawImage ()
	: width(-1), height(-1), filter(0), allocation(NULL), data(NULL),
	  profileLength(0), profileData(NULL), fujiWidth(0) {
}

RawImage::~RawImage() {

	delete [] allocation;
	delete [] data;
	delete [] profileData;
}

int RawImage::load (const Glib::ustring& fname) {

	LibRaw rawproc;
	
	// open raw file
    int res = rawproc.open_file (fname.c_str());
    
    if (res)
		return res;
		
	// unpack it 
	res = rawproc.unpack ();

	if (res)
		return res;
		
	// read out important parameters
	width = rawproc.imgdata.sizes.width;
	height = rawproc.imgdata.sizes.height;
	filter = rawproc.imgdata.idata.filters;
	
	// cam <-> srgb matrices
	cam_srgb = Matrix33 (rawproc.imgdata.color.rgb_cam);
    srgb_cam = cam_srgb.inverse ();
	
	// determine camera white balance
	float red_multiplier = rawproc.imgdata.color.pre_mul[0];
	float green_multiplier = rawproc.imgdata.color.pre_mul[1];
	float blue_multiplier = rawproc.imgdata.color.pre_mul[2];
	float pre_mul[4];
	int c;
	for (c=0; c<4; c++)
		pre_mul[c] = rawproc.imgdata.color.pre_mul[c];
	if (rawproc.imgdata.color.cam_mul[0] != -1) {
		unsigned sum[8]; int val;
		memset (sum, 0, sizeof(sum));
		for (int row=0; row < 8; row++)
			for (int col=0; col < 8; col++) {
				c = FC(row,col);
				if ((val = rawproc.imgdata.color.white[row][col] - rawproc.imgdata.color.cblack[c]) > 0)
					sum[c] += val;
				sum[c+4]++;
			}
		if (sum[0] && sum[1] && sum[2] && sum[3])
			for (c=0; c<4; c++)
				pre_mul[c] = (float) sum[c+4] / sum[c];
		else if (rawproc.imgdata.color.cam_mul[0] && rawproc.imgdata.color.cam_mul[2])
			memcpy (pre_mul, rawproc.imgdata.color.cam_mul, sizeof(pre_mul));
	}
	
	// fix zeros
	bool zero_is_bad = false;
	zero_is_bad = zero_is_bad || !strcmp(rawproc.imgdata.idata.make,"LEICA") || !strcmp(rawproc.imgdata.idata.make,"Panasonic");
	zero_is_bad = zero_is_bad ||  (!strcmp(rawproc.imgdata.idata.model,"PowerShot SD300") || !strcmp(rawproc.imgdata.idata.model,"PowerShot A460") || !strcmp(rawproc.imgdata.idata.model,"PowerShot A530")
		|| !strcmp(rawproc.imgdata.idata.model,"PowerShot A610") || !strcmp(rawproc.imgdata.idata.model,"PowerShot A620") || !strcmp(rawproc.imgdata.idata.model,"PowerShot A470")
		|| !strcmp(rawproc.imgdata.idata.model,"PowerShot A720") || !strcmp(rawproc.imgdata.idata.model,"PowerShot A630") || !strcmp(rawproc.imgdata.idata.model,"PowerShot A640")
		|| !strcmp(rawproc.imgdata.idata.model,"PowerShot A650") || !strcmp(rawproc.imgdata.idata.model,"PowerShot S3 IS"))  && width > 1600;
	zero_is_bad = zero_is_bad || !strcmp(rawproc.imgdata.idata.model,"PowerShot SX110 IS") || !strcmp(rawproc.imgdata.idata.model,"PowerShot SX20 IS");

	if (zero_is_bad) {
		for (int row=0; row < height; row++)
			for (int col=0; col < width; col++)
				if (rawproc.imgdata.image[row*rawproc.imgdata.sizes.iwidth+col][FC(row,col)] == 0) {
					int tot = 0, n = 0;
					for (int r = row-2; r <= row+2; r++)
						for (int c = col-2; c <= col+2; c++)
							if (r < height && c < width && FC(r,c) == FC(row,col) && rawproc.imgdata.image[r*rawproc.imgdata.sizes.iwidth+c][FC(r,c)])
								tot += (n++,rawproc.imgdata.image[r*rawproc.imgdata.sizes.iwidth+c][FC(r,c)]);
				if (n) 
					rawproc.imgdata.image[row*rawproc.imgdata.sizes.iwidth+col][FC(row,col)] = tot/n;
		  }
	}
	
	// scale colors
	if (pre_mul[3] == 0) 
		pre_mul[3] = pre_mul[1];
	unsigned max = rawproc.imgdata.color.maximum - rawproc.imgdata.color.black;
	double dmin, dmax;
	for (dmin=DBL_MAX, dmax=c=0; c < 4; c++) {
		if (dmin > pre_mul[c])
		dmin = pre_mul[c];
		if (dmax < pre_mul[c])
		dmax = pre_mul[c];
	}
	float scale_mul[4];
	for (c=0; c<4; c++) {
		pre_mul[c] /= dmax;
		scale_mul[c] = pre_mul[c] * 65535.0 / max;
	}
	int size = rawproc.imgdata.sizes.iheight * rawproc.imgdata.sizes.iwidth;
	for (int i=0; i < size*4; i++) {
		int val = rawproc.imgdata.image[0][i];
		val -= rawproc.imgdata.color.cblack[i & 3] + rawproc.imgdata.color.black;
		val *= scale_mul[i & 3];
		rawproc.imgdata.image[0][i] = CLIPTO(val,0,65535);
	}
	float camwb_red   = red_multiplier / pre_mul[0];
	float camwb_green = green_multiplier / pre_mul[1];
	float camwb_blue  = blue_multiplier / pre_mul[2];

	// calculate color temperatures
	camSpaceTemp = ColorTemp (camwb_red, camwb_green, camwb_blue);
    float cam_r, cam_g, cam_b;
    cam_srgb.transform (camwb_red, camwb_green, camwb_blue, cam_r, cam_g, cam_b);
    rgbSpaceTemp = ColorTemp (cam_r, cam_g, cam_b);

	// calculate defgain
	defgain = 1.0 / MIN(MIN(pre_mul[0],pre_mul[1]),pre_mul[2]);

	// pre-interpolate
	for (int row = FC(1,0) >> 1; row < height; row+=2)
		for (int col = FC(row,1) & 1; col < width; col+=2)
			rawproc.imgdata.image[row*width+col][1] = rawproc.imgdata.image[row*width+col][3];
	filter &= ~((filter & 0x55555555) << 1);

	// load image
	allocation = new unsigned short [height*width];
	data = new unsigned short* [height];
	for (int i = 0; i < height; i++)
		data[i] = allocation + i*width;
	for (int row = 0; row < height; row++)
		for (int col = 0; col < width; col++)
			data[row][col] = rawproc.imgdata.image[row*width+col][FC(row,col)];

	// copy embedded icc profile
	profileLength = 0;
	profileData = NULL;
	if (rawproc.imgdata.color.profile_length && rawproc.imgdata.color.profile) {
		profileLength = rawproc.imgdata.color.profile_length;
		profileData = new char [profileLength];
		memcpy (profileData, rawproc.imgdata.color.profile, profileLength);
	}
	return 0;
}

}
