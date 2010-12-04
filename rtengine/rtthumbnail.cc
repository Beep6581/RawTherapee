/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <rtengine.h>
#include <rtthumbnail.h>
#include <image8.h>
#include <lcms.h>
#include <curves.h>
#include <glibmm.h>
#include <improcfun.h>
#include <colortemp.h> 
#include <mytime.h>
#include <utils.h>
#include <iccstore.h>
#include <iccmatrices.h>
#include <rawimagesource.h>
#include <stdimagesource.h>
#include <glib/gstdio.h>
#include <setjmp.h>
#include <safekeyfile.h>
#include <rawimage.h>

extern "C" {
#include <jpeglib.h>
extern jmp_buf jpeg_jmp_buf;
extern GLOBAL(struct jpeg_error_mgr *)
my_jpeg_std_error (struct jpeg_error_mgr * err);
extern GLOBAL(void)
my_jpeg_stdio_src (j_decompress_ptr cinfo, FILE * infile);
}

#define MAXVAL  0xffff
#define CLIP(a) ((a)>0?((a)<MAXVAL?(a):MAXVAL):0)

namespace rtengine {

Thumbnail* Thumbnail::loadFromImage (const Glib::ustring& fname, int &w, int &h, int fixwh) {

    Image16* img = new Image16 ();
    int err = img->load (fname);
    if (err) {
        delete img;
        return NULL;
    }
    
    Thumbnail* tpp = new Thumbnail ();

    tpp->camwbRed = 1.0;
    tpp->camwbGreen = 1.0;
    tpp->camwbBlue = 1.0;

    tpp->embProfileLength = 0;
    unsigned char* data;
    img->getEmbeddedProfileData (tpp->embProfileLength, data);
    if (data && tpp->embProfileLength) {
        tpp->embProfileData = new unsigned char [tpp->embProfileLength];
        memcpy (tpp->embProfileData, data, tpp->embProfileLength);
    }
    else {
        tpp->embProfileLength = 0;
        tpp->embProfileData = NULL;
    }
    
    tpp->redMultiplier = 1.0;
    tpp->greenMultiplier = 1.0;
    tpp->blueMultiplier = 1.0;

    tpp->scaleForSave = 8192;
    tpp->defGain = 1.0;
    tpp->gammaCorrected = false;
    tpp->isRaw = 0;
    memset (tpp->colorMatrix, 0, sizeof(tpp->colorMatrix));
    tpp->colorMatrix[0][0] = 1.0;
    tpp->colorMatrix[1][1] = 1.0;
    tpp->colorMatrix[2][2] = 1.0;

    if (fixwh==1) {
        w = h * img->width / img->height;
        tpp->scale = (double)img->height / h;
    }
    else {
        h = w * img->height / img->width;
        tpp->scale = (double)img->width / w;
    }

    // bilinear interpolation
    tpp->thumbImg = img->resize (w, h, TI_Bilinear);

    // histogram computation
    tpp->aeHistCompression = 3;
    tpp->aeHistogram = new unsigned int[65536>>tpp->aeHistCompression];
	
	double avg_r = 0;
    double avg_g = 0;
    double avg_b = 0;
    int n = 0;
	
    memset (tpp->aeHistogram, 0, (65536>>tpp->aeHistCompression)*sizeof(int));
    int ix = 0;
    for (int i=0; i<img->height*img->width; i++) {
		int rtmp=CurveFactory::igamma_srgb (img->data[ix++]);
		int gtmp=CurveFactory::igamma_srgb (img->data[ix++]);
		int btmp=CurveFactory::igamma_srgb (img->data[ix++]);
		
		tpp->aeHistogram[rtmp>>tpp->aeHistCompression]++;
		tpp->aeHistogram[gtmp>>tpp->aeHistCompression]+=2;
		tpp->aeHistogram[btmp>>tpp->aeHistCompression]++;

		if (rtmp<64000 && gtmp<64000 && btmp<64000) {
			// autowb computation
			avg_r += rtmp;
            avg_g += gtmp;
            avg_b += btmp;
            n++;
		}
    }

    if (n>0)
        ColorTemp::mul2temp (avg_r/n, avg_g/n, avg_b/n, tpp->autowbTemp, tpp->autowbGreen);

    delete img;
    tpp->init ();
    return tpp;
}

Thumbnail* Thumbnail::loadQuickFromRaw (const Glib::ustring& fname, RawMetaDataLocation& rml, int &w, int &h, int fixwh)
{
	RawImage *ri= new RawImage(fname);
	int r = ri->loadRaw(false,false);
	if( r )
    {
        delete ri;
        return NULL;
    }

    rml.exifBase = ri->get_exifBase();
    rml.ciffBase = ri->get_ciffBase();
    rml.ciffLength = ri->get_ciffLen();

    Image16* img = new Image16 ();

    int err = 1;

    // see if it is something we support
    if ( ri->is_supportedThumb() )
    {
        const char* data((const char*)fdata(ri->get_thumbOffset(),ri->get_file()));
        if ( (unsigned char)data[1] == 0xd8 )
        {
            err = img->loadJPEGFromMemory(data,ri->get_thumbLength());
        }
        else
        {
            err = img->loadPPMFromMemory(data,ri->get_thumbWidth(),ri->get_thumbHeight(),ri->get_thumbSwap(),ri->get_thumbBPS());
        }
    }

    // did we succeed?
    if ( err ) 
    {
        printf("loadfromMemory: error\n");
        delete img;
        delete ri;
        return NULL;
    }

    Thumbnail* tpp = new Thumbnail ();

    tpp->camwbRed = 1.0;
    tpp->camwbGreen = 1.0;
    tpp->camwbBlue = 1.0;

    tpp->embProfileLength = 0;
	tpp->embProfile = NULL;
    tpp->embProfileData = NULL;

    tpp->redMultiplier = 1.0;
    tpp->greenMultiplier = 1.0;
    tpp->blueMultiplier = 1.0;

    tpp->scaleForSave = 8192;
    tpp->defGain = 1.0;
    tpp->gammaCorrected = false;
    tpp->isRaw = 1;
    memset (tpp->colorMatrix, 0, sizeof(tpp->colorMatrix));
    tpp->colorMatrix[0][0] = 1.0;
    tpp->colorMatrix[1][1] = 1.0;
    tpp->colorMatrix[2][2] = 1.0;

    if (fixwh==1) {
        w = h * img->width / img->height;
        tpp->scale = (double)img->height / h;
    }
    else {
        h = w * img->height / img->width;
        tpp->scale = (double)img->width / w;
    }

    tpp->thumbImg = img->resize (w, h, TI_Nearest);
    delete img;

    tpp->autowbTemp=2700;
    tpp->autowbGreen=1.0;

    if (ri->get_rotateDegree() > 0) {
        Image16* rot = tpp->thumbImg->rotate(ri->get_rotateDegree());
        delete tpp->thumbImg;
        tpp->thumbImg = rot;
    }

    tpp->init ();
    delete ri;

    return tpp;
}

#define FISRED(filter,row,col) \
	((filter >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==0 || !filter)
#define FISGREEN(filter,row,col) \
	((filter >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==1 || !filter)
#define FISBLUE(filter,row,col) \
	((filter >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)==2 || !filter)

Thumbnail* Thumbnail::loadFromRaw (const Glib::ustring& fname, RawMetaDataLocation& rml, int &w, int &h, int fixwh)
{
	RawImage *ri= new RawImage (fname);
	int r = ri->loadRaw();
	if( r ){
		delete ri;
		return NULL;
	}
	int width = ri->get_width();
	int height = ri->get_height();
	rtengine::Thumbnail* tpp = new rtengine::Thumbnail;

	tpp->isRaw = true;
	tpp->embProfile = NULL;
	tpp->embProfileData = NULL;
	tpp->embProfileLength = ri->get_profileLen();
	if (ri->get_profileLen())
		tpp->embProfile = cmsOpenProfileFromMem(ri->get_profile(),
				ri->get_profileLen()); //\ TODO check if mutex is needed

	tpp->redMultiplier = ri->get_pre_mul(0);
	tpp->greenMultiplier = ri->get_pre_mul(1);
	tpp->blueMultiplier = ri->get_pre_mul(2);

	float pre_mul[4], scale_mul[4];
	int cblack[4];

	ri->scale_colors();
	ri->pre_interpolate();

	rml.exifBase = ri->get_exifBase();
	rml.ciffBase = ri->get_ciffBase();
	rml.ciffLength = ri->get_ciffLen();

	tpp->camwbRed = tpp->redMultiplier / ri->get_pre_mul(0);
	tpp->camwbGreen = tpp->greenMultiplier / ri->get_pre_mul(1);
	tpp->camwbBlue = tpp->blueMultiplier / ri->get_pre_mul(2);
	tpp->defGain= 1.0/ MIN(MIN(ri->get_pre_mul(0),ri->get_pre_mul(1)),ri->get_pre_mul(2));
	tpp->gammaCorrected = true;

	unsigned filter = ri->get_filters();
	int firstgreen = 1;
	// locate first green location in the first row
	while (!FISGREEN(filter,1,firstgreen))
		firstgreen++;

	int skip = 1;
	if (fixwh == 1) // fix height, scale width
		skip = (ri->get_height() - firstgreen - 1) / h;
	else
		skip = (ri->get_width() - firstgreen - 1) / w;
	if (skip % 2)
		skip--;
	if (skip < 1)
		skip = 1;

	int hskip = skip, vskip = skip;
	if (!ri->get_model().compare("D1X"))
		hskip *= 2;

	int ix = 0;
	int rofs = 0;
	int tmpw = (width - 2) / hskip;
	int tmph = (height - 2) / vskip;

	DCraw::dcrawImage_t image = ri->get_image();

	Image16* tmpImg = new Image16(tmpw, tmph);
	if (ri->isBayer()) {
		for (int row = 1, y = 0; row < height - 1 && y < tmph; row += vskip, y++) {
			rofs = row * width;
			for (int col = firstgreen, x = 0; col < width - 1 && x < tmpw; col+= hskip, x++) {
				int ofs = rofs + col;
				int g = image[ofs][1];
				int r, b;
				if (FISRED(filter,row,col+1)) {
					r = (image[ofs + 1][0] + image[ofs - 1][0]) >> 1;
					b = (image[ofs + width][2] + image[ofs - width][2]) >> 1;
				} else {
					b = (image[ofs + 1][2] + image[ofs - 1][2]) >> 1;
					r = (image[ofs + width][0] + image[ofs - width][0]) >> 1;
				}
				tmpImg->r[y][x] = r;
				tmpImg->g[y][x] = g;
				tmpImg->b[y][x] = b;
			}
		}
	} else {
		for (int row = 1, y = 0; row < height - 1 && y < tmph; row += vskip, y++) {
			rofs = row * width;
			for (int col = firstgreen, x = 0; col < width - 1 && x < tmpw; col
					+= hskip, x++) {
				int ofs = rofs + col;
				tmpImg->r[y][x] = image[ofs][0];
				tmpImg->g[y][x] = image[ofs][1];
				tmpImg->b[y][x] = image[ofs][2];
			}
		}
	}

	if (ri->get_FujiWidth() != 0) {
		int fw = ri->get_FujiWidth() / hskip;
		double step = sqrt(0.5);
		int wide = fw / step;
		int high = (tmph - fw) / step;
		Image16* fImg = new Image16(wide, high);
		float r, c;

		for (int row = 0; row < high; row++)
			for (int col = 0; col < wide; col++) {
				unsigned ur = r = fw + (row - col) * step;
				unsigned uc = c = (row + col) * step;
				if (ur > tmph - 2 || uc > tmpw - 2)
					continue;
				double fr = r - ur;
				double fc = c - uc;
				int oofs = (ur * tmpw + uc) * 3;
				int fofs = (row * wide + col) * 3;
				fImg->r[row][col] = (tmpImg->r[ur][uc] * (1 - fc)
						+ tmpImg->r[ur][uc + 1] * fc) * (1 - fr)
						+ (tmpImg->r[ur + 1][uc] * (1 - fc)
								+ tmpImg->r[ur + 1][uc + 1] * fc) * fr;
				fImg->g[row][col] = (tmpImg->g[ur][uc] * (1 - fc)
						+ tmpImg->g[ur][uc + 1] * fc) * (1 - fr)
						+ (tmpImg->g[ur + 1][uc] * (1 - fc)
								+ tmpImg->g[ur + 1][uc + 1] * fc) * fr;
				fImg->b[row][col] = (tmpImg->b[ur][uc] * (1 - fc)
						+ tmpImg->b[ur][uc + 1] * fc) * (1 - fr)
						+ (tmpImg->b[ur + 1][uc] * (1 - fc)
								+ tmpImg->b[ur + 1][uc + 1] * fc) * fr;
			}
		delete tmpImg;
		tmpImg = fImg;
	}

	if (fixwh == 1) // fix height, scale width
		w = tmpw * h / tmph;
	else
		h = tmph * w / tmpw;

	tpp->thumbImg = tmpImg->resize(w, h, TI_Bilinear);
	delete tmpImg;

	if (ri->get_FujiWidth() != 0)
		tpp->scale = (double) (height - ri->get_FujiWidth()) / sqrt(0.5) / h;
	else
		tpp->scale = (double) height / h;

	// generate histogram for auto exposure
	tpp->aeHistCompression = 3;
	tpp->aeHistogram = new unsigned int[65536 >> tpp->aeHistCompression];
	memset(tpp->aeHistogram, 0, (65536 >> tpp->aeHistCompression) * sizeof(int));
	int radd = 4;
	int gadd = 4;
	int badd = 4;
	if (!filter)
		radd = gadd = badd = 1;
	for (int i = 8; i < height - 8; i++) {
		int start, end;
		if (ri->get_FujiWidth() != 0) {
			int fw = ri->get_FujiWidth();
			start = ABS(fw-i) + 8;
			end = MIN( height + width-fw-i, fw+i) - 8;
		} else {
			start = 8;
			end = width - 8;
		}
		for (int j = start; j < end; j++)
			if (FISGREEN(filter,i,j))
                tpp->aeHistogram[CLIP((int)(tpp->camwbGreen*image[i* width+j][1]))>>tpp->aeHistCompression]+=gadd;
			else if (FISRED(filter,i,j))
                tpp->aeHistogram[CLIP((int)(tpp->camwbRed * image[i* width+j][0]))>>tpp->aeHistCompression]+=radd;
			else if (FISBLUE(filter,i,j))
                tpp->aeHistogram[CLIP((int)(tpp->camwbBlue *image[i* width+j][2]))>>tpp->aeHistCompression]+=badd;
	}

	// generate autoWB
	double avg_r = 0;
	double avg_g = 0;
	double avg_b = 0;
	int rn = 0, gn = 0, bn = 0;

	for (int i = 32; i < height - 32; i++) {
		int start, end;
		if (ri->get_FujiWidth() != 0) {
			int fw = ri->get_FujiWidth();
			start = ABS(fw-i) + 32;
			end = MIN( height + width-fw-i, fw+i) - 32;
		} else {
			start = 32;
			end = width - 32;
		}
		for (int j = start; j < end; j++) {
			if (FISGREEN(filter,i,j)) {
				double d = tpp->defGain * image[i * width + j][1];
				if (d > 64000)
					continue;
				avg_g += d;
				gn++;
			}
			if (FISRED(filter,i,j)) {
				double d = tpp->defGain * image[i * width + j][0];
				if (d > 64000)
					continue;
				avg_r += d;
				rn++;
			}
			if (FISBLUE(filter,i,j)) {
				double d = tpp->defGain * image[i * width + j][2];
				if (d > 64000)
					continue;
				avg_b += d;
				bn++;
			}
		}
	}

	double reds = avg_r / rn * tpp->camwbRed;
	double greens = avg_g / gn * tpp->camwbGreen;
	double blues = avg_b / bn * tpp->camwbBlue;

	double rm = ri->get_rgb_cam(0, 0) * reds + ri->get_rgb_cam(0, 1) * greens + ri->get_rgb_cam(0, 2) * blues;
	double gm = ri->get_rgb_cam(1, 0) * reds + ri->get_rgb_cam(1, 1) * greens + ri->get_rgb_cam(1, 2) * blues;
	double bm = ri->get_rgb_cam(2, 0) * reds + ri->get_rgb_cam(2, 1) * greens + ri->get_rgb_cam(2, 2) * blues;

	ColorTemp::mul2temp(rm, gm, bm, tpp->autowbTemp, tpp->autowbGreen);

	if (ri->get_rotateDegree() > 0) {
		Image16* rot = tpp->thumbImg->rotate(ri->get_rotateDegree());
		delete tpp->thumbImg;
		tpp->thumbImg = rot;
	}

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			tpp->colorMatrix[a][b] = ri->get_rgb_cam(a, b);

	tpp->init();
	delete ri;
	return tpp;
}
#undef FISRED
#undef FISGREEN
#undef FISBLUE


void Thumbnail::init () {

    RawImageSource::inverse33 (colorMatrix, iColorMatrix);
    memset (camToD50, 0, sizeof(camToD50));
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            for (int k=0; k<3; k++)
                camToD50[i][j] += colorMatrix[k][i] * sRGB_d50[k][j];
    camProfile = iccStore->createFromMatrix (camToD50, false, "Camera");
}

bool Thumbnail::igammacomputed = false;
unsigned short Thumbnail::igammatab[256];
unsigned char Thumbnail::gammatab[65536];

Thumbnail::Thumbnail () :
    camProfile(NULL), thumbImg(NULL), aeHistogram(NULL), embProfileData(NULL), embProfile(NULL) {

    if (!igammacomputed) {
        for (int i=0; i<256; i++)
            igammatab[i] = (unsigned short)(255.0*pow(i/255.0,1.0/0.45));
        for (int i=0; i<65536; i++)
            gammatab[i] = (unsigned char)(255.0*pow(i/65535.0,0.45));
        igammacomputed = true;
    }
}

Thumbnail::~Thumbnail () {

    delete thumbImg;
    delete [] aeHistogram;
    delete [] embProfileData;
    if (embProfile)
        cmsCloseProfile(embProfile);
    if (camProfile)
        cmsCloseProfile(camProfile);
}

IImage8* Thumbnail::quickProcessImage (const procparams::ProcParams& params, int rheight, TypeInterpolation interp, double& myscale) {

    int rwidth;
    if (params.coarse.rotate==90 || params.coarse.rotate==270) {
        rwidth = rheight;
        rheight = thumbImg->height * rwidth / thumbImg->width;
    }
    else 
        rwidth = thumbImg->width * rheight / thumbImg->height;   
	Image16* baseImg = thumbImg->resize (rwidth, rheight, interp);

    if (params.coarse.rotate) {
        Image16* tmp = baseImg->rotate (params.coarse.rotate);
        rwidth = tmp->width;
        rheight = tmp->height;
        delete baseImg;
        baseImg = tmp;
    }
    if (params.coarse.hflip) {
        Image16* tmp = baseImg->hflip ();
        delete baseImg;
        baseImg = tmp;
    }
    if (params.coarse.vflip) {
        Image16* tmp = baseImg->vflip ();
        delete baseImg;
        baseImg = tmp;
    }
	Image8* img8 = baseImg->to8();
	delete baseImg;
	return img8;
}

IImage8* Thumbnail::processImage (const procparams::ProcParams& params, int rheight, TypeInterpolation interp, double& myscale) {

    // compute WB multipliers
    ColorTemp currWB = ColorTemp (params.wb.temperature, params.wb.green);
    if (params.wb.method=="Camera") {
        double cam_r = colorMatrix[0][0]*camwbRed + colorMatrix[0][1]*camwbGreen + colorMatrix[0][2]*camwbBlue;
        double cam_g = colorMatrix[1][0]*camwbRed + colorMatrix[1][1]*camwbGreen + colorMatrix[1][2]*camwbBlue;
        double cam_b = colorMatrix[2][0]*camwbRed + colorMatrix[2][1]*camwbGreen + colorMatrix[2][2]*camwbBlue;
        currWB = ColorTemp (cam_r, cam_g, cam_b);
    }
    else if (params.wb.method=="Auto")
        currWB = ColorTemp (autowbTemp, autowbGreen);
    double r, g, b;
    currWB.getMultipliers (r, g, b);
    double rm = iColorMatrix[0][0]*r + iColorMatrix[0][1]*g + iColorMatrix[0][2]*b;
    double gm = iColorMatrix[1][0]*r + iColorMatrix[1][1]*g + iColorMatrix[1][2]*b;
    double bm = iColorMatrix[2][0]*r + iColorMatrix[2][1]*g + iColorMatrix[2][2]*b;
    rm = camwbRed / rm;
    gm = camwbGreen / gm;
    bm = camwbBlue / bm;
    double mul_lum = 0.299*rm + 0.587*gm + 0.114*bm;
    double logDefGain = log(defGain) / log(2.0);
    int rmi, gmi, bmi;
    if (!isRaw || !params.hlrecovery.enabled) {
        logDefGain = 0.0;
        rmi = 1024.0 * rm * defGain / mul_lum;
        gmi = 1024.0 * gm * defGain / mul_lum;
        bmi = 1024.0 * bm * defGain / mul_lum;
    }
    else {
        rmi = 1024.0 * rm / mul_lum;
        gmi = 1024.0 * gm / mul_lum;
        bmi = 1024.0 * bm / mul_lum;  
    }
    // resize to requested width and perform coarse transformation
    int rwidth;
    if (params.coarse.rotate==90 || params.coarse.rotate==270) {
        rwidth = rheight;
        rheight = thumbImg->height * rwidth / thumbImg->width;
    }
    else 
        rwidth = thumbImg->width * rheight / thumbImg->height;   

    Image16* baseImg = thumbImg->resize (rwidth, rheight, interp);
    
    if (params.coarse.rotate) {
        Image16* tmp = baseImg->rotate (params.coarse.rotate);
        rwidth = tmp->width;
        rheight = tmp->height;
        delete baseImg;
        baseImg = tmp;
    }
    if (params.coarse.hflip) {
        Image16* tmp = baseImg->hflip ();
        delete baseImg;
        baseImg = tmp;
    }
    if (params.coarse.vflip) {
        Image16* tmp = baseImg->vflip ();
        delete baseImg;
        baseImg = tmp;
    }
    // apply white balance
    int val;
    for (int i=0; i<rheight; i++)
        for (int j=0; j<rwidth; j++) {
                val = baseImg->r[i][j]*rmi>>10;
                baseImg->r[i][j] = CLIP(val);
                val = baseImg->g[i][j]*gmi>>10;
                baseImg->g[i][j] = CLIP(val);
                val = baseImg->b[i][j]*bmi>>10;
                baseImg->b[i][j] = CLIP(val);
        }

    // apply highlight recovery, if needed
    if (isRaw && params.hlrecovery.enabled) {
        int maxval = 65535 / defGain;
        if (params.hlrecovery.method=="Luminance" || params.hlrecovery.method=="Color") 
            for (int i=0; i<rheight; i++)
                RawImageSource::HLRecovery_Luminance (baseImg->r[i], baseImg->g[i], baseImg->b[i], baseImg->r[i], baseImg->g[i], baseImg->b[i], rwidth, maxval);
        else if (params.hlrecovery.method=="CIELab blending") {
            double icamToD50[3][3];
            RawImageSource::inverse33 (camToD50, icamToD50);
            for (int i=0; i<rheight; i++)
                RawImageSource::HLRecovery_CIELab (baseImg->r[i], baseImg->g[i], baseImg->b[i], baseImg->r[i], baseImg->g[i], baseImg->b[i], rwidth, maxval, camToD50, icamToD50);
        }
    }

    // perform color space transformation
    if (isRaw)
        RawImageSource::colorSpaceConversion (baseImg, params.icm, embProfile, camProfile, camToD50, logDefGain);
    else
        StdImageSource::colorSpaceConversion (baseImg, params.icm, embProfile);
        
    int fw = baseImg->width;
    int fh = baseImg->height;

    ImProcFunctions ipf (&params, false);
    ipf.setScale (sqrt(double(fw*fw+fh*fh))/sqrt(double(thumbImg->width*thumbImg->width+thumbImg->height*thumbImg->height))*scale);

    unsigned int* hist16 = new unsigned int [65536];
    ipf.firstAnalysis (baseImg, &params, hist16, isRaw ? 2.2 : 0.0);

    // perform transform
    if (ipf.needsTransform()) {
        Image16* trImg = new Image16 (fw, fh);
        ipf.transform (baseImg, trImg, 0, 0, 0, 0, fw, fh);
        delete baseImg;
        baseImg = trImg;
    }
    
    // update blurmap
    SHMap* shmap = NULL;
    if (params.sh.enabled) {
        unsigned short** buffer = NULL;
        if (params.sh.hq) {
            buffer = new unsigned short*[fh];
            for (int i=0; i<fh; i++)
                buffer[i] = new unsigned short[fw];
        }
        shmap = new SHMap (fw, fh, false);
        double radius = sqrt (double(fw*fw+fh*fh)) / 2.0;
        double shradius = radius / 1800.0 * params.sh.radius;
        shmap->update (baseImg, buffer, shradius, ipf.lumimul, params.sh.hq);
        if (buffer) {
            for (int i=0; i<fh; i++)
                delete [] buffer[i];
            delete [] buffer;
        }
    }
    
    // RGB processing
    double br = params.toneCurve.expcomp;
    int    bl = params.toneCurve.black;

    if (params.toneCurve.autoexp && aeHistogram) 
        ipf.getAutoExp (aeHistogram, aeHistCompression, logDefGain, params.toneCurve.clip, br, bl);

	float* curve1 = new float [65536];
    float* curve2 = new float [65536];
	int* curve = new int [65536];
    CurveFactory::complexCurve (br, bl/65535.0, params.toneCurve.hlcompr, params.toneCurve.shcompr, params.toneCurve.brightness, params.toneCurve.contrast, logDefGain, isRaw ? 2.2 : 0, true, params.toneCurve.curve, hist16, curve1, curve2, curve, NULL, 16);

    LabImage* labView = new LabImage (baseImg);
    ipf.rgbProc (baseImg, labView, curve1, curve2, curve, shmap, params.toneCurve.saturation);

    if (shmap)
        delete shmap;

    // luminance histogram update
    memset (hist16, 0, 65536*sizeof(int));
    for (int i=0; i<fh; i++)
        for (int j=0; j<fw; j++)
            hist16[labView->L[i][j]]++;

    // luminance processing
    CurveFactory::complexCurve (0.0, 0.0, 0.0, 0.0, params.labCurve.brightness, params.labCurve.contrast, 0.0, 0.0, false, params.labCurve.lcurve, hist16, curve1, curve2, curve, NULL, 16);
    ipf.luminanceCurve (labView, labView, curve, 0, fh);
	CurveFactory::complexsgnCurve (0.0, 100.0, params.labCurve.saturation, 1.0, params.labCurve.acurve, curve, 16);
	ipf.chrominanceCurve (labView, labView, 0, curve, 0, fh);
	CurveFactory::complexsgnCurve (0.0, 100.0, params.labCurve.saturation, 1.0, params.labCurve.bcurve, curve, 16);
    ipf.chrominanceCurve (labView, labView, 1, curve, 0, fh);

	delete [] curve1;
    delete [] curve2;
    delete [] curve;
    delete [] hist16;

    // color processing
    ipf.colorCurve (labView, labView);

    // obtain final image
    Image8* readyImg = new Image8 (fw, fh);
    ipf.lab2rgb (labView, readyImg);
    delete labView;
    delete baseImg;

    // calculate scale
    if (params.coarse.rotate==90 || params.coarse.rotate==270) 
        myscale = scale * thumbImg->width / fh;
    else
        myscale = scale * thumbImg->height / fh;

    if (params.resize.enabled) {
        if (params.resize.dataspec==0)
            myscale *= params.resize.scale;
        else if (params.resize.dataspec==1)
            myscale *= (double)params.resize.width / (params.coarse.rotate==90 || params.coarse.rotate==270 ? thumbImg->height : thumbImg->width) / scale;
        else if (params.resize.dataspec==2)
            myscale *= (double)params.resize.height / (params.coarse.rotate==90 || params.coarse.rotate==270 ? thumbImg->width : thumbImg->height) / scale;
    }
    myscale = 1.0 / myscale;

/*    // apply crop
    if (params.crop.enabled) {
        int ix = 0;
        for (int i=0; i<fh; i++) 
            for (int j=0; j<fw; j++)
                if (i<params.crop.y/myscale || i>(params.crop.y+params.crop.h)/myscale || j<params.crop.x/myscale || j>(params.crop.x+params.crop.w)/myscale) {
                    readyImg->data[ix++] /= 3;
                    readyImg->data[ix++] /= 3;
                    readyImg->data[ix++] /= 3;
                }
                else
                    ix += 3;
    }*/
    return readyImg;
}

int Thumbnail::getImageWidth (const procparams::ProcParams& params, int rheight) {

    int rwidth;
    if (params.coarse.rotate==90 || params.coarse.rotate==270) 
        rwidth = thumbImg->height * rheight / thumbImg->width;
    else 
        rwidth = thumbImg->width * rheight / thumbImg->height;   

    return rwidth;
}

void Thumbnail::getFinalSize (const rtengine::procparams::ProcParams& params, int& fullw, int& fullh) {

    double fw = thumbImg->width*scale;
    double fh = thumbImg->height*scale;
    
    if (params.coarse.rotate==90 || params.coarse.rotate==270) {
        fh = thumbImg->width*scale;
        fw = thumbImg->height*scale;
    }
    if (!params.resize.enabled) {
        fullw = fw;
        fullh = fh;
    }
    else if (params.resize.dataspec==0) {
        fullw = fw*params.resize.scale;
        fullh = fh*params.resize.scale;
    }
    else if (params.resize.dataspec==1) {
        fullw = params.resize.width;
        fullh = (double)fh*params.resize.width/(params.coarse.rotate==90 || params.coarse.rotate==270 ? fh : fw);
    }
    else if (params.resize.dataspec==2) {
        fullw = (double)fw*params.resize.height/(params.coarse.rotate==90 || params.coarse.rotate==270 ? fw : fh);
        fullh = params.resize.height;
    }
}

void Thumbnail::getCamWB (double& temp, double& green) {

    double cam_r = colorMatrix[0][0]*camwbRed + colorMatrix[0][1]*camwbGreen + colorMatrix[0][2]*camwbBlue;
    double cam_g = colorMatrix[1][0]*camwbRed + colorMatrix[1][1]*camwbGreen + colorMatrix[1][2]*camwbBlue;
    double cam_b = colorMatrix[2][0]*camwbRed + colorMatrix[2][1]*camwbGreen + colorMatrix[2][2]*camwbBlue;
    ColorTemp currWB = ColorTemp (cam_r, cam_g, cam_b);
    temp = currWB.getTemp ();
    green = currWB.getGreen ();
}

void Thumbnail::getAutoWB (double& temp, double& green) {
    
    temp = autowbTemp;
    green = autowbGreen;
}

void Thumbnail::applyAutoExp (procparams::ProcParams& params) {

    if (params.toneCurve.autoexp && aeHistogram) {
        ImProcFunctions ipf (&params, false);
        ipf.getAutoExp (aeHistogram, aeHistCompression, log(defGain) / log(2.0), params.toneCurve.clip, params.toneCurve.expcomp, params.toneCurve.black);
    }
}

void Thumbnail::getSpotWB (const procparams::ProcParams& params, int xp, int yp, int rect, double& rtemp, double& rgreen) {

    std::vector<Coord2D> points, red, green, blue;
    for (int i=yp-rect; i<=yp+rect; i++)
        for (int j=xp-rect; j<=xp+rect; j++) 
            points.push_back (Coord2D (j, i));

    int fw = thumbImg->width, fh = thumbImg->height;
    if (params.coarse.rotate==90 || params.coarse.rotate==270) {
        fw = thumbImg->height;
        fh = thumbImg->width;
    }
    ImProcFunctions ipf (&params, false);
    ipf.transCoord (fw, fh, points, red, green, blue);
    int tr = TR_NONE;
    if (params.coarse.rotate==90)  tr |= TR_R90;
    if (params.coarse.rotate==180) tr |= TR_R180;
    if (params.coarse.rotate==270) tr |= TR_R270;
    if (params.coarse.hflip)       tr |= TR_HFLIP;
    if (params.coarse.vflip)       tr |= TR_VFLIP;

    // calculate spot wb (copy & pasted from stdimagesource)
    unsigned short igammatab[256];
    for (int i=0; i<256; i++)
        igammatab[i] = (unsigned short)(255.0*pow(i/255.0,1.0/0.45));
    int x; int y;
    double reds = 0, greens = 0, blues = 0;
    int rn = 0, gn = 0, bn = 0;
    for (int i=0; i<red.size(); i++) {
        transformPixel (red[i].x, red[i].y, tr, x, y);
        if (x>=0 && y>=0 && x<thumbImg->width && y<thumbImg->height) {
            reds += thumbImg->r[y][x];
            rn++;
        }
        transformPixel (green[i].x, green[i].y, tr, x, y);
        if (x>=0 && y>=0 && x<thumbImg->width && y<thumbImg->height) {
            greens += thumbImg->g[y][x];
            gn++;
        }
        transformPixel (blue[i].x, blue[i].y, tr, x, y);
        if (x>=0 && y>=0 && x<thumbImg->width && y<thumbImg->height) {
            blues += thumbImg->b[y][x];
            bn++;
        }
    }
    reds = reds/rn * camwbRed;
    greens = greens/gn * camwbGreen;
    blues = blues/bn * camwbBlue;

    double rm = colorMatrix[0][0]*reds + colorMatrix[0][1]*greens + colorMatrix[0][2]*blues;
    double gm = colorMatrix[1][0]*reds + colorMatrix[1][1]*greens + colorMatrix[1][2]*blues;
    double bm = colorMatrix[2][0]*reds + colorMatrix[2][1]*greens + colorMatrix[2][2]*blues;

    ColorTemp ct (rm, gm, bm);
    rtemp = ct.getTemp ();
    rgreen = ct.getGreen ();
}
void Thumbnail::transformPixel (int x, int y, int tran, int& tx, int& ty) {
    
    int W = thumbImg->width;
    int H = thumbImg->height;
    int sw = W, sh = H;  
    if ((tran & TR_ROT) == TR_R90 || (tran & TR_ROT) == TR_R270) {
        sw = H;
        sh = W;
    }

    int ppx = x, ppy = y;
    if (tran & TR_HFLIP) 
        ppx = sw - 1 - x ;
    if (tran & TR_VFLIP) 
        ppy = sh - 1 - y;
    
    tx = ppx;
    ty = ppy;
    
    if ((tran & TR_ROT) == TR_R180) {
        tx = W - 1 - ppx;
        ty = H - 1 - ppy;
    }
    else if ((tran & TR_ROT) == TR_R90) {
        tx = ppy;
        ty = H - 1 - ppx;
    }
    else if ((tran & TR_ROT) == TR_R270) {
        tx = W - 1 - ppy;
        ty = ppx;
    }
    tx/=scale;
    ty/=scale;
}

// format: 1=8bit direct, 2=16bit direct, 3=JPG
bool Thumbnail::writeImage (const Glib::ustring& fname, int format) {

    if (!thumbImg)
        return false;
    
    if (format==1 || format==3) {
        // to utilize the 8 bit color range of the thumbnail we brighten it and apply gamma correction
        unsigned char* tmpdata = new unsigned char[thumbImg->height*thumbImg->width*3];
        int ix = 0,max;

        if (gammaCorrected) {
            // if it's gamma correct (usually a RAW), we have the problem that there is a lot noise etc. that makes the maximum way too high.
            // Strategy is limit a certain percent of pixels so the overall picture quality when scaling to 8 bit is way better
            const double BurnOffPct=0.03;  // *100 = percent pixels that may be clipped

            // Calc the histogram
            unsigned int* hist16 = new unsigned int [65536];
            memset(hist16,0,sizeof(int)*65536);

            for (int row=0; row<thumbImg->height; row++)
                for (int col=0; col<thumbImg->width; col++) {
                    hist16[thumbImg->r[row][col]]++;
                    hist16[thumbImg->g[row][col]]+=2;  // Bayer 2x green correction
                    hist16[thumbImg->b[row][col]]++;
                }

            // Go down till we cut off that many pixels
            unsigned long cutoff = thumbImg->height * thumbImg->height * 4 * BurnOffPct;

            int max; unsigned long sum=0;
            for (max=65535; max>16384 && sum<cutoff; max--) sum+=hist16[max];

            delete[] hist16;

            scaleForSave = 65535*8192 / max;

            // Correction and gamma to 8 Bit
            for (int i=0; i<thumbImg->height; i++)
                for (int j=0; j<thumbImg->width; j++) {
                    tmpdata[ix++] = gammatab[MIN(thumbImg->r[i][j],max) * scaleForSave >> 13];
                    tmpdata[ix++] = gammatab[MIN(thumbImg->g[i][j],max) * scaleForSave >> 13];
                    tmpdata[ix++] = gammatab[MIN(thumbImg->b[i][j],max) * scaleForSave >> 13];
                }
        }
        else {
            // If it's not gamma corrected (usually a JPG) we take the normal maximum
            max=0;

            for (int row=0; row<thumbImg->height; row++)
                for (int col=0; col<thumbImg->width; col++) {
                    if (thumbImg->r[row][col]>max) max = thumbImg->r[row][col];
                    if (thumbImg->g[row][col]>max) max = thumbImg->g[row][col];
                    if (thumbImg->b[row][col]>max) max = thumbImg->b[row][col];
                }
            
            if (max < 16384) max = 16384;
            scaleForSave = 65535*8192 / max;

            // Correction and gamma to 8 Bit
            for (int i=0; i<thumbImg->height; i++)
                for (int j=0; j<thumbImg->width; j++) {
                    tmpdata[ix++] = thumbImg->r[i][j]*scaleForSave >> 21;
                    tmpdata[ix++] = thumbImg->g[i][j]*scaleForSave >> 21;
                    tmpdata[ix++] = thumbImg->b[i][j]*scaleForSave >> 21;
                }
        }

        if (format==1) {
            FILE* f = g_fopen (fname.c_str(), "wb");
            if (!f) {
                delete [] tmpdata;
                return false;
            }
            fwrite (&thumbImg->width, 1, sizeof (int), f);
            fwrite (&thumbImg->height, 1, sizeof (int), f);
            fwrite (tmpdata, thumbImg->width*thumbImg->height, 3, f);
            fclose (f);
        }
        else if (format==3) {
            FILE* f = g_fopen (fname.c_str(), "wb");
            if (!f) {
                delete [] tmpdata;
                return false;
            }
        	jpeg_compress_struct cinfo;
        	jpeg_error_mgr jerr;
	        cinfo.err = jpeg_std_error (&jerr);
	        jpeg_create_compress (&cinfo);
        	jpeg_stdio_dest (&cinfo, f);
            cinfo.image_width  = thumbImg->width;
	        cinfo.image_height = thumbImg->height;
	        cinfo.in_color_space = JCS_RGB;
	        cinfo.input_components = 3;
	        jpeg_set_defaults (&cinfo);
            cinfo.write_JFIF_header = FALSE;
			
            // compute optimal Huffman coding tables for the image. Bit slower to generate, but size of result image is a bit less (default was FALSE)
            cinfo.optimize_coding = TRUE;

            // Since math coprocessors are common these days, FLOAT should be a bit more accurate AND fast (default is ISLOW)
            // (machine dependency is not really an issue, since we all run on x86 and having exactly the same file is not a requirement)
            cinfo.dct_method = JDCT_FLOAT;

            jpeg_set_quality (&cinfo, 87, true);
        	jpeg_start_compress(&cinfo, TRUE);
            int rowlen = thumbImg->width*3;
        	while (cinfo.next_scanline < cinfo.image_height) {
                unsigned char* row = tmpdata + cinfo.next_scanline*thumbImg->width*3;
        		if (jpeg_write_scanlines (&cinfo, &row, 1) < 1) {
                    jpeg_finish_compress (&cinfo);
	                jpeg_destroy_compress (&cinfo);
        	        fclose (f);
                    delete [] tmpdata;
                    return false;
                }
            }
        	jpeg_finish_compress (&cinfo);
        	jpeg_destroy_compress (&cinfo);
        	fclose (f);
        }
        delete [] tmpdata;
        return true;
    }
    else if (format==2) {
        FILE* f = g_fopen (fname.c_str(), "wb");
        if (!f)
            return false;
        fwrite (&thumbImg->width, 1, sizeof (int), f);
        fwrite (&thumbImg->height, 1, sizeof (int), f);
        for (int i=0; i<thumbImg->height; i++)
            fwrite (thumbImg->r[i], thumbImg->width, 2, f);
        for (int i=0; i<thumbImg->height; i++)
            fwrite (thumbImg->g[i], thumbImg->width, 2, f);
        for (int i=0; i<thumbImg->height; i++)
            fwrite (thumbImg->b[i], thumbImg->width, 2, f);
        fclose (f);
        return true;
    }
    else
        return false;
}

bool Thumbnail::readImage (const Glib::ustring& fname) {
    
    delete thumbImg;
    thumbImg = NULL;

    int imgType = 0;
    if (Glib::file_test (fname+".cust16", Glib::FILE_TEST_EXISTS))
        imgType = 2;
    if (Glib::file_test (fname+".cust", Glib::FILE_TEST_EXISTS))
        imgType = 1;
    else if (Glib::file_test (fname+".jpg", Glib::FILE_TEST_EXISTS))
        imgType = 3;

    if (!imgType) 
        return false;
    else if (imgType==1) {
        FILE* f = g_fopen ((fname+".cust").c_str(), "rb");
        if (!f)
            return false;
        int width, height;
        fread (&width, 1, sizeof (int), f);
        fread (&height, 1, sizeof (int), f);
        unsigned char* tmpdata = new unsigned char [width*height*3];
        fread (tmpdata, width*height, 3, f);
        fclose (f);
        thumbImg = new Image16 (width, height);
        int ix = 0, val;
        for (int i=0; i<height; i++)
            for (int j=0; j<width; j++)
                if (gammaCorrected) {
                    val = igammatab[tmpdata[ix++]]*256*8192/scaleForSave;
                    thumbImg->r[i][j] = CLIP(val);
                    val = igammatab[tmpdata[ix++]]*256*8192/scaleForSave;
                    thumbImg->g[i][j] = CLIP(val);
                    val = igammatab[tmpdata[ix++]]*256*8192/scaleForSave;
                    thumbImg->b[i][j] = CLIP(val);
                }
                else {
                    val = tmpdata[ix++]*256*8192/scaleForSave;
                    thumbImg->r[i][j] = CLIP(val);
                    val = tmpdata[ix++]*256*8192/scaleForSave;
                    thumbImg->g[i][j] = CLIP(val);
                    val = tmpdata[ix++]*256*8192/scaleForSave;
                    thumbImg->b[i][j] = CLIP(val);
                }
        delete [] tmpdata;
        return true;
    }
    else if (imgType==2) {
        FILE* f = g_fopen ((fname+".cust16").c_str(), "rb");
        if (!f)
            return false;
        int width, height;
        fread (&width, 1, sizeof (int), f);
        fread (&height, 1, sizeof (int), f);
        thumbImg = new Image16 (width, height);
        for (int i=0; i<height; i++)
            fread (thumbImg->r[i], width, 2, f);
        for (int i=0; i<height; i++)
            fread (thumbImg->g[i], width, 2, f);
        for (int i=0; i<height; i++)
            fread (thumbImg->b[i], width, 2, f);
        fclose (f);
        return true;
    }
    else if (imgType==3) {
        FILE* f = g_fopen ((fname+".jpg").c_str(), "rb");
        if (!f) 
            return false;
        struct jpeg_decompress_struct cinfo;
        struct jpeg_error_mgr jerr;
        if (!setjmp(jpeg_jmp_buf)) {
            cinfo.err = my_jpeg_std_error (&jerr);
            jpeg_create_decompress (&cinfo);
            my_jpeg_stdio_src (&cinfo,f);
            jpeg_read_header (&cinfo, TRUE);
            int width, height;
            width = cinfo.image_width;
            height = cinfo.image_height;
            cinfo.dct_method = JDCT_FASTEST;
            cinfo.do_fancy_upsampling = 1;
            jpeg_start_decompress(&cinfo);
            thumbImg = new Image16 (width, height);
            unsigned char* row = new unsigned char [width*3];
            while (cinfo.output_scanline < cinfo.output_height) {
                jpeg_read_scanlines (&cinfo, &row, 1);
                int ix = 0, val;
                for (int j=0; j<width; j++) {
                    if (gammaCorrected) {
                        val = igammatab[row[ix++]]*256*8192/scaleForSave;
                        thumbImg->r[cinfo.output_scanline-1][j] = CLIP(val);
                        val = igammatab[row[ix++]]*256*8192/scaleForSave;
                        thumbImg->g[cinfo.output_scanline-1][j] = CLIP(val);
                        val = igammatab[row[ix++]]*256*8192/scaleForSave;
                        thumbImg->b[cinfo.output_scanline-1][j] = CLIP(val);
                    }
                    else {
                        val = row[ix++]*256*8192/scaleForSave;
                        thumbImg->r[cinfo.output_scanline-1][j] = CLIP(val);
                        val = row[ix++]*256*8192/scaleForSave;
                        thumbImg->g[cinfo.output_scanline-1][j] = CLIP(val);
                        val = row[ix++]*256*8192/scaleForSave;
                        thumbImg->b[cinfo.output_scanline-1][j] = CLIP(val);
                    }
                }
            }
            jpeg_finish_decompress (&cinfo);
            jpeg_destroy_decompress (&cinfo);
            fclose (f);
            delete [] row;
            return true;
        }
        else {
            fclose (f);
            return false;
        }
        return true;
    }
    return false;
}

bool Thumbnail::readData  (const Glib::ustring& fname) {

    SafeKeyFile keyFile;
    
    try {
        if (!keyFile.load_from_file (fname)) 
            return false;

        if (keyFile.has_group ("LiveThumbData")) { 
            if (keyFile.has_key ("LiveThumbData", "CamWBRed"))          camwbRed            = keyFile.get_double ("LiveThumbData", "CamWBRed");
            if (keyFile.has_key ("LiveThumbData", "CamWBGreen"))        camwbGreen          = keyFile.get_double ("LiveThumbData", "CamWBGreen");
            if (keyFile.has_key ("LiveThumbData", "CamWBBlue"))         camwbBlue           = keyFile.get_double ("LiveThumbData", "CamWBBlue");
            if (keyFile.has_key ("LiveThumbData", "AutoWBTemp"))        autowbTemp          = keyFile.get_double ("LiveThumbData", "AutoWBTemp");
            if (keyFile.has_key ("LiveThumbData", "AutoWBGreen"))       autowbGreen         = keyFile.get_double ("LiveThumbData", "AutoWBGreen");
            if (keyFile.has_key ("LiveThumbData", "AEHistCompression")) aeHistCompression   = keyFile.get_integer ("LiveThumbData", "AEHistCompression");
            if (keyFile.has_key ("LiveThumbData", "RedMultiplier"))     redMultiplier       = keyFile.get_double ("LiveThumbData", "RedMultiplier");
            if (keyFile.has_key ("LiveThumbData", "GreenMultiplier"))   greenMultiplier     = keyFile.get_double ("LiveThumbData", "GreenMultiplier");
            if (keyFile.has_key ("LiveThumbData", "BlueMultiplier"))    blueMultiplier      = keyFile.get_double ("LiveThumbData", "BlueMultiplier");
            if (keyFile.has_key ("LiveThumbData", "Scale"))             scale               = keyFile.get_double ("LiveThumbData", "Scale");
            if (keyFile.has_key ("LiveThumbData", "DefaultGain"))       defGain             = keyFile.get_double ("LiveThumbData", "DefaultGain");
            if (keyFile.has_key ("LiveThumbData", "ScaleForSave"))      scaleForSave        = keyFile.get_integer ("LiveThumbData", "ScaleForSave");
            if (keyFile.has_key ("LiveThumbData", "GammaCorrected"))    gammaCorrected      = keyFile.get_boolean ("LiveThumbData", "GammaCorrected");
            if (keyFile.has_key ("LiveThumbData", "ColorMatrix")) {
                std::vector<double> cm = keyFile.get_double_list ("LiveThumbData", "ColorMatrix");
                int ix = 0;
                for (int i=0; i<3; i++)
                    for (int j=0; j<3; j++)
                        colorMatrix[i][j] = cm[ix++];
            }
        }
        return true;
    }
    catch (Glib::Error) {
        return false;
    }
}

bool Thumbnail::writeData  (const Glib::ustring& fname) {

    SafeKeyFile keyFile;

    try {
    if( Glib::file_test(fname,Glib::FILE_TEST_EXISTS) )
        keyFile.load_from_file (fname); 
    } catch (...) {}

    keyFile.set_double  ("LiveThumbData", "CamWBRed", camwbRed);
    keyFile.set_double  ("LiveThumbData", "CamWBGreen", camwbGreen);
    keyFile.set_double  ("LiveThumbData", "CamWBBlue", camwbBlue);
    keyFile.set_double  ("LiveThumbData", "AutoWBTemp", autowbTemp);
    keyFile.set_double  ("LiveThumbData", "AutoWBGreen", autowbGreen);
    keyFile.set_integer ("LiveThumbData", "AEHistCompression", aeHistCompression);
    keyFile.set_double  ("LiveThumbData", "RedMultiplier", redMultiplier);
    keyFile.set_double  ("LiveThumbData", "GreenMultiplier", greenMultiplier);
    keyFile.set_double  ("LiveThumbData", "BlueMultiplier", blueMultiplier);
    keyFile.set_double  ("LiveThumbData", "Scale", scale);
    keyFile.set_double  ("LiveThumbData", "DefaultGain", defGain);
    keyFile.set_integer ("LiveThumbData", "ScaleForSave", scaleForSave);
    keyFile.set_boolean ("LiveThumbData", "GammaCorrected", gammaCorrected);
    Glib::ArrayHandle<double> cm ((double*)colorMatrix, 9, Glib::OWNERSHIP_NONE);
    keyFile.set_double_list ("LiveThumbData", "ColorMatrix", cm);

    FILE *f = g_fopen (fname.c_str(), "wt");
    if (!f)
        return false;
    else {
        fprintf (f, "%s", keyFile.to_data().c_str());
        fclose (f);
        return true;
    }
}

bool Thumbnail::readEmbProfile  (const Glib::ustring& fname) {

    FILE* f = fopen (fname.c_str(), "rb");
    if (!f) {
        embProfileData = NULL;
        embProfile = NULL;
        embProfileLength = 0;
    }
    else {
        fseek (f, 0, SEEK_END);
        embProfileLength = ftell (f);
        fseek (f, 0, SEEK_SET);
        embProfileData = new unsigned char[embProfileLength];
        fread (embProfileData, 1, embProfileLength, f);
        fclose (f);
        embProfile = cmsOpenProfileFromMem (embProfileData, embProfileLength);
        return true;
    }
    return false;
}

bool Thumbnail::writeEmbProfile (const Glib::ustring& fname) {
    
    if (embProfileData) {
        FILE* f = fopen (fname.c_str(), "wb");
        if (f) {
            fwrite (embProfileData, 1, embProfileLength, f);
            fclose (f);
            return true;
        }
    }
    return false;
}

bool Thumbnail::readAEHistogram  (const Glib::ustring& fname) {

    FILE* f = fopen (fname.c_str(), "rb");
    if (!f) 
        aeHistogram = NULL;
    else {
        aeHistogram = new unsigned int[65536>>aeHistCompression];
        fread (aeHistogram, 1, (65536>>aeHistCompression)*sizeof(int), f);
        fclose (f);
        return true;
    }
    return false;
}

bool Thumbnail::writeAEHistogram (const Glib::ustring& fname) {

    if (aeHistogram) {
        FILE* f = fopen (fname.c_str(), "wb");
        if (f) {
            fwrite (aeHistogram, 1, (65536>>aeHistCompression)*sizeof(int), f);
            fclose (f);
            return true;
        }
    }
    return false;
}

}
