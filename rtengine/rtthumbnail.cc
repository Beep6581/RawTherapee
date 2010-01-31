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
    tpp->aeHistogram = new int[65536>>tpp->aeHistCompression];
    memset (tpp->aeHistogram, 0, (65536>>tpp->aeHistCompression)*sizeof(int));
    int ix = 0;
    for (int i=0; i<img->height*img->width; i++) {
            tpp->aeHistogram[CurveFactory::igamma_srgb (img->data[ix++])>>tpp->aeHistCompression]++;
            tpp->aeHistogram[CurveFactory::igamma_srgb (img->data[ix++])>>tpp->aeHistCompression]++;
            tpp->aeHistogram[CurveFactory::igamma_srgb (img->data[ix++])>>tpp->aeHistCompression]++;
    }

    // autowb computation
    double avg_r = 0;
    double avg_g = 0;
    double avg_b = 0;
    int n = 0;
    int p = 6;
    for (int i=1; i<img->height-1; i++)
        for (int j=1; j<img->width-1; j++) {
            int ofs = 3*(i*img->width + j);
            if (img->data[ofs]>250 || img->data[ofs+1]>250 || img->data[ofs+2]>250)
                continue;
            avg_r += StdImageSource::intpow((double)img->data[ofs]*256, p);
            avg_g += StdImageSource::intpow((double)img->data[ofs+1]*256, p);
            avg_b += StdImageSource::intpow((double)img->data[ofs+2]*256, p);
            n++;
        }
    ColorTemp::mul2temp (pow(avg_r/n, 1.0/p), pow(avg_g/n, 1.0/p), pow(avg_b/n, 1.0/p), tpp->autowbTemp, tpp->autowbGreen);

    delete img;
    tpp->init ();
    return tpp;
}

void Thumbnail::init () {

    RawImageSource::inverse33 (colorMatrix, iColorMatrix);
    memset (camToD50, 0, sizeof(camToD50));
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            for (int k=0; k<3; k++)
                camToD50[i][j] += colorMatrix[k][i] * sRGB_d50[k][j];
    camProfile = iccStore.createFromMatrix (camToD50, false, "Camera");
}

bool Thumbnail::igammacomputed = false;
unsigned short Thumbnail::igammatab[256];
unsigned char Thumbnail::gammatab[65536];

Thumbnail::Thumbnail () :
    embProfile(NULL), camProfile(NULL), aeHistogram(NULL), thumbImg(NULL), embProfileData(NULL) {

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
    // resize to requested with and perform coarse transformation
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

    // appy highlight recovery, if needed
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

    ImProcFunctions ipf;
    int* hist16 = new int [65536];
    ipf.firstAnalysis (baseImg, &params, hist16, isRaw ? 2.2 : 0.0);

    // perform transform
    bool needstransform  = fabs(params.rotate.degree)>1e-15 || fabs(params.distortion.amount)>1e-15 || fabs(params.cacorrection.red)>1e-15 || fabs(params.cacorrection.blue)>1e-15;
    bool needsvignetting = params.vignetting.amount!=0;

    if (!needstransform && needsvignetting) {
        Image16* trImg = new Image16 (fw, fh);
        ipf.vignetting (baseImg, trImg, &params, 0, 0, fw, fh);
        delete baseImg;
        baseImg = trImg;
    }
    else if (needstransform) {
        Image16* trImg = new Image16 (fw, fh);
        ipf.simpltransform (baseImg, trImg, &params, 0, 0, 0, 0, fw, fh);
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
        shmap = new SHMap (fw, fh);
        double radius = sqrt (fw*fw+fh*fh) / 2.0;
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

    int* curve = new int [65536];
    CurveFactory::updateCurve3 (curve, hist16, params.toneCurve.curve, logDefGain, br, bl, params.toneCurve.hlcompr, params.toneCurve.shcompr, params.toneCurve.brightness, params.toneCurve.contrast,  isRaw ? 2.2 : 0, true);

    LabImage* labView = new LabImage (baseImg);
    ipf.rgbProc (baseImg, labView, &params, curve, shmap);

    if (shmap)
        delete shmap;

    // luminance histogram update
    memset (hist16, 0, 65536*sizeof(int));
    for (int i=0; i<fh; i++)
        for (int j=0; j<fw; j++)
            hist16[labView->L[i][j]]++;

    // luminance processing
    CurveFactory::updateCurve2 (curve, hist16, params.lumaCurve.curve, 0, params.lumaCurve.brightness, params.lumaCurve.black, params.lumaCurve.hlcompr, params.lumaCurve.shcompr, params.lumaCurve.contrast, 0.0, false);
    ipf.luminanceCurve (labView, labView, curve, 0, fh);

    delete [] curve;
    delete [] hist16;

    // color processing
    ipf.colorCurve (labView, labView, &params);

    // obtain final image
    Image8* readyImg = new Image8 (fw, fh);
    ipf.lab2rgb (labView, readyImg);
    ipf.release ();
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

    if (params.toneCurve.autoexp && aeHistogram) 
        ImProcFunctions::getAutoExp (aeHistogram, aeHistCompression, log(defGain) / log(2.0), params.toneCurve.clip, params.toneCurve.expcomp, params.toneCurve.black);
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
    ImProcFunctions ipf;  
    ipf.transCoord (&params, fw, fh, points, red, green, blue);
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

bool Thumbnail::writeImage (const Glib::ustring& fname, int format) {

    if (!thumbImg)
        return false;
    
    if (format==1 || format==3) {
        // to utilize the 8 bit color range of the thumbnail we brighten it and apply gamma correction
        int max = 0;
        for (int row=0; row<thumbImg->height; row++)
            for (int col=0; col<thumbImg->width; col++) {
                if (thumbImg->r[row][col]>max)
                    max = thumbImg->r[row][col];
                if (thumbImg->g[row][col]>max)
                    max = thumbImg->r[row][col];
                if (thumbImg->b[row][col]>max)
                    max = thumbImg->r[row][col];
            }
        if (max < 16384)
            max = 16384;
        scaleForSave = 65535*8192 / max;
        unsigned char* tmpdata = new unsigned char[thumbImg->height*thumbImg->width*3];
        int ix = 0;
        if (gammaCorrected) {
            for (int i=0; i<thumbImg->height; i++)
                for (int j=0; j<thumbImg->width; j++) {
                    tmpdata[ix++] = gammatab[thumbImg->r[i][j]*scaleForSave >> 13];
                    tmpdata[ix++] = gammatab[thumbImg->g[i][j]*scaleForSave >> 13];
                    tmpdata[ix++] = gammatab[thumbImg->b[i][j]*scaleForSave >> 13];
                }
        }
        else {
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
            jpeg_set_quality (&cinfo, 85, true);
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
}

bool Thumbnail::readData  (const Glib::ustring& fname) {

    Glib::KeyFile keyFile;
    
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

    Glib::KeyFile keyFile;

    try {
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
    
    if (embProfileLength) {
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
        aeHistogram = new int[65536>>aeHistCompression];
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
