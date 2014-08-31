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
#include <cmath>
#include <iostream>

#include "rtengine.h"
#include "rawimagesource.h"
#include "rawimagesource_i.h"
#include "median.h"
#include "rawimage.h"
#include "mytime.h"
#include "iccmatrices.h"
#include "iccstore.h"
#include "image8.h"
#include "curves.h"
#include "dfmanager.h"
#include "ffmanager.h"
#include "slicer.h"
#include "../rtgui/options.h"
#include "dcp.h"
#include "rt_math.h"
#include "improcfun.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "opthelper.h"

namespace rtengine {

extern const Settings* settings;
#undef ABS
#undef DIST

#define ABS(a) ((a)<0?-(a):(a))
#define DIST(a,b) (ABS(a-b))
	
#define PIX_SORT(a,b) { if ((a)>(b)) {temp=(a);(a)=(b);(b)=temp;} }
	
#define med3x3(a0,a1,a2,a3,a4,a5,a6,a7,a8,median) { \
p[0]=a0; p[1]=a1; p[2]=a2; p[3]=a3; p[4]=a4; p[5]=a5; p[6]=a6; p[7]=a7; p[8]=a8; \
PIX_SORT(p[1],p[2]); PIX_SORT(p[4],p[5]); PIX_SORT(p[7],p[8]); \
PIX_SORT(p[0],p[1]); PIX_SORT(p[3],p[4]); PIX_SORT(p[6],p[7]); \
PIX_SORT(p[1],p[2]); PIX_SORT(p[4],p[5]); PIX_SORT(p[7],p[8]); \
PIX_SORT(p[0],p[3]); PIX_SORT(p[5],p[8]); PIX_SORT(p[4],p[7]); \
PIX_SORT(p[3],p[6]); PIX_SORT(p[1],p[4]); PIX_SORT(p[2],p[5]); \
PIX_SORT(p[4],p[7]); PIX_SORT(p[4],p[2]); PIX_SORT(p[6],p[4]); \
PIX_SORT(p[4],p[2]); median=p[4];} //a4 is the median
	
#define med5(a0,a1,a2,a3,a4,median) { \
p[0]=a0; p[1]=a1; p[2]=a2; p[3]=a3; p[4]=a4; \
PIX_SORT(p[0],p[1]) ; PIX_SORT(p[3],p[4]) ; PIX_SORT(p[0],p[3]) ; \
PIX_SORT(p[1],p[4]) ; PIX_SORT(p[1],p[2]) ; PIX_SORT(p[2],p[3]) ; \
PIX_SORT(p[1],p[2]) ; median=p[2] ;}
	

RawImageSource::RawImageSource ()
:ImageSource()
,plistener(NULL)
,border(4)
,ri(NULL)
,cache(NULL)
,rawData(0,0)
,green(0,0)
,red(0,0)
,blue(0,0)
{
    hrmap[0] = NULL;
    hrmap[1] = NULL;
    hrmap[2] = NULL;
	//needhr = NULL;
    //hpmap = NULL;
	camProfile = NULL;
	embProfile = NULL;
	rgbSourceModified = false;
	hlmax[0] = hlmax[1] = hlmax[2] = hlmax[3] = 0.f;
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RawImageSource::~RawImageSource () {

    delete idata;
    if (ri) {
        delete ri;
    }

    flushRGB();
    flushRawData();

    if( cache )
        delete [] cache;
    if (hrmap[0]!=NULL) {
        int dh = H/HR_SCALE;
        freeArray<float>(hrmap[0], dh);
        freeArray<float>(hrmap[1], dh);
        freeArray<float>(hrmap[2], dh);
    }
    //if (needhr)
    //    freeArray<char>(needhr, H);
    //if (hpmap)
    //    freeArray<char>(hpmap, H);
    if (camProfile)
        cmsCloseProfile (camProfile);
    if (embProfile)
        cmsCloseProfile (embProfile);
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::transformRect (PreviewProps pp, int tran, int &ssx1, int &ssy1, int &width, int &height, int &fw) {

    pp.x += border;
    pp.y += border;

    if (d1x) {
        if ((tran & TR_ROT) == TR_R90 || (tran & TR_ROT) == TR_R270) {
            pp.x /= 2;
            pp.w = pp.w/2+1;
        }
        else {
            pp.y /= 2;
            pp.h = pp.h/2+1;
        }
    }

    int w = W, h = H;
    if (fuji) {
        w = ri->get_FujiWidth() * 2 + 1;
        h = (H - ri->get_FujiWidth())*2 + 1;
    }
    
    int sw = w, sh = h;  
    if ((tran & TR_ROT) == TR_R90 || (tran & TR_ROT) == TR_R270) {
        sw = h;
        sh = w;
    }
    if( pp.w > sw-2*border) pp.w = sw-2*border;
    if( pp.h > sh-2*border) pp.h = sh-2*border;

    int ppx = pp.x, ppy = pp.y;
    if (tran & TR_HFLIP) 
        ppx = sw - pp.x - pp.w;
    if (tran & TR_VFLIP) 
        ppy = sh - pp.y - pp.h;

    int sx1 = ppx;        // assuming it's >=0
    int sy1 = ppy;        // assuming it's >=0
    int sx2 = max(ppx + pp.w, w-1);
    int sy2 = max(ppy + pp.h, h-1);

    if ((tran & TR_ROT) == TR_R180) {
        sx1 = max(w - ppx - pp.w, 0);
        sy1 = max(h - ppy - pp.h, 0);
        sx2 = min(sx1 + pp.w, w-1);
        sy2 = min(sy1 + pp.h, h-1);
    }
    else if ((tran & TR_ROT) == TR_R90) {
        sx1 = ppy;
        sy1 = max(h - ppx - pp.w, 0);
        sx2 = min(sx1 + pp.h, w-1);
        sy2 = min(sy1 + pp.w, h-1);
    }
    else if ((tran & TR_ROT) == TR_R270) {
        sx1 = max(w - ppy - pp.h, 0);
        sy1 = ppx;
        sx2 = min(sx1 + pp.h, w-1);
        sy2 = min(sy1 + pp.w, h-1);
    }

    if (fuji) {
        // atszamoljuk a koordinatakat fuji-ra:
        // recalculate the coordinates fuji-ra:
        ssx1 = (sx1+sy1) / 2;
        ssy1 = (sy1 - sx2 ) / 2 + ri->get_FujiWidth();
        int ssx2 = (sx2+sy2) / 2 + 1;
        int ssy2 = (sy2 - sx1) / 2 + ri->get_FujiWidth();
        fw   = (sx2 - sx1) / 2 / pp.skip;
        width  = (ssx2 - ssx1) / pp.skip + ((ssx2 - ssx1) % pp.skip > 0);
        height = (ssy2 - ssy1) / pp.skip + ((ssy2 - ssy1) % pp.skip > 0); 
    }
    else {
        ssx1 = sx1;
        ssy1 = sy1;
        width  = (sx2 - sx1) / pp.skip + ((sx2 - sx1) % pp.skip > 0);
        height = (sy2 - sy1) / pp.skip + ((sy2 - sy1) % pp.skip > 0); 
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
static float
calculate_scale_mul(float scale_mul[4], const float pre_mul_[4], const float c_white[4], const float c_black[4], bool isMono, int colors)
{
    if (isMono || colors == 1) {
        for (int c = 0; c < 4; c++) {
            scale_mul[c] = 65535.0 / (c_white[c] - c_black[c]);
        }
    } else {
        float pre_mul[4];
        for (int c = 0; c < 4; c++) {
            pre_mul[c] = pre_mul_[c];
        }
        if (pre_mul[3] == 0) {
            pre_mul[3] = pre_mul[1]; // G2 == G1
        }
        float maxpremul = max(pre_mul[0], pre_mul[1], pre_mul[2], pre_mul[3]);
        for (int c = 0; c < 4; c++) {
            scale_mul[c] = (pre_mul[c] / maxpremul) * 65535.0 / (c_white[c] - c_black[c]);
        }
    }
    float gain = max(scale_mul[0], scale_mul[1], scale_mul[2], scale_mul[3]) / min(scale_mul[0], scale_mul[1], scale_mul[2], scale_mul[3]);
    return gain;
}

void RawImageSource::getImage (ColorTemp ctemp, int tran, Imagefloat* image, PreviewProps pp, ToneCurveParams  hrp, ColorManagementParams cmp, RAWParams raw )
{
    MyMutex::MyLock lock(getImageMutex);

    tran = defTransform (tran);

    // compute channel multipliers
    double r, g, b;
    float rm, gm, bm;
    ctemp.getMultipliers (r, g, b);
    rm = imatrices.cam_rgb[0][0]*r + imatrices.cam_rgb[0][1]*g + imatrices.cam_rgb[0][2]*b;
    gm = imatrices.cam_rgb[1][0]*r + imatrices.cam_rgb[1][1]*g + imatrices.cam_rgb[1][2]*b;
    bm = imatrices.cam_rgb[2][0]*r + imatrices.cam_rgb[2][1]*g + imatrices.cam_rgb[2][2]*b;

    if (true) {
        // adjust gain so the maximum raw value of the least scaled channel just hits max
        const float new_pre_mul[4] = { ri->get_pre_mul(0) / rm, ri->get_pre_mul(1) / gm, ri->get_pre_mul(2) / bm, ri->get_pre_mul(3) / gm };
        float new_scale_mul[4];

        bool isMono = (ri->getSensorType()==ST_FUJI_XTRANS && raw.xtranssensor.method == RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::mono])
                   || (ri->getSensorType()==ST_BAYER && raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::mono]);
        float gain = calculate_scale_mul(new_scale_mul, new_pre_mul, c_white, cblacksom, isMono, ri->get_colors());
        rm = new_scale_mul[0] / scale_mul[0] * gain;
        gm = new_scale_mul[1] / scale_mul[1] * gain;
        bm = new_scale_mul[2] / scale_mul[2] * gain;
        //fprintf(stderr, "camera gain: %f, current wb gain: %f, diff in stops %f\n", camInitialGain, gain, log2(camInitialGain) - log2(gain));
    } else {
        // old scaling: used a fixed reference gain based on camera (as-shot) white balance

        // how much we need to scale each channel to get our new white balance
        rm = refwb_red / rm;
        gm = refwb_green / gm;
        bm = refwb_blue / bm;
        // normalize so larger multiplier becomes 1.0
        float minval = min(rm, gm, bm);
        rm /= minval;
        gm /= minval;
        bm /= minval;
        // multiply with reference gain, ie as-shot WB
        rm *= camInitialGain;
        gm *= camInitialGain;
        bm *= camInitialGain;
    }

    defGain=0.0;
    // compute image area to render in order to provide the requested part of the image
    int sx1, sy1, imwidth, imheight, fw;
    transformRect (pp, tran, sx1, sy1, imwidth, imheight, fw);

    // check possible overflows
    int maximwidth, maximheight;
    if ((tran & TR_ROT) == TR_R90 || (tran & TR_ROT) == TR_R270) {
        maximwidth = image->height;
        maximheight = image->width;
    }
    else {
        maximwidth = image->width;
        maximheight = image->height;
    }
    if (d1x)
        maximheight /= 2;
        
    // correct if overflow (very rare), but not fuji because it is corrected in transline
    if (!fuji && imwidth>maximwidth) 
        imwidth = maximwidth;
    if (!fuji && imheight>maximheight)
        imheight = maximheight;
	
    int maxx=this->W,maxy=this->H,skip=pp.skip;

    //if (sx1+skip*imwidth>maxx) imwidth --; // very hard to fix this situation without an 'if' in the loop.
    float area=skip*skip;
    rm/=area;
    gm/=area;
    bm/=area;
	
	hlmax[0]=chmax[0]*rm*area;
	hlmax[1]=chmax[1]*gm*area;
	hlmax[2]=chmax[2]*bm*area;


#ifdef _OPENMP
#pragma omp parallel if(!d1x)		// omp disabled for D1x to avoid race conditions (see Issue 1088 http://code.google.com/p/rawtherapee/issues/detail?id=1088)
    {
#endif
    // render the requested image part
    float* line_red  = new float[imwidth];
    float* line_grn  = new float[imwidth];
    float* line_blue = new float[imwidth];
	//printf("clip[0]=%f  clip[1]=%f  clip[2]=%f\n",hlmax[0],hlmax[1],hlmax[2]);


#ifdef _OPENMP
#pragma omp for
#endif
	for (int ix=0; ix<imheight; ix++) { int i=sy1+skip*ix;if (i>=maxy-skip) i=maxy-skip-1; // avoid trouble
		if (ri->getSensorType()!=ST_NONE || ri->get_colors() == 1) {
            for (int j=0,jx=sx1; j<imwidth; j++,jx+=skip) {if (jx>=maxx-skip) jx=maxx-skip-1; // avoid trouble
            	float rtot,gtot,btot;
            	rtot=gtot=btot=0;
				for (int m=0; m<skip; m++)
					for (int n=0; n<skip; n++)
					{
						rtot += red[i+m][jx+n];
						gtot += green[i+m][jx+n];
						btot += blue[i+m][jx+n];
					}
				rtot*=rm;
				gtot*=gm;
				btot*=bm;
				if (!hrp.hrenabled)
				{
					rtot=CLIP(rtot);
					gtot=CLIP(gtot);
					btot=CLIP(btot);
				}
				line_red[j] = rtot;
				line_grn[j] = gtot;
				line_blue[j] = btot;
            }
        } else {
            for (int j=0,jx=sx1; j<imwidth; j++,jx+=skip) {if (jx>maxx-skip) jx=maxx-skip-1;
            	float rtot,gtot,btot;
            	rtot=gtot=btot=0;
				for (int m=0; m<skip; m++)
					for (int n=0; n<skip; n++)
					{
						rtot += rawData[i+m][(jx+n)*3+0];
						gtot += rawData[i+m][(jx+n)*3+1];
						btot += rawData[i+m][(jx+n)*3+2];
					}				
				rtot*=rm;
				gtot*=gm;
				btot*=bm;
				if (!hrp.hrenabled)
				{
					rtot=CLIP(rtot);
					gtot=CLIP(gtot);
					btot=CLIP(btot);
				}
				line_red[j] = rtot;
				line_grn[j] = gtot;
				line_blue[j] = btot;
				
            }
        }

		//process all highlight recovery other than "Color"
        if (hrp.hrenabled && hrp.method!="Color")
			hlRecovery (hrp.method, line_red, line_grn, line_blue, i, sx1, imwidth, skip, raw, hlmax);

        transLine (line_red, line_grn, line_blue, ix, image, tran, imwidth, imheight, fw);
		
    }

    delete [] line_red;
    delete [] line_grn;
    delete [] line_blue;
#ifdef _OPENMP
    }
#endif
    if (fuji) {   
        int a = ((tran & TR_ROT) == TR_R90 && image->width%2==0) || ((tran & TR_ROT) == TR_R180 && image->height%2+image->width%2==1) || ((tran & TR_ROT) == TR_R270 && image->height%2==0);
        // first row
        for (int j=1+a; j<image->width-1; j+=2) {
          image->r(0,j) = (image->r(1,j) + image->r(0,j+1) + image->r(0,j-1)) / 3;
          image->g(0,j) = (image->g(1,j) + image->g(0,j+1) + image->g(0,j-1)) / 3;
          image->b(0,j) = (image->b(1,j) + image->b(0,j+1) + image->b(0,j-1)) / 3;
        }
        // other rows
        for (int i=1; i<image->height-1; i++) {
          for (int j=2-(a+i+1)%2; j<image->width-1; j+=2) {
              // edge-adaptive interpolation
              double dh = (ABS(image->r(i,j+1) - image->r(i,j-1)) + ABS(image->g(i,j+1) - image->g(i,j-1)) + ABS(image->b(i,j+1) - image->b(i,j-1))) / 1.0;
              double dv = (ABS(image->r(i+1,j) - image->r(i-1,j)) + ABS(image->g(i+1,j) - image->g(i-1,j)) + ABS(image->b(i+1,j) - image->b(i-1,j))) / 1.0;
              double eh = 1.0 / (1.0 + dh);
              double ev = 1.0 / (1.0 + dv);
              image->r(i,j) = (eh * (image->r(i,j+1) + image->r(i,j-1)) + ev * (image->r(i+1,j) + image->r(i-1,j))) / (2.0 * (eh + ev));
              image->g(i,j) = (eh * (image->g(i,j+1) + image->g(i,j-1)) + ev * (image->g(i+1,j) + image->g(i-1,j))) / (2.0 * (eh + ev));
              image->b(i,j) = (eh * (image->b(i,j+1) + image->b(i,j-1)) + ev * (image->b(i+1,j) + image->b(i-1,j))) / (2.0 * (eh + ev));
          }
          // first pixel
          if (2-(a+i+1)%2==2) {
              image->r(i,0) = (image->r(i+1,0) + image->r(i-1,0) + image->r(i,1)) / 3;
              image->g(i,0) = (image->g(i+1,0) + image->g(i-1,0) + image->g(i,1)) / 3;
              image->b(i,0) = (image->b(i+1,0) + image->b(i-1,0) + image->b(i,1)) / 3;
          }
          // last pixel
          if (2-(a+i+image->width)%2==2) {
              image->r(i,image->width-1) = (image->r(i+1,image->width-1) + image->r(i-1,image->width-1) + image->r(i,image->width-2)) / 3;
              image->g(i,image->width-1) = (image->g(i+1,image->width-1) + image->g(i-1,image->width-1) + image->g(i,image->width-2)) / 3;
              image->b(i,image->width-1) = (image->b(i+1,image->width-1) + image->b(i-1,image->width-1) + image->b(i,image->width-2)) / 3;
          }
        }
        // last row
        int b = (a==1 && image->height%2) || (a==0 && image->height%2==0);
        for (int j=1+b; j<image->width-1; j+=2) {
          image->r(image->height-1,j) = (image->r(image->height-2,j) + image->r(image->height-1,j+1) + image->r(image->height-1,j-1)) / 3;
          image->g(image->height-1,j) = (image->g(image->height-2,j) + image->g(image->height-1,j+1) + image->g(image->height-1,j-1)) / 3;
          image->b(image->height-1,j) = (image->b(image->height-2,j) + image->b(image->height-1,j+1) + image->b(image->height-1,j-1)) / 3;
        }
    }


	
    // Flip if needed
    if (tran & TR_HFLIP)
        hflip (image);
    if (tran & TR_VFLIP)
        vflip (image);

    // Color correction (only when running on full resolution)
    if (ri->getSensorType()!=ST_NONE && pp.skip==1) {
        if (ri->getSensorType()==ST_BAYER)
            processFalseColorCorrection (image, raw.bayersensor.ccSteps);
        else if (ri->getSensorType()==ST_FUJI_XTRANS)
            processFalseColorCorrection (image, raw.xtranssensor.ccSteps);
    }
    // *** colorSpaceConversion was here ***
    //colorSpaceConversion (image, cmp, raw, embProfile, camProfile, xyz_cam, (static_cast<const ImageData*>(getMetaData()))->getCamera());
}

void RawImageSource::convertColorSpace(Imagefloat* image, ColorManagementParams cmp, ColorTemp &wb, RAWParams raw) {
    double pre_mul[3] = { ri->get_pre_mul(0), ri->get_pre_mul(1), ri->get_pre_mul(2) };
    colorSpaceConversion (image, cmp, wb, pre_mul, raw, embProfile, camProfile, imatrices.xyz_cam, (static_cast<const ImageData*>(getMetaData()))->getCamera());
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/* cfaCleanFromMap: correct raw pixels looking at the bitmap
 * takes into consideration if there are multiple bad pixels in the neighborhood
 */
int RawImageSource::cfaCleanFromMap( PixelsMap &bitmapBads )
{
	float eps=1.0;	
	int counter=0;
	for( int row = 2; row < H-2; row++ ){
		for(int col = 2; col <W-2; col++ ){
			int sk = bitmapBads.skipIfZero(col,row); //optimization for a stripe all zero
			if( sk ){
				col +=sk-1; //-1 is because of col++ in cycle
				continue;
			}
			if( ! bitmapBads.get(col,row ) )
				continue;

			double wtdsum=0,norm=0,sum=0,tot=0;
			for( int dy=-2;dy<=2;dy+=2){
				for( int dx=-2;dx<=2;dx+=2){
					if (dy==0 && dx==0) continue;
					if( bitmapBads.get(col+dx,row+dy) ) continue;
					sum += rawData[row+dy][col+dx];
					tot++;
					if (bitmapBads.get(col-dx,row-dy)) continue;

					double dirwt = 1/( fabs( rawData[row+dy][col+dx]- rawData[row-dy][col-dx])+eps);
					wtdsum += dirwt* rawData[row+dy][col+dx];
					norm += dirwt;
				}
			}
			if (norm > 0.0){
				rawData[row][col]= wtdsum / norm;//gradient weighted average
				counter++;
			} else {
				if (tot > 0.1) rawData[row][col] = sum/tot;//backup plan -- simple average
			}
		}
	}
	return counter;
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/*  Search for hot or dead pixels in the image and update the map
 *  For each pixel compare its value to the average of similar color surrounding
 *  (Taken from Emil Martinec idea)
 *  (Optimized by Ingo Weyrich 2013)
 */
int RawImageSource::findHotDeadPixel( PixelsMap &bpMap, float thresh)
{
	// counter for dead or hot pixels
	int counter=0;
	
	// allocate temporary buffer
	float (*cfablur);
	cfablur = (float (*)) malloc (H*W * sizeof *cfablur);
	
#pragma omp parallel
{
#pragma omp for
	for (int i=0; i<H; i++) {
		int iprev,inext,jprev,jnext;
		float p[9],temp;
		if (i<2) {iprev=i+2;} else {iprev=i-2;}
		if (i>H-3) {inext=i-2;} else {inext=i+2;}
		for (int j=0; j<W; j++) {
			if (j<2) {jprev=j+2;} else {jprev=j-2;}
			if (j>W-3) {jnext=j-2;} else {jnext=j+2;}
			med3x3(rawData[iprev][jprev],rawData[iprev][j],rawData[iprev][jnext],
				   rawData[i][jprev],rawData[i][j],rawData[i][jnext],
				   rawData[inext][jprev],rawData[inext][j],rawData[inext][jnext],temp);
			cfablur[i*W+j] = fabs(rawData[i][j]-temp);
		}
	}
#pragma omp for reduction(+:counter) schedule (dynamic,16)
	//cfa pixel heat/death evaluation
	for (int rr=0; rr < H; rr++) {
		int top=max(0,rr-2);
		int bottom=min(H-1,rr+2);
		int rrmWpcc = rr*W;
		for (int cc=0; cc < W; cc++,rrmWpcc++) {
			//evaluate pixel for heat/death
			float pixdev = cfablur[rrmWpcc];
			float hfnbrave = -pixdev;
			int left=max(0,cc-2);
			int right=min(W-1,cc+2);
			for (int mm=top; mm<=bottom; mm++) {
				int mmmWpnn = mm*W+left;
				for (int nn=left; nn<=right; nn++,mmmWpnn++) {
					hfnbrave += cfablur[mmmWpnn];
				}
			}
			if (pixdev * ((bottom-top+1)*(right-left+1)-1) > thresh*hfnbrave) {
				// mark the pixel as "bad"
				bpMap.set(cc,rr);
				counter++;
			}
		}//end of pixel evaluation
	}
}//end of parallel processing
	free (cfablur);
	return counter;
}

void RawImageSource::rotateLine (float* line, PlanarPtr<float> &channel, int tran, int i, int w, int h) {

    if ((tran & TR_ROT) == TR_R180) 
        for (int j=0; j<w; j++)
            channel(h-1-i,w-1-j) = line[j];

    else if ((tran & TR_ROT) == TR_R90)
        for (int j=0; j<w; j++) 
            channel(j,h-1-i) = line[j];

    else if ((tran & TR_ROT) == TR_R270)
        for (int j=0; j<w; j++) 
            channel(w-1-j,i) = line[j];
    else 
        for (int j=0; j<w; j++)
            channel(i,j) = line[j];
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::transLine (float* red, float* green, float* blue, int i, Imagefloat* image, int tran, int imwidth, int imheight, int fw) {

  // Fuji SuperCCD rotation + coarse rotation
  if (fuji) {
      int start = ABS(fw-i);
      int w = fw * 2 + 1;
      int h = (imheight - fw)*2 + 1;

      if ((tran & TR_ROT) == TR_R180) {
          int end = min(h+fw-i, w-fw+i);
          for (int j=start; j<end; j++) {
            int y = i+j-fw;
            int x = fw-i+j;
            if (x>=0 && y<image->height && y>=0 && x<image->width) {
                image->r(image->height-1-y,image->width-1-x) = red[j];
                image->g(image->height-1-y,image->width-1-x) = green[j];
                image->b(image->height-1-y,image->width-1-x) = blue[j];
            }
          }
      }
      else if ((tran & TR_ROT) == TR_R270) {
          int end = min(h+fw-i, w-fw+i);
          for (int j=start; j<end; j++) {
            int y = i+j-fw;
            int x = fw-i+j;
            if (x>=0 && x<image->height && y>=0 && y<image->width) {
                image->r(image->height-1-x,y) = red[j];
                image->g(image->height-1-x,y) = green[j];
                image->b(image->height-1-x,y) = blue[j];
            }
          }
      }
      else if ((tran & TR_ROT) == TR_R90) {
          int end = min(h+fw-i, w-fw+i);
          for (int j=start; j<end; j++) {
            int y = i+j-fw;
            int x = fw-i+j;
            if (x>=0 && y<image->width && y>=0 && x<image->height) {
                image->r(x,image->width-1-y) = red[j];
                image->g(x,image->width-1-y) = green[j];
                image->b(x,image->width-1-y) = blue[j];
            }
          }
      }
      else {
        int end = min(h+fw-i, w-fw+i);
        for (int j=start; j<end; j++) {
            int y = i+j-fw;
            int x = fw-i+j;
            if (x>=0 && y<image->height && y>=0 && x<image->width) {
                image->r(y,x) = red[j];
                image->g(y,x) = green[j];
                image->b(y,x) = blue[j];
            }
        }
      }
  }
  // Nikon D1X vertical interpolation + coarse rotation
  else if (d1x) {
    // copy new pixels
    if ((tran & TR_ROT) == TR_R180) {
      for (int j=0; j<imwidth; j++) {
        image->r(2*imheight-2-2*i,imwidth-1-j) = red[j];
        image->g(2*imheight-2-2*i,imwidth-1-j) = green[j];
        image->b(2*imheight-2-2*i,imwidth-1-j) = blue[j];
      }

      if (i==1 || i==2) { // linear interpolation
        int row = 2*imheight-1-2*i;
        for (int j=0; j<imwidth; j++) {
          int col = imwidth-1-j;
          image->r(row,col) = (red[j] + image->r(row+1,col)) /2;
          image->g(row,col) = (green[j] + image->g(row+1,col)) /2;
          image->b(row,col) = (blue[j] + image->b(row+1,col)) /2;
        }
      }
      else if (i==imheight-1) {
        int row = 2*imheight-1-2*i;
        for (int j=0; j<imwidth; j++) {
          int col = imwidth-1-j;
          image->r(row,col) = (red[j] + image->r(row+1,col)) /2;
          image->g(row,col) = (green[j] + image->g(row+1,col)) /2;
          image->b(row,col) = (blue[j] + image->b(row+1,col)) /2;
        }
        row = 2*imheight-1-2*i+2;
        for (int j=0; j<imwidth; j++) {
          int col = imwidth-1-j;
          image->r(row,col) = (red[j] + image->r(row+1,col)) /2;
          image->g(row,col) = (green[j] + image->g(row+1,col)) /2;
          image->b(row,col) = (blue[j] + image->b(row+1,col)) /2;
        }
      }
      else if (i>2 && i<imheight-1) { // vertical bicubic interpolation
        int row = 2*imheight-1-2*i+2;
        for (int j=0; j<imwidth; j++) {
          int col = imwidth-1-j;
          image->r(row,col) = CLIP((int)(-0.0625*red[j] + 0.5625*image->r(row-1,col) + 0.5625*image->r(row+1,col) - 0.0625*image->r(row+3,col)));
          image->g(row,col) = CLIP((int)(-0.0625*green[j] + 0.5625*image->g(row-1,col) + 0.5625*image->g(row+1,col) - 0.0625*image->g(row+3,col)));
          image->b(row,col) = CLIP((int)(-0.0625*blue[j] + 0.5625*image->b(row-1,col) + 0.5625*image->b(row+1,col) - 0.0625*image->b(row+3,col)));
        }
      }
    }
    else if ((tran & TR_ROT) == TR_R90) {
      for (int j=0; j<imwidth; j++) {
        image->r(j,2*imheight-2-2*i) = red[j];
        image->g(j,2*imheight-2-2*i) = green[j];
        image->b(j,2*imheight-2-2*i) = blue[j];
      }
      if (i==1 || i==2) { // linear interpolation
        int col = 2*imheight-1-2*i;
        for (int j=0; j<imwidth; j++) {
          image->r(j,col) = (red[j] + image->r(j,col+1)) /2;
          image->g(j,col) = (green[j] + image->g(j,col+1)) /2;
          image->b(j,col) = (blue[j] + image->b(j,col+1)) /2;
        }
      }
      else if (i==imheight-1) {
        int col = 2*imheight-1-2*i;
        for (int j=0; j<imwidth; j++) {
          image->r(j,col) = (red[j] + image->r(j,col+1)) /2;
          image->g(j,col) = (green[j] + image->g(j,col+1)) /2;
          image->b(j,col) = (blue[j] + image->b(j,col+1)) /2;
        }
        col = 2*imheight-1-2*i+2;
        for (int j=0; j<imwidth; j++) {
          image->r(j,col) = (red[j] + image->r(j,col+1)) /2;
          image->g(j,col) = (green[j] + image->g(j,col+1)) /2;
          image->b(j,col) = (blue[j] + image->b(j,col+1)) /2;
        }
      }
      else if (i>2 && i<imheight-1) { // vertical bicubic interpolation
        int col = 2*imheight-1-2*i+2;
        for (int j=0; j<imwidth; j++) {
          image->r(j,col) = CLIP((int)(-0.0625*red[j] + 0.5625*image->r(j,col-1) + 0.5625*image->r(j,col+1) - 0.0625*image->r(j,col+3)));
          image->g(j,col) = CLIP((int)(-0.0625*green[j] + 0.5625*image->g(j,col-1) + 0.5625*image->g(j,col+1) - 0.0625*image->g(j,col+3)));
          image->b(j,col) = CLIP((int)(-0.0625*blue[j] + 0.5625*image->b(j,col-1) + 0.5625*image->b(j,col+1) - 0.0625*image->b(j,col+3)));
        }
      }
    }
    else if ((tran & TR_ROT) == TR_R270) {
      for (int j=0; j<imwidth; j++) {
        image->r(imwidth-1-j,2*i) = red[j];
        image->g(imwidth-1-j,2*i) = green[j];
        image->b(imwidth-1-j,2*i) = blue[j];
      }
      if (i==1 || i==2) { // linear interpolation
        for (int j=0; j<imwidth; j++) {
          int row = imwidth-1-j;
          image->r(row,2*i-1) = (red[j] + image->r(row,2*i-2)) * 0.5f;
          image->g(row,2*i-1) = (green[j] + image->g(row,2*i-2)) * 0.5f;
          image->b(row,2*i-1) = (blue[j] + image->b(row,2*i-2)) * 0.5f;
        }
      }
      else if (i==imheight-1) {
        for (int j=0; j<imwidth; j++) {
          int row = imwidth-1-j;
          image->r(row,2*i-1) = (red[j] + image->r(row,2*i-2)) * 0.5f;
          image->g(row,2*i-1) = (green[j] + image->g(row,2*i-2)) * 0.5f;
          image->b(row,2*i-1) = (blue[j] + image->b(row,2*i-2)) * 0.5f;
          image->r(row,2*i-3) = (image->r(row,2*i-2) + image->r(row,2*i-4)) * 0.5f;
          image->g(row,2*i-3) = (image->g(row,2*i-2) + image->g(row,2*i-4)) * 0.5f;
          image->b(row,2*i-3) = (image->b(row,2*i-2) + image->b(row,2*i-4)) * 0.5f;
        }
      }
      else if (i>0 && i<imheight-1) { // vertical bicubic interpolationi
        for (int j=0; j<imwidth; j++) {
          int row = imwidth-1-j;
          image->r(row,2*i-3) = CLIP((int)(-0.0625*red[j] + 0.5625*image->r(row,2*i-2) + 0.5625*image->r(row,2*i-4) - 0.0625*image->r(row,2*i-6)));
          image->g(row,2*i-3) = CLIP((int)(-0.0625*green[j] + 0.5625*image->g(row,2*i-2) + 0.5625*image->g(row,2*i-4) - 0.0625*image->g(row,2*i-6)));
          image->b(row,2*i-3) = CLIP((int)(-0.0625*blue[j] + 0.5625*image->b(row,2*i-2) + 0.5625*image->b(row,2*i-4) - 0.0625*image->b(row,2*i-6)));
        }
      }    
    }
    else {
        rotateLine (red, image->r, tran, 2*i, imwidth, imheight);
        rotateLine (green, image->g, tran, 2*i, imwidth, imheight);
        rotateLine (blue, image->b, tran, 2*i, imwidth, imheight);

      if (i==1 || i==2) { // linear interpolation
        for (int j=0; j<imwidth; j++) {
          image->r(2*i-1,j) = (red[j] + image->r(2*i-2,j)) /2;
          image->g(2*i-1,j) = (green[j] + image->g(2*i-2,j)) /2;
          image->b(2*i-1,j) = (blue[j] + image->b(2*i-2,j)) /2;
        }
      }
      else if (i==imheight-1) {
            for (int j=0; j<imwidth; j++) {
              image->r(2*i-3,j) = (image->r(2*i-4,j) + image->r(2*i-2,j)) /2;
              image->g(2*i-3,j) = (image->g(2*i-4,j) + image->g(2*i-2,j)) /2;
              image->b(2*i-3,j) = (image->b(2*i-4,j) + image->b(2*i-2,j)) /2;
              image->r(2*i-1,j) = (red[j] + image->r(2*i-2,j)) /2;
              image->g(2*i-1,j) = (green[j] + image->g(2*i-2,j)) /2;
              image->b(2*i-1,j) = (blue[j] + image->b(2*i-2,j)) /2;
            }
      }
      else if (i>2 && i<imheight-1) { // vertical bicubic interpolationi
        for (int j=0; j<imwidth; j++) {
          image->r(2*i-3,j) = CLIP((int)(-0.0625*red[j] + 0.5625*image->r(2*i-2,j) + 0.5625*image->r(2*i-4,j) - 0.0625*image->r(2*i-6,j)));
          image->g(2*i-3,j) = CLIP((int)(-0.0625*green[j] + 0.5625*image->g(2*i-2,j) + 0.5625*image->g(2*i-4,j) - 0.0625*image->g(2*i-6,j)));
          image->b(2*i-3,j) = CLIP((int)(-0.0625*blue[j] + 0.5625*image->b(2*i-2,j) + 0.5625*image->b(2*i-4,j) - 0.0625*image->b(2*i-6,j)));
        }
      }
    }
  }  // if nikon dx1
  // other (conventional) CCD coarse rotation
  else {
    rotateLine (red, image->r, tran, i, imwidth, imheight);
    rotateLine (green, image->g, tran, i, imwidth, imheight);
    rotateLine (blue, image->b, tran, i, imwidth, imheight);
  }
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::getFullSize (int& w, int& h, int tr) {

    tr = defTransform (tr);

    if (fuji) {
        w = ri->get_FujiWidth() * 2 + 1;
        h = (H - ri->get_FujiWidth())*2 + 1;
    }
    else if (d1x) {
        w = W;
        h = 2*H-1;
    }
    else {
        w = W;
        h = H;
    }
   
    if ((tr & TR_ROT) == TR_R90 || (tr & TR_ROT) == TR_R270) {
        int tmp = w;
        w = h;
        h = tmp;
    }
    w -= 2 * border;
    h -= 2 * border;
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
void RawImageSource::getSize (int tran, PreviewProps pp, int& w, int& h) {

    tran = defTransform (tran);

//    if (fuji) {
//        return;
//    }
//    else if (d1x) {
//        return;
//    }
//    else {
        w = pp.w / pp.skip + (pp.w % pp.skip > 0);
        h = pp.h / pp.skip + (pp.h % pp.skip > 0);
//    }
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::hflip (Imagefloat* image) {
    image->hflip();
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::vflip (Imagefloat* image) {
    image->vflip();
}
	
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
int RawImageSource::load (Glib::ustring fname, bool batch) {

	MyTime t1,t2;
	t1.set();
    fileName = fname;

    if (plistener) {
        plistener->setProgressStr ("Decoding...");
        plistener->setProgress (0.0);
    }

    ri = new RawImage(fname);
    int errCode = ri->loadRaw (true, true, plistener, 0.8);
    if (errCode) return errCode;

    ri->compress_image();
    if (plistener) {
        plistener->setProgress (0.9);
    }
/***** Copy once constant data extracted from raw *******/
    W = ri->get_width();
    H = ri->get_height();
    fuji = ri->get_FujiWidth()!=0;
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
        	imatrices.rgb_cam[i][j] = ri->get_rgb_cam(i,j);
    // compute inverse of the color transformation matrix
	// first arg is matrix, second arg is inverse
    inverse33 (imatrices.rgb_cam, imatrices.cam_rgb);

    d1x  = ! ri->get_model().compare("D1X");
    if (d1x)
        border = 8;

    if(ri->getSensorType()==ST_FUJI_XTRANS)
		border = 7;

    if ( ri->get_profile() )
        embProfile = cmsOpenProfileFromMem (ri->get_profile(), ri->get_profileLen());

    // create profile
    memset (imatrices.xyz_cam, 0, sizeof(imatrices.xyz_cam));
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            for (int k=0; k<3; k++)
            	imatrices.xyz_cam[i][j] += xyz_sRGB[i][k] * imatrices.rgb_cam[k][j];
    camProfile = iccStore->createFromMatrix (imatrices.xyz_cam, false, "Camera");
    inverse33 (imatrices.xyz_cam, imatrices.cam_xyz);

        for (int c = 0; c < 4; c++) {
            c_white[c] = ri->get_white(c);
        }
        // First we get the "as shot" ("Camera") white balance and store it
	float pre_mul[4];
        // FIXME: get_colorsCoeff not so much used nowadays, when we have calculate_scale_mul() function here
        ri->get_colorsCoeff( pre_mul, scale_mul, c_black, false);//modify  for black level
        camInitialGain = max(scale_mul[0], scale_mul[1], scale_mul[2], scale_mul[3]) / min(scale_mul[0], scale_mul[1], scale_mul[2], scale_mul[3]);

	double camwb_red = ri->get_pre_mul(0) / pre_mul[0];
	double camwb_green = ri->get_pre_mul(1) / pre_mul[1];
	double camwb_blue = ri->get_pre_mul(2) / pre_mul[2];
	double cam_r = imatrices.rgb_cam[0][0]*camwb_red + imatrices.rgb_cam[0][1]*camwb_green + imatrices.rgb_cam[0][2]*camwb_blue;
	double cam_g = imatrices.rgb_cam[1][0]*camwb_red + imatrices.rgb_cam[1][1]*camwb_green + imatrices.rgb_cam[1][2]*camwb_blue;
	double cam_b = imatrices.rgb_cam[2][0]*camwb_red + imatrices.rgb_cam[2][1]*camwb_green + imatrices.rgb_cam[2][2]*camwb_blue;
	camera_wb = ColorTemp (cam_r, cam_g, cam_b, 1.); // as shot WB

	ColorTemp ReferenceWB;
        double ref_r, ref_g, ref_b;
	{
		// ...then we re-get the constants but now with auto which gives us better demosaicing and CA auto-correct
		// performance for strange white balance settings (such as UniWB)
		ri->get_colorsCoeff( ref_pre_mul, scale_mul, c_black, true);
		refwb_red = ri->get_pre_mul(0) / ref_pre_mul[0];
		refwb_green = ri->get_pre_mul(1) / ref_pre_mul[1];
		refwb_blue = ri->get_pre_mul(2) / ref_pre_mul[2];
                initialGain = max(scale_mul[0], scale_mul[1], scale_mul[2], scale_mul[3]) / min(scale_mul[0], scale_mul[1], scale_mul[2], scale_mul[3]);
		ref_r = imatrices.rgb_cam[0][0]*refwb_red + imatrices.rgb_cam[0][1]*refwb_green + imatrices.rgb_cam[0][2]*refwb_blue;
		ref_g = imatrices.rgb_cam[1][0]*refwb_red + imatrices.rgb_cam[1][1]*refwb_green + imatrices.rgb_cam[1][2]*refwb_blue;
		ref_b = imatrices.rgb_cam[2][0]*refwb_red + imatrices.rgb_cam[2][1]*refwb_green + imatrices.rgb_cam[2][2]*refwb_blue;
		ReferenceWB = ColorTemp (ref_r, ref_g, ref_b, 1.);
	}
	if (settings->verbose) {
		printf("Raw As Shot White balance: temp %f, tint %f\n", camera_wb.getTemp(), camera_wb.getGreen());
		printf("Raw Reference (auto) white balance: temp %f, tint %f, multipliers [%f %f %f | %f %f %f]\n", ReferenceWB.getTemp(), ReferenceWB.getGreen(), ref_r, ref_g, ref_b, refwb_red, refwb_blue, refwb_green);
	}

	/*{
	        // Test code: if you want to test a specific white balance
		ColorTemp d50wb = ColorTemp(5000.0, 1.0, 1.0, "Custom");
		double rm,gm,bm,r,g,b;
		d50wb.getMultipliers(r, g, b);
		camwb_red   = imatrices.cam_rgb[0][0]*r + imatrices.cam_rgb[0][1]*g + imatrices.cam_rgb[0][2]*b;
		camwb_green = imatrices.cam_rgb[1][0]*r + imatrices.cam_rgb[1][1]*g + imatrices.cam_rgb[1][2]*b;
		camwb_blue  = imatrices.cam_rgb[2][0]*r + imatrices.cam_rgb[2][1]*g + imatrices.cam_rgb[2][2]*b;
		double pre_mul[3], dmax = 0;
		pre_mul[0] = ri->get_pre_mul(0) / camwb_red;
		pre_mul[1] = ri->get_pre_mul(1) / camwb_green;
		pre_mul[2] = ri->get_pre_mul(2) / camwb_blue;
		for (int c = 0; c < 3; c++) {
			if (dmax < pre_mul[c])
				dmax = pre_mul[c];
                }
                for (int c = 0; c < 3; c++) {
			pre_mul[c] /= dmax;
                }
                camwb_red *= dmax;
                camwb_green *= dmax;
                camwb_blue *= dmax;
                for (int c = 0; c < 3; c++) {
			int sat = ri->get_white(c) - ri->get_cblack(c);
			scale_mul[c] = pre_mul[c] * 65535.0 / sat;
                }
                scale_mul[3] = pre_mul[1] * 65535.0 / (ri->get_white(3) - ri->get_cblack(3));
                initialGain = 1.0 / min(pre_mul[0], pre_mul[1], pre_mul[2]);
	}*/
            

    ri->set_prefilters();

    //Load complete Exif informations
    RawMetaDataLocation rml;
    rml.exifBase = ri->get_exifBase();
    rml.ciffBase = ri->get_ciffBase();
    rml.ciffLength = ri->get_ciffLen();
    idata = new ImageData (fname, &rml);

    green(W,H);
    red(W,H);
    blue(W,H);
    //hpmap = allocArray<char>(W, H);

    if (plistener) {
        plistener->setProgress (1.0);
    }
    plistener=NULL; // This must be reset, because only load() is called through progressConnector
    t2.set();
    if( settings->verbose )
       printf("Load %s: %d usec\n",fname.c_str(), t2.etime(t1));

    return 0; // OK!
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::preprocess  (const RAWParams &raw, const LensProfParams &lensProf, const CoarseTransformParams& coarse)
{
	MyTime t1,t2;
	t1.set();
	Glib::ustring newDF = raw.dark_frame;
	Glib::ustring makerstring = ((Glib::ustring)ri->get_maker()).uppercase();
	Glib::ustring modelstring = ((Glib::ustring)ri->get_model()).uppercase();
	RawImage *rid=NULL;
	if (!raw.df_autoselect) {
		if( !raw.dark_frame.empty())
			rid = dfm.searchDarkFrame( raw.dark_frame );
	} else {
		rid = dfm.searchDarkFrame( makerstring, modelstring, ri->get_ISOspeed(), ri->get_shutter(), ri->get_timestamp());
	}
	if( rid && settings->verbose){
		printf( "Subtracting Darkframe:%s\n",rid->get_filename().c_str());
	}
	//copyOriginalPixels(ri, rid);

	//FLATFIELD start
	Glib::ustring newFF = raw.ff_file;
	RawImage *rif=NULL;
	if (!raw.ff_AutoSelect) {
		if( !raw.ff_file.empty())
			rif = ffm.searchFlatField( raw.ff_file );
	} else {
		rif = ffm.searchFlatField( idata->getMake(), idata->getModel(),idata->getLens(),idata->getFocalLen(), idata->getFNumber(), idata->getDateTimeAsTS());
	}

    bool hasFlatField = (rif!=NULL);
	if( hasFlatField && settings->verbose) {
		printf( "Flat Field Correction:%s\n",rif->get_filename().c_str());
	}
	copyOriginalPixels(raw, ri, rid, rif);
	//FLATFIELD end
	
	
	
	PixelsMap bitmapBads(W,H);
	int totBP=0; // Hold count of bad pixels to correct

	// Always correct camera badpixels
	std::list<badPix> *bp = dfm.getBadPixels( ri->get_maker(), ri->get_model(), std::string("") );
	if( bp ){
		totBP+=bitmapBads.set( *bp );
		if( settings->verbose ){
			std::cout << "Correcting " << bp->size() << " pixels from .badpixels" << std::endl;
		}
	}

	// If darkframe selected, correct hotpixels found on darkframe
	bp = 0;
	if( raw.df_autoselect ){
		bp = dfm.getHotPixels( makerstring, modelstring, ri->get_ISOspeed(), ri->get_shutter(), ri->get_timestamp());
	}else if( !raw.dark_frame.empty() )
		bp = dfm.getHotPixels( raw.dark_frame );
	if(bp){
		totBP+=bitmapBads.set( *bp );
		if( settings->verbose && !bp->empty()){
			std::cout << "Correcting " << bp->size() << " hotpixels from darkframe" << std::endl;
		}
	}

    scaleColors( 0,0, W, H, raw);//+ + raw parameters for black level(raw.blackxx)

    // Correct vignetting of lens profile
    if (!hasFlatField && lensProf.useVign) {
        LCPProfile *pLCPProf=lcpStore->getProfile(lensProf.lcpFile);

        if (pLCPProf) {
            LCPMapper map(pLCPProf, idata->getFocalLen(), idata->getFocalLen35mm(), idata->getFocusDist(), idata->getFNumber(), true, false, W, H, coarse, -1);
        
            #pragma omp parallel for
            for (int y=0; y<H; y++) {
                for (int x=0; x<W; x++) {
                    if (rawData[y][x]>0) rawData[y][x] *= map.calcVignetteFac(x,y);
                }
            }
        }
    }

    defGain = 0.0;//log(initialGain) / log(2.0);

	if ( raw.hotdeadpix_filt>0 ) {
		if (plistener) {
			plistener->setProgressStr ("Hot/Dead Pixel Filter...");
			plistener->setProgress (0.0);
		}
		float varthresh = (20.0*((float)raw.hotdeadpix_thresh/100.0) + 1.0 ); 
		int nFound =findHotDeadPixel( bitmapBads, varthresh );
		totBP += nFound;
		if( settings->verbose && nFound>0){
			printf( "Correcting %d hot/dead pixels found inside image\n",nFound );
		}
	}
	if( totBP )
	   cfaCleanFromMap( bitmapBads );

    // check if it is an olympus E camera, if yes, compute G channel pre-compensation factors
    if ( ri->getSensorType()==ST_BAYER && (raw.bayersensor.greenthresh || (((idata->getMake().size()>=7 && idata->getMake().substr(0,7)=="OLYMPUS" && idata->getModel()[0]=='E') || (idata->getMake().size()>=9 && idata->getMake().substr(0,9)=="Panasonic")) && raw.bayersensor.method != RAWParams::BayerSensor::methodstring[ RAWParams::BayerSensor::vng4])) ) {
        // global correction
        int ng1=0, ng2=0, i=0;
        double avgg1=0., avgg2=0.;

#pragma omp parallel for default(shared) private(i) reduction(+: ng1, ng2, avgg1, avgg2)
        for (i=border; i<H-border; i++)
            for (int j=border; j<W-border; j++)
                if (ri->ISGREEN(i,j)) {
                    if (i&1) {
						avgg2 += rawData[i][j];
                        ng2++;
                    }
                    else {
                        avgg1 += rawData[i][j];
                        ng1++;
                    }
                }
        double corrg1 = ((double)avgg1/ng1 + (double)avgg2/ng2) / 2.0 / ((double)avgg1/ng1);
        double corrg2 = ((double)avgg1/ng1 + (double)avgg2/ng2) / 2.0 / ((double)avgg2/ng2);

#pragma omp parallel for default(shared)
        for (int i=border; i<H-border; i++)
            for (int j=border; j<W-border; j++)
                if (ri->ISGREEN(i,j)) {
                    float currData;
                    currData = (float)(rawData[i][j] * ((i&1) ? corrg2 : corrg1));
                    rawData[i][j] = (currData);
                }
	}

	if ( ri->getSensorType()==ST_BAYER && raw.bayersensor.greenthresh >0) {
		if (plistener) {
			plistener->setProgressStr ("Green equilibrate...");
			plistener->setProgress (0.0);
		}
		green_equilibrate(0.01*(raw.bayersensor.greenthresh));
    }

	
	if ( ri->getSensorType()==ST_BAYER && raw.bayersensor.linenoise >0 ) {
		if (plistener) {
			plistener->setProgressStr ("Line Denoise...");
			plistener->setProgress (0.0);
		}

		cfa_linedn(0.00002*(raw.bayersensor.linenoise));
	}
	
	if ( (raw.ca_autocorrect || fabs(raw.cared)>0.001 || fabs(raw.cablue)>0.001) && ri->getSensorType()!=ST_FUJI_XTRANS ) {  // Auto CA correction disabled for X-Trans, for now...
		if (plistener) {
			plistener->setProgressStr ("CA Auto Correction...");
			plistener->setProgress (0.0);
		}
		
		CA_correct_RT(raw.cared, raw.cablue);
	}
	
	if ( raw.expos !=1 ) processRawWhitepoint(raw.expos, raw.preser);
	
	if(dirpyrdenoiseExpComp == INFINITY) {
        LUTu aehist; int aehistcompr;
        double clip=0;
        int brightness, contrast, black, hlcompr, hlcomprthresh;
        getAutoExpHistogram (aehist, aehistcompr);
        ImProcFunctions::getAutoExp (aehist, aehistcompr, getDefGain(), clip, dirpyrdenoiseExpComp, brightness, contrast, black, hlcompr, hlcomprthresh);
	}
	
    t2.set();
    if( settings->verbose )
       printf("Preprocessing: %d usec\n", t2.etime(t1));
    return;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
void RawImageSource::demosaic(const RAWParams &raw)
{
	MyTime t1,t2;
	t1.set();

	if (ri->getSensorType()==ST_BAYER) {
		if ( raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::hphd] )
			   hphd_demosaic ();
		else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::vng4] )
			vng4_demosaic ();
		else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::ahd] )
			ahd_demosaic (0,0,W,H);
		else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::amaze] )
			amaze_demosaic_RT (0,0,W,H);
		else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::dcb] )
			dcb_demosaic(raw.bayersensor.dcb_iterations, raw.bayersensor.dcb_enhance);
		else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::eahd])
			eahd_demosaic ();
		else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::igv])
			igv_interpolate(W,H);
		else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::lmmse])
			lmmse_interpolate_omp(W,H,raw.bayersensor.lmmse_iterations);
		else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::fast] )
			fast_demosaic (0,0,W,H);
		else if (raw.bayersensor.method == RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::mono] )
			nodemosaic(true);
		else
			nodemosaic(false);
		
		//if (raw.all_enhance) refinement_lassus();

	} else if (ri->getSensorType()==ST_FUJI_XTRANS) {
		if (raw.xtranssensor.method == RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::fast] )
			fast_xtrans_interpolate();
		else if (raw.xtranssensor.method == RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::onePass])
			xtrans_interpolate(1,false);
		else if (raw.xtranssensor.method == RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::threePass] )
			xtrans_interpolate(3,true);
		else if(raw.xtranssensor.method == RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::mono] )
			nodemosaic(true);
		else
			nodemosaic(false);
    } else if (ri->get_colors() == 1) {
        // Monochrome
        nodemosaic(true);
    }
    t2.set();

    rgbSourceModified = false;
    if( settings->verbose ) {
        if (getSensorType() == ST_BAYER)
            printf("Demosaicing Bayer data: %s - %d usec\n",raw.bayersensor.method.c_str(), t2.etime(t1));
        else if (getSensorType() == ST_FUJI_XTRANS)
            printf("Demosaicing X-Trans data: %s - %d usec\n",raw.xtranssensor.method.c_str(), t2.etime(t1));
    }

}

void RawImageSource::flushRawData() {
    if(cache) {
        delete [] cache;
        cache = 0;
    }
    if (rawData) {
        rawData(0,0);
    }
}

void RawImageSource::flushRGB() {
    if (green) {
        green(0,0);
    }
    if (red) {
        red(0,0);
    }
    if (blue) {
        blue(0,0);
    }
}

void RawImageSource::HLRecovery_Global(ToneCurveParams hrp )
{
	//color propagation highlight recovery 
	if (hrp.hrenabled && hrp.method=="Color"){
		if (settings->verbose) printf ("Applying Highlight Recovery: Color propagation...\n");
		HLRecovery_inpaint (red,green,blue);
		rgbSourceModified = true;
	}
	else{
		rgbSourceModified = false;
	}
}


void RawImageSource::processFlatField(const RAWParams &raw, RawImage *riFlatFile, unsigned short black[4])
{
	float (*cfablur);
	cfablur = (float (*)) calloc (H*W, sizeof *cfablur);
	int BS = raw.ff_BlurRadius;
	BS += BS&1;
	
	//function call to cfabloxblur 
	if (raw.ff_BlurType == RAWParams::ff_BlurTypestring[RAWParams::v_ff])
		cfaboxblur(riFlatFile, cfablur, 2*BS, 0);
	else if (raw.ff_BlurType == RAWParams::ff_BlurTypestring[RAWParams::h_ff])
		cfaboxblur(riFlatFile, cfablur, 0, 2*BS);
	else if (raw.ff_BlurType == RAWParams::ff_BlurTypestring[RAWParams::vh_ff])
		//slightly more complicated blur if trying to correct both vertical and horizontal anomalies
		cfaboxblur(riFlatFile, cfablur, BS, BS);//first do area blur to correct vignette
	else //(raw.ff_BlurType == RAWParams::ff_BlurTypestring[RAWParams::area_ff])
		cfaboxblur(riFlatFile, cfablur, BS, BS);
	
	if(ri->getSensorType()==ST_BAYER) {
		float refcolor[2][2];
		//find center ave values by channel
		for (int m=0; m<2; m++)
			for (int n=0; n<2; n++) {
				int row = 2*(H>>2)+m;
				int col = 2*(W>>2)+n;
				int c  = FC(row, col);
				int c4 = ( c == 1 && !(row&1) ) ? 3 : c;
				refcolor[m][n] = max(0.0f,cfablur[row*W+col] - black[c4]);
			}

		float limitFactor = 1.f;
		if(raw.ff_AutoClipControl) {
			int clipControlGui = 0;
			for (int m=0; m<2; m++)
				for (int n=0; n<2; n++) {
					float maxval = 0.f;
					int c  = FC(m, n);
					int c4 = ( c == 1 && !(m&1) ) ? 3 : c;
#pragma omp parallel
{
					float maxvalthr = 0.f;
#pragma omp for
					for (int row = 0; row< H-m; row+=2) {
						for (int col = 0; col < W-n; col+=2) {
							float tempval = (rawData[row+m][col+n]-black[c4]) * ( refcolor[m][n]/max(1e-5f,cfablur[(row+m)*W+col+n]-black[c4]) );
							if(tempval > maxvalthr)
								maxvalthr = tempval;
						}
					}
#pragma omp critical
{

					if(maxvalthr>maxval)
						maxval = maxvalthr;
					
}
}
					// now we have the max value for the channel
					// if it clips, calculate factor to avoid clipping
					if(maxval + black[c4] >= ri->get_white(c4))
						limitFactor = min(limitFactor,ri->get_white(c4) / (maxval + black[c4]));
				}
			clipControlGui = (1.f - limitFactor) * 100.f;			// this value can be used to set the clip control slider in gui
		} else {
			limitFactor = max((float)(100 - raw.ff_clipControl)/100.f,0.01f);
		}
		for (int m=0; m<2; m++)
			for (int n=0; n<2; n++)
				refcolor[m][n] *= limitFactor;
	
			
		for (int m=0; m<2; m++)
			for (int n=0; n<2; n++) {
#pragma omp parallel
{
				int c  = FC(m, n);
				int c4 = ( c == 1 && !(m&1) ) ? 3 : c;
#pragma omp for
				for (int row = 0; row< H-m; row+=2) {
					for (int col = 0; col < W-n; col+=2) {
						float vignettecorr = ( refcolor[m][n]/max(1e-5f,cfablur[(row+m)*W+col+n]-black[c4]) );
						rawData[row+m][col+n] = (rawData[row+m][col+n]-black[c4]) * vignettecorr + black[c4]; 	
					}
				}
}
			}
	} else if(ri->getSensorType()==ST_FUJI_XTRANS) {
		float refcolor[3] = {0.f};
		int cCount[3] = {0};
		//find center ave values by channel
		for (int m=-3; m<3; m++)
			for (int n=-3; n<3; n++) {
				int row = 2*(H>>2)+m;
				int col = 2*(W>>2)+n;
				int c  = riFlatFile->XTRANSFC(row, col);
				refcolor[c] += max(0.0f,cfablur[row*W+col] - black[c]);
				cCount[c] ++;
			}
		for(int c=0;c<3;c++)
			refcolor[c] = refcolor[c] / cCount[c];

		float limitFactor = 1.f;

		if(raw.ff_AutoClipControl) {
			// determine maximum calculated value to avoid clipping
			int clipControlGui = 0;
			float maxval = 0.f;
			// xtrans files have only one black level actually, so we can simplify the code a bit
#pragma omp parallel
{
			float maxvalthr = 0.f;
#pragma omp for schedule(dynamic,16) nowait
			for (int row = 0; row< H; row++) {
				for (int col = 0; col < W; col++) {
					float tempval = (rawData[row][col]-black[0]) * ( refcolor[ri->XTRANSFC(row, col)]/max(1e-5f,cfablur[(row)*W+col]-black[0]) );
					if(tempval > maxvalthr)
						maxvalthr = tempval;
				}
			}
#pragma omp critical
{
			if(maxvalthr>maxval)
				maxval = maxvalthr;
}
}
			// there's only one white level for xtrans
			if(maxval + black[0] > ri->get_white(0)) {
				limitFactor = ri->get_white(0) / (maxval + black[0]);
				clipControlGui = (1.f - limitFactor) * 100.f;			// this value can be used to set the clip control slider in gui
			}
		} else { 
			limitFactor = max((float)(100 - raw.ff_clipControl)/100.f,0.01f);
		}


		for(int c=0;c<3;c++)
			refcolor[c] *= limitFactor;
		
#pragma omp parallel for
		for (int row = 0; row< H; row++) {
			for (int col = 0; col < W; col++) {
				int c  = ri->XTRANSFC(row, col);
				float vignettecorr = ( refcolor[c]/max(1e-5f,cfablur[(row)*W+col]-black[c]) );
				rawData[row][col] = (rawData[row][col]-black[c]) * vignettecorr + black[c]; 	
			}
		}
	}
	if (raw.ff_BlurType == RAWParams::ff_BlurTypestring[RAWParams::vh_ff]) {
		float (*cfablur1);
		cfablur1 = (float (*)) calloc (H*W, sizeof *cfablur1);
		float (*cfablur2);
		cfablur2 = (float (*)) calloc (H*W, sizeof *cfablur2);
		//slightly more complicated blur if trying to correct both vertical and horizontal anomalies
		cfaboxblur(riFlatFile, cfablur1, 0, 2*BS);//now do horizontal blur
		cfaboxblur(riFlatFile, cfablur2, 2*BS, 0);//now do vertical blur
		
		if(ri->getSensorType()==ST_BAYER) {
			for (int m=0; m<2; m++)
				for (int n=0; n<2; n++) {
#pragma omp parallel for
					for (int row = 0; row< H-m; row+=2) {
						int c  = FC(row, 0);
						int c4 = ( c == 1 && !(row&1) ) ? 3 : c;
						for (int col = 0; col < W-n; col+=2) {
							float hlinecorr = (max(1e-5f,cfablur[(row+m)*W+col+n]-black[c4])/max(1e-5f,cfablur1[(row+m)*W+col+n]-black[c4]) );
							float vlinecorr = (max(1e-5f,cfablur[(row+m)*W+col+n]-black[c4])/max(1e-5f,cfablur2[(row+m)*W+col+n]-black[c4]) );
							rawData[row+m][col+n] = ((rawData[row+m][col+n]-black[c4]) * hlinecorr * vlinecorr + black[c4]); 
						}
					}
				}
		} else if(ri->getSensorType()==ST_FUJI_XTRANS) {
#pragma omp parallel for
			for (int row = 0; row< H; row++) {
				for (int col = 0; col < W; col++) {
					int c  = ri->XTRANSFC(row, col);
					float hlinecorr = (max(1e-5f,cfablur[(row)*W+col]-black[c])/max(1e-5f,cfablur1[(row)*W+col]-black[c]) );
					float vlinecorr = (max(1e-5f,cfablur[(row)*W+col]-black[c])/max(1e-5f,cfablur2[(row)*W+col]-black[c]) );
					rawData[row][col] = ((rawData[row][col]-black[c]) * hlinecorr * vlinecorr + black[c]); 
				}
			}
			
		}
		free (cfablur1);
		free (cfablur2);
	}
	
	free (cfablur);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/* Copy original pixel data and
 * subtract dark frame (if present) from current image and apply flat field correction (if present)
 */
void RawImageSource::copyOriginalPixels(const RAWParams &raw, RawImage *src, RawImage *riDark, RawImage *riFlatFile )
{
	unsigned short black[4]={ri->get_cblack(0),ri->get_cblack(1),ri->get_cblack(2),ri->get_cblack(3)};

	if (ri->getSensorType()!=ST_NONE) {
		if (!rawData)
			rawData(W,H);
		if (riDark && W == riDark->get_width() && H == riDark->get_height()) { // This works also for xtrans-sensors, because black[0] to black[4] are equal for these
			for (int row = 0; row < H; row++) {
				for (int col = 0; col < W; col++) {
					int c  = FC(row, col);
					int c4 = ( c == 1 && !(row&1) ) ? 3 : c;
					rawData[row][col]	= max(src->data[row][col]+black[c4] - riDark->data[row][col], 0.0f);
				}
			}
		}else{
			for (int row = 0; row < H; row++) {
				for (int col = 0; col < W; col++) {
					rawData[row][col]	= src->data[row][col];
				}
			}
		}
		
		
		if (riFlatFile && W == riFlatFile->get_width() && H == riFlatFile->get_height()) {
			processFlatField(raw, riFlatFile, black);			
		}  // flatfield
	} else if (ri->get_colors() == 1) {
		// Monochrome
		if (!rawData) rawData(W,H);

		if (riDark && W == riDark->get_width() && H == riDark->get_height()) {
			for (int row = 0; row < H; row++) {
				for (int col = 0; col < W; col++) {
					rawData[row][col] = max(src->data[row][col]+black[0] - riDark->data[row][col], 0.0f);
				}
			}
		} else {
			for (int row = 0; row < H; row++) {
				for (int col = 0; col < W; col++) {
					rawData[row][col] = src->data[row][col];
				}
			}
		}
	} else {
        // No bayer pattern
        // TODO: Is there a flat field correction possible?
		if (!rawData) rawData(3*W,H);

		if (riDark && W == riDark->get_width() && H == riDark->get_height()) {
			for (int row = 0; row < H; row++) {
				for (int col = 0; col < W; col++) {
					int c  = FC(row, col);
					int c4 = ( c == 1 && !(row&1) ) ? 3 : c;
					rawData[row][3*col+0] = max(src->data[row][3*col+0]+black[c4] - riDark->data[row][3*col+0], 0.0f);
					rawData[row][3*col+1] = max(src->data[row][3*col+1]+black[c4] - riDark->data[row][3*col+1], 0.0f);
					rawData[row][3*col+2] = max(src->data[row][3*col+2]+black[c4] - riDark->data[row][3*col+2], 0.0f);
				}
			}
		} else {
			for (int row = 0; row < H; row++) {
				for (int col = 0; col < W; col++) {
					rawData[row][3*col+0] = src->data[row][3*col+0];
					rawData[row][3*col+1] = src->data[row][3*col+1];
					rawData[row][3*col+2] = src->data[row][3*col+2];
				}
			}
		}
	}
}

SSEFUNCTION void RawImageSource::cfaboxblur(RawImage *riFlatFile, float* cfablur, const int boxH, const int boxW ) {

	float (*cfatmp);
	cfatmp = (float (*)) calloc (H*W, sizeof *cfatmp);
//	const float hotdeadthresh = 0.5;

#pragma omp parallel
{
#pragma omp for
	for (int i=0; i<H; i++) {
		int iprev,inext,jprev,jnext;
		int p[5],temp, median;
		if (i<2) {iprev=i+2;} else {iprev=i-2;}
		if (i>H-3) {inext=i-2;} else {inext=i+2;}
		for (int j=0; j<W; j++) {
			if (j<2) {jprev=j+2;} else {jprev=j-2;}
			if (j>W-3) {jnext=j-2;} else {jnext=j+2;}
			//med3x3(riFlatFile->data[iprev][jprev], riFlatFile->data[iprev][j], riFlatFile->data[iprev][jnext],
			//	   riFlatFile->data[i][jprev], riFlatFile->data[i][j], riFlatFile->data[i][jnext],
			//	   riFlatFile->data[inext][jprev], riFlatFile->data[inext][j], riFlatFile->data[inext][jnext], cfatmp[i*W+j]);
			med5(riFlatFile->data[iprev][j], riFlatFile->data[i][jprev],riFlatFile->data[i][j],
				 riFlatFile->data[i][jnext], riFlatFile->data[inext][j],median);
				 
//			if (riFlatFile->data[i][j]>hotdeadthresh*median || median>hotdeadthresh*riFlatFile->data[i][j]) {
			if (((int)riFlatFile->data[i][j]<<1)>median || (median<<1) > riFlatFile->data[i][j]) {
				cfatmp[i*W+j] = median;
			} else {
				cfatmp[i*W+j] = riFlatFile->data[i][j];
			}

		}
	}
	
	//box blur cfa image; box size = BS
	//horizontal blur
#pragma omp for
	for (int row = 0; row < H; row++) {
		int len = boxW/2 + 1;
		cfatmp[row*W+0] = cfatmp[row*W+0]/len;
		cfatmp[row*W+1] = cfatmp[row*W+1]/len;
		for (int j=2; j<=boxW; j+=2) {
			cfatmp[row*W+0] += cfatmp[row*W+j]/len;
			cfatmp[row*W+1] += cfatmp[row*W+j+1]/len;
		}
		for (int col=2; col<=boxW; col+=2) {
			cfatmp[row*W+col] = (cfatmp[row*W+col-2]*len + cfatmp[row*W+boxW+col])/(len+1);
			cfatmp[row*W+col+1] = (cfatmp[row*W+col-1]*len + cfatmp[row*W+boxW+col+1])/(len+1);
			len ++;
		}
		for (int col = boxW+2; col < W-boxW; col++) {
			cfatmp[row*W+col] = cfatmp[row*W+col-2] + (cfatmp[row*W+boxW+col]-cfatmp[row*W+col-boxW-2])/len;
		}
		for (int col=W-boxW; col<W; col+=2) {
			cfatmp[row*W+col] = (cfatmp[row*W+col-2]*len - cfatmp[row*W+col-boxW-2])/(len-1);
			if (col+1<W) {
				cfatmp[row*W+col+1] = (cfatmp[row*W+col-1]*len - cfatmp[row*W+col-boxW-1])/(len-1);
			}
			len --;
		}
	}

	//vertical blur
#ifdef __SSE2__
	__m128	leninitv = _mm_set1_ps( (float)((int)(boxH/2 + 1)));
	__m128 	onev = _mm_set1_ps( 1.0f );
	__m128	temp1v,temp2v,lenv,lenp1v,lenm1v;
	int row;
#pragma omp for
	for (int col = 0; col < W-3; col+=4) {
		lenv = leninitv;
		temp1v = LVFU(cfatmp[0*W+col]) / lenv;
		temp2v = LVFU(cfatmp[1*W+col]) / lenv;
		
		for (int i=2; i<boxH+2; i+=2) {
			temp1v += LVFU(cfatmp[i*W+col]) / lenv;
			temp2v += LVFU(cfatmp[(i+1)*W+col]) / lenv;
		}
		_mm_storeu_ps(&cfablur[0*W+col], temp1v);
		_mm_storeu_ps(&cfablur[1*W+col], temp2v);
		for (row=2; row<boxH+2; row+=2) {
			lenp1v = lenv + onev;
			temp1v = (temp1v * lenv + LVFU(cfatmp[(row+boxH)*W+col])) / lenp1v;
			temp2v = (temp2v * lenv + LVFU(cfatmp[(row+boxH+1)*W+col])) / lenp1v;
			_mm_storeu_ps( &cfablur[row*W+col], temp1v);
			_mm_storeu_ps( &cfablur[(row+1)*W+col], temp2v);
			lenv = lenp1v;
		}
		for (; row < H-boxH-1; row+=2) {
			temp1v = temp1v + (LVFU(cfatmp[(row+boxH)*W+col]) - LVFU(cfatmp[(row-boxH-2)*W+col]))/lenv;
			temp2v = temp2v + (LVFU(cfatmp[(row+1+boxH)*W+col]) - LVFU(cfatmp[(row+1-boxH-2)*W+col]))/lenv;
			_mm_storeu_ps(&cfablur[row*W+col], temp1v);
			_mm_storeu_ps(&cfablur[(row+1)*W+col], temp2v);
		}
		for(; row < H-boxH; row++) {
			temp1v = temp1v + (LVFU(cfatmp[(row+boxH)*W+col]) - LVFU(cfatmp[(row-boxH-2)*W+col]))/lenv;
			_mm_storeu_ps(&cfablur[row*W+col], temp1v);
			__m128 swapv = temp1v;
			temp1v = temp2v;
			temp2v = swapv;
			
		}
		for (; row<H-1; row+=2) {
			lenm1v = lenv - onev;
			temp1v = (temp1v * lenv - LVFU(cfatmp[(row-boxH-2)*W+col])) / lenm1v;
			temp2v = (temp2v * lenv - LVFU(cfatmp[(row-boxH-1)*W+col])) / lenm1v;
			_mm_storeu_ps(&cfablur[row*W+col], temp1v);
			_mm_storeu_ps(&cfablur[(row+1)*W+col], temp2v);
			lenv = lenm1v;
		}
		for(; row < H; row++) {
			lenm1v = lenv - onev;
			temp1v = (temp1v * lenv - LVFU(cfatmp[(row-boxH-2)*W+col])) / lenm1v;
			_mm_storeu_ps(&cfablur[(row)*W+col], temp1v);
		}

	}
	for (int col = W-(W%4); col < W; col++) {
		int len = boxH/2 + 1;
		cfablur[0*W+col] = cfatmp[0*W+col]/len;
		cfablur[1*W+col] = cfatmp[1*W+col]/len;
		for (int i=2; i<boxH+2; i+=2) {
			cfablur[0*W+col] += cfatmp[i*W+col]/len;
			cfablur[1*W+col] += cfatmp[(i+1)*W+col]/len;
		}
		for (int row=2; row<boxH+2; row+=2) {
			cfablur[row*W+col] = (cfablur[(row-2)*W+col]*len + cfatmp[(row+boxH)*W+col])/(len+1);
			cfablur[(row+1)*W+col] = (cfablur[(row-1)*W+col]*len + cfatmp[(row+boxH+1)*W+col])/(len+1);
			len ++;
		}
		for (int row = boxH+2; row < H-boxH; row++) {
			cfablur[row*W+col] = cfablur[(row-2)*W+col] + (cfatmp[(row+boxH)*W+col] - cfatmp[(row-boxH-2)*W+col])/len;
		}
		for (int row=H-boxH; row<H; row+=2) {
			cfablur[row*W+col] = (cfablur[(row-2)*W+col]*len - cfatmp[(row-boxH-2)*W+col])/(len-1);
			if (row+1<H) 
				cfablur[(row+1)*W+col] = (cfablur[(row-1)*W+col]*len - cfatmp[(row-boxH-1)*W+col])/(len-1);
			len --;
		}
	}

#else
#pragma omp for
	for (int col = 0; col < W; col++) {
		int len = boxH/2 + 1;
		cfablur[0*W+col] = cfatmp[0*W+col]/len;
		cfablur[1*W+col] = cfatmp[1*W+col]/len;
		for (int i=2; i<boxH+2; i+=2) {
			cfablur[0*W+col] += cfatmp[i*W+col]/len;
			cfablur[1*W+col] += cfatmp[(i+1)*W+col]/len;
		}
		for (int row=2; row<boxH+2; row+=2) {
			cfablur[row*W+col] = (cfablur[(row-2)*W+col]*len + cfatmp[(row+boxH)*W+col])/(len+1);
			cfablur[(row+1)*W+col] = (cfablur[(row-1)*W+col]*len + cfatmp[(row+boxH+1)*W+col])/(len+1);
			len ++;
		}
		for (int row = boxH+2; row < H-boxH; row++) {
			cfablur[row*W+col] = cfablur[(row-2)*W+col] + (cfatmp[(row+boxH)*W+col] - cfatmp[(row-boxH-2)*W+col])/len;
		}
		for (int row=H-boxH; row<H; row+=2) {
			cfablur[row*W+col] = (cfablur[(row-2)*W+col]*len - cfatmp[(row-boxH-2)*W+col])/(len-1);
			if (row+1<H) 
				cfablur[(row+1)*W+col] = (cfablur[(row-1)*W+col]*len - cfatmp[(row-boxH-1)*W+col])/(len-1);
			len --;
		}
	}
#endif
}
	free (cfatmp);
}
	

// Scale original pixels into the range 0 65535 using black offsets and multipliers 
void RawImageSource::scaleColors(int winx,int winy,int winw,int winh, const RAWParams &raw)
{
	chmax[0]=chmax[1]=chmax[2]=chmax[3]=0;//channel maxima
	float black_lev[4];//black level

	//adjust black level  (eg Canon)
	bool isMono = false;
	if (getSensorType()==ST_BAYER) {

		black_lev[0]=raw.bayersensor.black1;//R
		black_lev[1]=raw.bayersensor.black0;//G1
		black_lev[2]=raw.bayersensor.black2;//B
		black_lev[3]=raw.bayersensor.black3;//G2

		isMono = RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::mono] == raw.bayersensor.method;
	}
	else if (getSensorType()==ST_FUJI_XTRANS) {

		black_lev[0]=raw.xtranssensor.blackred;//R
		black_lev[1]=raw.xtranssensor.blackgreen;//G1
		black_lev[2]=raw.xtranssensor.blackblue;//B
		black_lev[3]=raw.xtranssensor.blackgreen;//G2  (set, only used with a Bayer filter)

		isMono = RAWParams::XTransSensor::methodstring[RAWParams::XTransSensor::mono] == raw.xtranssensor.method;
	}

	for(int i=0; i<4 ;i++) cblacksom[i] = max( c_black[i]+black_lev[i], 0.0f ); // adjust black level
        initialGain = calculate_scale_mul(scale_mul, ref_pre_mul, c_white, cblacksom, isMono, ri->get_colors()); // recalculate scale colors with adjusted levels
        //fprintf(stderr, "recalc: %f [%f %f %f %f]\n", initialGain, scale_mul[0], scale_mul[1], scale_mul[2], scale_mul[3]);
        
	// this seems strange, but it works

		// scale image colors

	if( ri->getSensorType()==ST_BAYER){
#pragma omp parallel
{
		float tmpchmax[3];
		tmpchmax[0] = tmpchmax[1] = tmpchmax[2] = 0.0f;

#pragma omp for nowait
		for (int row = winy; row < winy+winh; row ++){
			for (int col = winx; col < winx+winw; col++) {
				float val = rawData[row][col];
				int c  = FC(row, col);                        // three colors,  0=R, 1=G,  2=B
				int c4 = ( c == 1 && !(row&1) ) ? 3 : c;      // four  colors,  0=R, 1=G1, 2=B, 3=G2
				val-=cblacksom[c4];
				val*=scale_mul[c4];
				
				rawData[row][col] = (val);
				tmpchmax[c] = max(tmpchmax[c],val);
			}
		}
#pragma omp critical
{
		chmax[0] = max(tmpchmax[0],chmax[0]);
		chmax[1] = max(tmpchmax[1],chmax[1]);
		chmax[2] = max(tmpchmax[2],chmax[2]);
}
}
	} else if ( ri->get_colors() == 1 ) {
#pragma omp parallel
{
		float tmpchmax = 0.0f;

#pragma omp for nowait
		for (int row = winy; row < winy+winh; row ++){
			for (int col = winx; col < winx+winw; col++) {
				float val = rawData[row][col];
				val -= cblacksom[0];
				val *= scale_mul[0];
				rawData[row][col] = (val);
				tmpchmax = max(tmpchmax,val);
			}
		}
#pragma omp critical
{
		chmax[0] = chmax[1] = chmax[2] = chmax[3] = max(tmpchmax,chmax[0]);
}
}
	} else if(ri->getSensorType()==ST_FUJI_XTRANS) {
#pragma omp parallel
{
		float tmpchmax[3];
		tmpchmax[0] = tmpchmax[1] = tmpchmax[2] = 0.0f;

#pragma omp for nowait
		for (int row = winy; row < winy+winh; row ++){
			for (int col = winx; col < winx+winw; col++) {
				float val = rawData[row][col];
				int c = ri->XTRANSFC(row, col);
				val-=cblacksom[c];
				val*=scale_mul[c];
				
				rawData[row][col] = (val);
				tmpchmax[c] = max(tmpchmax[c],val);
			}
		}
#pragma omp critical
{
		chmax[0] = max(tmpchmax[0],chmax[0]);
		chmax[1] = max(tmpchmax[1],chmax[1]);
		chmax[2] = max(tmpchmax[2],chmax[2]);
}
}
	} else {
#pragma omp parallel
{
		float tmpchmax[3];
		tmpchmax[0] = tmpchmax[1] = tmpchmax[2] = 0.0f;

#pragma omp for nowait
		for (int row = winy; row < winy+winh; row ++){
			for (int col = winx; col < winx+winw; col++) {
				for (int c=0; c<3; c++) {                     // three colors,  0=R, 1=G,  2=B
					float val = rawData[row][3*col+c];
					val -= cblacksom[c];
					val *= scale_mul[c];
					rawData[row][3*col+c] = (val);
					tmpchmax[c] = max(tmpchmax[c],val);
				}
			}
		}
#pragma omp critical
{
		chmax[0] = max(tmpchmax[0],chmax[0]);
		chmax[1] = max(tmpchmax[1],chmax[1]);
		chmax[2] = max(tmpchmax[2],chmax[2]);
}
}
		chmax[3]=chmax[1];
	}

}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int RawImageSource::defTransform (int tran) {

    int deg = ri->get_rotateDegree();
    if ((tran & TR_ROT) == TR_R180)
        deg += 180;
    else if ((tran & TR_ROT) == TR_R90)
        deg += 90;
    else if ((tran & TR_ROT) == TR_R270)
        deg += 270;
    deg %= 360;

    int ret = 0;
    if (deg==90)
        ret |= TR_R90;
    else if (deg==180)
        ret |= TR_R180;
    else if (deg==270)
        ret |= TR_R270;
    if (tran & TR_HFLIP)
        ret |= TR_HFLIP;
    if (tran & TR_VFLIP)
        ret |= TR_VFLIP;
    return ret;
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Thread called part
void RawImageSource::processFalseColorCorrectionThread  (Imagefloat* im, int row_from, int row_to) {
 
  int W = im->width;

  array2D<float> rbconv_Y (W,3);
  array2D<float> rbconv_I (W,3);
  array2D<float> rbconv_Q (W,3);
  array2D<float> rbout_I (W,3);
  array2D<float> rbout_Q (W,3);

  float* row_I = new float[W];
  float* row_Q = new float[W];

  float* pre1_I = new float[3];
  float* pre2_I = new float[3];
  float* post1_I = new float[3];
  float* post2_I = new float[3];
  float middle_I[6];
  float* pre1_Q = new float[3];
  float* pre2_Q = new float[3];
  float* post1_Q = new float[3];
  float* post2_Q = new float[3];
  float middle_Q[6];
  float* tmp;

  int px=(row_from-1)%3, cx=row_from%3, nx=0;
  
    convert_row_to_YIQ (im->r(row_from-1), im->g(row_from-1), im->b(row_from-1), rbconv_Y[px], rbconv_I[px], rbconv_Q[px], W);
    convert_row_to_YIQ (im->r(row_from), im->g(row_from), im->b(row_from), rbconv_Y[cx], rbconv_I[cx], rbconv_Q[cx], W);

    for (int j=0; j<W; j++) {
      rbout_I[px][j] = rbconv_I[px][j];
      rbout_Q[px][j] = rbconv_Q[px][j];
    }

    for (int i=row_from; i<row_to; i++) {
       
      px = (i-1)%3;
      cx = i%3;
      nx = (i+1)%3;

      convert_row_to_YIQ (im->r(i+1), im->g(i+1), im->b(i+1), rbconv_Y[nx], rbconv_I[nx], rbconv_Q[nx], W);

      SORT3(rbconv_I[px][0],rbconv_I[cx][0],rbconv_I[nx][0],pre1_I[0],pre1_I[1],pre1_I[2]);
      SORT3(rbconv_I[px][1],rbconv_I[cx][1],rbconv_I[nx][1],pre2_I[0],pre2_I[1],pre2_I[2]);
      SORT3(rbconv_Q[px][0],rbconv_Q[cx][0],rbconv_Q[nx][0],pre1_Q[0],pre1_Q[1],pre1_Q[2]);
      SORT3(rbconv_Q[px][1],rbconv_Q[cx][1],rbconv_Q[nx][1],pre2_Q[0],pre2_Q[1],pre2_Q[2]);

      // median I channel
      for (int j=1; j<W-2; j+=2) {
        SORT3(rbconv_I[px][j+1],rbconv_I[cx][j+1],rbconv_I[nx][j+1],post1_I[0],post1_I[1],post1_I[2]);
        SORT3(rbconv_I[px][j+2],rbconv_I[cx][j+2],rbconv_I[nx][j+2],post2_I[0],post2_I[1],post2_I[2]);
        MERGESORT(pre2_I[0],pre2_I[1],pre2_I[2],post1_I[0],post1_I[1],post1_I[2],middle_I[0],middle_I[1],middle_I[2],middle_I[3],middle_I[4],middle_I[5]);
        MEDIAN7(pre1_I[0],pre1_I[1],pre1_I[2],middle_I[1],middle_I[2],middle_I[3],middle_I[4],rbout_I[cx][j]);
        MEDIAN7(post2_I[0],post2_I[1],post2_I[2],middle_I[1],middle_I[2],middle_I[3],middle_I[4],rbout_I[cx][j+1]);
        tmp = pre1_I;
        pre1_I = post1_I;
        post1_I = tmp;
        tmp = pre2_I;
        pre2_I = post2_I;
        post2_I = tmp;

      }
      // median Q channel
      for (int j=1; j<W-2; j+=2) {
        SORT3(rbconv_Q[px][j+1],rbconv_Q[cx][j+1],rbconv_Q[nx][j+1],post1_Q[0],post1_Q[1],post1_Q[2]);
        SORT3(rbconv_Q[px][j+2],rbconv_Q[cx][j+2],rbconv_Q[nx][j+2],post2_Q[0],post2_Q[1],post2_Q[2]);
        MERGESORT(pre2_Q[0],pre2_Q[1],pre2_Q[2],post1_Q[0],post1_Q[1],post1_Q[2],middle_Q[0],middle_Q[1],middle_Q[2],middle_Q[3],middle_Q[4],middle_Q[5]);
        MEDIAN7(pre1_Q[0],pre1_Q[1],pre1_Q[2],middle_Q[1],middle_Q[2],middle_Q[3],middle_Q[4],rbout_Q[cx][j]);
        MEDIAN7(post2_Q[0],post2_Q[1],post2_Q[2],middle_Q[1],middle_Q[2],middle_Q[3],middle_Q[4],rbout_Q[cx][j+1]);
        tmp = pre1_Q;
        pre1_Q = post1_Q;
        post1_Q = tmp;
        tmp = pre2_Q;
        pre2_Q = post2_Q;
        post2_Q = tmp;
      }
      // fill first and last element in rbout
      rbout_I[cx][0] = rbconv_I[cx][0];
      rbout_I[cx][W-1] = rbconv_I[cx][W-1];
      rbout_I[cx][W-2] = rbconv_I[cx][W-2];
      rbout_Q[cx][0] = rbconv_Q[cx][0];
      rbout_Q[cx][W-1] = rbconv_Q[cx][W-1];
      rbout_Q[cx][W-2] = rbconv_Q[cx][W-2];

      // blur i-1th row
      if (i>row_from) {
        for (int j=1; j<W-1; j++) {
          row_I[j] = (rbout_I[px][j-1]+rbout_I[px][j]+rbout_I[px][j+1]+rbout_I[cx][j-1]+rbout_I[cx][j]+rbout_I[cx][j+1]+rbout_I[nx][j-1]+rbout_I[nx][j]+rbout_I[nx][j+1])/9;
          row_Q[j] = (rbout_Q[px][j-1]+rbout_Q[px][j]+rbout_Q[px][j+1]+rbout_Q[cx][j-1]+rbout_Q[cx][j]+rbout_Q[cx][j+1]+rbout_Q[nx][j-1]+rbout_Q[nx][j]+rbout_Q[nx][j+1])/9;
        }
        row_I[0] = rbout_I[px][0];
        row_Q[0] = rbout_Q[px][0];
        row_I[W-1] = rbout_I[px][W-1];
        row_Q[W-1] = rbout_Q[px][W-1];
        convert_row_to_RGB (im->r(i-1), im->g(i-1), im->b(i-1), rbconv_Y[px], row_I, row_Q, W);
      }
    }
    // blur last 3 row and finalize H-1th row
    for (int j=1; j<W-1; j++) {
      row_I[j] = (rbout_I[px][j-1]+rbout_I[px][j]+rbout_I[px][j+1]+rbout_I[cx][j-1]+rbout_I[cx][j]+rbout_I[cx][j+1]+rbconv_I[nx][j-1]+rbconv_I[nx][j]+rbconv_I[nx][j+1])/9;
      row_Q[j] = (rbout_Q[px][j-1]+rbout_Q[px][j]+rbout_Q[px][j+1]+rbout_Q[cx][j-1]+rbout_Q[cx][j]+rbout_Q[cx][j+1]+rbconv_Q[nx][j-1]+rbconv_Q[nx][j]+rbconv_Q[nx][j+1])/9;
    }
    row_I[0] = rbout_I[cx][0];
    row_Q[0] = rbout_Q[cx][0];
    row_I[W-1] = rbout_I[cx][W-1];
    row_Q[W-1] = rbout_Q[cx][W-1];
    convert_row_to_RGB (im->r(row_to-1), im->g(row_to-1), im->b(row_to-1), rbconv_Y[cx], row_I, row_Q, W);

  delete [] row_I;
  delete [] row_Q;
  delete [] pre1_I;
  delete [] pre2_I;
  delete [] post1_I;
  delete [] post2_I;
  delete [] pre1_Q;
  delete [] pre2_Q;
  delete [] post1_Q;
  delete [] post2_Q;
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// correction_YIQ_LQ
void RawImageSource::processFalseColorCorrection  (Imagefloat* im, int steps) {

    if (im->height<4)
        return;

    for (int t=0; t<steps; t++) {
#ifdef _OPENMP
	    #pragma omp parallel
    	{
    		int tid = omp_get_thread_num();
    		int nthreads = omp_get_num_threads();
    		int blk = (im->height-2)/nthreads;

    		if (tid<nthreads-1)
    			processFalseColorCorrectionThread (im, 1 + tid*blk, 1 + (tid+1)*blk);
    		else
    			processFalseColorCorrectionThread (im, 1 + tid*blk, im->height - 1);
    	}
#else
    	processFalseColorCorrectionThread (im, 1 , im->height - 1);
#endif
    }
}
	
// Some camera input profiles need gamma preprocessing
// gamma is applied before the CMS, correct line fac=lineFac*rawPixel+LineSum after the CMS
void RawImageSource::getProfilePreprocParams(cmsHPROFILE in, float& gammaFac, float& lineFac, float& lineSum) {
    gammaFac=0; lineFac=1; lineSum=0;

    char copyright[256];
    copyright[0]=0;

    if (cmsGetProfileInfoASCII(in, cmsInfoCopyright, cmsNoLanguage, cmsNoCountry, copyright, 256)>0) {
        if (strstr(copyright,"Phase One")!=NULL)
            gammaFac=0.55556;  // 1.8
        else if (strstr(copyright,"Nikon Corporation")!=NULL) {
            gammaFac=0.5; lineFac=-0.4; lineSum=1.35;  // determined in reverse by measuring NX an RT developed colorchecker PNGs
        }
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

static void
lab2ProphotoRgbD50(float L, float A, float B, float& r, float& g, float& b)
{
    float X;
    float Y;
    float Z;
#define CLIP01(a) ((a)>0?((a)<1?(a):1):0)
    { // convert from Lab to XYZ
        float x, y, z, fx, fy, fz;

        fy = (L + 16.0f)/116.0f;
        fx = A/500.0f + fy;
        fz = fy - B/200.0f;

        if (fy > 24.0f/116.0f) {
            y = fy*fy*fy;
        } else {
            y = (fy - 16.0f/116.0f)/7.787036979f;
        }
        if (fx > 24.0f/116.0f) {
            x = fx*fx*fx;
        } else {
            x = (fx - 16.0/116.0)/7.787036979f;
        }
        if (fz > 24.0f/116.0f) {
            z = fz*fz*fz;
        } else {
            z = (fz - 16.0f/116.0f)/7.787036979f;
        }
        //0.9642, 1.0000, 0.8249 D50
        X = x * 0.9642;
        Y = y;
        Z = z * 0.8249;
    }
    r = prophoto_xyz[0][0]*X + prophoto_xyz[0][1]*Y + prophoto_xyz[0][2]*Z;
    g = prophoto_xyz[1][0]*X + prophoto_xyz[1][1]*Y + prophoto_xyz[1][2]*Z;
    b = prophoto_xyz[2][0]*X + prophoto_xyz[2][1]*Y + prophoto_xyz[2][2]*Z;
    r = CLIP01(r);
    g = CLIP01(g);
    b = CLIP01(b);
}

// Converts raw image including ICC input profile to working space - floating point version
void RawImageSource::colorSpaceConversion_ (Imagefloat* im, ColorManagementParams &cmp, ColorTemp &wb, double pre_mul[3], const RAWParams &raw, cmsHPROFILE embedded, cmsHPROFILE camprofile, double camMatrix[3][3], const std::string &camName) {

//    MyTime t1, t2, t3;
//    t1.set ();
    cmsHPROFILE in;
    DCPProfile *dcpProf;

    if (!findInputProfile(cmp.input, embedded, camName, &dcpProf, in)) {
        return;
    }

    if (dcpProf!=NULL) {
        // DCP processing
        dcpProf->Apply(im, cmp.dcpIlluminant, cmp.working, wb, pre_mul, camMatrix, raw.expos, cmp.toneCurve);
        return;
    }

    if (in==NULL) {
        // use default camprofile, supplied by dcraw
        // in this case we avoid using the slllllooooooowwww lcms

        // Calculate matrix for direct conversion raw>working space
        TMatrix work = iccStore->workingSpaceInverseMatrix (cmp.working);
        double mat[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (int i=0; i<3; i++)
            for (int j=0; j<3; j++)
                for (int k=0; k<3; k++)
                    mat[i][j] += work[i][k] * camMatrix[k][j]; // rgb_xyz * imatrices.xyz_cam


        #pragma omp parallel for
        for (int i=0; i<im->height; i++)
            for (int j=0; j<im->width; j++) {

                float newr = mat[0][0]*im->r(i,j) + mat[0][1]*im->g(i,j) + mat[0][2]*im->b(i,j);
                float newg = mat[1][0]*im->r(i,j) + mat[1][1]*im->g(i,j) + mat[1][2]*im->b(i,j);
                float newb = mat[2][0]*im->r(i,j) + mat[2][1]*im->g(i,j) + mat[2][2]*im->b(i,j);

                im->r(i,j) = newr;
                im->g(i,j) = newg;
                im->b(i,j) = newb;

            }
    } else {
        const bool working_space_is_prophoto = (cmp.working == "ProPhoto");

        // use supplied input profile

        /*
          The goal here is to in addition to user-made custom ICC profiles also support profiles
          supplied with other popular raw converters. As curves affect color rendering and
          different raw converters deal with them differently (and few if any is as flexible
          as RawTherapee) we cannot really expect to get the *exact* same color rendering here.
          However we try hard to make the best out of it.

          Third-party input profiles that contain a LUT (usually A2B0 tag) often needs some preprocessing,
          as ICC LUTs are not really designed for dealing with linear camera data. Generally one
          must apply some sort of curve to get efficient use of the LUTs. Unfortunately how you
          should preprocess is not standardized so there are almost as many ways as there are
          software makers, and for each one we have to reverse engineer to find out how it has
          been done. (The ICC files made for RT has linear LUTs)

          ICC profiles which only contain the <r,g,b>XYZ tags (ie only a color matrix) should
          (hopefully) not require any pre-processing.

          Some LUT ICC profiles apply a contrast curve and desaturate highlights (to give a "film-like"
          behavior. These will generally work with RawTherapee, but will not produce good results when
          you enable highlight recovery/reconstruction, as that data is added linearly on top of the
          original range. RawTherapee works best with linear ICC profiles.
        */

        enum camera_icc_type {
            CAMERA_ICC_TYPE_GENERIC, // Generic, no special pre-processing required, RTs own is this way
            CAMERA_ICC_TYPE_PHASE_ONE, // Capture One profiles
            CAMERA_ICC_TYPE_LEAF, // Leaf profiles, former Leaf Capture now in Capture One, made for Leaf digital backs
            CAMERA_ICC_TYPE_NIKON // Nikon NX profiles
        } camera_icc_type = CAMERA_ICC_TYPE_GENERIC;

        float leaf_prophoto_mat[3][3];
        { // identify ICC type
            char copyright[256] = "";
            char description[256] = "";

            cmsGetProfileInfoASCII(in, cmsInfoCopyright, cmsNoLanguage, cmsNoCountry, copyright, 256);
            cmsGetProfileInfoASCII(in, cmsInfoDescription, cmsNoLanguage, cmsNoCountry, description, 256);
            camera_icc_type = CAMERA_ICC_TYPE_GENERIC;
            // Note: order the identification with the most detailed matching first since the more general ones may also match the more detailed
            if ((strstr(copyright, "Leaf") != NULL ||
                 strstr(copyright, "Phase One A/S") != NULL ||
                 strstr(copyright, "Kodak") != NULL ||
                 strstr(copyright, "Creo") != NULL) &&
                (strstr(description,"LF2 ") == description ||
                 strstr(description,"LF3 ") == description ||
                 strstr(description,"LeafLF2") == description ||
                 strstr(description,"LeafLF3") == description ||
                 strstr(description,"MamiyaLF2") == description ||
                 strstr(description,"MamiyaLF3") == description))
            {
                camera_icc_type = CAMERA_ICC_TYPE_LEAF;
            } else if (strstr(copyright, "Phase One A/S") != NULL) {
                camera_icc_type = CAMERA_ICC_TYPE_PHASE_ONE;
            } else if (strstr(copyright,"Nikon Corporation")!=NULL) {
                camera_icc_type = CAMERA_ICC_TYPE_NIKON;
            }
        }

        // Initialize transform
        cmsHTRANSFORM hTransform;
        cmsHPROFILE prophoto = iccStore->workingSpace("ProPhoto"); // We always use Prophoto to apply the ICC profile to minimize problems with clipping in LUT conversion.
        bool transform_via_pcs_lab = false;
        bool separate_pcs_lab_highlights = false;
        lcmsMutex->lock ();
        switch (camera_icc_type) {
        case CAMERA_ICC_TYPE_PHASE_ONE:
        case CAMERA_ICC_TYPE_LEAF: {
            // These profiles have a RGB to Lab cLUT, gives gamma 1.8 output, and expects a "film-like" curve on input
            transform_via_pcs_lab = true;
            separate_pcs_lab_highlights = true;
            // We transform to Lab because we can and that we avoid getting an unnecessary unmatched gamma conversion which we would need to revert.
            hTransform = cmsCreateTransform (in, TYPE_RGB_FLT, NULL, TYPE_Lab_FLT, INTENT_RELATIVE_COLORIMETRIC, cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE );
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    leaf_prophoto_mat[i][j] = 0;
                    for (int k=0; k<3; k++)
                        leaf_prophoto_mat[i][j] += prophoto_xyz[i][k] * camMatrix[k][j];
                }
            }
            break;
        }
        case CAMERA_ICC_TYPE_NIKON:
        case CAMERA_ICC_TYPE_GENERIC:
        default:
            hTransform = cmsCreateTransform (in, TYPE_RGB_FLT, prophoto, TYPE_RGB_FLT, INTENT_RELATIVE_COLORIMETRIC, cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE );  // NOCACHE is important for thread safety
            break;
        }

        lcmsMutex->unlock ();
        if (hTransform == NULL) {
            // Fallback: create transform from camera profile. Should not happen normally.
            lcmsMutex->lock ();
            hTransform = cmsCreateTransform (camprofile, TYPE_RGB_FLT, prophoto, TYPE_RGB_FLT, INTENT_RELATIVE_COLORIMETRIC, cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE );
            lcmsMutex->unlock ();
        }
        TMatrix toxyz, torgb;
        if (!working_space_is_prophoto) {
            toxyz = iccStore->workingSpaceMatrix ("ProPhoto");
            torgb = iccStore->workingSpaceInverseMatrix (cmp.working); //sRGB .. Adobe...Wide...
        }

        #ifdef _OPENMP
        #pragma omp parallel
        #endif
        {
            AlignedBuffer<float> buffer(im->width*3);
            AlignedBuffer<float> hl_buffer(im->width*3);
            AlignedBuffer<float> hl_scale(im->width);
            #ifdef _OPENMP
            #pragma omp for schedule(static)
            #endif
            for ( int h = 0; h < im->height; ++h ) {
                float *p=buffer.data, *pR=im->r(h), *pG=im->g(h), *pB=im->b(h);

                // Apply pre-processing
                for ( int w = 0; w < im->width; ++w ) {
                    float r = *(pR++);
                    float g = *(pG++);
                    float b = *(pB++);

                    // convert to 0-1 range as LCMS expects that
                    r /= 65535.0f;
                    g /= 65535.0f;
                    b /= 65535.0f;

                    float maxc = max(r,g,b);
                    if (maxc <= 1.0) {
                        hl_scale.data[w] = 1.0;
                    } else {
                        // highlight recovery extend the range past the clip point, which means we can get values larger than 1.0 here.
                        // LUT ICC profiles only work in the 0-1 range so we scale down to fit and restore after conversion.
                        hl_scale.data[w] = 1.0 / maxc;
                        r *= hl_scale.data[w];
                        g *= hl_scale.data[w];
                        b *= hl_scale.data[w];
                    }

                    switch (camera_icc_type) {
                    case CAMERA_ICC_TYPE_PHASE_ONE:
                        // Here we apply a curve similar to Capture One's "Film Standard" + gamma, the reason is that the LUTs embedded in the
                        // ICCs are designed to work on such input, and if you provide it with a different curve you don't get as good result.
                        // We will revert this curve after we've made the color transform. However when we revert the curve, we'll notice that
                        // highlight rendering suffers due to that the LUT transform don't expand well, therefore we do a less compressed
                        // conversion too and mix them, this gives us the highest quality and most flexible result.
                        hl_buffer.data[3*w+0] = pow_F(r, 1.0/1.8);
                        hl_buffer.data[3*w+1] = pow_F(g, 1.0/1.8);
                        hl_buffer.data[3*w+2] = pow_F(b, 1.0/1.8);
                        r = phaseOneIccCurveInv->getVal(r);
                        g = phaseOneIccCurveInv->getVal(g);
                        b = phaseOneIccCurveInv->getVal(b);
                        break;
                    case CAMERA_ICC_TYPE_LEAF: {
                        // Leaf profiles expect that the camera native RGB has been converted to Prophoto RGB
                        float newr = leaf_prophoto_mat[0][0]*r + leaf_prophoto_mat[0][1]*g + leaf_prophoto_mat[0][2]*b;
                        float newg = leaf_prophoto_mat[1][0]*r + leaf_prophoto_mat[1][1]*g + leaf_prophoto_mat[1][2]*b;
                        float newb = leaf_prophoto_mat[2][0]*r + leaf_prophoto_mat[2][1]*g + leaf_prophoto_mat[2][2]*b;
                        hl_buffer.data[3*w+0] = pow_F(newr, 1.0/1.8);
                        hl_buffer.data[3*w+1] = pow_F(newg, 1.0/1.8);
                        hl_buffer.data[3*w+2] = pow_F(newb, 1.0/1.8);
                        r = phaseOneIccCurveInv->getVal(newr);
                        g = phaseOneIccCurveInv->getVal(newg);
                        b = phaseOneIccCurveInv->getVal(newb);
                        break;
                    }
                    case CAMERA_ICC_TYPE_NIKON:
                        // gamma 0.5
                        r = sqrtf(r);
                        g = sqrtf(g);
                        b = sqrtf(b);
                        break;
                    case CAMERA_ICC_TYPE_GENERIC:
                    default:
                        // do nothing
                        break;
                    }

                    *(p++) = r; *(p++) = g; *(p++) = b;
                }

                // Run icc transform
                cmsDoTransform (hTransform, buffer.data, buffer.data, im->width);
                if (separate_pcs_lab_highlights) {
                    cmsDoTransform (hTransform, hl_buffer.data, hl_buffer.data, im->width);
                }

                // Apply post-processing
                p=buffer.data; pR=im->r(h); pG=im->g(h); pB=im->b(h);
                for ( int w = 0; w < im->width; ++w ) {

                    float r, g, b, hr, hg, hb;

                    if (transform_via_pcs_lab) {
                        float L = *(p++);
                        float A = *(p++);
                        float B = *(p++);
                        // profile connection space CIELAB should have D50 illuminant
                        lab2ProphotoRgbD50(L, A, B, r, g, b);
                        if (separate_pcs_lab_highlights) {
                            lab2ProphotoRgbD50(hl_buffer.data[3*w+0], hl_buffer.data[3*w+1], hl_buffer.data[3*w+2], hr, hg, hb);
                        }
                    } else {
                        r = *(p++);
                        g = *(p++);
                        b = *(p++);
                    }

                    // restore pre-processing and/or add post-processing for the various ICC types
                    switch (camera_icc_type) {
                    default:
                        break;
                    case CAMERA_ICC_TYPE_PHASE_ONE:
                    case CAMERA_ICC_TYPE_LEAF: {
                        // note the 1/1.8 gamma, it's the gamma that the profile has applied, which we must revert before we can revert the curve
                        r = phaseOneIccCurve->getVal(pow_F(r, 1.0/1.8));
                        g = phaseOneIccCurve->getVal(pow_F(g, 1.0/1.8));
                        b = phaseOneIccCurve->getVal(pow_F(b, 1.0/1.8));
                        const float mix = 0.25; // may seem a low number, but remember this is linear space, mixing starts 2 stops from clipping
                        const float maxc = max(r, g, b);
                        if (maxc > mix) {
                            float fac = (maxc - mix) / (1.0 - mix);
                            fac = sqrtf(sqrtf(fac)); // gamma 0.25 to mix in highlight render relatively quick
                            r = (1.0-fac) * r + fac * hr;
                            g = (1.0-fac) * g + fac * hg;
                            b = (1.0-fac) * b + fac * hb;
                        }
                        break;
                    }
                    case CAMERA_ICC_TYPE_NIKON: {
                        const float lineFac = -0.4;
                        const float lineSum = 1.35;
                        r *= r * lineFac + lineSum;
                        g *= g * lineFac + lineSum;
                        b *= b * lineFac + lineSum;
                        break;
                    }
                    }

                    // restore highlight scaling if any
                    if (hl_scale.data[w] != 1.0) {
                        float fac = 1.0 / hl_scale.data[w];
                        r *= fac;
                        g *= fac;
                        b *= fac;
                    }

                    // If we don't have ProPhoto as chosen working profile, convert. This conversion is clipless, ie if we convert
                    // to a small space such as sRGB we may end up with negative values and values larger than max.
                    if (!working_space_is_prophoto) {
                        //convert from Prophoto to XYZ
                        float x = (toxyz[0][0] * r + toxyz[0][1] * g + toxyz[0][2] * b ) ;
                        float y = (toxyz[1][0] * r + toxyz[1][1] * g + toxyz[1][2] * b ) ;
                        float z = (toxyz[2][0] * r + toxyz[2][1] * g + toxyz[2][2] * b ) ;
                        //convert from XYZ to cmp.working  (sRGB...Adobe...Wide..)
                        r = ((torgb[0][0]*x + torgb[0][1]*y + torgb[0][2]*z)) ;
                        g = ((torgb[1][0]*x + torgb[1][1]*y + torgb[1][2]*z)) ;
                        b = ((torgb[2][0]*x + torgb[2][1]*y + torgb[2][2]*z)) ;
                    }

                    // return to the 0.0 - 65535.0 range (with possible negative and > max values present)
                    r *= 65535.0;
                    g *= 65535.0;
                    b *= 65535.0;

                    *(pR++) = r; *(pG++) = g; *(pB++) = b;
                }
            }
        } // End of parallelization
        cmsDeleteTransform(hTransform);
    }
    
//t3.set ();
//        printf ("ICM TIME: %d usec\n", t3.etime(t1));
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Converts raw image including ICC input profile to working space - 16bit int version
/*void RawImageSource::colorSpaceConversion16 (Image16* im, ColorManagementParams cmp, cmsHPROFILE embedded, cmsHPROFILE camprofile, double camMatrix[3][3], std::string camName) {
    cmsHPROFILE in;
    DCPProfile *dcpProf;

    if (!findInputProfile(cmp.input, embedded, camName, &dcpProf, in)) return;

    if (dcpProf!=NULL) {
        dcpProf->Apply(im, (DCPLightType)cmp.preferredProfile, cmp.working, cmp.toneCurve);
    } else {
        if (in==NULL) {
            // Take camprofile from DCRAW
            // in this case we avoid using the slllllooooooowwww lcms
            TMatrix work = iccStore->workingSpaceInverseMatrix (cmp.working);
            double mat[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
            for (int i=0; i<3; i++)
                for (int j=0; j<3; j++)
                    for (int k=0; k<3; k++)
                        mat[i][j] += work[i][k] * camMatrix[k][j]; // rgb_xyz * imatrices.xyz_cam

#pragma omp parallel for
            for (int i=0; i<im->height; i++)
                for (int j=0; j<im->width; j++) {

                    float newr = mat[0][0]*im->r(i,j) + mat[0][1]*im->g(i,j) + mat[0][2]*im->b(i,j);
                    float newg = mat[1][0]*im->r(i,j) + mat[1][1]*im->g(i,j) + mat[1][2]*im->b(i,j);
                    float newb = mat[2][0]*im->r(i,j) + mat[2][1]*im->g(i,j) + mat[2][2]*im->b(i,j);

                    im->r(i,j) = CLIP((int)newr);
                    im->g(i,j) = CLIP((int)newg);
                    im->b(i,j) = CLIP((int)newb);
                }
        } else {
            // Gamma preprocessing
            float gammaFac, lineFac, lineSum;
            getProfilePreprocParams(in, gammaFac, lineFac, lineSum);

            if (gammaFac>0) {
#pragma omp parallel for
                for ( int h = 0; h < im->height; ++h )
                    for ( int w = 0; w < im->width; ++w ) {
                        im->r(h,w)=  (int) (pow ((double)(im->r(h,w) / 65535.0), (double)gammaFac) * 65535.0);
                        im->g(h,w)=  (int) (pow ((double)(im->g(h,w) / 65535.0), (double)gammaFac) * 65535.0);
                        im->b(h,w)=  (int) (pow ((double)(im->b(h,w) / 65535.0), (double)gammaFac) * 65535.0);
                    }
            }

            cmsHPROFILE out = iccStore->workingSpace (cmp.working);
            //out = iccStore->workingSpaceGamma (wProfile);
            lcmsMutex->lock ();
            cmsHTRANSFORM hTransform = cmsCreateTransform (in, TYPE_RGB_16, out, TYPE_RGB_16, settings->colorimetricIntent,
            cmsFLAGS_NOCACHE);  // NOCACHE is important for thread safety
            lcmsMutex->unlock ();

            if (hTransform) {
                im->ExecCMSTransform(hTransform);

                // There might be Nikon postprocessings
                if (lineSum>0) {
#pragma omp parallel for
                    for ( int h = 0; h < im->height; ++h )
                        for ( int w = 0; w < im->width; ++w ) {
                            im->r(h,w) *= im->r(h,w) * lineFac / 65535.0 + lineSum;
                            im->g(h,w) *= im->g(h,w) * lineFac / 65535.0 + lineSum;
                            im->b(h,w) *= im->b(h,w) * lineFac / 65535.0 + lineSum;
                        }
                }
            }
            else {
                lcmsMutex->lock ();
                hTransform = cmsCreateTransform (camprofile, TYPE_RGB_16, out, TYPE_RGB_16,
                              settings->colorimetricIntent, cmsFLAGS_NOCACHE);
                lcmsMutex->unlock ();

                im->ExecCMSTransform(hTransform);
            }

            cmsDeleteTransform(hTransform);
        }
    }
    //t3.set ();
    //printf ("ICM TIME: %d\n", t3.etime(t1));
}*/

// Determine RAW input and output profiles. Returns TRUE on success
bool RawImageSource::findInputProfile(Glib::ustring inProfile, cmsHPROFILE embedded, std::string camName, DCPProfile **dcpProf, cmsHPROFILE& in) {
    in=NULL; // cam will be taken on NULL
    *dcpProf=NULL;

    if (inProfile == "(none)") return false;

    if (inProfile == "(embedded)" && embedded) {
		in = embedded;
	} else if (inProfile=="(cameraICC)") {
        // DCPs have higher quality, so use them first
        *dcpProf=dcpStore->getStdProfile(camName);
        if (*dcpProf==NULL)  in = iccStore->getStdProfile(camName);
    } else if (inProfile!="(camera)" && inProfile!="") {
        Glib::ustring normalName=inProfile;
        if (!inProfile.compare (0, 5, "file:")) normalName=inProfile.substr(5);

        if (dcpStore->isValidDCPFileName(normalName)) *dcpProf=dcpStore->getProfile(normalName);
        if (*dcpProf==NULL) in = iccStore->getProfile (inProfile);
    }
    
    // "in" might be NULL because of "not found". That's ok, we take the cam profile then
	
    return true;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// derived from Dcraw "blend_highlights()"
//  very effective to reduce (or remove) the magenta, but with levels of grey !
void RawImageSource::HLRecovery_blend(float* rin, float* gin, float* bin, int width, float maxval, float* pre_mul, const RAWParams &raw, float* hlmax)
	{
		const int ColorCount=3;

		// Transform matrixes rgb>lab and back
		static const float trans[2][ColorCount][ColorCount] =
		{ { { 1,1,1 }, { 1.7320508,-1.7320508,0 }, { -1,-1,2 } },
			{ { 1,1,1 }, { 1,-1,1 }, { 1,1,-1 } } };
		static const float itrans[2][ColorCount][ColorCount] =
		{ { { 1,0.8660254,-0.5 }, { 1,-0.8660254,-0.5 }, { 1,0,1 } },
			{ { 1,1,1 }, { 1,-1,1 }, { 1,1,-1 } } };
		
#define FOREACHCOLOR for (int c=0; c < ColorCount; c++)
		
		float minpt=min(hlmax[0],hlmax[1],hlmax[2]);//min of the raw clip points
		//float maxpt=max(hlmax[0],hlmax[1],hlmax[2]);//max of the raw clip points
		//float medpt=hlmax[0]+hlmax[1]+hlmax[2]-minpt-maxpt;//median of the raw clip points
		float maxave=(hlmax[0]+hlmax[1]+hlmax[2])/3;//ave of the raw clip points
		//some thresholds:
		const float clipthresh = 0.95;
		const float fixthresh = 0.5;
		const float satthresh = 0.5;

		float clip[3];
		FOREACHCOLOR clip[c]=min(maxave,hlmax[c]);
		
		// Determine the maximum level (clip) of all channels
		const float clippt = clipthresh*maxval;
		const float fixpt = fixthresh*minpt;
		const float desatpt = satthresh*maxave+(1-satthresh)*maxval;


#pragma omp parallel for
		for (int col=0; col<width; col++) {
			float rgb[ColorCount], cam[2][ColorCount], lab[2][ColorCount], sum[2], chratio, lratio=0;
			float L,C,H,Lfrac;
			
			// Copy input pixel to rgb so it's easier to access in loops
			rgb[0] = rin[col]; rgb[1] = gin[col]; rgb[2] = bin[col];
			
			// If no channel is clipped, do nothing on pixel
			int c;
			for (c=0; c<ColorCount; c++) { if (rgb[c] > clippt) break; }
			if (c == ColorCount) continue;
			
			// Initialize cam with raw input [0] and potentially clipped input [1]
			FOREACHCOLOR {
				lratio += min(rgb[c],clip[c]);
				cam[0][c] = rgb[c];
				cam[1][c] = min(cam[0][c],maxval);
			}
			
			// Calculate the lightness correction ratio (chratio)
			for (int i=0; i<2; i++) {
				FOREACHCOLOR {
					lab[i][c]=0;
					for (int j=0; j < ColorCount; j++)
						lab[i][c] += trans[ColorCount-3][c][j] * cam[i][j];
				}
				
				sum[i]=0;
				for (int c=1; c < ColorCount; c++)
					sum[i] += SQR(lab[i][c]);
			}
			chratio = (sqrt(sum[1]/sum[0]));
			
			// Apply ratio to lightness in LCH space
			for (int c=1; c < ColorCount; c++) 
				lab[0][c] *= chratio;
			
			// Transform back from LCH to RGB
			FOREACHCOLOR {
				cam[0][c]=0;
				for (int j=0; j < ColorCount; j++) {
					cam[0][c] += itrans[ColorCount-3][c][j] * lab[0][j];
				}
			}
			FOREACHCOLOR rgb[c] = cam[0][c] / ColorCount;

			// Copy converted pixel back			
			if (rin[col] > fixpt) {
				float rfrac = SQR((min(clip[0],rin[col])-fixpt)/(clip[0]-fixpt));
				rin[col]= min(maxave,rfrac*rgb[0]+(1-rfrac)*rin[col]);
			}
			if (gin[col] > fixpt) {
				float gfrac = SQR((min(clip[1],gin[col])-fixpt)/(clip[1]-fixpt));
				gin[col]= min(maxave,gfrac*rgb[1]+(1-gfrac)*gin[col]);
			}
			if (bin[col] > fixpt) {
				float bfrac = SQR((min(clip[2],bin[col])-fixpt)/(clip[2]-fixpt));
				bin[col]= min(maxave,bfrac*rgb[2]+(1-bfrac)*bin[col]);
			}
			
			lratio /= (rin[col]+gin[col]+bin[col]);
			L = (rin[col]+gin[col]+bin[col])/3;
			C = lratio * 1.732050808 * (rin[col] - gin[col]);
			H = lratio * (2 * bin[col] - rin[col] - gin[col]);
			rin[col] = L - H / 6.0 + C / 3.464101615;
			gin[col] = L - H / 6.0 - C / 3.464101615;
			bin[col] = L + H / 3.0;
			
			if ((L=(rin[col]+gin[col]+bin[col])/3) > desatpt) {
				Lfrac = max(0.0f,(maxave-L)/(maxave-desatpt));
				C = Lfrac * 1.732050808 * (rin[col] - gin[col]);
				H = Lfrac * (2 * bin[col] - rin[col] - gin[col]);
				rin[col] = L - H / 6.0 + C / 3.464101615;
				gin[col] = L - H / 6.0 - C / 3.464101615;
				bin[col] = L + H / 3.0;
			}
		}
	}

void RawImageSource::HLRecovery_Luminance (float* rin, float* gin, float* bin, float* rout, float* gout, float* bout, int width, float maxval) {

    for (int i=0; i<width; i++) {
        float r = rin[i], g = gin[i], b = bin[i];
		if (r>maxval || g>maxval || b>maxval) {
		    float ro = min(r, maxval);
		    float go = min(g, maxval);
		    float bo = min(b, maxval);
            double L = r + g + b;
            double C = 1.732050808 * (r - g);
            double H = 2 * b - r - g;
            double Co = 1.732050808 * (ro - go);
            double Ho = 2 * bo - ro - go;
            if (r!=g && g!=b) {
                double ratio = sqrt ((Co*Co+Ho*Ho) / (C*C+H*H));
                C *= ratio;
                H *= ratio;
            }
            float rr = L / 3.0 - H / 6.0 + C / 3.464101615;
            float gr = L / 3.0 - H / 6.0 - C / 3.464101615;
            float br = L / 3.0 + H / 3.0;
			rout[i] = rr;
			gout[i] = gr;
			bout[i] = br;
		}
        else {
            rout[i] = rin[i];
            gout[i] = gin[i];
            bout[i] = bin[i];
        }
    }
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::HLRecovery_CIELab (float* rin, float* gin, float* bin, float* rout, float* gout, float* bout,
										int width, float maxval, double xyz_cam[3][3], double cam_xyz[3][3]) {

    //static bool crTableReady = false;
	
	// lookup table for Lab conversion
	// perhaps should be centralized, universally defined so we don't keep remaking it???
    /*for (int ix=0; ix < 0x10000; ix++) {
    	    float rx = ix / 65535.0;
        	fv[ix] = rx > 0.008856 ? exp(1.0/3 * log(rx)) : 7.787*rx + 16/116.0;
    	}*/
    	//crTableReady = true;


    for (int i=0; i<width; i++) {
        float r = rin[i], g = gin[i], b = bin[i];
		if (r>maxval || g>maxval || b>maxval) {
		    float ro = min(r, maxval);
		    float go = min(g, maxval);
		    float bo = min(b, maxval);
            float yy = xyz_cam[1][0]*r + xyz_cam[1][1]*g + xyz_cam[1][2]*b;
            float fy = (yy<65535.0 ? ImProcFunctions::cachef[yy]/327.68 : (exp(log(yy/MAXVALD)/3.0 )));
            // compute LCH decompostion of the clipped pixel (only color information, thus C and H will be used)
            float x = xyz_cam[0][0]*ro + xyz_cam[0][1]*go + xyz_cam[0][2]*bo;
            float y = xyz_cam[1][0]*ro + xyz_cam[1][1]*go + xyz_cam[1][2]*bo;
            float z = xyz_cam[2][0]*ro + xyz_cam[2][1]*go + xyz_cam[2][2]*bo;
			x = (x<65535.0 ? ImProcFunctions::cachef[x]/327.68 : (exp(log(x/MAXVALD)/3.0 )));
			y = (y<65535.0 ? ImProcFunctions::cachef[y]/327.68 : (exp(log(y/MAXVALD)/3.0 )));
			z = (z<65535.0 ? ImProcFunctions::cachef[z]/327.68 : (exp(log(z/MAXVALD)/3.0 )));
            // convert back to rgb
            double fz = fy - y + z;
            double fx = fy + x - y;

			double zr = Color::f2xyz(fz);
            double xr = Color::f2xyz(fx);

            x = xr*65535.0 ;
            y = yy;
            z = zr*65535.0 ;
            float rr = cam_xyz[0][0]*x + cam_xyz[0][1]*y + cam_xyz[0][2]*z;
            float gr = cam_xyz[1][0]*x + cam_xyz[1][1]*y + cam_xyz[1][2]*z;
            float br = cam_xyz[2][0]*x + cam_xyz[2][1]*y + cam_xyz[2][2]*z;
			rout[i] = (rr);
			gout[i] = (gr);
			bout[i] = (br);
		}
        else {
            rout[i] = (rin[i]);
            gout[i] = (gin[i]);
            bout[i] = (bin[i]);
        }
    }
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::hlRecovery (std::string method, float* red, float* green, float* blue, int i, int sx1, int width, int skip,const RAWParams &raw, float* hlmax ) {
	
    if (method=="Luminance")
        HLRecovery_Luminance (red, green, blue, red, green, blue, width, 65535.0);
    else if (method=="CIELab blending")
        HLRecovery_CIELab (red, green, blue, red, green, blue, width, 65535.0, imatrices.xyz_cam, imatrices.cam_xyz);
    /*else if (method=="Color")
        HLRecovery_ColorPropagation (red, green, blue, i, sx1, width, skip);*/
	else if (method=="Blend")	// derived from Dcraw
			{	float pre_mul[4];
				for(int c=0;c<4;c++) pre_mul[c]=ri->get_pre_mul(c);
				HLRecovery_blend(red, green, blue, width, 65535.0, pre_mul, raw, hlmax );}

}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::getAutoExpHistogram (LUTu & histogram, int& histcompr) {

    histcompr = 3;

    histogram(65536>>histcompr);
    histogram.clear();

#pragma omp parallel
{
	LUTu tmphistogram(65536>>histcompr);
	tmphistogram.clear();
#pragma omp for nowait
    for (int i=border; i<H-border; i++) {
        int start, end;
        getRowStartEnd (i, start, end);

        if (ri->getSensorType()==ST_BAYER) {
            for (int j=start; j<end; j++) {
				if (ri->ISGREEN(i,j))     tmphistogram[CLIP((int)(refwb_green*rawData[i][j]))>>histcompr]+=4;
				else if (ri->ISRED(i,j))  tmphistogram[CLIP((int)(refwb_red*  rawData[i][j]))>>histcompr]+=4;
				else if (ri->ISBLUE(i,j)) tmphistogram[CLIP((int)(refwb_blue* rawData[i][j]))>>histcompr]+=4;
			} 
        } else if (ri->getSensorType()==ST_FUJI_XTRANS) {
            for (int j=start; j<end; j++) {
				if (ri->ISXTRANSGREEN(i,j))     tmphistogram[CLIP((int)(refwb_green*rawData[i][j]))>>histcompr]+=4;
				else if (ri->ISXTRANSRED(i,j))  tmphistogram[CLIP((int)(refwb_red*  rawData[i][j]))>>histcompr]+=4;
				else if (ri->ISXTRANSBLUE(i,j)) tmphistogram[CLIP((int)(refwb_blue* rawData[i][j]))>>histcompr]+=4;
			} 
		} else if (ri->get_colors() == 1) {
			for (int j=start; j<end; j++) {
				tmphistogram[CLIP((int)(refwb_red*  rawData[i][j]))>>histcompr]++;
			}
		} else {
			for (int j=start; j<end; j++) {
				tmphistogram[CLIP((int)(refwb_red*  rawData[i][3*j+0]))>>histcompr]++;
				tmphistogram[CLIP((int)(refwb_green*rawData[i][3*j+1]))>>histcompr]+=2;
				tmphistogram[CLIP((int)(refwb_blue* rawData[i][3*j+2]))>>histcompr]++;
			}
		}
    }
#pragma omp critical
{
	for(int i=0; i<(65536>>histcompr);i++)
		histogram[i] += tmphistogram[i];
}
}
}
		
// Histogram MUST be 256 in size; gamma is applied, blackpoint and gain also
void RawImageSource::getRAWHistogram (LUTu & histRedRaw, LUTu & histGreenRaw, LUTu & histBlueRaw) {

	histRedRaw.clear(); histGreenRaw.clear(); histBlueRaw.clear();
	const float mult[4] = { 65535.0 / ri->get_white(0), 65535.0 / ri->get_white(1), 65535.0 / ri->get_white(2), 65535.0 / ri->get_white(3) };
	
#ifdef _OPENMP
	int numThreads;
	// reduce the number of threads under certain conditions to avoid overhaed of too many critical regions
	numThreads = sqrt((((H-2*border)*(W-2*border))/262144.f));
	numThreads = std::min(std::max(numThreads,1), omp_get_max_threads());

#pragma omp parallel num_threads(numThreads)
#endif
{
	// we need one LUT per color and thread, which corresponds to 1 MB per thread
	LUTu tmphist[4];
	tmphist[0](65536);tmphist[0].clear();
	tmphist[1](65536);tmphist[1].clear();
	tmphist[2](65536);tmphist[2].clear();
	tmphist[3](65536);tmphist[3].clear();
	
#ifdef _OPENMP	
#pragma omp for nowait
#endif
	for (int i=border; i<H-border; i++) {
		int start, end;
		getRowStartEnd (i, start, end);

		if (ri->getSensorType()==ST_BAYER) {
			int j;
			int c1 = FC(i,start);
			c1 = ( c1 == 1 && !(i&1) ) ? 3 : c1;
			int c2 = FC(i,start+1);
			c2 = ( c2 == 1 && !(i&1) ) ? 3 : c2;
			for (j=start; j<end-1; j+=2) {
				tmphist[c1][(int)ri->data[i][j]]++;
				tmphist[c2][(int)ri->data[i][j+1]]++;
			}
			if(j<end) { // last pixel of row if width is odd
				tmphist[c1][(int)ri->data[i][j]]++;
			}
		} else if (ri->get_colors() == 1) {
			for (int j=start; j<end; j++) {
				for (int c=0; c<3; c++){
					tmphist[c][(int)ri->data[i][j]]++;
				}
			}
		} else if(ri->getSensorType()==ST_FUJI_XTRANS) {
			for (int j=start; j<end-1; j+=2) {
				int c = ri->XTRANSFC(i,j);
				tmphist[c][(int)ri->data[i][j]]++;
			}
		} else {
			for (int j=start; j<end; j++) {
				for (int c=0; c<3; c++){
					tmphist[c][(int)ri->data[i][3*j+c]]++;
				}
			}
		}
	}
#ifdef _OPENMP
#pragma omp critical
#endif
{
	for(int i=0;i<65536;i++){
		int idx;
		idx = CLIP((int)Color::gamma(mult[0]*(i-(cblacksom[0]/*+black_lev[0]*/))));
		histRedRaw[idx>>8] += tmphist[0][i];
		idx = CLIP((int)Color::gamma(mult[1]*(i-(cblacksom[1]/*+black_lev[1]*/))));
		histGreenRaw[idx>>8] += tmphist[1][i];
		idx = CLIP((int)Color::gamma(mult[3]*(i-(cblacksom[3]/*+black_lev[3]*/))));
		histGreenRaw[idx>>8] += tmphist[3][i];
		idx = CLIP((int)Color::gamma(mult[2]*(i-(cblacksom[2]/*+black_lev[2]*/))));
		histBlueRaw[idx>>8] += tmphist[2][i];
	}
} // end of critical region
} // end of parallel region

    
    if (ri->getSensorType()==ST_BAYER)		// since there are twice as many greens, correct for it
		for (int i=0;i<256;i++)
			histGreenRaw[i]>>=1;
	else if(ri->getSensorType()==ST_FUJI_XTRANS)	// since Xtrans has 2.5 as many greens, correct for it
		for (int i=0;i<256;i++)
			histGreenRaw[i] = (histGreenRaw[i]*2)/5;
			

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

void RawImageSource::getRowStartEnd (int x, int &start, int &end) {
    if (fuji) {
        int fw = ri->get_FujiWidth();
        start = ABS(fw-x) + border;
        end = min(H+ W-fw-x, fw+x) - border;
    }
    else {
        start = border;
        end = W-border;
    }
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	void RawImageSource::getAutoWBMultipliers (double &rm, double &gm, double &bm) {

		if (ri->get_colors() == 1) {
			rm = gm = bm = 1;
			return;
		}
		if (redAWBMul != -1.) {
			rm = redAWBMul;
			gm = greenAWBMul;
			bm = blueAWBMul;
			return;
		}

		if (!isWBProviderReady()) {
			rm = -1.0;
			gm = -1.0;
			bm = -1.0;
			return;
		}

		double avg_r = 0;
		double avg_g = 0;
		double avg_b = 0;
		int rn = 0, gn = 0, bn = 0;
		
		if (fuji) {
			for (int i=32; i<H-32; i++) {
				int fw = ri->get_FujiWidth();
				int start = ABS(fw-i) + 32;
				int end = min(H+W-fw-i, fw+i) - 32;
				for (int j=start; j<end; j++) {
					if (ri->getSensorType()!=ST_BAYER) {
						double dr = CLIP(initialGain*(rawData[i][3*j]  ));
						double dg = CLIP(initialGain*(rawData[i][3*j+1]));
						double db = CLIP(initialGain*(rawData[i][3*j+2]));
						if (dr>64000. || dg>64000. || db>64000.) continue;
						avg_r += dr;
						avg_g += dg;
						avg_b += db;
						rn = gn = ++bn;
					}
					else {
						int c = FC( i, j);
						double d = CLIP(initialGain*(rawData[i][j]));
						if (d>64000.)
							continue;
						// Let's test green first, because they are more numerous
						if (c==1) {
							avg_g += d;
							gn++;
						}
						else if (c==0) {
							avg_r += d;
							rn++;
						}
						else /*if (c==2)*/ {
							avg_b += d;
							bn++;
						}
					}
				}
			}
		}
		else {
			if (ri->getSensorType()!=ST_BAYER) {
				if(ri->getSensorType()==ST_FUJI_XTRANS) {
					for (int i=32; i<H-32; i++)
						for (int j=32; j<W-32; j++) {
							// each loop read 1 rgb triplet value
							if(ri->ISXTRANSRED(i,j)) {
								float dr = CLIP(initialGain*(rawData[i][j]));
								if (dr>64000.f)
									continue;
								avg_r += dr;
								rn ++;
							}
							if(ri->ISXTRANSGREEN(i,j)) {
								float dg = CLIP(initialGain*(rawData[i][j]));
								if (dg>64000.f)
									continue;
								avg_g += dg;
								gn ++;
							}
							if(ri->ISXTRANSBLUE(i,j)) {
								float db = CLIP(initialGain*(rawData[i][j]));
								if (db>64000.f)
									continue;
								avg_b += db;
								bn ++;
							}
						}
				} else {
					for (int i=32; i<H-32; i++)
						for (int j=32; j<W-32; j++) {
							// each loop read 1 rgb triplet value

							double dr = CLIP(initialGain*(rawData[i][3*j]  ));
							double dg = CLIP(initialGain*(rawData[i][3*j+1]));
							double db = CLIP(initialGain*(rawData[i][3*j+2]));
							if (dr>64000. || dg>64000. || db>64000.) continue;
							avg_r += dr; rn++;
							avg_g += dg; 
							avg_b += db; 
						}
					gn = rn; bn=rn;
				}
			} else {
				//determine GRBG coset; (ey,ex) is the offset of the R subarray
				int ey, ex;
				if (ri->ISGREEN(0,0)) {//first pixel is G
					if (ri->ISRED(0,1)) {ey=0; ex=1;} else {ey=1; ex=0;}
				} else {//first pixel is R or B
					if (ri->ISRED(0,0)) {ey=0; ex=0;} else {ey=1; ex=1;}
				}
				double d[2][2];
				for (int i=32; i<H-32; i+=2)
					for (int j=32; j<W-32; j+=2) {
						//average each Bayer quartet component individually if non-clipped
						d[0][0] = CLIP(initialGain*(rawData[i][j]    ));
						d[0][1] = CLIP(initialGain*(rawData[i][j+1]  ));
						d[1][0] = CLIP(initialGain*(rawData[i+1][j]  ));
						d[1][1] = CLIP(initialGain*(rawData[i+1][j+1]));
						if (d[ey][ex] <= 64000.) {
							avg_r += d[ey][ex];
							rn++;
						}
						if (d[1-ey][ex] <= 64000.) {
							avg_g += d[1-ey][ex];
							gn++;
						}
						if (d[ey][1-ex] <= 64000.) {
							avg_g += d[ey][1-ex];
							gn++;
						}
						if (d[1-ey][1-ex] <= 64000.) {
							avg_b += d[1-ey][1-ex];
							bn++;
						}
					}
			}
		}
		if( settings->verbose )
			printf ("AVG: %g %g %g\n", avg_r/rn, avg_g/gn, avg_b/bn);
		
		//    return ColorTemp (pow(avg_r/rn, 1.0/6.0)*img_r, pow(avg_g/gn, 1.0/6.0)*img_g, pow(avg_b/bn, 1.0/6.0)*img_b);
		
		double reds   = avg_r/rn * refwb_red;
		double greens = avg_g/gn * refwb_green;
		double blues  = avg_b/bn * refwb_blue;
		
		redAWBMul   = rm = imatrices.rgb_cam[0][0]*reds + imatrices.rgb_cam[0][1]*greens + imatrices.rgb_cam[0][2]*blues;
		greenAWBMul = gm = imatrices.rgb_cam[1][0]*reds + imatrices.rgb_cam[1][1]*greens + imatrices.rgb_cam[1][2]*blues;
		blueAWBMul  = bm = imatrices.rgb_cam[2][0]*reds + imatrices.rgb_cam[2][1]*greens + imatrices.rgb_cam[2][2]*blues;
	}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	ColorTemp RawImageSource::getSpotWB (std::vector<Coord2D> &red, std::vector<Coord2D> &green, std::vector<Coord2D> &blue, int tran, double equal) {
		
		int x; int y;
		double reds = 0, greens = 0, blues = 0;
		int rn = 0;
		
		if (ri->getSensorType()!=ST_BAYER) {
			if(ri->getSensorType()==ST_FUJI_XTRANS) {
				int d[9][2] = {{0,0}, {-1,-1}, {-1,0}, {-1,1}, {0,-1}, {0,1}, {1,-1}, {1,0}, {1,1}};
				double rloc, gloc, bloc;
				int rnbrs, gnbrs, bnbrs;
				for (size_t i=0; i<red.size(); i++) {
					transformPosition (red[i].x, red[i].y, tran, x, y);
					rloc=gloc=bloc=rnbrs=gnbrs=bnbrs=0;
					for (int k=0; k<9; k++) {
						int xv = x + d[k][0];
						int yv = y + d[k][1];
						if(xv>=0 && yv>=0 && xv<W && yv<H) {
							if (ri->ISXTRANSRED(yv,xv)) { //RED
								rloc += (rawData[yv][xv]);
								rnbrs++;
								continue;
							} else if (ri->ISXTRANSBLUE(yv,xv)) { //BLUE
								bloc += (rawData[yv][xv]);
								bnbrs++;
								continue;
							} else { // GREEN
								gloc += (rawData[yv][xv]);
								gnbrs++;
								continue;
							}
						}
					}
					rloc /= rnbrs; gloc /= gnbrs; bloc /= bnbrs;
					if (rloc*initialGain<64000. && gloc*initialGain<64000. && bloc*initialGain<64000.) {
						reds += rloc; greens += gloc; blues += bloc; rn++;
					}
				}
				
			} else {
				int xmin, xmax, ymin, ymax;
				int xr, xg, xb, yr, yg, yb;
				for (size_t i=0; i<red.size(); i++) {
					transformPosition (red[i].x, red[i].y, tran, xr, yr);
					transformPosition (green[i].x, green[i].y, tran, xg, yg);
					transformPosition (blue[i].x, blue[i].y, tran, xb, yb);
					if (initialGain*(rawData[yr][3*xr]  )>52500 ||      
						initialGain*(rawData[yg][3*xg+1])>52500 ||        
						initialGain*(rawData[yb][3*xb+2])>52500) continue;
					xmin = min(xr,xg,xb);
					xmax = max(xr,xg,xb);
					ymin = min(yr,yg,yb);
					ymax = max(yr,yg,yb);
					if (xmin>=0 && ymin>=0 && xmax<W && ymax<H) {
						reds	+= (rawData[yr][3*xr]  );  
						greens	+= (rawData[yg][3*xg+1]);
						blues	+= (rawData[yb][3*xb+2]);  
						rn++;
					}
				}
			}
			
		} else {
			
			int d[9][2] = {{0,0}, {-1,-1}, {-1,0}, {-1,1}, {0,-1}, {0,1}, {1,-1}, {1,0}, {1,1}};
			double rloc, gloc, bloc;
			int rnbrs, gnbrs, bnbrs;
			for (size_t i=0; i<red.size(); i++) {
				transformPosition (red[i].x, red[i].y, tran, x, y);
				rloc=gloc=bloc=rnbrs=gnbrs=bnbrs=0;
				for (int k=0; k<9; k++) {
					int xv = x + d[k][0];
					int yv = y + d[k][1];
					int c = FC(yv,xv);
					if(xv>=0 && yv>=0 && xv<W && yv<H) {
						if (c==0) { //RED
							rloc += (rawData[yv][xv]);
							rnbrs++;
							continue;
						}else if (c==2) { //BLUE
							bloc += (rawData[yv][xv]);
							bnbrs++;
							continue;
						} else { // GREEN
							gloc += (rawData[yv][xv]);
							gnbrs++;
							continue;
						}
					}
				}
				rloc /= rnbrs; gloc /= gnbrs; bloc /= bnbrs;
				if (rloc*initialGain<64000. && gloc*initialGain<64000. && bloc*initialGain<64000.) {
					reds += rloc; greens += gloc; blues += bloc; rn++;
				}
				transformPosition (green[i].x, green[i].y, tran, x, y);//these are redundant now ??? if not, repeat for these blocks same as for red[]
				rloc=gloc=bloc=rnbrs=gnbrs=bnbrs=0;
				for (int k=0; k<9; k++) {
					int xv = x + d[k][0];
					int yv = y + d[k][1];
					int c = FC(yv,xv);
					if(xv>=0 && yv>=0 && xv<W && yv<H) {
						if (c==0) { //RED
							rloc += (rawData[yv][xv]);
							rnbrs++;
							continue;
						}else if (c==2) { //BLUE
							bloc += (rawData[yv][xv]);
							bnbrs++;
							continue;
						} else { // GREEN
							gloc += (rawData[yv][xv]);
							gnbrs++;
							continue;
						}
					}
				}
				rloc /= rnbrs; gloc /= gnbrs; bloc /= bnbrs;
				if (rloc*initialGain<64000. && gloc*initialGain<64000. && bloc*initialGain<64000.) {
					reds += rloc; greens += gloc; blues += bloc; rn++;
				}

				transformPosition (blue[i].x, blue[i].y, tran, x, y);
				rloc=gloc=bloc=rnbrs=gnbrs=bnbrs=0;
				for (int k=0; k<9; k++) {
					int xv = x + d[k][0];
					int yv = y + d[k][1];
					int c = FC(yv,xv);
					if(xv>=0 && yv>=0 && xv<W && yv<H) {
						if (c==0) { //RED
							rloc += (rawData[yv][xv]);
							rnbrs++;
							continue;
						} else if (c==2) { //BLUE
							bloc += (rawData[yv][xv]);
							bnbrs++;
							continue;
						} else { // GREEN
							gloc += (rawData[yv][xv]);
							gnbrs++;
							continue;
						}
					}
				}
				rloc /= rnbrs; gloc /= gnbrs; bloc /= bnbrs;
				if (rloc*initialGain<64000. && gloc*initialGain<64000. && bloc*initialGain<64000.) {
					reds += rloc; greens += gloc; blues += bloc; rn++;
				}
			}
		}
		
		if (2*rn < red.size()) {
			return ColorTemp (equal);
		}
		else {
			reds = reds/rn * refwb_red;
			greens = greens/rn * refwb_green;
			blues = blues/rn * refwb_blue;
			
			double rm = imatrices.rgb_cam[0][0]*reds + imatrices.rgb_cam[0][1]*greens + imatrices.rgb_cam[0][2]*blues;
			double gm = imatrices.rgb_cam[1][0]*reds + imatrices.rgb_cam[1][1]*greens + imatrices.rgb_cam[1][2]*blues;
			double bm = imatrices.rgb_cam[2][0]*reds + imatrices.rgb_cam[2][1]*greens + imatrices.rgb_cam[2][2]*blues;
			
			return ColorTemp (rm, gm, bm, equal);
		}
	}
	

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::transformPosition (int x, int y, int tran, int& ttx, int& tty) {

    tran = defTransform (tran);

    x += border;
    y += border;

    if (d1x) {
        if ((tran & TR_ROT) == TR_R90 || (tran & TR_ROT) == TR_R270) 
            x /= 2;
        else
            y /= 2;
    }

    int w = W, h = H;  
    if (fuji) {
        w = ri->get_FujiWidth() * 2 + 1;
        h = (H - ri->get_FujiWidth())*2 + 1;
    }
    int sw = w, sh = h;  
    if ((tran & TR_ROT) == TR_R90 || (tran & TR_ROT) == TR_R270) {
        sw = h;
        sh = w;
    }
    
    int ppx = x, ppy = y;
    if (tran & TR_HFLIP) 
        ppx = sw - 1 - x ;
    if (tran & TR_VFLIP) 
        ppy = sh - 1 - y;
    
    int tx = ppx;
    int ty = ppy;
    
    if ((tran & TR_ROT) == TR_R180) {
        tx = w - 1 - ppx;
        ty = h - 1 - ppy;
    }
    else if ((tran & TR_ROT) == TR_R90) {
        tx = ppy;
        ty = h - 1 - ppx;
    }
    else if ((tran & TR_ROT) == TR_R270) {
        tx = w - 1 - ppy;
        ty = ppx;
    }   

    if (fuji) {
        ttx = (tx+ty) / 2;
        tty = (ty-tx) / 2 + ri->get_FujiWidth();
    }
    else {
        ttx = tx;
        tty = ty;
    }
}
		
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::inverse33 (const double (*rgb_cam)[3], double (*cam_rgb)[3]) {
	double nom = (rgb_cam[0][2]*rgb_cam[1][1]*rgb_cam[2][0] - rgb_cam[0][1]*rgb_cam[1][2]*rgb_cam[2][0] -
				  rgb_cam[0][2]*rgb_cam[1][0]*rgb_cam[2][1] + rgb_cam[0][0]*rgb_cam[1][2]*rgb_cam[2][1] +
				  rgb_cam[0][1]*rgb_cam[1][0]*rgb_cam[2][2] - rgb_cam[0][0]*rgb_cam[1][1]*rgb_cam[2][2] );
	cam_rgb[0][0] = (rgb_cam[1][2]*rgb_cam[2][1]-rgb_cam[1][1]*rgb_cam[2][2]) / nom;
	cam_rgb[0][1] = -(rgb_cam[0][2]*rgb_cam[2][1]-rgb_cam[0][1]*rgb_cam[2][2]) / nom;
	cam_rgb[0][2] = (rgb_cam[0][2]*rgb_cam[1][1]-rgb_cam[0][1]*rgb_cam[1][2]) / nom;
	cam_rgb[1][0] = -(rgb_cam[1][2]*rgb_cam[2][0]-rgb_cam[1][0]*rgb_cam[2][2]) / nom;
	cam_rgb[1][1] = (rgb_cam[0][2]*rgb_cam[2][0]-rgb_cam[0][0]*rgb_cam[2][2]) / nom;
	cam_rgb[1][2] = -(rgb_cam[0][2]*rgb_cam[1][0]-rgb_cam[0][0]*rgb_cam[1][2]) / nom;
	cam_rgb[2][0] = (rgb_cam[1][1]*rgb_cam[2][0]-rgb_cam[1][0]*rgb_cam[2][1]) / nom;
	cam_rgb[2][1] = -(rgb_cam[0][1]*rgb_cam[2][0]-rgb_cam[0][0]*rgb_cam[2][1]) / nom;
	cam_rgb[2][2] = (rgb_cam[0][1]*rgb_cam[1][0]-rgb_cam[0][0]*rgb_cam[1][1]) / nom;
}

DiagonalCurve* RawImageSource::phaseOneIccCurve;
DiagonalCurve* RawImageSource::phaseOneIccCurveInv;

void RawImageSource::init () {

    { // Initialize Phase One ICC curves

        /* This curve is derived from TIFFTAG_TRANSFERFUNCTION of a Capture One P25+ image with applied film curve,
           exported to TIFF with embedded camera ICC. It's assumed to be similar to most standard curves in
           Capture One. It's not necessary to be exactly the same, it's just to be close to a typical curve to
           give the Phase One ICC files a good working space. */
        const double phase_one_forward[] = {
            0.0000000000, 0.0000000000, 0.0152590219, 0.0029602502, 0.0305180438, 0.0058899825, 0.0457770657, 0.0087739376, 0.0610360876, 0.0115968566,
            0.0762951095, 0.0143587396, 0.0915541314, 0.0171969177, 0.1068131533, 0.0201876860, 0.1220721752, 0.0232852674, 0.1373311971, 0.0264744030,
            0.1525902190, 0.0297245747, 0.1678492409, 0.0330205234, 0.1831082628, 0.0363775082, 0.1983672847, 0.0397802701, 0.2136263066, 0.0432593271,
            0.2288853285, 0.0467841611, 0.2441443503, 0.0503700313, 0.2594033722, 0.0540474556, 0.2746623941, 0.0577859159, 0.2899214160, 0.0616159304,
            0.3051804379, 0.0655222400, 0.3204394598, 0.0695353628, 0.3356984817, 0.0736552987, 0.3509575036, 0.0778973068, 0.3662165255, 0.0822461280,
            0.3814755474, 0.0867170214, 0.3967345693, 0.0913252461, 0.4119935912, 0.0960860609, 0.4272526131, 0.1009994659, 0.4425116350, 0.1060654612,
            0.4577706569, 0.1113298238, 0.4730296788, 0.1167925536, 0.4882887007, 0.1224841688, 0.5035477226, 0.1284046693, 0.5188067445, 0.1345540551,
            0.5340657664, 0.1409781033, 0.5493247883, 0.1476615549, 0.5645838102, 0.1546501869, 0.5798428321, 0.1619287404, 0.5951018540, 0.1695277333,
            0.6103608759, 0.1774776837, 0.6256198978, 0.1858091096, 0.6408789197, 0.1945525292, 0.6561379416, 0.2037384604, 0.6713969635, 0.2134279393,
            0.6866559854, 0.2236667430, 0.7019150072, 0.2345159075, 0.7171740291, 0.2460517281, 0.7324330510, 0.2583047227, 0.7476920729, 0.2714122225,
            0.7629510948, 0.2854352636, 0.7782101167, 0.3004959182, 0.7934691386, 0.3167620356, 0.8087281605, 0.3343862058, 0.8239871824, 0.3535820554,
            0.8392462043, 0.3745937285, 0.8545052262, 0.3977111467, 0.8697642481, 0.4232547494, 0.8850232700, 0.4515754940, 0.9002822919, 0.4830701152,
            0.9155413138, 0.5190966659, 0.9308003357, 0.5615320058, 0.9460593576, 0.6136263066, 0.9613183795, 0.6807965209, 0.9765774014, 0.7717402914,
            0.9918364233, 0.9052109560, 1.0000000000, 1.0000000000
        };
        std::vector<double> cForwardPoints;
        cForwardPoints.push_back(double(DCT_Spline));  // The first value is the curve type
        std::vector<double> cInversePoints;
        cInversePoints.push_back(double(DCT_Spline));  // The first value is the curve type
        for (int i = 0; i < sizeof(phase_one_forward)/sizeof(phase_one_forward[0]); i += 2) {
            cForwardPoints.push_back(phase_one_forward[i+0]);
            cForwardPoints.push_back(phase_one_forward[i+1]);
            cInversePoints.push_back(phase_one_forward[i+1]);
            cInversePoints.push_back(phase_one_forward[i+0]);
        }
        phaseOneIccCurve = new DiagonalCurve(cForwardPoints, CURVES_MIN_POLY_POINTS);
        phaseOneIccCurveInv = new DiagonalCurve(cInversePoints, CURVES_MIN_POLY_POINTS);
    }
}

void RawImageSource::cleanup () {
    delete phaseOneIccCurve;
    delete phaseOneIccCurveInv;
}
	
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//#include "demosaic_algos.cc"
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Emil's code 
/*
 * Now compiled separately
 *
#include "fast_demo.cc"//fast demosaic	
#include "amaze_demosaic_RT.cc"//AMaZE demosaic	
#include "CA_correct_RT.cc"//Emil's CA auto correction
#include "cfa_linedn_RT.cc"//Emil's line denoise
#include "green_equil_RT.cc"//Emil's green channel equilibration
#include "hilite_recon.cc"//Emil's highlight reconstruction

#include "expo_before_b.cc"//Jacques's exposure before interpolation
*/
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PIX_SORT
#undef med3x3

} /* namespace */

#undef PIX_SORT
#undef med3x3

