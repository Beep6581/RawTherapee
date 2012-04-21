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
#include <algorithm>
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

namespace rtengine {

extern const Settings* settings;
#undef ABS
#undef DIST
#undef CLIP

#define ABS(a) ((a)<0?-(a):(a))
#define DIST(a,b) (ABS(a-b))
#define MAXVAL 0xffff
#define CLIP(a) ((a)>0?((a)<MAXVAL?(a):MAXVAL):0)
	
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
,rawData(NULL)
,green(NULL)
,red(NULL)
,blue(NULL)
{
    hrmap[0] = NULL;
    hrmap[1] = NULL;
    hrmap[2] = NULL;
	//needhr = NULL;
    //hpmap = NULL;
	camProfile = NULL;
	embProfile = NULL;
	rgbSourceModified = false;
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RawImageSource::~RawImageSource () {

    delete idata;
    if (ri) {
        delete ri;
    }

    if (green)
        freeArray<float>(green, H);
    if (red)
        freeArray<float>(red, H);
    if (blue)
        freeArray<float>(blue, H);
    if(rawData)
    	freeArray<float>(rawData, H);
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
       
    int sx1 = ppx;
    int sy1 = ppy;
    int sx2 = ppx + pp.w;
    int sy2 = ppy + pp.h;
    
    if ((tran & TR_ROT) == TR_R180) {
        sx1 = w - ppx - pp.w;
        sy1 = h - ppy - pp.h;
        sx2 = sx1 + pp.w;
        sy2 = sy1 + pp.h;
    }
    else if ((tran & TR_ROT) == TR_R90) {
        sx1 = ppy;
        sy1 = h - ppx - pp.w;
        sx2 = sx1 + pp.h;
        sy2 = sy1 + pp.w;
    }
    else if ((tran & TR_ROT) == TR_R270) {
        sx1 = w - ppy - pp.h;
        sy1 = ppx;
        sx2 = sx1 + pp.h;
        sy2 = sy1 + pp.w;
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

void RawImageSource::getImage (ColorTemp ctemp, int tran, Imagefloat* image, PreviewProps pp, HRecParams hrp, ColorManagementParams cmp, RAWParams raw )
{
    Glib::Mutex::Lock lock(getImageMutex);

    tran = defTransform (tran);

    // compute channel multipliers
    double r, g, b;
    float rm, gm, bm;
    ctemp.getMultipliers (r, g, b);
    rm = cam_rgb[0][0]*r + cam_rgb[0][1]*g + cam_rgb[0][2]*b;
    gm = cam_rgb[1][0]*r + cam_rgb[1][1]*g + cam_rgb[1][2]*b;
    bm = cam_rgb[2][0]*r + cam_rgb[2][1]*g + cam_rgb[2][2]*b;
    rm = camwb_red / rm;
    gm = camwb_green / gm;
    bm = camwb_blue / bm;

    /*float mul_lum = 0.299*rm + 0.587*gm + 0.114*bm;
    rm /= mul_lum;
    gm /= mul_lum;
    bm /= mul_lum;*/    

    /*//initialGain=1.0;
	// in floating point, should keep white point fixed and recover higher values with exposure slider
    //if (hrp.enabled) */
    float min = rm;
    if (min>gm) min = gm;
    if (min>bm) min = bm;
        defGain=0.0;// = log(initialGain) / log(2.0);
        //printf(" Initial gain=%f defgain=%f min=%f\n",initialGain,defGain,min);
        //printf(" rm=%f gm=%f bm=%f\n",rm,gm,bm);
        min/=initialGain;
   //min=(float)1.0/min;
    //else {
        //defGain = 0.0;
        rm /= min;
        gm /= min;
        bm /= min;
    //}
	//defGain = 0.0;//no need now for making headroom for highlights???
    //printf("initial gain= %e\n",initialGain);
    //TODO: normalize the gain control





	//if (hrp.enabled==true && hrp.method=="Color" && hrmap[0]==NULL) 
    //    updateHLRecoveryMap_ColorPropagation ();

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
#pragma omp parallel
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
		if (ri->isBayer()) {
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
				if (!hrp.enabled)
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
				if (!hrp.enabled)
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
        if (hrp.enabled && hrp.method!="Color")
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
          image->r[0][j] = (image->r[1][j] + image->r[0][j+1] + image->r[0][j-1]) / 3;
          image->g[0][j] = (image->g[1][j] + image->g[0][j+1] + image->g[0][j-1]) / 3;
          image->b[0][j] = (image->b[1][j] + image->b[0][j+1] + image->b[0][j-1]) / 3;
        }
        // other rows
        for (int i=1; i<image->height-1; i++) {
          for (int j=2-(a+i+1)%2; j<image->width-1; j+=2) {
              // edge-adaptive interpolation
              double dh = (ABS(image->r[i][j+1] - image->r[i][j-1]) + ABS(image->g[i][j+1] - image->g[i][j-1]) + ABS(image->b[i][j+1] - image->b[i][j-1])) / 1.0;
              double dv = (ABS(image->r[i+1][j] - image->r[i-1][j]) + ABS(image->g[i+1][j] - image->g[i-1][j]) + ABS(image->b[i+1][j] - image->b[i-1][j])) / 1.0;
              double eh = 1.0 / (1.0 + dh);
              double ev = 1.0 / (1.0 + dv);
              image->r[i][j] = (eh * (image->r[i][j+1] + image->r[i][j-1]) + ev * (image->r[i+1][j] + image->r[i-1][j])) / (2.0 * (eh + ev));
              image->g[i][j] = (eh * (image->g[i][j+1] + image->g[i][j-1]) + ev * (image->g[i+1][j] + image->g[i-1][j])) / (2.0 * (eh + ev));
              image->b[i][j] = (eh * (image->b[i][j+1] + image->b[i][j-1]) + ev * (image->b[i+1][j] + image->b[i-1][j])) / (2.0 * (eh + ev));
          }
          // first pixel
          if (2-(a+i+1)%2==2) {
              image->r[i][0] = (image->r[i+1][0] + image->r[i-1][0] + image->r[i][1]) / 3;
              image->g[i][0] = (image->g[i+1][0] + image->g[i-1][0] + image->g[i][1]) / 3;
              image->b[i][0] = (image->b[i+1][0] + image->b[i-1][0] + image->b[i][1]) / 3;
          }
          // last pixel
          if (2-(a+i+image->width)%2==2) {
              image->r[i][image->width-1] = (image->r[i+1][image->width-1] + image->r[i-1][image->width-1] + image->r[i][image->width-2]) / 3;
              image->g[i][image->width-1] = (image->g[i+1][image->width-1] + image->g[i-1][image->width-1] + image->g[i][image->width-2]) / 3;
              image->b[i][image->width-1] = (image->b[i+1][image->width-1] + image->b[i-1][image->width-1] + image->b[i][image->width-2]) / 3;
          }
        }
        // last row
        int b = (a==1 && image->height%2) || (a==0 && image->height%2==0);
        for (int j=1+b; j<image->width-1; j+=2) {
          image->r[image->height-1][j] = (image->r[image->height-2][j] + image->r[image->height-1][j+1] + image->r[image->height-1][j-1]) / 3;
          image->g[image->height-1][j] = (image->g[image->height-2][j] + image->g[image->height-1][j+1] + image->g[image->height-1][j-1]) / 3;
          image->b[image->height-1][j] = (image->b[image->height-2][j] + image->b[image->height-1][j+1] + image->b[image->height-1][j-1]) / 3;
        }
    }


	
    // Flip if needed
    if (tran & TR_HFLIP)
        hflip (image);
    if (tran & TR_VFLIP)
        vflip (image);
        
    // Color correction (only when running on full resolution)
    if (ri->isBayer() && pp.skip==1)
        processFalseColorCorrection (image, raw.ccSteps);
    colorSpaceConversion (image, cmp, embProfile, camProfile, xyz_cam, (static_cast<const ImageData*>(getMetaData()))->getCamera(), defGain);
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
 */
int RawImageSource::findHotDeadPixel( PixelsMap &bpMap, float thresh)
{
	volatile int counter=0;
	
	float (*cfablur);
	cfablur = (float (*)) calloc (H*W, sizeof *cfablur);
	
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
				   rawData[inext][jprev],rawData[inext][j],rawData[inext][jnext],cfablur[i*W+j]);
		}
	}
	
#pragma omp  for
	//cfa pixel heat/death evaluation
	for (int rr=0; rr < H; rr++) {
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		for (int cc=0; cc < W; cc++) {
			//rawData[rr][cc] = cfablur[rr*W+cc];//diagnostic

			//evaluate pixel for heat/death
			float pixdev = fabs(rawData[rr][cc]-cfablur[rr*W+cc]);
			float hfnbrave=0;
			int top=max(0,rr-2);
			int bottom=min(H-1,rr+2);
			int left=max(0,cc-2);
			int right=min(W-1,cc+2);
			for (int mm=top; mm<=bottom; mm++)
				for (int nn=left; nn<=right; nn++) {
					hfnbrave += fabs(rawData[mm][nn]-cfablur[mm*W+nn]);
				}
			hfnbrave = (hfnbrave-pixdev)/((bottom-top+1)*(right-left+1)-1);
			
			if (pixdev > thresh*hfnbrave) {
				// mark the pixel as "bad"
				bpMap.set(cc,rr );
				counter++;
			}
		}//end of pixel evaluation
		
		
	}
	}//end pragma
	free (cfablur);
	//printf ("counter %d \n",counter);
	return counter;
}

void RawImageSource::rotateLine (float* line, float** channel, int tran, int i, int w, int h) {

    if ((tran & TR_ROT) == TR_R180) 
        for (int j=0; j<w; j++) 
            channel[h-1-i][w-1-j] = line[j];

    else if ((tran & TR_ROT) == TR_R90) 
        for (int j=0; j<w; j++) 
            channel[j][h-1-i] = line[j];

    else if ((tran & TR_ROT) == TR_R270) 
        for (int j=0; j<w; j++) 
            channel[w-1-j][i] = line[j];
    else 
		for (int j=0; j<w; j++) 
			channel[i][j] = line[j];
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
                image->r[image->height-1-y][image->width-1-x] = red[j];
                image->g[image->height-1-y][image->width-1-x] = green[j];
                image->b[image->height-1-y][image->width-1-x] = blue[j];
            }
          }
      }
      else if ((tran & TR_ROT) == TR_R270) {
          int end = min(h+fw-i, w-fw+i);
          for (int j=start; j<end; j++) {
            int y = i+j-fw;
            int x = fw-i+j;
            if (x>=0 && x<image->height && y>=0 && y<image->width) {
                image->r[image->height-1-x][y] = red[j];
                image->g[image->height-1-x][y] = green[j];
                image->b[image->height-1-x][y] = blue[j];
            }
          }
      }
      else if ((tran & TR_ROT) == TR_R90) {
          int end = min(h+fw-i, w-fw+i);
          for (int j=start; j<end; j++) {
            int y = i+j-fw;
            int x = fw-i+j;
            if (x>=0 && y<image->width && y>=0 && x<image->height) {
                image->r[x][image->width-1-y] = red[j];
                image->g[x][image->width-1-y] = green[j];
                image->b[x][image->width-1-y] = blue[j];
            }
          }
      }
      else {
        int end = min(h+fw-i, w-fw+i);
        for (int j=start; j<end; j++) {
            int y = i+j-fw;
            int x = fw-i+j;
            if (x>=0 && y<image->height && y>=0 && x<image->width) {
                image->r[y][x] = red[j];
                image->g[y][x] = green[j];
                image->b[y][x] = blue[j];
            }
        }
      }
  }
  // Nikon D1X vertical interpolation + coarse rotation
  else if (d1x) {
    // copy new pixels
    if ((tran & TR_ROT) == TR_R180) {
      for (int j=0; j<imwidth; j++) {
        image->r[2*imheight-2-2*i][imwidth-1-j] = red[j];
        image->g[2*imheight-2-2*i][imwidth-1-j] = green[j];
        image->b[2*imheight-2-2*i][imwidth-1-j] = blue[j];
      }

      if (i==1 || i==2) { // linear interpolation
        int row = 2*imheight-1-2*i;
        for (int j=0; j<imwidth; j++) {
          int col = imwidth-1-j;
          image->r[row][col] = (red[j] + image->r[row+1][col]) /2;
          image->g[row][col] = (green[j] + image->g[row+1][col]) /2;
          image->b[row][col] = (blue[j] + image->b[row+1][col]) /2;
        }
      }
      else if (i==imheight-1) {
        int row = 2*imheight-1-2*i;
        for (int j=0; j<imwidth; j++) {
          int col = imwidth-1-j;
          image->r[row][col] = (red[j] + image->r[row+1][col]) /2;
          image->g[row][col] = (green[j] + image->g[row+1][col]) /2;
          image->b[row][col] = (blue[j] + image->b[row+1][col]) /2;
        }
        row = 2*imheight-1-2*i+2;
        for (int j=0; j<imwidth; j++) {
          int col = imwidth-1-j;
          image->r[row][col] = (red[j] + image->r[row+1][col]) /2;
          image->g[row][col] = (green[j] + image->g[row+1][col]) /2;
          image->b[row][col] = (blue[j] + image->b[row+1][col]) /2;
        }
      }
      else if (i>2 && i<imheight-1) { // vertical bicubic interpolationi
        int row = 2*imheight-1-2*i+2;
        for (int j=0; j<imwidth; j++) {
          int col = imwidth-1-j;
          image->r[row][col] = CLIP((int)(-0.0625*red[j] + 0.5625*image->r[row-1][col] + 0.5625*image->r[row+1][col] - 0.0625*image->r[row+3][col]));
          image->g[row][col] = CLIP((int)(-0.0625*green[j] + 0.5625*image->g[row-1][col] + 0.5625*image->g[row+1][col] - 0.0625*image->g[row+3][col]));
          image->b[row][col] = CLIP((int)(-0.0625*blue[j] + 0.5625*image->b[row-1][col] + 0.5625*image->b[row+1][col] - 0.0625*image->b[row+3][col]));
        }
      }
    }
    else if ((tran & TR_ROT) == TR_R90) {
      for (int j=0; j<imwidth; j++) {
        image->r[j][2*imheight-2-2*i] = red[j];
        image->g[j][2*imheight-2-2*i] = green[j];
        image->b[j][2*imheight-2-2*i] = blue[j];
      }
      if (i==1 || i==2) { // linear interpolation
        int col = 2*imheight-1-2*i;
        for (int j=0; j<imwidth; j++) {
          image->r[j][col] = (red[j] + image->r[j][col+1]) /2;
          image->g[j][col] = (green[j] + image->g[j][col+1]) /2;
          image->b[j][col] = (blue[j] + image->b[j][col+1]) /2;
        }
      }
      else if (i==imheight-1) {
        int col = 2*imheight-1-2*i;
        for (int j=0; j<imwidth; j++) {
          image->r[j][col] = (red[j] + image->r[j][col+1]) /2;
          image->g[j][col] = (green[j] + image->g[j][col+1]) /2;
          image->b[j][col] = (blue[j] + image->b[j][col+1]) /2;
        }
        col = 2*imheight-1-2*i+2;
        for (int j=0; j<imwidth; j++) {
          image->r[j][col] = (red[j] + image->r[j][col+1]) /2;
          image->g[j][col] = (green[j] + image->g[j][col+1]) /2;
          image->b[j][col] = (blue[j] + image->b[j][col+1]) /2;
        }
      }
      else if (i>2 && i<imheight-1) { // vertical bicubic interpolation
        int col = 2*imheight-1-2*i+2;
        for (int j=0; j<imwidth; j++) {
          image->r[j][col] = CLIP((int)(-0.0625*red[j] + 0.5625*image->r[j][col-1] + 0.5625*image->r[j][col+1] - 0.0625*image->r[j][col+3]));
          image->g[j][col] = CLIP((int)(-0.0625*green[j] + 0.5625*image->g[j][col-1] + 0.5625*image->g[j][col+1] - 0.0625*image->g[j][col+3]));
          image->b[j][col] = CLIP((int)(-0.0625*blue[j] + 0.5625*image->b[j][col-1] + 0.5625*image->b[j][col+1] - 0.0625*image->b[j][col+3]));
        }
      }
    }
    else if ((tran & TR_ROT) == TR_R270) {
      for (int j=0; j<imwidth; j++) {
        image->r[imwidth-1-j][2*i] = red[j];
        image->g[imwidth-1-j][2*i] = green[j];
        image->b[imwidth-1-j][2*i] = blue[j];
      }
      if (i==1 || i==2) { // linear interpolation
        for (int j=0; j<imwidth; j++) {
          int row = imwidth-1-j;
          image->r[row][2*i-1] = (red[j] + image->r[row][2*i-2]) * 0.5f;
          image->g[row][2*i-1] = (green[j] + image->g[row][2*i-2]) * 0.5f;
          image->b[row][2*i-1] = (blue[j] + image->b[row][2*i-2]) * 0.5f;
        }
      }
      else if (i==imheight-1) {
        for (int j=0; j<imwidth; j++) {
          int row = imwidth-1-j;
          image->r[row][2*i-1] = (red[j] + image->r[row][2*i-2]) * 0.5f;
          image->g[row][2*i-1] = (green[j] + image->g[row][2*i-2]) * 0.5f;
          image->b[row][2*i-1] = (blue[j] + image->b[row][2*i-2]) * 0.5f;
          image->r[row][2*i-3] = (image->r[row][2*i-2] + image->r[row][2*i-4]) * 0.5f;
          image->g[row][2*i-3] = (image->g[row][2*i-2] + image->g[row][2*i-4]) * 0.5f;
          image->b[row][2*i-3] = (image->b[row][2*i-2] + image->b[row][2*i-4]) * 0.5f;
        }
      }
      else if (i>0 && i<imheight-1) { // vertical bicubic interpolationi
        for (int j=0; j<imwidth; j++) {
          int row = imwidth-1-j;
          image->r[row][2*i-3] = CLIP((int)(-0.0625*red[j] + 0.5625*image->r[row][2*i-2] + 0.5625*image->r[row][2*i-4] - 0.0625*image->r[row][2*i-6]));
          image->g[row][2*i-3] = CLIP((int)(-0.0625*green[j] + 0.5625*image->g[row][2*i-2] + 0.5625*image->g[row][2*i-4] - 0.0625*image->g[row][2*i-6]));
          image->b[row][2*i-3] = CLIP((int)(-0.0625*blue[j] + 0.5625*image->b[row][2*i-2] + 0.5625*image->b[row][2*i-4] - 0.0625*image->b[row][2*i-6]));
        }
      }    
    }
    else {
        rotateLine (red, image->r, tran, 2*i, imwidth, imheight);
        rotateLine (green, image->g, tran, 2*i, imwidth, imheight);
        rotateLine (blue, image->b, tran, 2*i, imwidth, imheight);

      if (i==1 || i==2) { // linear interpolation
        for (int j=0; j<imwidth; j++) {
          image->r[2*i-1][j] = (red[j] + image->r[2*i-2][j]) /2;
          image->g[2*i-1][j] = (green[j] + image->g[2*i-2][j]) /2;
          image->b[2*i-1][j] = (blue[j] + image->b[2*i-2][j]) /2;
        }
      }
      else if (i==imheight-1) {
            for (int j=0; j<imwidth; j++) {
              image->r[2*i-3][j] = (image->r[2*i-4][j] + image->r[2*i-2][j]) /2;
              image->g[2*i-3][j] = (image->g[2*i-4][j] + image->g[2*i-2][j]) /2;
              image->b[2*i-3][j] = (image->b[2*i-4][j] + image->b[2*i-2][j]) /2;
              image->r[2*i-1][j] = (red[j] + image->r[2*i-2][j]) /2;
              image->g[2*i-1][j] = (green[j] + image->g[2*i-2][j]) /2;
              image->b[2*i-1][j] = (blue[j] + image->b[2*i-2][j]) /2;
            }
      }
      else if (i>2 && i<imheight-1) { // vertical bicubic interpolationi
        for (int j=0; j<imwidth; j++) {
          image->r[2*i-3][j] = CLIP((int)(-0.0625*red[j] + 0.5625*image->r[2*i-2][j] + 0.5625*image->r[2*i-4][j] - 0.0625*image->r[2*i-6][j]));
          image->g[2*i-3][j] = CLIP((int)(-0.0625*green[j] + 0.5625*image->g[2*i-2][j] + 0.5625*image->g[2*i-4][j] - 0.0625*image->g[2*i-6][j]));
          image->b[2*i-3][j] = CLIP((int)(-0.0625*blue[j] + 0.5625*image->b[2*i-2][j] + 0.5625*image->b[2*i-4][j] - 0.0625*image->b[2*i-6][j]));
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
    int width  = image->width;
    int height = image->height;

    float* rowr = new float[width];
    float* rowg = new float[width];
    float* rowb = new float[width];
    for (int i=0; i<height; i++) {
      for (int j=0; j<width; j++) {
        rowr[j] = image->r[i][width-1-j];
        rowg[j] = image->g[i][width-1-j];
        rowb[j] = image->b[i][width-1-j];
      }
      memcpy (image->r[i], rowr, width*sizeof(float));
      memcpy (image->g[i], rowg, width*sizeof(float));
      memcpy (image->b[i], rowb, width*sizeof(float));
    }
    delete [] rowr;
    delete [] rowg;
    delete [] rowb;
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::vflip (Imagefloat* image) {
    int width  = image->width;
    int height = image->height;

    float tmp;
    for (int i=0; i<height/2; i++) 
      for (int j=0; j<width; j++) {
        tmp = image->r[i][j]; 
        image->r[i][j] = image->r[height-1-i][j];
        image->r[height-1-i][j] = tmp;
        tmp = image->g[i][j]; 
        image->g[i][j] = image->g[height-1-i][j];
        image->g[height-1-i][j] = tmp;
        tmp = image->b[i][j]; 
        image->b[i][j] = image->b[height-1-i][j];
        image->b[height-1-i][j] = tmp;
      }
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
    int errCode = ri->loadRaw ();
    if (errCode) return errCode;

    ri->compress_image();
    if (plistener) {
        plistener->setProgress (0.8);
    }
/***** Copy once constant data extracted from raw *******/
    W = ri->get_width();
    H = ri->get_height();
    fuji = ri->get_FujiWidth()!=0;
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            rgb_cam[i][j] = ri->get_rgb_cam(i,j);
    // compute inverse of the color transformation matrix
	// first arg is matrix, second arg is inverse
    inverse33 (rgb_cam, cam_rgb);

    d1x  = ! ri->get_model().compare("D1X");
    if (d1x)
        border = 8;
    if ( ri->get_profile() )
        embProfile = cmsOpenProfileFromMem (ri->get_profile(), ri->get_profileLen());

    // create profile
    memset (xyz_cam, 0, sizeof(xyz_cam));
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            for (int k=0; k<3; k++)
                xyz_cam[i][j] += xyz_sRGB[i][k] * rgb_cam[k][j]; 
    camProfile = iccStore->createFromMatrix (xyz_cam, false, "Camera");
    inverse33 (xyz_cam, cam_xyz);

	float pre_mul[4];
	ri->get_colorsCoeff( pre_mul, scale_mu_l, c_black);//modify  for black level

	camwb_red = ri->get_pre_mul(0) / pre_mul[0];
	camwb_green = ri->get_pre_mul(1) / pre_mul[1];
	camwb_blue = ri->get_pre_mul(2) / pre_mul[2];
	initialGain = 1.0 / min(pre_mul[0], pre_mul[1], pre_mul[2]);

    double cam_r = rgb_cam[0][0]*camwb_red + rgb_cam[0][1]*camwb_green + rgb_cam[0][2]*camwb_blue;
    double cam_g = rgb_cam[1][0]*camwb_red + rgb_cam[1][1]*camwb_green + rgb_cam[1][2]*camwb_blue;
    double cam_b = rgb_cam[2][0]*camwb_red + rgb_cam[2][1]*camwb_green + rgb_cam[2][2]*camwb_blue;

    wb = ColorTemp (cam_r, cam_g, cam_b);

    ri->set_prefilters();

    //Load complete Exif informations
    RawMetaDataLocation rml;
    rml.exifBase = ri->get_exifBase();
    rml.ciffBase = ri->get_ciffBase();
    rml.ciffLength = ri->get_ciffLen();
    idata = new ImageData (fname, &rml);

    green = allocArray<float>(W,H);
    red   = allocArray<float>(W,H);
    blue  = allocArray<float>(W,H);
    //hpmap = allocArray<char>(W, H);

    if (plistener) {
        plistener->setProgress (1.0);
    }
    plistener=NULL; // This must be reset, because only load() is called through progressConnector
    t2.set();
    if( settings->verbose )
       printf("Load %s: %d µsec\n",fname.c_str(), t2.etime(t1));

    return 0; // OK!
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::preprocess  (const RAWParams &raw)
{
	MyTime t1,t2;
	t1.set();
	Glib::ustring newDF = raw.dark_frame;
	RawImage *rid=NULL;
	if (!raw.df_autoselect) {
		if( !raw.dark_frame.empty())
			rid = dfm.searchDarkFrame( raw.dark_frame );
	} else {
		rid = dfm.searchDarkFrame( ri->get_maker(), ri->get_model(), ri->get_ISOspeed(), ri->get_shutter(), ri->get_timestamp());
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
	if( rif && settings->verbose) {
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
		bp = dfm.getHotPixels( ri->get_maker(), ri->get_model(), ri->get_ISOspeed(), ri->get_shutter(), ri->get_timestamp());
	}else if( !raw.dark_frame.empty() )
		bp = dfm.getHotPixels( raw.dark_frame );
	if(bp){
		totBP+=bitmapBads.set( *bp );
		if( settings->verbose && !bp->empty()){
			std::cout << "Correcting " << bp->size() << " hotpixels from darkframe" << std::endl;
		}
	}

    scaleColors( 0,0, W, H, raw);//+ + raw parameters for black level(raw.blackxx)

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
    if ( raw.greenthresh || (((idata->getMake().size()>=7 && idata->getMake().substr(0,7)=="OLYMPUS" && idata->getModel()[0]=='E') || (idata->getMake().size()>=9 && idata->getMake().substr(0,9)=="Panasonic")) && raw.dmethod != RAWParams::methodstring[ RAWParams::vng4] && ri->isBayer()) ) {
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

	if ( raw.greenthresh >0) {
		if (plistener) {
			plistener->setProgressStr ("Green equilibrate...");
			plistener->setProgress (0.0);
		}
		green_equilibrate(0.01*(raw.greenthresh));
    }

	
	if ( raw.linenoise >0 ) {
		if (plistener) {
			plistener->setProgressStr ("Line Denoise...");
			plistener->setProgress (0.0);
		}

		cfa_linedn(0.00002*(raw.linenoise));
	}
	
	if ( raw.ca_autocorrect || fabs(raw.cared)>0.001 || fabs(raw.cablue)>0.001 ) {
		if (plistener) {
			plistener->setProgressStr ("CA Auto Correction...");
			plistener->setProgress (0.0);
		}
		
		CA_correct_RT(raw.cared, raw.cablue);
	}
	
	if ( raw.expos !=1 ) processRawWhitepoint(raw.expos, raw.preser);
	
    t2.set();
    if( settings->verbose )
       printf("Preprocessing: %d usec\n", t2.etime(t1));
    return;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
void RawImageSource::demosaic(const RAWParams &raw)
{
    if (ri->isBayer()) {
    	MyTime t1,t2;
    	t1.set();
        if ( raw.dmethod == RAWParams::methodstring[RAWParams::hphd] )
                hphd_demosaic ();
        else if (raw.dmethod == RAWParams::methodstring[RAWParams::vng4] )
            vng4_demosaic ();
        else if (raw.dmethod == RAWParams::methodstring[RAWParams::ahd] )
            ahd_demosaic (0,0,W,H);
	    else if (raw.dmethod == RAWParams::methodstring[RAWParams::amaze] )
            amaze_demosaic_RT (0,0,W,H);
        else if (raw.dmethod == RAWParams::methodstring[RAWParams::dcb] )
            dcb_demosaic(raw.dcb_iterations, raw.dcb_enhance);
        else if (raw.dmethod == RAWParams::methodstring[RAWParams::eahd])
            eahd_demosaic ();
        else if (raw.dmethod == RAWParams::methodstring[RAWParams::fast] )
            fast_demosaic (0,0,W,H);
			//nodemosaic();//for testing
		else
        	nodemosaic();
        t2.set();
        if( settings->verbose )
           printf("Demosaicing: %s - %d usec\n",raw.dmethod.c_str(), t2.etime(t1));
		
        if (raw.all_enhance) refinement_lassus();
        
        rgbSourceModified = false;
    }
}

void RawImageSource::flushRawData() {
    if(cache) {
        delete [] cache;
        cache = 0;
    }
    if (rawData) {
        freeArray<float>(rawData, H);
        rawData = 0;
    }
}

void RawImageSource::flushRGB() {
    if (green) {
        freeArray<float>(green, H);
        green = 0;
    }
    if (red) {
        freeArray<float>(red, H);
        red = 0;
    }
    if (blue) {
        freeArray<float>(blue, H);
        blue = 0;
    }
}

void RawImageSource::HLRecovery_Global(HRecParams hrp )
{
	//color propagation highlight recovery 
	if (hrp.enabled && hrp.method=="Color"){
		if (settings->verbose) printf ("Applying Highlight Recovery: Color propagation...\n");
		HLRecovery_inpaint (red,green,blue);
		rgbSourceModified = true;
	}
	else{
		rgbSourceModified = false;
	}
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/* Copy original pixel data and
 * subtract dark frame (if present) from current image and apply flat field correction (if present)
 */
void RawImageSource::copyOriginalPixels(const RAWParams &raw, RawImage *src, RawImage *riDark, RawImage *riFlatFile )
{
	if (ri->isBayer()) {
		if (!rawData)
			rawData = allocArray<float>(W,H);
		if (riDark && W == riDark->get_width() && H == riDark->get_height()) {
			for (int row = 0; row < H; row++) {
				for (int col = 0; col < W; col++) {
					rawData[row][col]	= max(src->data[row][col]+ri->get_black() - riDark->data[row][col], 0);
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
			
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			float (*cfablur);
			cfablur = (float (*)) calloc (H*W, sizeof *cfablur);
//#define BS 32	
			int BS = raw.ff_BlurRadius;
			if (BS&1) BS++;
			
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
			
			float refcolor[2][2],vignettecorr;
			//find center ave values by channel
			for (int m=0; m<2; m++)
				for (int n=0; n<2; n++) {
					refcolor[m][n] = max(0.0f,cfablur[(2*(H>>2)+m)*W+2*(W>>2)+n] - ri->get_black());
				}
			
			for (int m=0; m<2; m++)
				for (int n=0; n<2; n++) {
					for (int row = 0; row+m < H; row+=2) 
						for (int col = 0; col+n < W; col+=2) {
							vignettecorr = ( refcolor[m][n]/max(1e-5f,cfablur[(row+m)*W+col+n]-ri->get_black()) );
							rawData[row+m][col+n] = (rawData[row+m][col+n]-ri->get_black()) * vignettecorr + ri->get_black(); 	
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
				
				float vlinecorr, hlinecorr;
				
				for (int m=0; m<2; m++)
					for (int n=0; n<2; n++) {
						for (int row = 0; row+m < H; row+=2) 
							for (int col = 0; col+n < W; col+=2) {
								hlinecorr = (max(1e-5f,cfablur[(row+m)*W+col+n]-ri->get_black())/max(1e-5f,cfablur1[(row+m)*W+col+n]-ri->get_black()) );
								vlinecorr = (max(1e-5f,cfablur[(row+m)*W+col+n]-ri->get_black())/max(1e-5f,cfablur2[(row+m)*W+col+n]-ri->get_black()) );
								rawData[row+m][col+n] = ((rawData[row+m][col+n]-ri->get_black()) * hlinecorr * vlinecorr + ri->get_black()); 
							}
					}
				free (cfablur1);
				free (cfablur2);
			}
			
			free (cfablur);
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//#undef BS
			

		}  // flatfield
	} else {
        // No bayer pattern
        // TODO: Is there a flat field correction possible?
		if (!rawData) rawData = allocArray<float>(3*W,H);

		if (riDark && W == riDark->get_width() && H == riDark->get_height()) {
			for (int row = 0; row < H; row++) {
				for (int col = 0; col < W; col++) {
					rawData[row][3*col+0] = max(src->data[row][3*col+0]+ri->get_black() - riDark->data[row][3*col+0], 0);
					rawData[row][3*col+1] = max(src->data[row][3*col+1]+ri->get_black() - riDark->data[row][3*col+1], 0);
					rawData[row][3*col+2] = max(src->data[row][3*col+2]+ri->get_black() - riDark->data[row][3*col+2], 0);
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

void RawImageSource::cfaboxblur(RawImage *riFlatFile, float* cfablur, int boxH, int boxW ) {

	float (*cfatmp);
	cfatmp = (float (*)) calloc (H*W, sizeof *cfatmp);
	float hotdeadthresh = 0.5;
		
	for (int i=0; i<H; i++) {
		int iprev,inext,jprev,jnext;
		float p[9],temp, median;
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
			if (riFlatFile->data[i][j]>hotdeadthresh*median || median>hotdeadthresh*riFlatFile->data[i][j]) {
				cfatmp[i*W+j] = median;
			} else {
				cfatmp[i*W+j] = riFlatFile->data[i][j];
			}

		}
	}
	
	//box blur cfa image; box size = BS
	//horizontal blur
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
	free (cfatmp);
	
}
	

// Scale original pixels into the range 0 65535 using black offsets and multipliers 
void RawImageSource::scaleColors(int winx,int winy,int winw,int winh, const RAWParams &raw)
{
float black_lev[4];//black level
	
	chmax[0]=chmax[1]=chmax[2]=chmax[3]=0;//channel maxima

//adjust black level  (eg Canon)
// cblack Bayer
//0 green
//1 red
//2 blue
// cblack no Bayer
//0 red
//1 green
//2 blue
if( ri->isBayer() ) {
black_lev[0]=raw.blackzero;//G1
black_lev[1]=raw.blackone;//R
black_lev[2]=raw.blacktwo;//B
black_lev[3]=raw.blackthree;//G2

}
else {
black_lev[0]=raw.blackone;//R
black_lev[1]=raw.blackzero;//G
black_lev[2]=raw.blacktwo;//B
black_lev[3]= raw.blackzero;
}
  for(int i=0; i<4; i++) {
   scale_mul[i]=scale_mu_l[i];}
   
  if( c_black[0]+black_lev[1] >0) cblacksom[0]=c_black[0]+black_lev[1]; else cblacksom[0]=0;// adjust black level
  if( c_black[3]+black_lev[3] >0) cblacksom[3]=c_black[3]+black_lev[3]; else cblacksom[3]=0;// adjust black level
  if( c_black[2]+black_lev[2] >0) cblacksom[2]=c_black[2]+black_lev[2]; else cblacksom[2]=0;// adjust black level
  if( c_black[1]+black_lev[0] >0) cblacksom[1]=c_black[1]+black_lev[0]; else cblacksom[1]=0;// adjust black level
// this seems strange, but it works

		// scale image colors
		
	if( ri->isBayer() ){
		for (int row = winy; row < winy+winh; row ++){
			for (int col = winx; col < winx+winw; col++) {
				float val = rawData[row][col];				
				int c = FC(row, col);
				if (ri->ISGREEN(row,col)) {
                    if (row&1) {
						val-=cblacksom[1]; 
						val *= scale_mul[1];
                    }
                    else {
                        val-=cblacksom[3];
						val *= scale_mul[3];
                    }
                }			
				else if (ri->ISRED(row,col)) {
				val-=cblacksom[0];
				val*=scale_mul[0];}	
				else if (ri->ISBLUE(row,col)) {
				val-=cblacksom[2];
				val*=scale_mul[2];}	
				
				rawData[row][col] = (val);
				chmax[c] = max(chmax[c],val);
			}
		}
	}else{
	// i don't know how it's run...
		for (int row = winy; row < winy+winh; row ++){
			for (int col = winx; col < winx+winw; col++) {
				float val = rawData[row][3*col+0];
				if (val){
					val -= cblack[0];
					val *= scale_mul[0];
					rawData[row][3*col+0] = (val);
					chmax[0] = max(chmax[0],val);
				}
				val = rawData[row][3*col+1];
				if (val){
					val -= cblack[1];
					val *= scale_mul[1];
					rawData[row][3*col+1] = (val);
					chmax[1] = max(chmax[1],val);
				}
				val = rawData[row][3*col+2];
				if (val){
					val -= cblack[2];
					val *= scale_mul[2];
					rawData[row][3*col+2] = (val);
					chmax[2] = max(chmax[2],val);
				}
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

  int ppx=0, px=(row_from-1)%3, cx=row_from%3, nx=0;
  
    convert_row_to_YIQ (im->r[row_from-1], im->g[row_from-1], im->b[row_from-1], rbconv_Y[px], rbconv_I[px], rbconv_Q[px], W);
    convert_row_to_YIQ (im->r[row_from], im->g[row_from], im->b[row_from], rbconv_Y[cx], rbconv_I[cx], rbconv_Q[cx], W);

    for (int j=0; j<W; j++) {
      rbout_I[px][j] = rbconv_I[px][j];
      rbout_Q[px][j] = rbconv_Q[px][j];
    }

    for (int i=row_from; i<row_to; i++) {
       
      ppx = (i-2)%3;
      px = (i-1)%3;
      cx = i%3;
      nx = (i+1)%3;

      convert_row_to_YIQ (im->r[i+1], im->g[i+1], im->b[i+1], rbconv_Y[nx], rbconv_I[nx], rbconv_Q[nx], W);

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
        convert_row_to_RGB (im->r[i-1], im->g[i-1], im->b[i-1], rbconv_Y[px], row_I, row_Q, W);
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
    convert_row_to_RGB (im->r[row_to-1], im->g[row_to-1], im->b[row_to-1], rbconv_Y[cx], row_I, row_Q, W);

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

// Converts raw image including ICC input profile to working space - floating point version
void RawImageSource::colorSpaceConversion (Imagefloat* im, ColorManagementParams cmp, cmsHPROFILE embedded, cmsHPROFILE camprofile, double camMatrix[3][3], std::string camName, double& defgain) {

    //MyTime t1, t2, t3;
    //t1.set ();
    cmsHPROFILE in;
    DCPProfile *dcpProf;

    if (!findInputProfile(cmp.input, embedded, camName, &dcpProf, in)) return;

    if (dcpProf!=NULL) {
        dcpProf->Apply(im, cmp.working);
    } else {
    // Calculate matrix for direct conversion raw>working space
        TMatrix work = iccStore->workingSpaceInverseMatrix (cmp.working);
        double mat[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (int i=0; i<3; i++)
            for (int j=0; j<3; j++) 
                for (int k=0; k<3; k++) 
                    mat[i][j] += work[i][k] * camMatrix[k][j]; // rgb_xyz * xyz_cam


    if (in==NULL) {
		// use default camprofile, supplied by dcraw
        // in this case we avoid using the slllllooooooowwww lcms

        #pragma omp parallel for
        for (int i=0; i<im->height; i++)
            for (int j=0; j<im->width; j++) {

                float newr = mat[0][0]*im->r[i][j] + mat[0][1]*im->g[i][j] + mat[0][2]*im->b[i][j];
                float newg = mat[1][0]*im->r[i][j] + mat[1][1]*im->g[i][j] + mat[1][2]*im->b[i][j];
                float newb = mat[2][0]*im->r[i][j] + mat[2][1]*im->g[i][j] + mat[2][2]*im->b[i][j];

                im->r[i][j] = (newr);
                im->g[i][j] = (newg);
                im->b[i][j] = (newb);

            }
    } else {
        Imagefloat* imgPreLCMS=NULL;
        if (cmp.blendCMSMatrix) imgPreLCMS=im->copy();

        // use supplied input profile
		// color space transform is expecting data in the range (0,1)
        #pragma omp parallel for
		for ( int h = 0; h < im->height; ++h )
			for ( int w = 0; w < im->width; ++w ) {
				im->r[h][w] /= 65535.0;
				im->g[h][w] /= 65535.0;
				im->b[h][w] /= 65535.0;
			}


        // Gamma preprocessing
        float gammaFac, lineFac, lineSum;
        getProfilePreprocParams(in, gammaFac, lineFac, lineSum);

        if (gammaFac>0) {
            #pragma omp parallel for
		    for ( int h = 0; h < im->height; ++h )
			    for ( int w = 0; w < im->width; ++w ) {
				    im->r[h][w] = pow (max(im->r[h][w],0.0f), gammaFac);
				    im->g[h][w] = pow (max(im->g[h][w],0.0f), gammaFac);
				    im->b[h][w] = pow (max(im->b[h][w],0.0f), gammaFac);
			    }
        }


	if(settings->gamutICC)  
	// use Prophoto to apply profil ICC, then converted to cmp.working (sRGB..., Adobe.., Wide..) to avoid color shifts due to relative colorimetric
	// LCMS use intent for applying profil => suppression of negatives values and > 65535
	//useful for correcting Munsell Lch

	{ if( settings->verbose ) printf("With Gamut ICC correction float\n");
		Glib::ustring profi ="ProPhoto";
        cmsHPROFILE out = iccStore->workingSpace (profi);//Prophoto	
	TMatrix wprof = iccStore->workingSpaceMatrix (profi);
	TMatrix wiprof = iccStore->workingSpaceInverseMatrix (cmp.working);//sRGB .. Adobe...Wide...
	double toxyz[3][3] = {
        {
        	( wprof[0][0]),
        	( wprof[0][1]),
        	( wprof[0][2])
        },{
			( wprof[1][0]		),
			( wprof[1][1]		),
			( wprof[1][2]		)
        },{
			( wprof[2][0]),
			( wprof[2][1]),
			( wprof[2][2])
        }
		};
        lcmsMutex->lock ();
        cmsHTRANSFORM hTransform = cmsCreateTransform (in, TYPE_RGB_FLT, out, TYPE_RGB_FLT, 
            INTENT_RELATIVE_COLORIMETRIC,  // float is clipless, so don't trim it
            cmsFLAGS_NOCACHE );  // NOCACHE is important for thread safety
        lcmsMutex->unlock ();
		if (hTransform) {
            im->ExecCMSTransform(hTransform);
			}
			else {
          // create the profile from camera
          lcmsMutex->lock ();
          hTransform = cmsCreateTransform (camprofile, TYPE_RGB_FLT, out, TYPE_RGB_FLT, settings->colorimetricIntent,
              cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE );  // NOCACHE is important for thread safety    
          lcmsMutex->unlock ();
				
          im->ExecCMSTransform(hTransform);
				}
		float x, y,z;
		Glib::ustring choiceprofile;
		choiceprofile=cmp.working;
		if(choiceprofile!="ProPhoto") {
		for ( int h = 0; h < im->height; ++h )
			for ( int w = 0; w < im->width; ++w ) {//convert from Prophoto to XYZ
		     x = (toxyz[0][0] * im->r[h][w] + toxyz[0][1] * im->g[h][w]  + toxyz[0][2] * im->b[h][w] ) ;
             y = (toxyz[1][0] * im->r[h][w] + toxyz[1][1] * im->g[h][w]  + toxyz[1][2] * im->b[h][w] ) ;
             z = (toxyz[2][0] * im->r[h][w] + toxyz[2][1] * im->g[h][w]  + toxyz[2][2] * im->b[h][w] ) ;
			 //convert from XYZ to cmp.working  (sRGB...Adobe...Wide..)
			im->r[h][w] = ((wiprof[0][0]*x + wiprof[0][1]*y + wiprof[0][2]*z)) ;
			im->g[h][w] = ((wiprof[1][0]*x + wiprof[1][1]*y + wiprof[1][2]*z)) ;
			im->b[h][w] = ((wiprof[2][0]*x + wiprof[2][1]*y + wiprof[2][2]*z)) ;			
			}	
			}
			
		        cmsDeleteTransform(hTransform);
	
	}
     else {	
	 if( settings->verbose ) printf("Without Gamut ICC correction float\n");
		cmsHPROFILE out = iccStore->workingSpace (cmp.working);	

//        out = iccStore->workingSpaceGamma (wProfile);
        lcmsMutex->lock ();
        cmsHTRANSFORM hTransform = cmsCreateTransform (in, TYPE_RGB_FLT, out, TYPE_RGB_FLT, 
            INTENT_RELATIVE_COLORIMETRIC,  // float is clipless, so don't trim it
            cmsFLAGS_NOCACHE );  // NOCACHE is important for thread safety
        lcmsMutex->unlock ();

        if (hTransform) {
            // there is an input profile
            im->ExecCMSTransform(hTransform);
        } else {
          // create the profile from camera
          lcmsMutex->lock ();
          hTransform = cmsCreateTransform (camprofile, TYPE_RGB_FLT, out, TYPE_RGB_FLT, settings->colorimetricIntent,
              cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE );  // NOCACHE is important for thread safety    
          lcmsMutex->unlock ();

          im->ExecCMSTransform(hTransform);
        }

        cmsDeleteTransform(hTransform);
		}
		
		// restore normalization to the range (0,65535) and blend matrix colors if LCMS is clipping
        const float RecoverTresh = 65535.0 * 0.98;  // just the last few percent highlights are merged
        #pragma omp parallel for
		for ( int h = 0; h < im->height; ++h )
			for ( int w = 0; w < im->width; ++w ) {
               
                // There might be Nikon postprocessings
                if (lineSum>0) {
                    im->r[h][w] *= im->r[h][w] * lineFac + lineSum;
                    im->g[h][w] *= im->g[h][w] * lineFac + lineSum;
                    im->b[h][w] *= im->b[h][w] * lineFac + lineSum;
                }

				im->r[h][w] *= 65535.0 ;
				im->g[h][w] *= 65535.0 ;
				im->b[h][w] *= 65535.0 ;

                if (cmp.blendCMSMatrix) {
                    // Red
                    float red=im->r[h][w];
                    if (red>RecoverTresh) {
                        float matrixRed = mat[0][0]*imgPreLCMS->r[h][w] + mat[0][1]*imgPreLCMS->g[h][w] + mat[0][2]*imgPreLCMS->b[h][w];

                        if (red>=65535.0)
                            im->r[h][w] = matrixRed;
                        else {
                            float fac = (red - RecoverTresh) / (65535.0 - RecoverTresh);
                            im->r[h][w] = (1.0-fac) * red + fac * matrixRed;
                        }
                    }

                    // Green
                    float green=im->g[h][w];
                    if (green>RecoverTresh) {
                        float matrixGreen = mat[1][0]*imgPreLCMS->r[h][w] + mat[1][1]*imgPreLCMS->g[h][w] + mat[1][2]*imgPreLCMS->b[h][w];

                        if (green>=65535.0)
                            im->g[h][w] = matrixGreen;
                        else {
                            float fac = (green - RecoverTresh) / (65535.0 - RecoverTresh);
                            im->g[h][w] = (1.0-fac) * green + fac * matrixGreen;
                        }
                    }


                    // Blue
                    float blue=im->b[h][w];
                    if (blue>RecoverTresh) {
                        float matrixBlue = mat[2][0]*imgPreLCMS->r[h][w] + mat[2][1]*imgPreLCMS->g[h][w] + mat[2][2]*imgPreLCMS->b[h][w];

                        if (blue>=65535.0)
                            im->b[h][w] = matrixBlue;
                        else {
                            float fac = (blue - RecoverTresh) / (65535.0 - RecoverTresh);
                            im->b[h][w] = (1.0-fac) * blue + fac * matrixBlue;
                        }
                    }
                }
			}

        if (imgPreLCMS!=NULL) delete imgPreLCMS;
    }
    }
        //t3.set ();
//        printf ("ICM TIME: %d\n", t3.etime(t1));
}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
// Converts raw image including ICC input profile to working space - 16bit int version
void RawImageSource::colorSpaceConversion16 (Image16* im, ColorManagementParams cmp, cmsHPROFILE embedded, cmsHPROFILE camprofile, double camMatrix[3][3], std::string camName, double& defgain) {
	cmsHPROFILE in;
    DCPProfile *dcpProf;
	
    if (!findInputProfile(cmp.input, embedded, camName, &dcpProf, in)) return;
	
    if (dcpProf!=NULL) {
        dcpProf->Apply(im, cmp.working);
    } else {
	if (in==NULL) {
        // Take camprofile from DCRAW
        // in this case we avoid using the slllllooooooowwww lcms
TMatrix work = iccStore->workingSpaceInverseMatrix (cmp.working);
		double mat[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++) 
				for (int k=0; k<3; k++) 
					mat[i][j] += work[i][k] * camMatrix[k][j]; // rgb_xyz * xyz_cam
		
#pragma omp parallel for
		for (int i=0; i<im->height; i++)
			for (int j=0; j<im->width; j++) {
				
				float newr = mat[0][0]*im->r[i][j] + mat[0][1]*im->g[i][j] + mat[0][2]*im->b[i][j];
				float newg = mat[1][0]*im->r[i][j] + mat[1][1]*im->g[i][j] + mat[1][2]*im->b[i][j];
				float newb = mat[2][0]*im->r[i][j] + mat[2][1]*im->g[i][j] + mat[2][2]*im->b[i][j];
				
				im->r[i][j] = CLIP((int)newr);
				im->g[i][j] = CLIP((int)newg);
				im->b[i][j] = CLIP((int)newb);
			}
	}
	else
	{
        // Gamma preprocessing
        float gammaFac, lineFac, lineSum;
        getProfilePreprocParams(in, gammaFac, lineFac, lineSum);

        if (gammaFac>0) {
            #pragma omp parallel for
		    for ( int h = 0; h < im->height; ++h )
			    for ( int w = 0; w < im->width; ++w ) {
				    im->r[h][w]=  (int) (pow ((double)(im->r[h][w] / 65535.0), (double)gammaFac) * 65535.0);
				    im->g[h][w]=  (int) (pow ((double)(im->g[h][w] / 65535.0), (double)gammaFac) * 65535.0);
				    im->b[h][w]=  (int) (pow ((double)(im->b[h][w] / 65535.0), (double)gammaFac) * 65535.0);
			    }
        }

		cmsHPROFILE out = iccStore->workingSpace (cmp.working);
		//        out = iccStore->workingSpaceGamma (wProfile);
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
                        im->r[h][w] *= im->r[h][w] * lineFac / 65535.0 + lineSum;
                        im->g[h][w] *= im->g[h][w] * lineFac / 65535.0 + lineSum;
                        im->b[h][w] *= im->b[h][w] * lineFac / 65535.0 + lineSum;
                    }
            }
		}
		else {
			lcmsMutex->lock ();
			hTransform = cmsCreateTransform (camprofile, TYPE_RGB_16, out, TYPE_RGB_16, settings->colorimetricIntent,
                cmsFLAGS_NOCACHE);   
			lcmsMutex->unlock ();

			im->ExecCMSTransform(hTransform);
		}

		cmsDeleteTransform(hTransform);
	}
    }
	//t3.set ();
	//        printf ("ICM TIME: %d\n", t3.etime(t1));
}

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
#include "rt_math.h"
		
		
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
			
			float clipfrac[3];
			FOREACHCOLOR clipfrac[c] = min(1.0f,rgb[c]/maxave);
			
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
            float fy = (yy<65535.0 ? ImProcFunctions::cachef[yy]/327.68 : (exp(log(yy/MAXVAL)/3.0 )));
            // compute LCH decompostion of the clipped pixel (only color information, thus C and H will be used)
            float x = xyz_cam[0][0]*ro + xyz_cam[0][1]*go + xyz_cam[0][2]*bo;
            float y = xyz_cam[1][0]*ro + xyz_cam[1][1]*go + xyz_cam[1][2]*bo;
            float z = xyz_cam[2][0]*ro + xyz_cam[2][1]*go + xyz_cam[2][2]*bo;
			x = (x<65535.0 ? ImProcFunctions::cachef[x]/327.68 : (exp(log(x/MAXVAL)/3.0 )));
			y = (y<65535.0 ? ImProcFunctions::cachef[y]/327.68 : (exp(log(y/MAXVAL)/3.0 )));
			z = (z<65535.0 ? ImProcFunctions::cachef[z]/327.68 : (exp(log(z/MAXVAL)/3.0 )));
            // convert back to rgb
            double fz = fy - y + z;
            double fx = fy + x - y;

			double zr = ImProcFunctions::f2xyz(fz);
            double xr = ImProcFunctions::f2xyz(fx);

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
        HLRecovery_CIELab (red, green, blue, red, green, blue, width, 65535.0, xyz_cam, cam_xyz);
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

    for (int i=border; i<H-border; i++) {
        int start, end;
        getRowStartEnd (i, start, end);

        if (ri->isBayer()) {
            for (int j=start; j<end; j++) {
                if (ri->ISGREEN(i,j))
                    histogram[CLIP((int)(camwb_green*rawData[i][j]))>>histcompr]+=4;
                else if (ri->ISRED(i,j))
					histogram[CLIP((int)(camwb_red*rawData[i][j]))>>histcompr]+=4;
				else if (ri->ISBLUE(i,j))
					histogram[CLIP((int)(camwb_blue*rawData[i][j]))>>histcompr]+=4;
			} 
			} else {
				for (int j=start; j<3*end; j++) {
                    histogram[CLIP((int)(camwb_red*rawData[i][j+0]))>>histcompr]++;
                    histogram[CLIP((int)(camwb_green*rawData[i][j+1]))>>histcompr]+=2;
                    histogram[CLIP((int)(camwb_blue*rawData[i][j+2]))>>histcompr]++;
				}
			}
    }
}
		
// Histogram MUST be 256 in size; gamma is applied, blackpoint and gain also
void RawImageSource::getRAWHistogram (LUTu & histRedRaw, LUTu & histGreenRaw, LUTu & histBlueRaw) {

    histRedRaw.clear(); histGreenRaw.clear(); histBlueRaw.clear();
	float mult = 65535.0 / ri->get_white();
    #pragma omp parallel for
    for (int i=border; i<H-border; i++) {
        int start, end, idx;
        getRowStartEnd (i, start, end);

	
        if (ri->isBayer()) {
            for (int j=start; j<end; j++) {
                if (ri->ISGREEN(i,j)) {
                    if(i &1) idx = CLIP((int)CurveFactory::gamma(mult*(ri->data[i][j]-(cblacksom[1]/*+black_lev[1]*/))));// green 1
					else 
					idx = CLIP((int)CurveFactory::gamma(mult*(ri->data[i][j]-(cblacksom[3]/*+black_lev[3]*/))));//green 2
                    histGreenRaw[idx>>8]++;
                } else if (ri->ISRED(i,j)) {
                    idx = CLIP((int)CurveFactory::gamma(mult*(ri->data[i][j]-(cblacksom[0]/*+black_lev[0]*/))));
					
                    histRedRaw[idx>>8]++;
                } else if (ri->ISBLUE(i,j)) {
                    idx = CLIP((int)CurveFactory::gamma(mult*(ri->data[i][j]-(cblacksom[2]/*+black_lev[2]*/))));
					
                    histBlueRaw[idx>>8]++;
    }
            }
        } else {
			for (int j=start; j<3*end; j++) {
				idx = CLIP((int)CurveFactory::gamma(mult*(ri->data[i][j]-cblack[0])));
                histRedRaw[idx>>8]++;

				idx = CLIP((int)CurveFactory::gamma(mult*(ri->data[i][j+1]-cblack[1])));
                histGreenRaw[idx>>8]++;

				idx = CLIP((int)CurveFactory::gamma(mult*(ri->data[i][j+2]-cblack[2])));
                histBlueRaw[idx>>8]++;
			}
		}
    }

    // since there are twice as many greens, correct for it
    if (ri->isBayer()) for (int i=0;i<256;i++) histGreenRaw[i]>>=1;
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
	ColorTemp RawImageSource::getAutoWB () {
		
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
					if (!ri->isBayer()) {
						double d = CLIP(initialGain*(rawData[i][3*j]));
						if (d>64000)
							continue;
						avg_r += d; rn++;
						d = CLIP(initialGain*(rawData[i][3*j+1]));
						if (d>64000)
							continue;
						avg_g += d; gn++;
						d = CLIP(initialGain*(rawData[i][3*j+2]));
						if (d>64000)
							continue;
						avg_b += d; bn++;
					}
					else {
						int c = FC( i, j);
						double d = CLIP(initialGain*(rawData[i][j]));
						if (d>64000)
							continue;
						double dp = d;
						if (c==0) {
							avg_r += dp;
							rn++;
						}
						else if (c==1) {
							avg_g += dp;
							gn++;
						}
						else if (c==2) {
							avg_b += dp;
							bn++;
						}
					}
				}
			}
		}
		else {
			if (!ri->isBayer()) {
				for (int i=32; i<H-32; i++)
					for (int j=32; j<W-32; j++) {
						double dr = CLIP(initialGain*(rawData[i][3*j]  ));
						double dg = CLIP(initialGain*(rawData[i][3*j+1]));
						double db = CLIP(initialGain*(rawData[i][3*j+2]));
						if (dr>64000 || dg>64000 || db>64000) continue;
						avg_r += dr; rn++;
						avg_g += dg; 
						avg_b += db; 
					}
				gn = rn; bn=rn;
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
						//average a Bayer quartet if nobody is clipped
						d[0][0] = CLIP(initialGain*(rawData[i][j]    ));
						d[0][1] = CLIP(initialGain*(rawData[i][j+1]  ));
						d[1][0] = CLIP(initialGain*(rawData[i+1][j]  ));
						d[1][1] = CLIP(initialGain*(rawData[i+1][j+1]));
						if ( d[0][0]>64000 || d[0][1]>64000 || d[1][0]>64000 || d[1][1]>64000 ) continue;
						avg_r += d[ey][ex];
						avg_g += d[1-ey][ex] + d[ey][1-ex];
						avg_b += d[1-ey][1-ex];
						rn++;
					}
				gn = 2*rn;
				bn = rn;
			}
		}
		if( settings->verbose )
			printf ("AVG: %g %g %g\n", avg_r/rn, avg_g/gn, avg_b/bn);
		
		//    return ColorTemp (pow(avg_r/rn, 1.0/6.0)*img_r, pow(avg_g/gn, 1.0/6.0)*img_g, pow(avg_b/bn, 1.0/6.0)*img_b);
		
		double reds   = avg_r/rn * camwb_red;
		double greens = avg_g/gn * camwb_green;
		double blues  = avg_b/bn * camwb_blue;
		
		double rm = rgb_cam[0][0]*reds + rgb_cam[0][1]*greens + rgb_cam[0][2]*blues;
		double gm = rgb_cam[1][0]*reds + rgb_cam[1][1]*greens + rgb_cam[1][2]*blues;
		double bm = rgb_cam[2][0]*reds + rgb_cam[2][1]*greens + rgb_cam[2][2]*blues;
		
		return ColorTemp (rm, gm, bm);
	}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	ColorTemp RawImageSource::getSpotWB (std::vector<Coord2D> red, std::vector<Coord2D> green, std::vector<Coord2D>& blue, int tran) {
		
		int x; int y;
		double reds = 0, greens = 0, blues = 0;
		int rn = 0;
		
		if (!ri->isBayer()) {
			int xmin, xmax, ymin, ymax;
			int xr, xg, xb, yr, yg, yb;
			for (int i=0; i<red.size(); i++) {
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
			
		} else {
			
			int d[9][2] = {{0,0}, {-1,-1}, {-1,0}, {-1,1}, {0,-1}, {0,1}, {1,-1}, {1,0}, {1,1}};
			int rloc, gloc, bloc, rnbrs, gnbrs, bnbrs;
			for (int i=0; i<red.size(); i++) {
				transformPosition (red[i].x, red[i].y, tran, x, y);
				rloc=gloc=bloc=rnbrs=gnbrs=bnbrs=0;
				for (int k=0; k<9; k++) {
					int xv = x + d[k][0];
					int yv = y + d[k][1];
					int c = FC(yv,xv);
					if (c==0 && xv>=0 && yv>=0 && xv<W && yv<H) { //RED
						rloc += (rawData[yv][xv]);
						rnbrs++;
						continue;
					}else if (c==2 && xv>=0 && yv>=0 && xv<W && yv<H) { //BLUE
						bloc += (rawData[yv][xv]);
						bnbrs++;
						continue;
					} else { // GREEN
						gloc += (rawData[yv][xv]);
						gnbrs++;
						continue;
					}
					
				}
				rloc /= rnbrs; gloc /= gnbrs; bloc /= bnbrs;
				if (rloc*initialGain<64000 && gloc*initialGain<64000 && bloc*initialGain<64000) {
					reds += rloc; greens += gloc; blues += bloc; rn++;
				}
				//transformPosition (green[i].x, green[i].y, tran, x, y);//these are redundant now ??? if not, repeat for these blocks same as for red[]
				//transformPosition (blue[i].x, blue[i].y, tran, x, y);
			}
		}
		
		if (2*rn < red.size()) {
			return ColorTemp ();
		}
		else {
			reds = reds/rn * camwb_red;
			greens = greens/rn * camwb_green;
			blues = blues/rn * camwb_blue;
			
			double rm = rgb_cam[0][0]*reds + rgb_cam[0][1]*greens + rgb_cam[0][2]*blues;
			double gm = rgb_cam[1][0]*reds + rgb_cam[1][1]*greens + rgb_cam[1][2]*blues;
			double bm = rgb_cam[2][0]*reds + rgb_cam[2][1]*greens + rgb_cam[2][2]*blues;
			
			return ColorTemp (rm, gm, bm);
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
	
	
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//#include "demosaic_algos.cc"
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Emil's code 
#include "fast_demo.cc"//fast demosaic	
#include "amaze_demosaic_RT.cc"//AMaZE demosaic	
#include "CA_correct_RT.cc"//Emil's CA auto correction
#include "cfa_linedn_RT.cc"//Emil's line denoise
#include "green_equil_RT.cc"//Emil's green channel equilibration
#include "hilite_recon.cc"//Emil's highlight reconstruction

#include "expo_before_b.cc"//Jacques's exposure before interpolation
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef PIX_SORT
#undef med3x3

} /* namespace */

#undef PIX_SORT
#undef med3x3

