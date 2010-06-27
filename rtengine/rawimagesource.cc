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
#include <rawimagesource.h>
#include <median.h>
#include <common.h>
#include <math.h>
#include <mytime.h>
#include <iccmatrices.h>
#include <iccstore.h>
#include <image8.h>
#include <curves.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace rtengine {

int loadRaw (const char* fname, struct RawImage* ri);

extern const Settings* settings;

#undef ABS
#undef MAX
#undef MIN
#undef DIST
#undef MAXVAL
#undef CLIP

#define ABS(a) ((a)<0?-(a):(a))
#define MAX(a,b) ((a)<(b)?(b):(a))
#define MIN(a,b) ((a)>(b)?(b):(a))
#define DIST(a,b) (ABS(a-b))
#define MAXVAL  0xffff
#define CLIP(a) ((a)>0?((a)<MAXVAL?(a):MAXVAL):0)

RawImageSource::RawImageSource () : ImageSource(), plistener(NULL), green(NULL), red(NULL), blue(NULL), cache(NULL), border(4) {

    hrmap[0] = NULL;
    hrmap[1] = NULL;
    hrmap[2] = NULL;
	needhr = NULL;
    hpmap = NULL;
	oldmethod = "None";
	camProfile = NULL;
	embProfile = NULL;
}

RawImageSource::~RawImageSource () {

    delete idata;
    if (ri) {
        if (ri->allocation)
            free (ri->allocation);
        if (ri->data)
            free (ri->data);
        if (ri->profile_data)
            free (ri->profile_data);
        delete ri;
    }
    if (green)
        freeArray<unsigned short>(green, H);
    if (red)
        freeArray<unsigned short>(red, H);
    if (blue)
        freeArray<unsigned short>(blue, H);
    
    delete [] cache;
    if (hrmap[0]!=NULL) {
        int dh = H/HR_SCALE;
        freeArray<float>(hrmap[0], dh);
        freeArray<float>(hrmap[1], dh);
        freeArray<float>(hrmap[2], dh);
    }
    if (needhr)
        freeArray<char>(needhr, H);
    if (hpmap)
        freeArray<char>(hpmap, H);
    if (camProfile)
        cmsCloseProfile (camProfile);
    if (embProfile)
        cmsCloseProfile (embProfile);
}

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
        w = ri->fuji_width * 2 + 1;
        h = (H - ri->fuji_width)*2 + 1;
    }
    
    int sw = w, sh = h;  
    if ((tran & TR_ROT) == TR_R90 || (tran & TR_ROT) == TR_R270) {
        sw = h;
        sh = w;
    }
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
        ssx1 = (sx1+sy1) / 2;
        ssy1 = (sy1 - sx2 ) / 2 + ri->fuji_width;
        int ssx2 = (sx2+sy2) / 2 + 1;
        int ssy2 = (sy2 - sx1) / 2 + ri->fuji_width;
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

void RawImageSource::getImage (ColorTemp ctemp, int tran, Image16* image, PreviewProps pp, HRecParams hrp, ColorManagementParams cmp) {

    isrcMutex.lock ();

    tran = defTransform (tran);

    // compute channel multipliers
    double r, g, b, rm, gm, bm;
    ctemp.getMultipliers (r, g, b);
    rm = icoeff[0][0]*r + icoeff[0][1]*g + icoeff[0][2]*b;
    gm = icoeff[1][0]*r + icoeff[1][1]*g + icoeff[1][2]*b;
    bm = icoeff[2][0]*r + icoeff[2][1]*g + icoeff[2][2]*b;
    rm = ri->camwb_red / rm;
    gm = ri->camwb_green / gm;
    bm = ri->camwb_blue / bm;
    double mul_lum = 0.299*rm + 0.587*gm + 0.114*bm;
    rm /= mul_lum;
    gm /= mul_lum;
    bm /= mul_lum;    

    if (hrp.enabled) 
        defGain = log(ri->defgain) / log(2.0);
    else {
        defGain = 0.0;
        rm *= ri->defgain;
        gm *= ri->defgain;
        bm *= ri->defgain;
    }

	if (hrp.enabled==true && hrp.method=="Color" && hrmap[0]==NULL) 
        updateHLRecoveryMap_ColorPropagation ();

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
       
    // render the requested image part
    unsigned short* red  = new unsigned short[imwidth];
    unsigned short* grn  = new unsigned short[imwidth];
    unsigned short* blue = new unsigned short[imwidth];

    for (int i=sy1,ix=0; ix<imheight; i+=pp.skip, ix++) {
        if (ri->filters && this->red && this->blue) {
            for (int j=0,jx=sx1; j<imwidth; j++,jx+=pp.skip) {
                red[j] = CLIP(rm*this->red[i][jx]);
                grn[j] = CLIP(gm*green[i][jx]);
                blue[j] = CLIP(bm*this->blue[i][jx]);
	    }
	} else if(ri->filters) {
            if (i==0)
                interpolate_row_rb_mul_pp (red, blue, NULL, green[i], green[i+1], i, rm, gm, bm, sx1, imwidth, pp.skip);
            else if (i==H-1)
                interpolate_row_rb_mul_pp (red, blue, green[i-1], green[i], NULL, i, rm, gm, bm, sx1, imwidth, pp.skip);
            else 
                interpolate_row_rb_mul_pp (red, blue, green[i-1], green[i], green[i+1], i, rm, gm, bm, sx1, imwidth, pp.skip);
                   
            for (int j=0,jx=sx1; j<imwidth; j++,jx+=pp.skip) 
                grn[j] = CLIP(gm*green[i][jx]);
        } else {
            for (int j=0,jx=sx1; j<imwidth; j++,jx+=pp.skip) {
                red[j]  = CLIP(rm*ri->data[i][jx*3+0]);
                grn[j]  = CLIP(gm*ri->data[i][jx*3+1]);
                blue[j] = CLIP(bm*ri->data[i][jx*3+2]);
            }
        }
                
        if (hrp.enabled)
            hlRecovery (hrp.method, red, grn, blue, i, sx1, imwidth, pp.skip);

        transLine (red, grn, blue, ix, image, tran, imwidth, imheight, fw);
    }

    delete [] red;
    delete [] grn;
    delete [] blue;
          
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
        
    // Color correction
    if (ri->filters && pp.skip==1)
        correction_YIQ_LQ (image, settings->colorCorrectionSteps);
 
    // Applying postmul
    colorSpaceConversion (image, cmp, embProfile, camProfile, cam, defGain);

    isrcMutex.unlock ();
}
	
	
void RawImageSource::cfa_clean(float thresh) //Emil's hot/dead pixel removal -- only filters egregiously impulsive pixels
	{  
		// local variables
		int rr, cc;
		int gin, g[8];
		
		float eps=1e-10;//tolerance to avoid dividing by zero
		float p[8];
		float pixave, pixratio;

		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//The cleaning algorithm starts here
				
		for (rr=4; rr < H-4; rr++)
			for (cc=4; cc < W-4; cc++) {
				
				//pixel neighbor average
				gin=ri->data[rr][cc];
				g[0]=ri->data[rr-2][cc-2];
				g[1]=ri->data[rr-2][cc];
				g[2]=ri->data[rr-2][cc+2];
				g[3]=ri->data[rr][cc-2];
				g[4]=ri->data[rr][cc+2];
				g[5]=ri->data[rr+2][cc-2];
				g[6]=ri->data[rr+2][cc];
				g[7]=ri->data[rr+2][cc+2];
				
				pixave=(float)(g[0]+g[1]+g[2]+g[3]+g[4]+g[5]+g[6]+g[7])/8;
				pixratio=MIN(gin,pixave)/(eps+MAX(gin,pixave));
				
				if (pixratio > thresh) continue;
				
				p[0]=1/(eps+fabs(g[0]-gin)+fabs(g[0]-ri->data[rr-4][cc-4])+fabs(ri->data[rr-1][cc-1]-ri->data[rr-3][cc-3]));
				p[1]=1/(eps+fabs(g[1]-gin)+fabs(g[1]-ri->data[rr-4][cc])+fabs(ri->data[rr-1][cc]-ri->data[rr-3][cc]));
				p[2]=1/(eps+fabs(g[2]-gin)+fabs(g[2]-ri->data[rr-4][cc+4])+fabs(ri->data[rr-1][cc+1]-ri->data[rr-3][cc+3]));
				p[3]=1/(eps+fabs(g[3]-gin)+fabs(g[3]-ri->data[rr][cc-4])+fabs(ri->data[rr][cc-1]-ri->data[rr][cc-3]));
				p[4]=1/(eps+fabs(g[4]-gin)+fabs(g[4]-ri->data[rr][cc+4])+fabs(ri->data[rr][cc+1]-ri->data[rr][cc+3]));
				p[5]=1/(eps+fabs(g[5]-gin)+fabs(g[5]-ri->data[rr+4][cc-4])+fabs(ri->data[rr+1][cc-1]-ri->data[rr+3][cc-3]));
				p[6]=1/(eps+fabs(g[6]-gin)+fabs(g[6]-ri->data[rr+4][cc])+fabs(ri->data[rr+1][cc]-ri->data[rr+3][cc]));
				p[7]=1/(eps+fabs(g[7]-gin)+fabs(g[7]-ri->data[rr+4][cc+4])+fabs(ri->data[rr+1][cc+1]-ri->data[rr+3][cc+3]));
				
				ri->data[rr][cc] = (int)((g[0]*p[0]+g[1]*p[1]+g[2]*p[2]+g[3]*p[3]+g[4]*p[4]+ \
										g[5]*p[5]+g[6]*p[6]+g[7]*p[7])/(p[0]+p[1]+p[2]+p[3]+p[4]+p[5]+p[6]+p[7]));
				
			}//now impulsive values have been corrected
		
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	}
	

void RawImageSource::rotateLine (unsigned short* line, unsigned short** channel, int tran, int i, int w, int h) {

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
        memcpy (channel[i], line,  w*sizeof(unsigned short));
}

void RawImageSource::transLine (unsigned short* red, unsigned short* green, unsigned short* blue, int i, Image16* image, int tran, int imwidth, int imheight, int fw) {

  // Fuji SuperCCD rotation + coarse rotation
  if (fuji) {
      int start = ABS(fw-i);
      int w = fw * 2 + 1;
      int h = (imheight - fw)*2 + 1;

      if ((tran & TR_ROT) == TR_R180) {
          int end = MIN(h+fw-i, w-fw+i);
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
          int end = MIN(h+fw-i, w-fw+i);
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
          int end = MIN(h+fw-i, w-fw+i);
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
        int end = MIN(h+fw-i, w-fw+i);
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
          image->r[row][col] = (red[j] + image->r[row+1][col]) >> 1;
          image->g[row][col] = (green[j] + image->g[row+1][col]) >> 1;
          image->b[row][col] = (blue[j] + image->b[row+1][col]) >> 1;
        }
      }
      else if (i==imheight-1) {
        int row = 2*imheight-1-2*i;
        for (int j=0; j<imwidth; j++) {
          int col = imwidth-1-j;
          image->r[row][col] = (red[j] + image->r[row+1][col]) >> 1;
          image->g[row][col] = (green[j] + image->g[row+1][col]) >> 1;
          image->b[row][col] = (blue[j] + image->b[row+1][col]) >> 1;
        }
        row = 2*imheight-1-2*i+2;
        for (int j=0; j<imwidth; j++) {
          int col = imwidth-1-j;
          image->r[row][col] = (red[j] + image->r[row+1][col]) >> 1;
          image->g[row][col] = (green[j] + image->g[row+1][col]) >> 1;
          image->b[row][col] = (blue[j] + image->b[row+1][col]) >> 1;
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
          image->r[j][col] = (red[j] + image->r[j][col+1]) >> 1;
          image->g[j][col] = (green[j] + image->g[j][col+1]) >> 1;
          image->b[j][col] = (blue[j] + image->b[j][col+1]) >> 1;
        }
      }
      else if (i==imheight-1) {
        int col = 2*imheight-1-2*i;
        for (int j=0; j<imwidth; j++) {
          image->r[j][col] = (red[j] + image->r[j][col+1]) >> 1;
          image->g[j][col] = (green[j] + image->g[j][col+1]) >> 1;
          image->b[j][col] = (blue[j] + image->b[j][col+1]) >> 1;
        }
        col = 2*imheight-1-2*i+2;
        for (int j=0; j<imwidth; j++) {
          image->r[j][col] = (red[j] + image->r[j][col+1]) >> 1;
          image->g[j][col] = (green[j] + image->g[j][col+1]) >> 1;
          image->b[j][col] = (blue[j] + image->b[j][col+1]) >> 1;
        }
      }
      else if (i>2 && i<imheight-1) { // vertical bicubic interpolationi
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
          image->r[row][2*i-1] = (red[j] + image->r[row][2*i-2]) >> 1;
          image->g[row][2*i-1] = (green[j] + image->g[row][2*i-2]) >> 1;
          image->b[row][2*i-1] = (blue[j] + image->b[row][2*i-2]) >> 1;
        }
      }
      else if (i==imheight-1) {
        for (int j=0; j<imwidth; j++) {
          int row = imwidth-1-j;
          image->r[row][2*i-1] = (red[j] + image->r[row][2*i-2]) >> 1;
          image->g[row][2*i-1] = (green[j] + image->g[row][2*i-2]) >> 1;
          image->b[row][2*i-1] = (blue[j] + image->b[row][2*i-2]) >> 1;
          image->r[row][2*i-3] = (image->r[row][2*i-2] + image->r[row][2*i-4]) >> 1;
          image->g[row][2*i-3] = (image->g[row][2*i-2] + image->g[row][2*i-4]) >> 1;
          image->b[row][2*i-3] = (image->b[row][2*i-2] + image->b[row][2*i-4]) >> 1;
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
          image->r[2*i-1][j] = (red[j] + image->r[2*i-2][j]) >> 1;
          image->g[2*i-1][j] = (green[j] + image->g[2*i-2][j]) >> 1;
          image->b[2*i-1][j] = (blue[j] + image->b[2*i-2][j]) >> 1;
        }
      }
      else if (i==imheight-1) {
            for (int j=0; j<imwidth; j++) {
              image->r[2*i-3][j] = (image->r[2*i-4][j] + image->r[2*i-2][j]) >> 1;
              image->g[2*i-3][j] = (image->g[2*i-4][j] + image->g[2*i-2][j]) >> 1;
              image->b[2*i-3][j] = (image->b[2*i-4][j] + image->b[2*i-2][j]) >> 1;
              image->r[2*i-1][j] = (red[j] + image->r[2*i-2][j]) >> 1;
              image->g[2*i-1][j] = (green[j] + image->g[2*i-2][j]) >> 1;
              image->b[2*i-1][j] = (blue[j] + image->b[2*i-2][j]) >> 1;
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
  }
  // other (conventional) CCD coarse rotation
  else {
    rotateLine (red, image->r, tran, i, imwidth, imheight);
    rotateLine (green, image->g, tran, i, imwidth, imheight);
    rotateLine (blue, image->b, tran, i, imwidth, imheight);
  }
}

void RawImageSource::getFullSize (int& w, int& h, int tr) {

    tr = defTransform (tr);

    if (fuji) {
        w = ri->fuji_width * 2 + 1;
        h = (H - ri->fuji_width)*2 + 1;
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

void RawImageSource::hflip (Image16* image) {
    int width  = image->width;
    int height = image->height;

    unsigned short* rowr = new unsigned short[width];
    unsigned short* rowg = new unsigned short[width];
    unsigned short* rowb = new unsigned short[width];
    for (int i=0; i<height; i++) {
      for (int j=0; j<width; j++) {
        rowr[j] = image->r[i][width-1-j];
        rowg[j] = image->g[i][width-1-j];
        rowb[j] = image->b[i][width-1-j];
      }
      memcpy (image->r[i], rowr, width*sizeof(unsigned short));
      memcpy (image->g[i], rowg, width*sizeof(unsigned short));
      memcpy (image->b[i], rowb, width*sizeof(unsigned short));
    }
    delete [] rowr;
    delete [] rowg;
    delete [] rowb;
}

void RawImageSource::vflip (Image16* image) {
    int width  = image->width;
    int height = image->height;

    register unsigned short tmp;
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

void RawImageSource::inverse33 (double (*coeff)[3], double (*icoeff)[3]) {
    double nom = coeff[0][2]*coeff[1][1]*coeff[2][0] - coeff[0][1]*coeff[1][2]*coeff[2][0] - coeff[0][2]*coeff[1][0]*coeff[2][1] + coeff[0][0]*coeff[1][2]*coeff[2][1] + coeff[0][1]*coeff[1][0]*coeff[2][2] - coeff[0][0]*coeff[1][1]*coeff[2][2];
    icoeff[0][0] = (coeff[1][2]*coeff[2][1]-coeff[1][1]*coeff[2][2]) / nom;
    icoeff[0][1] = -(coeff[0][2]*coeff[2][1]-coeff[0][1]*coeff[2][2]) / nom;
    icoeff[0][2] = (coeff[0][2]*coeff[1][1]-coeff[0][1]*coeff[1][2]) / nom;
    icoeff[1][0] = -(coeff[1][2]*coeff[2][0]-coeff[1][0]*coeff[2][2]) / nom;
    icoeff[1][1] = (coeff[0][2]*coeff[2][0]-coeff[0][0]*coeff[2][2]) / nom;
    icoeff[1][2] = -(coeff[0][2]*coeff[1][0]-coeff[0][0]*coeff[1][2]) / nom;
    icoeff[2][0] = (coeff[1][1]*coeff[2][0]-coeff[1][0]*coeff[2][1]) / nom;
    icoeff[2][1] = -(coeff[0][1]*coeff[2][0]-coeff[0][0]*coeff[2][1]) / nom;
    icoeff[2][2] = (coeff[0][1]*coeff[1][0]-coeff[0][0]*coeff[1][1]) / nom;
}
    
int RawImageSource::load (Glib::ustring fname) {

    fileName = fname;

    if (plistener) {
        plistener->setProgressStr ("Decoding...");
        plistener->setProgress (0.0);
    }

    ri = new RawImage;
    int res = loadRaw (fname.c_str(), ri);
    if (res)
        return res;

   if(red) {
	   delete red;
	   red = 0;
   }
   if(green) {
	   delete green;
	   green = 0;
   }
   if(blue) {
	   delete blue;
	   blue = 0;
   }

    W = ri->width;
    H = ri->height;

    d1x  = !strcmp(ri->model, "D1X");
    fuji = ri->fuji_width;
    if (d1x)
        border = 8;

    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            coeff[i][j] = ri->coeff[i][j];

    // compute inverse of the color transformation matrix
    inverse33 (coeff, icoeff);

    double cam_r = coeff[0][0]*ri->camwb_red + coeff[0][1]*ri->camwb_green + coeff[0][2]*ri->camwb_blue;
    double cam_g = coeff[1][0]*ri->camwb_red + coeff[1][1]*ri->camwb_green + coeff[1][2]*ri->camwb_blue;
    double cam_b = coeff[2][0]*ri->camwb_red + coeff[2][1]*ri->camwb_green + coeff[2][2]*ri->camwb_blue;

    wb = ColorTemp (cam_r, cam_g, cam_b);

    double tr = icoeff[0][0] * cam_r + icoeff[0][1] * cam_g + icoeff[0][2] * cam_b;
    double tg = icoeff[1][0] * cam_r + icoeff[1][1] * cam_g + icoeff[1][2] * cam_b;
    double tb = icoeff[2][0] * cam_r + icoeff[2][1] * cam_g + icoeff[2][2] * cam_b;

    // create profile
    memset (cam, 0, sizeof(cam));
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            for (int k=0; k<3; k++)
                cam[i][j] += coeff[k][i] * sRGB_d50[k][j];
    camProfile = iccStore.createFromMatrix (cam, false, "Camera");
    inverse33 (cam, icam);

    if (ri->profile_data)
        embProfile = cmsOpenProfileFromMem (ri->profile_data, ri->profile_len);

    defGain = log(ri->defgain) / log(2.0); 
    
    RawMetaDataLocation rml;
    rml.exifBase = ri->exifbase;
    rml.ciffBase = ri->ciff_base;
    rml.ciffLength = ri->ciff_len;

    idata = new ImageData (fname, &rml); 

    // check if it is an olympus E camera, if yes, compute G channel pre-compensation factors
    if (settings->greenthresh || (((idata->getMake().size()>=7 && idata->getMake().substr(0,7)=="OLYMPUS" && idata->getModel()[0]=='E') || (idata->getMake().size()>=9 && idata->getMake().substr(0,7)=="Panasonic")) && settings->demosaicMethod!="vng4" && ri->filters) ) {
        // global correction
        int ng1=0, ng2=0;
        double avgg1=0, avgg2=0;
        for (int i=border; i<H-border; i++)
            for (int j=border; j<W-border; j++)
                if (ISGREEN(ri,i,j)) {
                    if (i%2==0) {
                        avgg1 += ri->data[i][j];
                        ng1++;
                    }
                    else {
                        avgg2 += ri->data[i][j];
                        ng2++;
                    }
                }
        double corrg1 = ((double)avgg1/ng1 + (double)avgg2/ng2) / 2.0 / ((double)avgg1/ng1);
        double corrg2 = ((double)avgg1/ng1 + (double)avgg2/ng2) / 2.0 / ((double)avgg2/ng2);
        for (int i=border; i<H-border; i++)
            for (int j=border; j<W-border; j++)
                if (ISGREEN(ri,i,j)) 
                        ri->data[i][j] = CLIP(ri->data[i][j] * (i%2 ? corrg2 : corrg1));
	}
	

	// local correction in a 9x9 box
	if (settings->greenthresh) {
//Emil's green equilbration		
		if (plistener) {
			plistener->setProgressStr ("Green equilibrate...");
			plistener->setProgress (0.0);
		}
		green_equilibrate(0.01*(settings->greenthresh));

/*        unsigned short* corr_alloc = new unsigned short[W*H];
        unsigned short** corr_data = new unsigned short* [H];
        for (int i=0; i<H; i++) 
            corr_data[i] = corr_alloc + i*W;
        memcpy (corr_alloc, ri->allocation, W*H*sizeof(unsigned short));
        for (int i=border; i<H-border; i++)
            for (int j=border; j<W-border; j++)
                if (ISGREEN(ri,i,j)) {
                    unsigned int ag1 = ri->data[i-4][j-4] + ri->data[i-4][j-2] + ri->data[i-4][j] + ri->data[i-4][j+2] +
                              ri->data[i-2][j-4] + ri->data[i-2][j-2] + ri->data[i-2][j] + ri->data[i-2][j+2] +
                              ri->data[i][j-4] + ri->data[i][j-2] + ri->data[i][j] + ri->data[i][j+2] +
                              ri->data[i+2][j-4] + ri->data[i+2][j-2] + ri->data[i+2][j] + ri->data[i+2][j+2];
                    unsigned int ag2 = ri->data[i-3][j-3] + ri->data[i-3][j-1] + ri->data[i-3][j+1] + ri->data[i-3][j+1] +
                              ri->data[i-1][j-3] + ri->data[i-1][j-1] + ri->data[i-1][j+1] + ri->data[i-1][j+1] +
                              ri->data[i+1][j-3] + ri->data[i+1][j-1] + ri->data[i+1][j+1] + ri->data[i+1][j+1] +
                              ri->data[i+3][j-3] + ri->data[i+3][j-1] + ri->data[i+3][j+1] + ri->data[i+3][j+1];
                    unsigned int val = (ri->data[i][j] + ri->data[i][j] * ag2 / ag1) / 2;
                    corr_data[i][j] = CLIP (val);
                }
        memcpy (ri->allocation, corr_alloc, W*H*sizeof(unsigned short));
        delete corr_alloc;
        delete corr_data; */ 
                       
    }
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Emil's CFA hot/dead pixel filter
	
	if (settings->hotdeadpix_filt) {
		if (plistener) {
			plistener->setProgressStr ("Hot/Dead Pixel Filter...");
			plistener->setProgress (0.0);
		}
		
		cfa_clean(0.1);
	}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Emil's line noise filter
	
	if (settings->linenoise) {
		if (plistener) {
			plistener->setProgressStr ("Line Denoise...");
			plistener->setProgress (0.0);
		}
		
		cfa_linedn(0.00002*(settings->linenoise));
	}
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Emil's CA auto correction
	
	if (settings->ca_autocorrect) {
		if (plistener) {
			plistener->setProgressStr ("CA Auto Correction...");
			plistener->setProgress (0.0);
		}
		
		CA_correct_RT();
	}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (ri->filters) {
    	MyTime t1,t2;
    	t1.set();
        // demosaic
        if (settings->demosaicMethod=="hphd")
            hphd_demosaic ();
        else if (settings->demosaicMethod=="vng4")
            vng4_demosaic ();
        else if (settings->demosaicMethod=="ahd")
            ahd_demosaic ();
        //else if (settings->demosaicMethod=="ppg")
        //    ppg_demosaic ();
		else if (settings->demosaicMethod=="amaze")
            amaze_demosaic_RT ();//Emil's code for AMaZE
        else if (settings->demosaicMethod=="dcb")
            dcb_demosaic(settings->dcb_iterations, settings->dcb_enhance? 1:0);
        else
            eahd_demosaic ();
        t2.set();
        printf("Totale demosaicing:%d \n",t2.etime(t1));
    }


    if (plistener) {
        plistener->setProgressStr ("Ready.");
        plistener->setProgress (1.0);
    }

    return 0;
}

int RawImageSource::defTransform (int tran) {

    int deg = ri->rotate_deg;
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

void RawImageSource::correction_YIQ_LQ_  (Image16* im, int row_from, int row_to) {
 
  int W = im->width;

  int** rbconv_Y = new int*[3];
  int** rbconv_I = new int*[3];
  int** rbconv_Q = new int*[3];
  int** rbout_I = new int*[3];
  int** rbout_Q = new int*[3];
  for (int i=0; i<3; i++) {
    rbconv_Y[i] = new int[W];
    rbconv_I[i] = new int[W];
    rbconv_Q[i] = new int[W];
    rbout_I[i] = new int[W];
    rbout_Q[i] = new int[W];
  }

  int* row_I = new int[W];
  int* row_Q = new int[W];

  int* pre1_I = new int[3];
  int* pre2_I = new int[3];
  int* post1_I = new int[3];
  int* post2_I = new int[3];
  int middle_I[6];
  int* pre1_Q = new int[3];
  int* pre2_Q = new int[3];
  int* post1_Q = new int[3];
  int* post2_Q = new int[3];
  int middle_Q[6];
  int* tmp;

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

  freeArray<int>(rbconv_Y, 3);
  freeArray<int>(rbconv_I, 3);
  freeArray<int>(rbconv_Q, 3);
  freeArray<int>(rbout_I, 3);
  freeArray<int>(rbout_Q, 3);
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


void RawImageSource::correction_YIQ_LQ  (Image16* im, int times) {

    if (im->height<4)
        return;

    for (int t=0; t<times; t++) {
#ifdef _OPENMP
	    #pragma omp parallel
    	{
    		int tid = omp_get_thread_num();
    		int nthreads = omp_get_num_threads();
    		int blk = (im->height-2)/nthreads;

    		if (tid<nthreads-1)
    			correction_YIQ_LQ_ (im, 1 + tid*blk, 1 + (tid+1)*blk);
    		else
    			correction_YIQ_LQ_ (im, 1 + tid*blk, im->height - 1);
    	}
#else
    	correction_YIQ_LQ_ (im, 1 , im->height - 1);
#endif
    }
}



void RawImageSource::convert_row_to_YIQ (unsigned short* r, unsigned short* g, unsigned short* b, int* Y, int* I, int* Q, int W) {
  for (int j=0; j<W; j++) {
    Y[j] = 299 * r[j] + 587 * g[j] + 114 * b[j];
    I[j] = 596 * r[j] - 275 * g[j] - 321 * b[j];
    Q[j] = 212 * r[j] - 523 * g[j] + 311 * b[j];
  }
}


void RawImageSource::convert_row_to_RGB (unsigned short* r, unsigned short* g, unsigned short* b, int* Y, int* I, int* Q, int W) {
  for (int j=1; j<W-1; j++) {
    int ir = Y[j]/1000 + 0.956*I[j]/1000 + 0.621*Q[j]/1000;
    int ig = Y[j]/1000 - 0.272*I[j]/1000 - 0.647*Q[j]/1000;
    int ib = Y[j]/1000 - 1.105*I[j]/1000 + 1.702*Q[j]/1000;
    r[j] = CLIP(ir);
    g[j] = CLIP(ig);
    b[j] = CLIP(ib);
  }
}
#include <curves.h>
void RawImageSource::convert_to_cielab_row (unsigned short* ar, unsigned short* ag, unsigned short* ab, short* oL, short* oa, short* ob) {

  for (int j=0; j<W; j++) {
    double r = ar[j];
    double g = ag[j];
    double b = ab[j];
  
    double x = lc00 * r + lc01 * g + lc02 * b;
    double y = lc10 * r + lc11 * g + lc12 * b;
    double z = lc20 * r + lc21 * g + lc22 * b;

	if (y>threshold)
      oL[j] = 300.0*cache[(int)y];
    else
      oL[j] = 300.0 * 903.3 * y / CMAXVAL;

    oa[j] = 32.0 * 500.0 * ((x>threshold ? cache[(int)x] : 7.787*x/CMAXVAL+16.0/116.0) - (y>threshold ? cache[(int)y] : 7.787*y/CMAXVAL+16.0/116.0));
    ob[j] = 32.0 * 200.0 * ((y>threshold ? cache[(int)y] : 7.787*y/CMAXVAL+16.0/116.0) - (z>threshold ? cache[(int)z] : 7.787*z/CMAXVAL+16.0/116.0));
  }
}

void RawImageSource::interpolate_row_g (unsigned short* agh, unsigned short* agv, int i) {

  for (int j=0; j<W; j++) {
    if (ISGREEN(ri,i,j)) {
      agh[j] = ri->data[i][j];
      agv[j] = ri->data[i][j];
    }
    else {
      int gh=0;
      int gv=0;
      if (j>1 && j<W-2) {
        gh = (-ri->data[i][j-2] + 2*ri->data[i][j-1] + 2*ri->data[i][j] + 2*ri->data[i][j+1] -ri->data[i][j+2]) / 4;
        int maxgh = MAX(ri->data[i][j-1], ri->data[i][j+1]);
        int mingh = MIN(ri->data[i][j-1], ri->data[i][j+1]);
        if (gh>maxgh)
            gh = maxgh;
        else if (gh<mingh)
            gh = mingh;
      }
      else if (j==0)
        gh = ri->data[i][1];
      else if (j==1)
        gh = (ri->data[i][0] + ri->data[i][2]) / 2;
      else if (j==W-1)
        gh = ri->data[i][W-2];
      else if (j==W-2)
        gh = (ri->data[i][W-1] + ri->data[i][W-3]) / 2;
 
     if (i>1 && i<H-2) {
        gv = (-ri->data[i-2][j] + 2*ri->data[i-1][j] + 2*ri->data[i][j] + 2*ri->data[i+1][j] - ri->data[i+2][j]) / 4;
        int maxgv = MAX(ri->data[i-1][j], ri->data[i+1][j]);
        int mingv = MIN(ri->data[i-1][j], ri->data[i+1][j]);
        if (gv>maxgv)
          gv = maxgv;
        else if (gv<mingv)
          gv = mingv;
      }
      else if (i==0)
        gv = ri->data[1][j];
      else if (i==1)
        gv = (ri->data[0][j] + ri->data[2][j]) / 2;
      else if (i==H-1)
        gv = ri->data[H-2][j];
      else if (i==H-2)
        gv = (ri->data[H-1][j] + ri->data[H-3][j]) / 2;

      agh[j] = CLIP(gh);
      agv[j] = CLIP(gv);
    }
  }
}

void RawImageSource::interpolate_row_rb (unsigned short* ar, unsigned short* ab, unsigned short* pg, unsigned short* cg, unsigned short* ng, int i) {
  if (ISRED(ri,i,0) || ISRED(ri,i,1)) {
    // RGRGR or GRGRGR line
    for (int j=0; j<W; j++) {
      if (ISRED(ri,i,j)) {
         // red is simple
         ar[j] = ri->data[i][j];
         // blue: cross interpolation
         int b = 0;
         int n = 0;
         if (i>0 && j>0) {
           b += ri->data[i-1][j-1] - pg[j-1];
           n++;
         }
         if (i>0 && j<W-1) {
           b += ri->data[i-1][j+1] - pg[j+1];
           n++;
         }
         if (i<H-1 && j>0) {
           b += ri->data[i+1][j-1] - ng[j-1];
           n++;
         }
         if (i<H-1 && j<W-1) {
           b += ri->data[i+1][j+1] - ng[j+1];
           n++;
         }
         b = cg[j] + b / n;
         ab[j] = CLIP(b);
      }
      else {
        // linear R-G interp. horizontally
        int r;
        if (j==0)
          r = cg[0] + ri->data[i][1] - cg[1];
        else if (j==W-1)
          r = cg[W-1] + ri->data[i][W-2] - cg[W-2];
        else
          r = cg[j] + (ri->data[i][j-1] - cg[j-1] + ri->data[i][j+1] - cg[j+1]) / 2;
        ar[j] = CLIP(r);
        // linear B-G interp. vertically
        int b;
        if (i==0)
          b = ng[j] + ri->data[1][j] - cg[j];
        else if (i==H-1)
          b = pg[j] + ri->data[H-2][j] - cg[j];
        else
          b = cg[j] + (ri->data[i-1][j] - pg[j] + ri->data[i+1][j] - ng[j]) / 2;
        ab[j] = CLIP(b);
      }
    }
  }
  else {
    // BGBGB or GBGBGB line
    for (int j=0; j<W; j++) {
      if (ISBLUE(ri,i,j)) {
         // red is simple
         ab[j] = ri->data[i][j];
         // blue: cross interpolation
         int r = 0;
         int n = 0;
         if (i>0 && j>0) {
           r += ri->data[i-1][j-1] - pg[j-1];
           n++;
         }
         if (i>0 && j<W-1) {
           r += ri->data[i-1][j+1] - pg[j+1];
           n++;
         }
         if (i<H-1 && j>0) {
           r += ri->data[i+1][j-1] - ng[j-1];
           n++;
         }
         if (i<H-1 && j<W-1) {
           r += ri->data[i+1][j+1] - ng[j+1];
           n++;
         }
         r = cg[j] + r / n;

         ar[j] = CLIP(r);
      }
      else {
        // linear B-G interp. horizontally
        int b;
        if (j==0)
          b = cg[0] + ri->data[i][1] - cg[1];
        else if (j==W-1)
          b = cg[W-1] + ri->data[i][W-2] - cg[W-2];
        else
          b = cg[j] + (ri->data[i][j-1] - cg[j-1] + ri->data[i][j+1] - cg[j+1]) / 2;
        ab[j] = CLIP(b);
        // linear R-G interp. vertically
        int r;
        if (i==0)
          r = ng[j] + ri->data[1][j] - cg[j];
        else if (i==H-1)
          r = pg[j] + ri->data[H-2][j] - cg[j];
        else
          r = cg[j] + (ri->data[i-1][j] - pg[j] + ri->data[i+1][j] - ng[j]) / 2;
        ar[j] = CLIP(r);
      }
    }
  }
}

void RawImageSource::interpolate_row_rb_mul_pp (unsigned short* ar, unsigned short* ab, unsigned short* pg, unsigned short* cg, unsigned short* ng, int i, double r_mul, double g_mul, double b_mul, int x1, int width, int skip) {

  if (ISRED(ri,i,0) || ISRED(ri,i,1)) {
    // RGRGR or GRGRGR line
    for (int j=x1, jx=0; jx<width; j+=skip, jx++) {
      if (ISRED(ri,i,j)) {
         // red is simple
         ar[jx] = CLIP(r_mul * ri->data[i][j]);
         // blue: cross interpolation
         int b = 0;
         int n = 0;
         if (i>0 && j>0) {
           b += b_mul*ri->data[i-1][j-1] - g_mul*pg[j-1];
           n++;
         }
         if (i>0 && j<W-1) {
           b += b_mul*ri->data[i-1][j+1] - g_mul*pg[j+1];
           n++;
         }
         if (i<H-1 && j>0) {
           b += b_mul*ri->data[i+1][j-1] - g_mul*ng[j-1];
           n++;
         }
         if (i<H-1 && j<W-1) {
           b += b_mul*ri->data[i+1][j+1] - g_mul*ng[j+1];
           n++;
         }
         b = g_mul*cg[j] + b / n;
         ab[jx] = CLIP(b);
      }
      else {
        // linear R-G interp. horizontally
        int r;
        if (j==0)
          r = g_mul*cg[0] + r_mul*ri->data[i][1] - g_mul*cg[1];
        else if (j==W-1)
          r = g_mul*cg[W-1] + r_mul*ri->data[i][W-2] - g_mul*cg[W-2];
        else
          r = g_mul*cg[j] + (r_mul*ri->data[i][j-1] - g_mul*cg[j-1] + r_mul*ri->data[i][j+1] - g_mul*cg[j+1]) / 2;
        ar[jx] = CLIP(r);
        // linear B-G interp. vertically
        int b;
        if (i==0)
          b = g_mul*ng[j] + b_mul*ri->data[1][j] - g_mul*cg[j];
        else if (i==H-1)
          b = g_mul*pg[j] + b_mul*ri->data[H-2][j] - g_mul*cg[j];
        else
          b = g_mul*cg[j] + (b_mul*ri->data[i-1][j] - g_mul*pg[j] + b_mul*ri->data[i+1][j] - g_mul*ng[j]) / 2;
        ab[jx] = CLIP(b);
      }
    }
  }
  else {
    // BGBGB or GBGBGB line
    for (int j=x1, jx=0; jx<width; j+=skip, jx++) {
      if (ISBLUE(ri,i,j)) {
         // red is simple
         ab[jx] = CLIP(b_mul*ri->data[i][j]);
         // blue: cross interpolation
         int r = 0;
         int n = 0;
         if (i>0 && j>0) {
           r += r_mul*ri->data[i-1][j-1] - g_mul*pg[j-1];
           n++;
         }
         if (i>0 && j<W-1) {
           r += r_mul*ri->data[i-1][j+1] - g_mul*pg[j+1];
           n++;
         }
         if (i<H-1 && j>0) {
           r += r_mul*ri->data[i+1][j-1] - g_mul*ng[j-1];
           n++;
         }
         if (i<H-1 && j<W-1) {
           r += r_mul*ri->data[i+1][j+1] - g_mul*ng[j+1];
           n++;
         }
         r = g_mul*cg[j] + r / n;

         ar[jx] = CLIP(r);
      }
      else {
        // linear B-G interp. horizontally
        int b;
        if (j==0)
          b = g_mul*cg[0] + b_mul*ri->data[i][1] - g_mul*cg[1];
        else if (j==W-1)
          b = g_mul*cg[W-1] + b_mul*ri->data[i][W-2] - g_mul*cg[W-2];
        else
          b = g_mul*cg[j] + (b_mul*ri->data[i][j-1] - g_mul*cg[j-1] + b_mul*ri->data[i][j+1] - g_mul*cg[j+1]) / 2;
        ab[jx] = CLIP(b);
        // linear R-G interp. vertically
        int r;
        if (i==0)
          r = g_mul*ng[j] + r_mul*ri->data[1][j] - g_mul*cg[j];
        else if (i==H-1)
          r = g_mul*pg[j] + r_mul*ri->data[H-2][j] - g_mul*cg[j];
        else
          r = g_mul*cg[j] + (r_mul*ri->data[i-1][j] - g_mul*pg[j] + r_mul*ri->data[i+1][j] - g_mul*ng[j]) / 2;
        ar[jx] = CLIP(r);
      }
    }
  }
}


void RawImageSource::colorSpaceConversion (Image16* im, ColorManagementParams cmp, cmsHPROFILE embedded, cmsHPROFILE camprofile, double camMatrix[3][3], double& defgain) {

    if (cmp.input == "(none)")
        return;

    MyTime t1, t2, t3;

    t1.set ();

    cmsHPROFILE in;
    cmsHPROFILE out;
    
    Glib::ustring inProfile = cmp.input;

    if (inProfile=="(embedded)") {
        if (embedded)
            in = embedded;
        else
            in = camprofile;
    }
    else if (inProfile=="(camera)" || inProfile=="")
        in = camprofile;
    else {
        in = iccStore.getProfile (inProfile);
        if (in==NULL)
            inProfile = "(camera)";
    }

    
    if (inProfile=="(camera)" || inProfile=="" || (inProfile=="(embedded)" && !embedded)) {
        // in this case we avoid using the slllllooooooowwww lcms
    
//        out = iccStore.workingSpace (wProfile);
//        hTransform = cmsCreateTransform (in, TYPE_RGB_16_PLANAR, out, TYPE_RGB_16_PLANAR, settings->colorimetricIntent, cmsFLAGS_MATRIXINPUT | cmsFLAGS_MATRIXOUTPUT);//cmsFLAGS_MATRIXINPUT | cmsFLAGS_MATRIXOUTPUT);
//        cmsDoTransform (hTransform, im->data, im->data, im->planestride/2);
//        cmsDeleteTransform(hTransform);
        TMatrix work = iccStore.workingSpaceInverseMatrix (cmp.working);
        double mat[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        for (int i=0; i<3; i++)
            for (int j=0; j<3; j++) 
                for (int k=0; k<3; k++) 
                    mat[i][j] += camMatrix[i][k] * work[k][j];
                    
        for (int i=0; i<im->height; i++)
            for (int j=0; j<im->width; j++) {

                int newr = mat[0][0]*im->r[i][j] + mat[1][0]*im->g[i][j] + mat[2][0]*im->b[i][j];
                int newg = mat[0][1]*im->r[i][j] + mat[1][1]*im->g[i][j] + mat[2][1]*im->b[i][j];
                int newb = mat[0][2]*im->r[i][j] + mat[1][2]*im->g[i][j] + mat[2][2]*im->b[i][j];

                im->r[i][j] = CLIP(newr);
                im->g[i][j] = CLIP(newg);
                im->b[i][j] = CLIP(newb);
            }
    }
    else {
        out = iccStore.workingSpace (cmp.working);
//        out = iccStore.workingSpaceGamma (wProfile);
        lcmsMutex->lock ();
        cmsHTRANSFORM hTransform = cmsCreateTransform (in, TYPE_RGB_16_PLANAR, out, TYPE_RGB_16_PLANAR, settings->colorimetricIntent, 0);    
        lcmsMutex->unlock ();
        if (hTransform) {
            if (cmp.gammaOnInput) {
                double gd = pow (2.0, defgain);
                defgain = 0.0;                
                for (int i=0; i<im->height; i++)
                    for (int j=0; j<im->width; j++) {
                        im->r[i][j] = CurveFactory::gamma (CLIP(defgain*im->r[i][j]));
                        im->g[i][j] = CurveFactory::gamma (CLIP(defgain*im->g[i][j]));
                        im->b[i][j] = CurveFactory::gamma (CLIP(defgain*im->b[i][j]));
                    }
            }
            cmsDoTransform (hTransform, im->data, im->data, im->planestride/2);
        }
        else {
          lcmsMutex->lock ();
          hTransform = cmsCreateTransform (camprofile, TYPE_RGB_16_PLANAR, out, TYPE_RGB_16_PLANAR, settings->colorimetricIntent, 0);    
          lcmsMutex->unlock ();
          cmsDoTransform (hTransform, im->data, im->data, im->planestride/2);
        }
        cmsDeleteTransform(hTransform);
    }
        t3.set ();
//        printf ("ICM TIME: %d\n", t3.etime(t1));
}

void RawImageSource::eahd_demosaic () {

  if (plistener) {
    plistener->setProgressStr ("Demosaicing...");
    plistener->setProgress (0.0);
  }

  // prepare chache and constants for cielab conversion
  lc00 = (0.412453 * coeff[0][0] + 0.357580 * coeff[1][0] + 0.180423 * coeff[2][0]) / 0.950456;
  lc01 = (0.412453 * coeff[0][1] + 0.357580 * coeff[1][1] + 0.180423 * coeff[2][1]) / 0.950456;
  lc02 = (0.412453 * coeff[0][2] + 0.357580 * coeff[1][2] + 0.180423 * coeff[2][2]) / 0.950456;

  lc10 = 0.212671 * coeff[0][0] + 0.715160 * coeff[1][0] + 0.072169 * coeff[2][0];
  lc11 = 0.212671 * coeff[0][1] + 0.715160 * coeff[1][1] + 0.072169 * coeff[2][1];
  lc12 = 0.212671 * coeff[0][2] + 0.715160 * coeff[1][2] + 0.072169 * coeff[2][2];

  lc20 = (0.019334 * coeff[0][0] + 0.119193 * coeff[1][0] + 0.950227 * coeff[2][0]) / 1.088754;
  lc21 = (0.019334 * coeff[0][1] + 0.119193 * coeff[1][1] + 0.950227 * coeff[2][1]) / 1.088754;
  lc22 = (0.019334 * coeff[0][2] + 0.119193 * coeff[1][2] + 0.950227 * coeff[2][2]) / 1.088754;

  int maxindex = 2*65536;
  cache = new double[maxindex];
  threshold = (int)(0.008856*CMAXVAL);
  for (int i=0; i<maxindex; i++)
    cache[i] = exp(1.0/3.0 * log((double)i / CMAXVAL));
  
  // end of cielab preparation

  unsigned short* rh[3];
  unsigned short* gh[4];
  unsigned short* bh[3];
  unsigned short* rv[3];
  unsigned short* gv[4];
  unsigned short* bv[3];
  short* lLh[3];
  short* lah[3];
  short* lbh[3];
  short* lLv[3];
  short* lav[3];
  short* lbv[3];
  unsigned short* homh[3];
  unsigned short* homv[3];

  green = new unsigned short*[H];
  for (int i=0; i<H; i++)
    green[i] = new unsigned short[W];

  for (int i=0; i<4; i++) {
    gh[i] = new unsigned short[W];  
    gv[i] = new unsigned short[W];  
  }

  for (int i=0; i<3; i++) {
    rh[i] = new unsigned short[W];  
    bh[i] = new unsigned short[W];  
    rv[i] = new unsigned short[W];  
    bv[i] = new unsigned short[W];  
    lLh[i] = new short[W];  
    lah[i] = new short[W];  
    lbh[i] = new short[W];  
    lLv[i] = new short[W];  
    lav[i] = new short[W];  
    lbv[i] = new short[W];  
    homh[i] = new unsigned short[W];
    homv[i] = new unsigned short[W];
  }   

  // interpolate first two lines
  interpolate_row_g (gh[0], gv[0], 0);
  interpolate_row_g (gh[1], gv[1], 1);
  interpolate_row_g (gh[2], gv[2], 2);
  interpolate_row_rb (rh[0], bh[0], NULL, gh[0], gh[1], 0);
  interpolate_row_rb (rv[0], bv[0], NULL, gv[0], gv[1], 0);
  interpolate_row_rb (rh[1], bh[1], gh[0], gh[1], gh[2], 1);
  interpolate_row_rb (rv[1], bv[1], gv[0], gv[1], gv[2], 1);

  convert_to_cielab_row (rh[0], gh[0], bh[0], lLh[0], lah[0], lbh[0]);
  convert_to_cielab_row (rv[0], gv[0], bv[0], lLv[0], lav[0], lbv[0]);
  convert_to_cielab_row (rh[1], gh[1], bh[1], lLh[1], lah[1], lbh[1]);
  convert_to_cielab_row (rv[1], gv[1], bv[1], lLv[1], lav[1], lbv[1]);


  for (int j=0; j<W; j++) {
    homh[0][j] = 0;
    homv[0][j] = 0;
    homh[1][j] = 0;
    homv[1][j] = 0;
  }

  const static int delta = 1;

  int dLmaph[9];
  int dLmapv[9];
  int dCamaph[9];
  int dCamapv[9];
  int dCbmaph[9];
  int dCbmapv[9];

  for (int i=1; i<H-1; i++) {

    int ix = i%3;
    int imx = (i-1)%3;
    int ipx = (i+1)%3;

    if (i<H-2) {
      interpolate_row_g  (gh[(i+2)%4], gv[(i+2)%4], i+2);
      interpolate_row_rb (rh[(i+1)%3], bh[(i+1)%3], gh[i%4], gh[(i+1)%4], gh[(i+2)%4], i+1);
      interpolate_row_rb (rv[(i+1)%3], bv[(i+1)%3], gv[i%4], gv[(i+1)%4], gv[(i+2)%4], i+1);
    }
    else {
      interpolate_row_rb (rh[(i+1)%3], bh[(i+1)%3], gh[i%4], gh[(i+1)%4], NULL, i+1);
      interpolate_row_rb (rv[(i+1)%3], bv[(i+1)%3], gv[i%4], gv[(i+1)%4], NULL, i+1);
    }
 
    convert_to_cielab_row (rh[(i+1)%3], gh[(i+1)%4], bh[(i+1)%3], lLh[(i+1)%3], lah[(i+1)%3], lbh[(i+1)%3]);
    convert_to_cielab_row (rv[(i+1)%3], gv[(i+1)%4], bv[(i+1)%3], lLv[(i+1)%3], lav[(i+1)%3], lbv[(i+1)%3]);

    for (int j=0; j<W; j++) {
      homh[ipx][j] = 0;
      homv[ipx][j] = 0;
    }
    int sh, sv, idx, idx2;
    for (int j=1; j<W-1; j++) {

      int dmi = 0;
      for (int x=-1; x<=1; x++) {
        idx = (i+x)%3;
        idx2 = (i+x)%2;
        for (int y=-1; y<=1; y++) {
          // compute distance in a, b, and L
          if (dmi<4) {
            sh=homh[idx][j+y];
            sv=homv[idx][j+y];
            if (sh>sv) { // fixate horizontal pixel
              dLmaph[dmi]  = DIST(lLh[ix][j], lLh[idx][j+y]);
              dCamaph[dmi] = DIST(lah[ix][j], lah[idx][j+y]);
              dCbmaph[dmi] = DIST(lbh[ix][j], lbh[idx][j+y]);
              dLmapv[dmi]  = DIST(lLv[ix][j], lLh[idx][j+y]);
              dCamapv[dmi] = DIST(lav[ix][j], lah[idx][j+y]);
              dCbmapv[dmi] = DIST(lbv[ix][j], lbh[idx][j+y]);
            }
            else if (sh<sv) {
              dLmaph[dmi]  = DIST(lLh[ix][j], lLv[idx][j+y]);
              dCamaph[dmi] = DIST(lah[ix][j], lav[idx][j+y]);
              dCbmaph[dmi] = DIST(lbh[ix][j], lbv[idx][j+y]);
              dLmapv[dmi]  = DIST(lLv[ix][j], lLv[idx][j+y]);
              dCamapv[dmi] = DIST(lav[ix][j], lav[idx][j+y]);
              dCbmapv[dmi] = DIST(lbv[ix][j], lbv[idx][j+y]);
            }
            else {
              dLmaph[dmi]  = DIST(lLh[ix][j], lLh[idx][j+y]);
              dCamaph[dmi] = DIST(lah[ix][j], lah[idx][j+y]);
              dCbmaph[dmi] = DIST(lbh[ix][j], lbh[idx][j+y]);
              dLmapv[dmi]  = DIST(lLv[ix][j], lLv[idx][j+y]);
              dCamapv[dmi] = DIST(lav[ix][j], lav[idx][j+y]);
              dCbmapv[dmi] = DIST(lbv[ix][j], lbv[idx][j+y]);
            }
          }
          else {
            dLmaph[dmi]  = DIST(lLh[ix][j], lLh[idx][j+y]);
            dCamaph[dmi] = DIST(lah[ix][j], lah[idx][j+y]);
            dCbmaph[dmi] = DIST(lbh[ix][j], lbh[idx][j+y]);
            dLmapv[dmi]  = DIST(lLv[ix][j], lLv[idx][j+y]);
            dCamapv[dmi] = DIST(lav[ix][j], lav[idx][j+y]);
            dCbmapv[dmi] = DIST(lbv[ix][j], lbv[idx][j+y]);
          }
          dmi++;
        }
      }
      // compute eL & eC
      int eL = MIN(MAX(dLmaph[3],dLmaph[5]),MAX(dLmapv[1],dLmapv[7]));
      int eCa = MIN(MAX(dCamaph[3],dCamaph[5]),MAX(dCamapv[1],dCamapv[7]));
      int eCb = MIN(MAX(dCbmaph[3],dCbmaph[5]),MAX(dCbmapv[1],dCbmapv[7]));

      int wh = 0;
      for (int dmi=0; dmi<9; dmi++) 
          if (dLmaph[dmi]<=eL && dCamaph[dmi]<=eCa && dCbmaph[dmi]<=eCb) 
               wh++;

      int wv = 0;
      for (int dmi=0; dmi<9; dmi++) 
          if (dLmapv[dmi]<=eL && dCamapv[dmi]<=eCa && dCbmapv[dmi]<=eCb) 
               wv++;
 
      homh[imx][j-1]+=wh;
      homh[imx][j]  +=wh;
      homh[imx][j+1]+=wh;
      homh[ix][j-1] +=wh;
      homh[ix][j]   +=wh;
      homh[ix][j+1] +=wh;
      homh[ipx][j-1]+=wh;
      homh[ipx][j]  +=wh;
      homh[ipx][j+1]+=wh;

      homv[imx][j-1]+=wv;
      homv[imx][j]  +=wv;
      homv[imx][j+1]+=wv;
      homv[ix][j-1] +=wv;
      homv[ix][j]   +=wv;
      homv[ix][j+1] +=wv;
      homv[ipx][j-1]+=wv;
      homv[ipx][j]  +=wv;
      homv[ipx][j+1]+=wv;
    }
//}
    // finalize image
    int hc, vc;
    for (int j=0; j<W; j++) {
      if (ISGREEN(ri,i-1,j))
        green[i-1][j] = ri->data[i-1][j];
      else { 
        hc = homh[imx][j];
        vc = homv[imx][j];
        if (hc > vc) 
          green[i-1][j] = gh[(i-1)%4][j];
        else if (hc < vc) 
          green[i-1][j] = gv[(i-1)%4][j];
        else 
          green[i-1][j] = (gh[(i-1)%4][j] + gv[(i-1)%4][j]) / 2;
      }
    }

    if (!(i%20) && plistener) 
      plistener->setProgress ((double)i / (H-2));
  }
  // finish H-2th and H-1th row, homogenity value is still valailable
  int hc, vc;
  for (int i=H-1; i<H+1; i++)
    for (int j=0; j<W; j++) {
      hc = homh[(i-1)%3][j];
      vc = homv[(i-1)%3][j];
      if (hc > vc)
        green[i-1][j] = gh[(i-1)%4][j];
      else if (hc < vc)
        green[i-1][j] = gv[(i-1)%4][j];
      else 
        green[i-1][j] = (gh[(i-1)%4][j] + gv[(i-1)%4][j]) / 2;
    }
    
    freeArray2<unsigned short>(rh, 3);
    freeArray2<unsigned short>(gh, 4);
    freeArray2<unsigned short>(bh, 3);
    freeArray2<unsigned short>(rv, 3);
    freeArray2<unsigned short>(gv, 4);
    freeArray2<unsigned short>(bv, 3);
    freeArray2<short>(lLh, 3);
    freeArray2<short>(lah, 3);
    freeArray2<short>(lbh, 3);
    freeArray2<unsigned short>(homh, 3);
    freeArray2<short>(lLv, 3);
    freeArray2<short>(lav, 3);
    freeArray2<short>(lbv, 3);
    freeArray2<unsigned short>(homv, 3);
}

void RawImageSource::hphd_vertical (float** hpmap, int col_from, int col_to) {

  float* temp = new float[MAX(W,H)];
  float* avg = new float[MAX(W,H)];
  float* dev = new float[MAX(W,H)];
  
  memset (temp, 0, MAX(W,H)*sizeof(float));
  memset (avg, 0, MAX(W,H)*sizeof(float));
  memset (dev, 0, MAX(W,H)*sizeof(float));
  
  for (int k=col_from; k<col_to; k++) {
    for (int i=5; i<H-5; i++) {
      temp[i] = (ri->data[i-5][k] - 8*ri->data[i-4][k] + 27*ri->data[i-3][k] - 48*ri->data[i-2][k] + 42*ri->data[i-1][k] - 
                (ri->data[i+5][k] - 8*ri->data[i+4][k] + 27*ri->data[i+3][k] - 48*ri->data[i+2][k] + 42*ri->data[i+1][k])) / 100.0;
      temp[i] = ABS(temp[i]);
    }
    for (int j=4; j<H-4; j++) {
        float avgL = (temp[j-4] + temp[j-3] + temp[j-2] + temp[j-1] + temp[j] + temp[j+1] + temp[j+2] + temp[j+3] + temp[j+4]) / 9.0;
        avg[j] = avgL;
        float devL = ((temp[j-4]-avgL)*(temp[j-4]-avgL) + (temp[j-3]-avgL)*(temp[j-3]-avgL) + (temp[j-2]-avgL)*(temp[j-2]-avgL) + (temp[j-1]-avgL)*(temp[j-1]-avgL) + (temp[j]-avgL)*(temp[j]-avgL) + (temp[j+1]-avgL)*(temp[j+1]-avgL) + (temp[j+2]-avgL)*(temp[j+2]-avgL) + (temp[j+3]-avgL)*(temp[j+3]-avgL) + (temp[j+4]-avgL)*(temp[j+4]-avgL)) / 9.0;
        if (devL<0.001) devL = 0.001;
        dev[j] = devL;
    }       
    for (int j=5; j<H-5; j++) {
        float avgL = avg[j-1];
        float avgR = avg[j+1];
        float devL = dev[j-1];
        float devR = dev[j+1];
        hpmap[j][k] = avgL + (avgR - avgL) * devL / (devL + devR);
    }
  }
  delete [] temp;
  delete [] avg;
  delete [] dev;
}

void RawImageSource::hphd_horizontal (float** hpmap, int row_from, int row_to) {

  float* temp = new float[MAX(W,H)];
  float* avg = new float[MAX(W,H)];
  float* dev = new float[MAX(W,H)];
  
  memset (temp, 0, MAX(W,H)*sizeof(float));
  memset (avg, 0, MAX(W,H)*sizeof(float));
  memset (dev, 0, MAX(W,H)*sizeof(float));

  for (int i=row_from; i<row_to; i++) {
    for (int j=5; j<W-5; j++) {
      temp[j] = (ri->data[i][j-5] - 8*ri->data[i][j-4] + 27*ri->data[i][j-3] - 48*ri->data[i][j-2] + 42*ri->data[i][j-1] - 
                (ri->data[i][j+5] - 8*ri->data[i][j+4] + 27*ri->data[i][j+3] - 48*ri->data[i][j+2] + 42*ri->data[i][j+1])) / 100;
      temp[j] = ABS(temp[j]);
    }
    for (int j=4; j<W-4; j++) {
        float avgL = (temp[j-4] + temp[j-3] + temp[j-2] + temp[j-1] + temp[j] + temp[j+1] + temp[j+2] + temp[j+3] + temp[j+4]) / 9.0;
        avg[j] = avgL;
        float devL = ((temp[j-4]-avgL)*(temp[j-4]-avgL) + (temp[j-3]-avgL)*(temp[j-3]-avgL) + (temp[j-2]-avgL)*(temp[j-2]-avgL) + (temp[j-1]-avgL)*(temp[j-1]-avgL) + (temp[j]-avgL)*(temp[j]-avgL) + (temp[j+1]-avgL)*(temp[j+1]-avgL) + (temp[j+2]-avgL)*(temp[j+2]-avgL) + (temp[j+3]-avgL)*(temp[j+3]-avgL) + (temp[j+4]-avgL)*(temp[j+4]-avgL)) / 9.0;
        if (devL<0.001) devL = 0.001;
        dev[j] = devL;
    }    
    for (int j=5; j<W-5; j++) {
        float avgL = avg[j-1];
        float avgR = avg[j+1];
        float devL = dev[j-1];
        float devR = dev[j+1];
        float hpv = avgL + (avgR - avgL) * devL / (devL + devR);
        if (hpmap[i][j] < 0.8*hpv) 
            this->hpmap[i][j] = 2;
        else if (hpv < 0.8*hpmap[i][j]) 
            this->hpmap[i][j] = 1;
        else
            this->hpmap[i][j] = 0;
    }        
  }
  delete [] temp;
  delete [] avg;
  delete [] dev;
}

void RawImageSource::hphd_green () {

  #pragma omp parallel for
  for (int i=3; i<H-3; i++) {
    for (int j=3; j<W-3; j++) {
      if (ISGREEN(ri,i,j))
        green[i][j] = ri->data[i][j];
      else {
        if (this->hpmap[i][j]==1) { 
            int g2 = ri->data[i][j+1] + ((ri->data[i][j] - ri->data[i][j+2]) >> 1);
            int g4 = ri->data[i][j-1] + ((ri->data[i][j] - ri->data[i][j-2]) >> 1);

            int dx = ri->data[i][j+1] - ri->data[i][j-1];
            int d1 = ri->data[i][j+3] - ri->data[i][j+1];
            int d2 = ri->data[i][j+2] - ri->data[i][j];
            int d3 = (ri->data[i-1][j+2] - ri->data[i-1][j]) >> 1;
            int d4 = (ri->data[i+1][j+2] - ri->data[i+1][j]) >> 1;
        
            double e2 = 1.0 / (1.0 + ABS(dx) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));
        
            d1 = ri->data[i][j-3] - ri->data[i][j-1];
            d2 = ri->data[i][j-2] - ri->data[i][j];
            d3 = (ri->data[i-1][j-2] - ri->data[i-1][j]) >> 1;
            d4 = (ri->data[i+1][j-2] - ri->data[i+1][j]) >> 1;
        
            double e4 = 1.0 / (1.0 + ABS(dx) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

            green[i][j] = CLIP((e2 * g2 + e4 * g4) / (e2 + e4));
        }
        else if (this->hpmap[i][j]==2) { 
            int g1 = ri->data[i-1][j] + ((ri->data[i][j] - ri->data[i-2][j]) >> 1);
            int g3 = ri->data[i+1][j] + ((ri->data[i][j] - ri->data[i+2][j]) >> 1);

            int dy = ri->data[i+1][j] - ri->data[i-1][j];
            int d1 = ri->data[i-1][j] - ri->data[i-3][j];
            int d2 = ri->data[i][j] - ri->data[i-2][j];
            int d3 = (ri->data[i][j-1] - ri->data[i-2][j-1]) >> 1;
            int d4 = (ri->data[i][j+1] - ri->data[i-2][j+1]) >> 1;
        
            double e1 = 1.0 / (1.0 + ABS(dy) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));
        
            d1 = ri->data[i+1][j] - ri->data[i+3][j];
            d2 = ri->data[i][j] - ri->data[i+2][j];
            d3 = (ri->data[i][j-1] - ri->data[i+2][j-1]) >> 1;
            d4 = (ri->data[i][j+1] - ri->data[i+2][j+1]) >> 1;
        
            double e3 = 1.0 / (1.0 + ABS(dy) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));
        
            green[i][j] = CLIP((e1 * g1 + e3 * g3) / (e1 + e3));
        }
        else {
            int g1 = ri->data[i-1][j] + ((ri->data[i][j] - ri->data[i-2][j]) >> 1);
            int g2 = ri->data[i][j+1] + ((ri->data[i][j] - ri->data[i][j+2]) >> 1);
            int g3 = ri->data[i+1][j] + ((ri->data[i][j] - ri->data[i+2][j]) >> 1);
            int g4 = ri->data[i][j-1] + ((ri->data[i][j] - ri->data[i][j-2]) >> 1);
        
            int dx = ri->data[i][j+1] - ri->data[i][j-1];
            int dy = ri->data[i+1][j] - ri->data[i-1][j];

            int d1 = ri->data[i-1][j] - ri->data[i-3][j];
            int d2 = ri->data[i][j] - ri->data[i-2][j];
            int d3 = (ri->data[i][j-1] - ri->data[i-2][j-1]) >> 1;
            int d4 = (ri->data[i][j+1] - ri->data[i-2][j+1]) >> 1;
        
            double e1 = 1.0 / (1.0 + ABS(dy) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

            d1 = ri->data[i][j+3] - ri->data[i][j+1];
            d2 = ri->data[i][j+2] - ri->data[i][j];
            d3 = (ri->data[i-1][j+2] - ri->data[i-1][j]) >> 1;
            d4 = (ri->data[i+1][j+2] - ri->data[i+1][j]) >> 1;
        
            double e2 = 1.0 / (1.0 + ABS(dx) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

            d1 = ri->data[i+1][j] - ri->data[i+3][j];
            d2 = ri->data[i][j] - ri->data[i+2][j];
            d3 = (ri->data[i][j-1] - ri->data[i+2][j-1]) >> 1;
            d4 = (ri->data[i][j+1] - ri->data[i+2][j+1]) >> 1;
        
            double e3 = 1.0 / (1.0 + ABS(dy) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));
            
            d1 = ri->data[i][j-3] - ri->data[i][j-1];
            d2 = ri->data[i][j-2] - ri->data[i][j];
            d3 = (ri->data[i-1][j-2] - ri->data[i-1][j]) >> 1;
            d4 = (ri->data[i+1][j-2] - ri->data[i+1][j]) >> 1;
        
            double e4 = 1.0 / (1.0 + ABS(dx) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));            

            green[i][j] = CLIP((e1*g1 + e2*g2 + e3*g3 + e4*g4) / (e1 + e2 + e3 + e4));
        }
      }
    }
  }
}

void RawImageSource::hphd_demosaic () {

  if (plistener) {
    plistener->setProgressStr ("Demosaicing...");
    plistener->setProgress (0.0);
  }

  float** hpmap = new float*[H];
  for (int i=0; i<H; i++) {
    hpmap[i] = new float[W];
    memset(hpmap[i], 0, W*sizeof(float));
  }
  
#ifdef _OPENMP
  #pragma omp parallel
  {
		int tid = omp_get_thread_num();
		int nthreads = omp_get_num_threads();
		int blk = W/nthreads;

		if (tid<nthreads-1)
			hphd_vertical (hpmap, tid*blk, (tid+1)*blk);
		else
			hphd_vertical (hpmap, tid*blk, W);
  }
#else
  hphd_vertical (hpmap, 0, W);
#endif
  if (plistener) 
    plistener->setProgress (0.33);

  this->hpmap = allocArray<char>(W, H);
  for (int i=0; i<H; i++)
    memset(this->hpmap[i], 0, W*sizeof(char));

#ifdef _OPENMP
  #pragma omp parallel
  {
		int tid = omp_get_thread_num();
		int nthreads = omp_get_num_threads();
		int blk = H/nthreads;

		if (tid<nthreads-1)
			hphd_horizontal (hpmap, tid*blk, (tid+1)*blk);
		else
			hphd_horizontal (hpmap, tid*blk, H);
  }
#else
  hphd_horizontal (hpmap, 0, H);
#endif
  freeArray<float>(hpmap, H);

  if (plistener) 
    plistener->setProgress (0.66);
 
  green = new unsigned short*[H];
  for (int i=0; i<H; i++)
    green[i] = new unsigned short[W];
    
  hphd_green ();

  if (plistener) 
    plistener->setProgress (1.0);
}

void RawImageSource::HLRecovery_Luminance (unsigned short* rin, unsigned short* gin, unsigned short* bin, unsigned short* rout, unsigned short* gout, unsigned short* bout, int width, int maxval) {

    for (int i=0; i<width; i++) {
        int r = rin[i], g = gin[i], b = bin[i];
		if (r>maxval || g>maxval || b>maxval) {
		    int ro = MIN (r, maxval);
		    int go = MIN (g, maxval);
		    int bo = MIN (b, maxval);
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
            int rr = L / 3.0 - H / 6.0 + C / 3.464101615;
            int gr = L / 3.0 - H / 6.0 - C / 3.464101615;
            int br = L / 3.0 + H / 3.0;
			rout[i] = CLIP(rr);
			gout[i] = CLIP(gr);
			bout[i] = CLIP(br);
		}
        else {
            rout[i] = rin[i];
            gout[i] = gin[i];
            bout[i] = bin[i];
        }
    }
}

void RawImageSource::HLRecovery_CIELab (unsigned short* rin, unsigned short* gin, unsigned short* bin, unsigned short* rout, unsigned short* gout, unsigned short* bout, int width, int maxval, double cam[3][3], double icam[3][3]) {

    static bool crTableReady = false;
    static double fv[0x10000];
    if (!crTableReady) {
    	for (int ix=0; ix < 0x10000; ix++) {
    	    double rx = ix / 65535.0;
        	fv[ix] = rx > 0.008856 ? exp(1.0/3 * log(rx)) : 7.787*rx + 16/116.0;
    	}
    	crTableReady = true;
    }

    for (int i=0; i<width; i++) {
        int r = rin[i], g = gin[i], b = bin[i];
		if (r>maxval || g>maxval || b>maxval) {
		    int ro = MIN (r, maxval);
		    int go = MIN (g, maxval);
		    int bo = MIN (b, maxval);
            double yy = cam[0][1]*r + cam[1][1]*g + cam[2][1]*b;
            double fy = fv[CLIP((int)yy)];
            // compute LCH decompostion of the clipped pixel (only color information, thus C and H will be used)
            double x = cam[0][0]*ro + cam[1][0]*go + cam[2][0]*bo;
            double y = cam[0][1]*ro + cam[1][1]*go + cam[2][1]*bo;
            double z = cam[0][2]*ro + cam[1][2]*go + cam[2][2]*bo;
            x = fv[CLIP((int)x)];
            y = fv[CLIP((int)y)];
            z = fv[CLIP((int)z)];
            // convert back to rgb
            double fz = fy - y + z;
            double fx = fy + x - y;
            double zr = (fz<=0.206893) ? ((116.0*fz-16.0)/903.3) : (fz * fz * fz);
            double xr = (fx<=0.206893) ? ((116.0*fx-16.0)/903.3) : (fx * fx * fx);
            x = xr*65535.0 - 0.5;
            y = yy;
            z = zr*65535.0 - 0.5;
            int rr = icam[0][0]*x + icam[1][0]*y + icam[2][0]*z;
            int gr = icam[0][1]*x + icam[1][1]*y + icam[2][1]*z;
            int br = icam[0][2]*x + icam[1][2]*y + icam[2][2]*z;
			rout[i] = CLIP(rr);
			gout[i] = CLIP(gr);
			bout[i] = CLIP(br);
		}
        else {
            rout[i] = rin[i];
            gout[i] = gin[i];
            bout[i] = bin[i];
        }
    }
}

void RawImageSource::hlRecovery (std::string method, unsigned short* red, unsigned short* green, unsigned short* blue, int i, int sx1, int width, int skip) {

    if (method=="Luminance")
        HLRecovery_Luminance (red, green, blue, red, green, blue, width, 65535 / ri->defgain);
    else if (method=="CIELab blending")
        HLRecovery_CIELab (red, green, blue, red, green, blue, width, 65535 / ri->defgain, cam, icam);
    else if (method=="Color")
        HLRecovery_ColorPropagation (red, green, blue, i, sx1, width, skip);
}

int RawImageSource::getAEHistogram (unsigned int* histogram, int& histcompr) {

    histcompr = 3;

    memset (histogram, 0, (65536>>histcompr)*sizeof(int));

    for (int i=border; i<ri->height-border; i++) {
        int start, end;
        if (fuji) {
            int fw = ri->fuji_width;
            start = ABS(fw-i) + border;
            end = MIN( ri->height+ ri->width-fw-i, fw+i) - border;
        }
        else {
            start = border;
            end = ri->width-border;
        }
        if (ri->filters)
            for (int j=start; j<end; j++)
                if (ISGREEN(ri,i,j))
                    histogram[ri->data[i][j]>>histcompr]+=2;
                else
                    histogram[ri->data[i][j]>>histcompr]+=4;
        else
            for (int j=start; j<3*end; j++) {
                    histogram[ri->data[i][j+0]>>histcompr]++;
                    histogram[ri->data[i][j+1]>>histcompr]++;
                    histogram[ri->data[i][j+2]>>histcompr]++;
            }
    }
    return 1;
}

ColorTemp RawImageSource::getAutoWB () {

    double avg_r = 0;
    double avg_g = 0;
    double avg_b = 0;
    int rn = 0, gn = 0, bn = 0;

    if (fuji) {
        for (int i=32; i<ri->height-32; i++) {
            int fw = ri->fuji_width;
            int start = ABS(fw-i) + 32;
            int end = MIN(ri->height+ri->width-fw-i, fw+i) - 32;
            for (int j=start; j<end; j++) {
                if (!ri->filters) {
                    double d = CLIP(ri->defgain*ri->data[i][3*j]);
                    if (d>64000)
                        continue;
                    avg_r += d*d*d*d*d*d; rn++;
                    d = CLIP(ri->defgain*ri->data[i][3*j+1]);
                    if (d>64000)
                        continue;
                    avg_g += d*d*d*d*d*d; gn++;
                    d = CLIP(ri->defgain*ri->data[i][3*j+2]);
                    if (d>64000)
                        continue;
                    avg_b += d*d*d*d*d*d; bn++;
                }
                else {
                    double d = CLIP(ri->defgain*ri->data[i][j]);
                    if (d>64000)
                        continue;
                    double dp = d*d*d*d*d*d;
                    if (ISRED(ri,i,j)) {
                        avg_r += dp;
                        rn++;
                    }
                    else if (ISGREEN(ri,i,j)) {
                        avg_g += dp;
                        gn++;
                    }
                    else if (ISBLUE(ri,i,j)) {
                        avg_b += dp;
                        bn++;
                    }
                }
            }
        }
    }
    else {
        for (int i=32; i<ri->height-32; i++)
            for (int j=32; j<ri->width-32; j++) {
                if (!ri->filters) {
                    double d = CLIP(ri->defgain*ri->data[i][3*j]);
                    if (d>64000)
                        continue;
                    avg_r += d*d*d*d*d*d; rn++;
                    d = CLIP(ri->defgain*ri->data[i][3*j+1]);
                    if (d>64000)
                        continue;
                    avg_g += d*d*d*d*d*d; gn++;
                    d = CLIP(ri->defgain*ri->data[i][3*j+2]);
                    if (d>64000)
                        continue;
                    avg_b += d*d*d*d*d*d; bn++;
                }
                else {
                    double d = CLIP(ri->defgain*ri->data[i][j]);
                    if (d>64000)
                        continue;
                    double dp = d*d*d*d*d*d;
                    if (ISRED(ri,i,j)) {
                        avg_r += dp;
                        rn++;
                    }
                    else if (ISGREEN(ri,i,j)) {
                        avg_g += dp;
                        gn++;
                    }
                    else if (ISBLUE(ri,i,j)) {
                        avg_b += dp;
                        bn++;
                    }
                }
            }
    }

    printf ("AVG: %g %g %g\n", avg_r/rn, avg_g/gn, avg_b/bn);

//    double img_r, img_g, img_b;
//    wb.getMultipliers (img_r, img_g, img_b);

//    return ColorTemp (pow(avg_r/rn, 1.0/6.0)*img_r, pow(avg_g/gn, 1.0/6.0)*img_g, pow(avg_b/bn, 1.0/6.0)*img_b);

    double reds   = pow (avg_r/rn, 1.0/6.0) * ri->camwb_red;
    double greens = pow (avg_g/gn, 1.0/6.0) * ri->camwb_green;
    double blues  = pow (avg_b/bn, 1.0/6.0) * ri->camwb_blue;
    
    double rm = coeff[0][0]*reds + coeff[0][1]*greens + coeff[0][2]*blues;
    double gm = coeff[1][0]*reds + coeff[1][1]*greens + coeff[1][2]*blues;
    double bm = coeff[2][0]*reds + coeff[2][1]*greens + coeff[2][2]*blues;

    return ColorTemp (rm, gm, bm);
}

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
        w = ri->fuji_width * 2 + 1;
        h = (H - ri->fuji_width)*2 + 1;
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
        tty = (ty-tx) / 2 + ri->fuji_width;
    }
    else {
        ttx = tx;
        tty = ty;
    }
}

ColorTemp RawImageSource::getSpotWB (std::vector<Coord2D> red, std::vector<Coord2D> green, std::vector<Coord2D>& blue, int tran) {

    int x; int y;
    int d[9][2] = {0,0, -1,-1, -1,0, -1,1, 0,-1, 0,1, 1,-1, 1,0, 1,1};
    double reds = 0, greens = 0, blues = 0;
    int rn = 0, gn = 0, bn = 0;
    
    if (!ri->filters) {
        for (int i=0; i<red.size(); i++) {
            transformPosition (red[i].x, red[i].y, tran, x, y);
            if (x>=0 && y>=0 && x<W && y<H) {
                reds += ri->data[y][3*x];
                rn++;
            }
            transformPosition (green[i].x, green[i].y, tran, x, y);
            if (x>=0 && y>=0 && x<W && y<H) {
                greens += ri->data[y][3*x+1];
                gn++;
            }
            transformPosition (blue[i].x, blue[i].y, tran, x, y);
            if (x>=0 && y>=0 && x<W && y<H) {
                blues += ri->data[y][3*x+2];
                bn++;
            }
        }
    }
    else {
        for (int i=0; i<red.size(); i++) {
            transformPosition (red[i].x, red[i].y, tran, x, y);
            for (int k=0; k<9; k++) {
                int xv = x + d[k][0];
                int yv = y + d[k][1];
                if (ISRED(ri,yv,xv) && xv>=0 && yv>=0 && xv<W && yv<H) {
                    reds += ri->data[yv][xv];
                    rn++;
                    break;
                }
            }
            transformPosition (green[i].x, green[i].y, tran, x, y);
            for (int k=0; k<9; k++) {
                int xv = x + d[k][0];
                int yv = y + d[k][1];
                if (ISGREEN(ri,yv,xv) && xv>=0 && yv>=0 && xv<W && yv<H) {
                    greens += ri->data[yv][xv];
                    gn++;
                    break;
                }
            }
            transformPosition (blue[i].x, blue[i].y, tran, x, y);
            for (int k=0; k<9; k++) {
                int xv = x + d[k][0];
                int yv = y + d[k][1];
                if (ISBLUE(ri,yv,xv) && xv>=0 && yv>=0 && xv<W && yv<H) {
                    blues += ri->data[yv][xv];
                    bn++;
                    break;
                }
            }
        }
    }

    reds = reds/rn * ri->camwb_red;
    greens = greens/gn * ri->camwb_green;
    blues = blues/bn * ri->camwb_blue;

    double rm = coeff[0][0]*reds + coeff[0][1]*greens + coeff[0][2]*blues;
    double gm = coeff[1][0]*reds + coeff[1][1]*greens + coeff[1][2]*blues;
    double bm = coeff[2][0]*reds + coeff[2][1]*greens + coeff[2][2]*blues;

    return ColorTemp (rm, gm, bm);
}

#define FORCC for (c=0; c < colors; c++)
#define fc(row,col) \
	(ri->prefilters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)
typedef unsigned short ushort;
void RawImageSource::vng4_demosaic () {

  static const signed char *cp, terms[] = {
    -2,-2,+0,-1,0,0x01, -2,-2,+0,+0,1,0x01, -2,-1,-1,+0,0,0x01,
    -2,-1,+0,-1,0,0x02, -2,-1,+0,+0,0,0x03, -2,-1,+0,+1,1,0x01,
    -2,+0,+0,-1,0,0x06, -2,+0,+0,+0,1,0x02, -2,+0,+0,+1,0,0x03,
    -2,+1,-1,+0,0,0x04, -2,+1,+0,-1,1,0x04, -2,+1,+0,+0,0,0x06,
    -2,+1,+0,+1,0,0x02, -2,+2,+0,+0,1,0x04, -2,+2,+0,+1,0,0x04,
    -1,-2,-1,+0,0,0x80, -1,-2,+0,-1,0,0x01, -1,-2,+1,-1,0,0x01,
    -1,-2,+1,+0,1,0x01, -1,-1,-1,+1,0,0x88, -1,-1,+1,-2,0,0x40,
    -1,-1,+1,-1,0,0x22, -1,-1,+1,+0,0,0x33, -1,-1,+1,+1,1,0x11,
    -1,+0,-1,+2,0,0x08, -1,+0,+0,-1,0,0x44, -1,+0,+0,+1,0,0x11,
    -1,+0,+1,-2,1,0x40, -1,+0,+1,-1,0,0x66, -1,+0,+1,+0,1,0x22,
    -1,+0,+1,+1,0,0x33, -1,+0,+1,+2,1,0x10, -1,+1,+1,-1,1,0x44,
    -1,+1,+1,+0,0,0x66, -1,+1,+1,+1,0,0x22, -1,+1,+1,+2,0,0x10,
    -1,+2,+0,+1,0,0x04, -1,+2,+1,+0,1,0x04, -1,+2,+1,+1,0,0x04,
    +0,-2,+0,+0,1,0x80, +0,-1,+0,+1,1,0x88, +0,-1,+1,-2,0,0x40,
    +0,-1,+1,+0,0,0x11, +0,-1,+2,-2,0,0x40, +0,-1,+2,-1,0,0x20,
    +0,-1,+2,+0,0,0x30, +0,-1,+2,+1,1,0x10, +0,+0,+0,+2,1,0x08,
    +0,+0,+2,-2,1,0x40, +0,+0,+2,-1,0,0x60, +0,+0,+2,+0,1,0x20,
    +0,+0,+2,+1,0,0x30, +0,+0,+2,+2,1,0x10, +0,+1,+1,+0,0,0x44,
    +0,+1,+1,+2,0,0x10, +0,+1,+2,-1,1,0x40, +0,+1,+2,+0,0,0x60,
    +0,+1,+2,+1,0,0x20, +0,+1,+2,+2,0,0x10, +1,-2,+1,+0,0,0x80,
    +1,-1,+1,+1,0,0x88, +1,+0,+1,+2,0,0x08, +1,+0,+2,-1,0,0x40,
    +1,+0,+2,+1,0,0x10
  }, chood[] = { -1,-1, -1,0, -1,+1, 0,+1, +1,+1, +1,0, +1,-1, 0,-1 };

  if (plistener) {
    plistener->setProgressStr ("Demosaicing...");
    plistener->setProgress (0.0);
  }

  ushort (*brow[5])[4], *pix;
  int prow=7, pcol=1, *ip, *code[16][16], gval[8], gmin, gmax, sum[4];
  int row, col, x, y, x1, x2, y1, y2, t, weight, grads, color, diag;
  int g, diff, thold, num, c, width=W, height=H, colors=4;
  ushort (*image)[4];
  int lcode[16][16][32], shift, i, j;

  image = (ushort (*)[4]) calloc (H*W, sizeof *image);
  for (int ii=0; ii<H; ii++)
    for (int jj=0; jj<W; jj++)
        image[ii*W+jj][fc(ii,jj)] = ri->data[ii][jj];

// first linear interpolation
  for (row=0; row < 16; row++)
    for (col=0; col < 16; col++) {
      ip = lcode[row][col];
      memset (sum, 0, sizeof sum);
      for (y=-1; y <= 1; y++)
	for (x=-1; x <= 1; x++) {
	  shift = (y==0) + (x==0);
	  if (shift == 2) continue;
	  color = fc(row+y,col+x);
	  *ip++ = (width*y + x)*4 + color;
	  *ip++ = shift;
	  *ip++ = color;
	  sum[color] += 1 << shift;
	}
      FORCC
	if (c != fc(row,col)) {
	  *ip++ = c;
	  *ip++ = 256 / sum[c];
	}
    }

  for (row=1; row < height-1; row++)
    for (col=1; col < width-1; col++) {
      pix = image[row*width+col];
      ip = lcode[row & 15][col & 15];
      memset (sum, 0, sizeof sum);
      for (i=8; i--; ip+=3)
	sum[ip[2]] += pix[ip[0]] << ip[1];
      for (i=colors; --i; ip+=2)
	pix[ip[0]] = sum[ip[0]] * ip[1] >> 8;
    }

//  lin_interpolate();


  ip = (int *) calloc ((prow+1)*(pcol+1), 1280);
  for (row=0; row <= prow; row++)		/* Precalculate for VNG */
    for (col=0; col <= pcol; col++) {
      code[row][col] = ip;
      for (cp=terms, t=0; t < 64; t++) {
	y1 = *cp++;  x1 = *cp++;
	y2 = *cp++;  x2 = *cp++;
	weight = *cp++;
	grads = *cp++;
	color = fc(row+y1,col+x1);
	if (fc(row+y2,col+x2) != color) continue;
	diag = (fc(row,col+1) == color && fc(row+1,col) == color) ? 2:1;
	if (abs(y1-y2) == diag && abs(x1-x2) == diag) continue;
	*ip++ = (y1*width + x1)*4 + color;
	*ip++ = (y2*width + x2)*4 + color;
	*ip++ = weight;
	for (g=0; g < 8; g++)
	  if (grads & 1<<g) *ip++ = g;
	*ip++ = -1;
      }
      *ip++ = INT_MAX;
      for (cp=chood, g=0; g < 8; g++) {
	y = *cp++;  x = *cp++;
	*ip++ = (y*width + x) * 4;
	color = fc(row,col);
	if (fc(row+y,col+x) != color && fc(row+y*2,col+x*2) == color)
	  *ip++ = (y*width + x) * 8 + color;
	else
	  *ip++ = 0;
      }
    }
  brow[4] = (ushort (*)[4]) calloc (width*3, sizeof **brow);
  for (row=0; row < 3; row++)
    brow[row] = brow[4] + row*width;
  for (row=2; row < height-2; row++) {		/* Do VNG interpolation */
    for (col=2; col < width-2; col++) {
      color = fc(row,col);
      pix = image[row*width+col];
      ip = code[row & prow][col & pcol];
      memset (gval, 0, sizeof gval);
      while ((g = ip[0]) != INT_MAX) {		/* Calculate gradients */
	diff = ABS(pix[g] - pix[ip[1]]) << ip[2];
	gval[ip[3]] += diff;
	ip += 5;
	if ((g = ip[-1]) == -1) continue;
	gval[g] += diff;
	while ((g = *ip++) != -1)
	  gval[g] += diff;
      }
      ip++;
      gmin = gmax = gval[0];			/* Choose a threshold */
      for (g=1; g < 8; g++) {
	if (gmin > gval[g]) gmin = gval[g];
	if (gmax < gval[g]) gmax = gval[g];
      }
      if (gmax == 0) {
	memcpy (brow[2][col], pix, sizeof *image);
	continue;
      }
      thold = gmin + (gmax >> 1);
      memset (sum, 0, sizeof sum);
      for (num=g=0; g < 8; g++,ip+=2) {		/* Average the neighbors */
	if (gval[g] <= thold) {
	  FORCC
	    if (c == color && ip[1])
	      sum[c] += (pix[c] + pix[ip[1]]) >> 1;
	    else
	      sum[c] += pix[ip[0] + c];
	  num++;
	}
      }
      FORCC {					/* Save to buffer */
	t = pix[color];
	if (c != color)
	  t += (sum[c] - sum[color]) / num;
	brow[2][col][c] = CLIP(t);
      }
    }
    if (row > 3)				/* Write buffer to image */
      memcpy (image[(row-2)*width+2], brow[0]+2, (width-4)*sizeof *image);
    for (g=0; g < 4; g++)
      brow[(g-1) & 3] = brow[g];
    if (!(row%20) && plistener) 
      plistener->setProgress ((double)row / (H-2));
  }
  memcpy (image[(row-2)*width+2], brow[0]+2, (width-4)*sizeof *image);
  memcpy (image[(row-1)*width+2], brow[1]+2, (width-4)*sizeof *image);
  free (brow[4]);
  free (code[0][0]);

  green = new unsigned short*[H];
  for (int i=0; i<H; i++) {
    green[i] = new unsigned short[W];
    for (int j=0; j<W; j++)
        green[i][j] = (image[i*W+j][1] + image[i*W+j][3]) >> 1;
  }
  free (image);
}

//#define ABS(x) (((int)(x) ^ ((int)(x) >> 31)) - ((int)(x) >> 31))
//#define MIN(a,b) ((a) < (b) ? (a) : (b))
//#define MAX(a,b) ((a) > (b) ? (a) : (b))
//#define LIM(x,min,max) MAX(min,MIN(x,max))
//#define CLIP(x) LIM(x,0,65535)
#undef fc
#define fc(row,col) \
	(ri->filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)
#define FC(x,y) fc(x,y)
#define LIM(x,min,max) MAX(min,MIN(x,max))
#define ULIM(x,y,z) ((y) < (z) ? LIM(x,y,z) : LIM(x,z,y))

/*
   Patterned Pixel Grouping Interpolation by Alain Desbiolles
*/
void RawImageSource::ppg_demosaic()
{
  int width=W, height=H;
  int dir[5] = { 1, width, -1, -width, 1 };
  int row, col, diff[2], guess[2], c, d, i;
  ushort (*pix)[4];

  ushort (*image)[4];
  int colors = 3;

  if (plistener) {
    plistener->setProgressStr ("Demosaicing...");
    plistener->setProgress (0.0);
  }
  
  image = (ushort (*)[4]) calloc (H*W, sizeof *image);
  for (int ii=0; ii<H; ii++)
    for (int jj=0; jj<W; jj++)
        image[ii*W+jj][fc(ii,jj)] = ri->data[ii][jj];

  border_interpolate(3, image);

/*  Fill in the green layer with gradients and pattern recognition: */
  for (row=3; row < height-3; row++) {
    for (col=3+(FC(row,3) & 1), c=FC(row,col); col < width-3; col+=2) {
      pix = image + row*width+col;
      for (i=0; (d=dir[i]) > 0; i++) {
	guess[i] = (pix[-d][1] + pix[0][c] + pix[d][1]) * 2
		      - pix[-2*d][c] - pix[2*d][c];
	diff[i] = ( ABS(pix[-2*d][c] - pix[ 0][c]) +
		    ABS(pix[ 2*d][c] - pix[ 0][c]) +
		    ABS(pix[  -d][1] - pix[ d][1]) ) * 3 +
		  ( ABS(pix[ 3*d][1] - pix[ d][1]) +
		    ABS(pix[-3*d][1] - pix[-d][1]) ) * 2;
      }
      d = dir[i = diff[0] > diff[1]];
      pix[0][1] = ULIM(guess[i] >> 2, pix[d][1], pix[-d][1]);
    }
    if(plistener) plistener->setProgress(0.33*row/(height-3));
  }
/*  Calculate red and blue for each green pixel:		*/
  for (row=1; row < height-1; row++) {
    for (col=1+(FC(row,2) & 1), c=FC(row,col+1); col < width-1; col+=2) {
      pix = image + row*width+col;
      for (i=0; (d=dir[i]) > 0; c=2-c, i++)
	pix[0][c] = CLIP((pix[-d][c] + pix[d][c] + 2*pix[0][1]
			- pix[-d][1] - pix[d][1]) >> 1);
    }
    if(plistener) plistener->setProgress(0.33 + 0.33*row/(height-1));
  }
/*  Calculate blue for red pixels and vice versa:		*/
  for (row=1; row < height-1; row++) {
    for (col=1+(FC(row,1) & 1), c=2-FC(row,col); col < width-1; col+=2) {
      pix = image + row*width+col;
      for (i=0; (d=dir[i]+dir[i+1]) > 0; i++) {
	diff[i] = ABS(pix[-d][c] - pix[d][c]) +
		  ABS(pix[-d][1] - pix[0][1]) +
		  ABS(pix[ d][1] - pix[0][1]);
	guess[i] = pix[-d][c] + pix[d][c] + 2*pix[0][1]
		 - pix[-d][1] - pix[d][1];
      }
      if (diff[0] != diff[1])
	pix[0][c] = CLIP(guess[diff[0] > diff[1]] >> 1);
      else
	pix[0][c] = CLIP((guess[0]+guess[1]) >> 2);
    }
    if(plistener) plistener->setProgress(0.67 + 0.33*row/(height-1));
  }

  red = new unsigned short*[H];
  for (int i=0; i<H; i++) {
    red[i] = new unsigned short[W];
    for (int j=0; j<W; j++)
        red[i][j] = image[i*W+j][0];
  }
  green = new unsigned short*[H];
  for (int i=0; i<H; i++) {
    green[i] = new unsigned short[W];
    for (int j=0; j<W; j++)
        green[i][j] = image[i*W+j][1];
  }
  blue = new unsigned short*[H];
  for (int i=0; i<H; i++) {
    blue[i] = new unsigned short[W];
    for (int j=0; j<W; j++)
        blue[i][j] = image[i*W+j][2];
  }
  free (image);
}

void RawImageSource::border_interpolate(int border, ushort (*image)[4])
{
  unsigned row, col, y, x, f, c, sum[8];
  int width=W, height=H;
  int colors = 3;

  for (row=0; row < height; row++)
    for (col=0; col < width; col++) {
      if (col==border && row >= border && row < height-border)
	col = width-border;
      memset (sum, 0, sizeof sum);
      for (y=row-1; y != row+2; y++)
	for (x=col-1; x != col+2; x++)
	  if (y < height && x < width) {
	    f = fc(y,x);
	    sum[f] += image[y*width+x][f];
	    sum[f+4]++;
	  }
      f = fc(row,col);
      FORCC if (c != f && sum[c+4])
	image[row*width+col][c] = sum[c] / sum[c+4];
    }
}

/*
   Adaptive Homogeneity-Directed interpolation is based on
   the work of Keigo Hirakawa, Thomas Parks, and Paul Lee.
 */
#define TS 256		/* Tile Size */
#define FORC(cnt) for (c=0; c < cnt; c++)
#define FORC3 FORC(3)
#define SQR(x) ((x)*(x))

void RawImageSource::ahd_demosaic()
{
  int i, j, k, top, left, row, col, tr, tc, c, d, val, hm[2];
  ushort (*pix)[4], (*rix)[3];
  static const int dir[4] = { -1, 1, -TS, TS };
  unsigned ldiff[2][4], abdiff[2][4], leps, abeps;
  float r, cbrt[0x10000], xyz[3], xyz_cam[3][4];
  ushort (*rgb)[TS][TS][3];
   short (*lab)[TS][TS][3], (*lix)[3];
   char (*homo)[TS][TS], *buffer;

  int width=W, height=H;
  ushort (*image)[4];
  int colors = 3;

  const double xyz_rgb[3][3] = {			/* XYZ from RGB */
    { 0.412453, 0.357580, 0.180423 },
    { 0.212671, 0.715160, 0.072169 },
    { 0.019334, 0.119193, 0.950227 } };
  const float d65_white[3] = { 0.950456, 1, 1.088754 };

  if (plistener) {
    plistener->setProgressStr ("Demosaicing...");
    plistener->setProgress (0.0);
  }
  
  image = (ushort (*)[4]) calloc (H*W, sizeof *image);
  for (int ii=0; ii<H; ii++)
    for (int jj=0; jj<W; jj++)
        image[ii*W+jj][fc(ii,jj)] = ri->data[ii][jj];

  for (i=0; i < 0x10000; i++) {
    r = i / 65535.0;
    cbrt[i] = r > 0.008856 ? pow(r,1/3.0) : 7.787*r + 16/116.0;
  }
  for (i=0; i < 3; i++)
    for (j=0; j < colors; j++)
      for (xyz_cam[i][j] = k=0; k < 3; k++)
	xyz_cam[i][j] += xyz_rgb[i][k] * coeff[k][j] / d65_white[i];

  border_interpolate(5, image);
  buffer = (char *) malloc (26*TS*TS);		/* 1664 kB */
  //merror (buffer, "ahd_interpolate()");
  rgb  = (ushort(*)[TS][TS][3]) buffer;
  lab  = (short (*)[TS][TS][3])(buffer + 12*TS*TS);
  homo = (char  (*)[TS][TS])   (buffer + 24*TS*TS);

  for (top=2; top < height-5; top += TS-6)
    for (left=2; left < width-5; left += TS-6) {

/*  Interpolate green horizontally and vertically:		*/
      for (row = top; row < top+TS && row < height-2; row++) {
	col = left + (FC(row,left) & 1);
	for (c = FC(row,col); col < left+TS && col < width-2; col+=2) {
	  pix = image + row*width+col;
	  val = ((pix[-1][1] + pix[0][c] + pix[1][1]) * 2
		- pix[-2][c] - pix[2][c]) >> 2;
	  rgb[0][row-top][col-left][1] = ULIM(val,pix[-1][1],pix[1][1]);
	  val = ((pix[-width][1] + pix[0][c] + pix[width][1]) * 2
		- pix[-2*width][c] - pix[2*width][c]) >> 2;
	  rgb[1][row-top][col-left][1] = ULIM(val,pix[-width][1],pix[width][1]);
	}
      }

    if(plistener) plistener->setProgress (0.33);
/*  Interpolate red and blue, and convert to CIELab:		*/
      for (d=0; d < 2; d++)
	for (row=top+1; row < top+TS-1 && row < height-3; row++)
	  for (col=left+1; col < left+TS-1 && col < width-3; col++) {
	    pix = image + row*width+col;
	    rix = &rgb[d][row-top][col-left];
	    lix = &lab[d][row-top][col-left];
	    if ((c = 2 - FC(row,col)) == 1) {
	      c = FC(row+1,col);
	      val = pix[0][1] + (( pix[-1][2-c] + pix[1][2-c]
				 - rix[-1][1] - rix[1][1] ) >> 1);
	      rix[0][2-c] = CLIP(val);
	      val = pix[0][1] + (( pix[-width][c] + pix[width][c]
				 - rix[-TS][1] - rix[TS][1] ) >> 1);
	    } else
	      val = rix[0][1] + (( pix[-width-1][c] + pix[-width+1][c]
				 + pix[+width-1][c] + pix[+width+1][c]
				 - rix[-TS-1][1] - rix[-TS+1][1]
				 - rix[+TS-1][1] - rix[+TS+1][1] + 1) >> 2);
	    rix[0][c] = CLIP(val);
	    c = FC(row,col);
	    rix[0][c] = pix[0][c];
	    xyz[0] = xyz[1] = xyz[2] = 0.5;
	    FORCC {
	      xyz[0] += xyz_cam[0][c] * rix[0][c];
	      xyz[1] += xyz_cam[1][c] * rix[0][c];
	      xyz[2] += xyz_cam[2][c] * rix[0][c];
	    }
	    xyz[0] = cbrt[CLIP((int) xyz[0])];
	    xyz[1] = cbrt[CLIP((int) xyz[1])];
	    xyz[2] = cbrt[CLIP((int) xyz[2])];
	    lix[0][0] = 64 * (116 * xyz[1] - 16);
	    lix[0][1] = 64 * 500 * (xyz[0] - xyz[1]);
	    lix[0][2] = 64 * 200 * (xyz[1] - xyz[2]);
	  }

      if(plistener) plistener->setProgress (0.5);
/*  Build homogeneity maps from the CIELab images:		*/
      memset (homo, 0, 2*TS*TS);
      for (row=top+2; row < top+TS-2 && row < height-4; row++) {
	tr = row-top;
	for (col=left+2; col < left+TS-2 && col < width-4; col++) {
	  tc = col-left;
	  for (d=0; d < 2; d++) {
	    lix = &lab[d][tr][tc];
	    for (i=0; i < 4; i++) {
	       ldiff[d][i] = ABS(lix[0][0]-lix[dir[i]][0]);
	      abdiff[d][i] = SQR(lix[0][1]-lix[dir[i]][1])
			   + SQR(lix[0][2]-lix[dir[i]][2]);
	    }
	  }
	  leps = MIN(MAX(ldiff[0][0],ldiff[0][1]),
		     MAX(ldiff[1][2],ldiff[1][3]));
	  abeps = MIN(MAX(abdiff[0][0],abdiff[0][1]),
		      MAX(abdiff[1][2],abdiff[1][3]));
	  for (d=0; d < 2; d++)
	    for (i=0; i < 4; i++)
	      if (ldiff[d][i] <= leps && abdiff[d][i] <= abeps)
		homo[d][tr][tc]++;
	}
      }
      if(plistener) plistener->setProgress (0.8);
/*  Combine the most homogenous pixels for the final result:	*/
      for (row=top+3; row < top+TS-3 && row < height-5; row++) {
	tr = row-top;
	for (col=left+3; col < left+TS-3 && col < width-5; col++) {
	  tc = col-left;
	  for (d=0; d < 2; d++)
	    for (hm[d]=0, i=tr-1; i <= tr+1; i++)
	      for (j=tc-1; j <= tc+1; j++)
		hm[d] += homo[d][i][j];
	  if (hm[0] != hm[1])
	    FORC3 image[row*width+col][c] = rgb[hm[1] > hm[0]][tr][tc][c];
	  else
	    FORC3 image[row*width+col][c] =
		(rgb[0][tr][tc][c] + rgb[1][tr][tc][c]) >> 1;
	}
      }
    }
  if(plistener) plistener->setProgress (1.0);
  free (buffer);
  red = new unsigned short*[H];
  for (int i=0; i<H; i++) {
    red[i] = new unsigned short[W];
    for (int j=0; j<W; j++)
        red[i][j] = image[i*W+j][0];
  }
  green = new unsigned short*[H];
  for (int i=0; i<H; i++) {
    green[i] = new unsigned short[W];
    for (int j=0; j<W; j++)
        green[i][j] = image[i*W+j][1];
  }
  blue = new unsigned short*[H];
  for (int i=0; i<H; i++) {
    blue[i] = new unsigned short[W];
    for (int j=0; j<W; j++)
        blue[i][j] = image[i*W+j][2];
  }
  free (image);
}
#undef TS

/*
 *      Redistribution and use in source and binary forms, with or without
 *      modification, are permitted provided that the following conditions are
 *      met:
 *      
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above
 *        copyright notice, this list of conditions and the following disclaimer
 *        in the documentation and/or other materials provided with the
 *        distribution.
 *      * Neither the name of the author nor the names of its
 *        contributors may be used to endorse or promote products derived from
 *        this software without specific prior written permission.
 *      
 *      THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *      "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *      LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *      A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *      OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *      SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *      LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *      DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *      THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *      (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *      OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// If you want to use the code, you need to display name of the original authors in
// your software!

/* DCB demosaicing by Jacek Gozdz (cuniek@kft.umcs.lublin.pl)
 * the code is open source (BSD licence)
*/

#define TILESIZE 256
#define TILEBORDER 10
#define CACHESIZE (TILESIZE+2*TILEBORDER)

inline void RawImageSource::dcb_initTileLimits(int &colMin, int &rowMin, int &colMax, int &rowMax, int x0, int y0, int border)
{
	rowMin = border;
	colMin = border;
	rowMax = CACHESIZE-border;
	colMax = CACHESIZE-border;
	if(!y0 ) rowMin = TILEBORDER+border;
	if(!x0 ) colMin = TILEBORDER+border;
	if( y0+TILESIZE > H-border) rowMax = TILEBORDER+H-border-y0;
	if( x0+TILESIZE > W-border) colMax = TILEBORDER+W-border-x0;
}

void RawImageSource::fill_raw( ushort (*cache )[4], int x0, int y0, ushort** rawData)
{
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,0);

    for (int row=rowMin,y=y0-TILEBORDER+rowMin; row<rowMax; row++,y++)
    	for (int col=colMin,x=x0-TILEBORDER+colMin; col<colMax; col++,x++){
    		cache[row*CACHESIZE+col][fc(y,x)] = rawData[y][x];
    	}
}

void RawImageSource::fill_border( ushort (*cache )[4], int border, int x0, int y0)
{
	unsigned row, col, y, x, f, c, sum[8];
	int colors = 3;

	for (row = y0; row < y0+TILESIZE+TILEBORDER && row<H; row++){
		for (col = x0; col < x0+TILESIZE+TILEBORDER && col<W; col++) {
			if (col >= border && col < W - border && row >= border && row < H - border){
				col = W - border;
				if(col >= x0+TILESIZE+TILEBORDER )
					break;
			}
			memset(sum, 0, sizeof sum);
			for (y = row - 1; y != row + 2; y++)
				for (x = col - 1; x != col + 2; x++)
					if (y < H && y< y0+TILESIZE+TILEBORDER && x < W && x<x0+TILESIZE+TILEBORDER) {
						f = fc(y,x);
						sum[f] += cache[(y-y0 +TILEBORDER)* CACHESIZE +TILEBORDER+ x-x0][f];
						sum[f + 4]++;
					}
			f = fc(row,col);
			FORCC
				if (c != f && sum[c + 4])
					cache[(row-y0+TILEBORDER) * CACHESIZE +TILEBORDER + col-x0][c] = sum[c] / sum[c + 4];
		}
	}
}
// saves red and blue
void RawImageSource::copy_to_buffer( ushort (*buffer)[2], ushort (*image)[4])
{
	for (int indx=0; indx < CACHESIZE*CACHESIZE; indx++) {
		buffer[indx][0]=image[indx][0]; //R
		buffer[indx][1]=image[indx][2]; //B
	}
}

// restores red and blue
void RawImageSource::restore_from_buffer(ushort (*image)[4], ushort (*buffer)[2])
{
	for (int indx=0; indx < CACHESIZE*CACHESIZE; indx++) {
		image[indx][0]=buffer[indx][0]; //R
		image[indx][2]=buffer[indx][1]; //B
	}
}

// fast green interpolation
void RawImageSource::hid(ushort (*image)[4], int x0, int y0)
{
	int u=CACHESIZE, v=2*CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,2);
	
	for (int row = rowMin; row < rowMax; row++) {
		for (int col=colMin, indx=row*CACHESIZE+col; col < colMax; col++, indx++) {
			int c = fc(y0-TILEBORDER+row,x0-TILEBORDER+col);
			if( c != 1)
			{
				int current = (image[indx+u][1] + image[indx-u][1] + image[indx-1][1] + image[indx+1][1])/4 +
						      (image[indx][c] - ( image[indx+v][c] + image[indx-v][c] + image[indx-2][c] + image[indx+2][c])/4)/2;
				image[indx][1] = CLIP(current);
			}
		}
	}	
	
}

// missing colors are interpolated
void RawImageSource::dcb_color(ushort (*image)[4], int x0, int y0)
{
	int u=CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,1);

	// red in blue pixel, blue in red pixel
	for (int row=rowMin; row < rowMax; row++)
		for (int col=colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin) & 1), indx=row*CACHESIZE+col, c=2-FC(y0-TILEBORDER+row,x0-TILEBORDER+col); col < colMax; col+=2, indx+=2) {
			int current = ( 4*image[indx][1]
			                - image[indx+u+1][1] - image[indx+u-1][1] - image[indx-u+1][1] - image[indx-u-1][1]
			                + image[indx+u+1][c] + image[indx+u-1][c] + image[indx-u+1][c] + image[indx-u-1][c] )/4;
			image[indx][c] = CLIP(current);
		}

	// red or blue in green pixels
	for (int row=rowMin; row<rowMax; row++)
		for (int col=colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin+1)&1), indx=row*CACHESIZE+col,c=FC(y0-TILEBORDER+row,x0-TILEBORDER+col+1),d=2-c; col<colMax; col+=2, indx+=2) {
			int current = (2*image[indx][1] - image[indx+1][1] - image[indx-1][1] + image[indx+1][c] + image[indx-1][c])/2;
			image[indx][c] = CLIP( current );
			current = (2*image[indx][1] - image[indx+u][1] - image[indx-u][1] + image[indx+u][d] + image[indx-u][d])/2;
			image[indx][d] = CLIP( current );
		}	
}

// green correction
void RawImageSource::hid2(ushort (*image)[4], int x0, int y0)
{
	int u=CACHESIZE, v=2*CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,4);
	
	for (int row=rowMin; row < rowMax; row++) {
		for (int col=colMin, indx=row*CACHESIZE+col; col < colMax; col++, indx++) {
			int c = fc(y0-TILEBORDER+row,x0-TILEBORDER+col);
			if ( c != 1)
			{
				int current = (image[indx+v][1] + image[indx-v][1] + image[indx-2][1] + image[indx+2][1])/4 +
						   image[indx][c] - ( image[indx+v][c] + image[indx-v][c] + image[indx-2][c] + image[indx+2][c])/4;
				image[indx][1]=CLIP(current);
			}
		}
	}	
}

// green is used to create
// an interpolation direction map 
// 1 = vertical
// 0 = horizontal
// saved in image[][3]
void RawImageSource::dcb_map(ushort (*image)[4], int x0, int y0)
{	
	int u=CACHESIZE, v=2*CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,2);

	for (int row=rowMin; row < rowMax; row++) {
		for (int col=colMin, indx=row*CACHESIZE+col; col < colMax; col++, indx++) {
			if (image[indx][1] > ( image[indx-1][1] + image[indx+1][1] + image[indx-u][1] + image[indx+u][1])/4 )
				image[indx][3] = ((MIN( image[indx-1][1], image[indx+1][1]) + image[indx-1][1] + image[indx+1][1] ) < (MIN( image[indx-u][1], image[indx+u][1]) + image[indx-u][1] + image[indx+u][1]));
			else
				image[indx][3] = ((MAX( image[indx-1][1], image[indx+1][1]) + image[indx-1][1] + image[indx+1][1] ) > (MAX( image[indx-u][1], image[indx+u][1]) + image[indx-u][1] + image[indx+u][1]));
		}
	}
}


// interpolated green pixels are corrected using the map
void RawImageSource::dcb_correction(ushort (*image)[4], int x0, int y0)
{
	int u=CACHESIZE, v=2*CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,4);
	
	for (int row=rowMin; row < rowMax; row++) {
		for (int col=colMin, indx=row*CACHESIZE+col; col < colMax; col++, indx++) {
			if ( fc(y0-TILEBORDER+row,x0-TILEBORDER+col) != 1)
			{
				int current = 4*image[indx][3] +
						  2*(image[indx+u][3] + image[indx-u][3] + image[indx+1][3] + image[indx-1][3]) +
							image[indx+v][3] + image[indx-v][3] + image[indx+2][3] + image[indx-2][3];

				image[indx][1] = ((16-current)*(image[indx-1][1] + image[indx+1][1])/2 + current*(image[indx-u][1] + image[indx+u][1])/2)/16;
			}
		}
	}
}

// R and B smoothing using green contrast, all pixels except 2 pixel wide border
void RawImageSource::dcb_pp(ushort (*image)[4], int x0, int y0)
{
	int u=CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,2);
	
	for (int row=rowMin; row < rowMax; row++)
		for (int col=colMin, indx=row*CACHESIZE+col; col < colMax; col++, indx++) {
			int r1 = ( image[indx-1][0] + image[indx+1][0] + image[indx-u][0] + image[indx+u][0] + image[indx-u-1][0] + image[indx+u+1][0] + image[indx-u+1][0] + image[indx+u-1][0])/8;
			int g1 = ( image[indx-1][1] + image[indx+1][1] + image[indx-u][1] + image[indx+u][1] + image[indx-u-1][1] + image[indx+u+1][1] + image[indx-u+1][1] + image[indx+u-1][1])/8;
			int b1 = ( image[indx-1][2] + image[indx+1][2] + image[indx-u][2] + image[indx+u][2] + image[indx-u-1][2] + image[indx+u+1][2] + image[indx-u+1][2] + image[indx+u-1][2])/8;

			r1 = r1 + ( image[indx][1] - g1 );
			b1 = b1 + ( image[indx][1] - g1 );
			image[indx][0] = CLIP(r1);
			image[indx][2] = CLIP(b1);
		}
}

// interpolated green pixels are corrected using the map
// with correction
void RawImageSource::dcb_correction2(ushort (*image)[4], int x0, int y0)
{
	int u=CACHESIZE, v=2*CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,4);
	
	for (int row=rowMin; row < rowMax; row++) {
		for (int col=colMin, indx=row*CACHESIZE+col; col < colMax; col++, indx++) {
			int c =  fc(y0-TILEBORDER+row,x0-TILEBORDER+col);
			if (c != 1)
			{
				int current = 4*image[indx][3] +
						  2*(image[indx+u][3] + image[indx-u][3] + image[indx+1][3] + image[indx-1][3]) +
							image[indx+v][3] + image[indx-v][3] + image[indx+2][3] + image[indx-2][3];
				current = ((16-current)*((image[indx-1][1] + image[indx+1][1])/2 + image[indx][c] - (image[indx+2][c] + image[indx-2][c])/2) + current*((image[indx-u][1] + image[indx+u][1])/2 + image[indx][c] - (image[indx+v][c] + image[indx-v][c])/2))/16;
				image[indx][1] = CLIP(current);
			}
		}
	}
}

// image refinement
void RawImageSource::dcb_refinement(ushort (*image)[4], int x0, int y0)
{
	int u=CACHESIZE, v=2*CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,5);

	float f[4];
	int g[4];
	
	for (int row=rowMin; row < rowMax; row++)
		for (int col=colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin)&1),indx=row*CACHESIZE+col,c=FC(y0-TILEBORDER+row,x0-TILEBORDER+col); col < colMax; col+=2,indx+=2){

// Cubic Spline Interpolation by Li and Randhawa, modified by Jacek Gozdz and Luis Sanz Rodríguez
		f[0]=1.0/(1.0+abs(image[indx-u][c]-image[indx][c])+abs(image[indx-u][1]-image[indx][1]));
		f[1]=1.0/(1.0+abs(image[indx+1][c]-image[indx][c])+abs(image[indx+1][1]-image[indx][1]));
		f[2]=1.0/(1.0+abs(image[indx-1][c]-image[indx][c])+abs(image[indx-1][1]-image[indx][1]));
		f[3]=1.0/(1.0+abs(image[indx+u][c]-image[indx][c])+abs(image[indx+u][1]-image[indx][1]));

		g[0]=CLIP(image[indx-u][1]+(image[indx][c]-image[indx-u][c])/2 + (image[indx][c]-image[indx-v][c])/4);
		g[1]=CLIP(image[indx+1][1]+(image[indx][c]-image[indx+1][c])/2 + (image[indx][c]-image[indx+2][c])/4);
		g[2]=CLIP(image[indx-1][1]+(image[indx][c]-image[indx-1][c])/2 + (image[indx][c]-image[indx-2][c])/4);
		g[3]=CLIP(image[indx+u][1]+(image[indx][c]-image[indx+u][c])/2 + (image[indx][c]-image[indx+v][c])/4);

		image[indx][1]=CLIP(((f[0]*g[0]+f[1]*g[1]+f[2]*g[2]+f[3]*g[3])/(f[0]+f[1]+f[2]+f[3]) ));

// get rid of the overshooted pixels
		int min = MIN(image[indx+1+u][1], MIN(image[indx+1-u][1], MIN(image[indx-1+u][1], MIN(image[indx-1-u][1], MIN(image[indx-1][1], MIN(image[indx+1][1], MIN(image[indx-u][1], image[indx+u][1])))))));
		int max = MAX(image[indx+1+u][1], MAX(image[indx+1-u][1], MAX(image[indx-1+u][1], MAX(image[indx-1-u][1], MAX(image[indx-1][1], MAX(image[indx+1][1], MAX(image[indx-u][1], image[indx+u][1])))))));

		image[indx][1] =  LIM(image[indx][1], min, max);

			
		}
}

// missing colors are interpolated using high quality algorithm by Luis Sanz Rodríguez
void RawImageSource::dcb_color_full(ushort (*image)[4], int x0, int y0, float (*chroma)[2])
{
	int u=CACHESIZE, v=2*CACHESIZE, w=3*CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,3);

	int i,j;
	float f[4],g[4];

	for (int row=1; row < CACHESIZE-1; row++)
		for (int col=1+(FC(y0-TILEBORDER+row,x0-TILEBORDER+1)&1),indx=row*CACHESIZE+col,c=FC(y0-TILEBORDER+row,x0-TILEBORDER+col),d=c/2; col < CACHESIZE-1; col+=2,indx+=2)
			chroma[indx][d]=image[indx][c]-image[indx][1];

	for (int row=rowMin; row<rowMax; row++)
		for (int col=colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin)&1),indx=row*CACHESIZE+col,c=1-FC(y0-TILEBORDER+row,x0-TILEBORDER+col)/2,d=1-c; col<colMax; col+=2,indx+=2) {
			f[0]=1.0/(float)(1.0+fabs(chroma[indx-u-1][c]-chroma[indx+u+1][c])+fabs(chroma[indx-u-1][c]-chroma[indx-w-3][c])+fabs(chroma[indx+u+1][c]-chroma[indx-w-3][c]));
			f[1]=1.0/(float)(1.0+fabs(chroma[indx-u+1][c]-chroma[indx+u-1][c])+fabs(chroma[indx-u+1][c]-chroma[indx-w+3][c])+fabs(chroma[indx+u-1][c]-chroma[indx-w+3][c]));
			f[2]=1.0/(float)(1.0+fabs(chroma[indx+u-1][c]-chroma[indx-u+1][c])+fabs(chroma[indx+u-1][c]-chroma[indx+w+3][c])+fabs(chroma[indx-u+1][c]-chroma[indx+w-3][c]));
			f[3]=1.0/(float)(1.0+fabs(chroma[indx+u+1][c]-chroma[indx-u-1][c])+fabs(chroma[indx+u+1][c]-chroma[indx+w-3][c])+fabs(chroma[indx-u-1][c]-chroma[indx+w+3][c]));
			g[0]=1.325*chroma[indx-u-1][c]-0.175*chroma[indx-w-3][c]-0.075*chroma[indx-w-1][c]-0.075*chroma[indx-u-3][c];
			g[1]=1.325*chroma[indx-u+1][c]-0.175*chroma[indx-w+3][c]-0.075*chroma[indx-w+1][c]-0.075*chroma[indx-u+3][c];
			g[2]=1.325*chroma[indx+u-1][c]-0.175*chroma[indx+w-3][c]-0.075*chroma[indx+w-1][c]-0.075*chroma[indx+u-3][c];
			g[3]=1.325*chroma[indx+u+1][c]-0.175*chroma[indx+w+3][c]-0.075*chroma[indx+w+1][c]-0.075*chroma[indx+u+3][c];
			chroma[indx][c]=(f[0]*g[0]+f[1]*g[1]+f[2]*g[2]+f[3]*g[3])/(f[0]+f[1]+f[2]+f[3]);
		}
	for (int row=rowMin; row<rowMax; row++)
		for (int col=colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin+1)&1),indx=row*CACHESIZE+col,c=FC(y0-TILEBORDER+row,x0-TILEBORDER+col+1)/2; col<colMax; col+=2,indx+=2)
			for(int d=0;d<=1;c=1-c,d++){
				f[0]=1.0/(float)(1.0+fabs(chroma[indx-u][c]-chroma[indx+u][c])+fabs(chroma[indx-u][c]-chroma[indx-w][c])+fabs(chroma[indx+u][c]-chroma[indx-w][c]));
				f[1]=1.0/(float)(1.0+fabs(chroma[indx+1][c]-chroma[indx-1][c])+fabs(chroma[indx+1][c]-chroma[indx+3][c])+fabs(chroma[indx-1][c]-chroma[indx+3][c]));
				f[2]=1.0/(float)(1.0+fabs(chroma[indx-1][c]-chroma[indx+1][c])+fabs(chroma[indx-1][c]-chroma[indx-3][c])+fabs(chroma[indx+1][c]-chroma[indx-3][c]));
				f[3]=1.0/(float)(1.0+fabs(chroma[indx+u][c]-chroma[indx-u][c])+fabs(chroma[indx+u][c]-chroma[indx+w][c])+fabs(chroma[indx-u][c]-chroma[indx+w][c]));
			
				g[0]=0.875*chroma[indx-u][c]+0.125*chroma[indx-w][c];
				g[1]=0.875*chroma[indx+1][c]+0.125*chroma[indx+3][c];
				g[2]=0.875*chroma[indx-1][c]+0.125*chroma[indx-3][c];
				g[3]=0.875*chroma[indx+u][c]+0.125*chroma[indx+w][c];				

				chroma[indx][c]=(f[0]*g[0]+f[1]*g[1]+f[2]*g[2]+f[3]*g[3])/(f[0]+f[1]+f[2]+f[3]);
			}

	for(int row=rowMin; row<rowMax; row++)
		for(int col=colMin,indx=row*CACHESIZE+col; col<colMax; col++,indx++){
			image[indx][0]=CLIP(chroma[indx][0]+image[indx][1]);
			image[indx][2]=CLIP(chroma[indx][1]+image[indx][1]);
		}
}

// DCB demosaicing main routine (sharp version)
void RawImageSource::dcb_demosaic(int iterations, int dcb_enhance)
{

    red = new unsigned short*[H];
    green = new unsigned short*[H];
    blue = new unsigned short*[H];
    for (int i=0; i<H; i++) {
        red[i] = new unsigned short[W];
        green[i] = new unsigned short[W];
        blue[i] = new unsigned short[W];
        memset(red[i],0,W*2);
        memset(green[i],0,W*2);
        memset(blue[i],0,W*2);
    }
    if(plistener) {
        plistener->setProgressStr ("DCB Demosaicing...");
        plistener->setProgress (0.0);
    }

    int wTiles = W/TILESIZE + (W%TILESIZE?1:0);
    int hTiles = H/TILESIZE + (H%TILESIZE?1:0);
    int numTiles = wTiles * hTiles;

#ifdef _OPENMP
	int nthreads = omp_get_max_threads();
 	ushort (**image)[4]  =  (ushort(**)[4]) calloc( nthreads,sizeof( void*) );
	ushort (**image2)[2] =	(ushort(**)[2]) calloc( nthreads,sizeof( void*) );
	float  (**chroma)[2] =  (float (**)[2]) calloc( nthreads,sizeof( void*) );
	for(int i=0; i<nthreads; i++){
		image[i] = (ushort(*)[4]) calloc( CACHESIZE*CACHESIZE, sizeof **image);
		image2[i]= (ushort(*)[2]) calloc( CACHESIZE*CACHESIZE, sizeof **image2);
		chroma[i]= (float (*)[2]) calloc( CACHESIZE*CACHESIZE, sizeof **chroma);
	}
#else
	ushort (*image)[4]  = (ushort(*)[4]) calloc( CACHESIZE*CACHESIZE, sizeof *image);
	ushort (*image2)[2] = (ushort(*)[2]) calloc( CACHESIZE*CACHESIZE, sizeof *image2);
	float  (*chroma)[2] = (float (*)[2]) calloc( CACHESIZE*CACHESIZE, sizeof *chroma);
#endif

#pragma omp parallel for
    for( int iTile=0; iTile < numTiles; iTile++){
    	int xTile = iTile % wTiles;
    	int yTile = iTile / wTiles;
    	int x0 = xTile*TILESIZE;
    	int y0 = yTile*TILESIZE;
#ifdef _OPENMP
    	int tid = omp_get_thread_num();
    	ushort (*tile)[4]   = image[tid];
    	ushort (*buffer)[2] = image2[tid];
    	float  (*chrm)[2]   = chroma[tid];
#else
    	ushort (*tile)[4]   = image;
    	ushort (*buffer)[2] = image2;
    	float  (*chrm)[2]   = chroma;
#endif

		fill_raw( tile, x0,y0,ri->data );
		if( !xTile || !yTile || xTile==wTiles-1 || yTile==hTiles-1)
		   fill_border(tile,2, x0, y0);
        copy_to_buffer(buffer, tile);
        hid(tile,x0,y0);

       	dcb_color(tile,x0,y0);
        for (int i=iterations; i>0;i--) {
			hid2(tile,x0,y0);
			hid2(tile,x0,y0);
			hid2(tile,x0,y0);
			dcb_map(tile,x0,y0);
			dcb_correction(tile,x0,y0);
        }
        if( iterations )
        	dcb_color(tile,x0,y0);
        dcb_pp(tile,x0,y0);
        hid2(tile,x0,y0);
        hid2(tile,x0,y0);
        hid2(tile,x0,y0);
        dcb_map(tile,x0,y0);
        dcb_correction2(tile,x0,y0);
        restore_from_buffer(tile, buffer);
        dcb_map(tile,x0,y0);
        dcb_correction(tile,x0,y0);
        dcb_color(tile,x0,y0);
        dcb_pp(tile,x0,y0);
        dcb_map(tile,x0,y0);
        dcb_correction(tile,x0,y0);
        dcb_map(tile,x0,y0);
        dcb_correction(tile,x0,y0);
        restore_from_buffer(tile, buffer);
        dcb_color(tile,x0,y0);
        if (dcb_enhance) {
			dcb_refinement(tile,x0,y0);
			dcb_color_full(tile,x0,y0,chrm);
        }

        for(int y=0;y<TILESIZE && y0+y<H;y++){
			for (int j=0; j<TILESIZE && x0+j<W; j++){
				red[y0+y][x0+j]   = tile[(y+TILEBORDER)*CACHESIZE+TILEBORDER+j][0];
				green[y0+y][x0+j] = tile[(y+TILEBORDER)*CACHESIZE+TILEBORDER+j][1];
				blue[y0+y][x0+j]  = tile[(y+TILEBORDER)*CACHESIZE+TILEBORDER+j][2];
			}
        }
    }

#ifdef _OPENMP
	for(int i=0; i<nthreads; i++){
		free(image[i]);
		free(image2[i]);
		free(chroma[i])
	}
#endif
	free(image);
    free(image2);
    free(chroma);

    if(plistener) plistener->setProgress (1.0);

}
#undef TILEBORDER
#undef TILESIZE
	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Emil's code for AMaZE
#include "amaze_interpolate_RT.cc"//AMaZE demosaic	
#include "CA_correct_RT.cc"//Emil's CA auto correction
#include "cfa_linedn_RT.cc"//Emil's CA auto correction
#include "green_equil_RT.cc"//Emil's green channel equilibration
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


} /* namespace */

