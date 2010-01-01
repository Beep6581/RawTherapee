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

namespace rtengine {

int loadRaw (const char* fname, struct RawImage* ri);

extern const Settings* settings;

#undef ABS
#undef MAX
#undef MIN
#undef DIST
#undef MAXVAL
#undef CLIP
#undef THREAD_PRIORITY_NORMAL

#define ABS(a) ((a)<0?-(a):(a))
#define MAX(a,b) ((a)<(b)?(b):(a))
#define MIN(a,b) ((a)>(b)?(b):(a))
#define DIST(a,b) (ABS(a-b))
#define MAXVAL  0xffff
#define CLIP(a) ((a)>0?((a)<MAXVAL?(a):MAXVAL):0)

RawImageSource::RawImageSource () : ImageSource(), plistener(NULL), green(NULL), cache(NULL), border(4) {

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
        if (ri->filters) {
            if (i==0)
                interpolate_row_rb_mul_pp (red, blue, NULL, green[i], green[i+1], i, rm, gm, bm, sx1, imwidth, pp.skip);
            else if (i==H-1)
                interpolate_row_rb_mul_pp (red, blue, green[i-1], green[i], NULL, i, rm, gm, bm, sx1, imwidth, pp.skip);
            else 
                interpolate_row_rb_mul_pp (red, blue, green[i-1], green[i], green[i+1], i, rm, gm, bm, sx1, imwidth, pp.skip);
                   
            for (int j=0,jx=sx1; j<imwidth; j++,jx+=pp.skip) 
                grn[j] = CLIP(gm*green[i][jx]);
        }
        else {
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
    if (((idata->getMake().size()>=7 && idata->getMake().substr(0,7)=="OLYMPUS" && idata->getModel()[0]=='E') || (idata->getMake().size()>=9 && idata->getMake().substr(0,7)=="Panasonic")) && settings->demosaicMethod!="vng4" && ri->filters) {
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

        // local correction in a 9x9 box
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
        delete corr_data;
*/                        
    }
    
    if (ri->filters) {
        // demosaic
        if (settings->demosaicMethod=="hphd")
            hphd_demosaic ();
        else if (settings->demosaicMethod=="vng4")
            vng4_demosaic ();
        else
            eahd_demosaic ();
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

    MyTime t1, t2;

    t1.set ();
    for (int t=0; t<times; t++) {
        if (settings->dualThreadEnabled) {
            Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &RawImageSource::correction_YIQ_LQ_), im, 1, im->height/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &RawImageSource::correction_YIQ_LQ_), im, im->height/2, im->height-1), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            thread1->join ();
            thread2->join ();
        }
        else
            correction_YIQ_LQ_ (im, 1, im->height-1);
    }
    t2.set ();
//    printf ("Corrected. (%d)\n", t2.etime(t1));
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

void RawImageSource::hphd_green (int row_from, int row_to) {

  for (int i=row_from; i<row_to; i++) {
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
  
//  MyTime t1, t2, t3, t4;
//  t1.set ();

  // vertical 
  if (settings->dualThreadEnabled) {
      Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &RawImageSource::hphd_vertical), hpmap, 0, W/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
      Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &RawImageSource::hphd_vertical), hpmap, W/2, W), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
      thread1->join ();
      thread2->join ();
  }
  else
      hphd_vertical (hpmap, 0, W);

  if (plistener) 
    plistener->setProgress (0.33);

  // horizontal
  this->hpmap = allocArray<char>(W, H);
  for (int i=0; i<H; i++)
    memset(this->hpmap[i], 0, W*sizeof(char));
  
  if (settings->dualThreadEnabled) {
      Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &RawImageSource::hphd_horizontal), hpmap, 0, H/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
      Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &RawImageSource::hphd_horizontal), hpmap, H/2, H), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
      thread1->join ();
      thread2->join ();
  }
  else
      hphd_horizontal (hpmap, 0, H);

  freeArray<float>(hpmap, H);

  if (plistener) 
    plistener->setProgress (0.66);
 
// reconstruct G
  green = new unsigned short*[H];
  for (int i=0; i<H; i++)
    green[i] = new unsigned short[W];
    
  if (settings->dualThreadEnabled) {
      Glib::Thread *thread1 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &RawImageSource::hphd_green), 3, H/2), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
      Glib::Thread *thread2 = Glib::Thread::create(sigc::bind(sigc::mem_fun(*this, &RawImageSource::hphd_green), H/2, H-3), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
      thread1->join ();
      thread2->join ();
  }
  else
      hphd_green (3, H-3);

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

int RawImageSource::getAEHistogram (int* histogram, int& histcompr) {

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
}

