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
#include <rawimagesource_i.h>
#include <median.h>
#include <rawimage.h>
#include <math.h>
#include <mytime.h>
#include <iccmatrices.h>
#include <iccstore.h>
#include <image8.h>
#include <curves.h>
#include <dfmanager.h>
#include <slicer.h>
#include <assert.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace rtengine {
#undef ABS
#undef MAX
#undef MIN
#undef DIST

#define ABS(a) ((a)<0?-(a):(a))
#define MAX(a,b) ((a)<(b)?(b):(a))
#define MIN(a,b) ((a)>(b)?(b):(a))
#define DIST(a,b) (ABS(a-b))
#define CLIREF(x) LIM(x,-200000.0,200000.0) // avoid overflow : do not act directly on image[] or pix[]

#define PIX_SORT(a,b) { if ((a)>(b)) {temp=(a);(a)=(b);(b)=temp;} }
extern Settings* settings;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::eahd_demosaic () {
  if (plistener) {
    plistener->setProgressStr ("Demosaicing...");
    plistener->setProgress (0.0);
  }

  // prepare cache and constants for cielab conversion
	//TODO: revisit after conversion to D50 illuminant
	lc00 = (0.412453 * rgb_cam[0][0] + 0.357580 * rgb_cam[0][1] + 0.180423 * rgb_cam[0][2]) ;// / 0.950456;
  lc01 = (0.412453 * rgb_cam[1][0] + 0.357580 * rgb_cam[1][1] + 0.180423 * rgb_cam[1][2]) ;// / 0.950456;
  lc02 = (0.412453 * rgb_cam[2][0] + 0.357580 * rgb_cam[2][1] + 0.180423 * rgb_cam[2][2]) ;// / 0.950456;

  lc10 = 0.212671 * rgb_cam[0][0] + 0.715160 * rgb_cam[0][1] + 0.072169 * rgb_cam[0][2];
  lc11 = 0.212671 * rgb_cam[1][0] + 0.715160 * rgb_cam[1][1] + 0.072169 * rgb_cam[1][2];
  lc12 = 0.212671 * rgb_cam[2][0] + 0.715160 * rgb_cam[2][1] + 0.072169 * rgb_cam[2][2];

  lc20 = (0.019334 * rgb_cam[0][0] + 0.119193 * rgb_cam[0][1] + 0.950227 * rgb_cam[0][2]) ;// / 1.088754;
  lc21 = (0.019334 * rgb_cam[1][0] + 0.119193 * rgb_cam[1][1] + 0.950227 * rgb_cam[1][2]) ;// / 1.088754;
  lc22 = (0.019334 * rgb_cam[2][0] + 0.119193 * rgb_cam[2][1] + 0.950227 * rgb_cam[2][2]) ;// / 1.088754;

  int maxindex = 2*65536;
  cache = new double[maxindex];
  threshold = (int)(0.008856*MAXVAL);
  for (int i=0; i<maxindex; i++)
    cache[i] = exp(1.0/3.0 * log((double)i / MAXVAL));

  // end of cielab preparation

  float* rh[3];
  float* gh[4];
  float* bh[3];
  float* rv[3];
  float* gv[4];
  float* bv[3];
  float* lLh[3];
  float* lah[3];
  float* lbh[3];
  float* lLv[3];
  float* lav[3];
  float* lbv[3];
  float* homh[3];
  float* homv[3];

  for (int i=0; i<4; i++) {
    gh[i] = new float[W];
    gv[i] = new float[W];
  }

  for (int i=0; i<3; i++) {
    rh[i] = new float[W];
    bh[i] = new float[W];
    rv[i] = new float[W];
    bv[i] = new float[W];
    lLh[i] = new float[W];
    lah[i] = new float[W];
    lbh[i] = new float[W];
    lLv[i] = new float[W];
    lav[i] = new float[W];
    lbv[i] = new float[W];
    homh[i] = new float[W];
    homv[i] = new float[W];
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
      if (ri->ISGREEN(i-1,j))
        green[i-1][j] = rawData[i-1][j];
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

    freeArray2<float>(rh, 3);
    freeArray2<float>(gh, 4);
    freeArray2<float>(bh, 3);
    freeArray2<float>(rv, 3);
    freeArray2<float>(gv, 4);
    freeArray2<float>(bv, 3);
    freeArray2<float>(lLh, 3);
    freeArray2<float>(lah, 3);
    freeArray2<float>(lbh, 3);
    freeArray2<float>(homh, 3);
    freeArray2<float>(lLv, 3);
    freeArray2<float>(lav, 3);
    freeArray2<float>(lbv, 3);
    freeArray2<float>(homv, 3);

    // Interpolate R and B
    for (int i=0; i<H; i++) {
  	  if (i==0)
  		  // rm, gm, bm must be recovered
  		  //interpolate_row_rb_mul_pp (red, blue, NULL, green[i], green[i+1], i, rm, gm, bm, 0, W, 1);
  		  interpolate_row_rb_mul_pp (red[i], blue[i], NULL, green[i], green[i+1], i, 1.0, 1.0, 1.0, 0, W, 1);
  	  else if (i==H-1)
  		  interpolate_row_rb_mul_pp (red[i], blue[i], green[i-1], green[i], NULL, i, 1.0, 1.0, 1.0, 0, W, 1);
  	  else
  		  interpolate_row_rb_mul_pp (red[i], blue[i], green[i-1], green[i], green[i+1], i, 1.0, 1.0, 1.0, 0, W, 1);
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::hphd_vertical (float** hpmap, int col_from, int col_to) {
  float* temp = new float[MAX(W,H)];
  float* avg = new float[MAX(W,H)];
  float* dev = new float[MAX(W,H)];

  memset (temp, 0, MAX(W,H)*sizeof(float));
  memset (avg, 0, MAX(W,H)*sizeof(float));
  memset (dev, 0, MAX(W,H)*sizeof(float));

  for (int k=col_from; k<col_to; k++) {
    for (int i=5; i<H-5; i++) {
      temp[i] = (rawData[i-5][k] - 8*rawData[i-4][k] + 27*rawData[i-3][k] - 48*rawData[i-2][k] + 42*rawData[i-1][k] -
                (rawData[i+5][k] - 8*rawData[i+4][k] + 27*rawData[i+3][k] - 48*rawData[i+2][k] + 42*rawData[i+1][k])) / 100.0;
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

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::hphd_horizontal (float** hpmap, int row_from, int row_to) {
  float* temp = new float[MAX(W,H)];
  float* avg = new float[MAX(W,H)];
  float* dev = new float[MAX(W,H)];

  memset (temp, 0, MAX(W,H)*sizeof(float));
  memset (avg, 0, MAX(W,H)*sizeof(float));
  memset (dev, 0, MAX(W,H)*sizeof(float));

  for (int i=row_from; i<row_to; i++) {
    for (int j=5; j<W-5; j++) {
      temp[j] = (rawData[i][j-5] - 8*rawData[i][j-4] + 27*rawData[i][j-3] - 48*rawData[i][j-2] + 42*rawData[i][j-1] -
                (rawData[i][j+5] - 8*rawData[i][j+4] + 27*rawData[i][j+3] - 48*rawData[i][j+2] + 42*rawData[i][j+1])) / 100;
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
            hpmap[i][j] = 2;
        else if (hpv < 0.8*hpmap[i][j])
            hpmap[i][j] = 1;
        else
            hpmap[i][j] = 0;
    }
  }
  delete [] temp;
  delete [] avg;
  delete [] dev;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::hphd_green (float** hpmap) {
  #pragma omp parallel for
  for (int i=3; i<H-3; i++) {
    for (int j=3; j<W-3; j++) {
      if (ri->ISGREEN(i,j))
        green[i][j] = rawData[i][j];
      else {
        if (hpmap[i][j]==1) {
            int g2 = rawData[i][j+1] + ((rawData[i][j] - rawData[i][j+2]) /2);
            int g4 = rawData[i][j-1] + ((rawData[i][j] - rawData[i][j-2]) /2);

            int dx = rawData[i][j+1] - rawData[i][j-1];
            int d1 = rawData[i][j+3] - rawData[i][j+1];
            int d2 = rawData[i][j+2] - rawData[i][j];
            int d3 = (rawData[i-1][j+2] - rawData[i-1][j]) /2;
            int d4 = (rawData[i+1][j+2] - rawData[i+1][j]) /2;

            double e2 = 1.0 / (1.0 + ABS(dx) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

            d1 = rawData[i][j-3] - rawData[i][j-1];
            d2 = rawData[i][j-2] - rawData[i][j];
            d3 = (rawData[i-1][j-2] - rawData[i-1][j]) /2;
            d4 = (rawData[i+1][j-2] - rawData[i+1][j]) /2;

            double e4 = 1.0 / (1.0 + ABS(dx) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

            green[i][j] = (e2 * g2 + e4 * g4) / (e2 + e4);
        }
        else if (hpmap[i][j]==2) {
            int g1 = rawData[i-1][j] + ((rawData[i][j] - rawData[i-2][j]) /2);
            int g3 = rawData[i+1][j] + ((rawData[i][j] - rawData[i+2][j]) /2);

            int dy = rawData[i+1][j] - rawData[i-1][j];
            int d1 = rawData[i-1][j] - rawData[i-3][j];
            int d2 = rawData[i][j] - rawData[i-2][j];
            int d3 = (rawData[i][j-1] - rawData[i-2][j-1]) /2;
            int d4 = (rawData[i][j+1] - rawData[i-2][j+1]) /2;

            double e1 = 1.0 / (1.0 + ABS(dy) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

            d1 = rawData[i+1][j] - rawData[i+3][j];
            d2 = rawData[i][j] - rawData[i+2][j];
            d3 = (rawData[i][j-1] - rawData[i+2][j-1]) /2;
            d4 = (rawData[i][j+1] - rawData[i+2][j+1]) /2;

            double e3 = 1.0 / (1.0 + ABS(dy) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

            green[i][j] = (e1 * g1 + e3 * g3) / (e1 + e3);
        }
        else {
            int g1 = rawData[i-1][j] + ((rawData[i][j] - rawData[i-2][j]) /2);
            int g2 = rawData[i][j+1] + ((rawData[i][j] - rawData[i][j+2]) /2);
            int g3 = rawData[i+1][j] + ((rawData[i][j] - rawData[i+2][j]) /2);
            int g4 = rawData[i][j-1] + ((rawData[i][j] - rawData[i][j-2]) /2);

            int dx = rawData[i][j+1] - rawData[i][j-1];
            int dy = rawData[i+1][j] - rawData[i-1][j];

            int d1 = rawData[i-1][j] - rawData[i-3][j];
            int d2 = rawData[i][j] - rawData[i-2][j];
            int d3 = (rawData[i][j-1] - rawData[i-2][j-1]) /2;
            int d4 = (rawData[i][j+1] - rawData[i-2][j+1]) /2;

            double e1 = 1.0 / (1.0 + ABS(dy) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

            d1 = rawData[i][j+3] - rawData[i][j+1];
            d2 = rawData[i][j+2] - rawData[i][j];
            d3 = (rawData[i-1][j+2] - rawData[i-1][j]) /2;
            d4 = (rawData[i+1][j+2] - rawData[i+1][j]) /2;

            double e2 = 1.0 / (1.0 + ABS(dx) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

            d1 = rawData[i+1][j] - rawData[i+3][j];
            d2 = rawData[i][j] - rawData[i+2][j];
            d3 = (rawData[i][j-1] - rawData[i+2][j-1]) /2;
            d4 = (rawData[i][j+1] - rawData[i+2][j+1]) /2;

            double e3 = 1.0 / (1.0 + ABS(dy) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

            d1 = rawData[i][j-3] - rawData[i][j-1];
            d2 = rawData[i][j-2] - rawData[i][j];
            d3 = (rawData[i-1][j-2] - rawData[i-1][j]) /2;
            d4 = (rawData[i+1][j-2] - rawData[i+1][j]) /2;

            double e4 = 1.0 / (1.0 + ABS(dx) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

            green[i][j] = (e1*g1 + e2*g2 + e3*g3 + e4*g4) / (e1 + e2 + e3 + e4);
        }
      }
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::hphd_demosaic () {
  if (plistener) {
    plistener->setProgressStr ("Demosaicing...");
    plistener->setProgress (0.0);
  }

  float** hpmap = allocArray< float >(W,H, true);

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

  for (int i=0; i<H; i++)
    memset(hpmap[i], 0, W*sizeof(char));

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

  hphd_green (hpmap);
  freeArray<float>(hpmap, H);

  if (plistener)
    plistener->setProgress (0.66);

  for (int i=0; i<H; i++) {
	  if (i==0)
		  // rm, gm, bm must be recovered
		  //interpolate_row_rb_mul_pp (red, blue, NULL, green[i], green[i+1], i, rm, gm, bm, 0, W, 1);
		  interpolate_row_rb_mul_pp (red[i], blue[i], NULL, green[i], green[i+1], i, 1.0, 1.0, 1.0, 0, W, 1);
	  else if (i==H-1)
		  interpolate_row_rb_mul_pp (red[i], blue[i], green[i-1], green[i], NULL, i, 1.0, 1.0, 1.0, 0, W, 1);
	  else
		  interpolate_row_rb_mul_pp (red[i], blue[i], green[i-1], green[i], green[i+1], i, 1.0, 1.0, 1.0, 0, W, 1);
  }
  if (plistener)
    plistener->setProgress (1.0);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

  float (*brow[5])[4], *pix;
  int prow=7, pcol=1, *ip, *code[16][16], gval[8], gmin, gmax, sum[4];
  int row, col, x, y, x1, x2, y1, y2, t, weight, grads, color, diag;
  int g, diff, thold, num, c, width=W, height=H, colors=4;
  float (*image)[4];
  int lcode[16][16][32], shift, i, j;

  image = (float (*)[4]) calloc (H*W, sizeof *image);
  for (int ii=0; ii<H; ii++)
    for (int jj=0; jj<W; jj++)
        image[ii*W+jj][fc(ii,jj)] = rawData[ii][jj];

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
	  sum[color] += (1 << shift);
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
	sum[ip[2]] += pix[ip[0]] * (1 << ip[1]);
      for (i=colors; --i; ip+=2)
	pix[ip[0]] = sum[ip[0]] * ip[1] /256;
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
	  if (grads & (1<<g)) *ip++ = g;
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
  brow[4] = (float (*)[4]) calloc (width*3, sizeof **brow);
  for (row=0; row < 3; row++)
    brow[row] = brow[4] + row*width;
  for (row=2; row < height-2; row++) {		/* Do VNG interpolation */
    for (col=2; col < width-2; col++) {
      color = fc(row,col);
      pix = image[row*width+col];
      ip = code[row & prow][col & pcol];
      memset (gval, 0, sizeof gval);
      while ((g = ip[0]) != INT_MAX) {		/* Calculate gradients */
	diff = ABS(pix[g] - pix[ip[1]]) * (1<< ip[2]);
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
      thold = gmin + (gmax /2);
      memset (sum, 0, sizeof sum);
      for (num=g=0; g < 8; g++,ip+=2) {		/* Average the neighbors */
	if (gval[g] <= thold) {
	  FORCC
	    if (c == color && ip[1])
	      sum[c] += (pix[c] + pix[ip[1]]) /2;
	    else
	      sum[c] += pix[ip[0] + c];
	  num++;
	}
      }
      FORCC {					/* Save to buffer */
	t = pix[color];
	if (c != color)
	  t += (sum[c] - sum[color]) / num;
	brow[2][col][c] = t;
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

  for (int i=0; i<H; i++) {
    for (int j=0; j<W; j++)
        green[i][j] = (image[i*W+j][1] + image[i*W+j][3]) /2;
  }
  // Interpolate R and B
  for (int i=0; i<H; i++) {
	  if (i==0)
		  // rm, gm, bm must be recovered
		  //interpolate_row_rb_mul_pp (red, blue, NULL, green[i], green[i+1], i, rm, gm, bm, 0, W, 1);
		  interpolate_row_rb_mul_pp (red[i], blue[i], NULL, green[i], green[i+1], i, 1.0, 1.0, 1.0, 0, W, 1);
	  else if (i==H-1)
		  interpolate_row_rb_mul_pp (red[i], blue[i], green[i-1], green[i], NULL, i, 1.0, 1.0, 1.0, 0, W, 1);
	  else
		  interpolate_row_rb_mul_pp (red[i], blue[i], green[i-1], green[i], green[i+1], i, 1.0, 1.0, 1.0, 0, W, 1);
  }
  free (image);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef fc
#define fc(row,col) \
	(ri->get_filters() >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)
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
  float (*pix)[4];

  float (*image)[4];
  int colors = 3;

  if (plistener) {
    plistener->setProgressStr ("Demosaicing...");
    plistener->setProgress (0.0);
  }

  image = (float (*)[4]) calloc (H*W, sizeof *image);
  for (int ii=0; ii<H; ii++)
    for (int jj=0; jj<W; jj++)
        image[ii*W+jj][fc(ii,jj)] = rawData[ii][jj];

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
	pix[0][c] = CLIP(0.5*(pix[-d][c] + pix[d][c] + 2*pix[0][1]
			- pix[-d][1] - pix[d][1]) );
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
	pix[0][c] = CLIP(guess[diff[0] > diff[1]] /2);
      else
	pix[0][c] = CLIP((guess[0]+guess[1]) /4);
    }
    if(plistener) plistener->setProgress(0.67 + 0.33*row/(height-1));
  }

  red = new float*[H];
  for (int i=0; i<H; i++) {
    red[i] = new float[W];
    for (int j=0; j<W; j++)
        red[i][j] = image[i*W+j][0];
  }
  green = new float*[H];
  for (int i=0; i<H; i++) {
    green[i] = new float[W];
    for (int j=0; j<W; j++)
        green[i][j] = image[i*W+j][1];
  }
  blue = new float*[H];
  for (int i=0; i<H; i++) {
    blue[i] = new float[W];
    for (int j=0; j<W; j++)
        blue[i][j] = image[i*W+j][2];
  }
  free (image);
}

void RawImageSource::border_interpolate(int border, float (*image)[4], int start, int end)
{
  unsigned row, col, y, x, f, c, sum[8];
  int width=W, height=H;
  int colors = 3;

  if (end == 0 )end = H;
  for (row=start; row < end; row++)
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

void RawImageSource::ahd_demosaic(int winx, int winy, int winw, int winh)
{
    int i, j, k, top, left, row, col, tr, tc, c, d, val, hm[2];
    float (*pix)[4], (*rix)[3];
    static const int dir[4] = { -1, 1, -TS, TS };
    float ldiff[2][4], abdiff[2][4], leps, abeps;
    float r, xyz[3], xyz_cam[3][4];
	float (*cbrt);
    float (*rgb)[TS][TS][3];
    float (*lab)[TS][TS][3];
	float (*lix)[3];
    char (*homo)[TS][TS], *buffer;

    int width=W, height=H;
    float (*image)[4];
    int colors = 3;

    const double xyz_rgb[3][3] = {			/* XYZ from RGB */
        { 0.412453, 0.357580, 0.180423 },
        { 0.212671, 0.715160, 0.072169 },
        { 0.019334, 0.119193, 0.950227 }
    };

    const float d65_white[3] = { 0.950456, 1, 1.088754 };
	const float d50_white[3] = { 0.96422,  1, 0.82521  };

    if (plistener) {
        plistener->setProgressStr ("AHD Demosaicing...");
        plistener->setProgress (0.0);
    }

    image = (float (*)[4]) calloc (H*W, sizeof *image);
    for (int ii=0; ii<H; ii++)
        for (int jj=0; jj<W; jj++)
            image[ii*W+jj][fc(ii,jj)] = rawData[ii][jj];

	cbrt = (float (*)) calloc (0x10000, sizeof *cbrt);
    for (i=0; i < 0x10000; i++) {
        r = (double)i / 65535.0;
        cbrt[i] = r > 0.008856 ? pow(r,1/3.0) : 7.787*r + 16/116.0;
    }

    for (i=0; i < 3; i++)
        for (j=0; j < colors; j++)
            for (xyz_cam[i][j] = k=0; k < 3; k++)
	            xyz_cam[i][j] += xyz_rgb[i][k] * rgb_cam[k][j] / d65_white[i];

    border_interpolate(5, image);
    buffer = (char *) malloc (13*TS*TS*sizeof(float));		/* 1664 kB */
    //merror (buffer, "ahd_interpolate()");
    rgb  = (float(*)[TS][TS][3]) buffer;
    lab  = (float(*)[TS][TS][3])(buffer + 6*TS*TS*sizeof(float));
    homo = (char (*)[TS][TS])   (buffer + 12*TS*TS*sizeof(float));

    // helper variables for progress indication
    int n_tiles = ((height-7 + (TS-7))/(TS-6)) * ((width-7 + (TS-7))/(TS-6));
    int tile = 0;

    for (top=2; top < height-5; top += TS-6)
        for (left=2; left < width-5; left += TS-6) {
            /*  Interpolate green horizontally and vertically:		*/
            for (row = top; row < top+TS && row < height-2; row++) {
	            col = left + (FC(row,left) & 1);
	            for (c = FC(row,col); col < left+TS && col < width-2; col+=2) {
	                pix = image + (row*width+col);
	                val = 0.25*((pix[-1][1] + pix[0][c] + pix[1][1]) * 2
		                  - pix[-2][c] - pix[2][c]) ;
	                rgb[0][row-top][col-left][1] = ULIM(val,pix[-1][1],pix[1][1]);
	                val = 0.25*((pix[-width][1] + pix[0][c] + pix[width][1]) * 2
		                  - pix[-2*width][c] - pix[2*width][c]) ;
	                rgb[1][row-top][col-left][1] = ULIM(val,pix[-width][1],pix[width][1]);
	            }
            }

            /*  Interpolate red and blue, and convert to CIELab:		*/
            for (d=0; d < 2; d++)
	            for (row=top+1; row < top+TS-1 && row < height-3; row++)
	                for (col=left+1; col < left+TS-1 && col < width-3; col++) {
	                    pix = image + (row*width+col);
	                    rix = &rgb[d][row-top][col-left];
	                    lix = &lab[d][row-top][col-left];
	                    if ((c = 2 - FC(row,col)) == 1) {
	                        c = FC(row+1,col);
	                        val = pix[0][1] + (0.5*( pix[-1][2-c] + pix[1][2-c]
				                  - rix[-1][1] - rix[1][1] ) );
	                        rix[0][2-c] = CLIP(val);
	                        val = pix[0][1] + (0.5*( pix[-width][c] + pix[width][c]
				                  - rix[-TS][1] - rix[TS][1] ) );
	                    } else
	                        val = rix[0][1] + (0.25*( pix[-width-1][c] + pix[-width+1][c]
				                  + pix[+width-1][c] + pix[+width+1][c]
				                  - rix[-TS-1][1] - rix[-TS+1][1]
				                  - rix[+TS-1][1] - rix[+TS+1][1] + 1) );
	                    rix[0][c] = CLIP(val);
	                    c = FC(row,col);
	                    rix[0][c] = pix[0][c];
	                    xyz[0] = xyz[1] = xyz[2] = 0.0;
	                    FORCC {
	                        xyz[0] += xyz_cam[0][c] * rix[0][c];
	                        xyz[1] += xyz_cam[1][c] * rix[0][c];
	                        xyz[2] += xyz_cam[2][c] * rix[0][c];
	                    }

						xyz[0] = CurveFactory::flinterp(cbrt,xyz[0]);
	                    xyz[1] = CurveFactory::flinterp(cbrt,xyz[1]);
	                    xyz[2] = CurveFactory::flinterp(cbrt,xyz[2]);

						//xyz[0] = xyz[0] > 0.008856 ? pow(xyz[0]/65535,1/3.0) : 7.787*xyz[0] + 16/116.0;
						//xyz[1] = xyz[1] > 0.008856 ? pow(xyz[1]/65535,1/3.0) : 7.787*xyz[1] + 16/116.0;
						//xyz[2] = xyz[2] > 0.008856 ? pow(xyz[2]/65535,1/3.0) : 7.787*xyz[2] + 16/116.0;

	                    lix[0][0] = (116 * xyz[1] - 16);
	                    lix[0][1] = 500 * (xyz[0] - xyz[1]);
	                    lix[0][2] = 200 * (xyz[1] - xyz[2]);
	                }

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
                            0.5*(rgb[0][tr][tc][c] + rgb[1][tr][tc][c]) ;
                }
            }

            tile++;
            if(plistener) {
                plistener->setProgress((double)tile / n_tiles);
            }
        }

    if(plistener) plistener->setProgress (1.0);
    free (buffer);
    for (int i=0; i<H; i++) {
        for (int j=0; j<W; j++){
            red[i][j] = image[i*W+j][0];
            green[i][j] = image[i*W+j][1];
            blue[i][j] = image[i*W+j][2];
        }
    }

    free (image);
	free (cbrt);
}
#undef TS

void RawImageSource::nodemosaic()
{
    red = new float*[H];
    green = new float*[H];
    blue = new float*[H];
    for (int i=0; i<H; i++) {
        red[i] = new float[W];
        green[i] = new float[W];
        blue[i] = new float[W];
        for (int j=0; j<W; j++){
        	switch( FC(i,j)){
				case 0: red[i][j] = rawData[i][j]; green[i][j]=blue[i][j]=0; break;
				case 1: green[i][j] = rawData[i][j]; red[i][j]=blue[i][j]=0; break;
				case 2: blue[i][j] = rawData[i][j]; red[i][j]=green[i][j]=0; break;
					//red[i][j] = rawData[i][j];
					//green[i][j] = rawData[i][j];
					//blue[i][j] = rawData[i][j];
        	}
        }
    }
}

// Refinement based on EECI demosaicing algorithm by L. Chang and Y.P. Tan
// from "Lassus" : Luis Sanz, adapted by Jacques Desmis - JDC - and Oliver Duis for RawTherapee
// increases the signal to noise ratio (PSNR) # +1 to +2 dB : tested with Dcraw : eg: Lighthouse + AMaZE : whitout refinement:39.96dB, with refinement:41.86 dB
// reduce color artifacts, improves the interpolation
// but it's relatively slow
void RawImageSource::refinement_lassus()
{  
    const int PassCount=2;  // two passes EECI refine, slow but best results...

    if (settings->verbose) printf("Refinement Lassus\n");

    MyTime t1e,t2e;
    t1e.set();
    int u=W, v=2*u, w=3*u, x=4*u, y=5*u;
    float (*image)[4];
    image = (float(*)[4]) calloc(W*H, sizeof *image);
#pragma omp parallel shared(image)
    {
        // convert red, blue, green to image
#pragma omp for
        for (int i=0;i<H;i++) {
            for (int j=0;j<W;j++) {
                image[i*W+j][0] = red  [i][j];
                image[i*W+j][1] = green[i][j];
                image[i*W+j][2] = blue [i][j];			 
            }
        }

        for (int b=0; b<PassCount; b++) {
            if (plistener) {
                plistener->setProgressStr ("Refinement...");
                plistener->setProgress ((float)b/PassCount);
            }

            // Reinforce interpolated green pixels on RED/BLUE pixel locations
#pragma omp for
            for (int row=6; row<H-6; row++) {
                int c,d;
                for (int col=6+(FC(row,2)&1),c=FC(row,col); col<W-6; col+=2) {
                    float (*pix)[4]=image+row*W+col;

                    // Cubic Spline Interpolation by Li and Randhawa, modified by Luis Sanz Rodríguez

                    float f[4];
                    f[0]=1.0/(1.0+2.0*fabs(1.125*pix[-v][c]-0.875*pix[0][c]-0.250*pix[-x][c])+fabs(0.875*pix[u][1]-1.125*pix[-u][1]+0.250*pix[-w][1])+fabs(0.875*pix[-w][1]-1.125*pix[-u][1]+0.250*pix[-y][1]));   
                    f[1]=1.0/(1.0+2.0*fabs(1.125*pix[+2][c]-0.875*pix[0][c]-0.250*pix[+4][c])+fabs(0.875*pix[1][1]-1.125*pix[-1][1]+0.250*pix[+3][1])+fabs(0.875*pix[+3][1]-1.125*pix[+1][1]+0.250*pix[+5][1]));
                    f[2]=1.0/(1.0+2.0*fabs(1.125*pix[-2][c]-0.875*pix[0][c]-0.250*pix[-4][c])+fabs(0.875*pix[1][1]-1.125*pix[-1][1]+0.250*pix[-3][1])+fabs(0.875*pix[-3][1]-1.125*pix[-1][1]+0.250*pix[-5][1]));
                    f[3]=1.0/(1.0+2.0*fabs(1.125*pix[+v][c]-0.875*pix[0][c]-0.250*pix[+x][c])+fabs(0.875*pix[u][1]-1.125*pix[-u][1]+0.250*pix[+w][1])+fabs(0.875*pix[+w][1]-1.125*pix[+u][1]+0.250*pix[+y][1])); 

                    float g[4];//CLIREF avoid overflow
                    g[0]=pix[0][c]+(0.875*CLIREF(pix[-u][1]-pix[-u][c])+0.125*CLIREF(pix[+u][1]-pix[+u][c]));
                    g[1]=pix[0][c]+(0.875*CLIREF(pix[+1][1]-pix[+1][c])+0.125*CLIREF(pix[-1][1]-pix[-1][c]));
                    g[2]=pix[0][c]+(0.875*CLIREF(pix[-1][1]-pix[-1][c])+0.125*CLIREF(pix[+1][1]-pix[+1][c]));
                    g[3]=pix[0][c]+(0.875*CLIREF(pix[+u][1]-pix[+u][c])+0.125*CLIREF(pix[-u][1]-pix[-u][c]));
					
					
                    pix[0][1]=(f[0]*g[0]+f[1]*g[1]+f[2]*g[2]+f[3]*g[3]) / (f[0]+f[1]+f[2]+f[3]);

                }
            }
            // Reinforce interpolated red/blue pixels on GREEN pixel locations
#pragma omp for 
            for (int row=6; row<H-6; row++) {
                int c,d;
                for (int col=6+(FC(row,3)&1),c=FC(row,col+1); col<W-6; col+=2) {
                    float (*pix)[4]=image+row*W+col;
                    for (int i=0; i<2; c=2-c,i++) {
                        float f[4];
                        f[0]=1.0/(1.0+2.0*fabs(0.875*pix[-v][1]-1.125*pix[0][1]+0.250*pix[-x][1])+fabs(pix[u] [c]-pix[-u][c])+fabs(pix[-w][c]-pix[-u][c]));
                        f[1]=1.0/(1.0+2.0*fabs(0.875*pix[+2][1]-1.125*pix[0][1]+0.250*pix[+4][1])+fabs(pix[+1][c]-pix[-1][c])+fabs(pix[+3][c]-pix[+1][c]));
                        f[2]=1.0/(1.0+2.0*fabs(0.875*pix[-2][1]-1.125*pix[0][1]+0.250*pix[-4][1])+fabs(pix[+1][c]-pix[-1][c])+fabs(pix[-3][c]-pix[-1][c]));
                        f[3]=1.0/(1.0+2.0*fabs(0.875*pix[+v][1]-1.125*pix[0][1]+0.250*pix[+x][1])+fabs(pix[u] [c]-pix[-u][c])+fabs(pix[+w][c]-pix[+u][c]));

                        float g[5];//CLIREF avoid overflow
                        g[0]=CLIREF(pix[-u][1]-pix[-u][c]);
                        g[1]=CLIREF(pix[+1][1]-pix[+1][c]);
                        g[2]=CLIREF(pix[-1][1]-pix[-1][c]);
                        g[3]=CLIREF(pix[+u][1]-pix[+u][c]);
                        g[4]=((f[0]*g[0]+f[1]*g[1]+f[2]*g[2]+f[3]*g[3]) / (f[0]+f[1]+f[2]+f[3]));
                        pix[0][c]= pix[0][1]-(0.65*g[4]+0.35*CLIREF(pix[0][1]-pix[0][c]));
                    }
                }
				}
            // Reinforce integrated red/blue pixels on BLUE/RED pixel locations
#pragma omp for
            for (int row=6; row<H-6; row++) {
                int c,d;
                for (int col=6+(FC(row,2)&1),c=2-FC(row,col),d=2-c; col<W-6; col+=2) {
                    float (*pix)[4]=image+row*W+col;

                    float f[4];
                    f[0]=1.0/(1.0+2.0*fabs(1.125*pix[-v][d]-0.875*pix[0][d]-0.250*pix[-x][d])+fabs(0.875*pix[u][1]-1.125*pix[-u][1]+0.250*pix[-w][1])+fabs(0.875*pix[-w][1]-1.125*pix[-u][1]+0.250*pix[-y][1]));   
                    f[1]=1.0/(1.0+2.0*fabs(1.125*pix[+2][d]-0.875*pix[0][d]-0.250*pix[+4][d])+fabs(0.875*pix[1][1]-1.125*pix[-1][1]+0.250*pix[+3][1])+fabs(0.875*pix[+3][1]-1.125*pix[+1][1]+0.250*pix[+5][1]));
                    f[2]=1.0/(1.0+2.0*fabs(1.125*pix[-2][d]-0.875*pix[0][d]-0.250*pix[-4][d])+fabs(0.875*pix[1][1]-1.125*pix[-1][1]+0.250*pix[-3][1])+fabs(0.875*pix[-3][1]-1.125*pix[-1][1]+0.250*pix[-5][1]));
                    f[3]=1.0/(1.0+2.0*fabs(1.125*pix[+v][d]-0.875*pix[0][d]-0.250*pix[+x][d])+fabs(0.875*pix[u][1]-1.125*pix[-u][1]+0.250*pix[+w][1])+fabs(0.875*pix[+w][1]-1.125*pix[+u][1]+0.250*pix[+y][1])); 

                    float g[5];
                    g[0]=(0.875*(pix[-u][1]-pix[-u][c])+0.125*(pix[-v][1]-pix[-v][c]));
                    g[1]=(0.875*(pix[+1][1]-pix[+1][c])+0.125*(pix[+2][1]-pix[+2][c]));
                    g[2]=(0.875*(pix[-1][1]-pix[-1][c])+0.125*(pix[-2][1]-pix[-2][c]));
                    g[3]=(0.875*(pix[+u][1]-pix[+u][c])+0.125*(pix[+v][1]-pix[+v][c]));
					
                    g[4]=(f[0]*g[0]+f[1]*g[1]+f[2]*g[2]+f[3]*g[3]) / (f[0]+f[1]+f[2]+f[3]);

                    float p[9];
                    p[0]=(pix[-u-1][1]-pix[-u-1][c]);
                    p[1]=(pix[-u+0][1]-pix[-u+0][c]);
                    p[2]=(pix[-u+1][1]-pix[-u+1][c]);
                    p[3]=(pix[+0-1][1]-pix[+0-1][c]);
                    p[4]=(pix[+0+0][1]-pix[+0+0][c]);
                    p[5]=(pix[+0+1][1]-pix[+0+1][c]);
                    p[6]=(pix[+u-1][1]-pix[+u-1][c]);
                    p[7]=(pix[+u+0][1]-pix[+u+0][c]);
                    p[8]=(pix[+u+1][1]-pix[+u+1][c]);

                    // sort p[]
                    float temp;  // used in PIX_SORT macro;
                    PIX_SORT(p[1],p[2]); PIX_SORT(p[4],p[5]); PIX_SORT(p[7],p[8]); PIX_SORT(p[0],p[1]); PIX_SORT(p[3],p[4]); PIX_SORT(p[6],p[7]); PIX_SORT(p[1],p[2]); PIX_SORT(p[4],p[5]); PIX_SORT(p[7],p[8]); PIX_SORT(p[0],p[3]); PIX_SORT(p[5],p[8]); PIX_SORT(p[4],p[7]); PIX_SORT(p[3],p[6]); PIX_SORT(p[1],p[4]); PIX_SORT(p[2],p[5]); PIX_SORT(p[4],p[7]); PIX_SORT(p[4],p[2]); PIX_SORT(p[6],p[4]); PIX_SORT(p[4],p[2]);
                    pix[0][c]=LIM(pix[0][1]-(1.30*g[4]-0.30*(pix[0][1]-pix[0][c])), 0.99*(pix[0][1]-p[4]), 1.01*(pix[0][1]-p[4]));

                }
            }

        }

        // put modified values to red, green, blue
#pragma omp for
        for (int i=0;i<H;i++) {
            for (int j=0; j<W; j++) {
                red  [i][j] =image[i*W+j][0];
                green[i][j] =image[i*W+j][1];
                blue [i][j] =image[i*W+j][2];			 
            }
        }  
    }

    free(image);

    t2e.set();
    if (settings->verbose) printf("Refinement %d usec\n", t2e.etime(t1e));
}

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
	if( y0+TILESIZE+TILEBORDER >= H-border) rowMax = TILEBORDER+H-border-y0;
	if( x0+TILESIZE+TILEBORDER >= W-border) colMax = TILEBORDER+W-border-x0;
}

void RawImageSource::fill_raw( float (*cache )[4], int x0, int y0, float** rawData)
{
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,0);

    for (int row=rowMin,y=y0-TILEBORDER+rowMin; row<rowMax; row++,y++)
    	for (int col=colMin,x=x0-TILEBORDER+colMin,indx=row*CACHESIZE+col; col<colMax; col++,x++,indx++){
    		cache[indx][fc(y,x)] = rawData[y][x];
    	}
}

void RawImageSource::fill_border( float (*cache )[4], int border, int x0, int y0)
{
	unsigned row, col, y, x, f, c;
    float sum[8];
	const int colors = 3;  // used in FORCC

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
				if (c != f && sum[c + 4]>0)
					cache[(row-y0+TILEBORDER) * CACHESIZE +TILEBORDER + col-x0][c] = sum[c] / sum[c + 4];
		}
	}
}
// saves red and blue
void RawImageSource::copy_to_buffer( float (*buffer)[3], float (*image)[4])
{
	for (int indx=0; indx < CACHESIZE*CACHESIZE; indx++) {
		buffer[indx][0]=image[indx][0]; //R
		buffer[indx][2]=image[indx][2]; //B
	}
}

// restores red and blue
void RawImageSource::restore_from_buffer(float (*image)[4], float (*buffer)[3])
{
	for (int indx=0; indx < CACHESIZE*CACHESIZE; indx++) {
		image[indx][0]=buffer[indx][0]; //R
		image[indx][2]=buffer[indx][2]; //B
	}
}

// First pass green interpolation
void RawImageSource::dcb_hid(float (*image)[4],float (*bufferH)[3], float (*bufferV)[3], int x0, int y0)
{
	const int u=CACHESIZE, v=2*CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,2);

	// green pixels
	for (int row = rowMin; row < rowMax; row++) {
		for (int col = colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin)&1),indx=row*CACHESIZE+col,c=FC(y0-TILEBORDER+row,x0-TILEBORDER+col); col < colMax; col+=2, indx+=2) {
            assert(indx>=0 && indx<u*u);
			bufferH[indx][1] = (image[indx-1][1] + image[indx+1][1]) * 0.5f;
			bufferV[indx][1] = (image[indx+u][1] + image[indx-u][1]) * 0.5f;
		}
	}
	// red in blue pixel, blue in red pixel
	for (int row=rowMin; row < rowMax; row++)
		for (int col=colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin) & 1), indx=row*CACHESIZE+col, c=2-FC(y0-TILEBORDER+row,x0-TILEBORDER+col); col < colMax; col+=2, indx+=2) {
            assert(indx>=0 && indx<u*u && c>=0 && c<3);

            bufferH[indx][c] = ( 4.f * bufferH[indx][1]
			                - bufferH[indx+u+1][1] - bufferH[indx+u-1][1] - bufferH[indx-u+1][1] - bufferH[indx-u-1][1]
			                + image[indx+u+1][c] + image[indx+u-1][c] + image[indx-u+1][c] + image[indx-u-1][c] ) * 0.25f;
			bufferV[indx][c] = ( 4.f * bufferV[indx][1]
						    - bufferV[indx+u+1][1] - bufferV[indx+u-1][1] - bufferV[indx-u+1][1] - bufferV[indx-u-1][1]
						    + image[indx+u+1][c] + image[indx+u-1][c] + image[indx-u+1][c] + image[indx-u-1][c] ) * 0.25f;
		}

	// red or blue in green pixels
	for (int row=rowMin; row<rowMax; row++)
		for (int col=colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin+1)&1), indx=row*CACHESIZE+col,c=FC(y0-TILEBORDER+row,x0-TILEBORDER+col+1),d=2-c; col<colMax; col+=2, indx+=2) {
            assert(indx>=0 && indx<u*u && c>=0 && c<3 && d>=0 && d<3);
			bufferH[indx][c] = (image[indx+1][c] + image[indx-1][c]) * 0.5f;
			bufferH[indx][d] = (2.f * bufferH[indx][1] - bufferH[indx+u][1] - bufferH[indx-u][1] + image[indx+u][d] + image[indx-u][d]) * 0.5f;
			bufferV[indx][c] = (2.f * bufferV[indx][1] - bufferV[indx+1][1] - bufferV[indx-1][1] + image[indx+1][c] + image[indx-1][c]) * 0.5f;
			bufferV[indx][d] = (image[indx+u][d] + image[indx-u][d]) * 0.5f;
		}

    // Decide green pixels
    for (int row = rowMin; row < rowMax; row++)
        for (int col = colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin)&1),indx=row*CACHESIZE+col,c=FC(y0-TILEBORDER+row,x0-TILEBORDER+col),d=2-c; col < colMax; col+=2, indx+=2) {
            float current =  MAX(image[indx+v][c], MAX(image[indx-v][c], MAX(image[indx-2][c], image[indx+2][c]))) -
                            MIN(image[indx+v][c], MIN(image[indx-v][c], MIN(image[indx-2][c], image[indx+2][c]))) +
                            MAX(image[indx+1+u][d], MAX(image[indx+1-u][d], MAX(image[indx-1+u][d], image[indx-1-u][d]))) -
                            MIN(image[indx+1+u][d], MIN(image[indx+1-u][d], MIN(image[indx-1+u][d], image[indx-1-u][d])));

            float currentH = MAX(bufferH[indx+v][d], MAX(bufferH[indx-v][d], MAX(bufferH[indx-2][d], bufferH[indx+2][d]))) -
                            MIN(bufferH[indx+v][d], MIN(bufferH[indx-v][d], MIN(bufferH[indx-2][d], bufferH[indx+2][d]))) +
                            MAX(bufferH[indx+1+u][c], MAX(bufferH[indx+1-u][c], MAX(bufferH[indx-1+u][c], bufferH[indx-1-u][c]))) -
                            MIN(bufferH[indx+1+u][c], MIN(bufferH[indx+1-u][c], MIN(bufferH[indx-1+u][c], bufferH[indx-1-u][c])));

            float currentV = MAX(bufferV[indx+v][d], MAX(bufferV[indx-v][d], MAX(bufferV[indx-2][d], bufferV[indx+2][d]))) -
                            MIN(bufferV[indx+v][d], MIN(bufferV[indx-v][d], MIN(bufferV[indx-2][d], bufferV[indx+2][d]))) +
                            MAX(bufferV[indx+1+u][c], MAX(bufferV[indx+1-u][c], MAX(bufferV[indx-1+u][c], bufferV[indx-1-u][c]))) -
                            MIN(bufferV[indx+1+u][c], MIN(bufferV[indx+1-u][c], MIN(bufferV[indx-1+u][c], bufferV[indx-1-u][c])));

            assert(indx>=0 && indx<u*u);
            if (ABS(current-currentH) < ABS(current-currentV))
                image[indx][1] = bufferH[indx][1];
            else
                image[indx][1] = bufferV[indx][1];
        }
}

// missing colors are interpolated
void RawImageSource::dcb_color(float (*image)[4], int x0, int y0)
{
	const int u=CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,1);

	// red in blue pixel, blue in red pixel
	for (int row=rowMin; row < rowMax; row++)
		for (int col=colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin) & 1), indx=row*CACHESIZE+col, c=2-FC(y0-TILEBORDER+row,x0-TILEBORDER+col); col < colMax; col+=2, indx+=2) {
            assert(indx>=0 && indx<u*u && c>=0 && c<4);
			image[indx][c] = ( 4.f * image[indx][1]
			                - image[indx+u+1][1] - image[indx+u-1][1] - image[indx-u+1][1] - image[indx-u-1][1]
			                + image[indx+u+1][c] + image[indx+u-1][c] + image[indx-u+1][c] + image[indx-u-1][c] ) * 0.25f;
		}

	// red or blue in green pixels
	for (int row=rowMin; row<rowMax; row++)
		for (int col=colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin+1)&1), indx=row*CACHESIZE+col,c=FC(y0-TILEBORDER+row,x0-TILEBORDER+col+1),d=2-c; col<colMax; col+=2, indx+=2) {
            assert(indx>=0 && indx<u*u && c>=0 && c<4);
            image[indx][c] = (2.f * image[indx][1] - image[indx+1][1] - image[indx-1][1] + image[indx+1][c] + image[indx-1][c]) * 0.5f;
			image[indx][d] = (2.f * image[indx][1] - image[indx+u][1] - image[indx-u][1] + image[indx+u][d] + image[indx-u][d]) * 0.5f;
		}
}

// green correction
void RawImageSource::dcb_hid2(float (*image)[4], int x0, int y0)
{
	const int u=CACHESIZE, v=2*CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,2);

	for (int row=rowMin; row < rowMax; row++) {
		for (int col = colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin)&1),indx=row*CACHESIZE+col,c=FC(y0-TILEBORDER+row,x0-TILEBORDER+col); col < colMax; col+=2, indx+=2) {
            assert(indx>=0 && indx<u*u);
            image[indx][1] = (image[indx+v][1] + image[indx-v][1] + image[indx-2][1] + image[indx+2][1]) * 0.25f +
						   image[indx][c] - ( image[indx+v][c] + image[indx-v][c] + image[indx-2][c] + image[indx+2][c]) * 0.25f;
		}
	}
}

// green is used to create
// an interpolation direction map
// 1 = vertical
// 0 = horizontal
// saved in image[][3]
void RawImageSource::dcb_map(float (*image)[4], int x0, int y0)
{
	const int u=4*CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,2);

	for (int row=rowMin; row < rowMax; row++) {
		for (int col=colMin, indx=row*CACHESIZE+col; col < colMax; col++, indx++) {
			float *pix = &(image[indx][1]);

            assert(indx>=0 && indx<u*u);
			if ( *pix > ( pix[-4] + pix[+4] + pix[-u] + pix[+u])/4 )
				image[indx][3] = ((MIN( pix[-4], pix[+4]) + pix[-4] + pix[+4] ) < (MIN( pix[-u], pix[+u]) + pix[-u] + pix[+u]));
			else
				image[indx][3] = ((MAX( pix[-4], pix[+4]) + pix[-4] + pix[+4] ) > (MAX( pix[-u], pix[+u]) + pix[-u] + pix[+u]));
		}
	}
}

// interpolated green pixels are corrected using the map
void RawImageSource::dcb_correction(float (*image)[4], int x0, int y0)
{
	const int u=CACHESIZE, v=2*CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,2);

	for (int row=rowMin; row < rowMax; row++) {
		for (int col = colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin)&1),indx=row*CACHESIZE+col,c=FC(y0-TILEBORDER+row,x0-TILEBORDER+col); col < colMax; col+=2, indx+=2) {
			float current = 4.f * image[indx][3] +
						  2.f * (image[indx+u][3] + image[indx-u][3] + image[indx+1][3] + image[indx-1][3]) +
							image[indx+v][3] + image[indx-v][3] + image[indx+2][3] + image[indx-2][3];

            assert(indx>=0 && indx<u*u);
			image[indx][1] = ((16.f-current)*(image[indx-1][1] + image[indx+1][1]) * 0.5f + current*(image[indx-u][1] + image[indx+u][1]) * 0.5f ) * 0.0625f;
		}
	}
}

// R and B smoothing using green contrast, all pixels except 2 pixel wide border
void RawImageSource::dcb_pp(float (*image)[4], int x0, int y0)
{
	const int u=CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,2);

	for (int row=rowMin; row < rowMax; row++)
		for (int col=colMin, indx=row*CACHESIZE+col; col < colMax; col++, indx++) {
			//int r1 = ( image[indx-1][0] + image[indx+1][0] + image[indx-u][0] + image[indx+u][0] + image[indx-u-1][0] + image[indx+u+1][0] + image[indx-u+1][0] + image[indx+u-1][0])/8;
			//int g1 = ( image[indx-1][1] + image[indx+1][1] + image[indx-u][1] + image[indx+u][1] + image[indx-u-1][1] + image[indx+u+1][1] + image[indx-u+1][1] + image[indx+u-1][1])/8;
			//int b1 = ( image[indx-1][2] + image[indx+1][2] + image[indx-u][2] + image[indx+u][2] + image[indx-u-1][2] + image[indx+u+1][2] + image[indx-u+1][2] + image[indx+u-1][2])/8;
			float (*pix)[4] = image+(indx-u-1);
			float r1 = (*pix)[0];
			float g1 = (*pix)[1];
			float b1 = (*pix)[2];
			pix++;
			r1 += (*pix)[0];
			g1 += (*pix)[1];
			b1 += (*pix)[2];
			pix++;
			r1 += (*pix)[0];
			g1 += (*pix)[1];
			b1 += (*pix)[2];
			pix+=CACHESIZE-2;
			r1 += (*pix)[0];
			g1 += (*pix)[1];
			b1 += (*pix)[2];
			pix+=2;
			r1 += (*pix)[0];
			g1 += (*pix)[1];
			b1 += (*pix)[2];
			pix+=CACHESIZE-2;
			r1 += (*pix)[0];
			g1 += (*pix)[1];
			b1 += (*pix)[2];
			pix++;
			r1 += (*pix)[0];
			g1 += (*pix)[1];
			b1 += (*pix)[2];
			pix++;
			r1 += (*pix)[0];
			g1 += (*pix)[1];
			b1 += (*pix)[2];
			r1 *=0.125f;
			g1 *=0.125f;
			b1 *=0.125f;
			r1 = r1 + ( image[indx][1] - g1 );
			b1 = b1 + ( image[indx][1] - g1 );

            assert(indx>=0 && indx<u*u);
			image[indx][0] = r1;
			image[indx][2] = b1;
		}
}

// interpolated green pixels are corrected using the map
// with correction
void RawImageSource::dcb_correction2(float (*image)[4], int x0, int y0)
{
	const int u=CACHESIZE, v=2*CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,4);

	for (int row=rowMin; row < rowMax; row++) {
		for (int col = colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin)&1),indx=row*CACHESIZE+col,c=FC(y0-TILEBORDER+row,x0-TILEBORDER+col); col < colMax; col+=2, indx+=2) {
			float current = 4.f * image[indx][3] +
						  2.f * (image[indx+u][3] + image[indx-u][3] + image[indx+1][3] + image[indx-1][3]) +
							image[indx+v][3] + image[indx-v][3] + image[indx+2][3] + image[indx-2][3];

            assert(indx>=0 && indx<u*u);
			image[indx][1] = ((16.f-current)*((image[indx-1][1] + image[indx+1][1]) * 0.5f
                + image[indx][c] - (image[indx+2][c] + image[indx-2][c]) * 0.5f)
                + current*((image[indx-u][1] + image[indx+u][1]) * 0.5f + image[indx][c] - (image[indx+v][c] + image[indx-v][c]) * 0.5f)) * 0.0625f;
		}
	}
}

// image refinement
void RawImageSource::dcb_refinement(float (*image)[4], int x0, int y0)
{
	const int u=CACHESIZE, v=2*CACHESIZE, w=3*CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,4);

	float f[5],g1,g2;

	for (int row=rowMin; row < rowMax; row++)
		for (int col=colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin)&1),indx=row*CACHESIZE+col,c=FC(y0-TILEBORDER+row,x0-TILEBORDER+col); col < colMax; col+=2,indx+=2){
			float current = 4.f * image[indx][3] +
				      2.f * (image[indx+u][3] + image[indx-u][3] + image[indx+1][3] + image[indx-1][3])
				      +image[indx+v][3] + image[indx-v][3] + image[indx-2][3] + image[indx+2][3];

			f[0] = (float)(image[indx-u][1] + image[indx+u][1])/(2.f + 2.f * image[indx][c]);
			f[1] = 2.f * image[indx-u][1]/(2 + image[indx-v][c] + image[indx][c]);
			f[2] = (float)(image[indx-u][1] + image[indx-w][1])/(2.f + 2.f * image[indx-v][c]);
			f[3] = 2.f * image[indx+u][1]/(2 + image[indx+v][c] + image[indx][c]);
			f[4] = (float)(image[indx+u][1] + image[indx+w][1])/(2.f + 2.f * image[indx+v][c]);

			g1 = (f[0] + f[1] + f[2] + f[3] + f[4] - MAX(f[1], MAX(f[2], MAX(f[3], f[4]))) - MIN(f[1], MIN(f[2], MIN(f[3], f[4])))) / 3.f;

			f[0] = (float)(image[indx-1][1] + image[indx+1][1])/(2.f + 2.f * image[indx][c]);
			f[1] = 2.f * image[indx-1][1]/(2 + image[indx-2][c] + image[indx][c]);
			f[2] = (float)(image[indx-1][1] + image[indx-3][1])/(2.f + 2.f * image[indx-2][c]);
			f[3] = 2.f * image[indx+1][1]/(2 + image[indx+2][c] + image[indx][c]);
			f[4] = (float)(image[indx+1][1] + image[indx+3][1])/(2.f + 2.f * image[indx+2][c]);

			g2 = (f[0] + f[1] + f[2] + f[3] + f[4] - MAX(f[1], MAX(f[2], MAX(f[3], f[4]))) - MIN(f[1], MIN(f[2], MIN(f[3], f[4])))) / 3.f;

            assert(indx>=0 && indx<u*u);
			image[indx][1] = (2.f+image[indx][c]) * (current*g1 + (16.f-current)*g2) * 0.0625f;

            // get rid of the overshooted pixels
		    float min = MIN(image[indx+1+u][1], MIN(image[indx+1-u][1], MIN(image[indx-1+u][1], MIN(image[indx-1-u][1], MIN(image[indx-1][1], MIN(image[indx+1][1], MIN(image[indx-u][1], image[indx+u][1])))))));
		    float max = MAX(image[indx+1+u][1], MAX(image[indx+1-u][1], MAX(image[indx-1+u][1], MAX(image[indx-1-u][1], MAX(image[indx-1][1], MAX(image[indx+1][1], MAX(image[indx-u][1], image[indx+u][1])))))));

		    image[indx][1] =  LIM(image[indx][1], min, max);
		}
}

// missing colors are interpolated using high quality algorithm by Luis Sanz Rodr‚àö‚â†guez
void RawImageSource::dcb_color_full(float (*image)[4], int x0, int y0, float (*chroma)[2])
{
	const int u=CACHESIZE, v=2*CACHESIZE, w=3*CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,3);

	int i,j;
	float f[4],g[4];

	for (int row=1; row < CACHESIZE-1; row++)
		for (int col=1+(FC(y0-TILEBORDER+row,x0-TILEBORDER+1)&1),indx=row*CACHESIZE+col,c=FC(y0-TILEBORDER+row,x0-TILEBORDER+col),d=c/2; col < CACHESIZE-1; col+=2,indx+=2) {
            assert(indx>=0 && indx<u*u && c>=0 && c<4);
			chroma[indx][d]=image[indx][c]-image[indx][1];
        }

	for (int row=rowMin; row<rowMax; row++)
		for (int col=colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin)&1),indx=row*CACHESIZE+col,c=1-FC(y0-TILEBORDER+row,x0-TILEBORDER+col)/2,d=1-c; col<colMax; col+=2,indx+=2) {
			f[0]=1.f/(float)(1.f+fabs(chroma[indx-u-1][c]-chroma[indx+u+1][c])+fabs(chroma[indx-u-1][c]-chroma[indx-w-3][c])+fabs(chroma[indx+u+1][c]-chroma[indx-w-3][c]));
			f[1]=1.f/(float)(1.f+fabs(chroma[indx-u+1][c]-chroma[indx+u-1][c])+fabs(chroma[indx-u+1][c]-chroma[indx-w+3][c])+fabs(chroma[indx+u-1][c]-chroma[indx-w+3][c]));
			f[2]=1.f/(float)(1.f+fabs(chroma[indx+u-1][c]-chroma[indx-u+1][c])+fabs(chroma[indx+u-1][c]-chroma[indx+w+3][c])+fabs(chroma[indx-u+1][c]-chroma[indx+w-3][c]));
			f[3]=1.f/(float)(1.f+fabs(chroma[indx+u+1][c]-chroma[indx-u-1][c])+fabs(chroma[indx+u+1][c]-chroma[indx+w-3][c])+fabs(chroma[indx-u-1][c]-chroma[indx+w+3][c]));
			g[0]=1.325f * chroma[indx-u-1][c] - 0.175f*chroma[indx-w-3][c] - 0.075f*chroma[indx-w-1][c] - 0.075f*chroma[indx-u-3][c];
			g[1]=1.325f * chroma[indx-u+1][c] - 0.175f*chroma[indx-w+3][c] - 0.075f*chroma[indx-w+1][c] - 0.075f*chroma[indx-u+3][c];
			g[2]=1.325f * chroma[indx+u-1][c] - 0.175f*chroma[indx+w-3][c] - 0.075f*chroma[indx+w-1][c] - 0.075f*chroma[indx+u-3][c];
			g[3]=1.325f * chroma[indx+u+1][c] - 0.175f*chroma[indx+w+3][c] - 0.075f*chroma[indx+w+1][c] - 0.075f*chroma[indx+u+3][c];

            assert(indx>=0 && indx<u*u && c>=0 && c<2);
			chroma[indx][c]=(f[0]*g[0]+f[1]*g[1]+f[2]*g[2]+f[3]*g[3])/(f[0]+f[1]+f[2]+f[3]);
		}
	for (int row=rowMin; row<rowMax; row++)
		for (int col=colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin+1)&1),indx=row*CACHESIZE+col,c=FC(y0-TILEBORDER+row,x0-TILEBORDER+col+1)/2; col<colMax; col+=2,indx+=2)
			for(int d=0;d<=1;c=1-c,d++){
				f[0]= 1.f/(float)(1.f + fabs(chroma[indx-u][c]-chroma[indx+u][c])+fabs(chroma[indx-u][c]-chroma[indx-w][c])+fabs(chroma[indx+u][c]-chroma[indx-w][c]));
				f[1]= 1.f/(float)(1.f + fabs(chroma[indx+1][c]-chroma[indx-1][c])+fabs(chroma[indx+1][c]-chroma[indx+3][c])+fabs(chroma[indx-1][c]-chroma[indx+3][c]));
				f[2]= 1.f/(float)(1.f + fabs(chroma[indx-1][c]-chroma[indx+1][c])+fabs(chroma[indx-1][c]-chroma[indx-3][c])+fabs(chroma[indx+1][c]-chroma[indx-3][c]));
				f[3]= 1.f/(float)(1.f + fabs(chroma[indx+u][c]-chroma[indx-u][c])+fabs(chroma[indx+u][c]-chroma[indx+w][c])+fabs(chroma[indx-u][c]-chroma[indx+w][c]));

				g[0]= 0.875f * chroma[indx-u][c] + 0.125f * chroma[indx-w][c];
				g[1]= 0.875f * chroma[indx+1][c] + 0.125f * chroma[indx+3][c];
				g[2]= 0.875f * chroma[indx-1][c] + 0.125f * chroma[indx-3][c];
				g[3]= 0.875f * chroma[indx+u][c] + 0.125f * chroma[indx+w][c];

                assert(indx>=0 && indx<u*u && c>=0 && c<2);
				chroma[indx][c]=(f[0]*g[0]+f[1]*g[1]+f[2]*g[2]+f[3]*g[3])/(f[0]+f[1]+f[2]+f[3]);
			}

	for(int row=rowMin; row<rowMax; row++)
		for(int col=colMin,indx=row*CACHESIZE+col; col<colMax; col++,indx++){
            assert(indx>=0 && indx<u*u);

			image[indx][0] = chroma[indx][0] + image[indx][1];
			image[indx][2] = chroma[indx][1] + image[indx][1];
		}
}

// DCB demosaicing main routine (sharp version)
void RawImageSource::dcb_demosaic(int iterations, bool dcb_enhance)
{
    double currentProgress=0.0;
    if(plistener) {
        plistener->setProgressStr ("DCB Demosaicing...");
        plistener->setProgress (currentProgress);
    }

    int wTiles = W/TILESIZE + (W%TILESIZE?1:0);
    int hTiles = H/TILESIZE + (H%TILESIZE?1:0);
    int numTiles = wTiles * hTiles;
    int tilesDone=0;
#ifdef _OPENMP
	int nthreads = omp_get_max_threads();
 	float (**image)[4]  =  (float(**)[4]) calloc( nthreads,sizeof( void*) );
	float (**image2)[3] =	(float(**)[3]) calloc( nthreads,sizeof( void*) );
	float (**image3)[3] =	(float(**)[3]) calloc( nthreads,sizeof( void*) );
	float  (**chroma)[2] =  (float (**)[2]) calloc( nthreads,sizeof( void*) );
	for(int i=0; i<nthreads; i++){
		image[i] = (float(*)[4]) calloc( CACHESIZE*CACHESIZE, sizeof **image);
		image2[i]= (float(*)[3]) calloc( CACHESIZE*CACHESIZE, sizeof **image2);
		image3[i]= (float(*)[3]) calloc( CACHESIZE*CACHESIZE, sizeof **image3);
		chroma[i]= (float (*)[2]) calloc( CACHESIZE*CACHESIZE, sizeof **chroma);
	}
#else
	float (*image)[4]  = (float(*)[4]) calloc( CACHESIZE*CACHESIZE, sizeof *image);
	float (*image2)[3] = (float(*)[3]) calloc( CACHESIZE*CACHESIZE, sizeof *image2);
	float (*image3)[3] = (float(*)[3]) calloc( CACHESIZE*CACHESIZE, sizeof *image3);
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
    	float (*tile)[4]   = image[tid];
    	float (*buffer)[3] = image2[tid];
    	float (*buffer2)[3]= image3[tid];
    	float  (*chrm)[2]   = chroma[tid];
#else
    	float (*tile)[4]   = image;
    	float (*buffer)[3] = image2;
    	float (*buffer2)[3]= image3;
    	float  (*chrm)[2]   = chroma;
#endif

		fill_raw( tile, x0,y0,rawData );
		if( !xTile || !yTile || xTile==wTiles-1 || yTile==hTiles-1)
		   fill_border(tile,6, x0, y0);
		dcb_hid(tile,buffer,buffer2,x0,y0);
		copy_to_buffer(buffer, tile);
        for (int i=iterations; i>0;i--) {
			dcb_hid2(tile,x0,y0);
			dcb_hid2(tile,x0,y0);
			dcb_hid2(tile,x0,y0);
			dcb_map(tile,x0,y0);
			dcb_correction(tile,x0,y0);
        }
        dcb_color(tile,x0,y0);
        dcb_pp(tile,x0,y0);
        dcb_map(tile,x0,y0);
        dcb_correction2(tile,x0,y0);
        dcb_map(tile,x0,y0);
        dcb_correction(tile,x0,y0);
        dcb_color(tile,x0,y0);
        dcb_map(tile,x0,y0);
        dcb_correction(tile,x0,y0);
        dcb_map(tile,x0,y0);
        dcb_correction(tile,x0,y0);
        dcb_map(tile,x0,y0);
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

#ifdef _OPENMP
        if(omp_get_thread_num()==0)
#endif
        {
    		if( plistener && double(tilesDone)/numTiles > currentProgress){
    			currentProgress+=0.1; // Show progress each 10%
    			plistener->setProgress (currentProgress);
    		}
        }
#pragma omp atomic
        tilesDone++;
    }

#ifdef _OPENMP
	for(int i=0; i<nthreads; i++){
		free(image[i]);
		free(image2[i]);
		free(image3[i]);
		free(chroma[i]);
	}
#endif
	free(image);
    free(image2);
    free(image3);
    free(chroma);

    if(plistener) plistener->setProgress (1.0);
}
#undef TILEBORDER
#undef TILESIZE
#undef CACHESIZE
} /* namespace */