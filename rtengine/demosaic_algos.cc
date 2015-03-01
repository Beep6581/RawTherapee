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
#include <cassert>

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
#include "slicer.h"
#include "rt_math.h"
#include "color.h"
#include "../rtgui/multilangmgr.h"
#include "procparams.h"
#include "sleef.c"
#include "opthelper.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace rtengine {
#undef ABS
#undef DIST

#define ABS(a) ((a)<0?-(a):(a))
#define DIST(a,b) (ABS(a-b))
#define CLIREF(x) LIM(x,-200000.0f,200000.0f) // avoid overflow : do not act directly on image[] or pix[]
#define x1125(a) (a + xdivf(a, 3))
#define x0875(a) (a - xdivf(a, 3))
#define x0250(a) xdivf(a, 2)
#define x00625(a) xdivf(a, 4)
#define x0125(a) xdivf(a, 3)


#define PIX_SORT(a,b) { if ((a)>(b)) {temp=(a);(a)=(b);(b)=temp;} }
#define PIX_SORTV(av,bv)  tempv = _mm_min_ps(av,bv); bv = _mm_max_ps(av,bv); av = tempv; 
extern const Settings* settings;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void RawImageSource::eahd_demosaic () {
  if (plistener) {
    plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::eahd]));
    plistener->setProgress (0.0);
  }

  // prepare cache and constants for cielab conversion
  //TODO: revisit after conversion to D50 illuminant
  lc00 = (0.412453 * imatrices.rgb_cam[0][0] + 0.357580 * imatrices.rgb_cam[0][1] + 0.180423 * imatrices.rgb_cam[0][2]) ;// / 0.950456;
  lc01 = (0.412453 * imatrices.rgb_cam[1][0] + 0.357580 * imatrices.rgb_cam[1][1] + 0.180423 * imatrices.rgb_cam[1][2]) ;// / 0.950456;
  lc02 = (0.412453 * imatrices.rgb_cam[2][0] + 0.357580 * imatrices.rgb_cam[2][1] + 0.180423 * imatrices.rgb_cam[2][2]) ;// / 0.950456;

  lc10 = 0.212671 * imatrices.rgb_cam[0][0] + 0.715160 * imatrices.rgb_cam[0][1] + 0.072169 * imatrices.rgb_cam[0][2];
  lc11 = 0.212671 * imatrices.rgb_cam[1][0] + 0.715160 * imatrices.rgb_cam[1][1] + 0.072169 * imatrices.rgb_cam[1][2];
  lc12 = 0.212671 * imatrices.rgb_cam[2][0] + 0.715160 * imatrices.rgb_cam[2][1] + 0.072169 * imatrices.rgb_cam[2][2];

  lc20 = (0.019334 * imatrices.rgb_cam[0][0] + 0.119193 * imatrices.rgb_cam[0][1] + 0.950227 * imatrices.rgb_cam[0][2]) ;// / 1.088754;
  lc21 = (0.019334 * imatrices.rgb_cam[1][0] + 0.119193 * imatrices.rgb_cam[1][1] + 0.950227 * imatrices.rgb_cam[1][2]) ;// / 1.088754;
  lc22 = (0.019334 * imatrices.rgb_cam[2][0] + 0.119193 * imatrices.rgb_cam[2][1] + 0.950227 * imatrices.rgb_cam[2][2]) ;// / 1.088754;

  int maxindex = 3*65536;//2*65536 3 = avoid crash 3/2013 J.Desmis
  cache = new double[maxindex];
  threshold = (int)(0.008856*MAXVALD);
  for (int i=0; i<maxindex; i++)
    cache[i] = exp(1.0/3.0 * log(double(i) / MAXVALD));

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
    int sh, sv, idx;
    for (int j=1; j<W-1; j++) {
      int dmi = 0;
      for (int x=-1; x<=1; x++) {
        idx = (i+x)%3;
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
      int eL = min(max(dLmaph[3],dLmaph[5]),max(dLmapv[1],dLmapv[7]));
      int eCa = min(max(dCamaph[3],dCamaph[5]),max(dCamapv[1],dCamapv[7]));
      int eCb = min(max(dCbmaph[3],dCbmaph[5]),max(dCbmapv[1],dCbmapv[7]));

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
  float* temp = new float[max(W,H)];
  float* avg = new float[max(W,H)];
  float* dev = new float[max(W,H)];

  memset (temp, 0, max(W,H)*sizeof(float));
  memset (avg, 0, max(W,H)*sizeof(float));
  memset (dev, 0, max(W,H)*sizeof(float));

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
  float* temp = new float[max(W,H)];
  float* avg = new float[max(W,H)];
  float* dev = new float[max(W,H)];

  memset (temp, 0, max(W,H)*sizeof(float));
  memset (avg, 0, max(W,H)*sizeof(float));
  memset (dev, 0, max(W,H)*sizeof(float));

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
#ifdef _OPENMP
#pragma omp parallel for
#endif
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
    plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::hphd]));
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
  freeArray<float>(hpmap, H);//TODO: seems to cause sigabrt ???  why???

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
#define fc(row,col) (prefilters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)
typedef unsigned short ushort;
void RawImageSource::vng4_demosaic () {
	const signed short int *cp, terms[] = {
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
	},
	chood[] = { -1,-1, -1,0, -1,+1, 0,+1, +1,+1, +1,0, +1,-1, 0,-1 };

	double progress = 0.0;
	const bool plistenerActive = plistener;

	if (plistenerActive) {
		plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::vng4]));
		plistener->setProgress (progress);
	}

	const unsigned prefilters = ri->prefilters;
	const int width=W, height=H;
	const int colors=4;
	float (*image)[4];

	image = (float (*)[4]) calloc (height*width, sizeof *image);
	
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int ii=0; ii<H; ii++)
		for (int jj=0; jj<W; jj++)
			image[ii*W+jj][fc(ii,jj)] = rawData[ii][jj];

	{
	int lcode[16][16][32];
	float mul[16][16][8];
	float csum[16][16][4];

// first linear interpolation
	for (int row=0; row < 16; row++)
		for (int col=0; col < 16; col++) {
			int * ip = lcode[row][col];
			int mulcount=0;
			float sum[4];
			memset (sum, 0, sizeof sum);
			for (int y=-1; y <= 1; y++)
				for (int x=-1; x <= 1; x++) {
					int shift = (y==0) + (x==0);
					if (shift == 2)
						continue;
					int color = fc(row+y,col+x);
					*ip++ = (width*y + x)*4 + color;
					
					mul[row][col][mulcount] = (1 << shift);
					*ip++ = color;
					sum[color] += (1 << shift);
					mulcount++;
				}
			int colcount = 0;
			for (int c=0; c < colors; c++)
				if (c != fc(row,col)) {
					*ip++ = c;
					csum[row][col][colcount] = sum[c];
					colcount ++;
				}
		}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int row=1; row < height-1; row++) {
		for (int col=1; col < width-1; col++) {
			float * pix = image[row*width+col];
			int * ip = lcode[row & 15][col & 15];
			float sum[4];
			memset (sum, 0, sizeof sum);
			for (int i=0; i<8;i++, ip+=2)
				sum[ip[1]] += pix[ip[0]] * mul[row & 15][col & 15][i];
			for (int i=0; i<colors-1;i++, ip++)
				pix[ip[0]] = sum[ip[0]] / csum[row & 15][col & 15][i];
		}
	}
  }

	const int prow=7, pcol=1;
	int *code[8][2];
	int t, g;
	int * ip = (int *) calloc ((prow+1)*(pcol+1), 1280);
	for (int row=0; row <= prow; row++)		/* Precalculate for VNG */
		for (int col=0; col <= pcol; col++) {
			code[row][col] = ip;
			for (cp=terms, t=0; t < 64; t++) {
				int y1 = *cp++;
				int x1 = *cp++;
				int y2 = *cp++;
				int x2 = *cp++;
				int weight = *cp++;
				int grads = *cp++;
				int color = fc(row+y1,col+x1);
				if (fc(row+y2,col+x2) != color)
					continue;
				int diag = (fc(row,col+1) == color && fc(row+1,col) == color) ? 2:1;
				if (abs(y1-y2) == diag && abs(x1-x2) == diag)
					continue;
				*ip++ = (y1*width + x1)*4 + color;
				*ip++ = (y2*width + x2)*4 + color;
				*ip++ = weight;
				for (g=0; g < 8; g++)
					if (grads & (1<<g))
						*ip++ = g;
				*ip++ = -1;
			}
			*ip++ = INT_MAX;
			for (cp=chood, g=0; g < 8; g++) {
				int y = *cp++;
				int x = *cp++;
				*ip++ = (y*width + x) * 4;
				int color = fc(row,col);
				if (fc(row+y,col+x) != color && fc(row+y*2,col+x*2) == color)
					*ip++ = (y*width + x) * 8 + color;
				else
					*ip++ = 0;
			}
		}

	if(plistenerActive) {
		progress = 0.1;
		plistener->setProgress (progress);
	}


#ifdef _OPENMP
#pragma omp parallel
#endif
{
	float gval[8], thold, sum[3];
	int g;
	const int progressStep = 64;
	const double progressInc = (0.98-progress)/((height-2)/progressStep);
#ifdef _OPENMP
#pragma omp for nowait
#endif
	for (int row=2; row < height-2; row++) {		/* Do VNG interpolation */
		for (int col=2; col < width-2; col++) {
			float * pix = image[row*width+col];
			int * ip = code[row & prow][col & pcol];
			memset (gval, 0, sizeof gval);
			while ((g = ip[0]) != INT_MAX) {		/* Calculate gradients */
				float diff = fabsf(pix[g] - pix[ip[1]]) * (1<< ip[2]);
				gval[ip[3]] += diff;
				ip+=4;
				while ((g = *ip++) != -1)
					gval[g] += diff;
			}
			ip++;
			{
			float gmin, gmax;
			gmin = gmax = gval[0];			/* Choose a threshold */
			for (g=1; g < 8; g++) {
				if (gmin > gval[g])
					gmin = gval[g];
				if (gmax < gval[g])
					gmax = gval[g];
			}
			thold = gmin + (gmax /2);
			}
			memset (sum, 0, sizeof sum);
			float t1,t2;
			int color = fc(row,col);
			t1 = t2 = pix[color];
			if(color&1) {
				int num = 0;
				for (g=0; g < 8; g++,ip+=2) {		/* Average the neighbors */
					if (gval[g] <= thold) {
						if(ip[1])
							sum[0] += (t1 + pix[ip[1]])*0.5f;
						sum[1] += pix[ip[0] + (color ^ 2)];
						num++;
					}
				}
				t1 += (sum[1] - sum[0]) / num;
			} else {
				int num = 0;
				for (g=0; g < 8; g++,ip+=2) {		/* Average the neighbors */
					if (gval[g] <= thold) {
						sum[1] += pix[ip[0] + 1];
						sum[2] += pix[ip[0] + 3];
						if(ip[1])
							sum[0] += (t1 + pix[ip[1]])*0.5f;
						num++;
					}
				}
				t1 += (sum[1] - sum[0]) / num;
				t2 += (sum[2] - sum[0]) / num;
			}
			green[row][col] = 0.5f * (t1 + t2);
		}
	if(plistenerActive) {
		if((row % progressStep) == 0)
#ifdef _OPENMP
#pragma omp critical (updateprogress)
#endif
{
			progress += progressInc;
			plistener->setProgress (progress);
}
	}
  }

}
  free (code[0][0]);
  free (image);
	if(plistenerActive)
		plistener->setProgress (0.98);

  // Interpolate R and B
#ifdef _OPENMP
#pragma omp parallel for
#endif
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
	if(plistenerActive)
		plistener->setProgress (1.0);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef fc
#define fc(row,col) \
	(ri->get_filters() >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)

#define FORCC for (int c=0; c < colors; c++)

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

  if (plistener) {
    // looks like ppg isn't supported anymore
    //plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::ppg]));
    plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), "xxx"));
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
      pix[0][1] = ULIM(static_cast<float>(guess[i] >> 2), pix[d][1], pix[-d][1]);
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

  red(W,H);
  for (int i=0; i<H; i++)
	  for (int j=0; j<W; j++) {
        red[i][j] = image[i*W+j][0];
  }
  green(W,H);
  for (int i=0; i<H; i++)
	  for (int j=0; j<W; j++) {
        green[i][j] = image[i*W+j][1];
  }
  blue(W,H);
  for (int i=0; i<H; i++)
	  for (int j=0; j<W; j++) {
        blue[i][j] = image[i*W+j][2];
  }
  free (image);
}

void RawImageSource::border_interpolate(unsigned int border, float (*image)[4], unsigned int start, unsigned int end)
{
	unsigned row, col, y, x, f, c, sum[8];
	unsigned int width=W, height=H;
	unsigned int colors = 3;

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

void RawImageSource::border_interpolate2( int winw, int winh, int lborders)
{
int bord=lborders;
int width=winw;
int height=winh;
	for (int i=0; i<height; i++) {

        float sum[6];

		for (int j=0; j<bord; j++) {//first few columns
			for (int c=0; c<6; c++) sum[c]=0;
			for (int i1=i-1; i1<i+2; i1++)
				for (int j1=j-1; j1<j+2; j1++) {
					if ((i1 > -1) && (i1 < height) && (j1 > -1)) {
						int c = FC(i1,j1);
						sum[c] += rawData[i1][j1];
						sum[c+3]++;
					}
				}
			int c=FC(i,j);
			if (c==1) {
				red[i][j]=sum[0]/sum[3];
				green[i][j]=rawData[i][j];
				blue[i][j]=sum[2]/sum[5];
			} else {
				green[i][j]=sum[1]/sum[4];
				if (c==0) {
					red[i][j]=rawData[i][j];
					blue[i][j]=sum[2]/sum[5];
				} else {
					red[i][j]=sum[0]/sum[3];
					blue[i][j]=rawData[i][j];
				}
			}
		}//j

		for (int j=width-bord; j<width; j++) {//last few columns
			for (int c=0; c<6; c++) sum[c]=0;
			for (int i1=i-1; i1<i+2; i1++)
				for (int j1=j-1; j1<j+2; j1++) {
					if ((i1 > -1) && (i1 < height ) && (j1 < width)) {
						int c = FC(i1,j1);
						sum[c] += rawData[i1][j1];
						sum[c+3]++;
					}
				}
			int c=FC(i,j);
			if (c==1) {
				red[i][j]=sum[0]/sum[3];
				green[i][j]=rawData[i][j];
				blue[i][j]=sum[2]/sum[5];
			} else {
				green[i][j]=sum[1]/sum[4];
				if (c==0) {
					red[i][j]=rawData[i][j];
					blue[i][j]=sum[2]/sum[5];
				} else {
					red[i][j]=sum[0]/sum[3];
					blue[i][j]=rawData[i][j];
				}
			}
		}//j
	}//i
	for (int i=0; i<bord; i++) {

        float sum[6];

		for (int j=bord; j<width-bord; j++) {//first few rows
			for (int c=0; c<6; c++) sum[c]=0;
			for (int i1=i-1; i1<i+2; i1++)
				for (int j1=j-1; j1<j+2; j1++) {
					if ((i1 > -1) && (i1 < height) && (j1 > -1)) {
						int c = FC(i1,j1);
						sum[c] += rawData[i1][j1];
						sum[c+3]++;
					}
				}
			int c=FC(i,j);
			if (c==1) {
				red[i][j]=sum[0]/sum[3];
				green[i][j]=rawData[i][j];
				blue[i][j]=sum[2]/sum[5];
			} else {
				green[i][j]=sum[1]/sum[4];
				if (c==0) {
					red[i][j]=rawData[i][j];
					blue[i][j]=sum[2]/sum[5];
				} else {
					red[i][j]=sum[0]/sum[3];
					blue[i][j]=rawData[i][j];
				}
			}
		}//j
	}

	for (int i=height-bord; i<height; i++) {

        float sum[6];

		for (int j=bord; j<width-bord; j++) {//last few rows
			for (int c=0; c<6; c++) sum[c]=0;
			for (int i1=i-1; i1<i+2; i1++)
				for (int j1=j-1; j1<j+2; j1++) {
					if ((i1 > -1) && (i1 < height) && (j1 < width)) {
						int c = FC(i1,j1);
						sum[c] += rawData[i1][j1];
						sum[c+3]++;
					}
				}
			int c=FC(i,j);
			if (c==1) {
				red[i][j]=sum[0]/sum[3];
				green[i][j]=rawData[i][j];
				blue[i][j]=sum[2]/sum[5];
			} else {
				green[i][j]=sum[1]/sum[4];
				if (c==0) {
					red[i][j]=rawData[i][j];
					blue[i][j]=sum[2]/sum[5];
				} else {
					red[i][j]=sum[0]/sum[3];
					blue[i][j]=rawData[i][j];
				}
			}
		}//j
	}

}

// Joint Demosaicing and Denoising using High Order Interpolation Techniques
// Revision 0.9.1a - 09/02/2010 - Contact info: luis.sanz.rodriguez@gmail.com
// Copyright Luis Sanz Rodriguez 2010
// Adapted to RT by Jacques Desmis 3/2013

void RawImageSource::jdl_interpolate_omp()  // from "Lassus"
{
	int width=W, height=H;
	int row,col,c,d,i,u=width,v=2*u,w=3*u,x=4*u,y=5*u,z=6*u,indx,(*dif)[2],(*chr)[2];
	float f[4],g[4];
	float (*image)[4];
	image = (float (*)[4]) calloc (width*height, sizeof *image);
	dif = (int (*)[2]) calloc(width*height, sizeof *dif);
	chr = (int (*)[2]) calloc(width*height, sizeof *chr);
	if (plistener) {
		// this function seems to be unused
		//plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::jdl]));
		plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), "xxx"));
		plistener->setProgress (0.0);
	}

#ifdef _OPENMP
#pragma omp parallel default(none) shared(image,width,height,u,w,v,y,x,z,dif,chr) private(row,col,f,g,indx,c,d,i)
#endif
{
#ifdef _OPENMP
#pragma omp for
#endif
	for (int ii=0; ii<height; ii++)
		for (int jj=0; jj<width; jj++)
			image[ii*width+jj][fc(ii,jj)] = rawData[ii][jj];

	border_interpolate(6, image);

#ifdef _OPENMP
#pragma omp for
#endif
	for (row=5; row < height-5; row++)
		for (col=5+(FC(row,1)&1),indx=row*width+col,c=FC(row,col); col < u-5; col+=2,indx+=2) {
			f[0]=1.f+abs(image[indx-u][1]-image[indx-w][1])+abs(image[indx-u][1]-image[indx+u][1])+abs(image[indx][c]-image[indx-v][c])+abs(image[indx-v][c]-image[indx-x][c]);
			f[1]=1.f+abs(image[indx+1][1]-image[indx+3][1])+abs(image[indx+1][1]-image[indx-1][1])+abs(image[indx][c]-image[indx+2][c])+abs(image[indx+2][c]-image[indx+4][c]);
			f[2]=1.f+abs(image[indx-1][1]-image[indx-3][1])+abs(image[indx-1][1]-image[indx+1][1])+abs(image[indx][c]-image[indx-2][c])+abs(image[indx-2][c]-image[indx-4][c]);
			f[3]=1.f+abs(image[indx+u][1]-image[indx+w][1])+abs(image[indx+u][1]-image[indx-u][1])+abs(image[indx][c]-image[indx+v][c])+abs(image[indx+v][c]-image[indx+x][c]);
			g[0]=CLIP((22.f*image[indx-u][1]+22.f*image[indx-w][1]+2.f*image[indx-y][1]+2.f*image[indx+u][1]+40.f*image[indx][c]-32.f*image[indx-v][c]-8.f*image[indx-x][c])/48.f);
			g[1]=CLIP((22.f*image[indx+1][1]+22.f*image[indx+3][1]+2.f*image[indx+5][1]+2.f*image[indx-1][1]+40.f*image[indx][c]-32.f*image[indx+2][c]-8.f*image[indx+4][c])/48.f);
			g[2]=CLIP((22.f*image[indx-1][1]+22.f*image[indx-3][1]+2.f*image[indx-5][1]+2.f*image[indx+1][1]+40.f*image[indx][c]-32.f*image[indx-2][c]-8.f*image[indx-4][c])/48.f);
			g[3]=CLIP((22.f*image[indx+u][1]+22.f*image[indx+w][1]+2.f*image[indx+y][1]+2.f*image[indx-u][1]+40.f*image[indx][c]-32.f*image[indx+v][c]-8.f*image[indx+x][c])/48.f);
			dif[indx][0]=CLIP((f[3]*g[0]+f[0]*g[3])/(f[0]+f[3]))-image[indx][c];
			dif[indx][1]=CLIP((f[2]*g[1]+f[1]*g[2])/(f[1]+f[2]))-image[indx][c];
		}
#ifdef _OPENMP
#pragma omp for
#endif
	for (row=6; row < height-6; row++)
		for (col=6+(FC(row,2)&1),indx=row*width+col,c=FC(row,col)/2; col < u-6; col+=2,indx+=2) {
			f[0]=1.f+78.f*SQR((float)dif[indx][0])+69.f*(SQR((float) dif[indx-v][0])+SQR((float)dif[indx+v][0]))+51.f*(SQR((float)dif[indx-x][0])+SQR((float)dif[indx+x][0]))+21.f*(SQR((float)dif[indx-z][0])+SQR((float)dif[indx+z][0]))-6.f*SQR((float)dif[indx-v][0]+dif[indx][0]+dif[indx+v][0])-10.f*(SQR((float)dif[indx-x][0]+dif[indx-v][0]+dif[indx][0])+SQR((float)dif[indx][0]+dif[indx+v][0]+dif[indx+x][0]))-7.f*(SQR((float)dif[indx-z][0]+dif[indx-x][0]+dif[indx-v][0])+SQR((float)dif[indx+v][0]+dif[indx+x][0]+dif[indx+z][0]));
			f[1]=1.f+78.f*SQR((float)dif[indx][1])+69.f*(SQR((float)dif[indx-2][1])+SQR((float)dif[indx+2][1]))+51.f*(SQR((float)dif[indx-4][1])+SQR((float)dif[indx+4][1]))+21.f*(SQR((float)dif[indx-6][1])+SQR((float)dif[indx+6][1]))-6.f*SQR((float)dif[indx-2][1]+dif[indx][1]+dif[indx+2][1])-10.f*(SQR((float)dif[indx-4][1]+dif[indx-2][1]+dif[indx][1])+SQR((float)dif[indx][1]+dif[indx+2][1]+dif[indx+4][1]))-7.f*(SQR((float)dif[indx-6][1]+dif[indx-4][1]+dif[indx-2][1])+SQR((float)dif[indx+2][1]+dif[indx+4][1]+dif[indx+6][1]));
			g[0]=ULIM<float>(0.725f*dif[indx][0]+0.1375f*dif[indx-v][0]+0.1375f*dif[indx+v][0],dif[indx-v][0],dif[indx+v][0]);
			g[1]=ULIM<float>(0.725f*dif[indx][1]+0.1375f*dif[indx-2][1]+0.1375f*dif[indx+2][1],dif[indx-2][1],dif[indx+2][1]);
			chr[indx][c]=(f[1]*g[0]+f[0]*g[1])/(f[0]+f[1]);
		}
#ifdef _OPENMP
#pragma omp for
#endif
	for (row=6; row<height-6; row++)
		for (col=6+(FC(row,2)&1),indx=row*width+col,c=1-FC(row,col)/2,d=2*c; col<u-6; col+=2,indx+=2) {
			f[0]=1.f/(float)(1.f+fabs((float)chr[indx-u-1][c]-chr[indx+u+1][c])+fabs((float)chr[indx-u-1][c]-chr[indx-w-3][c])+fabs((float)chr[indx+u+1][c]-chr[indx-w-3][c]));
			f[1]=1.f/(float)(1.f+fabs((float)chr[indx-u+1][c]-chr[indx+u-1][c])+fabs((float)chr[indx-u+1][c]-chr[indx-w+3][c])+fabs((float)chr[indx+u-1][c]-chr[indx-w+3][c]));
			f[2]=1.f/(float)(1.f+fabs((float)chr[indx+u-1][c]-chr[indx-u+1][c])+fabs((float)chr[indx+u-1][c]-chr[indx+w+3][c])+fabs((float)chr[indx-u+1][c]-chr[indx+w-3][c]));
			f[3]=1.f/(float)(1.f+fabs((float)chr[indx+u+1][c]-chr[indx-u-1][c])+fabs((float)chr[indx+u+1][c]-chr[indx+w-3][c])+fabs((float)chr[indx-u-1][c]-chr[indx+w+3][c]));
			g[0]=ULIM(chr[indx-u-1][c],chr[indx-w-1][c],chr[indx-u-3][c]);
			g[1]=ULIM(chr[indx-u+1][c],chr[indx-w+1][c],chr[indx-u+3][c]);
			g[2]=ULIM(chr[indx+u-1][c],chr[indx+w-1][c],chr[indx+u-3][c]);
			g[3]=ULIM(chr[indx+u+1][c],chr[indx+w+1][c],chr[indx+u+3][c]);
			chr[indx][c]=(f[0]*g[0]+f[1]*g[1]+f[2]*g[2]+f[3]*g[3])/(f[0]+f[1]+f[2]+f[3]);
			image[indx][1]=CLIP(image[indx][2-d]+chr[indx][1-c]);
			image[indx][d]=CLIP(image[indx][1]-chr[indx][c]);
		}
#ifdef _OPENMP
#pragma omp for
#endif
	for (row=6; row<height-6; row++)
		for (col=6+(FC(row,1)&1),indx=row*width+col,c=FC(row,col+1)/2,d=2*c; col<u-6; col+=2,indx+=2)
			for(i=0;i<=1;c=1-c,d=2*c,i++){
				f[0]=1.f/(float)(1.f+fabs((float)chr[indx-u][c]-chr[indx+u][c])+fabs((float)chr[indx-u][c]-chr[indx-w][c])+fabs((float)chr[indx+u][c]-chr[indx-w][c]));
				f[1]=1.f/(float)(1.f+fabs((float)chr[indx+1][c]-chr[indx-1][c])+fabs((float)chr[indx+1][c]-chr[indx+3][c])+fabs((float)chr[indx-1][c]-chr[indx+3][c]));
				f[2]=1.f/(float)(1.f+fabs((float)chr[indx-1][c]-chr[indx+1][c])+fabs((float)chr[indx-1][c]-chr[indx-3][c])+fabs((float)chr[indx+1][c]-chr[indx-3][c]));
				f[3]=1.f/(float)(1.f+fabs((float)chr[indx+u][c]-chr[indx-u][c])+fabs((float)chr[indx+u][c]-chr[indx+w][c])+fabs((float)chr[indx-u][c]-chr[indx+w][c]));
				g[0]=0.875f*chr[indx-u][c]+0.125f*chr[indx-w][c];
				g[1]=0.875f*chr[indx+1][c]+0.125f*chr[indx+3][c];
				g[2]=0.875f*chr[indx-1][c]+0.125f*chr[indx-3][c];
				g[3]=0.875f*chr[indx+u][c]+0.125f*chr[indx+w][c];
				image[indx][d]=CLIP(image[indx][1]-(f[0]*g[0]+f[1]*g[1]+f[2]*g[2]+f[3]*g[3])/(f[0]+f[1]+f[2]+f[3]));
			}
#ifdef _OPENMP
#pragma omp for
#endif
	for (int ii=0; ii<height; ii++) {
		for (int jj=0; jj<width; jj++){
			red[ii][jj]	 = CLIP(image[ii*width+jj][0]);
			green[ii][jj]	 = CLIP(image[ii*width+jj][1]);
			blue[ii][jj] 	 = CLIP(image[ii*width+jj][2]);
		}
	}
} // End of parallelization
	free (image);
	free(dif);	free(chr);
	//RawImageSource::refinement_lassus();
}

// LSMME demosaicing algorithm
// L. Zhang and X. Wu,
// Color demozaicing via directional Linear Minimum Mean Square-error Estimation,
// IEEE Trans. on Image Processing, vol. 14, pp. 2167-2178,
// Dec. 2005.
// Adapted to RT by Jacques Desmis 3/2013
// Improved speed and reduced memory consumption by Ingo Weyrich 2/2015
//TODO Tiles to reduce memory consumption
SSEFUNCTION void RawImageSource::lmmse_interpolate_omp(int winw, int winh, int iterations)
{
	const int width=winw, height=winh;
	const int ba = 10;
	const int rr1 = height + 2*ba;
	const int cc1 = width + 2*ba;
	const int w1 = cc1;
	const int w2 = 2*w1;
	const int w3 = 3*w1;
	const int w4 = 4*w1;
	float h0, h1, h2, h3, h4, hs;
	h0 = 1.0f;
	h1 = exp( -1.0f/8.0f);
	h2 = exp( -4.0f/8.0f);
	h3 = exp( -9.0f/8.0f);
	h4 = exp(-16.0f/8.0f);
	hs = h0 + 2.0f*(h1 + h2 + h3 + h4);
	h0 /= hs;
	h1 /= hs;
	h2 /= hs;
	h3 /= hs;
	h4 /= hs;
	int passref;
	int iter;
	if(iterations <=4) {
		iter = iterations-1;
		passref=0;
	} else if (iterations <=6) {
		iter=3;
		passref=iterations-4;
	} else if (iterations <=8) {
		iter=3;
		passref=iterations-6;
	}
	bool applyGamma=true;
	if(iterations==0) {
		applyGamma=false;
		iter=0;
	} else
		applyGamma=true;

	LUTf *gamtab;
	if(applyGamma)
		gamtab = &(Color::gammatab_24_17a);
	else {
		gamtab = new LUTf(65536, LUT_CLIP_ABOVE | LUT_CLIP_BELOW);
		for(int i=0;i<65536;i++)
			(*gamtab)[i] = (float)i / 65535.f;
	}

	if (plistener) {
		plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::lmmse]));
		plistener->setProgress (0.0);
	}

	float *rix[5];
	float *qix[5];
	float *buffer = (float *)calloc(rr1*cc1*5*sizeof(float),1);
	if(buffer == NULL) { // allocation of big block of memory failed, try to get 6 smaller ones
		printf("lmmse_interpolate_omp: allocation of big memory block failed, try to get 5 smaller ones now...\n");
		for(int i=0;i<5;i++)
			qix[i] = (float *)calloc(rr1*cc1*sizeof(float),1);
	} else {
		qix[0] = buffer;
		for(int i=1;i<5;i++)
			qix[i] = qix[i-1] + rr1*cc1;
	}

#ifdef _OPENMP
#pragma omp parallel private(rix)
#endif
{
#ifdef _OPENMP
#pragma omp for
#endif
	for (int rrr=ba; rrr < rr1-ba; rrr++) {
		for (int ccc=ba, row=rrr-ba; ccc < cc1-ba; ccc++) {
			int col = ccc - ba;
			float *rix = qix[4] + rrr*cc1 + ccc;
			rix[0] = (*gamtab)[rawData[row][col]];
		}
	}

#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.1);
}

	// G-R(B)
#ifdef _OPENMP
#pragma omp for schedule(dynamic,16)
#endif
	for (int rr=2; rr < rr1-2; rr++) {
		// G-R(B) at R(B) location
		for (int cc=2+(FC(rr,2)&1); cc < cc1-2; cc+=2) {
			rix[4] = qix[4] + rr*cc1 + cc;
			float v0 = x00625(rix[4][-w1-1]+rix[4][-w1+1]+rix[4][w1-1]+rix[4][w1+1]) +x0250(rix[4][0]);
			// horizontal
			rix[0] = qix[0] + rr*cc1 + cc;
			rix[0][0] = - x0250(rix[4][ -2] + rix[4][ 2])+ xdiv2f(rix[4][ -1] + rix[4][0] + rix[4][ 1]);
			float Y = v0 + xdiv2f(rix[0][0]);
			if (rix[4][0] > 1.75f*Y)
				rix[0][0] = ULIM(rix[0][0],rix[4][ -1],rix[4][ 1]);
			else
				rix[0][0] = LIM(rix[0][0],0.0f,1.0f);
			rix[0][0] -= rix[4][0];
			// vertical
			rix[1] = qix[1] + rr*cc1 + cc;
			rix[1][0] = -x0250(rix[4][-w2] + rix[4][w2])+ xdiv2f(rix[4][-w1] + rix[4][0] + rix[4][w1]);
			Y = v0 + xdiv2f(rix[1][0]);
			if (rix[4][0] > 1.75f*Y)
				rix[1][0] = ULIM(rix[1][0],rix[4][-w1],rix[4][w1]);
			else
				rix[1][0] = LIM(rix[1][0],0.0f,1.0f);
			rix[1][0] -= rix[4][0];
		}
		// G-R(B) at G location
		for (int ccc=2+(FC(rr,3)&1); ccc < cc1-2; ccc+=2) {
			rix[0] = qix[0] + rr*cc1 + ccc;
			rix[1] = qix[1] + rr*cc1 + ccc;
			rix[4] = qix[4] + rr*cc1 + ccc;
			rix[0][0] = x0250(rix[4][ -2] + rix[4][ 2])- xdiv2f(rix[4][ -1] + rix[4][0] + rix[4][ 1]);
			rix[1][0] = x0250(rix[4][-w2] + rix[4][w2])- xdiv2f(rix[4][-w1] + rix[4][0] + rix[4][w1]);
			rix[0][0] = LIM(rix[0][0],-1.0f,0.0f) + rix[4][0];
			rix[1][0] = LIM(rix[1][0],-1.0f,0.0f) + rix[4][0];
		}
	}
#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.2);
}


	// apply low pass filter on differential colors
#ifdef _OPENMP
#pragma omp  for
#endif
	for (int rr=4; rr < rr1-4; rr++)
		for (int cc=4; cc < cc1-4; cc++) {
			rix[0] = qix[0] + rr*cc1 + cc;
			rix[2] = qix[2] + rr*cc1 + cc;
			rix[2][0] = h0*rix[0][0] + h1*(rix[0][ -1] + rix[0][ 1]) + h2*(rix[0][ -2] + rix[0][ 2]) + h3*(rix[0][ -3] + rix[0][ 3]) + h4*(rix[0][ -4] + rix[0][ 4]);
			rix[1] = qix[1] + rr*cc1 + cc;
			rix[3] = qix[3] + rr*cc1 + cc;
			rix[3][0] = h0*rix[1][0] + h1*(rix[1][-w1] + rix[1][w1]) + h2*(rix[1][-w2] + rix[1][w2]) + h3*(rix[1][-w3] + rix[1][w3]) + h4*(rix[1][-w4] + rix[1][w4]);
		}
#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.3);
}

	// interpolate G-R(B) at R(B)
#ifdef _OPENMP
#pragma omp for
#endif
	for (int rr=4; rr < rr1-4; rr++) {
		int cc=4+(FC(rr,4)&1);
#ifdef __SSE2__
		__m128 p1v, p2v, p3v, p4v, p5v, p6v, p7v, p8v, p9v, muv, vxv, vnv, xhv, vhv, xvv, vvv;
		__m128 epsv = _mm_set1_ps(1e-7);
		__m128 ninev = _mm_set1_ps(9.f);
		for (; cc < cc1-10; cc+=8) {
			rix[0] = qix[0] + rr*cc1 + cc;
			rix[1] = qix[1] + rr*cc1 + cc;
			rix[2] = qix[2] + rr*cc1 + cc;
			rix[3] = qix[3] + rr*cc1 + cc;
			rix[4] = qix[4] + rr*cc1 + cc;
			// horizontal
			p1v = LC2VFU(rix[2][-4]);
			p2v = LC2VFU(rix[2][-3]);
			p3v = LC2VFU(rix[2][-2]);
			p4v = LC2VFU(rix[2][-1]);
			p5v = LC2VFU(rix[2][ 0]);
			p6v = LC2VFU(rix[2][ 1]);
			p7v = LC2VFU(rix[2][ 2]);
			p8v = LC2VFU(rix[2][ 3]);
			p9v = LC2VFU(rix[2][ 4]);
			muv = (p1v + p2v + p3v + p4v + p5v + p6v + p7v + p8v + p9v) / ninev;
			vxv = epsv+SQRV(p1v-muv)+SQRV(p2v-muv)+SQRV(p3v-muv)+SQRV(p4v-muv)+SQRV(p5v-muv)+SQRV(p6v-muv)+SQRV(p7v-muv)+SQRV(p8v-muv)+SQRV(p9v-muv);
			p1v -= LC2VFU(rix[0][-4]);
			p2v -= LC2VFU(rix[0][-3]);
			p3v -= LC2VFU(rix[0][-2]);
			p4v -= LC2VFU(rix[0][-1]);
			p5v -= LC2VFU(rix[0][ 0]);
			p6v -= LC2VFU(rix[0][ 1]);
			p7v -= LC2VFU(rix[0][ 2]);
			p8v -= LC2VFU(rix[0][ 3]);
			p9v -= LC2VFU(rix[0][ 4]);
			vnv = epsv+SQRV(p1v)+SQRV(p2v)+SQRV(p3v)+SQRV(p4v)+SQRV(p5v)+SQRV(p6v)+SQRV(p7v)+SQRV(p8v)+SQRV(p9v);
			xhv = (LC2VFU(rix[0][0])*vxv + LC2VFU(rix[2][0])*vnv)/(vxv + vnv);
			vhv = vxv*vnv/(vxv + vnv);

			// vertical
			p1v = LC2VFU(rix[3][-w4]);
			p2v = LC2VFU(rix[3][-w3]);
			p3v = LC2VFU(rix[3][-w2]);
			p4v = LC2VFU(rix[3][-w1]);
			p5v = LC2VFU(rix[3][  0]);
			p6v = LC2VFU(rix[3][ w1]);
			p7v = LC2VFU(rix[3][ w2]);
			p8v = LC2VFU(rix[3][ w3]);
			p9v = LC2VFU(rix[3][ w4]);
			muv = (p1v + p2v + p3v + p4v + p5v + p6v + p7v + p8v + p9v) / ninev;
			vxv = epsv+SQRV(p1v-muv)+SQRV(p2v-muv)+SQRV(p3v-muv)+SQRV(p4v-muv)+SQRV(p5v-muv)+SQRV(p6v-muv)+SQRV(p7v-muv)+SQRV(p8v-muv)+SQRV(p9v-muv);
			p1v -= LC2VFU(rix[1][-w4]);
			p2v -= LC2VFU(rix[1][-w3]);
			p3v -= LC2VFU(rix[1][-w2]);
			p4v -= LC2VFU(rix[1][-w1]);
			p5v -= LC2VFU(rix[1][  0]);
			p6v -= LC2VFU(rix[1][ w1]);
			p7v -= LC2VFU(rix[1][ w2]);
			p8v -= LC2VFU(rix[1][ w3]);
			p9v -= LC2VFU(rix[1][ w4]);
			vnv = epsv+SQRV(p1v)+SQRV(p2v)+SQRV(p3v)+SQRV(p4v)+SQRV(p5v)+SQRV(p6v)+SQRV(p7v)+SQRV(p8v)+SQRV(p9v);
			xvv = (LC2VFU(rix[1][0])*vxv + LC2VFU(rix[3][0])*vnv)/(vxv + vnv);
			vvv = vxv*vnv/(vxv + vnv);
			// interpolated G-R(B)
			muv = (xhv*vvv + xvv*vhv)/(vhv + vvv);
			STC2VFU(rix[4][0], muv);
		}
#endif
		for (; cc < cc1-4; cc+=2) {
			rix[0] = qix[0] + rr*cc1 + cc;
			rix[1] = qix[1] + rr*cc1 + cc;
			rix[2] = qix[2] + rr*cc1 + cc;
			rix[3] = qix[3] + rr*cc1 + cc;
			rix[4] = qix[4] + rr*cc1 + cc;
			// horizontal
			float p1 = rix[2][-4];
			float p2 = rix[2][-3];
			float p3 = rix[2][-2];
			float p4 = rix[2][-1];
			float p5 = rix[2][ 0];
			float p6 = rix[2][ 1];
			float p7 = rix[2][ 2];
			float p8 = rix[2][ 3];
			float p9 = rix[2][ 4];
			float mu = (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) / 9.f;
			float vx = 1e-7+SQR(p1-mu)+SQR(p2-mu)+SQR(p3-mu)+SQR(p4-mu)+SQR(p5-mu)+SQR(p6-mu)+SQR(p7-mu)+SQR(p8-mu)+SQR(p9-mu);
			p1 -= rix[0][-4];
			p2 -= rix[0][-3];
			p3 -= rix[0][-2];
			p4 -= rix[0][-1];
			p5 -= rix[0][ 0];
			p6 -= rix[0][ 1];
			p7 -= rix[0][ 2];
			p8 -= rix[0][ 3];
			p9 -= rix[0][ 4];
			float vn = 1e-7+SQR(p1)+SQR(p2)+SQR(p3)+SQR(p4)+SQR(p5)+SQR(p6)+SQR(p7)+SQR(p8)+SQR(p9);
			float xh = (rix[0][0]*vx + rix[2][0]*vn)/(vx + vn);
			float vh = vx*vn/(vx + vn);

			// vertical
			p1 = rix[3][-w4];
			p2 = rix[3][-w3];
			p3 = rix[3][-w2];
			p4 = rix[3][-w1];
			p5 = rix[3][  0];
			p6 = rix[3][ w1];
			p7 = rix[3][ w2];
			p8 = rix[3][ w3];
			p9 = rix[3][ w4];
			mu = (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) / 9.f;
			vx = 1e-7+SQR(p1-mu)+SQR(p2-mu)+SQR(p3-mu)+SQR(p4-mu)+SQR(p5-mu)+SQR(p6-mu)+SQR(p7-mu)+SQR(p8-mu)+SQR(p9-mu);
			p1 -= rix[1][-w4];
			p2 -= rix[1][-w3];
			p3 -= rix[1][-w2];
			p4 -= rix[1][-w1];
			p5 -= rix[1][  0];
			p6 -= rix[1][ w1];
			p7 -= rix[1][ w2];
			p8 -= rix[1][ w3];
			p9 -= rix[1][ w4];
			vn = 1e-7+SQR(p1)+SQR(p2)+SQR(p3)+SQR(p4)+SQR(p5)+SQR(p6)+SQR(p7)+SQR(p8)+SQR(p9);
			float xv = (rix[1][0]*vx + rix[3][0]*vn)/(vx + vn);
			float vv = vx*vn/(vx + vn);
			// interpolated G-R(B)
			rix[4][0] = (xh*vv + xv*vh)/(vh + vv);
		}
	}

#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.4);
}

	// copy CFA values
#ifdef _OPENMP
#pragma omp for
#endif
	for (int rr=0; rr < rr1; rr++)
		for (int cc=0, row=rr-ba; cc < cc1; cc++) {
			int col=cc-ba;
			int c = FC(rr,cc);
			rix[c] = qix[c] + rr*cc1 + cc;
			if ((row >= 0) & (row < height) & (col >= 0) & (col < width)) {
				rix[c][0] = (*gamtab)[rawData[row][col]];
			}
			else
				rix[c][0] = 0.f;
			if (c != 1) {
				rix[1] = qix[1] + rr*cc1 + cc;
				rix[4] = qix[4] + rr*cc1 + cc;
				rix[1][0] = rix[c][0] + rix[4][0];
			}
		}

#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.5);
}

	// bilinear interpolation for R/B
	// interpolate R/B at G location
#ifdef _OPENMP
#pragma omp for
#endif
	for (int rr=1; rr < rr1-1; rr++)
		for (int cc=1+(FC(rr,2)&1), c=FC(rr,cc+1); cc < cc1-1; cc+=2) {
			rix[c] = qix[c] + rr*cc1 + cc;
			rix[1] = qix[1] + rr*cc1 + cc;
			rix[c][0] = rix[1][0] + xdiv2f(rix[c][ -1] - rix[1][ -1] + rix[c][ 1] - rix[1][ 1]);
			c = 2 - c;
			rix[c] = qix[c] + rr*cc1 + cc;
			rix[c][0] = rix[1][0]+ xdiv2f(rix[c][-w1] - rix[1][-w1] + rix[c][w1] - rix[1][w1]);
			c = 2 - c;
		}

#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.6);
}

	// interpolate R/B at B/R location
#ifdef _OPENMP
#pragma omp  for
#endif
	for (int rr=1; rr < rr1-1; rr++)
		for (int cc=1+(FC(rr,1)&1), c=2-FC(rr,cc); cc < cc1-1; cc+=2) {
			rix[c] = qix[c] + rr*cc1 + cc;
			rix[1] = qix[1] + rr*cc1 + cc;
			rix[c][0] = rix[1][0]+ x0250(rix[c][-w1] - rix[1][-w1] + rix[c][ -1] - rix[1][ -1]+ rix[c][  1] - rix[1][  1] + rix[c][ w1] - rix[1][ w1]);
		}

#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.7);
}

}// End of parallelization 1

	// median filter/
	for (int pass=0; pass < iter; pass++) {
		// Apply 3x3 median filter
		// Compute median(R-G) and median(B-G)

#ifdef _OPENMP
#pragma omp parallel for private(rix)
#endif
		for (int rr=1; rr < rr1-1; rr++) {
			for (int c=0; c < 3; c+=2) {
				int d = c + 3 - (c == 0 ? 0 : 1);
				int cc=1;
#ifdef __SSE2__
				__m128 p1v,p2v,p3v,p4v,p5v,p6v,p7v,p8v,p9v,tempv;
				for (; cc < cc1-4; cc+=4) {
					rix[d] = qix[d] + rr*cc1 + cc;
					rix[c] = qix[c] + rr*cc1 + cc;
					rix[1] = qix[1] + rr*cc1 + cc;
					// Assign 3x3 differential color values
					p1v = LVFU(rix[c][-w1-1])-LVFU(rix[1][-w1-1]); p2v = LVFU(rix[c][-w1])-LVFU(rix[1][-w1]); p3v = LVFU(rix[c][-w1+1])-LVFU(rix[1][-w1+1]);
					p4v = LVFU(rix[c][   -1])-LVFU(rix[1][   -1]); p5v = LVFU(rix[c][  0])-LVFU(rix[1][  0]); p6v = LVFU(rix[c][    1])-LVFU(rix[1][    1]);
					p7v = LVFU(rix[c][ w1-1])-LVFU(rix[1][ w1-1]); p8v = LVFU(rix[c][ w1])-LVFU(rix[1][ w1]); p9v = LVFU(rix[c][ w1+1])-LVFU(rix[1][ w1+1]);
					// Sort for median of 9 values
					PIX_SORTV(p2v,p3v); PIX_SORTV(p5v,p6v); PIX_SORTV(p8v,p9v);
					PIX_SORTV(p1v,p2v); PIX_SORTV(p4v,p5v); PIX_SORTV(p7v,p8v);
					PIX_SORTV(p2v,p3v); PIX_SORTV(p5v,p6v); PIX_SORTV(p8v,p9v);
					p4v = _mm_max_ps(p1v,p4v); p6v = _mm_min_ps(p6v,p9v); PIX_SORTV(p5v,p8v);
					p7v = _mm_max_ps(p4v,p7v); p5v = _mm_max_ps(p5v,p2v); p3v = _mm_min_ps(p3v,p6v);
					p5v = _mm_min_ps(p5v,p8v); PIX_SORTV(p5v,p3v); p5v = _mm_max_ps(p7v,p5v);
					p5v = _mm_min_ps(p3v,p5v);
					_mm_storeu_ps(&rix[d][0], p5v);
				}
#endif
				for (; cc < cc1-1; cc++) {
					float temp;
					rix[d] = qix[d] + rr*cc1 + cc;
					rix[c] = qix[c] + rr*cc1 + cc;
					rix[1] = qix[1] + rr*cc1 + cc;
					// Assign 3x3 differential color values
					float p1 = rix[c][-w1-1]-rix[1][-w1-1]; float p2 = rix[c][-w1]-rix[1][-w1]; float p3 = rix[c][-w1+1]-rix[1][-w1+1];
					float p4 = rix[c][   -1]-rix[1][   -1]; float p5 = rix[c][  0]-rix[1][  0]; float p6 = rix[c][    1]-rix[1][    1];
					float p7 = rix[c][ w1-1]-rix[1][ w1-1]; float p8 = rix[c][ w1]-rix[1][ w1]; float p9 = rix[c][ w1+1]-rix[1][ w1+1];
					// Sort for median of 9 values
					PIX_SORT(p2,p3); PIX_SORT(p5,p6); PIX_SORT(p8,p9);
					PIX_SORT(p1,p2); PIX_SORT(p4,p5); PIX_SORT(p7,p8);
					PIX_SORT(p2,p3); PIX_SORT(p5,p6); PIX_SORT(p8,p9);
					PIX_SORT(p1,p4); PIX_SORT(p6,p9); PIX_SORT(p5,p8);
					PIX_SORT(p4,p7); PIX_SORT(p2,p5); PIX_SORT(p3,p6);
					PIX_SORT(p5,p8); PIX_SORT(p5,p3); PIX_SORT(p7,p5);
					PIX_SORT(p5,p3);
					rix[d][0] = p5;
				}
			}
		}

		// red/blue at GREEN pixel locations & red/blue and green at BLUE/RED pixel locations
#ifdef _OPENMP
#pragma omp parallel for private (rix)
#endif
		for (int rr=0; rr < rr1; rr++) {
			rix[0] = qix[0] + rr*cc1;
			rix[1] = qix[1] + rr*cc1;
			rix[2] = qix[2] + rr*cc1;
			rix[3] = qix[3] + rr*cc1;
			rix[4] = qix[4] + rr*cc1;
			int c0 = FC(rr,0);
			int c1 = FC(rr,1);
			if(c0 == 1){
				c1 = 2 - c1;
				int d = c1 + 3 - (c1 == 0 ? 0 : 1);
				int cc;
				for (cc=0; cc < cc1-1; cc+=2) {
					rix[0][0] = rix[1][0] + rix[3][0];
					rix[2][0] = rix[1][0] + rix[4][0];
					rix[0]++; rix[1]++; rix[2]++; rix[3]++; rix[4]++;
					rix[c1][0] = rix[1][0] + rix[d][0];
					rix[1][0] = 0.5f * (rix[0][0] - rix[3][0] + rix[2][0] - rix[4][0]);
					rix[0]++; rix[1]++; rix[2]++; rix[3]++; rix[4]++;
				}
				if(cc < cc1) { // remaining pixel, only if width is odd
					rix[0][0] = rix[1][0] + rix[3][0];
					rix[2][0] = rix[1][0] + rix[4][0];
				}
			} else {
				c0 = 2 - c0;
				int d = c0 + 3 - (c0 == 0 ? 0 : 1);
				int cc;
				for (cc=0; cc < cc1-1; cc+=2) {
					rix[c0][0] = rix[1][0] + rix[d][0];
					rix[1][0] = 0.5f * (rix[0][0] - rix[3][0] + rix[2][0] - rix[4][0]);
					rix[0]++; rix[1]++; rix[2]++; rix[3]++; rix[4]++;
					rix[0][0] = rix[1][0] + rix[3][0];
					rix[2][0] = rix[1][0] + rix[4][0];
					rix[0]++; rix[1]++; rix[2]++; rix[3]++; rix[4]++;
				}
				if(cc < cc1) { // remaining pixel, only if width is odd
					rix[c0][0] = rix[1][0] + rix[d][0];
					rix[1][0] = 0.5f * (rix[0][0] - rix[3][0] + rix[2][0] - rix[4][0]);
				}
			}
		}
	}

	if (plistener) plistener->setProgress (0.8);

	if(applyGamma)
		gamtab = &(Color::igammatab_24_17);
	else {
		for(int i=0;i<65536;i++)
			(*gamtab)[i] = (float)i + 0.5f;
	}

	array2D<float> (*rgb[3]);
	rgb[0] = &red;
	rgb[1] = &green;
	rgb[2] = &blue;

	// copy result back to image matrix
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int row=0; row < height; row++) {
		for (int col=0, rr=row+ba; col < width; col++) {
			int cc = col+ba;
			int c = FC(row,col);
			for (int ii=0; ii < 3; ii++)
				if (ii != c) {
					float *rix = qix[ii] + rr*cc1 + cc;
					(*(rgb[ii]))[row][col] = (*gamtab)[65535.f*rix[0]];
				} else {
					(*(rgb[ii]))[row][col] = CLIP(rawData[row][col]);
				}
		}
	}

	if (plistener) plistener->setProgress (1.0);
	if(buffer)
		free(buffer);
	else
		for(int i=0;i<5;i++)
			free(qix[i]);

	if(!applyGamma)
		delete gamtab;

	if(iterations > 4 && iterations <=6) refinement(passref);
	else if(iterations > 6) refinement_lassus(passref);

}

/***
*
*   Bayer CFA Demosaicing using Integrated Gaussian Vector on Color Differences
*   Revision 1.0 - 2013/02/28
*
*   Copyright (c) 2007-2013 Luis Sanz Rodriguez
*   Using High Order Interpolation technique by Jim S, Jimmy Li, and Sharmil Randhawa
*
*   Contact info: luis.sanz.rodriguez@gmail.com
*
*   This code is distributed under a GNU General Public License, version 3.
*   Visit <http://www.gnu.org/licenses/> for more information.
*
***/
// Adapted to RT by Jacques Desmis 3/2013
// SSE version by Ingo Weyrich 5/2013
#ifdef __SSE2__
#define CLIPV(a) LIMV(a,zerov,c65535v)
SSEFUNCTION void RawImageSource::igv_interpolate(int winw, int winh)
{
	static const float eps=1e-5f, epssq=1e-5f;//mod epssq -10f =>-5f Jacques 3/2013 to prevent artifact (divide by zero)

	static const int h1=1, h2=2, h3=3, h5=5;
	const int width=winw, height=winh;
	const int v1=1*width, v2=2*width, v3=3*width, v5=5*width;
	float* rgb[2];
	float* chr[4];
	float *rgbarray, *vdif, *hdif, *chrarray;
	rgbarray	= (float (*)) malloc((width*height) * sizeof( float ) );
	rgb[0] = rgbarray;
	rgb[1] = rgbarray + (width*height)/2;

	vdif  = (float (*))	calloc( width*height/2, sizeof *vdif );
	hdif  = (float (*))	calloc( width*height/2, sizeof *hdif );

	chrarray	= (float (*)) calloc( width*height, sizeof( float ) );
	chr[0] = chrarray;
	chr[1] = chrarray + (width*height)/2;

	// mapped chr[2] and chr[3] to hdif and hdif, because these are out of use, when chr[2] and chr[3] are used
	chr[2] = hdif;
	chr[3] = vdif;

	border_interpolate2(winw,winh,7);

	if (plistener) {
		plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::igv]));
		plistener->setProgress (0.0);
	}
#ifdef _OPENMP
#pragma omp parallel default(none) shared(rgb,vdif,hdif,chr)
#endif
{
	__m128 ngv, egv, wgv, sgv, nvv, evv, wvv, svv, nwgv, negv, swgv, segv, nwvv, nevv, swvv, sevv,tempv,temp1v,temp2v,temp3v,temp4v,temp5v,temp6v,temp7v,temp8v;
	__m128 epsv = _mm_set1_ps( eps );
	__m128 epssqv = _mm_set1_ps( epssq );
	__m128 c65535v = _mm_set1_ps( 65535.f );
	__m128 c23v = _mm_set1_ps( 23.f );
	__m128 c40v = _mm_set1_ps( 40.f );
	__m128 c51v = _mm_set1_ps( 51.f );
	__m128 c32v = _mm_set1_ps( 32.f );
	__m128 c8v = _mm_set1_ps( 8.f );
	__m128 c7v = _mm_set1_ps( 7.f );
	__m128 c6v = _mm_set1_ps( 6.f );
	__m128 c10v = _mm_set1_ps( 10.f );
	__m128 c21v = _mm_set1_ps( 21.f );
	__m128 c78v = _mm_set1_ps( 78.f );
	__m128 c69v = _mm_set1_ps( 69.f );
	__m128 c3145680v = _mm_set1_ps( 3145680.f );
	__m128 onev = _mm_set1_ps ( 1.f );
	__m128 zerov = _mm_set1_ps ( 0.f );
	__m128 d725v = _mm_set1_ps ( 0.725f );
	__m128 d1375v = _mm_set1_ps ( 0.1375f );

	float *dest1, *dest2;
	float ng, eg, wg, sg, nv, ev, wv, sv, nwg, neg, swg, seg, nwv, nev, swv, sev;
#ifdef _OPENMP
#pragma omp for
#endif
	for (int row=0; row<height-0; row++) {
		dest1 = rgb[FC(row,0)&1];
		dest2 = rgb[FC(row,1)&1];
		int col, indx;
		for (col=0, indx=row*width+col; col<width-7; col+=8, indx+=8) {
			temp1v = LVFU( rawData[row][col] );
			temp1v = CLIPV( temp1v );
			temp2v = LVFU( rawData[row][col+4] );
			temp2v = CLIPV( temp2v );
			tempv = _mm_shuffle_ps( temp1v, temp2v, _MM_SHUFFLE( 2,0,2,0 ) );
			_mm_storeu_ps( &dest1[indx>>1], tempv );
			tempv = _mm_shuffle_ps( temp1v, temp2v, _MM_SHUFFLE( 3,1,3,1 ) );
			_mm_storeu_ps( &dest2[indx>>1], tempv );
		}
		for (; col<width; col++, indx+=2) {
			dest1[indx>>1] = CLIP(rawData[row][col]);	//rawData = RT datas
			col++;
			dest2[indx>>1] = CLIP(rawData[row][col]);	//rawData = RT datas
		}
	}
#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.13);
}

#ifdef _OPENMP
#pragma omp for
#endif
	for (int row=5; row<height-5; row++) {
		int col,indx,indx1;
		for (col=5+(FC(row,1)&1),indx=row*width+col, indx1=indx>>1; col<width-12; col+=8, indx+=8, indx1+=4) {
			//N,E,W,S Gradients
			ngv=(epsv+(vabsf(LVFU(rgb[1][(indx-v1)>>1])-LVFU(rgb[1][(indx-v3)>>1]))+vabsf(LVFU(rgb[0][indx1])-LVFU(rgb[0][(indx1-v1)])))/c65535v);
			egv=(epsv+(vabsf(LVFU(rgb[1][(indx+h1)>>1])-LVFU(rgb[1][(indx+h3)>>1]))+vabsf(LVFU(rgb[0][indx1])-LVFU(rgb[0][(indx1+h1)])))/c65535v);
			wgv=(epsv+(vabsf(LVFU(rgb[1][(indx-h1)>>1])-LVFU(rgb[1][(indx-h3)>>1]))+vabsf(LVFU(rgb[0][indx1])-LVFU(rgb[0][(indx1-h1)])))/c65535v);
			sgv=(epsv+(vabsf(LVFU(rgb[1][(indx+v1)>>1])-LVFU(rgb[1][(indx+v3)>>1]))+vabsf(LVFU(rgb[0][indx1])-LVFU(rgb[0][(indx1+v1)])))/c65535v);
			//N,E,W,S High Order Interpolation (Li & Randhawa)
			//N,E,W,S Hamilton Adams Interpolation
			// (48.f * 65535.f) = 3145680.f
			tempv = c40v*LVFU(rgb[0][indx1]);
			nvv=LIMV(((c23v*LVFU(rgb[1][(indx-v1)>>1])+c23v*LVFU(rgb[1][(indx-v3)>>1])+LVFU(rgb[1][(indx-v5)>>1])+LVFU(rgb[1][(indx+v1)>>1])+tempv-c32v*LVFU(rgb[0][(indx1-v1)])-c8v*LVFU(rgb[0][(indx1-v2)])))/c3145680v, zerov, onev);
			evv=LIMV(((c23v*LVFU(rgb[1][(indx+h1)>>1])+c23v*LVFU(rgb[1][(indx+h3)>>1])+LVFU(rgb[1][(indx+h5)>>1])+LVFU(rgb[1][(indx-h1)>>1])+tempv-c32v*LVFU(rgb[0][(indx1+h1)])-c8v*LVFU(rgb[0][(indx1+h2)])))/c3145680v, zerov, onev);
			wvv=LIMV(((c23v*LVFU(rgb[1][(indx-h1)>>1])+c23v*LVFU(rgb[1][(indx-h3)>>1])+LVFU(rgb[1][(indx-h5)>>1])+LVFU(rgb[1][(indx+h1)>>1])+tempv-c32v*LVFU(rgb[0][(indx1-h1)])-c8v*LVFU(rgb[0][(indx1-h2)])))/c3145680v, zerov, onev);
			svv=LIMV(((c23v*LVFU(rgb[1][(indx+v1)>>1])+c23v*LVFU(rgb[1][(indx+v3)>>1])+LVFU(rgb[1][(indx+v5)>>1])+LVFU(rgb[1][(indx-v1)>>1])+tempv-c32v*LVFU(rgb[0][(indx1+v1)])-c8v*LVFU(rgb[0][(indx1+v2)])))/c3145680v, zerov, onev);
			//Horizontal and vertical color differences
			tempv = LVFU( rgb[0][indx1] ) / c65535v;
			_mm_storeu_ps( &vdif[indx1], (sgv*nvv+ngv*svv)/(ngv+sgv)- tempv );
			_mm_storeu_ps( &hdif[indx1], (wgv*evv+egv*wvv)/(egv+wgv)- tempv );
		}
		// borders without SSE
		for (; col<width-5; col+=2, indx+=2, indx1++) {
			//N,E,W,S Gradients
			ng=(eps+(fabsf(rgb[1][(indx-v1)>>1]-rgb[1][(indx-v3)>>1])+fabsf(rgb[0][indx1]-rgb[0][(indx1-v1)]))/65535.f);;
			eg=(eps+(fabsf(rgb[1][(indx+h1)>>1]-rgb[1][(indx+h3)>>1])+fabsf(rgb[0][indx1]-rgb[0][(indx1+h1)]))/65535.f);
			wg=(eps+(fabsf(rgb[1][(indx-h1)>>1]-rgb[1][(indx-h3)>>1])+fabsf(rgb[0][indx1]-rgb[0][(indx1-h1)]))/65535.f);
			sg=(eps+(fabsf(rgb[1][(indx+v1)>>1]-rgb[1][(indx+v3)>>1])+fabsf(rgb[0][indx1]-rgb[0][(indx1+v1)]))/65535.f);
			//N,E,W,S High Order Interpolation (Li & Randhawa)
			//N,E,W,S Hamilton Adams Interpolation
			// (48.f * 65535.f) = 3145680.f
			nv=LIM(((23.0f*rgb[1][(indx-v1)>>1]+23.0f*rgb[1][(indx-v3)>>1]+rgb[1][(indx-v5)>>1]+rgb[1][(indx+v1)>>1]+40.0f*rgb[0][indx1]-32.0f*rgb[0][(indx1-v1)]-8.0f*rgb[0][(indx1-v2)]))/3145680.f, 0.0f, 1.0f);
			ev=LIM(((23.0f*rgb[1][(indx+h1)>>1]+23.0f*rgb[1][(indx+h3)>>1]+rgb[1][(indx+h5)>>1]+rgb[1][(indx-h1)>>1]+40.0f*rgb[0][indx1]-32.0f*rgb[0][(indx1+h1)]-8.0f*rgb[0][(indx1+h2)]))/3145680.f, 0.0f, 1.0f);
			wv=LIM(((23.0f*rgb[1][(indx-h1)>>1]+23.0f*rgb[1][(indx-h3)>>1]+rgb[1][(indx-h5)>>1]+rgb[1][(indx+h1)>>1]+40.0f*rgb[0][indx1]-32.0f*rgb[0][(indx1-h1)]-8.0f*rgb[0][(indx1-h2)]))/3145680.f, 0.0f, 1.0f);
			sv=LIM(((23.0f*rgb[1][(indx+v1)>>1]+23.0f*rgb[1][(indx+v3)>>1]+rgb[1][(indx+v5)>>1]+rgb[1][(indx-v1)>>1]+40.0f*rgb[0][indx1]-32.0f*rgb[0][(indx1+v1)]-8.0f*rgb[0][(indx1+v2)]))/3145680.f, 0.0f, 1.0f);
			//Horizontal and vertical color differences
			vdif[indx1]=(sg*nv+ng*sv)/(ng+sg)-(rgb[0][indx1])/65535.f;
			hdif[indx1]=(wg*ev+eg*wv)/(eg+wg)-(rgb[0][indx1])/65535.f;
		}
	}
#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.26);
}
#ifdef _OPENMP
#pragma omp for
#endif
	for (int row=7; row<height-7; row++) {
		int col,d,indx1;
		for (col=7+(FC(row,1)&1), indx1=(row*width+col)>>1, d=FC(row,col)/2; col<width-14; col+=8, indx1+=4) {
			//H&V integrated gaussian vector over variance on color differences
			//Mod Jacques 3/2013
			ngv=LIMV(epssqv+c78v*SQRV(LVFU(vdif[indx1]))+c69v*(SQRV(LVFU(vdif[indx1-v1]))+SQRV(LVFU(vdif[indx1+v1])))+c51v*(SQRV(LVFU(vdif[indx1-v2]))+SQRV(LVFU(vdif[indx1+v2])))+c21v*(SQRV(LVFU(vdif[indx1-v3]))+SQRV(LVFU(vdif[indx1+v3])))-c6v*SQRV(LVFU(vdif[indx1-v1])+LVFU(vdif[indx1])+LVFU(vdif[indx1+v1]))
			  -c10v*(SQRV(LVFU(vdif[indx1-v2])+LVFU(vdif[indx1-v1])+LVFU(vdif[indx1]))+SQRV(LVFU(vdif[indx1])+LVFU(vdif[indx1+v1])+LVFU(vdif[indx1+v2])))-c7v*(SQRV(LVFU(vdif[indx1-v3])+LVFU(vdif[indx1-v2])+LVFU(vdif[indx1-v1]))+SQRV(LVFU(vdif[indx1+v1])+LVFU(vdif[indx1+v2])+LVFU(vdif[indx1+v3]))),zerov,onev);
			egv=LIMV(epssqv+c78v*SQRV(LVFU(hdif[indx1]))+c69v*(SQRV(LVFU(hdif[indx1-h1]))+SQRV(LVFU(hdif[indx1+h1])))+c51v*(SQRV(LVFU(hdif[indx1-h2]))+SQRV(LVFU(hdif[indx1+h2])))+c21v*(SQRV(LVFU(hdif[indx1-h3]))+SQRV(LVFU(hdif[indx1+h3])))-c6v*SQRV(LVFU(hdif[indx1-h1])+LVFU(hdif[indx1])+LVFU(hdif[indx1+h1]))
			  -c10v*(SQRV(LVFU(hdif[indx1-h2])+LVFU(hdif[indx1-h1])+LVFU(hdif[indx1]))+SQRV(LVFU(hdif[indx1])+LVFU(hdif[indx1+h1])+LVFU(hdif[indx1+h2])))-c7v*(SQRV(LVFU(hdif[indx1-h3])+LVFU(hdif[indx1-h2])+LVFU(hdif[indx1-h1]))+SQRV(LVFU(hdif[indx1+h1])+LVFU(hdif[indx1+h2])+LVFU(hdif[indx1+h3]))),zerov,onev);
			//Limit chrominance using H/V neighbourhood
			nvv=ULIMV(d725v*LVFU(vdif[indx1])+d1375v*LVFU(vdif[indx1-v1])+d1375v*LVFU(vdif[indx1+v1]),LVFU(vdif[indx1-v1]),LVFU(vdif[indx1+v1]));
			evv=ULIMV(d725v*LVFU(hdif[indx1])+d1375v*LVFU(hdif[indx1-h1])+d1375v*LVFU(hdif[indx1+h1]),LVFU(hdif[indx1-h1]),LVFU(hdif[indx1+h1]));
			//Chrominance estimation
			tempv = (egv*nvv+ngv*evv)/(ngv+egv);
			_mm_storeu_ps(&(chr[d][indx1]), tempv);
			//Green channel population
			temp1v = c65535v * tempv + LVFU(rgb[0][indx1]);
			_mm_storeu_ps( &(rgb[0][indx1]), temp1v );
		}
		for (; col<width-7; col+=2, indx1++) {
			//H&V integrated gaussian vector over variance on color differences
			//Mod Jacques 3/2013
			ng=LIM(epssq+78.0f*SQR(vdif[indx1])+69.0f*(SQR(vdif[indx1-v1])+SQR(vdif[indx1+v1]))+51.0f*(SQR(vdif[indx1-v2])+SQR(vdif[indx1+v2]))+21.0f*(SQR(vdif[indx1-v3])+SQR(vdif[indx1+v3]))-6.0f*SQR(vdif[indx1-v1]+vdif[indx1]+vdif[indx1+v1])
			  -10.0f*(SQR(vdif[indx1-v2]+vdif[indx1-v1]+vdif[indx1])+SQR(vdif[indx1]+vdif[indx1+v1]+vdif[indx1+v2]))-7.0f*(SQR(vdif[indx1-v3]+vdif[indx1-v2]+vdif[indx1-v1])+SQR(vdif[indx1+v1]+vdif[indx1+v2]+vdif[indx1+v3])),0.f,1.f);
			eg=LIM(epssq+78.0f*SQR(hdif[indx1])+69.0f*(SQR(hdif[indx1-h1])+SQR(hdif[indx1+h1]))+51.0f*(SQR(hdif[indx1-h2])+SQR(hdif[indx1+h2]))+21.0f*(SQR(hdif[indx1-h3])+SQR(hdif[indx1+h3]))-6.0f*SQR(hdif[indx1-h1]+hdif[indx1]+hdif[indx1+h1])
			  -10.0f*(SQR(hdif[indx1-h2]+hdif[indx1-h1]+hdif[indx1])+SQR(hdif[indx1]+hdif[indx1+h1]+hdif[indx1+h2]))-7.0f*(SQR(hdif[indx1-h3]+hdif[indx1-h2]+hdif[indx1-h1])+SQR(hdif[indx1+h1]+hdif[indx1+h2]+hdif[indx1+h3])),0.f,1.f);
			//Limit chrominance using H/V neighbourhood
			nv=ULIM(0.725f*vdif[indx1]+0.1375f*vdif[indx1-v1]+0.1375f*vdif[indx1+v1],vdif[indx1-v1],vdif[indx1+v1]);
			ev=ULIM(0.725f*hdif[indx1]+0.1375f*hdif[indx1-h1]+0.1375f*hdif[indx1+h1],hdif[indx1-h1],hdif[indx1+h1]);
			//Chrominance estimation
			chr[d][indx1]=(eg*nv+ng*ev)/(ng+eg);
			//Green channel population
			rgb[0][indx1]=rgb[0][indx1]+65535.f*chr[d][indx1];
		}
	}
#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.39);
}
#ifdef _OPENMP
#pragma omp for
#endif
	for (int row=7; row<height-7; row++) {
		int col,indx,c;
		for (col=7+(FC(row,1)&1), indx=row*width+col, c=1-FC(row,col)/2; col<width-14; col+=8, indx+=8) {
			//NW,NE,SW,SE Gradients
			nwgv=onev/(epsv+vabsf(LVFU(chr[c][(indx-v1-h1)>>1])-LVFU(chr[c][(indx-v3-h3)>>1]))+vabsf(LVFU(chr[c][(indx+v1+h1)>>1])-LVFU(chr[c][(indx-v3-h3)>>1])));
			negv=onev/(epsv+vabsf(LVFU(chr[c][(indx-v1+h1)>>1])-LVFU(chr[c][(indx-v3+h3)>>1]))+vabsf(LVFU(chr[c][(indx+v1-h1)>>1])-LVFU(chr[c][(indx-v3+h3)>>1])));
			swgv=onev/(epsv+vabsf(LVFU(chr[c][(indx+v1-h1)>>1])-LVFU(chr[c][(indx+v3+h3)>>1]))+vabsf(LVFU(chr[c][(indx-v1+h1)>>1])-LVFU(chr[c][(indx+v3-h3)>>1])));
			segv=onev/(epsv+vabsf(LVFU(chr[c][(indx+v1+h1)>>1])-LVFU(chr[c][(indx+v3-h3)>>1]))+vabsf(LVFU(chr[c][(indx-v1-h1)>>1])-LVFU(chr[c][(indx+v3+h3)>>1])));
			//Limit NW,NE,SW,SE Color differences
			nwvv=ULIMV(LVFU(chr[c][(indx-v1-h1)>>1]),LVFU(chr[c][(indx-v3-h1)>>1]),LVFU(chr[c][(indx-v1-h3)>>1]));
			nevv=ULIMV(LVFU(chr[c][(indx-v1+h1)>>1]),LVFU(chr[c][(indx-v3+h1)>>1]),LVFU(chr[c][(indx-v1+h3)>>1]));
			swvv=ULIMV(LVFU(chr[c][(indx+v1-h1)>>1]),LVFU(chr[c][(indx+v3-h1)>>1]),LVFU(chr[c][(indx+v1-h3)>>1]));
			sevv=ULIMV(LVFU(chr[c][(indx+v1+h1)>>1]),LVFU(chr[c][(indx+v3+h1)>>1]),LVFU(chr[c][(indx+v1+h3)>>1]));
			//Interpolate chrominance: R@B and B@R
			tempv = (nwgv*nwvv+negv*nevv+swgv*swvv+segv*sevv)/(nwgv+negv+swgv+segv);
			_mm_storeu_ps( &(chr[c][indx>>1]), tempv);
		}
		for (; col<width-7; col+=2, indx+=2) {
			//NW,NE,SW,SE Gradients
			nwg=1.0f/(eps+fabsf(chr[c][(indx-v1-h1)>>1]-chr[c][(indx-v3-h3)>>1])+fabsf(chr[c][(indx+v1+h1)>>1]-chr[c][(indx-v3-h3)>>1]));
			neg=1.0f/(eps+fabsf(chr[c][(indx-v1+h1)>>1]-chr[c][(indx-v3+h3)>>1])+fabsf(chr[c][(indx+v1-h1)>>1]-chr[c][(indx-v3+h3)>>1]));
			swg=1.0f/(eps+fabsf(chr[c][(indx+v1-h1)>>1]-chr[c][(indx+v3+h3)>>1])+fabsf(chr[c][(indx-v1+h1)>>1]-chr[c][(indx+v3-h3)>>1]));
			seg=1.0f/(eps+fabsf(chr[c][(indx+v1+h1)>>1]-chr[c][(indx+v3-h3)>>1])+fabsf(chr[c][(indx-v1-h1)>>1]-chr[c][(indx+v3+h3)>>1]));
			//Limit NW,NE,SW,SE Color differences
			nwv=ULIM(chr[c][(indx-v1-h1)>>1],chr[c][(indx-v3-h1)>>1],chr[c][(indx-v1-h3)>>1]);
			nev=ULIM(chr[c][(indx-v1+h1)>>1],chr[c][(indx-v3+h1)>>1],chr[c][(indx-v1+h3)>>1]);
			swv=ULIM(chr[c][(indx+v1-h1)>>1],chr[c][(indx+v3-h1)>>1],chr[c][(indx+v1-h3)>>1]);
			sev=ULIM(chr[c][(indx+v1+h1)>>1],chr[c][(indx+v3+h1)>>1],chr[c][(indx+v1+h3)>>1]);
			//Interpolate chrominance: R@B and B@R
			chr[c][indx>>1]=(nwg*nwv+neg*nev+swg*swv+seg*sev)/(nwg+neg+swg+seg);
		}
	}
#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.65);
}
#ifdef _OPENMP
#pragma omp for
#endif
	for (int row=7; row<height-7; row++) {
		int col,indx;
		for (col=7+(FC(row,0)&1), indx=row*width+col; col<width-14; col+=8, indx+=8) {
			//N,E,W,S Gradients
			ngv=onev/(epsv+vabsf(LVFU(chr[0][(indx-v1)>>1])-LVFU(chr[0][(indx-v3)>>1]))+vabsf(LVFU(chr[0][(indx+v1)>>1])-LVFU(chr[0][(indx-v3)>>1])));
			egv=onev/(epsv+vabsf(LVFU(chr[0][(indx+h1)>>1])-LVFU(chr[0][(indx+h3)>>1]))+vabsf(LVFU(chr[0][(indx-h1)>>1])-LVFU(chr[0][(indx+h3)>>1])));
			wgv=onev/(epsv+vabsf(LVFU(chr[0][(indx-h1)>>1])-LVFU(chr[0][(indx-h3)>>1]))+vabsf(LVFU(chr[0][(indx+h1)>>1])-LVFU(chr[0][(indx-h3)>>1])));
			sgv=onev/(epsv+vabsf(LVFU(chr[0][(indx+v1)>>1])-LVFU(chr[0][(indx+v3)>>1]))+vabsf(LVFU(chr[0][(indx-v1)>>1])-LVFU(chr[0][(indx+v3)>>1])));
			//Interpolate chrominance: R@G and B@G
			tempv = ((ngv*LVFU(chr[0][(indx-v1)>>1])+egv*LVFU(chr[0][(indx+h1)>>1])+wgv*LVFU(chr[0][(indx-h1)>>1])+sgv*LVFU(chr[0][(indx+v1)>>1]))/(ngv+egv+wgv+sgv));
			_mm_storeu_ps( &chr[0+2][indx>>1], tempv);
			}
		for (; col<width-7; col+=2, indx+=2) {
			//N,E,W,S Gradients
			ng=1.0f/(eps+fabsf(chr[0][(indx-v1)>>1]-chr[0][(indx-v3)>>1])+fabsf(chr[0][(indx+v1)>>1]-chr[0][(indx-v3)>>1]));
			eg=1.0f/(eps+fabsf(chr[0][(indx+h1)>>1]-chr[0][(indx+h3)>>1])+fabsf(chr[0][(indx-h1)>>1]-chr[0][(indx+h3)>>1]));
			wg=1.0f/(eps+fabsf(chr[0][(indx-h1)>>1]-chr[0][(indx-h3)>>1])+fabsf(chr[0][(indx+h1)>>1]-chr[0][(indx-h3)>>1]));
			sg=1.0f/(eps+fabsf(chr[0][(indx+v1)>>1]-chr[0][(indx+v3)>>1])+fabsf(chr[0][(indx-v1)>>1]-chr[0][(indx+v3)>>1]));
			//Interpolate chrominance: R@G and B@G
			chr[0+2][indx>>1]=((ng*chr[0][(indx-v1)>>1]+eg*chr[0][(indx+h1)>>1]+wg*chr[0][(indx-h1)>>1]+sg*chr[0][(indx+v1)>>1])/(ng+eg+wg+sg));
			}
	}
#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.78);
}
#ifdef _OPENMP
#pragma omp for
#endif
	for (int row=7; row<height-7; row++) {
		int col,indx;
		for (col=7+(FC(row,0)&1), indx=row*width+col; col<width-14; col+=8, indx+=8) {
			//N,E,W,S Gradients
			ngv=onev/(epsv+vabsf(LVFU(chr[1][(indx-v1)>>1])-LVFU(chr[1][(indx-v3)>>1]))+vabsf(LVFU(chr[1][(indx+v1)>>1])-LVFU(chr[1][(indx-v3)>>1])));
			egv=onev/(epsv+vabsf(LVFU(chr[1][(indx+h1)>>1])-LVFU(chr[1][(indx+h3)>>1]))+vabsf(LVFU(chr[1][(indx-h1)>>1])-LVFU(chr[1][(indx+h3)>>1])));
			wgv=onev/(epsv+vabsf(LVFU(chr[1][(indx-h1)>>1])-LVFU(chr[1][(indx-h3)>>1]))+vabsf(LVFU(chr[1][(indx+h1)>>1])-LVFU(chr[1][(indx-h3)>>1])));
			sgv=onev/(epsv+vabsf(LVFU(chr[1][(indx+v1)>>1])-LVFU(chr[1][(indx+v3)>>1]))+vabsf(LVFU(chr[1][(indx-v1)>>1])-LVFU(chr[1][(indx+v3)>>1])));
			//Interpolate chrominance: R@G and B@G
			tempv = ((ngv*LVFU(chr[1][(indx-v1)>>1])+egv*LVFU(chr[1][(indx+h1)>>1])+wgv*LVFU(chr[1][(indx-h1)>>1])+sgv*LVFU(chr[1][(indx+v1)>>1]))/(ngv+egv+wgv+sgv));
			_mm_storeu_ps( &chr[1+2][indx>>1], tempv);
			}
		for (; col<width-7; col+=2, indx+=2) {
			//N,E,W,S Gradients
			ng=1.0f/(eps+fabsf(chr[1][(indx-v1)>>1]-chr[1][(indx-v3)>>1])+fabsf(chr[1][(indx+v1)>>1]-chr[1][(indx-v3)>>1]));
			eg=1.0f/(eps+fabsf(chr[1][(indx+h1)>>1]-chr[1][(indx+h3)>>1])+fabsf(chr[1][(indx-h1)>>1]-chr[1][(indx+h3)>>1]));
			wg=1.0f/(eps+fabsf(chr[1][(indx-h1)>>1]-chr[1][(indx-h3)>>1])+fabsf(chr[1][(indx+h1)>>1]-chr[1][(indx-h3)>>1]));
			sg=1.0f/(eps+fabsf(chr[1][(indx+v1)>>1]-chr[1][(indx+v3)>>1])+fabsf(chr[1][(indx-v1)>>1]-chr[1][(indx+v3)>>1]));
			//Interpolate chrominance: R@G and B@G
			chr[1+2][indx>>1]=((ng*chr[1][(indx-v1)>>1]+eg*chr[1][(indx+h1)>>1]+wg*chr[1][(indx-h1)>>1]+sg*chr[1][(indx+v1)>>1])/(ng+eg+wg+sg));
			}
	}
#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.91);
}
	float *src1, *src2, *redsrc0, *redsrc1, *bluesrc0, *bluesrc1;
#ifdef _OPENMP
#pragma omp for
#endif
	for(int row=7; row<height-7; row++){
		int col,indx,fc;
		fc = FC(row,7)&1;
		src1 = rgb[fc];
		src2 = rgb[fc^1];
		redsrc0 = chr[fc<<1];
		redsrc1 = chr[(fc^1)<<1];
		bluesrc0 = chr[(fc<<1)+1];
		bluesrc1 = chr[((fc^1)<<1)+1];
		for(col=7, indx=row*width+col; col<width-14; col+=8, indx+=8) {
			temp1v = LVFU( src1[indx>>1] );
			temp2v = LVFU( src2[(indx+1)>>1] );
			tempv = _mm_shuffle_ps( temp1v, temp2v, _MM_SHUFFLE( 1,0,1,0 ) );
			tempv = _mm_shuffle_ps( tempv, tempv, _MM_SHUFFLE( 3,1,2,0 ) );
			_mm_storeu_ps( &green[row][col], CLIPV( tempv ));
			temp5v = LVFU(redsrc0[indx>>1]);
			temp6v = LVFU(redsrc1[(indx+1)>>1]);
			temp3v = _mm_shuffle_ps( temp5v, temp6v, _MM_SHUFFLE( 1,0,1,0 ) );
			temp3v = _mm_shuffle_ps( temp3v, temp3v, _MM_SHUFFLE( 3,1,2,0 ) );
			temp3v = CLIPV( tempv - c65535v * temp3v );
			_mm_storeu_ps( &red[row][col], temp3v);
			temp7v = LVFU(bluesrc0[indx>>1]);
			temp8v = LVFU(bluesrc1[(indx+1)>>1]);
			temp4v = _mm_shuffle_ps( temp7v, temp8v, _MM_SHUFFLE( 1,0,1,0 ) );
			temp4v = _mm_shuffle_ps( temp4v, temp4v, _MM_SHUFFLE( 3,1,2,0 ) );
			temp4v = CLIPV( tempv - c65535v * temp4v );
			_mm_storeu_ps( &blue[row][col], temp4v);

			tempv = _mm_shuffle_ps( temp1v, temp2v, _MM_SHUFFLE( 3,2,3,2 ) );
			tempv = _mm_shuffle_ps( tempv, tempv, _MM_SHUFFLE( 3,1,2,0 ) );
			_mm_storeu_ps( &green[row][col+4], CLIPV( tempv ));

			temp3v = _mm_shuffle_ps( temp5v, temp6v, _MM_SHUFFLE( 3,2,3,2 ) );
			temp3v = _mm_shuffle_ps( temp3v, temp3v, _MM_SHUFFLE( 3,1,2,0 ) );
			temp3v = CLIPV( tempv - c65535v * temp3v );
			_mm_storeu_ps( &red[row][col+4], temp3v);
			temp4v = _mm_shuffle_ps( temp7v, temp8v, _MM_SHUFFLE( 3,2,3,2 ) );
			temp4v = _mm_shuffle_ps( temp4v, temp4v, _MM_SHUFFLE( 3,1,2,0 ) );
			temp4v = CLIPV( tempv - c65535v * temp4v );
			_mm_storeu_ps( &blue[row][col+4], temp4v);
		}

		for(; col<width-7; col++, indx+=2) {
			red  [row][col] = CLIP(src1[indx>>1]-65535.f*redsrc0[indx>>1]);
			green[row][col] = CLIP(src1[indx>>1]);
			blue [row][col] = CLIP(src1[indx>>1]-65535.f*bluesrc0[indx>>1]);
			col++;
			red  [row][col] = CLIP(src2[(indx+1)>>1]-65535.f*redsrc1[(indx+1)>>1]);
			green[row][col] = CLIP(src2[(indx+1)>>1]);
			blue [row][col] = CLIP(src2[(indx+1)>>1]-65535.f*bluesrc1[(indx+1)>>1]);
		}
	}
}// End of parallelization
	if (plistener) plistener->setProgress (1.0);

	free(chrarray); free(rgbarray);
	free(vdif); free(hdif);
}
#undef CLIPV
#else
void RawImageSource::igv_interpolate(int winw, int winh)
{
	static const float eps=1e-5f, epssq=1e-5f;//mod epssq -10f =>-5f Jacques 3/2013 to prevent artifact (divide by zero)
	static const int h1=1, h2=2, h3=3, h4=4, h5=5, h6=6;
	const int width=winw, height=winh;
	const int v1=1*width, v2=2*width, v3=3*width, v4=4*width, v5=5*width, v6=6*width;
	float* rgb[3];
	float* chr[2];
	float (*rgbarray), *vdif, *hdif, (*chrarray);

	rgbarray	= (float (*)) calloc(width*height*3, sizeof( float));
	rgb[0] = rgbarray;
	rgb[1] = rgbarray + (width*height);
	rgb[2] = rgbarray + 2*(width*height);

	chrarray	= (float (*)) calloc(width*height*2, sizeof( float));
	chr[0] = chrarray;
	chr[1] = chrarray + (width*height);

	vdif  = (float (*))    calloc(width*height/2, sizeof *vdif);
	hdif  = (float (*))    calloc(width*height/2, sizeof *hdif);

	border_interpolate2(winw,winh,7);

	if (plistener) {
		plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::igv]));
		plistener->setProgress (0.0);
	}
#ifdef _OPENMP
#pragma omp parallel default(none) shared(rgb,vdif,hdif,chr)
#endif
{

	float ng, eg, wg, sg, nv, ev, wv, sv, nwg, neg, swg, seg, nwv, nev, swv, sev;

#ifdef _OPENMP
#pragma omp for
#endif
	for (int row=0; row<height-0; row++)
		for (int col=0, indx=row*width+col; col<width-0; col++, indx++) {
			int c=FC(row,col);
			rgb[c][indx]=CLIP(rawData[row][col]);	//rawData = RT datas
		}
//	border_interpolate2(7, rgb);

#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.13);
}

#ifdef _OPENMP
#pragma omp for
#endif
	for (int row=5; row<height-5; row++)
		for (int col=5+(FC(row,1)&1), indx=row*width+col, c=FC(row,col); col<width-5; col+=2, indx+=2) {
			//N,E,W,S Gradients
			ng=(eps+(fabsf(rgb[1][indx-v1]-rgb[1][indx-v3])+fabsf(rgb[c][indx]-rgb[c][indx-v2]))/65535.f);;
			eg=(eps+(fabsf(rgb[1][indx+h1]-rgb[1][indx+h3])+fabsf(rgb[c][indx]-rgb[c][indx+h2]))/65535.f);
			wg=(eps+(fabsf(rgb[1][indx-h1]-rgb[1][indx-h3])+fabsf(rgb[c][indx]-rgb[c][indx-h2]))/65535.f);
			sg=(eps+(fabsf(rgb[1][indx+v1]-rgb[1][indx+v3])+fabsf(rgb[c][indx]-rgb[c][indx+v2]))/65535.f);
			//N,E,W,S High Order Interpolation (Li & Randhawa)
			//N,E,W,S Hamilton Adams Interpolation
			// (48.f * 65535.f) = 3145680.f
			nv=LIM(((23.0f*rgb[1][indx-v1]+23.0f*rgb[1][indx-v3]+rgb[1][indx-v5]+rgb[1][indx+v1]+40.0f*rgb[c][indx]-32.0f*rgb[c][indx-v2]-8.0f*rgb[c][indx-v4]))/3145680.f, 0.0f, 1.0f);
			ev=LIM(((23.0f*rgb[1][indx+h1]+23.0f*rgb[1][indx+h3]+rgb[1][indx+h5]+rgb[1][indx-h1]+40.0f*rgb[c][indx]-32.0f*rgb[c][indx+h2]-8.0f*rgb[c][indx+h4]))/3145680.f, 0.0f, 1.0f);
			wv=LIM(((23.0f*rgb[1][indx-h1]+23.0f*rgb[1][indx-h3]+rgb[1][indx-h5]+rgb[1][indx+h1]+40.0f*rgb[c][indx]-32.0f*rgb[c][indx-h2]-8.0f*rgb[c][indx-h4]))/3145680.f, 0.0f, 1.0f);
			sv=LIM(((23.0f*rgb[1][indx+v1]+23.0f*rgb[1][indx+v3]+rgb[1][indx+v5]+rgb[1][indx-v1]+40.0f*rgb[c][indx]-32.0f*rgb[c][indx+v2]-8.0f*rgb[c][indx+v4]))/3145680.f, 0.0f, 1.0f);
			//Horizontal and vertical color differences
			vdif[indx>>1]=(sg*nv+ng*sv)/(ng+sg)-(rgb[c][indx])/65535.f;
			hdif[indx>>1]=(wg*ev+eg*wv)/(eg+wg)-(rgb[c][indx])/65535.f;
		}

#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.26);
}

#ifdef _OPENMP
#pragma omp for
#endif
	for (int row=7; row<height-7; row++)
		for (int col=7+(FC(row,1)&1), indx=row*width+col, c=FC(row,col), d=c/2; col<width-7; col+=2, indx+=2) {
			//H&V integrated gaussian vector over variance on color differences
			//Mod Jacques 3/2013
			ng=LIM(epssq+78.0f*SQR(vdif[indx>>1])+69.0f*(SQR(vdif[(indx-v2)>>1])+SQR(vdif[(indx+v2)>>1]))+51.0f*(SQR(vdif[(indx-v4)>>1])+SQR(vdif[(indx+v4)>>1]))+21.0f*(SQR(vdif[(indx-v6)>>1])+SQR(vdif[(indx+v6)>>1]))-6.0f*SQR(vdif[(indx-v2)>>1]+vdif[indx>>1]+vdif[(indx+v2)>>1])
			  -10.0f*(SQR(vdif[(indx-v4)>>1]+vdif[(indx-v2)>>1]+vdif[indx>>1])+SQR(vdif[indx>>1]+vdif[(indx+v2)>>1]+vdif[(indx+v4)>>1]))-7.0f*(SQR(vdif[(indx-v6)>>1]+vdif[(indx-v4)>>1]+vdif[(indx-v2)>>1])+SQR(vdif[(indx+v2)>>1]+vdif[(indx+v4)>>1]+vdif[(indx+v6)>>1])),0.f,1.f);
			eg=LIM(epssq+78.0f*SQR(hdif[indx>>1])+69.0f*(SQR(hdif[(indx-h2)>>1])+SQR(hdif[(indx+h2)>>1]))+51.0f*(SQR(hdif[(indx-h4)>>1])+SQR(hdif[(indx+h4)>>1]))+21.0f*(SQR(hdif[(indx-h6)>>1])+SQR(hdif[(indx+h6)>>1]))-6.0f*SQR(hdif[(indx-h2)>>1]+hdif[indx>>1]+hdif[(indx+h2)>>1])
			  -10.0f*(SQR(hdif[(indx-h4)>>1]+hdif[(indx-h2)>>1]+hdif[indx>>1])+SQR(hdif[indx>>1]+hdif[(indx+h2)>>1]+hdif[(indx+h4)>>1]))-7.0f*(SQR(hdif[(indx-h6)>>1]+hdif[(indx-h4)>>1]+hdif[(indx-h2)>>1])+SQR(hdif[(indx+h2)>>1]+hdif[(indx+h4)>>1]+hdif[(indx+h6)>>1])),0.f,1.f);
			//Limit chrominance using H/V neighbourhood
			nv=ULIM(0.725f*vdif[indx>>1]+0.1375f*vdif[(indx-v2)>>1]+0.1375f*vdif[(indx+v2)>>1],vdif[(indx-v2)>>1],vdif[(indx+v2)>>1]);
			ev=ULIM(0.725f*hdif[indx>>1]+0.1375f*hdif[(indx-h2)>>1]+0.1375f*hdif[(indx+h2)>>1],hdif[(indx-h2)>>1],hdif[(indx+h2)>>1]);
			//Chrominance estimation
			chr[d][indx]=(eg*nv+ng*ev)/(ng+eg);
			//Green channel population
			rgb[1][indx]=rgb[c][indx]+65535.f*chr[d][indx];
		}

#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.39);
}

//	free(vdif); free(hdif);
#ifdef _OPENMP
#pragma omp for
#endif
	for (int row=7; row<height-7; row+=2)
		for (int col=7+(FC(row,1)&1), indx=row*width+col, c=1-FC(row,col)/2; col<width-7; col+=2, indx+=2) {
			//NW,NE,SW,SE Gradients
			nwg=1.0f/(eps+fabsf(chr[c][indx-v1-h1]-chr[c][indx-v3-h3])+fabsf(chr[c][indx+v1+h1]-chr[c][indx-v3-h3]));
			neg=1.0f/(eps+fabsf(chr[c][indx-v1+h1]-chr[c][indx-v3+h3])+fabsf(chr[c][indx+v1-h1]-chr[c][indx-v3+h3]));
			swg=1.0f/(eps+fabsf(chr[c][indx+v1-h1]-chr[c][indx+v3+h3])+fabsf(chr[c][indx-v1+h1]-chr[c][indx+v3-h3]));
			seg=1.0f/(eps+fabsf(chr[c][indx+v1+h1]-chr[c][indx+v3-h3])+fabsf(chr[c][indx-v1-h1]-chr[c][indx+v3+h3]));
			//Limit NW,NE,SW,SE Color differences
			nwv=ULIM(chr[c][indx-v1-h1],chr[c][indx-v3-h1],chr[c][indx-v1-h3]);
			nev=ULIM(chr[c][indx-v1+h1],chr[c][indx-v3+h1],chr[c][indx-v1+h3]);
			swv=ULIM(chr[c][indx+v1-h1],chr[c][indx+v3-h1],chr[c][indx+v1-h3]);
			sev=ULIM(chr[c][indx+v1+h1],chr[c][indx+v3+h1],chr[c][indx+v1+h3]);
			//Interpolate chrominance: R@B and B@R
			chr[c][indx]=(nwg*nwv+neg*nev+swg*swv+seg*sev)/(nwg+neg+swg+seg);
		}
#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.52);
}
#ifdef _OPENMP
#pragma omp for
#endif
	for (int row=8; row<height-7; row+=2)
		for (int col=7+(FC(row,1)&1), indx=row*width+col, c=1-FC(row,col)/2; col<width-7; col+=2, indx+=2) {
			//NW,NE,SW,SE Gradients
			nwg=1.0f/(eps+fabsf(chr[c][indx-v1-h1]-chr[c][indx-v3-h3])+fabsf(chr[c][indx+v1+h1]-chr[c][indx-v3-h3]));
			neg=1.0f/(eps+fabsf(chr[c][indx-v1+h1]-chr[c][indx-v3+h3])+fabsf(chr[c][indx+v1-h1]-chr[c][indx-v3+h3]));
			swg=1.0f/(eps+fabsf(chr[c][indx+v1-h1]-chr[c][indx+v3+h3])+fabsf(chr[c][indx-v1+h1]-chr[c][indx+v3-h3]));
			seg=1.0f/(eps+fabsf(chr[c][indx+v1+h1]-chr[c][indx+v3-h3])+fabsf(chr[c][indx-v1-h1]-chr[c][indx+v3+h3]));
			//Limit NW,NE,SW,SE Color differences
			nwv=ULIM(chr[c][indx-v1-h1],chr[c][indx-v3-h1],chr[c][indx-v1-h3]);
			nev=ULIM(chr[c][indx-v1+h1],chr[c][indx-v3+h1],chr[c][indx-v1+h3]);
			swv=ULIM(chr[c][indx+v1-h1],chr[c][indx+v3-h1],chr[c][indx+v1-h3]);
			sev=ULIM(chr[c][indx+v1+h1],chr[c][indx+v3+h1],chr[c][indx+v1+h3]);
			//Interpolate chrominance: R@B and B@R
			chr[c][indx]=(nwg*nwv+neg*nev+swg*swv+seg*sev)/(nwg+neg+swg+seg);
		}
#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.65);
}
#ifdef _OPENMP
#pragma omp for
#endif
	for (int row=7; row<height-7; row++)
		for (int col=7+(FC(row,0)&1), indx=row*width+col; col<width-7; col+=2, indx+=2) {
			//N,E,W,S Gradients
			ng=1.0f/(eps+fabsf(chr[0][indx-v1]-chr[0][indx-v3])+fabsf(chr[0][indx+v1]-chr[0][indx-v3]));
			eg=1.0f/(eps+fabsf(chr[0][indx+h1]-chr[0][indx+h3])+fabsf(chr[0][indx-h1]-chr[0][indx+h3]));
			wg=1.0f/(eps+fabsf(chr[0][indx-h1]-chr[0][indx-h3])+fabsf(chr[0][indx+h1]-chr[0][indx-h3]));
			sg=1.0f/(eps+fabsf(chr[0][indx+v1]-chr[0][indx+v3])+fabsf(chr[0][indx-v1]-chr[0][indx+v3]));
			//Interpolate chrominance: R@G and B@G
			chr[0][indx]=((ng*chr[0][indx-v1]+eg*chr[0][indx+h1]+wg*chr[0][indx-h1]+sg*chr[0][indx+v1])/(ng+eg+wg+sg));
		}
#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.78);
}
#ifdef _OPENMP
#pragma omp for
#endif
	for (int row=7; row<height-7; row++)
		for (int col=7+(FC(row,0)&1), indx=row*width+col; col<width-7; col+=2, indx+=2) {

			//N,E,W,S Gradients
			ng=1.0f/(eps+fabsf(chr[1][indx-v1]-chr[1][indx-v3])+fabsf(chr[1][indx+v1]-chr[1][indx-v3]));
			eg=1.0f/(eps+fabsf(chr[1][indx+h1]-chr[1][indx+h3])+fabsf(chr[1][indx-h1]-chr[1][indx+h3]));
			wg=1.0f/(eps+fabsf(chr[1][indx-h1]-chr[1][indx-h3])+fabsf(chr[1][indx+h1]-chr[1][indx-h3]));
			sg=1.0f/(eps+fabsf(chr[1][indx+v1]-chr[1][indx+v3])+fabsf(chr[1][indx-v1]-chr[1][indx+v3]));
			//Interpolate chrominance: R@G and B@G
			chr[1][indx]=((ng*chr[1][indx-v1]+eg*chr[1][indx+h1]+wg*chr[1][indx-h1]+sg*chr[1][indx+v1])/(ng+eg+wg+sg));
			}
#ifdef _OPENMP
#pragma omp single
#endif
{
	if (plistener) plistener->setProgress (0.91);

	//Interpolate borders
//	border_interpolate2(7, rgb);
}
/*
#ifdef _OPENMP
#pragma omp for
#endif
	for (int row=0; row < height; row++)  //borders
		for (int col=0; col < width; col++) {
			if (col==7 && row >= 7 && row < height-7)
				col = width-7;
			int indxc=row*width+col;
			red  [row][col] = rgb[indxc][0];
			green[row][col] = rgb[indxc][1];
			blue [row][col] = rgb[indxc][2];
		}
*/

#ifdef _OPENMP
#pragma omp for
#endif
	for(int row=7; row<height-7; row++)
		for(int col=7, indx=row*width+col; col<width-7; col++, indx++) {
			red  [row][col] = CLIP(rgb[1][indx]-65535.f*chr[0][indx]);
			green[row][col] = CLIP(rgb[1][indx]);
			blue [row][col] = CLIP(rgb[1][indx]-65535.f*chr[1][indx]);
		}
}// End of parallelization

	if (plistener) plistener->setProgress (1.0);

	free(chrarray); free(rgbarray);
	free(vdif); free(hdif);
}
#endif


/*
   Adaptive Homogeneity-Directed interpolation is based on
   the work of Keigo Hirakawa, Thomas Parks, and Paul Lee.
 */
#define TS 256		/* Tile Size */
#define FORC(cnt) for (c=0; c < cnt; c++)
#define FORC3 FORC(3)

void RawImageSource::ahd_demosaic(int winx, int winy, int winw, int winh)
{
    int i, j, k, top, left, row, col, tr, tc, c, d, val, hm[2];
    float (*pix)[4], (*rix)[3];
    static const int dir[4] = { -1, 1, -TS, TS };
    float ldiff[2][4], abdiff[2][4], leps, abeps;
    float xyz[3], xyz_cam[3][4];
    float (*cbrt);
    float (*rgb)[TS][TS][3];
    float (*lab)[TS][TS][3];
    float (*lix)[3];
    char (*homo)[TS][TS], *buffer;
    double r;

    int width=W, height=H;
    float (*image)[4];
    int colors = 3;

    const double xyz_rgb[3][3] = {        /* XYZ from RGB */
        { 0.412453, 0.357580, 0.180423 },
        { 0.212671, 0.715160, 0.072169 },
        { 0.019334, 0.119193, 0.950227 }
    };

    const float d65_white[3] = { 0.950456, 1, 1.088754 };

    if (plistener) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::ahd]));
        plistener->setProgress (0.0);
    }

    image = (float (*)[4]) calloc (H*W, sizeof *image);
    for (int ii=0; ii<H; ii++)
        for (int jj=0; jj<W; jj++)
            image[ii*W+jj][fc(ii,jj)] = rawData[ii][jj];

	cbrt = (float (*)) calloc (0x10000, sizeof *cbrt);
    for (i=0; i < 0x10000; i++) {
        r = (double)i / 65535.0;
        cbrt[i] = r > 0.008856 ? pow(r,0.333333333) : 7.787*r + 16/116.0;
    }

    for (i=0; i < 3; i++)
        for (j=0; j < colors; j++)
            for (xyz_cam[i][j] = k=0; k < 3; k++)
	            xyz_cam[i][j] += xyz_rgb[i][k] * imatrices.rgb_cam[k][j] / d65_white[i];

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
                    rgb[0][row-top][col-left][1] = ULIM(static_cast<float>(val),pix[-1][1],pix[1][1]);
                    val = 0.25*((pix[-width][1] + pix[0][c] + pix[width][1]) * 2
                          - pix[-2*width][c] - pix[2*width][c]) ;
                    rgb[1][row-top][col-left][1] = ULIM(static_cast<float>(val),pix[-width][1],pix[width][1]);
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
                                  - rix[+TS-1][1] - rix[+TS+1][1]) );
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

            /*  Build homogeneity maps from the CIELab images: */
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
                    leps = min(max(ldiff[0][0],ldiff[0][1]),
                               max(ldiff[1][2],ldiff[1][3]));
                    abeps = min(max(abdiff[0][0],abdiff[0][1]),
                                max(abdiff[1][2],abdiff[1][3]));
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

void RawImageSource::nodemosaic(bool bw)
{
    red(W,H);
    green(W,H);
    blue(W,H);
#pragma omp parallel for
    for (int i=0; i<H; i++) {
        for (int j=0; j<W; j++){
            if (bw) {
                red[i][j] = green[i][j] = blue[i][j] = rawData[i][j];
            } else if(ri->getSensorType()!=ST_FUJI_XTRANS){
                switch( FC(i,j)) {
                case 0: red[i][j] = rawData[i][j]; green[i][j]=blue[i][j]=0; break;
                case 1: green[i][j] = rawData[i][j]; red[i][j]=blue[i][j]=0; break;
                case 2: blue[i][j] = rawData[i][j]; red[i][j]=green[i][j]=0; break;
                }
            } else {
                switch( ri->XTRANSFC(i,j)) {
                case 0: red[i][j] = rawData[i][j]; green[i][j]=blue[i][j]=0; break;
                case 1: green[i][j] = rawData[i][j]; red[i][j]=blue[i][j]=0; break;
                case 2: blue[i][j] = rawData[i][j]; red[i][j]=green[i][j]=0; break;
                }
            }
        }
    }
}

/*
   Refinement based on EECI demosaicing algorithm by L. Chang and Y.P. Tan
   Paul Lee
   Adapted for Rawtherapee - Jacques Desmis 04/2013
*/

#ifdef __SSE2__
#define CLIPV(a) LIMV(a,ZEROV,c65535v)
#endif
SSEFUNCTION void RawImageSource::refinement(int PassCount)
{
	MyTime t1e,t2e;
    t1e.set();

	int width=W;
	int height=H;
	int w1 = width;
	int w2 = 2*w1;
	if (plistener) {
		plistener->setProgressStr (M("TP_RAW_DMETHOD_PROGRESSBAR_REFINE"));
	}

	array2D<float> *rgb[3];
	rgb[0] = &red;
	rgb[1] = &green;
	rgb[2] = &blue;

	for (int b=0; b<PassCount; b++) {
		if (plistener) {
			plistener->setProgress ((float)b/PassCount);
		}


#ifdef _OPENMP
#pragma omp parallel
#endif
{
		float *pix[3];

 /* Reinforce interpolated green pixels on RED/BLUE pixel locations */
#ifdef _OPENMP
#pragma omp for
#endif
		for (int row=2; row < height-2; row++) {
			int col = 2+(FC(row,2) & 1);
			int c = FC(row,col);
#ifdef __SSE2__
			__m128 dLv, dRv, dUv, dDv, v0v;
			__m128 onev = _mm_set1_ps(1.f);
			__m128 zd5v = _mm_set1_ps(0.5f);
			__m128 c65535v = _mm_set1_ps(65535.f);
			for (; col < width-8; col+=8) {
				int indx = row*width+col;
				pix[c] = (float*)(*rgb[c]) + indx;
				pix[1] = (float*)(*rgb[1]) + indx;
				dLv = onev/(onev+vabsf(LC2VFU(pix[c][ -2])-LC2VFU(pix[c][0]))+vabsf(LC2VFU(pix[1][ 1])-LC2VFU(pix[1][ -1])));
				dRv = onev/(onev+vabsf(LC2VFU(pix[c][  2])-LC2VFU(pix[c][0]))+vabsf(LC2VFU(pix[1][ 1])-LC2VFU(pix[1][ -1])));
				dUv = onev/(onev+vabsf(LC2VFU(pix[c][-w2])-LC2VFU(pix[c][0]))+vabsf(LC2VFU(pix[1][w1])-LC2VFU(pix[1][-w1])));
				dDv = onev/(onev+vabsf(LC2VFU(pix[c][ w2])-LC2VFU(pix[c][0]))+vabsf(LC2VFU(pix[1][w1])-LC2VFU(pix[1][-w1])));
				v0v = CLIPV(LC2VFU(pix[c][0]) + zd5v +((LC2VFU(pix[1][-1])-LC2VFU(pix[c][-1]))*dLv +(LC2VFU(pix[1][1])-LC2VFU(pix[c][1]))*dRv +(LC2VFU(pix[1][-w1])-LC2VFU(pix[c][-w1]))*dUv + (LC2VFU(pix[1][w1])-LC2VFU(pix[c][w1]))*dDv ) / (dLv+dRv+dUv+dDv));
				STC2VFU(pix[1][0],v0v);
			}
#endif
			for (; col < width-2; col+=2) {
				int indx = row*width+col;
				pix[c] = (float*)(*rgb[c]) + indx;
				pix[1] = (float*)(*rgb[1]) + indx;
				float dL = 1.f/(1.f+fabsf(pix[c][ -2]-pix[c][0])+fabsf(pix[1][ 1]-pix[1][ -1]));
				float dR = 1.f/(1.f+fabsf(pix[c][  2]-pix[c][0])+fabsf(pix[1][ 1]-pix[1][ -1]));
				float dU = 1.f/(1.f+fabsf(pix[c][-w2]-pix[c][0])+fabsf(pix[1][w1]-pix[1][-w1]));
				float dD = 1.f/(1.f+fabsf(pix[c][ w2]-pix[c][0])+fabsf(pix[1][w1]-pix[1][-w1]));
				float v0 = (pix[c][0] + 0.5f +((pix[1][ -1]-pix[c][ -1])*dL +(pix[1][  1]-pix[c][  1])*dR +(pix[1][-w1]-pix[c][-w1])*dU + (pix[1][ w1]-pix[c][ w1])*dD ) / (dL+dR+dU+dD));
				pix[1][0] = CLIP(v0);
			}
		}
  /* Reinforce interpolated red/blue pixels on GREEN pixel locations */
#ifdef _OPENMP
#pragma omp for
#endif
		for (int row=2; row < height-2; row++) {
			int col = 2+(FC(row,3) & 1);
			int c = FC(row,col+1);
#ifdef __SSE2__
			__m128 dLv, dRv, dUv, dDv, v0v;
			__m128 onev = _mm_set1_ps(1.f);
			__m128 zd5v = _mm_set1_ps(0.5f);
			__m128 c65535v = _mm_set1_ps(65535.f);
			for (; col < width-8; col+=8) {
				int indx = row*width+col;
				pix[1] = (float*)(*rgb[1]) + indx;
				for (int i=0; i < 2; c=2-c, i++) {
					pix[c] = (float*)(*rgb[c]) + indx;
					dLv = onev/(onev+vabsf(LC2VFU(pix[1][ -2])-LC2VFU(pix[1][0]))+vabsf(LC2VFU(pix[c][ 1])-LC2VFU(pix[c][ -1])));
					dRv = onev/(onev+vabsf(LC2VFU(pix[1][  2])-LC2VFU(pix[1][0]))+vabsf(LC2VFU(pix[c][ 1])-LC2VFU(pix[c][ -1])));
					dUv = onev/(onev+vabsf(LC2VFU(pix[1][-w2])-LC2VFU(pix[1][0]))+vabsf(LC2VFU(pix[c][w1])-LC2VFU(pix[c][-w1])));
					dDv = onev/(onev+vabsf(LC2VFU(pix[1][ w2])-LC2VFU(pix[1][0]))+vabsf(LC2VFU(pix[c][w1])-LC2VFU(pix[c][-w1])));
					v0v = CLIPV(LC2VFU(pix[1][0]) + zd5v -((LC2VFU(pix[1][-1])-LC2VFU(pix[c][-1]))*dLv + (LC2VFU(pix[1][1])-LC2VFU(pix[c][1]))*dRv +(LC2VFU(pix[1][-w1])-LC2VFU(pix[c][-w1]))*dUv + (LC2VFU(pix[1][w1])-LC2VFU(pix[c][w1]))*dDv ) / (dLv+dRv+dUv+dDv));
					STC2VFU(pix[c][0],v0v);
				}
			}
#endif
			for (; col < width-2; col+=2) {
				int indx = row*width+col;
				pix[1] = (float*)(*rgb[1]) + indx;
				for (int i=0; i < 2; c=2-c, i++) {
					pix[c] = (float*)(*rgb[c]) + indx;
					float dL = 1.f/(1.f+fabsf(pix[1][ -2]-pix[1][0])+fabsf(pix[c][ 1]-pix[c][ -1]));
					float dR = 1.f/(1.f+fabsf(pix[1][  2]-pix[1][0])+fabsf(pix[c][ 1]-pix[c][ -1]));
					float dU = 1.f/(1.f+fabsf(pix[1][-w2]-pix[1][0])+fabsf(pix[c][w1]-pix[c][-w1]));
					float dD = 1.f/(1.f+fabsf(pix[1][ w2]-pix[1][0])+fabsf(pix[c][w1]-pix[c][-w1]));
					float v0 = (pix[1][0] + 0.5f -((pix[1][ -1]-pix[c][ -1])*dL + (pix[1][  1]-pix[c][  1])*dR +(pix[1][-w1]-pix[c][-w1])*dU + (pix[1][ w1]-pix[c][ w1])*dD ) / (dL+dR+dU+dD));
					pix[c][0] = CLIP(v0);
				}
			}
		}
 /* Reinforce integrated red/blue pixels on BLUE/RED pixel locations */
#ifdef _OPENMP
#pragma omp for
#endif
		for (int row=2; row < height-2; row++) {
			int col = 2+(FC(row,2) & 1);
			int c = 2-FC(row,col);
#ifdef __SSE2__
			__m128 dLv, dRv, dUv, dDv, v0v;
			__m128 onev = _mm_set1_ps(1.f);
			__m128 zd5v = _mm_set1_ps(0.5f);
			__m128 c65535v = _mm_set1_ps(65535.f);
			for (; col < width-8; col+=8) {
				int indx = row*width+col;
				pix[0] = (float*)(*rgb[0]) + indx;
				pix[1] = (float*)(*rgb[1]) + indx;
				pix[2] = (float*)(*rgb[2]) + indx;
				int d = 2 - c;
				dLv = onev/(onev+vabsf(LC2VFU(pix[d][ -2])-LC2VFU(pix[d][0]))+vabsf(LC2VFU(pix[1][ 1])-LC2VFU(pix[1][ -1])));
				dRv = onev/(onev+vabsf(LC2VFU(pix[d][  2])-LC2VFU(pix[d][0]))+vabsf(LC2VFU(pix[1][ 1])-LC2VFU(pix[1][ -1])));
				dUv = onev/(onev+vabsf(LC2VFU(pix[d][-w2])-LC2VFU(pix[d][0]))+vabsf(LC2VFU(pix[1][w1])-LC2VFU(pix[1][-w1])));
				dDv = onev/(onev+vabsf(LC2VFU(pix[d][ w2])-LC2VFU(pix[d][0]))+vabsf(LC2VFU(pix[1][w1])-LC2VFU(pix[1][-w1])));
				v0v = CLIPV(LC2VFU(pix[1][0]) + zd5v -((LC2VFU(pix[1][-1])-LC2VFU(pix[c][-1]))*dLv +(LC2VFU(pix[1][1])-LC2VFU(pix[c][1]))*dRv +(LC2VFU(pix[1][-w1])-LC2VFU(pix[c][-w1]))*dUv + (LC2VFU(pix[1][w1])-LC2VFU(pix[c][w1]))*dDv ) / (dLv+dRv+dUv+dDv));
				STC2VFU(pix[c][0],v0v);
			}
#endif
			for (; col < width-2; col+=2) {
				int indx = row*width+col;
				pix[0] = (float*)(*rgb[0]) + indx;
				pix[1] = (float*)(*rgb[1]) + indx;
				pix[2] = (float*)(*rgb[2]) + indx;
				int d = 2 - c;
				float dL = 1.f/(1.f+fabsf(pix[d][ -2]-pix[d][0])+fabsf(pix[1][ 1]-pix[1][ -1]));
				float dR = 1.f/(1.f+fabsf(pix[d][  2]-pix[d][0])+fabsf(pix[1][ 1]-pix[1][ -1]));
				float dU = 1.f/(1.f+fabsf(pix[d][-w2]-pix[d][0])+fabsf(pix[1][w1]-pix[1][-w1]));
				float dD = 1.f/(1.f+fabsf(pix[d][ w2]-pix[d][0])+fabsf(pix[1][w1]-pix[1][-w1]));
				float v0 = (pix[1][0] + 0.5f -((pix[1][ -1]-pix[c][ -1])*dL +(pix[1][  1]-pix[c][  1])*dR +(pix[1][-w1]-pix[c][-w1])*dU + (pix[1][ w1]-pix[c][ w1])*dD ) / (dL+dR+dU+dD));
				pix[c][0] = CLIP(v0);
			}
		}
} // end parallel
	}
    t2e.set();
    if (settings->verbose) printf("Refinement Lee %d usec\n", t2e.etime(t1e));
}
#ifdef __SSE2__
#undef CLIPV
#endif


// Refinement based on EECI demozaicing algorithm by L. Chang and Y.P. Tan
// from "Lassus" : Luis Sanz Rodriguez, adapted by Jacques Desmis - JDC - and Oliver Duis for RawTherapee
// increases the signal to noise ratio (PSNR) # +1 to +2 dB : tested with Dcraw : eg: Lighthouse + AMaZE : whitout refinement:39.96dB, with refinement:41.86 dB
// reduce color artifacts, improves the interpolation
// but it's relatively slow
//
// Should be DISABLED if it decreases image quality by increases some image noise and generates blocky edges
void RawImageSource::refinement_lassus(int PassCount)
{
   // const int PassCount=1;

   // if (settings->verbose) printf("Refinement\n");

    MyTime t1e,t2e;
    t1e.set();
    int u=W, v=2*u, w=3*u, x=4*u, y=5*u;
    float (*image)[3];
    image = (float(*)[3]) calloc(W*H, sizeof *image);
#ifdef _OPENMP
#pragma omp parallel shared(image)
#endif
    {
        // convert red, blue, green to image
#ifdef _OPENMP
#pragma omp for
#endif
        for (int i=0;i<H;i++) {
            for (int j=0;j<W;j++) {
                image[i*W+j][0] = red  [i][j];
                image[i*W+j][1] = green[i][j];
                image[i*W+j][2] = blue [i][j];
            }
        }

        for (int b=0; b<PassCount; b++) {
            if (plistener) {
                plistener->setProgressStr (M("TP_RAW_DMETHOD_PROGRESSBAR_REFINE"));
                plistener->setProgress ((float)b/PassCount);
            }

            // Reinforce interpolated green pixels on RED/BLUE pixel locations
#ifdef _OPENMP
#pragma omp for
#endif
            for (int row=6; row<H-6; row++) {
                for (int col=6+(FC(row,2)&1),c=FC(row,col); col<W-6; col+=2) {
                    float (*pix)[3]=image+row*W+col;

                    // Cubic Spline Interpolation by Li and Randhawa, modified by Luis Sanz Rodriguez

                    float f[4];
                    f[0]=1.0f/(1.0f+xmul2f(fabs(x1125(pix[-v][c])-x0875(pix[0][c])-x0250(pix[-x][c])))+fabs(x0875(pix[u][1])-x1125(pix[-u][1])+x0250(pix[-w][1]))+fabs(x0875(pix[-w][1])-x1125(pix[-u][1])+x0250(pix[-y][1])));
                    f[1]=1.0f/(1.0f+xmul2f(fabs(x1125(pix[+2][c])-x0875(pix[0][c])-x0250(pix[+4][c])))+fabs(x0875(pix[1][1])-x1125(pix[-1][1])+x0250(pix[+3][1]))+fabs(x0875(pix[+3][1])-x1125(pix[+1][1])+x0250(pix[+5][1])));
                    f[2]=1.0f/(1.0f+xmul2f(fabs(x1125(pix[-2][c])-x0875(pix[0][c])-x0250(pix[-4][c])))+fabs(x0875(pix[1][1])-x1125(pix[-1][1])+x0250(pix[-3][1]))+fabs(x0875(pix[-3][1])-x1125(pix[-1][1])+x0250(pix[-5][1])));
                    f[3]=1.0f/(1.0f+xmul2f(fabs(x1125(pix[+v][c])-x0875(pix[0][c])-x0250(pix[+x][c])))+fabs(x0875(pix[u][1])-x1125(pix[-u][1])+x0250(pix[+w][1]))+fabs(x0875(pix[+w][1])-x1125(pix[+u][1])+x0250(pix[+y][1])));

                    float g[4];//CLIREF avoid overflow
                    g[0]=pix[0][c]+(x0875(CLIREF(pix[-u][1]-pix[-u][c]))+x0125(CLIREF(pix[+u][1]-pix[+u][c])));
                    g[1]=pix[0][c]+(x0875(CLIREF(pix[+1][1]-pix[+1][c]))+x0125(CLIREF(pix[-1][1]-pix[-1][c])));
                    g[2]=pix[0][c]+(x0875(CLIREF(pix[-1][1]-pix[-1][c]))+x0125(CLIREF(pix[+1][1]-pix[+1][c])));
                    g[3]=pix[0][c]+(x0875(CLIREF(pix[+u][1]-pix[+u][c]))+x0125(CLIREF(pix[-u][1]-pix[-u][c])));

                    pix[0][1]=(f[0]*g[0]+f[1]*g[1]+f[2]*g[2]+f[3]*g[3]) / (f[0]+f[1]+f[2]+f[3]);

                }
            }
            // Reinforce interpolated red/blue pixels on GREEN pixel locations
#ifdef _OPENMP
#pragma omp for
#endif
            for (int row=6; row<H-6; row++) {
                for (int col=6+(FC(row,3)&1),c=FC(row,col+1); col<W-6; col+=2) {
                    float (*pix)[3]=image+row*W+col;
                    for (int i=0; i<2; c=2-c,i++) {
                        float f[4];
                        f[0]=1.0f/(1.0f+xmul2f(fabs(x0875(pix[-v][1])-x1125(pix[0][1])+x0250(pix[-x][1])))+fabs(pix[u] [c]-pix[-u][c])+fabs(pix[-w][c]-pix[-u][c]));
                        f[1]=1.0f/(1.0f+xmul2f(fabs(x0875(pix[+2][1])-x1125(pix[0][1])+x0250(pix[+4][1])))+fabs(pix[+1][c]-pix[-1][c])+fabs(pix[+3][c]-pix[+1][c]));
                        f[2]=1.0f/(1.0f+xmul2f(fabs(x0875(pix[-2][1])-x1125(pix[0][1])+x0250(pix[-4][1])))+fabs(pix[+1][c]-pix[-1][c])+fabs(pix[-3][c]-pix[-1][c]));
                        f[3]=1.0f/(1.0f+xmul2f(fabs(x0875(pix[+v][1])-x1125(pix[0][1])+x0250(pix[+x][1])))+fabs(pix[u] [c]-pix[-u][c])+fabs(pix[+w][c]-pix[+u][c]));

                        float g[5];//CLIREF avoid overflow
                        g[0]=CLIREF(pix[-u][1]-pix[-u][c]);
                        g[1]=CLIREF(pix[+1][1]-pix[+1][c]);
                        g[2]=CLIREF(pix[-1][1]-pix[-1][c]);
                        g[3]=CLIREF(pix[+u][1]-pix[+u][c]);
                        g[4]=((f[0]*g[0]+f[1]*g[1]+f[2]*g[2]+f[3]*g[3]) / (f[0]+f[1]+f[2]+f[3]));
                        pix[0][c]= pix[0][1]-(0.65f*g[4]+0.35f*CLIREF(pix[0][1]-pix[0][c]));
                    }
                }
            }
            // Reinforce integrated red/blue pixels on BLUE/RED pixel locations
#ifdef _OPENMP
#pragma omp for
#endif
            for (int row=6; row<H-6; row++) {
                for (int col=6+(FC(row,2)&1),c=2-FC(row,col),d=2-c; col<W-6; col+=2) {
                    float (*pix)[3]=image+row*W+col;

                    float f[4];
                    f[0]=1.0f/(1.0f+xmul2f(fabs(x1125(pix[-v][d])-x0875(pix[0][d])-x0250(pix[-x][d])))+fabs(x0875(pix[u][1])-x1125(pix[-u][1])+x0250(pix[-w][1]))+fabs(x0875(pix[-w][1])-x1125(pix[-u][1])+x0250(pix[-y][1])));
                    f[1]=1.0f/(1.0f+xmul2f(fabs(x1125(pix[+2][d])-x0875(pix[0][d])-x0250(pix[+4][d])))+fabs(x0875(pix[1][1])-x1125(pix[-1][1])+x0250(pix[+3][1]))+fabs(x0875(pix[+3][1])-x1125(pix[+1][1])+x0250(pix[+5][1])));
                    f[2]=1.0f/(1.0f+xmul2f(fabs(x1125(pix[-2][d])-x0875(pix[0][d])-x0250(pix[-4][d])))+fabs(x0875(pix[1][1])-x1125(pix[-1][1])+x0250(pix[-3][1]))+fabs(x0875(pix[-3][1])-x1125(pix[-1][1])+x0250(pix[-5][1])));
                    f[3]=1.0f/(1.0f+xmul2f(fabs(x1125(pix[+v][d])-x0875(pix[0][d])-x0250(pix[+x][d])))+fabs(x0875(pix[u][1])-x1125(pix[-u][1])+x0250(pix[+w][1]))+fabs(x0875(pix[+w][1])-x1125(pix[+u][1])+x0250(pix[+y][1])));

                    float g[5];
                    g[0]=(x0875((pix[-u][1]-pix[-u][c]))+x0125((pix[-v][1]-pix[-v][c])));
                    g[1]=(x0875((pix[+1][1]-pix[+1][c]))+x0125((pix[+2][1]-pix[+2][c])));
                    g[2]=(x0875((pix[-1][1]-pix[-1][c]))+x0125((pix[-2][1]-pix[-2][c])));
                    g[3]=(x0875((pix[+u][1]-pix[+u][c]))+x0125((pix[+v][1]-pix[+v][c])));

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
                    pix[0][c]=LIM(pix[0][1]-(1.30f*g[4]-0.30f*(pix[0][1]-pix[0][c])), 0.99f*(pix[0][1]-p[4]), 1.01f*(pix[0][1]-p[4]));

                }
            }

        }

        // put modified values to red, green, blue
#ifdef _OPENMP
#pragma omp for
#endif
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
    if (settings->verbose) printf("Refinement Lassus %d usec\n", t2e.etime(t1e));
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
	const unsigned int colors = 3;  // used in FORCC

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
		for (int col = colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin)&1),indx=row*CACHESIZE+col; col < colMax; col+=2, indx+=2) {
            assert(indx-u>=0 && indx+u<u*u);
			bufferH[indx][1] = (image[indx-1][1] + image[indx+1][1]) * 0.5f;
			bufferV[indx][1] = (image[indx+u][1] + image[indx-u][1]) * 0.5f;
		}
	}
	// red in blue pixel, blue in red pixel
	for (int row=rowMin; row < rowMax; row++)
		for (int col=colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin) & 1), indx=row*CACHESIZE+col, c=2-FC(y0-TILEBORDER+row,x0-TILEBORDER+col); col < colMax; col+=2, indx+=2) {
            assert(indx-u-1>=0 && indx+u+1<u*u && c>=0 && c<3);

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
            assert(indx-u>=0 && indx+u<u*u && c>=0 && c<3 && d>=0 && d<3);
			bufferH[indx][c] = (image[indx+1][c] + image[indx-1][c]) * 0.5f;
			bufferH[indx][d] = (2.f * bufferH[indx][1] - bufferH[indx+u][1] - bufferH[indx-u][1] + image[indx+u][d] + image[indx-u][d]) * 0.5f;
			bufferV[indx][c] = (2.f * bufferV[indx][1] - bufferV[indx+1][1] - bufferV[indx-1][1] + image[indx+1][c] + image[indx-1][c]) * 0.5f;
			bufferV[indx][d] = (image[indx+u][d] + image[indx-u][d]) * 0.5f;
		}

    // Decide green pixels
    for (int row = rowMin; row < rowMax; row++)
        for (int col = colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin)&1),indx=row*CACHESIZE+col,c=FC(y0-TILEBORDER+row,x0-TILEBORDER+col),d=2-c; col < colMax; col+=2, indx+=2) {
            float current =  max(image[indx+v][c], image[indx-v][c], image[indx-2][c], image[indx+2][c]) -
                             min(image[indx+v][c], image[indx-v][c], image[indx-2][c], image[indx+2][c]) +
                             max(image[indx+1+u][d], image[indx+1-u][d], image[indx-1+u][d], image[indx-1-u][d]) -
                             min(image[indx+1+u][d], image[indx+1-u][d], image[indx-1+u][d], image[indx-1-u][d]);

            float currentH = max(bufferH[indx+v][d], bufferH[indx-v][d], bufferH[indx-2][d], bufferH[indx+2][d]) -
                             min(bufferH[indx+v][d], bufferH[indx-v][d], bufferH[indx-2][d], bufferH[indx+2][d]) +
                             max(bufferH[indx+1+u][c], bufferH[indx+1-u][c], bufferH[indx-1+u][c], bufferH[indx-1-u][c]) -
                             min(bufferH[indx+1+u][c], bufferH[indx+1-u][c], bufferH[indx-1+u][c], bufferH[indx-1-u][c]);

            float currentV = max(bufferV[indx+v][d], bufferV[indx-v][d], bufferV[indx-2][d], bufferV[indx+2][d]) -
                             min(bufferV[indx+v][d], bufferV[indx-v][d], bufferV[indx-2][d], bufferV[indx+2][d]) +
                             max(bufferV[indx+1+u][c], bufferV[indx+1-u][c], bufferV[indx-1+u][c], bufferV[indx-1-u][c]) -
                             min(bufferV[indx+1+u][c], bufferV[indx+1-u][c], bufferV[indx-1+u][c], bufferV[indx-1-u][c]);

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
            assert(indx-v>=0 && indx+v<u*u);
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
				image[indx][3] = ((min(pix[-4], pix[+4]) + pix[-4] + pix[+4] ) < (min(pix[-u], pix[+u]) + pix[-u] + pix[+u]));
			else
				image[indx][3] = ((max(pix[-4], pix[+4]) + pix[-4] + pix[+4] ) > (max(pix[-u], pix[+u]) + pix[-u] + pix[+u]));
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
		for (int col = colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin)&1),indx=row*CACHESIZE+col; col < colMax; col+=2, indx+=2) {
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

			g1 = (f[0] + f[1] + f[2] + f[3] + f[4] - max(f[1], f[2], f[3], f[4]) - min(f[1], f[2], f[3], f[4])) / 3.f;

			f[0] = (float)(image[indx-1][1] + image[indx+1][1])/(2.f + 2.f * image[indx][c]);
			f[1] = 2.f * image[indx-1][1]/(2 + image[indx-2][c] + image[indx][c]);
			f[2] = (float)(image[indx-1][1] + image[indx-3][1])/(2.f + 2.f * image[indx-2][c]);
			f[3] = 2.f * image[indx+1][1]/(2 + image[indx+2][c] + image[indx][c]);
			f[4] = (float)(image[indx+1][1] + image[indx+3][1])/(2.f + 2.f * image[indx+2][c]);

			g2 = (f[0] + f[1] + f[2] + f[3] + f[4] - max(f[1], f[2], f[3], f[4]) - min(f[1], f[2], f[3], f[4])) / 3.f;

            assert(indx>=0 && indx<u*u);
			image[indx][1] = (2.f+image[indx][c]) * (current*g1 + (16.f-current)*g2) * 0.0625f;

            // get rid of the overshooted pixels
		    float min_f = min(image[indx+1+u][1], min(image[indx+1-u][1], min(image[indx-1+u][1], min(image[indx-1-u][1], min(image[indx-1][1], min(image[indx+1][1], min(image[indx-u][1], image[indx+u][1])))))));
		    float max_f = max(image[indx+1+u][1], max(image[indx+1-u][1], max(image[indx-1+u][1], max(image[indx-1-u][1], max(image[indx-1][1], max(image[indx+1][1], max(image[indx-u][1], image[indx+u][1])))))));

		    image[indx][1] =  LIM(image[indx][1], min_f, max_f);
		}
}

// missing colors are interpolated using high quality algorithm by Luis Sanz Rodriguez
void RawImageSource::dcb_color_full(float (*image)[4], int x0, int y0, float (*chroma)[2])
{
	const int u=CACHESIZE, w=3*CACHESIZE;
	int rowMin,colMin,rowMax,colMax;
	dcb_initTileLimits(colMin,rowMin,colMax,rowMax,x0,y0,3);

	float f[4],g[4];

	for (int row=1; row < CACHESIZE-1; row++)
		for (int col=1+(FC(y0-TILEBORDER+row,x0-TILEBORDER+1)&1),indx=row*CACHESIZE+col,c=FC(y0-TILEBORDER+row,x0-TILEBORDER+col),d=c/2; col < CACHESIZE-1; col+=2,indx+=2) {
            assert(indx>=0 && indx<u*u && c>=0 && c<4);
			chroma[indx][d]=image[indx][c]-image[indx][1];
        }

	for (int row=rowMin; row<rowMax; row++)
		for (int col=colMin+(FC(y0-TILEBORDER+row,x0-TILEBORDER+colMin)&1),indx=row*CACHESIZE+col,c=1-FC(y0-TILEBORDER+row,x0-TILEBORDER+col)/2; col<colMax; col+=2,indx+=2) {
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
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::dcb]));
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

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for( int iTile=0; iTile < numTiles; iTile++){
    	int xTile = iTile % wTiles;
    	int yTile = iTile / wTiles;
    	int x0 = xTile*TILESIZE;
    	int y0 = yTile*TILESIZE;

#ifdef _OPENMP
    	int tid = omp_get_thread_num();
        assert(tid<nthreads);
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
#ifdef _OPENMP
#pragma omp atomic
#endif
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

const double xyz_rgb[3][3] = {			// XYZ from RGB
  { 0.412453, 0.357580, 0.180423 },
  { 0.212671, 0.715160, 0.072169 },
  { 0.019334, 0.119193, 0.950227 } };
const float d65_white[3] = { 0.950456, 1, 1.088754 };

void RawImageSource::cielab (const float (*rgb)[3], float* l, float* a, float *b, const int width, const int height, const int labWidth, const float xyz_cam[3][3])
{
	static float cbrt[0x10000];
	static bool cbrtinit = false;
	if (!rgb) { 
		int i, j, k;
		float r;
		if(!cbrtinit) {
			for (i=0; i < 0x10000; i++) {
				r = i / 65535.f;
				cbrt[i] = r > 0.008856f ? xcbrtf(r) : 7.787f*r + 16.f/116.f;
			}
			cbrtinit = true;
		}
		return;
	}

	int rgbOffset = (width - labWidth);
	for(int i=0;i<height;i++) {
		for(int j=0;j<labWidth;j++) {
			float xyz[3] = {0.5f};
			int c;
			FORC3 {
				xyz[0] += xyz_cam[0][c] * rgb[i*width+j][c];
				xyz[1] += xyz_cam[1][c] * rgb[i*width+j][c];
				xyz[2] += xyz_cam[2][c] * rgb[i*width+j][c];
			}
			xyz[0] = cbrt[CLIP((int) xyz[0])];
			xyz[1] = cbrt[CLIP((int) xyz[1])];
			xyz[2] = cbrt[CLIP((int) xyz[2])];

			l[i*labWidth+j] = 116 * xyz[1] - 16;
			a[i*labWidth+j] = 500 * (xyz[0] - xyz[1]);
			b[i*labWidth+j] = 200 * (xyz[1] - xyz[2]);
		}
	}
}

#define fcol(row,col) xtrans[(row)%6][(col)%6]

void RawImageSource::xtransborder_interpolate (int border)
{
	const int height = H, width = W;
	
	char xtrans[6][6];
	ri->getXtransMatrix(xtrans);

	for (int row=0; row < height; row++)
		for (int col=0; col < width; col++) {
			if (col==border && row >= border && row < height-border)
				col = width-border;
			float sum[6] = {0.f};
			for (int y=MAX(0,row-1); y <= MIN(row+1,height-1); y++)
				for (int x=MAX(0,col-1); x <= MIN(col+1,width-1); x++) {
					int f = fcol(y,x);
					sum[f] += rawData[y][x];
					sum[f+3]++;
				}

			switch(fcol(row,col)) {
				case 0:	red[row][col] = rawData[row][col];
						green[row][col] = (sum[1]/sum[4]);
						blue[row][col] = (sum[2]/sum[5]);
						break;
				case 1:	if(sum[3]==0.f) { // at the 4 corner pixels it can happen, that we have only green pixels in 2x2 area
							red[row][col] = green[row][col] = blue[row][col] = rawData[row][col];
						} else {
							red[row][col] = (sum[0]/sum[3]);
							green[row][col] = rawData[row][col];
							blue[row][col] = (sum[2]/sum[5]);
						}
						break;
				case 2:	red[row][col] = (sum[0]/sum[3]);
						green[row][col] = (sum[1]/sum[4]);
						blue[row][col] = rawData[row][col];
			}
		}
}

/*
   Frank Markesteijn's algorithm for Fuji X-Trans sensors
   adapted to RT by Ingo Weyrich 2014
*/

#define TS 122		/* Tile Size */

void RawImageSource::xtrans_interpolate (int passes, bool useCieLab)
{
	double progress = 0.0;
	const bool plistenerActive = plistener;
  
 	if (plistenerActive) {
		plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), "Xtrans"));
		plistener->setProgress (progress);
	}

	char xtrans[6][6];
	ri->getXtransMatrix(xtrans);

	static const short  orth[12] = { 1,0,0,1,-1,0,0,-1,1,0,0,1 },
						patt[2][16] = { { 0,1,0,-1,2,0,-1,0,1,1,1,-1,0,0,0,0 },
									 { 0,1,0,-2,1,0,-2,0,1,1,-2,-2,1,-1,-1,1 } },
						dir[4] = { 1,TS,TS+1,TS-1 };

	short allhex[2][3][3][8];

	// sgrow/sgcol is the offset in the sensor matrix of the solitary
	// green pixels 
	ushort sgrow, sgcol;
	
	const int height = H, width = W;

	if (settings->verbose)
		printf("%d-pass X-Trans interpolation using %s conversion...\n", passes, useCieLab ? "lab" : "yuv");

	xtransborder_interpolate(6);

	float xyz_cam[3][3];
	{
	float rgb_cam[3][4];
	ri->getRgbCam(rgb_cam);
	int k;
	for (int i=0; i < 3; i++)
		for (int j=0; j < 3; j++)
			for (xyz_cam[i][j] = k=0; k < 3; k++)
				xyz_cam[i][j] += xyz_rgb[i][k] * rgb_cam[k][j] / d65_white[i];
	}

/* Map a green hexagon around each non-green pixel and vice versa:	*/
{
	int gint, d, h, v, ng, row, col, c;

	for (row=0; row < 3; row++)
		for (col=0; col < 3; col++) {
			gint = fcol(row,col) == 1;
			for (ng=d=0; d < 10; d+=2) {
				if (fcol(row+orth[d]+6,col+orth[d+2]+6) == 1)
					ng=0;
				else
					ng++;
				if (ng == 4) {
					// if there are four non-green pixels adjacent in cardinal
					// directions, this is the solitary green pixel
					sgrow = row;
					sgcol = col;
				}
				if (ng == gint+1)
					FORC(8) {
						v = orth[d]*patt[gint][c*2] + orth[d+1]*patt[gint][c*2+1];
						h = orth[d+2]*patt[gint][c*2] + orth[d+3]*patt[gint][c*2+1];
						allhex[0][row][col][c^(gint*2 & d)] = h + v*width;
						allhex[1][row][col][c^(gint*2 & d)] = h + v*TS;
					}
			}
		}

}
	if(plistenerActive) {
		progress += 0.05;
		plistener->setProgress(progress);
	} 


	double progressInc = 36.0*(1.0-progress)/((H*W)/((TS-16)*(TS-16)));
	const int ndir = 4 << (passes > 1);
	cielab (0,0,0,0,0,0,0,0);
	struct s_minmaxgreen {
		float min;
		float max;
	};

	int	RightShift[6];
	for(int row=0;row<6;row++) {
		// count number of green pixels in three cols
		int greencount = 0;
		for(int col=0;col<3;col++)
			greencount += (fcol(row,col) == 1);
		RightShift[row] = (greencount == 2);
	}

	
#pragma omp parallel
{
	int progressCounter = 0;
	short *hex;
    int c, d, f, h, i, v, mrow, mcol;
    int pass;
	float color[3][8], g, val;
	float (*rgb)[TS][TS][3], (*rix)[3];
	float (*lab)[TS-8][TS-8];
	float (*drv)[TS-10][TS-10], diff[6], tr;
	s_minmaxgreen  (*greenminmaxtile)[TS];
	uint8_t (*homo)[TS][TS];
	uint8_t (*homosum)[TS][TS];
	float *buffer;
	buffer = (float *) malloc ((TS*TS*(ndir*3+11)+128)*sizeof(float));
	rgb  = (float(*)[TS][TS][3]) buffer;
	lab  = (float (*)    [TS-8][TS-8])(buffer + TS*TS*(ndir*3));
	drv  = (float (*)[TS-10][TS-10])   (buffer + TS*TS*(ndir*3+3));
	homo = (uint8_t  (*)[TS][TS])   (lab); // we can reuse the lab-buffer because they are not used together
	greenminmaxtile = (s_minmaxgreen(*)[TS]) (lab); // we can reuse the lab-buffer because they are not used together
	homosum = (uint8_t (*)[TS][TS]) (drv); // we can reuse the drv-buffer because they are not used together

#pragma omp for collapse(2)	schedule(dynamic) nowait
	for (int top=3; top < height-19; top += TS-16)
		for (int left=3; left < width-19; left += TS-16) {
			int mrow = MIN (top+TS, height-3);
			int mcol = MIN (left+TS, width-3);
			memset(rgb,0,TS*TS*3*sizeof(float));
			for (int row=top; row < mrow; row++)
				for (int col=left; col < mcol; col++) {
					rgb[0][row-top][col-left][fcol(row,col)] = rawData[row][col];
				}
			FORC3 memcpy (rgb[c+1], rgb[0], sizeof *rgb);

			/* Set green1 and green3 to the minimum and maximum allowed values:	*/
			for (int row=top; row < mrow; row++) {
				float minval=FLT_MAX;
				float maxval=0.f;
				int shiftindex = RightShift[(row)%6];
				for (int col=left; col < mcol; col++) {
					if (fcol(row,col) == 1) {
						minval=FLT_MAX;
						maxval=0.f;
						continue;
					}
					float *pix = &rawData[row][col];
					hex = allhex[0][row % 3][col % 3];
					if (maxval==0.f)
						FORC(6) {
							val = pix[hex[c]];
							if (minval > val)
								minval = val;
							if (maxval < val)
								maxval = val;
						}
					greenminmaxtile[row-top][(col-left)>>shiftindex].min = minval;
					greenminmaxtile[row-top][(col-left)>>shiftindex].max = maxval;
					switch ((row-sgrow) % 3) {
						case 1: if (row < mrow-1) {
									row++;
									shiftindex = RightShift[(row)%6];
									col--;
								}
								break;
						case 2: minval=FLT_MAX;
								maxval=0.f;
								if ((col+=2) < mcol-1 && row > top+1) {
									row--;
									shiftindex = RightShift[(row)%6];
								}
					}
				}
			}

			/* Interpolate green horizontally, vertically, and along both diagonals: */
			for (int row=top; row < mrow; row++) {
				// find first non-green pixel
				int leftstart = left;
				for(;leftstart<mcol;leftstart++)
					if(fcol(row,leftstart)!=1)
						break;
				const int shiftindex = RightShift[(row)%6];
				const int coloffset = (shiftindex == 1 ? 3:1);
				for (int col=leftstart; col < mcol; col+=coloffset) {
					if (fcol(row,col) == 1)
						continue;
					float *pix = &rawData[row][col];
					hex = allhex[0][row % 3][col % 3];
					color[1][0] = 0.6796875f * (pix[hex[1]] + pix[hex[0]]) -
						 0.1796875f * (pix[2*hex[1]] + pix[2*hex[0]]);
					color[1][1] = 0.87109375f *  pix[hex[3]] + pix[hex[2]] * 0.12890625f +
						 0.359375f * (pix[0] - pix[-hex[2]]);
					FORC(2)
						color[1][2+c] = 0.640625f * pix[hex[4+c]] + 0.359375f * pix[-2*hex[4+c]] + 0.12890625f *
										(2.f*pix[0] - pix[3*hex[4+c]] - pix[-3*hex[4+c]]);
					FORC(4)
						rgb[c^!((row-sgrow) % 3)][row-top][col-left][1] = LIM(color[1][c],greenminmaxtile[row-top][(col-left)>>shiftindex].min,greenminmaxtile[row-top][(col-left)>>shiftindex].max);
				}
			}
			for (pass=0; pass < passes; pass++) {
				if (pass == 1)
					memcpy (rgb+=4, buffer, 4*sizeof *rgb);

/* Recalculate green from interpolated values of closer pixels:	*/
				if (pass) {
					for (int row=top+2; row < mrow-2; row++) {
						int leftstart = left+2;
						for(;leftstart<mcol-2;leftstart++)
							if(fcol(row,leftstart)!=1)
								break;
						const int shiftindex = RightShift[(row)%6];
						const int coloffset = (shiftindex == 1 ? 3:1);
						for (int col=leftstart; col < mcol-2; col+=coloffset) {
							if ((f = fcol(row,col)) == 1)
								continue;
							hex = allhex[1][row % 3][col % 3];
							for (d=3; d < 6; d++) {
								rix = &rgb[(d-2)^!((row-sgrow) % 3)][row-top][col-left];
								val = rix[-2*hex[d]][1] + 2*(rix[hex[d]][1] - rix[hex[d]][f])
									- rix[-2*hex[d]][f] + 3*rix[0][f];
								rix[0][1] = LIM((float)(val*.33333333f),greenminmaxtile[row-top][(col-left)>>shiftindex].min,greenminmaxtile[row-top][(col-left)>>shiftindex].max);
							}
						}
					}
				}

/* Interpolate red and blue values for solitary green pixels:	*/
				for (int row=(top-sgrow+4)/3*3+sgrow; row < mrow-2; row+=3)
					for (int col=(left-sgcol+4)/3*3+sgcol; col < mcol-2; col+=3) {
						rix = &rgb[0][row-top][col-left];
						h = fcol(row,col+1);
						memset (diff, 0, sizeof diff);
						for (i=1, d=0; d < 6; d++, i^=TS^1, h^=2) {
							for (c=0; c < 2; c++, h^=2) {
								g = rix[0][1] + rix[0][1] - rix[i<<c][1] - rix[-i<<c][1];
								color[h][d] = g + rix[i<<c][h] + rix[-i<<c][h];
								if (d > 1)
									diff[d] += SQR (rix[i<<c][1] - rix[-i<<c][1]
											- rix[i<<c][h] + rix[-i<<c][h]) + SQR(g);
							}
							if (d > 2 && (d & 1))	 // 3, 5
								if (diff[d-1] < diff[d])
									FORC(2)
										color[c*2][d] = color[c*2][d-1];
							if ((d & 1) || d < 2) { // d: 0, 1, 3, 5
								FORC(2)
									rix[0][c*2] = CLIP(0.5f*color[c*2][d]);
								rix += TS*TS;
							}
						}
					}

/* Interpolate red for blue pixels and vice versa:		*/
				for (int row=top+3; row < mrow-3; row++) {
					int leftstart = left+3;
					for(;leftstart<mcol-1;leftstart++)
						if(fcol(row,leftstart)!=1)
							break;
					const int coloffset = (RightShift[(row)%6] == 1 ? 3:1);
					c = (row-sgrow) % 3 ? TS:1;
					h = 3 * (c ^ TS ^ 1);
					for (int col=leftstart; col < mcol-3; col+=coloffset) {
						if ((f = 2-fcol(row,col)) == 1)
							continue;
						rix = &rgb[0][row-top][col-left];
						for (d=0; d < 4; d++, rix += TS*TS) {
							i = d > 1 || ((d ^ c) & 1) ||
								((fabsf(rix[0][1]-rix[c][1])+fabsf(rix[0][1]-rix[-c][1])) <	2.f*(fabsf(rix[0][1]-rix[h][1])+fabsf(rix[0][1]-rix[-h][1]))) ? c:h;

							rix[0][f] = CLIP(0.5f*(rix[i][f] + rix[-i][f] +
										rix[0][1] + rix[0][1] - rix[i][1] - rix[-i][1]));
						}
					}
				}

/* Fill in red and blue for 2x2 blocks of green:		*/
				for (int row=top+2; row < mrow-2; row++)
					if ((row-sgrow) % 3) {
						for (int col=left+2; col < mcol-2; col++)
							if ((col-sgcol) % 3) {
								rix = &rgb[0][row-top][col-left];
								hex = allhex[1][row%3][col % 3];
								for (d=0; d < ndir; d+=2, rix += TS*TS)
									if (hex[d] + hex[d+1]) {
										g = 3*rix[0][1] - 2*rix[hex[d]][1] - rix[hex[d+1]][1];
										for (c=0; c < 4; c+=2)
											rix[0][c] =	CLIP((g + 2*rix[hex[d]][c] + rix[hex[d+1]][c])*0.33333333f);
									} else {
										g = 2*rix[0][1] - rix[hex[d]][1] - rix[hex[d+1]][1];
										for (c=0; c < 4; c+=2)
											rix[0][c] =	CLIP((g + rix[hex[d]][c] + rix[hex[d+1]][c])*0.5f);
									}
							}
					}
			}
			
// end of multipass part
			rgb = (float(*)[TS][TS][3]) buffer;
			mrow -= top;
			mcol -= left;

			if(useCieLab) {
				/* Convert to CIELab and differentiate in all directions:	*/
				// Original dcraw algorithm uses CIELab as perceptual space
				// (presumably coming from original AHD) and converts taking
				// camera matrix into account.  We use this in RT.
				for (d=0; d < ndir; d++) {
					float *l = &lab[0][0][0];
					float *a = &lab[1][0][0];
					float *b = &lab[2][0][0];
					cielab(&rgb[d][4][4],l,a,b,TS,mrow-8,TS-8,xyz_cam);
					int f = dir[d & 3];
					f = f == 1 ? 1 : f-8;
					for (int row=5; row < mrow-5; row++)
						for (int col=5; col < mcol-5; col++) {
							float *l = &lab[0][row-4][col-4];
							float *a = &lab[1][row-4][col-4];
							float *b = &lab[2][row-4][col-4];
							
							g = 2*l[0] - l[f] - l[-f];
							drv[d][row-5][col-5] =	SQR(g)
												+ SQR((2*a[0] - a[f] - a[-f] + g*2.1551724f))
												+ SQR((2*b[0] - b[f] - b[-f] - g*0.86206896f));
						}

				}
			} else {
				// Now use YPbPr which requires much
				// less code and is nearly indistinguishable. It assumes the
				// camera RGB is roughly linear.
				// 
				for (d=0; d < ndir; d++) {
					float (*yuv)[TS-8][TS-8] = lab; // we use the lab buffer, which has the same dimensions
					for (int row=4; row < mrow-4; row++)
						for (int col=4; col < mcol-4; col++) {
							// use ITU-R BT.2020 YPbPr, which is great, but could use
							// a better/simpler choice? note that imageop.h provides
							// dt_iop_RGB_to_YCbCr which uses Rec. 601 conversion,
							// which appears less good with specular highlights
							float y = 0.2627f * rgb[d][row][col][0] + 0.6780f * rgb[d][row][col][1] + 0.0593f * rgb[d][row][col][2];
							yuv[0][row-4][col-4] = y;
							yuv[1][row-4][col-4] = (rgb[d][row][col][2]-y)*0.56433f;
							yuv[2][row-4][col-4] = (rgb[d][row][col][0]-y)*0.67815f;
						}
					int f=dir[d & 3];
					f = f == 1 ? 1 : f-8;
					for (int row=5; row < mrow-5; row++)
						for (int col=5; col < mcol-5; col++) {
							float *y = &yuv[0][row-4][col-4];
							float *u = &yuv[1][row-4][col-4];
							float *v = &yuv[2][row-4][col-4];
							drv[d][row-5][col-5] = SQR(2*y[0] - y[f] - y[-f])
												 + SQR(2*u[0] - u[f] - u[-f])
												 + SQR(2*v[0] - v[f] - v[-f]);
					  }
				}
			}

/* Build homogeneity maps from the derivatives:			*/
			memset(homo, 0, ndir*TS*TS*sizeof(uint8_t));
			for (int row=6; row < mrow-6; row++)
				for (int col=6; col < mcol-6; col++) {
					for (tr=FLT_MAX, d=0; d < ndir; d++)
						tr = (drv[d][row-5][col-5] < tr ? drv[d][row-5][col-5] : tr);
					tr *= 8;
					for (d=0; d < ndir; d++)
						for (v=-1; v <= 1; v++)
							for (h=-1; h <= 1; h++)
								homo[d][row][col] += (drv[d][row+v-5][col+h-5] <= tr ? 1:0) ;
				}

			if (height-top < TS+4)
				mrow = height-top+2;
			if (width-left < TS+4)
				mcol = width-left+2;


/* Build 5x5 sum of homogeneity maps */
			for(d=0;d<ndir;d++) {
				for (int row = MIN(top,8); row < mrow-8; row++) {
					int v5sum[5] = {0};
					const int startcol = MIN(left,8);
					for(v=-2;v<=2;v++)
						for(h=-2;h<=2;h++)
							v5sum[2+h] += homo[d][row+v][startcol+h];
					int blocksum = v5sum[0] + v5sum[1] + v5sum[2] + v5sum[3] + v5sum[4];
					homosum[d][row][startcol] = blocksum;
					int voffset = -1;
					// now we can subtract a column of five from blocksum and get new colsum of 5
					for (int col = startcol+1; col < mcol-8; col++) {
						int colsum = homo[d][row-2][col+2];
						for(v=-1;v<=2;v++)
							colsum += homo[d][row+v][col+2];
						voffset ++;
						voffset = voffset == 5 ? 0 : voffset;  // faster than voffset %= 5;
						blocksum -= v5sum[voffset];
						blocksum += colsum;
						v5sum[voffset] = colsum;
						homosum[d][row][col] = blocksum;
					}
				}
			}

/* Average the most homogenous pixels for the final result:	*/
			for (int row = MIN(top,8); row < mrow-8; row++)
				for (int col = MIN(left,8); col < mcol-8; col++) {
					uint8_t hm[8];
					uint8_t maxval = 0;
					for (d=0; d < 4; d++) {
						hm[d] = homosum[d][row][col];
						maxval = (maxval < hm[d] ? hm[d] : maxval);
					}
					for (; d < ndir; d++) {
						hm[d] = homosum[d][row][col];
						maxval = (maxval < hm[d] ? hm[d] : maxval);
						if (hm[d-4] < hm[d])
							hm[d-4] = 0;
						else if (hm[d-4] > hm[d])
							hm[d] = 0;
					}
					maxval -= maxval >> 3;
					float avg[4] = {0.f};
					for (d=0; d < ndir; d++)
						if (hm[d] >= maxval) {
							FORC3 avg[c] += rgb[d][row][col][c];
							avg[3]++;
						}

					red[row+top][col+left] = (avg[0]/avg[3]);
					green[row+top][col+left] = (avg[1]/avg[3]);
					blue[row+top][col+left] = (avg[2]/avg[3]);
			}
			
			if(plistenerActive && ((++progressCounter) % 32 == 0)) {
#ifdef _OPENMP
#pragma omp critical (xtransdemosaic)
#endif
{
				progress += progressInc;
				progress = min(1.0,progress);
				plistener->setProgress (progress);
}
			}


		}
	free(buffer);
}

}

#undef TS

void RawImageSource::fast_xtrans_interpolate ()
{
	if (settings->verbose)
		printf("fast X-Trans interpolation...\n");

	double progress = 0.0;
	const bool plistenerActive = plistener;
  
 	if (plistenerActive) {
		plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), "fast Xtrans"));
		plistener->setProgress (progress);
	}

	const int height = H, width = W;

	xtransborder_interpolate (1);
	char xtrans[6][6];
	ri->getXtransMatrix(xtrans);

#pragma omp parallel for
	for(int row=1;row<height-1;row++) {
		for(int col=1;col<width-1;col++) {
			float sum[3] = {0.f};
			for(int v=-1;v<=1;v++) {
				for(int h=-1;h<=1;h++) {
					sum[fcol(row+v,col+h)] += rawData[row+v][(col+h)];
				}
			}
			switch(fcol(row,col)) {
				case 0: red[row][col] = rawData[row][col];
						green[row][col] = sum[1] *0.2f;
						blue[row][col] = sum[2] *0.33333333f;
						break;
				case 1: red[row][col] = sum[0] * 0.5f;
						green[row][col] = rawData[row][col];
						blue[row][col] = sum[2] * 0.5f;
						break;
				case 2:	red[row][col] = sum[0] * 0.33333333f;
						green[row][col] = sum[1] *0.2f;
						blue[row][col] = rawData[row][col];
						break;
			}
		}
	}
	
 	if (plistenerActive) {
		plistener->setProgress (1.0);
 	}
}
#undef fcol



#undef TILEBORDER
#undef TILESIZE
#undef CACHESIZE
} /* namespace */
