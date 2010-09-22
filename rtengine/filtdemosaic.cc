/*
 * filtdemosaic.cc
 *
 *  Created on: Sep 17, 2010
 *      Author: gabor
 */

#include "filtdemosaic.h"
#include "rtengine.h"
#include "macros.h"
#include <string.h>
#include <omp.h>
#include "median.h"

namespace rtengine {

DemosaicFilterDescriptor demosaicFilterDescriptor;

DemosaicFilterDescriptor::DemosaicFilterDescriptor ()
	: FilterDescriptor ("Demosaicing", MultiImage::Raw, MultiImage::RGB) {

	addTriggerEvent (EvDMMethod);
	addTriggerEvent (EvDMColorCorrSteps);
}

void DemosaicFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new DemosaicFilter ());
}

DemosaicFilter::DemosaicFilter ()
	: Filter (&demosaicFilterDescriptor), border (4) {
}

ImageView DemosaicFilter::calculateTargetImageView (const ImageView& requestedImView) {

    ImageView result;
    result.skip = 1;
    result.x = 0;
    result.y = 0;
    getPreviousFilter()->getFullImageSize (result.w, result.h);
    result.w -= 2 * border;
    result.h -= 2 * border;

    return result;
}

ImageView DemosaicFilter::calculateSourceImageView (const ImageView& requestedImView) {

    return calculateTargetImageView (requestedImView);
}

void DemosaicFilter::getFullImageSize (int& w, int& h) {

    getPreviousFilter()->getFullImageSize (w, h);
    w -= 2 * border;
    h -= 2 * border;
}

void DemosaicFilter::getReqiredBufferSize (int& w, int& h) {

    getPreviousFilter()->getFullImageSize (w, h);
}

void DemosaicFilter::reverseTransPoint (int x, int y, int& xv, int& yv) {

    xv = x + border;
    yv = y + border;
}

void DemosaicFilter::hphd_vertical (MultiImage* si, Buffer<float>* hpmap, int col_from, int col_to) {

    int W = si->width;
    int H = si->height;

    float* temp = new float[std::max(W,H)];
    float* avg  = new float[std::max(W,H)];
    float* dev  = new float[std::max(W,H)];

    memset (temp, 0, std::max(W,H)*sizeof(float));
    memset (avg,  0, std::max(W,H)*sizeof(float));
    memset (dev,  0, std::max(W,H)*sizeof(float));

    for (int k=col_from; k<col_to; k++) {
        for (int i=5; i<H-5; i++) {
            temp[i] = (si->raw[i-5][k] - 8*si->raw[i-4][k] + 27*si->raw[i-3][k] - 48*si->raw[i-2][k] + 42*si->raw[i-1][k] -
                    (si->raw[i+5][k] - 8*si->raw[i+4][k] + 27*si->raw[i+3][k] - 48*si->raw[i+2][k] + 42*si->raw[i+1][k])) / 100.0;
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
            hpmap->rows[j][k] = avgL + (avgR - avgL) * devL / (devL + devR);
        }
    }

    delete [] temp;
    delete [] avg;
    delete [] dev;
}

void DemosaicFilter::hphd_horizontal (MultiImage* si, Buffer<float>* hpmap, int row_from, int row_to) {

    int W = si->width;
    int H = si->height;

    float* temp = new float[std::max(W,H)];
    float* avg  = new float[std::max(W,H)];
    float* dev  = new float[std::max(W,H)];

    memset (temp, 0, std::max(W,H)*sizeof(float));
    memset (avg,  0, std::max(W,H)*sizeof(float));
    memset (dev,  0, std::max(W,H)*sizeof(float));

    for (int i=row_from; i<row_to; i++) {
        for (int j=5; j<W-5; j++) {
            temp[j] = (si->raw[i][j-5] - 8*si->raw[i][j-4] + 27*si->raw[i][j-3] - 48*si->raw[i][j-2] + 42*si->raw[i][j-1] -
                    (si->raw[i][j+5] - 8*si->raw[i][j+4] + 27*si->raw[i][j+3] - 48*si->raw[i][j+2] + 42*si->raw[i][j+1])) / 100;
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
            if (hpmap->rows[i][j] < 0.8*hpv)
                *(int*)(&(hpmap->rows[i][j])) = 2;
            else if (hpv < 0.8*hpmap->rows[i][j])
                *(int*)(&(hpmap->rows[i][j])) = 1;
            else
                *(int*)(&(hpmap->rows[i][j])) = 0;
        }
    }

    delete [] temp;
    delete [] avg;
    delete [] dev;
}

void DemosaicFilter::hphd_green (MultiImage* si, MultiImage* ti, Buffer<int>* hpmap) {

    #pragma omp parallel for if (multiThread)
    for (int i=border; i<si->height-border; i++) {
        for (int j=border; j<si->width-border; j++) {
            if (si->raw_isGreen(i,j))
                ti->g[i-border][j-border] = si->raw[i][j];
            else {
                if (hpmap->rows[i][j]==1) {
                    int g2 = si->raw[i][j+1] + ((si->raw[i][j] - si->raw[i][j+2]) >> 1);
                    int g4 = si->raw[i][j-1] + ((si->raw[i][j] - si->raw[i][j-2]) >> 1);

                    int dx = si->raw[i][j+1] - si->raw[i][j-1];
                    int d1 = si->raw[i][j+3] - si->raw[i][j+1];
                    int d2 = si->raw[i][j+2] - si->raw[i][j];
                    int d3 = (si->raw[i-1][j+2] - si->raw[i-1][j]) >> 1;
                    int d4 = (si->raw[i+1][j+2] - si->raw[i+1][j]) >> 1;

                    double e2 = 1.0 / (1.0 + ABS(dx) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

                    d1 = si->raw[i][j-3] - si->raw[i][j-1];
                    d2 = si->raw[i][j-2] - si->raw[i][j];
                    d3 = (si->raw[i-1][j-2] - si->raw[i-1][j]) >> 1;
                    d4 = (si->raw[i+1][j-2] - si->raw[i+1][j]) >> 1;

                    double e4 = 1.0 / (1.0 + ABS(dx) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

                    ti->g[i-border][j-border] = CLIP((e2 * g2 + e4 * g4) / (e2 + e4));
                }
                else if (hpmap->rows[i][j]==2) {
                    int g1 = si->raw[i-1][j] + ((si->raw[i][j] - si->raw[i-2][j]) >> 1);
                    int g3 = si->raw[i+1][j] + ((si->raw[i][j] - si->raw[i+2][j]) >> 1);

                    int dy = si->raw[i+1][j] - si->raw[i-1][j];
                    int d1 = si->raw[i-1][j] - si->raw[i-3][j];
                    int d2 = si->raw[i][j] - si->raw[i-2][j];
                    int d3 = (si->raw[i][j-1] - si->raw[i-2][j-1]) >> 1;
                    int d4 = (si->raw[i][j+1] - si->raw[i-2][j+1]) >> 1;

                    double e1 = 1.0 / (1.0 + ABS(dy) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

                    d1 = si->raw[i+1][j] - si->raw[i+3][j];
                    d2 = si->raw[i][j] - si->raw[i+2][j];
                    d3 = (si->raw[i][j-1] - si->raw[i+2][j-1]) >> 1;
                    d4 = (si->raw[i][j+1] - si->raw[i+2][j+1]) >> 1;

                    double e3 = 1.0 / (1.0 + ABS(dy) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

                    ti->g[i-border][j-border] = CLIP((e1 * g1 + e3 * g3) / (e1 + e3));
                }
                else {
                    int g1 = si->raw[i-1][j] + ((si->raw[i][j] - si->raw[i-2][j]) >> 1);
                    int g2 = si->raw[i][j+1] + ((si->raw[i][j] - si->raw[i][j+2]) >> 1);
                    int g3 = si->raw[i+1][j] + ((si->raw[i][j] - si->raw[i+2][j]) >> 1);
                    int g4 = si->raw[i][j-1] + ((si->raw[i][j] - si->raw[i][j-2]) >> 1);

                    int dx = si->raw[i][j+1] - si->raw[i][j-1];
                    int dy = si->raw[i+1][j] - si->raw[i-1][j];

                    int d1 = si->raw[i-1][j] - si->raw[i-3][j];
                    int d2 = si->raw[i][j] - si->raw[i-2][j];
                    int d3 = (si->raw[i][j-1] - si->raw[i-2][j-1]) >> 1;
                    int d4 = (si->raw[i][j+1] - si->raw[i-2][j+1]) >> 1;

                    double e1 = 1.0 / (1.0 + ABS(dy) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

                    d1 = si->raw[i][j+3] - si->raw[i][j+1];
                    d2 = si->raw[i][j+2] - si->raw[i][j];
                    d3 = (si->raw[i-1][j+2] - si->raw[i-1][j]) >> 1;
                    d4 = (si->raw[i+1][j+2] - si->raw[i+1][j]) >> 1;

                    double e2 = 1.0 / (1.0 + ABS(dx) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

                    d1 = si->raw[i+1][j] - si->raw[i+3][j];
                    d2 = si->raw[i][j] - si->raw[i+2][j];
                    d3 = (si->raw[i][j-1] - si->raw[i+2][j-1]) >> 1;
                    d4 = (si->raw[i][j+1] - si->raw[i+2][j+1]) >> 1;

                    double e3 = 1.0 / (1.0 + ABS(dy) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

                    d1 = si->raw[i][j-3] - si->raw[i][j-1];
                    d2 = si->raw[i][j-2] - si->raw[i][j];
                    d3 = (si->raw[i-1][j-2] - si->raw[i-1][j]) >> 1;
                    d4 = (si->raw[i+1][j-2] - si->raw[i+1][j]) >> 1;

                    double e4 = 1.0 / (1.0 + ABS(dx) + ABS(d1) + ABS(d2) + ABS(d3) + ABS(d4));

                    ti->g[i-border][j-border] = CLIP((e1*g1 + e2*g2 + e3*g3 + e4*g4) / (e1 + e2 + e3 + e4));
                }
            }
        }
    }
}

void DemosaicFilter::hphd_demosaic (MultiImage* si, MultiImage* ti, Buffer<int>* hpmap) {

    if (getProgressListener()) {
        getProgressListener()->startTimeConsumingOperation ();
        getProgressListener()->setProgressStr ("Demosaicing...");
        getProgressListener()->setProgress (0.0);
    }

    memset(hpmap->data, 0, hpmap->width*hpmap->height*sizeof(float));

    #pragma omp parallel if (multiThread)
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        int blk = si->width/nthreads;

        if (tid<nthreads-1)
            hphd_vertical (si, (Buffer<float>*)hpmap, tid*blk, (tid+1)*blk);
        else
            hphd_vertical (si, (Buffer<float>*)hpmap, tid*blk, si->width);
    }

    if (getProgressListener())
        getProgressListener()->setProgress (0.33);

    #pragma omp parallel if (multiThread)
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        int blk = si->height/nthreads;

        if (tid<nthreads-1)
            hphd_horizontal (si, (Buffer<float>*)hpmap, tid*blk, (tid+1)*blk);
        else
            hphd_horizontal (si, (Buffer<float>*)hpmap, tid*blk, si->height);
    }

    if (getProgressListener())
        getProgressListener()->setProgress (0.66);

    hphd_green (si, ti, hpmap);

    if (getProgressListener()) {
        getProgressListener()->setProgress (1.0);
        getProgressListener()->progressReady ();
    }
}

void DemosaicFilter::interpolate_rb_bilinear (MultiImage* si, MultiImage* ti) {

    int W = si->width;
    int H = si->height;
    #pragma omp parallel for if (multiThread)
    for (int i=border; i<si->height-border; i++) {
        int ix = i-border;
        if (si->raw_isRed (i,0) || si->raw_isRed(i,1)) {
            // RGRGR or GRGRGR line
            for (int j=border, jx=0; j<si->width-border; j++, jx++) {
                if (si->raw_isRed (i,j)) {
                    // red is simple
                    ti->r[ix][jx] = si->raw[i][j];
                    // blue: cross interpolation
                    int b = 0;
                    int n = 0;
                    if (ix>0 && jx>0) {
                        b += si->raw[i-1][j-1] - ti->g[ix-1][jx-1];
                        n++;
                    }
                    if (ix>0 && jx<W-1) {
                        b += si->raw[i-1][j+1] - ti->g[ix-1][jx+1];
                        n++;
                    }
                    if (ix<H-1 && jx>0) {
                        b += si->raw[i+1][j-1] - ti->g[ix+1][jx-1];
                        n++;
                    }
                    if (ix<H-1 && jx<W-1) {
                        b += si->raw[i+1][j+1] - ti->g[ix+1][jx+1];
                        n++;
                    }
                    b = ti->g[ix][jx] + b / n;
                    ti->b[ix][jx] = CLIP(b);
                }
                else {
                    // linear R-G interp. horizontally
                    int r;
                    if (jx==0)
                        r = ti->g[ix][0] + si->raw[i][1] - ti->g[ix][1];
                    else if (jx==W-1)
                        r = ti->g[ix][W-1] + si->raw[i][W-2] - ti->g[ix][W-2];
                    else
                        r = ti->g[ix][jx] + (si->raw[i][j-1] - ti->g[ix][jx-1] + si->raw[i][j+1] - ti->g[ix][jx+1]) / 2;
                    ti->r[ix][jx] = CLIP(r);
                    // linear B-G interp. vertically
                    int b;
                    if (ix==0)
                        b = ti->g[ix+1][jx] + si->raw[1][j] - ti->g[ix][jx];
                    else if (ix==H-1)
                        b = ti->g[ix-1][jx] + si->raw[H-2][j] - ti->g[ix][jx];
                    else
                        b = ti->g[ix][jx] + (si->raw[i-1][j] - ti->g[ix-1][jx] + si->raw[i+1][j] - ti->g[ix+1][jx]) / 2;
                    ti->b[ix][jx] = CLIP(b);
                }
            }
        }
        else {
            // BGBGB or GBGBGB line
            for (int j=border, jx=0; j<si->width-border; j++, jx++) {
                if (si->raw_isBlue (i,j)) {
                    // blue is simple
                    ti->b[ix][jx] = si->raw[i][j];
                    // red: cross interpolation
                    int r = 0;
                    int n = 0;
                    if (ix>0 && jx>0) {
                        r += si->raw[i-1][j-1] - ti->g[ix-1][jx-1];
                        n++;
                    }
                    if (ix>0 && jx<W-1) {
                        r += si->raw[i-1][j+1] - ti->g[ix-1][jx+1];
                        n++;
                    }
                    if (ix<H-1 && jx>0) {
                        r += si->raw[i+1][j-1] - ti->g[ix+1][jx-1];
                        n++;
                    }
                    if (ix<H-1 && jx<W-1) {
                        r += si->raw[i+1][j+1] - ti->g[ix+1][jx+1];
                        n++;
                    }
                    r = ti->g[ix][jx]; + r / n;
                    ti->r[ix][jx] = CLIP(r);
                }
                else {
                    // linear B-G interp. horizontally
                    int b;
                    if (jx==0)
                        b = ti->g[ix][0] + si->raw[i][1] - ti->g[ix][1];
                    else if (jx==W-1)
                        b = ti->g[ix][W-1] + si->raw[i][W-2] - ti->g[ix][W-2];
                    else
                        b = ti->g[ix][jx] + (si->raw[i][j-1] - ti->g[ix][jx-1] + si->raw[i][j+1] - ti->g[ix][jx+1]) / 2;
                    ti->b[ix][jx] = CLIP(b);
                    // linear R-G interp. vertically
                    int r;
                    if (ix==0)
                        r = ti->g[ix+1][jx] + si->raw[1][j] - ti->g[ix][jx];
                    else if (ix==H-1)
                        r = ti->g[ix-1][jx] + si->raw[H-2][j] - ti->g[ix][jx];
                    else
                        r = ti->g[ix][jx] + (si->raw[i-1][j] - ti->g[ix-1][jx] + si->raw[i+1][j] - ti->g[ix+1][jx]) / 2;
                    ti->r[ix][jx] = CLIP(r);
                }
            }
        }
    }
}

void DemosaicFilter::correction_YIQ_LQ_  (MultiImage* im, int row_from, int row_to) {

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

  for (int i=0; i<3; i++)
      delete [] rbconv_Y[i];
  delete [] rbconv_Y;
  for (int i=0; i<3; i++)
      delete [] rbconv_I[i];
  delete [] rbconv_I;
  for (int i=0; i<3; i++)
      delete [] rbconv_Q[i];
  delete [] rbconv_Q;
  for (int i=0; i<3; i++)
      delete [] rbout_I[i];
  delete [] rbout_I;
  for (int i=0; i<3; i++)
      delete [] rbout_Q[i];
  delete [] rbout_Q;
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

void DemosaicFilter::convert_row_to_YIQ (unsigned short* r, unsigned short* g, unsigned short* b, int* Y, int* I, int* Q, int W) {
  for (int j=0; j<W; j++) {
    Y[j] = 299 * r[j] + 587 * g[j] + 114 * b[j];
    I[j] = 596 * r[j] - 275 * g[j] - 321 * b[j];
    Q[j] = 212 * r[j] - 523 * g[j] + 311 * b[j];
  }
}


void DemosaicFilter::convert_row_to_RGB (unsigned short* r, unsigned short* g, unsigned short* b, int* Y, int* I, int* Q, int W) {
  for (int j=1; j<W-1; j++) {
    int ir = Y[j]/1000 + 0.956*I[j]/1000 + 0.621*Q[j]/1000;
    int ig = Y[j]/1000 - 0.272*I[j]/1000 - 0.647*Q[j]/1000;
    int ib = Y[j]/1000 - 1.105*I[j]/1000 + 1.702*Q[j]/1000;
    r[j] = CLIP(ir);
    g[j] = CLIP(ig);
    b[j] = CLIP(ib);
  }
}

void DemosaicFilter::correction_YIQ_LQ  (MultiImage* im, int times) {

    if (im->height < 4)
        return;

    for (int t=0; t<times; t++) {
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
    }
}

void DemosaicFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

    if (procParams->demosaicing.method == "hphd")
        hphd_demosaic (sourceImage, targetImage, buffer);

    interpolate_rb_bilinear (sourceImage, targetImage);
    correction_YIQ_LQ (targetImage, procParams->demosaicing.colorcorrectionsteps);
}

}
