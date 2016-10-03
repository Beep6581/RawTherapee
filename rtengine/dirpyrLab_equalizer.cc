/*
 *  This file is part of RawTherapee.
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
 *
 *  (C) 2010 Emil Martinec <ejmartin@uchicago.edu>
 *
 */

#include <cstddef>
#include <cmath>
#include "curves.h"
#include "labimage.h"
#include "improcfun.h"
#include "rawimagesource.h"
#include "rt_math.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define CLIPC(a) ((a)>-32000?((a)<32000?(a):32000):-32000)

#define DIRWT(i1,j1,i,j) (rangefn[abs((int)data_fine->L[i1][j1]-data_fine->L[i][j])+abs((int)data_fine->a[i1][j1]-data_fine->a[i][j])+abs((int)data_fine->b[i1][j1]-data_fine->b[i][j])] )

namespace rtengine
{

static const int maxlevel = 4;

//sequence of scales
static const int scales[8] = {1, 2, 4, 8, 16, 32, 64, 128};
//sequence of pitches
static const int pitches[8] = {1, 1, 1, 1, 1, 1, 1, 1};

//sequence of scales
//static const int scales[8] = {1,1,1,1,1,1,1,1};
//sequence of pitches
//static const int pitches[8] = {2,2,2,2,2,2,2,2};

//sequence of scales
//static const int scales[8] = {1,3,6,10,15,21,28,36};
//sequence of pitches
//static const int pitches[8] = {1,1,1,1,1,1,1,1};

//sequence of scales
//static const int scales[8] = {1,1,2,4,8,16,32,64};
//sequence of pitches
//static const int pitches[8] = {2,1,1,1,1,1,1,1};

//pitch is spacing of subsampling
//scale is spacing of directional averaging weights
//example 1: no subsampling at any level -- pitch=1, scale=2^n
//example 2: subsampling by 2 every level -- pitch=2, scale=1 at each level
//example 3: no subsampling at first level, subsampling by 2 thereafter --
//  pitch =1, scale=1 at first level; pitch=2, scale=2 thereafter




void ImProcFunctions :: dirpyrLab_equalizer(LabImage * src, LabImage * dst, /*float luma, float chroma, float gamma*/ const double * mult )
{

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    LUTf rangefn(0x20000);


    //set up weights
    float noise = 1500;


    //set up range functions

    for (int i = 0; i < 0x20000; i++) {
        rangefn[i] = (int)((noise / ((double)i + noise)));
    }


    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    int level;
    int ** buffer[3];


    LabImage * dirpyrLablo[maxlevel];

    int w = src->W;
    int h = src->H;

    buffer[0] = allocArray<int> (w + 128, h + 128);
    buffer[1] = allocArray<int> (w + 128, h + 128);
    buffer[2] = allocArray<int> (w + 128, h + 128);

    for (int i = 0; i < h + 128; i++)
        for (int j = 0; j < w + 128; j++) {
            for (int c = 0; c < 3; c++) {
                buffer[c][i][j] = 0;
            }
        }

    w = (int)((w - 1) / pitches[0]) + 1;
    h = (int)((h - 1) / pitches[0]) + 1;

    dirpyrLablo[0] = new LabImage(w, h);

    for (level = 1; level < maxlevel; level++) {
        w = (int)((w - 1) / pitches[level]) + 1;
        h = (int)((h - 1) / pitches[level]) + 1;
        dirpyrLablo[level] = new LabImage(w, h);
    };

    //////////////////////////////////////////////////////////////////////////////


    // c[0] = luma = noise_L
    // c[1] = chroma = noise_ab
    // c[2] decrease of noise var with scale
    // c[3] radius of domain blur at each level
    // c[4] shadow smoothing
    // c[5] edge preservation

    level = 0;

    int scale = scales[level];

    int pitch = pitches[level];

    //int thresh = 10 * c[8];
    //impulse_nr (src, src, m_w1, m_h1, thresh, noisevar);

    dirpyr_eq(src, dirpyrLablo[0], rangefn, 0, pitch, scale, mult );

    level = 1;

    int totalpitch = pitches[0];

    while(level < maxlevel) {
        scale = scales[level];
        pitch = pitches[level];

        dirpyr_eq(dirpyrLablo[level - 1], dirpyrLablo[level], rangefn, level, pitch, scale, mult );

        level ++;
        totalpitch *= pitch;
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //initiate buffer for final image
    for(int i = 0, i1 = 0; i < src->H; i += totalpitch, i1++)
        for(int j = 0, j1 = 0; j < src->W; j += totalpitch, j1++) {

            //copy pixels
            buffer[0][i][j] = dirpyrLablo[maxlevel - 1]->L[i1][j1];
            buffer[1][i][j] = dirpyrLablo[maxlevel - 1]->a[i1][j1];
            buffer[2][i][j] = dirpyrLablo[maxlevel - 1]->b[i1][j1];

        }

    //if we are not subsampling, this is lots faster but does the typecasting work???
    //memcpy(buffer[0],dirpyrLablo[maxlevel-1]->L,sizeof(buffer[0]));
    //memcpy(buffer[1],dirpyrLablo[maxlevel-1]->a,sizeof(buffer[1]));
    //memcpy(buffer[2],dirpyrLablo[maxlevel-1]->b,sizeof(buffer[2]));

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    for(int level = maxlevel - 1; level > 0; level--) {

        //int scale = scales[level];
        int pitch = pitches[level];

        totalpitch /= pitch;

        idirpyr_eq(dirpyrLablo[level], dirpyrLablo[level - 1], buffer, level, pitch, totalpitch, mult );

    }


    scale = scales[0];
    pitch = pitches[0];
    totalpitch /= pitch;

    idirpyr_eq(dirpyrLablo[0], dst, buffer, 0, pitch, totalpitch, mult );


    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    for (int i = 0; i < dst->H; i++)
        for (int j = 0; j < dst->W; j++) {

            // TODO: Is integer cast necessary here?
            dst->L[i][j] = CLIP((int)(  buffer[0][i][j]  ));
            dst->a[i][j] = CLIPC((int)( buffer[1][i][j]  ));
            dst->b[i][j] = CLIPC((int)( buffer[2][i][j]  ));

        }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    for(int i = 0; i < maxlevel; i++) {
        delete dirpyrLablo[i];
    }

    for (int c = 0; c < 3; c++) {
        freeArray<int>(buffer[c], h + 128);
    }


    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
}

void ImProcFunctions::dirpyr_eq(LabImage* data_fine, LabImage* data_coarse, LUTf & rangefn, int level, int pitch, int scale, const double * mult  )
{

    //pitch is spacing of subsampling
    //scale is spacing of directional averaging weights
    //example 1: no subsampling at any level -- pitch=1, scale=2^n
    //example 2: subsampling by 2 every level -- pitch=2, scale=1 at each level
    //example 3: no subsampling at first level, subsampling by 2 thereafter --
    //  pitch =1, scale=1 at first level; pitch=2, scale=2 thereafter


    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // calculate weights, compute directionally weighted average

    int width = data_fine->W;
    int height = data_fine->H;



    //generate domain kernel
    int halfwin = 1;//min(ceil(2*sig),3);
    int scalewin = halfwin * scale;


#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for(int i = 0; i < height; i += pitch) {
        int i1 = i / pitch;

        for(int j = 0, j1 = 0; j < width; j += pitch, j1++) {
            float Lout, aout, bout;
            float norm;
            norm = 0;//if we do want to include the input pixel in the sum
            Lout = 0;
            aout = 0;
            bout = 0;

            for(int inbr = max(0, i - scalewin); inbr <= min(height - 1, i + scalewin); inbr += scale) {
                for (int jnbr = max(0, j - scalewin); jnbr <= min(width - 1, j + scalewin); jnbr += scale) {
                    float dirwt = DIRWT(inbr, jnbr, i, j);
                    Lout += dirwt * data_fine->L[inbr][jnbr];
                    aout += dirwt * data_fine->a[inbr][jnbr];
                    bout += dirwt * data_fine->b[inbr][jnbr];
                    norm += dirwt;
                }
            }

            data_coarse->L[i1][j1] = Lout / norm; //low pass filter
            data_coarse->a[i1][j1] = aout / norm;
            data_coarse->b[i1][j1] = bout / norm;
        }
    }




}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void ImProcFunctions::idirpyr_eq(LabImage* data_coarse, LabImage* data_fine, int *** buffer, int level, int pitch, int scale, const double * mult )
{

    int width = data_fine->W;
    int height = data_fine->H;

    float lumamult[4], chromamult[4];

    for (int i = 0; i < 4; i++) {
        lumamult[i] = mult[i];
        chromamult[i] = mult[i + 4];
    }

    float wtdsum[6], norm, dirwt;
    float hipass[3];
    int i1, j1;


    // for coarsest level, take non-subsampled lopass image and subtract from lopass_fine to generate hipass image

    // denoise hipass image, add back into lopass_fine to generate denoised image at fine scale

    // now iterate:
    // (1) take denoised image at level n, expand and smooth using gradient weights from lopass image at level n-1
    //     the result is the smoothed image at level n-1
    // (2) subtract smoothed image at level n-1 from lopass image at level n-1 to make hipass image at level n-1
    // (3) denoise the hipass image at level n-1
    // (4) add the denoised image at level n-1 to the smoothed image at level n-1 to make the denoised image at level n-1

    // note that the coarsest level amounts to skipping step (1) and doing (2,3,4).
    // in other words, skip step one if pitch=1



    if (pitch == 1) {
        // step (1-2-3-4)
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for(int i = 0; i < height; i++)
            for(int j = 0; j < width; j++) {

                //luma
                float hipass0 = (float)data_fine->L[i][j] - data_coarse->L[i][j];
                buffer[0][i * scale][j * scale] += hipass0 * lumamult[level]; //*luma;

                //chroma
                float hipass1 = data_fine->a[i][j] - data_coarse->a[i][j];
                float hipass2 = data_fine->b[i][j] - data_coarse->b[i][j];
                buffer[1][i * scale][j * scale] += hipass1 * chromamult[level]; //*chroma;
                buffer[2][i * scale][j * scale] += hipass2 * chromamult[level]; //*chroma;
            }

    } else {

        // step (1)
        //if (pitch>1), pitch=2; expand coarse image, fill in missing data

        LabImage* smooth;
        smooth = new LabImage(width, height);
#ifdef _OPENMP
        #pragma omp parallel
#endif

        {
#ifdef _OPENMP
            #pragma omp for
#endif

            for(int i = 0; i < height; i += pitch) {
                int i2 = i / pitch;

                for(int j = 0, j2 = 0; j < width; j += pitch, j2++) {

                    //copy common pixels
                    smooth->L[i][j] = data_coarse->L[i2][j2];
                    smooth->a[i][j] = data_coarse->a[i2][j2];
                    smooth->b[i][j] = data_coarse->b[i2][j2];
                }
            }

#ifdef _OPENMP
            #pragma omp for
#endif

            for(int i = 0; i < height - 1; i += 2)
                for(int j = 0; j < width - 1; j += 2) {
                    //do midpoint first
                    norm = dirwt = 0;
                    wtdsum[0] = wtdsum[1] = wtdsum[2] = wtdsum[3] = wtdsum[4] = wtdsum[5] = 0.0;

                    for(i1 = i; i1 < min(height, i + 3); i1 += 2)
                        for (j1 = j; j1 < min(width, j + 3); j1 += 2) {
                            dirwt = 1;//IDIRWT(i1, j1, i, j);
                            wtdsum[0] += dirwt * smooth->L[i1][j1];
                            wtdsum[1] += dirwt * smooth->a[i1][j1];
                            wtdsum[2] += dirwt * smooth->b[i1][j1];
                            wtdsum[3] += dirwt * buffer[0][i1 * scale][j1 * scale]; // not completely right if j1*scale or i1*scale is out of bounds of original image ???
                            wtdsum[4] += dirwt * buffer[1][i1 * scale][j1 * scale]; // also should we use directional average?
                            wtdsum[5] += dirwt * buffer[2][i1 * scale][j1 * scale];
                            norm += dirwt;
                        }

                    norm = 1 / norm;
                    smooth->L[i + 1][j + 1] = wtdsum[0] * norm;
                    smooth->a[i + 1][j + 1] = wtdsum[1] * norm;
                    smooth->b[i + 1][j + 1] = wtdsum[2] * norm;
                    buffer[0][(i + 1)*scale][(j + 1)*scale] = wtdsum[3] * norm;
                    buffer[1][(i + 1)*scale][(j + 1)*scale] = wtdsum[4] * norm;
                    buffer[2][(i + 1)*scale][(j + 1)*scale] = wtdsum[5] * norm;
                }

#ifdef _OPENMP
            #pragma omp for
#endif

            for(int i = 0; i < height - 1; i += 2)
                for(int j = 0; j < width - 1; j += 2) {
                    //now right neighbor
                    if (j + 1 == width) {
                        continue;
                    }

                    norm = dirwt = 0;
                    wtdsum[0] = wtdsum[1] = wtdsum[2] = wtdsum[3] = wtdsum[4] = wtdsum[5] = 0.0;

                    for (j1 = j; j1 < min(width, j + 3); j1 += 2) {
                        dirwt = 1;//IDIRWT(i, j1, i, j);
                        wtdsum[0] += dirwt * smooth->L[i][j1];
                        wtdsum[1] += dirwt * smooth->a[i][j1];
                        wtdsum[2] += dirwt * smooth->b[i][j1];
                        wtdsum[3] += dirwt * buffer[0][i * scale][j1 * scale];
                        wtdsum[4] += dirwt * buffer[1][i * scale][j1 * scale];
                        wtdsum[5] += dirwt * buffer[2][i * scale][j1 * scale];
                        norm += dirwt;
                    }

                    for (i1 = max(0, i - 1); i1 < min(height, i + 2); i1 += 2) {
                        dirwt = 1;//IDIRWT(i1, j+1, i, j);
                        wtdsum[0] += dirwt * smooth->L[i1][j + 1];
                        wtdsum[1] += dirwt * smooth->a[i1][j + 1];
                        wtdsum[2] += dirwt * smooth->b[i1][j + 1];
                        wtdsum[3] += dirwt * buffer[0][i1 * scale][(j + 1) * scale];
                        wtdsum[4] += dirwt * buffer[1][i1 * scale][(j + 1) * scale];
                        wtdsum[5] += dirwt * buffer[2][i1 * scale][(j + 1) * scale];
                        norm += dirwt;
                    }

                    norm = 1 / norm;
                    smooth->L[i][j + 1] = wtdsum[0] * norm;
                    smooth->a[i][j + 1] = wtdsum[1] * norm;
                    smooth->b[i][j + 1] = wtdsum[2] * norm;
                    buffer[0][i][(j + 1)*scale] = wtdsum[3] * norm;
                    buffer[1][i][(j + 1)*scale] = wtdsum[4] * norm;
                    buffer[2][i][(j + 1)*scale] = wtdsum[5] * norm;

                    //now down neighbor
                    if (i + 1 == height) {
                        continue;
                    }

                    norm = 0;
                    wtdsum[0] = wtdsum[1] = wtdsum[2] = wtdsum[3] = wtdsum[4] = wtdsum[5] = 0.0;

                    for (i1 = i; i1 < min(height, i + 3); i1 += 2) {
                        dirwt = 1;//IDIRWT(i1, j, i, j);
                        wtdsum[0] += dirwt * smooth->L[i1][j];
                        wtdsum[1] += dirwt * smooth->a[i1][j];
                        wtdsum[2] += dirwt * smooth->b[i1][j];
                        wtdsum[3] += dirwt * buffer[0][i1 * scale][j * scale];
                        wtdsum[4] += dirwt * buffer[1][i1 * scale][j * scale];
                        wtdsum[5] += dirwt * buffer[2][i1 * scale][j * scale];
                        norm += dirwt;
                    }

                    for (j1 = max(0, j - 1); j1 < min(width, j + 2); j1 += 2) {
                        dirwt = 1;//IDIRWT(i+1, j1, i, j);
                        wtdsum[0] += dirwt * smooth->L[i + 1][j1];
                        wtdsum[1] += dirwt * smooth->a[i + 1][j1];
                        wtdsum[2] += dirwt * smooth->b[i + 1][j1];
                        wtdsum[3] += dirwt * buffer[0][(i + 1) * scale][j1 * scale];
                        wtdsum[4] += dirwt * buffer[1][(i + 1) * scale][j1 * scale];
                        wtdsum[5] += dirwt * buffer[2][(i + 1) * scale][j1 * scale];
                        norm += dirwt;
                    }

                    norm = 1 / norm;
                    smooth->L[i + 1][j] = wtdsum[0] * norm;
                    smooth->a[i + 1][j] = wtdsum[1] * norm;
                    smooth->b[i + 1][j] = wtdsum[2] * norm;
                    buffer[0][(i + 1)*scale][j * scale] = wtdsum[3] * norm;
                    buffer[1][(i + 1)*scale][j * scale] = wtdsum[4] * norm;
                    buffer[2][(i + 1)*scale][j * scale] = wtdsum[5] * norm;

                }


            // step (2-3-4)
#ifdef _OPENMP
            #pragma omp for
#endif

            for(int i = 0; i < height; i++)
                for(int j = 0; j < width; j++) {

                    //luma
                    hipass[0] = (float)data_fine->L[i][j] - smooth->L[i][j];
                    buffer[0][i * scale][j * scale] += hipass[0] * lumamult[level]; //*luma;

                    //chroma
                    hipass[1] = data_fine->a[i][j] - smooth->a[i][j];
                    hipass[2] = data_fine->b[i][j] - smooth->b[i][j];
                    buffer[1][i * scale][j * scale] += hipass[1] * chromamult[level]; //*chroma;
                    buffer[2][i * scale][j * scale] += hipass[2] * chromamult[level]; //*chroma;
                }
        }   // end parallel
        delete smooth;

    }
}


#undef DIRWT_L
#undef DIRWT_AB

#undef NRWT_L
#undef NRWT_AB

}

