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
 *  © 2010 Emil Martinec <ejmartin@uchicago.edu>
 *
 */

#include <cstddef>
#include <cmath>
#include "curves.h"
#include "labimage.h"
#include "color.h"
#include "mytime.h"
#include "improcfun.h"
#include "rawimagesource.h"
#include "array2D.h"
#include "rt_math.h"
#include "opthelper.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#define CLIPI(a) ((a)>0 ?((a)<32768 ?(a):32768):0)

#define RANGEFN(i) ((1000.0f / (i + 1000.0f)))
#define CLIPC(a) ((a)>-32000?((a)<32000?(a):32000):-32000)
#define DIRWT(i1,j1,i,j) ( domker[(i1-i)/scale+halfwin][(j1-j)/scale+halfwin] * RANGEFN(fabsf((data_fine[i1][j1]-data_fine[i][j]))) )

namespace rtengine
{

static const int maxlevel = 6;
static const float noise = 2000;

//sequence of scales
static const int scales[6] = {1, 2, 4, 8, 16, 32};
extern const Settings* settings;

//sequence of scales


SSEFUNCTION void ImProcFunctions :: dirpyr_equalizer(float ** src, float ** dst, int srcwidth, int srcheight, float ** l_a, float ** l_b, float ** dest_a, float ** dest_b, const double * mult, const double dirpyrThreshold, const double skinprot, const bool gamutlab, float b_l, float t_l, float t_r, float b_r, int choice, int scaleprev)
{
    int lastlevel = maxlevel;

    if(settings->verbose) {
        printf("Dirpyr scaleprev=%i\n", scaleprev);
    }

    float atten123 = (float) settings->level123_cbdl;

    if(atten123 > 50.f) {
        atten123 = 50.f;
    }

    if(atten123 < 0.f) {
        atten123 = 0.f;
    }

    float atten0 = (float) settings->level0_cbdl;

    if(atten0 > 40.f) {
        atten123 = 40.f;
    }

    if(atten0 < 0.f) {
        atten0 = 0.f;
    }

    if((t_r - t_l) < 0.55f) {
        t_l = t_r + 0.55f;    //avoid too small range
    }


    while (lastlevel > 0 && fabs(mult[lastlevel - 1] - 1) < 0.001) {
        lastlevel--;
        //printf("last level to process %d \n",lastlevel);
    }

    if (lastlevel == 0) {
        return;
    }

    int level;
    float multi[6] = {1.f, 1.f, 1.f, 1.f, 1.f, 1.f};
    float scalefl[6];

    for(int lv = 0; lv < 6; lv++) {
        scalefl[lv] = ((float) scales[lv]) / (float) scaleprev;

        if(lv >= 1) {
            if(scalefl[lv] < 1.f) {
                multi[lv] = (atten123 * ((float) mult[lv] - 1.f) / 100.f) + 1.f;    //modulate action if zoom < 100%
            } else {
                multi[lv] = (float) mult[lv];
            }
        } else  {
            if(scalefl[lv] < 1.f) {
                multi[lv] = (atten0 * ((float) mult[lv] - 1.f) / 100.f) + 1.f;    //modulate action if zoom < 100%
            } else {
                multi[lv] = (float) mult[lv];
            }
        }

    }

    if(settings->verbose) {
        printf("CbDL mult0=%f  1=%f 2=%f 3=%f 4=%f 5=%f\n", multi[0], multi[1], multi[2], multi[3], multi[4], multi[5]);
    }

    multi_array2D<float, maxlevel> dirpyrlo (srcwidth, srcheight);

    level = 0;

    //int thresh = 100 * mult[5];
    int scale = (int)(scales[level]) / scaleprev;

    if(scale < 1) {
        scale = 1;
    }


    dirpyr_channel(src, dirpyrlo[0], srcwidth, srcheight, 0, scale);

    level = 1;

    while(level < lastlevel) {

        scale = (int)(scales[level]) / scaleprev;

        if(scale < 1) {
            scale = 1;
        }

        dirpyr_channel(dirpyrlo[level - 1], dirpyrlo[level], srcwidth, srcheight, level, scale);

        level ++;
    }

    float **tmpHue, **tmpChr;

    if(skinprot != 0.f) {
        // precalculate hue and chroma, use SSE, if available
        // by precalculating these values we can greatly reduce the number of calculations in idirpyr_eq_channel()
        // but we need two additional buffers for this preprocessing
        tmpHue = new float*[srcheight];

        for (int i = 0; i < srcheight; i++) {
            tmpHue[i] = new float[srcwidth];
        }

#ifdef __SSE2__
        #pragma omp parallel for

        for(int i = 0; i < srcheight; i++) {
            int j;

            for(j = 0; j < srcwidth - 3; j += 4) {
                _mm_storeu_ps(&tmpHue[i][j], xatan2f(LVFU(l_b[i][j]), LVFU(l_a[i][j])));
            }

            for(; j < srcwidth; j++) {
                tmpHue[i][j] = xatan2f(l_b[i][j], l_a[i][j]);
            }
        }

#else
        #pragma omp parallel for

        for(int i = 0; i < srcheight; i++) {
            for(int j = 0; j < srcwidth; j++) {
                tmpHue[i][j] = xatan2f(l_b[i][j], l_a[i][j]);
            }
        }

#endif
        tmpChr = new float*[srcheight];

        for (int i = 0; i < srcheight; i++) {
            tmpChr[i] = new float[srcwidth];
        }

#ifdef __SSE2__
        #pragma omp parallel
        {
            __m128 div = _mm_set1_ps(327.68f);
            #pragma omp for

            for(int i = 0; i < srcheight; i++) {
                int j;

                for(j = 0; j < srcwidth - 3; j += 4) {
                    _mm_storeu_ps(&tmpChr[i][j], _mm_sqrt_ps(SQRV(LVFU(l_b[i][j])) + SQRV(LVFU(l_a[i][j]))) / div);
                }

                for(; j < srcwidth; j++) {
                    tmpChr[i][j] = sqrtf(SQR((l_b[i][j])) + SQR((l_a[i][j]))) / 327.68f;
                }
            }
        }
#else
        #pragma omp parallel for

        for(int i = 0; i < srcheight; i++) {
            for(int j = 0; j < srcwidth; j++) {
                tmpChr[i][j] = sqrtf(SQR((l_b[i][j])) + SQR((l_a[i][j]))) / 327.68f;
            }
        }

#endif
    }

    // with the current implementation of idirpyr_eq_channel we can safely use the buffer from last level as buffer, saves some memory
    float ** buffer = dirpyrlo[lastlevel - 1];

    for(int level = lastlevel - 1; level > 0; level--) {
        idirpyr_eq_channel(dirpyrlo[level], dirpyrlo[level - 1], buffer, srcwidth, srcheight, level, multi, dirpyrThreshold, tmpHue, tmpChr, skinprot, gamutlab, b_l, t_l, t_r, b_r, choice );
    }

    scale = scales[0];

    idirpyr_eq_channel(dirpyrlo[0], dst, buffer, srcwidth, srcheight, 0, multi, dirpyrThreshold, tmpHue, tmpChr, skinprot, gamutlab, b_l, t_l, t_r, b_r, choice );

    if(skinprot != 0.f) {
        for (int i = 0; i < srcheight; i++) {
            delete [] tmpChr[i];
        }

        delete [] tmpChr;

        for (int i = 0; i < srcheight; i++) {
            delete [] tmpHue[i];
        }

        delete [] tmpHue;
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #pragma omp parallel for

    for (int i = 0; i < srcheight; i++)
        for (int j = 0; j < srcwidth; j++) {
            dst[i][j] = CLIP(buffer[i][j]);  // TODO: Really a clip necessary?
        }

}



void ImProcFunctions :: dirpyr_equalizercam (CieImage *ncie, float ** src, float ** dst, int srcwidth, int srcheight, float ** h_p, float ** C_p, const double * mult, const double dirpyrThreshold, const double skinprot, bool execdir,  const bool gamutlab, float b_l, float t_l, float t_r, float b_r, int choice, int scaleprev)
{
    int lastlevel = maxlevel;

    if(settings->verbose) {
        printf("CAM dirpyr scaleprev=%i\n", scaleprev);
    }

    float atten123 = (float) settings->level123_cbdl;

    if(atten123 > 50.f) {
        atten123 = 50.f;
    }

    if(atten123 < 0.f) {
        atten123 = 0.f;
    }

//      printf("atten=%f\n",atten);
    float atten0 = (float) settings->level0_cbdl;

    if(atten0 > 40.f) {
        atten123 = 40.f;
    }

    if(atten0 < 0.f) {
        atten0 = 0.f;
    }

    if((t_r - t_l) < 0.55f) {
        t_l = t_r + 0.55f;    //avoid too small range
    }

    while (fabs(mult[lastlevel - 1] - 1) < 0.001 && lastlevel > 0) {
        lastlevel--;
        //printf("last level to process %d \n",lastlevel);
    }

    if (lastlevel == 0) {
        return;
    }

    int level;

    float multi[6] = {1.f, 1.f, 1.f, 1.f, 1.f, 1.f};
    float scalefl[6];

    for(int lv = 0; lv < 6; lv++) {
        scalefl[lv] = ((float) scales[lv]) / (float) scaleprev;

        //  if(scalefl[lv] < 1.f) multi[lv] = 1.f; else  multi[lv]=(float) mult[lv];
        if (lv >= 1) {
            if(scalefl[lv] < 1.f) {
                multi[lv] = (atten123 * ((float) mult[lv] - 1.f) / 100.f) + 1.f;
            } else {
                multi[lv] = (float) mult[lv];
            }
        } else {
            if(scalefl[lv] < 1.f) {
                multi[lv] = (atten0 * ((float) mult[lv] - 1.f) / 100.f) + 1.f;
            } else {
                multi[lv] = (float) mult[lv];
            }
        }


    }

    if(settings->verbose) {
        printf("CAM CbDL mult0=%f  1=%f 2=%f 3=%f 4=%f 5=%f\n", multi[0], multi[1], multi[2], multi[3], multi[4], multi[5]);
    }




    multi_array2D<float, maxlevel> dirpyrlo (srcwidth, srcheight);

    level = 0;

    int scale = (int)(scales[level]) / scaleprev;

    if(scale < 1) {
        scale = 1;
    }

    dirpyr_channel(src, dirpyrlo[0], srcwidth, srcheight, 0, scale);

    level = 1;

    while(level < lastlevel) {
        scale = (int)(scales[level]) / scaleprev;

        if(scale < 1) {
            scale = 1;
        }

        dirpyr_channel(dirpyrlo[level - 1], dirpyrlo[level], srcwidth, srcheight, level, scale);

        level ++;
    }


    // with the current implementation of idirpyr_eq_channel we can safely use the buffer from last level as buffer, saves some memory
    float ** buffer = dirpyrlo[lastlevel - 1];

    for(int level = lastlevel - 1; level > 0; level--) {
        idirpyr_eq_channelcam(dirpyrlo[level], dirpyrlo[level - 1], buffer, srcwidth, srcheight, level, multi, dirpyrThreshold , h_p, C_p, skinprot, b_l, t_l, t_r);
    }


    scale = scales[0];

    idirpyr_eq_channelcam(dirpyrlo[0], dst, buffer, srcwidth, srcheight, 0, multi, dirpyrThreshold,  h_p, C_p, skinprot, b_l, t_l, t_r);


    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(execdir) {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

        for (int i = 0; i < srcheight; i++)
            for (int j = 0; j < srcwidth; j++) {
                if(ncie->J_p[i][j] > 8.f && ncie->J_p[i][j] < 92.f) {
                    dst[i][j] = CLIP( buffer[i][j] );    // TODO: Really a clip necessary?
                } else {
                    dst[i][j] = src[i][j];
                }
            }
    } else
        for (int i = 0; i < srcheight; i++)
            for (int j = 0; j < srcwidth; j++) {
                dst[i][j] = CLIP( buffer[i][j] );  // TODO: Really a clip necessary?
            }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
}


SSEFUNCTION void ImProcFunctions::dirpyr_channel(float ** data_fine, float ** data_coarse, int width, int height, int level, int scale)
{
    //scale is spacing of directional averaging weights
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // calculate weights, compute directionally weighted average

    if(level > 1) {
        //generate domain kernel
        int domker[5][5] = {{1, 1, 1, 1, 1}, {1, 2, 2, 2, 1}, {1, 2, 2, 2, 1}, {1, 2, 2, 2, 1}, {1, 1, 1, 1, 1}};
        //  int domker[5][5] = {{1,1,1,1,1},{1,1,1,1,1},{1,1,1,1,1},{1,1,1,1,1},{1,1,1,1,1}};
        static const int halfwin = 2;
        const int scalewin = halfwin * scale;
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
#ifdef __SSE2__
            __m128 thousandv = _mm_set1_ps( 1000.0f );
            __m128 dirwtv, valv, normv, dftemp1v, dftemp2v;
//  multiplied each value of domkerv by 1000 to avoid multiplication by 1000 inside the loop
            float domkerv[5][5][4] ALIGNED16 = {{{1000, 1000, 1000, 1000}, {1000, 1000, 1000, 1000}, {1000, 1000, 1000, 1000}, {1000, 1000, 1000, 1000}, {1000, 1000, 1000, 1000}}, {{1000, 1000, 1000, 1000}, {2000, 2000, 2000, 2000}, {2000, 2000, 2000, 2000}, {2000, 2000, 2000, 2000}, {1000, 1000, 1000, 1000}}, {{1000, 1000, 1000, 1000}, {2000, 2000, 2000, 2000}, {2000, 2000, 2000, 2000}, {2000, 2000, 2000, 2000}, {1000, 1000, 1000, 1000}}, {{1000, 1000, 1000, 1000}, {2000, 2000, 2000, 2000}, {2000, 2000, 2000, 2000}, {2000, 2000, 2000, 2000}, {1000, 1000, 1000, 1000}}, {{1000, 1000, 1000, 1000}, {1000, 1000, 1000, 1000}, {1000, 1000, 1000, 1000}, {1000, 1000, 1000, 1000}, {1000, 1000, 1000, 1000}}};
#endif // __SSE2__

            int j;
#ifdef _OPENMP
            #pragma omp for //schedule (dynamic,8)
#endif

            for(int i = 0; i < height; i++) {
                float dirwt;

                for(j = 0; j < scalewin; j++) {
                    float val = 0.f;
                    float norm = 0.f;


                    for(int inbr = max(0, i - scalewin); inbr <= min(height - 1, i + scalewin); inbr += scale) {
                        for (int jnbr = max(0, j - scalewin); jnbr <= j + scalewin; jnbr += scale) {
                            //printf("i=%d ",(inbr-i)/scale+halfwin);
                            dirwt = DIRWT(inbr, jnbr, i, j);
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }

                    data_coarse[i][j] = val / norm; //low pass filter
                }

#ifdef __SSE2__

                for(; j < width - scalewin - 3; j += 4) {
                    valv = _mm_setzero_ps();
                    normv = _mm_setzero_ps();
                    dftemp1v = LVFU(data_fine[i][j]);

                    for(int inbr = MAX(0, i - scalewin); inbr <= MIN(height - 1, i + scalewin); inbr += scale) {
                        int indexihlp = (inbr - i) / scale + halfwin;

                        for (int jnbr = j - scalewin, indexjhlp = 0; jnbr <= j + scalewin; jnbr += scale, indexjhlp++) {
                            dftemp2v = LVFU(data_fine[inbr][jnbr]);
                            dirwtv = LVF(domkerv[indexihlp][indexjhlp]) / (vabsf(dftemp1v - dftemp2v) + thousandv);
                            valv += dirwtv * dftemp2v;
                            normv += dirwtv;
                        }
                    }

                    _mm_storeu_ps( &data_coarse[i][j], valv / normv); //low pass filter
                }

                for(; j < width - scalewin; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for(int inbr = max(0, i - scalewin); inbr <= min(height - 1, i + scalewin); inbr += scale) {
                        for (int jnbr = j - scalewin; jnbr <= j + scalewin; jnbr += scale) {
                            dirwt = DIRWT(inbr, jnbr, i, j);
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }

                    data_coarse[i][j] = val / norm; //low pass filter
                }

#else

                for(; j < width - scalewin; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for(int inbr = max(0, i - scalewin); inbr <= min(height - 1, i + scalewin); inbr += scale) {
                        for (int jnbr = j - scalewin; jnbr <= j + scalewin; jnbr += scale) {
                            dirwt = DIRWT(inbr, jnbr, i, j);
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }

                    data_coarse[i][j] = val / norm; //low pass filter
                }

#endif

                for(; j < width; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for(int inbr = max(0, i - scalewin); inbr <= min(height - 1, i + scalewin); inbr += scale) {
                        for (int jnbr = j - scalewin; jnbr <= min(width - 1, j + scalewin); jnbr += scale) {
                            dirwt = DIRWT(inbr, jnbr, i, j);
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }

                    data_coarse[i][j] = val / norm; //low pass filter
                }
            }
        }
    } else {    // level <=1 means that all values of domker would be 1.0f, so no need for multiplication
//      const int scalewin = scale;
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
#ifdef __SSE2__
            __m128 thousandv = _mm_set1_ps( 1000.0f );
            __m128 dirwtv, valv, normv, dftemp1v, dftemp2v;
#endif // __SSE2__
            int j;
#ifdef _OPENMP
            #pragma omp for schedule(dynamic,16)
#endif

            for(int i = 0; i < height; i++)
            {
                float dirwt;

                for(j = 0; j < scale; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for(int inbr = max(0, i - scale); inbr <= min(height - 1, i + scale); inbr += scale) {
                        for (int jnbr = max(0, j - scale); jnbr <= j + scale; jnbr += scale) {
                            dirwt = RANGEFN(fabsf(data_fine[inbr][jnbr] - data_fine[i][j]));
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }

                    data_coarse[i][j] = val / norm; //low pass filter
                }

#ifdef __SSE2__

                for(; j < width - scale - 3; j += 4) {
                    valv = _mm_setzero_ps();
                    normv = _mm_setzero_ps();
                    dftemp1v = LVFU(data_fine[i][j]);

                    for(int inbr = MAX(0, i - scale); inbr <= MIN(height - 1, i + scale); inbr += scale) {
                        for (int jnbr = j - scale; jnbr <= j + scale; jnbr += scale) {
                            dftemp2v = LVFU(data_fine[inbr][jnbr]);
                            dirwtv = thousandv / (vabsf(dftemp2v - dftemp1v) + thousandv);
                            valv += dirwtv * dftemp2v;
                            normv += dirwtv;
                        }
                    }

                    _mm_storeu_ps( &data_coarse[i][j], valv / normv); //low pass filter
                }

                for(; j < width - scale; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for(int inbr = max(0, i - scale); inbr <= min(height - 1, i + scale); inbr += scale) {
                        for (int jnbr = j - scale; jnbr <= j + scale; jnbr += scale) {
                            dirwt = RANGEFN(fabsf(data_fine[inbr][jnbr] - data_fine[i][j]));
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }

                    data_coarse[i][j] = val / norm; //low pass filter
                }

#else

                for(; j < width - scale; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for(int inbr = max(0, i - scale); inbr <= min(height - 1, i + scale); inbr += scale) {
                        for (int jnbr = j - scale; jnbr <= j + scale; jnbr += scale) {
                            dirwt = RANGEFN(fabsf(data_fine[inbr][jnbr] - data_fine[i][j]));
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }

                    data_coarse[i][j] = val / norm; //low pass filter
                }

#endif

                for(; j < width; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for(int inbr = max(0, i - scale); inbr <= min(height - 1, i + scale); inbr += scale) {
                        for (int jnbr = j - scale; jnbr <= min(width - 1, j + scale); jnbr += scale) {
                            dirwt = RANGEFN(fabsf(data_fine[inbr][jnbr] - data_fine[i][j]));
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }

                    data_coarse[i][j] = val / norm; //low pass filter
                }
            }
        }
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void ImProcFunctions::idirpyr_eq_channel(float ** data_coarse, float ** data_fine, float ** buffer, int width, int height, int level, float mult[5], const double dirpyrThreshold, float ** hue, float ** chrom, const double skinprot, const bool gamutlab, float b_l, float t_l, float t_r, float b_r , int choice)
{
    const float skinprotneg = -skinprot;
    const float factorHard = (1.f - skinprotneg / 100.f);

    float offs;

    if(skinprot == 0.f) {
        offs = 0.f;
    } else {
        offs = -1.f;
    }

    float multbis[6];

    multbis[level] = mult[level]; //multbis to reduce artifacts for high values mult

    if(level == 4 && mult[level] > 1.f) {
        multbis[level] = 1.f + 0.65f * (mult[level] - 1.f);
    }

    if(level == 5 && mult[level] > 1.f) {
        multbis[level] = 1.f + 0.45f * (mult[level] - 1.f);
    }

    LUTf irangefn (0x20000);
    {
        const float noisehi = 1.33f * noise * dirpyrThreshold / expf(level * log(3.0)), noiselo = 0.66f * noise * dirpyrThreshold / expf(level * log(3.0));
        //printf("level=%i multlev=%f noisehi=%f noiselo=%f skinprot=%f\n",level,mult[level], noisehi, noiselo, skinprot);

        for (int i = 0; i < 0x20000; i++) {
            if (abs(i - 0x10000) > noisehi || multbis[level] < 1.0) {
                irangefn[i] = multbis[level] + offs;
            } else {
                if (abs(i - 0x10000) < noiselo) {
                    irangefn[i] = 1.f + offs ;
                } else {
                    irangefn[i] = 1.f + offs + (multbis[level] - 1.f) * (noisehi - abs(i - 0x10000)) / (noisehi - noiselo + 0.01f) ;
                }
            }
        }
    }

    if(skinprot == 0.f)
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif
        for(int i = 0; i < height; i++) {
            for(int j = 0; j < width; j++) {
                float hipass = (data_fine[i][j] - data_coarse[i][j]);
                buffer[i][j] += irangefn[hipass + 0x10000] * hipass;
            }
        }
    else if(skinprot > 0.f)
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif
        for(int i = 0; i < height; i++) {
            for(int j = 0; j < width; j++) {
                float scale = 1.f;
                float hipass = (data_fine[i][j] - data_coarse[i][j]);
                // These values are precalculated now
                float modhue = hue[i][j];
                float modchro = chrom[i][j];
                Color::SkinSatCbdl ((data_fine[i][j]) / 327.68f, modhue, modchro, skinprot, scale, true, b_l, t_l, t_r);
                buffer[i][j] += (1.f + (irangefn[hipass + 0x10000]) * scale) * hipass ;
            }
        }
    else
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif
        for(int i = 0; i < height; i++) {
            for(int j = 0; j < width; j++) {
                float scale = 1.f;
                float hipass = (data_fine[i][j] - data_coarse[i][j]);
                // These values are precalculated now
                float modhue = hue[i][j];
                float modchro = chrom[i][j];
                Color::SkinSatCbdl ((data_fine[i][j]) / 327.68f, modhue, modchro, skinprotneg, scale, false, b_l, t_l, t_r);
                float correct = irangefn[hipass + 0x10000];

                if (scale == 1.f) {//image hard
                    buffer[i][j] += (1.f + (correct) * (factorHard)) * hipass ;
                } else { //image soft with scale < 1 ==> skin
                    buffer[i][j] += (1.f + (correct)) * hipass ;
                }
            }
        }
}


void ImProcFunctions::idirpyr_eq_channelcam(float ** data_coarse, float ** data_fine, float ** buffer, int width, int height, int level, float mult[5], const double dirpyrThreshold, float ** l_a_h, float ** l_b_c, const double skinprot, float b_l, float t_l, float t_r)
{

    const float skinprotneg = -skinprot;
    const float factorHard = (1.f - skinprotneg / 100.f);

    float offs;

    if(skinprot == 0.f) {
        offs = 0.f;
    } else {
        offs = -1.f;
    }

    float multbis[6];

    multbis[level] = mult[level]; //multbis to reduce artifacts for high values mult

    if(level == 4 && mult[level] > 1.f) {
        multbis[level] = 1.f + 0.65f * (mult[level] - 1.f);
    }

    if(level == 5 && mult[level] > 1.f) {
        multbis[level] = 1.f + 0.45f * (mult[level] - 1.f);
    }

    LUTf irangefn (0x20000);
    {
        const float noisehi = 1.33f * noise * dirpyrThreshold / expf(level * log(3.0)), noiselo = 0.66f * noise * dirpyrThreshold / expf(level * log(3.0));

        //printf("level=%i multlev=%f noisehi=%f noiselo=%f skinprot=%f\n",level,mult[level], noisehi, noiselo, skinprot);
        for (int i = 0; i < 0x20000; i++) {
            if (abs(i - 0x10000) > noisehi || multbis[level] < 1.0) {
                irangefn[i] = multbis[level] + offs;
            } else {
                if (abs(i - 0x10000) < noiselo) {
                    irangefn[i] = 1.f + offs ;
                } else {
                    irangefn[i] = 1.f + offs + (multbis[level] - 1.f) * (noisehi - abs(i - 0x10000)) / (noisehi - noiselo + 0.01f) ;
                }
            }
        }
    }

    if(skinprot == 0.f)
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif
        for(int i = 0; i < height; i++) {
            for(int j = 0; j < width; j++) {
                float hipass = (data_fine[i][j] - data_coarse[i][j]);
                buffer[i][j] += irangefn[hipass + 0x10000] * hipass ;
            }
        }
    else if(skinprot > 0.f)
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif
        for(int i = 0; i < height; i++) {
            for(int j = 0; j < width; j++) {
                float hipass = (data_fine[i][j] - data_coarse[i][j]);
                float scale = 1.f;
                Color::SkinSatCbdlCam ((data_fine[i][j]) / 327.68f, l_a_h[i][j] , l_b_c[i][j], skinprot, scale, true, b_l, t_l, t_r);
                buffer[i][j] += (1.f + (irangefn[hipass + 0x10000]) * scale) * hipass ;
            }
        }
    else
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif
        for(int i = 0; i < height; i++) {
            for(int j = 0; j < width; j++) {
                float hipass = (data_fine[i][j] - data_coarse[i][j]);
                float scale = 1.f;
                float correct;
                correct = irangefn[hipass + 0x10000];
                Color::SkinSatCbdlCam ((data_fine[i][j]) / 327.68f, l_a_h[i][j], l_b_c[i][j] , skinprotneg, scale, false, b_l, t_l, t_r);

                if (scale == 1.f) {//image hard
                    buffer[i][j] += (1.f + (correct) * factorHard) * hipass ;

                } else { //image soft
                    buffer[i][j] += (1.f + (correct)) * hipass ;
                }
            }
        }

    //  if(gamutlab) {
    //    ImProcFunctions::badpixcam (buffer[i][j], 6.0, 10, 2);//for bad pixels
    //  }

    /*      if(gamutlab) {//disabled
            float Lprov1=(buffer[i][j])/327.68f;
            float R,G,B;
    #ifdef _DEBUG
                bool neg=false;
                bool more_rgb=false;
                //gamut control : Lab values are in gamut
                Color::gamutLchonly(modhue,Lprov1,modchro, R, G, B, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
    #else
                //gamut control : Lab values are in gamut
                Color::gamutLchonly(modhue,Lprov1,modchro, R, G, B, wip, highlight, 0.15f, 0.96f);
    #endif
    //      Color::gamutLchonly(modhue,Lprov1,modchro, R, G, B, wip, highlight, 0.15f, 0.96f);//gamut control in Lab mode ..not in CIECAM
                buffer[i][j]=Lprov1*327.68f;
                float2 sincosval = xsincosf(modhue);
                l_a_h[i][j]=327.68f*modchro*sincosval.y;
                l_b_c[i][j]=327.68f*modchro*sincosval.x;
            }
            */
}

//  float hipass = (data_fine[i][j]-data_coarse[i][j]);
//  buffer[i][j] += irangefn[hipass+0x10000] * hipass ;

#undef DIRWT_L
#undef DIRWT_AB

#undef NRWT_L
#undef NRWT_AB

}

