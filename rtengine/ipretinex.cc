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

    *   adaptation to RawTherapee
    *   2015 Jacques Desmis <jdesmis@gmail.com>
    *   2015 Ingo Weyrich <heckflosse@i-weyrich.de>

    * D. J. Jobson, Z. Rahman, and G. A. Woodell. A multi-scale
    * Retinex for bridging the gap between color images and the
    * human observation of scenes. IEEE Transactions on Image Processing,
    * 1997, 6(7): 965-976

    * Fan Guo Zixing Cai Bin Xie Jin Tang
    * School of Information Science and Engineering, Central South University Changsha, China

    * Weixing Wang and Lian Xu
    * College of Physics and Information Engineering, Fuzhou University, Fuzhou, China

    * inspired from 2003 Fabien Pelisson <Fabien.Pelisson@inrialpes.fr>

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "rtengine.h"
#include "gauss.h"
#include "rawimagesource.h"
#include "improcfun.h"
#include "opthelper.h"
#include "StopWatch.h"
#define MAX_RETINEX_SCALES   8
#define clipretinex( val, minv, maxv )    (( val = (val < minv ? minv : val ) ) > maxv ? maxv : val )

#define med3(a0,a1,a2,a3,a4,a5,a6,a7,a8,median) { \
pp[0]=a0; pp[1]=a1; pp[2]=a2; pp[3]=a3; pp[4]=a4; pp[5]=a5; pp[6]=a6; pp[7]=a7; pp[8]=a8; \
PIX_SORT(pp[1],pp[2]); PIX_SORT(pp[4],pp[5]); PIX_SORT(pp[7],pp[8]); \
PIX_SORT(pp[0],pp[1]); PIX_SORT(pp[3],pp[4]); PIX_SORT(pp[6],pp[7]); \
PIX_SORT(pp[1],pp[2]); PIX_SORT(pp[4],pp[5]); PIX_SORT(pp[7],pp[8]); \
PIX_SORT(pp[0],pp[3]); PIX_SORT(pp[5],pp[8]); PIX_SORT(pp[4],pp[7]); \
PIX_SORT(pp[3],pp[6]); PIX_SORT(pp[1],pp[4]); PIX_SORT(pp[2],pp[5]); \
PIX_SORT(pp[4],pp[7]); PIX_SORT(pp[4],pp[2]); PIX_SORT(pp[6],pp[4]); \
PIX_SORT(pp[4],pp[2]); median=pp[4];} //pp4 = median

namespace rtengine
{

extern const Settings* settings;

static float RetinexScales[MAX_RETINEX_SCALES];

void retinex_scales( float* scales, int nscales, int mode, int s)
{
    if ( nscales == 1 ) {
        scales[0] =  (float)s / 2.f;
    } else if (nscales == 2) {
        scales[0] = (float) s / 2.f;
        scales[1] = (float) s;
    } else {
        float size_step = (float) s / (float) nscales;

        if (mode == 0) {
            for (int i = 0; i < nscales; ++i ) {
                scales[i] = 2.0f + (float)i * size_step;
            }
        } else if (mode == 1) {
            size_step = (float)log(s - 2.0f) / (float) nscales;

            for (int i = 0; i < nscales; ++i ) {
                scales[i] = 2.0f + (float)pow (10.f, (i * size_step) / log (10.f));
            }
        } else if (mode == 2) {
            size_step = (float) log(s - 2.0f) / (float) nscales;

            for ( int i = 0; i < nscales; ++i ) {
                scales[i] = s - (float)pow (10.f, (i * size_step) / log (10.f));
            }
        }
    }
}
void mean_stddv2( float **dst, float &mean, float &stddv, int W_L, int H_L, float &maxtr, float &mintr)
{
    // summation using double precision to avoid too large summation error for large pictures
    double vsquared = 0.f;
    double sum = 0.f;
    maxtr = -999999.f;
    mintr = 999999.f;
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        float lmax = -999999.f, lmin = 999999.f;
#ifdef _OPENMP
        #pragma omp for reduction(+:sum,vsquared) nowait // this leads to differences, but parallel summation is more accurate
#endif

        for (int i = 0; i < H_L; i++ )
            for (int j = 0; j < W_L; j++) {
                sum += dst[i][j];
                vsquared += (dst[i][j] * dst[i][j]);

                if ( dst[i][j] > lmax) {
                    lmax = dst[i][j] ;
                }

                if ( dst[i][j] < lmin) {
                    lmin = dst[i][j] ;
                }

            }

#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            maxtr = maxtr > lmax ? maxtr : lmax;
            mintr = mintr < lmin ? mintr : lmin;
        }

    }
    mean = sum / (double) (W_L * H_L);
    vsquared /= (double) W_L * H_L;
    stddv = ( vsquared - (mean * mean) );
    stddv = (float)sqrt(stddv);
}






void mean_stddv( float **dst, float &mean, float &stddv, int W_L, int H_L, const float factor, float &maxtr, float &mintr)

{
    // summation using double precision to avoid too large summation error for large pictures
    double vsquared = 0.f;
    double sum = 0.f;
    maxtr = 0.f;
    mintr = 0.f;
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        float lmax = 0.f, lmin = 0.f;

#ifdef _OPENMP
        #pragma omp for reduction(+:sum,vsquared) // this can lead to differences, but parallel summation is more accurate
#endif

        for (int i = 0; i < H_L; i++ )
            for (int j = 0; j < W_L; j++) {
                sum += dst[i][j];
                vsquared += (dst[i][j] * dst[i][j]);

                if ( dst[i][j] > lmax) {
                    lmax = dst[i][j] ;
                }

                if ( dst[i][j] < lmin) {
                    lmin = dst[i][j] ;
                }

            }

#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            maxtr = maxtr > lmax ? maxtr : lmax;
            mintr = mintr < lmin ? mintr : lmin;
        }

    }

    sum *= factor;
    maxtr *= factor;
    mintr *= factor;
    vsquared *= (factor * factor);
    mean = sum / (float) (W_L * H_L);
    vsquared /= (float) W_L * H_L;
    stddv = ( vsquared - (mean * mean) );
    stddv = (float)sqrt(stddv);
}

void RawImageSource::MSR(float** luminance, float** originalLuminance, int width, int height, RetinexParams deh, const RetinextransmissionCurve & dehatransmissionCurve, float &minCD, float &maxCD, float &mini, float &maxi, float &Tmean, float &Tsigma, float &Tmin, float &Tmax)
{
    if (deh.enabled) {//enabled
        StopWatch Stop1("MSR");
        float         mean, stddv, maxtr, mintr;
        //  float         mini, delta, maxi;
        float         delta;
        float eps = 2.f;
        bool useHsl = deh.retinexcolorspace == "HSLLOG";
        bool useHslLin = deh.retinexcolorspace == "HSLLIN";
        float gain2 = (float) deh.gain / 100.f; //def =1  not use
        gain2 = useHslLin ? gain2 * 0.5f : gain2;
        float offse = (float) deh.offs; //def = 0  not use
        int scal =  deh.scal; //def=3
        int nei = (int) 2.8f * deh.neigh; //def = 220
        float vart = (float)deh.vart / 100.f;//variance
        float strength = (float) deh.str / 100.f; // Blend with original L channel data
        float limD = (float) deh.limd;
        limD = pow(limD, 1.7f);//about 2500 enough
        limD *= useHslLin ? 10.f : 1.f;
        float ilimD = 1.f / limD;
        int moderetinex = 2; // default to 2 ( deh.retinexMethod == "high" )
        bool execcur = false;

        if (deh.retinexMethod == "uni") {
            moderetinex = 0;
        }

        if (deh.retinexMethod == "low") {
            moderetinex = 1;
        }

        retinex_scales( RetinexScales, scal, moderetinex, nei );

        int H_L = height;
        int W_L = width;

        float *src[H_L] ALIGNED16;
        float *srcBuffer = new float[H_L * W_L];

        for (int i = 0; i < H_L; i++) {
            src[i] = &srcBuffer[i * W_L];
        }

#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = 0; i < H_L; i++)
            for (int j = 0; j < W_L; j++) {
                src[i][j] = luminance[i][j] + eps;
                luminance[i][j] = 0.f;
            }

        float *out[H_L] ALIGNED16;
        float *outBuffer = new float[H_L * W_L];

        for (int i = 0; i < H_L; i++) {
            out[i] = &outBuffer[i * W_L];
        }

        float logBetaGain = xlogf(16384.f);
        float pond = logBetaGain / (float) scal;

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            AlignedBufferMP<double>* pBuffer = new AlignedBufferMP<double> (max(W_L, H_L));

            for ( int scale = 0; scale < scal; scale++ ) {
                gaussHorizontal<float> (src, out, *pBuffer, W_L, H_L, RetinexScales[scale]);
                gaussVertical<float>   (out, out, *pBuffer, W_L, H_L, RetinexScales[scale]);

#ifdef __SSE2__
                vfloat pondv = F2V(pond);
                vfloat limMinv = F2V(ilimD);
                vfloat limMaxv = F2V(limD);

#endif
#ifdef _OPENMP
                #pragma omp for
#endif

                for (int i = 0; i < H_L; i++) {
                    int j = 0;
#ifdef __SSE2__

                    if(useHslLin) {
                        for (; j < W_L - 3; j += 4) {
                            _mm_storeu_ps(&luminance[i][j], LVFU(luminance[i][j]) + pondv *  (LIMV(LVFU(src[i][j]) / LVFU(out[i][j]), limMinv, limMaxv) ));
                        }
                    } else {
                        for (; j < W_L - 3; j += 4) {
                            _mm_storeu_ps(&luminance[i][j], LVFU(luminance[i][j]) + pondv *  xlogf(LIMV(LVFU(src[i][j]) / LVFU(out[i][j]), limMinv, limMaxv) ));
                        }
                    }

#endif

                    if(useHslLin) {
                        for (; j < W_L; j++) {
                            luminance[i][j] +=  pond * (LIM(src[i][j] / out[i][j], ilimD, limD));
                        }
                    } else {
                        for (; j < W_L; j++) {
                            luminance[i][j] +=  pond * xlogf(LIM(src[i][j] / out[i][j], ilimD, limD));
                        }
                    }
                }
            }

            delete pBuffer;
        }

        delete [] outBuffer;
        delete [] srcBuffer;


        if (dehatransmissionCurve) {
            execcur = true;
        }

        mean = 0.f;
        stddv = 0.f;
        // I call mean_stddv2 instead of mean_stddv ==> logBetaGain

        mean_stddv2( luminance, mean, stddv, W_L, H_L, maxtr, mintr);
//        printf("mean=%f std=%f delta=%f maxtr=%f mintr=%f\n", mean, stddv, delta, maxtr, mintr);

        //  mean_stddv( luminance, mean, stddv, W_L, H_L, logBetaGain, maxtr, mintr);
        if (execcur) { //if curve
            float asig = 0.166666f / stddv;
            float bsig = 0.5f - asig * mean;
            //float insigma = 0.66666f; //SD
            float amean = 0.5f / mean;
            float asign = 0.166666f / stddv;
            float bsign = 0.5f - asign * mean;
            float amax = 0.333333f / (maxtr - mean - stddv);
            float bmax = 1.f - amax * maxtr;
            float amin = 0.333333f / (mean - stddv - mintr);
            float bmin = -amin * mintr;

#ifdef _OPENMP
            #pragma omp parallel
#endif
            {
                float absciss;
#ifdef _OPENMP
                #pragma omp for schedule(dynamic,16)
#endif

                for (int i = 0; i < H_L; i++ )
                    for (int j = 0; j < W_L; j++) { //for mintr to maxtr evalate absciss in function of original transmission
                        if (luminance[i][j] >= mean && luminance[i][j] < mean + stddv) {
                            absciss = asig * luminance[i][j]  + bsig;
                        } else if (luminance[i][j] >= mean + stddv) {
                            absciss = amax * luminance[i][j]  + bmax;
                        } else if (/*luminance[i][j] < mean && */luminance[i][j] > mean - stddv) {
                            absciss = asign * luminance[i][j]  + bsign;
                        } else { /*if(luminance[i][j] <= mean - stddv)*/
                            absciss = amin * luminance[i][j]  + bmin;
                        }

                        float kmul = 2.5f;
                        float kinterm = 1.f + kmul * (dehatransmissionCurve[absciss * 500.f] - 0.5f); //new transmission
                        luminance[i][j] *= kinterm;
//                        luminance[i][j] *= 1.000001f;
                    }
            }

            // median filter on transmission  ==> reduce artifacts
            if (deh.medianmap) {
                int wid = W_L;
                int hei = H_L;
                float *tmL[hei] ALIGNED16;
                float *tmLBuffer = new float[wid * hei];
                int borderL = 1;

                for (int i = 0; i < hei; i++) {
                    tmL[i] = &tmLBuffer[i * wid];
                }

                /*
                                    for(int i = borderL; i < hei - borderL; i++ ) {
                                        for(int j = borderL; j < wid - borderL; j++) {
                                            tmL[i][j] = luminance[i][j];
                                        }
                                    }
                */
#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int i = borderL; i < hei - borderL; i++) {
                    float pp[9], temp;

                    for (int j = borderL; j < wid - borderL; j++) {
                        med3(luminance[i][j], luminance[i - 1][j], luminance[i + 1][j], luminance[i][j + 1], luminance[i][j - 1], luminance[i - 1][j - 1], luminance[i - 1][j + 1], luminance[i + 1][j - 1], luminance[i + 1][j + 1], tmL[i][j]); //3x3
                    }
                }

#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int i = borderL; i < hei - borderL; i++ ) {
                    for (int j = borderL; j < wid - borderL; j++) {
                        luminance[i][j] = tmL[i][j];
                    }
                }

                delete [] tmLBuffer;

            }

            // I call mean_stddv2 instead of mean_stddv ==> logBetaGain
            //  mean_stddv( luminance, mean, stddv, W_L, H_L, 1.f, maxtr, mintr);
            mean_stddv2( luminance, mean, stddv, W_L, H_L, maxtr, mintr);

        }

        float epsil = 0.1f;

        mini = mean - vart * stddv;

        if (mini < mintr) {
            mini = mintr + epsil;
        }

        maxi = mean + vart * stddv;

        if (maxi > maxtr) {
            maxi = maxtr - epsil;
        }

        delta = maxi - mini;
        printf("maxi=%f mini=%f mean=%f std=%f delta=%f maxtr=%f mintr=%f\n", maxi, mini, mean, stddv, delta, maxtr, mintr);

        if ( !delta ) {
            delta = 1.0f;
        }

        float cdfactor = gain2 * 32768.f / delta;
        maxCD = -9999999.f;
        minCD = 9999999.f;

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            float cdmax = -999999.f, cdmin = 999999.f;

#ifdef _OPENMP
            #pragma omp for
#endif

            for ( int i = 0; i < H_L; i ++ )
                for (int j = 0; j < W_L; j++) {
                    //   float cd = cdfactor * ( luminance[i][j] * logBetaGain - mini ) + offse;
                    float cd = cdfactor * ( luminance[i][j]  - mini ) + offse;

                    if(cd > cdmax) {
                        cdmax = cd;
                    }

                    if(cd < cdmin) {
                        cdmin = cd;
                    }


                    luminance[i][j] = clipretinex( cd, 0.f, 32768.f ) * strength + (1.f - strength) * originalLuminance[i][j];
                }

#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                maxCD = maxCD > cdmax ? maxCD : cdmax;
                minCD = minCD < cdmin ? minCD : cdmin;
            }

        }
        //    printf("cdmin=%f cdmax=%f\n",minCD, maxCD);
        Tmean = mean;
        Tsigma = stddv;
        Tmin = mintr;
        Tmax = maxtr;
    }
}

}
