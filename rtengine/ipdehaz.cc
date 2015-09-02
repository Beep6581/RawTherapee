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

    * D. J. Jobson, Z. Rahman, and G. A. Woodell. A multi-scale
    * Retinex for bridging the gap between color images and the
    * human observation of scenes. IEEE Transactions on Image Processing,
    * 1997, 6(7): 965-976
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
#define MAX_DEHAZE_SCALES   6
#define clipdehaz( val, minv, maxv )    (( val = (val < minv ? minv : val ) ) > maxv ? maxv : val )

namespace rtengine
{

extern const Settings* settings;

static float DehazeScales[MAX_DEHAZE_SCALES];

void  dehaze_scales( float* scales, int nscales, int mode, int s)
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

void mean_stddv( float **dst, float &mean, float &stddv, int W_L, int H_L, const float factor )
{
    // summation using double precision to avoid too large summation error for large pictures
    double vsquared = 0.f;
    double sum = 0.f;

#ifdef _OPENMP
    #pragma omp parallel for reduction(+:sum,vsquared) // this can lead to differences, but parallel summation is more accurate
#endif

    for (int i = 0; i < H_L; i++ )
        for (int j = 0; j < W_L; j++) {
            sum += dst[i][j];
            vsquared += (dst[i][j] * dst[i][j]);
        }

    sum *= factor;
    vsquared *= (factor * factor);
    mean = sum / (float) (W_L * H_L);
    vsquared /= (float) W_L * H_L;
    stddv = ( vsquared - (mean * mean) );
    stddv = (float)sqrt(stddv);
}

void RawImageSource::MSR(float** luminance, float** originalLuminance, int width, int height, DehazParams lcur)
{
    if (lcur.enabled) {//enabled
        StopWatch Stop1("MSR");
        float         mean, stddv;
        float         mini, delta, maxi;
        float eps = 2.f;
        float gain2 = (float) lcur.gain / 100.f; //def =1  not use
        float offse = (float) lcur.offs; //def = 0  not use
        int scal =  lcur.scal; //def=3
        int nei = (int) 2.5f * lcur.neigh; //def = 200
        float vart = (float)lcur.vart / 100.f;//variance
        float strength = (float) lcur.str / 100.f; // Blend with original L channel data

        int modedehaz = 0; // default to 0 ( lcur.dehazmet == "uni" )

        if (lcur.dehazmet == "low") {
            modedehaz = 1;
        }

        if (lcur.dehazmet == "high") {
            modedehaz = 2;
        }

        dehaze_scales( DehazeScales, scal, modedehaz, nei );

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

        float pond = 1.0f / (float) scal;

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
        AlignedBufferMP<double>* pBuffer = new AlignedBufferMP<double> (max(W_L, H_L));

        for ( int scale = 0; scale < scal; scale++ ) {
                gaussHorizontal<float> (src, out, *pBuffer, W_L, H_L, DehazeScales[scale]);
                gaussVertical<float>   (out, out, *pBuffer, W_L, H_L, DehazeScales[scale]);
#ifdef __SSE2__
                vfloat pondv = F2V(pond);
                vfloat limMinv = F2V(0.0001f);
                vfloat limMaxv = F2V(10000.f);
#endif
#ifdef _OPENMP
                #pragma omp for
#endif

                for (int i = 0; i < H_L; i++)
                {
                    int j = 0;
#ifdef __SSE2__

                    for (; j < W_L - 3; j += 4) {
                        _mm_storeu_ps(&luminance[i][j], LVFU(luminance[i][j]) + pondv *  xlogf(LIMV(LVFU(src[i][j]) / LVFU(out[i][j]), limMinv, limMaxv) ));
                    }

#endif

                    for (; j < W_L; j++) {
                        luminance[i][j] +=  pond * xlogf(LIM(src[i][j] / out[i][j], 0.0001f, 10000.f));
                    }
                }
            }
            delete pBuffer;
        }

        delete [] outBuffer;
        delete [] srcBuffer;

        float logBetaGain = xlogf(16384.f);

        mean = 0.f;
        stddv = 0.f;
        mean_stddv( luminance, mean, stddv, W_L, H_L, logBetaGain);

        mini = mean - vart * stddv;
        maxi = mean + vart * stddv;
        delta = maxi - mini;
        printf("maxi=%f mini=%f mean=%f std=%f delta=%f\n", maxi, mini, mean, stddv, delta);

        if ( !delta ) {
            delta = 1.0f;
        }

        float cdfactor = gain2 * 32768.f / delta;
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for ( int i = 0; i < H_L; i ++ )
            for (int j = 0; j < W_L; j++) {
                float cd = cdfactor * ( luminance[i][j] * logBetaGain - mini ) + offse;
                luminance[i][j] = clipdehaz( cd, 0.f, 32768.f ) * strength + (1.f - strength) * originalLuminance[i][j];
            }
    }
}

}
