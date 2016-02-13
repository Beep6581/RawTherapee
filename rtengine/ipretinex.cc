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
    *   2015 Ingo Weyrich <heckflosse67@gmx.de>

    * D. J. Jobson, Z. Rahman, and G. A. Woodell. A multi-scale
    * Retinex for bridging the gap between color images and the
    * human observation of scenes. IEEE Transactions on Image Processing,
    * 1997, 6(7): 965-976

    * Fan Guo Zixing Cai Bin Xie Jin Tang
    * School of Information Science and Engineering, Central South University Changsha, China

    * Weixing Wang and Lian Xu
    * College of Physics and Information Engineering, Fuzhou University, Fuzhou, China

    * inspired from 2003 Fabien Pelisson <Fabien.Pelisson@inrialpes.fr>
    * some ideas taken (use of mask) Russell Cottrell - The Retinex .8bf Plugin

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

void retinex_scales( float* scales, int nscales, int mode, int s, float high)
{
    if ( nscales == 1 ) {
        scales[0] =  (float)s / 2.f;
    } else if (nscales == 2) {
        scales[1] = (float) s / 2.f;
        scales[0] = (float) s;
    } else {
        float size_step = (float) s / (float) nscales;

        if (mode == 0) {
            for (int i = 0; i < nscales; ++i ) {
                scales[nscales - i - 1] = 2.0f + (float)i * size_step;
            }
        } else if (mode == 1) {
            size_step = (float)log(s - 2.0f) / (float) nscales;

            for (int i = 0; i < nscales; ++i ) {
                scales[nscales - i - 1] = 2.0f + (float)pow (10.f, (i * size_step) / log (10.f));
            }
        } else if (mode == 2) {
            size_step = (float) log(s - 2.0f) / (float) nscales;

            for ( int i = 0; i < nscales; ++i ) {
                scales[i] = s - (float)pow (10.f, (i * size_step) / log (10.f));
            }
        } else if (mode == 3) {
            size_step = (float) log(s - 2.0f) / (float) nscales;

            for ( int i = 0; i < nscales; ++i ) {
                scales[i] = high * s - (float)pow (10.f, (i * size_step) / log (10.f));
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

void RawImageSource::MSR(float** luminance, float** originalLuminance, float **exLuminance,  LUTf & mapcurve, bool &mapcontlutili, int width, int height, RetinexParams deh, const RetinextransmissionCurve & dehatransmissionCurve, const RetinexgaintransmissionCurve & dehagaintransmissionCurve, float &minCD, float &maxCD, float &mini, float &maxi, float &Tmean, float &Tsigma, float &Tmin, float &Tmax)
{

    if (deh.enabled) {//enabled
        float         mean, stddv, maxtr, mintr;
        //float         mini, delta, maxi;
        float         delta;
        float eps = 2.f;
        bool useHsl = deh.retinexcolorspace == "HSLLOG";
        bool useHslLin = deh.retinexcolorspace == "HSLLIN";
        float gain2 = (float) deh.gain / 100.f; //def =1  not use
        gain2 = useHslLin ? gain2 * 0.5f : gain2;
        float offse = (float) deh.offs; //def = 0  not use
        int iter = deh.iter;
        int gradient  =  deh.scal;
        int scal = 3;//disabled scal
        int nei = (int) 2.8f * deh.neigh; //def = 220
        float vart = (float)deh.vart / 100.f;//variance
        float gradvart = (float)deh.grad;
        float gradstr = (float)deh.grads;
        float strength = (float) deh.str / 100.f; // Blend with original L channel data
        float limD = (float) deh.limd;
        limD = pow(limD, 1.7f);//about 2500 enough
        limD *= useHslLin ? 10.f : 1.f;
        float ilimD = 1.f / limD;
        int moderetinex = 2; // default to 2 ( deh.retinexMethod == "high" )
        float hig = ((float) deh.highl) / 100.f;
        bool higplus = false ;
        float elogt;
        float hl = deh.baselog;
        scal = deh.skal;

        if(hl >= 2.71828f) {
            elogt = 2.71828f + SQR(SQR(hl - 2.71828f));
        } else {
            elogt = hl;
        }

        int H_L = height;
        int W_L = width;

        float *tran[H_L] ALIGNED16;
        float *tranBuffer;
        int viewmet = 0;

        elogt = 2.71828f;//disabled baselog
        FlatCurve* shcurve = NULL;//curve L=f(H)
        bool lhutili = false;

        if (deh.enabled) {
            shcurve = new FlatCurve(deh.lhcurve);

            if (!shcurve || shcurve->isIdentity()) {
                if (shcurve) {
                    delete shcurve;
                    shcurve = NULL;
                }
            } else {
                lhutili = true;
            }
        }


        if(deh.retinexMethod == "highliplus") {
            higplus = true;
        }

        if (deh.retinexMethod == "uni") {
            moderetinex = 0;
        }

        if (deh.retinexMethod == "low") {
            moderetinex = 1;
        }

        if (deh.retinexMethod == "highli" || deh.retinexMethod == "highliplus") {
            moderetinex = 3;
        }

        for(int it = 1; it < iter + 1; it++) { //iter nb max of iterations
            float aahi = 49.f / 99.f; ////reduce sensibility 50%
            float bbhi = 1.f - aahi;
            float high;
            high =  bbhi + aahi * (float) deh.highl;

            float grads;
            float grad = 1.f;
            float sc = scal;

            if(gradient == 0) {
                grad = 1.f;
                sc = 3.f;
            } else if(gradient == 1) {
                grad = 0.25f * it + 0.75f;
                sc = -0.5f * it + 4.5f;
            } else if(gradient == 2) {
                grad = 0.5f * it + 0.5f;
                sc = -0.75f * it + 5.75f;
            } else if(gradient == 3) {
                grad = 0.666f * it + 0.333f;
                sc = -0.75f * it + 5.75f;
            } else if(gradient == 4) {
                grad = 0.8f * it + 0.2f;
                sc = -0.75f * it + 5.75f;
            } else if(gradient == 5) {
                if(moderetinex != 3) {
                    grad = 2.5f * it - 1.5f;
                } else {
                    float aa = (11.f * high - 1.f) / 4.f;
                    float bb = 1.f - aa;
                    grad = aa * it + bb;
                }

                sc = -0.75f * it + 5.75f;
            } else if(gradient == 6) {
                if(moderetinex != 3) {
                    grad = 5.f * it - 4.f;
                } else {
                    float aa = (21.f * high - 1.f) / 4.f;
                    float bb = 1.f - aa;
                    grad = aa * it + bb;
                }

                sc = -0.75f * it + 5.75f;
            }

            else if(gradient == -1) {
                grad = -0.125f * it + 1.125f;
                sc = 3.f;
            }

            if(iter == 1) {
                sc = scal;
            } else {
                //adjust sc in function of choice of scale by user if iterations
                if(scal < 3) {
                    sc -= 1;

                    if(sc < 1.f) {//avoid 0
                        sc = 1.f;
                    }
                }

                if(scal > 4) {
                    sc += 1;
                }
            }

            float varx;
            float limdx, ilimdx;

            if(gradvart != 0) {
                if(gradvart == 1) {
                    varx = vart * (-0.125f * it + 1.125f);
                    limdx = limD * (-0.125f * it + 1.125f);
                    ilimdx = 1.f / limdx;
                } else if(gradvart == 2) {
                    varx = vart * (-0.2f * it + 1.2f);
                    limdx = limD * (-0.2f * it + 1.2f);
                    ilimdx = 1.f / limdx;
                } else if(gradvart == -1) {
                    varx = vart * (0.125f * it + 0.875f);
                    limdx = limD * (0.125f * it + 0.875f);
                    ilimdx = 1.f / limdx;
                } else if(gradvart == -2) {
                    varx = vart * (0.4f * it + 0.6f);
                    limdx = limD * (0.4f * it + 0.6f);
                    ilimdx = 1.f / limdx;
                }
            } else {
                varx = vart;
                limdx = limD;
                ilimdx = ilimD;
            }

            scal = round(sc);
            float strengthx;
            float ks = 1.f;

            if(gradstr != 0) {
                if(gradstr == 1) {
                    if(it <= 3) {
                        ks = -0.3f * it + 1.6f;
                    } else {
                        ks = 0.5f;
                    }
                } else if(gradstr == 2) {
                    if(it <= 3) {
                        ks = -0.6f * it + 2.2f;
                    } else {
                        ks = 0.3f;
                    }
                } else if(gradstr == -1) {
                    if(it <= 3) {
                        ks = 0.2f * it + 0.6f;
                    } else {
                        ks = 1.2f;
                    }
                } else if(gradstr == -2) {
                    if(it <= 3) {
                        ks = 0.4f * it + 0.2f;
                    } else {
                        ks = 1.5f;
                    }
                }
            }

            strengthx = ks * strength;
            //printf("scale=%d\n", scal);
            retinex_scales( RetinexScales, scal, moderetinex, nei / grad, high );

            float *src[H_L] ALIGNED16;
            float *srcBuffer = new float[H_L * W_L];

            for (int i = 0; i < H_L; i++) {
                src[i] = &srcBuffer[i * W_L];
            }

            int h_th, s_th;

            int shHighlights = deh.highlights;
            int shShadows = deh.shadows;
            int mapmet = 0;

            if(deh.mapMethod == "map") {
                mapmet = 2;
            }

            if(deh.mapMethod == "mapT") {
                mapmet = 3;
            }

            /*if(deh.mapMethod == "curv") {
                mapmet = 1;
            }*/

            if(deh.mapMethod == "gaus") {
                mapmet = 4;
            }

            const double shradius = mapmet == 4 ? (double) deh.radius : 40.;

            if(deh.viewMethod == "mask") {
                viewmet = 1;
            }

            if(deh.viewMethod == "tran") {
                viewmet = 2;
            }

            if(deh.viewMethod == "tran2") {
                viewmet = 3;
            }

            if(deh.viewMethod == "unsharp") {
                viewmet = 4;
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

            if(viewmet == 3  || viewmet == 2) {
                tranBuffer = new float[H_L * W_L];

                for (int i = 0; i < H_L; i++) {
                    tran[i] = &tranBuffer[i * W_L];
                }
            }

            const float logBetaGain = xlogf(16384.f);
            float pond = logBetaGain / (float) scal;

            if(!useHslLin) {
                pond /= log(elogt);
            }

            auto shmap = ((mapmet == 2 || mapmet == 3 || mapmet == 4) && it == 1) ? new SHMap (W_L, H_L, true) : nullptr;

            float *buffer = new float[W_L * H_L];;

            for ( int scale = scal - 1; scale >= 0; scale-- ) {
#ifdef _OPENMP
                #pragma omp parallel
#endif
                {
                    if(scale == scal - 1)
                    {
                        gaussianBlur (src, out, W_L, H_L, RetinexScales[scale], buffer);
                    } else { // reuse result of last iteration
                        // out was modified in last iteration => restore it
                        if((((mapmet == 2 && scale > 1) || mapmet == 3 || mapmet == 4) || (mapmet > 0 && mapcontlutili)) && it == 1)
                        {
#ifdef _OPENMP
                            #pragma omp for
#endif

                            for (int i = 0; i < H_L; i++) {
                                for (int j = 0; j < W_L; j++) {
                                    out[i][j] = buffer[i * W_L + j];
                                }
                            }
                        }

                        gaussianBlur (out, out, W_L, H_L, sqrtf(SQR(RetinexScales[scale]) - SQR(RetinexScales[scale + 1])), buffer);
                    }
                    if((((mapmet == 2 && scale > 2) || mapmet == 3 || mapmet == 4) || (mapmet > 0 && mapcontlutili)) && it == 1 && scale > 0)
                    {
                        // out will be modified => store it for use in next iteration. We even don't need a new buffer because 'buffer' is free after gaussianBlur :)
#ifdef _OPENMP
                        #pragma omp for
#endif

                        for (int i = 0; i < H_L; i++) {
                            for (int j = 0; j < W_L; j++) {
                                buffer[i * W_L + j] = out[i][j];
                            }
                        }
                    }
                }

                if(((mapmet == 2 && scale > 2) || mapmet == 3 || mapmet == 4) && it == 1) {
                    shmap->updateL (out, shradius, true, 1);
                    h_th = shmap->max_f - deh.htonalwidth * (shmap->max_f - shmap->avg) / 100;
                    s_th = deh.stonalwidth * (shmap->avg - shmap->min_f) / 100;
                }

#ifdef __SSE2__
                vfloat pondv = F2V(pond);
                vfloat limMinv = F2V(ilimdx);
                vfloat limMaxv = F2V(limdx);

#endif

                if(mapmet > 0 && mapcontlutili && it == 1) {
#ifdef _OPENMP
                    #pragma omp parallel for
#endif

                    for (int i = 0; i < H_L; i++) {
                        for (int j = 0; j < W_L; j++) {
                            out[i][j] = mapcurve[2.f * out[i][j]] / 2.f;
                        }
                    }

                }

                if(((mapmet == 2 && scale > 2) || mapmet == 3 || mapmet == 4) && it == 1) {


#ifdef _OPENMP
                    #pragma omp parallel for
#endif

                    for (int i = 0; i < H_L; i++) {
                        for (int j = 0; j < W_L; j++) {
                            double mapval = 1.0 + shmap->map[i][j];
                            double factor = 1.0;

                            if (mapval > h_th) {
                                factor = (h_th + (100.0 - shHighlights) * (mapval - h_th) / 100.0) / mapval;
                            } else if (mapval < s_th) {
                                factor = (s_th - (100.0 - shShadows) * (s_th - mapval) / 100.0) / mapval;
                            }

                            out[i][j] *= factor;

                        }
                    }
                }

#ifdef _OPENMP
                #pragma omp parallel for
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
                            luminance[i][j] +=  pond * (LIM(src[i][j] / out[i][j], ilimdx, limdx));
                        }
                    } else {
                        for (; j < W_L; j++) {
                            luminance[i][j] +=  pond * xlogf(LIM(src[i][j] / out[i][j], ilimdx, limdx)); //  /logt ?
                        }
                    }
                }
            }

            if(mapmet > 1) {
                if(shmap) {
                    delete shmap;
                }
            }

            shmap = NULL;

            delete [] buffer;
            //delete [] outBuffer;
            delete [] srcBuffer;

            mean = 0.f;
            stddv = 0.f;
            // I call mean_stddv2 instead of mean_stddv ==> logBetaGain

            mean_stddv2( luminance, mean, stddv, W_L, H_L, maxtr, mintr);
            //printf("mean=%f std=%f delta=%f maxtr=%f mintr=%f\n", mean, stddv, delta, maxtr, mintr);

            //mean_stddv( luminance, mean, stddv, W_L, H_L, logBetaGain, maxtr, mintr);
            if (dehatransmissionCurve && mean != 0.f && stddv != 0.f) { //if curve
                float asig = 0.166666f / stddv;
                float bsig = 0.5f - asig * mean;
                float amax = 0.333333f / (maxtr - mean - stddv);
                float bmax = 1.f - amax * maxtr;
                float amin = 0.333333f / (mean - stddv - mintr);
                float bmin = -amin * mintr;

                asig *= 500.f;
                bsig *= 500.f;
                amax *= 500.f;
                bmax *= 500.f;
                amin *= 500.f;
                bmin *= 500.f;
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
                            if (LIKELY(fabsf(luminance[i][j] - mean) < stddv)) {
                                absciss = asig * luminance[i][j] + bsig;
                            } else if (luminance[i][j] >= mean) {
                                absciss = amax * luminance[i][j] + bmax;
                            } else { /*if(luminance[i][j] <= mean - stddv)*/
                                absciss = amin * luminance[i][j] + bmin;
                            }

                            luminance[i][j] *= (-1.f + 4.f * dehatransmissionCurve[absciss]); //new transmission

                            if(viewmet == 3 || viewmet == 2) {
                                tran[i][j] = luminance[i][j];
                            }
                        }
                }

                // median filter on transmission  ==> reduce artifacts
                if (deh.medianmap && it == 1) { //only one time
                    int wid = W_L;
                    int hei = H_L;
                    float *tmL[hei] ALIGNED16;
                    float *tmLBuffer = new float[wid * hei];
                    int borderL = 1;

                    for (int i = 0; i < hei; i++) {
                        tmL[i] = &tmLBuffer[i * wid];
                    }

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
                //mean_stddv( luminance, mean, stddv, W_L, H_L, 1.f, maxtr, mintr);
                mean_stddv2( luminance, mean, stddv, W_L, H_L, maxtr, mintr);

            }

            float epsil = 0.1f;

            mini = mean - varx * stddv;

            if (mini < mintr) {
                mini = mintr + epsil;
            }

            maxi = mean + varx * stddv;

            if (maxi > maxtr) {
                maxi = maxtr - epsil;
            }

            delta = maxi - mini;
            //printf("maxi=%f mini=%f mean=%f std=%f delta=%f maxtr=%f mintr=%f\n", maxi, mini, mean, stddv, delta, maxtr, mintr);

            if ( !delta ) {
                delta = 1.0f;
            }

            //  float cdfactor = gain2 * 32768.f / delta;
            float cdfactor = 32768.f / delta;
            maxCD = -9999999.f;
            minCD = 9999999.f;
            // coeff for auto    "transmission" with 2 sigma #95% datas
            float aza = 16300.f / (2.f * stddv);
            float azb = -aza * (mean - 2.f * stddv);
            float bza = 16300.f / (2.f * stddv);
            float bzb = 16300.f - bza * (mean);

//prepare work for curve gain
#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int i = 0; i < H_L; i++) {
                for (int j = 0; j < W_L; j++) {
                    luminance[i][j] = luminance[i][j] - mini;
                }
            }

            mean = 0.f;
            stddv = 0.f;
            // I call mean_stddv2 instead of mean_stddv ==> logBetaGain

            mean_stddv2( luminance, mean, stddv, W_L, H_L, maxtr, mintr);
            float asig, bsig, amax, bmax, amin, bmin;

            if (dehagaintransmissionCurve && mean != 0.f && stddv != 0.f) { //if curve
                asig = 0.166666f / stddv;
                bsig = 0.5f - asig * mean;
                amax = 0.333333f / (maxtr - mean - stddv);
                bmax = 1.f - amax * maxtr;
                amin = 0.333333f / (mean - stddv - mintr);
                bmin = -amin * mintr;

                asig *= 500.f;
                bsig *= 500.f;
                amax *= 500.f;
                bmax *= 500.f;
                amin *= 500.f;
                bmin *= 500.f;
            }

#ifdef _OPENMP
            #pragma omp parallel
#endif
            {
                float absciss;
                float cdmax = -999999.f, cdmin = 999999.f;
#ifdef _OPENMP
                #pragma omp for
#endif

                for ( int i = 0; i < H_L; i ++ )
                    for (int j = 0; j < W_L; j++) {
                        float gan;

                        if (dehagaintransmissionCurve && mean != 0.f && stddv != 0.f) {

                            if (LIKELY(fabsf(luminance[i][j] - mean) < stddv)) {
                                absciss = asig * luminance[i][j] + bsig;
                            } else if (luminance[i][j] >= mean) {
                                absciss = amax * luminance[i][j] + bmax;
                            } else { /*if(luminance[i][j] <= mean - stddv)*/
                                absciss = amin * luminance[i][j] + bmin;
                            }


                            //    float cd = cdfactor * ( luminance[i][j]  - mini ) + offse;
                            gan = 2.f * (dehagaintransmissionCurve[absciss]); //new gain function transmission
                        } else {
                            gan = 0.5f;
                        }

                        float cd = gan * cdfactor * ( luminance[i][j] ) + offse;

                        if(cd > cdmax) {
                            cdmax = cd;
                        }

                        if(cd < cdmin) {
                            cdmin = cd;
                        }

                        float str = strengthx;

                        if(lhutili  && it == 1) { // S=f(H)
                            {
                                float HH = exLuminance[i][j];
                                float valparam;

                                if(useHsl || useHslLin) {
                                    valparam = float((shcurve->getVal(HH) - 0.5f));
                                } else {
                                    valparam = float((shcurve->getVal(Color::huelab_to_huehsv2(HH)) - 0.5f));
                                }

                                str *= (1.f + 2.f * valparam);
                            }
                        }

                        if(exLuminance[i][j] > 65535.f * hig && higplus) {
                            str *= hig;
                        }

                        if(viewmet == 0) {
                            luminance[i][j] = clipretinex( cd, 0.f, 32768.f ) * str + (1.f - str) * originalLuminance[i][j];
                        }

                        if(viewmet == 1) {
                            luminance[i][j] = out[i][j];
                        }

                        if(viewmet == 4) {
                            luminance[i][j] = (1.f + str) * originalLuminance[i][j] -  str * out[i][j];    //unsharp
                        }

                        if(viewmet == 2) {
                            if(tran[i][j] <= mean) {
                                luminance[i][j] = azb + aza * tran[i][j];    //auto values
                            } else {
                                luminance[i][j] = bzb + bza * tran[i][j];
                            }
                        }

                        if(viewmet == 3) {
                            luminance[i][j] = 1000.f + tran[i][j] * 700.f;    //arbitrary values to help display log values which are between -20 to + 30 - usage values -4 + 5
                        }

                    }

#ifdef _OPENMP
                #pragma omp critical
#endif
                {
                    maxCD = maxCD > cdmax ? maxCD : cdmax;
                    minCD = minCD < cdmin ? minCD : cdmin;
                }

            }
            delete [] outBuffer;
            outBuffer = NULL;
            //printf("cdmin=%f cdmax=%f\n",minCD, maxCD);
            Tmean = mean;
            Tsigma = stddv;
            Tmin = mintr;
            Tmax = maxtr;


            if (shcurve && it == 1) {
                delete shcurve;
            }
        }

        if(viewmet == 3 || viewmet == 2) {
            delete [] tranBuffer;
            tranBuffer = NULL;
        }

    }
}

}
