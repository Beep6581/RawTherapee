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

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include "rtengine.h"
#include "gauss.h"
#include "rawimagesource.h"
#include "improcfun.h"
#include "opthelper.h"
#include "median.h"
#include "StopWatch.h"

#define clipretinex( val, minv, maxv )    (( val = (val < minv ? minv : val ) ) > maxv ? maxv : val )

namespace
{
void retinex_scales ( float* scales, int nscales, int mode, int s, float high)
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
            size_step = (float)log (s - 2.0f) / (float) nscales;

            for (int i = 0; i < nscales; ++i ) {
                scales[nscales - i - 1] = 2.0f + (float)pow (10.f, (i * size_step) / log (10.f));
            }
        } else if (mode == 2) {
            size_step = (float) log (s - 2.0f) / (float) nscales;

            for ( int i = 0; i < nscales; ++i ) {
                scales[i] = s - (float)pow (10.f, (i * size_step) / log (10.f));
            }
        } else if (mode == 3) {
            size_step = (float) log (s - 2.0f) / (float) nscales;

            for ( int i = 0; i < nscales; ++i ) {
                scales[i] = high * s - (float)pow (10.f, (i * size_step) / log (10.f));
            }
        }
    }
}

void mean_stddv2 ( float **dst, float &mean, float &stddv, int W_L, int H_L, float &maxtr, float &mintr)
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

                lmax = dst[i][j] > lmax ? dst[i][j] : lmax;
                lmin = dst[i][j] < lmin ? dst[i][j] : lmin;
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
    stddv = (float)sqrt (stddv);
}

}


namespace rtengine
{

extern const Settings* settings;

void RawImageSource::MSR(float** luminance, float** originalLuminance, float **exLuminance,  LUTf & mapcurve, bool &mapcontlutili, int width, int height, const RetinexParams &deh, const RetinextransmissionCurve & dehatransmissionCurve, const RetinexgaintransmissionCurve & dehagaintransmissionCurve, float &minCD, float &maxCD, float &mini, float &maxi, float &Tmean, float &Tsigma, float &Tmin, float &Tmax)
{

    if (deh.enabled) {//enabled
        float maxtr, mintr;
        constexpr float eps = 2.f;
        bool useHsl = deh.retinexcolorspace == "HSLLOG";
        bool useHslLin = deh.retinexcolorspace == "HSLLIN";
        float offse = (float) deh.offs; //def = 0  not use
        int iter = deh.iter;
        int gradient = deh.scal;
        int scal = 3;//disabled scal
        int nei = (int) (2.8f * deh.neigh); //def = 220
        float vart = (float)deh.vart / 100.f;//variance
        float gradvart = (float)deh.grad;
        float gradstr = (float)deh.grads;
        float strength = (float) deh.str / 100.f; // Blend with original L channel data
        float limD = (float) deh.limd;
        limD = pow (limD, 1.7f); //about 2500 enough
        limD *= useHslLin ? 10.f : 1.f;
        float ilimD = 1.f / limD;
        float hig = ((float) deh.highl) / 100.f;
        scal = deh.skal;

        int H_L = height;
        int W_L = width;

        float *tran[H_L] ALIGNED16;
        float *tranBuffer = nullptr;

        constexpr float elogt = 2.71828f;
        bool lhutili = false;

        FlatCurve* shcurve = new FlatCurve (deh.lhcurve); //curve L=f(H)

        if (!shcurve || shcurve->isIdentity()) {
            if (shcurve) {
                delete shcurve;
                shcurve = nullptr;
            }
        } else {
            lhutili = true;
        }


        bool higplus = false ;
        int moderetinex = 2; // default to 2 ( deh.retinexMethod == "high" )

        if (deh.retinexMethod == "highliplus") {
            higplus = true;
            moderetinex = 3;
        } else if (deh.retinexMethod == "uni") {
            moderetinex = 0;
        } else if (deh.retinexMethod == "low") {
            moderetinex = 1;
        } else { /*if (deh.retinexMethod == "highli") */
            moderetinex = 3;
        }

        constexpr float aahi = 49.f / 99.f; ////reduce sensibility 50%
        constexpr float bbhi = 1.f - aahi;

        for (int it = 1; it < iter + 1; it++) { //iter nb max of iterations
            float high = bbhi + aahi * (float) deh.highl;

            float grad = 1.f;
            float sc = scal;

            if (gradient == 0) {
                grad = 1.f;
                sc = 3.f;
            } else if (gradient == 1) {
                grad = 0.25f * it + 0.75f;
                sc = -0.5f * it + 4.5f;
            } else if (gradient == 2) {
                grad = 0.5f * it + 0.5f;
                sc = -0.75f * it + 5.75f;
            } else if (gradient == 3) {
                grad = 0.666f * it + 0.333f;
                sc = -0.75f * it + 5.75f;
            } else if (gradient == 4) {
                grad = 0.8f * it + 0.2f;
                sc = -0.75f * it + 5.75f;
            } else if (gradient == 5) {
                if (moderetinex != 3) {
                    grad = 2.5f * it - 1.5f;
                } else {
                    float aa = (11.f * high - 1.f) / 4.f;
                    float bb = 1.f - aa;
                    grad = aa * it + bb;
                }

                sc = -0.75f * it + 5.75f;
            } else if (gradient == 6) {
                if (moderetinex != 3) {
                    grad = 5.f * it - 4.f;
                } else {
                    float aa = (21.f * high - 1.f) / 4.f;
                    float bb = 1.f - aa;
                    grad = aa * it + bb;
                }

                sc = -0.75f * it + 5.75f;
            }

            else if (gradient == -1) {
                grad = -0.125f * it + 1.125f;
                sc = 3.f;
            }

            if (iter == 1) {
                sc = scal;
            } else {
                //adjust sc in function of choice of scale by user if iterations
                if (scal < 3) {
                    sc -= 1;

                    if (sc < 1.f) { //avoid 0
                        sc = 1.f;
                    }
                }

                if (scal > 4) {
                    sc += 1;
                }
            }

            float varx = vart;
            float limdx = limD;
            float ilimdx = ilimD;

            if (gradvart != 0) {
                if (gradvart == 1) {
                    varx = vart * (-0.125f * it + 1.125f);
                    limdx = limD * (-0.125f * it + 1.125f);
                    ilimdx = 1.f / limdx;
                } else if (gradvart == 2) {
                    varx = vart * (-0.2f * it + 1.2f);
                    limdx = limD * (-0.2f * it + 1.2f);
                    ilimdx = 1.f / limdx;
                } else if (gradvart == -1) {
                    varx = vart * (0.125f * it + 0.875f);
                    limdx = limD * (0.125f * it + 0.875f);
                    ilimdx = 1.f / limdx;
                } else if (gradvart == -2) {
                    varx = vart * (0.4f * it + 0.6f);
                    limdx = limD * (0.4f * it + 0.6f);
                    ilimdx = 1.f / limdx;
                }
            }

            scal = round (sc);
            float ks = 1.f;

            if (gradstr != 0) {
                if (gradstr == 1) {
                    if (it <= 3) {
                        ks = -0.3f * it + 1.6f;
                    } else {
                        ks = 0.5f;
                    }
                } else if (gradstr == 2) {
                    if (it <= 3) {
                        ks = -0.6f * it + 2.2f;
                    } else {
                        ks = 0.3f;
                    }
                } else if (gradstr == -1) {
                    if (it <= 3) {
                        ks = 0.2f * it + 0.6f;
                    } else {
                        ks = 1.2f;
                    }
                } else if (gradstr == -2) {
                    if (it <= 3) {
                        ks = 0.4f * it + 0.2f;
                    } else {
                        ks = 1.5f;
                    }
                }
            }

            float strengthx = ks * strength;

            constexpr auto maxRetinexScales = 8;
            float RetinexScales[maxRetinexScales];

            retinex_scales ( RetinexScales, scal, moderetinex, nei / grad, high );

            float *src[H_L] ALIGNED16;
            float *srcBuffer = new float[H_L * W_L];

            for (int i = 0; i < H_L; i++) {
                src[i] = &srcBuffer[i * W_L];
            }

            int h_th = 0, s_th = 0;

            int shHighlights = deh.highlights;
            int shShadows = deh.shadows;

            int mapmet = 0;

            if (deh.mapMethod == "map") {
                mapmet = 2;
            } else if (deh.mapMethod == "mapT") {
                mapmet = 3;
            } else if (deh.mapMethod == "gaus") {
                mapmet = 4;
            }

            const double shradius = mapmet == 4 ? (double) deh.radius : 40.;

            int viewmet = 0;

            if (deh.viewMethod == "mask") {
                viewmet = 1;
            } else if (deh.viewMethod == "tran") {
                viewmet = 2;
            } else if (deh.viewMethod == "tran2") {
                viewmet = 3;
            } else if (deh.viewMethod == "unsharp") {
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

            if (viewmet == 3  || viewmet == 2) {
                tranBuffer = new float[H_L * W_L];

                for (int i = 0; i < H_L; i++) {
                    tran[i] = &tranBuffer[i * W_L];
                }
            }

            const float logBetaGain = xlogf (16384.f);
            float pond = logBetaGain / (float) scal;

            if (!useHslLin) {
                pond /= log (elogt);
            }

            auto shmap = ((mapmet == 2 || mapmet == 3 || mapmet == 4) && it == 1) ? new SHMap (W_L, H_L, true) : nullptr;

            float *buffer = new float[W_L * H_L];;

            for ( int scale = scal - 1; scale >= 0; scale-- ) {
#ifdef _OPENMP
                #pragma omp parallel
#endif
                {
                    if (scale == scal - 1)
                    {
                        gaussianBlur (src, out, W_L, H_L, RetinexScales[scale], buffer);
                    } else { // reuse result of last iteration
                        // out was modified in last iteration => restore it
                        if ((((mapmet == 2 && scale > 1) || mapmet == 3 || mapmet == 4) || (mapmet > 0 && mapcontlutili)) && it == 1)
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

                        gaussianBlur (out, out, W_L, H_L, sqrtf (SQR (RetinexScales[scale]) - SQR (RetinexScales[scale + 1])), buffer);
                    }
                    if ((((mapmet == 2 && scale > 2) || mapmet == 3 || mapmet == 4) || (mapmet > 0 && mapcontlutili)) && it == 1 && scale > 0)
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

                if (((mapmet == 2 && scale > 2) || mapmet == 3 || mapmet == 4) && it == 1) {
                    shmap->updateL (out, shradius, true, 1);
                    h_th = shmap->max_f - deh.htonalwidth * (shmap->max_f - shmap->avg) / 100;
                    s_th = deh.stonalwidth * (shmap->avg - shmap->min_f) / 100;
                }

#ifdef __SSE2__
                vfloat pondv = F2V (pond);
                vfloat limMinv = F2V (ilimdx);
                vfloat limMaxv = F2V (limdx);

#endif

                if (mapmet > 0 && mapcontlutili && it == 1) {
                    // TODO: When rgbcurvespeedup branch is merged into master we can simplify the code by
                    // 1) in rawimagesource.retinexPrepareCurves() insert
                    //    mapcurve *= 0.5f;
                    //    after
                    //    CurveFactory::mapcurve (mapcontlutili, retinexParams.mapcurve, mapcurve, 1, lhist16RETI, histLRETI);
                    // 2) remove the division by 2.f from the code 7 lines below this line
#ifdef _OPENMP
                    #pragma omp parallel for
#endif

                    for (int i = 0; i < H_L; i++) {
                        for (int j = 0; j < W_L; j++) {
                            out[i][j] = mapcurve[2.f * out[i][j]] / 2.f;
                        }
                    }

                }

                if (((mapmet == 2 && scale > 2) || mapmet == 3 || mapmet == 4) && it == 1) {
                    float hWeight = (100.f - shHighlights) / 100.f;
                    float sWeight = (100.f - shShadows) / 100.f;
#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16)
#endif

                    for (int i = 0; i < H_L; i++) {
                        for (int j = 0; j < W_L; j++) {
                            float mapval = 1.f + shmap->map[i][j];
                            float factor = 1.f;

                            if (mapval > h_th) {
                                factor = (h_th + hWeight * (mapval - h_th)) / mapval;
                            } else if (mapval < s_th) {
                                factor = (s_th - sWeight * (s_th - mapval)) / mapval;
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

                    if (useHslLin) {
                        for (; j < W_L - 3; j += 4) {
                            _mm_storeu_ps (&luminance[i][j], LVFU (luminance[i][j]) + pondv *  (LIMV (LVFU (src[i][j]) / LVFU (out[i][j]), limMinv, limMaxv) ));
                        }
                    } else {
                        for (; j < W_L - 3; j += 4) {
                            _mm_storeu_ps (&luminance[i][j], LVFU (luminance[i][j]) + pondv *  xlogf (LIMV (LVFU (src[i][j]) / LVFU (out[i][j]), limMinv, limMaxv) ));
                        }
                    }

#endif

                    if (useHslLin) {
                        for (; j < W_L; j++) {
                            luminance[i][j] +=  pond * (LIM (src[i][j] / out[i][j], ilimdx, limdx));
                        }
                    } else {
                        for (; j < W_L; j++) {
                            luminance[i][j] +=  pond * xlogf (LIM (src[i][j] / out[i][j], ilimdx, limdx)); //  /logt ?
                        }
                    }
                }
            }

            if (mapmet > 1) {
                if (shmap) {
                    delete shmap;
                }
            }

            shmap = nullptr;

            delete [] buffer;
            delete [] srcBuffer;

            float mean = 0.f;
            float stddv = 0.f;
            // I call mean_stddv2 instead of mean_stddv ==> logBetaGain

            mean_stddv2 ( luminance, mean, stddv, W_L, H_L, maxtr, mintr);
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
                            if (LIKELY (fabsf (luminance[i][j] - mean) < stddv)) {
                                absciss = asig * luminance[i][j] + bsig;
                            } else if (luminance[i][j] >= mean) {
                                absciss = amax * luminance[i][j] + bmax;
                            } else { /*if(luminance[i][j] <= mean - stddv)*/
                                absciss = amin * luminance[i][j] + bmin;
                            }

                            //TODO : move multiplication by 4.f and subtraction of 1.f inside the curve
                            luminance[i][j] *= (-1.f + 4.f * dehatransmissionCurve[absciss]); //new transmission

                            if (viewmet == 3 || viewmet == 2) {
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
                        for (int j = borderL; j < wid - borderL; j++) {
                            tmL[i][j] = median (luminance[i][j], luminance[i - 1][j], luminance[i + 1][j], luminance[i][j + 1], luminance[i][j - 1], luminance[i - 1][j - 1], luminance[i - 1][j + 1], luminance[i + 1][j - 1], luminance[i + 1][j + 1]); //3x3
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
                mean_stddv2 ( luminance, mean, stddv, W_L, H_L, maxtr, mintr);

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

            float delta = maxi - mini;
            //printf("maxi=%f mini=%f mean=%f std=%f delta=%f maxtr=%f mintr=%f\n", maxi, mini, mean, stddv, delta, maxtr, mintr);

            if ( !delta ) {
                delta = 1.0f;
            }

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

            mean_stddv2 ( luminance, mean, stddv, W_L, H_L, maxtr, mintr);
            float asig = 0.f, bsig = 0.f, amax = 0.f, bmax = 0.f, amin = 0.f, bmin = 0.f;

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
                float cdmax = -999999.f, cdmin = 999999.f;
#ifdef _OPENMP
                #pragma omp for schedule(dynamic,16) nowait
#endif

                for ( int i = 0; i < H_L; i ++ )
                    for (int j = 0; j < W_L; j++) {
                        float gan;

                        if (dehagaintransmissionCurve && mean != 0.f && stddv != 0.f) {
                            float absciss;

                            if (LIKELY (fabsf (luminance[i][j] - mean) < stddv)) {
                                absciss = asig * luminance[i][j] + bsig;
                            } else if (luminance[i][j] >= mean) {
                                absciss = amax * luminance[i][j] + bmax;
                            } else { /*if(luminance[i][j] <= mean - stddv)*/
                                absciss = amin * luminance[i][j] + bmin;
                            }


                            //    float cd = cdfactor * ( luminance[i][j]  - mini ) + offse;
                            // TODO : move multiplication by 2.f inside the curve
                            gan = 2.f * (dehagaintransmissionCurve[absciss]); //new gain function transmission
                        } else {
                            gan = 0.5f;
                        }

                        float cd = gan * cdfactor * ( luminance[i][j] ) + offse;

                        cdmax = cd > cdmax ? cd : cdmax;
                        cdmin = cd < cdmin ? cd : cdmin;

                        float str = strengthx;

                        if (lhutili && it == 1) { // S=f(H)
                            {
                                float HH = exLuminance[i][j];
                                float valparam;

                                if (useHsl || useHslLin) {
                                    valparam = float ((shcurve->getVal (HH) - 0.5f));
                                } else {
                                    valparam = float ((shcurve->getVal (Color::huelab_to_huehsv2 (HH)) - 0.5f));
                                }

                                str *= (1.f + 2.f * valparam);
                            }
                        }

                        if (higplus && exLuminance[i][j] > 65535.f * hig) {
                            str *= hig;
                        }

                        if (viewmet == 0) {
                            luminance[i][j] = intp (str, clipretinex ( cd, 0.f, 32768.f ), originalLuminance[i][j]);
                        } else if (viewmet == 1) {
                            luminance[i][j] = out[i][j];
                        } else if (viewmet == 4) {
                            luminance[i][j] = originalLuminance[i][j] + str * (originalLuminance[i][j] - out[i][j]);//unsharp
                        } else if (viewmet == 2) {
                            if (tran[i][j] <= mean) {
                                luminance[i][j] = azb + aza * tran[i][j];    //auto values
                            } else {
                                luminance[i][j] = bzb + bza * tran[i][j];
                            }
                        } else { /*if(viewmet == 3) */
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
            outBuffer = nullptr;
            //printf("cdmin=%f cdmax=%f\n",minCD, maxCD);
            Tmean = mean;
            Tsigma = stddv;
            Tmin = mintr;
            Tmax = maxtr;

            if (shcurve) {
                delete shcurve;
                shcurve = nullptr;
            }
        }

        if (tranBuffer) {
            delete [] tranBuffer;
        }

    }
}

void ImProcFunctions::MSRLocal (float** luminance, float** templ, const float* const *originalLuminance, const int width, const int height, const LocallabParams &loc, const int skip, const LocretigainCurve &locRETgainCcurve, const int chrome, const int scall, const float krad, float &minCD, float &maxCD, float &mini, float &maxi, float &Tmean, float &Tsigma, float &Tmin, float &Tmax)
{
    BENCHFUN
    bool py = true;

    if (py) {//enabled
        float         mean, stddv, maxtr, mintr;
        float         delta;
        constexpr float eps = 2.f;
        //   constexpr bool useHsl = false; //never used
        constexpr bool useHslLin = false;//never used
        const float offse = 0.f; //loc.offs;
        const float chrT = (float) (loc.chrrt) / 100.f;
        const int scal = scall;//3;//loc.scale;;
        const float vart = loc.vart / 100.f;//variance
        const float strength = loc.str / 100.f; // Blend with original L channel data
        float limD = 10.f;//(float) loc.limd;
        limD = pow (limD, 1.7f); //about 2500 enough
        //limD *= useHslLin ? 10.f : 1.f;
        float ilimD = 1.f / limD;
        const float elogt = 2.71828f;

        //empirical skip evaluation : very difficult  because quasi all parameters interfere
        //to test on several images
        int nei = (int) (krad * loc.neigh);

        if (skip >= 4) {
            nei = (int) (0.1f * nei + 2.f);    //not too bad
        } else if (skip > 1 && skip < 4) {
            nei = (int) (0.3f * nei + 2.f);
        }

        int moderetinex = 0;

        if (loc.retinexMethod == "uni") {
            moderetinex = 0;
        } else if (loc.retinexMethod == "low") {
            moderetinex = 1;
        } else {
            if (loc.retinexMethod == "high") { // default to 2 ( deh.retinexMethod == "high" )
                moderetinex = 2;
            }
        }

        const float high = 0.f; // Dummy to pass to retinex_scales(...)

        constexpr auto maxRetinexScales = 8;
        float RetinexScales[maxRetinexScales];

        retinex_scales ( RetinexScales, scal, moderetinex, nei, high );


        const int H_L = height;
        const int W_L = width;
        float *src[H_L] ALIGNED16;
        float *srcBuffer = new float[H_L * W_L];

        for (int i = 0; i < H_L; i++) {
            src[i] = &srcBuffer[i * W_L];
        }

        float hi = 0.f;
        float lo = 0.f;
        const float shHighlights = (100.f - 0.f) / 100.f; //loc.highlights
        const float shShadows = (100.f - 0.f) / 100.f; //loc.shadows
        const int mapmet = (hi > 0 || lo > 0) ? 4 : 0;
        //  const double shradius = mapmet == 4 ? loc.radius : 40.;
        const double shradius = mapmet == 4 ? 40 : 40.;
        constexpr int it = 1;//in case of !!

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

        const float logBetaGain = xlogf (16384.f);
        float pond = logBetaGain / (float) scal;

        if (!useHslLin) {
            pond /= log (elogt);
        }

        auto shmap = mapmet == 4 ? new SHMap (W_L, H_L, true) : nullptr;

        float *buffer = new float[W_L * H_L];

        for ( int scale = scal - 1; scale >= 0; scale-- ) {
#ifdef _OPENMP
            #pragma omp parallel
#endif
            {

                if (scale == scal - 1)
                {
                    gaussianBlur (src, out, W_L, H_L, RetinexScales[scale], buffer);
                } else { // reuse result of last iteration
                    // out was modified in last iteration => restore it
                    if (((mapmet == 4)) && it == 1)
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

                    gaussianBlur (out, out, W_L, H_L, sqrtf (SQR (RetinexScales[scale]) - SQR (RetinexScales[scale + 1])), buffer);
                }
                if ((mapmet == 4) && it == 1 && scale > 0)
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

            float h_th, s_th;
            float h_thcomp, s_thcomp;

            if ((mapmet == 4) && it == 1) {
                shmap->updateL (out, shradius, true, 1);
                h_thcomp = 0.f;//shmap->max_f - loc.htonalwidth * (shmap->max_f - shmap->avg) / 100.f;
                h_th = h_thcomp - (shHighlights * h_thcomp);
                s_thcomp = 0.f;//loc.stonalwidth * (shmap->avg - shmap->min_f) / 100.f;
                s_th = s_thcomp - (shShadows * s_thcomp);
            }

            if ((mapmet == 4) && it == 1) {

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int i = 0; i < H_L; i++) {
                    for (int j = 0; j < W_L; j++) {
                        float mapval = 1.f + shmap->map[i][j];

                        if (mapval > h_thcomp) {
                            out[i][j] *= h_th / mapval + shHighlights;
                        } else if (mapval < s_thcomp) {
                            out[i][j] *= s_th / mapval + shShadows;
                        }
                    }
                }

            }

#ifdef __SSE2__
            vfloat pondv = F2V (pond);
            vfloat limMinv = F2V (ilimD);
            vfloat limMaxv = F2V (limD);

#endif
#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int i = 0; i < H_L; i++) {
                int j = 0;

#ifdef __SSE2__

                if (useHslLin) { //keep in case of ??
                    for (; j < W_L - 3; j += 4) {
                        _mm_storeu_ps (&luminance[i][j], LVFU (luminance[i][j]) + pondv *  (LIMV (LVFU (src[i][j]) / LVFU (out[i][j]), limMinv, limMaxv) ));
                    }
                } else {//always Lab mode due to Wavelet
                    for (; j < W_L - 3; j += 4) {
                        _mm_storeu_ps (&luminance[i][j], LVFU (luminance[i][j]) + pondv *  xlogf (LIMV (LVFU (src[i][j]) / LVFU (out[i][j]), limMinv, limMaxv) ));
                    }
                }

#endif

                if (useHslLin) {
                    for (; j < W_L; j++) {
                        luminance[i][j] +=  pond * (LIM (src[i][j] / out[i][j], ilimD, limD));
                    }
                } else {
                    for (; j < W_L; j++) {
                        luminance[i][j] +=  pond * xlogf (LIM (src[i][j] / out[i][j], ilimD, limD)); //  /logt ?
                    }
                }
            }
        }

        if (shmap) {
            delete shmap;
        }

        shmap = nullptr;

        delete [] buffer;
        delete [] outBuffer;
        outBuffer = nullptr;
        delete [] srcBuffer;

        mean = 0.f;
        stddv = 0.f;

        mean_stddv2 ( luminance, mean, stddv, W_L, H_L, maxtr, mintr);
        //printf("mean=%f std=%f delta=%f maxtr=%f mintr=%f\n", mean, stddv, delta, maxtr, mintr);

        //  mean_stddv( luminance, mean, stddv, W_L, H_L, logBetaGain, maxtr, mintr);

        bool retcurve = false;//wavRETCcurve

        if (retcurve && mean != 0.f && stddv != 0.f) { //if curve
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
            /*
            #ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16)
            #endif

                        for (int i = 0; i < H_L; i++ )
                            for (int j = 0; j < W_L; j++) { //for mintr to maxtr evalate absciss in function of original transmission
                                float absciss;

                                if (LIKELY (fabsf (luminance[i][j] - mean) < stddv)) {
                                    absciss = asig * luminance[i][j] + bsig;
                                } else if (luminance[i][j] >= mean) {
                                    absciss = amax * luminance[i][j] + bmax;
                                } else {
                                    absciss = amin * luminance[i][j] + bmin;
                                }

                                // luminance[i][j] *= (-1.f + 4.f * wavRETCcurve[absciss]); //new transmission
                            }
            */
            // median filter on transmission  ==> reduce artifacts
            bool ty = false;

            if (ty) {//not used here to simplify interface
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
                    // float pp[9], temp;

                    for (int j = borderL; j < wid - borderL; j++) {
                        tmL[i][j] = median (luminance[i][j], luminance[i - 1][j], luminance[i + 1][j], luminance[i][j + 1], luminance[i][j - 1], luminance[i - 1][j - 1], luminance[i - 1][j + 1], luminance[i + 1][j - 1], luminance[i + 1][j + 1]); //3x3
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
            mean_stddv2 ( luminance, mean, stddv, W_L, H_L, maxtr, mintr);

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
        //printf("maxi=%f mini=%f mean=%f std=%f delta=%f maxtr=%f mintr=%f\n", maxi, mini, mean, stddv, delta, maxtr, mintr);

        if ( !delta ) {
            delta = 1.0f;
        }

        float cdfactor = 32768.f / delta;
        maxCD = -9999999.f;
        minCD = 9999999.f;

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

        mean_stddv2 ( luminance, mean, stddv, W_L, H_L, maxtr, mintr);
        float asig = 0.f, bsig = 0.f, amax = 0.f, bmax = 0.f, amin = 0.f, bmin = 0.f;
        //   bool gaincurve = false; //wavRETgainCcurve
        const bool hasWavRetGainCurve =  locRETgainCcurve && mean != 0.f && stddv != 0.f;

        if (hasWavRetGainCurve) { //if curve
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
            cdfactor *= 2.f;
        }


        const float maxclip = (chrome == 0 ? 32768.f : 50000.f);
        float str = strength * (chrome == 0 ? 1.f : chrT);
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            //            float absciss;
            float cdmax = -999999.f, cdmin = 999999.f;
            float gan = 0.5f;
#ifdef _OPENMP
            #pragma omp for schedule(dynamic,16)
#endif

            for ( int i = 0; i < H_L; i ++ )
                for (int j = 0; j < W_L; j++) {
                    if (hasWavRetGainCurve) {
                        float absciss;

                        if (LIKELY (fabsf (luminance[i][j] - mean) < stddv)) {
                            absciss = asig * luminance[i][j] + bsig;
                        } else if (luminance[i][j] >= mean) {
                            absciss = amax * luminance[i][j] + bmax;
                        } else {
                            absciss = amin * luminance[i][j] + bmin;
                        }

                        gan = locRETgainCcurve[absciss]; //new gain function transmission
                    }

                    float cd = gan * cdfactor * luminance[i][j] + offse;

                    cdmax = cd > cdmax ? cd : cdmax;
                    cdmin = cd < cdmin ? cd : cdmin;
                    luminance[i][j] = LIM ( cd, 0.f, maxclip ) * str + (1.f - str) * originalLuminance[i][j];
                    //   templ[i][j] = LIM( cd, 0.f, maxclip ) * str + (1.f - str) * originalLuminance[i][j];
                    //  luminance[i][j] = LIM( cd, 0.f, maxclip );
                }

#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                maxCD = maxCD > cdmax ? maxCD : cdmax;
                minCD = minCD < cdmin ? minCD : cdmin;
            }

        }

        Tmean = mean;
        Tsigma = stddv;
        Tmin = mintr;
        Tmax = maxtr;


    }



}
}