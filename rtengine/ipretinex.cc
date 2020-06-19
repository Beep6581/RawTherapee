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
*  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.

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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "color.h"
#include "curves.h"
#include "gauss.h"
#include "improcfun.h"
#include "labimage.h"
#include "median.h"
#include "opthelper.h"
#include "procparams.h"
#include "rawimagesource.h"
#include "rtengine.h"
#include "shmap.h"
#define BENCHMARK
#include "StopWatch.h"
#include "guidedfilter.h"

#define clipretinex( val, minv, maxv )    (( val = (val < minv ? minv : val ) ) > maxv ? maxv : val )
#define CLIPLOC(x) LIM(x,0.f,32767.f)
#define CLIPC(a) LIM(a, -42000.f, 42000.f)  // limit a and b  to 130 probably enough ?

namespace
{

void calcGammaLut(double gamma, double ts, LUTf &gammaLut)
{
    double pwr = 1.0 / gamma;
    double gamm = gamma;
    const double gamm2 = gamma;
    rtengine::GammaValues g_a;

    if (gamm2 < 1.0) {
        std::swap(pwr, gamm);
    }

    rtengine::Color::calcGamma(pwr, ts, 0, g_a); // call to calcGamma with selected gamma and slope

    const double start = gamm2 < 1. ? g_a[2] : g_a[3];
    const double add = g_a[4];
    const double mul = 1.0 + g_a[4];

    if (gamm2 < 1.) {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1024)
#endif
        for (int i = 0; i < 65536; i++) {
            const double x = rtengine::Color::igammareti(i / 65535.0, gamm, start, ts, mul, add);
            gammaLut[i] = 0.5 * rtengine::CLIP(x * 65535.0);  // CLIP avoid in some case extra values
        }
    } else {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1024)
#endif
        for (int i = 0; i < 65536; i++) {
            const double x = rtengine::Color::gammareti(i / 65535.0, gamm, start, ts, mul, add);
            gammaLut[i] = 0.5 * rtengine::CLIP(x * 65535.0);  // CLIP avoid in some case extra values
        }
    }
}

void retinex_scales(float* scales, int nscales, int mode, int s, float high)
{
    if (s < 3) {
        s = 3;    //to avoid crash in MSRlocal if nei small
    }

    if (nscales == 1) {
        scales[0] = (float)s / 2.f;
    } else if (nscales == 2) {
        scales[1] = (float) s / 2.f;
        scales[0] = (float) s;
    } else {
        float size_step = (float) s / (float) nscales;

        if (mode == 0) {
            for (int i = 0; i < nscales; ++i) {
                scales[nscales - i - 1] = 2.0f + (float)i * size_step;
            }
        } else if (mode == 1) {
            size_step = (float)log(s - 2.0f) / (float) nscales;

            for (int i = 0; i < nscales; ++i) {
                scales[nscales - i - 1] = 2.0f + (float)pow(10.f, (i * size_step) / log(10.f));
            }
        } else if (mode == 2) {
            size_step = (float) log(s - 2.0f) / (float) nscales;

            for (int i = 0; i < nscales; ++i) {
                scales[i] = s - (float)pow(10.f, (i * size_step) / log(10.f));
            }
        } else if (mode == 3) {
            size_step = (float) log(s - 2.0f) / (float) nscales;

            for (int i = 0; i < nscales; ++i) {
                scales[i] = high * s - (float)pow(10.f, (i * size_step) / log(10.f));
            }
        }
    }
}

void mean_stddv2(float **dst, float &mean, float &stddv, int W_L, int H_L, float &maxtr, float &mintr)
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

        for (int i = 0; i < H_L; i++)
            for (int j = 0; j < W_L; j++) {
                sum += static_cast<double>(dst[i][j]);
                vsquared += rtengine::SQR<double>(dst[i][j]);

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
    mean = sum / (double)(W_L * H_L);
    vsquared /= (double) W_L * H_L;
    stddv = vsquared - rtengine::SQR<double>(mean);
    stddv = std::sqrt(stddv);
}

}


namespace rtengine
{

void RawImageSource::MSR(float** luminance, float** originalLuminance, float **exLuminance, const LUTf& mapcurve, bool mapcontlutili, int width, int height, const procparams::RetinexParams &deh, const RetinextransmissionCurve & dehatransmissionCurve, const RetinexgaintransmissionCurve & dehagaintransmissionCurve, float &minCD, float &maxCD, float &mini, float &maxi, float &Tmean, float &Tsigma, float &Tmin, float &Tmax)
{
    BENCHFUN

    if (!deh.enabled) {
        return;
    }

    constexpr float eps = 2.f;
    const bool useHsl = deh.retinexcolorspace == "HSLLOG";
    const bool useHslLin = deh.retinexcolorspace == "HSLLIN";
    const float offse = deh.offs; //def = 0  not use
    const int iter = deh.iter;
    const int gradient = deh.scal;
    int scal = deh.skal;
    const int nei = 2.8f * deh.neigh; //def = 220
    const float vart = deh.vart / 100.f;//variance
    const float gradvart = deh.grad;
    const float gradstr = deh.grads;
    const float strength = deh.str / 100.f; // Blend with original L channel data
    float limD = deh.limd;
    limD = pow(limD, 1.7f);//about 2500 enough
    limD *= useHslLin ? 10.f : 1.f;
    const float ilimD = 1.f / limD;
    const float hig = deh.highl / 100.f;

    const int H_L = height;
    const int W_L = width;

    constexpr float elogt = 2.71828f;
    bool lhutili = false;

    FlatCurve* shcurve = new FlatCurve(deh.lhcurve); //curve L=f(H)

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
    float high = bbhi + aahi * (float) deh.highl;

    for (int it = 1; it < iter + 1; it++) { //iter nb max of iterations
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
        } else if (gradient == -1) {
            grad = -0.125f * it + 1.125f;
            sc = 3.f;
        }

        if (iter == 1) {
            sc = scal;
        } else {
            //adjust sc in function of choice of scale by user if iterations
            if (scal < 3) {
                sc -= 1;

                if (sc < 1.f) {//avoid 0
                    sc = 1.f;
                }
            } else if (scal > 4) {
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

        scal = round(sc);
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

        const float strengthx = ks * strength;

        constexpr auto maxRetinexScales = 8;
        float RetinexScales[maxRetinexScales];

        retinex_scales(RetinexScales, scal, moderetinex, nei / grad, high);

        const int shHighlights = deh.highlights;
        const int shShadows = deh.shadows;

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

        std::unique_ptr<JaggedArray<float>> srcBuffer(new JaggedArray<float>(W_L, H_L));
        float** src = *(srcBuffer.get());

#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = 0; i < H_L; i++)
            for (int j = 0; j < W_L; j++) {
                src[i][j] = luminance[i][j] + eps;
                luminance[i][j] = 0.f;
            }

        JaggedArray<float> out(W_L, H_L);
        JaggedArray<float>& tran = out; // tran and out can safely use the same buffer

        const float logBetaGain = xlogf(16384.f);
        float pond = logBetaGain / (float) scal;

        if (!useHslLin) {
            pond /= log(elogt);
        }

        std::unique_ptr<SHMap> shmap;

        if (((mapmet == 2 || mapmet == 3 || mapmet == 4) && it == 1)) {
            shmap.reset(new SHMap(W_L, H_L));
        }

        std::unique_ptr<float[]> buffer;

        if (mapmet > 0) {
            buffer.reset(new float[W_L * H_L]);
        }

        for (int scale = scal - 1; scale >= 0; --scale) {
            if (scale == scal - 1) {
                gaussianBlur(src, out, W_L, H_L, RetinexScales[scale], true);
            } else { // reuse result of last iteration
                // out was modified in last iteration => restore it
                if ((((mapmet == 2 && scale > 1) || mapmet == 3 || mapmet == 4) || (mapmet > 0 && mapcontlutili)) && it == 1) {
#ifdef _OPENMP
                    #pragma omp parallel for
#endif

                    for (int i = 0; i < H_L; i++) {
                        for (int j = 0; j < W_L; j++) {
                            out[i][j] = buffer[i * W_L + j];
                        }
                    }
                }

                gaussianBlur(out, out, W_L, H_L, sqrtf(SQR(RetinexScales[scale]) - SQR(RetinexScales[scale + 1])), true);
            }

            if ((((mapmet == 2 && scale > 2) || mapmet == 3 || mapmet == 4) || (mapmet > 0 && mapcontlutili)) && it == 1 && scale > 0) {
                // out will be modified => store it for use in next iteration.
#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int i = 0; i < H_L; i++) {
                    for (int j = 0; j < W_L; j++) {
                        buffer[i * W_L + j] = out[i][j];
                    }
                }
            }

            int h_th = 0;
            int s_th = 0;

            if (((mapmet == 2 && scale > 2) || mapmet == 3 || mapmet == 4) && it == 1) {
                shmap->updateL(out, shradius, true, 1);
                h_th = shmap->max_f - deh.htonalwidth * (shmap->max_f - shmap->avg) / 100;
                s_th = deh.stonalwidth * (shmap->avg - shmap->min_f) / 100;
            }

            if (mapmet > 0 && mapcontlutili && it == 1) {
#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int i = 0; i < H_L; i++) {
                    for (int j = 0; j < W_L; j++) {
                        out[i][j] = mapcurve[2.f * out[i][j]];
                    }
                }

            }

            if (((mapmet == 2 && scale > 2) || mapmet == 3 || mapmet == 4) && it == 1) {
                const float hWeight = (100.f - shHighlights) / 100.f;
                const float sWeight = (100.f - shShadows) / 100.f;
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16)
#endif

                for (int i = 0; i < H_L; i++) {
                    for (int j = 0; j < W_L; j++) {
                        const float mapval = 1.f + shmap->map[i][j];
                        float factor;

                        if (mapval > h_th) {
                            factor = (h_th + hWeight * (mapval - h_th)) / mapval;
                        } else if (mapval < s_th) {
                            factor = (s_th - sWeight * (s_th - mapval)) / mapval;
                        } else {
                            factor = 1.f;
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
                const vfloat pondv = F2V(pond);
                const vfloat limMinv = F2V(ilimdx);
                const vfloat limMaxv = F2V(limdx);

                if (useHslLin) {
                    for (; j < W_L - 3; j += 4) {
                        STVFU(luminance[i][j], LVFU(luminance[i][j]) + pondv * vclampf(LVFU(src[i][j]) / LVFU(out[i][j]), limMinv, limMaxv));
                    }
                } else {
                    for (; j < W_L - 3; j += 4) {
                        STVFU(luminance[i][j], LVFU(luminance[i][j]) + pondv * xlogf(vclampf(LVFU(src[i][j]) / LVFU(out[i][j]), limMinv, limMaxv)));
                    }
                }

#endif

                if (useHslLin) {
                    for (; j < W_L; j++) {
                        luminance[i][j] +=  pond * LIM(src[i][j] / out[i][j], ilimdx, limdx);
                    }
                } else {
                    for (; j < W_L; j++) {
                        luminance[i][j] +=  pond * xlogf(LIM(src[i][j] / out[i][j], ilimdx, limdx)); //  /logt ?
                    }
                }
            }
        }

        srcBuffer.reset();

        float mean = 0.f;
        float stddv = 0.f;
        // I call mean_stddv2 instead of mean_stddv ==> logBetaGain

        float maxtr, mintr;
        mean_stddv2(luminance, mean, stddv, W_L, H_L, maxtr, mintr);
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
            #pragma omp parallel for schedule(dynamic,16)
#endif

            for (int i = 0; i < H_L; i++) {
                for (int j = 0; j < W_L; j++) { //for mintr to maxtr evalate absciss in function of original transmission
                    float absciss;

                    if (LIKELY(fabsf(luminance[i][j] - mean) < stddv)) {
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
                JaggedArray<float> tmL(W_L, H_L);
                constexpr int borderL = 1;

#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int i = borderL; i < H_L - borderL; i++) {
                    for (int j = borderL; j < W_L - borderL; j++) {
                        tmL[i][j] = median(luminance[i][j], luminance[i - 1][j], luminance[i + 1][j], luminance[i][j + 1], luminance[i][j - 1], luminance[i - 1][j - 1], luminance[i - 1][j + 1], luminance[i + 1][j - 1], luminance[i + 1][j + 1]); //3x3
                    }
                }

#ifdef _OPENMP
                #pragma omp parallel for
#endif

                for (int i = borderL; i < H_L - borderL; i++) {
                    for (int j = borderL; j < W_L - borderL; j++) {
                        luminance[i][j] = tmL[i][j];
                    }
                }
            }

            // I call mean_stddv2 instead of mean_stddv ==> logBetaGain
            //mean_stddv( luminance, mean, stddv, W_L, H_L, 1.f, maxtr, mintr);
            mean_stddv2(luminance, mean, stddv, W_L, H_L, maxtr, mintr);
        }

        constexpr float epsil = 0.1f;

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

        if (!delta) {
            delta = 1.0f;
        }

        // coeff for auto    "transmission" with 2 sigma #95% data
        const float aza = 16300.f / (2.f * stddv);
        const float azb = -aza * (mean - 2.f * stddv);
        const float bza = 16300.f / (2.f * stddv);
        const float bzb = 16300.f - bza * (mean);

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

        mean_stddv2(luminance, mean, stddv, W_L, H_L, maxtr, mintr);
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

        const float cdfactor = 32768.f / delta;
        maxCD = -9999999.f;
        minCD = 9999999.f;

#ifdef _OPENMP
        #pragma omp parallel for reduction(max:maxCD) reduction(min:minCD) schedule(dynamic, 16)
#endif

        for (int i = 0; i < H_L; i ++) {
            for (int j = 0; j < W_L; j++) {
                float gan;

                if (dehagaintransmissionCurve && mean != 0.f && stddv != 0.f) {
                    float absciss;

                    if (LIKELY(fabsf(luminance[i][j] - mean) < stddv)) {
                        absciss = asig * luminance[i][j] + bsig;
                    } else if (luminance[i][j] >= mean) {
                        absciss = amax * luminance[i][j] + bmax;
                    } else { /*if(luminance[i][j] <= mean - stddv)*/
                        absciss = amin * luminance[i][j] + bmin;
                    }


                    //    float cd = cdfactor * ( luminance[i][j]  - mini ) + offse;
                    // TODO : move multiplication by 2.f inside the curve
                    gan = 2.f * dehagaintransmissionCurve[absciss]; //new gain function transmission
                } else {
                    gan = 0.5f;
                }

                const float cd = gan * cdfactor * luminance[i][j] + offse;

                maxCD = cd > maxCD ? cd : maxCD;
                minCD = cd < minCD ? cd : minCD;

                float str = strengthx;

                if (lhutili && it == 1) { // S=f(H)
                    {
                        const float HH = exLuminance[i][j];
                        float valparam;

                        if (useHsl || useHslLin) {
                            valparam = shcurve->getVal(HH) - 0.5;
                        } else {
                            valparam = shcurve->getVal(Color::huelab_to_huehsv2(HH)) - 0.5;
                        }

                        str *= (1.f + 2.f * valparam);
                    }
                }

                if (higplus && exLuminance[i][j] > 65535.f * hig) {
                    str *= hig;
                }

                if (viewmet == 0) {
                    luminance[i][j] = intp(str, LIM(cd, 0.f, 32768.f), originalLuminance[i][j]);
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
                } else { /*if (viewmet == 3) */
                    luminance[i][j] = 1000.f + tran[i][j] * 700.f;    //arbitrary values to help display log values which are between -20 to + 30 - usage values -4 + 5
                }

            }
        }

        Tmean = mean;
        Tsigma = stddv;
        Tmin = mintr;
        Tmax = maxtr;

        if (shcurve) {
            delete shcurve;
            shcurve = nullptr;
        }
    }
}


void ImProcFunctions::maskforretinex(int sp, int before, float ** luminance, float ** out, int W_L, int H_L, int skip,
                                     const LocCCmaskCurve & locccmasretiCurve, bool &lcmasretiutili, const  LocLLmaskCurve & locllmasretiCurve, bool &llmasretiutili, const  LocHHmaskCurve & lochhmasretiCurve, bool & lhmasretiutili,
                                     int llretiMask, bool retiMasktmap, bool retiMask, float rad, float lap, bool pde, float gamm, float slop, float chro, float blend,
                                     LUTf & lmaskretilocalcurve, bool & localmaskretiutili,
                                     LabImage * bufreti, LabImage * bufmask, LabImage * buforig, LabImage * buforigmas, bool multiThread,
                                     bool delt, const float hueref, const float chromaref, const float lumaref,
                                     float maxdE, float mindE, float maxdElim,  float mindElim, float iterat, float limscope, int scope, float balance, float balanceh, float lumask)
{
    array2D<float> loctemp(W_L, H_L);
    array2D<float> ble(W_L, H_L);
    array2D<float> blechro(W_L, H_L);
    array2D<float> hue(W_L, H_L);
    array2D<float> guid(W_L, H_L);
    std::unique_ptr<LabImage> bufmaskblurreti;
    bufmaskblurreti.reset(new LabImage(W_L, H_L));
    std::unique_ptr<LabImage> bufmaskorigreti;
    bufmaskorigreti.reset(new LabImage(W_L, H_L));
    std::unique_ptr<LabImage> bufprov;
    bufprov.reset(new LabImage(W_L, H_L));

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif

    for (int y = 0; y < H_L; y++) {
        for (int x = 0; x < W_L; x++) {
            if (before == 1 && retiMasktmap) {
                loctemp[y][x] = LIM(luminance[y][x], 0.f, 32768.f);
            } else if (before == 0 && retiMasktmap) {
                loctemp[y][x] = out[y][x];
            } else {
                loctemp[y][x] = bufreti->L[y][x];
            }
        }
    }

    float chromult = 1.f - 0.01f * chro;
//fab
    float fab = 50.f;
    float meanfab = 0.f;
    const int nbfab = W_L * H_L;

    double sumab = 0.0;

#ifdef _OPENMP
        #pragma omp parallel for reduction(+:sumab)
#endif

    for (int y = 0; y < H_L; y++) {
        for (int x = 0; x < W_L; x++) {
            sumab += fabs(bufreti->a[y][x]);
            sumab += fabs(bufreti->b[y][x]);
        }
    }

    meanfab = sumab / (2.f * nbfab);
    double som = 0.0;

#ifdef _OPENMP
        #pragma omp parallel for reduction(+:som)
#endif

    for (int y = 0; y < H_L; y++) {
       for (int x = 0; x < W_L; x++) {
            som += SQR(fabs(bufreti->a[y][x]) - meanfab) + SQR(fabs(bufreti->b[y][x]) - meanfab);
       }
    }
    const float multsigma = (chro >= 0.f ? 0.035f : 0.018f) * chro + 1.f;
    const float stddv = sqrt(som / nbfab);
    fab = meanfab + multsigma * stddv;

    if (fab <= 0.f) {
        fab = 50.f;
    }
//end fab

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif

    for (int ir = 0; ir < H_L; ir++) {
        for (int jr = 0; jr < W_L; jr++) {
            float kmaskLexp = 0;
            float kmaskCH = 0;

            if (locllmasretiCurve  && llmasretiutili) {
                float ligh = loctemp[ir][jr] / 32768.f;
                kmaskLexp = 32768.f * LIM01(1.f - locllmasretiCurve[500.f * ligh]);

            }


            if (locllmasretiCurve  && llmasretiutili && retiMasktmap) {
            }

            if (llretiMask != 4) {
                if (locccmasretiCurve && lcmasretiutili) {
                    float chromask = 0.0001f + sqrt(SQR((bufreti->a[ir][jr]) / fab) + SQR((bufreti->b[ir][jr]) / fab));
                    kmaskCH = LIM01(1.f - locccmasretiCurve[500.f *  chromask]);
                }
            }

            if (lochhmasretiCurve && lhmasretiutili) {
                float huema = xatan2f(bufreti->b[ir][jr], bufreti->a[ir][jr]);
                float h = Color::huelab_to_huehsv2(huema);
                h += 1.f / 6.f;

                if (h > 1.f) {
                    h -= 1.f;
                }

                float valHH = LIM01(1.f - lochhmasretiCurve[500.f *  h]);

                if (llretiMask != 4) {
                    kmaskCH += chromult * valHH;
                }

                kmaskLexp += 32768.f * valHH;
            }

            bufmaskblurreti->L[ir][jr] = kmaskLexp;
            bufmaskblurreti->a[ir][jr] = kmaskCH;
            bufmaskblurreti->b[ir][jr] = kmaskCH;
            ble[ir][jr] = bufmaskblurreti->L[ir][jr] / 32768.f;
            hue[ir][jr] = xatan2f(bufmaskblurreti->b[ir][jr], bufmaskblurreti->a[ir][jr]);
            float chromah = sqrt(SQR(bufmaskblurreti->b[ir][jr]) + SQR(bufmaskblurreti->a[ir][jr]));
            blechro[ir][jr] = chromah / 32768.f;
            guid[ir][jr] = bufreti->L[ir][jr] / 32768.f;
            bufprov->L[ir][jr] = bufmaskblurreti->L[ir][jr];
            bufprov->a[ir][jr] = bufmaskblurreti->a[ir][jr];
            bufprov->b[ir][jr] = bufmaskblurreti->b[ir][jr];

        }
    }

    if (rad != 0.f) {
//        guidedFilter(guid, ble, ble, rad * 10.f / skip, 0.001, multiThread, 4);
            float blur = rad;
            blur = blur < 0.f ? -1.f / blur : 1.f + blur;
            int r1 = max(int(4 / skip * blur + 0.5), 1);
            int r2 = max(int(25 / skip * blur + 0.5), 1);

            double epsilmax = 0.0005;
            double epsilmin = 0.00001;

            double aepsil = (epsilmax - epsilmin) / 90.f;
            double bepsil = epsilmax - 100.f * aepsil;
            double epsil = aepsil * 0.1 * rad + bepsil;
            if (rad < 0.f) {
                epsil = 0.001;
            }
            rtengine::guidedFilter(guid, blechro, blechro, r1, epsil, multiThread);
            rtengine::guidedFilter(guid, ble, ble, r2, 0.2 * epsil, multiThread);

    }

    LUTf lutTonemaskreti(65536);
    calcGammaLut(gamm, slop, lutTonemaskreti);
    float radiusb = 1.f / skip;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif

    for (int ir = 0; ir < H_L; ir++)
        for (int jr = 0; jr < W_L; jr++) {
            float L_;
            float2 sincosval = xsincosf(hue[ir][jr]);
            bufmaskblurreti->L[ir][jr] = LIM01(ble[ir][jr]) * 32768.f;
            L_ = 2.f * bufmaskblurreti->L[ir][jr];
            bufmaskblurreti->L[ir][jr] = lutTonemaskreti[L_];
            bufmaskblurreti->a[ir][jr] = 32768.f * sincosval.y * blechro[ir][jr];
            bufmaskblurreti->b[ir][jr] = 32768.f * sincosval.x * blechro[ir][jr];
        }

    if (lmaskretilocalcurve && localmaskretiutili) {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

        for (int ir = 0; ir < H_L; ir++)
            for (int jr = 0; jr < W_L; jr++) {
                bufmaskblurreti->L[ir][jr] = 0.5f * lmaskretilocalcurve[2.f * bufmaskblurreti->L[ir][jr]];
            }
    }


    if (delt) {
        float *rdE[H_L] ALIGNED16;
        float *rdEBuffer = new float[H_L * W_L];

        for (int i = 0; i < H_L; i++) {
            rdE[i] = &rdEBuffer[i * W_L];
        }

        deltaEforMask(rdE, W_L, H_L, bufreti, hueref, chromaref, lumaref, maxdE, mindE, maxdElim, mindElim, iterat, limscope, scope, balance, balanceh);
        // printf("rde1=%f rde2=%f\n", rdE[1][1], rdE[100][100]);
        std::unique_ptr<LabImage> delta(new LabImage(W_L, H_L));
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

        for (int ir = 0; ir < H_L; ir++)
            for (int jr = 0; jr < W_L; jr++) {
                delta->L[ir][jr] = bufmaskblurreti->L[ir][jr] - bufprov->L[ir][jr];
                delta->a[ir][jr] = bufmaskblurreti->a[ir][jr] - bufprov->a[ir][jr];
                delta->b[ir][jr] = bufmaskblurreti->b[ir][jr] - bufprov->b[ir][jr];

                bufmaskblurreti->L[ir][jr] = bufprov->L[ir][jr] + rdE[ir][jr] * delta->L[ir][jr];
                bufmaskblurreti->a[ir][jr] = bufprov->a[ir][jr] + rdE[ir][jr] * delta->a[ir][jr];
                bufmaskblurreti->b[ir][jr] = bufprov->b[ir][jr] + rdE[ir][jr] * delta->b[ir][jr];
            }

        delete [] rdEBuffer;

    }


    if (lap > 0.f) {
        float *datain = new float[H_L * W_L];
        float *data_tmp = new float[H_L * W_L];

#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int y = 0; y < H_L; y++) {
            for (int x = 0; x < W_L; x++) {
                datain[y * W_L + x] =  bufmaskblurreti->L[y][x];
            }
        }

        if (!pde) {
            ImProcFunctions::discrete_laplacian_threshold(data_tmp, datain, W_L, H_L, 200.f * lap);
        } else {
            ImProcFunctions::retinex_pde(datain, data_tmp, W_L, H_L, 12.f * lap, 1.f, nullptr, 0, 0, 1);
        }

#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int y = 0; y < H_L; y++) {
            for (int x = 0; x < W_L; x++) {
                bufmaskblurreti->L[y][x] = data_tmp[y * W_L + x];
            }
        }

        delete [] datain;
        delete [] data_tmp;

    }

//blend
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        gaussianBlur(bufmaskblurreti->L, bufmaskorigreti->L, W_L, H_L, radiusb);
        gaussianBlur(bufmaskblurreti->a, bufmaskorigreti->a, W_L, H_L, 1.f + (0.5f * rad) / skip);
        gaussianBlur(bufmaskblurreti->b, bufmaskorigreti->b, W_L, H_L, 1.f + (0.5f * rad) / skip);
    }

    float modr = 0.01f * (float) blend;

    if (llretiMask != 3 && retiMask) {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

        for (int y = 0; y < H_L; y++) {
            for (int x = 0; x < W_L; x++) {
                if (before == 0 && retiMasktmap) {
                    out[y][x] += fabs(modr) * bufmaskorigreti->L[y][x];
                    out[y][x] = LIM(out[y][x], 0.f, 100000.f);
                } else {
                    bufreti->L[y][x] +=  bufmaskorigreti->L[y][x] * modr;
                    bufreti->L[y][x] = CLIPLOC(bufreti->L[y][x]);

                }

                bufreti->a[y][x] *= (1.f + bufmaskorigreti->a[y][x] * modr * (1.f + 0.01f * chro));
                bufreti->b[y][x] *= (1.f + bufmaskorigreti->b[y][x] * modr * (1.f + 0.01f * chro));
                bufreti->a[y][x] = CLIPC(bufreti->a[y][x]);
                bufreti->b[y][x] = CLIPC(bufreti->b[y][x]);
            }
        }


    }

    if (!retiMasktmap  && retiMask) { //new original blur mask for deltaE
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

        for (int y = 0; y < H_L; y++) {
            for (int x = 0; x < W_L; x++) {

                buforig->L[y][x] += (modr * bufmaskorigreti->L[y][x]);
                buforig->a[y][x] *= (1.f + modr * bufmaskorigreti->a[y][x]);
                buforig->b[y][x] *= (1.f + modr * bufmaskorigreti->b[y][x]);

                buforig->L[y][x] = CLIP(buforig->L[y][x]);
                buforig->a[y][x] = CLIPC(buforig->a[y][x]);
                buforig->b[y][x] = CLIPC(buforig->b[y][x]);

                buforig->L[y][x] = CLIP(buforig->L[y][x] - bufmaskorigreti->L[y][x]);
                buforig->a[y][x] = CLIPC(buforig->a[y][x] * (1.f - bufmaskorigreti->a[y][x]));
                buforig->b[y][x] = CLIPC(buforig->b[y][x] * (1.f - bufmaskorigreti->b[y][x]));
            }
        }

        float radius = 3.f / skip;

#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            gaussianBlur(buforig->L, buforigmas->L, W_L, H_L, radius);
            gaussianBlur(buforig->a, buforigmas->a, W_L, H_L, radius);
            gaussianBlur(buforig->b, buforigmas->b, W_L, H_L, radius);
        }

    }


    if (llretiMask == 3) {

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

        for (int y = 0; y < H_L; y++) {
            for (int x = 0; x < W_L; x++) {
                bufmask->L[y][x] = (lumask * 400.f) + CLIPLOC(bufmaskorigreti->L[y][x]);
                bufmask->a[y][x] = CLIPC(bufreti->a[y][x] * bufmaskorigreti->a[y][x]);
                bufmask->b[y][x] = CLIPC(bufreti->b[y][x] * bufmaskorigreti->b[y][x]);
            }
        }
    }

}



void ImProcFunctions::MSRLocal(int call, int sp, bool fftw, int lum, float** reducDE, LabImage * bufreti, LabImage * bufmask, LabImage * buforig, LabImage * buforigmas, float** luminance, const float* const *originalLuminance,
                               const int width, const int height, int bfwr, int bfhr, const procparams::LocallabParams &loc, const int skip, const LocretigainCurve &locRETgainCcurve, const LocretitransCurve &locRETtransCcurve,
                               const int chrome, const int scall, const float krad, float &minCD, float &maxCD, float &mini, float &maxi, float &Tmean, float &Tsigma, float &Tmin, float &Tmax,
                               const LocCCmaskCurve & locccmasretiCurve, bool &lcmasretiutili, const  LocLLmaskCurve & locllmasretiCurve, bool &llmasretiutili, const  LocHHmaskCurve & lochhmasretiCurve, bool & lhmasretiutili, int llretiMask,
                               LUTf & lmaskretilocalcurve, bool & localmaskretiutili,
                               LabImage * transformed, bool retiMasktmap, bool retiMask,
                               bool delt, const float hueref, const float chromaref, const float lumaref,
                               float maxdE, float mindE, float maxdElim,  float mindElim, float iterat, float limscope, int scope, float balance, float balanceh, float lumask)

{
    BENCHFUN
    bool py = true;

    if (py) {//enabled
        float mean, stddv, maxtr, mintr;
        mean = 0.f;
        stddv = 0.f;
        maxtr = 0.f;
        mintr = 0.f;
        float delta;
        constexpr float eps = 2.f;
        bool useHslLin = false;
        const float offse = loc.spots.at(sp).offs;
        const float chrT = (float)(loc.spots.at(sp).chrrt) / 100.f;
        const int scal = (loc.spots.at(sp).scalereti);
        float vart = loc.spots.at(sp).vart / 100.f;//variance
        const float strength = loc.spots.at(sp).str / 100.f; // Blend with original L channel data
        const float dar = loc.spots.at(sp).darkness;
        const float lig = loc.spots.at(sp).lightnessreti;
        float value = pow(strength, 0.4f);
        float value_1 = pow(strength, 0.3f);
        bool logli = loc.spots.at(sp).loglin;
        float limD = loc.spots.at(sp).limd;//10.f
        limD = pow(limD, 1.7f);  //about 2500 enough
        float ilimD = 1.f / limD;
        float threslum = loc.spots.at(sp).limd;
        const float elogt = 2.71828f;

        if (!logli) {
            useHslLin = true;
        }

        //empirical skip evaluation : very difficult  because quasi all parameters interfere
        //to test on several images
        int nei = (int)(krad * loc.spots.at(sp).neigh);
        // printf("neigh=%i\n", nei);
        //several test to find good values ???!!!
        //very difficult to do because 4 factor are correlate with skip and cannot been solved easily
        // size of spots
        // radius - neigh
        // scal
        // variance vart
        //not too bad proposition
        float divsca = 1.f;

        if (scal >= 3) {
            divsca = sqrt(scal / 3.f);
        }

        if (skip >= 4) {
            //nei = (int)(0.1f * nei + 2.f);     //not too bad
            nei = (int)(nei / (1.5f * skip)) / divsca;
            vart *= sqrt(skip);
        } else if (skip > 1 && skip < 4) {
            //nei = (int)(0.3f * nei + 2.f);
            nei = (int)(nei / skip) / divsca;
            vart *= sqrt(skip);
        }

        int moderetinex = 0;

        if (loc.spots.at(sp).retinexMethod == "uni") {
            moderetinex = 0;
        } else if (loc.spots.at(sp).retinexMethod == "low") {
            moderetinex = 1;
        } else if (loc.spots.at(sp).retinexMethod == "high") {
            moderetinex = 2;
        }

        const float high = 0.f; // Dummy to pass to retinex_scales(...)

        constexpr auto maxRetinexScales = 10;
        float RetinexScales[maxRetinexScales];

        retinex_scales(RetinexScales, scal, moderetinex, nei, high);


        const int H_L = height;
        const int W_L = width;
        std::unique_ptr<JaggedArray<float>> srcBuffer(new JaggedArray<float>(W_L, H_L));
        float** src = *(srcBuffer.get());


#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = 0; i < H_L; i++)
            for (int j = 0; j < W_L; j++) {
                src[i][j] = luminance[i][j] + eps;
                luminance[i][j] = 0.f;
            }

        JaggedArray<float> out(W_L, H_L);

        float clipt = loc.spots.at(sp).cliptm;

        const float logBetaGain = xlogf(16384.f);
        float pond = logBetaGain / (float) scal;

        if (!useHslLin) {
            pond /= log(elogt);
        }

        float kr = 1.f;//on FFTW
        float kg = 1.f;//on Gaussianblur
        std::unique_ptr<float[]> buffer;
        buffer.reset(new float[W_L * H_L]);

        for (int scale = scal - 1; scale >= 0; --scale) {
            //    printf("retscale=%f scale=%i \n", mulradiusfftw * RetinexScales[scale], scale);
            //emprical adjustment between FFTW radius and Gaussainblur
            //under 50 ==> 10.f
            // 400 ==> 1.f
            float sigm = 1.f;

            if (settings->fftwsigma == false) { //empirical formula
                sigm = RetinexScales[scale];
                float ak = -9.f / 350.f;
                float bk = 10.f - 50.f * ak;
                kr = ak * sigm + bk;

                if (sigm < 50.f) {
                    kr = 10.f;
                }

                //above 400 at 5000 ==> 20.f
                if (sigm > 400.f) { //increase ==> 5000
                    float ka = 19.f / 4600.f;
                    float kb = 1.f - 400 * ka;
                    kr = ka * sigm + kb;
                    float kga = -0.14f / 4600.f;//decrease
                    float kgb = 1.f - 400.f * kga;
                    kg = kga * sigm + kgb;

                    if (sigm > 5000.f) {
                        kr = ka * 5000.f + kb;
                        kg = kga * 5000.f + kgb;
                    }

                }
            } else {//sigma *= sigma
                kg = 1.f;
                kr = sigm;
            }
            printf("call=%i\n", call);
            if (!fftw) { // || (fftw && call != 2)) {
                if (scale == scal - 1) {
                    gaussianBlur(src, out, W_L, H_L, kg * RetinexScales[scale], true);
                } else { // reuse result of last iteration
                    // out was modified in last iteration => restore it
                    gaussianBlur(out, out, W_L, H_L, sqrtf(SQR(kg * RetinexScales[scale]) - SQR(kg * RetinexScales[scale + 1])), true);
                }
            } else {
                if (scale == scal - 1) {
                    if (settings->fftwsigma == false) { //empirical formula
                        ImProcFunctions::fftw_convol_blur2(src, out, bfwr, bfhr, (kr * RetinexScales[scale]), 0, 0);
                    } else {
                        ImProcFunctions::fftw_convol_blur2(src, out, bfwr, bfhr, (SQR(RetinexScales[scale])), 0, 0);
                    }
                } else { // reuse result of last iteration
                    // out was modified in last iteration => restore it
                    if (settings->fftwsigma == false) { //empirical formula
                        ImProcFunctions::fftw_convol_blur2(out, out, bfwr, bfhr, sqrtf(SQR(kr * RetinexScales[scale]) - SQR(kr * RetinexScales[scale + 1])), 0, 0);
                    } else {
                        ImProcFunctions::fftw_convol_blur2(out, out, bfwr, bfhr, (SQR(RetinexScales[scale]) - SQR(RetinexScales[scale + 1])), 0, 0);
                    }
                }
            }

            if (scale == 1) { //equalize last scale with darkness and lightness of course acts on TM!
                if (dar != 1.f || lig != 1.f) {

#ifdef _OPENMP
                    #pragma omp parallel for
#endif

                    for (int y = 0; y < H_L; ++y) {
                        for (int x = 0; x < W_L; ++x) {
                            float buf = (src[y][x] - out[y][x]) * value;
                            buf *= (buf > 0.f) ? lig : dar;
                            out[y][x] = LIM(out[y][x] + buf, -100000.f, 100000.f);
                        }
                    }
                }
                
            }

#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int i = 0; i < H_L; i++) {
                int j = 0;

#ifdef __SSE2__
                const vfloat pondv = F2V(pond);
                const vfloat limMinv = F2V(ilimD);
                const vfloat limMaxv = F2V(limD);

                if (useHslLin) {
                    for (; j < W_L - 3; j += 4) {
                        STVFU(luminance[i][j], LVFU(luminance[i][j]) + pondv * vclampf(LVFU(src[i][j]) / LVFU(out[i][j]), limMinv, limMaxv));
                    }
                } else {
                    for (; j < W_L - 3; j += 4) {
                        STVFU(luminance[i][j], LVFU(luminance[i][j]) + pondv * xlogf(vclampf(LVFU(src[i][j]) / LVFU(out[i][j]), limMinv, limMaxv)));
                    }
                }

#endif

                if (useHslLin) {
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

//        srcBuffer.reset();


        if (scal == 1) {//only if user select scal = 1

            float kval = 1.f;
#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int y = 0; y < H_L; ++y) {
                for (int x = 0; x < W_L; ++x) {
                    float threslow = threslum * 163.f;

                    if (src[y][x] < threslow) {
                        kval = src[y][x] / threslow;
                    }
                }
            }


#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int y = 0; y < H_L; ++y) {
                for (int x = 0; x < W_L; ++x) {
                    float buf = (src[y][x] - out[y][x]) * value_1;
                    buf *= (buf > 0.f) ? lig : dar;
                    luminance[y][x] = LIM(src[y][x] + (1.f + kval) * buf, -32768.f, 32768.f);
                }
            }

            double avg = 0.f;
            int ng = 0;

#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int i = 0; i < H_L; i++) {
                for (int j = 0; j < W_L; j++) {
                    avg += luminance[i][j];
                    ng++;
                }
            }

            avg /= ng;
            avg /= 32768.f;
            avg = LIM01(avg);
            float contreal = 0.5f * vart;
            DiagonalCurve reti_contrast({
                DCT_NURBS,
                0, 0,
                avg - avg * (0.6 - contreal / 250.0), avg - avg * (0.6 + contreal / 250.0),
                avg + (1 - avg) * (0.6 - contreal / 250.0), avg + (1 - avg) * (0.6 + contreal / 250.0),
                1, 1
            });

#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int i = 0; i < H_L; i++)
                for (int j = 0; j < W_L; j++) {
                    float buf = LIM01(luminance[i][j] / 32768.f);
                    buf = reti_contrast.getVal(buf);
                    buf *= 32768.f;
                    luminance[i][j] = buf;
                }
      
        }

        srcBuffer.reset();

        float str = strength * (chrome == 0 ? 1.f : 0.8f * (chrT - 0.4f));
        const float maxclip = (chrome == 0 ? 32768.f : 50000.f);

        if (scal != 1) {
            mean = 0.f;
            stddv = 0.f;

            mean_stddv2(luminance, mean, stddv, W_L, H_L, maxtr, mintr);
            // printf("mean=%f std=%f delta=%f maxtr=%f mintr=%f\n", mean, stddv, delta, maxtr, mintr);

            if (locRETtransCcurve && mean != 0.f && stddv != 0.f) { //if curve
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

                    for (int i = 0; i < H_L; i++)
                        for (int j = 0; j < W_L; j++) { //for mintr to maxtr evalate absciss in function of original transmission
                            if (LIKELY(fabsf(luminance[i][j] - mean) < stddv)) {
                                absciss = asig * luminance[i][j] + bsig;
                            } else if (luminance[i][j] >= mean) {
                                absciss = amax * luminance[i][j] + bmax;
                            } else { /*if(luminance[i][j] <= mean - stddv)*/
                                absciss = amin * luminance[i][j] + bmin;
                            }

                            //TODO : move multiplication by 4.f and subtraction of 1.f inside the curve
                            luminance[i][j] *= (-1.f + 4.f * locRETtransCcurve[absciss]); //new transmission

                        }
                }
            }

            mean = 0.f;
            stddv = 0.f;
            mean_stddv2(luminance, mean, stddv, W_L, H_L, maxtr, mintr);//new calculation of mean...

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
            //   printf("maxi=%f mini=%f mean=%f std=%f delta=%f maxtr=%f mintr=%f\n", maxi, mini, mean, stddv, delta, maxtr, mintr);

            if (!delta) {
                delta = 1.0f;
            }


            float *copylum[H_L] ALIGNED16;
            float *copylumBuffer = new float[H_L * W_L];

            for (int i = 0; i < H_L; i++) {
                copylum[i] = &copylumBuffer[i * W_L];
            }

            float cdfactor = (clipt * 32768.f) / delta;
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

            mean_stddv2(luminance, mean, stddv, W_L, H_L, maxtr, mintr);
//         printf("meanun=%f stdun=%f maxtr=%f mintr=%f\n", mean, stddv, maxtr, mintr);

            float asig = 0.f, bsig = 0.f, amax = 0.f, bmax = 0.f, amin = 0.f, bmin = 0.f;
            const bool hasRetGainCurve =  locRETgainCcurve && mean != 0.f && stddv != 0.f;

            if (hasRetGainCurve) { //if curve
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

                for (int i = 0; i < H_L; i ++)
                    for (int j = 0; j < W_L; j++) {

                        if (hasRetGainCurve) {
                            float absciss;

                            if (LIKELY(fabsf(luminance[i][j] - mean) < stddv)) {
                                absciss = asig * luminance[i][j] + bsig;
                            } else if (luminance[i][j] >= mean) {
                                absciss = amax * luminance[i][j] + bmax;
                            } else {
                                absciss = amin * luminance[i][j] + bmin;
                            }

                            gan = locRETgainCcurve[absciss]; //new gain function transmission
                        }

                        //but we don't update mean stddv for display only...
                        copylum[i][j] = gan * luminance[i][j];//update data for display
                        float cd = gan * cdfactor * luminance[i][j] + offse;

                        cdmax = cd > cdmax ? cd : cdmax;
                        cdmin = cd < cdmin ? cd : cdmin;
                        luminance[i][j] = intp(str * reducDE[i][j], clipretinex(cd, 0.f, maxclip), originalLuminance[i][j]);
                    }



#ifdef _OPENMP
                #pragma omp critical
#endif
                {
                    maxCD = maxCD > cdmax ? maxCD : cdmax;
                    minCD = minCD < cdmin ? minCD : cdmin;
                }
            }
            mean = 0.f;
            stddv = 0.f;

            mean_stddv2(copylum, mean, stddv, W_L, H_L, maxtr, mintr);
            delete [] copylumBuffer;
            copylumBuffer = nullptr;

//         printf("mean=%f std=%f maxtr=%f mintr=%f\n", mean, stddv, maxtr, mintr);

        } else {
#ifdef _OPENMP
            #pragma omp for schedule(dynamic,16)
#endif

            for (int i = 0; i < H_L; i ++)
                for (int j = 0; j < W_L; j++) {
                    luminance[i][j] =   LIM(luminance[i][j], 0.f, maxclip) * str + (1.f - str) * originalLuminance[i][j];

                }

        }

        float rad = loc.spots.at(sp).radmaskreti;
        float slop = loc.spots.at(sp).slomaskreti;
        float gamm = loc.spots.at(sp).gammaskreti;
        float blend = loc.spots.at(sp).blendmaskreti;
        float chro = loc.spots.at(sp).chromaskreti;
        float lap = loc.spots.at(sp).lapmaskreti;
        bool pde = params->locallab.spots.at(sp).laplac;

        if (lum == 1  && (llretiMask == 3 || llretiMask == 0 || llretiMask == 2 || llretiMask == 4)) { //only mask with luminance on last scale
            int before = 1;
            maskforretinex(sp, before, luminance, nullptr, W_L, H_L, skip,
                           locccmasretiCurve, lcmasretiutili, locllmasretiCurve, llmasretiutili, lochhmasretiCurve, lhmasretiutili, llretiMask, retiMasktmap, retiMask,
                           rad, lap, pde, gamm, slop, chro, blend,
                           lmaskretilocalcurve, localmaskretiutili,
                           bufreti, bufmask, buforig, buforigmas, multiThread,
                           delt, hueref, chromaref, lumaref,
                           maxdE, mindE, maxdElim, mindElim, iterat, limscope, scope, balance, balanceh, lumask
                          );
        }

        //mask does not interfered with data displayed

        Tmean = mean;
        Tsigma = stddv;
        Tmin = mintr;
        Tmax = maxtr;
    }
}
}
