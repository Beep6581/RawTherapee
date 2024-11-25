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
 */

#include "bilateral2.h"
#include "cieimage.h"
#include "gauss.h"
#include "improcfun.h"
#include "jaggedarray.h"
#include "labimage.h"
#include "opthelper.h"
#include "procparams.h"
#include "rt_algo.h"
#include "rt_math.h"
#include "settings.h"
#include "sleef.h"
#include "rtengine.h"
#include "color.h"

//#define BENCHMARK
#include "StopWatch.h"
#include "imagefloat.h"
#include "iccstore.h"

using namespace std;

namespace {

void sharpenHaloCtrl (float** luminance, float** blurmap, float** base, float** blend, int W, int H, const procparams::SharpeningParams &sharpenParam)
{

    const float scale = (100.f - sharpenParam.halocontrol_amount) * 0.01f;
    const float sharpFac = sharpenParam.amount * 0.01f;
    float** nL = base;

#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for (int i = 2; i < H - 2; i++) {
        float max1 = 0, max2 = 0, min1 = 0, min2 = 0;

        for (int j = 2; j < W - 2; j++) {
            // compute 3 iterations, only forward
            float np1 = 2.f * (nL[i - 2][j] + nL[i - 2][j + 1] + nL[i - 2][j + 2] + nL[i - 1][j] + nL[i - 1][j + 1] + nL[i - 1][j + 2] + nL[i]  [j] + nL[i]  [j + 1] + nL[i]  [j + 2]) / 27.f + nL[i - 1][j + 1] / 3.f;
            float np2 = 2.f * (nL[i - 1][j] + nL[i - 1][j + 1] + nL[i - 1][j + 2] + nL[i]  [j] + nL[i]  [j + 1] + nL[i]  [j + 2] + nL[i + 1][j] + nL[i + 1][j + 1] + nL[i + 1][j + 2]) / 27.f + nL[i]  [j + 1] / 3.f;
            float np3 = 2.f * (nL[i]  [j] + nL[i]  [j + 1] + nL[i]  [j + 2] + nL[i + 1][j] + nL[i + 1][j + 1] + nL[i + 1][j + 2] + nL[i + 2][j] + nL[i + 2][j + 1] + nL[i + 2][j + 2]) / 27.f + nL[i + 1][j + 1] / 3.f;

            // Max/Min of all these deltas and the last two max/min
            float maxn = rtengine::max(np1, np2, np3);
            float minn = rtengine::min(np1, np2, np3);
            float max_ = rtengine::max(max1, max2, maxn);
            float min_ = rtengine::min(min1, min2, minn);

            // Shift the queue
            max1 = max2;
            max2 = maxn;
            min1 = min2;
            min2 = minn;
            float labL = luminance[i][j];

            if (max_ < labL) {
                max_ = labL;
            }

            if (min_ > labL) {
                min_ = labL;
            }

            // deviation from the environment as measurement
            float diff = nL[i][j] - blurmap[i][j];

            constexpr float upperBound = 2000.f;  // WARNING: Duplicated value, it's baaaaaad !
            float delta = sharpenParam.threshold.multiply<float, float, float>(
                              rtengine::min(fabsf(diff), upperBound),   // X axis value = absolute value of the difference
                              sharpFac * diff               // Y axis max value = sharpening.amount * signed difference
                          );
            float newL = labL + delta;

            // applying halo control
            if (newL > max_) {
                newL = max_ + (newL - max_) * scale;
            } else if (newL < min_) {
                newL = min_ - (min_ - newL) * scale;
            }

            luminance[i][j] = intp(blend[i][j], newL, luminance[i][j]);
        }
    }
}

void dcdamping (float** aI, float** aO, float damping, int W, int H)
{

    const float dampingFac = -2.f / (damping * damping);

#ifdef __SSE2__
    vfloat Iv, Ov, Uv, zerov, onev, fourv, fivev, dampingFacv, Tv, Wv, Lv;
    zerov = _mm_setzero_ps();
    onev = F2V(1.f);
    fourv = F2V(4.f);
    fivev = F2V(5.f);
    dampingFacv = F2V(dampingFac);
#endif
#ifdef _OPENMP
    #pragma omp for
#endif

    for (int i = 0; i < H; i++) {
        int j = 0;
#ifdef __SSE2__

        for (; j < W - 3; j += 4) {
            Iv = LVFU(aI[i][j]);
            Ov = LVFU(aO[i][j]);
            Lv = xlogf(Iv / Ov);
            Wv = Ov - Iv;
            Uv = (Ov * Lv + Wv) * dampingFacv;
            Uv = vminf(Uv, onev);
            Tv = Uv * Uv;
            Tv = Tv * Tv;
            Uv = Tv * (fivev - Uv * fourv);
            Uv = (Wv / Iv) * Uv + onev;
            Uv = vselfzero(vmaskf_gt(Iv, zerov), Uv);
            Uv = vselfzero(vmaskf_gt(Ov, zerov), Uv);
            STVFU(aI[i][j], Uv);
        }

#endif

        for(; j < W; j++) {
            float I = aI[i][j];
            float O = aO[i][j];

            if (O <= 0.f || I <= 0.f) {
                aI[i][j] = 0.f;
                continue;
            }

            float U = (O * xlogf(I / O) - I + O) * dampingFac;
            U = rtengine::min(U, 1.0f);
            U = U * U * U * U * (5.f - U * 4.f);
            aI[i][j] = (O - I) / I * U + 1.f;
        }
    }
}

}

namespace rtengine
{

void ImProcFunctions::deconvsharpening (float** luminance, float** tmp, const float * const * blend, int W, int H, const procparams::SharpeningParams &sharpenParam, double Scale)
{
    if (sharpenParam.deconvamount == 0 && sharpenParam.blurradius < 0.25) {
        return;
    }
BENCHFUN
    JaggedArray<float> tmpI(W, H);

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < H; i++) {
        for(int j = 0; j < W; j++) {
            tmpI[i][j] = max(luminance[i][j], 0.f);
        }
    }

    JaggedArray<float>* blurbuffer = nullptr;

    if (sharpenParam.blurradius >= 0.25) {
        blurbuffer = new JaggedArray<float>(W, H);
        JaggedArray<float> &blur = *blurbuffer;
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            gaussianBlur(tmpI, blur, W, H, sharpenParam.blurradius);
#ifdef _OPENMP
            #pragma omp for
#endif
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    blur[i][j] = intp(blend[i][j], luminance[i][j], std::max(blur[i][j], 0.0f));
                }
            }
        }
    }
    const float damping = sharpenParam.deconvdamping / 5.0;
    const bool needdamp = sharpenParam.deconvdamping > 0;
    const double sigma = sharpenParam.deconvradius / Scale;
    const float amount = sharpenParam.deconvamount / 100.f;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        for (int k = 0; k < sharpenParam.deconviter; k++) {
            if (!needdamp) {
                // apply gaussian blur and divide luminance by result of gaussian blur
                gaussianBlur(tmpI, tmp, W, H, sigma, false, GAUSS_DIV, luminance);
            } else {
                // apply gaussian blur + damping
                gaussianBlur(tmpI, tmp, W, H, sigma);
                dcdamping(tmp, luminance, damping, W, H);
            }
            gaussianBlur(tmp, tmpI, W, H, sigma, false, GAUSS_MULT);
        } // end for

#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                luminance[i][j] = intp(blend[i][j] * amount, max(tmpI[i][j], 0.0f), luminance[i][j]);
            }
        }

        if (sharpenParam.blurradius >= 0.25) {
            JaggedArray<float> &blur = *blurbuffer;
#ifdef _OPENMP
        #pragma omp for
#endif
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    luminance[i][j] = intp(blend[i][j], luminance[i][j], max(blur[i][j], 0.0f));
                }
            }
        }
    } // end parallel
    delete blurbuffer;
}

void ImProcFunctions::deconvsharpeningloc (float** luminance, float** tmp, int W, int H, float** loctemp, int damp, double radi, int ite, int amo, int contrast, double blurrad, int sk)
{
    // BENCHFUN

    if (amo < 1) {
        return;
    }
    JaggedArray<float> blend(W, H);
    float contras = contrast / 100.f;
    buildBlendMask(luminance, blend, W, H, contras);


    JaggedArray<float> tmpI(W, H);


#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            tmpI[i][j] = max(luminance[i][j], 0.f);
        }
    }

    // calculate contrast based blend factors to reduce sharpening in regions with low contrast

    JaggedArray<float>* blurbuffer = nullptr;

    if (blurrad >= 0.25) {
        blurbuffer = new JaggedArray<float>(W, H);
        JaggedArray<float> &blur = *blurbuffer;
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            gaussianBlur(tmpI, blur, W, H, blurrad);
#ifdef _OPENMP
            #pragma omp for
#endif
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    blur[i][j] = intp(blend[i][j], luminance[i][j], std::max(blur[i][j], 0.0f));
                }
            }
        }
    }

    float damping = (float) damp / 5.0;
    bool needdamp = damp > 0;
    double sigma = radi / sk;
    const float amount = (float) amo / 100.f;

    if (sigma < 0.26f) {
        sigma = 0.26f;
    }

    int itera = ite;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        for (int k = 0; k < itera; k++) {
            if (!needdamp) {
                // apply gaussian blur and divide luminance by result of gaussian blur
            //    gaussianBlur (tmpI, tmp, W, H, sigma, nullptr, GAUSS_DIV, luminance);
                gaussianBlur(tmpI, tmp, W, H, sigma, false, GAUSS_DIV, luminance);
            } else {
                // apply gaussian blur + damping
                gaussianBlur (tmpI, tmp, W, H, sigma);
                dcdamping (tmp, luminance, damping, W, H);
            }

            gaussianBlur (tmp, tmpI, W, H, sigma, false, GAUSS_MULT);
        } // end for


#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < H; i++)
            for (int j = 0; j < W; j++) {
                loctemp[i][j] = intp(blend[i][j] * amount, max(tmpI[i][j], 0.0f), luminance[i][j]);
            }
            
        if (blurrad >= 0.25) {
            JaggedArray<float> &blur = *blurbuffer;
#ifdef _OPENMP
        #pragma omp for
#endif
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    loctemp[i][j] = intp(blend[i][j], loctemp[i][j], max(blur[i][j], 0.0f));
                }
            }
        }
            
    } // end parallel
    delete blurbuffer;


}

class CornerBoostMask {
public:
    CornerBoostMask(int ox, int oy, int width, int height, int latitude):
        ox_(ox), oy_(oy), w2_(width / 2), h2_(height / 2)
    {
        float radius = std::max(w2_, h2_);
        r2_ = (radius - radius * LIM01(float(latitude)/150.f)) / 2.f;
        sigma_ = 2.f * SQR(radius * 0.3f);
    }

    float operator()(int x, int y) const
    {
        int xx = x + ox_ - w2_;
        int yy = y + oy_ - h2_;
        float distance = std::sqrt(float(SQR(xx) + SQR(yy)));
        return 1.f - LIM01(xexpf((-SQR(std::max(distance - r2_, 0.f)) / sigma_)));
    }

private:
    int ox_;
    int oy_;
    int w2_;
    int h2_;
    float r2_;
    float sigma_;
};



void deconvsharpeningrgbloc(float **luminance, float **blend, char **impulse, int W, int H, double sigma, float amount, bool multiThread)
{
    if (amount <= 0) {
        return;
    }
BENCHFUN

    const int maxiter = 20;
    const float delta_factor = 0.2f;

    if (sigma < 0.2f) {
        return;
    }
    
    JaggedArray<float> tmp(W, H);
    JaggedArray<float> tmpI(W, H);
    JaggedArray<float> out(W, H);

    constexpr float offset = 1000.f;

#ifdef _OPENMP
#   pragma omp parallel for if (multiThread)
#endif
    for (int i = 0; i < H; i++) {
        for(int j = 0; j < W; j++) {
            luminance[i][j] += offset;
            tmpI[i][j] = std::max(luminance[i][j], 0.f);
            assert(std::isfinite(tmpI[i][j]));
            out[i][j] = RT_NAN;
        }
    }

    const auto get_output =
        [&](int i, int j) -> float
        {
            if (UNLIKELY(std::isnan(tmpI[i][j]))) {
                return luminance[i][j];
            }
            float b = impulse[i][j] ? 0.f : blend[i][j] * amount;
            return intp(b, std::max(tmpI[i][j], 0.0f), luminance[i][j]);
        };

    const auto check_stop =
        [&](int y, int x) -> void
        {
            if (LIKELY(std::isnan(out[y][x]))) {
                float l = luminance[y][x];
                float delta = l * delta_factor;
                if (UNLIKELY(std::abs(tmpI[y][x] - l) > delta)) {
                    out[y][x] = get_output(y, x);
                }
            }
        };

#ifdef _OPENMP
#   pragma omp parallel if (multiThread)
#endif
    {
        for (int k = 0; k < maxiter; k++) {
            gaussianBlur(tmpI, tmp, W, H, sigma);
            gaussianBlur(tmp, tmpI, W, H, sigma);
#ifdef _OPENMP
#           pragma omp for
#endif
            for (int y = 0; y < H; ++y) {
                for (int x = 0; x < W; ++x) {
                    check_stop(y, x);
                }
            }
        }

#ifdef _OPENMP
#       pragma omp for
#endif
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                float l = out[i][j];
                if (std::isnan(l)) {
                    l = get_output(i, j);
                }
                assert(std::isfinite(l));
                luminance[i][j] = std::max(l - offset, 0.f);
            }
        }
    }
}

bool checkForStop(float** tmpIThr, float** iterCheck, int fullTileSize, int border)
{
    for (int ii = border; ii < fullTileSize - border; ++ii) {
#ifdef __SSE2__
        for (int jj = border; jj < fullTileSize - border; jj += 4) {
            if (UNLIKELY(_mm_movemask_ps((vfloat)vmaskf_lt(LVFU(tmpIThr[ii][jj]), LVFU(iterCheck[ii - border][jj - border]))))) {
                return true;
            }
        }
#else
        for (int jj = border; jj < fullTileSize - border; ++jj) {
            if (tmpIThr[ii][jj] < iterCheck[ii - border][jj - border]) {
                return true;
            }
        }
#endif
    }
    return false;
}

void compute13x13kernel(float sigma, float kernel[13][13]) {

    const double temp = -2.f * rtengine::SQR(sigma);
    float sum = 0.f;
    for (int i = -6; i <= 6; ++i) {
        for (int j = -6; j <= 6; ++j) {
            if((rtengine::SQR(i) + rtengine::SQR(j)) <= rtengine::SQR(3.0 * 2.0)) {
                kernel[i + 6][j + 6] = std::exp((rtengine::SQR(i) + rtengine::SQR(j)) / temp);
                sum += kernel[i + 6][j + 6];
            } else {
                kernel[i + 6][j + 6] = 0.f;
            }
        }
    }

    for (int i = 0; i < 13; ++i) {
        for (int j = 0; j < 13; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

void compute9x9kernel(float sigma, float kernel[9][9]) {

    const double temp = -2.f * rtengine::SQR(sigma);
    float sum = 0.f;
    for (int i = -4; i <= 4; ++i) {
        for (int j = -4; j <= 4; ++j) {
            if((rtengine::SQR(i) + rtengine::SQR(j)) <= rtengine::SQR(3.0 * 1.5)) {
                kernel[i + 4][j + 4] = std::exp((rtengine::SQR(i) + rtengine::SQR(j)) / temp);
                sum += kernel[i + 4][j + 4];
            } else {
                kernel[i + 4][j + 4] = 0.f;
            }
        }
    }

    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 9; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

void compute7x7kernel(float sigma, float kernel[7][7]) {

    const double temp = -2.f * rtengine::SQR(sigma);
    float sum = 0.f;
    for (int i = -3; i <= 3; ++i) {
        for (int j = -3; j <= 3; ++j) {
            if((rtengine::SQR(i) + rtengine::SQR(j)) <= rtengine::SQR(3.0 * 1.15)) {
                kernel[i + 3][j + 3] = std::exp((rtengine::SQR(i) + rtengine::SQR(j)) / temp);
                sum += kernel[i + 3][j + 3];
            } else {
                kernel[i + 3][j + 3] = 0.f;
            }
        }
    }

    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

void compute5x5kernel(float sigma, float kernel[5][5]) {

    const double temp = -2.f * rtengine::SQR(sigma);
    float sum = 0.f;
    for (int i = -2; i <= 2; ++i) {
        for (int j = -2; j <= 2; ++j) {
            if((rtengine::SQR(i) + rtengine::SQR(j)) <= rtengine::SQR(3.0 * 0.84)) {
                kernel[i + 2][j + 2] = std::exp((rtengine::SQR(i) + rtengine::SQR(j)) / temp);
                sum += kernel[i + 2][j + 2];
            } else {
                kernel[i + 2][j + 2] = 0.f;
            }
        }
    }

    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

void compute3x3kernel(float sigma, float kernel[3][3]) {

    const double temp = -2.f * rtengine::SQR(sigma);
    float sum = 0.f;
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            kernel[i + 1][j + 1] = std::exp((rtengine::SQR(i) + rtengine::SQR(j)) / temp);
            sum += kernel[i + 1][j + 1];
        }
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            kernel[i][j] /= sum;
        }
    }
}

void gauss3x3div (float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[3][3])
{

    const float c11 = kernel[0][0];
    const float c10 = kernel[0][1];
    const float c00 = kernel[1][1];

    for (int i = 1; i < tileSize - 1; i++) {
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 1; j < tileSize - 1; j++) {
            const float val = c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) + 
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) + 
                              c00 * src[i][j];
            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

void gauss5x5div (float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[5][5])
{

    const float c21 = kernel[0][1];
    const float c20 = kernel[0][2];
    const float c11 = kernel[1][1];
    const float c10 = kernel[1][2];
    const float c00 = kernel[2][2];

    for (int i = 2; i < tileSize - 2; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 2; j < tileSize - 2; ++j) {
            const float val = c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

void gauss7x7div(float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[7][7])
{

    const float c31 = kernel[0][2];
    const float c30 = kernel[0][3];
    const float c22 = kernel[1][1];
    const float c21 = kernel[1][2];
    const float c20 = kernel[1][3];
    const float c11 = kernel[2][2];
    const float c10 = kernel[2][3];
    const float c00 = kernel[3][3];

    for (int i = 3; i < tileSize - 3; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 3; j < tileSize - 3; ++j) {
            const float val = c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

void gauss9x9div(float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[9][9])
{

    const float c42 = kernel[0][2];
    const float c41 = kernel[0][3];
    const float c40 = kernel[0][4];
    const float c33 = kernel[1][1];
    const float c32 = kernel[1][2];
    const float c31 = kernel[1][3];
    const float c30 = kernel[1][4];
    const float c22 = kernel[2][2];
    const float c21 = kernel[2][3];
    const float c20 = kernel[2][4];
    const float c11 = kernel[3][3];
    const float c10 = kernel[3][4];
    const float c00 = kernel[4][4];

    for (int i = 4; i < tileSize - 4; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 4; j < tileSize - 4; ++j) {
            const float val = c42 * ((src[i - 4][j - 2] + src[i - 4][j + 2]) + (src[i - 2][j - 4] + src[i - 2][j + 4]) + (src[i + 2][j - 4] + src[i + 2][j + 4]) + (src[i + 4][j - 2] + src[i + 4][j + 2])) +
                              c41 * ((src[i - 4][j - 1] + src[i - 4][j + 1]) + (src[i - 1][j - 4] + src[i - 1][j + 4]) + (src[i + 1][j - 4] + src[i + 1][j + 4]) + (src[i + 4][j - 1] + src[i + 4][j + 1])) +
                              c40 * (src[i - 4][j] + src[i][j - 4] + src[i][j + 4] + src[i + 4][j]) +
                              c33 * (src[i - 3][j - 3] + src[i - 3][j + 3] + src[i + 3][j - 3] + src[i + 3][j + 3]) +
                              c32 * ((src[i - 3][j - 2] + src[i - 3][j + 2]) + (src[i - 2][j - 3] + src[i - 2][j + 3]) + (src[i + 2][j - 3] + src[i + 2][j + 3]) + (src[i + 3][j - 2] + src[i + 3][j + 2])) +
                              c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

void gauss13x13div(float** RESTRICT src, float** RESTRICT dst, float** RESTRICT divBuffer, const int tileSize, const float kernel[13][13])
{
    const float c60 = kernel[0][6];
    const float c53 = kernel[1][3];
    const float c52 = kernel[1][4];
    const float c51 = kernel[1][5];
    const float c50 = kernel[1][6];
    const float c44 = kernel[2][2];
    const float c42 = kernel[2][4];
    const float c41 = kernel[2][5];
    const float c40 = kernel[2][6];
    const float c33 = kernel[3][3];
    const float c32 = kernel[3][4];
    const float c31 = kernel[3][5];
    const float c30 = kernel[3][6];
    const float c22 = kernel[4][4];
    const float c21 = kernel[4][5];
    const float c20 = kernel[4][6];
    const float c11 = kernel[5][5];
    const float c10 = kernel[5][6];
    const float c00 = kernel[6][6];

    for (int i = 6; i < tileSize - 6; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 6; j < tileSize - 6; ++j) {
            const float val = c60 * (src[i - 6][j] + src[i][j - 6] + src[i][j + 6] + src[i + 6][j]) +
                              c53 * ((src[i - 5][j - 3] + src[i - 5][j + 3]) + (src[i - 3][j - 5] + src[i - 3][j + 5]) + (src[i + 3][j - 5] + src[i + 3][j + 5]) + (src[i + 5][j - 3] + src[i + 5][j + 3])) +
                              c52 * ((src[i - 5][j - 2] + src[i - 5][j + 2]) + (src[i - 2][j - 5] + src[i - 2][j + 5]) + (src[i + 2][j - 5] + src[i + 2][j + 5]) + (src[i + 5][j - 2] + src[i + 5][j + 2])) +
                              c51 * ((src[i - 5][j - 1] + src[i - 5][j + 1]) + (src[i - 1][j - 5] + src[i - 1][j + 5]) + (src[i + 1][j - 5] + src[i + 1][j + 5]) + (src[i + 5][j - 1] + src[i + 5][j + 1])) +
                              c50 * ((src[i - 5][j] + src[i][j - 5] + src[i][j + 5] + src[i + 5][j]) + ((src[i - 4][j - 3] + src[i - 4][j + 3]) + (src[i - 3][j - 4] + src[i - 3][j + 4]) + (src[i + 3][j - 4] + src[i + 3][j + 4]) + (src[i + 4][j - 3] + src[i + 4][j + 3]))) +
                              c44 * (src[i - 4][j - 4] + src[i - 4][j + 4] + src[i + 4][j - 4] + src[i + 4][j + 4]) +
                              c42 * ((src[i - 4][j - 2] + src[i - 4][j + 2]) + (src[i - 2][j - 4] + src[i - 2][j + 4]) + (src[i + 2][j - 4] + src[i + 2][j + 4]) + (src[i + 4][j - 2] + src[i + 4][j + 2])) +
                              c41 * ((src[i - 4][j - 1] + src[i - 4][j + 1]) + (src[i - 1][j - 4] + src[i - 1][j + 4]) + (src[i + 1][j - 4] + src[i + 1][j + 4]) + (src[i + 4][j - 1] + src[i + 4][j + 1])) +
                              c40 * (src[i - 4][j] + src[i][j - 4] + src[i][j + 4] + src[i + 4][j]) +
                              c33 * (src[i - 3][j - 3] + src[i - 3][j + 3] + src[i + 3][j - 3] + src[i + 3][j + 3]) +
                              c32 * ((src[i - 3][j - 2] + src[i - 3][j + 2]) + (src[i - 2][j - 3] + src[i - 2][j + 3]) + (src[i + 2][j - 3] + src[i + 2][j + 3]) + (src[i + 3][j - 2] + src[i + 3][j + 2])) +
                              c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] = divBuffer[i][j] / std::max(val, 0.00001f);
        }
    }
}

void gauss3x3mult(float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[3][3])
{
    const float c11 = kernel[0][0];
    const float c10 = kernel[0][1];
    const float c00 = kernel[1][1];

    for (int i = 1; i < tileSize - 1; i++) {
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 1; j < tileSize - 1; j++) {
            const float val = c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) + 
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) + 
                              c00 * src[i][j];
            dst[i][j] *= val;
        }
    }

}

void gauss5x5mult (float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[5][5])
{

    const float c21 = kernel[0][1];
    const float c20 = kernel[0][2];
    const float c11 = kernel[1][1];
    const float c10 = kernel[1][2];
    const float c00 = kernel[2][2];

    for (int i = 2; i < tileSize - 2; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 2; j < tileSize - 2; ++j) {
            const float val = c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] *= val;
        }
    }
}

void gauss7x7mult(float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[7][7])
{

    const float c31 = kernel[0][2];
    const float c30 = kernel[0][3];
    const float c22 = kernel[1][1];
    const float c21 = kernel[1][2];
    const float c20 = kernel[1][3];
    const float c11 = kernel[2][2];
    const float c10 = kernel[2][3];
    const float c00 = kernel[3][3];

    for (int i = 3; i < tileSize - 3; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 3; j < tileSize - 3; ++j) {
            const float val = c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] *= val;
        }
    }
}

void gauss9x9mult(float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[9][9])
{

    const float c42 = kernel[0][2];
    const float c41 = kernel[0][3];
    const float c40 = kernel[0][4];
    const float c33 = kernel[1][1];
    const float c32 = kernel[1][2];
    const float c31 = kernel[1][3];
    const float c30 = kernel[1][4];
    const float c22 = kernel[2][2];
    const float c21 = kernel[2][3];
    const float c20 = kernel[2][4];
    const float c11 = kernel[3][3];
    const float c10 = kernel[3][4];
    const float c00 = kernel[4][4];

    for (int i = 4; i < tileSize - 4; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 4; j < tileSize - 4; ++j) {
            const float val = c42 * ((src[i - 4][j - 2] + src[i - 4][j + 2]) + (src[i - 2][j - 4] + src[i - 2][j + 4]) + (src[i + 2][j - 4] + src[i + 2][j + 4]) + (src[i + 4][j - 2] + src[i + 4][j + 2])) +
                              c41 * ((src[i - 4][j - 1] + src[i - 4][j + 1]) + (src[i - 1][j - 4] + src[i - 1][j + 4]) + (src[i + 1][j - 4] + src[i + 1][j + 4]) + (src[i + 4][j - 1] + src[i + 4][j + 1])) +
                              c40 * (src[i - 4][j] + src[i][j - 4] + src[i][j + 4] + src[i + 4][j]) +
                              c33 * (src[i - 3][j - 3] + src[i - 3][j + 3] + src[i + 3][j - 3] + src[i + 3][j + 3]) +
                              c32 * ((src[i - 3][j - 2] + src[i - 3][j + 2]) + (src[i - 2][j - 3] + src[i - 2][j + 3]) + (src[i + 2][j - 3] + src[i + 2][j + 3]) + (src[i + 3][j - 2] + src[i + 3][j + 2])) +
                              c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];
            dst[i][j] *= val;
        }
    }
}

void gauss13x13mult(float** RESTRICT src, float** RESTRICT dst, const int tileSize, const float kernel[13][13])
{

    const float c60 = kernel[0][6];
    const float c53 = kernel[1][3];
    const float c52 = kernel[1][4];
    const float c51 = kernel[1][5];
    const float c50 = kernel[1][6];
    const float c44 = kernel[2][2];
    const float c42 = kernel[2][4];
    const float c41 = kernel[2][5];
    const float c40 = kernel[2][6];
    const float c33 = kernel[3][3];
    const float c32 = kernel[3][4];
    const float c31 = kernel[3][5];
    const float c30 = kernel[3][6];
    const float c22 = kernel[4][4];
    const float c21 = kernel[4][5];
    const float c20 = kernel[4][6];
    const float c11 = kernel[5][5];
    const float c10 = kernel[5][6];
    const float c00 = kernel[6][6];

    for (int i = 6; i < tileSize - 6; ++i) {
        // I tried hand written SSE code but gcc vectorizes better
#if defined(__clang__)
        #pragma clang loop vectorize(assume_safety)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 6; j < tileSize - 6; ++j) {
            const float val = c60 * (src[i - 6][j] + src[i][j - 6] + src[i][j + 6] + src[i + 6][j]) +
                              c53 * ((src[i - 5][j - 3] + src[i - 5][j + 3]) + (src[i - 3][j - 5] + src[i - 3][j + 5]) + (src[i + 3][j - 5] + src[i + 3][j + 5]) + (src[i + 5][j - 3] + src[i + 5][j + 3])) +
                              c52 * ((src[i - 5][j - 2] + src[i - 5][j + 2]) + (src[i - 2][j - 5] + src[i - 2][j + 5]) + (src[i + 2][j - 5] + src[i + 2][j + 5]) + (src[i + 5][j - 2] + src[i + 5][j + 2])) +
                              c51 * ((src[i - 5][j - 1] + src[i - 5][j + 1]) + (src[i - 1][j - 5] + src[i - 1][j + 5]) + (src[i + 1][j - 5] + src[i + 1][j + 5]) + (src[i + 5][j - 1] + src[i + 5][j + 1])) +
                              c50 * ((src[i - 5][j] + src[i][j - 5] + src[i][j + 5] + src[i + 5][j]) + ((src[i - 4][j - 3] + src[i - 4][j + 3]) + (src[i - 3][j - 4] + src[i - 3][j + 4]) + (src[i + 3][j - 4] + src[i + 3][j + 4]) + (src[i + 4][j - 3] + src[i + 4][j + 3]))) +
                              c44 * (src[i - 4][j - 4] + src[i - 4][j + 4] + src[i + 4][j - 4] + src[i + 4][j + 4]) +
                              c42 * ((src[i - 4][j - 2] + src[i - 4][j + 2]) + (src[i - 2][j - 4] + src[i - 2][j + 4]) + (src[i + 2][j - 4] + src[i + 2][j + 4]) + (src[i + 4][j - 2] + src[i + 4][j + 2])) +
                              c41 * ((src[i - 4][j - 1] + src[i - 4][j + 1]) + (src[i - 1][j - 4] + src[i - 1][j + 4]) + (src[i + 1][j - 4] + src[i + 1][j + 4]) + (src[i + 4][j - 1] + src[i + 4][j + 1])) +
                              c40 * (src[i - 4][j] + src[i][j - 4] + src[i][j + 4] + src[i + 4][j]) +
                              c33 * (src[i - 3][j - 3] + src[i - 3][j + 3] + src[i + 3][j - 3] + src[i + 3][j + 3]) +
                              c32 * ((src[i - 3][j - 2] + src[i - 3][j + 2]) + (src[i - 2][j - 3] + src[i - 2][j + 3]) + (src[i + 2][j - 3] + src[i + 2][j + 3]) + (src[i + 3][j - 2] + src[i + 3][j + 2])) +
                              c31 * ((src[i - 3][j - 1] + src[i - 3][j + 1]) + (src[i - 1][j - 3] + src[i - 1][j + 3]) + (src[i + 1][j - 3] + src[i + 1][j + 3]) + (src[i + 3][j - 1] + src[i + 3][j + 1])) +
                              c30 * (src[i - 3][j] + src[i][j - 3] + src[i][j + 3] + src[i + 3][j]) +
                              c22 * (src[i - 2][j - 2] + src[i - 2][j + 2] + src[i + 2][j - 2] + src[i + 2][j + 2]) +
                              c21 * ((src[i - 2][j - 1] + src[i - 2][j + 1]) + (src[i - 1][j - 2] + src[i - 1][j + 2]) + (src[i + 1][j - 2] + src[i + 1][j + 2]) + (src[i + 2][j - 1] + src[i + 2][j + 1])) +
                              c20 * (src[i - 2][j] + src[i][j - 2] + src[i][j + 2] + src[i + 2][j]) +
                              c11 * (src[i - 1][j - 1] + src[i - 1][j + 1] + src[i + 1][j - 1] + src[i + 1][j + 1]) +
                              c10 * (src[i - 1][j] + src[i][j - 1] + src[i][j + 1] + src[i + 1][j]) +
                              c00 * src[i][j];

            dst[i][j] *= val;
        }
    }
}


void CaptureDeconvSharpening2 (float** luminance, const float* const * oldLuminance, const float * const * blend, int W, int H, float sigma, float sigmaCornerOffset, int iterations, bool checkIterStop, double startVal, double endVal)
{
BENCHFUN
printf("Sigma=%f \n", (double) sigma);
    const bool is9x9 = (sigma <= 1.5f && sigmaCornerOffset == 0.f);
    const bool is7x7 = (sigma <= 1.15f && sigmaCornerOffset == 0.f);
    const bool is5x5 = (sigma <= 0.84f && sigmaCornerOffset == 0.f);
    const bool is3x3 = (sigma < 0.6f && sigmaCornerOffset == 0.f);
    float kernel13[13][13];
    float kernel9[9][9];
    float kernel7[7][7];
    float kernel5[5][5];
    float kernel3[3][3];
    if (is3x3) {
        compute3x3kernel(sigma, kernel3);
    } else if (is5x5) {
        compute5x5kernel(sigma, kernel5);
    } else if (is7x7) {
        compute7x7kernel(sigma, kernel7);
    } else if (is9x9) {
        compute9x9kernel(sigma, kernel9);
    } else {
        compute13x13kernel(sigma, kernel13);
    }

    constexpr int tileSize = 32;
    const int border = (is3x3 || is5x5 || is7x7) ? iterations <= 30 ? 5 : 7 : 8;
    const int fullTileSize = tileSize + 2 * border;
    const float cornerRadius = std::min<float>(2.f, sigma + sigmaCornerOffset);
    const float cornerDistance = sqrt(rtengine::SQR(W * 0.5f) + rtengine::SQR(H * 0.5f));
    const float distanceFactor = (cornerRadius - sigma) / cornerDistance;

printf("fu=%i co=%f cod=%f df=%f\n", fullTileSize, cornerRadius, cornerDistance, distanceFactor);
   // double progress = startVal;
   // const double progressStep = (endVal - startVal) * rtengine::SQR(tileSize) / (W * H);

    constexpr float minBlend = 0.01f;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        int progresscounter = 0;
        array2D<float> tmpIThr(fullTileSize, fullTileSize);
        array2D<float> tmpThr(fullTileSize, fullTileSize);
        tmpThr.fill(1.f);
        array2D<float> lumThr(fullTileSize, fullTileSize);
        array2D<float> iterCheck(tileSize, tileSize);
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16) collapse(2)
#endif
        for (int i = border; i < H - border; i+= tileSize) {
            for(int j = border; j < W - border; j+= tileSize) {
                const bool endOfCol = (i + tileSize + border) >= H;
                const bool endOfRow = (j + tileSize + border) >= W;
                // fill tiles
                if (endOfRow || endOfCol) {
                    // special handling for small tiles at end of row or column
                    float maxVal = 0.f;
                    if (checkIterStop) {
                        for (int k = 0, ii = endOfCol ? H - fullTileSize + border : i; k < tileSize; ++k, ++ii) {
                            for (int l = 0, jj = endOfRow ? W - fullTileSize + border : j; l < tileSize; ++l, ++jj) {
                              //  if(blend[ii][jj] != 1.f) printf("z");
                                iterCheck[k][l] = oldLuminance[ii][jj] * blend[ii][jj] * 0.5f;
                                maxVal = std::max(maxVal, blend[ii][jj]);
                            }
                        }
                    } else {
                        for (int k = 0, ii = endOfCol ? H - fullTileSize + border : i; k < tileSize; ++k, ++ii) {
                            for (int l = 0, jj = endOfRow ? W - fullTileSize + border : j; l < tileSize; ++l, ++jj) {
                                maxVal = std::max(maxVal, blend[ii][jj]);
                            }
                        }
                    }
                    if (maxVal < minBlend) {
                        // no pixel of the tile has a blend factor >= minBlend => skip the tile
                        continue;
                    }
                    for (int k = 0, ii = endOfCol ? H - fullTileSize : i - border; k < fullTileSize; ++k, ++ii) {
                        for (int l = 0, jj = endOfRow ? W - fullTileSize : j - border; l < fullTileSize; ++l, ++jj) {
                            tmpIThr[k][l] = oldLuminance[ii][jj];
                            lumThr[k][l] = oldLuminance[ii][jj];
                        }
                    }
                } else {
                    float maxVal = 0.f;
                    if (checkIterStop) {
                        for (int ii = 0; ii < tileSize; ++ii) {
                            for (int jj = 0; jj < tileSize; ++jj) {
                                iterCheck[ii][jj] = oldLuminance[i + ii][j + jj] * blend[i + ii][j + jj] * 0.5f;
                                maxVal = std::max(maxVal, blend[i + ii][j + jj]);
                            }
                        }
                    } else {
                        for (int ii = 0; ii < tileSize; ++ii) {
                            for (int jj = 0; jj < tileSize; ++jj) {
                                maxVal = std::max(maxVal, blend[i + ii][j + jj]);
                            }
                        }
                    }
                    if (maxVal < minBlend) {
                        // no pixel of the tile has a blend factor >= minBlend => skip the tile
                        continue;
                    }
                    for (int ii = i; ii < i + fullTileSize; ++ii) {
                        for (int jj = j; jj < j + fullTileSize; ++jj) {
                            tmpIThr[ii - i][jj - j] = oldLuminance[ii - border][jj - border];
                            lumThr[ii - i][jj - j] = oldLuminance[ii - border][jj - border];
                        }
                    }
                }
                if (is3x3) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 3x3 gaussian blur and divide luminance by result of gaussian blur
                        gauss3x3div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel3);
                        gauss3x3mult(tmpThr, tmpIThr, fullTileSize, kernel3);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else if (is5x5) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 5x5 gaussian blur and divide luminance by result of gaussian blur
                        gauss5x5div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel5);
                        gauss5x5mult(tmpThr, tmpIThr, fullTileSize, kernel5);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else if (is7x7) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 5x5 gaussian blur and divide luminance by result of gaussian blur
                        gauss7x7div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel7);
                        gauss7x7mult(tmpThr, tmpIThr, fullTileSize, kernel7);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else if (is9x9) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 5x5 gaussian blur and divide luminance by result of gaussian blur
                        gauss9x9div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel9);
                        gauss9x9mult(tmpThr, tmpIThr, fullTileSize, kernel9);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else {
                    if (sigmaCornerOffset != 0.f) {
                        const float distance = sqrt(rtengine::SQR(i + tileSize / 2 - H / 2) + rtengine::SQR(j + tileSize / 2 - W / 2));
                        const float sigmaTile = static_cast<float>(sigma) + distanceFactor * distance;
                        if (sigmaTile >= 0.4f) {
                            if (sigmaTile > 1.5f) { // have to use 13x13 kernel
                                float lkernel13[13][13];
                                compute13x13kernel(static_cast<float>(sigma) + distanceFactor * distance, lkernel13);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 13x13 gaussian blur and divide luminance by result of gaussian blur
                                    gauss13x13div(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel13);
                                    gauss13x13mult(tmpThr, tmpIThr, fullTileSize, lkernel13);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            } else if (sigmaTile > 1.15f) { // have to use 9x9 kernel
                                float lkernel9[9][9];
                                compute9x9kernel(static_cast<float>(sigma) + distanceFactor * distance, lkernel9);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 9x9 gaussian blur and divide luminance by result of gaussian blur
                                    gauss9x9div(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel9);
                                    gauss9x9mult(tmpThr, tmpIThr, fullTileSize, lkernel9);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            } else if (sigmaTile > 0.84f) { // have to use 7x7 kernel
                                float lkernel7[7][7];
                                compute7x7kernel(static_cast<float>(sigma) + distanceFactor * distance, lkernel7);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 7x7 gaussian blur and divide luminance by result of gaussian blur
                                    gauss7x7div(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel7);
                                    gauss7x7mult(tmpThr, tmpIThr, fullTileSize, lkernel7);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            } else { // can use 5x5 kernel
                                float lkernel5[5][5];
                                compute5x5kernel(static_cast<float>(sigma) + distanceFactor * distance, lkernel5);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 7x7 gaussian blur and divide luminance by result of gaussian blur
                                    gauss5x5div(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel5);
                                    gauss5x5mult(tmpThr, tmpIThr, fullTileSize, lkernel5);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            }
                        }
                    } else {
                        for (int k = 0; k < iterations; ++k) {
                            // apply 13x13 gaussian blur and divide luminance by result of gaussian blur
                            gauss13x13div(tmpIThr, tmpThr, lumThr, fullTileSize, kernel13);
                            gauss13x13mult(tmpThr, tmpIThr, fullTileSize, kernel13);
                            if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                break;
                            }
                        }
                    }
                }
                if (endOfRow || endOfCol) {
                    // special handling for small tiles at end of row or column
                    for (int k = border, ii = endOfCol ? H - fullTileSize : i - border; k < fullTileSize - border; ++k) {
                        for (int l = border, jj = endOfRow ? W - fullTileSize : j - border; l < fullTileSize - border; ++l) {
                            luminance[ii + k][jj + l] = rtengine::intp(blend[ii + k][jj + l], tmpIThr[k][l], luminance[ii + k][jj + l]);
                        }
                    }
                } else {
                    for (int ii = border; ii < fullTileSize - border; ++ii) {
                        for (int jj = border; jj < fullTileSize - border; ++jj) {
                            luminance[i + ii - border][j + jj - border] = rtengine::intp(blend[i + ii - border][j + jj - border], tmpIThr[ii][jj], luminance[i + ii - border][j + jj - border]);
                        }
                    }
                }
                /*
                if (plistener) {
                    if (++progresscounter % 32 == 0) {
#ifdef _OPENMP
                        #pragma omp critical(csprogress)
#endif
                        {
                            progress += 32.0 * progressStep;
                            progress = rtengine::min(progress, endVal);
                            plistener->setProgress(progress);
                        }
                    }
                }
                */
            }
        }
    }
}



void ImProcFunctions::doSharpening(Imagefloat *rgb, int bfw, int bfh, int sk, float &sharpc, bool autoshar, float capradiu,  float deconvCo, float deconvLat, bool showMask)

{
    
    

    array2D<float> Y (bfw, bfh);

    get_luminance(rgb, Y, multiThread);

    float s_scale = std::sqrt(sk);
    float contrast = pow_F(sharpc / 100.f, 1.2f) * s_scale;
    JaggedArray<float> blend(bfw, bfh);
    buildBlendMask(Y, blend, bfw, bfh, contrast, autoshar);

    
    sharpc = 100.f * pow_F(contrast, 0.84f);
    printf("CONtrastSHAR=%f \n", (double) sharpc / s_scale);
    
  

    if (showMask) {
#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int i = 0; i < bfh; ++i) {
            for (int j = 0; j < bfw; ++j) {
                rgb->r(i, j)= rgb->g(i, j)= rgb->b(i, j) =  blend[i][j] * 65536.f;
                
              //  r[i][j] = g[i][j] = b[i][j] = blend[i][j] * 65536.f;
            }
        }

        return;
    }
   
    constexpr float xyz_rgb[3][3] = {          // XYZ from RGB
                                    { 0.412453, 0.357580, 0.180423 },
                                    { 0.212671, 0.715160, 0.072169 },
                                    { 0.019334, 0.119193, 0.950227 }
                                };
 
    array2D<float> clipMask(bfw, bfh);
  //  constexpr float clipLimit = 0.95f;
  //  constexpr float maxSigma = 2.f;
    
 //   float contrast = contra / 100.0;
    array2D<float> redVals (bfw, bfh);
    array2D<float> greenVals(bfw, bfh);
    array2D<float> blueVals(bfw, bfh);
         for (int i = 0; i < bfh; ++i) {
            for (int j = 0; j < bfw; ++j) {
                redVals[i][j] = rgb->r(i,j);
                greenVals[i][j] = rgb->g(i,j);
                blueVals[i][j] = rgb->b(i,j);
            }
        }

    array2D<float> L (bfw, bfh);
    array2D<float> YOld(bfw, bfh);
    array2D<float> YNew(bfw, bfh);
         for (int i = 0; i < bfh; ++i) {
            for (int j = 0; j < bfw; ++j) {
                L[i][j] = rgb->r(i,j);
                YOld[i][j] = rgb->g(i,j);
                YNew[i][j] = rgb->b(i,j);
            }
        }

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int i = 0; i < bfh; ++i) {
        Color::RGB2L(redVals[i], greenVals[i], blueVals[i], L[i], xyz_rgb, bfw);
        Color::RGB2Y(redVals[i], greenVals[i], blueVals[i], YOld[i], YNew[i], bfw);
    }
    buildBlendMask(L, clipMask, bfw, bfh, contrast, autoshar);//, clipMask);

    sharpc = contrast * 100.f;
    printf("SHAPC=%f \n", (double) sharpc);
    CaptureDeconvSharpening2(YNew, YOld, clipMask, bfw, bfh, capradiu, deconvCo, deconvLat, true, 0.2, 0.9);
 
          for (int i = 0; i < bfh; ++i) {
            for (int j = 0; j < bfw; ++j) {
                const float factor = YNew[i][j] / std::max(YOld[i][j], 0.00001f);
                rgb->r(i,j)= redVals[i][j] * factor;
                rgb->g(i,j)= greenVals[i][j] * factor;
                rgb->b(i,j)=  blueVals[i][j] * factor;
            }
        }
  
}

void ImProcFunctions::sharpening (LabImage* lab, const procparams::SharpeningParams &sharpenParam, bool showMask)
{

    if ((!sharpenParam.enabled) || sharpenParam.amount < 1 || lab->W < 8 || lab->H < 8) {
        return;
    }

    int W = lab->W, H = lab->H;

    // calculate contrast based blend factors to reduce sharpening in regions with low contrast
    JaggedArray<float> blend(W, H);
    float contrast = sharpenParam.contrast / 100.0;
    buildBlendMask(lab->L, blend, W, H, contrast);

    if(showMask) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                lab->L[i][j] = blend[i][j] * 32768.f;
            }
        }
        return;
    }

    JaggedArray<float> b2(W, H);

    if (sharpenParam.method == "rld") {
        deconvsharpening (lab->L, b2, blend, lab->W, lab->H, sharpenParam, scale);
        return;
    }
BENCHFUN

    // Rest is UNSHARP MASK
    float** b3 = nullptr;

    if (sharpenParam.edgesonly) {
        b3 = new float*[H];

        for (int i = 0; i < H; i++) {
            b3[i] = new float[W];
        }
    }

    JaggedArray<float> blur(W, H);

    if (sharpenParam.blurradius >= 0.25) {
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            gaussianBlur(lab->L, blur, W, H, sharpenParam.blurradius);
#ifdef _OPENMP
            #pragma omp for
#endif
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    blur[i][j] = intp(blend[i][j], lab->L[i][j], std::max(blur[i][j], 0.0f));
                }
            }
        }
    }


#ifdef _OPENMP
    #pragma omp parallel
#endif
    {

        if (!sharpenParam.edgesonly) {
            gaussianBlur (lab->L, b2, W, H, sharpenParam.radius / scale);
        } else {
            bilateral<float, float> (lab->L, (float**)b3, b2, W, H, sharpenParam.edges_radius / scale, sharpenParam.edges_tolerance, multiThread);
            gaussianBlur (b3, b2, W, H, sharpenParam.radius / scale);
        }
    }
    float** base = lab->L;

    if (sharpenParam.edgesonly) {
        base = b3;
    }

    if (!sharpenParam.halocontrol) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = 0; i < H; i++)
            for (int j = 0; j < W; j++) {
                constexpr float upperBound = 2000.f;  // WARNING: Duplicated value, it's baaaaaad !
                float diff = base[i][j] - b2[i][j];
                float delta = sharpenParam.threshold.multiply<float, float, float>(
                                  min(fabsf(diff), upperBound),                   // X axis value = absolute value of the difference, truncated to the max value of this field
                                  sharpenParam.amount * diff * 0.01f        // Y axis max value
                              );
                lab->L[i][j] = intp(blend[i][j], lab->L[i][j] + delta, lab->L[i][j]);
            }
    } else {
        if (!sharpenParam.edgesonly) {
            // make a deep copy of lab->L
            JaggedArray<float> labCopy(W, H);

#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for( int i = 0; i < H; i++ )
                for( int j = 0; j < W; j++ ) {
                    labCopy[i][j] = lab->L[i][j];
                }

            sharpenHaloCtrl (lab->L, b2, labCopy, blend, W, H, sharpenParam);
        } else {
            sharpenHaloCtrl (lab->L, b2, base, blend, W, H, sharpenParam);
        }

    }

    if (sharpenParam.edgesonly) {
        for (int i = 0; i < H; i++) {
            delete [] b3[i];
        }

        delete [] b3;
    }

    if (sharpenParam.blurradius >= 0.25) {
#ifdef _OPENMP
    #pragma omp parallel for
#endif
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                lab->L[i][j] = intp(blend[i][j], lab->L[i][j], max(blur[i][j], 0.0f));
            }
        }
    }

}

// To the extent possible under law, Manuel Llorens <manuelllorens@gmail.com>
// has waived all copyright and related or neighboring rights to this work.
// This code is licensed under CC0 v1.0, see license information at
// http://creativecommons.org/publicdomain/zero/1.0/

//! MicroContrast is a sharpening method developed by Manuel Llorens and documented here: http://www.rawness.es/sharpening/?lang=en
//! <BR>The purpose is maximize clarity of the image without creating halo's.
//! <BR>Addition from JD : pyramid  + pondered contrast with matrix 5x5
//! <BR>2017 Ingo Weyrich : reduced processing time
//! \param luminance : Luminance channel of image
void ImProcFunctions::MLmicrocontrast(float** luminance, int W, int H)
{
    if (!params->sharpenMicro.enabled || params->sharpenMicro.contrast == 100 || params->sharpenMicro.amount < 1.0) {
        return;
    }
BENCHFUN
    const int k = params->sharpenMicro.matrix ? 1 : 2;
    // k=2 matrix 5x5  k=1 matrix 3x3
    const int width = W, height = H;
    const int unif = params->sharpenMicro.uniformity;
    const float amount = (k == 1 ? 2.7 : 1.) * params->sharpenMicro.amount / 1500.0; //amount 2000.0 quasi no artifacts ==> 1500 = maximum, after artifacts, 25/9 if 3x3

    if (settings->verbose) {
        printf ("Micro-contrast amount %f\n", static_cast<double>(amount));
        printf ("Micro-contrast uniformity %i\n", unif);
    }

    //modulation uniformity in function of luminance
    const float L98[11] = {0.001f, 0.0015f, 0.002f, 0.004f, 0.006f, 0.008f, 0.01f, 0.03f, 0.05f, 0.1f, 0.1f};
    const float L95[11] = {0.0012f, 0.002f, 0.005f, 0.01f, 0.02f, 0.05f, 0.1f, 0.12f, 0.15f, 0.2f, 0.25f};
    const float L92[11] = {0.01f, 0.015f, 0.02f, 0.06f, 0.10f, 0.13f, 0.17f, 0.25f, 0.3f, 0.32f, 0.35f};
    const float L90[11] = {0.015f, 0.02f, 0.04f, 0.08f, 0.12f, 0.15f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f};
    const float L87[11] = {0.025f, 0.03f, 0.05f, 0.1f, 0.15f, 0.25f, 0.3f, 0.4f, 0.5f, 0.63f, 0.75f};
    const float L83[11] = {0.055f, 0.08f, 0.1f, 0.15f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.75f, 0.85f};
    const float L80[11] = {0.15f, 0.2f, 0.25f, 0.3f, 0.35f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f};
    const float L75[11] = {0.22f, 0.25f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.85f, 0.9f, 0.95f};
    const float L70[11] = {0.35f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.97f, 1.0f, 1.0f, 1.0f, 1.0f};
    const float L63[11] = {0.55f, 0.6f, 0.7f, 0.8f, 0.85f, 0.9f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    const float L58[11] = {0.75f, 0.77f, 0.8f, 0.9f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    //default 5
    //modulation contrast
    const float Cont0[11] = {0.05f, 0.1f, 0.2f, 0.25f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f};
    const float Cont1[11] = {0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 0.95f, 1.0f};
    const float Cont2[11] = {0.2f, 0.40f, 0.6f, 0.7f, 0.8f, 0.85f, 0.90f, 0.95f, 1.0f, 1.05f, 1.10f};
    const float Cont3[11] = {0.5f, 0.6f, 0.7f, 0.8f, 0.85f, 0.9f, 1.0f, 1.0f, 1.05f, 1.10f, 1.20f};
    const float Cont4[11] = {0.8f, 0.85f, 0.9f, 0.95f, 1.0f, 1.05f, 1.10f, 1.150f, 1.2f, 1.25f, 1.40f};
    const float Cont5[11] = {1.0f, 1.1f, 1.2f, 1.25f, 1.3f, 1.4f, 1.45f, 1.50f, 1.6f, 1.65f, 1.80f};

    const float sqrt2 = sqrt(2.0);
    const float sqrt1d25 = sqrt(1.25);
    float *LM = new float[width * height]; //allocation for Luminance

    // calculate contrast based blend factors to reduce sharpening in regions with low contrast
    JaggedArray<float> blend(W, H);
    float contrastThreshold = params->sharpenMicro.contrast / 100.0;
    buildBlendMask(luminance, blend, W, H, contrastThreshold);

#ifdef _OPENMP
    #pragma omp parallel
#endif
{

#ifdef _OPENMP
    #pragma omp for schedule(dynamic,16)
#endif

    for(int j = 0; j < height; j++)
        for(int i = 0, offset = j * width + i; i < width; i++, offset++) {
            LM[offset] = luminance[j][i] / 327.68f; // adjust to [0;100] and to RT variables
        }

#ifdef _OPENMP
    #pragma omp for schedule(dynamic,16)
#endif

    for(int j = k; j < height - k; j++)
        for(int i = k, offset = j * width + i; i < width - k; i++, offset++) {
            float v = LM[offset];

            float contrast;
            if (k == 1) {
                contrast = sqrtf(SQR(LM[offset + 1] - LM[offset - 1]) + SQR(LM[offset + width] - LM[offset - width])) * 0.125f;    //for 3x3
            } else /* if (k==2) */ contrast = sqrtf(SQR(LM[offset + 1] - LM[offset - 1]) + SQR(LM[offset + width] - LM[offset - width])
                                                       + SQR(LM[offset + 2] - LM[offset - 2]) + SQR(LM[offset + 2 * width] - LM[offset - 2 * width])) * 0.0625f; //for 5x5

            contrast = std::min(contrast, 1.f);

            //matrix 5x5
            float temp = v + 4.f *( v * (amount + sqrt2 * amount)); //begin 3x3
            float temp1 = sqrt2 * amount *(LM[offset - width - 1] + LM[offset - width + 1] + LM[offset + width - 1] + LM[offset + width + 1]);
            temp1 += amount * (LM[offset - width] + LM[offset - 1] + LM[offset + 1] + LM[offset + width]);

            temp -= temp1;

            // add JD continue 5x5
            if (k == 2) {
                float temp2 = -(LM[offset + 2 * width] + LM[offset - 2 * width] + LM[offset - 2] + LM[offset + 2]);

                temp2 -= sqrt1d25 * (LM[offset + 2 * width - 1] + LM[offset + 2 * width + 1] + LM[offset +  width + 2] + LM[offset +  width - 2] + 
                                     LM[offset - 2 * width - 1] + LM[offset - 2 * width + 1] + LM[offset -  width + 2] + LM[offset -  width - 2]);

                temp2 -= sqrt2 * (LM[offset + 2 * width - 2] + LM[offset + 2 * width + 2] + LM[offset - 2 * width - 2] + LM[offset - 2 * width + 2]);
                temp2 += 18.601126159f * v ; // 18.601126159 = 4 + 4 * sqrt(2) + 8 * sqrt(1.25)
                temp2 *= 2.f * amount;
                temp += temp2;
            }

            temp = std::max(temp, 0.f);

            for(int row = j - k; row <= j + k; ++row) {
                for(int offset2 = row * width + i - k; offset2 <= row * width + i + k; ++offset2) {
                    if((LM[offset2] - temp) * (v - LM[offset2]) > 0.f) {
                        temp = intp(0.75f, temp, LM[offset2]);
                        goto breakout;
                    }
                }
            }
            breakout:

            if (LM[offset] > 95.0f || LM[offset] < 5.0f) {
                contrast *= Cont0[unif];    //+ JD : luminance  pyramid to adjust contrast by evaluation of LM[offset]
            } else if (LM[offset] > 90.0f || LM[offset] < 10.0f) {
                contrast *= Cont1[unif];
            } else if (LM[offset] > 80.0f || LM[offset] < 20.0f) {
                contrast *= Cont2[unif];
            } else if (LM[offset] > 70.0f || LM[offset] < 30.0f) {
                contrast *= Cont3[unif];
            } else if (LM[offset] > 60.0f || LM[offset] < 40.0f) {
                contrast *= Cont4[unif];
            } else {
                contrast *= Cont5[unif];    //(2.0f/k)*Cont5[unif];
            }

            contrast = std::min(contrast, 1.f);

            float tempL = intp(contrast, LM[offset], temp);
            // JD: modulation of microcontrast in function of original Luminance and modulation of luminance
            if (tempL > LM[offset]) {
                float temp2 = tempL / LM[offset]; //for highlights
                temp2 = std::min(temp2, 1.7f); //limit action
                temp2 -= 1.f;
                if (LM[offset] > 98.0f) {
                    temp = 0.f;
                } else if (LM[offset] > 95.0f) {
                    temp = L95[unif];
                } else if (LM[offset] > 92.0f) {
                    temp = L92[unif];
                } else if (LM[offset] > 90.0f) {
                    temp = L90[unif];
                } else if (LM[offset] > 87.0f) {
                    temp = L87[unif];
                } else if (LM[offset] > 83.0f) {
                    temp = L83[unif];
                } else if (LM[offset] > 80.0f) {
                    temp = L80[unif];
                } else if (LM[offset] > 75.0f) {
                    temp = L75[unif];
                } else if (LM[offset] > 70.0f) {
                    temp = L70[unif];
                } else if (LM[offset] > 63.0f) {
                    temp = L63[unif];
                } else if (LM[offset] > 58.0f) {
                    temp = L58[unif];
                } else if (LM[offset] > 42.0f) {
                    temp = L58[unif];
                } else if (LM[offset] > 37.0f) {
                    temp = L63[unif];
                } else if (LM[offset] > 30.0f) {
                    temp = L70[unif];
                } else if (LM[offset] > 25.0f) {
                    temp = L75[unif];
                } else if (LM[offset] > 20.0f) {
                    temp = L80[unif];
                } else if (LM[offset] > 17.0f) {
                    temp = L83[unif];
                } else if (LM[offset] > 13.0f) {
                    temp = L87[unif];
                } else if (LM[offset] > 10.0f) {
                    temp = L90[unif];
                } else if (LM[offset] > 5.0f) {
                    temp = L95[unif];
                } else {
                    temp = 0.f;
                }
                luminance[j][i] = intp(blend[j][i], luminance[j][i] * (temp * temp2 + 1.f), luminance[j][i]);
            } else {

                float temp4 = LM[offset] / tempL; //

                if (temp4 > 1.0f) {
                    temp4 = std::min(temp4, 1.7f); //limit action
                    temp4 -= 1.f;
                    if (LM[offset] < 2.0f) {
                        temp = L98[unif];
                    } else if (LM[offset] < 5.0f) {
                        temp = L95[unif];
                    } else if (LM[offset] < 8.0f) {
                        temp = L92[unif];
                    } else if (LM[offset] < 10.0f) {
                        temp = L90[unif];
                    } else if (LM[offset] < 13.0f) {
                        temp = L87[unif];
                    } else if (LM[offset] < 17.0f) {
                        temp = L83[unif];
                    } else if (LM[offset] < 20.0f) {
                        temp = L80[unif];
                    } else if (LM[offset] < 25.0f) {
                        temp = L75[unif];
                    } else if (LM[offset] < 30.0f) {
                        temp = L70[unif];
                    } else if (LM[offset] < 37.0f) {
                        temp = L63[unif];
                    } else if (LM[offset] < 42.0f) {
                        temp = L58[unif];
                    } else if (LM[offset] < 58.0f) {
                        temp = L58[unif];
                    } else if (LM[offset] < 63.0f) {
                        temp = L63[unif];
                    } else if (LM[offset] < 70.0f) {
                        temp = L70[unif];
                    } else if (LM[offset] < 75.0f) {
                        temp = L75[unif];
                    } else if (LM[offset] < 80.0f) {
                        temp = L80[unif];
                    } else if (LM[offset] < 83.0f) {
                        temp = L83[unif];
                    } else if (LM[offset] < 87.0f) {
                        temp = L87[unif];
                    } else if (LM[offset] < 90.0f) {
                        temp = L90[unif];
                    } else if (LM[offset] < 95.0f) {
                        temp = L95[unif];
                    } else {
                        temp = 0.f;
                    }
                    luminance[j][i] = intp(blend[j][i], luminance[j][i] / (temp * temp4 + 1.f), luminance[j][i]);
                }
            }
        }
}
    delete [] LM;
}

void ImProcFunctions::MLmicrocontrast(LabImage* lab)
{
    MLmicrocontrast(lab->L, lab->W, lab->H);
}

void ImProcFunctions::MLmicrocontrastcam(CieImage* ncie)
{
    MLmicrocontrast(ncie->sh_p, ncie->W, ncie->H);
}

void ImProcFunctions::sharpeningcam (CieImage* ncie, float** b2, bool showMask)
{
    if ((!params->sharpening.enabled) || params->sharpening.amount < 1 || ncie->W < 8 || ncie->H < 8) {
        return;
    }

    int W = ncie->W, H = ncie->H;

    // calculate contrast based blend factors to reduce sharpening in regions with low contrast
    JaggedArray<float> blend(W, H);
    float contrast = params->sharpening.contrast / 100.0;
    buildBlendMask(ncie->sh_p, blend, W, H, contrast);
    if(showMask) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                ncie->sh_p[i][j] = blend[i][j] * 32768.f;
            }
        }
        return;
    }

    if (params->sharpening.method == "rld") {
        deconvsharpening (ncie->sh_p, b2, blend, ncie->W, ncie->H, params->sharpening, scale);
        return;
    }

    // Rest is UNSHARP MASK

    float** b3 = nullptr;

    if (params->sharpening.edgesonly) {
        b3 = new float*[H];

        for (int i = 0; i < H; i++) {
            b3[i] = new float[W];
        }
    }

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {

        if (!params->sharpening.edgesonly) {
            gaussianBlur (ncie->sh_p, b2, W, H, params->sharpening.radius / scale);
        } else {
            bilateral<float, float> (ncie->sh_p, (float**)b3, b2, W, H, params->sharpening.edges_radius / scale, params->sharpening.edges_tolerance, multiThread);
            gaussianBlur (b3, b2, W, H, params->sharpening.radius / scale);
        }
    }

    float** base = ncie->sh_p;

    if (params->sharpening.edgesonly) {
        base = b3;
    }

    if (!params->sharpening.halocontrol) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = 0; i < H; i++)
            for (int j = 0; j < W; j++) {
                constexpr float upperBound = 2000.f;  // WARNING: Duplicated value, it's baaaaaad !
                float diff = base[i][j] - b2[i][j];
                float delta = params->sharpening.threshold.multiply<float, float, float>(
                                  min(fabsf(diff), upperBound),                   // X axis value = absolute value of the difference, truncated to the max value of this field
                                  params->sharpening.amount * diff * 0.01f      // Y axis max value
                              );

                if(ncie->J_p[i][j] > 8.0f && ncie->J_p[i][j] < 92.0f) {
                    ncie->sh_p[i][j] = intp(blend[i][j], ncie->sh_p[i][j] + delta, ncie->sh_p[i][j]);
                }
            }
    } else {
        float** ncieCopy = nullptr;

        if (!params->sharpening.edgesonly) {
            // make deep copy of ncie->sh_p
            ncieCopy = new float*[H];

            for( int i = 0; i < H; i++ ) {
                ncieCopy[i] = new float[W];
            }

#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for( int i = 0; i < H; i++ )
                for( int j = 0; j < W; j++ ) {
                    ncieCopy[i][j] = ncie->sh_p[i][j];
                }

            base = ncieCopy;
        }

        sharpenHaloCtrl (ncie->sh_p, b2, base, blend, W, H, params->sharpening);

        if(ncieCopy) {
            for( int i = 0; i < H; i++ ) {
                delete[] ncieCopy[i];
            }

            delete[] ncieCopy;
        }
    }

    if (params->sharpening.edgesonly) {
        for (int i = 0; i < H; i++) {
            delete [] b3[i];
        }

        delete [] b3;
    }
}

}
