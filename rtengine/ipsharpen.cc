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
#include "boxblur.h"

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
            tmpI[i][j] = luminance[i][j] = max(luminance[i][j], 0.f);
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
            //tmpI[i][j] = max(luminance[i][j], 0.f)
            tmpI[i][j] = luminance[i][j] = max(luminance[i][j], 0.f);
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

void ImProcFunctions::CaptureDeconvSharpening_SE (float** luminance, const float* const * oldLuminance, const float * const * blend, int bfw, int bfh, struct localpass &locp, float sigma, float sigmaCornerOffset, int iterations, bool checkIterStop, double startVal, double endVal)
{
 // Copyright (c) 2019 Ingo Weyrich (heckflosse67@gmx.de)
 // adaptation november 2024 - Jacques Desmis - for Selective Editing 
BENCHFUN
   
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
        rtengine::compute3x3kernel2(sigma, kernel3);
    } else if (is5x5) {
        rtengine::compute5x5kernel2(sigma, kernel5);
    } else if (is7x7) {
        rtengine::compute7x7kernel2(sigma, kernel7);
    } else if (is9x9) {
        rtengine::compute9x9kernel2(sigma, kernel9);
    } else {
        rtengine::compute13x13kernel2(sigma, kernel13);
    }
    
    const float nocoboot = locp.sizenocoboot;
    constexpr int tileSize = 32;
    const int border = (is3x3 || is5x5 || is7x7) ? iterations <= 30 ? 5 : 7 : 8;
    const int fullTileSize = tileSize + 2 * border;
    const float cornerRadius = std::min<float>(2.f, sigma + sigmaCornerOffset);
    const float cornerDistance = sqrt(rtengine::SQR(bfw * nocoboot) + rtengine::SQR(bfh * nocoboot));
    const float distanceFactor = (cornerRadius - sigma) / cornerDistance;

    if (settings->verbose) {
        printf("sigma=%f Tilesize=%i cornerrad=%f cornerDist=%f distanceFactor=%f\n", sigma, fullTileSize, cornerRadius, cornerDistance, distanceFactor);
    }
    //I don't use startval and endval to show progress... but ...
    constexpr float minBlend = 0.01f;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        array2D<float> tmpIThr(fullTileSize, fullTileSize);
        array2D<float> tmpThr(fullTileSize, fullTileSize);
        tmpThr.fill(1.f);
        array2D<float> lumThr(fullTileSize, fullTileSize);
        array2D<float> iterCheck(tileSize, tileSize);
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16) collapse(2)
#endif
        for (int i = border; i < bfh - border; i+= tileSize) {
            for(int j = border; j < bfw - border; j+= tileSize) {
                const bool endOfCol = (i + tileSize + border) >= bfh;
                const bool endOfRow = (j + tileSize + border) >= bfw;
                // fill tiles
                if (endOfRow || endOfCol) {
                    // special handling for small tiles at end of row or column
                    float maxVal = 0.f;
                    if (checkIterStop) {
                        for (int k = 0, ii = endOfCol ? bfh - fullTileSize + border : i; k < tileSize; ++k, ++ii) {
                            for (int l = 0, jj = endOfRow ? bfw - fullTileSize + border : j; l < tileSize; ++l, ++jj) {
                                iterCheck[k][l] = oldLuminance[ii][jj] * blend[ii][jj] * 0.5f;
                                maxVal = std::max(maxVal, blend[ii][jj]);
                            }
                        }
                    } else {
                        for (int k = 0, ii = endOfCol ? bfh - fullTileSize + border : i; k < tileSize; ++k, ++ii) {
                            for (int l = 0, jj = endOfRow ? bfw - fullTileSize + border : j; l < tileSize; ++l, ++jj) {
                                maxVal = std::max(maxVal, blend[ii][jj]);
                            }
                        }
                    }
                    if (maxVal < minBlend) {
                        // no pixel of the tile has a blend factor >= minBlend => skip the tile
                        continue;
                    }
                    for (int k = 0, ii = endOfCol ? bfh - fullTileSize : i - border; k < fullTileSize; ++k, ++ii) {
                        for (int l = 0, jj = endOfRow ? bfw - fullTileSize : j - border; l < fullTileSize; ++l, ++jj) {
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
                        rtengine::gauss3x3div2(tmpIThr, tmpThr, lumThr, fullTileSize, kernel3);
                        rtengine::gauss3x3mult2(tmpThr, tmpIThr, fullTileSize, kernel3);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else if (is5x5) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 5x5 gaussian blur and divide luminance by result of gaussian blur
                        rtengine::gauss5x5div2(tmpIThr, tmpThr, lumThr, fullTileSize, kernel5);
                        rtengine::gauss5x5mult2(tmpThr, tmpIThr, fullTileSize, kernel5);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else if (is7x7) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 5x5 gaussian blur and divide luminance by result of gaussian blur
                        rtengine::gauss7x7div2(tmpIThr, tmpThr, lumThr, fullTileSize, kernel7);
                        rtengine::gauss7x7mult2(tmpThr, tmpIThr, fullTileSize, kernel7);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else if (is9x9) {
                    for (int k = 0; k < iterations; ++k) {
                        // apply 5x5 gaussian blur and divide luminance by result of gaussian blur
                        rtengine::gauss9x9div2(tmpIThr, tmpThr, lumThr, fullTileSize, kernel9);
                        rtengine::gauss9x9mult2(tmpThr, tmpIThr, fullTileSize, kernel9);
                        if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                            break;
                        }
                    }
                } else {
                    if (sigmaCornerOffset != 0.f) {
                        const float distance = sqrt(rtengine::SQR(i + tileSize / 2 - bfh / 2) + rtengine::SQR(j + tileSize / 2 - bfw / 2));
                        const float sigmaTile = static_cast<float>(sigma) + distanceFactor * distance;
                        if (sigmaTile >= 0.4f) {
                            if (sigmaTile > 1.5f) { // have to use 13x13 kernel
                                float lkernel13[13][13];
                                rtengine::compute13x13kernel2(static_cast<float>(sigma) + distanceFactor * distance, lkernel13);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 13x13 gaussian blur and divide luminance by result of gaussian blur
                                    rtengine::gauss13x13div2(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel13);
                                    rtengine::gauss13x13mult2(tmpThr, tmpIThr, fullTileSize, lkernel13);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            } else if (sigmaTile > 1.15f) { // have to use 9x9 kernel
                                float lkernel9[9][9];
                                rtengine::compute9x9kernel2(static_cast<float>(sigma) + distanceFactor * distance, lkernel9);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 9x9 gaussian blur and divide luminance by result of gaussian blur
                                    rtengine::gauss9x9div2(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel9);
                                    rtengine::gauss9x9mult2(tmpThr, tmpIThr, fullTileSize, lkernel9);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            } else if (sigmaTile > 0.84f) { // have to use 7x7 kernel
                                float lkernel7[7][7];
                                rtengine::compute7x7kernel2(static_cast<float>(sigma) + distanceFactor * distance, lkernel7);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 7x7 gaussian blur and divide luminance by result of gaussian blur
                                    rtengine::gauss7x7div2(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel7);
                                    rtengine::gauss7x7mult2(tmpThr, tmpIThr, fullTileSize, lkernel7);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            } else { // can use 5x5 kernel
                                float lkernel5[5][5];
                                rtengine::compute5x5kernel2(static_cast<float>(sigma) + distanceFactor * distance, lkernel5);
                                for (int k = 0; k < iterations; ++k) {
                                    // apply 7x7 gaussian blur and divide luminance by result of gaussian blur
                                    rtengine::gauss5x5div2(tmpIThr, tmpThr, lumThr, fullTileSize, lkernel5);
                                    rtengine::gauss5x5mult2(tmpThr, tmpIThr, fullTileSize, lkernel5);
                                    if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                        break;
                                    }
                                }
                            }
                        }
                    } else {
                        for (int k = 0; k < iterations; ++k) {
                            // apply 13x13 gaussian blur and divide luminance by result of gaussian blur
                            rtengine::gauss13x13div2(tmpIThr, tmpThr, lumThr, fullTileSize, kernel13);
                            rtengine::gauss13x13mult2(tmpThr, tmpIThr, fullTileSize, kernel13);
                            if (checkIterStop && k < iterations - 1 && checkForStop(tmpIThr, iterCheck, fullTileSize, border)) {
                                break;
                            }
                        }
                    }
                }
                if (endOfRow || endOfCol) {
                    // special handling for small tiles at end of row or column
                    for (int k = border, ii = endOfCol ? bfh - fullTileSize : i - border; k < fullTileSize - border; ++k) {
                        for (int l = border, jj = endOfRow ? bfw - fullTileSize : j - border; l < fullTileSize - border; ++l) {
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
            }
        }
    }
}

float igammalog2(float x, float p, float s, float g2, float g4) // same function as in iplocallab.cc
{
    return x <= g2 ? x / s : pow_F((x + g4) / (1.f + g4), p);//continuous
}


float gammalog2(float x, float p, float s, float g3, float g4) // same function as in iplocallab.cc
{
    return x <= g3 ? x * s : (1.f + g4) * xexpf(xlogf(x) / p) - g4;
}


void ImProcFunctions::doCapture_Sharpening_SE(Imagefloat *rgb, int bfw, int bfh, struct localpass &locp, int sk, float &sharpc, bool autoshar, float capradiu,  float deconvCo, float deconvLat, bool itcheck, bool showMask, float deconvgam)
//Jacques Desmis - 
{
    
    TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);

    const float wip[3][3] = {
        {(float) wprof[0][0], (float) wprof[0][1], (float) wprof[0][2]},
        {(float) wprof[1][0], (float) wprof[1][1], (float) wprof[1][2]},
        {(float) wprof[2][0], (float) wprof[2][1], (float) wprof[2][2]}
    };
        //tempory variables for gamma
        array2D<float> redgam (bfw, bfh);
        array2D<float> greengam(bfw, bfh);
        array2D<float> bluegam(bfw, bfh);
        
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
        for (int i = 0; i < bfh; ++i) {//save temp gamma
            for (int j = 0; j < bfw; ++j) {
                redgam[i][j] = rgb->r(i,j);
                greengam[i][j] = rgb->g(i,j);
                bluegam[i][j] = rgb->b(i,j);
            }
        }
        //calculate gamma as it was with Lab calculation but in RGB linear mode to have the same 'behavior' as other gamma in Selective Editing
        //This is actually not of major importance, the main thing is to use the inverse function afterwards
        float gamma1 = deconvgam;
        rtengine::GammaValues g_a; //gamma parameters
        double pwr1 = 1.0 / (double) gamma1;//default 3.0 - gamma Lab
        double ts1 = 9.03296;//always the same 'slope' in the extreme shadows - slope Lab , I can choose also 5 or 13!
        rtengine::Color::calcGamma(pwr1, ts1, g_a); // call to calcGamma with selected gamma and slope
        
        if (gamma1 != 1.f) {//calculate new values with gamma for R, G, B of course with 65535 instead of 32768
#ifdef _OPENMP
            #   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int i = 0; i < bfh; ++i) {
                for (int j = 0; j < bfw; ++j) {
                    redgam[i][j] = 65535.f * igammalog2(redgam[i][j] / 65535.f, gamma1, ts1, g_a[2], g_a[4]);
                    greengam[i][j] = 65535.f * igammalog2(greengam[i][j] / 65535.f, gamma1, ts1, g_a[2], g_a[4]);
                    bluegam[i][j] = 65535.f * igammalog2(bluegam[i][j] / 65535.f, gamma1, ts1, g_a[2], g_a[4]);
                }
            }
        }
        
        
    float s_scale = std::sqrt(sk);
   
    if (showMask) {
        array2D<float> clipMask(bfw, bfh);

        array2D<float> redVals (bfw, bfh);
        array2D<float> greenVals(bfw, bfh);
        array2D<float> blueVals(bfw, bfh);

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif

        for (int i = 0; i < bfh; ++i) {
            for (int j = 0; j < bfw; ++j) {
                redVals[i][j] = redgam[i][j];
                greenVals[i][j] = greengam[i][j];
                blueVals[i][j] = bluegam[i][j];
            }
        }
        
        array2D<float> Y (bfw, bfh);
        float contrastsh = pow_F(sharpc / 100.f, 1.f) * s_scale;
        
         for (int i = 0; i < bfh; ++i) {
            Color::RGB2L(redVals[i], greenVals[i], blueVals[i], Y[i], wip, bfw);//mask values with gamma
        }
        float reducautocontrast = 1.f;//to take noise into account

        buildBlendMask2(Y, clipMask, bfw, bfh, contrastsh, 1.f, autoshar, 2.f / s_scale, 1.f, reducautocontrast);//build with gamma if need

        sharpc = 100.f * pow_F(contrastsh, 1.f)/ s_scale;
       
        
#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int i = 0; i < bfh; ++i) {
            for (int j = 0; j < bfw; ++j) {
                rgb->r(i, j) = rgb->g(i, j) = rgb->b(i, j) =  clipMask[i][j] * 65536.f; //same values R G B to black and white image, just for mask             
            }
        }
        if (settings->verbose) {
            int autos = 0;
            if(autoshar) {autos = 1;}
            printf("Contrast threshold Selective Editing Capture Show mask=%f auto=%i\n", (double) sharpc, autos);
        }

    } else {//general case
        array2D<float> clipMask2(bfw, bfh); 
        array2D<float> redVal (bfw, bfh);
        array2D<float> greenVal(bfw, bfh);
        array2D<float> blueVal(bfw, bfh);

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif   
         for (int i = 0; i < bfh; ++i) {//values with gamma
            for (int j = 0; j < bfw; ++j) {
                redVal[i][j] = redgam[i][j];
                greenVal[i][j] = greengam[i][j];
                blueVal[i][j] = bluegam[i][j];
            }
        }

        array2D<float> L (bfw, bfh);
        array2D<float> YOld(bfw, bfh);
        array2D<float> YNew(bfw, bfh);

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif        
         for (int i = 0; i < bfh; ++i) {
            for (int j = 0; j < bfw; ++j) {
                L[i][j] = redgam[i][j];
                YOld[i][j] = greengam[i][j];
                YNew[i][j] = bluegam[i][j];
            }
        }

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 16)
#endif
        for (int i = 0; i < bfh; ++i) {
            Color::RGB2L(redVal[i], greenVal[i], blueVal[i], L[i], wip, bfw);
            Color::RGB2Y(redVal[i], greenVal[i], blueVal[i], YOld[i], YNew[i], bfw);
        }
        float contrast = pow_F(sharpc / 100.f, 1.f) * s_scale;
        float reducautocontrast = 1.f;//to take noise into account

        buildBlendMask2(L, clipMask2, bfw, bfh, contrast, 1.f, autoshar, 2.f / s_scale, 1.f, reducautocontrast);
        sharpc = 100.f * pow_F(contrast, 1.f) / s_scale;

        if (settings->verbose) {
            printf("Contrast threshold SE Captur=%f \n", (double) sharpc);
        }

        CaptureDeconvSharpening_SE(YNew, YOld, clipMask2, bfw, bfh, locp, capradiu, deconvCo, deconvLat, itcheck, 0.2, 0.9);// 0.2 and 0.9 not used but in case.
 
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 16)
#endif
        for (int i = 0; i < bfh; ++i) {
            for (int j = 0; j < bfw; ++j) {
                const float factor = YNew[i][j] / std::max(YOld[i][j], 0.00001f);
                redgam[i][j] = redVal[i][j] * factor;
                greengam[i][j] = greenVal[i][j] * factor;
                bluegam[i][j] = blueVal[i][j] * factor;
            }
        }
    
        if (gamma1 != 1.f) {//inverse gamma 
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int i = 0; i < bfh; ++i) {
                for (int j = 0; j < bfw; ++j) {
                    redgam[i][j] = 65535.f * gammalog2(redgam[i][j] / 65535.f, gamma1, ts1, g_a[3], g_a[4]);
                    greengam[i][j] = 65535.f * gammalog2(greengam[i][j] / 65535.f, gamma1, ts1, g_a[3], g_a[4]);
                    bluegam[i][j] = 65535.f * gammalog2(bluegam[i][j] / 65535.f, gamma1, ts1, g_a[3], g_a[4]);
                }
            }
        }

#ifdef _OPENMP
         #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif               
        for (int i = 0; i < bfh; ++i) {
            for (int j = 0; j < bfw; ++j) {
                rgb->r(i,j) = redgam[i][j]; 
                rgb->g(i,j) = greengam[i][j]; 
                rgb->b(i,j) = bluegam[i][j]; 
            }
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
