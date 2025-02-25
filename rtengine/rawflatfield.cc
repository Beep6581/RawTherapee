/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2019 Gabor Horvath <hgabor@rawtherapee.com>
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

#include <cstring>
#include <memory>
#include <new>

#include "rawimagesource.h"
#include "procparams.h"
#include "rawimage.h"
//#define BENCHMARK
//#include "StopWatch.h"
#include "opthelper.h"

namespace {

void cfaboxblur(const float* const * riFlatFile, float* cfablur, int boxH, int boxW, int H, int W)
{
    if (boxW < 0 || boxH < 0 || (boxW == 0 && boxH == 0)) { // nothing to blur or negative values
        memcpy(cfablur, riFlatFile[0], static_cast<unsigned long>(W) * H * sizeof(float));
        return;
    }

    std::unique_ptr<float []> tmpBuffer;
    float *cfatmp = cfablur;


    if (boxH > 0 && boxW > 0) {
        // we need a temporary buffer if we have to blur both directions
        tmpBuffer.reset(new float [H * W]);
        cfatmp = tmpBuffer.get();
    }

    // if boxW == 0 we can skip the horizontal blur and process the vertical blur from riFlatFile to cfablur without using a temporary buffer
    const float* srcVertical = boxW == 0 ? riFlatFile[0] : cfatmp;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {

        if (boxW > 0) {
            //box blur cfa image; box size = BS
            //horizontal blur
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int row = 0; row < H; ++row) {
                int len = boxW / 2 + 1;
                cfatmp[row * W] = riFlatFile[row][0] / len;
                cfatmp[row * W + 1] = riFlatFile[row][1] / len;

                for (int j = 2; j <= boxW; j += 2) {
                    cfatmp[row * W] += riFlatFile[row][j] / len;
                    cfatmp[row * W + 1] += riFlatFile[row][j + 1] / len;
                }

                for (int col = 2; col <= boxW; col += 2) {
                    cfatmp[row * W + col] = (cfatmp[row * W + col - 2] * len + riFlatFile[row][boxW + col]) / (len + 1);
                    cfatmp[row * W + col + 1] = (cfatmp[row * W + col - 1] * len + riFlatFile[row][boxW + col + 1]) / (len + 1);
                    len ++;
                }

                const float rlen = 1.f / len;
                for (int col = boxW + 2; col < W - boxW; col++) {
                    cfatmp[row * W + col] = cfatmp[row * W + col - 2] + (riFlatFile[row][boxW + col] - cfatmp[row * W + col - boxW - 2]) * rlen;
                }

                for (int col = W - boxW; col < W; col += 2) {
                    cfatmp[row * W + col] = (cfatmp[row * W + col - 2] * len - cfatmp[row * W + col - boxW - 2]) / (len - 1);

                    if (col + 1 < W) {
                        cfatmp[row * W + col + 1] = (cfatmp[row * W + col - 1] * len - cfatmp[row * W + col - boxW - 1]) / (len - 1);
                    }

                    len --;
                }
            }
        }

        if (boxH > 0) {
            //vertical blur
#ifdef __SSE2__
            const vfloat leninitv = F2V(boxH / 2 + 1);
            const vfloat onev = F2V(1.f);
#ifdef _OPENMP
            #pragma omp for nowait
#endif

            for (int col = 0; col < W - 7; col += 8) {
                vfloat lenv = leninitv;
                vfloat temp1v = LVFU(srcVertical[col]) / lenv;
                vfloat temp2v = LVFU(srcVertical[W + col]) / lenv;
                vfloat temp3v = LVFU(srcVertical[col + 4]) / lenv;
                vfloat temp4v = LVFU(srcVertical[W + col + 4]) / lenv;

                for (int i = 2; i < boxH + 2; i += 2) {
                    temp1v += LVFU(srcVertical[i * W + col]) / lenv;
                    temp2v += LVFU(srcVertical[(i + 1) * W + col]) / lenv;
                    temp3v += LVFU(srcVertical[i * W + col + 4]) / lenv;
                    temp4v += LVFU(srcVertical[(i + 1) * W + col + 4]) / lenv;
                }

                STVFU(cfablur[col], temp1v);
                STVFU(cfablur[W + col], temp2v);
                STVFU(cfablur[col + 4], temp3v);
                STVFU(cfablur[W + col + 4], temp4v);

                int row;
                for (row = 2; row < boxH + 2; row += 2) {
                    const vfloat lenp1v = lenv + onev;
                    temp1v = (temp1v * lenv + LVFU(srcVertical[(row + boxH) * W + col])) / lenp1v;
                    temp2v = (temp2v * lenv + LVFU(srcVertical[(row + boxH + 1) * W + col])) / lenp1v;
                    temp3v = (temp3v * lenv + LVFU(srcVertical[(row + boxH) * W + col + 4])) / lenp1v;
                    temp4v = (temp4v * lenv + LVFU(srcVertical[(row + boxH + 1) * W + col + 4])) / lenp1v;
                    STVFU(cfablur[row * W + col], temp1v);
                    STVFU(cfablur[(row + 1) * W + col], temp2v);
                    STVFU(cfablur[row * W + col + 4], temp3v);
                    STVFU(cfablur[(row + 1) * W + col + 4], temp4v);
                    lenv = lenp1v;
                }

                for (; row < H - boxH - 1; row += 2) {
                    temp1v = temp1v + (LVFU(srcVertical[(row + boxH) * W + col]) - LVFU(srcVertical[(row - boxH - 2) * W + col])) / lenv;
                    temp2v = temp2v + (LVFU(srcVertical[(row + 1 + boxH) * W + col]) - LVFU(srcVertical[(row + 1 - boxH - 2) * W + col])) / lenv;
                    temp3v = temp3v + (LVFU(srcVertical[(row + boxH) * W + col + 4]) - LVFU(srcVertical[(row - boxH - 2) * W + col + 4])) / lenv;
                    temp4v = temp4v + (LVFU(srcVertical[(row + 1 + boxH) * W + col + 4]) - LVFU(srcVertical[(row + 1 - boxH - 2) * W + col + 4])) / lenv;
                    STVFU(cfablur[row * W + col], temp1v);
                    STVFU(cfablur[(row + 1) * W + col], temp2v);
                    STVFU(cfablur[row * W + col + 4], temp3v);
                    STVFU(cfablur[(row + 1) * W + col + 4], temp4v);
                }

                if (row < H - boxH) {
                    temp1v = temp1v + (LVFU(srcVertical[(row + boxH) * W + col]) - LVFU(srcVertical[(row - boxH - 2) * W + col])) / lenv;
                    temp3v = temp3v + (LVFU(srcVertical[(row + boxH) * W + col + 4]) - LVFU(srcVertical[(row - boxH - 2) * W + col + 4])) / lenv;
                    STVFU(cfablur[row * W + col], temp1v);
                    STVFU(cfablur[row * W + col + 4], temp3v);
                    vfloat swapv = temp1v;
                    temp1v = temp2v;
                    temp2v = swapv;
                    swapv = temp3v;
                    temp3v = temp4v;
                    temp4v = swapv;
                    ++row;
                }

                for (; row < H - 1; row += 2) {
                    const vfloat lenm1v = lenv - onev;
                    temp1v = (temp1v * lenv - LVFU(srcVertical[(row - boxH - 2) * W + col])) / lenm1v;
                    temp2v = (temp2v * lenv - LVFU(srcVertical[(row - boxH - 1) * W + col])) / lenm1v;
                    temp3v = (temp3v * lenv - LVFU(srcVertical[(row - boxH - 2) * W + col + 4])) / lenm1v;
                    temp4v = (temp4v * lenv - LVFU(srcVertical[(row - boxH - 1) * W + col + 4])) / lenm1v;
                    STVFU(cfablur[row * W + col], temp1v);
                    STVFU(cfablur[(row + 1) * W + col], temp2v);
                    STVFU(cfablur[row * W + col + 4], temp3v);
                    STVFU(cfablur[(row + 1) * W + col + 4], temp4v);
                    lenv = lenm1v;
                }

                if (row < H) {
                    vfloat lenm1v = lenv - onev;
                    temp1v = (temp1v * lenv - LVFU(srcVertical[(row - boxH - 2) * W + col])) / lenm1v;
                    temp3v = (temp3v * lenv - LVFU(srcVertical[(row - boxH - 2) * W + col + 4])) / lenm1v;
                    STVFU(cfablur[row * W + col], temp1v);
                    STVFU(cfablur[row * W + col + 4], temp3v);
                }

            }

#ifdef _OPENMP
            #pragma omp single
#endif

            for (int col = W - (W % 8); col < W; ++col) {
                int len = boxH / 2 + 1;
                cfablur[col] = srcVertical[col] / len;
                cfablur[W + col] = srcVertical[W + col] / len;

                for (int i = 2; i < boxH + 2; i += 2) {
                    cfablur[col] += srcVertical[i * W + col] / len;
                    cfablur[W + col] += srcVertical[(i + 1) * W + col] / len;
                }

                for (int row = 2; row < boxH + 2; row += 2) {
                    cfablur[row * W + col] = (cfablur[(row - 2) * W + col] * len + srcVertical[(row + boxH) * W + col]) / (len + 1);
                    cfablur[(row + 1) * W + col] = (cfablur[(row - 1) * W + col] * len + srcVertical[(row + boxH + 1) * W + col]) / (len + 1);
                    ++len;
                }

                for (int row = boxH + 2; row < H - boxH; ++row) {
                    cfablur[row * W + col] = cfablur[(row - 2) * W + col] + (srcVertical[(row + boxH) * W + col] - srcVertical[(row - boxH - 2) * W + col]) / len;
                }

                for (int row = H - boxH; row < H; row += 2) {
                    cfablur[row * W + col] = (cfablur[(row - 2) * W + col] * len - srcVertical[(row - boxH - 2) * W + col]) / (len - 1);

                    if (row + 1 < H) {
                        cfablur[(row + 1) * W + col] = (cfablur[(row - 1) * W + col] * len - srcVertical[(row - boxH - 1) * W + col]) / (len - 1);
                    }
                    --len;
                }
            }

#else
#ifdef _OPENMP
            #pragma omp for
#endif

            for (int col = 0; col < W; ++col) {
                int len = boxH / 2 + 1;
                cfablur[col] = srcVertical[col] / len;
                cfablur[W + col] = srcVertical[W + col] / len;

                for (int i = 2; i < boxH + 2; i += 2) {
                    cfablur[col] += srcVertical[i * W + col] / len;
                    cfablur[W + col] += srcVertical[(i + 1) * W + col] / len;
                }

                for (int row = 2; row < boxH + 2; row += 2) {
                    cfablur[row * W + col] = (cfablur[(row - 2) * W + col] * len + srcVertical[(row + boxH) * W + col]) / (len + 1);
                    cfablur[(row + 1) * W + col] = (cfablur[(row - 1) * W + col] * len + srcVertical[(row + boxH + 1) * W + col]) / (len + 1);
                    ++len;
                }

                for (int row = boxH + 2; row < H - boxH; ++row) {
                    cfablur[row * W + col] = cfablur[(row - 2) * W + col] + (srcVertical[(row + boxH) * W + col] - srcVertical[(row - boxH - 2) * W + col]) / len;
                }

                for (int row = H - boxH; row < H; row += 2) {
                    cfablur[row * W + col] = (cfablur[(row - 2) * W + col] * len - srcVertical[(row - boxH - 2) * W + col]) / (len - 1);

                    if (row + 1 < H) {
                        cfablur[(row + 1) * W + col] = (cfablur[(row - 1) * W + col] * len - srcVertical[(row - boxH - 1) * W + col]) / (len - 1);
                    }
                    --len;
                }
            }
#endif
        }
    }
}

}

namespace rtengine
{

void RawImageSource::processFlatField(const procparams::RAWParams &raw, RawImage *riFlatFile, array2D<float> &rawData, const float black[4])
{
//    BENCHFUN
    std::unique_ptr<float[]> cfablur(new float[H * W]);

    const int BS = raw.ff_BlurRadius + (raw.ff_BlurRadius & 1);

    std::array<float, 4> ffblack;
    {
        const auto tmpfilters = riFlatFile->get_filters();
        riFlatFile->set_filters(riFlatFile->prefilters); // we need 4 blacks for bayer processing
        riFlatFile->get_colorsCoeff(nullptr, nullptr, ffblack.data(), false);
        riFlatFile->set_filters(tmpfilters);
    }

    if (raw.ff_BlurType == procparams::RAWParams::getFlatFieldBlurTypeString(procparams::RAWParams::FlatFieldBlurType::V)) {
        cfaboxblur(riFlatFile->data, cfablur.get(), 2 * BS, 0, H, W);
    } else if (raw.ff_BlurType == procparams::RAWParams::getFlatFieldBlurTypeString(procparams::RAWParams::FlatFieldBlurType::H)) {
        cfaboxblur(riFlatFile->data, cfablur.get(), 0, 2 * BS, H, W);
    } else if (raw.ff_BlurType == procparams::RAWParams::getFlatFieldBlurTypeString(procparams::RAWParams::FlatFieldBlurType::VH)) {
        //slightly more complicated blur if trying to correct both vertical and horizontal anomalies
        cfaboxblur(riFlatFile->data, cfablur.get(), BS, BS, H, W);    //first do area blur to correct vignette
    } else { //(raw.ff_BlurType == RAWParams::getFlatFieldBlurTypeString(RAWParams::area_ff))
        cfaboxblur(riFlatFile->data, cfablur.get(), BS, BS, H, W);
    }

    if (ri->getSensorType() == ST_BAYER || ri->get_colors() == 1) {
        float refcolor[2][2];

        // find center values by channel
        for (int m = 0; m < 2; ++m)
            for (int n = 0; n < 2; ++n) {
                const int row = 2 * (H >> 2) + m;
                const int col = 2 * (W >> 2) + n;
                const int c  = ri->get_colors() != 1 ? FC(row, col) : 0;
                const int c4 = ri->get_colors() != 1 ? ((c == 1 && !(row & 1)) ? 3 : c) : 0;
                refcolor[m][n] = std::max(0.0f, cfablur[row * W + col] - ffblack[c4]);
            }

        float limitFactor = 1.f;

        if (raw.ff_AutoClipControl) {
            bool clippedBefore = false;
            for (int m = 0; m < 2 && !clippedBefore; ++m) {
                for (int n = 0; n < 2 && !clippedBefore; ++n) {
                    float maxval = 0.f;
                    const int c  = ri->get_colors() != 1 ? FC(m, n) : 0;
                    const int c4 = ri->get_colors() != 1 ? ((c == 1 && !(m & 1)) ? 3 : c) : 0;
                    const float clipVal = ri->get_white(c4);
#ifdef _OPENMP
                    #pragma omp parallel for reduction(max:maxval) schedule(dynamic, 16)
#endif
                    for (int row = 0; row < H - m; row += 2) {
                        for (int col = 0; col < W - n && !clippedBefore; col += 2) {
                            const float rawVal = rawData[row + m][col + n];
                            if (rawVal >= clipVal) {
                                clippedBefore = true;
                                break;
                            }
                            const float tempval = (rawVal - black[c4]) * (refcolor[m][n] / std::max(1e-5f, cfablur[(row + m) * W + col + n] - ffblack[c4]));
                            maxval = std::max(maxval, tempval);
                        }
                    }

                    // now we have the max value for the channel
                    // if it clips, calculate factor to avoid clipping
                    if (maxval + black[c4] >= ri->get_white(c4)) {
                        if (!clippedBefore) {
                            limitFactor = std::min(limitFactor, ri->get_white(c4) / (maxval + black[c4]));
                        } else {
                            limitFactor = 1.f;
                        }
                    }
                }
            }
            flatFieldAutoClipValue = (1.f - limitFactor) * 100.f;           // this value can be used to set the clip control slider in gui
        } else {
            limitFactor = std::max((100 - raw.ff_clipControl) / 100.f, 0.01f);
        }

        for (int m = 0; m < 2; ++m)
            for (int n = 0; n < 2; ++n) {
                refcolor[m][n] *= limitFactor;
            }

        unsigned int c[2][2] {};
        unsigned int c4[2][2] {};
        if (ri->get_colors() != 1) {
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    c[i][j] = FC(i, j);
                }
            }
            c4[0][0] = (c[0][0] == 1) ? 3 : c[0][0];
            c4[0][1] = (c[0][1] == 1) ? 3 : c[0][1];
            c4[1][0] = c[1][0];
            c4[1][1] = c[1][1];
        }

        constexpr float minValue = 1.f; // if the pixel value in the flat field is less or equal this value, no correction will be applied.

#ifdef __SSE2__
        const vfloat refcolorv[2] = {_mm_set_ps(refcolor[0][1], refcolor[0][0], refcolor[0][1], refcolor[0][0]),
                                     _mm_set_ps(refcolor[1][1], refcolor[1][0], refcolor[1][1], refcolor[1][0])
                                    };
        const vfloat blackv[2] = {_mm_set_ps(black[c4[0][1]], black[c4[0][0]], black[c4[0][1]], black[c4[0][0]]),
                                  _mm_set_ps(black[c4[1][1]], black[c4[1][0]], black[c4[1][1]], black[c4[1][0]])
                                 };
        const vfloat ffblackv[2] = {_mm_set_ps(ffblack[c4[0][1]], ffblack[c4[0][0]], ffblack[c4[0][1]], ffblack[c4[0][0]]),
                                    _mm_set_ps(ffblack[c4[1][1]], ffblack[c4[1][0]], ffblack[c4[1][1]], ffblack[c4[1][0]])
                                   };

        const vfloat onev = F2V(1.f);
        const vfloat minValuev = F2V(minValue);
#endif
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

        for (int row = 0; row < H; ++row) {
            int col = 0;
#ifdef __SSE2__
            const vfloat rowBlackv = blackv[row & 1];
            const vfloat ffrowBlackv = ffblackv[row & 1];
            const vfloat rowRefcolorv = refcolorv[row & 1];

            for (; col < W - 3; col += 4) {
                const vfloat blurv = LVFU(cfablur[row * W + col]) - ffrowBlackv;
                vfloat vignettecorrv = rowRefcolorv / blurv;
                vignettecorrv = vself(vmaskf_le(blurv, minValuev), onev, vignettecorrv);
                const vfloat valv = LVFU(rawData[row][col]) - rowBlackv;
                STVFU(rawData[row][col], valv * vignettecorrv + rowBlackv);
            }

#endif

            for (; col < W; ++col) {
                const float blur = cfablur[row * W + col] - ffblack[c4[row & 1][col & 1]];
                const float vignettecorr = blur <= minValue ? 1.f : refcolor[row & 1][col & 1] / blur;
                rawData[row][col] = (rawData[row][col] - black[c4[row & 1][col & 1]]) * vignettecorr + black[c4[row & 1][col & 1]];
            }
        }
    } else if (ri->getSensorType() == ST_FUJI_XTRANS) {
        float refcolor[3] = {0.f};
        int cCount[3] = {0};

        // find center average values by channel
        for (int m = -3; m < 3; ++m)
            for (int n = -3; n < 3; ++n) {
                const int row = 2 * (H >> 2) + m;
                const int col = 2 * (W >> 2) + n;
                const int c  = riFlatFile->XTRANSFC(row, col);
                refcolor[c] += std::max(0.0f, cfablur[row * W + col] - black[c]);
                cCount[c] ++;
            }

        for (int c = 0; c < 3; ++c) {
            refcolor[c] = refcolor[c] / cCount[c];
        }

        float limitFactor = 1.f;

        if (raw.ff_AutoClipControl) {
            // determine maximum calculated value to avoid clipping
            bool clippedBefore = false;
            const float clipVal = ri->get_white(0);
            float maxval = 0.f;
            // xtrans files have only one black level actually, so we can simplify the code a bit
#ifdef _OPENMP
            #pragma omp parallel for reduction(max:maxval) schedule(dynamic,16)
#endif
            for (int row = 0; row < H; ++row) {
                for (int col = 0; col < W && !clippedBefore; ++col) {
                    const float rawVal = rawData[row][col];
                    if (rawVal >= clipVal) {
                        clippedBefore = true;
                        break;
                    }
                    const float tempval = (rawVal - black[0]) * (refcolor[ri->XTRANSFC(row, col)] / std::max(1e-5f, cfablur[(row) * W + col] - black[0]));
                    maxval = std::max(maxval, tempval);
                }
            }

            // there's only one white level for xtrans
            if (!clippedBefore && maxval + black[0] > ri->get_white(0)) {
                limitFactor = ri->get_white(0) / (maxval + black[0]);
                flatFieldAutoClipValue = (1.f - limitFactor) * 100.f;           // this value can be used to set the clip control slider in gui
            }
        } else {
            limitFactor = std::max((float)(100 - raw.ff_clipControl) / 100.f, 0.01f);
        }


        for (int c = 0; c < 3; ++c) {
            refcolor[c] *= limitFactor;
        }

        constexpr float minValue = 1.f; // if the pixel value in the flat field is less or equal this value, no correction will be applied.

#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int row = 0; row < H; ++row) {
            for (int col = 0; col < W; ++col) {
                const int c = ri->XTRANSFC(row, col);
                const float blur = cfablur[(row) * W + col] - black[c];
                const float vignettecorr = blur <= minValue ? 1.f : refcolor[c] / blur;
                rawData[row][col] = (rawData[row][col] - black[c]) * vignettecorr + black[c];
            }
        }
    }

    if (raw.ff_BlurType == procparams::RAWParams::getFlatFieldBlurTypeString(procparams::RAWParams::FlatFieldBlurType::VH)) {
        std::unique_ptr<float []> cfablur1(new float[H * W]);
        std::unique_ptr<float []> cfablur2(new float[H * W]);
        //slightly more complicated blur if trying to correct both vertical and horizontal anomalies
        cfaboxblur(riFlatFile->data, cfablur1.get(), 0, 2 * BS, H, W); //now do horizontal blur
        cfaboxblur(riFlatFile->data, cfablur2.get(), 2 * BS, 0, H, W); //now do vertical blur

        if (ri->getSensorType() == ST_BAYER || ri->get_colors() == 1) {
            unsigned int c[2][2] {};
            unsigned int c4[2][2] {};
            if (ri->get_colors() != 1) {
                for (int i = 0; i < 2; ++i) {
                    for (int j = 0; j < 2; ++j) {
                        c[i][j] = FC(i, j);
                    }
                }
                c4[0][0] = (c[0][0] == 1) ? 3 : c[0][0];
                c4[0][1] = (c[0][1] == 1) ? 3 : c[0][1];
                c4[1][0] = c[1][0];
                c4[1][1] = c[1][1];
            }

#ifdef __SSE2__
            const vfloat blackv[2] = {_mm_set_ps(black[c4[0][1]], black[c4[0][0]], black[c4[0][1]], black[c4[0][0]]),
                                      _mm_set_ps(black[c4[1][1]], black[c4[1][0]], black[c4[1][1]], black[c4[1][0]])
                                     };
            const vfloat ffblackv[2] = {_mm_set_ps(ffblack[c4[0][1]], ffblack[c4[0][0]], ffblack[c4[0][1]], ffblack[c4[0][0]]),
                                        _mm_set_ps(ffblack[c4[1][1]], ffblack[c4[1][0]], ffblack[c4[1][1]], ffblack[c4[1][0]])
                                       };


            const vfloat epsv = F2V(1e-5f);
#endif
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16)
#endif

            for (int row = 0; row < H; ++row) {
                int col = 0;
#ifdef __SSE2__
                const vfloat rowBlackv = blackv[row & 1];
                const vfloat ffrowBlackv = ffblackv[row & 1];

                for (; col < W - 3; col += 4) {
                    const vfloat linecorrv = SQRV(vmaxf(LVFU(cfablur[row * W + col]) - ffrowBlackv, epsv)) /
                                             (vmaxf(LVFU(cfablur1[row * W + col]) - ffrowBlackv, epsv) * vmaxf(LVFU(cfablur2[row * W + col]) - ffrowBlackv, epsv));
                    const vfloat valv = LVFU(rawData[row][col]) - rowBlackv;
                    STVFU(rawData[row][col], valv * linecorrv + rowBlackv);
                }

#endif

                for (; col < W; ++col) {
                    const float linecorr = SQR(std::max(1e-5f, cfablur[row * W + col] - ffblack[c4[row & 1][col & 1]])) /
                                           (std::max(1e-5f, cfablur1[row * W + col] - ffblack[c4[row & 1][col & 1]]) * std::max(1e-5f, cfablur2[row * W + col] - ffblack[c4[row & 1][col & 1]]));
                    rawData[row][col] = (rawData[row][col] - black[c4[row & 1][col & 1]]) * linecorr + black[c4[row & 1][col & 1]];
                }
            }
        } else if (ri->getSensorType() == ST_FUJI_XTRANS) {
#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int row = 0; row < H; ++row) {
                for (int col = 0; col < W; ++col) {
                    const int c  = ri->XTRANSFC(row, col);
                    const float hlinecorr = std::max(1e-5f, cfablur[(row) * W + col] - black[c]) / std::max(1e-5f, cfablur1[(row) * W + col] - black[c]);
                    const float vlinecorr = std::max(1e-5f, cfablur[(row) * W + col] - black[c]) / std::max(1e-5f, cfablur2[(row) * W + col] - black[c]);
                    rawData[row][col] = (rawData[row][col] - black[c]) * hlinecorr * vlinecorr + black[c];
                }
            }
        }
    }
}
} /* namespace */
