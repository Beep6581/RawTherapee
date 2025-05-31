////////////////////////////////////////////////////////////////
//
//      Highlight reconstruction
//
//      copyright (c) 2008-2011  Emil Martinec <ejmartin@uchicago.edu>
//      copyright (c) 2019 Ingo Weyrich <heckflosse67@gmx.de>
//
//
// code dated: June 16, 2011
// code dated: July 09, 2019, speedups by Ingo Weyrich <heckflosse67@gmx.de>
//
//  hilite_recon.cc is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////

#include <cassert>
#include <cmath>
#include <cstddef>

#include "array2D.h"
#include "opthelper.h"
#include "rawimagesource.h"
#include "rt_math.h"
//#define BENCHMARK
#include "StopWatch.h"
#include "guidedfilter.h"
#include "settings.h"
#include "gauss.h"
#include "rescale.h"
#include "iccstore.h"
#include "color.h"
#include "linalgebra.h"

namespace
{

void boxblur2(const float* const* src, float** dst, float** temp, int startY, int startX, int H, int W, int bufferW, int box)
{
    constexpr int numCols = 16;
    assert((bufferW % numCols) == 0);

    //box blur image channel; box size = 2*box+1
    //horizontal blur
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int row = 0; row < H; ++row) {
        int len = box + 1;
        temp[row][0] = src[row + startY][startX] / len;

        for (int j = 1; j <= box; ++j) {
            temp[row][0] += src[row + startY][j + startX] / len;
        }

        for (int col = 1; col <= box; ++col, ++len) {
            temp[row][col] = (temp[row][col - 1] * len + src[row + startY][col + box + startX]) / (len + 1);
        }

        for (int col = box + 1; col < W - box; ++col) {
            temp[row][col] = temp[row][col - 1] + (src[row + startY][col + box + startX] - src[row + startY][col - box - 1 + startX]) / len;
        }

        for (int col = W - box; col < W; ++col, --len) {
            temp[row][col] = (temp[row][col - 1] * len - src[row + startY][col - box - 1 + startX]) / (len - 1);
        }
    }

    //vertical blur
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        float tempvalN[numCols] ALIGNED64;
#ifdef _OPENMP
    #pragma omp for
#endif
        for (int col = 0; col < bufferW - numCols + 1; col += numCols) {
            float len = box + 1;

            for (int n = 0; n < numCols; ++n) {
                tempvalN[n] = temp[0][col + n] / len;
            }

            for (int i = 1; i <= box; ++i) {
                for (int n = 0; n < numCols; ++n) {
                    tempvalN[n] += temp[i][col + n] / len;
                }
            }

            for (int n = 0; n < numCols; ++n) {
                dst[0][col + n] = tempvalN[n];
            }

            for (int row = 1; row <= box; ++row, ++len) {
                for (int n = 0; n < numCols; ++n) {
                    tempvalN[n] = (tempvalN[n] * len + temp[(row + box)][col + n]) / (len + 1);
                    dst[row][col + n] = tempvalN[n];
                }
            }

            const float rlen = 1.f / len;

            for (int row = box + 1; row < H - box; ++row) {
                for (int n = 0; n < numCols; ++n) {
                    tempvalN[n] += (temp[(row + box)][col + n] - temp[(row - box - 1)][col + n]) * rlen;
                    dst[row][col + n] = tempvalN[n];
                }
            }

            for (int row = H - box; row < H; ++row, --len) {
                for (int n = 0; n < numCols; ++n) {
                    tempvalN[n] = (dst[(row - 1)][col + n] * len - temp[(row - box - 1)][col + n]) / (len - 1);
                    dst[row][col + n] = tempvalN[n];
                }
            }
        }
    }
}

void boxblur_resamp(const float* const* src, float** dst, float** temp, int H, int W, int box, int samp)
{
    assert(samp != 0);

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef _OPENMP
        #pragma omp for
#endif
        //box blur image channel; box size = 2*box+1
        //horizontal blur
        for (int row = 0; row < H; ++row) {
            int len = box + 1;
            float tempval = src[row][0] / len;

            for (int j = 1; j <= box; ++j) {
                tempval += src[row][j] / len;
            }

            temp[row][0] = tempval;

            for (int col = 1; col <= box; ++col, ++len) {
                tempval = (tempval * len + src[row][col + box]) / (len + 1);

                if (col % samp == 0) {
                    temp[row][col / samp] = tempval;
                }
            }

            const float oneByLen = 1.f / static_cast<float>(len);

            for (int col = box + 1; col < W - box; ++col) {
                tempval = tempval + (src[row][col + box] - src[row][col - box - 1]) * oneByLen;

                if (col % samp == 0) {
                    temp[row][col / samp] = tempval;
                }
            }

            for (int col = W - box; col < W; ++col, --len) {
                tempval = (tempval * len - src[row][col - box - 1]) / (len - 1);

                if (col % samp == 0) {
                    temp[row][col / samp] = tempval;
                }
            }
        }
    }

    constexpr int numCols = 8;   // process numCols columns at once for better L1 CPU cache usage

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        float tempvalN[numCols] ALIGNED64;

#ifdef _OPENMP
        #pragma omp for nowait
#endif
        //vertical blur
        for (int col = 0; col < (W / samp) - (numCols - 1); col += numCols) {
            float len = box + 1;

            for (int n = 0; n < numCols; ++n) {
                tempvalN[n] = temp[0][col + n] / len;
            }

            for (int i = 1; i <= box; ++i) {
                for (int n = 0; n < numCols; ++n) {
                    tempvalN[n] += temp[i][col + n] / len;
                }
            }

            for (int n = 0; n < numCols; ++n) {
                dst[0][col + n] = tempvalN[n];
            }

            for (int row = 1; row <= box; ++row, ++len) {
                for (int n = 0; n < numCols; ++n) {
                    tempvalN[n] = (tempvalN[n] * len + temp[(row + box)][col + n]) / (len + 1);
                }

                if (row % samp == 0) {
                    for (int n = 0; n < numCols; ++n) {
                        dst[row / samp][col + n] = tempvalN[n];
                    }
                }
            }

            const float rlen = 1.f / len;

            for (int row = box + 1; row < H - box; ++row) {
                for (int n = 0; n < numCols; ++n) {
                    tempvalN[n] += (temp[(row + box)][col + n] - temp[(row - box - 1)][col + n]) * rlen;
                }

                if (row % samp == 0) {
                    for (int n = 0; n < numCols; ++n) {
                        dst[row / samp][col + n] = tempvalN[n];
                    }
                }
            }

            for (int row = H - box; row < H; ++row, --len) {
                for (int n = 0; n < numCols; ++n) {
                    tempvalN[n] = (tempvalN[n] * len - temp[(row - box - 1)][col + n]) / (len - 1);
                }

                if (row % samp == 0) {
                    for (int n = 0; n < numCols; ++n) {
                        dst[row / samp][col + n] = tempvalN[n];
                    }
                }
            }
        }

        // process remaining columns
#ifdef _OPENMP
        #pragma omp single
#endif
        {

            //vertical blur
            for (int col = (W / samp) - ((W / samp) % numCols); col < W / samp; ++col) {
                int len = box + 1;
                float tempval = temp[0][col] / len;

                for (int i = 1; i <= box; ++i) {
                    tempval += temp[i][col] / len;
                }

                dst[0][col] = tempval;

                for (int row = 1; row <= box; ++row, ++len) {
                    tempval = (tempval * len + temp[(row + box)][col]) / (len + 1);

                    if (row % samp == 0) {
                        dst[row / samp][col] = tempval;
                    }
                }

                for (int row = box + 1; row < H - box; ++row) {
                    tempval += (temp[(row + box)][col] - temp[(row - box - 1)][col]) / len;

                    if (row % samp == 0) {
                        dst[row / samp][col] = tempval;
                    }
                }

                for (int row = H - box; row < H; ++row, --len) {
                    tempval = (tempval * len - temp[(row - box - 1)][col]) / (len - 1);

                    if (row % samp == 0) {
                        dst[row / samp][col] = tempval;
                    }
                }
            }
        }
    }
}

}

namespace rtengine
{
extern const Settings *settings;
using namespace procparams;
    const ProcParams params;

void RawImageSource::HLRecovery_inpaint(float** red, float** green, float** blue, int blur)
{  
    //BENCHFUN
    double progress = 0.0;

    if (plistener) {
        plistener->setProgressStr("PROGRESSBAR_HLREC");
        plistener->setProgress(progress);
    }

    const int height = H;
    const int width = W;

    constexpr int range = 2;
    constexpr int pitch = 4;

    constexpr float threshpct = 0.25f;
    constexpr float maxpct = 0.95f;
    constexpr float epsilon = 0.00001f;

    //for blend algorithm:
    constexpr float blendthresh = 1.0;
    // Transform matrixes rgb>lab and back
    constexpr float trans[3][3] = {
        {1.f, 1.f, 1.f},
        {1.7320508f, -1.7320508f, 0.f},
        {-1.f, -1.f, 2.f}
    };
    constexpr float itrans[3][3] = {
        {1.f, 0.8660254f, -0.5f},
        {1.f, -0.8660254f, -0.5f},
        {1.f, 0.f, 1.f}
    };

    if (settings->verbose) {
        for (int c = 0; c < 3; ++c) {
            printf("chmax[%d] : %f\tclmax[%d] : %f\tratio[%d] : %f\n", c, static_cast<double>(chmax[c]), c, static_cast<double>(clmax[c]), c, static_cast<double>(chmax[c] / clmax[c]));
        }
    }

    float factor[3];

    for (int c = 0; c < 3; ++c) {
        factor[c] = chmax[c] / clmax[c];
    }

    const float minFactor = min(factor[0], factor[1], factor[2]);

    if (minFactor > 1.f) { // all 3 channels clipped
        // calculate clip factor per channel
        for (int c = 0; c < 3; ++c) {
            factor[c] /= minFactor;
        }

        // get max clip factor
        int maxpos = 0;
        float maxValNew = 0.f;

        for (int c = 0; c < 3; ++c) {
            if (chmax[c] / factor[c] > maxValNew) {
                maxValNew = chmax[c] / factor[c];
                maxpos = c;
            }
        }

        const float clipFactor = clmax[maxpos] / maxValNew;

        if (clipFactor < maxpct) {
            // if max clipFactor < maxpct (0.95) adjust per channel factors
            for (int c = 0; c < 3; ++c) {
                factor[c] *= (maxpct / clipFactor);
            }
        }
    } else {
        factor[0] = factor[1] = factor[2] = 1.f;
    }

    if (settings->verbose) {
        for (int c = 0; c < 3; ++c) {
            printf("correction factor[%d] : %f\n", c, static_cast<double>(factor[c]));
        }
    }

    float max_f[3];
    float thresh[3];

    for (int c = 0; c < 3; ++c) {
        thresh[c] = chmax[c] * threshpct / factor[c];
        max_f[c] = chmax[c] * maxpct / factor[c];
    }

    const float whitept = max(max_f[0], max_f[1], max_f[2]);
    const float clippt  = min(max_f[0], max_f[1], max_f[2]);
    const float medpt   = max_f[0] + max_f[1] + max_f[2] - whitept - clippt;
    const float blendpt = blendthresh * clippt;

    float medFactor[3];

    for (int c = 0; c < 3; ++c) {
        medFactor[c] = max(1.0f, max_f[c] / medpt) / -blendpt;
    }

    int minx = width - 1;
    int maxx = 0;
    int miny = height - 1;
    int maxy = 0;

#ifdef _OPENMP
    #pragma omp parallel for reduction(min:minx,miny) reduction(max:maxx,maxy) schedule(dynamic, 16)
#endif
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j< width; ++j) {
            if (red[i][j] >= max_f[0] || green[i][j] >= max_f[1] || blue[i][j] >= max_f[2]) {
                minx = std::min(minx, j);
                maxx = std::max(maxx, j);
                miny = std::min(miny, i);
                maxy = std::max(maxy, i);
            }
        }
    }

    if (minx > maxx || miny > maxy) { // nothing to reconstruct
        return;
    }

    constexpr int blurBorder = 256;
    minx = std::max(0, minx - blurBorder);
    miny = std::max(0, miny - blurBorder);
    maxx = std::min(width - 1, maxx + blurBorder);
    maxy = std::min(height - 1, maxy + blurBorder);
    const int blurWidth = maxx - minx + 1;
    const int blurHeight = maxy - miny + 1;
    const int bufferWidth = blurWidth + ((16 - (blurWidth % 16)) & 15);

    multi_array2D<float, 3> channelblur(bufferWidth, blurHeight, 0, 48);
    array2D<float> temp(bufferWidth, blurHeight); // allocate temporary buffer

    // blur RGB channels
    boxblur2(red, channelblur[0], temp, miny, minx, blurHeight, blurWidth, bufferWidth, 4);
    boxblur2(green, channelblur[1], temp, miny, minx, blurHeight, blurWidth, bufferWidth, 4);
    boxblur2(blue, channelblur[2], temp, miny, minx, blurHeight, blurWidth, bufferWidth, 4);
 
    if (plistener) {
        progress += 0.07;
        plistener->setProgress(progress);
    }

    // reduce channel blur to one array
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < blurHeight; ++i) {
        for (int j = 0; j < blurWidth; ++j) {
            channelblur[0][i][j] = std::fabs(channelblur[0][i][j] - red[i + miny][j + minx]) + std::fabs(channelblur[1][i][j] - green[i + miny][j + minx]) + std::fabs(channelblur[2][i][j] - blue[i + miny][j + minx]);
        }
    }

    for (int c = 1; c < 3; ++c) {
        channelblur[c].free();    //free up some memory
    }

    if (plistener) {
        progress += 0.05;
        plistener->setProgress(progress);
    }

    multi_array2D<float, 4> hilite_full(bufferWidth, blurHeight, ARRAY2D_CLEAR_DATA, 32);

    if (plistener) {
        progress += 0.05;
        plistener->setProgress(progress);
    }

    double hipass_sum = 0.0;
    int hipass_norm = 0;

    // set up which pixels are clipped or near clipping
#ifdef _OPENMP
    #pragma omp parallel for reduction(+:hipass_sum,hipass_norm) schedule(dynamic,16)
#endif
    for (int i = 0; i < blurHeight; ++i) {
        for (int j = 0; j < blurWidth; ++j) {
            if (
                (
                    red[i + miny][j + minx] > thresh[0]
                    || green[i + miny][j + minx] > thresh[1]
                    || blue[i + miny][j + minx] > thresh[2]
                )
                && red[i + miny][j + minx] < max_f[0]
                && green[i + miny][j + minx] < max_f[1]
                && blue[i + miny][j + minx] < max_f[2]
            ) {
                // if one or more channels is highlight but none are blown, add to highlight accumulator
                hipass_sum += static_cast<double>(channelblur[0][i][j]);
                ++hipass_norm;

                hilite_full[0][i][j] = red[i + miny][j + minx];
                hilite_full[1][i][j] = green[i + miny][j + minx];
                hilite_full[2][i][j] = blue[i + miny][j + minx];
                hilite_full[3][i][j] = 1.f;
            }
        }
    }

    const float hipass_ave = 2.0 * hipass_sum / (hipass_norm + static_cast<double>(epsilon));

    if (plistener) {
        progress += 0.05;
        plistener->setProgress(progress);
    }

    array2D<float> hilite_full4(bufferWidth, blurHeight);

    //blur highlight data
    boxblur2(hilite_full[3], hilite_full4, temp, 0, 0, blurHeight, blurWidth, bufferWidth, 1);

    temp.free(); // free temporary buffer

    if (plistener) {
        progress += 0.07;
        plistener->setProgress(progress);
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif
    for (int i = 0; i < blurHeight; ++i) {
        for (int j = 0; j < blurWidth; ++j) {
            if (channelblur[0][i][j] > hipass_ave) {
                //too much variation
                hilite_full[0][i][j] = hilite_full[1][i][j] = hilite_full[2][i][j] = hilite_full[3][i][j] = 0.f;
                continue;
            }

            if (hilite_full4[i][j] > epsilon && hilite_full4[i][j] < 0.95f) {
                //too near an edge, could risk using CA affected pixels, therefore omit
                hilite_full[0][i][j] = hilite_full[1][i][j] = hilite_full[2][i][j] = hilite_full[3][i][j] = 0.f;
            }
        }
    }

    channelblur[0].free();    //free up some memory
    hilite_full4.free();    //free up some memory

    const int hfh = (blurHeight - blurHeight % pitch) / pitch;
    const int hfw = (blurWidth - blurWidth % pitch) / pitch;

    multi_array2D<float, 4> hilite(hfw + 1, hfh + 1, ARRAY2D_CLEAR_DATA, 48);

    // blur and resample highlight data; range=size of blur, pitch=sample spacing
    array2D<float> temp2(blurWidth / pitch + (blurWidth % pitch == 0 ? 0 : 1), blurHeight);

    for (int m = 0; m < 4; ++m) {
        boxblur_resamp(hilite_full[m], hilite[m], temp2, blurHeight, blurWidth, range, pitch);

        if (plistener) {
            progress += 0.05;
            plistener->setProgress(progress);
        }
    }

    temp2.free();

    for (int c = 0; c < 4; ++c) {
        hilite_full[c].free();    //free up some memory
    }

    multi_array2D<float, 8> hilite_dir(hfw, hfh, ARRAY2D_CLEAR_DATA, 64);
    // for faster processing we create two buffers using (height,width) instead of (width,height)
    multi_array2D<float, 4> hilite_dir0(hfh, hfw, ARRAY2D_CLEAR_DATA, 64);
    multi_array2D<float, 4> hilite_dir4(hfh, hfw, ARRAY2D_CLEAR_DATA, 64);

    if (plistener) {
        progress += 0.05;
        plistener->setProgress(progress);
    }

    //fill gaps in highlight map by directional extension
    //raster scan from four corners
    for (int j = 1; j < hfw - 1; ++j) {
        for (int i = 2; i < hfh - 2; ++i) {
            //from left
            if (hilite[3][i][j] > epsilon) {
                hilite_dir0[3][j][i] = 1.f;
            } else {
                hilite_dir0[3][j][i] = (hilite_dir0[0 + 3][j - 1][i - 2] + hilite_dir0[0 + 3][j - 1][i - 1] + hilite_dir0[0 + 3][j - 1][i] + hilite_dir0[0 + 3][j - 1][i + 1] + hilite_dir0[0 + 3][j - 1][i + 2]) == 0.f ? 0.f : 0.1f;
            }
        }

        if (hilite[3][2][j] <= epsilon) {
            hilite_dir[0 + 3][0][j]  = hilite_dir0[3][j][2];
        }

        if (hilite[3][3][j] <= epsilon) {
            hilite_dir[0 + 3][1][j]  = hilite_dir0[3][j][3];
        }

        if (hilite[3][hfh - 3][j] <= epsilon) {
            hilite_dir[4 + 3][hfh - 1][j] = hilite_dir0[3][j][hfh - 3];
        }

        if (hilite[3][hfh - 4][j] <= epsilon) {
            hilite_dir[4 + 3][hfh - 2][j] = hilite_dir0[3][j][hfh - 4];
        }
    }

    for (int i = 2; i < hfh - 2; ++i) {
        if (hilite[3][i][hfw - 2] <= epsilon) {
            hilite_dir4[3][hfw - 1][i] = hilite_dir0[3][hfw - 2][i];
        }
    }

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef _OPENMP
        #pragma omp for nowait
#endif
        for (int c = 0; c < 3; ++c) {
            for (int j = 1; j < hfw - 1; ++j) {
                for (int i = 2; i < hfh - 2; ++i) {
                    //from left
                    if (hilite[3][i][j] > epsilon) {
                        hilite_dir0[c][j][i] = hilite[c][i][j] / hilite[3][i][j];
                    } else {
                        hilite_dir0[c][j][i] = 0.1f * ((hilite_dir0[0 + c][j - 1][i - 2] + hilite_dir0[0 + c][j - 1][i - 1] + hilite_dir0[0 + c][j - 1][i] + hilite_dir0[0 + c][j - 1][i + 1] + hilite_dir0[0 + c][j - 1][i + 2]) /
                                                       (hilite_dir0[0 + 3][j - 1][i - 2] + hilite_dir0[0 + 3][j - 1][i - 1] + hilite_dir0[0 + 3][j - 1][i] + hilite_dir0[0 + 3][j - 1][i + 1] + hilite_dir0[0 + 3][j - 1][i + 2] + epsilon));
                    }
                }

                if (hilite[3][2][j] <= epsilon) {
                    hilite_dir[0 + c][0][j] = hilite_dir0[c][j][2];
                }

                if (hilite[3][3][j] <= epsilon) {
                    hilite_dir[0 + c][1][j] = hilite_dir0[c][j][3];
                }

                if (hilite[3][hfh - 3][j] <= epsilon) {
                    hilite_dir[4 + c][hfh - 1][j] = hilite_dir0[c][j][hfh - 3];
                }

                if (hilite[3][hfh - 4][j] <= epsilon) {
                    hilite_dir[4 + c][hfh - 2][j] = hilite_dir0[c][j][hfh - 4];
                }
            }

            for (int i = 2; i < hfh - 2; ++i) {
                if (hilite[3][i][hfw - 2] <= epsilon) {
                    hilite_dir4[c][hfw - 1][i] = hilite_dir0[c][hfw - 2][i];
                }
            }
        }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            for (int j = hfw - 2; j > 0; --j) {
                for (int i = 2; i < hfh - 2; ++i) {
                    //from right
                    if (hilite[3][i][j] > epsilon) {
                        hilite_dir4[3][j][i] = 1.f;
                    } else {
                        hilite_dir4[3][j][i] = (hilite_dir4[3][(j + 1)][(i - 2)] + hilite_dir4[3][(j + 1)][(i - 1)] + hilite_dir4[3][(j + 1)][(i)] + hilite_dir4[3][(j + 1)][(i + 1)] + hilite_dir4[3][(j + 1)][(i + 2)]) == 0.f ? 0.f : 0.1f;
                    }
                }

                if (hilite[3][2][j] <= epsilon) {
                    hilite_dir[0 + 3][0][j] += hilite_dir4[3][j][2];
                }

                if (hilite[3][hfh - 3][j] <= epsilon) {
                    hilite_dir[4 + 3][hfh - 1][j] += hilite_dir4[3][j][hfh - 3];
                }
            }

            for (int i = 2; i < hfh - 2; ++i) {
                if (hilite[3][i][0] <= epsilon) {
                    hilite_dir[0 + 3][i - 2][0] += hilite_dir4[3][0][i];
                    hilite_dir[4 + 3][i + 2][0] += hilite_dir4[3][0][i];
                }

                if (hilite[3][i][1] <= epsilon) {
                    hilite_dir[0 + 3][i - 2][1] += hilite_dir4[3][1][i];
                    hilite_dir[4 + 3][i + 2][1] += hilite_dir4[3][1][i];
                }

                if (hilite[3][i][hfw - 2] <= epsilon) {
                    hilite_dir[0 + 3][i - 2][hfw - 2] += hilite_dir4[3][hfw - 2][i];
                    hilite_dir[4 + 3][i + 2][hfw - 2] += hilite_dir4[3][hfw - 2][i];
                }
            }
        }
    }
    if (plistener) {
        progress += 0.05;
        plistener->setProgress(progress);
    }

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef _OPENMP
        #pragma omp for nowait
#endif
        for (int c = 0; c < 3; ++c) {
            for (int j = hfw - 2; j > 0; --j) {
                for (int i = 2; i < hfh - 2; ++i) {
                    //from right
                    if (hilite[3][i][j] > epsilon) {
                        hilite_dir4[c][j][i] = hilite[c][i][j] / hilite[3][i][j];
                    } else {
                        hilite_dir4[c][j][i] = 0.1f * ((hilite_dir4[c][(j + 1)][(i - 2)] + hilite_dir4[c][(j + 1)][(i - 1)] + hilite_dir4[c][(j + 1)][(i)] + hilite_dir4[c][(j + 1)][(i + 1)] + hilite_dir4[c][(j + 1)][(i + 2)]) /
                                                      (hilite_dir4[3][(j + 1)][(i - 2)] + hilite_dir4[3][(j + 1)][(i - 1)] + hilite_dir4[3][(j + 1)][(i)] + hilite_dir4[3][(j + 1)][(i + 1)] + hilite_dir4[3][(j + 1)][(i + 2)] + epsilon));
                    }
                }

                if (hilite[3][2][j] <= epsilon) {
                    hilite_dir[0 + c][0][j] += hilite_dir4[c][j][2];
                }

                if (hilite[3][hfh - 3][j] <= epsilon) {
                    hilite_dir[4 + c][hfh - 1][j] += hilite_dir4[c][j][hfh - 3];
                }
            }

            for (int i = 2; i < hfh - 2; ++i) {
                if (hilite[3][i][0] <= epsilon) {
                    hilite_dir[0 + c][i - 2][0] += hilite_dir4[c][0][i];
                    hilite_dir[4 + c][i + 2][0] += hilite_dir4[c][0][i];
                }

                if (hilite[3][i][1] <= epsilon) {
                    hilite_dir[0 + c][i - 2][1] += hilite_dir4[c][1][i];
                    hilite_dir[4 + c][i + 2][1] += hilite_dir4[c][1][i];
                }

                if (hilite[3][i][hfw - 2] <= epsilon) {
                    hilite_dir[0 + c][i - 2][hfw - 2] += hilite_dir4[c][hfw - 2][i];
                    hilite_dir[4 + c][i + 2][hfw - 2] += hilite_dir4[c][hfw - 2][i];
                }
            }
        }

#ifdef _OPENMP
        #pragma omp single
#endif
        {
            for (int i = 1; i < hfh - 1; ++i)
                for (int j = 2; j < hfw - 2; ++j) {
                    //from top
                    if (hilite[3][i][j] > epsilon) {
                        hilite_dir[0 + 3][i][j] = 1.f;
                    } else {
                        hilite_dir[0 + 3][i][j] = (hilite_dir[0 + 3][i - 1][j - 2] + hilite_dir[0 + 3][i - 1][j - 1] + hilite_dir[0 + 3][i - 1][j] + hilite_dir[0 + 3][i - 1][j + 1] + hilite_dir[0 + 3][i - 1][j + 2]) == 0.f ? 0.f : 0.1f;
                    }
                }

            for (int j = 2; j < hfw - 2; ++j) {
                if (hilite[3][hfh - 2][j] <= epsilon) {
                    hilite_dir[4 + 3][hfh - 1][j] += hilite_dir[0 + 3][hfh - 2][j];
                }
            }
        }
    }
    if (plistener) {
        progress += 0.05;
        plistener->setProgress(progress);
    }

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef _OPENMP
        #pragma omp for nowait
#endif
        for (int c = 0; c < 3; ++c) {
            for (int i = 1; i < hfh - 1; ++i) {
                for (int j = 2; j < hfw - 2; ++j) {
                    //from top
                    if (hilite[3][i][j] > epsilon) {
                        hilite_dir[0 + c][i][j] = hilite[c][i][j] / hilite[3][i][j];
                    } else {
                        hilite_dir[0 + c][i][j] = 0.1f * ((hilite_dir[0 + c][i - 1][j - 2] + hilite_dir[0 + c][i - 1][j - 1] + hilite_dir[0 + c][i - 1][j] + hilite_dir[0 + c][i - 1][j + 1] + hilite_dir[0 + c][i - 1][j + 2]) /
                                                         (hilite_dir[0 + 3][i - 1][j - 2] + hilite_dir[0 + 3][i - 1][j - 1] + hilite_dir[0 + 3][i - 1][j] + hilite_dir[0 + 3][i - 1][j + 1] + hilite_dir[0 + 3][i - 1][j + 2] + epsilon));
                    }
                }
            }

            for (int j = 2; j < hfw - 2; ++j) {
                if (hilite[3][hfh - 2][j] <= epsilon) {
                    hilite_dir[4 + c][hfh - 1][j] += hilite_dir[0 + c][hfh - 2][j];
                }
            }
        }


#ifdef _OPENMP
        #pragma omp single
#endif
        for (int i = hfh - 2; i > 0; --i) {
            for (int j = 2; j < hfw - 2; ++j) {
                //from bottom
                if (hilite[3][i][j] > epsilon) {
                    hilite_dir[4 + 3][i][j] = 1.f;
                } else {
                    hilite_dir[4 + 3][i][j] = (hilite_dir[4 + 3][(i + 1)][(j - 2)] + hilite_dir[4 + 3][(i + 1)][(j - 1)] + hilite_dir[4 + 3][(i + 1)][(j)] + hilite_dir[4 + 3][(i + 1)][(j + 1)] + hilite_dir[4 + 3][(i + 1)][(j + 2)]) == 0.f ? 0.f : 0.1f;
                }
            }
        }
    }

    if (plistener) {
        progress += 0.05;
        plistener->setProgress(progress);
    }

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int c = 0; c < 4; ++c) {
        for (int i = hfh - 2; i > 0; --i) {
            for (int j = 2; j < hfw - 2; ++j) {
                //from bottom
                if (hilite[3][i][j] > epsilon) {
                    hilite_dir[4 + c][i][j] = hilite[c][i][j] / hilite[3][i][j];
                } else {
                    hilite_dir[4 + c][i][j] = 0.1f * ((hilite_dir[4 + c][(i + 1)][(j - 2)] + hilite_dir[4 + c][(i + 1)][(j - 1)] + hilite_dir[4 + c][(i + 1)][(j)] + hilite_dir[4 + c][(i + 1)][(j + 1)] + hilite_dir[4 + c][(i + 1)][(j + 2)]) /
                                                     (hilite_dir[4 + 3][(i + 1)][(j - 2)] + hilite_dir[4 + 3][(i + 1)][(j - 1)] + hilite_dir[4 + 3][(i + 1)][(j)] + hilite_dir[4 + 3][(i + 1)][(j + 1)] + hilite_dir[4 + 3][(i + 1)][(j + 2)] + epsilon));
                }
            }
        }
    }

    if (plistener) {
        progress += 0.05;
        plistener->setProgress(progress);
    }

    //fill in edges
    for (int dir = 0; dir < 2; ++dir) {
        for (int i = 1; i < hfh - 1; ++i) {
            for (int c = 0; c < 4; ++c) {
                hilite_dir[dir * 4 + c][i][0] = hilite_dir[dir * 4 + c][i][1];
                hilite_dir[dir * 4 + c][i][hfw - 1] = hilite_dir[dir * 4 + c][i][hfw - 2];
            }
        }

        for (int j = 1; j < hfw - 1; ++j) {
            for (int c = 0; c < 4; ++c) {
                hilite_dir[dir * 4 + c][0][j] = hilite_dir[dir * 4 + c][1][j];
                hilite_dir[dir * 4 + c][hfh - 1][j] = hilite_dir[dir * 4 + c][hfh - 2][j];
            }
        }

        for (int c = 0; c < 4; ++c) {
            hilite_dir[dir * 4 + c][0][0] = hilite_dir[dir * 4 + c][1][0] = hilite_dir[dir * 4 + c][0][1] = hilite_dir[dir * 4 + c][1][1] = hilite_dir[dir * 4 + c][2][2];
            hilite_dir[dir * 4 + c][0][hfw - 1] = hilite_dir[dir * 4 + c][1][hfw - 1] = hilite_dir[dir * 4 + c][0][hfw - 2] = hilite_dir[dir * 4 + c][1][hfw - 2] = hilite_dir[dir * 4 + c][2][hfw - 3];
            hilite_dir[dir * 4 + c][hfh - 1][0] = hilite_dir[dir * 4 + c][hfh - 2][0] = hilite_dir[dir * 4 + c][hfh - 1][1] = hilite_dir[dir * 4 + c][hfh - 2][1] = hilite_dir[dir * 4 + c][hfh - 3][2];
            hilite_dir[dir * 4 + c][hfh - 1][hfw - 1] = hilite_dir[dir * 4 + c][hfh - 2][hfw - 1] = hilite_dir[dir * 4 + c][hfh - 1][hfw - 2] = hilite_dir[dir * 4 + c][hfh - 2][hfw - 2] = hilite_dir[dir * 4 + c][hfh - 3][hfw - 3];
        }
    }

    for (int i = 1; i < hfh - 1; ++i) {
        for (int c = 0; c < 4; ++c) {
            hilite_dir0[c][0][i] = hilite_dir0[c][1][i];
            hilite_dir0[c][hfw - 1][i] = hilite_dir0[c][hfw - 2][i];
        }
    }

    for (int j = 1; j < hfw - 1; ++j) {
        for (int c = 0; c < 4; ++c) {
            hilite_dir0[c][j][0] = hilite_dir0[c][j][1];
            hilite_dir0[c][j][hfh - 1] = hilite_dir0[c][j][hfh - 2];
        }
    }

    for (int c = 0; c < 4; ++c) {
        hilite_dir0[c][0][0] = hilite_dir0[c][0][1] = hilite_dir0[c][1][0] = hilite_dir0[c][1][1] = hilite_dir0[c][2][2];
        hilite_dir0[c][hfw - 1][0] = hilite_dir0[c][hfw - 1][1] = hilite_dir0[c][hfw - 2][0] = hilite_dir0[c][hfw - 2][1] = hilite_dir0[c][hfw - 3][2];
        hilite_dir0[c][0][hfh - 1] = hilite_dir0[c][0][hfh - 2] = hilite_dir0[c][1][hfh - 1] = hilite_dir0[c][1][hfh - 2] = hilite_dir0[c][2][hfh - 3];
        hilite_dir0[c][hfw - 1][hfh - 1] = hilite_dir0[c][hfw - 1][hfh - 2] = hilite_dir0[c][hfw - 2][hfh - 1] = hilite_dir0[c][hfw - 2][hfh - 2] = hilite_dir0[c][hfw - 3][hfh - 3];
    }

    for (int i = 1; i < hfh - 1; ++i) {
        for (int c = 0; c < 4; ++c) {
            hilite_dir4[c][0][i] = hilite_dir4[c][1][i];
            hilite_dir4[c][hfw - 1][i] = hilite_dir4[c][hfw - 2][i];
        }
    }

    for (int j = 1; j < hfw - 1; ++j) {
        for (int c = 0; c < 4; ++c) {
            hilite_dir4[c][j][0] = hilite_dir4[c][j][1];
            hilite_dir4[c][j][hfh - 1] = hilite_dir4[c][j][hfh - 2];
        }
    }

    for (int c = 0; c < 4; ++c) {
        hilite_dir4[c][0][0] = hilite_dir4[c][0][1] = hilite_dir4[c][1][0] = hilite_dir4[c][1][1] = hilite_dir4[c][2][2];
        hilite_dir4[c][hfw - 1][0] = hilite_dir4[c][hfw - 1][1] = hilite_dir4[c][hfw - 2][0] = hilite_dir4[c][hfw - 2][1] = hilite_dir4[c][hfw - 3][2];
        hilite_dir4[c][0][hfh - 1] = hilite_dir4[c][0][hfh - 2] = hilite_dir4[c][1][hfh - 1] = hilite_dir4[c][1][hfh - 2] = hilite_dir4[c][2][hfh - 3];
        hilite_dir4[c][hfw - 1][hfh - 1] = hilite_dir4[c][hfw - 1][hfh - 2] = hilite_dir4[c][hfw - 2][hfh - 1] = hilite_dir4[c][hfw - 2][hfh - 2] = hilite_dir4[c][hfw - 3][hfh - 3];
    }

    if (plistener) {
        progress += 0.05;
        plistener->setProgress(progress);
    }

    //free up some memory
    for (int c = 0; c < 4; ++c) {
        hilite[c].free();
    }

    // now reconstruct clipped channels using color ratios
    //using code from ART - thanks to Alberto Griggio
    const int W2 = blur > 0 ? blurWidth / 2.f + 0.5f : 0;
    const int H2 = blur > 0 ? blurHeight / 2.f + 0.5f : 0;
    array2D<float> mask(W2, H2, ARRAY2D_CLEAR_DATA);
    array2D<float> rbuf(W2, H2);
    array2D<float> gbuf(W2, H2);
    array2D<float> bbuf(W2, H2);
    array2D<float> guide(W2, H2);
   
    if (blur > 0) {
        array2D<float> rbuffer(blurWidth, blurHeight, minx, miny, red, ARRAY2D_BYREFERENCE);
        rescaleNearest(rbuffer, rbuf, true);
        array2D<float> gbuffer(blurWidth, blurHeight, minx, miny, green, ARRAY2D_BYREFERENCE);
        rescaleNearest(gbuffer, gbuf, true);
        array2D<float> bbuffer(blurWidth, blurHeight, minx, miny, blue, ARRAY2D_BYREFERENCE);
        rescaleNearest(bbuffer, bbuf, true);

        LUTf gamma(65536);
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int i = 0; i < 65536; ++i) {
            gamma[i] = pow_F(i / 65535.f, 2.2f);
        }

        const float xyzcam[3] = {static_cast<float>(imatrices.xyz_cam[1][0]), static_cast<float>(imatrices.xyz_cam[1][1]), static_cast<float>(imatrices.xyz_cam[1][2])};
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int y = 0; y < H2; ++y) {
            for (int x = 0; x < W2; ++x) {
                guide[y][x] = gamma[Color::rgbLuminance(rbuf[y][x], gbuf[y][x], bbuf[y][x], xyzcam)];
            }
        }
    }
//end adding code ART

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif
    for (int i = 0; i < blurHeight; ++i) {
        const int i1 = min((i - i % pitch) / pitch, hfh - 1);

        for (int j = 0; j < blurWidth; ++j) {
            const float pixel[3] = {
                red[i + miny][j + minx],
                green[i + miny][j + minx],
                blue[i + miny][j + minx]
            };

            if (pixel[0] < max_f[0] && pixel[1] < max_f[1] && pixel[2] < max_f[2]) {
                continue;    //pixel not clipped
            }

            const int j1 = min((j - j % pitch) / pitch, hfw - 1);

            //estimate recovered values using modified HLRecovery_blend algorithm
            float rgb[3] = {
                pixel[0],
                pixel[1],
                pixel[2]
            };// Copy input pixel to rgb so it's easier to access in loops
            float rgb_blend[3] = {};
            float cam[2][3];
            float lab[2][3];
            float sum[2];

            // Initialize cam with raw input [0] and potentially clipped input [1]
            for (int c = 0; c < 3; ++c) {
                cam[0][c] = rgb[c];
                cam[1][c] = min(cam[0][c], clippt);
            }

            // Calculate the lightness correction ratio (chratio)
            for (int i2 = 0; i2 < 2; ++i2) {
                for (int c = 0; c < 3; ++c) {
                    lab[i2][c] = 0;

                    for (int j2 = 0; j2 < 3; ++j2) {
                        lab[i2][c] += trans[c][j2] * cam[i2][j2];
                    }
                }

                sum[i2] = 0.f;

                for (int c = 1; c < 3; ++c) {
                    sum[i2] += SQR(lab[i2][c]);
                }
            }

            // avoid division by zero
            sum[0] = std::max(sum[0], epsilon);

            const float chratio = sqrtf(sum[1] / sum[0]);

            // Apply ratio to lightness in lab space
            for (int c = 1; c < 3; ++c) {
                lab[0][c] *= chratio;
            }

            // Transform back from lab to RGB
            for (int c = 0; c < 3; ++c) {
                cam[0][c] = 0.f;

                for (int j2 = 0; j2 < 3; ++j2) {
                    cam[0][c] += itrans[c][j2] * lab[0][j2];
                }
            }

            for (int c = 0; c < 3; ++c) {
                rgb[c] = cam[0][c] / 3;
            }

            // Copy converted pixel back
            if (pixel[0] > blendpt) {
                const float rfrac = LIM01(medFactor[0] * (pixel[0] - blendpt));
                rgb_blend[0] = intp(rfrac, rgb[0], pixel[0]);
            }

            if (pixel[1] > blendpt) {
                const float gfrac = LIM01(medFactor[1] * (pixel[1] - blendpt));
                rgb_blend[1] = intp(gfrac, rgb[1], pixel[1]);
            }

            if (pixel[2] > blendpt) {
                const float bfrac = LIM01(medFactor[2] * (pixel[2] - blendpt));
                rgb_blend[2] = intp(bfrac, rgb[2], pixel[2]);
            }

            //end of HLRecovery_blend estimation

            //there are clipped highlights
            //first, determine weighted average of unclipped extensions (weighting is by 'hue' proximity)
            bool totwt = false;
            float clipfix[3] = {0.f, 0.f, 0.f};

            float Y = epsilon + rgb_blend[0] + rgb_blend[1] + rgb_blend[2];

            for (int c = 0; c < 3; ++c) {
                rgb_blend[c] /= Y;
            }

            float Yhi = 1.f / (hilite_dir0[0][j1][i1] + hilite_dir0[1][j1][i1] + hilite_dir0[2][j1][i1]);

            if (Yhi < 2.f) {
                const float dirwt = 1.f / ((1.f + 65535.f * (SQR(rgb_blend[0] - hilite_dir0[0][j1][i1] * Yhi) +
                                                      SQR(rgb_blend[1] - hilite_dir0[1][j1][i1] * Yhi) +
                                                      SQR(rgb_blend[2] - hilite_dir0[2][j1][i1] * Yhi))) * (hilite_dir0[3][j1][i1] + epsilon));
                totwt = true;
                clipfix[0] = dirwt * hilite_dir0[0][j1][i1];
                clipfix[1] = dirwt * hilite_dir0[1][j1][i1];
                clipfix[2] = dirwt * hilite_dir0[2][j1][i1];
            }

            for (int dir = 0; dir < 2; ++dir) {
                const float Yhi2 = 1.f / (hilite_dir[dir * 4 + 0][i1][j1] + hilite_dir[dir * 4 + 1][i1][j1] + hilite_dir[dir * 4 + 2][i1][j1]);

                if (Yhi2 < 2.f) {
                    const float dirwt = 1.f / ((1.f + 65535.f * (SQR(rgb_blend[0] - hilite_dir[dir * 4 + 0][i1][j1] * Yhi2) +
                                                          SQR(rgb_blend[1] - hilite_dir[dir * 4 + 1][i1][j1] * Yhi2) +
                                                          SQR(rgb_blend[2] - hilite_dir[dir * 4 + 2][i1][j1] * Yhi2))) * (hilite_dir[dir * 4 + 3][i1][j1] + epsilon));
                    totwt = true;
                    clipfix[0] += dirwt * hilite_dir[dir * 4 + 0][i1][j1];
                    clipfix[1] += dirwt * hilite_dir[dir * 4 + 1][i1][j1];
                    clipfix[2] += dirwt * hilite_dir[dir * 4 + 2][i1][j1];
                }
            }


            Yhi = 1.f / (hilite_dir4[0][j1][i1] + hilite_dir4[1][j1][i1] + hilite_dir4[2][j1][i1]);

            if (Yhi < 2.f) {
                const float dirwt = 1.f / ((1.f + 65535.f * (SQR(rgb_blend[0] - hilite_dir4[0][j1][i1] * Yhi) +
                                                      SQR(rgb_blend[1] - hilite_dir4[1][j1][i1] * Yhi) +
                                                      SQR(rgb_blend[2] - hilite_dir4[2][j1][i1] * Yhi))) * (hilite_dir4[3][j1][i1] + epsilon));
                totwt = true;
                clipfix[0] += dirwt * hilite_dir4[0][j1][i1];
                clipfix[1] += dirwt * hilite_dir4[1][j1][i1];
                clipfix[2] += dirwt * hilite_dir4[2][j1][i1];
            }

            if (UNLIKELY(!totwt)) {
                continue;
            }

           //using code from ART - thanks to Alberto Griggio
            float maskval = 1.f;
            const int yy = i + miny;
            const int xx = j + minx;

            //now correct clipped channels
            if (pixel[0] > max_f[0] && pixel[1] > max_f[1] && pixel[2] > max_f[2]) {
                //all channels clipped
                const float mult = whitept / (0.299f * clipfix[0] + 0.587f * clipfix[1] + 0.114f * clipfix[2]);
                red[yy][xx]   = clipfix[0] * mult;
                green[yy][xx] = clipfix[1] * mult;
                blue[yy][xx]  = clipfix[2] * mult;
            } else {//some channels clipped
                const float notclipped[3] = {
                    pixel[0] <= max_f[0] ? 1.f : 0.f,
                    pixel[1] <= max_f[1] ? 1.f : 0.f,
                    pixel[2] <= max_f[2] ? 1.f : 0.f
                };

                if (notclipped[0] == 0.f) { //red clipped
                    red[yy][xx]  = max(pixel[0], clipfix[0] * ((notclipped[1] * pixel[1] + notclipped[2] * pixel[2]) /
                                                 (notclipped[1] * clipfix[1] + notclipped[2] * clipfix[2] + epsilon)));
                }

                if (notclipped[1] == 0.f) { //green clipped
                    green[yy][xx] = max(pixel[1], clipfix[1] * ((notclipped[2] * pixel[2] + notclipped[0] * pixel[0]) /
                                                    (notclipped[2] * clipfix[2] + notclipped[0] * clipfix[0] + epsilon)));
                }

                if (notclipped[2] == 0.f) { //blue clipped
                    blue[yy][xx]  = max(pixel[2], clipfix[2] * ((notclipped[0] * pixel[0] + notclipped[1] * pixel[1]) /
                                                   (notclipped[0] * clipfix[0] + notclipped[1] * clipfix[1] + epsilon)));
                }

                maskval = 1.f - (notclipped[0] + notclipped[1] + notclipped[2]) / 5.f;
            }

            Y = 0.299f * red[yy][xx] + 0.587f * green[yy][xx] + 0.114f * blue[yy][xx];

            if (Y > whitept) {
                const float mult = whitept / Y;
                red[yy][xx] *= mult;
                green[yy][xx] *= mult;
                blue[yy][xx] *= mult;
            }

            if (blur > 0) {
                const int ii = i / 2;
                const int jj = j / 2;
                rbuf[ii][jj] = red[yy][xx];
                gbuf[ii][jj] = green[yy][xx];
                bbuf[ii][jj] = blue[yy][xx];
                mask[ii][jj] = maskval;
            }
        }
    }

    if (blur > 0) {
        if (plistener) {
            progress += 0.05;
            plistener->setProgress(progress);
        }
        blur = rtengine::LIM(blur - 1, 0, 3);

        constexpr float vals[4][3] = {{4.0f, 0.3f, 0.3f},
                                  //    {3.5f, 0.5f, 0.2f},
                                      {3.0f, 1.0f, 0.1f},
                                      {3.0f, 2.0f, 0.01f},
                                      {2.0f, 3.0f, 0.001f}
                                     };

        const float radius1 = vals[blur][0];
        const float radius2 = vals[blur][1];
        const float th = vals[blur][2];

        guidedFilter(guide, mask, mask, radius1, th, true, 1);
        if (plistener) {
            progress += 0.03;
            plistener->setProgress(progress);
        }
        if (blur > 0) { //no use of 2nd guidedFilter if Blur = 0 (slider to 1)..speed-up and very small differences.
            guidedFilter(guide, rbuf, rbuf, radius2, 0.01f * 65535.f, true, 1);
            if (plistener) {
                progress += 0.03;
                plistener->setProgress(progress);
            }
            guidedFilter(guide, gbuf, gbuf, radius2, 0.01f * 65535.f, true, 1);
            if (plistener) {
                progress += 0.03;
                plistener->setProgress(progress);
            }
            guidedFilter(guide, bbuf, bbuf, radius2, 0.01f * 65535.f, true, 1);
            if (plistener) {
                progress += 0.03;
                plistener->setProgress(progress);
            }
        }
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif
        for (int y = 0; y < blurHeight; ++y) {
            const float fy = y * 0.5f;
            const int yy = y / 2;
            for (int x = 0; x < blurWidth; ++x) {
                const int xx = x / 2;
                const float m = mask[yy][xx];
                if (m > 0.f) {
                    const float fx = x * 0.5f;
                    red[y + miny][x + minx] = intp(m, getBilinearValue(rbuf, fx, fy), red[y + miny][x + minx]);
                    green[y + miny][x + minx] = intp(m, getBilinearValue(gbuf, fx, fy), green[y + miny][x + minx]);
                    blue[y + miny][x + minx] = intp(m, getBilinearValue(bbuf, fx, fy), blue[y + miny][x + minx]);
                }
            }
        }
    }

    if (plistener) {
        plistener->setProgress(1.00);
    }

}// end of HLReconstruction


//-----------------------------------------------------------------------------
// "inpaint opposed" algorithm taken from darktable
//
// (Very effective, very simple, very neat)
//
// Kudos to the original authors (@jenshannoschwalm from dt, in collaboration
// with @garagecoder and @Iain from gmic).
//
// Copyright and description of the original code follows
//
/*
    Copyright (C) 2022 darktable developers.

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/

/* The refavg values are calculated in raw-RGB-cube3 space
   We calculate all color channels in the 3x3 photosite area, this can be understaood as a "superpixel",
   the "asking" location is in the centre.
   As this works for bayer and xtrans sensors we don't have a fixed ratio but calculate the average
   for every color channel first.
   refavg for one of red, green or blue is defined as means of both other color channels (opposing).
   
   The basic idea / observation for the _process_opposed algorithm is, the refavg is a good estimate
   for any clipped color channel in the vast majority of images, working mostly fine both for small specular
   highlighted spots and large areas.
   
   The correction via some sort of global chrominance further helps to correct color casts.
   The chrominace data are taken from the areas morphologically very close to clipped data.
   Failures of the algorithm (color casts) are in most cases related to
    a) very large differences between optimal white balance coefficients vs what we have as D65 in the darktable pipeline
    b) complicated lightings so the gradients are not well related
    c) a wrong whitepoint setting in the rawprepare module. 
    d) the maths might not be best
*/ 
//-----------------------------------------------------------------------------

namespace {

constexpr int HL_BORDER = 8;
constexpr float HL_POWERF = 3.0f;

// void border_fill_zero(int *d, int width, int height)
// {
//     for (int i = 0; i < HL_BORDER * width; i++) {
//         d[i] = 0;
//     }
//     for (int i = (height - HL_BORDER - 1) * width; i < width*height; i++) {
//         d[i] = 0;
//     }
//     for (int row = HL_BORDER; row < height - HL_BORDER; row++) {
//         int *p1 = d + row*width;
//         int *p2 = d + (row+1)*width - HL_BORDER;
//         for(int i = 0; i < HL_BORDER; i++) {
//             p1[i] = p2[i] = 0;
//         }
//     }
// }


int test_dilate(const int *img, int i, int w1)
{
    int retval = 0;
    retval = img[i-w1-1] | img[i-w1] | img[i-w1+1] |
        img[i-1]    | img[i]    | img[i+1] |
        img[i+w1-1] | img[i+w1] | img[i+w1+1];
    if (retval) {
        return retval;
    }

    const size_t w2 = 2*w1;
    retval = img[i-w2-1] | img[i-w2]   | img[i-w2+1] |
        img[i-w1-2] | img[i-w1+2] | 
        img[i-2]    | img[i+2] |
        img[i+w1-2] | img[i+w1+2] |
        img[i+w2-1] | img[i+w2]   | img[i+w2+1];
    if (retval) {
        return retval;
    }

    const size_t w3 = 3*w1;
    retval = img[i-w3-2] | img[i-w3-1] | img[i-w3] | img[i-w3+1] | img[i-w3+2] |
        img[i-w2-3] | img[i-w2-2] | img[i-w2+2] | img[i-w2+3] |
        img[i-w1-3] | img[i-w1+3] | 
        img[i-3]    | img[i+3]    | 
        img[i+w1-3] | img[i+w1+3] | 
        img[i+w2-3] | img[i+w2-2] | img[i+w2+2] | img[i+w2+3] |
        img[i+w3-2] | img[i+w3-1] | img[i+w3] | img[i+w3+1] | img[i+w3+2]; 
    return retval;
}


void dilating(const int *img, int *o, int w1, int height)
{
#ifdef _OPENMP
#   pragma omp parallel for
#endif
    for (int row = HL_BORDER; row < height - HL_BORDER; row++) {
        for (int col = HL_BORDER, i = row*w1 + col; col < w1 - HL_BORDER; col++, i++) {
            o[i] = test_dilate(img, i, w1);
        }
    }
}

} // namespace

void RawImageSource::highlight_recovery_opposed(float scale_mul[3], const ColorTemp &wb, float gainth)
{
    //BENCHFUN

    if (settings->verbose) {
        std::cout << "Applying Highlight Recovery: Inpaint opposed" << std::endl;
    }

    if (plistener) {
        plistener->setProgressStr("PROGRESSBAR_HLREC");
        plistener->setProgress(0);
    }
    
    double rr, gg, bb;
    wb.getMultipliers(rr, gg, bb);
    wbMul2Camera(rr, gg, bb);
    
	float gain = 1.2f * gainth;
	
	float clipval = 0.987f / gain;
    const float scalecoeffs[3] = {
        scale_mul[0] * float(rr) / 65535.f,
        scale_mul[1] * float(gg) / 65535.f,
        scale_mul[2] * float(bb) / 65535.f,
    };
    const float clips[3] = {
        clipval * float(rr),
        clipval * float(gg),
        clipval * float(bb)
    };
    const float clipdark[3] = {
        0.03f * clips[0],
        0.125f * clips[1],
        0.03f * clips[2]
    };

    bool anyclipped = false;
    float **chan[3] = { red, green, blue };

    const float clipscale[3] = {
        clips[0] / scalecoeffs[0],
        clips[1] / scalecoeffs[1],
        clips[2] / scalecoeffs[2]
    };

    int x1 = W, y1 = H, x2 = 0, y2 = 0;
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            for (int c = 0; c < 3; ++c) {
                if (chan[c][y][x] >= clipscale[c]) {
                    anyclipped = true;
                    x1 = std::min(x, x1);
                    x2 = std::max(x, x2);
                    y1 = std::min(y, y1);
                    y2 = std::max(y, y2);
                }
            }
        }
    }

    if (!anyclipped) {
        if (plistener) {
            plistener->setProgress(1.0);
        }
        return;
    }

    x1 = std::max(x1-1, 0);
    x2 = std::min(x2+1, W-1);
    y1 = std::max(y1-1, 0);
    y2 = std::min(y2+1, H-1);

    const int cW = x2 - x1 + 1;
    const int cH = y2 - y1 + 1;
    
#ifdef _OPENMP
#   pragma omp parallel for
#endif
    for (int y = 0; y < cH; ++y) {
        const int yy = y + y1;
        for (int x = 0; x < cW; ++x) {
            const int xx = x + x1;
            for (int c = 0; c < 3; ++c) {
                chan[c][yy][xx] *= scalecoeffs[c];
            }
        }
    }

    if (plistener) {
        plistener->setProgress(0.1);
    }

    multi_array2D<float, 3> tmp(cW, cH);

    const int pwidth = cW + 2 * HL_BORDER;
    const int pheight = cH + 2 * HL_BORDER;
    const int p_size = pwidth * pheight;
    AlignedBuffer<int> mask_vec(4 * p_size);
    int *mask_buffer = mask_vec.data;

    const auto mask_val =
        [&](int c, int y, int x) -> int &
        {
            return mask_buffer[c * p_size + (HL_BORDER + y) * pwidth + x + HL_BORDER];
        };
    
    const auto set_refavg =
        [&](int y, int x) -> bool
        {
            const int yy = y + y1;
            const int xx = x + x1;
            bool found = false;
            for (int c = 0; c < 3 && !found; ++c) {
                if (chan[c][yy][xx] >= clips[c]) {
                    found = true;
                }
            }
            if (!found) {
                return false;
            }
            
            float mean[3] = { 0.0f, 0.0f, 0.0f };
            for (int dy = -1; dy < 2; dy++) {
                for (int dx = -1; dx < 2; dx++) {
                    for (int c = 0; c < 3; ++c) {
                        mean[c] += std::max(0.0f, chan[c][yy+dy][xx+dx]);
                    }
                }
            }
            for (int c = 0; c < 3; ++c) {
                mean[c] = pow_F(mean[c] / 9.0f, 1.0f / HL_POWERF);
            }

            const float croot_refavg[3] = {
                0.5f * (mean[1] + mean[2]),
                0.5f * (mean[0] + mean[2]),
                0.5f * (mean[0] + mean[1])
            };
            
            for (int c = 0; c < 3; ++c) {
                if (chan[c][yy][xx] >= clips[c]) {
                    tmp[c][y][x] = pow_F(croot_refavg[c], HL_POWERF);
                    mask_val(c, y, x) = 1;
                }
            }
            return true;
        };

#ifdef _OPENMP
#   pragma omp parallel for
#endif
    for (int y = 0; y < cH; ++y) {
        const int yy = y + y1;
        for (int x = 0; x < cW; ++x) {
            const int xx = x + x1;
            for (int c = 0; c < 3; ++c) {
                tmp[c][y][x] = std::max(0.f, chan[c][yy][xx]);
            }

            if ((x > 0) && (x < cW - 1) && (y > 0) && (y < cH - 1)) {
                set_refavg(y, x);
            }
        }
    }

    if (plistener) {
        plistener->setProgress(0.3);
    }

    for (size_t i = 0; i < 3; i++) {
        int *mask = mask_buffer + i * p_size;
        int *tmp = mask_buffer + 3 * p_size;
        //border_fill_zero(mask, pwidth, pheight);
        dilating(mask, tmp, pwidth, pheight);
        memcpy(mask, tmp, p_size * sizeof(int));
    }

    float cr_sum[3] = { 0.f, 0.f, 0.f };
    int cr_cnt[3] = { 0, 0, 0 };

#ifdef _OPENMP
#   pragma omp parallel for reduction(+ : cr_sum, cr_cnt)
#endif
    for (int y = 1; y < cH-1; ++y) {
        const int yy = y + y1;
        for (int x = 1; x < cW-1; ++x) {
            const int xx = x + x1;
            for (int c = 0; c < 3; ++c) {
                const float inval = std::max(0.0f, chan[c][yy][xx]); 
                if (mask_val(c, y, x) && (inval > clipdark[c]) && (inval < clips[c])) {
                    cr_sum[c] += inval - tmp[c][y][x];
                    ++cr_cnt[c];
                }
            }
        }
    }

    if (plistener) {
        plistener->setProgress(0.6);
    }

    float chrominance[3] = {
        cr_sum[0] / std::max(1.f, float(cr_cnt[0])),
        cr_sum[1] / std::max(1.f, float(cr_cnt[1])),
        cr_sum[2] / std::max(1.f, float(cr_cnt[2]))
    };

#ifdef _OPENMP
#   pragma omp parallel for 
#endif
    for (int y = 0; y < cH; ++y) {
        const int yy = y + y1;
        for (int x = 0; x < cW; ++x) {
            const int xx = x + x1;
            for (int c = 0; c < 3; ++c) {
                const float inval = std::max(0.0f, chan[c][yy][xx]);
                if (inval >= clips[c]) {
                    chan[c][yy][xx] = std::max(inval, tmp[c][y][x] + chrominance[c]);
                }
            }
        }
    }

    if (plistener) {
        plistener->setProgress(0.9);
    }

#ifdef _OPENMP
#   pragma omp parallel for
#endif
    for (int y = 0; y < cH; ++y) {
        const int yy = y + y1;
        for (int x = 0; x < cW; ++x) {
            const int xx = x + x1;
            for (int c = 0; c < 3; ++c) {
                chan[c][yy][xx] /= scalecoeffs[c];
            }
        }
    }

    if (plistener) {
        plistener->setProgress(1.0);
    }
}




}

