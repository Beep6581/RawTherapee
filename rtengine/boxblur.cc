/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (C) 2010 Emil Martinec <ejmartin@uchicago.edu>
 *  Copyright (C) 2019 Ingo Weyrich <heckflosse67@gmx.de>
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

#include <memory>
#include <cmath>

#include "boxblur.h"

#include "rt_math.h"
#include "opthelper.h"

namespace rtengine
{

void boxblur(float** src, float** dst, int radius, int W, int H, bool multiThread)
{
    //box blur using rowbuffers and linebuffers instead of a full size buffer

    radius = rtengine::min(radius, W - 1, H - 1);
    if (radius == 0) {
        if (src != dst) {
#ifdef _OPENMP
            #pragma omp parallel for if (multiThread)
#endif

            for (int row = 0; row < H; ++row) {
                for (int col = 0; col < W; ++col) {
                    dst[row][col] = src[row][col];
                }
            }
        }
        return;
    }

    constexpr int numCols = 8; // process numCols columns at once for better usage of L1 cpu cache
#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        std::unique_ptr<float[]> buffer(new float[numCols * (radius + 1)]);

        //horizontal blur
        float* const lineBuffer = buffer.get();
#ifdef _OPENMP
        #pragma omp for
#endif
        for (int row = 0; row < H; ++row) {
            float len = radius + 1;
            float tempval = src[row][0];
            lineBuffer[0] = tempval;
            for (int j = 1; j <= radius; j++) {
                tempval += src[row][j];
            }

            tempval /= len;
            dst[row][0] = tempval;

            for (int col = 1; col <= radius; ++col) {
                lineBuffer[col] = src[row][col];
                tempval = (tempval * len + src[row][col + radius]) / (len + 1);
                dst[row][col] = tempval;
                ++len;
            }
            int pos = 0;
            for (int col = radius + 1; col < W - radius; ++col) {
                const float oldVal = lineBuffer[pos];
                lineBuffer[pos] = src[row][col];
                tempval = tempval + (src[row][col + radius] - oldVal) / len;
                dst[row][col] = tempval;
                ++pos;
                pos = pos <= radius ? pos : 0;
            }

            for (int col = W - radius; col < W; ++col) {
                tempval = (tempval * len - lineBuffer[pos]) / (len - 1);
                dst[row][col] = tempval;
                --len;
                ++pos;
                pos = pos <= radius ? pos : 0;
            }
        }

        //vertical blur
#ifdef __SSE2__
        vfloat (* const rowBuffer)[2] = (vfloat(*)[2]) buffer.get();
        const vfloat leninitv = F2V(radius + 1);
        const vfloat onev = F2V(1.f);
        vfloat tempv, temp1v, lenv, lenp1v, lenm1v, rlenv;

#ifdef _OPENMP
        #pragma omp for nowait
#endif

        for (int col = 0; col < W - 7; col += 8) {
            lenv = leninitv;
            tempv = LVFU(dst[0][col]);
            temp1v = LVFU(dst[0][col + 4]);
            rowBuffer[0][0] = tempv;
            rowBuffer[0][1] = temp1v;

            for (int i = 1; i <= radius; ++i) {
                tempv = tempv + LVFU(dst[i][col]);
                temp1v = temp1v + LVFU(dst[i][col + 4]);
            }

            tempv = tempv / lenv;
            temp1v = temp1v / lenv;
            STVFU(dst[0][col], tempv);
            STVFU(dst[0][col + 4], temp1v);

            for (int row = 1; row <= radius; ++row) {
                rowBuffer[row][0] = LVFU(dst[row][col]);
                rowBuffer[row][1] = LVFU(dst[row][col + 4]);
                lenp1v = lenv + onev;
                tempv = (tempv * lenv + LVFU(dst[row + radius][col])) / lenp1v;
                temp1v = (temp1v * lenv + LVFU(dst[row + radius][col + 4])) / lenp1v;
                STVFU(dst[row][col], tempv);
                STVFU(dst[row][col + 4], temp1v);
                lenv = lenp1v;
            }

            rlenv = onev / lenv;
            int pos = 0;
            for (int row = radius + 1; row < H - radius; ++row) {
                vfloat oldVal0 = rowBuffer[pos][0];
                vfloat oldVal1 = rowBuffer[pos][1];
                rowBuffer[pos][0] = LVFU(dst[row][col]);
                rowBuffer[pos][1] = LVFU(dst[row][col + 4]);
                tempv = tempv + (LVFU(dst[row + radius][col]) - oldVal0) * rlenv ;
                temp1v = temp1v + (LVFU(dst[row + radius][col + 4]) - oldVal1) * rlenv ;
                STVFU(dst[row][col], tempv);
                STVFU(dst[row][col + 4], temp1v);
                ++pos;
                pos = pos <= radius ? pos : 0;
            }

            for (int row = H - radius; row < H; ++row) {
                lenm1v = lenv - onev;
                tempv = (tempv * lenv - rowBuffer[pos][0]) / lenm1v;
                temp1v = (temp1v * lenv - rowBuffer[pos][1]) / lenm1v;
                STVFU(dst[row][col], tempv);
                STVFU(dst[row][col + 4], temp1v);
                lenv = lenm1v;
                ++pos;
                pos = pos <= radius ? pos : 0;
            }
        }

#else
        float (* const rowBuffer)[8] = (float(*)[8]) buffer.get();
#ifdef _OPENMP
        #pragma omp for nowait
#endif

        for (int col = 0; col < W - numCols + 1; col += 8) {
            float len = radius + 1;

            for (int k = 0; k < numCols; ++k) {
                rowBuffer[0][k] = dst[0][col + k];
            }

            for (int i = 1; i <= radius; ++i) {
                for (int k = 0; k < numCols; ++k) {
                    dst[0][col + k] += dst[i][col + k];
                }
            }

            for(int k = 0; k < numCols; ++k) {
                dst[0][col + k] /= len;
            }

            for (int row = 1; row <= radius; ++row) {
                for(int k = 0; k < numCols; ++k) {
                    rowBuffer[row][k] = dst[row][col + k];
                    dst[row][col + k] = (dst[row - 1][col + k] * len + dst[row + radius][col + k]) / (len + 1);
                }

                len ++;
            }

            int pos = 0;
            for (int row = radius + 1; row < H - radius; ++row) {
                for(int k = 0; k < numCols; ++k) {
                    float oldVal = rowBuffer[pos][k];
                    rowBuffer[pos][k] = dst[row][col + k];
                    dst[row][col + k] = dst[row - 1][col + k] + (dst[row + radius][col + k] - oldVal) / len;
                }
                ++pos;
                pos = pos <= radius ? pos : 0;
            }

            for (int row = H - radius; row < H; ++row) {
                for(int k = 0; k < numCols; ++k) {
                    dst[row][col + k] = (dst[row - 1][col + k] * len - rowBuffer[pos][k]) / (len - 1);
                }
                len --;
                ++pos;
                pos = pos <= radius ? pos : 0;
            }
        }

#endif
        //vertical blur, remaining columns
#ifdef _OPENMP
        #pragma omp single
#endif
        {
            const int remaining = W % numCols;

            if (remaining > 0) {
                float (* const rowBuffer)[8] = (float(*)[8]) buffer.get();
                const int col = W - remaining;

                float len = radius + 1;
                for(int k = 0; k < remaining; ++k) {
                    rowBuffer[0][k] = dst[0][col + k];
                }
                for (int row = 1; row <= radius; ++row) {
                    for(int k = 0; k < remaining; ++k) {
                        dst[0][col + k] += dst[row][col + k];
                    }
                }
                for(int k = 0; k < remaining; ++k) {
                    dst[0][col + k] /= len;
                }
                for (int row = 1; row <= radius; ++row) {
                    for(int k = 0; k < remaining; ++k) {
                        rowBuffer[row][k] = dst[row][col + k];
                        dst[row][col + k] = (dst[row - 1][col + k] * len + dst[row + radius][col + k]) / (len + 1);
                    }
                    len ++;
                }
                const float rlen = 1.f / len;
                int pos = 0;
                for (int row = radius + 1; row < H - radius; ++row) {
                    for(int k = 0; k < remaining; ++k) {
                        float oldVal = rowBuffer[pos][k];
                        rowBuffer[pos][k] = dst[row][col + k];
                        dst[row][col + k] = dst[row - 1][col + k] + (dst[row + radius][col + k] - oldVal) * rlen;
                    }
                    ++pos;
                    pos = pos <= radius ? pos : 0;
                }
                for (int row = H - radius; row < H; ++row) {
                    for(int k = 0; k < remaining; ++k) {
                        dst[row][col + k] = (dst[(row - 1)][col + k] * len - rowBuffer[pos][k]) / (len - 1);
                    }
                    len --;
                    ++pos;
                    pos = pos <= radius ? pos : 0;
                }
            }
        }
    }
}

void boxabsblur(float** src, float** dst, int radius, int W, int H, bool multiThread)
{
    //abs box blur using rowbuffers and linebuffers instead of a full size buffer, W should be a multiple of 16

    if (radius == 0) {
        if (src != dst) {
#ifdef _OPENMP
            #pragma omp parallel for if (multiThread)
#endif

            for (int row = 0; row < H; ++row) {
                for (int col = 0; col < W; ++col) {
                    dst[row][col] = std::fabs(src[row][col]);
                }
            }
        }
        return;
    }

    constexpr int numCols = 16; // process numCols columns at once for better usage of L1 cpu cache
#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        float buffer[numCols * (radius + 1)] ALIGNED64;

        //horizontal blur
        float* const lineBuffer = buffer;
#ifdef _OPENMP
        #pragma omp for
#endif
        for (int row = 0; row < H; ++row) {
            float len = radius + 1;
            float tempval = std::fabs(src[row][0]);
            lineBuffer[0] = tempval;
            for (int j = 1; j <= radius; j++) {
                tempval += std::fabs(src[row][j]);
            }

            tempval /= len;
            dst[row][0] = tempval;

            for (int col = 1; col <= radius; ++col) {
                lineBuffer[col] = std::fabs(src[row][col]);
                tempval = (tempval * len + std::fabs(src[row][col + radius])) / (len + 1);
                dst[row][col] = tempval;
                ++len;
            }

            const float rlen = 1.f / len;
            int pos = 0;
            for (int col = radius + 1; col < W - radius; ++col) {
                const float oldVal = lineBuffer[pos];
                lineBuffer[pos] = std::fabs(src[row][col]);
                tempval = tempval + (std::fabs(src[row][col + radius]) - oldVal) * rlen;
                dst[row][col] = tempval;
                ++pos;
                pos = pos <= radius ? pos : 0;
            }

            for (int col = W - radius; col < W; ++col) {
                tempval = (tempval * len - lineBuffer[pos]) / (len - 1);
                dst[row][col] = tempval;
                --len;
                ++pos;
                pos = pos <= radius ? pos : 0;
            }
        }

        //vertical blur
        float (* const rowBuffer)[numCols] = (float(*)[numCols]) buffer;
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int col = 0; col < W; col += numCols) {
            float len = radius + 1;

            for (int k = 0; k < numCols; ++k) {
                rowBuffer[0][k] = dst[0][col + k];
            }

            for (int i = 1; i <= radius; ++i) {
                for (int k = 0; k < numCols; ++k) {
                    dst[0][col + k] += dst[i][col + k];
                }
            }

            for(int k = 0; k < numCols; ++k) {
                dst[0][col + k] /= len;
            }

            for (int row = 1; row <= radius; ++row) {
                for(int k = 0; k < numCols; ++k) {
                    rowBuffer[row][k] = dst[row][col + k];
                    dst[row][col + k] = (dst[row - 1][col + k] * len + dst[row + radius][col + k]) / (len + 1);
                }

                ++len;
            }

            const float rlen = 1.f / len;
            int pos = 0;
            for (int row = radius + 1; row < H - radius; ++row) {
                for(int k = 0; k < numCols; ++k) {
                    float oldVal = rowBuffer[pos][k];
                    rowBuffer[pos][k] = dst[row][col + k];
                    dst[row][col + k] = dst[row - 1][col + k] + (dst[row + radius][col + k] - oldVal) * rlen;
                }
                ++pos;
                pos = pos <= radius ? pos : 0;
            }

            for (int row = H - radius; row < H; ++row) {
                for(int k = 0; k < numCols; ++k) {
                    dst[row][col + k] = (dst[row - 1][col + k] * len - rowBuffer[pos][k]) / (len - 1);
                }
                --len;
                ++pos;
                pos = pos <= radius ? pos : 0;
            }
        }
    }
}

void boxblur(float* src, float* dst, int radius, int W, int H, bool multiThread)
{
    float* srcp[H];
    float* dstp[H];
    for (int i = 0; i < H; ++i) {
        srcp[i] = src + i * W;
        dstp[i] = dst + i * W;
    }
    boxblur(srcp, dstp, radius, W, H, multiThread);
}

void boxabsblur(float* src, float* dst, int radius, int W, int H, bool multiThread)
{
    float* srcp[H];
    float* dstp[H];
    for (int i = 0; i < H; ++i) {
        srcp[i] = src + i * W;
        dstp[i] = dst + i * W;
    }
    boxabsblur(srcp, dstp, radius, W, H, multiThread);
}

}
