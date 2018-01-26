/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (C) 2010 Emil Martinec <ejmartin@uchicago.edu>
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
 */
#ifndef _BOXBLUR_H_
#define _BOXBLUR_H_

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "alignedbuffer.h"
#include "rt_math.h"
#include "opthelper.h"


namespace rtengine
{

// classical filtering if the support window is small:

template<class T, class A> void boxblur (T** src, A** dst, int radx, int rady, int W, int H)
{
    //box blur image; box range = (radx,rady)

    AlignedBuffer<float>* buffer = new AlignedBuffer<float> (W * H);
    float* temp = buffer->data;

    if (radx == 0) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int row = 0; row < H; row++)
            for (int col = 0; col < W; col++) {
                temp[row * W + col] = (float)src[row][col];
            }
    } else {
        //horizontal blur
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int row = 0; row < H; row++) {
            int len = radx + 1;
            temp[row * W + 0] = (float)src[row][0] / len;

            for (int j = 1; j <= radx; j++) {
                temp[row * W + 0] += (float)src[row][j] / len;
            }

            for (int col = 1; col <= radx; col++) {
                temp[row * W + col] = (temp[row * W + col - 1] * len + (float)src[row][col + radx]) / (len + 1);
                len ++;
            }

            for (int col = radx + 1; col < W - radx; col++) {
                temp[row * W + col] = temp[row * W + col - 1] + ((float)(src[row][col + radx] - src[row][col - radx - 1])) / len;
            }

            for (int col = W - radx; col < W; col++) {
                temp[row * W + col] = (temp[row * W + col - 1] * len - src[row][col - radx - 1]) / (len - 1);
                len --;
            }
        }
    }

    if (rady == 0) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int row = 0; row < H; row++)
            for (int col = 0; col < W; col++) {
                dst[row][col] = temp[row * W + col];
            }
    } else {
        //vertical blur
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int col = 0; col < W; col++) {
            int len = rady + 1;
            dst[0][col] = temp[0 * W + col] / len;

            for (int i = 1; i <= rady; i++) {
                dst[0][col] += temp[i * W + col] / len;
            }

            for (int row = 1; row <= rady; row++) {
                dst[row][col] = (dst[(row - 1)][col] * len + temp[(row + rady) * W + col]) / (len + 1);
                len ++;
            }

            for (int row = rady + 1; row < H - rady; row++) {
                dst[row][col] = dst[(row - 1)][col] + (temp[(row + rady) * W + col] - temp[(row - rady - 1) * W + col]) / len;
            }

            for (int row = H - rady; row < H; row++) {
                dst[row][col] = (dst[(row - 1)][col] * len - temp[(row - rady - 1) * W + col]) / (len - 1);
                len --;
            }
        }
    }

    delete buffer;

}

template<class T, class A> void boxblur (T** src, A** dst, T* buffer, int radx, int rady, int W, int H)
{
    //box blur image; box range = (radx,rady)

    float* temp = buffer;

    if (radx == 0) {
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int row = 0; row < H; row++)
            for (int col = 0; col < W; col++) {
                temp[row * W + col] = (float)src[row][col];
            }
    } else {
        //horizontal blur
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int row = 0; row < H; row++) {
            float len = radx + 1;
            float tempval = (float)src[row][0];

            for (int j = 1; j <= radx; j++) {
                tempval += (float)src[row][j];
            }

            tempval /= len;
            temp[row * W + 0] = tempval;

            for (int col = 1; col <= radx; col++) {
                temp[row * W + col] = tempval = (tempval * len + (float)src[row][col + radx]) / (len + 1);
                len ++;
            }

            for (int col = radx + 1; col < W - radx; col++) {
                temp[row * W + col] = tempval = tempval + ((float)(src[row][col + radx] - src[row][col - radx - 1])) / len;
            }

            for (int col = W - radx; col < W; col++) {
                temp[row * W + col] = tempval = (tempval * len - src[row][col - radx - 1]) / (len - 1);
                len --;
            }
        }
    }

    if (rady == 0) {
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int row = 0; row < H; row++)
            for (int col = 0; col < W; col++) {
                dst[row][col] = temp[row * W + col];
            }
    } else {
        const int numCols = 8; // process numCols columns at once for better usage of L1 cpu cache
#ifdef __SSE2__
        vfloat  leninitv = F2V( (float)(rady + 1));
        vfloat  onev = F2V( 1.f );
        vfloat  tempv, temp1v, lenv, lenp1v, lenm1v, rlenv;

#ifdef _OPENMP
        #pragma omp for
#endif

        for (int col = 0; col < W - 7; col += 8) {
            lenv = leninitv;
            tempv = LVFU(temp[0 * W + col]);
            temp1v = LVFU(temp[0 * W + col + 4]);

            for (int i = 1; i <= rady; i++) {
                tempv = tempv + LVFU(temp[i * W + col]);
                temp1v = temp1v + LVFU(temp[i * W + col + 4]);
            }

            tempv = tempv / lenv;
            temp1v = temp1v / lenv;
            STVFU( dst[0][col], tempv);
            STVFU( dst[0][col + 4], temp1v);

            for (int row = 1; row <= rady; row++) {
                lenp1v = lenv + onev;
                tempv = (tempv * lenv + LVFU(temp[(row + rady) * W + col])) / lenp1v;
                temp1v = (temp1v * lenv + LVFU(temp[(row + rady) * W + col + 4])) / lenp1v;
                STVFU( dst[row][col], tempv);
                STVFU( dst[row][col + 4], temp1v);
                lenv = lenp1v;
            }

            rlenv = onev / lenv;

            for (int row = rady + 1; row < H - rady; row++) {
                tempv = tempv + (LVFU(temp[(row + rady) * W + col]) - LVFU(temp[(row - rady - 1) * W + col])) * rlenv ;
                temp1v = temp1v + (LVFU(temp[(row + rady) * W + col + 4]) - LVFU(temp[(row - rady - 1) * W + col + 4])) * rlenv ;
                STVFU( dst[row][col], tempv);
                STVFU( dst[row][col + 4], temp1v);
            }

            for (int row = H - rady; row < H; row++) {
                lenm1v = lenv - onev;
                tempv = (tempv * lenv - LVFU(temp[(row - rady - 1) * W + col])) / lenm1v;
                temp1v = (temp1v * lenv - LVFU(temp[(row - rady - 1) * W + col + 4])) / lenm1v;
                STVFU( dst[row][col], tempv);
                STVFU( dst[row][col + 4], temp1v);
                lenv = lenm1v;
            }
        }

#else
        //vertical blur
#ifdef _OPENMP
        #pragma omp for
#endif

        for (int col = 0; col < W - numCols + 1; col += 8) {
            float len = rady + 1;

            for(int k = 0; k < numCols; k++) {
                dst[0][col + k] = temp[0 * W + col + k];
            }

            for (int i = 1; i <= rady; i++) {
                for(int k = 0; k < numCols; k++) {
                    dst[0][col + k] += temp[i * W + col + k];
                }
            }

            for(int k = 0; k < numCols; k++) {
                dst[0][col + k] /= len;
            }

            for (int row = 1; row <= rady; row++) {
                for(int k = 0; k < numCols; k++) {
                    dst[row][col + k] = (dst[(row - 1)][col + k] * len + temp[(row + rady) * W + col + k]) / (len + 1);
                }

                len ++;
            }

            for (int row = rady + 1; row < H - rady; row++) {
                for(int k = 0; k < numCols; k++) {
                    dst[row][col + k] = dst[(row - 1)][col + k] + (temp[(row + rady) * W + col + k] - temp[(row - rady - 1) * W + col + k]) / len;
                }
            }

            for (int row = H - rady; row < H; row++) {
                for(int k = 0; k < numCols; k++) {
                    dst[row][col + k] = (dst[(row - 1)][col + k] * len - temp[(row - rady - 1) * W + col + k]) / (len - 1);
                }

                len --;
            }
        }

#endif
#ifdef _OPENMP
        #pragma omp single
#endif

        for (int col = W - (W % numCols); col < W; col++) {
            float len = rady + 1;
            dst[0][col] = temp[0 * W + col] / len;

            for (int i = 1; i <= rady; i++) {
                dst[0][col] += temp[i * W + col] / len;
            }

            for (int row = 1; row <= rady; row++) {
                dst[row][col] = (dst[(row - 1)][col] * len + temp[(row + rady) * W + col]) / (len + 1);
                len ++;
            }

            for (int row = rady + 1; row < H - rady; row++) {
                dst[row][col] = dst[(row - 1)][col] + (temp[(row + rady) * W + col] - temp[(row - rady - 1) * W + col]) / len;
            }

            for (int row = H - rady; row < H; row++) {
                dst[row][col] = (dst[(row - 1)][col] * len - temp[(row - rady - 1) * W + col]) / (len - 1);
                len --;
            }
        }
    }

}

template<class T, class A> void boxblur (T* src, A* dst, A* buffer, int radx, int rady, int W, int H)
{
    //box blur image; box range = (radx,rady) i.e. box size is (2*radx+1)x(2*rady+1)

    float* temp = buffer;

    if (radx == 0) {
        for (int row = 0; row < H; row++)
            for (int col = 0; col < W; col++) {
                temp[row * W + col] = src[row * W + col];
            }
    } else {
        //horizontal blur
        for (int row = H - 1; row >= 0; row--) {
            int len = radx + 1;
            float tempval = (float)src[row * W];

            for (int j = 1; j <= radx; j++) {
                tempval += (float)src[row * W + j];
            }

            tempval = tempval / len;
            temp[row * W] = tempval;

            for (int col = 1; col <= radx; col++) {
                tempval = (tempval * len + src[row * W + col + radx]) / (len + 1);
                temp[row * W + col] = tempval;
                len ++;
            }

            float reclen = 1.f / len;

            for (int col = radx + 1; col < W - radx; col++) {
                tempval = tempval + ((float)(src[row * W + col + radx] - src[row * W + col - radx - 1])) * reclen;
                temp[row * W + col] = tempval;
            }

            for (int col = W - radx; col < W; col++) {
                tempval = (tempval * len - src[row * W + col - radx - 1]) / (len - 1);
                temp[row * W + col] = tempval;
                len --;
            }
        }
    }

    if (rady == 0) {
        for (int row = 0; row < H; row++)
            for (int col = 0; col < W; col++) {
                dst[row * W + col] = temp[row * W + col];
            }
    } else {
        //vertical blur
#ifdef __SSE2__
        vfloat  leninitv = F2V( (float)(rady + 1));
        vfloat  onev = F2V( 1.f );
        vfloat  tempv, temp1v, lenv, lenp1v, lenm1v, rlenv;
        int col;

        for (col = 0; col < W - 7; col += 8) {
            lenv = leninitv;
            tempv = LVFU(temp[0 * W + col]);
            temp1v = LVFU(temp[0 * W + col + 4]);

            for (int i = 1; i <= rady; i++) {
                tempv = tempv + LVFU(temp[i * W + col]);
                temp1v = temp1v + LVFU(temp[i * W + col + 4]);
            }

            tempv = tempv / lenv;
            temp1v = temp1v / lenv;
            STVFU( dst[0 * W + col], tempv);
            STVFU( dst[0 * W + col + 4], temp1v);

            for (int row = 1; row <= rady; row++) {
                lenp1v = lenv + onev;
                tempv = (tempv * lenv + LVFU(temp[(row + rady) * W + col])) / lenp1v;
                temp1v = (temp1v * lenv + LVFU(temp[(row + rady) * W + col + 4])) / lenp1v;
                STVFU( dst[row * W + col], tempv);
                STVFU( dst[row * W + col + 4], temp1v);
                lenv = lenp1v;
            }

            rlenv = onev / lenv;

            for (int row = rady + 1; row < H - rady; row++) {
                tempv = tempv + (LVFU(temp[(row + rady) * W + col]) - LVFU(temp[(row - rady - 1) * W + col])) * rlenv ;
                temp1v = temp1v + (LVFU(temp[(row + rady) * W + col + 4]) - LVFU(temp[(row - rady - 1) * W + col + 4])) * rlenv ;
                STVFU( dst[row * W + col], tempv);
                STVFU( dst[row * W + col + 4], temp1v);
            }

            for (int row = H - rady; row < H; row++) {
                lenm1v = lenv - onev;
                tempv = (tempv * lenv - LVFU(temp[(row - rady - 1) * W + col])) / lenm1v;
                temp1v = (temp1v * lenv - LVFU(temp[(row - rady - 1) * W + col + 4])) / lenm1v;
                STVFU( dst[row * W + col], tempv);
                STVFU( dst[row * W + col + 4], temp1v);
                lenv = lenm1v;
            }
        }

        for (; col < W - 3; col += 4) {
            lenv = leninitv;
            tempv = LVFU(temp[0 * W + col]);

            for (int i = 1; i <= rady; i++) {
                tempv = tempv + LVFU(temp[i * W + col]);
            }

            tempv = tempv / lenv;
            STVFU( dst[0 * W + col], tempv);

            for (int row = 1; row <= rady; row++) {
                lenp1v = lenv + onev;
                tempv = (tempv * lenv + LVFU(temp[(row + rady) * W + col])) / lenp1v;
                STVFU( dst[row * W + col], tempv);
                lenv = lenp1v;
            }

            rlenv = onev / lenv;

            for (int row = rady + 1; row < H - rady; row++) {
                tempv = tempv + (LVFU(temp[(row + rady) * W + col]) - LVFU(temp[(row - rady - 1) * W + col])) * rlenv ;
                STVFU( dst[row * W + col], tempv);
            }

            for (int row = H - rady; row < H; row++) {
                lenm1v = lenv - onev;
                tempv = (tempv * lenv - LVFU(temp[(row - rady - 1) * W + col])) / lenm1v;
                STVFU( dst[row * W + col], tempv);
                lenv = lenm1v;
            }
        }

        for (; col < W; col++) {
            int len = rady + 1;
            dst[0 * W + col] = temp[0 * W + col] / len;

            for (int i = 1; i <= rady; i++) {
                dst[0 * W + col] += temp[i * W + col] / len;
            }

            for (int row = 1; row <= rady; row++) {
                dst[row * W + col] = (dst[(row - 1) * W + col] * len + temp[(row + rady) * W + col]) / (len + 1);
                len ++;
            }

            for (int row = rady + 1; row < H - rady; row++) {
                dst[row * W + col] = dst[(row - 1) * W + col] + (temp[(row + rady) * W + col] - temp[(row - rady - 1) * W + col]) / len;
            }

            for (int row = H - rady; row < H; row++) {
                dst[row * W + col] = (dst[(row - 1) * W + col] * len - temp[(row - rady - 1) * W + col]) / (len - 1);
                len --;
            }
        }

#else

        for (int col = 0; col < W; col++) {
            int len = rady + 1;
            dst[0 * W + col] = temp[0 * W + col] / len;

            for (int i = 1; i <= rady; i++) {
                dst[0 * W + col] += temp[i * W + col] / len;
            }

            for (int row = 1; row <= rady; row++) {
                dst[row * W + col] = (dst[(row - 1) * W + col] * len + temp[(row + rady) * W + col]) / (len + 1);
                len ++;
            }

            for (int row = rady + 1; row < H - rady; row++) {
                dst[row * W + col] = dst[(row - 1) * W + col] + (temp[(row + rady) * W + col] - temp[(row - rady - 1) * W + col]) / len;
            }

            for (int row = H - rady; row < H; row++) {
                dst[row * W + col] = (dst[(row - 1) * W + col] * len - temp[(row - rady - 1) * W + col]) / (len - 1);
                len --;
            }
        }

#endif
    }

}

template<class T, class A> void boxabsblur (T* src, A* dst, int radx, int rady, int W, int H, float * temp)
{

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //box blur image; box range = (radx,rady) i.e. box size is (2*radx+1)x(2*rady+1)

    if (radx == 0) {
        for (int row = 0; row < H; row++)
            for (int col = 0; col < W; col++) {
                temp[row * W + col] = fabs(src[row * W + col]);
            }
    } else {
        //horizontal blur
        for (int row = 0; row < H; row++) {
            int len = radx + 1;
            float tempval = fabsf((float)src[row * W + 0]);

            for (int j = 1; j <= radx; j++) {
                tempval += fabsf((float)src[row * W + j]);
            }

            tempval /= len;
            temp[row * W + 0] = tempval;

            for (int col = 1; col <= radx; col++) {
                tempval = (tempval * len + fabsf(src[row * W + col + radx])) / (len + 1);
                temp[row * W + col] = tempval;
                len ++;
            }

            float rlen = 1.f / (float)len;

            for (int col = radx + 1; col < W - radx; col++) {
                tempval = tempval + ((float)(fabsf(src[row * W + col + radx]) - fabsf(src[row * W + col - radx - 1]))) * rlen;
                temp[row * W + col] = tempval;
            }

            for (int col = W - radx; col < W; col++) {
                tempval = (tempval * len - fabsf(src[row * W + col - radx - 1])) / (len - 1);
                temp[row * W + col] = tempval;
                len --;
            }
        }
    }

    if (rady == 0) {
        for (int row = 0; row < H; row++)
            for (int col = 0; col < W; col++) {
                dst[row * W + col] = temp[row * W + col];
            }
    } else {
        //vertical blur
#ifdef __SSE2__
        vfloat  leninitv = F2V( (float)(rady + 1));
        vfloat  onev = F2V( 1.f );
        vfloat  tempv, lenv, lenp1v, lenm1v, rlenv;

        for (int col = 0; col < W - 3; col += 4) {
            lenv = leninitv;
            tempv = LVF(temp[0 * W + col]);

            for (int i = 1; i <= rady; i++) {
                tempv = tempv + LVF(temp[i * W + col]);
            }

            tempv = tempv / lenv;
            STVF(dst[0 * W + col], tempv);

            for (int row = 1; row <= rady; row++) {
                lenp1v = lenv + onev;
                tempv = (tempv * lenv + LVF(temp[(row + rady) * W + col])) / lenp1v;
                STVF(dst[row * W + col], tempv);
                lenv = lenp1v;
            }

            rlenv = onev / lenv;

            for (int row = rady + 1; row < H - rady; row++) {
                tempv = tempv + (LVF(temp[(row + rady) * W + col]) - LVF(temp[(row - rady - 1) * W + col])) * rlenv;
                STVF(dst[row * W + col], tempv);
            }

            for (int row = H - rady; row < H; row++) {
                lenm1v = lenv - onev;
                tempv = (tempv * lenv - LVF(temp[(row - rady - 1) * W + col])) / lenm1v;
                STVF(dst[row * W + col], tempv);
                lenv = lenm1v;
            }
        }

        for (int col = W - (W % 4); col < W; col++) {
            int len = rady + 1;
            dst[0 * W + col] = temp[0 * W + col] / len;

            for (int i = 1; i <= rady; i++) {
                dst[0 * W + col] += temp[i * W + col] / len;
            }

            for (int row = 1; row <= rady; row++) {
                dst[row * W + col] = (dst[(row - 1) * W + col] * len + temp[(row + rady) * W + col]) / (len + 1);
                len ++;
            }

            for (int row = rady + 1; row < H - rady; row++) {
                dst[row * W + col] = dst[(row - 1) * W + col] + (temp[(row + rady) * W + col] - temp[(row - rady - 1) * W + col]) / len;
            }

            for (int row = H - rady; row < H; row++) {
                dst[row * W + col] = (dst[(row - 1) * W + col] * len - temp[(row - rady - 1) * W + col]) / (len - 1);
                len --;
            }
        }

#else

        for (int col = 0; col < W; col++) {
            int len = rady + 1;
            dst[0 * W + col] = temp[0 * W + col] / len;

            for (int i = 1; i <= rady; i++) {
                dst[0 * W + col] += temp[i * W + col] / len;
            }

            for (int row = 1; row <= rady; row++) {
                dst[row * W + col] = (dst[(row - 1) * W + col] * len + temp[(row + rady) * W + col]) / (len + 1);
                len ++;
            }

            for (int row = rady + 1; row < H - rady; row++) {
                dst[row * W + col] = dst[(row - 1) * W + col] + (temp[(row + rady) * W + col] - temp[(row - rady - 1) * W + col]) / len;
            }

            for (int row = H - rady; row < H; row++) {
                dst[row * W + col] = (dst[(row - 1) * W + col] * len - temp[(row - rady - 1) * W + col]) / (len - 1);
                len --;
            }
        }

#endif
    }

}

}
#endif /* _BOXBLUR_H_ */
