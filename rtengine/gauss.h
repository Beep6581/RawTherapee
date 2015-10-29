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
 */
#ifndef _GAUSS_H_
#define _GAUSS_H_

#include <cstdlib>
#include <cstring>
#include <cmath>
#include "opthelper.h"
#include "stdio.h"
// classical filtering if the support window is small:

template<class T> void gaussHorizontal3 (T** src, T** dst, int W, int H, const float c0, const float c1)
{

#ifdef _OPENMP
    #pragma omp for
#endif

    for (int i = 0; i < H; i++) {
        T temp[W] ALIGNED16;

        for (int j = 1; j < W - 1; j++) {
            temp[j] = (T)(c1 * (src[i][j - 1] + src[i][j + 1]) + c0 * src[i][j]);
        }

        dst[i][0] = src[i][0];
        memcpy (dst[i] + 1, temp + 1, (W - 2)*sizeof(T));

        dst[i][W - 1] = src[i][W - 1];
    }
}

template<class T> void gaussVertical3 (T** src, T** dst, int W, int H, const float c0, const float c1)
{

#ifdef _OPENMP
    #pragma omp for
#endif

    for (int i = 0; i < W; i++) {
        T temp[H] ALIGNED16;

        for (int j = 1; j < H - 1; j++) {
            temp[j] = (T)(c1 * (src[j - 1][i] + src[j + 1][i]) + c0 * src[j][i]);
        }

        dst[0][i] = src[0][i];

        for (int j = 1; j < H - 1; j++) {
            dst[j][i] = temp[j];
        }

        dst[H - 1][i] = src[H - 1][i];
    }
}

#ifdef __SSE2__
template<class T> SSEFUNCTION void gaussVertical3Sse (T** src, T** dst, int W, int H, const float c0, const float c1)
{
    __m128 Tv, Tm1v, Tp1v;
    __m128 c0v, c1v;
    c0v = F2V(c0);
    c1v = F2V(c1);
#ifdef _OPENMP
    #pragma omp for
#endif

    for (int i = 0; i < W - 3; i += 4) {
        Tm1v = LVFU( src[0][i] );
        STVFU( dst[0][i], Tm1v);

        if (H > 1) {
            Tv = LVFU( src[1][i]);
        }

        for (int j = 1; j < H - 1; j++) {
            Tp1v = LVFU( src[j + 1][i]);
            STVFU( dst[j][i], c1v * (Tp1v + Tm1v) + Tv * c0v);
            Tm1v = Tv;
            Tv = Tp1v;
        }

        STVFU( dst[H - 1][i], LVFU( src[H - 1][i]));
    }

// Borders are done without SSE
#ifdef _OPENMP
    #pragma omp for
#endif

    for (int i = W - (W % 4); i < W; i++) {
        dst[0][i] = src[0][i];

        for (int j = 1; j < H - 1; j++) {
            dst[j][i] = c1 * (src[j - 1][i] + src[j + 1][i]) + c0 * src[j][i];
        }

        dst[H - 1][i] = src[H - 1][i];
    }
}


template<class T> SSEFUNCTION void gaussHorizontal3Sse (T** src, T** dst, int W, int H, const float c0, const float c1)
{
    float tmp[W][4] ALIGNED16;

    __m128 Tv, Tm1v, Tp1v;
    __m128 c0v, c1v;
    c0v = F2V(c0);
    c1v = F2V(c1);
#ifdef _OPENMP
    #pragma omp for
#endif

    for (int i = 0; i < H - 3; i += 4) {
        dst[i][0] = src[i][0];
        dst[i + 1][0] = src[i + 1][0];
        dst[i + 2][0] = src[i + 2][0];
        dst[i + 3][0] = src[i + 3][0];
        Tm1v = _mm_set_ps( src[i][0], src[i + 1][0], src[i + 2][0], src[i + 3][0] );

        if (W > 1) {
            Tv = _mm_set_ps( src[i][1], src[i + 1][1], src[i + 2][1], src[i + 3][1] );
        }

        for (int j = 1; j < W - 1; j++) {
            Tp1v = _mm_set_ps( src[i][j + 1], src[i + 1][j + 1], src[i + 2][j + 1], src[i + 3][j + 1] );
            STVF( tmp[j][0], c1v * (Tp1v + Tm1v) + Tv * c0v);
            Tm1v = Tv;
            Tv = Tp1v;
        }

        for (int j = 1; j < W - 1; j++) {
            dst[i + 3][j] = tmp[j][0];
            dst[i + 2][j] = tmp[j][1];
            dst[i + 1][j] = tmp[j][2];
            dst[i][j] = tmp[j][3];
        }

        dst[i][W - 1] = src[i][W - 1];
        dst[i + 1][W - 1] = src[i + 1][W - 1];
        dst[i + 2][W - 1] = src[i + 2][W - 1];
        dst[i + 3][W - 1] = src[i + 3][W - 1];
    }

// Borders are done without SSE
#ifdef _OPENMP
    #pragma omp for
#endif

    for (int i = H - (H % 4); i < H; i++) {
        dst[i][0] = src[i][0];

        for (int j = 1; j < W - 1; j++) {
            dst[i][j] = c1 * (src[i][j - 1] + src[i][j + 1]) + c0 * src[i][j];
        }

        dst[i][W - 1] = src[i][W - 1];
    }
}



// fast gaussian approximation if the support window is large
template<class T> SSEFUNCTION void gaussHorizontalSse (T** src, T** dst, int W, int H, float sigma)
{

    if (sigma < 0.25) {
        // dont perform filtering
        if (src != dst)
#ifdef _OPENMP
            #pragma omp for
#endif
            for (int i = 0; i < H; i++) {
                memcpy (dst[i], src[i], W * sizeof(T));
            }

        return;
    }

    if (sigma < 0.6) {
        // compute 3x3 kernel
        float c1 = exp (-1.0 / (2.0 * sigma * sigma));
        float csum = 2.0 * c1 + 1.0;
        c1 /= csum;
        float c0 = 1.0 / csum;
        gaussHorizontal3Sse<T> (src, dst, W, H, c0, c1);
        return;
    }

    // coefficient calculation
    float q = 0.98711 * sigma - 0.96330;

    if (sigma < 2.5) {
        q = 3.97156 - 4.14554 * sqrt (1.0 - 0.26891 * sigma);
    }

    float b0 = 1.57825 + 2.44413 * q + 1.4281 * q * q + 0.422205 * q * q * q;
    float b1 = 2.44413 * q + 2.85619 * q * q + 1.26661 * q * q * q;
    float b2 = -1.4281 * q * q - 1.26661 * q * q * q;
    float b3 = 0.422205 * q * q * q;
    float B = 1.0 - (b1 + b2 + b3) / b0;

    b1 /= b0;
    b2 /= b0;
    b3 /= b0;

    // From: Bill Triggs, Michael Sdika: Boundary Conditions for Young-van Vliet Recursive Filtering
    float M[3][3];
    M[0][0] = -b3 * b1 + 1.0 - b3 * b3 - b2;
    M[0][1] = (b3 + b1) * (b2 + b3 * b1);
    M[0][2] = b3 * (b1 + b3 * b2);
    M[1][0] = b1 + b3 * b2;
    M[1][1] = -(b2 - 1.0) * (b2 + b3 * b1);
    M[1][2] = -(b3 * b1 + b3 * b3 + b2 - 1.0) * b3;
    M[2][0] = b3 * b1 + b2 + b1 * b1 - b2 * b2;
    M[2][1] = b1 * b2 + b3 * b2 * b2 - b1 * b3 * b3 - b3 * b3 * b3 - b3 * b2 + b3;
    M[2][2] = b3 * (b1 + b3 * b2);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            M[i][j] *= (1.0 + b2 + (b1 - b3) * b3);
            M[i][j] /= (1.0 + b1 - b2 + b3) * (1.0 - b1 - b2 - b3);
        }

    float tmp[W][4] ALIGNED16;
    float tmpV[4] ALIGNED16;
    __m128 Rv;
    __m128 Tv, Tm2v, Tm3v;
    __m128 Bv, b1v, b2v, b3v;
    __m128 temp2W, temp2Wp1;
    Bv = F2V(B);
    b1v = F2V(b1);
    b2v = F2V(b2);
    b3v = F2V(b3);

#ifdef _OPENMP
    #pragma omp for
#endif

    for (int i = 0; i < H - 3; i += 4) {
        tmpV[0] = src[i + 3][0];
        tmpV[1] = src[i + 2][0];
        tmpV[2] = src[i + 1][0];
        tmpV[3] = src[i][0];
        Tv = LVF(tmpV[0]);
        Rv = Tv * (Bv + b1v + b2v + b3v);
        Tm3v = Rv;
        STVF( tmp[0][0], Rv );

        tmpV[0] = src[i + 3][1];
        tmpV[1] = src[i + 2][1];
        tmpV[2] = src[i + 1][1];
        tmpV[3] = src[i][1];
        Rv = LVF(tmpV[0]) * Bv + Rv * b1v + Tv * (b2v + b3v);
        Tm2v = Rv;
        STVF( tmp[1][0], Rv );

        tmpV[0] = src[i + 3][2];
        tmpV[1] = src[i + 2][2];
        tmpV[2] = src[i + 1][2];
        tmpV[3] = src[i][2];
        Rv = LVF(tmpV[0]) * Bv + Rv * b1v + Tm3v * b2v + Tv * b3v;
        STVF( tmp[2][0], Rv );

        for (int j = 3; j < W; j++) {
            Tv = Rv;
            Rv = _mm_set_ps(src[i][j], src[i + 1][j], src[i + 2][j], src[i + 3][j]) * Bv + Tv * b1v + Tm2v * b2v + Tm3v * b3v;
            STVF( tmp[j][0], Rv );
            Tm3v = Tm2v;
            Tm2v = Tv;
        }

        Tv = _mm_set_ps(src[i][W - 1], src[i + 1][W - 1], src[i + 2][W - 1], src[i + 3][W - 1]);

        temp2Wp1 = Tv + F2V(M[2][0]) * (Rv - Tv) + F2V(M[2][1]) * ( Tm2v - Tv ) +  F2V(M[2][2]) * (Tm3v - Tv);
        temp2W = Tv + F2V(M[1][0]) * (Rv - Tv) + F2V(M[1][1]) * (Tm2v - Tv) + F2V(M[1][2]) * (Tm3v - Tv);

        Rv = Tv + F2V(M[0][0]) * (Rv - Tv) + F2V(M[0][1]) * (Tm2v - Tv) + F2V(M[0][2]) * (Tm3v - Tv);
        STVF( tmp[W - 1][0], Rv );

        Tm2v = Bv * Tm2v + b1v * Rv + b2v * temp2W + b3v * temp2Wp1;
        STVF( tmp[W - 2][0], Tm2v );

        Tm3v = Bv * Tm3v + b1v * Tm2v + b2v * Rv + b3v * temp2W;
        STVF( tmp[W - 3][0], Tm3v );

        Tv = Rv;
        Rv = Tm3v;
        Tm3v = Tv;

        for (int j = W - 4; j >= 0; j--) {
            Tv = Rv;
            Rv = LVF(tmp[j][0]) * Bv + Tv * b1v + Tm2v * b2v + Tm3v * b3v;
            STVF( tmp[j][0], Rv );
            Tm3v = Tm2v;
            Tm2v = Tv;
        }

        for (int j = 0; j < W; j++) {
            dst[i + 3][j] = tmp[j][0];
            dst[i + 2][j] = tmp[j][1];
            dst[i + 1][j] = tmp[j][2];
            dst[i][j] = tmp[j][3];
        }


    }

// Borders are done without SSE
#ifdef _OPENMP
    #pragma omp for
#endif

    for (int i = H - (H % 4); i < H; i++) {
        tmp[0][0] = B * src[i][0] + b1 * src[i][0] + b2 * src[i][0] + b3 * src[i][0];
        tmp[1][0] = B * src[i][1] + b1 * tmp[0][0]  + b2 * src[i][0] + b3 * src[i][0];
        tmp[2][0] = B * src[i][2] + b1 * tmp[1][0]  + b2 * tmp[0][0]  + b3 * src[i][0];

        for (int j = 3; j < W; j++) {
            tmp[j][0] = B * src[i][j] + b1 * tmp[j - 1][0] + b2 * tmp[j - 2][0] + b3 * tmp[j - 3][0];
        }

        float temp2Wm1 = src[i][W - 1] + M[0][0] * (tmp[W - 1][0] - src[i][W - 1]) + M[0][1] * (tmp[W - 2][0] - src[i][W - 1]) + M[0][2] * (tmp[W - 3][0] - src[i][W - 1]);
        float temp2W   = src[i][W - 1] + M[1][0] * (tmp[W - 1][0] - src[i][W - 1]) + M[1][1] * (tmp[W - 2][0] - src[i][W - 1]) + M[1][2] * (tmp[W - 3][0] - src[i][W - 1]);
        float temp2Wp1 = src[i][W - 1] + M[2][0] * (tmp[W - 1][0] - src[i][W - 1]) + M[2][1] * (tmp[W - 2][0] - src[i][W - 1]) + M[2][2] * (tmp[W - 3][0] - src[i][W - 1]);

        tmp[W - 1][0] = temp2Wm1;
        tmp[W - 2][0] = B * tmp[W - 2][0] + b1 * tmp[W - 1][0] + b2 * temp2W + b3 * temp2Wp1;
        tmp[W - 3][0] = B * tmp[W - 3][0] + b1 * tmp[W - 2][0] + b2 * tmp[W - 1][0] + b3 * temp2W;

        for (int j = W - 4; j >= 0; j--) {
            tmp[j][0] = B * tmp[j][0] + b1 * tmp[j + 1][0] + b2 * tmp[j + 2][0] + b3 * tmp[j + 3][0];
        }

        for (int j = 0; j < W; j++) {
            dst[i][j] = tmp[j][0];
        }

    }
}
#endif

// fast gaussian approximation if the support window is large

template<class T> void gaussHorizontal (T** src, T** dst, int W, int H, double sigma)
{

#ifdef __SSE2__

    if (sigma < 70) { // bigger sigma only with double precision
        gaussHorizontalSse<T> (src, dst, W, H, sigma);
        return;
    }

#endif

    if (sigma < 0.25) {
        // dont perform filtering
        if (src != dst)
#ifdef _OPENMP
            #pragma omp for
#endif
            for (int i = 0; i < H; i++) {
                memcpy (dst[i], src[i], W * sizeof(T));
            }

        return;
    }

    if (sigma < 0.6) {
        // compute 3x3 kernel
        double c1 = exp (-1.0 / (2.0 * sigma * sigma));
        double csum = 2.0 * c1 + 1.0;
        c1 /= csum;
        double c0 = 1.0 / csum;
        gaussHorizontal3<T> (src, dst, W, H, c0, c1);
        return;
    }

    // coefficient calculation
    double q = 0.98711 * sigma - 0.96330;

    if (sigma < 2.5) {
        q = 3.97156 - 4.14554 * sqrt (1.0 - 0.26891 * sigma);
    }

    double b0 = 1.57825 + 2.44413 * q + 1.4281 * q * q + 0.422205 * q * q * q;
    double b1 = 2.44413 * q + 2.85619 * q * q + 1.26661 * q * q * q;
    double b2 = -1.4281 * q * q - 1.26661 * q * q * q;
    double b3 = 0.422205 * q * q * q;
    double B = 1.0 - (b1 + b2 + b3) / b0;

    b1 /= b0;
    b2 /= b0;
    b3 /= b0;

    // From: Bill Triggs, Michael Sdika: Boundary Conditions for Young-van Vliet Recursive Filtering
    double M[3][3];
    M[0][0] = -b3 * b1 + 1.0 - b3 * b3 - b2;
    M[0][1] = (b3 + b1) * (b2 + b3 * b1);
    M[0][2] = b3 * (b1 + b3 * b2);
    M[1][0] = b1 + b3 * b2;
    M[1][1] = -(b2 - 1.0) * (b2 + b3 * b1);
    M[1][2] = -(b3 * b1 + b3 * b3 + b2 - 1.0) * b3;
    M[2][0] = b3 * b1 + b2 + b1 * b1 - b2 * b2;
    M[2][1] = b1 * b2 + b3 * b2 * b2 - b1 * b3 * b3 - b3 * b3 * b3 - b3 * b2 + b3;
    M[2][2] = b3 * (b1 + b3 * b2);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            M[i][j] /= (1.0 + b1 - b2 + b3) * (1.0 + b2 + (b1 - b3) * b3);
        }

    double temp2[W] ALIGNED16;

#ifdef _OPENMP
    #pragma omp for
#endif

    for (int i = 0; i < H; i++) {

        temp2[0] = B * src[i][0] + b1 * src[i][0] + b2 * src[i][0] + b3 * src[i][0];
        temp2[1] = B * src[i][1] + b1 * temp2[0]  + b2 * src[i][0] + b3 * src[i][0];
        temp2[2] = B * src[i][2] + b1 * temp2[1]  + b2 * temp2[0]  + b3 * src[i][0];

        for (int j = 3; j < W; j++) {
            temp2[j] = B * src[i][j] + b1 * temp2[j - 1] + b2 * temp2[j - 2] + b3 * temp2[j - 3];
        }

        double temp2Wm1 = src[i][W - 1] + M[0][0] * (temp2[W - 1] - src[i][W - 1]) + M[0][1] * (temp2[W - 2] - src[i][W - 1]) + M[0][2] * (temp2[W - 3] - src[i][W - 1]);
        double temp2W   = src[i][W - 1] + M[1][0] * (temp2[W - 1] - src[i][W - 1]) + M[1][1] * (temp2[W - 2] - src[i][W - 1]) + M[1][2] * (temp2[W - 3] - src[i][W - 1]);
        double temp2Wp1 = src[i][W - 1] + M[2][0] * (temp2[W - 1] - src[i][W - 1]) + M[2][1] * (temp2[W - 2] - src[i][W - 1]) + M[2][2] * (temp2[W - 3] - src[i][W - 1]);

        temp2[W - 1] = temp2Wm1;
        temp2[W - 2] = B * temp2[W - 2] + b1 * temp2[W - 1] + b2 * temp2W + b3 * temp2Wp1;
        temp2[W - 3] = B * temp2[W - 3] + b1 * temp2[W - 2] + b2 * temp2[W - 1] + b3 * temp2W;

        for (int j = W - 4; j >= 0; j--) {
            temp2[j] = B * temp2[j] + b1 * temp2[j + 1] + b2 * temp2[j + 2] + b3 * temp2[j + 3];
        }

        for (int j = 0; j < W; j++) {
            dst[i][j] = (T)temp2[j];
        }

    }
}

#ifdef __SSE2__
template<class T> SSEFUNCTION void gaussVerticalSse (T** src, T** dst, int W, int H, float sigma)
{

    if (sigma < 0.25) {
        // dont perform filtering
        if (src != dst)
#ifdef _OPENMP
            #pragma omp for
#endif
            for (int i = 0; i < H; i++) {
                memcpy (dst[i], src[i], W * sizeof(T));
            }

        return;
    }

    if (sigma < 0.6) {
        // compute 3x3 kernel
        double c1 = exp (-1.0 / (2.0 * sigma * sigma));
        double csum = 2.0 * c1 + 1.0;
        c1 /= csum;
        double c0 = 1.0 / csum;
        gaussVertical3Sse<T> (src, dst, W, H, c0, c1);
        return;
    }

    // coefficient calculation
    double q = 0.98711 * sigma - 0.96330;

    if (sigma < 2.5) {
        q = 3.97156 - 4.14554 * sqrt (1.0 - 0.26891 * sigma);
    }

    double b0 = 1.57825 + 2.44413 * q + 1.4281 * q * q + 0.422205 * q * q * q;
    double b1 = 2.44413 * q + 2.85619 * q * q + 1.26661 * q * q * q;
    double b2 = -1.4281 * q * q - 1.26661 * q * q * q;
    double b3 = 0.422205 * q * q * q;
    double B = 1.0 - (b1 + b2 + b3) / b0;

    b1 /= b0;
    b2 /= b0;
    b3 /= b0;

    // From: Bill Triggs, Michael Sdika: Boundary Conditions for Young-van Vliet Recursive Filtering
    double M[3][3];
    M[0][0] = -b3 * b1 + 1.0 - b3 * b3 - b2;
    M[0][1] = (b3 + b1) * (b2 + b3 * b1);
    M[0][2] = b3 * (b1 + b3 * b2);
    M[1][0] = b1 + b3 * b2;
    M[1][1] = -(b2 - 1.0) * (b2 + b3 * b1);
    M[1][2] = -(b3 * b1 + b3 * b3 + b2 - 1.0) * b3;
    M[2][0] = b3 * b1 + b2 + b1 * b1 - b2 * b2;
    M[2][1] = b1 * b2 + b3 * b2 * b2 - b1 * b3 * b3 - b3 * b3 * b3 - b3 * b2 + b3;
    M[2][2] = b3 * (b1 + b3 * b2);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            M[i][j] *= (1.0 + b2 + (b1 - b3) * b3);
            M[i][j] /= (1.0 + b1 - b2 + b3) * (1.0 - b1 - b2 - b3);
        }

    float tmp[H][4] ALIGNED16;
    __m128 Rv;
    __m128 Tv, Tm2v, Tm3v;
    __m128 Bv, b1v, b2v, b3v;
    __m128 temp2W, temp2Wp1;
    Bv = F2V(B);
    b1v = F2V(b1);
    b2v = F2V(b2);
    b3v = F2V(b3);


#ifdef _OPENMP
    #pragma omp for
#endif

    for (int i = 0; i < W - 3; i += 4) {
        Tv = LVFU( src[0][i]);
        Rv = Tv * (Bv + b1v + b2v + b3v);
        Tm3v = Rv;
        STVF( tmp[0][0], Rv );

        Rv = LVFU(src[1][i]) * Bv + Rv * b1v + Tv * (b2v + b3v);
        Tm2v = Rv;
        STVF( tmp[1][0], Rv );

        Rv = LVFU(src[2][i]) * Bv + Rv * b1v + Tm3v * b2v + Tv * b3v;
        STVF( tmp[2][0], Rv );

        for (int j = 3; j < H; j++) {
            Tv = Rv;
            Rv = LVFU(src[j][i]) * Bv +  Tv * b1v + Tm2v * b2v + Tm3v * b3v;
            STVF( tmp[j][0], Rv );
            Tm3v = Tm2v;
            Tm2v = Tv;
        }

        Tv = LVFU(src[H - 1][i]);

        temp2Wp1 = Tv + F2V(M[2][0]) * (Rv - Tv) + F2V(M[2][1]) * (Tm2v - Tv) + F2V(M[2][2]) * (Tm3v - Tv);
        temp2W = Tv + F2V(M[1][0]) * (Rv - Tv) + F2V(M[1][1]) * (Tm2v - Tv) + F2V(M[1][2]) * (Tm3v - Tv);

        Rv = Tv + F2V(M[0][0]) * (Rv - Tv) + F2V(M[0][1]) * (Tm2v - Tv) + F2V(M[0][2]) * (Tm3v - Tv);
        STVFU( dst[H - 1][i], Rv );

        Tm2v = Bv * Tm2v + b1v * Rv + b2v * temp2W + b3v * temp2Wp1;
        STVFU( dst[H - 2][i], Tm2v );

        Tm3v = Bv * Tm3v + b1v * Tm2v + b2v * Rv + b3v * temp2W;
        STVFU( dst[H - 3][i], Tm3v );

        Tv = Rv;
        Rv = Tm3v;
        Tm3v = Tv;

        for (int j = H - 4; j >= 0; j--) {
            Tv = Rv;
            Rv = LVF(tmp[j][0]) * Bv +  Tv * b1v + Tm2v * b2v + Tm3v * b3v;
            STVFU( dst[j][i], Rv );
            Tm3v = Tm2v;
            Tm2v = Tv;
        }
    }

// Borders are done without SSE
#ifdef _OPENMP
    #pragma omp for
#endif

    for (int i = W - (W % 4); i < W; i++) {
        tmp[0][0] = B * src[0][i] + b1 * src[0][i] + b2 * src[0][i] + b3 * src[0][i];
        tmp[1][0] = B * src[1][i] + b1 * tmp[0][0] + b2 * src[0][i] + b3 * src[0][i];
        tmp[2][0] = B * src[2][i] + b1 * tmp[1][0] + b2 * tmp[0][0] + b3 * src[0][i];

        for (int j = 3; j < H; j++) {
            tmp[j][0] = B * src[j][i] + b1 * tmp[j - 1][0] + b2 * tmp[j - 2][0] + b3 * tmp[j - 3][0];
        }

        float temp2Hm1 = src[H - 1][i] + M[0][0] * (tmp[H - 1][0] - src[H - 1][i]) + M[0][1] * (tmp[H - 2][0] - src[H - 1][i]) + M[0][2] * (tmp[H - 3][0] - src[H - 1][i]);
        float temp2H   = src[H - 1][i] + M[1][0] * (tmp[H - 1][0] - src[H - 1][i]) + M[1][1] * (tmp[H - 2][0] - src[H - 1][i]) + M[1][2] * (tmp[H - 3][0] - src[H - 1][i]);
        float temp2Hp1 = src[H - 1][i] + M[2][0] * (tmp[H - 1][0] - src[H - 1][i]) + M[2][1] * (tmp[H - 2][0] - src[H - 1][i]) + M[2][2] * (tmp[H - 3][0] - src[H - 1][i]);

        tmp[H - 1][0] = temp2Hm1;
        tmp[H - 2][0] = B * tmp[H - 2][0] + b1 * tmp[H - 1][0] + b2 * temp2H + b3 * temp2Hp1;
        tmp[H - 3][0] = B * tmp[H - 3][0] + b1 * tmp[H - 2][0] + b2 * tmp[H - 1][0] + b3 * temp2H;

        for (int j = H - 4; j >= 0; j--) {
            tmp[j][0] = B * tmp[j][0] + b1 * tmp[j + 1][0] + b2 * tmp[j + 2][0] + b3 * tmp[j + 3][0];
        }

        for (int j = 0; j < H; j++) {
            dst[j][i] = tmp[j][0];
        }

    }
}

#endif

template<class T> void gaussVertical (T** src, T** dst, int W, int H, double sigma)
{

#ifdef __SSE2__

    if (sigma < 70) { // bigger sigma only with double precision
        gaussVerticalSse<T> (src, dst, W, H, sigma);
        return;
    }

#endif

    if (sigma < 0.25) {
        // don't perform filtering
        if (src != dst)
#ifdef _OPENMP
            #pragma omp for
#endif
            for (int i = 0; i < H; i++) {
                memcpy (dst[i], src[i], W * sizeof(T));
            }

        return;
    }

    if (sigma < 0.6) {
        // compute 3x3 kernel
        double c1 = exp (-1.0 / (2.0 * sigma * sigma));
        double csum = 2.0 * c1 + 1.0;
        c1 /= csum;
        double c0 = 1.0 / csum;
        gaussVertical3<T> (src, dst, W, H, c0, c1);
        return;
    }

    // coefficient calculation
    double q = 0.98711 * sigma - 0.96330;

    if (sigma < 2.5) {
        q = 3.97156 - 4.14554 * sqrt (1.0 - 0.26891 * sigma);
    }

    double b0 = 1.57825 + 2.44413 * q + 1.4281 * q * q + 0.422205 * q * q * q;
    double b1 = 2.44413 * q + 2.85619 * q * q + 1.26661 * q * q * q;
    double b2 = -1.4281 * q * q - 1.26661 * q * q * q;
    double b3 = 0.422205 * q * q * q;
    double B = 1.0 - (b1 + b2 + b3) / b0;

    b1 /= b0;
    b2 /= b0;
    b3 /= b0;

    // From: Bill Triggs, Michael Sdika: Boundary Conditions for Young-van Vliet Recursive Filtering
    double M[3][3];
    M[0][0] = -b3 * b1 + 1.0 - b3 * b3 - b2;
    M[0][1] = (b3 + b1) * (b2 + b3 * b1);
    M[0][2] = b3 * (b1 + b3 * b2);
    M[1][0] = b1 + b3 * b2;
    M[1][1] = -(b2 - 1.0) * (b2 + b3 * b1);
    M[1][2] = -(b3 * b1 + b3 * b3 + b2 - 1.0) * b3;
    M[2][0] = b3 * b1 + b2 + b1 * b1 - b2 * b2;
    M[2][1] = b1 * b2 + b3 * b2 * b2 - b1 * b3 * b3 - b3 * b3 * b3 - b3 * b2 + b3;
    M[2][2] = b3 * (b1 + b3 * b2);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            M[i][j] /= (1.0 + b1 - b2 + b3) * (1.0 + b2 + (b1 - b3) * b3);
        }

    // process 'numcols' columns for better usage of L1 cpu cache (especially faster for large values of H)
    static const int numcols = 8;
    double temp2[H][numcols] ALIGNED16;
    double temp2Hm1[numcols], temp2H[numcols], temp2Hp1[numcols];
#ifdef _OPENMP
    #pragma omp for nowait
#endif

    for (int i = 0; i < W - numcols + 1; i += numcols) {
        for (int k = 0; k < numcols; k++) {
            temp2[0][k] = B * src[0][i + k] + b1 * src[0][i + k] + b2 * src[0][i + k] + b3 * src[0][i + k];
            temp2[1][k] = B * src[1][i + k] + b1 * temp2[0][k] + b2 * src[0][i + k] + b3 * src[0][i + k];
            temp2[2][k] = B * src[2][i + k] + b1 * temp2[1][k] + b2 * temp2[0][k] + b3 * src[0][i + k];
        }

        for (int j = 3; j < H; j++) {
            for (int k = 0; k < numcols; k++) {
                temp2[j][k] = B * src[j][i + k] + b1 * temp2[j - 1][k] + b2 * temp2[j - 2][k] + b3 * temp2[j - 3][k];
            }
        }

        for (int k = 0; k < numcols; k++) {
            temp2Hm1[k] = src[H - 1][i + k] + M[0][0] * (temp2[H - 1][k] - src[H - 1][i + k]) + M[0][1] * (temp2[H - 2][k] - src[H - 1][i + k]) + M[0][2] * (temp2[H - 3][k] - src[H - 1][i + k]);
            temp2H[k]   = src[H - 1][i + k] + M[1][0] * (temp2[H - 1][k] - src[H - 1][i + k]) + M[1][1] * (temp2[H - 2][k] - src[H - 1][i + k]) + M[1][2] * (temp2[H - 3][k] - src[H - 1][i + k]);
            temp2Hp1[k] = src[H - 1][i + k] + M[2][0] * (temp2[H - 1][k] - src[H - 1][i + k]) + M[2][1] * (temp2[H - 2][k] - src[H - 1][i + k]) + M[2][2] * (temp2[H - 3][k] - src[H - 1][i + k]);
        }

        for (int k = 0; k < numcols; k++) {
            dst[H - 1][i + k] = temp2[H - 1][k] = temp2Hm1[k];
            dst[H - 2][i + k] = temp2[H - 2][k] = B * temp2[H - 2][k] + b1 * temp2[H - 1][k] + b2 * temp2H[k] + b3 * temp2Hp1[k];
            dst[H - 3][i + k] = temp2[H - 3][k] = B * temp2[H - 3][k] + b1 * temp2[H - 2][k] + b2 * temp2[H - 1][k] + b3 * temp2H[k];
        }

        for (int j = H - 4; j >= 0; j--) {
            for (int k = 0; k < numcols; k++) {
                dst[j][i + k] = temp2[j][k] = B * temp2[j][k] + b1 * temp2[j + 1][k] + b2 * temp2[j + 2][k] + b3 * temp2[j + 3][k];
            }
        }
    }

#ifdef _OPENMP
    #pragma omp single
#endif

    // process remaining column
    for (int i = W - (W % numcols); i < W; i++) {
        temp2[0][0] = B * src[0][i] + b1 * src[0][i] + b2 * src[0][i] + b3 * src[0][i];
        temp2[1][0] = B * src[1][i] + b1 * temp2[0][0]  + b2 * src[0][i] + b3 * src[0][i];
        temp2[2][0] = B * src[2][i] + b1 * temp2[1][0]  + b2 * temp2[0][0]  + b3 * src[0][i];

        for (int j = 3; j < H; j++) {
            temp2[j][0] = B * src[j][i] + b1 * temp2[j - 1][0] + b2 * temp2[j - 2][0] + b3 * temp2[j - 3][0];
        }

        double temp2Hm1 = src[H - 1][i] + M[0][0] * (temp2[H - 1][0] - src[H - 1][i]) + M[0][1] * (temp2[H - 2][0] - src[H - 1][i]) + M[0][2] * (temp2[H - 3][0] - src[H - 1][i]);
        double temp2H   = src[H - 1][i] + M[1][0] * (temp2[H - 1][0] - src[H - 1][i]) + M[1][1] * (temp2[H - 2][0] - src[H - 1][i]) + M[1][2] * (temp2[H - 3][0] - src[H - 1][i]);
        double temp2Hp1 = src[H - 1][i] + M[2][0] * (temp2[H - 1][0] - src[H - 1][i]) + M[2][1] * (temp2[H - 2][0] - src[H - 1][i]) + M[2][2] * (temp2[H - 3][0] - src[H - 1][i]);

        dst[H - 1][i] = temp2[H - 1][0] = temp2Hm1;
        dst[H - 2][i] = temp2[H - 2][0] = B * temp2[H - 2][0] + b1 * temp2[H - 1][0] + b2 * temp2H + b3 * temp2Hp1;
        dst[H - 3][i] = temp2[H - 3][0] = B * temp2[H - 3][0] + b1 * temp2[H - 2][0] + b2 * temp2[H - 1][0] + b3 * temp2H;

        for (int j = H - 4; j >= 0; j--) {
            dst[j][i] = temp2[j][0] = B * temp2[j][0] + b1 * temp2[j + 1][0] + b2 * temp2[j + 2][0] + b3 * temp2[j + 3][0];
        }
    }
}

template<class T> void gaussianBlur(T** src, T** dst, const int W, const int H, const double sigma, bool forceLowSigma = false)
{
    double newSigma = sigma;

    if(forceLowSigma && newSigma > 170.f) {
        newSigma /= sqrt(2.0);

        if(newSigma < 0.6) { // barrier to avoid using simple gauss version for higher radius
            newSigma = sigma;
            forceLowSigma = false;
            gaussianBlur(src,dst,W,H,newSigma,forceLowSigma);
        } else {
            gaussianBlur(src,dst,W,H,newSigma,forceLowSigma);
            gaussianBlur(dst,dst,W,H,newSigma,forceLowSigma);
        }
    } else {
        gaussHorizontal<T> (src, dst, W, H, newSigma);
        gaussVertical<T>   (dst, dst, W, H, newSigma);

    }
//    #pragma omp critical
//    printf("gauss sigma : %f / %f\n",sigma,newSigma);

/*
    if(forceLowSigma && newSigma > 170.f) {
        gaussianBlur(dst,dst,W,H,newSigma,forceLowSigma);
//        gaussHorizontal<T> (dst, dst, W, H, newSigma);
//        gaussVertical<T>   (dst, dst, W, H, newSigma);
    }
*/
}

#endif
