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
#include "alignedbuffer.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __SSE__
#if defined( WIN32 ) && defined(__x86_64__)
    #include <intrin.h>
#else
    #include <xmmintrin.h>
#endif
#endif

// classical filtering if the support window is small:

template<class T> void gaussHorizontal3 (T** src, T** dst, AlignedBufferMP<double> &buffer, int W, int H, const float c0, const float c1) {

#ifdef _OPENMP
#pragma omp for
#endif
    for (int i=0; i<H; i++) {
    	AlignedBuffer<double>* pBuf = buffer.acquire();
        T* temp=(T*)pBuf->data;

        for (int j=1; j<W-1; j++)
            temp[j] = (T)(c1 * (src[i][j-1] + src[i][j+1]) + c0 * src[i][j]);
        dst[i][0] = src[i][0];
        memcpy (dst[i]+1, temp+1, (W-2)*sizeof(T));

        buffer.release(pBuf);

        dst[i][W-1] = src[i][W-1];
    }
}

template<class T> void gaussVertical3 (T** src, T** dst, AlignedBufferMP<double> &buffer, int W, int H, const float c0, const float c1) {

#ifdef _OPENMP
#pragma omp for
#endif
    for (int i=0; i<W; i++) {
        AlignedBuffer<double>* pBuf = buffer.acquire();
    	T* temp = (T*)pBuf->data;

        for (int j = 1; j<H-1; j++)
        	temp[j] = (T)(c1 * (src[j-1][i] + src[j+1][i]) + c0 * src[j][i]);
        dst[0][i] = src[0][i];
	    for (int j=1; j<H-1; j++)
            dst[j][i] = temp[j];

        buffer.release(pBuf);

        dst[H-1][i] = src[H-1][i];
    }
}

#ifdef __SSE__
#ifdef WIN32
template<class T> __attribute__((force_align_arg_pointer)) void gaussVertical3Sse (T** src, T** dst, int W, int H, const float c0, const float c1) {
#else
template<class T> void gaussVertical3Sse (T** src, T** dst, int W, int H, const float c0, const float c1) {
#endif
    __m128 Tv,Tm1v,Tp1v;
    __m128 c0v,c1v;
    c0v = _mm_set1_ps(c0);
    c1v = _mm_set1_ps(c1);
#ifdef _OPENMP
#pragma omp for
#endif
    for (int i=0; i<W-3; i+=4) {
        Tm1v = _mm_loadu_ps( &src[0][i] );
        _mm_storeu_ps( &dst[0][i], Tm1v);
        if(H>1)
            Tv = _mm_loadu_ps( &src[1][i]);
        for (int j=1; j<H-1; j++){
            Tp1v = _mm_loadu_ps( &src[j+1][i]);
            _mm_storeu_ps( &dst[j][i], c1v * (Tp1v + Tm1v) + Tv * c0v);
            Tm1v = Tv;
            Tv = Tp1v;
        }
        _mm_storeu_ps( &dst[H-1][i], _mm_loadu_ps( &src[H-1][i]));
    }

// Borders are done without SSE
#ifdef _OPENMP
#pragma omp for
#endif
    for(int i=W-(W%4);i<W;i++)
        {
        dst[0][i] = src[0][i];
        for (int j = 1; j<H-1; j++)
        	dst[j][i] = c1 * (src[j-1][i] + src[j+1][i]) + c0 * src[j][i];
        dst[H-1][i] = src[H-1][i];
        }
}


#ifdef WIN32
template<class T> __attribute__((force_align_arg_pointer)) void gaussHorizontal3Sse (T** src, T** dst, int W, int H, const float c0, const float c1) {
#else
template<class T> void gaussHorizontal3Sse (T** src, T** dst, int W, int H, const float c0, const float c1) {
#endif
    float tmp[W][4] __attribute__ ((aligned (16)));

    __m128 Tv,Tm1v,Tp1v;
    __m128 c0v,c1v;
    c0v = _mm_set1_ps(c0);
    c1v = _mm_set1_ps(c1);
#ifdef _OPENMP
#pragma omp for
#endif
    for (int i=0; i<H-3; i+=4) {
        dst[i][0] = src[i][0];
        dst[i+1][0] = src[i+1][0];
        dst[i+2][0] = src[i+2][0];
        dst[i+3][0] = src[i+3][0];
        Tm1v = _mm_set_ps( src[i][0], src[i+1][0], src[i+2][0], src[i+3][0] );
        if(W>1)
            Tv = _mm_set_ps( src[i][1], src[i+1][1], src[i+2][1], src[i+3][1] );
        for (int j=1; j<W-1; j++){
            Tp1v = _mm_set_ps( src[i][j+1], src[i+1][j+1], src[i+2][j+1], src[i+3][j+1] );
            _mm_store_ps( &tmp[j][0], c1v * (Tp1v + Tm1v) + Tv * c0v);
            Tm1v = Tv;
            Tv = Tp1v;
        }

        for (int j=1; j<W-1; j++) {
            dst[i+3][j] = tmp[j][0];
            dst[i+2][j] = tmp[j][1];
            dst[i+1][j] = tmp[j][2];
            dst[i][j] = tmp[j][3];
        }

        dst[i][W-1] = src[i][W-1];
        dst[i+1][W-1] = src[i+1][W-1];
        dst[i+2][W-1] = src[i+2][W-1];
        dst[i+3][W-1] = src[i+3][W-1];
    }
// Borders are done without SSE
#ifdef _OPENMP
#pragma omp for
#endif
    for(int i=H-(H%4);i<H;i++)
        {
        dst[i][0] = src[i][0];
        for (int j = 1; j<W-1; j++)
        	dst[i][j] = c1 * (src[i][j-1] + src[i][j+1]) + c0 * src[i][j];
        dst[i][W-1] = src[i][W-1];
        }
}



// fast gaussian approximation if the support window is large
#ifdef WIN32
template<class T> __attribute__((force_align_arg_pointer)) void gaussHorizontalSse (T** src, T** dst, int W, int H, float sigma) {
#else
template<class T> void gaussHorizontalSse (T** src, T** dst, int W, int H, float sigma) {
#endif
    if (sigma<0.25) {
        // dont perform filtering
        if (src!=dst)
#pragma omp for
            for (int i = 0; i<H; i++)
                memcpy (dst[i], src[i], W*sizeof(T));
        return;
    }

    if (sigma<0.6) {
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
    if (sigma<2.5)
        q = 3.97156 - 4.14554 * sqrt (1.0 - 0.26891 * sigma);
    float b0 = 1.57825 + 2.44413*q + 1.4281*q*q + 0.422205*q*q*q;
    float b1 = 2.44413*q + 2.85619*q*q + 1.26661*q*q*q;
    float b2 = -1.4281*q*q - 1.26661*q*q*q;
    float b3 = 0.422205*q*q*q;
    float B = 1.0 - (b1+b2+b3) / b0;

    b1 /= b0;
    b2 /= b0;
    b3 /= b0;

    // From: Bill Triggs, Michael Sdika: Boundary Conditions for Young-van Vliet Recursive Filtering
    float M[3][3];
    M[0][0] = -b3*b1+1.0-b3*b3-b2;
    M[0][1] = (b3+b1)*(b2+b3*b1);
    M[0][2] = b3*(b1+b3*b2);
    M[1][0] = b1+b3*b2;
    M[1][1] = -(b2-1.0)*(b2+b3*b1);
    M[1][2] = -(b3*b1+b3*b3+b2-1.0)*b3;
    M[2][0] = b3*b1+b2+b1*b1-b2*b2;
    M[2][1] = b1*b2+b3*b2*b2-b1*b3*b3-b3*b3*b3-b3*b2+b3;
    M[2][2] = b3*(b1+b3*b2);
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++) {
            M[i][j] *= (1.0+b2+(b1-b3)*b3);
            M[i][j] /= (1.0+b1-b2+b3)*(1.0-b1-b2-b3);
        }
    float tmp[W][4] __attribute__ ((aligned (16)));
    float tmpV[4] __attribute__ ((aligned (16)));
    __m128 Rv;
    __m128 Tv,Tm2v,Tm3v;
    __m128 Bv,b1v,b2v,b3v;
    __m128 temp2W,temp2Wp1;
    Bv = _mm_set1_ps(B);
    b1v = _mm_set1_ps(b1);
    b2v = _mm_set1_ps(b2);
    b3v = _mm_set1_ps(b3);
#pragma omp for
    for (int i=0; i<H-3; i+=4) {
        tmpV[0] = src[i+3][0]; tmpV[1] = src[i+2][0]; tmpV[2] = src[i+1][0]; tmpV[3] = src[i][0];
        Tv = _mm_load_ps(tmpV);
        Rv = Tv * (Bv + b1v + b2v + b3v);
        Tm3v = Rv;
        _mm_store_ps( &tmp[0][0], Rv );

        tmpV[0] = src[i+3][1]; tmpV[1] = src[i+2][1]; tmpV[2] = src[i+1][1]; tmpV[3] = src[i][1];
        Rv = _mm_load_ps(tmpV) * Bv + Rv * b1v + Tv * (b2v + b3v);
        Tm2v = Rv;
        _mm_store_ps( &tmp[1][0], Rv );

        tmpV[0] = src[i+3][2]; tmpV[1] = src[i+2][2]; tmpV[2] = src[i+1][2]; tmpV[3] = src[i][2];
        Rv = _mm_load_ps(tmpV) * Bv + Rv * b1v + Tm3v * b2v + Tv * b3v;
        _mm_store_ps( &tmp[2][0], Rv );

        for (int j=3; j<W; j++) {
            Tv = Rv;
            Rv = _mm_set_ps(src[i][j],src[i+1][j],src[i+2][j],src[i+3][j]) * Bv + Tv * b1v + Tm2v * b2v + Tm3v * b3v;
            _mm_store_ps( &tmp[j][0], Rv );
            Tm3v = Tm2v;
            Tm2v = Tv;
        }

        Tv = _mm_set_ps(src[i][W-1],src[i+1][W-1],src[i+2][W-1],src[i+3][W-1]);

        temp2Wp1 = Tv + _mm_set1_ps(M[2][0]) * (Rv - Tv) + _mm_set1_ps(M[2][1]) * ( Tm2v - Tv ) +  _mm_set1_ps(M[2][2]) * (Tm3v - Tv);
        temp2W = Tv + _mm_set1_ps(M[1][0]) * (Rv - Tv) + _mm_set1_ps(M[1][1]) * (Tm2v - Tv) + _mm_set1_ps(M[1][2]) * (Tm3v - Tv);

        Rv = Tv + _mm_set1_ps(M[0][0]) * (Rv - Tv) + _mm_set1_ps(M[0][1]) * (Tm2v - Tv) + _mm_set1_ps(M[0][2]) * (Tm3v - Tv);
        _mm_store_ps( &tmp[W-1][0], Rv );

        Tm2v = Bv * Tm2v + b1v * Rv + b2v * temp2W + b3v * temp2Wp1;
        _mm_store_ps( &tmp[W-2][0], Tm2v );

        Tm3v = Bv * Tm3v + b1v * Tm2v + b2v * Rv + b3v * temp2W;
        _mm_store_ps( &tmp[W-3][0], Tm3v );

        Tv = Rv;
        Rv = Tm3v;
        Tm3v = Tv;

        for (int j=W-4; j>=0; j--) {
            Tv = Rv;
            Rv = _mm_load_ps(&tmp[j][0]) * Bv + Tv * b1v + Tm2v * b2v + Tm3v * b3v;
            _mm_store_ps( &tmp[j][0], Rv );
            Tm3v = Tm2v;
            Tm2v = Tv;
        }

        for (int j=0; j<W; j++) {
            dst[i+3][j] = tmp[j][0];
            dst[i+2][j] = tmp[j][1];
            dst[i+1][j] = tmp[j][2];
            dst[i][j] = tmp[j][3];
        }


    }
// Borders are done without SSE
#pragma omp for
        for(int i=H-(H%4);i<H;i++)
        {
        tmp[0][0] = B * src[i][0] + b1*src[i][0] + b2*src[i][0] + b3*src[i][0];
        tmp[1][0] = B * src[i][1] + b1*tmp[0][0]  + b2*src[i][0] + b3*src[i][0];
        tmp[2][0] = B * src[i][2] + b1*tmp[1][0]  + b2*tmp[0][0]  + b3*src[i][0];

        for (int j=3; j<W; j++)
            tmp[j][0] = B * src[i][j] + b1*tmp[j-1][0] + b2*tmp[j-2][0] + b3*tmp[j-3][0];

        float temp2Wm1 = src[i][W-1] + M[0][0]*(tmp[W-1][0] - src[i][W-1]) + M[0][1]*(tmp[W-2][0] - src[i][W-1]) + M[0][2]*(tmp[W-3][0] - src[i][W-1]);
        float temp2W   = src[i][W-1] + M[1][0]*(tmp[W-1][0] - src[i][W-1]) + M[1][1]*(tmp[W-2][0] - src[i][W-1]) + M[1][2]*(tmp[W-3][0] - src[i][W-1]);
        float temp2Wp1 = src[i][W-1] + M[2][0]*(tmp[W-1][0] - src[i][W-1]) + M[2][1]*(tmp[W-2][0] - src[i][W-1]) + M[2][2]*(tmp[W-3][0] - src[i][W-1]);

        tmp[W-1][0] = temp2Wm1;
        tmp[W-2][0] = B * tmp[W-2][0] + b1*tmp[W-1][0] + b2*temp2W + b3*temp2Wp1;
        tmp[W-3][0] = B * tmp[W-3][0] + b1*tmp[W-2][0] + b2*tmp[W-1][0] + b3*temp2W;

        for (int j=W-4; j>=0; j--)
            tmp[j][0] = B * tmp[j][0] + b1*tmp[j+1][0] + b2*tmp[j+2][0] + b3*tmp[j+3][0];

        for (int j=0; j<W; j++)
            dst[i][j] = tmp[j][0];

        }
}
#endif

// fast gaussian approximation if the support window is large

template<class T> void gaussHorizontal (T** src, T** dst, AlignedBufferMP<double> &buffer, int W, int H, double sigma) {

#ifdef __SSE__
	if(sigma < 70) { // bigger sigma only with double precision
		gaussHorizontalSse<T> (src, dst, W, H, sigma);
		return;
	}
#endif
    if (sigma<0.25) {
        // dont perform filtering
        if (src!=dst)
#pragma omp for
            for (int i = 0; i<H; i++)
                memcpy (dst[i], src[i], W*sizeof(T));
        return;
    }

    if (sigma<0.6) {
        // compute 3x3 kernel
        double c1 = exp (-1.0 / (2.0 * sigma * sigma));
        double csum = 2.0 * c1 + 1.0;
        c1 /= csum;
        double c0 = 1.0 / csum;
        gaussHorizontal3<T> (src, dst, buffer, W, H, c0, c1);
        return;
    }

    // coefficient calculation
    double q = 0.98711 * sigma - 0.96330;
    if (sigma<2.5)
        q = 3.97156 - 4.14554 * sqrt (1.0 - 0.26891 * sigma);
    double b0 = 1.57825 + 2.44413*q + 1.4281*q*q + 0.422205*q*q*q;
    double b1 = 2.44413*q + 2.85619*q*q + 1.26661*q*q*q;
    double b2 = -1.4281*q*q - 1.26661*q*q*q;
    double b3 = 0.422205*q*q*q;
    double B = 1.0 - (b1+b2+b3) / b0;

    b1 /= b0;
    b2 /= b0;
    b3 /= b0;

    // From: Bill Triggs, Michael Sdika: Boundary Conditions for Young-van Vliet Recursive Filtering
    double M[3][3];
    M[0][0] = -b3*b1+1.0-b3*b3-b2;
    M[0][1] = (b3+b1)*(b2+b3*b1);
    M[0][2] = b3*(b1+b3*b2);
    M[1][0] = b1+b3*b2;
    M[1][1] = -(b2-1.0)*(b2+b3*b1);
    M[1][2] = -(b3*b1+b3*b3+b2-1.0)*b3;
    M[2][0] = b3*b1+b2+b1*b1-b2*b2;
    M[2][1] = b1*b2+b3*b2*b2-b1*b3*b3-b3*b3*b3-b3*b2+b3;
    M[2][2] = b3*(b1+b3*b2);
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            M[i][j] /= (1.0+b1-b2+b3)*(1.0+b2+(b1-b3)*b3);

	#pragma omp for
    for (int i=0; i<H; i++) {
        AlignedBuffer<double>* pBuf = buffer.acquire();
        double* temp2 = pBuf->data;

        temp2[0] = B * src[i][0] + b1*src[i][0] + b2*src[i][0] + b3*src[i][0];
        temp2[1] = B * src[i][1] + b1*temp2[0]  + b2*src[i][0] + b3*src[i][0];
        temp2[2] = B * src[i][2] + b1*temp2[1]  + b2*temp2[0]  + b3*src[i][0];

        for (int j=3; j<W; j++)
            temp2[j] = B * src[i][j] + b1*temp2[j-1] + b2*temp2[j-2] + b3*temp2[j-3];

        double temp2Wm1 = src[i][W-1] + M[0][0]*(temp2[W-1] - src[i][W-1]) + M[0][1]*(temp2[W-2] - src[i][W-1]) + M[0][2]*(temp2[W-3] - src[i][W-1]);
        double temp2W   = src[i][W-1] + M[1][0]*(temp2[W-1] - src[i][W-1]) + M[1][1]*(temp2[W-2] - src[i][W-1]) + M[1][2]*(temp2[W-3] - src[i][W-1]);
        double temp2Wp1 = src[i][W-1] + M[2][0]*(temp2[W-1] - src[i][W-1]) + M[2][1]*(temp2[W-2] - src[i][W-1]) + M[2][2]*(temp2[W-3] - src[i][W-1]);

        temp2[W-1] = temp2Wm1;
        temp2[W-2] = B * temp2[W-2] + b1*temp2[W-1] + b2*temp2W + b3*temp2Wp1;
        temp2[W-3] = B * temp2[W-3] + b1*temp2[W-2] + b2*temp2[W-1] + b3*temp2W;

        for (int j=W-4; j>=0; j--)
            temp2[j] = B * temp2[j] + b1*temp2[j+1] + b2*temp2[j+2] + b3*temp2[j+3];
        for (int j=0; j<W; j++)
            dst[i][j] = (T)temp2[j];

        buffer.release(pBuf);
    }
}

#ifdef __SSE__
#ifdef WIN32
template<class T> __attribute__((force_align_arg_pointer)) void gaussVerticalSse (T** src, T** dst, int W, int H, float sigma) {
#else
template<class T> void gaussVerticalSse (T** src, T** dst, int W, int H, float sigma) {
#endif
    if (sigma<0.25) {
        // dont perform filtering
        if (src!=dst)
#pragma omp for
            for (int i = 0; i<H; i++)
                memcpy (dst[i], src[i], W*sizeof(T));
        return;
    }

    if (sigma<0.6) {
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
    if (sigma<2.5)
        q = 3.97156 - 4.14554 * sqrt (1.0 - 0.26891 * sigma);
    double b0 = 1.57825 + 2.44413*q + 1.4281*q*q + 0.422205*q*q*q;
    double b1 = 2.44413*q + 2.85619*q*q + 1.26661*q*q*q;
    double b2 = -1.4281*q*q - 1.26661*q*q*q;
    double b3 = 0.422205*q*q*q;
    double B = 1.0 - (b1+b2+b3) / b0;

    b1 /= b0;
    b2 /= b0;
    b3 /= b0;

    // From: Bill Triggs, Michael Sdika: Boundary Conditions for Young-van Vliet Recursive Filtering
    double M[3][3];
    M[0][0] = -b3*b1+1.0-b3*b3-b2;
    M[0][1] = (b3+b1)*(b2+b3*b1);
    M[0][2] = b3*(b1+b3*b2);
    M[1][0] = b1+b3*b2;
    M[1][1] = -(b2-1.0)*(b2+b3*b1);
    M[1][2] = -(b3*b1+b3*b3+b2-1.0)*b3;
    M[2][0] = b3*b1+b2+b1*b1-b2*b2;
    M[2][1] = b1*b2+b3*b2*b2-b1*b3*b3-b3*b3*b3-b3*b2+b3;
    M[2][2] = b3*(b1+b3*b2);
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++) {
            M[i][j] *= (1.0+b2+(b1-b3)*b3);
            M[i][j] /= (1.0+b1-b2+b3)*(1.0-b1-b2-b3);
        }
    float tmp[H][4] __attribute__ ((aligned (16)));
    __m128 Rv;
    __m128 Tv,Tm2v,Tm3v;
    __m128 Bv,b1v,b2v,b3v;
    __m128 temp2W,temp2Wp1;
    Bv = _mm_set1_ps(B);
    b1v = _mm_set1_ps(b1);
    b2v = _mm_set1_ps(b2);
    b3v = _mm_set1_ps(b3);


#ifdef _OPENMP
#pragma omp for
#endif
    for (int i=0; i<W-3; i+=4) {
        Tv = _mm_loadu_ps( &src[0][i]);
        Rv = Tv * (Bv + b1v + b2v + b3v);
        Tm3v = Rv;
        _mm_store_ps( &tmp[0][0], Rv );

        Rv = _mm_loadu_ps(&src[1][i]) * Bv + Rv * b1v + Tv * (b2v + b3v);
        Tm2v = Rv;
        _mm_store_ps( &tmp[1][0], Rv );

        Rv = _mm_loadu_ps(&src[2][i]) * Bv + Rv * b1v + Tm3v * b2v + Tv * b3v;
        _mm_store_ps( &tmp[2][0], Rv );

        for (int j=3; j<H; j++) {
            Tv = Rv;
            Rv = _mm_loadu_ps(&src[j][i]) * Bv +  Tv * b1v + Tm2v * b2v + Tm3v * b3v;
            _mm_store_ps( &tmp[j][0], Rv );
            Tm3v = Tm2v;
            Tm2v = Tv;
        }
        Tv = _mm_loadu_ps(&src[H-1][i]);

        temp2Wp1 = Tv + _mm_set1_ps(M[2][0]) * (Rv - Tv) + _mm_set1_ps(M[2][1]) * (Tm2v - Tv) + _mm_set1_ps(M[2][2]) * (Tm3v - Tv);
        temp2W = Tv + _mm_set1_ps(M[1][0]) * (Rv - Tv) + _mm_set1_ps(M[1][1]) * (Tm2v - Tv) + _mm_set1_ps(M[1][2]) * (Tm3v - Tv);

        Rv = Tv + _mm_set1_ps(M[0][0]) * (Rv - Tv) + _mm_set1_ps(M[0][1]) * (Tm2v - Tv) + _mm_set1_ps(M[0][2]) * (Tm3v - Tv);
        _mm_storeu_ps( &dst[H-1][i], Rv );

        Tm2v = Bv * Tm2v + b1v * Rv + b2v * temp2W + b3v * temp2Wp1;
        _mm_storeu_ps( &dst[H-2][i], Tm2v );

        Tm3v = Bv * Tm3v + b1v * Tm2v + b2v * Rv + b3v * temp2W;
        _mm_storeu_ps( &dst[H-3][i], Tm3v );

        Tv = Rv;
        Rv = Tm3v;
        Tm3v = Tv;

        for (int j=H-4; j>=0; j--) {
            Tv = Rv;
            Rv = _mm_load_ps(&tmp[j][0]) * Bv +  Tv * b1v + Tm2v * b2v + Tm3v * b3v;
            _mm_storeu_ps( &dst[j][i], Rv );
            Tm3v = Tm2v;
            Tm2v = Tv;
        }
    }
// Borders are done without SSE
#pragma omp for
    for(int i=W-(W%4);i<W;i++)
        {
    	tmp[0][0] = B * src[0][i] + b1*src[0][i] + b2*src[0][i] + b3*src[0][i];
        tmp[1][0] = B * src[1][i] + b1*tmp[0][0] + b2*src[0][i] + b3*src[0][i];
        tmp[2][0] = B * src[2][i] + b1*tmp[1][0] + b2*tmp[0][0] + b3*src[0][i];

        for (int j=3; j<H; j++)
            tmp[j][0] = B * src[j][i] + b1*tmp[j-1][0] + b2*tmp[j-2][0] + b3*tmp[j-3][0];

        float temp2Hm1 = src[H-1][i] + M[0][0]*(tmp[H-1][0] - src[H-1][i]) + M[0][1]*(tmp[H-2][0] - src[H-1][i]) + M[0][2]*(tmp[H-3][0] - src[H-1][i]);
        float temp2H   = src[H-1][i] + M[1][0]*(tmp[H-1][0] - src[H-1][i]) + M[1][1]*(tmp[H-2][0] - src[H-1][i]) + M[1][2]*(tmp[H-3][0] - src[H-1][i]);
        float temp2Hp1 = src[H-1][i] + M[2][0]*(tmp[H-1][0] - src[H-1][i]) + M[2][1]*(tmp[H-2][0] - src[H-1][i]) + M[2][2]*(tmp[H-3][0] - src[H-1][i]);

        tmp[H-1][0] = temp2Hm1;
        tmp[H-2][0] = B * tmp[H-2][0] + b1*tmp[H-1][0] + b2*temp2H + b3*temp2Hp1;
        tmp[H-3][0] = B * tmp[H-3][0] + b1*tmp[H-2][0] + b2*tmp[H-1][0] + b3*temp2H;

        for (int j=H-4; j>=0; j--)
            tmp[j][0] = B * tmp[j][0] + b1*tmp[j+1][0] + b2*tmp[j+2][0] + b3*tmp[j+3][0];

        for (int j=0; j<H; j++)
            dst[j][i] = tmp[j][0];

        }
}

#endif

template<class T> void gaussVertical (T** src, T** dst, AlignedBufferMP<double> &buffer, int W, int H, double sigma) {

#ifdef __SSE__
	if(sigma < 70) { // bigger sigma only with double precision
		gaussVerticalSse<T> (src, dst, W, H, sigma);
		return;
	}
#endif

    if (sigma<0.25) {
        // dont perform filtering
        if (src!=dst)
#pragma omp for
            for (int i = 0; i<H; i++)
                memcpy (dst[i], src[i], W*sizeof(T));
        return;
    }

    if (sigma<0.6) {
        // compute 3x3 kernel
        double c1 = exp (-1.0 / (2.0 * sigma * sigma));
        double csum = 2.0 * c1 + 1.0;
        c1 /= csum;
        double c0 = 1.0 / csum;
        gaussVertical3<T> (src, dst, buffer, W, H, c0, c1);
        return;
    }

    // coefficient calculation
    double q = 0.98711 * sigma - 0.96330;
    if (sigma<2.5)
        q = 3.97156 - 4.14554 * sqrt (1.0 - 0.26891 * sigma);
    double b0 = 1.57825 + 2.44413*q + 1.4281*q*q + 0.422205*q*q*q;
    double b1 = 2.44413*q + 2.85619*q*q + 1.26661*q*q*q;
    double b2 = -1.4281*q*q - 1.26661*q*q*q;
    double b3 = 0.422205*q*q*q;
    double B = 1.0 - (b1+b2+b3) / b0;

    b1 /= b0;
    b2 /= b0;
    b3 /= b0;

    // From: Bill Triggs, Michael Sdika: Boundary Conditions for Young-van Vliet Recursive Filtering
    double M[3][3];
    M[0][0] = -b3*b1+1.0-b3*b3-b2;
    M[0][1] = (b3+b1)*(b2+b3*b1);
    M[0][2] = b3*(b1+b3*b2);
    M[1][0] = b1+b3*b2;
    M[1][1] = -(b2-1.0)*(b2+b3*b1);
    M[1][2] = -(b3*b1+b3*b3+b2-1.0)*b3;
    M[2][0] = b3*b1+b2+b1*b1-b2*b2;
    M[2][1] = b1*b2+b3*b2*b2-b1*b3*b3-b3*b3*b3-b3*b2+b3;
    M[2][2] = b3*(b1+b3*b2);
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            M[i][j] /= (1.0+b1-b2+b3)*(1.0+b2+(b1-b3)*b3);
#ifdef _OPENMP
#pragma omp for
#endif
    for (int i=0; i<W; i++) {
        AlignedBuffer<double>* pBuf = buffer.acquire();
        double* temp2 = pBuf->data;
    	temp2[0] = B * src[0][i] + b1*src[0][i] + b2*src[0][i] + b3*src[0][i];
        temp2[1] = B * src[1][i] + b1*temp2[0]  + b2*src[0][i] + b3*src[0][i];
        temp2[2] = B * src[2][i] + b1*temp2[1]  + b2*temp2[0]  + b3*src[0][i];

        for (int j=3; j<H; j++)
            temp2[j] = B * src[j][i] + b1*temp2[j-1] + b2*temp2[j-2] + b3*temp2[j-3];

        double temp2Hm1 = src[H-1][i] + M[0][0]*(temp2[H-1] - src[H-1][i]) + M[0][1]*(temp2[H-2] - src[H-1][i]) + M[0][2]*(temp2[H-3] - src[H-1][i]);
        double temp2H   = src[H-1][i] + M[1][0]*(temp2[H-1] - src[H-1][i]) + M[1][1]*(temp2[H-2] - src[H-1][i]) + M[1][2]*(temp2[H-3] - src[H-1][i]);
        double temp2Hp1 = src[H-1][i] + M[2][0]*(temp2[H-1] - src[H-1][i]) + M[2][1]*(temp2[H-2] - src[H-1][i]) + M[2][2]*(temp2[H-3] - src[H-1][i]);

        temp2[H-1] = temp2Hm1;
        temp2[H-2] = B * temp2[H-2] + b1*temp2[H-1] + b2*temp2H + b3*temp2Hp1;
        temp2[H-3] = B * temp2[H-3] + b1*temp2[H-2] + b2*temp2[H-1] + b3*temp2H;

        for (int j=H-4; j>=0; j--)
            temp2[j] = B * temp2[j] + b1*temp2[j+1] + b2*temp2[j+2] + b3*temp2[j+3];

        for (int j=0; j<H; j++)
            dst[j][i] = (T)temp2[j];

        buffer.release(pBuf);
    }
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template<class T> void gaussDerivH (T** src, T** dst, AlignedBufferMP<double> &buffer, int W, int H, double sigma) {


    if (sigma<0.6) {
        // apply symmetric derivative
#ifdef _OPENMP
#pragma omp for
#endif
		for (int i=0; i<H; i++) {
			AlignedBuffer<double>* pBuf = buffer.acquire();
			T* temp = (T*)pBuf->data;
			// double* temp = buffer->data;// replaced by 2 lines above
			for (int j=1; j<W-1; j++)
				temp[j] = (0.5 * (src[i][j+1] - src[i][j-1]) );
			dst[i][0] = (src[i][1]-src[i][0]);
			//memcpy (dst[i]+1, temp+1, (W-2)*sizeof(T));
			for (int j=1; j<W-1; j++)
				dst[i][j] = temp[j];

			buffer.release(pBuf);
			dst[i][W-1] = (src[i][W-1]-src[i][W-2]);
		}
        return;
    }

    // coefficient calculation
    double q = 0.98711 * sigma - 0.96330;
    if (sigma<2.5)
        q = 3.97156 - 4.14554 * sqrt (1.0 - 0.26891 * sigma);
    double b0 = 1.57825 + 2.44413*q + 1.4281*q*q + 0.422205*q*q*q;
    double b1 = 2.44413*q + 2.85619*q*q + 1.26661*q*q*q;
    double b2 = -1.4281*q*q - 1.26661*q*q*q;
    double b3 = 0.422205*q*q*q;
    double B = 1.0 - (b1+b2+b3) / b0;

    b1 /= b0;
    b2 /= b0;
    b3 /= b0;

    // From: Bill Triggs, Michael Sdika: Boundary Conditions for Young-van Vliet Recursive Filtering
    double M[3][3];
    M[0][0] = -b3*b1+1.0-b3*b3-b2;
    M[0][1] = (b3+b1)*(b2+b3*b1);
    M[0][2] = b3*(b1+b3*b2);
    M[1][0] = b1+b3*b2;
    M[1][1] = -(b2-1.0)*(b2+b3*b1);
    M[1][2] = -(b3*b1+b3*b3+b2-1.0)*b3;
    M[2][0] = b3*b1+b2+b1*b1-b2*b2;
    M[2][1] = b1*b2+b3*b2*b2-b1*b3*b3-b3*b3*b3-b3*b2+b3;
    M[2][2] = b3*(b1+b3*b2);
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            M[i][j] /= (1.0+b1-b2+b3)*(1.0+b2+(b1-b3)*b3);

#pragma omp for
    for (int i=0; i<H; i++) {
		AlignedBuffer<double>* pBuf = buffer.acquire();
		T* temp2 = (T*)pBuf->data;
		// double* temp2 = buffer->data;// replaced by 2 lines above

		double src0 = (src[i][1]-src[i][0]);

        temp2[0] = B * src0 + b1*src0 + b2*src0 + b3*src0;
        temp2[1] = B * 0.5*(src[i][2]-src[i][0]) + b1*temp2[0]  + b2*src0 + b3*src0;
        temp2[2] = B * 0.5*(src[i][3]-src[i][1]) + b1*temp2[1]  + b2*temp2[0]  + b3*src0;

        for (int j=3; j<W-1; j++)
            temp2[j] = B * 0.5*(src[i][j+1]-src[i][j-1]) + b1*temp2[j-1] + b2*temp2[j-2] + b3*temp2[j-3];

		double srcWm1 = (src[i][W-1]-src[i][W-2]);

		temp2[W-1] = B * srcWm1 + b1*temp2[W-2] + b2*temp2[W-3] + b3*temp2[W-4];

        double temp2Wm1 = srcWm1 + M[0][0]*(temp2[W-1] - srcWm1) + M[0][1]*(temp2[W-2] - srcWm1) + M[0][2]*(temp2[W-3] - srcWm1);
        double temp2W   = srcWm1 + M[1][0]*(temp2[W-1] - srcWm1) + M[1][1]*(temp2[W-2] - srcWm1) + M[1][2]*(temp2[W-3] - srcWm1);
        double temp2Wp1 = srcWm1 + M[2][0]*(temp2[W-1] - srcWm1) + M[2][1]*(temp2[W-2] - srcWm1) + M[2][2]*(temp2[W-3] - srcWm1);

        temp2[W-1] = temp2Wm1;
        temp2[W-2] = B * temp2[W-2] + b1*temp2[W-1] + b2*temp2W + b3*temp2Wp1;
        temp2[W-3] = B * temp2[W-3] + b1*temp2[W-2] + b2*temp2[W-1] + b3*temp2W;

        for (int j=W-4; j>=0; j--)
            temp2[j] = B * temp2[j] + b1*temp2[j+1] + b2*temp2[j+2] + b3*temp2[j+3];
        for (int j=0; j<W; j++)
            dst[i][j] = (T)temp2[j];

        buffer.release(pBuf);
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template<class T> void gaussDerivV (T** src, T** dst, AlignedBufferMP<double> &buffer, int W, int H, double sigma) {

	if (sigma<0.6) {
        // apply symmetric derivative
#ifdef _OPENMP
#pragma omp for
#endif
		for (int j=0; j<W; j++) {
			AlignedBuffer<double>* pBuf = buffer.acquire();
			T* temp = (T*)pBuf->data;
			// double* temp = buffer->data;// replaced by 2 lines above
			for (int i = 1; i<H-1; i++)
				temp[i] = (0.5 * (src[i+1][j] - src[i-1][j]) );
			dst[0][j] = (src[1][j]-src[0][j]);
			for (int i=1; i<H-1; i++)
				dst[i][j] = temp[i];

			buffer.release(pBuf);

			dst[H-1][j] = (src[H-1][j]-src[H-2][j]);
		}
        return;
    }

    // coefficient calculation
    double q = 0.98711 * sigma - 0.96330;
    if (sigma<2.5)
        q = 3.97156 - 4.14554 * sqrt (1.0 - 0.26891 * sigma);
    double b0 = 1.57825 + 2.44413*q + 1.4281*q*q + 0.422205*q*q*q;
    double b1 = 2.44413*q + 2.85619*q*q + 1.26661*q*q*q;
    double b2 = -1.4281*q*q - 1.26661*q*q*q;
    double b3 = 0.422205*q*q*q;
    double B = 1.0 - (b1+b2+b3) / b0;

    b1 /= b0;
    b2 /= b0;
    b3 /= b0;

    // From: Bill Triggs, Michael Sdika: Boundary Conditions for Young-van Vliet Recursive Filtering
    double M[3][3];
    M[0][0] = -b3*b1+1.0-b3*b3-b2;
    M[0][1] = (b3+b1)*(b2+b3*b1);
    M[0][2] = b3*(b1+b3*b2);
    M[1][0] = b1+b3*b2;
    M[1][1] = -(b2-1.0)*(b2+b3*b1);
    M[1][2] = -(b3*b1+b3*b3+b2-1.0)*b3;
    M[2][0] = b3*b1+b2+b1*b1-b2*b2;
    M[2][1] = b1*b2+b3*b2*b2-b1*b3*b3-b3*b3*b3-b3*b2+b3;
    M[2][2] = b3*(b1+b3*b2);
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            M[i][j] /= (1.0+b1-b2+b3)*(1.0+b2+(b1-b3)*b3);
#ifdef _OPENMP
#pragma omp for
#endif
    for (int i=0; i<W; i++) {
		AlignedBuffer<double>* pBuf = buffer.acquire();
		T* temp2 = (T*)pBuf->data;
		// double* temp2 = buffer->data;// replaced by 2 lines above

		double src0 = 0.5*(src[1][i]-src[0][i]);

    	temp2[0] = B * src0 + b1*src0 + b2*src0 + b3*src0;
        temp2[1] = B * 0.5*(src[2][i]-src[0][i]) + b1*temp2[0]  + b2*src0 + b3*src0;
        temp2[2] = B * 0.5*(src[3][i]-src[1][i]) + b1*temp2[1]  + b2*temp2[0]  + b3*src0;

        for (int j=3; j<H-1; j++)
            temp2[j] = B * 0.5*(src[j+1][i]-src[j-1][i]) + b1*temp2[j-1] + b2*temp2[j-2] + b3*temp2[j-3];

		double srcHm1 = 0.5*(src[H-1][i]-src[H-2][i]);

		temp2[H-1] = B * srcHm1 + b1*temp2[H-2] + b2*temp2[H-3] + b3*temp2[H-4];

        double temp2Hm1 = srcHm1 + M[0][0]*(temp2[H-1] - srcHm1) + M[0][1]*(temp2[H-2] - srcHm1) + M[0][2]*(temp2[H-3] - srcHm1);
        double temp2H   = srcHm1 + M[1][0]*(temp2[H-1] - srcHm1) + M[1][1]*(temp2[H-2] - srcHm1) + M[1][2]*(temp2[H-3] - srcHm1);
        double temp2Hp1 = srcHm1 + M[2][0]*(temp2[H-1] - srcHm1) + M[2][1]*(temp2[H-2] - srcHm1) + M[2][2]*(temp2[H-3] - srcHm1);

        temp2[H-1] = temp2Hm1;
        temp2[H-2] = B * temp2[H-2] + b1*temp2[H-1] + b2*temp2H + b3*temp2Hp1;
        temp2[H-3] = B * temp2[H-3] + b1*temp2[H-2] + b2*temp2[H-1] + b3*temp2H;

        for (int j=H-4; j>=0; j--)
            temp2[j] = B * temp2[j] + b1*temp2[j+1] + b2*temp2[j+2] + b3*temp2[j+3];

        for (int j=0; j<H; j++)
            dst[j][i] = (T)temp2[j];

        buffer.release(pBuf);
    }
}

#endif
