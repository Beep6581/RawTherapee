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

#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <alignedbuffer.h>

#define NOSSE 1

#ifndef NOSSE
#include <xmmintrin.h>
template<class T> void gaussHorizontal3 (T** src, T** dst, T* buffer, int W, int row_from, int row_to, const float c0, const float c1) {
    
    typedef float pfloat[4];
    pfloat* temp = (pfloat*)buffer + 1;
    pfloat* tmp = (pfloat*)buffer;

    __m128 xmm1; xmm1 = _mm_load1_ps (&c0);
    __m128 xmm2; xmm2 = _mm_load1_ps (&c1);
    __m128 xmm3; __m128 xmm4; __m128 xmm5; __m128 xmm6;
    __m128 xmm0;

    int i;
    for (i = row_from; i<row_to-3; i+=4) {
        T* row1 = src[i];
        T* row2 = src[i+1];
        T* row3 = src[i+2];
        T* row4 = src[i+3];
        (*tmp)[0] = (float)row1[0];
        (*tmp)[1] = (float)row2[0];
        (*tmp)[2] = (float)row3[0];
        (*tmp)[3] = (float)row4[0];
        xmm3 = _mm_load_ps ((float*)tmp);   // previous element
        (*tmp)[0] = (float)row1[1];
        (*tmp)[1] = (float)row2[1];
        (*tmp)[2] = (float)row3[1];
        (*tmp)[3] = (float)row4[1];
        xmm4 = _mm_load_ps ((float*)tmp);   // current element

        for (int j=1; j<W-1; j++) {
            (*tmp)[0] = (float)row1[j+1];
            (*tmp)[1] = (float)row2[j+1];
            (*tmp)[2] = (float)row3[j+1];
            (*tmp)[3] = (float)row4[j+1];
            xmm5 = _mm_load_ps ((float*)tmp);   // next element
            // compute blured pixel
            xmm0 = _mm_mul_ps (xmm3, xmm2);
            xmm6 = _mm_mul_ps (xmm5, xmm2);
            xmm0 = _mm_add_ps (xmm6, xmm0);
            xmm6 = _mm_mul_ps (xmm4, xmm1);
            xmm0 = _mm_add_ps (xmm6, xmm0);
            // store blured pixel
            _mm_store_ps ((float*)&temp[j], xmm0);
            // update prev and current
            xmm3 = xmm4;
            xmm4 = xmm5;
        }
        for (int j=1; j<W-1; j++) {
          dst[i][j]   = (T)temp[j][0];
          dst[i+1][j] = (T)temp[j][1];
          dst[i+2][j] = (T)temp[j][2];
          dst[i+3][j] = (T)temp[j][3];
        }                    
        dst[i][0] = src[i][0];
        dst[i+1][0] = src[i+1][0];
        dst[i+2][0] = src[i+2][0];
        dst[i+3][0] = src[i+3][0];
        dst[i][W-1] = src[i][W-1];
        dst[i+1][W-1] = src[i+1][W-1];
        dst[i+2][W-1] = src[i+2][W-1];
        dst[i+3][W-1] = src[i+3][W-1];
    }
    for (; i<row_to; i++) {
        for (int j=1; j<W-1; j++)
            buffer[j] = (T)(c1 * (src[i][j-1] + src[i][j+1]) + c0 * src[i][j]);
        dst[i][0] = src[i][0];
        memcpy (dst[i]+1, buffer+1, (W-2)*sizeof(T));
        dst[i][W-1] = src[i][W-1];
    }
}

/*template<class T> void gaussHorizontal3 (T** src, T** dst, T* buffer, int W, int row_from, int row_to, const int c0, const int c1) {
    
    time_t t1 = clock ();

    const int csum = c0 + 2 * c1;    
    for (int i = row_from; i<row_to; i++) {
        for (int j=1; j<W-1; j++)
            buffer[j] = (c1 * (src[i][j-1] + src[i][j+1]) + c0 * src[i][j]) / csum;
        dst[i][0] = src[i][0];
        memcpy (dst[i]+1, buffer+1, (W-2)*sizeof(T));
        dst[i][W-1] = src[i][W-1];
    }
    printf ("horizontal time = %d\n", clock()-t1);
}
*/

template<class T> void gaussVertical3 (T** src, T** dst, T* buffer, int H, int col_from, int col_to, const float c0, const float c1) {
    
    typedef float pfloat[4];
    pfloat* temp = (pfloat*)buffer + 1;
    pfloat* tmp = (pfloat*)buffer;

    __m128 xmm1; xmm1 = _mm_load1_ps (&c0);
    __m128 xmm2; xmm2 = _mm_load1_ps (&c1);
    __m128 xmm3; __m128 xmm4; __m128 xmm5; __m128 xmm6;
    __m128 xmm0;

    int i;
    for (i = col_from; i<col_to-3; i+=4) {
        (*tmp)[0] = (float)src[0][i];
        (*tmp)[1] = (float)src[0][i+1];
        (*tmp)[2] = (float)src[0][i+2];
        (*tmp)[3] = (float)src[0][i+3];
        xmm3 = _mm_load_ps ((float*)tmp);   // previous element
        (*tmp)[0] = (float)src[1][i];
        (*tmp)[1] = (float)src[1][i+1];
        (*tmp)[2] = (float)src[1][i+2];
        (*tmp)[3] = (float)src[1][i+3];
        xmm4 = _mm_load_ps ((float*)tmp);   // current element

        for (int j=1; j<H-1; j++) {
            (*tmp)[0] = (float)src[j+1][i];
            (*tmp)[1] = (float)src[j+1][i+1];
            (*tmp)[2] = (float)src[j+1][i+2];
            (*tmp)[3] = (float)src[j+1][i+3];
            xmm5 = _mm_load_ps ((float*)tmp);   // next element
            // compute blured pixel
            xmm0 = _mm_mul_ps (xmm3, xmm2);
            xmm6 = _mm_mul_ps (xmm5, xmm2);
            xmm0 = _mm_add_ps (xmm6, xmm0);
            xmm6 = _mm_mul_ps (xmm4, xmm1);
            xmm0 = _mm_add_ps (xmm6, xmm0);
            // store blured pixel
            _mm_store_ps ((float*)&temp[j], xmm0);
            // update prev and current
            xmm3 = xmm4;
            xmm4 = xmm5;
        }
        for (int j=1; j<H-1; j++) {
            dst[j][i]   = (T)temp[j][0];
            dst[j][i+1] = (T)temp[j][1];
            dst[j][i+2] = (T)temp[j][2];
            dst[j][i+3] = (T)temp[j][3];
        }        
        dst[0][i] = src[0][i];
        dst[0][i+1] = src[0][i+1];
        dst[0][i+2] = src[0][i+2];
        dst[0][i+3] = src[0][i+3];
        dst[H-1][i] = src[H-1][i];
        dst[H-1][i+1] = src[H-1][i+1];
        dst[H-1][i+2] = src[H-1][i+2];
        dst[H-1][i+3] = src[H-1][i+3];
    }
//    int stride = dst[1] - dst[0];
    for (; i<col_to; i++) {
        for (int j = 1; j<H-1; j++) 
            buffer[j] = (T)(c1 * (src[j-1][i] + src[j+1][i]) + c0 * src[j][i]);
        dst[0][i] = src[0][i];
	for (int j=1; j<H-1; j++)
          dst[j][i] = buffer[j];
//        T* crow = &dst[1][i];
//        T* cbuff = &buffer[1];
//        for (int j = 1; j<H-1; j++, crow += stride, cbuff++)
//            (*crow) = (*cbuff);
        dst[H-1][i] = src[H-1][i];
    }
}
/*
template<class T> void gaussVertical3 (T** src, T** dst, T* buffer, int H, int col_from, int col_to, const float c0, const float c1) {

    time_t t1 = clock ();
    
    int stride = dst[1] - dst[0];
    for (int j=col_from; j<col_to; j++) {
        for (int i = 1; i<H-1; i++) 
            buffer[i] = (c1 * (src[i-1][j] + src[i+1][j]) + c0 * src[i][j]);
        dst[0][j] = src[0][j];
	for (int i=1; i<H-1; i++)
          dst[i][j] = buffer[i];
//        T* crow = &dst[1][j];
//        T* cbuff = &buffer[1];
//        for (int i = 1; i<H-1; i++, crow += stride, cbuff++)
//            (*crow) = (*cbuff);
        dst[H-1][j] = src[H-1][j];
    }
    
    printf ("vertical time = %d\n", clock()-t1);
}
*/
template<class T> void gaussHorizontal (T** src, T** dst, AlignedBuffer<float>* buffer, int W, int row_from, int row_to, double sigma) {

    if (sigma<0.25) {
        // dont perform filtering
        if (src!=dst)
            for (int i = row_from; i<row_to; i++) 
                memcpy (dst[i], src[i], W*sizeof(T));
        return;
    }

    if (sigma<0.6) {
        // compute 3x3 kernel
//        double c0 = 2.0 / exp (-1.0 / (2.0 * sigma * sigma));
//        printf ("c0=%g\n", c0);
//        gaussHorizontal3<T> (src, dst, (T*)(buffer->data), W, row_from, row_to, (int) round(c0), 2);
        double c1 = exp (-1.0 / (2.0 * sigma * sigma));
        double csum = 2.0 * c1 + 1.0;
        c1 /= csum;
        double c0 = 1.0 / csum;
        gaussHorizontal3<T> (src, dst, (T*)(buffer->data), W, row_from, row_to, c0, c1);
        return;
    }

    // horizontal
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
    
    // SSE optimized version:
    __m128 xmm1; xmm1 = _mm_load1_ps (&B);
    __m128 xmm2; xmm2 = _mm_load1_ps (&b1);
    __m128 xmm3; xmm3 = _mm_load1_ps (&b2);
    __m128 xmm4; xmm4 = _mm_load1_ps (&b3);
    __m128 xmm5;
    __m128 xmm6;
      
    typedef float pfloat[4];
    pfloat* temp = (pfloat*)buffer->data + 1;
    pfloat* tmp = (pfloat*)buffer->data;

    memset (temp, 0, W*sizeof(pfloat));

    int i;
    for (i=row_from; i<row_to-3; i+=4) {
        T* row1 = src[i];
        T* row2 = src[i+1];
        T* row3 = src[i+2];
        T* row4 = src[i+3];
        for (int j=0; j<3; j++) {
            (*tmp)[0] = (float)row1[3];
            (*tmp)[1] = (float)row2[3];
            (*tmp)[2] = (float)row3[3];
            (*tmp)[3] = (float)row4[3];
            xmm5 = _mm_load_ps ((float*)tmp);
            _mm_store_ps ((float*)&temp[j], xmm5);
        }
        for (int j=3; j<W; j++) {
            (*tmp)[0] = (float)row1[j];
            (*tmp)[1] = (float)row2[j];
            (*tmp)[2] = (float)row3[j];
            (*tmp)[3] = (float)row4[j];
            xmm5 = _mm_load_ps ((float*)tmp);
            xmm5 =_mm_mul_ps (xmm1, xmm5);
            xmm6 = _mm_load_ps ((float*)&temp[j-1]);
            xmm6 = _mm_mul_ps (xmm2, xmm6);
            xmm5 = _mm_add_ps (xmm6, xmm5);
            xmm6 = _mm_load_ps ((float*)&temp[j-2]);
            xmm6 = _mm_mul_ps (xmm6, xmm3);
            xmm5 = _mm_add_ps (xmm5, xmm6);
            xmm6 = _mm_load_ps ((float*)&temp[j-3]);
            xmm6 = _mm_mul_ps (xmm6, xmm4);
            xmm5 = _mm_add_ps (xmm5, xmm6);
            _mm_store_ps ((float*)&temp[j], xmm5);
        }
        for (int j=W-4; j>=0; j--) {
          xmm5 = _mm_load_ps ((float*)&temp[j]);
          xmm5 =_mm_mul_ps (xmm1, xmm5);
          xmm6 = _mm_load_ps ((float*)&temp[j+1]);
          xmm6 = _mm_mul_ps (xmm2, xmm6);
          xmm5 = _mm_add_ps (xmm6, xmm5);
          xmm6 = _mm_load_ps ((float*)&temp[j+2]);
          xmm6 = _mm_mul_ps (xmm6, xmm3);
          xmm5 = _mm_add_ps (xmm5, xmm6);
          xmm6 = _mm_load_ps ((float*)&temp[j+3]);
          xmm6 = _mm_mul_ps (xmm6, xmm4);
          xmm5 = _mm_add_ps (xmm5, xmm6);
          _mm_store_ps ((float*)&temp[j], xmm5);
        }
        for (int j=0; j<W; j++) {
          dst[i][j]   = (T)temp[j][0];
          dst[i+1][j] = (T)temp[j][1];
          dst[i+2][j] = (T)temp[j][2];
          dst[i+3][j] = (T)temp[j][3];
        }
    }
    // blur remaining rows
    float* temp2 = (float*)buffer->data;
    for (; i<row_to; i++) {
        for (int j=0; j<3; j++)
            temp2[j] = src[i][3];
        for (int j=3; j<W; j++)
            temp2[j] = B * src[i][j] + b1*temp2[j-1] + b2*temp2[j-2] + b3*temp2[j-3];
        for (int j=W-4; j>=0; j--)
            temp2[j] = B * temp2[j] + b1*temp2[j+1] + b2*temp2[j+2] + b3*temp2[j+3];
        for (int j=0; j<W; j++)
            dst[i][j] = (T)temp2[j];
    }
}

template<class T> void gaussVertical (T** src, T** dst, AlignedBuffer<float>* buffer, int H, int col_from, int col_to, double sigma) {

    if (sigma<0.25) {
        // dont perform filtering
        if (src!=dst)
            for (int i = 0; i<H; i++) 
                for (int j=col_from; j<col_to; j++)
                    dst[i][j] = src[i][j];
        return;
    }

    if (sigma<0.6) {
        // compute 3x3 kernel
        double c1 = exp (-1.0 / (2.0 * sigma * sigma));
        double csum = 2.0 * c1 + 1.0;
        c1 /= csum;
        double c0 = 1.0 / csum;
        gaussVertical3<T> (src, dst, (T*)(buffer->data), H, col_from, col_to, c0, c1);
//        double c0 = 2.0 / exp (-1.0 / (2.0 * sigma * sigma));
//        gaussVertical3<T> (src, dst, (T*)(buffer->data), H, col_from, col_to, (int) round(c0), 2);
        return;
    }

    // vertical
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
    
    // SSE optimized version:
    __m128 xmm1; xmm1 = _mm_load1_ps (&B);
    __m128 xmm2; xmm2 = _mm_load1_ps (&b1);
    __m128 xmm3; xmm3 = _mm_load1_ps (&b2);
    __m128 xmm4; xmm4 = _mm_load1_ps (&b3);
    __m128 xmm5;
    __m128 xmm6;
      
    typedef float pfloat[4];
    pfloat* temp = (pfloat*)buffer->data + 1;
    pfloat* tmp = (pfloat*)buffer->data;

    memset (temp, 0, H*sizeof(pfloat));

    int i;
    for (i=col_from; i<col_to-3; i+=4) {
        for (int j=0; j<3; j++) {
            (*tmp)[0] = (float)src[3][i];
            (*tmp)[1] = (float)src[3][i+1];
            (*tmp)[2] = (float)src[3][i+2];
            (*tmp)[3] = (float)src[3][i+3];
            xmm5 = _mm_load_ps ((float*)tmp);
            _mm_store_ps ((float*)&temp[j], xmm5);
        }
        for (int j=3; j<H; j++) {
            (*tmp)[0] = (float)src[j][i];
            (*tmp)[1] = (float)src[j][i+1];
            (*tmp)[2] = (float)src[j][i+2];
            (*tmp)[3] = (float)src[j][i+3];
            xmm5 = _mm_load_ps ((float*)tmp);
            xmm5 =_mm_mul_ps (xmm1, xmm5);
            xmm6 = _mm_load_ps ((float*)&temp[j-1]);
            xmm6 = _mm_mul_ps (xmm2, xmm6);
            xmm5 = _mm_add_ps (xmm6, xmm5);
            xmm6 = _mm_load_ps ((float*)&temp[j-2]);
            xmm6 = _mm_mul_ps (xmm6, xmm3);
            xmm5 = _mm_add_ps (xmm5, xmm6);
            xmm6 = _mm_load_ps ((float*)&temp[j-3]);
            xmm6 = _mm_mul_ps (xmm6, xmm4);
            xmm5 = _mm_add_ps (xmm5, xmm6);
            _mm_store_ps ((float*)&temp[j], xmm5);
        }
        for (int j=H-4; j>=0; j--) {
            xmm5 = _mm_load_ps ((float*)&temp[j]);
            xmm5 =_mm_mul_ps (xmm1, xmm5);
            xmm6 = _mm_load_ps ((float*)&temp[j+1]);
            xmm6 = _mm_mul_ps (xmm2, xmm6);
            xmm5 = _mm_add_ps (xmm6, xmm5);
            xmm6 = _mm_load_ps ((float*)&temp[j+2]);
            xmm6 = _mm_mul_ps (xmm6, xmm3);
            xmm5 = _mm_add_ps (xmm5, xmm6);
            xmm6 = _mm_load_ps ((float*)&temp[j+3]);
            xmm6 = _mm_mul_ps (xmm6, xmm4);
            xmm5 = _mm_add_ps (xmm5, xmm6);
            _mm_store_ps ((float*)&temp[j], xmm5);
        }
        for (int j=0; j<H; j++) {
            dst[j][i]   = (T)temp[j][0];
            dst[j][i+1] = (T)temp[j][1];
            dst[j][i+2] = (T)temp[j][2];
            dst[j][i+3] = (T)temp[j][3];
        }
    }
    // blur remaining columns
    float* temp2 = (float*)buffer->data;
    for (; i<col_to; i++) {
        for (int j=0; j<3; j++)
            temp2[j] = src[3][i];
        for (int j=3; j<H; j++)
            temp2[j] = B * src[j][i] + b1*temp2[j-1] + b2*temp2[j-2] + b3*temp2[j-3];
        for (int j=H-4; j>=0; j--)
            temp2[j] = B * temp2[j] + b1*temp2[j+1] + b2*temp2[j+2] + b3*temp2[j+3];
        for (int j=0; j<H; j++)
            dst[j][i] = (T)temp2[j];
    }
}

#else

template<class T> void gaussHorizontal3 (T** src, T** dst, T* buffer, int W, int row_from, int row_to, const float c0, const float c1) {
    
    for (int i=row_from; i<row_to; i++) {
        for (int j=1; j<W-1; j++)
            buffer[j] = (T)(c1 * (src[i][j-1] + src[i][j+1]) + c0 * src[i][j]);
        dst[i][0] = src[i][0];
        memcpy (dst[i]+1, buffer+1, (W-2)*sizeof(T));
        dst[i][W-1] = src[i][W-1];
    }
}

template<class T> void gaussVertical3 (T** src, T** dst, T* buffer, int H, int col_from, int col_to, const float c0, const float c1) {
    
    for (int i=col_from; i<col_to; i++) {
        for (int j = 1; j<H-1; j++) 
            buffer[j] = (T)(c1 * (src[j-1][i] + src[j+1][i]) + c0 * src[j][i]);
        dst[0][i] = src[0][i];
	    for (int j=1; j<H-1; j++)
            dst[j][i] = buffer[j];
        dst[H-1][i] = src[H-1][i];
    }
}

template<class T> void gaussHorizontal (T** src, T** dst, AlignedBuffer<double>* buffer, int W, int row_from, int row_to, double sigma) {

    if (sigma<0.25) {
        // dont perform filtering
        if (src!=dst)
            for (int i = row_from; i<row_to; i++) 
                memcpy (dst[i], src[i], W*sizeof(T));
        return;
    }

    if (sigma<0.6) {
        // compute 3x3 kernel
        double c1 = exp (-1.0 / (2.0 * sigma * sigma));
        double csum = 2.0 * c1 + 1.0;
        c1 /= csum;
        double c0 = 1.0 / csum;
        gaussHorizontal3<T> (src, dst, (T*)(buffer->data), W, row_from, row_to, c0, c1);
        return;
    }

    // horizontal
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
    
    double* temp2 = (double*)buffer->data;
    for (int i=row_from; i<row_to; i++) {
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
    }
}

template<class T> void gaussVertical (T** src, T** dst, AlignedBuffer<double>* buffer, int H, int col_from, int col_to, double sigma) {

    if (sigma<0.25) {
        // dont perform filtering
        if (src!=dst)
            for (int i = 0; i<H; i++) 
                for (int j=col_from; j<col_to; j++)
                    dst[i][j] = src[i][j];
        return;
    }

    if (sigma<0.6) {
        // compute 3x3 kernel
        double c1 = exp (-1.0 / (2.0 * sigma * sigma));
        double csum = 2.0 * c1 + 1.0;
        c1 /= csum;
        double c0 = 1.0 / csum;
        gaussVertical3<T> (src, dst, (T*)(buffer->data), H, col_from, col_to, c0, c1);
        return;
    }

    // vertical
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

    double* temp2 = (double*)buffer->data;
    for (int i=col_from; i<col_to; i++) {
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
    }
}   
#endif


void gaussHorizontal_unsigned (unsigned short** src, unsigned short** dst, AlignedBuffer<double>* buffer, int W, int row_from, int row_to, double sigma);
void gaussVertical_unsigned (unsigned short** src, unsigned short** dst, AlignedBuffer<double>* buffer, int H, int col_from, int col_to, double sigma);
void gaussHorizontal_signed (short** src, short** dst, AlignedBuffer<double>* buffer, int W, int row_from, int row_to, double sigma);
void gaussVertical_signed (short** src, short** dst, AlignedBuffer<double>* buffer, int H, int col_from, int col_to, double sigma);
void gaussHorizontal_float (float** src, float** dst, AlignedBuffer<double>* buffer, int W, int row_from, int row_to, double sigma);
void gaussVertical_float (float** src, float** dst, AlignedBuffer<double>* buffer, int H, int col_from, int col_to, double sigma);

#endif
