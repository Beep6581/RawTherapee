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
#include <omp.h>
#include <buffer.h>

namespace rtengine {

// classical filtering if the support window is small:

template<class T> void gaussHorizontal3 (Buffer<T>* src, Buffer<T>* dst, Dim size, T* buffer, const float c0, const float c1, bool multiThread) {

	int W = size.width;
    #pragma omp parallel for if (multiThread)
    for (int i=0; i<size.height; i++) {
    	T* temp = buffer + omp_get_thread_num() * W;
        for (int j=1; j<W-1; j++)
            temp[j] = (T)(c1 * (src->rows[i][j-1] + src->rows[i][j+1]) + c0 * src->rows[i][j]);
        dst->rows[i][0] = src->rows[i][0];
        memcpy (dst->rows[i]+1, temp+1, (W-2)*sizeof(T));
        dst->rows[i][W-1] = src->rows[i][W-1];
    }
}

template<class T> void gaussVertical3 (Buffer<T>* src, Buffer<T>* dst, Dim size, T* buffer, const float c0, const float c1, bool multiThread) {
    
    int H = size.height;
	#pragma omp parallel for if (multiThread)
    for (int i=0; i<size.width; i++) {
    	T* temp = buffer + omp_get_thread_num() * H;
        for (int j = 1; j<H-1; j++) 
        	temp[j] = (T)(c1 * (src->rows[j-1][i] + src->rows[j+1][i]) + c0 * src->rows[j][i]);
        dst->rows[0][i] = src->rows[0][i];
	    for (int j=1; j<H-1; j++)
            dst->rows[j][i] = temp[j];
        dst->rows[H-1][i] = src->rows[H-1][i];
    }
}

// fast gaussian approximation if the support window is large

template<class T> void gaussHorizontal (Buffer<T>* src, Buffer<T>* dst, Dim size, double* buffer, double sigma, bool multiThread) {

    if (sigma<0.25) {
        // dont perform filtering
        if (src!=dst)
        	for (int i=0; i<size.height; i++)
        		for (int j=0; j<size.width; j++)
        			dst->rows[i][j] = src->rows[i][j];
        return;
    }

    if (sigma<0.6) {
        // compute 3x3 kernel
        double c1 = exp (-1.0 / (2.0 * sigma * sigma));
        double csum = 2.0 * c1 + 1.0;
        c1 /= csum;
        double c0 = 1.0 / csum;
        gaussHorizontal3<T> (src, dst, size, (T*)buffer, c0, c1, multiThread);
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
    
    int W = size.width;
	#pragma omp parallel for if (multiThread)
    for (int i=0; i<size.height; i++) {

        double* temp2 = buffer + omp_get_thread_num() * W;

        temp2[0] = B * src->rows[i][0] + b1*src->rows[i][0] + b2*src->rows[i][0] + b3*src->rows[i][0];
        temp2[1] = B * src->rows[i][1] + b1*temp2[0]  + b2*src->rows[i][0] + b3*src->rows[i][0];
        temp2[2] = B * src->rows[i][2] + b1*temp2[1]  + b2*temp2[0]  + b3*src->rows[i][0];

        for (int j=3; j<W; j++)
            temp2[j] = B * src->rows[i][j] + b1*temp2[j-1] + b2*temp2[j-2] + b3*temp2[j-3];

        double temp2Wm1 = src->rows[i][W-1] + M[0][0]*(temp2[W-1] - src->rows[i][W-1]) + M[0][1]*(temp2[W-2] - src->rows[i][W-1]) + M[0][2]*(temp2[W-3] - src->rows[i][W-1]);
        double temp2W   = src->rows[i][W-1] + M[1][0]*(temp2[W-1] - src->rows[i][W-1]) + M[1][1]*(temp2[W-2] - src->rows[i][W-1]) + M[1][2]*(temp2[W-3] - src->rows[i][W-1]);
        double temp2Wp1 = src->rows[i][W-1] + M[2][0]*(temp2[W-1] - src->rows[i][W-1]) + M[2][1]*(temp2[W-2] - src->rows[i][W-1]) + M[2][2]*(temp2[W-3] - src->rows[i][W-1]);

        temp2[W-1] = temp2Wm1;
        temp2[W-2] = B * temp2[W-2] + b1*temp2[W-1] + b2*temp2W + b3*temp2Wp1;
        temp2[W-3] = B * temp2[W-3] + b1*temp2[W-2] + b2*temp2[W-1] + b3*temp2W;

        for (int j=W-4; j>=0; j--)
            temp2[j] = B * temp2[j] + b1*temp2[j+1] + b2*temp2[j+2] + b3*temp2[j+3];
        for (int j=0; j<W; j++)
            dst->rows[i][j] = (T)temp2[j];
    }
}

template<class T> void gaussVertical (Buffer<T>* src, Buffer<T>* dst, Dim size, double* buffer, double sigma, bool multiThread) {

    if (sigma<0.25) {
        // dont perform filtering
        if (src!=dst)
        	for (int i=0; i<size.height; i++)
        		for (int j=0; j<size.width; j++)
        			dst->rows[i][j] = src->rows[i][j];
        return;
    }

    if (sigma<0.6) {
        // compute 3x3 kernel
        double c1 = exp (-1.0 / (2.0 * sigma * sigma));
        double csum = 2.0 * c1 + 1.0;
        c1 /= csum;
        double c0 = 1.0 / csum;
        gaussVertical3<T> (src, dst, size, (T*)buffer, c0, c1, multiThread);
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

    int H = size.height;
    #pragma omp parallel for if (multiThread)
    for (int i=0; i<size.width; i++) {

        double* temp2 = buffer + omp_get_thread_num() * H;

    	temp2[0] = B * src->rows[0][i] + b1*src->rows[0][i] + b2*src->rows[0][i] + b3*src->rows[0][i];
        temp2[1] = B * src->rows[1][i] + b1*temp2[0]  + b2*src->rows[0][i] + b3*src->rows[0][i];
        temp2[2] = B * src->rows[2][i] + b1*temp2[1]  + b2*temp2[0]  + b3*src->rows[0][i];

        for (int j=3; j<H; j++)
            temp2[j] = B * src->rows[j][i] + b1*temp2[j-1] + b2*temp2[j-2] + b3*temp2[j-3];

        double temp2Hm1 = src->rows[H-1][i] + M[0][0]*(temp2[H-1] - src->rows[H-1][i]) + M[0][1]*(temp2[H-2] - src->rows[H-1][i]) + M[0][2]*(temp2[H-3] - src->rows[H-1][i]);
        double temp2H   = src->rows[H-1][i] + M[1][0]*(temp2[H-1] - src->rows[H-1][i]) + M[1][1]*(temp2[H-2] - src->rows[H-1][i]) + M[1][2]*(temp2[H-3] - src->rows[H-1][i]);
        double temp2Hp1 = src->rows[H-1][i] + M[2][0]*(temp2[H-1] - src->rows[H-1][i]) + M[2][1]*(temp2[H-2] - src->rows[H-1][i]) + M[2][2]*(temp2[H-3] - src->rows[H-1][i]);

        temp2[H-1] = temp2Hm1;
        temp2[H-2] = B * temp2[H-2] + b1*temp2[H-1] + b2*temp2H + b3*temp2Hp1;
        temp2[H-3] = B * temp2[H-3] + b1*temp2[H-2] + b2*temp2[H-1] + b3*temp2H;
        
        for (int j=H-4; j>=0; j--)
            temp2[j] = B * temp2[j] + b1*temp2[j+1] + b2*temp2[j+2] + b3*temp2[j+3];
        
        for (int j=0; j<H; j++)
            dst->rows[j][i] = (T)temp2[j];
    }
}

}
#endif
