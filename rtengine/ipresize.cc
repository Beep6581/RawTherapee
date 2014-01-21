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

#include "improcfun.h"
#include "rt_math.h"

#ifdef _OPENMP
#include <omp.h>
#endif

//#define PROFILE

#ifdef PROFILE
#   include <iostream>
#endif

namespace rtengine {

static inline float Lanc(float x, float a)
{
    if (x * x < 1e-6f)
        return 1.0f;
    else if (x * x > a * a)
        return 0.0f;
    else {
        x = static_cast<float>(M_PI) * x;
        return sinf(x) * sinf(x / a) / (x * x / a);
    }
}

static void Lanczos(const Image16* src, Image16* dst, float scale)
{
    const float delta = 1.0f / scale;
    const float a = 3.0f;
    const float sc = min(scale, 1.0f);
    const int support = static_cast<int>(2.0f * a / sc) + 1;
    
    // storage for precomputed parameters for horisontal interpolation
    float * wwh = new float[support * dst->width];
    int * jj0 = new int[dst->width];
    int * jj1 = new int[dst->width];
    
    // temporal storage for vertically-interpolated row of pixels
    float * lr = new float[src->width];
    float * lg = new float[src->width];
    float * lb = new float[src->width];

    // Phase 1: precompute coefficients for horisontal interpolation
    
    for (int j = 0; j < dst->width; j++) {
    
        // x coord of the center of pixel on src image
        float x0 = (static_cast<float>(j) + 0.5f) * delta - 0.5f;

        // weights for interpolation in horisontal direction
        float * w = wwh + j * support;
        
        // sum of weights used for normalization
        float ws = 0.0f;

        jj0[j] = max(0, static_cast<int>(floorf(x0 - a / sc)) + 1);
        jj1[j] = min(src->width, static_cast<int>(floorf(x0 + a / sc)) + 1);

        // calculate weights
        for (int jj = jj0[j]; jj < jj1[j]; jj++) {
            int k = jj - jj0[j];
            float z = sc * (x0 - static_cast<float>(jj));
            w[k] = Lanc(z, a);
            ws += w[k];
        }
        
        // normalize weights
        for (int k = 0; k < support; k++) {
            w[k] /= ws;
        }
    }
    
    // Phase 2: do actual interpolation
    
    for (int i = 0; i < dst->height; i++) {
        
        // y coord of the center of pixel on src image
        float y0 = (static_cast<float>(i) + 0.5f) * delta - 0.5f;
        
        // weights for interpolation in y direction
        float w[support];
        
        // sum of weights used for normalization
        float ws= 0.0f;

        int ii0 = max(0, static_cast<int>(floorf(y0 - a / sc)) + 1);
        int ii1 = min(src->height, static_cast<int>(floorf(y0 + a / sc)) + 1);
        
        // calculate weights for vertical interpolation
        for (int ii = ii0; ii < ii1; ii++) {
            int k = ii - ii0;
            float z = sc * (y0 - static_cast<float>(ii));
            w[k] = Lanc(z, a);
            ws += w[k];
        }

        // normalize weights
        for (int k = 0; k < support; k++) {
            w[k] /= ws;
        }
        
        // Do vertical interpolation. Store results.
        for (int j = 0; j < src->width; j++) {
            
            float r = 0.0f, g = 0.0f, b = 0.0f;
            
            for (int ii = ii0; ii < ii1; ii++) {
                int k = ii - ii0;
            
                r += w[k] * src->r(ii,j);
                g += w[k] * src->g(ii,j);
                b += w[k] * src->b(ii,j);
            }
            
            lr[j] = r;
            lg[j] = g;
            lb[j] = b;
        }
        
        // Do horizontal interpolation
        for(int j = 0; j < dst->width; j++) {

            float * wh = wwh + support * j;
            
            float r = 0.0f, g = 0.0f, b = 0.0f;
            
            for (int jj = jj0[j]; jj < jj1[j]; jj++) {
                int k = jj - jj0[j];
            
                r += wh[k] * lr[jj];
                g += wh[k] * lg[jj];
                b += wh[k] * lb[jj];
            }
            
            dst->r(i,j) = CLIP(static_cast<int>(r));
            dst->g(i,j) = CLIP(static_cast<int>(g));
            dst->b(i,j) = CLIP(static_cast<int>(b));
        }
    }
    
    delete[] wwh;
    delete[] jj0;
    delete[] jj1;
    delete[] lr;
    delete[] lg;
    delete[] lb;
}

void ImProcFunctions::resize (Image16* src, Image16* dst, float dScale) {

#ifdef PROFILE
    time_t t1 = clock();
#endif

    if(params->resize.method == "Lanczos" ||
       params->resize.method == "Downscale (Better)" ||
       params->resize.method == "Downscale (Faster)"
      ) {
        Lanczos(src, dst, dScale);
    }
    else if (params->resize.method.substr(0,7)=="Bicubic") {
        float Av = -0.5f;
        if (params->resize.method=="Bicubic (Sharper)")
            Av = -0.75f;
        else if (params->resize.method=="Bicubic (Softer)")
            Av = -0.25f;
		#pragma omp parallel for if (multiThread)
        for (int i=0; i<dst->height; i++) {
            float wx[4], wy[4];
            float Dy = i / dScale;
            int yc  =  (int) Dy;
            Dy -= (float)yc;
            int ys = yc - 1; // smallest y-index used for interpolation
            // compute vertical weights
            float t1y = -Av*(Dy-1.0f)*Dy;
            float t2y = (3.0f - 2.0f*Dy)*Dy*Dy;
            wy[3] = t1y*Dy;
            wy[2] = t1y*(Dy - 1.0f) + t2y;
            wy[1] = -t1y*Dy + 1.0f - t2y;
            wy[0] = -t1y*(Dy - 1.0f);
            for (int j = 0; j < dst->width; j++) {
                float Dx = j / dScale;
                int xc  =  (int) Dx;
                Dx -= (float)xc;
                int xs = xc - 1; // smallest x-index used for interpolation
                if (ys >= 0 && ys < src->height-3 && xs >= 0 && xs <= src->width-3) {
                    // compute horizontal weights
                    float t1 = -Av*(Dx-1.0f)*Dx;
                    float t2 = (3.0f - 2.0f*Dx)*Dx*Dx;
                    wx[3] = t1*Dx;
                    wx[2] = t1*(Dx - 1.0f) + t2;
                    wx[1] = -t1*Dx + 1.0f - t2;
                    wx[0] = -t1*(Dx - 1.0f);
                    // compute weighted sum
                    int r = 0;
                    int g = 0;
                    int b = 0;
                    for (int x=0; x<4; x++)
                        for (int y=0; y<4; y++) {
                            float w = wx[x]*wy[y];
                            r += w*src->r(ys+y,xs+x);
                            g += w*src->g(ys+y,xs+x);
                            b += w*src->b(ys+y,xs+x);
                        }
                    dst->r(i,j) = CLIP(r);
                    dst->g(i,j) = CLIP(g);
                    dst->b(i,j) = CLIP(b);
                }
                else {
                    xc = LIM(xc, 0, src->width-1);
                    yc = LIM(yc, 0, src->height-1);
                    int nx = xc + 1;
                    if (nx >= src->width)
                        nx = xc;
                    int ny = yc + 1;
                    if (ny >= src->height)
                        ny = yc;
                    dst->r(i,j) = (1-Dx)*(1-Dy)*src->r(yc,xc) + (1-Dx)*Dy*src->r(ny,xc) + Dx*(1-Dy)*src->r(yc,nx) + Dx*Dy*src->r(ny,nx);
                    dst->g(i,j) = (1-Dx)*(1-Dy)*src->g(yc,xc) + (1-Dx)*Dy*src->g(ny,xc) + Dx*(1-Dy)*src->g(yc,nx) + Dx*Dy*src->g(ny,nx);
                    dst->b(i,j) = (1-Dx)*(1-Dy)*src->b(yc,xc) + (1-Dx)*Dy*src->b(ny,xc) + Dx*(1-Dy)*src->b(yc,nx) + Dx*Dy*src->b(ny,nx);
                }
            }
        }
    }
    else if (params->resize.method=="Bilinear") {
		#pragma omp parallel for if (multiThread)
        for (int i=0; i<dst->height; i++) {
            int sy = i/dScale;
            sy = LIM(sy, 0, src->height-1);
            float dy = i/dScale - sy;
            int ny = sy+1;
            if (ny>=src->height)
                ny = sy;
            for (int j=0; j<dst->width; j++) {
                int sx = j/dScale;
                sx = LIM(sx, 0, src->width-1);
                float dx = j/dScale - sx;
                int nx = sx+1;
                if (nx>=src->width)
                    nx = sx;
                dst->r(i,j) = (1-dx)*(1-dy)*src->r(sy,sx) + (1-dx)*dy*src->r(ny,sx) + dx*(1-dy)*src->r(sy,nx) + dx*dy*src->r(ny,nx);
                dst->g(i,j) = (1-dx)*(1-dy)*src->g(sy,sx) + (1-dx)*dy*src->g(ny,sx) + dx*(1-dy)*src->g(sy,nx) + dx*dy*src->g(ny,nx);
                dst->b(i,j) = (1-dx)*(1-dy)*src->b(sy,sx) + (1-dx)*dy*src->b(ny,sx) + dx*(1-dy)*src->b(sy,nx) + dx*dy*src->b(ny,nx);
            }
        }
    }
    else {
        // Nearest neighbour algorithm
		#pragma omp parallel for if (multiThread)
        for (int i=0; i<dst->height; i++) {
            int sy = i/dScale;
            sy = LIM(sy, 0, src->height-1);
            for (int j=0; j<dst->width; j++) {
                int sx = j/dScale;
                sx = LIM(sx, 0, src->width-1);
                dst->r(i,j) = src->r(sy,sx);
                dst->g(i,j) = src->g(sy,sx);
                dst->b(i,j) = src->b(sy,sx);
            }
        }
    }

#ifdef PROFILE    
    time_t t2 = clock();
    std::cout << "Resize: " << params->resize.method << ": "
        << (float)(t2 - t1) / CLOCKS_PER_SEC << std::endl;
#endif
}

}
