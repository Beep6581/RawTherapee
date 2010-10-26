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
#include <rtengine.h>
#include <improcfun.h>
#include <glibmm.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>

namespace rtengine {

#undef CLIP
#undef CLIPTO
#undef CMAXVAL

#define CMAXVAL 0xffff
#define CLIP(a) ((a)>0?((a)<CMAXVAL?(a):CMAXVAL):0)
#define CLIPTO(a,b,c) ((a)>(b)?((a)<(c)?(a):(c)):(b))

void ImProcFunctions::resize (Image16* src, Image16* dst) {

    if(true) {
    //if(params->resize.method == "Lanczos") {
        double delta = 1.0 / params->resize.scale;
        const double a = 3.0;
        const int support = 6;
        const int kc = 2;
        
        Image16 * tmp = new Image16(src->width, dst->height);
        
        for (int i = 0; i < tmp->height; i++) {
            // y coord of the center of pixel on src image
            double y0 = (i + 0.5) * delta - 0.5;
            int i0 = floor(y0);
            
            // weights for interpolation in y direction
            double w[support];
            
            // sum of weights used for normalization
            double ww = 0.0;

            int ii0 = std::max(0, i0 - kc);
            int ii1 = std::min(src->height, i0 - kc + support);
            
            // calculate weights
            for (int ii = ii0; ii < ii1; ii++) {
                int k = ii - i0 + kc;
                double z = M_PI * (y0 - (i0 + k - kc));
                w[k] = sin(z) * sin(a*z) / (a * z * z);
                ww += w[k];
            }
            
            // normalize weights
            for (int k = 0; k < support; k++) {
                w[k] /= ww;
            }
            
            for (int j = 0; j < tmp->width; j++) {
                
                double r = 0.0, g = 0.0, b = 0.0;
                
                for (int ii = ii0; ii < ii1; ii++) {
                    int k = ii - i0 + kc;
                
                    r += w[k] * src->r[ii][j];
                    g += w[k] * src->g[ii][j];
                    b += w[k] * src->b[ii][j];
                }
                
                tmp->r[i][j] = CLIP((int)r);
                tmp->g[i][j] = CLIP((int)g);
                tmp->b[i][j] = CLIP((int)b);
            }
        }
        
        for (int j = 0; j < dst->width; j++) {
            // y coord of the center of pixel on src image
            double x0 = (j + 0.5) * delta - 0.5;
            int j0 = floor(x0);

            // weights for interpolation in y direction
            double w[support];
            
            // sum of weights used for normalization
            double ww = 0.0;

            int jj0 = std::max(0, j0 - kc);
            int jj1 = std::min(tmp->width, j0 - kc + support);

            // calculate weights
            for (int jj = jj0; jj < jj1; jj++) {
                int k = jj - j0 + kc;
                double z = M_PI * (x0 - (j0 + k - kc));
                w[k] = sin(z) * sin(a*z) / (a * z * z);
                ww += w[k];
            }
            
            // normalize weights
            for (int k = 0; k < support; k++) {
                w[k] /= ww;
            }
            
            for (int i = 0; i < dst->height; i++) {
                
                double r = 0.0, g = 0.0, b = 0.0;
                
                for (int jj = jj0; jj < jj1; jj++) {
                    int k = jj - j0 + kc;
                
                    r += w[k] * tmp->r[i][jj];
                    g += w[k] * tmp->g[i][jj];
                    b += w[k] * tmp->b[i][jj];
                }
                
                dst->r[i][j] = CLIP((int)r);
                dst->g[i][j] = CLIP((int)g);
                dst->b[i][j] = CLIP((int)b);
            }
        }
        
        delete tmp;
    }
    
    return;

	if(params->resize.method == "Downscale (Better)") {
        // small-scale algorithm by Ilia
        // provides much better quality on small scales
        // calculates mean value over source pixels which current destination pixel covers
        // works only for scales < 1
        // for scales ~1 it is analogous to bilinear
        // possibly, for even less scale factors (< 0.2 possibly) boundary pixels are not needed, omitting them can give a speedup
        // this algorithm is much slower on small factors than others, because it uses all pixels of the SOURCE image
        // Ilia Popov ilia_popov@rambler.ru 2010

        double delta = 1.0 / params->resize.scale;
        double k = params->resize.scale * params->resize.scale;

		#pragma omp parallel for if (multiThread)
        for(int i = 0; i < dst->height; i++) {
            // top and bottom boundary coordinates
            double y0 = i * delta;
            double y1 = (i + 1) * delta;

            int m0 = y0;
            m0 = CLIPTO(m0, 0, src->height-1);

            int m1 = y1;
            m1 = CLIPTO(m1, 0, src->height-1);

            // weights of boundary pixels
            double wy0 = 1.0 - (y0 - m0);
            double wy1 = y1 - m1;

            for(int j = 0; j < dst->width; j++) {
                // left and right boundary coordinates
                double x0 = j * delta;
                double x1 = (j + 1) * delta;

                int n0 = x0;
                n0 = CLIPTO(n0, 0, src->width-1);
                int n1 = x1;
                n1 = CLIPTO(n1, 0, src->width-1);

                double wx0 = 1.0 - (x0 - n0);
                double wx1 = x1 - n1;

                double r = 0;
                double g = 0;
                double b = 0;

                // integration
                // corners
                r += wy0 * wx0 * src->r[m0][n0] + wy0 * wx1 * src->r[m0][n1] + wy1 * wx0 * src->r[m1][n0] + wy1 * wx1 * src->r[m1][n1];
                g += wy0 * wx0 * src->g[m0][n0] + wy0 * wx1 * src->g[m0][n1] + wy1 * wx0 * src->g[m1][n0] + wy1 * wx1 * src->g[m1][n1];
                b += wy0 * wx0 * src->b[m0][n0] + wy0 * wx1 * src->b[m0][n1] + wy1 * wx0 * src->b[m1][n0] + wy1 * wx1 * src->b[m1][n1];

                // top and bottom boundaries
                for(int n = n0 + 1; n < n1; n++) {
                    r += wy0 * src->r[m0][n] + wy1 * src->r[m1][n];
                    g += wy0 * src->g[m0][n] + wy1 * src->g[m1][n];
                    b += wy0 * src->b[m0][n] + wy1 * src->b[m1][n];
                }

                // inner rows
                for(int m = m0 + 1; m < m1; m++) {
                    // left and right boundaries
                    r += wx0 * src->r[m][n0] + wx1 * src->r[m][n1];
                    g += wx0 * src->g[m][n0] + wx1 * src->g[m][n1];
                    b += wx0 * src->b[m][n0] + wx1 * src->b[m][n1];
                    // inner pixels
                    for(int n = n0 + 1; n < n1; n++) {
                        r += src->r[m][n];
                        g += src->g[m][n];
                        b += src->b[m][n];
                    }
                }
                // overall weight is equal to the DST pixel area in SRC coordinates
                r *= k;
                g *= k;
                b *= k;

                dst->r[i][j] = CLIP((int)r);
                dst->g[i][j] = CLIP((int)g);
                dst->b[i][j] = CLIP((int)b);
            }
        }
        return;
    }

    if(params->resize.method == "Downscale (Faster)")
	{
        // faster version of algo above, does not take into account border pixels,
        // which are summed with non-unity weights in slow algo. So, no need
        // for weights at all
        // Ilia Popov ilia_popov@rambler.ru 5.04.2010

        double delta = 1.0 / params->resize.scale;

        int p = (int) delta;

        // if actually we are doing upscaling, behave like Nearest
        if(p == 0)
            p = 1;

        int q = p/2;

        // may cause problems on 32-bit systems on extremely small factors.
        // In that case change 1024 to smth less
        const int divider = 1024;

        // scaling factor after summation
        int k = divider / (p * p);

		#pragma omp parallel for if (multiThread)
        for(int i = 0; i < dst->height; i++) {
            // y coordinate of center of destination pixel
            double y = (i + 0.5) * delta;

            int m0 = (int) (y) - q;
            m0 = CLIPTO(m0, 0, src->height-1);

            int m1 = m0 + p;
            if(m1 > src->height) {
                m1 = src->height;
                m0 = m1 - p;
            }
            m1 = CLIPTO(m1, 0, src->height);

            for(int j = 0; j < dst->width; j++) {
                // x coordinate of center of destination pixel
                double x = (j + 0.5) * delta;

                int n0 = (int) (x) - q;
                n0 = CLIPTO(n0, 0, src->width-1);

                int n1 = n0 + p;
                if(n1 > src->width) {
                    n1 = src->width;
                    n0 = n1 - p;
                }
                n1 = CLIPTO(n1, 0, src->width);

                int r = 0;
                int g = 0;
                int b = 0;

                // integration
                for(int m = m0; m < m1; m++) {
                    for(int n = n0; n < n1; n++) {
                        r += src->r[m][n];
                        g += src->g[m][n];
                        b += src->b[m][n];
                    }
                }
                dst->r[i][j] = CLIP( r * k / divider);
                dst->g[i][j] = CLIP( g * k / divider);
                dst->b[i][j] = CLIP( b * k / divider);
            }
        }
        return;
    }
    if (params->resize.method.substr(0,7)=="Bicubic") {
        double Av = -0.5;
        if (params->resize.method=="Bicubic (Sharper)")
            Av = -0.75;
        else if (params->resize.method=="Bicubic (Softer)")
            Av = -0.25;
		#pragma omp parallel for if (multiThread)
        for (int i=0; i<dst->height; i++) {
            double wx[4], wy[4];
            double Dy = i / params->resize.scale;
            int yc  =  (int) Dy; Dy -= (double)yc;
            int ys = yc - 1; // smallest y-index used for interpolation
            // compute vertical weights
            double t1y = -Av*(Dy-1.0)*Dy;
            double t2y = (3.0-2.0*Dy)*Dy*Dy;
            wy[3] = t1y*Dy;
            wy[2] = t1y*(Dy-1.0) + t2y;
            wy[1] = -t1y*Dy + 1.0 - t2y;
            wy[0] = -t1y*(Dy-1.0);
            for (int j=0; j<dst->width; j++) {
                double Dx = j / params->resize.scale;
                int xc  =  (int) Dx; Dx -= (double)xc;
                int xs = xc - 1; // smallest x-index used for interpolation
                if (ys >= 0 && ys <src->height-3 && xs >= 0 && xs <= src->width-3) {
                    // compute horizontal weights
                    double t1 = -Av*(Dx-1.0)*Dx;
                    double t2 = (3.0-2.0*Dx)*Dx*Dx;
                    wx[3] = t1*Dx;
                    wx[2] = t1*(Dx-1.0) + t2;
                    wx[1] = -t1*Dx + 1.0 - t2;
                    wx[0] = -t1*(Dx-1.0);
                    // compute weighted sum
                    int r = 0;
                    int g = 0;
                    int b = 0;
                    for (int x=0; x<4; x++)
                        for (int y=0; y<4; y++) {
                            double w = wx[x]*wy[y];
                            r += w*src->r[ys+y][xs+x];
                            g += w*src->g[ys+y][xs+x];
                            b += w*src->b[ys+y][xs+x];
                        }
                    dst->r[i][j] = CLIP(r);
                    dst->g[i][j] = CLIP(g);
                    dst->b[i][j] = CLIP(b);
                }
                else {
                    xc = CLIPTO(xc, 0, src->width-1);
                    yc = CLIPTO(yc, 0, src->height-1);
                    int nx = xc + 1;
                    if (nx>=src->width)
                        nx = xc;
                    int ny = yc + 1;
                    if (ny>=src->height)
                        ny = yc;
                    dst->r[i][j] = (1-Dx)*(1-Dy)*src->r[yc][xc] + (1-Dx)*Dy*src->r[ny][xc] + Dx*(1-Dy)*src->r[yc][nx] + Dx*Dy*src->r[ny][nx];
                    dst->g[i][j] = (1-Dx)*(1-Dy)*src->g[yc][xc] + (1-Dx)*Dy*src->g[ny][xc] + Dx*(1-Dy)*src->g[yc][nx] + Dx*Dy*src->g[ny][nx];
                    dst->b[i][j] = (1-Dx)*(1-Dy)*src->b[yc][xc] + (1-Dx)*Dy*src->b[ny][xc] + Dx*(1-Dy)*src->b[yc][nx] + Dx*Dy*src->b[ny][nx];
                }
            }
        }
    }
    else if (params->resize.method=="Bilinear") {
		#pragma omp parallel for if (multiThread)
        for (int i=0; i<dst->height; i++) {
            int sy = i/params->resize.scale;
            sy = CLIPTO(sy, 0, src->height-1);
            double dy = i/params->resize.scale - sy;
            int ny = sy+1;
            if (ny>=src->height)
                ny = sy;
            for (int j=0; j<dst->width; j++) {
                int sx = j/params->resize.scale;
                sx = CLIPTO(sx, 0, src->width-1);
                double dx = j/params->resize.scale - sx;
                int nx = sx+1;
                if (nx>=src->width)
                    nx = sx;
                dst->r[i][j] = (1-dx)*(1-dy)*src->r[sy][sx] + (1-dx)*dy*src->r[ny][sx] + dx*(1-dy)*src->r[sy][nx] + dx*dy*src->r[ny][nx];
                dst->g[i][j] = (1-dx)*(1-dy)*src->g[sy][sx] + (1-dx)*dy*src->g[ny][sx] + dx*(1-dy)*src->g[sy][nx] + dx*dy*src->g[ny][nx];
                dst->b[i][j] = (1-dx)*(1-dy)*src->b[sy][sx] + (1-dx)*dy*src->b[ny][sx] + dx*(1-dy)*src->b[sy][nx] + dx*dy*src->b[ny][nx];
            }
        }
    }
    else {
		#pragma omp parallel for if (multiThread)
        for (int i=0; i<dst->height; i++) {
            int sy = i/params->resize.scale;
            sy = CLIPTO(sy, 0, src->height-1);
            for (int j=0; j<dst->width; j++) {
                int sx = j/params->resize.scale;
                sx = CLIPTO(sx, 0, src->width-1);
                dst->r[i][j] = src->r[sy][sx];
                dst->g[i][j] = src->g[sy][sx];
                dst->b[i][j] = src->b[sy][sx];
            }
        }
    }
}

}
