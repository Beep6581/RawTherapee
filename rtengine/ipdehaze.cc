/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2018 Alberto Griggio <alberto.griggio@gmail.com>
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

/*
 * Haze removal using the algorithm described in the paper:
 *
 * Single Image Haze Removal Using Dark Channel Prior
 * by He, Sun and Tang
 *
 * using a guided filter for the "soft matting" of the transmission map
 *
 */  

#include "improcfun.h"
#include "guidedfilter.h"
#include "rt_math.h"
#include "rt_algo.h"
#include <iostream>
#include <queue>

extern Options options;

namespace rtengine {

namespace {

#if 0
#  define DEBUG_DUMP(arr)                                                 \
    do {                                                                \
        Imagefloat im(arr.width(), arr.height());                      \
        const char *out = "/tmp/" #arr ".tif";                     \
        for (int y = 0; y < im.getHeight(); ++y) {                      \
            for (int x = 0; x < im.getWidth(); ++x) {                   \
                im.r(y, x) = im.g(y, x) = im.b(y, x) = arr[y][x] * 65535.f; \
            }                                                           \
        }                                                               \
        im.saveTIFF(out, 16);                                           \
    } while (false)
#else
#  define DEBUG_DUMP(arr)
#endif


int get_dark_channel(const Imagefloat &src, array2D<float> &dst,
                     int patchsize, float *ambient, bool multithread)
{
    const int w = src.getWidth();
    const int h = src.getHeight();

    int npatches = 0;

#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < src.getHeight(); y += patchsize) {
        int pH = std::min(y+patchsize, h);
        for (int x = 0; x < src.getWidth(); x += patchsize, ++npatches) {
            float val = RT_INFINITY_F;
            int pW = std::min(x+patchsize, w);
            for (int yy = y; yy < pH; ++yy) {
                float yval = RT_INFINITY_F;
                for (int xx = x; xx < pW; ++xx) {
                    float r = src.r(yy, xx);
                    float g = src.g(yy, xx);
                    float b = src.b(yy, xx);
                    if (ambient) {
                        r /= ambient[0];
                        g /= ambient[1];
                        b /= ambient[2];
                    }
                    yval = min(yval, r, g, b);
                }
                val = min(val, yval);
            }
            for (int yy = y; yy < pH; ++yy) {
                std::fill(dst[yy]+x, dst[yy]+pW, val);
            }
            for (int yy = y; yy < pH; ++yy) {
                for (int xx = x; xx < pW; ++xx) {
                    float r = src.r(yy, xx);
                    float g = src.g(yy, xx);
                    float b = src.b(yy, xx);
                    if (ambient) {
                        r /= ambient[0];
                        g /= ambient[1];
                        b /= ambient[2];
                    }
                    float l = min(r, g, b);
                    if (l >= 2.f * val) {
                        dst[yy][xx] = l;
                    }
                }
            }
        }
    }

    return npatches;
}


int estimate_ambient_light(const Imagefloat *img, const array2D<float> &dark, const array2D<float> &Y, int patchsize, int npatches, float ambient[3])
{
    const int W = img->getWidth();
    const int H = img->getHeight();

    const auto get_percentile =
        [](std::priority_queue<float> &q, float prcnt) -> float
        {
            size_t n = LIM<size_t>(q.size() * prcnt, 1, q.size());
            while (q.size() > n) {
                q.pop();
            }
            return q.top();
        };
    
    float lim = RT_INFINITY_F;
    {
        std::priority_queue<float> p;
        for (int y = 0; y < H; y += patchsize) {
            for (int x = 0; x < W; x += patchsize) {
                p.push(dark[y][x]);
            }
        }
        lim = get_percentile(p, 0.95);
    }

    std::vector<std::pair<int, int>> patches;
    patches.reserve(npatches);

    for (int y = 0; y < H; y += patchsize) {
        for (int x = 0; x < W; x += patchsize) {
            if (dark[y][x] >= lim) {
                patches.push_back(std::make_pair(x, y));
            }
        }
    }

    if (options.rtSettings.verbose) {
        std::cout << "dehaze: computing ambient light from " << patches.size()
                  << " patches" << std::endl;
    }

    {
        std::priority_queue<float> l;
        
        for (auto &p : patches) {
            const int pW = std::min(p.first+patchsize, W);
            const int pH = std::min(p.second+patchsize, H);
            
            for (int y = p.second; y < pH; ++y) {
                for (int x = p.first; x < pW; ++x) {
                    l.push(Y[y][x]);
                }
            }
        }

        lim = get_percentile(l, 0.95);
    }

    double rr = 0, gg = 0, bb = 0;
    int n = 0;
    for (auto &p : patches) {
        const int pW = std::min(p.first+patchsize, W);
        const int pH = std::min(p.second+patchsize, H);
            
        for (int y = p.second; y < pH; ++y) {
            for (int x = p.first; x < pW; ++x) {
                if (Y[y][x] >= lim) {
                    float r = img->r(y, x);
                    float g = img->g(y, x);
                    float b = img->b(y, x);
                    rr += r;
                    gg += g;
                    bb += b;
                    ++n;
                }
            }
        }
    }
    ambient[0] = rr / n;
    ambient[1] = gg / n;
    ambient[2] = bb / n;

    return n;
}


void get_luminance(Imagefloat *img, array2D<float> &Y, TMatrix ws, bool multithread)
{
    const int W = img->getWidth();
    const int H = img->getHeight();
    
#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            Y[y][x] = Color::rgbLuminance(img->r(y, x), img->g(y, x), img->b(y, x), ws);
        }
    }
}


} // namespace


void ImProcFunctions::dehaze(Imagefloat *img)
{
    if (!params->dehaze.enabled) {
        return;
    }

    img->normalizeFloatTo1();
    
    const int W = img->getWidth();
    const int H = img->getHeight();
    const float strength = LIM01(float(params->dehaze.strength) / 100.f * 0.9f);

    if (options.rtSettings.verbose) {
        std::cout << "dehaze: strength = " << strength << std::endl;
    }
    
    array2D<float> dark(W, H);
    const int patchsize = std::max(W / 200, 2);
    int npatches = get_dark_channel(*img, dark, patchsize, nullptr, multiThread);
    DEBUG_DUMP(dark);

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
    array2D<float> Y(W, H);
    get_luminance(img, Y, ws, multiThread);
    
    float ambient[3];
    int n = estimate_ambient_light(img, dark, Y, patchsize, npatches, ambient);
    float ambient_Y = Color::rgbLuminance(ambient[0], ambient[1], ambient[2], ws);

    if (options.rtSettings.verbose) {
        std::cout << "dehaze: ambient light is "
                  << ambient[0] << ", " << ambient[1] << ", " << ambient[2]
                  << " (average of " << n << ")"
                  << std::endl;
        std::cout << "        ambient luminance is " << ambient_Y << std::endl;
    }

    if (min(ambient[0], ambient[1], ambient[2]) < 0.01f) {
        if (options.rtSettings.verbose) {
            std::cout << "dehaze: no haze detected" << std::endl;
        }
        img->normalizeFloatTo65535();
        return; // probably no haze at all
    }

    array2D<float> &t_tilde = dark;
    get_dark_channel(*img, dark, patchsize, ambient, multiThread);
    DEBUG_DUMP(t_tilde);
    
#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            dark[y][x] = 1.f - strength * dark[y][x];
        }
    }

    const int radius = patchsize * 2;
    const float epsilon = 2.5e-4;
    array2D<float> &t = t_tilde;
    
    guidedFilter(Y, t_tilde, t, radius, epsilon, multiThread);

    DEBUG_DUMP(t);

    const float t0 = 0.01;
#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            float mt = std::max(t[y][x], t0);
            float r = (img->r(y, x) - ambient[0]) / mt + ambient[0];
            float g = (img->g(y, x) - ambient[1]) / mt + ambient[1];
            float b = (img->b(y, x) - ambient[2]) / mt + ambient[2];
            img->r(y, x) = r;
            img->g(y, x) = g;
            img->b(y, x) = b;
        }
    }

    float oldmed;
    findMinMaxPercentile(Y, Y.width() * Y.height(), 0.5, oldmed, 0.5, oldmed, multiThread);

    get_luminance(img, Y, ws, multiThread);
    float newmed;

    findMinMaxPercentile(Y, Y.width() * Y.height(), 0.5, newmed, 0.5, newmed, multiThread);

    if (newmed > 1e-5f) {
        const float f1 = oldmed / newmed;
        const float f = f1 * 65535.f;
#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                float r = img->r(y, x);
                float g = img->g(y, x);
                float b = img->b(y, x);
                float h, s, l;
                Color::rgb2hslfloat(r * f, g * f, b * f, h, s, l);
                s = LIM01(s / f1);
                Color::hsl2rgbfloat(h, s, l, img->r(y, x), img->g(y, x), img->b(y, x));
            }
        }
    }
}


} // namespace rtengine
