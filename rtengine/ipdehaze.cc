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


int get_dark_channel(const array2D<float> &R, const array2D<float> &G, const array2D<float> &B, array2D<float> &dst,
                     int patchsize, float *ambient, bool multithread)
{
    const int W = R.width();
    const int H = R.height();

    int npatches = 0;

#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; y += patchsize) {
        int pH = std::min(y+patchsize, H);
        for (int x = 0; x < W; x += patchsize, ++npatches) {
            float val = RT_INFINITY_F;
            int pW = std::min(x+patchsize, W);
            for (int yy = y; yy < pH; ++yy) {
                float yval = RT_INFINITY_F;
                for (int xx = x; xx < pW; ++xx) {
                    float r = R[yy][xx];
                    float g = G[yy][xx];
                    float b = B[yy][xx];
                    if (ambient) {
                        r /= ambient[0];
                        g /= ambient[1];
                        b /= ambient[2];
                    }
                    yval = min(yval, r, g, b);
                }
                val = min(val, yval);
            }
            val = LIM01(val);
            for (int yy = y; yy < pH; ++yy) {
                std::fill(dst[yy]+x, dst[yy]+pW, val);
            }
            float val2 = RT_INFINITY_F;
            for (int yy = y; yy < pH; ++yy) {
                for (int xx = x; xx < pW; ++xx) {
                    float r = R[yy][xx];
                    float g = G[yy][xx];
                    float b = B[yy][xx];
                    if (ambient) {
                        r /= ambient[0];
                        g /= ambient[1];
                        b /= ambient[2];
                    }
                    float l = min(r, g, b);
                    if (l >= 2.f * val) {
                        val2 = min(val2, l);
                        dst[yy][xx] = -1;
                    }
                }
            }
            if (val2 < RT_INFINITY_F) {
                val2 = LIM01(val2);
                for (int yy = y; yy < pH; ++yy) {
                    for (int xx = x; xx < pW; ++xx) {
                        if (dst[yy][xx] < 0.f) {
                            dst[yy][xx] = val2;
                        }
                    }
                }
            }
        }
    }

    return npatches;
}


int estimate_ambient_light(const array2D<float> &R, const array2D<float> &G, const array2D<float> &B, const array2D<float> &dark, const array2D<float> &Y, int patchsize, int npatches, float ambient[3])
{
    const int W = R.width();
    const int H = R.height();

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
                    float r = R[y][x];
                    float g = G[y][x];
                    float b = B[y][x];
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


void apply_contrast(array2D<float> &dark, float ambient, int contrast, double scale, bool multithread)
{
    if (contrast) {
        const int W = dark.width();
        const int H = dark.height();
        
        float avg = ambient * 0.25f;
        float c = contrast * 0.3f;

        std::vector<double> pts = {
            DCT_NURBS,
            0, //black point.  Value in [0 ; 1] range
            0, //black point.  Value in [0 ; 1] range

            avg - avg * (0.6 - c / 250.0), //toe point
            avg - avg * (0.6 + c / 250.0), //value at toe point

            avg + (1 - avg) * (0.6 - c / 250.0), //shoulder point
            avg + (1 - avg) * (0.6 + c / 250.0), //value at shoulder point

            1., // white point
            1. // value at white point
        };

        const DiagonalCurve curve(pts, CURVES_MIN_POLY_POINTS / scale);

#ifdef _OPENMP
        #pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                dark[y][x] = curve.getVal(dark[y][x]);
            }
        }
    }
}


void extract_channels(Imagefloat *img, const array2D<float> &Y, array2D<float> &r, array2D<float> &g, array2D<float> &b, int radius, float epsilon, bool multithread)
{
    const int W = img->getWidth();
    const int H = img->getHeight();

#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            r[y][x] = img->r(y, x);
            g[y][x] = img->g(y, x);
            b[y][x] = img->b(y, x);
        }
    }

    guidedFilter(Y, r, r, radius, epsilon, multithread);
    guidedFilter(Y, g, g, radius, epsilon, multithread);
    guidedFilter(Y, b, b, radius, epsilon, multithread);
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
    float strength = LIM01(float(params->dehaze.strength) / 100.f * 0.9f);

    if (options.rtSettings.verbose) {
        std::cout << "dehaze: strength = " << strength << std::endl;
    }

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
    array2D<float> Y(W, H);
    get_luminance(img, Y, ws, multiThread);

    array2D<float> R(W, H);
    array2D<float> G(W, H);
    array2D<float> B(W, H);
    int patchsize = max(int(20 / scale), 2);
    extract_channels(img, Y, R, G, B, patchsize, 1e-1, multiThread);
    
    array2D<float> dark(W, H);
    patchsize = std::max(W / (200 + params->dehaze.detail * (SGN(params->dehaze.detail) > 0 ? 4 : 1)), 2);
    int npatches = get_dark_channel(R, G, B, dark, patchsize, nullptr, multiThread);
    DEBUG_DUMP(dark);

    float ambient[3];
    int n = estimate_ambient_light(R, G, B, dark, Y, patchsize, npatches, ambient);
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
    get_dark_channel(R, G, B, dark, patchsize, ambient, multiThread);
    apply_contrast(dark, ambient_Y, params->dehaze.depth, scale, multiThread);
    DEBUG_DUMP(t_tilde);

    if (!params->dehaze.showDepthMap) {
#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                dark[y][x] = 1.f - strength * dark[y][x];
            }
        }
    }

    float mult = 2.f;
    if (params->dehaze.detail > 0) {
        mult -= (params->dehaze.detail / 100.f) * 1.9f;
    } else {
        mult -= params->dehaze.detail / 10.f;
    }
    const int radius = max(int(patchsize * mult), 1);
    const float epsilon = 2.5e-4;
    array2D<float> &t = t_tilde;

    if (!params->dehaze.showDepthMap)
    guidedFilter(Y, t_tilde, t, radius, epsilon, multiThread);

    DEBUG_DUMP(t);

    
    if (params->dehaze.showDepthMap) {
#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                img->r(y, x) = img->g(y, x) = img->b(y, x) = t[y][x] * 65535.f;
            }
        }
        return;
    }

    const float t0 = 0.1;
    const float teps = 1e-3;
#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            float rgb[3] = { img->r(y, x), img->g(y, x), img->b(y, x) };
            float tl = 1.f - min(rgb[0]/ambient[0], rgb[1]/ambient[1], rgb[2]/ambient[2]);
            float tu = t0 - teps;
            for (int c = 0; c < 3; ++c) {
                if (ambient[c] < 1) {
                    tu = max(tu, (rgb[c] - ambient[c])/(1.f - ambient[c]));
                }
            }
            float mt = max(t[y][x], t0, tl + teps, tu + teps);
            float r = (rgb[0] - ambient[0]) / mt + ambient[0];
            float g = (rgb[1] - ambient[1]) / mt + ambient[1];
            float b = (rgb[2] - ambient[2]) / mt + ambient[2];

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
    } else {
        img->normalizeFloatTo65535();
    }
}


} // namespace rtengine
