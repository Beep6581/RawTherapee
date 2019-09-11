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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
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

#include <iostream>
#include <queue>

#include "guidedfilter.h"
#include "improcfun.h"
#include "procparams.h"
#include "rt_algo.h"
#include "rt_math.h"

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


int get_dark_channel(const array2D<float> &R, const array2D<float> &G, const array2D<float> &B, array2D<float> &dst, int patchsize, const float ambient[3], bool clip, bool multithread)
{
    const int W = R.width();
    const int H = R.height();

#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; y += patchsize) {
        const int pH = min(y + patchsize, H);
        for (int x = 0; x < W; x += patchsize) {
            float val = RT_INFINITY_F;
            const int pW = min(x + patchsize, W);
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
                    val = min(val, r, g, b);
                }
            }
            if (clip) {
                val = LIM01(val);
            }
            for (int yy = y; yy < pH; ++yy) {
                std::fill(dst[yy] + x, dst[yy] + pW, val);
            }
        }
    }

    return (W / patchsize + ((W % patchsize) > 0)) *  (H / patchsize + ((H % patchsize) > 0));
}


float estimate_ambient_light(const array2D<float> &R, const array2D<float> &G, const array2D<float> &B, const array2D<float> &dark, int patchsize, int npatches, float ambient[3])
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
    
    float darklim = RT_INFINITY_F;
    {
        std::priority_queue<float> p;
        for (int y = 0; y < H; y += patchsize) {
            for (int x = 0; x < W; x += patchsize) {
                if (!OOG(dark[y][x], 1.f - 1e-5f)) {
                    p.push(dark[y][x]);
                }
            }
        }
        darklim = get_percentile(p, 0.95);
    }

    std::vector<std::pair<int, int>> patches;
    patches.reserve(npatches);

    for (int y = 0; y < H; y += patchsize) {
        for (int x = 0; x < W; x += patchsize) {
            if (dark[y][x] >= darklim && !OOG(dark[y][x], 1.f)) {
                patches.push_back(std::make_pair(x, y));
            }
        }
    }

    if (options.rtSettings.verbose) {
        std::cout << "dehaze: computing ambient light from " << patches.size()
                  << " patches" << std::endl;
    }

    float bright_lim = RT_INFINITY_F;
    {
        std::priority_queue<float> l;
        
        for (auto &p : patches) {
            const int pW = min(p.first+patchsize, W);
            const int pH = min(p.second+patchsize, H);
            
            for (int y = p.second; y < pH; ++y) {
                for (int x = p.first; x < pW; ++x) {
                    l.push(R[y][x] + G[y][x] + B[y][x]);
                }
            }
        }

        bright_lim = get_percentile(l, 0.95);
    }

    double rr = 0, gg = 0, bb = 0;
    int n = 0;
    for (auto &p : patches) {
        const int pW = min(p.first+patchsize, W);
        const int pH = min(p.second+patchsize, H);
            
        for (int y = p.second; y < pH; ++y) {
            for (int x = p.first; x < pW; ++x) {
                float r = R[y][x];
                float g = G[y][x];
                float b = B[y][x];
                if (r + g + b >= bright_lim) {
                    rr += r;
                    gg += g;
                    bb += b;
                    ++n;
                }
            }
        }
    }
    n = std::max(n, 1);
    ambient[0] = rr / n;
    ambient[1] = gg / n;
    ambient[2] = bb / n;

    // taken from darktable
    return darklim > 0 ? -1.125f * std::log(darklim) : std::log(std::numeric_limits<float>::max()) / 2;
}


void extract_channels(Imagefloat *img, array2D<float> &r, array2D<float> &g, array2D<float> &b, int radius, float epsilon, bool multithread)
{
    const int W = img->getWidth();
    const int H = img->getHeight();

    array2D<float> imgR(W, H, img->r.ptrs, ARRAY2D_BYREFERENCE);
    guidedFilter(imgR, imgR, r, radius, epsilon, multithread);

    array2D<float> imgG(W, H, img->g.ptrs, ARRAY2D_BYREFERENCE);
    guidedFilter(imgG, imgG, g, radius, epsilon, multithread);

    array2D<float> imgB(W, H, img->b.ptrs, ARRAY2D_BYREFERENCE);
    guidedFilter(imgB, imgB, b, radius, epsilon, multithread);
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

    int patchsize = max(int(5 / scale), 2);
    float ambient[3];
    array2D<float> &t_tilde = dark;
    float max_t = 0.f;

    {
        int npatches = 0;
        array2D<float> R(W, H);
        array2D<float> G(W, H);
        array2D<float> B(W, H);
        extract_channels(img, R, G, B, patchsize, 1e-1, multiThread);
    
        patchsize = max(max(W, H) / 600, 2);
        npatches = get_dark_channel(R, G, B, dark, patchsize, nullptr, false, multiThread);
        DEBUG_DUMP(dark);

        max_t = estimate_ambient_light(R, G, B, dark, patchsize, npatches, ambient);

        if (options.rtSettings.verbose) {
            std::cout << "dehaze: ambient light is "
                      << ambient[0] << ", " << ambient[1] << ", " << ambient[2]
                      << std::endl;
        }

        get_dark_channel(R, G, B, dark, patchsize, ambient, true, multiThread);
    }

    if (min(ambient[0], ambient[1], ambient[2]) < 0.01f) {
        if (options.rtSettings.verbose) {
            std::cout << "dehaze: no haze detected" << std::endl;
        }
        img->normalizeFloatTo65535();
        return; // probably no haze at all
    }

    DEBUG_DUMP(t_tilde);

#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            dark[y][x] = 1.f - strength * dark[y][x];
        }
    }

    const int radius = patchsize * 4;
    const float epsilon = 1e-5;
    array2D<float> &t = t_tilde;

    {
        array2D<float> guideB(W, H, img->b.ptrs, ARRAY2D_BYREFERENCE);
        guidedFilter(guideB, t_tilde, t, radius, epsilon, multiThread);
    }
        
    DEBUG_DUMP(t);

    if (options.rtSettings.verbose) {
        std::cout << "dehaze: max distance is " << max_t << std::endl;
    }

    float depth = -float(params->dehaze.depth) / 100.f;
    const float t0 = max(1e-3f, std::exp(depth * max_t));
    const float teps = 1e-3f;
#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            // ensure that the transmission is such that to avoid clipping...
            float rgb[3] = { img->r(y, x), img->g(y, x), img->b(y, x) };
            // ... t >= tl to avoid negative values
            float tl = 1.f - min(rgb[0]/ambient[0], rgb[1]/ambient[1], rgb[2]/ambient[2]);
            // ... t >= tu to avoid values > 1
            float tu = t0 - teps;
            for (int c = 0; c < 3; ++c) {
                if (ambient[c] < 1) {
                    tu = max(tu, (rgb[c] - ambient[c])/(1.f - ambient[c]));
                }
            }
            float mt = max(t[y][x], t0, tl + teps, tu + teps);
            if (params->dehaze.showDepthMap) {
                img->r(y, x) = img->g(y, x) = img->b(y, x) = LIM01(1.f - mt);
            } else {
                float r = (rgb[0] - ambient[0]) / mt + ambient[0];
                float g = (rgb[1] - ambient[1]) / mt + ambient[1];
                float b = (rgb[2] - ambient[2]) / mt + ambient[2];

                img->r(y, x) = r;
                img->g(y, x) = g;
                img->b(y, x) = b;
            }
        }
    }

    img->normalizeFloatTo65535();
}


} // namespace rtengine
