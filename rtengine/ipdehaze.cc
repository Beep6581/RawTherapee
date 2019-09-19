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
#include "rt_algo.h"
#include "rt_math.h"
#define BENCHMARK
#include "StopWatch.h"

extern Options options;

namespace rtengine {

namespace {

int get_dark_channel(const array2D<float> &R, const array2D<float> &G, const array2D<float> &B, array2D<float> &dst, int patchsize, const float ambient[3], bool clip, bool multithread, float strength)
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
            for (int xx = x; xx < pW; ++xx) {
                for (int yy = y; yy < pH; ++yy) {
                    val = min(val, R[yy][xx] / ambient[0], G[yy][xx] / ambient[1], B[yy][xx] / ambient[2]);
                }
            }
            val = 1.f - strength * LIM01(val);
            for (int yy = y; yy < pH; ++yy) {
                std::fill(dst[yy] + x, dst[yy] + pW, val);
            }
        }
    }

    return (W / patchsize + ((W % patchsize) > 0)) *  (H / patchsize + ((H % patchsize) > 0));
}

int get_dark_channel_downsized(const array2D<float> &R, const array2D<float> &G, const array2D<float> &B, array2D<float> &dst, int patchsize, bool multithread)
{
    const int W = R.width();
    const int H = R.height();

#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; y += patchsize) {
        int yy = y / patchsize;
        const int pH = min(y + patchsize, H);
        for (int x = 0, xx = 0; x < W; x += patchsize, ++xx) {
            float val = RT_INFINITY_F;
            const int pW = min(x + patchsize, W);
            for (int xp = x; xp < pW; ++xp) {
                for (int yp = y; yp < pH; ++yp) {
                    val = min(val, R[yp][xp], G[yp][xp], B[yp][xp]);
                }
            }
            dst[yy][xx] = val;
        }
    }

    return (W / patchsize + ((W % patchsize) > 0)) *  (H / patchsize + ((H % patchsize) > 0));
}


float estimate_ambient_light(const array2D<float> &R, const array2D<float> &G, const array2D<float> &B, const array2D<float> &dark, int patchsize, int npatches, float ambient[3])
{
    const int W = R.width();
    const int H = R.height();

    float darklim = RT_INFINITY_F;
    {
        std::vector<float> p;
        for (int y = 0, yy = 0; y < H; y += patchsize, ++yy) {
            for (int x = 0, xx = 0; x < W; x += patchsize, ++xx) {
                if (!OOG(dark[yy][xx], 1.f - 1e-5f)) {
                    p.push_back(dark[yy][xx]);
                }
            }
        }
        const int pos = p.size() * 0.95;
        std::nth_element(p.begin(), p.begin() + pos, p.end());
        darklim = p[pos];
    }

    std::vector<std::pair<int, int>> patches;
    patches.reserve(npatches);

    for (int y = 0, yy = 0; y < H; y += patchsize, ++yy) {
        for (int x = 0, xx = 0; x < W; x += patchsize, ++xx) {
            if (dark[yy][xx] >= darklim && !OOG(dark[yy][xx], 1.f)) {
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
        std::vector<float> l;
        l.reserve(patches.size() * patchsize * patchsize);
        
        for (const auto &p : patches) {
            const int pW = min(p.first + patchsize, W);
            const int pH = min(p.second + patchsize, H);
            
            for (int y = p.second; y < pH; ++y) {
                for (int x = p.first; x < pW; ++x) {
                    l.push_back(R[y][x] + G[y][x] + B[y][x]);
                }
            }
        }
        const int pos = l.size() * 0.95;
        std::nth_element(l.begin(), l.begin() + pos, l.end());
        bright_lim = l[pos];
    }

    double rr = 0, gg = 0, bb = 0;
    int n = 0;
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) reduction(+:rr,gg,bb,n)
#endif
    for (const auto &p : patches) {
        const int pW = min(p.first + patchsize, W);
        const int pH = min(p.second + patchsize, H);
            
        for (int y = p.second; y < pH; ++y) {
            for (int x = p.first; x < pW; ++x) {
                const float r = R[y][x];
                const float g = G[y][x];
                const float b = B[y][x];
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
    if (!params->dehaze.enabled || params->dehaze.strength == 0.0) {
        return;
    }
BENCHFUN
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
    float max_t = 0.f;

    {
        array2D<float> R(W, H);
        array2D<float> G(W, H);
        array2D<float> B(W, H);
        extract_channels(img, R, G, B, patchsize, 1e-1, multiThread);

        patchsize = max(max(W, H) / 600, 2);
        array2D<float> darkDownsized(W / patchsize + 1, H / patchsize + 1);
        const int npatches = get_dark_channel_downsized(R, G, B, darkDownsized, patchsize, multiThread);

        max_t = estimate_ambient_light(R, G, B, darkDownsized, patchsize, npatches, ambient);

        if (options.rtSettings.verbose) {
            std::cout << "dehaze: ambient light is "
                      << ambient[0] << ", " << ambient[1] << ", " << ambient[2]
                      << std::endl;
        }

        if (min(ambient[0], ambient[1], ambient[2]) < 0.01f) {
            if (options.rtSettings.verbose) {
                std::cout << "dehaze: no haze detected" << std::endl;
            }
            img->normalizeFloatTo65535();
            return; // probably no haze at all
        }

        get_dark_channel(R, G, B, dark, patchsize, ambient, true, multiThread, strength);
    }

    const int radius = patchsize * 4;
    constexpr float epsilon = 1e-5f;

    {
        array2D<float> guideB(W, H, img->b.ptrs, ARRAY2D_BYREFERENCE);
        guidedFilter(guideB, dark, dark, radius, epsilon, multiThread);
    }
        
    if (options.rtSettings.verbose) {
        std::cout << "dehaze: max distance is " << max_t << std::endl;
    }

    const float depth = -float(params->dehaze.depth) / 100.f;
    const float t0 = max(1e-3f, std::exp(depth * max_t));
    const float teps = 1e-3f;
#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < H; ++y) {
        int x = 0;
#ifdef __SSE2__
        const vfloat onev = F2V(1.f);
        const vfloat ambient0v = F2V(ambient[0]);
        const vfloat ambient1v = F2V(ambient[1]);
        const vfloat ambient2v = F2V(ambient[2]);
        const vfloat t0v = F2V(t0);
        const vfloat tepsv = F2V(teps);
        const vfloat c65535v = F2V(65535.f);
        for (; x < W - 3; x += 4) {
            // ensure that the transmission is such that to avoid clipping...
            vfloat r = LVFU(img->r(y, x));
            vfloat g = LVFU(img->g(y, x));
            vfloat b = LVFU(img->b(y, x));
            // ... t >= tl to avoid negative values
            const vfloat tlv = onev - vminf(r / ambient0v, vminf(g / ambient1v, b / ambient2v));
            // ... t >= tu to avoid values > 1
            r -= ambient0v;
            g -= ambient1v;
            b -= ambient2v;

            vfloat tuv = t0v - tepsv;
            tuv = vself(vmaskf_lt(ambient0v, onev), vmaxf(tuv, r / (onev - ambient0v)), tuv);
            tuv = vself(vmaskf_lt(ambient1v, onev), vmaxf(tuv, g / (onev - ambient1v)), tuv);
            tuv = vself(vmaskf_lt(ambient2v, onev), vmaxf(tuv, b / (onev - ambient2v)), tuv);

            const vfloat mtv = vmaxf(LVFU(dark[y][x]), vmaxf(tlv, tuv) + tepsv);
            if (params->dehaze.showDepthMap) {
                const vfloat valv = vclampf(onev - mtv, ZEROV, onev) * c65535v;
                STVFU(img->r(y, x), valv);
                STVFU(img->g(y, x), valv);
                STVFU(img->b(y, x), valv);
            } else {
                STVFU(img->r(y, x), (r / mtv + ambient0v) * c65535v);
                STVFU(img->g(y, x), (g / mtv + ambient1v) * c65535v);
                STVFU(img->b(y, x), (b / mtv + ambient2v) * c65535v);
            }
        }
#endif
        for (; x < W; ++x) {
            // ensure that the transmission is such that to avoid clipping...
            float r = img->r(y, x);
            float g = img->g(y, x);
            float b = img->b(y, x);
            // ... t >= tl to avoid negative values
            const float tl = 1.f - min(r / ambient[0], g / ambient[1], b / ambient[2]);
            // ... t >= tu to avoid values > 1
            r -= ambient[0];
            g -= ambient[1];
            b -= ambient[2];

            float tu = t0 - teps;
            tu = ambient[0] < 1.f ? max(tu, r / (1.f - ambient[0])) : tu;
            tu = ambient[1] < 1.f ? max(tu, g / (1.f - ambient[1])) : tu;
            tu = ambient[2] < 1.f ? max(tu, b / (1.f - ambient[2])) : tu;

            const float mt = max(dark[y][x], tl + teps, tu + teps);
            if (params->dehaze.showDepthMap) {
                img->r(y, x) = img->g(y, x) = img->b(y, x) = LIM01(1.f - mt) * 65535.f;
            } else {
                img->r(y, x) = (r / mt + ambient[0]) * 65535.f;
                img->g(y, x) = (g / mt + ambient[1]) * 65535.f;
                img->b(y, x) = (b / mt + ambient[2]) * 65535.f;
            }
        }
    }
}


} // namespace rtengine
