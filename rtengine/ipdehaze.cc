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

#include <algorithm>
#include <iostream>
#include <vector>

#include "array2D.h"
#include "color.h"
#include "guidedfilter.h"
#include "iccstore.h"
#include "imagefloat.h"
#include "improcfun.h"
#include "procparams.h"
#include "rescale.h"
#include "rt_math.h"
//#define BENCHMARK
#include "StopWatch.h"

#include "rtgui/options.h"

namespace rtengine
{

namespace
{

float normalize(Imagefloat *rgb, bool multithread)
{
    float maxval = 0.f;
    const int W = rgb->getWidth();
    const int H = rgb->getHeight();
#ifdef _OPENMP
    #pragma omp parallel for reduction(max:maxval) schedule(dynamic, 16) if (multithread)
#endif

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            maxval = max(maxval, rgb->r(y, x), rgb->g(y, x), rgb->b(y, x));
        }
    }

    maxval = max(maxval * 2.f, 65535.f);
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16) if (multithread)
#endif

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            rgb->r(y, x) /= maxval;
            rgb->g(y, x) /= maxval;
            rgb->b(y, x) /= maxval;
        }
    }

    return maxval;
}

void restore(Imagefloat *rgb, float maxval, bool multithread)
{
    const int W = rgb->getWidth();
    const int H = rgb->getHeight();

    if (maxval > 0.f && maxval != 1.f) {
#ifdef _OPENMP
        #       pragma omp parallel for if (multithread)
#endif

        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                rgb->r(y, x) *= maxval;
                rgb->g(y, x) *= maxval;
                rgb->b(y, x) *= maxval;
            }
        }
    }
}

int get_dark_channel(const array2D<float> &R, const array2D<float> &G, const array2D<float> &B, array2D<float> &dst, int patchsize, const float ambient[3], bool clip, bool multithread, float strength)
{
    const int W = R.getWidth();
    const int H = R.getHeight();

#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif

    for (int y = 0; y < H; y += patchsize) {
        const int pH = min(y + patchsize, H);

        for (int x = 0; x < W; x += patchsize) {
            float minR = RT_INFINITY_F;
            float minG = RT_INFINITY_F;
            float minB = RT_INFINITY_F;
#ifdef __SSE2__
            vfloat minRv = F2V(minR);
            vfloat minGv = F2V(minG);
            vfloat minBv = F2V(minB);
#endif
            const int pW = min(x + patchsize, W);

            for (int yy = y; yy < pH; ++yy) {
                int xx = x;
#ifdef __SSE2__

                for (; xx < pW - 3; xx += 4) {
                    minRv = vminf(minRv, LVFU(R[yy][xx]));
                    minGv = vminf(minGv, LVFU(G[yy][xx]));
                    minBv = vminf(minBv, LVFU(B[yy][xx]));
                }

#endif

                for (; xx < pW; ++xx) {
                    minR = min(minR, R[yy][xx]);
                    minG = min(minG, G[yy][xx]);
                    minB = min(minB, B[yy][xx]);
                }
            }

#ifdef __SSE2__
            minR = min(minR, vhmin(minRv));
            minG = min(minG, vhmin(minGv));
            minB = min(minB, vhmin(minBv));
#endif
            float val = min(minR / ambient[0], minG / ambient[1], minB / ambient[2]);
            val = 1.f - strength * LIM01(val);

            for (int yy = y; yy < pH; ++yy) {
                std::fill(dst[yy] + x, dst[yy] + pW, val);
            }
        }
    }

    return (W / patchsize + ((W % patchsize) > 0)) * (H / patchsize + ((H % patchsize) > 0));
}

int get_dark_channel_downsized(const array2D<float> &R, const array2D<float> &G, const array2D<float> &B, array2D<float> &dst, int patchsize, bool multithread)
{
    const int W = R.getWidth();
    const int H = R.getHeight();

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
                    val = min(val, R[yy][xx], G[yy][xx], B[yy][xx]);
                }
            }

            for (int yy = y; yy < pH; ++yy) {
                std::fill(dst[yy] + x, dst[yy] + pW, val);
            }
        }
    }

    return (W / patchsize + ((W % patchsize) > 0)) * (H / patchsize + ((H % patchsize) > 0));
}

float estimate_ambient_light(const array2D<float> &R, const array2D<float> &G, const array2D<float> &B, const array2D<float> &dark, int patchsize, int npatches, float ambient[3])
{
    const int W = R.getWidth();
    const int H = R.getHeight();

    float darklim = RT_INFINITY_F;
    {
        std::vector<float> p;

        for (int y = 0; y < H; y += patchsize) {
            for (int x = 0; x < W; x += patchsize) {
                if (!OOG(dark[y][x], 1.f - 1e-5f)) {
                    p.push_back(dark[y][x]);
                }
            }
        }

        if (p.empty()) {
            return 0.f;
        }

        const int pos = p.size() * 0.95;
        std::nth_element(p.begin(), p.begin() + pos, p.end());
        darklim = p[pos];
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

    if (settings->verbose) {
        std::cout << "dehaze: computing ambient light from " << patches.size()
                  << " patches" << std::endl;
    }

    float bright_lim = RT_INFINITY_F;
    {
        std::vector<float> l;
        l.reserve(patches.size() * patchsize * patchsize);

        for (auto &p : patches) {
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

    for (auto &p : patches) {
        const int pW = min(p.first + patchsize, W);
        const int pH = min(p.second + patchsize, H);

        for (int y = p.second; y < pH; ++y) {
            for (int x = p.first; x < pW; ++x) {
                float r = R[y][x];
                float g = G[y][x];
                float b = B[y][x];

                if (r + g + b >= bright_lim) {
                    rr += static_cast<double>(r);
                    gg += static_cast<double>(g);
                    bb += static_cast<double>(b);
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

void ImProcFunctions::dehaze(Imagefloat *img, const DehazeParams &dehazeParams)
{
    if (!dehazeParams.enabled || dehazeParams.strength == 0.0) {
        return;
    }

    const float maxChannel = normalize(img, multiThread);

    const int W = img->getWidth();
    const int H = img->getHeight();
    const float strength = LIM01(float(dehazeParams.strength) / 100.f * 0.9f);

    array2D<float> dark(W, H);

    int patchsize = max(int(5 / scale), 2);
    float ambient[3];
    float maxDistance = 0.f;

    {
        array2D<float>& R = dark; // R and dark can safely use the same buffer, which is faster and reduces memory allocations/deallocations
        array2D<float> G(W, H);
        array2D<float> B(W, H);
        extract_channels(img, R, G, B, patchsize, 1e-1, multiThread);

        {
            constexpr int sizecap = 200;
            const float r = static_cast<float>(W) / static_cast<float>(H);
            const int hh = r >= 1.f ? sizecap : sizecap / r;
            const int ww = r >= 1.f ? sizecap * r : sizecap;

            if (W <= ww && H <= hh) {
                // don't rescale small thumbs
                array2D<float> D(W, H);
                const int npatches = get_dark_channel_downsized(R, G, B, D, 2, multiThread);
                maxDistance = estimate_ambient_light(R, G, B, D, patchsize, npatches, ambient);
            } else {
                array2D<float> RR(ww, hh);
                array2D<float> GG(ww, hh);
                array2D<float> BB(ww, hh);
                rescaleNearest(R, RR, multiThread);
                rescaleNearest(G, GG, multiThread);
                rescaleNearest(B, BB, multiThread);
                array2D<float> D(ww, hh);

                const int npatches = get_dark_channel_downsized(RR, GG, BB, D, 2, multiThread);
                maxDistance = estimate_ambient_light(RR, GG, BB, D, patchsize, npatches, ambient);
            }
        }

        if (min(ambient[0], ambient[1], ambient[2]) < 0.01f) {
            if (settings->verbose) {
                std::cout << "dehaze: no haze detected" << std::endl;
            }
            restore(img, maxChannel, multiThread);
            return; // probably no haze at all
        }
        patchsize = max(max(W, H) / 600, 2);

        if (settings->verbose) {
            std::cout << "dehaze: ambient light is "
                      << ambient[0] << ", " << ambient[1] << ", " << ambient[2]
                      << std::endl;
        }

        get_dark_channel(R, G, B, dark, patchsize, ambient, true, multiThread, strength);
    }

    const int radius = patchsize * 4;
    constexpr float epsilon = 1e-5f;

    array2D<float> guideB(W, H, img->b.ptrs, ARRAY2D_BYREFERENCE);
    guidedFilter(guideB, dark, dark, radius, epsilon, multiThread);
        
    if (settings->verbose) {
        std::cout << "dehaze: max distance is " << maxDistance << std::endl;
    }

    const float depth = -float(dehazeParams.depth) / 100.f;
    const float t0 = max(1e-3f, std::exp(depth * maxDistance));
    constexpr float teps = 1.f + 1e-3f;

    const float satBlend = dehazeParams.saturation / 100.f;
    const TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
#ifdef __SSE2__
    const vfloat wsv[3] = {F2V(ws[1][0]), F2V(ws[1][1]),F2V(ws[1][2])};
#endif
    const float ambientY = Color::rgbLuminance(ambient[0], ambient[1], ambient[2], ws);

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
        const vfloat ambientYv = F2V(ambientY);
        const vfloat epsYv = F2V(1e-5f);
        const vfloat t0v = F2V(t0);
        const vfloat tepsv = F2V(teps);
        const vfloat cmaxChannelv = F2V(maxChannel);
        const vfloat satBlendv = F2V(satBlend);
        for (; x < W - 3; x += 4) {
            // ensure that the transmission is such that to avoid clipping...
            const vfloat r = LVFU(img->r(y, x));
            const vfloat g = LVFU(img->g(y, x));
            const vfloat b = LVFU(img->b(y, x));
            // ... t >= tl to avoid negative values
            const vfloat tlv = tepsv - vminf(r / ambient0v, vminf(g / ambient1v, b / ambient2v));
            const vfloat mtv = vmaxf(LVFU(dark[y][x]), vmaxf(tlv, t0v));
            if (dehazeParams.showDepthMap) {
                const vfloat valv = vclampf(onev - mtv, ZEROV, onev) * cmaxChannelv;
                STVFU(img->r(y, x), valv);
                STVFU(img->g(y, x), valv);
                STVFU(img->b(y, x), valv);
            } else {
                const vfloat Yv = Color::rgbLuminance(r, g, b, wsv);
                const vfloat YYv = (Yv - ambientYv) / mtv + ambientYv;
                const vfloat fv = vself(vmaskf_gt(Yv, epsYv), cmaxChannelv * YYv / Yv, cmaxChannelv);
                STVFU(img->r(y, x), vintpf(satBlendv, ((r - ambient0v) / mtv + ambient0v) * cmaxChannelv, r * fv));
                STVFU(img->g(y, x), vintpf(satBlendv, ((g - ambient1v) / mtv + ambient1v) * cmaxChannelv, g * fv));
                STVFU(img->b(y, x), vintpf(satBlendv, ((b - ambient2v) / mtv + ambient2v) * cmaxChannelv, b * fv));
            }
        }
#endif
        for (; x < W; ++x) {
            // ensure that the transmission is such that to avoid clipping...
            const float r = img->r(y, x);
            const float g = img->g(y, x);
            const float b = img->b(y, x);
            // ... t >= tl to avoid negative values
            const float tl = teps - min(r / ambient[0], g / ambient[1], b / ambient[2]);
            const float mt = max(dark[y][x], t0, tl);
            if (dehazeParams.showDepthMap) {
                img->r(y, x) = img->g(y, x) = img->b(y, x) = LIM01(1.f - mt) * maxChannel;
            } else {
                const float Y = Color::rgbLuminance(img->r(y, x), img->g(y, x), img->b(y, x), ws);
                const float YY = (Y - ambientY) / mt + ambientY;
                const float f = Y > 1e-5f ? maxChannel * YY / Y : maxChannel;
                img->r(y, x) = intp(satBlend, ((r - ambient[0]) / mt + ambient[0]) * maxChannel, r * f);
                img->g(y, x) = intp(satBlend, ((g - ambient[1]) / mt + ambient[1]) * maxChannel, g * f);
                img->b(y, x) = intp(satBlend, ((b - ambient[2]) / mt + ambient[2]) * maxChannel, b * f);
            }
        }
    }
}

void ImProcFunctions::dehazeloc(Imagefloat *img, const DehazeParams &dehazeParams, int sk, int sp)
{
    //J.Desmis 12 2019 - this version derived from ART, is slower than the main from maximum 10% - probably use of SSE
    //Probably Ingo could solved this problem in some times

    if (!dehazeParams.enabled || dehazeParams.strength == 0.0) {
        return;
    }

    const float maxChannel = normalize(img, multiThread);

    const int W = img->getWidth();
    const int H = img->getHeight();
    const float strength = LIM01(float(std::abs(dehazeParams.strength)) / 100.f * 0.9f);
    const bool add_haze = dehazeParams.strength < 0;

    array2D<float> dark(W, H);

    int patchsize = max(int(5 / scale), 2);
    float ambient[3];
    float maxDistance = 0.f;

    int whit = 0;
    int blac = params->locallab.spots.at(sp).dehazeblack;

    if(blac != 0) {
        ImProcFunctions::tone_eqdehaz(this, img, whit, blac, params->icm.workingProfile, sk, multiThread);
    }

    {
        array2D<float>& R = dark; // R and dark can safely use the same buffer, which is faster and reduces memory allocations/deallocations
        array2D<float> G(W, H);
        array2D<float> B(W, H);
        extract_channels(img, R, G, B, patchsize, 1e-1, multiThread);

        {
            constexpr int sizecap = 200;
            const float r = static_cast<float>(W) / static_cast<float>(H);
            const int hh = r >= 1.f ? sizecap : sizecap / r;
            const int ww = r >= 1.f ? sizecap * r : sizecap;

            if (W <= ww && H <= hh) {
                // don't rescale small thumbs
                array2D<float> D(W, H);
                const int npatches = get_dark_channel_downsized(R, G, B, D, 2, multiThread);
                maxDistance = estimate_ambient_light(R, G, B, D, patchsize, npatches, ambient);
            } else {
                array2D<float> RR(ww, hh);
                array2D<float> GG(ww, hh);
                array2D<float> BB(ww, hh);
                rescaleNearest(R, RR, multiThread);
                rescaleNearest(G, GG, multiThread);
                rescaleNearest(B, BB, multiThread);
                array2D<float> D(ww, hh);

                const int npatches = get_dark_channel_downsized(RR, GG, BB, D, 2, multiThread);
                maxDistance = estimate_ambient_light(RR, GG, BB, D, patchsize, npatches, ambient);
            }
        }

        if (min(ambient[0], ambient[1], ambient[2]) < 0.01f) {
            if (settings->verbose) {
                std::cout << "dehaze: no haze detected" << std::endl;
            }

            restore(img, maxChannel, multiThread);
            return; // probably no haze at all
        }

        patchsize = max(max(W, H) / 600, 2);

        if (settings->verbose) {
            std::cout << "dehaze: ambient light is "
                      << ambient[0] << ", " << ambient[1] << ", " << ambient[2]
                      << std::endl;
        }

        get_dark_channel(R, G, B, dark, patchsize, ambient, true, multiThread, strength);
    }


    const int radius = patchsize * 4;
    constexpr float epsilon = 1e-5f;

    array2D<float> guideB(W, H, img->b.ptrs, ARRAY2D_BYREFERENCE);
    guidedFilter(guideB, dark, dark, radius, epsilon, multiThread);

    if (settings->verbose) {
        std::cout << "dehaze: max distance is " << maxDistance << std::endl;
    }

    const float depth = -float(dehazeParams.depth) / 100.f;
    constexpr float teps = 1e-6f;
    const float t0 = max(teps, std::exp(depth * maxDistance));

    const float satBlend = dehazeParams.saturation / 100.f;

    const TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);

    const float ambientY = Color::rgbLuminance(ambient[0], ambient[1], ambient[2], ws);
#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            // ensure that the transmission is such that to avoid clipping...
            const float rIn = img->r(y, x);
            const float gIn = img->g(y, x);
            const float bIn = img->b(y, x);
            // ... t >= tl to avoid negative values
            const float tl = 1.f + teps - min(rIn / ambient[0], gIn / ambient[1], bIn / ambient[2]);
            const float mt = max(dark[y][x], t0, tl);

            if (dehazeParams.showDepthMap) {
                img->r(y, x) = img->g(y, x) = img->b(y, x) = LIM01(1.f - mt) * maxChannel;
            } else {
                float f = 1.f;
                const float Y = Color::rgbLuminance(rIn, gIn, bIn, ws);
                if (Y > 1e-5f) {
                    float YY = (Y - ambientY) / mt + ambientY;
                    if (add_haze) {
                        YY = Y + Y - YY;
                    }
                    f = YY / Y;
                }
                const float r1 = rIn * f;
                const float g1 = gIn * f;
                const float b1 = bIn * f;

                float r2 = ((rIn - ambient[0]) / mt + ambient[0]);
                float g2 = ((gIn - ambient[1]) / mt + ambient[1]);
                float b2 = ((bIn - ambient[2]) / mt + ambient[2]);

                if (add_haze) {
                    r2 = rIn + rIn - r2;
                    g2 = gIn + gIn - g2;
                    b2 = bIn + bIn - b2;
                }
                img->r(y, x) = intp(satBlend, r2, r1);
                img->g(y, x) = intp(satBlend, g2, g1);
                img->b(y, x) = intp(satBlend, b2, b1);
            }
        }
    }

    restore(img, maxChannel, multiThread);

}

} // namespace rtengine
