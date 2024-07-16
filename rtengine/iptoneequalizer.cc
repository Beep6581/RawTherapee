#include "color.h"
#include "guidedfilter.h"
#include "iccstore.h"
#include "imagefloat.h"
#include "improcfun.h"
#include "sleef.h"
#include "StopWatch.h"


namespace
{

const std::vector<std::array<float, 3>> colormap = {
    {0.5f, 0.f, 0.5f},
    {0.5f, 0.f, 0.5f},
    {0.5f, 0.f, 0.5f},
    {0.5f, 0.f, 0.5f},
    {0.5f, 0.f, 0.5f}, // blacks
    {0.f, 0.f, 1.f}, // shadows
    {0.5f, 0.5f, 0.5f}, // midtones
    {1.f, 1.f, 0.f}, // highlights
    {1.f, 0.f, 0.f}, // whites
    {1.f, 0.f, 0.f},
    {1.f, 0.f, 0.f},
    {1.f, 0.f, 0.f},
    {1.f, 0.f, 0.f},
    {1.f, 0.f, 0.f},
    {1.f, 0.f, 0.f}

};


void toneEqualizer(
    array2D<float> &R, array2D<float> &G, array2D<float> &B,
    const rtengine::ToneEqualizerParams &params,
    const Glib::ustring &workingProfile,
    double scale,
    bool multithread)
// adapted from the tone equalizer of darktable
/*
    Copyright 2019 Alberto Griggio <alberto.griggio@gmail.com>
    Small adaptation to Local Adjustment 10 2019 Jacques Desmis <jdesmis@gmail.com>
    This file is part of darktable,
    copyright (c) 2018 Aurelien Pierre.

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/

{
   // BENCHFUN

    const int W = R.getWidth();
    const int H = R.getHeight();
    array2D<float> Y(W, H);

    const auto log2 =
    [](float x) -> float {
        static const float l2 = xlogf(2);
        return xlogf(x) / l2;
    };

    const auto exp2 =
    [](float x) -> float {
        return pow_F(2.f, x);
    };
    // Build the luma channels: band-pass filters with gaussian windows of
    // std 2 EV, spaced by 2 EV
    const float centers[15] = {
        -16.0f, -14.0f, -12.0f, -10.0f, -8.0f, -6.0f,
        -4.0f, -2.0f, 0.0f, 2.0f, 4.0f, 6.0f, 8.0f, 10.0f, 12.0f
    };

    const auto conv = [&](int v, float lo, float hi) -> float {
        const float f = v < 0 ? lo : hi;
        return exp2(float(v) / 100.f * f);
    };
    const float factors[15] = {
        conv(params.bands[0], 2.f, 3.f), // -16 EV
        conv(params.bands[0], 2.f, 3.f), // -14 EV
        conv(params.bands[0], 2.f, 3.f), // -12 EV
        conv(params.bands[0], 2.f, 3.f), // -10 EV
        conv(params.bands[0], 2.f, 3.f), //  -8 EV
        conv(params.bands[1], 2.f, 3.f), //  -6 EV
        conv(params.bands[2], 2.5f, 2.5f), //  -4 EV
        conv(params.bands[3], 3.f, 2.f), //  -2 EV
        conv(params.bands[4], 3.f, 2.f), //   0 EV
        conv(params.bands[4], 3.f, 2.f), //   2 EV
        conv(params.bands[4], 3.f, 2.f), //   4 EV
        conv(params.bands[4], 3.f, 2.f),  //   6 EV
        // this settings under are very rarely used...images with very high DR - I add a slider, but it's not the goal: the goal is white distribution (in very rare case)
        // I have not change "main" Tone Equalizer.
        conv(params.bands[5], 3.f, 2.f),  //   8 EV  Added for white distribution (Cam16 and Log encode) and images with very high DR
        conv(params.bands[5], 3.f, 2.f),  //   10 EV Added for white distribution(Cam16 and Log encode) and images with very high DR
        conv(params.bands[5], 3.f, 2.f)  //   12 EV Added for white distribution(Cam16 and Log encode) and images with very high DR
    };

    rtengine::TMatrix ws = rtengine::ICCStore::getInstance()->workingSpaceMatrix(workingProfile);

#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            Y[y][x] = rtengine::Color::rgbLuminance(R[y][x], G[y][x], B[y][x], ws);
        }
    }

    int detail = rtengine::LIM(params.regularization + 5, 0, 5);
    int radius = detail / scale + 0.5;
    float epsilon2 = 0.01f + 0.002f * rtengine::max(detail - 3, 0);

    if (radius > 0) {
        rtengine::guidedFilterLog(10.f, Y, radius, epsilon2, multithread);
    }

    if (params.regularization > 0) {
        array2D<float> Y2(W, H);
        constexpr float base_epsilon = 0.02f;
        constexpr float base_posterization = 5.f;

#ifdef _OPENMP
        #pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                float l = rtengine::LIM(log2(rtengine::max(Y[y][x], 1e-9f)), centers[0], centers[13]);
                float ll = round(l * base_posterization) / base_posterization;
                Y2[y][x] = Y[y][x];
                Y[y][x] = exp2(ll);
            }
        }

        radius = 350.0 / scale;
        epsilon2 = base_epsilon / float(6 - rtengine::min(params.regularization, 5));
        rtengine::guidedFilter(Y2, Y, Y, radius, epsilon2, multithread);
    }

    const auto gauss =
    [](float b, float x) -> float {
        return xexpf((-rtengine::SQR(x - b) / 4.0f));
    };

    // For every pixel luminance, the sum of the gaussian masks
    float w_sum = 0.f;

    for (int i = 0; i < 15; ++i) {
        w_sum += gauss(centers[i], 0.f);
    }

    constexpr float luma_lo = -14.f;
    constexpr float luma_hi = 6.f;

    const auto process_pixel =
    [&](float y) -> float {
        // convert to log space
        const float luma = rtengine::LIM(log2(rtengine::max(y, 0.f)), luma_lo, luma_hi);

        // build the correction as the sum of the contribution of each
        // luminance channel to current pixel
        float correction = 0.0f;

        for (int c = 0; c < 15; ++c)
        {
            correction += gauss(centers[c], luma) * factors[c];
        }

        correction /= w_sum;

        return correction;
    };

    std::vector<std::array<float, 3>> cur_colormap;
    if (params.show_colormap) {
        rtengine::lcmsMutex->lock();
        cmsHPROFILE in = rtengine::ICCStore::getInstance()->getsRGBProfile();
        cmsHPROFILE out = rtengine::ICCStore::getInstance()->workingSpace(workingProfile);
        cmsHTRANSFORM xform = cmsCreateTransform(in, TYPE_RGB_FLT, out, TYPE_RGB_FLT, INTENT_RELATIVE_COLORIMETRIC, cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE);
        rtengine::lcmsMutex->unlock();

        for (auto &c : colormap) {
            cur_colormap.push_back(c);
            auto &cc = cur_colormap.back();
            cmsDoTransform(xform, &cc[0], &cc[0], 1);
        }

        cmsDeleteTransform(xform);
    }

    const auto process_colormap =
        [&](float y) -> std::array<float, 3>
        {
            std::array<float, 3> ret = { 0.f, 0.f, 0.f };

            // convert to log space
            const float luma = rtengine::LIM(log2(rtengine::max(y, 0.f)), luma_lo, luma_hi);

            // build the correction as the sum of the contribution of each
            // luminance channel to current pixel
            for (int c = 0; c < 15; ++c) {
                float w = gauss(centers[c], luma);
                for (int i = 0; i < 3; ++i) {
                    ret[i] += w * cur_colormap[c][i];
                }
            }
            for (int i = 0; i < 3; ++i) {
                ret[i] = rtengine::LIM01(ret[i] / w_sum);
            }

            return ret;
        };


#ifdef __SSE2__
    vfloat vfactors[15];
    vfloat vcenters[15];

    for (int i = 0; i < 15; ++i) {
        vfactors[i] = F2V(factors[i]);
        vcenters[i] = F2V(centers[i]);
    }

    const auto vgauss =
    [](vfloat b, vfloat x) -> vfloat {
        static const vfloat fourv = F2V(4.f);
        return xexpf((-rtengine::SQR(x - b) / fourv));
    };

    vfloat zerov = F2V(0.f);
    vfloat vw_sum = F2V(w_sum);

    const vfloat vluma_lo = F2V(luma_lo);
    const vfloat vluma_hi = F2V(luma_hi);
    const vfloat xlog2v = F2V(xlogf(2.f));

    const auto vprocess_pixel =
    [&](vfloat y) -> vfloat {
        const vfloat luma = vminf(vmaxf(xlogf(vmaxf(y, zerov)) / xlog2v, vluma_lo), vluma_hi);

        vfloat correction = zerov;

        for (int c = 0; c < 15; ++c)
        {
            correction += vgauss(vcenters[c], luma) * vfactors[c];
        }

        correction /= vw_sum;

        return correction;
    };


    vfloat v1 = F2V(1.f);
    vfloat v65535 = F2V(65535.f);
#endif // __SSE2__


    if (params.show_colormap) {
        LUTf lut_r(65537), lut_g(65537), lut_b(65537);
        for (int i = 0; i < 65536; ++i) {
            float y = float(i)/65535.f;
            auto rgb = process_colormap(y);
            lut_r[i] = rgb[0];
            lut_g[i] = rgb[1];
            lut_b[i] = rgb[2];
        }
        lut_r[65536] = cur_colormap.back()[0];
        lut_g[65536] = cur_colormap.back()[1];
        lut_b[65536] = cur_colormap.back()[2];

#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                float cY = Y[y][x] * 65535.f;
                R[y][x] = lut_r[cY];
                G[y][x] = lut_g[cY];
                B[y][x] = lut_b[cY];
            }
        }
        return;
    }


    LUTf lut(65536);

    for (int i = 0; i < 65536; ++i) {
        float y = float(i) / 65535.f;
        float c = process_pixel(y);
        lut[i] = c;
    }

#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        int x = 0;


#ifdef __SSE2__

        for (; x < W - 3; x += 4) {
            vfloat cY = LVFU(Y[y][x]);
            vmask m = vmaskf_gt(cY, v1);
            vfloat corr;

            if (_mm_movemask_ps((vfloat)m)) {
                corr = vprocess_pixel(cY);
            } else {
                corr = lut[cY * v65535];
            }

            STVF(R[y][x], LVF(R[y][x]) * corr);
            STVF(G[y][x], LVF(G[y][x]) * corr);
            STVF(B[y][x], LVF(B[y][x]) * corr);
        }

#endif // __SSE2__

        for (; x < W; ++x) {
            float cY = Y[y][x];
            float corr = cY > 1.f ? process_pixel(cY) : lut[cY * 65535.f];
            R[y][x] *= corr;
            G[y][x] *= corr;
            B[y][x] *= corr;
        }
    }

}

}


namespace rtengine
{

void ImProcFunctions::toneEqualizer(
    Imagefloat *rgb,
    const ToneEqualizerParams &params,
    const Glib::ustring &workingProfile,
    double scale,
    bool multiThread)
{
    if (!params.enabled) {
        return;
    }

    BENCHFUN

    const float gain = 1.f / 65535.f * std::pow(2.f, -params.pivot);

    rgb->multiply(gain, multiThread);

    const int W = rgb->getWidth();
    const int H = rgb->getHeight();

    array2D<float> R(W, H, rgb->r.ptrs, ARRAY2D_BYREFERENCE);
    array2D<float> G(W, H, rgb->g.ptrs, ARRAY2D_BYREFERENCE);
    array2D<float> B(W, H, rgb->b.ptrs, ARRAY2D_BYREFERENCE);

    ::toneEqualizer(R, G, B, params, workingProfile, scale, multiThread);

    rgb->multiply(params.show_colormap ? 65535.f : 1.f/gain, multiThread);
}

void ImProcFunctions::toneEqualizer(Imagefloat *rgb)
{
    toneEqualizer(rgb, params->toneEqualizer, params->icm.workingProfile, scale, multiThread);
}

}
