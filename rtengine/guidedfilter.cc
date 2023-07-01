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

/**
 * This is a Fast Guided Filter implementation, derived directly from the
 * pseudo-code of the paper:
 *
 * Fast Guided Filter
 * by Kaiming He, Jian Sun
 *
 * available at https://arxiv.org/abs/1505.00996
 */

#include "array2D.h"
#include "boxblur.h"
#include "guidedfilter.h"
#include "boxblur.h"
#include "sleef.h"
#include "rescale.h"
#include "imagefloat.h"

namespace rtengine {

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


namespace {

int calculate_subsampling(int w, int h, int r)
{
    if (r == 1) {
        return 1;
    }
    
    if (max(w, h) <= 600) {
        return 1;
    }
    
    for (int s = 5; s > 0; --s) {
        if (r % s == 0) {
            return s;
        }
    }

    return LIM(r / 2, 2, 4);
}

} // namespace


void guidedFilter(const array2D<float> &guide, const array2D<float> &src, array2D<float> &dst, int r, float epsilon, bool multithread, int subsampling)
{

    const int W = src.getWidth();
    const int H = src.getHeight();

    if (subsampling <= 0) {
        subsampling = calculate_subsampling(W, H, r);
    }

    enum Op { MUL, DIVEPSILON, ADD, SUB, ADDMUL, SUBMUL };

    const auto apply =
        [=](Op op, array2D<float> &res, const array2D<float> &a, const array2D<float> &b, const array2D<float> &c=array2D<float>()) -> void
        {
            const int w = res.getWidth();
            const int h = res.getHeight();
            
#ifdef _OPENMP
            #pragma omp parallel for if (multithread)
#endif
            for (int y = 0; y < h; ++y) {
                for (int x = 0; x < w; ++x) {
                    float r;
                    float aa = a[y][x];
                    float bb = b[y][x];
                    switch (op) {
                    case MUL:
                        r = aa * bb;
                        break;
                    case DIVEPSILON:
                        r = aa / (bb + epsilon);
                        break;
                    case ADD:
                        r = aa + bb;
                        break;
                    case SUB:
                        r = aa - bb;
                        break;
                    case ADDMUL:
                        r = aa * bb + c[y][x];
                        break;
                    case SUBMUL:
                        r = c[y][x] - (aa * bb);
                        break;
                    default:
                        assert(false);
                        r = 0;
                        break;
                    }
                    res[y][x] = r;
                }
            }
        };

    // use the terminology of the paper (Algorithm 2)
    const array2D<float> &I = guide;
    const array2D<float> &p = src;
    array2D<float> &q = dst;

    const auto f_subsample =
        [=](array2D<float> &d, const array2D<float> &s) -> void
        {
            if (d.getWidth() == s.getWidth() && d.getHeight() == s.getHeight()) {
#ifdef _OPENMP
#               pragma omp parallel for if (multithread)
#endif
                for (int y = 0; y < s.getHeight(); ++y) {
                    for (int x = 0; x < s.getWidth(); ++x) {
                        d[y][x] = s[y][x];
                    }
                }
            } else {
                rescaleBilinear(s, d, multithread);
            }
        };

    // const auto f_upsample = f_subsample;
    
    const size_t w = W / subsampling;
    const size_t h = H / subsampling;

    const auto f_mean =
        [multithread](array2D<float> &d, array2D<float> &s, int rad) -> void
        {
            rad = LIM(rad, 0, (min(s.getWidth(), s.getHeight()) - 1) / 2 - 1);
           // boxblur(s, d, rad, s.getWidth(), s.getHeight(), multithread);
            boxblur(static_cast<float**>(s), static_cast<float**>(d), rad, s.getWidth(), s.getHeight(), multithread);
        };

    array2D<float> I1(w, h);
    array2D<float> p1(w, h);

    f_subsample(I1, I);
    f_subsample(p1, p);

    DEBUG_DUMP(I);
    DEBUG_DUMP(p);
    DEBUG_DUMP(I1);
    DEBUG_DUMP(p1);

    float r1 = float(r) / subsampling;

    array2D<float> meanI(w, h);
    f_mean(meanI, I1, r1);
    DEBUG_DUMP(meanI);

    array2D<float> meanp(w, h);
    f_mean(meanp, p1, r1);
    DEBUG_DUMP(meanp);

    array2D<float> &corrIp = p1;
    apply(MUL, corrIp, I1, p1);
    f_mean(corrIp, corrIp, r1);
    DEBUG_DUMP(corrIp);

    array2D<float> &corrI = I1;
    apply(MUL, corrI, I1, I1);
    f_mean(corrI, corrI, r1);
    DEBUG_DUMP(corrI);

    array2D<float> &varI = corrI;
    apply(SUBMUL, varI, meanI, meanI, corrI);
    DEBUG_DUMP(varI);

    array2D<float> &covIp = corrIp;
    apply(SUBMUL, covIp, meanI, meanp, corrIp);
    DEBUG_DUMP(covIp);

    array2D<float> &a = varI;
    apply(DIVEPSILON, a, covIp, varI);
    DEBUG_DUMP(a);

    array2D<float> &b = covIp;
    apply(SUBMUL, b, a, meanI, meanp);
    DEBUG_DUMP(b);

    array2D<float> &meana = a;
    f_mean(meana, a, r1);
    DEBUG_DUMP(meana);

    array2D<float> &meanb = b;
    f_mean(meanb, b, r1);
    DEBUG_DUMP(meanb);

    // speedup by heckflosse67
    const int Ws = meana.getWidth();
    const int Hs = meana.getHeight();
    const int Wd = q.getWidth();
    const int Hd = q.getHeight();
    const float col_scale = float(Ws) / float(Wd);
    const float row_scale = float(Hs) / float(Hd);

#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < Hd; ++y) {
        float ymrs = y * row_scale; 
        for (int x = 0; x < Wd; ++x) {
            q[y][x] = getBilinearValue(meana, x * col_scale, ymrs) * I[y][x] + getBilinearValue(meanb, x * col_scale, ymrs);
        }
    }
}


void guidedFilterLog(const array2D<float> &guide, float base, array2D<float> &chan, int r, float eps, bool multithread, int subsampling)
{
#ifdef _OPENMP
#    pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < chan.getHeight(); ++y) {
        for (int x = 0; x < chan.getWidth(); ++x) {
            chan[y][x] = xlin2log(max(chan[y][x], 0.f), base);
        }
    }

    guidedFilter(guide, chan, chan, r, eps, multithread, subsampling);

#ifdef _OPENMP
#    pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < chan.getHeight(); ++y) {
        for (int x = 0; x < chan.getWidth(); ++x) {
            chan[y][x] = xlog2lin(max(chan[y][x], 0.f), base);
        }
    }
}


void guidedFilterLog(float base, array2D<float> &chan, int r, float eps, bool multithread, int subsampling)
{
    guidedFilterLog(chan, base, chan, r, eps, multithread, subsampling);
}

} // namespace rtengine




