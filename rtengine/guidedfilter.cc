/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2018 Alberto Griggio <alberto.griggio@gmail.com>
 *  Optimized 2019 Ingo Weyrich <heckflosse67@gmx.de>
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
 * This is a Fast Guided Filter implementation, derived directly from the
 * pseudo-code of the paper:
 *
 * Fast Guided Filter
 * by Kaiming He, Jian Sun
 *
 * available at https://arxiv.org/abs/1505.00996
*/

#include "guidedfilter.h"
#include "boxblur.h"
#include "rescale.h"
#include "imagefloat.h"
#define BENCHMARK
#include "StopWatch.h"
namespace rtengine {

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
    enum Op {MUL, DIVEPSILON, SUBMUL};

    const auto apply =
        [multithread, epsilon](Op op, array2D<float> &res, const array2D<float> &a, const array2D<float> &b, const array2D<float> &c=array2D<float>()) -> void
        {
            const int w = res.width();
            const int h = res.height();
            
#ifdef _OPENMP
            #pragma omp parallel for if (multithread)
#endif
            for (int y = 0; y < h; ++y) {
                for (int x = 0; x < w; ++x) {
                    switch (op) {
                        case MUL:
                            res[y][x] = a[y][x] * b[y][x];
                            break;
                        case DIVEPSILON:
                            res[y][x] = a[y][x] / (b[y][x] + epsilon); // note: the value of epsilon intentionally has an impact on the result. It is not only to avoid divisions by zero
                            break;
                        case SUBMUL:
                            res[y][x] = c[y][x] - (a[y][x] * b[y][x]);
                            break;
                        default:
                            assert(false);
                            res[y][x] = 0;
                            break;
                    }
                }
            }
        };

    const auto f_subsample =
        [multithread](array2D<float> &d, const array2D<float> &s) -> void
        {
            rescaleBilinear(s, d, multithread);
        };

    const auto f_mean =
        [multithread](array2D<float> &d, array2D<float> &s, int rad) -> void
        {
            rad = LIM(rad, 0, (min(s.width(), s.height()) - 1) / 2 - 1);
            boxblur(s, d, rad, s.width(), s.height(), multithread);
        };

    const int W = src.width();
    const int H = src.height();

    if (subsampling <= 0) {
        subsampling = calculate_subsampling(W, H, r);
    }

    const size_t w = W / subsampling;
    const size_t h = H / subsampling;
    const float r1 = float(r) / subsampling;

    array2D<float> I1(w, h);
    array2D<float> p1(w, h);

    f_subsample(I1, guide);

    if (&guide == &src) {
        f_mean(p1, I1, r1);

        apply(MUL, I1, I1, I1);        // I1 = I1 * I1

        f_mean(I1, I1, r1);

        apply(SUBMUL, I1, p1, p1, I1); // I1 = I1 - p1 * p1
        apply(DIVEPSILON, I1, I1, I1); // I1 = I1 / (I1 + epsilon)
        apply(SUBMUL, p1, I1, p1, p1); // p1 = p1 - I1 * p1

    } else {
        f_subsample(p1, src);

        array2D<float> meanI(w, h);
        f_mean(meanI, I1, r1);

        array2D<float> meanp(w, h);
        f_mean(meanp, p1, r1);

        apply(MUL, p1, I1, p1);

        f_mean(p1, p1, r1);

        apply(MUL, I1, I1, I1);

        f_mean(I1, I1, r1);

        apply(SUBMUL, I1, meanI, meanI, I1);
        apply(SUBMUL, p1, meanI, meanp, p1);
        apply(DIVEPSILON, I1, p1, I1);
        apply(SUBMUL, p1, I1, meanI, meanp);
    }

    f_mean(I1, I1, r1);
    f_mean(p1, p1, r1);

    const int Ws = I1.width();
    const int Hs = I1.height();
    const int Wd = dst.width();
    const int Hd = dst.height();

    const float col_scale = static_cast<float>(Ws) / static_cast<float>(Wd);
    const float row_scale = static_cast<float>(Hs) / static_cast<float>(Hd);

#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif

    for (int y = 0; y < Hd; ++y) {
        const float ymrs = y * row_scale;
        for (int x = 0; x < Wd; ++x) {
            dst[y][x] = getBilinearValue(I1, x * col_scale, ymrs) * guide[y][x] + getBilinearValue(p1, x * col_scale, ymrs);
        }
    }
}

} // namespace rtengine
