/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2017 Gabor Horvath <hgabor@rawtherapee.com>
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

#include <memory>

#include "labimage.h"

namespace rtengine
{

LabImage::LabImage (int w, int h, bool initZero, bool multiThread) : W(w), H(h)
{
    allocLab(w, h);
    if (initZero) {
        clear(multiThread);
    }
}

LabImage::~LabImage ()
{
    deleteLab();
}

void LabImage::CopyFrom(LabImage *Img, bool multiThread)
{
#ifdef _OPENMP
    #pragma omp parallel sections if(multiThread)
    {
        #pragma omp section
        memcpy(L[0], Img->L[0], static_cast<std::size_t>(W) * H * sizeof(float));
        #pragma omp section
        memcpy(a[0], Img->a[0], static_cast<std::size_t>(W) * H * sizeof(float));
        #pragma omp section
        memcpy(b[0], Img->b[0], static_cast<std::size_t>(W) * H * sizeof(float));
    }
#else
    memcpy(data, Img->data, static_cast<std::size_t>(W) * H * 3 * sizeof(float));
#endif
}

void LabImage::getPipetteData (float &v1, float &v2, float &v3, int posX, int posY, int squareSize)
{
    float accumulator_L = 0.f;
    float accumulator_a = 0.f;
    float accumulator_b = 0.f;
    unsigned long int n = 0;
    int halfSquare = squareSize / 2;

    for (int iy = posY - halfSquare; iy < posY - halfSquare + squareSize; ++iy) {
        for (int ix = posX - halfSquare; ix < posX - halfSquare + squareSize; ++ix) {
            if (ix >= 0 && iy >= 0 && ix < W && iy < H) {
                accumulator_L += L[iy][ix];
                accumulator_a += a[iy][ix];
                accumulator_b += b[iy][ix];
                ++n;
            }
        }
    }

    v1 = n ? accumulator_L / float(n) : 0.f;
    v2 = n ? accumulator_a / float(n) : 0.f;
    v3 = n ? accumulator_b / float(n) : 0.f;
}

void LabImage::allocLab(size_t w, size_t h)
{
    L = new float*[h];
    a = new float*[h];
    b = new float*[h];

    data = new float [w * h * 3];
    float * index = data;

    for (size_t i = 0; i < h; i++) {
        L[i] = index + i * w;
    }

    index += w * h;

    for (size_t i = 0; i < h; i++) {
        a[i] = index + i * w;
    }

    index += w * h;

    for (size_t i = 0; i < h; i++) {
        b[i] = index + i * w;
    }
}

void LabImage::deleteLab()
{
    delete [] L;
    delete [] a;
    delete [] b;
    delete [] data;
}

void LabImage::reallocLab()
{
    allocLab(W, H);
}

void LabImage::clear(bool multiThread) {
#ifdef _OPENMP
        #pragma omp parallel for if(multiThread)
#endif
        for(size_t i = 0; i < static_cast<size_t>(H) * W * 3; ++i) {
            data[i] = 0.f;
        }
    }
}
