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
#include "shmap.h"
#include "gauss.h"
#include "rtengine.h"
#include "rt_math.h"
#include "rawimagesource.h"
#include "jaggedarray.h"
#undef THREAD_PRIORITY_NORMAL
#include "opthelper.h"

namespace rtengine
{

extern const Settings* settings;

SHMap::SHMap (int w, int h, bool multiThread) : W(w), H(h), multiThread(multiThread)
{

    map = new float*[H];

    for (int i = 0; i < H; i++) {
        map[i] = new float[W];
    }

}

SHMap::~SHMap ()
{

    for (int i = 0; i < H; i++) {
        delete [] map[i];
    }

    delete [] map;
}

void SHMap::fillLuminance( Imagefloat * img, float **luminance, double lumi[3] ) // fill with luminance
{

#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++) {
            luminance[i][j] = lumi[0] * std::max(img->r(i, j), 0.f) + lumi[1] * std::max(img->g(i, j), 0.f) + lumi[2] * std::max(img->b(i, j), 0.f);
        }

}

void SHMap::fillLuminanceL( float ** L, float **luminance) // fill with luminance
{

#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++) {
            luminance[i][j] = std::max(L[i][j], 0.f) ;//we can put here some enhancements Gamma, compression data,...
        }

}

void SHMap::update (Imagefloat* img, double radius, double lumi[3], bool hq, int skip)
{

    if (!hq) {
        fillLuminance( img, map, lumi);

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            gaussianBlur (map, map, W, H, radius);
        }
    }

    else {
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //experimental dirpyr shmap

        float thresh = (100.f * radius); //1000;

        // set up range function
        // calculate size of Lookup table. That's possible because from a value k for all i>=k rangefn[i] will be exp(-10)
        // So we use this fact and the automatic clip of lut to reduce the size of lut and the number of calculations to fill the lut
        // In past this lut had only integer precision with rangefn[i] = 0 for all i>=k
        // We set the last element to a small epsilon 1e-15 instead of zero to avoid divisions by zero
        const int lutSize = thresh * sqrtf(10.f) + 1;
        thresh *= thresh;
        LUTf rangefn(lutSize);

        for (int i = 0; i < lutSize - 1; i++) {
            rangefn[i] = xexpf(-min(10.f, (static_cast<float>(i) * i) / thresh )); //*intfactor;
        }

        rangefn[lutSize - 1] = 1e-15f;

        // We need one temporary buffer
        const JaggedArray<float> buffer (W, H);

        // the final result has to be in map
        // for an even number of levels that means: map => buffer, buffer => map
        // for an odd number of levels that means: buffer => map, map => buffer, buffer => map
        // so let's calculate the number of levels first
        // There are at least two levels
        int numLevels = 2;
        int scale = 2;

        while (skip * scale < 16) {
            scale *= 2;
            numLevels++;
        }

        float ** dirpyrlo[2];

        if(numLevels & 1) { // odd number of levels, start with buffer
            dirpyrlo[0] = buffer;
            dirpyrlo[1] = map;
        } else { // even number of levels, start with map
            dirpyrlo[0] = map;
            dirpyrlo[1] = buffer;
        }

        fillLuminance( img, dirpyrlo[0], lumi);

        scale = 1;
        int level = 0;
        int indx = 0;
        dirpyr_shmap(dirpyrlo[indx], dirpyrlo[1 - indx], W, H, rangefn, level, scale );
        scale *= 2;
        level ++;
        indx = 1 - indx;

        while (skip * scale < 16) {
            dirpyr_shmap(dirpyrlo[indx], dirpyrlo[1 - indx], W, H, rangefn, level, scale );
            scale *= 2;
            level ++;
            indx = 1 - indx;
        }

        dirpyr_shmap(dirpyrlo[indx], dirpyrlo[1 - indx], W, H, rangefn, level, scale );
    }

    // update average, minimum, maximum
    double _avg = 0.0f; // use double precision to gain precision especially at systems with few cores and big pictures (error for 36 MPixel on single core was about 8% with float)
    min_f = 65535;
    max_f = 0;
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        float _min_f = 65535.0f;
        float _max_f = 0.0f;
        float _val;
#ifdef _OPENMP
        #pragma omp for reduction(+:_avg) schedule(dynamic,16) nowait
#endif

        for (int i = 0; i < H; i++)
            for (int j = 0; j < W; j++) {
                _val = map[i][j];

                if (_val < _min_f) {
                    _min_f = _val;
                }

                if (_val > _max_f) {
                    _max_f = _val;
                }

                _avg += _val;
            }

#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            if(_min_f < min_f ) {
                min_f = _min_f;
            }

            if(_max_f > max_f ) {
                max_f = _max_f;
            }
        }
    }
    _avg /= ((H) * (W));
    avg = _avg;

}

void SHMap::updateL (float** L, double radius, bool hq, int skip)
{

    if (!hq) {
        fillLuminanceL( L, map);
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            gaussianBlur (map, map, W, H, radius);
        }
    }

    else

    {
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //experimental dirpyr shmap
        float thresh = (100.f * radius); //1000;
        int levrad = 16;
        levrad = 2; //for retinex - otherwise levrad = 16
        // set up range function
        // calculate size of Lookup table. That's possible because from a value k for all i>=k rangefn[i] will be exp(-10)
        // So we use this fact and the automatic clip of lut to reduce the size of lut and the number of calculations to fill the lut
        // In past this lut had only integer precision with rangefn[i] = 0 for all i>=k
        // We set the last element to a small epsilon 1e-15 instead of zero to avoid divisions by zero
        const int lutSize = (int) thresh * sqrtf(10.f) + 1;
        thresh *= thresh;
        LUTf rangefn(lutSize);

        for (int i = 0; i < lutSize - 1; i++) {
            rangefn[i] = xexpf(-min(10.f, (static_cast<float>(i) * i) / thresh )); //*intfactor;
        }

        rangefn[lutSize - 1] = 1e-15f;
        //printf("lut=%d rf5=%f rfm=%f\n thre=%f",lutSize, rangefn[5],rangefn[lutSize-10],thresh );

        // We need one temporary buffer
        const JaggedArray<float> buffer (W, H);

        // the final result has to be in map
        // for an even number of levels that means: map => buffer, buffer => map
        // for an odd number of levels that means: buffer => map, map => buffer, buffer => map
        // so let's calculate the number of levels first
        // There are at least two levels
        int numLevels = 2;
        int scale = 2;

        while (skip * scale < levrad) {
            scale *= 2;
            numLevels++;
        }

        //printf("numlev=%d\n",numLevels);
        float ** dirpyrlo[2];

        if(numLevels & 1) { // odd number of levels, start with buffer
            dirpyrlo[0] = buffer;
            dirpyrlo[1] = map;
        } else { // even number of levels, start with map
            dirpyrlo[0] = map;
            dirpyrlo[1] = buffer;
        }

        fillLuminanceL( L, dirpyrlo[0]);

        scale = 1;
        int level = 0;
        int indx = 0;
        dirpyr_shmap(dirpyrlo[indx], dirpyrlo[1 - indx], W, H, rangefn, level, scale );
        scale *= 2;
        level ++;
        indx = 1 - indx;

        while (skip * scale < levrad) {
            dirpyr_shmap(dirpyrlo[indx], dirpyrlo[1 - indx], W, H, rangefn, level, scale );
            scale *= 2;
            level ++;
            indx = 1 - indx;
        }

        dirpyr_shmap(dirpyrlo[indx], dirpyrlo[1 - indx], W, H, rangefn, level, scale );
    }

    // update average, minimum, maximum
    double _avg = 0.0f; // use double precision to gain precision especially at systems with few cores and big pictures (error for 36 MPixel on single core was about 8% with float)
    min_f = 65535;
    max_f = 0;
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        float _min_f = 65535.0f;
        float _max_f = 0.0f;
        float _val;
#ifdef _OPENMP
        #pragma omp for reduction(+:_avg) schedule(dynamic,16) nowait
#endif

        for (int i = 0; i < H; i++)
            for (int j = 0; j < W; j++) {
                _val = map[i][j];

                if (_val < _min_f) {
                    _min_f = _val;
                }

                if (_val > _max_f) {
                    _max_f = _val;
                }

                _avg += _val;
            }

#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            if(_min_f < min_f ) {
                min_f = _min_f;
            }

            if(_max_f > max_f ) {
                max_f = _max_f;
            }
        }
    }
    _avg /= ((H) * (W));
    avg = _avg;

}


void SHMap::forceStat (float max_, float min_, float avg_)
{

    max_f = max_;
    min_f = min_;
    avg = avg_;
}

SSEFUNCTION void SHMap::dirpyr_shmap(float ** data_fine, float ** data_coarse, int width, int height, LUTf & rangefn, int level, int scale)
{
    //scale is spacing of directional averaging weights

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // calculate weights, compute directionally weighted average

    int scalewin, halfwin;

    if(level < 2) {
        halfwin = 1;
        scalewin = halfwin * scale;

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
#if defined( __SSE2__ ) && defined( __x86_64__ )
            __m128 dirwtv, valv, normv, dftemp1v, dftemp2v, fg;
#endif // __SSE2__
            int j;
#ifdef _OPENMP
            #pragma omp for
#endif

            for(int i = 0; i < height; i++) {
                float dirwt;

                for(j = 0; j < scalewin; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for(int inbr = max(i - scalewin, i % scale); inbr <= min(i + scalewin, height - 1); inbr += scale) {
                        for (int jnbr = j % scale; jnbr <= j + scalewin; jnbr += scale) {
                            //printf("dat=%f ",abs(data_fine[inbr][jnbr] - data_fine[i][j]));
                            dirwt = ( rangefn[abs(data_fine[inbr][jnbr] - data_fine[i][j])] );
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }

                    data_coarse[i][j] = val / norm; // low pass filter
                }

#if defined( __SSE2__ ) && defined( __x86_64__ )
                int inbrMin = max(i - scalewin, i % scale);

                for(; j < (width - scalewin) - 3; j += 4) {
                    valv = _mm_setzero_ps();
                    normv = _mm_setzero_ps();
                    dftemp1v = LVFU(data_fine[i][j]);

                    for(int inbr = inbrMin; inbr <= min(i + scalewin, height - 1); inbr += scale) {
                        for (int jnbr = j - scalewin; jnbr <= j + scalewin; jnbr += scale) {
                            dftemp2v = LVFU(data_fine[inbr][jnbr]);
                            dirwtv = ( rangefn[_mm_cvttps_epi32(vabsf(dftemp2v - dftemp1v))] );
                            valv += dirwtv * dftemp2v;
                            normv += dirwtv;
                        }
                    }

                    _mm_storeu_ps( &data_coarse[i][j], valv / normv);
                }

                for(; j < width - scalewin; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for(int inbr = inbrMin; inbr <= min(i + scalewin, height - 1); inbr += scale) {
                        for (int jnbr = j - scalewin; jnbr <= j + scalewin; jnbr += scale) {
                            dirwt = ( rangefn[abs(data_fine[inbr][jnbr] - data_fine[i][j])] );
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }

                    data_coarse[i][j] = val / norm; // low pass filter
                }

#else

                for(; j < width - scalewin; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for(int inbr = max(i - scalewin, i % scale); inbr <= min(i + scalewin, height - 1); inbr += scale) {
                        for (int jnbr = j - scalewin; jnbr <= j + scalewin; jnbr += scale) {
                            dirwt = ( rangefn[abs(data_fine[inbr][jnbr] - data_fine[i][j])] );
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }

                    data_coarse[i][j] = val / norm; // low pass filter
                }

#endif

                for(; j < width; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for(int inbr = max(i - scalewin, i % scale); inbr <= min(i + scalewin, height - 1); inbr += scale) {
                        for (int jnbr = j - scalewin; jnbr < width; jnbr += scale) {
                            dirwt = ( rangefn[abs(data_fine[inbr][jnbr] - data_fine[i][j])] );
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }

                    data_coarse[i][j] = val / norm; // low pass filter
                }
            }
        }
    } else {
        halfwin = 2;
        scalewin = halfwin * scale;
        int domker[5][5] = {{1, 1, 1, 1, 1}, {1, 2, 2, 2, 1}, {1, 2, 2, 2, 1}, {1, 2, 2, 2, 1}, {1, 1, 1, 1, 1}};
        //generate domain kernel

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
#if defined( __SSE2__ ) && defined( __x86_64__ )
            __m128 dirwtv, valv, normv, dftemp1v, dftemp2v, fgg;
            float domkerv[5][5][4] __attribute__ ((aligned (16))) = {{{1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}}, {{1, 1, 1, 1}, {2, 2, 2, 2}, {2, 2, 2, 2}, {2, 2, 2, 2}, {1, 1, 1, 1}}, {{1, 1, 1, 1}, {2, 2, 2, 2}, {2, 2, 2, 2}, {2, 2, 2, 2}, {1, 1, 1, 1}}, {{1, 1, 1, 1}, {2, 2, 2, 2}, {2, 2, 2, 2}, {2, 2, 2, 2}, {1, 1, 1, 1}}, {{1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}}};

#endif // __SSE2__
            int j;
#ifdef _OPENMP
            #pragma omp for schedule(dynamic,16)
#endif

            for(int i = 0; i < height; i++) {
                float dirwt;

                for(j = 0; j < scalewin; j++) {
                    float val = 0.f;
                    float norm = 0.f;

                    for(int inbr = max(i - scalewin, i % scale); inbr <= min(i + scalewin, height - 1); inbr += scale) {
                        for (int jnbr = j % scale; jnbr <= j + scalewin; jnbr += scale) {
                            dirwt = ( domker[(inbr - i) / scale + halfwin][(jnbr - j) / scale + halfwin] * rangefn[abs(data_fine[inbr][jnbr] - data_fine[i][j])] );
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }

                    data_coarse[i][j] = val / norm; // low pass filter
                }

#if defined( __SSE2__ ) && defined( __x86_64__ )

                for(; j < width - scalewin - 3; j += 4) {
                    valv = _mm_setzero_ps();
                    normv = _mm_setzero_ps();
                    dftemp1v = LVFU(data_fine[i][j]);

                    for(int inbr = max(i - scalewin, i % scale); inbr <= MIN(i + scalewin, height - 1); inbr += scale) {
                        int indexihlp = (inbr - i) / scale + halfwin;

                        for (int jnbr = j - scalewin, indexjhlp = 0; jnbr <= j + scalewin; jnbr += scale, indexjhlp++) {
                            dftemp2v = LVFU(data_fine[inbr][jnbr]);
                            dirwtv = ( _mm_load_ps((float*)&domkerv[indexihlp][indexjhlp]) * rangefn[_mm_cvttps_epi32(vabsf(dftemp2v - dftemp1v))] );
                            valv += dirwtv * dftemp2v;
                            normv += dirwtv;
                        }
                    }

                    _mm_storeu_ps( &data_coarse[i][j], valv / normv);
                }

                for(; j < width - scalewin; j++) {
                    float val = 0;
                    float norm = 0;

                    for(int inbr = max(i - scalewin, i % scale); inbr <= min(i + scalewin, height - 1); inbr += scale) {
                        for (int jnbr = j - scalewin; jnbr <= j + scalewin; jnbr += scale) {
                            dirwt = ( domker[(inbr - i) / scale + halfwin][(jnbr - j) / scale + halfwin] * rangefn[abs(data_fine[inbr][jnbr] - data_fine[i][j])] );
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }

                    data_coarse[i][j] = val / norm; // low pass filter
                }

#else

                for(; j < width - scalewin; j++) {
                    float val = 0;
                    float norm = 0;

                    for(int inbr = max(i - scalewin, i % scale); inbr <= min(i + scalewin, height - 1); inbr += scale) {
                        for (int jnbr = j - scalewin; jnbr <= j + scalewin; jnbr += scale) {
                            dirwt = ( domker[(inbr - i) / scale + halfwin][(jnbr - j) / scale + halfwin] * rangefn[abs(data_fine[inbr][jnbr] - data_fine[i][j])] );
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }

                    data_coarse[i][j] = val / norm; // low pass filter
                }

#endif

                for(; j < width; j++) {
                    float val = 0;
                    float norm = 0;

                    for(int inbr = max(i - scalewin, i % scale); inbr <= min(i + scalewin, height - 1); inbr += scale) {
                        for (int jnbr = j - scalewin; jnbr < width; jnbr += scale) {
                            dirwt = ( domker[(inbr - i) / scale + halfwin][(jnbr - j) / scale + halfwin] * rangefn[abs(data_fine[inbr][jnbr] - data_fine[i][j])] );
                            val += dirwt * data_fine[inbr][jnbr];
                            norm += dirwt;
                        }
                    }

                    data_coarse[i][j] = val / norm; // low pass filter
                }
            }
        }

    }

}

}//end of SHMap
