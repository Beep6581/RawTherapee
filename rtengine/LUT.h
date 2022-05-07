/*
 * LUT.h
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2011 Jan Rinze Peterzon (janrinze@gmail.com)
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
 *  Declaration of flexible Lookup Tables
 *
 *  Usage:
 *
 *      LUT<type> name (size);
 *      LUT<type> name (size, flags);
 *
 *      creates an array which is valid within the normal C/C++ scope "{ ... }"
 *
 *      access to elements is a simple as:
 *
 *          LUT<float> my_lut (10);
 *          float value = my_lut[3];
 *          float value = my_lut[2.5]; // this will interpolate
 *
 *      when using a float type index it will interpolate the lookup values
 *
 *      extra setting in flags: (clipping is set by default)
 *      LUT_CLIP_ABOVE
 *      LUT_CLIP_BELOW
 *
 *      example:
 *          LUT<float> my_lut (10,LUT_CLIP_BELOW);
 *          float value = my_lut[22.5];  // this will extrapolate
 *          float value = my_lut[-22.5]; // this will not extrapolate
 *
 *          LUT<float> my_lut (10,0); // this will extrapolate on either side
 *
 *      shotcuts:
 *
 *          LUTf stands for LUT<float>
 *          LUTi stands for LUT<int>
 *          LUTu stands for LUT<unsigned int>
 *          LUTd stands for LUT<double>
 *          LUTuc stands for LUT<unsigned char>
 */

#pragma once

#include <algorithm>
#include <cstring>
#include <cstdint>
#include <cassert>
#include <vector>

#ifndef NDEBUG
#include <fstream>
#endif

#include "opthelper.h"
#include "rt_math.h"

// Bit representations of flags
enum {
    LUT_CLIP_OFF,   // LUT does not clip input values
    LUT_CLIP_BELOW, // LUT clips input values at lower bound
    LUT_CLIP_ABOVE  // LUT clips input values at upper bound
};

template<typename T>
class LUT;

using LUTf = LUT<float>;
using LUTi = LUT<int32_t>;
using LUTu = LUT<uint32_t>;
using LUTd = LUT<double>;
using LUTuc = LUT<uint8_t>;

template<typename T>
class LUT
{
protected:
    // list of variables ordered to improve cache speed
    int maxs;
    float maxsf;
    T * data;
    unsigned int clip;
    unsigned int size;
    unsigned int upperBound;  // always equals size-1, parameter created for performance reason
private:
    unsigned int owner;
#ifdef __SSE2__
    alignas(16) vfloat maxsv;
    alignas(16) vfloat sizev;
    alignas(16) vint sizeiv;
#endif
public:
    /// convenience flag! If one doesn't want to delete the buffer but want to flag it to be recomputed...
    /// The user have to handle it itself, even if some method can (re)initialize it
    bool dirty;

    explicit LUT(int s, int flags = LUT_CLIP_BELOW | LUT_CLIP_ABOVE, bool initZero = false)
    {
#ifndef NDEBUG

        if (s <= 0) {
            printf("s<=0!\n");
        }

        assert (s > 0);
#endif
        dirty = true;
        clip = flags;
        // Add a few extra elements so [](vfloat) won't access out-of-bounds memory.
        // The routine would still produce the right answer, but might cause issues
        // with address/heap checking programs.
        data = new T[s + 3];
        owner = 1;
        size = s;
        upperBound = size - 1;
        maxs = size - 2;
        maxsf = (float)maxs;
#ifdef __SSE2__
        maxsv =  F2V( maxs );
        sizeiv =  _mm_set1_epi32( (int)(size - 1) );
        sizev = F2V( size - 1 );
#endif
        if (initZero) {
            clear();
        }
    }

    explicit LUT(const std::vector<T>& input, int flags = LUT_CLIP_BELOW | LUT_CLIP_ABOVE) :
        maxs(input.size() - 2),
        maxsf(maxs),
        data(new T[input.size() + 3]), // Add a few extra elements so [](vfloat) won't access out-of-bounds memory.
        clip(flags),
        size(input.size()),
        upperBound(size - 1),
        owner(1),
#ifdef __SSE2__
        maxsv(F2V(maxs)),
        sizev(F2V(size - 1)),
        sizeiv(_mm_set1_epi32(size - 1)),
#endif
        dirty(true)
    {
#ifndef NDEBUG

        if (input.empty()) {
            printf("s=0!\n");
        }

        assert(!input.empty());
#endif
        std::copy_n(input.begin(), input.size(), data);
    }

    void operator ()(int s, int flags = LUT_CLIP_BELOW | LUT_CLIP_ABOVE, bool initZero = false)
    {
#ifndef NDEBUG

        if (s <= 0) {
            printf("s<=0!\n");
        }

        assert (s > 0);
#endif

        if (owner && data) {
            delete[] data;
        }

        dirty = true; // Assumption!
        clip = flags;
        // See comment in constructor.
        data = new T[s + 3];
        owner = 1;
        size = s;
        upperBound = size - 1;
        maxs = size - 2;
        maxsf = (float)maxs;
#ifdef __SSE2__
        maxsv =  F2V( maxs );
        sizeiv =  _mm_set1_epi32( (int)(size - 1) );
        sizev = F2V( size - 1 );
#endif
        if (initZero) {
            clear();
        }

    }

    LUT()
    {
        data = nullptr;
        reset();
#ifdef __SSE2__
        maxsv = ZEROV;
        sizev = ZEROV;
        sizeiv = _mm_setzero_si128();
#endif
    }

    ~LUT()
    {
        if (owner) {
            delete[] data;
#ifndef NDEBUG
            data = (T*)0xBAADF00D;
#endif
        }
    }

    explicit LUT(const LUT&) = delete;

    void setClip(int flags)
    {
        clip = flags;
    }

    int getClip() const {
        return clip;
    }

    /** @brief Get the number of element in the LUT (i.e. dimension of the array)
     *  For a LUT(500), it will return 500
     *  @return number of element in the array
     */
    unsigned int getSize() const
    {
        return size;
    }

    /** @brief Get the highest value possible (i.e. dimension of the array)
     *  For a LUT(500), it will return 499, because 500 elements, starting from 0, goes up to 499
     *  @return number of element in the array
     */
    unsigned int getUpperBound() const
    {
        return size > 0 ? upperBound : 0;
    }

    LUT<T>& operator=(const LUT<T>& rhs)
    {
        if (this != &rhs) {
            if (rhs.size > this->size) {
                delete [] this->data;
                this->data = nullptr;
            }

            if (this->data == nullptr) {
                // See comment in constructor.
                this->data = new T[rhs.size + 3];
            }

            this->clip = rhs.clip;
            this->owner = 1;
            memcpy(this->data, rhs.data, rhs.size * sizeof(T));
            this->size = rhs.size;
            this->upperBound = rhs.upperBound;
            this->maxs = this->size - 2;
            this->maxsf = (float)this->maxs;
#ifdef __SSE2__
            this->maxsv =  F2V( this->size - 2);
            this->sizeiv =  _mm_set1_epi32( (int)(this->size - 1) );
            this->sizev = F2V( this->size - 1 );
#endif
        }

        return *this;
    }

    // handy to sum up per thread histograms. #pragma omp simd speeds up the loop by about factor 3 for LUTu (uint32_t).
    LUT<T>& operator+=(const LUT<T>& rhs)
    {
        if (rhs.size == this->size) {
#ifdef _OPENMP
            #pragma omp simd
#endif

            for(unsigned int i = 0; i < this->size; i++) {
                data[i] += rhs.data[i];
            }
        }

        return *this;
    }

    // multiply all elements of LUT<float> with a constant float value
    template<typename U = T, typename = typename std::enable_if<std::is_same<U, float>::value>::type>
    LUT<float>& operator*=(float factor)
    {
#ifdef _OPENMP
        #pragma omp simd
#endif

        for(unsigned int i = 0; i < this->size; i++) {
            data[i] *= factor;
        }

        return *this;
    }

    // divide all elements of LUT<float> by a constant float value
    template<typename U = T, typename = typename std::enable_if<std::is_same<U, float>::value>::type>
    LUT<float>& operator/=(float divisor)
    {
#ifdef _OPENMP
        #pragma omp simd
#endif

        for(unsigned int i = 0; i < this->size; i++) {
            data[i] /= divisor;
        }

        return *this;
    }


    // use with integer indices
    T& operator[](int index) const
    {
        return data[ rtengine::LIM<int>(index, 0, upperBound) ];
    }

#ifdef __SSE2__


    // NOTE: This function requires LUTs which clips only at lower bound
    vfloat cb(vfloat indexv) const
    {
        static_assert(std::is_same<T, float>::value, "This method only works for float LUTs");

        // Clamp and convert to integer values. Extract out of SSE register because all
        // lookup operations use regular addresses.
        vfloat clampedIndexes = vclampf(indexv, ZEROV, maxsv); // this automagically uses ZEROV in case indexv is NaN
        vint indexes = _mm_cvttps_epi32(clampedIndexes);
        int indexArray[4];
        _mm_storeu_si128(reinterpret_cast<__m128i*>(&indexArray[0]), indexes);

        // Load data from the table. This reads more than necessary, but there don't seem
        // to exist more granular operations (though we could try non-SSE).
        // Cast to int for convenience in the next operation (partial transpose).
        vint values[4];
        for (int i = 0; i < 4; ++i) {
            values[i] = _mm_castps_si128(LVFU(data[indexArray[i]]));
        }

        // Partial 4x4 transpose operation. We want two new vectors, the first consisting
        // of [values[0][0] ... values[3][0]] and the second [values[0][1] ... values[3][1]].
        __m128i temp0 = _mm_unpacklo_epi32(values[0], values[1]);
        __m128i temp1 = _mm_unpacklo_epi32(values[2], values[3]);
        vfloat lowerVal = _mm_castsi128_ps(_mm_unpacklo_epi64(temp0, temp1));
        vfloat upperVal = _mm_castsi128_ps(_mm_unpackhi_epi64(temp0, temp1));

        vfloat diff = vmaxf(ZEROV, indexv) - _mm_cvtepi32_ps(indexes);
        return vintpf(diff, upperVal, lowerVal);
    }

    // NOTE: This version requires LUTs which clip at upper and lower bounds
    // (which is the default).
    vfloat operator[](vfloat indexv) const
    {
        static_assert(std::is_same<T, float>::value, "This method only works for float LUTs");

        // Clamp and convert to integer values. Extract out of SSE register because all
        // lookup operations use regular addresses.
        vfloat clampedIndexes = vclampf(indexv, ZEROV, maxsv); // this automagically uses ZEROV in case indexv is NaN
        vint indexes = _mm_cvttps_epi32(clampedIndexes);
        int indexArray[4];
        _mm_storeu_si128(reinterpret_cast<__m128i*>(&indexArray[0]), indexes);

        // Load data from the table. This reads more than necessary, but there don't seem
        // to exist more granular operations (though we could try non-SSE).
        // Cast to int for convenience in the next operation (partial transpose).
        vint values[4];
        for (int i = 0; i < 4; ++i) {
            values[i] = _mm_castps_si128(LVFU(data[indexArray[i]]));
        }

        // Partial 4x4 transpose operation. We want two new vectors, the first consisting
        // of [values[0][0] ... values[3][0]] and the second [values[0][1] ... values[3][1]].
        __m128i temp0 = _mm_unpacklo_epi32(values[0], values[1]);
        __m128i temp1 = _mm_unpacklo_epi32(values[2], values[3]);
        vfloat lowerVal = _mm_castsi128_ps(_mm_unpacklo_epi64(temp0, temp1));
        vfloat upperVal = _mm_castsi128_ps(_mm_unpackhi_epi64(temp0, temp1));

        vfloat diff = vclampf(indexv, ZEROV, sizev) - _mm_cvtepi32_ps(indexes); // this automagically uses ZEROV in case indexv is NaN
        return vintpf(diff, upperVal, lowerVal);
    }

    // NOTE: This version requires LUTs which do not clip at upper and lower bounds
    vfloat operator()(vfloat indexv) const
    {
        static_assert(std::is_same<T, float>::value, "This method only works for float LUTs");

        // Clamp and convert to integer values. Extract out of SSE register because all
        // lookup operations use regular addresses.
        vfloat clampedIndexes = vclampf(indexv, ZEROV, maxsv); // this automagically uses ZEROV in case indexv is NaN
        vint indexes = _mm_cvttps_epi32(clampedIndexes);
        int indexArray[4];
        _mm_storeu_si128(reinterpret_cast<__m128i*>(&indexArray[0]), indexes);

        // Load data from the table. This reads more than necessary, but there don't seem
        // to exist more granular operations (though we could try non-SSE).
        // Cast to int for convenience in the next operation (partial transpose).
        vint values[4];
        for (int i = 0; i < 4; ++i) {
            values[i] = _mm_castps_si128(LVFU(data[indexArray[i]]));
        }

        // Partial 4x4 transpose operation. We want two new vectors, the first consisting
        // of [values[0][0] ... values[3][0]] and the second [values[0][1] ... values[3][1]].
        __m128i temp0 = _mm_unpacklo_epi32(values[0], values[1]);
        __m128i temp1 = _mm_unpacklo_epi32(values[2], values[3]);
        vfloat lowerVal = _mm_castsi128_ps(_mm_unpacklo_epi64(temp0, temp1));
        vfloat upperVal = _mm_castsi128_ps(_mm_unpackhi_epi64(temp0, temp1));

        vfloat diff = indexv - _mm_cvtepi32_ps(indexes);
        return vintpf(diff, upperVal, lowerVal);
    }

    // vectorized LUT access with integer indices. Clips at lower and upper bounds
#ifdef __SSE4_1__
    template<typename U = T, typename = typename std::enable_if<std::is_same<U, float>::value>::type>
    vfloat operator[](vint idxv) const
    {
        idxv = _mm_max_epi32( _mm_setzero_si128(), _mm_min_epi32(idxv, sizeiv));
        // access the LUT 4 times. Trust the compiler. It generates good code here, better than hand written SSE code
        return _mm_setr_ps(data[_mm_extract_epi32(idxv,0)], data[_mm_extract_epi32(idxv,1)], data[_mm_extract_epi32(idxv,2)], data[_mm_extract_epi32(idxv,3)]);
    }
#else
    template<typename U = T, typename = typename std::enable_if<std::is_same<U, float>::value>::type>
    vfloat operator[](vint idxv) const
    {
        // convert to float because SSE2 has no min/max for 32bit integers
        vfloat tempv = vclampf(_mm_cvtepi32_ps(idxv), ZEROV, sizev); // this automagically uses ZEROV in case idxv is NaN (which will never happen because it is a vector of int)
        idxv = _mm_cvttps_epi32(tempv);
        // access the LUT 4 times. Trust the compiler. It generates good code here, better than hand written SSE code
        return _mm_setr_ps(data[_mm_cvtsi128_si32(idxv)],
                           data[_mm_cvtsi128_si32(_mm_shuffle_epi32(idxv, _MM_SHUFFLE(1, 1, 1, 1)))],
                           data[_mm_cvtsi128_si32(_mm_shuffle_epi32(idxv, _MM_SHUFFLE(2, 2, 2, 2)))],
                           data[_mm_cvtsi128_si32(_mm_shuffle_epi32(idxv, _MM_SHUFFLE(3, 3, 3, 3)))]);
    }
#endif
#endif

    // use with float indices
    template<typename U = T, typename V, typename = typename std::enable_if<std::is_floating_point<V>::value && std::is_same<U, float>::value>::type>
    T operator[](V index) const
    {
        int idx = (int)index;  // don't use floor! The difference in negative space is no problems here

        if (index < 0.f) {
            if (clip & LUT_CLIP_BELOW) {
                return data[0];
            }

            idx = 0;
        } else if (index > maxsf) {
            if (clip & LUT_CLIP_ABOVE) {
                return data[upperBound];
            }

            idx = maxs;
        }

        float diff = index - (float) idx;
        T p1 = data[idx];
        T p2 = data[idx + 1] - p1;
        return (p1 + p2 * diff);
    }

    // Return the value for "index" that is in the [0-1] range.
    template<typename U = T, typename = typename std::enable_if<std::is_same<U, float>::value>::type>
    T getVal01(float index) const
    {
        index *= (float)upperBound;
        int idx = (int)index;  // don't use floor! The difference in negative space is no problems here

        if (index < 0.f) {
            if (clip & LUT_CLIP_BELOW) {
                return data[0];
            }

            idx = 0;
        } else if (index > maxsf) {
            if (clip & LUT_CLIP_ABOVE) {
                return data[upperBound];
            }

            idx = maxs;
        }

        float diff = index - (float) idx;
        T p1 = data[idx];
        T p2 = data[idx + 1] - p1;
        return (p1 + p2 * diff);
    }

    operator bool() const // FIXME: Should be explicit
    {
        return size > 0;
    }

    void clear()
    {
        if (data && size) {
            memset(data, 0, size * sizeof(T));
        }
    }

    void reset()
    {
        if (data) {
            delete[] data;
        }

        dirty = true;
        data = nullptr;
        owner = 1;
        size = 0;
        upperBound = 0;
        maxs = 0;
        maxsf = 0.f;
        clip = 0;
    }

    // create an identity LUT (LUT(x) = x) or a scaled identity LUT (LUT(x) = x / divisor)
    template<typename U = T, typename = typename std::enable_if<std::is_same<U, float>::value>::type>
    void makeIdentity(float divisor = 1.f)
    {
        if(divisor == 1.f) {
            for(unsigned int i = 0; i < size; i++) {
                data[i] = i;
            }
        } else {
            for(unsigned int i = 0; i < size; i++) {
                data[i] = i / divisor;
            }
        }
    }

    // compress a LUT<uint32_t> with size y into a LUT<uint32_t> with size x (y>x)
    template<typename U = T, typename = typename std::enable_if<std::is_same<U, std::uint32_t>::value>::type>
    void compressTo(LUT<T> &dest, unsigned int numVals = 0) const
    {
        numVals = numVals == 0 ? size : numVals;
        numVals = std::min(numVals, size);
        float divisor = numVals - 1;
        float mult = (dest.size - 1) / divisor;

        for (unsigned int i = 0; i < numVals; i++) {
            int hi = (int)(mult * i);
            dest.data[hi] += this->data[i] ;
        }
    }

    // compress a LUT<uint32_t> with size y into a LUT<uint32_t> with size x (y>x) by using the passThrough LUT to calculate indexes
    template<typename U = T, typename = typename std::enable_if<std::is_same<U, std::uint32_t>::value>::type>
    void compressTo(LUT<T> &dest, unsigned int numVals, const LUT<float> &passThrough) const
    {
        if(passThrough) {
            numVals = std::min(numVals, size);
            numVals = std::min(numVals, passThrough.getSize());
            float mult = dest.size - 1;

            for (unsigned int i = 0; i < numVals; i++) {
                int hi = (int)(mult * passThrough[i]);
                dest[hi] += this->data[i] ;
            }
        }
    }

    // compute sum and average of a LUT<uint32_t>
    template<typename U = T, typename = typename std::enable_if<std::is_same<U, std::uint32_t>::value>::type>
    void getSumAndAverage(float &sum, float &avg) const
    {
        sum = 0.f;
        avg = 0.f;
        int i = 0;
#ifdef __SSE2__
        vfloat iv = _mm_set_ps(3.f, 2.f, 1.f, 0.f);
        vfloat fourv = F2V(4.f);
        vint sumv = (vint)ZEROV;
        vfloat avgv = ZEROV;

        for(; i < static_cast<int>(size) - 3; i += 4) {
            vint datav = _mm_loadu_si128((__m128i*)&data[i]);
            sumv += datav;
            avgv += iv * _mm_cvtepi32_ps(datav);
            iv += fourv;

        }

        sum = vhadd(_mm_cvtepi32_ps(sumv));
        avg = vhadd(avgv);
#endif

        for (; i < static_cast<int>(size); i++) {
            T val = data[i];
            sum += val;
            avg += i * val;
        }

        avg /= sum;
    }


    template<typename U = T, typename = typename std::enable_if<std::is_same<U, float>::value>::type>
    void makeConstant(float value, unsigned int numVals = 0)
    {
        numVals = numVals == 0 ? size : numVals;
        numVals = std::min(numVals, size);

        for(unsigned int i = 0; i < numVals; i++) {
            data[i] = value;
        }
    }

    // share the buffer with another LUT, handy for same data but different clip flags
    void share(const LUT<T> &source, int flags = LUT_CLIP_BELOW | LUT_CLIP_ABOVE)
    {
        if (owner && data) {
            delete[] data;
        }

        dirty = false;  // Assumption
        clip = flags;
        data = source.data;
        owner = 0;
        size = source.getSize();
        upperBound = size - 1;
        maxs = size - 2;
        maxsf = (float)maxs;
#ifdef __SSE2__
        maxsv =  F2V( size - 2);
        sizeiv =  _mm_set1_epi32( (int)(size - 1) );
        sizev = F2V( size - 1 );
#endif
    }


};
