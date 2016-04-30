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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
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

#ifndef LUT_H_
#define LUT_H_

// bit representations of flags
#define LUT_CLIP_BELOW 1
#define LUT_CLIP_ABOVE 2

#define LUTf LUT<float>
#define LUTi LUT<int>
#define LUTu LUT<unsigned int>
#define LUTd LUT<double>
#define LUTuc LUT<unsigned char>

#include <cstring>
#ifndef NDEBUG
#include <glibmm.h>
#include <fstream>
#endif
#include "opthelper.h"
#include <assert.h>
#include "rt_math.h"

template<typename T>
class LUT
{
protected:
    // list of variables ordered to improve cache speed
    unsigned int maxs;
    float maxsf;
    T * data;
    unsigned int clip;
    unsigned int size;
    unsigned int upperBound;  // always equals size-1, parameter created for performance reason
private:
    unsigned int owner;
#if defined( __SSE2__ ) && defined( __x86_64__ )
    vfloat maxsv ALIGNED16;
    vfloat sizev ALIGNED16;
    vint sizeiv ALIGNED16;
#endif
public:
    /// convenience flag! If one doesn't want to delete the buffer but want to flag it to be recomputed...
    /// The user have to handle it itself, even if some method can (re)initialize it
    bool dirty;

    LUT(int s, int flags = 0xfffffff)
    {
#ifndef NDEBUG

        if (s <= 0) {
            printf("s<=0!\n");
        }

        assert (s > 0);
#endif
        dirty = true;
        clip = flags;
        data = new T[s];
        owner = 1;
        size = s;
        upperBound = size - 1;
        maxs = size - 2;
        maxsf = (float)maxs;
#if defined( __SSE2__ ) && defined( __x86_64__ )
        maxsv =  F2V( maxs );
        sizeiv =  _mm_set1_epi32( (int)(size - 1) );
        sizev = F2V( size - 1 );
#endif
    }
    void operator ()(int s, int flags = 0xfffffff)
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
        data = new T[s];
        owner = 1;
        size = s;
        upperBound = size - 1;
        maxs = size - 2;
        maxsf = (float)maxs;
#if defined( __SSE2__ ) && defined( __x86_64__ )
        maxsv =  F2V( maxs );
        sizeiv =  _mm_set1_epi32( (int)(size - 1) );
        sizev = F2V( size - 1 );
#endif
    }

    LUT(int s, T * source, int flags = 0xfffffff)
    {
#ifndef NDEBUG

        if (s <= 0) {
            printf("s<=0!\n");
        }

        assert (s > 0);

        if (!source) {
            printf("source is NULL!\n");
        }

        assert (source != nullptr);
#endif
        dirty = false;  // Assumption
        clip = flags;
        data = new T[s];
        owner = 1;
        size = s;
        upperBound = size - 1;
        maxs = size - 2;
        maxsf = (float)maxs;
#if defined( __SSE2__ ) && defined( __x86_64__ )
        maxsv =  F2V( size - 2);
        sizeiv =  _mm_set1_epi32( (int)(size - 1) );
        sizev = F2V( size - 1 );
#endif

        for (int i = 0; i < s; i++) {
            data[i] = source[i];
        }
    }

    LUT()
    {
        data = nullptr;
        reset();
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

    void setClip(int flags)
    {
        clip = flags;
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
    int getUpperBound()
    {
        return size > 0 ? upperBound : 0;
    }

    LUT<T> & operator=(LUT<T> &rhs)
    {
        if (this != &rhs) {
            if (rhs.size > this->size) {
                delete [] this->data;
                this->data = nullptr;
            }

            if (this->data == nullptr) {
                this->data = new T[rhs.size];
            }

            this->clip = rhs.clip;
            this->owner = 1;
            memcpy(this->data, rhs.data, rhs.size * sizeof(T));
            this->size = rhs.size;
            this->upperBound = rhs.upperBound;
            this->maxs = this->size - 2;
            this->maxsf = (float)this->maxs;
#if defined( __SSE2__ ) && defined( __x86_64__ )
            this->maxsv =  F2V( this->size - 2);
            this->sizeiv =  _mm_set1_epi32( (int)(this->size - 1) );
            this->sizev = F2V( this->size - 1 );
#endif
        }

        return *this;
    }

    // handy to sum up per thread histograms. #pragma omp simd speeds up the loop by about factor 3 for LUTu (unsigned int).
    LUT<T> & operator+=(LUT<T> &rhs)
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

    // mutliply all elements of LUT with a constant float value
    LUT<float> & operator*=(float factor)
    {
#ifdef _OPENMP
        #pragma omp simd
#endif

        for(unsigned int i = 0; i < this->size; i++) {
            data[i] *= factor;
        }

        return *this;
    }

    // use with integer indices
    T& operator[](int index) const
    {
        return data[ rtengine::LIM<int>(index, 0, upperBound) ];
    }

#if defined( __SSE2__ ) && defined( __x86_64__ )
/*
    vfloat operator[](vfloat indexv ) const
    {
//      printf("don't use this operator. It's not ready for production");
        return _mm_setzero_ps();

        // convert floats to ints
        vint idxv =  _mm_cvttps_epi32( indexv );
        vfloat tempv, resultv, p1v, p2v;
        vmask maxmask = vmaskf_gt(indexv, maxsv);
        idxv = _mm_castps_si128(vself(maxmask, maxsv, _mm_castsi128_ps(idxv)));
        vmask minmask = vmaskf_lt(indexv, _mm_setzero_ps());
        idxv = _mm_castps_si128(vself(minmask, _mm_setzero_ps(), _mm_castsi128_ps(idxv)));
        // access the LUT 4 times and shuffle the values into p1v and p2v

        int idx;

        // get 4th value
        idx = _mm_cvtsi128_si32 (_mm_shuffle_epi32(idxv, _MM_SHUFFLE(3, 3, 3, 3)));
        tempv = LVFU(data[idx]);
        p1v = _mm_shuffle_ps(tempv, tempv, _MM_SHUFFLE(0, 0, 0, 0));
        p2v = _mm_shuffle_ps(tempv, tempv, _MM_SHUFFLE(1, 1, 1, 1));
        // now p1v is 3 3 3 3
        //     p2v is 3 3 3 3

        // get 3rd value
        idx = _mm_cvtsi128_si32 (_mm_shuffle_epi32(idxv, _MM_SHUFFLE(2, 2, 2, 2)));
        tempv = LVFU(data[idx]);
        p1v = _mm_move_ss( p1v, tempv);
        tempv = _mm_shuffle_ps(tempv, tempv, _MM_SHUFFLE(1, 1, 1, 1));
        p2v = _mm_move_ss( p2v, tempv);
        // now p1v is 3 3 3 2
        //     p2v is 3 3 3 2

        // get 2nd value
        idx = _mm_cvtsi128_si32 (_mm_shuffle_epi32(idxv, _MM_SHUFFLE(1, 1, 1, 1)));
        tempv = LVFU(data[idx]);
        p1v = _mm_shuffle_ps( p1v, p1v, _MM_SHUFFLE(1, 0, 1, 0));
        p2v = _mm_shuffle_ps( p2v, p2v, _MM_SHUFFLE(1, 0, 1, 0));
        // now p1v is 3 2 3 2
        // now p2v is 3 2 3 2
        p1v = _mm_move_ss( p1v, tempv );
        // now p1v is 3 2 3 1
        tempv = _mm_shuffle_ps(tempv, tempv, _MM_SHUFFLE(1, 1, 1, 1));
        p2v = _mm_move_ss( p2v, tempv);
        // now p1v is 3 2 3 1

        // get 1st value
        idx = _mm_cvtsi128_si32 (_mm_shuffle_epi32(idxv, _MM_SHUFFLE(0, 0, 0, 0)));
        tempv = LVFU(data[idx]);
        p1v = _mm_shuffle_ps( p1v, p1v, _MM_SHUFFLE(3, 2, 0, 0));
        // now p1v is 3 2 1 1
        p2v = _mm_shuffle_ps( p2v, p2v, _MM_SHUFFLE(3, 2, 0, 0));
        // now p2v is 3 2 1 1
        p1v = _mm_move_ss( p1v, tempv );
        // now p1v is 3 2 1 0
        tempv = _mm_shuffle_ps(tempv, tempv, _MM_SHUFFLE(1, 1, 1, 1));
        p2v = _mm_move_ss( p2v, tempv);
        // now p2v is 3 2 1 0

        vfloat diffv = indexv - _mm_cvtepi32_ps ( idxv );
        diffv = vself(vorm(maxmask, minmask), _mm_setzero_ps(), diffv);
        resultv = p1v + p2v * diffv;
        return resultv  ;
    }
*/
#ifdef __SSE4_1__
    vfloat operator[](vint idxv ) const
    {
        vfloat tempv, p1v;
        idxv = _mm_max_epi32( _mm_setzero_si128(), _mm_min_epi32(idxv, sizeiv));
        // access the LUT 4 times and shuffle the values into p1v

        int idx;

        // get 4th value
        idx = _mm_extract_epi32(idxv, 3);
        tempv = _mm_load_ss(&data[idx]);
        p1v = PERMUTEPS(tempv, _MM_SHUFFLE(0, 0, 0, 0));
        // now p1v is 3 3 3 3

        // get 3rd value
        idx = _mm_extract_epi32(idxv, 2);
        tempv = _mm_load_ss(&data[idx]);
        p1v = _mm_move_ss( p1v, tempv);
        // now p1v is 3 3 3 2

        // get 2nd value
        idx = _mm_extract_epi32(idxv, 1);
        tempv = _mm_load_ss(&data[idx]);
        p1v = PERMUTEPS( p1v, _MM_SHUFFLE(1, 0, 1, 0));
        // now p1v is 3 2 3 2
        p1v = _mm_move_ss( p1v, tempv );
        // now p1v is 3 2 3 1

        // get 1st value
        idx = _mm_cvtsi128_si32(idxv);
        tempv = _mm_load_ss(&data[idx]);
        p1v = PERMUTEPS( p1v, _MM_SHUFFLE(3, 2, 0, 0));
        // now p1v is 3 2 1 1
        p1v = _mm_move_ss( p1v, tempv );
        // now p1v is 3 2 1 0

        return p1v;
    }
#else
    vfloat operator[](vint idxv ) const
    {
        vfloat tempv, p1v;
        tempv = _mm_cvtepi32_ps(idxv);
        tempv = _mm_min_ps( tempv, sizev );
        idxv = _mm_cvttps_epi32(_mm_max_ps( tempv, _mm_setzero_ps( )  ));
        // access the LUT 4 times and shuffle the values into p1v

        int idx;

        // get 4th value
        idx = _mm_cvtsi128_si32 (_mm_shuffle_epi32(idxv, _MM_SHUFFLE(3, 3, 3, 3)));
        tempv = _mm_load_ss(&data[idx]);
        p1v = PERMUTEPS(tempv, _MM_SHUFFLE(0, 0, 0, 0));
        // now p1v is 3 3 3 3

        // get 3rd value
        idx = _mm_cvtsi128_si32 (_mm_shuffle_epi32(idxv, _MM_SHUFFLE(2, 2, 2, 2)));
        tempv = _mm_load_ss(&data[idx]);
        p1v = _mm_move_ss( p1v, tempv);
        // now p1v is 3 3 3 2

        // get 2nd value
        idx = _mm_cvtsi128_si32 (_mm_shuffle_epi32(idxv, _MM_SHUFFLE(1, 1, 1, 1)));
        tempv = _mm_load_ss(&data[idx]);
        p1v = PERMUTEPS( p1v, _MM_SHUFFLE(1, 0, 1, 0));
        // now p1v is 3 2 3 2
        p1v = _mm_move_ss( p1v, tempv );
        // now p1v is 3 2 3 1

        // get 1st value
        idx = _mm_cvtsi128_si32 (idxv);
        tempv = _mm_load_ss(&data[idx]);
        p1v = PERMUTEPS( p1v, _MM_SHUFFLE(3, 2, 0, 0));
        // now p1v is 3 2 1 1
        p1v = _mm_move_ss( p1v, tempv );
        // now p1v is 3 2 1 0

        return p1v;
    }
#endif
#endif

    // use with float indices
    T operator[](float index) const
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
    T getVal01 (float index) const
    {
        index *= float(upperBound);
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

#ifndef NDEBUG
    // Debug facility ; dump the content of the LUT in a file. No control of the filename is done
    void dump(Glib::ustring fname)
    {
        if (size) {
            Glib::ustring fname_ = fname + ".xyz"; // TopSolid'Design "plot" file format
            std::ofstream f (fname_.c_str());
            f << "$" << std::endl;

            for (unsigned int iter = 0; iter < size; iter++) {
                f << iter << ", " << data[iter] << ", 0." << std::endl;
            }

            f << "$" << std::endl;
            f.close ();
        }
    }
#endif


    operator bool (void) const
    {
        return size > 0;
    }

    void clear(void)
    {
        if (data && size) {
            memset(data, 0, size * sizeof(T));
        }
    }

    void reset(void)
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
    }

    // create an identity LUT (LUT(x) = x) or a scaled identity LUT (LUT(x) = x / divisor)
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

    // compress a LUT with size y into a LUT with size y (y<x)
    void compressTo(LUT<T> &dest, unsigned int numVals) const
    {
        numVals = std::min(numVals, size);
        float divisor = numVals - 1;
        float mult = (dest.size - 1) / divisor;

        for (unsigned int i = 0; i < numVals; i++) {
            int hi = (int)(mult * i);
            dest.data[hi] += this->data[i] ;
        }
    }

    // compress a LUT with size y into a LUT with size y (y<x) by using the passTrough LUT to calculate indexes
    void compressTo(LUT<T> &dest, unsigned int numVals, const LUT<float> &passThrough) const
    {
        if(passThrough) {
            numVals = std::min(numVals, size);
            numVals = std::min(numVals, passThrough.getSize());
            float mult = dest.size - 1;

            for (int i = 0; i < numVals; i++) {
                int hi = (int)(mult * passThrough[i]);
                dest[hi] += this->data[i] ;
            }
        }
    }

    void makeConstant(float value, unsigned int numVals = 0)
    {
        numVals = numVals == 0 ? size : numVals;
        numVals = std::min(numVals, size);
        for(unsigned int i = 0; i < numVals; i++) {
            data[i] = value;
        }
    }


};

#endif /* LUT_H_ */
