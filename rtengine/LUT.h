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
 *  	LUT<type> name (size);
 *		LUT<type> name (size, flags);
 *
 *		creates an array which is valid within the normal C/C++ scope "{ ... }"
 *
 *      access to elements is a simple as:
 *
 *      	LUT<float> my_lut (10);
 *      	float value = my_lut[3];
 *          float value = my_lut[2.5]; // this will interpolate
 *
 *      when using a float type index it will interpolate the lookup values
 *
 *      extra setting in flags: (clipping is set by default)
 *      LUT_CLIP_ABOVE
 *      LUT_CLIP_BELOW
 *
 *      example:
 *      	LUT<float> my_lut (10,LUT_CLIP_BELOW);
 *          float value = my_lut[22.5];  // this will extrapolate
 *          float value = my_lut[-22.5]; // this will not extrapolate
 *
 *          LUT<float> my_lut (10,0); // this will extrapolate on either side
 *
 *      shotcuts:
 *
 *      	LUTf stands for LUT<float>
 *          LUTi stands for LUT<int>
 *          LUTu stands for LUT<unsigned int>
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

#include <cstring>
#ifndef NDEBUG
#include <glibmm.h>
#include <fstream>
#endif
#ifdef __SSE2__
#include "sleefsseavx.c"
#endif
#include <assert.h>
#include "rt_math.h"

template<typename T>
class LUT {
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
	__m128 maxsv __attribute__ ((aligned (16)));
	__m128 sizev __attribute__ ((aligned (16)));
	__m128i maxsiv __attribute__ ((aligned (16)));
	__m128i sizeiv __attribute__ ((aligned (16)));
#endif
public:
	/// convenience flag! If one doesn't want to delete the buffer but want to flag it to be recomputed...
	/// The user have to handle it itself, even if some method can (re)initialize it
	bool dirty;

	LUT(int s, int flags = 0xfffffff) {
		#ifndef NDEBUG
		if (s<=0)
			printf("s<=0!\n");
		assert (s>0);
		#endif
		dirty = true;
		clip = flags;
		data = new T[s];
		owner = 1;
		size = s;
		upperBound = size-1;
		maxs=size-2;
		maxsf = (float)maxs;
#if defined( __SSE2__ ) && defined( __x86_64__ )
		maxsv =  _mm_set1_ps( maxs );
		maxsiv = _mm_cvttps_epi32( maxsv );
		sizeiv =  _mm_set1_epi32( (int)(size-1) );
		sizev = _mm_set1_ps( size-1 );
#endif
	}
	void operator ()(int s, int flags = 0xfffffff) {
		#ifndef NDEBUG
		if (s<=0)
			printf("s<=0!\n");
		assert (s>0);
		#endif
		if (owner&&data)
			delete[] data;
		dirty = true; // Assumption!
		clip = flags;
		data = new T[s];
		owner = 1;
		size = s;
		upperBound = size-1;
		maxs=size-2;
		maxsf = (float)maxs;
#if defined( __SSE2__ ) && defined( __x86_64__ )
		maxsv =  _mm_set1_ps( maxs );
		maxsiv = _mm_cvttps_epi32( maxsv );
		sizeiv =  _mm_set1_epi32( (int)(size-1) );
		sizev = _mm_set1_ps( size-1 );
#endif
	}

	LUT(int s, T * source, int flags = 0xfffffff) {
		#ifndef NDEBUG
		if (s<=0)
			printf("s<=0!\n");
		assert (s>0);
		if (source==NULL)
			printf("source is NULL!\n");
		assert (source != NULL);
		#endif
		dirty = false;  // Assumption
		clip = flags;
		data = new T[s];
		owner = 1;
		size = s;
		upperBound = size-1;
		maxs=size-2;
		maxsf = (float)maxs;
#if defined( __SSE2__ ) && defined( __x86_64__ )
		maxsv =  _mm_set1_ps( size - 2);
		maxsiv = _mm_cvttps_epi32( maxsv );
		sizeiv =  _mm_set1_epi32( (int)(size-1) );
		sizev = _mm_set1_ps( size-1 );
#endif
		for (int i = 0; i < s; i++) {
			data[i] = source[i];
		}
	}

	LUT() {
		data = NULL;
		reset();
	}

	~LUT() {
		if (owner) {
			delete[] data;
			#ifndef NDEBUG
			data=(T*)0xBAADF00D;
			#endif
		}
	}

	void setClip(int flags) {
		clip = flags;
	}

	/** @brief Get the number of element in the LUT (i.e. dimension of the array)
	 *  For a LUT(500), it will return 500
	 *  @return number of element in the array
	 */
	int getSize() {
		return size;
	}

	/** @brief Get the highest value possible (i.e. dimension of the array)
	 *  For a LUT(500), it will return 499, because 500 elements, starting from 0, goes up to 499
	 *  @return number of element in the array
	 */
	int getUpperBound() {
		return size>0 ? upperBound : 0;
	}

	LUT<T> & operator=(LUT<T> &rhs) {
	    if (this != &rhs) {
	      if (rhs.size>this->size)
	      {
	    	delete [] this->data;
	    	this->data=NULL;
	      }
	      if (this->data==NULL) this->data=new T[rhs.size];
	      this->clip=rhs.clip;
	      this->owner=1;
	      memcpy(this->data,rhs.data,rhs.size*sizeof(T));
	      this->size=rhs.size;
	      this->upperBound=rhs.upperBound;
	      this->maxs=this->size-2;
		  this->maxsf = (float)this->maxs;
#if defined( __SSE2__ ) && defined( __x86_64__ )
		  this->maxsv =  _mm_set1_ps( this->size - 2);
		  this->maxsiv = _mm_cvttps_epi32( this->maxsv );
		  this->sizeiv =  _mm_set1_epi32( (int)(this->size-1) );
		  this->sizev = _mm_set1_ps( this->size-1 );
#endif
	    }

	    return *this;
	  }
	// use with integer indices
	T& operator[](int index) const {
		return data[ rtengine::LIM<int>(index, 0, upperBound) ];
	}

#if defined( __SSE2__ ) && defined( __x86_64__ )
	__m128 operator[](__m128 indexv ) const {
		printf("don't use this operator. It's not ready for production");
		return _mm_setzero_ps();

		// convert floats to ints
		__m128i	idxv =  _mm_cvttps_epi32( indexv );
		__m128 tempv, resultv, p1v, p2v;
		vmask maxmask = vmaskf_gt(indexv, maxsv);
		idxv = _mm_castps_si128(vself(maxmask, maxsv, _mm_castsi128_ps(idxv)));
		vmask minmask = vmaskf_lt(indexv, _mm_setzero_ps());
		idxv = _mm_castps_si128(vself(minmask, _mm_setzero_ps(), _mm_castsi128_ps(idxv)));
		// access the LUT 4 times and shuffle the values into p1v and p2v

		int idx;

		// get 4th value
		idx = _mm_cvtsi128_si32 (_mm_shuffle_epi32(idxv,_MM_SHUFFLE(3,3,3,3)));
		tempv = LVFU(data[idx]);
		p1v = _mm_shuffle_ps(tempv, tempv, _MM_SHUFFLE(0,0,0,0));
		p2v = _mm_shuffle_ps(tempv, tempv, _MM_SHUFFLE(1,1,1,1));
		// now p1v is 3 3 3 3
		//     p2v is 3 3 3 3

		// get 3rd value
		idx = _mm_cvtsi128_si32 (_mm_shuffle_epi32(idxv,_MM_SHUFFLE(2,2,2,2)));
		tempv = LVFU(data[idx]);
		p1v = _mm_move_ss( p1v, tempv);
		tempv = _mm_shuffle_ps(tempv, tempv, _MM_SHUFFLE(1,1,1,1));
		p2v = _mm_move_ss( p2v, tempv);
		// now p1v is 3 3 3 2
		//     p2v is 3 3 3 2

		// get 2nd value
		idx = _mm_cvtsi128_si32 (_mm_shuffle_epi32(idxv,_MM_SHUFFLE(1,1,1,1)));
		tempv = LVFU(data[idx]);
		p1v = _mm_shuffle_ps( p1v, p1v, _MM_SHUFFLE(1,0,1,0));
		p2v = _mm_shuffle_ps( p2v, p2v, _MM_SHUFFLE(1,0,1,0));
		// now p1v is 3 2 3 2
		// now p2v is 3 2 3 2
		p1v = _mm_move_ss( p1v, tempv );
		// now p1v is 3 2 3 1
		tempv = _mm_shuffle_ps(tempv, tempv, _MM_SHUFFLE(1,1,1,1));
		p2v = _mm_move_ss( p2v, tempv);
		// now p1v is 3 2 3 1

		// get 1st value
		idx = _mm_cvtsi128_si32 (_mm_shuffle_epi32(idxv,_MM_SHUFFLE(0,0,0,0)));
		tempv = LVFU(data[idx]);
		p1v = _mm_shuffle_ps( p1v, p1v, _MM_SHUFFLE(3,2,0,0));
		// now p1v is 3 2 1 1
		p2v = _mm_shuffle_ps( p2v, p2v, _MM_SHUFFLE(3,2,0,0));
		// now p2v is 3 2 1 1
		p1v = _mm_move_ss( p1v, tempv );
		// now p1v is 3 2 1 0
		tempv = _mm_shuffle_ps(tempv, tempv, _MM_SHUFFLE(1,1,1,1));
		p2v = _mm_move_ss( p2v, tempv);
		// now p2v is 3 2 1 0

		__m128 diffv = indexv - _mm_cvtepi32_ps ( idxv );
		diffv = vself(vorm(maxmask,minmask), _mm_setzero_ps(), diffv);
		resultv = p1v + p2v * diffv;
		return resultv	;
	}

	__m128 operator[](__m128i idxv ) const
	 {
		__m128 tempv, p1v;
		tempv = _mm_cvtepi32_ps(idxv);
		tempv = _mm_min_ps( tempv, sizev );
		idxv = _mm_cvttps_epi32(_mm_max_ps( tempv, _mm_setzero_ps( )  ));
		// access the LUT 4 times and shuffle the values into p1v

		int idx;

		// get 4th value
		idx = _mm_cvtsi128_si32 (_mm_shuffle_epi32(idxv,_MM_SHUFFLE(3,3,3,3)));
		tempv = _mm_load_ss(&data[idx]);
		p1v = _mm_shuffle_ps(tempv, tempv, _MM_SHUFFLE(0,0,0,0));
		// now p1v is 3 3 3 3

		// get 3rd value
		idx = _mm_cvtsi128_si32 (_mm_shuffle_epi32(idxv,_MM_SHUFFLE(2,2,2,2)));
		tempv = _mm_load_ss(&data[idx]);
		p1v = _mm_move_ss( p1v, tempv);
		// now p1v is 3 3 3 2

		// get 2nd value
		idx = _mm_cvtsi128_si32 (_mm_shuffle_epi32(idxv,_MM_SHUFFLE(1,1,1,1)));
		tempv = _mm_load_ss(&data[idx]);
		p1v = _mm_shuffle_ps( p1v, p1v, _MM_SHUFFLE(1,0,1,0));
		// now p1v is 3 2 3 2
		p1v = _mm_move_ss( p1v, tempv );
		// now p1v is 3 2 3 1

		// get 1st value
		idx = _mm_cvtsi128_si32 (idxv);
		tempv = _mm_load_ss(&data[idx]);
		p1v = _mm_shuffle_ps( p1v, p1v, _MM_SHUFFLE(3,2,0,0));
		// now p1v is 3 2 1 1
		p1v = _mm_move_ss( p1v, tempv );
		// now p1v is 3 2 1 0

		return p1v;
	}
#endif

	// use with float indices
	T operator[](float index) const {
		int idx = (int)index;  // don't use floor! The difference in negative space is no problems here
		if (index<0.f)
		{
			if (clip & LUT_CLIP_BELOW)
				return data[0];
			idx=0;
		}
		else if (index > maxsf)
		{
			if (clip & LUT_CLIP_ABOVE)
				return data[upperBound];
			idx =maxs;
		}
		float diff = index - (float) idx;
		T p1 = data[idx];
		T p2 = data[idx + 1]-p1;
		return (p1 + p2*diff);
	}

	// Return the value for "index" that is in the [0-1] range.
	T getVal01 (float index) const {
		index *= float(upperBound);
		int idx = (int)index;  // don't use floor! The difference in negative space is no problems here
		if (index<0.f)
		{
			if (clip & LUT_CLIP_BELOW)
				return data[0];
			idx=0;
		}
		else if (index > maxsf)
		{
			if (clip & LUT_CLIP_ABOVE)
				return data[upperBound];
			idx =maxs;
		}
		float diff = index - (float) idx;
		T p1 = data[idx];
		T p2 = data[idx + 1]-p1;
		return (p1 + p2*diff);
	}

#ifndef NDEBUG
	// Debug facility ; dump the content of the LUT in a file. No control of the filename is done
	void dump(Glib::ustring fname) {
		if (size) {
		    Glib::ustring fname_ = fname + ".xyz"; // TopSolid'Design "plot" file format
			std::ofstream f (fname_.c_str());
			f << "$" << std::endl;
			for (unsigned int iter=0; iter<size; iter++) {
				f << iter << ", " << data[iter] << ", 0." << std::endl;
			}
			f << "$" << std::endl;
			f.close ();
		}
	}
#endif


	operator bool (void) const
		{
			return size>0;
		}

	void clear(void) {
		if (data && size)
			memset(data, 0, size * sizeof(T));
	}

	void reset(void) {
		if (data) delete[] data;
		dirty = true;
		data = NULL;
		owner = 1;
		size = 0;
		upperBound=0;
		maxs=0;
    }
};



// TODO: HOMBRE: HueLUT is actually unused, could we delete this class now that LUT::getVal01 has been created?


/** @brief LUT subclass handling hue values specifically.
    The array has a fixed size of float values and have to be in the [0.; 1.] range in both axis (no error checking implemented) */
class HueLUT : public LUTf {
	public:
		HueLUT() : LUTf() {}
		HueLUT(bool createArray) : LUTf() {
			if (createArray)
				this->operator () (501, LUT_CLIP_BELOW|LUT_CLIP_ABOVE);
		}

		void create() {
			this->operator () (501, LUT_CLIP_BELOW|LUT_CLIP_ABOVE);
		}

		// use with integer indices
		float& operator[](int index) const {
			return data[ rtengine::LIM<int>(index, 0, upperBound) ];
		}

		// use with float indices in the [0.;1.] range
		float operator[](float index) const {
			int idx = int(index*500.f);  // don't use floor! The difference in negative space is no problems here
			if (index<0.f)
				return data[0];
			else if (index > 1.f)
				return data[upperBound];

			float balance = index - float(idx/500.f);
			float h1 = data[idx];
			float h2 = data[idx + 1];

			if (h1==h2)
				return h1;
			if ((h1 > h2) && (h1-h2 > 0.5f)){
				h1 -= 1.f;
				float value = h1 + balance * (h2-h1);
				if (value < 0.f)
					value += 1.f;
				return value;
			}
			else if (h2-h1 > 0.5f) {
				h2 -= 1.f;
				float value = h1 + balance * (h2-h1);
				if (value < 0.f)
					value += 1.f;
				return value;
			}
			else
				return h1 + balance * (h2-h1);
		}
};


#endif /* LUT_H_ */
