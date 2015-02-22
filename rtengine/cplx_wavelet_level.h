/*
 *  This file is part of RawTherapee.
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
 *
 *  2010 Ilya Popov <ilia_popov@rambler.ru>
 *  2012 Emil Martinec <ejmartin@uchicago.edu>
 *  2014 Ingo Weyrich <heckflosse@i-weyrich.de>
 */

#ifndef CPLX_WAVELET_LEVEL_H_INCLUDED
#define CPLX_WAVELET_LEVEL_H_INCLUDED

#include <cstddef>
#include "rt_math.h"
#include "opthelper.h"
#include "stdio.h"
namespace rtengine {

	template<typename T>
	class wavelet_level
	{

		// level of decomposition
		int lvl;

		// whether to subsample the output
		bool subsamp_out;
		
		int numThreads;

		// spacing of filter taps
		int skip;

		bool bigBlockOfMemory;
		// allocation and destruction of data storage
		T ** create(int n);
		void destroy(T ** subbands);

		// load a row/column of input data, possibly with padding

		void AnalysisFilterHaarVertical (const T * const srcbuffer, T * dstLo, T * dstHi, const int width, const int height, const int row);
		void AnalysisFilterHaarHorizontal (const T * const srcbuffer, T * dstLo, T * dstHi, const int width, const int row);
		void SynthesisFilterHaarHorizontal (const T * const srcLo, const T * const srcHi, T * dst, const int width, const int height);
		void SynthesisFilterHaarVertical (const T * const srcLo, const T * const srcHi, T * dst, const int width, const int height);

		void AnalysisFilterSubsampHorizontal (T * srcbuffer, T * dstLo, T * dstHi, float *filterLo, float *filterHi, 
							 const int taps, const int offset, const int srcwidth, const int dstwidth, const int row);
#ifdef __SSE2__
		void AnalysisFilterSubsampVertical (T * srcbuffer, T * dstLo, T * dstHi, float (*filterLo)[4], float (*filterHi)[4],
							 const int taps, const int offset, const int width, const int height, const int row);
#else
		void AnalysisFilterSubsampVertical (T * srcbuffer, T * dstLo, T * dstHi, float *filterLo, float *filterHi, 
							 int const taps, const int offset, const int width, const int height, const int row);
#endif
		void SynthesisFilterSubsampHorizontal (T * srcLo, T * srcHi, T * dst,
							  float *filterLo, float *filterHi, const int taps, const int offset, const int scrwidth, const int dstwidth, const int height);
#ifdef __SSE2__
		void SynthesisFilterSubsampVertical (T * srcLo, T * srcHi, T * dst, float (*filterLo)[4], float (*filterHi)[4], const int taps, const int offset, const int width, const int srcheight, const int dstheight, const float blend);
#else
		void SynthesisFilterSubsampVertical (T * srcLo, T * srcHi, T * dst, float *filterLo, float *filterHi, const int taps, const int offset, const int width, const int srcheight, const int dstheight, const float blend);
#endif
	public:
		bool memoryAllocationFailed;

		T ** wavcoeffs;
		// full size
		int m_w, m_h;

		// size of low frequency part
		int m_w2, m_h2;

		template<typename E>
		wavelet_level(E * src, E * dst, int level, int subsamp, int w, int h, float *filterV, float *filterH, int len, int offset, int skipcrop, int numThreads)
		: lvl(level), subsamp_out((subsamp>>level)&1), numThreads(numThreads), skip(1<<level), bigBlockOfMemory(true), memoryAllocationFailed(false), wavcoeffs(NULL), m_w(w), m_h(h), m_w2(w), m_h2(h)
		{
			if (subsamp) {
				skip = 1;
				for (int n=0; n<level; n++) {
					skip *= 2-((subsamp>>n)&1);
				}
				skip /= skipcrop;
				if(skip < 1) skip=1;

			}
			m_w2 = (subsamp_out ? (w+1)/2 : w);
			m_h2 = (subsamp_out ? (h+1)/2 : h);
			
			wavcoeffs = create((m_w2)*(m_h2));
			if(!memoryAllocationFailed)
				decompose_level(src, dst, filterV, filterH, len, offset);
			
		}

		~wavelet_level() {
			destroy(wavcoeffs);
		}

		T ** subbands() const {
			return wavcoeffs;
		}

		T * lopass() const {
			return wavcoeffs[0];
		}

		int width() const {
			return m_w2;
		}

		int height() const {
			return m_h2;
		}

		int stride() const {
			return skip;
		}
		
		bool bigBlockOfMemoryUsed() const {
			return bigBlockOfMemory;
		}

		template<typename E>
		void decompose_level(E *src, E *dst, float *filterV, float *filterH, int len, int offset);

		template<typename E>
		void reconstruct_level(E* tmpLo, E* tmpHi, E *src, E *dst, float *filterV, float *filterH, int taps, int offset, const float blend = 1.f);
	};

	template<typename T>
	T ** wavelet_level<T>::create(int n)	{
		T * data = new (std::nothrow) T[3*n];
		if(data == NULL) {
			bigBlockOfMemory = false;
		}
		T ** subbands = new T*[4];
		for(int j = 1; j < 4; j++) {
			if(bigBlockOfMemory)
				subbands[j] = data + n * (j-1);
			else {
				subbands[j] = new (std::nothrow) T[n];
				if(subbands[j] == NULL) {
					printf("Couldn't allocate memory in level %d of wavelet\n",lvl);
					memoryAllocationFailed = true;
				}
			}
		}
		return subbands;
	}

	template<typename T>
	void wavelet_level<T>::destroy(T ** subbands) {
		if(subbands) {
			if(bigBlockOfMemory)
				delete[] subbands[1];
			else {
				for(int j = 1; j < 4; j++) {
					if(subbands[j] != NULL)
						delete[] subbands[j];
				}
			}
			delete[] subbands;
		}
	}

	template<typename T>
	void wavelet_level<T>::AnalysisFilterHaarHorizontal (const T * const RESTRICT srcbuffer, T * RESTRICT dstLo, T * RESTRICT dstHi, const int width, const int row) {
		/* Basic convolution code
		 * Applies a Haar filter 
		*/							
			for(int i = 0; i < (width - skip); i++) {
				dstLo[row*width+i] = (srcbuffer[i] + srcbuffer[i+skip]);
				dstHi[row*width+i] = (srcbuffer[i] - srcbuffer[i+skip]);
			}
			for(int i = max(width-skip,skip); i < (width); i++) {
				dstLo[row*width+i] = (srcbuffer[i] + srcbuffer[i-skip]);
				dstHi[row*width+i] = (srcbuffer[i] - srcbuffer[i-skip]);
			}
	}

	template<typename T> void wavelet_level<T>::AnalysisFilterHaarVertical (const T * const RESTRICT srcbuffer, T * RESTRICT dstLo, T * RESTRICT dstHi, const int width, const int height, const int row) {
	/* Basic convolution code
	 * Applies a Haar filter 
	*/
		if(row < (height - skip)) {
			for(int j=0;j<width;j++) {
				dstLo[j] = 0.25f*(srcbuffer[row*width+j] + srcbuffer[(row+skip)*width+j]);
				dstHi[j] = 0.25f*(srcbuffer[row*width+j] - srcbuffer[(row+skip)*width+j]);
			}
		} else if(row>=max(height-skip,skip)) {
			for(int j=0;j<width;j++) {
				dstLo[j] = 0.25f*(srcbuffer[row*width+j] + srcbuffer[(row-skip)*width+j]);
				dstHi[j] = 0.25f*(srcbuffer[row*width+j] - srcbuffer[(row-skip)*width+j]);
			}
		}
	}

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	
	template<typename T> void wavelet_level<T>::SynthesisFilterHaarHorizontal (const T * const RESTRICT srcLo, const T * const RESTRICT srcHi, T * RESTRICT dst, const int width, const int height) {

		/* Basic convolution code
		 * Applies a Haar filter 
		 *
		 */
#ifdef _OPENMP
#pragma omp parallel for num_threads(numThreads) if(numThreads>1)
#endif
		for (int k=0; k<height; k++) {
			for(int i = 0; i < skip; i++) {
				dst[k*width+i] = (srcLo[k*width+i] + srcHi[k*width+i]);			
			}
			for(int i = skip; i < width; i++) {
				dst[k*width+i] = 0.5f*(srcLo[k*width+i] + srcHi[k*width+i] + srcLo[k*width+i-skip] - srcHi[k*width+i-skip]);			
			}
		}
	}

	template<typename T> void wavelet_level<T>::SynthesisFilterHaarVertical (const T * const RESTRICT srcLo, const T * const RESTRICT srcHi, T * RESTRICT dst, const int width, const int height) {

		/* Basic convolution code
		 * Applies a Haar filter 
		 *
		 */
#ifdef _OPENMP
#pragma omp parallel num_threads(numThreads) if(numThreads>1)
#endif
{
#ifdef _OPENMP
#pragma omp for nowait
#endif
		for(int i = 0; i < skip; i++) {
			for(int j=0;j<width;j++)
				dst[width*i+j] = (srcLo[i*width+j] + srcHi[i*width+j]);			
		}
#ifdef _OPENMP
#pragma omp for
#endif
		for(int i = skip; i < height; i++) {
			for(int j=0;j<width;j++)
				dst[width*i+j] = 0.5f*(srcLo[i*width+j] + srcHi[i*width+j] + srcLo[(i-skip)*width+j] - srcHi[(i-skip)*width+j]);			
		}
}
	}

	template<typename T>
	void wavelet_level<T>::AnalysisFilterSubsampHorizontal (T * RESTRICT srcbuffer, T * RESTRICT dstLo, T * RESTRICT dstHi, float * RESTRICT filterLo, float *RESTRICT filterHi, 
												const int taps, const int offset, const int srcwidth, const int dstwidth, const int row) {
		/* Basic convolution code
		 * Applies an FIR filter 'filter' with filter length 'taps', 
		 * aligning the 'offset' element of the filter with
		 * the input pixel, and skipping 'skip' pixels between taps 
		 * Output is subsampled by two
		 */
		// calculate coefficients
		for(int i = 0; i < srcwidth; i+=2) {
			float lo = 0.f, hi = 0.f;
			if (LIKELY(i>skip*taps && i<srcwidth-skip*taps)) {//bulk
				for (int j=0, l=-skip*offset; j<taps; j++, l+=skip) {
					float src = srcbuffer[i-l];
					lo += filterLo[j] * src;//lopass channel
					hi += filterHi[j] * src;//hipass channel
				}
			} else {
				for (int j=0; j<taps; j++) {
					int arg = max(0,min(i+skip*(offset-j),srcwidth-1));//clamped BC's
					lo += filterLo[j] * srcbuffer[arg];//lopass channel
					hi += filterHi[j] * srcbuffer[arg];//hipass channel
				}
			}
			dstLo[row*dstwidth+((i/2))] = lo;
			dstHi[row*dstwidth+((i/2))] = hi;
		}
	}

#ifdef __SSE2__
	template<typename T> SSEFUNCTION void wavelet_level<T>::AnalysisFilterSubsampVertical (T * RESTRICT srcbuffer, T * RESTRICT dstLo, T * RESTRICT dstHi, float (* RESTRICT filterLo)[4], float (* RESTRICT filterHi)[4],
													const int taps, const int offset, const int width, const int height, const int row) {

		/* Basic convolution code
		 * Applies an FIR filter 'filter' with filter length 'taps', 
		 * aligning the 'offset' element of the filter with
		 * the input pixel, and skipping 'skip' pixels between taps 
		 * Output is subsampled by two
		 */

		// calculate coefficients
		if (LIKELY(row>skip*taps && row<height-skip*taps)) {//bulk
			int k;
			for (k=0; k<width-3; k+=4) {
				__m128 lov = _mm_setzero_ps();
				__m128 hiv = _mm_setzero_ps();
				for (int j=0, l=-skip*offset; j<taps; j++, l+=skip) {
					__m128 srcv = LVFU(srcbuffer[(row-l)*width+k]);
					lov += LVF(filterLo[j][0]) * srcv;//lopass channel
					hiv += LVF(filterHi[j][0]) * srcv;//hipass channel
				}
				STVF(dstLo[k], lov);
				STVF(dstHi[k], hiv);
			}
			for (; k<width; k++) {
				float lo = 0.f, hi = 0.f;
				for (int j=0, l=-skip*offset; j<taps; j++, l+=skip) {
					lo += filterLo[j][0] * srcbuffer[(row-l)*width+k];//lopass channel
					hi += filterHi[j][0] * srcbuffer[(row-l)*width+k];//hipass channel
				}
				dstLo[k] = lo;
				dstHi[k] = hi;
			}
		} else {//boundary
			int k;
			for (k=0; k<width-3; k+=4) {
				__m128 lov = _mm_setzero_ps();
				__m128 hiv = _mm_setzero_ps();
				for (int j=0; j<taps; j++) {
					int arg = max(0,min(row+skip*(offset-j),height-1))*width+k;//clamped BC's
					__m128 srcv = LVFU(srcbuffer[arg]);
					lov += LVF(filterLo[j][0]) * srcv;//lopass channel
					hiv += LVF(filterHi[j][0]) * srcv;//hipass channel
				}
				STVF(dstLo[k], lov);
				STVF(dstHi[k], hiv);
			}
			for (; k<width; k++) {
				float lo = 0.f, hi = 0.f;
				for (int j=0; j<taps; j++) {
					int arg = max(0,min(row+skip*(offset-j),height-1))*width+k;//clamped BC's
					lo += filterLo[j][0] * srcbuffer[arg];//lopass channel
					hi += filterHi[j][0] * srcbuffer[arg];//hipass channel
				}
				dstLo[k] = lo;
				dstHi[k] = hi;
			}
		}
	}
#else
	template<typename T> void wavelet_level<T>::AnalysisFilterSubsampVertical (T * RESTRICT srcbuffer, T * RESTRICT dstLo, T * RESTRICT dstHi, float * RESTRICT filterLo, float * RESTRICT filterHi, 
													const int taps, const int offset, const int width, const int height, const int row) {

		/* Basic convolution code
		 * Applies an FIR filter 'filter' with filter length 'taps', 
		 * aligning the 'offset' element of the filter with
		 * the input pixel, and skipping 'skip' pixels between taps 
		 * Output is subsampled by two
		 */

		// calculate coefficients
		if (LIKELY(row>skip*taps && row<height-skip*taps)) {//bulk
			for (int k=0; k<width; k++) {
				float lo = 0.f, hi = 0.f;
				for (int j=0, l=-skip*offset; j<taps; j++, l+=skip) {
					lo += filterLo[j] * srcbuffer[(row-l)*width+k];//lopass channel
					hi += filterHi[j] * srcbuffer[(row-l)*width+k];//hipass channel
				}
				dstLo[k] = lo;
				dstHi[k] = hi;
			}
		} else {//boundary
			for (int k=0; k<width; k++) {
				float lo = 0.f, hi = 0.f;
				for (int j=0; j<taps; j++) {
					int arg = max(0,min(row+skip*(offset-j),height-1))*width+k;//clamped BC's
					lo += filterLo[j] * srcbuffer[arg];//lopass channel
					hi += filterHi[j] * srcbuffer[arg];//hipass channel
				}
				dstLo[k] = lo;
				dstHi[k] = hi;
			}
		}
	}
#endif


	template<typename T> void wavelet_level<T>::SynthesisFilterSubsampHorizontal (T * RESTRICT srcLo, T * RESTRICT srcHi, T * RESTRICT dst, float * RESTRICT filterLo, float * RESTRICT filterHi, const int taps, const int offset, const int srcwidth, const int dstwidth, const int height) {

		/* Basic convolution code
		 * Applies an FIR filter 'filter' with filter length 'taps', 
		 * aligning the 'offset' element of the filter with
		 * the input pixel, and skipping 'skip' pixels between taps 
		 * Output is subsampled by two
		 */

		// calculate coefficients
		int shift = skip*(taps-offset-1);//align filter with data
#ifdef _OPENMP
#pragma omp parallel for num_threads(numThreads) if(numThreads>1)
#endif
		for (int k=0; k<height; k++) {
			int i;	
			for(i=0; i<=min(skip*taps,dstwidth); i++) {
				float tot=0.f;
				//TODO: this is correct only if skip=1; otherwise, want to work with cosets of length 'skip'
				int i_src = (i+shift)/2;
				int begin = (i+shift)%2;
				for (int j=begin, l=0; j<taps; j+=2, l+=skip) {
					int arg = max(0,min((i_src-l),srcwidth-1));//clamped BC's
					tot += ((filterLo[j] * srcLo[k*srcwidth+arg] + filterHi[j] * srcHi[k*srcwidth+arg]));
				}
				dst[k*dstwidth+i] = tot;
			}
			for(; i<min(dstwidth-skip*taps,dstwidth); i++) {
				float tot=0.f;
				//TODO: this is correct only if skip=1; otherwise, want to work with cosets of length 'skip'
				int i_src = (i+shift)/2;
				int begin = (i+shift)%2;
				for (int j=begin, l=0; j<taps; j+=2, l+=skip) {
					tot += ((filterLo[j] * srcLo[k*srcwidth+i_src-l] + filterHi[j] * srcHi[k*srcwidth+i_src-l]));
				}
				dst[k*dstwidth+i] = tot;
			}
			for(; i < dstwidth; i++) {
				float tot=0.f;
				//TODO: this is correct only if skip=1; otherwise, want to work with cosets of length 'skip'
				int i_src = (i+shift)/2;
				int begin = (i+shift)%2;
				for (int j=begin, l=0; j<taps; j+=2, l+=skip) {
					int arg = max(0,min((i_src-l),srcwidth-1));//clamped BC's
					tot += ((filterLo[j] * srcLo[k*srcwidth+arg] + filterHi[j] * srcHi[k*srcwidth+arg]));
				}
				dst[k*dstwidth+i] = tot;
			}
		}
	}

#ifdef __SSE2__
	template<typename T> SSEFUNCTION void wavelet_level<T>::SynthesisFilterSubsampVertical (T * RESTRICT srcLo, T * RESTRICT srcHi, T * RESTRICT dst, float (* RESTRICT filterLo)[4], float (* RESTRICT filterHi)[4], const int taps, const int offset, const int width, const int srcheight, const int dstheight, const float blend)
	 {

		/* Basic convolution code
		 * Applies an FIR filter 'filter' with filter length 'taps', 
		 * aligning the 'offset' element of the filter with
		 * the input pixel, and skipping 'skip' pixels between taps 
		 * Output is subsampled by two
		 */
		const float srcFactor = 1.f - blend;
		// calculate coefficients
		int shift=skip*(taps-offset-1);//align filter with data
		__m128 fourv = _mm_set1_ps(4.f);
		__m128 srcFactorv = _mm_set1_ps(srcFactor);
		__m128 dstFactorv = _mm_set1_ps(blend);
#ifdef _OPENMP
#pragma omp parallel for num_threads(numThreads) if(numThreads>1)
#endif
		for(int i = 0; i < dstheight; i++) {
			int i_src = (i+shift)/2;
			int begin = (i+shift)%2;
			//TODO: this is correct only if skip=1; otherwise, want to work with cosets of length 'skip'
			if (LIKELY(i>skip*taps && i<(dstheight-skip*taps))) {//bulk
				int k;
				for (k=0; k<width-3; k+=4) {
					__m128 totv = _mm_setzero_ps();
					for (int j=begin, l=0; j<taps; j+=2, l+=skip) {
						totv += ((LVF(filterLo[j][0]) * LVFU(srcLo[(i_src-l)*width+k]) + LVF(filterHi[j][0]) * LVFU(srcHi[(i_src-l)*width+k])));
					}
					_mm_storeu_ps(&dst[width*i+k], LVFU(dst[width*i+k]) * srcFactorv + dstFactorv * fourv * totv);
				}
				for (; k<width; k++) {
					float tot = 0.f;
					for (int j=begin, l=0; j<taps; j+=2, l+=skip) {
						tot += ((filterLo[j][0] * srcLo[(i_src-l)*width+k] + filterHi[j][0] * srcHi[(i_src-l)*width+k]));
					}
					dst[width*i+k] = dst[width*i+k] * srcFactor + blend * 4.f * tot;
				}
			} else {//boundary
				int k;
				for (k=0; k<width-3; k+=4) {
					__m128 totv = _mm_setzero_ps();
					for (int j=begin, l=0; j<taps; j+=2, l+=skip) {
						int arg = max(0,min((i_src-l),srcheight-1))*width+k;//clamped BC's
						totv += ((LVF(filterLo[j][0]) * LVFU(srcLo[arg]) + LVF(filterHi[j][0]) * LVFU(srcHi[arg])));
					}
					_mm_storeu_ps(&dst[width*i+k], LVFU(dst[width*i+k]) * srcFactorv + dstFactorv * fourv * totv);
				}
				for (; k<width; k++) {
					float tot = 0.f;
					for (int j=begin, l=0; j<taps; j+=2, l+=skip) {
						int arg = max(0,min((i_src-l),srcheight-1))*width+k;//clamped BC's
						tot += ((filterLo[j][0] * srcLo[arg] + filterHi[j][0] * srcHi[arg]));
					}
					dst[width*i+k] = dst[width*i+k] * srcFactor + blend * 4.f * tot;
				}
			}
		}
	}
#else
	template<typename T> void wavelet_level<T>::SynthesisFilterSubsampVertical (T * RESTRICT srcLo, T * RESTRICT srcHi, T * RESTRICT dst, float * RESTRICT filterLo, float * RESTRICT filterHi, const int taps, const int offset, const int width, const int srcheight, const int dstheight, const float blend)
	 {

		/* Basic convolution code
		 * Applies an FIR filter 'filter' with filter length 'taps', 
		 * aligning the 'offset' element of the filter with
		 * the input pixel, and skipping 'skip' pixels between taps 
		 * Output is subsampled by two
		 */

		const float srcFactor = 1.f - blend;
		// calculate coefficients
		int shift=skip*(taps-offset-1);//align filter with data

#ifdef _OPENMP
#pragma omp parallel for num_threads(numThreads) if(numThreads>1)
#endif
		for(int i = 0; i < dstheight; i++) {
			int i_src = (i+shift)/2;
			int begin = (i+shift)%2;
			//TODO: this is correct only if skip=1; otherwise, want to work with cosets of length 'skip'
			if (LIKELY(i>skip*taps && i<(dstheight-skip*taps))) {//bulk
				for (int k=0; k<width; k++) {
					float tot = 0.f;
					for (int j=begin, l=0; j<taps; j+=2, l+=skip) {
						tot += ((filterLo[j] * srcLo[(i_src-l)*width+k] + filterHi[j] * srcHi[(i_src-l)*width+k]));
					}
					dst[width*i+k] = dst[width*i+k] * srcFactor + blend * 4.f * tot;
				}
			} else {//boundary
				for (int k=0; k<width; k++) {
					float tot = 0.f;
					for (int j=begin, l=0; j<taps; j+=2, l+=skip) {
						int arg = max(0,min((i_src-l),srcheight-1))*width+k;//clamped BC's
						tot += ((filterLo[j] * srcLo[arg] + filterHi[j] * srcHi[arg]));
					}
					dst[width*i+k] = dst[width*i+k] * srcFactor + blend * 4.f * tot;
				}
			}
		}
	}
#endif

#ifdef __SSE2__
	template<typename T> template<typename E> SSEFUNCTION void wavelet_level<T>::decompose_level(E *src, E *dst, float *filterV, float *filterH, int taps, int offset) { 

		/* filter along rows and columns */
		float filterVarray[2*taps][4] ALIGNED64;
		if(subsamp_out) {
			for(int i=0;i<2*taps;i++) {
				for(int j=0;j<4;j++) {
					filterVarray[i][j] = filterV[i];
				}
			}
		}
#ifdef _OPENMP
#pragma omp parallel num_threads(numThreads) if(numThreads>1)
#endif
{
		T tmpLo[m_w] ALIGNED64;
		T tmpHi[m_w] ALIGNED64;
		if(subsamp_out) {
#ifdef _OPENMP
#pragma omp for
#endif
			for(int row=0;row<m_h;row+=2) {
				AnalysisFilterSubsampVertical (src, tmpLo, tmpHi, filterVarray, filterVarray+taps, taps, offset, m_w, m_h, row);
				AnalysisFilterSubsampHorizontal (tmpLo, dst, wavcoeffs[1], filterH, filterH+taps, taps, offset, m_w, m_w2, row/2);
				AnalysisFilterSubsampHorizontal (tmpHi, wavcoeffs[2], wavcoeffs[3], filterH, filterH+taps, taps, offset, m_w, m_w2, row/2);
			}
		} else {
#ifdef _OPENMP
#pragma omp for
#endif
			for(int row=0;row<m_h;row++) {
				AnalysisFilterHaarVertical (src, tmpLo, tmpHi, m_w, m_h, row);
				AnalysisFilterHaarHorizontal (tmpLo, dst, wavcoeffs[1], m_w, row);
				AnalysisFilterHaarHorizontal (tmpHi, wavcoeffs[2], wavcoeffs[3], m_w, row);
			}
		}
}
	}
#else
	template<typename T> template<typename E> void wavelet_level<T>::decompose_level(E *src, E *dst, float *filterV, float *filterH, int taps, int offset) { 

#ifdef _OPENMP
#pragma omp parallel num_threads(numThreads) if(numThreads>1)
#endif
{
		T tmpLo[m_w] ALIGNED64;
		T tmpHi[m_w] ALIGNED64;
		/* filter along rows and columns */
		if(subsamp_out) {
#ifdef _OPENMP
#pragma omp for
#endif
			for(int row=0;row<m_h;row+=2) {
				AnalysisFilterSubsampVertical (src, tmpLo, tmpHi, filterV, filterV+taps, taps, offset, m_w, m_h, row);
				AnalysisFilterSubsampHorizontal (tmpLo, dst, wavcoeffs[1], filterH, filterH+taps, taps, offset, m_w, m_w2, row/2);
				AnalysisFilterSubsampHorizontal (tmpHi, wavcoeffs[2], wavcoeffs[3], filterH, filterH+taps, taps, offset, m_w, m_w2, row/2);
			}
		} else {
#ifdef _OPENMP
#pragma omp for
#endif
			for(int row=0;row<m_h;row++) {
				AnalysisFilterHaarVertical (src, tmpLo, tmpHi, m_w, m_h, row);
				AnalysisFilterHaarHorizontal (tmpLo, dst, wavcoeffs[1], m_w, row);
				AnalysisFilterHaarHorizontal (tmpHi, wavcoeffs[2], wavcoeffs[3], m_w, row);
			}
		}
}
	}
#endif

#ifdef __SSE2__

	template<typename T> template<typename E> SSEFUNCTION void wavelet_level<T>::reconstruct_level(E* tmpLo, E* tmpHi, E * src, E *dst, float *filterV, float *filterH, int taps, int offset, const float blend) { 
		if(memoryAllocationFailed)
			return;

		/* filter along rows and columns */
		if (subsamp_out) {
			float filterVarray[2*taps][4] ALIGNED64;
			for(int i=0;i<2*taps;i++) {
				for(int j=0;j<4;j++) {
					filterVarray[i][j] = filterV[i];
				}
			}
			SynthesisFilterSubsampHorizontal (wavcoeffs[2], wavcoeffs[3], tmpHi, filterH, filterH+taps, taps, offset, m_w2, m_w, m_h2);
			SynthesisFilterSubsampHorizontal (src, wavcoeffs[1], tmpLo, filterH, filterH+taps, taps, offset, m_w2, m_w, m_h2);
			SynthesisFilterSubsampVertical (tmpLo, tmpHi, dst, filterVarray, filterVarray+taps, taps, offset, m_w, m_h2, m_h, blend);
		} else {
			SynthesisFilterHaarHorizontal (wavcoeffs[2], wavcoeffs[3], tmpHi, m_w, m_h2);
			SynthesisFilterHaarHorizontal (src, wavcoeffs[1], tmpLo, m_w, m_h2);
			SynthesisFilterHaarVertical (tmpLo, tmpHi, dst, m_w, m_h);
		}
	}
#else
	template<typename T> template<typename E> void wavelet_level<T>::reconstruct_level(E* tmpLo, E* tmpHi, E * src, E *dst, float *filterV, float *filterH, int taps, int offset, const float blend) { 
		if(memoryAllocationFailed)
			return;
		/* filter along rows and columns */
		if (subsamp_out) {
			SynthesisFilterSubsampHorizontal (wavcoeffs[2], wavcoeffs[3], tmpHi, filterH, filterH+taps, taps, offset, m_w2, m_w, m_h2);
			SynthesisFilterSubsampHorizontal (src, wavcoeffs[1], tmpLo, filterH, filterH+taps, taps, offset, m_w2, m_w, m_h2);
			SynthesisFilterSubsampVertical (tmpLo, tmpHi, dst, filterV, filterV+taps, taps, offset, m_w, m_h2, m_h, blend);
		} else {
			SynthesisFilterHaarHorizontal (wavcoeffs[2], wavcoeffs[3], tmpHi, m_w, m_h2);
			SynthesisFilterHaarHorizontal (src, wavcoeffs[1], tmpLo, m_w, m_h2);
			SynthesisFilterHaarVertical (tmpLo, tmpHi, dst, m_w, m_h);
		}
	}
#endif
};

#endif
