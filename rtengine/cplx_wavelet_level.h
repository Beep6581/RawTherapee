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
 */

#ifndef CPLX_WAVELET_LEVEL_H_INCLUDED
#define CPLX_WAVELET_LEVEL_H_INCLUDED

#include <cstddef>
#include <algorithm>

#include "array2D.h"

namespace rtengine {
	
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) > (b) ? (b) : (a))
	
	
	//////////////////////////////////////////////////////////////////////////////
	
	template<typename T>
	class cplx_wavelet_level
	{
		// full size
		size_t m_w, m_h;
		
		// size of low frequency part
		size_t m_w2, m_h2;
		
		// size of padded border
		size_t m_pad;
		
		// array of pointers to lines of coeffs
		// actually is a single contiguous data array pointed by m_coeffs[0]
		//T ** m_coeffs;
		//array2D<float> wavcoeffs(4,1);
		//data structure: first label is output channel (LL,LH,HL,HH), second is pixel location in flattened array
		
		// weights storage
		//T ** m_weights_rows;
		//T ** m_weights_cols;
		
		// allocation and destruction of data storage
		T ** create(size_t n);
		void destroy(T ** subbands);
		
		// load a row/column of input data, possibly with padding
		template<typename E>
		void loadbuffer(E * src, E * dst, int srclen, int pitch);
		
		//void dwt_2d(size_t w, size_t h);
		//void idwt_2d(size_t w, size_t h, int alpha);
		
		void AnalysisFilter (T * src, T * dstLo, T * dstHi, float *filterLo, float *filterHi, 
							 int taps, int offset, int pitch, int srclen);
		void SynthesisFilter (T * srcLo, T * srcHi, T * dst, T *bufferLo, T *bufferHi, 
							  float *filterLo, float *filterHi, int taps, int offset, int pitch, int dstlen);
		
	public:
		
		T ** wavcoeffs;
		
		template<typename E>
		cplx_wavelet_level(E * src, int padding, size_t w, size_t h, float *filterV, float *filterH, int len, int offset)
		: m_w(w), m_h(h), m_w2((w+1+2*padding)/2), m_h2((h+1+2*padding)/2), m_pad(padding), wavcoeffs(NULL)
		{
			
			//m_coeffs = create(w, h);
			//m_weights_rows = create(w + 4, h);
			//m_weights_cols = create(h + 4, w);
			
			//decompose_level(src, w, h, wavcoeffs, float **filterV, float **filterH, int len, int offset);
			
			wavcoeffs = create((m_w2)*(m_h2));
			decompose_level(src, filterV, filterH, len, offset);
			
		}
		
		~cplx_wavelet_level()
		{
			destroy(wavcoeffs);
		}
		
		T ** subbands() const
		{
			return wavcoeffs;
		}
		
		T * lopass() const
		{
			return wavcoeffs[0];
		}
		
		size_t width() const
		{
			return m_w2;
		}
		
		size_t height() const
		{
			return m_h2;
		}
		
		template<typename E>
		void decompose_level(E *src, float *filterV, float *filterH, int len, int offset);
		
		template<typename E>
		void reconstruct_level(E *dst, float *filterV, float *filterH, int len, int offset);
		
	};
	
	//////////////////////////////////////////////////////////////////////////////
	
	
	
	template<typename T>
	T ** cplx_wavelet_level<T>::create(size_t n)
	{
		T * data = new T[4*n];
		T ** subbands = new T*[4];
		for(size_t j = 0; j < 4; j++)
		{
			subbands[j] = data + n * j;
		}
		return subbands;
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	
	
	template<typename T>
	void cplx_wavelet_level<T>::destroy(T ** subbands)
	{
		if(subbands)
		{
			delete[] subbands[0];
			delete[] subbands;
		}
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	
	template<typename T> template<typename E>
	void cplx_wavelet_level<T>::loadbuffer(E * src, E * dst, int pitch, int srclen)
	{
		E * tmp = dst + m_pad;
		memset(dst, 0, (srclen+2*m_pad)*sizeof(E));

		for(size_t i = 0, j = 0; i<srclen; i++, j += pitch)
		{
			tmp[i] = src[j];
		}	
		
		// extend mirror-like
		
		for (size_t i=1; i<=MIN(srclen-1,m_pad); i++) {
			tmp[-i] = tmp[i];
			tmp[srclen+i-1] = tmp[srclen-i-1];
		}
	}
	
	
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	template<typename T>
	void cplx_wavelet_level<T>::AnalysisFilter (T * src, T * dstLo, T * dstHi, float *filterLo, float *filterHi, 
												int taps, int offset, int pitch, int srclen) {
		
		/* Basic convolution code
		 * Applies an FIR filter 'filter' with filter length 'taps', 
		 * aligning the 'offset' element of the filter
		 * with the input pixel, and skipping 'pitch' pixels
		 * between taps (eg pitch=1 for horizontal filtering, 
		 * pitch=W for vertical, pitch=W+1,W-1 for diagonals.
		 * Currently diagonal filtering is not supported
		 * for the full source array, until a more sophisticated 
		 * treatment of mirror BC's is implemented.
		 *
		 */
				
		// calculate coefficients
		
		for(int i = 0; i < (srclen); i+=2) {
			float lo=0,hi=0;
			if (i>taps && i<srclen-taps) {//bulk
				for (int j=0; j<taps; j++) {
					lo += filterLo[j] * src[(i+offset-j)];//lopass channel
					hi += filterHi[j] * src[(i+offset-j)];//hipass channel
				}
			} else {//boundary
				for (int j=0; j<taps; j++) {
					int arg = MAX(0,MIN(i+offset-j,srclen-1));//clamped BC's
					lo += filterLo[j] * src[arg];//lopass channel
					hi += filterHi[j] * src[arg];//hipass channel
				}
			}

			dstLo[(pitch*(i/2))] = lo;
			dstHi[(pitch*(i/2))] = hi;
		}
		
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	template<typename T>
	void cplx_wavelet_level<T>::SynthesisFilter (T * srcLo, T * srcHi, T * dst, T *bufferLo, T *bufferHi, 
												 float *filterLo, float *filterHi, int taps, int offset, int pitch, int dstlen) {
		
		/* Basic convolution code
		 * Applies an FIR filter 'filter' with 'len' taps, 
		 * aligning the 'offset' element of the filter
		 * with the input pixel, and skipping 'pitch' pixels
		 * between taps (eg pitch=1 for horizontal filtering, 
		 * pitch=W for vertical, pitch=W+1,W-1 for diagonals.
		 * Currently diagonal filtering is not supported
		 * for the full source array, until a more sophisticated 
		 * treatment of mirror BC's is implemented.
		 *
		 */

		
		
		// calculate coefficients
		
		int srclen=(dstlen+1+2*m_pad)/2;
		for (int i=0; i<srclen; i++) {
			bufferLo[i]=srcLo[i*pitch];
			bufferHi[i]=srcHi[i*pitch];
		}
		
		int shift=taps-offset-1;
		for(int i = m_pad; i < (dstlen-m_pad); i++) {
			if (bufferLo[i]!=0) {
				float xxx=bufferLo[i];
			}
			float tot=0;
			int i_src = (i+shift)/2;
			int begin = (i+shift)%2;
			if (i>taps && i<(srclen-taps)) {//bulk
				for (int j=begin, l=0; j<taps; j+=2, l++) {
					tot += (filterLo[j] * bufferLo[i_src-l] + filterHi[j] * bufferHi[i_src-l]);
				}
			} else {//boundary
				for (int j=begin, l=0; j<taps; j+=2, l++) {
					int arg = MAX(0,MIN((i_src-l),srclen-1));//clamped BC's
					tot += (filterLo[j] * bufferLo[arg] + filterHi[j] * bufferHi[arg]);
				}
			}

			dst[pitch*(i-m_pad)] = tot;
		}
		
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	template<typename T> template<typename E>
	void cplx_wavelet_level<T>::decompose_level(E *src, float *filterV, float *filterH, int taps, int offset) { 
		
		T *tmpLo = new T[m_w*m_h2];
		T *tmpHi = new T[m_w*m_h2];
		
		T *buffer = new T[MAX(m_w,m_h)+2*m_pad];
		
		/* filter along columns */
		for (int j=0; j<m_w; j++) {
			loadbuffer(src+j, buffer, m_w/*pitch*/, m_h/*srclen*/);//pad a column of data and load it to buffer
			AnalysisFilter (buffer, tmpLo+j, tmpHi+j, filterV, filterV+taps, taps, offset, m_w/*output_pitch*/, m_h+2*m_pad/*srclen*/);
		}
		
		/* filter along rows */
		for (int i=0; i<m_h2; i++) {
			loadbuffer(tmpLo+i*m_w, buffer, 1/*pitch*/, m_w/*srclen*/);//pad a row of data and load it to buffer
			AnalysisFilter (buffer, wavcoeffs[0]+i*m_w2, wavcoeffs[1]+i*m_w2, filterH, filterH+taps, taps, offset, 1/*output_pitch*/, m_w+2*m_pad/*srclen*/);
			loadbuffer(tmpHi+i*m_w, buffer, 1/*pitch*/, m_w/*srclen*/);
			AnalysisFilter (buffer, wavcoeffs[2]+i*m_w2, wavcoeffs[3]+i*m_w2, filterH, filterH+taps, taps, offset, 1/*output_pitch*/, m_w+2*m_pad/*srclen*/);
		}

		delete[] tmpLo;
		delete[] tmpHi;
		delete[] buffer;
	}
	
	/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
	
	template<typename T> template<typename E>
	void cplx_wavelet_level<T>::reconstruct_level(E *dst, float *filterV, float *filterH, int taps, int offset) { 
		
		
		//int hfw = (W+1)/2;
		//int hfh = (H+1)/2;
		T *tmpLo = new T[m_w*m_h2];
		T *tmpHi = new T[m_w*m_h2];
		
		int buflen = MAX(m_w,m_h);
		float *bufferLo = new float[buflen];
		float *bufferHi = new float[buflen];
		
		/* filter along rows */
		for (int i=0; i<m_h2; i++) {
			
			SynthesisFilter (wavcoeffs[0]+i*m_w2, wavcoeffs[1]+i*m_w2, tmpLo+i*m_w, bufferLo, bufferHi, \ 
							 filterH, filterH+taps, taps, offset, 1/*pitch*/, m_w/*dstlen*/);
			SynthesisFilter (wavcoeffs[2]+i*m_w2, wavcoeffs[3]+i*m_w2, tmpHi+i*m_w, bufferLo, bufferHi, \ 
							 filterH, filterH+taps, taps, offset, 1/*pitch*/, m_w/*dstlen*/);
		}
		
		/* filter along columns */
		for (int j=0; j<m_w; j++) {
			SynthesisFilter (tmpLo+j, tmpHi+j, dst+j, bufferLo, bufferHi, 
							 filterV, filterV+taps, taps, offset, m_w/*pitch*/, m_h/*dstlen*/);
		}

		
		delete[] tmpLo;
		delete[] tmpHi;
		delete[] bufferLo;
		delete[] bufferHi;
		
	}
	
	/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
	/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
	
	
};

#endif
