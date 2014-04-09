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

#include "gauss.h"
#include "rt_math.h"

namespace rtengine {

	
	//////////////////////////////////////////////////////////////////////////////
	
	template<typename T>
	class wavelet_level
	{
		// full size
		size_t m_w, m_h;
		
		// size of low frequency part
		size_t m_w2, m_h2;
		
		// size of padded border
		size_t m_pad;
		
		// level of decomposition
		int lvl;
		
		// whether to subsample the output
		bool subsamp_out;
		
		// spacing of filter taps
		size_t skip;
		
		// allocation and destruction of data storage
		T ** create(size_t n);
		void destroy(T ** subbands);
		
		// load a row/column of input data, possibly with padding
		template<typename E>
		void loadbuffer(E * src, E * dst, int srclen, int pitch);
		
		void AnalysisFilter (T * srcbuffer, T * dstLo, T * dstHi, float *filterLo, float *filterHi, 
							 int taps, int offset, int pitch, int srclen);
		void SynthesisFilter (T * srcLo, T * srcHi, T * dst, T *bufferLo, T *bufferHi, 
							  float *filterLo, float *filterHi, int taps, int offset, int pitch, int dstlen);
		
		void AnalysisFilterHaar (T * srcbuffer, T * dstLo, T * dstHi, int pitch, int srclen);
		void SynthesisFilterHaar (T * srcLo, T * srcHi, T * dst, T *bufferLo, T *bufferHi, int pitch, int dstlen);
		
		void AnalysisFilterSubsamp (T * srcbuffer, T * dstLo, T * dstHi, float *filterLo, float *filterHi, 
							 int taps, int offset, int pitch, int srclen);
		void SynthesisFilterSubsamp (T * srcLo, T * srcHi, T * dst, T *bufferLo, T *bufferHi, 
							  float *filterLo, float *filterHi, int taps, int offset, int pitch, int dstlen);
		
		void AnalysisFilterSubsampHaar (T * srcbuffer, T * dstLo, T * dstHi, int pitch, int srclen);
		void SynthesisFilterSubsampHaar (T * srcLo, T * srcHi, T * dst, int pitch, int dstlen);
		
		void imp_nr (T* src, int width, int height, double thresh);

		
	public:
		
		T ** wavcoeffs;
		
		template<typename E>
		wavelet_level(E * src, int level, int subsamp, int padding, size_t w, size_t h, float *filterV, float *filterH, int len, int offset)
		: m_w(w), m_h(h), m_w2(w), m_h2(h), m_pad(padding), wavcoeffs(NULL), lvl(level), skip(1<<level), subsamp_out((subsamp>>level)&1)
		{
			if (subsamp) {
				skip = 1;
				for (int n=0; n<level; n++) {
					skip *= 2-((subsamp>>n)&1);
				}
			}
			m_w2 = (subsamp_out ? ((w+1+2*skip*padding)/2) : (w+2*skip*padding));
			m_h2 = (subsamp_out ? ((h+1+2*skip*padding)/2) : (h+2*skip*padding));
			m_pad= skip*padding;
			
			wavcoeffs = create((m_w2)*(m_h2));
			decompose_level(src, filterV, filterH, len, offset);
			
		}
		
		~wavelet_level()
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
		
		size_t padding() const
		{
			return m_pad/skip;
		}
		
		size_t stride() const
		{
			return skip;
		}
		
		template<typename E>
		void decompose_level(E *src, float *filterV, float *filterH, int len, int offset);
		
		template<typename E>
		void reconstruct_level(E *dst, float *filterV, float *filterH, int len, int offset);
		
	};
	
	//////////////////////////////////////////////////////////////////////////////
	
	
	
	template<typename T>
	T ** wavelet_level<T>::create(size_t n)
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
	void wavelet_level<T>::destroy(T ** subbands)
	{
		if(subbands)
		{
			delete[] subbands[0];
			delete[] subbands;
		}
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

	template<typename T> template<typename E>
	void wavelet_level<T>::loadbuffer(E * src, E * dst, int pitch, int srclen)
	{
		E * tmp = dst + m_pad;
		memset(dst, 0, (MAX(m_w2,m_h2))*sizeof(E));
		
		//create padded buffer from src data
		for (size_t i = 0, j = 0; i<srclen; i++, j += pitch)
		{
			tmp[i] = src[j];
		}	
		
		// extend each coset mirror-like by padding amount 'm_pad' and to a multiple of 'skip'
		for (size_t i=1; i<=MIN(srclen-1,m_pad); i++) {
			tmp[-i] = tmp[i];
			tmp[srclen+i-1] = tmp[srclen-i-1];
		}
		for (size_t i=0; i<srclen%skip; i++) {
			tmp[srclen+m_pad+i] = tmp[srclen+m_pad-i-2];
		}
		
	}
	
	
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	template<typename T>
	void wavelet_level<T>::AnalysisFilter (T * srcbuffer, T * dstLo, T * dstHi, float *filterLo, float *filterHi, 
												int taps, int offset, int pitch, int srclen) {
		
		/* Basic convolution code
		 * Applies an FIR filter 'filter' with filter length 'taps', 
		 * aligning the 'offset' element of the filter with
		 * the input pixel, and skipping 'skip' pixels between taps 
		 *
		 */
		
				
		for (size_t i = 0; i < (srclen); i++) {
			float lo=0,hi=0;
			if (i>skip*taps && i<srclen-skip*taps) {//bulk
				for (int j=0, l=-skip*offset; j<taps; j++, l+=skip) {
					lo += filterLo[j] * srcbuffer[i-l];//lopass channel
					hi += filterHi[j] * srcbuffer[i-l];//hipass channel
				}
			} else {//boundary
				for (int j=0; j<taps; j++) {
					int arg = MAX(0,MIN(i+skip*(offset-j),srclen-1));//clamped BC's
					lo += filterLo[j] * srcbuffer[arg];//lopass channel
					hi += filterHi[j] * srcbuffer[arg];//hipass channel
				}
			}
			
			dstLo[(pitch*(i))] = lo;
			dstHi[(pitch*(i))] = hi;
		}
		
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

	template<typename T>
	void wavelet_level<T>::SynthesisFilter (T * srcLo, T * srcHi, T * dst, T *bufferLo, T *bufferHi, float *filterLo, 
												 float *filterHi, int taps, int offset, int pitch, int dstlen) {
		
		/* Basic convolution code
		 * Applies an FIR filter 'filter' with filter length 'taps', 
		 * aligning the 'offset' element of the filter with
		 * the input pixel, and skipping 'skip' pixels between taps 
		 *
		 */
		
		
		// load into buffer
		
		int srclen = (dstlen==m_w ? m_w2 : m_h2);//length of row/col in src (coarser level)
		
		for (size_t i=0, j=0; i<srclen; i++, j+=pitch) {
			bufferLo[i]=srcLo[j];
			bufferHi[i]=srcHi[j];
		}
		
		int shift=skip*(taps-offset-1);
		for(size_t i = m_pad; i < (dstlen+m_pad); i++) {
			float tot=0;
			if (i>skip*taps && i<(srclen-skip*taps)) {//bulk
				for (int j=0, l=-shift; j<taps; j++, l+=skip) {
					tot += (filterLo[j] * bufferLo[i-l] + filterHi[j] * bufferHi[i-l]);
				}
			} else {//boundary
				for (int j=0, l=-shift; j<taps; j++, l+=skip) {
					int arg = MAX(0,MIN((i-l),srclen-1));//clamped BC's
					tot += (filterLo[j] * bufferLo[arg] + filterHi[j] * bufferHi[arg]);
				}
			}
			
			dst[pitch*(i-m_pad)] = tot;

		}
		
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	template<typename T>
	void wavelet_level<T>::AnalysisFilterHaar (T * srcbuffer, T * dstLo, T * dstHi, int pitch, int srclen) {
		
		/* Basic convolution code
		 * Applies a Haar filter 
		 *
		 */										
		
		for(size_t i = 0; i < (srclen - skip); i++) {
			dstLo[(pitch*(i))] = 0.5*(srcbuffer[i] + srcbuffer[i+skip]);
			dstHi[(pitch*(i))] = 0.5*(srcbuffer[i] - srcbuffer[i+skip]);
		}
		// Start the loop at max(srclen-skip,skip) to avoid buffer underrun
		for(size_t i = max(srclen-skip,skip); i < (srclen); i++) {
			dstLo[(pitch*(i))] = 0.5*(srcbuffer[i] + srcbuffer[i-skip]);
			dstHi[(pitch*(i))] = 0.5*(srcbuffer[i] - srcbuffer[i-skip]);
		}
																		 
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	
	template<typename T>
	void wavelet_level<T>::SynthesisFilterHaar (T * srcLo, T * srcHi, T * dst, T *bufferLo, T *bufferHi, int pitch, int dstlen) {
		
		/* Basic convolution code
		 * Applies a Haar filter 
		 *
		 */
		
		int srclen = (dstlen==m_w ? m_w2 : m_h2);//length of row/col in src (coarser level)
		
		for (size_t i=0, j=0; i<srclen; i++, j+=pitch) {
			bufferLo[i]=srcLo[j];
			bufferHi[i]=srcHi[j];
		}
		
		for(size_t i = m_pad+skip; i < (dstlen+m_pad); i++) {
			dst[pitch*(i-m_pad)] = 0.5*(bufferLo[i] + bufferHi[i] + bufferLo[i-skip] - bufferHi[i-skip]);			
		}
		
		for(size_t i = (m_pad); i < (m_pad+skip); i++) {
			dst[pitch*(i-m_pad)] = (bufferLo[i] + bufferHi[i]);			
		}
	}

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


	template<typename T>
	void wavelet_level<T>::AnalysisFilterSubsamp (T * srcbuffer, T * dstLo, T * dstHi, float *filterLo, float *filterHi, 
												int taps, int offset, int pitch, int srclen) {
		
		/* Basic convolution code
		 * Applies an FIR filter 'filter' with filter length 'taps', 
		 * aligning the 'offset' element of the filter with
		 * the input pixel, and skipping 'skip' pixels between taps 
		 * Output is subsampled by two
		 */
		
		// calculate coefficients
		
		for(int i = 0; i < (srclen); i+=2) {
			float lo=0,hi=0;
			if (i>skip*taps && i<srclen-skip*taps) {//bulk
				for (int j=0, l=-skip*offset; j<taps; j++, l+=skip) {
					lo += filterLo[j] * srcbuffer[i-l];//lopass channel
					hi += filterHi[j] * srcbuffer[i-l];//hipass channel
				}
			} else {//boundary
				for (int j=0; j<taps; j++) {
					int arg = MAX(0,MIN(i+skip*(offset-j),srclen-1));//clamped BC's
					lo += filterLo[j] * srcbuffer[arg];//lopass channel
					hi += filterHi[j] * srcbuffer[arg];//hipass channel
				}
			}
			
			dstLo[(pitch*(i/2))] = lo;
			dstHi[(pitch*(i/2))] = hi;
		}
		
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	template<typename T>
	void wavelet_level<T>::SynthesisFilterSubsamp (T * srcLo, T * srcHi, T * dst, T *bufferLo, T *bufferHi, 
												 float *filterLo, float *filterHi, int taps, int offset, int pitch, int dstlen) {
		
		/* Basic convolution code
		 * Applies an FIR filter 'filter' with filter length 'taps', 
		 * aligning the 'offset' element of the filter with
		 * the input pixel, and skipping 'skip' pixels between taps 
		 * Output is subsampled by two
		 */
		
		
		
		// calculate coefficients
		
		int srclen = (dstlen==m_w ? m_w2 : m_h2);//length of row/col in src (coarser level)
		
		//fill a buffer with a given row/column of data
		for (size_t i=0, j=0; i<srclen; i++, j+=pitch) {
			bufferLo[i]=srcLo[j];
			bufferHi[i]=srcHi[j];
		}
				
		int shift=skip*(taps-offset-1);//align filter with data
		for(size_t i = m_pad; i < (dstlen+m_pad); i++) {
			
			float tot=0;
			//TODO: this is correct only if skip=1; otherwise, want to work with cosets of length 'skip'
			int i_src = (i+shift)/2;
			int begin = (i+shift)%2;
			if (i>skip*taps && i<(srclen-skip*taps)) {//bulk
				for (int j=begin, l=0; j<taps; j+=2, l+=skip) {
					tot += 2*((filterLo[j] * bufferLo[i_src-l] + filterHi[j] * bufferHi[i_src-l]));
				}
			} else {//boundary
				for (int j=begin, l=0; j<taps; j+=2, l+=skip) {
					int arg = MAX(0,MIN((i_src-l),srclen-1));//clamped BC's
					tot += 2*((filterLo[j] * bufferLo[arg] + filterHi[j] * bufferHi[arg]));
				}
			}
			
			dst[pitch*(i-m_pad)] = tot;
			
		}
		
	}
	
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	
	
	template<typename T>
	void wavelet_level<T>::AnalysisFilterSubsampHaar (T * srcbuffer, T * dstLo, T * dstHi, int pitch, int srclen) {
		
		/* Basic convolution code
		 * Applies a Haar filter
		 * Output is subsampled by two
		 */
		
		// calculate coefficients
		
		for(size_t i = 0; i < (srclen - skip); i+=2) {
			dstLo[(pitch*(i/2))] = 0.5*(srcbuffer[i] + srcbuffer[i+skip]);
			dstHi[(pitch*(i/2))] = 0.5*(srcbuffer[i] - srcbuffer[i+skip]);
		}
		
		for(size_t i = (srclen-skip)-((srclen-skip)&1); i < (srclen); i+=2) {
			dstLo[(pitch*(i/2))] = 0.5*(srcbuffer[i] + srcbuffer[i-skip]);
			dstHi[(pitch*(i/2))] = 0.5*(srcbuffer[i] - srcbuffer[i-skip]);
		}
		
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	template<typename T>
	void wavelet_level<T>::SynthesisFilterSubsampHaar (T * srcLo, T * srcHi, T * dst, int pitch, int dstlen) {
		
		/* Basic convolution code
		 * Applies a Haar filter 
		 * Input was subsampled by two
		 */
		
		
		// calculate coefficients
				
		//TODO: this code is buggy...
		for (int n=0; n<skip; n++) {
			for (size_t i = m_pad; i < (dstlen+m_pad-2*skip); i+=2*skip) {
				dst[pitch*(i-m_pad+n)] = srcLo[i/2+n]+srcHi[i/2+n];
				dst[pitch*(i-m_pad+skip+n)] = srcLo[i/2+n]-srcHi[i/2+n];
			}
		}
		
		if ((dstlen+m_pad-2*skip)<dstlen-1) {
			for (int n=0; n<skip; n++) {
				for (size_t i=(dstlen+m_pad-2*skip); i<dstlen+m_pad-1; i++) {
					dst[pitch*(dstlen-m_pad+n)] = srcLo[i/2+n]+srcHi[i/2+n];
				}
			}
		}
		
		
	}
	


	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	
	
	template<typename T> template<typename E>
	void wavelet_level<T>::decompose_level(E *src, float *filterV, float *filterH, int taps, int offset) { 
		
		T *tmpLo = new T[m_w*m_h2];
		T *tmpHi = new T[m_w*m_h2];
		
		T *buffer = new T[MAX(m_w,m_h)+2*m_pad+skip];
		
		/* filter along columns */
//OpenMP here		
		for (int j=0; j<m_w; j++) {
			loadbuffer(src+j, buffer, m_w/*pitch*/, m_h/*srclen*/);//pad a column of data and load it to buffer
			if (subsamp_out) {
				AnalysisFilterSubsamp (buffer, tmpLo+j, tmpHi+j, filterV, filterV+taps, taps, offset, m_w/*output_pitch*/, m_h/*srclen*/);
				//AnalysisFilterSubsampHaar (buffer, tmpLo+j, tmpHi+j, m_w, m_h);
			} else {
				//AnalysisFilter (buffer, tmpLo+j, tmpHi+j, filterV, filterV+taps, taps, offset, m_w/*output_pitch*/, m_h/*srclen*/);
				AnalysisFilterHaar (buffer, tmpLo+j, tmpHi+j, m_w, m_h);
			}
		}
		
		/* filter along rows */
//OpenMP here		
		for (int i=0; i<m_h2; i++) {
			loadbuffer(tmpLo+i*m_w, buffer, 1/*pitch*/, m_w/*srclen*/);//pad a row of data and load it to buffer
			if (subsamp_out) {
				AnalysisFilterSubsamp (buffer, wavcoeffs[0]+i*m_w2, wavcoeffs[1]+i*m_w2, filterH, filterH+taps, taps, offset, 1/*output_pitch*/, m_w/*srclen*/);
				//AnalysisFilterSubsampHaar (buffer, wavcoeffs[0]+i*m_w2, wavcoeffs[1]+i*m_w2, 1, m_w);
			} else {
				//AnalysisFilter (buffer, wavcoeffs[0]+i*m_w2, wavcoeffs[1]+i*m_w2, filterH, filterH+taps, taps, offset, 1/*output_pitch*/, m_w/*srclen*/);
				AnalysisFilterHaar (buffer, wavcoeffs[0]+i*m_w2, wavcoeffs[1]+i*m_w2, 1, m_w);
			}
			loadbuffer(tmpHi+i*m_w, buffer, 1/*pitch*/, m_w/*srclen*/);
			if (subsamp_out) {
				AnalysisFilterSubsamp (buffer, wavcoeffs[2]+i*m_w2, wavcoeffs[3]+i*m_w2, filterH, filterH+taps, taps, offset, 1/*output_pitch*/, m_w/*srclen*/);
				//AnalysisFilterSubsampHaar (buffer, wavcoeffs[2]+i*m_w2, wavcoeffs[3]+i*m_w2, 1, m_w);
			} else {
				//AnalysisFilter (buffer, wavcoeffs[2]+i*m_w2, wavcoeffs[3]+i*m_w2, filterH, filterH+taps, taps, offset, 1/*output_pitch*/, m_w/*srclen*/);
				AnalysisFilterHaar (buffer, wavcoeffs[2]+i*m_w2, wavcoeffs[3]+i*m_w2, 1, m_w);
			}
		}
				
		delete[] tmpLo;
		delete[] tmpHi;
		delete[] buffer;
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	
	template<typename T> template<typename E>
	void wavelet_level<T>::reconstruct_level(E *dst, float *filterV, float *filterH, int taps, int offset) { 
		
		T *tmpLo = new T[m_w*m_h2];
		T *tmpHi = new T[m_w*m_h2];
		
		int buflen = MAX(m_w2,m_h2);
		float *bufferLo = new float[buflen];
		float *bufferHi = new float[buflen];
		
		/* filter along rows */
//OpenMP here		
		for (int i=0; i<m_h2; i++) {
			
			if (subsamp_out) {
				SynthesisFilterSubsamp (wavcoeffs[0]+i*m_w2, wavcoeffs[1]+i*m_w2, tmpLo+i*m_w, bufferLo, bufferHi,
										filterH, filterH+taps, taps, offset, 1/*pitch*/, m_w/*dstlen*/);
				SynthesisFilterSubsamp (wavcoeffs[2]+i*m_w2, wavcoeffs[3]+i*m_w2, tmpHi+i*m_w, bufferLo, bufferHi,
										filterH, filterH+taps, taps, offset, 1/*pitch*/, m_w/*dstlen*/);
				//SynthesisFilterSubsampHaar (wavcoeffs[0]+i*m_w2, wavcoeffs[1]+i*m_w2, tmpLo+i*m_w, 1, m_w);//TODO: this is buggy
				//SynthesisFilterSubsampHaar (wavcoeffs[2]+i*m_w2, wavcoeffs[3]+i*m_w2, tmpHi+i*m_w, 1, m_w);
			} else {
				//SynthesisFilter (wavcoeffs[0]+i*m_w2, wavcoeffs[1]+i*m_w2, tmpLo+i*m_w, bufferLo, bufferHi,
				//				 filterH, filterH+taps, taps, offset, 1/*pitch*/, m_w/*dstlen*/);
				//SynthesisFilter (wavcoeffs[2]+i*m_w2, wavcoeffs[3]+i*m_w2, tmpHi+i*m_w, bufferLo, bufferHi,
				//				 filterH, filterH+taps, taps, offset, 1/*pitch*/, m_w/*dstlen*/);
				SynthesisFilterHaar (wavcoeffs[0]+i*m_w2, wavcoeffs[1]+i*m_w2, tmpLo+i*m_w, bufferLo, bufferHi, 1, m_w);
				SynthesisFilterHaar (wavcoeffs[2]+i*m_w2, wavcoeffs[3]+i*m_w2, tmpHi+i*m_w, bufferLo, bufferHi, 1, m_w);
			}
		}
		
		/* filter along columns */
//OpenMP here		
		for (int j=0; j<m_w; j++) {
			if (subsamp_out) {
				SynthesisFilterSubsamp (tmpLo+j, tmpHi+j, dst+j, bufferLo, bufferHi,
										filterV, filterV+taps, taps, offset, m_w/*pitch*/, m_h/*dstlen*/);
				//SynthesisFilterSubsampHaar (tmpLo+j, tmpHi+j, dst+j, m_w, m_h);
			} else {
				//SynthesisFilter (tmpLo+j, tmpHi+j, dst+j, bufferLo, bufferHi,
				//				 filterV, filterV+taps, taps, offset, m_w/*pitch*/, m_h/*dstlen*/);
				SynthesisFilterHaar (tmpLo+j, tmpHi+j, dst+j, bufferLo, bufferHi, m_w, m_h);
			}
		}
		
		
		delete[] tmpLo;
		delete[] tmpHi;
		delete[] bufferLo;
		delete[] bufferHi;
		
	}

	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
};

#endif
