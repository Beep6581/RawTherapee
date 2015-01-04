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
namespace rtengine {

	template<typename T>
	class wavelet_level
	{

		// size of padded border
		size_t m_pad;

		// level of decomposition
		int lvl;

		// whether to subsample the output
		bool subsamp_out;

		// spacing of filter taps
		int skip;

		// allocation and destruction of data storage
		T ** create(size_t n);
		void destroy(T ** subbands);

		// load a row/column of input data, possibly with padding

		void AnalysisFilterHaarVertical (T * srcbuffer, T * dstLo, T * dstHi, int pitch, int srclen, int row);
		void AnalysisFilterHaarHorizontal (T * srcbuffer, T * dstLo, T * dstHi, int srclen, int row);
		void SynthesisFilterHaarHorizontal (T * srcLo, T * srcHi, T * dst, int dstlen);
		void SynthesisFilterHaarVertical (T * srcLo, T * srcHi, T * dst, int pitch, int dstlen);

		void AnalysisFilterSubsampHorizontal (T * srcbuffer, T * dstLo, T * dstHi, float *filterLo, float *filterHi, 
							 int taps, int offset, int pitch, int srclen, int m_w2, int row);
		void AnalysisFilterSubsampVertical (T * srcbuffer, T * dstLo, T * dstHi, float *filterLo, float *filterHi, 
							 int taps, int offset, int pitch, int srclen, int row);
		void SynthesisFilterSubsampHorizontal (T * srcLo, T * srcHi, T * dst,
							  float *filterLo, float *filterHi, int taps, int offset, int dstlen);
		void SynthesisFilterSubsampVertical (T * srcLo, T * srcHi, T * dst, float *filterLo, float *filterHi, int taps, int offset, int pitch, int dstlen);

	public:

		T ** wavcoeffs;
		// full size
		size_t m_w, m_h;

		// size of low frequency part
		size_t m_w2, m_h2;

		template<typename E>
		wavelet_level(E * src, E * dst, int level, int subsamp, int padding, size_t w, size_t h, float *filterV, float *filterH, int len, int offset)
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
			decompose_level(src, dst, filterV, filterH, len, offset);
			
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
		void decompose_level(E *src, E *dst, float *filterV, float *filterH, int len, int offset);

		template<typename E>
		void reconstruct_level(E* tmpLo, E* tmpHi, E *src, E *dst, float *filterV, float *filterH, int taps, int offset);
	};

	template<typename T>
	T ** wavelet_level<T>::create(size_t n)	{
		T * data = new T[3*n];
		T ** subbands = new T*[4];
		for(size_t j = 1; j < 4; j++) {
			subbands[j] = data + n * (j-1);
		}
		return subbands;
	}

	template<typename T>
	void wavelet_level<T>::destroy(T ** subbands) {
		if(subbands) {
			delete[] subbands[1];
			delete[] subbands;
		}
	}

	template<typename T>
	void wavelet_level<T>::AnalysisFilterHaarHorizontal (T * RESTRICT srcbuffer, T * RESTRICT dstLo, T * RESTRICT dstHi, int srclen, int row) {
		/* Basic convolution code
		 * Applies a Haar filter 
		*/							
			for(int i = 0; i < (srclen - skip); i++) {
				dstLo[row*srclen+i] = (srcbuffer[i] + srcbuffer[i+skip]);
				dstHi[row*srclen+i] = (srcbuffer[i] - srcbuffer[i+skip]);
			}
			for(size_t i = max(srclen-skip,skip); i < (srclen); i++) {
				dstLo[row*srclen+i] = (srcbuffer[i] + srcbuffer[i-skip]);
				dstHi[row*srclen+i] = (srcbuffer[i] - srcbuffer[i-skip]);
			}
	}

	template<typename T> void wavelet_level<T>::AnalysisFilterHaarVertical (T * RESTRICT srcbuffer, T * RESTRICT dstLo, T * RESTRICT dstHi, int pitch, int srclen, int row) {
	/* Basic convolution code
	 * Applies a Haar filter 
	*/
		if(row < (srclen - skip)) {
			for(int j=0;j<pitch;j++) {
				dstLo[j] = 0.25f*(srcbuffer[row*pitch+j] + srcbuffer[(row+skip)*pitch+j]);
				dstHi[j] = 0.25f*(srcbuffer[row*pitch+j] - srcbuffer[(row+skip)*pitch+j]);
			}
		} else if(row>=max(srclen-skip,skip)) {
			for(int j=0;j<pitch;j++) {
				dstLo[j] = 0.25f*(srcbuffer[row*pitch+j] + srcbuffer[(row-skip)*pitch+j]);
				dstHi[j] = 0.25f*(srcbuffer[row*pitch+j] - srcbuffer[(row-skip)*pitch+j]);
			}
		}
	}

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	
	template<typename T> void wavelet_level<T>::SynthesisFilterHaarHorizontal (T * RESTRICT srcLo, T * RESTRICT srcHi, T * RESTRICT dst, int dstlen) {

		/* Basic convolution code
		 * Applies a Haar filter 
		 *
		 */

		for (int k=0; k<m_h2; k++) {
			for(size_t i = (m_pad); i < (m_pad+skip); i++) {
				dst[k*m_w+(i-m_pad)] = (srcLo[k*m_w2+i] + srcHi[k*m_w2+i]);			
			}
			for(size_t i = m_pad+skip; i < (dstlen+m_pad); i++) {
				dst[k*m_w+(i-m_pad)] = 0.5f*(srcLo[k*m_w2+i] + srcHi[k*m_w2+i] + srcLo[k*m_w2+i-skip] - srcHi[k*m_w2+i-skip]);			
			}
		}
	}

	template<typename T> void wavelet_level<T>::SynthesisFilterHaarVertical (T * RESTRICT srcLo, T * RESTRICT srcHi, T * RESTRICT dst, int pitch, int dstlen) {

		/* Basic convolution code
		 * Applies a Haar filter 
		 *
		 */

		for(size_t i = (m_pad); i < (m_pad+skip); i++) {
			for(int j=0;j<pitch;j++)
				dst[pitch*(i-m_pad)+j] = (srcLo[i*pitch+j] + srcHi[i*pitch+j]);			
		}
		for(size_t i = m_pad+skip; i < (dstlen+m_pad); i++) {
			for(int j=0;j<pitch;j++)
				dst[pitch*(i-m_pad)+j] = 0.5f*(srcLo[i*pitch+j] + srcHi[i*pitch+j] + srcLo[(i-skip)*pitch+j] - srcHi[(i-skip)*pitch+j]);			
		}
	}

	template<typename T>
	void wavelet_level<T>::AnalysisFilterSubsampHorizontal (T * RESTRICT srcbuffer, T * RESTRICT dstLo, T * RESTRICT dstHi, float * RESTRICT filterLo, float *filterHi, 
												int taps, int offset, int pitch, int srclen, int m_w2, int row) {
		/* Basic convolution code
		 * Applies an FIR filter 'filter' with filter length 'taps', 
		 * aligning the 'offset' element of the filter with
		 * the input pixel, and skipping 'skip' pixels between taps 
		 * Output is subsampled by two
		 */
		// calculate coefficients
			for(int i = 0; i < srclen; i+=2) {
				float lo = 0.f, hi = 0.f;
				if (LIKELY(i>skip*taps && i<srclen-skip*taps)) {//bulk
					for (int j=0, l=-skip*offset; j<taps; j++, l+=skip) {
						float src = srcbuffer[i-l];
						lo += filterLo[j] * src;//lopass channel
						hi += filterHi[j] * src;//hipass channel
					}
				} else {
					for (int j=0; j<taps; j++) {
						int arg = max(0,min(i+skip*(offset-j),srclen-1));//clamped BC's
						lo += filterLo[j] * srcbuffer[arg];//lopass channel
						hi += filterHi[j] * srcbuffer[arg];//hipass channel
					}
				}
				dstLo[row*m_w2+((i/2))] = lo;
				dstHi[row*m_w2+((i/2))] = hi;
			}
	}

	template<typename T> void wavelet_level<T>::AnalysisFilterSubsampVertical (T * RESTRICT srcbuffer, T * RESTRICT dstLo, T * RESTRICT dstHi, float * RESTRICT filterLo, float * RESTRICT filterHi, 
													int taps, int offset, int pitch, int srclen, int row) {

			/* Basic convolution code
			 * Applies an FIR filter 'filter' with filter length 'taps', 
			 * aligning the 'offset' element of the filter with
			 * the input pixel, and skipping 'skip' pixels between taps 
			 * Output is subsampled by two
			 */

			// calculate coefficients
			if (LIKELY(row>skip*taps && row<srclen-skip*taps)) {//bulk
				for (int k=0; k<pitch; k++) {
					float lo = 0.f, hi = 0.f;
					for (int j=0, l=-skip*offset; j<taps; j++, l+=skip) {
						lo += filterLo[j] * srcbuffer[(row-l)*pitch+k];//lopass channel
						hi += filterHi[j] * srcbuffer[(row-l)*pitch+k];//hipass channel
					}
					dstLo[k] = lo;
					dstHi[k] = hi;
				}
			} else {//boundary
				for (int k=0; k<pitch; k++) {
					float lo = 0.f, hi = 0.f;
					for (int j=0; j<taps; j++) {
						int arg = max(0,min(row+skip*(offset-j),srclen-1))*pitch+k;//clamped BC's
						lo += filterLo[j] * srcbuffer[arg];//lopass channel
						hi += filterHi[j] * srcbuffer[arg];//hipass channel
					}
					dstLo[k] = lo;
					dstHi[k] = hi;
				}
			}
	}

	template<typename T> void wavelet_level<T>::SynthesisFilterSubsampHorizontal (T * RESTRICT srcLo, T * RESTRICT srcHi, T * RESTRICT dst, float * RESTRICT filterLo, float * RESTRICT filterHi, int taps, int offset, int dstlen) {

		/* Basic convolution code
		 * Applies an FIR filter 'filter' with filter length 'taps', 
		 * aligning the 'offset' element of the filter with
		 * the input pixel, and skipping 'skip' pixels between taps 
		 * Output is subsampled by two
		 */

		// calculate coefficients
		int srclen = (dstlen==m_w ? m_w2 : m_h2);//length of row/col in src (coarser level)
		int shift = skip*(taps-offset-1);//align filter with data

		for (int k=0; k<m_h2; k++) {
			for(size_t i = m_pad; i < (dstlen+m_pad); i++) {
				float tot=0.f;
				//TODO: this is correct only if skip=1; otherwise, want to work with cosets of length 'skip'
				int i_src = (i+shift)/2;
				int begin = (i+shift)%2;
				if (LIKELY(i>skip*taps && i<(srclen-skip*taps))) {//bulk
					for (int j=begin, l=0; j<taps; j+=2, l+=skip) {
						tot += ((filterLo[j] * srcLo[k*m_w2+i_src-l] + filterHi[j] * srcHi[k*m_w2+i_src-l]));
					}
				} else {//boundary
					for (int j=begin, l=0; j<taps; j+=2, l+=skip) {
						int arg = max(0,min((i_src-l),srclen-1));//clamped BC's
						tot += ((filterLo[j] * srcLo[k*m_w2+arg] + filterHi[j] * srcHi[k*m_w2+arg]));
					}
				}
				dst[k*m_w+(i-m_pad)] = tot;
			}
		}
	}

	template<typename T> void wavelet_level<T>::SynthesisFilterSubsampVertical (T * RESTRICT srcLo, T * RESTRICT srcHi, T * RESTRICT dst, float * RESTRICT filterLo, float * RESTRICT filterHi, int taps, int offset, int pitch, int dstlen) {

			/* Basic convolution code
			 * Applies an FIR filter 'filter' with filter length 'taps', 
			 * aligning the 'offset' element of the filter with
			 * the input pixel, and skipping 'skip' pixels between taps 
			 * Output is subsampled by two
			 */

		// calculate coefficients
		int srclen = (dstlen==m_w ? m_w2 : m_h2);//length of row/col in src (coarser level)
		int shift=skip*(taps-offset-1);//align filter with data

		for(size_t i = m_pad; i < (dstlen+m_pad); i++) {
			int i_src = (i+shift)/2;
			int begin = (i+shift)%2;
			//TODO: this is correct only if skip=1; otherwise, want to work with cosets of length 'skip'
			if (LIKELY(i>skip*taps && i<(srclen-skip*taps))) {//bulk
				for (int k=0; k<pitch; k++) {
					float tot = 0.f;
					for (int j=begin, l=0; j<taps; j+=2, l+=skip) {
						tot += ((filterLo[j] * srcLo[(i_src-l)*pitch+k] + filterHi[j] * srcHi[(i_src-l)*pitch+k]));
					}
					dst[pitch*(i-m_pad)+k] = 4.f * tot;
				}
			} else {//boundary
				for (int k=0; k<pitch; k++) {
					float tot = 0.f;
					for (int j=begin, l=0; j<taps; j+=2, l+=skip) {
						int arg = max(0,min((i_src-l),srclen-1))*pitch+k;//clamped BC's
						tot += ((filterLo[j] * srcLo[arg] + filterHi[j] * srcHi[arg]));
					}
					dst[pitch*(i-m_pad)+k] = 4.f * tot;
				}
			}
		}
	}

	template<typename T> template<typename E> void wavelet_level<T>::decompose_level(E *src, E *dst, float *filterV, float *filterH, int taps, int offset) { 

		T tmpLo[m_w] ALIGNED64;
		T tmpHi[m_w] ALIGNED64;
		/* filter along rows and columns */
		if(subsamp_out) {
			for(int row=0;row<m_h;row+=2) {
				AnalysisFilterSubsampVertical (src, tmpLo, tmpHi, filterV, filterV+taps, taps, offset, m_w/*output_pitch*/, m_h/*srclen*/, row);
				AnalysisFilterSubsampHorizontal (tmpLo, dst, wavcoeffs[1], filterH, filterH+taps, taps, offset, m_h2/*output_pitch*/, m_w/*srclen*/, m_w2, row/2);
				AnalysisFilterSubsampHorizontal (tmpHi, wavcoeffs[2], wavcoeffs[3], filterH, filterH+taps, taps, offset, m_h2/*output_pitch*/, m_w/*srclen*/, m_w2, row/2);
			}
		} else {
			for(int row=0;row<m_h;row++) {
				AnalysisFilterHaarVertical (src, tmpLo, tmpHi, m_w, m_h, row);
				AnalysisFilterHaarHorizontal (tmpLo, dst, wavcoeffs[1], m_w, row);
				AnalysisFilterHaarHorizontal (tmpHi, wavcoeffs[2], wavcoeffs[3], m_w, row);
			}
		}
	}

	template<typename T> template<typename E> void wavelet_level<T>::reconstruct_level(E* tmpLo, E* tmpHi, E * src, E *dst, float *filterV, float *filterH, int taps, int offset) { 

		/* filter along rows and columns */
		if (subsamp_out) {
			SynthesisFilterSubsampHorizontal (src, wavcoeffs[1], tmpLo, filterH, filterH+taps, taps, offset, m_w/*dstlen*/);
			SynthesisFilterSubsampHorizontal (wavcoeffs[2], wavcoeffs[3], tmpHi, filterH, filterH+taps, taps, offset, m_w/*dstlen*/);
			SynthesisFilterSubsampVertical (tmpLo, tmpHi, dst, filterV, filterV+taps, taps, offset, m_w/*pitch*/, m_h/*dstlen*/);
		} else {
			SynthesisFilterHaarHorizontal (src, wavcoeffs[1], tmpLo, m_w);
			SynthesisFilterHaarHorizontal (wavcoeffs[2], wavcoeffs[3], tmpHi, m_w);
			SynthesisFilterHaarVertical (tmpLo, tmpHi, dst, m_w, m_h);
		}
	}
	
};

#endif
