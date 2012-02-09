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
#include "gauss.h"

namespace rtengine {
	
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) > (b) ? (b) : (a))
#define SQR(x) ((x)*(x))

	
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
		
		// spacing of filter taps
		size_t skip;
		
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
		
		void AnalysisFilter (T * srcbuffer, T * dstLo, T * dstHi, float *filterLo, float *filterHi, 
							 int taps, int offset, int pitch, int srclen);
		void SynthesisFilter (T * srcLo, T * srcHi, T * dst, T *bufferLo, T *bufferHi, 
							  float *filterLo, float *filterHi, int taps, int offset, int pitch, int dstlen);
		
		void imp_nr (T* src, int width, int height, double thresh);

		
	public:
		
		T ** wavcoeffs;
		
		template<typename E>
		wavelet_level(E * src, int level, int padding, size_t w, size_t h, float *filterV, float *filterH, int len, int offset)
		: m_w(w), m_h(h), m_w2(w), m_h2(h), m_pad(padding), wavcoeffs(NULL), lvl(level), skip(1<<level)
		{
			m_w2 = (w+2*skip*padding);
			m_h2 = (h+2*skip*padding);
			m_pad= skip*padding;
			
			wavcoeffs = create((m_w2)*(m_h2));
			decompose_level(src, filterV, filterH, len, offset, skip);
			
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
		
		template<typename E>
		void decompose_level(E *src, float *filterV, float *filterH, int len, int offset, int skip);
		
		template<typename E>
		void reconstruct_level(E *dst, float *filterV, float *filterH, int len, int offset, int skip);
		
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
	
	template<typename T> template<typename E>
	void wavelet_level<T>::loadbuffer(E * src, E * dst, int pitch, int srclen)
	{
		E * tmp = dst + m_pad;
		memset(dst, 0, (MAX(m_w2,m_h2))*sizeof(E));
		
		/*int cosetlen = (srclen+1)/skip;
		
		//create buffer with 'skip' rows and 'cosetlen' columns from src data
		//'skip' is the spacing of taps on the wavelet filter to be applied to src rows/columns
		//therefore there are 'skip' cosets of the row/column data, each of length 'cosetlen'
		//'pitch' is 1 for rows, W for columns
		for (size_t i = 0, j = 0; i<srclen; i++, j += pitch)
		{
			int coset = i%skip;
			int indx  = i/skip;
			tmp[coset*cosetlen + indx] = src[j];
		}	
		
		//even up last row/column if srclen is not a multiple of 'skip'
		for (size_t i=srclen; i<srclen+(srclen%skip); i++) {
			tmp[i] = tmp[i-skip];
		}
		
		// extend each coset mirror-like by padding amount 'm_pad'
		for (size_t coset=0; coset<skip*cosetlen; coset+=cosetlen) {
			for (size_t i=1; i<=MIN(cosetlen-1,m_pad); i++) {
				tmp[coset-i] = tmp[coset+i];
				tmp[coset+cosetlen+i-1] = tmp[coset+cosetlen-i-1];
			}
		}*/
		
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
		 * aligning the 'offset' element of the filter
		 * with the input pixel, and skipping 'pitch' pixels
		 * between taps (eg pitch=1 for horizontal filtering, 
		 * pitch=W for vertical, pitch=W+1,W-1 for diagonals.
		 * Currently diagonal filtering is not supported
		 * for the full source array, until a more sophisticated 
		 * treatment of mirror BC's is implemented.
		 *
		 */
		
		//input data is 'skip' rows and cosetlen=srclen/skip columns (which includes padding at either and)
		
		/*int cosetlen = srclen/skip;
		
		for (size_t coset=0; coset<srclen; coset+=cosetlen) {
			for (size_t i = 0; i < (cosetlen); i++) {
				float lo=0,hi=0;
				if (i>taps && i<cosetlen-taps) {//bulk
					for (int j=0, l=-offset; j<taps; j++, l++) {
						lo += filterLo[j] * src[i-l];//lopass channel
						hi += filterHi[j] * src[i-l];//hipass channel
					}
				} else {//boundary
					for (int j=0; j<taps; j++) {
						int arg = MAX(0,MIN(i+(offset-j),srclen-1));//clamped BC's
						lo += filterLo[j] * src[arg];//lopass channel
						hi += filterHi[j] * src[arg];//hipass channel
					}
				}
			
				dstLo[(pitch*(coset+i))] = lo;
				dstHi[(pitch*(coset+i))] = hi;
			}
		}*/
				
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
	
	template<typename T>
	void wavelet_level<T>::SynthesisFilter (T * srcLo, T * srcHi, T * dst, T *bufferLo, T *bufferHi, float *filterLo, 
												 float *filterHi, int taps, int offset, int pitch, int dstlen) {
		
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

		
		
		// load into buffer
		/*
		int srclen=(dstlen+(dstlen%skip)+2*m_pad);	//length of row/col in src (coarser level)
		int cosetlen = srclen/skip;						//length of coset (skip is spacing of taps in filter)

		for (size_t i=0, j=0; i<srclen; i++, j+=pitch) {
			int indx = (i%skip)*cosetlen + i/skip;
			bufferLo[indx]=srcLo[j];
			bufferHi[indx]=srcHi[j];
		}
		
		for (size_t coset=0; coset<srclen; coset+=cosetlen) {
			for (size_t i = m_pad; i < (cosetlen-m_pad); i++) {
				float tot=0;
				if (i>taps && i<(cosetlen-taps)) {//bulk
					for (int j=0, l=-shift; j<taps; j++, l++) {
						tot += (filterLo[j] * bufferLo[i-l] + filterHi[j] * bufferHi[i-l]);
					}
				} else {//boundary
					if (coset+i-m_pad == srclen) return;
					for (int j=0, l=-shift; j<taps; j++, l++) {
						int arg = MAX(0,MIN((i-l),srclen-1));//clamped BC's
						tot += (filterLo[j] * bufferLo[arg] + filterHi[j] * bufferHi[arg]);
					}
				}
				
				dst[pitch*(coset+i-m_pad)] = tot;
			}
		}*/
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		// load into buffer
		
		int srclen = (dstlen==m_w ? m_w2 : m_h2);//length of row/col in src (coarser level)
		
		for (size_t i=0, j=0; i<srclen; i++, j+=pitch) {
			bufferLo[i]=srcLo[j];
			bufferHi[i]=srcHi[j];
		}
		
		int shift=(taps-offset-1);
		for(size_t i = m_pad; i < (dstlen+m_pad); i++) {
			float tot=0;
			if (i>skip*taps && i<(srclen-skip*taps)) {//bulk
				for (int j=0, l=-skip*shift; j<taps; j++, l+=skip) {
					tot += (filterLo[j] * bufferLo[i-l] + filterHi[j] * bufferHi[i-l]);
				}
			} else {//boundary
				for (int j=0, l=-skip*shift; j<taps; j++, l+=skip) {
					int arg = MAX(0,MIN((i-l),srclen-1));//clamped BC's
					tot += (filterLo[j] * bufferLo[arg] + filterHi[j] * bufferHi[arg]);
				}
			}
			
			dst[pitch*(i-m_pad)] = tot;
			if (tot<0.0f || tot>65535.0f) {
				float xxx=tot;
				float yyy=1.0f;
			}
		}
		
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	template<typename T> template<typename E>
	void wavelet_level<T>::decompose_level(E *src, float *filterV, float *filterH, int taps, int offset, int skip) { 
		
		T *tmpLo = new T[m_w*m_h2];
		T *tmpHi = new T[m_w*m_h2];
		
		T *buffer = new T[MAX(m_w2,m_h2)];
		
		/* filter along columns */
		for (int j=0; j<m_w; j++) {
			loadbuffer(src+j, buffer, m_w/*pitch*/, m_h/*srclen*/);//pad a column of data and load it to buffer
			AnalysisFilter (buffer, tmpLo+j, tmpHi+j, filterV, filterV+taps, taps, offset, m_w/*output_pitch*/, m_h/*srclen*/);
		}
		
		/* filter along rows */
		for (int i=0; i<m_h2; i++) {
			loadbuffer(tmpLo+i*m_w, buffer, 1/*pitch*/, m_w/*srclen*/);//pad a row of data and load it to buffer
			AnalysisFilter (buffer, wavcoeffs[0]+i*m_w2, wavcoeffs[1]+i*m_w2, filterH, filterH+taps, taps, offset, 1/*output_pitch*/, m_w/*srclen*/);
			loadbuffer(tmpHi+i*m_w, buffer, 1/*pitch*/, m_w/*srclen*/);
			AnalysisFilter (buffer, wavcoeffs[2]+i*m_w2, wavcoeffs[3]+i*m_w2, filterH, filterH+taps, taps, offset, 1/*output_pitch*/, m_w/*srclen*/);
		}
		
		//imp_nr (wavcoeffs[0], m_w2, m_h2, 50.0f/20.0f);

		delete[] tmpLo;
		delete[] tmpHi;
		delete[] buffer;
	}
	
	/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
	
	template<typename T> template<typename E>
	void wavelet_level<T>::reconstruct_level(E *dst, float *filterV, float *filterH, int taps, int offset, int skip) { 
		
		T *tmpLo = new T[m_w*m_h2];
		T *tmpHi = new T[m_w*m_h2];
		
		int buflen = MAX(m_w2,m_h2);
		float *bufferLo = new float[buflen];
		float *bufferHi = new float[buflen];
		
		/* filter along rows */
		for (int i=0; i<m_h2; i++) {
			
			SynthesisFilter (wavcoeffs[0]+i*m_w2, wavcoeffs[1]+i*m_w2, tmpLo+i*m_w, bufferLo, bufferHi,  
							 filterH, filterH+taps, taps, offset, 1/*pitch*/, m_w/*dstlen*/);
			SynthesisFilter (wavcoeffs[2]+i*m_w2, wavcoeffs[3]+i*m_w2, tmpHi+i*m_w, bufferLo, bufferHi,  
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
	
	template<typename T> 
	void wavelet_level<T>::imp_nr (T* src, int width, int height, double thresh) {
		
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// impulse noise removal
		// local variables
		
		float hpfabs, hfnbrave;
		const float eps = 0.01;
		
		// buffer for the lowpass image
		float * lpf = new float[width*height];
		// buffer for the highpass image
		float * impish = new float[width*height];
		
		//The cleaning algorithm starts here
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// modified bilateral filter for lowpass image, omitting input pixel; or Gaussian blur
		/*
		static float eps = 1.0;
		float wtdsum[3], dirwt, norm;
		int i1, j1;	
				
		AlignedBuffer<double>* buffer = new AlignedBuffer<double> (MAX(width,height));
		
		gaussHorizontal<float> (src, lpf, buffer, width, height, MAX(2.0,thresh-1.0), false);
		gaussVertical<float>   (lpf, lpf, buffer, width, height, MAX(2.0,thresh-1.0), false);
		
		delete buffer;
		*/
		
		boxblur(src, lpf, 2, 2, width, height);
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		float impthr = MAX(1.0,5.5-thresh);
		
		for (int i=0; i < height; i++) 
			for (int j=0; j < width; j++) {
				hpfabs = fabs(src[i*width+j]-lpf[i*width+j]);
				//block average of high pass data
				for (int i1=MAX(0,i-2), hfnbrave=0; i1<=MIN(i+2,height-1); i1++ )
					for (int j1=MAX(0,j-2); j1<=MIN(j+2,width-1); j1++ ) {
						hfnbrave += fabs(src[i1*width+j1]-lpf[i1*width+j1]);
					}
				hfnbrave = (hfnbrave-hpfabs)/24;
				hpfabs>(hfnbrave*impthr) ? impish[i*width+j]=1 : impish[i*width+j]=0;
				
			}//now impulsive values have been identified
		
		for (int i=0; i < height; i++)
			for (int j=0; j < width; j++) {
				if (!impish[i*width+j]) continue;
				float norm=0.0;
				float wtdsum=0.0;
				for (int i1=MAX(0,i-2), hfnbrave=0; i1<=MIN(i+2,height-1); i1++ )
					for (int j1=MAX(0,j-2); j1<=MIN(j+2,width-1); j1++ ) {
						if (i1==i && j1==j) continue;
						if (impish[i1*width+j1]) continue;
						float dirwt = 1/(SQR(src[i1*width+j1]-src[i*width+j])+eps);//use more sophisticated rangefn???
						wtdsum += dirwt*src[i1*width+j1];
						norm += dirwt;
					}
				//wtdsum /= norm;
				if (norm) {
					src[i*width+j]=wtdsum/norm;//low pass filter
				} 
				
			}//now impulsive values have been corrected
		
	delete [] lpf;
	delete [] impish;
	
	}

	
};

#endif
