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


/*template<typename T>
class limiter //for limiting output between specified bounds
{
    T min_value, max_value;
public:
    limiter(T min, T max)
    : min_value(min), max_value(max)
    {}
    
    T operator()(T x)
    {
        if(x < min_value)
            return min_value;
        if(x > max_value)
            return max_value;
        return x;
    }
};*/

/*template<typename T>
class noop
{
public:
    T operator()(T x)
    {
        return x;
    }
};*/

/*template<typename T>
inline T clip(T x, T min_value, T max_value)
{
    if(x < min_value)
        return min_value;
    if(x > max_value)
        return max_value;
    return x;
}*/

/*template <typename A, typename B>
void plane_copy(A ** a, B * b, size_t datalen)
{
	for (size_t i=0; i<datalen; i++) {
		b[j] = static_cast<B> (0.25f*(a[0][j]+a[1][j]+a[2][j]+a[3][j]))
	}
}*/

//////////////////////////////////////////////////////////////////////////////

template<typename T>
class cplx_wavelet_level
{
    // full size
    size_t m_w, m_h;

    // size of low frequency part
    size_t m_w2, m_h2;

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

    //void dwt_2d(size_t w, size_t h);
    //void idwt_2d(size_t w, size_t h, int alpha);
	
	void AnalysisFilter (T * src, T * dstLo, T * dstHi, T *buffer, float *filterLo, float *filterHi, 
						 int taps, int offset, int pitch, int srclen);
	void SynthesisFilter (T * srcLo, T * srcHi, T * dst, T *bufferLo, T *bufferHi, 
						  float *filterLo, float *filterHi, int taps, int offset, int pitch, int dstlen);
		
public:
	
	T ** wavcoeffs;

    template<typename E>
    cplx_wavelet_level(E * src, size_t w, size_t h, float *filterV, float *filterH, int len, int offset)
    : m_w(w), m_h(h), m_w2((w+1)/2), m_h2((h+1)/2), 
      wavcoeffs(NULL)//,m_coeffs(NULL), m_weights_rows(NULL), m_weights_cols(NULL)
    {

        //m_coeffs = create(w, h);
        //m_weights_rows = create(w + 4, h);
        //m_weights_cols = create(h + 4, w);
        
        //decompose_level(src, w, h, wavcoeffs, float **filterV, float **filterH, int len, int offset);
		
		wavcoeffs = create(m_w2*m_h2);
		decompose_level(src, filterV, filterH, len, offset);

    }
    
    ~cplx_wavelet_level()
    {
        //destroy(m_coeffs);
        //destroy(m_weights_rows);
        //destroy(m_weights_cols);
		destroy(wavcoeffs);
    }
    
    T ** subbands() const
    {
        return wavcoeffs;//m_coeffs;
    }
	
	T * lopass() const
    {
        return wavcoeffs[0];//m_coeffs;
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
/*
template<typename T>
void wavelet_level<T>::dwt_2d(size_t w, size_t h)
{
    T * buffer = new T[std::max(w, h) + 4];
    
    for(size_t j = 0; j < h; j++)
    {
        //dwt_haar(m_coeffs[j], 1, buffer, w);
        //dwt_53(m_coeffs[j], 1, buffer, w);
        dwt_wcdf(m_coeffs[j], 1, buffer, w, m_weights_rows[j]);
    }
    
    for(size_t i = 0; i < w; i++)
    {
        //dwt_haar(&m_coeffs[0][i], m_pitch, buffer, h);
        //dwt_53(&m_coeffs[0][i], w, buffer, h);
        dwt_wcdf(&m_coeffs[0][i], w, buffer, h, m_weights_cols[i]);
    }
    
    delete[] buffer;
}

template<typename T>
void wavelet_level<T>::idwt_2d(size_t w, size_t h, int alpha)
{
    T * buffer = new T[std::max(w, h) + 4];
    
    for(size_t i = 0; i < w; i++)
    {
        //idwt_haar(&m_coeffs[0][i], m_pitch, buffer, h, alpha);
        //idwt_53(&m_coeffs[0][i], w, buffer, h, alpha);
        idwt_wcdf(&m_coeffs[0][i], w, buffer, h, alpha, m_weights_cols[i]);
        //idwt_noop(&m_coeffs[0][i], w, buffer, h, alpha);
    }
    
    for(size_t j = 0; j < h; j++)
    {
        //idwt_haar(m_coeffs[j], 1, buffer, w, alpha);
        //idwt_53(m_coeffs[j], 1, buffer, w, alpha);
        idwt_wcdf(m_coeffs[j], 1, buffer, w, alpha, m_weights_rows[j]);
        //idwt_noop(m_coeffs[j], 1, buffer, w, alpha);
    }

    delete[] buffer;
}
*/
	

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

template<typename T>
void cplx_wavelet_level<T>::destroy(T ** subbands)
{
    if(subbands)
    {
        delete[] subbands[0];
    
        delete[] subbands;
    }
}

	
/*template<typename T> template<typename E>
void wavelet_level<T>::decompose(E ** src)
{
    noop<T> l;

    plane_copy(src, m_coeffs, m_w, m_h, l);

    dwt_2d(m_w, m_h);
}

template<typename T> template<typename E, typename L>
void wavelet_level<T>::reconstruct(E ** dst, int alpha, L & l)
{
    idwt_2d(m_w, m_h, alpha);

    plane_copy(m_coeffs, dst, m_w, m_h, l);
}*/
	
	
	
	
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	template<typename T>
	void cplx_wavelet_level<T>::AnalysisFilter (T * src, T * dstLo, T * dstHi, T *buffer, float *filterLo, float *filterHi, 
						 int taps, int offset, int pitch, int srclen) {
		
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
		 * Destination arrays must be initialized to zero.
		 */
		
		T * tmp = buffer + taps;//offset
		
		// copy data
		
		for(size_t i = 0, j = 0; i < srclen; i++, j += pitch)
		{
			tmp[i] = src[j];
		}
		
		// extend mirror-like
		
		for (size_t i=-1; i!=-offset; i--) {
			tmp[i] = tmp[-i];
		}
		for (size_t i=0; i<taps-offset; i++) {
			tmp[srclen+i] = tmp[srclen-i-2];
		}
		
		// calculate coefficients
		
		for(ptrdiff_t i = 0; i < (ptrdiff_t)srclen; i+=2) {
			float lo=0,hi=0;
			for (int j=0; j<taps; j++) {
				lo += filterLo[j] * tmp[i-offset+j];//lopass channel
				hi += filterHi[j] * tmp[i-offset+j];//hipass channel
			}
			dstLo[i] = lo;
			dstHi[i] = hi;
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
		 * Destination arrays must be initialized to zero.
		 */
		
		T * tmpLo = bufferLo + taps;//offset
		T * tmpHi = bufferHi + taps;//offset
		
		// copy data
		
		for(size_t i = 0, j = 0; i < dstlen; i++, j += pitch)
		{
			tmpLo[2*i] = srcLo[j];
			tmpHi[2*i] = srcHi[j];
		}
		
		// extend mirror-like
		
		for (size_t i=-1; i!=-offset; i--) {
			tmpLo[2*i] = tmpLo[-i];
			tmpHi[2*i] = tmpHi[-i];
		}
		for (size_t i=0; i<taps-offset; i++) {
			tmpLo[2*(dstlen+i)] = tmpLo[dstlen-i-2];
			tmpHi[2*(dstlen+i)] = tmpHi[dstlen-i-2];
		}
		
		// calculate coefficients
		
		for(ptrdiff_t i = 0; i < (ptrdiff_t)dstlen; i++) {
			float tot=0;
			for (int j=0; j<taps; j++) {
				tot += 0.5*(filterLo[j] * tmpLo[i-offset+j] + filterHi[j] * tmpHi[i-offset+j]);//lopass channel
			}
			dst[i] = tot;
		}
		
		
	}
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	template<typename T> template<typename E>
	void cplx_wavelet_level<T>::decompose_level(E *src, float *filterV, float *filterH, int taps, int offset) { 
		
		//int hfw = (W+1)/2;
		//int hfh = (H+1)/2;
		T *tmpLo = new T(m_w*m_h2);
		T *tmpHi = new T(m_w*m_h2);
		
		T *buffer = new T[MAX(m_w,m_h)+taps];
		
		/* filter along columns */
		for (int j=0; j<m_w; j++) {
			AnalysisFilter (src+j, tmpLo+j, tmpHi+j, buffer, filterV, filterV+taps, taps, offset, m_w/*pitch*/, m_h/*srclen*/);
		}
		
		/* filter along rows */
		for (int i=0; i<m_h2; i++) {
			AnalysisFilter (tmpLo+i*m_w, wavcoeffs[0]+i*m_w2, wavcoeffs[1]+i*m_w2, buffer, filterH, filterH+taps, taps, offset, 1/*pitch*/, m_w/*srclen*/);
			AnalysisFilter (tmpHi+i*m_w, wavcoeffs[2]+i*m_w2, wavcoeffs[3]+i*m_w2, buffer, filterH, filterH+taps, taps, offset, 1/*pitch*/, m_w/*srclen*/);
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
		array2D<float> tmpLo(m_w2,m_h);
		array2D<float> tmpHi(m_w2,m_h);
		
		float *bufferLo = new float[MAX(m_w,m_h)+taps];
		float *bufferHi = new float[MAX(m_w,m_h)+taps];
		//bufferLo = (float (*)) calloc (MAX(m_w,m_h)+taps, sizeof *bufferLo);
		//bufferHi = (float (*)) calloc (MAX(m_w,m_h)+taps, sizeof *bufferHi);
		
		/* filter along columns */
		for (int j=0; j<m_w2; j++) {
			SynthesisFilter (wavcoeffs[0], wavcoeffs[1], tmpLo, bufferLo, bufferHi, 
							 filterV, filterV+taps, taps, offset, m_w2/*pitch*/, m_h2/*srclen*/);
			SynthesisFilter (wavcoeffs[2], wavcoeffs[3], tmpLo, bufferLo, bufferHi, 
							 filterV, filterV+taps, taps, offset, m_w2/*pitch*/, m_h2/*srclen*/);
		}
		
		/* filter along rows */
		for (int i=0; i<m_h2; i++) {
			SynthesisFilter (tmpLo, tmpHi, dst, bufferLo, bufferHi, 
							 filterH, filterH+taps, taps, offset, 1/*pitch*/, m_w2/*srclen*/);
		}
		
		//free (bufferLo);
		//free (bufferHi);
		
		delete[] bufferLo;
		delete[] bufferHi;
		
	}
	
	/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
	/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
	

};

#endif
