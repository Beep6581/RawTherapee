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
 */

#ifndef WAVELET_LEVEL_H_INCLUDED
#define WAVELET_LEVEL_H_INCLUDED

#include <cstddef>
#include <algorithm>

#include "wavelet.h"

namespace rtengine
{

template<typename T>
class limiter
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
};

template<typename T>
class noop
{
public:
    T operator()(T x)
    {
        return x;
    }
};

template<typename T>
inline T clip(T x, T min_value, T max_value)
{
    if(x < min_value)
        return min_value;
    if(x > max_value)
        return max_value;
    return x;
}

template <typename A, typename B, typename L>
void plane_copy(A ** a, B ** b, size_t w, size_t h, L & l)
{
    for(size_t j = 0; j < h; j++)
        for(size_t i = 0; i < w; i++)
            b[j][i] = static_cast<B>(l(a[j][i]));
}

//////////////////////////////////////////////////////////////////////////////

template<typename T>
class wavelet_level
{
    // full size
    size_t m_w, m_h;

    // size of low frequency part
    size_t m_w2, m_h2;

    // array of pointers to lines of coeffs
    // actually is a single contiguous data array pointed by m_coeffs[0]
    T ** m_coeffs;
    
    // weights storage
    T ** m_weights_rows;
    T ** m_weights_cols;
    
    // allocation and destruction of data storage
    T ** create(size_t w, size_t h);
    void destroy(T ** p);

    void dwt_2d(size_t w, size_t h);
    void idwt_2d(size_t w, size_t h, int alpha);

public:

    template<typename E>
    wavelet_level(E ** src, size_t w, size_t h)
    : m_w(w), m_h(h), m_w2((w+1)/2), m_h2((h+1)/2), 
      m_coeffs(NULL), m_weights_rows(NULL), m_weights_cols(NULL)
    {
        m_coeffs = create(w, h);
        m_weights_rows = create(w + 4, h);
        m_weights_cols = create(h + 4, w);
        
        decompose(src);
    }
    
    ~wavelet_level()
    {
        destroy(m_coeffs);
        destroy(m_weights_rows);
        destroy(m_weights_cols);
    }
    
    T ** lowfreq() const
    {
        return m_coeffs;
    }
    
    size_t lfw() const
    {
        return m_w2;
    }

    size_t lfh() const
    {
        return m_h2;
    }

    template<typename E>
    void decompose(E ** src);

    template<typename E, typename L>
    void reconstruct(E ** dst, int alpha, L & limiter);
};

//////////////////////////////////////////////////////////////////////////////

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

template<typename T>
T ** wavelet_level<T>::create(size_t w, size_t h)
{
    T * data = new T[w * h];
    T ** p = new T*[h];
    for(size_t j = 0; j < h; j++)
    {
        p[j] = data + w * j;
    }
    return p;
}

template<typename T>
void wavelet_level<T>::destroy(T ** p)
{
    if(p)
    {
        delete[] p[0];
    
        delete[] p;
    }
}

template<typename T> template<typename E>
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
}

};

#endif
