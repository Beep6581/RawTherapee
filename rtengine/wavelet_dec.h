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

#ifndef WAVELET_DEC_H_INCLUDED
#define WAVELET_DEC_H_INCLUDED

#include <cstddef>

#include "wavelet_level.h"

namespace rtengine {

class wavelet_decomposition
{
public:

    typedef int internal_type;

private:

    static const int maxlevels = 8;
    
    int nlevels;
    size_t m_w, m_h;
    size_t m_w1, m_h1;
    
    wavelet_level<internal_type> * m_c[maxlevels];
    
public:

    template<typename E>
    wavelet_decomposition(E ** src, size_t w, size_t h);
    
    ~wavelet_decomposition();
    
    template<typename E, typename L>
    void reconstruct(E ** dst, const int * c, L & limiter);
};

//////////////////////////////////////////////////////////////////////////////

template<typename E>
wavelet_decomposition::wavelet_decomposition(E ** src, size_t w, size_t h)
: nlevels(0), m_w(w), m_h(h), m_w1(0), m_h1(0)
{
    m_w1 = w;
    m_h1 = h;
    
    m_c[0] = new wavelet_level<internal_type>(src, m_w1, m_h1);
    nlevels = 1;
    
    while(nlevels < maxlevels)
    {
        m_c[nlevels] = new wavelet_level<internal_type>(m_c[nlevels - 1]->lowfreq(), m_c[nlevels-1]->lfw(), m_c[nlevels-1]->lfh());
        nlevels ++;
    }
}


template<typename E, typename L>
void wavelet_decomposition::reconstruct(E ** dst, const int * c, L & l)
{
    noop<internal_type> n;

    for(int level = nlevels - 1; level > 0; level--)
    {
        int alpha = 1024 + 10 * c[level];
        m_c[level]->reconstruct(m_c[level-1]->lowfreq(), alpha, n);
    }
    
    int alpha = 1024 + 10 * c[0];
    m_c[0]->reconstruct(dst, alpha, l);
}


};

#endif
