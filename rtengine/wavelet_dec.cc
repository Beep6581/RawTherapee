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

#include "wavelet_dec.h"

namespace rtengine
{

wavelet_decomposition::wavelet_decomposition(unsigned short ** src, size_t w, size_t h)
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

wavelet_decomposition::~wavelet_decomposition()
{
    for(int i = 0; i < nlevels; i++)
    {
        delete m_c[i];
    }
}

void wavelet_decomposition::reconstruct(unsigned short ** dst, const int * c)
{
    for(int level = nlevels - 1; level > 0; level--)
    {
        int alpha = 1024 + 10 * c[level];
        m_c[level]->reconstruct(m_c[level-1]->lowfreq(), alpha);
    }
    
    int alpha = 1024 + 10 * c[0];
    m_c[0]->reconstruct(dst, alpha, wavelet_level<internal_type>::CLIP);
}

};

