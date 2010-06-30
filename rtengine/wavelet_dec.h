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
    typedef int internal_type;

    static const int maxlevels = 8;
    
    int nlevels;
    size_t m_w, m_h;
    size_t m_w1, m_h1;
    
    wavelet_level<internal_type> * m_c[maxlevels];
    
public:
    wavelet_decomposition(unsigned short ** src, size_t w, size_t h);
    
    ~wavelet_decomposition();
    
    void reconstruct(unsigned short ** dst, const int * c);
};

//////////////////////////////////////////////////////////////////////////////

};

#endif
