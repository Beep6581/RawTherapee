/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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
 */
#ifndef _ALIGNEDBUFFER_
#define _ALIGNEDBUFFER_

template <class T> class AlignedBuffer {

    private:
      T* real ;
      
    public:
      T* data ;

        AlignedBuffer (size_t size, size_t align=16) {
            real = new T[size+2*align];
            data = (T*)((uintptr_t)real + (align-((uintptr_t)real)%align));
        }

        ~AlignedBuffer () {
            delete [] real;
        }
};

#endif
