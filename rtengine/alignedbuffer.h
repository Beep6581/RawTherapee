/*
 *  This file is part of RawTherapee.
 *
*  Copyright (c) 2004-2012 Gabor Horvath <hgabor@rawtherapee.com>, Oliver Duis <oduis@oliverduis.de>
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
#include <stdint.h>
#include <vector>
#include <glibmm.h>

// Aligned buffer that should be faster
template <class T> class AlignedBuffer {

    private:
      T* real ;
      
    public:
      T* data ;
    bool inUse;

        AlignedBuffer (size_t size, size_t align=16) {
            real = new T[size+2*align];
            data = (T*)((uintptr_t)real + (align-((uintptr_t)real)%align));
        inUse=true;
        }

        ~AlignedBuffer () {
            delete [] real;
        }
};

// Multi processor version, use with OpenMP
template <class T> class AlignedBufferMP {
private:
    Glib::Mutex mtx;
    std::vector<AlignedBuffer<T>*> buffers;
    size_t size;

public:
    AlignedBufferMP(size_t sizeP) {
        size=sizeP;
    }

    ~AlignedBufferMP() {
        for (int i=0;i<buffers.size();i++) delete buffers[i];
    }

    AlignedBuffer<T>* acquire() {
        Glib::Mutex::Lock lock(mtx);

        // Find available buffer
        for (int i;i<buffers.size();i++) {
            if (!buffers[i]->inUse) {
                buffers[i]->inUse=true;
                return buffers[i];
            }
        }

        // Add new buffer if nothing is free
        AlignedBuffer<T>* buffer=new AlignedBuffer<T>(size);
        buffers.push_back(buffer);

        return buffer;
    }

    void release(AlignedBuffer<T>* buffer) {
        Glib::Mutex::Lock lock(mtx);

        buffer->inUse=false;
    }
};
#endif
