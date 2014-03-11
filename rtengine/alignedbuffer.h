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
#include <cstdlib>
#include <vector>
#include <glibmm.h>
#include "../rtgui/threadutils.h"

// Aligned buffer that should be faster
template <class T> class AlignedBuffer {

private:
    void* real ;
    char alignment;
    size_t allocatedSize;

public:
    T* data ;
    bool inUse;

    /* size is the number of elements of size T, i.e. allocated size will be sizeof(T)*size ; set it to 0 if you want to defer the allocation
     * align is expressed in bytes; SSE instructions need 128 bits alignment, which mean 16 bytes, which is the default value
     */
    AlignedBuffer (size_t size=0, size_t align=16) : real(NULL), alignment(align), allocatedSize(0), data(NULL), inUse(false) {
        if (size)
            resize(size);
    }

    ~AlignedBuffer () {
        if (real) free(real);
    }

    /* Allocate the the "size" amount of elements of "structSize" length each
     * params:
     * @size: number of elements to allocate
     * @structSize: if non null, will let you override the default struct's size (unit: byte)
     */
    bool resize(size_t size, int structSize=0) {
        if (allocatedSize != size) {
            if (!size) {
                // The user want to free the memory
                if (real) free(real);
                real = NULL;
                data = NULL;
                inUse = false;
            }
            else {
                int sSize = structSize ? structSize : sizeof(T);
                allocatedSize = size*sSize;
                real = realloc(real, allocatedSize+alignment);
                if (real) {
                    //data = (T*)( (uintptr_t)real + (alignment-((uintptr_t)real)%alignment) );
                    data = (T*)( ( uintptr_t(real) + uintptr_t(alignment-1)) / alignment * alignment);
                    inUse = true;
                }
                else {
                    allocatedSize = 0;
                    data = NULL;
                    inUse = false;
                    return false;
                }
            }
        }
        return true;
    }

    void swap(AlignedBuffer<T> &other) {
        void *tmpReal = other.real;
        other.real = real;
        real = tmpReal;

        char tmpAlignt = other.alignment;
        other.alignment = alignment;
        alignment = tmpAlignt;

        size_t tmpAllocSize = other.allocatedSize;
        other.allocatedSize = allocatedSize;
        allocatedSize = tmpAllocSize;

        T* tmpData = other.data;
        other.data = data;
        data = tmpData;

        bool tmpInUse = other.inUse;
        other.inUse = inUse;
        inUse = tmpInUse;
    }
};

// Multi processor version, use with OpenMP
template <class T> class AlignedBufferMP {
private:
    MyMutex mtx;
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
        MyMutex::MyLock lock(mtx);

        // Find available buffer
        for (int i=0;i<buffers.size();i++) {
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
    	MyMutex::MyLock lock(mtx);

        buffer->inUse=false;
    }
};
#endif
