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

     /** @brief Allocate aligned memory
     * @param size Number of elements of size T to allocate, i.e. allocated size will be sizeof(T)*size ; set it to 0 if you want to defer the allocation
     * @param align Expressed in bytes; SSE instructions need 128 bits alignment, which mean 16 bytes, which is the default value
     */
    AlignedBuffer (size_t size=0, size_t align=16) : real(NULL), alignment(align), allocatedSize(0), data(NULL), inUse(false) {
        if (size)
            resize(size);
    }

    ~AlignedBuffer () {
        if (real) free(real);
    }

    /** @brief Return true if there's no memory allocated
    */
    bool isEmpty() {
        return allocatedSize==0;
    }

    /** @brief Allocate the the "size" amount of elements of "structSize" length each
    * @param size number of elements to allocate
    * @param structSize if non null, will let you override the default struct's size (unit: byte)
    * @return True is everything went fine, including freeing memory when size==0, false if the allocation failed
    */
    bool resize(size_t size, int structSize=0) {
        if (allocatedSize != size) {
            if (!size) {
                // The user want to free the memory
                if (real) free(real);
                real = NULL;
                data = NULL;
                inUse = false;
                allocatedSize = 0;
            }
            else {
                int sSize = structSize ? structSize : sizeof(T);
                size_t oldAllocatedSize = allocatedSize;
                allocatedSize = size*sSize;

                // realloc were used here to limit memory fragmentation, specially when the size was smaller than the previous one.
                // But realloc copies the content to the eventually new location, which is unnecessary. To avoid this performance penalty,
                // we're freeing the memory and allocate it again if the new size is bigger.

                if (allocatedSize < oldAllocatedSize)
                    real = realloc(real, allocatedSize+alignment);
                else {
                    if (real) free (real);
                    real = malloc(allocatedSize+alignment);
                }

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
    	MyMutex::MyLock lock(mtx);

        buffer->inUse=false;
    }
};
#endif
