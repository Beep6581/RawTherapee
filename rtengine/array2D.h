/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2011 Jan Rinze Peterzon (janrinze@gmail.com)
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

/*
 *  Declaration of flexible 2D arrays
 *
 *  Usage:
 *
 *      array2D<type> name (X-size,Y-size);
 *      array2D<type> name (X-size,Y-size,type ** data);
 *
 *      creates an array which is valid within the normal C/C++ scope "{ ... }"
 *
 *      access to elements is a simple as:
 *
 *          array2D<float> my_array (10,10); // creates 10x10 array of floats
 *          value =  my_array[3][5];
 *          my_array[4][6]=value;
 *
 *      or copy an existing 2D array
 *
 *          float ** mydata;
 *          array2D<float> my_array (10,10,mydata);
 *
 *
 *      Useful extra pointers
 *
 *          <type> ** my_array      gives access to the pointer for access with [][]
 *          <type> *  my_array      gives access to the flat stored data.
 *
 *      Advanced usage:
 *          array2D<float> my_array             ; // empty container.
 *          my_array(10,10)                     ; // resize to 10x10 array
 *          my_array(10,10,ARRAY2D_CLEAR_DATA)  ; // resize to 10x10 and clear data
 *          my_array(10,10,ARRAY2D_CLEAR_DATA|ARRAY2D_LOCK_DATA)  ; same but set a lock on changes
 *
 *          !! locked arrays cannot be resized and cannot be unlocked again !!
 */
#ifndef ARRAY2D_H_
#define ARRAY2D_H_
#include <csignal>  // for raise()
#include <cassert>

// flags for use
#define ARRAY2D_LOCK_DATA   1
#define ARRAY2D_CLEAR_DATA  2
#define ARRAY2D_BYREFERENCE 4
#define ARRAY2D_VERBOSE     8

#include <cstring>
#include <cstdio>

template<typename T>
class array2D
{

private:
    int x, y, owner, flags;
    T ** ptr;
    T * data;
    bool lock; // useful lock to ensure data is not changed anymore.
    void ar_realloc(int w, int h, int offset = 0)
    {
        if ((ptr) && ((h > y) || (4 * h < y))) {
            delete[] ptr;
            ptr = NULL;
        }

        if ((data) && (((h * w) > (x * y)) || ((h * w) < ((x * y) / 4)))) {
            delete[] data;
            data = NULL;
        }

        if (ptr == NULL) {
            ptr = new T*[h];
        }

        if (data == NULL) {
            data = new T[h * w + offset];
        }

        x = w;
        y = h;

        for (int i = 0; i < h; i++) {
            ptr[i] = data + offset + w * i;
        }

        owner = 1;
    }
public:

    // use as empty declaration, resize before use!
    // very useful as a member object
    array2D() :
        x(0), y(0), owner(0), ptr(NULL), data(NULL), lock(0)
    {
        //printf("got empty array2D init\n");
    }

    // creator type1
    array2D(int w, int h, unsigned int flgs = 0)
    {
        flags = flgs;
        lock = flags & ARRAY2D_LOCK_DATA;
        data = new T[h * w];
        owner = 1;
        x = w;
        y = h;
        ptr = new T*[h];

        for (int i = 0; i < h; i++) {
            ptr[i] = data + i * w;
        }

        if (flags & ARRAY2D_CLEAR_DATA) {
            memset(data, 0, w * h * sizeof(T));
        }
    }

    // creator type 2
    array2D(int w, int h, T ** source, unsigned int flgs = 0)
    {
        flags = flgs;
        //if (lock) { printf("array2D attempt to overwrite data\n");raise(SIGSEGV);}
        lock |= flags & ARRAY2D_LOCK_DATA;
        // when by reference
        // TODO: improve this code with ar_realloc()
        owner = (flags & ARRAY2D_BYREFERENCE) ? 0 : 1;

        if (owner) {
            data = new T[h * w];
        } else {
            data = NULL;
        }

        x = w;
        y = h;
        ptr = new T*[h];

        for (int i = 0; i < h; i++) {
            if (owner) {
                ptr[i] = data + i * w;

                for (int j = 0; j < w; j++) {
                    ptr[i][j] = source[i][j];
                }
            } else {
                ptr[i] = source[i];
            }
        }
    }

    array2D(const array2D&) = delete;

    // destructor
    ~array2D()
    {

        if (flags & ARRAY2D_VERBOSE) {
            printf(" deleting array2D size %dx%d \n", x, y);
        }

        if ((owner) && (data)) {
            delete[] data;
        }

        if (ptr) {
            delete[] ptr;
        }
    }

    void free()
    {
        if ((owner) && (data)) {
            delete[] data;
            data = NULL;
        }

        if (ptr) {
            delete [] ptr;
            ptr = NULL;
        }
    }

    // use with indices
    T * operator[](int index)
    {
        assert((index >= 0) && (index < y));
        return ptr[index];
    }

    // use as pointer to T**
    operator T**()
    {
        return ptr;
    }

    // use as pointer to data
    operator T*()
    {
        // only if owner this will return a valid pointer
        return data;
    }


    // useful within init of parent object
    // or use as resize of 2D array
    void operator()(int w, int h, unsigned int flgs = 0, int offset = 0)
    {
        flags = flgs;

        if (flags & ARRAY2D_VERBOSE) {
            printf("got init request %dx%d flags=%u\n", w, h, flags);
            printf("previous was data %p ptr %p \n", data, ptr);
        }

        if (lock) { // our object was locked so don't allow a change.
            printf("got init request but object was locked!\n");
            raise( SIGSEGV);
        }

        lock = flags & ARRAY2D_LOCK_DATA;

        ar_realloc(w, h, offset);

        if (flags & ARRAY2D_CLEAR_DATA) {
            memset(data + offset, 0, w * h * sizeof(T));
        }
    }

    // import from flat data
    void operator()(int w, int h, T* copy, unsigned int flgs = 0)
    {
        flags = flgs;

        if (flags & ARRAY2D_VERBOSE) {
            printf("got init request %dx%d flags=%u\n", w, h, flags);
            printf("previous was data %p ptr %p \n", data, ptr);
        }

        if (lock) { // our object was locked so don't allow a change.
            printf("got init request but object was locked!\n");
            raise( SIGSEGV);
        }

        lock = flags & ARRAY2D_LOCK_DATA;

        ar_realloc(w, h);
        memcpy(data, copy, w * h * sizeof(T));
    }
    int width()
    {
        return x;
    }
    int height()
    {
        return y;
    }

    operator bool()
    {
        return (x > 0 && y > 0);
    }

    array2D<T> & operator=( array2D<T> & rhs)
    {
        if (this != &rhs)

        {
            if (!owner) { // we can only copy same size data
                if ((x != rhs.x) || (y != rhs.y)) {
                    printf(" assignment error in array2D\n");
                    printf(" sizes differ and not owner\n");
                    raise( SIGSEGV);
                }

            } else {
                ar_realloc(rhs.x, rhs.y);
            }

            // we could have been created from a different
            // array format where each row is created by 'new'
            for (int i = 0; i < y; i++) {
                memcpy(ptr[i], rhs.ptr[i], x * sizeof(T));
            }
        }

        return *this;
    }

};
template<typename T, const size_t num>
class multi_array2D
{
private:
    array2D<T> list[num];

public:
    multi_array2D(int x, int y, int flags = 0, int offset = 0)
    {
        for (size_t i = 0; i < num; i++) {
            list[i](x, y, flags, (i + 1) * offset);
        }
    }

    ~multi_array2D()
    {
        //printf("trying to delete the list of array2D objects\n");
    }

    array2D<T> & operator[](int index)
    {
        assert(static_cast<size_t>(index) < num);
        return list[index];
    }
};
#endif /* array2D_H_ */
