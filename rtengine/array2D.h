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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
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
 *      access to elements is as simple as:
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
 *
 */
#pragma once

#include <cassert>
#include <cstring>
#include <sys/types.h>
#include <vector>
#include "noncopyable.h"

// flags for use
constexpr unsigned int ARRAY2D_CLEAR_DATA = 1;
constexpr unsigned int ARRAY2D_BYREFERENCE = 2;


template<typename T>
class array2D
{

private:
    ssize_t width;
    std::vector<T*> rows;
    std::vector<T> buffer;

    void initRows(ssize_t h, int offset = 0)
    {
        rows.resize(h);
        T* start = buffer.data() + offset;
        for (ssize_t i = 0; i < h; ++i) {
            rows[i] = start + width * i;
        }
    }

    void ar_realloc(ssize_t w, ssize_t h, int offset = 0)
    {
        width = w;
        buffer.resize(h * width + offset);
        initRows(h, offset);
    }
public:

    // use as empty declaration, resize before use!
    // very useful as a member object
    array2D() : width(0) {}

    // creator type1
    array2D(int w, int h, unsigned int flags = 0) : width(w)
    {
        if (flags & ARRAY2D_CLEAR_DATA) {
            buffer.resize(h * width, 0);
        } else {
            buffer.resize(h * width);
        }
        initRows(h);
    }

    // creator type 2
    array2D(int w, int h, T ** source, unsigned int flags = 0) : width(w)
    {
        rows.resize(h);
        if (!(flags & ARRAY2D_BYREFERENCE)) {
            buffer.resize(h * width);
            T* start = buffer.data();
            for (ssize_t i = 0; i < h; ++i) {
                rows[i] = start + i * width;
                for (ssize_t j = 0; j < width; ++j) {
                    rows[i][j] = source[i][j];
                }
            }
        } else {
            for (ssize_t i = 0; i < h; ++i) {
                rows[i] = source[i];
            }
        }
    }

    // creator type 3
    array2D(int w, int h, int startx, int starty, T ** source, unsigned int flags = 0) : width(w)
    {
        rows.resize(h);
        if (!(flags & ARRAY2D_BYREFERENCE)) {
            buffer.resize(h * width);
            T* start = buffer.data();
            for (ssize_t i = 0; i < h; ++i) {
                rows[i] = start + i * width;
                for (ssize_t j = 0; j < width; ++j) {
                    rows[i][j] = source[i + starty][j + startx];
                }
            }
        } else {
            for (ssize_t i = 0; i < h; ++i) {
                rows[i] = source[i + starty] + startx;
            }
        }
    }

    array2D(const array2D& other) :
        width(other.width),
        buffer(other.buffer)
    {
        initRows(other.rows.size());
    }

    array2D& operator =(const array2D& other)
    {
        if (this != &other) {
            free();
            width = other.width;
            buffer = other.buffer;
            initRows(other.rows.size());
        }

        return *this;
    }

    void fill(const T val, bool multiThread = false)
    {
        const ssize_t height = rows.size();
#ifdef _OPENMP
        #pragma omp parallel for if(multiThread)
#endif
        for (ssize_t i = 0; i < width * height; ++i) {
            buffer[i] = val;
        }
    }

    void free()
    {
        buffer.clear();
        rows.clear();
        width = 0;
    }

    // use with indices
    T * operator[](int index)
    {
        assert((index >= 0) && (std::size_t(index) < rows.size()));
        return rows[index];
    }

    const T * operator[](int index) const
    {
        assert((index >= 0) && (std::size_t(index) < rows.size()));
        return rows[index];
    }

    // use as pointer to T**
    operator T**()
    {
        return rows.data();
    }

    // use as pointer to T**
    operator const T* const *() const
    {
        return rows.data();
    }

    // use as pointer to buffer
    operator T*()
    {
        // only if owner this will return a valid pointer
        return buffer.data();
    }

    operator const T*() const
    {
        // only if owner this will return a valid pointer
        return buffer.data();
    }


    // useful within init of parent object
    // or use as resize of 2D array
    void operator()(int w, int h, unsigned int flags = 0, int offset = 0)
    {
        ar_realloc(w, h, offset);

        if (flags & ARRAY2D_CLEAR_DATA) {
            fill(0);
        }
    }

    array2D<T>& operator+=(const array2D<T>& rhs)
    {
        if (rhs.getWidth() == this->getWidth() && rhs.getHeight() == this->getHeight()) {
            for (int i = 0; i < getHeight(); ++i) {
#ifdef _OPENMP
                #pragma omp simd
#endif

                for (int j = 0; j < getWidth(); ++j) {
                    rows[i][j] += rhs[i][j];
                }
            }
        }

        return *this;
    }

    // import from flat data
    void operator()(std::size_t w, std::size_t h, const T* const copy)
    {
        ar_realloc(w, h);
        for (std::size_t y = 0; y < h; ++y) {
            std::copy(copy + y * w, copy + y * w + w, rows.data()[y]);
        }
    }

    int getWidth() const
    {
        return width;
    }
    int getHeight() const
    {
        return rows.size();
    }

    operator bool()
    {
        return (width > 0 && !rows.empty());
    }

};
template<typename T, const size_t num>
class multi_array2D : public rtengine::NonCopyable
{
private:
    array2D<T> list[num];

public:
    multi_array2D(int width, int height, int flags = 0, int offset = 0)
    {
        for (size_t i = 0; i < num; ++i) {
            list[i](width, height, flags, (i + 1) * offset);
        }
    }

    array2D<T> & operator[](int index)
    {
        assert(static_cast<size_t>(index) < num);
        return list[index];
    }
};
