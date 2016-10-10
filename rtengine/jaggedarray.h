/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2016 Ingo Weyrich <heckflosse67@gmx.de>
 *  Copyright (c) 2016 Adam Reichold <adam.reichold@t-online.de>
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
#pragma once

#include "noncopyable.h"

namespace rtengine
{

// These emulate a jagged array, but use only 2 allocations instead of 1 + H.

template<class T>
inline T** const allocJaggedArray (const int W, const int H, const bool initZero = false)
{
    T** const a = new T*[H];
    a[0] = new T[H * W];

    for (int i = 1; i < H; ++i) {
        a[i] = a[i - 1] + W;
    }

    if (initZero) {
        std::memset(a[0], 0, sizeof(T) * W * H);
    }

    return a;
}

template<class T>
inline void freeJaggedArray (T** const a)
{
    delete [] a[0];
    delete [] a;
}

template<class T>
class JaggedArray :
    public NonCopyable
{
public:
    JaggedArray (const int W, const int H, const bool initZero = false)
    {
        a = allocJaggedArray<T> (W, H, initZero);
    }
    ~JaggedArray ()
    {
        if (a) {
            freeJaggedArray<T> (a);
            a = nullptr;
        }
    }

    operator T** const () const
    {
        return a;
    }

private:
    T** a;

};

// Declared but not defined to prevent
// explicitly freeing a JaggedArray<T> implicitly cast to T**.
template<class T>
void freeJaggedArray (JaggedArray<T>&);

} // rtengine
