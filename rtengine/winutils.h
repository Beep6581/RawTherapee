/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2021 Lawrence Lee
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
#pragma once

#ifdef _WIN32

#include <aclapi.h>
#include <windows.h>

#include "noncopyable.h"


/**
 * Wrapper for pointers to memory allocated by HeapAlloc.
 *
 * Memory is automatically freed when the object goes out of scope.
 */
template <typename T>
class WinHeapPtr : public rtengine::NonCopyable
{
private:
    const T ptr;

public:
    WinHeapPtr() = delete;

    /** Allocates the specified number of bytes in the process heap. */
    explicit WinHeapPtr(SIZE_T bytes): ptr(static_cast<T>(HeapAlloc(GetProcessHeap(), 0, bytes))) {};

    ~WinHeapPtr()
    {
        // HeapFree does a null check.
        HeapFree(GetProcessHeap(), 0, static_cast<LPVOID>(ptr));
    }

    T operator ->() const
    {
        return ptr;
    }

    operator T() const
    {
        return ptr;
    }
};

/**
 * Wrapper for HLOCAL pointers to memory allocated by LocalAlloc.
 *
 * Memory is automatically freed when the object goes out of scope.
 */
template <typename T>
class WinLocalPtr : public rtengine::NonCopyable
{
private:
    const T ptr;

public:
    WinLocalPtr() = delete;

    /** Wraps a raw pointer. */
    WinLocalPtr(T pointer): ptr(pointer) {};

    ~WinLocalPtr()
    {
        // LocalFree does a null check.
        LocalFree(static_cast<HLOCAL>(ptr));
    }

    T operator ->() const
    {
        return ptr;
    }

    operator T() const
    {
        return ptr;
    }
};

/**
 * Wrapper for HANDLEs.
 *
 * Handles are automatically closed when the object goes out of scope.
 */
class WinHandle : public rtengine::NonCopyable
{
private:
    const HANDLE handle;

public:
    WinHandle() = delete;

    /** Wraps a HANDLE. */
    WinHandle(HANDLE handle): handle(handle) {};

    ~WinHandle()
    {
        CloseHandle(handle);
    }

    operator HANDLE() const
    {
        return handle;
    }
};

#endif
