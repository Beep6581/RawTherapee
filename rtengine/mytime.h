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
#ifndef _MYTIME_
#define _MYTIME_

#ifdef WIN32
#include <windows.h>
#elif defined __APPLE__
#include <sys/time.h>
#else
#include <ctime>
#endif

class MyTime
{

public:
#ifndef WIN32
    timespec t;
#else
    LONGLONG t;
    LONGLONG baseFrequency;
    MyTime()
    {
        LARGE_INTEGER ulf;
        QueryPerformanceFrequency(&ulf);
        baseFrequency = ulf.QuadPart;
        QueryPerformanceCounter(&ulf);
        t = ulf.QuadPart;
    }
#endif

    void set ()
    {
#ifdef WIN32
        LARGE_INTEGER ulf;
        QueryPerformanceCounter(&ulf);
        t = ulf.QuadPart;
#elif defined __APPLE__
        struct timeval tv;
        gettimeofday(&tv, NULL);
        t.tv_sec = tv.tv_sec;
        t.tv_nsec = tv.tv_usec * 1000;
#else
        clock_gettime (CLOCK_REALTIME, &t);
#endif
    }

    int etime (MyTime a)
    {
#ifndef WIN32
        return (t.tv_sec - a.t.tv_sec) * 1000000 + (t.tv_nsec - a.t.tv_nsec) / 1000;
#else
        return (t - a.t) * 1000 / (baseFrequency / 1000);
#endif
    }
};


#endif
