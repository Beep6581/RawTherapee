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
 *
 *  Author: reine
 */

#ifndef STOPWATCH_H
#define STOPWATCH_H
#include <iostream>
#include "mytime.h"

#ifdef BENCHMARK
    #define BENCHFUN StopWatch StopFun(__func__);
    #define BENCHFUNMICRO StopWatch StopFun(__func__, true);
#else
    #define BENCHFUN
    #define BENCHFUNMICRO
#endif

class StopWatch
{
public:

    explicit StopWatch( const char* msg, bool microseconds = false ) : microseconds(microseconds)
    {
        message = msg;
        start();
        stopped = false;
    }
    ~StopWatch()
    {
        if(!stopped) {
            stop();
        }
    }
    void start()
    {
        startTime.set();
    };
    void stop()
    {
        stopTime.set();
        if(!microseconds) {
            long elapsedTime = stopTime.etime(startTime) / 1000;
            std::cout << message << " took " << elapsedTime << " ms" << std::endl;
        } else {
            long elapsedTime = stopTime.etime(startTime);
            std::cout << message << " took " << elapsedTime << " us" << std::endl;
        }
        stopped = true;
    }
    void stop(const char *msg)
    {
        message = msg;
        stop();
    };
private:
    bool microseconds;
    MyTime startTime;
    MyTime stopTime;
    const char *message;
    bool stopped;
};

#endif  /* STOPWATCH_H */
