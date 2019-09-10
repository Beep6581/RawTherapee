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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Author: reine
 */
#pragma once

#include <iostream>
#include <string>
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

    explicit StopWatch(const char* msg, bool microSeconds = false) : message(msg), unit(microSeconds ? " us" : " ms"), divisor(microSeconds ? 1 : 1000)
    {
        start();
        stopped = false;
    }
    ~StopWatch()
    {
        if (!stopped) {
            stop();
        }
    }
    void start()
    {
        startTime.set();
    }
    void stop()
    {
        stopTime.set();
        const long elapsedTime = stopTime.etime(startTime) / divisor;
        std::cout << message << " took " << elapsedTime << unit << std::endl;
        stopped = true;
    }

private:
    MyTime startTime;
    MyTime stopTime;
    const std::string message;
    const std::string unit;
    const int divisor;
    bool stopped;
};
