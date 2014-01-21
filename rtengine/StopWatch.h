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
#define	STOPWATCH_H
#include <iostream>
#include <sys/time.h>

class StopWatch {
public:
    StopWatch( ) { stopped = false; }
    StopWatch( const char* msg) { message = msg; start(); stopped = false; }
    ~StopWatch() { if(!stopped) stop(); }
    void start()
    {
        gettimeofday(&startStruct,NULL);
    };
    void stop()
    {
        gettimeofday(&stopStruct,NULL);
        long elapsedTime = (stopStruct.tv_sec - startStruct.tv_sec) * 1000.0;      // sec to ms
        elapsedTime += (stopStruct.tv_usec - startStruct.tv_usec)/1000;
        std::cout << message << " took " << elapsedTime << "ms" <<std::endl;
        stopped = true;
    }
    void stop(const char *msg)
    {
    	message = msg;
    	stop();
    };
private:
    struct timeval startStruct;
    struct timeval stopStruct;
    const char *message;
    bool stopped;
};

#endif	/* STOPWATCH_H */
