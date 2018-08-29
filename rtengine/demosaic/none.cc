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
#include <cmath>

#include "../rawimagesource.h"
#include "../rt_math.h"
#define BENCHMARK
#include "../StopWatch.h"

using namespace std;

namespace rtengine
{

extern const Settings* settings;

void RawImageSource::nodemosaic(bool bw)
{
    red(W, H);
    green(W, H);
    blue(W, H);
    #pragma omp parallel for

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            if (bw) {
                red[i][j] = green[i][j] = blue[i][j] = rawData[i][j];
            } else if(ri->getSensorType() != ST_FUJI_XTRANS) {
                switch( FC(i, j)) {
                case 0:
                    red[i][j] = rawData[i][j];
                    green[i][j] = blue[i][j] = 0;
                    break;

                case 1:
                    green[i][j] = rawData[i][j];
                    red[i][j] = blue[i][j] = 0;
                    break;

                case 2:
                    blue[i][j] = rawData[i][j];
                    red[i][j] = green[i][j] = 0;
                    break;
                }
            } else {
                switch( ri->XTRANSFC(i, j)) {
                case 0:
                    red[i][j] = rawData[i][j];
                    green[i][j] = blue[i][j] = 0;
                    break;

                case 1:
                    green[i][j] = rawData[i][j];
                    red[i][j] = blue[i][j] = 0;
                    break;

                case 2:
                    blue[i][j] = rawData[i][j];
                    red[i][j] = green[i][j] = 0;
                    break;
                }
            }
        }
    }
}

} /* namespace */
