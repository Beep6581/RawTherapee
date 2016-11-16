////////////////////////////////////////////////////////////////
//
//  pentax pixelshift algorithm with motion detection
//
//  derived from dcrawps (https://github.com/tomtor/dcrawps), but with additional motion correction methods and adapted for RawTherapee data structures
//
//  If motion correction is enabled only the pixels which are not detected as motion are set
//  That means for a complete image you have to demosaic one of the frames with a bayer demosaicer to fill red, green and blue
//  before calling pixelshift in case motion correction is enabled.
//
//  copyright (c) Ingo Weyrich 2016
//
//
//  pixelshift.cc is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////

#include <cmath>
#include "rawimagesource.h"
#include "../rtgui/multilangmgr.h"
#include "procparams.h"
#define BENCHMARK
#include "StopWatch.h"

namespace
{

float greenDiff(float a, float b)
{
    // calculate the difference between to green samples
    // add a small epsilon to avoid division by zero
    return std::fabs(a - b) / (std::max(a, b) + 0.0001f);

}

}

using namespace std;
using namespace rtengine;

void RawImageSource::pixelshift_simple(int winx, int winy, int winw, int winh, bool detectMotion, int motion, bool showMotion, unsigned int frame, unsigned int gridSize)
{

    BENCHFUN

    if (plistener) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::pixelshift_simple]));
        plistener->setProgress(0.0);
    }

    // If the values of two corresponding green pixels differ my more then motionThreshold %, the pixel will be treated as a badGreen pixel
    const float motionThreshold = 1.f - (motion / 100.f);

    unsigned int offsX = 0, offsY = 0;
    if(detectMotion && !showMotion) {
        // if motion correction is enabled we have to adjust the offsets for the selected subframe we use for areas with motion
        switch (frame) {
            case 0:
                offsX = offsY = 0;
                break;

            case 1:
                offsX = 0;
                offsY = 1;
                break;

            case 2:
                offsX = offsY = 1;
                break;

            case 3:
                offsX = 1;
                offsY = 0;
        }
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif

    for(int i = winy + border - offsY; i < winh - (border + offsY); ++i) {
        float *greenDest = green[i + offsY];
        float *nonGreenDest0 = red[i + offsY];
        float *nonGreenDest1 = blue[i + offsY];
        int j = winx + border - offsX;
        int c = FC(i, j);

        if (c == 2 || ((c & 1) && FC(i, j + 1) == 2)) {
            // row with blue pixels => swap destination pointers for non green pixels
            std::swap(nonGreenDest0, nonGreenDest1);
        }

        // offset to keep the code short. It changes its value between 0 and 1 for each iteration of the loop
        unsigned int offset = (c & 1);

        
        float greenDifMax[gridSize];
        // motion detection checks the grid around the pixel for differences in green channels
        if(detectMotion) {
            if(gridSize == 3) {
                // compute maximum of differences for first two columns of 3x3 grid
                greenDifMax[0] =  max(greenDiff(riFrames[0 + offset]->data[i + offset][j - 1], riFrames[2 + offset]->data[i - offset + 1][j]),
                                      greenDiff(riFrames[1 - offset]->data[i - offset][j - 1], riFrames[3 - offset]->data[i + offset - 1][j]),
                                      greenDiff(riFrames[1 - offset]->data[i - offset + 2][j - 1], riFrames[3 - offset]->data[i + offset + 1][j])
                                     );
                greenDifMax[1] =  max(greenDiff(riFrames[1 - offset]->data[i - offset + 1][j], riFrames[3 - offset]->data[i + offset][j + 1]),
                                      greenDiff(riFrames[0 + offset]->data[i + offset - 1][j], riFrames[2 + offset]->data[i - offset][j + 1]),
                                      greenDiff(riFrames[0 + offset]->data[i + offset + 1][j], riFrames[2 + offset]->data[i - offset + 2][j + 1])
                                     );
            } else if(gridSize == 5) {
                // compute maximum of differences for first four columns of 5x5 grid
                greenDifMax[0] =  max(greenDiff(riFrames[1 - offset]->data[i - offset - 1][j-2], riFrames[3 - offset]->data[i + offset -2][j - 1]),
                                      greenDiff(riFrames[1 - offset]->data[i - offset + 1][j-2], riFrames[3 - offset]->data[i + offset][j - 1]),
                                      greenDiff(riFrames[1 - offset]->data[i - offset + 3][j-2], riFrames[3 - offset]->data[i + offset +2][j - 1]),
                                      greenDiff(riFrames[0 + offset]->data[i + offset - 1][j-2], riFrames[2 + offset]->data[i - offset][j - 1]),
                                      greenDiff(riFrames[0 + offset]->data[i + offset + 1][j-2], riFrames[2 + offset]->data[i - offset + 2][j - 1])
                                     );
                greenDifMax[1] =  max(greenDiff(riFrames[0 + offset]->data[i + offset-2][j - 1], riFrames[2 + offset]->data[i - offset - 1][j]),
                                      greenDiff(riFrames[0 + offset]->data[i + offset][j - 1], riFrames[2 + offset]->data[i - offset + 1][j]),
                                      greenDiff(riFrames[0 + offset]->data[i + offset+2][j - 1], riFrames[2 + offset]->data[i - offset + 3][j]),
                                      greenDiff(riFrames[1 - offset]->data[i - offset][j - 1], riFrames[3 - offset]->data[i + offset - 1][j]),
                                      greenDiff(riFrames[1 - offset]->data[i - offset + 2][j - 1], riFrames[3 - offset]->data[i + offset + 1][j])
                                     );
                greenDifMax[2] =  max(greenDiff(riFrames[1 - offset]->data[i - offset - 1][j], riFrames[3 - offset]->data[i + offset -2][j + 1]),
                                      greenDiff(riFrames[1 - offset]->data[i - offset + 1][j], riFrames[3 - offset]->data[i + offset][j + 1]),
                                      greenDiff(riFrames[1 - offset]->data[i - offset + 3][j], riFrames[3 - offset]->data[i + offset +2][j + 1]),
                                      greenDiff(riFrames[0 + offset]->data[i + offset - 1][j], riFrames[2 + offset]->data[i - offset][j + 1]),
                                      greenDiff(riFrames[0 + offset]->data[i + offset + 1][j], riFrames[2 + offset]->data[i - offset + 2][j + 1])
                                     );
                greenDifMax[3] =  max(greenDiff(riFrames[0 + offset]->data[i + offset-2][j + 1], riFrames[2 + offset]->data[i - offset - 1][j+2]),
                                      greenDiff(riFrames[0 + offset]->data[i + offset][j + 1], riFrames[2 + offset]->data[i - offset + 1][j+2]),
                                      greenDiff(riFrames[0 + offset]->data[i + offset+2][j + 1], riFrames[2 + offset]->data[i - offset + 3][j+2]),
                                      greenDiff(riFrames[1 - offset]->data[i - offset][j + 1], riFrames[3 - offset]->data[i + offset - 1][j+2]),
                                      greenDiff(riFrames[1 - offset]->data[i - offset + 2][j +- 1], riFrames[3 - offset]->data[i + offset + 1][j+2])
                                     );
            }
        }

        offset ^= 1; // 0 => 1 or 1 => 0

        // this is the index for the last column of the grid. Obviously we have to start with gridSize - 1
        int lastIndex = gridSize - 1;

        for(; j < winw - (border + offsX); ++j) {
            offset ^= 1; // 0 => 1 or 1 => 0

            if(detectMotion) {
                bool skipNext = false;
                float gridMax;
                if(gridSize == 1) {
                    // compute difference for current pixel and skip next pixel, that's the method from dcrawps
                    gridMax = greenDiff(riFrames[1 - offset]->data[i - offset + 1][j], riFrames[3 - offset]->data[i + offset][j + 1]);
                    skipNext = !showMotion;
                } else if(gridSize == 3) {
                    // compute maximum of differences for third column of 3x3 grid and save at position lastIndex
                    greenDifMax[lastIndex] = max(greenDiff(riFrames[0 + offset]->data[i + offset][j + 1], riFrames[2 + offset]->data[i - offset + 1][j + 2]),
                                                 greenDiff(riFrames[1 - offset]->data[i - offset][j + 1], riFrames[3 - offset]->data[i + offset - 1][j + 2]),
                                                 greenDiff(riFrames[1 - offset]->data[i - offset + 2][j + 1], riFrames[3 - offset]->data[i + offset + 1][j + 2])
                                                );
                    gridMax = max(greenDifMax[0],greenDifMax[1],greenDifMax[2]);
                } else if(gridSize == 5) {
                    // compute maximum of differences for fifth column of 5x5 grid and save at position lastIndex
                    greenDifMax[lastIndex] = max(greenDiff(riFrames[1 - offset]->data[i - offset - 1][j+2], riFrames[3 - offset]->data[i + offset -2][j + 3]),
                                                 greenDiff(riFrames[1 - offset]->data[i - offset + 1][j+2], riFrames[3 - offset]->data[i + offset][j + 3]),
                                                 greenDiff(riFrames[1 - offset]->data[i - offset + 3][j+2], riFrames[3 - offset]->data[i + offset +2][j + 3]),
                                                 greenDiff(riFrames[0 + offset]->data[i + offset - 1][j+2], riFrames[2 + offset]->data[i - offset][j + 3]),
                                                 greenDiff(riFrames[0 + offset]->data[i + offset + 1][j+2], riFrames[2 + offset]->data[i - offset + 2][j + 3])
                                                );
                    gridMax = max(greenDifMax[0],greenDifMax[1],greenDifMax[2],greenDifMax[3],greenDifMax[4]);
                }
                // adjust index for next column
                lastIndex ++;
                lastIndex = lastIndex == gridSize ? 0 : lastIndex;

                if (gridMax > motionThreshold) {
                    // at least one of the tested pixels of the grid is detected as motion
                    if(showMotion) {
                        // if showMotion is enabled make the pixel green
                        greenDest[j] = 10000.f;
                        nonGreenDest0[j] = nonGreenDest1[j] = 0.f;
                    }
                    if(skipNext) {
                        // treat the horizontally next pixel also as motion
                        j++;
                        offset ^= 1;
                    }
                    // do not set the motion pixel values. They have already been set by demosaicer or showMotion
                    continue;
                }
            }

            // motion correction disabled or no motion detected => combine the values from the four pixelshift frames
            greenDest[j + offsX] = (riFrames[1 - offset]->data[i - offset + 1][j] + riFrames[3 - offset]->data[i + offset][j + 1]) / 2.f;
            nonGreenDest0[j + offsX] = riFrames[(offset << 1) + offset]->data[i][j + offset];
            nonGreenDest1[j + offsX] = riFrames[2 - offset]->data[i + 1][j - offset + 1];
        }
    }

    if(plistener) {
        plistener->setProgress(1.0);
    }
}
