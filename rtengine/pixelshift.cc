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

float greenDiff(float a, float b, bool adaptive, float scale, float stddevFactor, float eperIso, float nreadIso, float prnu)
{
    // calculate the difference between to green samples
    // add a small epsilon to avoid division by zero
    float diff = std::fabs(a - b) / (std::max(a, b) + 0.01f);
    if(adaptive) {
        float avg = (a+b)/2.f;
        avg *= scale;
        float stddev = sqrtf(avg * eperIso + nreadIso * nreadIso + prnu * prnu);
        float korr = stddevFactor * stddev / (a * scale);
        diff -= korr;
    }
    return diff;
}

}

using namespace std;
using namespace rtengine;

void RawImageSource::pixelshift_simple(int winx, int winy, int winw, int winh, bool detectMotion, int motion, bool showMotion, unsigned int frame, unsigned int gridSize, bool adaptive, float stddevFactor, float eperIso, float nreadIso, float prnu)
{

    BENCHFUN

    if (plistener) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::pixelshift_simple]));
        plistener->setProgress(0.0);
    }


    gridSize += ((gridSize & 1) == 0 ? 1 : 0);
    // Lookup table for non adaptive (slider) mode
    LUTf log2Lut(32768, LUT_CLIP_BELOW | LUT_CLIP_ABOVE);
    log2Lut[0] = 0;
    const float lutStrength = 2.f;
    const float scaleGreen = 1.f / scale_mul[1];
    for(int i=2; i < 65536; i+=2)
        log2Lut[i>>1] = 2.f * log2(i) / 100.f;

//    const float eperIso = 0.75f * idata->getISOSpeed() / 100;
    eperIso *= (idata->getISOSpeed() / 100);
    nreadIso *= (idata->getISOSpeed() / 100);
//    const float nreadIso = 5.f * idata->getISOSpeed() / 100;
//    const float prnu = 1.f;
//    const float stddevFactor = 4.f;
    // If the values of two corresponding green pixels differ my more then motionThreshold %, the pixel will be treated as a badGreen pixel
    float motionThreshold = 1.f - (motion / 100.f);
    // For shades of green motion indicators 
    const float blendFactor = (motion == 0.f ? 1.f : 1.f / (1.f - motionThreshold));

//    bool checkRedBlue = (gridSize == 5);
//    bool checkRedBlue = false;
    unsigned int offsX = 0, offsY = 0;
    // We have to adjust the offsets for the selected subframe we use for areas with motion
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
        if(detectMotion || adaptive) {
            if(gridSize == 3) {
                // compute maximum of differences for first two columns of 3x3 grid
                greenDifMax[0] =  max(greenDiff(riFrames[0 + offset]->data[i + offset][j - 1], riFrames[2 + offset]->data[i - offset + 1][j], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[1 - offset]->data[i - offset][j - 1], riFrames[3 - offset]->data[i + offset - 1][j], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[1 - offset]->data[i - offset + 2][j - 1], riFrames[3 - offset]->data[i + offset + 1][j], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu)
                                     );
                greenDifMax[1] =  max(greenDiff(riFrames[1 - offset]->data[i - offset + 1][j], riFrames[3 - offset]->data[i + offset][j + 1], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[0 + offset]->data[i + offset - 1][j], riFrames[2 + offset]->data[i - offset][j + 1], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[0 + offset]->data[i + offset + 1][j], riFrames[2 + offset]->data[i - offset + 2][j + 1], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu)
                                     );
            } else if(gridSize == 5) {
                // compute maximum of differences for first four columns of 5x5 grid
                greenDifMax[0] =  max(greenDiff(riFrames[1 - offset]->data[i - offset - 1][j-2], riFrames[3 - offset]->data[i + offset -2][j - 1], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[1 - offset]->data[i - offset + 1][j-2], riFrames[3 - offset]->data[i + offset][j - 1], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[1 - offset]->data[i - offset + 3][j-2], riFrames[3 - offset]->data[i + offset +2][j - 1], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[0 + offset]->data[i + offset - 1][j-2], riFrames[2 + offset]->data[i - offset][j - 1], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[0 + offset]->data[i + offset + 1][j-2], riFrames[2 + offset]->data[i - offset + 2][j - 1], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu)
                                     );
                greenDifMax[1] =  max(greenDiff(riFrames[0 + offset]->data[i + offset-2][j - 1], riFrames[2 + offset]->data[i - offset - 1][j], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[0 + offset]->data[i + offset][j - 1], riFrames[2 + offset]->data[i - offset + 1][j], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[0 + offset]->data[i + offset+2][j - 1], riFrames[2 + offset]->data[i - offset + 3][j], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[1 - offset]->data[i - offset][j - 1], riFrames[3 - offset]->data[i + offset - 1][j], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[1 - offset]->data[i - offset + 2][j - 1], riFrames[3 - offset]->data[i + offset + 1][j], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu)
                                     );
                greenDifMax[2] =  max(greenDiff(riFrames[1 - offset]->data[i - offset - 1][j], riFrames[3 - offset]->data[i + offset -2][j + 1], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[1 - offset]->data[i - offset + 1][j], riFrames[3 - offset]->data[i + offset][j + 1], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[1 - offset]->data[i - offset + 3][j], riFrames[3 - offset]->data[i + offset +2][j + 1], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[0 + offset]->data[i + offset - 1][j], riFrames[2 + offset]->data[i - offset][j + 1], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[0 + offset]->data[i + offset + 1][j], riFrames[2 + offset]->data[i - offset + 2][j + 1], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu)
                                     );
                greenDifMax[3] =  max(greenDiff(riFrames[0 + offset]->data[i + offset-2][j + 1], riFrames[2 + offset]->data[i - offset - 1][j+2], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[0 + offset]->data[i + offset][j + 1], riFrames[2 + offset]->data[i - offset + 1][j+2], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[0 + offset]->data[i + offset+2][j + 1], riFrames[2 + offset]->data[i - offset + 3][j+2], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[1 - offset]->data[i - offset][j + 1], riFrames[3 - offset]->data[i + offset - 1][j+2], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                      greenDiff(riFrames[1 - offset]->data[i - offset + 2][j +- 1], riFrames[3 - offset]->data[i + offset + 1][j+2], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu)
                                     );
            }
        }

        offset ^= 1; // 0 => 1 or 1 => 0

        // this is the index for the last column of the grid. Obviously we have to start with gridSize - 1
        int lastIndex = gridSize - 1;

        for(; j < winw - (border + offsX); ++j) {
            offset ^= 1; // 0 => 1 or 1 => 0

            if(detectMotion || adaptive) {
                bool skipNext = false;
                float gridMax;
                if(gridSize == 1) {
                    // compute difference for current pixel and skip next pixel, that's the method from dcrawps
                    gridMax = greenDiff(riFrames[1 - offset]->data[i - offset + 1][j], riFrames[3 - offset]->data[i + offset][j + 1], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu);
                    skipNext = !showMotion;
                } else if(gridSize == 3) {
                    // compute maximum of differences for third column of 3x3 grid and save at position lastIndex
                    greenDifMax[lastIndex] = max(greenDiff(riFrames[0 + offset]->data[i + offset][j + 1], riFrames[2 + offset]->data[i - offset + 1][j + 2], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                                 greenDiff(riFrames[1 - offset]->data[i - offset][j + 1], riFrames[3 - offset]->data[i + offset - 1][j + 2], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                                 greenDiff(riFrames[1 - offset]->data[i - offset + 2][j + 1], riFrames[3 - offset]->data[i + offset + 1][j + 2], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu)
                                                );
                    gridMax = max(greenDifMax[0],greenDifMax[1],greenDifMax[2]);
                } else if(gridSize == 5) {
                    // compute maximum of differences for fifth column of 5x5 grid and save at position lastIndex
                    greenDifMax[lastIndex] = max(greenDiff(riFrames[1 - offset]->data[i - offset - 1][j+2], riFrames[3 - offset]->data[i + offset -2][j + 3], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                                 greenDiff(riFrames[1 - offset]->data[i - offset + 1][j+2], riFrames[3 - offset]->data[i + offset][j + 3], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                                 greenDiff(riFrames[1 - offset]->data[i - offset + 3][j+2], riFrames[3 - offset]->data[i + offset +2][j + 3], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                                 greenDiff(riFrames[0 + offset]->data[i + offset - 1][j+2], riFrames[2 + offset]->data[i - offset][j + 3], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu),
                                                 greenDiff(riFrames[0 + offset]->data[i + offset + 1][j+2], riFrames[2 + offset]->data[i - offset + 2][j + 3], adaptive, scaleGreen, stddevFactor, eperIso, nreadIso, prnu)
                                                );
                    gridMax = max(greenDifMax[0],greenDifMax[1],greenDifMax[2],greenDifMax[3],greenDifMax[4]);
                }
                // adjust index for next column
                lastIndex ++;
                lastIndex = lastIndex == gridSize ? 0 : lastIndex;

                // increase motion detection dependent on brightness
                float korr;
                float thresh;
                if(adaptive) {
                    korr = 0.f;
                    thresh = 0;
                } else {
                    korr = log2Lut[((int)(riFrames[1 - offset]->data[i - offset + 1][j] * scaleGreen))>>1];
                    thresh = motionThreshold;
                }

                if (gridMax > thresh - korr) {
                    float blend = (gridMax - thresh + korr) * blendFactor;
                    // at least one of the tested pixels of the grid is detected as motion
                    if(showMotion) {
                        // if showMotion is enabled make the pixel green
                        greenDest[j + offsX] = 1000.f + 25000.f * blend;
                        nonGreenDest0[j + offsX] = nonGreenDest1[j + offsX] = 0.f;
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

//            if(false && detectMotion && checkRedBlue) {
//                float ng1 = riFrames[(offset << 1) + offset]->data[i][j + offset];
//                float ng0 = riFrames[((offset^1) << 1) + (offset^1)]->data[i][j + (offset^1)+1];
//                float ng2 = riFrames[((offset^1) << 1) + (offset^1)]->data[i][j + (offset^1)-1];
//                float diff0 = ng1 - ng0;
//                float diff2 = ng1 - ng2;
//                    float gridMax;
//                if(diff0 * diff2 > 0.f) {
////                    if(greenDiff(ng1, fabsf(diff0) < fabsf(diff2) ? ng0 : ng2) > motionThreshold ) {
//                    gridMax = greenDiff(ng1, std::max(ng0, ng2));
////                    gridMax = greenDiff(ng1, ((ng0 + ng2) / 2.f));
//                    if(gridMax > motionThreshold ) {
//                    float factor = 1.f / (1.f - motionThreshold);
//                    float blend = (gridMax - motionThreshold) * factor;
//                    if(showMotion) {
//                            // if showMotion is enabled make the pixel green
//                            greenDest[j + offsX] = nonGreenDest1[j + offsX] = 0.f;
//                            nonGreenDest0[j + offsX] = 20000.f;
//                            continue;
////                            greenDest[j + offsX+1] = nonGreenDest1[j + offsX+1] = 0.f;
////                            nonGreenDest0[j + offsX+1] = 20000.f;
//                        }
//                        greenDest[j + offsX] = (riFrames[1 - offset]->data[i - offset + 1][j] + riFrames[3 - offset]->data[i + offset][j + 1]) / 2.f;
//
////                        greenDest[j + offsX] = intp(blend, greenDest[j + offsX],(riFrames[1 - offset]->data[i - offset + 1][j] + riFrames[3 - offset]->data[i + offset][j + 1]) / 2.f);
////                        nonGreenDest0[j + offsX] = (ng1 + (diff0 < diff2 ? ng0 : ng2)) / 2.f;
//                        nonGreenDest0[j + offsX] = intp(blend, nonGreenDest0[j + offsX], riFrames[(offset << 1) + offset]->data[i][j + offset]);
////                        nonGreenDest0[j + offsX] = intp(blend, nonGreenDest0[j + offsX], riFrames[(offset << 1) + offset]->data[i][j + offset]);
//                        nonGreenDest1[j + offsX] = riFrames[2 - offset]->data[i + 1][j - offset + 1];
//
////                        nonGreenDest1[j + offsX] = intp(blend, nonGreenDest1[j + offsX], riFrames[2 - offset]->data[i + 1][j - offset + 1]);
//
////                        if(skipNext) {
////                            // treat the horizontally next pixel also as motion
////                            j++;
////                            offset ^= 1;
////                        }
//                        // do not set the motion pixel values. They have already been set by demosaicer or showMotion
//                        continue;
//                    }
//                }
//                ng1 = riFrames[2 - offset]->data[i + 1][j - offset + 1];
//                ng0 = riFrames[2 - (offset^1)]->data[i + 1][j - (offset^1) + 2];
//                ng2 = riFrames[2 - (offset^1)]->data[i + 1][j - (offset^1)];
//                diff0 = ng1 - ng0;
//                diff2 = ng1 - ng2;
//                if(signbit(diff0) == signbit(diff2)) {
////                    if(greenDiff(ng1, fabsf(diff0) < fabsf(diff2) ? ng0 : ng2) > motionThreshold ) {
//                    gridMax = greenDiff(ng1, std::max(ng0, ng2));
////                    gridMax = greenDiff(ng1, ((ng0 + ng2) / 2.f));
//                    if(gridMax > motionThreshold ) {
//                    float factor = 1.f / (1.f - motionThreshold);
//                    float blend = (gridMax - motionThreshold) * factor;
//                    if(showMotion) {
//                            // if showMotion is enabled make the pixel green
//                            greenDest[j + offsX] = nonGreenDest0[j + offsX] = 0.f;
//                            nonGreenDest1[j + offsX] = 20000.f;
////                            greenDest[j + offsX+1] = nonGreenDest0[j + offsX+1] = 0.f;
////                            nonGreenDest1[j + offsX+1] = 20000.f;
//continue;
//                        }
//                        greenDest[j + offsX] = (riFrames[1 - offset]->data[i - offset + 1][j] + riFrames[3 - offset]->data[i + offset][j + 1]) / 2.f;
////                                            greenDest[j + offsX] = intp(blend, greenDest[j + offsX],(riFrames[1 - offset]->data[i - offset + 1][j] + riFrames[3 - offset]->data[i + offset][j + 1]) / 2.f);
////                        nonGreenDest0[j + offsX] = intp(blend, nonGreenDest0[j + offsX], riFrames[(offset << 1) + offset]->data[i][j + offset]);
//                        nonGreenDest0[j + offsX] = riFrames[(offset << 1) + offset]->data[i][j + offset];
//
////                        nonGreenDest1[j + offsX] = (ng1 + (diff0 < diff2 ? ng0 : ng2)) / 2.f;
//                        nonGreenDest1[j + offsX] = intp(blend, nonGreenDest1[j + offsX], riFrames[2 - offset]->data[i + 1][j - offset + 1]);
////                        if(skipNext) {
////                            // treat the horizontally next pixel also as motion
////                            j++;
////                            offset ^= 1;
////                        }
//                        // do not set the motion pixel values. They have already been set by demosaicer or showMotion
//                        continue;
//                    }
//                }
//            }
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
