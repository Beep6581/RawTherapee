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

float colourDiff(float a, float b, bool adaptive, float stddevFactor, float eperIso, float nreadIso, float prnu, bool showMotion)
{
    // calculate the difference between to green samples
    if(adaptive) {
        float gDiff = a - b;
        gDiff *= eperIso;
        gDiff *= gDiff;
        float avg = (a + b) / 2.f;
        avg *= eperIso;
        prnu *= avg;
        float stddev = stddevFactor * (avg + nreadIso + prnu * prnu);
        float result = gDiff - stddev;
        if(!showMotion) {
            return result;
        } else if(result > 0.f) { // for the motion mask
            return std::fabs(a - b) / std::max(a, b) + 0.01f;
        } else {
            return 0.f;
        }
    } else {
        float gDiff = std::fabs(a - b);
        // add a small epsilon to avoid division by zero
        float maxVal = std::max(a, b) + 0.01f;
        return gDiff / maxVal;
    }
}

}

using namespace std;
using namespace rtengine;

void RawImageSource::pixelshift(int winx, int winy, int winw, int winh, bool detectMotion, int motion, bool showMotion, bool showOnlyMask, unsigned int frame, unsigned int gridSize, bool adaptive, float stddevFactor, float eperIso, float nreadIso, float prnu, bool checkNonGreenHorizontal)
{

    BENCHFUN

    if (plistener) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::pixelshift_simple]));
        plistener->setProgress(0.0);
    }


    printf("Pixelshift parameters : gridSize %d\tadaptive %d\tstdDevFactor %f\telectrons %f\tnread %f\tprnu %f\n",gridSize, adaptive, stddevFactor, eperIso, nreadIso, prnu);
    gridSize += ((gridSize & 1) == 0 ? 1 : 0);
    // Lookup table for non adaptive (slider) mode
    LUTf log2Lut(32768, LUT_CLIP_BELOW | LUT_CLIP_ABOVE);
    
    if(detectMotion && !adaptive) {
        const float lutStrength = 2.f;
        log2Lut[0] = 0;
        for(int i=2; i < 65536; i+=2)
            log2Lut[i>>1] = lutStrength * log2(i) / 100.f;
    }
    const float scaleGreen = 1.f / scale_mul[1];

    eperIso *= (100.f / idata->getISOSpeed());
    
    float eperIsoGreen = eperIso * scaleGreen;

    prnu /= 100.f;
    stddevFactor *= stddevFactor;
    nreadIso *= nreadIso;
    
    // If the values of two corresponding green pixels differ my more then motionThreshold %, the pixel will be treated as a badGreen pixel
    float motionThreshold = 1.f - (motion / 100.f);
    // For shades of green motion indicators 
    const float blendFactor = (motion == 0.f ? 1.f : 1.f / (1.f - motionThreshold));

    bool checkNonGreen = true;
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

    const float thresh = adaptive ? 0.f : motionThreshold;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif

    for(int i = winy + border - offsY; i < winh - (border + offsY); ++i) {
        float *greenDest = green[i + offsY];
        float *nonGreenDest0 = red[i + offsY];
        float *nonGreenDest1 = blue[i + offsY];
        int j = winx + border - offsX;
        int c = FC(i, j);
        float scaleNonGreen0 = 1.f / scale_mul[0];
        float scaleNonGreen2 = 1.f / scale_mul[2];
        float eperIsoNonGreen0 = eperIso / scale_mul[0];
        float eperIsoNonGreen2 = eperIso / scale_mul[2];

        if (c == 2 || ((c & 1) && FC(i, j + 1) == 2)) {
            // row with blue pixels => swap destination pointers for non green pixels
            std::swap(nonGreenDest0, nonGreenDest1);
            std::swap(scaleNonGreen0, scaleNonGreen2);
            std::swap(eperIsoNonGreen0, eperIsoNonGreen2);
        }

        // offset to keep the code short. It changes its value between 0 and 1 for each iteration of the loop
        unsigned int offset = (c & 1);

        float greenDifMax[gridSize];
        // motion detection checks the grid around the pixel for differences in green channels
        if(detectMotion || adaptive) {
            if(gridSize == 3) {
                // compute maximum of differences for first two columns of 3x3 grid
                greenDifMax[0] =  max(colourDiff(riFrames[0 + offset]->data[i + offset][j - 1], riFrames[2 + offset]->data[i - offset + 1][j], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[1 - offset]->data[i - offset][j - 1], riFrames[3 - offset]->data[i + offset - 1][j], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[1 - offset]->data[i - offset + 2][j - 1], riFrames[3 - offset]->data[i + offset + 1][j], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion)
                                     );
                greenDifMax[1] =  max(colourDiff(riFrames[1 - offset]->data[i - offset + 1][j], riFrames[3 - offset]->data[i + offset][j + 1], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[0 + offset]->data[i + offset - 1][j], riFrames[2 + offset]->data[i - offset][j + 1], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[0 + offset]->data[i + offset + 1][j], riFrames[2 + offset]->data[i - offset + 2][j + 1], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion)
                                     );
            } else if(gridSize == 5) {
                // compute maximum of differences for first four columns of 5x5 grid
                greenDifMax[0] =  max(colourDiff(riFrames[1 - offset]->data[i - offset - 1][j-2], riFrames[3 - offset]->data[i + offset -2][j - 1], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[1 - offset]->data[i - offset + 1][j-2], riFrames[3 - offset]->data[i + offset][j - 1], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[1 - offset]->data[i - offset + 3][j-2], riFrames[3 - offset]->data[i + offset +2][j - 1], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[0 + offset]->data[i + offset - 1][j-2], riFrames[2 + offset]->data[i - offset][j - 1], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[0 + offset]->data[i + offset + 1][j-2], riFrames[2 + offset]->data[i - offset + 2][j - 1], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion)
                                     );
                greenDifMax[1] =  max(colourDiff(riFrames[0 + offset]->data[i + offset-2][j - 1], riFrames[2 + offset]->data[i - offset - 1][j], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[0 + offset]->data[i + offset][j - 1], riFrames[2 + offset]->data[i - offset + 1][j], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[0 + offset]->data[i + offset+2][j - 1], riFrames[2 + offset]->data[i - offset + 3][j], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[1 - offset]->data[i - offset][j - 1], riFrames[3 - offset]->data[i + offset - 1][j], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[1 - offset]->data[i - offset + 2][j - 1], riFrames[3 - offset]->data[i + offset + 1][j], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion)
                                     );
                greenDifMax[2] =  max(colourDiff(riFrames[1 - offset]->data[i - offset - 1][j], riFrames[3 - offset]->data[i + offset -2][j + 1], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[1 - offset]->data[i - offset + 1][j], riFrames[3 - offset]->data[i + offset][j + 1], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[1 - offset]->data[i - offset + 3][j], riFrames[3 - offset]->data[i + offset +2][j + 1], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[0 + offset]->data[i + offset - 1][j], riFrames[2 + offset]->data[i - offset][j + 1], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[0 + offset]->data[i + offset + 1][j], riFrames[2 + offset]->data[i - offset + 2][j + 1], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion)
                                     );
                greenDifMax[3] =  max(colourDiff(riFrames[0 + offset]->data[i + offset-2][j + 1], riFrames[2 + offset]->data[i - offset - 1][j+2], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[0 + offset]->data[i + offset][j + 1], riFrames[2 + offset]->data[i - offset + 1][j+2], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[0 + offset]->data[i + offset+2][j + 1], riFrames[2 + offset]->data[i - offset + 3][j+2], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[1 - offset]->data[i - offset][j + 1], riFrames[3 - offset]->data[i + offset - 1][j+2], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                      colourDiff(riFrames[1 - offset]->data[i - offset + 2][j +- 1], riFrames[3 - offset]->data[i + offset + 1][j+2], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion)
                                     );
            }
        }

        offset ^= 1; // 0 => 1 or 1 => 0

        // this is the index for the last column of the grid. Obviously we have to start with gridSize - 1
        int lastIndex = gridSize - 1;
        float korr = 0.f;

        for(; j < winw - (border + offsX); ++j) {
            offset ^= 1; // 0 => 1 or 1 => 0

            if(detectMotion || adaptive) {
                bool skipNext = false;
                float gridMax;
                if(gridSize == 1) {
                    // compute difference for current pixel and skip next pixel, that's the method from dcrawps
                    gridMax = colourDiff(riFrames[1 - offset]->data[i - offset + 1][j], riFrames[3 - offset]->data[i + offset][j + 1], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion);
                    skipNext = !showMotion;
                } else if(gridSize == 3) {
                    // compute maximum of differences for third column of 3x3 grid and save at position lastIndex
                    greenDifMax[lastIndex] = max(colourDiff(riFrames[0 + offset]->data[i + offset][j + 1], riFrames[2 + offset]->data[i - offset + 1][j + 2], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                                 colourDiff(riFrames[1 - offset]->data[i - offset][j + 1], riFrames[3 - offset]->data[i + offset - 1][j + 2], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                                 colourDiff(riFrames[1 - offset]->data[i - offset + 2][j + 1], riFrames[3 - offset]->data[i + offset + 1][j + 2], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion)
                                                );
                    gridMax = max(greenDifMax[0],greenDifMax[1],greenDifMax[2]);
                } else if(gridSize == 5) {
                    // compute maximum of differences for fifth column of 5x5 grid and save at position lastIndex
                    greenDifMax[lastIndex] = max(colourDiff(riFrames[1 - offset]->data[i - offset - 1][j+2], riFrames[3 - offset]->data[i + offset -2][j + 3], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                                 colourDiff(riFrames[1 - offset]->data[i - offset + 1][j+2], riFrames[3 - offset]->data[i + offset][j + 3], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                                 colourDiff(riFrames[1 - offset]->data[i - offset + 3][j+2], riFrames[3 - offset]->data[i + offset +2][j + 3], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                                 colourDiff(riFrames[0 + offset]->data[i + offset - 1][j+2], riFrames[2 + offset]->data[i - offset][j + 3], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion),
                                                 colourDiff(riFrames[0 + offset]->data[i + offset + 1][j+2], riFrames[2 + offset]->data[i - offset + 2][j + 3], adaptive, stddevFactor, eperIsoGreen, nreadIso, prnu, showMotion)
                                                );
                    gridMax = max(greenDifMax[0],greenDifMax[1],greenDifMax[2],greenDifMax[3],greenDifMax[4]);
                }
                // adjust index for next column
                lastIndex ++;
                lastIndex = lastIndex == gridSize ? 0 : lastIndex;

                // increase motion detection dependent on brightness
                if(!adaptive) {
                    korr = log2Lut[((int)(riFrames[1 - offset]->data[i - offset + 1][j] * scaleGreen))>>1];
                }

                if (gridMax > thresh - korr) {
                    // at least one of the tested pixels of the grid is detected as motion
                    if(showMotion) {
                        float blend = (gridMax - thresh + korr) * blendFactor;
                        if(!showOnlyMask) {
                            // if showMotion is enabled make the pixel green
                            greenDest[j + offsX] = 1000.f + 25000.f * blend;
                            nonGreenDest0[j + offsX] = nonGreenDest1[j + offsX] = 0.f;
                        } else {
                            greenDest[j + offsX] = nonGreenDest0[j + offsX] = nonGreenDest1[j + offsX] = 1000.f + 25000.f * blend;
                        }
                    }

                    if(skipNext) {
                        // treat the horizontally next pixel also as motion
                        j++;
                        offset ^= 1;
                    }
                    // do not set the motion pixel values. They have already been set by demosaicer or showMotion
                    continue;
                }
//                } else if(showMotion && showOnlyMask) {
//                    greenDest[j + offsX] = nonGreenDest0[j + offsX] = nonGreenDest1[j + offsX] = 0.f;
//                    continue;
//                }

                if(adaptive && checkNonGreenHorizontal) {
                    float ng1 = riFrames[(offset << 1) + offset]->data[i][j + offset];
                    float ng0 = riFrames[((offset^1) << 1) + (offset^1)]->data[i][j + (offset^1)+1];
                    float ng2 = riFrames[((offset^1) << 1) + (offset^1)]->data[i][j + (offset^1)-1];
                    float diff0 = ng0 - ng1;
                    float diff2 = ng2 - ng1;
                    if(diff0 * diff2 >= 0.f) {
//                        float val = std::abs(diff0) > std::abs(diff2) ? ng0 : ng2;
                        float val = (ng0 + ng2) / 2.f;
                        float gridMax = colourDiff(ng1, val, true, stddevFactor, eperIsoNonGreen0, nreadIso, prnu, showMotion);
                        if(gridMax > 0.f) {
                            if(showMotion) {
                                float blend = gridMax * blendFactor;
                                if(!showOnlyMask) {
                                    // if showMotion is enabled make the pixel green
                                    nonGreenDest0[j + offsX] = 1000.f + 25000.f * blend;
                                    nonGreenDest1[j + offsX] = greenDest[j + offsX] = 0.f;
                                } else {
                                    greenDest[j + offsX] = nonGreenDest0[j + offsX] = nonGreenDest1[j + offsX] = 1000.f + 25000.f * blend;
                                }
                            }
                            continue;
                        }
                    }
                    ng1 = riFrames[2 - offset]->data[i + 1][j - offset + 1];
                    ng0 = riFrames[2 - (offset^1)]->data[i + 1][j - (offset^1) + 2];
                    ng2 = riFrames[2 - (offset^1)]->data[i + 1][j - (offset^1)];
                    diff0 = ng0 - ng1;
                    diff2 = ng2 - ng1;
                    if(diff0 * diff2 >= 0.f) {
//                        float val = std::abs(diff0) > std::abs(diff2) ? ng0 : ng2;
                        float val = (ng0 + ng2) / 2.f;
                        float gridMax = colourDiff(ng1, val, true, stddevFactor, eperIsoNonGreen2, nreadIso, prnu, showMotion);
                        if(gridMax > 0.f) {
                            if(showMotion) {
                                float blend = gridMax * blendFactor;
                                if(!showOnlyMask) {
                                    nonGreenDest1[j + offsX] = 1000.f + 25000.f * blend;
                                    nonGreenDest0[j + offsX] = greenDest[j + offsX] = 0.f;
                                } else {
                                    greenDest[j + offsX] = nonGreenDest0[j + offsX] = nonGreenDest1[j + offsX] = 1000.f + 25000.f * blend;
                                }
                            }
                            continue;
                        }
                    }

//                    float ngDiff = nonGreenDiff(ng0,ng1,ng2);
//                    float korr = log2Lut[((int)(riFrames[(offset << 1) + offset]->data[i][j + offset] * scaleNonGreen0))>>1];
//                    if(ngDiff > 0.5f - korr) {
//                        if(showMotion) {
//                            if(!showOnlyMask) {
//                                // if showMotion is enabled make the pixel green
//                                nonGreenDest0[j + offsX] = 1000.f + 25000.f;
//                                nonGreenDest1[j + offsX] = greenDest[j + offsX] = 0.f;
//                            } else {
//                                greenDest[j + offsX] = nonGreenDest0[j + offsX] = nonGreenDest1[j + offsX] = 1000.f + 25000.f;
//                            }
//                        }
//                        continue;
//                    }
//
//                    ng1 = riFrames[2 - offset]->data[i + 1][j - offset + 1];
//                    ng0 = riFrames[2 - (offset^1)]->data[i + 1][j - (offset^1) + 2];
//                    ng2 = riFrames[2 - (offset^1)]->data[i + 1][j - (offset^1)];
//                    ngDiff = nonGreenDiff(ng0,ng1,ng2);
//                    korr = log2Lut[((int)(riFrames[2 - offset]->data[i + 1][j - offset + 1] * scaleNonGreen2))>>1];
//                    if(ngDiff > 0.5f - korr) {
//                        if(showMotion) {
//                            if(!showOnlyMask) {
//                                // if showMotion is enabled make the pixel green
//                                nonGreenDest0[j + offsX] = 1000.f + 25000.f;
//                                nonGreenDest1[j + offsX] = greenDest[j + offsX] = 0.f;
//                            } else {
//                                greenDest[j + offsX] = nonGreenDest0[j + offsX] = nonGreenDest1[j + offsX] = 1000.f + 25000.f;
//                            }
//                        }
//                        continue;
//                    }
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
////                    if(colourDiff(ng1, fabsf(diff0) < fabsf(diff2) ? ng0 : ng2) > motionThreshold ) {
//                    gridMax = colourDiff(ng1, std::max(ng0, ng2));
////                    gridMax = colourDiff(ng1, ((ng0 + ng2) / 2.f));
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
////                    if(colourDiff(ng1, fabsf(diff0) < fabsf(diff2) ? ng0 : ng2) > motionThreshold ) {
//                    gridMax = colourDiff(ng1, std::max(ng0, ng2));
////                    gridMax = colourDiff(ng1, ((ng0 + ng2) / 2.f));
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

            if(showMotion && showOnlyMask) {
                greenDest[j + offsX] = nonGreenDest0[j + offsX] = nonGreenDest1[j + offsX] = 0.f;
                continue;
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
