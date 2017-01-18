////////////////////////////////////////////////////////////////
//
//  pentax pixelshift algorithm with motion detection
//
//  non adaptive mode is derived from dcrawps (https://github.com/tomtor/dcrawps), but with additional motion correction methods and adapted for RawTherapee data structures
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
#include "gauss.h"
#include "median.h"
#define BENCHMARK
#include "StopWatch.h"

namespace
{

float greenDiff(float a, float b, bool adaptive, float stddevFactor, float eperIso, float nreadIso, float prnu, bool showMotion, int x, int y)
{
    // calculate the difference between two green samples
    if(adaptive) {
        float gDiff = a - b;
        gDiff *= eperIso;
        gDiff *= gDiff;
        float avg = (a + b) / 2.f;
        avg *= eperIso;
        prnu *= avg;
        float stddev = stddevFactor * (avg + nreadIso + prnu * prnu);

//        if(x >= 4294 && x <= 4303 && y >= 3056 && y <= 3058) {
//            #pragma omp critical
//            std::cout << "x : " << x << " y : " << y << " stddev : " << stddev << " avg : " << avg << " gDiff  : " << gDiff << std::endl;
//        }

        float result = gDiff - stddev;

        if(!showMotion) {
            return result;
        } else if(result > 0.f) { // for the motion mask
            return std::fabs(a - b) / (std::max(a, b) + 0.01f);
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

float nonGreenDiff(float a, float b, float stddevFactor, float eperIso, float nreadIso, float prnu, bool showMotion)
{
    // calculate the difference between two nongreen samples
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
        return std::fabs(a - b) / (std::max(a, b) + 0.01f);
    } else {
        return 0.f;
    }
}

float nonGreenDiffCross(float right, float left, float top, float bottom, float centre, float stddevFactor, float eperIso, float nreadIso, float prnu, bool showMotion)
{
    // check non green cross
    float hDiff = (right + left) / 2.f - centre;
    hDiff *= eperIso;
    hDiff *= hDiff;
    float vDiff = (top + bottom) / 2.f - centre;
    vDiff *= eperIso;
    vDiff *= vDiff;
    float avg = (right + left + top + bottom) / 4.f;
    avg *= eperIso;
    prnu *= avg;
    float stddev = stddevFactor * (avg + nreadIso + prnu * prnu);
    float result = std::min(hDiff - stddev, vDiff - stddev);

    if(!showMotion) {
        return result;
    } else if(result > 0.f) { // for the motion mask
        return std::sqrt((result / (stddev + result + 0.01f)));
    } else {
        return 0.f;
    }
}

void paintMotionMask(int index, bool showMotion, float gridMax, bool showOnlyMask, float *maskDest, float *nonMaskDest0, float *nonMaskDest1)
{
    if(showMotion) {
        if(!showOnlyMask) {
            // if showMotion is enabled colourize the pixel
            maskDest[index] = 1000.f + 25000.f * gridMax;
            nonMaskDest1[index] = nonMaskDest0[index] = 0.f;
        } else {
            maskDest[index] = nonMaskDest0[index] = nonMaskDest1[index] = 1000.f + 25000.f * gridMax;
        }
    }
}

void invertMask(int xStart, int xEnd, int yStart, int yEnd, array2D<uint8_t> &maskIn, array2D<uint8_t> &maskOut)
{
    #pragma omp parallel for schedule(dynamic,16)

    for(int i = yStart; i < yEnd; ++i) {
        #pragma omp simd

        for(int j = xStart; j < xEnd; ++j) {
            maskOut[i][j] = ~maskIn[i][j];
        }
    }
}

void xorMasks(int xStart, int xEnd, int yStart, int yEnd, array2D<uint8_t> &maskIn, array2D<uint8_t> &maskOut)
{
    #pragma omp parallel for schedule(dynamic,16)

    for(int i = yStart; i < yEnd; ++i) {
        #pragma omp simd

        for(int j = xStart; j < xEnd; ++j) {
            maskOut[i][j] ^= maskIn[i][j];
        }
    }
}

void floodFill4Impl(int y, int x, int xStart, int xEnd, int yStart, int yEnd, array2D<uint8_t> &mask, std::stack<std::pair<uint16_t, uint16_t>, std::vector<std::pair<uint16_t, uint16_t>>> &coordStack)
{
    coordStack.emplace(x, y);

    while(!coordStack.empty()) {
        auto coord = coordStack.top();
        coordStack.pop();
        auto x = coord.first, y = coord.second;

        if (mask[y][x] == 255) {
            auto yUp = y - 1, yDown = y + 1;
            bool lastXUp = false, lastXDown = false, firstXUp = false, firstXDown = false;
            mask[y][x] = 0;
            if(yUp >= yStart && mask[yUp][x] == 255) {
                coordStack.emplace(x, yUp);
                firstXUp = lastXUp = true;
            }
            if(yDown < yEnd && mask[yDown][x] == 255) {
                coordStack.emplace(x, yDown);
                firstXDown = lastXDown = true;
            }
            auto xr = x + 1;
            
            while(xr < xEnd && mask[y][xr] == 255) {
                mask[y][xr] = 0;
                if(yUp >= yStart && mask[yUp][xr] == 255) {
                    if(!lastXUp) {
                        coordStack.emplace(xr, yUp);
                        lastXUp = true;
                    }
                } else {
                    lastXUp = false;
                }
                if(yDown < yEnd && mask[yDown][xr] == 255) {
                    if(!lastXDown) {
                        coordStack.emplace(xr, yDown);
                        lastXDown = true;
                    }
                } else {
                    lastXDown = false;
                }
                xr++;
            }

            auto xl = x - 1;
            lastXUp = firstXUp;
            lastXDown = firstXDown;
            while(xl >= xStart && mask[y][xl] == 255) {
                mask[y][xl] = 0;
                if(yUp >= yStart && mask[yUp][xl] == 255) {
                    if(!lastXUp) {
                        coordStack.emplace(xl, yUp);
                        lastXUp = true;
                    }
                } else {
                    lastXUp = false;
                }
                if(yDown < yEnd && mask[yDown][xl] == 255) {
                    if(!lastXDown) {
                        coordStack.emplace(xl, yDown);
                        lastXDown = true;
                    }
                } else {
                    lastXDown = false;
                }
                xl--;
            }
            mask[y][x] = 0;
        }
    }
}

void floodFill4(int xStart, int xEnd, int yStart, int yEnd, array2D<uint8_t> &mask)
{
    #pragma omp parallel
    {
    std::stack<std::pair<uint16_t, uint16_t>, std::vector<std::pair<uint16_t, uint16_t>>> coordStack;

    #pragma omp for schedule(dynamic,128) nowait
    for(uint16_t i = yStart;i<yEnd;i++)
        floodFill4Impl(i, xStart, xStart, xEnd, yStart, yEnd, mask, coordStack);

    #pragma omp for schedule(dynamic,128) nowait
    for(int16_t i = yEnd-1; i >= 0 ;i--)
        floodFill4Impl(i, xEnd - 1, xStart, xEnd, yStart, yEnd, mask, coordStack);

    #pragma omp sections nowait
    {
        #pragma omp section
        {
            uint16_t i = yStart;

            for(uint16_t j = xStart; j < xEnd; j++)
            {
                floodFill4Impl(i, j, xStart, xEnd, yStart, yEnd, mask, coordStack);
            }
        }
        #pragma omp section
        {
            uint16_t i = yStart;

            for(uint16_t j = xEnd - 1; j >= xStart; j--)
            {
                floodFill4Impl(i, j, xStart, xEnd, yStart, yEnd, mask, coordStack);
            }
        }
        #pragma omp section
        {
            uint16_t i = yEnd;

            for(uint16_t j = xStart; j < xEnd; j++)
            {
                floodFill4Impl(i, j, xStart, xEnd, yStart, yEnd, mask, coordStack);
            }
        }
        #pragma omp section
        {
            uint16_t i = yEnd;

            for(uint16_t j = xEnd - 1; j >= xStart; j--)
            {
                floodFill4Impl(i, j, xStart, xEnd, yStart, yEnd, mask, coordStack);
            }
        }
    }
    }
}


}

using namespace std;
using namespace rtengine;
#ifdef __OLDPS__
void RawImageSource::pixelshift(int winx, int winy, int winw, int winh, bool detectMotion, int motion, bool showMotion, bool showOnlyMask, unsigned int frame, RAWParams::BayerSensor::ePSMotionCorrection gridSize_, bool adaptive, float stddevFactorGreen, float stddevFactorRed, float stddevFactorBlue, float eperIso, float nreadIso, float prnu, const std::string &model, float rawWpCorrection, bool checkNonGreenHorizontal, bool checkNonGreenVertical, bool checkNonGreenCross)
{

    BENCHFUN

    static const float nReadK3II[] = { 3.4f,  // ISO 100
                                       3.1f,  // ISO 125
                                       2.5f,  // ISO 160
                                       2.5f,  // ISO 200
                                       2.5f,  // ISO 250
                                       2.5f,  // ISO 320
                                       2.3f,  // ISO 400
                                       2.5f,  // ISO 500
                                       2.3f,  // ISO 640
                                       2.3f,  // ISO 800
                                       2.4f,  // ISO 1000
                                       2.3f,  // ISO 1250
                                       1.75f, // ISO 1600
                                       1.75f, // ISO 2000
                                       1.75f, // ISO 2500
                                       1.75f, // ISO 3200
                                       1.75f, // ISO 4000
                                       1.75f, // ISO 5000
                                       1.75f, // ISO 6400
                                       1.75f, // ISO 8000
                                       1.75f, // ISO 10000
                                       1.5f,  // ISO 12800
                                       1.5f,  // ISO 16000
                                       1.5f,  // ISO 20000
                                       1.5f,  // ISO 25600
                                       1.5f,  // ISO 32000
                                       1.5f,  // ISO 40000
                                       1.5f,  // ISO 51200
                                       1.5f   // ISO > 51200 (we get a max ISO value of 65535 from dcraw)
                                     };

    static const float ePerIsoK3II = 0.35f;

    static const float nReadK1[] = {   3.45f, // ISO 100
                                       3.15f, // ISO 125
                                       3.45f, // ISO 160
                                       3.0f,  // ISO 200
                                       3.0f,  // ISO 250
                                       3.0f,  // ISO 320
                                       2.7f,  // ISO 400
                                       2.7f,  // ISO 500
                                       2.7f,  // ISO 640
                                       2.5f,  // ISO 800
                                       2.5f,  // ISO 1000
                                       2.5f,  // ISO 1250
                                       2.4f,  // ISO 1600
                                       2.4f,  // ISO 2000
                                       2.4f,  // ISO 2500
                                       2.4f,  // ISO 3200
                                       2.4f,  // ISO 4000
                                       2.4f,  // ISO 5000
                                       2.4f,  // ISO 6400
                                       2.4f,  // ISO 8000
                                       2.4f,  // ISO 10000
                                       2.4f,  // ISO 12800
                                       2.4f,  // ISO 16000
                                       2.4f,  // ISO 20000
                                       2.4f,  // ISO 25600
                                       2.4f,  // ISO 32000
                                       2.4f,  // ISO 40000
                                       2.4f,  // ISO 51200
                                       2.4f   // ISO > 51200 (we get a max ISO value of 65535 from dcraw)
                                   };

    static const float ePerIsoK1 = 0.75f;

    static const float nReadK70[] = {  3.0f, // ISO 100
                                       3.0f, // ISO 125
                                       3.0f, // ISO 160
                                       3.0f,  // ISO 200
                                       3.0f,  // ISO 250
                                       3.0f,  // ISO 320
                                       3.0f,  // ISO 400
                                       3.0f,  // ISO 500
                                       3.0f,  // ISO 640
                                       3.0f,  // ISO 800
                                       3.0f,  // ISO 1000
                                       3.0f,  // ISO 1250
                                       3.0f,  // ISO 1600
                                       3.0f,  // ISO 2000
                                       3.0f,  // ISO 2500
                                       3.0f,  // ISO 3200
                                       3.0f,  // ISO 4000
                                       3.0f,  // ISO 5000
                                       3.0f,  // ISO 6400
                                       3.0f,  // ISO 8000
                                       3.0f,  // ISO 10000
                                       3.0f,  // ISO 12800
                                       3.0f,  // ISO 16000
                                       3.0f,  // ISO 20000
                                       3.0f,  // ISO 25600
                                       3.0f,  // ISO 32000
                                       3.0f,  // ISO 40000
                                       3.0f,  // ISO 51200
                                       3.0f   // ISO > 51200 (we get a max ISO value of 65535 from dcraw)
                                    };

    static const float ePerIsoK70 = 0.5f;

    if (plistener) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::pixelshift]));
        plistener->setProgress(0.0);
    }


    const bool skip = (gridSize_ == RAWParams::BayerSensor::ePSMotionCorrection::Grid1x2);
    int gridSize = 1;

    switch (gridSize_) {
        case RAWParams::BayerSensor::ePSMotionCorrection::Grid1x1:
        case RAWParams::BayerSensor::ePSMotionCorrection::Grid1x2:
            gridSize = 1;
            break;

        case RAWParams::BayerSensor::ePSMotionCorrection::Grid3x3:
            gridSize = 3;
            break;

        case RAWParams::BayerSensor::ePSMotionCorrection::Grid5x5:
            gridSize = 5;
            break;

        case RAWParams::BayerSensor::ePSMotionCorrection::Grid7x7:
            gridSize = 7;
    }

    // Lookup table for non adaptive (slider) mode
    LUTf log2Lut(32768, LUT_CLIP_BELOW | LUT_CLIP_ABOVE);

    if(detectMotion && !adaptive) {
        const float lutStrength = 2.f;
        log2Lut[0] = 0;

        for(int i = 2; i < 65536; i += 2) {
            log2Lut[i >> 1] = lutStrength * log2(i) / 100.f;
        }
    }

    const float scaleGreen = 1.f / scale_mul[1];

    float nRead;
    float eperIsoModel;

    int nReadIndex = static_cast<int>(round(log2(idata->getISOSpeed() /  100.f) * 3.f));

    if(model.find("K-3") != string::npos) {
        nRead = nReadK3II[nReadIndex];
        eperIsoModel = ePerIsoK3II;
    } else if(model.find("K-1") != string::npos) {
        nRead = nReadK1[nReadIndex];
        eperIsoModel = ePerIsoK1;
    } else {
        nRead = nReadK70[nReadIndex];
        eperIsoModel = ePerIsoK70;
    }

    nRead *= pow(2.f, nreadIso);
    eperIsoModel *= pow(2.f, eperIso);
    eperIso = eperIsoModel * (100.f / (rawWpCorrection * idata->getISOSpeed()));
    float eperIsoGreen = eperIso * scaleGreen;

//    printf("Pixelshift parameters : gridSize %d\tadaptive %d\tstdDevFactorGreen %f\telectrons %1.8f\tnread %f\tprnu %1.1f%%\n", gridSize, adaptive, stddevFactorGreen, eperIso, nRead, prnu);

    prnu /= 100.f;
    stddevFactorGreen *= stddevFactorGreen;
    stddevFactorRed *= stddevFactorRed;
    stddevFactorBlue *= stddevFactorBlue;


    nRead *= nRead;

    // If the values of two corresponding green pixels differ my more then motionThreshold %, the pixel will be treated as a badGreen pixel
    float motionThreshold = 1.f - (motion / 100.f);
    // For shades of green motion indicators
    const float blendFactor = ((adaptive || motion == 0.f) ? 1.f : 1.f / (1.f - motionThreshold));

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
        float stddevFactorNonGreen0 = stddevFactorRed;
        float stddevFactorNonGreen2 = stddevFactorBlue;
        bool blueRow = false;

        if (c == 2 || ((c & 1) && FC(i, j + 1) == 2)) {
            // row with blue pixels => swap destination pointers for non green pixels
            blueRow = true;
            std::swap(nonGreenDest0, nonGreenDest1);
            std::swap(scaleNonGreen0, scaleNonGreen2);
            std::swap(eperIsoNonGreen0, eperIsoNonGreen2);
            std::swap(stddevFactorNonGreen0, stddevFactorNonGreen2);
        }

        // offset to keep the code short. It changes its value between 0 and 1 for each iteration of the loop
        unsigned int offset = (c & 1);

        float greenDifMax[gridSize];

        // green channel motion detection checks the grid around the pixel for differences in green channels
        if(detectMotion || adaptive) {
            if(gridSize == 3) {
                // compute maximum of differences for first two columns of 3x3 grid
                greenDifMax[0] =  std::max({greenDiff((*rawDataFrames[1 - offset])[i - offset][j - 1], (*rawDataFrames[3 - offset])[i + offset - 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset][j - 1], (*rawDataFrames[2 + offset])[i - offset + 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 2][j - 1], (*rawDataFrames[3 - offset])[i + offset + 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[1] =  std::max({greenDiff((*rawDataFrames[0 + offset])[i + offset - 1][j], (*rawDataFrames[2 + offset])[i - offset][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 1][j], (*rawDataFrames[3 - offset])[i + offset][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset + 1][j], (*rawDataFrames[2 + offset])[i - offset + 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
            } else if(gridSize == 5) {
                // compute maximum of differences for first four columns of 5x5 grid
                greenDifMax[0] =  std::max({greenDiff((*rawDataFrames[1 - offset])[i - offset - 1][j - 2], (*rawDataFrames[3 - offset])[i + offset - 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset - 1][j - 2], (*rawDataFrames[2 + offset])[i - offset][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 1][j - 2], (*rawDataFrames[3 - offset])[i + offset][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset + 1][j - 2], (*rawDataFrames[2 + offset])[i - offset + 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 3][j - 2], (*rawDataFrames[3 - offset])[i + offset + 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[1] =  std::max({greenDiff((*rawDataFrames[0 + offset])[i + offset - 2][j - 1], (*rawDataFrames[2 + offset])[i - offset - 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset][j - 1], (*rawDataFrames[3 - offset])[i + offset - 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset][j - 1], (*rawDataFrames[2 + offset])[i - offset + 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 2][j - 1], (*rawDataFrames[3 - offset])[i + offset + 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset + 2][j - 1], (*rawDataFrames[2 + offset])[i - offset + 3][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[2] =  std::max({greenDiff((*rawDataFrames[1 - offset])[i - offset - 1][j], (*rawDataFrames[3 - offset])[i + offset - 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset - 1][j], (*rawDataFrames[2 + offset])[i - offset][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 1][j], (*rawDataFrames[3 - offset])[i + offset][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset + 1][j], (*rawDataFrames[2 + offset])[i - offset + 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 3][j], (*rawDataFrames[3 - offset])[i + offset + 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[3] =  std::max({greenDiff((*rawDataFrames[0 + offset])[i + offset - 2][j + 1], (*rawDataFrames[2 + offset])[i - offset - 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset][j + 1], (*rawDataFrames[3 - offset])[i + offset - 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset][j + 1], (*rawDataFrames[2 + offset])[i - offset + 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 2][j + 1], (*rawDataFrames[3 - offset])[i + offset + 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset + 2][j + 1], (*rawDataFrames[2 + offset])[i - offset + 3][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
            } else if(gridSize == 7) {
                // compute maximum of differences for first six columns of 7x7 grid
                greenDifMax[0] =  std::max({greenDiff((*rawDataFrames[1 - offset])[i - offset - 2][j - 3], (*rawDataFrames[3 - offset])[i + offset - 3][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset - 2][j - 3], (*rawDataFrames[2 + offset])[i - offset - 1][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset][j - 3], (*rawDataFrames[3 - offset])[i + offset - 1][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset][j - 3], (*rawDataFrames[2 + offset])[i - offset + 1][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 2][j - 3], (*rawDataFrames[3 - offset])[i + offset + 1][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset + 2][j - 3], (*rawDataFrames[2 + offset])[i - offset + 3][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 4][j - 3], (*rawDataFrames[3 - offset])[i + offset + 3][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[1] =  std::max({greenDiff((*rawDataFrames[0 + offset])[i + offset - 3][j - 2], (*rawDataFrames[2 + offset])[i - offset - 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset - 1][j - 2], (*rawDataFrames[3 - offset])[i + offset - 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset - 1][j - 2], (*rawDataFrames[2 + offset])[i - offset][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 1][j - 2], (*rawDataFrames[3 - offset])[i + offset][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset + 1][j - 2], (*rawDataFrames[2 + offset])[i - offset + 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 3][j - 2], (*rawDataFrames[3 - offset])[i + offset + 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset + 3][j - 2], (*rawDataFrames[2 + offset])[i - offset + 4][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[2] =  std::max({greenDiff((*rawDataFrames[1 - offset])[i - offset - 2][j - 1], (*rawDataFrames[3 - offset])[i + offset - 3][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset - 2][j - 1], (*rawDataFrames[2 + offset])[i - offset - 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset][j - 1], (*rawDataFrames[3 - offset])[i + offset - 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset][j - 1], (*rawDataFrames[2 + offset])[i - offset + 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 2][j - 1], (*rawDataFrames[3 - offset])[i + offset + 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset + 2][j - 1], (*rawDataFrames[2 + offset])[i - offset + 3][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 4][j - 1], (*rawDataFrames[3 - offset])[i + offset + 3][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[3] =  std::max({greenDiff((*rawDataFrames[0 + offset])[i + offset - 3][j], (*rawDataFrames[2 + offset])[i - offset - 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset - 1][j], (*rawDataFrames[3 - offset])[i + offset - 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset - 1][j], (*rawDataFrames[2 + offset])[i - offset][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 1][j], (*rawDataFrames[3 - offset])[i + offset][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset + 1][j], (*rawDataFrames[2 + offset])[i - offset + 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 3][j], (*rawDataFrames[3 - offset])[i + offset + 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset + 3][j], (*rawDataFrames[2 + offset])[i - offset + 4][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[4] =  std::max({greenDiff((*rawDataFrames[1 - offset])[i - offset - 2][j + 1], (*rawDataFrames[3 - offset])[i + offset - 3][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset - 2][j + 1], (*rawDataFrames[2 + offset])[i - offset - 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset][j + 1], (*rawDataFrames[3 - offset])[i + offset - 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset][j + 1], (*rawDataFrames[2 + offset])[i - offset + 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 2][j + 1], (*rawDataFrames[3 - offset])[i + offset + 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset + 2][j + 1], (*rawDataFrames[2 + offset])[i - offset + 3][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 4][j + 1], (*rawDataFrames[3 - offset])[i + offset + 3][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[5] =  std::max({greenDiff((*rawDataFrames[0 + offset])[i + offset - 3][j + 2], (*rawDataFrames[2 + offset])[i - offset - 2][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset - 1][j + 2], (*rawDataFrames[3 - offset])[i + offset - 2][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset - 1][j + 2], (*rawDataFrames[2 + offset])[i - offset][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 1][j + 2], (*rawDataFrames[3 - offset])[i + offset][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset + 1][j + 2], (*rawDataFrames[2 + offset])[i - offset + 2][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[1 - offset])[i - offset + 3][j + 2], (*rawDataFrames[3 - offset])[i + offset + 2][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff((*rawDataFrames[0 + offset])[i + offset + 3][j + 2], (*rawDataFrames[2 + offset])[i - offset + 4][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
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

                if(gridSize < 2) {
                    // compute difference for current pixel and skip next pixel, that's the method from dcrawps
                    gridMax = greenDiff((*rawDataFrames[1 - offset])[i - offset + 1][j], (*rawDataFrames[3 - offset])[i + offset][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i);
                    skipNext = skip;
                } else if(gridSize == 3) {
                    // compute maximum of differences for third column of 3x3 grid and save at position lastIndex
                    greenDifMax[lastIndex] = std::max({greenDiff((*rawDataFrames[1 - offset])[i - offset][j + 1], (*rawDataFrames[3 - offset])[i + offset - 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff((*rawDataFrames[0 + offset])[i + offset][j + 1], (*rawDataFrames[2 + offset])[i - offset + 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff((*rawDataFrames[1 - offset])[i - offset + 2][j + 1], (*rawDataFrames[3 - offset])[i + offset + 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                                      });
                    gridMax = std::max({greenDifMax[0], greenDifMax[1], greenDifMax[2]});
                } else if(gridSize == 5) {
                    // compute maximum of differences for fifth column of 5x5 grid and save at position lastIndex
                    greenDifMax[lastIndex] = std::max({greenDiff((*rawDataFrames[1 - offset])[i - offset - 1][j + 2], (*rawDataFrames[3 - offset])[i + offset - 2][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff((*rawDataFrames[0 + offset])[i + offset - 1][j + 2], (*rawDataFrames[2 + offset])[i - offset][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff((*rawDataFrames[1 - offset])[i - offset + 1][j + 2], (*rawDataFrames[3 - offset])[i + offset][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff((*rawDataFrames[0 + offset])[i + offset + 1][j + 2], (*rawDataFrames[2 + offset])[i - offset + 2][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff((*rawDataFrames[1 - offset])[i - offset + 3][j + 2], (*rawDataFrames[3 - offset])[i + offset + 2][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                                      });
                    gridMax = std::max({greenDifMax[0], greenDifMax[1], greenDifMax[2], greenDifMax[3], greenDifMax[4]});
                } else if(gridSize == 7) {
                    // compute maximum of differences for 7th column of 7x7 grid and save at position lastIndex
                    greenDifMax[lastIndex] = std::max({greenDiff((*rawDataFrames[1 - offset])[i - offset - 2][j + 3], (*rawDataFrames[3 - offset])[i + offset - 3][j + 4], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff((*rawDataFrames[0 + offset])[i + offset - 2][j + 3], (*rawDataFrames[2 + offset])[i - offset - 1][j + 4], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff((*rawDataFrames[1 - offset])[i - offset][j + 3], (*rawDataFrames[3 - offset])[i + offset - 1][j + 4], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff((*rawDataFrames[0 + offset])[i + offset][j + 3], (*rawDataFrames[2 + offset])[i - offset + 1][j + 4], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff((*rawDataFrames[1 - offset])[i - offset + 2][j + 3], (*rawDataFrames[3 - offset])[i + offset + 1][j + 4], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff((*rawDataFrames[0 + offset])[i + offset + 2][j + 3], (*rawDataFrames[2 + offset])[i - offset + 3][j + 4], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff((*rawDataFrames[1 - offset])[i - offset + 4][j + 3], (*rawDataFrames[3 - offset])[i + offset + 3][j + 4], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                                      });
                    gridMax = std::max({greenDifMax[0], greenDifMax[1], greenDifMax[2], greenDifMax[3], greenDifMax[4], greenDifMax[5], greenDifMax[6]});
                }


                // adjust index for next column
                lastIndex ++;
                lastIndex = lastIndex == gridSize ? 0 : lastIndex;

                // increase motion detection dependent on brightness
                if(!adaptive) {
                    korr = log2Lut[((int)((*rawDataFrames[1 - offset])[i - offset + 1][j] * scaleGreen)) >> 1];
                }

                if (gridMax > thresh - korr) {
                    // at least one of the tested pixels of the grid is detected as motion
                    paintMotionMask(j + offsX, showMotion, (gridMax - thresh + korr) * blendFactor, showOnlyMask, greenDest, nonGreenDest0, nonGreenDest1);

                    if(skipNext) {
                        // treat the horizontally next pixel also as motion
                        j++;
                        paintMotionMask(j + offsX, showMotion, (gridMax - thresh + korr) * blendFactor, showOnlyMask, greenDest, nonGreenDest0, nonGreenDest1);
                        offset ^= 1;
                    }

                    // do not set the motion pixel values. They have already been set by demosaicer or showMotion
                    continue;
                }

                if(adaptive && checkNonGreenCross) {
                    float ngCentre = (*rawDataFrames[(offset << 1) + offset])[i][j + offset];
                    float ngRight = (*rawDataFrames[((offset ^ 1) << 1) + (offset ^ 1)])[i][j + (offset ^ 1) + 1];
                    float ngLeft = (*rawDataFrames[((offset ^ 1) << 1) + (offset ^ 1)])[i][j + (offset ^ 1) - 1];
                    float diffRight = ngRight - ngCentre;
                    float diffLeft = ngLeft - ngCentre;
                    float diffHorNg0 = -1.f;

                    if(diffRight * diffLeft >= 0.f) {
                        float avg = (ngRight + ngLeft) / 2.f;
                        diffHorNg0 = nonGreenDiff(ngCentre, avg, stddevFactorNonGreen0, eperIsoNonGreen0, nRead, prnu, showMotion);

//                        if(diff > 0.f) {
//                            paintMotionMask(j + offsX, showMotion, diff, showOnlyMask, nonGreenDest0, nonGreenDest1, greenDest);
//                            continue;
//                        }
                    }

                    float diffHorNg1 = -1.f;
                    ngCentre = (*rawDataFrames[2 - offset])[i + 1][j - offset + 1];
                    ngRight = (*rawDataFrames[2 - (offset ^ 1)])[i + 1][j - (offset ^ 1) + 2];
                    ngLeft = (*rawDataFrames[2 - (offset ^ 1)])[i + 1][j - (offset ^ 1)];
                    diffRight = ngRight - ngCentre;
                    diffLeft = ngLeft - ngCentre;

                    if(diffRight * diffLeft >= 0.f) {
                        float avg = (ngRight + ngLeft) / 2.f;
                        float diffHorNg1 = nonGreenDiff(ngCentre, avg, stddevFactorNonGreen2, eperIsoNonGreen2, nRead, prnu, showMotion);

//                        if(diff > 0.f) {
//                            paintMotionMask(j + offsX, showMotion, diff, showOnlyMask, nonGreenDest1, nonGreenDest0, greenDest);
//                            continue;
//                        }
                    }

                    if( diffHorNg0 * diffHorNg1 < 0.f) {
                        paintMotionMask(j + offsX, showMotion, 1.f, showOnlyMask, nonGreenDest0, nonGreenDest1, greenDest);
                        continue;

                    }

//                    bool motion = false;
//                    float ngCentre = (*rawDataFrames[(offset << 1) + offset])[i][j + offset];
//                    float ngRight = (*rawDataFrames[((offset ^ 1) << 1) + (offset ^ 1)])[i][j + (offset ^ 1) + 1];
//                    float ngLeft = (*rawDataFrames[((offset ^ 1) << 1) + (offset ^ 1)])[i][j + (offset ^ 1) - 1];
//                    float ngTop = (*rawDataFrames[((offset << 1) + offset) ^ 1])[i][j + offset];
//                    float ngBottom = (*rawDataFrames[((offset << 1) + offset) ^ 1])[i + 2][j + offset];
//                    float diff = nonGreenDiffCross(ngRight, ngLeft, ngTop, ngBottom, ngCentre, stddevFactorNonGreen0, eperIsoNonGreen0, nRead, prnu, showMotion);
//
//                    if(diff > 0.f) {
//                        motion = true;
//                        paintMotionMask(j + offsX, showMotion, diff, showOnlyMask, nonGreenDest0, nonGreenDest1, greenDest);
//                        continue;
//                    }
//
//                    ngCentre = (*rawDataFrames[2 - offset])[i + 1][j - offset + 1];
//                    ngRight = (*rawDataFrames[2 - (offset ^ 1)])[i + 1][j - (offset ^ 1) + 2];
//                    ngLeft = (*rawDataFrames[2 - (offset ^ 1)])[i + 1][j - (offset ^ 1)];
//                    ngTop = (*rawDataFrames[3 - ((offset << 1) + offset)])[i - 1][j - offset + 1];
//                    ngBottom = (*rawDataFrames[3 - ((offset << 1) + offset)])[i + 1][j - offset + 1];
//                    diff = nonGreenDiffCross(ngRight, ngLeft, ngTop, ngBottom, ngCentre, stddevFactorNonGreen2, eperIsoNonGreen2, nRead, prnu, showMotion);
//
////                    if(diff > 0.f) {
//                    if((diff > 0.f && !motion) || (diff <= 0.f && motion) ) {
//                        paintMotionMask(j + offsX, showMotion, diff, showOnlyMask, nonGreenDest1, nonGreenDest0, greenDest);
//                        continue;
//                    }
                }

                if(adaptive && checkNonGreenHorizontal) {
//                    float lg = ((*rawDataFrames[1 - (offset^1)])[i - (offset^1) + 1][j - 1] + (*rawDataFrames[3 - (offset^1)])[i + (offset^1)][j]) / 2.f;
//                    float cg = ((*rawDataFrames[1 - offset])[i - offset + 1][j] + (*rawDataFrames[3 - offset])[i + offset][j + 1]) / 2.f;
//                    float rg = ((*rawDataFrames[1 - (offset^1)])[i - (offset^1) + 1][j + 1] + (*rawDataFrames[3 - (offset^1)])[i + (offset^1)][j + 2]) / 2.f;
//
//                    float lr = (*rawDataFrames[((offset^1) << 1) + (offset^1)])[i][j + (offset^1) - 1];
//                    float cr = (*rawDataFrames[(offset << 1) + offset])[i][j + offset];
//                    float rr = (*rawDataFrames[((offset^1) << 1) + (offset^1)])[i][j + (offset^1) + 1];
//
//                    float lb = (*rawDataFrames[2 - (offset^1)])[i + 1][j - (offset^1)];
//                    float cb = (*rawDataFrames[2 - offset])[i + 1][j - offset + 1];
//                    float rb = (*rawDataFrames[2 - (offset^1)])[i + 1][j - (offset^1) + 2];
//
//                    if(blueRow) {
//                        std::swap(lr, lb);
//                        std::swap(cr, cb);
//                        std::swap(rr, rb);
//                    }
//
//                    float lh = Color::rgb2h(lr, lg, lb);
//                    float ch = Color::rgb2h(cr, cg, cb);
//                    float rh = Color::rgb2h(rr, rg, rb);
//
//                    float lHueDiff = lh - ch;
//                    float rHueDiff = rh - ch;
//                    if(lHueDiff * rHueDiff > 0.f) {
//                        if(std::fabs(lHueDiff) > 0.5f && std::fabs(rHueDiff) > 0.5f/* && std::fabs(lHueDiff) < 3.f && std::fabs(rHueDiff) < 3.f*/) {
//                            paintMotionMask(j + offsX, showMotion, 1.f, showOnlyMask, nonGreenDest0, nonGreenDest1, greenDest);
//                            continue;
//                        }
//                    }
                    float ngCentre = (*rawDataFrames[(offset << 1) + offset])[i][j + offset];
                    float ngRight = (*rawDataFrames[((offset ^ 1) << 1) + (offset ^ 1)])[i][j + (offset ^ 1) + 1];
                    float ngLeft = (*rawDataFrames[((offset ^ 1) << 1) + (offset ^ 1)])[i][j + (offset ^ 1) - 1];
                    float diffRight = ngRight - ngCentre;
                    float diffLeft = ngLeft - ngCentre;

                    if(diffRight * diffLeft >= 0.f) {
                        float avg = (ngRight + ngLeft) / 2.f;
                        float diff = nonGreenDiff(ngCentre, avg, stddevFactorNonGreen0, eperIsoNonGreen0, nRead, prnu, showMotion);

                        if(diff > 0.f) {
                            paintMotionMask(j + offsX, showMotion, diff, showOnlyMask, nonGreenDest0, nonGreenDest1, greenDest);
                            continue;
                        }
                    }

                    ngCentre = (*rawDataFrames[2 - offset])[i + 1][j - offset + 1];
                    ngRight = (*rawDataFrames[2 - (offset ^ 1)])[i + 1][j - (offset ^ 1) + 2];
                    ngLeft = (*rawDataFrames[2 - (offset ^ 1)])[i + 1][j - (offset ^ 1)];
                    diffRight = ngRight - ngCentre;
                    diffLeft = ngLeft - ngCentre;

                    if(diffRight * diffLeft >= 0.f) {
                        float avg = (ngRight + ngLeft) / 2.f;
                        float diff = nonGreenDiff(ngCentre, avg, stddevFactorNonGreen2, eperIsoNonGreen2, nRead, prnu, showMotion);

                        if(diff > 0.f) {
                            paintMotionMask(j + offsX, showMotion, diff, showOnlyMask, nonGreenDest1, nonGreenDest0, greenDest);
                            continue;
                        }
                    }
                }

                if(adaptive && checkNonGreenVertical) {
                    float ngCentre = (*rawDataFrames[(offset << 1) + offset])[i][j + offset];
                    float ngTop = (*rawDataFrames[((offset << 1) + offset) ^ 1])[i][j + offset];
                    float ngBottom = (*rawDataFrames[((offset << 1) + offset) ^ 1])[i + 2][j + offset];

                    float diffTop = ngTop - ngCentre;
                    float diffBottom = ngBottom - ngCentre;

                    if(diffTop * diffBottom >= 0.f) {
                        float avg = (ngTop + ngBottom) / 2.f;
                        float diff = nonGreenDiff(ngCentre, avg, stddevFactorNonGreen0, eperIsoNonGreen0, nRead, prnu, showMotion);

                        if(diff > 0.f) {
                            paintMotionMask(j + offsX, showMotion, diff, showOnlyMask, nonGreenDest0, nonGreenDest1, greenDest);
                            continue;
                        }
                    }

                    ngCentre = (*rawDataFrames[2 - offset])[i + 1][j - offset + 1];
                    ngTop = (*rawDataFrames[3 - ((offset << 1) + offset)])[i - 1][j - offset + 1];
                    ngBottom = (*rawDataFrames[3 - ((offset << 1) + offset)])[i + 1][j - offset + 1];

                    diffTop = ngTop - ngCentre;
                    diffBottom = ngBottom - ngCentre;

                    if(diffTop * diffBottom >= 0.f) {
                        float avg = (ngTop + ngBottom) / 2.f;
                        float diff = nonGreenDiff(ngCentre, avg, stddevFactorNonGreen2, eperIsoNonGreen2, nRead, prnu, showMotion);

                        if(diff > 0.f) {
                            paintMotionMask(j + offsX, showMotion, diff, showOnlyMask, nonGreenDest1, nonGreenDest0, greenDest);
                            continue;
                        }
                    }
                }
            }

            if(showMotion && showOnlyMask) {
                greenDest[j + offsX] = nonGreenDest0[j + offsX] = nonGreenDest1[j + offsX] = 0.f;
                continue;
            }

            // motion correction disabled or no motion detected => combine the values from the four pixelshift frames
            greenDest[j + offsX] = ((*rawDataFrames[1 - offset])[i - offset + 1][j] + (*rawDataFrames[3 - offset])[i + offset][j + 1]) / 2.f;
            nonGreenDest0[j + offsX] = (*rawDataFrames[(offset << 1) + offset])[i][j + offset];
            nonGreenDest1[j + offsX] = (*rawDataFrames[2 - offset])[i + 1][j - offset + 1];
        }
    }

    if(plistener) {
        plistener->setProgress(1.0);
    }
}
#else
void RawImageSource::pixelshift(int winx, int winy, int winw, int winh, const RAWParams::BayerSensor &bayerParams, unsigned int frame, const std::string &model, float rawWpCorrection)
{

    BENCHFUN

    const bool detectMotion = bayerParams.pixelShiftMotion > 0;
    const int motion = bayerParams.pixelShiftMotion;
    const bool showMotion = bayerParams.pixelshiftShowMotion;
    const bool showOnlyMask = bayerParams.pixelshiftShowMotionMaskOnly;
    const RAWParams::BayerSensor::ePSMotionCorrection gridSize_ = bayerParams.pixelShiftMotionCorrection;
    const bool adaptive = bayerParams.pixelShiftAutomatic;
    float stddevFactorGreen = bayerParams.pixelShiftStddevFactorGreen;
    float stddevFactorRed = bayerParams.pixelShiftStddevFactorRed;
    float stddevFactorBlue = bayerParams.pixelShiftStddevFactorBlue;
    float eperIso = bayerParams.pixelShiftEperIso;
    float nreadIso =  bayerParams.pixelShiftNreadIso;
    float prnu = bayerParams.pixelShiftPrnu;
    const bool checkNonGreenHorizontal = bayerParams.pixelShiftNonGreenHorizontal;
    const bool checkNonGreenVertical = bayerParams.pixelShiftNonGreenVertical;
    const bool checkNonGreenCross = bayerParams.pixelShiftNonGreenCross;
    const bool checkNonGreenAmaze = bayerParams.pixelShiftNonGreenAmaze;
    const bool checkNonGreenCross2 = bayerParams.pixelShiftNonGreenCross2;
    const bool checkGreen = bayerParams.pixelShiftGreen;
    const float redBlueWeight = bayerParams.pixelShiftRedBlueWeight;
    const bool blurMap = bayerParams.pixelShiftBlur;
    const float sigma = bayerParams.pixelShiftSigma;
    const float threshold = bayerParams.pixelShiftSum;
    const bool experimental0 = bayerParams.pixelShiftExp0;
    const bool holeFill = bayerParams.pixelShiftHoleFill;

    static const float nReadK3II[] = { 3.4f,  // ISO 100
                                       3.1f,  // ISO 125
                                       2.5f,  // ISO 160
                                       2.5f,  // ISO 200
                                       2.5f,  // ISO 250
                                       2.5f,  // ISO 320
                                       2.3f,  // ISO 400
                                       2.5f,  // ISO 500
                                       2.3f,  // ISO 640
                                       2.3f,  // ISO 800
                                       2.4f,  // ISO 1000
                                       2.3f,  // ISO 1250
                                       1.75f, // ISO 1600
                                       1.75f, // ISO 2000
                                       1.75f, // ISO 2500
                                       1.75f, // ISO 3200
                                       1.75f, // ISO 4000
                                       1.75f, // ISO 5000
                                       1.75f, // ISO 6400
                                       1.75f, // ISO 8000
                                       1.75f, // ISO 10000
                                       1.5f,  // ISO 12800
                                       1.5f,  // ISO 16000
                                       1.5f,  // ISO 20000
                                       1.5f,  // ISO 25600
                                       1.5f,  // ISO 32000
                                       1.5f,  // ISO 40000
                                       1.5f,  // ISO 51200
                                       1.5f   // ISO > 51200 (we get a max ISO value of 65535 from dcraw)
                                     };

    static const float ePerIsoK3II = 4 * 0.35f;

    static const float nReadK1[] = {   3.45f, // ISO 100
                                       3.15f, // ISO 125
                                       3.45f, // ISO 160
                                       3.0f,  // ISO 200
                                       3.0f,  // ISO 250
                                       3.0f,  // ISO 320
                                       2.7f,  // ISO 400
                                       2.7f,  // ISO 500
                                       2.7f,  // ISO 640
                                       2.5f,  // ISO 800
                                       2.5f,  // ISO 1000
                                       2.5f,  // ISO 1250
                                       2.4f,  // ISO 1600
                                       2.4f,  // ISO 2000
                                       2.4f,  // ISO 2500
                                       2.4f,  // ISO 3200
                                       2.4f,  // ISO 4000
                                       2.4f,  // ISO 5000
                                       2.4f,  // ISO 6400
                                       2.4f,  // ISO 8000
                                       2.4f,  // ISO 10000
                                       2.4f,  // ISO 12800
                                       2.4f,  // ISO 16000
                                       2.4f,  // ISO 20000
                                       2.4f,  // ISO 25600
                                       2.4f,  // ISO 32000
                                       2.4f,  // ISO 40000
                                       2.4f,  // ISO 51200
                                       2.4f,  // ISO 64000
                                       2.4f,  // ISO 80000
                                       2.4f,  // ISO 102400
                                       2.4f,  // ISO 128000
                                       2.4f,  // ISO 160000
                                       2.4f   // ISO 204800
                                   };

    static const float ePerIsoK1 = 4 * 0.75f;

    static const float nReadK70[] = {  4.0f, // ISO 100
                                       4.0f, // ISO 125
                                       4.0f, // ISO 160
                                       4.0f,  // ISO 200
                                       4.0f,  // ISO 250
                                       4.0f,  // ISO 320
                                       4.0f,  // ISO 400
                                       4.0f,  // ISO 500
                                       4.0f,  // ISO 640
                                       3.0f,  // ISO 800
                                       3.0f,  // ISO 1000
                                       3.0f,  // ISO 1250
                                       3.0f,  // ISO 1600
                                       3.0f,  // ISO 2000
                                       3.0f,  // ISO 2500
                                       3.0f,  // ISO 3200
                                       3.0f,  // ISO 4000
                                       3.0f,  // ISO 5000
                                       3.0f,  // ISO 6400
                                       3.0f,  // ISO 8000
                                       3.0f,  // ISO 10000
                                       3.0f,  // ISO 12800
                                       3.0f,  // ISO 16000
                                       3.0f,  // ISO 20000
                                       3.0f,  // ISO 25600
                                       3.0f,  // ISO 32000
                                       3.0f,  // ISO 40000
                                       3.0f,  // ISO 51200
                                       3.0f   // ISO > 51200 (we get a max ISO value of 65535 from dcraw)
                                    };

    static const float ePerIsoK70 = 4 * 0.5f;

    if (plistener) {
        plistener->setProgressStr (Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), RAWParams::BayerSensor::methodstring[RAWParams::BayerSensor::pixelshift]));
        plistener->setProgress(0.0);
    }


    const bool skip = (gridSize_ == RAWParams::BayerSensor::ePSMotionCorrection::Grid1x2);
    int gridSize = 1;
    bool nOf3x3 = false;

    switch (gridSize_) {
        case RAWParams::BayerSensor::ePSMotionCorrection::Grid1x1:
        case RAWParams::BayerSensor::ePSMotionCorrection::Grid1x2:
            gridSize = 1;
            break;

        case RAWParams::BayerSensor::ePSMotionCorrection::Grid3x3:
            gridSize = 3;
            break;

        case RAWParams::BayerSensor::ePSMotionCorrection::Grid5x5:
            gridSize = 5;
            break;

        case RAWParams::BayerSensor::ePSMotionCorrection::Grid7x7:
            gridSize = 7;
            break;

        case RAWParams::BayerSensor::ePSMotionCorrection::Grid3x3New:
            gridSize = 1;
            nOf3x3 = true;
    }

    // Lookup table for non adaptive (slider) mode
    LUTf log2Lut(32768, LUT_CLIP_BELOW | LUT_CLIP_ABOVE);

    if(detectMotion && !adaptive) {
        const float lutStrength = 2.f;
        log2Lut[0] = 0;

        for(int i = 2; i < 65536; i += 2) {
            log2Lut[i >> 1] = lutStrength * log2(i) / 100.f;
        }
    }

    const float scaleGreen = 1.f / scale_mul[1];

    float nRead;
    float eperIsoModel;

    int nReadIndex = static_cast<int>(round(log2(idata->getISOSpeed() /  100.f) * 3.f));

    if(model.find("K-3") != string::npos) {
        nRead = nReadK3II[nReadIndex];
        eperIsoModel = ePerIsoK3II;
    } else if(model.find("K-1") != string::npos) {
        nRead = nReadK1[nReadIndex];
        eperIsoModel = ePerIsoK1;
    } else {
        nRead = nReadK70[nReadIndex];
        eperIsoModel = ePerIsoK70;
    }

    nRead *= pow(2.f, nreadIso);
    eperIsoModel *= pow(2.f, eperIso);
    eperIso = eperIsoModel * (100.f / (rawWpCorrection * idata->getISOSpeed()));

    float eperIsoRed = eperIso / scale_mul[0];
    float eperIsoGreen = eperIso * scaleGreen;
    float eperIsoBlue = eperIso / scale_mul[2];

//    printf("Pixelshift parameters : gridSize %d\tadaptive %d\tstdDevFactorGreen %f\telectrons %1.8f\tnread %f\tprnu %1.1f%%\n", gridSize, adaptive, stddevFactorGreen, eperIso, nRead, prnu);

    prnu /= 100.f;
    stddevFactorGreen *= stddevFactorGreen;
    stddevFactorRed *= stddevFactorRed;
    stddevFactorBlue *= stddevFactorBlue;


    nRead *= nRead;

    // If the values of two corresponding green pixels differ my more then motionThreshold %, the pixel will be treated as a badGreen pixel
    float motionThreshold = 1.f - (motion / 100.f);
    // For shades of green motion indicators
    const float blendFactor = ((adaptive || motion == 0.f) ? 1.f : 1.f / (1.f - motionThreshold));

    unsigned int offsX = 0, offsY = 0;

    if(!bayerParams.pixelShiftMedian || !adaptive) {
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
    }

    const float thresh = adaptive ? 0.f : motionThreshold;
    array2D<float> psRed(winw + 32, winh); // increase width to avoid cache conflicts
    array2D<float> psG1(winw + 32, winh);
    array2D<float> psG2(winw + 32, winh);
    array2D<float> psBlue(winw + 32, winh);

    array2D<float> psMask(winw, winh, ARRAY2D_CLEAR_DATA);
// fill channels psRed, psG1, psG2 and psBlue
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif

    for(int i = winy + 1; i < winh - 1; ++i) {
        float *greenDest1 = psG1[i];
        float *greenDest2 = psG2[i];
        float *nonGreenDest0 = psRed[i];
        float *nonGreenDest1 = psBlue[i];
        int j = winx + 1;
        int c = FC(i, j);

        if (c == 2 || ((c & 1) && FC(i, j + 1) == 2)) {
            // row with blue pixels => swap destination pointers for non green pixels
            std::swap(nonGreenDest0, nonGreenDest1);
            std::swap(greenDest1, greenDest2);
        }

        // offset to keep the code short. It changes its value between 0 and 1 for each iteration of the loop
        unsigned int offset = (c & 1);
        offset ^= 1; // 0 => 1 or 1 => 0

        for(; j < winw - 1; ++j) {
            offset ^= 1; // 0 => 1 or 1 => 0

            // store the values from the 4 frames into 4 different temporary planes
            greenDest1[j] = (*rawDataFrames[1 - offset])[i - offset + 1][j];
            greenDest2[j] = (*rawDataFrames[3 - offset])[i + offset][j + 1];
            nonGreenDest0[j] = (*rawDataFrames[(offset << 1) + offset])[i][j + offset];
            nonGreenDest1[j] = (*rawDataFrames[2 - offset])[i + 1][j - offset + 1];
        }
    }

// now that the temporary planes are filled for easy access we do the motion detection
    int sum0 = 0;
    int sum1 = 0;
    float pixelcount = ((winh - (border + offsY) - (winy + border - offsY)) * (winw - (border + offsX) - (winx + border - offsX))) / 2.f;
#ifdef _OPENMP
    #pragma omp parallel for reduction(+:sum0,sum1) schedule(dynamic,16)
#endif

    for(int i = winy + border - offsY; i < winh - (border + offsY); ++i) {
        float *greenDest = green[i + offsY];
        float *redDest = red[i + offsY];
        float *blueDest = blue[i + offsY];
        int j = winx + border - offsX;

        float greenDifMax[gridSize]; // Here we store the maximum differences per Column

        // green channel motion detection checks the grid around the pixel for differences in green channels
        if(detectMotion || (adaptive && checkGreen)) {
            if(gridSize == 3) {
                // compute maximum of differences for first two columns of 3x3 grid
                greenDifMax[0] =  std::max({greenDiff(psG1[i - 1][j - 1], psG2[i - 1][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[ i ][j - 1], psG2[ i ][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 1][j - 1], psG2[i + 1][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[1] =  std::max({greenDiff(psG1[i - 1][ j ], psG2[i - 1][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[ i ][ j ], psG2[ i ][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 1][ j ], psG2[i + 1][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
            } else if(gridSize == 5) {
                // compute maximum of differences for first four columns of 5x5 grid
                greenDifMax[0] =  std::max({greenDiff(psG1[i - 2][j - 2], psG2[i - 2][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i - 1][j - 2], psG2[i - 1][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[ i ][j - 2], psG2[ i ][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 1][j - 2], psG2[i + 1][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 2][j - 2], psG2[i + 2][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[1] =  std::max({greenDiff(psG1[i - 2][j - 1], psG2[i - 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i - 1][j - 1], psG2[i - 1][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[ i ][j - 1], psG2[ i ][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 1][j - 1], psG2[i + 1][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 2][j - 1], psG2[i + 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[2] =  std::max({greenDiff(psG1[i - 2][ j ], psG2[i - 2][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i - 1][ j ], psG2[i - 1][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[ i ][ j ], psG2[ i ][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 1][ j ], psG2[i + 1][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 2][ j ], psG2[i + 2][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[3] =  std::max({greenDiff(psG1[i - 2][j + 1], psG2[i - 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i - 1][j + 1], psG2[i - 1][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[ i ][j + 1], psG2[ i ][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 1][j + 1], psG2[i + 1][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 2][j + 1], psG2[i + 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
            } else if(gridSize == 7) {
                // compute maximum of differences for first six columns of 7x7 grid
                greenDifMax[0] =  std::max({greenDiff(psG1[i - 3][j - 3], psG2[i - 3][j - 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i - 2][j - 3], psG2[i - 2][j - 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i - 1][j - 3], psG2[i - 1][j - 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[ i ][j - 3], psG2[ i ][j - 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 1][j - 3], psG2[i + 1][j - 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 2][j - 3], psG2[i + 2][j - 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 3][j - 3], psG2[i + 3][j - 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[1] =  std::max({greenDiff(psG1[i - 3][j - 2], psG2[i - 3][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i - 2][j - 2], psG2[i - 2][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i - 1][j - 2], psG2[i - 1][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[ i ][j - 2], psG2[ i ][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 1][j - 2], psG2[i + 1][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 2][j - 2], psG2[i + 2][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 3][j - 2], psG2[i + 3][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[2] =  std::max({greenDiff(psG1[i - 3][j - 1], psG2[i - 3][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i - 2][j - 1], psG2[i - 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i - 1][j - 1], psG2[i - 1][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[ i ][j - 1], psG2[ i ][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 1][j - 1], psG2[i + 1][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 2][j - 1], psG2[i + 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 3][j - 1], psG2[i + 3][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[3] =  std::max({greenDiff(psG1[i - 3][ j ], psG2[i - 3][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i - 2][ j ], psG2[i - 2][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i - 1][ j ], psG2[i - 1][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[ i ][ j ], psG2[ i ][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 1][ j ], psG2[i + 1][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 2][ j ], psG2[i + 2][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 3][ j ], psG2[i + 3][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[4] =  std::max({greenDiff(psG1[i - 3][j + 1], psG2[i - 3][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i - 2][j + 1], psG2[i - 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i - 1][j + 1], psG2[i - 1][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[ i ][j + 1], psG2[ i ][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 1][j + 1], psG2[i + 1][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 2][j + 1], psG2[i + 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 3][j + 1], psG2[i + 3][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
                greenDifMax[5] =  std::max({greenDiff(psG1[i - 3][j + 2], psG2[i - 3][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i - 2][j + 2], psG2[i - 2][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i - 1][j + 2], psG2[i - 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[ i ][j + 2], psG2[ i ][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 1][j + 2], psG2[i + 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 2][j + 2], psG2[i + 2][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                            greenDiff(psG1[i + 3][j + 2], psG2[i + 3][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                           });
            }

        }


        // this is the index for the last column of the grid. Obviously we have to start with gridSize - 1
        int lastIndex = gridSize - 1;
        float korr = 0.f;
        int c = FC(i, j);
        bool blueRow = false;

        if (c == 2 || ((c & 1) && FC(i, j + 1) == 2)) {
            // row with blue pixels => swap destination pointers for non green pixels
            blueRow = true;
        }

        // offset to keep the code short. It changes its value between 0 and 1 for each iteration of the loop
        unsigned int offset = (c & 1);
//        offset ^= 1; // 0 => 1 or 1 => 0

        for(; j < winw - (border + offsX); ++j) {
            bool greenFromPs = false;
            offset ^= 1; // 0 => 1 or 1 => 0

            if(detectMotion || (adaptive && checkGreen)) {
                bool skipNext = false;
                float gridMax;

                if(gridSize < 2) {
                    // compute difference for current pixel and skip next pixel, that's roughly the method from dcrawps
                    gridMax = greenDiff(psG1[i][j], psG2[i][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i);
                    skipNext = skip;
                } else if(gridSize == 3) {
                    // compute maximum of differences for third column of 3x3 grid and save at position lastIndex
                    greenDifMax[lastIndex] = std::max({greenDiff(psG1[i - 1][j + 1], psG2[i - 1][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff(psG1[ i ][j + 1], psG2[ i ][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff(psG1[i + 1][j + 1], psG2[i + 1][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                                      });
                    // calculate maximum of whole grid by calculating maximum of grid column max values
                    gridMax = std::max({greenDifMax[0], greenDifMax[1], greenDifMax[2]});
                } else if(gridSize == 5) {
                    // compute maximum of differences for fifth column of 5x5 grid and save at position lastIndex
                    greenDifMax[lastIndex] = std::max({greenDiff(psG1[i - 2][j + 2], psG2[i - 2][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff(psG1[i - 1][j + 2], psG2[i - 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff(psG1[ i ][j + 2], psG2[ i ][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff(psG1[i + 1][j + 2], psG2[i + 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff(psG1[i + 2][j + 2], psG2[i + 2][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i)
                                                      });
                    // calculate maximum of whole grid by calculating maximum of grid column max values
                    gridMax = std::max({greenDifMax[0], greenDifMax[1], greenDifMax[2], greenDifMax[3], greenDifMax[4]});
                } else if(gridSize == 7) {
                    // compute maximum of differences for 7th column of 7x7 grid and save at position lastIndex
                    greenDifMax[lastIndex] = std::max({greenDiff(psG1[i - 3][j + 3], psG2[i - 3][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff(psG1[i - 2][j + 3], psG2[i - 2][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff(psG1[i - 1][j + 3], psG2[i - 1][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff(psG1[ i ][j + 3], psG2[ i ][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff(psG1[i + 1][j + 3], psG2[i + 1][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff(psG1[i + 2][j + 3], psG2[i + 2][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                       greenDiff(psG1[i + 3][j + 3], psG2[i + 3][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion, j, i),
                                                      });
                    // calculate maximum of whole grid by calculating maximum of grid column max values
                    gridMax = std::max({greenDifMax[0], greenDifMax[1], greenDifMax[2], greenDifMax[3], greenDifMax[4], greenDifMax[5], greenDifMax[6]});
                }


                // adjust index for next column
                lastIndex ++;
                lastIndex = lastIndex == gridSize ? 0 : lastIndex;

                // increase motion detection dependent on brightness
                if(!adaptive) {
                    korr = log2Lut[((int)(psG1[i][j] * scaleGreen)) >> 1];
                }

                if (gridMax > thresh - korr) {
                    if(offset == 0) {
                        sum0 ++;
                    } else {
                        sum1 ++;
                    }

                    if(nOf3x3) {
                        psMask[i][j] = 1.f;
                    } else if((offset == (frame & 1)) && checkNonGreenVertical) {
                        if(frame > 1) {
                            green[i + offsY][j + offsX] = blueRow ? psG1[i][j] : psG2[i][j];
                        } else {
                            green[i + offsY][j + offsX] = blueRow ? psG2[i][j] : psG1[i][j];;
                        }

                        continue;
                    } else {
                        // at least one of the tested green pixels of the grid is detected as motion
                        paintMotionMask(j + offsX, showMotion, (gridMax - thresh + korr) * blendFactor, showOnlyMask, greenDest, redDest, blueDest);

                        if(skipNext) {
                            // treat the horizontally next pixel also as motion
                            j++;
                            paintMotionMask(j + offsX, showMotion, (gridMax - thresh + korr) * blendFactor, showOnlyMask, greenDest, redDest, blueDest);
                        }
                        // do not set the motion pixel values. They have already been set by demosaicer or showMotion
                        continue;
                    }
                }
            }

            if(adaptive && checkNonGreenCross) {
                // check red cross
                float redTop    = psRed[i - 1][ j ];
                float redLeft   = psRed[ i ][j - 1];
                float redCentre = psRed[ i ][ j ];
                float redRight  = psRed[ i ][j + 1];
                float redBottom = psRed[i + 1][ j ];
                float redDiff = nonGreenDiffCross(redRight, redLeft, redTop, redBottom, redCentre, stddevFactorRed, eperIsoRed, nRead, prnu, showMotion);

                if(redDiff > 0.f) {
                    if(nOf3x3) {
                        psMask[i][j] = redBlueWeight;
                    } else {
                        paintMotionMask(j + offsX, showMotion, redDiff, showOnlyMask, redDest, blueDest, greenDest);
                    }

                    continue;
                }

                // check blue cross
                float blueTop    = psBlue[i - 1][ j ];
                float blueLeft   = psBlue[ i ][j - 1];
                float blueCentre = psBlue[ i ][ j ];
                float blueRight  = psBlue[ i ][j + 1];
                float blueBottom = psBlue[i + 1][ j ];
                float blueDiff = nonGreenDiffCross(blueRight, blueLeft, blueTop, blueBottom, blueCentre, stddevFactorBlue, eperIsoBlue, nRead, prnu, showMotion);

                if(blueDiff > 0.f) {
                    if(nOf3x3) {
                        psMask[i][j] = redBlueWeight;
                    } else {
                        paintMotionMask(j + offsX, showMotion, blueDiff, showOnlyMask, blueDest, redDest, greenDest);
                    }

                    continue;

                }
            }

            if(adaptive && checkNonGreenHorizontal) {
                float redLeft    = psRed[ i ][j - 1];
                float redCentre  = psRed[ i ][ j ];
                float redRight   = psRed[ i ][j + 1];

                float redDiffLeft = redLeft - redCentre;
                float redDiffRight = redRight - redCentre;

                if(redDiffLeft * redDiffRight >= 0.f) {
                    float redAvg = (redRight + redLeft) / 2.f;
                    float redDiffHor = nonGreenDiff(redCentre, redAvg, stddevFactorRed, eperIsoRed, nRead, prnu, showMotion);

                    if(redDiffHor > 0.f) {
                        if(nOf3x3) {
                            psMask[i][j] = redBlueWeight;
                        } else {
                            paintMotionMask(j + offsX, showMotion, redDiffHor, showOnlyMask, redDest, blueDest, greenDest);
                        }

                        continue;
                    }
                }

                float blueLeft    = psBlue[ i ][j - 1];
                float blueCentre  = psBlue[ i ][ j ];
                float blueRight   = psBlue[ i ][j + 1];

                float blueDiffLeft = blueLeft - blueCentre;
                float blueDiffRight = blueRight - blueCentre;

                if(blueDiffLeft * blueDiffRight >= 0.f) {
                    float blueAvg = (blueRight + blueLeft) / 2.f;
                    float blueDiffHor = nonGreenDiff(blueCentre, blueAvg, stddevFactorBlue, eperIsoBlue, nRead, prnu, showMotion);

                    if(blueDiffHor > 0.f) {
                        if(nOf3x3) {
                            psMask[i][j] = redBlueWeight;
                        } else {
                            paintMotionMask(j + offsX, showMotion, blueDiffHor, showOnlyMask, blueDest, redDest, greenDest);
                        }

                        continue;
                    }
                }
            }

            if(adaptive && checkNonGreenVertical) {
                // check red vertically
                float redTop    = psRed[i - 1][ j ];
                float redCentre = psRed[ i ][ j ];
                float redBottom = psRed[i + 1][ j ];

                float redDiffTop = redTop - redCentre;
                float redDiffBottom = redBottom - redCentre;

                if(redDiffTop * redDiffBottom >= 0.f) {
                    float redAvg = (redTop + redBottom) / 2.f;
                    float redDiff = nonGreenDiff(redCentre, redAvg, stddevFactorRed, eperIsoRed, nRead, prnu, showMotion);

                    if(redDiff > 0.f) {
                        if(nOf3x3) {
                            psMask[i][j] = redBlueWeight;
                        } else {
                            paintMotionMask(j + offsX, showMotion, redDiff, showOnlyMask, redDest, blueDest, greenDest);
                        }

                        continue;
                    }
                }

                // check blue vertically
                float blueTop    = psBlue[i - 1][ j ];
                float blueCentre = psBlue[ i ][ j ];
                float blueBottom = psBlue[i + 1][ j ];

                float blueDiffTop = blueTop - blueCentre;
                float blueDiffBottom = blueBottom - blueCentre;

                if(blueDiffTop * blueDiffBottom >= 0.f) {
                    float blueAvg = (blueTop + blueBottom) / 2.f;
                    float blueDiff = nonGreenDiff(blueCentre, blueAvg, stddevFactorBlue, eperIsoBlue, nRead, prnu, showMotion);

                    if(blueDiff > 0.f) {
                        if(nOf3x3) {
                            psMask[i][j] = redBlueWeight;
                        } else {
                            paintMotionMask(j + offsX, showMotion, blueDiff, showOnlyMask, blueDest, redDest, greenDest);
                        }

                        continue;
                    }
                }
            }

            if(adaptive && checkNonGreenAmaze) {
                // check current pixel against amaze
                float redCentre  = psRed[ i ][ j ];
                float redAmaze   = red[i + offsY][j + offsX];

                float redDiffAmaze = nonGreenDiff(redCentre, redAmaze, stddevFactorRed, eperIsoRed, nRead, prnu, showMotion);

                if(redDiffAmaze > 0.f) {
                    if(nOf3x3) {
                        psMask[i][j] = redBlueWeight;
                    } else {
                        paintMotionMask(j + offsX, showMotion, redDiffAmaze, showOnlyMask, redDest, blueDest, greenDest);
                    }

                    continue;
                }

                float blueCentre  = psBlue[ i ][ j ];
                float blueAmaze   = blue[i + offsY][j + offsX];

                float blueDiffAmaze = nonGreenDiff(blueCentre, blueAmaze, stddevFactorBlue, eperIsoBlue, nRead, prnu, showMotion);

                if(blueDiffAmaze > 0.f) {
                    if(nOf3x3) {
                        psMask[i][j] = redBlueWeight;
                    } else {
                        paintMotionMask(j + offsX, showMotion, blueDiffAmaze, showOnlyMask, blueDest, redDest, greenDest);
                    }

                    continue;
                }
            }

            if(adaptive && checkNonGreenCross2) { // for green amaze
                float greenCentre = (psG1[ i ][ j ] + psG2[ i ][ j ]) / 2.f;
                float greenAmaze = green[i + offsY][j + offsX];
                float greenDiffAmaze = nonGreenDiff(greenCentre, greenAmaze, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion);

                if(greenDiffAmaze > 0.f) {
                    if(nOf3x3) {
                        psMask[i][j] = 1.f;
                    } else {
                        paintMotionMask(j + offsX, showMotion, greenDiffAmaze, showOnlyMask, greenDest, redDest, blueDest);
                    }

                    continue;
                }
            }

            if(adaptive && experimental0) { // for experiments
//                float green1Median, green2Median;
//                green1Median = median(psG1[ i - 1 ][ j - 1 ],psG1[ i - 1 ][ j + 1 ],psG1[ i ][ j ],psG1[ i + 1 ][ j -1 ],psG1[ i + 1 ][ j + 1 ]);
//                green2Median = median(psG2[ i - 1 ][ j - 1 ],psG2[ i - 1 ][ j + 1 ],psG2[ i ][ j ],psG2[ i + 1 ][ j -1 ],psG2[ i + 1 ][ j + 1 ]);
//                float greenDiffMedian = nonGreenDiff(green1Median, green2Median, stddevFactorGreen * 0.36f, eperIsoGreen, nRead, prnu, showMotion);
//
//                if(greenDiffMedian > 0.f) {
//                    if(nOf3x3) {
//                        psMask[i][j] = 1.f;
//                    } else {
//                        paintMotionMask(j + offsX, showMotion, greenDiffMedian, showOnlyMask, greenDest, redDest, blueDest);
//                    }
//
//                    continue;
//                }

            }

            if(showMotion && showOnlyMask) { // we want only motion mask => paint areas without motion in pure black
                red[i + offsY][j + offsX] = green[i + offsY][j + offsX] = blue[i + offsY][j + offsX] = 0.f;
            } else if(!(adaptive && nOf3x3)) {
                // no motion detected, replace the a priori demosaiced values by the pixelshift combined values
                red[i + offsY][j + offsX] = psRed[i][j];
                green[i + offsY][j + offsX] = (psG1[i][j] + psG2[i][j]) / 2.f;
                blue[i + offsY][j + offsX] = psBlue[i][j];
            }
        }
    }

    float percent0 = 100.f * sum0 / pixelcount;
    float percent1 = 100.f * sum1 / pixelcount;

    std::cout << fileName <<  " : Green detections at stddev " << std::setprecision( 2 ) << bayerParams.pixelShiftStddevFactorGreen << " : Frame 1/3 : " << std::setprecision( 6 ) << sum0 << " (" << percent0 << "%)" << " Frame 2/4 : " << sum1 << " (" << percent1 << "%)" << std::endl;

    if(adaptive && nOf3x3) {
        if(blurMap) {
            #pragma omp parallel
            {
                gaussianBlur(psMask, psMask, winw, winh, sigma);
            }
        }

        array2D<uint8_t> mask(W, H, ARRAY2D_CLEAR_DATA);
        array2D<uint8_t> maskInv(W, H, ARRAY2D_CLEAR_DATA);
        #pragma omp parallel for schedule(dynamic,16)

        for(int i = winy + border - offsY; i < winh - (border + offsY); ++i) {
            float *greenDest = green[i + offsY];
            float *redDest = red[i + offsY];
            float *blueDest = blue[i + offsY];
            int j = winx + border - offsX;
            float v3sum[3] = {0.f};

            for(int v = -1; v <= 1; v++) {
                for(int h = -1; h < 1; h++) {
                    v3sum[1 + h] += psMask[i + v][j + h];
                }
            }

            float blocksum = v3sum[0] + v3sum[1];

            for(int voffset = 2; j < winw - (border + offsX); ++j, ++voffset) {
                float colSum = psMask[i - 1][j + 1] + psMask[i][j + 1] + psMask[i + 1][j + 1];
                voffset = voffset == 3 ? 0 : voffset;  // faster than voffset %= 3;
                blocksum -= v3sum[voffset];
                blocksum += colSum;
                v3sum[voffset] = colSum;

                if(blocksum >= threshold) {
                    mask[i][j] = 255;
                }
            }
        }

        if(holeFill) {
            invertMask(winx + border - offsX, winw - (border + offsX), winy + border - offsY, winh - (border + offsY), mask, maskInv);
            floodFill4(winx + border - offsX, winw - (border + offsX), winy + border - offsY, winh - (border + offsY), maskInv);
            xorMasks(winx + border - offsX, winw - (border + offsX), winy + border - offsY, winh - (border + offsY), maskInv, mask);
        }

        #pragma omp parallel for schedule(dynamic,16)

        for(int i = winy + border - offsY; i < winh - (border + offsY); ++i) {
            float *greenDest = green[i + offsY];
            float *redDest = red[i + offsY];
            float *blueDest = blue[i + offsY];

            for(int j = winx + border - offsX; j < winw - (border + offsX); ++j) {
                if(mask[i][j] == 255 ) {
                    paintMotionMask(j + offsX, showMotion, 0.5f, showOnlyMask, greenDest, redDest, blueDest);
                    continue;
                }

                if(showMotion && showOnlyMask) { // we want only motion mask => paint areas without motion in pure black
                    red[i + offsY][j + offsX] = green[i + offsY][j + offsX] = blue[i + offsY][j + offsX] = 0.f;
                } else {
                    if(blurMap && experimental0) {
                        red[i + offsY][j + offsX] = intp(psMask[i][j], red[i + offsY][j + offsX], psRed[i][j] );
                        green[i + offsY][j + offsX] = intp(psMask[i][j],green[i + offsY][j + offsX],(psG1[i][j] + psG2[i][j]) / 2.f);
                        blue[i + offsY][j + offsX] = intp(psMask[i][j],blue[i + offsY][j + offsX], psBlue[i][j]);
                    } else {
                        red[i + offsY][j + offsX] = psRed[i][j];
                        green[i + offsY][j + offsX] = (psG1[i][j] + psG2[i][j]) / 2.f;
                        blue[i + offsY][j + offsX] = psBlue[i][j];
                    }
                }
            }
        }
    }

    if(plistener) {
        plistener->setProgress(1.0);
    }
}
#endif