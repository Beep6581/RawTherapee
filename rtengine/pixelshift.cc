////////////////////////////////////////////////////////////////
//
//  pentax pixelshift algorithm with motion detection
//
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

float greenDiff(float a, float b, bool adaptive, float stddevFactor, float eperIso, float nreadIso, float prnu, bool showMotion)
{
    // calculate the difference between two green samples
#ifdef PIXELSHIFTDEV
    if(adaptive) {
#endif
        float gDiff = a - b;
        gDiff *= eperIso;
        gDiff *= gDiff;
        float avg = (a + b) * 0.5f;
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
#ifdef PIXELSHIFTDEV

    } else {
        float gDiff = std::fabs(a - b);
        // add a small epsilon to avoid division by zero
        float maxVal = std::max(a, b) + 0.01f;
        return gDiff / maxVal;
    }
#endif
}

#ifdef PIXELSHIFTDEV
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
#endif

float nonGreenDiffCross(float right, float left, float top, float bottom, float centre, float clippedVal, float stddevFactor, float eperIso, float nreadIso, float prnu, bool showMotion)
{
    if(rtengine::max(right, left, top, bottom, centre) > clippedVal) {
        return 0.f;
    }
    // check non green cross
    float hDiff = (right + left) * 0.5f - centre;
    hDiff *= eperIso;
    hDiff *= hDiff;
    float vDiff = (top + bottom) * 0.5f - centre;
    vDiff *= eperIso;
    vDiff *= vDiff;
    float avg = ((right + left) + (top + bottom)) * 0.25f;
    avg *= eperIso;
    prnu *= avg;
    float stddev = stddevFactor * (avg + nreadIso + prnu * prnu);
    float result = std::min(hDiff, vDiff) - stddev;

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

void invertMask(int xStart, int xEnd, int yStart, int yEnd, const array2D<uint8_t> &maskIn, array2D<uint8_t> &maskOut)
{
    #pragma omp parallel for schedule(dynamic,16)

    for(int i = yStart; i < yEnd; ++i) {
        #pragma omp simd

        for(int j = xStart; j < xEnd; ++j) {
            maskOut[i][j] = ~maskIn[i][j];
        }
    }
}

void xorMasks(int xStart, int xEnd, int yStart, int yEnd, const array2D<uint8_t> &maskIn, array2D<uint8_t> &maskOut)
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

        for(uint16_t i = yStart; i < yEnd; i++)
        {
            floodFill4Impl(i, xStart, xStart, xEnd, yStart, yEnd, mask, coordStack);
        }

        #pragma omp for schedule(dynamic,128) nowait

        for(int16_t i = yEnd - 1; i >= 0 ; i--)
        {
            floodFill4Impl(i, xEnd - 1, xStart, xEnd, yStart, yEnd, mask, coordStack);
        }

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
void RawImageSource::pixelshift(int winx, int winy, int winw, int winh, const RAWParams::BayerSensor &bayerParamsIn, unsigned int frame, const std::string &model, float rawWpCorrection)
{
#ifdef PIXELSHIFTDEV
    BENCHFUN
#endif

    if(numFrames != 4) { // fallback for non pixelshift files
        amaze_demosaic_RT (0, 0, winw, winh, rawData, red, green, blue);
        return;
    }

    RAWParams::BayerSensor bayerParams = bayerParamsIn;

#ifndef PIXELSHIFTDEV
    bayerParams.pixelShiftAutomatic = true;
#endif

    if(bayerParams.pixelShiftMotionCorrectionMethod == RAWParams::BayerSensor::Automatic) {
        bool pixelShiftEqualBright = bayerParams.pixelShiftEqualBright;
        bayerParams.setPixelShiftDefaults();
        bayerParams.pixelShiftEqualBright = pixelShiftEqualBright;
    } else if(bayerParams.pixelShiftMotionCorrectionMethod == RAWParams::BayerSensor::Off) {
        bayerParams.pixelShiftMotion = 0;
        bayerParams.pixelShiftAutomatic = false;
        bayerParams.pixelshiftShowMotion = false;
    }

    if((bayerParams.pixelShiftMotion > 0 || bayerParams.pixelShiftAutomatic)) {
        if(bayerParams.pixelShiftMedian) { // We need the amaze demosaiced frames for motion correction
#ifdef PIXELSHIFTDEV
            if(!bayerParams.pixelShiftMedian3) {
#endif
                if(bayerParams.pixelShiftLmmse) {
                    lmmse_interpolate_omp(winw, winh, *(rawDataFrames[0]), red, green, blue, bayerParams.lmmse_iterations);
                } else {
                    amaze_demosaic_RT (0, 0, winw, winh, *(rawDataFrames[0]), red, green, blue);
                }
                multi_array2D<float,3> redTmp(W,H);
                multi_array2D<float,3> greenTmp(W,H);
                multi_array2D<float,3> blueTmp(W,H);
                for(int i=0;i<3;i++) {
                    if(bayerParams.pixelShiftLmmse) {
                        lmmse_interpolate_omp(winw, winh, *(rawDataFrames[i+1]), redTmp[i], greenTmp[i], blueTmp[i], bayerParams.lmmse_iterations);
                    } else {
                        amaze_demosaic_RT (0, 0, winw, winh, *(rawDataFrames[i+1]), redTmp[i], greenTmp[i], blueTmp[i]);
                    }
                }
                #pragma omp parallel for schedule(dynamic,16)
                for(int i=border;i<H-border;i++) {
                    for(int j=border;j<W-border;j++) {
                        red[i][j] = median(red[i][j],redTmp[0][i+1][j],redTmp[1][i+1][j+1],redTmp[2][i][j+1]);
                    }
                    for(int j=border;j<W-border;j++) {
                        green[i][j] = median(green[i][j],greenTmp[0][i+1][j],greenTmp[1][i+1][j+1],greenTmp[2][i][j+1]);
                    }
                    for(int j=border;j<W-border;j++) {
                        blue[i][j] = median(blue[i][j],blueTmp[0][i+1][j],blueTmp[1][i+1][j+1],blueTmp[2][i][j+1]);
                    }
                }
#ifdef PIXELSHIFTDEV
            } else {
                multi_array2D<float,3> redTmp(W,H);
                multi_array2D<float,3> greenTmp(W,H);
                multi_array2D<float,3> blueTmp(W,H);
                for(int i=0, frameIndex = 0;i<4;++i) {
                    if(i != currFrame) {
                        if(bayerParams.pixelShiftLmmse) {
                            lmmse_interpolate_omp(winw, winh, *(rawDataFrames[i]), redTmp[frameIndex], greenTmp[frameIndex], blueTmp[frameIndex], bayerParams.lmmse_iterations);
                        } else {
                            amaze_demosaic_RT (0, 0, winw, winh, *(rawDataFrames[i]), redTmp[frameIndex], greenTmp[frameIndex], blueTmp[frameIndex]);
                        }
                        ++frameIndex;
                    }
                }
                unsigned int offsX0 = 0, offsY0 = 0;
                unsigned int offsX1 = 0, offsY1 = 0;
                unsigned int offsX2 = 0, offsY2 = 0;

                // We have to adjust the offsets for the selected subframe we exclude from median
                switch (currFrame) {
                    case 0:
                        offsY0 = 1;
                        offsX0 = 0;
                        offsY1 = 1;
                        offsX1 = 1;
                        offsY2 = 0;
                        offsX2 = 1;
                        break;

                    case 1:
                        offsY0 = 0;
                        offsX0 = 0;
                        offsY1 = 1;
                        offsX1 = 1;
                        offsY2 = 0;
                        offsX2 = 1;
                        break;

                    case 2:
                        offsY0 = 0;
                        offsX0 = 0;
                        offsY1 = 1;
                        offsX1 = 0;
                        offsY2 = 0;
                        offsX2 = 1;
                        break;

                    case 3:
                        offsY0 = 0;
                        offsX0 = 0;
                        offsY1 = 1;
                        offsX1 = 0;
                        offsY2 = 1;
                        offsX2 = 1;
                }

                #pragma omp parallel for schedule(dynamic,16)
                for(int i=border;i<H-border;i++) {
                    for(int j=border;j<W-border;j++) {
                        red[i][j] = median(redTmp[0][i+offsY0][j+offsX0],redTmp[1][i+offsY1][j+offsX1],redTmp[2][i+offsY2][j+offsX2]);
                    }
                    for(int j=border;j<W-border;j++) {
                        green[i][j] = median(greenTmp[0][i+offsY0][j+offsX0],greenTmp[1][i+offsY1][j+offsX1],greenTmp[2][i+offsY2][j+offsX2]);
                    }
                    for(int j=border;j<W-border;j++) {
                        blue[i][j] = median(blueTmp[0][i+offsY0][j+offsX0],blueTmp[1][i+offsY1][j+offsX1],blueTmp[2][i+offsY2][j+offsX2]);
                    }
                }
            }
#endif
        } else {
            if(bayerParams.pixelShiftLmmse) {
                lmmse_interpolate_omp(winw, winh, rawData, red, green, blue, bayerParams.lmmse_iterations);
            } else {
                amaze_demosaic_RT (0, 0, winw, winh, rawData, red, green, blue);
            }
        }
    } else if(bayerParams.pixelShiftMotionCorrectionMethod != RAWParams::BayerSensor::Off) {
        if(bayerParams.pixelShiftLmmse) {
            lmmse_interpolate_omp(winw, winh, rawData, red, green, blue, bayerParams.lmmse_iterations);
        } else {
            amaze_demosaic_RT (0, 0, winw, winh, rawData, red, green, blue);
        }
    }

    const int motion = bayerParams.pixelShiftMotion;
    const bool showMotion = bayerParams.pixelshiftShowMotion;
    const bool showOnlyMask = bayerParams.pixelshiftShowMotionMaskOnly && showMotion;
    const RAWParams::BayerSensor::ePSMotionCorrection gridSize_ = bayerParams.pixelShiftMotionCorrection;
    const bool adaptive = bayerParams.pixelShiftAutomatic;
#ifdef PIXELSHIFTDEV
    const bool detectMotion = bayerParams.pixelShiftMotion > 0;
    float stddevFactorGreen = bayerParams.pixelShiftStddevFactorGreen;
    float stddevFactorRed = bayerParams.pixelShiftStddevFactorRed;
    float stddevFactorBlue = bayerParams.pixelShiftStddevFactorBlue;
    float nreadIso =  bayerParams.pixelShiftNreadIso;
    float prnu = bayerParams.pixelShiftPrnu;
    const float redBlueWeight = bayerParams.pixelShiftRedBlueWeight + 1.f;
#else
    float stddevFactorGreen = 5.f;
    float stddevFactorRed = 5.f;
    float stddevFactorBlue = 5.f;
    float nreadIso =  0.f;
    float prnu = 1.f;
    const float redBlueWeight = 0.7f + 1.f;
#endif
    float eperIso = bayerParams.pixelShiftEperIso;
    const bool checkNonGreenHorizontal = bayerParams.pixelShiftNonGreenHorizontal;
    const bool checkNonGreenVertical = bayerParams.pixelShiftNonGreenVertical;
    const bool checkNonGreenCross = bayerParams.pixelShiftNonGreenCross;
    const bool checkNonGreenAmaze = bayerParams.pixelShiftNonGreenAmaze;
    const bool checkNonGreenCross2 = bayerParams.pixelShiftNonGreenCross2;
    const bool checkGreen = bayerParams.pixelShiftGreen;
    const float greenWeight = 2.f;
    const bool blurMap = bayerParams.pixelShiftBlur;
    const float sigma = bayerParams.pixelShiftSigma;
#ifdef PIXELSHIFTDEV
    const float threshold = bayerParams.pixelShiftSum + 9.f;
#else
    constexpr float threshold = 3.f + 9.f;
#endif
    const bool experimental0 = bayerParams.pixelShiftExp0;
    const bool holeFill = bayerParams.pixelShiftHoleFill;
    const bool equalBrightness = bayerParams.pixelShiftEqualBright;
    const bool smoothTransitions = blurMap && bayerParams.pixelShiftSmoothFactor > 0. && !showOnlyMask;
    const bool automatic = bayerParams.pixelShiftMotionCorrectionMethod == RAWParams::BayerSensor::Automatic;
    const float smoothFactor = 1.0 - bayerParams.pixelShiftSmoothFactor;

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
                                       2.4f,  // ISO 64000
                                       2.4f,  // ISO 80000
                                       2.4f,  // ISO 102400
                                       2.4f,  // ISO 128000
                                       2.4f,  // ISO 160000
                                       2.4f   // ISO 204800
                                   };

    static const float ePerIsoK1 = 0.75f;

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

    static const float ePerIsoK70 = 0.5f;

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

    if(adaptive && blurMap && nOf3x3 && smoothFactor == 0.f && !showMotion) {
        if(plistener) {
            plistener->setProgress(1.0);
        }

        return;
    }

#ifdef PIXELSHIFTDEV
    // Lookup table for non adaptive (slider) mode
    LUTf log2Lut(32768, LUT_CLIP_BELOW | LUT_CLIP_ABOVE);

    if(detectMotion && !adaptive) {
        const float lutStrength = 2.f;
        log2Lut[0] = 0;

        for(int i = 2; i < 65536; i += 2) {
            log2Lut[i >> 1] = lutStrength * log2(i) / 100.f;
        }
    }
#endif

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

#ifdef PIXELSHIFTDEV
    if(adaptive && experimental0) {
        eperIso = eperIsoModel * sqrtf(100.f / (rawWpCorrection * idata->getISOSpeed()));
    } else {
        eperIso = eperIsoModel * (100.f / (rawWpCorrection * idata->getISOSpeed()));
    }
#else
    eperIso = eperIsoModel * (100.f / (rawWpCorrection * idata->getISOSpeed()));
#endif

#ifdef PIXELSHIFTDEV
    std::cout << "WL: " << c_white[0] << " BL: " << c_black[0] << " ePerIso multiplicator: " << (65535.f / (c_white[0] - c_black[0])) << std::endl;
#endif
    const float eperIsoRed = (eperIso / scale_mul[0]) * (65535.f / (c_white[0] - c_black[0]));
    const float eperIsoGreen = (eperIso * scaleGreen) * (65535.f / (c_white[1] - c_black[1]));
    const float eperIsoBlue = (eperIso / scale_mul[2]) * (65535.f / (c_white[2] - c_black[2]));

    const float clippedRed = 65535.f / scale_mul[0];
    const float clippedBlue = 65535.f / scale_mul[2];

    prnu /= 100.f;
    stddevFactorGreen *= stddevFactorGreen;
    stddevFactorRed *= stddevFactorRed;
    stddevFactorBlue *= stddevFactorBlue;


    nRead *= nRead;

    // If the values of two corresponding green pixels differ my more then motionThreshold %, the pixel will be treated as a badGreen pixel
    const float motionThreshold = 1.f - (motion / 100.f);
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

    // calculate average green brightness for each frame
    float greenBrightness[4] = {1.f, 1.f, 1.f, 1.f};

    if(equalBrightness) {
        LUT<uint32_t> *histo[4];

        for(int i = 0; i < 4; ++i) {
            histo[i] = new LUT<uint32_t>(65536);
            histo[i]->clear();
        }

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            LUT<uint32_t> *histoThr[4];

            for(int i = 0; i < 4; ++i) {
                histoThr[i] = new LUT<uint32_t>(65536);
                histoThr[i]->clear();
            }

#ifdef _OPENMP
            #pragma omp for schedule(dynamic,16) nowait
#endif

            for(int i = winy + 1; i < winh - 1; ++i) {
                for(int j = winx + 1, offset = FC(i, j) & 1; j < winw - 1; ++j, offset ^= 1) {
                    (*histoThr[1 - offset])[(*rawDataFrames[1 - offset])[i - offset + 1][j]]++;
                    (*histoThr[3 - offset])[(*rawDataFrames[3 - offset])[i + offset][j + 1]]++;
                }
            }

            #pragma omp critical
            {
                for(int i = 0; i < 4; ++i) {
                    (*histo[i]) += (*histoThr[i]);
                    delete histoThr[i];
                }
            }
        }

        float medians[4];

        for(int i = 0; i < 4; ++i) {
            //find median of histogram
            uint32_t median = 0, count = 0;
            uint32_t datalen = (winh - 2) * (winw - 2) / 2;

            while (count < datalen / 2) {
                count += (*histo[i])[median];
                ++median;
            }

            const float weight = (count - datalen / 2.f) / (*histo[i])[median - 1];
            medians[i] = intp(weight, (float)(median - 2), (float)(median - 1));
            delete histo[i];
        }

        for(int i = 0; i < 4; ++i) {
            greenBrightness[i] = medians[frame] / medians[i];
        }

#ifdef PIXELSHIFTDEV
        std::cout << "brightness factors by median : " << greenBrightness[0] << " " << greenBrightness[1] << " " << greenBrightness[2] << " " << greenBrightness[3] << std::endl;
#endif

    }

    const float thresh = adaptive ? 0.f : motionThreshold;
    array2D<float> psRed(winw + 32, winh); // increase width to avoid cache conflicts
    array2D<float> psG1(winw + 32, winh);
    array2D<float> psG2(winw + 32, winh);
    array2D<float> psBlue(winw + 32, winh);

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

        if ((c + FC(i, j + 1)) == 3) {
            // row with blue pixels => swap destination pointers for non green pixels
            std::swap(nonGreenDest0, nonGreenDest1);
            std::swap(greenDest1, greenDest2);
        }

        // offset to keep the code short. It changes its value between 0 and 1 for each iteration of the loop
        unsigned int offset = c & 1;

        for(; j < winw - 1; ++j) {
            // store the values from the 4 frames into 4 different temporary planes
            greenDest1[j] = (*rawDataFrames[1 - offset])[i - offset + 1][j] * greenBrightness[1 - offset];
            greenDest2[j] = (*rawDataFrames[3 - offset])[i + offset][j + 1] * greenBrightness[3 - offset];
            nonGreenDest0[j] = (*rawDataFrames[(offset << 1) + offset])[i][j + offset] * greenBrightness[(offset << 1) + offset];
            nonGreenDest1[j] = (*rawDataFrames[2 - offset])[i + 1][j - offset + 1] * greenBrightness[2 - offset];
            offset ^= 1; // 0 => 1 or 1 => 0
        }
    }

// now that the temporary planes are filled for easy access we do the motion detection
#ifdef PIXELSHIFTDEV
    int sum[2] = {0};
    float pixelcount = ((winh - (border + offsY) - (winy + border - offsY)) * (winw - (border + offsX) - (winx + border - offsX))) / 2.f;
#endif

    array2D<float> psMask(winw, winh);

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
#ifdef PIXELSHIFTDEV
        int sumThr[2] = {0};
#endif
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16) nowait
#endif

        for(int i = winy + border - offsY; i < winh - (border + offsY); ++i) {
            float *greenDest = green[i + offsY];
            float *redDest = red[i + offsY];
            float *blueDest = blue[i + offsY];
            int j = winx + border - offsX;

#ifdef PIXELSHIFTDEV
            float greenDifMax[gridSize]; // Here we store the maximum differences per Column

            // green channel motion detection checks the grid around the pixel for differences in green channels

            if(detectMotion || (adaptive && checkGreen)) {
                if(gridSize == 3) {
                    // compute maximum of differences for first two columns of 3x3 grid
                    greenDifMax[0] =  std::max({greenDiff(psG1[i - 1][j - 1], psG2[i - 1][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[ i ][j - 1], psG2[ i ][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 1][j - 1], psG2[i + 1][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                               });
                    greenDifMax[1] =  std::max({greenDiff(psG1[i - 1][ j ], psG2[i - 1][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[ i ][ j ], psG2[ i ][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 1][ j ], psG2[i + 1][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                               });
                } else if(gridSize == 5) {
                    // compute maximum of differences for first four columns of 5x5 grid
                    greenDifMax[0] =  std::max({greenDiff(psG1[i - 2][j - 2], psG2[i - 2][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i - 1][j - 2], psG2[i - 1][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[ i ][j - 2], psG2[ i ][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 1][j - 2], psG2[i + 1][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 2][j - 2], psG2[i + 2][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                               });
                    greenDifMax[1] =  std::max({greenDiff(psG1[i - 2][j - 1], psG2[i - 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i - 1][j - 1], psG2[i - 1][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[ i ][j - 1], psG2[ i ][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 1][j - 1], psG2[i + 1][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 2][j - 1], psG2[i + 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                               });
                    greenDifMax[2] =  std::max({greenDiff(psG1[i - 2][ j ], psG2[i - 2][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i - 1][ j ], psG2[i - 1][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[ i ][ j ], psG2[ i ][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 1][ j ], psG2[i + 1][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 2][ j ], psG2[i + 2][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                               });
                    greenDifMax[3] =  std::max({greenDiff(psG1[i - 2][j + 1], psG2[i - 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i - 1][j + 1], psG2[i - 1][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[ i ][j + 1], psG2[ i ][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 1][j + 1], psG2[i + 1][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 2][j + 1], psG2[i + 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                               });
                } else if(gridSize == 7) {
                    // compute maximum of differences for first six columns of 7x7 grid
                    greenDifMax[0] =  std::max({greenDiff(psG1[i - 3][j - 3], psG2[i - 3][j - 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i - 2][j - 3], psG2[i - 2][j - 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i - 1][j - 3], psG2[i - 1][j - 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[ i ][j - 3], psG2[ i ][j - 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 1][j - 3], psG2[i + 1][j - 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 2][j - 3], psG2[i + 2][j - 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 3][j - 3], psG2[i + 3][j - 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                               });
                    greenDifMax[1] =  std::max({greenDiff(psG1[i - 3][j - 2], psG2[i - 3][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i - 2][j - 2], psG2[i - 2][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i - 1][j - 2], psG2[i - 1][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[ i ][j - 2], psG2[ i ][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 1][j - 2], psG2[i + 1][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 2][j - 2], psG2[i + 2][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 3][j - 2], psG2[i + 3][j - 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                               });
                    greenDifMax[2] =  std::max({greenDiff(psG1[i - 3][j - 1], psG2[i - 3][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i - 2][j - 1], psG2[i - 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i - 1][j - 1], psG2[i - 1][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[ i ][j - 1], psG2[ i ][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 1][j - 1], psG2[i + 1][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 2][j - 1], psG2[i + 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 3][j - 1], psG2[i + 3][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                               });
                    greenDifMax[3] =  std::max({greenDiff(psG1[i - 3][ j ], psG2[i - 3][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i - 2][ j ], psG2[i - 2][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i - 1][ j ], psG2[i - 1][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[ i ][ j ], psG2[ i ][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 1][ j ], psG2[i + 1][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 2][ j ], psG2[i + 2][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 3][ j ], psG2[i + 3][ j ], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                               });
                    greenDifMax[4] =  std::max({greenDiff(psG1[i - 3][j + 1], psG2[i - 3][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i - 2][j + 1], psG2[i - 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i - 1][j + 1], psG2[i - 1][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[ i ][j + 1], psG2[ i ][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 1][j + 1], psG2[i + 1][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 2][j + 1], psG2[i + 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 3][j + 1], psG2[i + 3][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                               });
                    greenDifMax[5] =  std::max({greenDiff(psG1[i - 3][j + 2], psG2[i - 3][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i - 2][j + 2], psG2[i - 2][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i - 1][j + 2], psG2[i - 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[ i ][j + 2], psG2[ i ][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 1][j + 2], psG2[i + 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 2][j + 2], psG2[i + 2][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                greenDiff(psG1[i + 3][j + 2], psG2[i + 3][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                               });
                }

            }

            // this is the index for the last column of the grid. Obviously we have to start with gridSize - 1
            int lastIndex = gridSize - 1;
            float korr = 0.f;
            bool blueRow = false;
#endif

            int c = FC(i, j);

#ifdef PIXELSHIFTDEV
            if (c == 2 || ((c & 1) && FC(i, j + 1) == 2)) {
                // row with blue pixels => swap destination pointers for non green pixels
                blueRow = true;
            }
#endif
            // offset to keep the code short. It changes its value between 0 and 1 for each iteration of the loop
            unsigned int offset = c & 1;

            for(; j < winw - (border + offsX); ++j) {
                psMask[i][j] = 1.f;

                offset ^= 1; // 0 => 1 or 1 => 0

#ifdef PIXELSHIFTDEV
                if(detectMotion || (adaptive && checkGreen)) {
                    bool skipNext = false;
#else
                if(adaptive && checkGreen) {
#endif
                    float gridMax;

#ifdef PIXELSHIFTDEV

                    if(gridSize < 2) {
                        // compute difference for current pixel and skip next pixel, that's roughly the method from dcrawps
#endif
                        gridMax = greenDiff(psG1[i][j], psG2[i][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion);
#ifdef PIXELSHIFTDEV

                        skipNext = skip;
                    } else if(gridSize == 3) {
                        // compute maximum of differences for third column of 3x3 grid and save at position lastIndex
                        greenDifMax[lastIndex] = std::max({greenDiff(psG1[i - 1][j + 1], psG2[i - 1][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                           greenDiff(psG1[ i ][j + 1], psG2[ i ][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                           greenDiff(psG1[i + 1][j + 1], psG2[i + 1][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                                          });
                        // calculate maximum of whole grid by calculating maximum of grid column max values
                        gridMax = std::max({greenDifMax[0], greenDifMax[1], greenDifMax[2]});
                        // adjust index for next column
                        lastIndex ++;
                        lastIndex = lastIndex == gridSize ? 0 : lastIndex;
                    } else if(gridSize == 5) {
                        // compute maximum of differences for fifth column of 5x5 grid and save at position lastIndex
                        greenDifMax[lastIndex] = std::max({greenDiff(psG1[i - 2][j + 2], psG2[i - 2][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                           greenDiff(psG1[i - 1][j + 2], psG2[i - 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                           greenDiff(psG1[ i ][j + 2], psG2[ i ][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                           greenDiff(psG1[i + 1][j + 2], psG2[i + 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                           greenDiff(psG1[i + 2][j + 2], psG2[i + 2][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                                          });
                        // calculate maximum of whole grid by calculating maximum of grid column max values
                        gridMax = std::max({greenDifMax[0], greenDifMax[1], greenDifMax[2], greenDifMax[3], greenDifMax[4]});
                        // adjust index for next column
                        lastIndex ++;
                        lastIndex = lastIndex == gridSize ? 0 : lastIndex;
                    } else if(gridSize == 7) {
                        // compute maximum of differences for 7th column of 7x7 grid and save at position lastIndex
                        greenDifMax[lastIndex] = std::max({greenDiff(psG1[i - 3][j + 3], psG2[i - 3][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                           greenDiff(psG1[i - 2][j + 3], psG2[i - 2][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                           greenDiff(psG1[i - 1][j + 3], psG2[i - 1][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                           greenDiff(psG1[ i ][j + 3], psG2[ i ][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                           greenDiff(psG1[i + 1][j + 3], psG2[i + 1][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                           greenDiff(psG1[i + 2][j + 3], psG2[i + 2][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                           greenDiff(psG1[i + 3][j + 3], psG2[i + 3][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                          });
                        // calculate maximum of whole grid by calculating maximum of grid column max values
                        gridMax = std::max({greenDifMax[0], greenDifMax[1], greenDifMax[2], greenDifMax[3], greenDifMax[4], greenDifMax[5], greenDifMax[6]});
                        // adjust index for next column
                        lastIndex ++;
                        lastIndex = lastIndex == gridSize ? 0 : lastIndex;
                    }

                    if(!adaptive) {
                        // increase motion detection dependent on brightness
                        korr = log2Lut[((int)(psG1[i][j] * scaleGreen)) >> 1];
                    }

                    if (gridMax > thresh - korr) {
#else
                    if (gridMax > thresh) {

#endif

#ifdef PIXELSHIFTDEV
                        sumThr[offset] ++;

                        if(nOf3x3) {
#endif
                            psMask[i][j] = greenWeight;
#ifdef PIXELSHIFTDEV
                        }

                        else if((offset == (frame & 1)) && checkNonGreenVertical) {
                            if(frame > 1) {
                                green[i + offsY][j + offsX] = blueRow ? psG1[i][j] : psG2[i][j];
                            } else {
                                green[i + offsY][j + offsX] = blueRow ? psG2[i][j] : psG1[i][j];;
                            }

                        } else {
                            // at least one of the tested green pixels of the grid is detected as motion
                            paintMotionMask(j + offsX, showMotion, (gridMax - thresh + korr) * blendFactor, showOnlyMask, greenDest, redDest, blueDest);

                            if(skipNext) {
                                // treat the horizontally next pixel also as motion
                                j++;
                                paintMotionMask(j + offsX, showMotion, (gridMax - thresh + korr) * blendFactor, showOnlyMask, greenDest, redDest, blueDest);
                            }
                        }

#endif
                        // do not set the motion pixel values. They have already been set by demosaicer or showMotion
                        continue;
                    }
                }

                if(adaptive) {
                    if(checkNonGreenCross) {
                        // check red cross
                        float redTop    = psRed[i - 1][ j ];
                        float redLeft   = psRed[ i ][j - 1];
                        float redCentre = psRed[ i ][ j ];
                        float redRight  = psRed[ i ][j + 1];
                        float redBottom = psRed[i + 1][ j ];
                        float redDiff = nonGreenDiffCross(redRight, redLeft, redTop, redBottom, redCentre, clippedRed, stddevFactorRed, eperIsoRed, nRead, prnu, showMotion);

                        if(redDiff > 0.f) {
#ifdef PIXELSHIFTDEV

                            if(nOf3x3) {
#endif
                                psMask[i][j] = redBlueWeight;
#ifdef PIXELSHIFTDEV
                            } else {
                                paintMotionMask(j + offsX, showMotion, redDiff, showOnlyMask, redDest, blueDest, greenDest);
                            }

#endif
                            continue;
                        }

                        // check blue cross
                        float blueTop    = psBlue[i - 1][ j ];
                        float blueLeft   = psBlue[ i ][j - 1];
                        float blueCentre = psBlue[ i ][ j ];
                        float blueRight  = psBlue[ i ][j + 1];
                        float blueBottom = psBlue[i + 1][ j ];
                        float blueDiff = nonGreenDiffCross(blueRight, blueLeft, blueTop, blueBottom, blueCentre, clippedBlue, stddevFactorBlue, eperIsoBlue, nRead, prnu, showMotion);

                        if(blueDiff > 0.f) {
#ifdef PIXELSHIFTDEV

                            if(nOf3x3) {
#endif
                                psMask[i][j] = redBlueWeight;
#ifdef PIXELSHIFTDEV
                            } else {
                                paintMotionMask(j + offsX, showMotion, blueDiff, showOnlyMask, blueDest, redDest, greenDest);
                            }

#endif
                            continue;

                        }
                    }

#ifdef PIXELSHIFTDEV

                    if(checkNonGreenHorizontal) {
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

                    if(checkNonGreenVertical) {
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

                    if(checkNonGreenAmaze) {
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

                    if(checkNonGreenCross2) { // for green amaze
                        float greenCentre = (psG1[ i ][ j ] + psG2[ i ][ j ]) / 2.f;
                        float greenAmaze = green[i + offsY][j + offsX];
                        float greenDiffAmaze = nonGreenDiff(greenCentre, greenAmaze, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion);

                        if(greenDiffAmaze > 0.f) {
                            if(nOf3x3) {
                                psMask[i][j] = greenWeight;
                            } else {
                                paintMotionMask(j + offsX, showMotion, greenDiffAmaze, showOnlyMask, greenDest, redDest, blueDest);
                            }

                            continue;
                        }
                    }

                    if(experimental0) { // for experiments

                    }

#endif
                }

                if(showOnlyMask) { // we want only motion mask => paint areas without motion in pure black
                    red[i + offsY][j + offsX] = green[i + offsY][j + offsX] = blue[i + offsY][j + offsX] = 0.f;
                } else if(!(adaptive && nOf3x3)) {
                    // no motion detected, replace the a priori demosaiced values by the pixelshift combined values
                    red[i + offsY][j + offsX] = psRed[i][j];
                    green[i + offsY][j + offsX] = (psG1[i][j] + psG2[i][j]) / 2.f;
                    blue[i + offsY][j + offsX] = psBlue[i][j];
                }
            }
        }

#ifdef PIXELSHIFTDEV

#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            sum[0] += sumThr[0];
            sum[1] += sumThr[1];
        }
#endif
    }


#ifdef PIXELSHIFTDEV
    float percent0 = 100.f * sum[0] / pixelcount;
    float percent1 = 100.f * sum[1] / pixelcount;

    std::cout << fileName <<  " : Green detections at stddev " << std::setprecision( 2 ) << bayerParams.pixelShiftStddevFactorGreen << " : Frame 1/3 : " << std::setprecision( 6 ) << sum[0] << " (" << percent0 << "%)" << " Frame 2/4 : " << sum[1] << " (" << percent1 << "%)" << std::endl;
#endif

    if(adaptive && nOf3x3) {
        if(blurMap) {
            #pragma omp parallel
            {
                gaussianBlur(psMask, psMask, winw, winh, sigma);
            }
        }

        array2D<uint8_t> mask(W, H, ARRAY2D_CLEAR_DATA);

        #pragma omp parallel for schedule(dynamic,16)

        for(int i = winy + border - offsY; i < winh - (border + offsY); ++i) {
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
            array2D<uint8_t> maskInv(W, H);
            invertMask(winx + border - offsX, winw - (border + offsX), winy + border - offsY, winh - (border + offsY), mask, maskInv);
            floodFill4(winx + border - offsX, winw - (border + offsX), winy + border - offsY, winh - (border + offsY), maskInv);
            xorMasks(winx + border - offsX, winw - (border + offsX), winy + border - offsY, winh - (border + offsY), maskInv, mask);
        }


        #pragma omp parallel for schedule(dynamic,16)

        for(int i = winy + border - offsY; i < winh - (border + offsY); ++i) {
#ifdef __SSE2__

            if(smoothTransitions) { //
                vfloat onev = F2V(1.f);
                vfloat smoothv = F2V(smoothFactor);
                int j = winx + border - offsX;

                for(; j < winw - (border + offsX) - 3; j += 4) {
                    vfloat blendv = vmaxf(LVFU(psMask[i][j]), onev) - onev;
                    blendv = pow_F(blendv, smoothv);
                    blendv = vself(vmaskf_eq(smoothv, ZEROV), onev, blendv);
                    STVFU(psMask[i][j], blendv);
                }

                for(; j < winw - (border + offsX); ++j) {
                    psMask[i][j] = smoothFactor == 0.f ? 1.f : pow_F(std::max(psMask[i][j] - 1.f, 0.f), smoothFactor);
                }
            }

#endif
            float *greenDest = green[i + offsY];
            float *redDest = red[i + offsY];
            float *blueDest = blue[i + offsY];

            for(int j = winx + border - offsX; j < winw - (border + offsX); ++j) {
                if(mask[i][j] == 255) {
                    paintMotionMask(j + offsX, showMotion, 0.5f, showOnlyMask, greenDest, redDest, blueDest);
                } else if(showOnlyMask) { // we want only motion mask => paint areas without motion in pure black
                    redDest[j + offsX] = greenDest[j + offsX] = blueDest[j + offsX] = 0.f;
                } else {
                    if(smoothTransitions) {
#ifdef __SSE2__
                        const float blend = psMask[i][j];
#else
                        const float blend = smoothFactor == 0.f ? 1.f : pow_F(std::max(psMask[i][j] - 1.f, 0.f), smoothFactor);
#endif
                        redDest[j + offsX] = intp(blend, redDest[j + offsX], psRed[i][j] );
                        greenDest[j + offsX] = intp(blend, greenDest[j + offsX], (psG1[i][j] + psG2[i][j]) * 0.5f);
                        blueDest[j + offsX] = intp(blend, blueDest[j + offsX], psBlue[i][j]);
                    } else {
                        redDest[j + offsX] = psRed[i][j];
                        greenDest[j + offsX] = (psG1[i][j] + psG2[i][j]) * 0.5f;
                        blueDest[j + offsX] = psBlue[i][j];
                    }
                }
            }
        }
    }

    if(plistener) {
        plistener->setProgress(1.0);
    }
}
