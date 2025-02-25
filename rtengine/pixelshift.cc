////////////////////////////////////////////////////////////////
//
//  Algorithm for Pentax/Sony Pixel Shift raw files with motion detection
//
//  Copyright (C) 2016 - 2018 Ingo Weyrich <heckflosse67@gmx.de>
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
//  along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////

#include <cmath>
#include <stack>

#include "array2D.h"
#include "gauss.h"
#include "median.h"
#include "procparams.h"
#include "rawimagesource.h"
#include "sleef.h"
#include "rtgui/multilangmgr.h"
#include "rtgui/options.h"

//#define BENCHMARK
#include "StopWatch.h"

namespace
{

unsigned fc(const unsigned int cfa[2][2], int r, int c) {
    return cfa[r & 1][c & 1];
}

float greenDiff(float a, float b, float stddevFactor, float eperIso, float nreadIso, float prnu)
{
    // calculate the difference between two green samples
    float gDiff = a - b;
    gDiff *= eperIso;
    gDiff *= gDiff;
    float avg = (a + b) * 0.5f;
    avg *= eperIso;
    prnu *= avg;
    float stddev = stddevFactor * (avg + nreadIso + prnu * prnu);
    return gDiff - stddev;
}

float nonGreenDiffCross(float right, float left, float top, float bottom, float centre, float clippedVal, float stddevFactor, float eperIso, float nreadIso, float prnu)
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
    return std::min(hDiff, vDiff) - stddev;
}

void paintMotionMask(int index, bool showMotion, float *maskDest, float *nonMaskDest0, float *nonMaskDest1)
{
    if(showMotion) {
        // if showMotion is enabled colourize the pixel
        maskDest[index] = 13500.f;
        nonMaskDest1[index] = nonMaskDest0[index] = 0.f;
    }
}

void invertMask(int xStart, int xEnd, int yStart, int yEnd, const array2D<uint8_t> &maskIn, array2D<uint8_t> &maskOut)
{
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif

    for(int i = yStart; i < yEnd; ++i) {
#ifdef _OPENMP
        #pragma omp simd
#endif

        for(int j = xStart; j < xEnd; ++j) {
            maskOut[i][j] = ~maskIn[i][j];
        }
    }
}

void xorMasks(int xStart, int xEnd, int yStart, int yEnd, const array2D<uint8_t> &maskIn, array2D<uint8_t> &maskOut)
{
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif

    for(int i = yStart; i < yEnd; ++i) {
#ifdef _OPENMP
        #pragma omp simd
#endif

        for(int j = xStart; j < xEnd; ++j) {
            maskOut[i][j] ^= maskIn[i][j];
        }
    }
}

void floodFill4Impl(int yin, int xin, int xStart, int xEnd, int yStart, int yEnd, array2D<uint8_t> &mask, std::stack<std::pair<uint16_t, uint16_t>, std::vector<std::pair<uint16_t, uint16_t>>> &coordStack)
{
    coordStack.emplace(xin, yin);

    while(!coordStack.empty()) {
        auto coord = coordStack.top();
        coordStack.pop();
        auto x = coord.first, y = coord.second;

        if(mask[y][x] == 255) {
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
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        std::stack<std::pair<uint16_t, uint16_t>, std::vector<std::pair<uint16_t, uint16_t>>> coordStack;

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,128) nowait
#endif

        for(uint16_t i = yStart; i < yEnd; i++)
        {
            floodFill4Impl(i, xStart, xStart, xEnd, yStart, yEnd, mask, coordStack);
        }

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,128) nowait
#endif

        for(int16_t i = yEnd - 1; i >= 0 ; i--)
        {
            floodFill4Impl(i, xEnd - 1, xStart, xEnd, yStart, yEnd, mask, coordStack);
        }

#ifdef _OPENMP
        #pragma omp sections nowait
#endif
        {
#ifdef _OPENMP
            #pragma omp section
#endif
            {
                uint16_t i = yStart;

                for(uint16_t j = xStart; j < xEnd; j++)
                {
                    floodFill4Impl(i, j, xStart, xEnd, yStart, yEnd, mask, coordStack);
                }
            }
#ifdef _OPENMP
            #pragma omp section
#endif
            {
                uint16_t i = yStart;

                for(uint16_t j = xEnd - 1; j >= xStart; j--)
                {
                    floodFill4Impl(i, j, xStart, xEnd, yStart, yEnd, mask, coordStack);
                }
            }
#ifdef _OPENMP
            #pragma omp section
#endif
            {
                uint16_t i = yEnd;

                for(uint16_t j = xStart; j < xEnd; j++)
                {
                    floodFill4Impl(i, j, xStart, xEnd, yStart, yEnd, mask, coordStack);
                }
            }
#ifdef _OPENMP
            #pragma omp section
#endif
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

void calcFrameBrightnessFactor(unsigned int frame, uint32_t datalen, LUTu *histo[4], float brightnessFactor[4])
{
    float medians[4];

    for (int i = 0; i < 4; ++i) {
        //find median of histogram
        uint32_t lmedian = 0, count = 0;

        while (count < datalen / 2) {
            count += (*histo[i])[lmedian];
            ++lmedian;
        }

        const float weight = (count - datalen / 2.f) / (*histo[i])[lmedian - 1];
        medians[i] = rtengine::intp(weight, (float)(lmedian - 2), (float)(lmedian - 1));
    }

    for (int i = 0; i < 4; ++i) {
        brightnessFactor[i] = medians[frame] / medians[i];
    }

}

}

using namespace std;
using namespace rtengine;
void RawImageSource::pixelshift(int winx, int winy, int winw, int winh, const procparams::RAWParams &rawParamsIn, unsigned int frame, const std::string &make, const std::string &model, float rawWpCorrection)
{
BENCHFUN
    if(numFrames != 4) { // fallback for non pixelshift files
        amaze_demosaic_RT(winx, winy, winw, winh, rawData, red, green, blue, options.chunkSizeAMAZE, options.measure);
        return;
    }

    procparams::RAWParams::BayerSensor bayerParams = rawParamsIn.bayersensor;

    bool motionDetection = true;

    if(bayerParams.pixelShiftMotionCorrectionMethod == procparams::RAWParams::BayerSensor::PSMotionCorrectionMethod::AUTO) {
        bool pixelShiftEqualBright = bayerParams.pixelShiftEqualBright;
        bayerParams.setPixelShiftDefaults();
        bayerParams.pixelShiftEqualBright = pixelShiftEqualBright;
    } else if(bayerParams.pixelShiftMotionCorrectionMethod == procparams::RAWParams::BayerSensor::PSMotionCorrectionMethod::OFF) {
        motionDetection = false;
        bayerParams.pixelShiftShowMotion = false;
    }

    const unsigned int cfarray[2][2] = {{FC(0,0), FC(0,1)}, {FC(1,0), FC(1,1)}};
    const bool showMotion = bayerParams.pixelShiftShowMotion;
    const bool showOnlyMask = bayerParams.pixelShiftShowMotionMaskOnly && showMotion;
    const float smoothFactor = 1.0 - bayerParams.pixelShiftSmoothFactor;

    if(motionDetection) {
        if(!showOnlyMask) {
            if(bayerParams.pixelShiftMedian || bayerParams.pixelShiftAverage) { // We need the demosaiced frames for motion correction
                if (bayerParams.pixelShiftDemosaicMethod == bayerParams.getPSDemosaicMethodString(procparams::RAWParams::BayerSensor::PSDemosaicMethod::LMMSE)) {
                    lmmse_interpolate_omp(winw, winh, *(rawDataFrames[0]), red, green, blue, bayerParams.lmmse_iterations);
                } else if (bayerParams.pixelShiftDemosaicMethod == bayerParams.getPSDemosaicMethodString(procparams::RAWParams::BayerSensor::PSDemosaicMethod::AMAZEVNG4)) {
                    dual_demosaic_RT (true, rawParamsIn, winw, winh, *(rawDataFrames[0]), red, green, blue, bayerParams.dualDemosaicContrast, true);
                } else if (bayerParams.pixelShiftDemosaicMethod == bayerParams.getPSDemosaicMethodString(procparams::RAWParams::BayerSensor::PSDemosaicMethod::RCDVNG4)) {
                    dual_demosaic_RT (true, rawParamsIn, winw, winh, *(rawDataFrames[0]), red, green, blue, bayerParams.dualDemosaicContrast, true);
                } else {
                    amaze_demosaic_RT(winx, winy, winw, winh, *(rawDataFrames[0]), red, green, blue, options.chunkSizeAMAZE, options.measure);
                }
                multi_array2D<float, 3> redTmp(winw, winh);
                multi_array2D<float, 3> greenTmp(winw, winh);
                multi_array2D<float, 3> blueTmp(winw, winh);

                for(int i = 0; i < 3; i++) {
                    if (bayerParams.pixelShiftDemosaicMethod == bayerParams.getPSDemosaicMethodString(procparams::RAWParams::BayerSensor::PSDemosaicMethod::LMMSE)) {
                        lmmse_interpolate_omp(winw, winh, *(rawDataFrames[i + 1]), redTmp[i], greenTmp[i], blueTmp[i], bayerParams.lmmse_iterations);
                    } else if (bayerParams.pixelShiftDemosaicMethod == bayerParams.getPSDemosaicMethodString(procparams::RAWParams::BayerSensor::PSDemosaicMethod::AMAZEVNG4)) {
                        dual_demosaic_RT (true, rawParamsIn, winw, winh, *(rawDataFrames[i + 1]), redTmp[i], greenTmp[i], blueTmp[i], bayerParams.dualDemosaicContrast, true);
                    } else if (bayerParams.pixelShiftDemosaicMethod == bayerParams.getPSDemosaicMethodString(procparams::RAWParams::BayerSensor::PSDemosaicMethod::RCDVNG4)) {
                        dual_demosaic_RT (true, rawParamsIn, winw, winh, *(rawDataFrames[i + 1]), redTmp[i], greenTmp[i], blueTmp[i], bayerParams.dualDemosaicContrast, true);
                    } else {
                        amaze_demosaic_RT(winx, winy, winw, winh, *(rawDataFrames[i + 1]), redTmp[i], greenTmp[i], blueTmp[i], options.chunkSizeAMAZE, options.measure);
                    }
                }
                if(bayerParams.pixelShiftMedian) {

#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16)
#endif

                    for(int i = winy + border; i < winh - border; i++) {
                        for(int j = winx + border; j < winw - border; j++) {
                            red[i][j] = median(red[i][j], redTmp[0][i + 1][j], redTmp[1][i + 1][j + 1], redTmp[2][i][j + 1]);
                        }

                        for(int j = winx + border; j < winw - border; j++) {
                            green[i][j] = median(green[i][j], greenTmp[0][i + 1][j], greenTmp[1][i + 1][j + 1], greenTmp[2][i][j + 1]);
                        }

                        for(int j = winx + border; j < winw - border; j++) {
                            blue[i][j] = median(blue[i][j], blueTmp[0][i + 1][j], blueTmp[1][i + 1][j + 1], blueTmp[2][i][j + 1]);
                        }
                    }
                } else {
#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16)
#endif

                    for(int i = winy + border; i < winh - border; i++) {
                        for(int j = winx + border; j < winw - border; j++) {
                            red[i][j] = 0.25f * ((red[i][j] + redTmp[0][i + 1][j]) + (redTmp[1][i + 1][j + 1] + redTmp[2][i][j + 1]));
                        }

                        for(int j = winx + border; j < winw - border; j++) {
                            green[i][j] = 0.25f * ((green[i][j] + greenTmp[0][i + 1][j]) + (greenTmp[1][i + 1][j + 1] + greenTmp[2][i][j + 1]));
                        }

                        for(int j = winx + border; j < winw - border; j++) {
                            blue[i][j] = 0.25f * ((blue[i][j] + blueTmp[0][i + 1][j]) + (blueTmp[1][i + 1][j + 1] + blueTmp[2][i][j + 1]));
                        }
                    }
                }
            } else {
                if (bayerParams.pixelShiftDemosaicMethod == bayerParams.getPSDemosaicMethodString(procparams::RAWParams::BayerSensor::PSDemosaicMethod::LMMSE)) {
                    lmmse_interpolate_omp(winw, winh, rawData, red, green, blue, bayerParams.lmmse_iterations);
                } else if (bayerParams.pixelShiftDemosaicMethod == bayerParams.getPSDemosaicMethodString(procparams::RAWParams::BayerSensor::PSDemosaicMethod::AMAZEVNG4)) {
                    procparams::RAWParams rawParamsTmp = rawParamsIn;
                    rawParamsTmp.bayersensor.method = procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::AMAZEVNG4);
                    dual_demosaic_RT (true, rawParamsTmp, winw, winh, rawData, red, green, blue, bayerParams.dualDemosaicContrast, true);
                } else if (bayerParams.pixelShiftDemosaicMethod == bayerParams.getPSDemosaicMethodString(procparams::RAWParams::BayerSensor::PSDemosaicMethod::RCDVNG4)) {
                    procparams::RAWParams rawParamsTmp = rawParamsIn;
                    rawParamsTmp.bayersensor.method = procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::RCDVNG4);
                    dual_demosaic_RT (true, rawParamsTmp, winw, winh, rawData, red, green, blue, bayerParams.dualDemosaicContrast, true);
                } else {
                    amaze_demosaic_RT(winx, winy, winw, winh, rawData, red, green, blue, options.chunkSizeAMAZE, options.measure);
                }
            }
        }
    }

    constexpr float stddevFactorGreen = 25.f;
    constexpr float stddevFactorRed = 25.f;
    constexpr float stddevFactorBlue = 25.f;
    constexpr float prnu = 0.01f;
    constexpr float redBlueWeight = 0.7f + 1.f;
    float eperIso = bayerParams.pixelShiftEperIso;
    const bool checkNonGreenCross = bayerParams.pixelShiftNonGreenCross;
    const bool checkGreen = bayerParams.pixelShiftGreen;
    constexpr float greenWeight = 2.f;
    const bool blurMap = bayerParams.pixelShiftBlur;
    const float sigma = bayerParams.pixelShiftSigma;
    constexpr float noMotion = 0.99f;
    constexpr float threshold = 3.f + 9 * noMotion;
    const bool holeFill = bayerParams.pixelShiftHoleFill;
    const bool equalBrightness = bayerParams.pixelShiftEqualBright;
    const bool equalChannel = bayerParams.pixelShiftEqualBrightChannel;
    const bool smoothTransitions = blurMap && bayerParams.pixelShiftSmoothFactor > 0.;

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
                                       1.5f   // ISO 51200
                                     };

    static const float ePerIsoK3II = 0.35f;

    // currently we use the same values for K-1 and K-1 Mark II, though for the K-1 Mark II the values seem a bit aggressive
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
                                       2.4f,  // ISO 204800
                                       2.4f,  // ISO 256000  // these are for K-1 Mark II to avoid crashes when using high-ISO files
                                       2.4f,  // ISO 320000
                                       2.4f,  // ISO 409600
                                       2.4f,  // ISO 512000
                                       2.4f,  // ISO 640000
                                       2.4f   // ISO 819200
                                   };

    static const float ePerIsoK1 = 0.75f;

    // currently nReadK70 is used for K-70 and KP
    static const float nReadK70[] = {  4.0f,  // ISO 100
                                       4.0f,  // ISO 125
                                       4.0f,  // ISO 160
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
                                       3.0f,  // ISO 64000
                                       3.0f,  // ISO 80000
                                       3.0f,  // ISO 102400
                                       3.0f,  // ISO 128000
                                       3.0f,  // ISO 160000
                                       3.0f,  // ISO 204800
                                       3.0f,  // ISO 256000
                                       3.0f,  // ISO 320000
                                       3.0f,  // ISO 409600
                                       3.0f,  // ISO 512000
                                       3.0f,  // ISO 640000
                                       3.0f   // ISO 819200
                                    };

    static const float ePerIsoK70 = 0.5f;

    // preliminary ILCE-7RM3 data, good fidelity except from A) small inaccuracy at places
    // due to integer scaling quantization, B) much different noise behaviour of PDAF pixels
    static const float nReadILCE7RM3[] = { 4.2f,  // ISO 100
                                           3.9f,  // ISO 125
                                           3.6f,  // ISO 160
                                           3.55f, // ISO 200
                                           3.5f,  // ISO 250
                                           3.45f, // ISO 320
                                           3.35f, // ISO 400
                                           3.3f,  // ISO 500
                                           1.3f,  // ISO 640
                                           1.2f,  // ISO 800
                                           1.2f,  // ISO 1000
                                           1.2f,  // ISO 1250
                                           1.15f, // ISO 1600
                                           1.2f,  // ISO 2000
                                           1.15f, // ISO 2500
                                           1.15f, // ISO 3200
                                           1.1f,  // ISO 4000
                                           1.1f,  // ISO 5000
                                           1.05f, // ISO 6400
                                           1.05f, // ISO 8000
                                           1.05f, // ISO 10000
                                           1.0f,  // ISO 12800
                                           1.0f,  // ISO 16000
                                           1.0f,  // ISO 20000
                                           1.0f,  // ISO 25600
                                           1.0f,  // ISO 32000
                                           1.0f,  // ISO 40000
                                           1.0f,  // ISO 51200
                                           1.1f,  // ISO 64000
                                           1.1f,  // ISO 80000
                                           1.1f,  // ISO 102400
                                        };

    static const float ePerIsoILCE7RM3 = 0.8f;

    //TODO: Add data for ILCE-7RM4, and ILCE-1
    if(plistener) {
        plistener->setProgressStr(Glib::ustring::compose(M("TP_RAW_DMETHOD_PROGRESSBAR"), M("TP_RAW_PIXELSHIFT")));
        plistener->setProgress(0.0);
    }

    if(motionDetection && blurMap && smoothFactor == 0.f && !showMotion) {
        if(plistener) {
            plistener->setProgress(1.0);
        }

        return;
    }

    float nRead;
    float eperIsoModel;

    int nReadIndex = static_cast<int>(round(log2(idata->getISOSpeed() /  100.f) * 3.f));

    if(make.find("Sony") != string::npos) {
        if(model.find("ILCE-7RM3") != string::npos) {
            nRead = nReadILCE7RM3[nReadIndex];
            eperIsoModel = ePerIsoILCE7RM3;
        /* TODO: When we have data for missing ILCE-7RM4, and ILCE-1, add it here
        } else if(model.find("ILCE-7RM4") != string::npos) {
        } else if(model.find("ILCE-1") != string::npos) {
        */
        } else { // default to ILCE-7RM3 for Sony cameras without data
           nRead = nReadILCE7RM3[nReadIndex];
           eperIsoModel = ePerIsoILCE7RM3;
        }
    } else { // Pentax
        if(model.find("K-3") != string::npos) {
            nRead = nReadK3II[nReadIndex];
            eperIsoModel = ePerIsoK3II;
        } else if(model.find("K-1") != string::npos) { // this also matches K-1 Mark II
            nRead = nReadK1[nReadIndex];
            eperIsoModel = ePerIsoK1;
        } else { // as long as we don't have values for Pentax KP, we use the values from K-70
            nRead = nReadK70[nReadIndex];
            eperIsoModel = ePerIsoK70;
        }
    }

    eperIsoModel *= pow(2.f, eperIso - 1.f);

    eperIso = eperIsoModel * (100.f / (rawWpCorrection * idata->getISOSpeed()));

    const float eperIsoRed = (eperIso / scale_mul[0]) * (65535.f / (c_white[0] - c_black[0]));
    const float eperIsoGreen = (eperIso / scale_mul[1]) * (65535.f / (c_white[1] - c_black[1]));
    const float eperIsoBlue = (eperIso / scale_mul[2]) * (65535.f / (c_white[2] - c_black[2]));

    const float clippedRed = 65535.f / scale_mul[0];
    const float clippedBlue = 65535.f / scale_mul[2];

    nRead *= nRead;

    // calculate channel median brightness for each frame
    float greenBrightness[4] = {1.f, 1.f, 1.f, 1.f};
    float redBrightness[4] = {1.f, 1.f, 1.f, 1.f};
    float blueBrightness[4] = {1.f, 1.f, 1.f, 1.f};

    if(equalBrightness) {
        if(rawDirty) {
            LUTu *histogreen[4];
            LUTu *histored[4];
            LUTu *histoblue[4];

            for(int i = 0; i < 4; ++i) {
                histogreen[i] = new LUTu(65536, LUT_CLIP_BELOW | LUT_CLIP_ABOVE, true);
                histored[i] = new LUTu(65536, LUT_CLIP_BELOW | LUT_CLIP_ABOVE, true);
                histoblue[i] = new LUTu(65536, LUT_CLIP_BELOW | LUT_CLIP_ABOVE, true);
            }

#ifdef _OPENMP
            #pragma omp parallel
#endif
            {
                LUTu *histogreenThr[4];
                LUTu *historedThr[4];
                LUTu *histoblueThr[4];

                for(int i = 0; i < 4; ++i) {
                    histogreenThr[i] = new LUTu(65536, LUT_CLIP_BELOW | LUT_CLIP_ABOVE, true);
                    historedThr[i] = new LUTu(65536, LUT_CLIP_BELOW | LUT_CLIP_ABOVE, true);
                    histoblueThr[i] = new LUTu(65536, LUT_CLIP_BELOW | LUT_CLIP_ABOVE, true);
                }

#ifdef _OPENMP
                #pragma omp for schedule(dynamic,16) nowait
#endif

                for(int i = winy + 1; i < winh - 1; ++i) {
                    int j = winx + 1;
                    int c = fc(cfarray, i, j);

                    bool bluerow = (c + fc(cfarray, i, j + 1)) == 3;

                    for(int j = winx + 1, offset = fc(cfarray, i, j) & 1; j < winw - 1; ++j, offset ^= 1) {
                        (*histogreenThr[1 - offset])[(*rawDataFrames[1 - offset])[i - offset + 1][j]]++;
                        (*histogreenThr[3 - offset])[(*rawDataFrames[3 - offset])[i + offset][j + 1]]++;

                        if(bluerow) {
                            (*historedThr[2 - offset])[(*rawDataFrames[2 - offset])[i + 1][j - offset + 1]]++;
                            (*histoblueThr[(offset << 1) + offset])[(*rawDataFrames[(offset << 1) + offset])[i][j + offset]]++;
                        } else {
                            (*historedThr[(offset << 1) + offset])[(*rawDataFrames[(offset << 1) + offset])[i][j + offset]]++;
                            (*histoblueThr[2 - offset])[(*rawDataFrames[2 - offset])[i + 1][j - offset + 1]]++;
                        }
                    }
                }

#ifdef _OPENMP
                #pragma omp critical
#endif
                {
                    for(int i = 0; i < 4; ++i) {
                        (*histogreen[i]) += (*histogreenThr[i]);
                        delete histogreenThr[i];
                        (*histored[i]) += (*historedThr[i]);
                        delete historedThr[i];
                        (*histoblue[i]) += (*histoblueThr[i]);
                        delete histoblueThr[i];
                    }
                }
            }

            calcFrameBrightnessFactor(frame, (winh - 2) * (winw - 2) / 4, histored, redBrightness);
            calcFrameBrightnessFactor(frame, (winh - 2) * (winw - 2) / 4, histoblue, blueBrightness);
            calcFrameBrightnessFactor(frame, (winh - 2) * (winw - 2) / 2, histogreen, greenBrightness);

            for(int i = 0; i < 4; ++i) {
                psRedBrightness[i] = redBrightness[i];
                psGreenBrightness[i] = greenBrightness[i];
                psBlueBrightness[i] = blueBrightness[i];
            }
            rawDirty = false;

            for(int i = 0; i < 4; ++i) {
                delete histored[i];
                delete histoblue[i];
                delete histogreen[i];
            }
            if(plistener) {
                plistener->setProgress(0.15);
            }

        } else {
            for(int i = 0; i < 4; ++i) {
                redBrightness[i] = psRedBrightness[i];
                greenBrightness[i] = psGreenBrightness[i];
                blueBrightness[i] = psBlueBrightness[i];
            }
        }
        if(!equalChannel) {
            for(int i = 0; i < 4; ++i) {
                redBrightness[i] = blueBrightness[i] = greenBrightness[i];
            }
        }
    }


    if(motionDetection) {
        // fill channels psRed and psBlue
        array2D<float> psRed(winw + 32, winh); // increase width to avoid cache conflicts
        array2D<float> psBlue(winw + 32, winh);

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

        for(int i = winy + 1; i < winh - 1; ++i) {
            float *nonGreenDest0 = psRed[i];
            float *nonGreenDest1 = psBlue[i];
            float ngbright[2][4] = {{redBrightness[0], redBrightness[1], redBrightness[2], redBrightness[3]},
                                    {blueBrightness[0], blueBrightness[1], blueBrightness[2], blueBrightness[3]}
                                   };
            int ng = 0;
            int j = winx + 1;
            int c = fc(cfarray, i, j);

            if((c + fc(cfarray, i, j + 1)) == 3) {
                // row with blue pixels => swap destination pointers for non green pixels
                std::swap(nonGreenDest0, nonGreenDest1);
                ng ^= 1;
            }

            // offset to keep the code short. It changes its value between 0 and 1 for each iteration of the loop
            unsigned int offset = c & 1;

            for(; j < winw - 1; ++j) {
                // store the non green values from the 4 frames into 2 temporary planes
                nonGreenDest0[j] = (*rawDataFrames[(offset << 1) + offset])[i][j + offset] * ngbright[ng][(offset << 1) + offset];
                nonGreenDest1[j] = (*rawDataFrames[2 - offset])[i + 1][j - offset + 1] * ngbright[ng ^ 1][2 - offset];
                offset ^= 1; // 0 => 1 or 1 => 0
            }
        }

        if(plistener) {
            plistener->setProgress(0.3);
        }

        // now that the temporary planes are filled for easy access we do the motion detection
        array2D<float> psMask(winw, winh);

        int offsX = 0, offsY = 0;

        if(!bayerParams.pixelShiftMedian) {
            // We have to adjust the offsets for the selected subframe we use for areas with motion
            switch(frame) {
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
            // offset to keep the code short. It changes its value between 0 and 1 for each iteration of the loop
            unsigned int offset = fc(cfarray, i, winx + border - offsX) & 1;

            for(int j = winx + border - offsX; j < winw - (border + offsX); ++j, offset ^= 1) {
                psMask[i][j] = noMotion;

                if(checkGreen) {
                    if(greenDiff((*rawDataFrames[1 - offset])[i - offset + 1][j] * greenBrightness[1 - offset], (*rawDataFrames[3 - offset])[i + offset][j + 1] * greenBrightness[3 - offset], stddevFactorGreen, eperIsoGreen, nRead, prnu) > 0.f) {
                        psMask[i][j] = greenWeight;
                        // do not set the motion pixel values. They have already been set by demosaicer
                        continue;
                    }
                }

                if(checkNonGreenCross) {
                    // check red cross
                    float redTop    = psRed[i - 1][j];
                    float redLeft   = psRed[i][j - 1];
                    float redCentre = psRed[i][j];
                    float redRight  = psRed[i][j + 1];
                    float redBottom = psRed[i + 1][j];
                    float redDiff   = nonGreenDiffCross(redRight, redLeft, redTop, redBottom, redCentre, clippedRed, stddevFactorRed, eperIsoRed, nRead, prnu);

                    if(redDiff > 0.f) {
                        psMask[i][j] = redBlueWeight;
                        continue;
                    }

                    // check blue cross
                    float blueTop    = psBlue[i - 1][j];
                    float blueLeft   = psBlue[i][j - 1];
                    float blueCentre = psBlue[i][j];
                    float blueRight  = psBlue[i][j + 1];
                    float blueBottom = psBlue[i + 1][j];
                    float blueDiff   = nonGreenDiffCross(blueRight, blueLeft, blueTop, blueBottom, blueCentre, clippedBlue, stddevFactorBlue, eperIsoBlue, nRead, prnu);

                    if(blueDiff > 0.f) {
                        psMask[i][j] = redBlueWeight;
                        continue;
                    }
                }
            }
        }

        if(plistener) {
            plistener->setProgress(0.45);
        }

        if(blurMap) {
#ifdef _OPENMP
            #pragma omp parallel
#endif
            {
                gaussianBlur(psMask, psMask, winw, winh, sigma);
            }
            if(plistener) {
                plistener->setProgress(0.6);
            }
        }

        array2D<uint8_t> mask(winw, winh, ARRAY2D_CLEAR_DATA);

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

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

        if(plistener) {
            plistener->setProgress(0.75);
        }

        if(holeFill) {
            array2D<uint8_t> maskInv(winw, winh);
            invertMask(winx + border - offsX, winw - (border + offsX), winy + border - offsY, winh - (border + offsY), mask, maskInv);
            floodFill4(winx + border - offsX, winw - (border + offsX), winy + border - offsY, winh - (border + offsY), maskInv);
            xorMasks(winx + border - offsX, winw - (border + offsX), winy + border - offsY, winh - (border + offsY), maskInv, mask);
        }

        if(plistener) {
            plistener->setProgress(0.9);
        }

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

        for(int i = winy + border - offsY; i < winh - (border + offsY); ++i) {
#ifdef __SSE2__

            // pow() is expensive => pre calculate blend factor using SSE
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

            // offset to keep the code short. It changes its value between 0 and 1 for each iteration of the loop
            unsigned int offset = fc(cfarray, i, winx + border - offsX) & 1;

            for(int j = winx + border - offsX; j < winw - (border + offsX); ++j, offset ^= 1) {
                if(showOnlyMask) {
                    if(smoothTransitions) { // we want only motion mask => paint areas according to their motion (dark = no motion, bright = motion)
#ifdef __SSE2__
                        // use pre calculated blend factor
                        const float blend = psMask[i][j];
#else
                        const float blend = smoothFactor == 0.f ? 1.f : pow_F(std::max(psMask[i][j] - 1.f, 0.f), smoothFactor);
#endif
                        redDest[j + offsX] = greenDest[j + offsX] = blueDest[j + offsX] = blend * 32768.f;
                    } else {
                        redDest[j + offsX] = greenDest[j + offsX] = blueDest[j + offsX] = mask[i][j] == 255 ? 65535.f : 0.f;
                    }
                } else if(mask[i][j] == 255) {
                    paintMotionMask(j + offsX, showMotion, greenDest, redDest, blueDest);
                } else {
                    if(smoothTransitions) {
#ifdef __SSE2__
                        // use pre calculated blend factor
                        const float blend = psMask[i][j];
#else
                        const float blend = smoothFactor == 0.f ? 1.f : pow_F(std::max(psMask[i][j] - 1.f, 0.f), smoothFactor);
#endif
                        redDest[j + offsX] = intp(blend, showMotion ? 0.f : redDest[j + offsX], psRed[i][j] );
                        greenDest[j + offsX] = intp(blend, showMotion ? 13500.f : greenDest[j + offsX], ((*rawDataFrames[1 - offset])[i - offset + 1][j] * greenBrightness[1 - offset] + (*rawDataFrames[3 - offset])[i + offset][j + 1] * greenBrightness[3 - offset]) * 0.5f);
                        blueDest[j + offsX] = intp(blend, showMotion ? 0.f : blueDest[j + offsX], psBlue[i][j]);
                    } else {
                        redDest[j + offsX] = psRed[i][j];
                        greenDest[j + offsX] = ((*rawDataFrames[1 - offset])[i - offset + 1][j] * greenBrightness[1 - offset] + (*rawDataFrames[3 - offset])[i + offset][j + 1] * greenBrightness[3 - offset]) * 0.5f;
                        blueDest[j + offsX] = psBlue[i][j];
                    }
                }
            }
        }
    } else {
        // motion detection off => combine the 4 raw frames
        float ngbright[2][4] = {{redBrightness[0], redBrightness[1], redBrightness[2], redBrightness[3]},
                                {blueBrightness[0], blueBrightness[1], blueBrightness[2], blueBrightness[3]}
        };
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

        for(int i = winy + 1; i < winh - 1; ++i) {
            float *nonGreenDest0 = red[i];
            float *nonGreenDest1 = blue[i];
            int ng = 0;
            int j = winx + 1;
            int c = fc(cfarray, i, j);

            if((c + fc(cfarray, i, j + 1)) == 3) {
                // row with blue pixels => swap destination pointers for non green pixels
                std::swap(nonGreenDest0, nonGreenDest1);
                ng ^= 1;
            }

            // offset to keep the code short. It changes its value between 0 and 1 for each iteration of the loop
            unsigned int offset = c & 1;

            for(; j < winw - 1; ++j) {
                // set red, green and blue values
                green[i][j] = ((*rawDataFrames[1 - offset])[i - offset + 1][j] * greenBrightness[1 - offset] + (*rawDataFrames[3 - offset])[i + offset][j + 1] * greenBrightness[3 - offset]) * 0.5f;
                nonGreenDest0[j] = (*rawDataFrames[(offset << 1) + offset])[i][j + offset] * ngbright[ng][(offset << 1) + offset];
                nonGreenDest1[j] = (*rawDataFrames[2 - offset])[i + 1][j - offset + 1] * ngbright[ng ^ 1][2 - offset];
                offset ^= 1; // 0 => 1 or 1 => 0
            }
        }
    }

    if(plistener) {
        plistener->setProgress(1.0);
    }
}
