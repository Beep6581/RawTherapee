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
#define BENCHMARK
#include "StopWatch.h"

namespace
{

float greenDiff(float a, float b, bool adaptive, float stddevFactor, float eperIso, float nreadIso, float prnu, bool showMotion)
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

}

using namespace std;
using namespace rtengine;

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

    if(model.find("K-3") != string::npos) {
        nRead = nReadK3II[static_cast<int>(round(log2(idata->getISOSpeed() /  100.f) * 3.f))];
        eperIsoModel = ePerIsoK3II;
    } else if(model.find("K-1") != string::npos) {
        nRead = nReadK1[static_cast<int>(round(log2(idata->getISOSpeed() /  100.f) * 3.f))];
        eperIsoModel = ePerIsoK1;
    } else {
        nRead = nReadK70[static_cast<int>(round(log2(idata->getISOSpeed() /  100.f) * 3.f))];
        eperIsoModel = ePerIsoK70;
    }

    nRead *= pow(2.f, nreadIso);
    eperIsoModel *= pow(2.f, eperIso);
    eperIso = eperIsoModel * (100.f / (rawWpCorrection * idata->getISOSpeed()));
    float eperIsoGreen = eperIso * scaleGreen;

    printf("Pixelshift parameters : gridSize %d\tadaptive %d\tstdDevFactorGreen %f\telectrons %1.8f\tnread %f\tprnu %1.1f%%\n", gridSize, adaptive, stddevFactorGreen, eperIso, nRead, prnu);

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

        if (c == 2 || ((c & 1) && FC(i, j + 1) == 2)) {
            // row with blue pixels => swap destination pointers for non green pixels
            std::swap(nonGreenDest0, nonGreenDest1);
            std::swap(scaleNonGreen0, scaleNonGreen2);
            std::swap(eperIsoNonGreen0, eperIsoNonGreen2);
            std::swap(stddevFactorNonGreen0, stddevFactorNonGreen2);
        }

        // offset to keep the code short. It changes its value between 0 and 1 for each iteration of the loop
        unsigned int offset = (c & 1);

        float greenDifMax[gridSize];

        // motion detection checks the grid around the pixel for differences in green channels
        if(detectMotion || adaptive) {
            if(gridSize == 3) {
                // compute maximum of differences for first two columns of 3x3 grid
                greenDifMax[0] =  max(greenDiff((*rawDataFrames[1 - offset])[i - offset][j - 1], (*rawDataFrames[3 - offset])[i + offset - 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[0 + offset])[i + offset][j - 1], (*rawDataFrames[2 + offset])[i - offset + 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[1 - offset])[i - offset + 2][j - 1], (*rawDataFrames[3 - offset])[i + offset + 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                     );
                greenDifMax[1] =  max(greenDiff((*rawDataFrames[0 + offset])[i + offset - 1][j], (*rawDataFrames[2 + offset])[i - offset][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[1 - offset])[i - offset + 1][j], (*rawDataFrames[3 - offset])[i + offset][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[0 + offset])[i + offset + 1][j], (*rawDataFrames[2 + offset])[i - offset + 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                     );
            } else if(gridSize == 5) {
                // compute maximum of differences for first four columns of 5x5 grid
                greenDifMax[0] =  max(greenDiff((*rawDataFrames[1 - offset])[i - offset - 1][j - 2], (*rawDataFrames[3 - offset])[i + offset - 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[0 + offset])[i + offset - 1][j - 2], (*rawDataFrames[2 + offset])[i - offset][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[1 - offset])[i - offset + 1][j - 2], (*rawDataFrames[3 - offset])[i + offset][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[0 + offset])[i + offset + 1][j - 2], (*rawDataFrames[2 + offset])[i - offset + 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[1 - offset])[i - offset + 3][j - 2], (*rawDataFrames[3 - offset])[i + offset + 2][j - 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                     );
                greenDifMax[1] =  max(greenDiff((*rawDataFrames[0 + offset])[i + offset - 2][j - 1], (*rawDataFrames[2 + offset])[i - offset - 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[1 - offset])[i - offset][j - 1], (*rawDataFrames[3 - offset])[i + offset - 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[0 + offset])[i + offset][j - 1], (*rawDataFrames[2 + offset])[i - offset + 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[1 - offset])[i - offset + 2][j - 1], (*rawDataFrames[3 - offset])[i + offset + 1][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[0 + offset])[i + offset + 2][j - 1], (*rawDataFrames[2 + offset])[i - offset + 3][j], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                     );
                greenDifMax[2] =  max(greenDiff((*rawDataFrames[1 - offset])[i - offset - 1][j], (*rawDataFrames[3 - offset])[i + offset - 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[0 + offset])[i + offset - 1][j], (*rawDataFrames[2 + offset])[i - offset][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[1 - offset])[i - offset + 1][j], (*rawDataFrames[3 - offset])[i + offset][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[0 + offset])[i + offset + 1][j], (*rawDataFrames[2 + offset])[i - offset + 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[1 - offset])[i - offset + 3][j], (*rawDataFrames[3 - offset])[i + offset + 2][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                     );
                greenDifMax[3] =  max(greenDiff((*rawDataFrames[0 + offset])[i + offset - 2][j + 1], (*rawDataFrames[2 + offset])[i - offset - 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[1 - offset])[i - offset][j + 1], (*rawDataFrames[3 - offset])[i + offset - 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[0 + offset])[i + offset][j + 1], (*rawDataFrames[2 + offset])[i - offset + 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[1 - offset])[i - offset + 2][j + 1], (*rawDataFrames[3 - offset])[i + offset + 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                      greenDiff((*rawDataFrames[0 + offset])[i + offset + 2][j + 1], (*rawDataFrames[2 + offset])[i - offset + 3][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
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

                if(gridSize < 2) {
                    // compute difference for current pixel and skip next pixel, that's the method from dcrawps
                    gridMax = greenDiff((*rawDataFrames[1 - offset])[i - offset + 1][j], (*rawDataFrames[3 - offset])[i + offset][j + 1], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion);
                    skipNext = skip;
                } else if(gridSize == 3) {
                    // compute maximum of differences for third column of 3x3 grid and save at position lastIndex
                    greenDifMax[lastIndex] = max(greenDiff((*rawDataFrames[1 - offset])[i - offset][j + 1], (*rawDataFrames[3 - offset])[i + offset - 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                 greenDiff((*rawDataFrames[0 + offset])[i + offset][j + 1], (*rawDataFrames[2 + offset])[i - offset + 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                 greenDiff((*rawDataFrames[1 - offset])[i - offset + 2][j + 1], (*rawDataFrames[3 - offset])[i + offset + 1][j + 2], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                                );
                    gridMax = max(greenDifMax[0], greenDifMax[1], greenDifMax[2]);
                } else if(gridSize == 5) {
                    // compute maximum of differences for fifth column of 5x5 grid and save at position lastIndex
                    greenDifMax[lastIndex] = max(greenDiff((*rawDataFrames[1 - offset])[i - offset - 1][j + 2], (*rawDataFrames[3 - offset])[i + offset - 2][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                 greenDiff((*rawDataFrames[0 + offset])[i + offset - 1][j + 2], (*rawDataFrames[2 + offset])[i - offset][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                 greenDiff((*rawDataFrames[1 - offset])[i - offset + 1][j + 2], (*rawDataFrames[3 - offset])[i + offset][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                 greenDiff((*rawDataFrames[0 + offset])[i + offset + 1][j + 2], (*rawDataFrames[2 + offset])[i - offset + 2][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion),
                                                 greenDiff((*rawDataFrames[1 - offset])[i - offset + 3][j + 2], (*rawDataFrames[3 - offset])[i + offset + 2][j + 3], adaptive, stddevFactorGreen, eperIsoGreen, nRead, prnu, showMotion)
                                                );
                    gridMax = max(greenDifMax[0], greenDifMax[1], greenDifMax[2], greenDifMax[3], greenDifMax[4]);
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
                        offset ^= 1;
                    }

                    // do not set the motion pixel values. They have already been set by demosaicer or showMotion
                    continue;
                }

                if(adaptive && checkNonGreenCross) {
                    float ngCentre = (*rawDataFrames[(offset << 1) + offset])[i][j + offset];
                    float ngRight = (*rawDataFrames[((offset ^ 1) << 1) + (offset ^ 1)])[i][j + (offset ^ 1) + 1];
                    float ngLeft = (*rawDataFrames[((offset ^ 1) << 1) + (offset ^ 1)])[i][j + (offset ^ 1) - 1];
                    float ngTop = (*rawDataFrames[((offset << 1) + offset) ^ 1])[i][j + offset];
                    float ngBottom = (*rawDataFrames[((offset << 1) + offset) ^ 1])[i + 2][j + offset];
                    float gridMax = nonGreenDiffCross(ngRight, ngLeft, ngTop, ngBottom, ngCentre, stddevFactorNonGreen0, eperIsoNonGreen0, nRead, prnu, showMotion);

                    if(gridMax > 0.f) {
                        if(showMotion) {
                            float blend = gridMax * blendFactor;

                            if(!showOnlyMask) {
                                // if showMotion is enabled colourize the pixel
                                nonGreenDest0[j + offsX] = 1000.f + 25000.f * blend;
                                nonGreenDest1[j + offsX] = greenDest[j + offsX] = 0.f;
                            } else {
                                greenDest[j + offsX] = nonGreenDest0[j + offsX] = nonGreenDest1[j + offsX] = 1000.f + 25000.f * blend;
                            }
                        }

                        continue;
                    }

                    ngCentre = (*rawDataFrames[2 - offset])[i + 1][j - offset + 1];
                    ngRight = (*rawDataFrames[2 - (offset ^ 1)])[i + 1][j - (offset ^ 1) + 2];
                    ngLeft = (*rawDataFrames[2 - (offset ^ 1)])[i + 1][j - (offset ^ 1)];
                    ngTop = (*rawDataFrames[3 - ((offset << 1) + offset)])[i - 1][j - offset + 1];
                    ngBottom = (*rawDataFrames[3 - ((offset << 1) + offset)])[i + 1][j - offset + 1];
                    gridMax = nonGreenDiffCross(ngRight, ngLeft, ngTop, ngBottom, ngCentre, stddevFactorNonGreen2, eperIsoNonGreen2, nRead, prnu, showMotion);

                    if(gridMax > 0.f) {
                        if(showMotion) {
                            float blend = gridMax * blendFactor;

                            if(!showOnlyMask) {
                                // if showMotion is enabled colourize the pixel
                                nonGreenDest1[j + offsX] = 1000.f + 25000.f * blend;
                                nonGreenDest0[j + offsX] = greenDest[j + offsX] = 0.f;
                            } else {
                                greenDest[j + offsX] = nonGreenDest0[j + offsX] = nonGreenDest1[j + offsX] = 1000.f + 25000.f * blend;
                            }
                        }

                        continue;
                    }
                }

                if(adaptive && checkNonGreenHorizontal) {
                    float ngCentre = (*rawDataFrames[(offset << 1) + offset])[i][j + offset];
                    float ngRight = (*rawDataFrames[((offset ^ 1) << 1) + (offset ^ 1)])[i][j + (offset ^ 1) + 1];
                    float ngLeft = (*rawDataFrames[((offset ^ 1) << 1) + (offset ^ 1)])[i][j + (offset ^ 1) - 1];
                    float diffRight = ngRight - ngCentre;
                    float diff2 = ngLeft - ngCentre;

                    if(diffRight * diff2 >= 0.f) {
                        float val = (ngRight + ngLeft) / 2.f;
                        float gridMax = nonGreenDiff(ngCentre, val, stddevFactorNonGreen0, eperIsoNonGreen0, nRead, prnu, showMotion);

                        if(gridMax > 0.f) {
                            if(showMotion) {
                                float blend = gridMax * blendFactor;

                                if(!showOnlyMask) {
                                    // if showMotion is enabled colourize the pixel
                                    nonGreenDest0[j + offsX] = 1000.f + 25000.f * blend;
                                    nonGreenDest1[j + offsX] = greenDest[j + offsX] = 0.f;
                                } else {
                                    greenDest[j + offsX] = nonGreenDest0[j + offsX] = nonGreenDest1[j + offsX] = 1000.f + 25000.f * blend;
                                }
                            }

                            continue;
                        }
                    }

                    ngCentre = (*rawDataFrames[2 - offset])[i + 1][j - offset + 1];
                    ngRight = (*rawDataFrames[2 - (offset ^ 1)])[i + 1][j - (offset ^ 1) + 2];
                    ngLeft = (*rawDataFrames[2 - (offset ^ 1)])[i + 1][j - (offset ^ 1)];
                    diffRight = ngRight - ngCentre;
                    diff2 = ngLeft - ngCentre;

                    if(diffRight * diff2 >= 0.f) {
                        float val = (ngRight + ngLeft) / 2.f;
                        float gridMax = nonGreenDiff(ngCentre, val, stddevFactorNonGreen2, eperIsoNonGreen2, nRead, prnu, showMotion);

                        if(gridMax > 0.f) {
                            if(showMotion) {
                                float blend = gridMax * blendFactor;

                                if(!showOnlyMask) {
                                    // if showMotion is enabled colourize the pixel
                                    nonGreenDest1[j + offsX] = 1000.f + 25000.f * blend;
                                    nonGreenDest0[j + offsX] = greenDest[j + offsX] = 0.f;
                                } else {
                                    greenDest[j + offsX] = nonGreenDest0[j + offsX] = nonGreenDest1[j + offsX] = 1000.f + 25000.f * blend;
                                }
                            }

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
                        float val = (ngTop + ngBottom) / 2.f;
                        float gridMax = nonGreenDiff(ngCentre, val, stddevFactorNonGreen0, eperIsoNonGreen0, nRead, prnu, showMotion);

                        if(gridMax > 0.f) {
                            if(showMotion) {
                                float blend = gridMax * blendFactor;

                                if(!showOnlyMask) {
                                    // if showMotion is enabled colourize the pixel
                                    nonGreenDest0[j + offsX] = 1000.f + 25000.f * blend;
                                    nonGreenDest1[j + offsX] = greenDest[j + offsX] = 0.f;
                                } else {
                                    greenDest[j + offsX] = nonGreenDest0[j + offsX] = nonGreenDest1[j + offsX] = 1000.f + 25000.f * blend;
                                }
                            }

                            continue;
                        }
                    }

                    ngCentre = (*rawDataFrames[2 - offset])[i + 1][j - offset + 1];
                    ngTop = (*rawDataFrames[3 - ((offset << 1) + offset)])[i - 1][j - offset + 1];
                    ngBottom = (*rawDataFrames[3 - ((offset << 1) + offset)])[i + 1][j - offset + 1];

                    diffTop = ngTop - ngCentre;
                    diffBottom = ngBottom - ngCentre;

                    if(diffTop * diffBottom >= 0.f) {
                        float val = (ngTop + ngBottom) / 2.f;
                        float gridMax = nonGreenDiff(ngCentre, val, stddevFactorNonGreen2, eperIsoNonGreen2, nRead, prnu, showMotion);

                        if(gridMax > 0.f) {
                            if(showMotion) {
                                float blend = gridMax * blendFactor;

                                if(!showOnlyMask) {
                                    // if showMotion is enabled colourize the pixel
                                    nonGreenDest1[j + offsX] = 1000.f + 25000.f * blend;
                                    nonGreenDest0[j + offsX] = greenDest[j + offsX] = 0.f;
                                } else {
                                    greenDest[j + offsX] = nonGreenDest0[j + offsX] = nonGreenDest1[j + offsX] = 1000.f + 25000.f * blend;
                                }
                            }

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
