/*
*  This file is part of RawTherapee.
*
*  Copyright (c) 2012 Oliver Duis <www.oliverduis.de>
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
#include <cstring>

#include "dcp.h"
#include "safegtk.h"
#include "iccmatrices.h"
#include "iccstore.h"
#include "rawimagesource.h"
#include "improcfun.h"
#include "rt_math.h"

using namespace std;
using namespace rtengine;
using namespace rtexif;

static void Invert3x3(const double (*A)[3], double (*B)[3]) {

    double a00 = A[0][0];
    double a01 = A[0][1];
    double a02 = A[0][2];
    double a10 = A[1][0];
    double a11 = A[1][1];
    double a12 = A[1][2];
    double a20 = A[2][0];
    double a21 = A[2][1];
    double a22 = A[2][2];
    double temp [3][3];

    temp[0][0] = a11 * a22 - a21 * a12;
    temp[0][1] = a21 * a02 - a01 * a22;
    temp[0][2] = a01 * a12 - a11 * a02;
    temp[1][0] = a20 * a12 - a10 * a22;
    temp[1][1] = a00 * a22 - a20 * a02;
    temp[1][2] = a10 * a02 - a00 * a12;
    temp[2][0] = a10 * a21 - a20 * a11;
    temp[2][1] = a20 * a01 - a00 * a21;
    temp[2][2] = a00 * a11 - a10 * a01;

    double det = a00 * temp[0][0] + a01 * temp[1][0] + a02 * temp[2][0];

    if (fabs(det) < 1.0E-10) {
        abort(); // can't be inverted, we shouldn't be dealing with such matrices
    }

    for (int j = 0; j < 3; j++)	{
        for (int k = 0; k < 3; k++) {
            B[j][k] = temp[j][k] / det;
        }
    }
}

static void Multiply3x3(const double (*A)[3], const double (*B)[3], double (*C)[3]) {

    // use temp to support having output same as input
    double M[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            M[i][j] = 0;
            for (int k = 0; k < 3; k++) {
                M[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    memcpy(C, M, 3 * 3 * sizeof(double));
}

static void Multiply3x3_v3(const double (*A)[3], const double B[3], double C[3]) {

    // use temp to support having output same as input
    double M[3] = { 0, 0, 0 };
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            M[i] += A[i][j] * B[j];
        }
    }
    memcpy(C, M, 3 * sizeof(double));
}

static void Mix3x3(const double (*A)[3], double mulA, const double (*B)[3], double mulB, double (*C)[3]) {

    double M[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            M[i][j] = A[i][j] * mulA + B[i][j] * mulB;
        }
    }
    memcpy(C, M, 3 * 3 * sizeof(double));
}

static void MapWhiteMatrix(const double white1[3], const double white2[3], double (*B)[3]) {

    // code adapted from dng_color_spec::MapWhiteMatrix

    // Use the linearized Bradford adaptation matrix.
    double Mb[3][3] = { { 0.8951,  0.2664, -0.1614 }, { -0.7502,  1.7135,  0.0367 }, { 0.0389, -0.0685,  1.0296 }};

    double w1[3];
    Multiply3x3_v3(Mb, white1, w1);
    double w2[3];
    Multiply3x3_v3(Mb, white2, w2);

    // Negative white coordinates are kind of meaningless.
    w1[0] = std::max(w1[0], 0.0);
    w1[1] = std::max(w1[1], 0.0);
    w1[2] = std::max(w1[2], 0.0);
    w2[0] = std::max(w2[0], 0.0);
    w2[1] = std::max(w2[1], 0.0);
    w2[2] = std::max(w2[2], 0.0);

    // Limit scaling to something reasonable.
    double A[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    A[0][0] = std::max(0.1, std::min(w1[0] > 0.0 ? w2[0] / w1[0] : 10.0, 10.0));
    A[1][1] = std::max(0.1, std::min(w1[1] > 0.0 ? w2[1] / w1[1] : 10.0, 10.0));
    A[2][2] = std::max(0.1, std::min(w1[2] > 0.0 ? w2[2] / w1[2] : 10.0, 10.0));

    double temp[3][3];
    Invert3x3(Mb, temp);
    Multiply3x3(temp, A, temp);
    Multiply3x3(temp, Mb, B);
}

static void XYZtoXY(const double XYZ[3], double XY[2]) {
    double X = XYZ[0];
    double Y = XYZ[1];
    double Z = XYZ[2];
    double total = X + Y + Z;
    if (total > 0.0) {
        XY[0] = X / total;
        XY[1] = Y / total;
    } else {
        XY[0] = 0.3457;
        XY[1] = 0.3585;
    }
}

static void XYtoXYZ(const double XY[2], double XYZ[3]) {
    double temp[2] = { XY[0], XY[1] };
    // Restrict xy coord to someplace inside the range of real xy coordinates.
    // This prevents math from doing strange things when users specify
    // extreme temperature/tint coordinates.
    temp[0] = std::max(0.000001, std::min(temp[0], 0.999999));
    temp[1] = std::max(0.000001, std::min(temp[1], 0.999999));
    if (temp[0] + temp[1] > 0.999999) {
        double scale = 0.999999 / (temp[0] + temp[1]);
        temp[0] *= scale;
        temp[1] *= scale;
    }
    XYZ[0] = temp[0] / temp[1];
    XYZ[1] = 1.0;
    XYZ[2] = (1.0 - temp[0] - temp[1]) / temp[1];
}

enum dngCalibrationIlluminant {
    lsUnknown = 0,
    lsDaylight = 1,
    lsFluorescent = 2,
    lsTungsten = 3,
    lsFlash = 4,
    lsFineWeather = 9,
    lsCloudyWeather = 10,
    lsShade = 11,
    lsDaylightFluorescent = 12,	// D  5700 - 7100K
    lsDayWhiteFluorescent = 13,	// N  4600 - 5500K
    lsCoolWhiteFluorescent = 14, // W  3800 - 4500K
    lsWhiteFluorescent = 15, // WW 3250 - 3800K
    lsWarmWhiteFluorescent = 16, // L  2600 - 3250K
    lsStandardLightA = 17,
    lsStandardLightB = 18,
    lsStandardLightC = 19,
    lsD55 = 20,
    lsD65 = 21,
    lsD75 = 22,
    lsD50 = 23,
    lsISOStudioTungsten = 24,
    lsOther = 255
};

// should probably be moved to colortemp.cc
static double calibrationIlluminantToTemperature(int light) {

    // these temperatures are those found in DNG SDK reference code.
    switch (light) {
    case lsStandardLightA:
    case lsTungsten:
        return 2850.0;
    case lsISOStudioTungsten:
        return 3200.0;
    case lsD50:
        return 5000.0;
    case lsD55:
    case lsDaylight:
    case lsFineWeather:
    case lsFlash:
    case lsStandardLightB:
        return 5500.0;
    case lsD65:
    case lsStandardLightC:
    case lsCloudyWeather:
        return 6500.0;
    case lsD75:
    case lsShade:
        return 7500.0;
    case lsDaylightFluorescent:
        return (5700.0 + 7100.0) * 0.5;
    case lsDayWhiteFluorescent:
        return (4600.0 + 5500.0) * 0.5;
    case lsCoolWhiteFluorescent:
    case lsFluorescent:
        return (3800.0 + 4500.0) * 0.5;
    case lsWhiteFluorescent:
        return (3250.0 + 3800.0) * 0.5;
    case lsWarmWhiteFluorescent:
        return (2600.0 + 3250.0) * 0.5;
    default:
        return 0.0;
    }
}
void DCPProfile::MakeXYZCAM(ColorTemp &wb, double pre_mul[3], double camWbMatrix[3][3], int preferredIlluminant, double (*mXYZCAM)[3]) const
{
    // code adapted from dng_color_spec::FindXYZtoCamera
    // note that we do not support monochrome or colorplanes > 3 (no reductionMatrix support)
    // we do not support cameracalibration either

    double neutral[3]; // same as the DNG "AsShotNeutral" tag if white balance is Camera's own
    {
        /* A bit messy matrixing and conversions to get the neutral[] array from RT's own white balance which is stored in
           sRGB space, while the DCP code needs multipliers in CameraRGB space */
        double r, g, b;
        wb.getMultipliers(r, g, b);

        // camWbMatrix == imatrices.xyz_cam
        double cam_xyz[3][3];
        Invert3x3(camWbMatrix, cam_xyz);
        double cam_rgb[3][3];
        Multiply3x3(cam_xyz, xyz_sRGB, cam_rgb);
        double camwb_red   = cam_rgb[0][0]*r + cam_rgb[0][1]*g + cam_rgb[0][2]*b;
        double camwb_green = cam_rgb[1][0]*r + cam_rgb[1][1]*g + cam_rgb[1][2]*b;
        double camwb_blue  = cam_rgb[2][0]*r + cam_rgb[2][1]*g + cam_rgb[2][2]*b;
        neutral[0] = camwb_red / pre_mul[0];
        neutral[1] = camwb_green / pre_mul[1];
        neutral[2] = camwb_blue / pre_mul[2];
        double maxentry = 0;
        for (int i = 0; i < 3; i++) {
            if (neutral[i] > maxentry) {
                maxentry = neutral[i];
            }
        }
        for (int i = 0; i < 3; i++) {
            neutral[i] /= maxentry;
        }
    }

    /* Calculate what the RGB multipliers corresponds to as a white XY coordinate, based on the
       DCP ColorMatrix or ColorMatrices if dual-illuminant. This is the DNG reference code way to
       do it, which is a bit different from RT's own white balance model at the time of writing.
       When RT's white balance can make use of the DCP color matrices we could use that instead. */
    double white_xy[2];
    dngref_NeutralToXY(neutral, preferredIlluminant, white_xy);

    bool hasFwd1 = hasForwardMatrix1;
    bool hasFwd2 = hasForwardMatrix2;
    bool hasCol1 = hasColorMatrix1;
    bool hasCol2 = hasColorMatrix2;
    if (preferredIlluminant == 1) {
        if (hasFwd1) hasFwd2 = false;
        if (hasCol1) hasCol2 = false;
    } else if (preferredIlluminant == 2) {
        if (hasFwd2) hasFwd1 = false;
        if (hasCol2) hasCol1 = false;
    }

    // mix if we have two matrices
    double mix = 1.0;
    if ((hasCol1 && hasCol2) || (hasFwd1 && hasFwd2)) {
        double wbtemp;
        /* DNG ref way to convert XY to temperature, which affect matrix mixing. A different model here
           typically does not affect the result too much, ie it's probably not strictly necessary to
           use the DNG reference code here, but we do it for now. */
        dngref_XYCoord2Temperature(white_xy, &wbtemp, NULL);
        if (wbtemp <= temperature1) {
            mix = 1.0;
        } else if (wbtemp >= temperature2) {
            mix = 0.0;
        } else {
            double invT = 1.0 / wbtemp;
            mix = (invT - (1.0 / temperature2)) / ((1.0 / temperature1) - (1.0 / temperature2));
        }
    }

    // Colormatrix
    double mCol[3][3];
    if (hasCol1 && hasCol2) {
        // interpolate
        if (mix >= 1.0) {
            memcpy(mCol, mColorMatrix1, sizeof(mCol));
        } else if (mix <= 0.0) {
            memcpy(mCol, mColorMatrix2, sizeof(mCol));
        } else {
            Mix3x3(mColorMatrix1, mix, mColorMatrix2, 1.0 - mix, mCol);
        }
    } else if (hasCol1) {
        memcpy(mCol, mColorMatrix1, sizeof(mCol));
    } else {
        memcpy(mCol, mColorMatrix2, sizeof(mCol));
    }

    /*
      The exact position of the white XY coordinate affects the result very much, thus
      it's important that the result is very similar or the same as DNG reference code.
      Especially important is it that the raw-embedded "AsShot" multipliers is translated
      to the same white XY coordinate as the DNG reference code, or else third party DCPs
      will show incorrect color.
    */

    double white_xyz[3];
    XYtoXYZ(white_xy, white_xyz);

    double cam_xyz[3][3];
    if (hasFwd1 || hasFwd2) {
        // always prefer ForwardMatrix ahead of ColorMatrix
        double mFwd[3][3];
        if (hasFwd1 && hasFwd2) {
            // interpolate
            if (mix >= 1.0) {
                memcpy(mFwd, mForwardMatrix1, sizeof(mFwd));
            } else if (mix <= 0.0) {
                memcpy(mFwd, mForwardMatrix2, sizeof(mFwd));
            } else {
                Mix3x3(mForwardMatrix1, mix, mForwardMatrix2, 1.0 - mix, mFwd);
            }
        } else if (hasFwd1) {
            memcpy(mFwd, mForwardMatrix1, sizeof(mFwd));
        } else {
            memcpy(mFwd, mForwardMatrix2, sizeof(mFwd));
        }
        // adapted from dng_color_spec::SetWhiteXY
        double CameraWhite[3];
        Multiply3x3_v3(mCol, white_xyz, CameraWhite);

        double whiteDiag[3][3] = {{CameraWhite[0], 0, 0}, {0, CameraWhite[1], 0}, {0, 0, CameraWhite[2]}};
        double whiteDiagInv[3][3];
        Invert3x3(whiteDiag, whiteDiagInv);

        double xyz_cam[3][3];
        Multiply3x3(mFwd, whiteDiagInv, xyz_cam);
        Invert3x3(xyz_cam, cam_xyz);
    } else {
        double whiteMatrix[3][3];
        const double white_d50[3] = { 0.3457, 0.3585, 0.2958 }; // D50
        MapWhiteMatrix(white_d50, white_xyz, whiteMatrix);
        Multiply3x3(mCol, whiteMatrix, cam_xyz);
    }

    // convert cam_xyz (XYZ D50 to CameraRGB, "PCS to Camera" in DNG terminology) to mXYZCAM

    {
        // This block can probably be simplified, seems unnecessary to pass through the sRGB matrix
        // (probably dcraw legacy), it does no harm though as we don't clip anything.
        int i,j,k;

        // Multiply out XYZ colorspace
        double cam_rgb[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (i=0; i < 3; i++)
            for (j=0; j < 3; j++)
                for (k=0; k < 3; k++)
                    cam_rgb[i][j] += cam_xyz[i][k] * xyz_sRGB[k][j];

        // Normalize cam_rgb so that:  cam_rgb * (1,1,1) is (1,1,1,1)
        double num;
        for (i=0; i<3; i++) {
            for (num=j=0; j<3; j++) num += cam_rgb[i][j];
            for (j=0; j<3; j++) cam_rgb[i][j] /= num;
        }

        double rgb_cam[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        RawImageSource::inverse33 (cam_rgb, rgb_cam);

        for (i=0; i<3; i++)
            for (j=0; j<3; j++) mXYZCAM[i][j]=0;
        for (i=0; i<3; i++)
            for (j=0; j<3; j++)
                for (k=0; k<3; k++)
                    mXYZCAM[i][j] += xyz_sRGB[i][k] * rgb_cam[k][j];
    }
}

const DCPProfile::HSBModify* DCPProfile::MakeHueSatMap(ColorTemp &wb, int preferredIlluminant, HSBModify **deleteHandle) const {

    *deleteHandle = NULL;
    if (!aDeltas1) {
        return NULL;
    }
    if (!aDeltas2) {
        return aDeltas1;
    }

    if (preferredIlluminant == 1) {
        return aDeltas1;
    } else if (preferredIlluminant == 2) {
        return aDeltas2;
    }

    // Interpolate based on color temperature.
    if (temperature1 <= 0.0 || temperature2 <= 0.0 || temperature1 == temperature2) {
        return aDeltas1;
    }
    bool reverseOrder = temperature1 > temperature2;
    double t1, t2;
    if (reverseOrder) {
        t1 = temperature2;
        t2 = temperature1;
    } else {
        t1 = temperature1;
        t2 = temperature2;
    }

    double mix;
    if (wb.getTemp() <= t1) {
        mix = 1.0;
    } else if (wb.getTemp() >= t2) {
        mix = 0.0;
    } else {
        double invT = 1.0 / wb.getTemp();
        mix = (invT - (1.0 / t2)) / ((1.0 / t1) - (1.0 / t2));
    }

    if (reverseOrder) {
        mix = 1.0 - mix;
    }

    if (mix >= 1.0) {
        return aDeltas1;
    } else if (mix <= 0.0) {
        return aDeltas2;
    }

    // Interpolate between the tables.
    HSBModify *aDeltas = new HSBModify[DeltaInfo.iArrayCount];
    *deleteHandle = aDeltas;
    float w1 = (float)mix;
    float w2 = 1.0f - (float)mix;
    for (int i = 0; i < DeltaInfo.iArrayCount; i++) {
        aDeltas[i].fHueShift = w1 * aDeltas1[i].fHueShift + w2 * aDeltas2[i].fHueShift;
        aDeltas[i].fSatScale = w1 * aDeltas1[i].fSatScale + w2 * aDeltas2[i].fSatScale;
        aDeltas[i].fValScale = w1 * aDeltas1[i].fValScale + w2 * aDeltas2[i].fValScale;
    }
    return aDeltas;
}

DCPProfile::DCPProfile(Glib::ustring fname, bool isRTProfile) {
    const int TIFFFloatSize=4;
    const int TagColorMatrix1=50721, TagColorMatrix2=50722, TagProfileHueSatMapDims=50937;
    const int TagForwardMatrix1=50964, TagForwardMatrix2=50965;
    const int TagProfileHueSatMapData1=50938, TagProfileHueSatMapData2=50939;
    const int TagCalibrationIlluminant1=50778, TagCalibrationIlluminant2=50779;
    const int TagProfileLookTableData=50982, TagProfileLookTableDims=50981;  // ProfileLookup is the low quality variant
    const int TagProfileToneCurve=50940;

    aDeltas1=aDeltas2=aLookTable=NULL;

    FILE *pFile = safe_g_fopen(fname, "rb");

    TagDirectory *tagDir=ExifManager::parseTIFF(pFile, false);

    Tag* tag = tagDir->getTag(TagCalibrationIlluminant1); iLightSource1 = (tag!=NULL ? tag->toInt(0,rtexif::SHORT) : -1);
    tag = tagDir->getTag(TagCalibrationIlluminant2); iLightSource2 = (tag!=NULL ? tag->toInt(0,rtexif::SHORT) : -1);
    temperature1 = calibrationIlluminantToTemperature(iLightSource1);
    temperature2 = calibrationIlluminantToTemperature(iLightSource2);

    bool hasSecondHueSat = tagDir->getTag(TagProfileHueSatMapData2)!=NULL;  // some profiles have two matrices, but just one huesat

    // Fetch Forward Matrices, if any
    hasForwardMatrix1 = false;
    hasForwardMatrix2 = false;
    hasColorMatrix1 = false;
    hasColorMatrix2 = false;
    hasToneCurve = false;
    tag = tagDir->getTag(TagForwardMatrix1);
    if (tag) {
        hasForwardMatrix1 = true;
        for (int row=0;row<3;row++) { 
            for (int col=0;col<3;col++) {
                mForwardMatrix1[row][col]=(float)tag->toDouble((col+row*3)*8);
            }
        }
    }
    tag = tagDir->getTag(TagForwardMatrix2);
    if (tag) {
        hasForwardMatrix2 = true;
        for (int row=0;row<3;row++) { 
            for (int col=0;col<3;col++) {
                mForwardMatrix2[row][col]=(float)tag->toDouble((col+row*3)*8);
            }
        }
    }

    // Color Matrix (1 is always there)
    tag = tagDir->getTag(TagColorMatrix1);
    if (!tag) {
        // FIXME: better error handling
        fprintf(stderr, "Bad DCP, no ColorMatrix1\n");
        abort();
    }
    hasColorMatrix1 = true;

    for (int row=0;row<3;row++) { 
        for (int col=0;col<3;col++) {
            mColorMatrix1[row][col]=(float)tag->toDouble((col+row*3)*8);
        }
    }

    tag=tagDir->getTag(TagProfileLookTableDims);
    if (tag!=NULL) {
        LookInfo.iHueDivisions=tag->toInt(0); LookInfo.iSatDivisions=tag->toInt(4); LookInfo.iValDivisions=tag->toInt(8);

        tag = tagDir->getTag(TagProfileLookTableData);
        LookInfo.iArrayCount = tag->getCount()/3;

        aLookTable =new HSBModify[LookInfo.iArrayCount];

        for (int i=0;i<LookInfo.iArrayCount;i++) {
            aLookTable[i].fHueShift=tag->toDouble((i*3)*TIFFFloatSize);
            aLookTable[i].fSatScale=tag->toDouble((i*3+1)*TIFFFloatSize);
            aLookTable[i].fValScale=tag->toDouble((i*3+2)*TIFFFloatSize);
        }

        // precalculated constants for table application
        LookInfo.pc.hScale = (LookInfo.iHueDivisions < 2) ? 0.0f : (LookInfo.iHueDivisions * (1.0f / 6.0f));
        LookInfo.pc.sScale = (float) (LookInfo.iSatDivisions - 1);
        LookInfo.pc.vScale = (float) (LookInfo.iValDivisions - 1);
        LookInfo.pc.maxHueIndex0 = LookInfo.iHueDivisions - 1;
        LookInfo.pc.maxSatIndex0 = LookInfo.iSatDivisions - 2;
        LookInfo.pc.maxValIndex0 = LookInfo.iValDivisions - 2;
        LookInfo.pc.hueStep = LookInfo.iSatDivisions;
        LookInfo.pc.valStep = LookInfo.iHueDivisions * LookInfo.pc.hueStep;
    }

    tag = tagDir->getTag(TagProfileHueSatMapDims);
    if (tag!=NULL) {
        DeltaInfo.iHueDivisions=tag->toInt(0); DeltaInfo.iSatDivisions=tag->toInt(4); DeltaInfo.iValDivisions=tag->toInt(8);

        tag = tagDir->getTag(TagProfileHueSatMapData1);
        DeltaInfo.iArrayCount = tag->getCount()/3;

        aDeltas1=new HSBModify[DeltaInfo.iArrayCount];

        for (int i=0;i<DeltaInfo.iArrayCount;i++) {
            aDeltas1[i].fHueShift=tag->toDouble((i*3)*TIFFFloatSize);
            aDeltas1[i].fSatScale=tag->toDouble((i*3+1)*TIFFFloatSize);
            aDeltas1[i].fValScale=tag->toDouble((i*3+2)*TIFFFloatSize);
        }

        DeltaInfo.pc.hScale = (DeltaInfo.iHueDivisions < 2) ? 0.0f : (DeltaInfo.iHueDivisions * (1.0f / 6.0f));
        DeltaInfo.pc.sScale = (float) (DeltaInfo.iSatDivisions - 1);
        DeltaInfo.pc.vScale = (float) (DeltaInfo.iValDivisions - 1);
        DeltaInfo.pc.maxHueIndex0 = DeltaInfo.iHueDivisions - 1;
        DeltaInfo.pc.maxSatIndex0 = DeltaInfo.iSatDivisions - 2;
        DeltaInfo.pc.maxValIndex0 = DeltaInfo.iValDivisions - 2;
        DeltaInfo.pc.hueStep = DeltaInfo.iSatDivisions;
        DeltaInfo.pc.valStep = DeltaInfo.iHueDivisions * DeltaInfo.pc.hueStep;
    }

    if (iLightSource2!=-1) {
        // Second matrix
        tag = tagDir->getTag(TagColorMatrix2);
        hasColorMatrix2 = true;

        for (int row=0;row<3;row++) { 
            for (int col=0;col<3;col++) {
                mColorMatrix2[row][col]= (tag!=NULL ? (float)tag->toDouble((col+row*3)*8) : mColorMatrix1[row][col]);
            }
        }

        // Second huesatmap
        if (hasSecondHueSat) {
            aDeltas2=new HSBModify[DeltaInfo.iArrayCount];

            // Saturation maps. Need to be unwinded.
            tag = tagDir->getTag(TagProfileHueSatMapData2);

            for (int i=0;i<DeltaInfo.iArrayCount;i++) {
                aDeltas2[i].fHueShift=tag->toDouble((i*3)*TIFFFloatSize);
                aDeltas2[i].fSatScale=tag->toDouble((i*3+1)*TIFFFloatSize);
                aDeltas2[i].fValScale=tag->toDouble((i*3+2)*TIFFFloatSize);
            }
        }
    }

    // Read tone curve points, if any, but disable to RTs own profiles
    // the DCP tone curve is subjective and of low quality in comparison to RTs tone curves
    tag = tagDir->getTag(TagProfileToneCurve);
    if (tag!=NULL && !isRTProfile) {
        std::vector<double> cPoints;
        cPoints.push_back(double(DCT_Spline));  // The first value is the curve type

        // push back each X/Y coordinates in a loop
        bool curve_is_linear = true;
        for (int i=0;i<tag->getCount(); i+= 2) {
            double x = tag->toDouble((i+0)*TIFFFloatSize);
            double y = tag->toDouble((i+1)*TIFFFloatSize);
            if (x != y) {
                curve_is_linear = false;
            }
            cPoints.push_back( x );
            cPoints.push_back( y );
        }

        if (!curve_is_linear) {
            // Create the curve
            DiagonalCurve rawCurve(cPoints, CURVES_MIN_POLY_POINTS);

            toneCurve.Set((Curve*)&rawCurve);
            hasToneCurve = true;
        }
    }

    willInterpolate = false;
    if (hasForwardMatrix1) {
        if (hasForwardMatrix2) {
            if (memcmp(mForwardMatrix1, mForwardMatrix2, sizeof(mForwardMatrix1)) != 0) {
                // common that forward matrices are the same!
                willInterpolate = true;
            }
            if (aDeltas1 && aDeltas2) {
                // we assume tables are different
                willInterpolate = true;
            }
        }
    }
    if (hasColorMatrix1 && hasColorMatrix2) {
        if (memcmp(mColorMatrix1, mColorMatrix2, sizeof(mColorMatrix1)) != 0) {
            willInterpolate = true;
        }
        if (aDeltas1 && aDeltas2) {
            willInterpolate = true;
        }
    }

    if (pFile!=NULL) fclose(pFile);
    delete tagDir;
}

DCPProfile::~DCPProfile() {
    delete[] aDeltas1; delete[] aDeltas2;
}

void DCPProfile::HSDApply(const HSDTableInfo &ti, const HSBModify *tableBase, const float hs, const float ss, const float vs, float &h, float &s, float &v) const {

    // Apply the HueSatMap. Ported from Adobes reference implementation
    float hueShift, satScale, valScale;

    if (ti.iValDivisions < 2)  // Optimize most common case of "2.5D" table.
    {
        float hScaled = hs * ti.pc.hScale;
        float sScaled = ss * ti.pc.sScale;

        int hIndex0 = max((int)hScaled, 0);
        int sIndex0 = max(min((int)sScaled,ti.pc.maxSatIndex0),0);

        int hIndex1 = hIndex0 + 1;

        if (hIndex0 >= ti.pc.maxHueIndex0)
        {
            hIndex0 = ti.pc.maxHueIndex0;
            hIndex1 = 0;
        }

        float hFract1 = hScaled - (float) hIndex0;
        float sFract1 = sScaled - (float) sIndex0;

        float hFract0 = 1.0f - hFract1;
        float sFract0 = 1.0f - sFract1;

        const HSBModify *entry00 = tableBase + hIndex0 * ti.pc.hueStep + sIndex0;
        const HSBModify *entry01 = entry00 + (hIndex1 - hIndex0) * ti.pc.hueStep;

        float hueShift0 = hFract0 * entry00->fHueShift + hFract1 * entry01->fHueShift;
        float satScale0 = hFract0 * entry00->fSatScale + hFract1 * entry01->fSatScale;
        float valScale0 = hFract0 * entry00->fValScale + hFract1 * entry01->fValScale;

        entry00++;
        entry01++;

        float hueShift1 = hFract0 * entry00->fHueShift +
            hFract1 * entry01->fHueShift;

        float satScale1 = hFract0 * entry00->fSatScale +
            hFract1 * entry01->fSatScale;

        float valScale1 = hFract0 * entry00->fValScale +
            hFract1 * entry01->fValScale;

        hueShift = sFract0 * hueShift0 + sFract1 * hueShift1;
        satScale = sFract0 * satScale0 + sFract1 * satScale1;
        valScale = sFract0 * valScale0 + sFract1 * valScale1;

    } else {

        float hScaled = hs * ti.pc.hScale;
        float sScaled = ss * ti.pc.sScale;
        float vScaled = vs * ti.pc.vScale;

        int hIndex0 = (int) hScaled;
        int sIndex0 = max(min((int)sScaled,ti.pc.maxSatIndex0),0);
        int vIndex0 = max(min((int)vScaled,ti.pc.maxValIndex0),0);

        int hIndex1 = hIndex0 + 1;

        if (hIndex0 >= ti.pc.maxHueIndex0)
        {
            hIndex0 = ti.pc.maxHueIndex0;
            hIndex1 = 0;
        }

        float hFract1 = hScaled - (float) hIndex0;
        float sFract1 = sScaled - (float) sIndex0;
        float vFract1 = vScaled - (float) vIndex0;

        float hFract0 = 1.0f - hFract1;
        float sFract0 = 1.0f - sFract1;
        float vFract0 = 1.0f - vFract1;

        const HSBModify *entry00 = tableBase + vIndex0 * ti.pc.valStep + hIndex0 * ti.pc.hueStep + sIndex0;

        const HSBModify *entry01 = entry00 + (hIndex1 - hIndex0) * ti.pc.hueStep;

        const HSBModify *entry10 = entry00 + ti.pc.valStep;
        const HSBModify *entry11 = entry01 + ti.pc.valStep;

        float hueShift0 = vFract0 * (hFract0 * entry00->fHueShift +
                                     hFract1 * entry01->fHueShift) +
            vFract1 * (hFract0 * entry10->fHueShift +
                       hFract1 * entry11->fHueShift);

        float satScale0 = vFract0 * (hFract0 * entry00->fSatScale +
                                     hFract1 * entry01->fSatScale) +
            vFract1 * (hFract0 * entry10->fSatScale +
                       hFract1 * entry11->fSatScale);

        float valScale0 = vFract0 * (hFract0 * entry00->fValScale +
                                     hFract1 * entry01->fValScale) +
            vFract1 * (hFract0 * entry10->fValScale +
                       hFract1 * entry11->fValScale);

        entry00++;
        entry01++;
        entry10++;
        entry11++;

        float hueShift1 = vFract0 * (hFract0 * entry00->fHueShift +
                                     hFract1 * entry01->fHueShift) +
            vFract1 * (hFract0 * entry10->fHueShift +
                       hFract1 * entry11->fHueShift);

        float satScale1 = vFract0 * (hFract0 * entry00->fSatScale +
                                     hFract1 * entry01->fSatScale) +
            vFract1 * (hFract0 * entry10->fSatScale +
                       hFract1 * entry11->fSatScale);

        float valScale1 = vFract0 * (hFract0 * entry00->fValScale +
                                     hFract1 * entry01->fValScale) +
            vFract1 * (hFract0 * entry10->fValScale +
                       hFract1 * entry11->fValScale);

        hueShift = sFract0 * hueShift0 + sFract1 * hueShift1;
        satScale = sFract0 * satScale0 + sFract1 * satScale1;
        valScale = sFract0 * valScale0 + sFract1 * valScale1;
    }

    hueShift *= (6.0f / 360.0f);	// Convert to internal hue range.

    h += hueShift;
    s *= satScale;  // no clipping here, we are RT float :-)
    v *= valScale;
}

struct ruvt {
    double r;
    double u;
    double v;
    double t;
};

static const double kTintScale = -3000.0;
static const ruvt kTempTable [] =
	{
	{   0, 0.18006, 0.26352, -0.24341 },
	{  10, 0.18066, 0.26589, -0.25479 },
	{  20, 0.18133, 0.26846, -0.26876 },
	{  30, 0.18208, 0.27119, -0.28539 },
	{  40, 0.18293, 0.27407, -0.30470 },
	{  50, 0.18388, 0.27709, -0.32675 },
	{  60, 0.18494, 0.28021, -0.35156 },
	{  70, 0.18611, 0.28342, -0.37915 },
	{  80, 0.18740, 0.28668, -0.40955 },
	{  90, 0.18880, 0.28997, -0.44278 },
	{ 100, 0.19032, 0.29326, -0.47888 },
	{ 125, 0.19462, 0.30141, -0.58204 },
	{ 150, 0.19962, 0.30921, -0.70471 },
	{ 175, 0.20525, 0.31647, -0.84901 },
	{ 200, 0.21142, 0.32312, -1.0182 },
	{ 225, 0.21807, 0.32909, -1.2168 },
	{ 250, 0.22511, 0.33439, -1.4512 },
	{ 275, 0.23247, 0.33904, -1.7298 },
	{ 300, 0.24010, 0.34308, -2.0637 },
	{ 325, 0.24702, 0.34655, -2.4681 },
	{ 350, 0.25591, 0.34951, -2.9641 },
	{ 375, 0.26400, 0.35200, -3.5814 },
	{ 400, 0.27218, 0.35407, -4.3633 },
	{ 425, 0.28039, 0.35577, -5.3762 },
	{ 450, 0.28863, 0.35714, -6.7262 },
	{ 475, 0.29685, 0.35823, -8.5955 },
	{ 500, 0.30505, 0.35907, -11.324 },
	{ 525, 0.31320, 0.35968, -15.628 },
	{ 550, 0.32129, 0.36011, -23.325 },
	{ 575, 0.32931, 0.36038, -40.770 },
	{ 600, 0.33724, 0.36051, -116.45 }
	};

void DCPProfile::dngref_XYCoord2Temperature(const double whiteXY[2], double *temp, double *tint) const {
    double fTemperature = 0;
    double fTint = 0;

    // Convert to uv space.
    double u = 2.0 * whiteXY[0] / (1.5 - whiteXY[0] + 6.0 * whiteXY[1]);
    double v = 3.0 * whiteXY[1] / (1.5 - whiteXY[0] + 6.0 * whiteXY[1]);

    // Search for line pair coordinate is between.
    double last_dt = 0.0;
    double last_dv = 0.0;
    double last_du = 0.0;

    for (uint32_t index = 1; index <= 30; index++) {
        // Convert slope to delta-u and delta-v, with length 1.
        double du = 1.0;
        double dv = kTempTable [index] . t;
        double len = sqrt (1.0 + dv * dv);
        du /= len;
        dv /= len;

        // Find delta from black body point to test coordinate.
        double uu = u - kTempTable [index] . u;
        double vv = v - kTempTable [index] . v;

        // Find distance above or below line.
        double dt = - uu * dv + vv * du;

        // If below line, we have found line pair.
        if (dt <= 0.0 || index == 30) {
            // Find fractional weight of two lines.
            if (dt > 0.0)
                dt = 0.0;
            dt = -dt;
            double f;
            if (index == 1)
            {
                f = 0.0;
            }
            else
            {
                f = dt / (last_dt + dt);
            }

            // Interpolate the temperature.
            fTemperature = 1.0E6 / (kTempTable [index - 1] . r * f +
                                    kTempTable [index    ] . r * (1.0 - f));

            // Find delta from black body point to test coordinate.
            uu = u - (kTempTable [index - 1] . u * f +
                      kTempTable [index    ] . u * (1.0 - f));
            vv = v - (kTempTable [index - 1] . v * f +
                      kTempTable [index    ] . v * (1.0 - f));
            // Interpolate vectors along slope.
            du = du * (1.0 - f) + last_du * f;
            dv = dv * (1.0 - f) + last_dv * f;
            len = sqrt (du * du + dv * dv);
            du /= len;
            dv /= len;

            // Find distance along slope.
            fTint = (uu * du + vv * dv) * kTintScale;
            break;
        }
        // Try next line pair.
        last_dt = dt;
        last_du = du;
        last_dv = dv;
    }
    if (temp != NULL)
        *temp = fTemperature;
    if (tint != NULL)
        *tint = fTint;
}

void DCPProfile::dngref_FindXYZtoCamera(const double whiteXY[2], int preferredIlluminant, double (*xyzToCamera)[3]) const {

    bool hasCol1 = hasColorMatrix1;
    bool hasCol2 = hasColorMatrix2;
    if (preferredIlluminant == 1) {
        if (hasCol1) hasCol2 = false;
    } else if (preferredIlluminant == 2) {
        if (hasCol2) hasCol1 = false;
    }

    // mix if we have two matrices
    double mix;
    if (hasCol1 && hasCol2) {
        double wbtemp;
        /*
          Note: we're using DNG SDK reference code for XY to temperature translation to get the exact same mix as
          the reference code does.
        */
        dngref_XYCoord2Temperature(whiteXY, &wbtemp, NULL);
        if (wbtemp <= temperature1) {
            mix = 1.0;
        } else if (wbtemp >= temperature2) {
            mix = 0.0;
        } else {
            double invT = 1.0 / wbtemp;
            mix = (invT - (1.0 / temperature2)) / ((1.0 / temperature1) - (1.0 / temperature2));
        }
    }

    // Interpolate the color matrix.
    double mCol[3][3];
    if (hasCol1 && hasCol2) {
        // interpolate
        if (mix >= 1.0) {
            memcpy(mCol, mColorMatrix1, sizeof(mCol));
        } else if (mix <= 0.0) {
            memcpy(mCol, mColorMatrix2, sizeof(mCol));
        } else {
            Mix3x3(mColorMatrix1, mix, mColorMatrix2, 1.0 - mix, mCol);
        }
    } else if (hasCol1) {
        memcpy(mCol, mColorMatrix1, sizeof(mCol));
    } else {
        memcpy(mCol, mColorMatrix2, sizeof(mCol));
    }
    memcpy(xyzToCamera, mCol, sizeof(mCol));
}

void DCPProfile::dngref_NeutralToXY(double neutral[3], int preferredIlluminant, double XY[2]) const {
    const int kMaxPasses = 30;
    double lastXY[2] = { 0.3457, 0.3585 }; // D50
    for (int pass = 0; pass < kMaxPasses; pass++) {
        double xyzToCamera[3][3];
        dngref_FindXYZtoCamera(lastXY, preferredIlluminant, xyzToCamera);

        double invM[3][3], nextXYZ[3], nextXY[2];
        Invert3x3(xyzToCamera, invM);
        Multiply3x3_v3(invM, neutral, nextXYZ);
        XYZtoXY(nextXYZ, nextXY);

        if (fabs(nextXY[0] - lastXY[0]) +
            fabs(nextXY[1] - lastXY[1]) < 0.0000001)
        {
            XY[0] = nextXY[0];
            XY[1] = nextXY[1];
            return;
        }
        // If we reach the limit without converging, we are most likely
        // in a two value oscillation.  So take the average of the last
        // two estimates and give up.
        if (pass == kMaxPasses - 1) {
            nextXY[0] = (lastXY[0] + nextXY[0]) * 0.5;
            nextXY[1] = (lastXY[1] + nextXY[1]) * 0.5;
        }
        lastXY[0] = nextXY[0];
        lastXY[1] = nextXY[1];
    }
    XY[0] = lastXY[0];
    XY[1] = lastXY[1];
}

void DCPProfile::Apply(Imagefloat *pImg, int preferredIlluminant, Glib::ustring workingSpace, ColorTemp &wb, double pre_mul[3], double camWbMatrix[3][3], float rawWhiteFac, bool useToneCurve) const {

    TMatrix mWork = iccStore->workingSpaceInverseMatrix (workingSpace);

    double mXYZCAM[3][3]; // Camera RGB to XYZ D50 matrix
    MakeXYZCAM(wb, pre_mul, camWbMatrix, preferredIlluminant, mXYZCAM);
    HSBModify *deleteTableHandle;
    const HSBModify *deltaBase = MakeHueSatMap(wb, preferredIlluminant, &deleteTableHandle);

    useToneCurve&=toneCurve;

    if (deltaBase == NULL && aLookTable == NULL && !useToneCurve) {
        //===== The fast path: no LUT and not tone curve- Calculate matrix for direct conversion raw>working space
        double mat[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (int i=0; i<3; i++) 
            for (int j=0; j<3; j++)
                for (int k=0; k<3; k++)
                    mat[i][j] += mWork[i][k] * mXYZCAM[k][j];

        // Apply the matrix part
#pragma omp parallel for
        for (int y=0; y<pImg->height; y++) {
            float newr, newg, newb;
            for (int x=0; x<pImg->width; x++) {
                newr = mat[0][0]*pImg->r(y,x) + mat[0][1]*pImg->g(y,x) + mat[0][2]*pImg->b(y,x);
                newg = mat[1][0]*pImg->r(y,x) + mat[1][1]*pImg->g(y,x) + mat[1][2]*pImg->b(y,x);
                newb = mat[2][0]*pImg->r(y,x) + mat[2][1]*pImg->g(y,x) + mat[2][2]*pImg->b(y,x);

                pImg->r(y,x) = newr; pImg->g(y,x) = newg; pImg->b(y,x) = newb;
            }
        }
    }
    else {
        //===== LUT available- Calculate matrix for conversion raw>ProPhoto
        double m2ProPhoto[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (int i=0; i<3; i++) 
            for (int j=0; j<3; j++)
                for (int k=0; k<3; k++)
                    m2ProPhoto[i][j] += prophoto_xyz[i][k] * mXYZCAM[k][j];

        double m2Work[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (int i=0; i<3; i++) 
            for (int j=0; j<3; j++)
                for (int k=0; k<3; k++)
                    m2Work[i][j] += mWork[i][k] * xyz_prophoto[k][j];

        bool useRawWhite=fabs(rawWhiteFac)>0.001;

        // Convert to prophoto and apply LUT
#pragma omp parallel for
        for (int y=0; y<pImg->height; y++) {
            float newr, newg, newb, h,s,v,hs,ss,vs;
            for (int x=0; x<pImg->width; x++) {
                newr = m2ProPhoto[0][0]*pImg->r(y,x) + m2ProPhoto[0][1]*pImg->g(y,x) + m2ProPhoto[0][2]*pImg->b(y,x);
                newg = m2ProPhoto[1][0]*pImg->r(y,x) + m2ProPhoto[1][1]*pImg->g(y,x) + m2ProPhoto[1][2]*pImg->b(y,x);
                newb = m2ProPhoto[2][0]*pImg->r(y,x) + m2ProPhoto[2][1]*pImg->g(y,x) + m2ProPhoto[2][2]*pImg->b(y,x);

                // if point is in negative area, just the matrix, but not the LUT
                if ((deltaBase || aLookTable) && newr>=0 && newg>=0 && newb>=0) {
                    Color::rgb2hsv(newr, newg, newb, h , s, v);
                    h*=6.f;  // RT calculates in [0,1]

                    if (useRawWhite) {
                        // Retro-calculate what the point was like before RAW white came in
                        Color::rgb2hsv(newr/rawWhiteFac, newg/rawWhiteFac, newb/rawWhiteFac, hs, ss, vs);
                        hs*=6.f;  // RT calculates in [0,1]
                    } else {
                        hs=h; ss=s; vs=v;
                    }

                    if (deltaBase) {
                        HSDApply(DeltaInfo, deltaBase, hs, ss, vs, h, s, v);
                    }
                    if (aLookTable) {
                        HSDApply(LookInfo, aLookTable, hs, ss, vs, h, s, v);
                    }

                    // RT range correction
                    if (h < 0.0f) h += 6.0f;
                    if (h >= 6.0f) h -= 6.0f;
                    h/=6.f;  
                    Color::hsv2rgb( h, s, v, newr, newg, newb);
                }
                // tone curve
                if (useToneCurve) toneCurve.Apply(newr, newg, newb);

                pImg->r(y,x) = m2Work[0][0]*newr + m2Work[0][1]*newg + m2Work[0][2]*newb;
                pImg->g(y,x) = m2Work[1][0]*newr + m2Work[1][1]*newg + m2Work[1][2]*newb;
                pImg->b(y,x) = m2Work[2][0]*newr + m2Work[2][1]*newg + m2Work[2][2]*newb;
            }
        }
    }

    if (deleteTableHandle) delete[] deleteTableHandle;
}

// Integer variant is legacy, only used for thumbs. Simply take the matrix here
void DCPProfile::Apply(Image16 *pImg, int preferredIlluminant, Glib::ustring workingSpace, ColorTemp &wb, double pre_mul[3], double camWbMatrix[3][3], bool useToneCurve) const {
    TMatrix mWork = iccStore->workingSpaceInverseMatrix (workingSpace);

    double mXYZCAM[3][3];
    MakeXYZCAM(wb, pre_mul, camWbMatrix, preferredIlluminant, mXYZCAM);

    useToneCurve&=toneCurve;

    // Calculate matrix for direct conversion raw>working space
    double mat[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int i=0; i<3; i++) 
        for (int j=0; j<3; j++)
            for (int k=0; k<3; k++)
                mat[i][j] += mWork[i][k] * mXYZCAM[k][j];

    // Apply the matrix part
#pragma omp parallel for
    for (int y=0; y<pImg->height; y++) {
        float newr, newg, newb;
        for (int x=0; x<pImg->width; x++) {
            newr = mat[0][0]*pImg->r(y,x) + mat[0][1]*pImg->g(y,x) + mat[0][2]*pImg->b(y,x);
            newg = mat[1][0]*pImg->r(y,x) + mat[1][1]*pImg->g(y,x) + mat[1][2]*pImg->b(y,x);
            newb = mat[2][0]*pImg->r(y,x) + mat[2][1]*pImg->g(y,x) + mat[2][2]*pImg->b(y,x);

            // tone curve
            if (useToneCurve) toneCurve.Apply(newr, newg, newb);

            pImg->r(y,x) = CLIP((int)newr); pImg->g(y,x) = CLIP((int)newg); pImg->b(y,x) = CLIP((int)newb);
        }
    }
}


// Generates as singleton
DCPStore* DCPStore::getInstance()
{
    static DCPStore* instance_ = 0;
    if ( instance_ == 0 )
    {
        static MyMutex smutex_;
        MyMutex::MyLock lock(smutex_);
        if ( instance_ == 0 )
        {
            instance_ = new DCPStore();
        }
    }
    return instance_;
}

// Reads all profiles from the given profiles dir
void DCPStore::init (Glib::ustring rtProfileDir) {
    MyMutex::MyLock lock(mtx);

    fileStdProfiles.clear();

    Glib::ustring rootDirName=rtProfileDir;

    if (rootDirName!="") {
        std::deque<Glib::ustring> qDirs;

        qDirs.push_front(rootDirName);

        while (!qDirs.empty()) {
            // process directory
            Glib::ustring dirname = qDirs.back();
            qDirs.pop_back();

            Glib::Dir* dir = NULL;
            try {
                if (!safe_file_test (dirname, Glib::FILE_TEST_IS_DIR)) return;
                dir = new Glib::Dir (dirname);
            }
            catch (Glib::Exception& fe) {
                return;
            }
            dirname = dirname + "/";
            for (Glib::DirIterator i = dir->begin(); i!=dir->end(); ++i) {
                Glib::ustring fname = dirname + *i;
                Glib::ustring sname = *i;
                // ignore directories
                if (!safe_file_test (fname, Glib::FILE_TEST_IS_DIR)) {
                    size_t lastdot = sname.find_last_of ('.');
                    if (lastdot!=Glib::ustring::npos && lastdot<=sname.size()-4 && (!sname.casefold().compare (lastdot, 4, ".dcp"))) {
                        Glib::ustring camShortName = sname.substr(0,lastdot).uppercase();
                        fileStdProfiles[camShortName]=fname;  // they will be loaded and cached on demand
                    }
                } else qDirs.push_front(fname);  // for later scanning
            }
            delete dir;
        }
    }
}

DCPProfile* DCPStore::getProfile (Glib::ustring filename, bool isRTProfile) {
    MyMutex::MyLock lock(mtx);

    std::map<Glib::ustring, DCPProfile*>::iterator r = profileCache.find (filename);
    if (r!=profileCache.end()) return r->second;

    // Add profile
    profileCache[filename]=new DCPProfile(filename, isRTProfile);

    return profileCache[filename];
}

DCPProfile* DCPStore::getStdProfile(Glib::ustring camShortName) {
    Glib::ustring name2=camShortName.uppercase();

    // Warning: do NOT use map.find(), since it does not seem to work reliably here
    for (std::map<Glib::ustring, Glib::ustring>::iterator i=fileStdProfiles.begin();i!=fileStdProfiles.end();i++)
        if (name2==(*i).first) return getProfile((*i).second, true);

    return NULL;
}

bool DCPStore::isValidDCPFileName(Glib::ustring filename) const {
    if (!safe_file_test (filename, Glib::FILE_TEST_EXISTS) || safe_file_test (filename, Glib::FILE_TEST_IS_DIR)) return false;
    size_t pos=filename.find_last_of ('.');
    return pos>0 && (!filename.casefold().compare (pos, 4, ".dcp") || !filename.casefold().compare (pos, 4, ".dng"));
}
