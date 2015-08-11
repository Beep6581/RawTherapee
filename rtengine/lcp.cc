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

#include "lcp.h"
#include "safegtk.h"
#include "iccmatrices.h"
#include "iccstore.h"
#include "rawimagesource.h"
#include "improcfun.h"
#include "rt_math.h"

#ifdef WIN32
#include <windows.h>
// for GCC32
#ifndef _WIN32_IE
#define _WIN32_IE 0x0600
#endif
#include <shlobj.h>
#endif


using namespace std;
using namespace rtengine;
using namespace rtexif;


LCPModelCommon::LCPModelCommon()
{
    focLenX = focLenY = -1;
    imgXCenter = imgYCenter = 0.5;
    x0 = y0 = fx = fy = meanErr = 0;
    badErr = false;

    for (int i = 0; i < 5; i++) {
        param[i] = 0;
    }

    scaleFac = 1;
}

bool LCPModelCommon::empty() const
{
    return param[0] == 0 && param[1] == 0 && param[2] == 0;
}

void LCPModelCommon::print() const
{
    printf("focLen %g/%g; imgCenter %g/%g; scale %g; err %g\n", focLenX, focLenY, imgXCenter, imgYCenter, scaleFac, meanErr);
    printf("xy0 %g/%g  fxy %g/%g\n", x0, y0, fx, fy);
    printf("param: %g/%g/%g/%g/%g\n", param[0], param[1], param[2], param[3], param[4]);
}

// weightened merge two parameters
void LCPModelCommon::merge(const LCPModelCommon& a, const LCPModelCommon& b, float facA)
{
    float facB = 1 - facA;

    focLenX    = facA * a.focLenX    + facB * b.focLenX;
    focLenY    = facA * a.focLenY    + facB * b.focLenY;
    imgXCenter = facA * a.imgXCenter + facB * b.imgXCenter;
    imgYCenter = facA * a.imgYCenter + facB * b.imgYCenter;
    scaleFac   = facA * a.scaleFac   + facB * b.scaleFac;
    meanErr    = facA * a.meanErr    + facB * b.meanErr;

    for (int i = 0; i < 5; i++) {
        param[i] = facA * a.param[i] + facB * b.param[i];
    }
}

void LCPModelCommon::prepareParams(int fullWidth, int fullHeight, float focalLength, float focalLength35mm, float sensorFormatFactor, bool swapXY, bool mirrorX, bool mirrorY)
{
    // Mention that the Adobe technical paper has a bug here, the DMAX is handled differently for focLen and imgCenter
    int Dmax = fullWidth;

    if (fullHeight > fullWidth) {
        Dmax = fullHeight;
    }

    // correct focLens
    if (focLenX < 0) { // they may not be given
        // and 35mm may not be given either
        if (focalLength35mm < 1) {
            focalLength35mm = focalLength * sensorFormatFactor;
        }

        focLenX = focLenY = focalLength / ( 35 * focalLength / focalLength35mm); // focLen must be calculated in pixels
    }

    if (swapXY) {
        x0 = (mirrorX ? 1. - imgYCenter : imgYCenter) * fullWidth;
        y0 = (mirrorY ? 1. - imgXCenter : imgXCenter) * fullHeight;
        fx = focLenY * Dmax;
        fy = focLenX * Dmax;
    } else {
        x0 = (mirrorX ? 1. - imgXCenter : imgXCenter) * fullWidth;
        y0 = (mirrorY ? 1. - imgYCenter : imgYCenter) * fullHeight;
        fx = focLenX * Dmax;
        fy = focLenY * Dmax;
    }

    //printf("FW %i /X0 %g   FH %i /Y0 %g  %g\n",fullWidth,x0,fullHeight,y0, imgYCenter);
}

LCPPersModel::LCPPersModel()
{
    focLen = focDist = aperture = 0;
}

// mode: 0=distortion, 1=vignette, 2=CA
bool LCPPersModel::hasModeData(int mode) const
{
    return (mode == 0 && !vignette.empty() && !vignette.badErr) || (mode == 1 && !base.empty() && !base.badErr)
           || (mode == 2 && !chromRG.empty() && !chromG.empty() && !chromBG.empty() &&
               !chromRG.badErr && !chromG.badErr && !chromBG.badErr);
}

void LCPPersModel::print() const
{
    printf("--- PersModel focLen %g; focDist %g; aperture %g\n", focLen, focDist, aperture);
    printf("Base:\n");
    base.print();

    if (!chromRG.empty()) {
        printf("ChromRG:\n");
        chromRG.print();
    }

    if (!chromG.empty()) {
        printf("ChromG:\n");
        chromG.print();
    }

    if (!chromBG.empty()) {
        printf("ChromBG:\n");
        chromBG.print();
    }

    if (!vignette.empty()) {
        printf("Vignette:\n");
        vignette.print();
    }

    printf("\n");
}

// if !vignette then geometric and CA
LCPMapper::LCPMapper(LCPProfile* pProf, float focalLength, float focalLength35mm, float focusDist, float aperture, bool vignette, bool useCADistP,
                     int fullWidth, int fullHeight, const CoarseTransformParams& coarse, int rawRotationDeg)
{
    if (pProf == NULL) {
        return;
    }

    useCADist = useCADistP;

    // determine in what the image with the RAW landscape in comparison (calibration target)
    // in vignetting, the rotation has not taken place yet
    int rot = 0;

    if (rawRotationDeg >= 0) {
        rot = (coarse.rotate + rawRotationDeg) % 360;
    }

    swapXY  = (rot == 90  || rot == 270);
    bool mirrorX = (rot == 90  || rot == 180);
    bool mirrorY = (rot == 180 || rot == 270);
    //printf("Vign: %i, fullWidth: %i/%i, focLen %g SwapXY: %i / MirX/Y %i / %i on rot:%i from %i\n",vignette, fullWidth, fullHeight, focalLength, swapXY, mirrorX, mirrorY, rot, rawRotationDeg);

    pProf->calcParams(vignette ? 0 : 1, focalLength, focusDist, aperture, &mc, NULL, NULL);
    mc.prepareParams(fullWidth, fullHeight, focalLength, focalLength35mm, pProf->sensorFormatFactor, swapXY, mirrorX, mirrorY);

    if (!vignette) {
        pProf->calcParams(2, focalLength, focusDist, aperture, &chrom[0], &chrom[1], &chrom[2]);

        for (int i = 0; i < 3; i++) {
            chrom[i].prepareParams(fullWidth, fullHeight, focalLength, focalLength35mm, pProf->sensorFormatFactor, swapXY, mirrorX, mirrorY);
        }
    }

    enableCA = !vignette && focusDist > 0;
}

void LCPMapper::correctDistortion(double& x, double& y) const
{
    double xd = (x - mc.x0) / mc.fx, yd = (y - mc.y0) / mc.fy;

    const float* aDist = mc.param;
    double rsqr      = xd * xd + yd * yd;
    double xfac = aDist[swapXY ? 3 : 4], yfac = aDist[swapXY ? 4 : 3];

    double commonFac = (((aDist[2] * rsqr + aDist[1]) * rsqr + aDist[0]) * rsqr + 1.)
                       + 2. * (yfac * yd + xfac * xd);

    double xnew = xd * commonFac + xfac * rsqr;
    double ynew = yd * commonFac + yfac * rsqr;

    x = xnew * mc.fx + mc.x0;
    y = ynew * mc.fy + mc.y0;
}

void LCPMapper::correctCA(double& x, double& y, int channel) const
{
    if (!enableCA) {
        return;
    }

    double rsqr, xgreen, ygreen;

    // First calc the green channel like normal distortion
    // the other are just deviations from it
    double xd = (x - chrom[1].x0) / chrom[1].fx, yd = (y - chrom[1].y0) / chrom[1].fy;

    // Green contains main distortion, just like base
    if (useCADist) {
        const float* aDist = chrom[1].param;
        double rsqr      = xd * xd + yd * yd;
        double xfac = aDist[swapXY ? 3 : 4], yfac = aDist[swapXY ? 4 : 3];

        double commonFac = (((aDist[2] * rsqr + aDist[1]) * rsqr + aDist[0]) * rsqr + 1.)
                           + 2. * (yfac * yd + xfac * xd);

        xgreen = xd * commonFac + aDist[4] * rsqr;
        ygreen = yd * commonFac + aDist[3] * rsqr;
    } else {
        xgreen = xd;
        ygreen = yd;
    }

    if (channel == 1) {
        // green goes directly
        x = xgreen * chrom[1].fx + chrom[1].x0;
        y = ygreen * chrom[1].fy + chrom[1].y0;
    } else {
        // others are diffs from green
        xd = xgreen;
        yd = ygreen;
        rsqr = xd * xd + yd * yd;

        const float* aCA = chrom[channel].param;
        double xfac = aCA[swapXY ? 3 : 4], yfac = aCA[swapXY ? 4 : 3];
        double commonSum = 1. + rsqr * (aCA[0] + rsqr * (aCA[1] + aCA[2] * rsqr)) + 2. * (yfac * yd + xfac * xd);

        x = (chrom[channel].scaleFac * ( xd * commonSum + xfac * rsqr )) * chrom[channel].fx + chrom[channel].x0;
        y = (chrom[channel].scaleFac * ( yd * commonSum + yfac * rsqr )) * chrom[channel].fy + chrom[channel].y0;
    }
}

float LCPMapper::calcVignetteFac(int x, int y) const
{
    // No need for swapXY, since vignette is in RAW and always before rotation
    double xd = ((double)x - mc.x0) / mc.fx, yd = ((double)y - mc.y0) / mc.fy;

    const float* aVig = mc.param;
    double rsqr      = xd * xd + yd * yd;
    double param0Sqr = aVig[0] * aVig[0];

    return 1. + rsqr * (-aVig[0] + rsqr * ((param0Sqr - aVig[1])
                                           - (param0Sqr * aVig[0] - 2.*aVig[0] * aVig[1] + aVig[2]) * rsqr
                                           + (param0Sqr * param0Sqr + aVig[1] * aVig[1]
                                                   + 2.*aVig[0] * aVig[2] - 3.*param0Sqr * aVig[1]) * rsqr * rsqr));
}

LCPProfile::LCPProfile(Glib::ustring fname)
{
    const int BufferSize = 8192;
    char buf[BufferSize];

    XML_Parser parser = XML_ParserCreate(NULL);

    if (!parser) {
        throw "Couldn't allocate memory for XML parser";
    }

    XML_SetElementHandler(parser, XmlStartHandler, XmlEndHandler);
    XML_SetCharacterDataHandler(parser, XmlTextHandler);
    XML_SetUserData(parser, (void *)this);


    isFisheye = inCamProfiles = firstLIDone = inPerspect = inAlternateLensID = inAlternateLensNames = false;
    sensorFormatFactor = 1;

    for (int i = 0; i < MaxPersModelCount; i++) {
        aPersModel[i] = NULL;
    }

    persModelCount = 0;
    *inInvalidTag = 0;

    FILE *pFile = safe_g_fopen(fname, "rb");

    bool done;

    do {
        int bytesRead = (int)fread(buf, 1, BufferSize, pFile);
        done = feof(pFile);

        if (XML_Parse(parser, buf, bytesRead, done) == XML_STATUS_ERROR) {
            throw "Invalid XML in LCP file";
        }
    } while (!done);

    fclose(pFile);

    XML_ParserFree(parser);

    //printf("Parsing %s\n", fname.c_str());
    // Two phase filter: first filter out the very rough ones, that distord the average a lot
    // force it, even if there are few frames (community profiles)
    filterBadFrames(2.0, 0);
    // from the non-distorded, filter again on new average basis, but only if there are enough frames left
    filterBadFrames(1.5, 100);
}

// from all frames not marked as bad already, take average and filter out frames with higher deviation than this if there are enough values
int LCPProfile::filterBadFrames(double maxAvgDevFac, int minFramesLeft)
{
    // take average error per type, then calculated the maximum deviation allowed
    double errBase = 0, errChrom = 0, errVignette = 0;
    int baseCount = 0, chromCount = 0, vignetteCount = 0;

    for (int pm = 0; pm < MaxPersModelCount && aPersModel[pm]; pm++) {
        if (aPersModel[pm]->hasModeData(0)) {
            errVignette += aPersModel[pm]->vignette.meanErr;
            vignetteCount++;
        }

        if (aPersModel[pm]->hasModeData(1)) {
            errBase += aPersModel[pm]->base.meanErr;
            baseCount++;
        }

        if (aPersModel[pm]->hasModeData(2)) {
            errChrom += std::max(std::max(aPersModel[pm]->chromRG.meanErr, aPersModel[pm]->chromG.meanErr), aPersModel[pm]->chromBG.meanErr);
            chromCount++;
        }
    }

    // Only if we have enough frames, filter out errors
    int filtered = 0;

    if (baseCount + chromCount + vignetteCount >= minFramesLeft) {
        if (baseCount > 0) {
            errBase /= (double)baseCount;
        }

        if (chromCount > 0) {
            errChrom /= (double)chromCount;
        }

        if (vignetteCount > 0) {
            errVignette /= (double)vignetteCount;
        }

        // Now mark all the bad ones as bad, and hasModeData will return false;
        for (int pm = 0; pm < MaxPersModelCount && aPersModel[pm]; pm++) {
            if (aPersModel[pm]->hasModeData(0) && aPersModel[pm]->vignette.meanErr > maxAvgDevFac * errVignette) {
                aPersModel[pm]->vignette.badErr = true;
                filtered++;
            }

            if (aPersModel[pm]->hasModeData(1) && aPersModel[pm]->base.meanErr > maxAvgDevFac * errBase) {
                aPersModel[pm]->base.badErr = true;
                filtered++;
            }

            if (aPersModel[pm]->hasModeData(2) &&
                    (aPersModel[pm]->chromRG.meanErr > maxAvgDevFac * errChrom || aPersModel[pm]->chromG.meanErr > maxAvgDevFac * errChrom
                     || aPersModel[pm]->chromBG.meanErr > maxAvgDevFac * errChrom)) {
                aPersModel[pm]->chromRG.badErr = aPersModel[pm]->chromG.badErr = aPersModel[pm]->chromBG.badErr = true;
                filtered++;
            }
        }

        //printf("Filtered %.1f%% frames for maxAvgDevFac %g leaving %i\n", filtered*100./(baseCount+chromCount+vignetteCount), maxAvgDevFac, baseCount+chromCount+vignetteCount-filtered);
    }

    return filtered;
}


// mode: 0=vignette, 1=distortion, 2=CA
void LCPProfile::calcParams(int mode, float focalLength, float focusDist, float aperture, LCPModelCommon *pCorr1, LCPModelCommon *pCorr2, LCPModelCommon *pCorr3) const
{
    float euler = exp(1.0);

    // find the frames with the least distance, focal length wise
    LCPPersModel *pLow = NULL, *pHigh = NULL;

    float focalLengthLog = log(focalLength); //, apertureLog=aperture>0 ? log(aperture) : 0;
    float focusDistLog = focusDist > 0 ? log(focusDist) + euler : 0;

    // Pass 1: determining best focal length, if possible different focusDistances (for the focDist is not given case)
    for (int pm = 0; pm < persModelCount; pm++) {
        float f = aPersModel[pm]->focLen;

        if (aPersModel[pm]->hasModeData(mode)) {
            if (f <= focalLength && (pLow == NULL || f > pLow->focLen || (focusDist == 0 && f == pLow->focLen && pLow->focDist > aPersModel[pm]->focDist))) {
                pLow = aPersModel[pm];
            }

            if (f >= focalLength && (pHigh == NULL || f < pHigh->focLen || (focusDist == 0 && f == pHigh->focLen && pHigh->focDist < aPersModel[pm]->focDist))) {
                pHigh = aPersModel[pm];
            }
        }
    }

    if (!pLow) {
        pLow = pHigh;
    } else if (!pHigh) {
        pHigh = pLow;
    } else {
        // Pass 2: We have some, so take the best aperture for vignette and best focus for CA and distortion
        // there are usually several frame per focal length. In the end pLow will have both flen and apterure/focdis below the target,
        // and vice versa pHigh
        float bestFocLenLow = pLow->focLen, bestFocLenHigh = pHigh->focLen;

        for (int pm = 0; pm < persModelCount; pm++) {
            float aper = aPersModel[pm]->aperture; // float aperLog=log(aper);
            float focDist = aPersModel[pm]->focDist;
            float focDistLog = log(focDist) + euler;
            double meanErr;

            if (aPersModel[pm]->hasModeData(mode)) {
                if (mode == 0) {
                    meanErr = aPersModel[pm]->vignette.meanErr;

                    // by aperture (vignette), and max out focus distance
                    // tests showed doing this by log(aperture) is not as advisable
                    if (aPersModel[pm]->focLen == bestFocLenLow && (
                                (aper == aperture && pLow->vignette.meanErr > meanErr)
                                || (aper >= aperture && aper < pLow->aperture && pLow->aperture > aperture)
                                || (aper <= aperture && (pLow->aperture > aperture || fabs(aperture - aper) < fabs(aperture - pLow->aperture))))) {
                        pLow = aPersModel[pm];
                    }

                    if (aPersModel[pm]->focLen == bestFocLenHigh && (
                                (aper == aperture && pHigh->vignette.meanErr > meanErr)
                                || (aper <= aperture && aper > pHigh->aperture && pHigh->aperture < aperture)
                                || (aper >= aperture && (pHigh->aperture < aperture || fabs(aperture - aper) < fabs(aperture - pHigh->aperture))))) {
                        pHigh = aPersModel[pm];
                    }
                } else {
                    meanErr = (mode == 1 ? aPersModel[pm]->base.meanErr : aPersModel[pm]->chromG.meanErr);

                    if (focusDist > 0) {
                        // by focus distance
                        if (aPersModel[pm]->focLen == bestFocLenLow && (
                                    (focDist == focusDist && (mode == 1 ? pLow->base.meanErr : pLow->chromG.meanErr) > meanErr)
                                    || (focDist >= focusDist && focDist < pLow->focDist && pLow->focDist > focusDist)
                                    || (focDist <= focusDist && (pLow->focDist > focusDist || fabs(focusDistLog - focDistLog) < fabs(focusDistLog - (log(pLow->focDist) + euler)))))) {
                            pLow = aPersModel[pm];
                        }

                        if (aPersModel[pm]->focLen == bestFocLenHigh && (
                                    (focDist == focusDist && (mode == 1 ? pHigh->base.meanErr : pHigh->chromG.meanErr) > meanErr)
                                    || (focDist <= focusDist && focDist > pHigh->focDist && pHigh->focDist < focusDist)
                                    || (focDist >= focusDist && (pHigh->focDist < focusDist || fabs(focusDistLog - focDistLog) < fabs(focusDistLog - (log(pHigh->focDist) + euler)))))) {
                            pHigh = aPersModel[pm];
                        }
                    } else {
                        // no focus distance available, just error
                        if (aPersModel[pm]->focLen == bestFocLenLow && (mode == 1 ? pLow->base.meanErr : pLow->chromG.meanErr) > meanErr) {
                            pLow = aPersModel[pm];
                        }

                        if (aPersModel[pm]->focLen == bestFocLenHigh && (mode == 1 ? pHigh->base.meanErr : pHigh->chromG.meanErr) > meanErr) {
                            pHigh = aPersModel[pm];
                        }
                    }
                }
            }
        }
    }

    if (pLow != NULL && pHigh != NULL) {
        // average out the factors, linear interpolation in logarithmic scale
        float facLow = 0.5;
        bool focLenOnSpot = false; // pretty often, since max/min are often as frames in LCP

        // There is as foclen range, take that as basis
        if (pLow->focLen < pHigh->focLen) {
            facLow = (log(pHigh->focLen) - focalLengthLog) / (log(pHigh->focLen) - log(pLow->focLen));
        } else {
            focLenOnSpot = pLow->focLen == pHigh->focLen && pLow->focLen == focalLength;
        }

        // and average the other factor if available
        if (mode == 0 && pLow->aperture < aperture && pHigh->aperture > aperture) {
            // Mix in aperture
            float facAperLow = (pHigh->aperture - aperture) / (pHigh->aperture - pLow->aperture);
            facLow = focLenOnSpot ? facAperLow : (0.5 * facLow + 0.5 * facAperLow);
        } else if (mode != 0 && focusDist > 0 && pLow->focDist < focusDist && pHigh->focDist > focusDist) {
            // focus distance for all else (if focus distance is given)
            float facDistLow = (log(pHigh->focDist) + euler - focusDistLog) / (log(pHigh->focDist) - log(pLow->focDist));
            facLow = focLenOnSpot ? facDistLow : (0.8 * facLow + 0.2 * facDistLow);
        }

        switch (mode) {
        case 0:  // vignette
            pCorr1->merge(pLow->vignette, pHigh->vignette, facLow);
            break;

        case 1:  // distortion
            pCorr1->merge(pLow->base, pHigh->base, facLow);
            break;

        case 2:  // CA
            pCorr1->merge(pLow->chromRG, pHigh->chromRG, facLow);
            pCorr2->merge(pLow->chromG,  pHigh->chromG,  facLow);
            pCorr3->merge(pLow->chromBG, pHigh->chromBG, facLow);
            break;
        }

        //printf("LCP mode=%i, dist: %g found frames: Fno %g-%g; FocLen %g-%g; Dist %g-%g with weight %g\n", mode, focusDist, pLow->aperture, pHigh->aperture, pLow->focLen, pHigh->focLen, pLow->focDist, pHigh->focDist, facLow);
    } else {
        printf("Error: LCP file contained no %s parameters\n", mode == 0 ? "vignette" : mode == 1 ? "distortion" : "CA" );
    }
}

void LCPProfile::print() const
{
    printf("=== Profile %s\n", profileName.c_str());
    printf("Frames: %i, RAW: %i; Fisheye: %i; Sensorformat: %f\n", persModelCount, isRaw, isFisheye, sensorFormatFactor);

    for (int pm = 0; pm < persModelCount; pm++) {
        aPersModel[pm]->print();
    }
}

void XMLCALL LCPProfile::XmlStartHandler(void *pLCPProfile, const char *el, const char **attr)
{
    LCPProfile *pProf = static_cast<LCPProfile*>(pLCPProfile);
    bool parseAttr = false;

    if (*pProf->inInvalidTag) {
        return;    // We ignore everything in dirty tag till it's gone
    }

    // clean up tagname
    const char* src = strrchr(el, ':');

    if (src == NULL) {
        src = const_cast<char*>(el);
    } else {
        src++;
    }

    strcpy(pProf->lastTag, src);

    if (!strcmp("VignetteModelPiecewiseParam", src)) {
        strcpy(pProf->inInvalidTag, src);
    }

    if (!strcmp("CameraProfiles", src)) {
        pProf->inCamProfiles = true;
    }

    if (!strcmp("AlternateLensIDs", src)) {
        pProf->inAlternateLensID = true;
    }

    if (!strcmp("AlternateLensNames", src)) {
        pProf->inAlternateLensNames = true;
    }

    if (!pProf->inCamProfiles || pProf->inAlternateLensID || pProf->inAlternateLensNames) {
        return;
    }

    if (!strcmp("li", src)) {
        pProf->pCurPersModel = new LCPPersModel();
        pProf->pCurCommon = &pProf->pCurPersModel->base; // iterated to next tags within persModel
        return;
    }

    if (!strcmp("PerspectiveModel", src)) {
        pProf->firstLIDone = true;
        pProf->inPerspect = true;
        return;
    } else if (!strcmp("FisheyeModel", src)) {
        pProf->firstLIDone = true;
        pProf->inPerspect = true;
        pProf->isFisheye = true; // just misses third param, and different path, rest is the same
        return;
    } else if (!strcmp("Description", src)) {
        parseAttr = true;
    }

    // Move pointer to general section
    if (pProf->inPerspect) {
        if (!strcmp("ChromaticRedGreenModel", src)) {
            pProf->pCurCommon = &pProf->pCurPersModel->chromRG;
            parseAttr = true;
        } else if (!strcmp("ChromaticGreenModel", src)) {
            pProf->pCurCommon = &pProf->pCurPersModel->chromG;
            parseAttr = true;
        } else if (!strcmp("ChromaticBlueGreenModel", src)) {
            pProf->pCurCommon = &pProf->pCurPersModel->chromBG;
            parseAttr = true;
        } else if (!strcmp("VignetteModel", src)) {
            pProf->pCurCommon = &pProf->pCurPersModel->vignette;
            parseAttr = true;
        }
    }

    // some profiles (espc. Pentax) have a different structure that is attributes based
    // simulate tags by feeding them in
    if (parseAttr && attr != NULL) {
        for (int i = 0; attr[i]; i += 2) {
            const char* nameStart = strrchr(attr[i], ':');

            if (nameStart == NULL) {
                nameStart = const_cast<char*>(attr[i]);
            } else {
                nameStart++;
            }

            strcpy(pProf->lastTag, nameStart);
            XmlTextHandler(pLCPProfile, attr[i + 1], strlen(attr[i + 1]));
        }
    }
}

void XMLCALL LCPProfile::XmlTextHandler(void *pLCPProfile, const XML_Char *s, int len)
{
    LCPProfile *pProf = static_cast<LCPProfile*>(pLCPProfile);

    if (!pProf->inCamProfiles || pProf->inAlternateLensID || pProf->inAlternateLensNames || *pProf->inInvalidTag) {
        return;
    }

    // Check if it contains non-whitespaces (there are several calls to this for one tag unfortunately)
    bool onlyWhiteSpace = true;
    int i = 0;

    while (i < len && onlyWhiteSpace) {
        onlyWhiteSpace = isspace(s[i]);
        i++;
    }

    if (onlyWhiteSpace) {
        return;
    }

    // convert to null terminated
    char raw[len + 1];
    memcpy(raw, s, len);
    raw[len] = 0;
    char* tag = pProf->lastTag;

    // Common data section
    if (!pProf->firstLIDone) {
        // Generic tags are the same for all
        if (!strcmp("ProfileName", tag)) {
            pProf->profileName = Glib::ustring(raw);
        } else if (!strcmp("Model", tag)) {
            pProf->camera = Glib::ustring(raw);
        } else if (!strcmp("Lens", tag)) {
            pProf->lens = Glib::ustring(raw);
        } else if (!strcmp("CameraPrettyName", tag)) {
            pProf->cameraPrettyName = Glib::ustring(raw);
        } else if (!strcmp("LensPrettyName", tag)) {
            pProf->lensPrettyName = Glib::ustring(raw);
        } else if (!strcmp("CameraRawProfile", tag)) {
            pProf->isRaw = !strcmp("True", raw);
        }
    }

    // --- Now all floating points. Must replace local dot characters
    // WARNING: called by different threads, that may run on different local settings,
    // so don't use system params
    if (atof("1,2345") == 1.2345) {
        char* p = raw;

        while (*p) {
            if (*p == '.') {
                *p = ',';
            }

            p++;
        }
    }

    if (!pProf->firstLIDone) {
        if (!strcmp("SensorFormatFactor", tag)) {
            pProf->sensorFormatFactor = atof(raw);
        }
    }

    // Perspective model base data
    if (!strcmp("FocalLength", tag)) {
        pProf->pCurPersModel->focLen = atof(raw);
    } else if (!strcmp("FocusDistance", tag)) {
        double focDist = atof(raw);
        pProf->pCurPersModel->focDist = focDist < 10000 ? focDist : 10000;
    } else if (!strcmp("ApertureValue", tag)) {
        pProf->pCurPersModel->aperture = atof(raw);
    }

    // Section depended
    if (!strcmp("FocalLengthX", tag)) {
        pProf->pCurCommon->focLenX = atof(raw);
    } else if (!strcmp("FocalLengthY", tag)) {
        pProf->pCurCommon->focLenY = atof(raw);
    } else if (!strcmp("ImageXCenter", tag)) {
        pProf->pCurCommon->imgXCenter = atof(raw);
    } else if (!strcmp("ImageYCenter", tag)) {
        pProf->pCurCommon->imgYCenter = atof(raw);
    } else if (!strcmp("ScaleFactor", tag)) {
        pProf->pCurCommon->scaleFac = atof(raw);
    } else if (!strcmp("ResidualMeanError", tag)) {
        pProf->pCurCommon->meanErr = atof(raw);
    } else if (!strcmp("RadialDistortParam1", tag) || !strcmp("VignetteModelParam1", tag)) {
        pProf->pCurCommon->param[0] = atof(raw);
    } else if (!strcmp("RadialDistortParam2", tag) || !strcmp("VignetteModelParam2", tag)) {
        pProf->pCurCommon->param[1] = atof(raw);
    } else if (!strcmp("RadialDistortParam3", tag) || !strcmp("VignetteModelParam3", tag)) {
        pProf->pCurCommon->param[2] = atof(raw);
    } else if (!strcmp("RadialDistortParam4", tag) || !strcmp("TangentialDistortParam1", tag)) {
        pProf->pCurCommon->param[3] = atof(raw);
    } else if (!strcmp("RadialDistortParam5", tag) || !strcmp("TangentialDistortParam2", tag)) {
        pProf->pCurCommon->param[4] = atof(raw);
    }
}

void XMLCALL LCPProfile::XmlEndHandler(void *pLCPProfile, const char *el)
{
    LCPProfile *pProf = static_cast<LCPProfile*>(pLCPProfile);

    // We ignore everything in dirty tag till it's gone
    if (*pProf->inInvalidTag) {
        if (strstr(el, pProf->inInvalidTag)) {
            *pProf->inInvalidTag = 0;
        }

        return;
    }

    if (strstr(el, ":CameraProfiles")) {
        pProf->inCamProfiles = false;
    }

    if (strstr(el, ":AlternateLensIDs")) {
        pProf->inAlternateLensID = false;
    }

    if (strstr(el, ":AlternateLensNames")) {
        pProf->inAlternateLensNames = false;
    }

    if (!pProf->inCamProfiles || pProf->inAlternateLensID || pProf->inAlternateLensNames) {
        return;
    }

    if (strstr(el, ":PerspectiveModel") || strstr(el, ":FisheyeModel")) {
        pProf->inPerspect = false;
    } else if (strstr(el, ":li")) {
        pProf->aPersModel[pProf->persModelCount] = pProf->pCurPersModel;
        pProf->pCurPersModel = NULL;
        pProf->persModelCount++;
    }
}

// Generates as singleton
LCPStore* LCPStore::getInstance()
{
    static LCPStore* instance_ = 0;

    if ( instance_ == 0 ) {
        static MyMutex smutex_;
        MyMutex::MyLock lock(smutex_);

        if ( instance_ == 0 ) {
            instance_ = new LCPStore();
        }
    }

    return instance_;
}

LCPProfile* LCPStore::getProfile (Glib::ustring filename)
{
    if (filename.length() == 0 || !isValidLCPFileName(filename)) {
        return NULL;
    }

    MyMutex::MyLock lock(mtx);

    std::map<Glib::ustring, LCPProfile*>::iterator r = profileCache.find (filename);

    if (r != profileCache.end()) {
        return r->second;
    }

    // Add profile (if exists)
    profileCache[filename] = new LCPProfile(filename);
    //profileCache[filename]->print();
    return profileCache[filename];
}

bool LCPStore::isValidLCPFileName(Glib::ustring filename) const
{
    if (!safe_file_test (filename, Glib::FILE_TEST_EXISTS) || safe_file_test (filename, Glib::FILE_TEST_IS_DIR)) {
        return false;
    }

    size_t pos = filename.find_last_of ('.');
    return pos > 0 && !filename.casefold().compare (pos, 4, ".lcp");
}

Glib::ustring LCPStore::getDefaultCommonDirectory() const
{
    Glib::ustring dir;

#ifdef WIN32
    WCHAR pathW[MAX_PATH] = {0};
    char pathA[MAX_PATH];

    if (SHGetSpecialFolderPathW(NULL, pathW, CSIDL_COMMON_APPDATA, false)) {
        char pathA[MAX_PATH];
        WideCharToMultiByte(CP_UTF8, 0, pathW, -1, pathA, MAX_PATH, 0, 0);
        Glib::ustring fullDir = Glib::ustring(pathA) + Glib::ustring("\\Adobe\\CameraRaw\\LensProfiles\\1.0");

        if (safe_file_test(fullDir, Glib::FILE_TEST_IS_DIR)) {
            dir = fullDir;
        }
    }

#endif

    // TODO: Add Mac paths here

    return dir;
}
