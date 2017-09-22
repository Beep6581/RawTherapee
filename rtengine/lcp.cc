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

#include <algorithm>
#include <cstring>

#include <glib/gstdio.h>

#ifdef WIN32
#include <shlobj.h>
#include <windows.h>
#endif

#include "lcp.h"

#include "settings.h"

namespace rtengine
{

extern const Settings* settings;

}

class rtengine::LCPProfile::LCPPersModel
{
public:
    LCPPersModel();
    bool hasModeData(LCPCorrectionMode mode) const;
    void print() const;

    float focLen;
    float focDist;
    float aperture;  // this is what it refers to

    LCPModelCommon base;  // base perspective correction
    LCPModelCommon chromRG;
    LCPModelCommon chromG;
    LCPModelCommon chromBG;  // red/green, green, blue/green (may be empty)
    LCPModelCommon vignette;  // vignette (may be empty)
};

rtengine::LCPModelCommon::LCPModelCommon() :
    foc_len_x(-1.0f),
    foc_len_y(-1.0f),
    img_center_x(0.5f),
    img_center_y(0.5f),
    param{{}},
    scale_factor(1.0f),
    mean_error(0.0),
    bad_error(false),
    x0(0.0f),
    y0(0.0f),
    fx(0.0f),
    fy(0.0f),
    rfx(0.0f),
    rfy(0.0f),
    vign_param{{}}
{
}

bool rtengine::LCPModelCommon::empty() const
{
    return
        param[0] == 0.0f
        && param[1] == 0.0f
        && param[2] == 0.0f;
}

void rtengine::LCPModelCommon::print() const
{
    std::printf("focLen %g/%g; imgCenter %g/%g; scale %g; err %g\n", foc_len_x, foc_len_y, img_center_x, img_center_y, scale_factor, mean_error);
    std::printf("xy0 %g/%g  fxy %g/%g\n", x0, y0, fx, fy);
    std::printf("param: %g/%g/%g/%g/%g\n", param[0], param[1], param[2], param[3], param[4]);
}

// weighted merge two parameters
void rtengine::LCPModelCommon::merge(const LCPModelCommon& a, const LCPModelCommon& b, float facA)
{
    const float facB = 1.0f - facA;

    foc_len_x    = facA * a.foc_len_x    + facB * b.foc_len_x;
    foc_len_y    = facA * a.foc_len_y    + facB * b.foc_len_y;
    img_center_x = facA * a.img_center_x + facB * b.img_center_x;
    img_center_y = facA * a.img_center_y + facB * b.img_center_y;
    scale_factor = facA * a.scale_factor + facB * b.scale_factor;
    mean_error   = facA * a.mean_error   + facB * b.mean_error;

    for (int i = 0; i < 5; ++i) {
        param[i] = facA * a.param[i] + facB * b.param[i];
    }

    const float param0Sqr = param[0] * param[0];

    vign_param[0] = -param[0];
    vign_param[1] = param0Sqr - param[1];
    vign_param[2] = param0Sqr * param[0] - 2.0f * param[0] * param[1] + param[2];
    vign_param[3] = param0Sqr * param0Sqr + param[1] * param[1] + 2.0f * param[0] * param[2] - 3.0f * param0Sqr * param[1];

}

void rtengine::LCPModelCommon::prepareParams(
    int fullWidth,
    int fullHeight,
    float focalLength,
    float focalLength35mm,
    float sensorFormatFactor,
    bool swapXY,
    bool mirrorX,
    bool mirrorY
)
{
    // Mention that the Adobe technical paper has a bug here, the DMAX is handled differently for focLen and imgCenter
    const int Dmax = std::max(fullWidth, fullHeight);

    // correct focLens
    if (foc_len_x < 0.0f) { // they may not be given
        // and 35mm may not be given either
        if (focalLength35mm < 1.0f) {
            focalLength35mm = focalLength * sensorFormatFactor;
        }

        foc_len_x = foc_len_y = focalLength / (35.0f * focalLength / focalLength35mm); // focLen must be calculated in pixels
    }

    if (swapXY) {
        x0 = (mirrorX ? 1.0f - img_center_y : img_center_y) * fullWidth;
        y0 = (mirrorY ? 1.0f - img_center_x : img_center_x) * fullHeight;
        fx = foc_len_y * Dmax;
        fy = foc_len_x * Dmax;
    } else {
        x0 = (mirrorX ? 1.0f - img_center_x : img_center_x) * fullWidth;
        y0 = (mirrorY ? 1.0f - img_center_y : img_center_y) * fullHeight;
        fx = foc_len_x * Dmax;
        fy = foc_len_y * Dmax;
    }
    rfx = 1.0f / fx;
    rfy = 1.0f / fy;

    //std::printf("FW %i /X0 %g   FH %i /Y0 %g  %g\n",fullWidth,x0,fullHeight,y0, imgYCenter);
}

rtengine::LCPProfile::LCPPersModel::LCPPersModel() :
    focLen(0.f),
    focDist(0.f),
    aperture(0.f)
{
}

bool rtengine::LCPProfile::LCPPersModel::hasModeData(LCPCorrectionMode mode) const
{
    switch (mode) {
        case LCPCorrectionMode::VIGNETTE: {
            return !vignette.empty() && !vignette.bad_error;
        }

        case LCPCorrectionMode::DISTORTION: {
            return !base.empty() && !base.bad_error;
        }

        case LCPCorrectionMode::CA: {
            return
                !chromRG.empty()
                && !chromG.empty()
                && !chromBG.empty()
                && !chromRG.bad_error
                && !chromG.bad_error
                && !chromBG.bad_error;
        }
    }

    assert(false);
    return false;
}

void rtengine::LCPProfile::LCPPersModel::print() const
{
    std::printf("--- PersModel focLen %g; focDist %g; aperture %g\n", focLen, focDist, aperture);
    std::printf("Base:\n");
    base.print();

    if (!chromRG.empty()) {
        std::printf("ChromRG:\n");
        chromRG.print();
    }

    if (!chromG.empty()) {
        std::printf("ChromG:\n");
        chromG.print();
    }

    if (!chromBG.empty()) {
        std::printf("ChromBG:\n");
        chromBG.print();
    }

    if (!vignette.empty()) {
        std::printf("Vignette:\n");
        vignette.print();
    }

    std::printf("\n");
}

rtengine::LCPProfile::LCPProfile(const Glib::ustring& fname) :
    isFisheye(false),
    sensorFormatFactor(1.f),
    persModelCount(0),
    inCamProfiles(false),
    firstLIDone(false),
    inPerspect(false),
    inAlternateLensID(false),
    inAlternateLensNames(false),
    lastTag{},
    inInvalidTag{},
    pCurPersModel(nullptr),
    pCurCommon(nullptr),
    aPersModel{}
{
    const int BufferSize = 8192;
    char buf[BufferSize];

    XML_Parser parser = XML_ParserCreate(nullptr);

    if (!parser) {
        throw "Couldn't allocate memory for XML parser";
    }

    XML_SetElementHandler(parser, XmlStartHandler, XmlEndHandler);
    XML_SetCharacterDataHandler(parser, XmlTextHandler);
    XML_SetUserData(parser, static_cast<void*>(this));

    FILE* const pFile = g_fopen(fname.c_str (), "rb");

    if (pFile) {
        bool done;

        do {
            int bytesRead = fread(buf, 1, BufferSize, pFile);
            done = feof(pFile);

            if (XML_Parse(parser, buf, bytesRead, done) == XML_STATUS_ERROR) {
                XML_ParserFree(parser);
                throw "Invalid XML in LCP file";
            }
        } while (!done);

        fclose(pFile);
    }

    XML_ParserFree(parser);

    if (settings->verbose) {
        std::printf("Parsing %s\n", fname.c_str());
    }
    // Two phase filter: first filter out the very rough ones, that distord the average a lot
    // force it, even if there are few frames (community profiles)
    filterBadFrames(LCPCorrectionMode::VIGNETTE, 2.0, 0);
    filterBadFrames(LCPCorrectionMode::CA, 2.0, 0);
    // from the non-distorded, filter again on new average basis, but only if there are enough frames left
    filterBadFrames(LCPCorrectionMode::VIGNETTE, 1.5, 50);
    filterBadFrames(LCPCorrectionMode::CA, 1.5, 50);
}

rtengine::LCPProfile::~LCPProfile()
{
    delete pCurPersModel;

    for (int i = 0; i < MaxPersModelCount; ++i) {
        delete aPersModel[i];
    }
}

void rtengine::LCPProfile::calcParams(
    LCPCorrectionMode mode,
    float focalLength,
    float focusDist,
    float aperture,
    LCPModelCommon* pCorr1,
    LCPModelCommon* pCorr2,
    LCPModelCommon* pCorr3
) const
{
    const float euler = std::exp(1.0);

    // find the frames with the least distance, focal length wise
    LCPPersModel* pLow = nullptr;
    LCPPersModel* pHigh = nullptr;

    const float focalLengthLog = std::log(focalLength); //, apertureLog=aperture>0 ? std::log(aperture) : 0;
    const float focusDistLog = focusDist > 0 ? std::log(focusDist) + euler : 0;

    // Pass 1: determining best focal length, if possible different focusDistances (for the focDist is not given case)
    for (int pm = 0; pm < persModelCount; ++pm) {
        const float f = aPersModel[pm]->focLen;

        if (aPersModel[pm]->hasModeData(mode)) {
            if (
                f <= focalLength
                && (
                    pLow == nullptr
                    || f > pLow->focLen
                    || (
                        focusDist == 0
                        && f == pLow->focLen
                        && pLow->focDist > aPersModel[pm]->focDist
                    )
                )
            ) {
                pLow = aPersModel[pm];
            }

            if (
                f >= focalLength
                && (
                    pHigh == nullptr
                    || f < pHigh->focLen
                    || (
                        focusDist == 0
                        && f == pHigh->focLen
                        && pHigh->focDist < aPersModel[pm]->focDist
                    )
                )
            ) {
                pHigh = aPersModel[pm];
            }
        }
    }

    if (!pLow) {
        pLow = pHigh;
    }
    else if (!pHigh) {
        pHigh = pLow;
    }
    else {
        // Pass 2: We have some, so take the best aperture for vignette and best focus for CA and distortion
        // there are usually several frame per focal length. In the end pLow will have both flen and apterure/focdis below the target,
        // and vice versa pHigh
        const float bestFocLenLow = pLow->focLen;
        const float bestFocLenHigh = pHigh->focLen;

        for (int pm = 0; pm < persModelCount; ++pm) {
            const float aper = aPersModel[pm]->aperture; // float aperLog=std::log(aper);
            const float focDist = aPersModel[pm]->focDist;
            const float focDistLog = std::log(focDist) + euler;

            double meanErr;

            if (aPersModel[pm]->hasModeData(mode)) {
                double lowMeanErr = 0.0;
                double highMeanErr = 0.0;

                switch (mode) {
                    case LCPCorrectionMode::VIGNETTE: {
                        meanErr = aPersModel[pm]->vignette.mean_error;
                        lowMeanErr = pLow->vignette.mean_error;
                        highMeanErr = pHigh->vignette.mean_error;
                        break;
                    }

                    case LCPCorrectionMode::DISTORTION: {
                        meanErr = aPersModel[pm]->base.mean_error;
                        lowMeanErr = pLow->base.mean_error;
                        highMeanErr = pHigh->base.mean_error;
                        break;
                    }

                    case LCPCorrectionMode::CA: {
                        meanErr = aPersModel[pm]->chromG.mean_error;
                        lowMeanErr = pLow->chromG.mean_error;
                        highMeanErr = pHigh->chromG.mean_error;
                        break;
                    }
                }

                if (aperture > 0 && mode != LCPCorrectionMode::CA) {
                    if (
                        aPersModel[pm]->focLen == bestFocLenLow
                        && (
                            (
                                aper == aperture
                                && lowMeanErr > meanErr
                            )
                            || (
                                aper >= aperture
                                && aper < pLow->aperture
                                && pLow->aperture > aperture
                            )
                            || (
                                aper <= aperture
                                && (
                                    pLow->aperture > aperture
                                    || fabs(aperture - aper) < fabs(aperture - pLow->aperture)
                                )
                            )
                        )
                    ) {
                        pLow = aPersModel[pm];
                    }

                    if (
                        aPersModel[pm]->focLen == bestFocLenHigh
                        && (
                            (
                                aper == aperture
                                && highMeanErr > meanErr
                            )
                            || (
                                aper <= aperture
                                && aper > pHigh->aperture
                                && pHigh->aperture < aperture
                            )
                            || (
                                aper >= aperture
                                && (
                                    pHigh->aperture < aperture
                                    || fabs(aperture - aper) < fabs(aperture - pHigh->aperture)
                                )
                            )
                        )
                    ) {
                        pHigh = aPersModel[pm];
                    }
                }
                else if (focusDist > 0 && mode != LCPCorrectionMode::VIGNETTE) {
                    // by focus distance
                    if (
                        aPersModel[pm]->focLen == bestFocLenLow
                        && (
                            (
                                focDist == focusDist
                                && lowMeanErr > meanErr
                            )
                            || (
                                focDist >= focusDist
                                && focDist < pLow->focDist
                                && pLow->focDist > focusDist
                            )
                            || (
                                focDist <= focusDist
                                && (
                                    pLow->focDist > focusDist
                                    || fabs(focusDistLog - focDistLog) < fabs(focusDistLog - (std::log(pLow->focDist) + euler))
                                )
                            )
                        )
                    ) {
                        pLow = aPersModel[pm];
                    }

                    if (
                        aPersModel[pm]->focLen == bestFocLenHigh
                        && (
                            (
                                focDist == focusDist
                                && highMeanErr > meanErr
                            )
                            || (
                                focDist <= focusDist
                                && focDist > pHigh->focDist
                                && pHigh->focDist < focusDist
                            )
                            || (
                                focDist >= focusDist
                                && (
                                    pHigh->focDist < focusDist
                                    || fabs(focusDistLog - focDistLog) < fabs(focusDistLog - (std::log(pHigh->focDist) + euler))
                                )
                            )
                        )
                    ) {
                        pHigh = aPersModel[pm];
                    }
                }
                else {
                    // no focus distance available, just error
                    if (aPersModel[pm]->focLen == bestFocLenLow && lowMeanErr > meanErr) {
                        pLow = aPersModel[pm];
                    }

                    if (aPersModel[pm]->focLen == bestFocLenHigh && highMeanErr > meanErr) {
                        pHigh = aPersModel[pm];
                    }
                }

            }
        }
    }

    if (pLow != nullptr && pHigh != nullptr) {
        // average out the factors, linear interpolation in logarithmic scale
        float facLow = 0.5f;
        bool focLenOnSpot = false; // pretty often, since max/min are often as frames in LCP

        // There is as foclen range, take that as basis
        if (pLow->focLen < pHigh->focLen) {
            facLow = (std::log(pHigh->focLen) - focalLengthLog) / (std::log(pHigh->focLen) - std::log(pLow->focLen));
        } else {
            focLenOnSpot = pLow->focLen == pHigh->focLen && pLow->focLen == focalLength;
        }

        // and average the other factor if available
        if (
            mode == LCPCorrectionMode::VIGNETTE
            && pLow->aperture < aperture
            && pHigh->aperture > aperture
        ) {
            // Mix in aperture
            const float facAperLow = (pHigh->aperture - aperture) / (pHigh->aperture - pLow->aperture);
            facLow = focLenOnSpot ? facAperLow : (0.5 * facLow + 0.5 * facAperLow);
        }
        else if (
            mode != LCPCorrectionMode::VIGNETTE
            && focusDist > 0
            && pLow->focDist < focusDist
            && pHigh->focDist > focusDist
        ) {
            // focus distance for all else (if focus distance is given)
            const float facDistLow = (std::log(pHigh->focDist) + euler - focusDistLog) / (std::log(pHigh->focDist) - std::log(pLow->focDist));
            facLow = focLenOnSpot ? facDistLow : (0.8 * facLow + 0.2 * facDistLow);
        }

        switch (mode) {
            case LCPCorrectionMode::VIGNETTE: {
                pCorr1->merge(pLow->vignette, pHigh->vignette, facLow);
                break;
            }

            case LCPCorrectionMode::DISTORTION: {
                pCorr1->merge(pLow->base, pHigh->base, facLow);
                break;
            }

            case LCPCorrectionMode::CA: {
                pCorr1->merge(pLow->chromRG, pHigh->chromRG, facLow);
                pCorr2->merge(pLow->chromG,  pHigh->chromG,  facLow);
                pCorr3->merge(pLow->chromBG, pHigh->chromBG, facLow);
                break;
            }
        }

        if (settings->verbose) {
            std::printf("LCP mode=%i, dist: %g found frames: Fno %g-%g; FocLen %g-%g; Dist %g-%g with weight %g\n", toUnderlying(mode), focusDist, pLow->aperture, pHigh->aperture, pLow->focLen, pHigh->focLen, pLow->focDist, pHigh->focDist, facLow);
        }
    } else {
        if (settings->verbose) {
            std::printf("Error: LCP file contained no %s parameters\n", mode == LCPCorrectionMode::VIGNETTE ? "vignette" : mode == LCPCorrectionMode::DISTORTION ? "distortion" : "CA" );
        }
    }
}

void rtengine::LCPProfile::print() const
{
    std::printf("=== Profile %s\n", profileName.c_str());
    std::printf("Frames: %i, RAW: %i; Fisheye: %i; Sensorformat: %f\n", persModelCount, isRaw, isFisheye, sensorFormatFactor);

    for (int pm = 0; pm < persModelCount; ++pm) {
        aPersModel[pm]->print();
    }
}

// from all frames not marked as bad already, take average and filter out frames with higher deviation than this if there are enough values
int rtengine::LCPProfile::filterBadFrames(LCPCorrectionMode mode, double maxAvgDevFac, int minFramesLeft)
{
    // take average error, then calculated the maximum deviation allowed
    double err = 0.0;
    int count = 0;

    for (int pm = 0; pm < MaxPersModelCount && aPersModel[pm]; ++pm) {
        if (aPersModel[pm]->hasModeData(mode)) {
            ++count;
            switch (mode) {
                case LCPCorrectionMode::VIGNETTE: {
                    err += aPersModel[pm]->vignette.mean_error;
                    break;
                }

                case LCPCorrectionMode::DISTORTION: {
                    err += aPersModel[pm]->base.mean_error;
                    break;
                }

                case LCPCorrectionMode::CA: {
                    err += rtengine::max(aPersModel[pm]->chromRG.mean_error, aPersModel[pm]->chromG.mean_error, aPersModel[pm]->chromBG.mean_error);
                    break;
                }
            }
        }
    }

    // Only if we have enough frames, filter out errors
    int filtered = 0;

    if (count >= minFramesLeft) {
        if (count > 0) {
            err /= count;
        }

        // Now mark all the bad ones as bad, and hasModeData will return false;
        for (int pm = 0; pm < MaxPersModelCount && aPersModel[pm]; ++pm) {
            if (aPersModel[pm]->hasModeData(mode)) {
                switch (mode) {
                    case LCPCorrectionMode::VIGNETTE: {
                        if (aPersModel[pm]->vignette.mean_error > maxAvgDevFac * err) {
                            aPersModel[pm]->vignette.bad_error = true;
                            filtered++;
                        }
                        break;
                    }

                    case LCPCorrectionMode::DISTORTION: {
                        if (aPersModel[pm]->base.mean_error > maxAvgDevFac * err) {
                            aPersModel[pm]->base.bad_error = true;
                            filtered++;
                        }
                        break;
                    }

                    case LCPCorrectionMode::CA: {
                        if (
                            aPersModel[pm]->chromRG.mean_error > maxAvgDevFac * err
                            || aPersModel[pm]->chromG.mean_error > maxAvgDevFac * err
                            || aPersModel[pm]->chromBG.mean_error > maxAvgDevFac * err
                        ) {
                            aPersModel[pm]->chromRG.bad_error = true;
                            aPersModel[pm]->chromG.bad_error = true;
                            aPersModel[pm]->chromBG.bad_error = true;
                            ++filtered;
                        }
                        break;
                    }
                }
            }
        }

        if (settings->verbose && count) {
            std::printf("Filtered %.1f%% frames for maxAvgDevFac %g leaving %i\n", filtered * 100.f / count, maxAvgDevFac, count - filtered);
        }
    }

    return filtered;
}

void rtengine::LCPProfile::handle_text(const std::string& text)
{
    // Check if it contains non-whitespaces (there are several calls to this for one tag unfortunately)
    bool onlyWhiteSpace = true;
    for (auto c : text) {
        if (!std::isspace(c)) {
            onlyWhiteSpace = false;
            break;
        }
    }

    if (onlyWhiteSpace) {
        return;
    }

    LCPProfile* const pProf = this;

    // convert to null terminated
    const std::string tag = pProf->lastTag;

    // Common data section
    if (!pProf->firstLIDone) {
        // Generic tags are the same for all
        if (tag == "ProfileName") {
            pProf->profileName = text;
        } else if (tag == "Model") {
            pProf->camera = text;
        } else if (tag == "Lens") {
            pProf->lens = text;
        } else if (tag == "CameraPrettyName") {
            pProf->cameraPrettyName = text;
        } else if (tag == "LensPrettyName") {
            pProf->lensPrettyName = text;
        } else if (tag == "CameraRawProfile") {
            pProf->isRaw = text == "True";
        }
    }

    // Locale should be already set
    assert(std::atof("1.2345") == 1.2345);

    if (!pProf->firstLIDone) {
        if (tag == "SensorFormatFactor") {
            pProf->sensorFormatFactor = std::atof(text.c_str());
        }
    }

    // Perspective model base data
    if (tag == "FocalLength") {
        pProf->pCurPersModel->focLen = std::atof(text.c_str());
    } else if (tag == "FocusDistance") {
        double focDist = std::atof(text.c_str());
        pProf->pCurPersModel->focDist = focDist < 10000 ? focDist : 10000;
    } else if (tag == "ApertureValue") {
        pProf->pCurPersModel->aperture = std::atof(text.c_str());
    }

    // Section depended
    if (tag == "FocalLengthX") {
        pProf->pCurCommon->foc_len_x = std::atof(text.c_str());
    } else if (tag == "FocalLengthY") {
        pProf->pCurCommon->foc_len_y = std::atof(text.c_str());
    } else if (tag == "ImageXCenter") {
        pProf->pCurCommon->img_center_x = std::atof(text.c_str());
    } else if (tag == "ImageYCenter") {
        pProf->pCurCommon->img_center_y = std::atof(text.c_str());
    } else if (tag == "ScaleFactor") {
        pProf->pCurCommon->scale_factor = std::atof(text.c_str());
    } else if (tag == "ResidualMeanError") {
        pProf->pCurCommon->mean_error = std::atof(text.c_str());
    } else if (tag == "RadialDistortParam1" || tag == "VignetteModelParam1") {
        pProf->pCurCommon->param[0] = std::atof(text.c_str());
    } else if (tag == "RadialDistortParam2" || tag == "VignetteModelParam2") {
        pProf->pCurCommon->param[1] = std::atof(text.c_str());
    } else if (tag == "RadialDistortParam3" || tag == "VignetteModelParam3") {
        pProf->pCurCommon->param[2] = std::atof(text.c_str());
    } else if (tag == "RadialDistortParam4" || tag == "TangentialDistortParam1") {
        pProf->pCurCommon->param[3] = std::atof(text.c_str());
    } else if (tag == "RadialDistortParam5" || tag == "TangentialDistortParam2") {
        pProf->pCurCommon->param[4] = std::atof(text.c_str());
    }
}

void XMLCALL rtengine::LCPProfile::XmlStartHandler(void* pLCPProfile, const char* el, const char** attr)
{
    LCPProfile* const pProf = static_cast<LCPProfile*>(pLCPProfile);

    bool parseAttr = false;

    if (*pProf->inInvalidTag) {
        return;    // We ignore everything in dirty tag till it's gone
    }

    // clean up tagname
    const char* src = strrchr(el, ':');

    if (src == nullptr) {
        src = el;
    } else {
        ++src;
    }

    strcpy(pProf->lastTag, src);

    const std::string src_str = src;

    if (src_str == "VignetteModelPiecewiseParam") {
        strcpy(pProf->inInvalidTag, src);
    }

    if (src_str == "CameraProfiles") {
        pProf->inCamProfiles = true;
    }

    if (src_str == "AlternateLensIDs") {
        pProf->inAlternateLensID = true;
    }

    if (src_str == "AlternateLensNames") {
        pProf->inAlternateLensNames = true;
    }

    if (
        !pProf->inCamProfiles
        || pProf->inAlternateLensID
        || pProf->inAlternateLensNames
    ) {
        return;
    }

    if (src_str == "li") {
        pProf->pCurPersModel = new LCPPersModel();
        pProf->pCurCommon = &pProf->pCurPersModel->base; // iterated to next tags within persModel
        return;
    }

    if (src_str == "PerspectiveModel") {
        pProf->firstLIDone = true;
        pProf->inPerspect = true;
        return;
    } else if (src_str == "FisheyeModel") {
        pProf->firstLIDone = true;
        pProf->inPerspect = true;
        pProf->isFisheye = true; // just misses third param, and different path, rest is the same
        return;
    } else if (src_str == "Description") {
        parseAttr = true;
    }

    // Move pointer to general section
    if (pProf->inPerspect) {
        if (src_str == "ChromaticRedGreenModel") {
            pProf->pCurCommon = &pProf->pCurPersModel->chromRG;
            parseAttr = true;
        } else if (src_str == "ChromaticGreenModel") {
            pProf->pCurCommon = &pProf->pCurPersModel->chromG;
            parseAttr = true;
        } else if (src_str == "ChromaticBlueGreenModel") {
            pProf->pCurCommon = &pProf->pCurPersModel->chromBG;
            parseAttr = true;
        } else if (src_str == "VignetteModel") {
            pProf->pCurCommon = &pProf->pCurPersModel->vignette;
            parseAttr = true;
        }
    }

    // some profiles (espc. Pentax) have a different structure that is attributes based
    // simulate tags by feeding them in
    if (parseAttr && attr != nullptr) {
        for (int i = 0; attr[i]; i += 2) {
            const char* nameStart = strrchr(attr[i], ':');

            if (nameStart == nullptr) {
                nameStart = attr[i];
            } else {
                ++nameStart;
            }

            strncpy(pProf->lastTag, nameStart, 255);

            pProf->handle_text(attr[i + 1]);
        }
    }
}

void XMLCALL rtengine::LCPProfile::XmlTextHandler(void* pLCPProfile, const XML_Char* s, int len)
{
    LCPProfile* const pProf = static_cast<LCPProfile*>(pLCPProfile);

    if (
        !pProf->inCamProfiles
        || pProf->inAlternateLensID
        || pProf->inAlternateLensNames
        || *pProf->inInvalidTag
    ) {
        return;
    }

    for (int i = 0; i < len; ++i) {
        pProf->textbuf << s[i];
    }
}

void XMLCALL rtengine::LCPProfile::XmlEndHandler(void* pLCPProfile, const char* el)
{
    LCPProfile* const pProf = static_cast<LCPProfile*>(pLCPProfile);

    pProf->handle_text(pProf->textbuf.str());
    pProf->textbuf.str("");

    // We ignore everything in dirty tag till it's gone
    if (*pProf->inInvalidTag) {
        if (std::strstr(el, pProf->inInvalidTag)) {
            *pProf->inInvalidTag = 0;
        }

        return;
    }

    if (std::strstr(el, ":CameraProfiles")) {
        pProf->inCamProfiles = false;
    }

    if (std::strstr(el, ":AlternateLensIDs")) {
        pProf->inAlternateLensID = false;
    }

    if (std::strstr(el, ":AlternateLensNames")) {
        pProf->inAlternateLensNames = false;
    }

    if (
        !pProf->inCamProfiles
        || pProf->inAlternateLensID
        || pProf->inAlternateLensNames
    ) {
        return;
    }

    if (std::strstr(el, ":PerspectiveModel") || std::strstr(el, ":FisheyeModel")) {
        pProf->inPerspect = false;
    } else if (std::strstr(el, ":li")) {
        pProf->aPersModel[pProf->persModelCount] = pProf->pCurPersModel;
        pProf->pCurPersModel = nullptr;
        ++pProf->persModelCount;
    }
}

// Generates as singleton
rtengine::LCPStore* rtengine::LCPStore::getInstance()
{
    static LCPStore instance_;
    return &instance_;
}

bool rtengine::LCPStore::isValidLCPFileName(const Glib::ustring& filename) const
{
    if (!Glib::file_test(filename, Glib::FILE_TEST_EXISTS) || Glib::file_test (filename, Glib::FILE_TEST_IS_DIR)) {
        return false;
    }

    const size_t pos = filename.find_last_of ('.');
    return pos > 0 && !filename.casefold().compare(pos, 4, ".lcp");
}

std::shared_ptr<rtengine::LCPProfile> rtengine::LCPStore::getProfile(const Glib::ustring& filename) const
{
    if (filename.length() == 0 || !isValidLCPFileName(filename)) {
        return nullptr;
    }

    std::shared_ptr<LCPProfile> res;
    if (!cache.get(filename, res)) {
        res.reset(new LCPProfile(filename));
        cache.set(filename, res);
    }

    return res;
}

Glib::ustring rtengine::LCPStore::getDefaultCommonDirectory() const
{
    Glib::ustring dir;

#ifdef WIN32
    WCHAR pathW[MAX_PATH] = {0};

    if (SHGetSpecialFolderPathW(NULL, pathW, CSIDL_COMMON_APPDATA, false)) {
        char pathA[MAX_PATH];
        WideCharToMultiByte(CP_UTF8, 0, pathW, -1, pathA, MAX_PATH, 0, 0);
        Glib::ustring fullDir = Glib::ustring(pathA) + Glib::ustring("\\Adobe\\CameraRaw\\LensProfiles\\1.0");

        if (Glib::file_test (fullDir, Glib::FILE_TEST_IS_DIR)) {
            dir = fullDir;
        }
    }

#endif

    // TODO: Add Mac paths here

    return dir;
}

rtengine::LCPStore::LCPStore(unsigned int _cache_size) :
    cache(_cache_size)
{
}

// if !vignette then geometric and CA
rtengine::LCPMapper::LCPMapper(
    const std::shared_ptr<LCPProfile>& pProf,
    float focalLength,
    float focalLength35mm,
    float focusDist,
    float aperture,
    bool vignette,
    bool useCADistP,
    int fullWidth,
    int fullHeight,
    const CoarseTransformParams& coarse,
    int rawRotationDeg
) :
    enableCA(false),
    useCADist(useCADistP),
    swapXY(false),
    isFisheye(false)
{
    if (!pProf) {
        return;
    }

    // determine in what the image with the RAW landscape in comparison (calibration target)
    // in vignetting, the rotation has not taken place yet
    int rot = 0;

    if (rawRotationDeg >= 0) {
        rot = (coarse.rotate + rawRotationDeg) % 360;
    }

    swapXY = (rot == 90  || rot == 270);

    const bool mirrorX = (rot == 90  || rot == 180);
    const bool mirrorY = (rot == 180 || rot == 270);
    if (settings->verbose) {
        std::printf("Vign: %i, fullWidth: %i/%i, focLen %g SwapXY: %i / MirX/Y %i / %i on rot:%i from %i\n",vignette, fullWidth, fullHeight, focalLength, swapXY, mirrorX, mirrorY, rot, rawRotationDeg);
    }

    pProf->calcParams(vignette ? LCPCorrectionMode::VIGNETTE : LCPCorrectionMode::DISTORTION, focalLength, focusDist, aperture, &mc, nullptr, nullptr);
    mc.prepareParams(fullWidth, fullHeight, focalLength, focalLength35mm, pProf->sensorFormatFactor, swapXY, mirrorX, mirrorY);

    if (!vignette) {
        pProf->calcParams(LCPCorrectionMode::CA, focalLength, focusDist, aperture, &chrom[0], &chrom[1], &chrom[2]);

        for (int i = 0; i < 3; ++i) {
            chrom[i].prepareParams(fullWidth, fullHeight, focalLength, focalLength35mm, pProf->sensorFormatFactor, swapXY, mirrorX, mirrorY);
        }
    }

    enableCA = !vignette && focusDist > 0.f;
    isFisheye = pProf->isFisheye;
}

bool rtengine::LCPMapper::isCACorrectionAvailable() const
{
    return enableCA;
}

void rtengine::LCPMapper::correctDistortion(double &x, double &y, int cx, int cy, double scale) const
{
    x += cx;
    y += cy;

    if (isFisheye) {
        const double u = x * scale;
        const double v = y * scale;
        const double u0 = mc.x0 * scale;
        const double v0 = mc.y0 * scale;
        const double du = (u - u0);
        const double dv = (v - v0);
        const double fx = mc.fx;
        const double fy = mc.fy;
        const double k1 = mc.param[0];
        const double k2 = mc.param[1];
        const double r = sqrt(du * du + dv * dv);
        const double f = sqrt(fx*fy / (scale * scale));
        const double th = atan2(r, f);
        const double th2 = th * th;
        const double cfact = (((k2 * th2 + k1) * th2 + 1) * th) / r;
        const double ud = cfact * fx * du + u0;
        const double vd = cfact * fy * dv + v0;

        x = ud;
        y = vd;
    } else {
        x *= scale;
        y *= scale;
        const double x0 = mc.x0 * scale;
        const double y0 = mc.y0 * scale;
        const double xd = (x - x0) / mc.fx, yd = (y - y0) / mc.fy;

        const LCPModelCommon::Param aDist = mc.param;
        const double rsqr      = xd * xd + yd * yd;
        const double xfac = aDist[swapXY ? 3 : 4], yfac = aDist[swapXY ? 4 : 3];

        const double commonFac = (((aDist[2] * rsqr + aDist[1]) * rsqr + aDist[0]) * rsqr + 1.)
            + 2. * (yfac * yd + xfac * xd);

        const double xnew = xd * commonFac + xfac * rsqr;
        const double ynew = yd * commonFac + yfac * rsqr;

        x = xnew * mc.fx + x0;
        y = ynew * mc.fy + y0;
    }

    x -= cx * scale;
    y -= cy * scale;
}

void rtengine::LCPMapper::correctCA(double& x, double& y, int channel) const
{
    if (!enableCA) {
        return;
    }

    double xgreen, ygreen;

    // First calc the green channel like normal distortion
    // the other are just deviations from it
    double xd = (x - chrom[1].x0) / chrom[1].fx;
    double yd = (y - chrom[1].y0) / chrom[1].fy;

    // Green contains main distortion, just like base
    if (useCADist) {
        const LCPModelCommon::Param aDist = chrom[1].param;
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
        const double rsqr = xd * xd + yd * yd;

        const LCPModelCommon::Param aCA = chrom[channel].param;
        const double xfac = aCA[swapXY ? 3 : 4], yfac = aCA[swapXY ? 4 : 3];
        const double commonSum = 1. + rsqr * (aCA[0] + rsqr * (aCA[1] + aCA[2] * rsqr)) + 2. * (yfac * yd + xfac * xd);

        x = (chrom[channel].scale_factor * ( xd * commonSum + xfac * rsqr )) * chrom[channel].fx + chrom[channel].x0;
        y = (chrom[channel].scale_factor * ( yd * commonSum + yfac * rsqr )) * chrom[channel].fy + chrom[channel].y0;
    }
}

SSEFUNCTION void rtengine::LCPMapper::processVignetteLine(int width, int y, float* line) const
{
    // No need for swapXY, since vignette is in RAW and always before rotation
    float yd = ((float)y - mc.y0) * mc.rfy;
    yd *= yd;
    int x = 0;
#ifdef __SSE2__
    const vfloat fourv = F2V(4.f);
    const vfloat zerov = F2V(0.f);
    const vfloat ydv = F2V(yd);
    const vfloat p0 = F2V(mc.vign_param[0]);
    const vfloat p1 = F2V(mc.vign_param[1]);
    const vfloat p2 = F2V(mc.vign_param[2]);
    const vfloat p3 = F2V(mc.vign_param[3]);
    const vfloat x0v = F2V(mc.x0);
    const vfloat rfxv = F2V(mc.rfx);

    vfloat xv = _mm_setr_ps(0.f, 1.f, 2.f, 3.f);
    for (; x < width-3; x+=4) {
        const vfloat xdv = (xv - x0v) * rfxv;
        const vfloat rsqr = xdv * xdv + ydv;
        const vfloat vignFactorv = rsqr * (p0 + rsqr * (p1 - p2 * rsqr + p3 * rsqr * rsqr));
        vfloat valv = LVFU(line[x]);
        valv += valv * vselfzero(vmaskf_gt(valv, zerov), vignFactorv);
        STVFU(line[x], valv);
        xv += fourv;
    }
#endif // __SSE2__
    for (; x < width; x++) {
        if (line[x] > 0) {
            const float xd = ((float)x - mc.x0) * mc.rfx;
            const LCPModelCommon::VignParam vignParam = mc.vign_param;
            const float rsqr = xd * xd + yd;
            line[x] += line[x] * rsqr * (vignParam[0] + rsqr * ((vignParam[1]) - (vignParam[2]) * rsqr + (vignParam[3]) * rsqr * rsqr));
        }
    }
}

SSEFUNCTION void rtengine::LCPMapper::processVignetteLine3Channels(int width, int y, float* line) const
{
    // No need for swapXY, since vignette is in RAW and always before rotation
    float yd = ((float)y - mc.y0) * mc.rfy;
    yd *= yd;
    const LCPModelCommon::VignParam vignParam = mc.vign_param;
    for (int x = 0; x < width; x++) {
        const float xd = ((float)x - mc.x0) * mc.rfx;
        const float rsqr = xd * xd + yd;
        const float vignetteFactor = rsqr * (vignParam[0] + rsqr * ((vignParam[1]) - (vignParam[2]) * rsqr + (vignParam[3]) * rsqr * rsqr));
        for(int c = 0;c < 3; ++c) {
            if (line[3*x+c] > 0) {
                line[3*x+c] += line[3*x+c] * vignetteFactor;
            }
        }
    }
}
