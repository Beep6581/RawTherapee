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

#pragma once

#include <array>
#include <map>
#include <memory>
#include <string>
#include <sstream>

#include <glibmm.h>
#include <expat.h>

#include "cache.h"
#include "imagefloat.h"
#include "opthelper.h"

namespace rtengine
{

enum class LCPCorrectionMode {
    VIGNETTE,
    DISTORTION,
    CA
};

// Perspective model common data, also used for Vignette and Fisheye
class LCPModelCommon final
{
public:
    LCPModelCommon();

    bool empty() const;  // is it empty
    void print() const;  // printf all values
    void merge(const LCPModelCommon& a, const LCPModelCommon& b, float facA);
    void prepareParams(
        int fullWidth,
        int fullHeight,
        float focalLength,
        float focalLength35mm,
        float sensorFormatFactor,
        bool swapXY,
        bool mirrorX,
        bool mirrorY
    );

//private:
    using Param = std::array<float, 5>;
    using VignParam = std::array<float, 4>;

    float foc_len_x;
    float foc_len_y;
    float img_center_x;
    float img_center_y;
    Param param;  // k1..k5, resp. alpha1..5
    float scale_factor;  // alpha0
    double mean_error;
    bool bad_error;

    // prepared params
    float x0;
    float y0;
    float fx;
    float fy;
    float rfx;
    float rfy;
    VignParam vign_param;
};

class LCPProfile
{
public:
    explicit LCPProfile(const Glib::ustring& fname);
    ~LCPProfile();

    void calcParams(
        LCPCorrectionMode mode,
        float focalLength,
        float focusDist,
        float aperture,
        LCPModelCommon* pCorr1,
        LCPModelCommon* pCorr2,
        LCPModelCommon *pCorr3
    ) const; // Interpolates between the persModels frames

    void print() const;

//private:
    // Common data
    Glib::ustring profileName;
    Glib::ustring lensPrettyName;
    Glib::ustring cameraPrettyName;
    Glib::ustring lens;
    Glib::ustring camera;  // lens/camera(=model) can be auto-matched with DNG
    bool isRaw;
    bool isFisheye;
    float sensorFormatFactor;
    int persModelCount;

private:
    class LCPPersModel;

    int filterBadFrames(LCPCorrectionMode mode, double maxAvgDevFac, int minFramesLeft);

    void handle_text(const std::string& text);

    static void XMLCALL XmlStartHandler(void* pLCPProfile, const char* el, const char** attr);
    static void XMLCALL XmlTextHandler(void* pLCPProfile, const XML_Char* s, int len);
    static void XMLCALL XmlEndHandler(void* pLCPProfile, const char* el);

    // Temporary data for parsing
    bool inCamProfiles;
    bool firstLIDone;
    bool inPerspect;
    bool inAlternateLensID;
    bool inAlternateLensNames;
    char lastTag[256];
    char inInvalidTag[256];
    LCPPersModel* pCurPersModel;
    LCPModelCommon* pCurCommon;

    std::ostringstream textbuf;

    // The correction frames
    static constexpr int MaxPersModelCount = 3000;
    LCPPersModel* aPersModel[MaxPersModelCount];  // Do NOT use std::list or something, it's buggy in GCC!
};

class LCPStore
{
public:
    static LCPStore* getInstance();

    bool isValidLCPFileName(const Glib::ustring& filename) const;
    std::shared_ptr<LCPProfile> getProfile(const Glib::ustring& filename) const;
    Glib::ustring getDefaultCommonDirectory() const;

private:
    LCPStore(unsigned int _cache_size = 32);

    // Maps file name to profile as cache
    mutable Cache<Glib::ustring, std::shared_ptr<LCPProfile>> cache;
};

class LensCorrection {
public:
    virtual ~LensCorrection() {}
    virtual void correctDistortion(double &x, double &y, int cx, int cy, double scale) const = 0;
    virtual bool isCACorrectionAvailable() const = 0;
    virtual void correctCA(double &x, double &y, int channel) const = 0;
    virtual void processVignetteLine(int width, int y, int xoffs, int yoffs, float *line) const = 0;
    virtual void processVignetteLine3Channels(int width, int y, float *line) const = 0;
};


// Once precalculated class to correct a point
class LCPMapper: public LensCorrection
{
public:
    // Precalculates the mapper
    LCPMapper(
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
    );


    void correctDistortion(double &x, double &y, int cx, int cy, double scale) const;  // MUST be the first stage
    bool isCACorrectionAvailable() const;
    void correctCA(double& x, double& y, int channel) const;
    void processVignetteLine(int width, int y, int xoffs, int yoffs, float* line) const;
    void processVignetteLine3Channels(int width, int y, float* line) const;

private:
    bool enableCA;  // is the mapper capable if CA correction?
    bool useCADist;  // should the distortion in the CA info be used?
    bool swapXY;
    LCPModelCommon mc;
    LCPModelCommon chrom[3];  // in order RedGreen/Green/BlueGreen
    bool isFisheye;
};

}
