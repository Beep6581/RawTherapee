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

#include <map>

#include <glibmm.h>

#include "imagefloat.h"
#include "curves.h"
#include "colortemp.h"

#include "../rtgui/threadutils.h"

namespace rtengine
{

class DCPProfile final
{
public:
    struct dcpApplyState{
        double m2ProPhoto[3][3];
        double m2Work[3][3];
        bool alreadyProPhoto;
        bool useToneCurve;
        bool applyLookTable;
        float blScale;
    };

    DCPProfile(const Glib::ustring &fname);
    ~DCPProfile();

    bool getHasToneCurve() const;
    bool getHasLookTable() const;
    bool getHasHueSatMap() const;
    bool getHasBaselineExposureOffset() const;
    void getIlluminants(int &i1, double &temp1, int &i2, double &temp2, bool &willInterpolate_) const;

    void Apply(Imagefloat *pImg, int preferredIlluminant, const Glib::ustring &workingSpace, const ColorTemp &wb, double pre_mul[3], double camMatrix[3][3], bool useToneCurve = false, bool applyHueSatMap = true, bool applyLookTable = false) const;
    void setStep2ApplyState(const Glib::ustring &workingSpace, bool useToneCurve, bool applyLookTable, bool applyBaselineExposure, dcpApplyState &asOut);
    void step2ApplyTile(float *r, float *g, float *b, int width, int height, int tileWidth, const dcpApplyState &asIn) const;

private:
    struct HSBModify {
        float fHueShift;
        float fSatScale;
        float fValScale;
    };
    struct HSDTableInfo {
        int iHueDivisions, iSatDivisions, iValDivisions;
        int iHueStep, iValStep, iArrayCount;
        bool sRGBGamma;
        struct {
            float hScale, sScale, vScale;
            int maxHueIndex0, maxSatIndex0, maxValIndex0;
            int hueStep, valStep;
        } pc;
    };

    void dngref_XYCoord2Temperature(const double whiteXY[2], double *temp, double *tint) const;
    void dngref_FindXYZtoCamera(const double whiteXY[2], int preferredIlluminant, double (*xyzToCamera)[3]) const;
    void dngref_NeutralToXY(double neutral[3], int preferredIlluminant, double XY[2]) const;
    void MakeXYZCAM(const ColorTemp &wb, double pre_mul[3], double camWbMatrix[3][3], int preferredIlluminant, double (*mXYZCAM)[3]) const;
    const HSBModify* MakeHueSatMap(const ColorTemp &wb, int preferredIlluminant, HSBModify **deleteHandle) const;
    void HSDApply(const HSDTableInfo &ti, const HSBModify *tableBase, float &h, float &s, float &v) const;

    double mColorMatrix1[3][3], mColorMatrix2[3][3];
    bool hasColorMatrix1, hasColorMatrix2, hasForwardMatrix1, hasForwardMatrix2, hasToneCurve, hasBaselineExposureOffset, willInterpolate;
    double mForwardMatrix1[3][3], mForwardMatrix2[3][3];
    double temperature1, temperature2;
    double baselineExposureOffset;
    HSBModify *aDeltas1, *aDeltas2, *aLookTable;
    HSDTableInfo DeltaInfo, LookInfo;
    short iLightSource1, iLightSource2;

    AdobeToneCurve toneCurve;
};

class DCPStore final
{
public:
    static DCPStore* getInstance();

    DCPStore(const DCPStore& other) = delete;
    DCPStore& operator =(const DCPStore& other) = delete;

    void init(const Glib::ustring& rt_profile_dir);

    bool isValidDCPFileName(const Glib::ustring& filename) const;

    DCPProfile* getProfile(const Glib::ustring& filename) const;
    DCPProfile* getStdProfile(const Glib::ustring& camShortName) const;

private:
    DCPStore() = default;

    mutable MyMutex mutex;

    // these contain standard profiles from RT. keys are all in uppercase, file path is value
    std::map<Glib::ustring, Glib::ustring> file_std_profiles;

    // Maps file name to profile as cache
    mutable std::map<Glib::ustring, DCPProfile*> profile_cache;
};

}
