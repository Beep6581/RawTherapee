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

#ifndef _DCP_
#define _DCP_

#include "imagefloat.h"
#include <glibmm.h>
#include <map>
#include <string>

namespace rtengine {
    enum DCPLightType {
        Daylight=1, Tungsten=2, Fluorescent=3, Flash=4
    };

    class DCPProfile {
        struct HSBModify
        {
            float fHueShift;
            float fSatScale;
            float fValScale;
        };

        double mColorMatrix1[3][3],mColorMatrix2[3][3];
        double mXYZCAM1[3][3],mXYZCAM2[3][3];  // compatible to RTs xyz_cam 
        HSBModify *aDeltas1,*aDeltas2;
        short iLightSource1,iLightSource2;

        int iHueDivisions, iSatDivisions, iValDivisions;

        int iHueStep, iValStep, iArrayCount;

        LUTf lutToneCurve;  // 0..0xffff values to 0..1
        void ApplyToneCurve(float& r, float& g, float& b) const;
        void RGBTone(float& r, float& g, float& b) const;  // helper for tone curve

        void ConvertDNGMatrix2XYZCAM(const double (*mColorMatrix)[3], double (*mXYZCAM)[3]);

        const HSBModify* GetBestProfile(DCPLightType preferredProfile,  double (*mXYZCAM)[3]) const;

        DCPLightType GetLightType(short iLightSource) const;

    public:
        DCPProfile(Glib::ustring fname, bool isRTProfile);
        ~DCPProfile();

        void Apply(Imagefloat *pImg, DCPLightType preferredProfile, Glib::ustring workingSpace, float rawWhiteFac=1, bool useToneCurve=false) const;
        void Apply(Image16 *pImg, DCPLightType preferredProfile, Glib::ustring workingSpace, bool useToneCurve) const;
    };

    class DCPStore {
        Glib::Mutex mtx;

        // these contain standard profiles from RT. keys are all in uppercase, file path is value
        std::map<Glib::ustring, Glib::ustring> fileStdProfiles;
        
        // Maps file name to profile as cache
        std::map<Glib::ustring, DCPProfile*> profileCache;

    public:
        void init(Glib::ustring rtProfileDir);

        bool isValidDCPFileName(Glib::ustring filename) const;

        DCPProfile* getProfile(Glib::ustring filename, bool isRTProfile=false);
        DCPProfile* getStdProfile(Glib::ustring camShortName);
        
        static DCPStore* getInstance();
    };

    #define dcpStore DCPStore::getInstance()
}
#endif
