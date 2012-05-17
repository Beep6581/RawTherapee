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

#ifndef _LCP_
#define _LCP_

#include "imagefloat.h"
#include <glibmm.h>
#include <map>
#include <list>
#include <string>
#include <expat.h>

namespace rtengine {
    // Perspective model common data, also used for Vignette and Fisheye
    class LCPModelCommon {
    public:
       float focLenX, focLenY, imgXCenter, imgYCenter;
       float param[5];  // k1..k5
       float scaleFac;

       LCPModelCommon();
       bool empty() const;  // is it empty?
       void print() const;  // printf all values
    };

    class LCPPersModel {
    public:
        float focLen, focDist, aperture;  // this is what it refers to

        LCPModelCommon base;  // base perspective correction
        LCPModelCommon chromRG, chromG, chromBG;  // red/green, green, blue/green (may be empty)
        LCPModelCommon vignette;  // vignette (may be empty)

        LCPPersModel();
        void print() const;
    };


    class LCPProfile {
        // Temporary data for parsing
        bool inCamProfiles,firstLIDone,inPerspect;
        char lastTag[256];
        LCPPersModel* pCurPersModel;
        LCPModelCommon* pCurCommon;

        static void XMLCALL XmlStartHandler(void *pLCPProfile, const char *el, const char **attr);
        static void XMLCALL XmlTextHandler (void *pLCPProfile, const XML_Char *s, int len);
        static void XMLCALL XmlEndHandler  (void *pLCPProfile, const char *el);

    public:
        // Common data
        Glib::ustring profileName, lensPrettyName, cameraPrettyName, lens, camera;  // lens/camera(=model) can be auto-matched with DNG
        bool isRaw,isFisheye;

        // The correction frames
        std::list<LCPPersModel*> persModels;

        LCPProfile(Glib::ustring fname);

        void calcBasePerspectiveParams(float focalLength, bool vignette, LCPModelCommon& corr);  // Interpolates between the persModels frames

        void print() const;
    };

    class LCPStore {
        Glib::Mutex mtx;


        // Maps file name to profile as cache
        std::map<Glib::ustring, LCPProfile*> profileCache;

    public:
        Glib::ustring getDefaultCommonDirectory() const;
        bool isValidLCPFileName(Glib::ustring filename) const;
        LCPProfile* getProfile(Glib::ustring filename);

        static LCPStore* getInstance();
    };

#define lcpStore LCPStore::getInstance()


    // Once precalculated class to correct a point
    class LCPMapper {
        double x0,y0,fx,fy;
        LCPModelCommon mc;
        bool swapXY;
        bool mirrorX, mirrorY;

    public:
        // precalculates the mapper.
        LCPMapper(LCPProfile* pProf, float focalLength, bool vignette, int fullWidth, int fullHeight,
            const CoarseTransformParams& coarse, int rawRotationDeg);

        void  correctDistortion(double& x, double& y) const;  // MUST be the first stage
        float correctVignette  (double x, double y)   const;  // MUST be in RAW
    };
}
#endif
