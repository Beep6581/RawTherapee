/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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
#ifndef __ICCSTORE__
#define __ICCSTORE__

#include <lcms2.h>
#include <glibmm.h>
#include <map>
#include <string>
#include "../rtgui/threadutils.h"

namespace rtengine {

typedef const double (*TMatrix)[3];

class ProfileContent {

    public:
        char* data;
        int   length;

        ProfileContent (): data(NULL), length(0) {}
        ProfileContent (Glib::ustring fileName);
        ProfileContent (const ProfileContent& other);
        ProfileContent (cmsHPROFILE hProfile);
        ~ProfileContent ();
        ProfileContent& operator= (const rtengine::ProfileContent& other);
        cmsHPROFILE toProfile ();
};

class ICCStore {

        std::map<Glib::ustring, cmsHPROFILE> wProfiles;
        std::map<Glib::ustring, cmsHPROFILE> wProfilesGamma;
        std::map<Glib::ustring, TMatrix> wMatrices;
        std::map<Glib::ustring, TMatrix> iwMatrices;
        
        // these contain profiles from user/system directory (supplied on init)
        std::map<Glib::ustring, cmsHPROFILE>    fileProfiles;
        std::map<Glib::ustring, ProfileContent> fileProfileContents;

        // these contain standard profiles from RT. keys are all in uppercase
        std::map<Glib::ustring, cmsHPROFILE>    fileStdProfiles;
        std::map<Glib::ustring, ProfileContent> fileStdProfileContents;

        cmsHPROFILE xyz;
        cmsHPROFILE srgb;

        MyMutex mutex_;

        ICCStore (); 
        void             loadICCs(Glib::ustring rootDirName, bool nameUpper, std::map<Glib::ustring, cmsHPROFILE>& resultProfiles, std::map<Glib::ustring, ProfileContent> &resultProfileContents);
        
    public:

        static ICCStore* getInstance(void);
        static cmsHPROFILE makeStdGammaProfile(cmsHPROFILE iprof);

        Glib::ustring    defaultMonitorProfile;  // Main monitors standard profile name, from OS
        void             findDefaultMonitorProfile();

        int              numOfWProfiles ();
        cmsHPROFILE      createFromMatrix (const double matrix[3][3], bool gamma=false, Glib::ustring name="");
        cmsHPROFILE      workingSpace (Glib::ustring name);
        cmsHPROFILE      workingSpaceGamma (Glib::ustring name);
        TMatrix          workingSpaceMatrix (Glib::ustring name);
        TMatrix          workingSpaceInverseMatrix (Glib::ustring name);
        
        cmsHPROFILE      getProfile   (Glib::ustring name);
        cmsHPROFILE      getStdProfile(Glib::ustring name);

        void             init         (Glib::ustring usrICCDir, Glib::ustring stdICCDir);
        ProfileContent   getContent   (Glib::ustring name);

        cmsHPROFILE      getXYZProfile ()  { return xyz;  }
        cmsHPROFILE      getsRGBProfile () { return srgb; }
        std::vector<Glib::ustring> getOutputProfiles ();
};

#define iccStore ICCStore::getInstance()

//extern const char* wpnames[];
}
#endif
                                                      
