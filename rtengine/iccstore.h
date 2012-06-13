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

        std::map<std::string, cmsHPROFILE> wProfiles;
        std::map<std::string, cmsHPROFILE> wProfilesGamma;
        std::map<std::string, TMatrix> wMatrices;
        std::map<std::string, TMatrix> iwMatrices;
        
        // these contain profiles from user/system directory (supplied on init)
        std::map<std::string, cmsHPROFILE>    fileProfiles;
        std::map<std::string, ProfileContent> fileProfileContents;

        // these contain standard profiles from RT. keys are all in uppercase
        std::map<std::string, cmsHPROFILE>    fileStdProfiles;
        std::map<std::string, ProfileContent> fileStdProfileContents;

        cmsHPROFILE xyz;
        cmsHPROFILE srgb;

        Glib::Mutex mutex_;

        ICCStore (); 
        void loadICCs(Glib::ustring rootDirName, bool nameUpper, std::map<std::string, cmsHPROFILE>& resultProfiles, std::map<std::string, ProfileContent> &resultProfileContents);
        
    public:

        static ICCStore* getInstance(void);

        Glib::ustring defaultMonitorProfile;  // Main monitors standard profile name, from OS
        void findDefaultMonitorProfile();

        int     numOfWProfiles ();
        cmsHPROFILE createFromMatrix (const double matrix[3][3], bool gamma=false, Glib::ustring name="");
        cmsHPROFILE workingSpace (Glib::ustring name);
        cmsHPROFILE workingSpaceGamma (Glib::ustring name);
        TMatrix workingSpaceMatrix (Glib::ustring name);
        TMatrix workingSpaceInverseMatrix (Glib::ustring name);
        
        cmsHPROFILE              getProfile   (Glib::ustring name);
        cmsHPROFILE      getStdProfile(Glib::ustring name);

        void             init         (Glib::ustring usrICCDir, Glib::ustring stdICCDir);
        ProfileContent           getContent   (Glib::ustring name);

        cmsHPROFILE getXYZProfile ()  { return xyz;  }
        cmsHPROFILE getsRGBProfile () { return srgb; }
        std::vector<std::string> getOutputProfiles ();
};

#define iccStore ICCStore::getInstance()

//extern const char* wpnames[];
}
#endif
                                                      
