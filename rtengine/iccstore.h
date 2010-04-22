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

#include <lcms.h>
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
        ~ProfileContent ();
        ProfileContent& operator= (const ProfileContent other);
        cmsHPROFILE toProfile ();
};

class ICCStore {

        std::map<std::string, cmsHPROFILE> wProfiles;
        std::map<std::string, cmsHPROFILE> wProfilesGamma;
        std::map<std::string, TMatrix> wMatrices;
        std::map<std::string, TMatrix> iwMatrices;
        
        std::map<std::string, cmsHPROFILE>    fileProfiles;
        std::map<std::string, ProfileContent> fileProfileContents;

        cmsHPROFILE xyz;
        cmsHPROFILE srgb;
        
    public:
        ICCStore (); 
        
        int     numOfWProfiles ();
        cmsHPROFILE createFromMatrix (const double matrix[3][3], bool gamma=false, Glib::ustring name="");
        cmsHPROFILE workingSpace (Glib::ustring name);
        cmsHPROFILE workingSpaceGamma (Glib::ustring name);
        TMatrix workingSpaceMatrix (Glib::ustring name);
        TMatrix workingSpaceInverseMatrix (Glib::ustring name);
        
        cmsHPROFILE              getProfile   (Glib::ustring name);
        std::vector<std::string> parseDir     (Glib::ustring pdir);
        ProfileContent           getContent   (Glib::ustring name);

        cmsHPROFILE getXYZProfile ()  { return xyz;  }
        cmsHPROFILE getsRGBProfile () { return srgb; }
        std::vector<std::string> getOutputProfiles ();
};

extern ICCStore iccStore;

//extern const char* wpnames[];
}
#endif
                                                      
