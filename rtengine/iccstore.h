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
#include <rtengine.h>
#include <map>
#include <string>
#include "matrix33.h"

namespace rtengine {

//typedef const float (*TMatrix)[3];

class ProfileContent {

    public:
        char* data;
        int   length;

        ProfileContent (): data(NULL), length(0) {}
        ProfileContent (const String& fileName);
        ProfileContent (const ProfileContent& other);       
        ~ProfileContent ();
        ProfileContent& operator= (const ProfileContent other);
        cmsHPROFILE toProfile ();
};

class ICCStore {

        std::map<std::string, cmsHPROFILE> wProfiles;
        std::map<std::string, cmsHPROFILE> wProfilesGamma;
        std::map<std::string, Matrix33> wMatrices;
        std::map<std::string, Matrix33> iwMatrices;
        
        std::map<std::string, cmsHPROFILE>    fileProfiles;
        std::map<std::string, ProfileContent> fileProfileContents;

        cmsHPROFILE xyz;
        cmsHPROFILE srgb;

        
    public:
        ICCStore (); 
        
        std::vector<std::string> getWorkingProfiles ();			/// returns the vector of names of the working profiles available
        cmsHPROFILE workingSpace (const String& name);			/// returns the profile corresponding to the given working color space
        cmsHPROFILE workingSpaceGamma (const String& name); 	/// returns the profile corresponding to the given working color space, with gamma correction
        Matrix33 workingSpaceMatrix (const String& name);		/// returns the transformation matrix belonging to the given working profile
        Matrix33 workingSpaceInverseMatrix (const String& name); /// returns the inverse transformation matrix belonging to the given working profile

        cmsHPROFILE createFromMatrix (const Matrix33& matrix, bool gamma=false, const String& name="");	/// create profile from matrix

        std::vector<std::string> parseDir     (const String& pdir);		/// parse the given directory and load all icc files found
        std::vector<std::string> getOutputProfiles ();					/// returns the list of names of all icc profiles loaded
        cmsHPROFILE              getProfile   (const String& name);		/// returns profile corresponding to the given name
        ProfileContent           getContent   (const String& name);		/// returns the content corresponding to profile "name" (can be used to attach to image files)

        cmsHPROFILE getXYZProfile ()  { return xyz;  }	/// return identity transform
        cmsHPROFILE getsRGBProfile () { return srgb; }  /// return standard srgb transform
};

extern ICCStore* iccStore;
}
#endif
                                                      
