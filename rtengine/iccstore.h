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

namespace rtengine
{

typedef const double (*TMatrix)[3];

class ProfileContent
{

public:
    char* data;
    int   length;

    ProfileContent ();
    ~ProfileContent ();

    ProfileContent (const ProfileContent& other);
    ProfileContent& operator= (const rtengine::ProfileContent& other);

    ProfileContent (const Glib::ustring& fileName);
    ProfileContent (cmsHPROFILE hProfile);
    cmsHPROFILE toProfile () const;
};

class ICCStore
{
    typedef std::map<Glib::ustring, cmsHPROFILE> ProfileMap;
    typedef std::map<Glib::ustring, TMatrix> MatrixMap;
    typedef std::map<Glib::ustring, ProfileContent> ContentMap;
    typedef std::map<Glib::ustring, Glib::ustring> NameMap;

    ProfileMap wProfiles;
    ProfileMap wProfilesGamma;
    MatrixMap wMatrices;
    MatrixMap iwMatrices;

    // these contain profiles from user/system directory (supplied on init)
    Glib::ustring profilesDir;
    ProfileMap fileProfiles;
    ContentMap fileProfileContents;

    // these contain standard profiles from RT. keys are all in uppercase
    Glib::ustring stdProfilesDir;
    NameMap fileStdProfilesFileNames;
    ProfileMap fileStdProfiles;

    Glib::ustring defaultMonitorProfile;

    const cmsHPROFILE xyz;
    const cmsHPROFILE srgb;

    mutable MyMutex mutex_;

    ICCStore ();

public:

    static ICCStore* getInstance ();

    void init (const Glib::ustring& usrICCDir, const Glib::ustring& stdICCDir);

    static cmsHPROFILE makeStdGammaProfile (cmsHPROFILE iprof);
    static cmsHPROFILE createFromMatrix (const double matrix[3][3], bool gamma = false, const Glib::ustring& name = Glib::ustring());

    // Main monitors standard profile name, from OS
    void findDefaultMonitorProfile ();
    cmsHPROFILE getDefaultMonitorProfile () const;

    cmsHPROFILE      workingSpace (const Glib::ustring& name) const;
    cmsHPROFILE      workingSpaceGamma (const Glib::ustring& name) const;
    TMatrix          workingSpaceMatrix (const Glib::ustring& name) const;
    TMatrix          workingSpaceInverseMatrix (const Glib::ustring& name) const;

    cmsHPROFILE      getProfile    (const Glib::ustring& name) const;
    cmsHPROFILE      getStdProfile (const Glib::ustring& name) const;
    ProfileContent   getContent    (const Glib::ustring& name) const;

    cmsHPROFILE      getXYZProfile  () const;
    cmsHPROFILE      getsRGBProfile () const;

    std::vector<Glib::ustring> getProfiles () const;
    std::vector<Glib::ustring> getProfilesFromDir (const Glib::ustring& dirName) const;

    std::uint8_t     getInputIntents  (cmsHPROFILE profile) const;
    std::uint8_t     getOutputIntents (cmsHPROFILE profile) const;
    std::uint8_t     getProofIntents  (cmsHPROFILE profile) const;

    std::uint8_t     getInputIntents  (const Glib::ustring& name) const;
    std::uint8_t     getOutputIntents (const Glib::ustring& name) const;
    std::uint8_t     getProofIntents  (const Glib::ustring& name) const;
};

#define iccStore ICCStore::getInstance()

inline ProfileContent::ProfileContent () :
    data(NULL),
    length(0)
{
}

inline ProfileContent::~ProfileContent ()
{
    delete [] data;
}

inline cmsHPROFILE ICCStore::getDefaultMonitorProfile () const
{
    return getProfile (defaultMonitorProfile);
}

inline std::uint8_t ICCStore::getInputIntents (const Glib::ustring &name) const
{
    return getInputIntents (getProfile (name));
}

inline std::uint8_t ICCStore::getOutputIntents (const Glib::ustring &name) const
{
    return getOutputIntents (getProfile (name));
}

inline std::uint8_t ICCStore::getProofIntents (const Glib::ustring &name) const
{
    return getProofIntents (getProfile (name));
}

inline cmsHPROFILE ICCStore::getXYZProfile () const
{
    return xyz;
}

inline cmsHPROFILE ICCStore::getsRGBProfile () const
{
    return srgb;
}

}

#endif

