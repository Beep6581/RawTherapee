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
#include <cstdint>
#include "procparams.h"
#include "color.h"
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

    explicit ProfileContent (const Glib::ustring& fileName);
    explicit ProfileContent (cmsHPROFILE hProfile);
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
    Glib::ustring userICCDir;
    ProfileMap fileProfiles;
    ContentMap fileProfileContents;

    // these contain standard profiles from RT. keys are all in uppercase
    Glib::ustring stdProfilesDir;
    NameMap fileStdProfilesFileNames;
    ProfileMap fileStdProfiles;

    Glib::ustring defaultMonitorProfile;

    bool loadAll;

    const cmsHPROFILE xyz;
    const cmsHPROFILE srgb;

    mutable MyMutex mutex_;

    ICCStore ();

public:

    enum class ProfileType {
        MONITOR,
        PRINTER,
        OUTPUT  // (actually correspond to the same profiles than with MONITOR)
    };

    static ICCStore* getInstance ();

    void init (const Glib::ustring& usrICCDir, const Glib::ustring& stdICCDir, bool loadAll);

    static void getGammaArray(const procparams::ColorManagementParams &icm, GammaValues &ga);
    static cmsHPROFILE makeStdGammaProfile (cmsHPROFILE iprof);
    static cmsHPROFILE createFromMatrix (const double matrix[3][3], bool gamma = false, const Glib::ustring& name = Glib::ustring());
    static cmsHPROFILE createGammaProfile (const procparams::ColorManagementParams &icm, GammaValues &ga);
    static cmsHPROFILE createCustomGammaOutputProfile (const procparams::ColorManagementParams &icm, GammaValues &ga);

    // Main monitors standard profile name, from OS
    void findDefaultMonitorProfile ();
    cmsHPROFILE getDefaultMonitorProfile ();
    Glib::ustring getDefaultMonitorProfileName ();

    cmsHPROFILE      workingSpace (const Glib::ustring& name) const;
    cmsHPROFILE      workingSpaceGamma (const Glib::ustring& name) const;
    TMatrix          workingSpaceMatrix (const Glib::ustring& name) const;
    TMatrix          workingSpaceInverseMatrix (const Glib::ustring& name) const;

    bool             outputProfileExist (const Glib::ustring& name) const;
    cmsHPROFILE      getProfile         (const Glib::ustring& name);
    cmsHPROFILE      getStdProfile      (const Glib::ustring& name);
    ProfileContent   getContent         (const Glib::ustring& name) const;

    cmsHPROFILE      getXYZProfile  () const;
    cmsHPROFILE      getsRGBProfile () const;

    std::vector<Glib::ustring> getProfiles (const ProfileType type = ProfileType::MONITOR) const;
    std::vector<Glib::ustring> getProfilesFromDir (const Glib::ustring& dirName) const;

    uint8_t     getInputIntents  (cmsHPROFILE profile);
    uint8_t     getOutputIntents (cmsHPROFILE profile);
    uint8_t     getProofIntents  (cmsHPROFILE profile);

    uint8_t     getInputIntents  (const Glib::ustring& name);
    uint8_t     getOutputIntents (const Glib::ustring& name);
    uint8_t     getProofIntents  (const Glib::ustring& name);
};

#define iccStore ICCStore::getInstance()

inline ProfileContent::ProfileContent () :
    data(nullptr),
    length(0)
{
}

inline ProfileContent::~ProfileContent ()
{
    delete [] data;
}

inline cmsHPROFILE ICCStore::getDefaultMonitorProfile ()
{
    return getProfile (defaultMonitorProfile);
}

inline Glib::ustring ICCStore::getDefaultMonitorProfileName ()
{
    return defaultMonitorProfile;
}

inline uint8_t ICCStore::getInputIntents (const Glib::ustring &name)
{
    return getInputIntents (getProfile (name));
}

inline uint8_t ICCStore::getOutputIntents (const Glib::ustring &name)
{
    return getOutputIntents (getProfile (name));
}

inline uint8_t ICCStore::getProofIntents (const Glib::ustring &name)
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

