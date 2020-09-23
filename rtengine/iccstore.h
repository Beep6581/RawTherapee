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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <glibmm/ustring.h>

#include <lcms2.h>

namespace rtengine
{

namespace procparams
{

    struct ColorManagementParams;

}

typedef const double(*TMatrix)[3];

class ProfileContent final
{
public:
    ProfileContent();

    explicit ProfileContent(const Glib::ustring& fileName);
    explicit ProfileContent(cmsHPROFILE hProfile);
    cmsHPROFILE toProfile() const;

    const std::string& getData() const;

private:
    std::string data;
};

class ICCStore final
{
public:
    enum class ProfileType {
        MONITOR,
        PRINTER,
        OUTPUT  //(actually correspond to the same profiles than with MONITOR)
    };

    static ICCStore* getInstance();

    void init(const Glib::ustring& usrICCDir, const Glib::ustring& stdICCDir, bool loadAll);

    cmsHPROFILE      workingSpace(const Glib::ustring& name) const;
    // cmsHPROFILE      workingSpaceGamma(const Glib::ustring& name) const;
    TMatrix          workingSpaceMatrix(const Glib::ustring& name) const;
    TMatrix          workingSpaceInverseMatrix(const Glib::ustring& name) const;

    bool             outputProfileExist(const Glib::ustring& name) const;
    cmsHPROFILE      getProfile(const Glib::ustring& name) const;
    cmsHPROFILE      getStdProfile(const Glib::ustring& name) const;
    ProfileContent   getContent(const Glib::ustring& name) const;

    Glib::ustring getDefaultMonitorProfileName() const;
    void setDefaultMonitorProfileName(const Glib::ustring &name);

    cmsHPROFILE      getXYZProfile() const;
    cmsHPROFILE      getsRGBProfile() const;

    std::vector<Glib::ustring> getProfiles(ProfileType type = ProfileType::MONITOR) const;
    std::vector<Glib::ustring> getProfilesFromDir(const Glib::ustring& dirName) const;

    std::uint8_t     getInputIntents(cmsHPROFILE profile) const;
    std::uint8_t     getOutputIntents(cmsHPROFILE profile) const;
    std::uint8_t     getProofIntents(cmsHPROFILE profile) const;

    std::uint8_t     getInputIntents(const Glib::ustring& name) const;
    std::uint8_t     getOutputIntents(const Glib::ustring& name) const;
    std::uint8_t     getProofIntents(const Glib::ustring& name) const;

    /*static*/ std::vector<Glib::ustring> getWorkingProfiles();

    static cmsHPROFILE createFromMatrix(const double matrix[3][3], bool gamma = false, const Glib::ustring& name = Glib::ustring());

private:
    class Implementation;

    ICCStore();
    ~ICCStore();

    const std::unique_ptr<Implementation> implementation;
};

}
