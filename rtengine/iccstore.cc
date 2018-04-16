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
#include <cstring>

#include <glibmm.h>
#include <glib/gstdio.h>

#ifdef WIN32
#include <winsock2.h>
#else
#include <netinet/in.h>
#endif

#include <iostream>

#include "iccstore.h"

#include "iccmatrices.h"
#include "procparams.h"

#include "../rtgui/options.h"
#include "../rtgui/threadutils.h"

#include "cJSON.h"
namespace rtengine
{

extern const Settings* settings;

}

namespace
{

// Not recursive
void loadProfiles(
    const Glib::ustring& dirName,
    std::map<Glib::ustring, cmsHPROFILE>* profiles,
    std::map<Glib::ustring, rtengine::ProfileContent>* profileContents,
    std::map<Glib::ustring, Glib::ustring>* profileNames,
    bool nameUpper
)
{
    if (dirName.empty()) {
        return;
    }

    try {
        Glib::Dir dir(dirName);

        for (Glib::DirIterator entry = dir.begin(); entry != dir.end(); ++entry) {
            const Glib::ustring fileName = *entry;

            if (fileName.size() < 4) {
                continue;
            }

            const Glib::ustring extension = rtengine::getFileExtension(fileName);

            if (extension != "icc" && extension != "icm") {
                continue;
            }

            const Glib::ustring filePath = Glib::build_filename(dirName, fileName);

            if (!Glib::file_test(filePath, Glib::FILE_TEST_IS_REGULAR)) {
                continue;
            }

            Glib::ustring name = fileName.substr(0, fileName.size() - 4);

            if (nameUpper) {
                name = name.uppercase();
            }

            if (profiles) {
                const rtengine::ProfileContent content(filePath);
                const cmsHPROFILE profile = content.toProfile();

                if (profile) {
                    profiles->emplace(name, profile);

                    if (profileContents) {
                        profileContents->emplace(name, content);
                    }
                }
            }

            if (profileNames) {
                profileNames->emplace(name, filePath);
            }
        }
    } catch (Glib::Exception&) {
    }
}

// Version dedicated to single profile load when loadAll==false (cli version "-q" mode)
bool loadProfile(
    const Glib::ustring& profile,
    const Glib::ustring& dirName,
    std::map<Glib::ustring, cmsHPROFILE>* profiles,
    std::map<Glib::ustring, rtengine::ProfileContent>* profileContents
)
{
    if (dirName.empty() || profiles == nullptr) {
        return false;
    }

    try {
        Glib::Dir dir(dirName);

        for (Glib::DirIterator entry = dir.begin(); entry != dir.end(); ++entry) {
            const Glib::ustring fileName = *entry;

            if (fileName.size() < 4) {
                continue;
            }

            const Glib::ustring extension = rtengine::getFileExtension(fileName);

            if (extension != "icc" && extension != "icm") {
                continue;
            }

            const Glib::ustring filePath = Glib::build_filename(dirName, fileName);

            if (!Glib::file_test(filePath, Glib::FILE_TEST_IS_REGULAR)) {
                continue;
            }

            const Glib::ustring name = fileName.substr(0, fileName.size() - 4);

            if (name == profile) {
                const rtengine::ProfileContent content(filePath);
                const cmsHPROFILE profile = content.toProfile();

                if (profile) {
                    profiles->emplace(name, profile);

                    if (profileContents) {
                        profileContents->emplace(name, content);
                    }

                    return true;
                }
            }
        }
    } catch (Glib::Exception&) {
    }

    return false;
}

void getSupportedIntent(cmsHPROFILE profile, cmsUInt32Number intent, cmsUInt32Number direction, uint8_t& result)
{
    if (cmsIsIntentSupported(profile, intent, direction)) {
        result |= 1 << intent;
    }
}

uint8_t getSupportedIntents(cmsHPROFILE profile, cmsUInt32Number direction)
{
    if (!profile) {
        return 0;
    }

    uint8_t result = 0;

    getSupportedIntent(profile, INTENT_PERCEPTUAL, direction, result);
    getSupportedIntent(profile, INTENT_RELATIVE_COLORIMETRIC, direction, result);
    getSupportedIntent(profile, INTENT_SATURATION, direction, result);
    getSupportedIntent(profile, INTENT_ABSOLUTE_COLORIMETRIC, direction, result);

    return result;
}

cmsHPROFILE createXYZProfile()
{
    double mat[3][3] = { {1.0, 0, 0}, {0, 1.0, 0}, {0, 0, 1.0} };
    return rtengine::ICCStore::createFromMatrix(mat, false, "XYZ");
}

const double(*wprofiles[])[3]  = {xyz_sRGB, xyz_adobe, xyz_prophoto, xyz_widegamut, xyz_bruce, xyz_beta, xyz_best, xyz_rec2020, xyz_ACESp0, xyz_ACESp1};//
const double(*iwprofiles[])[3] = {sRGB_xyz, adobe_xyz, prophoto_xyz, widegamut_xyz, bruce_xyz, beta_xyz, best_xyz, rec2020_xyz, ACESp0_xyz, ACESp1_xyz};//
const char* wpnames[] = {"sRGB", "Adobe RGB", "ProPhoto", "WideGamut", "BruceRGB", "Beta RGB", "BestRGB", "Rec2020", "ACESp0", "ACESp1"};//
const char* wpgamma[] = {"Free", "BT709_g2.2_s4.5", "sRGB_g2.4_s12.92", "linear_g1.0", "standard_g2.2", "standard_g1.8", "High_g1.3_s3.35", "Low_g2.6_s6.9", "Lab_g3.0s9.03296"}; //gamma free
//default = gamma inside profile
//BT709 g=2.22 s=4.5  sRGB g=2.4 s=12.92
//linear g=1.0
//std22 g=2.2   std18 g=1.8
// high  g=1.3 s=3.35  for high dynamic images
//low  g=2.6 s=6.9  for low contrast images

}

rtengine::ProfileContent::ProfileContent() = default;

rtengine::ProfileContent::ProfileContent(const Glib::ustring& fileName)
{
    FILE* const f = g_fopen(fileName.c_str(), "rb");

    if (!f) {
        return;
    }

    fseek(f, 0, SEEK_END);
    long length = ftell(f);

    if (length > 0) {
        char* d = new char[length + 1];
        fseek(f, 0, SEEK_SET);
        length = fread(d, 1, length, f);
        d[length] = 0;
        data.assign(d, length);
        delete[] d;
    } else {
        data.clear();
    }

    fclose(f);

}

rtengine::ProfileContent::ProfileContent(cmsHPROFILE hProfile)
{
    if (hProfile != nullptr) {
        cmsUInt32Number bytesNeeded = 0;
        cmsSaveProfileToMem(hProfile, nullptr, &bytesNeeded);

        if (bytesNeeded > 0) {
            char* d = new char[bytesNeeded + 1];
            cmsSaveProfileToMem(hProfile, d, &bytesNeeded);
            data.assign(d, bytesNeeded);
            delete[] d;
        }
    }
}

cmsHPROFILE rtengine::ProfileContent::toProfile() const
{

    return
        !data.empty()
        ? cmsOpenProfileFromMem(data.c_str(), data.size())
        : nullptr;
}

const std::string& rtengine::ProfileContent::getData() const
{
    return data;
}

class rtengine::ICCStore::Implementation
{
public:
    Implementation() :
        loadAll(true),
        xyz(createXYZProfile()),
        srgb(cmsCreate_sRGBProfile())
    {
        //cmsErrorAction(LCMS_ERROR_SHOW);

        constexpr int N = sizeof(wpnames) / sizeof(wpnames[0]);

        for (int i = 0; i < N; ++i) {
            wProfiles[wpnames[i]] = createFromMatrix(wprofiles[i]);
            // wProfilesGamma[wpnames[i]] = createFromMatrix(wprofiles[i], true);
            wMatrices[wpnames[i]] = wprofiles[i];
            iwMatrices[wpnames[i]] = iwprofiles[i];
        }
    }

    ~Implementation()
    {
        for (auto &p : wProfiles) {
            if (p.second) {
                cmsCloseProfile(p.second);
            }
        }

        // for (auto &p : wProfilesGamma) {
        //     if (p.second) {
        //         cmsCloseProfile(p.second);
        //     }
        // }
        for (auto &p : fileProfiles) {
            if (p.second) {
                cmsCloseProfile(p.second);
            }
        }

        if (srgb) {
            cmsCloseProfile(srgb);
        }

        if (xyz) {
            cmsCloseProfile(xyz);
        }
    }

    void init(const Glib::ustring& usrICCDir, const Glib::ustring& rtICCDir, bool loadAll)
    {
        // Reads all profiles from the given profiles dir

        MyMutex::MyLock lock(mutex);

        this->loadAll = loadAll;

        // RawTherapee's profiles take precedence if a user's profile of the same name exists
        profilesDir = Glib::build_filename(rtICCDir, "output");
        userICCDir = usrICCDir;
        fileProfiles.clear();
        fileProfileContents.clear();

        if (loadAll) {
            loadProfiles(profilesDir, &fileProfiles, &fileProfileContents, nullptr, false);
            loadProfiles(userICCDir, &fileProfiles, &fileProfileContents, nullptr, false);
        }

        // Input profiles
        // Load these to different areas, since the short name(e.g. "NIKON D700" may overlap between system/user and RT dir)
        stdProfilesDir = Glib::build_filename(rtICCDir, "input");
        fileStdProfiles.clear();
        fileStdProfilesFileNames.clear();

        if (loadAll) {
            loadProfiles(stdProfilesDir, nullptr, nullptr, &fileStdProfilesFileNames, true);
        }

        defaultMonitorProfile = settings->monitorProfile;

        loadWorkingSpaces(rtICCDir);
        loadWorkingSpaces(userICCDir);

        // initialize the alarm colours for lcms gamut checking -- we use bright green
        cmsUInt16Number cms_alarm_codes[cmsMAXCHANNELS] = { 0, 65535, 65535 };
        cmsSetAlarmCodes(cms_alarm_codes);
    }

    cmsHPROFILE workingSpace(const Glib::ustring& name) const
    {
        const ProfileMap::const_iterator r = wProfiles.find(name);

        return
            r != wProfiles.end()
            ? r->second
            : wProfiles.find("sRGB")->second;
    }

    // cmsHPROFILE workingSpaceGamma(const Glib::ustring& name) const
    // {

    //     const ProfileMap::const_iterator r = wProfilesGamma.find(name);

    //     return
    //         r != wProfilesGamma.end()
    //             ? r->second
    //             : wProfilesGamma.find("sRGB")->second;
    // }

    TMatrix workingSpaceMatrix(const Glib::ustring& name) const
    {
        const MatrixMap::const_iterator r = wMatrices.find(name);

        return
            r != wMatrices.end()
            ? r->second
            : wMatrices.find("sRGB")->second;
    }

    TMatrix workingSpaceInverseMatrix(const Glib::ustring& name) const
    {

        const MatrixMap::const_iterator r = iwMatrices.find(name);

        return
            r != iwMatrices.end()
            ? r->second
            : iwMatrices.find("sRGB")->second;
    }

    bool outputProfileExist(const Glib::ustring& name) const
    {
        MyMutex::MyLock lock(mutex);
        return fileProfiles.find(name) != fileProfiles.end();
    }

    cmsHPROFILE getProfile(const Glib::ustring& name)
    {
        MyMutex::MyLock lock(mutex);

        const ProfileMap::const_iterator r = fileProfiles.find(name);

        if (r != fileProfiles.end()) {
            return r->second;
        }

        if (!name.compare(0, 5, "file:")) {
            const ProfileContent content(name.substr(5));
            const cmsHPROFILE profile = content.toProfile();

            if (profile) {
                fileProfiles.emplace(name, profile);
                fileProfileContents.emplace(name, content);

                return profile;
            }
        } else if (!loadAll) {
            // Look for a standard profile
            if (!loadProfile(name, profilesDir, &fileProfiles, &fileProfileContents)) {
                loadProfile(name, userICCDir, &fileProfiles, &fileProfileContents);
            }

            const ProfileMap::const_iterator r = fileProfiles.find(name);

            if (r != fileProfiles.end()) {
                return r->second;
            }
        }

        return nullptr;
    }

    cmsHPROFILE getStdProfile(const Glib::ustring& name)
    {
        const Glib::ustring nameUpper = name.uppercase();

        MyMutex::MyLock lock(mutex);

        const ProfileMap::const_iterator r = fileStdProfiles.find(nameUpper);

        // Return profile from store
        if (r != fileStdProfiles.end()) {
            return r->second;
        } else if (!loadAll) {
            // Directory not scanned, so looking and adding now...
            if (!loadProfile(name, profilesDir, &fileProfiles, &fileProfileContents)) {
                loadProfile(name, userICCDir, &fileProfiles, &fileProfileContents);
            }

            const ProfileMap::const_iterator r = fileProfiles.find(name);

            if (r != fileProfiles.end()) {
                return r->second;
            }
        }

        // Profile is not yet in store
        const NameMap::const_iterator f = fileStdProfilesFileNames.find(nameUpper);

        // Profile does not exist
        if (f == fileStdProfilesFileNames.end()) {
            return nullptr;
        }

        // But there exists one --> load it
        const ProfileContent content(f->second);
        const cmsHPROFILE profile = content.toProfile();

        if (profile) {
            fileStdProfiles.emplace(f->first, profile);
        }

        // Profile invalid or stored now --> remove entry from fileStdProfilesFileNames
        fileStdProfilesFileNames.erase(f);
        return profile;
    }

    ProfileContent getContent(const Glib::ustring& name) const
    {
        MyMutex::MyLock lock(mutex);

        const ContentMap::const_iterator r = fileProfileContents.find(name);

        return
            r != fileProfileContents.end()
            ? r->second
            : ProfileContent();
    }

    cmsHPROFILE getXYZProfile() const
    {
        return xyz;
    }

    cmsHPROFILE getsRGBProfile() const
    {
        return srgb;
    }

    std::vector<Glib::ustring> getProfiles(ProfileType type) const
    {
        std::vector<Glib::ustring> res;

        MyMutex::MyLock lock(mutex);

        for (const auto profile : fileProfiles) {
            if (
                (
                    type == ICCStore::ProfileType::MONITOR
                    && cmsGetDeviceClass(profile.second) == cmsSigDisplayClass
                    && cmsGetColorSpace(profile.second) == cmsSigRgbData
                )
                || (
                    type == ICCStore::ProfileType::PRINTER
                    && cmsGetDeviceClass(profile.second) == cmsSigOutputClass
                )
                || (
                    type == ICCStore::ProfileType::OUTPUT
                    && (cmsGetDeviceClass(profile.second) == cmsSigDisplayClass
                        || cmsGetDeviceClass(profile.second) == cmsSigInputClass
                        || cmsGetDeviceClass(profile.second) == cmsSigOutputClass)
                    && cmsGetColorSpace(profile.second) == cmsSigRgbData
                )
            ) {
                res.push_back(profile.first);
            }
        }

        return res;
    }

    std::vector<Glib::ustring> getProfilesFromDir(const Glib::ustring& dirName) const
    {
        std::vector<Glib::ustring> res;
        ProfileMap profiles;

        MyMutex::MyLock lock(mutex);

        loadProfiles(profilesDir, &profiles, nullptr, nullptr, false);
        loadProfiles(dirName, &profiles, nullptr, nullptr, false);

        for (const auto& profile : profiles) {
            res.push_back(profile.first);
        }

        return res;
    }

    std::uint8_t getInputIntents(cmsHPROFILE profile)
    {
        MyMutex::MyLock lock(mutex);

        return getSupportedIntents(profile, LCMS_USED_AS_INPUT);
    }

    std::uint8_t getOutputIntents(cmsHPROFILE profile)
    {
        MyMutex::MyLock lock(mutex);

        return getSupportedIntents(profile, LCMS_USED_AS_OUTPUT);
    }

    std::uint8_t getProofIntents(cmsHPROFILE profile)
    {
        MyMutex::MyLock lock(mutex);

        return getSupportedIntents(profile, LCMS_USED_AS_PROOF);
    }

    std::uint8_t getInputIntents(const Glib::ustring &name)
    {
        return getInputIntents(getProfile(name));
    }

    std::uint8_t getOutputIntents(const Glib::ustring &name)
    {
        return getOutputIntents(getProfile(name));
    }

    std::uint8_t getProofIntents(const Glib::ustring &name)
    {
        return getProofIntents(getProfile(name));
    }

    Glib::ustring getDefaultMonitorProfileName() const
    {
        return defaultMonitorProfile;
    }

    void setDefaultMonitorProfileName(const Glib::ustring &name)
    {
        MyMutex::MyLock lock(mutex);
        defaultMonitorProfile = name;
    }

    std::vector<Glib::ustring> getWorkingProfiles()
    {
        std::vector<Glib::ustring> res;

        // for (unsigned int i = 0; i < sizeof(wpnames) / sizeof(wpnames[0]); i++) {
        //     res.push_back(wpnames[i]);
        // }
        for (const auto &p : wProfiles) {
            res.push_back(p.first);
        }

        return res;
    }

private:
    using CVector = std::array<double, 3>;
    using CMatrix = std::array<CVector, 3>;
    struct PMatrix {
        double matrix[3][3];
        PMatrix(): matrix{} {}
        PMatrix(const CMatrix &m)
        {
            set(m);
        }

        CMatrix toMatrix() const
        {
            CMatrix ret;

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    ret[i][j] = matrix[i][j];
                }
            }

            return ret;
        }

        void set(const CMatrix &m)
        {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    matrix[i][j] = m[i][j];
                }
            }
        }
    };

    bool computeWorkingSpaceMatrix(const Glib::ustring &path, const Glib::ustring &filename, PMatrix &out)
    {
        Glib::ustring fullpath = filename;

        if (!Glib::path_is_absolute(fullpath)) {
            fullpath = Glib::build_filename(path, filename);
        }

        ProfileContent content(fullpath);
        cmsHPROFILE prof = content.toProfile();

        if (!prof) {
            return false;
        }

        if (cmsGetColorSpace(prof) != cmsSigRgbData) {
            cmsCloseProfile(prof);
            return false;
        }

        if (!cmsIsMatrixShaper(prof)) {
            cmsCloseProfile(prof);
            return false;
        }

        cmsCIEXYZ *red = static_cast<cmsCIEXYZ *>(cmsReadTag(prof, cmsSigRedMatrixColumnTag));
        cmsCIEXYZ *green  = static_cast<cmsCIEXYZ *>(cmsReadTag(prof, cmsSigGreenMatrixColumnTag));
        cmsCIEXYZ *blue  = static_cast<cmsCIEXYZ *>(cmsReadTag(prof, cmsSigBlueMatrixColumnTag));

        if (!red || !green || !blue) {
            cmsCloseProfile(prof);
            return false;
        }

        CMatrix m = {
            CVector({ red->X, green->X, blue->X }),
            CVector({ red->Y, green->Y, blue->Y }),
            CVector({ red->Z, green->Z, blue->Z })
        };
        m[1][0] = red->Y;
        m[1][1] = green->Y;
        m[1][2] = blue->Y;
        m[2][0] = red->Z;
        m[2][1] = green->Z;
        m[2][2] = blue->Z;
        out.set(m);

        cmsCloseProfile(prof);
        return true;
    }

    bool loadWorkingSpaces(const Glib::ustring &path)
    {
        Glib::ustring fileName = Glib::build_filename(path, "workingspaces.json");
        FILE* const f = g_fopen(fileName.c_str(), "r");

        if (settings->verbose) {
            std::cout << "trying to load extra working spaces from " << fileName << std::flush;
        }

        if (!f) {
            if (settings->verbose) {
                std::cout << " FAIL" << std::endl;
            }

            return false;
        }

        fseek(f, 0, SEEK_END);
        long length = ftell(f);

        if (length <= 0) {
            if (settings->verbose) {
                std::cout << " FAIL" << std::endl;
            }

            fclose(f);
            return false;
        }

        char *buf = new char[length + 1];
        fseek(f, 0, SEEK_SET);
        length = fread(buf, 1, length, f);
        buf[length] = 0;

        fclose(f);

        cJSON_Minify(buf);
        cJSON *root = cJSON_Parse(buf);

        if (!root) {
            if (settings->verbose) {
                std::cout << " FAIL" << std::endl;
            }

            return false;
        }

        delete[] buf;

        cJSON *js = cJSON_GetObjectItem(root, "working_spaces");

        if (!js) {
            goto parse_error;
        }

        for (js = js->child; js != nullptr; js = js->next) {
            cJSON *ji = cJSON_GetObjectItem(js, "name");
            std::unique_ptr<PMatrix> m(new PMatrix);
            std::string name;

            if (!ji || ji->type != cJSON_String) {
                goto parse_error;
            }

            name = ji->valuestring;

            if (wProfiles.find(name) != wProfiles.end()) {
                continue; // already there -- ignore
            }

            bool found_matrix = false;

            ji = cJSON_GetObjectItem(js, "matrix");

            if (ji) {
                if (ji->type != cJSON_Array) {
                    goto parse_error;
                }

                ji = ji->child;

                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j, ji = ji->next) {
                        if (!ji || ji->type != cJSON_Number) {
                            goto parse_error;
                        }

                        m->matrix[i][j] = ji->valuedouble;
                    }
                }

                if (ji) {
                    goto parse_error;
                }

                found_matrix = true;
            } else {
                ji = cJSON_GetObjectItem(js, "file");

                if (!ji || ji->type != cJSON_String) {
                    goto parse_error;
                }

                found_matrix = computeWorkingSpaceMatrix(path, ji->valuestring, *m);
            }

            if (!found_matrix) {
                if (settings->verbose) {
                    std::cout << "Could not find suitable matrix for working space: " << name << std::endl;
                }

                continue;
            }

            pMatrices.emplace_back(std::move(m));
            TMatrix w = pMatrices.back()->matrix;

            CMatrix b = {};

            if (!rtengine::invertMatrix(pMatrices.back()->toMatrix(), b)) {
                if (settings->verbose) {
                    std::cout << "Matrix for working space: " << name << " is not invertible, skipping" << std::endl;
                }

                pMatrices.pop_back();
            } else {
                wMatrices[name] = w;
                pMatrices.emplace_back(new PMatrix(b));
                TMatrix iw = pMatrices.back()->matrix;
                iwMatrices[name] = iw;
                wProfiles[name] = createFromMatrix(w);

                if (settings->verbose) {
                    std::cout << "Added working space: " << name << std::endl;
                    std::cout << "  matrix: [";

                    for (int i = 0; i < 3; ++i) {
                        std::cout << " [";

                        for (int j = 0; j < 3; ++j) {
                            std::cout << " " << w[i][j];
                        }

                        std::cout << "]";
                    }

                    std::cout << " ]" << std::endl;
                }
            }
        }

        cJSON_Delete(root);

        if (settings->verbose) {
            std::cout << " OK" << std::endl;
        }

        return true;

parse_error:

        if (settings->verbose) {
            std::cout << " ERROR in parsing " << fileName << std::endl;
        }

        cJSON_Delete(root);
        return false;
    }

    using ProfileMap = std::map<Glib::ustring, cmsHPROFILE>;
    using MatrixMap = std::map<Glib::ustring, TMatrix>;
    using ContentMap = std::map<Glib::ustring, ProfileContent>;
    using NameMap = std::map<Glib::ustring, Glib::ustring>;

    ProfileMap wProfiles;
    // ProfileMap wProfilesGamma;
    MatrixMap wMatrices;
    MatrixMap iwMatrices;

    std::vector<std::unique_ptr<PMatrix>> pMatrices;

    // These contain profiles from user/system directory(supplied on init)
    Glib::ustring profilesDir;
    Glib::ustring userICCDir;
    ProfileMap fileProfiles;
    ContentMap fileProfileContents;

    //These contain standard profiles from RT. Keys are all in uppercase.
    Glib::ustring stdProfilesDir;
    NameMap fileStdProfilesFileNames;
    ProfileMap fileStdProfiles;

    Glib::ustring defaultMonitorProfile;

    bool loadAll;

    const cmsHPROFILE xyz;
    const cmsHPROFILE srgb;

    mutable MyMutex mutex;
};

rtengine::ICCStore* rtengine::ICCStore::getInstance()
{
    static rtengine::ICCStore instance;
    return &instance;
}

void rtengine::ICCStore::init(const Glib::ustring& usrICCDir, const Glib::ustring& stdICCDir, bool loadAll)
{
    implementation->init(usrICCDir, stdICCDir, loadAll);
}

cmsHPROFILE rtengine::ICCStore::workingSpace(const Glib::ustring& name) const
{
    return implementation->workingSpace(name);
}

// cmsHPROFILE rtengine::ICCStore::workingSpaceGamma(const Glib::ustring& name) const
// {
//     return implementation->workingSpaceGamma(name);
// }

rtengine::TMatrix rtengine::ICCStore::workingSpaceMatrix(const Glib::ustring& name) const
{
    return implementation->workingSpaceMatrix(name);
}

rtengine::TMatrix rtengine::ICCStore::workingSpaceInverseMatrix(const Glib::ustring& name) const
{
    return implementation->workingSpaceInverseMatrix(name);
}

bool rtengine::ICCStore::outputProfileExist(const Glib::ustring& name) const
{
    return implementation->outputProfileExist(name);
}

cmsHPROFILE rtengine::ICCStore::getProfile(const Glib::ustring& name) const
{
    return implementation->getProfile(name);
}

cmsHPROFILE rtengine::ICCStore::getStdProfile(const Glib::ustring& name) const
{
    return implementation->getStdProfile(name);
}

rtengine::ProfileContent rtengine::ICCStore::getContent(const Glib::ustring& name) const
{
    return implementation->getContent(name);
}


Glib::ustring rtengine::ICCStore::getDefaultMonitorProfileName() const
{
    return implementation->getDefaultMonitorProfileName();
}


void rtengine::ICCStore::setDefaultMonitorProfileName(const Glib::ustring &name)
{
    implementation->setDefaultMonitorProfileName(name);
}

cmsHPROFILE rtengine::ICCStore::getXYZProfile() const
{
    return implementation->getXYZProfile();
}

cmsHPROFILE rtengine::ICCStore::getsRGBProfile() const
{
    return implementation->getsRGBProfile();
}

std::vector<Glib::ustring> rtengine::ICCStore::getProfiles(ProfileType type) const
{
    return implementation->getProfiles(type);
}

std::vector<Glib::ustring> rtengine::ICCStore::getProfilesFromDir(const Glib::ustring& dirName) const
{
    return implementation->getProfilesFromDir(dirName);
}

std::uint8_t rtengine::ICCStore::getInputIntents(cmsHPROFILE profile) const
{
    return implementation->getInputIntents(profile);
}

std::uint8_t rtengine::ICCStore::getOutputIntents(cmsHPROFILE profile) const
{
    return implementation->getOutputIntents(profile);
}

std::uint8_t rtengine::ICCStore::getProofIntents(cmsHPROFILE profile) const
{
    return implementation->getProofIntents(profile);
}

std::uint8_t rtengine::ICCStore::getInputIntents(const Glib::ustring& name) const
{
    return implementation->getInputIntents(name);
}

std::uint8_t rtengine::ICCStore::getOutputIntents(const Glib::ustring& name) const
{
    return implementation->getOutputIntents(name);
}

std::uint8_t rtengine::ICCStore::getProofIntents(const Glib::ustring& name) const
{
    return implementation->getProofIntents(name);
}

rtengine::ICCStore::ICCStore() :
    implementation(new Implementation)
{
}

rtengine::ICCStore::~ICCStore() = default;

std::vector<Glib::ustring> rtengine::ICCStore::getWorkingProfiles()
{
    return implementation->getWorkingProfiles();
}

std::vector<Glib::ustring> rtengine::ICCStore::getGamma()
{

    std::vector<Glib::ustring> res;

    for (unsigned int i = 0; i < sizeof(wpgamma) / sizeof(wpgamma[0]); i++) {
        res.push_back(wpgamma[i]);
    }

    return res;
}

void rtengine::ICCStore::getGammaArray(const procparams::ColorManagementParams &icm, GammaValues &ga)
{
    const double eps = 0.000000001; // not divide by zero

    if (icm.freegamma  && icm.gamma != "Free") { //if Free gamma selected with other than Free
        // gamma : ga[0],ga[1],ga[2],ga[3],ga[4],ga[5] by calcul
        if (icm.gamma == "BT709_g2.2_s4.5")      {
            ga[0] = 2.22;    //BT709  2.2  4.5  - my preferred as D.Coffin
            ga[1] = 0.909995;
            ga[2] = 0.090005;
            ga[3] = 0.222222;
            ga[4] = 0.081071;
        } else if (icm.gamma == "sRGB_g2.4_s12.92")   {
            ga[0] = 2.40;    //sRGB 2.4 12.92  - RT default as Lightroom
            ga[1] = 0.947858;
            ga[2] = 0.052142;
            ga[3] = 0.077399;
            ga[4] = 0.039293;
        } else if (icm.gamma == "High_g1.3_s3.35")    {
            ga[0] = 1.3 ;    //for high dynamic images
            ga[1] = 0.998279;
            ga[2] = 0.001721;
            ga[3] = 0.298507;
            ga[4] = 0.005746;
        } else if (icm.gamma == "Low_g2.6_s6.9")   {
            ga[0] = 2.6 ;    //gamma 2.6 variable : for low contrast images
            ga[1] = 0.891161;
            ga[2] = 0.108839;
            ga[3] = 0.144928;
            ga[4] = 0.076332;
        } else if (icm.gamma == "standard_g2.2")   {
            ga[0] = 2.2;    //gamma=2.2(as gamma of Adobe, Widegamut...)
            ga[1] = 1.;
            ga[2] = 0.;
            ga[3] = 1. / eps;
            ga[4] = 0.;
        } else if (icm.gamma == "standard_g1.8")   {
            ga[0] = 1.8;    //gamma=1.8(as gamma of Prophoto)
            ga[1] = 1.;
            ga[2] = 0.;
            ga[3] = 1. / eps;
            ga[4] = 0.;
        } else if (icm.gamma == "Lab_g3.0s9.03296")   {
            ga[0] = 3.0;    //Lab gamma =3 slope=9.03296
            ga[1] = 0.8621;
            ga[2] = 0.1379;
            ga[3] = 0.1107;
            ga[4] = 0.08;

        } else { /* if (icm.gamma == "linear_g1.0") */
            ga[0] = 1.0;    //gamma=1 linear : for high dynamic images(cf : D.Coffin...)
            ga[1] = 1.;
            ga[2] = 0.;
            ga[3] = 1. / eps;
            ga[4] = 0.;
        }

        ga[5] = 0.0;
        ga[6] = 0.0;
    } else { //free gamma selected
        GammaValues g_a; //gamma parameters
        double pwr = 1.0 / icm.gampos;
        double ts = icm.slpos;
        double slope = icm.slpos == 0 ? eps : icm.slpos;

        int mode = 0;
        Color::calcGamma(pwr, ts, mode, g_a); // call to calcGamma with selected gamma and slope : return parameters for LCMS2
        ga[4] = g_a[3] * ts;
        //printf("g_a.gamma0=%f g_a.gamma1=%f g_a.gamma2=%f g_a.gamma3=%f g_a.gamma4=%f\n", g_a.gamma0,g_a.gamma1,g_a.gamma2,g_a.gamma3,g_a.gamma4);
        ga[0] = icm.gampos;
        ga[1] = 1. / (1.0 + g_a[4]);
        ga[2] = g_a[4] / (1.0 + g_a[4]);
        ga[3] = 1. / slope;
        ga[5] = 0.0;
        ga[6] = 0.0;
        //   printf("ga[0]=%f ga[1]=%f ga[2]=%f ga[3]=%f ga[4]=%f\n", ga[0],ga[1],ga[2],ga[3],ga[4]);
    }
}

// WARNING: the caller must lock lcmsMutex
cmsHPROFILE rtengine::ICCStore::makeStdGammaProfile(cmsHPROFILE iprof)
{
    // forgive me for the messy code, quick hack to change gamma of an ICC profile to the RT standard gamma
    if (!iprof) {
        return nullptr;
    }

    cmsUInt32Number bytesNeeded = 0;
    cmsSaveProfileToMem(iprof, nullptr, &bytesNeeded);

    if (bytesNeeded == 0) {
        return nullptr;
    }

    uint8_t *data = new uint8_t[bytesNeeded + 1];
    cmsSaveProfileToMem(iprof, data, &bytesNeeded);
    const uint8_t *p = &data[128]; // skip 128 byte header
    uint32_t tag_count;
    memcpy(&tag_count, p, 4);
    p += 4;
    tag_count = ntohl(tag_count);

    struct icctag {
        uint32_t sig;
        uint32_t offset;
        uint32_t size;
    } tags[tag_count];

    const uint32_t gamma = 0x239;
    int gamma_size = 14;
    int data_size = (gamma_size + 3) & ~3;

    for (uint32_t i = 0; i < tag_count; i++) {
        memcpy(&tags[i], p, 12);
        tags[i].sig = ntohl(tags[i].sig);
        tags[i].offset = ntohl(tags[i].offset);
        tags[i].size = ntohl(tags[i].size);
        p += 12;

        if (tags[i].sig != 0x62545243 && // bTRC
                tags[i].sig != 0x67545243 && // gTRC
                tags[i].sig != 0x72545243 && // rTRC
                tags[i].sig != 0x6B545243) { // kTRC
            data_size += (tags[i].size + 3) & ~3;
        }
    }

    uint32_t sz = 128 + 4 + tag_count * 12 + data_size;
    uint8_t *nd = new uint8_t[sz];
    memset(nd, 0, sz);
    memcpy(nd, data, 128 + 4);
    sz = htonl(sz);
    memcpy(nd, &sz, 4);
    uint32_t offset = 128 + 4 + tag_count * 12;
    uint32_t gamma_offset = 0;

    for (uint32_t i = 0; i < tag_count; i++) {
        struct icctag tag;
        tag.sig = htonl(tags[i].sig);

        if (tags[i].sig == 0x62545243 || // bTRC
                tags[i].sig == 0x67545243 || // gTRC
                tags[i].sig == 0x72545243 || // rTRC
                tags[i].sig == 0x6B545243) { // kTRC
            if (gamma_offset == 0) {
                gamma_offset = offset;
                uint32_t pcurve[] = { htonl(0x63757276), htonl(0), htonl(gamma_size == 12 ? 0U : 1U) };
                memcpy(&nd[offset], pcurve, 12);

                if (gamma_size == 14) {
                    uint16_t gm = htons(gamma);
                    memcpy(&nd[offset + 12], &gm, 2);
                }

                offset += (gamma_size + 3) & ~3;
            }

            tag.offset = htonl(gamma_offset);
            tag.size = htonl(gamma_size);
        } else {
            tag.offset = htonl(offset);
            tag.size = htonl(tags[i].size);
            memcpy(&nd[offset], &data[tags[i].offset], tags[i].size);
            offset += (tags[i].size + 3) & ~3;
        }

        memcpy(&nd[128 + 4 + i * 12], &tag, 12);
    }

    cmsHPROFILE oprof = cmsOpenProfileFromMem(nd, ntohl(sz));
    delete [] nd;
    delete [] data;
    return oprof;
}

cmsHPROFILE rtengine::ICCStore::createFromMatrix(const double matrix[3][3], bool gamma, const Glib::ustring& name)
{

    static const unsigned phead[] = {
        1024, 0, 0x2100000, 0x6d6e7472, 0x52474220, 0x58595a20, 0, 0, 0,
        0x61637370, 0, 0, 0, 0, 0, 0, 0, 0xf6d6, 0x10000, 0xd32d
    };
    unsigned pbody[] = {
        10, 0x63707274, 0, 36,  /* cprt */
        0x64657363, 0, 40,  /* desc */
        0x77747074, 0, 20,  /* wtpt */
        0x626b7074, 0, 20,  /* bkpt */
        0x72545243, 0, 14,  /* rTRC */
        0x67545243, 0, 14,  /* gTRC */
        0x62545243, 0, 14,  /* bTRC */
        0x7258595a, 0, 20,  /* rXYZ */
        0x6758595a, 0, 20,  /* gXYZ */
        0x6258595a, 0, 20
    };    /* bXYZ */
    static const unsigned pwhite[] = { 0xf351, 0x10000, 0x116cc };//D65
    //static const unsigned pwhite[] = { 0xf6d6, 0x10000, 0xd340 };//D50

    // 0x63757276 : curveType, 0 : reserved, 1 : entries(1=gamma, 0=identity), 0x1000000=1.0
    unsigned pcurve[] = { 0x63757276, 0, 0, 0x1000000 };
//    unsigned pcurve[] = { 0x63757276, 0, 1, 0x1000000 };

    if (gamma) {
        pcurve[2] = 1;
        // pcurve[3] = 0x1f00000;// pcurve for gamma BT709 : g=2.22 s=4.5
        // normalize gamma in RT, default(Emil's choice = sRGB)
        pcurve[3] = 0x2390000;//pcurve for gamma sRGB : g:2.4 s=12.92

    } else {
        // lcms2 up to 2.4 has a bug with linear gamma causing precision loss(banding)
        // of floating point data when a normal icc encoding of linear gamma is used
        //(i e 0 table entries), but by encoding a gamma curve which is 1.0 the
        // floating point path is taken within lcms2 so no precision loss occurs and
        // gamma is still 1.0.
        pcurve[2] = 1;
        pcurve[3] = 0x1000000; //pcurve for gamma 1
    }

    // constructing profile header
    unsigned* oprof = new unsigned [phead[0] / sizeof(unsigned)];
    memset(oprof, 0, phead[0]);
    memcpy(oprof, phead, sizeof(phead));

    oprof[0] = 132 + 12 * pbody[0];

    // constructing tag directory(pointers inside the file), and types
    // 0x74657874 : text
    // 0x64657363 : description tag
    for (unsigned int i = 0; i < pbody[0]; i++) {
        oprof[oprof[0] / 4] = i ? (i > 1 ? 0x58595a20 : 0x64657363) : 0x74657874;
        pbody[i * 3 + 2] = oprof[0];
        oprof[0] += (pbody[i * 3 + 3] + 3) & -4;
    }

    memcpy(oprof + 32, pbody, sizeof(pbody));

    // wtpt
    memcpy((char *)oprof + pbody[8] + 8, pwhite, sizeof(pwhite));

    // r/g/b TRC
    for (int i = 4; i < 7; i++) {
        memcpy((char *)oprof + pbody[i * 3 + 2], pcurve, sizeof(pcurve));
    }

    // r/g/b XYZ
//    pseudoinverse((double(*)[3]) out_rgb[output_color-1], inverse, 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            oprof[pbody[j * 3 + 23] / 4 + i + 2] = matrix[i][j] * 0x10000 + 0.5;
//      for (num = k=0; k < 3; k++)
//        num += xyzd50_srgb[i][k] * inverse[j][k];
        }

    // convert to network byte order
    for (unsigned int i = 0; i < phead[0] / 4; i++) {
        oprof[i] = htonl(oprof[i]);
    }

    // cprt
    strcpy((char *)oprof + pbody[2] + 8, "--rawtherapee profile--");

    // desc
    oprof[pbody[5] / 4 + 2] = name.size() + 1;
    strcpy((char *)oprof + pbody[5] + 12, name.c_str());


    cmsHPROFILE p = cmsOpenProfileFromMem(oprof, ntohl(oprof[0]));
    delete [] oprof;
    return p;
}

cmsHPROFILE rtengine::ICCStore::createGammaProfile(const procparams::ColorManagementParams &icm, GammaValues &ga)
{
    float p[6]; //primaries
    ga[6] = 0.0;

    enum class ColorTemp {
        D50 = 5003,  // for Widegamut, Prophoto Best, Beta -> D50
        D65 = 6504,   // for sRGB, AdobeRGB, Bruce Rec2020  -> D65
        D60 = 6005        //for ACESP0 and AcesP1

    };
    ColorTemp temp = ColorTemp::D50;

    //primaries for 10 working profiles ==> output profiles
    if (icm.wprimari == "wideg") {
        p[0] = 0.7350;    //Widegamut primaries
        p[1] = 0.2650;
        p[2] = 0.1150;
        p[3] = 0.8260;
        p[4] = 0.1570;
        p[5] = 0.0180;
    } else if (icm.wprimari == "adob") {
        p[0] = 0.6400;    //Adobe primaries
        p[1] = 0.3300;
        p[2] = 0.2100;
        p[3] = 0.7100;
        p[4] = 0.1500;
        p[5] = 0.0600;
        temp = ColorTemp::D65;
    } else if (icm.wprimari == "srgb") {
        p[0] = 0.6400;    // sRGB primaries
        p[1] = 0.3300;
        p[2] = 0.3000;
        p[3] = 0.6000;
        p[4] = 0.1500;
        p[5] = 0.0600;
        temp = ColorTemp::D65;
    } else if (icm.wprimari == "BruceRGB") {
        p[0] = 0.6400;    // Bruce primaries
        p[1] = 0.3300;
        p[2] = 0.2800;
        p[3] = 0.6500;
        p[4] = 0.1500;
        p[5] = 0.0600;
        temp = ColorTemp::D65;
    } else if (icm.wprimari == "BetaRGB") {
        p[0] = 0.6888;    // Beta primaries
        p[1] = 0.3112;
        p[2] = 0.1986;
        p[3] = 0.7551;
        p[4] = 0.1265;
        p[5] = 0.0352;
    } else if (icm.wprimari == "BestRGB") {
        p[0] = 0.7347;    // Best primaries
        p[1] = 0.2653;
        p[2] = 0.2150;
        p[3] = 0.7750;
        p[4] = 0.1300;
        p[5] = 0.0350;
    } else if (icm.wprimari == "rec2020") {
        p[0] = 0.7080;    // Rec2020 primaries
        p[1] = 0.2920;
        p[2] = 0.1700;
        p[3] = 0.7970;
        p[4] = 0.1310;
        p[5] = 0.0460;
        temp = ColorTemp::D65;
    } else if (icm.wprimari == "acesp0") {
        p[0] = 0.7347;    // ACES P0 primaries
        p[1] = 0.2653;
        p[2] = 0.0000;
        p[3] = 1.0;
        p[4] = 0.0001;
        p[5] = -0.0770;
        temp = ColorTemp::D60;
    } else if (icm.wprimari == "acesp1") {
        p[0] = 0.713;    // ACES P1 primaries
        p[1] = 0.293;
        p[2] = 0.165;
        p[3] = 0.830;
        p[4] = 0.128;
        p[5] = 0.044;
        temp = ColorTemp::D60;
    } else if (icm.wprimari == "proph") {
        p[0] = 0.7347;    //ProPhoto and default primaries
        p[1] = 0.2653;
        p[2] = 0.1596;
        p[3] = 0.8404;
        p[4] = 0.0366;
        p[5] = 0.0001;
    } else {
        p[0] = 0.7347;    //default primaries
        p[1] = 0.2653;
        p[2] = 0.1596;
        p[3] = 0.8404;
        p[4] = 0.0366;
        p[5] = 0.0001;
    }


    cmsCIExyY xyD;
    cmsCIExyYTRIPLE Primaries = {
        {p[0], p[1], 1.0}, // red
        {p[2], p[3], 1.0}, // green
        {p[4], p[5], 1.0}  // blue
    };
    cmsToneCurve* GammaTRC[3];

    // 7 parameters for smoother curves
    cmsFloat64Number Parameters[7] = { ga[0],  ga[1], ga[2], ga[3], ga[4], ga[5], ga[6] } ;

    //lcmsMutex->lock();  Mutex acquired by the caller
    cmsWhitePointFromTemp(&xyD, (double)temp);

    GammaTRC[0] = GammaTRC[1] = GammaTRC[2] = cmsBuildParametricToneCurve(nullptr, 5, Parameters); //5 = smoother than 4
    cmsHPROFILE oprofdef = cmsCreateRGBProfile(&xyD, &Primaries, GammaTRC); //oprofdef  become Outputprofile

    cmsFreeToneCurve(GammaTRC[0]);
    //lcmsMutex->unlock();

    return oprofdef;
}

// WARNING: the caller must lock lcmsMutex
cmsHPROFILE rtengine::ICCStore::createCustomGammaOutputProfile(const procparams::ColorManagementParams &icm, GammaValues &ga)
{
    bool pro = false;
    Glib::ustring outProfile;
    cmsHPROFILE outputProfile = nullptr;
    Glib::ustring outPr;
    Glib::ustring gammaStr;


    if (icm.freegamma && icm.gampos < 1.35) {
        pro = true;    //select profil with gammaTRC modified :
    } else if (icm.gamma == "linear_g1.0" || (icm.gamma == "High_g1.3_s3.35")) {
        pro = true;    //pro=0  RT_sRGB || Prophoto
    }

    //necessary for V2 profile
    if (icm.wprimari == "proph"    && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.prophoto)   && !pro) {
        outProfile = options.rtSettings.prophoto;
        outPr = "RT_large";

    } else if (icm.wprimari == "adob" && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.adobe)) {
        outProfile = options.rtSettings.adobe;
        outPr = "RT_adob";
    } else if (icm.wprimari == "wideg" && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.widegamut)) {
        outProfile = options.rtSettings.widegamut;
        outPr = "RT_wide";
    } else if (icm.wprimari == "BetaRGB"  && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.beta)) {
        outProfile = options.rtSettings.beta;
        outPr = "RT_beta";
    } else if (icm.wprimari == "BestRGB"   && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.best)) {
        outProfile = options.rtSettings.best;
        outPr = "RT_best";
    } else if (icm.wprimari == "BruceRGB"  && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.bruce)) {
        outProfile = options.rtSettings.bruce;
        outPr = "RT_bruce";
    } else if (icm.wprimari == "srgb"      && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.srgb)       && !pro) {
        outProfile = options.rtSettings.srgb;
        outPr = "RT_srgb";
    } else if (icm.wprimari == "sRGB"      && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.srgb10)     &&  pro) {
        outProfile = options.rtSettings.srgb10;
        outPr = "RT_srgb";
    } else if (icm.wprimari == "proph"  && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.prophoto10) &&  pro) {
        outProfile = options.rtSettings.prophoto10;
        outPr = "RT_large";
    } else if (icm.wprimari == "rec2020"   && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.rec2020)) {
        outProfile = options.rtSettings.rec2020;
        outPr = "RT_rec2020";
    } else if (icm.wprimari == "acesp0"   && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.ACESp0)) {
        outProfile = options.rtSettings.ACESp0;
        outPr = "RT_acesp0";
    } else if (icm.wprimari == "acesp1"   && rtengine::ICCStore::getInstance()->outputProfileExist(options.rtSettings.ACESp1)) {
        outProfile = options.rtSettings.ACESp1;
        outPr = "RT_acesp1";
    } else {
        // Should not occurs
        if (settings->verbose) {
            printf("\"%s\": unknown working profile! - use LCMS2 substitution\n", icm.working.c_str());
        }

        return nullptr;
    }


    //begin adaptation rTRC gTRC bTRC
    //"outputProfile" profile has the same characteristics than RGB values, but TRC are adapted... for applying profile
    if (settings->verbose) {
        printf("Output Gamma - profile Primaries as RT profile: \"%s\"\n", outProfile.c_str());      //c_str()
    }

    outputProfile = ICCStore::getInstance()->getProfile(outProfile); //get output profile

    if (outputProfile == nullptr) {

        if (settings->verbose) {
            printf("\"%s\" ICC output profile not found!\n", outProfile.c_str());
        }

        return nullptr;
    }

    // 7 parameters for smoother curves
    cmsFloat64Number Parameters[7] = { ga[0], ga[1], ga[2], ga[3], ga[4], ga[5], ga[6] };

    //change desc Tag , to "free gamma", or "BT709", etc.
    cmsMLU *mlu;
    cmsContext ContextID = cmsGetProfileContextID(outputProfile); // create context to modify some TAGs
    mlu = cmsMLUalloc(ContextID, 1);
    Glib::ustring outPro;

    if (icm.gamma == "High_g1.3_s3.35") {
        gammaStr = "_High_g=1.3_s=3.35";
    } else if (icm.gamma == "Low_g2.6_s6.9") {
        gammaStr = "_Low_g=2.6_s=6.9";
    } else if (icm.gamma == "sRGB_g2.4_s12.92") {
        gammaStr = "_sRGB_g=2.4_s=12.92";
    } else if (icm.gamma == "BT709_g2.2_s4.5") {
        gammaStr = "_BT709_g=2.2_s=4.5";
    } else if (icm.gamma == "linear_g1.0") {
        gammaStr = "_Linear_g=1.0";
    } else if (icm.gamma == "standard_g2.2") {
        gammaStr = "_g=2.2";
    } else if (icm.gamma == "standard_g1.8") {
        gammaStr = "_g=1.8";
    } else if (icm.gamma == "Lab_g3.0s9.03296") {
        gammaStr = "_LAB_g3.0_s9.03296";
    }

    // create description with gamma + slope + primaries
    std::wostringstream gammaWs;
    std::wstring gammaStrICC;

    gammaWs.precision(3);

    if (icm.gamma == "Free") {
        if (icm.wprofile == "v4") {
            outPro = outPr + "_V4_" + std::to_string((float)icm.gampos) + " " + std::to_string((float)icm.slpos) + ".icc";
        } else if (icm.wprofile == "v2" || icm.wprofile == "none" ) { 
            outPro = outPr + "_V2_" + std::to_string((float)icm.gampos) + " " + std::to_string((float)icm.slpos) + ".icc";
        }

        gammaWs.precision(2);
        gammaWs << outPr << (float)icm.gampos << " s=" << (float)icm.slpos;


        cmsMLUsetWide(mlu,  "en", "US", gammaWs.str().c_str());

    } else {

        if (icm.wprofile == "v4") {
            outPro = outPr + "_V4_" + gammaStr + ".icc";

        } else if (icm.wprofile == "v2" || icm.wprofile == "none" ) {
            outPro = outPr + "_V2_" + gammaStr + ".icc";
        }

        gammaWs << outPro.c_str()  << " s=";


    }



    cmsMLUsetWide(mlu,  "en", "US", gammaWs.str().c_str());
    cmsMLU *copyright = cmsMLUalloc(NULL, 1);

    cmsMLUsetASCII(copyright, "en", "US", "No copyright Rawtherapee");
    cmsWriteTag(outputProfile, cmsSigCopyrightTag, copyright);
    cmsMLUfree(copyright);

    cmsWriteTag(outputProfile, cmsSigProfileDescriptionTag,  mlu);//desc changed

    cmsMLU *descrip = cmsMLUalloc(NULL, 1);

    cmsMLUsetASCII(descrip, "en", "US", "Rawtherapee");
    cmsWriteTag(outputProfile, cmsSigDeviceModelDescTag, descrip);
    cmsMLUfree(descrip);



    // instruction with //ICC are used to generate ICC profile
    if (mlu == nullptr) {
        printf("Description error\n");
    } else {

        if (icm.wprofile == "v4") {
            cmsSetProfileVersion(outputProfile, 4.3);
        } else {
            cmsSetProfileVersion(outputProfile, 2.0);
        }

//change
        float p[6]; //primaries
        ga[6] = 0.0;

        enum class ColorTemp {
            D50 = 5003,  // for Widegamut, Prophoto Best, Beta -> D50
            D65 = 6504,   // for sRGB, AdobeRGB, Bruce Rec2020  -> D65
            D60 = 6005        //for ACESc->D60
        };
        ColorTemp temp = ColorTemp::D50;

        if (icm.wprimari == "wideg") {
            p[0] = 0.7350;    //Widegamut primaries
            p[1] = 0.2650;
            p[2] = 0.1150;
            p[3] = 0.8260;
            p[4] = 0.1570;
            p[5] = 0.0180;

        } else if (icm.wprimari == "adob") {
            p[0] = 0.6400;    //Adobe primaries
            p[1] = 0.3300;
            p[2] = 0.2100;
            p[3] = 0.7100;
            p[4] = 0.1500;
            p[5] = 0.0600;
            temp = ColorTemp::D65;
        } else if (icm.wprimari == "srgb") {
            p[0] = 0.6400;    // sRGB primaries
            p[1] = 0.3300;
            p[2] = 0.3000;
            p[3] = 0.6000;
            p[4] = 0.1500;
            p[5] = 0.0600;
            temp = ColorTemp::D65;
        } else if (icm.wprimari == "BruceRGB") {
            p[0] = 0.6400;    // Bruce primaries
            p[1] = 0.3300;
            p[2] = 0.2800;
            p[3] = 0.6500;
            p[4] = 0.1500;
            p[5] = 0.0600;
            temp = ColorTemp::D65;
        } else if (icm.wprimari == "BetaRGB") {
            p[0] = 0.6888;    // Beta primaries
            p[1] = 0.3112;
            p[2] = 0.1986;
            p[3] = 0.7551;
            p[4] = 0.1265;
            p[5] = 0.0352;
        } else if (icm.wprimari == "BestRGB") {
            p[0] = 0.7347;    // Best primaries
            p[1] = 0.2653;
            p[2] = 0.2150;
            p[3] = 0.7750;
            p[4] = 0.1300;
            p[5] = 0.0350;
        } else if (icm.wprimari == "rec2020") {
            p[0] = 0.7080;    // Rec2020 primaries
            p[1] = 0.2920;
            p[2] = 0.1700;
            p[3] = 0.7970;
            p[4] = 0.1310;
            p[5] = 0.0460;
            temp = ColorTemp::D65;
        } else if (icm.wprimari == "acesp0") {
            p[0] = 0.7347;    // ACES P0 primaries
            p[1] = 0.2653;
            p[2] = 0.0000;
            p[3] = 1.0;
            p[4] = 0.0001;
            p[5] = -0.0770;
            temp = ColorTemp::D60;
        } else if (icm.wprimari == "acesp1") {
            p[0] = 0.713;    // ACES P1 primaries
            p[1] = 0.293;
            p[2] = 0.165;
            p[3] = 0.830;
            p[4] = 0.128;
            p[5] = 0.044;
            temp = ColorTemp::D60;
        } else if (icm.wprimari == "proph") {
            p[0] = 0.7347;    //ProPhoto and default primaries
            p[1] = 0.2653;
            p[2] = 0.1596;
            p[3] = 0.8404;
            p[4] = 0.0366;
            p[5] = 0.0001;
        } else {
            p[0] = 0.7347;    //default primaries
            p[1] = 0.2653;
            p[2] = 0.1596;
            p[3] = 0.8404;
            p[4] = 0.0366;
            p[5] = 0.0001;
        }

        cmsCIExyY xyD;
        cmsCIExyYTRIPLE Primaries = {
            {p[0], p[1], 1.0}, // red
            {p[2], p[3], 1.0}, // green
            {p[4], p[5], 1.0}  // blue
        };

        cmsWhitePointFromTemp(&xyD, (double)temp);
        cmsToneCurve* GammaTRC[3];


        // Calculate output profile's rTRC gTRC bTRC
        GammaTRC[0] = GammaTRC[1] = GammaTRC[2] = cmsBuildParametricToneCurve(nullptr, 5, Parameters);

        if (icm.wprofile == "v4") {
            outputProfile = cmsCreateRGBProfile(&xyD, &Primaries, GammaTRC);
        }

        cmsWriteTag(outputProfile, cmsSigRedTRCTag, GammaTRC[0]);
        cmsWriteTag(outputProfile, cmsSigGreenTRCTag, GammaTRC[1]);
        cmsWriteTag(outputProfile, cmsSigBlueTRCTag, GammaTRC[2]);
        cmsWriteTag(outputProfile, cmsSigProfileDescriptionTag,  mlu);//desc changed

        /*  //to read XYZ values
            cmsCIEXYZ *redT = static_cast<cmsCIEXYZ*>(cmsReadTag(outputProfile, cmsSigRedMatrixColumnTag));
            cmsCIEXYZ *greenT  = static_cast<cmsCIEXYZ*>(cmsReadTag(outputProfile, cmsSigGreenMatrixColumnTag));
            cmsCIEXYZ *blueT  = static_cast<cmsCIEXYZ*>(cmsReadTag(outputProfile, cmsSigBlueMatrixColumnTag));
            printf("rx=%f gx=%f bx=%f ry=%f gy=%f by=%f rz=%f gz=%f bz=%f\n", redT->X, greenT->X, blueT->X, redT->Y, greenT->Y, blueT->Y, redT->Z, greenT->Z, blueT->Z);
        */

        cmsMLUfree(mlu);



        if (icm.wprofile == "v2" || icm.wprofile == "v4") {
            cmsSaveProfileToFile(outputProfile,  outPro.c_str());

        }

        if (GammaTRC) {
            cmsFreeToneCurve(GammaTRC[0]);
        }
    }

    return outputProfile;

}
