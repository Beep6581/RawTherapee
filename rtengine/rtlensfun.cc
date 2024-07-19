/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2017 Alberto Griggio <alberto.griggio@gmail.com>
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

#include <iostream>
#include <regex>

#include "imagedata.h"
#include "procparams.h"
#include "rtlensfun.h"
#include "settings.h"


namespace
{

bool isCStringIn(const char *str, const char *const *list)
{
    for (auto element_ptr = list; *element_ptr; element_ptr++) {
        if (!strcmp(str, *element_ptr)) {
            return true;
        }
    }
    return false;
}

bool isNextLensCropFactorBetter(const lfLens *current_lens, const lfCamera *camera, float next_lens_crop_factor)
{
    if (!current_lens) {
        // No current lens, so next lens's crop factor is
        // automatically better.
        return true;
    }

    const float current_lens_crop_factor = current_lens->CropFactor;

    if (!camera) {
        // Favor the smaller crop factor for maximum coverage.
        return current_lens_crop_factor > next_lens_crop_factor;
    }

    const float camera_crop_factor = camera->CropFactor;

    if (current_lens_crop_factor > camera_crop_factor) {
        // Current lens's data does not cover the entire camera
        // sensor. Any lens's data with a smaller crop factor is
        // better.
        return current_lens->CropFactor > next_lens_crop_factor;
    }

    // Current lens's data covers the entire camera sensor. A lens
    // with data from a larger crop factor will be more precise, but
    // also must not be larger than the camera sensor's crop factor
    // to maintain full coverage.
    return current_lens->CropFactor < next_lens_crop_factor &&
           next_lens_crop_factor <= camera_crop_factor;
}

bool isNextLensBetter(const lfCamera *camera, const lfLens *current_lens, const lfLens &next_lens, const Glib::ustring &lens_name, const Glib::ustring &next_lens_name)
{
    return isNextLensCropFactorBetter(current_lens, camera, next_lens.CropFactor) &&
           lens_name == next_lens_name &&
           (!camera || isCStringIn(camera->Mount, next_lens.Mounts));
}

} // namespace

namespace rtengine
{

//-----------------------------------------------------------------------------
// LFModifier
//-----------------------------------------------------------------------------

LFModifier::~LFModifier()
{
    if (data_) {
        data_->Destroy();
    }
}


LFModifier::operator bool() const
{
    return data_;
}


bool LFModifier::hasDistortionCorrection() const
{
    return (flags_ & LF_MODIFY_DISTORTION);
}

bool LFModifier::hasCACorrection() const
{
    return (flags_ & LF_MODIFY_TCA);
}

bool LFModifier::hasVignettingCorrection() const
{
    return (flags_ & LF_MODIFY_VIGNETTING);
}

void LFModifier::correctDistortion(double &x, double &y, int cx, int cy) const
{
    if (!data_) {
        return;
    }

    float pos[2];
    float xx = x + cx;
    float yy = y + cy;
    if (swap_xy_) {
        std::swap(xx, yy);
    }
    if (data_->ApplyGeometryDistortion(xx, yy, 1, 1, pos)) {  // This is thread-safe
        x = pos[0];
        y = pos[1];
        if (swap_xy_) {
            std::swap(x, y);
        }
        x -= cx;
        y -= cy;
    }
}

void LFModifier::correctCA(double &x, double &y, int cx, int cy, int channel) const
{
    assert(channel >= 0 && channel <= 2);

    // agriggio: RT currently applies the CA correction per channel, whereas
    // lensfun applies it to all the three channels simultaneously. This means
    // we do the work 3 times, because each time we discard 2 of the 3
    // channels. We could consider caching the info to speed this up
    x += cx;
    y += cy;

    float pos[6];
    if (swap_xy_) {
        std::swap(x, y);
    }
    data_->ApplySubpixelDistortion(x, y, 1, 1, pos);  // This is thread-safe
    x = pos[2*channel];
    y = pos[2*channel+1];
    if (swap_xy_) {
        std::swap(x, y);
    }
    x -= cx;
    y -= cy;
}

void LFModifier::correctDistortionAndCA(double &x, double &y, int cx, int cy, int channel) const
{
    assert(channel >= 0 && channel <= 2);

    // RT currently applies the CA correction per channel, whereas
    // lensfun applies it to all the three channels simultaneously. This means
    // we do the work 3 times, because each time we discard 2 of the 3
    // channels. We could consider caching the info to speed this up
    x += cx;
    y += cy;

    float pos[6];
    if (swap_xy_) {
        std::swap(x, y);
    }
    data_->ApplySubpixelGeometryDistortion(x, y, 1, 1, pos);  // This is thread-safe
    x = pos[2*channel];
    y = pos[2*channel+1];
    if (swap_xy_) {
        std::swap(x, y);
    }
    x -= cx;
    y -= cy;
}

#ifdef _OPENMP
void LFModifier::processVignette(int width, int height, float** rawData) const
{
    #pragma omp parallel for schedule(dynamic,16)

    for (int y = 0; y < height; ++y) {
        data_->ApplyColorModification(rawData[y], 0, y, width, 1, LF_CR_1(INTENSITY), 0);
    }
}
#else
void LFModifier::processVignette(int width, int height, float** rawData) const
{
    data_->ApplyColorModification(rawData[0], 0, 0, width, height, LF_CR_1(INTENSITY), width * sizeof(float));
}

#endif

#ifdef _OPENMP
void LFModifier::processVignette3Channels(int width, int height, float** rawData) const
{
    #pragma omp parallel for schedule(dynamic,16)

    for (int y = 0; y < height; ++y) {
        data_->ApplyColorModification(rawData[y], 0, y, width, 1, LF_CR_3(RED, GREEN, BLUE), 0);
    }
}
#else
void LFModifier::processVignette3Channels(int width, int height, float** rawData) const
{
    data_->ApplyColorModification(rawData[0], 0, 0, width, height, LF_CR_3(RED, GREEN, BLUE), width * 3 * sizeof(float));
}

#endif
Glib::ustring LFModifier::getDisplayString() const
{
    if (!data_) {
        return "NONE";
    } else {
        Glib::ustring ret;
        Glib::ustring sep;
        if (flags_ & LF_MODIFY_DISTORTION) {
            ret += "distortion";
            sep = ", ";
        }
        if (flags_ & LF_MODIFY_VIGNETTING) {
            ret += sep;
            ret += "vignetting";
            sep = ", ";
        }
        if (flags_ & LF_MODIFY_TCA) {
            ret += sep;
            ret += "CA";
            sep = ", ";
        }
        if (flags_ & LF_MODIFY_SCALE) {
            ret += sep;
            ret += "autoscaling";
        }
        return ret;
    }
}


LFModifier::LFModifier(lfModifier *m, bool swap_xy, int flags):
    data_(m),
    swap_xy_(swap_xy),
    flags_(flags)
{
}


//-----------------------------------------------------------------------------
// LFCamera
//-----------------------------------------------------------------------------

LFCamera::LFCamera():
    data_(nullptr)
{
}


LFCamera::operator bool() const
{
    return data_;
}


Glib::ustring LFCamera::getMake() const
{
    if (data_) {
        return data_->Maker;
    } else {
        return "";
    }
}


Glib::ustring LFCamera::getModel() const
{
    if (data_) {
        return data_->Model;
    } else {
        return "";
    }
}


float LFCamera::getCropFactor() const
{
    if (data_) {
        return data_->CropFactor;
    } else {
        return 0;
    }
}


bool LFCamera::isFixedLens() const
{
    // per lensfun's main developer Torsten Bronger:
    // "Compact camera mounts can be identified by the fact that the mount
    // starts with a lowercase letter"
    return data_ && data_->Mount && std::islower(data_->Mount[0]);
}


Glib::ustring LFCamera::getDisplayString() const
{
    if (data_) {
        return getMake() + ' ' + getModel();
    } else {
        return "---";
    }
}


//-----------------------------------------------------------------------------
// LFLens
//-----------------------------------------------------------------------------

LFLens::LFLens():
    data_(nullptr)
{
}


LFLens::operator bool() const
{
    return data_;
}


Glib::ustring LFLens::getMake() const
{
    if (data_) {
        return data_->Maker;
    } else {
        return "";
    }
}


Glib::ustring LFLens::getLens() const
{
    if (data_) {
        return Glib::ustring(data_->Maker) + ' ' + data_->Model;
    } else {
        return "---";
    }
}


float LFLens::getCropFactor() const
{
    if (data_) {
        return data_->CropFactor;
    } else {
        return 0;
    }
}

bool LFLens::hasVignettingCorrection() const
{
    if (data_) {
        return data_->CalibVignetting;
    } else {
        return false;
    }
}

bool LFLens::hasDistortionCorrection() const
{
    if (data_) {
        return data_->CalibDistortion;
    } else {
        return false;
    }
}

bool LFLens::hasCACorrection() const
{
    if (data_) {
        return data_->CalibTCA;
    } else {
        return false;
    }
}


//-----------------------------------------------------------------------------
// LFDatabase
//-----------------------------------------------------------------------------

LFDatabase LFDatabase::instance_;


bool LFDatabase::init(const Glib::ustring &dbdir)
{
    instance_.data_ = lfDatabase::Create();

    if (settings->verbose) {
        std::cout << "Loading lensfun database from ";
        if (dbdir.empty()) {
            std::cout << "the default directories";
        } else {
            std::cout << "'" << dbdir << "'";
        }
        std::cout << "..." << std::flush;
    }

    bool ok = false;
    if (dbdir.empty()) {
        ok = (instance_.data_->Load() ==  LF_NO_ERROR);
    } else {
        ok = instance_.LoadDirectory(dbdir.c_str());
    }

    if (settings->verbose) {
        std::cout << (ok ? "OK" : "FAIL") << std::endl;
    }

    return ok;
}


bool LFDatabase::LoadDirectory(const char *dirname)
{
#if RT_LENSFUN_HAS_LOAD_DIRECTORY
    return instance_.data_->LoadDirectory(dirname);
#else
    // backported from lensfun 0.3.x
    bool database_found = false;

    GDir *dir = g_dir_open (dirname, 0, NULL);
    if (dir)
    {
        GPatternSpec *ps = g_pattern_spec_new ("*.xml");
        if (ps)
        {
            const gchar *fn;
            while ((fn = g_dir_read_name (dir)))
            {
                size_t sl = strlen (fn);
                if (g_pattern_match (ps, sl, fn, NULL))
                {
                    gchar *ffn = g_build_filename (dirname, fn, NULL);
                    /* Ignore errors */
                    if (data_->Load (ffn) == LF_NO_ERROR)
                        database_found = true;
                    g_free (ffn);
                }
            }
            g_pattern_spec_free (ps);
        }
        g_dir_close (dir);
    }

    return database_found;
#endif
}


LFDatabase::LFDatabase():
    data_(nullptr)
{
}


LFDatabase::~LFDatabase()
{
    if (data_) {
        MyMutex::MyLock lock(lfDBMutex);
        data_->Destroy();
    }
}


const LFDatabase *LFDatabase::getInstance()
{
    return &instance_;
}


std::vector<LFCamera> LFDatabase::getCameras() const
{
    std::vector<LFCamera> ret;
    if (data_) {
        MyMutex::MyLock lock(lfDBMutex);
        auto cams = data_->GetCameras();
        while (*cams) {
            ret.emplace_back();
            ret.back().data_ = *cams;
            ++cams;
        }
    }
    return ret;
}


std::vector<LFLens> LFDatabase::getLenses() const
{
    std::vector<LFLens> ret;
    if (data_) {
        MyMutex::MyLock lock(lfDBMutex);
        auto lenses = data_->GetLenses();
        while (*lenses) {
            ret.emplace_back();
            ret.back().data_ = *lenses;
            ++lenses;
        }
    }
    return ret;
}


LFCamera LFDatabase::findCamera(const Glib::ustring &make, const Glib::ustring &model, bool autoMatch) const
{
    LFCamera ret;
    if (data_ && !make.empty()) {
        MyMutex::MyLock lock(lfDBMutex);
        if (!autoMatch) {
            // Try to find exact match by name.
            for (auto camera_list = data_->GetCameras(); camera_list[0]; camera_list++) {
                const auto camera = camera_list[0];
                if (make == camera->Maker && model == camera->Model) {
                    ret.data_ = camera;
                    return ret;
                }
            }
        }
        auto found = data_->FindCamerasExt(make.c_str(), model.c_str());
        if (found) {
            ret.data_ = found[0];
            lf_free(found);
        }
    }
    return ret;
}


LFLens LFDatabase::findLens(const LFCamera &camera, const Glib::ustring &name, bool autoMatch) const
{
    LFLens ret;
    if (data_ && !name.empty()) {
        MyMutex::MyLock lock(lfDBMutex);
        if (!autoMatch) {
            // Only the lens name provided. Try to find exact match by name.
            LFLens candidate;
            LFLens bestCandidate;

            for (auto lens_list = data_->GetLenses(); lens_list[0]; lens_list++) {
                candidate.data_ = lens_list[0];
                if (isNextLensBetter(camera.data_, bestCandidate.data_, *(candidate.data_), name, candidate.getLens())) {
                    bestCandidate.data_ = candidate.data_;
                }
            }
            if (bestCandidate.data_) {
                return bestCandidate;
            }
        }
        const auto find_lens_from_name = [](const lfDatabase *database, const lfCamera *cam, const Glib::ustring &lens_name) {
            auto found = database->FindLenses(cam, nullptr, lens_name.c_str());
            for (size_t pos = 0; !found && pos < lens_name.size(); ) {
                // try to split the maker from the model of the lens -- we have to
                // guess a bit here, since there are makers with a multi-word name
                // (e.g. "Leica Camera AG")
                if (lens_name.find("f/", pos) == 0) {
                    break; // no need to search further
                }
                Glib::ustring make, model;
                auto i = lens_name.find(' ', pos);
                if (i != Glib::ustring::npos) {
                    make = lens_name.substr(0, i);
                    model = lens_name.substr(i+1);
                    found = database->FindLenses(cam, make.c_str(), model.c_str());
                    pos = i+1;
                } else {
                    break;
                }
            }
            return found;
        };
        auto found = find_lens_from_name(data_, camera.data_, name);
        if (!found) {
            // Some names have white-space around the dash(s) while Lensfun does
            // not have any.
            const std::regex pattern("\\s*-\\s*");
            const auto formatted_name = std::regex_replace(name.raw(), pattern, "-");
            if (name != formatted_name) {
                found = find_lens_from_name(data_, camera.data_, formatted_name);
            }
        }
        if (!found && camera && camera.isFixedLens()) {
            found = data_->FindLenses(camera.data_, nullptr, "");
        }
        if (found) {
            ret.data_ = found[0];
            lf_free(found);
        }
    }
    return ret;
}


std::unique_ptr<LFModifier> LFDatabase::getModifier(const LFCamera &camera, const LFLens &lens,
                                    float focalLen, float aperture, float focusDist,
                                    int width, int height, bool swap_xy) const
{
    std::unique_ptr<LFModifier> ret;
    if (data_) {
        MyMutex::MyLock lock(lfDBMutex);
        if (camera && lens) {
            lfModifier *mod = lfModifier::Create(lens.data_, camera.getCropFactor(), width, height);
            int flags = LF_MODIFY_DISTORTION | LF_MODIFY_SCALE | LF_MODIFY_TCA;
            if (aperture > 0) {
                flags |= LF_MODIFY_VIGNETTING;
            }
            flags = mod->Initialize(lens.data_, LF_PF_F32, focalLen, aperture, focusDist > 0 ? focusDist : 1000, 0.0, LF_RECTILINEAR, flags, false);
            ret.reset(new LFModifier(mod, swap_xy, flags));
        }
    }
    return ret;
}

std::unique_ptr<LFModifier> LFDatabase::findModifier(
    const procparams::LensProfParams &lensProf,
    const FramesMetaData *idata,
    int width,
    int height,
    const procparams::CoarseTransformParams &coarse,
    int rawRotationDeg
) const
{
    const float focallen = idata->getFocalLen();

    Glib::ustring make, model, lens;
    if (lensProf.lfAutoMatch()) {
        if (focallen <= 0.f) {
            return nullptr;
        }
        make = idata->getMake();
        model = idata->getModel();
        lens = idata->getLens();
    } else {
        make = lensProf.lfCameraMake;
        model = lensProf.lfCameraModel;
        lens = lensProf.lfLens;
    }

    if (make.empty() || model.empty() || lens.empty()) {
        return nullptr;
    }

    const std::string key = (make + model + lens).collate_key();
    if (notFound.find(key) != notFound.end()) {
        // This combination was not found => do not search again
        return nullptr;
    }

    const LFCamera c = findCamera(make, model, lensProf.lfAutoMatch());
    const LFLens l = findLens(c, lens, lensProf.lfAutoMatch());

    bool swap_xy = false;
    if (rawRotationDeg >= 0) {
        int rot = (coarse.rotate + rawRotationDeg) % 360;
        swap_xy = (rot == 90 || rot == 270);
        if (swap_xy) {
            std::swap(width, height);
        }
    }

    std::unique_ptr<LFModifier> ret = getModifier(
        c,
        l,
        idata->getFocalLen(),
        idata->getFNumber(),
        idata->getFocusDist(),
        width,
        height,
        swap_xy
    );

    if (settings->verbose) {
        std::cout << "LENSFUN:\n"
                  << "  camera: " << c.getDisplayString() << "\n"
                  << "  lens: " << l.getDisplayString() << "\n"
                  << "  correction: "
                  << (ret ? ret->getDisplayString() : "NONE")
                  << std::endl;
    }

    if (!ret) {
        notFound.insert(key);
    }

    return ret;
}


} // namespace rtengine
