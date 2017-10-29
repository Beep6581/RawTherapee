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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "rtlensfun.h"
#include "settings.h"
#include <iostream>

namespace rtengine {

extern const Settings *settings;

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


void LFModifier::correctDistortion(double &x, double &y, int cx, int cy, double scale) const
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
    if (data_->ApplyGeometryDistortion(xx, yy, 1, 1, pos)) {
        x = pos[0];
        y = pos[1];
        if (swap_xy_) {
            std::swap(x, y);
        }
        x -= cx;
        y -= cy;
    }
    x *= scale;
    y *= scale;
}


bool LFModifier::isCACorrectionAvailable() const
{
    return (flags_ & LF_MODIFY_TCA);
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
    data_->ApplySubpixelDistortion(x, y, 1, 1, pos);
    x = pos[2*channel];
    y = pos[2*channel+1];
    if (swap_xy_) {
        std::swap(x, y);
    }
    x -= cx;
    y -= cy;
}


void LFModifier::processVignetteLine(int width, int y, float *line) const
{
    data_->ApplyColorModification(line, 0, y, width, 1, LF_CR_1(INTENSITY), 0);
}


void LFModifier::processVignetteLine3Channels(int width, int y, float *line) const
{
    data_->ApplyColorModification(line, 0, y, width, 1, LF_CR_3(RED, GREEN, BLUE), 0);
}


Glib::ustring LFModifier::getDisplayString() const
{
    if (!data_) {
        return "NONE";
    } else {
        Glib::ustring ret;
        Glib::ustring sep = "";
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
        auto lenses = data_->GetLenses();
        while (*lenses) {
            ret.emplace_back();
            ret.back().data_ = *lenses;
            ++lenses;
        }
    }
    return ret;
}


LFCamera LFDatabase::findCamera(const Glib::ustring &make, const Glib::ustring &model) const
{
    LFCamera ret;
    if (data_) {
        auto found = data_->FindCamerasExt(make.c_str(), model.c_str());
        if (found) {
            ret.data_ = found[0];
            lf_free(found);
        }
    }
    return ret;
}


LFLens LFDatabase::findLens(const LFCamera &camera, const Glib::ustring &name) const
{
    LFLens ret;
    if (data_) {
        Glib::ustring lname = name;
        bool stdlens = camera && (name.empty() || name.find("Unknown") == 0);
        if (stdlens) {
            lname = camera.getModel(); // "Standard"
        }
        auto found = data_->FindLenses(camera.data_, nullptr, lname.c_str());
        for (size_t pos = 0; !found && pos < name.size(); ) {
            // try to split the maker from the model of the lens -- we have to
            // guess a bit here, since there are makers with a multi-word name
            // (e.g. "Leica Camera AG")
            if (name.find("f/", pos) == 0) {
                break; // no need to search further
            }
            Glib::ustring make, model;
            auto i = name.find(' ', pos);
            if (i != Glib::ustring::npos) {
                make = name.substr(0, i);
                model = name.substr(i+1);
                found = data_->FindLenses(camera.data_, make.c_str(), model.c_str());
                pos = i+1;
            } else {
                break;
            }
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


std::unique_ptr<LFModifier> LFDatabase::findModifier(const LensProfParams &lensProf, const FramesMetaData *idata, int width, int height, const CoarseTransformParams &coarse, int rawRotationDeg)
{
    const LFDatabase *db = getInstance();
    Glib::ustring make, model, lens;
    float focallen = idata->getFocalLen();
    if (lensProf.lfAutoMatch()) {
        if (focallen <= 0) {
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
    LFCamera c = db->findCamera(make, model);
    LFLens l = db->findLens(lensProf.lfAutoMatch() ? c : LFCamera(), lens);
    if (focallen <= 0 && l.data_ && l.data_->MinFocal == l.data_->MaxFocal) {
        focallen = l.data_->MinFocal;
    }
    if (focallen <= 0) {
        return nullptr;
    }
    bool swap_xy = false;
    if (rawRotationDeg >= 0) {
        int rot = (coarse.rotate + rawRotationDeg) % 360;
        swap_xy = (rot == 90 || rot == 270);
        if (swap_xy) {
            std::swap(width, height);
        }
    }

    std::unique_ptr<LFModifier> ret = db->getModifier(c, l, idata->getFocalLen(), idata->getFNumber(), idata->getFocusDist(), width, height, swap_xy);

    if (settings->verbose) {
        std::cout << "LENSFUN:\n"
                  << "  camera: " << c.getDisplayString() << "\n"
                  << "  lens: " << l.getDisplayString() << "\n"
                  << "  correction: "
                  << (ret ? ret->getDisplayString() : "NONE") << std::endl;
    }

    return ret;
}


} // namespace rtengine
