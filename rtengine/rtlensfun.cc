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

LFModifier::LFModifier(lfModifier *m, bool swap_xy, int flags):
    data_(m),
    swap_xy_(swap_xy),
    flags_(flags)
{
}


LFModifier::~LFModifier()
{
    if (data_) {
        data_->Destroy();
    }
}

bool LFModifier::ok() const
{
    return data_;
}


void LFModifier::correctDistortion(double &x, double &y, int cx, int cy, double scale) const
{
    if (!ok()) {
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
        if (flags_ & LF_MODIFY_SCALE) {
            ret += sep;
            ret += "autoscaling";
        }
        return ret;
    }
}


//-----------------------------------------------------------------------------
// LFCamera
//-----------------------------------------------------------------------------

LFCamera::LFCamera():
    data_(nullptr)
{
}


bool LFCamera::ok() const
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
        return Glib::ustring::compose("%1 %2", getMake(), getModel());
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


bool LFLens::ok() const
{
    return data_;
}


Glib::ustring LFLens::getLens() const
{
    if (data_) {
        return Glib::ustring::compose("%1 %2", data_->Maker, data_->Model);
    } else {
        return "---";
    }
}


//-----------------------------------------------------------------------------
// LFDatabase
//-----------------------------------------------------------------------------

LFDatabase LFDatabase::instance_;


bool LFDatabase::init()
{
    instance_.data_ = lfDatabase::Create();
    return instance_.data_->Load() != LF_NO_ERROR;
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
            ret.emplace_back(LFCamera());
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
            ret.emplace_back(LFLens());
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
        bool stdlens = camera.ok() && (name.empty() || name.find("Unknown ") == 0);
        if (stdlens) {
            lname = camera.getModel(); // "Standard"
        }
        auto found = data_->FindLenses(camera.data_, nullptr, lname.c_str());
        if (!found) {
            // try to split the maker from the model of the lens
            Glib::ustring make, model;
            auto i = name.find_first_of(' ');
            if (i != Glib::ustring::npos) {
                make = name.substr(0, i);
                model = name.substr(i+1);
                found = data_->FindLenses(camera.data_, make.c_str(), model.c_str());
            }
        }
        if (found) {
            ret.data_ = found[0];
            lf_free(found);
        }
    }
    return ret;
}


LFModifier *LFDatabase::getModifier(const LFCamera &camera, const LFLens &lens,
                                    float focalLen, float aperture, float focusDist,
                                    int width, int height, bool swap_xy) const
{
    LFModifier *ret = nullptr;
    if (data_) {
        if (camera.ok() && lens.ok()) {
            lfModifier *mod = lfModifier::Create(lens.data_, camera.getCropFactor(), width, height);
            int flags = LF_MODIFY_DISTORTION | LF_MODIFY_SCALE;
            if (aperture > 0) {
                flags |= LF_MODIFY_VIGNETTING;
            }
            flags = mod->Initialize(lens.data_, LF_PF_F32, focalLen, aperture, focusDist > 0 ? focusDist : 1000, 0.0, LF_RECTILINEAR, flags, false);
            ret = new LFModifier(mod, swap_xy, flags);
        }
    }
    return ret;
}


LFModifier *LFDatabase::findModifier(const LensProfParams &lensProf, const ImageMetaData *idata, int width, int height, const CoarseTransformParams &coarse, int rawRotationDeg)
{
    const LFDatabase *db = getInstance();
    Glib::ustring make, model, lens;
    float focallen = idata->getFocalLen();
    if (lensProf.lfAutoMatch) {
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
    LFLens l = db->findLens(lensProf.lfAutoMatch ? c : LFCamera(), lens);
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
    
    LFModifier *ret = db->getModifier(c, l, idata->getFocalLen(), idata->getFNumber(), idata->getFocusDist(), width, height, swap_xy);

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
