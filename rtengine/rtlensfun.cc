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

namespace rtengine {

//-----------------------------------------------------------------------------
// LFModifier
//-----------------------------------------------------------------------------

LFModifier::LFModifier(lfModifier *m):
    data_(m)
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
    if (data_->ApplyGeometryDistortion(x+cx, y+cy, 1, 1, pos)) {
        x = pos[0] - cx;
        y = pos[1] - cy;
    }
}


void LFModifier::processVignetteLine(int width, int y, float *line) const
{
    // TODO
}


void LFModifier::processVignetteLine3Channels(int width, int y, float *line) const
{
    // TODO
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


Glib::ustring LFLens::getDisplayString() const
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


std::vector<LFLens> LFDatabase::getLenses(const LFCamera &camera) const
{
    std::vector<LFLens> ret;
    if (data_) {
        auto lenses = data_->FindLenses(camera.data_, NULL, "", LF_SEARCH_LOOSE /*| LF_SEARCH_SORT_AND_UNIQUIFY*/);
        while (*lenses) {
            ret.emplace_back(LFLens());
            ret.back().data_ = *lenses;
            ++lenses;
        }
        lf_free(lenses);
    }
    return ret;
}


LFCamera LFDatabase::findCamera(const Glib::ustring &make, const Glib::ustring &model) const
{
    LFCamera ret;
    if (data_) {
        auto found = data_->FindCamerasExt(make.c_str(), model.c_str(), LF_SEARCH_LOOSE);
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
        auto found = data_->FindLenses(camera.data_, NULL, name.c_str(), LF_SEARCH_LOOSE);
        if (!found) {
            // try to split the maker from the model of the lens
            Glib::ustring make, model;
            auto i = name.find_first_of(' ');
            if (i != Glib::ustring::npos) {
                make = name.substr(0, i);
                model = name.substr(i+1);
                found = data_->FindLenses(camera.data_, make.c_str(), model.c_str(), LF_SEARCH_LOOSE);
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
                                    int width, int height, float focalLen,
                                    float aperture, float focusDist) const
{
    LFModifier *ret = nullptr;
    if (data_) {
        if (camera.ok() && lens.ok()) {
            lfModifier *mod = lfModifier::Create(lens.data_, camera.getCropFactor(), width, height);
            mod->Initialize(lens.data_, LF_PF_F32, focalLen, aperture, focusDist > 0 ? focusDist : 1000, 0.0, LF_RECTILINEAR, LF_MODIFY_VIGNETTING | LF_MODIFY_DISTORTION, false);
            ret = new LFModifier(mod);
        }
    }
    return ret;
}
    

} // namespace rtengine
