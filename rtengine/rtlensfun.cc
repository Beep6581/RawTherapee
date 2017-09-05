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

bool LFModifier::ok() const
{
    return data_.get();
}


void LFModifier::correctDistortion(double &x, double &y)
{
    if (!ok()) {
        return;
    }

    float pos[2];
    data_->ApplyGeometryDistortion(x, y, 1, 1, pos);
    x = pos[0];
    y = pos[1];
}


void LFModifier::processVignetteLine(int width, int y, float *line)
{
    // TODO
}


void LFModifier::processVignetteLine3Channels(int width, int y, float *line)
{
    // TODO
}


//-----------------------------------------------------------------------------
// LFCamera
//-----------------------------------------------------------------------------

bool LFCamera::ok() const
{
    return data_.get();
}


Glib::ustring LFCamera::getMake() const
{
    if (ok()) {
        return data_->Maker;
    } else {
        return "";
    }
}


Glib::ustring LFCamera::getModel() const
{
    if (ok()) {
        return data_->Model;
    } else {
        return "";
    }
}


float LFCamera::getCropFactor() const
{
    if (ok()) {
        return data_->CropFactor;
    } else {
        return 0;
    }
}


Glib::ustring LFCamera::getDisplayString() const
{
    if (ok()) {
        return Glib::ustring::compose("%1 %2", getMake(), getModel());
    } else {
        return "---";
    }
}


//-----------------------------------------------------------------------------
// LFLens
//-----------------------------------------------------------------------------

bool LFLens::ok() const
{
    return data_->get();
}


Glib::ustring LFLens::getDisplayString() const
{
    if (ok()) {
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
    instance_.data_.reset(new lfDatabase());
    return instance_.data_->Load() != LF_NO_ERROR;
}


LFDatabase *LFDatabase::getInstance()
{
    return &instance_;
}


std::vector<LFCamera> LFDatabase::getCameras()
{
    auto cams = data_->GetCameras();
    std::vector<LFCamera> ret;
    while (*cams) {
        ret.emplace_back(LFCamera());
        ret.back().data_.reset(new lfCamera(**cams));
        ++cams;
    }
    return ret;
}


std::vector<LFLens> getLenses(const LFCamera &camera)
{
    auto lenses = data_->FindLenses(*camera.data_->get(), NULL, "", LF_SEARCH_LOOSE | LF_SEARCH_SORT_AND_UNIQUIFY);
    std::vector<LFLens> ret;
    while (*lenses) {
        ret.emplace_back(LFLens());
        ret.back().data_.reset(new lfLens(**lenses));
        ++lenses;
    }
    lf_free(lenses);
    return ret;
}


LFCamera LFDatabase::findCamera(const Glib::ustring &make, const Glib::ustring &model)
{
    LFCamera ret;
    auto found = data_->FindCamerasExt(make.c_str(), model.c_str(), LF_SEARCH_LOOSE);
    if (found) {
        ret.data_.reset(new lfCamera(*found[0]));
        lf_free(found);
    }
    return ret;
}


LFLens LFDatabase::findLens(const LFCamera &camera, const Glib::ustring &name)
{
    LFLens ret;
    auto found = data_->FindLenses(camera.data_.get(), NULL, name.c_str(), LF_SEARCH_LOOSE);
    if (!found) {
        // try to split the maker from the model of the lens
        Glib::ustring make, model;
        auto i = name.find_first_of(' ');
        if (i != Glib::ustring::npos) {
            make = name.substr(0, i);
            model = name.substr(i+1);
            found = data_->FindLenses(camera.data_.get(), make.c_str(), model.c_str(), LF_SEARCH_LOOSE);
        }
    }
    if (found) {
        ret.data_.reset(new lfLens(*found[0]));
        lf_free(found);
    }
    return ret;
}


LFModifier LFDatabase::getModifier(const LFCamera &camera, const LFLens &lens,
                                   int width, int height, float focalLen,
                                   float aperture)
{
    LFModifier ret;
    if (camera.ok() && lens.ok()) {
        lfModifier *mod = lfModifier::Create(lens.data_.get(), camera.getCropFactor(), width, height);
        mod->Initialize(lens.data_.get(), LF_PF_F32, focalLen, aperture, 1000, 1, LF_RECTILINEAR, LF_MODIFY_VIGNETTING | LF_MODIFY_DISTORTION, false);
        ret.data_.reset(mod);
    }
    return ret;
}
    

} // namespace rtengine
