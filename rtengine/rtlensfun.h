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

#pragma once

#include <lensfun.h>
#include <glibmm.h>
#include <memory>

namespace rtengine {

class LFModifier {
public:
    bool ok() const;
    
    void correctDistortion(double &x, double &y);
    void processVignetteLine(int width, int y, float *line);
    void processVignetteLine3Channels(int width, int y, float *line);
    
private:
    friend class LFDatabase;
    std::shared_ptr<lfModifier> data_;
};

class LFCamera {
public:
    bool ok() const;
    
    Glib::ustring getMake() const;
    Glib::ustring getModel() const;
    float getCropFactor() const;

    Glib::ustring getDisplayString() const;

private:
    friend class LFDatabase;
    std::shared_ptr<lfCamera> data_;
};

class LFLens {
public:
    bool ok() const;
    
    Glib::ustring getDisplayString() const;
private:
    friend class LFDatabase;
    std::shared_ptr<lfLens> data_;
};

class LFDatabase {
public:
    static bool init();
    static LFDatabase *getInstance();

    std::vector<LFCamera> getCameras();
    std::vector<LFLens> getLenses(const LFCamera &camera);
    LFCamera findCamera(const Glib::ustring &make, const Glib::ustring &model);
    LFLens findLens(const LFCamera &camera, const Glib::ustring &name);
    LFModifier getModifier(const LFCamera &camera, const LFLens &lens,
                           int width, int height,
                           float focalLen, float aperture);

private:
    static LFDatabase instance_;
    std::shared_ptr<lfDatabase> data_;
};

} // namespace rtengine
