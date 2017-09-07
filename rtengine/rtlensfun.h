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
#include "lcp.h"
#include "procparams.h"

namespace rtengine {

class LFModifier: public LensCorrection {
public:
    ~LFModifier();
    bool ok() const;
    
    void correctDistortion(double &x, double &y, int cx, int cy, double scale) const;
    bool supportsAutoFill() const { return false; }
    bool supportsCA() const { return false; }
    void correctCA(double &x, double &y, int channel) const {}
    void processVignetteLine(int width, int y, float *line) const;
    void processVignetteLine3Channels(int width, int y, float *line) const;
    
private:
    explicit LFModifier(lfModifier *m, bool rotateXY);
    LFModifier(const LFModifier &);
    LFModifier &operator=(const LFModifier &);
    
    friend class LFDatabase;
    lfModifier *data_;
    bool swap_xy_;
};

class LFCamera {
public:
    LFCamera();
    bool ok() const;
    
    Glib::ustring getMake() const;
    Glib::ustring getModel() const;
    float getCropFactor() const;

    Glib::ustring getDisplayString() const;

private:
    friend class LFDatabase;
    const lfCamera *data_;
};

class LFLens {
public:
    LFLens();
    bool ok() const;
    
    Glib::ustring getDisplayString() const;
private:
    friend class LFDatabase;
    const lfLens *data_;
};

class LFDatabase {
public:
    static bool init();
    static const LFDatabase *getInstance();

    ~LFDatabase();
    
    std::vector<LFCamera> getCameras() const;
    std::vector<LFLens> getLenses() const;
    LFCamera findCamera(const Glib::ustring &make, const Glib::ustring &model) const;
    LFLens findLens(const LFCamera &camera, const Glib::ustring &name) const;
    LFModifier *getModifier(const LFCamera &camera, const LFLens &lens,
                            float focalLen, float aperture, float focusDist,
                            int width, int height, bool swap_xy) const;

    static LFModifier *findModifier(const LensProfParams &lensProf, const ImageMetaData *idata, int width, int height, const CoarseTransformParams &coarse, int rawRotationDeg);

private:
    LFDatabase();
    LFDatabase(const LFDatabase &);
    LFDatabase &operator=(const LFDatabase &);
    static LFDatabase instance_;
    lfDatabase *data_;
};

} // namespace rtengine
