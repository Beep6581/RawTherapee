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

#pragma once

#include <memory>
#include <set>
#include <vector>

#include <glibmm/ustring.h>

#include <lensfun.h>

#include "lcp.h"
#include "noncopyable.h"

namespace rtengine
{

class FramesMetaData;

namespace procparams
{

struct CoarseTransformParams;
struct LensProfParams;

}

class LFModifier final :
    public LensCorrection,
    public NonCopyable
{
public:
    ~LFModifier() override;

    explicit operator bool() const;

    void correctDistortion(double &x, double &y, int cx, int cy) const override;
    bool isCACorrectionAvailable() const override;
    void correctCA(double &x, double &y, int cx, int cy, int channel) const override;
    void processVignette(int width, int height, float** rawData) const override;
    void processVignette3Channels(int width, int height, float** rawData) const override;

    Glib::ustring getDisplayString() const;

private:
    LFModifier(lfModifier *m, bool swap_xy, int flags);

    friend class LFDatabase;
    lfModifier *data_;
    bool swap_xy_;
    int flags_;
};

class LFCamera final
{
public:
    LFCamera();

    explicit operator bool() const;

    Glib::ustring getMake() const;
    Glib::ustring getModel() const;
    float getCropFactor() const;
    bool isFixedLens() const;

    Glib::ustring getDisplayString() const;

private:
    friend class LFDatabase;
    const lfCamera *data_;
};

class LFLens final
{
public:
    LFLens();

    explicit operator bool() const;

    Glib::ustring getMake() const;
    Glib::ustring getLens() const;
    Glib::ustring getDisplayString() const { return getLens(); }
    float getCropFactor() const;
    bool hasVignettingCorrection() const;
    bool hasDistortionCorrection() const;
    bool hasCACorrection() const;

private:
    friend class LFDatabase;
    const lfLens *data_;
};

class LFDatabase final :
    public NonCopyable
{
public:
    static bool init(const Glib::ustring &dbdir);
    static const LFDatabase *getInstance();

    ~LFDatabase();

    std::vector<LFCamera> getCameras() const;
    std::vector<LFLens> getLenses() const;
    LFCamera findCamera(const Glib::ustring &make, const Glib::ustring &model) const;
    LFLens findLens(const LFCamera &camera, const Glib::ustring &name) const;

    std::unique_ptr<LFModifier> findModifier(
        const procparams::LensProfParams &lensProf,
        const FramesMetaData *idata,
        int width, int height,
        const procparams::CoarseTransformParams &coarse,
        int rawRotationDeg
    ) const;

private:
    std::unique_ptr<LFModifier> getModifier(const LFCamera &camera, const LFLens &lens,
                                            float focalLen, float aperture, float focusDist,
                                            int width, int height, bool swap_xy) const;
    LFDatabase();
    bool LoadDirectory(const char *dirname);

    mutable MyMutex lfDBMutex;
    static LFDatabase instance_;
    lfDatabase *data_;
    mutable std::set<std::string> notFound;
};

} // namespace rtengine
