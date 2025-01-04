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

#include <memory>
#include <vector>

#include "rtengine/diagonalcurvetypes.h"
#include "rtengine/flatcurvetypes.h"

struct ParamsEdited;

namespace rtengine
{

namespace procparams
{
class ProcParams;
class PartialProfile;
class IPTCPairs;

}

}

class Clipboard
{
public:
    Clipboard ();
    ~Clipboard ();

    bool hasIPTC() const;
    const rtengine::procparams::IPTCPairs& getIPTC() const;
    void setIPTC(const rtengine::procparams::IPTCPairs& iptcc);

    const rtengine::procparams::PartialProfile& getPartialProfile() const;
    void setPartialProfile(const rtengine::procparams::PartialProfile& pprofile);

    const rtengine::procparams::ProcParams& getProcParams() const;
    void setProcParams(const rtengine::procparams::ProcParams& pparams);

    const ParamsEdited& getParamsEdited() const;

    bool hasProcParams() const;
    bool hasPEdited() const;

    DiagonalCurveType hasDiagonalCurveData() const;
    const std::vector<double>& getDiagonalCurveData() const;
    void setDiagonalCurveData(const std::vector<double>& p, DiagonalCurveType type);

    void setFlatCurveData(const std::vector<double>& p, FlatCurveType type);
    const std::vector<double>& getFlatCurveData() const;
    FlatCurveType hasFlatCurveData() const;

private:
    bool _hasIPTC;
    const std::unique_ptr<rtengine::procparams::IPTCPairs> iptc;
    const std::unique_ptr<rtengine::procparams::PartialProfile> partProfile;
    DiagonalCurveType hasDiagonalCurveDataType;
    FlatCurveType hasFlatCurveDataType;
    std::vector<double> diagonalCurve;
    std::vector<double> flatCurve;
};

extern Clipboard clipboard;
