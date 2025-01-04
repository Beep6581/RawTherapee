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
#include "clipboard.h"

#include "paramsedited.h"
#include "rtengine/procparams.h"

Clipboard clipboard;

Clipboard::Clipboard () :
    _hasIPTC(false),
    iptc(new rtengine::procparams::IPTCPairs),
    partProfile(new rtengine::procparams::PartialProfile(false)),
    hasDiagonalCurveDataType(DCT_Empty),
    hasFlatCurveDataType(FCT_Empty)
{
}

Clipboard::~Clipboard ()
{
    partProfile->deleteInstance();
}

bool Clipboard::hasIPTC() const
{
    return _hasIPTC;
}

const rtengine::procparams::IPTCPairs& Clipboard::getIPTC() const
{
    return *iptc;
}

void Clipboard::setIPTC(const rtengine::procparams::IPTCPairs& iptcc)
{
    *iptc = iptcc;
    _hasIPTC = true;
}

const rtengine::procparams::PartialProfile& Clipboard::getPartialProfile() const
{
    return *partProfile;
}

/*
 * set both the "pparams" and "pedited" field of the PartialProfile; each one can be NULL
 */
void Clipboard::setPartialProfile(const rtengine::procparams::PartialProfile& pprofile)
{
    if (pprofile.pparams) {
        if (!partProfile->pparams) {
            partProfile->pparams = new rtengine::procparams::ProcParams();
        }

        *partProfile->pparams = *pprofile.pparams;
    } else {
        if (partProfile->pparams) {
            delete partProfile->pparams;
            partProfile->pparams = nullptr;
        }
    }

    if (pprofile.pedited) {
        if (!partProfile->pedited) {
            partProfile->pedited = new ParamsEdited();
        }

        *partProfile->pedited = *pprofile.pedited;
    } else {
        if (partProfile->pedited) {
            delete partProfile->pedited;
            partProfile->pedited = nullptr;
        }
    }
}

const rtengine::procparams::ProcParams& Clipboard::getProcParams() const
{
    return *partProfile->pparams;
}

/*
 * this method copy the procparams to "pparams" and delete "pedited"
 */
void Clipboard::setProcParams(const rtengine::procparams::ProcParams& pparams)
{
    // copy procparams
    if (!partProfile->pparams) {
        partProfile->pparams = new rtengine::procparams::ProcParams();
    }

    *partProfile->pparams = pparams;

    // delete pedited
    if (partProfile->pedited) {
        delete partProfile->pedited;
        partProfile->pedited = nullptr;
    }
}

const ParamsEdited& Clipboard::getParamsEdited() const
{
    return *partProfile->pedited;
}

bool Clipboard::hasProcParams() const
{
    return partProfile->pparams;
}

bool Clipboard::hasPEdited() const
{
    return partProfile->pedited;
}

DiagonalCurveType Clipboard::hasDiagonalCurveData() const
{
    return hasDiagonalCurveDataType;
}

const std::vector<double>& Clipboard::getDiagonalCurveData() const
{
    return diagonalCurve;
}

void Clipboard::setDiagonalCurveData(const std::vector<double>& p, DiagonalCurveType type)
{
    diagonalCurve = p;
    hasDiagonalCurveDataType = type;
    return;
}

FlatCurveType Clipboard::hasFlatCurveData() const
{
    return hasFlatCurveDataType;
}

const std::vector<double>& Clipboard:: getFlatCurveData() const
{
    return flatCurve;
}

void Clipboard::setFlatCurveData(const std::vector<double>& p, FlatCurveType type)
{
    flatCurve = p;
    hasFlatCurveDataType = type;
    return;
}
