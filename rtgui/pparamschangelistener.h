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

struct ParamsEdited;

namespace Glib
{

class ustring;

}
namespace rtengine
{

class ProcEvent;

namespace procparams
{

class ProcParams;


}

}

class PParamsChangeListener
{
public:
    virtual ~PParamsChangeListener() = default;
    virtual void procParamsChanged(
        const rtengine::procparams::ProcParams* params,
        const rtengine::ProcEvent& ev,
        const Glib::ustring& descr,
        const ParamsEdited* paramsEdited = nullptr
    ) = 0;
    virtual void clearParamChanges() = 0;
};

class BatchPParamsChangeListener
{
public:
    virtual ~BatchPParamsChangeListener() = default;
    virtual void beginBatchPParamsChange(int numberOfEntries) = 0;
    virtual void endBatchPParamsChange() = 0;
};
