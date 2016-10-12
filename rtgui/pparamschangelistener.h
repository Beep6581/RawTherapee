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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _PPARAMSCHANGELISTENER_
#define _PPARAMSCHANGELISTENER_

#include "../rtengine/rtengine.h"
#include <glibmm.h>
#include "paramsedited.h"

class PParamsChangeListener
{

public:
    virtual     ~PParamsChangeListener() {}
    virtual void procParamsChanged (rtengine::procparams::ProcParams* params, rtengine::ProcEvent ev, Glib::ustring descr, ParamsEdited* paramsEdited = nullptr) {}
    virtual void clearParamChanges () {}
};

class BatchPParamsChangeListener
{

public:
    virtual     ~BatchPParamsChangeListener() {}
    virtual void beginBatchPParamsChange(int numberOfEntries) {}
    virtual void endBatchPParamsChange() {}
};

#endif

