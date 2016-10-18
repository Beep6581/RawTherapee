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
#ifndef _PROFILECHANGELISTENER_
#define _PROFILECHANGELISTENER_

#include "../rtengine/rtengine.h"
#include <glibmm.h>

class ProfileChangeListener
{

public:
    virtual     ~ProfileChangeListener() {}
    virtual void profileChange  (const rtengine::procparams::PartialProfile* nparams, rtengine::ProcEvent event, const Glib::ustring& descr, const ParamsEdited* paramsEdited = nullptr) {}
    virtual void setDefaults    (rtengine::procparams::ProcParams* defparams) {}
};

class BatchProfileChangeListener
{

public:
    virtual     ~BatchProfileChangeListener() {}
    virtual void beginBatchProfileChange(int numberOfEntries) {}
    virtual void endBatchProfileChange() {}
};


#endif

