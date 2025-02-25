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

#include <gtkmm.h>

#include "lwbuttonset.h"

class BatchQueueEntry;
class RTSurface;

class BatchQueueButtonSet : public LWButtonSet
{

    static bool iconsLoaded;

public:
    static std::shared_ptr<RTSurface> cancelIcon;
    static std::shared_ptr<RTSurface> headIcon;
    static std::shared_ptr<RTSurface> tailIcon;

    static Glib::ustring moveHeadToolTip;
    static Glib::ustring moveEndToolTip;
    static Glib::ustring cancelJobToolTip;

    explicit BatchQueueButtonSet (BatchQueueEntry* myEntry);
};
