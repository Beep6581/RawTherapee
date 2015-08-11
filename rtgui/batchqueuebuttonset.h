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
#ifndef _BATCHQUEUEBUTTONSET_
#define _BATCHQUEUEBUTTONSET_

#include "lwbuttonset.h"
#include <gtkmm.h>

class BatchQueueEntry;
class BatchQueueButtonSet : public LWButtonSet
{

    static bool iconsLoaded;

public:
    static Cairo::RefPtr<Cairo::ImageSurface> cancelIcon;
    static Cairo::RefPtr<Cairo::ImageSurface> headIcon;
    static Cairo::RefPtr<Cairo::ImageSurface> tailIcon;

    BatchQueueButtonSet (BatchQueueEntry* myEntry);
};

#endif
