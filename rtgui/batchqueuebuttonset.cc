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
#include "batchqueuebuttonset.h"
#include "multilangmgr.h"
#include "../rtengine/safegtk.h"

extern Glib::ustring argv0;

bool BatchQueueButtonSet::iconsLoaded = false;

Cairo::RefPtr<Cairo::ImageSurface> BatchQueueButtonSet::cancelIcon;
Cairo::RefPtr<Cairo::ImageSurface> BatchQueueButtonSet::headIcon;
Cairo::RefPtr<Cairo::ImageSurface> BatchQueueButtonSet::tailIcon;

BatchQueueButtonSet::BatchQueueButtonSet (BatchQueueEntry* myEntry)
{

    if (!iconsLoaded) {
        cancelIcon = safe_create_from_png ("gtk-close.png");
        headIcon   = safe_create_from_png ("toleftend.png");
        tailIcon   = safe_create_from_png ("torightend.png");
        iconsLoaded = true;
    }

    add (new LWButton (headIcon, 8, myEntry, LWButton::Left, LWButton::Center, M("FILEBROWSER_POPUPMOVEHEAD")));
    add (new LWButton (tailIcon, 9, myEntry, LWButton::Left, LWButton::Center, M("FILEBROWSER_POPUPMOVEEND")));
    add (new LWButton (cancelIcon, 10, myEntry, LWButton::Right, LWButton::Center, M("FILEBROWSER_POPUPCANCELJOB")));
}
