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
#include "batchqueuebuttonset.h"

#include "lwbutton.h"
#include "multilangmgr.h"
#include "rtimage.h"
#include "rtsurface.h"

bool BatchQueueButtonSet::iconsLoaded = false;

Cairo::RefPtr<RTSurface> BatchQueueButtonSet::cancelIcon;
Cairo::RefPtr<RTSurface> BatchQueueButtonSet::headIcon;
Cairo::RefPtr<RTSurface> BatchQueueButtonSet::tailIcon;

Glib::ustring BatchQueueButtonSet::moveHeadToolTip;
Glib::ustring BatchQueueButtonSet::moveEndToolTip;
Glib::ustring BatchQueueButtonSet::cancelJobToolTip;

BatchQueueButtonSet::BatchQueueButtonSet (BatchQueueEntry* myEntry)
{

    if (!iconsLoaded) {
        cancelIcon = Cairo::RefPtr<RTSurface>(new RTSurface("cancel-small.png"));
        headIcon = Cairo::RefPtr<RTSurface>(new RTSurface("goto-start-small.png"));
        tailIcon = Cairo::RefPtr<RTSurface>(new RTSurface("goto-end-small.png"));
        moveHeadToolTip = M("FILEBROWSER_POPUPMOVEHEAD");
        moveEndToolTip = M("FILEBROWSER_POPUPMOVEEND");
        cancelJobToolTip = M("FILEBROWSER_POPUPCANCELJOB");
        iconsLoaded = true;
    }

    add(new LWButton(headIcon, 8, myEntry, LWButton::Left, LWButton::Center, &moveHeadToolTip));
    add(new LWButton(tailIcon, 9, myEntry, LWButton::Left, LWButton::Center, &moveEndToolTip));
    add(new LWButton(cancelIcon, 10, myEntry, LWButton::Right, LWButton::Center, &cancelJobToolTip));
}
