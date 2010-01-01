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
#include <filethumbnailbuttonset.h>
#include <multilangmgr.h>

extern Glib::ustring argv0;

bool FileThumbnailButtonSet::iconsLoaded = false;

Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::rankIcon;
Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::gRankIcon;
Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::unRankIcon;
Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::trashIcon;
Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::unTrashIcon;
Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::processIcon;

FileThumbnailButtonSet::FileThumbnailButtonSet (FileBrowserEntry* myEntry) {

    if (!iconsLoaded) {
        unRankIcon  = Cairo::ImageSurface::create_from_png (argv0+"/images/unrated.png");
        rankIcon    = Cairo::ImageSurface::create_from_png (argv0+"/images/rated.png");
        gRankIcon   = Cairo::ImageSurface::create_from_png (argv0+"/images/grayrated.png");
        trashIcon   = Cairo::ImageSurface::create_from_png (argv0+"/images/trash.png");
        unTrashIcon = Cairo::ImageSurface::create_from_png (argv0+"/images/undelete.png");
        processIcon = Cairo::ImageSurface::create_from_png (argv0+"/images/processing.png");
        iconsLoaded = true;
    }

    add (new LWButton (unRankIcon, 0, myEntry, LWButton::Left, LWButton::Center, M("FILEBROWSER_POPUPUNRANK")));
    for (int i=0; i<5; i++)
        add (new LWButton (rankIcon, i+1, myEntry, LWButton::Left));
    add (new LWButton (processIcon, 6, myEntry, LWButton::Right, LWButton::Center, M("FILEBROWSER_POPUPPROCESS")));
    add (new LWButton (trashIcon, 7, myEntry, LWButton::Right, LWButton::Center, M("FILEBROWSER_POPUPTRASH")));

    buttons[1]->setToolTip (M("FILEBROWSER_POPUPRANK1"));
    buttons[2]->setToolTip (M("FILEBROWSER_POPUPRANK2"));
    buttons[3]->setToolTip (M("FILEBROWSER_POPUPRANK3"));
    buttons[4]->setToolTip (M("FILEBROWSER_POPUPRANK4"));
    buttons[5]->setToolTip (M("FILEBROWSER_POPUPRANK5"));
}

void FileThumbnailButtonSet::setRank (int stars) {

    for (int i=1; i<=5; i++)
        buttons[i]->setIcon (i<=stars ? rankIcon : gRankIcon);
}

void FileThumbnailButtonSet::setInTrash (bool inTrash) {

    buttons[7]->setIcon (inTrash ? unTrashIcon : trashIcon);
    buttons[7]->setToolTip (inTrash ? M("FILEBROWSER_POPUPUNTRASH") : M("FILEBROWSER_POPUPTRASH"));
}
