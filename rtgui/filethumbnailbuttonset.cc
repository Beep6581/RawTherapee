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
#include "filethumbnailbuttonset.h"
#include "multilangmgr.h"
#include "../rtengine/safegtk.h"

extern Glib::ustring argv0;

bool FileThumbnailButtonSet::iconsLoaded = false;

Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::rankIcon;
Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::gRankIcon;
Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::unRankIcon;
Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::trashIcon;
Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::unTrashIcon;
Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::processIcon;
Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::colorLabelIcon_0;
Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::colorLabelIcon_1;
Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::colorLabelIcon_2;
Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::colorLabelIcon_3;
Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::colorLabelIcon_4;
Cairo::RefPtr<Cairo::ImageSurface> FileThumbnailButtonSet::colorLabelIcon_5;

FileThumbnailButtonSet::FileThumbnailButtonSet (FileBrowserEntry* myEntry)
{

    if (!iconsLoaded) {
        unRankIcon  = safe_create_from_png ("ratednotg.png");
        rankIcon    = safe_create_from_png ("rated.png");
        gRankIcon   = safe_create_from_png ("grayrated.png");
        trashIcon   = safe_create_from_png ("trash-thumbnail.png");
        unTrashIcon = safe_create_from_png ("undelete-thumbnail.png");
        processIcon = safe_create_from_png ("processing-thumbnail.png");

        colorLabelIcon_0 = safe_create_from_png ("cglabel0.png"); //("nocolorlabel.png");
        colorLabelIcon_1 = safe_create_from_png ("clabel1.png");
        colorLabelIcon_2 = safe_create_from_png ("clabel2.png");
        colorLabelIcon_3 = safe_create_from_png ("clabel3.png");
        colorLabelIcon_4 = safe_create_from_png ("clabel4.png");
        colorLabelIcon_5 = safe_create_from_png ("clabel5.png");
        iconsLoaded = true;
    }

    add (new LWButton (processIcon, 6, myEntry, LWButton::Left, LWButton::Center, M("FILEBROWSER_POPUPPROCESS")));
    add (new LWButton (unRankIcon, 0, myEntry, LWButton::Left, LWButton::Center, M("FILEBROWSER_UNRANK_TOOLTIP")));

    for (int i = 0; i < 5; i++) {
        add (new LWButton (rankIcon, i + 1, myEntry, LWButton::Left));
    }

    add (new LWButton (trashIcon, 7, myEntry, LWButton::Right, LWButton::Center, M("FILEBROWSER_POPUPTRASH")));

    add (new LWButton (colorLabelIcon_0, 8, myEntry, LWButton::Right, LWButton::Center, M("FILEBROWSER_COLORLABEL_TOOLTIP")));

    buttons[2]->setToolTip (M("FILEBROWSER_RANK1_TOOLTIP"));
    buttons[3]->setToolTip (M("FILEBROWSER_RANK2_TOOLTIP"));
    buttons[4]->setToolTip (M("FILEBROWSER_RANK3_TOOLTIP"));
    buttons[5]->setToolTip (M("FILEBROWSER_RANK4_TOOLTIP"));
    buttons[6]->setToolTip (M("FILEBROWSER_RANK5_TOOLTIP"));
}

void FileThumbnailButtonSet::setRank (int stars)
{

    for (int i = 1; i <= 5; i++) {
        buttons[i + 1]->setIcon (i <= stars ? rankIcon : gRankIcon);
    }
}

void FileThumbnailButtonSet::setColorLabel (int colorLabel)
{

    if (colorLabel == 0) {
        buttons[8]->setIcon (colorLabelIcon_0);    //transparent label
    }

    if (colorLabel == 1) {
        buttons[8]->setIcon (colorLabelIcon_1);
    }

    if (colorLabel == 2) {
        buttons[8]->setIcon (colorLabelIcon_2);
    }

    if (colorLabel == 3) {
        buttons[8]->setIcon (colorLabelIcon_3);
    }

    if (colorLabel == 4) {
        buttons[8]->setIcon (colorLabelIcon_4);
    }

    if (colorLabel == 5) {
        buttons[8]->setIcon (colorLabelIcon_5);
    }
}

void FileThumbnailButtonSet::setInTrash (bool inTrash)
{

    buttons[7]->setIcon (inTrash ? unTrashIcon : trashIcon);
    buttons[7]->setToolTip (inTrash ? M("FILEBROWSER_POPUPUNTRASH") : M("FILEBROWSER_POPUPTRASH"));
}
