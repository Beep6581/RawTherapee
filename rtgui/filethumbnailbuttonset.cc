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

#include "rtimage.h"
#include "multilangmgr.h"

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
        unRankIcon  = RTImage::createFromPng ("star-hollow-small.png");
        rankIcon    = RTImage::createFromPng ("star-gold-small.png");
        gRankIcon   = RTImage::createFromPng ("star-small.png");
        trashIcon   = RTImage::createFromPng ("trash-small.png");
        unTrashIcon = RTImage::createFromPng ("trash-remove-small.png");
        processIcon = RTImage::createFromPng ("gears-small.png");

        colorLabelIcon_0 = RTImage::createFromPng ("circle-empty-gray-small.png"); //("nocolorlabel.png");
        colorLabelIcon_1 = RTImage::createFromPng ("circle-red-small.png");
        colorLabelIcon_2 = RTImage::createFromPng ("circle-yellow-small.png");
        colorLabelIcon_3 = RTImage::createFromPng ("circle-green-small.png");
        colorLabelIcon_4 = RTImage::createFromPng ("circle-blue-small.png");
        colorLabelIcon_5 = RTImage::createFromPng ("circle-purple-small.png");
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
