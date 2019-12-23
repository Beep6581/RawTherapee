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
#include "filethumbnailbuttonset.h"

#include "rtimage.h"
#include "multilangmgr.h"
#include "lwbutton.h"
#include "rtsurface.h"

bool FileThumbnailButtonSet::iconsLoaded = false;

Cairo::RefPtr<RTSurface> FileThumbnailButtonSet::rankIcon;
Cairo::RefPtr<RTSurface> FileThumbnailButtonSet::gRankIcon;
Cairo::RefPtr<RTSurface> FileThumbnailButtonSet::unRankIcon;
Cairo::RefPtr<RTSurface> FileThumbnailButtonSet::trashIcon;
Cairo::RefPtr<RTSurface> FileThumbnailButtonSet::unTrashIcon;
Cairo::RefPtr<RTSurface> FileThumbnailButtonSet::processIcon;
std::array<Cairo::RefPtr<RTSurface>, 6> FileThumbnailButtonSet::colorLabelIcon;

Glib::ustring FileThumbnailButtonSet::processToolTip;
Glib::ustring FileThumbnailButtonSet::unrankToolTip;
Glib::ustring FileThumbnailButtonSet::trashToolTip;
Glib::ustring FileThumbnailButtonSet::untrashToolTip;
Glib::ustring FileThumbnailButtonSet::colorLabelToolTip;
std::array<Glib::ustring, 5> FileThumbnailButtonSet::rankToolTip;

FileThumbnailButtonSet::FileThumbnailButtonSet (FileBrowserEntry* myEntry)
{

    if (!iconsLoaded) {
        unRankIcon  = Cairo::RefPtr<RTSurface>(new RTSurface("star-hollow-narrow.png"));
        rankIcon    = Cairo::RefPtr<RTSurface>(new RTSurface("star-gold-narrow.png"));
        gRankIcon   = Cairo::RefPtr<RTSurface>(new RTSurface("star-narrow.png"));
        trashIcon   = Cairo::RefPtr<RTSurface>(new RTSurface("trash-small.png"));
        unTrashIcon = Cairo::RefPtr<RTSurface>(new RTSurface("trash-remove-small.png"));
        processIcon = Cairo::RefPtr<RTSurface>(new RTSurface("gears-small.png"));
        colorLabelIcon[0] = Cairo::RefPtr<RTSurface>(new RTSurface("circle-empty-gray-small.png"));
        colorLabelIcon[1] = Cairo::RefPtr<RTSurface>(new RTSurface("circle-red-small.png"));
        colorLabelIcon[2] = Cairo::RefPtr<RTSurface>(new RTSurface("circle-yellow-small.png"));
        colorLabelIcon[3] = Cairo::RefPtr<RTSurface>(new RTSurface("circle-green-small.png"));
        colorLabelIcon[4] = Cairo::RefPtr<RTSurface>(new RTSurface("circle-blue-small.png"));
        colorLabelIcon[5] = Cairo::RefPtr<RTSurface>(new RTSurface("circle-purple-small.png"));

        processToolTip = M("FILEBROWSER_POPUPPROCESS");
        unrankToolTip = M("FILEBROWSER_UNRANK_TOOLTIP");
        trashToolTip = M("FILEBROWSER_POPUPTRASH");
        untrashToolTip = M("FILEBROWSER_POPUPUNTRASH");
        colorLabelToolTip = M("FILEBROWSER_COLORLABEL_TOOLTIP");
        rankToolTip[0] = M("FILEBROWSER_RANK1_TOOLTIP");
        rankToolTip[1] = M("FILEBROWSER_RANK2_TOOLTIP");
        rankToolTip[2] = M("FILEBROWSER_RANK3_TOOLTIP");
        rankToolTip[3] = M("FILEBROWSER_RANK4_TOOLTIP");
        rankToolTip[4] = M("FILEBROWSER_RANK5_TOOLTIP");

        iconsLoaded = true;
    }

    add(new LWButton(processIcon, 6, myEntry, LWButton::Left, LWButton::Center, &processToolTip));
    add(new LWButton(unRankIcon, 0, myEntry, LWButton::Left, LWButton::Center, &unrankToolTip));

    for (int i = 0; i < 5; i++) {
        add(new LWButton(rankIcon, i + 1, myEntry, LWButton::Left, LWButton::Center, &rankToolTip[i]));
    }

    add(new LWButton(trashIcon, 7, myEntry, LWButton::Right, LWButton::Center, &trashToolTip));
    add(new LWButton(colorLabelIcon[0], 8, myEntry, LWButton::Right, LWButton::Center, &colorLabelToolTip));
}

void FileThumbnailButtonSet::setRank (int stars)
{

    for (int i = 1; i <= 5; i++) {
        buttons[i + 1]->setIcon(i <= stars ? rankIcon : gRankIcon);
    }
}

void FileThumbnailButtonSet::setColorLabel (int colorLabel)
{

    if (colorLabel >= 0 && colorLabel <= 5) {
        buttons[8]->setIcon(colorLabelIcon[colorLabel]);
    }
}

void FileThumbnailButtonSet::setInTrash (bool inTrash)
{

    buttons[7]->setIcon(inTrash ? unTrashIcon : trashIcon);
    buttons[7]->setToolTip(inTrash ? &untrashToolTip : &trashToolTip);
}
