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
#ifndef _FILETHUMBNAILBUTTONSET_
#define _FILETHUMBNAILBUTTONSET_

#include "lwbuttonset.h"
#include <gtkmm.h>
#include "filebrowserentry.h"

class FileBrowserEntry;
class FileThumbnailButtonSet : public LWButtonSet
{

    static bool iconsLoaded;

public:
    static Cairo::RefPtr<Cairo::ImageSurface> rankIcon;
    static Cairo::RefPtr<Cairo::ImageSurface> gRankIcon;
    static Cairo::RefPtr<Cairo::ImageSurface> unRankIcon;
    static Cairo::RefPtr<Cairo::ImageSurface> trashIcon;
    static Cairo::RefPtr<Cairo::ImageSurface> unTrashIcon;
    static Cairo::RefPtr<Cairo::ImageSurface> processIcon;

    static Cairo::RefPtr<Cairo::ImageSurface> colorLabelIcon_0;
    static Cairo::RefPtr<Cairo::ImageSurface> colorLabelIcon_1;
    static Cairo::RefPtr<Cairo::ImageSurface> colorLabelIcon_2;
    static Cairo::RefPtr<Cairo::ImageSurface> colorLabelIcon_3;
    static Cairo::RefPtr<Cairo::ImageSurface> colorLabelIcon_4;
    static Cairo::RefPtr<Cairo::ImageSurface> colorLabelIcon_5;

    explicit FileThumbnailButtonSet (FileBrowserEntry* myEntry);
    void    setRank (int stars);
    void    setColorLabel (int colorlabel);
    void    setInTrash (bool inTrash);

};

#endif
