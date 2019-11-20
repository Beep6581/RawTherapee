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

class CropGUIListener;
class Thumbnail;
class ToolBar;
class ImageAreaToolListener
{

public:
    virtual ~ImageAreaToolListener() = default;
    virtual void spotWBselected(int x, int y, Thumbnail* thm = nullptr) = 0;
    virtual void sharpMaskSelected(bool sharpMask) = 0;
    virtual int getSpotWBRectSize() const = 0;
    virtual void cropSelectionReady() = 0;
    virtual void rotateSelectionReady(double rotate_deg, Thumbnail* thm = nullptr) = 0;
    virtual ToolBar* getToolBar() const = 0;
    virtual CropGUIListener* startCropEditing(Thumbnail* thm = nullptr) = 0;
};
