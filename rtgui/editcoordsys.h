/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Jean-Christophe FRISCH <natureh.510@gmail.com>
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

/** @file
 * See editwidgets.h for documentation
 */

 /** @brief Coordinate system where the widgets will be drawn
  *
  * The EditCoordSystem is used to define a screen and an image coordinate system.
  */
#pragma once


 class EditCoordSystem
 {
 public:
     virtual ~EditCoordSystem() {}

     /// Convert the widget's DrawingArea (i.e. preview area) coords to the edit buffer coords
     virtual void screenCoordToCropBuffer (int phyx, int phyy, int& cropx, int& cropy) = 0;
     /// Convert the widget's DrawingArea (i.e. preview area) coords to the full image coords
     virtual void screenCoordToImage (int phyx, int phyy, int& imgx, int& imgy) = 0;
     /// Convert the image coords to the widget's DrawingArea (i.e. preview area) coords
     virtual void imageCoordToScreen (int imgx, int imgy, int& phyx, int& phyy) = 0;
     /// Convert the image coords to the crop's canvas coords (full image + padding)
     virtual void imageCoordToCropCanvas (int imgx, int imgy, int& phyx, int& phyy) = 0;
     /// Convert the image coords to the edit buffer coords  (includes borders)
     virtual void imageCoordToCropBuffer (int imgx, int imgy, int& phyx, int& phyy) = 0;
     /// Convert the image coords to the displayed image coords  (no borders here)
     virtual void imageCoordToCropImage (int imgx, int imgy, int& phyx, int& phyy) = 0;
     /// Convert a size value from the preview's scale to the image's scale
     virtual int scaleValueToImage (int value) = 0;
     /// Convert a size value from the preview's scale to the image's scale
     virtual float scaleValueToImage (float value) = 0;
     /// Convert a size value from the preview's scale to the image's scale
     virtual double scaleValueToImage (double value) = 0;
     /// Convert a size value from the image's scale to the preview's scale
     virtual int scaleValueToCanvas (int value) = 0;
     /// Convert a size value from the image's scale to the preview's scale
     virtual float scaleValueToCanvas (float value) = 0;
     /// Convert a size value from the image's scale to the preview's scale
     virtual double scaleValueToCanvas (double value) = 0;
 };

