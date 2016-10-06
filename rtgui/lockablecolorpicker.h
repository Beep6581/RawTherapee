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

#ifndef __COLORPICKER__
#define __COLORPICKER__


#include "../rtengine/coord.h"
#include "guiutils.h"
#include "edit.h"

class CropWindow;

class LockablePickerToolListener {
public:
    virtual ~LockablePickerToolListener () {}

    /// Callback on Color Picker's visibility switch
    virtual void switchPickerVisibility (bool isVisible) {}
};

class LockableColorPicker : BackBuffer
{
public:
    enum class Size {
        S5=5,
        S10=10,
        S15=15,
        S20=20,
        S25=25,
        S30=30
    };
    enum class Validity {
        INSIDE,
        CROSSING,
        OUTSIDE
    };
private:
    enum class ColorPickerType {
        RGB,
        HSV,
        LAB
    };
    CropWindow* cropWindow;  // the color picker is displayed in a single cropWindow, the one that the user has clicked in
    ColorPickerType displayedValues;
    rtengine::Coord position;  // Coordinate in image space
    rtengine::Coord anchorOffset;
    Size size;
    Glib::ustring *outputProfile;
    Glib::ustring *workingProfile;
    Validity validity;
    float r, g, b;  // red green blue in [0;1] range
    float rpreview, gpreview, bpreview;
    float h, s, v;  // hue saturation value in [0;1] range
    float L, a, bb;  // L*a*b value in [0;1] range

    void updateBackBuffer ();

public:

    LockableColorPicker (CropWindow* cropWindow, Glib::ustring *oProfile, Glib::ustring *wProfile);
    LockableColorPicker (int x, int y, Size size, const float R, const float G, const float B, CropWindow* cropWindow, Glib::ustring *oProfile, Glib::ustring *wProfile);

    void draw (Cairo::RefPtr<Cairo::Context> &cr);

    // Used to update the RGB color, the HSV values will be updated accordingly
    void setPosition (const rtengine::Coord &newPos, const float R, const float G, float B, const float previewR, const float previewG, const float previewB);
    void setRGB (const float R, const float G, const float B, const float previewR, const float previewG, const float previewB);
    void getImagePosition (rtengine::Coord &imgPos);
    void getScreenPosition (rtengine::Coord &screenPos);
    Size getSize ();
    bool isOver (int x, int y);
    void setValidity (Validity isValid);
    void setSize (Size newSize);
    void rollDisplayedValues ();
    void incSize ();
    void decSize ();
    bool cycleRGB ();
    bool cycleHSV ();
};

#endif
