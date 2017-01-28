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

#include "lockablecolorpicker.h"
#include "options.h"
#include "../rtengine/color.h"
#include "../rtengine/rt_math.h"
#include "imagearea.h"
#include "multilangmgr.h"
#include "navigator.h"

extern Options options;

LockableColorPicker::LockableColorPicker (CropWindow* cropWindow, Glib::ustring *oProfile, Glib::ustring *wProfile)
: cropWindow(cropWindow), displayedValues(ColorPickerType::RGB), position(0, 0), size(Size::S15),
  outputProfile(oProfile), workingProfile(wProfile), validity(Validity::OUTSIDE),
  r(0.f), g(0.f), b(0.f), rpreview(0.f), gpreview(0.f), bpreview(0.f), hue(0.f), sat(0.f), val(0.f), L(0.f), a(0.f), bb(0.f)
{}

void LockableColorPicker::updateBackBuffer ()
{
    int newW, newH;

    // -------------------- setting some key constants ---------------------
    constexpr float circlePadding = 3.f;  // keep this value odd
    constexpr double opacity = 0.62;
    // ---------------------------------------------------------------------

    if (validity == Validity::INSIDE) {
        Gtk::DrawingArea *iArea = cropWindow->getImageArea();

        Glib::RefPtr<Pango::Context> pangoContext = iArea->get_pango_context ();
        Pango::FontDescription fontd = pangoContext->get_font_description();
        // set font family and size
        fontd.set_family(options.CPFontFamily == "default" ? "sans" : options.CPFontFamily);
        fontd.set_size((options.CPFontFamily == "default" ? 8 : options.CPFontSize) * Pango::SCALE);
        fontd.set_weight(Pango::WEIGHT_NORMAL);
        pangoContext->set_font_description (fontd);

        Glib::RefPtr<Pango::Layout> layout[3][2];
        Glib::ustring s1, s2, s3;
        PointerMotionListener* navigator = cropWindow->getPointerMotionListener ();

        switch (displayedValues) {
        case ColorPickerType::RGB:
            navigator->getRGBText ((int)(r*255.f), (int)(g*255.f), (int)(b*255.f), s1, s2, s3);
            layout[0][0] = iArea->create_pango_layout(M("NAVIGATOR_R"));
            layout[0][1] = iArea->create_pango_layout(s1);
            layout[1][0] = iArea->create_pango_layout(M("NAVIGATOR_G"));
            layout[1][1] = iArea->create_pango_layout(s2);
            layout[2][0] = iArea->create_pango_layout(M("NAVIGATOR_B"));
            layout[2][1] = iArea->create_pango_layout(s3);
            break;
        case ColorPickerType::HSV:
            navigator->getHSVText (hue, sat, val, s1, s2, s3);
            layout[0][0] = iArea->create_pango_layout(M("NAVIGATOR_H"));
            layout[0][1] = iArea->create_pango_layout(s1);
            layout[1][0] = iArea->create_pango_layout(M("NAVIGATOR_S"));
            layout[1][1] = iArea->create_pango_layout(s2);
            layout[2][0] = iArea->create_pango_layout(M("NAVIGATOR_V"));
            layout[2][1] = iArea->create_pango_layout(s3);
            break;
        case ColorPickerType::LAB:
        default:
            navigator->getLABText (L, a, bb, s1, s2, s3);
            layout[0][0] = iArea->create_pango_layout(M("NAVIGATOR_LAB_L"));
            layout[0][1] = iArea->create_pango_layout(s1);
            layout[1][0] = iArea->create_pango_layout(M("NAVIGATOR_LAB_A"));
            layout[1][1] = iArea->create_pango_layout(s2);
            layout[2][0] = iArea->create_pango_layout(M("NAVIGATOR_LAB_B"));
            layout[2][1] = iArea->create_pango_layout(s3);
        }

        int w00, w01, w10, w11, w20, w21, h00, h01, h10, h11, h20, h21;
        layout[0][0]->get_pixel_size(w00, h00);
        layout[1][0]->get_pixel_size(w10, h10);
        layout[2][0]->get_pixel_size(w20, h20);
        layout[0][1]->get_pixel_size(w01, h01);
        layout[1][1]->get_pixel_size(w11, h11);
        layout[2][1]->get_pixel_size(w21, h21);
        int maxWCol0 = rtengine::max(w00, w10, w20);
        int maxWCol1 = rtengine::max(w01, w11, w21);
        int maxHRow0 = rtengine::max(h00, h01);
        int maxHRow1 = rtengine::max(h10, h11);
        int maxHRow2 = rtengine::max(h20, h21);

        // -------------------- setting some key constants ---------------------
        constexpr int textPadding = 3;
        const int textWidth = maxWCol0 + maxWCol1 + textPadding;
        const int textHeight = maxHRow0 + maxHRow1 + maxHRow2 + 2*textPadding;
        // ---------------------------------------------------------------------

        newW = rtengine::max<int>((int)size + 2 * circlePadding, textWidth + 2 * textPadding);
        newH = (int)size + 2 * circlePadding + textHeight + 2 * textPadding;

        setDrawRectangle(Cairo::FORMAT_ARGB32, 0, 0, newW, newH, true);

        Cairo::RefPtr<Cairo::Context> bbcr = BackBuffer::getContext();

        // cleaning the back buffer
        bbcr->set_source_rgba (0., 0., 0., 0.);
        bbcr->set_operator (Cairo::OPERATOR_CLEAR);
        bbcr->paint ();
        bbcr->set_operator (Cairo::OPERATOR_OVER);

        bbcr->set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
        bbcr->set_line_width (0.);

        float center = (float)size / 2.f + circlePadding;

        // black background of the whole color picker
        bbcr->set_line_width (0.);
        bbcr->set_source_rgba (0., 0., 0., opacity);
        bbcr->arc_negative (center, center, center, 0., (double)rtengine::RT_PI);
        bbcr->line_to (0, 2. * center + textHeight);
        bbcr->arc_negative (2. * textPadding, 2. * center + textHeight, 2. * textPadding, (double)rtengine::RT_PI, (double)rtengine::RT_PI / 2.);
        bbcr->line_to (textWidth, 2. * center + textHeight + 2. * textPadding);
        bbcr->arc_negative (textWidth, 2. * center + textHeight, 2. * textPadding, (double)rtengine::RT_PI / 2., 0.);
        bbcr->line_to (textWidth + 2. * textPadding, 2. * center + 2. * textPadding);
        bbcr->arc_negative (textWidth, 2. * center + 2. * textPadding, 2. * textPadding, 0., (double)rtengine::RT_PI * 1.5);
        bbcr->line_to (2. * center, 2. * center);
        bbcr->close_path();
        bbcr->set_line_join (Cairo::LINE_JOIN_BEVEL);
        bbcr->set_line_cap (Cairo::LINE_CAP_SQUARE);
        bbcr->fill ();

        // light grey circle around the color mark
        bbcr->arc (center, center, center - circlePadding / 2., 0., 2. * (double)rtengine::RT_PI);
        bbcr->set_source_rgb (0.75, 0.75, 0.75);
        bbcr->set_line_width (circlePadding - 2.);
        bbcr->stroke ();

        // spot disc with picked color
        bbcr->arc (center, center, center - circlePadding, 0., 2. * (double)rtengine::RT_PI);
        bbcr->set_source_rgb (rpreview, gpreview, bpreview);  // <- set the picker color here
        bbcr->set_line_width (0.);
        bbcr->fill();

        // adding the font
        bbcr->set_line_width (0.);
        bbcr->set_line_join (Cairo::LINE_JOIN_ROUND);
        bbcr->set_line_cap (Cairo::LINE_CAP_ROUND);
        bbcr->set_source_rgb (1., 1., 1.);
        double txtOffsetX = textPadding;
        double txtOffsetY = (double)size + 2. * circlePadding + textPadding;
        switch (iArea->get_direction()) {
        case Gtk::TEXT_DIR_RTL:
            bbcr->move_to (txtOffsetX                         , txtOffsetY);
            layout[0][1]->add_to_cairo_context (bbcr);
            bbcr->fill ();
            bbcr->move_to (txtOffsetX + maxWCol1 + textPadding, txtOffsetY);
            layout[0][0]->add_to_cairo_context (bbcr);
            bbcr->fill ();

            bbcr->move_to (txtOffsetX                         , txtOffsetY + maxHRow0 + textPadding);
            layout[1][1]->add_to_cairo_context (bbcr);
            bbcr->fill ();
            bbcr->move_to (txtOffsetX + maxWCol1 + textPadding, txtOffsetY + maxHRow0 + textPadding);
            layout[1][0]->add_to_cairo_context (bbcr);
            bbcr->fill ();

            bbcr->move_to (txtOffsetX                         , txtOffsetY + maxHRow0 + maxHRow1 + 2*textPadding);
            layout[2][1]->add_to_cairo_context (bbcr);
            bbcr->fill ();
            bbcr->move_to (txtOffsetX + maxWCol1 + textPadding, txtOffsetY + maxHRow0 + maxHRow1 + 2*textPadding);
            layout[2][0]->add_to_cairo_context (bbcr);
            bbcr->fill ();
            break;
        default:
            bbcr->move_to (txtOffsetX                         , txtOffsetY);
            layout[0][0]->add_to_cairo_context (bbcr);
            bbcr->fill ();
            bbcr->move_to (txtOffsetX + maxWCol0 + textPadding, txtOffsetY);
            layout[0][1]->add_to_cairo_context (bbcr);
            bbcr->fill ();

            bbcr->move_to (txtOffsetX                         , txtOffsetY + maxHRow0 + textPadding);
            layout[1][0]->add_to_cairo_context (bbcr);
            bbcr->fill ();
            bbcr->move_to (txtOffsetX + maxWCol0 + textPadding, txtOffsetY + maxHRow0 + textPadding);
            layout[1][1]->add_to_cairo_context (bbcr);
            bbcr->fill ();

            bbcr->move_to (txtOffsetX                         , txtOffsetY + maxHRow0 + maxHRow1 + 2*textPadding);
            layout[2][0]->add_to_cairo_context (bbcr);
            bbcr->fill ();
            bbcr->move_to (txtOffsetX + maxWCol0 + textPadding, txtOffsetY + maxHRow0 + maxHRow1 + 2*textPadding);
            layout[2][1]->add_to_cairo_context (bbcr);
            bbcr->fill ();
        }

        anchorOffset.set (center, center);

        setDirty (false);
    } else if (validity == Validity::CROSSING) {
        newH = newW = (int)size + 2 * circlePadding;

        setDrawRectangle(Cairo::FORMAT_ARGB32, 0, 0, newW, newH, true);

        Cairo::RefPtr<Cairo::Context> bbcr = BackBuffer::getContext();

        // cleaning the back buffer
        bbcr->set_source_rgba (0., 0., 0., 0.);
        bbcr->set_operator (Cairo::OPERATOR_CLEAR);
        bbcr->paint ();
        bbcr->set_operator (Cairo::OPERATOR_OVER);

        bbcr->set_antialias(Cairo::ANTIALIAS_SUBPIXEL);

        float center = (float)size / 2.f + circlePadding;

        // light grey circle around the color mark
        bbcr->arc (center, center, center - circlePadding / 2., 0., 2. * (double)rtengine::RT_PI);
        bbcr->set_source_rgba (0., 0., 0., opacity);
        bbcr->set_line_width(circlePadding);
        bbcr->stroke_preserve();
        bbcr->set_source_rgb (0.75, 0.75, 0.75);
        bbcr->set_line_width (circlePadding - 2.);
        bbcr->stroke ();

        anchorOffset.set (center, center);

        setDirty (false);
    }

}

void LockableColorPicker::draw (Cairo::RefPtr<Cairo::Context> &cr)
{
    if (validity == Validity::OUTSIDE) {
        return;
    }

    if (isDirty()) {
        updateBackBuffer();
    }
    int px, py;
    cropWindow->imageCoordToScreen(position.x, position.y, px, py);
    setDestPosition(px - anchorOffset.x, py - anchorOffset.y);

    copySurface(cr);
}

void LockableColorPicker::setPosition (const rtengine::Coord &newPos)
{
    // we're not checking bounds here, this will be done at rendering time
    position = newPos;
}

void LockableColorPicker::setRGB (const float R, const float G, const float B, const float previewR, const float previewG, const float previewB)
{
    if (r==R && g==G && b==B) {
        return;
    }

    r = R;
    g = G;
    b = B;

    rpreview = previewR;
    gpreview = previewG;
    bpreview = previewB;

    rtengine::Color::rgb2hsv(r*65535.f, g*65535.f, b*65535.f, hue, sat, val);
    rtengine::Color::rgb2lab (*outputProfile, *workingProfile, r * 65535.f, g * 65535.f, b * 65535.f, L, a, bb, options.rtSettings.HistogramWorking);  // TODO: Really sure this function works?

    if (validity != Validity::OUTSIDE) {
        setDirty(true);
    }
}

void LockableColorPicker::getImagePosition (rtengine::Coord &imgPos)
{
    imgPos = position;
}

void LockableColorPicker::getScreenPosition (rtengine::Coord &screenPos)
{
    if (cropWindow) {
        cropWindow->imageCoordToScreen(position.x, position.y, screenPos.x, screenPos.y);
    }
}

bool LockableColorPicker::isOver (int x, int y)
{
    if (!cropWindow) {
        return false;
    }
    rtengine::Coord pickerScreenPos;
    cropWindow->imageCoordToScreen(position.x, position.y, pickerScreenPos.x, pickerScreenPos.y);

    rtengine::Coord mousePos(x, y);
    rtengine::Coord wh(getWidth(), getHeight());
    rtengine::Coord tl(pickerScreenPos - anchorOffset);
    rtengine::Coord br(tl + wh);
    return mousePos >= tl && mousePos <= br;
}

void LockableColorPicker::setValidity (Validity validity)
{
    if (this->validity != validity) {
        setDirty(true);
    }
    this->validity = validity;
}

void LockableColorPicker::setSize (Size newSize)
{
    if (size != newSize)
    {
        size = newSize;
        setDirty(true);
    }
}

LockableColorPicker::Size LockableColorPicker::getSize ()
{
    return size;
}

void LockableColorPicker::rollDisplayedValues ()
{
    if (displayedValues < ColorPickerType::LAB) {
        displayedValues = (ColorPickerType)(rtengine::toUnderlying(displayedValues) + 1);
    } else {
        displayedValues = ColorPickerType::RGB;
    }
    setDirty(true);

}

bool LockableColorPicker::incSize ()
{
    if (size < Size::S30) {
        size = (Size)(rtengine::toUnderlying(size) + 5);
        setDirty(true);
        return true;
    }
    return false;
}

bool LockableColorPicker::decSize ()
{
    if (size > Size::S5) {
        size = (Size)(rtengine::toUnderlying(size) - 5);
        setDirty(true);
        return true;
    }
    return false;
}

// return true if the picker has to be redrawn
bool LockableColorPicker::cycleRGB ()
{
    if (displayedValues == ColorPickerType::RGB) {
        setDirty (true);
        return true;
    }
    return false;
}

// return true if the picker has to be redrawn
bool LockableColorPicker::cycleHSV ()
{
    if (displayedValues == ColorPickerType::HSV) {
        setDirty (true);
        return true;
    }
    return false;
}
