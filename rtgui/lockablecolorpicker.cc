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

extern Options options;

LockableColorPicker::LockableColorPicker (CropWindow* cropWindow)
: cropWindow(cropWindow), displayedValues(ColorPickerType::RGB), position(0, 0), size(PickerSize::S20),
  isValid(false), r(0.f), g(0.f), b(0.f), h(0.f), s(0.f), v(0.f)
{
}

LockableColorPicker::LockableColorPicker (int x, int y, PickerSize size, const float R, const float G, const float B, CropWindow* cropWindow)
: cropWindow(cropWindow), displayedValues(ColorPickerType::RGB), position(x, y), size(size),
  isValid(false), r(R), g(G), b(B)
{
    float h_, s_, v_;
    rtengine::Color::rgb2hsv((float)r*65535.f, (float)g*65535.f, (float)b*65535.f, h_, s_, v_);
    h = (int)(h_*255.f);
    s = (int)(s_*255.f);
    v = (int)(v_*255.f);
}

void LockableColorPicker::updateBackBuffer ()
{
    int newW, newH;

    // -------------------- setting some key constants ---------------------
    const float circlePadding = 3.f;
    // ---------------------------------------------------------------------

    if (isValid) {
        Gtk::DrawingArea *iArea = cropWindow->getImageArea();
        Cairo::RefPtr<Cairo::Context> bbcr = BackBuffer::getContext();

        Glib::RefPtr<Pango::Context> pangoContext = iArea->get_pango_context ();
        Pango::FontDescription fontd(options.colorPickerFont);
        pangoContext->set_font_description (fontd);

        // for drawing text
        bbcr->set_line_width (0.);

        Glib::RefPtr<Pango::Layout> layout[3][2];

        switch (displayedValues) {
        case ColorPickerType::RGB:
            layout[0][0] = iArea->create_pango_layout(M("NAVIGATOR_R") + " ");
            layout[0][1] = iArea->create_pango_layout(Glib::ustring::format(std::fixed, std::setprecision(0), (int)(r*255.f)));
            layout[1][0] = iArea->create_pango_layout(M("NAVIGATOR_G") + " ");
            layout[1][1] = iArea->create_pango_layout(Glib::ustring::format(std::fixed, std::setprecision(0), (int)(g*255.f)));
            layout[2][0] = iArea->create_pango_layout(M("NAVIGATOR_B") + " ");
            layout[2][1] = iArea->create_pango_layout(Glib::ustring::format(std::fixed, std::setprecision(0), (int)(b*255.f)));
            break;
        case ColorPickerType::HSV:
        default:
            layout[0][0] = iArea->create_pango_layout(M("NAVIGATOR_H") + " ");
            layout[0][1] = iArea->create_pango_layout(Glib::ustring::format(std::fixed, std::setprecision(0), (int)(h*255.f)));
            layout[1][0] = iArea->create_pango_layout(M("NAVIGATOR_S") + " ");
            layout[1][1] = iArea->create_pango_layout(Glib::ustring::format(std::fixed, std::setprecision(0), (int)(s*255.f)));
            layout[2][0] = iArea->create_pango_layout(M("NAVIGATOR_V") + " ");
            layout[2][1] = iArea->create_pango_layout(Glib::ustring::format(std::fixed, std::setprecision(0), (int)(v*255.f)));
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
        const int textPadding = 2;
        const int textWidth = maxWCol0 + maxWCol1 + textPadding;
        const int textHeight = maxHRow0 + maxHRow1 + maxHRow2 + 2*textPadding;
        const double opacity = 0.75;
        // ---------------------------------------------------------------------

        newW = rtengine::max<int>((int)size + 2 * circlePadding, textWidth + 2 * textPadding);
        newH = (int)size + 2 * circlePadding + textHeight + 2 * textPadding;

        setDrawRectangle(Cairo::FORMAT_ARGB32, 0, 0, newW, newH, true);

        float center = (float)size/2.f + circlePadding;

        // black background of the whole color picker
        bbcr->set_line_width (0.);
        bbcr->arc_negative (center, center, center, 0.f, (float)M_PI);
        bbcr->line_to (0, 2.f * center + textHeight + 2*textPadding);
        bbcr->line_to (textWidth + 2*textPadding, 2.f * center + textHeight + 2*textPadding);
        bbcr->line_to (textWidth + 2*textPadding, 2.f * center);
        bbcr->line_to (2.f * center, 2.f * center);
        bbcr->close_path();
        bbcr->set_source_rgba (0., 0., 0., opacity);
        bbcr->set_line_join (Cairo::LINE_JOIN_BEVEL);
        bbcr->set_line_cap (Cairo::LINE_CAP_SQUARE);
        bbcr->fill ();

        // light grey circle around the color mark
        bbcr->arc (center, center, center - circlePadding / 2.f - 0.5f, 0.f, 2.f * (float)M_PI);
        bbcr->set_source_rgb (0.7, 0.7, 0.7);
        bbcr->set_line_width (circlePadding-2.f);
        bbcr->stroke ();

        // spot disc with picked color
        bbcr->arc (center, center, center - circlePadding - 0.5f, 0.f, 2.f * (float)M_PI);
        bbcr->set_source_rgb (r, g, b);  // <- set the picker color here
        bbcr->set_line_width (0.);
        bbcr->fill();

        // adding the font
        bbcr->set_line_width (0.);
        bbcr->set_line_join (Cairo::LINE_JOIN_ROUND);
        bbcr->set_line_cap (Cairo::LINE_CAP_ROUND);
        bbcr->set_source_rgb (1., 1., 1.);
        double txtOffsetX = textPadding;
        double txtOffsetY = (double)size + 2.f * circlePadding + textPadding;
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

    } else {
        newH = newW = (int)size + 2 * circlePadding;

        setDrawRectangle(Cairo::FORMAT_ARGB32, 0, 0, newW, newH, true);

        Cairo::RefPtr<Cairo::Context> bbcr = BackBuffer::getContext();

        float center = (float)size/2.f + circlePadding;
        bbcr->arc (center, center, center - circlePadding/2.f, 0.f, 2.f * (float)M_PI);
        bbcr->set_source_rgb (0., 0., 0.);
        bbcr->set_line_width(circlePadding);
        bbcr->stroke_preserve();
        bbcr->set_source_rgb (1., 1., 1.);
        bbcr->set_line_width(circlePadding-2.f);
        bbcr->stroke ();

        anchorOffset.set (center, center);

    }

    setDirty (false);
}

void LockableColorPicker::draw (Cairo::RefPtr<Cairo::Context> &cr)
{
    if (isDirty()) {
        updateBackBuffer();
    }
    int px, py;
    cropWindow->imageCoordToScreen(position.x, position.y, px, py);
    setDestPosition(px - anchorOffset.x, py - anchorOffset.y);

    copySurface(cr);
}

void LockableColorPicker::setPosition (const rtengine::Coord &newPos, const float R, const float G, const float B)
{
    // we're not checking bounds here, this will be done at rendering time
    position = newPos;

    r = R;
    g = G;
    b = B;

    float h_, s_, v_;
    rtengine::Color::rgb2hsv((float)r*65535.f, (float)g*65535.f, (float)b*65535.f, h_, s_, v_);
    h = (float)h_;
    s = (float)s_;
    v = (float)v_;

    if (isValid) {
        setDirty(true);
    }
}

void LockableColorPicker::setRGB (const float R, const float G, const float B)
{
    if (r==R && g==G && b==B) {
        return;
    }

    r = R;
    g = G;
    b = B;

    float h_, s_, v_;
    rtengine::Color::rgb2hsv((float)r*65535.f, (float)g*65535.f, (float)b*65535.f, h_, s_, v_);
    h = (float)h_;
    s = (float)s_;
    v = (float)v_;

    if (isValid) {
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

void LockableColorPicker::setValidity (bool isValid)
{
    if (this->isValid != isValid) {
        setDirty(true);
    }
    this->isValid = isValid;
}

void LockableColorPicker::setSize (PickerSize newSize)
{
    if (size != newSize)
    {
        size = newSize;
        setDirty(true);
    }
}

LockableColorPicker::PickerSize LockableColorPicker::getSize ()
{
    return size;
}

void LockableColorPicker::incSize ()
{
    if (size < PickerSize::S30) {
        size = (PickerSize)((int)size + 5);
        setDirty(true);
    }
}

void LockableColorPicker::decSize ()
{
    if (size > PickerSize::S5) {
        size = (PickerSize)((int)size - 5);
        setDirty(true);
    }
}
