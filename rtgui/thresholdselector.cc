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

#include <cassert>
#include <cmath>

#include "thresholdselector.h"
#include "multilangmgr.h"
#include "mycurve.h"

ThresholdSelector::ThresholdSelector(double minValueBottom, double maxValueBottom, double defBottom, Glib::ustring labelBottom, unsigned int precisionBottom,
                                     double minValueTop,    double maxValueTop,    double defTop,    Glib::ustring labelTop,    unsigned int precisionTop,
                                     ThresholdCurveProvider* curveProvider)
    : coloredBar(RTO_Left2Right)
{
    positions[TS_BOTTOMLEFT]  = defPos[TS_BOTTOMLEFT]  = defBottom;
    positions[TS_TOPLEFT]     = defPos[TS_TOPLEFT]     = defTop;
    positions[TS_BOTTOMRIGHT] = defPos[TS_BOTTOMRIGHT] = 0;  // unused
    positions[TS_TOPRIGHT]    = defPos[TS_TOPRIGHT]    = 0;  // unused
    this->precisionTop = precisionTop;
    this->precisionBottom = precisionBottom;
    doubleThresh = false;

    separatedLabelBottom = labelBottom;
    separatedLabelTop = labelTop;

    bgCurveProvider = curveProvider;
    separatedSliders = true;
    initalEq1 = false;  // unused
    minValBottom = minValueBottom;
    maxValBottom = maxValueBottom;
    minValTop = minValueTop;
    maxValTop = maxValueTop;

    initValues ();
}

ThresholdSelector::ThresholdSelector(double minValue, double maxValue, double defBottom,
                                     double defTop, unsigned int precision, bool startAtOne)
    : coloredBar(RTO_Left2Right)
{
    positions[TS_BOTTOMLEFT]  = defPos[TS_BOTTOMLEFT]  = defBottom;
    positions[TS_TOPLEFT]     = defPos[TS_TOPLEFT]     = defTop;
    positions[TS_BOTTOMRIGHT] = defPos[TS_BOTTOMRIGHT] = maxValue;
    positions[TS_TOPRIGHT]    = defPos[TS_TOPRIGHT]    = maxValue;
    this->precisionTop = precision;
    this->precisionBottom = precision;
    doubleThresh = false;

    separatedLabelBottom = "";
    separatedLabelTop = "";

#ifndef NDEBUG

    if (startAtOne) {
        assert (defBottom >= defTop);
        assert (defTop >= minValue);
        assert (defBottom <= maxValue);
    } else {
        assert (defTop >= defBottom);
        assert (defBottom >= minValue);
        assert (defTop <= maxValue);
    }

    assert(minValue < maxValue);
#endif

    bgCurveProvider = nullptr;
    separatedSliders = false;
    initalEq1 = startAtOne;
    minValTop = minValBottom = minValue;
    maxValTop = maxValBottom = maxValue;

    initValues ();

}

ThresholdSelector::ThresholdSelector(double minValue, double maxValue, double defBottomLeft, double defTopLeft,
                                     double defBottomRight, double defTopRight, unsigned int precision, bool startAtOne)
    : coloredBar(RTO_Left2Right)
{
    positions[TS_BOTTOMLEFT]  = defPos[TS_BOTTOMLEFT]  = defBottomLeft;
    positions[TS_TOPLEFT]     = defPos[TS_TOPLEFT]     = defTopLeft;
    positions[TS_BOTTOMRIGHT] = defPos[TS_BOTTOMRIGHT] = defBottomRight;
    positions[TS_TOPRIGHT]    = defPos[TS_TOPRIGHT]    = defTopRight;
    this->precisionTop = precision;
    this->precisionBottom = precision;
    doubleThresh = true;

    separatedLabelBottom = "";
    separatedLabelTop = "";

#ifndef NDEBUG

    if (startAtOne) {
        assert (minValue <= defTopLeft);
        assert (defTopLeft <= defBottomLeft);
        assert (defBottomLeft <= defBottomRight);
        assert (defBottomRight <= defTopRight);
        assert (defTopRight <= maxValue);
    } else {
        assert (minValue <= defBottomLeft);
        assert (defBottomLeft <= defTopLeft);
        assert (defTopLeft <= defTopRight);
        assert (defTopRight <= defBottomRight);
        assert (defBottomRight <= maxValue);
    }

    assert(minValue < maxValue);
#endif

    bgCurveProvider = nullptr;
    separatedSliders = false;
    initalEq1 = startAtOne;
    minValTop = minValBottom = minValue;
    maxValTop = maxValBottom = maxValue;

    initValues ();
}

void ThresholdSelector::initValues ()
{

    updatePolicy = RTUP_STATIC;
    additionalTTip = "";
    oldLitCursor = litCursor = TS_UNDEFINED;
    movedCursor = TS_UNDEFINED;
    secondaryMovedCursor = TS_UNDEFINED;
    Glib::RefPtr<Gtk::StyleContext> style = get_style_context();

    style->add_class("drawingarea");
    style->add_class(GTK_STYLE_CLASS_TROUGH);
    //style->add_class(GTK_STYLE_CLASS_SCALE);
    style->add_class(GTK_STYLE_CLASS_SLIDER);

    set_name("ThresholdSelector");
    set_can_focus(false);
    set_app_paintable(true);
    setDirty(true);
    updateTooltip();
}

Gtk::SizeRequestMode ThresholdSelector::get_request_mode_vfunc () const
{
    return Gtk::SIZE_REQUEST_CONSTANT_SIZE;
}

void ThresholdSelector::get_preferred_height_vfunc (int &minimum_height, int &natural_height) const
{
    int minimumWidth = 0;
    int naturalWidth = 0;
    get_preferred_width_vfunc (minimumWidth, naturalWidth);
    get_preferred_height_for_width_vfunc (minimumWidth, minimum_height, natural_height);
}

void ThresholdSelector::get_preferred_width_vfunc (int &minimum_width, int &natural_width) const
{
    minimum_width = 60;
    natural_width = 150;
}

void ThresholdSelector::get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const
{
    natural_height = minimum_height = 23;
}

void ThresholdSelector::get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const
{
    get_preferred_width_vfunc (minimum_width, natural_width);
}

/*
 * Set the position of the sliders without telling it to the listener
 */
void ThresholdSelector::setPositions (double bottom, double top)
{

    setPositions(bottom, top, maxValBottom, maxValTop);

    if (updatePolicy == RTUP_DYNAMIC) {
        setDirty(true);
    }
}

/*
 * Set the position of the sliders without telling it to the listener
 */
void ThresholdSelector::setPositions (double bottomLeft, double topLeft, double bottomRight, double topRight)
{

    bool different = (  (positions[TS_TOPLEFT]    != topLeft)    || (positions[TS_TOPRIGHT]    != topRight)    ||
                        (positions[TS_BOTTOMLEFT] != bottomLeft) || (positions[TS_BOTTOMRIGHT] != bottomRight) );
    positions[TS_BOTTOMLEFT]  = bottomLeft;
    positions[TS_TOPLEFT]     = topLeft;
    positions[TS_BOTTOMRIGHT] = bottomRight;
    positions[TS_TOPRIGHT]    = topRight;

    if (different) {
        if (updatePolicy == RTUP_DYNAMIC) {
            setDirty(true);
        }

        sig_val_changed.emit();
        updateTooltip();
        queue_draw ();
    }
}

void ThresholdSelector::setDefaults (double bottom, double top)
{

    setDefaults(bottom, top, maxValBottom, maxValTop);
}

void ThresholdSelector::setDefaults (double bottomLeft, double topLeft, double bottomRight, double topRight)
{

    defPos[TS_BOTTOMLEFT] = bottomLeft;
    defPos[TS_TOPLEFT]    = topLeft;

    if (doubleThresh) {
        defPos[TS_BOTTOMRIGHT] = bottomRight;
        defPos[TS_TOPRIGHT]    = topRight;
    }
}

void ThresholdSelector::on_realize()
{

    Gtk::DrawingArea::on_realize();

    add_events(Gdk::POINTER_MOTION_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK | Gdk::LEAVE_NOTIFY_MASK);
}

void ThresholdSelector::updateBackBuffer()
{

    // This will create or update the size of the BackBuffer::surface
    setDrawRectangle(Cairo::FORMAT_ARGB32, 0, 0, get_width(), get_height(), true);

    if (!surface) {
        return;
    }

    Cairo::RefPtr<Cairo::Context> cr = Cairo::Context::create(surface);
    Glib::RefPtr<Gtk::StyleContext> style = get_style_context();

    cr->set_source_rgba (0., 0., 0., 0.);
    cr->set_operator (Cairo::OPERATOR_CLEAR);
    cr->paint ();
    cr->set_operator (Cairo::OPERATOR_OVER);

    double positions01[4];
    int w = get_width ();
    int h = get_height ();

    int wslider = 10;
    int hwslider = wslider / 2;
    int iw = w - wslider - 2 * hb; // inner width  (excluding padding for sliders)

    positions01[TS_BOTTOMLEFT]  = to01(TS_BOTTOMLEFT);
    positions01[TS_TOPLEFT]     = to01(TS_TOPLEFT);
    positions01[TS_BOTTOMRIGHT] = to01(TS_BOTTOMRIGHT);
    positions01[TS_TOPRIGHT]    = to01(TS_TOPRIGHT);

    // set the box's colors
    cr->set_line_width (1.0);
    cr->set_line_cap(Cairo::LINE_CAP_BUTT);

    if (is_sensitive() && coloredBar.canGetColors()) {
        if (updatePolicy == RTUP_DYNAMIC) {
            coloredBar.setDirty(true);
        }
        // this will eventually create/update the off-screen Surface for the gradient area only !
        coloredBar.setDrawRectangle(hb + hwslider, int(float(h) * 1.5f / 7.f + 0.5f), iw + 1, int(float(h) * 4.f / 7.f + 0.5f));
        // that we're displaying here
        coloredBar.expose(*this, cr);
    } else {
        style->render_background(cr, hb + hwslider, int(float(h) * 1.5f / 7.f + 0.5f), iw + 1, int(float(h) * 4.f / 7.f + 0.5f));
    }

    // draw the box's borders
    style->render_frame(cr, hb + hwslider - 0.5, double(int(float(h) * 1.5f / 7.f)) + 0.5, iw + 1, double(int(float(h) * 4.f / 7.f)));

    cr->set_line_width (1.);
    cr->set_antialias(Cairo::ANTIALIAS_NONE);

    // draw curve

    if (bgCurveProvider) {
        double yStart = double(int(float(h) * 5.5f / 7.f)) - 0.5;
        double yEnd   = double(int(float(h) * 1.5f / 7.f)) + 1.5;

        std::vector<double> pts = bgCurveProvider->getCurvePoints(this);  // the values sent by the provider are not checked (assumed to be correct)

        if (pts.size() >= 4) {
            std::vector<double>::iterator i = pts.begin();
            double x = *i;
            ++i;
            double y = *i;
            ++i;
            cr->move_to (hb + hwslider + iw * x + 0.5, (yEnd - yStart)*y + yStart);

            for (; i < pts.end(); ) {
                x = *i;
                ++i;
                y = *i;
                ++i;
                cr->line_to (hb + hwslider + iw * x + 0.5, (yEnd - yStart)*y + yStart);
            }
        } else {
            // Draw a straight line because not enough points has been sent
            double yStart = double(int(float(h) * 5.5f / 7.f)) - 0.5;
            cr->move_to (hb + hwslider + 0.5, yStart);
            cr->line_to (hb + hwslider + iw + 0.5, yStart);
        }

    } else {
        if (!separatedSliders) {
            double yStart = initalEq1 ? double(int(float(h) * 1.5f / 7.f)) + 1.5 : double(int(float(h) * 5.5f / 7.f)) - 0.5;
            double yEnd   = initalEq1 ? double(int(float(h) * 5.5f / 7.f)) - 0.5 : double(int(float(h) * 1.5f / 7.f)) + 1.5;
            ThreshCursorId p[4];

            if (initalEq1) {
                p[0] = TS_TOPLEFT;
                p[1] = TS_BOTTOMLEFT;
                p[2] = TS_BOTTOMRIGHT;
                p[3] = TS_TOPRIGHT;
            } else           {
                p[0] = TS_BOTTOMLEFT;
                p[1] = TS_TOPLEFT;
                p[2] = TS_TOPRIGHT;
                p[3] = TS_BOTTOMRIGHT;
            }

            if (positions[p[1]] > minValTop) { // we use minValTop since if this block is executed, it means that we are in a simple Threshold where both bottom and top range are the same
                cr->move_to (hb + hwslider, yStart);
            } else {
                cr->move_to (hb + hwslider, yEnd);
            }

            if (positions[p[0]] > minValTop) {
                cr->line_to (hb + hwslider + iw * positions01[p[0]] + 0.5, yStart);
            }

            if (positions[p[1]] > minValTop) {
                cr->line_to (hb + hwslider + iw * positions01[p[1]] + 0.5, yEnd);
            }

            cr->line_to (hb + hwslider + iw * positions01[p[2]] + 0.5, yEnd);

            if (doubleThresh && positions[p[2]] < maxValTop) {
                cr->line_to (hb + hwslider + iw * positions01[p[3]] + 0.5, yStart);

                if (positions[p[3]] < maxValTop) {
                    cr->line_to (hb + hwslider + iw + 0.5, yStart);
                }
            }
        }
    }

    cr->set_antialias(Cairo::ANTIALIAS_SUBPIXEL);

    if (is_sensitive() && coloredBar.canGetColors()) {
        // draw surrounding curve
        cr->set_source_rgb (0., 0., 0.);
        cr->set_line_width (5.0);
        cr->stroke_preserve();
    }

    // draw curve
    if (is_sensitive()) {
        cr->set_source_rgb (1., 1., 1.);
    } else {
        cr->set_source_rgba (0., 0., 0., 0.5);
    }

    cr->set_line_width (1.5);
    cr->stroke ();

    // draw sliders
    //if (!(litCursor == TS_UNDEFINED && movedCursor == TS_UNDEFINED)) {
    Gtk::StateFlags currState = style->get_state();

    for (int i = 0; i < (doubleThresh ? 4 : 2); ++i) {
        if (!is_sensitive()) {
            style->set_state(Gtk::STATE_FLAG_INSENSITIVE);
        } else if (i == movedCursor) {
            style->set_state(Gtk::STATE_FLAG_ACTIVE);
        } else if (i == litCursor) {
            style->set_state(Gtk::STATE_FLAG_PRELIGHT);
        } else {
            style->set_state(Gtk::STATE_FLAG_NORMAL);
        }

        double posX = hb + iw * positions01[i] + 0.5;
        double arrowY = i == 0 || i == 2 ? h - (h * 3. / 7. - 0.5) - vb : h * 3. / 7. - 0.5 + vb;
        double baseY = i == 0 || i == 2 ? h - 0.5 - vb : 0.5 + vb;

        style->render_slider(cr, posX, i == 0 || i == 2 ? arrowY : baseY, wslider, i == 0 || i == 2 ? baseY - arrowY : arrowY - baseY, Gtk::ORIENTATION_HORIZONTAL);
    }

    style->set_state(currState);

    /*
     *
    Gtk::StateFlags state = !is_sensitive() ? Gtk::STATE_FLAG_INSENSITIVE : Gtk::STATE_FLAG_NORMAL;

    cr->set_line_width (1.);
    for (int i=0; i<(doubleThresh?4:2); i++) {
        double posX = hb+hwslider+iw*positions01[i]+0.5;
        double arrowY = i==0 || i==2 ? h-(h*2.5/7.-0.5)-vb : h*2.5/7.-0.5+vb;
        double baseY = i==0 || i==2 ? h-0.5-vb : 0.5+vb;
        double centerY = (arrowY+baseY)/2.;
        cr->move_to (posX, arrowY);
        cr->line_to (posX+hwslider, centerY);
        cr->line_to (posX+hwslider, baseY);
        cr->line_to (posX-hwslider, baseY);
        cr->line_to (posX-hwslider, centerY);
        cr->close_path();
        if (i==movedCursor) {
            // moved (selected)
            c = style->get_background_color(Gtk::STATE_FLAG_SELECTED);
            cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
            cr->fill_preserve ();
            c = style->get_border_color (Gtk::STATE_FLAG_SELECTED);
            cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
            cr->stroke ();
        }
        else if (i==secondaryMovedCursor || (movedCursor==TS_UNDEFINED && i==litCursor)) {
            // prelight
            c = style->get_background_color(Gtk::STATE_FLAG_PRELIGHT);
            cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
            cr->fill_preserve ();
            c = style->get_border_color (Gtk::STATE_FLAG_PRELIGHT);
            cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
            cr->stroke ();
        }
        else {
            // normal
            c = style->get_background_color(is_sensitive() ? Gtk::STATE_FLAG_ACTIVE : Gtk::STATE_FLAG_INSENSITIVE);
            cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
            cr->fill_preserve ();
            c = style->get_border_color (is_sensitive() ? Gtk::STATE_FLAG_ACTIVE : Gtk::STATE_FLAG_INSENSITIVE);
            cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
            cr->stroke ();
        }
    }
    */
    //}
}

bool ThresholdSelector::on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr)
{

    // on_realize & updateBackBuffer have to be called before
    if (get_realized() && get_width() && get_height()) {
        if (isDirty()) {
            updateBackBuffer();
        }

        if (surface) {
            copySurface(cr);
        }
    }

    return true;
}

bool ThresholdSelector::on_button_press_event (GdkEventButton* event)
{

    if (event->button == 1)  {
        movedCursor = litCursor;
        findSecondaryMovedCursor(event->state);
        tmpX = event->x;

        queue_draw ();
    }

    grab_focus();
    return true;
}

bool ThresholdSelector::on_button_release_event (GdkEventButton* event)
{

    if (event->button == 1)  {
        findLitCursor(event->x, event->y);
        movedCursor = TS_UNDEFINED;
        secondaryMovedCursor = TS_UNDEFINED;
        queue_draw ();
    }

    return true;
}

bool ThresholdSelector::on_leave_notify_event (GdkEventCrossing* event)
{
    if (movedCursor == TS_UNDEFINED) {
        litCursor = TS_UNDEFINED;
        oldLitCursor = TS_UNDEFINED;
        queue_draw();
    }

    return true;
}

bool ThresholdSelector::on_motion_notify_event (GdkEventMotion* event)
{

    int w = get_width ();

    findLitCursor(event->x, event->y);

    if (movedCursor != TS_UNDEFINED) {
        // user is moving a cursor or two
        double minBound, maxBound, dRange;

        findSecondaryMovedCursor(event->state);

        // computing the boundaries
        findBoundaries(minBound, maxBound);

        if (movedCursor == TS_BOTTOMLEFT || movedCursor == TS_BOTTOMRIGHT) {
            dRange = maxValBottom - minValBottom;
        } else {
            dRange = maxValTop - minValTop;
        }

        double dX = ( (event->x - tmpX) * dRange ) / ( w - 2 * hb );

        // slow motion if CTRL is pressed
        if (event->state & Gdk::CONTROL_MASK) {
            dX *= 0.05;
        }

        // get the new X value, inside bounds
        double newX = positions[movedCursor] + dX;

        if (newX > maxBound) {
            newX = maxBound;
        } else if (newX < minBound) {
            newX = minBound;
        }

        // compute the effective dX
        dX = newX - positions[movedCursor];
        // set the new position of the moved cursor
        positions[movedCursor] = newX;

        // apply the decay to the secondary moved cursor, if necessary
        if (secondaryMovedCursor != TS_UNDEFINED) {
            positions[secondaryMovedCursor] += dX;
        }

        // set the new reference value for the next move
        tmpX = event->x;

        // ask to redraw the background
        if (updatePolicy == RTUP_DYNAMIC) {
            setDirty(true);
        }

        // update the tooltip
        updateTooltip();

        sig_val_changed.emit();

        queue_draw ();
    } else {
        if (litCursor != oldLitCursor) {
            queue_draw ();
        }

        oldLitCursor = litCursor;
    }


    return true;
}

void ThresholdSelector::findLitCursor(int posX, int posY)
{
    int w = get_width ();
    int h = get_height ();

    litCursor = TS_UNDEFINED;

    if (posY >= 0 && posY <= h / 2) {
        if (posX > 0 && posX < w) {
            litCursor = TS_TOPLEFT;

            if (doubleThresh) {
                // we use minValTop since if this block is executed, it means that we are in a simple Threshold where both bottom and top range are the same
                double cursorX = (posX - hb) * (maxValTop - minValTop) / (w - 2 * hb) + minValTop;

                if (cursorX > positions[TS_TOPRIGHT] || std::fabs(cursorX - positions[TS_TOPRIGHT]) < std::fabs(cursorX - positions[TS_TOPLEFT])) {
                    litCursor = TS_TOPRIGHT;
                }
            }
        }
    } else if (posY > h / 2 && posY < h) {
        if (posX > 0 && posX < w) {
            litCursor = TS_BOTTOMLEFT;

            if (doubleThresh) {
                // we use minValTop since if this block is executed, it means that we are in a simple Threshold where both bottom and top range are the same
                double cursorX = (posX - hb) * (maxValTop - minValTop) / (w - 2 * hb) + minValTop;

                if (cursorX > positions[TS_BOTTOMRIGHT] || std::fabs(cursorX - positions[TS_BOTTOMRIGHT]) < std::fabs(cursorX - positions[TS_BOTTOMLEFT])) {
                    litCursor = TS_BOTTOMRIGHT;
                }
            }
        }
    }
}

void ThresholdSelector::findBoundaries(double &min, double &max)
{

    switch (movedCursor) {
    case (TS_BOTTOMLEFT):
        if (separatedSliders) {
            if (movedCursor == TS_BOTTOMLEFT) {
                min = minValBottom;
                max = maxValBottom;
            } else if (movedCursor == TS_TOPLEFT) {
                min = minValTop;
                max = maxValTop;
            }
        } else if (initalEq1) {
            min = secondaryMovedCursor == TS_UNDEFINED ? positions[TS_TOPLEFT] : minValTop + (positions[TS_BOTTOMLEFT] - positions[TS_TOPLEFT]);
            max = positions[TS_BOTTOMRIGHT];
        } else {
            min = minValTop;
            max = secondaryMovedCursor == TS_UNDEFINED ? positions[TS_TOPLEFT] : positions[TS_TOPRIGHT] - (positions[TS_TOPLEFT] - positions[TS_BOTTOMLEFT]);
        }

        break;

    case (TS_TOPLEFT):
        if (separatedSliders) {
            if (movedCursor == TS_BOTTOMLEFT) {
                min = minValBottom;
                max = maxValBottom;
            } else if (movedCursor == TS_TOPLEFT) {
                min = minValTop;
                max = maxValTop;
            }
        } else if (initalEq1) {
            min = minValTop;
            max = secondaryMovedCursor == TS_UNDEFINED ? positions[TS_BOTTOMLEFT] : positions[TS_BOTTOMRIGHT] - (positions[TS_BOTTOMLEFT] - positions[TS_TOPLEFT]);
        } else {
            min = secondaryMovedCursor == TS_UNDEFINED ? positions[TS_BOTTOMLEFT] : minValTop + (positions[TS_TOPLEFT] - positions[TS_BOTTOMLEFT]);
            max = positions[TS_TOPRIGHT];
        }

        break;

    case (TS_BOTTOMRIGHT):
        if (initalEq1) {
            min = positions[TS_BOTTOMLEFT];
            max = secondaryMovedCursor == TS_UNDEFINED ? positions[TS_TOPRIGHT] : maxValTop - (positions[TS_TOPRIGHT] - positions[TS_BOTTOMRIGHT]);
        } else {
            min = secondaryMovedCursor == TS_UNDEFINED ? positions[TS_TOPRIGHT] : positions[TS_TOPLEFT] + (positions[TS_BOTTOMRIGHT] - positions[TS_TOPRIGHT]);
            max = maxValTop;
        }

        break;

    case (TS_TOPRIGHT):
        if (initalEq1) {
            min = secondaryMovedCursor == TS_UNDEFINED ? positions[TS_BOTTOMRIGHT] : positions[TS_BOTTOMLEFT] + (positions[TS_TOPRIGHT] - positions[TS_BOTTOMRIGHT]);
            max = maxValTop;
        } else {
            min = positions[TS_TOPLEFT];
            max = secondaryMovedCursor == TS_UNDEFINED ? positions[TS_BOTTOMRIGHT] : maxValTop - (positions[TS_BOTTOMRIGHT] - positions[TS_TOPRIGHT]);
        }

        break;

    default:
        min = minValTop;
        max = maxValTop;
        break;
    }
}

void ThresholdSelector::findSecondaryMovedCursor(guint state)
{
    secondaryMovedCursor = TS_UNDEFINED;

    if (!separatedSliders && !(state & Gdk::SHIFT_MASK)) {
        switch (movedCursor) {
        case (TS_BOTTOMLEFT):
            secondaryMovedCursor = TS_TOPLEFT;
            break;

        case (TS_TOPLEFT):
            secondaryMovedCursor = TS_BOTTOMLEFT;
            break;

        case (TS_BOTTOMRIGHT):
            secondaryMovedCursor = TS_TOPRIGHT;
            break;

        case (TS_TOPRIGHT):
            secondaryMovedCursor = TS_BOTTOMRIGHT;
            break;

        default:
            secondaryMovedCursor = TS_UNDEFINED;
            break;
        }
    }
}

void ThresholdSelector::styleChanged (const Glib::RefPtr<Gtk::Style>& style)
{

    queue_draw ();
}

void ThresholdSelector::reset ()
{

    positions[0] = defPos[0];
    positions[1] = defPos[1];
    positions[2] = defPos[2];
    positions[3] = defPos[3];

    if (updatePolicy == RTUP_DYNAMIC) {
        setDirty(true);
    }

    updateTooltip();
    queue_draw ();
}

double ThresholdSelector::to01(ThreshCursorId cursorId)
{

    double rVal;

    if (cursorId == TS_BOTTOMLEFT || cursorId == TS_BOTTOMRIGHT) {
        rVal = (positions[cursorId] - minValBottom) / (maxValBottom - minValBottom);
    } else {
        rVal = (positions[cursorId] - minValTop) / (maxValTop - minValTop);
    }

    if (rVal < 0.) {
        rVal = 0.;
    } else if (rVal > 1.) {
        rVal = 1.;
    }

    return rVal;
}

void ThresholdSelector::setBgCurveProvider (ThresholdCurveProvider* provider)
{
    bgCurveProvider = provider;
}

void ThresholdSelector::setSeparatedSliders(bool separated)
{
    separatedSliders = separated;
}

bool ThresholdSelector::getSeparatedSliders()
{
    return separatedSliders;
}

void ThresholdSelector::updateTooltip()
{

    Glib::ustring tTip;

    if (doubleThresh) {
        tTip  = Glib::ustring::compose("<b>%1:</b> %2     <b>%3:</b> %4\n<b>%5:</b> %6     <b>%7:</b> %8",
                                       M("THRESHOLDSELECTOR_TL"), Glib::ustring::format(std::fixed, std::setprecision(precisionTop), positions[TS_TOPLEFT]),
                                       M("THRESHOLDSELECTOR_TR"), Glib::ustring::format(std::fixed, std::setprecision(precisionTop), positions[TS_TOPRIGHT]),
                                       M("THRESHOLDSELECTOR_BL"), Glib::ustring::format(std::fixed, std::setprecision(precisionBottom), positions[TS_BOTTOMLEFT]),
                                       M("THRESHOLDSELECTOR_BR"), Glib::ustring::format(std::fixed, std::setprecision(precisionBottom), positions[TS_BOTTOMRIGHT])
                                      );

        if (!additionalTTip.empty()) {
            tTip += Glib::ustring::compose("\n\n%1", additionalTTip);
        }

        tTip += Glib::ustring::compose("\n\n%1", M("THRESHOLDSELECTOR_HINT"));
    } else if (separatedSliders) {
        tTip  = Glib::ustring::compose("<b>%1:</b> %2\n<b>%3:</b> %4",
                                       separatedLabelTop,    Glib::ustring::format(std::fixed, std::setprecision(precisionTop),    positions[TS_TOPLEFT]),
                                       separatedLabelBottom, Glib::ustring::format(std::fixed, std::setprecision(precisionBottom), positions[TS_BOTTOMLEFT])
                                      );

        if (!additionalTTip.empty()) {
            tTip += Glib::ustring::compose("\n\n%1", additionalTTip);
        }
    } else {
        tTip  = Glib::ustring::compose("<b>%1:</b> %2\n<b>%3:</b> %4",
                                       M("THRESHOLDSELECTOR_T"), Glib::ustring::format(std::fixed, std::setprecision(precisionTop), positions[TS_TOPLEFT]),
                                       M("THRESHOLDSELECTOR_B"), Glib::ustring::format(std::fixed, std::setprecision(precisionBottom), positions[TS_BOTTOMLEFT])
                                      );

        if (!additionalTTip.empty()) {
            tTip += Glib::ustring::compose("\n\n%1", additionalTTip);
        }

        tTip += Glib::ustring::compose("\n\n%1", M("THRESHOLDSELECTOR_HINT"));
    }

    Gtk::Widget::set_tooltip_markup(tTip);
}

sigc::signal<void> ThresholdSelector::signal_value_changed()
{
    return sig_val_changed;
}

double ThresholdSelector::shapePositionValue (ThreshCursorId cursorId)
{
    unsigned int precision = (cursorId == TS_BOTTOMLEFT || cursorId == TS_BOTTOMRIGHT) ? precisionBottom : precisionTop;
    return round(positions[cursorId] * pow(double(10), precision)) / pow(double(10), precision);
}

void ThresholdSelector::set_tooltip_markup(const Glib::ustring& markup)
{
    additionalTTip = markup;
    updateTooltip();
}

void ThresholdSelector::set_tooltip_text(const Glib::ustring& text)
{
    additionalTTip = text;
    updateTooltip();
}
