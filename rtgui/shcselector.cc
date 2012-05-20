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

#include "shcselector.h"
#include <iomanip>
#include "mycurve.h"

SHCSelector::SHCSelector() : movingPosition(-1), cl(NULL) {

    positions[0] = 0.25;
    positions[1] = 0.5;
    positions[2] = 0.75;
}

void SHCSelector::setPositions (double spos, double cpos, double hpos) {
    
    positions[0] = spos;
    positions[1] = cpos;
    positions[2] = hpos;
    
    queue_draw ();
}

void SHCSelector::getPositions (double& spos, double& cpos, double& hpos) {

    spos = positions[0];
    cpos = positions[1];
    hpos = positions[2];
}

void SHCSelector::on_realize() {

  Gtk::DrawingArea::on_realize();

  add_events(Gdk::EXPOSURE_MASK | Gdk::POINTER_MOTION_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK);
}

bool SHCSelector::on_expose_event(GdkEventExpose* event) {
    
    Cairo::RefPtr<Cairo::Context> cr = get_window()->create_cairo_context();

    int w = get_width () - RADIUS*2;
    int h = get_height ();

    wslider = h *2.0 / 5.0;

    Gdk::Color bgc = get_style()->get_bg (Gtk::STATE_NORMAL);
    Gdk::Color fgc = get_style()->get_text (Gtk::STATE_NORMAL);

    // clear bg
    cr->set_source_rgb (bgc.get_red_p(), bgc.get_green_p(), bgc.get_blue_p());
    cr->rectangle (0, 0, w, h);
    cr->fill();

    // draw gradient background
    Cairo::RefPtr< Cairo::LinearGradient > bggradient = Cairo::LinearGradient::create (0, 0, w, 0);
    bggradient->add_color_stop_rgb (0, 0, 0, 0);
    bggradient->add_color_stop_rgb (1, 1, 1, 1);

    cr->set_line_width (1.0);
    cr->set_source (bggradient);
    cr->rectangle (0.5+RADIUS, h*2.0/7.0 + 0.5, w-0.5, h*3.0/7.0-0.5);
    cr->fill_preserve();
    cr->set_source_rgb (fgc.get_red_p(), fgc.get_green_p(), fgc.get_blue_p());
    cr->stroke ();

    // draw sliders
    cr->set_line_width (1.0);
    for (int i=0; i<3; i++) {
        cr->move_to (RADIUS+w*positions[i]-wslider/2+0.5, h-0.5);
        cr->line_to (RADIUS+w*positions[i]-wslider/2+0.5, wslider/2 + 0.5);
        cr->line_to (RADIUS+w*positions[i], 0.5);
        cr->line_to (RADIUS+w*positions[i]+wslider/2-0.5, wslider/2 + 0.5);
        cr->line_to (RADIUS+w*positions[i]+wslider/2-0.5, h-0.5);
        cr->line_to (RADIUS+w*positions[i]-wslider/2+0.5, h-0.5);
        cr->set_source_rgb (bgc.get_red_p(), bgc.get_green_p(), bgc.get_blue_p());
        cr->fill_preserve ();
        cr->set_source_rgb (fgc.get_red_p(), fgc.get_green_p(), fgc.get_blue_p());
        cr->stroke ();
    }
    
    // draw text for the slider that is being moved
    Glib::RefPtr<Pango::Context> context = get_pango_context () ;
    cr->set_line_width (0.5);
    if (movingPosition >= 0) {
        int i = movingPosition;
        int offset;
        int layout_width, layout_height;
        Glib::RefPtr<Pango::Layout> layout = create_pango_layout(Glib::ustring::format(std::setprecision(2), positions[i]));
        layout->get_pixel_size(layout_width, layout_height);
        offset = positions[i] > 0.5 ? -layout_width-1-wslider/2 : 1+wslider/2;
        cr->move_to (RADIUS+w*positions[i]+offset-0.5, 0);
        cr->set_source_rgb (bgc.get_red_p(), bgc.get_green_p(), bgc.get_blue_p());
        layout->add_to_cairo_context (cr);
        cr->fill_preserve ();
        cr->stroke ();
        cr->move_to (RADIUS+w*positions[i]+offset+0.5, 1);
        layout->add_to_cairo_context (cr);
        cr->fill_preserve ();
        cr->stroke ();
        cr->set_source_rgb (fgc.get_red_p(), fgc.get_green_p(), fgc.get_blue_p());
        cr->move_to (RADIUS+w*positions[i]+offset, 0.5);
        layout->add_to_cairo_context (cr);
        cr->fill_preserve ();
        cr->stroke ();
    }
    return true;
}

bool SHCSelector::on_button_press_event (GdkEventButton* event) {
    
    // check if a slider is under the cursor
    int w = get_width ();
    movingPosition = -1;
    for (int i=0; i<3; i++)
        if (event->x > w*positions[i]-wslider/2 && event->x < w*positions[i]+wslider/2) {
            movingPosition = i;
            tmpX = event->x;
            tmpPos = positions[i];
            break;
        }
    queue_draw ();
    return true;
}

bool SHCSelector::on_button_release_event (GdkEventButton* event) {
    
    if (movingPosition >= 0) {
        movingPosition = -1;
        queue_draw ();
    }
    return true;
}

bool SHCSelector::on_motion_notify_event (GdkEventMotion* event) {
    
    if (movingPosition >= 0) {
        int w = get_width ();
        positions[movingPosition] = tmpPos + (event->x - tmpX) / w;
        if (positions[movingPosition] < 0)
            positions[movingPosition] = 0.0;
        else if (movingPosition > 0 && positions[movingPosition] < positions[movingPosition-1]+wslider/w)
            positions[movingPosition] = positions[movingPosition-1]+wslider/w;
        if (positions[movingPosition] > 1.0)
            positions[movingPosition] = 1.0;
        else if (movingPosition <3 && positions[movingPosition] > positions[movingPosition+1]-wslider/w)
            positions[movingPosition] = positions[movingPosition+1]-wslider/w;

        if (cl)
            cl->shcChanged ();
        queue_draw ();
    }
    return true;
}

void SHCSelector::styleChanged (const Glib::RefPtr<Gtk::Style>& style) {
    
    queue_draw ();
}

void SHCSelector::reset () {	//  : movingPosition(-1), cl(NULL) {
	positions[0] = 0.25;
	positions[1] = 0.5;
	positions[2] = 0.75;
    queue_draw ();
}
