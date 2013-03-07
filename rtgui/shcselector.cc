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
#include "multilangmgr.h"
#include <iomanip>
#include "mycurve.h"

SHCSelector::SHCSelector() : ColoredBar(RTO_Left2Right), movingPosition(-1), cl(NULL) {

    positions[0] = defaults[0] = 0.25;
    positions[1] = defaults[1] = 0.5;
    positions[2] = defaults[2] = 0.75;
	leftMargin = RADIUS;
	rightMargin = RADIUS;

    // TODO: This is a hack :) ; change this name to a specific one and create a new entry in all gtkrc theme files
	set_name("ThresholdSelector");
	set_can_focus(false);
	set_size_request (-1, 12);
	set_tooltip_text(M("SHCSELECTOR_TOOLTIP"));
}

void SHCSelector::setMargins(int left, int right) {
	leftMargin = left;
	rightMargin = right;
}

void SHCSelector::setDefaults (double spos, double cpos, double hpos) {
	defaults[0] = spos;
	defaults[1] = cpos;
	defaults[2] = hpos;
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
    
	Gdk::Color c;
    Cairo::RefPtr<Cairo::Context> cr = get_window()->create_cairo_context();

    int w = get_width () - leftMargin - rightMargin;
    int h = get_height ();

	wslider = std::max(int(h / 5), 10);
	double hwslider = double(wslider)/2.;

	Gtk::StateType state = !is_sensitive() ? Gtk::STATE_INSENSITIVE : Gtk::STATE_NORMAL;
	Glib::RefPtr<Gtk::Style> style = get_style();

    // clear bg

	// set the box's colors
	cr->set_line_width (1.0);
	cr->set_line_cap(Cairo::LINE_CAP_BUTT);
	if (is_sensitive() && canGetColors()) {
		// gradient background
		Glib::RefPtr<Gdk::Window> win = get_window();
		// this will eventually create/update the off-screen pixmap
		setDrawRectangle(win, leftMargin+1, 1, w-2, int(float(h)*5.5f/7.f+0.5f));
		// that we're displaying here
		ColoredBar::expose(win);
	}
	else {
		// solid background
		c = style->get_bg (state);
		cr->set_source_rgb (c.get_red_p()*0.85, c.get_green_p()*0.85, c.get_blue_p()*0.85);

		// draw the box's background
		cr->rectangle (leftMargin+1, 1, w-2, int(float(h)*5.5f/7.f+0.5f));
		cr->fill();
	}

	// draw the box's borders
	cr->set_line_width (1.);
	cr->rectangle (leftMargin+0.5, 0.5, w-1, int(float(h)*5.5f/7.f+0.5f)+1);
    c = style->get_bg (state);
    cr->set_source_rgb (c.get_red_p()*0.7, c.get_green_p()*0.7, c.get_blue_p()*0.7);
	cr->stroke ();

    // draw sliders
    //cr->set_line_width (1.0);
    for (int i=0; i<3; i++) {
        cr->move_to (leftMargin+0.5+(w-1)*positions[i]+hwslider, double(h)-0.5);
        cr->rel_line_to (0., double(-h/3));
        cr->rel_line_to (-hwslider, double(-h/3));
        cr->rel_line_to (-hwslider, double(h/3));
        cr->rel_line_to (0., double(h/3));
        cr->close_path();
		// normal
		c = style->get_bg (is_sensitive() ? Gtk::STATE_ACTIVE : Gtk::STATE_INSENSITIVE);
		cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
		cr->fill_preserve ();
	    c = style->get_bg (state);
	    cr->set_source_rgb (c.get_red_p()*0.7, c.get_green_p()*0.7, c.get_blue_p()*0.7);
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
        offset = positions[i] > 0.5 ? -layout_width-1-hwslider : 1+hwslider;
		c = style->get_bg (state);
		cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());

		cr->set_line_width(3.);
		cr->set_line_join(Cairo::LINE_JOIN_ROUND);
		cr->set_line_cap(Cairo::LINE_CAP_ROUND);

        cr->move_to (leftMargin+w*positions[i]+offset, 0.);
		layout->add_to_cairo_context (cr);
        cr->stroke_preserve();
		c = style->get_fg (Gtk::STATE_PRELIGHT);
		cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
        cr->fill ();
    }
    return true;
}

bool SHCSelector::on_button_press_event (GdkEventButton* event) {
    
    // check if a slider is under the cursor
    double w = double(get_width ()-leftMargin-rightMargin);
    movingPosition = -1;
    for (int i=0; i<3; i++)
        if (event->x > double(leftMargin)+w*positions[i]-wslider/2. && event->x < double(leftMargin)+w*positions[i]+wslider/2) {
            movingPosition = i;
            tmpX = event->x;
            tmpPos = positions[i];
            break;
        }
    queue_draw ();
    return true;
}

bool SHCSelector::on_button_release_event (GdkEventButton* event) {
    
    if (event->button == 1) {
		if (movingPosition >= 0) {
			movingPosition = -1;
			queue_draw ();
		}
    }
    else if (event->button == 3) {
		if (movingPosition >= 0)
			movingPosition = -1;
    	// right mouse button reset the selector to the stored default values
    	if (reset()) {
    		// rest has modified the values
            if (cl)
                cl->shcChanged ();
    	}
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

bool SHCSelector::reset () {	//  : movingPosition(-1), cl(NULL) {
	if ( positions[0] != defaults[0] ||
	     positions[1] != defaults[1] ||
	     positions[2] != defaults[2]
	) {

		positions[0] = defaults[0];
		positions[1] = defaults[1];
		positions[2] = defaults[2];
		queue_draw ();
		return true;
	}
	return false;
}

void SHCSelector::refresh() {
	setDirty(true);
	Glib::RefPtr<Gdk::Window> win = get_window();
	if (win)
		win->invalidate(true);
}
