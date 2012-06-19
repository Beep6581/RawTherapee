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

#include "thresholdselector.h"
#include "multilangmgr.h"
#include <cassert>
#include <iomanip>
#include "mycurve.h"

ThresholdSelector::ThresholdSelector(double minValue, double maxValue, double defBottom, double defTop, unsigned int precision, bool startAtOne) {
	positions[TS_BOTTOMLEFT]  = defPos[TS_BOTTOMLEFT]  = defBottom;
	positions[TS_TOPLEFT]     = defPos[TS_TOPLEFT]     = defTop;
	positions[TS_BOTTOMRIGHT] = defPos[TS_BOTTOMRIGHT] = maxValue;
	positions[TS_TOPRIGHT]    = defPos[TS_TOPRIGHT]    = maxValue;
	this->precision = precision;
	doubleThresh = false;

#ifndef NDEBUG
	if (startAtOne) {
		assert (defBottom >= defTop);
		assert (defTop >= minValue);
		assert (defBottom <= maxValue);
	}
	else {
		assert (defTop >= defBottom);
		assert (defBottom >= minValue);
		assert (defTop <= maxValue);
	}
#endif

	initValues (minValue, maxValue, startAtOne);
}

ThresholdSelector::ThresholdSelector(double minValue, double maxValue, double defBottomLeft, double defTopLeft, double defBottomRight, double defTopRight, unsigned int precision, bool startAtOne) {
	positions[TS_BOTTOMLEFT]  = defPos[TS_BOTTOMLEFT]  = defBottomLeft;
	positions[TS_TOPLEFT]     = defPos[TS_TOPLEFT]     = defTopLeft;
	positions[TS_BOTTOMRIGHT] = defPos[TS_BOTTOMRIGHT] = defBottomRight;
	positions[TS_TOPRIGHT]    = defPos[TS_TOPRIGHT]    = defTopRight;
	this->precision = precision;
	doubleThresh = true;

#ifndef NDEBUG
	if (startAtOne) {
		assert (minValue <= defTopLeft);
		assert (defTopLeft <= defBottomLeft);
		assert (defBottomLeft <= defBottomRight);
		assert (defBottomRight <= defTopRight);
		assert (defTopRight <= maxValue);
	}
	else {
		assert (minValue <= defBottomLeft);
		assert (defBottomLeft <= defTopLeft);
		assert (defTopLeft <= defTopRight);
		assert (defTopRight <= defBottomRight);
		assert (defBottomRight <= maxValue);
	}
#endif

	initValues (minValue, maxValue, startAtOne);
}

void ThresholdSelector::initValues (double minValue, double maxValue, bool startAtOne) {
	assert(minValue <= maxValue);
	initalEq1 = startAtOne;
	minVal = minValue;
	maxVal = maxValue;
	oldLitCursor = litCursor = TS_UNDEFINED;
	movedCursor = TS_UNDEFINED;
	secondaryMovedCursor = TS_UNDEFINED;
	set_size_request (-1, 20);
	add_events(Gdk::LEAVE_NOTIFY_MASK);
	set_name("ThresholdSelector");
	set_can_focus(false);
	set_app_paintable(true);
	updateTooltip();
}

/*
 * Set the position of the sliders without telling it to the listener
 */
void ThresholdSelector::setPositions (double bottom, double top) {

	setPositions(bottom, top, maxVal, maxVal);
}

/*
 * Set the position of the sliders without telling it to the listener
 */
void ThresholdSelector::setPositions (double bottomLeft, double topLeft, double bottomRight, double topRight) {

	bool different = (  (positions[TS_TOPLEFT]    != topLeft)    || (positions[TS_TOPRIGHT]    != topRight)    ||
						(positions[TS_BOTTOMLEFT] != bottomLeft) || (positions[TS_BOTTOMRIGHT] != bottomRight) );
	positions[TS_BOTTOMLEFT]  = bottomLeft;
	positions[TS_TOPLEFT]     = topLeft;
	positions[TS_BOTTOMRIGHT] = bottomRight;
	positions[TS_TOPRIGHT]    = topRight;

	if (different) {
		sig_val_changed.emit();
		updateTooltip();
		queue_draw ();
	}
}

void ThresholdSelector::setDefaults (double bottom, double top) {

	setDefaults(bottom, top, maxVal, maxVal);
}

void ThresholdSelector::setDefaults (double bottomLeft, double topLeft, double bottomRight, double topRight) {

	defPos[TS_BOTTOMLEFT] = bottomLeft;
	defPos[TS_TOPLEFT]    = topLeft;
	if (doubleThresh) {
		defPos[TS_BOTTOMRIGHT] = bottomRight;
		defPos[TS_TOPRIGHT]    = topRight;
	}
}

void ThresholdSelector::getPositions (Glib::ustring& bottom, Glib::ustring& top) {


	bottom = Glib::ustring::format(std::fixed, std::setprecision(precision),positions[TS_BOTTOMLEFT]);
	top    = Glib::ustring::format(std::fixed, std::setprecision(precision),positions[TS_TOPLEFT]);
}

void ThresholdSelector::getPositions (Glib::ustring& bottomLeft, Glib::ustring& topLeft, Glib::ustring& bottomRight, Glib::ustring& topRight) {

	bottomLeft  = Glib::ustring::format(std::fixed, std::setprecision(precision),positions[TS_BOTTOMLEFT]);
	topLeft     = Glib::ustring::format(std::fixed, std::setprecision(precision),positions[TS_TOPLEFT]);
	bottomRight = Glib::ustring::format(std::fixed, std::setprecision(precision),positions[TS_BOTTOMRIGHT]);
	topRight    = Glib::ustring::format(std::fixed, std::setprecision(precision),positions[TS_TOPRIGHT]);
}

void ThresholdSelector::setBgGradient (const std::vector<GradientMilestone> &milestones) {
	bgGradient.clear();
	bgGradient = milestones;
}

void ThresholdSelector::on_realize() {

	Gtk::DrawingArea::on_realize();

	add_events(Gdk::EXPOSURE_MASK | Gdk::POINTER_MOTION_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK);
}

bool ThresholdSelector::on_expose_event(GdkEventExpose* event) {

	Gdk::Color c;
	Glib::RefPtr<Gdk::Window> win = get_window();
	Cairo::RefPtr<Cairo::Context> cr = win->create_cairo_context();

	double positions01[4];
	int w = get_width ();
	int h = get_height ();

	wslider = std::max(int(h / 5), 10);
	int hwslider = wslider/2;

	int iw = w-wslider-2*hb;  // inner width  (excluding padding for tabs)

	for (int i=0; i<4; i++) {
		positions01[i] = to01(positions[i]);
	}

	Gtk::StateType state = !is_sensitive() ? Gtk::STATE_INSENSITIVE : Gtk::STATE_NORMAL;
	Glib::RefPtr<Gtk::Style> style = get_style();

	// set the box's colors
	cr->set_line_width (1.0);
	cr->set_line_cap(Cairo::LINE_CAP_BUTT);
	if (is_sensitive() && bgGradient.size()>1) {
		// gradient background
		Cairo::RefPtr< Cairo::LinearGradient > bggradient = Cairo::LinearGradient::create (hwslider, 0, hwslider+iw, 0);
		for (std::vector<GradientMilestone>::iterator i=bgGradient.begin(); i!=bgGradient.end(); i++) {
			bggradient->add_color_stop_rgb (i->position, i->r, i->g, i->b);
		}
		cr->set_source (bggradient);

		// draw the box's background
		cr->rectangle (hb+hwslider-0.5, double(int(float(h)*1.5f/7.f))+0.5, iw+1, double(int(float(h)*4.f/7.f)));
		cr->fill();
	}
	else if (is_sensitive()) {
		// solid background
		c = style->get_bg (state);
		cr->set_source_rgb (c.get_red_p()*0.85, c.get_green_p()*0.85, c.get_blue_p()*0.85);

		// draw the box's background
		cr->rectangle (hb+hwslider-0.5, double(int(float(h)*1.5f/7.f))+0.5, iw+1, double(int(float(h)*4.f/7.f)));
		cr->fill();
	}

	// draw curve
	double yStart = initalEq1 ? double(int(float(h)*1.5f/7.f))+1.5 : double(int(float(h)*5.5f/7.f))-0.5;
	double yEnd   = initalEq1 ? double(int(float(h)*5.5f/7.f))-0.5 : double(int(float(h)*1.5f/7.f))+1.5;
	ThreshCursorId p[4];
	if (initalEq1) { p[0] = TS_TOPLEFT;    p[1] = TS_BOTTOMLEFT; p[2] = TS_BOTTOMRIGHT; p[3] = TS_TOPRIGHT; }
	else           { p[0] = TS_BOTTOMLEFT; p[1] = TS_TOPLEFT;    p[2] = TS_TOPRIGHT;    p[3] = TS_BOTTOMRIGHT; }
	if (positions[p[1]] > minVal)
		cr->move_to (hb+hwslider, yStart);
	else
		cr->move_to (hb+hwslider, yEnd);
	if (positions[p[0]] > minVal)
		cr->line_to (hb+hwslider+iw*positions01[p[0]]+0.5, yStart);
	if (positions[p[1]] > minVal)
		cr->line_to (hb+hwslider+iw*positions01[p[1]]+0.5, yEnd);
	cr->line_to (hb+hwslider+iw*positions01[p[2]]+0.5, yEnd);
	if (doubleThresh && positions[p[2]] < maxVal) {
		cr->line_to (hb+hwslider+iw*positions01[p[3]]+0.5, yStart);
		if (positions[p[3]] < maxVal)
			cr->line_to (hb+hwslider+iw+0.5, yStart);
	}
	if (is_sensitive() && bgGradient.size()>1) {
		// draw surrounding curve
		c = style->get_bg (state);
		cr->set_source_rgb (c.get_red_p()*0.85, c.get_green_p()*0.85, c.get_blue_p()*0.85);
		cr->set_line_width (5.0);
		cr->stroke_preserve();
	}
	// draw curve
	if (is_sensitive()) {
		c = style->get_fg (movedCursor!=TS_UNDEFINED || litCursor!=TS_UNDEFINED ? Gtk::STATE_PRELIGHT : Gtk::STATE_ACTIVE);
		cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
	}
	else {
		c = style->get_bg (Gtk::STATE_INSENSITIVE);
		cr->set_source_rgb (c.get_red_p()*0.7, c.get_green_p()*0.7, c.get_blue_p()*0.7);
	}
	cr->set_line_width (1.5);
	cr->stroke ();

	// draw the box's borders
	cr->set_line_width (1.);
	cr->rectangle (hb+hwslider-0.5, double(int(float(h)*1.5f/7.f))+0.5, iw+1, double(int(float(h)*4.f/7.f)));
    c = style->get_bg (state);
    cr->set_source_rgb (c.get_red_p()*0.7, c.get_green_p()*0.7, c.get_blue_p()*0.7);
	cr->stroke ();

	// draw sliders
	//if (!(litCursor == TS_UNDEFINED && movedCursor == TS_UNDEFINED)) {
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
			c = style->get_bg (Gtk::STATE_SELECTED);
			cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
			cr->fill_preserve ();
			//c = style->get_dark (Gtk::STATE_SELECTED);
			//cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
			c = style->get_bg (state);
			cr->set_source_rgb (c.get_red_p()*0.55, c.get_green_p()*0.55, c.get_blue_p()*0.55);
			cr->stroke ();
		}
		else if (i==secondaryMovedCursor || (movedCursor==TS_UNDEFINED && i==litCursor)) {
			// prelight
			c = style->get_bg (Gtk::STATE_PRELIGHT);
			cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
			cr->fill_preserve ();
			c = style->get_bg (state);
			cr->set_source_rgb (c.get_red_p()*0.55, c.get_green_p()*0.55, c.get_blue_p()*0.55);
			cr->stroke ();
		}
		else {
			// normal
			c = style->get_bg (is_sensitive() ? Gtk::STATE_ACTIVE : Gtk::STATE_INSENSITIVE);
			cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
			cr->fill_preserve ();
		    c = style->get_bg (state);
		    cr->set_source_rgb (c.get_red_p()*0.7, c.get_green_p()*0.7, c.get_blue_p()*0.7);
			cr->stroke ();
		}
	}
	//}
	//printf("\n\n");

	// draw text for the slider that is being moved
	/*
	 * Original code from shcselector.cc
	 *
	Glib::RefPtr<Pango::Context> context = get_pango_context () ;
	cr->set_line_width (0.5);
	if (litCursor != TS_UNDEFINED) {
		int offset;
		int layout_width, layout_height;
		Glib::RefPtr<Pango::Layout> layout = create_pango_layout(Glib::ustring::format(std::setprecision(2), positions01[litCursor]));
		layout->get_pixel_size(layout_width, layout_height);
		offset = positions01[litCursor] > 0.5 ? -layout_width-1-wslider/2 : 1+wslider/2;
		cr->move_to (w*positions01[litCursor]+offset-0.5, 0);
		cr->set_source_rgb (bgnc.get_red_p(), bgnc.get_green_p(), bgnc.get_blue_p());
		layout->add_to_cairo_context (cr);
		cr->fill_preserve ();
		cr->stroke ();
		cr->move_to (w*positions01[litCursor]+offset+0.5, 1);
		layout->add_to_cairo_context (cr);
		cr->fill_preserve ();
		cr->stroke ();
		cr->set_source_rgb (fgnc.get_red_p(), fgnc.get_green_p(), fgnc.get_blue_p());
		cr->move_to (w*positions01[litCursor]+offset, 0.5);
		layout->add_to_cairo_context (cr);
		cr->fill_preserve ();
		cr->stroke ();
	}*/
	return true;
}

bool ThresholdSelector::on_button_press_event (GdkEventButton* event) {

	if (event->button == 1)  {
		movedCursor = litCursor;
		findSecondaryMovedCursor(event->state);
		tmpX = event->x;

		queue_draw ();
	}
	grab_focus();
	return true;
}

bool ThresholdSelector::on_button_release_event (GdkEventButton* event) {

	if (event->button == 1)  {
		findLitCursor(event->x, event->y);
		movedCursor = TS_UNDEFINED;
		secondaryMovedCursor = TS_UNDEFINED;
		queue_draw ();
	}
	return true;
}

bool ThresholdSelector::on_leave_notify_event (GdkEventCrossing* event) {
	if (movedCursor == TS_UNDEFINED) {
		litCursor = TS_UNDEFINED;
		oldLitCursor = TS_UNDEFINED;
		queue_draw();
	}
	return true;
}

bool ThresholdSelector::on_motion_notify_event (GdkEventMotion* event) {

	int w = get_width ();

	findLitCursor(event->x, event->y);

	if (movedCursor != TS_UNDEFINED) {
		// user is moving a cursor or two
		double minBound, maxBound;

		findSecondaryMovedCursor(event->state);

		// computing the boundaries
		findBoundaries(minBound, maxBound);

		double dX = ( (event->x-tmpX)*(maxVal-minVal) )/( w-2*hb );
		// slow motion if CTRL is pressed
		if (event->state & Gdk::CONTROL_MASK)
			dX *= 0.05;

		// get the new X value, inside bounds
		double newX = positions[movedCursor] + dX;

		if (newX > maxBound) newX = maxBound;
		else if (newX < minBound) newX = minBound;

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

		// update the tooltip
		updateTooltip();

		sig_val_changed.emit();

		queue_draw ();
	}
	else {
		if (litCursor != oldLitCursor)
			queue_draw ();
		oldLitCursor = litCursor;
	}


	return true;
}

void ThresholdSelector::findLitCursor(int posX, int posY) {
	int w = get_width ();
	int h = get_height ();

	litCursor = TS_UNDEFINED;
	if (posY >=0 && posY <= h/2) {
		if (posX > 0 && posX < w) {
			litCursor = TS_TOPLEFT;

			if (doubleThresh) {
				double cursorX = (posX-hb)*(maxVal-minVal)/(w-2*hb)+minVal;

				if (cursorX>positions[TS_TOPRIGHT] || abs(cursorX-positions[TS_TOPRIGHT]) < abs(cursorX-positions[TS_TOPLEFT]))
					litCursor = TS_TOPRIGHT;
			}
		}
	}
	else if (posY > h/2 && posY < h) {
		if (posX > 0 && posX < w) {
			litCursor = TS_BOTTOMLEFT;
			if (doubleThresh) {
				double cursorX = (posX-hb)*(maxVal-minVal)/(w-2*hb)+minVal;

				if (cursorX>positions[TS_BOTTOMRIGHT] || abs(cursorX-positions[TS_BOTTOMRIGHT]) < abs(cursorX-positions[TS_BOTTOMLEFT]))
					litCursor = TS_BOTTOMRIGHT;
			}
		}
	}
}

void ThresholdSelector::findBoundaries(double &min, double &max) {

	switch (movedCursor) {
	case (TS_BOTTOMLEFT):
		if (initalEq1) {
			min = secondaryMovedCursor == TS_UNDEFINED ? positions[TS_TOPLEFT] : minVal+(positions[TS_BOTTOMLEFT]-positions[TS_TOPLEFT]);
			max = positions[TS_BOTTOMRIGHT];
		}
		else {
			min = minVal;
			max = secondaryMovedCursor == TS_UNDEFINED ? positions[TS_TOPLEFT] : positions[TS_TOPRIGHT]-(positions[TS_TOPLEFT]-positions[TS_BOTTOMLEFT]);
		}
		break;
	case (TS_TOPLEFT):
		if (initalEq1) {
			min = minVal;
			max = secondaryMovedCursor == TS_UNDEFINED ? positions[TS_BOTTOMLEFT] : positions[TS_BOTTOMRIGHT]-(positions[TS_BOTTOMLEFT]-positions[TS_TOPLEFT]);
		}
		else {
			min = secondaryMovedCursor == TS_UNDEFINED ? positions[TS_BOTTOMLEFT] : minVal+(positions[TS_TOPLEFT]-positions[TS_BOTTOMLEFT]);
			max = positions[TS_TOPRIGHT];
		}
		break;
	case (TS_BOTTOMRIGHT):
		if (initalEq1) {
			min = positions[TS_BOTTOMLEFT];
			max = secondaryMovedCursor == TS_UNDEFINED ? positions[TS_TOPRIGHT] : maxVal-(positions[TS_TOPRIGHT]-positions[TS_BOTTOMRIGHT]);
		}
		else {
			min = secondaryMovedCursor == TS_UNDEFINED ? positions[TS_TOPRIGHT] : positions[TS_TOPLEFT]+(positions[TS_BOTTOMRIGHT]-positions[TS_TOPRIGHT]);
			max = maxVal;
		}
		break;
	case (TS_TOPRIGHT):
		if (initalEq1) {
			min = secondaryMovedCursor == TS_UNDEFINED ? positions[TS_BOTTOMRIGHT] : positions[TS_BOTTOMLEFT]+(positions[TS_TOPRIGHT]-positions[TS_BOTTOMRIGHT]);
			max = maxVal;
		}
		else {
			min = positions[TS_TOPLEFT];
			max = secondaryMovedCursor == TS_UNDEFINED ? positions[TS_BOTTOMRIGHT] : maxVal-(positions[TS_BOTTOMRIGHT]-positions[TS_TOPRIGHT]);
		}
		break;
	default:
		min = minVal;
		max = maxVal;
		break;
	}
}

void ThresholdSelector::findSecondaryMovedCursor(guint state) {
	secondaryMovedCursor = TS_UNDEFINED;
	if (!(state & Gdk::SHIFT_MASK)) {
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

void ThresholdSelector::styleChanged (const Glib::RefPtr<Gtk::Style>& style) {

	queue_draw ();
}

void ThresholdSelector::reset () {

	positions[0] = defPos[0];
	positions[1] = defPos[1];
	positions[2] = defPos[2];
	positions[2] = defPos[3];
	updateTooltip();
	queue_draw ();
}

inline double ThresholdSelector::to01(double value) {

	double rVal = (value-minVal)/(maxVal-minVal);
	if (rVal < 0.) rVal = 0.;
	else if (rVal > 1.) rVal = 1.;
	return rVal;
}

void ThresholdSelector::updateTooltip() {

	Glib::ustring tTip;
	if (doubleThresh)
		tTip  = Glib::ustring::compose("<b>%1:</b> %2     <b>%3:</b> %4\n<b>%5:</b> %6     <b>%7:</b> %8\n%9",
										M("THRESHOLDSELECTOR_TL"), Glib::ustring::format(std::fixed, std::setprecision(precision), positions[TS_TOPLEFT]),
										M("THRESHOLDSELECTOR_TR"), Glib::ustring::format(std::fixed, std::setprecision(precision), positions[TS_TOPRIGHT]),
										M("THRESHOLDSELECTOR_BL"), Glib::ustring::format(std::fixed, std::setprecision(precision), positions[TS_BOTTOMLEFT]),
										M("THRESHOLDSELECTOR_BR"), Glib::ustring::format(std::fixed, std::setprecision(precision), positions[TS_BOTTOMRIGHT]),
										M("THRESHOLDSELECTOR_HINT")
		);
	else
		tTip  = Glib::ustring::compose("<b>%1:</b> %2\n<b>%3:</b> %4\n%5",
										M("THRESHOLDSELECTOR_T"), Glib::ustring::format(std::fixed, std::setprecision(precision), positions[TS_TOPLEFT]),
										M("THRESHOLDSELECTOR_B"), Glib::ustring::format(std::fixed, std::setprecision(precision), positions[TS_BOTTOMLEFT]),
										M("THRESHOLDSELECTOR_HINT")
		);
	set_tooltip_markup(tTip);
}

sigc::signal<void> ThresholdSelector::signal_value_changed() {
	return sig_val_changed;
}
