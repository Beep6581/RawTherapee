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

#include "coloredbar.h"

ColoredBar::ColoredBar (eRTOrientation orient) {
	orientation = orient;
	dirty = true;
	cp = NULL;
	this->x = this->y = this->w = this->h = 0;
}

void ColoredBar::setColorProvider (ColorProvider* p) {
	cp = p;
}

/*
 * Redraw the bar to a Cairo::Surface
 */
void ColoredBar::expose(Cairo::RefPtr<Cairo::Surface> destSurface) {
	// look out if the Surface has to be redrawn
	if (!surfaceCreated() || !destSurface)
		return;
	draw();
	copySurface(destSurface);
}

/*
 * Redraw the bar to a Gdk::Window
 */
void ColoredBar::expose(Glib::RefPtr<Gdk::Window> destWindow) {
	// look out if the Surface has to be redrawn
	if (!surfaceCreated() || !destWindow)
		return;
	draw();
	copySurface(destWindow);
}

/*
 * Redraw the bar to a Gdk::Window
 */
void ColoredBar::expose(BackBuffer *backBuffer) {
	// look out if the Surface has to be redrawn
	if (!surfaceCreated() || !backBuffer)
		return;
	draw();
	copySurface(backBuffer);
}

void ColoredBar::draw() {
	if (isDirty()) {
		Cairo::RefPtr<Cairo::Context> cr = getContext();
		// the bar has to be drawn to the Surface first
		if (!bgGradient.empty()) {
			// a gradient has been set, we use it
			cr->set_line_width(0.);

			// gradient background
			Cairo::RefPtr< Cairo::LinearGradient > bggradient;
			switch (orientation) {
				case (RTO_Left2Right):
					bggradient = Cairo::LinearGradient::create (0., 0., double(w), 0.);
					break;
				case (RTO_Right2Left):
					bggradient = Cairo::LinearGradient::create (double(w), 0., 0., 0.);
					break;
				case (RTO_Bottom2Top):
					bggradient = Cairo::LinearGradient::create (0., double(h), 0., 0.);
					break;
				case (RTO_Top2Bottom):
				default:
					bggradient = Cairo::LinearGradient::create (0., 0., 0., double(h));
					break;
			}

			for (std::vector<GradientMilestone>::iterator i=bgGradient.begin(); i!=bgGradient.end(); i++) {
				bggradient->add_color_stop_rgb (i->position, i->r, i->g, i->b);
			}
			cr->set_source (bggradient);
			cr->rectangle(0, 0, w, h);
			cr->fill();
		}
		else {
			// ask the ColorProvider to provide colors :) for each pixels
			if (cp) {
				cr->set_antialias(Cairo::ANTIALIAS_NONE);
				cr->set_line_width(1.);
				for (int x=0; x<w; x++) {
					for (int y=0; y<h; y++) {
						double x2 = double(x)+0.5;
						double y2 = double(y)+0.5;
						double x01 = x2/(w-1);
						double y01 = y2/h;
						cp->colorForValue (x01, y01);
						cr->set_source_rgb(cp->red, cp->green, cp->blue);
						cr->move_to(x2, y2);
						cr->stroke();
					}
				}
			}
		}
		// has it been updated or not, we assume that the Surface has been correctly set (we don't handle allocation error)
		setDirty(false);
	}
}

void ColoredBar::setBgGradient (const std::vector<GradientMilestone> &milestones) {
	bgGradient = milestones;
	setDirty(true);
}

void ColoredBar::clearBgGradient () {
	bgGradient.clear();
	setDirty(true);
}

bool ColoredBar::canGetColors() {
	return cp!=NULL || bgGradient.size()>0;
}
