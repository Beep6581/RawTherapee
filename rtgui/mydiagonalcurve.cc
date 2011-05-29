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
#include <mydiagonalcurve.h>
#include <curves.h>
#include <string.h>
#include <gdkmm/types.h>

MyDiagonalCurve::MyDiagonalCurve () : activeParam(-1), bghistvalid(false) {

    innerWidth = get_allocation().get_width() - RADIUS * 2;
    innerHeight = get_allocation().get_height() - RADIUS * 2;
    prevInnerHeight = innerHeight;
    grab_point = -1;
    lit_point = -1;
    buttonPressed = false;

    bghist = new unsigned int[256];

    signal_event().connect( sigc::mem_fun(*this, &MyDiagonalCurve::handleEvents) );

    curve.type = DCT_Spline;

    curve.x.push_back(0.);
    curve.y.push_back(0.);
    curve.x.push_back(1.);
    curve.y.push_back(1.);
}

MyDiagonalCurve::~MyDiagonalCurve () {

    delete [] bghist;
}

std::vector<double> MyDiagonalCurve::get_vector (int veclen) {

    std::vector<double> vector;
    vector.resize (veclen);

    if (curve.type != DCT_Parametric) {
        // count active points:
        double prev =- 1.0;
        int active = 0;
        int firstact = -1;
        for (int i = 0; i < (int)curve.x.size(); ++i)
            if (curve.x[i] > prev) {
                if (firstact < 0)
                  firstact = i;
                prev = curve.x[i];
                ++active;
            }
        // handle degenerate case:
        if (active < 2) {
            double ry;
            if (active > 0)
                ry = curve.y[firstact];
            else
                ry = 0.0;
            if (ry < 0.0) ry = 0.0;
            if (ry > 1.0) ry = 1.0;
            for (int x = 0; x < veclen; ++x)
                vector[x] = ry;
            return vector;
        }
    }

    // calculate remaining points
    std::vector<double> curveDescr = getPoints ();
    rtengine::DiagonalCurve* rtcurve = new rtengine::DiagonalCurve (curveDescr, veclen*1.5);
    std::vector<double> t;
    t.resize (veclen);
    for (int i = 0; i < veclen; i++)
        t[i] = (double) i / (veclen - 1.0);
    rtcurve->getVal (t, vector);
    delete rtcurve;
    return vector;
}

void MyDiagonalCurve::interpolate () {

    prevInnerHeight = innerHeight;
    point.resize (innerWidth);
    std::vector<double> vector = get_vector (innerWidth);
    prevInnerHeight = innerHeight;
    for (int i = 0; i < innerWidth; ++i)
        point[i] = Gdk::Point (RADIUS + i, RADIUS + innerHeight - (int)((innerHeight-1) * vector[i] + 0.5));
    upoint.clear ();
    lpoint.clear ();

    if (curve.type==DCT_Parametric && activeParam>0) {
        double tmp = curve.x[activeParam-1];
        if (activeParam>=4) {
            upoint.resize(innerWidth);
            lpoint.resize(innerWidth);
            curve.x[activeParam-1] = 100;
            vector = get_vector (innerWidth);
            for (int i = 0; i < innerWidth; ++i)
                upoint[i] = Gdk::Point (RADIUS + i, RADIUS + innerHeight - (int)((innerHeight-1) * vector[i] + 0.5));
            curve.x[activeParam-1] = -100;
            vector = get_vector (innerWidth);
            for (int i = 0; i < innerWidth; ++i)
                lpoint[i] = Gdk::Point (RADIUS + i, RADIUS + innerHeight - (int)((innerHeight-1) * vector[i] + 0.5));
            curve.x[activeParam-1] = tmp;
        }
    }
}

void MyDiagonalCurve::draw (int handle) {
    if (!pixmap)
        return;

    // re-calculate curve if dimensions changed
    if (prevInnerHeight != innerHeight || (int)point.size() != innerWidth)
        interpolate ();

    Gtk::StateType state = Gtk::STATE_NORMAL;
    if (!is_sensitive())
        state = Gtk::STATE_INSENSITIVE;

    Glib::RefPtr<Gtk::Style> style = get_style ();
    Cairo::RefPtr<Cairo::Context> cr = pixmap->create_cairo_context();

    // bounding rectangle
    Gdk::Color c = style->get_bg (state);
    cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
    cr->rectangle (0, 0, innerWidth + RADIUS*2, innerHeight + RADIUS*2);
    cr->fill ();

    // histogram in the background
    if (bghistvalid) {
        // find highest bin
        unsigned int histheight = 0;
        for (int i=0; i<256; i++)
            if (bghist[i]>histheight)
	            histheight = bghist[i];
        // draw histogram
        cr->set_line_width (1.0);
        double stepSize = (innerWidth-1) / 256.0;
        cr->move_to (RADIUS, innerHeight-1+RADIUS);
        cr->set_source_rgb (0.75, 0.75, 0.75);
        for (int i=0; i<256; i++) {
            double val = bghist[i] * (double)(innerHeight-2) / (double)histheight;
            if (val>innerHeight-1)
                val = innerHeight-1;
      	    if (i>0)
       	        cr->line_to (i*stepSize+RADIUS, innerHeight-1+RADIUS-val);
    	}
        cr->line_to (innerWidth-1+RADIUS, innerHeight-1+RADIUS);
    	cr->fill ();
    }

    // draw the grid lines:
    cr->set_line_width (1.0);
    c = style->get_dark (state);
    cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
    cr->set_antialias (Cairo::ANTIALIAS_NONE);
    for (int i = 0; i < 5; i++) { // + 0.5 to align well with f(x)=x so it will cut through the center
        cr->move_to (RADIUS, MAX(0,i * (innerHeight + 0.5) / 4) + RADIUS);
        cr->line_to (innerWidth + RADIUS, MAX(0,i * (innerHeight + 0.5) / 4) + RADIUS);
        cr->move_to (MAX(0,i * innerWidth / 4) + RADIUS, RADIUS);
        cr->line_to (MAX(0,i * innerWidth / 4) + RADIUS, innerHeight + RADIUS);
    }
    cr->stroke ();

    // draw f(x)=x line
    if (snapToElmt == -2)
        cr->set_source_rgb (1.0, 0.0, 0.0);
    else
        cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
    std::valarray<double> ds (1);
    ds[0] = 4;
    cr->set_dash (ds, 0);
    cr->move_to (RADIUS, innerHeight + RADIUS);
    cr->line_to (innerWidth + RADIUS, RADIUS);
    cr->stroke ();
    cr->unset_dash ();

    cr->set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
    cr->set_line_width (1.0);

    // draw upper and lower bounds
    if (curve.type==DCT_Parametric && activeParam>0 && lpoint.size()>1 && upoint.size()>1) {
        cr->set_source_rgba (0.0, 0.0, 0.0, 0.15);
        cr->move_to (upoint[0].get_x(), upoint[0].get_y());
        for (int i=1; i<(int)upoint.size(); i++)
            cr->line_to (upoint[i].get_x(), upoint[i].get_y());
        cr->line_to (lpoint[lpoint.size()-1].get_x(), lpoint[lpoint.size()-1].get_y());
        for (int i=(int)lpoint.size()-2; i>=0; i--)
            cr->line_to (lpoint[i].get_x(), lpoint[i].get_y());
        cr->line_to (upoint[0].get_x(), upoint[0].get_y());
        cr->fill ();
    }

    // draw the cage of the NURBS curve
    if (curve.type==DCT_NURBS) {
    	unsigned int nbPoints;
        std::valarray<double> ch_ds (1);
        ch_ds[0] = 2;
        cr->set_dash (ch_ds, 0);
        cr->set_source_rgb (0.0, 0.0, 0.0);
        std::vector<double> points = getPoints();
        nbPoints = ((int)points.size()-1)/2;
        for (unsigned int i = 1; i < nbPoints; i++) {
        	int pos = i*2+1;

        	double x1 = ((innerWidth-1) * points[pos-2] + 0.5)+RADIUS;    // project (curve.x[i], 0, 1, innerWidth);
			double y1 = innerHeight - ((innerHeight-1) * points[pos-1] + 0.5)+RADIUS; // project (curve.y[i], 0, 1, innerHeight);
			double x2 = ((innerWidth-1) * points[pos] + 0.5)+RADIUS;    // project (curve.x[i], 0, 1, innerWidth);
			double y2 = innerHeight - ((innerHeight-1) * points[pos+1] + 0.5)+RADIUS; // project (curve.y[i], 0, 1, innerHeight);

			// set the color of the line when the point is snapped to the cage
			if (curve.x.size() == nbPoints && snapToElmt >= 1000 && ((i == (snapToElmt-1000)) || (i == (snapToElmt-999))))
				cr->set_source_rgb (1.0, 0.0, 0.0);
			else
				cr->set_source_rgb (0.0, 0.0, 0.0);
			cr->move_to (x1, y1);
			cr->line_to (x2, y2);
			cr->stroke ();
        }
        cr->unset_dash ();
    }

    // draw curve
    cr->set_source_rgb (0.0, 0.0, 0.0);
    cr->move_to (point[0].get_x(), point[0].get_y());
    for (int i=1; i<(int)point.size(); i++)
        cr->line_to (point[i].get_x(), point[i].get_y());
    cr->stroke ();

    // draw bullets
    if (curve.type!=DCT_Parametric)
        for (int i = 0; i < (int)curve.x.size(); ++i) {
        	if (curve.x[i] == -1) continue;
        	if (snapToElmt >= 1000) {
        		int pt = snapToElmt-1000;
        		if (i >= (pt-1) && i <= (pt+1))
        			cr->set_source_rgb(1.0, 0.0, 0.0);
        		else
        			cr->set_source_rgb(0.0, 0.0, 0.0);
        	}
        	else
        		cr->set_source_rgb ((i == handle || i == snapToElmt ? 1.0 : 0.0), 0.0, 0.0);

            double x = ((innerWidth-1) * curve.x[i] + 0.5)+RADIUS;    // project (curve.x[i], 0, 1, innerWidth);
            double y = innerHeight - ((innerHeight-1) * curve.y[i] + 0.5)+RADIUS; // project (curve.y[i], 0, 1, innerHeight);

            cr->arc (x, y, RADIUS+0.5, 0, 2*M_PI);
            cr->fill ();
        }

    get_window()->draw_drawable (style->get_fg_gc (state), pixmap, 0, 0, 0, 0, innerWidth + RADIUS * 2, innerHeight + RADIUS * 2);
}

/*void MyDiagonalCurve::graphSizeRequest (Gtk::Requisition* req) {
	req->width = getGraphMinSize();
	// The real height request should take care of the presence of the vertical
	// scroll bar and its width
	req->height = sized ? getGraphMinSize() : get_allocation().get_width();
}*/

bool MyDiagonalCurve::handleEvents (GdkEvent* event) {

	CursorShape new_type = cursor_type;
	int src, dst;
	std::vector<double>::iterator itx, ity;

	bool retval = false;
	int num = (int)curve.x.size();

	/* innerWidth and innerHeight are the size of the graph */
	innerWidth = get_allocation().get_width() - RADIUS * 2;
	innerHeight = get_allocation().get_height() - RADIUS * 2;

	double minDistanceX = (double)(MIN_DISTANCE) / (double)(innerWidth-1);
	double minDistanceY = (double)(MIN_DISTANCE) / (double)(innerHeight-1);

	if ((innerWidth < 0) || (innerHeight < 0))
		return false;

	switch (event->type) {
	case Gdk::CONFIGURE: {
		// Happen when the the window is resized
		if (sized & (RS_Pending | RS_Force)) {
			int size = get_allocation().get_width();
			set_size_request(-1, size);
			sized = RS_Done;
		}
		if (pixmap)
			pixmap.clear ();
		retval = true;
		break;
	}
	case Gdk::EXPOSE:
		if (sized & (RS_Pending | RS_Force)) {
			int size = get_allocation().get_width();
			set_size_request(-1, size);
		}
		sized = RS_Pending;
		if (!pixmap) {
			pixmap = Gdk::Pixmap::create (get_window(), get_allocation().get_width(),  get_allocation().get_height());
			interpolate ();
		}
		draw (lit_point);
		retval = true;
		break;

	case Gdk::BUTTON_PRESS:
		snapToElmt = -100;
		if (curve.type!=DCT_Parametric) {
			if (event->button.button == 1) {
				buttonPressed = true;
				add_modal_grab ();

				// get the pointer position
				getCursorPosition(event);
				findClosestPoint();

				new_type = CSMove;
				if (distanceX > minDistanceX) {
					/* insert a new control point */
					if (num > 0) {
						if (clampedX > curve.x[closest_point])
							++closest_point;
					}
					itx = curve.x.begin();
					ity = curve.y.begin();
					for (int i=0; i<closest_point; i++) { itx++; ity++; }
					curve.x.insert (itx, 0);
					curve.y.insert (ity, 0);
					num++;

					// the graph is refreshed only if a new point is created (snaped to a pixel)
					curve.x[closest_point] = clampedX;
					curve.y[closest_point] = clampedY;

					interpolate ();
					draw (closest_point);
					notifyListener ();
				}
				grab_point = closest_point;
				lit_point = closest_point;
				ugpX = curve.x[closest_point];
				ugpY = curve.y[closest_point];
			}
			if (buttonPressed) retval = true;
		}
		break;

	case Gdk::BUTTON_RELEASE:
		snapToElmt = -100;
		if (curve.type!=DCT_Parametric) {
			if (buttonPressed && event->button.button == 1) {
				buttonPressed = false;
				/*  get the pointer position  */
				getCursorPosition(event);
				findClosestPoint();

				remove_modal_grab ();
				int previous_lit_point = lit_point;
				/* delete inactive points: */
				itx = curve.x.begin();
				ity = curve.y.begin();
				for (src = dst = 0; src < num; ++src)
					if (curve.x[src] >= 0.0) {
						curve.x[dst] = curve.x[src];
						curve.y[dst] = curve.y[src];
						++dst;
						++itx;
						++ity;
					}
				if (dst < src) {
					curve.x.erase (itx, curve.x.end());
					curve.y.erase (ity, curve.y.end());
					if (!curve.x.size()) {
						curve.x.push_back (0);
						curve.y.push_back (0);
						interpolate ();
						draw (lit_point);
					}
				}
				if (distanceX <= minDistanceX) {
					new_type = CSMove;
					lit_point = closest_point;
				}
				else {
					new_type = CSPlus;
					lit_point = -1;
				}
				if (lit_point != previous_lit_point)
					draw (lit_point);
				grab_point = -1;
				retval = true;
				notifyListener ();
			}
		}
		break;

	case Gdk::LEAVE_NOTIFY:
		// Pointer can LEAVE even when dragging the point, so we don't modify the cursor in this case
		// The cursor will have to LEAVE another time after the drag...
		if (!buttonPressed)
			if (grab_point == -1) {
				new_type = CSArrow;
				lit_point = -1;
				draw (lit_point);
			}
		retval = true;
		break;

	case Gdk::MOTION_NOTIFY:
		snapToElmt = -100;
		if (curve.type == DCT_Linear || curve.type == DCT_Spline || curve.type == DCT_NURBS) {

			snapToMinDist = 10.;
			snapToVal = 0.;
			snapToElmt = -100;

			// get the pointer position
			getCursorPosition(event);

			if (grab_point == -1) {
				// there's no point currently being moved
				int previous_lit_point = lit_point;
				findClosestPoint();
				if (distanceX <= minDistanceX) {
					new_type = CSMove;
					lit_point = closest_point;
				}
				else {
					new_type = CSPlus;
					lit_point = -1;
				}
				if (lit_point != previous_lit_point)
					draw (lit_point);
			}
			else {
				// a point is being moved

				// bounds of the grabbed point
				double leftBound         = (grab_point == 0    ) ? 0. : curve.x[grab_point-1];
				double rightBound        = (grab_point == num-1) ? 1. : curve.x[grab_point+1];
				double const bottomBound = 0.;
				double const topBound    = 1.;

				double leftDeletionBound   = leftBound   - minDistanceX;
				double rightDeletionBound  = rightBound  + minDistanceX;
				double bottomDeletionBound = bottomBound - minDistanceY;
				double topDeletionBound    = topBound    + minDistanceY;

				// we memorize the previous position of the point, for optimization purpose
				double prevPosX = curve.x[grab_point];
				double prevPosY = curve.y[grab_point];

				// we memorize the previous position of the point, for optimization purpose
				ugpX += deltaX;
				ugpY += deltaY;

				// the unclamped grabbed point is brought back in the range when snapTo is active
				if (snapTo) ugpY = CLAMP(ugpY, 0.0, 1.0);

				// handling limitations along X axis
				if (ugpX >= rightDeletionBound && (grab_point > 0 && grab_point < (num-1))) {
					curve.x[grab_point] = -1.;
				}
				else if (ugpX <= leftDeletionBound && (grab_point > 0 && grab_point < (num-1))) {
					curve.x[grab_point] = -1.;
				}
				else
					// nextPosX is in bounds
					curve.x[grab_point] = CLAMP(ugpX, leftBound, rightBound);

				// Handling limitations along Y axis
				if (ugpY >= topDeletionBound && grab_point != 0 && grab_point != num-1) {
					curve.x[grab_point] = -1.;
				}
				else if (ugpY <= bottomDeletionBound && grab_point != 0 && grab_point != num-1) {
					curve.x[grab_point] = -1.;
				}
				else {
					// snapping point to specific values
					if (snapTo && curve.x[grab_point] != -1.) {
						if (grab_point > 0 && grab_point < (curve.y.size()-1)) {
							double prevX = curve.x[grab_point-1];
							double prevY = curve.y[grab_point-1];
							double nextX = curve.x[grab_point+1];
							double nextY = curve.y[grab_point+1];

							double ratio = (curve.x[grab_point]-prevX)/(nextX-prevX);
							double y = (nextY-prevY) * ratio + prevY;

							if (snapCoordinate(y, ugpY)) snapToElmt = 1000+grab_point;
						}
						if (grab_point > 0) {
							int prevP = grab_point-1;
							if (snapCoordinate(curve.y[prevP], ugpY)) snapToElmt = prevP;
						}
						if (grab_point < (curve.y.size()-1)) {
							int nextP = grab_point+1;
							if (snapCoordinate(curve.y[nextP], ugpY)) snapToElmt = nextP;
						}
						if (snapCoordinate(1.0,                 ugpY)) snapToElmt = -3;
						if (snapCoordinate(curve.x[grab_point], ugpY)) snapToElmt = -2;
						if (snapCoordinate(0.0,                 ugpY)) snapToElmt = -1;

						curve.y[grab_point] = snapToVal;
					}
					else {
						// nextPosY is in the bounds
						curve.y[grab_point] = CLAMP(ugpY, 0.0, 1.0);
					}
				}

				if (curve.x[grab_point] != prevPosX || curve.y[grab_point] != prevPosY) {
					// we recalculate the curve only if we have to
					interpolate ();
					draw (lit_point);
					notifyListener ();
				}
			}
		}

		retval = true;
		break;

	default:
		break;
	}
	if (new_type != cursor_type) {
		cursor_type = new_type;
		cursorManager.setCursor(cursor_type);
	}
	return retval;
}

void MyDiagonalCurve::getCursorPosition(GdkEvent* event) {
	int tx, ty;
	int prevCursorX, prevCursorY;
	double incrementX = 1. / (double)(innerWidth-1);
	double incrementY = 1. / (double)(innerHeight-1);

	// getting the cursor position
	switch (event->type) {
	case (Gdk::MOTION_NOTIFY) :
		if (event->motion.is_hint) {
			get_window()->get_pointer (tx, ty, mod_type);
		}
		else {
			tx = (int)event->button.x;
			ty = (int)event->button.y;
			mod_type = (Gdk::ModifierType)event->button.state;
		}
		break;
	case (Gdk::BUTTON_PRESS) :
	case (Gdk::BUTTON_RELEASE) :
		tx = (int)event->button.x;
		ty = (int)event->button.y;
		mod_type = (Gdk::ModifierType)event->button.state;
		break;
	default :
		// The cursor position is not available
		return;
		break;
	}

	if (grab_point != -1) {
		prevCursorX = cursorX;
		prevCursorY = cursorY;
    }
    cursorX = tx - RADIUS;
    cursorY = (innerHeight-1) - (ty - RADIUS);

    snapTo = ST_None;

    // update deltaX/Y if the user drags a point
    if (grab_point != -1) {
    	// set the dragging factor
    	int control_key = mod_type & GDK_CONTROL_MASK;
    	int shift_key = mod_type & GDK_SHIFT_MASK;

    	// the increment get smaller if modifier key are used, and "snap to" may be enabled
    	if (control_key) { incrementX *= 0.05; incrementY *= 0.05; }
    	if (shift_key)   { snapTo = true; }

    	deltaX = (double)(cursorX - prevCursorX) * incrementX;
    	deltaY = (double)(cursorY - prevCursorY) * incrementY;
    }
    // otherwise set the position of the new point (modifier keys has no effect here)
    else {
		double tempCursorX = cursorX * incrementX;
		double tempCursorY = cursorY * incrementY;
		clampedX = CLAMP (tempCursorX, 0., 1.);  // X position of the pointer from the origin of the graph
		clampedY = CLAMP (tempCursorY, 0., 1.); // Y position of the pointer from the origin of the graph
    }

}

void MyDiagonalCurve::findClosestPoint() {
    distanceX = 10.0;  distanceY = 10.0;
    closest_point = -1;

    if (curve.type!=DCT_Parametric) {
        for (int i = 0; i < (int)curve.x.size(); i++) {
            double dX = curve.x[i] - clampedX;
            double dY = curve.y[i] - clampedY;
            double currDistX = dX < 0. ? -dX : dX; //abs (dX);
            double currDistY = dY < 0. ? -dY : dY; //abs (dY);
            if (currDistX < distanceX) {
                distanceX = currDistX;
                distanceY = currDistY;
                closest_point = i;
            }
            else if (currDistX == distanceX && currDistY < distanceY) {
            	// there is more than 1 point for that X coordinate, we select the closest point to the cursor
				distanceY = currDistY;
				closest_point = i;
            }
        }
    }
}

std::vector<double> MyDiagonalCurve::getPoints () {
    std::vector<double> result;
    if (curve.type==DCT_Parametric) {
        result.push_back ((double)(DCT_Parametric));
        for (int i=0; i<(int)curve.x.size(); i++) {
            result.push_back (curve.x[i]);
        }
    }
    else {
    	// the first value gives the type of the curve
        if (curve.type==DCT_Linear)
            result.push_back ((double)(DCT_Linear));
        else if (curve.type==DCT_Spline)
            result.push_back ((double)(DCT_Spline));
        else if (curve.type==DCT_NURBS)
            result.push_back ((double)(DCT_NURBS));
        // then we push all the points coordinate
        for (int i=0; i<(int)curve.x.size(); i++) {
            if (curve.x[i]>=0) {
                result.push_back (curve.x[i]);
                result.push_back (curve.y[i]);
            }
        }
    }
    return result;
}

void MyDiagonalCurve::setPoints (const std::vector<double>& p) {
    int ix = 0;
    DiagonalCurveType t = (DiagonalCurveType)p[ix++];
    curve.type = t;
    if (t==DCT_Parametric) {
        curve.x.clear ();
        curve.y.clear ();
        for (int i=1; i<(int)p.size(); i++)
            curve.x.push_back (p[ix++]);
    }
    else {
        curve.x.clear ();
        curve.y.clear ();
        for (int i=0; i<(int)p.size()/2; i++) {
            curve.x.push_back (p[ix++]);
            curve.y.push_back (p[ix++]);
        }
        activeParam = -1;
    }
    pixmap.clear ();
    queue_draw ();
}

void MyDiagonalCurve::setType (DiagonalCurveType t) {

    curve.type = t;
    pixmap.clear ();
}

void MyDiagonalCurve::setActiveParam (int ac) {
    
    activeParam = ac;
    pixmap.clear ();
    queue_draw ();
}

int diagonalmchistupdate (void* data) {

    gdk_threads_enter ();

    MyCurveIdleHelper* mcih = (MyCurveIdleHelper*)data;

    if (mcih->destroyed) {
        if (mcih->pending == 1)
            delete mcih;
        else    
            mcih->pending--;
        gdk_threads_leave ();
        return 0;
    }

    mcih->clearPixmap ();
    mcih->myCurve->queue_draw ();

    mcih->pending--;
    gdk_threads_leave ();
    return 0;
}

void MyDiagonalCurve::updateBackgroundHistogram (LUTu & hist) {
	
    if (hist!=NULL) {
        //memcpy (bghist, hist, 256*sizeof(unsigned int));
        for (int i=0; i<256; i++) bghist[i]=hist[i];
        //hist = bghist;
        bghistvalid = true;
    }
    else
        bghistvalid = false;

    mcih->pending++;
    g_idle_add (diagonalmchistupdate, mcih);

}

void MyDiagonalCurve::reset() {
	innerWidth = get_allocation().get_width() - RADIUS * 2;
	innerHeight = get_allocation().get_height() - RADIUS * 2;

	switch (curve.type) {
	case DCT_Spline :
	case  DCT_NURBS :
		curve.x.clear();
		curve.y.clear();
	    curve.x.push_back(0.);
	    curve.y.push_back(0.);
	    curve.x.push_back(1.);
	    curve.y.push_back(1.);
	    grab_point = -1;
	    lit_point = -1;
        interpolate ();
		break;
	case DCT_Parametric :
		// Nothing to do (?)
	default:
		break;
	}
	draw(-1);
}
