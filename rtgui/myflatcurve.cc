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
#include <myflatcurve.h>
#include <curves.h>
#include <string.h>
#include <gdkmm/types.h>

MyFlatCurve::MyFlatCurve () {

    innerWidth = get_allocation().get_width() - RADIUS * 2;
    innerHeight = get_allocation().get_height() - RADIUS * 2;
    prevInnerHeight = innerHeight;
    lit_point = -1;
    closest_point = 0;
    buttonPressed = false;
    editedHandle = FCT_EditedHandle_None;
    area = FCT_Area_None;
    tanHandlesDisplayed = false;
    periodic = true;

    //bghist = new unsigned int[256];

    signal_event().connect( sigc::mem_fun(*this, &MyFlatCurve::handleEvents) );

    // By default, we create a curve with 8 control points
    curve.type = FCT_MinMaxCPoints;

    defaultCurve();
}

/*MyFlatCurve::~MyFlatCurve () {
}*/

std::vector<double> MyFlatCurve::get_vector (int veclen) {

	// Create the output variable
    std::vector<double> convertedValues;

    // Get the curve control points
    std::vector<double> curveDescr = getPoints ();
    rtengine::FlatCurve* rtcurve = new rtengine::FlatCurve (curveDescr, periodic, veclen*1.5 > 5000 ? 5000 : veclen*1.5);

    // Create the sample values that will be converted
    std::vector<double> samples;
    samples.resize (veclen);
    for (int i = 0; i < veclen; i++)
    	samples[i] = (double) i / (veclen - 1.0);

    // Converting the values
    rtcurve->getVal (samples, convertedValues);

    // Cleanup and return
    delete rtcurve;
    return convertedValues;
}

void MyFlatCurve::interpolate () {

    prevInnerHeight = innerHeight;
    point.resize (innerWidth+1);
    std::vector<double> vector = get_vector (innerWidth+1);
    for (int i = 0; i <= innerWidth; ++i)
        point[i] = Gdk::Point ((double)RADIUS+0.5 + i, (double)RADIUS+0.5 + (double)innerHeight*(1.-vector[i]));
    upoint.clear ();
    lpoint.clear ();

    /*if (curve.type==FCT_Parametric && activeParam>0) {
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
    }*/
}

void MyFlatCurve::draw () {
    if (!pixmap)
        return;

    // re-calculate curve if dimensions changed
    if (prevInnerHeight != innerHeight || (int)point.size() != (innerWidth+1)) {
        interpolate ();

    }

    Gtk::StateType state = !is_sensitive() ? Gtk::STATE_INSENSITIVE : Gtk::STATE_NORMAL;

    Glib::RefPtr<Gtk::Style> style = get_style ();
    Cairo::RefPtr<Cairo::Context> cr = pixmap->create_cairo_context();

    // bounding rectangle
    Gdk::Color c = style->get_bg (Gtk::STATE_NORMAL);
    cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
    cr->rectangle (0, 0, innerWidth+RADIUS*2+1.5, innerHeight+RADIUS*2+1.5);
    cr->fill ();

    // histogram in the background
    /*if (bghistvalid) {
        // find highest bin
        unsigned int histheight = 0;
        for (int i=0; i<256; i++)
            if (bghist[i]>histheight)
	            histheight = bghist[i];
        // draw histogram
        cr->set_line_width (1.0);
        double stepSize = (innerWidth-1) / 256.0;
        cr->move_to (RADIUS, innerHeight-1+RADIUS);
        c = style->get_fg (Gtk::STATE_INSENSITIVE);
        cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
        for (int i=0; i<256; i++) {
            double val = bghist[i] * (double)(innerHeight-2) / (double)histheight;
            if (val>innerHeight-1)
                val = innerHeight-1;
      	    if (i>0)
       	        cr->line_to (i*stepSize+RADIUS, innerHeight-1+RADIUS-val);
    	}
        cr->line_to (innerWidth-1+RADIUS, innerHeight-1+RADIUS);
    	cr->fill ();
    }*/

    cr->set_line_cap(Cairo::LINE_CAP_SQUARE);

    // draw the grid lines:
    cr->set_line_width (1.0);
    cr->set_antialias (Cairo::ANTIALIAS_NONE);
    c = style->get_dark (state);
    cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());

    double x0 = (double)RADIUS-0.5;
    double x1 = (double)RADIUS-0.5 + (double)innerWidth + 2.;
    double y0 = (double)RADIUS-0.5;
    double y1 = (double)RADIUS-0.5 + (double)innerHeight + 2.;
    for (int i = 0; i < 5; i++) {
        cr->move_to (x0, y0);
        cr->line_to (x0, y1);
        cr->line_to (x1, y1);
        cr->line_to (x1, y0);
        cr->line_to (x0, y0);
    }
    /*for (int i = 0; i < 5; i++) {
        double currX = (double)RADIUS-0.5 + (double)i*((double)innerWidth + 2.)/4.;
        double currY = (double)RADIUS-0.5 + (double)i*((double)innerHeight + 2.)/4.;
        cr->move_to (x0, currY);
        cr->line_to (x1, currY);
        cr->move_to (currX, y0);
        cr->line_to (currX, y1);
    }*/
    cr->stroke ();

    // draw f(x)=0.5 line
    c = style->get_fg (state);
    cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
    std::valarray<double> ds (1);
    ds[0] = 4;
    cr->set_dash (ds, 0);
    cr->move_to (x0, (double)RADIUS+0.5 + (double)innerHeight/2.);
    cr->line_to (x1, (double)RADIUS+0.5 + (double)innerHeight/2.);
    cr->stroke ();

    cr->unset_dash ();

    cr->set_antialias (Cairo::ANTIALIAS_SUBPIXEL);

    cr->set_line_width (1.0);

    // draw the color feedback of the control points
    if (colorProvider) {

        //if (curve.type!=FCT_Parametric)
            for (int i=0; i<(int)curve.x.size(); ++i) {

                if (curve.x[i] != -1.) {

                    cr->set_line_width (1.0);
                    colorProvider->colorForValue(curve.x[i], 0.5);
                    cr->set_source_rgb (colorProvider->red, colorProvider->green, colorProvider->blue);

                    double x = (double)RADIUS+0.5 + innerWidth*curve.x[i];

                    if (i == lit_point && editedHandle&(FCT_EditedHandle_CPointUD|FCT_EditedHandle_CPoint|FCT_EditedHandle_CPointX)) {
                        cr->set_line_width (4.0);
                    }
                    cr->move_to (x, (double)RADIUS+0.5);
                    cr->line_to (x, (double)RADIUS+0.5 + innerHeight);
                    cr->stroke ();
                    cr->set_line_width (1.0);

                    // draw the lit_point's horizontal line
                    if (i == lit_point) {

                        if (area&(FCT_Area_H|FCT_Area_V|FCT_Area_Point) || editedHandle==FCT_EditedHandle_CPointUD) {

                            if (editedHandle&(FCT_EditedHandle_CPointUD|FCT_EditedHandle_CPoint|FCT_EditedHandle_CPointY)) {
                                cr->set_line_width (4.0);
                            }

                            colorProvider->colorForValue(curve.x[i], curve.y[i]);
                            cr->set_source_rgb (colorProvider->red, colorProvider->green, colorProvider->blue);

                            double y = (double)RADIUS+0.5 + (double)innerHeight*(1.-curve.y[lit_point]);
                            cr->move_to (                 RADIUS, y);
                            cr->line_to ((double)RADIUS+0.5 + (double)innerWidth, y);
                            cr->stroke ();
                        }
                    }
                }
            }
        // endif
	    cr->set_line_width (1.0);
    }
    else {
        cr->set_source_rgb (0.5, 0.0, 0.0);

        if (area==(FCT_Area_H|FCT_Area_V|FCT_Area_Point) || editedHandle==FCT_EditedHandle_CPointUD) {
        	double position;

            // draw the lit_point's vertical line
            if (editedHandle==(FCT_EditedHandle_CPointUD|FCT_EditedHandle_CPoint|FCT_EditedHandle_CPointY)) {
                cr->set_line_width (2.0);
            }
            position = (double)RADIUS+0.5 + (double)innerWidth*curve.x[lit_point];
            cr->move_to (position, (double)RADIUS+0.5);
            cr->line_to (position, (double)RADIUS+0.5 + (double)innerHeight);
            cr->stroke ();
            cr->set_line_width (1.0);

            // draw the lit_point's horizontal line
            if (editedHandle==(FCT_EditedHandle_CPointUD|FCT_EditedHandle_CPoint|FCT_EditedHandle_CPointY)) {
                cr->set_line_width (2.0);
            }
            position = (double)RADIUS+0.5 + (double)innerHeight*(1.-curve.y[lit_point]);
            cr->move_to ((double)RADIUS+0.5                     , position);
            cr->line_to ((double)RADIUS+0.5 + (double)innerWidth, position);
            cr->stroke ();
            cr->set_line_width (1.0);
        }
    }

    double lineMinLength = 1. / innerWidth * SQUARE * 0.9;
    if (lit_point!=-1 && getHandles(lit_point) && curve.x[lit_point]!=-1.) {
        double x = (double)RADIUS+0.5 + (double)innerWidth * curve.x[lit_point];
        double y = (double)RADIUS+0.5 + (double)innerHeight * (1.-curve.y[lit_point]);
        double x2;
        double square;
        bool crossingTheFrame;

        // left handle is yellow
        // TODO: finding a way to set the left handle color for flat curve editor
        cr->set_source_rgb (1.0, 1.0, 0.0);

        // draw tangential vectors

        crossingTheFrame = false;
        // We display the line only if it's longer than the handle knot half-size...
        if (leftTanX < -0.00001) {
            leftTanX += 1.0;
            crossingTheFrame = true;
        }
        x2 = (double)RADIUS+0.5 + (double)innerWidth * leftTanX;
        if (curve.x[lit_point] - leftTanX > lineMinLength || crossingTheFrame) {
            // The left tangential vector reappear on the right side
            // draw the line
            cr->move_to (x, y);
            if (crossingTheFrame) {
                cr->line_to ((double)RADIUS+0.5, y);
                cr->stroke ();
                cr->move_to ((double)RADIUS+0.5 + (double)innerWidth, y);
            }
            cr->line_to (x2, y);
            cr->stroke ();
        }
        // draw tangential knot
        square = area == FCT_Area_LeftTan ? SQUARE*2. : SQUARE;
        cr->move_to(x2-square, y+square);
        cr->line_to(x2+square, y+square);
        cr->line_to(x2+square, y-square);
        cr->line_to(x2-square, y-square);
        cr->line_to(x2-square, y+square);
        cr->fill();

        // right handle is blue
        // TODO: finding a way to set the right handle color for flat curve editor
        cr->set_source_rgb (0.0, 0.0, 1.0);

        // draw tangential vectors

        crossingTheFrame = false;
        // We display the line only if it's longer than the handle knot half-size...
        if (rightTanX > 1.00001) {
            rightTanX -= 1.0;
            crossingTheFrame = true;
        }
        x2 = (double)RADIUS+0.5 + (double)innerWidth * rightTanX;
        if (rightTanX - curve.x[lit_point] > lineMinLength || crossingTheFrame) {
            // The left tangential vector reappear on the right side
            // draw the line
            cr->move_to (x, y);
            if (crossingTheFrame) {
                cr->line_to ((double)RADIUS+0.5 + (double)innerWidth, y);
                cr->stroke ();
                cr->move_to ((double)RADIUS+0.5, y);
            }
            cr->line_to (x2, y);
            cr->stroke ();
        }
        // draw tangential knot
        square = area == FCT_Area_RightTan ? SQUARE*2. : SQUARE;
        cr->move_to(x2-square, y+square);
        cr->line_to(x2+square, y+square);
        cr->line_to(x2+square, y-square);
        cr->line_to(x2-square, y-square);
        cr->line_to(x2-square, y+square);
        cr->fill();
    }

    // draw curve
    cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
    cr->move_to (point[0].get_x(), point[0].get_y());
    for (int i=1; i<(int)point.size(); i++)
        cr->line_to (point[i].get_x(), point[i].get_y());
    cr->stroke ();

    // draw bullets
    //if (curve.type!=FCT_Parametric)
        for (int i = 0; i < (int)curve.x.size(); ++i) {
            if (curve.x[i] != -1.) {
                if (i == lit_point) {
                	if (colorProvider) {
                        colorProvider->colorForValue(curve.x[i], curve.y[i]);
                        cr->set_source_rgb (colorProvider->red, colorProvider->green, colorProvider->blue);
                	}
                	else
                		cr->set_source_rgb (1.0, 0.0, 0.0);
                }
                else if (i == snapToElmt) {
                    cr->set_source_rgb (1.0, 0.0, 0.0);
                }
                else if (curve.y[i] == 0.5)
                    cr->set_source_rgb (0.0, 0.5, 0.0);
                else
                    cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
                double x = (double)RADIUS+0.5 + (double)innerWidth * curve.x[i];    // project (curve.x[i], 0, 1, innerWidth);
                double y = (double)RADIUS+0.5 + (double)innerHeight * (1.-curve.y[i]); // project (curve.y[i], 0, 1, innerHeight);

                cr->arc (x, y, (double)RADIUS, 0, 2*M_PI);
                cr->fill ();
            }
        }
    // endif

    // draw the left and right tangent handles
    if (tanHandlesDisplayed) {
    	double top, bottom, left, right;
    	double halfSquareSizeX, halfSquareSizeY;

    	// LEFT handle
    	halfSquareSizeX = minDistanceX/2.;
    	halfSquareSizeY = minDistanceY/2.;
    	//halfSquareSizeX = area == FCT_Area_LeftTan ? minDistanceX : minDistanceX/2.;
    	//halfSquareSizeY = area == FCT_Area_LeftTan ? minDistanceY : minDistanceY/2.;
		top    = leftTanHandle.centerY + halfSquareSizeY;
		bottom = leftTanHandle.centerY - halfSquareSizeY;
		left   = leftTanHandle.centerX - halfSquareSizeX;
		right  = leftTanHandle.centerX + halfSquareSizeX;

		// yellow
        cr->set_source_rgb (1.0, 1.0, 0.0);
        cr->move_to((double)RADIUS+0.5 + (double)innerWidth * left,  (double)RADIUS+0.5 + (double)innerHeight * (1.-top));
        cr->line_to((double)RADIUS+0.5 + (double)innerWidth * right, (double)RADIUS+0.5 + (double)innerHeight * (1.-top));
        cr->line_to((double)RADIUS+0.5 + (double)innerWidth * right, (double)RADIUS+0.5 + (double)innerHeight * (1.-bottom));
        cr->line_to((double)RADIUS+0.5 + (double)innerWidth * left,  (double)RADIUS+0.5 + (double)innerHeight * (1.-bottom));
        cr->line_to((double)RADIUS+0.5 + (double)innerWidth * left,  (double)RADIUS+0.5 + (double)innerHeight * (1.-top));
        cr->fill();

        // RIGHT handle
    	//halfSquareSizeX = area == FCT_Area_RightTan ? minDistanceX : minDistanceX/2.;
    	//halfSquareSizeY = area == FCT_Area_RightTan ? minDistanceY : minDistanceY/2.;
		top    = rightTanHandle.centerY + halfSquareSizeY;
		bottom = rightTanHandle.centerY - halfSquareSizeY;
		left   = rightTanHandle.centerX - halfSquareSizeX;
		right  = rightTanHandle.centerX + halfSquareSizeX;

        // blue
        cr->set_source_rgb (0.0, 0.0, 1.0);
        cr->move_to((double)RADIUS+0.5 + (double)innerWidth * left,  (double)RADIUS+0.5 + (double)innerHeight * (1.-top));
        cr->line_to((double)RADIUS+0.5 + (double)innerWidth * right, (double)RADIUS+0.5 + (double)innerHeight * (1.-top));
        cr->line_to((double)RADIUS+0.5 + (double)innerWidth * right, (double)RADIUS+0.5 + (double)innerHeight * (1.-bottom));
        cr->line_to((double)RADIUS+0.5 + (double)innerWidth * left,  (double)RADIUS+0.5 + (double)innerHeight * (1.-bottom));
        cr->line_to((double)RADIUS+0.5 + (double)innerWidth * left,  (double)RADIUS+0.5 + (double)innerHeight * (1.-top));
        cr->fill();
    }

    get_window()->draw_drawable (style->get_fg_gc (state), pixmap, 0, 0, 0, 0, innerWidth + RADIUS * 2 + 1, innerHeight + RADIUS * 2 + 1);
}

/*
 * Return the X1, X2, Y position of the tangential handles.
 */
bool MyFlatCurve::getHandles(int n) {
	int N = curve.x.size();
	double prevX, nextX;
	double prevY, nextY;
	double prevTan, nextTan;
	double x, y, leftTan, rightTan;

	if (n == -1) return false;

	x = curve.x[n];
	y = curve.y[n];
	leftTan = curve.leftTangent[n];
	rightTan = curve.rightTangent[n];

	if (!n) {
		// first point, the left handle is then computed with the last point's right handle
		prevX = curve.x[N-1]-1.0;
		prevY = curve.y[N-1];
		prevTan = curve.rightTangent[N-1];

		nextX = curve.x[n+1];
		nextY = curve.y[n+1];
		nextTan = curve.leftTangent[n+1];
	}
	else if (n == N-1) {
		// last point, the right handle is then computed with the first point's left handle
		prevX = curve.x[n-1];
		prevY = curve.y[n-1];
		prevTan = curve.rightTangent[n-1];

		nextX = curve.x[0]+1.0;
		nextY = curve.y[0];
		nextTan = curve.leftTangent[0];
	}
	else {
		// last point, the right handle is then computed with the first point's left handle
		prevX = curve.x[n-1];
		prevY = curve.y[n-1];
		prevTan = curve.rightTangent[n-1];

		nextX = curve.x[n+1];
		nextY = curve.y[n+1];
		nextTan = curve.leftTangent[n+1];
	}

	if (leftTan == 0.0) leftTanX = x;
	else if (leftTan == 1.0) leftTanX = prevX;
	else {
		leftTanX = (prevX - x) * leftTan + x;
	}

	if (rightTan == 0.0) rightTanX = x;
	else if (rightTan == 1.0) rightTanX = nextX;
	else {
		rightTanX = (nextX - x) * rightTan + x;
	}
	return true;
}

bool MyFlatCurve::handleEvents (GdkEvent* event) {

	CursorShape new_type = cursor_type;
	int src, dst;
	std::vector<double>::iterator itx, ity, itlt, itrt;

	snapToElmt = -100;
	bool retval = false;
	int num = (int)curve.x.size();

	/* innerWidth and innerHeight are the size of the graph */
	innerWidth = get_allocation().get_width() - 2*RADIUS - 1;
	innerHeight = get_allocation().get_height() - 2*RADIUS - 1;

	minDistanceX = (double)(MIN_DISTANCE) / (double)innerWidth;
	minDistanceY = (double)(MIN_DISTANCE) / (double)innerHeight;

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
		draw ();
		retval = true;
		break;

	case Gdk::BUTTON_PRESS:
		//if (curve.type!=FCT_Parametric) {
			if (event->button.button == 1) {
				buttonPressed = true;
				add_modal_grab ();

				// get the pointer position
				getCursorPosition(event);
				getMouseOverArea();

				// hide the tangent handles
				tanHandlesDisplayed = false;

				// Action on BUTTON_PRESS and no edited point
				switch (area) {

				case (FCT_Area_Insertion):
					new_type = CSMove;

					/* insert a new control point */
					if (num > 0) {
						if (clampedX > curve.x[closest_point])
							++closest_point;
					}
					itx = curve.x.begin();
					ity = curve.y.begin();
					itlt = curve.leftTangent.begin();
					itrt = curve.rightTangent.begin();
					for (int i=0; i<closest_point; i++) { itx++; ity++; itlt++; itrt++; }
					curve.x.insert (itx, 0);
					curve.y.insert (ity, 0);
					curve.leftTangent.insert (itlt, 0);
					curve.rightTangent.insert (itrt, 0);
					num++;

					// the graph is refreshed only if a new point is created
					curve.x[closest_point] = clampedX;
					curve.y[closest_point] = clampedY;
					curve.leftTangent[closest_point] = 0.35;
					curve.rightTangent[closest_point] = 0.35;

					interpolate ();
					draw ();
					notifyListener ();

					lit_point = closest_point;

					// point automatically activated
					editedHandle = FCT_EditedHandle_CPoint;
					ugpX = curve.x[lit_point];
					ugpY = curve.y[lit_point];
					break;

				case (FCT_Area_Point):
					new_type = CSMove;
					editedHandle = FCT_EditedHandle_CPoint;
					ugpX = curve.x[lit_point];
					ugpY = curve.y[lit_point];
					break;

				case (FCT_Area_H):
				case (FCT_Area_V):
					new_type = CSMove;
					editedHandle = FCT_EditedHandle_CPointUD;
					ugpX = curve.x[lit_point];
					ugpY = curve.y[lit_point];
					break;

				case (FCT_Area_LeftTan):
					new_type = CSEmpty;
					editedHandle = FCT_EditedHandle_LeftTan;
					ugpX = curve.leftTangent[lit_point];
					break;

				case (FCT_Area_RightTan):
					new_type = CSEmpty;
					editedHandle = FCT_EditedHandle_RightTan;
					ugpX = curve.rightTangent[lit_point];
					break;

				default:
					break;
				}

			}
			if (buttonPressed) retval = true;
		//}
		break;

	case Gdk::BUTTON_RELEASE:
		//if (curve.type!=FCT_Parametric) {
			if (buttonPressed && event->button.button == 1) {
				buttonPressed = false;
				enum MouseOverAreas prevArea = area;
				remove_modal_grab ();

				int previous_lit_point = lit_point;

				// Removing any deleted point if we were previously modifying the point position
				if (editedHandle & (FCT_EditedHandle_CPoint|FCT_EditedHandle_CPointX|FCT_EditedHandle_CPointY)) {
					/* delete inactive points: */
					itx = curve.x.begin();
					ity = curve.y.begin();
					itlt = curve.leftTangent.begin();
					itrt = curve.rightTangent.begin();
					for (src = dst = 0; src < num; ++src)
						if (curve.x[src] >= 0.0) {
							curve.x[dst] = curve.x[src];
							curve.y[dst] = curve.y[src];
							curve.leftTangent[dst] = curve.leftTangent[src];
							curve.rightTangent[dst] = curve.rightTangent[src];
							++dst;
							++itx;
							++ity;
							++itlt;
							++itrt;
						}
					if (dst < src) {

						// curve cleanup
						curve.x.erase (itx, curve.x.end());
						curve.y.erase (ity, curve.y.end());
						curve.leftTangent.erase (itlt, curve.leftTangent.end());
						curve.rightTangent.erase (itrt, curve.rightTangent.end());
						if (!curve.x.size()) {
							curve.x.push_back (0.5);
							curve.y.push_back (0.5);
							curve.leftTangent.push_back (0.3);
							curve.rightTangent.push_back (0.3);
							interpolate ();
						}
					}
				}

				editedHandle = FCT_EditedHandle_None;
				lit_point = -1;

				// get the pointer position
				getCursorPosition(event);
				getMouseOverArea();

				switch (area) {

				case (FCT_Area_Insertion):
					new_type = CSArrow;
					break;

				case (FCT_Area_Point):
					new_type = CSMove;
					break;

				case (FCT_Area_H):
					new_type = CSResizeHeight;
					break;

				case (FCT_Area_V):
					new_type = CSMove;
					break;

				case (FCT_Area_LeftTan):
					new_type = CSEmpty;
					break;

				case (FCT_Area_RightTan):
					new_type = CSEmpty;
					break;

				default:
					break;
				}

				if ((lit_point != previous_lit_point) || (prevArea != area))
					draw ();
				retval = true;
				//notifyListener ();
			}
		//}
		break;

	case Gdk::MOTION_NOTIFY:
		if (curve.type == FCT_Linear || curve.type == FCT_MinMaxCPoints) {

			int previous_lit_point = lit_point;
			enum MouseOverAreas prevArea = area;

			snapToMinDist = 10.;
			snapToVal = 0.;
			snapToElmt = -100;

			// get the pointer position
			getCursorPosition(event);
			getMouseOverArea();

			if (editedHandle == FCT_EditedHandle_CPointUD) {
				double dX = deltaX;
				double dY = deltaY;
				if (dX < 0.) dX = -dX;
				if (dY < 0.) dY = -dY;
				if (dX > dY) {
					editedHandle = FCT_EditedHandle_CPointX;
					area = FCT_Area_V;
					new_type = CSResizeWidth;
				}
				else {
					editedHandle = FCT_EditedHandle_CPointY;
					area = FCT_Area_H;
					new_type = CSResizeHeight;
				}
			}

			switch (editedHandle) {

			case (FCT_EditedHandle_None): {

				if ((lit_point != -1 && previous_lit_point != lit_point) && area & (FCT_Area_V|FCT_Area_Point)) {

					bool sameSide = false;

					// display the handles
					tanHandlesDisplayed = true;

					if (curve.x[lit_point] < 3.*minDistanceX) {
						// lit_point too near of the left border -> both handles are displayed on the right of the vertical line
						rightTanHandle.centerX  = leftTanHandle.centerX  = curve.x[lit_point] + 2.*minDistanceX;
						sameSide = true;
					}
					else if (curve.x[lit_point] > 1.-3.*minDistanceX) {
						// lit_point too near of the left border -> both handles are displayed on the right of the vertical line
						rightTanHandle.centerX = leftTanHandle.centerX = curve.x[lit_point]  - 2.*minDistanceX;
						sameSide = true;
					}
					else {
						leftTanHandle.centerX = curve.x[lit_point] - 2.*minDistanceX;
						rightTanHandle.centerX = curve.x[lit_point] + 2.*minDistanceX;
					}
					if (sameSide) {
						if (clampedY > 1.-2.*minDistanceY) {
							// lit_point too near of the top border
							leftTanHandle.centerY = 1. - minDistanceY;
							rightTanHandle.centerY = 1. - 3.*minDistanceY;
						}
						else if (clampedY < 2.*minDistanceY) {
							// lit_point too near of the bottom border
							leftTanHandle.centerY = 3.*minDistanceY;
							rightTanHandle.centerY = minDistanceY;
						}
						else {
							leftTanHandle.centerY = clampedY + minDistanceY;
							rightTanHandle.centerY = clampedY - minDistanceY;
						}
					}
					else {
						if (clampedY > 1.-minDistanceY) {
							rightTanHandle.centerY = leftTanHandle.centerY = 1. - minDistanceY;
						}
						else if (clampedY < minDistanceY) {
							rightTanHandle.centerY = leftTanHandle.centerY = minDistanceY;
						}
						else {
							rightTanHandle.centerY = leftTanHandle.centerY = clampedY;
						}
					}
				}
				else if (lit_point == -1) {
					tanHandlesDisplayed = false;
				}


				switch (area) {

				case (FCT_Area_Insertion):
					new_type = CSPlus;
					break;
				case (FCT_Area_Point):
					//new_type = CSMove;
					//break;
				case (FCT_Area_V):
					new_type = CSMove;
					break;
				case (FCT_Area_H):
					new_type = CSResizeHeight;
					break;
				case (FCT_Area_LeftTan):
					new_type = CSMoveLeft;
					break;
				case (FCT_Area_RightTan):
					new_type = CSMoveRight;
					break;
				case (FCT_Area_None):
				default:
					new_type = CSArrow;
					break;
				}

				if ((lit_point != previous_lit_point) || (prevArea != area))
					draw ();
				break;
			}

			case (FCT_EditedHandle_CPoint):
				movePoint(true, true);
				break;

			case (FCT_EditedHandle_CPointX):
				movePoint(true, false);
				break;

			case (FCT_EditedHandle_CPointY):
				movePoint(false, true);
				break;

			case (FCT_EditedHandle_LeftTan): {
				double prevValue = curve.leftTangent[lit_point];

				ugpX -= deltaX*3;
				ugpX = CLAMP(ugpX, 0., 1.);
				if (snapTo) {
					snapCoordinate(0.0,  ugpX);
					snapCoordinate(0.35, ugpX);
					snapCoordinate(0.5,  ugpX);
					snapCoordinate(1.0,  ugpX);
					curve.leftTangent[lit_point] = snapToVal;
				}
				else {
					curve.leftTangent[lit_point] = ugpX;
				}

				if (curve.leftTangent[lit_point] != prevValue) {
					interpolate ();
					draw ();
					notifyListener ();
				}
				break;
			}

			case (FCT_EditedHandle_RightTan): {
				double prevValue = curve.rightTangent[lit_point];

				ugpX += deltaX*3;
				ugpX = CLAMP(ugpX, 0., 1.);
				if (snapTo) {
					snapCoordinate(0.0,  ugpX);
					snapCoordinate(0.35, ugpX);
					snapCoordinate(0.5,  ugpX);
					snapCoordinate(1.0,  ugpX);
					curve.rightTangent[lit_point] = snapToVal;
				}
				else {
					curve.rightTangent[lit_point] = ugpX;
				}

				if (curve.rightTangent[lit_point] != prevValue) {
					interpolate ();
					draw ();
					notifyListener ();
				}
				break;
			}

			// already processed before the "switch" instruction
			//case (FCT_EditedHandle_CPointUD):

			default:
				break;
			}
		}

		retval = true;
		break;

	case Gdk::LEAVE_NOTIFY:
		// Pointer can LEAVE even when dragging the point, so we don't modify the cursor in this case
		// The cursor will have to LEAVE another time after the drag...
		if (editedHandle == FCT_EditedHandle_None) {
			new_type = CSArrow;
			lit_point = -1;
			tanHandlesDisplayed = false;
			draw ();
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

void MyFlatCurve::movePoint(bool moveX, bool moveY) {

	// bounds of the grabbed point
	double leftBound;
	double rightBound;
	double const bottomBound = 0.;
	double const topBound    = 1.;

	// we memorize the previous position of the point, for optimization purpose
	double prevPosX = curve.x[lit_point];
	double prevPosY = curve.y[lit_point];

	int nbPoints = (int)curve.x.size();

	// left and right bound rely on curve periodicity
	leftBound         = (lit_point == 0         ) ? (periodic ? curve.x[nbPoints-1]-1. : 0.) : curve.x[lit_point-1];
	rightBound        = (lit_point == nbPoints-1) ? (periodic ? curve.x[0]+1.          : 1.) : curve.x[lit_point+1];

	double leftDeletionBound   = leftBound   - minDistanceX;
	double rightDeletionBound  = rightBound  + minDistanceX;
	double bottomDeletionBound = bottomBound - minDistanceY;
	double topDeletionBound    = topBound    + minDistanceY;

	if (moveX) {
		// we memorize the previous position of the point, for optimization purpose
		ugpX += deltaX;

		// handling periodicity (the first and last point can reappear at the other side of the X range)
		if (periodic) {
			if (lit_point==0 && ugpX<0.) {
				// the first point has to be placed at the tail of the point list
				std::vector<double>::iterator itx, ity, itlt, itrt;

				ugpX += 1.;
				leftBound += 1.;
				rightBound += 1.;
				leftDeletionBound += 1.;
				rightDeletionBound += 1.;

				// adding a copy of the first point to the tail of the list
				curve.x.push_back(curve.x[0]);
				curve.y.push_back(curve.y[0]);
				curve.leftTangent.push_back(curve.leftTangent[0]);
				curve.rightTangent.push_back(curve.rightTangent[0]);

				// deleting the first point
				itx = curve.x.begin();
				ity = curve.y.begin();
				itlt = curve.leftTangent.begin();
				itrt = curve.rightTangent.begin();

				curve.x.erase(itx);
				curve.y.erase(ity);
				curve.leftTangent.erase(itlt);
				curve.rightTangent.erase(itrt);

				lit_point = nbPoints-1;
			}
			else if (lit_point==(nbPoints-1) && ugpX>1.) {
				// the last point has to be placed at the head of the point list
				std::vector<double>::iterator itx, ity, itlt, itrt;

				ugpX -= 1.;
				leftBound -= 1.;
				rightBound -= 1.;
				leftDeletionBound -= 1.;
				rightDeletionBound -= 1.;

				// adding a copy of the last point to the head of the list
				itx = curve.x.begin();
				ity = curve.y.begin();
				itlt = curve.leftTangent.begin();
				itrt = curve.rightTangent.begin();

				curve.x.insert(itx,0);
				curve.y.insert(ity,0);
				curve.leftTangent.insert(itlt,0);
				curve.rightTangent.insert(itrt,0);

				curve.x[0] = curve.x[nbPoints];
				curve.y[0] = curve.y[nbPoints];
				curve.leftTangent[0] = curve.leftTangent[nbPoints];
				curve.rightTangent[0] = curve.rightTangent[nbPoints];

				// deleting the last point
				curve.x.pop_back();
				curve.y.pop_back();
				curve.leftTangent.pop_back();
				curve.rightTangent.pop_back();

				lit_point = 0;
			}
		}

		// handling limitations along X axis
		if (ugpX >= rightDeletionBound && nbPoints>2) {
			curve.x[lit_point] = -1.;
		}
		else if (ugpX <= leftDeletionBound && nbPoints>2) {
			curve.x[lit_point] = -1.;
		}
		else
			// nextPosX is in bounds
			curve.x[lit_point] = CLAMP(ugpX, leftBound, rightBound);
	}

	if (moveY) {

		// we memorize the previous position of the point, for optimization purpose
		ugpY += deltaY;

		// snapping point to specific values
		if (snapTo && curve.x[lit_point] != -1) {

			// the unclamped grabbed point is brought back in the range
			ugpY = CLAMP(ugpY, 0.0, 1.0);

			if (lit_point == 0) {
				int prevP = curve.y.size()-1;
				if (snapCoordinate(curve.y[prevP], ugpY)) snapToElmt = prevP;
			}
			else {
				int prevP = lit_point-1;
				if (snapCoordinate(curve.y[prevP], ugpY)) snapToElmt = prevP;
			}

			if (curve.y.size() > 2) {
				if (lit_point == (curve.y.size()-1)) {
					if (snapCoordinate(curve.y[0], ugpY)) snapToElmt = 0;
				}
				else {
					int nextP = lit_point+1;
					if (snapCoordinate(curve.y[nextP], ugpY)) snapToElmt = nextP;
				}
			}
			if (snapCoordinate(1.0, ugpY)) snapToElmt = -3;
			if (snapCoordinate(0.5, ugpY)) snapToElmt = -2;
			if (snapCoordinate(0.0, ugpY)) snapToElmt = -1;

			curve.y[lit_point] = snapToVal;
		}

		// Handling limitations along Y axis
		if (ugpY >= topDeletionBound && nbPoints>2) {
			if (curve.x[lit_point] != -1.) {
				deletedPointX = curve.x[lit_point];
				curve.x[lit_point] = -1.;
				curve.y[lit_point] = ugpY; // This is only to force the redraw of the curve
			}
		}
		else if (ugpY <= bottomDeletionBound  && nbPoints>2) {
			if (curve.x[lit_point] != -1.) {
				deletedPointX = curve.x[lit_point];
				curve.x[lit_point] = -1.;
				curve.y[lit_point] = ugpY; // This is only to force the redraw of the curve
			}
		}
		else {
			// nextPosY is in the bounds
			if (!snapTo)  curve.y[lit_point] = CLAMP(ugpY, 0.0, 1.0);
			if (!moveX && curve.x[lit_point] == -1.) {
				// bring back the X value of the point if it reappear
				curve.x[lit_point] = deletedPointX;
			}
		}
	}

	if (curve.x[lit_point] != prevPosX || curve.y[lit_point] != prevPosY) {
		// we recompute the curve only if we have to
		interpolate ();
		draw ();
		notifyListener ();
	}
}

// Set datas relative to cursor position
void MyFlatCurve::getCursorPosition(GdkEvent* event) {
	int tx, ty;
	int prevCursorX, prevCursorY;
	double incrementX = 1. / (double)innerWidth;
	double incrementY = 1. / (double)innerHeight;

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

	if (editedHandle != FCT_EditedHandle_None) {
		prevCursorX = cursorX;
		prevCursorY = cursorY;
    }
    cursorX = tx - RADIUS;
    cursorY = innerHeight - (ty - RADIUS);

    preciseCursorX = cursorX * incrementX;
    preciseCursorY = cursorY * incrementY;

    snapTo = false;

    // update deltaX/Y if the user drags a point
    if (editedHandle != FCT_EditedHandle_None) {
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
		clampedX = CLAMP (preciseCursorX, 0., 1.); // X position of the pointer from the origin of the graph
		clampedY = CLAMP (preciseCursorY, 0., 1.); // Y position of the pointer from the origin of the graph
    }

}

// Find out the active area under the cursor
void MyFlatCurve::getMouseOverArea () {

	// When dragging an element, editedHandle keep its value
	if (editedHandle == FCT_EditedHandle_None) { // && curve.type!=Parametric

		double minDist = 1000;	// used to find out the point pointed by the cursor (over it)
		double minDistX = 1000;	// used to find out the closest point
		double dX, dY;
		double absDX;
		double dist;
		bool aboveVLine = false;

		// NB: this function assume that the graph's shape is a square

		// Check if the cursor is over a tangent handle
		if (tanHandlesDisplayed) {

			if (preciseCursorX>=(leftTanHandle.centerX-minDistanceX) && preciseCursorX<=(leftTanHandle.centerX+minDistanceX+0.00001)
			&&  preciseCursorY>=(leftTanHandle.centerY-minDistanceY) && preciseCursorY<=(leftTanHandle.centerY+minDistanceY)) {
					area = FCT_Area_LeftTan;
					return;
			}
			if (preciseCursorX>=(rightTanHandle.centerX-minDistanceX-0.00001) && preciseCursorX<=(rightTanHandle.centerX+minDistanceX)
			&&  preciseCursorY>=(rightTanHandle.centerY-minDistanceY        ) && preciseCursorY<=(rightTanHandle.centerY+minDistanceY)) {
					area = FCT_Area_RightTan;
					return;
			}
		}

		area = FCT_Area_None;
		closest_point = 0;
		lit_point = -1;

		for (int i = 0; i < (int)curve.x.size(); i++) {
			if (curve.x[i] != -1) {
				dX = curve.x[i] - preciseCursorX;
				absDX = dX>0 ? dX : -dX;
				if (absDX < minDistX) {
					minDistX = absDX;
					closest_point = i;
					lit_point = i;
				}
				if (absDX <= minDistanceX) {
					aboveVLine = true;
					dY = curve.y[i] - preciseCursorY;
					dist = sqrt(dX*dX + dY*dY);
					if (dist < minDist) {
						minDist = dist;
					}
				}
			}
		}
		if (minDist <= minDistanceX) {
			// the cursor is over the point
			area = FCT_Area_Point;
		}
		else if (aboveVLine) {
			area = FCT_Area_V;
		}
		else {
			// Check if the cursor is in an insertion area
			lit_point = -1;
			if (preciseCursorX>=0.0 && preciseCursorX<=1.0 && preciseCursorY>=0.0 && preciseCursorY<=1.0)
				area = FCT_Area_Insertion;
		}
	}
}

std::vector<double> MyFlatCurve::getPoints () {
    std::vector<double> result;
    /*if (curve.type==FCT_Parametric) {
        result.push_back ((double)(Parametric));
        for (int i=0; i<(int)curve.x.size(); i++) {
            result.push_back (curve.x[i]);
        }
    }
    else {*/
    	// the first value gives the type of the curve
        if (curve.type==FCT_Linear)
            result.push_back ((double)(FCT_Linear));
        else if (curve.type==FCT_MinMaxCPoints)
            result.push_back ((double)(FCT_MinMaxCPoints));
        // then we push all the points coordinate
        for (int i=0; i<(int)curve.x.size(); i++) {
            if (curve.x[i]>=0) {
                result.push_back (curve.x[i]);
                result.push_back (curve.y[i]);
                result.push_back (curve.leftTangent[i]);
                result.push_back (curve.rightTangent[i]);
            }
        }
    //}
    return result;
}

void MyFlatCurve::setPoints (const std::vector<double>& p) {
    int ix = 0;
    FlatCurveType t = (FlatCurveType)p[ix++];
    curve.type = t;
    if (t==FCT_MinMaxCPoints) {
        curve.x.clear ();
        curve.y.clear ();
        curve.leftTangent.clear();
        curve.rightTangent.clear();
        for (int i=0; i<(int)p.size()/4; i++) {
            curve.x.push_back (p[ix++]);
            curve.y.push_back (p[ix++]);
            curve.leftTangent.push_back (p[ix++]);
            curve.rightTangent.push_back (p[ix++]);
        }
    }
    pixmap.clear ();
    queue_draw ();
}

void MyFlatCurve::setType (FlatCurveType t) {

    curve.type = t;
    pixmap.clear ();
}

/*int flatmchistupdate (void* data) {

    gdk_threads_enter ();

    MyFlatCurveIdleHelper* mcih = (MyFlatCurveIdleHelper*)data;

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
}*/

/*void MyFlatCurve::updateBackgroundHistogram (unsigned int* hist) {

    if (hist!=NULL) {
        memcpy (bghist, hist, 256*sizeof(unsigned int));
        bghistvalid = true;
    }
    else
        bghistvalid = false;

    mcih->pending++;
    g_idle_add (flatmchistupdate, mcih);

}*/

void MyFlatCurve::reset() {
	innerWidth = get_allocation().get_width() - RADIUS * 2;
	innerHeight = get_allocation().get_height() - RADIUS * 2;

	switch (curve.type) {
	case  FCT_MinMaxCPoints :
		defaultCurve();
	    lit_point = -1;
        interpolate ();
		break;
	//case Parametric :
		// Nothing to do (?)
	default:
		break;
	}
	draw();
}

void MyFlatCurve::defaultCurve () {

	curve.x.clear();
	curve.y.clear();
    curve.leftTangent.clear();
    curve.rightTangent.clear();

    // Point for RGBCMY colors
    for (int i=0; i<6; i++) {
		curve.x.push_back((1./6.)*i);
		curve.y.push_back(0.5);
		curve.leftTangent.push_back(0.35);
		curve.rightTangent.push_back(0.35);
    }

}
