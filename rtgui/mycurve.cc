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
#include <mycurve.h>
#include <curves.h>
#include <string.h>
#include <gdkmm/types.h>

MyCurve::MyCurve () : listener(NULL), activeParam(-1), bghistvalid(false) {

    cursor_type = CSArrow;
    curve.type = Spline;
    height = 0;
    grab_point = -1;
    lit_point = -1;

    add_events(Gdk::EXPOSURE_MASK |	Gdk::POINTER_MOTION_MASK |	Gdk::POINTER_MOTION_HINT_MASK |	Gdk::ENTER_NOTIFY_MASK | Gdk::LEAVE_NOTIFY_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK | Gdk::BUTTON1_MOTION_MASK);
    signal_event().connect( sigc::mem_fun(*this, &MyCurve::handleEvents) );
    
    curve.x.push_back(0);
    curve.y.push_back(0);
    curve.x.push_back(1);
    curve.y.push_back(1);
    curve.type = Spline;

    mcih = new MyCurveIdleHelper;
    mcih->myCurve = this;
    mcih->destroyed = false;
    mcih->pending = 0;
}

MyCurve::~MyCurve () {
    
    if (mcih->pending)
        mcih->destroyed = true;
    else
        delete mcih;
}

std::vector<double> MyCurve::get_vector (int veclen) {

    std::vector<double> vector;
    vector.resize (veclen);

    if (curve.type != Parametric) {
        // count active points:
        double prev =- 1.0;
        int active = 0;
        int firstact = -1;
        for (int i = 0; i < curve.x.size(); ++i)
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
    rtengine::Curve* rtcurve = new rtengine::Curve (curveDescr, veclen*1.5);
    std::vector<double> t;
    t.resize (veclen);
    for (int i = 0; i < veclen; i++)
        t[i] = (double) i / (veclen - 1.0);
    rtcurve->getVal (t, vector);
    delete rtcurve;
    return vector;
}

void MyCurve::interpolate (int width, int height) {

    this->height = height;
    point.resize (width);
    std::vector<double> vector = get_vector (width);
    this->height = height;
    for (int i = 0; i < width; ++i)
        point[i] = Gdk::Point (RADIUS + i, RADIUS + height - (int)((height-1) * vector[i] + 0.5));
    upoint.clear ();
    lpoint.clear ();

    if (curve.type==Parametric && activeParam>0) {
        double tmp = curve.x[activeParam-1];
        if (activeParam>=4) {
            upoint.resize(width);
            lpoint.resize(width);
            curve.x[activeParam-1] = 100;
            vector = get_vector (width);
            for (int i = 0; i < width; ++i)
                upoint[i] = Gdk::Point (RADIUS + i, RADIUS + height - (int)((height-1) * vector[i] + 0.5));
            curve.x[activeParam-1] = -100;
            vector = get_vector (width);
            for (int i = 0; i < width; ++i)
                lpoint[i] = Gdk::Point (RADIUS + i, RADIUS + height - (int)((height-1) * vector[i] + 0.5));
            curve.x[activeParam-1] = tmp;
        }
    }
}

void MyCurve::draw (int width, int height, int handle) {
	// width and heigth are the size of the graph

    if (!pixmap)
        return;

    // re-calculate curve if dimensions changed
    if (this->height != height || point.size() != width)
        interpolate (width, height);

    Gtk::StateType state = Gtk::STATE_NORMAL;
    if (!is_sensitive())
        state = Gtk::STATE_INSENSITIVE;

    Glib::RefPtr<Gtk::Style> style = get_style ();
    Cairo::RefPtr<Cairo::Context> cr = pixmap->create_cairo_context();

    // bounding rectangle
    Gdk::Color c = style->get_bg (state);
    cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
    cr->rectangle (0, 0, width + RADIUS*2, height + RADIUS*2);
    cr->fill ();

    // histogram in the background
    if (bghistvalid) {
        // find heighest bin
        int histheight = 0;
        for (int i=0; i<256; i++)
            if (bghist[i]>histheight)
	            histheight = bghist[i];
        // draw histogram
        cr->set_line_width (1.0);
        double stepSize = (width-1) / 256.0;
        cr->move_to (RADIUS, height-1+RADIUS);
        cr->set_source_rgb (0.75, 0.75, 0.75);
        for (int i=0; i<256; i++) {
            double val = bghist[i] * (double)(height-2) / histheight;
            if (val>height-1)
                val = height-1;
      	    if (i>0)
       	        cr->line_to (i*stepSize+RADIUS, height-1+RADIUS-val);
    	}
        cr->line_to (width-1+RADIUS, height-1+RADIUS);
    	cr->fill ();
    }

    // draw the grid lines:
    cr->set_line_width (1.0);
    c = style->get_dark (state);
    cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
    cr->set_antialias (Cairo::ANTIALIAS_NONE);
    for (int i = 0; i < 5; i++) {
        cr->move_to (RADIUS, MAX(0,i * height / 4 - 1) + RADIUS);
        cr->line_to (width + RADIUS, MAX(0,i * height / 4 - 1) + RADIUS);
        cr->move_to (MAX(0,i * width / 4 - 1) + RADIUS, RADIUS);
        cr->line_to (MAX(0,i * width / 4 - 1) + RADIUS, height + RADIUS);
    }
    cr->stroke ();

    // draw f(x)=x line
    cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
    std::valarray<double> ds (1);
    ds[0] = 4;
    cr->set_dash (ds, 0);
    cr->move_to (RADIUS, height + RADIUS);
    cr->line_to (width + RADIUS, RADIUS);
    cr->stroke ();
    cr->unset_dash ();


    cr->set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
    cr->set_line_width (1.0);

    // draw upper and lower bounds
    if (curve.type==Parametric && activeParam>0 && lpoint.size()>1 && upoint.size()>1) {
        cr->set_source_rgba (0.0, 0.0, 0.0, 0.15);
        cr->move_to (upoint[0].get_x(), upoint[0].get_y());
        for (int i=1; i<upoint.size(); i++)
            cr->line_to (upoint[i].get_x(), upoint[i].get_y());
        cr->line_to (lpoint[lpoint.size()-1].get_x(), lpoint[lpoint.size()-1].get_y());
        for (int i=lpoint.size()-2; i>=0; i--)
            cr->line_to (lpoint[i].get_x(), lpoint[i].get_y());
        cr->line_to (upoint[0].get_x(), upoint[0].get_y());
        cr->fill ();
    }

    // draw the cage of the NURBS curve
    if (curve.type==NURBS) {
        std::valarray<double> ch_ds (1);
        ch_ds[0] = 2;
        cr->set_dash (ch_ds, 0);
        cr->set_source_rgb (0.0, 0.0, 0.0);
        std::vector<double> points = getPoints();
        for (int i = 1; i < points.size(); ) {
			double x = ((width-1) * points[i++] + 0.5)+RADIUS;    // project (curve.x[i], 0, 1, width);
			double y = height - ((height-1) * points[i++] + 0.5)+RADIUS; // project (curve.y[i], 0, 1, height);
			if (i==3)
				cr->move_to (x, y);
			else
				cr->line_to (x, y);
        }
        cr->stroke ();
        cr->unset_dash ();
    }

    // draw curve
    cr->set_source_rgb (0.0, 0.0, 0.0);
    cr->move_to (point[0].get_x(), point[0].get_y());
    for (int i=1; i<point.size(); i++)
        cr->line_to (point[i].get_x(), point[i].get_y());
    cr->stroke ();

    // draw bullets
    if (curve.type!=Parametric) 
        for (int i = 0; i < curve.x.size(); ++i) {
            cr->set_source_rgb ((i == handle ? 1.0 : 0.0), 0.0, 0.0);
            double x = ((width-1) * curve.x[i] + 0.5)+RADIUS;    // project (curve.x[i], 0, 1, width);
            double y = height - ((height-1) * curve.y[i] + 0.5)+RADIUS; // project (curve.y[i], 0, 1, height);

            cr->arc (x, y, RADIUS+0.5, 0, 2*M_PI);
            cr->fill ();
        }

    get_window()->draw_drawable (style->get_fg_gc (state), pixmap, 0, 0, 0, 0, width + RADIUS * 2, height + RADIUS * 2);
}

bool MyCurve::handleEvents (GdkEvent* event) {

    CursorShape new_type = cursor_type;
    int src, dst;
    GdkEventMotion *mevent;
    std::vector<double>::iterator itx, ity;

    Glib::RefPtr<Gdk::Display> rt_display = Gtk::Widget::get_display();
    Glib::RefPtr<Gdk::Screen> rt_screen = Gtk::Widget::get_screen();

    bool retval = false;

    /* width and height are the size of the graph */
    int width = get_allocation().get_width() - RADIUS * 2;
    int height = get_allocation().get_height() - RADIUS * 2;

    if ((width < 0) || (height < 0))
        return false;

    /*  get the pointer position  */
    int tx, ty;
    Gdk::ModifierType gm;
    get_window()->get_pointer (tx, ty, gm);
    int x =            CLAMP ((tx - RADIUS), 0, width-1);  // X position of the pointer from the origin of the graph
    int y = height-1 - CLAMP ((ty - RADIUS), 0, height-1); // Y position of the pointer from the origin of the graph

    unsigned int distance_x = ~0U,  distance_y = ~0U;
    int num = curve.x.size();
    int closest_point = -1;

    if (curve.type!=Parametric) {
        for (int i = 0; i < num; ++i) {
            int cx = (int)((width-1) * curve.x[i] + 0.5); //project (c->ctlpoint[i][0], min_x, c->max_x, width);
            int cy = (int)((height-1) * curve.y[i] + 0.5); //project (c->ctlpoint[i][0], min_x, c->max_x, width);
            unsigned int curr_dist_x = abs (x - cx);
            unsigned int curr_dist_y = abs (y - cy);
            if (curr_dist_x < distance_x) {
                distance_x = curr_dist_x;
                distance_y = curr_dist_y;
                closest_point = i;
            }
            else if (curr_dist_x == distance_x && curr_dist_y < distance_y) {
			// there is mode than 1 point for that X coordinate, we select the point closest to the cursor
				distance_y = curr_dist_y;
				closest_point = i;
            }
        }
    }

    switch (event->type) {
        case Gdk::CONFIGURE:
            if (pixmap)
                pixmap.clear ();

        case Gdk::EXPOSE:
        	// When does this event occurs ???
            if (!pixmap) {
	            pixmap = Gdk::Pixmap::create (get_window(), get_allocation().get_width(),  get_allocation().get_height());
	            interpolate (width, height);
	        }
            draw (width, height, lit_point);
            break;

        case Gdk::BUTTON_PRESS:
            if (curve.type!=Parametric) {
                add_modal_grab ();

                // get cursor position
            	Gdk::ModifierType mod_type;
                rt_display->get_pointer(cursor_x, cursor_y, mod_type);

                new_type = CSEmpty;
    	        if (distance_x > MIN_DISTANCE) {
	                /* insert a new control point */
	                if (num > 0) {
		                int cx = (int)((width-1)*curve.x[closest_point]+0.5);
		                if (x > cx)
		                    ++closest_point;
		            }
		            itx = curve.x.begin();
		            ity = curve.y.begin();
	                for (int i=0; i<closest_point; i++) { itx++; ity++; }
					curve.x.insert (itx, 0);
					curve.y.insert (ity, 0);
					num++;

					// the graph is refreshed only if a new point is created (snaped to a pixel)
					curve.x[closest_point] =  (double) x / (width-1);
					curve.y[closest_point] =  (double) y / (height-1);
					interpolate (width, height);
		            draw (width, height, closest_point);
    	        }
				grab_point = closest_point;
				lit_point = closest_point;
				ugp_x = curve.x[closest_point];
				ugp_y = curve.y[closest_point];
				notifyListener ();
				break;
	        }
            retval = true;
            break;

        case Gdk::BUTTON_RELEASE:
            if (curve.type!=Parametric) {
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
	                if (curve.x.size() <= 0) {
    		            curve.x.push_back (0);
    		            curve.y.push_back (0);
    		            interpolate (width, height);
    		            draw (width, height, lit_point);
                    }
	            }
	            if (distance_x <= MIN_DISTANCE) {
		            new_type = CSMove;
					lit_point = closest_point;
				}
	            else {
		            new_type = CSPlus;
		            lit_point = -1;
	            }
	            if (lit_point != previous_lit_point)
	            	draw (width, height, lit_point);
                grab_point = -1;
                retval = true;
                notifyListener ();
            }
            break;

        case Gdk::LEAVE_NOTIFY:
        	// Pointer can LEAVE even when dragging the point, so we don't modify the cursor in this case
        	// The cursor will have to LEAVE another time after the drag...
        	if (grab_point == -1)
        		new_type = CSArrow;
        	break;

        case Gdk::MOTION_NOTIFY:
            mevent = (GdkEventMotion *) event;

            if (curve.type == Linear || curve.type == Spline || curve.type == NURBS) {
    	        if (grab_point == -1) {
    	        	int previous_lit_point = lit_point;
    	            /* if no point is grabbed...  */
    	            if (distance_x <= MIN_DISTANCE) {
    		            new_type = CSMove;
						lit_point = closest_point;
					}
    	            else {
    		            new_type = CSPlus;
    		            lit_point = -1;
    	            }
    	            if ((new_type != cursor_type) || (lit_point != previous_lit_point))
    	            	draw (width, height, lit_point);
    	        }
    	        else  {
    	            int new_cursor_x, new_cursor_y;
    	        	double factor = 0.5;

    	            // get cursor position
    	        	Gdk::ModifierType mod_type;
    	            rt_display->get_pointer(new_cursor_x, new_cursor_y, mod_type);

    	            // set the dragging factor
    	        	int control_key = gm & GDK_CONTROL_MASK;
    	        	int shift_key = gm & GDK_SHIFT_MASK;

    	            // what is the speed factor
    	        	if (control_key && shift_key) factor = 0.005;
    	            else if (shift_key)           factor = 0.02;
    	            else if (control_key)         factor = 0.1;

    	            // calculate the delta in [0.0 ; 1.0] range
    	            double delta_x = (double)(new_cursor_x - cursor_x) * factor / (double)(width-1);
    	            double delta_y = (double)(cursor_y - new_cursor_y) * factor / (double)(height-1);

    	            // bounds of the grabed point
	            	double leftbound = (grab_point == 0) ? 0. : curve.x[grab_point-1];
	            	double rightbound = (grab_point == num-1) ? 1. : curve.x[grab_point+1];
	            	double bottombound = (double)(-MIN_DISTANCE) * factor / (double)(height-1);
	            	double topbound = (double)1.0 + (double)(MIN_DISTANCE) * factor / (double)(height-1);

    	            // modification of the unclamped grabed point
					bool delete_me = false;
					// Handling limitations along X axis
    	            if (ugp_x >= leftbound && ugp_x <= rightbound) {
    	            	ugp_x += delta_x;
    	            	if (ugp_x > rightbound) {
    	            		if (grab_point == num-1)
    	            			curve.x[grab_point] = 1.;
    	            		else
    	            			if (num == 2)
    	            				curve.x[grab_point] = rightbound;
    	            			else
    	            				curve.x[grab_point] = -1.;
    	            	}
    	            	else if (ugp_x < leftbound) {
    	            		if (grab_point == 0)
    	            			curve.x[grab_point] = 0.;
    	            		else
    	            			if (num == 2)
    	            				curve.x[grab_point] = leftbound;
    	            			else
    	            				curve.x[grab_point] = -1.;
    	            	}
    	            	else
    	            		curve.x[grab_point] = ugp_x;
    	            }
    	            else if (ugp_x > rightbound && delta_x < 0.)
						curve.x[grab_point] = ugp_x = rightbound;
					else if (ugp_x < leftbound && delta_x > 0.)
						curve.x[grab_point] = ugp_x = leftbound;

					// Handling limitations along Y axis
    	            if (ugp_y >= bottombound && ugp_y <= topbound) {
    	            	ugp_y += delta_y;
    	            	if (ugp_y > topbound) {
    	            		if (grab_point == 0 || grab_point == num-1)
    	            			curve.y[grab_point] = 1.;
    	            		else
   	            				curve.x[grab_point] = -1.;
    	            	}
    	            	else if (ugp_y < bottombound) {
    	            		if (grab_point == 0 || grab_point == num-1)
    	            			curve.y[grab_point] = 0.;
    	            		else
   	            				curve.x[grab_point] = -1.;
    	            	}
    	            	else
       	            		curve.y[grab_point] = CLAMP(ugp_y, 0.0, 1.0);
    	            }
    	            else if (ugp_y > 1. && delta_y < 0.)
						curve.y[grab_point] = ugp_y = 1.0;
					else if (ugp_y < 0. && delta_y > 0.)
						curve.y[grab_point] = ugp_y = 0.;
    	            else if ((grab_point > 0 && grab_point < num-1) && (ugp_y > topbound || ugp_y < bottombound))
           				curve.x[grab_point] = -1.;

    	            interpolate (width, height);

    	            // move the cursor back (to avoid being limited by the screen)
    	        	rt_display->warp_pointer(rt_screen, cursor_x, cursor_y);

    	            draw (width, height, lit_point);
                    notifyListener ();
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

std::vector<double> MyCurve::getPoints () {

    std::vector<double> result;
    if (curve.type==Parametric) {
        result.push_back ((double)(Parametric));
        for (int i=0; i<curve.x.size(); i++)
            result.push_back (curve.x[i]);
    }
    else {
    	// the first value gives the type of the curve
        if (curve.type==Linear)
            result.push_back ((double)(Linear));
        else if (curve.type==Spline)
            result.push_back ((double)(Spline));
        else if (curve.type==NURBS)
            result.push_back ((double)(NURBS));
        // then we push all the points coordinate
        for (int i=0; i<curve.x.size(); i++)
            if (curve.x[i]>=0) {
                result.push_back (curve.x[i]);
                result.push_back (curve.y[i]);
            }
    }
    return result;
}

void MyCurve::setPoints (const std::vector<double>& p) {

    int ix = 0;
    CurveType t = (CurveType)p[ix++];
    curve.type = t;
    if (t==Parametric) {
        curve.x.clear ();
        curve.y.clear ();
        for (int i=1; i<p.size(); i++)
            curve.x.push_back (p[ix++]);
    }
    else {
        curve.x.clear ();
        curve.y.clear ();
        for (int i=0; i<p.size()/2; i++) {
            curve.x.push_back (p[ix++]);
            curve.y.push_back (p[ix++]);
        }
        activeParam = -1;
    }
    pixmap.clear ();
    queue_draw ();
}

void MyCurve::setType (CurveType t) {

    curve.type = t;
    pixmap.clear ();
}

void MyCurve::notifyListener () {

    if (listener)
        listener->curveChanged ();    	      
}

void MyCurve::setActiveParam (int ac) {
    
    activeParam = ac;
    pixmap.clear ();
    queue_draw ();
}

int mchistupdate (void* data) {

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
    
    mcih->myCurve->pixmap.clear ();
    mcih->myCurve->queue_draw ();

    mcih->pending--;
    gdk_threads_leave ();
    return 0;
}

void MyCurve::updateBackgroundHistogram (unsigned int* hist) {

    if (hist!=NULL) {
        memcpy (bghist, hist, 256*sizeof(unsigned int));
        bghistvalid = true;
    }
    else
        bghistvalid = false;

    mcih->pending++;
    g_idle_add (mchistupdate, mcih);

}

void MyCurve::reset() {
	int width = get_allocation().get_width() - RADIUS * 2;
	int height = get_allocation().get_height() - RADIUS * 2;

	switch (curve.type) {
	case Spline :
	case  NURBS :
		curve.x.clear();
		curve.y.clear();
	    curve.x.push_back(0.);
	    curve.y.push_back(0.);
	    curve.x.push_back(1.);
	    curve.y.push_back(1.);
	    grab_point = -1;
	    lit_point = -1;
        interpolate (width, height);
		break;
	case Parametric :
		// Nothing to do (?)
	default:
		break;
	}
	draw(width, height, -1);
}
