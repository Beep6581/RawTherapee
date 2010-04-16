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

#define RADIUS		3	/* radius of the control points */
#define MIN_DISTANCE	8	/* min distance between control points */

MyCurve::MyCurve () : listener(NULL), activeParam(-1), bghistvalid(false) {

    cursor_type = Gdk::TOP_LEFT_ARROW;
    curve.type = Spline;
    height = 0;
    grab_point = -1;

    add_events(Gdk::EXPOSURE_MASK |	Gdk::POINTER_MOTION_MASK |	Gdk::POINTER_MOTION_HINT_MASK |	Gdk::ENTER_NOTIFY_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK | Gdk::BUTTON1_MOTION_MASK);
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
    rtengine::Curve* rtcurve = new rtengine::Curve (curveDescr);
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

void MyCurve::draw (int width, int height) {

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
        cr->move_to (0, height-1);
        cr->set_source_rgb (0.75, 0.75, 0.75);
        for (int i=0; i<256; i++) {
            double val = bghist[i] * (double)(height-2) / histheight;
            if (val>height-1)
                val = height-1;
      	    if (i>0)
       	        cr->line_to (i*stepSize, height-1-val);
    	}
        cr->line_to (width-1, height-1);
    	cr->fill ();
    }

    // draw the grid lines:
    cr->set_line_width (1.0);  
    c = style->get_dark (state);
    cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
    cr->set_antialias (Cairo::ANTIALIAS_NONE);
    for (int i = 0; i < 5; i++) {
        cr->move_to (RADIUS, i * height / 4 + RADIUS);
        cr->line_to (width + RADIUS, i * height / 4 + RADIUS);
        cr->move_to (i * width / 4 + RADIUS, RADIUS);
        cr->line_to (i * width / 4 + RADIUS, height + RADIUS);  
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

    // draw curve
    cr->set_source_rgb (0.0, 0.0, 0.0);
    cr->move_to (point[0].get_x(), point[0].get_y());
    for (int i=1; i<point.size(); i++)
        cr->line_to (point[i].get_x(), point[i].get_y());
    cr->stroke ();

    // draw bullets
    if (curve.type!=Parametric) 
        for (int i = 0; i < curve.x.size(); ++i) {
            double x = ((width-1) * curve.x[i] + 0.5)+RADIUS;    // project (curve.x[i], 0, 1, width);
            double y = height - ((height-1) * curve.y[i] + 0.5)+RADIUS; // project (curve.y[i], 0, 1, height);

            cr->arc (x, y, RADIUS, 0, 2*M_PI);
            cr->fill ();
        }

    get_window()->draw_drawable (style->get_fg_gc (state), pixmap, 0, 0, 0, 0, width + RADIUS * 2, height + RADIUS * 2);
}

bool MyCurve::handleEvents (GdkEvent* event) {

    Gdk::CursorType new_type = cursor_type;
    int src, dst;
    GdkEventMotion *mevent;
    std::vector<double>::iterator itx, ity;

    bool retval = false;

    int width = get_allocation().get_width() - RADIUS * 2;
    int height = get_allocation().get_height() - RADIUS * 2;

    if ((width < 0) || (height < 0))
        return false;

    /*  get the pointer position  */
    int tx, ty;
    Gdk::ModifierType gm;
    get_window()->get_pointer (tx, ty, gm);
    int x = CLAMP ((tx - RADIUS), 0, width-1);
    int y = CLAMP ((ty - RADIUS), 0, height-1);

    unsigned int distance = ~0U;
    int num = curve.x.size();
    int closest_point = 0;
    
    if (curve.type!=Parametric) {
        for (int i = 0; i < num; ++i) {
            int cx = (int)((width-1) * curve.x[i] + 0.5); //project (c->ctlpoint[i][0], min_x, c->max_x, width);
            if ((unsigned int) abs (x - cx) < distance) {
                distance = abs (x - cx);
                closest_point = i;
            }
        }
    }
    
    switch (event->type) {
        case Gdk::CONFIGURE:
            if (pixmap)
                pixmap.clear ();

        case Gdk::EXPOSE:
            if (!pixmap) {
	            pixmap = Gdk::Pixmap::create (get_window(), get_allocation().get_width(),  get_allocation().get_height());
	            interpolate (width, height);
	        }
            draw (width, height);
            break;

        case Gdk::BUTTON_PRESS:
            if (curve.type!=Parametric) {
                add_modal_grab ();
                new_type = Gdk::PLUS;
    	        if (distance > MIN_DISTANCE) {
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
	                }
	                grab_point = closest_point;
	                curve.x[grab_point] =  (double) x / (width-1);
	                curve.y[grab_point] =  (double) (height-y) / (height-1);
    	            interpolate (width, height);
                    notifyListener ();
	                break;
	        }
            draw (width, height);
            retval = true;
            break;

        case Gdk::BUTTON_RELEASE:
            if (curve.type!=Parametric) {
                remove_modal_grab ();
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
    		            draw (width, height);
                    }
	            }
                new_type = Gdk::FLEUR;
                grab_point = -1;
                retval = true;
                notifyListener ();
            }
            break;

        case Gdk::MOTION_NOTIFY:
            mevent = (GdkEventMotion *) event;
            if (curve.type == Linear || curve.type == Spline) {
    	        if (grab_point == -1) {
    	            /* if no point is grabbed...  */
    	            if (distance <= MIN_DISTANCE)
    		            new_type = Gdk::FLEUR;
    	            else
    		            new_type = Gdk::PLUS;
    	        }
    	        else  {
    	            /* drag the grabbed point  */
    	            new_type = Gdk::FLEUR;
    	            int leftbound = -MIN_DISTANCE;
    	            if (grab_point > 0)
    		            leftbound = (int)((width-1)*curve.x[grab_point-1]+0.5);
    
    	            int rightbound = width + RADIUS * 2 + MIN_DISTANCE;
    	            if (grab_point + 1 < num)
    		            rightbound = (int)((width-1)*curve.x[grab_point+1]+0.5);
    
    	            if (tx <= leftbound || tx >= rightbound || ty > height + RADIUS * 2 + MIN_DISTANCE || ty < -MIN_DISTANCE)
    		            curve.x[grab_point] = -1.0;
    	            else	{
    		            curve.x[grab_point] = (double) x / (width-1); 
    		            curve.y[grab_point] = (double) (height-y) / (height-1);
                    }
    	            interpolate (width, height);
    	            draw (width, height);
                    notifyListener ();
    	        }
    	    }
	
        if (new_type != cursor_type) {
	        cursor_type = new_type;
            Gdk::Cursor* cursor = new Gdk::Cursor (get_display(), cursor_type);
	        get_window ()->set_cursor (*cursor);
            delete cursor;
	    }
        retval = true;
        break;

    default:
      break;
  }
  return retval;
}

std::vector<double> MyCurve::getPoints () {

    std::vector<double> result;
    if (curve.type==Parametric) {
        result.push_back (+2.0);
        for (int i=0; i<curve.x.size(); i++)
            result.push_back (curve.x[i]);
    }
    else {
        if (curve.type==Linear)
            result.push_back (-1.0);
        else
            result.push_back (+1.0);
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
    int t = p[ix++];
    if (t==2) {
        curve.type = Parametric;
        curve.x.clear ();
        curve.y.clear ();
        for (int i=1; i<p.size(); i++)
            curve.x.push_back (p[ix++]);
    }
    else {
        if (t==1)
            curve.type = Spline;
        else
            curve.type = Linear;
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

