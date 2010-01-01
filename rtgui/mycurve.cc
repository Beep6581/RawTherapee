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

#define RADIUS		3	/* radius of the control points */
#define MIN_DISTANCE	8	/* min distance between control points */

MyCurve::MyCurve () : listener(NULL) {

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
}

void MyCurve::spline_solve (int n, double x[], double y[], double y2[]) {

  double* u = new double[n-1];

  y2[0] = u[0] = 0.0;	/* set lower boundary condition to "natural" */

  for (int i = 1; i < n - 1; ++i)
    {
      double sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
      double p = sig * y2[i - 1] + 2.0;
      y2[i] = (sig - 1.0) / p;
      u[i] = ((y[i + 1] - y[i])
	      / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]));
      u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }

  y2[n - 1] = 0.0;
  for (int k = n - 2; k >= 0; --k)
    y2[k] = y2[k] * y2[k + 1] + u[k];

  delete [] u;
}

double MyCurve::spline_eval (int n, double x[], double y[], double y2[], double val) {

  if (val>x[n-1])
      return y[n-1];
  else if (val<x[0])
      return y[0];

  /* do a binary search for the right interval: */
  int k_lo = 0, k_hi = n - 1;
  while (k_hi - k_lo > 1){
      int k = (k_hi + k_lo) / 2;
      if (x[k] > val)
	    k_hi = k;
      else
	    k_lo = k;
  }

  double h = x[k_hi] - x[k_lo];
  double a = (x[k_hi] - val) / h;
  double b = (val - x[k_lo]) / h;
  return a*y[k_lo] + b*y[k_hi] + ((a*a*a - a)*y2[k_lo] + (b*b*b - b)*y2[k_hi]) * (h*h)/6.0;
}


std::vector<double> MyCurve::get_vector (int veclen) {

  std::vector<double> vector;
  
  int num = curve.x.size(); 

  /* count active points: */
  double prev =- 1.0;
  int active = 0;
  int firstact = -1;
  for (int i = 0; i < num; ++i)
	if (curve.x[i] > prev) {
	    if (firstact < 0)
	      firstact = i;
	    prev = curve.x[i];
	    ++active;
	}
  /* handle degenerate case: */
  if (active < 2) {
      double ry;
	  if (active > 0)
	    ry = curve.y[firstact];
	  else
	    ry = 0.0;
	  if (ry < 0.0) ry = 0.0;
	  if (ry > 1.0) ry = 1.0;
	  for (int x = 0; x < veclen; ++x)
	    vector.push_back(ry);
	  return vector;
  }

  if (curve.type==Spline) {

      double* mem = new double [3*active];
      double* xv  = mem;
      double* yv  = mem + active;
      double* y2v = mem + 2*active;

      prev = -1.0;
      int dst = 0;
      for (int i = 0; i < num; ++i) {
    	if (curve.x[i] > prev) {
  	      prev  = curve.x[i];
 	      xv[dst] = curve.x[i];
	      yv[dst] = curve.y[i];
	      dst++;
	    }
	  }
      spline_solve (active, xv, yv, y2v);

      double dx = 1.0 / (veclen - 1);
      double rx = 0.0;
      for (int x = 0; x < veclen; ++x, rx += dx) {
    	  double ry = spline_eval (active, xv, yv, y2v, rx);
	      if (ry < 0.0) ry = 0;
	      if (ry > 1.0) ry = 1.0;
    	  vector.push_back (ry);
      }
      delete [] mem;
  }
  else if (curve.type==Linear) {
      double dx = 1.0 / (veclen - 1);
      double rx = 0;
      double ry = 0;
      double dy = 0.0;
      int i  = firstact;
      for (int x = 0; x < veclen; ++x, rx += dx) {
	    if (rx >= curve.x[i]) {
	      if (rx > curve.x[i])
		    ry = 0.0;
	      dy = 0.0;
	      int next = i + 1;
	      while (next < num && curve.x[next] <= curve.x[i])
	        ++next;
	      if (next < num) {
		    double delta_x = curve.x[next] - curve.x[i];
		    dy = (curve.y[next] - curve.y[i]) / delta_x;
		    dy *= dx;
		    ry = curve.y[i];
		    i = next;
		  }
	    }
	    if (rx<curve.x[0])
	        vector.push_back (curve.y[0]);
	    else if (rx>curve.x[num-1])
	        vector.push_back (curve.y[num-1]);
	    else
    	    vector.push_back (ry);
	    ry += dy;
	  }
  }
  return vector;
}

void MyCurve::interpolate (int width, int height) {

  this->height = height;
  point.clear ();
  std::vector<double> vector = get_vector (width);
  this->height = height;
  for (int i = 0; i < width; ++i) {
    Gdk::Point p (RADIUS + i, RADIUS + height - (int)((height-1) * vector[i] + 0.5));
    point.push_back (p);
  }
}


void MyCurve::draw (int width, int height) {

  if (!pixmap)
    return;

  if (this->height != height || point.size() != width)
    interpolate (width, height);

  Gtk::StateType state = Gtk::STATE_NORMAL;
  if (!is_sensitive())
    state = Gtk::STATE_INSENSITIVE;

  Glib::RefPtr<Gtk::Style> style = get_style ();

  Cairo::RefPtr<Cairo::Context> cr = pixmap->create_cairo_context();

  /* clear the pixmap: */
//  gtk_paint_flat_box (style->gobj(), pixmap->gobj(), GTK_STATE_NORMAL, GTK_SHADOW_NONE,
//		      NULL, (GtkWidget*)gobj(), "curve_bg", 0, 0, , height + RADIUS * 2);

//  pixmap->draw_rectangle (style->get_bg_gc (state), false, 0, 0, width + RADIUS*2 - 1, height + RADIUS*2 - 1);

  Gdk::Color c = style->get_bg (state);
  cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
  cr->rectangle (0, 0, width + RADIUS*2, height + RADIUS*2);
  cr->fill ();

  /* draw the grid lines: (XXX make more meaningful) */
  cr->set_line_width (1.0);  
  c = style->get_dark (state);
  cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
  cr->set_antialias (Cairo::ANTIALIAS_NONE);
  for (int i = 0; i < 5; i++) {
  
      cr->move_to (RADIUS, i * height / 4 + RADIUS);
      cr->line_to (width + RADIUS, i * height / 4 + RADIUS);
      cr->move_to (i * width / 4 + RADIUS, RADIUS);
      cr->line_to (i * width / 4 + RADIUS, height + RADIUS);  
//      pixmap->draw_line (style->get_dark_gc (state), RADIUS, i * height / 4 + RADIUS, width + RADIUS, i * height / 4 + RADIUS);
//      pixmap->draw_line (style->get_dark_gc (state), i * width / 4 + RADIUS, RADIUS, i * width / 4 + RADIUS, height + RADIUS);
  }
  cr->stroke ();

  cr->set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
  cr->set_line_width (1.0);
  cr->set_source_rgb (0.0, 0.0, 0.0);
  cr->move_to (point[0].get_x(), point[0].get_y());
  for (int i=1; i<point.size(); i++)
      cr->line_to (point[i].get_x(), point[i].get_y());
  cr->stroke ();

  for (int i = 0; i < curve.x.size(); ++i) {

	double x = ((width-1) * curve.x[i] + 0.5)+RADIUS;    // project (curve.x[i], 0, 1, width);
	double y = height - ((height-1) * curve.y[i] + 0.5)+RADIUS; // project (curve.y[i], 0, 1, height);

	/* draw a bullet: */
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
  for (int i = 0; i < num; ++i) {
      int cx = (int)((width-1) * curve.x[i] + 0.5); //project (c->ctlpoint[i][0], min_x, c->max_x, width);
      if ((unsigned int) abs (x - cx) < distance) {
	    distance = abs (x - cx);
	    closest_point = i;
	  }
  }

  switch (event->type) {
    case Gdk::CONFIGURE:
      if (pixmap)
        pixmap.clear ();
      /* fall through */

    case Gdk::EXPOSE:
      if (!pixmap) {
	    pixmap = Gdk::Pixmap::create (get_window(), get_allocation().get_width(),  get_allocation().get_height());
	    interpolate (width, height);
	  }
      draw (width, height);
	  
      break;

    case Gdk::BUTTON_PRESS:
      add_modal_grab ();
      new_type = Gdk::PLUS;
      switch (curve.type) {
    	case Linear:
    	case Spline:
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
      remove_modal_grab ();

      /* delete inactive points: */
      itx = curve.x.begin();
      ity = curve.y.begin();
	  for (src = dst = 0; src < num; ++src) {
	      if (curve.x[src] >= 0.0) {
            curve.x[dst] = curve.x[src];
            curve.y[dst] = curve.y[src];
  		    ++dst;
  		    ++itx;
  		    ++ity;
		  }
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
      break;

    case Gdk::MOTION_NOTIFY:
      mevent = (GdkEventMotion *) event;

      switch (curve.type) {
	    case Linear:
    	case Spline:
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
    	  break;
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
    if (curve.type==Linear)
        result.push_back (-1.0);
    else
        result.push_back (+1.0);
    for (int i=0; i<curve.x.size(); i++)
        if (curve.x[i]>=0) {
            result.push_back (curve.x[i]);
            result.push_back (curve.y[i]);
        }
    return result;
}

void MyCurve::setPoints (const std::vector<double>& p) {

    int ix = 0;
    if (p[ix++]>0)
        curve.type = Spline;
    else
        curve.type = Linear;

    curve.x.clear ();
    curve.y.clear ();
    for (int i=0; i<p.size()/2; i++) {
        curve.x.push_back (p[ix++]);
        curve.y.push_back (p[ix++]);
    }
    pixmap.clear ();
    bool pi = pixmap;
}

void MyCurve::setType (CurveType t) {

    curve.type = t;
    pixmap.clear ();
}

void MyCurve::notifyListener () {
    if (listener)
        listener->curveChanged ();    	      
}
