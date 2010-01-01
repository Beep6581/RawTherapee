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
#include <histogrampanel.h>
#include <multilangmgr.h>

HistogramPanel::HistogramPanel () {

   histogramArea = Gtk::manage (new HistogramArea ());
   showRed   = Gtk::manage (new Gtk::ToggleButton ("R"));
   showGreen = Gtk::manage (new Gtk::ToggleButton ("G"));
   showBlue  = Gtk::manage (new Gtk::ToggleButton ("B"));
   showValue = Gtk::manage (new Gtk::ToggleButton ("L"));
   Gtk::VBox* vbox = Gtk::manage (new Gtk::VBox (false, 0));

   showRed->set_active (true);
   showGreen->set_active (true);
   showBlue->set_active (true);
   showValue->set_active (true);
   vbox->pack_start (*showRed, Gtk::PACK_SHRINK, 2);
   vbox->pack_start (*showGreen, Gtk::PACK_SHRINK, 2);
   vbox->pack_start (*showBlue, Gtk::PACK_SHRINK, 2);
   vbox->pack_start (*showValue, Gtk::PACK_SHRINK, 2);
   pack_start (*histogramArea);
   pack_end (*vbox, Gtk::PACK_SHRINK, 2);

   showRed->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::rgbv_toggled) );
   showGreen->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::rgbv_toggled) );
   showBlue->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::rgbv_toggled) );
   showValue->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::rgbv_toggled) );

   show_all ();

   showRed->set_tooltip_text   (M("HISTOGRAM_TOOLTIP_R"));
   showGreen->set_tooltip_text (M("HISTOGRAM_TOOLTIP_G"));
   showBlue->set_tooltip_text  (M("HISTOGRAM_TOOLTIP_B"));
   showValue->set_tooltip_text (M("HISTOGRAM_TOOLTIP_L"));

    rconn = signal_size_allocate().connect( sigc::mem_fun(*this, &HistogramPanel::resized) );
}

void HistogramPanel::resized (Gtk::Allocation& req) {

    rconn.block (true);
    
    if (req.get_width()/2>150)
        set_size_request (req.get_width(), 150);
    else
        set_size_request (req.get_width(), req.get_width()/2);
    rconn.block (false);
    histogramArea->renderHistogram ();
    histogramArea->queue_draw ();
}

void HistogramPanel::rgbv_toggled () {

  histogramArea->updateOptions (showRed->get_active(), showGreen->get_active(), showBlue->get_active(), showValue->get_active());
  histogramArea->queue_draw ();
}

HistogramArea::HistogramArea () : 
      needVal(true), needRed(true), needGreen(true), needBlue(true), oldwidth(-1), valid(false), showFull(true) {

    haih = new HistogramAreaIdleHelper;
    haih->harea = this;
    haih->destroyed = false;
    haih->pending = 0;
    
    signal_style_changed().connect( sigc::mem_fun(*this, &HistogramArea::styleChanged) );
}

HistogramArea::~HistogramArea () {

    if (haih->pending)
        haih->destroyed = true;
    else
        delete haih;
}

void HistogramArea::updateOptions (bool r, bool g, bool b, bool v) {

    needRed   = r;
    needGreen = g;
    needBlue  = b;
    needVal   = v;

    renderHistogram ();
}

int histupdate (void* data) {

    gdk_threads_enter ();

    HistogramAreaIdleHelper* haih = (HistogramAreaIdleHelper*)data;

    if (haih->destroyed) {
        if (haih->pending == 1)
            delete haih;
        else    
            haih->pending--;
        gdk_threads_leave ();
        return 0;
    }
    
    haih->harea->renderHistogram ();
    haih->harea->queue_draw ();

    haih->pending--;
    gdk_threads_leave ();
    return 0;
}

void HistogramArea::update (unsigned int* rh, unsigned int* gh, unsigned int* bh, unsigned int* lh) {

    if (rh!=NULL) {
        memcpy (lhist, lh, 256*sizeof(unsigned int));
        memcpy (rhist, rh, 256*sizeof(unsigned int));
        memcpy (ghist, gh, 256*sizeof(unsigned int));
        memcpy (bhist, bh, 256*sizeof(unsigned int));
        valid = true;
    }
    else
        valid = false;

    haih->pending++;
    g_idle_add (histupdate, haih);
}

void HistogramArea::renderHistogram () {

  if (!is_realized ())
        return;

  Glib::RefPtr<Gdk::Window> window = get_window();
  int winx, winy, winw, winh, wind;
  window->get_geometry(winx, winy, winw, winh, wind);

  backBuffer = Gdk::Pixmap::create (window, winw, winh, -1);

  Glib::RefPtr<Gdk::GC> bgc = Gdk::GC::create(backBuffer);
  
  bgc->set_foreground (white);
  backBuffer->draw_rectangle (bgc, true, 0, 0, winw, winh);

  if (valid) {

    // compute height of the full histogram (realheight) and

    int fullhistheight = 0;
    for (int i=0; i<256; i++) {
        if (needVal && lhist[i]>fullhistheight)
	        fullhistheight = lhist[i];
        if (needRed && rhist[i]>fullhistheight)
	        fullhistheight = rhist[i];
        if (needGreen && ghist[i]>fullhistheight)
	        fullhistheight = ghist[i];
        if (needBlue && bhist[i]>fullhistheight)
	        fullhistheight = bhist[i];
    }
    
    // compute two hights, one for the magnified view and one for the threshold
    int realhistheight = fullhistheight;

    if (!showFull) {    
        int area1thres = 0;
        int area2thres = 0;
        int area = 0;
        for (int i=0; i<fullhistheight; i++) {
            for (int j=0; j<256; j++)
                if ((needVal && lhist[j]>i) || (needRed && rhist[j]>i) || (needGreen && ghist[j]>i) || (needBlue && bhist[j]>i))
                    area++;
            if (area1thres==0 && (double)area / (256*(i+1)) < 0.3)
                area1thres = i;
            if (area2thres==0 && (double)area / (256*(i+1)) < 0.3)
                area2thres = i;
            if (area1thres && area2thres)
                break;
        }
        if (area1thres>0 && area2thres>0 && area1thres<fullhistheight)
            realhistheight = area2thres;
    }
    
    if (realhistheight<winh-2)
        realhistheight = winh-2;

    Cairo::RefPtr<Cairo::Context> cr = backBuffer->create_cairo_context();
    cr->set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
    cr->set_line_width (1.0);
    double stepSize = (winw-1) / 256.0;
    if (needVal) {
        cr->move_to (0, winh-1);
        cr->set_source_rgb (0.75, 0.75, 0.75);
        for (int i=0; i<256; i++) {
            double val = lhist[i] * (double)(winh-2) / realhistheight;
            if (val>winh-1)
                val = winh-1;
      	    if (i>0)
       	        cr->line_to (i*stepSize, winh-1-val);
    	}
        cr->save ();
        cr->line_to (winw-1, winh-1);    	cr->fill_preserve ();
        cr->restore ();
        cr->set_source_rgb (0.5, 0.5, 0.5);
    	cr->stroke ();
    }
    if (needRed) {
        cr->move_to (0, winh-1);
        cr->set_source_rgb (1.0, 0.0, 0.0);
        for (int i=0; i<256; i++) {
            double val = rhist[i] * (double)(winh-2) / realhistheight;
            if (val>winh-1)
                val = winh-1;
  	        if (i>0)
   	            cr->line_to (i*stepSize, winh-1-val);	
    	}
    	cr->stroke ();
    }
    if (needGreen) {        cr->move_to (0, winh-1);
        cr->set_source_rgb (0.0, 1.0, 0.0);
        for (int i=0; i<256; i++) {
            double val = ghist[i] * (double)(winh-2) / realhistheight;
            if (val>winh-1)
                val = winh-1;
  	        if (i>0)
   	            cr->line_to (i*stepSize, winh-1-val);	
    	}
    	cr->stroke ();
    }
    if (needBlue) {
        cr->move_to (0, winh-1);
        cr->set_source_rgb (0.0, 0.0, 1.0);
        for (int i=0; i<256; i++) {
            int val = bhist[i] * (double)(winh-2) / realhistheight;
            if (val>winh-1)
                val = winh-1;
      	    if (i>0)
       	        cr->line_to (i*stepSize, winh-1-val);
    	}
    	cr->stroke ();
    }
  }

/*

    // scale histogram to width winw-1

    int* vval   = new int[winw-1];
    int* vred   = new int[winw-1];
    int* vgreen = new int[winw-1];
    int* vblue  = new int[winw-1];

    memset (vval, 0, sizeof(int)*(winw-1));
    memset (vred, 0, sizeof(int)*(winw-1));
    memset (vgreen, 0, sizeof(int)*(winw-1));
    memset (vblue, 0, sizeof(int)*(winw-1));

    int index = 0;
    double scale = 256.0 / (winw-2);
    for (int i=0; i<=winw-2; i++) {
        int samples = 0;
        while (index < 256 && (int)(index/scale) == i) {
	        vval[i]    += lhist[index];
	        vred[i]    += rhist[index];
	        vgreen[i]  += ghist[index];
	        vblue[i]   += bhist[index];
	        index++;
	        samples++;
        }
        if (samples>0) {
            vval[i] /= samples;
            vred[i] /= samples;
            vgreen[i] /= samples;
            vblue[i] /= samples;
        }
    }

    // compute height of the full histogram (realheight) and

    int fullhistheight = 0;
    for (int i=0; i<=winw-2; i++) {
        if (needVal && vval[i]>fullhistheight)
	        fullhistheight = vval[i];
        if (needRed && vred[i]>fullhistheight)
	        fullhistheight = vred[i];
        if (needGreen && vgreen[i]>fullhistheight)
	        fullhistheight = vgreen[i];
        if (needBlue && vblue[i]>fullhistheight)
	        fullhistheight = vblue[i];
    }
    
    // compute two hights, one for the magnified view and one for the threshold

    int realhistheight = fullhistheight;

    if (!showFull) {    
        int area1thres = 0;
        int area2thres = 0;
        int area = 0;
        for (int i=0; i<fullhistheight; i++) {
            for (int j=0; j<winw-1; j++)
                if ((needVal && vval[j]>i) || (needRed && vred[j]>i) || (needGreen && vgreen[j]>i) || (needBlue && vblue[j]>i))
                    area++;
            if (area1thres==0 && (double)area / ((winw-1)*(i+1)) < 0.3)
                area1thres = i;
            if (area2thres==0 && (double)area / ((winw-1)*(i+1)) < 0.3)
                area2thres = i;
            if (area1thres && area2thres)
                break;
        }
        if (area1thres>0 && area2thres>0 && area1thres<fullhistheight)
            realhistheight = area2thres;
    }
    
    if (realhistheight<winh-2)
        realhistheight = winh-2;

    Cairo::RefPtr<Cairo::Context> cr = backBuffer->create_cairo_context();
    cr->set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
    cr->set_line_width (1.0);
    if (needVal) {
        cr->move_to (0, winh-1);
        cr->set_source_rgb (0.75, 0.75, 0.75);
        for (int i=0; i<=winw-2; i++) {
            int val = (int)(vval[i] * (double)(winh-2) / realhistheight);
            if (val>winh-1)
                val = winh-1;
      	    if (i>0)
       	        cr->line_to (i+1, winh-1-val);
    	}
    	cr->fill_preserve ();
        cr->set_source_rgb (0.5, 0.5, 0.5);
    	cr->stroke ();
    }
    if (needRed) {
        cr->move_to (0, winh-1);
        cr->set_source_rgb (1.0, 0.0, 0.0);
        for (int i=0; i<=winw-2; i++) {
            int val = (int)(vred[i] * (double)(winh-2) / realhistheight);
            if (val>winh-1)
                val = winh-1;
  	        if (i>0)
   	            cr->line_to (i+1, winh-1-val);	
    	}
    	cr->stroke ();
    }
    if (needGreen) {
        cr->move_to (0, winh-1);
        cr->set_source_rgb (0.0, 1.0, 0.0);
        for (int i=0; i<=winw-2; i++) {
            int val = (int)(vgreen[i] * (double)(winh-2) / realhistheight);
            if (val>winh-1)
                val = winh-1;
  	        if (i>0)
   	            cr->line_to (i+1, winh-1-val);	
    	}
    	cr->stroke ();
    }
    if (needBlue) {
        cr->move_to (0, winh-1);
        cr->set_source_rgb (0.0, 0.0, 1.0);
        for (int i=0; i<=winw-2; i++) {
            int val = (int)(vblue[i] * (double)(winh-2) / realhistheight);
            if (val>winh-1)
                val = winh-1;
      	    if (i>0)
       	        cr->line_to (i+1, winh-1-val);
    	}
    	cr->stroke ();
    }
    
    delete [] vval;
    delete [] vred;
    delete [] vgreen;
    delete [] vblue;
  }
*/
  bgc->set_foreground (mgray);
  backBuffer->draw_rectangle (bgc, false, 0, 0, winw-1, winh-1);

  bgc->set_line_attributes (1, Gdk::LINE_ON_OFF_DASH, Gdk::CAP_NOT_LAST, Gdk::JOIN_MITER);

  backBuffer->draw_line (bgc, winw/3, 0, winw/3, winh-1);
  backBuffer->draw_line (bgc, 2*winw/3, 0, 2*winw/3, winh-1);
  backBuffer->draw_line (bgc, 0, winh/3, winw-1, winh/3);
  backBuffer->draw_line (bgc, 0, 2*winh/3, winw-1, 2*winh/3);

  bgc->set_line_attributes (1, Gdk::LINE_SOLID, Gdk::CAP_NOT_LAST, Gdk::JOIN_MITER);

  oldwidth = winw;
  oldheight = winh;
}

void HistogramArea::on_realize () {

  Gtk::DrawingArea::on_realize();
  Glib::RefPtr<Gdk::Window> window = get_window();
  gc_ = Gdk::GC::create(window);
  add_events(Gdk::EXPOSURE_MASK | Gdk::BUTTON_PRESS_MASK);
  Glib::RefPtr<Gdk::Colormap> colormap = get_default_colormap();

  black = Gdk::Color ("black");
  red = Gdk::Color ("red");
  green = Gdk::Color ("green");
  blue = Gdk::Color ("blue");
  lgray = Gdk::Color ("gray75");
  mgray = Gdk::Color ("gray50");
  dgray = Gdk::Color ("gray25");
  colormap->alloc_color(black);
  colormap->alloc_color(white);
  colormap->alloc_color(red);
  colormap->alloc_color(green);
  colormap->alloc_color(blue);
  colormap->alloc_color(lgray);
  colormap->alloc_color(mgray);
  colormap->alloc_color(dgray);

}

void HistogramArea::styleChanged (const Glib::RefPtr<Gtk::Style>& style) {

  white = get_style()->get_base(Gtk::STATE_NORMAL);
  queue_draw ();
}

bool HistogramArea::on_expose_event(GdkEventExpose* event) {

  Glib::RefPtr<Gdk::Window> window = get_window();

  int winx, winy, winw, winh, wind;
  window->get_geometry(winx, winy, winw, winh, wind);

  if (winw!=oldwidth && winh!=oldheight)
    renderHistogram ();
  window->draw_drawable (gc_, backBuffer, 0, 0, 0, 0, -1, -1);

  return true;
}


bool HistogramArea::on_button_press_event (GdkEventButton* event) {

    if (event->type==GDK_2BUTTON_PRESS && event->button==1) {
        showFull = !showFull;
        renderHistogram ();
        queue_draw ();
    }
}


