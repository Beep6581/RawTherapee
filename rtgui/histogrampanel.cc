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
#include <string.h>
#include <LUT.h>

//
//
// HistogramPanel
HistogramPanel::HistogramPanel () {

   // Gtk::HBox* outer_hbox = Gtk::manage (new Gtk::HBox (false, 0));
	   
   histogramArea = Gtk::manage (new HistogramArea ());
   histogramRGBArea = Gtk::manage (new HistogramRGBArea ());

   Gtk::VBox* left_vbox = Gtk::manage (new Gtk::VBox (false, 0));
   left_vbox->pack_start (*histogramArea, Gtk::PACK_EXPAND_WIDGET, 2);
   left_vbox->pack_start (*histogramRGBArea, Gtk::PACK_SHRINK, 2);
   set_size_request (-1, 170);
   histogramArea->set_size_request (-1, 150);
   histogramRGBArea->set_size_request(-1, 10);

   showRed   = Gtk::manage (new Gtk::ToggleButton (M("HISTOGRAM_BUTTON_R")));
   showGreen = Gtk::manage (new Gtk::ToggleButton (M("HISTOGRAM_BUTTON_G")));
   showBlue  = Gtk::manage (new Gtk::ToggleButton (M("HISTOGRAM_BUTTON_B")));
   showValue = Gtk::manage (new Gtk::ToggleButton (M("HISTOGRAM_BUTTON_L")));
   showRAW   = Gtk::manage (new Gtk::ToggleButton (M("HISTOGRAM_BUTTON_RAW")));
   showBAR   = Gtk::manage (new Gtk::ToggleButton (M("HISTOGRAM_BUTTON_BAR")));

   showRed->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::rgbv_toggled) );
   showGreen->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::rgbv_toggled) );
   showBlue->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::rgbv_toggled) );
   showValue->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::rgbv_toggled) );
   showRAW->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::rgbv_toggled) );
   showBAR->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::rgbv_toggled) );

   showRed->set_tooltip_text   (M("HISTOGRAM_TOOLTIP_R"));
   showGreen->set_tooltip_text (M("HISTOGRAM_TOOLTIP_G"));
   showBlue->set_tooltip_text  (M("HISTOGRAM_TOOLTIP_B"));
   showValue->set_tooltip_text (M("HISTOGRAM_TOOLTIP_L"));
   showRAW->set_tooltip_text   (M("HISTOGRAM_TOOLTIP_RAW"));
   showBAR->set_tooltip_text   (M("HISTOGRAM_TOOLTIP_BAR"));

   Gtk::VBox* vbox = Gtk::manage (new Gtk::VBox (false, 0));
   showRed->set_active (true);
   showGreen->set_active (true);
   showBlue->set_active (true);
   showValue->set_active (true);
   showRAW->set_active (false);
   showBAR->set_active (true);
   vbox->pack_start (*showRed, Gtk::PACK_SHRINK, 2);
   vbox->pack_start (*showGreen, Gtk::PACK_SHRINK, 2);
   vbox->pack_start (*showBlue, Gtk::PACK_SHRINK, 2);
   vbox->pack_start (*showValue, Gtk::PACK_SHRINK, 2);
   vbox->pack_start (*showRAW, Gtk::PACK_SHRINK, 2);
   vbox->pack_start (*showBAR, Gtk::PACK_SHRINK, 2);

   pack_start (*left_vbox,Gtk::PACK_EXPAND_WIDGET, 2);
   pack_start (*vbox, Gtk::PACK_SHRINK, 2);

   show_all ();

   rconn = signal_size_allocate().connect( sigc::mem_fun(*this, &HistogramPanel::resized) );
}

void HistogramPanel::resized (Gtk::Allocation& req) {

    rconn.block (true);
    
    if (req.get_width()/2>150) {
        set_size_request (req.get_width(), 170);
        histogramArea->set_size_request (req.get_width(), 150);
        histogramRGBArea->set_size_request (req.get_width(), 10);
	
	// Probably set R,G,B,V, RAW, and BAR button here to their original size if resizing below is implemented.
    }
    else {
        set_size_request (req.get_width(), req.get_width()/2);
        histogramArea->set_size_request (req.get_width(), req.get_width()/2);
        histogramRGBArea->set_size_request(req.get_width(),5);

	// Probably reduce R,G,B,V, RAW, and BAR button here a to a half-sized version.
    }

    rconn.block (false);

    histogramArea->renderHistogram ();
    histogramArea->queue_draw ();

    if (histogramRGBArea->getFreeze()==true) {
        histogramRGBArea->updateFreeze(false);
        // set histogramRGBArea invalid;
        histogramRGBArea->renderRGBMarks(-1, -1, -1);
        // re-set freeze to old state
        histogramRGBArea->updateFreeze(true);
        histogramRGBArea->queue_draw ();
    }
    else {
        // set histogramRGBArea invalid;
        histogramRGBArea->renderRGBMarks(-1, -1, -1);
        histogramRGBArea->queue_draw ();
    }
}

void HistogramPanel::rgbv_toggled () {

  // Update Display
  histogramArea->updateOptions (showRed->get_active(), showGreen->get_active(), showBlue->get_active(), showValue->get_active(), showRAW->get_active());
  histogramArea->queue_draw ();

  histogramRGBArea->updateOptions (showRed->get_active(), showGreen->get_active(), showBlue->get_active(), showValue->get_active(), showRAW->get_active(), showBAR->get_active());
  histogramRGBArea->renderRGBMarks (0,0,0);
  histogramArea->queue_draw ();
}


void HistogramPanel::setHistRGBInvalid () {
        // do something to un-show vertical bars
        histogramRGBArea->renderRGBMarks(-1, -1, -1);
        histogramRGBArea->queue_draw (); 
        }
        
// "Freeze" is not a button, but a RMB-click, so this is not in the RGBV-Toggle method 
void HistogramPanel::toggleFreeze () {
        if (histogramRGBArea->getFreeze()==true) { histogramRGBArea->updateFreeze(false); }
        else {
	   if (histogramRGBArea->getShow()==true) {
		histogramRGBArea->updateFreeze(true);
	   }
	} 
        return;
}

void HistogramPanel::pointerMoved (bool validPos, Glib::ustring profile, int x, int y, int r, int g, int b) {

        if (!validPos) {
                // do something to un-show vertical bars
               histogramRGBArea->renderRGBMarks(-1, -1, -1);
               histogramRGBArea->queue_draw ();
        }
        else {
               // do something to show vertical bars
              histogramRGBArea->renderRGBMarks(r, g, b);
              histogramRGBArea->queue_draw ();
        }
}

//
//
//
// HistogramRGBArea
HistogramRGBArea::HistogramRGBArea () :
  frozen(false), valid(false), showMode(true), rawMode(false), needLuma(true), needRed(true), needGreen(true), needBlue(true) {

  harih = new HistogramRGBAreaIdleHelper;
  harih->harea = this;
  harih->destroyed = false;
  harih->pending = 0;

  signal_style_changed().connect( sigc::mem_fun(*this, &HistogramRGBArea::styleChanged) );
}

HistogramRGBArea::~HistogramRGBArea () {

    if (harih->pending)
         harih->destroyed = true;
    else
         delete harih;
}

bool HistogramRGBArea::getFreeze() {
    return(frozen);
}

bool HistogramRGBArea::getShow() {
    return(showMode);
}

void HistogramRGBArea::updateFreeze (bool f) {
    frozen = f;
    return;
}

void HistogramRGBArea::renderRGBMarks (int r, int g, int b) {

  if (!is_realized ())
        return;
  
  if (frozen)  {
        return;
  }

  Glib::RefPtr<Gdk::Window> window = get_window();
  int winx, winy, winw, winh, wind;
  window->get_geometry(winx, winy, winw, winh, wind);

  overlay = Gdk::Pixmap::create (window, winw, winh, -1);
  Glib::RefPtr<Gdk::GC> ovrl = Gdk::GC::create(overlay);

  Glib::RefPtr<Gtk::Style> style = get_style ();

  if (!showMode)  {
    ovrl->set_foreground (style->get_bg (Gtk::STATE_NORMAL));
    overlay->draw_rectangle (ovrl, true, 0, 0, winw, winh);
    if (rgbgc_ && overlay) {
      window->draw_drawable (rgbgc_, overlay, 0, 0, 0, 0, -1, -1); }
    return; }
  else {
    ovrl->set_foreground (style->get_fg (Gtk::STATE_NORMAL));
    overlay->draw_rectangle (ovrl, true, 0, 0, winw, winh);
    if (rgbgc_ && overlay) {
      window->draw_drawable (rgbgc_, overlay, 0, 0, 0, 0, -1, -1); }
  }

  Cairo::RefPtr<Cairo::Context> cr = overlay->create_cairo_context();
  cr->set_line_width (1.0);

  if ( r != -1 && g != -1 && b != -1 ) {
    if (needRed) {
      // Red
      cr->set_source_rgb(1.0, 0.0, 0.0);
      cr->move_to((int)(r*(winw/256.0)), 0);
      cr->line_to((int)(r*(winw/256.0)), winh-0);
      cr->stroke();
    }
    if (needGreen) {
      // Green
      cr->set_source_rgb(0.0, 1.0, 0.0);
      cr->move_to((int)(g*(winw/256.0)), 0);
      cr->line_to((int)(g*(winw/256.0)), winh-0);
      cr->stroke();
  }
    if (needBlue) {
      // Blue
      cr->set_source_rgb(0.0, 0.0, 1.0);
      cr->move_to((int)(b*(winw/256.0)), 0);
      cr->line_to((int)(b*(winw/256.0)), winh-0);
      cr->stroke();
  }
    if (needLuma) {
      // Luma
      cr->set_source_rgb(1.0, 1.0, 1.0);
      cr->move_to((int)((r+g+b)/3*(winw/256.0)), 0);
      cr->line_to((int)((r+g+b)/3*(winw/256.0)), winh-0);
      cr->stroke();
    }
  }
}

int histrgbupdate (void* data) {

    gdk_threads_enter ();

    HistogramRGBAreaIdleHelper* harih = (HistogramRGBAreaIdleHelper*)data;

    if (harih->destroyed) {
        if (harih->pending == 1)
            delete harih;
        else    
            harih->pending--;
        gdk_threads_leave ();
        return 0;
    }
    
    harih->harea->renderRGBMarks(-1,-1,-1);
    harih->harea->queue_draw ();

    harih->pending--;
    gdk_threads_leave ();

    return 0;
}

void HistogramRGBArea::update (int valh, int rh, int  gh, int bh) {
	
    if (valh) {
        val=valh;
        r=rh;
        g=gh;
        b=bh;
        valid = true;
    }
    else
        valid = false;
	
    harih->pending++;
    g_idle_add (histrgbupdate, harih);
}

void HistogramRGBArea::updateOptions (bool r, bool g, bool b, bool l, bool raw, bool show) {

    needRed   = r;
    needGreen = g;
    needBlue  = b;
    needLuma  = l;
    rawMode   = raw;
    showMode   = show;

    // Histogram RGB BAR button logic goes here

    // Disable bar button when RAW histogram is displayed
    if ( rawMode && showMode) {
        showMode = false;
   }

    // When un-showing the bar, set the freeze state to off
    if (!showMode) {
        updateFreeze(false); 
    }
}

void HistogramRGBArea::on_realize () {

  Gtk::DrawingArea::on_realize();
  Glib::RefPtr<Gdk::Window> window = get_window();
  rgbgc_ = Gdk::GC::create(window);
  add_events(Gdk::EXPOSURE_MASK | Gdk::BUTTON_PRESS_MASK);

  Glib::RefPtr<Gdk::Colormap> rgbcolormap = get_default_colormap();
  black = Gdk::Color ("black");
  red = Gdk::Color ("red");
  green = Gdk::Color ("green");
  blue = Gdk::Color ("blue");
  lgray = Gdk::Color ("gray75");
  mgray = Gdk::Color ("gray50");
  dgray = Gdk::Color ("gray25");
  rgbcolormap->alloc_color(black);
  rgbcolormap->alloc_color(white);
  rgbcolormap->alloc_color(red);
  rgbcolormap->alloc_color(green);
  rgbcolormap->alloc_color(blue);
  rgbcolormap->alloc_color(lgray);
  rgbcolormap->alloc_color(mgray);
  rgbcolormap->alloc_color(dgray);

}

bool HistogramRGBArea::on_expose_event(GdkEventExpose* event) {

  Glib::RefPtr<Gdk::Window> window = get_window();

  // on_realize & RenderRGBMarks have to be called before 
  if (rgbgc_ && overlay) {
      window->draw_drawable (rgbgc_, overlay, 0, 0, 0, 0, -1, -1);
  }

  return true;
}

bool HistogramRGBArea::on_button_press_event (GdkEventButton* event) {

  if (event->type==GDK_2BUTTON_PRESS && event->button==1) {
  // do something. Maybe un-freeze ? 
  }
  return true;
}

void HistogramRGBArea::styleChanged (const Glib::RefPtr<Gtk::Style>& style) {

  white = get_style()->get_base(Gtk::STATE_NORMAL);
  queue_draw ();
}

//
//
//
// HistogramArea
HistogramArea::HistogramArea () : 
      valid(false), showFull(true), oldwidth(-1), needLuma(true), needRed(true), needGreen(true), needBlue(true), rawMode(false) {

         lhist(256);
	 rhist(256);
	 ghist(256);
	 bhist(256);

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

void HistogramArea::updateOptions (bool r, bool g, bool b, bool l, bool raw) {

    needRed   = r;
    needGreen = g;
    needBlue  = b;
    needLuma  = l;
    rawMode   = raw;

    renderHistogram ();
}

int histupdateUI (void* data) {

    HistogramAreaIdleHelper* haih = (HistogramAreaIdleHelper*)data;

    if (haih->destroyed) {
        if (haih->pending == 1)
            delete haih;
        else    
            haih->pending--;

        return 0;
    }
    
    haih->harea->renderHistogram ();
    haih->harea->queue_draw ();

    haih->pending--;

    return 0;
}

void HistogramArea::update (LUTu &histRed, LUTu &histGreen, LUTu &histBlue, LUTu &histLuma, LUTu &histRedRaw, LUTu &histGreenRaw, LUTu &histBlueRaw) {
	
    if (histRed) {
        lhist=histLuma;
        rhist=histRed; ghist=histGreen; bhist=histBlue;
        rhistRaw=histRedRaw; ghistRaw=histGreenRaw; bhistRaw=histBlueRaw;

        valid = true;
    }
    else
        valid = false;
	
    haih->pending++;
    g_idle_add (histupdateUI, haih);
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
    // For RAW mode use the other hists
    LUTu& rh = rawMode ? rhistRaw : rhist;
    LUTu& gh = rawMode ? ghistRaw : ghist;
    LUTu& bh = rawMode ? bhistRaw : bhist;

    // compute height of the full histogram (realheight) and
    // does not take into account 0 and 255 values
    // them are handled separately

    int fullhistheight = 0;
    for (int i=1; i<255; i++) {
        if (needLuma && lhist[i]>fullhistheight)
	        fullhistheight = lhist[i];
        if (needRed && (rawMode?rhistRaw:rhist)[i]>fullhistheight)
	        fullhistheight = rh[i];
        if (needGreen && (rawMode?ghistRaw:ghist)[i]>fullhistheight)
	        fullhistheight = gh[i];
        if (needBlue && (rawMode?bhistRaw:bhist)[i]>fullhistheight)
	        fullhistheight = bh[i];
    }

    int realhistheight = fullhistheight;

    if (!showFull) {    
        int area1thres = 0;
        int area2thres = 0;
        int area = 0;
        for (int i=0; i<fullhistheight; i++) {
            for (int j=0; j<256; j++)
                if ((needLuma && !rawMode && lhist[j]>i) || (needRed && rh[j]>i) || (needGreen && gh[j]>i) || (needBlue && bh[j]>i))
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

    int ui = 0, oi = 0;

    if (needLuma && !rawMode) {
        drawCurve(cr, lhist, realhistheight, winw, winh);
        cr->set_source_rgb (0.75, 0.75, 0.75);
        cr->fill_preserve ();
        cr->set_source_rgb (0.5, 0.5, 0.5);
    	cr->stroke ();

        drawMarks(cr, lhist, realhistheight, winw, ui, oi);
    }
    if (needRed) {
        drawCurve(cr, rh, realhistheight, winw, winh);
        cr->set_source_rgb (1.0, 0.0, 0.0);
    	cr->stroke ();

        drawMarks(cr, rh, realhistheight, winw, ui, oi);
    }
    if (needGreen) {
        drawCurve(cr, gh, realhistheight, winw, winh);
        cr->set_source_rgb (0.0, 1.0, 0.0);
    	cr->stroke ();

        drawMarks(cr, gh, realhistheight, winw, ui, oi);
    }
    if (needBlue) {
        drawCurve(cr, bh, realhistheight, winw, winh);
        cr->set_source_rgb (0.0, 0.0, 1.0);
    	cr->stroke ();

        drawMarks(cr, bh, realhistheight, winw, ui, oi);
    }
  }

  bgc->set_foreground (mgray);
  backBuffer->draw_rectangle (bgc, false, 0, 0, winw-1, winh-1);

  bgc->set_line_attributes (1, Gdk::LINE_ON_OFF_DASH, Gdk::CAP_NOT_LAST, Gdk::JOIN_MITER);

  backBuffer->draw_line (bgc, winw/4, 0, winw/4, winh);
  backBuffer->draw_line (bgc, 2*winw/4, 0, 2*winw/4, winh);
  backBuffer->draw_line (bgc, 3*winw/4, 0, 3*winw/4, winh);
  backBuffer->draw_line (bgc, 0, winh/4, winw, winh/4);
  backBuffer->draw_line (bgc, 0, 2*winh/4, winw, 2*winh/4);
  backBuffer->draw_line (bgc, 0, 3*winh/4, winw, 3*winh/4);

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

void HistogramArea::drawCurve(Cairo::RefPtr<Cairo::Context> &cr,
    LUTu & data, double scale, int hsize, int vsize)
{
    cr->move_to (0, vsize-1);
    for (int i = 0; i < 256; i++) {
        double val = data[i] * (double)(vsize-2) / scale;
        if (val > vsize - 1)
            val = vsize - 1;
        cr->line_to ((i/255.0)*(hsize - 1), vsize - 1 - val);
	}
    cr->line_to (hsize - 1, vsize - 1);
}

void HistogramArea::drawMarks(Cairo::RefPtr<Cairo::Context> &cr,
    LUTu & data, double scale, int hsize, int & ui, int & oi)
{
    int s = 8;
    
    if(data[0] > scale) {
            cr->rectangle(0, (ui++)*s, s, s);
    }
    if(data[255] > scale) {
            cr->rectangle(hsize - s, (oi++)*s, s, s);
    }
    cr->fill();
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
    return true;
}
