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


HistogramPanel::HistogramPanel () {

   histogramArea = Gtk::manage (new HistogramArea ());
   showRed   = Gtk::manage (new Gtk::ToggleButton (M("HISTOGRAM_BUTTON_R")));
   showGreen = Gtk::manage (new Gtk::ToggleButton (M("HISTOGRAM_BUTTON_G")));
   showBlue  = Gtk::manage (new Gtk::ToggleButton (M("HISTOGRAM_BUTTON_B")));
   showValue = Gtk::manage (new Gtk::ToggleButton (M("HISTOGRAM_BUTTON_L")));
   showRAW   = Gtk::manage (new Gtk::ToggleButton (M("HISTOGRAM_BUTTON_RAW")));
   Gtk::VBox* vbox = Gtk::manage (new Gtk::VBox (false, 0));

   showRed->set_active (true);
   showGreen->set_active (true);
   showBlue->set_active (true);
   showValue->set_active (true);
   showRAW->set_active (false);
   vbox->pack_start (*showRed, Gtk::PACK_SHRINK, 2);
   vbox->pack_start (*showGreen, Gtk::PACK_SHRINK, 2);
   vbox->pack_start (*showBlue, Gtk::PACK_SHRINK, 2);
   vbox->pack_start (*showValue, Gtk::PACK_SHRINK, 2);
   vbox->pack_end (*showRAW, Gtk::PACK_SHRINK, 2);
   pack_start (*histogramArea);
   pack_end (*vbox, Gtk::PACK_SHRINK, 2);

   showRed->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::rgbv_toggled) );
   showGreen->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::rgbv_toggled) );
   showBlue->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::rgbv_toggled) );
   showValue->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::rgbv_toggled) );
   showRAW->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::rgbv_toggled) );

   show_all ();

   showRed->set_tooltip_text   (M("HISTOGRAM_TOOLTIP_R"));
   showGreen->set_tooltip_text (M("HISTOGRAM_TOOLTIP_G"));
   showBlue->set_tooltip_text  (M("HISTOGRAM_TOOLTIP_B"));
   showValue->set_tooltip_text (M("HISTOGRAM_TOOLTIP_L"));
   showRAW->set_tooltip_text   (M("HISTOGRAM_TOOLTIP_RAW"));

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

  histogramArea->updateOptions (showRed->get_active(), showGreen->get_active(), showBlue->get_active(), showValue->get_active(), showRAW->get_active());
  histogramArea->queue_draw ();
}

HistogramArea::HistogramArea () : 
      valid(false), showFull(true), oldwidth(-1), needLuma(true), needRed(true), needGreen(true), needBlue(true), rawMode(false) {

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

