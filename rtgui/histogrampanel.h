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
#ifndef _HISTOGRAMPANEL_
#define _HISTOGRAMPANEL_

#include <gtkmm.h>
#include <glibmm.h>
#include <LUT.h>

class HistogramArea;
struct HistogramAreaIdleHelper {
    HistogramArea* harea;
    bool destroyed;
    int pending;
};

class HistogramArea : public Gtk::DrawingArea {  

  protected:

    Glib::RefPtr<Gdk::GC> gc_;
    Glib::RefPtr<Gdk::Pixmap> backBuffer;

    Gdk::Color black;
    Gdk::Color white;
    Gdk::Color red;
    Gdk::Color green;
    Gdk::Color blue;
    Gdk::Color lgray;
    Gdk::Color mgray;
    Gdk::Color dgray;
    LUTu lhist, rhist, ghist, bhist;
    LUTu lhistRaw, rhistRaw, ghistRaw, bhistRaw;

    bool valid;
    bool showFull;
    int oldwidth, oldheight;

    bool needLuma, needRed, needGreen, needBlue, rawMode;

    HistogramAreaIdleHelper* haih;

  public:
         
    HistogramArea();
    ~HistogramArea();

    void renderHistogram ();
    void update (LUTu &histRed, LUTu &histGreen, LUTu &histBlue, LUTu &histLuma, LUTu &histRedRaw, LUTu &histGreenRaw, LUTu &histBlueRaw);
    void updateOptions (bool r, bool g, bool b, bool l, bool raw);
    void on_realize();
    bool on_expose_event(GdkEventExpose* event);
    bool on_button_press_event (GdkEventButton* event);
    void styleChanged (const Glib::RefPtr<Gtk::Style>& style);
  private:
    void drawCurve(Cairo::RefPtr<Cairo::Context> &cr,
				   LUTu & data, double scale, int hsize, int vsize);
    void drawMarks(Cairo::RefPtr<Cairo::Context> &cr,
				   LUTu & data, double scale, int hsize, int & ui, int & oi);
};

class HistogramPanel : public Gtk::HBox {

  protected:

    HistogramArea* histogramArea;
    Gtk::ToggleButton* showRed;
    Gtk::ToggleButton* showGreen;
    Gtk::ToggleButton* showBlue;
    Gtk::ToggleButton* showValue;
    Gtk::ToggleButton* showRAW;
    
    sigc::connection rconn;
    
  public:

    HistogramPanel ();

    void histogramChanged (LUTu &histRed, LUTu &histGreen, LUTu &histBlue, LUTu &histLuma, LUTu &histRedRaw, LUTu &histGreenRaw, LUTu &histBlueRaw) { 
        histogramArea->update (histRed, histGreen, histBlue, histLuma, histRedRaw, histGreenRaw, histBlueRaw);
    }
    void rgbv_toggled ();
    void resized (Gtk::Allocation& req);
};

#endif
