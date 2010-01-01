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
#include <rtengine.h>

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
    unsigned int lhist[256];
    unsigned int rhist[256];
    unsigned int ghist[256];
    unsigned int bhist[256];
    bool valid;
    bool showFull;
    int oldwidth, oldheight;

    bool needVal;
    bool needRed;
    bool needGreen;
    bool needBlue;

    HistogramAreaIdleHelper* haih;

  public:
         
    HistogramArea();
    ~HistogramArea();

    void renderHistogram ();
    void update (unsigned int* rh, unsigned int* gh, unsigned int* bh, unsigned int* lh);
    void updateOptions (bool r, bool g, bool b, bool v);
    void on_realize();
    bool on_expose_event(GdkEventExpose* event);
    bool on_button_press_event (GdkEventButton* event);
    void styleChanged (const Glib::RefPtr<Gtk::Style>& style);
};

class HistogramPanel : public Gtk::HBox, public rtengine::HistogramListener {

  protected:

    HistogramArea* histogramArea;
    Gtk::ToggleButton* showRed;
    Gtk::ToggleButton* showGreen;
    Gtk::ToggleButton* showBlue;
    Gtk::ToggleButton* showValue;
    
    sigc::connection rconn;
    
  public:

    HistogramPanel ();

    void histogramChanged (unsigned int* rh, unsigned int* gh, unsigned int* bh, unsigned int* lh) { histogramArea->update (rh, gh, bh, lh); }
    void rgbv_toggled ();
    void resized (Gtk::Allocation& req);
};

#endif
