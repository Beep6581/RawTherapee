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
#include "../rtengine/LUT.h"

#include "pointermotionlistener.h"

class HistogramArea;
struct HistogramAreaIdleHelper {
    HistogramArea* harea;
    bool destroyed;
    int pending;
};

class HistogramRGBArea;
struct HistogramRGBAreaIdleHelper {
    HistogramRGBArea* harea;
    bool destroyed;
    int pending;
};

class HistogramRGBArea : public Gtk::DrawingArea {

  protected:

    Glib::RefPtr<Gdk::GC> rgbgc_;
    Glib::RefPtr<Gdk::Pixmap> overlay;

    Gdk::Color black;
    Gdk::Color white;
    Gdk::Color red;
    Gdk::Color green;
    Gdk::Color blue;
    Gdk::Color lgray;
    Gdk::Color mgray;
    Gdk::Color dgray;

    int val;
    int r;
    int g;
    int b;

    bool frozen;
    bool valid;

    bool needRed;
    bool needGreen;
    bool needBlue;
    bool needLuma;
    bool rawMode;
    bool showMode;
    bool barDisplayed;

    Gtk::VBox* parent;

    HistogramRGBAreaIdleHelper* harih;

  public:

    HistogramRGBArea();
    ~HistogramRGBArea();

    void renderRGBMarks (int r, int g, int b);
    void updateFreeze (bool f);
    bool getFreeze ();
    bool getShow ();
    void setParent (Gtk::VBox* p) { parent = p; };

    void update (int val, int rh, int gh, int bh);
    void updateOptions (bool r, bool g, bool b, bool l, bool raw, bool show);

    void on_realize();
    bool on_expose_event(GdkEventExpose* event);
    bool on_button_press_event (GdkEventButton* event);
    void styleChanged (const Glib::RefPtr<Gtk::Style>& style);
  private:
    // Some ...
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

class HistogramPanel : public Gtk::HBox, public PointerMotionListener  {

  protected:

	Gtk::VBox* gfxVBox;
	Gtk::VBox* buttonVBox;
	HistogramArea* histogramArea;
    HistogramRGBArea* histogramRGBArea;
    Gtk::ToggleButton* showRed;
    Gtk::ToggleButton* showGreen;
    Gtk::ToggleButton* showBlue;
    Gtk::ToggleButton* showValue;
    Gtk::ToggleButton* showRAW;
    Gtk::ToggleButton* showBAR;
    
    Gtk::Image *redImage;
    Gtk::Image *greenImage;
    Gtk::Image *blueImage;
    Gtk::Image *valueImage;
    Gtk::Image *rawImage;
    Gtk::Image *barImage;
    Gtk::Image *redImage_g;
    Gtk::Image *greenImage_g;
    Gtk::Image *blueImage_g;
    Gtk::Image *valueImage_g;
    Gtk::Image *rawImage_g;
    Gtk::Image *barImage_g;


    sigc::connection rconn;
    void setHistInvalid ();
    
  public:

    HistogramPanel ();
    ~HistogramPanel ();

    void histogramChanged (LUTu &histRed, LUTu &histGreen, LUTu &histBlue, LUTu &histLuma, LUTu &histRedRaw, LUTu &histGreenRaw, LUTu &histBlueRaw) { 
        histogramArea->update (histRed, histGreen, histBlue, histLuma, histRedRaw, histGreenRaw, histBlueRaw);
    }
    // pointermotionlistener interface
    void pointerMoved (bool validPos, Glib::ustring profile, int x, int y, int r, int g, int b);
    // added pointermotionlistener interface
    void toggleFreeze();
    // TODO should be protected
    void setHistRGBInvalid ();

    void reorder (Gtk::AlignmentEnum align);
    void red_toggled ();
    void green_toggled ();
    void blue_toggled ();
    void value_toggled ();
    void raw_toggled ();
    void bar_toggled ();
    void rgbv_toggled ();
    void resized (Gtk::Allocation& req);
};

#endif
