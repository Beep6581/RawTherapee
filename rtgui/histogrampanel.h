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
#include <cairomm/cairomm.h>
#include "../rtengine/LUT.h"
#include "../rtengine/improccoordinator.h"
#include "guiutils.h"

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

class HistogramRGBArea : public Gtk::DrawingArea, public BackBuffer
{

    typedef const double (*TMatrix)[3];

protected:

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
    bool needChroma;

    Gtk::Grid* parent;

    HistogramRGBAreaIdleHelper* harih;

public:

    HistogramRGBArea();
    ~HistogramRGBArea();

    void updateBackBuffer (int r, int g, int b, Glib::ustring profile = "", Glib::ustring profileW = "");
    void updateFreeze (bool f);
    bool getFreeze ();
    bool getShow ();
    void setParent (Gtk::Grid* p)
    {
        parent = p;
    };

    void update (int val, int rh, int gh, int bh);
    void updateOptions (bool r, bool g, bool b, bool l, bool raw, bool show, bool c);

    void on_realize();
    bool on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr);
    bool on_button_press_event (GdkEventButton* event);
private:
    void rgb2lab (Glib::ustring profile, Glib::ustring profileW, int r, int g, int b, float &LAB_l, float &LAB_a, float &LAB_b);
    Gtk::SizeRequestMode get_request_mode_vfunc () const;
    void get_preferred_height_vfunc (int& minimum_height, int& natural_height) const;
    void get_preferred_width_vfunc (int &minimum_width, int &natural_width) const;
    void get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const;
    void get_preferred_width_for_height_vfunc (int width, int &minimum_width, int &natural_width) const;
    // Some ...
};


class FullModeListener
{
public:
    virtual ~FullModeListener() {}
    virtual void toggle_button_full () {}
};

class HistogramArea : public Gtk::DrawingArea, public BackBuffer
{

protected:

    LUTu lhist, rhist, ghist, bhist, chist;
    LUTu lhistRaw, rhistRaw, ghistRaw, bhistRaw;

    bool valid;
    bool fullMode;
    FullModeListener *myFullModeListener;
    int oldwidth, oldheight;

    bool needLuma, needRed, needGreen, needBlue, rawMode, needChroma;

    HistogramAreaIdleHelper* haih;

public:

    explicit HistogramArea(FullModeListener *fml = nullptr);
    ~HistogramArea();

    void updateBackBuffer ();
    void update (LUTu &histRed, LUTu &histGreen, LUTu &histBlue, LUTu &histLuma, LUTu &histRedRaw, LUTu &histGreenRaw, LUTu &histBlueRaw, LUTu &histChroma);
    void updateOptions (bool r, bool g, bool b, bool l, bool raw, bool full , bool c);
    void on_realize();
    bool on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr);
    bool on_button_press_event (GdkEventButton* event);

private:
    void drawCurve(Cairo::RefPtr<Cairo::Context> &cr, LUTu & data, double scale, int hsize, int vsize);
    void drawMarks(Cairo::RefPtr<Cairo::Context> &cr, LUTu & data, double scale, int hsize, int & ui, int & oi);
    Gtk::SizeRequestMode get_request_mode_vfunc () const;
    void get_preferred_height_vfunc (int& minimum_height, int& natural_height) const;
    void get_preferred_width_vfunc (int &minimum_width, int &natural_width) const;
    void get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const;
    void get_preferred_width_for_height_vfunc (int width, int &minimum_width, int &natural_width) const;
};

class HistogramPanel : public Gtk::Grid, public PointerMotionListener, public FullModeListener
{

protected:

    Gtk::Grid* gfxGrid;
    Gtk::Grid* buttonGrid;
    HistogramArea* histogramArea;
    HistogramRGBArea* histogramRGBArea;
    Gtk::ToggleButton* showRed;
    Gtk::ToggleButton* showGreen;
    Gtk::ToggleButton* showBlue;
    Gtk::ToggleButton* showValue;
    Gtk::ToggleButton* showRAW;
    Gtk::ToggleButton* showFull;
    Gtk::ToggleButton* showBAR;
    Gtk::ToggleButton* showChro;

    Gtk::Image *redImage;
    Gtk::Image *greenImage;
    Gtk::Image *blueImage;
    Gtk::Image *valueImage;
    Gtk::Image *rawImage;
    Gtk::Image *fullImage;
    Gtk::Image *barImage;
    Gtk::Image *chroImage;

    Gtk::Image *redImage_g;
    Gtk::Image *greenImage_g;
    Gtk::Image *blueImage_g;
    Gtk::Image *valueImage_g;
    Gtk::Image *rawImage_g;
    Gtk::Image *fullImage_g;
    Gtk::Image *barImage_g;
    Gtk::Image *chroImage_g;


    sigc::connection rconn;
    void setHistInvalid ();

public:

    HistogramPanel ();
    ~HistogramPanel ();

    void histogramChanged (LUTu &histRed, LUTu &histGreen, LUTu &histBlue, LUTu &histLuma, LUTu &histRedRaw, LUTu &histGreenRaw, LUTu &histBlueRaw, LUTu &histChroma)
    {
        histogramArea->update (histRed, histGreen, histBlue, histLuma, histRedRaw, histGreenRaw, histBlueRaw, histChroma);
    }
    // pointermotionlistener interface
    void pointerMoved (bool validPos, Glib::ustring profile, Glib::ustring profileW, int x, int y, int r, int g, int b);
    // added pointermotionlistener interface
    void toggleFreeze();
    // TODO should be protected
    void setHistRGBInvalid ();

    void reorder (Gtk::PositionType position);
    void red_toggled ();
    void green_toggled ();
    void blue_toggled ();
    void value_toggled ();
    void raw_toggled ();
    void full_toggled ();
    void chro_toggled ();
    void bar_toggled ();
    void rgbv_toggled ();
    void resized (Gtk::Allocation& req);

    // fullModeListener interface
    void toggle_button_full ();
};

#endif
