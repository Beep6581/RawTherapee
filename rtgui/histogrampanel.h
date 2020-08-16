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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <gtkmm.h>

#include <glibmm/ustring.h>

#include <cairomm/cairomm.h>

#include "delayed.h"
#include "guiutils.h"
#include "pointermotionlistener.h"

#include "../rtengine/LUT.h"
#include "../rtengine/noncopyable.h"

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

class HistogramScaling
{
public:
    double factor;
    HistogramScaling() : factor(10.0) {}
    double log (double vsize, double val);
};

class HistogramRGBArea final : public Gtk::DrawingArea, public BackBuffer, private HistogramScaling, public rtengine::NonCopyable
{
private:
    typedef const double (*TMatrix)[3];

    IdleRegister idle_register;

protected:
    int val;
    int r;
    int g;
    int b;

    bool valid;

    bool needRed;
    bool needGreen;
    bool needBlue;
    bool needLuma;
    bool needChroma;
    bool rawMode;
    bool showMode;
    bool barDisplayed;

    Gtk::Grid* parent;
    
    double padding = 5.0;

    HistogramRGBAreaIdleHelper* harih;

public:
    HistogramRGBArea();
    ~HistogramRGBArea() override;

    void updateBackBuffer (int r, int g, int b, const Glib::ustring &profile = "", const Glib::ustring &profileW = "");
    bool getShow ();
    void setParent (Gtk::Grid* p)
    {
        parent = p;
    };

    void update (int val, int rh, int gh, int bh);
    void updateOptions (bool r, bool g, bool b, bool l, bool c, bool raw, bool show);

    void on_realize() override;
    bool on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr) override;
    bool on_button_press_event (GdkEventButton* event) override;
    void factorChanged (double newFactor);

private:
    Gtk::SizeRequestMode get_request_mode_vfunc () const override;
    void get_preferred_height_vfunc (int& minimum_height, int& natural_height) const override;
    void get_preferred_width_vfunc (int &minimum_width, int &natural_width) const override;
    void get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const override;
    void get_preferred_width_for_height_vfunc (int h, int &minimum_width, int &natural_width) const override;

};

class DrawModeListener
{
public:
    virtual ~DrawModeListener() = default;
    virtual void toggleButtonMode() = 0;
};

class HistogramArea final : public Gtk::DrawingArea, public BackBuffer, private HistogramScaling, public rtengine::NonCopyable
{
public:
    typedef sigc::signal<void, double> type_signal_factor_changed;

private:
    IdleRegister idle_register;
    type_signal_factor_changed sigFactorChanged;

protected:
    LUTu rhist, ghist, bhist, lhist, chist;
    LUTu rhistRaw, ghistRaw, bhistRaw, lhistRaw; //lhistRaw is unused?

    bool valid;
    int drawMode;
    DrawModeListener *myDrawModeListener;
    int oldwidth, oldheight;

    bool needRed, needGreen, needBlue, needLuma, needChroma;
    bool rawMode;
    bool isPressed;
    double movingPosition;
    
    double padding = 5.0;

    HistogramAreaIdleHelper* haih;

public:
    explicit HistogramArea(DrawModeListener *fml = nullptr);
    ~HistogramArea() override;

    void updateBackBuffer ();
    void update(
        const LUTu& histRed,
        const LUTu& histGreen,
        const LUTu& histBlue,
        const LUTu& histLuma,
        const LUTu& histChroma,
        const LUTu& histRedRaw,
        const LUTu& histGreenRaw,
        const LUTu& histBlueRaw
    );
    void updateOptions (bool r, bool g, bool b, bool l, bool c, bool raw, int mode);
    void on_realize() override;
    bool on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr) override;
    bool on_button_press_event (GdkEventButton* event) override;
    bool on_button_release_event (GdkEventButton* event) override;
    bool on_motion_notify_event (GdkEventMotion* event) override;
    type_signal_factor_changed signal_factor_changed();

private:
    void drawCurve(Cairo::RefPtr<Cairo::Context> &cr, const LUTu & data, double scale, int hsize, int vsize);
    void drawMarks(Cairo::RefPtr<Cairo::Context> &cr, const LUTu & data, double scale, int hsize, int & ui, int & oi);
    Gtk::SizeRequestMode get_request_mode_vfunc () const override;
    void get_preferred_height_vfunc (int& minimum_height, int& natural_height) const override;
    void get_preferred_width_vfunc (int &minimum_width, int &natural_width) const override;
    void get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const override;
    void get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const override;
};

class HistogramPanel final : public Gtk::Grid, public PointerMotionListener, public DrawModeListener, public rtengine::NonCopyable
{
private:
    DelayedCall<bool, Glib::ustring, Glib::ustring, int, int, int> pointer_moved_delayed_call;

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
    Gtk::ToggleButton* showBAR;
    Gtk::ToggleButton* showChro;
    Gtk::Button* showMode;

    Gtk::Image *redImage;
    Gtk::Image *greenImage;
    Gtk::Image *blueImage;
    Gtk::Image *valueImage;
    Gtk::Image *rawImage;
    Gtk::Image *barImage;
    Gtk::Image *chroImage;

    Gtk::Image *redImage_g;
    Gtk::Image *greenImage_g;
    Gtk::Image *blueImage_g;
    Gtk::Image *valueImage_g;
    Gtk::Image *rawImage_g;
    Gtk::Image *barImage_g;
    Gtk::Image *chroImage_g;

    Gtk::Image *mode0Image;
    Gtk::Image *mode1Image;
    Gtk::Image *mode2Image;

    sigc::connection rconn;
    void setHistInvalid ();

public:

    HistogramPanel ();
    ~HistogramPanel () override;

    void histogramChanged(
        const LUTu& histRed,
        const LUTu& histGreen,
        const LUTu& histBlue,
        const LUTu& histLuma,
        const LUTu& histChroma,
        const LUTu& histRedRaw,
        const LUTu& histGreenRaw,
        const LUTu& histBlueRaw)
    {
        histogramArea->update(histRed, histGreen, histBlue, histLuma, histChroma, histRedRaw, histGreenRaw, histBlueRaw);
    }
    // pointermotionlistener interface
    void pointerMoved (bool validPos, const Glib::ustring &profile, const Glib::ustring &profileW, int x, int y, int r, int g, int b, bool isRaw = false) override;

    // TODO should be protected
    void setHistRGBInvalid ();

    void reorder (Gtk::PositionType position);
    void red_toggled ();
    void green_toggled ();
    void blue_toggled ();
    void value_toggled ();
    void raw_toggled ();
    void chro_toggled ();
    void bar_toggled ();
    void mode_released ();
    void rgbv_toggled ();
    void resized (Gtk::Allocation& req);

    // drawModeListener interface
    void toggleButtonMode () override;
};
