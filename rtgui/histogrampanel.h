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

#include <vector>

#include <gtkmm.h>

#include <glibmm/ustring.h>

#include <cairomm/cairomm.h>

#include "delayed.h"
#include "guiutils.h"
#include "options.h"
#include "pointermotionlistener.h"

#include "rtengine/array2D.h"
#include "rtengine/LUT.h"
#include "rtengine/noncopyable.h"


class HistogramScaling
{
public:
    double factor;
    HistogramScaling() : factor(10.0) {}
    double log (double vsize, double val);
};

class HistogramRGBArea : public Gtk::DrawingArea, protected HistogramScaling, public rtengine::NonCopyable
{
protected:
    // Saved pointer parameters
    int r;
    int g;
    int b;
    float lab_L;
    float lab_a;
    float lab_b;
    bool pointerValid;

    // Drawing options
    bool needRed;
    bool needGreen;
    bool needBlue;
    bool needLuma;
    bool needChroma;
    Options::ScopeType scopeType;
    int scaleMode;
    bool showBar;

    const double padding = 5.0;

    // Internal drawing functions
    void updateDrawingArea (const ::Cairo::RefPtr< Cairo::Context> &cc);
    virtual void drawBar(const Cairo::RefPtr<Cairo::Context> &cc, const double value, const double max_value, const int winw, const int winh) = 0;

    // GtkDrawingArea override functions
    void on_realize() override;
    bool on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr) override;

    // Widget size management functions
    void getPreferredThickness(int& min_thickness, int& natural_length) const;
    void getPreferredLength(int& min_length, int& natural_length) const;
    void getPreferredThicknessForLength(int length, int& min_thickness, int& natural_length) const;
    void getPreferredLengthForThickness(int thickness, int& min_length, int& natural_length) const;

public:
    HistogramRGBArea();
    ~HistogramRGBArea() override;

    // Update pointer values: returns true if widget needs redrawing
    bool updatePointer (const int new_r, const int new_g, const int new_b, const rtengine::procparams::ColorManagementParams *cmp = nullptr);
    // Update internal parameters from options
    void updateFromOptions ();

    // Event management functions
    bool on_button_press_event (GdkEventButton* event) override;
    void factorChanged (double newFactor);
};

class HistogramRGBAreaHori final : public HistogramRGBArea
{
private:
    // Internal drawing function
    void drawBar(const Cairo::RefPtr<Cairo::Context> &cc, const double value, const double max_value, const int winw, const int winh) override;

    // Widget size management functions
    Gtk::SizeRequestMode get_request_mode_vfunc () const override;
    void get_preferred_height_vfunc (int& minimum_height, int& natural_height) const override;
    void get_preferred_width_vfunc (int &minimum_width, int &natural_width) const override;
    void get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const override;
    void get_preferred_width_for_height_vfunc (int h, int &minimum_width, int &natural_width) const override;
};

class HistogramRGBAreaVert final : public HistogramRGBArea
{
private:
    // Internal drawing function
    void drawBar(const Cairo::RefPtr<Cairo::Context> &cc, const double value, const double max_value, const int winw, const int winh) override;

    // Widget size management functions
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

class HistogramArea final : public Gtk::DrawingArea, private HistogramScaling, public rtengine::NonCopyable
{
public:
    typedef sigc::signal<void, double> type_signal_factor_changed;
    typedef sigc::signal<void, float> SignalBrightnessChanged;

    static constexpr float MIN_BRIGHT = 0.1;
    static constexpr float MAX_BRIGHT = 3;
private:
    IdleRegister idle_register;
    type_signal_factor_changed sigFactorChanged;

protected:
    // Histogram parameters
    LUTu rhist, ghist, bhist, lhist, chist;
    LUTu rhistRaw, ghistRaw, bhistRaw;
    int vectorscope_scale;
    array2D<int> vect_hc, vect_hs;
    std::vector<unsigned char> vect_hc_buffer, vect_hs_buffer;
    bool vect_hc_buffer_dirty, vect_hs_buffer_dirty;
    int waveform_scale;
    array2D<int> rwave, gwave, bwave, lwave;
    MyRWMutex wave_mutex;
    std::vector<unsigned char> parade_buffer_r;
    std::vector<unsigned char> parade_buffer_g;
    std::vector<unsigned char> parade_buffer_b;
    bool parade_buffer_r_dirty, parade_buffer_g_dirty, parade_buffer_b_dirty;
    std::vector<unsigned char> wave_buffer;
    std::vector<unsigned char> wave_buffer_luma;
    bool wave_buffer_dirty, wave_buffer_luma_dirty;
    bool LUT_valid;

    // Intensity of waveform and vectorscope trace
    float trace_brightness;
    SignalBrightnessChanged signal_brightness_changed;

    // Saved pointer parameters
    int pointer_red;
    int pointer_green;
    int pointer_blue;
    float pointer_a;
    float pointer_b;
    bool pointer_valid;

    // Drawing options
    bool needRed;
    bool needGreen;
    bool needBlue;
    bool needLuma;
    bool needChroma;
    bool needPointer;
    Options::ScopeType scopeType;
    int drawMode;
    DrawModeListener *myDrawModeListener;

    // Motion event management
    bool isPressed;
    double movingPosition;

    double padding = 5.0;

public:
    explicit HistogramArea(DrawModeListener *fml = nullptr);
    ~HistogramArea() override;

    // Update pointer values: returns true if widget needs redrawing
    bool updatePointer(const int r, const int g, const int b, const rtengine::procparams::ColorManagementParams *cmp = nullptr);
    // Update histogram data
    void update(
        const LUTu& histRed,
        const LUTu& histGreen,
        const LUTu& histBlue,
        const LUTu& histLuma,
        const LUTu& histChroma,
        const LUTu& histRedRaw,
        const LUTu& histGreenRaw,
        const LUTu& histBlueRaw,
        int vectorscopeScale,
        const array2D<int>& vectorscopeHC,
        const array2D<int>& vectorscopeHS,
        int waveformScale,
        const array2D<int>& waveformRed,
        const array2D<int>& waveformGreen,
        const array2D<int>& waveformBlue,
        const array2D<int>& waveformLuma
    );
    // Update internal parameters from options
    void updateFromOptions();

    // Event management functions
    bool on_button_press_event (GdkEventButton* event) override;
    bool on_button_release_event (GdkEventButton* event) override;
    bool on_motion_notify_event (GdkEventMotion* event) override;

    // Brightness management functions (nominal = 1)
    float getBrightness(void);
    void setBrightness(float brightness);
    SignalBrightnessChanged getBrighnessChangedSignal(void);

    // Factor management function
    type_signal_factor_changed signal_factor_changed();

private:
    // GtkDrawingArea override functions
    void on_realize() override;
    bool on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr) override;

    // Internal drawing functions
    void updateDrawingArea (const ::Cairo::RefPtr< Cairo::Context> &cr);
    void drawCurve(const Cairo::RefPtr<Cairo::Context> &cr, const LUTu & data, const double scale, const int hsize, const int vsize);
    void drawMarks(const Cairo::RefPtr<Cairo::Context> &cr, const LUTu & data, const double scale, const int wsize, int & ui, int & oi);
    void drawParade(const Cairo::RefPtr<Cairo::Context> &cr, const int hsize, const int vsize);
    void drawVectorscope(const Cairo::RefPtr<Cairo::Context> &cr, const int hsize, const int vsize);
    void drawWaveform(const Cairo::RefPtr<Cairo::Context> &cr, const int hsize, const int vsize);

    // Widget size management functions
    Gtk::SizeRequestMode get_request_mode_vfunc () const override;
    void get_preferred_height_vfunc (int& minimum_height, int& natural_height) const override;
    void get_preferred_width_vfunc (int &minimum_width, int &natural_width) const override;
    void get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const override;
    void get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const override;
};

class HistogramPanelListener
{
public:
    virtual void scopeTypeChanged(Options::ScopeType new_type) = 0;
};

class HistogramPanel final : public Gtk::Grid, public PointerMotionListener, public DrawModeListener, public rtengine::NonCopyable
{
private:
    DelayedCall<bool, const rtengine::procparams::ColorManagementParams *, int, int, int> pointer_moved_delayed_call;

protected:

    Gtk::Grid* gfxGrid;
    Gtk::Grid* buttonGrid;
    Gtk::Box* persistentButtons;
    Gtk::Box* optionButtons;
    HistogramArea* histogramArea;
    HistogramRGBArea* histogramRGBArea;
    std::unique_ptr<HistogramRGBAreaHori> histogramRGBAreaHori;
    std::unique_ptr<HistogramRGBAreaVert> histogramRGBAreaVert;
    Gtk::ToggleButton* showRed;
    Gtk::ToggleButton* showGreen;
    Gtk::ToggleButton* showBlue;
    Gtk::ToggleButton* showValue;
    Gtk::ToggleButton* showBAR;
    Gtk::ToggleButton* showChro;
    Gtk::Button* showMode;
    Gtk::ToggleButton* scopeOptions;
    Gtk::Scale* brightnessWidget;

    Gtk::RadioButton* scopeHistBtn;
    Gtk::RadioButton* scopeHistRawBtn;
    Gtk::RadioButton* scopeParadeBtn;
    Gtk::RadioButton* scopeWaveBtn;
    Gtk::RadioButton* scopeVectHcBtn;
    Gtk::RadioButton* scopeVectHsBtn;

    Gtk::Image *redImage;
    Gtk::Image *greenImage;
    Gtk::Image *blueImage;
    Gtk::Image *valueImage;
    Gtk::Image *barImage;
    Gtk::Image *chroImage;

    Gtk::Image *redImage_g;
    Gtk::Image *greenImage_g;
    Gtk::Image *blueImage_g;
    Gtk::Image *valueImage_g;
    Gtk::Image *barImage_g;
    Gtk::Image *chroImage_g;

    Gtk::Image *mode0Image;
    Gtk::Image *mode1Image;
    Gtk::Image *mode2Image;

    HistogramPanelListener* panel_listener;

    sigc::connection brightness_changed_connection;
    sigc::connection rconn;

    void setHistInvalid ();
    void showRGBBar();

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
        const LUTu& histBlueRaw,
        int vectorscopeScale,
        const array2D<int>& vectorscopeHC,
        const array2D<int>& vectorscopeHS,
        int waveformScale,
        const array2D<int>& waveformRed,
        const array2D<int>& waveformGreen,
        const array2D<int>& waveformBlue,
        const array2D<int>& waveformLuma
    )
    {
        histogramArea->update(histRed, histGreen, histBlue, histLuma, histChroma, histRedRaw, histGreenRaw, histBlueRaw, vectorscopeScale, vectorscopeHC, vectorscopeHS, waveformScale, waveformRed, waveformGreen, waveformBlue, waveformLuma);
    }
    // pointermotionlistener interface
    void pointerMoved (bool validPos, const rtengine::procparams::ColorManagementParams &cmp, int x, int y, int r, int g, int b, bool isRaw = false) override;

    void reorder (Gtk::PositionType position);
    void red_toggled ();
    void green_toggled ();
    void blue_toggled ();
    void value_toggled ();
    void chro_toggled ();
    void bar_toggled ();
    void mode_released ();
    void brightnessWidgetValueChanged();
    void brightnessUpdated(float brightness);
    void scopeOptionsToggled();
    void type_selected(Gtk::RadioButton* button);
    void type_changed ();
    void rgbv_toggled ();
    void resized (Gtk::Allocation& req);

    // drawModeListener interface
    void toggleButtonMode () override;

    void setPanelListener(HistogramPanelListener* listener);
};
