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
#include "histogrampanel.h"
#include "multilangmgr.h"
#include "guiutils.h"
#include "options.h"
#include <cstring>
#include "../rtengine/LUT.h"
#include "rtimage.h"
#include "../rtengine/improccoordinator.h"
#include "../rtengine/color.h"
#include "../rtengine/opthelper.h"
#include "../rtengine/iccstore.h"
using namespace rtengine;

extern Options options;


//
//
// HistogramPanel
HistogramPanel::HistogramPanel ()
{

    set_vexpand(false);
    set_hexpand(true);
    set_valign(Gtk::ALIGN_START);
    set_halign(Gtk::ALIGN_FILL);
    set_name("HistogramPanel");

    histogramArea = Gtk::manage (new HistogramArea (this));
    histogramArea->set_hexpand(true);
    histogramArea->set_vexpand(true);
    histogramRGBArea = Gtk::manage (new HistogramRGBArea ());
    histogramRGBArea->set_hexpand(true);
    histogramRGBArea->set_vexpand(false);
    histogramRGBArea->set_halign(Gtk::ALIGN_FILL);
    histogramRGBArea->set_valign(Gtk::ALIGN_END);
    histogramRGBArea->show();

    gfxGrid = Gtk::manage (new Gtk::Grid ());
    gfxGrid->set_orientation(Gtk::ORIENTATION_VERTICAL);
    gfxGrid->set_row_spacing(1);
    gfxGrid->set_column_spacing(1);
    histogramRGBArea->setParent(gfxGrid);
    gfxGrid->add(*histogramArea);

    if (options.histogramBar) {
        gfxGrid->add (*histogramRGBArea);
    }

    redImage   = new RTImage ("histRed.png");
    greenImage = new RTImage ("histGreen.png");
    blueImage  = new RTImage ("histBlue.png");
    valueImage = new RTImage ("histValue.png");
    chroImage  = new RTImage ("histChro.png");
    rawImage   = new RTImage ("histRaw.png");
    fullImage  = new RTImage ("histFull.png");
    barImage   = new RTImage ("histBar.png");

    redImage_g   = new RTImage ("histRedg.png");
    greenImage_g = new RTImage ("histGreeng.png");
    blueImage_g  = new RTImage ("histBlueg.png");
    valueImage_g = new RTImage ("histValueg.png");
    chroImage_g  = new RTImage ("histChrog.png");
    rawImage_g   = new RTImage ("histRawg.png");
    fullImage_g  = new RTImage ("histFullg.png");
    barImage_g   = new RTImage ("histBarg.png");

    showRed   = Gtk::manage (new Gtk::ToggleButton ());
    showGreen = Gtk::manage (new Gtk::ToggleButton ());
    showBlue  = Gtk::manage (new Gtk::ToggleButton ());
    showValue = Gtk::manage (new Gtk::ToggleButton ());
    showChro  = Gtk::manage (new Gtk::ToggleButton ());
    showRAW   = Gtk::manage (new Gtk::ToggleButton ());
    showFull  = Gtk::manage (new Gtk::ToggleButton ());
    showBAR   = Gtk::manage (new Gtk::ToggleButton ());

    showRed->set_name("histButton");
    showRed->set_can_focus(false);
    showGreen->set_name("histButton");
    showGreen->set_can_focus(false);
    showBlue->set_name("histButton");
    showBlue->set_can_focus(false);
    showValue->set_name("histButton");
    showValue->set_can_focus(false);
    showChro->set_name("histButton");
    showChro->set_can_focus(false);
    showRAW->set_name("histButton");
    showRAW->set_can_focus(false);
    showFull->set_name("fullButton");
    showFull->set_can_focus(false);
    showBAR->set_name("histButton");
    showBAR->set_can_focus(false);

    showRed->set_relief (Gtk::RELIEF_NONE);
    showGreen->set_relief (Gtk::RELIEF_NONE);
    showBlue->set_relief (Gtk::RELIEF_NONE);
    showValue->set_relief (Gtk::RELIEF_NONE);
    showChro->set_relief (Gtk::RELIEF_NONE);
    showRAW->set_relief (Gtk::RELIEF_NONE);
    showFull->set_relief (Gtk::RELIEF_NONE);
    showBAR->set_relief (Gtk::RELIEF_NONE);

    showRed->set_tooltip_text   (M("HISTOGRAM_TOOLTIP_R"));
    showGreen->set_tooltip_text (M("HISTOGRAM_TOOLTIP_G"));
    showBlue->set_tooltip_text  (M("HISTOGRAM_TOOLTIP_B"));
    showValue->set_tooltip_text (M("HISTOGRAM_TOOLTIP_L"));
    showChro->set_tooltip_text  (M("HISTOGRAM_TOOLTIP_CHRO"));
    showRAW->set_tooltip_text   (M("HISTOGRAM_TOOLTIP_RAW"));
    showFull->set_tooltip_text  (M("HISTOGRAM_TOOLTIP_FULL"));
    showBAR->set_tooltip_text   (M("HISTOGRAM_TOOLTIP_BAR"));

    buttonGrid = Gtk::manage (new Gtk::Grid ());
    buttonGrid->set_orientation(Gtk::ORIENTATION_VERTICAL);
    showRed->set_active (true);
    showGreen->set_active (true);
    showBlue->set_active (true);
    showValue->set_active (true);
    showChro->set_active (false);//unactive by default

    showRAW->set_active (false);
    showFull->set_active (!options.histogramFullMode);
    showBAR->set_active (options.histogramBar);

    showRed->set_image   (showRed->get_active()   ? *redImage   : *redImage_g);
    showGreen->set_image (showGreen->get_active() ? *greenImage : *greenImage_g);
    showBlue->set_image  (showBlue->get_active()  ? *blueImage  : *blueImage_g);
    showValue->set_image (showValue->get_active() ? *valueImage : *valueImage_g);
    showChro->set_image  (showChro->get_active()   ? *chroImage : *chroImage_g);
    showRAW->set_image   (showRAW->get_active()   ? *rawImage   : *rawImage_g);
    showFull->set_image  (showFull->get_active()  ? *fullImage  : *fullImage_g);
    showBAR->set_image   (showBAR->get_active()   ? *barImage   : *barImage_g);

    showRed->set_hexpand(false);
    showRed->set_vexpand(false);
    showRed->set_halign(Gtk::ALIGN_CENTER);
    showRed->set_valign(Gtk::ALIGN_START);
    showGreen->set_hexpand(false);
    showGreen->set_vexpand(false);
    showGreen->set_halign(Gtk::ALIGN_CENTER);
    showGreen->set_valign(Gtk::ALIGN_START);
    showBlue->set_hexpand(false);
    showRed->set_vexpand(false);
    showBlue->set_halign(Gtk::ALIGN_CENTER);
    showBlue->set_valign(Gtk::ALIGN_START);
    showValue->set_hexpand(false);
    showValue->set_vexpand(false);
    showValue->set_halign(Gtk::ALIGN_CENTER);
    showValue->set_valign(Gtk::ALIGN_START);
    showChro->set_hexpand(false);
    showChro->set_vexpand(false);
    showChro->set_halign(Gtk::ALIGN_CENTER);
    showChro->set_valign(Gtk::ALIGN_START);
    showRAW->set_hexpand(false);
    showRAW->set_vexpand(false);
    showRAW->set_halign(Gtk::ALIGN_CENTER);
    showRAW->set_valign(Gtk::ALIGN_START);
    showFull->set_hexpand(false);
    showFull->set_vexpand(false);
    showFull->set_halign(Gtk::ALIGN_CENTER);
    showFull->set_valign(Gtk::ALIGN_START);
    showBAR->set_hexpand(false);
    showBAR->set_vexpand(false);
    showBAR->set_halign(Gtk::ALIGN_CENTER);
    showBAR->set_valign(Gtk::ALIGN_START);

    showRed->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::red_toggled), showRed );
    showGreen->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::green_toggled), showGreen );
    showBlue->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::blue_toggled), showBlue );
    showValue->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::value_toggled), showValue );
    showChro->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::chro_toggled), showChro );
    showRAW->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::raw_toggled), showRAW );
    showFull->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::full_toggled), showFull );
    showBAR->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::bar_toggled), showBAR );

    buttonGrid->add (*showRed);
    buttonGrid->add (*showGreen);
    buttonGrid->add (*showBlue);
    buttonGrid->add (*showValue);
    buttonGrid->add (*showChro);
    buttonGrid->add (*showRAW);
    buttonGrid->add (*showFull);
    buttonGrid->add (*showBAR);

    // Put the button vbox next to the window's border to be less disturbing
    if (options.histogramPosition == 1) {
        add (*buttonGrid);
        add (*gfxGrid);
    } else {
        add (*gfxGrid);
        add (*buttonGrid);
    }

    show_all ();

    rconn = signal_size_allocate().connect( sigc::mem_fun(*this, &HistogramPanel::resized) );
}

HistogramPanel::~HistogramPanel ()
{
    delete redImage;
    delete greenImage;
    delete blueImage;
    delete valueImage;
    delete chroImage;
    delete rawImage;
    delete fullImage;
    delete barImage;

    delete redImage_g;
    delete greenImage_g;
    delete blueImage_g;
    delete valueImage_g;
    delete chroImage_g;
    delete rawImage_g;
    delete fullImage_g;
    delete barImage_g;

}

void HistogramPanel::resized (Gtk::Allocation& req)
{

    /*
    rconn.block (true);

    int gHeight = req.get_width()/2;
    if (gHeight > 150) gHeight = 150; else if (gHeight < 100) gHeight = 100;
    int bHeight = req.get_width()/30;
    if (bHeight >  10) bHeight =  10; else if (bHeight < 5  ) bHeight = 5;
    histogramArea->set_size_request (req.get_width(), gHeight);
    histogramRGBArea->set_size_request (req.get_width(), bHeight);

    rconn.block (false);
    */

    histogramArea->updateBackBuffer ();
    histogramArea->queue_draw ();

    if (histogramRGBArea->getFreeze()) {
        histogramRGBArea->updateFreeze(false);
        // set histogramRGBArea invalid;
        histogramRGBArea->updateBackBuffer(-1, -1, -1);
        // re-set freeze to old state
        histogramRGBArea->updateFreeze(true);
        histogramRGBArea->queue_draw ();
    } else {
        // set histogramRGBArea invalid;
        histogramRGBArea->updateBackBuffer(-1, -1, -1);
        histogramRGBArea->queue_draw ();
    }
}

void HistogramPanel::red_toggled ()
{
    showRed->set_image(showRed->get_active() ? *redImage : *redImage_g);
    rgbv_toggled();
}
void HistogramPanel::green_toggled ()
{
    showGreen->set_image(showGreen->get_active() ? *greenImage : *greenImage_g);
    rgbv_toggled();
}
void HistogramPanel::blue_toggled ()
{
    showBlue->set_image(showBlue->get_active() ? *blueImage : *blueImage_g);
    rgbv_toggled();
}
void HistogramPanel::value_toggled ()
{
    removeIfThere(showValue, valueImage, false);
    removeIfThere(showValue, valueImage_g, false);
    showValue->set_image(showValue->get_active() ? *valueImage : *valueImage_g);
    rgbv_toggled();
}
void HistogramPanel::chro_toggled ()
{
    removeIfThere(showChro, chroImage, false);
    removeIfThere(showChro, chroImage_g, false);
    showChro->set_image(showChro->get_active() ? *chroImage : *chroImage_g);
    rgbv_toggled();
}

void HistogramPanel::raw_toggled ()
{
    if (showRAW->get_active()) {
        showRAW->set_image(*rawImage);
        showValue->set_sensitive(false);
        showChro->set_sensitive(false);
    } else {
        showRAW->set_image(*rawImage_g);
        showValue->set_sensitive(true);
        showChro->set_sensitive(true);
    }

    rgbv_toggled();
}
void HistogramPanel::full_toggled ()
{
    options.histogramFullMode = !showFull->get_active();
    showFull->set_image(showFull->get_active() ? *fullImage : *fullImage_g);
    rgbv_toggled();
}
void HistogramPanel::bar_toggled ()
{
    showBAR->set_image(showBAR->get_active() ? *barImage : *barImage_g);
    rgbv_toggled();
}
void HistogramPanel::rgbv_toggled ()
{
    // Update Display
    histogramArea->updateOptions (showRed->get_active(), showGreen->get_active(), showBlue->get_active(), showValue->get_active(), showRAW->get_active(), showFull->get_active(), showChro->get_active());
    histogramArea->queue_draw ();

    histogramRGBArea->updateOptions (showRed->get_active(), showGreen->get_active(), showBlue->get_active(), showValue->get_active(), showRAW->get_active(), showBAR->get_active(), showChro->get_active());
    histogramRGBArea->updateBackBuffer (0, 0, 0);
    histogramArea->queue_draw ();
}

void HistogramPanel::setHistRGBInvalid ()
{
    // do something to un-show vertical bars
    histogramRGBArea->updateBackBuffer(-1, -1, -1);
    histogramRGBArea->queue_draw ();
}

// "Freeze" is not a button, but a RMB-click, so this is not in the RGBV-Toggle method
void HistogramPanel::toggleFreeze ()
{
    if (histogramRGBArea->getFreeze()) {
        histogramRGBArea->updateFreeze(false);
    } else if (histogramRGBArea->getShow()) {
        histogramRGBArea->updateFreeze(true);
    }

    return;
}

void HistogramPanel::pointerMoved (bool validPos, Glib::ustring profile, Glib::ustring profileW, int x, int y, int r, int g, int b)
{

    if (!validPos) {
        // do something to un-show vertical bars
        histogramRGBArea->updateBackBuffer(-1, -1, -1);
        histogramRGBArea->queue_draw ();
    } else {
        // do something to show vertical bars
        histogramRGBArea->updateBackBuffer(r, g, b, profile, profileW);
        histogramRGBArea->queue_draw ();
    }
}

/*
 * Move the vertical button bar to the right side
 * only allowed values for align are Gtk::POS_LEFT and Gtk::POS_RIGHT
 */
void HistogramPanel::reorder (Gtk::PositionType align)
{
    if (align == Gtk::POS_LEFT) {
        gfxGrid->reference();
        removeIfThere(this, gfxGrid, false);
        add (*gfxGrid);
        gfxGrid->unreference();
    } else {
        buttonGrid->reference();
        removeIfThere(this, buttonGrid, false);
        add (*buttonGrid);
        buttonGrid->unreference();
    }
}

// FullModeListener interface:
void HistogramPanel::toggle_button_full ()
{
    showFull->set_active (!showFull->get_active ());
    showFull->set_image(showFull->get_active() ? *fullImage : *fullImage_g);
}

//
//
//
// HistogramRGBArea
HistogramRGBArea::HistogramRGBArea () ://needChroma unactive by default
    val(0), r(0), g(0), b(0), frozen(false), valid(false), needRed(true), needGreen(true), needBlue(true), needLuma(true), rawMode(false), showMode(options.histogramBar), barDisplayed(options.histogramBar), needChroma(false), parent(nullptr)
{

    get_style_context()->add_class("drawingarea");
    set_name("HistogramRGBArea");

    harih = new HistogramRGBAreaIdleHelper;
    harih->harea = this;
    harih->destroyed = false;
    harih->pending = 0;
}

HistogramRGBArea::~HistogramRGBArea ()
{
    idle_register.destroy();

    if (harih->pending) {
        harih->destroyed = true;
    } else {
        delete harih;
    }
}


Gtk::SizeRequestMode HistogramRGBArea::get_request_mode_vfunc () const
{
    return Gtk::SIZE_REQUEST_HEIGHT_FOR_WIDTH;
}

void HistogramRGBArea::get_preferred_height_vfunc (int &minimum_height, int &natural_height) const
{
    int minimumWidth = 0;
    int naturalWidth = 0;
    get_preferred_width_vfunc(minimumWidth, naturalWidth);
    get_preferred_height_for_width_vfunc (minimumWidth, minimum_height, natural_height);
}

void HistogramRGBArea::get_preferred_width_vfunc (int &minimum_width, int &natural_width) const
{
    minimum_width = 60;
    natural_width = 200;
}

void HistogramRGBArea::get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const
{
    int bHeight = width / 30;

    if (bHeight >  10) {
        bHeight =  10;
    } else if (bHeight < 5  ) {
        bHeight = 5;
    }

    minimum_height = bHeight;
    natural_height = bHeight;
}

void HistogramRGBArea::get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const
{
    get_preferred_width_vfunc (minimum_width, natural_width);
}

bool HistogramRGBArea::getFreeze()
{
    return(frozen);
}

bool HistogramRGBArea::getShow()
{
    return(showMode);
}

void HistogramRGBArea::updateFreeze (bool f)
{
    frozen = f;
    return;
}

void HistogramRGBArea::updateBackBuffer (int r, int g, int b, const Glib::ustring &profile, const Glib::ustring &profileW)
{
    if (!get_realized () || frozen || !showMode) {
        return;
    }

    // Mostly not necessary, but should be in some case
    GThreadLock lock; // All GUI acces from idle_add callbacks or separate thread HAVE to be protected

    Glib::RefPtr<Gdk::Window> window = get_window();
    int winx, winy, winw, winh;
    window->get_geometry(winx, winy, winw, winh);

    // This will create or update the size of the BackBuffer::surface
    setDrawRectangle(Cairo::FORMAT_ARGB32, 0, 0, winw, winh, true);

    if (surface)  {
        Cairo::RefPtr<Cairo::Context> cc = Cairo::Context::create(surface);
        Glib::RefPtr<Gtk::StyleContext> style = get_style_context();

        cc->set_source_rgba (0., 0., 0., 0.);
        cc->set_operator (Cairo::OPERATOR_CLEAR);
        cc->paint ();
        cc->set_operator (Cairo::OPERATOR_OVER);

        cc->set_antialias(Cairo::ANTIALIAS_NONE);
        cc->set_line_width (1.0);

        if ( r != -1 && g != -1 && b != -1 ) {
            if (needRed) {
                // Red
                cc->set_source_rgb(1.0, 0.0, 0.0);
                cc->move_to((int)(r * (winw / 256.0)), 0);
                cc->line_to((int)(r * (winw / 256.0)), winh - 0);
                cc->stroke();
            }

            if (needGreen) {
                // Green
                cc->set_source_rgb(0.0, 1.0, 0.0);
                cc->move_to((int)(g * (winw / 256.0)), 0);
                cc->line_to((int)(g * (winw / 256.0)), winh - 0);
                cc->stroke();
            }

            if (needBlue) {
                // Blue
                cc->set_source_rgb(0.0, 0.0, 1.0);
                cc->move_to((int)(b * (winw / 256.0)), 0);
                cc->line_to((int)(b * (winw / 256.0)), winh - 0);
                cc->stroke();
            }

            if(needLuma || needChroma) {
                float Lab_L, Lab_a, Lab_b;
                rtengine::Color::rgb2lab01(profile, profileW, r / 255.f, g / 255.f, b / 255.f, Lab_L, Lab_a, Lab_b, options.rtSettings.HistogramWorking);

                if (needLuma) {
                    // Luma
                    cc->set_source_rgb(1.0, 1.0, 1.0);
                    cc->move_to((int)((Lab_L) * (winw / 100.0)), 0);
                    cc->line_to((int)((Lab_L) * (winw / 100.0)), winh - 0);
                    cc->stroke();
                }

                if (needChroma) {
                    // Chroma
                    float chromaval = sqrt(Lab_a * Lab_a + Lab_b * Lab_b) / 1.8;
                    //  float chromaval = sqrt(Lab_a*Lab_a + Lab_b*Lab_b);
                    cc->set_source_rgb(0.0, 0.0, 0.0);
                    cc->move_to((int)(chromaval * (winw / 100.0)), 0);
                    cc->line_to((int)(chromaval * (winw / 100.0)), winh - 0);
                    cc->stroke();
                }
            }
        }
    }

    setDirty(false);
}

void HistogramRGBArea::update (int valh, int rh, int  gh, int bh)
{

    if (valh) {
        val = valh;
        r = rh;
        g = gh;
        b = bh;
        valid = true;
    } else {
        valid = false;
    }

    harih->pending++;

    const auto func = [](gpointer data) -> gboolean {
        HistogramRGBAreaIdleHelper* const harih = static_cast<HistogramRGBAreaIdleHelper*>(data);

        if (harih->destroyed) {
            if (harih->pending == 1) {
                delete harih;
            } else {
                harih->pending--;
            }

            return 0;
        }

        harih->harea->updateBackBuffer(-1, -1, -1);
        harih->harea->queue_draw ();

        harih->pending--;

        return FALSE;
    };

    idle_register.add(func, harih);
}

void HistogramRGBArea::updateOptions (bool r, bool g, bool b, bool l, bool raw, bool bar, bool c)
{

    needRed   = r;
    needGreen = g;
    needBlue  = b;
    needLuma  = l;
    rawMode   = raw;
    showMode  = bar;
    needChroma = c;

    // Histogram RGB BAR button logic goes here

    if (bar && !barDisplayed) {
        // Toggled on, add (show) the widget
        parent->add(*this);
        options.histogramBar = true;
        barDisplayed = true;
    } else if (!bar && barDisplayed) {
        // Toggled off, remove (hide) the widget
        removeIfThere(parent, this, false);
        options.histogramBar = false;
        barDisplayed = false;
        // unfreeze
        updateFreeze(false);
    }

    // Disable (but don't hide it) the bar button when RAW histogram is displayed
    if (rawMode) {
        showMode = false;
    }
}

void HistogramRGBArea::on_realize ()
{

    Gtk::DrawingArea::on_realize();
    Glib::RefPtr<Gdk::Window> window = get_window();
    add_events(Gdk::BUTTON_PRESS_MASK);
}

bool HistogramRGBArea::on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr)
{

    const Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
    style->render_background(cr, 0, 0, get_width(), get_height());

    // on_realize & updateBackBuffer have to be called before
    if (surface) {
        if (isDirty()) { // not sure this could happen...
            updateBackBuffer(-1, -1, -1);
        }

        copySurface(cr, NULL);
    }

    style->render_frame (cr, 0, 0, get_width(), get_height());

    return true;
}

bool HistogramRGBArea::on_button_press_event (GdkEventButton* event)
{

    if (event->type == GDK_2BUTTON_PRESS && event->button == 1) {
        // do something. Maybe un-freeze ?
    }

    return true;
}

//
//
//
// HistogramArea
HistogramArea::HistogramArea (FullModeListener *fml) : //needChroma unactive by default
    valid(false), fullMode(options.histogramFullMode), myFullModeListener(fml), oldwidth(-1), oldheight(-1), needLuma(true), needRed(true), needGreen(true), needBlue(true), rawMode(false), needChroma(false)
{

    lhist(256);
    rhist(256);
    ghist(256);
    bhist(256);
    chist(256);

    get_style_context()->add_class("drawingarea");
    set_name("HistogramArea");

    haih = new HistogramAreaIdleHelper;
    haih->harea = this;
    haih->destroyed = false;
    haih->pending = 0;
}

HistogramArea::~HistogramArea ()
{
    idle_register.destroy();

    if (haih->pending) {
        haih->destroyed = true;
    } else {
        delete haih;
    }
}

Gtk::SizeRequestMode HistogramArea::get_request_mode_vfunc () const
{
    return Gtk::SIZE_REQUEST_HEIGHT_FOR_WIDTH;
}

void HistogramArea::get_preferred_height_vfunc (int &minimum_height, int &natural_height) const
{
    int minimumWidth = 0;
    int naturalWidth = 0;
    get_preferred_width_vfunc (minimumWidth, naturalWidth);
    get_preferred_height_for_width_vfunc (minimumWidth, minimum_height, natural_height);
}

void HistogramArea::get_preferred_width_vfunc (int &minimum_width, int &natural_width) const
{
    minimum_width = 60;
    natural_width = 200;
}

void HistogramArea::get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const
{
    int gHeight = width / 2;

    if (gHeight > 150) {
        gHeight = 150;
    } else if (gHeight < 100) {
        gHeight = 100;
    }

    minimum_height = gHeight * 0.7;
    natural_height = gHeight;
}

void HistogramArea::get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const
{
    get_preferred_width_vfunc (minimum_width, natural_width);
}

void HistogramArea::updateOptions (bool r, bool g, bool b, bool l, bool raw, bool full, bool c)
{

    needRed   = r;
    needGreen = g;
    needBlue  = b;
    needLuma  = l;
    rawMode   = raw;
    fullMode  = !full;
    needChroma = c;

    updateBackBuffer ();
}

void HistogramArea::update (LUTu &histRed, LUTu &histGreen, LUTu &histBlue, LUTu &histLuma, LUTu &histRedRaw, LUTu &histGreenRaw, LUTu &histBlueRaw, LUTu &histChroma)
{
    if (histRed) {
        lhist = histLuma;
        chist = histChroma;
        rhist = histRed;
        ghist = histGreen;
        bhist = histBlue;
        rhistRaw = histRedRaw;
        ghistRaw = histGreenRaw;
        bhistRaw = histBlueRaw;

        valid = true;
    } else {
        valid = false;
    }

    haih->pending++;
    // Can be done outside of the GUI thread
    const auto func = [](gpointer data) -> gboolean {
        HistogramAreaIdleHelper* const haih = static_cast<HistogramAreaIdleHelper*>(data);

        if (haih->destroyed) {
            if (haih->pending == 1) {
                delete haih;
            } else {
                haih->pending--;
            }

            return 0;
        }

        haih->harea->setDirty (true);
        haih->harea->updateBackBuffer ();
        haih->harea->queue_draw ();

        haih->pending--;

        return FALSE;
    };

    idle_register.add(func, haih);
}

void HistogramArea::updateBackBuffer ()
{

    if (!get_realized ()) {
        return;
    }

    Glib::RefPtr<Gdk::Window> window = get_window();
    int winx, winy, winw, winh;
    window->get_geometry(winx, winy, winw, winh);

    // This will create or update the size of the BackBuffer::surface
    setDrawRectangle(Cairo::FORMAT_ARGB32, 0, 0, winw, winh, true);

    Cairo::RefPtr<Cairo::Context> cr = Cairo::Context::create(surface);
    const Glib::RefPtr<Gtk::StyleContext> style = get_style_context();

    cr->set_source_rgba (0., 0., 0., 0.);
    cr->set_operator (Cairo::OPERATOR_CLEAR);
    cr->paint ();
    cr->set_operator (Cairo::OPERATOR_OVER);

    if (valid) {
        // For RAW mode use the other hists
        LUTu& rh = rawMode ? rhistRaw : rhist;
        LUTu& gh = rawMode ? ghistRaw : ghist;
        LUTu& bh = rawMode ? bhistRaw : bhist;

        // make double copies of LUT, one for faster access, another one to scale down the raw histos
        LUTu rhchanged(256), ghchanged(256), bhchanged(256);
        unsigned int lhisttemp[256] ALIGNED16 {0}, chisttemp[256] ALIGNED16 {0}, rhtemp[256] ALIGNED16 {0}, ghtemp[256] ALIGNED16 {0}, bhtemp[256] ALIGNED16 {0};
        const int scale = (rawMode ? 8 : 1);

        for(int i = 0; i < 256; i++) {
            if(needLuma) {
                lhisttemp[i] = lhist[i];
            }

            if(needChroma) {
                chisttemp[i] = chist[i];
            }

            if(needRed) {
                rhchanged[i] = rhtemp[i] = rh[i] / scale;
            }

            if(needGreen) {
                ghchanged[i] = ghtemp[i] = gh[i] / scale;
            }

            if(needBlue) {
                bhchanged[i] = bhtemp[i] = bh[i] / scale;
            }
        }

        // compute height of the full histogram (realheight) and
        // does not take into account 0 and 255 values
        // them are handled separately

        unsigned int fullhistheight = 0;

        for (int i = 1; i < 255; i++) {
            if (needLuma && lhisttemp[i] > fullhistheight) {
                fullhistheight = lhisttemp[i];
            }

            if (needChroma && chisttemp[i] > fullhistheight) {
                fullhistheight = chisttemp[i];
            }

            if (needRed && rhtemp[i] > fullhistheight) {
                fullhistheight = rhtemp[i];
            }

            if (needGreen && ghtemp[i] > fullhistheight) {
                fullhistheight = ghtemp[i];
            }

            if (needBlue && bhtemp[i] > fullhistheight) {
                fullhistheight = bhtemp[i];
            }
        }

        int realhistheight = fullhistheight;

        // though much faster than before, this still takes a lot of time especially for big files if rawMode is true
        if (!fullMode) {
            int area = 0;

#ifdef __SSE2__
            vint onev = _mm_set1_epi32(1);
            vint iv = (vint)ZEROV;
#endif

            for (unsigned i = 0; i < fullhistheight; i++) {
#ifdef __SSE2__
                vint areatempv = (vint)ZEROV;

                for (int j = 0; j < 256; j += 4) {
                    vmask mask1v = _mm_cmpgt_epi32(LVI(lhisttemp[j]), iv);
                    vmask mask2v = _mm_cmpgt_epi32(LVI(rhtemp[j]), iv);
                    vmask mask3v = _mm_cmpgt_epi32(LVI(ghtemp[j]), iv);
                    vmask mask4v = _mm_cmpgt_epi32(LVI(bhtemp[j]), iv);
                    mask1v = _mm_or_si128(mask1v, mask2v);
                    mask3v = _mm_or_si128(mask3v, mask4v);
                    mask2v = _mm_cmpgt_epi32(LVI(chisttemp[j]), iv);
                    mask1v = _mm_or_si128(mask1v, mask3v);
                    mask1v = _mm_or_si128(mask1v, mask2v);
                    areatempv = _mm_add_epi32(areatempv, _mm_and_si128(mask1v, onev));

                }

                areatempv = _mm_add_epi32(areatempv, (vint)_mm_movehl_ps((vfloat)areatempv, (vfloat)areatempv));
                areatempv = _mm_add_epi32(areatempv, _mm_shuffle_epi32(areatempv, 1));
                area += _mm_cvtsi128_si32(areatempv);
                iv = _mm_add_epi32(iv, onev);

#else

                for (int j = 0; j < 256; j++)
                    if (lhisttemp[j] > i || rhtemp[j] > i || ghtemp[j] > i || bhtemp[j] > i || chisttemp[j] > i) {
                        area++;
                    }

#endif

                if ((double)area / (256 * (i + 1)) < 0.3) {
                    realhistheight = i;
                    break;
                }
            }
        }

        if (realhistheight < winh - 2) {
            realhistheight = winh - 2;
        }

        cr->set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
        cr->set_line_width (1.0);

        int ui = 0, oi = 0;

        if (needLuma && !rawMode) {
            drawCurve(cr, lhist, realhistheight, w, h);
            cr->set_source_rgb (0.65, 0.65, 0.65);
            cr->fill ();

            drawMarks(cr, lhist, realhistheight, w, ui, oi);
        }

        if (needChroma && !rawMode) {
            drawCurve(cr, chist, realhistheight, w, h);
            cr->set_source_rgb (0., 0., 0.);
            cr->stroke ();

            drawMarks(cr, chist, realhistheight, w, ui, oi);
        }

        if (needRed) {
            drawCurve(cr, rhchanged, realhistheight, w, h);
            cr->set_source_rgb (1.0, 0.0, 0.0);
            cr->stroke ();

            drawMarks(cr, rhchanged, realhistheight, w, ui, oi);
        }

        if (needGreen) {
            drawCurve(cr, ghchanged, realhistheight, w, h);
            cr->set_source_rgb (0.0, 1.0, 0.0);
            cr->stroke ();

            drawMarks(cr, ghchanged, realhistheight, w, ui, oi);
        }

        if (needBlue) {
            drawCurve(cr, bhchanged, realhistheight, w, h);
            cr->set_source_rgb (0.0, 0.0, 1.0);
            cr->stroke ();

            drawMarks(cr, bhchanged, realhistheight, w, ui, oi);
        }
    }

    cr->set_source_rgba (1., 1., 1., 0.35);
    cr->set_line_width (1.0);
    cr->set_antialias(Cairo::ANTIALIAS_NONE);

    // Draw the content
    cr->set_line_join(Cairo::LINE_JOIN_MITER);
    std::valarray<double> ch_ds (1);
    ch_ds[0] = 4;
    cr->set_dash (ch_ds, 0);

    cr->move_to(w / 4 + 0.5, 1.5);
    cr->line_to(w / 4 + 0.5, h - 2);
    cr->stroke();
    cr->move_to(2 * w / 4 + 0.5, 1.5);
    cr->line_to(2 * w / 4 + 0.5, h - 2);
    cr->stroke();
    cr->move_to(3 * w / 4 + 0.5, 1.5);
    cr->line_to(3 * w / 4 + 0.5, h - 2);
    cr->stroke();
    cr->move_to(1.5, h / 4 + 0.5);
    cr->line_to(w - 2, h / 4 + 0.5);
    cr->stroke();
    cr->move_to(1.5, 2 * h / 4 + 0.5);
    cr->line_to(w - 2, 2 * h / 4 + 0.5);
    cr->stroke();
    cr->move_to(1.5, 3 * h / 4 + 0.5);
    cr->line_to(w - 2, 3 * h / 4 + 0.5);
    cr->stroke();

    cr->unset_dash();

    // Draw the frame's border
    style->render_frame(cr, 0, 0, surface->get_width(), surface->get_height());

    oldwidth = w;
    oldheight = h;

    setDirty(false);
}

void HistogramArea::on_realize ()
{

    Gtk::DrawingArea::on_realize();
    Glib::RefPtr<Gdk::Window> window = get_window();
    add_events(Gdk::BUTTON_PRESS_MASK);
}

void HistogramArea::drawCurve(Cairo::RefPtr<Cairo::Context> &cr,
                              LUTu & data, double scale, int hsize, int vsize)
{
    cr->move_to (0, vsize - 1);
    scale = scale <= 0.f ? 0.001f : scale; // avoid division by zero and negative values

    for (int i = 0; i < 256; i++) {
        double val = data[i] * (double)(vsize - 2) / scale;

        if (val > vsize - 1) {
            val = vsize - 1;
        }

        double posX = (i / 255.0) * (hsize - 1);
        double posY = vsize - 1 - val;
        cr->line_to (posX, posY);
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

bool HistogramArea::on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr)
{

    Glib::RefPtr<Gdk::Window> window = get_window();

    if (get_width() != oldwidth || get_height() != oldheight || isDirty ()) {
        updateBackBuffer ();
    }

    const Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
    style->render_background(cr, 0, 0, get_width(), get_height());
    copySurface(cr, NULL);
    style->render_frame (cr, 0, 0, get_width(), get_height());

    return true;
}

bool HistogramArea::on_button_press_event (GdkEventButton* event)
{

    if (event->type == GDK_2BUTTON_PRESS && event->button == 1) {
        fullMode = !fullMode;
        options.histogramFullMode = fullMode;

        if (myFullModeListener) {
            myFullModeListener->toggle_button_full ();
        }

        updateBackBuffer ();
        queue_draw ();
    }

    return true;
}

