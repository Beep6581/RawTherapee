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
#include <splash.h>
#include <multilangmgr.h>

extern Glib::ustring argv0;
extern Glib::ustring versionString;

SplashImage::SplashImage () {

    pixbuf = Gdk::Pixbuf::create_from_file (argv0+"/images/splash.png");
    set_size_request (pixbuf->get_width(), pixbuf->get_height());
}

void SplashImage::on_realize () {

    Gtk::DrawingArea::on_realize();
    add_events(Gdk::EXPOSURE_MASK);

    gc_ = Gdk::GC::create (get_window());
    Glib::RefPtr<Gdk::Colormap> colormap = get_default_colormap();
    Gdk::Color fontc = Gdk::Color ("white");
    colormap->alloc_color (fontc);
    gc_->set_foreground (fontc);
  
}

bool SplashImage::on_expose_event (GdkEventExpose* event) {

    Glib::RefPtr<Gdk::Window> window = get_window();
    pixbuf->render_to_drawable (window, gc_, 0, 0, 0, 0, pixbuf->get_width(), pixbuf->get_height(), Gdk::RGB_DITHER_NONE, 0, 0);

    Cairo::FontOptions cfo;
    cfo.set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
    Glib::RefPtr<Pango::Context> context = get_pango_context () ;
    context->set_cairo_font_options (cfo);
    Pango::FontDescription fontd = context->get_font_description ();
    fontd.set_weight (Pango::WEIGHT_SEMIBOLD);
    fontd.set_size (12*Pango::SCALE);
    context->set_font_description (fontd);

    version = create_pango_layout (versionString);
    int w, h;
    version->get_pixel_size (w, h);  
    window->draw_layout(gc_, pixbuf->get_width() - w - 28, 44-h, version);

    return true;
}

Splash::Splash (int maxtime) {

    set_title (M("GENERAL_ABOUT"));

    splashImage = new SplashImage ();
//    add (*splashImage);
    get_vbox()->pack_start (*splashImage);
    set_has_separator (false);
    splashImage->show ();

    if (maxtime>0)
        Glib::signal_timeout().connect (sigc::mem_fun(*this, &Splash::on_timer), maxtime);
    set_position (Gtk::WIN_POS_CENTER);
    if (maxtime>0)
        set_decorated (false);
    add_events(Gdk::BUTTON_RELEASE_MASK);
    set_resizable (false);

    set_keep_above (true);
}

bool Splash::on_timer () {

    hide ();
    return false;
}

bool Splash::on_button_release_event (GdkEventButton* event) {

    hide ();
    return true;
}
