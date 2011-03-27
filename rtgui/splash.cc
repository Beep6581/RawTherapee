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
#include <glib/gstdio.h>
#include <safegtk.h>

#ifndef WIN32
#include <config.h>
#endif

extern Glib::ustring argv0;
extern Glib::ustring versionString;

SplashImage::SplashImage () {

    pixbuf = safe_create_from_file (argv0+"/images/splash.png");
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
    Glib::RefPtr<Pango::Context> context = get_pango_context ();
    context->set_cairo_font_options (cfo);
    Pango::FontDescription fontd = context->get_font_description ();
    fontd.set_weight (Pango::WEIGHT_LIGHT);
    fontd.set_absolute_size (12*Pango::SCALE);
    context->set_font_description (fontd);
    Gdk::Color *textColor = new Gdk::Color();
    textColor->set_rgb(0, 0, 0);
    gc_->set_foreground(*textColor);

    int w, h;
    version = create_pango_layout (versionString);
    version->get_pixel_size (w, h);  
    window->draw_layout(gc_, pixbuf->get_width() - w - 4, pixbuf->get_height() - h - 4, version);

    return true;
}

Splash::Splash () {

    set_title (M("GENERAL_ABOUT"));
    set_border_width (4);

    Gtk::Notebook* nb = Gtk::manage (new Gtk::Notebook ());
    get_vbox()->pack_start (*nb);

    // Tab 1: the image
    splashImage = new SplashImage ();
    nb->append_page (*splashImage,  M("ABOUT_TAB_SPLASH"));
    splashImage->show ();

    // Tab 2: the informations about the current version
#if defined _WIN32 || defined __APPLE__
 	std::string buildFileName = Glib::build_filename (argv0, "AboutThisBuild.txt");
#else
	std::string buildFileName = Glib::build_filename (CREDITS_SEARCH_PATH, "AboutThisBuild.txt");
#endif
	if ( Glib::file_test(buildFileName, (Glib::FILE_TEST_EXISTS)) ) {
	    FILE *f = g_fopen (buildFileName.c_str(), "rt");
	    if (f != NULL) {
	        char* buffer = new char[1024];
	        std::ostringstream ostr;
	        while (fgets (buffer, 1024, f))
	            ostr << buffer;
	        delete [] buffer;
	        fclose (f);

	        Glib::RefPtr<Gtk::TextBuffer> textBuffer = Gtk::TextBuffer::create();
	        textBuffer->set_text((Glib::ustring)(ostr.str()));

	        Gtk::ScrolledWindow *buildSW = Gtk::manage (new Gtk::ScrolledWindow());
			Gtk::TextView *buildTV = Gtk::manage (new Gtk::TextView (textBuffer));
			buildTV->set_editable(false);
			buildSW->add(*buildTV);
		    nb->append_page (*buildSW, M("ABOUT_TAB_BUILD"));
	    }
	}

    // Tab 3: the credits
#if defined _WIN32 || defined __APPLE__
	std::string creditsFileName = Glib::build_filename (argv0, "AUTHORS.txt");
#else
	std::string creditsFileName = Glib::build_filename (CREDITS_SEARCH_PATH, "AUTHORS.txt");
#endif
	if ( Glib::file_test(creditsFileName, (Glib::FILE_TEST_EXISTS)) ) {
	    FILE *f = g_fopen (creditsFileName.c_str(), "rt");
	    if (f != NULL) {
	        char* buffer = new char[1024];
	        std::ostringstream ostr;
	        while (fgets (buffer, 1024, f))
	            ostr << buffer;
	        delete [] buffer;
	        fclose (f);

	        Glib::RefPtr<Gtk::TextBuffer> textBuffer = Gtk::TextBuffer::create();
	        textBuffer->set_text((Glib::ustring)(ostr.str()));

	        Gtk::ScrolledWindow *creditsSW = Gtk::manage (new Gtk::ScrolledWindow());
			Gtk::TextView *creditsTV = Gtk::manage (new Gtk::TextView (textBuffer));
			creditsTV->set_editable(false);
			creditsSW->add(*creditsTV);
		    nb->append_page (*creditsSW, M("ABOUT_TAB_CREDITS"));
	    }
	}

    // Tab 4: the license
#if defined _WIN32 || defined __APPLE__
	std::string licenseFileName = Glib::build_filename (argv0, "LICENSE.txt");
#else
	std::string licenseFileName = Glib::build_filename (LICENCE_SEARCH_PATH, "LICENSE.txt");
#endif
	if ( Glib::file_test(licenseFileName, (Glib::FILE_TEST_EXISTS)) ) {
	    FILE *f = g_fopen (licenseFileName.c_str(), "rt");
	    if (f != NULL) {
	        char* buffer = new char[1024];
	        std::ostringstream ostr;
	        while (fgets (buffer, 1024, f))
	            ostr << buffer;
	        delete [] buffer;
	        fclose (f);

	        Glib::RefPtr<Gtk::TextBuffer> textBuffer = Gtk::TextBuffer::create();
	        textBuffer->set_text((Glib::ustring)(ostr.str()));

	        Gtk::ScrolledWindow *licenseSW = Gtk::manage (new Gtk::ScrolledWindow());
			Gtk::TextView *creditsTV = Gtk::manage (new Gtk::TextView (textBuffer));
			creditsTV->set_editable(false);
			licenseSW->add(*creditsTV);
		    nb->append_page (*licenseSW, M("ABOUT_TAB_LICENSE"));
	    }
	}


    set_position (Gtk::WIN_POS_CENTER);
    //add_events(Gdk::BUTTON_RELEASE_MASK);
    set_resizable (true);

    nb->set_current_page (0);

    show_all_children ();
    set_modal (true);
    set_keep_above (true);
}

Splash::Splash (int maxtime) {

    set_title (M("GENERAL_ABOUT"));

    splashImage = new SplashImage ();
//    add (*splashImage);
    get_vbox()->pack_start (*splashImage);
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
