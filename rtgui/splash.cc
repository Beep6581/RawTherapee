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
#include "splash.h"
#include "multilangmgr.h"
#include <glib/gstdio.h>
#include "../rtengine/safegtk.h"

extern Glib::ustring creditsPath;
extern Glib::ustring licensePath;
extern Glib::ustring versionString;
extern Glib::ustring versionSuffixString;

SplashImage::SplashImage ()
{

    pixbuf = safe_create_from_file ("splash.png");
    set_size_request (pixbuf->get_width(), pixbuf->get_height());
}

bool SplashImage::on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr)
{

    Glib::RefPtr<Gdk::Window> window = get_window();
    Gdk::Cairo::set_source_pixbuf(cr, pixbuf);
    cr->rectangle(0, 0, pixbuf->get_width(), pixbuf->get_height());
    cr->fill();

    Cairo::FontOptions cfo;
    cfo.set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
    Glib::RefPtr<Pango::Context> context = get_pango_context ();
    context->set_cairo_font_options (cfo);
    Pango::FontDescription fontd = context->get_font_description ();
    fontd.set_weight (Pango::WEIGHT_LIGHT);
    fontd.set_absolute_size (12 * Pango::SCALE);
    context->set_font_description (fontd);

    int w, h;
    Glib::ustring versionStr(versionString);

    if (!versionSuffixString.empty()) {
        versionStr += " " + versionSuffixString;
    }

    version = create_pango_layout (versionStr);
    version->set_text(versionStr);
    version->get_pixel_size (w, h);
    cr->set_source_rgb (0., 0., 0.);
    cr->set_line_width(3.);
    cr->set_line_join(Cairo::LINE_JOIN_ROUND);
    cr->move_to (pixbuf->get_width() - w - 32, pixbuf->get_height() - h - 20);
    version->add_to_cairo_context (cr);
    cr->stroke_preserve();
    cr->set_source_rgb (1., 1., 1.);
    cr->set_line_width(0.5);
    cr->stroke_preserve();
    cr->fill();

    return true;
}

Splash::Splash (Gtk::Window& parent) : Gtk::Dialog(M("GENERAL_ABOUT"), parent, true)
{

    set_border_width (4);

    releaseNotesSW = NULL;

    nb = Gtk::manage (new Gtk::Notebook ());
    get_content_area()->pack_start (*nb);

    // Add close button to bottom of the notebook
    Gtk::Button* closeButton = Gtk::manage (new Gtk::Button (M("GENERAL_CLOSE")));
    closeButton->signal_clicked().connect( sigc::mem_fun(*this, &Splash::closePressed) );
    Gtk::HBox* bottomHBox = Gtk::manage (new Gtk::HBox ());
    bottomHBox->pack_end (*closeButton, Gtk::PACK_SHRINK, 0);
    get_content_area()->pack_start (*bottomHBox, Gtk::PACK_SHRINK, 0);

    Glib::RefPtr<Gtk::CssProvider> localCSS = Gtk::CssProvider::create();
    localCSS->load_from_data ("GtkTextView { font: monospace 8; }");

    // Tab 1: the image
    splashImage = Gtk::manage(new SplashImage ());
    nb->append_page (*splashImage,  M("ABOUT_TAB_SPLASH"));
    splashImage->show ();

    // Tab 2: the informations about the current version
    std::string buildFileName = Glib::build_filename (creditsPath, "AboutThisBuild.txt");

    if ( safe_file_test(buildFileName, (Glib::FILE_TEST_EXISTS)) ) {
        FILE *f = safe_g_fopen (buildFileName, "rt");

        if (f != NULL) {
            char* buffer = new char[1024];
            std::ostringstream ostr;

            while (fgets (buffer, 1024, f)) {
                ostr << buffer;
            }

            delete [] buffer;
            fclose (f);

            Glib::RefPtr<Gtk::TextBuffer> textBuffer = Gtk::TextBuffer::create();
            textBuffer->set_text((Glib::ustring)(ostr.str()));

            Gtk::ScrolledWindow *buildSW = Gtk::manage (new Gtk::ScrolledWindow());
            Gtk::TextView *buildTV = Gtk::manage (new Gtk::TextView (textBuffer));
            buildTV->get_style_context()->add_provider(localCSS, GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);
            buildTV->set_editable(false);
            buildTV->set_left_margin (10);
            buildTV->set_right_margin (5);
            buildSW->add(*buildTV);
            nb->append_page (*buildSW, M("ABOUT_TAB_BUILD"));
        }
    }

    // Tab 3: the credits
    std::string creditsFileName = Glib::build_filename (creditsPath, "AUTHORS.txt");

    if ( safe_file_test(creditsFileName, (Glib::FILE_TEST_EXISTS)) ) {
        FILE *f = safe_g_fopen (creditsFileName, "rt");

        if (f != NULL) {
            char* buffer = new char[1024];
            std::ostringstream ostr;

            while (fgets (buffer, 1024, f)) {
                ostr << buffer;
            }

            delete [] buffer;
            fclose (f);

            Glib::RefPtr<Gtk::TextBuffer> textBuffer = Gtk::TextBuffer::create();
            textBuffer->set_text((Glib::ustring)(ostr.str()));

            Gtk::ScrolledWindow *creditsSW = Gtk::manage (new Gtk::ScrolledWindow());
            Gtk::TextView *creditsTV = Gtk::manage (new Gtk::TextView (textBuffer));
            //creditsTV->get_style_context()->add_provider(localCSS, GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);
            creditsTV->set_left_margin (10);
            creditsTV->set_right_margin (5);
            creditsTV->set_wrap_mode(Gtk::WRAP_WORD);
            creditsTV->set_editable(false);
            creditsSW->add(*creditsTV);
            nb->append_page (*creditsSW, M("ABOUT_TAB_CREDITS"));
        }
    }

    // Tab 4: the license
    std::string licenseFileName = Glib::build_filename (licensePath, "LICENSE.txt");

    if ( safe_file_test(licenseFileName, (Glib::FILE_TEST_EXISTS)) ) {
        FILE *f = safe_g_fopen (licenseFileName, "rt");

        if (f != NULL) {
            char* buffer = new char[1024];
            std::ostringstream ostr;

            while (fgets (buffer, 1024, f)) {
                ostr << buffer;
            }

            delete [] buffer;
            fclose (f);

            Glib::RefPtr<Gtk::TextBuffer> textBuffer = Gtk::TextBuffer::create();
            textBuffer->set_text((Glib::ustring)(ostr.str()));

            Gtk::ScrolledWindow *licenseSW = Gtk::manage (new Gtk::ScrolledWindow());
            Gtk::TextView *licenseTV = Gtk::manage (new Gtk::TextView (textBuffer));
            //licenseTV->get_style_context()->add_provider(localCSS, GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);

            // set monospace font to enhance readability of formatted text
            licenseTV->set_left_margin (10);
            licenseTV->set_right_margin (5);
            licenseTV->set_editable(false);
            licenseSW->add(*licenseTV);
            nb->append_page (*licenseSW, M("ABOUT_TAB_LICENSE"));
        }
    }

    // Tab 5: the Release Notes
    std::string releaseNotesFileName = Glib::build_filename (creditsPath, "RELEASE_NOTES.txt");

    if ( safe_file_test(releaseNotesFileName, (Glib::FILE_TEST_EXISTS)) ) {
        FILE *f = safe_g_fopen (releaseNotesFileName, "rt");

        if (f != NULL) {
            char* buffer = new char[1024];
            std::ostringstream ostr;

            while (fgets (buffer, 1024, f)) {
                ostr << buffer;
            }

            delete [] buffer;
            fclose (f);

            Glib::RefPtr<Gtk::TextBuffer> textBuffer = Gtk::TextBuffer::create();
            textBuffer->set_text((Glib::ustring)(ostr.str()));

            releaseNotesSW = Gtk::manage (new Gtk::ScrolledWindow());
            Gtk::TextView *releaseNotesTV = Gtk::manage (new Gtk::TextView (textBuffer));
            releaseNotesTV->get_style_context()->add_provider(localCSS, GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);

            // set monospace font to enhance readability of formatted text
            releaseNotesTV->set_left_margin (10);
            releaseNotesTV->set_right_margin (3);
            releaseNotesTV->set_editable(false);
            releaseNotesTV->set_wrap_mode(Gtk::WRAP_WORD);
            releaseNotesSW->add(*releaseNotesTV);
            nb->append_page (*releaseNotesSW, M("ABOUT_TAB_RELEASENOTES"));
        }
    }


    set_position (Gtk::WIN_POS_CENTER);
    //add_events(Gdk::BUTTON_RELEASE_MASK);
    set_resizable (true);

    nb->set_current_page (0);

    show_all_children ();
    set_keep_above (true);
}

Splash::Splash (Gtk::Window& parent, int maxtime) : Gtk::Dialog(M("GENERAL_ABOUT"), parent, true)
{

    splashImage = Gtk::manage(new SplashImage ());
    get_content_area()->pack_start (*splashImage);
    splashImage->show ();

    if (maxtime > 0) {
        Glib::signal_timeout().connect (sigc::mem_fun(*this, &Splash::on_timer), maxtime);
    }

    set_position (Gtk::WIN_POS_CENTER);

    if (maxtime > 0) {
        set_decorated (false);
    }

    add_events(Gdk::BUTTON_RELEASE_MASK);
    set_resizable (false);

    set_keep_above (true);
}

bool Splash::on_timer ()
{

    hide ();
    return false;
}

/*
 * removed as it seem to be too sensitive in some OS
bool Splash::on_button_release_event (GdkEventButton* event) {

    hide ();
    return true;
}
*/

void Splash::showReleaseNotes()
{
    if (releaseNotesSW) {
        nb->set_current_page(nb->page_num(*releaseNotesSW));
    }
}

void Splash::closePressed()
{
    hide();
}
