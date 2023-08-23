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
#include "splash.h"

#include <glib/gstdio.h>

#include "multilangmgr.h"

extern Glib::ustring creditsPath;
extern Glib::ustring licensePath;
extern Glib::ustring versionString;

SplashImage::SplashImage () : surface(new RTSurface("splash.svg"))
{
}

bool SplashImage::on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr)
{
    if (surface->hasSurface()) {
        cr->set_source(surface->get(), 0., 0.);
        cr->rectangle(0, 0, surface->getWidth(), surface->getHeight());
        cr->fill();

        Cairo::FontOptions cfo;
        cfo.set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
        Glib::RefPtr<Pango::Context> context = get_pango_context ();
        context->set_cairo_font_options (cfo);
        Pango::FontDescription fontd = get_style_context()->get_font();
        fontd.set_weight (Pango::WEIGHT_LIGHT);
        const int fontSize = 12; // pt
        // Non-absolute size is defined in "Pango units" and shall be multiplied by
        // Pango::SCALE from "pt":
        fontd.set_size(fontSize * Pango::SCALE);
        context->set_font_description (fontd);

        int w, h;
        Glib::ustring versionStr(versionString);

        version = create_pango_layout (versionStr);
        version->set_text(versionStr);
        version->get_pixel_size (w, h);
        cr->set_source_rgb (0., 0., 0.);
        cr->set_line_width(3.);
        cr->set_line_join(Cairo::LINE_JOIN_ROUND);
        cr->move_to (surface->getWidth() - w - 32, surface->getHeight() - h - 20);
        version->add_to_cairo_context (cr);
        cr->stroke_preserve();
        cr->set_source_rgb (1., 1., 1.);
        cr->set_line_width(0.5);
        cr->stroke_preserve();
        cr->fill();
    }

    return true;
}

Gtk::SizeRequestMode SplashImage::get_request_mode_vfunc () const
{
    return Gtk::SIZE_REQUEST_CONSTANT_SIZE;
}

void SplashImage::get_preferred_height_vfunc (int &minimum_height, int &natural_height) const
{
    minimum_height = natural_height = surface ? surface->getHeight() : RTScalable::scalePixelSize(100);
}

void SplashImage::get_preferred_width_vfunc (int &minimum_width, int &natural_width) const
{
    minimum_width = natural_width = surface ? surface->getWidth() : RTScalable::scalePixelSize(100);
}

void SplashImage::get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const
{
    get_preferred_height_vfunc (minimum_height, natural_height);
}

void SplashImage::get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const
{
    get_preferred_width_vfunc (minimum_width, natural_width);
}

Splash::Splash (Gtk::Window& parent) : Gtk::Dialog(M("GENERAL_ABOUT"), parent, true)
{

    releaseNotesSW = nullptr;

    nb = Gtk::manage (new Gtk::Notebook ());
    nb->set_name ("AboutNotebook");
    get_content_area()->pack_start (*nb);

    // Add close button to bottom of the notebook
    Gtk::Button* closeButton = Gtk::manage (new Gtk::Button (M("GENERAL_CLOSE")));
    closeButton->signal_clicked().connect( sigc::mem_fun(*this, &Splash::closePressed) );
    get_action_area()->pack_start (*closeButton, Gtk::PACK_SHRINK, 0);

    Glib::RefPtr<Gtk::CssProvider> localCSS = Gtk::CssProvider::create();
    localCSS->load_from_data ("textview { font-family: monospace; }");

    // Tab 1: the image
    splashImage = Gtk::manage(new SplashImage ());
    splashImage->set_halign(Gtk::ALIGN_CENTER);
    splashImage->set_valign(Gtk::ALIGN_CENTER);
    nb->append_page (*splashImage,  M("ABOUT_TAB_SPLASH"));
    splashImage->show ();

    // Tab 2: the information about the current version
    std::string buildFileName = Glib::build_filename (creditsPath, "AboutThisBuild.txt");

    if ( Glib::file_test(buildFileName, (Glib::FILE_TEST_EXISTS)) ) {
        FILE *f = g_fopen (buildFileName.c_str (), "rt");

        if (f != nullptr) {
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

    if ( Glib::file_test(creditsFileName, (Glib::FILE_TEST_EXISTS)) ) {
        FILE *f = g_fopen (creditsFileName.c_str (), "rt");

        if (f != nullptr) {
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
    std::string licenseFileName = Glib::build_filename (licensePath, "LICENSE");

    if ( Glib::file_test(licenseFileName, (Glib::FILE_TEST_EXISTS)) ) {
        FILE *f = g_fopen (licenseFileName.c_str (), "rt");

        if (f != nullptr) {
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

    if ( Glib::file_test(releaseNotesFileName, (Glib::FILE_TEST_EXISTS)) ) {
        FILE *f = g_fopen (releaseNotesFileName.c_str (), "rt");

        if (f != nullptr) {
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
    close();
}
