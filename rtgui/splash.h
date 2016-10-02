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
#ifndef __SPLASH__
#define __SPLASH__

#include <gtkmm.h>

class SplashImage : public Gtk::DrawingArea
{

private:
    Glib::RefPtr<Gdk::GC> gc_;
    Glib::RefPtr<Gdk::Pixbuf> pixbuf;
    Glib::RefPtr<Pango::Layout> version;

public:
    SplashImage ();
    void on_realize ();
    bool on_expose_event (GdkEventExpose* event);
};

//class Splash : public Gtk::Window {
class Splash : public Gtk::Dialog
{

private:
    SplashImage* splashImage;
    Gtk::Notebook* nb;
    Gtk::ScrolledWindow* releaseNotesSW;

public:
    Splash (Gtk::Window& parent, int maxtime);
    explicit Splash (Gtk::Window& parent);

    bool hasReleaseNotes()
    {
        return releaseNotesSW != NULL;
    };
    void showReleaseNotes();
    bool on_timer ();
    //virtual bool on_button_release_event (GdkEventButton* event);
    void closePressed();
};

#endif
