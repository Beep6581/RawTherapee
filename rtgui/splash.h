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
#include "rtsurface.h"

class SplashImage final :
    public Gtk::DrawingArea
{

private:
    std::shared_ptr<RTSurface> surface;
    Glib::RefPtr<Pango::Layout> version;

public:
    SplashImage ();
    bool on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr) override;
    Gtk::SizeRequestMode get_request_mode_vfunc () const override;
    void get_preferred_height_vfunc (int &minimum_height, int &natural_height) const override;
    void get_preferred_width_vfunc (int &minimum_width, int &natural_width) const override;
    void get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const override;
    void get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const override;
};

//class Splash : public Gtk::Window {
class Splash final : public Gtk::Dialog
{

private:
    SplashImage* splashImage;
    Gtk::Notebook* nb;
    Gtk::ScrolledWindow* releaseNotesSW;

public:
    explicit Splash (Gtk::Window& parent);

    bool hasReleaseNotes()
    {
        return releaseNotesSW != nullptr;
    };
    void showReleaseNotes();
    bool on_timer ();
    //virtual bool on_button_release_event (GdkEventButton* event);
    void closePressed();
};
