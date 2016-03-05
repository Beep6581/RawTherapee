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
#ifndef _RTIMAGE_
#define _RTIMAGE_

#include <gtkmm/image.h>

class Options;

/**
 * @brief A derived class of Gtk::Image in order to handle theme-related icon sets.
 */
class RTImage : public Gtk::Image
{
public:
    RTImage (const Glib::ustring& fileName, const Glib::ustring& rtlFileName = Glib::ustring());

    void changeImage (const Glib::ustring& imageName);
    static void updateImages ();

    static Glib::ustring findIconAbsolutePath (const Glib::ustring& iconName);
    static void setPaths (const Options& options);

    static Glib::RefPtr<Gdk::Pixbuf> createFromFile (const Glib::ustring& fileName);
    static Cairo::RefPtr<Cairo::ImageSurface> createFromPng (const Glib::ustring& fileName);
};

#endif
