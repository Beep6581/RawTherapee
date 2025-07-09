/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2021 Lawrence Lee <billee@ucdavis.edu>
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

#include <giomm/appinfo.h>
#include <glibmm/refptr.h>
#include <glibmm/ustring.h>
#include <gtkmm/appchooserdialog.h>

/**
 * Custom version of gtkmm's Gtk::AppChooserDialog to work around crashes
 * (https://gitlab.gnome.org/GNOME/glib/-/issues/1104 and
 * https://gitlab.gnome.org/GNOME/glibmm/-/issues/94).
 */
class RTAppChooserDialog : public Gtk::AppChooserDialog
{
public:
    RTAppChooserDialog(const Glib::ustring &content_type);
    ~RTAppChooserDialog();

    Glib::RefPtr<Gio::AppInfo> get_app_info();
    Glib::RefPtr<const Gio::AppInfo> get_app_info() const;
};
