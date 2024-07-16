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
#include "rtappchooserdialog.h"

#if !(defined _WIN32 || defined __APPLE__)
#define GTKMM_APPCHOOSERDIALOG
#endif

RTAppChooserDialog::~RTAppChooserDialog() {}

#ifdef GTKMM_APPCHOOSERDIALOG // Use Gtk::AppChooserDialog directly.

RTAppChooserDialog::RTAppChooserDialog(const Glib::ustring &content_type) :
    Gtk::AppChooserDialog(content_type)
{
}

Glib::RefPtr<Gio::AppInfo> RTAppChooserDialog::get_app_info()
{
    return Gtk::AppChooserDialog::get_app_info();
}

Glib::RefPtr<const Gio::AppInfo> RTAppChooserDialog::get_app_info() const
{
    return Gtk::AppChooserDialog::get_app_info();
}

#else // Work around bugs with GLib and glibmm.

RTAppChooserDialog::RTAppChooserDialog(const Glib::ustring &content_type) :
    Gtk::AppChooserDialog(content_type)
{
    // GTK calls a faulty GLib function to update the most recently selected
    // application after an application is selected. This removes all signal
    // handlers to prevent the function call.
    auto signal_id = g_signal_lookup("response", GTK_TYPE_APP_CHOOSER_DIALOG);
    while (true) {
        auto handler_id = g_signal_handler_find(gobj(), G_SIGNAL_MATCH_ID, signal_id, GQuark(), nullptr, nullptr, nullptr);
        if (!handler_id) {
            break;
        }
        g_signal_handler_disconnect(gobj(), handler_id);
    }
}

Glib::RefPtr<Gio::AppInfo> RTAppChooserDialog::get_app_info()
{
    // glibmm wrapping of GAppInfo does not work on some platforms. Manually
    // wrap it here.
    GAppInfo *gAppInfo = gtk_app_chooser_get_app_info(GTK_APP_CHOOSER(gobj()));
    return Glib::wrap(gAppInfo, true);
}

Glib::RefPtr<const Gio::AppInfo> RTAppChooserDialog::get_app_info() const
{
    GAppInfo *gAppInfo = gtk_app_chooser_get_app_info(GTK_APP_CHOOSER(
        const_cast<GtkAppChooserDialog *>(gobj())));
    return Glib::wrap(gAppInfo, true);
}

#endif
