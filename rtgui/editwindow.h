/*
 *  This file is part of RawTherapee.
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

#include <set>

#include <gtkmm.h>

#include "editorpanel.h"
#include "filepanel.h"
#include "rtimage.h"

class EditWindow :
    public Gtk::Window
{

private:
    double resolution;
    RTWindow* parent;
    RTImage appIcon;

    Gtk::Notebook* mainNB;
    std::set<Glib::ustring> filesEdited;
    std::map<Glib::ustring, EditorPanel*> epanels;

    bool isFullscreen;
    bool isClosed;
    void toggleFullscreen ();
    void restoreWindow();
    bool updateResolution();
    void setAppIcon();

public:
    // Check if the system has more than one display and option is set
    static bool isMultiDisplayEnabled();

    // Should only be created once, auto-creates window on correct display
    static EditWindow* getInstance(RTWindow* p, bool restore = true);

    explicit EditWindow (RTWindow* p);

    void writeOptions();
    void addEditorPanel (EditorPanel* ep, const std::string &name);
    void remEditorPanel (EditorPanel* ep);
    bool selectEditorPanel(const std::string &name);
    bool closeOpenEditors();
    bool isProcessing();

    void toFront();
    bool keyPressed (GdkEventKey* event);
    bool on_configure_event(GdkEventConfigure* event) override;
    bool on_delete_event(GdkEventAny* event) override;
    //bool on_window_state_event(GdkEventWindowState* event);
    void on_mainNB_switch_page(Gtk::Widget* page, guint page_num);
    void set_title_decorated(Glib::ustring fname);
    void on_realize () override;
};
