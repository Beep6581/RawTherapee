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

#include "rtimage.h"
#include "guiutils.h"

class EditorPanel;
struct ExternalEditor;
class RTWindow;

class EditWindow :
    public Gtk::Window
{

private:
    RTWindow* parent;
    RTImage appIcon;

    Gtk::Notebook* mainNB;
    std::set<Glib::ustring> filesEdited;
    std::map<Glib::ustring, EditorPanel*> epanels;

    sigc::signal<void> externalEditorChangedSignal;

    bool isFullscreen;
    bool isClosed;
    bool isMinimized;
    sigc::connection onConfEventConn;
    void toggleFullscreen ();

    IdleRegister idle_register;

public:
    // Check if the system has more than one display and option is set
    static bool isMultiDisplayEnabled();

    // Should only be created once
    static EditWindow* getInstance(RTWindow* p);

    explicit EditWindow (RTWindow* p);

    void writeOptions();
    void addEditorPanel (EditorPanel* ep, const std::string &name);
    void remEditorPanel (EditorPanel* ep);
    bool selectEditorPanel(const std::string &name);
    bool closeOpenEditors();
    bool isProcessing();
    void updateExternalEditorWidget(int selectedIndex, const std::vector<ExternalEditor> &editors);
    void updateToolPanelToolLocations(
        const std::vector<Glib::ustring> &favorites, bool cloneFavoriteTools);

    void toFront();
    bool keyPressed (GdkEventKey* event);
    bool on_configure_event(GdkEventConfigure* event) override;
    bool on_delete_event(GdkEventAny* event) override;
    bool on_window_state_event(GdkEventWindowState* event) override;
    void on_mainNB_switch_page(Gtk::Widget* page, guint page_num);
    void set_title_decorated(Glib::ustring fname);
    void on_realize () override;
    void get_position(int& x, int& y) const;
    void restoreWindow();
};
