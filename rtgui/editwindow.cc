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
*  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "editwindow.h"
#include "options.h"
#include "preferences.h"
#include "cursormanager.h"
#include "rtwindow.h"
#include <gtk/gtk.h>
#include "rtimage.h"
#include "threadutils.h"
#include "../rtengine/icons.h"

// Check if the system has more than one display and option is set
bool EditWindow::isMultiDisplayEnabled()
{
    return options.multiDisplayMode > 0 && Gdk::Screen::get_default()->get_n_monitors() > 1;
}

// Should only be created once, auto-creates window on correct display
EditWindow* EditWindow::getInstance(RTWindow* p, bool restore)
{
    struct EditWindowInstance
    {
        EditWindow editWnd;

        explicit EditWindowInstance(RTWindow* p) : editWnd(p)
        {
        }
    };

    static EditWindowInstance instance_(p);
    if(restore) {
        instance_.editWnd.restoreWindow();
    }
    return &instance_.editWnd;
}

EditWindow::EditWindow (RTWindow* p) : parent(p) , isFullscreen(false), isClosed(true)
{

    Glib::ustring fName = "rt-logo-tiny.png";
    Glib::ustring fullPath = rtengine::findIconAbsolutePath(fName);

    try {
        set_default_icon_from_file (fullPath);
    } catch(Glib::Exception& ex) {
        printf ("%s\n", ex.what().c_str());
    }

    set_title_decorated("");
    set_modal(false);
    set_resizable(true);
    set_default_size(options.meowWidth, options.meowHeight);

    property_destroy_with_parent().set_value(false);

    mainNB = Gtk::manage (new Gtk::Notebook ());
    mainNB->set_scrollable (true);
    mainNB->signal_switch_page().connect_notify( sigc::mem_fun(*this, &EditWindow::on_mainNB_switch_page) );

    signal_key_press_event().connect( sigc::mem_fun(*this, &EditWindow::keyPressed) );

    Gtk::VBox* mainBox = Gtk::manage (new Gtk::VBox ());
    mainBox->pack_start (*mainNB);

    add (*mainBox);

}

void EditWindow::restoreWindow() {

    if(isClosed) {
        int meowMonitor = 0;
        if(isMultiDisplayEnabled()) {
            if(options.meowMonitor >= 0) { // use display from last session if available
                meowMonitor = std::min(options.meowMonitor, Gdk::Screen::get_default()->get_n_monitors());
            } else { // Determine the other display
                const Glib::RefPtr< Gdk::Window >& wnd = parent->get_window();
                meowMonitor = parent->get_screen()->get_monitor_at_window(wnd) == 0 ? 1 : 0;
            }
        }

        Gdk::Rectangle lMonitorRect;
        get_screen()->get_monitor_geometry(meowMonitor, lMonitorRect);
        if(options.meowMaximized) {
            move(lMonitorRect.get_x(), lMonitorRect.get_y());
            maximize();
        } else {
            resize(options.meowWidth, options.meowHeight);
            if(options.meowX <= lMonitorRect.get_x() + lMonitorRect.get_width() && options.meowY <= lMonitorRect.get_y() + lMonitorRect.get_height()) {
                move(options.meowX, options.meowY);
            } else {
                move(lMonitorRect.get_x(), lMonitorRect.get_y());
            }
        }
        show_all();

        isFullscreen = options.meowFullScreen;

        if(isFullscreen) {
            fullscreen();
        }

        isClosed = false;
    }

}

void EditWindow::on_realize ()
{
    Gtk::Window::on_realize ();

    editWindowCursorManager.init (get_window());
}

bool EditWindow::on_configure_event(GdkEventConfigure* event)
{
    if (get_realized() && is_visible()) {
        if(!is_maximized()) {
            get_position(options.meowX, options.meowY);
            get_size(options.meowWidth, options.meowHeight);
        }
        options.meowMaximized = is_maximized();
    }

    return Gtk::Widget::on_configure_event(event);
}

/*  HOMBRE: Disabling this since it's maximized when opened anyway.
 *  Someday, the EditorWindow migh save it own position and state, so it'll have to be uncommented
bool EditWindow::on_window_state_event(GdkEventWindowState* event)
{
    if (event->changed_mask & GDK_WINDOW_STATE_MAXIMIZED) {
        options.windowMaximized = event->new_window_state & GDK_WINDOW_STATE_MAXIMIZED;
    }

    return Gtk::Widget::on_window_state_event(event);
}*/

void EditWindow::on_mainNB_switch_page(Gtk::Widget* widget, guint page_num)
{
    //if (page_num > 1) {
    EditorPanel *ep = static_cast<EditorPanel*>(widget);

    if (mainNB->get_n_pages() > 1 && page_num <= (filesEdited.size() - 1)) {
        set_title_decorated(ep->getFileName());
    }

    ep->setAspect();
    //}
}

void EditWindow::addEditorPanel (EditorPanel* ep, const std::string &name)
{
    ep->setParent (parent);
    ep->setParentWindow(this);

    // construct closeable tab for the image
    Gtk::HBox* hb = Gtk::manage (new Gtk::HBox ());
    hb->pack_start (*Gtk::manage (new RTImage ("rtwindow.png")));
    hb->pack_start (*Gtk::manage (new Gtk::Label (Glib::path_get_basename (name))));
    hb->set_tooltip_markup (name);
    Gtk::Button* closeb = Gtk::manage (new Gtk::Button ());
    closeb->set_image (*Gtk::manage(new RTImage ("gtk-close.png")));
    closeb->set_relief (Gtk::RELIEF_NONE);
    closeb->set_focus_on_click (false);

    // make the button as small as possible thanks via css
    closeb->set_name("notebook_close_button");

    closeb->signal_clicked().connect( sigc::bind (sigc::mem_fun(*this, &EditWindow::remEditorPanel) , ep));
    hb->pack_end (*closeb);
    hb->set_spacing (2);
    hb->show_all ();

    mainNB->append_page (*ep, *hb);
    mainNB->set_current_page (mainNB->page_num (*ep));
    mainNB->set_tab_reorderable (*ep, true);

    set_title_decorated(name);

    epanels[ name ] = ep;
    filesEdited.insert ( name );
    parent->fpanel->refreshEditedState (filesEdited);
    ep->setAspect();
}

void EditWindow::remEditorPanel (EditorPanel* ep)
{
    if (ep->getIsProcessing()) {
        return;    // Will crash if destroyed while loading
    }

    epanels.erase (ep->getFileName());
    filesEdited.erase (ep->getFileName ());
    parent->fpanel->refreshEditedState (filesEdited);

    mainNB->remove_page (*ep);

    if (mainNB->get_n_pages() > 0) {
        EditorPanel* ep1 = static_cast<EditorPanel*>(mainNB->get_nth_page (mainNB->get_current_page()));
        set_title_decorated(ep1->getFileName());
    } else {
        set_title_decorated("");
    }

    // TODO: save options if wanted
}

bool EditWindow::selectEditorPanel(const std::string &name)
{
    std::map<Glib::ustring, EditorPanel*>::iterator iep = epanels.find(name);

    if (iep != epanels.end()) {
        mainNB->set_current_page (mainNB->page_num (*iep->second));
        set_title_decorated(name);
        return true;
    }

    return false;
}

void EditWindow::toFront ()
{
    // when using MEOW mode on a single monitor we need to present the secondary window.
    // If we don't it will stay in background when opening 2nd, 3rd... editor, which is annoying
    // It will also deiconify the window
    if(!isMultiDisplayEnabled()) {
        present();
    }
}

bool EditWindow::keyPressed (GdkEventKey* event)
{
    bool ctrl = event->state & GDK_CONTROL_MASK;

    if(event->keyval == GDK_KEY_F11) {
        toggleFullscreen();
        return true;
    } else {
        if(mainNB->get_n_pages () > 0) { //pass the handling for the editor panels, if there are any
            if (event->keyval == GDK_KEY_w && ctrl) { //remove editor panel
                EditorPanel* ep = static_cast<EditorPanel*>(mainNB->get_nth_page (mainNB->get_current_page()));
                remEditorPanel (ep);
                return true;
            } else if(mainNB->get_n_pages () > 0) {
                EditorPanel* ep = static_cast<EditorPanel*>(mainNB->get_nth_page (mainNB->get_current_page()));
                return ep->handleShortcutKey (event);
            }
        }

        return false;
    }

}

void EditWindow::toggleFullscreen ()
{
    isFullscreen ? unfullscreen() : fullscreen();
    options.meowFullScreen = isFullscreen = !isFullscreen;
}

void EditWindow::writeOptions() {

    if(is_visible()) {
        if(isMultiDisplayEnabled()) {
            options.meowMonitor = get_screen()->get_monitor_at_window(get_window());
        }

        options.meowMaximized = is_maximized();
        get_position(options.meowX, options.meowY);
        get_size(options.meowWidth,options.meowHeight);
    }
}
bool EditWindow::on_delete_event(GdkEventAny* event)
{

    if (!closeOpenEditors()) {
        return true;
    }

    writeOptions();
    hide();
    isClosed = true;

    return false;
}

bool EditWindow::isProcessing ()
{
    for ( std::set <Glib::ustring>::iterator iter = filesEdited.begin(); iter != filesEdited.end(); ++iter ) {
        if (epanels[*iter]->getIsProcessing()) {
            return true;
        }
    }

    return false;
}

bool EditWindow::closeOpenEditors()
{
    // Check if any editor is still processing, and do NOT quit if so. Otherwise crashes and inconsistent caches
    if (isProcessing()) {
        return false;
    }

    if (epanels.size()) {
        int page = mainNB->get_current_page();
        Gtk::Widget *w = mainNB->get_nth_page(page);
        bool optionsWritten = false;

        for (std::map<Glib::ustring, EditorPanel*>::iterator i = epanels.begin(); i != epanels.end(); ++i) {
            if (i->second == w) {
                i->second->writeOptions();
                optionsWritten = true;
            }
        }

        if (!optionsWritten) {
            // fallback solution: save the options of the first editor panel
            std::map<Glib::ustring, EditorPanel*>::iterator i = epanels.begin();
            i->second->writeOptions();
        }
    }

    for ( std::set <Glib::ustring>::iterator iter = filesEdited.begin(); iter != filesEdited.end(); ++iter ) {
        mainNB->remove_page (*epanels[*iter]);
    }

    epanels.clear();
    filesEdited.clear();
    parent->fpanel->refreshEditedState (filesEdited);

    return true;
}

void EditWindow::set_title_decorated(Glib::ustring fname)
{
    Glib::ustring subtitle;

    if (!fname.empty()) {
        subtitle = " - " + fname;
    }

    set_title("RawTherapee " + M("EDITWINDOW_TITLE") + subtitle);
}
