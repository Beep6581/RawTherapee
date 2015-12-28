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

// Check if the system has more than one display and option is set
bool EditWindow::isMultiDisplayEnabled()
{
    return options.multiDisplayMode > 0 && Gdk::Screen::get_default()->get_n_monitors () > 1;
}

// Should only be created once, auto-creates window on correct display
EditWindow* EditWindow::getInstance(RTWindow* p)
{
    struct EditWindowInstance
    {
        EditWindow editWnd;

        EditWindowInstance(RTWindow* p) : editWnd(p)
        {
            // Determine the other display and maximize the window on that
            const Glib::RefPtr< Gdk::Window >& wnd = p->get_window();
            int monNo = p->get_screen()->get_monitor_at_window (wnd);

            Gdk::Rectangle lMonitorRect;
            editWnd.get_screen()->get_monitor_geometry(isMultiDisplayEnabled() ? (monNo == 0 ?  1 : 0) : monNo, lMonitorRect);
            editWnd.move(lMonitorRect.get_x(), lMonitorRect.get_y());
            editWnd.maximize();
        }
    };

    static EditWindowInstance instance_(p);
    instance_.editWnd.show_all();
    return &instance_.editWnd;
}

EditWindow::EditWindow (RTWindow* p) : parent(p) , isFullscreen(false)
{

    Glib::ustring fName = "rt-logo.png";
    Glib::ustring fullPath = RTImage::findIconAbsolutePath(fName);

    try {
        set_default_icon_from_file (fullPath);
    } catch(Glib::Exception& ex) {
        printf ("%s\n", ex.what().c_str());
    }

    set_title_decorated("");
    property_allow_shrink() = true;
    set_modal(false);
    set_resizable(true);

    property_destroy_with_parent().set_value(false);
    signal_window_state_event().connect( sigc::mem_fun(*this, &EditWindow::on_window_state_event) );

    mainNB = Gtk::manage (new Gtk::Notebook ());
    mainNB->set_scrollable (true);
    mainNB->signal_switch_page().connect_notify( sigc::mem_fun(*this, &EditWindow::on_mainNB_switch_page) );

    signal_key_press_event().connect( sigc::mem_fun(*this, &EditWindow::keyPressed) );

    Gtk::VBox* mainBox = Gtk::manage (new Gtk::VBox ());
    mainBox->pack_start (*mainNB);

    add (*mainBox);
    show_all ();
}

void EditWindow::on_realize ()
{
    Gtk::Window::on_realize ();

    cursorManager.init (get_window());
}

bool EditWindow::on_window_state_event(GdkEventWindowState* event)
{
    if (!event->new_window_state) {
        // Window mode
        options.windowMaximized = false;
    } else if (event->new_window_state & (GDK_WINDOW_STATE_MAXIMIZED | GDK_WINDOW_STATE_FULLSCREEN)) {
        // Fullscreen mode
        options.windowMaximized = true;
    }

    return true;
}

void EditWindow::on_mainNB_switch_page(GtkNotebookPage* page, guint page_num)
{
    //if (page_num > 1) {
    EditorPanel *ep = static_cast<EditorPanel*>(mainNB->get_nth_page(page_num));

    if (mainNB->get_n_pages() > 1 && page_num <= (filesEdited.size() - 1)) {
        set_title_decorated(ep->getFileName());
    }

    ep->setAspect();
    //}
}

void EditWindow::addEditorPanel (EditorPanel* ep, const std::string &name)
{
    ep->setParent (parent);

    // construct closeable tab for the image
    Gtk::HBox* hb = Gtk::manage (new Gtk::HBox ());
    hb->pack_start (*Gtk::manage (new RTImage ("rtwindow.png")));
    hb->pack_start (*Gtk::manage (new Gtk::Label (Glib::path_get_basename (name))));
    hb->set_tooltip_markup (name);
    Gtk::Button* closeb = Gtk::manage (new Gtk::Button ());
    closeb->set_image (*Gtk::manage(new RTImage ("gtk-close.png")));
    closeb->set_relief (Gtk::RELIEF_NONE);
    closeb->set_focus_on_click (false);
    // make the button as small as possible
    Glib::RefPtr<Gtk::RcStyle> style = Gtk::RcStyle::create ();
    style->set_xthickness (0);
    style->set_ythickness (0);

    closeb->modify_style (style);
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

bool EditWindow::keyPressed (GdkEventKey* event)
{
    bool ctrl = event->state & GDK_CONTROL_MASK;

    if(event->keyval == GDK_F11) {
        toggleFullscreen();
        return true;
    } else {
        if(mainNB->get_n_pages () > 0) { //pass the handling for the editor panels, if there are any
            if (event->keyval == GDK_w && ctrl) { //remove editor panel
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
    isFullscreen = !isFullscreen;
}

bool EditWindow::on_delete_event(GdkEventAny* event)
{
    // Check if any editor is still processing, and do NOT quit if so. Otherwise crashes and inconsistent caches
    bool isProcessing = false;

    for ( std::set <Glib::ustring>::iterator iter = filesEdited.begin(); iter != filesEdited.end() && !isProcessing; iter++ ) {
        if (epanels[*iter]->getIsProcessing()) {
            isProcessing = true;
        }
    }

    if (isProcessing) {
        return true;
    }

    for ( std::set <Glib::ustring>::iterator iter = filesEdited.begin(); iter != filesEdited.end(); iter++ ) {
        mainNB->remove_page (*epanels[*iter]);
    }

    epanels.clear();

    filesEdited.clear();
    parent->fpanel->refreshEditedState (filesEdited);

    hide ();
    return false;
}

void EditWindow::set_title_decorated(Glib::ustring fname)
{
    Glib::ustring subtitle;

    if (!fname.empty()) {
        subtitle = " - " + fname;
    }

    set_title("RawTherapee " + M("EDITWINDOW_TITLE") + subtitle);
}
