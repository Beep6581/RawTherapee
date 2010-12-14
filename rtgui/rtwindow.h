/*
 *  This file is part of RawTherapee.
 *
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
#ifndef _RTWINDOW_
#define _RTWINDOW_

#include <gtkmm.h>
#include <filepanel.h>
#include <editorpanel.h>
#include <batchqueuepanel.h>
#include <set>
#include <progressconnector.h>
#include <editwindow.h>

class RTWindow : public Gtk::Window, public rtengine::ProgressListener{

    private:
        Gtk::Notebook* mainNB;
        BatchQueuePanel* bpanel;
        std::set<Glib::ustring> filesEdited;
        std::map<Glib::ustring, EditorPanel*> epanels;

        Gtk::Label prLabel;
        Gtk::ProgressBar prProgBar;
        PLDBridge* pldBridge;
        bool is_fullscreen;
        Gtk::Button * btn_fullscreen;
        
        bool isSingleTabMode() { return !options.tabbedUI && !EditWindow::isMultiDisplayEnabled(); };

        bool on_expose_event_epanel(GdkEventExpose* event);
        bool on_expose_event_fpanel(GdkEventExpose* event);

    public:
        RTWindow ();

        void addEditorPanel (EditorPanel* ep,const std::string &name);
        void remEditorPanel (EditorPanel* ep);

        void addBatchQueueJob       (BatchQueueEntry* bqe, bool head=false);
        void addBatchQueueJobs      (std::vector<BatchQueueEntry*> &entries);

        bool keyPressed (GdkEventKey* event);
        bool on_delete_event(GdkEventAny* event);
        bool on_window_state_event(GdkEventWindowState* event);
        void on_mainNB_switch_page(GtkNotebookPage* page, guint page_num);

        void imageDeveloped (Glib::ustring fname); // called by the batchqueue when it finishes an image
        void showPreferences ();
        void on_realize ();
        void toggle_fullscreen ();
        void setProgress (double p);
        void setProgressStr (Glib::ustring str);
        void setProgressState (int state);
        void error (Glib::ustring descr);
        rtengine::ProgressListener* getProgressListener () { return pldBridge; }
        
        EditorPanel*  epanel;
        FilePanel* fpanel;

        void SetEditorCurrent();
        void SetMainCurrent();
        void MoveFileBrowserToEditor();
        void MoveFileBrowserToMain();
};

#endif
