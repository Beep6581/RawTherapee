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

#include "exportpanel.h"
#include "filecatalog.h"
#include "fileselectionlistener.h"
#include "filterpanel.h"
#include "history.h"
#include "placesbrowser.h"
#include "pparamschangelistener.h"
#include "progressconnector.h"
#include "recentbrowser.h"

#include "rtengine/noncopyable.h"

class BatchToolPanelCoordinator;
class RTWindow;
class DirBrowser;

class FilePanel final :
    public Gtk::Paned,
    public FileSelectionListener,
    public rtengine::NonCopyable
{
public:
    FilePanel ();
    ~FilePanel () override;

    Gtk::Paned* placespaned;
    Gtk::Paned* dirpaned;

    Gtk::Box* rightBox;

    DirBrowser* dirBrowser;
    FilterPanel* filterPanel;
    ExportPanel* exportPanel;
    FileCatalog* fileCatalog;
    Gtk::Paned *ribbonPane;

    void setParent (RTWindow* p)
    {
        parent = p;
    }
    void init (); // don't call it directly, the constructor calls it as idle source
    void on_realize () override;
    void setAspect();
    void open (const Glib::ustring& d); // open a file or a directory
    void refreshEditedState (const std::set<Glib::ustring>& efiles)
    {
        fileCatalog->refreshEditedState (efiles);
    }
    void loadingThumbs(Glib::ustring str, double rate);

    // call this before closing RT: it saves file browser's related things into options
    void saveOptions ();

    // interface fileselectionlistener
    bool fileSelected(Thumbnail* thm) override;
    bool addBatchQueueJobs(const std::vector<BatchQueueEntry*>& entries) override;

    void optionsChanged         ();
    bool imageLoaded( Thumbnail* thm, ProgressConnector<rtengine::InitialImage*> * );

    bool handleShortcutKey (GdkEventKey* event);
    bool handleShortcutKeyRelease(GdkEventKey *event);
    void updateTPVScrollbar (bool hide);
    void updateToolPanelToolLocations(
        const std::vector<Glib::ustring> &favorites, bool cloneFavoriteTools);

private:
    void on_NB_switch_page(Gtk::Widget* page, guint page_num);

    PlacesBrowser* placesBrowser;
    RecentBrowser* recentBrowser;

    Inspector* inspectorPanel;
    Gtk::Paned* tpcPaned;
    BatchToolPanelCoordinator* tpc;
    History* history;
    RTWindow* parent;
    Gtk::Notebook* rightNotebook;
    sigc::connection rightNotebookSwitchConn;

    struct pendingLoad {
        bool complete;
        ProgressConnector<rtengine::InitialImage*> *pc;
        Thumbnail *thm;
    };
    MyMutex pendingLoadMutex;
    std::vector<struct pendingLoad*> pendingLoads;

    int error;

    IdleRegister idle_register;
};
