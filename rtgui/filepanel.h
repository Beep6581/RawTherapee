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
#ifndef _FILEPANEL_
#define _FILEPANEL_

#include <gtkmm.h>
#include <batchtoolpanelcoord.h>
#include <filecatalog.h>
#include <dirbrowser.h>
#include <fileselectionlistener.h>
#include <placesbrowser.h>
#include <recentbrowser.h>
#include <pparamschangelistener.h>
#include <history.h>
#include <filterpanel.h>

class RTWindow;
class FilePanel : public Gtk::HPaned,
                  public FileSelectionListener,  
                  public PParamsChangeListener
{

    protected:
        Gtk::VPaned* placespaned;
        Gtk::HPaned* dirpaned;
        DirBrowser* dirBrowser;
        PlacesBrowser* placesBrowser;
        RecentBrowser* recentBrowser;
        FileCatalog* fileCatalog;   // filecatalog is the file browser with the button bar above it
        Gtk::HBox* rightBox;
        BatchToolPanelCoordinator* tpc;
        History* history;
		FilterPanel* filterPanel;
        RTWindow* parent;      
        Gtk::Notebook* rightNotebook;

    public:
        FilePanel ();

        void on_realize ();

        void setParent (RTWindow* p) { parent = p; }
        void init (); // dont call it directly, the constructor calls it as idle source
        void open (const Glib::ustring& d); // open a file or a directory
        void refreshEditedState (const std::set<Glib::ustring>& efiles) { fileCatalog->refreshEditedState (efiles); }
        
        // call this before closeing rt: it saves file browser relating things into options
        void saveOptions ();
        
        // interface fileselectionlistener
        bool fileSelected           (Thumbnail* thm);
        bool addBatchQueueJob       (BatchQueueEntry* bqe);

        void optionsChanged         ();
};

#endif

