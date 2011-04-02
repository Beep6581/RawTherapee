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
#ifndef _FILECATALOG_
#define _FILECATALOG_

#ifdef _WIN32
#include <windirmonitor.h>
#endif
#include <dirbrowserremoteinterface.h>
#include <dirselectionlistener.h>
#include <filebrowser.h>
#include <exiffiltersettings.h>
#include <giomm.h>
#include <fileselectionlistener.h>
#include <set>
#include <fileselectionchangelistener.h>
#include <coarsepanel.h>
#include <toolbar.h>
#include <filterpanel.h>
#include <previewloader.h>
#include <multilangmgr.h>


class DirEntry {

    public:
        Glib::ustring fullName;
        
        DirEntry (const Glib::ustring& n) : fullName (n) {}
		
	bool operator< (DirEntry& other) {
		return fullName.casefold() < other.fullName.casefold();
	}
};
class FilePanel;
class FileCatalog : public Gtk::VBox,
                    public DirSelectionListener, 
                    public PreviewLoaderListener, 
					public FilterPanelListener,
                    public FileBrowserListener
#ifdef _WIN32
                  , public WinDirChangeListener
#endif
 {

        FilePanel* filepanel;
        Gtk::HBox* hBox;
        Glib::ustring selectedDirectory;
        int selectedDirectoryId;
        bool enabled;
        bool inTabMode;  // Tab mode has e.g. different progress bar handling

        FileSelectionListener* listener;
        FileSelectionChangeListener* fslistener;
        ImageAreaToolListener* iatlistener;
        DirBrowserRemoteInterface*   dirlistener;

        Gtk::HBox* buttonBar;
        Gtk::HBox* buttonBar2;
        Gtk::ToggleButton* tbLeftPanel_1;
        Gtk::ToggleButton* tbRightPanel_1;
        Gtk::ToggleButton* bDir;
        Gtk::ToggleButton* bUnRanked;
        Gtk::ToggleButton* bRank[5];
        Gtk::ToggleButton* bTrash;
        Gtk::ToggleButton* categoryButtons[8];
        Gtk::ToggleButton* exifInfo;
        sigc::connection bCateg[8];
        Gtk::Image* iranked[5], *igranked[5];
        Gtk::Image *iTrashEmpty, *iTrashFull;
        Gtk::Image *iRightArrow_red, *iRightArrow;
        Gtk::Image *iLeftPanel_1_Show, *iLeftPanel_1_Hide, *iRightPanel_1_Show, *iRightPanel_1_Hide;
        Gtk::Entry* BrowsePath;
        Gtk::Button* buttonBrowsePath;
        sigc::connection BrowsePathconn;
        
        double hScrollPos[8];
        double vScrollPos[8];
        int lastScrollPos;

        Gtk::VBox* trashButtonBox;
        
        Gtk::Button* zoomInButton;
        Gtk::Button* zoomOutButton;

        ExifFilterSettings dirEFS;
        ExifFilterSettings currentEFS;
        bool hasValidCurrentEFS;  

		FilterPanel* filterPanel;

        Glib::RefPtr<Gio::FileMonitor> dirMonitor;

        int previewsToLoad;
        int previewsLoaded;


#ifdef _WIN32
        WinDirMonitor* wdMonitor;
     public:
        void winDirChanged ();
     private:
#endif		
        std::vector<Glib::ustring> fileNameList;
        std::set<Glib::ustring> editedFiles;
        guint modifierKey; // any modifiers held when rank button was pressed

        void addAndOpenFile (const Glib::ustring& fname);
        void checkAndAddFile (Glib::RefPtr<Gio::File> info);
        std::vector<Glib::ustring> getFileList ();
        BrowserFilter getFilter ();
        void trashChanged ();

    public:
            // thumbnail browsers
            FileBrowser* fileBrowser;
            
            CoarsePanel* coarsePanel;
            ToolBar* toolBar;

                     FileCatalog (CoarsePanel* cp, ToolBar* tb, FilePanel* filepanel);
                void dirSelected (const Glib::ustring& dirname, const Glib::ustring& openfile="");
                void closeDir    ();
                void refreshEditedState (const std::set<Glib::ustring>& efiles);
                
                // previewloaderlistener interface
				void _previewReady (int dir_id, FileBrowserEntry* fdn);
				void previewReady (int dir_id, FileBrowserEntry* fdn);
				void previewsFinished (int dir_id);
                void _previewsFinished ();
                void _refreshProgressBar ();

				// filterpanel interface
				void exifFilterChanged ();
				
       Glib::ustring lastSelectedDir () { return selectedDirectory; }
                void setEnabled (bool e);   // if not enabled, it does not open image
                void enableTabMode(bool enable);  // sets progress bar

                void redrawAll ();
                void refreshAll ();
                void refreshHeight ();
                
                void openRequested          (std::vector<Thumbnail*> tbe);
                void deleteRequested        (std::vector<FileBrowserEntry*> tbe, bool inclBatchProcessed);
                void copyMoveRequested      (std::vector<FileBrowserEntry*> tbe, bool moveRequested);
                void developRequested       (std::vector<FileBrowserEntry*> tbe);
                void renameRequested        (std::vector<FileBrowserEntry*> tbe);
                void clearFromCacheRequested(std::vector<FileBrowserEntry*> tbe, bool leavenotrace);
                void selectionChanged       (std::vector<Thumbnail*> tbe);
                void emptyTrash ();
                bool trashIsEmpty ();
                
                void setFileSelectionListener (FileSelectionListener* l) { listener = l; }
                void setFileSelectionChangeListener (FileSelectionChangeListener* l) { fslistener = l; }
                void setImageAreaToolListener (ImageAreaToolListener* l) { iatlistener = l; }
                void setDirBrowserRemoteInterface (DirBrowserRemoteInterface* l) { dirlistener = l; }

				void setFilterPanel (FilterPanel* fpanel);
				void exifInfoButtonToggled();
                void categoryButtonToggled (Gtk::ToggleButton* b);
                bool capture_event(GdkEventButton* event);
                void filterChanged ();
                void runFilterDialog ();

                void on_realize();
                void on_dir_changed (const Glib::RefPtr<Gio::File>& file, const Glib::RefPtr<Gio::File>& other_file, Gio::FileMonitorEvent event_type, bool internal);
				int  reparseDirectory ();
                void _openImage (std::vector<Thumbnail*> tmb);
                
                void zoomIn ();
                void zoomOut ();

                void buttonBrowsePathPressed ();

                void tbLeftPanel_1_toggled ();
                void tbLeftPanel_1_visible (bool visible);
                void tbRightPanel_1_toggled ();
                void tbRightPanel_1_visible (bool visible);

                void openNextImage () { fileBrowser->openNextImage(); }
                void openPrevImage () { fileBrowser->openPrevImage(); }               

                bool handleShortcutKey (GdkEventKey* event);

                bool CheckSidePanelsVisibility();
                void toggleSidePanels();
};

#endif
