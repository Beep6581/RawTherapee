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
#include <dirselectionlistener.h>
#include <filebrowser.h>
#include <procthread.h>
#include <exiffiltersettings.h>
#include <giomm.h>
#include <fileselectionlistener.h>
#include <set>
#include <fileselectionchangelistener.h>
#include <coarsepanel.h>
#include <toolbar.h>
#include <filterpanel.h>

class PreviewLoaderListener {
  public:
    virtual void previewReady (FileBrowserEntry* fd) {}
    virtual void previewsFinished () {}
};

class DirEntry {

    public:
        Glib::ustring fullName;
        
        DirEntry (const Glib::ustring& n) : fullName (n) {}
		
	bool operator< (DirEntry& other) {
		return fullName.casefold() < other.fullName.casefold();
	}
};

class PreviewLoader : public ProcessingThread<DirEntry> {

  protected:
    PreviewLoaderListener* pl;

  public:
    PreviewLoader () : pl(NULL) { ProcessingThread<DirEntry>(); }
    void setPreviewLoaderListener (PreviewLoaderListener* p) { pl = p; }
    void start ();
    void process () { ProcessingThread<DirEntry>::process (); }
    void process (DirEntry& current);
	void remove  (Glib::ustring fname);
    void end ();
};

class FileCatalog : public Gtk::VBox, 
                    public DirSelectionListener, 
                    public PreviewLoaderListener, 
					public FilterPanelListener,
                    public FileBrowserListener
#ifdef _WIN32
                  , public WinDirChangeListener
#endif
 {

        // thumbnail browsers
        FileBrowser* fileBrowser;
        
        Gtk::HBox* hBox;
        Glib::ustring selectedDirectory;
        bool enabled;
        
        PreviewLoader previewLoader;
        FileSelectionListener* listener;
        FileSelectionChangeListener* fslistener;
        ImageAreaToolListener* iatlistener;

        Gtk::HBox* buttonBar;
        Gtk::HBox* buttonBar2;
        Gtk::ToggleButton* bDir;
        Gtk::ToggleButton* bUnRanked;
        Gtk::ToggleButton* bRank[5];
        Gtk::ToggleButton* bTrash;
        Gtk::ToggleButton* categoryButtons[8];
        Gtk::ToggleButton* exifInfo;
        sigc::connection bCateg[8];
        Gtk::Image* iranked[5], *igranked[5];
        
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

        Gtk::ProgressBar* progressBar;
        int previewsToLoad;
        int previewsLoaded;


#ifdef _WIN32
        WinDirMonitor* wdMonitor;
     public:
        int checkCounter;
        void winDirChanged ();
     private:
#endif		
        std::vector<Glib::ustring> fileNameList;
        std::set<Glib::ustring> editedFiles;

        void addAndOpenFile (const Glib::ustring& fname);
        void checkAndAddFile (Glib::RefPtr<Gio::File> info);
        std::vector<Glib::ustring> getFileList ();
        BrowserFilter getFilter ();
               
    public:
            CoarsePanel* coarsePanel;
            ToolBar* toolBar;

                     FileCatalog (CoarsePanel* cp, ToolBar* tb);
                void dirSelected (const Glib::ustring& dirname, const Glib::ustring& openfile="");
                void closeDir    ();
                void refreshEditedState (const std::set<Glib::ustring>& efiles);
                
                // previewloaderlistener interface
				void previewReady (FileBrowserEntry* fdn);
				void previewsFinished ();
                void _previewsFinished ();
                void _refreshProgressBar ();

				// filterpanel interface
				void exifFilterChanged ();
				
       Glib::ustring lastSelectedDir () { return selectedDirectory; }
                void setEnabled (bool e);   // if not enabled, it does not open image

                void redrawAll ();
                void refreshAll ();
                
                void openRequested          (std::vector<Thumbnail*> tbe);
                void deleteRequested        (std::vector<FileBrowserEntry*> tbe);
                void developRequested       (std::vector<FileBrowserEntry*> tbe);
                void renameRequested        (std::vector<FileBrowserEntry*> tbe);
                void selectionChanged       (std::vector<Thumbnail*> tbe);
                void emptyTrash ();
                
                void setFileSelectionListener (FileSelectionListener* l) { listener = l; }
                void setFileSelectionChangeListener (FileSelectionChangeListener* l) { fslistener = l; }
                void setImageAreaToolListener (ImageAreaToolListener* l) { iatlistener = l; }
				void setFilterPanel (FilterPanel* fpanel);
				void exifInfoButtonToggled();
                void categoryButtonToggled (Gtk::ToggleButton* b);
                void filterChanged ();
                void runFilterDialog ();

                void on_realize();
                void on_dir_changed (const Glib::RefPtr<Gio::File>& file, const Glib::RefPtr<Gio::File>& other_file, Gio::FileMonitorEvent event_type, bool internal);
				int  reparseDirectory ();
                void _openImage (std::vector<Thumbnail*> tmb);
                
                void zoomIn ();
                void zoomOut ();

                void openNextImage () { fileBrowser->openNextImage(); }
                void openPrevImage () { fileBrowser->openPrevImage(); }               
};

#endif
