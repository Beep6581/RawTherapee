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

#ifdef WIN32
#include "windirmonitor.h"
#endif
#include "filebrowser.h"
#include "exiffiltersettings.h"
#include <giomm.h>
#include "fileselectionlistener.h"
#include <set>
#include "fileselectionchangelistener.h"
#include "coarsepanel.h"
#include "toolbar.h"
#include "filterpanel.h"
#include "exportpanel.h"
#include "previewloader.h"
#include "multilangmgr.h"
#include "threadutils.h"


class DirEntry
{

public:
    Glib::ustring fullName;

    explicit DirEntry (const Glib::ustring& n) : fullName (n) {}

    bool operator< (DirEntry& other)
    {
        return fullName.casefold() < other.fullName.casefold();
    }
};
class FilePanel;
/*
 * Class:
 *   - handling the list of file (add/remove them)
 *   - handling the thumbnail toolbar,
 *   - monitoring the directory (for any change)
 */
class FileCatalog : public Gtk::VBox,
    public PreviewLoaderListener,
    public FilterPanelListener,
    public FileBrowserListener,
    public ExportPanelListener
#ifdef WIN32
    , public WinDirChangeListener
#endif
{
public:
    typedef sigc::slot<void, const Glib::ustring&> DirSelectionSlot;

private:
    FilePanel* filepanel;
    Gtk::HBox* hBox;
    Glib::ustring selectedDirectory;
    int selectedDirectoryId;
    bool enabled;
    bool inTabMode;  // Tab mode has e.g. different progress bar handling
    Glib::ustring imageToSelect_fname;
    Glib::ustring refImageForOpen_fname; // Next/previous for Editor's perspective
    eRTNav actionNextPrevious;

    FileSelectionListener* listener;
    FileSelectionChangeListener* fslistener;
    ImageAreaToolListener* iatlistener;
    DirSelectionSlot selectDir;

    Gtk::HBox* buttonBar;
    Gtk::HBox* hbToolBar1;

    Gtk::HBox* fltrRankbox;
    Gtk::HBox* fltrLabelbox;
    Gtk::VBox* fltrVbox1;

    Gtk::HBox* fltrEditedBox;
    Gtk::HBox* fltrRecentlySavedBox;
    Gtk::VBox* fltrVbox2;

    Gtk::VSeparator* vSepiLeftPanel;

    Gtk::ToggleButton* tbLeftPanel_1;
    Gtk::ToggleButton* tbRightPanel_1;
    Gtk::ToggleButton* bFilterClear;
    Gtk::ToggleButton* bUnRanked;
    Gtk::ToggleButton* bRank[5];
    Gtk::ToggleButton* bUnCLabeled;
    Gtk::ToggleButton* bCLabel[5];//color label
    Gtk::ToggleButton* bEdited[2];
    Gtk::ToggleButton* bRecentlySaved[2];
    Gtk::ToggleButton* bTrash;
    Gtk::ToggleButton* bNotTrash;
    Gtk::ToggleButton* bOriginal;
    Gtk::ToggleButton* categoryButtons[20];
    Gtk::ToggleButton* exifInfo;
    sigc::connection bCateg[20];
    Gtk::Image* iFilterClear, *igFilterClear;
    Gtk::Image* iranked[5], *igranked[5], *iUnRanked, *igUnRanked;
    Gtk::Image* iCLabeled[5], *igCLabeled[5], *iUnCLabeled, *igUnCLabeled;
    Gtk::Image* iEdited[2], *igEdited[2];
    Gtk::Image* iRecentlySaved[2], *igRecentlySaved[2];
    Gtk::Image *iTrashEmpty, *iTrashFull;
    Gtk::Image *iNotTrash, *iOriginal;
    Gtk::Image *iRefreshWhite, *iRefreshRed;
    Gtk::Image *iLeftPanel_1_Show, *iLeftPanel_1_Hide, *iRightPanel_1_Show, *iRightPanel_1_Hide;
    Gtk::Image *iQueryClear;

    Gtk::Entry* BrowsePath;
    Gtk::Button* buttonBrowsePath;

    Gtk::Entry* Query;
    Gtk::Button* buttonQueryClear;

    double hScrollPos[18];
    double vScrollPos[18];
    int lastScrollPos;

    Gtk::VBox* trashButtonBox;

    Gtk::Button* zoomInButton;
    Gtk::Button* zoomOutButton;

    MyMutex dirEFSMutex;
    ExifFilterSettings dirEFS;
    ExifFilterSettings currentEFS;
    bool hasValidCurrentEFS;

    FilterPanel* filterPanel;
    ExportPanel* exportPanel;

    int previewsToLoad;
    int previewsLoaded;


    std::vector<Glib::ustring> fileNameList;
    std::set<Glib::ustring> editedFiles;
    guint modifierKey; // any modifiers held when rank button was pressed

#ifndef _WIN32
    Glib::RefPtr<Gio::FileMonitor> dirMonitor;
#else
    WinDirMonitor* wdMonitor;
#endif

    IdleRegister idle_register;

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
    ~FileCatalog();
    void dirSelected (const Glib::ustring& dirname, const Glib::ustring& openfile);
    void closeDir    ();
    void refreshEditedState (const std::set<Glib::ustring>& efiles);

    // previewloaderlistener interface
    void previewReady (int dir_id, FileBrowserEntry* fdn);
    void previewsFinished (int dir_id);
    void previewsFinishedUI ();
    void _refreshProgressBar ();

    void setInspector(Inspector* inspector)
    {
        if (fileBrowser) {
            fileBrowser->setInspector(inspector);
        }
    }
    void disableInspector()
    {
        if (fileBrowser) {
            fileBrowser->disableInspector();
        }
    }
    void enableInspector()
    {
        if (fileBrowser) {
            fileBrowser->enableInspector();
        }
    }

    // filterpanel interface
    void exifFilterChanged ();

    // exportpanel interface
    void exportRequested();

    Glib::ustring lastSelectedDir ()
    {
        return selectedDirectory;
    }
    void setEnabled (bool e);   // if not enabled, it does not open image
    void enableTabMode(bool enable);  // sets progress bar

    // accessors for FileBrowser
    void redrawAll ();
    void refreshThumbImages ();
    void refreshHeight ();

    void openRequested          (std::vector<Thumbnail*> tbe);
    void deleteRequested        (std::vector<FileBrowserEntry*> tbe, bool inclBatchProcessed);
    void copyMoveRequested      (std::vector<FileBrowserEntry*> tbe, bool moveRequested);
    void developRequested       (std::vector<FileBrowserEntry*> tbe, bool fastmode);
    void renameRequested        (std::vector<FileBrowserEntry*> tbe);
    void clearFromCacheRequested(std::vector<FileBrowserEntry*> tbe, bool leavenotrace);
    void selectionChanged       (std::vector<Thumbnail*> tbe);
    void emptyTrash ();
    bool trashIsEmpty ();

    void setFileSelectionListener (FileSelectionListener* l)
    {
        listener = l;
    }
    void setFileSelectionChangeListener (FileSelectionChangeListener* l)
    {
        fslistener = l;
    }
    void setImageAreaToolListener (ImageAreaToolListener* l)
    {
        iatlistener = l;
    }
    void setDirSelector (const DirSelectionSlot& selectDir);

    void setFilterPanel (FilterPanel* fpanel);
    void setExportPanel (ExportPanel* expanel);
    void exifInfoButtonToggled();
    void categoryButtonToggled (Gtk::ToggleButton* b, bool isMouseClick);
    bool capture_event(GdkEventButton* event);
    void filterChanged ();
    void runFilterDialog ();

    void on_realize();
    void reparseDirectory ();
    void _openImage (std::vector<Thumbnail*> tmb);

    void zoomIn ();
    void zoomOut ();

    void buttonBrowsePathPressed ();
    bool BrowsePath_key_pressed (GdkEventKey *event);
    void buttonQueryClearPressed ();
    void executeQuery ();
    bool Query_key_pressed(GdkEventKey *event);
    void updateFBQueryTB (bool singleRow);
    void updateFBToolBarVisibility (bool showFilmStripToolBar);

    void tbLeftPanel_1_toggled ();
    void tbLeftPanel_1_visible (bool visible);
    void tbRightPanel_1_toggled ();
    void tbRightPanel_1_visible (bool visible);

    void openNextImage ()
    {
        fileBrowser->openNextImage();
    }
    void openPrevImage ()
    {
        fileBrowser->openPrevImage();
    }
    void selectImage (Glib::ustring fname, bool clearFilters);
    void openNextPreviousEditorImage (Glib::ustring fname, bool clearFilters, eRTNav nextPrevious);

    bool handleShortcutKey (GdkEventKey* event);

    bool isInTabMode()
    {
        return inTabMode;
    }

    bool CheckSidePanelsVisibility();
    void toggleSidePanels();
    void toggleLeftPanel();
    void toggleRightPanel();

    void showToolBar();
    void hideToolBar();
    void filterApplied();

#ifndef _WIN32
    void on_dir_changed (const Glib::RefPtr<Gio::File>& file, const Glib::RefPtr<Gio::File>& other_file, Gio::FileMonitorEvent event_type, bool internal);
#else
    void winDirChanged ();
#endif

};

inline void FileCatalog::setDirSelector (const FileCatalog::DirSelectionSlot& selectDir)
{
    this->selectDir = selectDir;
}

#endif
