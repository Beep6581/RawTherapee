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
#ifndef _FILEBROWSER_
#define _FILEBROWSER_

#include <gtkmm.h>
#include <map>
#include "thumbbrowserbase.h"
#include "exiffiltersettings.h"
#include "filebrowserentry.h"
#include "browserfilter.h"
#include "pparamschangelistener.h"
#include "partialpastedlg.h"
#include "exportpanel.h"
#include "extprog.h"
#include "profilestore.h"

class ProfileStoreLabel;
class FileBrowser;
class FileBrowserEntry;
class FileBrowserListener
{

public:
    virtual     ~FileBrowserListener    () {}
    virtual void filterApplied () {}
    virtual void openRequested          (std::vector<Thumbnail*> tbe) {}
    virtual void developRequested       (std::vector<FileBrowserEntry*> tbe, bool fastmode) {}
    virtual void renameRequested        (std::vector<FileBrowserEntry*> tbe) {}
    virtual void deleteRequested        (std::vector<FileBrowserEntry*> tbe, bool inclBatchProcessed) {}
    virtual void copyMoveRequested      (std::vector<FileBrowserEntry*> tbe, bool moveRequested) {}
    virtual void selectionChanged       (std::vector<Thumbnail*> tbe) {}
    virtual void clearFromCacheRequested(std::vector<FileBrowserEntry*> tbe, bool leavenotrace) {}
    virtual bool isInTabMode            ()
    {
        return false;
    }
};

struct FileBrowserIdleHelper {
    FileBrowser* fbrowser;
    bool destroyed;
    int pending;
};

/*
 * Class handling actions common to all thumbnails of the file browser
 */
class FileBrowser  : public ThumbBrowserBase,
    public LWButtonListener,
    public ExportPanelListener,
    public ProfileStoreListener
{
private:
    typedef sigc::signal<void> type_trash_changed;

    IdleRegister idle_register;

protected:
    Gtk::MenuItem* rank[6];
    MyImageMenuItem* colorlabel[6];
    Gtk::MenuItem* trash;
    Gtk::MenuItem* untrash;
    Gtk::MenuItem* develop;
    Gtk::MenuItem* developfast;
    Gtk::MenuItem* rename;
    Gtk::MenuItem* remove;
    Gtk::MenuItem* removeInclProc;
    Gtk::MenuItem* open;
    Gtk::MenuItem* selall;
    Gtk::MenuItem* copyTo;
    Gtk::MenuItem* moveTo;

    Gtk::MenuItem* menuRank;
    Gtk::MenuItem* menuLabel;
    Gtk::MenuItem* menuFileOperations;
    MyImageMenuItem* menuProfileOperations;
    Gtk::MenuItem* menuExtProg;
    Gtk::MenuItem** amiExtProg;
    Gtk::MenuItem* miOpenDefaultViewer;
    std::map<Glib::ustring, const ExtProgAction*> mMenuExtProgs;  // key is menuitem label

    Gtk::MenuItem* menuDF;
    Gtk::MenuItem* selectDF;
    Gtk::MenuItem* thisIsDF;
    Gtk::MenuItem* autoDF;

    Gtk::MenuItem* menuFF;
    Gtk::MenuItem* selectFF;
    Gtk::MenuItem* thisIsFF;
    Gtk::MenuItem* autoFF;

    Gtk::MenuItem* copyprof;
    Gtk::MenuItem* pasteprof;
    Gtk::MenuItem* partpasteprof;
    Gtk::MenuItem* applyprof;
    Gtk::MenuItem* applypartprof;
    Gtk::MenuItem* resetdefaultprof;
    Gtk::MenuItem* clearprof;
    Gtk::MenuItem* cachemenu;
    Gtk::MenuItem* clearFromCache;
    Gtk::MenuItem* clearFromCacheFull;
    Gtk::Menu* pmenu;

    MyImageMenuItem* colorlabel_pop[6];
    Gtk::Menu* pmenuColorLabels;
    void* colorLabel_actionData;
    void menuColorlabelActivated (Gtk::MenuItem* m); // use only when menu is invoked via FileBrowser::buttonPressed to pass actionData

    Glib::RefPtr<Gtk::AccelGroup> pmaccelgroup;

    BatchPParamsChangeListener* bppcl;
    FileBrowserListener* tbl;
    BrowserFilter filter;
    int numFiltered;
    FileBrowserIdleHelper* fbih;

    void toTrashRequested   (std::vector<FileBrowserEntry*> tbe);
    void fromTrashRequested (std::vector<FileBrowserEntry*> tbe);
    void rankingRequested   (std::vector<FileBrowserEntry*> tbe, int rank);
    void colorlabelRequested   (std::vector<FileBrowserEntry*> tbe, int colorlabel);
    void requestRanking (int rank);
    void requestColorLabel(int colorlabel);
    void notifySelectionListener ();
    void openRequested( std::vector<FileBrowserEntry*> mselected);
    ExportPanel* exportPanel;

    type_trash_changed m_trash_changed;

public:
    FileBrowser ();
    ~FileBrowser ();

    void addEntry (FileBrowserEntry* entry); // can be called from any thread
    void addEntry_ (FileBrowserEntry* entry); // this must be executed inside the gtk thread
    FileBrowserEntry*  delEntry (const Glib::ustring& fname);    // return the entry if found here return NULL otherwise
    void close ();

    void setBatchPParamsChangeListener (BatchPParamsChangeListener* l)
    {
        bppcl = l;
    }
    void setFileBrowserListener (FileBrowserListener* l)
    {
        tbl = l;
    }

    void menuItemActivated (Gtk::MenuItem* m);
    void applyMenuItemActivated (ProfileStoreLabel *label);
    void applyPartialMenuItemActivated (ProfileStoreLabel *label);

    void applyFilter (const BrowserFilter& filter);
    int getNumFiltered()
    {
        return numFiltered;
    }

    void buttonPressed (LWButton* button, int actionCode, void* actionData);
    void redrawNeeded  (LWButton* button);
    bool checkFilter (ThumbBrowserEntryBase* entry);
    void rightClicked (ThumbBrowserEntryBase* entry);
    void doubleClicked (ThumbBrowserEntryBase* entry);
    bool keyPressed (GdkEventKey* event);

    void saveThumbnailHeight (int height);
    int  getThumbnailHeight ();

    bool isInTabMode()
    {
        return tbl ? tbl->isInTabMode() : false;
    }

    void openNextImage ();
    void openPrevImage ();
    void copyProfile ();
    void pasteProfile ();
    void partPasteProfile ();
    void selectImage (Glib::ustring fname);
    void openNextPreviousEditorImage (Glib::ustring fname, eRTNav eNextPrevious);

    void openDefaultViewer (int destination);

    void thumbRearrangementNeeded ();
    void _thumbRearrangementNeeded ();

    void selectionChanged ();

    void setExportPanel (ExportPanel* expanel);
    // exportpanel interface
    void exportRequested();

    void updateProfileList ();

    type_trash_changed trash_changed();
};

#endif
