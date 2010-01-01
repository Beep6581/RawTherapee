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
#include <thumbbrowserbase.h>
#include <exiffiltersettings.h>
#include <filebrowserentry.h>
#include <browserfilter.h>
#include <partialpastedlg.h>

class FileBrowser;
class FileBrowserEntry;
class FileBrowserListener {

    public:
        virtual void openRequested          (std::vector<Thumbnail*> tbe) {}
        virtual void developRequested       (std::vector<FileBrowserEntry*> tbe) {}
        virtual void renameRequested        (std::vector<FileBrowserEntry*> tbe) {}
        virtual void deleteRequested        (std::vector<FileBrowserEntry*> tbe) {}
        virtual void selectionChanged       (std::vector<Thumbnail*> tbe) {}
};

struct FileBrowserIdleHelper {
    FileBrowser* fbrowser;
    bool destroyed;
    int pending;
};

class FileBrowser  : public ThumbBrowserBase, public LWButtonListener {  

  protected:

    Gtk::MenuItem* rank[6];
    Gtk::MenuItem* trash;
    Gtk::MenuItem* untrash;
    Gtk::MenuItem* develop;
    Gtk::MenuItem* rename;
    Gtk::MenuItem* remove;
    Gtk::MenuItem* open;
    Gtk::MenuItem* selall;
    Gtk::MenuItem* copyprof;
    Gtk::MenuItem* pasteprof;
    Gtk::MenuItem* partpasteprof;
    Gtk::MenuItem* applyprof;
    Gtk::MenuItem* clearprof;
    Gtk::Menu* pmenu;
    Gtk::Menu* profmenu;
    Glib::RefPtr<Gtk::AccelGroup> pmaccelgroup;

    FileBrowserListener* tbl;
    BrowserFilter filter;
    PartialPasteDlg partialPasteDlg;

    FileBrowserIdleHelper* fbih;

    void toTrashRequested   (std::vector<FileBrowserEntry*> tbe);
    void fromTrashRequested (std::vector<FileBrowserEntry*> tbe);
    void rankingRequested   (std::vector<FileBrowserEntry*> tbe, int rank);
    void notifySelectionListener ();
    
  public:
   
    FileBrowser ();

    void addEntry (FileBrowserEntry* entry); // can be called from any thread
    void addEntry_ (FileBrowserEntry* entry); // this must be executed inside the gtk thread 
    FileBrowserEntry*  delEntry (const Glib::ustring& fname);    // return the entry if found here return NULL otherwise
    FileBrowserEntry*  findEntry (const Glib::ustring& fname);    // return the entry if found here return NULL otherwise
    void close ();
    
    void setFileBrowserListener (FileBrowserListener* l) { tbl = l; }
     
    void menuItemActivated (Gtk::MenuItem* m);
    void applyMenuItemActivated (Glib::ustring ppname);

    void applyFilter (const BrowserFilter& filter);

    void buttonPressed (LWButton* button, int actionCode, void* actionData);
    void redrawNeeded  (LWButton* button);
    bool checkFilter (ThumbBrowserEntryBase* entry);
    void rightClicked (ThumbBrowserEntryBase* entry);
    void doubleClicked (ThumbBrowserEntryBase* entry);
    bool keyPressed (GdkEventKey* event);

    void openNextImage ();
    void openPrevImage ();
    void copyProfile ();
    void pasteProfile ();
    void partPasteProfile ();

    void redrawNeeded (ThumbBrowserEntryBase* entry);
    void thumbRearrangementNeeded ();
    void _thumbRearrangementNeeded ();
    
    void selectionChanged ();
};

#endif
