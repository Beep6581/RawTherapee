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
#ifndef _FAVORITBROWSER_
#define _FAVORITBROWSER_

#include <gtkmm.h>
#include "dirbrowserremoteinterface.h"
#include "dirselectionlistener.h"

class FavoritBrowser : public Gtk::VBox, public DirSelectionListener {

        class FavoritColumns : public Gtk::TreeModel::ColumnRecord {
            public:
                Gtk::TreeModelColumn<Glib::RefPtr<Gio::Icon> >   icon;
                Gtk::TreeModelColumn<Glib::ustring> shortdir;
                Gtk::TreeModelColumn<Glib::ustring> fulldir;
                FavoritColumns() { add(icon); add(shortdir), add(fulldir); }
        };
        
        FavoritColumns               favoritColumns;
        Gtk::ScrolledWindow*         scrollw;
        Gtk::TreeView*               treeView;
        Glib::RefPtr<Gtk::ListStore> favoritModel;
        DirBrowserRemoteInterface*   listener;
        Glib::ustring                lastSelectedDir;
        Gtk::Button*                 add;
        Gtk::Button*                 del;
    public:
    
        FavoritBrowser ();
        
        void setDirBrowserRemoteInterface (DirBrowserRemoteInterface* l) { listener = l; }
        void dirSelected (const Glib::ustring& dirname, const Glib::ustring& openfile="");
        
        void addPressed ();
        void delPressed ();
        void selectionChanged ();
};

#endif


