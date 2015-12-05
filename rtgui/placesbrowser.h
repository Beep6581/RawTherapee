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
#ifndef _PLACESBROWSER_
#define _PLACESBROWSER_

#include <gtkmm.h>
#include <giomm.h>
#include "multilangmgr.h"

class PlacesBrowser : public Gtk::VBox
{
public:
    typedef sigc::slot<void, const Glib::ustring&> DirSelectionSlot;

private:

    class PlacesColumns : public Gtk::TreeModel::ColumnRecord
    {
    public:
        Gtk::TreeModelColumn<Glib::RefPtr<Gio::Icon> >   icon;
        Gtk::TreeModelColumn<Glib::ustring>              label;
        Gtk::TreeModelColumn<Glib::ustring>              root;
        Gtk::TreeModelColumn<int>                        type;
        Gtk::TreeModelColumn<bool>                       rowSeparator;
        PlacesColumns()
        {
            add(icon);
            add(label);
            add(root);
            add(type);
            add(rowSeparator);
        }
    };
    PlacesColumns            placesColumns;
    Gtk::ScrolledWindow*    scrollw;
    Gtk::TreeView*          treeView;
    Glib::RefPtr<Gtk::ListStore> placesModel;
    Glib::RefPtr<Gio::VolumeMonitor> vm;
    DirSelectionSlot             selectDir;
    Glib::ustring                lastSelectedDir;
    Gtk::Button*                 add;
    Gtk::Button*                 del;

public:

    PlacesBrowser ();

    void setDirSelector (const DirSelectionSlot& selectDir);
    void dirSelected (const Glib::ustring& dirname, const Glib::ustring& openfile);

    void refreshPlacesList ();
    void mountChanged (const Glib::RefPtr<Gio::Mount>& m);
    void volumeChanged (const Glib::RefPtr<Gio::Volume>& v);
    void driveChanged (const Glib::RefPtr<Gio::Drive>& d);
    bool rowSeparatorFunc (const Glib::RefPtr<Gtk::TreeModel>& model, const Gtk::TreeModel::iterator& iter);
    void selectionChanged ();
    void addPressed ();
    void delPressed ();
};

inline void PlacesBrowser::setDirSelector (const PlacesBrowser::DirSelectionSlot& selectDir)
{
    this->selectDir = selectDir;
}

#endif


