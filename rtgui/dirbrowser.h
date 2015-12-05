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
#ifndef _DIRBROWSER_
#define _DIRBROWSER_

#include <gtkmm.h>
#include <giomm.h>
#ifdef WIN32
#include "windirmonitor.h"
#endif
#include "dirbrowserremoteinterface.h"

class DirBrowser : public Gtk::VBox, public DirBrowserRemoteInterface
#ifdef WIN32
    , public WinDirChangeListener
#endif
{
public:
    typedef sigc::signal<void, const Glib::ustring&, const Glib::ustring&> DirSelectionSignal;

private:

    Glib::RefPtr<Gtk::TreeStore> dirTreeModel;

    struct DirTreeColumns : public Gtk::TreeModelColumnRecord {
    public:
        Gtk::TreeModelColumn<Glib::ustring> filename;
        Gtk::TreeModelColumn<Glib::RefPtr<Gdk::Pixbuf> > icon1;
        Gtk::TreeModelColumn<Glib::RefPtr<Gdk::Pixbuf> > icon2;
        Gtk::TreeModelColumn<Glib::ustring> dirname;
#ifdef WIN32
        Gtk::TreeModelColumn<Glib::RefPtr<WinDirMonitor>  > monitor;
#else
        Gtk::TreeModelColumn<Glib::RefPtr<Gio::FileMonitor> > monitor;
#endif

        DirTreeColumns()
        {
            add(icon1);
            add(icon2);
            add(filename);
            add(dirname);
            add(monitor);
        }
    };

    DirTreeColumns dtColumns;
    Gtk::TreeViewColumn tvc;
    Gtk::CellRendererText crt;
    Gtk::CellRendererPixbuf crb;


    Gtk::TreeView *dirtree;
    Gtk::ScrolledWindow *scrolledwindow4;
    DirSelectionSignal dirSelectionSignal;

    void fillRoot ();

    Glib::RefPtr<Gdk::Pixbuf> openfolder;
    Glib::RefPtr<Gdk::Pixbuf> closedfolder;
    Glib::RefPtr<Gdk::Pixbuf> icdrom;
    Glib::RefPtr<Gdk::Pixbuf> ifloppy;
    Glib::RefPtr<Gdk::Pixbuf> ihdd;
    Glib::RefPtr<Gdk::Pixbuf> inetwork;
    Glib::RefPtr<Gdk::Pixbuf> iremovable;

    bool expandSuccess;

#ifdef WIN32
    int volumes;
public:
    void updateVolumes ();
    void updateDirTree  (const Gtk::TreeModel::iterator& iter);
    void updateDirTreeRoot  ();
    void winDirChanged ();
private:
    void addRoot (char letter);
#endif
    void addDir (const Gtk::TreeModel::iterator& iter, const Glib::ustring& dirname);
    Gtk::TreePath expandToDir (const Glib::ustring& dirName);
    void updateDir (const Gtk::TreeModel::iterator& iter);

public:
    DirBrowser ();

    void fillDirTree ();
    void on_sort_column_changed() const;
    void row_expanded   (const Gtk::TreeModel::iterator& iter, const Gtk::TreeModel::Path& path);
    void row_activated  (const Gtk::TreeModel::Path& path, Gtk::TreeViewColumn* column);
    void file_changed   (const Glib::RefPtr<Gio::File>& file, const Glib::RefPtr<Gio::File>& other_file, Gio::FileMonitorEvent event_type, const Gtk::TreeModel::iterator& iter, const Glib::ustring& dirName);
    void open           (const Glib::ustring& dirName, const Glib::ustring& fileName = ""); // goes to dir "dirName" and selects file "fileName"
    void selectDir      (Glib::ustring dir);

    DirSelectionSignal dirSelected () const;
};

inline DirBrowser::DirSelectionSignal DirBrowser::dirSelected () const
{
    return dirSelectionSignal;
}

#endif
