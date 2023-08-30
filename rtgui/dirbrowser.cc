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
#include "dirbrowser.h"

#include <iostream>
#include <cstring>

#ifdef _WIN32_WINNT
#undef _WIN32_WINNT
#define _WIN32_WINNT 0x0600
#include <windows.h>
#endif

#include "guiutils.h"
#include "rtimage.h"
#include "multilangmgr.h"
#include "options.h"

namespace
{

std::vector<Glib::ustring> listSubDirs (const Glib::RefPtr<Gio::File>& dir, bool addHidden)
{
    std::vector<Glib::ustring> subDirs;

    try {

        // CD-ROM with no disc inserted are reported, but do not exist.
        if (!Glib::file_test (dir->get_path (), Glib::FILE_TEST_EXISTS)) {
            return subDirs;
        }

        auto enumerator = dir->enumerate_children ("standard::name,standard::type,standard::is-hidden");

        while (true) {
            try {
                auto file = enumerator->next_file ();
                if (!file) {
                    break;
                }
                if (file->get_file_type () != Gio::FILE_TYPE_DIRECTORY) {
                    continue;
                }
                if (!addHidden && file->is_hidden ()) {
                    continue;
                }
                subDirs.push_back (file->get_name ());
            } catch (const Glib::Exception& exception) {

                if (rtengine::settings->verbose) {
                    std::cerr << exception.what().c_str() << std::endl;
                }

            }
        }

    } catch (const Glib::Exception& exception) {

        if (rtengine::settings->verbose) {
            std::cerr << "Failed to list subdirectories of \"" << dir->get_parse_name() << "\": " << exception.what () << std::endl;
        }

    }

    return subDirs;
}

}

DirBrowser::DirBrowser () : dirTreeModel(),
    dtColumns(),
    tvc(M("DIRBROWSER_FOLDERS")),
    expandSuccess(false)
#ifdef _WIN32
    , volumes(0)
#endif
{
    set_orientation(Gtk::ORIENTATION_VERTICAL);
    dirtree = Gtk::manage ( new Gtk::TreeView() );
    scrolledwindow4 = Gtk::manage ( new Gtk::ScrolledWindow() );
    crt.property_ellipsize() = Pango::ELLIPSIZE_END;

//   dirtree->set_flags(Gtk::CAN_FOCUS);
    dirtree->set_headers_visible();
    dirtree->set_headers_clickable();
    dirtree->set_rules_hint(false);
    dirtree->set_reorderable(false);
    dirtree->set_enable_search(false);
    scrolledwindow4->set_can_focus(true);
    scrolledwindow4->set_shadow_type(Gtk::SHADOW_NONE);
    scrolledwindow4->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
    scrolledwindow4->property_window_placement().set_value(Gtk::CORNER_TOP_LEFT);
    scrolledwindow4->add(*dirtree);

    pack_start (*scrolledwindow4);
    dirtree->show ();
    scrolledwindow4->show ();
}

DirBrowser::~DirBrowser()
{
    idle_register.destroy();
}

void DirBrowser::fillDirTree ()
{

    openfolder = RTImage::createPixbufFromFile ("folder-open-small.png");
    closedfolder = RTImage::createPixbufFromFile ("folder-closed-small.png");
    icdrom = RTImage::createPixbufFromFile ("device-optical.png");
    ifloppy = RTImage::createPixbufFromFile ("device-floppy.png");
    ihdd = RTImage::createPixbufFromFile ("device-hdd.png");
    iremovable = RTImage::createPixbufFromFile ("device-usb.png");
    inetwork = RTImage::createPixbufFromFile ("device-network.png");

    //Create the Tree model:
    dirTreeModel = Gtk::TreeStore::create(dtColumns);
    dirtree->set_model (dirTreeModel);

    fillRoot ();

    Gtk::CellRendererPixbuf* render_pb = Gtk::manage ( new Gtk::CellRendererPixbuf () );
    tvc.pack_start (*render_pb, false);
    tvc.add_attribute(*render_pb, "pixbuf-expander-closed", dtColumns.icon2);
    tvc.add_attribute(*render_pb, "pixbuf", dtColumns.icon2);
    tvc.add_attribute(*render_pb, "pixbuf-expander-open", dtColumns.icon1);
    tvc.pack_start (crt);
    tvc.add_attribute(crt, "text", dtColumns.filename);

    dirtree->append_column(tvc);

    tvc.set_sort_order(options.dirBrowserSortType);
    tvc.set_sort_column(dtColumns.filename);
    tvc.set_sort_indicator(true);
    tvc.set_clickable();

    dirTreeModel->set_sort_column(dtColumns.filename, options.dirBrowserSortType);

    crt.property_ypad() = 0;
    render_pb->property_ypad() = 0;

    dirtree->signal_row_expanded().connect(sigc::mem_fun(*this, &DirBrowser::row_expanded));
    dirtree->signal_row_activated().connect(sigc::mem_fun(*this, &DirBrowser::row_activated));
    dirTreeModel->signal_sort_column_changed().connect(sigc::mem_fun(*this, &DirBrowser::on_sort_column_changed));
}

#ifdef _WIN32
void DirBrowser::addRoot (char letter)
{

    char volume[4];
    volume[0] = letter;
    strcpy (volume + 1, ":\\");

    Gtk::TreeModel::iterator root = dirTreeModel->append();
    root->set_value (dtColumns.filename, Glib::ustring(volume));
    root->set_value (dtColumns.dirname, Glib::ustring(volume));

    int type = GetDriveType (volume);

    if (type == DRIVE_CDROM) {
        root->set_value (0, icdrom);
        root->set_value (1, icdrom);
    } else if (type == DRIVE_REMOVABLE) {
        if (letter - 'A' < 2) {
            root->set_value (0, ifloppy);
            root->set_value (1, ifloppy);
        } else {
            root->set_value (0, iremovable);
            root->set_value (1, iremovable);
        }
    } else if (type == DRIVE_REMOTE) {
        root->set_value (0, inetwork);
        root->set_value (1, inetwork);
    } else if (type == DRIVE_FIXED) {
        root->set_value (0, ihdd);
        root->set_value (1, ihdd);
    }

    Gtk::TreeModel::iterator child = dirTreeModel->append (root->children());
    child->set_value (dtColumns.filename, Glib::ustring("foo"));
}

void DirBrowser::updateDirTreeRoot ()
{

    for (Gtk::TreeModel::iterator i = dirTreeModel->children().begin(); i != dirTreeModel->children().end(); ++i) {
        updateDirTree (i);
    }
}

void DirBrowser::updateDirTree (const Gtk::TreeModel::iterator& iter)
{

    if (dirtree->row_expanded (dirTreeModel->get_path (iter))) {
        updateDir (iter);

        for (Gtk::TreeModel::iterator i = iter->children().begin(); i != iter->children().end(); ++i) {
            updateDirTree (i);
        }
    }
}

void DirBrowser::updateVolumes ()
{

    unsigned int nvolumes = GetLogicalDrives ();

    if (nvolumes != volumes) {
        GThreadLock lock;

        for (int i = 0; i < 32; i++)
            if (((volumes >> i) & 1) && !((nvolumes >> i) & 1)) { // volume i has been deleted
                for (Gtk::TreeModel::iterator iter = dirTreeModel->children().begin(); iter != dirTreeModel->children().end(); ++iter)
                    if (iter->get_value (dtColumns.filename).c_str()[0] - 'A' == i) {
                        dirTreeModel->erase (iter);
                        break;
                    }
            } else if (!((volumes >> i) & 1) && ((nvolumes >> i) & 1)) {
                addRoot ('A' + i);    // volume i has been added
            }

        volumes = nvolumes;
    }
}

int updateVolumesUI (void* br)
{
    (static_cast<DirBrowser*>(br))->updateVolumes ();
    return 1;
}

#endif

void DirBrowser::fillRoot ()
{

#ifdef _WIN32
    volumes = GetLogicalDrives ();

    for (int i = 0; i < 32; i++)
        if ((volumes >> i) & 1) {
            addRoot ('A' + i);
        }

    // since sigc++ is not thread safe, we have to use the glib function
    g_timeout_add (5000, updateVolumesUI, this);
#else
    Gtk::TreeModel::Row rootRow = *(dirTreeModel->append());
    rootRow[dtColumns.filename] = "/";
    rootRow[dtColumns.dirname] = "/";
    Gtk::TreeModel::Row childRow = *(dirTreeModel->append(rootRow.children()));
    childRow[dtColumns.filename] = "foo";
#endif
}

void DirBrowser::on_sort_column_changed() const
{
    options.dirBrowserSortType = tvc.get_sort_order();
}

void DirBrowser::row_expanded (const Gtk::TreeModel::iterator& iter, const Gtk::TreeModel::Path& path)
{

    expandSuccess = false;

    // We will disable model's sorting because it decreases speed of inserting new items
    // in list tree dramatically. Therefore will do:
    // 1) Disable sorting in model
    // 2) Manually sort data in the order determined by the options
    // 3) Enable sorting in model again for UI (sorting by click on header)
    int prevSortColumn;
    Gtk::SortType prevSortType;
    dirTreeModel->get_sort_column_id(prevSortColumn, prevSortType);
    dirTreeModel->set_sort_column(Gtk::TreeSortable::DEFAULT_UNSORTED_COLUMN_ID, Gtk::SORT_ASCENDING);

    auto dir = Gio::File::create_for_path (iter->get_value (dtColumns.dirname));
    auto subDirs = listSubDirs (dir, options.fbShowHidden);

    Gtk::TreeNodeChildren children = iter->children();
    std::list<Gtk::TreeIter> forErase(children.begin(), children.end());

    std::sort (subDirs.begin (), subDirs.end (), [] (const Glib::ustring& firstDir, const Glib::ustring& secondDir)
    {
        switch (options.dirBrowserSortType) {
        default:
        case Gtk::SORT_ASCENDING:
            return firstDir < secondDir;
        case Gtk::SORT_DESCENDING:
            return firstDir > secondDir;
        }
    });

    for (auto it = subDirs.begin(), end = subDirs.end(); it != end; ++it) {
        addDir(iter, *it);
    }

    for (auto it = forErase.begin(), end = forErase.end(); it != end; ++it) {
        dirTreeModel->erase(*it);
    }

    dirTreeModel->set_sort_column(prevSortColumn, prevSortType);

    expandSuccess = true;

    Glib::RefPtr<Gio::FileMonitor> monitor = dir->monitor_directory ();
    iter->set_value (dtColumns.monitor, monitor);
    monitor->signal_changed().connect (sigc::bind(sigc::mem_fun(*this, &DirBrowser::file_changed), iter, dir->get_parse_name()));
}

void DirBrowser::updateDir (const Gtk::TreeModel::iterator& iter)
{

    // first test if some files are deleted
    bool change = true;

    while (change) {
        change = false;

        for (Gtk::TreeModel::iterator it = iter->children().begin(); it != iter->children().end(); ++it)
            if (!Glib::file_test (it->get_value (dtColumns.dirname), Glib::FILE_TEST_EXISTS)
                    || !Glib::file_test (it->get_value (dtColumns.dirname), Glib::FILE_TEST_IS_DIR)) {
                GThreadLock lock;
                dirTreeModel->erase (it);
                change = true;
                break;
            }
    }

    // test if new files are created
    auto dir = Gio::File::create_for_path (iter->get_value (dtColumns.dirname));
    auto subDirs = listSubDirs (dir, options.fbShowHidden);

    for (size_t i = 0; i < subDirs.size(); i++) {
        bool found = false;

        for (Gtk::TreeModel::iterator it = iter->children().begin(); it != iter->children().end() && !found ; ++it) {
            found = (it->get_value (dtColumns.filename) == subDirs[i]);
        }

        if (!found) {
            GThreadLock lock;
            addDir (iter, subDirs[i]);
        }
    }
}

void DirBrowser::addDir (const Gtk::TreeModel::iterator& iter, const Glib::ustring& dirname)
{

    Gtk::TreeModel::iterator child = dirTreeModel->append(iter->children());
    child->set_value (dtColumns.filename, dirname);
    child->set_value (dtColumns.icon1, openfolder);
    child->set_value (dtColumns.icon2, closedfolder);
    Glib::ustring fullname = Glib::build_filename (iter->get_value (dtColumns.dirname), dirname);
    child->set_value (dtColumns.dirname, fullname);
    Gtk::TreeModel::iterator fooRow = dirTreeModel->append(child->children());
    fooRow->set_value (dtColumns.filename, Glib::ustring("foo"));
}

void DirBrowser::row_activated (const Gtk::TreeModel::Path& path, Gtk::TreeViewColumn* column)
{

    Glib::ustring dname = dirTreeModel->get_iter (path)->get_value (dtColumns.dirname);

    if (Glib::file_test (dname, Glib::FILE_TEST_IS_DIR)) {
        dirSelectionSignal (dname, Glib::ustring());
        dirtree->expand_row(path, false);
    }
}

Gtk::TreePath DirBrowser::expandToDir (const Glib::ustring& absDirPath)
{

    Gtk::TreeModel::Path path;
    path.push_back(0);

    char* dcpy = strdup (absDirPath.c_str());
    char* dir = strtok (dcpy, "/\\");
    int count = 0;
    expandSuccess = true;

#ifndef _WIN32
    Gtk::TreeModel::iterator j = dirTreeModel->get_iter (path);
    path.up ();
    path.push_back (0);
    row_expanded(j, path);
    path.push_back (0);
#endif

    while (dir) {
        Glib::ustring dirstr = dir;
#ifdef _WIN32

        if (count == 0) {
            dirstr = dirstr + "\\";
        }

#endif
        Gtk::TreeModel::iterator i = dirTreeModel->get_iter (path);
        int ix = 0;

        while (i && expandSuccess) {
            Gtk::TreeModel::Row crow = *i;
            Glib::ustring str = crow[dtColumns.filename];
#ifdef _WIN32

            if (str.casefold() == dirstr.casefold()) {
#else

            if (str == dirstr) {
#endif
                path.up ();
                path.push_back (ix);
                row_expanded(i, path);
                path.push_back (0);
                break;
            }

            ++ix;
            ++i;
        }

        count++;
        dir = strtok(nullptr, "/\\");
    }

    free(dcpy);

    path.up ();
    dirtree->expand_to_path (path);

    return path;
}

void DirBrowser::open (const Glib::ustring& dirname, const Glib::ustring& fileName)
{

    dirtree->collapse_all ();

    // WARNING & TODO: One should test here if the directory/file has R/W access permission to avoid crash

    Glib::RefPtr<Gio::File> dir = Gio::File::create_for_path(dirname);

    if( !dir->query_exists()) {
        return;
    }

    Glib::ustring absDirPath = dir->get_parse_name ();
    Gtk::TreePath path = expandToDir (absDirPath);
    dirtree->scroll_to_row (path);
    dirtree->get_selection()->select (path);
    Glib::ustring absFilePath;

    if (!fileName.empty()) {
        absFilePath = Glib::build_filename (absDirPath, fileName);
    }

    dirSelectionSignal (absDirPath, absFilePath);
}

void DirBrowser::file_changed (const Glib::RefPtr<Gio::File>& file, const Glib::RefPtr<Gio::File>& other_file, Gio::FileMonitorEvent event_type, const Gtk::TreeModel::iterator& iter, const Glib::ustring& dirName)
{

    if (!file || !Glib::file_test (dirName, Glib::FILE_TEST_IS_DIR) || event_type == Gio::FILE_MONITOR_EVENT_ATTRIBUTE_CHANGED) {
        return;
    }

    updateDir (iter);
}

void DirBrowser::selectDir (Glib::ustring dir)
{

    open (dir, "");
}

