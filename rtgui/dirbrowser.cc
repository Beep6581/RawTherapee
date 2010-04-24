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
#include <dirbrowser.h>
#ifdef WIN32
#define _WIN32_WINNT 0x0600
#include <windows.h>
#endif
#include <options.h>
#include <safegtk.h>

#define CHECKTIME 5000
extern Glib::ustring argv0;

DirBrowser::DirBrowser () {

   dirtree = new Gtk::TreeView();
   scrolledwindow4 = new Gtk::ScrolledWindow();

//   dirtree->set_flags(Gtk::CAN_FOCUS);
   dirtree->set_headers_visible(false);
   dirtree->set_rules_hint(false);
   dirtree->set_reorderable(false);
   dirtree->set_enable_search(false);
   scrolledwindow4->set_flags(Gtk::CAN_FOCUS);
   scrolledwindow4->set_border_width(2);
   scrolledwindow4->set_shadow_type(Gtk::SHADOW_NONE);
   scrolledwindow4->set_policy(Gtk::POLICY_ALWAYS, Gtk::POLICY_ALWAYS);
   scrolledwindow4->property_window_placement().set_value(Gtk::CORNER_TOP_LEFT);
   scrolledwindow4->add(*dirtree);

   pack_start (*scrolledwindow4);
   dirtree->show ();
   scrolledwindow4->show ();
}

void DirBrowser::fillDirTree () {

  openfolder = safe_create_from_file (argv0+"/images/folder_open.png");
  closedfolder = safe_create_from_file (argv0+"/images/folder.png");
  icdrom = safe_create_from_file (argv0+"/images/cdrom.png");
  ifloppy = safe_create_from_file (argv0+"/images/floppy.png");
  ihdd = safe_create_from_file (argv0+"/images/hdd.png");
  iremovable = safe_create_from_file (argv0+"/images/usbpendrive.png");
  inetwork = safe_create_from_file (argv0+"/images/network.png");

  //Create the Tree model:
  dirTreeModel = Gtk::TreeStore::create(dtColumns);
  dirtree->set_model (dirTreeModel);

  fillRoot ();

  Gtk::CellRendererPixbuf* render_pb = new Gtk::CellRendererPixbuf ();
  tvc.pack_start (*render_pb, false);
  tvc.add_attribute(*render_pb, "pixbuf-expander-closed", 1);
  tvc.add_attribute(*render_pb, "pixbuf", 1);
  tvc.add_attribute(*render_pb, "pixbuf-expander-open", 0);
  tvc.pack_start (crt);
  tvc.add_attribute(crt, "text", 2);
 
  crt.property_ypad() = 0;
  render_pb->property_ypad() = 0;
  
  dirtree->append_column(tvc); 

  dirtree->signal_row_expanded().connect(sigc::mem_fun(*this, &DirBrowser::row_expanded));
  dirtree->signal_row_activated().connect(sigc::mem_fun(*this, &DirBrowser::row_activated));
}

#ifdef WIN32
void DirBrowser::addRoot (char letter) {

    char volume[4];
    volume[0] = letter;
    strcpy (volume+1, ":\\");

    Gtk::TreeModel::iterator root = dirTreeModel->append();
    root->set_value (dtColumns.filename, Glib::ustring(volume));
    root->set_value (dtColumns.dirname, Glib::ustring(volume));

    int type = GetDriveType (volume);
    if (type==DRIVE_CDROM) {
        root->set_value (0, icdrom);
        root->set_value (1, icdrom);
    }
    else if (type==DRIVE_REMOVABLE) {
        if (letter-'A'<2) {
            root->set_value (0, ifloppy);
            root->set_value (1, ifloppy);
        }
        else {
            root->set_value (0, iremovable);
            root->set_value (1, iremovable);
        }
    }
    else if (type==DRIVE_REMOTE) {
        root->set_value (0, inetwork);
        root->set_value (1, inetwork);
    }
    else if (type==DRIVE_FIXED) {
        root->set_value (0, ihdd);
        root->set_value (1, ihdd);
    }

    Gtk::TreeModel::iterator child = dirTreeModel->append (root->children());
    child->set_value (dtColumns.filename, Glib::ustring("foo"));
}

void DirBrowser::updateDirTreeRoot () {

    for (Gtk::TreeModel::iterator i=dirTreeModel->children().begin(); i!=dirTreeModel->children().end(); i++)
		updateDirTree (i);
}

void DirBrowser::updateDirTree (const Gtk::TreeModel::iterator& iter) {
	
	if (dirtree->row_expanded (dirTreeModel->get_path (iter))) {
		updateDir (iter);
		for (Gtk::TreeModel::iterator i=iter->children().begin(); i!=iter->children().end(); i++)
			updateDirTree (i);
	}
}

void DirBrowser::updateVolumes () {

    int nvolumes = GetLogicalDrives ();
    if (nvolumes!=volumes) {
        for (int i=0; i<32; i++) 
            if (((volumes >> i) & 1) && !((nvolumes >> i) & 1)) { // volume i has been deleted
                for (Gtk::TreeModel::iterator iter = dirTreeModel->children().begin(); iter!=dirTreeModel->children().end(); iter++) 
                    if (iter->get_value (dtColumns.filename).c_str()[0]-'A' == i) {
                        dirTreeModel->erase (iter);
                        break;
                    }
            }
            else if (!((volumes >> i) & 1) && ((nvolumes >> i) & 1)) 
                addRoot ('A'+i); // volume i has been added
        volumes = nvolumes;
    }
}

int _updateVolumes (void* br) {

    gdk_threads_enter ();
    ((DirBrowser*)br)->updateVolumes ();
    gdk_threads_leave ();
    return 1;
}
int _updateDirTree (void* br) {

    gdk_threads_enter ();
    ((DirBrowser*)br)->updateDirTreeRoot ();
    gdk_threads_leave ();
    return 0;
}

void DirBrowser::winDirChanged () {

    g_idle_add (_updateDirTree, this);
}
#endif

void DirBrowser::fillRoot () {

#ifdef WIN32
  volumes = GetLogicalDrives ();
  for (int i=0; i<32; i++)
    if ((volumes >> i) & 1) 
        addRoot ('A'+i);
  // since sigc++ is not thread safe, we have to use the glib function
  g_timeout_add (CHECKTIME, _updateVolumes, this);
#else
  Gtk::TreeModel::Row rootRow = *(dirTreeModel->append());
  rootRow[dtColumns.filename] = "/";
  rootRow[dtColumns.dirname] = "/";
  Gtk::TreeModel::Row childRow = *(dirTreeModel->append(rootRow.children()));
  childRow[dtColumns.filename] = "foo";
#endif
}

void DirBrowser::row_expanded (const Gtk::TreeModel::iterator& iter, const Gtk::TreeModel::Path& path) {

  expandSuccess = false;

  int todel = iter->children().size();

	std::vector<Glib::ustring> subDirs;
  Glib::RefPtr<Gio::File> dir = Gio::File::create_for_path (iter->get_value (dtColumns.dirname));
    
  safe_build_subdir_list (dir, subDirs, options.fbShowHidden);

	if (subDirs.size() == 0)
			dirtree->collapse_row (path);
	else {
	
			std::sort (subDirs.begin(), subDirs.end());
			for (int i=0; i<subDirs.size(); i++) 
					addDir (iter, subDirs[i]);

			for (int i=0; i<todel; i++)
					dirTreeModel->erase (iter->children().begin());
			expandSuccess = true;
	}
#ifdef _WIN32
  Glib::RefPtr<WinDirMonitor> monitor = Glib::RefPtr<WinDirMonitor>(new WinDirMonitor (iter->get_value (dtColumns.dirname), this));
  iter->set_value (dtColumns.monitor, monitor);
#elif defined __APPLE__
  printf("TODO fix dir->monitor_directory () for OSX\n"); 
#else
  Glib::RefPtr<Gio::FileMonitor> monitor = dir->monitor_directory ();
  iter->set_value (dtColumns.monitor, monitor);
  monitor->signal_changed().connect (sigc::bind(sigc::mem_fun(*this, &DirBrowser::file_changed), iter, dir->get_parse_name()));
#endif
}

void DirBrowser::updateDir (const Gtk::TreeModel::iterator& iter) {

    // first test if some files are deleted
    bool change = true;
    while (change) {
        change = false;
        for (Gtk::TreeModel::iterator it=iter->children().begin(); it!=iter->children().end(); it++)
            if (!Glib::file_test (it->get_value (dtColumns.dirname), Glib::FILE_TEST_EXISTS) 
             || !Glib::file_test (it->get_value (dtColumns.dirname), Glib::FILE_TEST_IS_DIR)) {
                dirTreeModel->erase (it);
                change = true;
                break;
            }
    }
    // test if new files are created
		std::vector<Glib::ustring> subDirs;
    Glib::RefPtr<Gio::File> dir = Gio::File::create_for_path (iter->get_value (dtColumns.dirname));
		safe_build_subdir_list (dir, subDirs, options.fbShowHidden);

    for (int i=0; i<subDirs.size(); i++) {
        bool found = false;
        for (Gtk::TreeModel::iterator it=iter->children().begin(); it!=iter->children().end() && !found ; it++) 
            found = (it->get_value (dtColumns.filename)==subDirs[i]);

        if (!found)
            addDir (iter, subDirs[i]);
    }
}

void DirBrowser::addDir (const Gtk::TreeModel::iterator& iter, const Glib::ustring& dirname) {

    Gtk::TreeModel::iterator child = dirTreeModel->append(iter->children());
    child->set_value (dtColumns.filename, dirname);
    child->set_value (0, openfolder);
    child->set_value (1, closedfolder);
    Glib::ustring fullname = Glib::build_filename (iter->get_value (dtColumns.dirname), dirname);
    child->set_value (dtColumns.dirname, fullname);
    Glib::RefPtr<Gio::File> f = Gio::File::create_for_path (fullname);
    Gtk::TreeModel::iterator fooRow = dirTreeModel->append(child->children());
    fooRow->set_value (dtColumns.filename, Glib::ustring("foo"));
}

void DirBrowser::row_activated (const Gtk::TreeModel::Path& path, Gtk::TreeViewColumn* column) {

    Glib::ustring dname = dirTreeModel->get_iter (path)->get_value (dtColumns.dirname);
    if (Glib::file_test (dname, Glib::FILE_TEST_IS_DIR)) 
        for (int i=0; i<dllisteners.size(); i++)
            dllisteners[i]->dirSelected (dname);
}

Gtk::TreePath DirBrowser::expandToDir (const Glib::ustring& absDirPath) {

    Gtk::TreeModel::Path path;
    path.append_index(0);

    int end = 0;
    int beg = 0;
    char* dcpy = strdup (absDirPath.c_str());
    char* dir = strtok (dcpy, "/\\");
    int count = 0;
    expandSuccess = true;

#ifndef _WIN32
    Gtk::TreeModel::iterator j = dirTreeModel->get_iter (path);
    path.up ();
    path.append_index (0);
    row_expanded(j, path);
    path.append_index (0);
#endif

    while (dir) {
        Glib::ustring dirstr = dir;
#ifdef _WIN32
        if (count==0)
            dirstr = dirstr + "\\";
#endif
        Gtk::TreeModel::iterator i = dirTreeModel->get_iter (path);
        int ix = 0;
        while (i && expandSuccess) {
            Gtk::TreeModel::Row crow = *i;
            Glib::ustring str =crow[dtColumns.filename]; 
#ifdef _WIN32
            if (str.casefold()==dirstr.casefold()) {
#else
            if (str==dirstr) {
#endif
                path.up ();
                path.append_index (ix);
                row_expanded(i, path);
                path.append_index (0);
                break;
            }
            ix++;
            i++;
        }
        count++;
        dir = strtok(NULL, "/\\");
    }

    free(dcpy);

    path.up ();
    dirtree->expand_to_path (path);

    return path;
}

void DirBrowser::open (const Glib::ustring& dirname, const Glib::ustring& fileName) {
  
    dirtree->collapse_all ();

    Glib::ustring absDirPath = Gio::File::create_for_path(dirname)->get_parse_name ();
    Gtk::TreePath path = expandToDir (absDirPath);

    if (expandSuccess) {
        dirtree->scroll_to_row (path);
        dirtree->get_selection()->select (path);
        for (int i=0; i<dllisteners.size(); i++)
            dllisteners[i]->dirSelected (absDirPath, Glib::build_filename (absDirPath, fileName));
    }
}

void DirBrowser::file_changed (const Glib::RefPtr<Gio::File>& file, const Glib::RefPtr<Gio::File>& other_file, Gio::FileMonitorEvent event_type, const Gtk::TreeModel::iterator& iter, const Glib::ustring& dirName) {

    if (!file || !Glib::file_test (dirName, Glib::FILE_TEST_IS_DIR) || event_type==Gio::FILE_MONITOR_EVENT_ATTRIBUTE_CHANGED) 
        return;

    gdk_threads_enter();
    updateDir (iter);
    gdk_threads_leave();
}

void DirBrowser::selectDir (Glib::ustring dir) {

    open (dir, "");
}

