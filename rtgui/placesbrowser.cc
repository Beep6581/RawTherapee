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
#include "placesbrowser.h"
#include "options.h"
#include "toolpanel.h"
#include "../rtengine/safegtk.h"
#include "guiutils.h"
#include "rtimage.h"

PlacesBrowser::PlacesBrowser () : listener (NULL)
{

    scrollw = Gtk::manage (new Gtk::ScrolledWindow ());
    scrollw->set_policy (Gtk::POLICY_NEVER, Gtk::POLICY_AUTOMATIC);
    pack_start (*scrollw);

    // Since Gtk3, we can't have image+text buttons natively. We'll comply to the Gtk guidelines and choose one of them (icons here)
    add = Gtk::manage (new Gtk::Button ());
    add->set_tooltip_text(M("MAIN_FRAME_PLACES_ADD"));
    setExpandAlignProperties(add, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);
    //add->get_style_context()->set_junction_sides(Gtk::JUNCTION_RIGHT);
    add->get_style_context()->add_class("Left");
    add->set_image (*Gtk::manage (new RTImage ("list-add.png")));
    del = Gtk::manage (new Gtk::Button ());
    del->set_tooltip_text(M("MAIN_FRAME_PLACES_DEL"));
    setExpandAlignProperties(del, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);
    //del->get_style_context()->set_junction_sides(Gtk::JUNCTION_LEFT);
    del->get_style_context()->add_class("Right");
    del->set_image (*Gtk::manage (new RTImage ("list-remove.png")));
    Gtk::Grid* buttonBox = Gtk::manage (new Gtk::Grid ());
    buttonBox->set_orientation(Gtk::ORIENTATION_HORIZONTAL);
    buttonBox->attach_next_to(*add, Gtk::POS_LEFT, 1, 1);
    buttonBox->attach_next_to(*del, *add, Gtk::POS_RIGHT, 1, 1);

    pack_start (*buttonBox, Gtk::PACK_SHRINK, 2);

    treeView = Gtk::manage (new Gtk::TreeView ());
    treeView->set_can_focus(false);
    scrollw->add (*treeView);

    placesModel = Gtk::ListStore::create (placesColumns);
    treeView->set_model (placesModel);
    treeView->set_headers_visible (true);

    Gtk::TreeView::Column *iviewcol = Gtk::manage (new Gtk::TreeView::Column (M("MAIN_FRAME_PLACES")));
    Gtk::CellRendererPixbuf *iconCR  = Gtk::manage (new Gtk::CellRendererPixbuf());
    Gtk::CellRendererText *labelCR  = Gtk::manage (new Gtk::CellRendererText());
    iviewcol->pack_start (*iconCR, false);
    iviewcol->pack_start (*labelCR, true);
    iviewcol->add_attribute (*iconCR, "gicon", 0);
    iviewcol->add_attribute (*labelCR, "text", placesColumns.label);
    treeView->append_column (*iviewcol);

    treeView->set_row_separator_func (sigc::mem_fun(*this, &PlacesBrowser::rowSeparatorFunc));

    vm = Gio::VolumeMonitor::get();

    vm->signal_mount_changed().connect (sigc::mem_fun(*this, &PlacesBrowser::mountChanged));
    vm->signal_mount_added().connect (sigc::mem_fun(*this, &PlacesBrowser::mountChanged));
    vm->signal_mount_removed().connect (sigc::mem_fun(*this, &PlacesBrowser::mountChanged));
    vm->signal_volume_changed().connect (sigc::mem_fun(*this, &PlacesBrowser::volumeChanged));
    vm->signal_volume_added().connect (sigc::mem_fun(*this, &PlacesBrowser::volumeChanged));
    vm->signal_volume_removed().connect (sigc::mem_fun(*this, &PlacesBrowser::volumeChanged));
    vm->signal_drive_connected().connect (sigc::mem_fun(*this, &PlacesBrowser::driveChanged));
    vm->signal_drive_disconnected().connect (sigc::mem_fun(*this, &PlacesBrowser::driveChanged));
    vm->signal_drive_changed().connect (sigc::mem_fun(*this, &PlacesBrowser::driveChanged));

    treeView->get_selection()->signal_changed().connect(sigc::mem_fun(*this, &PlacesBrowser::selectionChanged));
    add->signal_clicked().connect(sigc::mem_fun(*this, &PlacesBrowser::addPressed));
    del->signal_clicked().connect(sigc::mem_fun(*this, &PlacesBrowser::delPressed));

    show_all ();
}

// For drive letter comparison
bool compareMountByRoot (Glib::RefPtr<Gio::Mount> a, Glib::RefPtr<Gio::Mount> b)
{
    return a->get_root()->get_parse_name() < b->get_root()->get_parse_name();
}

void PlacesBrowser::refreshPlacesList ()
{

    placesModel->clear ();

    // append home directory
    Glib::RefPtr<Gio::File> hfile = Gio::File::create_for_path (safe_get_user_home_dir());  // Will send back "My documents" on Windows now, which has no restricted access

    if (hfile && hfile->query_exists()) {
        try {
            Glib::RefPtr<Gio::FileInfo> info = safe_query_file_info (hfile);

            if (info) {
                Gtk::TreeModel::Row newrow = *(placesModel->append());
                newrow[placesColumns.label] = info->get_display_name ();
                newrow[placesColumns.icon]  = info->get_icon ();
                newrow[placesColumns.root]  = hfile->get_parse_name ();
                newrow[placesColumns.type]  = 4;
                newrow[placesColumns.rowSeparator] = false;
            }
        } catch (Gio::Error&) {
            /* This will be thrown if the path doesn't exist */
        }
    }

    // append pictures directory
    hfile = Gio::File::create_for_path (safe_get_user_picture_dir());

    if (hfile && hfile->query_exists()) {
        try {
            Glib::RefPtr<Gio::FileInfo> info = safe_query_file_info (hfile);

            if (info) {
                Gtk::TreeModel::Row newrow = *(placesModel->append());
                newrow[placesColumns.label] = info->get_display_name ();
                newrow[placesColumns.icon]  = info->get_icon ();
                newrow[placesColumns.root]  = hfile->get_parse_name ();
                newrow[placesColumns.type]  = 4;
                newrow[placesColumns.rowSeparator] = false;
            }
        } catch (Gio::Error&) {
            /* This will be thrown if the path doesn't exist */
        }
    }

    if (!placesModel->children().empty()) {
        Gtk::TreeModel::Row newrow = *(placesModel->append());
        newrow[placesColumns.rowSeparator] = true;
    }

    // scan all drives
    std::vector<Glib::RefPtr<Gio::Drive> > drives = vm->get_connected_drives ();

    for (size_t j = 0; j < drives.size (); j++) {
        std::vector<Glib::RefPtr<Gio::Volume> > volumes = drives[j]->get_volumes ();

        if (volumes.empty()) {
            Gtk::TreeModel::Row newrow = *(placesModel->append());
            newrow[placesColumns.label] = drives[j]->get_name ();
            newrow[placesColumns.icon]  = drives[j]->get_icon ();
            newrow[placesColumns.root]  = "";
            newrow[placesColumns.type]  = 3;
            newrow[placesColumns.rowSeparator] = false;
        }

        for (size_t i = 0; i < volumes.size (); i++) {
            Glib::RefPtr<Gio::Mount> mount = volumes[i]->get_mount ();

            if (mount) { // placesed volumes
                Gtk::TreeModel::Row newrow = *(placesModel->append());
                newrow[placesColumns.label] = mount->get_name ();
                newrow[placesColumns.icon]  = mount->get_icon ();
                newrow[placesColumns.root]  = mount->get_root ()->get_parse_name ();
                newrow[placesColumns.type]  = 1;
                newrow[placesColumns.rowSeparator] = false;
            } else { // unplacesed volumes
                Gtk::TreeModel::Row newrow = *(placesModel->append());
                newrow[placesColumns.label] = volumes[i]->get_name ();
                newrow[placesColumns.icon]  = volumes[i]->get_icon ();
                newrow[placesColumns.root]  = "";
                newrow[placesColumns.type]  = 2;
                newrow[placesColumns.rowSeparator] = false;
            }
        }
    }

    // volumes not belonging to drives
    std::vector<Glib::RefPtr<Gio::Volume> > volumes = vm->get_volumes ();

    for (size_t i = 0; i < volumes.size (); i++) {
        if (!volumes[i]->get_drive ()) {
            Glib::RefPtr<Gio::Mount> mount = volumes[i]->get_mount ();

            if (mount) { // placesed volumes
                Gtk::TreeModel::Row newrow = *(placesModel->append());
                newrow[placesColumns.label] = mount->get_name ();
                newrow[placesColumns.icon]  = mount->get_icon ();
                newrow[placesColumns.root]  = mount->get_root ()->get_parse_name ();
                newrow[placesColumns.type]  = 1;
                newrow[placesColumns.rowSeparator] = false;
            } else { // unplacesed volumes
                Gtk::TreeModel::Row newrow = *(placesModel->append());
                newrow[placesColumns.label] = volumes[i]->get_name ();
                newrow[placesColumns.icon]  = volumes[i]->get_icon ();
                newrow[placesColumns.root]  = "";
                newrow[placesColumns.type]  = 2;
                newrow[placesColumns.rowSeparator] = false;
            }
        }
    }

    // places not belonging to volumes
    // (Drives in Windows)
    std::vector<Glib::RefPtr<Gio::Mount> > mounts = vm->get_mounts ();

#ifdef WIN32
    // on Windows, it's usual to sort by drive letter, not by name
    std::sort (mounts.begin(), mounts.end(), compareMountByRoot);
#endif

    for (size_t i = 0; i < mounts.size (); i++) {
        if (!mounts[i]->get_volume ()) {
            Gtk::TreeModel::Row newrow = *(placesModel->append());
            newrow[placesColumns.label] = mounts[i]->get_name ();
            newrow[placesColumns.icon]  = mounts[i]->get_icon ();
            newrow[placesColumns.root]  = mounts[i]->get_root ()->get_parse_name ();
            newrow[placesColumns.type]  = 1;
            newrow[placesColumns.rowSeparator] = false;
        }
    }

    // append favorites
    if (!placesModel->children().empty()) {
        Gtk::TreeModel::Row newrow = *(placesModel->append());
        newrow[placesColumns.rowSeparator] = true;
    }

    for (size_t i = 0; i < options.favoriteDirs.size(); i++) {
        Glib::RefPtr<Gio::File> hfile = Gio::File::create_for_path (options.favoriteDirs[i]);

        if (hfile && hfile->query_exists()) {
            Glib::RefPtr<Gio::FileInfo> info = safe_query_file_info (hfile);

            if (info) {
                Gtk::TreeModel::Row newrow = *(placesModel->append());
                newrow[placesColumns.label] = info->get_display_name ();
                newrow[placesColumns.icon]  = info->get_icon ();
                newrow[placesColumns.root]  = hfile->get_parse_name ();
                newrow[placesColumns.type]  = 5;
                newrow[placesColumns.rowSeparator] = false;
            }
        }
    }
}

bool PlacesBrowser::rowSeparatorFunc (const Glib::RefPtr<Gtk::TreeModel>& model, const Gtk::TreeModel::iterator& iter)
{

    return iter->get_value (placesColumns.rowSeparator);
}

void PlacesBrowser::mountChanged (const Glib::RefPtr<Gio::Mount>& m)
{
    GThreadLock lock;
    refreshPlacesList ();
}

void PlacesBrowser::volumeChanged (const Glib::RefPtr<Gio::Volume>& m)
{
    GThreadLock lock;
    refreshPlacesList ();
}

void PlacesBrowser::driveChanged (const Glib::RefPtr<Gio::Drive>& m)
{
    GThreadLock lock;
    refreshPlacesList ();
}

void PlacesBrowser::selectionChanged ()
{

    Glib::RefPtr<Gtk::TreeSelection> selection = treeView->get_selection();
    Gtk::TreeModel::iterator iter = selection->get_selected();

    if (iter) {
        if (iter->get_value (placesColumns.type) == 2) {
            std::vector<Glib::RefPtr<Gio::Volume> > volumes = vm->get_volumes ();

            for (size_t i = 0; i < volumes.size(); i++)
                if (volumes[i]->get_name () == iter->get_value (placesColumns.label)) {
                    volumes[i]->mount ();
                    break;
                }
        } else if (iter->get_value (placesColumns.type) == 3) {
            std::vector<Glib::RefPtr<Gio::Drive> > drives = vm->get_connected_drives ();

            for (size_t i = 0; i < drives.size(); i++)
                if (drives[i]->get_name () == iter->get_value (placesColumns.label)) {
                    drives[i]->poll_for_media ();
                    break;
                }
        } else if (listener) {
            listener->selectDir (iter->get_value (placesColumns.root));
        }
    }
}

void PlacesBrowser::dirSelected (const Glib::ustring& dirname, const Glib::ustring& openfile)
{

    lastSelectedDir = dirname;
}

void PlacesBrowser::addPressed ()
{

    if (lastSelectedDir == "") {
        return;
    }

    // check if the dirname is already in the list. If yes, return.
    for (size_t i = 0; i < options.favoriteDirs.size(); i++)
        if (options.favoriteDirs[i] == lastSelectedDir) {
            return;
        }

    // append
    Glib::RefPtr<Gio::File> hfile = Gio::File::create_for_path (lastSelectedDir);

    if (hfile && hfile->query_exists()) {
        Glib::RefPtr<Gio::FileInfo> info = safe_query_file_info (hfile);

        if (info) {
            options.favoriteDirs.push_back (hfile->get_parse_name ());
            refreshPlacesList ();
        }
    }
}

void PlacesBrowser::delPressed ()
{

    // lookup the selected item in the bookmark
    Glib::RefPtr<Gtk::TreeSelection> selection = treeView->get_selection();
    Gtk::TreeModel::iterator iter = selection->get_selected();

    if (iter && iter->get_value (placesColumns.type) == 5) {
        std::vector<Glib::ustring>::iterator i = std::find (options.favoriteDirs.begin(), options.favoriteDirs.end(), iter->get_value (placesColumns.root));

        if (i != options.favoriteDirs.end()) {
            options.favoriteDirs.erase (i);
        }
    }

    refreshPlacesList ();
}

