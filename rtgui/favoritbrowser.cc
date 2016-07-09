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
#include <favoritbrowser.h>
#include "multilangmgr.h"
#include "rtimage.h"

FavoritBrowser::FavoritBrowser () : listener (NULL), lastSelectedDir ("")
{

    scrollw = Gtk::manage (new Gtk::ScrolledWindow ());
    scrollw->set_policy (Gtk::POLICY_NEVER, Gtk::POLICY_AUTOMATIC);

    Gtk::Frame* frame = Gtk::manage (new Gtk::Frame ("Favorite Folders"));
    frame->add (*scrollw);

    pack_start (*frame);

    treeView = Gtk::manage (new Gtk::TreeView ());
    scrollw->add (*treeView);

    favoritModel = Gtk::ListStore::create (favoritColumns);
    treeView->set_model (favoritModel);
    treeView->set_headers_visible (false);

    Gtk::TreeView::Column *iviewcol = Gtk::manage (new Gtk::TreeView::Column ("icon"));
    Gtk::CellRendererPixbuf *iconCR  = Gtk::manage (new Gtk::CellRendererPixbuf());
    iviewcol->pack_start (*iconCR, false);
    iviewcol->add_attribute (*iconCR, "gicon", 0);

    treeView->append_column (*iviewcol);
    treeView->append_column ("text", favoritColumns.shortdir);

    treeView->set_tooltip_column (2);
    treeView->get_selection()->signal_changed().connect(sigc::mem_fun(*this, &FavoritBrowser::selectionChanged));

    add = Gtk::manage (new Gtk::Button ());
    add->set_tooltip_text(M("MAIN_FRAME_PLACES_ADD"));
    add->set_hexpand(true);
    add->set_vexpand(false);
    add->set_halign(Gtk::ALIGN_FILL);
    add->set_valign(Gtk::ALIGN_START);
    add->set_image (*Gtk::manage (new RTImage ("gtk-add.png")));
    add->get_style_context()->add_class("Left");
    del = Gtk::manage (new Gtk::Button ());
    del->set_tooltip_text(M("MAIN_FRAME_PLACES_DEL"));
    del->set_hexpand(true);
    del->set_vexpand(false);
    del->set_halign(Gtk::ALIGN_FILL);
    del->set_valign(Gtk::ALIGN_START);
    del->set_image (*Gtk::manage (new RTImage ("list-remove.png")));
    del->get_style_context()->add_class("Right");
    Gtk::HBox* buttonBox = Gtk::manage (new Gtk::HBox ());
    buttonBox->pack_start (*add);
    buttonBox->pack_start (*del);

    pack_start (*buttonBox, Gtk::PACK_SHRINK, 2);

    add->signal_clicked().connect(sigc::mem_fun(*this, &FavoritBrowser::addPressed));
    del->signal_clicked().connect(sigc::mem_fun(*this, &FavoritBrowser::delPressed));

    show_all ();
}

void FavoritBrowser::selectionChanged ()
{

    Glib::RefPtr<Gtk::TreeSelection> selection = treeView->get_selection();
    Gtk::TreeModel::iterator iter = selection->get_selected();

    if (iter && listener) {
        listener->selectDir (iter->get_value (favoritColumns.fulldir));
    }
}

void FavoritBrowser::dirSelected (const Glib::ustring& dirname, const Glib::ustring& openfile)
{

    lastSelectedDir = dirname;
}

void FavoritBrowser::addPressed ()
{

    if (lastSelectedDir == "") {
        return;
    }

    // check if the dirname is already in the list. If yes, return.
    Gtk::TreeModel::iterator iter = favoritModel->children ().begin();

    while (iter != favoritModel->children().end()) {
        if (iter->get_value (favoritColumns.fulldir) == lastSelectedDir) {
            return;
        }

        ++iter;
    }

    Glib::RefPtr<Gio::File> hfile = Gio::File::create_for_parse_name (lastSelectedDir);

    if (hfile) {
        Glib::RefPtr<Gio::FileInfo> info = hfile->query_info ();

        if (info) {
            Gtk::TreeModel::Row newrow = *(favoritModel->append());
            newrow[favoritColumns.shortdir] = info->get_display_name ();
            newrow[favoritColumns.fulldir] = lastSelectedDir;
            newrow[favoritColumns.icon] = info->get_icon ();
        }
    }
}

void FavoritBrowser::delPressed ()
{

    // lookup the selected item in the bookmark
    Glib::RefPtr<Gtk::TreeSelection> selection = treeView->get_selection();
    Gtk::TreeModel::iterator iter = selection->get_selected();

    if (iter) {
        favoritModel->erase (iter);
    }
}

