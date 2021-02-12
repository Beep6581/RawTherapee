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
#include "renamedlg.h"
#include "cacheimagedata.h"
#include "multilangmgr.h"
#include "rtimage.h"

RenameDialog::RenameDialog (Gtk::Window* parent)
    : Gtk::Dialog (M("FILEBROWSER_RENAMEDLGLABEL"), *parent, true), p(parent), imageData(nullptr)
{

    Gtk::Table* names = Gtk::manage (new Gtk::Table (2, 2));
    Gtk::Label* onlab = Gtk::manage (new Gtk::Label (M("FILEBROWSER_CURRENT_NAME")));
    Gtk::Label* nnlab = Gtk::manage (new Gtk::Label (M("FILEBROWSER_NEW_NAME")));
    oldName = Gtk::manage (new Gtk::Label ("alma"));
    newName = Gtk::manage (new Gtk::Entry ());

    names->attach (*onlab, 0, 1, 0, 1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    names->attach (*oldName, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);
    names->attach (*nnlab, 0, 1, 1, 2, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    names->attach (*newName, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    get_content_area()->pack_start (*names, Gtk::PACK_SHRINK, 4);

// Issue 316
//    Gtk::Box* tbox = Gtk::manage (new Gtk::Box());
//    useTmpl = Gtk::manage (new Gtk::CheckButton (M("FILEBROWSER_USETEMPLATE")));
//    templates = Gtk::manage (new MyComboBox ());
//    templateModel = Gtk::ListStore::create (templateColumns);
//    templates->set_model (templateModel);
//    templates->pack_start (templateColumns.tmplName);

//    tbox->pack_start (*useTmpl, Gtk::PACK_SHRINK, 4);
//    tbox->pack_start (*templates);

//    get_content_area()->pack_start (*tbox, Gtk::PACK_SHRINK, 4);

    add_button (Gtk::Stock::OK, Gtk::RESPONSE_OK);
    add_button (Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
// Issue 316
//    all = add_button ("All", RESPONSE_ALL);

    newName->set_activates_default (true);
    set_default_response (Gtk::RESPONSE_OK);
// Issue 316
//    fillTemplateList ();

//    templates->set_row_separator_func (sigc::mem_fun(*this, &RenameDialog::rowSeparatorFunc));
//    templates->signal_changed().connect(sigc::mem_fun(*this, &RenameDialog::tmplSelectionChanged));
//    useTmpl->signal_toggled().connect( sigc::mem_fun(*this, &RenameDialog::useTemplToggled) );

//    useTmpl->set_active (options.renameUseTemplates);

    show_all_children ();
}

void RenameDialog::initName (const Glib::ustring& iname, const CacheImageData* cid)
{

    imageData = cid;
    oldName->set_text (iname);
    newName->set_text (iname);
// Issue 316
//    if (useTmpl->get_active () && isTemplSelected ())
//        newName->set_text (applyTemplate (iname, cid, getActiveTemplate()));
    newName->select_region (0, newName->get_text().size());
}

Glib::ustring RenameDialog::getNewName ()
{

    return newName->get_text ();
}

