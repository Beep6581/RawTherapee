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
#include "renamedlg.h"
#include "multilangmgr.h"
#include "options.h"
#include "rtimage.h"

RenameDialog::RenameDialog (Gtk::Window* parent) 
    : Gtk::Dialog (M("FILEBROWSER_RENAMEDLGLABEL"), *parent, true, true), p(parent), imageData(NULL) {
    
    Gtk::Table* names = Gtk::manage (new Gtk::Table (2, 2));
    Gtk::Label* onlab = Gtk::manage (new Gtk::Label (M("FILEBROWSER_CURRENT_NAME")));
      Gtk::Label* nnlab = Gtk::manage (new Gtk::Label (M("FILEBROWSER_NEW_NAME")));
    oldName = Gtk::manage (new Gtk::Label ("alma"));
    newName = Gtk::manage (new Gtk::Entry ());
    
    names->attach (*onlab, 0, 1, 0, 1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    names->attach (*oldName, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);
    names->attach (*nnlab, 0, 1, 1, 2, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    names->attach (*newName, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    get_vbox()->pack_start (*names, Gtk::PACK_SHRINK, 4);

// Issue 316
//    Gtk::HBox* tbox = Gtk::manage (new Gtk::HBox());
//    useTmpl = Gtk::manage (new Gtk::CheckButton (M("FILEBROWSER_USETEMPLATE")));
//    templates = Gtk::manage (new MyComboBox ());
//    templateModel = Gtk::ListStore::create (templateColumns);
//    templates->set_model (templateModel);
//    templates->pack_start (templateColumns.tmplName);

//    tbox->pack_start (*useTmpl, Gtk::PACK_SHRINK, 4);
//    tbox->pack_start (*templates);
    
//    get_vbox()->pack_start (*tbox, Gtk::PACK_SHRINK, 4);
    
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

void RenameDialog::initName (const Glib::ustring& iname, const CacheImageData* cid) {

    imageData = cid;
    oldName->set_text (iname);
    newName->set_text (iname);
// Issue 316
//    if (useTmpl->get_active () && isTemplSelected ()) 
//        newName->set_text (applyTemplate (iname, cid, getActiveTemplate()));
    newName->select_region (0, newName->get_text().size());
}

Glib::ustring RenameDialog::getNewName () {

    return newName->get_text ();
}

void RenameDialog::fillTemplateList () {

    templateModel->clear ();

    for (size_t i=0; i<options.renameTemplates.size(); i++) {
        Gtk::TreeModel::iterator iter = templateModel->append ();
        iter->set_value (templateColumns.tmplName, options.renameTemplates[i]);
        iter->set_value (templateColumns.rowSeparator, false);
    }
    // append separator and the manage... item
    Gtk::TreeModel::iterator iter = templateModel->append ();
    iter->set_value (templateColumns.tmplName, Glib::ustring(""));
    iter->set_value (templateColumns.rowSeparator, true);
    iter = templateModel->append ();
    iter->set_value (templateColumns.tmplName, Glib::ustring(M("FILEBROWSER_ADDDELTEMPLATE")));
    iter->set_value (templateColumns.rowSeparator, false);
}

bool RenameDialog::rowSeparatorFunc (const Glib::RefPtr<Gtk::TreeModel>& model, const Gtk::TreeModel::iterator& iter) {

    return iter->get_value (templateColumns.rowSeparator);
}

void RenameDialog::useTemplToggled () {

    templates->set_sensitive (useTmpl->get_active ());
    if (useTmpl->get_active () && isTemplSelected ()) {
        all->set_sensitive (true);
        newName->set_text (applyTemplate (oldName->get_text(), imageData, getActiveTemplate()));
    }
    else
        all->set_sensitive (false);
    newName->select_region (0, newName->get_text().size());
}

bool RenameDialog::isTemplSelected () {

    Gtk::TreeModel::iterator iter = templates->get_active();
    return iter && iter->get_value (templateColumns.tmplName)!=M("FILEBROWSER_ADDDELTEMPLATE");
}

Glib::ustring RenameDialog::getActiveTemplate () {

    Gtk::TreeModel::iterator iter = templates->get_active();
    if (iter && iter->get_value (templateColumns.tmplName)!=M("FILEBROWSER_ADDDELTEMPLATE"))
        return iter->get_value (templateColumns.tmplName);
    else
        return "";
}

void RenameDialog::tmplSelectionChanged () {

    Gtk::TreeModel::iterator iter = templates->get_active();
    if (iter && iter->get_value (templateColumns.tmplName)==M("FILEBROWSER_ADDDELTEMPLATE")) {
        RenameTemplateEditor* rte = new RenameTemplateEditor (p);
        if (rte->run()==Gtk::RESPONSE_OK) {
            fillTemplateList ();
        }
        delete rte;
        // show add/del template dialog
    }
    else
        useTemplToggled ();
}

RenameTemplateEditor::RenameTemplateEditor (Gtk::Window* parent) 
    : Gtk::Dialog ("Edit rename templates", *parent, true, true) {
    
    list = Gtk::manage (new Gtk::ListViewText (1, false, Gtk::SELECTION_MULTIPLE));
    list->set_headers_visible (false);
    get_vbox ()->pack_start (*list);
    
    Gtk::HBox* hb = Gtk::manage (new Gtk::HBox ());
    templ = Gtk::manage (new Gtk::Entry ());
    Gtk::Button* add = Gtk::manage (new Gtk::Button ());
    Gtk::Button* del = Gtk::manage (new Gtk::Button ());
    add->add (*Gtk::manage (new RTImage ("list-add-small.png")));
    del->add (*Gtk::manage (new RTImage ("list-remove-red-small.png")));
    hb->pack_start (*templ);
    hb->pack_start (*add, Gtk::PACK_SHRINK, 2);
    hb->pack_start (*del, Gtk::PACK_SHRINK, 2);
    
    get_vbox ()->pack_start (*hb, Gtk::PACK_SHRINK, 4);

    add_button (Gtk::Stock::OK, Gtk::RESPONSE_OK);

    refreshTemplateList ();

    add->signal_pressed().connect( sigc::mem_fun(*this, &RenameTemplateEditor::addPressed) );
    del->signal_pressed().connect( sigc::mem_fun(*this, &RenameTemplateEditor::delPressed) );
    
    show_all_children ();
    
    set_size_request (-1, 250);
}

void RenameTemplateEditor::refreshTemplateList () {

    list->clear_items ();

    for (size_t i=0; i<options.renameTemplates.size(); i++)
        list->append_text (options.renameTemplates[i]);
}

void RenameTemplateEditor::addPressed () {

    if (templ->get_text()!="") {
        options.renameTemplates.push_back (templ->get_text ());
        refreshTemplateList ();
        templ->set_text("");
    }
}

void RenameTemplateEditor::delPressed () {

    std::vector<int> sel = list->get_selected ();
    for (size_t i=0; i<sel.size(); i++) {
        Glib::ustring toDel = list->get_text (sel[i]);
        std::vector<Glib::ustring>::iterator f = std::find (options.renameTemplates.begin(), options.renameTemplates.end(), toDel);
        if (f!=options.renameTemplates.end())
            options.renameTemplates.erase (f);
    }
    refreshTemplateList ();
}

Glib::ustring RenameDialog::applyTemplate (const Glib::ustring& oName, const CacheImageData* cid, const Glib::ustring& templ) {

    return Glib::ustring ("szeva");
    
}


