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
#include "exifpanel.h"
#include "../rtengine/safegtk.h"
#include "guiutils.h"
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;
using namespace rtexif;

ExifPanel::ExifPanel () : idata(NULL) {

   recursiveOp = true;

   exifTree = Gtk::manage(new Gtk::TreeView());
   scrolledWindow = Gtk::manage(new Gtk::ScrolledWindow());

   exifTree->set_headers_visible(false);
   exifTree->set_rules_hint(false);
   exifTree->set_reorderable(false);
   exifTree->set_enable_search(true);
   exifTree->get_selection()->set_mode (Gtk::SELECTION_MULTIPLE);
   scrolledWindow->set_border_width(2);
   scrolledWindow->set_shadow_type(Gtk::SHADOW_NONE);
   scrolledWindow->set_policy(Gtk::POLICY_ALWAYS, Gtk::POLICY_ALWAYS);
   scrolledWindow->property_window_placement().set_value(Gtk::CORNER_TOP_LEFT);
   scrolledWindow->add(*exifTree);
   
   exifTreeModel = Gtk::TreeStore::create(exifColumns);
   exifTree->set_model (exifTreeModel);

   delicon = safe_create_from_file ("gtk-close.png");
   keepicon = safe_create_from_file ("gtk-apply.png");
   editicon = safe_create_from_file ("gtk-add.png");

   Gtk::TreeView::Column *viewcol = Gtk::manage(new Gtk::TreeView::Column ("Field Name"));
   Gtk::CellRendererPixbuf* render_pb = Gtk::manage(new Gtk::CellRendererPixbuf ());
   Gtk::CellRendererText *render_txt = Gtk::manage(new Gtk::CellRendererText());
   viewcol->pack_start (*render_pb, false);
   viewcol->pack_start (*render_txt, true);
   viewcol->add_attribute (*render_pb, "pixbuf", exifColumns.icon);
   viewcol->add_attribute (*render_txt, "markup", exifColumns.field);
 
   render_pb->property_ypad() = 0;
   render_txt->property_ypad() = 0;
   render_pb->property_yalign() = 0;
   render_txt->property_yalign() = 0;
  
   exifTree->append_column (*viewcol); 
   
   Gtk::TreeView::Column *viewcolv = Gtk::manage(new Gtk::TreeView::Column ("Value"));
   Gtk::CellRendererText *render_txtv = Gtk::manage(new Gtk::CellRendererText());
   viewcolv->pack_start (*render_txtv, true);
   viewcolv->add_attribute (*render_txtv, "markup", exifColumns.value);
   
   render_txtv->property_ypad() = 0;
  
   exifTree->append_column (*viewcolv); 
  
   pack_start (*scrolledWindow);
   
   Gtk::HBox* buttons1 = Gtk::manage(new Gtk::HBox ());
   Gtk::HBox* buttons2 = Gtk::manage(new Gtk::HBox ());
   
   remove = Gtk::manage(new Gtk::Button (M("EXIFPANEL_REMOVE")));
   remove->set_image (*Gtk::manage(new Gtk::Image (delicon)));
   remove->set_tooltip_text (M("EXIFPANEL_REMOVEHINT"));
   buttons1->pack_start (*remove);
   
   keep = Gtk::manage(new Gtk::Button (M("EXIFPANEL_KEEP")));
   keep->set_image (*Gtk::manage(new Gtk::Image (keepicon)));
   keep->set_tooltip_text (M("EXIFPANEL_KEEPHINT"));
   buttons1->pack_start (*keep);

   add = Gtk::manage(new Gtk::Button (M("EXIFPANEL_ADDEDIT")));
   add->set_image (*Gtk::manage(new Gtk::Image (editicon)));
   add->set_tooltip_text (M("EXIFPANEL_ADDEDITHINT"));
   buttons1->pack_start (*add);

   reset = Gtk::manage(new Gtk::Button (M("EXIFPANEL_RESET")));
   reset->set_image (*Gtk::manage(new RTImage ("gtk-undo-ltr.png", "gtk-undo-rtl.png")));
   reset->set_tooltip_text (M("EXIFPANEL_RESETHINT"));
   buttons2->pack_start (*reset);

   resetAll = Gtk::manage(new Gtk::Button (M("EXIFPANEL_RESETALL")));
   resetAll->set_image (*Gtk::manage(new RTImage ("gtk-undoall-ltr.png", "gtk-undoall-rtl.png")));
   resetAll->set_tooltip_text (M("EXIFPANEL_RESETALLHINT"));
   buttons2->pack_start (*resetAll);

   pack_end (*buttons2, Gtk::PACK_SHRINK);
   pack_end (*buttons1, Gtk::PACK_SHRINK);

   exifTree->get_selection()->signal_changed().connect(sigc::mem_fun(*this, &ExifPanel::exifSelectionChanged));
   exifTree->signal_row_activated().connect(sigc::mem_fun(*this, &ExifPanel::row_activated));

   remove->signal_clicked().connect( sigc::mem_fun(*this, &ExifPanel::removePressed) );
   keep->signal_clicked().connect( sigc::mem_fun(*this, &ExifPanel::keepPressed) );
   reset->signal_clicked().connect( sigc::mem_fun(*this, &ExifPanel::resetPressed) );
   resetAll->signal_clicked().connect( sigc::mem_fun(*this, &ExifPanel::resetAllPressed) );
   add->signal_clicked().connect( sigc::mem_fun(*this, &ExifPanel::addPressed) );
   
   show_all ();
}

ExifPanel::~ExifPanel () {
}

void ExifPanel::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    changeList = pp->exif;
    setImageData (idata);
    applyChangeList ();
    exifSelectionChanged ();
    
    enableListener ();
}

void ExifPanel::write (ProcParams* pp, ParamsEdited* pedited) {
    
//    updateChangeList ();
    pp->exif = changeList;
}

void ExifPanel::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    defChangeList = defParams->exif;
}
        
void ExifPanel::setImageData (const ImageMetaData* id) {
    
    idata = id; 
    exifTreeModel->clear ();

    const std::vector<Tag*>& defTags = ExifManager::getDefaultTIFFTags (NULL);
    for (size_t i=0; i<defTags.size(); i++)
        if (defTags[i]->nameToString() == "ImageWidth" || defTags[i]->nameToString() == "ImageHeight" || defTags[i]->nameToString() == "BitsPerSample")
            addTag (exifTreeModel->children(), defTags[i]->nameToString(), "?", AC_SYSTEM, false);
        else
            addTag (exifTreeModel->children(), defTags[i]->nameToString(), defTags[i]->valueToString(), AC_SYSTEM, false);
    
    if (id && id->getExifData ()) {
//        id->getExifData ()->printAll ();     
        addDirectory (id->getExifData (), exifTreeModel->children());
    }
}

Gtk::TreeModel::Children ExifPanel::addTag (const Gtk::TreeModel::Children& root, Glib::ustring field, Glib::ustring value, rtexif::ActionCode action, bool editable) {

    Gtk::TreeModel::Row row = *(exifTreeModel->append(root));
    row[exifColumns.action]   = action;
    row[exifColumns.editable] = editable;
    row[exifColumns.edited]   = false;
    row[exifColumns.field_nopango] = field;
    row[exifColumns.value_nopango] = value;
    row[exifColumns.orig_value]    = value;
    
    if (action==AC_WRITE)
        row[exifColumns.icon] = keepicon;
    else if (action==AC_DONTWRITE)
        row[exifColumns.icon] = delicon;
    
    if (editable) {
        row[exifColumns.field] = Glib::ustring("<b>") + escapeHtmlChars(field) + "</b>";
        row[exifColumns.value] = Glib::ustring("<b>") + escapeHtmlChars(value) + "</b>";
    }
    else if (action==AC_SYSTEM) {
        row[exifColumns.field] = Glib::ustring("<i>") + escapeHtmlChars(field) + "</i>";
        row[exifColumns.value] = Glib::ustring("<i>") + escapeHtmlChars(value) + "</i>";
    }
    else {
        row[exifColumns.field] = escapeHtmlChars(field);
        row[exifColumns.value] = escapeHtmlChars(value);
    }
    
    return row.children();
}

void ExifPanel::addDirectory (const TagDirectory* dir, Gtk::TreeModel::Children root) {

    for (int i=0; i<dir->getCount(); i++) {
	Tag* t = (const_cast<TagDirectory*>(dir))->getTagByIndex (i);
        if (t->getAttrib() && t->getAttrib()->action==AC_SYSTEM)
            continue;
        if (t->isDirectory()) 
            for (int j=0; t->getDirectory(j); j++) {
                Gtk::TreeModel::Children ch = addTag (root, t->nameToString (j), M("EXIFPANEL_SUBDIRECTORY"), t->getAttrib() ? t->getAttrib()->action : AC_DONTWRITE, t->getAttrib() && t->getAttrib()->editable);
                addDirectory (t->getDirectory(j), ch);
            }
        else 
            addTag (root, t->nameToString (), t->valueToString (), t->getAttrib() ? (t->getOwnMemory()?t->getAttrib()->action:AC_SYSTEM) : AC_DONTWRITE, t->getAttrib() && t->getAttrib()->editable);
    }
}

void ExifPanel::exifSelectionChanged () {

    Glib::RefPtr<Gtk::TreeSelection> selection = exifTree->get_selection();
    std::vector<Gtk::TreeModel::Path> sel = selection->get_selected_rows();
    if (sel.size()>1) {
        remove->set_sensitive (1);
        keep->set_sensitive (1);
        reset->set_sensitive (1);
    }
    else if (sel.size()==1) {
        Gtk::TreeModel::iterator iter = exifTreeModel->get_iter (sel[0]);
        if (iter->get_value (exifColumns.action)==AC_SYSTEM) {
            remove->set_sensitive (0);
            keep->set_sensitive (0);
            reset->set_sensitive (0);
        }
        else if (!iter->children().empty()) {
            remove->set_sensitive (1);
            keep->set_sensitive (1);
            reset->set_sensitive (1);
        }
        else if (iter->get_value(exifColumns.icon)==delicon) {
            remove->set_sensitive (0);
            keep->set_sensitive (1);
            reset->set_sensitive (1);
        }
        else if (iter->get_value(exifColumns.icon)==keepicon || iter->get_value(exifColumns.icon)==editicon) {
            keep->set_sensitive (0);
            remove->set_sensitive (1);
            reset->set_sensitive (1);
        }
    }
    else {
        remove->set_sensitive (0);
        keep->set_sensitive (0);
        reset->set_sensitive (0);
    }
}

void ExifPanel::delIt (Gtk::TreeModel::iterator iter) {
    
    if (!iter)
        return;

    if (iter->get_value (exifColumns.action) != AC_SYSTEM)
        iter->set_value (exifColumns.icon, delicon);
    if (recursiveOp)
        for (Gtk::TreeModel::iterator i=iter->children().begin(); i!=iter->children().end(); i++)
            delIt (i);
}

void ExifPanel::removePressed () {

    std::vector<Gtk::TreeModel::Path> sel = exifTree->get_selection()->get_selected_rows();
    for (size_t i=0; i<sel.size(); i++)
        delIt (exifTreeModel->get_iter (sel[i]));

    exifSelectionChanged ();
    updateChangeList ();
    notifyListener ();
}

void ExifPanel::keepIt (Gtk::TreeModel::iterator iter) {
    
    if (!iter)
        return;

    if (iter->get_value (exifColumns.action) != AC_SYSTEM)
        iter->set_value (exifColumns.icon, iter->get_value (exifColumns.edited) ? editicon : keepicon);
    if (recursiveOp)
        for (Gtk::TreeModel::iterator i=iter->children().begin(); i!=iter->children().end(); i++)
            keepIt (i);
}

void ExifPanel::keepPressed () {

    std::vector<Gtk::TreeModel::Path> sel = exifTree->get_selection()->get_selected_rows();
    for (size_t i=0; i<sel.size(); i++)
        keepIt (exifTreeModel->get_iter (sel[i]));

    exifSelectionChanged ();
    updateChangeList ();
    notifyListener ();
}

/*void ExifPanel::resetIt (Gtk::TreeModel::iterator  iter) {
    
    if (!iter)
        return;

    if (iter->get_value (exifColumns.action)!=AC_SYSTEM)
        iter->set_value (exifColumns.icon, iter->get_value (exifColumns.action) ? keepicon : delicon);
    if (iter->get_value (exifColumns.edited)) {
        iter->set_value (exifColumns.value, Glib::ustring("<b>") + iter->get_value(exifColumns.orig_value) + "</b>");
        iter->set_value (exifColumns.value_nopango, iter->get_value(exifColumns.orig_value));
        iter->set_value (exifColumns.edited, false);
    }
    if (iter->get_value (exifColumns.action)==AC_INVALID)
        exifTreeModel->erase (iter);
    else 
    if (recursiveOp)
        for (Gtk::TreeModel::iterator i=iter->children().begin(); i!=iter->children().end(); i++)
            resetIt (i);
}*/
Gtk::TreeModel::iterator ExifPanel::resetIt (Gtk::TreeModel::iterator  iter) {
    
    if (!iter)
        return iter;

    if (iter->get_value (exifColumns.action)!=AC_SYSTEM)
        iter->set_value (exifColumns.icon, iter->get_value (exifColumns.action) ? keepicon : delicon);
    if (iter->get_value (exifColumns.edited)) {
        iter->set_value (exifColumns.value, Glib::ustring("<b>") + iter->get_value(exifColumns.orig_value) + "</b>");
        iter->set_value (exifColumns.value_nopango, iter->get_value(exifColumns.orig_value));
        iter->set_value (exifColumns.edited, false);
    }
    if (iter->get_value (exifColumns.action)==AC_INVALID) {
        return exifTreeModel->erase (iter);
    }
    else 
    if (recursiveOp) {
        Gtk::TreeModel::iterator i = iter->children().begin();
        while (i && i != iter->children().end())
            i = resetIt (i);
    }
    return ++iter;
}
void ExifPanel::resetPressed () {

    std::vector<Gtk::TreeModel::Path> sel = exifTree->get_selection()->get_selected_rows();
    for (size_t i=0; i<sel.size(); i++)
        resetIt (exifTreeModel->get_iter (sel[i]));

    exifSelectionChanged ();
    updateChangeList ();
    notifyListener ();
}

void ExifPanel::resetAllPressed () {

    setImageData (idata);
    changeList = defChangeList;
    applyChangeList ();
    exifSelectionChanged ();
    notifyListener ();
}

void ExifPanel::addPressed () {

    Gtk::Dialog* dialog = new Gtk::Dialog (M("EXIFPANEL_ADDTAGDLG_TITLE"), *((Gtk::Window*)get_toplevel()), true, true);
    dialog->add_button (Gtk::Stock::OK, Gtk::RESPONSE_OK);   
    dialog->add_button (Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
    
    Gtk::HBox* hb1 = new Gtk::HBox ();
    Gtk::HBox* hb2 = new Gtk::HBox ();
    
    Gtk::Label* tlabel = new Gtk::Label (M("EXIFPANEL_ADDTAGDLG_SELECTTAG")+":");
    MyComboBoxText* tcombo = new MyComboBoxText ();
    
    tcombo->append_text ("Artist");
    tcombo->append_text ("Copyright");
    tcombo->append_text ("ImageDescription");
    tcombo->append_text ("Exif.UserComment");
    
    hb1->pack_start (*tlabel, Gtk::PACK_SHRINK, 4);
    hb1->pack_start (*tcombo);

    Gtk::Label* vlabel = new Gtk::Label (M("EXIFPANEL_ADDTAGDLG_ENTERVALUE")+":");
    Gtk::Entry* ventry = new Gtk::Entry ();
    hb2->pack_start (*vlabel, Gtk::PACK_SHRINK, 4);
    hb2->pack_start (*ventry);

    Glib::ustring sel = getSelection (true);
    if (sel=="")
        tcombo->set_active_text ("Exif.UserComment");
    else {
        tcombo->set_active_text (sel);
        if (tcombo->get_active ()<0) {
            tcombo->append_text (sel);
            tcombo->set_active_text (sel);
        }
        ventry->set_text (getSelectedValue ());
    }
    
    ventry->set_activates_default (true);
    dialog->set_default_response (Gtk::RESPONSE_OK);
    dialog->get_vbox()->pack_start (*hb1, Gtk::PACK_SHRINK);
    dialog->get_vbox()->pack_start (*hb2, Gtk::PACK_SHRINK, 4);
    tlabel->show ();
    tcombo->show ();
    vlabel->show ();
    ventry->show ();
    hb1->show ();
    hb2->show ();
    
    if (dialog->run ()== Gtk::RESPONSE_OK) {   
        editTag (exifTreeModel->children(), tcombo->get_active_text(), ventry->get_text());
        updateChangeList ();
        notifyListener ();
    }
    
    delete dialog;
    delete tlabel;
    delete tcombo;  
    delete vlabel;
    delete ventry;  
    delete hb1;
    delete hb2;  
}

void ExifPanel::editTag (Gtk::TreeModel::Children root, Glib::ustring name, Glib::ustring value) {

    Glib::ustring::size_type dp = name.find_first_of ('.');
    Glib::ustring fseg = name.substr (0,dp);
    // look up first segment of the path
    Gtk::TreeModel::iterator iter;
    for (iter = root.begin(); iter!=root.end(); iter++) 
        if (iter->get_value (exifColumns.field_nopango) == fseg)
            break;
    
    if (iter==root.end() && value!="#keep" && value!="#delete") {
        iter = exifTreeModel->append(root);
        iter->set_value (exifColumns.field_nopango, fseg);
        iter->set_value (exifColumns.action, AC_INVALID);
        if (dp==Glib::ustring::npos) {
            iter->set_value (exifColumns.value, Glib::ustring("<b>") + value + "</b>");
            iter->set_value (exifColumns.value_nopango, value);
            iter->set_value (exifColumns.orig_value, value);
            iter->set_value (exifColumns.field, Glib::ustring("<b>") + fseg + "</b>");
            iter->set_value (exifColumns.edited, true);
            iter->set_value (exifColumns.editable, true);
            iter->set_value (exifColumns.icon, editicon);
        }
        else {
            iter->set_value (exifColumns.value, Glib::ustring(M("EXIFPANEL_SUBDIRECTORY")));
            iter->set_value (exifColumns.value_nopango, Glib::ustring(M("EXIFPANEL_SUBDIRECTORY")));
            iter->set_value (exifColumns.field, fseg);
            iter->set_value (exifColumns.icon, keepicon);
            iter->set_value (exifColumns.orig_value, Glib::ustring(M("EXIFPANEL_SUBDIRECTORY")));
        }
    }
    
    if (dp==Glib::ustring::npos) {
        if (value=="#keep" && iter->get_value (exifColumns.action)!=AC_SYSTEM)
            iter->set_value (exifColumns.icon, iter->get_value (exifColumns.edited) ? editicon : keepicon);
        else if (value=="#delete" && iter->get_value (exifColumns.action)!=AC_SYSTEM)
            iter->set_value (exifColumns.icon, delicon);
        else {
            iter->set_value (exifColumns.value, Glib::ustring("<b>") + value + "</b>");
            iter->set_value (exifColumns.value_nopango, value);
            iter->set_value (exifColumns.edited, true);
            iter->set_value (exifColumns.icon, editicon);
        }
    }
    else 
        editTag (iter->children(), name.substr (dp+1, Glib::ustring::npos), value);
}

Glib::ustring ExifPanel::getSelectedValue () {

    Glib::RefPtr<Gtk::TreeSelection> selection = exifTree->get_selection();
    std::vector<Gtk::TreeModel::Path> rows = selection->get_selected_rows(); 
    if (rows.size()!=1)
        return "";   
    Gtk::TreeModel::iterator iter = exifTreeModel->get_iter (rows[0]);
    if (iter) 
        return iter->get_value (exifColumns.value_nopango);
    return "";
}

Glib::ustring ExifPanel::getSelection (bool onlyeditable) {

    Glib::RefPtr<Gtk::TreeSelection> selection = exifTree->get_selection();
    std::vector<Gtk::TreeModel::Path> rows = selection->get_selected_rows(); 

    if (rows.size()!=1)
        return "";   
    Gtk::TreeModel::iterator iter = exifTreeModel->get_iter (rows[0]);

    Glib::ustring ret = "";
    bool first = true;
    bool editable = false;
    while (iter) {
        if (first)
            ret = iter->get_value (exifColumns.field_nopango);
        else
            ret = iter->get_value (exifColumns.field_nopango) + "." + ret;
        editable = iter->get_value (exifColumns.editable);
        iter = iter->parent ();
        first = false;
    }
    if (!editable && onlyeditable)
        return "";
    return ret;
}

void ExifPanel::updateChangeList (Gtk::TreeModel::Children root, std::string prefix) {

    if (prefix!="")
        prefix = prefix + ".";

    Gtk::TreeModel::iterator iter;
    for (iter = root.begin(); iter!=root.end(); iter++)  {
        if (iter->get_value (exifColumns.edited) == true)
            changeList[ prefix+iter->get_value (exifColumns.field_nopango) ] = iter->get_value (exifColumns.value_nopango);
        else if (iter->get_value (exifColumns.action) == AC_WRITE && iter->get_value (exifColumns.icon) == delicon)
            changeList[ prefix+iter->get_value (exifColumns.field_nopango) ] = "#delete";
        else if (iter->get_value (exifColumns.action) == AC_DONTWRITE && iter->get_value (exifColumns.icon) == keepicon)
            changeList[ prefix+iter->get_value (exifColumns.field_nopango) ] = "#keep";
        if (iter->get_value (exifColumns.icon) == keepicon)
            updateChangeList (iter->children(), prefix + iter->get_value (exifColumns.field_nopango));
    }
}

void ExifPanel::updateChangeList () {

    changeList.clear ();
    updateChangeList (exifTreeModel->children(), "");
}

void ExifPanel::applyChangeList () {

    for (rtengine::procparams::ExifPairs::iterator i=changeList.begin(); i!=changeList.end(); i++)
        editTag (exifTreeModel->children(), i->first, i->second);
}

void ExifPanel::row_activated (const Gtk::TreeModel::Path& path, Gtk::TreeViewColumn* column) {

    Gtk::TreeModel::iterator iter = exifTreeModel->get_iter (path);
    if (iter) {
        if (!iter->children().empty())
            if (exifTree->row_expanded (path))
                exifTree->collapse_row (path);
            else
                exifTree->expand_row (path, false);
        else if (iter->get_value (exifColumns.editable))
            addPressed ();
    }
}


void ExifPanel::notifyListener () {
    
    if (listener)
        listener->panelChanged (EvExif, M("HISTORY_CHANGED"));
}
