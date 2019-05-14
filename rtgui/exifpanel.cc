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
#include <unordered_set>

#include "exifpanel.h"

#include "guiutils.h"
#include "rtimage.h"
#include "options.h"

#include "../rtengine/imagedata.h"
#include "../rtengine/metadata.h"
#include "../rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

ExifPanel::ExifPanel() :
    idata(nullptr),
    changeList(new rtengine::procparams::ExifPairs),
    defChangeList(new rtengine::procparams::ExifPairs),
    editableTags{
        {"Exif.Photo.UserComment", "User Comment"},
        {"Exif.Image.Artist", "Artist"},
        {"Exif.Image.Copyright", "Copyright"},
        {"Exif.Image.ImageDescription", "Image Description"}
    }
{
    set_orientation(Gtk::ORIENTATION_VERTICAL);

    exifTree = Gtk::manage (new Gtk::TreeView());
    scrolledWindow = Gtk::manage (new Gtk::ScrolledWindow());

    exifTree->set_headers_visible (false);
    exifTree->set_rules_hint (false);
    exifTree->set_reorderable (false);
    exifTree->set_enable_search (false);
    exifTree->get_selection()->set_mode (Gtk::SELECTION_MULTIPLE);
    scrolledWindow->set_shadow_type (Gtk::SHADOW_NONE);
    scrolledWindow->set_policy (Gtk::POLICY_ALWAYS, Gtk::POLICY_ALWAYS);
    scrolledWindow->property_window_placement().set_value (Gtk::CORNER_TOP_LEFT);
    scrolledWindow->add (*exifTree);

    exifTreeModel = Gtk::TreeStore::create (exifColumns);
    exifTree->set_model (exifTreeModel);
    exifTree->set_grid_lines (Gtk::TREE_VIEW_GRID_LINES_NONE);
    exifTree->set_show_expanders(false);

    keepicon = RTImage::createPixbufFromFile ("tick-small.png");
    editicon = RTImage::createPixbufFromFile("add-small.png");

    Gtk::TreeView::Column *viewcol = Gtk::manage (new Gtk::TreeView::Column ("Field Name"));
    Gtk::CellRendererPixbuf* render_pb = Gtk::manage (new Gtk::CellRendererPixbuf ());
    Gtk::CellRendererText *render_txt = Gtk::manage (new Gtk::CellRendererText());
    render_txt->property_ellipsize() = Pango::ELLIPSIZE_END;
    viewcol->pack_start (*render_pb, false);
    viewcol->pack_start (*render_txt, true);
    viewcol->add_attribute (*render_pb, "pixbuf", exifColumns.icon);
    viewcol->add_attribute (*render_txt, "markup", exifColumns.label);
    viewcol->set_expand (true);
    viewcol->set_resizable (true);
    viewcol->set_fixed_width (35);
    viewcol->set_min_width (35);
    viewcol->set_sizing (Gtk::TREE_VIEW_COLUMN_AUTOSIZE);

    render_pb->property_ypad() = 0;
    render_txt->property_ypad() = 0;
    render_pb->property_yalign() = 0;
    render_txt->property_yalign() = 0;

    exifTree->append_column (*viewcol);

    Gtk::TreeView::Column *viewcolv = Gtk::manage (new Gtk::TreeView::Column ("Value"));
    Gtk::CellRendererText *render_txtv = Gtk::manage (new Gtk::CellRendererText());
    render_txtv->property_ellipsize() = Pango::ELLIPSIZE_END;
    viewcolv->pack_start (*render_txtv, true);
    viewcolv->add_attribute (*render_txtv, "markup", exifColumns.value);
    viewcolv->set_expand (true);
    viewcolv->set_resizable (true);
    viewcol->set_fixed_width (35);
    viewcolv->set_min_width (35);
    viewcolv->set_sizing (Gtk::TREE_VIEW_COLUMN_AUTOSIZE);

    render_txtv->property_ypad() = 0;

    exifTree->append_column (*viewcolv);

    pack_start (*scrolledWindow);

    Gtk::Grid* buttons1 = Gtk::manage (new Gtk::Grid());
    buttons1->set_row_homogeneous (true);
    buttons1->set_column_homogeneous (true);
    setExpandAlignProperties (buttons1, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    add = Gtk::manage (new Gtk::Button ()); // M("EXIFPANEL_ADDEDIT")
    add->set_image (*Gtk::manage (new RTImage(editicon)));
    add->set_tooltip_text (M ("EXIFPANEL_ADDEDITHINT"));
    add->get_style_context()->add_class ("Right");
    setExpandAlignProperties (add, true, true, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    buttons1->attach_next_to (*add, Gtk::POS_RIGHT, 1, 1);

    reset = Gtk::manage (new Gtk::Button ()); // M("EXIFPANEL_RESET")
    reset->set_image (*Gtk::manage (new RTImage("undo.png", "redo.png")));
    reset->set_tooltip_text (M ("EXIFPANEL_RESETHINT"));
    reset->get_style_context()->add_class ("MiddleH");
    setExpandAlignProperties (reset, true, true, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    buttons1->attach_next_to (*reset, Gtk::POS_RIGHT, 1, 1);

    resetAll = Gtk::manage (new Gtk::Button ()); // M("EXIFPANEL_RESETALL")
    resetAll->set_image (*Gtk::manage (new RTImage ("undo-all.png", "redo-all.png")));
    resetAll->set_tooltip_text (M ("EXIFPANEL_RESETALLHINT"));
    resetAll->get_style_context()->add_class ("Right");
    setExpandAlignProperties (resetAll, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    buttons1->attach_next_to (*resetAll, Gtk::POS_RIGHT, 1, 1);

    pack_end (*buttons1, Gtk::PACK_SHRINK);

    exifTree->get_selection()->signal_changed().connect (sigc::mem_fun (*this, &ExifPanel::exifSelectionChanged));
    // exifTree->signal_row_activated().connect (sigc::mem_fun (*this, &ExifPanel::row_activated));

    reset->signal_clicked().connect ( sigc::mem_fun (*this, &ExifPanel::resetPressed) );
    resetAll->signal_clicked().connect ( sigc::mem_fun (*this, &ExifPanel::resetAllPressed) );
    add->signal_clicked().connect ( sigc::mem_fun (*this, &ExifPanel::addPressed) );

    show_all ();
}

ExifPanel::~ExifPanel ()
{
}

void ExifPanel::read (const ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();

    *changeList = pp->exif;
    setImageData (idata);
    refreshTags();

    enableListener ();
}


void ExifPanel::write (ProcParams* pp, ParamsEdited* pedited)
{
    pp->exif = *changeList;
}

void ExifPanel::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    *defChangeList = defParams->exif;
}

void ExifPanel::setImageData (const FramesMetaData* id)
{

    idata = id;
}

Gtk::TreeModel::Children ExifPanel::addTag(const std::string &key, const Glib::ustring &label, const Glib::ustring &value, bool editable, bool edited)
{

    // TODO Re-fix #5923 if necessary
    //if (!value.validate()) {
    //    value = "???";
    //}

    auto root = exifTreeModel->children();

    Gtk::TreeModel::Row row = * (exifTreeModel->append (root));
    row[exifColumns.editable] = editable;
    row[exifColumns.edited] = edited;
    row[exifColumns.key] = key;
    row[exifColumns.label] = label;
    row[exifColumns.value_nopango] = value;
    row[exifColumns.value] = value;

    row[exifColumns.label] = escapeHtmlChars(label);
    row[exifColumns.value] = escapeHtmlChars(value);

    if (edited) {
        row[exifColumns.icon] = editicon;
    } else if (editable) {
        row[exifColumns.icon] = keepicon;
    }

    return row.children();
}

void ExifPanel::refreshTags()
{
    Glib::RefPtr<Gtk::TreeSelection> selection = exifTree->get_selection();
    std::vector<Gtk::TreeModel::Path> sel = selection->get_selected_rows();

    exifTreeModel->clear();

    if (!idata) {
        return;
    }

    const Glib::ustring fn = idata->getFileName();
    if (fn.empty()) {
        return;
    }

    std::unordered_set<std::string> ed;
    for (const auto &p : editableTags) {
        ed.insert(p.first);
    }

    try {
        rtengine::Exiv2Metadata meta(fn);
        meta.load();
        auto& exif = meta.exifData();

        for (const auto& p : *changeList) {
            try {
                exif[p.first] = p.second;
            } catch (const Exiv2::AnyError& exc) {
            }
        }

        for (const auto& p : editableTags) {
            const auto pos = exif.findKey(Exiv2::ExifKey(p.first));
            if (pos != exif.end() && pos->size()) {
                const bool edited = changeList->find(pos->key()) != changeList->end();
                addTag(pos->key(), pos->tagLabel(), pos->print(&exif), true, edited);
            }
        }
        for (const auto& tag : exif) {
            const bool editable = ed.find(tag.key()) != ed.end();
            if (
                !editable
                && !tag.tagLabel().empty()
                && tag.typeId() != Exiv2::undefined
                && (
                    tag.typeId() == Exiv2::asciiString
                    || tag.size() < 256
                )
            ) {
                addTag(tag.key(), tag.tagLabel(), tag.print(&exif), false, false);
            }
        }
    } catch (const Exiv2::AnyError& exc) {
        return;
    }

    for (const auto& p : sel) {
        exifTree->get_selection()->select(p);
    }
}

void ExifPanel::exifSelectionChanged ()
{

    Glib::RefPtr<Gtk::TreeSelection> selection = exifTree->get_selection();
    std::vector<Gtk::TreeModel::Path> sel = selection->get_selected_rows();

    if (sel.size() > 1) {
        reset->set_sensitive (1);
    } else if (sel.size() == 1) {
        Gtk::TreeModel::iterator iter = exifTreeModel->get_iter (sel[0]);

        if (iter->get_value(exifColumns.icon) == editicon) {
            reset->set_sensitive (1);
        }
    } else {
        reset->set_sensitive (0);
    }
}

void ExifPanel::resetIt(const Gtk::TreeModel::const_iterator& iter)
{
    if (!iter) {
        return;
    }

    const auto key = iter->get_value(exifColumns.key);
    changeList->erase(key);
}

void ExifPanel::resetPressed ()
{

    std::vector<Gtk::TreeModel::Path> sel = exifTree->get_selection()->get_selected_rows();

    for (size_t i = 0; i < sel.size(); i++) {
        resetIt(exifTreeModel->get_iter(sel[i]));
    }

    refreshTags();
    notifyListener();
}

void ExifPanel::resetAllPressed ()
{
    setImageData(idata);
    *changeList = *defChangeList;
    refreshTags();
    notifyListener ();
}

void ExifPanel::addPressed ()
{

    Gtk::Dialog* dialog = new Gtk::Dialog (M ("EXIFPANEL_ADDTAGDLG_TITLE"), * ((Gtk::Window*)get_toplevel()), true);
    dialog->add_button ("_OK", Gtk::RESPONSE_OK);
    dialog->add_button ("_Cancel", Gtk::RESPONSE_CANCEL);

    Gtk::Box* hb1 = new Gtk::Box ();
    Gtk::Box* hb2 = new Gtk::Box ();

    Gtk::Label* tlabel = new Gtk::Label (M ("EXIFPANEL_ADDTAGDLG_SELECTTAG") + ":");
    MyComboBoxText* tcombo = new MyComboBoxText ();

    for (const auto& p : editableTags) {
        tcombo->append(p.second);
    }

    hb1->pack_start (*tlabel, Gtk::PACK_SHRINK, 4);
    hb1->pack_start (*tcombo);

    Gtk::Label* vlabel = new Gtk::Label (M ("EXIFPANEL_ADDTAGDLG_ENTERVALUE") + ":");
    Gtk::Entry* ventry = new Gtk::Entry ();
    hb2->pack_start (*vlabel, Gtk::PACK_SHRINK, 4);
    hb2->pack_start (*ventry);

    Glib::ustring sel;
    Glib::ustring val;
    {
        const Glib::RefPtr<Gtk::TreeSelection> selection = exifTree->get_selection();
        const std::vector<Gtk::TreeModel::Path> rows = selection->get_selected_rows();

        if (rows.size() == 1) {
            const Gtk::TreeModel::iterator iter = exifTreeModel->get_iter(rows[0]);
            if (iter->get_value(exifColumns.editable)) {
                sel = iter->get_value(exifColumns.key);
                val = iter->get_value(exifColumns.value_nopango);
            }
        }
    }

    if (sel.empty()) {
        tcombo->set_active(0);
    } else {
        for (size_t i = 0; i < editableTags.size(); ++i) {
            if (editableTags[i].first == sel) {
                tcombo->set_active(i);
                break;
            }
        }
    }

    ventry->set_text(val);
    ventry->set_activates_default (true);
    dialog->set_default_response (Gtk::RESPONSE_OK);
    dialog->get_content_area()->pack_start (*hb1, Gtk::PACK_SHRINK);
    dialog->get_content_area()->pack_start (*hb2, Gtk::PACK_SHRINK, 4);
    tlabel->show ();
    tcombo->show ();
    vlabel->show ();
    ventry->show ();
    hb1->show ();
    hb2->show ();

    if (dialog->run () == Gtk::RESPONSE_OK) {
        auto key = editableTags[tcombo->get_active_row_number()].first;
        const auto value = ventry->get_text();
        (*changeList)[key] = value;
        refreshTags();
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

void ExifPanel::notifyListener ()
{
    if (listener) {
        listener->panelChanged (EvExif, M ("HISTORY_CHANGED"));
    }
}
