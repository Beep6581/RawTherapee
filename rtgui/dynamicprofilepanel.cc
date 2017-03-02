/* -*- C++ -*-
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2017 Alberto Griggio
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

#include "dynamicprofilepanel.h"
#include "multilangmgr.h"
#include <sstream>
#include <iomanip>

namespace {

template <class V>
Glib::ustring to_str(V n)
{
    std::ostringstream buf;
    buf << std::setprecision(1) << std::fixed << n;
    return buf.str();
}


int to_int(const Glib::ustring &s)
{
    std::istringstream buf(s);
    int r = -1;
    buf >> r;
    return r;
}


double to_double(const Glib::ustring &s)
{
    std::istringstream buf(s);
    double r = 0.0;
    buf >> r;
    return r;
}

} // namespace


//-----------------------------------------------------------------------------
// DynamicProfilePanel::EditDialog
//-----------------------------------------------------------------------------

DynamicProfilePanel::EditDialog::EditDialog(const Glib::ustring &title,
                                            Gtk::Window &parent):
    Gtk::Dialog(title, parent)
{
    profilepath_ = Gtk::manage(new ProfileStoreComboBox());
    Gtk::HBox *hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("DYNPROFILEEDITOR_PROFILE"))),
                   false, false, 4);
    hb->pack_start(*profilepath_, true, true, 2);
    get_content_area()->pack_start(*hb, Gtk::PACK_SHRINK, 4);

    add_optional(M("DYNPROFILEEDITOR_CAMERA_MAKE"), has_make, make);
    add_optional(M("DYNPROFILEEDITOR_CAMERA_MODEL"), has_model, model);
    add_optional(M("EXIFFILTER_LENS"), has_lens, lens);
            
    add_range(M("EXIFFILTER_ISO"), iso_min_, iso_max_);
    add_range(M("EXIFFILTER_APERTURE"), fnumber_min_, fnumber_max_);
    add_range(M("EXIFFILTER_FOCALLEN"), focallen_min_, focallen_max_);
    add_range(M("EXIFFILTER_SHUTTER"), shutterspeed_min_, shutterspeed_max_);
    add_range(M("EXIFFILTER_EXPOSURECOMPENSATION"), expcomp_min_, expcomp_max_);

    add_button(M("GENERAL_OK"), 1);
    add_button(M("GENERAL_CANCEL"), 2);

    set_ranges();
            
    show_all_children();
}


void DynamicProfilePanel::EditDialog::set_entry(
    const DynamicProfileEntry &entry)
{
    iso_min_->set_value(entry.iso.min);
    iso_max_->set_value(entry.iso.max);

    fnumber_min_->set_value(entry.fnumber.min);
    fnumber_max_->set_value(entry.fnumber.max);

    focallen_min_->set_value(entry.focallen.min);
    focallen_max_->set_value(entry.focallen.max);

    shutterspeed_min_->set_value(entry.shutterspeed.min);
    shutterspeed_max_->set_value(entry.shutterspeed.max);

    expcomp_min_->set_value(entry.expcomp.min);
    expcomp_max_->set_value(entry.expcomp.max);

    has_make->set_active(entry.make.enabled);
    make->set_text(entry.make.value);

    has_model->set_active(entry.model.enabled);
    model->set_text(entry.model.value);

    has_lens->set_active(entry.lens.enabled);
    lens->set_text(entry.lens.value);

    profilepath_->updateProfileList();
    if (!profilepath_->setActiveRowFromFullPath(entry.profilepath)) {
        profilepath_->setInternalEntry();
    }
}
        

DynamicProfileEntry DynamicProfilePanel::EditDialog::get_entry()
{
    DynamicProfileEntry ret;
    ret.iso.min = iso_min_->get_value_as_int();
    ret.iso.max = iso_max_->get_value_as_int();

    ret.fnumber.min = fnumber_min_->get_value();
    ret.fnumber.max = fnumber_max_->get_value();
            
    ret.focallen.min = focallen_min_->get_value();
    ret.focallen.max = focallen_max_->get_value();

    ret.shutterspeed.min = shutterspeed_min_->get_value();
    ret.shutterspeed.max = shutterspeed_max_->get_value();

    ret.expcomp.min = expcomp_min_->get_value();
    ret.expcomp.max = expcomp_max_->get_value();

    ret.make.enabled = has_make->get_active();
    ret.make.value = make->get_text();

    ret.model.enabled = has_model->get_active();
    ret.model.value = model->get_text();

    ret.lens.enabled = has_lens->get_active();
    ret.lens.value = lens->get_text();
 
    ret.profilepath = profilepath_->getFullPathFromActiveRow();

    return ret;
}

void DynamicProfilePanel::EditDialog::set_ranges()
{
    DynamicProfileEntry default_entry;
    iso_min_->set_digits(0);
    iso_max_->set_digits(0);
    iso_min_->set_increments(1, 10);
    iso_max_->set_increments(1, 10);
    iso_min_->set_range(default_entry.iso.min, default_entry.iso.max);
    iso_max_->set_range(default_entry.iso.min, default_entry.iso.max);
    iso_min_->set_value(default_entry.iso.min);
    iso_max_->set_value(default_entry.iso.max);

#define DOIT_(name)                                     \
    name ## _min_->set_digits(1);                       \
    name ## _max_->set_digits(1);                       \
    name ## _min_->set_increments(0.1, 1);              \
    name ## _max_->set_increments(0.1, 1);              \
    name ## _min_->set_range(default_entry. name .min,  \
                             default_entry. name .max); \
    name ## _max_->set_range(default_entry. name .min,  \
                             default_entry. name .max); \
    name ## _min_->set_value(default_entry. name .min); \
    name ## _max_->set_value(default_entry. name .max)

    DOIT_(fnumber);
    DOIT_(focallen);
    DOIT_(shutterspeed);
    DOIT_(expcomp);
#undef DOIT_

    profilepath_->setInternalEntry();
}


void DynamicProfilePanel::EditDialog::add_range(const Glib::ustring &name,
        Gtk::SpinButton *&from, Gtk::SpinButton *&to)
{
    Gtk::HBox *hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*Gtk::manage(new Gtk::Label(name)), false, false, 4);
    from = Gtk::manage(new Gtk::SpinButton());
    to = Gtk::manage(new Gtk::SpinButton());
    from->set_numeric(true);
    to->set_numeric(true);
    hb->pack_start(*from, true, true, 2);
    hb->pack_start(*Gtk::manage(new Gtk::Label(" - ")),
        false, false, 4);
    hb->pack_start(*to, true, true, 2);
    get_content_area()->pack_start(*hb, Gtk::PACK_SHRINK, 4);
}


void DynamicProfilePanel::EditDialog::add_optional(const Glib::ustring &name,
        Gtk::CheckButton *&check, Gtk::Entry *&field)
{
    check = Gtk::manage (new Gtk::CheckButton(name));
    Gtk::HBox *hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*check, Gtk::PACK_SHRINK, 4);
    field = Gtk::manage(new Gtk::Entry());
    hb->pack_start(*field, true, true, 2);
    get_content_area()->pack_start(*hb, Gtk::PACK_SHRINK, 4);
}


//-----------------------------------------------------------------------------
// DynamicProfilePanel
//-----------------------------------------------------------------------------

DynamicProfilePanel::DynamicProfilePanel():
    vbox_(Gtk::ORIENTATION_VERTICAL),
    button_up_(M("DYNPROFILEEDITOR_MOVE_UP")),
    button_down_(M("DYNPROFILEEDITOR_MOVE_DOWN")),
    button_new_(M("DYNPROFILEEDITOR_NEW")),
    button_edit_(M("DYNPROFILEEDITOR_EDIT")),
    button_delete_(M("DYNPROFILEEDITOR_DELETE"))
{
    add(vbox_);

    //Add the TreeView, inside a ScrolledWindow, with the button underneath:
    scrolledwindow_.add(treeview_);

    //Only show the scrollbars when they are necessary:
    scrolledwindow_.set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);

    vbox_.pack_start(scrolledwindow_);
    vbox_.pack_start(buttonbox_, Gtk::PACK_SHRINK);

    buttonbox_.pack_start(button_new_, Gtk::PACK_SHRINK);
    buttonbox_.pack_start(button_edit_, Gtk::PACK_SHRINK);
    buttonbox_.pack_start(button_delete_, Gtk::PACK_SHRINK);
    buttonbox_.pack_start(button_up_, Gtk::PACK_SHRINK);
    buttonbox_.pack_start(button_down_, Gtk::PACK_SHRINK);
    buttonbox_.set_border_width(5);
    buttonbox_.set_layout(Gtk::BUTTONBOX_END);
    button_up_.signal_clicked().connect(
        sigc::mem_fun(*this, &DynamicProfilePanel::on_button_up));
    button_down_.signal_clicked().connect(
        sigc::mem_fun(*this, &DynamicProfilePanel::on_button_down));
    button_new_.signal_clicked().connect(
        sigc::mem_fun(*this, &DynamicProfilePanel::on_button_new));
    button_edit_.signal_clicked().connect(
        sigc::mem_fun(*this, &DynamicProfilePanel::on_button_edit));
    button_delete_.signal_clicked().connect(
        sigc::mem_fun(*this, &DynamicProfilePanel::on_button_delete));

    //Create the Tree model:
    treemodel_ = Gtk::ListStore::create(columns_);
    treeview_.set_model(treemodel_);

    //Add the TreeView's view columns:
    auto cell = Gtk::manage(new Gtk::CellRendererText());
    int cols_count = treeview_.append_column(
        M("DYNPROFILEEDITOR_PROFILE"), *cell);
    auto col = treeview_.get_column(cols_count - 1);
    if (col) {
        col->set_cell_data_func(
            *cell, sigc::mem_fun(
                *this, &DynamicProfilePanel::render_profilepath));
    }
    cell = Gtk::manage(new Gtk::CellRendererText());
    cols_count = treeview_.append_column(
        M("DYNPROFILEEDITOR_CAMERA_MAKE"), *cell);
    col = treeview_.get_column(cols_count - 1);
    if (col) {
        col->set_cell_data_func(
            *cell, sigc::mem_fun(
                *this, &DynamicProfilePanel::render_make));
    }
    cell = Gtk::manage(new Gtk::CellRendererText());
    cols_count = treeview_.append_column(
        M("DYNPROFILEEDITOR_CAMERA_MODEL"), *cell);
    col = treeview_.get_column(cols_count - 1);
    if (col) {
        col->set_cell_data_func(
            *cell, sigc::mem_fun(
                *this, &DynamicProfilePanel::render_model));
    }
    cell = Gtk::manage(new Gtk::CellRendererText());
    cols_count = treeview_.append_column(M("EXIFFILTER_LENS"), *cell);
    col = treeview_.get_column(cols_count - 1);
    if (col) {
        col->set_cell_data_func(
            *cell, sigc::mem_fun(
                *this, &DynamicProfilePanel::render_lens));
    }
    cell = Gtk::manage(new Gtk::CellRendererText());
    cols_count = treeview_.append_column(M("EXIFFILTER_ISO"), *cell);
    col = treeview_.get_column(cols_count - 1);
    if (col) {
        col->set_cell_data_func(
            *cell, sigc::mem_fun(
                *this, &DynamicProfilePanel::render_iso));
    }
    cell = Gtk::manage(new Gtk::CellRendererText());
    cols_count = treeview_.append_column(M("EXIFFILTER_APERTURE"), *cell);
    col = treeview_.get_column(cols_count - 1);
    if (col) {
        col->set_cell_data_func(
            *cell, sigc::mem_fun(
                *this, &DynamicProfilePanel::render_fnumber));
    }
    cell = Gtk::manage(new Gtk::CellRendererText());
    cols_count = treeview_.append_column(M("EXIFFILTER_FOCALLEN"), *cell);
    col = treeview_.get_column(cols_count - 1);
    if (col) {
        col->set_cell_data_func(
            *cell, sigc::mem_fun(
                *this, &DynamicProfilePanel::render_focallen));
    }
    cell = Gtk::manage(new Gtk::CellRendererText());
    cols_count = treeview_.append_column(M("EXIFFILTER_SHUTTER"), *cell);
    col = treeview_.get_column(cols_count - 1);
    if (col) {
        col->set_cell_data_func(
            *cell, sigc::mem_fun(
                *this, &DynamicProfilePanel::render_shutterspeed));
    }
    cell = Gtk::manage(new Gtk::CellRendererText());
    cols_count = treeview_.append_column(
        M("EXIFFILTER_EXPOSURECOMPENSATION"), *cell);
    col = treeview_.get_column(cols_count - 1);
    if (col) {
        col->set_cell_data_func(
            *cell, sigc::mem_fun(
                *this, &DynamicProfilePanel::render_expcomp));
    }
    
    show_all_children();

    std::vector<DynamicProfileEntry> entries;
    if (loadDynamicProfileEntries(entries)) {
        for (auto &e : entries) {
            add_entry(e);
        }
    }
}


void DynamicProfilePanel::update_entry(Gtk::TreeModel::Row row,
                                       const DynamicProfileEntry &entry)
{
    row[columns_.iso] = entry.iso;
    row[columns_.fnumber] = entry.fnumber;
    row[columns_.focallen] = entry.focallen;
    row[columns_.shutterspeed] = entry.shutterspeed;
    row[columns_.expcomp] = entry.expcomp;
    row[columns_.make] = entry.make;
    row[columns_.model] = entry.model;
    row[columns_.lens] = entry.lens;
    row[columns_.profilepath] = entry.profilepath;
}

void DynamicProfilePanel::add_entry(const DynamicProfileEntry &entry)
{
    auto row = *(treemodel_->append());
    update_entry(row, entry);
}


DynamicProfileEntry DynamicProfilePanel::to_entry(Gtk::TreeModel::Row row,
                                                  int serial)
{
    DynamicProfileEntry ret;
    ret.serial_number = serial;
    ret.iso = row[columns_.iso];
    ret.fnumber = row[columns_.fnumber];
    ret.focallen = row[columns_.focallen];
    ret.shutterspeed = row[columns_.shutterspeed];
    ret.expcomp = row[columns_.expcomp];
    ret.make = row[columns_.make];
    ret.model = row[columns_.model];
    ret.lens = row[columns_.lens];
    ret.profilepath = row[columns_.profilepath];
    return ret;
}


void DynamicProfilePanel::render_profilepath(
    Gtk::CellRenderer *cell, const Gtk::TreeModel::iterator &iter)
{
    auto row = *iter;
    Gtk::CellRendererText *ct = static_cast<Gtk::CellRendererText *>(cell);
    auto value = row[columns_.profilepath];
    ct->property_text() = value;
}


#define RENDER_RANGE_(tp, name)                                          \
    auto row = *iter;                                                   \
    Gtk::CellRendererText *ct = static_cast<Gtk::CellRendererText *>(cell); \
    DynamicProfileEntry::Range<tp> r = row[columns_. name];         \
    auto value = to_str(r.min) + " - " + to_str(r.max); \
    ct->property_text() = value;

void DynamicProfilePanel::render_iso(
    Gtk::CellRenderer *cell, const Gtk::TreeModel::iterator &iter)
{
    RENDER_RANGE_(int, iso);
}


void DynamicProfilePanel::render_fnumber(
    Gtk::CellRenderer *cell, const Gtk::TreeModel::iterator &iter)
{
    RENDER_RANGE_(double, fnumber);
}


void DynamicProfilePanel::render_focallen(
    Gtk::CellRenderer *cell, const Gtk::TreeModel::iterator &iter)
{
    RENDER_RANGE_(double, focallen);
}


void DynamicProfilePanel::render_shutterspeed(
    Gtk::CellRenderer *cell, const Gtk::TreeModel::iterator &iter)
{
    RENDER_RANGE_(double, shutterspeed);
}


void DynamicProfilePanel::render_expcomp(
    Gtk::CellRenderer *cell, const Gtk::TreeModel::iterator &iter)
{
    RENDER_RANGE_(double, expcomp);
}

#undef RENDER_RANGE_

#define RENDER_OPTIONAL_(name) \
    auto row = *iter; \
    Gtk::CellRendererText *ct = static_cast<Gtk::CellRendererText *>(cell); \
    DynamicProfileEntry::Optional<Glib::ustring> o = row[columns_. name]; \
    if (o.enabled) { \
        ct->property_text() = o.value; \
    } else { \
        ct->property_text() = ""; \
    }

void DynamicProfilePanel::render_make(
    Gtk::CellRenderer *cell, const Gtk::TreeModel::iterator &iter)
{
    RENDER_OPTIONAL_(make);
}


void DynamicProfilePanel::render_model(
    Gtk::CellRenderer *cell, const Gtk::TreeModel::iterator &iter)
{
    RENDER_OPTIONAL_(model);
}


void DynamicProfilePanel::render_lens(
    Gtk::CellRenderer *cell, const Gtk::TreeModel::iterator &iter)
{
    RENDER_OPTIONAL_(lens);
}

#undef RENDER_OPTIONAL_

void DynamicProfilePanel::on_button_up()
{
    auto s = treeview_.get_selection();
    if (!s->count_selected_rows()) {
        return;
    }
    auto it = s->get_selected();
    if (it != treemodel_->children().begin()) {
        auto it2 = it;
        --it2;
        treemodel_->iter_swap(it, it2);
    }
}

void DynamicProfilePanel::on_button_down()
{
    auto s = treeview_.get_selection();
    if (!s->count_selected_rows()) {
        return;
    }
    auto it = s->get_selected();
    auto it2 = it;
    ++it2;
    if (it2 != treemodel_->children().end()) {
        treemodel_->iter_swap(it, it2);        
    }
}


void DynamicProfilePanel::on_button_delete()
{
    auto s = treeview_.get_selection();
    if (!s->count_selected_rows()) {
        return;
    }
    auto it = s->get_selected();
    treemodel_->erase(it);
}


void DynamicProfilePanel::on_button_new()
{
    EditDialog d(M("DYNPROFILEEDITOR_NEW_RULE"),
                 static_cast<Gtk::Window &>(*get_toplevel())); 
    int status = d.run();
    if (status == 1) {
        DynamicProfileEntry entry = d.get_entry();
        add_entry(entry);
    }
}


void DynamicProfilePanel::on_button_edit()
{
    auto s = treeview_.get_selection();
    if (!s->count_selected_rows()) {
        return;
    }
    EditDialog d(M("DYNPROFILEEDITOR_EDIT_RULE"),
                 static_cast<Gtk::Window &>(*get_toplevel()));
    auto it = s->get_selected();
    Gtk::TreeModel::Row row = *(s->get_selected());
    d.set_entry(to_entry(row));
    int status = d.run();
    if (status == 1) {
        update_entry(row, d.get_entry());
    }
}


void DynamicProfilePanel::save()
{
    std::vector<DynamicProfileEntry> entries;
    int serial = 1;
    for (auto row : treemodel_->children()) {
        entries.emplace_back(to_entry(row, serial++));
    }
    if (!storeDynamicProfileEntries(entries)) {
        printf("Error in saving dynamic profile rules\n");
    } else {
        printf("Saved %d dynamic profile rules\n", int(entries.size()));
    }
}
