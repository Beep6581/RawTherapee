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

#include "rtengine/imagedata.h"
#include "rtengine/metadata.h"
#include "rtengine/procparams.h"

#include <glibmm/regex.h>

using namespace rtengine;
using namespace rtengine::procparams;

ExifPanel::ExifPanel() :
    idata(nullptr),
    changeList(new rtengine::procparams::ExifPairs),
    defChangeList(new rtengine::procparams::ExifPairs),
    //keepicon("tick-small"),
    editicon("edit-small"),
    open_icon_("expander-open-small"),
    closed_icon_("expander-closed-small"),
    pl_(nullptr)
{
    for (auto &k : MetaDataParams::basicExifKeys) {
        editableTags.push_back(std::make_pair(k, ""));
    }

    set_orientation(Gtk::ORIENTATION_VERTICAL);
    exifTree = Gtk::manage (new Gtk::TreeView());
    scrolledWindow = Gtk::manage (new Gtk::ScrolledWindow());

    exifTree->set_headers_visible (false);
    exifTree->set_rules_hint (false);
    exifTree->set_reorderable (false);
    exifTree->set_enable_search (false);
    exifTree->get_selection()->set_mode(Gtk::SELECTION_SINGLE);
    scrolledWindow->set_shadow_type(Gtk::SHADOW_NONE);
    scrolledWindow->set_policy(Gtk::POLICY_ALWAYS, Gtk::POLICY_ALWAYS);
    scrolledWindow->property_window_placement().set_value(Gtk::CORNER_TOP_LEFT);
    scrolledWindow->add(*exifTree);

    exifTreeModel = Gtk::TreeStore::create(exifColumns);
    exifTree->set_model(exifTreeModel);
    exifTree->set_grid_lines(Gtk::TREE_VIEW_GRID_LINES_NONE);
    exifTree->set_show_expanders(false);
    exifTree->set_tooltip_column(0);
    exifTree->set_enable_search(false);

    exif_active_renderer_.property_mode() = Gtk::CELL_RENDERER_MODE_ACTIVATABLE;
    exif_active_renderer_.signal_toggled().connect(sigc::mem_fun(this, &ExifPanel::onKeyActiveToggled));
    exif_active_column_.pack_start(exif_active_renderer_);
    exif_active_column_.set_cell_data_func(exif_active_renderer_, sigc::mem_fun(this, &ExifPanel::setKeyActive));

    exifTree->append_column(exif_active_column_);

    Gtk::TreeView::Column *viewcol = Gtk::manage (new Gtk::TreeView::Column ("Field Name"));
    Gtk::CellRendererPixbuf* render_pb = Gtk::manage (new Gtk::CellRendererPixbuf());
    render_pb->property_stock_size() = Gtk::ICON_SIZE_SMALL_TOOLBAR;
    Gtk::CellRendererText *render_txt = Gtk::manage (new Gtk::CellRendererText());
    render_txt->property_ellipsize() = Pango::ELLIPSIZE_END;
    viewcol->pack_start(*render_pb, false);
    viewcol->pack_start(*render_txt, true);
    viewcol->add_attribute(*render_pb, "icon-name", exifColumns.icon);
    viewcol->add_attribute(*render_txt, "markup", exifColumns.label);
    viewcol->set_expand(true);
    viewcol->set_resizable(true);
    viewcol->set_fixed_width(35);
    viewcol->set_min_width(35);
    viewcol->set_sizing(Gtk::TREE_VIEW_COLUMN_AUTOSIZE);

    render_pb->property_ypad() = 0;
    render_txt->property_ypad() = 0;
    render_pb->property_yalign() = 0;
    render_txt->property_yalign() = 0;

    exifTree->append_column(*viewcol);
    //exifTree->set_expander_column(*viewcol);

    Gtk::TreeView::Column *viewcolv = Gtk::manage(new Gtk::TreeView::Column ("Value"));
    Gtk::CellRendererText *render_txtv = Gtk::manage(new Gtk::CellRendererText());
    render_txtv->property_ellipsize() = Pango::ELLIPSIZE_END;
    viewcolv->pack_start (*render_txtv, true);
    // viewcolv->add_attribute (*render_txtv, "markup", exifColumns.value);
    viewcolv->set_expand (true);
    viewcolv->set_resizable (true);
    viewcol->set_fixed_width (35);
    viewcolv->set_min_width (35);
    viewcolv->set_sizing (Gtk::TREE_VIEW_COLUMN_AUTOSIZE);

    render_txtv->property_ypad() = 0;
    viewcolv->set_cell_data_func(*render_txtv, sigc::mem_fun(this, &ExifPanel::setExifTagValue));
    render_txtv->signal_edited().connect(sigc::mem_fun(this, &ExifPanel::onEditExifTagValue));

    exifTree->append_column(*viewcolv);

    pack_start (*scrolledWindow);

    Gtk::Grid* buttons1 = Gtk::manage (new Gtk::Grid());
    buttons1->set_row_homogeneous (true);
    buttons1->set_column_homogeneous (true);
    setExpandAlignProperties (buttons1, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    const auto addbtn =
        [&](const Glib::ustring &tip, const Glib::ustring &icon) -> Gtk::Button *
        {
            Gtk::Button *b = Gtk::manage(new Gtk::Button());
            b->set_image(*Gtk::manage(new RTImage(icon, Gtk::ICON_SIZE_BUTTON)));
            b->set_tooltip_text(M(tip));
            b->get_style_context()->add_class("Right");
            setExpandAlignProperties(b, true, true, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
            buttons1->attach_next_to(*b, Gtk::POS_RIGHT, 1, 1);
            return b;
        };

    activate_all_ = addbtn("EXIFPANEL_ACTIVATE_ALL_HINT", "tick");
    activate_none_ = addbtn("EXIFPANEL_ACTIVATE_NONE_HINT", "box");
    add = addbtn("EXIFPANEL_ADDEDIT", "edit");
    reset = addbtn("EXIFPANEL_RESETHINT", "undo");
    resetAll = addbtn("EXIFPANEL_RESETALLHINT", "undo-all");

    pack_end (*buttons1, Gtk::PACK_SHRINK);

    exifTree->get_selection()->signal_changed().connect (sigc::mem_fun (*this, &ExifPanel::exifSelectionChanged));

    reset->signal_clicked().connect ( sigc::mem_fun (*this, &ExifPanel::resetPressed) );
    resetAll->signal_clicked().connect ( sigc::mem_fun (*this, &ExifPanel::resetAllPressed) );
    add->signal_clicked().connect ( sigc::mem_fun (*this, &ExifPanel::addPressed) );
    activate_all_->signal_clicked().connect(sigc::mem_fun(*this, &ExifPanel::activateAllPressed));
    activate_none_->signal_clicked().connect(sigc::mem_fun(*this, &ExifPanel::activateNonePressed));

    exifTree->signal_button_press_event().connect_notify(sigc::mem_fun(*this, &ExifPanel::onExifTreeClick));
    exifTree->signal_row_expanded().connect(sigc::mem_fun(*this, &ExifPanel::onExifRowExpanded));
    exifTree->signal_row_collapsed().connect(sigc::mem_fun(*this, &ExifPanel::onExifRowCollapsed));

    show_all ();
}

ExifPanel::~ExifPanel ()
{
}

void ExifPanel::read (const ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener();

    *changeList = pp->metadata.exif;
    initial_active_keys_.clear();
    initial_active_keys_.insert(pp->metadata.exifKeys.begin(), pp->metadata.exifKeys.end());
    cur_active_keys_ = initial_active_keys_;
    setImageData(idata);
    refreshTags();

    enableListener();
}


void ExifPanel::write (ProcParams* pp, ParamsEdited* pedited)
{
    pp->metadata.exif = *changeList;
    cur_active_keys_ = get_active_keys();
    pp->metadata.exifKeys.assign(cur_active_keys_.begin(), cur_active_keys_.end());
    std::sort(pp->metadata.exifKeys.begin(), pp->metadata.exifKeys.end());
}

void ExifPanel::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{
    *defChangeList = defParams->metadata.exif;
}

void ExifPanel::setImageData (const FramesMetaData* id)
{
    idata = id;
}

void ExifPanel::addTag(const std::string &key, const std::pair<Glib::ustring, Glib::ustring> &label, const Glib::ustring &exifValue, bool editable, bool edited)
{

    const Glib::ustring& value = exifValue.validate() ? exifValue : "???";

//    auto root = exifTreeModel->children();

    const auto getgroup =
        [&]() -> Gtk::TreeNodeChildren
        {
            auto root = exifTreeModel->children();

            // for (auto it = root.rbegin(), end = root.rend(); it != end; ++it) {
            //     auto row = *it;
            for (auto &row : root) {
                // auto row = *it;
                std::string key = row[exifColumns.key];
                if (row[exifColumns.is_group] && key == label.first) {
                    return row./*it->*/children();
                }
            }
            auto it = exifTreeModel->append(root);
            auto row = *it;

            row[exifColumns.editable] = false;
            row[exifColumns.edited] = false;
            row[exifColumns.key] = label.first;
            row[exifColumns.label] = "<b>" + label.first + "</b>";
            row[exifColumns.value_nopango] = "";
            row[exifColumns.value] = "";
            row[exifColumns.is_group] = true;

            return it->children();
        };

    auto root = getgroup();

    Gtk::TreeModel::Row row = *(exifTreeModel->append(root));
    row[exifColumns.editable] = editable;
    row[exifColumns.edited] = edited;
    row[exifColumns.key] = key;
    row[exifColumns.is_group] = false;
    //row[exifColumns.label] = label.second;
    row[exifColumns.value_nopango] = value;
    //row[exifColumns.value] = value;

    row[exifColumns.label] = escapeHtmlChars(label.second);
    row[exifColumns.value] = value;//escapeHtmlChars(value);

    bool active = all_keys_active() || cur_active_keys_.find(key) != cur_active_keys_.end();
    row[exifColumns.active] = active;

    if (edited) {
        row[exifColumns.icon] = editicon;
    // } else if (editable) {
    //     row[exifColumns.icon] = keepicon;
    }
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

    const auto to_label =
        [](const Exiv2::Exifdatum &tag) -> std::pair<Glib::ustring, Glib::ustring>
        {
            auto s = tag.key();
            auto pos = s.find('.');
            if (pos != std::string::npos) {
                s = s.substr(pos+1);
            }
            Glib::ustring g = "";
            pos = s.find('.');
            if (pos != std::string::npos) {
                g = s.substr(0, pos);
                s = s.substr(pos+1);
            }
            auto l = tag.tagLabel();
            if (!l.empty()) {
                s = l;
            }
            return std::make_pair(g, Glib::ustring(s));
        };

    try {
        rtengine::Exiv2Metadata meta(fn);
        meta.load();
        Exiv2::ExifData exif = meta.exifData();

        const auto to_value =
            [&](Exiv2::Exifdatum &tag) -> Glib::ustring
            {
                try {
                    if (!tag.tagLabel().empty() && //tag.typeId() != Exiv2::undefined &&
                        (tag.typeId() == Exiv2::asciiString || tag.size() < 256)) {
                        return escapeHtmlChars(tag.print(&exif));
                    }
                } catch (std::exception &) {}
                return "<i>(" + M("EXIFPANEL_VALUE_NOT_SHOWN") + ")</i>";
            };

        if (const rtengine::FramesData *fd = dynamic_cast<const rtengine::FramesData *>(idata)) {
            fd->fillBasicTags(exif);
        }

        for (const auto& p : *changeList) {
            try {
                auto &datum = exif[p.first];
                if (datum.setValue(p.second) != 0) {
                    if (pl_) {
                        pl_->error(Glib::ustring::compose(M("ERROR_MSG_METADATA_VALUE"), p.first, p.second));
                    }
                }
            } catch (const std::exception& exc) {
                if (pl_) {
                    pl_->error(Glib::ustring::compose(M("ERROR_MSG_METADATA_VALUE"), p.first, p.second));
                }
            }
        }

        for (auto& p : editableTags) {
            Exiv2::ExifKey k(p.first);
            const auto pos = exif.findKey(k);
            bool edited = false;
            Glib::ustring value = "";
            auto lbl = std::make_pair(M("EXIFPANEL_BASIC_GROUP"), k.tagLabel());
            p.second = k.tagLabel();
            if (pos != exif.end() && pos->size()) {
                edited = changeList->find(pos->key()) != changeList->end();
                value = escapeHtmlChars(pos->print(&exif));
            }
            addTag(p.first, lbl, value, true, edited);
        }
        struct KeyLt {
            KeyLt():
                order_{
                    {"Exif.GPSInfo", 0},
                    {"Exif.Photo", 1},
                    {"Exif.Image", 2}
                }
            {}
            bool operator()(const std::string &a, const std::string &b) const
            {
                auto p1 = a.find_last_of('.');
                auto p2 = b.find_last_of('.');
                const char *sa = a.c_str();
                const char *sb = b.c_str();
                if (p1 != std::string::npos && p2 != std::string::npos) {
                    bool hex_a = strncmp(sa+p1+1, "0x", 2) == 0;
                    bool hex_b = strncmp(sb+p2+1, "0x", 2) == 0;
                    if (hex_a != hex_b) {
                        return !hex_a;
                    }
                }
                if (p1 != p2 || strncmp(sa, sb, p1) != 0) {
                    std::string ga(sa, sa+p1);
                    std::string gb(sb, sb+p2);
                    int ia = getorder(ga);
                    int ib = getorder(gb);
                    if (ia != ib) {
                        return ia < ib;
                    }
                }
                return strcmp(sa, sb) < 0;
            }

            int getorder(const std::string &key) const
            {
                auto it = order_.find(key);
                if (it == order_.end()) {
                    return 1000;
                }
                return it->second;
            }

            std::unordered_map<std::string, int> order_;
        };
        std::set<std::string, KeyLt> keyset;
        for (const auto& tag : exif) {
            const bool editable = ed.find(tag.key()) != ed.end();
            if (!editable) {
                keyset.insert(tag.key());
            }
        }
        for (auto &k : keyset) {
            auto &tag = *(exif.findKey(Exiv2::ExifKey(k)));
            addTag(tag.key(), to_label(tag), to_value(tag), false, false);
        }
    } catch (const std::exception& exc) {
        return;
    }

    for (const auto& p : sel) {
        exifTree->get_selection()->select(p);
    }

    exifTree->expand_all();
}

void ExifPanel::exifSelectionChanged ()
{
    Glib::RefPtr<Gtk::TreeSelection> selection = exifTree->get_selection();
    std::vector<Gtk::TreeModel::Path> sel = selection->get_selected_rows();

    if (sel.size() >= 1) {
        Gtk::TreeModel::iterator iter = exifTreeModel->get_iter (sel[0]);
        if (iter->get_value(exifColumns.editable)) {
            reset->set_sensitive(true);
            add->set_sensitive(true);
        }
    } else {
        reset->set_sensitive(false);
        add->set_sensitive(false);
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

void ExifPanel::resetPressed()
{
    cur_active_keys_ = get_active_keys();

    auto sel = exifTree->get_selection()->get_selected_rows();
    for (size_t i = 0; i < sel.size(); i++) {
        resetIt(exifTreeModel->get_iter(sel[i]));
    }

    refreshTags();

    if (!sel.empty()) {
        exifTree->get_selection()->select(sel[0]);
    }

    notifyListener();
}

void ExifPanel::resetAllPressed()
{
    auto sel = exifTree->get_selection()->get_selected_rows();
    setImageData(idata);
    *changeList = *defChangeList;
    cur_active_keys_ = initial_active_keys_;
    refreshTags();
    if (!sel.empty()) {
        exifTree->get_selection()->select(sel[0]);
    }
    notifyListener();
}

void ExifPanel::addPressed()
{
    Gtk::TreeModel::Path path;
    Gtk::TreeViewColumn *col;

    exifTree->get_cursor(path, col);
    const auto it = exifTreeModel->get_iter(path);
    if (it && it->get_value(exifColumns.editable)) {
        exifTree->set_cursor(path, *col, true);
    }
}

void ExifPanel::activateAllPressed()
{
    disableListener();
    auto root = exifTreeModel->children();
    for (auto &group : root->children()) {
        group[exifColumns.active] = true;
        for (auto &row : group.children()) {
            row[exifColumns.active] = true;
        }
    }
    enableListener();
    notifyListener();
}


void ExifPanel::activateNonePressed()
{
    disableListener();
    auto root = exifTreeModel->children();
    for (auto &group : root->children()) {
        group[exifColumns.active] = false;
        for (auto &row : group.children()) {
            row[exifColumns.active] = false;
        }
    }
    enableListener();
    notifyListener();
}


void ExifPanel::notifyListener()
{
    if (listener) {
        listener->panelChanged(EvExif, M("HISTORY_CHANGED"));
    }
}


void ExifPanel::onKeyActiveToggled(const Glib::ustring &path)
{
    auto it = exifTreeModel->get_iter(path);
    if (it) {
        auto row = *it;
        bool b = !row[exifColumns.active];
        row[exifColumns.active] = b;
        if (row[exifColumns.is_group]) {
            for (auto &c : row.children()) {
                c[exifColumns.active] = b;
            }
        }
        notifyListener();
    }
}


void ExifPanel::setKeyActive(Gtk::CellRenderer *renderer, const Gtk::TreeModel::iterator &it)
{
    auto row = *it;
    Gtk::CellRendererToggle *toggle = static_cast<Gtk::CellRendererToggle *>(renderer);
    if (row[exifColumns.is_group]) {
        bool all_true = true, all_false = true;
        for (auto &c : row.children()) {
            if (c[exifColumns.active]) {
                all_false = false;
            } else {
                all_true = false;
            }
        }
        if (!all_true && !all_false) {
            toggle->property_inconsistent() = true;
        } else {
            toggle->property_inconsistent() = false;
            toggle->set_active(all_true);
        }
    } else {
        toggle->property_inconsistent() = false;
        toggle->set_active(row[exifColumns.active]);
    }
}


bool ExifPanel::all_keys_active() const
{
    return (cur_active_keys_.size() == 1 && *(cur_active_keys_.begin()) == "*");
}


std::unordered_set<std::string> ExifPanel::get_active_keys() const
{
    bool all_active = true;
    std::unordered_set<std::string> ret;
    auto root = exifTreeModel->children();
    for (auto &group : root->children()) {
        for (auto &entry : group.children()) {
            std::string key = entry[exifColumns.key];
            if (entry[exifColumns.active]) {
                ret.insert(key);
            } else {
                all_active = false;
            }
        }
    }
    if (all_active) {
        ret.clear();
        ret.insert("*");
    }
    return ret;
}

void ExifPanel::onExifTreeClick(GdkEventButton *event)
{
    Gtk::TreeModel::Path pth;
    Gtk::TreeViewColumn *col;
    int cell_x;
    int cell_y;
    if (exifTree->get_path_at_pos(event->x, event->y, pth, col, cell_x, cell_y) && col == exifTree->get_column(1) && cell_x <= 22) {
        auto it = exifTreeModel->get_iter(pth);
        auto row = *it;
        if (row[exifColumns.is_group]) {
            if (exifTree->row_expanded(pth)) {
                exifTree->collapse_row(pth);
            } else {
                exifTree->expand_row(pth, false);
            }
        }
    }
}


void ExifPanel::onExifRowExpanded(const Gtk::TreeModel::iterator &it, const Gtk::TreeModel::Path &path)
{
    auto row = *it;
    if (row[exifColumns.is_group]) {
        row[exifColumns.icon] = open_icon_;
    }
}


void ExifPanel::onExifRowCollapsed(const Gtk::TreeModel::iterator &it, const Gtk::TreeModel::Path &path)
{
    auto row = *it;
    if (row[exifColumns.is_group]) {
        row[exifColumns.icon] = closed_icon_;
    }
}


void ExifPanel::setExifTagValue(Gtk::CellRenderer *renderer, const Gtk::TreeModel::iterator &it)
{
    auto row = *it;
    Gtk::CellRendererText *txt = static_cast<Gtk::CellRendererText *>(renderer);
    txt->property_editable() = row[exifColumns.editable];
    txt->property_markup() = row[exifColumns.value];
}


namespace {


Glib::ustring to_fraction(const Glib::ustring &s)
{
    auto i = s.find(".");
    if (i != Glib::ustring::npos) {
        return s.substr(0, i) + s.substr(i+1) + "/1" + Glib::ustring(s.size() - i - 1, '0');
    }
    return s;
}

typedef Glib::ustring (*validator_func)(const Glib::ustring &);

Glib::ustring get_fnumber(const Glib::ustring &val)
{
    Glib::MatchInfo m;
    auto re = Glib::Regex::create("f? *([0-9.]+) *", Glib::REGEX_CASELESS);
    if (re->match(val, m)) {
        auto s = m.fetch(1);
        return to_fraction(s);
    }
    return val;
}

Glib::ustring get_shutterspeed(const Glib::ustring &val)
{
    Glib::MatchInfo m;
    auto re = Glib::Regex::create(" *([0-9/]+) *s? *", Glib::REGEX_CASELESS);
    if (re->match(val, m)) {
        auto s = m.fetch(1);
        return s;
    }
    return val;
}

Glib::ustring get_focallen(const Glib::ustring &val)
{
    Glib::MatchInfo m;
    auto re = Glib::Regex::create(" *([0-9.]+) *(mm)? *", Glib::REGEX_CASELESS);
    if (re->match(val, m)) {
        auto s = m.fetch(1);
        return to_fraction(s);
    }
    return val;
}

Glib::ustring get_expcomp(const Glib::ustring &val)
{
    Glib::MatchInfo m;
    auto re = Glib::Regex::create(" *(-?[0-9.]+) *(EV)? *", Glib::REGEX_CASELESS);
    if (re->match(val, m)) {
        auto s = m.fetch(1);
        return to_fraction(s);
    }
    return val;
}

std::unordered_map<std::string, validator_func> validators = {
    {"Exif.Photo.FNumber", get_fnumber},
    {"Exif.Photo.ExposureTime", get_shutterspeed},
    {"Exif.Photo.FocalLength", get_focallen},
    {"Exif.Photo.ExposureBiasValue", get_expcomp}
};

} // namespace


void ExifPanel::onEditExifTagValue(const Glib::ustring &path, const Glib::ustring &val)
{
    auto it = exifTreeModel->get_iter(path);
    auto row = *it;
    std::string key = row[exifColumns.key];
    auto value = val;

    bool good = true;
    try {
        Exiv2::ExifData data;
        auto &datum = data[key];
        auto it = validators.find(key);
        if (it != validators.end()) {
            auto v = it->second(value);
            if (datum.setValue(v) == 0) {
                value = v;
            }
        }
        if (datum.setValue(value) != 0) {
            if ((datum.typeId() == Exiv2::signedRational || datum.typeId() == Exiv2::unsignedRational) && datum.setValue(value + "/1") == 0) {
                value += "/1";
            } else {
                good = false;
            }
        }
    } catch (std::exception &exc) {
        good = false;
    }

    if (good) {
        (*changeList)[key] = value;
        if (!all_keys_active()) {
            cur_active_keys_.insert(key);
        }
        refreshTags();

        it = exifTreeModel->get_iter(path);
        exifTree->get_selection()->select(it);
        notifyListener();
    } else if (pl_) {
        pl_->error(Glib::ustring::compose(M("ERROR_MSG_METADATA_VALUE"), key, value));
    }
}


void ExifPanel::setProgressListener(rtengine::ProgressListener *pl)
{
    pl_ = pl;
}
