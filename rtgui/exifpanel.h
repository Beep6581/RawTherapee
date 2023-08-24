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
#pragma once

#include <memory>

#include <gtkmm.h>
#include <unordered_set>

#include "toolpanel.h"

namespace rtengine
{

namespace procparams
{

class ExifPairs;

}

}

class ExifPanel final :
    public Gtk::Box,
    public ToolPanel
{

private:
    const rtengine::FramesMetaData* idata;
    const std::unique_ptr<rtengine::procparams::ExifPairs> changeList;
    const std::unique_ptr<rtengine::procparams::ExifPairs> defChangeList;

    class ExifColumns : public Gtk::TreeModelColumnRecord
    {
    public:
        Gtk::TreeModelColumn<Glib::ustring> icon;
        // Gtk::TreeModelColumn<Glib::RefPtr<Gdk::Pixbuf>> expander_icon;
        Gtk::TreeModelColumn<std::string> key;
        Gtk::TreeModelColumn<Glib::ustring> label;
        Gtk::TreeModelColumn<Glib::ustring> value;
        Gtk::TreeModelColumn<Glib::ustring> value_nopango;
        Gtk::TreeModelColumn<bool> editable;
        Gtk::TreeModelColumn<bool> edited;
        Gtk::TreeModelColumn<bool> active;
        Gtk::TreeModelColumn<bool> is_group;

        ExifColumns()
        {
            add(key);
            add(label);
            add(value);
            add(icon);
            add(edited);
            add(value_nopango);
            add(editable);
            add(active);
            add(is_group);
            // add(expander_icon);
        }
    };

    //Glib::ustring keepicon;
    Glib::ustring editicon;
    Glib::ustring open_icon_;
    Glib::ustring closed_icon_;

    ExifColumns exifColumns;
    Gtk::TreeView* exifTree;
    Gtk::ScrolledWindow* scrolledWindow;
    Glib::RefPtr<Gtk::TreeStore> exifTreeModel;

    Gtk::Button* add;
    Gtk::Button* reset;
    Gtk::Button* resetAll;
    Gtk::Button *activate_all_;
    Gtk::Button *activate_none_;

    Gtk::CellRendererToggle exif_active_renderer_;
    Gtk::TreeView::Column exif_active_column_;

    std::vector<std::pair<std::string, Glib::ustring>> editableTags;

    std::unordered_set<std::string> initial_active_keys_;
    std::unordered_set<std::string> cur_active_keys_;

    rtengine::ProgressListener *pl_;

    void addTag(const std::string &key, const std::pair<Glib::ustring, Glib::ustring> &label, const Glib::ustring &value, bool editable, bool edited);
    void refreshTags();
    void resetIt(const Gtk::TreeModel::const_iterator& iter);
    void resetPressed();
    void resetAllPressed();
    void addPressed();
    void activateAllPressed();
    void activateNonePressed();

    void setKeyActive(Gtk::CellRenderer *renderer, const Gtk::TreeModel::iterator &it);
    void onKeyActiveToggled(const Glib::ustring &path);

    bool all_keys_active() const;
    std::unordered_set<std::string> get_active_keys() const;

    void onExifTreeClick(GdkEventButton *event);
    void onExifRowExpanded(const Gtk::TreeModel::iterator &it, const Gtk::TreeModel::Path &path);
    void onExifRowCollapsed(const Gtk::TreeModel::iterator &it, const Gtk::TreeModel::Path &path);

    void setExifTagValue(Gtk::CellRenderer *renderer, const Gtk::TreeModel::iterator &it);
    void onEditExifTagValue(const Glib::ustring &path, const Glib::ustring &value);

public:
    ExifPanel ();
    ~ExifPanel() override;

    void read (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;

    void setImageData (const rtengine::FramesMetaData* id);

    void exifSelectionChanged();
    // void row_activated (const Gtk::TreeModel::Path& path, Gtk::TreeViewColumn* column);

    void notifyListener();

    void setProgressListener(rtengine::ProgressListener *pl);
};
