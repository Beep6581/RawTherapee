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
#ifndef _DYNAMICPROFILEPANEL_H_
#define _DYNAMICPROFILEPANEL_H_

#include <gtkmm.h>
#include "../rtengine/dynamicprofile.h"
#include "profilestorecombobox.h"

class DynamicProfilePanel: public Gtk::VBox
{
public:
    DynamicProfilePanel();
    void save();

private:
    void update_rule (Gtk::TreeModel::Row row, const DynamicProfileRule &rule);
    void add_rule (const DynamicProfileRule &rule);
    DynamicProfileRule to_rule (Gtk::TreeModel::Row row, int serial = 0);

    void on_button_quit();
    void on_button_up();
    void on_button_down();
    void on_button_new();
    void on_button_edit();
    void on_button_delete();

    class DynamicProfileColumns: public Gtk::TreeModel::ColumnRecord
    {
    public:
        DynamicProfileColumns()
        {
            add (iso);
            add (fnumber);
            add (focallen);
            add (shutterspeed);
            add (expcomp);
            add (camera);
            add (lens);
            add (profilepath);
        }

        Gtk::TreeModelColumn<DynamicProfileRule::Range<int>> iso;
        Gtk::TreeModelColumn<DynamicProfileRule::Range<double>> fnumber;
        Gtk::TreeModelColumn<DynamicProfileRule::Range<double>> focallen;
        Gtk::TreeModelColumn<DynamicProfileRule::Range<double>> shutterspeed;
        Gtk::TreeModelColumn<DynamicProfileRule::Range<double>> expcomp;
        Gtk::TreeModelColumn<DynamicProfileRule::Optional> camera;
        Gtk::TreeModelColumn<DynamicProfileRule::Optional> lens;
        Gtk::TreeModelColumn<Glib::ustring> profilepath;
    };

    // cell renderers
    void render_iso (Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);
    void render_fnumber (Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);
    void render_focallen (Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);
    void render_shutterspeed (Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);
    void render_expcomp (Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);
    void render_camera (Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);
    void render_lens (Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);
    void render_profilepath (Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);

    class EditDialog: public Gtk::Dialog
    {
    public:
        EditDialog (const Glib::ustring &title, Gtk::Window &parent);
        void set_rule (const DynamicProfileRule &rule);
        DynamicProfileRule get_rule();

    private:
        void set_ranges();
        void add_range (const Glib::ustring &name, Gtk::SpinButton *&from, Gtk::SpinButton *&to);
        void add_optional (const Glib::ustring &name, Gtk::CheckButton *&check, Gtk::Entry *&field);

        Gtk::SpinButton *iso_min_;
        Gtk::SpinButton *iso_max_;

        Gtk::SpinButton *fnumber_min_;
        Gtk::SpinButton *fnumber_max_;

        Gtk::SpinButton *focallen_min_;
        Gtk::SpinButton *focallen_max_;

        Gtk::SpinButton *shutterspeed_min_;
        Gtk::SpinButton *shutterspeed_max_;

        Gtk::SpinButton *expcomp_min_;
        Gtk::SpinButton *expcomp_max_;

        Gtk::CheckButton *has_camera_;
        Gtk::Entry *camera_;

        Gtk::CheckButton *has_lens_;
        Gtk::Entry *lens_;

        ProfileStoreComboBox *profilepath_;
    };

    DynamicProfileColumns columns_;

    //Child widgets:
    Gtk::Box vbox_;

    Gtk::ScrolledWindow scrolledwindow_;
    Gtk::TreeView treeview_;
    Glib::RefPtr<Gtk::ListStore> treemodel_;

    Gtk::ButtonBox buttonbox_;
    Gtk::Button button_up_;
    Gtk::Button button_down_;
    Gtk::Button button_new_;
    Gtk::Button button_edit_;
    Gtk::Button button_delete_;
};

#endif // _DYNAMICPROFILEPANEL_H_
