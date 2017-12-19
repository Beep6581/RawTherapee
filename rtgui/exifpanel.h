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
#ifndef _EXIFPANEL_
#define _EXIFPANEL_

#include <gtkmm.h>
#include "toolpanel.h"

class ExifPanel : public Gtk::VBox, public ToolPanel
{

private:
    const rtengine::FramesMetaData* idata;
    rtengine::procparams::ExifPairs changeList;
    rtengine::procparams::ExifPairs defChangeList;
    bool recursiveOp;

    class ExifColumns : public Gtk::TreeModelColumnRecord
    {
    public:
        Gtk::TreeModelColumn<Glib::RefPtr<Gdk::Pixbuf> > icon;
        Gtk::TreeModelColumn<Glib::ustring> field;
        Gtk::TreeModelColumn<Glib::ustring> field_nopango;
        Gtk::TreeModelColumn<Glib::ustring> value;
        Gtk::TreeModelColumn<Glib::ustring> value_nopango;
        Gtk::TreeModelColumn<Glib::ustring> orig_value;
        Gtk::TreeModelColumn<rtexif::ActionCode> action;
        Gtk::TreeModelColumn<bool> editable;
        Gtk::TreeModelColumn<bool> edited;
        Gtk::TreeModelColumn<bool> isSeparator;

        ExifColumns()
        {
            add (field);
            add (value);
            add (icon);
            add (action);
            add (edited);
            add (field_nopango);
            add (value_nopango);
            add (editable);
            add (orig_value);
            add (isSeparator);
        }
    };
    Glib::RefPtr<Gdk::Pixbuf> delicon;
    Glib::RefPtr<Gdk::Pixbuf> keepicon;
    Glib::RefPtr<Gdk::Pixbuf> editicon;

    ExifColumns exifColumns;
    Gtk::TreeView* exifTree;
    Gtk::ScrolledWindow* scrolledWindow;
    Glib::RefPtr<Gtk::TreeStore> exifTreeModel;

    Gtk::Button* remove;
    Gtk::Button* keep;
    Gtk::Button* add;
    Gtk::Button* reset;
    Gtk::Button* resetAll;
    Gtk::ToggleButton* showAll;

    Gtk::TreeModel::Children addTag (const Gtk::TreeModel::Children& root, Glib::ustring field, Glib::ustring value, rtexif::ActionCode action, bool editable);
    void editTag (Gtk::TreeModel::Children root, Glib::ustring name, Glib::ustring value);
    void updateChangeList (Gtk::TreeModel::Children root, std::string prefix);
    void addDirectory (const rtexif::TagDirectory* dir, Gtk::TreeModel::Children root, bool checkForSeparator = false);
    Gtk::TreeModel::Children addSeparator();
    Glib::ustring getSelection (bool onlyifeditable = false);
    Glib::ustring getSelectedValue();
    void updateChangeList();
    void applyChangeList();
    void keepIt (Gtk::TreeModel::iterator iter);
    void delIt (Gtk::TreeModel::iterator iter);
    Gtk::TreeModel::iterator resetIt (Gtk::TreeModel::iterator iter);
    void removePressed();
    void keepPressed();
    void resetPressed();
    void resetAllPressed();
    void addPressed();
    void showAlltoggled();
    bool rowSeperatorFunc(const Glib::RefPtr<Gtk::TreeModel>& model, const Gtk::TreeModel::iterator& iter);

public:
    ExifPanel ();
    virtual ~ExifPanel();

    void read (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);

    void setImageData (const rtengine::FramesMetaData* id);

    void exifSelectionChanged();
    void row_activated (const Gtk::TreeModel::Path& path, Gtk::TreeViewColumn* column);

    void notifyListener();

};

#endif
