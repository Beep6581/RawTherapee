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
#ifndef _RENAMEDLG_
#define _RENAMEDLG_

#include <gtkmm.h>
#include "cacheimagedata.h"
#include "guiutils.h"

#define RESPONSE_ALL 100

class RenameDialog : public Gtk::Dialog
{

protected:

    class TemplateColumns : public Gtk::TreeModel::ColumnRecord
    {
    public:
        Gtk::TreeModelColumn<Glib::ustring> tmplName;
        Gtk::TreeModelColumn<bool>          rowSeparator;
        TemplateColumns()
        {
            add(tmplName);
            add(rowSeparator);
        }
    };
    TemplateColumns         templateColumns;
    Glib::RefPtr<Gtk::ListStore> templateModel;

    Gtk::Window* p;
    Gtk::Label* oldName;
    Gtk::Entry* newName;
    Gtk::CheckButton* useTmpl;
    MyComboBox* templates;
    Gtk::Button* all;
    const CacheImageData* imageData;

    void fillTemplateList ();

public:
    RenameDialog (Gtk::Window* parent);

    void initName (const Glib::ustring& iname, const CacheImageData* cid);
    Glib::ustring getNewName ();

    bool rowSeparatorFunc (const Glib::RefPtr<Gtk::TreeModel>& model, const Gtk::TreeModel::iterator& iter);
    void tmplSelectionChanged ();
    void useTemplToggled ();

    bool isTemplSelected ();
    Glib::ustring getActiveTemplate ();

    static Glib::ustring applyTemplate (const Glib::ustring& oName, const CacheImageData* cid, const Glib::ustring& templ);
};

class RenameTemplateEditor : public Gtk::Dialog
{

protected:
    Gtk::ListViewText* list;
    Gtk::Entry* templ;

    void refreshTemplateList ();
public:
    RenameTemplateEditor (Gtk::Window* parent);

    Glib::ustring getSelectedTemplate ();

    void addPressed ();
    void delPressed ();
};

#endif

