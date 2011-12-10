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
#include <toolpanel.h>
#include <imagedata.h>

class ExifPanel : public Gtk::VBox, public ToolPanel {

    private: 
        const rtengine::ImageMetaData* idata;
        
        class ExifColumns : public Gtk::TreeModelColumnRecord {
            public:
                Gtk::TreeModelColumn<Glib::RefPtr<Gdk::Pixbuf> > icon;
                Gtk::TreeModelColumn<Glib::ustring> field;
                Gtk::TreeModelColumn<Glib::ustring> value;
                Gtk::TreeModelColumn<Glib::ustring> orig_value;

                ExifColumns() { add(field); add(value); add(icon); add(orig_value); }
        };
        
        ExifColumns exifColumns;
        Gtk::TreeView* exifTree;
        Gtk::ScrolledWindow* scrolledWindow;
        Glib::RefPtr<Gtk::TreeStore> exifTreeModel;
        
        Gtk::TreeModel::Children addTag  (const Gtk::TreeModel::Children& root, Glib::ustring field, Glib::ustring value );
        Gtk::TreeModel::Children addGroup (const Gtk::TreeModel::Children& root, Glib::ustring group );
     public:
        ExifPanel ();
        virtual ~ExifPanel ();
        
		void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
		void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
		void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
        
        void setImageData (const rtengine::ImageMetaData* id );

        void notifyListener ();
        
};

#endif
