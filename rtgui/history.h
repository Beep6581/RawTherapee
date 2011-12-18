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
#ifndef _HISTORY_
#define _HISTORY_

#include <gtkmm.h>
#include <rtengine.h>
#include <pparamschangelistener.h>
#include <profilechangelistener.h>
#include <paramsedited.h>
#include <thumbnail.h>
#include <imagedata.h>

class HistoryBeforeLineListener {

    public:
        virtual void historyBeforeLineChanged (const rtengine::procparams::ProcParams& params) {}
};



class History : public Gtk::VBox, public PParamsChangeListener {

    public:

        class HistoryColumns : public Gtk::TreeModel::ColumnRecord {
            public:
                Gtk::TreeModelColumn<Glib::ustring>  realText;
                Gtk::TreeModelColumn<Glib::ustring>  text;
                Gtk::TreeModelColumn<Glib::ustring>  value;
                Gtk::TreeModelColumn<rtengine::procparams::ProcParams>     params;
                Gtk::TreeModelColumn<rtengine::ProcEvent>    chev;
                Gtk::TreeModelColumn<ParamsEdited>     paramsEdited;
                HistoryColumns() { add(text); add(realText); add(value); add(chev); add(params); add(paramsEdited); }
        };
        HistoryColumns historyColumns;
        class BookmarkColumns : public Gtk::TreeModel::ColumnRecord {
            public:
                //Gtk::TreeModelColumn<Glib::ustring>  text;
        	    Gtk::TreeModelColumn<int> id;
                Gtk::TreeModelColumn<rtengine::procparams::ProcParams>     params;
                Gtk::TreeModelColumn<ParamsEdited>     paramsEdited;
                Gtk::TreeModelColumn<Glib::ustring> name;
                Gtk::TreeModelColumn< Glib::RefPtr<Gdk::Pixbuf> > status;
                BookmarkColumns() { add(id); add(params); add(paramsEdited); add(name); add(status);}
        };
        BookmarkColumns bookmarkColumns;

    protected:    
        Gtk::ScrolledWindow*    hscrollw;
        Gtk::TreeView*          hTreeView;
        Glib::RefPtr<Gtk::ListStore> historyModel;

        Gtk::ScrolledWindow*    bscrollw;
        Gtk::TreeView*          bTreeView;
        Glib::RefPtr<Gtk::ListStore> bookmarkModel;
        Gtk::CellRendererText m_cellrenderer_validated;
        Gtk::TreeView::Column m_treeviewcolumn_validated;
        bool m_validate_retry;
        Glib::ustring m_invalid_text_for_retry;


        Gtk::Button*            addBookmark;
        Gtk::Button*            delBookmark;

        sigc::connection        selchangehist;
        sigc::connection        selchangebm;
        
        HistoryBeforeLineListener * blistener;
        SnapshotListener *slistener;
        ProfileChangeListener* tpc;
        ParamsEdited defParamsEdited;
        int bmnum;
        static bool iconsLoaded;
        static Glib::RefPtr<Gdk::Pixbuf> recentlySavedIcon;
        static Glib::RefPtr<Gdk::Pixbuf> enqueuedIcon;
        static Glib::RefPtr<Gdk::Pixbuf> voidIcon;
    public:
        
        History (bool bookmarkSupport = true);

        void setProfileChangeListener     (ProfileChangeListener* tpc_) { tpc = tpc_; }
        void setHistoryBeforeLineListener (HistoryBeforeLineListener* bll) { blistener = bll; }
        void setSnapshotListener (SnapshotListener* sl) { slistener = sl; }
        
        // pparamschangelistener interface
        void procParamsChanged (rtengine::procparams::ProcParams* params, rtengine::ProcEvent ev, Glib::ustring descr, ParamsEdited* paramsEdited=NULL);
        void clearParamChanges ();

        void historySelectionChanged ();
        void bookmarkSelectionChanged ();
        void initHistory ();

        bool getBeforeLineParams (rtengine::procparams::ProcParams& params);
        
        void addBookmarkWithText (Glib::ustring text);
        void addBookmarkPressed ();
        void delBookmarkPressed ();
        void addSnapshot(const rtengine::SnapshotInfo &snapInfo );
        void updateSnapshot(const rtengine::SnapshotInfo &snapInfo );
        int  getSelectedSnapshot();
        bool findName( Glib::ustring text );

        void resized (Gtk::Allocation& req);

        void undo ();
        void redo ();
        void cellrenderer_validated_on_editing_started(Gtk::CellEditable* cell_editable, const Glib::ustring& path);
        void cellrenderer_validated_on_edited(const Glib::ustring& path_string, const Glib::ustring& new_text);
        void column_validated_on_cell_data( Gtk::CellRenderer* renderer , const Gtk::TreeModel::iterator& iter);
        bool blistenerLock;
};

#endif
