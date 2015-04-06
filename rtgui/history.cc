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
#include "history.h"
#include "multilangmgr.h"
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

Glib::ustring eventDescrArray[NUMOFEVENTS];
extern Glib::ustring argv0;

History::History (bool bookmarkSupport) : blistener(NULL), tpc (NULL), bmnum (1) {

	blistenerLock = false; // sets default that the Before preview will not be locked

    // fill history event message array
    for (int i=0; i<NUMOFEVENTS; i++) 
        eventDescrArray[i] = M(Glib::ustring::compose("HISTORY_MSG_%1", i+1));

    // History List
    // ~~~~~~~~~~~~
    hscrollw = Gtk::manage (new Gtk::ScrolledWindow ());
    hscrollw->set_policy (Gtk::POLICY_NEVER, Gtk::POLICY_AUTOMATIC);

    Gtk::Frame* histFrame = Gtk::manage (new Gtk::Frame (M("HISTORY_LABEL")));
    histFrame->add (*hscrollw);
   
    hTreeView = Gtk::manage (new Gtk::TreeView ());
    hscrollw->add (*hTreeView);
    
    historyModel = Gtk::ListStore::create (historyColumns);
    hTreeView->set_model (historyModel);
//    hTreeView->set_headers_visible (false);

    Gtk::CellRendererText *changecrt = Gtk::manage (new Gtk::CellRendererText());
    Gtk::CellRendererText *valuecrt  = Gtk::manage (new Gtk::CellRendererText());
    Gtk::TreeView::Column *hviewcol = Gtk::manage (new Gtk::TreeView::Column (""));
    hviewcol->pack_start (*changecrt, true);
    hviewcol->add_attribute (changecrt->property_markup (), historyColumns.text);
    hviewcol->set_resizable (true);

    Gtk::TreeView::Column *hviewcol2 = Gtk::manage (new Gtk::TreeView::Column (""));
    hviewcol2->pack_start (*valuecrt, true);
    hviewcol2->add_attribute (valuecrt->property_markup (), historyColumns.value);
    valuecrt->set_property ("xalign", 1.0);

    hTreeView->append_column (*hviewcol); 
    hTreeView->append_column (*hviewcol2); 

    hviewcol2->set_sizing (Gtk::TREE_VIEW_COLUMN_FIXED);

    selchangehist = hTreeView->get_selection()->signal_changed().connect(sigc::mem_fun(*this, &History::historySelectionChanged));

    // Bookmark List
    // ~~~~~~~~~~~~~

    Gtk::HSeparator* hsepb = Gtk::manage (new Gtk::HSeparator ());
    pack_end (*hsepb, Gtk::PACK_SHRINK, 0);

    Gtk::HBox* ahbox = Gtk::manage (new Gtk::HBox ());
    addBookmark = Gtk::manage (new Gtk::Button (M("HISTORY_NEWSNAPSHOT")));
    addBookmark->set_tooltip_markup (M("HISTORY_NEWSNAPSHOT_TOOLTIP"));
    Gtk::Image* addimg = Gtk::manage (new RTImage ("gtk-add.png"));
    addBookmark->set_image (*addimg);
    ahbox->pack_start (*addBookmark);

    delBookmark = Gtk::manage (new Gtk::Button (M("HISTORY_DELSNAPSHOT")));
    Gtk::Image* delimg = Gtk::manage (new RTImage ("list-remove.png"));
    delBookmark->set_image (*delimg);
    ahbox->pack_start (*delBookmark);

    bscrollw = Gtk::manage (new Gtk::ScrolledWindow ());
//    bscrollw->set_policy (Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
    bscrollw->set_policy (Gtk::POLICY_NEVER, Gtk::POLICY_AUTOMATIC);
    bscrollw->set_size_request (-1, 75);

    Gtk::Frame* bmFrame = Gtk::manage (new Gtk::Frame (M("HISTORY_SNAPSHOTS")));
    Gtk::VBox* bmBox = Gtk::manage (new Gtk::VBox ());
    bmFrame->add (*bmBox);
    bmBox->pack_start (*bscrollw, Gtk::PACK_EXPAND_WIDGET, 4);
    bmBox->pack_end (*ahbox, Gtk::PACK_SHRINK, 4);
    
    if (bookmarkSupport){
        historyVPaned = Gtk::manage ( new Gtk::VPaned () );
        historyVPaned->pack1 (*histFrame, true, true);
        historyVPaned->pack2 (*bmFrame, false, true);
        pack_start(*historyVPaned);
    }
    else{
        pack_start (*histFrame);
    }

    
    bTreeView = Gtk::manage (new Gtk::TreeView ());
    bscrollw->add (*bTreeView);
    
    bookmarkModel = Gtk::ListStore::create (bookmarkColumns);
    bTreeView->set_model (bookmarkModel);
    bTreeView->set_headers_visible (false);
    bTreeView->append_column_editable (M("HISTORY_SNAPSHOTS"), bookmarkColumns.text);

    selchangebm = bTreeView->get_selection()->signal_changed().connect(sigc::mem_fun(*this, &History::bookmarkSelectionChanged));

    addBookmark->signal_clicked().connect( sigc::mem_fun(*this, &History::addBookmarkPressed) );
    delBookmark->signal_clicked().connect( sigc::mem_fun(*this, &History::delBookmarkPressed) );

//    hTreeView->set_grid_lines (Gtk::TREE_VIEW_GRID_LINES_HORIZONTAL); 
    hTreeView->set_grid_lines (Gtk::TREE_VIEW_GRID_LINES_BOTH); 
    hTreeView->signal_size_allocate().connect( sigc::mem_fun(*this, &History::resized) );

    hTreeView->set_enable_search(false);
    bTreeView->set_enable_search(false);

    show_all_children ();
}

void History::initHistory () {

    historyModel->clear ();
    bookmarkModel->clear ();
}

void History::clearParamChanges () {

    initHistory ();
}

void History::historySelectionChanged () {

    Glib::RefPtr<Gtk::TreeSelection> selection = hTreeView->get_selection();
    Gtk::TreeModel::iterator iter = selection->get_selected();
    if (iter) {
        Gtk::TreeModel::Row row = *iter;
        if (row) 
            bTreeView->get_selection()->unselect_all ();
        if (row && tpc) {
            ProcParams pparams = row[historyColumns.params];
            ParamsEdited pe;
            PartialProfile pp(&pparams, &pe);
            ParamsEdited paramsEdited = row[historyColumns.paramsEdited];
            tpc->profileChange (&pp, EvHistoryBrowsed, row[historyColumns.text], &paramsEdited);
        }
        if (blistener && blistenerLock==false) {
            Gtk::TreeModel::Path path = historyModel->get_path (iter);
            path.prev ();
            iter = historyModel->get_iter (path);
            if (blistener && iter)
                blistener->historyBeforeLineChanged (iter->get_value (historyColumns.params));
        }
    }
}

void History::bookmarkSelectionChanged () {

    Glib::RefPtr<Gtk::TreeSelection> selection = bTreeView->get_selection();
    Gtk::TreeModel::iterator iter = selection->get_selected();
    if (iter) {
        Gtk::TreeModel::Row row = *iter;
        if (row) 
            hTreeView->get_selection()->unselect_all ();
        if (row && tpc) {
            ProcParams pparams = row[bookmarkColumns.params];
            ParamsEdited pe;
            PartialProfile pp(&pparams, &pe);
            ParamsEdited paramsEdited = row[bookmarkColumns.paramsEdited];
            tpc->profileChange (&pp, EvBookmarkSelected, row[bookmarkColumns.text], &paramsEdited);
        }
    }
}

void History::procParamsChanged (ProcParams* params, ProcEvent ev, Glib::ustring descr, ParamsEdited* paramsEdited) {

    // to prevent recursion, we filter out the events triggered by the history
    if (ev==EvHistoryBrowsed)
        return;

    selchangehist.block (true);

    if (ev==EvPhotoLoaded)
        initHistory ();
    // construct formatted list content
    Glib::ustring text = Glib::ustring::compose ("%1", eventDescrArray[ev]);

    Glib::RefPtr<Gtk::TreeSelection> selection = hTreeView->get_selection();
    Gtk::TreeModel::iterator iter = selection->get_selected();
    // remove all rows after the selection
    if (iter) {
        iter++;
        while (iter)
            iter = historyModel->erase (iter);
    }
    // lookup the last remaining item in the list
    int size = historyModel->children().size ();
    Gtk::TreeModel::Row row;
    if (size>0)
        row = historyModel->children()[size-1];
    // if there is no last item or its chev!=ev, create a new one
    if (size==0 || !row || row[historyColumns.chev]!=ev || ev==EvProfileChanged) {
        Gtk::TreeModel::Row newrow = *(historyModel->append());
        newrow[historyColumns.realText] = eventDescrArray[ev];
        newrow[historyColumns.text] = text;
        newrow[historyColumns.value] = descr;
        newrow[historyColumns.chev] = ev;
        newrow[historyColumns.params] = *params;
        newrow[historyColumns.paramsEdited] = paramsEdited ? *paramsEdited : defParamsEdited;
        if (ev!=EvBookmarkSelected)
            selection->select (newrow);
        if (blistener && row && blistenerLock==false)
            blistener->historyBeforeLineChanged (row[historyColumns.params]);
        else if (blistener && size==0 && blistenerLock==false)
            blistener->historyBeforeLineChanged (newrow[historyColumns.params]);
    }
    // else just update it
    else {
        row[historyColumns.realText] = eventDescrArray[ev];
        row[historyColumns.text] = text;
        row[historyColumns.value] = descr;
        row[historyColumns.chev] = ev;
        row[historyColumns.params] = *params;
        row[historyColumns.paramsEdited] = paramsEdited ? *paramsEdited : defParamsEdited;
        if (ev!=EvBookmarkSelected)
            selection->select (row);
    }
    if (ev!=EvBookmarkSelected)
        bTreeView->get_selection()->unselect_all ();


    if (!selection->get_selected_rows().empty()) {
        Gtk::TreeView::Selection::ListHandle_Path selp = selection->get_selected_rows();
        hTreeView->scroll_to_row (*selp.begin());
    }
    selchangehist.block (false);
}

void History::addBookmarkWithText (Glib::ustring text) {

    // lookup the selected item in the history
    Glib::RefPtr<Gtk::TreeSelection> selection = hTreeView->get_selection();
    Gtk::TreeModel::iterator iter = selection->get_selected();
    Gtk::TreeModel::Row row = *iter;

    if (!row) {
        return;
    }

    // append new row to bookmarks
    Gtk::TreeModel::Row newrow = *(bookmarkModel->append());
    newrow[bookmarkColumns.text] = text;
    ProcParams params = row[historyColumns.params];
    newrow[bookmarkColumns.params] = params;
    ParamsEdited paramsEdited = row[historyColumns.paramsEdited];
    newrow[bookmarkColumns.paramsEdited] = paramsEdited;
}

void History::addBookmarkPressed () {
    
    if (hTreeView->get_selection()->get_selected()) 
        addBookmarkWithText (Glib::ustring::compose ("%1 %2", M("HISTORY_SNAPSHOT"), bmnum++));
}

void History::delBookmarkPressed () {

    // lookup the selected item in the bookmark
    Glib::RefPtr<Gtk::TreeSelection> selection = bTreeView->get_selection();
    Gtk::TreeModel::iterator iter = selection->get_selected();

    if (!iter) {
        return;
    }
    // remove selected bookmark
    bookmarkModel->erase (iter);
    // select last item in history
    int size = historyModel->children().size ();
    Gtk::TreeModel::Row row = historyModel->children()[size-1];
    hTreeView->get_selection()->select (row);
}

void History::undo () {

    Glib::RefPtr<Gtk::TreeSelection> selection = hTreeView->get_selection();
    Gtk::TreeModel::iterator iter = selection->get_selected();

    if (iter && iter!=historyModel->children().begin()) 
        selection->select (--iter);
    else if (!iter) {
        int size = historyModel->children().size ();
        if (size>1) 
            selection->select (historyModel->children()[size-2]); 
    }
}

void History::redo () {

    Glib::RefPtr<Gtk::TreeSelection> selection = hTreeView->get_selection();
    Gtk::TreeModel::iterator iter = selection->get_selected();

    if (iter) {
        iter++;
        if (iter!=historyModel->children().end())
            selection->select (iter);
    }
    else {
        int size = historyModel->children().size ();
        if (size>1) 
            selection->select (historyModel->children()[size-2]); 
    }
}

void History::resized (Gtk::Allocation& req) {
}

bool History::getBeforeLineParams (rtengine::procparams::ProcParams& params) {

    int size = historyModel->children().size ();
    if (size==0 || !blistener)
        return false;

    Gtk::TreeModel::Row row;
    row = historyModel->children()[size==1 ? 0 : size-2];
    params = row[historyColumns.params];
    return true;
}

