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
#include <history.h>
#include <multilangmgr.h>
#include <rtimage.h>

using namespace rtengine;
using namespace rtengine::procparams;

Glib::ustring eventDescrArray[NUMOFEVENTS];
extern Glib::ustring argv0;

History::History (bool bookmarkSupport) : blistener(NULL), tpc (NULL), slistener(NULL), bmnum (1) {

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

    pack_start (*histFrame);
    
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
    bmBox->pack_start (*bscrollw, Gtk::PACK_SHRINK, 4);
    bmBox->pack_end (*ahbox, Gtk::PACK_SHRINK, 4);
    
    if (bookmarkSupport)
        pack_end (*bmFrame, Gtk::PACK_SHRINK, 4);
    
    bTreeView = Gtk::manage (new Gtk::TreeView ());
    bscrollw->add (*bTreeView);
    
    bookmarkModel = Gtk::ListStore::create (bookmarkColumns);
    bTreeView->set_model (bookmarkModel);
    bTreeView->set_headers_visible (false);
    //bTreeView->append_column_editable (M("HISTORY_SNAPSHOTS"), bookmarkColumns.text);

    m_treeviewcolumn_validated.set_title("validated");
    m_treeviewcolumn_validated.pack_start(m_cellrenderer_validated);
    bTreeView->append_column(m_treeviewcolumn_validated);
    m_cellrenderer_validated.property_editable() = true;
    m_cellrenderer_validated.signal_editing_started().connect(sigc::mem_fun(*this,  &History::cellrenderer_validated_on_editing_started) );
    m_cellrenderer_validated.signal_edited().connect( sigc::mem_fun(*this,          &History::cellrenderer_validated_on_edited) );
    m_treeviewcolumn_validated.set_cell_data_func(m_cellrenderer_validated,  sigc::mem_fun(*this, &History::column_validated_on_cell_data) );
    m_validate_retry=false;

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


void History::cellrenderer_validated_on_editing_started(Gtk::CellEditable* cell_editable, const Glib::ustring& path)
{
	  //Start editing with previously-entered (but invalid) text,
	  //if we are allowing the user to correct some invalid data.
	  if(m_validate_retry)
	  {
	    //This is the CellEditable inside the CellRenderer.
	    Gtk::CellEditable* celleditable_validated = cell_editable;

	    //It's usually an Entry, at least for a CellRendererText:
	    Gtk::Entry* pEntry = dynamic_cast<Gtk::Entry*>(celleditable_validated);
	    if(pEntry)
	    {
	      pEntry->set_text(m_invalid_text_for_retry);
	      m_validate_retry = false;
	      m_invalid_text_for_retry.clear();
	    }
	  }
}

void History::cellrenderer_validated_on_edited(const Glib::ustring& path_string, const Glib::ustring& new_text)
{
	 Gtk::TreePath path(path_string);
	 Gtk::TreeModel::iterator iter = bookmarkModel->get_iter(path);
	 //Prevent entry of the same name
	 if( new_text.compare( SnapshotInfo::kCurrentSnapshotName)==0 || findName( new_text ) ){

	    //Start editing again, with the bad text, so that the user can correct it.
	    //A real application should probably allow the user to revert to the
	    //previous text.

	    //Set the text to be used in the start_editing signal handler:
		if(iter)
	       m_invalid_text_for_retry = (*iter)[bookmarkColumns.name];
		else
			m_invalid_text_for_retry = new_text;
	    m_validate_retry = true;

	    //Start editing again:
	    bTreeView->set_cursor(path, m_treeviewcolumn_validated, m_cellrenderer_validated, true /* start_editing */);
	  }
	  else
	  {

	    if(iter)
	    {
	      Gtk::TreeModel::Row row = *iter;
	      Glib::ustring old_value(row[bookmarkColumns.name]);
	      //Put the new value in the model:
	      row[bookmarkColumns.name] = new_text;
	      int id = row[bookmarkColumns.id];
	      if( slistener )
	    	  slistener->renameSnapshot(id,new_text );
	    }
	  }
}

void History::column_validated_on_cell_data( Gtk::CellRenderer*  renderer , const Gtk::TreeModel::iterator& iter)
{
	  //Get the value from the model and show it appropriately in the view:
	  if(iter)
	  {
	    Gtk::TreeModel::Row row = *iter;
	    Glib::ustring value = row[bookmarkColumns.name];
	    m_cellrenderer_validated.property_text() = value;
	  }

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
            ProcParams params = row[historyColumns.params];
            ParamsEdited paramsEdited = row[historyColumns.paramsEdited];
            tpc->profileChange (&params, EvHistoryBrowsed, row[historyColumns.text], &paramsEdited);
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
        	int id = row[bookmarkColumns.id];
            ProcParams params = row[bookmarkColumns.params];
            ParamsEdited paramsEdited = row[bookmarkColumns.paramsEdited];
            tpc->profileChange (&params, EvBookmarkSelected, row[bookmarkColumns.name], &paramsEdited);
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
    Glib::ustring text = Glib::ustring::compose ("<b>%1</b>", eventDescrArray[ev]);

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
    if (size==0 || !row || row[historyColumns.chev]!=ev) {
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
    ProcParams params = row[historyColumns.params];
    ParamsEdited paramsEdited = row[historyColumns.paramsEdited];
    int newId = rand();
    if( slistener ){
    	newId = slistener->newSnapshot( text, params);
    	if( newId<0 )
    		return;
    }

    // append new row to bookmarks
    Gtk::TreeModel::Row newrow = *(bookmarkModel->append());
    //newrow[bookmarkColumns.text] = text;
    newrow[bookmarkColumns.id] = newId;
    newrow[bookmarkColumns.name] = text;
    newrow[bookmarkColumns.params] = params;
    newrow[bookmarkColumns.paramsEdited] = paramsEdited;

}

void History::addSnapshot( const rtengine::SnapshotInfo &snapInfo )
{
    // append new row to bookmarks
    Gtk::TreeModel::Row newrow = *(bookmarkModel->append());
    newrow[bookmarkColumns.id] = snapInfo.id;
    newrow[bookmarkColumns.name] = snapInfo.name;
    newrow[bookmarkColumns.params] = snapInfo.params;
    ParamsEdited paramsEdited;
    paramsEdited.set(true);
    newrow[bookmarkColumns.paramsEdited] = paramsEdited;

}

int  History::getSelectedSnapshot()
{
    // lookup the selected item in the bookmark
    Glib::RefPtr<Gtk::TreeSelection> selection = bTreeView->get_selection();
    Gtk::TreeModel::iterator iter = selection->get_selected();

    if (!iter) {
        return -1;
    }
    return (*iter)[bookmarkColumns.id];
}

bool History::findName( Glib::ustring text )
{
	Gtk::ListStore::Children l= bookmarkModel->children();

	for( Gtk::ListStore::Children::iterator iter = l.begin();iter!=l.end();iter++ ){
		Glib::ustring s = (*iter)[bookmarkColumns.name];
		if( text.compare( s )==0 )
			return true;
	}
	return false;
}

void History::addBookmarkPressed () {
    
    if (hTreeView->get_selection()->get_selected()) {
    	Glib::ustring newName;
    	do{
    		newName = Glib::ustring::compose ("%1 %2", M("HISTORY_SNAPSHOT"), bmnum++);
    	}while( findName(newName) );
        addBookmarkWithText (newName);
    }
}

void History::delBookmarkPressed () {

    // lookup the selected item in the bookmark
    Glib::RefPtr<Gtk::TreeSelection> selection = bTreeView->get_selection();
    Gtk::TreeModel::iterator iter = selection->get_selected();

    if (!iter) {
        return;
    }
    Glib::ustring text = (*iter)[bookmarkColumns.name];
    int id = (*iter)[bookmarkColumns.id];
    if( slistener ){
    	if( !slistener->deleteSnapshot( id ) )
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

