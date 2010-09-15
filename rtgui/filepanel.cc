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
#include <filepanel.h>
#include <rtwindow.h>

int fbinit (void* data) {

    gdk_threads_enter ();
    ((FilePanel*)data)->init ();
    gdk_threads_leave ();

    return 0;
}

FilePanel::FilePanel () : parent(NULL) {

    dirpaned = new Gtk::HPaned ();
    dirpaned->set_position (options.dirBrowserWidth);

    dirBrowser = new DirBrowser ();
    placesBrowser = new PlacesBrowser ();
    recentBrowser = new RecentBrowser ();

    placespaned = new Gtk::VPaned ();
    placespaned->set_position (options.dirBrowserHeight);

    Gtk::VBox* obox = Gtk::manage (new Gtk::VBox ());
    obox->pack_start (*recentBrowser, Gtk::PACK_SHRINK, 4);
    obox->pack_start (*dirBrowser);

    placespaned->pack1 (*placesBrowser, false, true);
    placespaned->pack2 (*obox, true, true);

    dirpaned->pack1 (*placespaned, true, true);

    tpc = new BatchToolPanelCoordinator (this);
    fileCatalog = new FileCatalog (tpc->coarse, tpc->getToolBar());
    dirpaned->pack2 (*fileCatalog, true, true);

    placesBrowser->setDirBrowserRemoteInterface (dirBrowser);
    recentBrowser->setDirBrowserRemoteInterface (dirBrowser);
    dirBrowser->addDirSelectionListener (fileCatalog);
    dirBrowser->addDirSelectionListener (recentBrowser);
    dirBrowser->addDirSelectionListener (placesBrowser);
    fileCatalog->setFileSelectionListener (this);
    
    rightBox = new Gtk::HBox ();
    rightNotebook = new Gtk::Notebook ();
    Gtk::VBox* taggingBox = new Gtk::VBox ();
    
    history = new History (false);

    tpc->addPParamsChangeListener (history);
    history->setProfileChangeListener (tpc);

    Gtk::ScrolledWindow* sFilterPanel = new Gtk::ScrolledWindow();
	filterPanel = new FilterPanel ();
	sFilterPanel->add (*filterPanel);
	sFilterPanel->set_policy(Gtk::POLICY_NEVER, Gtk::POLICY_AUTOMATIC);

	fileCatalog->setFilterPanel (filterPanel);
    fileCatalog->setImageAreaToolListener (tpc);
    
    //------------------

    rightNotebook->set_tab_pos (Gtk::POS_LEFT);
    
    Gtk::Label* devLab = new Gtk::Label (M("MAIN_TAB_DEVELOP"));
    devLab->set_angle (90);
    Gtk::Label* filtLab = new Gtk::Label (M("MAIN_TAB_FILTER"));
    filtLab->set_angle (90);
    Gtk::Label* tagLab = new Gtk::Label (M("MAIN_TAB_TAGGING"));
    tagLab->set_angle (90);

    Gtk::VPaned* tpcPaned = new Gtk::VPaned ();
    tpcPaned->pack1 (*tpc->toolPanelNotebook, true, true);
    tpcPaned->pack2 (*history, true, true);
    
    rightNotebook->append_page (*tpcPaned, *devLab);
    rightNotebook->append_page (*sFilterPanel, *filtLab);
    rightNotebook->append_page (*taggingBox, *tagLab);

    rightBox->pack_start (*rightNotebook);

    pack1(*dirpaned, true, true);
    pack2(*rightBox, true, true);

    set_position(options.browserToolPanelWidth);

    fileCatalog->setFileSelectionChangeListener (tpc);

    fileCatalog->setFileSelectionListener (this);
    g_idle_add (fbinit, this);

    show_all ();
}

void FilePanel::init () {
  
    dirBrowser->fillDirTree ();
    placesBrowser->refreshPlacesList ();

    if (argv1!="")
        dirBrowser->open (argv1);
    else {
        if (options.startupDir==STARTUPDIR_HOME) 
            dirBrowser->open (Glib::get_home_dir());
        else if (options.startupDir==STARTUPDIR_CURRENT)
            dirBrowser->open (argv0);
        else if (options.startupDir==STARTUPDIR_CUSTOM || options.startupDir==STARTUPDIR_LAST) 
            dirBrowser->open (options.startupPath);
    }
} 

bool FilePanel::fileSelected (Thumbnail* thm) {

    if (!parent)
        return false;

    // try to open the file
    fileCatalog->setEnabled (false);
    ProgressConnector<rtengine::InitialImage*> *ld = new ProgressConnector<rtengine::InitialImage*>();
    ld->startFunc (sigc::bind(sigc::ptr_fun(&rtengine::InitialImage::load), thm->getFileName (), thm->getType()==FT_Raw, &error, parent->getProgressListener()),
   		           sigc::bind(sigc::mem_fun(*this,&FilePanel::imageLoaded), thm, ld) );
    return true;
}
bool FilePanel::imageLoaded( Thumbnail* thm, ProgressConnector<rtengine::InitialImage*> *pc ){

	if (pc->returnValue() && thm) {
		EditorPanel* epanel = Gtk::manage (new EditorPanel ());
		parent->addEditorPanel (epanel,Glib::path_get_basename (thm->getFileName()));
		epanel->open(thm, pc->returnValue() );
	}else {
		Glib::ustring msg_ = Glib::ustring("<b>") + M("MAIN_MSG_CANNOTLOAD") + " \"" + thm->getFileName() + "\" .\n</b>";
		Gtk::MessageDialog msgd (msg_, true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
		msgd.run ();
	}
	delete pc;

	parent->setProgress(0.);
	parent->setProgressStr("");
	fileCatalog->setEnabled (true);
	return false; // MUST return false from idle function
}

void FilePanel::saveOptions () { 

    options.dirBrowserWidth = dirpaned->get_position ();
    options.dirBrowserHeight = placespaned->get_position ();
    options.browserToolPanelWidth = get_position();
    if (options.startupDir==STARTUPDIR_LAST && fileCatalog->lastSelectedDir ()!="")
        options.startupPath = fileCatalog->lastSelectedDir ();
    fileCatalog->closeDir (); 
}

void FilePanel::open (const Glib::ustring& d) {

    if (Glib::file_test (d, Glib::FILE_TEST_IS_DIR))
        dirBrowser->open (d.c_str());
    else if (Glib::file_test (d, Glib::FILE_TEST_EXISTS))
        dirBrowser->open (Glib::path_get_dirname(d), Glib::path_get_basename(d));
}

bool FilePanel::addBatchQueueJob (BatchQueueEntry* bqe) {

    if (parent)
        parent->addBatchQueueJob (bqe);
	return true;
}

void FilePanel::optionsChanged () {

    tpc->optionsChanged ();
    fileCatalog->refreshAll ();
}

bool FilePanel::handleShortcutKey (GdkEventKey* event) {

    bool ctrl = event->state & GDK_CONTROL_MASK;
    bool shift = event->state & GDK_SHIFT_MASK;
    
    if (!ctrl) {
        switch(event->keyval) {
        }
    }
    else {
        switch (event->keyval) {
        }
    }
    
    if(tpc->getToolBar()->handleShortcutKey(event))
        return true;
    
    if(fileCatalog->handleShortcutKey(event))
        return true;

    return false;
}

