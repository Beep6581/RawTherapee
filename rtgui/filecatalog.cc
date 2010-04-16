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
#include <filecatalog.h>
#include <options.h>
#include <cachemanager.h>
#include <multilangmgr.h>
#include <guiutils.h>
#include <glib/gstdio.h>
#include <iostream>
#include <renamedlg.h>
#include <thumbimageupdater.h>
#include <safegtk.h>

#define CHECKTIME 2000

extern Glib::ustring argv0;

#ifdef _WIN32
int _directoryUpdater (void* cat) {

    ((FileCatalog*)cat)->checkCounter++;
    if (((FileCatalog*)cat)->checkCounter==2) {
        gdk_threads_enter ();
        ((FileCatalog*)cat)->reparseDirectory ();
        gdk_threads_leave ();
    }
    return 1;
}
#endif

FileCatalog::FileCatalog (CoarsePanel* cp, ToolBar* tb) : listener(NULL), fslistener(NULL), hasValidCurrentEFS(false), coarsePanel(cp), filterPanel(NULL), toolBar(tb) {

    previewLoader.setPreviewLoaderListener (this);
    
    //  construct and initialize thumbnail browsers
    fileBrowser = new FileBrowser();
    fileBrowser->setFileBrowserListener (this);
    fileBrowser->setArrangement (ThumbBrowserBase::TB_Vertical);
    fileBrowser->show ();
       
    // construct trash panel with the extra "empty trash" button
    trashButtonBox = new Gtk::VBox;
    Gtk::Button* emptyT = new Gtk::Button (M("FILEBROWSER_EMPTYTRASH"));
    emptyT->set_tooltip_text (M("FILEBROWSER_EMPTYTRASHHINT"));
    emptyT->set_image (*(new Gtk::Image (Gtk::StockID("gtk-delete"), Gtk::ICON_SIZE_BUTTON)));
    emptyT->signal_pressed().connect (sigc::mem_fun(*this, &FileCatalog::emptyTrash));    
    trashButtonBox->pack_start (*emptyT, Gtk::PACK_SHRINK, 4);
    emptyT->show ();
    trashButtonBox->show ();
    
    // setup button bar
    buttonBar = new Gtk::HBox ();
    pack_start (*buttonBar, Gtk::PACK_SHRINK);
    
    buttonBar->pack_start (*(new Gtk::VSeparator), Gtk::PACK_SHRINK);
    bDir = new Gtk::ToggleButton ();
    bDir->set_active (true);
    bDir->set_image (*(new Gtk::Image (argv0+"/images/folder.png")));
    bDir->set_relief (Gtk::RELIEF_NONE);
    bDir->set_tooltip_text (M("FILEBROWSER_SHOWDIRHINT"));
    bCateg[0] = bDir->signal_toggled().connect (sigc::bind(sigc::mem_fun(*this, &FileCatalog::categoryButtonToggled), bDir));    
    buttonBar->pack_start (*bDir, Gtk::PACK_SHRINK);
    buttonBar->pack_start (*(new Gtk::VSeparator), Gtk::PACK_SHRINK);

    bUnRanked = new Gtk::ToggleButton ();
    bUnRanked->set_active (false);
    bUnRanked->set_image (*(new Gtk::Image (argv0+"/images/unrated.png")));
    bUnRanked->set_relief (Gtk::RELIEF_NONE);
    bUnRanked->set_tooltip_text (M("FILEBROWSER_SHOWUNRANKHINT"));
    bCateg[1] = bUnRanked->signal_toggled().connect (sigc::bind(sigc::mem_fun(*this, &FileCatalog::categoryButtonToggled), bUnRanked));    
    buttonBar->pack_start (*bUnRanked, Gtk::PACK_SHRINK);
    buttonBar->pack_start (*(new Gtk::VSeparator), Gtk::PACK_SHRINK);

    for (int i=0; i<5; i++) {
        iranked[i] = new Gtk::Image (argv0+"/images/rated.png");
        igranked[i] = new Gtk::Image (argv0+"/images/grayrated.png");
        iranked[i]->show ();
        igranked[i]->show ();
        bRank[i] = new Gtk::ToggleButton ();
        bRank[i]->set_image (*igranked[i]);
        bRank[i]->set_relief (Gtk::RELIEF_NONE);
        buttonBar->pack_start (*bRank[i], Gtk::PACK_SHRINK);
        bCateg[i+2] = bRank[i]->signal_toggled().connect (sigc::bind(sigc::mem_fun(*this, &FileCatalog::categoryButtonToggled), bRank[i]));    
    }  
    bRank[0]->set_tooltip_text (M("FILEBROWSER_SHOWRANK1HINT"));
    bRank[1]->set_tooltip_text (M("FILEBROWSER_SHOWRANK2HINT"));
    bRank[2]->set_tooltip_text (M("FILEBROWSER_SHOWRANK3HINT"));
    bRank[3]->set_tooltip_text (M("FILEBROWSER_SHOWRANK4HINT"));
    bRank[4]->set_tooltip_text (M("FILEBROWSER_SHOWRANK5HINT"));
    buttonBar->pack_start (*(new Gtk::VSeparator), Gtk::PACK_SHRINK);

    bTrash = new Gtk::ToggleButton ();
    bTrash->set_image (*(new Gtk::Image (Gtk::StockID("gtk-delete"), Gtk::ICON_SIZE_SMALL_TOOLBAR)));
    bTrash->set_relief (Gtk::RELIEF_NONE);
    bTrash->set_tooltip_text (M("FILEBROWSER_SHOWTRASHHINT"));
    bCateg[7] = bTrash->signal_toggled().connect (sigc::bind(sigc::mem_fun(*this, &FileCatalog::categoryButtonToggled), bTrash));    
    buttonBar->pack_start (*bTrash, Gtk::PACK_SHRINK);
    buttonBar->pack_start (*(new Gtk::VSeparator), Gtk::PACK_SHRINK);

    categoryButtons[0] = bDir;
    categoryButtons[1] = bUnRanked;
    for (int i=0; i<5; i++)
        categoryButtons[i+2] = bRank[i];
    categoryButtons[7] = bTrash;

    // thumbnail zoom
    Gtk::HBox* zoomBox = new Gtk::HBox ();
    zoomInButton  = new Gtk::Button ();
    zoomInButton->set_image (*(new Gtk::Image (Gtk::StockID("gtk-zoom-in"), Gtk::ICON_SIZE_SMALL_TOOLBAR)));
    zoomInButton->signal_pressed().connect (sigc::mem_fun(*this, &FileCatalog::zoomIn));    
    zoomInButton->set_relief (Gtk::RELIEF_NONE);
    zoomInButton->set_tooltip_text (M("FILEBROWSER_ZOOMINHINT"));
    zoomBox->pack_end (*zoomInButton, Gtk::PACK_SHRINK);
    zoomOutButton  = new Gtk::Button ();
    zoomOutButton->set_image (*(new Gtk::Image (Gtk::StockID("gtk-zoom-out"), Gtk::ICON_SIZE_SMALL_TOOLBAR)));
    zoomOutButton->signal_pressed().connect (sigc::mem_fun(*this, &FileCatalog::zoomOut));    
    zoomOutButton->set_relief (Gtk::RELIEF_NONE);
    zoomOutButton->set_tooltip_text (M("FILEBROWSER_ZOOMOUTHINT"));
    zoomBox->pack_end (*zoomOutButton, Gtk::PACK_SHRINK);   

    // add default panel 
    hBox = new Gtk::HBox ();
    hBox->show ();
    hBox->pack_end (*fileBrowser);
    fileBrowser->applyFilter (getFilter());
    pack_start (*hBox);

    buttonBar2 = new Gtk::HBox ();
    pack_end (*buttonBar2, Gtk::PACK_SHRINK);
    progressBar = new Gtk::ProgressBar ();
    buttonBar2->pack_start (*progressBar, Gtk::PACK_SHRINK, 4);
    progressBar->set_size_request (-1, 16);

    buttonBar->pack_start (*zoomBox, Gtk::PACK_SHRINK);   

    buttonBar->pack_end (*coarsePanel, Gtk::PACK_SHRINK);   
    buttonBar->pack_end (*Gtk::manage(new Gtk::VSeparator), Gtk::PACK_SHRINK, 4);
    buttonBar->pack_end (*toolBar, Gtk::PACK_SHRINK);   
    buttonBar->pack_end (*Gtk::manage(new Gtk::VSeparator), Gtk::PACK_SHRINK, 4);

    enabled = true;  

    lastScrollPos = 0;
    for (int i=0; i<8; i++) {
        hScrollPos[i] = 0;
        vScrollPos[i] = 0;
    }
    
    selectedDirectory = "";
#ifdef _WIN32
    wdMonitor = NULL;
    checkCounter = 2;
    g_timeout_add (CHECKTIME, _directoryUpdater, this);
#endif
}

void FileCatalog::on_realize() {

    Gtk::VBox::on_realize();
    Pango::FontDescription fontd = get_pango_context()->get_font_description ();  
    fileBrowser->get_pango_context()->set_font_description (fontd);
//    batchQueue->get_pango_context()->set_font_description (fontd);
}

void FileCatalog::closeDir () {

	if (filterPanel)
		filterPanel->set_sensitive (false);

#ifndef _WIN32
    if (dirMonitor)
        dirMonitor->cancel ();
#else
    if (wdMonitor) {
        delete wdMonitor;
        wdMonitor = NULL;
    }
#endif
    // terminate thumbnail preview loading
    previewLoader.terminate ();

    // terminate thumbnail updater
    thumbImageUpdater.terminate ();

    // remove entries
    fileBrowser->close ();  
	fileNameList.clear ();
	
    dirEFS.clear ();
    hasValidCurrentEFS = false;
    selectedDirectory = "";
    redrawAll ();
}

std::vector<Glib::ustring> FileCatalog::getFileList () {

    std::vector<Glib::ustring> names;
    Glib::RefPtr<Gio::File> dir = Gio::File::create_for_path (selectedDirectory);
		safe_build_file_list (dir, names, selectedDirectory);
    return names;
}

void FileCatalog::dirSelected (const Glib::ustring& dirname, const Glib::ustring& openfile) {

    try {
        Glib::RefPtr<Gio::File> dir = Gio::File::create_for_path (dirname);

        if (!dir)
            return;
        closeDir ();
        previewsToLoad = 0;
        previewsLoaded = 0;
        // if openfile exists, we have to open it first (it is a command line argument)
        if (openfile!="") 
            addAndOpenFile (openfile);

		selectedDirectory = dir->get_parse_name();
        fileNameList = getFileList ();

        for (int i=0; i<fileNameList.size(); i++) {
            Glib::RefPtr<Gio::File> f = Gio::File::create_for_path(fileNameList[i]);
            if (f->get_parse_name() != openfile) // if we opened a file at the beginning dont add it again
                checkAndAddFile (f);
        }

        _refreshProgressBar ();
        previewLoader.process ();
		
#ifdef _WIN32
      wdMonitor = new WinDirMonitor (selectedDirectory, this);
#elif defined __APPLE__
      printf("TODO fix dir->monitor_directory () for OSX\n");
#else  
        dirMonitor = dir->monitor_directory ();
        dirMonitor->signal_changed().connect (sigc::bind(sigc::mem_fun(*this, &FileCatalog::on_dir_changed), false));
#endif
    }
    catch (Glib::Exception& ex) {
        std::cout << ex.what();
    }
}

void FileCatalog::_refreshProgressBar () {

    // check if progress bar is visible
/*    Glib::ListHandle<Gtk::Widget*> list = buttonBar2->get_children ();
    Glib::ListHandle<Gtk::Widget*>::iterator i = list.begin ();
    for (; i!=list.end() && *i!=progressBar; i++);
    if (i==list.end()) {
        buttonBar2->pack_start (*progressBar, Gtk::PACK_SHRINK, 4);
        buttonBar2->reorder_child (*progressBar, 2);
    }
*/
    progressBar->show ();    
    if (previewsToLoad>0)
        progressBar->set_fraction ((double)previewsLoaded / previewsToLoad);
    else
        progressBar->set_fraction (1.0);
}

int refreshpb (void* data) {

    gdk_threads_enter ();
    ((FileCatalog*)data)->_refreshProgressBar ();
    gdk_threads_leave ();
    return 0;
}

void FileCatalog::previewReady (FileBrowserEntry* fdn) {

    // put it into the "full directory" browser
    fdn->setImageAreaToolListener (iatlistener);
    fileBrowser->addEntry (fdn);
    
    // update exif filter settings (minimal & maximal values of exif tags, cameras, lenses, etc...)
    const CacheImageData* cfs = fdn->thumbnail->getCacheImageData();
    if (cfs->exifValid) {
        if (cfs->fnumber < dirEFS.fnumberFrom)
            dirEFS.fnumberFrom = cfs->fnumber;
        if (cfs->fnumber > dirEFS.fnumberTo)
            dirEFS.fnumberTo = cfs->fnumber;
        if (cfs->shutter < dirEFS.shutterFrom)
            dirEFS.shutterFrom = cfs->shutter;
        if (cfs->shutter > dirEFS.shutterTo)
            dirEFS.shutterTo = cfs->shutter;
        if (cfs->iso>0 && cfs->iso < dirEFS.isoFrom)
            dirEFS.isoFrom = cfs->iso;
        if (cfs->iso>0 && cfs->iso > dirEFS.isoTo)
            dirEFS.isoTo = cfs->iso;
        if (cfs->focalLen < dirEFS.focalFrom)
            dirEFS.focalFrom = cfs->focalLen;
        if (cfs->focalLen > dirEFS.focalTo)
            dirEFS.focalTo = cfs->focalLen;
    }
    dirEFS.cameras.insert (cfs->camera);
    dirEFS.lenses.insert (cfs->lens);
    previewsLoaded++;
    g_idle_add (refreshpb, this);
}

int prevfinished (void* data) {

    gdk_threads_enter();
    ((FileCatalog*)data)->_previewsFinished ();
    gdk_threads_leave();
    return 0;
}

void FileCatalog::_previewsFinished () {

    redrawAll ();
    previewsToLoad = 0;
    previewsLoaded = 0;
//    removeIfThere (buttonBar2, progressBar);
    progressBar->hide ();
	if (filterPanel) {
		filterPanel->set_sensitive (true);
	    if ( !hasValidCurrentEFS ){
	        currentEFS = dirEFS;
		    filterPanel->setFilter ( dirEFS,true );
	    }else {
		    filterPanel->setFilter ( currentEFS,false );
	    }
	}
}

void FileCatalog::previewsFinished () {

    if (!hasValidCurrentEFS) 
        currentEFS = dirEFS;
    g_idle_add (prevfinished, this);
}

void PreviewLoader::remove (Glib::ustring fname) {
	std::list<DirEntry>::iterator i;
	for (i=jqueue.begin(); i!=jqueue.end(); i++) 
		if (i->fullName==fname)
			break;
	if (i!=jqueue.end())
		jqueue.erase (i);	
}

void PreviewLoader::start () {

    jqueue.sort ();
}

void PreviewLoader::process (DirEntry& current) {

	if (Glib::file_test (current.fullName, Glib::FILE_TEST_EXISTS)) {
	    Thumbnail* tmb = cacheMgr.getEntry (current.fullName);
	    if (tmb && pl) 
	  	    pl->previewReady (new FileBrowserEntry (tmb, current.fullName));
	}
}

void PreviewLoader::end () {

    if (pl)
        pl->previewsFinished ();
}

void FileCatalog::setEnabled (bool e) {

    enabled = e;
}

void FileCatalog::redrawAll () {

    fileBrowser->queue_draw ();
}

void FileCatalog::refreshAll () {

    fileBrowser->refreshThumbImages ();
}

void FileCatalog::_openImage (std::vector<Thumbnail*> tmb) {

    if (enabled && listener!=NULL) {
        previewLoader.stop ();
        thumbImageUpdater.stop ();
        for (int i=0; i<tmb.size(); i++) {
            if (editedFiles.find (tmb[i]->getFileName())==editedFiles.end())
                listener->fileSelected (tmb[i]);
            tmb[i]->decreaseRef ();
        }
        previewLoader.process ();
        thumbImageUpdater.process ();
    }
}

struct FCOIParams {
    FileCatalog* catalog;
    std::vector<Thumbnail*> tmb;
};

int fcopenimg (void* p) {

    gdk_threads_enter ();
    FCOIParams* params = (FCOIParams*)p;
    params->catalog->_openImage (params->tmb);
    delete params;
    gdk_threads_leave ();
    return 0;
}

void FileCatalog::openRequested  (std::vector<Thumbnail*> tmb) {

    FCOIParams* params = new FCOIParams;
    params->catalog = this;
    params->tmb = tmb;
    for (int i=0; i<tmb.size(); i++)
        tmb[i]->increaseRef ();
    g_idle_add (fcopenimg, params);
}

void FileCatalog::deleteRequested  (std::vector<FileBrowserEntry*> tbe) {

    if (tbe.size()==0)
        return;

    Gtk::MessageDialog msd (M("FILEBROWSER_DELETEDLGLABEL"), false, Gtk::MESSAGE_QUESTION, Gtk::BUTTONS_YES_NO, true);
    msd.set_secondary_text(Glib::ustring::compose (M("FILEBROWSER_DELETEDLGMSG"), tbe.size()));

    if (msd.run()==Gtk::RESPONSE_YES) {
        for (int i=0; i<tbe.size(); i++) {
            Glib::ustring fname = tbe[i]->filename;
            // remove from browser
            FileBrowserEntry* t = fileBrowser->delEntry (fname);
//            t->thumbnail->decreaseRef ();
            delete t;
            // remove from cache
            cacheMgr.deleteEntry (fname);
            // delete from file system
            ::g_remove (fname.c_str());
            // delete .pp2 if found
            ::g_remove (Glib::ustring(fname+".pp2").c_str());
            ::g_remove (Glib::ustring(removeExtension(fname)+".pp2").c_str());
            // delete .thm file
            ::g_remove (Glib::ustring(removeExtension(fname)+".thm").c_str());
            ::g_remove (Glib::ustring(removeExtension(fname)+".THM").c_str());
        }
        redrawAll ();    
    }
}

void FileCatalog::developRequested (std::vector<FileBrowserEntry*> tbe) {

    if (listener) {
        thumbImageUpdater.stop ();
        for (int i=0; i<tbe.size(); i++) {
            rtengine::procparams::ProcParams params = tbe[i]->thumbnail->getProcParams();
            rtengine::ProcessingJob* pjob = rtengine::ProcessingJob::create (tbe[i]->filename, tbe[i]->thumbnail->getType()==FT_Raw, params);
            double tmpscale;
            rtengine::IImage8* img = tbe[i]->thumbnail->processThumbImage (params, options.maxThumbnailHeight, tmpscale);
            if (img) {
                int pw = img->getWidth ();
                int ph = img->getHeight ();
                guint8* prev = new guint8 [pw*ph*3];
                memcpy (prev, img->getData (), pw*ph*3);
                listener->addBatchQueueJob (new BatchQueueEntry (pjob, params, tbe[i]->filename, prev, pw, ph, tbe[i]->thumbnail));
            }
            else {
                int pw, ph;
                tbe[i]->thumbnail->getThumbnailSize (pw, ph);
                listener->addBatchQueueJob (new BatchQueueEntry (pjob, params, tbe[i]->filename, NULL, pw, ph, tbe[i]->thumbnail));
            }
        }
        thumbImageUpdater.process ();
    }
}

void FileCatalog::renameRequested  (std::vector<FileBrowserEntry*> tbe) {

    RenameDialog* renameDlg = new RenameDialog ((Gtk::Window*)get_toplevel());

    for (int i=0; i<tbe.size(); i++) {
        renameDlg->initName (Glib::path_get_basename (tbe[i]->filename), tbe[i]->thumbnail->getCacheImageData());

        Glib::ustring ofname = tbe[i]->filename;
        Glib::ustring dirName = Glib::path_get_dirname (tbe[i]->filename);
        Glib::ustring baseName = Glib::path_get_basename (tbe[i]->filename);

        if (renameDlg->run ()== Gtk::RESPONSE_OK) {
            Glib::ustring nBaseName = renameDlg->getNewName ();
            // if path has directory components, exit
            if (Glib::path_get_dirname (nBaseName) != ".")
                continue;
            // if no extension is given, concatenate the extension of the original file
            Glib::ustring ext = getExtension (nBaseName);
            if (ext=="") 
                nBaseName += "." + getExtension (baseName);
            Glib::ustring nfname = Glib::build_filename (dirName, nBaseName);
            if (!::g_rename (ofname.c_str(), nfname.c_str())) {
				cacheMgr.renameEntry (ofname, tbe[i]->thumbnail->getMD5(), nfname);
				reparseDirectory ();
            }
            renameDlg->hide ();
        }
    }
    delete renameDlg;
/*    // ask for new file name
    Gtk::Dialog dialog (M("FILEBROWSER_RENAMEDLGLABEL"), *((Gtk::Window*)get_toplevel()), true, true);
    
    dialog.add_button (Gtk::Stock::OK, Gtk::RESPONSE_OK);
    dialog.add_button (Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);

    Gtk::Label l;
    dialog.get_vbox()->pack_start (l, Gtk::PACK_SHRINK);

    Gtk::Entry nfentry;

    dialog.get_vbox()->pack_start (nfentry, Gtk::PACK_SHRINK);
    dialog.get_vbox()->show_all ();

    nfentry.set_activates_default (true);
    dialog.set_default_response (Gtk::RESPONSE_OK);

    for (int i=0; i<tbe.size(); i++) {

        Glib::ustring ofname = tbe[i]->filename;
        Glib::ustring dirName = Glib::path_get_dirname (tbe[i]->filename);
        Glib::ustring baseName = Glib::path_get_basename (tbe[i]->filename);

        l.set_markup (Glib::ustring("<big><b>") + Glib::ustring::compose (M("FILEBROWSER_RENAMEDLGMSG"), baseName) + Glib::ustring("</b></big>"));
        nfentry.set_text (baseName);
        nfentry.select_region (0, baseName.size());

        if (dialog.run ()== Gtk::RESPONSE_OK) {
            Glib::ustring nBaseName = nfentry.get_text ();
            // if path has directory components, exit
            if (Glib::path_get_dirname (nBaseName) != ".")
                continue;
            // if no extension is given, concatenate the extension of the original file
            if (nBaseName.find ('.')==nBaseName.npos) {
                int lastdot = baseName.find_last_of ('.');
                nBaseName += "." + (lastdot!=Glib::ustring::npos ? baseName.substr (lastdot+1) : "");
            }
            Glib::ustring nfname = Glib::build_filename (dirName, nBaseName);
            if (!::g_rename (ofname.c_str(), nfname.c_str())) {
				cacheMgr.renameEntry (ofname, tbe[i]->thumbnail->getMD5(), nfname);
				// the remaining part (removing old and adding new entry) is done by the directory monitor
				reparseDirectory ();
//                on_dir_changed (Gio::File::create_for_path (nfname), Gio::File::create_for_path (nfname), Gio::FILE_MONITOR_EVENT_CHANGED, true);
            }
        }
    }
    */
}

void FileCatalog::categoryButtonToggled (Gtk::ToggleButton* b) {

    for (int i=0; i<8; i++)
        bCateg[i].block (true);

    fileBrowser->getScrollPosition (hScrollPos[lastScrollPos], vScrollPos[lastScrollPos]);

    // seek the one pressed
    for (int i=0; i<8; i++) {
        categoryButtons[i]->set_active (categoryButtons[i]==b);
        if (categoryButtons[i]==b)
            lastScrollPos = i;
    }
       
    // change the images of the buttons to reflect current ranking
    for (int i=0; i<5; i++) 
        bRank[i]->set_image (*igranked[i]);
    for (int i=0; i<5; i++)
        if (b==bRank[i]) 
            for (int j=0; j<=i; j++)
                bRank[j]->set_image (*iranked[j]);               

    fileBrowser->applyFilter (getFilter ());

    // rearrange panels according to the selected filter
    removeIfThere (hBox, trashButtonBox);
    if (bTrash->get_active ())
        hBox->pack_start (*trashButtonBox, Gtk::PACK_SHRINK, 4);
    hBox->queue_draw ();

    fileBrowser->setScrollPosition (hScrollPos[lastScrollPos], vScrollPos[lastScrollPos]);

    for (int i=0; i<8; i++)
        bCateg[i].block (false);
}

BrowserFilter FileCatalog::getFilter () {

    BrowserFilter filter;
    filter.showRanked[0] = bDir->get_active() || bUnRanked->get_active () || bTrash->get_active ();
    for (int i=1; i<=5; i++)
        filter.showRanked[i] = bDir->get_active() || bRank[i-1]->get_active () || bTrash->get_active ();
    filter.showTrash = bDir->get_active() || bTrash->get_active ();
    filter.showNotTrash = !bTrash->get_active ();
    if (!filterPanel)
		filter.exifFilterEnabled = false;
	else {
		if (!hasValidCurrentEFS)
			filter.exifFilter = dirEFS;
		else
			filter.exifFilter = currentEFS;
		filter.exifFilterEnabled = filterPanel->isEnabled ();
	}
    return filter;
}

void FileCatalog::filterChanged () {
    
    fileBrowser->applyFilter (getFilter());   
}

int FileCatalog::reparseDirectory () {

    if (selectedDirectory=="")
        return 0;

    if (!Glib::file_test (selectedDirectory, Glib::FILE_TEST_IS_DIR)) {
        closeDir ();
		return 0;
	}
	
	std::vector<Glib::ustring> nfileNameList = getFileList ();

	// check if a thumbnailed file has been deleted
	const std::vector<ThumbBrowserEntryBase*>& t = fileBrowser->getEntries ();
	std::vector<Glib::ustring> fileNamesToDel;
	for (int i=0; i<t.size(); i++) 
		if (!Glib::file_test (t[i]->filename, Glib::FILE_TEST_EXISTS))
			fileNamesToDel.push_back (t[i]->filename);
	for (int i=0; i<fileNamesToDel.size(); i++) {
		delete fileBrowser->delEntry (fileNamesToDel[i]);
		cacheMgr.deleteEntry (fileNamesToDel[i]);
	}

	// check if a new file has been added
	for (int i=0; i<nfileNameList.size(); i++) {
		bool found = false;
		for (int j=0; j<fileNameList.size(); j++)
			if (nfileNameList[i]==fileNameList[j]) {
				found = true;
				break;
			}
		if (!found) {
			previewLoader.stop ();
			checkAndAddFile (Gio::File::create_for_parse_name (nfileNameList[i]));
            _refreshProgressBar ();
			previewLoader.process ();
		}
	}

	fileNameList = nfileNameList;
	return 1;
}

#ifdef _WIN32
void FileCatalog::winDirChanged () {

    checkCounter = 0;
}
#endif

void FileCatalog::on_dir_changed (const Glib::RefPtr<Gio::File>& file, const Glib::RefPtr<Gio::File>& other_file, Gio::FileMonitorEvent event_type, bool internal) {

    if (!internal)
        gdk_threads_enter();

    if (event_type == Gio::FILE_MONITOR_EVENT_CREATED || event_type == Gio::FILE_MONITOR_EVENT_DELETED || event_type == Gio::FILE_MONITOR_EVENT_CHANGED) 
		reparseDirectory ();

	if (!internal)
        gdk_threads_leave();
}

void FileCatalog::checkAndAddFile (Glib::RefPtr<Gio::File> file) {

    if (!file)
        return;
    Glib::RefPtr<Gio::FileInfo> info = safe_query_file_info(file);
    if (info && info->get_file_type() != Gio::FILE_TYPE_DIRECTORY && (!info->is_hidden() || !options.fbShowHidden)) {
        int lastdot = info->get_name().find_last_of ('.');
        if (options.is_extention_enabled(lastdot!=Glib::ustring::npos ? info->get_name().substr (lastdot+1) : "")){
						previewLoader.add (DirEntry (file->get_parse_name()));
            previewsToLoad++;
				}
    }
}

void FileCatalog::addAndOpenFile (const Glib::ustring& fname) {

    Glib::RefPtr<Gio::File> file = Gio::File::create_for_path (fname);
    if (!file)
        return;
    Glib::RefPtr<Gio::FileInfo> info = safe_query_file_info(file);
    int lastdot = info->get_name().find_last_of ('.');
    if (options.is_extention_enabled(lastdot!=Glib::ustring::npos ? info->get_name().substr (lastdot+1) : "")){
        // if supported, load thumbnail first
        Thumbnail* tmb = cacheMgr.getEntry (file->get_parse_name());
        if (tmb) {
            FileBrowserEntry* entry = new FileBrowserEntry (tmb, file->get_parse_name());
  	        previewReady (entry);
            // open the file
            FCOIParams* params = new FCOIParams;
            params->catalog = this;
            params->tmb.push_back (tmb);
            tmb->increaseRef ();
            g_idle_add (fcopenimg, params);
        }
    }
}

void FileCatalog::emptyTrash () {
    
    const std::vector<ThumbBrowserEntryBase*> t = fileBrowser->getEntries ();
    std::vector<FileBrowserEntry*> toDel;
    for (int i=0; i<t.size(); i++)
        if (((FileBrowserEntry*)t[i])->thumbnail->getStage()==1)
            toDel.push_back (((FileBrowserEntry*)t[i]));
    deleteRequested (toDel);
}

void FileCatalog::zoomIn () {

    bool pLoad = previewLoader.runs();
    if (pLoad)
        previewLoader.stop ();
        
    fileBrowser->zoomIn ();
        
    if (pLoad)
        previewLoader.process ();
}
void FileCatalog::zoomOut () {

    bool pLoad = previewLoader.runs();
    if (pLoad)
        previewLoader.stop ();
        
    fileBrowser->zoomOut ();
        
    if (pLoad)
        previewLoader.process ();
}
void FileCatalog::refreshEditedState (const std::set<Glib::ustring>& efiles) {

    editedFiles = efiles;
    fileBrowser->refreshEditedState (efiles);
}

void FileCatalog::selectionChanged (std::vector<Thumbnail*> tbe) {

    if (fslistener)
        fslistener->selectionChanged (tbe);
}

void FileCatalog::exifFilterChanged () {

	currentEFS = filterPanel->getFilter ();
    hasValidCurrentEFS = true;
    fileBrowser->applyFilter (getFilter ());
}

void FileCatalog::setFilterPanel (FilterPanel* fpanel) { 

	filterPanel = fpanel; 
	filterPanel->set_sensitive (false);
	filterPanel->setFilterPanelListener (this);
}
