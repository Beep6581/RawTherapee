/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Copyright (c) 2011 Oliver Duis <www.oliverduis.de>
 *  Copyright (c) 2011 Michael Ezra <www.michaelezra.com>
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
#include <filebrowser.h>
#include <glibmm.h>
#include <options.h>
#include <multilangmgr.h>
#include <clipboard.h>
#include <profilestore.h>
#include <procparamchangers.h>
#include <dfmanager.h>
#include <ffmanager.h>

extern Options options;

FileBrowser::FileBrowser () 
    : tbl(NULL),numFiltered(0) {

    fbih = new FileBrowserIdleHelper;
    fbih->fbrowser = this;
    fbih->destroyed = false;
    fbih->pending = 0;

  //  profileStore.parseProfiles ();

    signal_style_changed().connect( sigc::mem_fun(*this, &FileBrowser::styleChanged) );
    
    int p = 0;
    pmenu = new Gtk::Menu ();
    pmenu->attach (*Gtk::manage(open = new Gtk::MenuItem (M("FILEBROWSER_POPUPOPEN"))), 0, 1, p, p+1); p++;
    pmenu->attach (*Gtk::manage(develop = new Gtk::ImageMenuItem (M("FILEBROWSER_POPUPPROCESS"))), 0, 1, p, p+1); p++;
    develop->set_image(*Gtk::manage(new Gtk::Image (argv0+"/images/processing.png")));

    pmenu->attach (*Gtk::manage(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;
    pmenu->attach (*Gtk::manage(selall = new Gtk::MenuItem (M("FILEBROWSER_POPUPSELECTALL"))), 0, 1, p, p+1); p++;

    /***********************
     * rank
     ***********************/
    if (options.menuGroupRank){
    	pmenu->attach (*Gtk::manage(menuRank = new Gtk::MenuItem (M("FILEBROWSER_POPUPRANK"))), 0, 1, p, p+1); p++;
    	Gtk::Menu* submenuRank = Gtk::manage (new Gtk::Menu ());
    	submenuRank->attach (*Gtk::manage(rank[0] = new Gtk::MenuItem (M("FILEBROWSER_POPUPUNRANK"))), 0, 1, p, p+1); p++;
    	for (int i=1; i<=5; i++){
    		submenuRank->attach (*Gtk::manage(rank[i] = new Gtk::MenuItem (M(Glib::ustring::compose("%1%2","FILEBROWSER_POPUPRANK",i)))), 0, 1, p, p+1); p++;
    	}
    	submenuRank->show_all ();
		menuRank->set_submenu (*submenuRank);
    }
    else{
        pmenu->attach (*Gtk::manage(rank[0] = new Gtk::MenuItem (M("FILEBROWSER_POPUPUNRANK"))), 0, 1, p, p+1); p++;
    	for (int i=1; i<=5; i++){
    		pmenu->attach (*Gtk::manage(rank[i] = new Gtk::MenuItem (M(Glib::ustring::compose("%1%2","FILEBROWSER_POPUPRANK",i)))), 0, 1, p, p+1); p++;
    	}
        pmenu->attach (*Gtk::manage(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;
    }

    if (!options.menuGroupRank || !options.menuGroupLabel) // separate Rank and Color Labels if either is not grouped
    	 pmenu->attach (*Gtk::manage(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;

    /***********************
     * color labels
     ***********************/
    if (options.menuGroupLabel){
    	pmenu->attach (*Gtk::manage(menuLabel = new Gtk::MenuItem (M("FILEBROWSER_POPUPCOLORLABEL"))), 0, 1, p, p+1); p++;
    	Gtk::Menu* submenuLabel = Gtk::manage (new Gtk::Menu ());

    	for (int i=0; i<=5; i++){
    		submenuLabel->attach (*Gtk::manage(colorlabel[i] = new Gtk::ImageMenuItem (M(Glib::ustring::compose("%1%2","FILEBROWSER_POPUPCOLORLABEL",i)))), 0, 1, p, p+1); p++;
    	}
    	submenuLabel->show_all ();
		menuLabel->set_submenu (*submenuLabel);
    }
    else{
    	for (int i=0; i<=5; i++){
    		pmenu->attach (*Gtk::manage(colorlabel[i] = new Gtk::ImageMenuItem (M(Glib::ustring::compose("%1%2","FILEBROWSER_POPUPCOLORLABEL",i)))), 0, 1, p, p+1); p++;
    	}
    }
    for (int i=1; i<=5; i++){//set color label images
    	colorlabel[i]->set_image(*Gtk::manage(new Gtk::Image (Glib::ustring::compose("%1%2%3%4",argv0,"/images/clabel",i,".png"))));
    }
        
    pmenu->attach (*Gtk::manage(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;

    /***********************
     * File Operations
     * *********************/
    if (options.menuGroupFileOperations){
    	pmenu->attach (*Gtk::manage(menuFileOperations = new Gtk::MenuItem (M("FILEBROWSER_POPUPFILEOPERATIONS"))), 0, 1, p, p+1); p++;
    	Gtk::Menu* submenuFileOperations = Gtk::manage (new Gtk::Menu ());

    	submenuFileOperations->attach (*Gtk::manage(trash = new Gtk::MenuItem (M("FILEBROWSER_POPUPTRASH"))), 0, 1, p, p+1); p++;
    	submenuFileOperations->attach (*Gtk::manage(untrash = new Gtk::MenuItem (M("FILEBROWSER_POPUPUNTRASH"))), 0, 1, p, p+1); p++;
    	submenuFileOperations->attach (*Gtk::manage(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;
    	submenuFileOperations->attach (*Gtk::manage(rename = new Gtk::MenuItem (M("FILEBROWSER_POPUPRENAME"))), 0, 1, p, p+1); p++;
    	submenuFileOperations->attach (*Gtk::manage(remove = new Gtk::MenuItem (M("FILEBROWSER_POPUPREMOVE"))), 0, 1, p, p+1); p++;
    	submenuFileOperations->attach (*Gtk::manage(removeInclProc = new Gtk::MenuItem (M("FILEBROWSER_POPUPREMOVEINCLPROC"))), 0, 1, p, p+1); p++;
    	submenuFileOperations->attach (*Gtk::manage(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;
    	submenuFileOperations->attach (*Gtk::manage(copyTo = new Gtk::MenuItem (M("FILEBROWSER_POPUPCOPYTO"))), 0, 1, p, p+1); p++;
    	submenuFileOperations->attach (*Gtk::manage(moveTo = new Gtk::MenuItem (M("FILEBROWSER_POPUPMOVETO"))), 0, 1, p, p+1); p++;

    	submenuFileOperations->show_all ();
    	menuFileOperations->set_submenu (*submenuFileOperations);
    }
    else{
        pmenu->attach (*Gtk::manage(trash = new Gtk::MenuItem (M("FILEBROWSER_POPUPTRASH"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(untrash = new Gtk::MenuItem (M("FILEBROWSER_POPUPUNTRASH"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(rename = new Gtk::MenuItem (M("FILEBROWSER_POPUPRENAME"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(remove = new Gtk::MenuItem (M("FILEBROWSER_POPUPREMOVE"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(removeInclProc = new Gtk::MenuItem (M("FILEBROWSER_POPUPREMOVEINCLPROC"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(copyTo = new Gtk::MenuItem (M("FILEBROWSER_POPUPCOPYTO"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(moveTo = new Gtk::MenuItem (M("FILEBROWSER_POPUPMOVETO"))), 0, 1, p, p+1); p++;
    }

    pmenu->attach (*Gtk::manage(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;

    /***********************
     * Profile Operations
     * *********************/
    if (options.menuGroupProfileOperations){
    	pmenu->attach (*Gtk::manage(menuProfileOperations = new Gtk::ImageMenuItem (M("FILEBROWSER_POPUPPROFILEOPERATIONS"))), 0, 1, p, p+1); p++;
    	menuProfileOperations->set_image(*Gtk::manage(new Gtk::Image (argv0+"/images/logoicon_wind_16.png")));

    	Gtk::Menu* submenuProfileOperations = Gtk::manage (new Gtk::Menu ());

    	submenuProfileOperations->attach (*Gtk::manage(copyprof = new Gtk::MenuItem (M("FILEBROWSER_COPYPROFILE"))), 0, 1, p, p+1); p++;
    	submenuProfileOperations->attach (*Gtk::manage(pasteprof = new Gtk::MenuItem (M("FILEBROWSER_PASTEPROFILE"))), 0, 1, p, p+1); p++;
    	submenuProfileOperations->attach (*Gtk::manage(partpasteprof = new Gtk::MenuItem (M("FILEBROWSER_PARTIALPASTEPROFILE"))), 0, 1, p, p+1); p++;
    	submenuProfileOperations->attach (*Gtk::manage(applyprof = new Gtk::MenuItem (M("FILEBROWSER_APPLYPROFILE"))), 0, 1, p, p+1); p++;
    	submenuProfileOperations->attach (*Gtk::manage(applypartprof = new Gtk::MenuItem (M("FILEBROWSER_APPLYPROFILE_PARTIAL"))), 0, 1, p, p+1); p++;
    	submenuProfileOperations->attach (*Gtk::manage(clearprof = new Gtk::MenuItem (M("FILEBROWSER_CLEARPROFILE"))), 0, 1, p, p+1); p++;

    	submenuProfileOperations->show_all ();
    	menuProfileOperations->set_submenu (*submenuProfileOperations);
    }
    else{
        pmenu->attach (*Gtk::manage(copyprof = new Gtk::MenuItem (M("FILEBROWSER_COPYPROFILE"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(pasteprof = new Gtk::MenuItem (M("FILEBROWSER_PASTEPROFILE"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(partpasteprof = new Gtk::MenuItem (M("FILEBROWSER_PARTIALPASTEPROFILE"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(applyprof = new Gtk::MenuItem (M("FILEBROWSER_APPLYPROFILE"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(applypartprof = new Gtk::MenuItem (M("FILEBROWSER_APPLYPROFILE_PARTIAL"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(clearprof = new Gtk::MenuItem (M("FILEBROWSER_CLEARPROFILE"))), 0, 1, p, p+1); p++;
    }


    pmenu->attach (*Gtk::manage(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;
    pmenu->attach (*Gtk::manage(menuDF = new Gtk::MenuItem (M("FILEBROWSER_DARKFRAME"))), 0, 1, p, p+1); p++;
    pmenu->attach (*Gtk::manage(menuFF = new Gtk::MenuItem (M("FILEBROWSER_FLATFIELD"))), 0, 1, p, p+1); p++;

    pmenu->attach (*Gtk::manage(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;
    pmenu->attach (*Gtk::manage(cachemenu = new Gtk::MenuItem (M("FILEBROWSER_CACHE"))), 0, 1, p, p+1); p++;

    pmenu->show_all ();

    /***********************
     * Accelerators
     * *********************/
    pmaccelgroup = Gtk::AccelGroup::create ();
    pmenu->set_accel_group (pmaccelgroup);
    selall->add_accelerator ("activate", pmenu->get_accel_group(), GDK_a, Gdk::CONTROL_MASK, Gtk::ACCEL_VISIBLE);
    trash->add_accelerator ("activate", pmenu->get_accel_group(), GDK_Delete, (Gdk::ModifierType)0, Gtk::ACCEL_VISIBLE);
    untrash->add_accelerator ("activate", pmenu->get_accel_group(), GDK_Delete, Gdk::SHIFT_MASK, Gtk::ACCEL_VISIBLE);
    develop->add_accelerator ("activate", pmenu->get_accel_group(), GDK_Q, Gdk::CONTROL_MASK, Gtk::ACCEL_VISIBLE);
    copyprof->add_accelerator ("activate", pmenu->get_accel_group(), GDK_C, Gdk::CONTROL_MASK, Gtk::ACCEL_VISIBLE);
    pasteprof->add_accelerator ("activate", pmenu->get_accel_group(), GDK_V, Gdk::CONTROL_MASK, Gtk::ACCEL_VISIBLE);
    partpasteprof->add_accelerator ("activate", pmenu->get_accel_group(), GDK_V, Gdk::CONTROL_MASK | Gdk::SHIFT_MASK, Gtk::ACCEL_VISIBLE);

    open->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), open));    
    for (int i=0; i<6; i++)
        rank[i]->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), rank[i]));
    for (int i=0; i<6; i++)
    	colorlabel[i]->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), colorlabel[i]));

    trash->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), trash));    
    untrash->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), untrash));    
    develop->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), develop));    
    rename->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), rename));    
    remove->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), remove));    
    removeInclProc->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), removeInclProc));
    selall->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), selall));
    copyTo->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), copyTo));
    moveTo->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), moveTo));
    copyprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), copyprof));    
    pasteprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), pasteprof));    
    partpasteprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), partpasteprof));    
    applyprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), applyprof));    
    applypartprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), applypartprof));
    clearprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), clearprof));
    cachemenu->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), cachemenu));
}

FileBrowser::~FileBrowser ()
{
    delete pmenu;
}

void FileBrowser::rightClicked (ThumbBrowserEntryBase* entry) {

    trash->set_sensitive (false);
    untrash->set_sensitive (false);
    for (int i=0; i<selected.size(); i++) 
        if (((FileBrowserEntry*)selected[i])->thumbnail->getStage()==1) {
            untrash->set_sensitive (true);
            break;
        }
    for (int i=0; i<selected.size(); i++) 
        if (((FileBrowserEntry*)selected[i])->thumbnail->getStage()==0) {
            trash->set_sensitive (true);
            break;
        }

    pasteprof->set_sensitive (clipboard.hasProcParams());
    partpasteprof->set_sensitive (clipboard.hasProcParams());
    copyprof->set_sensitive (selected.size()==1);
    clearprof->set_sensitive (selected.size()>0);

    // submenu applmenu
    int p = 0;
    Gtk::Menu* applmenu = Gtk::manage (new Gtk::Menu ());
    std::vector<Glib::ustring> profnames = profileStore.getProfileNames ();
    for (int i=0; i<profnames.size(); i++) {
        Gtk::MenuItem* mi = Gtk::manage (new Gtk::MenuItem (profnames[i]));
        applmenu->attach (*mi, 0, 1, p, p+1); p++;
        mi->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::applyMenuItemActivated), profnames[i]));    
        mi->show ();
    }
    applyprof->set_submenu (*applmenu);

    // submenu applpartmenu
    p = 0;
    Gtk::Menu* applpartmenu = Gtk::manage (new Gtk::Menu ());
    //std::vector<Glib::ustring> profnames = profileStore.getProfileNames (); // this is already created for submenu applmenu above
    for (int i=0; i<profnames.size(); i++) {
	    Gtk::MenuItem* mi = Gtk::manage (new Gtk::MenuItem (profnames[i]));
	    applpartmenu->attach (*mi, 0, 1, p, p+1); p++;
	    mi->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::applyPartialMenuItemActivated), profnames[i]));
	    mi->show ();
    }
    applypartprof->set_submenu (*applpartmenu);

    // submenuDF
    p = 0;
    Gtk::Menu* submenuDF = Gtk::manage (new Gtk::Menu ());
    submenuDF->attach (*Gtk::manage(selectDF = new Gtk::MenuItem (M("FILEBROWSER_SELECTDARKFRAME"))), 0, 1, p, p+1); p++;
    submenuDF->attach (*Gtk::manage(autoDF = new Gtk::MenuItem (M("FILEBROWSER_AUTODARKFRAME"))), 0, 1, p, p+1); p++;
    submenuDF->attach (*Gtk::manage(thisIsDF = new Gtk::MenuItem (M("FILEBROWSER_MOVETODARKFDIR"))), 0, 1, p, p+1); p++;
    selectDF->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), selectDF));
	autoDF->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), autoDF));
	thisIsDF->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated),thisIsDF ));
	submenuDF->show_all ();
	menuDF->set_submenu (*submenuDF);

	// submenuFF
	p = 0;
	Gtk::Menu* submenuFF = Gtk::manage (new Gtk::Menu ());
	submenuFF->attach (*Gtk::manage(selectFF = new Gtk::MenuItem (M("FILEBROWSER_SELECTFLATFIELD"))), 0, 1, p, p+1); p++;
	submenuFF->attach (*Gtk::manage(autoFF = new Gtk::MenuItem (M("FILEBROWSER_AUTOFLATFIELD"))), 0, 1, p, p+1); p++;
	submenuFF->attach (*Gtk::manage(thisIsFF = new Gtk::MenuItem (M("FILEBROWSER_MOVETOFLATFIELDDIR"))), 0, 1, p, p+1); p++;
	selectFF->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), selectFF));
	autoFF->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), autoFF));
	thisIsFF->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated),thisIsFF ));
	submenuFF->show_all ();
	menuFF->set_submenu (*submenuFF);

    // build cache sub menu
    p = 0;
	Gtk::Menu* cachesubmenu = Gtk::manage (new Gtk::Menu ());
	cachesubmenu->attach (*Gtk::manage(clearFromCache = new Gtk::MenuItem (M("FILEBROWSER_CACHECLEARFROMPARTIAL"))), 0, 1, p, p+1); p++;
	cachesubmenu->attach (*Gtk::manage(clearFromCacheFull = new Gtk::MenuItem (M("FILEBROWSER_CACHECLEARFROMFULL"))), 0, 1, p, p+1); p++;
    clearFromCache->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), clearFromCache));
    clearFromCacheFull->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), clearFromCacheFull));
    cachesubmenu->show_all ();
    cachemenu->set_submenu (*cachesubmenu);

    pmenu->popup (3, this->eventTime);
}

void FileBrowser::doubleClicked (ThumbBrowserEntryBase* entry) {

    if (tbl && entry) {
        std::vector<Thumbnail*> entries;
        entries.push_back (((FileBrowserEntry*)entry)->thumbnail);
        tbl->openRequested (entries);
    }
}

struct addparams {
    FileBrowserIdleHelper* fbih;
    FileBrowserEntry* entry;
};

int AddEntryUIThread (void* data) {
    
    addparams* ap = (addparams*) data;
    FileBrowserIdleHelper* fbih = ap->fbih;

    gdk_threads_enter();

    if (fbih->destroyed) {
        if (fbih->pending == 1)
            delete fbih;
        else    
            fbih->pending--;
        delete ap->entry;
        delete ap;
        gdk_threads_leave ();
        return 0;
    }

    ap->fbih->fbrowser->addEntry_ (ap->entry);
    delete ap;
    fbih->pending--;
    gdk_threads_leave();
    return 0;
}

void FileBrowser::addEntry (FileBrowserEntry* entry) {

    fbih->pending++;
    entry->setParent (this);
    addparams* ap = new addparams;
    ap->fbih = fbih;
    ap->entry = entry;
    g_idle_add (AddEntryUIThread, ap);
}

void FileBrowser::addEntry_ (FileBrowserEntry* entry) {

    entry->selected = false;
    entry->drawable = false;
    entry->framed = editedFiles.find (entry->filename)!=editedFiles.end();
        
    // add button set to the thumbbrowserentry
    entry->addButtonSet (new FileThumbnailButtonSet (entry));
    entry->getThumbButtonSet()->setRank (entry->thumbnail->getRank());
    entry->getThumbButtonSet()->setColorLabel (entry->thumbnail->getColorLabel());
    entry->getThumbButtonSet()->setInTrash (entry->thumbnail->getStage()==1);
    entry->getThumbButtonSet()->setButtonListener (this);
    entry->resize (getCurrentThumbSize());

    // find place in abc order
	{
		// TODO: Check for Linux
		#ifdef WIN32
		Glib::Mutex::Lock lock(entryMutex);
		#endif

    std::vector<ThumbBrowserEntryBase*>::iterator i = fd.begin();
    while (i!=fd.end() && *entry < *((FileBrowserEntry*)*i))
        i++;
        
    fd.insert (i, entry);

    initEntry (entry);
	}
    redraw ();
}

FileBrowserEntry* FileBrowser::delEntry (const Glib::ustring& fname) {
	// TODO: Check for Linux
	#ifdef WIN32
	Glib::Mutex::Lock lock(entryMutex);
	#endif

    for (std::vector<ThumbBrowserEntryBase*>::iterator i=fd.begin(); i!=fd.end(); i++) 
        if ((*i)->filename==fname) {
            ThumbBrowserEntryBase* entry = *i;
            entry->selected = false;
            fd.erase (i);
            std::vector<ThumbBrowserEntryBase*>::iterator j = std::find (selected.begin(), selected.end(), entry);
            if (j!=selected.end()) {
                selected.erase (j);
                notifySelectionListener ();
            }
            if (lastClicked==entry)
                lastClicked = NULL;
            redraw ();
            return (FileBrowserEntry*)entry;
        }
    return NULL;
}

FileBrowserEntry* FileBrowser::findEntry (const Glib::ustring& fname) {
	// TODO: Check for Linux
	#ifdef WIN32
	Glib::Mutex::Lock lock(entryMutex);
	#endif

    for (std::vector<ThumbBrowserEntryBase*>::iterator i=fd.begin(); i!=fd.end(); i++) 
        if ((*i)->filename==fname) 
            return (FileBrowserEntry*)*i;
    return NULL;
}

void FileBrowser::close () {
    if (fbih->pending)
        fbih->destroyed = true;
    else
        delete fbih;

    fbih = new FileBrowserIdleHelper;
    fbih->fbrowser = this;
    fbih->destroyed = false;
    fbih->pending = 0;

	{
		// TODO: Check for Linux
		#ifdef WIN32
		Glib::Mutex::Lock lock(entryMutex);
		#endif

        
		selected.clear ();
		notifySelectionListener ();

        // The listener merges parameters with old values, so delete afterwards
    for (int i=0; i<fd.size(); i++)
    {
        delete fd[i];
    }
    fd.clear ();
	}

    lastClicked = NULL;
}

void FileBrowser::menuItemActivated (Gtk::MenuItem* m) {

    std::vector<FileBrowserEntry*> mselected;
    for (int i=0; i<selected.size(); i++)
        mselected.push_back ((FileBrowserEntry*)selected[i]);

    if (!tbl || (m!=selall && mselected.size()==0) )
        return;

    for (int i=0; i<6; i++) 
        if (m==rank[i]) {
            rankingRequested (mselected, i);
            return;
        }
    for (int i=0; i<6; i++)
        if (m==colorlabel[i]) {
            colorlabelRequested (mselected, i);
            return;
        }
    if (m==open) {
        std::vector<Thumbnail*> entries;
        for (int i=0; i<mselected.size(); i++)
            entries.push_back (mselected[i]->thumbnail);
        tbl->openRequested (entries);
     }
    else if (m==remove)
        tbl->deleteRequested (mselected, false);
	else if (m==removeInclProc)
        tbl->deleteRequested (mselected, true);
    else if (m==trash) 
        toTrashRequested (mselected);
    else if (m==untrash) 
        fromTrashRequested (mselected);
    else if (m==develop)
        tbl->developRequested (mselected);
    else if (m==rename)
        tbl->renameRequested (mselected);
    else if (m==selall) {
        lastClicked = NULL;
		{
			// TODO: Check for Linux
			#ifdef WIN32
			Glib::Mutex::Lock lock(entryMutex);
			#endif

        selected.clear ();
        for (int i=0; i<fd.size(); i++)
            if (checkFilter (fd[i])) {
                fd[i]->selected = true;
                selected.push_back (fd[i]);
            }
		}
        queue_draw ();
        notifySelectionListener ();
    }
    else if( m==copyTo){
    	tbl->copyMoveRequested (mselected, false);
    }

    else if( m==moveTo){
       	tbl->copyMoveRequested (mselected, true);
    }

    else if (m==autoDF){
		for (int i=0; i<mselected.size(); i++){
			rtengine::procparams::ProcParams pp=mselected[i]->thumbnail->getProcParams();
			pp.raw.df_autoselect= true;
			pp.raw.dark_frame.clear();
			mselected[i]->thumbnail->setProcParams(pp,FILEBROWSER,false);
		}
    }else if (m==selectDF){
    	if( mselected.size() > 0 ){
    		rtengine::procparams::ProcParams pp=mselected[0]->thumbnail->getProcParams();
    		Gtk::FileChooserDialog fc("Dark Frame",Gtk::FILE_CHOOSER_ACTION_OPEN );
    		fc.add_button( Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    		fc.add_button( Gtk::StockID("gtk-apply"), Gtk::RESPONSE_APPLY);
    		if( pp.raw.dark_frame.empty())
    		   fc.set_current_folder( options.rtSettings.darkFramesPath );
    		else
    		   fc.set_filename( pp.raw.dark_frame );
    		if( fc.run() == Gtk::RESPONSE_APPLY ){
    			for (int i=0; i<mselected.size(); i++){
    				rtengine::procparams::ProcParams pp=mselected[i]->thumbnail->getProcParams();
					pp.raw.dark_frame= fc.get_filename();
					pp.raw.df_autoselect= false;
					mselected[i]->thumbnail->setProcParams(pp,FILEBROWSER,false);
				}
			}
    	}
    }else if( m==thisIsDF){
    	if( options.rtSettings.darkFramesPath.size() >0 && Gio::File::create_for_path(options.rtSettings.darkFramesPath)->query_exists() ){
			for (int i=0; i<mselected.size(); i++){
				Glib::RefPtr<Gio::File> file = Gio::File::create_for_path ( mselected[i]->filename );
				if( !file )continue;
				Glib::ustring destName = options.rtSettings.darkFramesPath+ "/" + file->get_basename();
				Glib::RefPtr<Gio::File> dest = Gio::File::create_for_path ( destName );
				file->move(  dest );
			}
			// Reinit cache
			rtengine::dfm.init( options.rtSettings.darkFramesPath );
    	}
    }
    else if (m==autoFF){
		for (int i=0; i<mselected.size(); i++){
			rtengine::procparams::ProcParams pp=mselected[i]->thumbnail->getProcParams();
			pp.raw.ff_AutoSelect= true;
			pp.raw.ff_file.clear();
			mselected[i]->thumbnail->setProcParams(pp,FILEBROWSER,false);
		}
    }
    else if (m==selectFF){
    	if( mselected.size() > 0 ){
    		rtengine::procparams::ProcParams pp=mselected[0]->thumbnail->getProcParams();
    		Gtk::FileChooserDialog fc("Flat Field",Gtk::FILE_CHOOSER_ACTION_OPEN );
    		fc.add_button( Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    		fc.add_button( Gtk::StockID("gtk-apply"), Gtk::RESPONSE_APPLY);
    		if( pp.raw.ff_file.empty())
    		   fc.set_current_folder( options.rtSettings.flatFieldsPath );
    		else
    		   fc.set_filename( pp.raw.ff_file );
    		if( fc.run() == Gtk::RESPONSE_APPLY ){
    			for (int i=0; i<mselected.size(); i++){
    				rtengine::procparams::ProcParams pp=mselected[i]->thumbnail->getProcParams();
					pp.raw.ff_file= fc.get_filename();
					pp.raw.ff_AutoSelect= false;
					mselected[i]->thumbnail->setProcParams(pp,FILEBROWSER,false);
				  }
			  }
    	}
    }
    else if( m==thisIsFF){
    	if( options.rtSettings.flatFieldsPath.size() >0 && Gio::File::create_for_path(options.rtSettings.flatFieldsPath)->query_exists() ){
			for (int i=0; i<mselected.size(); i++){
				Glib::RefPtr<Gio::File> file = Gio::File::create_for_path ( mselected[i]->filename );
				if( !file )continue;
				Glib::ustring destName = options.rtSettings.flatFieldsPath+ "/" + file->get_basename();
				Glib::RefPtr<Gio::File> dest = Gio::File::create_for_path ( destName );
				file->move(  dest );
			}
			// Reinit cache
			rtengine::ffm.init( options.rtSettings.flatFieldsPath );
    	}
    }
    else if (m==copyprof)
        copyProfile ();
    else if (m==pasteprof) 
        pasteProfile ();
    else if (m==partpasteprof) 
        partPasteProfile ();
    else if (m==clearprof) {
        for (int i=0; i<mselected.size(); i++) 
            mselected[i]->thumbnail->clearProcParams (FILEBROWSER);
        queue_draw ();
    }
	else if (m==clearFromCache) {
		for (int i=0; i<mselected.size(); i++)
			tbl->clearFromCacheRequested (mselected, false);
		//queue_draw ();
    }
	else if (m==clearFromCacheFull) {
		for (int i=0; i<mselected.size(); i++)
			tbl->clearFromCacheRequested (mselected, true);
		//queue_draw ();
    }
}

void FileBrowser::copyProfile () {

    if (selected.size()==1)
        clipboard.setProcParams (((FileBrowserEntry*)selected[0])->thumbnail->getProcParams());
}

void FileBrowser::pasteProfile () {

    std::vector<FileBrowserEntry*> mselected;
    for (int i=0; i<selected.size(); i++)
        mselected.push_back ((FileBrowserEntry*)selected[i]);

    if (!tbl || mselected.size()==0)
        return;

    for (int i=0; i<mselected.size(); i++) 
        mselected[i]->thumbnail->setProcParams (clipboard.getProcParams(), FILEBROWSER);
    
    queue_draw ();
}

void FileBrowser::partPasteProfile () {

    std::vector<FileBrowserEntry*> mselected;
    for (int i=0; i<selected.size(); i++)
        mselected.push_back ((FileBrowserEntry*)selected[i]);

    if (!tbl || mselected.size()==0)
        return;

    if (partialPasteDlg.run ()) {
    
        for (int i=0; i<mselected.size(); i++) {
            rtengine::procparams::ProcParams params = mselected[i]->thumbnail->getProcParams ();
            partialPasteDlg.applyPaste (&params, &clipboard.getProcParams());
            mselected[i]->thumbnail->setProcParams (params, FILEBROWSER);
        }

        queue_draw ();
    }
    partialPasteDlg.hide ();
}

void FileBrowser::openDefaultViewer (int destination) {
    if (selected.size()==1)
        ((FileBrowserEntry*)selected[0])->thumbnail->openDefaultViewer(destination);
}

bool FileBrowser::keyPressed (GdkEventKey* event) {

    if ((event->keyval==GDK_C || event->keyval==GDK_c) && event->state & GDK_CONTROL_MASK) {
        copyProfile ();
        return true;
    }
    else if ((event->keyval==GDK_V || event->keyval==GDK_v) && event->state & GDK_CONTROL_MASK && !(event->state & GDK_SHIFT_MASK)) {
        pasteProfile ();
        return true;
    }
    else if ((event->keyval==GDK_V || event->keyval==GDK_v) && event->state & GDK_CONTROL_MASK && event->state & GDK_SHIFT_MASK) {
        partPasteProfile ();
        return true;
    }
    else if (event->keyval==GDK_Delete && !(event->state & GDK_SHIFT_MASK)) {
        menuItemActivated (trash);
        return true;
    }
    else if (event->keyval==GDK_Delete && event->state & GDK_SHIFT_MASK) {
        menuItemActivated (untrash);
        return true;
    }
    else if ((event->keyval==GDK_Q || event->keyval==GDK_q) && event->state & GDK_CONTROL_MASK) {
        menuItemActivated (develop);
        return true;
    }
    else if ((event->keyval==GDK_A || event->keyval==GDK_a) && event->state & GDK_CONTROL_MASK) {
        menuItemActivated (selall);
        return true;
    }
    else if (event->keyval==GDK_F5) {
        int dest = 1;
        if (event->state & GDK_SHIFT_MASK)
            dest = 2;
        else if (event->state & GDK_CONTROL_MASK)
            dest = 3;

        openDefaultViewer (dest);
        return true;
    }
    else if (event->keyval==GDK_Page_Up) {
        scrollPage(GDK_SCROLL_UP);
        return true;
    }
    else if (event->keyval==GDK_Page_Down) {
        scrollPage(GDK_SCROLL_DOWN);
        return true;
    }
        
    return false;
}

void FileBrowser::applyMenuItemActivated (Glib::ustring ppname) {

    rtengine::procparams::ProcParams* pparams = profileStore.getProfile (ppname);
    if (pparams && selected.size()>0) {
        for (int i=0; i<selected.size(); i++) 
            ((FileBrowserEntry*)selected[i])->thumbnail->setProcParams (*pparams, FILEBROWSER);
        queue_draw ();
    }
}

void FileBrowser::applyPartialMenuItemActivated (Glib::ustring ppname) {

	if (!tbl || selected.size()==0)
		return;

	rtengine::procparams::ProcParams* pparams = profileStore.getProfile (ppname);

	if (pparams) {
		if (partialPasteDlg.run ()) {

			for (int i=0; i<selected.size(); i++) {
				rtengine::procparams::ProcParams params = ((FileBrowserEntry*)selected[i])->thumbnail->getProcParams ();
				partialPasteDlg.applyPaste (&params, pparams);
				((FileBrowserEntry*)selected[i])->thumbnail->setProcParams (params, FILEBROWSER);
			}
			queue_draw ();
		}
		partialPasteDlg.hide ();
	}
}

void FileBrowser::applyFilter (const BrowserFilter& filter) {

    this->filter = filter;

    // remove items not complying the filter from the selection
    bool selchanged = false;
    numFiltered=0;
	{
		// TODO: Check for Linux
		#ifdef WIN32
		Glib::Mutex::Lock lock(entryMutex);
		#endif

    for (int i=0; i<fd.size(); i++)
    	if(checkFilter (fd[i]))
    		numFiltered++;
    	else if (fd[i]->selected ) {
            fd[i]->selected = false;
            std::vector<ThumbBrowserEntryBase*>::iterator j = std::find (selected.begin(), selected.end(), fd[i]);
            selected.erase (j);
            if (lastClicked==fd[i])
                lastClicked = NULL;
            selchanged = true;
        }
	}

    if (selchanged)
        notifySelectionListener ();
    redraw ();
}

bool FileBrowser::checkFilter (ThumbBrowserEntryBase* entryb) { // true -> entry complies filter
    
    FileBrowserEntry* entry = (FileBrowserEntry*)entryb;
    // return false if basic filter settings are not satisfied
    if (filter.showRanked[entry->thumbnail->getRank()]==false || filter.showCLabeled[entry->thumbnail->getColorLabel()]==false || (entry->thumbnail->getStage()==1 && !filter.showTrash) || (entry->thumbnail->getStage()==0 && !filter.showNotTrash))
        return false;

    // return false is query is not satisfied
    if (filter.queryFileName.size()>0){
    	// check if image's FileName contains queryFileName (case insensitive)
    	// TODO should we provide case-sensitive search option via preferences?
    	Glib::ustring FileName;
    	FileName = Glib::path_get_basename (entry->thumbnail->getFileName());
    	FileName = FileName.uppercase();
    	//printf("FileBrowser::checkFilter FileName = '%s'; find() result= %i \n",FileName.c_str(), FileName.find(filter.queryFileName.uppercase()));
    	
    	if (FileName.find(filter.queryFileName.uppercase())==-1)
    		 return false;
    }

    // check exif filter
    const CacheImageData* cfs = entry->thumbnail->getCacheImageData();
    double tol = 0.01;
    double tol2 = 1e-8;
    
	if (!filter.exifFilterEnabled)
		return true;
	
	if (!cfs->exifValid)
		return (!filter.exifFilter.filterCamera || filter.exifFilter.cameras.count(cfs->camera)>0) 
			&& (!filter.exifFilter.filterLens || filter.exifFilter.lenses.count(cfs->lens)>0)
			&& (!filter.exifFilter.filterFiletype || filter.exifFilter.filetypes.count(cfs->filetype)>0);
		
    return 
         (!filter.exifFilter.filterShutter || (rtengine::ImageMetaData::shutterFromString(rtengine::ImageMetaData::shutterToString(cfs->shutter)) >= filter.exifFilter.shutterFrom-tol2 && rtengine::ImageMetaData::shutterFromString(rtengine::ImageMetaData::shutterToString(cfs->shutter)) <= filter.exifFilter.shutterTo+tol2))
      && (!filter.exifFilter.filterFNumber || (rtengine::ImageMetaData::apertureFromString(rtengine::ImageMetaData::apertureToString(cfs->fnumber)) >= filter.exifFilter.fnumberFrom-tol2 && rtengine::ImageMetaData::apertureFromString(rtengine::ImageMetaData::apertureToString(cfs->fnumber)) <= filter.exifFilter.fnumberTo+tol2))
      && (!filter.exifFilter.filterFocalLen || (cfs->focalLen >= filter.exifFilter.focalFrom-tol && cfs->focalLen <= filter.exifFilter.focalTo+tol))
	  && (!filter.exifFilter.filterISO     || (cfs->iso >= filter.exifFilter.isoFrom && cfs->iso <= filter.exifFilter.isoTo))
      && (!filter.exifFilter.filterCamera  || filter.exifFilter.cameras.count(cfs->camera)>0)
	  && (!filter.exifFilter.filterLens    || filter.exifFilter.lenses.count(cfs->lens)>0)
	    && (!filter.exifFilter.filterFiletype  || filter.exifFilter.filetypes.count(cfs->filetype)>0);
}

void FileBrowser::toTrashRequested (std::vector<FileBrowserEntry*> tbe) {

    for (int i=0; i<tbe.size(); i++) {
    	// try to load the last saved parameters from the cache or from the paramfile file
    	tbe[i]->thumbnail->createProcParamsForUpdate();  // this can execute customprofilebuilder to generate param file

    	// no need to notify listeners as item goes to trash, likely to be deleted

    	if (tbe[i]->thumbnail->getStage()==1)
            continue;               
        tbe[i]->thumbnail->setStage (1);
        if (tbe[i]->getThumbButtonSet()) {
            tbe[i]->getThumbButtonSet()->setRank (tbe[i]->thumbnail->getRank());
            tbe[i]->getThumbButtonSet()->setColorLabel (tbe[i]->thumbnail->getColorLabel());
            tbe[i]->getThumbButtonSet()->setInTrash (true);
            tbe[i]->thumbnail->updateCache(); // needed to save the rank to disk
        }
    }
    trash_changed().emit();
    applyFilter (filter);
}

void FileBrowser::fromTrashRequested (std::vector<FileBrowserEntry*> tbe) {

    for (int i=0; i<tbe.size(); i++) {
    	// if thumbnail was marked inTrash=true then param file must be there, no need to run customprofilebuilder

        if (tbe[i]->thumbnail->getStage()==0)
            continue;
        tbe[i]->thumbnail->setStage (0);
        if (tbe[i]->getThumbButtonSet()) {
            tbe[i]->getThumbButtonSet()->setRank (tbe[i]->thumbnail->getRank());
            tbe[i]->getThumbButtonSet()->setColorLabel (tbe[i]->thumbnail->getColorLabel());
            tbe[i]->getThumbButtonSet()->setInTrash (false);
            tbe[i]->thumbnail->updateCache(); // needed to save the rank to disk
        }
    }
    trash_changed().emit();
    applyFilter (filter);
}

void FileBrowser::rankingRequested (std::vector<FileBrowserEntry*> tbe, int rank) {

    for (int i=0; i<tbe.size(); i++) {
    	// try to load the last saved parameters from the cache or from the paramfile file
    	tbe[i]->thumbnail->createProcParamsForUpdate();  // this can execute customprofilebuilder to generate param file

    	// notify listeners TODO: should do this ONLY when params changed by customprofilebuilder?
    	tbe[i]->thumbnail->notifylisterners_procParamsChanged(FILEBROWSER);

        tbe[i]->thumbnail->setRank (rank);
        tbe[i]->thumbnail->updateCache(); // needed to save the rank to disk
        //TODO? - should update pparams instead?

        if (tbe[i]->getThumbButtonSet())
                tbe[i]->getThumbButtonSet()->setRank (tbe[i]->thumbnail->getRank());
    }
    applyFilter (filter);
}

void FileBrowser::colorlabelRequested (std::vector<FileBrowserEntry*> tbe, int colorlabel) {

    for (int i=0; i<tbe.size(); i++) {
    	// try to load the last saved parameters from the cache or from the paramfile file
    	tbe[i]->thumbnail->createProcParamsForUpdate();  // this can execute customprofilebuilder to generate param file

    	// notify listeners TODO: should do this ONLY when params changed by customprofilebuilder?
    	tbe[i]->thumbnail->notifylisterners_procParamsChanged(FILEBROWSER);

        tbe[i]->thumbnail->setColorLabel (colorlabel);
        tbe[i]->thumbnail->updateCache(); // needed to save the colorlabel to disk
        //TODO? - should update pparams instead?
        if (tbe[i]->getThumbButtonSet())
                tbe[i]->getThumbButtonSet()->setColorLabel (tbe[i]->thumbnail->getColorLabel());
    }
    applyFilter (filter);
}

void FileBrowser::buttonPressed (LWButton* button, int actionCode, void* actionData) {

    if (actionCode>=0 && actionCode<=5) { // rank
        std::vector<FileBrowserEntry*> tbe;
        tbe.push_back ((FileBrowserEntry*)actionData);
        rankingRequested (tbe, actionCode);
    }
    else if (actionCode==6 && tbl) { // to processing queue
        std::vector<FileBrowserEntry*> tbe;
        tbe.push_back ((FileBrowserEntry*)actionData);
        tbl->developRequested (tbe);
    }
    else if (actionCode==7) { // to trash / undelete
        std::vector<FileBrowserEntry*> tbe;
        FileBrowserEntry* entry = (FileBrowserEntry*)actionData;
        tbe.push_back (entry);
        if (entry->thumbnail->getStage()==0)
            toTrashRequested (tbe);
        else
            fromTrashRequested (tbe);
    }
}

void FileBrowser::openNextImage () {
	// TODO: Check for Linux
	#ifdef WIN32
	Glib::Mutex::Lock lock(entryMutex);
	#endif

    if (fd.size()>0) {
        for (int i=fd.size()-1; i>=0; i--)
            if (editedFiles.find (fd[i]->filename)!=editedFiles.end()) 
                if (i<fd.size()-1 && tbl) {
                    std::vector<Thumbnail*> entries;
                    entries.push_back (((FileBrowserEntry*)fd[i+1])->thumbnail);
                    tbl->openRequested (entries);
                    return;
                }
        if (tbl) {
            std::vector<Thumbnail*> entries;
            entries.push_back (((FileBrowserEntry*)fd[0])->thumbnail);
            tbl->openRequested (entries);
        }
    }
}

void FileBrowser::openPrevImage () {
	// TODO: Check for Linux
	#ifdef WIN32
	Glib::Mutex::Lock lock(entryMutex);
	#endif

    if (fd.size()>0) {
        for (int i=0; i<fd.size(); i++)
            if (editedFiles.find (fd[i]->filename)!=editedFiles.end()) 
                if (i>0 && tbl) {
                    std::vector<Thumbnail*> entries;
                    entries.push_back (((FileBrowserEntry*)fd[i-1])->thumbnail);
                    tbl->openRequested (entries);
                    return;
                }
        if (tbl) {
            std::vector<Thumbnail*> entries;
            entries.push_back (((FileBrowserEntry*)fd[fd.size()-1])->thumbnail);
            tbl->openRequested (entries);
        }
    }
}

int redrawtb (void* data) {

    ((FileBrowser*)data)->_thumbRearrangementNeeded ();
    return 0;
}

void FileBrowser::_thumbRearrangementNeeded () {

    refreshThumbImages ();
}

void FileBrowser::thumbRearrangementNeeded () {

    g_idle_add (redrawtb, this);
}
void FileBrowser::selectionChanged () {

    notifySelectionListener ();
}

void FileBrowser::notifySelectionListener () {

    if (tbl) {
        std::vector<Thumbnail*> thm;
        for (int i=0; i<selected.size(); i++)
            thm.push_back (((FileBrowserEntry*)selected[i])->thumbnail);
        tbl->selectionChanged (thm);
    }    
}

void FileBrowser::redrawNeeded (LWButton* button) {
    
    queue_draw ();
}
FileBrowser::type_trash_changed FileBrowser::trash_changed () {
    return m_trash_changed;
}
