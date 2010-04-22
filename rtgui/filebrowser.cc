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
#include <filebrowser.h>
#include <glibmm.h>
#include <options.h>
#include <multilangmgr.h>
#include <clipboard.h>
#include <profilestore.h>
#include <procparamchangers.h>

FileBrowser::FileBrowser () 
    : tbl(NULL) {

    fbih = new FileBrowserIdleHelper;
    fbih->fbrowser = this;
    fbih->destroyed = false;
    fbih->pending = 0;

    profileStore.parseProfiles ();

    signal_style_changed().connect( sigc::mem_fun(*this, &FileBrowser::styleChanged) );
    
    int p = 0;
    pmenu = new Gtk::Menu ();
    pmenu->attach (*(open = new Gtk::MenuItem (M("FILEBROWSER_POPUPOPEN"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(develop = new Gtk::MenuItem (M("FILEBROWSER_POPUPPROCESS"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;
    pmenu->attach (*(rank[0] = new Gtk::MenuItem (M("FILEBROWSER_POPUPUNRANK"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(rank[1] = new Gtk::MenuItem (M("FILEBROWSER_POPUPRANK1"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(rank[2] = new Gtk::MenuItem (M("FILEBROWSER_POPUPRANK2"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(rank[3] = new Gtk::MenuItem (M("FILEBROWSER_POPUPRANK3"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(rank[4] = new Gtk::MenuItem (M("FILEBROWSER_POPUPRANK4"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(rank[5] = new Gtk::MenuItem (M("FILEBROWSER_POPUPRANK5"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;
    pmenu->attach (*(trash = new Gtk::MenuItem (M("FILEBROWSER_POPUPTRASH"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(untrash = new Gtk::MenuItem (M("FILEBROWSER_POPUPUNTRASH"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;
    pmenu->attach (*(rename = new Gtk::MenuItem (M("FILEBROWSER_POPUPRENAME"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(remove = new Gtk::MenuItem (M("FILEBROWSER_POPUPREMOVE"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;
    pmenu->attach (*(selall = new Gtk::MenuItem (M("FILEBROWSER_POPUPSELECTALL"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;
    pmenu->attach (*(copyprof = new Gtk::MenuItem (M("FILEBROWSER_COPYPROFILE"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(pasteprof = new Gtk::MenuItem (M("FILEBROWSER_PASTEPROFILE"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(partpasteprof = new Gtk::MenuItem (M("FILEBROWSER_PARTIALPASTEPROFILE"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(applyprof = new Gtk::MenuItem (M("FILEBROWSER_APPLYPROFILE"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(clearprof = new Gtk::MenuItem (M("FILEBROWSER_CLEARPROFILE"))), 0, 1, p, p+1); p++;
    pmenu->show_all ();

    pmaccelgroup = Gtk::AccelGroup::create ();
    pmenu->set_accel_group (pmaccelgroup);
    selall->add_accelerator ("activate", pmenu->get_accel_group(), GDK_a, Gdk::CONTROL_MASK, Gtk::ACCEL_VISIBLE);
    trash->add_accelerator ("activate", pmenu->get_accel_group(), GDK_Delete, (Gdk::ModifierType)0, Gtk::ACCEL_VISIBLE);
    untrash->add_accelerator ("activate", pmenu->get_accel_group(), GDK_Delete, Gdk::SHIFT_MASK, Gtk::ACCEL_VISIBLE);
    develop->add_accelerator ("activate", pmenu->get_accel_group(), GDK_Q, Gdk::CONTROL_MASK, Gtk::ACCEL_VISIBLE);
    copyprof->add_accelerator ("activate", pmenu->get_accel_group(), GDK_C, Gdk::CONTROL_MASK, Gtk::ACCEL_VISIBLE);
    pasteprof->add_accelerator ("activate", pmenu->get_accel_group(), GDK_V, Gdk::CONTROL_MASK, Gtk::ACCEL_VISIBLE);
    partpasteprof->add_accelerator ("activate", pmenu->get_accel_group(), GDK_V, Gdk::CONTROL_MASK | Gdk::SHIFT_MASK, Gtk::ACCEL_VISIBLE);

    profmenu = new Gtk::Menu ();

    open->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), open));    
    for (int i=0; i<6; i++)
        rank[i]->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), rank[i]));    
    trash->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), trash));    
    untrash->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), untrash));    
    develop->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), develop));    
    rename->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), rename));    
    remove->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), remove));    
    selall->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), selall));    
    copyprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), copyprof));    
    pasteprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), pasteprof));    
    partpasteprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), partpasteprof));    
    applyprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), applyprof));    
    clearprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), clearprof));    
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

int addfl (void* data) {
    
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
    g_idle_add (addfl, ap);
}

void FileBrowser::addEntry_ (FileBrowserEntry* entry) {

    entry->selected = false;
    entry->drawable = false;
    entry->framed = editedFiles.find (entry->filename)!=editedFiles.end();
        
    // add button set to the thumbbrowserentry
    entry->addButtonSet (new FileThumbnailButtonSet (entry));
    entry->getThumbButtonSet()->setRank (entry->thumbnail->getRank());
    entry->getThumbButtonSet()->setInTrash (entry->thumbnail->getStage()==1);
    entry->getThumbButtonSet()->setButtonListener (this);
    entry->resize (options.thumbSize);

    // find place in abc order
    std::vector<ThumbBrowserEntryBase*>::iterator i = fd.begin();
    while (i!=fd.end() && *entry < *((FileBrowserEntry*)*i))
        i++;
        
    fd.insert (i, entry);    

    initEntry (entry);
    redraw ();
}

FileBrowserEntry* FileBrowser::delEntry (const Glib::ustring& fname) {

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

    for (int i=0; i<fd.size(); i++)
        delete fd[i];
    fd.clear ();
    selected.clear ();
    notifySelectionListener ();
    lastClicked = NULL;
    for (int i=0; i<fd.size(); i++)
        ((FileBrowserEntry*)fd[i])->thumbnail->decreaseRef ();
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
    if (m==open) {
        std::vector<Thumbnail*> entries;
        for (int i=0; i<mselected.size(); i++)
            entries.push_back (mselected[i]->thumbnail);
        tbl->openRequested (entries);
     }
    else if (m==remove)
        tbl->deleteRequested (mselected);
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
        selected.clear ();
        for (int i=0; i<fd.size(); i++)
            if (checkFilter (fd[i])) {
                fd[i]->selected = true;
                selected.push_back (fd[i]);
            }
        queue_draw ();
        notifySelectionListener ();
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

void FileBrowser::applyFilter (const BrowserFilter& filter) {

    this->filter = filter;

    // remove items not complying the filter from the selection
    bool selchanged = false;
    for (int i=0; i<fd.size(); i++) 
        if (fd[i]->selected && !checkFilter (fd[i])) {
            fd[i]->selected = false;
            std::vector<ThumbBrowserEntryBase*>::iterator j = std::find (selected.begin(), selected.end(), fd[i]);
            selected.erase (j);
            if (lastClicked==fd[i])
                lastClicked = NULL;
            selchanged = true;
        }
    if (selchanged)
        notifySelectionListener ();
    redraw ();
}

bool FileBrowser::checkFilter (ThumbBrowserEntryBase* entryb) { // true -> entry complies filter
    
    FileBrowserEntry* entry = (FileBrowserEntry*)entryb;
    // return false if basic filter settings are not satisfied
    if (filter.showRanked[entry->thumbnail->getRank()]==false || (entry->thumbnail->getStage()==1 && !filter.showTrash) || (entry->thumbnail->getStage()==0 && !filter.showNotTrash))
        return false;
    
    // check exif filter
    const CacheImageData* cfs = entry->thumbnail->getCacheImageData();
    double tol = 0.01;
    double tol2 = 1e-8;
    
	if (!filter.exifFilterEnabled)
		return true;
	
	if (!cfs->exifValid)
		return (!filter.exifFilter.filterCamera || filter.exifFilter.cameras.count(cfs->camera)>0) 
			&& (!filter.exifFilter.filterLens || filter.exifFilter.lenses.count(cfs->lens)>0);
		
    return 
         (!filter.exifFilter.filterShutter || (rtengine::ImageMetaData::shutterFromString(rtengine::ImageMetaData::shutterToString(cfs->shutter)) >= filter.exifFilter.shutterFrom-tol2 && rtengine::ImageMetaData::shutterFromString(rtengine::ImageMetaData::shutterToString(cfs->shutter)) <= filter.exifFilter.shutterTo+tol2))
      && (!filter.exifFilter.filterFNumber || (rtengine::ImageMetaData::apertureFromString(rtengine::ImageMetaData::apertureToString(cfs->fnumber)) >= filter.exifFilter.fnumberFrom-tol2 && rtengine::ImageMetaData::apertureFromString(rtengine::ImageMetaData::apertureToString(cfs->fnumber)) <= filter.exifFilter.fnumberTo+tol2))
      && (!filter.exifFilter.filterFocalLen || (cfs->focalLen >= filter.exifFilter.focalFrom-tol && cfs->focalLen <= filter.exifFilter.focalTo+tol))
	  && (!filter.exifFilter.filterISO     || (cfs->iso >= filter.exifFilter.isoFrom && cfs->iso <= filter.exifFilter.isoTo))
      && (!filter.exifFilter.filterCamera  || filter.exifFilter.cameras.count(cfs->camera)>0)
	  && (!filter.exifFilter.filterLens    || filter.exifFilter.lenses.count(cfs->lens)>0);
}

void FileBrowser::toTrashRequested (std::vector<FileBrowserEntry*> tbe) {

    for (int i=0; i<tbe.size(); i++) {
        if (tbe[i]->thumbnail->getStage()==1)
            continue;               
        tbe[i]->thumbnail->setStage (1);
        if (tbe[i]->getThumbButtonSet()) {
            tbe[i]->getThumbButtonSet()->setRank (tbe[i]->thumbnail->getRank());
            tbe[i]->getThumbButtonSet()->setInTrash (true);
        }
    }
    applyFilter (filter);
}

void FileBrowser::fromTrashRequested (std::vector<FileBrowserEntry*> tbe) {

    for (int i=0; i<tbe.size(); i++) {
        if (tbe[i]->thumbnail->getStage()==0)
            continue;
        tbe[i]->thumbnail->setStage (0);
        if (tbe[i]->getThumbButtonSet()) {
            tbe[i]->getThumbButtonSet()->setRank (tbe[i]->thumbnail->getRank());
            tbe[i]->getThumbButtonSet()->setInTrash (false);
        }
    }
    applyFilter (filter);
}

void FileBrowser::rankingRequested (std::vector<FileBrowserEntry*> tbe, int rank) {

    for (int i=0; i<tbe.size(); i++) {
        tbe[i]->thumbnail->setRank (rank);
        if (tbe[i]->getThumbButtonSet())
                tbe[i]->getThumbButtonSet()->setRank (tbe[i]->thumbnail->getRank());
    }
    applyFilter (filter);
}

void FileBrowser::buttonPressed (LWButton* button, int actionCode, void* actionData) {

    if (actionCode>=0 && actionCode<=5) { // rank
        std::vector<FileBrowserEntry*> tbe;
        tbe.push_back ((FileBrowserEntry*)actionData);
        rankingRequested (tbe, actionCode);
    }
    else if (actionCode==6 && tbl) { // to processin queue
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

void FileBrowser::redrawNeeded (ThumbBrowserEntryBase* entry) {

    if (entry->insideWindow (0, 0, internal.get_width(), internal.get_height())) {
        if (!internal.isDirty ()) {
            internal.setDirty ();
            internal.queue_draw ();
        }
    }        
}

void FileBrowser::redrawNeeded (LWButton* button) {
    
    queue_draw ();
}
