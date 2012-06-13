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
#include "filebrowser.h"
#include <glibmm.h>
#include "options.h"
#include "multilangmgr.h"
#include "clipboard.h"
#include "profilestore.h"
#include "procparamchangers.h"
#include "batchqueue.h"
#include "../rtengine/dfmanager.h"
#include "../rtengine/ffmanager.h"
#include "rtimage.h"

extern Options options;

FileBrowser::FileBrowser () 
    : tbl(NULL),numFiltered(0), partialPasteDlg(M("PARTIALPASTE_DIALOGLABEL")) {

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
    develop->set_image(*Gtk::manage(new RTImage ("processing.png")));
    pmenu->attach (*Gtk::manage(developfast = new Gtk::ImageMenuItem (M("FILEBROWSER_POPUPPROCESSFAST"))), 0, 1, p, p+1); p++;

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
    		submenuLabel->attach (*Gtk::manage(colorlabel[i] = new Gtk::ImageMenuItem (M("FILEBROWSER_POPUPCOLORLABEL")+":"+options.colorLabels[i])), 0, 1, p, p+1); p++;
    	}
    	submenuLabel->show_all ();
		menuLabel->set_submenu (*submenuLabel);
    }
    else{
    	for (int i=0; i<=5; i++){
    		pmenu->attach (*Gtk::manage(colorlabel[i] = new Gtk::ImageMenuItem (M("FILEBROWSER_POPUPCOLORLABEL")+":"+options.colorLabels[i])), 0, 1, p, p+1); p++;
    	}
    }
    for (int i=1; i<=5; i++){//set color label images
    	colorlabel[i]->set_image(*Gtk::manage(new RTImage (Glib::ustring::compose("%1%2%3","clabel",i,".png"))));
    }
        
    pmenu->attach (*Gtk::manage(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;

    /***********************
     * external programs
     * *********************/
#ifdef WIN32
    Gtk::manage(miOpenDefaultViewer=new Gtk::MenuItem (M("FILEBROWSER_OPENDEFAULTVIEWER")));
#else
    miOpenDefaultViewer=NULL;
#endif

    // Build a list of menu items
    mMenuExtProgs.clear(); amiExtProg=NULL;
    for (std::list<ExtProgAction*>::iterator it=extProgStore->lActions.begin();it!=extProgStore->lActions.end();it++) {
        ExtProgAction* pAct=*it;

        if (pAct->target==1 || pAct->target==2) mMenuExtProgs[pAct->GetFullName()]=pAct;
    }

    // Attach them to menu
    if (!mMenuExtProgs.empty() || miOpenDefaultViewer!=NULL) {
        amiExtProg=new Gtk::MenuItem*[mMenuExtProgs.size()];
        int itemNo=0;

        if (options.menuGroupExtProg) {
    	    pmenu->attach (*Gtk::manage(menuExtProg = new Gtk::MenuItem (M("FILEBROWSER_EXTPROGMENU"))), 0, 1, p, p+1); p++;
    	    Gtk::Menu* submenuExtProg = Gtk::manage (new Gtk::Menu());
            
            if (miOpenDefaultViewer!=NULL) {
                submenuExtProg->attach (*miOpenDefaultViewer, 0, 1, p, p+1); p++;
            }

            for (std::map<Glib::ustring, ExtProgAction*>::iterator it=mMenuExtProgs.begin();it!=mMenuExtProgs.end();it++,itemNo++) {
                submenuExtProg->attach (*Gtk::manage(amiExtProg[itemNo] = new Gtk::MenuItem ((*it).first)), 0, 1, p, p+1); p++;
            }

            submenuExtProg->show_all ();
		    menuExtProg->set_submenu (*submenuExtProg);
        } else {
            if (miOpenDefaultViewer!=NULL) {
                pmenu->attach (*miOpenDefaultViewer, 0, 1, p, p+1); p++;
            }

            for (std::map<Glib::ustring, ExtProgAction*>::iterator it=mMenuExtProgs.begin();it!=mMenuExtProgs.end();it++,itemNo++) {
                pmenu->attach (*Gtk::manage(amiExtProg[itemNo] = new Gtk::MenuItem ((*it).first)), 0, 1, p, p+1); p++;
            }
        }

        pmenu->attach (*Gtk::manage(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;
    }

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
    	menuProfileOperations->set_image(*Gtk::manage(new RTImage ("logoicon-wind.png")));

    	Gtk::Menu* submenuProfileOperations = Gtk::manage (new Gtk::Menu ());

    	submenuProfileOperations->attach (*Gtk::manage(copyprof = new Gtk::MenuItem (M("FILEBROWSER_COPYPROFILE"))), 0, 1, p, p+1); p++;
    	submenuProfileOperations->attach (*Gtk::manage(pasteprof = new Gtk::MenuItem (M("FILEBROWSER_PASTEPROFILE"))), 0, 1, p, p+1); p++;
    	submenuProfileOperations->attach (*Gtk::manage(partpasteprof = new Gtk::MenuItem (M("FILEBROWSER_PARTIALPASTEPROFILE"))), 0, 1, p, p+1); p++;
    	submenuProfileOperations->attach (*Gtk::manage(applyprof = new Gtk::MenuItem (M("FILEBROWSER_APPLYPROFILE"))), 0, 1, p, p+1); p++;
    	submenuProfileOperations->attach (*Gtk::manage(applypartprof = new Gtk::MenuItem (M("FILEBROWSER_APPLYPROFILE_PARTIAL"))), 0, 1, p, p+1); p++;
        submenuProfileOperations->attach (*Gtk::manage(execcustprof = new Gtk::MenuItem (M("FILEBROWSER_EXEC_CPB"))), 0, 1, p, p+1); p++;
    	submenuProfileOperations->attach (*Gtk::manage(clearprof = new Gtk::MenuItem (M("FILEBROWSER_CLEARPROFILE"))), 0, 1, p, p+1); p++;
    	submenuProfileOperations->show_all ();
    	menuProfileOperations->set_submenu (*submenuProfileOperations);

    	pmenu->attach (*Gtk::manage(menuIPTCOperations = new Gtk::ImageMenuItem (M("FILEBROWSER_POPUPIPTCOPERATIONS"))), 0, 1, p, p+1); p++;
    	menuIPTCOperations->set_image(*Gtk::manage(new Gtk::Image (argv0+"/images/barcode-16.png")));
    	Gtk::Menu* submenuIPTCOperations = Gtk::manage (new Gtk::Menu ());
    	submenuIPTCOperations->attach (*Gtk::manage(copyIPTC = new Gtk::MenuItem (M("FILEBROWSER_COPYIPTC"))), 0, 1, p, p+1); p++;
    	submenuIPTCOperations->attach (*Gtk::manage(pasteIPTC = new Gtk::MenuItem (M("FILEBROWSER_PASTEIPTC"))), 0, 1, p, p+1); p++;
    	submenuIPTCOperations->attach (*Gtk::manage(partpasteIPTC = new Gtk::MenuItem (M("FILEBROWSER_PARTIALPASTEIPTC"))), 0, 1, p, p+1); p++;
    	submenuIPTCOperations->attach (*Gtk::manage(resyncIPTC = new Gtk::MenuItem (M("FILEBROWSER_RESYNCIPTC"))), 0, 1, p, p+1); p++;
    	submenuIPTCOperations->show_all();
    	menuIPTCOperations->set_submenu( *submenuIPTCOperations );

    }
    else{
        pmenu->attach (*Gtk::manage(copyprof = new Gtk::MenuItem (M("FILEBROWSER_COPYPROFILE"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(pasteprof = new Gtk::MenuItem (M("FILEBROWSER_PASTEPROFILE"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(partpasteprof = new Gtk::MenuItem (M("FILEBROWSER_PARTIALPASTEPROFILE"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(applyprof = new Gtk::MenuItem (M("FILEBROWSER_APPLYPROFILE"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(applypartprof = new Gtk::MenuItem (M("FILEBROWSER_APPLYPROFILE_PARTIAL"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(execcustprof = new Gtk::MenuItem (M("FILEBROWSER_EXEC_CPB"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(clearprof = new Gtk::MenuItem (M("FILEBROWSER_CLEARPROFILE"))), 0, 1, p, p+1); p++;

        pmenu->attach (*Gtk::manage(copyIPTC = new Gtk::MenuItem (M("FILEBROWSER_COPYIPTC"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(pasteIPTC = new Gtk::MenuItem (M("FILEBROWSER_PASTEIPTC"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(partpasteIPTC = new Gtk::MenuItem (M("FILEBROWSER_PARTIALPASTEIPTC"))), 0, 1, p, p+1); p++;
        pmenu->attach (*Gtk::manage(resyncIPTC = new Gtk::MenuItem (M("FILEBROWSER_RESYNCIPTC"))), 0, 1, p, p+1); p++;

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

    for (int i=0; i<mMenuExtProgs.size(); i++)
    	amiExtProg[i]->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), amiExtProg[i]));

    if (miOpenDefaultViewer!=NULL) {
        miOpenDefaultViewer->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), miOpenDefaultViewer));
    }

    trash->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), trash));    
    untrash->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), untrash));    
    develop->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), develop));    
    developfast->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), developfast));
    rename->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), rename));    
    remove->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), remove));    
    removeInclProc->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), removeInclProc));
    selall->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), selall));
    copyTo->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), copyTo));
    moveTo->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), moveTo));
    copyprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), copyprof));    
    pasteprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), pasteprof));    
    partpasteprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), partpasteprof));
    copyIPTC->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), copyIPTC));
    pasteIPTC->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), pasteIPTC));
    partpasteIPTC->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), partpasteIPTC));
    resyncIPTC->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated),resyncIPTC ));
    applyprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), applyprof));    
    applypartprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), applypartprof));
    execcustprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), execcustprof));    
    clearprof->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), clearprof));
    cachemenu->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &FileBrowser::menuItemActivated), cachemenu));
}

FileBrowser::~FileBrowser ()
{
    delete pmenu;
    delete[] amiExtProg;
}

void FileBrowser::rightClicked (ThumbBrowserEntryBase* entry) {

    trash->set_sensitive (false);
    untrash->set_sensitive (false);
    for (size_t i=0; i<selected.size(); i++)
        if ((static_cast<FileBrowserEntry*>(selected[i]))->thumbnail->getRank()==-1) {
            untrash->set_sensitive (true);
            break;
        }
        for (size_t i=0; i<selected.size(); i++)
        if ((static_cast<FileBrowserEntry*>(selected[i]))->thumbnail->getRank()>=0) {
            trash->set_sensitive (true);
            break;
        }

    pasteprof->set_sensitive (clipboard.hasProcParams());
    partpasteprof->set_sensitive (clipboard.hasProcParams());
    copyprof->set_sensitive (selected.size()==1);
    clearprof->set_sensitive (!selected.empty());
    copyIPTC->set_sensitive (selected.size()==1);
    pasteIPTC->set_sensitive (clipboard.hasIPTC());
    partpasteIPTC->set_sensitive (clipboard.hasIPTC());
    resyncIPTC->set_sensitive (!selected.empty());

    // submenu applmenu
    int p = 0;
    Gtk::Menu* applmenu = Gtk::manage (new Gtk::Menu ());
    std::vector<Glib::ustring> profnames = profileStore.getProfileNames ();
    for (size_t i=0; i<profnames.size(); i++) {
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
    for (size_t i=0; i<profnames.size(); i++) {
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
        entries.push_back ((static_cast<FileBrowserEntry*>(entry))->thumbnail);
        tbl->openRequested (entries);
    }
}

struct addparams {
    FileBrowserIdleHelper* fbih;
    FileBrowserEntry* entry;
};

int AddEntryUIThread (void* data) {
    
    addparams* ap = static_cast<addparams*>(data);
    FileBrowserIdleHelper* fbih = ap->fbih;

    if (fbih->destroyed) {
        if (fbih->pending == 1)
            delete fbih;
        else    
            fbih->pending--;
        delete ap->entry;
        delete ap;

        return 0;
    }

    ap->fbih->fbrowser->addEntry_ (ap->entry);
    delete ap;
    fbih->pending--;

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
    entry->getThumbButtonSet()->setColorLabel ( options.getColorFromLabel( entry->thumbnail->getLabel() ));
    entry->getThumbButtonSet()->setInTrash (entry->thumbnail->getRank()==-1);
    entry->getThumbButtonSet()->setButtonListener (this);
    entry->resize (getCurrentThumbSize());

    // find place in abc order
	{
		// TODO: Check for Linux
		#ifdef WIN32
		Glib::RWLock::WriterLock l(entryRW);
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
	Glib::RWLock::WriterLock l(entryRW);
	#endif

    for (std::vector<ThumbBrowserEntryBase*>::iterator i=fd.begin(); i!=fd.end(); i++) 
        if ((*i)->filename==fname) {
            ThumbBrowserEntryBase* entry = *i;
            entry->selected = false;
            fd.erase (i);
            std::vector<ThumbBrowserEntryBase*>::iterator j = std::find (selected.begin(), selected.end(), entry);

            #ifdef WIN32
            l.release();
            #endif

            if (j!=selected.end()) {
                if (checkFilter (*j)) numFiltered--;
                selected.erase (j);
                notifySelectionListener ();
            }

            if (lastClicked==entry)
                lastClicked = NULL;
            redraw ();

            return (static_cast<FileBrowserEntry*>(entry));
        }
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
		Glib::RWLock::WriterLock l(entryRW);
		#endif
        
		selected.clear ();
		notifySelectionListener ();

        // The listener merges parameters with old values, so delete afterwards
		for (size_t i=0; i<fd.size(); i++)
        {
            delete fd[i];
        }
        fd.clear ();
	}

    lastClicked = NULL;
}

void FileBrowser::menuItemActivated (Gtk::MenuItem* m) {

    std::vector<FileBrowserEntry*> mselected;
    for (size_t i=0; i<selected.size(); i++)
        mselected.push_back (static_cast<FileBrowserEntry*>(selected[i]));

    if (!tbl || (m!=selall && mselected.empty()) )
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

    for (int j=0; j<mMenuExtProgs.size(); j++) {
        if (m==amiExtProg[j]) {
            ExtProgAction* pAct = mMenuExtProgs[m->get_label()];

            // Build vector of all file names
            std::vector<Glib::ustring> selFileNames;
            for (int i=0; i<selected.size(); i++) {
                Glib::ustring fn=selected[i]->thumbnail->getFileName();

                // Maybe batch processed version
                if (pAct->target==2) fn = Glib::ustring::compose ("%1.%2", BatchQueue::calcAutoFileNameBase(fn), options.saveFormatBatch.format);

                selFileNames.push_back(fn);
            }

            pAct->Execute(selFileNames);
            return;
        }
    }

    if (m==open) {
        std::vector<Thumbnail*> entries;
	for (size_t i=0; i<mselected.size(); i++)
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
        tbl->developRequested (mselected, false);
    else if (m==developfast)
        tbl->developRequested (mselected, true);

    else if (m==rename)
        tbl->renameRequested (mselected);
    else if (m==selall) {
        lastClicked = NULL;
		{
			// TODO: Check for Linux
			#ifdef WIN32
			Glib::RWLock::ReaderLock l(entryRW);
			#endif

            selected.clear ();
	    for (size_t i=0; i<fd.size(); i++)
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
	    for (size_t i=0; i<mselected.size(); i++){
			rtengine::procparams::ProcParams pp=mselected[i]->thumbnail->getProcParams();
			pp.raw.df_autoselect= true;
			pp.raw.dark_frame.clear();
			mselected[i]->thumbnail->setProcParams(pp,NULL,FILEBROWSER,false);
		}
    }else if (m==selectDF){
    	if( !mselected.empty() ){
    		rtengine::procparams::ProcParams pp=mselected[0]->thumbnail->getProcParams();
    		Gtk::FileChooserDialog fc("Dark Frame",Gtk::FILE_CHOOSER_ACTION_OPEN );
    		fc.add_button( Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    		fc.add_button( Gtk::StockID("gtk-apply"), Gtk::RESPONSE_APPLY);
    		if( pp.raw.dark_frame.empty())
    		   fc.set_current_folder( options.rtSettings.darkFramesPath );
    		else
    		   fc.set_filename( pp.raw.dark_frame );
    		if( fc.run() == Gtk::RESPONSE_APPLY ){
			for (size_t i=0; i<mselected.size(); i++){
    				rtengine::procparams::ProcParams pp=mselected[i]->thumbnail->getProcParams();
					pp.raw.dark_frame= fc.get_filename();
					pp.raw.df_autoselect= false;
					mselected[i]->thumbnail->setProcParams(pp,NULL,FILEBROWSER,false);
				}
			}
    	}
    }else if( m==thisIsDF){
    	if( !options.rtSettings.darkFramesPath.empty() && Gio::File::create_for_path(options.rtSettings.darkFramesPath)->query_exists() ){
		for (size_t i=0; i<mselected.size(); i++){
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
	    for (size_t i=0; i<mselected.size(); i++){
			rtengine::procparams::ProcParams pp=mselected[i]->thumbnail->getProcParams();
			pp.raw.ff_AutoSelect= true;
			pp.raw.ff_file.clear();
			mselected[i]->thumbnail->setProcParams(pp,NULL,FILEBROWSER,false);
		}
    }
    else if (m==selectFF){
    	if( !mselected.empty() ){
    		rtengine::procparams::ProcParams pp=mselected[0]->thumbnail->getProcParams();
    		Gtk::FileChooserDialog fc("Flat Field",Gtk::FILE_CHOOSER_ACTION_OPEN );
    		fc.add_button( Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    		fc.add_button( Gtk::StockID("gtk-apply"), Gtk::RESPONSE_APPLY);
    		if( pp.raw.ff_file.empty())
    		   fc.set_current_folder( options.rtSettings.flatFieldsPath );
    		else
    		   fc.set_filename( pp.raw.ff_file );
    		if( fc.run() == Gtk::RESPONSE_APPLY ){
			for (size_t i=0; i<mselected.size(); i++){
    				rtengine::procparams::ProcParams pp=mselected[i]->thumbnail->getProcParams();
					pp.raw.ff_file= fc.get_filename();
					pp.raw.ff_AutoSelect= false;
					mselected[i]->thumbnail->setProcParams(pp,NULL,FILEBROWSER,false);
				  }
			  }
    	}
    }
    else if( m==thisIsFF){
        if( !options.rtSettings.flatFieldsPath.empty() && Gio::File::create_for_path(options.rtSettings.flatFieldsPath)->query_exists() ){
            for (size_t i=0; i<mselected.size(); i++){
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
	    for (size_t i=0; i<mselected.size(); i++)
            mselected[i]->thumbnail->clearProcParams (FILEBROWSER);
        queue_draw ();
    }else if (m==copyIPTC)
		copyMetadata ();
	else if (m==pasteIPTC)
		pasteMetadata ();
	else if (m==partpasteIPTC)
		partPasteMetadata ();
	else if( m== resyncIPTC)
		resyncMetadata ();
    else if (m==execcustprof) {
	    for (size_t i=0; i<mselected.size(); i++)  {
            mselected[i]->thumbnail->createProcParamsForUpdate (false, true);

            // Empty run to update the thumb
            rtengine::procparams::ProcParams params = mselected[i]->thumbnail->getProcParams ();
            mselected[i]->thumbnail->setProcParams (params, NULL, FILEBROWSER);
    }
    } else if (m==clearFromCache) {
	    for (size_t i=0; i<mselected.size(); i++)
			tbl->clearFromCacheRequested (mselected, false);
		//queue_draw ();
    }
	else if (m==clearFromCacheFull) {
		for (size_t i=0; i<mselected.size(); i++)
			tbl->clearFromCacheRequested (mselected, true);
		//queue_draw ();
    } else if (miOpenDefaultViewer!=NULL && m==miOpenDefaultViewer) {
        openDefaultViewer(1);
    }
}

void FileBrowser::copyProfile () {

    if (selected.size()==1)
        clipboard.setProcParams ((static_cast<FileBrowserEntry*>(selected[0]))->thumbnail->getProcParams());
}

void FileBrowser::pasteProfile () {

    if (clipboard.hasProcParams()) {
        std::vector<FileBrowserEntry*> mselected;
        for (unsigned int i=0; i<selected.size(); i++)
            mselected.push_back (static_cast<FileBrowserEntry*>(selected[i]));

        if (!tbl || mselected.empty())
            return;

        for (unsigned int i=0; i<mselected.size(); i++) {
            // copying read only clipboard PartialProfile to a temporary one
            rtengine::procparams::PartialProfile cbPartProf = clipboard.getPartialProfile();
            rtengine::procparams::PartialProfile pastedPartProf(cbPartProf.pparams, cbPartProf.pedited, true);

            // applying the PartialProfile to the thumb's ProcParams
            mselected[i]->thumbnail->setProcParams (*pastedPartProf.pparams, pastedPartProf.pedited, FILEBROWSER);
            pastedPartProf.deleteInstance();
        }

        queue_draw ();
    }
}

void FileBrowser::partPasteProfile () {

    if (clipboard.hasProcParams()) {

        std::vector<FileBrowserEntry*> mselected;
        for (unsigned int i=0; i<selected.size(); i++)
            mselected.push_back (static_cast<FileBrowserEntry*>(selected[i]));

        if (!tbl || mselected.empty())
            return;

        int i = partialPasteDlg.run ();
        if (i == Gtk::RESPONSE_OK) {

            for (unsigned int i=0; i<mselected.size(); i++) {
                // copying read only clipboard PartialProfile to a temporary one, initialized to the thumb's ProcParams
                mselected[i]->thumbnail->createProcParamsForUpdate(false,false);  // this can execute customprofilebuilder to generate param file
                rtengine::procparams::PartialProfile cbPartProf = clipboard.getPartialProfile();
                rtengine::procparams::PartialProfile pastedPartProf(&mselected[i]->thumbnail->getProcParams (), NULL);

                // pushing the selected values of the clipboard PartialProfile to the temporary PartialProfile
                partialPasteDlg.applyPaste (pastedPartProf.pparams, pastedPartProf.pedited, cbPartProf.pparams, cbPartProf.pedited);

                // applying the temporary PartialProfile to the thumb's ProcParams
                mselected[i]->thumbnail->setProcParams (*pastedPartProf.pparams, pastedPartProf.pedited, FILEBROWSER);
                pastedPartProf.deleteInstance();
            }

            queue_draw ();
        }
        partialPasteDlg.hide ();
    }
}

void FileBrowser::copyMetadata () {

    if (selected.size()==1)
        clipboard.setIPTC( ((FileBrowserEntry*)selected[0])->thumbnail->getMetadata()->getIPTCData() );
}

void FileBrowser::pasteMetadata() {

    std::vector<FileBrowserEntry*> mselected;
    for (int i=0; i<selected.size(); i++)
        mselected.push_back ((FileBrowserEntry*)selected[i]);

    if (!tbl || mselected.size()==0)
        return;

    for (int i=0; i<mselected.size(); i++)
        mselected[i]->thumbnail->getMetadata()->setIPTCData( clipboard.getIPTC() );
}

void FileBrowser::partPasteMetadata () {

    std::vector<FileBrowserEntry*> mselected;
    for (int i=0; i<selected.size(); i++)
        mselected.push_back ((FileBrowserEntry*)selected[i]);

    if (!tbl || mselected.size()==0)
        return;

    PartialPasteIPTCDlg dlg( clipboard.getIPTC() );

    if ( dlg.run ()) {
    	rtengine::MetadataList iptc = dlg.getIPTC();
        for (int i=0; i<mselected.size(); i++) {
            mselected[i]->thumbnail->getMetadata()->setIPTCData( iptc );
        }
    }
}

void FileBrowser::resyncMetadata ()
{
    std::vector<FileBrowserEntry*> mselected;
    for (int i=0; i<selected.size(); i++)
        mselected.push_back ((FileBrowserEntry*)selected[i]);

    if (!tbl || mselected.size()==0)
        return;

    for (int i=0; i<mselected.size(); i++)
        mselected[i]->thumbnail->getMetadata()->resync();
}

void FileBrowser::openDefaultViewer (int destination) {
    bool success=true;
    if (selected.size()==1)
        success=(static_cast<FileBrowserEntry*>(selected[0]))->thumbnail->openDefaultViewer(destination);
    
    if (!success) {
        Gtk::MessageDialog msgd (M("MAIN_MSG_IMAGEUNPROCESSED"), true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
        msgd.run ();
    }
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
    else if (event->keyval==GDK_F2 && !(event->state & GDK_CONTROL_MASK)) {
        menuItemActivated (rename);
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

    rtengine::procparams::PartialProfile* partProfile = profileStore.getProfile (ppname);
    if (partProfile->pparams && !selected.empty()) {
	    for (size_t i=0; i<selected.size(); i++)
            (static_cast<FileBrowserEntry*>(selected[i]))->thumbnail->setProcParams (*partProfile->pparams, partProfile->pedited, FILEBROWSER);
        queue_draw ();
    }
}

void FileBrowser::applyPartialMenuItemActivated (Glib::ustring ppname) {

	if (!tbl || selected.empty())
		return;

	rtengine::procparams::PartialProfile* srcProfiles = profileStore.getProfile (ppname);

	if (srcProfiles->pparams) {
		if (partialPasteDlg.run()==Gtk::RESPONSE_OK) {

			for (size_t i=0; i<selected.size(); i++) {
				selected[i]->thumbnail->createProcParamsForUpdate(false, false);  // this can execute customprofilebuilder to generate param file

				rtengine::procparams::PartialProfile dstProfile(true);
				*dstProfile.pparams = (static_cast<FileBrowserEntry*>(selected[i]))->thumbnail->getProcParams ();
				dstProfile.set(true);
				partialPasteDlg.applyPaste (dstProfile.pparams, dstProfile.pedited, srcProfiles->pparams, srcProfiles->pedited);
				(static_cast<FileBrowserEntry*>(selected[i]))->thumbnail->setProcParams (*dstProfile.pparams, dstProfile.pedited, FILEBROWSER);
				dstProfile.deleteInstance();
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
		Glib::RWLock::ReaderLock l(entryRW);  // Don't make this a writer lock!
		#endif

		for (size_t i=0; i<fd.size(); i++) {
    	    if (checkFilter (fd[i]))
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
	}

    if (selchanged)
        notifySelectionListener ();
    redraw ();
}

bool FileBrowser::checkFilter (ThumbBrowserEntryBase* entryb) { // true -> entry complies filter
    
    FileBrowserEntry* entry = static_cast<FileBrowserEntry*>(entryb);
    // return false if basic filter settings are not satisfied
    if ((entry->thumbnail->getRank()>=0 && filter.showRanked[entry->thumbnail->getRank()]==false ) ||
        (filter.showCLabeled[options.getColorFromLabel( entry->thumbnail->getLabel() )]==false ) ||

        ((entry->thumbnail->hasProcParams() && filter.showEdited[0]) && !filter.showEdited[1]) ||
        ((!entry->thumbnail->hasProcParams() && filter.showEdited[1])&& !filter.showEdited[0]) ||

        ((entry->thumbnail->isRecentlySaved() && filter.showRecentlySaved[0]) && !filter.showRecentlySaved[1]) ||
        ((!entry->thumbnail->isRecentlySaved() && filter.showRecentlySaved[1]) && !filter.showRecentlySaved[0]) ||

        (entry->thumbnail->getRank()==-1 && !filter.showTrash) ||
        (entry->thumbnail->getRank()>=0 && !filter.showNotTrash))
        return false;

    // return false is query is not satisfied
    if (!filter.queryFileName.empty()){
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
			&& (!filter.exifFilter.filterFiletype || filter.exifFilter.filetypes.count(cfs->filetype)>0)
			&& (!filter.exifFilter.filterExpComp || filter.exifFilter.expcomp.count(cfs->expcomp)>0);

    return 
         (!filter.exifFilter.filterShutter || (rtengine::ImageMetaData::shutterFromString(rtengine::ImageMetaData::shutterToString(cfs->shutter)) >= filter.exifFilter.shutterFrom-tol2 && rtengine::ImageMetaData::shutterFromString(rtengine::ImageMetaData::shutterToString(cfs->shutter)) <= filter.exifFilter.shutterTo+tol2))
      && (!filter.exifFilter.filterFNumber || (rtengine::ImageMetaData::apertureFromString(rtengine::ImageMetaData::apertureToString(cfs->fnumber)) >= filter.exifFilter.fnumberFrom-tol2 && rtengine::ImageMetaData::apertureFromString(rtengine::ImageMetaData::apertureToString(cfs->fnumber)) <= filter.exifFilter.fnumberTo+tol2))
      && (!filter.exifFilter.filterFocalLen || (cfs->focalLen >= filter.exifFilter.focalFrom-tol && cfs->focalLen <= filter.exifFilter.focalTo+tol))
	  && (!filter.exifFilter.filterISO     || (cfs->iso >= filter.exifFilter.isoFrom && cfs->iso <= filter.exifFilter.isoTo))
	  && (!filter.exifFilter.filterExpComp || filter.exifFilter.expcomp.count(cfs->expcomp)>0)
      && (!filter.exifFilter.filterCamera  || filter.exifFilter.cameras.count(cfs->camera)>0)
	  && (!filter.exifFilter.filterLens    || filter.exifFilter.lenses.count(cfs->lens)>0)
	  && (!filter.exifFilter.filterFiletype  || filter.exifFilter.filetypes.count(cfs->filetype)>0);
}

void FileBrowser::toTrashRequested (std::vector<FileBrowserEntry*> tbe) {

	for (size_t i=0; i<tbe.size(); i++) {
    	// try to load the last saved parameters from the cache or from the paramfile file
    	tbe[i]->thumbnail->createProcParamsForUpdate(false, false);  // this can execute customprofilebuilder to generate param file

    	// no need to notify listeners as item goes to trash, likely to be deleted

    	if (tbe[i]->thumbnail->getRank()==-1)
            continue;               
        tbe[i]->thumbnail->setRank (-1);
        if (tbe[i]->getThumbButtonSet()) {
            tbe[i]->getThumbButtonSet()->setRank (tbe[i]->thumbnail->getRank());
            //tbe[i]->getThumbButtonSet()->setColorLabel (options.getColorFromLabel( tbe[i]->thumbnail->getLabel() ));
            tbe[i]->getThumbButtonSet()->setInTrash (true);
            //tbe[i]->thumbnail->updateCache (); // needed to save the colorlabel to disk in the procparam file(s) and the cache image data file
        }
    }
    trash_changed().emit();
    applyFilter (filter);
}

void FileBrowser::fromTrashRequested (std::vector<FileBrowserEntry*> tbe) {

	for (size_t i=0; i<tbe.size(); i++) {
    	// if thumbnail was marked inTrash=true then param file must be there, no need to run customprofilebuilder

        if (tbe[i]->thumbnail->getRank()>=0)
            continue;
        tbe[i]->thumbnail->setRank (0); // Out of trash: set unranked
        if (tbe[i]->getThumbButtonSet()) {
            tbe[i]->getThumbButtonSet()->setRank (tbe[i]->thumbnail->getRank());
            //tbe[i]->getThumbButtonSet()->setColorLabel (options.getColorFromLabel( tbe[i]->thumbnail->getLabel() ));
            tbe[i]->getThumbButtonSet()->setInTrash (false);
            //tbe[i]->thumbnail->updateCache (); // needed to save the colorlabel to disk in the procparam file(s) and the cache image data file
        }
    }
    trash_changed().emit();
    applyFilter (filter);
}

void FileBrowser::rankingRequested (std::vector<FileBrowserEntry*> tbe, int rank) {

	for (size_t i=0; i<tbe.size(); i++) {
    	// try to load the last saved parameters from the cache or from the paramfile file
    	tbe[i]->thumbnail->createProcParamsForUpdate(false, false);  // this can execute customprofilebuilder to generate param file

    	// notify listeners TODO: should do this ONLY when params changed by customprofilebuilder?
    	tbe[i]->thumbnail->notifylisterners_procParamsChanged(FILEBROWSER);

        tbe[i]->thumbnail->setRank (rank);

        if (tbe[i]->getThumbButtonSet()){
        	tbe[i]->getThumbButtonSet()->setRank (tbe[i]->thumbnail->getRank());
        	tbe[i]->getThumbButtonSet()->setInTrash (false);
        }
    }
    applyFilter (filter);
}

void FileBrowser::colorlabelRequested (std::vector<FileBrowserEntry*> tbe, int colorlabel) {

	for (size_t i=0; i<tbe.size(); i++) {
    	// try to load the last saved parameters from the cache or from the paramfile file
    	tbe[i]->thumbnail->createProcParamsForUpdate(false, false);  // this can execute customprofilebuilder to generate param file

    	// notify listeners TODO: should do this ONLY when params changed by customprofilebuilder?
    	tbe[i]->thumbnail->notifylisterners_procParamsChanged(FILEBROWSER);

    	if( colorlabel<=0 )
            tbe[i]->thumbnail->setLabel ( "" );
    	else if( colorlabel< options.colorLabels.size() )
    		tbe[i]->thumbnail->setLabel ( options.colorLabels[colorlabel] );

        if (tbe[i]->getThumbButtonSet())
                tbe[i]->getThumbButtonSet()->setColorLabel ( options.getColorFromLabel( tbe[i]->thumbnail->getLabel() ));
    }
    applyFilter (filter);
}

void FileBrowser::buttonPressed (LWButton* button, int actionCode, void* actionData) {

    if (actionCode>=0 && actionCode<=5) { // rank
        std::vector<FileBrowserEntry*> tbe;
	    tbe.push_back (static_cast<FileBrowserEntry*>(actionData));
        rankingRequested (tbe, actionCode);
    }
    else if (actionCode==6 && tbl) { // to processing queue
        std::vector<FileBrowserEntry*> tbe;
	    tbe.push_back (static_cast<FileBrowserEntry*>(actionData));
        tbl->developRequested (tbe, false); // not a fast, but a FULL mode
    }
    else if (actionCode==7) { // to trash / undelete
        std::vector<FileBrowserEntry*> tbe;
	    FileBrowserEntry* entry = static_cast<FileBrowserEntry*>(actionData);
        tbe.push_back (entry);
        if (entry->thumbnail->getRank()>=0)
            toTrashRequested (tbe);
        else
            fromTrashRequested (tbe);
    }
}

void FileBrowser::openNextImage () {
	// TODO: Check for Linux
	#ifdef WIN32
	Glib::RWLock::ReaderLock l(entryRW);
	#endif

    if (!fd.empty()) {
	    for (size_t i=fd.size()-1; i>0; i--)
            if (editedFiles.find (fd[i]->filename)!=editedFiles.end()) 
                if (i<fd.size()-1 && tbl) {
                    std::vector<Thumbnail*> entries;
                    entries.push_back ((static_cast<FileBrowserEntry*>(fd[i+1]))->thumbnail);
                    tbl->openRequested (entries);
                    return;
                }
        if (tbl) {
            std::vector<Thumbnail*> entries;
            entries.push_back ((static_cast<FileBrowserEntry*>(fd[0]))->thumbnail);
            tbl->openRequested (entries);
        }
    }
}

void FileBrowser::openPrevImage () {
	// TODO: Check for Linux
	#ifdef WIN32
	Glib::RWLock::ReaderLock l(entryRW);
	#endif

    if (!fd.empty()) {
	    for (size_t i=0; i<fd.size(); i++)
            if (editedFiles.find (fd[i]->filename)!=editedFiles.end()) 
                if (i>0 && tbl) {
                    std::vector<Thumbnail*> entries;
                    entries.push_back ((static_cast<FileBrowserEntry*>(fd[i-1]))->thumbnail);
                    tbl->openRequested (entries);
                    return;
                }
        if (tbl) {
            std::vector<Thumbnail*> entries;
            entries.push_back ((static_cast<FileBrowserEntry*>(fd[fd.size()-1]))->thumbnail);
            tbl->openRequested (entries);
        }
    }
}

int refreshThumbImagesUI (void* data) {
    (static_cast<FileBrowser*>(data))->_thumbRearrangementNeeded ();
    return 0;
}

void FileBrowser::_thumbRearrangementNeeded () {
    refreshThumbImages ();  // arrangeFiles is NOT enough
}

void FileBrowser::thumbRearrangementNeeded () {
    g_idle_add (refreshThumbImagesUI, this);
}

void FileBrowser::selectionChanged () {

    notifySelectionListener ();
}

void FileBrowser::notifySelectionListener () {

    if (tbl) {
        std::vector<Thumbnail*> thm;
        for (size_t i=0; i<selected.size(); i++)
            thm.push_back ((static_cast<FileBrowserEntry*>(selected[i]))->thumbnail);
        tbl->selectionChanged (thm);
    }    
}

void FileBrowser::redrawNeeded (LWButton* button) {
    
    queue_draw ();
}
FileBrowser::type_trash_changed FileBrowser::trash_changed () {
    return m_trash_changed;
}


// ExportPanel interface
void FileBrowser::exportRequested (){
	FileBrowser::menuItemActivated(developfast);
}

void FileBrowser::setExportPanel (ExportPanel* expanel) {

	exportPanel = expanel;
	exportPanel->set_sensitive (false);
	exportPanel->setExportPanelListener (this);
}
