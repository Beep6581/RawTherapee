/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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
#include <filecatalog.h>
#include <filepanel.h>
#include <options.h>
#include <cachemanager.h>
#include <multilangmgr.h>
#include <guiutils.h>
#include <glib/gstdio.h>
#include <iostream>
#include <iomanip>
#include <renamedlg.h>
#include <thumbimageupdater.h>
#include <safegtk.h>
#include <batchqueue.h>


#define CHECKTIME 2000

extern Glib::ustring argv0;

FileCatalog::FileCatalog (CoarsePanel* cp, ToolBar* tb, FilePanel* filepanel) :
		selectedDirectoryId(1),
		listener(NULL),
		fslistener(NULL),
		dirlistener(NULL),
		hasValidCurrentEFS(false),
		filterPanel(NULL),
		coarsePanel(cp),
		toolBar(tb),
		filepanel(filepanel),
		previewsToLoad(0),
		previewsLoaded(0)
		{

    inTabMode=false;

    //  construct and initialize thumbnail browsers
        fileBrowser = Gtk::manage( new FileBrowser() );
        fileBrowser->setFileBrowserListener (this);
        fileBrowser->setArrangement (ThumbBrowserBase::TB_Vertical);
        fileBrowser->show ();

    set_size_request(0,250);
    // construct trash panel with the extra "empty trash" button
    trashButtonBox = Gtk::manage( new Gtk::VBox );
    Gtk::Button* emptyT = Gtk::manage( new Gtk::Button (M("FILEBROWSER_EMPTYTRASH")));
    emptyT->set_tooltip_markup (M("FILEBROWSER_EMPTYTRASHHINT"));
    emptyT->set_image (*Gtk::manage(new Gtk::Image (Gtk::StockID("gtk-delete"), Gtk::ICON_SIZE_BUTTON)));
    emptyT->signal_pressed().connect (sigc::mem_fun(*this, &FileCatalog::emptyTrash));    
    trashButtonBox->pack_start (*emptyT, Gtk::PACK_SHRINK, 4);
    emptyT->show ();
    trashButtonBox->show ();
    
    //initialize hbToolBar1
    Gtk::HBox* hbToolBar1 = Gtk::manage(new Gtk::HBox ());

    //setup BrowsePath
    iRefreshWhite = new Gtk::Image(argv0+"/images/refresh_white.png");
    iRefreshRed = new Gtk::Image(argv0+"/images/refresh_red.png");

    BrowsePath = Gtk::manage(new Gtk::Entry ());
    BrowsePath->set_width_chars (50);
    BrowsePath->set_tooltip_markup (M("FILEBROWSER_BROWSEPATHHINT"));
    Gtk::HBox* hbBrowsePath = Gtk::manage(new Gtk::HBox ());
    buttonBrowsePath = Gtk::manage(new Gtk::Button ());
    buttonBrowsePath->set_image (*iRefreshWhite);
    buttonBrowsePath->set_tooltip_markup (M("FILEBROWSER_BROWSEPATHBUTTONHINT"));
    buttonBrowsePath->set_relief (Gtk::RELIEF_NONE);
    buttonBrowsePath->signal_clicked().connect( sigc::mem_fun(*this, &FileCatalog::buttonBrowsePathPressed) );
    hbBrowsePath->pack_start (*BrowsePath, Gtk::PACK_EXPAND_WIDGET,0);
    hbBrowsePath->pack_start (*buttonBrowsePath,Gtk::PACK_SHRINK, 0);
    hbToolBar1->pack_start (*hbBrowsePath, Gtk::PACK_EXPAND_WIDGET,0);

    BrowsePath->signal_activate().connect (sigc::mem_fun(*this, &FileCatalog::buttonBrowsePathPressed)); //respond to the Enter key

    //setup Query
    iQueryClear = new Gtk::Image(argv0+"/images/x_10.png");
    Gtk::Label* labelQuery = Gtk::manage(new Gtk::Label(M("FILEBROWSER_QUERYLABEL")));
    Query = Gtk::manage(new Gtk::Entry ()); // cannot use Gtk::manage here as FileCatalog::getFilter will fail on Query->get_text()
    Query->set_text("");
    Query->set_width_chars (20); // TODO !!! add this value to options?
    Query->set_tooltip_markup (M("FILEBROWSER_QUERYHINT"));
    Gtk::HBox* hbQuery = Gtk::manage(new Gtk::HBox ());
    buttonQueryClear = Gtk::manage(new Gtk::Button ());
    buttonQueryClear->set_image (*iQueryClear);
    buttonQueryClear->set_tooltip_markup (M("FILEBROWSER_QUERYBUTTONHINT"));
    buttonQueryClear->set_relief (Gtk::RELIEF_NONE);
    buttonQueryClear->signal_clicked().connect( sigc::mem_fun(*this, &FileCatalog::buttonQueryClearPressed) );
    hbQuery->pack_start (*labelQuery,Gtk::PACK_SHRINK, 0);
    hbQuery->pack_start (*Query,Gtk::PACK_SHRINK, 0);
    hbQuery->pack_start (*buttonQueryClear,Gtk::PACK_SHRINK, 0);
    hbToolBar1->pack_start (*hbQuery, Gtk::PACK_SHRINK,0);

    Query->signal_activate().connect (sigc::mem_fun(*this, &FileCatalog::executeQuery)); //respond to the Enter key

    // if NOT a single row toolbar
    if (!options.FileBrowserToolbarSingleRow) pack_start (*hbToolBar1, Gtk::PACK_SHRINK,0);

    // setup button bar
    buttonBar = Gtk::manage( new Gtk::HBox () );
    pack_start (*buttonBar, Gtk::PACK_SHRINK);
    
    buttonBar->pack_start (*Gtk::manage(new Gtk::VSeparator), Gtk::PACK_SHRINK);

    tbLeftPanel_1 = new Gtk::ToggleButton ();
    iLeftPanel_1_Show = new Gtk::Image(argv0+"/images/panel_to_right.png");
    iLeftPanel_1_Hide = new Gtk::Image(argv0+"/images/panel_to_left.png");

    tbLeftPanel_1->set_relief(Gtk::RELIEF_NONE);
    tbLeftPanel_1->set_active (true);
    tbLeftPanel_1->set_tooltip_markup (M("MAIN_TOOLTIP_SHOWHIDELP1"));
    tbLeftPanel_1->set_image (*iLeftPanel_1_Hide);
    tbLeftPanel_1->signal_toggled().connect( sigc::mem_fun(*this, &FileCatalog::tbLeftPanel_1_toggled) );
    buttonBar->pack_start (*tbLeftPanel_1, Gtk::PACK_SHRINK);

    buttonBar->pack_start (*(new Gtk::VSeparator), Gtk::PACK_SHRINK);


    iFilterClear = new Gtk::Image (argv0+"/images/filterclear_18.png");
    igFilterClear = new Gtk::Image (argv0+"/images/filterclearg_18.png");
    bFilterClear = Gtk::manage(new Gtk::ToggleButton ());
    bFilterClear->set_active (true);
    bFilterClear->set_image(*iFilterClear);//(*Gtk::manage(new Gtk::Image (argv0+"/images/filterclear.png")));
    bFilterClear->set_relief (Gtk::RELIEF_NONE);
    bFilterClear->set_tooltip_markup (M("FILEBROWSER_SHOWDIRHINT"));
    bFilterClear->signal_button_press_event().connect (sigc::mem_fun(*this, &FileCatalog::capture_event),false);
    bCateg[0] = bFilterClear->signal_toggled().connect (sigc::bind(sigc::mem_fun(*this, &FileCatalog::categoryButtonToggled), bFilterClear));
    buttonBar->pack_start (*bFilterClear, Gtk::PACK_SHRINK);
    buttonBar->pack_start (*Gtk::manage(new Gtk::VSeparator), Gtk::PACK_SHRINK);

    fltrVbox1 = Gtk::manage (new Gtk::VBox());
    fltrRankbox = Gtk::manage (new Gtk::HBox());
    fltrLabelbox = Gtk::manage (new Gtk::HBox());

    iUnRanked = new Gtk::Image (argv0+"/images/ratednot_10.png");
    igUnRanked = new Gtk::Image (argv0+"/images/ratednotg_10.png");
    bUnRanked = Gtk::manage( new Gtk::ToggleButton () );
    bUnRanked->set_active (false);
    bUnRanked->set_image (*igUnRanked);
    bUnRanked->set_relief (Gtk::RELIEF_NONE);
    bUnRanked->set_tooltip_markup (M("FILEBROWSER_SHOWUNRANKHINT"));
    bCateg[1] = bUnRanked->signal_toggled().connect (sigc::bind(sigc::mem_fun(*this, &FileCatalog::categoryButtonToggled), bUnRanked));    
    fltrRankbox->pack_start (*bUnRanked, Gtk::PACK_SHRINK);
    bUnRanked->signal_button_press_event().connect (sigc::mem_fun(*this, &FileCatalog::capture_event),false);

    for (int i=0; i<5; i++) {
        iranked[i] = new Gtk::Image (argv0+"/images/rated.png");
        igranked[i] = new Gtk::Image (argv0+"/images/grayrated.png");
        iranked[i]->show ();
        igranked[i]->show ();
        bRank[i] = Gtk::manage( new Gtk::ToggleButton () );
        bRank[i]->set_image (*igranked[i]);
        bRank[i]->set_relief (Gtk::RELIEF_NONE);
        fltrRankbox->pack_start (*bRank[i], Gtk::PACK_SHRINK);
        bCateg[i+2] = bRank[i]->signal_toggled().connect (sigc::bind(sigc::mem_fun(*this, &FileCatalog::categoryButtonToggled), bRank[i]));    
        bRank[i]->signal_button_press_event().connect (sigc::mem_fun(*this, &FileCatalog::capture_event),false);
    }  

    iUnCLabeled = new Gtk::Image (argv0+"/images/clabel0.png");
    igUnCLabeled = new Gtk::Image (argv0+"/images/cglabel0.png");
    bUnCLabeled = Gtk::manage(new Gtk::ToggleButton ());
    bUnCLabeled->set_active (false);
    bUnCLabeled->set_image (*igUnCLabeled);
    bUnCLabeled->set_relief (Gtk::RELIEF_NONE);
    bUnCLabeled->set_tooltip_markup (M("FILEBROWSER_SHOWUNCOLORHINT"));
    bCateg[7] = bUnCLabeled->signal_toggled().connect (sigc::bind(sigc::mem_fun(*this, &FileCatalog::categoryButtonToggled), bUnCLabeled));
    fltrLabelbox->pack_start (*bUnCLabeled, Gtk::PACK_SHRINK);
    bUnCLabeled->signal_button_press_event().connect (sigc::mem_fun(*this, &FileCatalog::capture_event),false);

    for (int i=0; i<5; i++) {
		iCLabeled[i] = new Gtk::Image (Glib::ustring::compose("%1%2%3%4",argv0,"/images/clabel",i+1,".png"));
		igCLabeled[i] = new Gtk::Image (Glib::ustring::compose("%1%2%3%4",argv0,"/images/cglabel",i+1,".png"));
		iCLabeled[i]->show ();
		igCLabeled[i]->show ();
		bCLabel[i] = Gtk::manage(new Gtk::ToggleButton ());
		bCLabel[i]->set_image (*igCLabeled[i]);
		bCLabel[i]->set_relief (Gtk::RELIEF_NONE);
		fltrLabelbox->pack_start (*bCLabel[i], Gtk::PACK_SHRINK);
		bCateg[i+8] = bCLabel[i]->signal_toggled().connect (sigc::bind(sigc::mem_fun(*this, &FileCatalog::categoryButtonToggled), bCLabel[i]));
		bCLabel[i]->signal_button_press_event().connect (sigc::mem_fun(*this, &FileCatalog::capture_event),false);
    }

    fltrVbox1->pack_start (*fltrRankbox, Gtk::PACK_SHRINK,0);
    fltrVbox1->pack_start (*fltrLabelbox, Gtk::PACK_SHRINK,0);
    buttonBar->pack_start (*fltrVbox1, Gtk::PACK_SHRINK);

    bRank[0]->set_tooltip_markup (M("FILEBROWSER_SHOWRANK1HINT"));
    bRank[1]->set_tooltip_markup (M("FILEBROWSER_SHOWRANK2HINT"));
    bRank[2]->set_tooltip_markup (M("FILEBROWSER_SHOWRANK3HINT"));
    bRank[3]->set_tooltip_markup (M("FILEBROWSER_SHOWRANK4HINT"));
    bRank[4]->set_tooltip_markup (M("FILEBROWSER_SHOWRANK5HINT"));
    
    bCLabel[0]->set_tooltip_markup (M("FILEBROWSER_SHOWCOLORLABEL1HINT"));
    bCLabel[1]->set_tooltip_markup (M("FILEBROWSER_SHOWCOLORLABEL2HINT"));
    bCLabel[2]->set_tooltip_markup (M("FILEBROWSER_SHOWCOLORLABEL3HINT"));
    bCLabel[3]->set_tooltip_markup (M("FILEBROWSER_SHOWCOLORLABEL4HINT"));
    bCLabel[4]->set_tooltip_markup (M("FILEBROWSER_SHOWCOLORLABEL5HINT"));

    buttonBar->pack_start (*Gtk::manage(new Gtk::VSeparator), Gtk::PACK_SHRINK);

    fltrVbox2 = Gtk::manage (new Gtk::VBox());
    fltrEditedBox = Gtk::manage (new Gtk::HBox());
    fltrRecentlySavedBox = Gtk::manage (new Gtk::HBox());

    // bEdited
    iEdited[0] = new Gtk::Image (argv0+"/images/editednot_10.png");
    igEdited[0] = new Gtk::Image (argv0+"/images/editednotg_10.png");
    iEdited[1] = new Gtk::Image (argv0+"/images/edited_10.png");
    igEdited[1] = new Gtk::Image (argv0+"/images/editedg_10.png");
    for (int i=0; i<2; i++) {
		iEdited[i]->show ();
		bEdited[i] = Gtk::manage(new Gtk::ToggleButton ());
		bEdited[i]->set_active (false);
		bEdited[i]->set_image (*igEdited[i]);
		bEdited[i]->set_relief (Gtk::RELIEF_NONE);
		fltrEditedBox->pack_start (*bEdited[i], Gtk::PACK_SHRINK);
		//13, 14
		bCateg[i+13] = bEdited[i]->signal_toggled().connect (sigc::bind(sigc::mem_fun(*this, &FileCatalog::categoryButtonToggled), bEdited[i]));
		bEdited[i]->signal_button_press_event().connect (sigc::mem_fun(*this, &FileCatalog::capture_event),false);
    }
    bEdited[0]->set_tooltip_markup (M("FILEBROWSER_SHOWEDITEDNOTHINT"));
    bEdited[1]->set_tooltip_markup (M("FILEBROWSER_SHOWEDITEDHINT"));

    // RecentlySaved
    iRecentlySaved[0] = new Gtk::Image (argv0+"/images/savednot_10.png");
    igRecentlySaved[0] = new Gtk::Image (argv0+"/images/savednotg_10.png");
    iRecentlySaved[1] = new Gtk::Image (argv0+"/images/saved_10.png");
    igRecentlySaved[1] = new Gtk::Image (argv0+"/images/savedg_10.png");
    for (int i=0; i<2; i++) {
    	iRecentlySaved[i]->show ();
		bRecentlySaved[i] = Gtk::manage(new Gtk::ToggleButton ());
		bRecentlySaved[i]->set_active (false);
		bRecentlySaved[i]->set_image (*igRecentlySaved[i]);
		bRecentlySaved[i]->set_relief (Gtk::RELIEF_NONE);
		fltrRecentlySavedBox->pack_start (*bRecentlySaved[i], Gtk::PACK_SHRINK);
		//15, 16
		bCateg[i+15] = bRecentlySaved[i]->signal_toggled().connect (sigc::bind(sigc::mem_fun(*this, &FileCatalog::categoryButtonToggled), bRecentlySaved[i]));
		bRecentlySaved[i]->signal_button_press_event().connect (sigc::mem_fun(*this, &FileCatalog::capture_event),false);
    }
    bRecentlySaved[0]->set_tooltip_markup (M("FILEBROWSER_SHOWRECENTLYSAVEDNOTHINT"));
    bRecentlySaved[1]->set_tooltip_markup (M("FILEBROWSER_SHOWRECENTLYSAVEDHINT"));

    fltrVbox2->pack_start (*fltrEditedBox, Gtk::PACK_SHRINK,0);
    fltrVbox2->pack_start (*fltrRecentlySavedBox, Gtk::PACK_SHRINK,0);
    buttonBar->pack_start (*fltrVbox2, Gtk::PACK_SHRINK);

    buttonBar->pack_start (*Gtk::manage(new Gtk::VSeparator), Gtk::PACK_SHRINK);

    // Trash
    iTrashEmpty = new Gtk::Image(argv0+"/images/trash-show-empty.png") ;
    iTrashFull  = new Gtk::Image(argv0+"/images/trash-show-full.png") ;

    bTrash = Gtk::manage(  new Gtk::ToggleButton () );
    bTrash->set_image (*iTrashEmpty);
    bTrash->set_relief (Gtk::RELIEF_NONE);
    bTrash->set_tooltip_markup (M("FILEBROWSER_SHOWTRASHHINT"));   
    bCateg[17] = bTrash->signal_toggled().connect (sigc::bind(sigc::mem_fun(*this, &FileCatalog::categoryButtonToggled), bTrash));
    bTrash->signal_button_press_event().connect (sigc::mem_fun(*this, &FileCatalog::capture_event),false);
    buttonBar->pack_start (*bTrash, Gtk::PACK_SHRINK);
    buttonBar->pack_start (*Gtk::manage(new Gtk::VSeparator), Gtk::PACK_SHRINK);
    fileBrowser->trash_changed().connect( sigc::mem_fun(*this, &FileCatalog::trashChanged) );

    // 0  - bFilterClear
    // 1  - bUnRanked
    // 2  - bRank[0]
    // 3  - bRank[1]
    // 4  - bRank[2]
    // 5  - bRank[3]
    // 6  - bRank[4]
    // 7  - bUnCLabeled
    // 8  - bCLabel[0]
    // 9  - bCLabel[1]
    // 10 - bCLabel[2]
    // 11 - bCLabel[3]
    // 12 - bCLabel[4]
    // 13 - bEdited[0]
    // 14 - bEdited[1]
    // 15 - bRecentlySaved[0]
    // 16 - bRecentlySaved[1]
    // 17 - bTrash

    categoryButtons[0] = bFilterClear;
    categoryButtons[1] = bUnRanked;
    for (int i=0; i<5; i++){ categoryButtons[i+2] = bRank[i];}
    categoryButtons[7] = bUnCLabeled;
    for (int i=0; i<5; i++){ categoryButtons[i+8] = bCLabel[i];}
    for (int i=0; i<2; i++){ categoryButtons[i+13] = bEdited[i];}
    for (int i=0; i<2; i++){ categoryButtons[i+15] = bRecentlySaved[i];}
    categoryButtons[17] = bTrash;

    exifInfo = Gtk::manage(new Gtk::ToggleButton ());
    exifInfo->set_image (*Gtk::manage(new Gtk::Image (argv0+"/images/info.png")));
    exifInfo->set_relief (Gtk::RELIEF_NONE);
    exifInfo->set_tooltip_markup (M("FILEBROWSER_SHOWEXIFINFO"));
    exifInfo->set_active( options.showFileNames );
    exifInfo->signal_toggled().connect(sigc::mem_fun(*this, &FileCatalog::exifInfoButtonToggled));
    buttonBar->pack_start (*exifInfo, Gtk::PACK_SHRINK);

    // thumbnail zoom
    Gtk::HBox* zoomBox = Gtk::manage( new Gtk::HBox () );
    zoomInButton  = Gtk::manage(  new Gtk::Button () );
    zoomInButton->set_image (*Gtk::manage(new Gtk::Image (Gtk::StockID("gtk-zoom-in"), Gtk::ICON_SIZE_SMALL_TOOLBAR)));
    zoomInButton->signal_pressed().connect (sigc::mem_fun(*this, &FileCatalog::zoomIn));
    zoomInButton->set_relief (Gtk::RELIEF_NONE);
    zoomInButton->set_tooltip_markup (M("FILEBROWSER_ZOOMINHINT"));
    zoomBox->pack_end (*zoomInButton, Gtk::PACK_SHRINK);
    zoomOutButton  = Gtk::manage( new Gtk::Button () );
    zoomOutButton->set_image (*Gtk::manage(new Gtk::Image (Gtk::StockID("gtk-zoom-out"), Gtk::ICON_SIZE_SMALL_TOOLBAR)));
    zoomOutButton->signal_pressed().connect (sigc::mem_fun(*this, &FileCatalog::zoomOut));
    zoomOutButton->set_relief (Gtk::RELIEF_NONE);
    zoomOutButton->set_tooltip_markup (M("FILEBROWSER_ZOOMOUTHINT"));
    zoomBox->pack_end (*zoomOutButton, Gtk::PACK_SHRINK);   

    buttonBar->pack_start (*zoomBox, Gtk::PACK_SHRINK);
    buttonBar->pack_start (*Gtk::manage(new Gtk::VSeparator), Gtk::PACK_SHRINK);

    //iRightArrow = new Gtk::Image(argv0+"/images/right.png");
    //iRightArrow_red = new Gtk::Image(argv0+"/images/right_red.png");

    // if it IS a single row toolbar
    if (options.FileBrowserToolbarSingleRow) buttonBar->pack_start (*hbToolBar1, Gtk::PACK_EXPAND_WIDGET,0);

    tbRightPanel_1 = new Gtk::ToggleButton ();
    iRightPanel_1_Show = new Gtk::Image(argv0+"/images/panel_to_left.png");
    iRightPanel_1_Hide = new Gtk::Image(argv0+"/images/panel_to_right.png");

    tbRightPanel_1->set_relief(Gtk::RELIEF_NONE);
    tbRightPanel_1->set_active (true);
    tbRightPanel_1->set_tooltip_markup (M("MAIN_TOOLTIP_SHOWHIDERP1"));
    tbRightPanel_1->set_image (*iRightPanel_1_Hide);
    tbRightPanel_1->signal_toggled().connect( sigc::mem_fun(*this, &FileCatalog::tbRightPanel_1_toggled) );
    buttonBar->pack_end (*tbRightPanel_1, Gtk::PACK_SHRINK);

    buttonBar->pack_end (*coarsePanel, Gtk::PACK_SHRINK);   
    buttonBar->pack_end (*Gtk::manage(new Gtk::VSeparator), Gtk::PACK_SHRINK, 4);
    buttonBar->pack_end (*toolBar, Gtk::PACK_SHRINK);   
    buttonBar->pack_end (*Gtk::manage(new Gtk::VSeparator), Gtk::PACK_SHRINK, 4);

    // add default panel
    hBox = Gtk::manage( new Gtk::HBox () );
    hBox->show ();
    hBox->pack_end (*fileBrowser);
    fileBrowser->applyFilter (getFilter()); // warning: can call this only after all objects used in getFilter (e.g. Query) are instantiated
    //printf("FileCatalog::FileCatalog  fileBrowser->applyFilter (getFilter())\n");
    pack_start (*hBox);

    enabled = true;  

    lastScrollPos = 0;
    for (int i=0; i<18; i++) {
        hScrollPos[i] = 0;
        vScrollPos[i] = 0;
    }
    
    selectedDirectory = "";
#ifdef WIN32
    wdMonitor = NULL;
#endif   
}

FileCatalog::~FileCatalog(){
    for (int i=0; i<5; i++) {
   	    delete iranked[i];
	    delete igranked[i];
		delete iCLabeled[i];
		delete igCLabeled[i];
	}
    for (int i=0; i<2; i++) {
		delete iEdited[i];
		delete igEdited[i];
		delete iRecentlySaved[i];
		delete igRecentlySaved[i];
    }
    delete iFilterClear;
    delete igFilterClear;
    delete iUnRanked;
    delete igUnRanked;
    delete iUnCLabeled;
    delete igUnCLabeled;
	delete iTrashEmpty;
	delete iTrashFull;
	delete iRefreshWhite;
	delete iRefreshRed;
	delete iQueryClear;
	delete iLeftPanel_1_Show;
	delete iLeftPanel_1_Hide;
	delete iRightPanel_1_Show;
	delete iRightPanel_1_Hide;
}

bool FileCatalog::capture_event(GdkEventButton* event){
    // need to record modifiers on the button press, because signal_toggled does not pass the event.
    modifierKey = event->state;
    return false;
}

void FileCatalog::exifInfoButtonToggled()
{
   options.showFileNames =  exifInfo->get_active();
   fileBrowser->refreshThumbImages ();
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

#ifndef WIN32
    if (dirMonitor)
        dirMonitor->cancel ();
#else
    if (wdMonitor) {
        delete wdMonitor;
        wdMonitor = NULL;
    }
#endif

	// ignore old requests
	++selectedDirectoryId;

    // terminate thumbnail preview loading
    previewLoader->removeAllJobs ();

    // terminate thumbnail updater
    thumbImageUpdater->removeAllJobs ();

    // remove entries
    selectedDirectory = "";
    fileBrowser->close ();  
	fileNameList.clear ();
	
    dirEFS.clear ();
    hasValidCurrentEFS = false;
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
		//printf("FileCatalog::dirSelected  selectedDirectory = %s\n",selectedDirectory.c_str());
		BrowsePath->set_text (selectedDirectory);
		buttonBrowsePath->set_image (*iRefreshWhite);
        fileNameList = getFileList ();

        for (unsigned int i=0; i<fileNameList.size(); i++) {
            Glib::RefPtr<Gio::File> f = Gio::File::create_for_path(fileNameList[i]);
            if (f->get_parse_name() != openfile) // if we opened a file at the beginning don't add it again
                checkAndAddFile (f);
        }

        _refreshProgressBar ();
        filepanel->loadingThumbs(M("PROGRESSBAR_LOADINGTHUMBS"),0);

#ifdef WIN32
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

void FileCatalog::enableTabMode(bool enable) {
    inTabMode = enable;

    fileBrowser->enableTabMode(inTabMode);

    redrawAll();
}

void FileCatalog::_refreshProgressBar () {
    // In tab mode, no progress bar at all
    // Also mention that this progress bar only measures the FIRST pass (quick thumbnails)
    // The second, usually longer pass is done multithreaded down in the single entries and is NOT measured by this
    if (!inTabMode) {
    	Gtk::Notebook *nb =(Gtk::Notebook *)(filepanel->get_parent());
    	Gtk::Box* hbb=NULL;
    	Gtk::Label *label=NULL;
    	if( options.mainNBVertical )
           hbb = Gtk::manage (new Gtk::VBox ());
    	else
    	   hbb = Gtk::manage (new Gtk::HBox ());
        if (!previewsToLoad ) {
            hbb->pack_start (*Gtk::manage (new Gtk::Image (Gtk::Stock::DIRECTORY, Gtk::ICON_SIZE_MENU)));
            label = Gtk::manage (new Gtk::Label (M("MAIN_FRAME_FILEBROWSER")+" ("+Glib::ustring::format(fileBrowser->getNumFiltered())+"/"+Glib::ustring::format(previewsLoaded)+")"));
        } else {
            hbb->pack_start (*Gtk::manage (new Gtk::Image (Gtk::Stock::FIND, Gtk::ICON_SIZE_MENU)));
            label = Gtk::manage (new Gtk::Label (M("MAIN_FRAME_FILEBROWSER")+" [" +Glib::ustring::format(std::fixed, std::setprecision(0), std::setw(3), (double)previewsLoaded / previewsToLoad*100 )+"%]" ));
        	filepanel->loadingThumbs("",(double)previewsLoaded / previewsToLoad);
        }
        if( options.mainNBVertical )
        	label->set_angle(90);
        hbb->pack_start (*label);
        hbb->set_spacing (2);
        hbb->set_tooltip_markup (M("MAIN_FRAME_FILEBROWSER_TOOLTIP"));
        hbb->show_all ();
        nb->set_tab_label(*filepanel,*hbb);

    }
}

int refreshpb (void* data) {

    gdk_threads_enter ();
    ((FileCatalog*)data)->_refreshProgressBar ();
    gdk_threads_leave ();
    return 0;
}

void FileCatalog::previewReady (int dir_id, FileBrowserEntry* fdn) {
    GThreadLock lock;
	_previewReady (dir_id,fdn);
}
void FileCatalog::_previewReady (int dir_id, FileBrowserEntry* fdn) {

	if ( dir_id != selectedDirectoryId )
	{
		return;
	}

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
        if (cfs->iso>0 && (int)cfs->iso < dirEFS.isoFrom)
            dirEFS.isoFrom = (int)cfs->iso;
        if (cfs->iso>0 && (int)cfs->iso > dirEFS.isoTo)
            dirEFS.isoTo = (int)cfs->iso;
        if (cfs->focalLen < dirEFS.focalFrom)
            dirEFS.focalFrom = cfs->focalLen;
        if (cfs->focalLen > dirEFS.focalTo)
            dirEFS.focalTo = cfs->focalLen;
    }
    dirEFS.filetypes.insert (cfs->filetype);
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

	if (filterPanel) {
		filterPanel->set_sensitive (true);
	    if ( !hasValidCurrentEFS ){
	        currentEFS = dirEFS;
		    filterPanel->setFilter ( dirEFS,true );
	    }else {
		    filterPanel->setFilter ( currentEFS,false );
	    }
	}
 	// restart anything that might have been loaded low quality
 	fileBrowser->refreshQuickThumbImages();
 	fileBrowser->applyFilter (getFilter());
    _refreshProgressBar();
    filepanel->loadingThumbs(M("PROGRESSBAR_READY"),0);
}

void FileCatalog::previewsFinished (int dir_id) {

 	GThreadLock lock;

	if ( dir_id != selectedDirectoryId )
	{
		return;
	}

    if (!hasValidCurrentEFS) 
        currentEFS = dirEFS;

    g_idle_add (prevfinished, this);
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

void FileCatalog::refreshHeight () {
    int newHeight=fileBrowser->getEffectiveHeight() + buttonBar->get_height();
    set_size_request(0, newHeight);
}

void FileCatalog::_openImage (std::vector<Thumbnail*> tmb) {

    if (enabled && listener!=NULL) {
    	bool continueToLoad=true;
        for (size_t i=0; i< tmb.size() && continueToLoad; i++) {
            if (editedFiles.find (tmb[i]->getFileName())==editedFiles.end()){
            	listener->fileSelected (tmb[i]);
                if( !options.tabbedUI )
                	continueToLoad = false;
            }
            tmb[i]->decreaseRef ();
        }
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
    for (size_t i=0; i<tmb.size(); i++)
        tmb[i]->increaseRef ();
    g_idle_add (fcopenimg, params);
}

void FileCatalog::deleteRequested  (std::vector<FileBrowserEntry*> tbe, bool inclBatchProcessed) {

    if (tbe.size()==0)
        return;

    Gtk::MessageDialog msd (M("FILEBROWSER_DELETEDLGLABEL"), false, Gtk::MESSAGE_QUESTION, Gtk::BUTTONS_YES_NO, true);
    msd.set_secondary_text(Glib::ustring::compose ( inclBatchProcessed ? M("FILEBROWSER_DELETEDLGMSGINCLPROC") : M("FILEBROWSER_DELETEDLGMSG"), tbe.size()));

    if (msd.run()==Gtk::RESPONSE_YES) {
        for (unsigned int i=0; i<tbe.size(); i++) {
            Glib::ustring fname = tbe[i]->filename;
            // remove from browser
            FileBrowserEntry* t = fileBrowser->delEntry (fname);
//            t->thumbnail->decreaseRef ();
            delete t;
            // remove from cache
            cacheMgr->deleteEntry (fname);
            // delete from file system
            safe_g_remove (fname);
            // delete paramfile if found
            safe_g_remove (Glib::ustring(fname+paramFileExtension));
            safe_g_remove (Glib::ustring(removeExtension(fname)+paramFileExtension));
            // delete .thm file
            safe_g_remove (Glib::ustring(removeExtension(fname)+".thm"));
            safe_g_remove (Glib::ustring(removeExtension(fname)+".THM"));

			if (inclBatchProcessed) {
			    Glib::ustring procfName = Glib::ustring::compose ("%1.%2", BatchQueue::calcAutoFileNameBase(fname), options.saveFormat.format);
				if (safe_file_test (procfName, Glib::FILE_TEST_EXISTS)) safe_g_remove (procfName);
			}
        }
        redrawAll ();    
    }
}


void FileCatalog::copyMoveRequested  (std::vector<FileBrowserEntry*> tbe, bool moveRequested) {

    if (tbe.size()==0)
        return;


    Glib::ustring fc_title;
    if (moveRequested) fc_title=M("FILEBROWSER_POPUPMOVETO");
    else fc_title=M("FILEBROWSER_POPUPCOPYTO");
    Gtk::FileChooserDialog fc(fc_title,Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER );
	fc.add_button( Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
	fc.add_button( Gtk::StockID("gtk-ok"), Gtk::RESPONSE_OK);
	// open dialog at the 1-st file's path
	fc.set_filename(tbe[0]->filename);
	//!!! TOFO prevent dialog closing on "enter" key press

	bool filecopymovecomplete;
	int i_copyindex;

	if( fc.run() == Gtk::RESPONSE_OK ){
		Glib::ustring dest_Dir = fc.get_current_folder();

		// iterate through selected files
		for (unsigned int i=0; i<tbe.size(); i++) {
			Glib::ustring src_fPath = tbe[i]->filename;
			Glib::ustring src_Dir = Glib::path_get_dirname(src_fPath);
			Glib::RefPtr<Gio::File> src_file = Gio::File::create_for_path ( src_fPath );
			if( !src_file ) continue; // if file is missing - skip it

			Glib::ustring fname=src_file->get_basename();
			Glib::ustring fname_noExt = removeExtension(fname);
			Glib::ustring fname_Ext = getExtension(fname);

			// construct  destination File Paths
			Glib::ustring dest_fPath = Glib::build_filename (dest_Dir, fname);
			Glib::ustring dest_fPath_param= dest_fPath + paramFileExtension;

			if (moveRequested && (src_Dir==dest_Dir)) continue;
			/* comparison of src_Dir and dest_Dir is done per image for compatibility with
			possible future use of Collections as source where each file's source path may be different.*/

			filecopymovecomplete = false;
			i_copyindex = 1;
			while(!filecopymovecomplete){
				// check for filename conflicts at destination - prevent overwriting (actually RT will crash on overwriting attempt)
				if (!safe_file_test(dest_fPath, Glib::FILE_TEST_EXISTS) && !safe_file_test(dest_fPath_param, Glib::FILE_TEST_EXISTS)){
					// copy/move file to destination
					Glib::RefPtr<Gio::File> dest_file = Gio::File::create_for_path ( dest_fPath );
					if (moveRequested) {
						// move file
						src_file->move(dest_file);
						// re-attach cache files
						cacheMgr->renameEntry (src_fPath, tbe[i]->thumbnail->getMD5(), dest_fPath);
						// remove from browser
						FileBrowserEntry* t = fileBrowser->delEntry (src_fPath);
					}
					else
						src_file->copy(dest_file);


					// attempt to copy/move paramFile only if it exist next to the src
					Glib::RefPtr<Gio::File> scr_param = Gio::File::create_for_path (  src_fPath + paramFileExtension );

					if (safe_file_test( src_fPath + paramFileExtension, Glib::FILE_TEST_EXISTS)){
						Glib::RefPtr<Gio::File> dest_param = Gio::File::create_for_path ( dest_fPath_param);
						// copy/move paramFile to destination
						if (moveRequested){
							if (safe_file_test( dest_fPath + paramFileExtension, Glib::FILE_TEST_EXISTS)){
								// profile already got copied to destination from cache after cacheMgr->renameEntry
								// delete source profile as cleanup
								safe_g_remove (src_fPath + paramFileExtension);
							}
							else
								scr_param->move(dest_param);
						}
						else
							scr_param->copy(dest_param);
					}
					filecopymovecomplete = true;
				}
				else{
					// adjust destination fname to avoid conflicts (append "_<index>", preserve extension)
					Glib::ustring dest_fname = Glib::ustring::compose("%1%2%3%4%5",fname_noExt,"_",i_copyindex,".",fname_Ext);
					// re-construct  destination File Paths
					dest_fPath = Glib::build_filename (dest_Dir, dest_fname);
					dest_fPath_param= dest_fPath + paramFileExtension;
					i_copyindex++;
				}
			}//while
		} // i<tbe.size() loop
		redrawAll ();
	} // Gtk::RESPONSE_OK
}
void FileCatalog::developRequested (std::vector<FileBrowserEntry*> tbe) {

    if (listener) {
    	std::vector<BatchQueueEntry*> entries;
        for (size_t i=0; i<tbe.size(); i++) {
            rtengine::procparams::ProcParams params = tbe[i]->thumbnail->getProcParams();
            rtengine::ProcessingJob* pjob = rtengine::ProcessingJob::create (tbe[i]->filename, tbe[i]->thumbnail->getType()==FT_Raw, params);
            double tmpscale;
            rtengine::IImage8* img = tbe[i]->thumbnail->processThumbImage (params, options.maxThumbnailHeight, tmpscale);
            if (img) {
                int pw = img->getWidth ();
                int ph = img->getHeight ();
                guint8* prev = new guint8 [pw*ph*3];
                memcpy (prev, img->getData (), pw*ph*3);
                img->free();
                entries.push_back(new BatchQueueEntry (pjob, params, tbe[i]->filename, prev, pw, ph, tbe[i]->thumbnail));
            }
            else {
                int pw, ph;
                tbe[i]->thumbnail->getThumbnailSize (pw, ph);
                entries.push_back(new BatchQueueEntry (pjob, params, tbe[i]->filename, NULL, pw, ph, tbe[i]->thumbnail));
            }
        }
        listener->addBatchQueueJobs( entries );
    }
}

void FileCatalog::renameRequested  (std::vector<FileBrowserEntry*> tbe) {

	bool success;

    RenameDialog* renameDlg = new RenameDialog ((Gtk::Window*)get_toplevel());

    for (size_t i=0; i<tbe.size(); i++) {
    	renameDlg->initName (Glib::path_get_basename (tbe[i]->filename), tbe[i]->thumbnail->getCacheImageData());

        Glib::ustring ofname = tbe[i]->filename;
        Glib::ustring dirName = Glib::path_get_dirname (tbe[i]->filename);
        Glib::ustring baseName = Glib::path_get_basename (tbe[i]->filename);

    	success = false;
        do {
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

				/* check if filename already exists*/
				if (safe_file_test (nfname, Glib::FILE_TEST_EXISTS)) {
					Glib::ustring msg_ = Glib::ustring("<b>") + nfname + ": " + M("MAIN_MSG_ALREADYEXISTS") + "</b>";
					Gtk::MessageDialog msgd (msg_, true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
					msgd.run ();
				}
				else {
					success = true;
            if (!safe_g_rename (ofname, nfname)) {
						cacheMgr->renameEntry (ofname, tbe[i]->thumbnail->getMD5(), nfname);
						reparseDirectory ();
					}
				}
			}
			else
				success = true;
        } while (!success);
		renameDlg->hide ();
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
            if (!safe_g_rename (ofname, nfname)) {
				cacheMgr->renameEntry (ofname, tbe[i]->thumbnail->getMD5(), nfname);
				// the remaining part (removing old and adding new entry) is done by the directory monitor
				reparseDirectory ();
//                on_dir_changed (Gio::File::create_for_path (nfname), Gio::File::create_for_path (nfname), Gio::FILE_MONITOR_EVENT_CHANGED, true);
            }
        }
    }
    */
}

void FileCatalog::clearFromCacheRequested  (std::vector<FileBrowserEntry*> tbe, bool leavenotrace) {

	 if (tbe.size()==0)
	        return;

	 for (unsigned int i=0; i<tbe.size(); i++) {
		 Glib::ustring fname = tbe[i]->filename;
		 // remove from cache
		 cacheMgr->clearFromCache (fname,leavenotrace);
	 }
}

void FileCatalog::categoryButtonToggled (Gtk::ToggleButton* b) {
    
    //was control key pressed
    bool control_down = modifierKey & GDK_CONTROL_MASK;

    //was shift key pressed
    bool shift_down   = modifierKey & GDK_SHIFT_MASK;

	for (int i=0; i<18; i++)
		bCateg[i].block (true);

	//button is already toggled when entering this function, so we switch it back to its initial state
    b->set_active(!b->get_active());

    //if both control and shift keys were pressed, do nothing
    if (!(control_down && shift_down)) {

		fileBrowser->getScrollPosition (hScrollPos[lastScrollPos], vScrollPos[lastScrollPos]);

		//we look how many stars are already toggled on, if any
		int toggled_stars_count=0, buttons=0, start_star=0, toggled_button=0;

		for (int i=0; i<18; i++) {
			if (categoryButtons[i]->get_active()) {
				if (i>0 && i<17) {
					toggled_stars_count ++;
					start_star = i;
				}
				buttons |= (1 << i);
			}
			if (categoryButtons[i] == b) toggled_button = i;
		}

		//if no modifier key were pressed, we can switch-on the button, and clear all others
		if (!(control_down || shift_down)) {
			for (int i=0; i<18; i++) {
				categoryButtons[i]->set_active (i==toggled_button);
			}
		}
		//modifier key allowed only for stars and color labels
		else if (toggled_button>0 && toggled_button<17) {
			if (control_down) {
				//control is pressed
				if (toggled_stars_count == 1 && (buttons & (1 << toggled_button))) {
					//we're deselecting the only star still active, so we activate the folder filter
					categoryButtons[0]->set_active(true);
					//and we deselect the toggled star
					categoryButtons[toggled_button]->set_active (false);
				}
				else if (toggled_stars_count >= 1) {
					//we toggle the state of a star (eventually another one than the only one selected)
					categoryButtons[toggled_button]->set_active(!categoryButtons[toggled_button]->get_active());
				}
				else {
					//no star selected
					//we deselect the 2 non star filters
					if (buttons &  1    ) categoryButtons[0]->set_active(false);
					if (buttons & (1 << 17)) categoryButtons[17]->set_active(false);
					//and we toggle on the star
					categoryButtons[toggled_button]->set_active (true);
				}
			}
			else {
				//shift is pressed, only allowed if 0 or 1 star & labels is selected
				if (!toggled_stars_count) {
					//we deselect the 2 non star filters
					if (buttons &  1      ) categoryButtons[0]->set_active(false);
					if (buttons & (1 << 7)) categoryButtons[7]->set_active(false);
					if (buttons & (1 << 13)) categoryButtons[13]->set_active(false);
					if (buttons & (1 << 17)) categoryButtons[17]->set_active(false);
					//and we set the start star to 1 (unrated images)
					start_star = 1;
					//we act as if one star were selected
					toggled_stars_count = 1;
				}
				if (toggled_stars_count == 1) {
					int current_star=MIN(start_star,toggled_button);
					int last_star   =MAX(start_star,toggled_button);
					//we permute the start and the end star for the next loop
					for (; current_star <= last_star; current_star++) {
						//we toggle on all the star in the range
						if (!(buttons & (1 << current_star))) categoryButtons[current_star]->set_active(true);
					}
				}
				//if more than one star & color label is selected, do nothing
			}
		}

		bool active_now, active_before;

		// FilterClear: set the right images
		// TODO: swapping FilterClear icon needs more work in categoryButtonToggled
		/*active_now = bFilterClear->get_active();
		active_before = buttons & (1 << (0)); // 0
		if      ( active_now && !active_before) bFilterClear->set_image (*iFilterClear);
		else if (!active_now &&  active_before) bFilterClear->set_image (*igFilterClear);*/

		// rank: set the right images
		for (int i=0; i<5; i++) {
			active_now = bRank[i]->get_active();
			active_before = buttons & (1 << (i+2)); // 2,3,4,5,6
			if      ( active_now && !active_before) bRank[i]->set_image (*iranked[i]);
			else if (!active_now &&  active_before) bRank[i]->set_image (*igranked[i]);
		}
		active_now = bUnRanked->get_active();
		active_before = buttons & (1 << (1)); // 1
		if      ( active_now && !active_before) bUnRanked->set_image (*iUnRanked);
		else if (!active_now &&  active_before) bUnRanked->set_image (*igUnRanked);

		// color labels: set the right images
		for (int i=0; i<5; i++) {
			active_now = bCLabel[i]->get_active();
			active_before = buttons & (1 << (i+8)); // 8,9,10,11,12
			if      ( active_now && !active_before) bCLabel[i]->set_image (*iCLabeled[i]);
			else if (!active_now &&  active_before) bCLabel[i]->set_image (*igCLabeled[i]);
		}
		active_now = bUnCLabeled->get_active();
		active_before = buttons & (1 << (7)); // 7
		if      ( active_now && !active_before) bUnCLabeled->set_image (*iUnCLabeled);
		else if (!active_now &&  active_before) bUnCLabeled->set_image (*igUnCLabeled);

		// Edited: set the right images
		for (int i=0; i<2; i++) {
			active_now = bEdited[i]->get_active();
			active_before = buttons & (1 << (i+13)); //13,14
			if      ( active_now && !active_before) bEdited[i]->set_image (*iEdited[i]);
			else if (!active_now &&  active_before) bEdited[i]->set_image (*igEdited[i]);
		}

		// RecentlySaved: set the right images
		for (int i=0; i<2; i++) {
			active_now = bRecentlySaved[i]->get_active();
			active_before = buttons & (1 << (i+15));//15,16
			if      ( active_now && !active_before) bRecentlySaved[i]->set_image (*iRecentlySaved[i]);
			else if (!active_now &&  active_before) bRecentlySaved[i]->set_image (*igRecentlySaved[i]);
		}

		fileBrowser->applyFilter (getFilter ());
		_refreshProgressBar();

		//rearrange panels according to the selected filter
		removeIfThere (hBox, trashButtonBox);
		if (bTrash->get_active ())
			hBox->pack_start (*trashButtonBox, Gtk::PACK_SHRINK, 4);
		hBox->queue_draw ();

		fileBrowser->setScrollPosition (hScrollPos[lastScrollPos], vScrollPos[lastScrollPos]);
    }

	for (int i=0; i<18; i++)
		bCateg[i].block (false);
}

BrowserFilter FileCatalog::getFilter () {

	BrowserFilter filter;

	bool anyRankFilterActive = bUnRanked->get_active () || bRank[0]->get_active () || bRank[1]->get_active () || bRank[2]->get_active () || bRank[3]->get_active () || bRank[4]->get_active ();
	bool anyCLabelFilterActive = bUnCLabeled->get_active () || bCLabel[0]->get_active ()|| bCLabel[1]->get_active ()|| bCLabel[2]->get_active ()|| bCLabel[3]->get_active ()|| bCLabel[4]->get_active ();
	bool anyEditedFilterActive = bEdited[0]->get_active() || bEdited[1]->get_active();
	bool anyRecentlySavedFilterActive = bRecentlySaved[0]->get_active() || bRecentlySaved[1]->get_active();
	/*
	 * filter is setup in 2 steps
	 * Step 1: handle individual filters
	*/
    filter.showRanked[0] = bFilterClear->get_active() || bUnRanked->get_active () || bTrash->get_active () || \
                           anyCLabelFilterActive || anyEditedFilterActive || anyRecentlySavedFilterActive;

    filter.showCLabeled[0] = bFilterClear->get_active() || bUnCLabeled->get_active () || bTrash->get_active ()  || \
        		             anyRankFilterActive || anyEditedFilterActive || anyRecentlySavedFilterActive;

    for (int i=1; i<=5; i++){
        filter.showRanked[i] = bFilterClear->get_active() || bRank[i-1]->get_active () || bTrash->get_active () || \
                               anyCLabelFilterActive || anyEditedFilterActive || anyRecentlySavedFilterActive;

        filter.showCLabeled[i] = bFilterClear->get_active() || bCLabel[i-1]->get_active () || bTrash->get_active ()  || \
                                 anyRankFilterActive || anyEditedFilterActive || anyRecentlySavedFilterActive;
    }

    for (int i=0; i<2; i++){
    	filter.showEdited[i] = bFilterClear->get_active() || bEdited[i]->get_active () || bTrash->get_active ()  || \
                               anyRankFilterActive || anyCLabelFilterActive || anyRecentlySavedFilterActive;

    	filter.showRecentlySaved[i] = bFilterClear->get_active() || bRecentlySaved[i]->get_active () || bTrash->get_active ()  || \
    	                          anyRankFilterActive || anyCLabelFilterActive || anyEditedFilterActive;
    }
    if( options.rtSettings.verbose ){
		printf ("\n**************** FileCatalog::getFilter *** AFTER STEP 1 \n");
		for (int i=0; i<=5; i++) printf ("filter.showRanked[%i] = %i\n",i,filter.showRanked[i]);
		for (int i=0; i<=5; i++) printf ("filter.showCLabeled[%i] = %i\n",i,filter.showCLabeled[i]);
		for (int i=0; i<2; i++) printf ("filter.showEdited[%i] = %i\n",i,filter.showEdited[i]);
		for (int i=0; i<2; i++)	printf ("filter.showRecentlySaved[%i] = %i\n",i,filter.showRecentlySaved[i]);
    }
	filter.multiselect = false;

	/*
	 * Step 2
	 * handle the case when more than 1 filter is selected. This overrides values set in Step
	 * if no filters in a group are active, filter.show for each member of that group will be set to true
	 * otherwise they are set based on UI input
	 */
    if ((anyRankFilterActive && anyCLabelFilterActive ) || \
    	(anyRankFilterActive && anyEditedFilterActive ) || \
    	(anyRankFilterActive && anyRecentlySavedFilterActive ) || \
    	(anyCLabelFilterActive && anyEditedFilterActive ) || \
    	(anyCLabelFilterActive && anyRecentlySavedFilterActive ) || \
        (anyEditedFilterActive && anyRecentlySavedFilterActive)){

    	filter.multiselect = true;
    	filter.showRanked[0] = anyRankFilterActive?bUnRanked->get_active ():true;
    	filter.showCLabeled[0] = anyCLabelFilterActive?bUnCLabeled->get_active ():true;

    	for (int i=1; i<=5; i++){
    		filter.showRanked[i] = anyRankFilterActive?bRank[i-1]->get_active ():true;
    		filter.showCLabeled[i] = anyCLabelFilterActive?bCLabel[i-1]->get_active ():true;
    	}
    	for (int i=0; i<2; i++){
    		filter.showEdited[i] = anyEditedFilterActive?bEdited[i]->get_active():true;
    	    filter.showRecentlySaved[i] = anyRecentlySavedFilterActive?bRecentlySaved[i]->get_active():true;
    	}
    	if( options.rtSettings.verbose ){
    		printf ("\n**************** FileCatalog::getFilter *** AFTER STEP 2 \n");
			for (int i=0; i<=5; i++) printf ("filter.showRanked[%i] = %i\n",i,filter.showRanked[i]);
			for (int i=0; i<=5; i++) printf ("filter.showCLabeled[%i] = %i\n",i,filter.showCLabeled[i]);
			for (int i=0; i<2; i++) printf ("filter.showEdited[%i] = %i\n",i,filter.showEdited[i]);
			for (int i=0; i<2; i++)	printf ("filter.showRecentlySaved[%i] = %i\n",i,filter.showRecentlySaved[i]);
			printf ("filter.multiselect = %i\n",filter.multiselect);
    	}
    }


    filter.showTrash = bFilterClear->get_active() || bTrash->get_active ();
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

    //TODO add support for more query options. e.g by date, iso, f-number, etc
    //TODO could use date:<value>;iso:<value>  etc
    // default will be filename

    /* // this is for safe execution if getFilter is called before Query object is instantiated
    Glib::ustring tempQuery;
    tempQuery="";
    if (Query) tempQuery = Query->get_text();
    */
    filter.queryString = Query->get_text(); // full query string from Query Entry
    filter.queryFileName = Query->get_text(); // for now Query is only by file name

    return filter;
}

void FileCatalog::filterChanged () {
    //TODO !!! there is too many repetitive and unnecessary executions of
	// " fileBrowser->applyFilter (getFilter()); " throughout the code
	// this needs further analysis and cleanup
    fileBrowser->applyFilter (getFilter());
    _refreshProgressBar();
}

int FileCatalog::reparseDirectory () {

    if (selectedDirectory=="")
        return 0;

    if (!safe_file_test (selectedDirectory, Glib::FILE_TEST_IS_DIR)) {
        closeDir ();
		return 0;
	}
	
	std::vector<Glib::ustring> nfileNameList = getFileList ();

	// check if a thumbnailed file has been deleted
	const std::vector<ThumbBrowserEntryBase*>& t = fileBrowser->getEntries ();
	std::vector<Glib::ustring> fileNamesToDel;
	for (size_t i=0; i<t.size(); i++)
		if (!safe_file_test (t[i]->filename, Glib::FILE_TEST_EXISTS))
			fileNamesToDel.push_back (t[i]->filename);
	for (size_t i=0; i<fileNamesToDel.size(); i++) {
		delete fileBrowser->delEntry (fileNamesToDel[i]);
		cacheMgr->deleteEntry (fileNamesToDel[i]);
	}

	// check if a new file has been added
	for (size_t i=0; i<nfileNameList.size(); i++) {
		bool found = false;
		for (int j=0; j<fileNameList.size(); j++)
			if (nfileNameList[i]==fileNameList[j]) {
				found = true;
				break;
			}
		if (!found) {
			checkAndAddFile (Gio::File::create_for_parse_name (nfileNameList[i]));
            _refreshProgressBar ();
		}
	}

	fileNameList = nfileNameList;
	return 1;
}

#ifdef WIN32
int winDirChangedUITread (void* cat) {
	gdk_threads_enter ();
    ((FileCatalog*)cat)->reparseDirectory ();
    gdk_threads_leave ();

    return 0;
}

void FileCatalog::winDirChanged () {
	g_idle_add(winDirChangedUITread, this);
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

    if (!file )
        return;
    if( !file->query_exists())
    	return;
    Glib::RefPtr<Gio::FileInfo> info = safe_query_file_info(file);
    if (info && info->get_file_type() != Gio::FILE_TYPE_DIRECTORY && (!info->is_hidden() || !options.fbShowHidden)) {
        int lastdot = info->get_name().find_last_of ('.');
        if (options.is_extention_enabled(lastdot!=(int)Glib::ustring::npos ? info->get_name().substr (lastdot+1) : "")){
						previewLoader->add (selectedDirectoryId,file->get_parse_name(),this);
            previewsToLoad++;
				}
    }
}

void FileCatalog::addAndOpenFile (const Glib::ustring& fname) {

    Glib::RefPtr<Gio::File> file = Gio::File::create_for_path (fname);
    if (!file )
        return;
    if( !file->query_exists())
    	return;
    Glib::RefPtr<Gio::FileInfo> info = safe_query_file_info(file);
    if( !info )
    	return;
    int lastdot = info->get_name().find_last_of ('.');
    if (options.is_extention_enabled(lastdot!=(int)Glib::ustring::npos ? info->get_name().substr (lastdot+1) : "")){
        // if supported, load thumbnail first
        Thumbnail* tmb = cacheMgr->getEntry (file->get_parse_name());
        if (tmb) {
            FileBrowserEntry* entry = new FileBrowserEntry (tmb, file->get_parse_name());
			_previewReady (selectedDirectoryId,entry);
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
    for (size_t i=0; i<t.size(); i++)
        if (((FileBrowserEntry*)t[i])->thumbnail->getStage()==1)
            toDel.push_back (((FileBrowserEntry*)t[i]));
    deleteRequested (toDel, false);
    trashChanged();
}

bool FileCatalog::trashIsEmpty () {
    const std::vector<ThumbBrowserEntryBase*> t = fileBrowser->getEntries ();
    for (size_t i=0; i<t.size(); i++)
        if (((FileBrowserEntry*)t[i])->thumbnail->getStage()==1)
            return false;

    return true;
}

void FileCatalog::zoomIn () {

    fileBrowser->zoomIn ();
    refreshHeight();
        
}
void FileCatalog::zoomOut () {

    fileBrowser->zoomOut ();
    refreshHeight();
        
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
    _refreshProgressBar();
}

void FileCatalog::setFilterPanel (FilterPanel* fpanel) { 

	filterPanel = fpanel; 
	filterPanel->set_sensitive (false);
	filterPanel->setFilterPanelListener (this);
}
void FileCatalog::trashChanged () {
    if (trashIsEmpty()) {
        bTrash->set_image(*iTrashEmpty);
    }
    else {
        bTrash->set_image(*iTrashFull);
    }
}
void FileCatalog::buttonQueryClearPressed () {
	Query->set_text("");
	FileCatalog::executeQuery ();
}

void FileCatalog::executeQuery(){
	// if BrowsePath text was changed, do a full browse;
	// otherwise filter only

	if (BrowsePath->get_text()!=selectedDirectory)
		buttonBrowsePathPressed ();
	else
		FileCatalog::filterChanged ();
}

void FileCatalog::buttonBrowsePathPressed () {
	Glib::ustring BrowsePathValue = BrowsePath->get_text();
	Glib::ustring DecodedPathPrefix="";
	Glib::ustring FirstChar;

	// handle shortcuts in the BrowsePath -- START
	// read the 1-st character from the path
	FirstChar = BrowsePathValue.substr (0,1);

	if (FirstChar=="~"){ // home directory
		DecodedPathPrefix = Glib::get_home_dir();
	}
	else if (FirstChar=="!"){ // user's pictures directory
		//DecodedPathPrefix = g_get_user_special_dir(G_USER_DIRECTORY_PICTURES);
		DecodedPathPrefix = Glib::get_user_special_dir (G_USER_DIRECTORY_PICTURES);
	}

	if (DecodedPathPrefix!=""){
		BrowsePathValue = Glib::ustring::compose ("%1%2",DecodedPathPrefix,BrowsePathValue.substr (1,BrowsePath->get_text_length()-1));
		BrowsePath->set_text(BrowsePathValue);
	}
	// handle shortcuts in the BrowsePath -- END

	// validate the path
	if (safe_file_test(BrowsePathValue, Glib::FILE_TEST_IS_DIR) && dirlistener){
		dirlistener->selectDir (BrowsePathValue);
	}
	else
		// error, likely path not found: show red arrow
		buttonBrowsePath->set_image (*iRefreshRed);
}

void FileCatalog::tbLeftPanel_1_visible (bool visible){
	if (visible)
		tbLeftPanel_1->show();
	else
		tbLeftPanel_1->hide();
}
void FileCatalog::tbRightPanel_1_visible (bool visible){
	if (visible)
		tbRightPanel_1->show();
	else
		tbRightPanel_1->hide();
}
void FileCatalog::tbLeftPanel_1_toggled () {
    removeIfThere (filepanel->dirpaned, filepanel->placespaned, false);
    if (tbLeftPanel_1->get_active()){
    	filepanel->dirpaned->pack1 (*filepanel->placespaned, false, true);
        tbLeftPanel_1->set_image (*iLeftPanel_1_Hide);
    }
    else {
    	tbLeftPanel_1->set_image (*iLeftPanel_1_Show);
    }
}

void FileCatalog::tbRightPanel_1_toggled () {
    if (tbRightPanel_1->get_active()){
    	filepanel->rightBox->show();
    	tbRightPanel_1->set_image (*iRightPanel_1_Hide);
    }
    else{
    	filepanel->rightBox->hide();
    	tbRightPanel_1->set_image (*iRightPanel_1_Show);
    }
}

bool FileCatalog::CheckSidePanelsVisibility(){
	if(tbLeftPanel_1->get_active()==false && tbRightPanel_1->get_active()==false)
		return false;
	else
		return true;
}
void FileCatalog::toggleSidePanels(){
	// toggle left AND right panels

	bool bAllSidePanelsVisible;
	bAllSidePanelsVisible= CheckSidePanelsVisibility();

	tbLeftPanel_1->set_active (!bAllSidePanelsVisible);
	tbRightPanel_1->set_active (!bAllSidePanelsVisible);
}

bool FileCatalog::handleShortcutKey (GdkEventKey* event) {

    bool ctrl = event->state & GDK_CONTROL_MASK;
    bool shift = event->state & GDK_SHIFT_MASK;
    bool alt = event->state & GDK_MOD1_MASK;
    
    modifierKey = event->state;
    
    // GUI Layout
    switch(event->keyval) {
        case GDK_l:
        	if (!alt)tbLeftPanel_1->set_active (!tbLeftPanel_1->get_active()); // toggle left panel
			if (alt && !ctrl) tbRightPanel_1->set_active (!tbRightPanel_1->get_active()); // toggle right panel
			if (alt && ctrl) {
				tbLeftPanel_1->set_active (!tbLeftPanel_1->get_active()); // toggle left panel
				tbRightPanel_1->set_active (!tbRightPanel_1->get_active()); // toggle right panel
			}
			return true;
        case GDK_m:
        	if (!ctrl && !alt) toggleSidePanels();
			return true;
    }

    if (!alt) {
        switch(event->keyval) {
			case GDK_grave:
				categoryButtonToggled(bUnRanked);
				return true;
            case GDK_1:
                categoryButtonToggled(bRank[0]);
                return true;
            case GDK_2:
                categoryButtonToggled(bRank[1]);
                return true;
            case GDK_3:
                categoryButtonToggled(bRank[2]);
                return true;
            case GDK_4:
                categoryButtonToggled(bRank[3]);
                return true;
            case GDK_5:
                categoryButtonToggled(bRank[4]);
                return true;
            case GDK_6:
                categoryButtonToggled(bEdited[0]);
                return true;
            case GDK_7:
                categoryButtonToggled(bEdited[1]);
                return true;

            case GDK_Return:
            case GDK_KP_Enter:
                FileCatalog::buttonBrowsePathPressed ();
                return true;
        }
    }
    
    if (alt) {
        switch(event->keyval) {
        	case GDK_grave:
				categoryButtonToggled(bUnCLabeled);
				return true;
            case GDK_1:
                categoryButtonToggled(bCLabel[0]);
                return true;
            case GDK_2:
                categoryButtonToggled(bCLabel[1]);
                return true;
            case GDK_3:
                categoryButtonToggled(bCLabel[2]);
                return true;
            case GDK_4:
                categoryButtonToggled(bCLabel[3]);
                return true;
            case GDK_5:
                categoryButtonToggled(bCLabel[4]);
                return true;
            case GDK_6:
                categoryButtonToggled(bRecentlySaved[0]);
                return true;
            case GDK_7:
                categoryButtonToggled(bRecentlySaved[1]);
                return true;
        }
    }
    
    if (!ctrl && !alt) {
        switch(event->keyval) {
			case GDK_d:
			case GDK_D:
				categoryButtonToggled(bFilterClear);
				return true;
			case GDK_t:
			case GDK_T:
				categoryButtonToggled(bTrash);
				return true;
        }
    }

    if (!ctrl) {
        switch(event->keyval) {

			case GDK_bracketright:
				coarsePanel->rotateRight();
				return true;
			case GDK_bracketleft:
				coarsePanel->rotateLeft();
				return true;
            case GDK_i:
            case GDK_I:
                exifInfo->set_active (!exifInfo->get_active());
                return true;

            case GDK_plus:
            case GDK_equal:
                zoomIn();
                return true;
            case GDK_minus:
            case GDK_underscore:
                zoomOut();
                return true;
        }
    }
    else { // with Ctrl
        switch (event->keyval) {
			case GDK_o:
				if (!alt){
					BrowsePath->select_region(0, BrowsePath->get_text_length());
					BrowsePath->grab_focus();
					return true;
				}
			case GDK_f:
				if (!alt){
					Query->select_region(0, Query->get_text_length());
					Query->grab_focus();
					return true;
				}
        }
    }

    return false;
}
