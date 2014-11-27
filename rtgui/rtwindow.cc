/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>, Oliver Duis <www.oliverduis.de>
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

#include "rtwindow.h"
#include "options.h"
#include "preferences.h"
#include "cursormanager.h"
#include "rtimage.h"
#include "whitebalance.h"

#if defined(__APPLE__)
static gboolean
osx_should_quit_cb (GtkosxApplication *app, gpointer data)
{
    RTWindow *rtWin = (RTWindow *)data;
    return rtWin->on_delete_event(0);
}

static void
osx_will_quit_cb (GtkosxApplication *app, gpointer data)
{
    RTWindow *rtWin = (RTWindow *)data;
    rtWin->on_delete_event(0);
    gtk_main_quit ();
}

bool RTWindow::osxFileOpenEvent(Glib::ustring path) {

    CacheManager* cm = CacheManager::getInstance();
    Thumbnail* thm= cm->getEntry( path );
    if(thm && fpanel){
        std::vector<Thumbnail*> entries;
	entries.push_back(thm);
        fpanel->fileCatalog->openRequested(entries);
        return true;
    }
    return false;
}

static gboolean
osx_open_file_cb (GtkosxApplication *app, gchar *path_, gpointer data)
{
    RTWindow *rtWin = (RTWindow *)data;
    if (!argv1.empty()) {
        // skip handling if we have a file argument or else we get double open of same file
        return false;
    }
    Glib::ustring path = Glib::ustring(path_);
    Glib::ustring suffix = path.length() > 4 ? path.substr(path.length()-3) : "";
    suffix = suffix.lowercase();
    if (suffix == "pp3")  {
        path = path.substr(0, path.length()-4);
    }
    return rtWin->osxFileOpenEvent(path);
}
#endif // __APPLE__

RTWindow::RTWindow ()
:mainNB(NULL)
,bpanel(NULL)
,splash(NULL)
,epanel(NULL)
,fpanel(NULL)
,btn_fullscreen(NULL)
{

    cacheMgr->init ();
    WhiteBalance::init();
    ProfilePanel::init();

    Glib::ustring fName = "rt-logo.png";
    Glib::ustring fullPath = RTImage::findIconAbsolutePath(fName);

#ifdef GLIBMM_EXCEPTIONS_ENABLED
		try { set_default_icon_from_file (fullPath);
		} catch(Glib::Exception& ex) {		printf ("%s\n", ex.what().c_str());	}
#else
 	  {		std::auto_ptr<Glib::Error> error;
				set_default_icon_from_file (fullPath, error);
		}
#endif //GLIBMM_EXCEPTIONS_ENABLED

#if defined(__APPLE__)
    {
        osxApp  = (GtkosxApplication *)g_object_new (GTKOSX_TYPE_APPLICATION, NULL);
        gboolean falseval = FALSE;
        gboolean trueval = TRUE;
        RTWindow *rtWin = this;
        g_signal_connect (osxApp, "NSApplicationBlockTermination", G_CALLBACK (osx_should_quit_cb), rtWin);
        g_signal_connect (osxApp, "NSApplicationWillTerminate",  G_CALLBACK (osx_will_quit_cb), rtWin);
        g_signal_connect (osxApp, "NSApplicationOpenFile", G_CALLBACK (osx_open_file_cb), rtWin);
        // RT don't have a menu, but we must create a dummy one to get the default OS X app menu working
        GtkWidget *menubar;
        menubar = gtk_menu_bar_new ();
        gtkosx_application_set_menu_bar (osxApp, GTK_MENU_SHELL (menubar));
        gtkosx_application_set_use_quartz_accelerators (osxApp, FALSE);
        gtkosx_application_ready (osxApp);
    }
#endif
    versionStr = "RawTherapee "+versionString;
    if (!versionSuffixString.empty())
        versionStr += " "+versionSuffixString;
    set_title_decorated("");
    property_allow_shrink() = true;
    set_default_size(options.windowWidth, options.windowHeight);
    set_modal(false);
    set_resizable(true);
    if (options.windowMaximized)
    	maximize();
    else {
    	unmaximize();
    	move(options.windowX,options.windowY);
    }

    on_delete_has_run = false;
    is_fullscreen = false;
    property_destroy_with_parent().set_value(false);
    signal_window_state_event().connect( sigc::mem_fun(*this, &RTWindow::on_window_state_event) );
	signal_key_press_event().connect( sigc::mem_fun(*this, &RTWindow::keyPressed) );

    if(simpleEditor)
    {
		epanel = Gtk::manage( new EditorPanel (NULL) );
		epanel->setParent (this);
    	add (*epanel);
    	show_all ();

        pldBridge=NULL; // No progress listener

    	CacheManager* cm = CacheManager::getInstance();
    	Thumbnail* thm= cm->getEntry( argv1 );
    	if(thm){
    		int error;
    		rtengine::InitialImage *ii= rtengine::InitialImage::load(argv1,thm->getType() == FT_Raw,&error,NULL);
    		epanel->open( thm, ii );
    	}
	} else {
		mainNB = Gtk::manage (new Gtk::Notebook ());
		mainNB->set_scrollable (true);
		mainNB->signal_switch_page().connect_notify( sigc::mem_fun(*this, &RTWindow::on_mainNB_switch_page) );

		fpanel =  new FilePanel () ;
		fpanel->setParent (this);

		// decorate tab
		if (options.mainNBVertical) {
			mainNB->set_tab_pos (Gtk::POS_LEFT);

			Gtk::VBox* vbf = Gtk::manage (new Gtk::VBox ());
			vbf->pack_start (*Gtk::manage (new Gtk::Image (Gtk::Stock::DIRECTORY, Gtk::ICON_SIZE_MENU)));
            Gtk::Label* l= Gtk::manage(new Gtk::Label (Glib::ustring(" ") + M("MAIN_FRAME_FILEBROWSER")));
			l->set_angle (90);
			vbf->pack_start (*l);
			vbf->set_spacing (2);
			vbf->set_tooltip_markup (M("MAIN_FRAME_FILEBROWSER_TOOLTIP"));
			vbf->show_all ();
			mainNB->append_page (*fpanel, *vbf);
		 } else {
			Gtk::HBox* hbf = Gtk::manage (new Gtk::HBox ());
			hbf->pack_start (*Gtk::manage (new Gtk::Image (Gtk::Stock::DIRECTORY, Gtk::ICON_SIZE_MENU)));
			hbf->pack_start (*Gtk::manage (new Gtk::Label (M("MAIN_FRAME_FILEBROWSER"))));
			hbf->set_spacing (2);
			hbf->set_tooltip_markup (M("MAIN_FRAME_FILEBROWSER_TOOLTIP"));
			hbf->show_all ();
			mainNB->append_page (*fpanel, *hbf);
		}

        bpanel = Gtk::manage ( new BatchQueuePanel () );
		bpanel->setParent (this);

		// decorate tab, the label is unimportant since its updated in batchqueuepanel anyway
        Gtk::Label* lbq = Gtk::manage ( new Gtk::Label (M("MAIN_FRAME_BATCHQUEUE")) );
        if (options.mainNBVertical) lbq->set_angle(90);
		mainNB->append_page (*bpanel, *lbq);

		// epanel is only for single tab mode
        epanel = Gtk::manage ( new EditorPanel (fpanel) );
		epanel->setParent (this);

		// decorate tab
		if (options.mainNBVertical) {
			Gtk::VBox* vbe = Gtk::manage (new Gtk::VBox ());
			vbe->pack_start (*Gtk::manage (new RTImage ("rt-logo.png")));
            Gtk::Label* l=Gtk::manage (new Gtk::Label( Glib::ustring(" ") + M("MAIN_FRAME_EDITOR") ));
			//l->set_markup(Glib::ustring("<b>Editor</b>"));  Bold difficult to read
			l->set_angle (90);
			vbe->pack_start (*l);
			vbe->set_spacing (2);
			vbe->set_tooltip_markup (M("MAIN_FRAME_EDITOR_TOOLTIP"));
			vbe->show_all ();
			epanel->tbTopPanel_1_visible(true); //show the toggle Top Panel button
			mainNB->append_page (*epanel, *vbe);
		} else {
			Gtk::HBox* hbe = Gtk::manage (new Gtk::HBox ());
			hbe->pack_start (*Gtk::manage (new RTImage ("rt-logo.png")));
			hbe->pack_start (*Gtk::manage (new Gtk::Label( Glib::ustring(" ") + M("MAIN_FRAME_EDITOR") )));
			hbe->set_spacing (2);
			hbe->set_tooltip_markup (M("MAIN_FRAME_EDITOR_TOOLTIP"));
			hbe->show_all ();
			epanel->tbTopPanel_1_visible(true); //show the toggle Top Panel button
			mainNB->append_page (*epanel, *hbe);
		}

		mainNB->set_current_page (mainNB->page_num (*fpanel));

		Gtk::VBox* mainBox = Gtk::manage (new Gtk::VBox ());
		mainBox->pack_start (*mainNB);

		// filling bottom box
		iFullscreen = new RTImage ("fullscreen.png");
		iFullscreen_exit = new RTImage ("fullscreen-exit.png");

		Gtk::LinkButton* rtWeb = Gtk::manage (new Gtk::LinkButton ("http://rawtherapee.com"));
		//Gtk::Button* preferences = Gtk::manage (new Gtk::Button (M("MAIN_BUTTON_PREFERENCES")+"..."));
		Gtk::Button* preferences = Gtk::manage (new Gtk::Button ());
		preferences->set_image (*Gtk::manage(new RTImage ("gtk-preferences.png")));
		preferences->set_tooltip_markup (M("MAIN_BUTTON_PREFERENCES"));
		preferences->signal_clicked().connect( sigc::mem_fun(*this, &RTWindow::showPreferences) );

		//btn_fullscreen = Gtk::manage( new Gtk::Button(M("MAIN_BUTTON_FULLSCREEN")));
		btn_fullscreen = Gtk::manage( new Gtk::Button());
		btn_fullscreen->set_tooltip_markup (M("MAIN_BUTTON_FULLSCREEN"));
		btn_fullscreen->set_image (*iFullscreen);
		btn_fullscreen->signal_clicked().connect( sigc::mem_fun(*this, &RTWindow::toggle_fullscreen) );
#if GTKMM_MINOR_VERSION >= 20
		if (options.mainNBVertical) {
			Gtk::VBox* bottomVBox = Gtk::manage (new Gtk::VBox ());
			bottomVBox->pack_start (prProgBar, Gtk::PACK_SHRINK, 1);
			bottomVBox->pack_end (*preferences, Gtk::PACK_SHRINK, 0);
			bottomVBox->pack_end (*btn_fullscreen, Gtk::PACK_EXPAND_WIDGET, 1);
			prProgBar.set_orientation(Gtk::PROGRESS_BOTTOM_TO_TOP);
			mainNB->set_action_widget( bottomVBox,Gtk::PACK_END);
			bottomVBox->show_all();
		}else{
			Gtk::HBox* bottomHBox = Gtk::manage (new Gtk::HBox ());
			bottomHBox->pack_end (*btn_fullscreen, Gtk::PACK_SHRINK, 1);
			bottomHBox->pack_end (*preferences, Gtk::PACK_SHRINK, 0);
			bottomHBox->pack_start (prProgBar, Gtk::PACK_EXPAND_WIDGET, 1);
			mainNB->set_action_widget( bottomHBox,Gtk::PACK_END);
			bottomHBox->show_all();
		}
#else
		Gtk::HBox* bottomBox = Gtk::manage (new Gtk::HBox ());
		bottomBox->pack_end (*btn_fullscreen, Gtk::PACK_SHRINK, 1);
		bottomBox->pack_end (*preferences, Gtk::PACK_SHRINK, 0);
		bottomBox->pack_start (prProgBar, Gtk::PACK_EXPAND_WIDGET, 1);
		mainBox->pack_start (*bottomBox, Gtk::PACK_SHRINK, 1);
#endif

		pldBridge = new PLDBridge(static_cast<rtengine::ProgressListener*>(this));

		Glib::RefPtr<Gtk::RcStyle> style = Gtk::RcStyle::create ();
		style->set_xthickness (0);
		style->set_ythickness (0);
		rtWeb->modify_style (style);

		add (*mainBox);
		show_all ();
    }
    if (!isSingleTabMode()&& !simpleEditor) epanel->hide_all();
}

RTWindow::~RTWindow()
{
    if(!simpleEditor)
        delete pldBridge;
    pldBridge = NULL;
#if defined(__APPLE__)
    g_object_unref (osxApp);
#endif
}

void RTWindow::findVerNumbers(int* numbers, Glib::ustring versionStr) {
	numbers[0] = numbers[1] = numbers[2] = numbers[3] = 0;
	int n=0;
	for (unsigned int i=0; i<versionStr.length(); i++) {
		char chr = (char)versionStr.at(i);
		if (chr >= '0' && chr <= '9') {
			numbers[n] *= 10;
			numbers[n] += (int)(chr - '0');
		}
		else {
			n++;
			if (n>4) {
				printf("Error: malformed version string; \"%s\" must follow this format: xx.xx.xx.xx. Admitting it's a developer version...\n", versionStr.c_str());
				// Reseting the already found numbers
				numbers[0] = numbers[1] = numbers[2] = numbers[3] = 100;
				return;
			}
		}
	}
}

void RTWindow::on_realize () {
    Gtk::Window::on_realize ();

    if( fpanel )
       fpanel->setAspect();

    cursorManager.init (get_window());

    // Check if first run of this version, then display the Release Notes text
    if (options.version != versionString) {
        int prevVerNbr[4];
        int currVerNbr[4];
        findVerNumbers(prevVerNbr, options.version);
        findVerNumbers(currVerNbr, versionString);

        // Now we can update the version parameter with the right value
        options.version = versionString;

        bool showReleaseNotes = false;
        // Check if the current version is newer
        if      (currVerNbr[0] > prevVerNbr[0]) showReleaseNotes = true;
        else if (currVerNbr[1] > prevVerNbr[1]) showReleaseNotes = true;
        else if (currVerNbr[2] > prevVerNbr[2]) showReleaseNotes = true;

        if (showReleaseNotes) {
        	// this is a first run!
            splash = new Splash (*this);
            splash->set_transient_for (*this);
            splash->signal_delete_event().connect( sigc::mem_fun(*this, &RTWindow::splashClosed) );
            if (splash->hasReleaseNotes()) {
            	splash->showReleaseNotes();
            	splash->show ();
            }
            else {
            	delete splash;
            	splash = NULL;
            }
        }
    }
}

bool RTWindow::on_window_state_event(GdkEventWindowState* event) {
	if (!event->new_window_state) {
		// Window mode
		options.windowMaximized = false;
	}
	else if (event->new_window_state & (GDK_WINDOW_STATE_MAXIMIZED|GDK_WINDOW_STATE_FULLSCREEN)) {
		// Fullscreen mode
		options.windowMaximized = true;
	}
	return true;
}

void RTWindow::on_mainNB_switch_page(GtkNotebookPage* page, guint page_num) {
	if(!on_delete_has_run) {
		if(isEditorPanel(page_num)) {
			if (isSingleTabMode() && epanel)
				MoveFileBrowserToEditor();
			
			EditorPanel *ep = static_cast<EditorPanel*>(mainNB->get_nth_page(page_num));
			ep->setAspect();
			
			if (!isSingleTabMode()){
				if (filesEdited.size()>0){
					set_title_decorated(ep->getFileName());
				}
			}
		} else {
			// in single tab mode with command line filename epanel does not exist yet
			if (isSingleTabMode() && epanel) {
				// Save profile on leaving the editor panel
				epanel->saveProfile();

				// Moving the FileBrowser only if the user has switched to the FileBrowser tab
				if (mainNB->get_nth_page(page_num)==fpanel)
					MoveFileBrowserToMain();
			}
		}
	}
}

void RTWindow::addEditorPanel (EditorPanel* ep, const std::string &name) {
    if (options.multiDisplayMode>0) {
        EditWindow * wndEdit = EditWindow::getInstance(this);
        wndEdit->show_all();
        wndEdit->addEditorPanel(ep,name);
    } else {
        ep->setParent (this);

        // construct closeable tab for the image
        Gtk::HBox* hb = Gtk::manage (new Gtk::HBox ());
        hb->pack_start (*Gtk::manage (new RTImage ("rtwindow.png")));
        hb->pack_start (*Gtk::manage (new Gtk::Label (Glib::path_get_basename (name))));
        hb->set_tooltip_markup (name);
        Gtk::Button* closeb = Gtk::manage (new Gtk::Button ());
        closeb->set_image (*Gtk::manage(new RTImage ("gtk-close.png")));
        closeb->set_relief (Gtk::RELIEF_NONE);
        closeb->set_focus_on_click (false);
        // make the button as small as possible
        Glib::RefPtr<Gtk::RcStyle> style = Gtk::RcStyle::create ();
        style->set_xthickness (0);
        style->set_ythickness (0);    

        closeb->modify_style (style);
        closeb->signal_clicked().connect( sigc::bind (sigc::mem_fun(*this, &RTWindow::remEditorPanel) , ep));
        hb->pack_end (*closeb);
        hb->set_spacing (2);
        hb->show_all ();

        mainNB->append_page (*ep, *hb);
        //ep->setAspect ();
        mainNB->set_current_page (mainNB->page_num (*ep));
        mainNB->set_tab_reorderable (*ep, true);

        set_title_decorated(name);
        epanels[ name ] = ep;
        filesEdited.insert ( name );
        fpanel->refreshEditedState (filesEdited);
        ep->tbTopPanel_1_visible(false); //hide the toggle Top Panel button
    }
}

void RTWindow::remEditorPanel (EditorPanel* ep) {
    if (ep->getIsProcessing()) return;  // Will crash if destroyed while loading

    if (options.multiDisplayMode>0) {
        EditWindow * wndEdit = EditWindow::getInstance(this);
        wndEdit->remEditorPanel(ep);
    } else {
    	bool queueHadFocus = (mainNB->get_current_page() == mainNB->page_num (*bpanel));
	    epanels.erase (ep->getFileName());
	    filesEdited.erase (ep->getFileName ());
	    fpanel->refreshEditedState (filesEdited);

	    mainNB->remove_page (*ep);

	    if (!isEditorPanel(mainNB->get_current_page())){
			if(!queueHadFocus)
				mainNB->set_current_page (mainNB->page_num (*fpanel));
		    set_title_decorated("");
	    } else {
		    EditorPanel* ep = static_cast<EditorPanel*>(mainNB->get_nth_page (mainNB->get_current_page()));
	        set_title_decorated(ep->getFileName());
	    }
        // TODO: ask what to do: close & apply, close & apply selection, close & revert, cancel
    }
}

bool RTWindow::selectEditorPanel(const std::string &name) {
    if (options.multiDisplayMode>0) {
        EditWindow * wndEdit = EditWindow::getInstance(this);
        if (wndEdit->selectEditorPanel(name)) {
        	set_title_decorated(name);
        	return true;
        }
    } else {
        std::map<Glib::ustring, EditorPanel*>::iterator iep = epanels.find(name);

        if (iep!=epanels.end()) {
            mainNB->set_current_page (mainNB->page_num (*iep->second));
            set_title_decorated(name);
            return true;
        }
        else{
        	//set_title_decorated(name);
        	//printf("RTWindow::selectEditorPanel - plain set\n");
        }
    }

    return false;
}

bool RTWindow::keyPressed (GdkEventKey* event) {

	bool ctrl = event->state & GDK_CONTROL_MASK;
	//bool shift = event->state & GDK_SHIFT_MASK;

        bool try_quit = false;
#if defined(__APPLE__)
        bool apple_cmd = event->state & GDK_MOD2_MASK;
        if (event->keyval == GDK_q && apple_cmd)
            try_quit = true;
#else
        if (event->keyval == GDK_q && ctrl)
            try_quit = true;
#endif
        if (try_quit) {
            if (!on_delete_event(0)) {
                gtk_main_quit();
            }
        }

    if(event->keyval == GDK_F11)
        toggle_fullscreen();

    if (simpleEditor)
    	// in simpleEditor mode, there's no other tab that can handle pressed keys, so we can send the event to editor panel then return
    	return epanel->handleShortcutKey (event);;

	if (ctrl) {
		switch(event->keyval) {
			case GDK_F2: // file browser panel
				mainNB->set_current_page (mainNB->page_num (*fpanel));
				return true;
			case GDK_F3: // batch queue panel
				mainNB->set_current_page (mainNB->page_num (*bpanel));
				return true;
			case GDK_F4: //single tab mode, editor panel
				if (isSingleTabMode() && epanel) {
				    mainNB->set_current_page (mainNB->page_num (*epanel));
				}
				return true;
			case GDK_w: //multi-tab mode, close editor panel
				if (!isSingleTabMode() && 
					mainNB->get_current_page()!=mainNB->page_num(*fpanel) &&
					mainNB->get_current_page()!=mainNB->page_num(*bpanel)) {
					
			    	EditorPanel* ep = static_cast<EditorPanel*>(mainNB->get_nth_page (mainNB->get_current_page()));
    				remEditorPanel (ep);
		    		return true;
	    		}
		}
	}

    if (mainNB->get_current_page() == mainNB->page_num(*fpanel)) {
        return fpanel->handleShortcutKey (event);
    }
    else if (mainNB->get_current_page() == mainNB->page_num(*bpanel)) {
        return bpanel->handleShortcutKey (event);
    }
    else {
        EditorPanel* ep = static_cast<EditorPanel*>(mainNB->get_nth_page (mainNB->get_current_page()));
        return ep->handleShortcutKey (event);
    }
    return false;
}

void RTWindow::addBatchQueueJob (BatchQueueEntry* bqe, bool head) {

	std::vector<BatchQueueEntry*> entries;
	entries.push_back(bqe);
    bpanel->addBatchQueueJobs (entries, head);
    fpanel->queue_draw ();
}

void RTWindow::addBatchQueueJobs (std::vector<BatchQueueEntry*> &entries) {

    bpanel->addBatchQueueJobs (entries, false);
    fpanel->queue_draw ();
}

bool RTWindow::on_delete_event(GdkEventAny* event) {

    if (on_delete_has_run) {
        // on Mac OSX we can get multiple events
        return false;
    }
    // Check if any editor is still processing, and do NOT quit if so. Otherwise crashes and inconsistent caches
    bool isProcessing=false;

    if (isSingleTabMode() || simpleEditor) 
        isProcessing=epanel->getIsProcessing();
    else {
        int pageCount=mainNB->get_n_pages();

        for (int i=0;i<pageCount && !isProcessing;i++) {
			if(isEditorPanel(i))
				isProcessing |= (static_cast<EditorPanel*>(mainNB->get_nth_page(i)))->getIsProcessing();
		}
    }

    if (isProcessing)
		return true;

    if( fpanel )
       fpanel->saveOptions ();
    if( bpanel )
       bpanel->saveOptions ();
  
    if (isSingleTabMode() || simpleEditor) {
       epanel->saveProfile();
       epanel->writeOptions ();
    }
    else {
        // Storing the options of the last EditorPanel before Gtk destroys everything
        // Look at the active panel first, if any, otherwise look at the first one (sorted on the filename)
        if (epanels.size()) {
            int page = mainNB->get_current_page();
            Gtk::Widget *w = mainNB->get_nth_page(page);
            bool optionsWritten = false;
            for (std::map<Glib::ustring, EditorPanel*>::iterator i=epanels.begin(); i!=epanels.end(); i++) {
                if (i->second == w) {
                    i->second->writeOptions();
                    optionsWritten = true;
                }
            }
            if (!optionsWritten) {
                // fallback solution: save the options of the first editor panel
                std::map<Glib::ustring, EditorPanel*>::iterator i=epanels.begin();
                i->second->writeOptions();
            }
    	}
    }

    cacheMgr->closeCache ();  // also makes cleanup if too large
    WhiteBalance::cleanup();
    ProfilePanel::cleanup();


    if (!options.windowMaximized) {
		options.windowWidth = get_width();
		options.windowHeight = get_height();
		get_position (options.windowX,options.windowY);
    }

    Options::save ();
    hide();
    on_delete_has_run = true;
    return false;
}

void RTWindow::showPreferences () {
  Preferences *pref = new Preferences (this);
  pref->run ();
  delete pref;
  
  fpanel->optionsChanged ();
}

void RTWindow::setProgress (double p) {
	prProgBar.set_fraction (p);
}

void RTWindow::setProgressStr (Glib::ustring str) {
	if (!options.mainNBVertical)
	   prProgBar.set_text ( str );
}

void RTWindow::setProgressState (bool inProcessing) {
	if (inProcessing)
		prProgBar.show();
	else
		prProgBar.hide();
	}
		
void RTWindow::toggle_fullscreen () {
    if (is_fullscreen) {
        unfullscreen();
        is_fullscreen = false;
        if (btn_fullscreen) {
            //btn_fullscreen->set_label(M("MAIN_BUTTON_FULLSCREEN"));
            btn_fullscreen->set_tooltip_markup(M("MAIN_BUTTON_FULLSCREEN"));
            btn_fullscreen->set_image (*iFullscreen);
        }
    } else {
        fullscreen();
        is_fullscreen = true;
        if (btn_fullscreen) {
            //btn_fullscreen->set_label(M("MAIN_BUTTON_UNFULLSCREEN"));
            btn_fullscreen->set_tooltip_markup(M("MAIN_BUTTON_UNFULLSCREEN"));
            btn_fullscreen->set_image (*iFullscreen_exit);
        }
    }
}

void RTWindow::error (Glib::ustring descr){
	prProgBar.set_text ( descr );
}

void RTWindow::SetEditorCurrent()
{
  mainNB->set_current_page (mainNB->page_num (*epanel));
}

void RTWindow::SetMainCurrent()
{
  mainNB->set_current_page (mainNB->page_num (*fpanel));
}

void RTWindow::MoveFileBrowserToMain()
{
    if( fpanel->ribbonPane->get_children().empty())
    {
        FileCatalog *fCatalog = fpanel->fileCatalog;
        epanel->catalogPane->remove(*fCatalog);
        fpanel->ribbonPane->add(*fCatalog);
        fCatalog->enableTabMode(false);
        fCatalog->tbLeftPanel_1_visible(true);
        fCatalog->tbRightPanel_1_visible(true);
    }
}

void RTWindow::MoveFileBrowserToEditor()
{
    if(epanel->catalogPane->get_children().empty() )
    {
        FileCatalog *fCatalog = fpanel->fileCatalog;
        fpanel->ribbonPane->remove(*fCatalog);
        epanel->catalogPane->add(*fCatalog);
        epanel->showTopPanel(options.editorFilmStripOpened);
        fCatalog->enableTabMode(true);
        fCatalog->refreshHeight();
        fCatalog->tbLeftPanel_1_visible(false);
        fCatalog->tbRightPanel_1_visible(false);
    }
}

void RTWindow::updateTPProfileSelector (bool showMe) {
	epanel->updateProfileSelector (showMe);

	std::map<Glib::ustring, EditorPanel*>::const_iterator itr;
	for(itr = epanels.begin(); itr != epanels.end(); ++itr){
		((*itr).second)->updateProfileSelector (showMe);
	}
}

void RTWindow::updateTPVScrollbar (bool hide) {
	fpanel->updateTPVScrollbar (hide);
	epanel->updateTPVScrollbar (hide);

	std::map<Glib::ustring, EditorPanel*>::const_iterator itr;
	for(itr = epanels.begin(); itr != epanels.end(); ++itr){
		((*itr).second)->updateTPVScrollbar (hide);
	}
}

void RTWindow::updateTabsUsesIcons (bool useIcons) {
	fpanel->updateTabsUsesIcons (useIcons);
	epanel->updateTabsUsesIcons (useIcons);

	std::map<Glib::ustring, EditorPanel*>::const_iterator itr;
	for(itr = epanels.begin(); itr != epanels.end(); ++itr){
		((*itr).second)->updateTabsUsesIcons (useIcons);
	}
}

void RTWindow::updateFBQueryTB (bool singleRow) {
	fpanel->fileCatalog->updateFBQueryTB (singleRow);
}

void RTWindow::updateHistogramPosition (int oldPosition, int newPosition) {
	epanel->updateHistogramPosition (oldPosition, newPosition);

	std::map<Glib::ustring, EditorPanel*>::const_iterator itr;
	for(itr = epanels.begin(); itr != epanels.end(); ++itr){
		((*itr).second)->updateHistogramPosition (oldPosition, newPosition);
	}
}

bool RTWindow::splashClosed(GdkEventAny* event) {
	delete splash;
	splash = NULL;
	return true;
}

void RTWindow::set_title_decorated(Glib::ustring fname){
	Glib::ustring subtitle;
	if (!fname.empty()) subtitle= " - " + fname;
	set_title(versionStr + subtitle);	
}

void RTWindow::CloseOpenEditors(){
	std::map<Glib::ustring, EditorPanel*>::const_iterator itr;
	itr = epanels.begin();
	while(itr != epanels.end()) {
		remEditorPanel((*itr).second);
		itr = epanels.begin();
	}
}

bool RTWindow::isEditorPanel(Widget* panel) {
	return (panel != bpanel) && (panel != fpanel);
}

bool RTWindow::isEditorPanel(guint pageNum) {
	return isEditorPanel(mainNB->get_nth_page(pageNum));
}