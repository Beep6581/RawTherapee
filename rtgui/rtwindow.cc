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
#include <rtwindow.h>
#include <options.h>
#include <preferences.h>
#include <cursormanager.h>


RTWindow::RTWindow () {

    epanel=NULL;  // to prevent eventing errors

    cacheMgr->init ();

#ifdef GLIBMM_EXCEPTIONS_ENABLED
		try { set_default_icon_from_file (argv0+"/images/logoicon16.png");
		} catch(Glib::Exception& ex) {		printf ("%s\n", ex.what().c_str());	}
#else
 	  {		std::auto_ptr<Glib::Error> error;
				set_default_icon_from_file (argv0+"/images/logoicon16.png", error);
		}
#endif //GLIBMM_EXCEPTIONS_ENABLED

    set_title("RawTherapee "+versionString);
    property_allow_shrink() = true;
    set_default_size(options.windowWidth, options.windowHeight);
    set_modal(false);
    set_resizable(true);
    if (options.windowMaximized)
    	maximize();
    else
    	unmaximize();
    property_destroy_with_parent().set_value(false);
    signal_window_state_event().connect( sigc::mem_fun(*this, &RTWindow::on_window_state_event) );

    mainNB = Gtk::manage (new Gtk::Notebook ());
    mainNB->set_scrollable (true);
    mainNB->signal_switch_page().connect_notify( sigc::mem_fun(*this, &RTWindow::on_mainNB_switch_page) );

    fpanel = Gtk::manage ( new FilePanel () );
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
        vbf->show_all ();
        mainNB->append_page (*fpanel, *vbf);
     } else {
        Gtk::HBox* hbf = Gtk::manage (new Gtk::HBox ());
        hbf->pack_start (*Gtk::manage (new Gtk::Image (Gtk::Stock::DIRECTORY, Gtk::ICON_SIZE_MENU)));
        hbf->pack_start (*Gtk::manage (new Gtk::Label (M("MAIN_FRAME_FILEBROWSER"))));
        hbf->set_spacing (2);
        hbf->show_all ();
        mainNB->append_page (*fpanel, *hbf);
    }

    bpanel = Gtk::manage ( new BatchQueuePanel () );
    bpanel->setParent (this);

    // decorate tab, the label is unimportant since its updated in batchqueuepanel anyway
    Gtk::Label* lbq = Gtk::manage ( new Gtk::Label (M("MAIN_FRAME_BATCHQUEUE")) );
    mainNB->append_page (*bpanel, *lbq);
    
    // epanel is only for single tab mode
    epanel = Gtk::manage ( new EditorPanel (fpanel) );
    epanel->setParent (this);

    // decorate tab
    if (options.mainNBVertical) {
        Gtk::VBox* vbe = Gtk::manage (new Gtk::VBox ());
        vbe->pack_start (*Gtk::manage (new Gtk::Image (argv0+"/images/logoicon16.png")));
        Gtk::Label* l=Gtk::manage (new Gtk::Label( Glib::ustring(" ") + M("MAIN_FRAME_EDITOR") ));
        //l->set_markup(Glib::ustring("<b>Editor</b>"));  Bold difficult to read
        l->set_angle (90);
        vbe->pack_start (*l);
        vbe->set_spacing (2);
        vbe->show_all ();
        mainNB->append_page (*epanel, *vbe);
    } else {
        Gtk::HBox* hbe = Gtk::manage (new Gtk::HBox ());
        hbe->pack_start (*Gtk::manage (new Gtk::Image (argv0+"/images/logoicon16.png")));
        hbe->pack_start (*Gtk::manage (new Gtk::Label(M("MAIN_FRAME_EDITOR"))));
        hbe->set_spacing (2);
        hbe->show_all ();
        mainNB->append_page (*epanel, *hbe);
    }

    mainNB->set_current_page (mainNB->page_num (*fpanel));

    signal_key_press_event().connect( sigc::mem_fun(*this, &RTWindow::keyPressed) );

    Gtk::VBox* mainBox = Gtk::manage (new Gtk::VBox ());
    mainBox->pack_start (*mainNB);   
    Gtk::HBox* bottomBox = Gtk::manage (new Gtk::HBox ());
    mainBox->pack_start (*bottomBox, Gtk::PACK_SHRINK, 1);

    // filling bottom box
    Gtk::LinkButton* rtWeb = Gtk::manage (new Gtk::LinkButton ("http://rawtherapee.com"));
    Gtk::Button* preferences = Gtk::manage (new Gtk::Button (M("MAIN_BUTTON_PREFERENCES")+"..."));
    preferences->set_image (*Gtk::manage(new Gtk::Image (Gtk::StockID("gtk-preferences"), Gtk::ICON_SIZE_BUTTON)));
    preferences->signal_clicked().connect( sigc::mem_fun(*this, &RTWindow::showPreferences) );
    is_fullscreen = false;
    btn_fullscreen = Gtk::manage( new Gtk::Button(M("MAIN_BUTTON_FULLSCREEN")));
    btn_fullscreen->signal_clicked().connect( sigc::mem_fun(*this, &RTWindow::toggle_fullscreen) );
    bottomBox->pack_start (*preferences, Gtk::PACK_SHRINK, 0);
    bottomBox->pack_end (*btn_fullscreen, Gtk::PACK_SHRINK, 4);
    bottomBox->pack_start (*rtWeb, Gtk::PACK_SHRINK, 4);
    bottomBox->pack_start (prLabel );
    prLabel.set_alignment(Gtk::ALIGN_RIGHT);
    bottomBox->pack_start (prProgBar, Gtk::PACK_SHRINK, 4);
    pldBridge = new PLDBridge(&prLabel,&prProgBar);

    Glib::RefPtr<Gtk::RcStyle> style = Gtk::RcStyle::create ();
    style->set_xthickness (0);
    style->set_ythickness (0);    
    rtWeb->modify_style (style);

    add (*mainBox);
    show_all ();

    if (!isSingleTabMode()) epanel->hide_all();
}

void RTWindow::on_realize () {
    Gtk::Window::on_realize ();

    fpanel->setAspect();

    cursorManager.init (get_window());
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
	if (page_num > 1) {
        if (isSingleTabMode()) MoveFileBrowserToEditor();

		EditorPanel *ep = (EditorPanel *)mainNB->get_nth_page(page_num);
		ep->setAspect();
	} else {
        // in single tab mode with command line filename epanel does not exist yet
        if (isSingleTabMode() && epanel) {
            // Save profile on leaving the editor panel
            epanel->saveProfile();

            MoveFileBrowserToMain();
        }
    }
}

void RTWindow::addEditorPanel (EditorPanel* ep, const std::string &name) {
    if (EditWindow::isMultiDisplayEnabled()) {
        EditWindow * wndEdit = EditWindow::getInstance(this);
        wndEdit->show_all();
        wndEdit->addEditorPanel(ep,name);
    } else {
        ep->setParent (this);

        // construct closeable tab for the image
        Gtk::HBox* hb = Gtk::manage (new Gtk::HBox ());
        hb->pack_start (*Gtk::manage (new Gtk::Image (Gtk::Stock::FILE, Gtk::ICON_SIZE_MENU)));
        hb->pack_start (*Gtk::manage (new Gtk::Label (name)));
        Gtk::Button* closeb = Gtk::manage (new Gtk::Button ());
        closeb->set_image (*Gtk::manage(new Gtk::Image (Gtk::Stock::CLOSE, Gtk::ICON_SIZE_MENU)));
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

        epanels[ name ] = ep;
        filesEdited.insert ( name );
        fpanel->refreshEditedState (filesEdited);
    }
}

void RTWindow::remEditorPanel (EditorPanel* ep) {
    if (EditWindow::isMultiDisplayEnabled()) {
        EditWindow * wndEdit = EditWindow::getInstance(this);
        wndEdit->remEditorPanel(ep);
    } else {
	    epanels.erase (ep->getShortName());
	    filesEdited.erase (ep->getShortName ());
	    fpanel->refreshEditedState (filesEdited);

	    mainNB->remove_page (*ep);

	    if (mainNB->get_current_page () == mainNB->page_num (*bpanel))
		    mainNB->set_current_page (mainNB->page_num (*fpanel));
        // TODO: ask what to do: close & apply, close & apply selection, close & revert, cancel
    }
}

bool RTWindow::keyPressed (GdkEventKey* event) {
    if(event->keyval == GDK_F11) {
        toggle_fullscreen();
    }

    if (mainNB->get_current_page() == mainNB->page_num(*fpanel)) {
        return fpanel->handleShortcutKey (event);
    }
    else if (mainNB->get_current_page() == mainNB->page_num(*bpanel)) {
        return false;
    }
    else {
        EditorPanel* ep = (EditorPanel*)mainNB->get_nth_page (mainNB->get_current_page());
        return ep->handleShortcutKey (event);
    }
    return false;
}

void RTWindow::imageDeveloped (Glib::ustring fname) {

//    fpanel->refreshThumbnail (fname);
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


    fpanel->saveOptions ();
    bpanel->saveOptions ();
  
    if (isSingleTabMode()) epanel->saveProfile();
    
    cacheMgr->closeCache ();  // also makes cleanup if too large

    options.firstRun = false;

    if (!options.windowMaximized) {
		options.windowWidth = get_width();
		options.windowHeight = get_height();
    }

    Options::save ();
    hide();
    return true;
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
	prLabel.set_text ( str );
}

void RTWindow::setProgressState (int state) {
	if (state) {
		prProgBar.show();
		prLabel.show();
	} else {
		prProgBar.hide();
		prLabel.hide();
	}
}
		
void RTWindow::toggle_fullscreen () {
    if (is_fullscreen) {
        unfullscreen();
        is_fullscreen = false;
        btn_fullscreen->set_label(M("MAIN_BUTTON_FULLSCREEN"));
    } else {
        fullscreen();
        is_fullscreen = true;
        btn_fullscreen->set_label(M("MAIN_BUTTON_UNFULLSCREEN"));
    }
}

void RTWindow::error (Glib::ustring descr){
	prLabel.set_text ( descr );
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
    if( fpanel->ribbonPane->get_children().size() ==0)
    {
        FileCatalog *fCatalog = fpanel->fileCatalog;
        epanel->catalogPane->remove(*fCatalog);
        fpanel->ribbonPane->add(*fCatalog);
        fCatalog->enableTabMode(false);
    }
}

void RTWindow::MoveFileBrowserToEditor()
{
    if(epanel->catalogPane->get_children().size() ==0 )
    {
        FileCatalog *fCatalog = fpanel->fileCatalog;
        fpanel->ribbonPane->remove(*fCatalog);
        epanel->catalogPane->add(*fCatalog);
        fCatalog->enableTabMode(true);
        fCatalog->refreshHeight();
    }
}

