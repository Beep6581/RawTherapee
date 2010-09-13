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
#include <rtwindow.h>
#include <options.h>
#include <preferences.h>
#include <cursormanager.h>

RTWindow::RTWindow () {

    cacheMgr.init ();

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
    property_destroy_with_parent().set_value(false);

    mainNB = Gtk::manage (new Gtk::Notebook ());
    mainNB->set_scrollable (true);

    fpanel = new FilePanel ();
    fpanel->setParent (this);

    // decorate tab
    Gtk::HBox* hbf = Gtk::manage (new Gtk::HBox ());
    hbf->pack_start (*Gtk::manage (new Gtk::Image (Gtk::Stock::DIRECTORY, Gtk::ICON_SIZE_MENU)));
    hbf->pack_start (*Gtk::manage (new Gtk::Label (M("MAIN_FRAME_FILEBROWSER"))));
    hbf->set_spacing (2);
    hbf->show_all ();
    mainNB->append_page (*fpanel, *hbf);

    bpanel = new BatchQueuePanel ();
    bpanel->setParent (this);

    // decorate tab
    Gtk::HBox* hbb = Gtk::manage (new Gtk::HBox ());
    hbb->pack_start (*Gtk::manage (new Gtk::Image (Gtk::Stock::EXECUTE, Gtk::ICON_SIZE_MENU)));
    hbb->pack_start (*Gtk::manage (new Gtk::Label (M("MAIN_FRAME_BATCHQUEUE"))));
    hbb->set_spacing (2);
    hbb->show_all ();
    mainNB->append_page (*bpanel, *hbb);

    signal_key_press_event().connect( sigc::mem_fun(*this, &RTWindow::keyPressed) );

    Gtk::VBox* mainBox = Gtk::manage (new Gtk::VBox ());
    mainBox->pack_start (*mainNB);   
    Gtk::HBox* bottomBox = Gtk::manage (new Gtk::HBox ());
    mainBox->pack_start (*bottomBox, Gtk::PACK_SHRINK, 1);

    // filling bottom box
    Gtk::LinkButton* rtWeb = Gtk::manage (new Gtk::LinkButton ("http://rawtherapee.com"));
    Gtk::Button* preferences = Gtk::manage (new Gtk::Button (M("MAIN_BUTTON_PREFERENCES")));
    preferences->set_image (*Gtk::manage(new Gtk::Image (Gtk::StockID("gtk-preferences"), Gtk::ICON_SIZE_BUTTON)));
    preferences->set_relief (Gtk::RELIEF_NONE);
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
    preferences->modify_style (style);

    add (*mainBox);
    show_all ();
}

void RTWindow::on_realize () {

    Gtk::Window::on_realize ();

    cursorManager.init (get_window());
}

void RTWindow::addEditorPanel (EditorPanel* ep, const std::string &name) {

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
    mainNB->set_current_page (mainNB->page_num (*ep));
    mainNB->set_tab_reorderable (*ep, true);

    epanels[ name ] = ep;
    filesEdited.insert ( name );
    fpanel->refreshEditedState (filesEdited);
}

void RTWindow::remEditorPanel (EditorPanel* ep) {

    if (ep->beforeClosing ()) {
        ep->saveOptions ();
        epanels.erase (ep->getFileName());
        filesEdited.erase (ep->getFileName ());
        fpanel->refreshEditedState (filesEdited);

        mainNB->remove_page (*ep);
       
        if (mainNB->get_current_page () == mainNB->page_num (*bpanel))
            mainNB->set_current_page (mainNB->page_num (*fpanel));
    }
    // TODO: ask what to do: close & apply, close & apply selection, close & revert, cancel
}

bool RTWindow::keyPressed (GdkEventKey* event) {
    if(event->keyval == GDK_F11) {
        toggle_fullscreen();
    }

    if (mainNB->get_nth_page (mainNB->get_current_page()) == fpanel) {
    }
//    else if (mainNB->get_nth_page (mainNB->get_current_page()) == bqpanel) {
//    }
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

    bpanel->addBatchQueueJob (bqe, head);
    fpanel->queue_draw ();
}

bool RTWindow::on_delete_event(GdkEventAny* event) {

    fpanel->saveOptions ();
    bpanel->saveOptions ();

/*    if (fileBrowser->getFileCatalog()->getBatchQueue()->hasJobs()) {
        Gtk::MessageDialog msgd (M("MAIN_MSG_EXITJOBSINQUEUEQUEST"), false, Gtk::MESSAGE_QUESTION, Gtk::BUTTONS_YES_NO, true);
        msgd.set_secondary_text (M("MAIN_MSG_EXITJOBSINQUEUEINFO"));
        int response = msgd.run ();
        if (response==Gtk::RESPONSE_NO)
            return true;
    }

    editCoord->close ();

    if (options.startupDir==STARTUPDIR_LAST && fileBrowser->lastSelectedDir ()!="")
        options.startupPath = fileBrowser->lastSelectedDir ();
    fileBrowser->close ();
    cacheMgr.closeCache ();
        
      
    options.lastScale = editorPanel->zoomBar->getScale ();
    options.lastCropSize = editorPanel->zoomBar->getCropSize ();
    if (options.showFilePanelState==0 || options.showFilePanelState==2)
        options.fileBrowserHeight = fileBrowser->get_height ();
    options.historyPanelWidth = ppaned->get_position ();
    options.toolPanelWidth = vboxright->get_width();//hpaned->get_width() - hpaned->get_position ();
    options.showHistory  = editorPanel->hidehp->get_active ();
    options.showInfo = editorPanel->info->get_active ();
    options.showClippedHighlights = editorPanel->indclippedh->get_active ();
    options.showClippedShadows = editorPanel->indclippeds->get_active ();
    options.bgcolor = editorPanel->iarea->imageArea->getBGColor ();
    options.lastSaveAsPath = saveAsDialog->getDirectory ();
    options.procQueueEnabled = fileBrowser->getFileCatalog()->getBatchQueue()->isEnabled();
    options.fbArrangement = fileBrowser->getFileCatalog()->getArrangement ();
    options.firstRun = false;
*/
    //options.windowWidth = get_width();
    //options.windowHeight = get_height();
   

    Options::save ();
    hide();
    return true;
}

void RTWindow::showPreferences () {

  Preferences *pref = new Preferences ();
  pref->run ();
  delete pref;
  
  fpanel->optionsChanged ();
}
void RTWindow::setProgress (double p){
	prProgBar.set_fraction (p);
}
void RTWindow::setProgressStr (Glib::ustring str){
	prLabel.set_text ( str );
}
void RTWindow::setProgressState (int state){
	if(state){
		prProgBar.show();
		prLabel.show();
	}else{
		prProgBar.hide();
		prLabel.hide();
	}
}
		
void RTWindow::toggle_fullscreen () {
    if(is_fullscreen){
        unfullscreen();
        is_fullscreen = false;
        btn_fullscreen->set_label(M("MAIN_BUTTON_FULLSCREEN"));
    }
    else {
        fullscreen();
        is_fullscreen = true;
        btn_fullscreen->set_label(M("MAIN_BUTTON_UNFULLSCREEN"));
    }
}

void RTWindow::error (Glib::ustring descr){
	prLabel.set_text ( descr );
}
