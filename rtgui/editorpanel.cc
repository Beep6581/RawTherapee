/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Copyright (c) 2010 Oliver Duis <www.oliverduis.de>
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
#include <editorpanel.h>
#include <options.h>
#include <progressconnector.h>
#include <rtwindow.h>
#include <guiutils.h>
#include <procparamchangers.h>
#include <safegtk.h>
#include <imagesource.h>
#include <soundman.h>

using namespace rtengine::procparams;

EditorPanel::EditorPanel (FilePanel* filePanel) 
    : beforePreviewHandler(NULL), beforeIarea(NULL), parent(NULL), beforeIpc(NULL), ipc(NULL), catalogPane(NULL), isProcessing(false) {

    epih = new EditorPanelIdleHelper;
    epih->epanel = this;
    epih->destroyed = false;
    epih->pending = 0;

    processingStartedTime = 0;
    firstProcessingDone = false;

// construct toolpanelcoordinator
    tpc = new ToolPanelCoordinator ();

// build GUI
    // build left side panel
    leftbox = new Gtk::VBox ();
    leftbox->set_border_width (2);
    leftbox->set_size_request(100,250);

    if (options.histogramPosition>0) {
    histogramPanel = Gtk::manage (new HistogramPanel ());
    histogramPanel->set_size_request (-1, 160);
        if (options.histogramPosition==1) leftbox->pack_start (*histogramPanel, Gtk::PACK_SHRINK, 4);
    } else histogramPanel = NULL;

    profilep = Gtk::manage (new ProfilePanel ());
    Gtk::Frame* ppframe = Gtk::manage (new Gtk::Frame ());
    ppframe->add (*profilep);
    ppframe->set_label (M("PROFILEPANEL_LABEL"));
//    leftbox->pack_start (*ppframe, Gtk::PACK_SHRINK, 4);

    navigator = Gtk::manage (new Navigator ());
    navigator->previewWindow->set_size_request (-1, 150);
    leftbox->pack_start (*navigator, Gtk::PACK_SHRINK, 2);

    history = Gtk::manage (new History ());
    leftbox->pack_start (*history);

    leftbox->show_all ();

    // build the middle of the screen
    Gtk::VBox* editbox = Gtk::manage (new Gtk::VBox ());

    info = Gtk::manage (new Gtk::ToggleButton ());
    Gtk::Image* infoimg = Gtk::manage (new Gtk::Image (argv0+"/images/info.png"));
    info->add (*infoimg);
    info->set_relief(Gtk::RELIEF_NONE);
    info->set_tooltip_markup (M("MAIN_TOOLTIP_QINFO"));

    beforeAfter = Gtk::manage (new Gtk::ToggleButton ());
    Gtk::Image* beforeAfterIcon = Gtk::manage (new Gtk::Image (argv0+"/images/beforeafter.png"));
    beforeAfter->add(*beforeAfterIcon);
    beforeAfter->set_relief(Gtk::RELIEF_NONE);
    beforeAfter->set_tooltip_markup (M("MAIN_TOOLTIP_TOGGLE"));


    Gtk::VSeparator* vsept = Gtk::manage (new Gtk::VSeparator ());
    Gtk::VSeparator* vsepz = Gtk::manage (new Gtk::VSeparator ());
    Gtk::VSeparator* vsepi = Gtk::manage (new Gtk::VSeparator ());
    Gtk::VSeparator* vseph = Gtk::manage (new Gtk::VSeparator ());
    Gtk::VSeparator* vsep1 = Gtk::manage (new Gtk::VSeparator ());

    hidehp = Gtk::manage (new Gtk::ToggleButton ());

    iHistoryShow = new Gtk::Image(argv0+"/images/panel_to_right.png");
    iHistoryHide = new Gtk::Image(argv0+"/images/panel_to_left.png");

    hidehp->set_relief(Gtk::RELIEF_NONE);
    hidehp->set_active (options.showHistory);
    hidehp->set_tooltip_markup (M("MAIN_TOOLTIP_HIDEHP"));
    if (options.showHistory){
    	hidehp->set_image (*iHistoryHide);
    }
    else {
    	hidehp->set_image (*iHistoryShow);
    }

    tbTopPanel_1 = new Gtk::ToggleButton ();
    iTopPanel_1_Show = new Gtk::Image(argv0+"/images/panel_to_bottom.png");
    iTopPanel_1_Hide = new Gtk::Image(argv0+"/images/panel_to_top.png");
    tbTopPanel_1->set_relief(Gtk::RELIEF_NONE);
    tbTopPanel_1->set_active (true);
    tbTopPanel_1->set_tooltip_markup (M("MAIN_TOOLTIP_SHOWHIDETP1"));
    tbTopPanel_1->set_image (*iTopPanel_1_Hide);

    tbRightPanel_1 = new Gtk::ToggleButton ();
    iRightPanel_1_Show = new Gtk::Image(argv0+"/images/panel_to_left.png");
    iRightPanel_1_Hide = new Gtk::Image(argv0+"/images/panel_to_right.png");
    tbRightPanel_1->set_relief(Gtk::RELIEF_NONE);
    tbRightPanel_1->set_active (true);
    tbRightPanel_1->set_tooltip_markup (M("MAIN_TOOLTIP_SHOWHIDERP1"));
    tbRightPanel_1->set_image (*iRightPanel_1_Hide);

    Gtk::VSeparator* vsepcl = Gtk::manage (new Gtk::VSeparator ());
    Gtk::VSeparator* vsepz2 = Gtk::manage (new Gtk::VSeparator ());
    Gtk::VSeparator* vsepz3 = Gtk::manage (new Gtk::VSeparator ());

    iarea = new ImageAreaPanel ();

    Gtk::HBox* toolBarPanel = Gtk::manage (new Gtk::HBox ());
    toolBarPanel->pack_start (*hidehp, Gtk::PACK_SHRINK, 1);
    toolBarPanel->pack_start (*vseph, Gtk::PACK_SHRINK, 2);
    toolBarPanel->pack_start (*info, Gtk::PACK_SHRINK, 1);
    toolBarPanel->pack_start (*beforeAfter, Gtk::PACK_SHRINK, 1);
    toolBarPanel->pack_start (*vsepi, Gtk::PACK_SHRINK, 2);
    toolBarPanel->pack_start (*tpc->getToolBar(), Gtk::PACK_SHRINK, 1);
    toolBarPanel->pack_start (*vsept, Gtk::PACK_SHRINK, 2);
    
    toolBarPanel->pack_end   (*tbTopPanel_1, Gtk::PACK_SHRINK, 1);
    toolBarPanel->pack_end   (*vsep1, Gtk::PACK_SHRINK, 2);
    toolBarPanel->pack_end   (*tpc->coarse, Gtk::PACK_SHRINK, 2);
    toolBarPanel->pack_end   (*vsepcl, Gtk::PACK_SHRINK, 2);
    toolBarPanel->pack_end   (*iarea->imageArea->indClippedPanel, Gtk::PACK_SHRINK, 0);
    toolBarPanel->pack_end   (*vsepz, Gtk::PACK_SHRINK, 2);

    afterBox = Gtk::manage (new Gtk::VBox ());
    afterBox->pack_start (*iarea);

    beforeAfterBox = Gtk::manage (new Gtk::HBox());
    beforeAfterBox->pack_start (*afterBox);

    editbox->pack_start (*toolBarPanel, Gtk::PACK_SHRINK,0);
    editbox->pack_start (*beforeAfterBox);

    // build right side panel
    vboxright = new Gtk::VBox (false, 0);
    vboxright->set_size_request(100,250);

    vboxright->set_border_width (2);

    if (options.histogramPosition==2) vboxright->pack_start (*histogramPanel, Gtk::PACK_SHRINK, 2);
    if (options.showProfileSelector) vboxright->pack_start (*ppframe, Gtk::PACK_SHRINK, 2);
    // main notebook
    vboxright->pack_start (*tpc->toolPanelNotebook);

    // Save buttons
    Gtk::HBox* iops = Gtk::manage (new Gtk::HBox ());

    //Gtk::Image *saveButtonImage = Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-save"), Gtk::ICON_SIZE_BUTTON));
    Gtk::Image *saveButtonImage =  Gtk::manage (new Gtk::Image (argv0+"/images/save_hdd_01.png"));
    saveimgas = Gtk::manage (new Gtk::Button ());
    saveimgas->add(*saveButtonImage);
    saveimgas->set_tooltip_markup(M("MAIN_BUTTON_SAVE_TOOLTIP"));

    Gtk::Image *queueButtonImage = Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-execute"), Gtk::ICON_SIZE_BUTTON));
    queueimg = Gtk::manage (new Gtk::Button ());
    queueimg->add(*queueButtonImage);
    queueimg->set_tooltip_markup(M("MAIN_BUTTON_PUTTOQUEUE_TOOLTIP"));

    Gtk::Image *sendToEditorButtonImage = Gtk::manage (new Gtk::Image (argv0+"/images/gimp.png"));
    sendtogimp = Gtk::manage (new Gtk::Button ());
    sendtogimp->add(*sendToEditorButtonImage);
    sendtogimp->set_tooltip_markup(M("MAIN_BUTTON_SENDTOEDITOR_TOOLTIP"));

    iops->pack_start (*saveimgas, Gtk::PACK_SHRINK);
    if(!simpleEditor)
       iops->pack_start (*queueimg, Gtk::PACK_SHRINK);
    iops->pack_start (*sendtogimp, Gtk::PACK_SHRINK);

    // Status box
    statusBox = Gtk::manage (new Gtk::HBox ());
    progressLabel = Gtk::manage (new Gtk::ProgressBar());
    progressLabel->set_fraction(0.0);
    //progressLabel->modify_bg( Gtk::STATE_NORMAL,Gdk::Color("grey") );  // Disable, because in single mode this is may be permanent red without processing

    statusBox->pack_start (*progressLabel);
    iops->pack_start(*statusBox, Gtk::PACK_SHRINK, 2);

    // tbRightPanel_1
    iops->pack_end (*tbRightPanel_1, Gtk::PACK_SHRINK,0);

	// ShowHideSidePanels
    tbShowHideSidePanels = new Gtk::ToggleButton ();
    iShowHideSidePanels = new Gtk::Image(argv0+"/images/crossed_arrows_out_45_02.png");
    tbShowHideSidePanels->set_relief(Gtk::RELIEF_NONE);
    tbShowHideSidePanels->set_active (false);
    tbShowHideSidePanels->set_tooltip_markup (M("MAIN_BUTTON_SHOWHIDESIDEPANELS_TOOLTIP"));
    tbShowHideSidePanels->set_image (*iShowHideSidePanels);
    iops->pack_end (*tbShowHideSidePanels, Gtk::PACK_SHRINK,0);
	iops->pack_end (*vsepz2, Gtk::PACK_SHRINK,1);

    // Zoom panel
    iops->pack_end (*iarea->imageArea->zoomPanel, Gtk::PACK_SHRINK, 1);
    iops->pack_end (*vsepz3, Gtk::PACK_SHRINK, 2);

    editbox->pack_start (*Gtk::manage(new Gtk::HSeparator()), Gtk::PACK_SHRINK, 0);
    editbox->pack_start (*iops, Gtk::PACK_SHRINK, 0);
    editbox->show_all ();

    // build screen
    hpanedl = Gtk::manage (new Gtk::HPaned());
    hpanedr = Gtk::manage (new Gtk::HPaned());
    leftbox->reference ();
    vboxright->reference ();
    if (options.showHistory) {
       hpanedl->pack1(*leftbox, false, true);
       hpanedl->set_position (options.historyPanelWidth);
    }

  
    Gtk::VPaned * viewpaned = Gtk::manage (new Gtk::VPaned());
    fPanel = filePanel;
    if(filePanel)
    {
        catalogPane = new Gtk::Paned();
        viewpaned->pack1(*catalogPane, false, true);                    
    }
    viewpaned->pack2(*editbox, true, true);


    Gtk::Frame* vbfr = Gtk::manage (new Gtk::Frame ());    
    vbfr->add (*viewpaned);
    vbfr->set_size_request(100,250);
    hpanedl->pack2(*vbfr, true, true);

    hpanedr->pack1(*hpanedl, true, true);
    hpanedr->pack2(*vboxright, false, true);
	hpanedl->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &EditorPanel::leftPaneButtonReleased) );
	hpanedr->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &EditorPanel::rightPaneButtonReleased) );

    pack_start (*hpanedr);
    show_all ();

    // save as dialog
    if (safe_file_test (options.lastSaveAsPath, Glib::FILE_TEST_IS_DIR))
        saveAsDialog = new SaveAsDialog (options.lastSaveAsPath);
    else
        saveAsDialog = new SaveAsDialog (Glib::get_user_special_dir (G_USER_DIRECTORY_PICTURES));

    saveAsDialog->set_default_size (options.saveAsDialogWidth, options.saveAsDialogHeight);

// connect listeners
    profilep->setProfileChangeListener (tpc);
    history->setProfileChangeListener (tpc);
    history->setHistoryBeforeLineListener (this);
    tpc->addPParamsChangeListener (profilep);
    tpc->addPParamsChangeListener (history);
    tpc->addPParamsChangeListener (this);
    iarea->imageArea->setCropGUIListener (tpc->getCropGUIListener());
    iarea->imageArea->setPointerMotionListener (navigator);
    iarea->imageArea->setPointerMotionHListener (histogramPanel);
    iarea->imageArea->setImageAreaToolListener (tpc);

// initialize components
    info->set_active (options.showInfo);
    tpc->readOptions ();

// connect event handlers
    info->signal_toggled().connect( sigc::mem_fun(*this, &EditorPanel::info_toggled) );
    beforeAfter->signal_toggled().connect( sigc::mem_fun(*this, &EditorPanel::beforeAfterToggled) );
    hidehp->signal_toggled().connect( sigc::mem_fun(*this, &EditorPanel::hideHistoryActivated) );
    tbTopPanel_1->signal_toggled().connect( sigc::mem_fun(*this, &EditorPanel::tbTopPanel_1_toggled) );
    tbRightPanel_1->signal_toggled().connect( sigc::mem_fun(*this, &EditorPanel::tbRightPanel_1_toggled) );
    saveimgas->signal_pressed().connect( sigc::mem_fun(*this, &EditorPanel::saveAsPressed) );
    queueimg->signal_pressed().connect( sigc::mem_fun(*this, &EditorPanel::queueImgPressed) );
    sendtogimp->signal_pressed().connect( sigc::mem_fun(*this, &EditorPanel::sendToGimpPressed) );
    ShowHideSidePanelsconn = tbShowHideSidePanels->signal_toggled().connect ( sigc::mem_fun(*this, &EditorPanel::toggleSidePanels), true);
}

EditorPanel::~EditorPanel () {

    history->setHistoryBeforeLineListener (NULL);
    // the order is important!
    iarea->setBeforeAfterViews (NULL, iarea);
    delete iarea;
    iarea = NULL;

    if (beforeIpc)
        beforeIpc->stopProcessing ();

    delete beforeIarea;
    beforeIarea = NULL;

    if (beforeIpc)
        beforeIpc->setPreviewImageListener (NULL);

    delete beforePreviewHandler;
    beforePreviewHandler = NULL;
    if (beforeIpc)
        rtengine::StagedImageProcessor::destroy (beforeIpc);
    beforeIpc = NULL;

    close ();

    if (epih->pending)
        epih->destroyed = true;
    else
        delete epih;

    delete tpc;

    delete leftbox;
    delete vboxright;    
    delete saveAsDialog;
    if(catalogPane)
        delete catalogPane;
}

void EditorPanel::leftPaneButtonReleased(GdkEventButton *event) {
	if (event->button == 1) {
		// Button 1 released : it's a resize
		options.historyPanelWidth = hpanedl->get_position();
	}
	/*else if (event->button == 3) {
	}*/
}

void EditorPanel::rightPaneButtonReleased(GdkEventButton *event) {
	if (event->button == 1) {
		int winW, winH;
		parent->get_size(winW, winH);
		// Button 1 released : it's a resize
		options.toolPanelWidth = winW - hpanedr->get_position();
	}
	/*else if (event->button == 3) {
	}*/
}

void EditorPanel::setAspect () {
	int winW, winH;
	parent->get_size(winW, winH);
	hpanedl->set_position(options.historyPanelWidth);
	hpanedr->set_position(winW - options.toolPanelWidth);
	// initialize components
	if (info->get_active() != options.showInfo)
		info->set_active (options.showInfo);
}

void EditorPanel::on_realize () {
    
    Gtk::VBox::on_realize ();
    // This line is needed to avoid autoexpansion of the window :-/
    vboxright->set_size_request (options.toolPanelWidth, -1);
}

void EditorPanel::open (Thumbnail* tmb, rtengine::InitialImage* isrc) {

    close();

    isProcessing=true;  // prevents closing-on-init

    // initialize everything
    openThm = tmb;
    openThm->increaseRef ();

    fname=openThm->getFileName();

    previewHandler = new PreviewHandler ();

    this->isrc = isrc;
    ipc = rtengine::StagedImageProcessor::create (isrc);
    ipc->setProgressListener (this);
    ipc->setPreviewImageListener (previewHandler);
    ipc->setPreviewScale (10);  // Important
    tpc->initImage (ipc, tmb->getType()==FT_Raw);
    ipc->setHistogramListener (this);

//    iarea->fitZoom ();   // tell to the editorPanel that the next image has to be fitted to the screen
    iarea->imageArea->setPreviewHandler (previewHandler);
    iarea->imageArea->setImProcCoordinator (ipc);
    navigator->previewWindow->setPreviewHandler (previewHandler);
    navigator->previewWindow->setImageArea (iarea->imageArea);

    rtengine::ImageSource* is=isrc->getImageSource();
    is->setProgressListener( this );

    // try to load the last saved parameters from the cache or from the paramfile file
    ProcParams* ldprof = openThm->createProcParamsForUpdate();  // will be freed by initProfile

    // initialize profile
    Glib::ustring defProf = openThm->getType()==FT_Raw ? options.defProfRaw : options.defProfImg;
    profilep->initProfile (defProf, ldprof, NULL);

    openThm->addThumbnailListener (this);
    info_toggled ();
    
    if (beforeIarea)
    {
        beforeAfterToggled();
        beforeAfterToggled();
    }

    // If in single tab mode, the main crop window is not constructed the very first time
    // since there was no resize event
    if (iarea->imageArea->mainCropWindow)
    {
    	iarea->imageArea->mainCropWindow->cropHandler.newImage(ipc);
        iarea->imageArea->mainCropWindow->initialImageArrived();

		// In single tab mode, the image is not always updated between switches
		// normal redraw don't work, so this is the hard way
		if (!options.tabbedUI) iarea->imageArea->mainCropWindow->cropHandler.update();
    } else {
        Gtk::Allocation alloc;
        iarea->imageArea->on_resized(alloc);
    }
}

void EditorPanel::close () {
    if (ipc)
    {
        saveProfile ();
        // close image processor and the current thumbnail
        tpc->closeImage ();    // this call stops image processing
        tpc->writeOptions ();
        rtengine::ImageSource* is=isrc->getImageSource();
        is->setProgressListener( NULL );

        if (ipc)
            ipc->setPreviewImageListener (NULL);

        if (beforeIpc)
            beforeIpc->setPreviewImageListener (NULL);

        delete previewHandler;
        previewHandler= NULL;

        if(iarea)
        {
            iarea->imageArea->setPreviewHandler (NULL);
            iarea->imageArea->setImProcCoordinator (NULL);
        }
        rtengine::StagedImageProcessor::destroy (ipc);
        ipc = NULL;
        navigator->previewWindow->setPreviewHandler (NULL);

        // If the file was deleted somewhere, the openThm.descreaseRef delete the object, but we don't know here
        if (safe_file_test(fname, Glib::FILE_TEST_EXISTS)) {
        openThm->removeThumbnailListener (this);
        openThm->decreaseRef ();
        }
    }
}

void EditorPanel::saveProfile () {
    if (!ipc || !openThm) return;

    // If the file was deleted, do not generate ghost entries
    if (safe_file_test(fname, Glib::FILE_TEST_EXISTS)) {
    ProcParams params;
    ipc->getParams (&params);

        // Will call updateCache, which will update both the cached and sidecar files if necessary
        openThm->setProcParams (params, EDITOR);
}
}

Glib::ustring EditorPanel::getShortName () {

    return Glib::path_get_basename (openThm->getFileName ());
}

Glib::ustring EditorPanel::getFileName () {

    return openThm->getFileName ();
}

// TODO!!!
void EditorPanel::procParamsChanged (rtengine::procparams::ProcParams* params, rtengine::ProcEvent ev, Glib::ustring descr, ParamsEdited* paramsEdited) {

//    if (ev!=EvPhotoLoaded)
//        saveLabel->set_markup (Glib::ustring("<span foreground=\"#AA0000\" weight=\"bold\">") + M("MAIN_BUTTON_SAVE") + "</span>");
}

struct spsparams {
    bool inProcessing;
    EditorPanelIdleHelper* epih;
};

int setProgressStateUIThread (void* data) {

    spsparams* p = (spsparams*)data;

    if (p->epih->destroyed) {
        if (p->epih->pending == 1)
            delete p->epih;
        else
            p->epih->pending--;
        delete p;

        return 0;
    }

    p->epih->epanel->refreshProcessingState (p->inProcessing);
    p->epih->pending--;
    delete p;

    return 0;
}

void EditorPanel::setProgressState (bool inProcessing) {

    epih->pending++;

    spsparams* p = new spsparams;
    p->inProcessing = inProcessing;
    p->epih = epih;
    g_idle_add (setProgressStateUIThread, p);
}

struct spparams {
    double val;
    Glib::ustring str;
    Gtk::ProgressBar *pProgress;
};

int setprogressStrUI( void *p )
{
	spparams *s= (spparams*)p;

	if( ! s->str.empty() )
	   s->pProgress->set_text( M(s->str) );
	if( s->val >=0 ){
	   s->pProgress->set_fraction( s->val );
		if( s->val <1.0 )
			s->pProgress->modify_bg( Gtk::STATE_NORMAL,Gdk::Color("red") );
		else
			s->pProgress->modify_bg( Gtk::STATE_NORMAL,Gdk::Color("grey") );
	}

	delete s;
	return 0;
}

void EditorPanel::setProgress (double p)
{
	spparams *s=new spparams;
	s->val = p;
	s->pProgress = progressLabel;
	g_idle_add (setprogressStrUI, s);
}

void EditorPanel::setProgressStr (Glib::ustring str)
{
	spparams *s=new spparams;
	s->str = str;
	s->val = -1;
	s->pProgress = progressLabel;
	g_idle_add (setprogressStrUI, s);
}

// This is only called from the ThreadUI, so within the gtk thread
void EditorPanel::refreshProcessingState (bool inProcessingP) {
	spparams *s=new spparams;
    s->pProgress = progressLabel;

    if (inProcessingP) {
        if (processingStartedTime==0) processingStartedTime = ::time(NULL);

    	s->str = "PROGRESSBAR_PROCESSING";
    	s->val = 0.0;
    } else {
    // Set proc params of thumbnail. It saves it into the cache and updates the file browser.
			if (ipc && openThm && tpc->getChangedState()) {
        rtengine::procparams::ProcParams pparams;
        ipc->getParams (&pparams);
        openThm->setProcParams (pparams, EDITOR, false);
    }

		// Ring a sound if it was a long event
        if (processingStartedTime!=0) {
            time_t curTime= ::time(NULL);
            if (::difftime(curTime, processingStartedTime) > options.sndLngEditProcDoneSecs) 
                SoundManager::playSoundAsync(options.sndLngEditProcDone);

            processingStartedTime = 0;
        }

		// Set progress bar "done"
    	s->str = "PROGRESSBAR_READY";
    	s->val = 1.0;

#ifdef WIN32
        if (!firstProcessingDone && (RTWindow*)parent->getIsFullscreen()) { parent->fullscreen(); }
#endif
        firstProcessingDone = true;
}

    isProcessing=inProcessingP;

	setprogressStrUI(s);
}

struct errparams {
    Glib::ustring descr;
    EditorPanelIdleHelper* epih;
};

void EditorPanel::displayError (Glib::ustring descr) {

    if (parent) {
        Gtk::MessageDialog* msgd = new Gtk::MessageDialog (*parent, descr, true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
        msgd->set_title (M("MAIN_MSG_CANNOTSAVE"));
        msgd->run ();
        delete msgd;
    }
}

int disperrorUI (void* data) {
    errparams* p = (errparams*)data;

    if (p->epih->destroyed) {
        if (p->epih->pending == 1)
            delete p->epih;
        else
            p->epih->pending--;
        delete p;

        return 0;
    }

    p->epih->epanel->displayError (p->descr);
    p->epih->pending--;
    delete p;

    return 0;
}

void EditorPanel::error (Glib::ustring descr) {

    epih->pending++;
    errparams* p = new errparams;
    p->descr = descr;
    p->epih = epih;
    g_idle_add (disperrorUI, p);
}

void EditorPanel::info_toggled () {

    Glib::ustring infoString;
    if (!ipc || !openThm) return;
    const rtengine::ImageMetaData* idata = ipc->getInitialImage()->getMetaData();
    if (idata && idata->hasExif())
//        infoString = Glib::ustring::compose ("%1 %2\nF/%3 %4 sec\n%5: %6\n%7: %8 mm\n",
//            Glib::ustring(idata->getMake()), Glib::ustring(idata->getModel()),
//            Glib::ustring(idata->apertureToString(idata->getFNumber())), Glib::ustring(idata->shutterToString(idata->getShutterSpeed())),
//            M("QINFO_ISO"), idata->getISOSpeed(),
//            M("QINFO_FOCALLENGTH"), idata->getFocalLen())
//            + Glib::ustring::compose ("%1: %2", M("QINFO_LENS"), Glib::ustring(idata->getLens()));
infoString = Glib::ustring::compose (
            "%1 + %2\n<span size=\"small\">f/</span><span size=\"large\">%3</span>  <span size=\"large\">%4</span><span size=\"small\">s</span>  <span size=\"small\">%5</span><span size=\"large\">%6</span>  <span size=\"large\">%7</span><span size=\"small\">mm</span>\n<span size=\"small\">%8</span><span>%9</span>",
            Glib::ustring(idata->getMake()+" "+idata->getModel()),
            Glib::ustring(idata->getLens()),
            Glib::ustring(idata->apertureToString(idata->getFNumber())),
            Glib::ustring(idata->shutterToString(idata->getShutterSpeed())),
            M("QINFO_ISO"), idata->getISOSpeed(),
            idata->getFocalLen(),
            Glib::path_get_dirname(openThm->getFileName()) + G_DIR_SEPARATOR_S,
            Glib::path_get_basename(openThm->getFileName())
            );
    else
        infoString = M("QINFO_NOEXIF");

    iarea->imageArea->setInfoText (infoString);
    iarea->imageArea->infoEnabled (info->get_active ());
}

void EditorPanel::hideHistoryActivated () {

    removeIfThere (hpanedl, leftbox, false);
    if (hidehp->get_active())
        hpanedl->pack1 (*leftbox, false, true);
    options.showHistory = hidehp->get_active();

    if (options.showHistory){
    	hidehp->set_image (*iHistoryHide);
    }
    else {
    	hidehp->set_image (*iHistoryShow);
    }

    tbShowHideSidePanels_managestate();
}


void EditorPanel::tbRightPanel_1_toggled () {
/*
    removeIfThere (hpanedr, vboxright, false);
    if (tbRightPanel_1->get_active()){
        hpanedr->pack2(*vboxright, false, true);
    	tbRightPanel_1->set_image (*iRightPanel_1_Hide);
    }
    else {
    	tbRightPanel_1->set_image (*iRightPanel_1_Show);
    }
    tbShowHideSidePanels_managestate();
    */
	if (vboxright){
		if (tbRightPanel_1->get_active()){
			vboxright->show();
			tbRightPanel_1->set_image (*iRightPanel_1_Hide);
		}
		else{
			vboxright->hide();
			tbRightPanel_1->set_image (*iRightPanel_1_Show);
		}
		tbShowHideSidePanels_managestate();
	}
}

void EditorPanel::tbTopPanel_1_visible (bool visible){
	if (visible)
		tbTopPanel_1->show();
	else
		tbTopPanel_1->hide();
}

void EditorPanel::tbTopPanel_1_toggled () {

	if (catalogPane){ // catalogPane does not exist in multitab mode
		tbTopPanel_1_Active = tbTopPanel_1->get_active();

		if (tbTopPanel_1_Active){
			catalogPane->show();
			tbTopPanel_1->set_image (*iTopPanel_1_Hide);
		}
		else {
			catalogPane->hide();
			tbTopPanel_1->set_image (*iTopPanel_1_Show);
		}

		tbShowHideSidePanels_managestate();
	}
}

bool EditorPanel::handleShortcutKey (GdkEventKey* event) {

    bool ctrl = event->state & GDK_CONTROL_MASK;
    bool shift = event->state & GDK_SHIFT_MASK;
    bool alt = event->state & GDK_MOD1_MASK;


    // Editor Layout
    switch(event->keyval) {
        case GDK_L:
		    tbTopPanel_1->set_active (!tbTopPanel_1->get_active()); // toggle top panel
		    if (ctrl) hidehp->set_active (!hidehp->get_active()); // toggle History (left panel)
			if (alt) tbRightPanel_1->set_active (!tbRightPanel_1->get_active()); // toggle right panel
			return true;
    	case GDK_l:
    		if (!shift && !alt /*&& !ctrl*/){
    			hidehp->set_active (!hidehp->get_active()); // toggle History (left panel)
    		    return true;
    		}
			if (alt && !ctrl){ // toggle right panel
				tbRightPanel_1->set_active (!tbRightPanel_1->get_active());
				return true;
			}
			if (alt && ctrl){ // toggle left and right panels
				hidehp->set_active (!hidehp->get_active());
				tbRightPanel_1->set_active (!tbRightPanel_1->get_active());
				return true;
			}
        case GDK_m: // Maximize preview panel: hide top AND right AND history panels
        	if (!ctrl && !alt) {
        		toggleSidePanels();
        		return true;
        	}
        case GDK_M: // Maximize preview panel: hide top AND right AND history panels AND (fit image preview)
        	if (!ctrl && !alt) {
        		toggleSidePanelsZoomFit();
        		return true;
        	}
    }

    if (!alt){
		if (!ctrl) {
			// Normal
			switch(event->keyval) {
				case GDK_bracketright:
					tpc->coarse->rotateRight();
					return true;
				case GDK_bracketleft:
					tpc->coarse->rotateLeft();
					return true;

				case GDK_i:
				case GDK_I:
					info->set_active (!info->get_active());
					return true;
				case GDK_b:
				case GDK_B:
					beforeAfter->set_active (!beforeAfter->get_active());
					return true;
				case GDK_plus:
				case GDK_equal:
					iarea->imageArea->zoomPanel->zoomInClicked();
					return true;
				case GDK_minus:
				case GDK_underscore:
					iarea->imageArea->zoomPanel->zoomOutClicked();
					return true;
				case GDK_1:
					iarea->imageArea->zoomPanel->zoom11Clicked();
					return true;

				case GDK_f:
					iarea->imageArea->zoomPanel->zoomFitClicked();
					return true;
				case GDK_less:
					iarea->imageArea->indClippedPanel->toggleClipped(true);
					return true;
				case GDK_greater:
					iarea->imageArea->indClippedPanel->toggleClipped(false);
					return true;

				case GDK_F5:
					openThm->openDefaultViewer(event->state & GDK_SHIFT_MASK ? 2 : 1);
					return true;
			}
		}
		else {
			// With control
			switch (event->keyval) {
				case GDK_s:
					saveAsPressed();
					return true;
				case GDK_q:
					queueImgPressed();
					return true;
				case GDK_e:
					sendToGimpPressed();
					return true;
				case GDK_z:
					history->undo ();
					return true;
				case GDK_Z:
					history->redo ();
					return true;
				case GDK_F5:
					openThm->openDefaultViewer(3);
					return true;
			}
		} //if (!ctrl)
    } //if (!alt)

    if(tpc->getToolBar()->handleShortcutKey(event))
        return true;
    if(tpc->handleShortcutKey(event))
        return true;

    return false;
}

void EditorPanel::procParamsChanged (Thumbnail* thm, int whoChangedIt) {

    if (whoChangedIt!=EDITOR)
      tpc->profileChange (&openThm->getProcParams(), rtengine::EvProfileChangeNotification, M("PROGRESSDLG_PROFILECHANGEDINBROWSER"));
}

bool EditorPanel::idle_saveImage (ProgressConnector<rtengine::IImage16*> *pc, Glib::ustring fname, SaveFormat sf) {
	rtengine::IImage16* img = pc->returnValue();
	delete pc;
	if( img ) {
        setProgressStr(M("GENERAL_SAVE")); setProgress(0.9f);

        ProgressConnector<int> *ld = new ProgressConnector<int>();
        img->setSaveProgressListener (parent->getProgressListener());
        if (sf.format=="tif")
            ld->startFunc (sigc::bind(sigc::mem_fun(img, &rtengine::IImage16::saveAsTIFF), fname, sf.tiffBits, sf.tiffUncompressed),
    			            sigc::bind(sigc::mem_fun(*this,&EditorPanel::idle_imageSaved), ld, img, fname, sf));
        else if (sf.format=="png")
            ld->startFunc (sigc::bind(sigc::mem_fun(img, &rtengine::IImage16::saveAsPNG), fname, sf.pngCompression, sf.pngBits),
    			            sigc::bind(sigc::mem_fun(*this,&EditorPanel::idle_imageSaved), ld, img, fname, sf));
        else if (sf.format=="jpg")
            ld->startFunc (sigc::bind(sigc::mem_fun(img, &rtengine::IImage16::saveAsJPEG), fname, sf.jpegQuality),
    			            sigc::bind(sigc::mem_fun(*this,&EditorPanel::idle_imageSaved), ld, img, fname, sf));
    } else {
		Glib::ustring msg_ = Glib::ustring("<b>") + fname + ": Error during image processing\n</b>";
		Gtk::MessageDialog msgd (*parent, msg_, true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
		msgd.run ();

        saveimgas->set_sensitive(true);
        sendtogimp->set_sensitive(true);

	}
	rtengine::ImageSource* imgsrc = isrc->getImageSource ();
    imgsrc->setProgressListener(this);
	return false;
}

bool EditorPanel::idle_imageSaved(ProgressConnector<int> *pc,rtengine::IImage16* img,Glib::ustring fname, SaveFormat sf){
	img->free ();

	if (! pc->returnValue() ) {
		openThm->imageDeveloped ();
		// save processing parameters, if needed
		if (sf.saveParams) {
			rtengine::procparams::ProcParams pparams;
			ipc->getParams (&pparams);
			// We keep the extension to avoid overwriting the profile when we have
			// the same output filename with different extension
			//pparams.save (removeExtension (fname) + ".out" + paramFileExtension);
			pparams.save (fname + ".out" + paramFileExtension);
		}
	} else {
		Glib::ustring msg_ = Glib::ustring("<b>") + fname + ": Error during image saving\n</b>";
		Gtk::MessageDialog msgd (*parent, msg_, true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
		msgd.run ();
    }

    saveimgas->set_sensitive(true);
    sendtogimp->set_sensitive(true);

	parent->setProgressStr("");
	parent->setProgress(0.);

    setProgressState(false);

	delete pc;
    return false;
}

BatchQueueEntry* EditorPanel::createBatchQueueEntry () {

    rtengine::procparams::ProcParams pparams;
    ipc->getParams (&pparams);
    //rtengine::ProcessingJob* job = rtengine::ProcessingJob::create (ipc->getInitialImage(), pparams);
    rtengine::ProcessingJob* job = rtengine::ProcessingJob::create (openThm->getFileName (), openThm->getType()==FT_Raw, pparams);
    int prevh = options.maxThumbnailHeight;
    int prevw = prevh;
    guint8* prev = NULL;//(guint8*) previewHandler->getImagePreview (prevw, prevh);
    double tmpscale;
    rtengine::IImage8* img = openThm->processThumbImage (pparams, options.maxThumbnailHeight, tmpscale);
    if (img) {
    	prevw = img->getWidth ();
    	prevh = img->getHeight ();
        prev = new guint8 [prevw*prevh*3];
        memcpy (prev, img->getData (), prevw*prevh*3);
        img->free();
    }
    return new BatchQueueEntry (job, pparams, openThm->getFileName(), prev, prevw, prevh, openThm);
}



void EditorPanel::saveAsPressed () {
	if (!ipc || !openThm) return;
	bool fnameOK = false;
	Glib::ustring fname;

	saveAsDialog->setInitialFileName (removeExtension (Glib::path_get_basename (openThm->getFileName())));
	do {
		saveAsDialog->run ();
		if (saveAsDialog->getResponse()==Gtk::RESPONSE_CANCEL)
			return;

		// The SaveAsDialog ensure that a filename has been specified
		fname = saveAsDialog->getFileName ();

		options.lastSaveAsPath = saveAsDialog->getDirectory ();
		options.saveAsDialogWidth = saveAsDialog->get_width();
		options.saveAsDialogHeight = saveAsDialog->get_height();

		SaveFormat sf = saveAsDialog->getFormat ();

		options.saveFormat = sf;
		options.autoSuffix = saveAsDialog->getAutoSuffix();

		if (saveAsDialog->getImmediately ()) {
			// separate filename and the path to the destination directory
			Glib::ustring dstdir = Glib::path_get_dirname (fname);
			Glib::ustring dstfname = Glib::path_get_basename (removeExtension(fname));

			if (saveAsDialog->getAutoSuffix()) {

				Glib::ustring fnameTemp;
				for (int tries=0; tries<100; tries++) {
					if (tries==0)
						fnameTemp = Glib::ustring::compose ("%1.%2", Glib::build_filename (dstdir,  dstfname), sf.format);
					else
						fnameTemp = Glib::ustring::compose ("%1-%2.%3", Glib::build_filename (dstdir,  dstfname), tries, sf.format);

					if (!safe_file_test (fnameTemp, Glib::FILE_TEST_EXISTS)) {
						fname = fnameTemp;
						fnameOK = true;
						break;
					}
				}
			}
			// check if it exists
			if (!fnameOK) {
				if (safe_file_test (fname, Glib::FILE_TEST_EXISTS)) {
					Glib::ustring msg_ = Glib::ustring("<b>") + fname + ": " + M("MAIN_MSG_ALREADYEXISTS") + "\n" + M("MAIN_MSG_QOVERWRITE") + "</b>";
					Gtk::MessageDialog msgd (*parent, msg_, true, Gtk::MESSAGE_WARNING, Gtk::BUTTONS_YES_NO, true);
					int response = msgd.run ();
					fnameOK = (response==Gtk::RESPONSE_YES);
				}
				else fnameOK = true;
			}

			if (fnameOK) {
				// save image
				rtengine::procparams::ProcParams pparams;
				ipc->getParams (&pparams);
				rtengine::ProcessingJob* job = rtengine::ProcessingJob::create (ipc->getInitialImage(), pparams);

				ProgressConnector<rtengine::IImage16*> *ld = new ProgressConnector<rtengine::IImage16*>();
				ld->startFunc(sigc::bind(sigc::ptr_fun(&rtengine::processImage), job, err, parent->getProgressListener(), options.tunnelMetaData ),
							  sigc::bind(sigc::mem_fun( *this,&EditorPanel::idle_saveImage ),ld,fname,sf ));
				saveimgas->set_sensitive(false);
				sendtogimp->set_sensitive(false);
			}
		}
		else {
			BatchQueueEntry* bqe = createBatchQueueEntry ();
			bqe->outFileName = fname;
			bqe->saveFormat = saveAsDialog->getFormat ();
			parent->addBatchQueueJob (bqe, saveAsDialog->getToHeadOfQueue ());
			fnameOK = true;
		}
		// ask parent to redraw file browser
		// ... or does it automatically when the tab is switched to it
	} while (!fnameOK);
}

void EditorPanel::queueImgPressed () {
	if (!ipc || !openThm) return;
    saveProfile ();
    parent->addBatchQueueJob (createBatchQueueEntry ());
}

void EditorPanel::sendToGimpPressed () {
	if (!ipc || !openThm) return;
    // develop image
    rtengine::procparams::ProcParams pparams;
    ipc->getParams (&pparams);
    rtengine::ProcessingJob* job = rtengine::ProcessingJob::create (ipc->getInitialImage(), pparams);
    ProgressConnector<rtengine::IImage16*> *ld = new ProgressConnector<rtengine::IImage16*>();
    ld->startFunc(sigc::bind(sigc::ptr_fun(&rtengine::processImage), job, err, parent->getProgressListener(), options.tunnelMetaData ),
    		      sigc::bind(sigc::mem_fun( *this,&EditorPanel::idle_sendToGimp ),ld ));
    saveimgas->set_sensitive(false);
    sendtogimp->set_sensitive(false);
}

bool EditorPanel::idle_sendToGimp( ProgressConnector<rtengine::IImage16*> *pc){

	rtengine::IImage16* img = pc->returnValue();
	delete pc;
    if (img) {
        // get file name base
        Glib::ustring shortname = removeExtension (Glib::path_get_basename (openThm->getFileName()));
        Glib::ustring dirname = Glib::get_tmp_dir ();
        Glib::ustring fname = Glib::build_filename (dirname, shortname);

        SaveFormat sf;
        sf.format = "tif";
        sf.tiffBits = 16;
        sf.tiffUncompressed = true;
        sf.saveParams = true;

        Glib::ustring fileName = Glib::ustring::compose ("%1.%2", fname, sf.format);

				int tries = 1;
				while (safe_file_test (fileName, Glib::FILE_TEST_EXISTS) && tries<1000) {
					fileName = Glib::ustring::compose("%1-%2.%3", fname, tries, sf.format);
					tries++;
				}
				if (tries==1000){
					img->free ();
					return false;
				}

        ProgressConnector<int> *ld = new ProgressConnector<int>();
        img->setSaveProgressListener (parent->getProgressListener());
       	ld->startFunc (sigc::bind(sigc::mem_fun(img, &rtengine::IImage16::saveAsTIFF), fileName, sf.tiffBits, sf.tiffUncompressed),
        			   sigc::bind(sigc::mem_fun(*this,&EditorPanel::idle_sentToGimp), ld, img, fileName));
    }else{
				Glib::ustring msg_ = Glib::ustring("<b> Error during image processing\n</b>");
				Gtk::MessageDialog msgd (*parent, msg_, true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
				msgd.run ();
        saveimgas->set_sensitive(true);
        sendtogimp->set_sensitive(true);
    }
    return false;
}

bool EditorPanel::idle_sentToGimp(ProgressConnector<int> *pc,rtengine::IImage16* img,Glib::ustring filename){
    img->free ();
    int errore = pc->returnValue();
    delete pc;
    if (!errore) {
                        saveimgas->set_sensitive(true);
                        sendtogimp->set_sensitive(true);
    	                parent->setProgressStr("");
    	                parent->setProgress(0.);
						bool success=false;
						Glib::ustring cmdLine;
						Glib::ustring executable;
						// start gimp
						if (options.editorToSendTo==1) {
#ifdef WIN32
								executable = Glib::build_filename (Glib::build_filename(options.gimpDir,"bin"), "gimp-win-remote");
								cmdLine = Glib::ustring("\"") + executable + Glib::ustring("\" gimp-2.4.exe ") + Glib::ustring("\"") + filename + Glib::ustring("\"");
								if ( safe_file_test(executable, (Glib::FILE_TEST_EXISTS|Glib::FILE_TEST_IS_EXECUTABLE)) ) {
									success = safe_spawn_command_line_async (cmdLine);
								}
#else
								cmdLine = Glib::ustring("gimp-remote ") + Glib::ustring(" \"") + filename + Glib::ustring("\"");
								success = safe_spawn_command_line_async (cmdLine);
#endif
								if (!success){
#ifdef WIN32
										int ver = 12;
										while (!success && ver) {
												executable = Glib::build_filename (Glib::build_filename(options.gimpDir,"bin"), Glib::ustring::compose(Glib::ustring("gimp-2.%1.exe"),ver));
												if ( safe_file_test(executable, (Glib::FILE_TEST_EXISTS|Glib::FILE_TEST_IS_EXECUTABLE)) ) {
													cmdLine = Glib::ustring("\"") + executable + Glib::ustring("\" \"") + filename + Glib::ustring("\"");
													success = safe_spawn_command_line_async (cmdLine);
												}
												ver--;
										}
#elif defined __APPLE__
										cmdLine = Glib::ustring("gimp \"") + filename + Glib::ustring("\"");
										success = safe_spawn_command_line_async (cmdLine);
#else
										cmdLine = Glib::ustring("gimp \"") + filename + Glib::ustring("\"");
										success = safe_spawn_command_line_async (cmdLine);
#endif
								}
						}
						else if (options.editorToSendTo==2) {
#ifdef WIN32
								executable = Glib::build_filename(options.psDir,"Photoshop.exe");
								if ( safe_file_test(executable, (Glib::FILE_TEST_EXISTS|Glib::FILE_TEST_IS_EXECUTABLE)) ) {
										cmdLine = Glib::ustring("\"") + executable + Glib::ustring("\" \"") + filename + Glib::ustring("\"");
										success = safe_spawn_command_line_async (cmdLine);
								}
#else
		#ifdef __APPLE__
								cmdLine = Glib::ustring("open -a \'") + Glib::build_filename(options.psDir,"Photoshop.app\' ")  + Glib::ustring("\'") + filename + Glib::ustring("\'");
		#else
								cmdLine = Glib::ustring("\"") + Glib::build_filename(options.psDir,"Photoshop.exe") + Glib::ustring("\" \"") + filename + Glib::ustring("\"");
		#endif
								success = safe_spawn_command_line_async (cmdLine);
#endif
						}
						else if (options.editorToSendTo==3) {
#ifdef WIN32
								if ( safe_file_test(options.customEditorProg, (Glib::FILE_TEST_EXISTS|Glib::FILE_TEST_IS_EXECUTABLE)) ) {
										cmdLine = Glib::ustring("\"") + options.customEditorProg + Glib::ustring("\" \"") + filename + Glib::ustring("\"");
										success = safe_spawn_command_line_async (cmdLine);
								}
#else
		#ifdef __APPLE__
								cmdLine = options.customEditorProg + filename;
		#else
								cmdLine = Glib::ustring("\"") + options.customEditorProg + Glib::ustring("\" \"") + filename + Glib::ustring("\"");
		#endif
								success = safe_spawn_command_line_async (cmdLine);
#endif
						}

						if (!success) {
								Gtk::MessageDialog* msgd = new Gtk::MessageDialog (*parent, M("MAIN_MSG_CANNOTSTARTEDITOR"), false, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
								msgd->set_secondary_text (M("MAIN_MSG_CANNOTSTARTEDITOR_SECONDARY"));
								msgd->set_title (M("MAIN_BUTTON_SENDTOEDITOR"));
								msgd->run ();
								delete msgd;
						}

    }

    return false;
}

void EditorPanel::historyBeforeLineChanged (const rtengine::procparams::ProcParams& params) {

    if (beforeIpc) {
        ProcParams* pparams = beforeIpc->getParamsForUpdate (rtengine::EvProfileChanged);
        *pparams = params;
        beforeIpc->paramsUpdateReady ();  // starts the IPC processinp
    }
}

void EditorPanel::beforeAfterToggled () {

    if(!ipc)
        return;

    removeIfThere (beforeAfterBox,  beforeBox, false);
    removeIfThere (afterBox,  afterLabel, false);

    if (beforeIarea) {
        if (beforeIpc)
            beforeIpc->stopProcessing ();
        iarea->setBeforeAfterViews (NULL, iarea);
        delete beforeIarea;
        beforeIarea = NULL;
        if (beforeIpc)
            beforeIpc->setPreviewImageListener (NULL);
        delete beforePreviewHandler;
        beforePreviewHandler = NULL;
        if (beforeIpc)
            rtengine::StagedImageProcessor::destroy (beforeIpc);
        beforeIpc = NULL;
    }

    if (beforeAfter->get_active ()) {
        int errorCode=0;
        rtengine::InitialImage *beforeImg = rtengine::InitialImage::load ( isrc->getImageSource ()->getFileName(),  openThm->getType()==FT_Raw , &errorCode, NULL);
        if( !beforeImg || errorCode )
        	return;

        beforeIarea = new ImageAreaPanel ();

        beforeLabel = Gtk::manage (new Gtk::Label ());
        beforeLabel->set_markup (Glib::ustring("<b>") + M("GENERAL_BEFORE") + "</b>");
        beforeBox = Gtk::manage (new Gtk::VBox ());
        beforeBox->pack_start (*beforeLabel, Gtk::PACK_SHRINK, 2);
        beforeBox->pack_start (*beforeIarea);

        afterLabel = Gtk::manage (new Gtk::Label ());
        afterLabel->set_markup (Glib::ustring("<b>") + M("GENERAL_AFTER") + "</b>");
        afterBox->pack_start (*afterLabel, Gtk::PACK_SHRINK, 2);
        afterBox->reorder_child (*afterLabel, 0);

        beforeAfterBox->pack_start (*beforeBox);
        beforeAfterBox->reorder_child (*beforeBox, 0);
        beforeAfterBox->show_all ();

        beforePreviewHandler = new PreviewHandler ();

        beforeIpc = rtengine::StagedImageProcessor::create (beforeImg);
        beforeIpc->setPreviewScale (10);
        beforeIpc->setPreviewImageListener (beforePreviewHandler);
        beforeIarea->imageArea->setPreviewHandler (beforePreviewHandler);
        beforeIarea->imageArea->setImProcCoordinator (beforeIpc);

        iarea->setBeforeAfterViews (beforeIarea, iarea);
        beforeIarea->setBeforeAfterViews (beforeIarea, iarea);

        rtengine::procparams::ProcParams params;
        if (history->getBeforeLineParams (params))
            historyBeforeLineChanged (params);
    }
}

void EditorPanel::histogramChanged (LUTu & histRed, LUTu & histGreen, LUTu & histBlue, LUTu & histLuma, LUTu & histToneCurve, LUTu & histLCurve,
    LUTu & histRedRaw, LUTu & histGreenRaw, LUTu & histBlueRaw) {

    if (histogramPanel) histogramPanel->histogramChanged (histRed, histGreen, histBlue, histLuma, histRedRaw, histGreenRaw, histBlueRaw);
    tpc->updateCurveBackgroundHistogram (histToneCurve, histLCurve);
}

bool EditorPanel::CheckSidePanelsVisibility(){
	if(tbTopPanel_1->get_active()==false && tbRightPanel_1->get_active()==false && hidehp->get_active()==false)
		return false;
	else
		return true;
}
void EditorPanel::toggleSidePanels(){
	// Maximize preview panel:
	// toggle top AND right AND history panels

	bool bAllSidePanelsVisible;
	bAllSidePanelsVisible= CheckSidePanelsVisibility();

	tbTopPanel_1->set_active (!bAllSidePanelsVisible);
	tbRightPanel_1->set_active (!bAllSidePanelsVisible);
	hidehp->set_active (!bAllSidePanelsVisible);
}

void EditorPanel::toggleSidePanelsZoomFit(){
	toggleSidePanels();

	// fit image preview
	// !!! TODO this does not want to work... seems to have an effect on a subsequent key press
	// iarea->imageArea->zoomPanel->zoomFitClicked();
}

void EditorPanel::tbShowHideSidePanels_managestate(){
	bool bAllSidePanelsVisible;
	bAllSidePanelsVisible = CheckSidePanelsVisibility();
	ShowHideSidePanelsconn.block (true);

	tbShowHideSidePanels->set_active (!bAllSidePanelsVisible);

	ShowHideSidePanelsconn.block (false);
}
