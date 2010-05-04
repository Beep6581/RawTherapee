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
#include <editorpanel.h>
#include <options.h>
#include <progressdialog.h>
#include <rtwindow.h>
#include <guiutils.h>
#include <procparamchangers.h>
#include <safegtk.h>

using namespace rtengine::procparams;

EditorPanel::EditorPanel () : parent(NULL), beforeIarea(NULL), beforePreviewHandler(NULL), beforeIpc(NULL) {

    epih = new EditorPanelIdleHelper;
    epih->epanel = this;
    epih->destroyed = false;
    epih->pending = 0;

// construct toolpanelcoordinator
    tpc = new ToolPanelCoordinator ();

// build GUI
    // build left side panel
    leftbox = Gtk::manage (new Gtk::VBox ());
    leftbox->set_border_width (4);

    histogramPanel = Gtk::manage (new HistogramPanel ());
    histogramPanel->set_size_request (-1, 150);
//    leftbox->pack_start (*histogramPanel, Gtk::PACK_SHRINK, 4);

    profilep = Gtk::manage (new ProfilePanel ());
    Gtk::Frame* ppframe = Gtk::manage (new Gtk::Frame ());
    ppframe->add (*profilep);
    ppframe->set_label (M("PROFILEPANEL_LABEL"));
//    leftbox->pack_start (*ppframe, Gtk::PACK_SHRINK, 4);

    navigator = Gtk::manage (new Navigator ());
    navigator->previewWindow->set_size_request (-1, 150);
    leftbox->pack_start (*navigator, Gtk::PACK_SHRINK, 4);

    history = Gtk::manage (new History ());
    leftbox->pack_start (*history);

    leftbox->show_all ();

    // build the middle of the screen
    Gtk::VBox* editbox = Gtk::manage (new Gtk::VBox ());

    info = Gtk::manage (new Gtk::ToggleButton ());
    Gtk::Image* infoimg = Gtk::manage (new Gtk::Image (argv0+"/images/info.png"));
    info->add (*infoimg);
    info->set_relief(Gtk::RELIEF_NONE);
    info->set_tooltip_text (M("MAIN_TOOLTIP_QINFO"));

    beforeAfter = Gtk::manage (new Gtk::ToggleButton (M("MAIN_TOGGLE_BEFORE_AFTER")));
    beforeAfter->set_tooltip_text (M("MAIN_TOOLTIP_TOGGLE"));


    Gtk::VSeparator* vsept = Gtk::manage (new Gtk::VSeparator ());
    Gtk::VSeparator* vsepz = Gtk::manage (new Gtk::VSeparator ());
    Gtk::VSeparator* vsepi = Gtk::manage (new Gtk::VSeparator ());
    Gtk::VSeparator* vseph = Gtk::manage (new Gtk::VSeparator ());

    hidehp = Gtk::manage (new Gtk::ToggleButton ());
    Gtk::Label* hidehpLabel = Gtk::manage (new Gtk::Label ());
    hidehpLabel->set_markup ("<b>H</b>");
    Gtk::Image* hpimg = Gtk::manage (new Gtk::Image (argv0+"/images/left.png"));
    Gtk::HBox* hidehpBox = Gtk::manage (new Gtk::HBox ());
    hidehpBox->pack_start (*hpimg, Gtk::PACK_SHRINK, 2);
    hidehpBox->pack_start (*hidehpLabel, Gtk::PACK_SHRINK, 2);
    hidehp->add (*hidehpBox);
    hidehp->set_relief(Gtk::RELIEF_NONE);
    hidehp->set_active (options.showHistory);
    hidehp->set_tooltip_text (M("MAIN_TOOLTIP_HIDEHP"));

    Gtk::VSeparator* vsepcl = Gtk::manage (new Gtk::VSeparator ());
    Gtk::VSeparator* vsepz2 = Gtk::manage (new Gtk::VSeparator ());

    iarea = new ImageAreaPanel ();

    Gtk::HBox* toolBarPanel = Gtk::manage (new Gtk::HBox ());
    toolBarPanel->pack_start (*hidehp, Gtk::PACK_SHRINK, 1);
    toolBarPanel->pack_start (*vseph, Gtk::PACK_SHRINK, 2);
    toolBarPanel->pack_start (*info, Gtk::PACK_SHRINK, 1);
    toolBarPanel->pack_start (*beforeAfter, Gtk::PACK_SHRINK, 1);
    toolBarPanel->pack_start (*vsepi, Gtk::PACK_SHRINK, 2);
    toolBarPanel->pack_start (*tpc->getToolBar(), Gtk::PACK_SHRINK, 1);
    toolBarPanel->pack_start (*vsept, Gtk::PACK_SHRINK, 2);
    toolBarPanel->pack_end   (*tpc->coarse, Gtk::PACK_SHRINK, 4);
    toolBarPanel->pack_end   (*vsepcl, Gtk::PACK_SHRINK, 4);
    toolBarPanel->pack_end   (*iarea->imageArea->indClippedPanel, Gtk::PACK_SHRINK, 0);
    toolBarPanel->pack_end   (*vsepz, Gtk::PACK_SHRINK, 2);
    
    afterBox = Gtk::manage (new Gtk::VBox ());
    afterBox->pack_start (*iarea);

    beforeAfterBox = Gtk::manage (new Gtk::HBox());
    beforeAfterBox->pack_start (*afterBox);
    
    editbox->pack_start (*toolBarPanel, Gtk::PACK_SHRINK);
    editbox->pack_start (*beforeAfterBox);

    // build right side panel
    vboxright = Gtk::manage (new Gtk::VBox (false, 0));
    vboxright->set_border_width (4);
    vboxright->pack_start (*histogramPanel, Gtk::PACK_SHRINK, 4);
    vboxright->pack_start (*ppframe, Gtk::PACK_SHRINK, 4);
    // main notebook
    vboxright->pack_start (*tpc->toolPanelNotebook);

    // buttons & status
    Gtk::HBox* iops = Gtk::manage (new Gtk::HBox ());
    saveimgas = Gtk::manage (new Gtk::Button (M("MAIN_BUTTON_SAVE")));
    saveimgas->set_image (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-save"), Gtk::ICON_SIZE_BUTTON)));
    queueimg = Gtk::manage (new Gtk::Button (M("MAIN_BUTTON_PUTTOQUEUE")));
    queueimg->set_image (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-execute"), Gtk::ICON_SIZE_BUTTON)));
    sendtogimp = Gtk::manage (new Gtk::Button (M("MAIN_BUTTON_SENDTOEDITOR")));
    sendtogimp->set_image (*Gtk::manage(new Gtk::Image (argv0+"/images/gimp.png")));
    iops->pack_start (*saveimgas, Gtk::PACK_SHRINK);
    iops->pack_start (*queueimg, Gtk::PACK_SHRINK);
    iops->pack_start (*sendtogimp, Gtk::PACK_SHRINK);

    statusBox = Gtk::manage (new Gtk::HBox ());
    progressLabel = Gtk::manage (new Gtk::Label(""));
    statusBox->pack_start (*progressLabel);
    red = new Gtk::Image (argv0+"/images/red.png");
    green = new Gtk::Image (argv0+"/images/green.png");
    red->show ();
    green->show ();
    statusBox->pack_end (*green, Gtk::PACK_SHRINK, 4);
    iops->pack_start(*statusBox, Gtk::PACK_SHRINK, 4);

    iops->pack_end (*iarea->imageArea->zoomPanel, Gtk::PACK_SHRINK, 1);
    iops->pack_end (*vsepz2, Gtk::PACK_SHRINK, 2);

    
    editbox->pack_start (*Gtk::manage(new Gtk::HSeparator()), Gtk::PACK_SHRINK, 4);
    editbox->pack_start (*iops, Gtk::PACK_SHRINK, 4);
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

    Gtk::Frame* vbfr = Gtk::manage (new Gtk::Frame ());
    vbfr->add (*editbox);
    hpanedl->pack2(*vbfr, true, true);

    hpanedr->pack1(*hpanedl, true, true);
    hpanedr->pack2(*vboxright, false, true);

    pack_start (*hpanedr);
    show_all ();

    // save as dialog
    if (Glib::file_test (options.lastSaveAsPath, Glib::FILE_TEST_IS_DIR)) 
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
	iarea->imageArea->setImageAreaToolListener (tpc);
    
// initialize components
    info->set_active (options.showInfo);
    tpc->readOptions ();

// connect event handlers
    info->signal_toggled().connect( sigc::mem_fun(*this, &EditorPanel::info_toggled) );
    beforeAfter->signal_toggled().connect( sigc::mem_fun(*this, &EditorPanel::beforeAfterToggled) );
    hidehp->signal_toggled().connect( sigc::mem_fun(*this, &EditorPanel::hideHistoryActivated) );
    saveimgas->signal_pressed().connect( sigc::mem_fun(*this, &EditorPanel::saveAsPressed) );
    queueimg->signal_pressed().connect( sigc::mem_fun(*this, &EditorPanel::queueImgPressed) );
    sendtogimp->signal_pressed().connect( sigc::mem_fun(*this, &EditorPanel::sendToGimpPressed) );

}

bool EditorPanel::beforeClosing () {

    options.toolPanelWidth = vboxright->get_width ();
    return true;
}

EditorPanel::~EditorPanel () {

    history->setHistoryBeforeLineListener (NULL);
    // the order is important!
    delete iarea;
    delete beforeIarea;

    if (ipc)
        ipc->setPreviewImageListener (NULL);
    if (beforeIpc)
        beforeIpc->setPreviewImageListener (NULL);

    delete previewHandler;
    delete beforePreviewHandler;

    if (ipc)
        close ();

    if (epih->pending)
        epih->destroyed = true;
    else
        delete epih;

    delete tpc;

    delete red;
    delete green;
    delete leftbox;
    delete vboxright;
    
    delete saveAsDialog;
}

void EditorPanel::on_realize () {
    
    Gtk::VBox::on_realize ();
    vboxright->set_size_request (options.toolPanelWidth, -1);
}

void EditorPanel::open (Thumbnail* tmb, rtengine::InitialImage* isrc) {

    // initialize everything
    openThm = tmb;
    openThm->increaseRef ();

    previewHandler = new PreviewHandler ();

    this->isrc = isrc;
    ipc = rtengine::StagedImageProcessor::create (isrc);
    ipc->setProgressListener (this);
    ipc->setPreviewImageListener (previewHandler);
    ipc->setPreviewScale (10);
    tpc->initImage (ipc, tmb->getType()==FT_Raw);
    ipc->setHistogramListener (this);

//    iarea->fitZoom ();   // tell to the editorPanel that the next image has to be fitted to the screen
    iarea->imageArea->setPreviewHandler (previewHandler);
    iarea->imageArea->setImProcCoordinator (ipc);
    navigator->previewWindow->setPreviewHandler (previewHandler);
    navigator->previewWindow->setImageArea (iarea->imageArea);

    // try to load the last saved parameters from the cache or from the pp2 file
    ProcParams* ldprof = NULL;
    if (openThm->hasProcParams()) {
        ldprof = new ProcParams ();
        *ldprof = openThm->getProcParams ();
    }

    // initialize profile
    if (openThm->getType()!=FT_Raw)
        profilep->initProfile (options.defProfImg, ldprof, NULL);
    else
        profilep->initProfile (options.defProfRaw, ldprof, NULL);      

    openThm->addThumbnailListener (this);
    info_toggled ();
}

void EditorPanel::close () {

    saveProfile ();

    // close image processor and the current thumbnail
    tpc->closeImage ();    // this call stops image processing
    tpc->writeOptions ();

    if (ipc)
        rtengine::StagedImageProcessor::destroy (ipc);
    if (beforeIpc)
        rtengine::StagedImageProcessor::destroy (beforeIpc);

    openThm->removeThumbnailListener (this);
    openThm->decreaseRef ();
}

void EditorPanel::saveProfile () {

    ProcParams params;
    ipc->getParams (&params);
    
    if (options.saveParamsFile)
        params.save (openThm->getFileName() + ".pp2");
    if (openThm && options.saveParamsCache) 
        openThm->setProcParams (params, EDITOR);
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
    bool state;
    EditorPanelIdleHelper* epih;
};

int setprocstate (void* data) {

    gdk_threads_enter ();
    spsparams* p = (spsparams*)data;

    if (p->epih->destroyed) {
        if (p->epih->pending == 1)
            delete p->epih;
        else    
            p->epih->pending--;
        delete p;
        gdk_threads_leave ();
        return 0;
    }

    p->epih->epanel->refreshProcessingState (p->state);
    p->epih->pending--;
    delete p;
    gdk_threads_leave ();
    return 0;
}

void EditorPanel::setProgressState (int state) {

    epih->pending++;

    spsparams* p = new spsparams;
    p->state = state;
    p->epih = epih;
    g_idle_add (setprocstate, p);
}

void EditorPanel::refreshProcessingState (bool state) {

    // Set proc params of thumbnail. It saves it into the cache and updates the file browser.
    if (ipc && openThm && !state && tpc->getChangedState()) {
        rtengine::procparams::ProcParams pparams;
        ipc->getParams (&pparams);
        openThm->setProcParams (pparams, EDITOR, false);
    }

    // change state of the led
    std::vector<Widget*> children = (std::vector<Widget*>) statusBox->get_children();
    if (children.size()>=1) {
        Gtk::Widget* wlast = children[children.size()-1];
        if (wlast)
            statusBox->remove (*wlast);
    }
    if (state) 
        statusBox->pack_end (*red, Gtk::PACK_SHRINK, 4);
    else
        statusBox->pack_end (*green, Gtk::PACK_SHRINK, 4);
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

int disperror (void* data) {

    gdk_threads_enter ();
    errparams* p = (errparams*)data;

    if (p->epih->destroyed) {
        if (p->epih->pending == 1)
            delete p->epih;
        else    
            p->epih->pending--;
        delete p;
        gdk_threads_leave ();
        return 0;
    }

    p->epih->epanel->displayError (p->descr);
    p->epih->pending--;
    delete p;
    gdk_threads_leave ();
    return 0;
}

void EditorPanel::error (Glib::ustring descr) {

    epih->pending++;
    errparams* p = new errparams;
    p->descr = descr;
    p->epih = epih;
    g_idle_add (disperror, p);
}

void EditorPanel::info_toggled () {

    Glib::ustring infoString;

    const rtengine::ImageMetaData* idata = ipc->getInitialImage()->getMetaData();
    if (idata && idata->hasExif())
        infoString = Glib::ustring::compose ("%1 %2\nF/%3 %4 sec\n%5: %6\n%7: %8 mm\n", 
            Glib::ustring(idata->getMake()), Glib::ustring(idata->getModel()),
            Glib::ustring(idata->apertureToString(idata->getFNumber())), Glib::ustring(idata->shutterToString(idata->getShutterSpeed())),
            M("QINFO_ISO"), idata->getISOSpeed(),
            M("QINFO_FOCALLENGTH"), idata->getFocalLen())
            + Glib::ustring::compose ("%1: %2", M("QINFO_LENS"), Glib::ustring(idata->getLens()));
    else 
        infoString = M("QINFO_NOEXIF");

    iarea->imageArea->setInfoText (infoString);
    iarea->imageArea->infoEnabled (info->get_active ());
}

void EditorPanel::hideHistoryActivated () {

    removeIfThere (hpanedl, leftbox, false);
    if (hidehp->get_active()) 
        hpanedl->pack1 (*leftbox, false, true);
}

bool EditorPanel::handleShortcutKey (GdkEventKey* event) {

    if (event->keyval==GDK_H || event->keyval==GDK_h) {
        hidehp->set_active (!hidehp->get_active());
        return true;
    }
    else if ((event->keyval==GDK_Z || event->keyval==GDK_z) && event->state & GDK_CONTROL_MASK && !(event->state & GDK_SHIFT_MASK)) {
        history->undo ();
        return true;
    }
    else if ((event->keyval==GDK_Z || event->keyval==GDK_z) && event->state & GDK_CONTROL_MASK && event->state & GDK_SHIFT_MASK) {
        history->redo ();
        return true;
    }
    else if (event->keyval==GDK_w || event->keyval==GDK_W) {
        tpc->getToolBar()->wb_pressed ();
        return true;
    }
    else if (event->keyval==GDK_c || event->keyval==GDK_C) {
        tpc->getToolBar()->crop_pressed ();
        return true;
    }
    else if (event->keyval==GDK_s || event->keyval==GDK_S) {
        tpc->getToolBar()->stra_pressed ();
        return true;
    }
    else if (event->keyval==GDK_n || event->keyval==GDK_N) {
        tpc->getToolBar()->hand_pressed ();
        return true;
    }
    else
        return false;
}

void EditorPanel::procParamsChanged (Thumbnail* thm, int whoChangedIt) {

    if (whoChangedIt!=EDITOR) 
      tpc->profileChange (&openThm->getProcParams(), rtengine::EvProfileChangeNotification, M("PROGRESSDLG_PROFILECHANGEDINBROWSER"));    
}

bool EditorPanel::idle_saveImage (ProgressConnector<rtengine::IImage16*> *pc, Glib::ustring fname, SaveFormat sf, bool findNewNameIfNeeded){
	rtengine::IImage16* img = pc->returnValue();
	delete pc;
	if( img )
	   saveImage( img, fname, sf, findNewNameIfNeeded);
	else{
		gdk_threads_enter ();
		Glib::ustring msg_ = Glib::ustring("<b>") + fname + ": Error during image processing\n</b>";
		Gtk::MessageDialog msgd (*parent, msg_, true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
		msgd.run ();
		gdk_threads_leave ();

        saveimgas->set_sensitive(true);
        sendtogimp->set_sensitive(true);

	}
	return false;
}

BatchQueueEntry* EditorPanel::createBatchQueueEntry () {

    rtengine::procparams::ProcParams pparams;
    ipc->getParams (&pparams);
    rtengine::ProcessingJob* job = rtengine::ProcessingJob::create (ipc->getInitialImage(), pparams);
    int prevh = options.maxThumbnailHeight;
    int prevw = prevh;
    guint8* prev = NULL;//(guint8*) previewHandler->getImagePreview (prevw, prevh);
    return new BatchQueueEntry (job, pparams, openThm->getFileName(), prev, prevw, prevh, openThm);
}

int EditorPanel::saveImage (rtengine::IImage16* img, Glib::ustring& fname, SaveFormat sf, bool findNewNameIfNeeded) {

    Glib::ustring fileName = Glib::ustring::compose ("%1.%2", fname, sf.format);
    if (findNewNameIfNeeded) {
        int tries = 1;
        while (Glib::file_test (fileName, Glib::FILE_TEST_EXISTS) && tries<1000) {
            fileName = Glib::ustring::compose("%1-%2.%3", fname, tries, sf.format);
            tries++;
        }
        if (tries==1000)
            return -1000;
    }
    ProgressConnector<int> *ld = new ProgressConnector<int>();
    img->setSaveProgressListener (parent->getProgressListener());
    if (sf.format=="tif")
    	ld->startFunc (sigc::bind(sigc::mem_fun(img, &rtengine::IImage16::saveAsTIFF), fileName, sf.tiffBits),
    			       sigc::bind(sigc::mem_fun(*this,&EditorPanel::idle_imageSaved), ld, img, fileName,sf));
    else if (sf.format=="png")
    	ld->startFunc (sigc::bind(sigc::mem_fun(img, &rtengine::IImage16::saveAsPNG), fileName, sf.pngCompression, sf.pngBits),
    			       sigc::bind(sigc::mem_fun(*this,&EditorPanel::idle_imageSaved), ld, img, fileName,sf));
    else if (sf.format=="jpg")
    	ld->startFunc (sigc::bind(sigc::mem_fun(img, &rtengine::IImage16::saveAsJPEG), fileName, sf.jpegQuality),
    			       sigc::bind(sigc::mem_fun(*this,&EditorPanel::idle_imageSaved), ld, img, fileName,sf));
    return 0;
}

bool EditorPanel::idle_imageSaved(ProgressConnector<int> *pc,rtengine::IImage16* img,Glib::ustring fname, SaveFormat sf){
	img->free ();
	if (! pc->returnValue() ) {
		openThm->imageDeveloped ();
		// save processing parameters, if needed
		if (sf.saveParams) {
			rtengine::procparams::ProcParams pparams;
			ipc->getParams (&pparams);
			pparams.save (removeExtension (fname) + ".out.pp2");
		}
	}else{
		gdk_threads_enter ();
		Glib::ustring msg_ = Glib::ustring("<b>") + fname + ": Error during image saving\n</b>";
		Gtk::MessageDialog msgd (*parent, msg_, true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
		msgd.run ();
		gdk_threads_leave ();
    }
    saveimgas->set_sensitive(true);
    sendtogimp->set_sensitive(true);
	parent->setProgressStr("");
	parent->setProgress(0.);
	delete pc;
    return false;
}

void EditorPanel::saveAsPressed () {

    // obtaining short name without extension
    saveAsDialog->setInitialFileName (removeExtension (Glib::path_get_basename (openThm->getFileName())));
    saveAsDialog->run ();
    Glib::ustring fname = saveAsDialog->getFileName ();
    if (fname=="")
        return;

    options.lastSaveAsPath = saveAsDialog->getDirectory ();
    options.saveAsDialogWidth = saveAsDialog->get_width();
    options.saveAsDialogHeight = saveAsDialog->get_height();

    SaveFormat sf = saveAsDialog->getFormat ();
    if (getExtension (fname)!=sf.format)
        fname = fname + "." + sf.format;

    if (saveAsDialog->getImmediately ()) {
        // check if it exists
        if (Glib::file_test (fname, Glib::FILE_TEST_EXISTS)) {
            Glib::ustring msg_ = Glib::ustring("<b>") + fname + ": " + M("MAIN_MSG_ALREADYEXISTS") + "\n" + M("MAIN_MSG_QOVERWRITE") + "</b>";
            Gtk::MessageDialog msgd (*parent, msg_, true, Gtk::MESSAGE_WARNING, Gtk::BUTTONS_YES_NO, true);
            int response = msgd.run ();
            if (response==Gtk::RESPONSE_NO)
                return;
        }
        // save image 
        rtengine::procparams::ProcParams pparams;
        ipc->getParams (&pparams);
        rtengine::ProcessingJob* job = rtengine::ProcessingJob::create (ipc->getInitialImage(), pparams);
        fname = removeExtension (fname);
        ProgressConnector<rtengine::IImage16*> *ld = new ProgressConnector<rtengine::IImage16*>();
        ld->startFunc(sigc::bind(sigc::ptr_fun(&rtengine::processImage), job, err, parent->getProgressListener() ),
        		      sigc::bind(sigc::mem_fun( *this,&EditorPanel::idle_saveImage ),ld,fname,sf,false ));
        saveimgas->set_sensitive(false);
        sendtogimp->set_sensitive(false);
    }
    else {
        BatchQueueEntry* bqe = createBatchQueueEntry ();
        bqe->outFileName = fname;
        bqe->saveFormat = saveAsDialog->getFormat ();
        parent->addBatchQueueJob (bqe, saveAsDialog->getToHeadOfQueue ());
    }
    // ask parent to redraw file browser
    // ... or does it automatically when the tab is switched to it
}

void EditorPanel::queueImgPressed () {

    saveProfile ();
    parent->addBatchQueueJob (createBatchQueueEntry ());
}

void EditorPanel::sendToGimpPressed () {
    // develop image
    rtengine::procparams::ProcParams pparams;
    ipc->getParams (&pparams);
    rtengine::ProcessingJob* job = rtengine::ProcessingJob::create (ipc->getInitialImage(), pparams);
    ProgressConnector<rtengine::IImage16*> *ld = new ProgressConnector<rtengine::IImage16*>();
    ld->startFunc(sigc::bind(sigc::ptr_fun(&rtengine::processImage), job, err, parent->getProgressListener() ),
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

        Glib::ustring fileName = Glib::ustring::compose ("%1.%2", fname, sf.format);

		int tries = 1;
		while (Glib::file_test (fileName, Glib::FILE_TEST_EXISTS) && tries<1000) {
			fileName = Glib::ustring::compose("%1-%2.%3", fname, tries, sf.format);
			tries++;
		}
		if (tries==1000){
			img->free ();
			return false;
		}

        ProgressConnector<int> *ld = new ProgressConnector<int>();
        img->setSaveProgressListener (parent->getProgressListener());
       	ld->startFunc (sigc::bind(sigc::mem_fun(img, &rtengine::IImage16::saveAsTIFF), fileName, sf.tiffBits),
        			   sigc::bind(sigc::mem_fun(*this,&EditorPanel::idle_sentToGimp), ld, img, fileName));
    }else{
    	gdk_threads_enter();
		Glib::ustring msg_ = Glib::ustring("<b> Error during image processing\n</b>");
		Gtk::MessageDialog msgd (*parent, msg_, true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
		msgd.run ();
		gdk_threads_leave ();
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
						// start gimp
						if (options.editorToSendTo==1) {
#ifdef _WIN32
								cmdLine = Glib::ustring("\"") + Glib::build_filename (Glib::build_filename(options.gimpDir,"bin"), "gimp-win-remote") + "\" gimp-2.4.exe" + " \"" + filename + "\"";
#else
								cmdLine = Glib::ustring("gimp-remote ") + " \"" + filename + "\"";
#endif
								success = safe_spawn_command_line_async (cmdLine);
								if (!success){
#ifdef _WIN32
										int ver = 12;
										while (!success && ver) {
												cmdLine = Glib::ustring("\"") + Glib::build_filename (Glib::build_filename(options.gimpDir,"bin"), Glib::ustring::compose("gimp-2.%1.exe",ver)) + "\" \"" + filename + "\"";
												ver--;
												success = safe_spawn_command_line_async (cmdLine);
										}
#elif defined __APPLE__
										cmdLine = Glib::ustring("gimp ") + " \"" + filename + "\"";
										success = safe_spawn_command_line_async (cmdLine);
#else
										cmdLine = Glib::ustring("gimp ") + " \"" + filename + "\"";
										success = safe_spawn_command_line_async (cmdLine);
#endif
								}
						}
						else if (options.editorToSendTo==2) {
#ifdef __APPLE__
								cmdLine = Glib::ustring("open -a \'") + Glib::build_filename(options.psDir,"Photoshop.app\' ")  + "\'" + filename + "\'";
#else
								cmdLine = Glib::ustring("\"") + Glib::build_filename(options.psDir,"Photoshop.exe") + "\" \"" + filename + "\"";
#endif
								success = safe_spawn_command_line_async (cmdLine);
						}
						else if (options.editorToSendTo==3) {
#ifdef __APPLE__
								cmdLine = Glib::ustring("") + options.customEditorProg + filename;
#else
								cmdLine = Glib::ustring("\"") + options.customEditorProg + "\" \"" + filename + "\"";
#endif
								success = safe_spawn_command_line_async (cmdLine);
						}

						if (!success) {
							    gdk_threads_enter ();
								Gtk::MessageDialog* msgd = new Gtk::MessageDialog (*parent, M("MAIN_MSG_CANNOTSTARTEDITOR"), false, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
								msgd->set_secondary_text (M("MAIN_MSG_CANNOTSTARTEDITOR_SECONDARY"));
								msgd->set_title (M("MAIN_BUTTON_SENDTOEDITOR"));
								msgd->run ();
								delete msgd;
								gdk_threads_leave ();
						}

    }

    return false;
}

void EditorPanel::saveOptions () {

    options.historyPanelWidth = hpanedl->get_position ();
    options.toolPanelWidth = vboxright->get_width ();
}

void EditorPanel::historyBeforeLineChanged (const rtengine::procparams::ProcParams& params) {

    if (beforeIpc) {
        ProcParams* pparams = beforeIpc->getParamsForUpdate (rtengine::EvProfileChanged);
        *pparams = params;
        beforeIpc->paramsUpdateReady ();
    }
}

void EditorPanel::beforeAfterToggled () {

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
        isrc->increaseRef ();
        beforeIpc = rtengine::StagedImageProcessor::create (isrc);
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

void EditorPanel::histogramChanged (unsigned int* rh, unsigned int* gh, unsigned int* bh, unsigned int* lh, unsigned int* bcrgb, unsigned int* bcl) {
    
    histogramPanel->histogramChanged (rh, gh, bh, lh); 
    tpc->updateCurveBackgroundHistogram (bcrgb, bcl);
}
