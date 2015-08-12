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
#include "editorpanel.h"
#include "options.h"
#include "progressconnector.h"
#include "rtwindow.h"
#include "guiutils.h"
#include "procparamchangers.h"
#include "../rtengine/safegtk.h"
#include "../rtengine/imagesource.h"
#include "soundman.h"
#include "rtimage.h"
#include <iostream>

using namespace rtengine::procparams;

EditorPanel::EditorPanel (FilePanel* filePanel)
    : realized(false), iHistoryShow(NULL), iHistoryHide(NULL), iTopPanel_1_Show(NULL), iTopPanel_1_Hide(NULL), iRightPanel_1_Show(NULL), iRightPanel_1_Hide(NULL), iBeforeLockON(NULL), iBeforeLockOFF(NULL), beforePreviewHandler(NULL), beforeIarea(NULL), beforeBox(NULL), afterBox(NULL), afterHeaderBox(NULL), parent(NULL), openThm(NULL), ipc(NULL), beforeIpc(NULL), isProcessing(false), catalogPane(NULL)
{

    epih = new EditorPanelIdleHelper;
    epih->epanel = this;
    epih->destroyed = false;
    epih->pending = 0;
    //rtengine::befaf=true;
    processingStartedTime = 0;
    firstProcessingDone = false;

    // construct toolpanelcoordinator
    tpc = new ToolPanelCoordinator ();

    // build GUI

    // build left side panel
    leftbox = new Gtk::VBox ();
    leftbox->set_border_width (2);
    leftbox->set_size_request(100, 250);

    histogramPanel = NULL;

    profilep = Gtk::manage (new ProfilePanel ());
    ppframe = new Gtk::Frame ();
    ppframe->add (*profilep);
    ppframe->set_label (M("PROFILEPANEL_LABEL"));
    //leftbox->pack_start (*ppframe, Gtk::PACK_SHRINK, 4);

    navigator = Gtk::manage (new Navigator ());
    navigator->previewWindow->set_size_request (-1, 150);
    leftbox->pack_start (*navigator, Gtk::PACK_SHRINK, 2);

    history = Gtk::manage (new History ());
    leftbox->pack_start (*history);

    leftbox->show_all ();

    // build the middle of the screen
    Gtk::VBox* editbox = Gtk::manage (new Gtk::VBox ());

    info = Gtk::manage (new Gtk::ToggleButton ());
    Gtk::Image* infoimg = Gtk::manage (new RTImage ("info.png"));
    info->add (*infoimg);
    info->set_relief(Gtk::RELIEF_NONE);
    info->set_tooltip_markup (M("MAIN_TOOLTIP_QINFO"));

    beforeAfter = Gtk::manage (new Gtk::ToggleButton ());
    Gtk::Image* beforeAfterIcon = Gtk::manage (new RTImage ("beforeafter.png"));
    beforeAfter->add(*beforeAfterIcon);
    beforeAfter->set_relief(Gtk::RELIEF_NONE);
    beforeAfter->set_tooltip_markup (M("MAIN_TOOLTIP_TOGGLE"));

    iBeforeLockON = new RTImage ("lock-on.png");
    iBeforeLockOFF = new RTImage ("lock-off.png");

    Gtk::VSeparator* vsept = Gtk::manage (new Gtk::VSeparator ());
    Gtk::VSeparator* vsepz = Gtk::manage (new Gtk::VSeparator ());
    Gtk::VSeparator* vsepi = Gtk::manage (new Gtk::VSeparator ());
    Gtk::VSeparator* vseph = Gtk::manage (new Gtk::VSeparator ());

    hidehp = Gtk::manage (new Gtk::ToggleButton ());

    iHistoryShow = new RTImage ("panel-to-right.png");
    iHistoryHide = new RTImage ("panel-to-left.png");

    hidehp->set_relief(Gtk::RELIEF_NONE);
    hidehp->set_active (options.showHistory);
    hidehp->set_tooltip_markup (M("MAIN_TOOLTIP_HIDEHP"));

    if (options.showHistory) {
        hidehp->set_image (*iHistoryHide);
    } else {
        hidehp->set_image (*iHistoryShow);
    }

    tbTopPanel_1 = NULL;

    if (!simpleEditor && filePanel) {
        tbTopPanel_1 = new Gtk::ToggleButton ();
        iTopPanel_1_Show = new RTImage ("panel-to-bottom.png");
        iTopPanel_1_Hide = new RTImage ("panel-to-top.png");
        tbTopPanel_1->set_relief(Gtk::RELIEF_NONE);
        tbTopPanel_1->set_active (true);
        tbTopPanel_1->set_tooltip_markup (M("MAIN_TOOLTIP_SHOWHIDETP1"));
        tbTopPanel_1->set_image (*iTopPanel_1_Hide);
    }

    tbRightPanel_1 = new Gtk::ToggleButton ();
    iRightPanel_1_Show = new RTImage ("panel-to-left.png");
    iRightPanel_1_Hide = new RTImage ("panel-to-right.png");
    tbRightPanel_1->set_relief(Gtk::RELIEF_NONE);
    tbRightPanel_1->set_active (true);
    tbRightPanel_1->set_tooltip_markup (M("MAIN_TOOLTIP_SHOWHIDERP1"));
    tbRightPanel_1->set_image (*iRightPanel_1_Hide);

    Gtk::VSeparator* vsepcl = Gtk::manage (new Gtk::VSeparator ());
    Gtk::VSeparator* vsepz2 = Gtk::manage (new Gtk::VSeparator ());
    Gtk::VSeparator* vsepz3 = Gtk::manage (new Gtk::VSeparator ());
    Gtk::VSeparator* vsepz4 = Gtk::manage (new Gtk::VSeparator ());

    iareapanel = new ImageAreaPanel ();
    tpc->setEditProvider(iareapanel->imageArea);

    Gtk::HBox* toolBarPanel = Gtk::manage (new Gtk::HBox ());
    toolBarPanel->pack_start (*hidehp, Gtk::PACK_SHRINK, 1);
    toolBarPanel->pack_start (*vseph, Gtk::PACK_SHRINK, 2);
    toolBarPanel->pack_start (*info, Gtk::PACK_SHRINK, 1);
    toolBarPanel->pack_start (*beforeAfter, Gtk::PACK_SHRINK, 1);
    toolBarPanel->pack_start (*vsepi, Gtk::PACK_SHRINK, 2);
    toolBarPanel->pack_start (*tpc->getToolBar(), Gtk::PACK_SHRINK, 1);
    toolBarPanel->pack_start (*vsept, Gtk::PACK_SHRINK, 2);

    if (tbTopPanel_1) {
        toolBarPanel->pack_end   (*tbTopPanel_1, Gtk::PACK_SHRINK, 1);
        Gtk::VSeparator* vsep1 = Gtk::manage (new Gtk::VSeparator ());
        toolBarPanel->pack_end   (*vsep1, Gtk::PACK_SHRINK, 2);
    }

    toolBarPanel->pack_end   (*tpc->coarse, Gtk::PACK_SHRINK, 2);
    toolBarPanel->pack_end   (*vsepcl, Gtk::PACK_SHRINK, 2);
    toolBarPanel->pack_end   (*iareapanel->imageArea->indClippedPanel, Gtk::PACK_SHRINK, 0);
    toolBarPanel->pack_end   (*vsepz, Gtk::PACK_SHRINK, 2);
    toolBarPanel->pack_end   (*iareapanel->imageArea->previewModePanel, Gtk::PACK_SHRINK, 0);
    toolBarPanel->pack_end   (*vsepz4, Gtk::PACK_SHRINK, 2);

    afterBox = Gtk::manage (new Gtk::VBox ());
    afterBox->pack_start (*iareapanel);

    beforeAfterBox = Gtk::manage (new Gtk::HBox());
    beforeAfterBox->pack_start (*afterBox);

    editbox->pack_start (*toolBarPanel, Gtk::PACK_SHRINK, 0);
    editbox->pack_start (*beforeAfterBox);

    // build right side panel
    vboxright = new Gtk::VBox (false, 0);
    vboxright->set_size_request(100, 250);

    vboxright->set_border_width (2);

    vboxright->pack_start (*ppframe, Gtk::PACK_SHRINK, 2);
    // main notebook
    vboxright->pack_start (*tpc->toolPanelNotebook);

    // Save buttons
    Gtk::HBox* iops = Gtk::manage (new Gtk::HBox ());

    Gtk::Image *saveButtonImage =  Gtk::manage (new RTImage ("gtk-save-large.png"));
    saveimgas = Gtk::manage (new Gtk::Button ());
    saveimgas->add(*saveButtonImage);
    saveimgas->set_tooltip_markup(M("MAIN_BUTTON_SAVE_TOOLTIP"));

    Gtk::Image *queueButtonImage = Gtk::manage (new RTImage ("processing.png"));
    queueimg = Gtk::manage (new Gtk::Button ());
    queueimg->add(*queueButtonImage);
    queueimg->set_tooltip_markup(M("MAIN_BUTTON_PUTTOQUEUE_TOOLTIP"));

    Gtk::Image *sendToEditorButtonImage = Gtk::manage (new RTImage ("image-editor.png"));
    sendtogimp = Gtk::manage (new Gtk::Button ());
    sendtogimp->add(*sendToEditorButtonImage);
    sendtogimp->set_tooltip_markup(M("MAIN_BUTTON_SENDTOEDITOR_TOOLTIP"));

    iops->pack_start (*saveimgas, Gtk::PACK_SHRINK);

    if(!simpleEditor) {
        iops->pack_start (*queueimg, Gtk::PACK_SHRINK);
    }

    iops->pack_start (*sendtogimp, Gtk::PACK_SHRINK);

    // Status box
    statusBox = Gtk::manage (new Gtk::HBox ());
    cssProvider = Gtk::CssProvider::create();
    progressLabel = Gtk::manage (new Gtk::ProgressBar());
    progressLabel->set_show_text(true);
    setExpandAlignProperties(progressLabel, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_CENTER);
    progressLabel->set_fraction(0.0);

    if (cssProvider) {
        progressLabel->get_style_context()->add_provider (cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);    // Can't find the gtkmm version of the enum!
    }

    statusBox->pack_start (*progressLabel);
    iops->pack_start(*statusBox, Gtk::PACK_SHRINK, 2);

    // tbRightPanel_1
    iops->pack_end (*tbRightPanel_1, Gtk::PACK_SHRINK, 0);

    // ShowHideSidePanels
    tbShowHideSidePanels = new Gtk::ToggleButton ();
    iShowHideSidePanels = new RTImage ("crossed-arrows-out.png");
    iShowHideSidePanels_exit = new RTImage ("crossed-arrows-in.png");
    tbShowHideSidePanels->set_relief(Gtk::RELIEF_NONE);
    tbShowHideSidePanels->set_active (false);
    tbShowHideSidePanels->set_tooltip_markup (M("MAIN_BUTTON_SHOWHIDESIDEPANELS_TOOLTIP"));
    tbShowHideSidePanels->set_image (*iShowHideSidePanels);
    iops->pack_end (*tbShowHideSidePanels, Gtk::PACK_SHRINK, 0);
    iops->pack_end (*vsepz2, Gtk::PACK_SHRINK, 1);

    // Zoom panel
    iops->pack_end (*iareapanel->imageArea->zoomPanel, Gtk::PACK_SHRINK, 1);
    iops->pack_end (*vsepz3, Gtk::PACK_SHRINK, 2);

    navPrev = navNext = navSync = NULL;

    if (!simpleEditor && !options.tabbedUI) {
        // Navigation buttons
        Gtk::Image *navPrevImage = Gtk::manage (new RTImage ("nav-prev.png"));
        navPrevImage->set_padding(0, 0);
        navPrev = Gtk::manage (new Gtk::Button ());
        navPrev->add(*navPrevImage);
        navPrev->set_relief(Gtk::RELIEF_NONE);
        navPrev->set_tooltip_markup(M("MAIN_BUTTON_NAVPREV_TOOLTIP"));

        Gtk::Image *navNextImage = Gtk::manage (new RTImage ("nav-next.png"));
        navNextImage->set_padding(0, 0);
        navNext = Gtk::manage (new Gtk::Button ());
        navNext->add(*navNextImage);
        navNext->set_relief(Gtk::RELIEF_NONE);
        navNext->set_tooltip_markup(M("MAIN_BUTTON_NAVNEXT_TOOLTIP"));

        Gtk::Image *navSyncImage = Gtk::manage (new RTImage ("nav-sync.png"));
        navSyncImage->set_padding(0, 0);
        navSync = Gtk::manage (new Gtk::Button ());
        navSync->add(*navSyncImage);
        navSync->set_relief(Gtk::RELIEF_NONE);
        navSync->set_tooltip_markup(M("MAIN_BUTTON_NAVSYNC_TOOLTIP"));

        iops->pack_end (*Gtk::manage(new Gtk::VSeparator()), Gtk::PACK_SHRINK, 0);
        iops->pack_end (*navNext, Gtk::PACK_SHRINK, 0);
        iops->pack_end (*navSync, Gtk::PACK_SHRINK, 0);
        iops->pack_end (*navPrev, Gtk::PACK_SHRINK, 0);
    }

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

    if(filePanel) {
        catalogPane = new Gtk::Paned();
        viewpaned->pack1(*catalogPane, false, true);
    }

    viewpaned->pack2(*editbox, true, true);


    Gtk::Frame* vbfr = Gtk::manage (new Gtk::Frame ());
    vbfr->add (*viewpaned);
    vbfr->set_size_request(100, 250);
    hpanedl->pack2(*vbfr, true, true);

    hpanedr->pack1(*hpanedl, true, true);
    hpanedr->pack2(*vboxright, false, true);
    hpanedl->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &EditorPanel::leftPaneButtonReleased) );
    hpanedr->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &EditorPanel::rightPaneButtonReleased) );

    pack_start (*hpanedr);

    updateHistogramPosition (0, options.histogramPosition);

    show_all ();
    /*
        // save as dialog
        if (safe_file_test (options.lastSaveAsPath, Glib::FILE_TEST_IS_DIR))
            saveAsDialog = new SaveAsDialog (options.lastSaveAsPath);
        else
            saveAsDialog = new SaveAsDialog (safe_get_user_picture_dir());

        saveAsDialog->set_default_size (options.saveAsDialogWidth, options.saveAsDialogHeight);
    */
    // connect listeners
    profilep->setProfileChangeListener (tpc);
    history->setProfileChangeListener (tpc);
    history->setHistoryBeforeLineListener (this);
    tpc->addPParamsChangeListener (profilep);
    tpc->addPParamsChangeListener (history);
    tpc->addPParamsChangeListener (this);
    iareapanel->imageArea->setCropGUIListener (tpc->getCropGUIListener());
    iareapanel->imageArea->setPointerMotionListener (navigator);
    iareapanel->imageArea->setImageAreaToolListener (tpc);

    // initialize components
    info->set_active (options.showInfo);
    tpc->readOptions ();

    // connect event handlers
    info->signal_toggled().connect( sigc::mem_fun(*this, &EditorPanel::info_toggled) );
    beforeAfter->signal_toggled().connect( sigc::mem_fun(*this, &EditorPanel::beforeAfterToggled) );
    hidehp->signal_toggled().connect( sigc::mem_fun(*this, &EditorPanel::hideHistoryActivated) );
    tbRightPanel_1->signal_toggled().connect( sigc::mem_fun(*this, &EditorPanel::tbRightPanel_1_toggled) );
    saveimgas->signal_pressed().connect( sigc::mem_fun(*this, &EditorPanel::saveAsPressed) );
    queueimg->signal_pressed().connect( sigc::mem_fun(*this, &EditorPanel::queueImgPressed) );
    sendtogimp->signal_pressed().connect( sigc::mem_fun(*this, &EditorPanel::sendToGimpPressed) );

    if(navPrev) {
        navPrev->signal_pressed().connect( sigc::mem_fun(*this, &EditorPanel::openPreviousEditorImage) );
    }

    if(navNext) {
        navNext->signal_pressed().connect( sigc::mem_fun(*this, &EditorPanel::openNextEditorImage) );
    }

    if(navSync) {
        navSync->signal_pressed().connect( sigc::mem_fun(*this, &EditorPanel::syncFileBrowser) );
    }

    ShowHideSidePanelsconn = tbShowHideSidePanels->signal_toggled().connect ( sigc::mem_fun(*this, &EditorPanel::toggleSidePanels), true);

    if (tbTopPanel_1) {
        tbTopPanel_1->signal_toggled().connect( sigc::mem_fun(*this, &EditorPanel::tbTopPanel_1_toggled) );
    }
}

EditorPanel::~EditorPanel ()
{

    history->setHistoryBeforeLineListener (NULL);
    // the order is important!
    iareapanel->setBeforeAfterViews (NULL, iareapanel);
    delete iareapanel;
    iareapanel = NULL;

    if (beforeIpc) {
        beforeIpc->stopProcessing ();
    }

    delete beforeIarea;
    beforeIarea = NULL;

    if (beforeIpc) {
        beforeIpc->setPreviewImageListener (NULL);
    }

    delete beforePreviewHandler;
    beforePreviewHandler = NULL;

    if (beforeIpc) {
        rtengine::StagedImageProcessor::destroy (beforeIpc);
    }

    beforeIpc = NULL;

    close ();

    if (epih->pending) {
        epih->destroyed = true;
    } else {
        delete epih;
    }

    delete tpc;

    delete ppframe;
    delete leftbox;
    delete vboxright;

    //delete saveAsDialog;
    if(catalogPane) {
        delete catalogPane;
    }

    if (iTopPanel_1_Show) {
        delete iTopPanel_1_Show;
    }

    if (iTopPanel_1_Hide) {
        delete iTopPanel_1_Hide;
    }

    if (iHistoryShow) {
        delete iHistoryShow;
    }

    if (iHistoryHide) {
        delete iHistoryHide;
    }

    if(iBeforeLockON) {
        delete iBeforeLockON;
    }

    if(iBeforeLockOFF) {
        delete iBeforeLockOFF;
    }

    if(iRightPanel_1_Show) {
        delete iRightPanel_1_Show;
    }

    if(iRightPanel_1_Hide) {
        delete iRightPanel_1_Hide;
    }
}

void EditorPanel::leftPaneButtonReleased(GdkEventButton *event)
{
    if (event->button == 1) {
        // Button 1 released : it's a resize
        options.historyPanelWidth = hpanedl->get_position();
    }

    /*else if (event->button == 3) {
    }*/
}

void EditorPanel::rightPaneButtonReleased(GdkEventButton *event)
{
    if (event->button == 1) {
        int winW, winH;
        parent->get_size(winW, winH);
        // Button 1 released : it's a resize
        options.toolPanelWidth = winW - hpanedr->get_position();
    }

    /*else if (event->button == 3) {
    }*/
}

void EditorPanel::writeOptions()
{
    if (profilep) {
        profilep->writeOptions();
    }

    if (tpc) {
        tpc->writeOptions();
    }
}

void EditorPanel::showTopPanel(bool show)
{
    if (tbTopPanel_1->get_active() != show) {
        tbTopPanel_1->set_active(show);
    }
}

void EditorPanel::setAspect ()
{
    int winW, winH;
    parent->get_size(winW, winH);
    hpanedl->set_position(options.historyPanelWidth);
    hpanedr->set_position(winW - options.toolPanelWidth);

    // initialize components
    if (info->get_active() != options.showInfo) {
        info->set_active (options.showInfo);
    }
}

void EditorPanel::on_realize ()
{
    realized = true;
    Gtk::VBox::on_realize ();
    // This line is needed to avoid autoexpansion of the window :-/
    //vboxright->set_size_request (options.toolPanelWidth, -1);
    tpc->updateToolState();
}

void EditorPanel::open (Thumbnail* tmb, rtengine::InitialImage* isrc)
{

    close();

    isProcessing = true; // prevents closing-on-init

    // initialize everything
    openThm = tmb;
    openThm->increaseRef ();

    fname = openThm->getFileName();
    lastSaveAsFileName = removeExtension (Glib::path_get_basename (fname));

    previewHandler = new PreviewHandler ();
    previewHandler2 = new PreviewHandler ();

    this->isrc = isrc;
    ipc = rtengine::StagedImageProcessor::create (isrc);
    ipc->setProgressListener (this);
    ipc->setPreviewImageListener (previewHandler);
    ipc->setPreviewScale (10);  // Important
    tpc->initImage (ipc, tmb->getType() == FT_Raw);
    ipc->setHistogramListener (this);

//    iarea->fitZoom ();   // tell to the editorPanel that the next image has to be fitted to the screen
    iareapanel->imageArea->setPreviewHandler (previewHandler);
    iareapanel->imageArea->setImProcCoordinator (ipc);
    navigator->previewWindow->setPreviewHandler (previewHandler);
    navigator->previewWindow->setImageArea (iareapanel->imageArea);

    rtengine::ImageSource* is = isrc->getImageSource();
    is->setProgressListener( this );

    // try to load the last saved parameters from the cache or from the paramfile file
    ProcParams* ldprof = openThm->createProcParamsForUpdate(true, false); // will be freed by initProfile

    // initialize profile
    Glib::ustring defProf = openThm->getType() == FT_Raw ? options.defProfRaw : options.defProfImg;
    profilep->initProfile (defProf, ldprof);
    profilep->setInitialFileName (fname);

    openThm->addThumbnailListener (this);
    info_toggled ();

    if (beforeIarea) {
        beforeAfterToggled();
        beforeAfterToggled();
    }

    // If in single tab mode, the main crop window is not constructed the very first time
    // since there was no resize event
    if (iareapanel->imageArea->mainCropWindow) {
        iareapanel->imageArea->mainCropWindow->cropHandler.newImage(ipc, false);
        iareapanel->imageArea->mainCropWindow->initialImageArrived();

        // In single tab mode, the image is not always updated between switches
        // normal redraw don't work, so this is the hard way
        // Disabled this with Issue 2435 because it seems to work fine now
//        if (!options.tabbedUI && iareapanel->imageArea->mainCropWindow->getZoomFitVal() == 1.0) {
//          iareapanel->imageArea->mainCropWindow->cropHandler.update();
//        }
    } else {
        Gtk::Allocation alloc;
        iareapanel->imageArea->on_resized(alloc);
    }
}

void EditorPanel::close ()
{
    if (ipc) {
        saveProfile ();
        // close image processor and the current thumbnail
        tpc->closeImage ();    // this call stops image processing
        tpc->writeOptions ();
        rtengine::ImageSource* is = isrc->getImageSource();
        is->setProgressListener( NULL );

        if (ipc) {
            ipc->setPreviewImageListener (NULL);
        }

        if (beforeIpc) {
            beforeIpc->setPreviewImageListener (NULL);
        }

        delete previewHandler;
        previewHandler = NULL;

        if(iareapanel) {
            iareapanel->imageArea->setPreviewHandler (NULL);
            iareapanel->imageArea->setImProcCoordinator (NULL);
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

void EditorPanel::saveProfile ()
{
    if (!ipc || !openThm) {
        return;
    }

    // If the file was deleted, do not generate ghost entries
    if (safe_file_test(fname, Glib::FILE_TEST_EXISTS)) {
        ProcParams params;
        ipc->getParams (&params);

        // Will call updateCache, which will update both the cached and sidecar files if necessary
        openThm->setProcParams (params, NULL, EDITOR);
    }
}

Glib::ustring EditorPanel::getShortName ()
{
    if (openThm) {
        return Glib::path_get_basename (openThm->getFileName ());
    } else {
        return "";
    }
}

Glib::ustring EditorPanel::getFileName ()
{
    if (openThm) {
        return openThm->getFileName ();
    } else {
        return "";
    }
}

// TODO!!!
void EditorPanel::procParamsChanged (rtengine::procparams::ProcParams* params, rtengine::ProcEvent ev, Glib::ustring descr, ParamsEdited* paramsEdited)
{

//    if (ev!=EvPhotoLoaded)
//        saveLabel->set_markup (Glib::ustring("<span foreground=\"#AA0000\" weight=\"bold\">") + M("MAIN_BUTTON_SAVE") + "</span>");
}

struct spsparams {
    bool inProcessing;
    EditorPanelIdleHelper* epih;
};

int setProgressStateUIThread (void* data)
{

    spsparams* p = static_cast<spsparams*>(data);

    if (p->epih->destroyed) {
        if (p->epih->pending == 1) {
            delete p->epih;
        } else {
            p->epih->pending--;
        }

        delete p;

        return 0;
    }

    p->epih->epanel->refreshProcessingState (p->inProcessing);
    p->epih->pending--;
    delete p;

    return 0;
}

void EditorPanel::setProgressState (bool inProcessing)
{

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
    Glib::RefPtr<Gtk::CssProvider> cssProvider;

};

int setprogressStrUI( void *p )
{
    spparams *s = static_cast<spparams*>(p);

    if( ! s->str.empty() ) {
        s->pProgress->set_text( M(s->str) );
    }

    if( s->val >= 0 ) {
        s->pProgress->set_fraction( s->val );

        if (s->cssProvider) {
            if( s->val < 1.0 ) {
                s->cssProvider->load_from_data("ProgressBar { background-color: red }");
            } else {
                s->cssProvider->load_from_data("ProgressBar { background-color: grey }");
            }

            s->pProgress->get_style_context()->set_background(s->pProgress->get_window());
        }
    }

    delete s;
    return 0;
}

void EditorPanel::setProgress (double p)
{
    spparams *s = new spparams;
    s->val = p;
    s->pProgress = progressLabel;
    add_idle (setprogressStrUI, s);
}

void EditorPanel::setProgressStr (Glib::ustring str)
{
    spparams *s = new spparams;
    s->str = str;
    s->val = -1;
    s->pProgress = progressLabel;
    add_idle (setprogressStrUI, s);
}

// This is only called from the ThreadUI, so within the gtk thread
void EditorPanel::refreshProcessingState (bool inProcessingP)
{
    spparams *s = new spparams;
    s->pProgress = progressLabel;

    if (inProcessingP) {
        if (processingStartedTime == 0) {
            processingStartedTime = ::time(NULL);
        }

        s->str = "PROGRESSBAR_PROCESSING";
        s->val = 0.0;
    } else {
        // Set proc params of thumbnail. It saves it into the cache and updates the file browser.
        if (ipc && openThm && tpc->getChangedState()) {
            rtengine::procparams::ProcParams pparams;
            ipc->getParams (&pparams);
            openThm->setProcParams (pparams, NULL, EDITOR, false);
        }

        // Ring a sound if it was a long event
        if (processingStartedTime != 0) {
            time_t curTime = ::time(NULL);

            if (::difftime(curTime, processingStartedTime) > options.sndLngEditProcDoneSecs) {
                SoundManager::playSoundAsync(options.sndLngEditProcDone);
            }

            processingStartedTime = 0;
        }

        // Set progress bar "done"
        s->str = "PROGRESSBAR_READY";
        s->val = 1.0;

#ifdef WIN32

        // Maybe accessing "parent", which is a Gtk object, can justify to get the Gtk lock...
        if (!firstProcessingDone && static_cast<RTWindow*>(parent)->getIsFullscreen()) {
            parent->fullscreen();
        }

#endif
        firstProcessingDone = true;
    }

    isProcessing = inProcessingP;

    setprogressStrUI(s);
}

struct errparams {
    Glib::ustring descr;
    Glib::ustring title;
    EditorPanelIdleHelper* epih;
};

void EditorPanel::displayError (Glib::ustring title, Glib::ustring descr)
{
    GtkWidget* msgd = gtk_message_dialog_new_with_markup (NULL,
                      GTK_DIALOG_DESTROY_WITH_PARENT,
                      GTK_MESSAGE_ERROR,
                      GTK_BUTTONS_OK,
                      "<b>%s</b>",
                      descr.data());
    gtk_window_set_title((GtkWindow*)msgd, title.data());
    g_signal_connect_swapped (msgd, "response",
                              G_CALLBACK (gtk_widget_destroy),
                              msgd);
    gtk_widget_show_all (msgd);
}

int disperrorUI (void* data)
{
    errparams* p = static_cast<errparams*>(data);

    if (p->epih->destroyed) {
        if (p->epih->pending == 1) {
            delete p->epih;
        } else {
            p->epih->pending--;
        }

        delete p;

        return 0;
    }

    p->epih->epanel->displayError (p->title, p->descr);
    p->epih->pending--;
    delete p;

    return 0;
}

void EditorPanel::error (Glib::ustring title, Glib::ustring descr)
{

    epih->pending++;
    errparams* p = new errparams;
    p->descr = descr;
    p->title = title;
    p->epih = epih;
    g_idle_add (disperrorUI, p);
}

void EditorPanel::info_toggled ()
{

    Glib::ustring infoString;
    Glib::ustring infoString1; //1-st line
    Glib::ustring infoString2; //2-nd line
    Glib::ustring infoString3; //3-rd line
    Glib::ustring infoString4; //4-th line
    Glib::ustring expcomp;

    if (!ipc || !openThm) {
        return;
    }

    const rtengine::ImageMetaData* idata = ipc->getInitialImage()->getMetaData();

    if (idata && idata->hasExif()) {
        infoString1 = Glib::ustring::compose ("%1 + %2",
                                              Glib::ustring(idata->getMake() + " " + idata->getModel()),
                                              Glib::ustring(idata->getLens()));

        infoString2 = Glib::ustring::compose ("<span size=\"small\">f/</span><span size=\"large\">%1</span>  <span size=\"large\">%2</span><span size=\"small\">s</span>  <span size=\"small\">%3</span><span size=\"large\">%4</span>  <span size=\"large\">%5</span><span size=\"small\">mm</span>",
                                              Glib::ustring(idata->apertureToString(idata->getFNumber())),
                                              Glib::ustring(idata->shutterToString(idata->getShutterSpeed())),
                                              M("QINFO_ISO"), idata->getISOSpeed(),
                                              Glib::ustring::format(std::setw(3), std::fixed, std::setprecision(2), idata->getFocalLen()));

        expcomp = Glib::ustring(idata->expcompToString(idata->getExpComp(), true)); // maskZeroexpcomp

        if (expcomp != "") {
            infoString2 = Glib::ustring::compose("%1  <span size=\"large\">%2</span><span size=\"small\">EV</span>",
                                                 infoString2,
                                                 expcomp /*Glib::ustring(idata->expcompToString(idata->getExpComp()))*/);
        }

        infoString3 = Glib::ustring::compose ("<span size=\"small\">%1</span><span>%2</span>",
                                              escapeHtmlChars(Glib::path_get_dirname(openThm->getFileName())) + G_DIR_SEPARATOR_S,
                                              escapeHtmlChars(Glib::path_get_basename(openThm->getFileName()))  );

        int ww = ipc->getFullWidth();
        int hh = ipc->getFullHeight();
        //megapixels
        infoString4 = Glib::ustring::compose ("<span size=\"small\">%1 MP (%2x%3)</span>", Glib::ustring::format(std::setw(4), std::fixed, std::setprecision(1), (float)ww * hh / 1000000), ww, hh);

        infoString = Glib::ustring::compose ("%1\n%2\n%3\n%4", infoString1, infoString2, infoString3, infoString4);
    } else {
        infoString = M("QINFO_NOEXIF");
    }

    iareapanel->imageArea->setInfoText (infoString);
    iareapanel->imageArea->infoEnabled (info->get_active ());
}

void EditorPanel::hideHistoryActivated ()
{

    removeIfThere (hpanedl, leftbox, false);

    if (hidehp->get_active()) {
        hpanedl->pack1 (*leftbox, false, true);
    }

    options.showHistory = hidehp->get_active();

    if (options.showHistory) {
        hidehp->set_image (*iHistoryHide);
    } else {
        hidehp->set_image (*iHistoryShow);
    }

    tbShowHideSidePanels_managestate();
}


void EditorPanel::tbRightPanel_1_toggled ()
{
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
    if (vboxright) {
        if (tbRightPanel_1->get_active()) {
            vboxright->show();
            tbRightPanel_1->set_image (*iRightPanel_1_Hide);
        } else {
            vboxright->hide();
            tbRightPanel_1->set_image (*iRightPanel_1_Show);
        }

        tbShowHideSidePanels_managestate();
    }
}

void EditorPanel::tbTopPanel_1_visible (bool visible)
{
    if (!tbTopPanel_1) {
        return;
    }

    if (visible) {
        tbTopPanel_1->show();
    } else {
        tbTopPanel_1->hide();
    }
}

void EditorPanel::tbTopPanel_1_toggled ()
{

    if (catalogPane) { // catalogPane does not exist in multitab mode

        if (tbTopPanel_1->get_active()) {
            catalogPane->show();
            tbTopPanel_1->set_image (*iTopPanel_1_Hide);
            options.editorFilmStripOpened = true;
        } else {
            catalogPane->hide();
            tbTopPanel_1->set_image (*iTopPanel_1_Show);
            options.editorFilmStripOpened = false;
        }

        tbShowHideSidePanels_managestate();
    }
}

/*
 * WARNING: Take care of the simpleEditor value when adding or modifying shortcut keys,
 *          since handleShortcutKey is now also triggered in simple editor mode
 */
bool EditorPanel::handleShortcutKey (GdkEventKey* event)
{

    bool ctrl = event->state & GDK_CONTROL_MASK;
    bool shift = event->state & GDK_SHIFT_MASK;
    bool alt = event->state & GDK_MOD1_MASK;
#ifdef __WIN32__
    bool altgr = event->state & GDK_MOD2_MASK;
#else
    bool altgr = event->state & GDK_MOD5_MASK;
#endif

    // Editor Layout
    switch(event->keyval) {
    case GDK_KEY_L:
        if (tbTopPanel_1) {
            tbTopPanel_1->set_active (!tbTopPanel_1->get_active());    // toggle top panel
        }

        if (ctrl) {
            hidehp->set_active (!hidehp->get_active());    // toggle History (left panel)
        }

        if (alt) {
            tbRightPanel_1->set_active (!tbRightPanel_1->get_active());    // toggle right panel
        }

        return true;
        break;

    case GDK_KEY_l:
        if (!shift && !alt /*&& !ctrl*/) {
            hidehp->set_active (!hidehp->get_active()); // toggle History (left panel)
            return true;
        }

        if (alt && !ctrl) { // toggle right panel
            tbRightPanel_1->set_active (!tbRightPanel_1->get_active());
            return true;
        }

        if (alt && ctrl) { // toggle left and right panels
            hidehp->set_active (!hidehp->get_active());
            tbRightPanel_1->set_active (!tbRightPanel_1->get_active());
            return true;
        }

        break;

    case GDK_KEY_m: // Maximize preview panel: hide top AND right AND history panels
        if (!ctrl && !alt) {
            toggleSidePanels();
            return true;
        }

        break;

    case GDK_KEY_M: // Maximize preview panel: hide top AND right AND history panels AND (fit image preview)
        if (!ctrl && !alt) {
            toggleSidePanelsZoomFit();
            return true;
        }

        break;
    }

#ifdef __WIN32__

    if (!alt && !ctrl && !altgr && event->hardware_keycode == 0x39 ) {
        iareapanel->imageArea->previewModePanel->togglebackColor();
        return true;
    }

#else

    if (!alt && !ctrl && !altgr && event->hardware_keycode == 0x12 ) {
        iareapanel->imageArea->previewModePanel->togglebackColor();
        return true;
    }

#endif

    if (!alt) {
        if (!ctrl) {
            // Normal
            switch(event->keyval) {
            case GDK_KEY_bracketright:
                tpc->coarse->rotateRight();
                return true;

            case GDK_KEY_bracketleft:
                tpc->coarse->rotateLeft();
                return true;

            case GDK_KEY_i:
            case GDK_KEY_I:
                info->set_active (!info->get_active());
                return true;

            case GDK_KEY_B:
                beforeAfter->set_active (!beforeAfter->get_active());
                return true;

            case GDK_KEY_plus:
            case GDK_KEY_equal:
            case GDK_KEY_KP_Add:
                iareapanel->imageArea->zoomPanel->zoomInClicked();
                return true;

            case GDK_KEY_minus:
            case GDK_KEY_underscore:
            case GDK_KEY_KP_Subtract:
                iareapanel->imageArea->zoomPanel->zoomOutClicked();
                return true;

            case GDK_KEY_z://GDK_1
                iareapanel->imageArea->zoomPanel->zoom11Clicked();
                return true;

            /*
            #ifndef __WIN32__
                            case GDK_KEY_9: // toggle background color of the preview
                                iareapanel->imageArea->previewModePanel->togglebackColor();
                                return true;
            #endif
            */
            case GDK_KEY_r: //preview mode Red
                iareapanel->imageArea->previewModePanel->toggleR();
                return true;

            case GDK_KEY_g: //preview mode Green
                iareapanel->imageArea->previewModePanel->toggleG();
                return true;

            case GDK_KEY_b: //preview mode Blue
                iareapanel->imageArea->previewModePanel->toggleB();
                return true;

            case GDK_KEY_v: //preview mode Luminosity
                iareapanel->imageArea->previewModePanel->toggleL();
                return true;

            case GDK_KEY_F: //preview mode Focus Mask
                iareapanel->imageArea->previewModePanel->toggleFocusMask();
                return true;

            case GDK_KEY_f:
                iareapanel->imageArea->zoomPanel->zoomFitClicked();
                return true;

            case GDK_KEY_less:
                iareapanel->imageArea->indClippedPanel->toggleClipped(true);
                return true;

            case GDK_KEY_greater:
                iareapanel->imageArea->indClippedPanel->toggleClipped(false);
                return true;

            case GDK_KEY_F5:
                openThm->openDefaultViewer((event->state & GDK_SHIFT_MASK) ? 2 : 1);
                return true;

            case GDK_KEY_y: // synchronize filebrowser with image in Editor
                if (!simpleEditor && fPanel && !fname.empty()) {
                    fPanel->fileCatalog->selectImage(fname, false);
                    return true;
                }

                break; // to avoid gcc complain

            case GDK_KEY_x: // clear filters and synchronize filebrowser with image in Editor
                if (!simpleEditor && fPanel && !fname.empty()) {
                    fPanel->fileCatalog->selectImage(fname, true);
                    return true;
                }

                break; // to avoid gcc complain
            }
        } else {
            // With control
            switch (event->keyval) {
            case GDK_KEY_S:
                saveProfile();
                setProgressStr(M("PROGRESSBAR_PROCESSING_PROFILESAVED"));
                return true;

            case GDK_KEY_s:
                saveAsPressed();
                return true;

            case GDK_KEY_b:
                if (!simpleEditor) {
                    queueImgPressed();
                }

                return true;

            case GDK_KEY_e:
                sendToGimpPressed();
                return true;

            case GDK_KEY_z:
                history->undo ();
                return true;

            case GDK_KEY_Z:
                history->redo ();
                return true;

            case GDK_KEY_F5:
                openThm->openDefaultViewer(3);
                return true;
            }
        } //if (!ctrl)
    } //if (!alt)

    if (alt) {
        switch (event->keyval) {
        case GDK_KEY_s:
            history->addBookmarkPressed ();
            setProgressStr(M("PROGRESSBAR_SNAPSHOT_ADDED"));
            return true;

        case GDK_KEY_f:
            iareapanel->imageArea->zoomPanel->zoomFitCropClicked();
            return true;
        }
    }

    if (shift) {
        switch (event->keyval) {
        case GDK_KEY_F3: // open Previous image from Editor's perspective
            if (!simpleEditor && fPanel && !fname.empty()) {
                EditorPanel::openPreviousEditorImage();
                return true;
            }

            break; // to avoid gcc complain

        case GDK_KEY_F4: // open next image from Editor's perspective
            if (!simpleEditor && fPanel && !fname.empty()) {
                EditorPanel::openNextEditorImage();
                return true;
            }

            break; // to avoid gcc complain
        }
    }

    if(tpc->getToolBar() && tpc->getToolBar()->handleShortcutKey(event)) {
        return true;
    }

    if(tpc->handleShortcutKey(event)) {
        return true;
    }

    if (!simpleEditor && fPanel) {
        if (fPanel->handleShortcutKey(event)) {
            return true;
        }
    }

    return false;
}

void EditorPanel::procParamsChanged (Thumbnail* thm, int whoChangedIt)
{

    if (whoChangedIt != EDITOR) {
        PartialProfile pp(true);
        pp.set(true);
        *(pp.pparams) = openThm->getProcParams();
        tpc->profileChange (&pp, rtengine::EvProfileChangeNotification, M("PROGRESSDLG_PROFILECHANGEDINBROWSER"));
        pp.deleteInstance();
    }
}

bool EditorPanel::idle_saveImage (ProgressConnector<rtengine::IImage16*> *pc, Glib::ustring fname, SaveFormat sf)
{
    rtengine::IImage16* img = pc->returnValue();
    delete pc;

    if( img ) {
        setProgressStr(M("GENERAL_SAVE"));
        setProgress(0.9f);

        ProgressConnector<int> *ld = new ProgressConnector<int>();
        img->setSaveProgressListener (parent->getProgressListener());

        if (sf.format == "tif")
            ld->startFunc (sigc::bind(sigc::mem_fun(img, &rtengine::IImage16::saveAsTIFF), fname, sf.tiffBits, sf.tiffUncompressed),
                           sigc::bind(sigc::mem_fun(*this, &EditorPanel::idle_imageSaved), ld, img, fname, sf));
        else if (sf.format == "png")
            ld->startFunc (sigc::bind(sigc::mem_fun(img, &rtengine::IImage16::saveAsPNG), fname, sf.pngCompression, sf.pngBits),
                           sigc::bind(sigc::mem_fun(*this, &EditorPanel::idle_imageSaved), ld, img, fname, sf));
        else if (sf.format == "jpg")
            ld->startFunc (sigc::bind(sigc::mem_fun(img, &rtengine::IImage16::saveAsJPEG), fname, sf.jpegQuality, sf.jpegSubSamp),
                           sigc::bind(sigc::mem_fun(*this, &EditorPanel::idle_imageSaved), ld, img, fname, sf));
    } else {
        Glib::ustring msg_ = Glib::ustring("<b>") + fname + ": Error during image processing\n</b>";
        Gtk::MessageDialog msgd (*parent, msg_, true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
        msgd.run ();

        saveimgas->set_sensitive(true);
        sendtogimp->set_sensitive(true);
        isProcessing = false;

    }

    rtengine::ImageSource* imgsrc = isrc->getImageSource ();
    imgsrc->setProgressListener(this);
    return false;
}

bool EditorPanel::idle_imageSaved(ProgressConnector<int> *pc, rtengine::IImage16* img, Glib::ustring fname, SaveFormat sf)
{
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
        error(M("MAIN_MSG_CANNOTSAVE"), fname);
    }

    saveimgas->set_sensitive(true);
    sendtogimp->set_sensitive(true);

    parent->setProgressStr("");
    parent->setProgress(0.);

    setProgressState(false);

    delete pc;
    SoundManager::playSoundAsync(options.sndBatchQueueDone);
    isProcessing = false;
    return false;
}

BatchQueueEntry* EditorPanel::createBatchQueueEntry ()
{

    rtengine::procparams::ProcParams pparams;
    ipc->getParams (&pparams);
    //rtengine::ProcessingJob* job = rtengine::ProcessingJob::create (ipc->getInitialImage(), pparams);
    rtengine::ProcessingJob* job = rtengine::ProcessingJob::create (openThm->getFileName (), openThm->getType() == FT_Raw, pparams);
    int fullW = 0, fullH = 0;
    isrc->getImageSource()->getFullSize(fullW, fullH, pparams.coarse.rotate == 90 || pparams.coarse.rotate == 270 ? TR_R90 : TR_NONE);
    int prevh = BatchQueue::calcMaxThumbnailHeight();
    int prevw = int((size_t)fullW * (size_t)prevh / (size_t)fullH);
    return new BatchQueueEntry (job, pparams, openThm->getFileName(), prevw, prevh, openThm);
}



void EditorPanel::saveAsPressed ()
{
    if (!ipc || !openThm) {
        return;
    }

    bool fnameOK = false;
    Glib::ustring fnameOut;

    SaveAsDialog* saveAsDialog;

    if (safe_file_test (options.lastSaveAsPath, Glib::FILE_TEST_IS_DIR)) {
        saveAsDialog = new SaveAsDialog (options.lastSaveAsPath);
    } else {
        saveAsDialog = new SaveAsDialog (safe_get_user_picture_dir());
    }

    saveAsDialog->set_default_size (options.saveAsDialogWidth, options.saveAsDialogHeight);
    saveAsDialog->setInitialFileName (lastSaveAsFileName);
    saveAsDialog->setImagePath (fname);

    do {
        int result = saveAsDialog->run ();

        // The SaveAsDialog ensure that a filename has been specified
        fnameOut = saveAsDialog->getFileName ();

        options.lastSaveAsPath = saveAsDialog->getDirectory ();
        options.saveAsDialogWidth = saveAsDialog->get_width ();
        options.saveAsDialogHeight = saveAsDialog->get_height ();
        options.autoSuffix = saveAsDialog->getAutoSuffix ();
        options.saveMethodNum = saveAsDialog->getSaveMethodNum ();
        lastSaveAsFileName = Glib::path_get_basename (removeExtension (fnameOut));
        SaveFormat sf = saveAsDialog->getFormat ();
        options.saveFormat = sf;
        options.forceFormatOpts = saveAsDialog->getForceFormatOpts ();

        if (result != Gtk::RESPONSE_OK) {
            break;
        }

        if (saveAsDialog->getImmediately ()) {
            // separate filename and the path to the destination directory
            Glib::ustring dstdir = Glib::path_get_dirname (fnameOut);
            Glib::ustring dstfname = Glib::path_get_basename (removeExtension(fnameOut));
            Glib::ustring dstext = getExtension (fnameOut);

            if (saveAsDialog->getAutoSuffix()) {

                Glib::ustring fnameTemp;

                for (int tries = 0; tries < 100; tries++) {
                    if (tries == 0) {
                        fnameTemp = Glib::ustring::compose ("%1.%2", Glib::build_filename (dstdir,  dstfname), dstext);
                    } else {
                        fnameTemp = Glib::ustring::compose ("%1-%2.%3", Glib::build_filename (dstdir,  dstfname), tries, dstext);
                    }

                    if (!safe_file_test (fnameTemp, Glib::FILE_TEST_EXISTS)) {
                        fnameOut = fnameTemp;
                        fnameOK = true;
                        break;
                    }
                }
            }

            // check if it exists
            if (!fnameOK) {
                fnameOK = confirmOverwrite (*saveAsDialog, fnameOut);
            }

            if (fnameOK) {
                isProcessing = true;
                // save image
                rtengine::procparams::ProcParams pparams;
                ipc->getParams (&pparams);
                rtengine::ProcessingJob* job = rtengine::ProcessingJob::create (ipc->getInitialImage(), pparams);

                ProgressConnector<rtengine::IImage16*> *ld = new ProgressConnector<rtengine::IImage16*>();
                ld->startFunc(sigc::bind(sigc::ptr_fun(&rtengine::processImage), job, err, parent->getProgressListener(), options.tunnelMetaData, false ),
                              sigc::bind(sigc::mem_fun( *this, &EditorPanel::idle_saveImage ), ld, fnameOut, sf ));
                saveimgas->set_sensitive(false);
                sendtogimp->set_sensitive(false);
            }
        } else {
            BatchQueueEntry* bqe = createBatchQueueEntry ();
            bqe->outFileName = fnameOut;
            bqe->saveFormat = saveAsDialog->getFormat ();
            bqe->forceFormatOpts = saveAsDialog->getForceFormatOpts ();
            parent->addBatchQueueJob (bqe, saveAsDialog->getToHeadOfQueue ());
            fnameOK = true;
        }

        // ask parent to redraw file browser
        // ... or does it automatically when the tab is switched to it
    } while (!fnameOK);

    saveAsDialog->hide();
}

void EditorPanel::queueImgPressed ()
{
    if (!ipc || !openThm) {
        return;
    }

    saveProfile ();
    parent->addBatchQueueJob (createBatchQueueEntry ());
}

void EditorPanel::sendToGimpPressed ()
{
    if (!ipc || !openThm) {
        return;
    }

    // develop image
    rtengine::procparams::ProcParams pparams;
    ipc->getParams (&pparams);
    rtengine::ProcessingJob* job = rtengine::ProcessingJob::create (ipc->getInitialImage(), pparams);
    ProgressConnector<rtengine::IImage16*> *ld = new ProgressConnector<rtengine::IImage16*>();
    ld->startFunc(sigc::bind(sigc::ptr_fun(&rtengine::processImage), job, err, parent->getProgressListener(), options.tunnelMetaData, false ),
                  sigc::bind(sigc::mem_fun( *this, &EditorPanel::idle_sendToGimp ), ld, openThm->getFileName() ));
    saveimgas->set_sensitive(false);
    sendtogimp->set_sensitive(false);
}


void EditorPanel::openPreviousEditorImage()
{
    if (!simpleEditor && fPanel && !fname.empty()) {
        fPanel->fileCatalog->openNextPreviousEditorImage(fname, false, NAV_PREVIOUS);
    }
}

void EditorPanel::openNextEditorImage()
{
    if (!simpleEditor && fPanel && !fname.empty()) {
        fPanel->fileCatalog->openNextPreviousEditorImage(fname, false, NAV_NEXT);
    }
}

void EditorPanel::syncFileBrowser()   // synchronize filebrowser with image in Editor
{
    if (!simpleEditor && fPanel && !fname.empty()) {
        fPanel->fileCatalog->selectImage(fname, false);
    }
}

bool EditorPanel::idle_sendToGimp( ProgressConnector<rtengine::IImage16*> *pc, Glib::ustring fname)
{

    rtengine::IImage16* img = pc->returnValue();
    delete pc;

    if (img) {
        // get file name base
        Glib::ustring shortname = removeExtension (Glib::path_get_basename (fname));
        Glib::ustring dirname = Glib::get_tmp_dir ();
        Glib::ustring fname = Glib::build_filename (dirname, shortname);

        SaveFormat sf;
        sf.format = "tif";
        sf.tiffBits = 16;
        sf.tiffUncompressed = true;
        sf.saveParams = true;

        Glib::ustring fileName = Glib::ustring::compose ("%1.%2", fname, sf.format);

        int tries = 1;

        while (safe_file_test (fileName, Glib::FILE_TEST_EXISTS) && tries < 1000) {
            fileName = Glib::ustring::compose("%1-%2.%3", fname, tries, sf.format);
            tries++;
        }

        if (tries == 1000) {
            img->free ();
            return false;
        }

        ProgressConnector<int> *ld = new ProgressConnector<int>();
        img->setSaveProgressListener (parent->getProgressListener());
        ld->startFunc (sigc::bind(sigc::mem_fun(img, &rtengine::IImage16::saveAsTIFF), fileName, sf.tiffBits, sf.tiffUncompressed),
                       sigc::bind(sigc::mem_fun(*this, &EditorPanel::idle_sentToGimp), ld, img, fileName));
    } else {
        Glib::ustring msg_ = Glib::ustring("<b> Error during image processing\n</b>");
        Gtk::MessageDialog msgd (*parent, msg_, true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
        msgd.run ();
        saveimgas->set_sensitive(true);
        sendtogimp->set_sensitive(true);
    }

    return false;
}

bool EditorPanel::idle_sentToGimp(ProgressConnector<int> *pc, rtengine::IImage16* img, Glib::ustring filename)
{
    img->free ();
    int errore = pc->returnValue();
    delete pc;

    if (!errore) {
        saveimgas->set_sensitive(true);
        sendtogimp->set_sensitive(true);
        parent->setProgressStr("");
        parent->setProgress(0.);
        bool success = false;
        Glib::ustring cmdLine;
        Glib::ustring executable;

        // start gimp
        if (options.editorToSendTo == 1) {
#ifdef WIN32
            executable = Glib::build_filename (Glib::build_filename(options.gimpDir, "bin"), "gimp-win-remote");
            cmdLine = Glib::ustring("\"") + executable + Glib::ustring("\" gimp-2.4.exe ") + Glib::ustring("\"") + filename + Glib::ustring("\"");

            if ( safe_file_test(executable, (Glib::FILE_TEST_EXISTS | Glib::FILE_TEST_IS_EXECUTABLE)) ) {
                success = safe_spawn_command_line_async (cmdLine);
            }

#elif defined __APPLE__
            cmdLine = Glib::ustring("open -a /Applications/GIMP.app \'") + filename + Glib::ustring("\'");
            success = safe_spawn_command_line_async (cmdLine);
            std::cout << cmdLine << std::endl;
#else
            cmdLine = Glib::ustring("gimp \"") + filename + Glib::ustring("\"");
            success = safe_spawn_command_line_async (cmdLine);
            std::cout << cmdLine << std::endl;
#endif

            if (!success) {
#ifdef WIN32
                int ver = 12;

                while (!success && ver) {
                    executable = Glib::build_filename (Glib::build_filename(options.gimpDir, "bin"), Glib::ustring::compose(Glib::ustring("gimp-2.%1.exe"), ver));

                    if ( safe_file_test(executable, (Glib::FILE_TEST_EXISTS | Glib::FILE_TEST_IS_EXECUTABLE)) ) {
                        cmdLine = Glib::ustring("\"") + executable + Glib::ustring("\" \"") + filename + Glib::ustring("\"");
                        success = safe_spawn_command_line_async (cmdLine);
                    }

                    ver--;
                }

#elif defined __APPLE__
                cmdLine = Glib::ustring("open -a /Applications/Gimp.app/Contents/Resources/start \'") + filename + Glib::ustring("\'");
                success = safe_spawn_command_line_async (cmdLine);
                std::cout << cmdLine << std::endl;
#else
                cmdLine = Glib::ustring("gimp-remote \"") + filename + Glib::ustring("\"");
                success = safe_spawn_command_line_async (cmdLine);
                std::cout << cmdLine << std::endl;
#endif
            }
        } else if (options.editorToSendTo == 2) {
#ifdef WIN32
            executable = Glib::build_filename(options.psDir, "Photoshop.exe");

            if ( safe_file_test(executable, (Glib::FILE_TEST_EXISTS | Glib::FILE_TEST_IS_EXECUTABLE)) ) {
                cmdLine = Glib::ustring("\"") + executable + Glib::ustring("\" \"") + filename + Glib::ustring("\"");
                success = safe_spawn_command_line_async (cmdLine);
            }

#else
#ifdef __APPLE__
            cmdLine = Glib::ustring("open -a \'") + Glib::build_filename(options.psDir, "Photoshop.app\' ")  + Glib::ustring("\'") + filename + Glib::ustring("\'");
#else
            cmdLine = Glib::ustring("\"") + Glib::build_filename(options.psDir, "Photoshop.exe") + Glib::ustring("\" \"") + filename + Glib::ustring("\"");
#endif
            success = safe_spawn_command_line_async (cmdLine);
            std::cout << cmdLine << std::endl;
#endif
        } else if (options.editorToSendTo == 3) {
#ifdef WIN32

            if ( safe_file_test(options.customEditorProg, (Glib::FILE_TEST_EXISTS | Glib::FILE_TEST_IS_EXECUTABLE)) ) {
                cmdLine = Glib::ustring("\"") + options.customEditorProg + Glib::ustring("\" \"") + filename + Glib::ustring("\"");
                success = safe_spawn_command_line_async (cmdLine);
            }

#else
#ifdef __APPLE__
            cmdLine = options.customEditorProg + Glib::ustring(" \"") + filename + Glib::ustring("\"");
#else
            cmdLine = Glib::ustring("\"") + options.customEditorProg + Glib::ustring("\" \"") + filename + Glib::ustring("\"");
#endif
            success = safe_spawn_command_line_async (cmdLine);
            std::cout << cmdLine << std::endl;
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

void EditorPanel::historyBeforeLineChanged (const rtengine::procparams::ProcParams& params)
{

    if (beforeIpc) {
        ProcParams* pparams = beforeIpc->beginUpdateParams ();
        *pparams = params;
        beforeIpc->endUpdateParams (rtengine::EvProfileChanged);  // starts the IPC processing
    }
}

void EditorPanel::beforeAfterToggled ()
{

    if(!ipc) {
        return;
    }

    removeIfThere (beforeAfterBox,  beforeBox, false);
    removeIfThere (afterBox,  afterHeaderBox, false);

    if (beforeIarea) {
        if (beforeIpc) {
            beforeIpc->stopProcessing ();
        }

        iareapanel->setBeforeAfterViews (NULL, iareapanel);
        iareapanel->imageArea->iLinkedImageArea = NULL;
        delete beforeIarea;
        beforeIarea = NULL;

        if (beforeIpc) {
            beforeIpc->setPreviewImageListener (NULL);
        }

        delete beforePreviewHandler;
        beforePreviewHandler = NULL;

        if (beforeIpc) {
            rtengine::StagedImageProcessor::destroy (beforeIpc);
        }

        beforeIpc = NULL;
    }

    if (beforeAfter->get_active ()) {

        int errorCode = 0;
        rtengine::InitialImage *beforeImg = rtengine::InitialImage::load ( isrc->getImageSource ()->getFileName(),  openThm->getType() == FT_Raw , &errorCode, NULL);

        if( !beforeImg || errorCode ) {
            return;
        }

        beforeIarea = new ImageAreaPanel ();

        int HeaderBoxHeight = 17;

        beforeLabel = Gtk::manage (new Gtk::Label ());
        beforeLabel->set_markup (Glib::ustring("<b>") + M("GENERAL_BEFORE") + "</b>");
        tbBeforeLock = Gtk::manage (new Gtk::ToggleButton ());
        tbBeforeLock->set_tooltip_markup (M("MAIN_TOOLTIP_BEFOREAFTERLOCK"));
        tbBeforeLock->signal_toggled().connect( sigc::mem_fun(*this, &EditorPanel::tbBeforeLock_toggled) );
        beforeHeaderBox = Gtk::manage (new Gtk::HBox ());
        beforeHeaderBox->pack_end (*tbBeforeLock, Gtk::PACK_SHRINK, 2);
        beforeHeaderBox->pack_end (*beforeLabel, Gtk::PACK_SHRINK, 2);
        beforeHeaderBox->set_size_request(0, HeaderBoxHeight);

        history->blistenerLock ? tbBeforeLock->set_image (*iBeforeLockON) : tbBeforeLock->set_image (*iBeforeLockOFF);
        tbBeforeLock->set_active(history->blistenerLock);

        beforeBox = Gtk::manage (new Gtk::VBox ());
        beforeBox->pack_start (*beforeHeaderBox, Gtk::PACK_SHRINK, 2);
        beforeBox->pack_start (*beforeIarea);

        afterLabel = Gtk::manage (new Gtk::Label ());
        afterLabel->set_markup (Glib::ustring("<b>") + M("GENERAL_AFTER") + "</b>");
        afterHeaderBox = Gtk::manage (new Gtk::HBox ());
        afterHeaderBox->set_size_request(0, HeaderBoxHeight);
        afterHeaderBox->pack_end (*afterLabel, Gtk::PACK_SHRINK, 2);
        afterBox->pack_start (*afterHeaderBox, Gtk::PACK_SHRINK, 2);
        afterBox->reorder_child (*afterHeaderBox, 0);

        beforeAfterBox->pack_start (*beforeBox);
        beforeAfterBox->reorder_child (*beforeBox, 0);
        beforeAfterBox->show_all ();

        beforePreviewHandler = new PreviewHandler ();

        beforeIpc = rtengine::StagedImageProcessor::create (beforeImg);
        beforeIpc->setPreviewScale (10);
        beforeIpc->setPreviewImageListener (beforePreviewHandler);
        beforeIarea->imageArea->setPreviewHandler (beforePreviewHandler);
        beforeIarea->imageArea->setImProcCoordinator (beforeIpc);

        beforeIarea->imageArea->setPreviewModePanel(iareapanel->imageArea->previewModePanel);
        beforeIarea->imageArea->setIndicateClippedPanel(iareapanel->imageArea->indClippedPanel);
        iareapanel->imageArea->iLinkedImageArea = beforeIarea->imageArea;

        iareapanel->setBeforeAfterViews (beforeIarea, iareapanel);
        beforeIarea->setBeforeAfterViews (beforeIarea, iareapanel);

        rtengine::procparams::ProcParams params;

        if (history->getBeforeLineParams (params)) {
            historyBeforeLineChanged (params);
        }
    }
}

void EditorPanel::tbBeforeLock_toggled ()
{
    history->blistenerLock = tbBeforeLock->get_active();
    tbBeforeLock->get_active() ? tbBeforeLock->set_image (*iBeforeLockON) : tbBeforeLock->set_image (*iBeforeLockOFF);
}

void EditorPanel::histogramChanged (LUTu & histRed, LUTu & histGreen, LUTu & histBlue, LUTu & histLuma, LUTu & histToneCurve, LUTu & histLCurve, LUTu & histCCurve, /*LUTu & histCLurve, LUTu & histLLCurve,*/ LUTu & histLCAM, LUTu & histCCAM,
                                    LUTu & histRedRaw, LUTu & histGreenRaw, LUTu & histBlueRaw , LUTu & histChroma)
{

    if (histogramPanel) {
        histogramPanel->histogramChanged (histRed, histGreen, histBlue, histLuma, histRedRaw, histGreenRaw, histBlueRaw, histChroma);
    }

    tpc->updateCurveBackgroundHistogram (histToneCurve, histLCurve, histCCurve,/*histCLurve,  histLLCurve,*/ histLCAM, histCCAM, histRed, histGreen, histBlue, histLuma);
}

bool EditorPanel::CheckSidePanelsVisibility()
{
    if (tbTopPanel_1) {
        if(tbTopPanel_1->get_active() == false && tbRightPanel_1->get_active() == false && hidehp->get_active() == false) {
            return false;
        }

        return true;
    }

    if(tbRightPanel_1->get_active() == false && hidehp->get_active() == false) {
        return false;
    }

    return true;
}
void EditorPanel::toggleSidePanels()
{
    // Maximize preview panel:
    // toggle top AND right AND history panels

    bool bAllSidePanelsVisible;
    bAllSidePanelsVisible = CheckSidePanelsVisibility();

    if (tbTopPanel_1) {
        tbTopPanel_1->set_active (!bAllSidePanelsVisible);
    }

    tbRightPanel_1->set_active (!bAllSidePanelsVisible);
    hidehp->set_active (!bAllSidePanelsVisible);

    if (bAllSidePanelsVisible == false) {
        tbShowHideSidePanels->set_image (*iShowHideSidePanels);
    } else {
        tbShowHideSidePanels->set_image (*iShowHideSidePanels_exit);
    }
}

void EditorPanel::toggleSidePanelsZoomFit()
{
    toggleSidePanels();

    // fit image preview
    // !!! TODO this does not want to work... seems to have an effect on a subsequent key press
    // iarea->imageArea->zoomPanel->zoomFitClicked();
}

void EditorPanel::tbShowHideSidePanels_managestate()
{
    bool bAllSidePanelsVisible;
    bAllSidePanelsVisible = CheckSidePanelsVisibility();
    ShowHideSidePanelsconn.block (true);

    tbShowHideSidePanels->set_active (!bAllSidePanelsVisible);

    ShowHideSidePanelsconn.block (false);
}

void EditorPanel::updateTPVScrollbar (bool hide)
{
    tpc->updateTPVScrollbar (hide);
}

void EditorPanel::updateTabsUsesIcons (bool useIcons)
{
    tpc->updateTabsUsesIcons (useIcons);
}

void EditorPanel::updateHistogramPosition (int oldPosition, int newPosition)
{

    switch (newPosition) {
    case 0:

        // No histogram
        if (!oldPosition) {
            // An histogram actually exist, we delete it
            if      (oldPosition == 1) {
                removeIfThere(leftbox, histogramPanel, false);
            } else if (oldPosition == 2) {
                removeIfThere(vboxright, histogramPanel, false);
            }

            delete histogramPanel;
            histogramPanel = NULL;
        }

        // else no need to create it
        break;

    case 1:

        // Histogram on the left pane
        if (oldPosition == 0) {
            // There was no Histogram before, so we create it
            histogramPanel = Gtk::manage (new HistogramPanel ());
            leftbox->pack_start (*histogramPanel, Gtk::PACK_SHRINK, 2);
        } else if (oldPosition == 2) {
            // The histogram was on the right side, so we move it to the left
            histogramPanel->reference();
            removeIfThere(vboxright, histogramPanel, false);
            leftbox->pack_start (*histogramPanel, Gtk::PACK_SHRINK, 2);
            histogramPanel->unreference();
        }

        histogramPanel->reorder(Gtk::POS_LEFT);
        leftbox->reorder_child(*histogramPanel, 0);
        break;

    case 2:
    default:

        // Histogram on the right pane
        if (oldPosition == 0) {
            // There was no Histogram before, so we create it
            histogramPanel = Gtk::manage (new HistogramPanel ());
            vboxright->pack_start (*histogramPanel, Gtk::PACK_SHRINK, 2);
        } else if (oldPosition == 1) {
            // The histogram was on the left side, so we move it to the right
            histogramPanel->reference();
            removeIfThere(leftbox, histogramPanel, false);
            vboxright->pack_start (*histogramPanel, Gtk::PACK_SHRINK, 2);
            histogramPanel->unreference();
        }

        histogramPanel->reorder(Gtk::POS_RIGHT);
        vboxright->reorder_child(*histogramPanel, 0);
        break;
    }

    iareapanel->imageArea->setPointerMotionHListener (histogramPanel);
}
