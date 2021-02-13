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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "batchqueuepanel.h"
#include "options.h"
#include "multilangmgr.h"
#include "rtwindow.h"
#include "soundman.h"
#include "rtimage.h"

static Glib::ustring makeFolderLabel(Glib::ustring path)
{
    if (!Glib::file_test (path, Glib::FILE_TEST_IS_DIR)) {
        return "(" + M("GENERAL_NONE") + ")";
    }

    if (path.size() > 40) {
        size_t last_ds = path.find_last_of (G_DIR_SEPARATOR);

        if (last_ds != Glib::ustring::npos && last_ds > 10) {
            path = "..." + path.substr(last_ds);
        }
    }

    return path;
}

BatchQueuePanel::BatchQueuePanel (FileCatalog* aFileCatalog) : parent(nullptr)
{
    set_orientation(Gtk::ORIENTATION_VERTICAL);
    batchQueue = Gtk::manage( new BatchQueue(aFileCatalog) );

    Gtk::Box* batchQueueButtonBox = Gtk::manage (new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    batchQueueButtonBox->set_name("BatchQueueButtons");

    qStartStop = Gtk::manage (new Gtk::Switch());
    qStartStop->set_tooltip_markup (M("QUEUE_STARTSTOP_TOOLTIP"));
    qStartStopConn = qStartStop->property_active().signal_changed().connect (sigc::mem_fun(*this, &BatchQueuePanel::startOrStopBatchProc));

    qAutoStart = Gtk::manage (new Gtk::CheckButton (M("QUEUE_AUTOSTART")));
    qAutoStart->set_tooltip_text (M("QUEUE_AUTOSTART_TOOLTIP"));
    qAutoStart->set_active (options.procQueueEnabled);

    queueShouldRun = false;

    batchQueueButtonBox->pack_start (*qStartStop, Gtk::PACK_SHRINK, 4);
    batchQueueButtonBox->pack_start (*qAutoStart, Gtk::PACK_SHRINK, 4);
    Gtk::Frame *bbox = Gtk::manage(new Gtk::Frame(M("MAIN_FRAME_QUEUE")));
    bbox->set_label_align(0.025, 0.5);
    bbox->add(*batchQueueButtonBox);

    // Output directory selection
    fdir = Gtk::manage (new Gtk::Frame (M("QUEUE_LOCATION_TITLE")));
    fdir->set_label_align(0.025, 0.5);
    Gtk::Box* odvb = Gtk::manage (new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    Gtk::Box* hb2 = Gtk::manage (new Gtk::Box ());
    useTemplate = Gtk::manage (new Gtk::RadioButton (M("QUEUE_LOCATION_TEMPLATE") + ":"));
    hb2->pack_start (*useTemplate, Gtk::PACK_SHRINK, 4);
    outdirTemplate = Gtk::manage (new Gtk::Entry ());
    hb2->pack_start (*outdirTemplate);
    odvb->pack_start (*hb2, Gtk::PACK_SHRINK, 4);
    outdirTemplate->set_tooltip_markup (M("QUEUE_LOCATION_TEMPLATE_TOOLTIP"));
    useTemplate->set_tooltip_markup (M("QUEUE_LOCATION_TEMPLATE_TOOLTIP"));
    Gtk::Box* hb3 = Gtk::manage (new Gtk::Box ());
    useFolder = Gtk::manage (new Gtk::RadioButton (M("QUEUE_LOCATION_FOLDER") + ":"));
    hb3->pack_start (*useFolder, Gtk::PACK_SHRINK, 4);

#if 0 //defined(__APPLE__) || defined(__linux__)
    // At the time of writing (2013-11-11) the gtkmm FileChooserButton with ACTION_SELECT_FOLDER
    // is so buggy on these platforms (OS X and Linux) that we rather employ this ugly button hack.
    // When/if GTKMM gets fixed we can go back to use the FileChooserButton, like we do on Windows.
    outdirFolderButton = Gtk::manage (new Gtk::Button("(" + M("GENERAL_NONE") + ")"));
    outdirFolderButton->set_alignment(0.0, 0.0);
    hb3->pack_start (*outdirFolderButton);
    outdirFolderButton->signal_pressed().connect( sigc::mem_fun(*this, &BatchQueuePanel::pathFolderButtonPressed) );
    outdirFolderButton->set_label(makeFolderLabel(options.savePathFolder));
    Gtk::Image* folderImg = Gtk::manage (new RTImage ("folder-closed.png"));
    folderImg->show ();
    outdirFolderButton->set_image (*folderImg);
    outdirFolder = nullptr;
#else
    outdirFolder = Gtk::manage (new MyFileChooserButton (M("QUEUE_LOCATION_FOLDER"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    hb3->pack_start (*outdirFolder);
    outdirFolder->signal_selection_changed().connect (sigc::mem_fun(*this, &BatchQueuePanel::pathFolderChanged));

    if (Glib::file_test (options.savePathFolder, Glib::FILE_TEST_IS_DIR)) {
        outdirFolder->set_current_folder (options.savePathFolder);
    } else {
        outdirFolder->set_current_folder (Glib::get_home_dir());
    }

    outdirFolderButton = 0;
#endif

    odvb->pack_start (*hb3, Gtk::PACK_SHRINK, 4);
    Gtk::RadioButton::Group g = useTemplate->get_group();
    useFolder->set_group (g);
    fdir->add (*odvb);

    // Output file format selection
    fformat = Gtk::manage (new Gtk::Frame (M("QUEUE_FORMAT_TITLE")));
    fformat->set_label_align(0.025, 0.5);
    saveFormatPanel = Gtk::manage (new SaveFormatPanel ());
    setExpandAlignProperties(saveFormatPanel, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    fformat->add (*saveFormatPanel);

    outdirTemplate->set_text (options.savePathTemplate);
    useTemplate->set_active (options.saveUsePathTemplate);
    useFolder->set_active (!options.saveUsePathTemplate);

    // setup signal handlers
    outdirTemplate->signal_changed().connect (sigc::mem_fun(*this, &BatchQueuePanel::saveOptions));
    useTemplate->signal_toggled().connect (sigc::mem_fun(*this, &BatchQueuePanel::saveOptions));
    useFolder->signal_toggled().connect (sigc::mem_fun(*this, &BatchQueuePanel::saveOptions));
    saveFormatPanel->setListener (this);

    // setup button bar
    topBox = Gtk::manage (new Gtk::Box ());
    pack_start (*topBox, Gtk::PACK_SHRINK);
    topBox->set_name("BatchQueueButtonsMainContainer");

    topBox->pack_start (*bbox, Gtk::PACK_SHRINK, 4);
    topBox->pack_start (*fdir, Gtk::PACK_EXPAND_WIDGET, 4);
    topBox->pack_start (*fformat, Gtk::PACK_EXPAND_WIDGET, 4);

    // add middle browser area
    pack_start (*batchQueue);

    // lower box with thumbnail zoom
    bottomBox = Gtk::manage (new Gtk::Box ());
    pack_start (*bottomBox, Gtk::PACK_SHRINK);

    // thumbnail zoom
    Gtk::Box* zoomBox = Gtk::manage (new Gtk::Box ());
    zoomBox->pack_start (*Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_VERTICAL)), Gtk::PACK_SHRINK, 4);
    Gtk::Label* zoomLabel = Gtk::manage (new Gtk::Label (Glib::ustring("<b>") + M("FILEBROWSER_THUMBSIZE") + ":</b>"));
    zoomLabel->set_use_markup (true);
    zoomBox->pack_start (*zoomLabel, Gtk::PACK_SHRINK, 4);
    zoomInButton  = Gtk::manage (new Gtk::Button ());
    zoomInButton->set_image (*Gtk::manage (new RTImage ("magnifier-plus.png")));
    zoomInButton->signal_pressed().connect (sigc::mem_fun(*batchQueue, &BatchQueue::zoomIn));
    zoomInButton->set_relief (Gtk::RELIEF_NONE);
    zoomInButton->set_tooltip_markup (M("FILEBROWSER_ZOOMINHINT"));
    zoomBox->pack_end (*zoomInButton, Gtk::PACK_SHRINK);
    zoomOutButton  = Gtk::manage (new Gtk::Button ());
    zoomOutButton->set_image (*Gtk::manage (new RTImage ("magnifier-minus.png")));
    zoomOutButton->signal_pressed().connect (sigc::mem_fun(*batchQueue, &BatchQueue::zoomOut));
    zoomOutButton->set_relief (Gtk::RELIEF_NONE);
    zoomOutButton->set_tooltip_markup (M("FILEBROWSER_ZOOMOUTHINT"));
    zoomBox->pack_end (*zoomOutButton, Gtk::PACK_SHRINK);
    bottomBox->pack_end (*zoomBox, Gtk::PACK_SHRINK);


    batchQueue->setBatchQueueListener (this);

    show_all ();

    if (batchQueue->loadBatchQueue()) {
        idle_register.add(
            [this]() -> bool
            {
                batchQueue->resizeLoadedQueue();
                return false;
            },
            G_PRIORITY_LOW
        );
    }
}

BatchQueuePanel::~BatchQueuePanel()
{
    idle_register.destroy();
}

void BatchQueuePanel::init (RTWindow *parent)
{
    this->parent = parent;

    saveFormatPanel->init (options.saveFormatBatch);
}

// it is expected to have a non null forceOrientation value on Preferences update only. In this case, qsize is ignored and computed automatically
void BatchQueuePanel::updateTab (int qsize, int forceOrientation)
{
    Gtk::Notebook *nb = (Gtk::Notebook *)(this->get_parent());

    if (forceOrientation > 0) {
        qsize = batchQueue->getEntries().size();
    }

    Gtk::Grid* grid = Gtk::manage (new Gtk::Grid ());
    if ((forceOrientation == 0 && options.mainNBVertical) || (forceOrientation == 2)) {
        Gtk::Label* l;

        if(!qsize ) {
            grid->attach_next_to(*Gtk::manage (new RTImage ("gears.png")), Gtk::POS_TOP, 1, 1);
            l = Gtk::manage (new Gtk::Label (Glib::ustring(" ") + M("MAIN_FRAME_QUEUE")) );
        } else if (qStartStop->get_active()) {
            grid->attach_next_to(*Gtk::manage (new RTImage ("gears-play.png")), Gtk::POS_TOP, 1, 1);
            l = Gtk::manage (new Gtk::Label (Glib::ustring(" ") + M("MAIN_FRAME_QUEUE") + " [" + Glib::ustring::format( qsize ) + "]"));
        } else {
            grid->attach_next_to(*Gtk::manage (new RTImage ("gears-pause.png")), Gtk::POS_TOP, 1, 1);
            l = Gtk::manage (new Gtk::Label (Glib::ustring(" ") + M("MAIN_FRAME_QUEUE") + " [" + Glib::ustring::format( qsize ) + "]" ));
        }

        l->set_angle (90);
        grid->attach_next_to(*l, Gtk::POS_TOP, 1, 1);
        grid->set_tooltip_markup (M("MAIN_FRAME_QUEUE_TOOLTIP"));
        grid->show_all ();

        if (nb) {
            nb->set_tab_label(*this, *grid);
        }
    } else {
        if (!qsize ) {
            grid->attach_next_to(*Gtk::manage (new RTImage ("gears.png")), Gtk::POS_RIGHT, 1, 1);
            grid->attach_next_to(*Gtk::manage (new Gtk::Label (M("MAIN_FRAME_QUEUE") )), Gtk::POS_RIGHT, 1, 1);
        } else if (qStartStop->get_active()) {
            grid->attach_next_to(*Gtk::manage (new RTImage ("gears-play.png")), Gtk::POS_RIGHT, 1, 1);
            grid->attach_next_to(*Gtk::manage (new Gtk::Label (M("MAIN_FRAME_QUEUE") + " [" + Glib::ustring::format( qsize ) + "]" )), Gtk::POS_RIGHT, 1, 1);
        } else {
            grid->attach_next_to(*Gtk::manage (new RTImage ("gears-pause.png")), Gtk::POS_RIGHT, 1, 1);
            grid->attach_next_to(*Gtk::manage (new Gtk::Label (M("MAIN_FRAME_QUEUE") + " [" + Glib::ustring::format( qsize ) + "]" )), Gtk::POS_RIGHT, 1, 1);
        }

        grid->set_tooltip_markup (M("MAIN_FRAME_QUEUE_TOOLTIP"));
        grid->show_all ();

        if (nb) {
            nb->set_tab_label(*this, *grid);
        }
    }
}

void BatchQueuePanel::queueSizeChanged(int qsize, bool queueRunning, bool queueError, const Glib::ustring& queueErrorMessage)
{
    setGuiFromBatchState(queueRunning, qsize);

    if (!queueRunning && qsize == 0 && queueShouldRun) {
        // There was work, but it is all done now.
        queueShouldRun = false;

        SoundManager::playSoundAsync(options.sndBatchQueueDone);
    }

    if (queueError) {
        Gtk::MessageDialog msgd (queueErrorMessage, true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
        msgd.run ();
    }
}

void BatchQueuePanel::startOrStopBatchProc()
{
    if (qStartStop->get_state()) {
        startBatchProc();
    } else {
        stopBatchProc();
    }
}

void BatchQueuePanel::startBatchProc ()
{
    if (batchQueue->hasJobs()) {
        // Update the *desired* state of the queue, then launch it.  The switch
        // state is not updated here; it is changed by the queueSizeChanged()
        // callback in response to the *reported* state.
        queueShouldRun = true;

        saveOptions();
        batchQueue->startProcessing ();
    }
}

void BatchQueuePanel::stopBatchProc ()
{
    // There is nothing much to do here except set the desired state, which the
    // background queue thread must check.  It will notify queueSizeChanged()
    // when it stops.
    queueShouldRun = false;
}

void BatchQueuePanel::setGuiFromBatchState(bool queueRunning, int qsize)
{
    // Change the GUI state in response to the reported queue state
    if (qsize == 0 || (qsize == 1 && queueRunning)) {
        qStartStop->set_sensitive(false);
    } else {
        qStartStop->set_sensitive(true);
    }

    qStartStopConn.block(true);
    qStartStop->set_active(queueRunning);
    qStartStopConn.block(false);

    fdir->set_sensitive (!queueRunning);
    fformat->set_sensitive (!queueRunning);

    updateTab(qsize);
}

void BatchQueuePanel::addBatchQueueJobs(const std::vector<BatchQueueEntry*>& entries, bool head)
{
    batchQueue->addEntries(entries, head);

    if (!qStartStop->get_active() && qAutoStart->get_active()) {
        // Auto-start as if the user had pressed the qStartStop switch
        startBatchProc ();
    }
}

void BatchQueuePanel::saveOptions ()
{

    options.savePathTemplate    = outdirTemplate->get_text();
    options.saveUsePathTemplate = useTemplate->get_active();
    options.procQueueEnabled    = qAutoStart->get_active();
}

bool BatchQueuePanel::handleShortcutKey (GdkEventKey* event)
{
    bool ctrl = event->state & GDK_CONTROL_MASK;

    if (ctrl) {
        switch(event->keyval) {
        case GDK_KEY_s:
            if (qStartStop->get_active()) {
                stopBatchProc();
            } else {
                startBatchProc();
            }

            return true;
        }
    }

    return batchQueue->keyPressed (event);
}

bool BatchQueuePanel::canStartNext ()
{
    // This function is called from the background BatchQueue thread.  It
    // cannot call UI functions; we keep the desired state in an atomic.
    return queueShouldRun;
}

void BatchQueuePanel::pathFolderButtonPressed ()
{

    Gtk::FileChooserDialog fc (getToplevelWindow (this), M("QUEUE_LOCATION_FOLDER"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER );
    fc.add_button( "_Cancel", Gtk::RESPONSE_CANCEL); // STOCKICON WAS THERE
    fc.add_button( "_OK", Gtk::RESPONSE_OK); // STOCKICON WAS THERE
    fc.set_filename(options.savePathFolder);
    fc.set_transient_for(*parent);
    int result = fc.run();

    if (result == Gtk::RESPONSE_OK) {
        if (Glib::file_test(fc.get_current_folder(), Glib::FILE_TEST_IS_DIR)) {
            options.savePathFolder = fc.get_filename ();
            outdirFolderButton->set_label(makeFolderLabel(options.savePathFolder));
        }
    }
}

// We only want to save the following when it changes,
// since these settings are shared with editorpanel :
void BatchQueuePanel::pathFolderChanged ()
{
    options.savePathFolder = outdirFolder->get_filename();
}

void BatchQueuePanel::formatChanged(const Glib::ustring& format)
{
    options.saveFormatBatch = saveFormatPanel->getFormat();
}
