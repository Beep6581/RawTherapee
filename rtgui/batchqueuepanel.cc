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
#include "batchqueuepanel.h"
#include "options.h"
#include "preferences.h"
#include "multilangmgr.h"
#include "rtwindow.h"
#include "soundman.h"
#include "rtimage.h"

struct BQProcessLoaded {
    BatchQueue* bq;
};

int processLoadedBatchQueueUIThread (void* data)
{

    BatchQueue* bq = static_cast<BatchQueue*>(data);
    bq->resizeLoadedQueue();
    return 0;
}

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

BatchQueuePanel::BatchQueuePanel (FileCatalog* aFileCatalog)
{

    batchQueue = Gtk::manage( new BatchQueue(aFileCatalog) );

    // construct batch queue panel with the extra "start" and "stop" button
    Gtk::VBox* batchQueueButtonBox = Gtk::manage (new Gtk::VBox);
    start = Gtk::manage (new Gtk::ToggleButton ());
    stop = Gtk::manage (new Gtk::ToggleButton ());
    autoStart = Gtk::manage (new Gtk::CheckButton (M("BATCHQUEUE_AUTOSTART")));
    start->set_tooltip_markup (M("FILEBROWSER_STARTPROCESSINGHINT"));
    stop->set_tooltip_markup (M("FILEBROWSER_STOPPROCESSINGHINT"));
    autoStart->set_tooltip_text (M("FILEBROWSER_TOOLTIP_STOPPROCESSING"));
    start->set_active (false);
    stop->set_active (true);
    autoStart->set_active (options.procQueueEnabled);

    start->set_image (*Gtk::manage (new RTImage ("gtk-media-play.png")));
    start->get_style_context()->add_class("BIG");
    startConnection = start->signal_toggled().connect (sigc::mem_fun(*this, &BatchQueuePanel::startBatchProc));
    stop->set_image (*Gtk::manage (new RTImage ("gtk-media-stop.png")));
    stop->get_style_context()->add_class("BIG");
    stopConnection = stop->signal_toggled().connect (sigc::mem_fun(*this, &BatchQueuePanel::stopBatchProc));
    batchQueueButtonBox->pack_start (*start, Gtk::PACK_SHRINK, 4);
    batchQueueButtonBox->pack_start (*stop, Gtk::PACK_SHRINK, 4);
    batchQueueButtonBox->pack_start (*autoStart, Gtk::PACK_SHRINK, 4);

    // Output directory selection
    fdir = Gtk::manage (new Gtk::Frame (M("PREFERENCES_OUTDIR")));
    Gtk::VBox* odvb = Gtk::manage (new Gtk::VBox ());
    Gtk::HBox* hb2 = Gtk::manage (new Gtk::HBox ());
    useTemplate = Gtk::manage (new Gtk::RadioButton (M("PREFERENCES_OUTDIRTEMPLATE") + ":"));
    hb2->pack_start (*useTemplate, Gtk::PACK_SHRINK, 4);
    outdirTemplate = Gtk::manage (new Gtk::Entry ());
    hb2->pack_start (*outdirTemplate);
    odvb->pack_start (*hb2, Gtk::PACK_SHRINK, 4);
    outdirTemplate->set_tooltip_markup (M("PREFERENCES_OUTDIRTEMPLATEHINT"));
    useTemplate->set_tooltip_markup (M("PREFERENCES_OUTDIRTEMPLATEHINT"));
    Gtk::HBox* hb3 = Gtk::manage (new Gtk::HBox ());
    useFolder = Gtk::manage (new Gtk::RadioButton (M("PREFERENCES_OUTDIRFOLDER") + ":"));
    hb3->pack_start (*useFolder, Gtk::PACK_SHRINK, 4);

#if 0 //defined(__APPLE__) || defined(__linux__)
    // At the time of writing (2013-11-11) the gtkmm FileChooserButton with ACTION_SELECT_FOLDER
    // is so buggy on these platforms (OS X and Linux) that we rather employ this ugly button hack.
    // When/if GTKMM gets fixed we can go back to use the FileChooserButton, like we do on Windows.
    outdirFolderButton = Gtk::manage (new Gtk::Button("(" + M("GENERAL_NONE") + ")"));
    outdirFolderButton->set_alignment(0.0, 0.0);
    hb3->pack_start (*outdirFolderButton);
    outdirFolderButton->signal_pressed().connect( sigc::mem_fun(*this, &BatchQueuePanel::pathFolderButtonPressed) );
    outdirFolderButton->set_tooltip_markup (M("PREFERENCES_OUTDIRFOLDERHINT"));
    outdirFolderButton->set_label(makeFolderLabel(options.savePathFolder));
    Gtk::Image* folderImg = Gtk::manage (new RTImage ("gtk-directory.png"));
    folderImg->show ();
    outdirFolderButton->set_image (*folderImg);
    outdirFolder = 0;
#else
    outdirFolder = Gtk::manage (new MyFileChooserButton (M("PREFERENCES_OUTDIRFOLDER"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    hb3->pack_start (*outdirFolder);
    outdirFolder->signal_current_folder_changed().connect (sigc::mem_fun(*this, &BatchQueuePanel::pathFolderChanged));
    outdirFolder->set_tooltip_markup (M("PREFERENCES_OUTDIRFOLDERHINT"));

    if (Glib::file_test (options.savePathFolder, Glib::FILE_TEST_IS_DIR)) {
        outdirFolder->set_current_folder (options.savePathFolder);
    }

    outdirFolderButton = 0;
#endif

    odvb->pack_start (*hb3, Gtk::PACK_SHRINK, 4);
    useFolder->set_tooltip_markup (M("PREFERENCES_OUTDIRFOLDERHINT"));
    Gtk::RadioButton::Group g = useTemplate->get_group();
    useFolder->set_group (g);
    fdir->add (*odvb);

    // Output file format selection
    fformat = Gtk::manage (new Gtk::Frame (M("PREFERENCES_FILEFORMAT")));
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
    topBox = Gtk::manage (new Gtk::HBox ());
    pack_start (*topBox, Gtk::PACK_SHRINK);

    topBox->pack_start (*batchQueueButtonBox, Gtk::PACK_SHRINK, 4);
    topBox->pack_start (*fdir, Gtk::PACK_EXPAND_WIDGET, 4);
    topBox->pack_start (*fformat, Gtk::PACK_EXPAND_WIDGET, 4);

    // add middle browser area
    pack_start (*batchQueue);

    // lower box with thumbnail zoom
    bottomBox = Gtk::manage (new Gtk::HBox ());
    pack_start (*bottomBox, Gtk::PACK_SHRINK);

    // thumbnail zoom
    Gtk::HBox* zoomBox = Gtk::manage (new Gtk::HBox ());
    zoomBox->pack_start (*Gtk::manage (new Gtk::VSeparator), Gtk::PACK_SHRINK, 4);
    Gtk::Label* zoomLabel = Gtk::manage (new Gtk::Label (Glib::ustring("<b>") + M("FILEBROWSER_THUMBSIZE") + ":</b>"));
    zoomLabel->set_use_markup (true);
    zoomBox->pack_start (*zoomLabel, Gtk::PACK_SHRINK, 4);
    zoomInButton  = Gtk::manage (new Gtk::Button ());
    zoomInButton->set_image (*Gtk::manage (new RTImage ("gtk-zoom-in.png")));
    zoomInButton->signal_pressed().connect (sigc::mem_fun(*batchQueue, &BatchQueue::zoomIn));
    zoomInButton->set_relief (Gtk::RELIEF_NONE);
    zoomInButton->set_tooltip_markup (M("FILEBROWSER_ZOOMINHINT"));
    zoomBox->pack_end (*zoomInButton, Gtk::PACK_SHRINK);
    zoomOutButton  = Gtk::manage (new Gtk::Button ());
    zoomOutButton->set_image (*Gtk::manage (new RTImage ("gtk-zoom-out.png")));
    zoomOutButton->signal_pressed().connect (sigc::mem_fun(*batchQueue, &BatchQueue::zoomOut));
    zoomOutButton->set_relief (Gtk::RELIEF_NONE);
    zoomOutButton->set_tooltip_markup (M("FILEBROWSER_ZOOMOUTHINT"));
    zoomBox->pack_end (*zoomOutButton, Gtk::PACK_SHRINK);
    bottomBox->pack_end (*zoomBox, Gtk::PACK_SHRINK);


    batchQueue->setBatchQueueListener (this);

    show_all ();

    if (batchQueue->loadBatchQueue ()) {
        g_idle_add_full (G_PRIORITY_LOW, processLoadedBatchQueueUIThread, batchQueue, NULL);
    }
}

void BatchQueuePanel::init (RTWindow *parent)
{
    this->parent = parent;

    saveFormatPanel->init (options.saveFormatBatch);
}

// it is expected to have a non null forceOrientation value on Preferences update only. In this case, qsize is ingored and computed automatically
void BatchQueuePanel::updateTab (int qsize, int forceOrientation)
{
    Gtk::Notebook *nb = (Gtk::Notebook *)(this->get_parent());

    if (forceOrientation > 0) {
        qsize = batchQueue->getEntries().size();
    }

    if ((forceOrientation == 0 && options.mainNBVertical) || (forceOrientation == 2)) {
        Gtk::VBox* vbb = Gtk::manage (new Gtk::VBox ());
        Gtk::Label* l;

        if(!qsize ) {
            vbb->pack_start (*Gtk::manage (new RTImage ("processing.png")));
            l = Gtk::manage (new Gtk::Label (Glib::ustring(" ") + M("MAIN_FRAME_BATCHQUEUE")) );
        } else if( start->get_active () ) {
            vbb->pack_start (*Gtk::manage (new RTImage ("processing-play.png")));
            l = Gtk::manage (new Gtk::Label (Glib::ustring(" ") + M("MAIN_FRAME_BATCHQUEUE") + " [" + Glib::ustring::format( qsize ) + "]"));
        } else {
            vbb->pack_start (*Gtk::manage (new RTImage ("processing-pause.png")));
            l = Gtk::manage (new Gtk::Label (Glib::ustring(" ") + M("MAIN_FRAME_BATCHQUEUE") + " [" + Glib::ustring::format( qsize ) + "]" ));
        }

        l->set_angle (90);
        vbb->pack_start (*l);
        vbb->set_spacing (2);
        vbb->set_tooltip_markup (M("MAIN_FRAME_BATCHQUEUE_TOOLTIP"));
        vbb->show_all ();

        if (nb) {
            nb->set_tab_label(*this, *vbb);
        }
    } else {
        Gtk::HBox* hbb = Gtk::manage (new Gtk::HBox ());

        if (!qsize ) {
            hbb->pack_start (*Gtk::manage (new RTImage ("processing.png")));
            hbb->pack_start (*Gtk::manage (new Gtk::Label (M("MAIN_FRAME_BATCHQUEUE") )));
        } else if ( start->get_active () ) {
            hbb->pack_start (*Gtk::manage (new RTImage ("processing-play.png")));
            hbb->pack_start (*Gtk::manage (new Gtk::Label (M("MAIN_FRAME_BATCHQUEUE") + " [" + Glib::ustring::format( qsize ) + "]" )));
        } else {
            hbb->pack_start (*Gtk::manage (new RTImage ("processing-pause.png")));
            hbb->pack_start (*Gtk::manage (new Gtk::Label (M("MAIN_FRAME_BATCHQUEUE") + " [" + Glib::ustring::format( qsize ) + "]" )));
        }

        hbb->set_spacing (2);
        hbb->set_tooltip_markup (M("MAIN_FRAME_BATCHQUEUE_TOOLTIP"));
        hbb->show_all ();

        if (nb) {
            nb->set_tab_label(*this, *hbb);
        }
    }
}

void BatchQueuePanel::queueSizeChanged (int qsize, bool queueEmptied, bool queueError, Glib::ustring queueErrorMessage)
{
    updateTab ( qsize);

    if (queueEmptied || queueError) {
        stopBatchProc ();
        fdir->set_sensitive (true);
        fformat->set_sensitive (true);

        SoundManager::playSoundAsync(options.sndBatchQueueDone);
    }

    if (queueError) {
        Gtk::MessageDialog msgd (queueErrorMessage, true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
        msgd.run ();
    }
}

void BatchQueuePanel::startBatchProc ()
{

    stopConnection.block (true);
    startConnection.block (true);
    stop->set_active (false);
    start->set_active (true);
    stopConnection.block (false);
    startConnection.block (false);

    if (batchQueue->hasJobs()) {
        fdir->set_sensitive (false);
        fformat->set_sensitive (false);
        saveOptions();
        batchQueue->startProcessing ();
    } else {
        stopBatchProc ();
    }

    updateTab (batchQueue->getEntries().size());
}

void BatchQueuePanel::stopBatchProc ()
{

    stopConnection.block (true);
    startConnection.block (true);
    stop->set_active (true);
    start->set_active (false);
    stopConnection.block (false);
    startConnection.block (false);
    updateTab (batchQueue->getEntries().size());
}

void BatchQueuePanel::addBatchQueueJobs ( std::vector<BatchQueueEntry*> &entries, bool head)
{

    batchQueue->addEntries (entries, head);

    if (stop->get_active () && autoStart->get_active ()) {
        startBatchProc ();
    }
}

bool BatchQueuePanel::canStartNext ()
{

    if (start->get_active ()) {
        return true;
    } else {
        fdir->set_sensitive (true);
        fformat->set_sensitive (true);
        return false;
    }
}

void BatchQueuePanel::saveOptions ()
{

    options.savePathTemplate    = outdirTemplate->get_text();
    options.saveUsePathTemplate = useTemplate->get_active();
    options.procQueueEnabled    = autoStart->get_active ();
}

void BatchQueuePanel::pathFolderButtonPressed ()
{

    Gtk::FileChooserDialog fc (getToplevelWindow (this), M("PREFERENCES_OUTDIRFOLDER"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER );
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

void BatchQueuePanel::formatChanged (Glib::ustring f)
{

    options.saveFormatBatch = saveFormatPanel->getFormat ();

}

bool BatchQueuePanel::handleShortcutKey (GdkEventKey* event)
{
    bool ctrl = event->state & GDK_CONTROL_MASK;

    if (ctrl) {
        switch(event->keyval) {
        case GDK_KEY_s:
            if (start->get_active()) {
                stopBatchProc();
            } else {
                startBatchProc();
            }

            return true;
        }
    }

    return batchQueue->keyPressed (event);
}
