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
#include <batchqueuepanel.h>
#include <options.h>
#include <preferences.h>
#include <multilangmgr.h>
#include <rtwindow.h>

BatchQueuePanel::BatchQueuePanel () {

    batchQueue = new BatchQueue();

    // construct batch queue panel with the extra "start" and "stop" button
    Gtk::VBox* batchQueueButtonBox = Gtk::manage (new Gtk::VBox);
    start = Gtk::manage (new Gtk::ToggleButton (M("FILEBROWSER_STARTPROCESSING")));
    stop = Gtk::manage (new Gtk::ToggleButton (M("FILEBROWSER_STOPPROCESSING")));
    autoStart = Gtk::manage (new Gtk::CheckButton ("Auto start"));
    start->set_tooltip_text (M("FILEBROWSER_STARTPROCESSINGHINT"));
    stop->set_tooltip_text (M("FILEBROWSER_STOPPROCESSINGHINT"));
    autoStart->set_tooltip_text ("Start processing automatically when a new job arrives");
    start->set_active (false);
    stop->set_active (true);
    autoStart->set_active (options.procQueueEnabled);
    
    start->set_image (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-media-play"), Gtk::ICON_SIZE_BUTTON)));
    startConnection = start->signal_toggled().connect (sigc::mem_fun(*this, &BatchQueuePanel::startBatchProc));    
    stop->set_image (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-media-stop"), Gtk::ICON_SIZE_BUTTON)));
    stopConnection = stop->signal_toggled().connect (sigc::mem_fun(*this, &BatchQueuePanel::stopBatchProc));    
    batchQueueButtonBox->pack_start (*start, Gtk::PACK_SHRINK, 4);
    batchQueueButtonBox->pack_start (*stop, Gtk::PACK_SHRINK, 4);
    batchQueueButtonBox->pack_start (*autoStart, Gtk::PACK_SHRINK, 4);

    // Output directory selection
    fdir = Gtk::manage (new Gtk::Frame (M("PREFERENCES_OUTDIR")));
    Gtk::VBox* odvb = Gtk::manage (new Gtk::VBox ());
    odvb->set_border_width (4);
    Gtk::HBox* hb2 = Gtk::manage (new Gtk::HBox ());
    useTemplate = Gtk::manage (new Gtk::RadioButton (M("PREFERENCES_OUTDIRTEMPLATE")+":"));
    hb2->pack_start (*useTemplate, Gtk::PACK_SHRINK,4);
    outdirTemplate = Gtk::manage (new Gtk::Entry ());
    hb2->pack_start (*outdirTemplate);
    odvb->pack_start (*hb2, Gtk::PACK_SHRINK, 4);
    outdirTemplate->set_tooltip_markup (M("PREFERENCES_OUTDIRTEMPLATEHINT"));
    useTemplate->set_tooltip_markup (M("PREFERENCES_OUTDIRTEMPLATEHINT"));
    Gtk::HBox* hb3 = Gtk::manage (new Gtk::HBox ());
    useFolder = Gtk::manage (new Gtk::RadioButton (M("PREFERENCES_OUTDIRFOLDER")+":"));
    hb3->pack_start (*useFolder, Gtk::PACK_SHRINK,4);
    outdirFolder = Gtk::manage (new Gtk::FileChooserButton (M("PREFERENCES_OUTDIRFOLDER"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    hb3->pack_start (*outdirFolder);
    odvb->pack_start (*hb3, Gtk::PACK_SHRINK, 4);
    outdirFolder->set_tooltip_markup (M("PREFERENCES_OUTDIRFOLDERHINT"));
    useFolder->set_tooltip_markup (M("PREFERENCES_OUTDIRFOLDERHINT"));
    Gtk::RadioButton::Group g = useTemplate->get_group();
    useFolder->set_group (g);
    fdir->add (*odvb);

    // Output file format selection
    fformat = Gtk::manage (new Gtk::Frame (M("PREFERENCES_FILEFORMAT")));
    saveFormatPanel = Gtk::manage (new SaveFormatPanel ());
    fformat->add (*saveFormatPanel);       

    saveFormatPanel->init (options.saveFormat);
    outdirTemplate->set_text (options.savePathTemplate);
    if (Glib::file_test (options.savePathFolder, Glib::FILE_TEST_IS_DIR)) 
        outdirFolder->set_filename (options.savePathFolder);
    useTemplate->set_active (options.saveUsePathTemplate);
    useFolder->set_active (!options.saveUsePathTemplate);

    // setup signal handlers
    outdirTemplate->signal_changed().connect (sigc::mem_fun(*this, &BatchQueuePanel::saveOptions));    
    outdirFolder->signal_current_folder_changed().connect (sigc::mem_fun(*this, &BatchQueuePanel::saveOptions));    
    useTemplate->signal_toggled().connect (sigc::mem_fun(*this, &BatchQueuePanel::saveOptions));    
    useFolder->signal_toggled().connect (sigc::mem_fun(*this, &BatchQueuePanel::saveOptions));    
    saveFormatPanel->setListener (this);

    // setup button bar
    topBox = Gtk::manage (new Gtk::HBox ());
    pack_start (*topBox, Gtk::PACK_SHRINK);

    topBox->pack_start (*batchQueueButtonBox, Gtk::PACK_SHRINK, 4);
    topBox->pack_start (*fdir);
    topBox->pack_start (*fformat, Gtk::PACK_SHRINK, 4);

    // add middle browser area 
    Gtk::HBox* hBox = Gtk::manage (new Gtk::HBox ());
    pack_start (*batchQueue);

    // lower box with thumbnail zoom
    bottomBox = Gtk::manage (new Gtk::HBox ());
    pack_start (*bottomBox, Gtk::PACK_SHRINK);

    // change thumbnail arrangement button
    hAlignIcon = new Gtk::Image (argv0+"/images/horizontals.png");
    vAlignIcon = new Gtk::Image (argv0+"/images/verticals.png");
    hAlignIcon->show ();
    vAlignIcon->show ();
    chAlign = Gtk::manage (new Gtk::Button ());
    chAlign->show ();
    bottomBox->pack_end (*chAlign, Gtk::PACK_SHRINK);
    chAlign->set_image (*hAlignIcon);
    chAlign->set_relief (Gtk::RELIEF_NONE);
    chAlign->signal_pressed().connect (sigc::mem_fun(*this, &BatchQueuePanel::arrangementButtonPressed));    
    chAlign->set_tooltip_text (M("FILEBROWSER_ARRANGEMENTHINT"));
    bottomBox->pack_end (*Gtk::manage (new Gtk::VSeparator), Gtk::PACK_SHRINK, 4);
    if (options.fbArrangement==1) 
        chAlign->set_image (*vAlignIcon);
    else 
        chAlign->set_image (*hAlignIcon);
    arrangementButtonPressed ();

    // thumbnail zoom
    Gtk::HBox* zoomBox = Gtk::manage (new Gtk::HBox ());
    zoomBox->pack_start (*Gtk::manage (new Gtk::VSeparator), Gtk::PACK_SHRINK, 4);
    Gtk::Label* zoomLabel = Gtk::manage (new Gtk::Label (Glib::ustring("<b>")+M("FILEBROWSER_THUMBSIZE")+":</b>"));
    zoomLabel->set_use_markup (true);
    zoomBox->pack_start (*zoomLabel, Gtk::PACK_SHRINK, 4);   
    zoomInButton  = Gtk::manage (new Gtk::Button ());
    zoomInButton->set_image (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-zoom-in"), Gtk::ICON_SIZE_SMALL_TOOLBAR)));
    zoomInButton->signal_pressed().connect (sigc::mem_fun(*batchQueue, &BatchQueue::zoomIn));    
    zoomInButton->set_relief (Gtk::RELIEF_NONE);
    zoomInButton->set_tooltip_text (M("FILEBROWSER_ZOOMINHINT"));
    zoomBox->pack_end (*zoomInButton, Gtk::PACK_SHRINK);
    zoomOutButton  = Gtk::manage (new Gtk::Button ());
    zoomOutButton->set_image (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-zoom-out"), Gtk::ICON_SIZE_SMALL_TOOLBAR)));
    zoomOutButton->signal_pressed().connect (sigc::mem_fun(*batchQueue, &BatchQueue::zoomOut));    
    zoomOutButton->set_relief (Gtk::RELIEF_NONE);
    zoomOutButton->set_tooltip_text (M("FILEBROWSER_ZOOMOUTHINT"));
    zoomBox->pack_end (*zoomOutButton, Gtk::PACK_SHRINK);   
    bottomBox->pack_end (*zoomBox, Gtk::PACK_SHRINK);   


    batchQueue->setBatchQueueListener (this);

    show_all ();
}

void BatchQueuePanel::arrangementButtonPressed () {

    if (chAlign->get_image()==hAlignIcon) {
        chAlign->set_image (*vAlignIcon);
        batchQueue->setArrangement (BatchQueue::TB_Vertical);
    }
    else {
        chAlign->set_image (*hAlignIcon);
        batchQueue->setArrangement (BatchQueue::TB_Horizontal);
    }
}

void BatchQueuePanel::queueSizeChanged (int qsize) {

    // TODO: write it somewhere
}

void BatchQueuePanel::imageProcessingReady (Glib::ustring fname) {

    parent->imageDeveloped (fname);
}

void BatchQueuePanel::startBatchProc () {

    stopConnection.block (true);
    startConnection.block (true);
    stop->set_active (false);
    start->set_active (true);
    stopConnection.block (false);
    startConnection.block (false);
    
    if (batchQueue->hasJobs()) {
        fdir->set_sensitive (false);
        fformat->set_sensitive (false);
        batchQueue->startProcessing ();
    }
    else 
        stopBatchProc ();
}

void BatchQueuePanel::stopBatchProc () {

    stopConnection.block (true);
    startConnection.block (true);
    stop->set_active (true);
    start->set_active (false);
    stopConnection.block (false);
    startConnection.block (false);
}

void BatchQueuePanel::addBatchQueueJob (BatchQueueEntry* bqe, bool head) {

    batchQueue->addEntry (bqe, head);
    
    if (stop->get_active () && autoStart->get_active ())
        startBatchProc ();
}

void BatchQueuePanel::queueEmpty () {

    stopBatchProc ();
    fdir->set_sensitive (true);
    fformat->set_sensitive (true);
}

bool BatchQueuePanel::canStartNext () {

    if (start->get_active ())
        return true;
    else {
        fdir->set_sensitive (true);
        fformat->set_sensitive (true);
        return false;
    }
}

void BatchQueuePanel::saveOptions () { 

    options.saveFormat          = saveFormatPanel->getFormat ();
    options.savePathTemplate    = outdirTemplate->get_text();
    options.savePathFolder      = outdirFolder->get_filename();
    options.saveUsePathTemplate = useTemplate->get_active();
    options.procQueueEnabled    = autoStart->get_active ();
}

void BatchQueuePanel::formatChanged () {
    
    saveOptions ();
}
