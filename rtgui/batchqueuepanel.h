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
#ifndef _BATCHQUEUEPANEL_
#define _BATCHQUEUEPANEL_

#include <gtkmm.h>
#include "batchqueue.h"
#include "saveformatpanel.h"
#include "guiutils.h"

class RTWindow;
class BatchQueuePanel : public Gtk::VBox,
                        public BatchQueueListener,
                        public FormatChangeListener {

        Gtk::Button* zoomInButton;
        Gtk::Button* zoomOutButton;
        Gtk::ToggleButton* start;
        Gtk::ToggleButton* stop;
        Gtk::CheckButton* autoStart;
        sigc::connection startConnection;
        sigc::connection stopConnection;

        Gtk::Entry* outdirTemplate;
        MyFileChooserButton* outdirFolder;
        Gtk::Button* outdirFolderButton;
        Gtk::RadioButton* useTemplate;
        Gtk::RadioButton* useFolder;
        SaveFormatPanel* saveFormatPanel;
        Gtk::Frame *fdir, *fformat;

        RTWindow* parent;
        BatchQueue* batchQueue;
        Gtk::HBox* bottomBox;
        Gtk::HBox* topBox;

    public:
        
        BatchQueuePanel ();

        void setParent (RTWindow* p) { parent = p; }

        void addBatchQueueJobs (std::vector<BatchQueueEntry*> &entries , bool head=false);

        // batchqueuelistener interface
        void queueSizeChanged     (int qsize, bool queueEmptied, bool queueError, Glib::ustring queueErrorMessage);
        bool canStartNext         ();
        
        void startBatchProc ();
        void stopBatchProc ();
        
        void saveOptions ();
        void pathFolderChanged ();
        void pathFolderButtonPressed ();
        void formatChanged (Glib::ustring f);
        void updateTab (int qsize);
        
        bool handleShortcutKey (GdkEventKey* event);
};
#endif

