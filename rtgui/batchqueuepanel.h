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
#include <batchqueue.h>
#include <saveformatpanel.h>

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
        Gtk::FileChooserButton* outdirFolder;
        Gtk::RadioButton* useTemplate;
        Gtk::RadioButton* useFolder;
        SaveFormatPanel* saveFormatPanel;
        Gtk::Frame *fdir, *fformat;

        Gtk::Image* hAlignIcon;
        Gtk::Image* vAlignIcon;
        Gtk::Button* chAlign;

        RTWindow* parent;
        BatchQueue* batchQueue;
        Gtk::HBox* bottomBox;
        Gtk::HBox* topBox;

    public:
        
        BatchQueuePanel ();

        void setParent (RTWindow* p) { parent = p; }
        void arrangementButtonPressed ();

        void addBatchQueueJob (BatchQueueEntry* bqe, bool head=false);

        // batchqueuelistener interface
        void queueSizeChanged     (int qsize);
        void imageProcessingReady (Glib::ustring fname);
        void queueEmpty           ();
        bool canStartNext         ();
        
        void startBatchProc ();
        void stopBatchProc ();
        
        void saveOptions ();
        void formatChanged ();       
};

#endif

