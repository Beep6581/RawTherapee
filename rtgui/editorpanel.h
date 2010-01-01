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
#ifndef _EDITORPANEL_
#define _EDITORPANEL_

#include <gtkmm.h>
#include <imageareapanel.h>
#include <toolpanelcoord.h>
#include <profilepanel.h>
#include <rtengine.h>
#include <history.h>
#include <histogrampanel.h>
#include <thumbnail.h>
#include <saveasdlg.h>
#include <batchqueueentry.h>
#include <thumbnaillistener.h>
#include <navigator.h>

class EditorPanel;
struct EditorPanelIdleHelper {
    EditorPanel* epanel;
    bool destroyed;
    int pending;
};

class RTWindow;
class EditorPanel : public Gtk::VBox, 
                    public PParamsChangeListener,
                    public rtengine::ProgressListener,
                    public ThumbnailListener,
                    public HistoryBeforeLineListener {

    protected:      
        Gtk::Label *progressLabel;
        Gtk::ToggleButton* info;
        Gtk::ToggleButton* hidehp;
        Gtk::ToggleButton* beforeAfter;
        Gtk::HPaned* hpanedl;
        Gtk::HPaned* hpanedr;
        Gtk::HBox* statusBox;
        Gtk::Image* red;
        Gtk::Image* green;
        Gtk::VBox* leftbox, *vboxright;

        Gtk::Button* queueimg;
        Gtk::Button* saveimgas;
        Gtk::Button* sendtogimp;

        ImageAreaPanel* iarea;
        PreviewHandler* previewHandler;
        PreviewHandler* beforePreviewHandler;   // for the before-after view
        Navigator* navigator;
        ImageAreaPanel* beforeIarea;    // for the before-after view
        Gtk::VBox* beforeBox;
        Gtk::VBox* afterBox;
        Gtk::Label* beforeLabel;
        Gtk::Label* afterLabel;
        Gtk::HBox* beforeAfterBox;
        
        ProfilePanel* profilep;
        History* history;
        HistogramPanel* histogramPanel;
        ToolPanelCoordinator* tpc;
        RTWindow* parent;
        SaveAsDialog* saveAsDialog;
    
        Thumbnail* openThm;
        rtengine::InitialImage* isrc;
        rtengine::StagedImageProcessor* ipc;
        rtengine::StagedImageProcessor* beforeIpc;    // for the before-after view

        EditorPanelIdleHelper* epih;

        void open (Thumbnail* tmb, rtengine::InitialImage* isrc);
        void close ();

        rtengine::IImage16* processImage ();
        BatchQueueEntry*    createBatchQueueEntry ();
        int                 saveImage (rtengine::IImage16* img, Glib::ustring& fname, SaveFormat sf, bool findNewNameIfNeeded);

    public:
    
        static rtengine::InitialImage* loadImage (Thumbnail* tmb);

        EditorPanel (Thumbnail* tmb, rtengine::InitialImage* isrc);
        virtual ~EditorPanel ();

        bool beforeClosing ();
        void on_realize ();

        void setParent (RTWindow* p) { parent = p; }

        // progresslistener interface
        void setProgressState (int state);
        void error (Glib::ustring descr);
        void refreshProcessingState (bool state); // this is called by setProcessingState in the gtk thread
        void displayError (Glib::ustring descr);  // this is called by error in the gtk thread
        
        // PParamsChangeListener interface
        void procParamsChanged (rtengine::procparams::ProcParams* params, rtengine::ProcEvent ev, Glib::ustring descr, ParamsEdited* paramsEdited=NULL);

        // thumbnaillistener interface
        void procParamsChanged (Thumbnail* thm, int whoChangedIt);
        
        // HistoryBeforeLineListener
        void historyBeforeLineChanged (const rtengine::procparams::ProcParams& params);
        
        // event handlers
        void info_toggled ();
        void hideHistoryActivated ();
        void beforeAfterToggled ();
        void saveAsPressed ();
        void queueImgPressed ();
        void sendToGimpPressed ();

        void saveProfile ();
        Glib::ustring getShortName ();
        Glib::ustring getFileName ();
        bool handleShortcutKey (GdkEventKey* event);
        
        void saveOptions ();
};

#endif

