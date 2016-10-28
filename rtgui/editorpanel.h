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
#ifndef _EDITORPANEL_
#define _EDITORPANEL_

#include <gtkmm.h>
#include "imageareapanel.h"
#include "toolpanelcoord.h"
#include "profilepanel.h"
#include "../rtengine/rtengine.h"
#include "history.h"
#include "histogrampanel.h"
#include "thumbnail.h"
#include "saveasdlg.h"
#include "batchqueueentry.h"
#include "thumbnaillistener.h"
#include "navigator.h"
#include "progressconnector.h"
#include "filepanel.h"

class EditorPanel;
class MyProgressBar;

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
    public HistoryBeforeLineListener,
    public rtengine::HistogramListener
{
private:

    Glib::ustring lastSaveAsFileName;
    bool realized;

protected:
    MyProgressBar  *progressLabel;
    Gtk::ToggleButton* info;
    Gtk::ToggleButton* hidehp;
    Gtk::ToggleButton* tbShowHideSidePanels;
    Gtk::ToggleButton* tbTopPanel_1;
    Gtk::ToggleButton* tbRightPanel_1;
    Gtk::ToggleButton* tbBeforeLock;
    //bool bAllSidePanelsVisible;
    Gtk::ToggleButton* beforeAfter;
    Gtk::HPaned* hpanedl;
    Gtk::HPaned* hpanedr;
    Gtk::Image *iHistoryShow, *iHistoryHide;
    Gtk::Image *iTopPanel_1_Show, *iTopPanel_1_Hide;
    Gtk::Image *iRightPanel_1_Show, *iRightPanel_1_Hide;
    Gtk::Image *iShowHideSidePanels;
    Gtk::Image *iShowHideSidePanels_exit;
    Gtk::Image *iBeforeLockON, *iBeforeLockOFF;
    Gtk::VBox *leftbox;
    Gtk::VBox *vboxright;

    Gtk::Button* queueimg;
    Gtk::Button* saveimgas;
    Gtk::Button* sendtogimp;
    Gtk::Button* navSync;
    Gtk::Button* navNext;
    Gtk::Button* navPrev;

    class ColorManagementToolbar;
    std::unique_ptr<ColorManagementToolbar> colorMgmtToolBar;

    ImageAreaPanel* iareapanel;
    PreviewHandler* previewHandler;
    PreviewHandler* beforePreviewHandler;   // for the before-after view
    Navigator* navigator;
    ImageAreaPanel* beforeIarea;    // for the before-after view
    Gtk::VBox* beforeBox;
    Gtk::VBox* afterBox;
    Gtk::Label* beforeLabel;
    Gtk::Label* afterLabel;
    Gtk::HBox* beforeAfterBox;
    Gtk::HBox* beforeHeaderBox;
    Gtk::HBox* afterHeaderBox;

    Gtk::Frame* ppframe;
    ProfilePanel* profilep;
    History* history;
    HistogramPanel* histogramPanel;
    ToolPanelCoordinator* tpc;
    RTWindow* parent;
    //SaveAsDialog* saveAsDialog;
    BatchToolPanelCoordinator* btpCoordinator;
    FilePanel* fPanel;

    bool firstProcessingDone;

    Thumbnail* openThm;  // may get invalid on external delete event
    Glib::ustring fname;  // must be saved separately

    rtengine::InitialImage* isrc;
    rtengine::StagedImageProcessor* ipc;
    rtengine::StagedImageProcessor* beforeIpc;    // for the before-after view

    EditorPanelIdleHelper* epih;

    void close ();

    BatchQueueEntry*    createBatchQueueEntry ();
    bool                idle_imageSaved (ProgressConnector<int> *pc, rtengine::IImage16* img, Glib::ustring fname, SaveFormat sf);
    bool                idle_saveImage (ProgressConnector<rtengine::IImage16*> *pc, Glib::ustring fname, SaveFormat sf);
    bool                idle_sendToGimp ( ProgressConnector<rtengine::IImage16*> *pc, Glib::ustring fname);
    bool                idle_sentToGimp (ProgressConnector<int> *pc, rtengine::IImage16* img, Glib::ustring filename);
    int err;

    time_t processingStartedTime;

    sigc::connection ShowHideSidePanelsconn;

    bool isProcessing;


public:
    explicit EditorPanel (FilePanel* filePanel = nullptr);
    virtual ~EditorPanel ();

    void open (Thumbnail* tmb, rtengine::InitialImage* isrc);
    void setAspect ();
    void on_realize ();
    void leftPaneButtonReleased (GdkEventButton *event);
    void rightPaneButtonReleased (GdkEventButton *event);

    void setParent (RTWindow* p)
    {
        parent = p;
    }
    void writeOptions();

    void showTopPanel (bool show);
    bool isRealized()
    {
        return realized;
    }
    // progresslistener interface
    void setProgress (double p);
    void setProgressStr (Glib::ustring str);
    void setProgressState (bool inProcessing);
    void error (Glib::ustring title, Glib::ustring descr);
    void displayError (Glib::ustring title, Glib::ustring descr);  // this is called by error in the gtk thread
    void refreshProcessingState (bool inProcessing); // this is called by setProcessingState in the gtk thread

    // PParamsChangeListener interface
    void procParamsChanged (rtengine::procparams::ProcParams* params, rtengine::ProcEvent ev, Glib::ustring descr, ParamsEdited* paramsEdited = nullptr);

    // thumbnaillistener interface
    void procParamsChanged (Thumbnail* thm, int whoChangedIt);

    // HistoryBeforeLineListener
    void historyBeforeLineChanged (const rtengine::procparams::ProcParams& params);

    // HistogramListener
    void histogramChanged (LUTu & histRed, LUTu & histGreen, LUTu & histBlue, LUTu & histLuma, LUTu & histToneCurve, LUTu & histLCurve, LUTu & histCCurve,/* LUTu & histCLurve, LUTu & histLLCurve,*/ LUTu & histLCAM, LUTu & histCCAM,
                           LUTu & histRedRaw, LUTu & histGreenRaw, LUTu & histBlueRaw, LUTu & histChroma, LUTu & histLRETI);

    // event handlers
    void info_toggled ();
    void hideHistoryActivated ();
    void tbRightPanel_1_toggled ();
    void tbTopPanel_1_toggled ();
    void beforeAfterToggled ();
    void tbBeforeLock_toggled();
    void saveAsPressed ();
    void queueImgPressed ();
    void sendToGimpPressed ();
    void openNextEditorImage ();
    void openPreviousEditorImage ();
    void syncFileBrowser ();

    void tbTopPanel_1_visible (bool visible);
    bool CheckSidePanelsVisibility();
    void tbShowHideSidePanels_managestate();
    void toggleSidePanels();
    void toggleSidePanelsZoomFit();

    void saveProfile ();
    Glib::ustring getShortName ();
    Glib::ustring getFileName ();
    bool handleShortcutKey (GdkEventKey* event);

    bool getIsProcessing() const
    {
        return isProcessing;
    }
    void updateTPVScrollbar (bool hide);
    void updateTabsUsesIcons (bool useIcons);
    void updateHistogramPosition (int oldPosition, int newPosition);

    Gtk::Paned *catalogPane;
};

#endif

