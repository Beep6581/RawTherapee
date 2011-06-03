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
#ifndef __TOOLPANELCCORD__
#define __TOOLPANELCCORD__

#include <rtengine.h>
#include <toolpanel.h>
#include <vector>
#include <pparamschangelistener.h>
#include <profilechangelistener.h>
#include <imageareatoollistener.h>
#include <gtkmm.h>
#include <whitebalance.h>
#include <colorboost.h>
#include <coarsepanel.h>
#include <tonecurve.h>
#include <shadowshighlights.h>
#include <lumadenoise.h>
#include <colordenoise.h>
#include <impulsedenoise.h>
#include <defringe.h>
#include <dirpyrdenoise.h>
#include <sharpening.h>
#include <labcurve.h>
#include <exifpanel.h>
#include <iptcpanel.h>
#include <crop.h>
#include <icmpanel.h>
#include <resize.h>
#include <chmixer.h>
#include <hlrec.h>
#include <colorshift.h>
#include <cacorrection.h>
#include <distortion.h>
#include <perspective.h>
#include <rotate.h>
#include <vignetting.h>
#include <toolbar.h>
#include <lensgeom.h>
#include <lensgeomlistener.h>
#include <equalizer.h>
#include <dirpyrequalizer.h>
#include <hsvequalizer.h>
#include <rawprocess.h>
#include <preprocess.h>
#include <darkframe.h>
#include <flatfield.h>
#include <rawcacorrection.h>
#include <rawexposure.h>

class ImageEditorCoordinator;

class ToolPanelCoordinator :    public ToolPanelListener, 
                                public ProfileChangeListener, 
                                public WBProvider,
                                public DFProvider,
                                public FFProvider,
                                public LensGeomListener,
                                public SpotWBListener,
                                public CropPanelListener,
                                public ICMPanelListener,
                                public ImageAreaToolListener {

    protected:
        
        WhiteBalance* whitebalance;
        Vignetting* vignetting;
        LensGeometry* lensgeom;
        Rotate* rotate;
        Distortion* distortion;
        PerspCorrection* perspective;
        CACorrection* cacorrection;
        //ColorShift* colorshift;
        HLRecovery* hlrecovery;
        ChMixer* chmixer;
        //ColorBoost* colorboost;
        Resize* resize;
        ICMPanel* icm;
        Crop* crop;
        ToneCurve* toneCurve;
        ShadowsHighlights* shadowshighlights;
        //LumaDenoise* lumadenoise;
        //ColorDenoise* colordenoise;
        Defringe* defringe;
        ImpulseDenoise* impulsedenoise;
        DirPyrDenoise* dirpyrdenoise;
        Sharpening* sharpening;
        LCurve* lcurve;
        //Equalizer * equalizer;
        DirPyrEqualizer * dirpyrequalizer;
        HSVEqualizer * hsvequalizer;
        RawProcess* rawprocess;
        PreProcess* preprocess;
        DarkFrame* darkframe;
        FlatField* flatfield;
        RAWCACorr* rawcacorrection;
        RAWExposure* rawexposure;

        std::vector<PParamsChangeListener*> paramcListeners;

        rtengine::StagedImageProcessor* ipc;

        std::vector<ToolPanel*> toolPanels;
        Gtk::VBox* exposurePanel;
        Gtk::VBox* detailsPanel;
        Gtk::VBox* colorPanel;
        Gtk::VBox* transformPanel;
        Gtk::VBox* rawPanel;
        Gtk::Notebook* metadataPanel;
        ExifPanel* exifpanel;
        IPTCPanel* iptcpanel;
        ToolBar* toolBar;

        Gtk::Image* imgPanelEnd[5];
        Gtk::VBox* vbPanelEnd[5];

        Gtk::ScrolledWindow* exposurePanelSW;
        Gtk::ScrolledWindow* detailsPanelSW;
        Gtk::ScrolledWindow* colorPanelSW;
        Gtk::ScrolledWindow* transformPanelSW;
        Gtk::ScrolledWindow* rawPanelSW;

        std::vector<Gtk::Expander*> expList;
        
        bool hasChanged;

        void addPanel (Gtk::Box* where, FoldableToolPanel* panel, Glib::ustring label);
        void foldThemAll (GdkEventButton* event);

    public:
    
        CoarsePanel* coarse;
        Gtk::Notebook* toolPanelNotebook;

        ToolPanelCoordinator ();
        ~ToolPanelCoordinator ();

        bool getChangedState                ()                                      { return hasChanged; }
		void updateCurveBackgroundHistogram (LUTu & histToneCurve, LUTu & histLCurves);
        void foldAllButOne (Gtk::Box* parent, FoldableToolPanel* openedSection);

        // multiple listeners can be added that are notified on changes (typical: profile panel and the history)
        void addPParamsChangeListener   (PParamsChangeListener* pp) { paramcListeners.push_back (pp); }
        
        // toolpanellistener interface
        void panelChanged   (rtengine::ProcEvent event, const Glib::ustring& descr);

        // profilechangelistener interface
        void profileChange  (const rtengine::procparams::ProcParams* nparams, rtengine::ProcEvent event, const Glib::ustring& descr, const ParamsEdited* paramsEdited=NULL);    
        void setDefaults    (rtengine::procparams::ProcParams* defparams);

        // to support the GUI:
        CropGUIListener* getCropGUIListener (); // through the CropGUIListener the editor area can notify the "crop" ToolPanel when the crop selection changes
        
        // init the toolpanelcoordinator with an image & close it        
        void initImage          (rtengine::StagedImageProcessor* ipc_, bool israw);
        void closeImage         ();

        // read/write the "expanded" state of the expanders & read/write the crop panel settings (ratio, guide type, etc.)
        void readOptions        ();
        void writeOptions       ();       

        // wbprovider interface
        void getAutoWB (double& temp, double& green) { if (ipc) ipc->getAutoWB (temp, green); }
        void getCamWB (double& temp, double& green)  { if (ipc) ipc->getCamWB (temp, green); }

        //DFProvider interface
        rtengine::RawImage* getDF();

        //FFProvider interface
        rtengine::RawImage* getFF();

        // rotatelistener interface
        void straightenRequested ();
        void autoCropRequested ();
        double autoDistorRequested ();

        // spotwblistener interface
        void spotWBRequested (int size);

        // croppanellistener interface
        void cropSelectRequested ();
        
        // icmpanellistener interface
        void saveInputICCReference (Glib::ustring fname);
    
        // imageareatoollistener interface
        void spotWBselected (int x, int y, Thumbnail* thm=NULL);
        void cropSelectionReady ();
        void rotateSelectionReady (double rotate_deg, Thumbnail* thm=NULL);
        ToolBar* getToolBar () { return toolBar; }
        int  getSpotWBRectSize ();
        CropGUIListener* startCropEditing (Thumbnail* thm=NULL) { return crop; }

        bool handleShortcutKey (GdkEventKey* event);
};

#endif
