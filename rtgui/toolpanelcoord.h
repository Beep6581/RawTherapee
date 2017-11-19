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

#include "../rtengine/rtengine.h"
#include "toolpanel.h"
#include <vector>
#include "pparamschangelistener.h"
#include "profilechangelistener.h"
#include "imageareatoollistener.h"
#include <gtkmm.h>
#include "whitebalance.h"
#include "coarsepanel.h"
#include "tonecurve.h"
#include "vibrance.h"
#include "colorappearance.h"
#include "shadowshighlights.h"
#include "impulsedenoise.h"
#include "defringe.h"
#include "dirpyrdenoise.h"
#include "epd.h"
#include "sharpening.h"
#include "labcurve.h"
#include "exifpanel.h"
#include "iptcpanel.h"
#include "crop.h"
#include "icmpanel.h"
#include "resize.h"
#include "chmixer.h"
#include "blackwhite.h"
#include "cacorrection.h"
#include "lensprofile.h"
#include "distortion.h"
#include "perspective.h"
#include "rotate.h"
#include "vignetting.h"
#include "retinex.h"
#include "gradient.h"
#include "locallab.h"
#include "pcvignette.h"
#include "toolbar.h"
#include "lensgeom.h"
#include "lensgeomlistener.h"
#include "wavelet.h"
#include "dirpyrequalizer.h"
#include "hsvequalizer.h"
#include "preprocess.h"
#include "bayerpreprocess.h"
#include "bayerprocess.h"
#include "xtransprocess.h"
#include "darkframe.h"
#include "flatfield.h"
#include "sensorbayer.h"
#include "sensorxtrans.h"
#include "rawcacorrection.h"
#include "rawexposure.h"
#include "bayerrawexposure.h"
#include "xtransrawexposure.h"
#include "sharpenmicro.h"
#include "sharpenedge.h"
#include "rgbcurves.h"
#include "colortoning.h"
#include "filmsimulation.h"
#include "prsharpening.h"
#include "fattaltonemap.h"
#include "guiutils.h"

class ImageEditorCoordinator;

class ToolPanelCoordinator :    public ToolPanelListener,
    public ToolBarListener,
    public ProfileChangeListener,
    public WBProvider,
    public DFProvider,
    public FFProvider,
    public LensGeomListener,
    public SpotWBListener,
    public CropPanelListener,
    public ICMPanelListener,
    public ImageAreaToolListener,
    public rtengine::ImageTypeListener
{

protected:

    WhiteBalance* whitebalance;
    Vignetting* vignetting;
    Gradient* gradient;
    Locallab* locallab;
    Retinex*  retinex;
    PCVignette* pcvignette;
    LensGeometry* lensgeom;
    LensProfilePanel* lensProf;
    Rotate* rotate;
    Distortion* distortion;
    PerspCorrection* perspective;
    CACorrection* cacorrection;
    ColorAppearance* colorappearance;
    Vibrance* vibrance;
    ChMixer* chmixer;
    BlackWhite* blackwhite;
    Resize* resize;
    PrSharpening* prsharpening;
    ICMPanel* icm;
    Crop* crop;
    ToneCurve* toneCurve;
    ShadowsHighlights* shadowshighlights;
    Defringe* defringe;
    ImpulseDenoise* impulsedenoise;
    DirPyrDenoise* dirpyrdenoise;
    EdgePreservingDecompositionUI *epd;
    Sharpening* sharpening;
    SharpenEdge* sharpenEdge;
    SharpenMicro* sharpenMicro;
    LCurve* lcurve;
    RGBCurves* rgbcurves;
    ColorToning* colortoning;
    Wavelet * wavelet;
    DirPyrEqualizer* dirpyrequalizer;
    HSVEqualizer* hsvequalizer;
    FilmSimulation *filmSimulation;
    SensorBayer * sensorbayer;
    SensorXTrans * sensorxtrans;
    BayerProcess* bayerprocess;
    XTransProcess* xtransprocess;
    BayerPreProcess* bayerpreprocess;
    PreProcess* preprocess;
    DarkFrame* darkframe;
    FlatField* flatfield;
    RAWCACorr* rawcacorrection;
    RAWExposure* rawexposure;
    BayerRAWExposure* bayerrawexposure;
    XTransRAWExposure* xtransrawexposure;
    FattalToneMapping *fattal;

    std::vector<PParamsChangeListener*> paramcListeners;

    rtengine::StagedImageProcessor* ipc;

    std::vector<ToolPanel*> toolPanels;
    ToolVBox* exposurePanel;
    ToolVBox* detailsPanel;
    ToolVBox* colorPanel;
    ToolVBox* transformPanel;
    ToolVBox* rawPanel;
    ToolVBox* waveletPanel;
    ToolVBox* locallabPanel;
    Gtk::Notebook* metadataPanel;
    ExifPanel* exifpanel;
    IPTCPanel* iptcpanel;
    ToolBar* toolBar;

    TextOrIcon* toiE;
    TextOrIcon* toiD;
    TextOrIcon* toiC;
    TextOrIcon* toiT;
    TextOrIcon* toiR;
    TextOrIcon* toiM;
    TextOrIcon* toiW;
    TextOrIcon* toiL;

    Gtk::Image* imgPanelEnd[7];
    Gtk::VBox* vbPanelEnd[7];

    Gtk::ScrolledWindow* exposurePanelSW;
    Gtk::ScrolledWindow* detailsPanelSW;
    Gtk::ScrolledWindow* colorPanelSW;
    Gtk::ScrolledWindow* transformPanelSW;
    Gtk::ScrolledWindow* rawPanelSW;
    Gtk::ScrolledWindow* waveletPanelSW;
    Gtk::ScrolledWindow* locallabPanelSW;

    std::vector<MyExpander*> expList;

    bool hasChanged;

    void addPanel (Gtk::Box* where, FoldableToolPanel* panel, int level = 1);
    void foldThemAll (GdkEventButton* event);
    void updateVScrollbars (bool hide);
    void updateTabsHeader (bool useIcons);

private:

    EditDataProvider *editDataProvider;

public:

    CoarsePanel* coarse;
    Gtk::Notebook* toolPanelNotebook;

    ToolPanelCoordinator (bool batch = false);
    virtual ~ToolPanelCoordinator ();

    bool getChangedState                ()
    {
        return hasChanged;
    }
    void updateCurveBackgroundHistogram (LUTu & histToneCurve, LUTu & histLCurve, LUTu & histCCurve, /*LUTu & histCLurve, LUTu & histLLCurve,*/ LUTu & histLCAM,  LUTu & histCCAM, LUTu & histRed, LUTu & histGreen, LUTu & histBlue, LUTu & histLuma, LUTu & histLRETI);
    void foldAllButOne (Gtk::Box* parent, FoldableToolPanel* openedSection);

    // multiple listeners can be added that are notified on changes (typical: profile panel and the history)
    void addPParamsChangeListener   (PParamsChangeListener* pp)
    {
        paramcListeners.push_back (pp);
    }

    // toolpanellistener interface
    void panelChanged   (rtengine::ProcEvent event, const Glib::ustring& descr);

    void imageTypeChanged (bool isRaw, bool isBayer, bool isXtrans);
    // profilechangelistener interface
    void profileChange  (const rtengine::procparams::PartialProfile* nparams, rtengine::ProcEvent event, const Glib::ustring& descr, const ParamsEdited* paramsEdited = nullptr);
    void setDefaults    (rtengine::procparams::ProcParams* defparams);

    // DirSelectionListener interface
    void dirSelected (const Glib::ustring& dirname, const Glib::ustring& openfile);

    // to support the GUI:
    CropGUIListener* getCropGUIListener (); // through the CropGUIListener the editor area can notify the "crop" ToolPanel when the crop selection changes

    // init the toolpanelcoordinator with an image & close it
    void initImage          (rtengine::StagedImageProcessor* ipc_, bool israw);
    void closeImage         ();

    // update the "expanded" state of the Tools
    void updateToolState    ();
    void openAllTools       ();
    void closeAllTools      ();
    // read/write the "expanded" state of the expanders & read/write the crop panel settings (ratio, guide type, etc.)
    void readOptions        ();
    void writeOptions       ();
    void writeToolExpandedStatus (std::vector<int> &tpOpen);


    // wbprovider interface
    void getAutoWB (double& temp, double& green, double equal, double tempBias)
    {
        if (ipc) {
            ipc->getAutoWB (temp, green, equal, tempBias);
        }
    }
    void getCamWB (double& temp, double& green)
    {
        if (ipc) {
            ipc->getCamWB (temp, green);
        }
    }

    //DFProvider interface
    rtengine::RawImage* getDF();

    //FFProvider interface
    rtengine::RawImage* getFF();
    Glib::ustring GetCurrentImageFilePath();

    // rotatelistener interface
    void straightenRequested ();
    void autoCropRequested ();
    double autoDistorRequested ();

    // spotwblistener interface
    void spotWBRequested (int size);

    // croppanellistener interface
    void cropSelectRequested ();

    // icmpanellistener interface
    void saveInputICCReference (Glib::ustring fname, bool apply_wb);

    // imageareatoollistener interface
    void spotWBselected (int x, int y, Thumbnail* thm = nullptr);
    void cropSelectionReady ();
    void rotateSelectionReady (double rotate_deg, Thumbnail* thm = nullptr);
    ToolBar* getToolBar ()
    {
        return toolBar;
    }
    int  getSpotWBRectSize ();
    CropGUIListener* startCropEditing (Thumbnail* thm = nullptr)
    {
        return crop;
    }

    void updateTPVScrollbar (bool hide);
    void updateTabsUsesIcons (bool useIcons);
    bool handleShortcutKey (GdkEventKey* event);

    // ToolBarListener interface
    void toolSelected (ToolMode tool);
    void editModeSwitchedOff ();

    void setEditProvider (EditDataProvider *provider);
};

#endif
