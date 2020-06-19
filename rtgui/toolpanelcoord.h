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
#pragma once

#include <vector>

#include <gtkmm.h>

#include "bayerpreprocess.h"
#include "bayerprocess.h"
#include "bayerrawexposure.h"
#include "blackwhite.h"
#include "cacorrection.h"
#include "chmixer.h"
#include "coarsepanel.h"
#include "colorappearance.h"
#include "colortoning.h"
#include "crop.h"
#include "darkframe.h"
#include "defringe.h"
#include "dehaze.h"
#include "dirpyrdenoise.h"
#include "dirpyrequalizer.h"
#include "distortion.h"
#include "epd.h"
#include "fattaltonemap.h"
#include "filmnegative.h"
#include "filmsimulation.h"
#include "flatfield.h"
#include "gradient.h"
#include "guiutils.h"
#include "hsvequalizer.h"
#include "icmpanel.h"
#include "imageareatoollistener.h"
#include "impulsedenoise.h"
#include "labcurve.h"
#include "lensgeom.h"
#include "lensgeomlistener.h"
#include "lensprofile.h"
#include "localcontrast.h"
#include "locallab.h"
#include "pcvignette.h"
#include "pdsharpening.h"
#include "perspective.h"
#include "pparamschangelistener.h"
#include "preprocess.h"
#include "preprocesswb.h"
#include "profilechangelistener.h"
#include "prsharpening.h"
#include "rawcacorrection.h"
#include "rawexposure.h"
#include "resize.h"
#include "retinex.h"
#include "rgbcurves.h"
#include "rotate.h"
#include "sensorbayer.h"
#include "sensorxtrans.h"
#include "shadowshighlights.h"
#include "sharpenedge.h"
#include "sharpening.h"
#include "sharpenmicro.h"
#include "softlight.h"
#include "tonecurve.h"
#include "toolbar.h"
#include "toolpanel.h"
#include "vibrance.h"
#include "vignetting.h"
#include "wavelet.h"
#include "whitebalance.h"
#include "xtransprocess.h"
#include "xtransrawexposure.h"

#include "../rtengine/noncopyable.h"
#include "../rtengine/rtengine.h"

class ImageEditorCoordinator;
class MetaDataPanel;

class ToolPanelCoordinator :
    public ToolPanelListener,
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
    public rtengine::ImageTypeListener,
    public FilmNegProvider,
    public rtengine::NonCopyable
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
    LocalContrast *localContrast;
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
    SoftLight *softlight;
    Dehaze *dehaze;
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
    PreprocessWB* preprocessWB;
    BayerRAWExposure* bayerrawexposure;
    XTransRAWExposure* xtransrawexposure;
    FattalToneMapping *fattal;
    MetaDataPanel* metadata;
    FilmNegative* filmNegative;
    PdSharpening* pdSharpening;
    std::vector<PParamsChangeListener*> paramcListeners;

    rtengine::StagedImageProcessor* ipc;

    std::vector<ToolPanel*> toolPanels;
    std::vector<FoldableToolPanel*> favorites;
    ToolVBox* favoritePanel;
    ToolVBox* exposurePanel;
    ToolVBox* detailsPanel;
    ToolVBox* colorPanel;
    ToolVBox* transformPanel;
    ToolVBox* rawPanel;
    ToolVBox* advancedPanel;
    ToolVBox* locallabPanel;
    ToolBar* toolBar;

    TextOrIcon* toiF;
    TextOrIcon* toiE;
    TextOrIcon* toiD;
    TextOrIcon* toiC;
    TextOrIcon* toiT;
    TextOrIcon* toiR;
    TextOrIcon* toiM;
    TextOrIcon* toiW;
    TextOrIcon* toiL;

    Gtk::Image* imgPanelEnd[8];
    Gtk::VBox* vbPanelEnd[8];

    Gtk::ScrolledWindow* favoritePanelSW;
    Gtk::ScrolledWindow* exposurePanelSW;
    Gtk::ScrolledWindow* detailsPanelSW;
    Gtk::ScrolledWindow* colorPanelSW;
    Gtk::ScrolledWindow* transformPanelSW;
    Gtk::ScrolledWindow* rawPanelSW;
    Gtk::ScrolledWindow* advancedPanelSW;
    Gtk::ScrolledWindow* locallabPanelSW;

    std::vector<MyExpander*> expList;

    bool hasChanged;

    void addPanel(Gtk::Box* where, FoldableToolPanel* panel, int level = 1);
    void foldThemAll(GdkEventButton* event);
    void updateVScrollbars(bool hide);
    void addfavoritePanel (Gtk::Box* where, FoldableToolPanel* panel, int level = 1);
    void notebookPageChanged(Gtk::Widget* page, guint page_num);

private:
    EditDataProvider *editDataProvider;
    sigc::connection notebookconn;
    bool photoLoadedOnce; // Used to indicated that a photo has been loaded yet
    Gtk::Widget* prevPage;

public:
    CoarsePanel* coarse;
    Gtk::Notebook* toolPanelNotebook;

    ToolPanelCoordinator(bool batch = false);
    ~ToolPanelCoordinator () override;

    bool getChangedState()
    {
        return hasChanged;
    }
    void updateCurveBackgroundHistogram(
        const LUTu& histToneCurve,
        const LUTu& histLCurve,
        const LUTu& histCCurve,
        const LUTu& histLCAM,
        const LUTu& histCCAM,
        const LUTu& histRed,
        const LUTu& histGreen,
        const LUTu& histBlue,
        const LUTu& histLuma,
        const LUTu& histLRETI
    );
    void foldAllButOne(Gtk::Box* parent, FoldableToolPanel* openedSection);

    // multiple listeners can be added that are notified on changes (typical: profile panel and the history)
    void addPParamsChangeListener(PParamsChangeListener* pp)
    {
        paramcListeners.push_back(pp);
    }

    // toolpanellistener interface
    void panelChanged(const rtengine::ProcEvent& event, const Glib::ustring& descr) override;

    void imageTypeChanged (bool isRaw, bool isBayer, bool isXtrans, bool isMono = false) override;

    // profilechangelistener interface
    void profileChange(
        const rtengine::procparams::PartialProfile* nparams,
        const rtengine::ProcEvent& event,
        const Glib::ustring& descr,
        const ParamsEdited* paramsEdited = nullptr,
        bool fromLastSave = false
    ) override;
    void setDefaults(const rtengine::procparams::ProcParams* defparams) override;

    // DirSelectionListener interface
    void dirSelected(const Glib::ustring& dirname, const Glib::ustring& openfile);

    // to support the GUI:
    CropGUIListener* getCropGUIListener();  // through the CropGUIListener the editor area can notify the "crop" ToolPanel when the crop selection changes

    // init the toolpanelcoordinator with an image & close it
    void initImage(rtengine::StagedImageProcessor* ipc_, bool israw);
    void closeImage();

    // update the "expanded" state of the Tools
    void updateToolState();
    void openAllTools();
    void closeAllTools();
    // read/write the "expanded" state of the expanders & read/write the crop panel settings (ratio, guide type, etc.)
    void readOptions();
    void writeOptions();
    void writeToolExpandedStatus(std::vector<int> &tpOpen);
    void updateShowtooltipVisibility (bool showtooltip);

    // wbprovider interface
    void getAutoWB (double& temp, double& green, double equal, double tempBias) override
    {
        if (ipc) {
            ipc->getAutoWB(temp, green, equal, tempBias);
        }
    }
    void getCamWB (double& temp, double& green) override
    {
        if (ipc) {
            ipc->getCamWB(temp, green);
        }
    }

    //DFProvider interface
    rtengine::RawImage* getDF() override;

    //FFProvider interface
    rtengine::RawImage* getFF() override;
    Glib::ustring GetCurrentImageFilePath() override;

    // FilmNegProvider interface
    bool getFilmNegativeExponents(rtengine::Coord spotA, rtengine::Coord spotB, std::array<float, 3>& newExps) override;
    bool getRawSpotValues(rtengine::Coord spot, int spotSize, std::array<float, 3>& rawValues) override;

    // rotatelistener interface
    void straightenRequested () override;
    void autoCropRequested () override;
    void autoPerspRequested (bool corr_pitch, bool corr_yaw, double& rot, double& pitch, double& yaw) override;
    double autoDistorRequested () override;

    // spotwblistener interface
    void spotWBRequested (int size) override;

    // croppanellistener interface
    void cropSelectRequested () override;

    // icmpanellistener interface
    void saveInputICCReference(const Glib::ustring& fname, bool apply_wb) override;

    // imageareatoollistener interface
    void spotWBselected(int x, int y, Thumbnail* thm = nullptr) override;
    void sharpMaskSelected(bool sharpMask) override final;
    int getSpotWBRectSize() const override;
    void cropSelectionReady() override;
    void rotateSelectionReady(double rotate_deg, Thumbnail* thm = nullptr) override;
    ToolBar* getToolBar() const final;
    CropGUIListener* startCropEditing(Thumbnail* thm = nullptr) override;

    void updateTPVScrollbar(bool hide);
    bool handleShortcutKey(GdkEventKey* event);

    // ToolBarListener interface
    void toolSelected (ToolMode tool) override;
    void editModeSwitchedOff () final;

    void setEditProvider(EditDataProvider *provider);

private:
    IdleRegister idle_register;
};
