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
#include "multilangmgr.h"
#include "toolpanelcoord.h"
#include "metadatapanel.h"
#include "options.h"
#include "rtimage.h"

#include "../rtengine/imagesource.h"
#include "../rtengine/dfmanager.h"
#include "../rtengine/ffmanager.h"
#include "../rtengine/improcfun.h"
#include "../rtengine/perspectivecorrection.h"
#include "../rtengine/procevents.h"
#include "../rtengine/refreshmap.h"

using namespace rtengine::procparams;

ToolPanelCoordinator::ToolPanelCoordinator (bool batch) : ipc (nullptr), favoritePanelSW(nullptr), hasChanged (false), editDataProvider (nullptr), photoLoadedOnce(false)
{

    favoritePanel   = Gtk::manage (new ToolVBox ());
    exposurePanel   = Gtk::manage (new ToolVBox ());
    detailsPanel    = Gtk::manage (new ToolVBox ());
    colorPanel      = Gtk::manage (new ToolVBox ());
    transformPanel  = Gtk::manage (new ToolVBox ());
    rawPanel        = Gtk::manage (new ToolVBox ());
    advancedPanel    = Gtk::manage (new ToolVBox ());
    locallabPanel    = Gtk::manage(new ToolVBox());

    coarse              = Gtk::manage (new CoarsePanel ());
    toneCurve           = Gtk::manage (new ToneCurve ());
    shadowshighlights   = Gtk::manage (new ShadowsHighlights ());
    impulsedenoise      = Gtk::manage (new ImpulseDenoise ());
    defringe            = Gtk::manage (new Defringe ());
    dirpyrdenoise       = Gtk::manage (new DirPyrDenoise ());
    epd                 = Gtk::manage (new EdgePreservingDecompositionUI ());
    sharpening          = Gtk::manage (new Sharpening ());
    localContrast       = Gtk::manage(new LocalContrast());
    sharpenEdge         = Gtk::manage(new SharpenEdge());
    sharpenMicro        = Gtk::manage(new SharpenMicro());
    lcurve              = Gtk::manage(new LCurve());
    rgbcurves           = Gtk::manage(new RGBCurves());
    colortoning         = Gtk::manage(new ColorToning());
    lensgeom            = Gtk::manage(new LensGeometry());
    lensProf            = Gtk::manage(new LensProfilePanel());
    distortion          = Gtk::manage(new Distortion());
    rotate              = Gtk::manage(new Rotate());
    vibrance            = Gtk::manage(new Vibrance());
    colorappearance     = Gtk::manage(new ColorAppearance());
    whitebalance        = Gtk::manage(new WhiteBalance());
    vignetting          = Gtk::manage(new Vignetting());
    retinex             = Gtk::manage(new Retinex());
    gradient            = Gtk::manage(new Gradient());
    locallab            = Gtk::manage(new Locallab());
    pcvignette          = Gtk::manage(new PCVignette());
    perspective         = Gtk::manage(new PerspCorrection());
    cacorrection        = Gtk::manage(new CACorrection());
    chmixer             = Gtk::manage(new ChMixer());
    blackwhite          = Gtk::manage(new BlackWhite());
    resize              = Gtk::manage(new Resize());
    prsharpening        = Gtk::manage(new PrSharpening());
    crop                = Gtk::manage(new Crop());
    icm                 = Gtk::manage(new ICMPanel());
    metadata            = Gtk::manage(new MetaDataPanel());
    wavelet             = Gtk::manage(new Wavelet());
    dirpyrequalizer     = Gtk::manage(new DirPyrEqualizer());
    hsvequalizer        = Gtk::manage(new HSVEqualizer());
    filmSimulation      = Gtk::manage(new FilmSimulation());
    softlight           = Gtk::manage(new SoftLight());
    dehaze              = Gtk::manage(new Dehaze());
    sensorbayer         = Gtk::manage(new SensorBayer());
    sensorxtrans        = Gtk::manage(new SensorXTrans());
    bayerprocess        = Gtk::manage(new BayerProcess());
    xtransprocess       = Gtk::manage(new XTransProcess());
    bayerpreprocess     = Gtk::manage(new BayerPreProcess());
    preprocess          = Gtk::manage(new PreProcess());
    darkframe           = Gtk::manage(new DarkFrame());
    flatfield           = Gtk::manage(new FlatField());
    rawcacorrection     = Gtk::manage(new RAWCACorr());
    rawexposure         = Gtk::manage(new RAWExposure());
    preprocessWB        = Gtk::manage (new PreprocessWB ());
    bayerrawexposure    = Gtk::manage(new BayerRAWExposure());
    xtransrawexposure   = Gtk::manage(new XTransRAWExposure());
    fattal              = Gtk::manage(new FattalToneMapping());
    filmNegative        = Gtk::manage (new FilmNegative ());
    pdSharpening        = Gtk::manage (new PdSharpening());
    // So Demosaic, Line noise filter, Green Equilibration, Ca-Correction (garder le nom de section identique!) and Black-Level will be moved in a "Bayer sensor" tool,
    // and a separate Demosaic and Black Level tool will be created in an "X-Trans sensor" tool

    // X-Trans demozaic methods: "3-pass (best), 1-pass (medium), fast"
    // Mettre  jour les profils fournis pour inclure les nouvelles section Raw, notamment pour "Default High ISO"
    // Valeurs par dfaut:
    //     Best -> low ISO
    //     Medium -> High ISO
    favorites.resize(options.favorites.size(), nullptr);

    addfavoritePanel (colorPanel, whitebalance);
    addfavoritePanel (exposurePanel, toneCurve);
    addfavoritePanel (colorPanel, vibrance);
    addfavoritePanel (colorPanel, chmixer);
    addfavoritePanel (colorPanel, blackwhite);
    addfavoritePanel (exposurePanel, shadowshighlights);
    addfavoritePanel (detailsPanel, sharpening);
    addfavoritePanel (detailsPanel, localContrast);
    addfavoritePanel (detailsPanel, sharpenEdge);
    addfavoritePanel (detailsPanel, sharpenMicro);
    addfavoritePanel (colorPanel, hsvequalizer);
    addfavoritePanel (colorPanel, filmSimulation);
    addfavoritePanel (colorPanel, softlight);
    addfavoritePanel (colorPanel, rgbcurves);
    addfavoritePanel (colorPanel, colortoning);
    addfavoritePanel (exposurePanel, epd);
    addfavoritePanel (exposurePanel, fattal);
    addfavoritePanel (advancedPanel, retinex);
    addfavoritePanel (exposurePanel, pcvignette);
    addfavoritePanel (exposurePanel, gradient);
    addfavoritePanel (exposurePanel, lcurve);
    addfavoritePanel (advancedPanel, colorappearance);
    addfavoritePanel (detailsPanel, impulsedenoise);
    addfavoritePanel (detailsPanel, dirpyrdenoise);
    addfavoritePanel (detailsPanel, defringe);
    addfavoritePanel (detailsPanel, dirpyrequalizer);
    addfavoritePanel (detailsPanel, dehaze);
    addfavoritePanel (advancedPanel, wavelet);
    addfavoritePanel(locallabPanel, locallab);
    
    addfavoritePanel (transformPanel, crop);
    addfavoritePanel (transformPanel, resize);
    addPanel (resize->getPackBox(), prsharpening, 2);
    addfavoritePanel (transformPanel, lensgeom);
    addfavoritePanel (lensgeom->getPackBox(), rotate, 2);
    addfavoritePanel (lensgeom->getPackBox(), perspective, 2);
    addfavoritePanel (lensgeom->getPackBox(), lensProf, 2);
    addfavoritePanel (lensgeom->getPackBox(), distortion, 2);
    addfavoritePanel (lensgeom->getPackBox(), cacorrection, 2);
    addfavoritePanel (lensgeom->getPackBox(), vignetting, 2);
    addfavoritePanel (colorPanel, icm);
    addfavoritePanel (rawPanel, sensorbayer);
    addfavoritePanel (sensorbayer->getPackBox(), bayerprocess, 2);
    addfavoritePanel (sensorbayer->getPackBox(), bayerrawexposure, 2);
    addfavoritePanel (sensorbayer->getPackBox(), bayerpreprocess, 2);
    addfavoritePanel (sensorbayer->getPackBox(), rawcacorrection, 2);
    addfavoritePanel (rawPanel, sensorxtrans);
    addfavoritePanel (sensorxtrans->getPackBox(), xtransprocess, 2);
    addfavoritePanel (sensorxtrans->getPackBox(), xtransrawexposure, 2);
    addfavoritePanel (rawPanel, rawexposure);
    addfavoritePanel (rawPanel, preprocessWB);
    addfavoritePanel (rawPanel, preprocess);
    addfavoritePanel (rawPanel, darkframe);
    addfavoritePanel (rawPanel, flatfield);
    addfavoritePanel (rawPanel, filmNegative);
    addfavoritePanel (rawPanel, pdSharpening);

    int favoriteCount = 0;
    for(auto it = favorites.begin(); it != favorites.end(); ++it) {
        if (*it) {
            addPanel(favoritePanel, *it);
            ++favoriteCount;
        }
    }

    toolPanels.push_back (coarse);
    toolPanels.push_back(metadata);

    toolPanelNotebook = new Gtk::Notebook();
    toolPanelNotebook->set_name("ToolPanelNotebook");
    exposurePanelSW    = Gtk::manage (new MyScrolledWindow ());
    detailsPanelSW     = Gtk::manage (new MyScrolledWindow ());
    colorPanelSW       = Gtk::manage (new MyScrolledWindow ());
    transformPanelSW   = Gtk::manage (new MyScrolledWindow ());
    rawPanelSW         = Gtk::manage (new MyScrolledWindow ());
    advancedPanelSW    = Gtk::manage (new MyScrolledWindow ());
    locallabPanelSW     = Gtk::manage(new MyScrolledWindow());    

    // load panel endings
    for (int i = 0; i < 8; i++) {
        vbPanelEnd[i] = Gtk::manage (new Gtk::VBox ());
        imgPanelEnd[i] = Gtk::manage (new RTImage ("ornament1.png"));
        imgPanelEnd[i]->show();
        vbPanelEnd[i]->pack_start(*imgPanelEnd[i], Gtk::PACK_SHRINK);
        vbPanelEnd[i]->show_all();
    }
    if(favoriteCount > 0) {
        favoritePanelSW = Gtk::manage(new MyScrolledWindow());
        favoritePanelSW->add(*favoritePanel);
        favoritePanel->pack_start(*Gtk::manage(new Gtk::HSeparator), Gtk::PACK_SHRINK, 0);
        favoritePanel->pack_start(*vbPanelEnd[0], Gtk::PACK_SHRINK, 4);
    }
    updateVScrollbars(options.hideTPVScrollbar);

    exposurePanelSW->add  (*exposurePanel);
    exposurePanel->pack_start (*Gtk::manage (new Gtk::HSeparator), Gtk::PACK_SHRINK, 0);
    exposurePanel->pack_start (*vbPanelEnd[1], Gtk::PACK_SHRINK, 4);

    detailsPanelSW->add   (*detailsPanel);
    detailsPanel->pack_start (*Gtk::manage (new Gtk::HSeparator), Gtk::PACK_SHRINK, 0);
    detailsPanel->pack_start (*vbPanelEnd[2], Gtk::PACK_SHRINK, 4);

    colorPanelSW->add     (*colorPanel);
    colorPanel->pack_start (*Gtk::manage (new Gtk::HSeparator), Gtk::PACK_SHRINK, 0);
    colorPanel->pack_start (*vbPanelEnd[3], Gtk::PACK_SHRINK, 4);

    advancedPanelSW->add       (*advancedPanel);
    advancedPanel->pack_start (*Gtk::manage (new Gtk::HSeparator), Gtk::PACK_SHRINK, 0);
    advancedPanel->pack_start (*vbPanelEnd[6], Gtk::PACK_SHRINK, 0);

    locallabPanelSW->add(*locallabPanel);
    locallabPanel->pack_start(*Gtk::manage(new Gtk::HSeparator), Gtk::PACK_SHRINK, 0);
    locallabPanel->pack_start(*vbPanelEnd[7], Gtk::PACK_SHRINK, 4);
    
    transformPanelSW->add (*transformPanel);
    transformPanel->pack_start (*Gtk::manage (new Gtk::HSeparator), Gtk::PACK_SHRINK, 0);
    transformPanel->pack_start (*vbPanelEnd[4], Gtk::PACK_SHRINK, 4);

    rawPanelSW->add       (*rawPanel);
    rawPanel->pack_start (*Gtk::manage (new Gtk::HSeparator), Gtk::PACK_SHRINK, 0);
    rawPanel->pack_start (*vbPanelEnd[5], Gtk::PACK_SHRINK, 0);

    toiF = Gtk::manage (new TextOrIcon ("star.png", M ("MAIN_TAB_FAVORITES"), M ("MAIN_TAB_FAVORITES_TOOLTIP")));
    toiE = Gtk::manage (new TextOrIcon ("exposure.png", M ("MAIN_TAB_EXPOSURE"), M ("MAIN_TAB_EXPOSURE_TOOLTIP")));
    toiD = Gtk::manage (new TextOrIcon ("detail.png", M ("MAIN_TAB_DETAIL"), M ("MAIN_TAB_DETAIL_TOOLTIP")));
    toiC = Gtk::manage (new TextOrIcon ("color-circles.png", M ("MAIN_TAB_COLOR"), M ("MAIN_TAB_COLOR_TOOLTIP")));
    toiW = Gtk::manage (new TextOrIcon ("atom.png", M ("MAIN_TAB_ADVANCED"), M ("MAIN_TAB_ADVANCED_TOOLTIP")));
    toiL = Gtk::manage(new TextOrIcon("hand-open.png", M("MAIN_TAB_LOCALLAB"), M("MAIN_TAB_LOCALLAB_TOOLTIP")));

    toiT = Gtk::manage (new TextOrIcon ("transform.png", M ("MAIN_TAB_TRANSFORM"), M ("MAIN_TAB_TRANSFORM_TOOLTIP")));
    toiR = Gtk::manage (new TextOrIcon ("bayer.png", M ("MAIN_TAB_RAW"), M ("MAIN_TAB_RAW_TOOLTIP")));
    toiM = Gtk::manage (new TextOrIcon ("metadata.png", M ("MAIN_TAB_METADATA"), M ("MAIN_TAB_METADATA_TOOLTIP")));
    if (favoritePanelSW) {
        toolPanelNotebook->append_page (*favoritePanelSW,  *toiF);
    }
    toolPanelNotebook->append_page (*exposurePanelSW,  *toiE);
    toolPanelNotebook->append_page (*detailsPanelSW,   *toiD);
    toolPanelNotebook->append_page (*colorPanelSW,     *toiC);
    toolPanelNotebook->append_page (*advancedPanelSW,   *toiW);

    // Locallab notebook is hidden in batch mode
    if (!batch) {
        toolPanelNotebook->append_page(*locallabPanelSW,   *toiL);
    }

    toolPanelNotebook->append_page (*transformPanelSW, *toiT);
    toolPanelNotebook->append_page (*rawPanelSW,       *toiR);
    toolPanelNotebook->append_page (*metadata,    *toiM);

    toolPanelNotebook->set_scrollable();
    toolPanelNotebook->show_all();

    notebookconn = toolPanelNotebook->signal_switch_page().connect(
                       sigc::mem_fun(*this, &ToolPanelCoordinator::notebookPageChanged));

    // In batch mode, notebookPageChanged method is blocked because it's useless to display spots
    if (batch) {
        notebookconn.block(true);
    }

    for (auto toolPanel : toolPanels) {
        toolPanel->setListener(this);
    }

    whitebalance->setWBProvider(this);
    whitebalance->setSpotWBListener(this);
    darkframe->setDFProvider(this);
    flatfield->setFFProvider(this);
    lensgeom->setLensGeomListener(this);
    rotate->setLensGeomListener(this);
    perspective->setLensGeomListener(this);
    distortion->setLensGeomListener(this);
    crop->setCropPanelListener(this);
    icm->setICMPanelListener(this);
    filmNegative->setFilmNegProvider (this);

    toolBar = new ToolBar();
    toolBar->setToolBarListener(this);

    prevPage = toolPanelNotebook->get_nth_page(0);
}

void ToolPanelCoordinator::notebookPageChanged(Gtk::Widget* page, guint page_num)
{
    // Locallab spot curves are set visible if at least one photo has been loaded (to avoid
    // segfault) and locallab panel is active
    if (photoLoadedOnce) {
        if (page == locallabPanelSW) {
            toolBar->blockEditDeactivation(); // Avoid edit tool deactivation when Locallab page is active (except if pressing other tools button)
            locallab->subscribe();
        }

        if (prevPage == locallabPanelSW) { // To deactivate Locallab only when switching from Locallab page
            toolBar->blockEditDeactivation(false);
            locallab->unsubscribe();
        }

        prevPage = page;
    }
}

void ToolPanelCoordinator::addPanel(Gtk::Box* where, FoldableToolPanel* panel, int level)
{

    panel->setParent(where);
    panel->setLevel(level);

    expList.push_back(panel->getExpander());
    where->pack_start(*panel->getExpander(), false, false);
    toolPanels.push_back(panel);
}
void ToolPanelCoordinator::addfavoritePanel (Gtk::Box* where, FoldableToolPanel* panel, int level)
{
    auto name = panel->getToolName();
    auto it = std::find(options.favorites.begin(), options.favorites.end(), name);
    if (it != options.favorites.end()) {
        int index = std::distance(options.favorites.begin(), it);
        favorites[index] = panel;
    } else {
        addPanel(where, panel, level);
    }
}

ToolPanelCoordinator::~ToolPanelCoordinator ()
{
    idle_register.destroy();

    closeImage();

    // When deleting toolPanelNotebook, pages removal activates notebookPageChanged function
    // which is responsible of segfault if listener isn't deactivated before
    notebookconn.block(true);

    delete toolPanelNotebook;
    delete toolBar;
}

void ToolPanelCoordinator::imageTypeChanged(bool isRaw, bool isBayer, bool isXtrans, bool isMono)
{
    if (isRaw) {
        if (isBayer) {
            idle_register.add(
                [this]() -> bool
                {
                    rawPanelSW->set_sensitive(true);
                    sensorxtrans->FoldableToolPanel::hide();
                    xtransprocess->FoldableToolPanel::hide();
                    xtransrawexposure->FoldableToolPanel::hide();
                    sensorbayer->FoldableToolPanel::show();
                    bayerprocess->FoldableToolPanel::show();
                    bayerpreprocess->FoldableToolPanel::show();
                    rawcacorrection->FoldableToolPanel::show();
                    preprocess->FoldableToolPanel::show();
                    flatfield->FoldableToolPanel::show();
                    filmNegative->FoldableToolPanel::show();
                    pdSharpening->FoldableToolPanel::show();
                    retinex->FoldableToolPanel::setGrayedOut(false);
                    return false;
                }
            );
        } else if (isXtrans) {
            idle_register.add(
                [this]() -> bool
                {
                    rawPanelSW->set_sensitive(true);
                    sensorxtrans->FoldableToolPanel::show();
                    xtransprocess->FoldableToolPanel::show();
                    xtransrawexposure->FoldableToolPanel::show();
                    sensorbayer->FoldableToolPanel::hide();
                    bayerprocess->FoldableToolPanel::hide();
                    bayerpreprocess->FoldableToolPanel::hide();
                    rawcacorrection->FoldableToolPanel::hide();
                    preprocess->FoldableToolPanel::show();
                    flatfield->FoldableToolPanel::show();
                    filmNegative->FoldableToolPanel::show();
                    pdSharpening->FoldableToolPanel::show();
                    retinex->FoldableToolPanel::setGrayedOut(false);
                    return false;
                }
            );
        } else if (isMono) {
            idle_register.add(
                [this]() -> bool
                {
                    rawPanelSW->set_sensitive(true);
                    sensorbayer->FoldableToolPanel::hide();
                    bayerprocess->FoldableToolPanel::hide();
                    bayerpreprocess->FoldableToolPanel::hide();
                    rawcacorrection->FoldableToolPanel::hide();
                    sensorxtrans->FoldableToolPanel::hide();
                    xtransprocess->FoldableToolPanel::hide();
                    xtransrawexposure->FoldableToolPanel::hide();
                    preprocessWB->FoldableToolPanel::hide();
                    preprocess->FoldableToolPanel::hide();
                    flatfield->FoldableToolPanel::show();
                    filmNegative->FoldableToolPanel::hide();
                    pdSharpening->FoldableToolPanel::show();
                    retinex->FoldableToolPanel::setGrayedOut(false);
                    return false;
                }
            );
        } else {
            idle_register.add(
                [this]() -> bool
                {
                    rawPanelSW->set_sensitive(true);
                    sensorbayer->FoldableToolPanel::hide();
                    bayerprocess->FoldableToolPanel::hide();
                    bayerpreprocess->FoldableToolPanel::hide();
                    rawcacorrection->FoldableToolPanel::hide();
                    sensorxtrans->FoldableToolPanel::hide();
                    xtransprocess->FoldableToolPanel::hide();
                    xtransrawexposure->FoldableToolPanel::hide();
                    preprocessWB->FoldableToolPanel::hide();
                    preprocess->FoldableToolPanel::hide();
                    flatfield->FoldableToolPanel::hide();
                    filmNegative->FoldableToolPanel::hide();
                    pdSharpening->FoldableToolPanel::hide();
                    retinex->FoldableToolPanel::setGrayedOut(false);
                    return false;
                }
            );
        }
    } else {
        idle_register.add(
            [this]() -> bool
            {
                rawPanelSW->set_sensitive(false);
                sensorbayer->FoldableToolPanel::hide();
                bayerprocess->FoldableToolPanel::hide();
                bayerpreprocess->FoldableToolPanel::hide();
                rawcacorrection->FoldableToolPanel::hide();
                sensorxtrans->FoldableToolPanel::hide();
                xtransprocess->FoldableToolPanel::hide();
                xtransrawexposure->FoldableToolPanel::hide();
                preprocessWB->FoldableToolPanel::hide();
                preprocess->FoldableToolPanel::hide();
                flatfield->FoldableToolPanel::hide();
                filmNegative->FoldableToolPanel::hide();
                pdSharpening->FoldableToolPanel::hide();
                retinex->FoldableToolPanel::setGrayedOut(true);
                return false;
            }
        );
    }

}


void ToolPanelCoordinator::panelChanged(const rtengine::ProcEvent& event, const Glib::ustring& descr)
{
    if (!ipc) {
        return;
    }

    int changeFlags = rtengine::RefreshMapper::getInstance()->getAction(event);

    ProcParams* params = ipc->beginUpdateParams();

    for (auto toolPanel : toolPanels) {
        toolPanel->write(params);
    }

    // Compensate rotation on flip
    if (event == rtengine::EvCTHFlip || event == rtengine::EvCTVFlip) {
        if (fabs(params->rotate.degree) > 0.001) {
            params->rotate.degree *= -1;
            changeFlags |= rtengine::RefreshMapper::getInstance()->getAction(rtengine::EvROTDegree);
            rotate->read(params);
        }
    }

    int tr = TR_NONE;

    if (params->coarse.rotate == 90) {
        tr = TR_R90;
    } else if (params->coarse.rotate == 180) {
        tr = TR_R180;
    } else if (params->coarse.rotate == 270) {
        tr = TR_R270;
    }

    // Update "on preview" geometry
    if (event == rtengine::EvPhotoLoaded || event == rtengine::EvProfileChanged || event == rtengine::EvHistoryBrowsed || event == rtengine::EvCTRotate) {
        // updating the "on preview" geometry
        int fw, fh;
        ipc->getInitialImage()->getImageSource()->getFullSize(fw, fh, tr);
        gradient->updateGeometry(params->gradient.centerX, params->gradient.centerY, params->gradient.feather, params->gradient.degree, fw, fh);
    }

    // some transformations make the crop change for convenience
    if (event == rtengine::EvCTHFlip) {
        crop->hFlipCrop();
        crop->write(params);
    } else if (event == rtengine::EvCTVFlip) {
        crop->vFlipCrop();
        crop->write(params);
    } else if (event == rtengine::EvCTRotate) {
        crop->rotateCrop(params->coarse.rotate, params->coarse.hflip, params->coarse.vflip);
        crop->write(params);
        resize->update(params->crop.enabled, params->crop.w, params->crop.h, ipc->getFullWidth(), ipc->getFullHeight());
        resize->write(params);
    } else if (event == rtengine::EvCrop) {
        resize->update(params->crop.enabled, params->crop.w, params->crop.h);
        resize->write(params);
    }

    /*
     * Manage Locallab mask visibility:
     * - Mask preview is updated when choosing a mask preview method
     * - Mask preview is also updated when modifying (to avoid hiding a potentially visible mask combobox):
     *   - Color&Light invers
     *   - Exposure inversex
     *   - Shadow Highlight inverssh
     *   - Soft Light softMethod
     * - Mask preview is stopped when creating, deleting or selecting a spot
     * - Mask preview is also stopped when removing a spot or resetting all mask visibility
     */
    if (event == rtengine::EvlocallabshowmaskMethod) {
        const Locallab::llMaskVisibility maskStruc = locallab->getMaskVisibility();
        ipc->setLocallabMaskVisibility(maskStruc.previewDeltaE, maskStruc.colorMask, maskStruc.colorMaskinv, maskStruc.expMask, maskStruc.expMaskinv,
                maskStruc.SHMask, maskStruc.SHMaskinv, maskStruc.vibMask, maskStruc.softMask,
                maskStruc.blMask, maskStruc.tmMask, maskStruc.retiMask, maskStruc.sharMask,
                maskStruc.lcMask, maskStruc.cbMask, maskStruc.maskMask);
    } else if (event == rtengine::EvLocallabSpotCreated || event == rtengine::EvLocallabSpotSelectedWithMask ||
            event == rtengine::EvLocallabSpotDeleted || event == rtengine::Evlocallabshowreset ||
            event == rtengine::EvlocallabToolRemovedWithRefresh) {
        locallab->resetMaskVisibility();
        ipc->setLocallabMaskVisibility(false, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    }

    ipc->endUpdateParams(changeFlags);    // starts the IPC processing

    hasChanged = true;

    for (auto paramcListener : paramcListeners) {
        paramcListener->procParamsChanged(params, event, descr);
    }

    // Locallab spot curves are set visible if at least one photo has been loaded (to avoid
    // segfault) and locallab panel is active
    // When a new photo is loaded, Locallab spot curves need to be set visible again
const auto func =
    [this]() -> bool
    {
        if (photoLoadedOnce && (toolPanelNotebook->get_nth_page(toolPanelNotebook->get_current_page()) == locallabPanelSW)) {
            locallab->subscribe();
       }

        return false;
    };

if (event == rtengine::EvPhotoLoaded) {
    idle_register.add(func);
}

    photoLoadedOnce = true;

}

void ToolPanelCoordinator::profileChange(
    const PartialProfile* nparams,
    const rtengine::ProcEvent& event,
    const Glib::ustring& descr,
    const ParamsEdited* paramsEdited,
    bool fromLastSave
)
{
    int fw, fh, tr;

    if (!ipc) {
        return;
    }

    ProcParams *params = ipc->beginUpdateParams();
    ProcParams *mergedParams = new ProcParams();

    // Copy the current params as default values for the fusion
    *mergedParams = *params;

    // Reset IPTC values when switching procparams from the History
    if (event == rtengine::EvHistoryBrowsed) {
        mergedParams->iptc.clear();
        mergedParams->exif.clear();
    }

    // And apply the partial profile nparams to mergedParams
    nparams->applyTo(mergedParams, fromLastSave);

    // Derive the effective changes, if it's a profile change, to prevent slow RAW rerendering if not necessary
    bool filterRawRefresh = false;

    if (event != rtengine::EvPhotoLoaded) {
        ParamsEdited pe(true);
        std::vector<rtengine::procparams::ProcParams> lParams(2);
        lParams[0] = *params;
        lParams[1] = *mergedParams;
        pe.initFrom(lParams);

        filterRawRefresh = pe.raw.isUnchanged() && pe.lensProf.isUnchanged() && pe.retinex.isUnchanged() && pe.filmNegative.isUnchanged() && pe.pdsharpening.isUnchanged();
    }

    *params = *mergedParams;
    delete mergedParams;

    tr = TR_NONE;

    if (params->coarse.rotate == 90) {
        tr = TR_R90;
    } else if (params->coarse.rotate == 180) {
        tr = TR_R180;
    } else if (params->coarse.rotate == 270) {
        tr = TR_R270;
    }

    // trimming overflowing cropped area
    ipc->getInitialImage()->getImageSource()->getFullSize(fw, fh, tr);
    crop->trim(params, fw, fh);

    // updating the GUI with updated values
    for (auto toolPanel : toolPanels) {
        toolPanel->read(params);

        if (event == rtengine::EvPhotoLoaded || event == rtengine::EvProfileChanged) {
            toolPanel->autoOpenCurve();

            // For Locallab, reset tool expanders visibility only when a photo or profile is loaded
            locallab->openAllTools();
        }
    }

    if (event == rtengine::EvPhotoLoaded || event == rtengine::EvProfileChanged || event == rtengine::EvHistoryBrowsed || event == rtengine::EvCTRotate) {
        // updating the "on preview" geometry
        gradient->updateGeometry(params->gradient.centerX, params->gradient.centerY, params->gradient.feather, params->gradient.degree, fw, fh);
    }

    // Reset Locallab mask visibility
    locallab->resetMaskVisibility();
    ipc->setLocallabMaskVisibility(false, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

    // start the IPC processing
    if (filterRawRefresh) {
        ipc->endUpdateParams(rtengine::RefreshMapper::getInstance()->getAction(event) & ALLNORAW);
    } else {
        ipc->endUpdateParams(event);
    }

    hasChanged = event != rtengine::EvProfileChangeNotification;

    for (auto paramcListener : paramcListeners) {
        paramcListener->procParamsChanged(params, event, descr);
    }

    // Locallab spot curves are set visible if at least one photo has been loaded (to avoid
    // segfault) and locallab panel is active
    // When a new photo is loaded, Locallab spot curves need to be set visible again
const auto func =
    [this]() -> bool
    {
        if (photoLoadedOnce && (toolPanelNotebook->get_nth_page(toolPanelNotebook->get_current_page()) == locallabPanelSW)) {
            locallab->subscribe();
        }

        return false;
    };

if (event == rtengine::EvPhotoLoaded) {
    idle_register.add(func);
}

    photoLoadedOnce = true;
}

void ToolPanelCoordinator::setDefaults(const ProcParams* defparams)
{
    if (defparams) {
        for (auto toolPanel : toolPanels) {
            toolPanel->setDefaults(defparams);
        }
    }
}

CropGUIListener* ToolPanelCoordinator::getCropGUIListener()
{

    return crop;
}

void ToolPanelCoordinator::initImage(rtengine::StagedImageProcessor* ipc_, bool raw)
{

    ipc = ipc_;
    toneCurve->disableListener();
    toneCurve->enableAll();
    toneCurve->enableListener();

    if (ipc) {
        const rtengine::FramesMetaData* pMetaData = ipc->getInitialImage()->getMetaData();
        metadata->setImageData(pMetaData);

        ipc->setAutoExpListener(toneCurve);
        ipc->setAutoCamListener(colorappearance);
        ipc->setAutoBWListener(blackwhite);
        ipc->setFrameCountListener(bayerprocess);
        ipc->setFlatFieldAutoClipListener (flatfield);
        ipc->setBayerAutoContrastListener (bayerprocess);
        ipc->setXtransAutoContrastListener (xtransprocess);
        ipc->setpdSharpenAutoContrastListener (pdSharpening);
        ipc->setpdSharpenAutoRadiusListener (pdSharpening);
        ipc->setAutoWBListener(whitebalance);
        ipc->setAutoColorTonListener(colortoning);
        ipc->setAutoChromaListener(dirpyrdenoise);
        ipc->setWaveletListener(wavelet);
        ipc->setRetinexListener(retinex);
        ipc->setSizeListener(crop);
        ipc->setSizeListener(resize);
        ipc->setLocallabListener(locallab);
        ipc->setImageTypeListener(this);
        ipc->setFilmNegListener (filmNegative);
        flatfield->setShortcutPath(Glib::path_get_dirname(ipc->getInitialImage()->getFileName()));

        icm->setRawMeta(raw, (const rtengine::FramesData*)pMetaData);
        lensProf->setRawMeta(raw, pMetaData);
        perspective->setMetadata(pMetaData);
    }


    toneCurve->setRaw(raw);
    hasChanged = true;
}


void ToolPanelCoordinator::closeImage()
{

    if (ipc) {
        ipc->stopProcessing();
        ipc = nullptr;
    }
}

void ToolPanelCoordinator::closeAllTools()
{
    for (size_t i = 0; i < options.tpOpen.size(); ++i) {
        if (i < expList.size()) {
            expList[i]->set_expanded(false);
        }
    }
}

void ToolPanelCoordinator::openAllTools()
{
    for (size_t i = 0; i < options.tpOpen.size(); ++i) {
        if (i < expList.size()) {
            expList[i]->set_expanded(true);
        }
    }
}

void ToolPanelCoordinator::updateToolState()
{
    if (options.tpOpen.empty()) {
        for (auto expander : expList) {
            expander->set_expanded(false);
        }

        wavelet->updateToolState({});
        retinex->updateToolState({});

        return;
    }

    for (size_t i = 0; i < options.tpOpen.size(); ++i) {
        if (i < expList.size()) {
            expList[i]->set_expanded(options.tpOpen[i]);
        }
    }

    if (options.tpOpen.size() > expList.size()) {
        const size_t sizeWavelet = options.tpOpen.size() - expList.size();

        std::vector<int> temp;

        for (size_t i = 0; i < sizeWavelet; ++i) {
            temp.push_back(options.tpOpen[i + expList.size()]);
        }

        wavelet->updateToolState(temp);
        retinex->updateToolState(temp);
    }
}

void ToolPanelCoordinator::readOptions()
{

    crop->readOptions();
}

void ToolPanelCoordinator::writeOptions()
{

    crop->writeOptions();

    if (options.autoSaveTpOpen) {
        writeToolExpandedStatus(options.tpOpen);
    }
}


void ToolPanelCoordinator::writeToolExpandedStatus(std::vector<int> &tpOpen)
{
    tpOpen.clear();

    for (size_t i = 0; i < expList.size(); i++) {
        tpOpen.push_back(expList.at(i)->get_expanded());
    }

    wavelet->writeOptions(tpOpen);
    retinex->writeOptions(tpOpen);

}


void ToolPanelCoordinator::updateShowtooltipVisibility (bool showtooltip)
{
    locallab->updateShowtooltipVisibility(showtooltip);
}


void ToolPanelCoordinator::spotWBselected(int x, int y, Thumbnail* thm)
{
    if (!ipc) {
        return;
    }

//    toolBar->setTool (TOOL_HAND);
    int rect = whitebalance->getSize();
    int ww = ipc->getFullWidth();
    int hh = ipc->getFullHeight();

    if (x - rect > 0 && y - rect > 0 && x + rect < ww && y + rect < hh) {
        double temp;
        double green;
        ipc->getSpotWB(x, y, rect, temp, green);
        whitebalance->setWB(temp, green);
    }
}

void ToolPanelCoordinator::sharpMaskSelected(bool sharpMask)
{
    if (!ipc) {
        return;
    }

    ipc->beginUpdateParams();
    ipc->endUpdateParams (ipc->setSharpMask(sharpMask));
}

int ToolPanelCoordinator::getSpotWBRectSize() const
{
    return whitebalance->getSize();
}

void ToolPanelCoordinator::cropSelectionReady()
{
    toolBar->setTool (TMHand);

    if (!ipc) {
        return;
    }
}

void ToolPanelCoordinator::rotateSelectionReady(double rotate_deg, Thumbnail* thm)
{
    toolBar->setTool (TMHand);

    if (!ipc) {
        return;
    }

    if (rotate_deg != 0.0) {
        rotate->straighten (rotate_deg);
    }
}

ToolBar* ToolPanelCoordinator::getToolBar() const
{
    return toolBar;
}

CropGUIListener* ToolPanelCoordinator::startCropEditing(Thumbnail* thm)
{
    return crop;
}

void ToolPanelCoordinator::autoCropRequested()
{

    if (!ipc) {
        return;
    }

    int x1, y1, x2, y2, w, h;
    ipc->getAutoCrop(crop->getRatio(), x1, y1, w, h);
    x2 = x1 + w - 1;
    y2 = y1 + h - 1;
    crop->cropInit(x1, y1, w, h);
    crop->cropResized(x1, y1, x2, y2);
    crop->cropManipReady();
}

rtengine::RawImage* ToolPanelCoordinator::getDF()
{
    if (!ipc) {
        return nullptr;
    }

    const rtengine::FramesMetaData *imd = ipc->getInitialImage()->getMetaData();

    if (imd) {
        int iso = imd->getISOSpeed();
        double shutter = imd->getShutterSpeed();
        std::string maker(imd->getMake());
        std::string model(imd->getModel());
        time_t timestamp = imd->getDateTimeAsTS();

        return rtengine::dfm.searchDarkFrame(maker, model, iso, shutter, timestamp);
    }

    return nullptr;
}

rtengine::RawImage* ToolPanelCoordinator::getFF()
{
    if (!ipc) {
        return nullptr;
    }

    const rtengine::FramesMetaData *imd = ipc->getInitialImage()->getMetaData();

    if (imd) {
        // int iso = imd->getISOSpeed();              temporarily removed because unused
        // double shutter = imd->getShutterSpeed();   temporarily removed because unused
        double aperture = imd->getFNumber();
        double focallength = imd->getFocalLen();
        std::string maker(imd->getMake());
        std::string model(imd->getModel());
        std::string lens(imd->getLens());
        time_t timestamp = imd->getDateTimeAsTS();

        return rtengine::ffm.searchFlatField(maker, model, lens, focallength, aperture, timestamp);
    }

    return nullptr;
}

Glib::ustring ToolPanelCoordinator::GetCurrentImageFilePath()
{
    if (!ipc) {
        return "";
    }

    return ipc->getInitialImage()->getFileName();
}

void ToolPanelCoordinator::straightenRequested()
{

    if (!ipc) {
        return;
    }

    toolBar->setTool(TMStraighten);
}

void ToolPanelCoordinator::autoPerspRequested (bool corr_pitch, bool corr_yaw, double& rot, double& pitch, double& yaw)
{
    if (!(ipc && (corr_pitch || corr_yaw))) {
        return;
    }

    rtengine::ImageSource *src = dynamic_cast<rtengine::ImageSource *>(ipc->getInitialImage());
    if (!src) {
        return;
    }

    rtengine::procparams::ProcParams params;
    ipc->getParams(&params);

    auto res = rtengine::PerspectiveCorrection::autocompute(src, corr_pitch, corr_yaw, &params, src->getMetaData());
    rot = res.angle;
    pitch = res.pitch;
    yaw = res.yaw;
}

double ToolPanelCoordinator::autoDistorRequested()
{
    if (!ipc) {
        return 0.0;
    }

    return rtengine::ImProcFunctions::getAutoDistor(ipc->getInitialImage()->getFileName(), 400);
}

void ToolPanelCoordinator::spotWBRequested(int size)
{

    if (!ipc) {
        return;
    }

    toolBar->setTool(TMSpotWB);
}

void ToolPanelCoordinator::cropSelectRequested()
{

    if (!ipc) {
        return;
    }

    toolBar->setTool(TMCropSelect);
}

void ToolPanelCoordinator::saveInputICCReference(const Glib::ustring& fname, bool apply_wb)
{
    if (ipc) {
        ipc->saveInputICCReference(fname, apply_wb);
    }
}

void ToolPanelCoordinator::updateCurveBackgroundHistogram(
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
)
{
    colorappearance->updateCurveBackgroundHistogram(histToneCurve, histLCurve, histCCurve, histLCAM,  histCCAM, histRed, histGreen, histBlue, histLuma, histLRETI);
    toneCurve->updateCurveBackgroundHistogram(histToneCurve, histLCurve, histCCurve,histLCAM,  histCCAM, histRed, histGreen, histBlue, histLuma, histLRETI);
    lcurve->updateCurveBackgroundHistogram(histToneCurve, histLCurve, histCCurve, histLCAM, histCCAM, histRed, histGreen, histBlue, histLuma, histLRETI);
    rgbcurves->updateCurveBackgroundHistogram(histToneCurve, histLCurve, histCCurve, histLCAM, histCCAM, histRed, histGreen, histBlue, histLuma, histLRETI);
    retinex->updateCurveBackgroundHistogram(histToneCurve, histLCurve, histCCurve, histLCAM, histCCAM, histRed, histGreen, histBlue, histLuma, histLRETI);
}

void ToolPanelCoordinator::foldAllButOne(Gtk::Box* parent, FoldableToolPanel* openedSection)
{

    for (auto toolPanel : toolPanels) {
        if (toolPanel->getParent() != nullptr) {
            ToolPanel* currentTP = toolPanel;

            if (currentTP->getParent() == parent) {
                // Section in the same tab, we unfold it if it's not the one that has been clicked
                if (currentTP != openedSection) {
                    currentTP->setExpanded(false);
                } else {
                    if (!currentTP->getExpanded()) {
                        currentTP->setExpanded(true);
                    }
                }
            }
        }
    }
}

bool ToolPanelCoordinator::handleShortcutKey(GdkEventKey* event)
{

    //bool ctrl = event->state & GDK_CONTROL_MASK;  temporarily removed because unused
    //bool shift = event->state & GDK_SHIFT_MASK;   temporarily removed because unused
    bool alt = event->state & GDK_MOD1_MASK;

    if (alt) {
        switch (event->keyval) {
            case GDK_KEY_u:
                if (favoritePanelSW) {
                    toolPanelNotebook->set_current_page (toolPanelNotebook->page_num (*favoritePanelSW));
                }
                return true;

            case GDK_KEY_e:
                toolPanelNotebook->set_current_page(toolPanelNotebook->page_num(*exposurePanelSW));
                return true;

            case GDK_KEY_d:
                toolPanelNotebook->set_current_page(toolPanelNotebook->page_num(*detailsPanelSW));
                return true;

            case GDK_KEY_c:
                toolPanelNotebook->set_current_page(toolPanelNotebook->page_num(*colorPanelSW));
                return true;

            case GDK_KEY_t:
                toolPanelNotebook->set_current_page(toolPanelNotebook->page_num(*transformPanelSW));
                return true;

            case GDK_KEY_r:
                toolPanelNotebook->set_current_page(toolPanelNotebook->page_num(*rawPanelSW));
                return true;

            case GDK_KEY_a:
                toolPanelNotebook->set_current_page(toolPanelNotebook->page_num(*advancedPanelSW));
                return true;

            case GDK_KEY_o:
                toolPanelNotebook->set_current_page(toolPanelNotebook->page_num(*locallabPanelSW));
                return true;

            case GDK_KEY_m:
                toolPanelNotebook->set_current_page(toolPanelNotebook->page_num(*metadata));
                return true;
        }
    }

    return false;
}

void ToolPanelCoordinator::updateVScrollbars(bool hide)
{
    GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected
    Gtk::PolicyType policy = hide ? Gtk::POLICY_NEVER : Gtk::POLICY_AUTOMATIC;
    if (favoritePanelSW) {
        favoritePanelSW->set_policy     (Gtk::POLICY_AUTOMATIC, policy);
    }
    exposurePanelSW->set_policy     (Gtk::POLICY_AUTOMATIC, policy);
    detailsPanelSW->set_policy      (Gtk::POLICY_AUTOMATIC, policy);
    colorPanelSW->set_policy        (Gtk::POLICY_AUTOMATIC, policy);
    transformPanelSW->set_policy    (Gtk::POLICY_AUTOMATIC, policy);
    rawPanelSW->set_policy          (Gtk::POLICY_AUTOMATIC, policy);
    advancedPanelSW->set_policy      (Gtk::POLICY_AUTOMATIC, policy);
    locallabPanelSW->set_policy(Gtk::POLICY_AUTOMATIC, policy);
    

    for (auto currExp : expList) {
        currExp->updateVScrollbars(hide);
    }
}


void ToolPanelCoordinator::updateTPVScrollbar(bool hide)
{
    updateVScrollbars(hide);
}

void ToolPanelCoordinator::toolSelected(ToolMode tool)
{
    GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected
    notebookconn.block(true); // "signal_switch_page" event is blocked to avoid unsubscribing Locallab (allows a correct behavior when switching to another tool using toolbar)

    auto checkFavorite = [this](FoldableToolPanel* tool) {
        for (auto fav : favorites) {
            if (fav == tool) {
                return true;
            }
        }
        return false;
    };

    switch (tool) {
        case TMCropSelect: {
            toolBar->blockEditDeactivation(false); // To allow deactivating Locallab when switching to another tool using toolbar
            crop->setExpanded(true);
            toolPanelNotebook->set_current_page(toolPanelNotebook->page_num(checkFavorite(crop) ? *favoritePanelSW : *transformPanelSW));
            prevPage = toolPanelNotebook->get_nth_page(toolPanelNotebook->get_current_page()); // Updating prevPage as "signal_switch_page" event
            break;
        }

        case TMSpotWB: {
            toolBar->blockEditDeactivation(false); // To allow deactivating Locallab when switching to another tool using toolbar
            whitebalance->setExpanded(true);
            toolPanelNotebook->set_current_page(toolPanelNotebook->page_num(checkFavorite(whitebalance) ? *favoritePanelSW : *colorPanelSW));
            prevPage = toolPanelNotebook->get_nth_page(toolPanelNotebook->get_current_page()); // Updating prevPage as "signal_switch_page" event
            break;
        }

        case TMStraighten: {
            toolBar->blockEditDeactivation(false); // To allow deactivating Locallab when switching to another tool using toolbar
            rotate->setExpanded(true);
            bool isFavorite = checkFavorite(rotate);
            if (!isFavorite) {
                isFavorite = checkFavorite(lensgeom);
                lensgeom->setExpanded(true);
            }
            toolPanelNotebook->set_current_page(toolPanelNotebook->page_num(isFavorite ? *favoritePanelSW : *transformPanelSW));
            prevPage = toolPanelNotebook->get_nth_page(toolPanelNotebook->get_current_page()); // Updating prevPage as "signal_switch_page" event
            break;
        }

        default:
            break;
    }

    notebookconn.block(false);
}

void ToolPanelCoordinator::editModeSwitchedOff()
{
    if (editDataProvider) {
        editDataProvider->switchOffEditMode();
    }
}

void ToolPanelCoordinator::dirSelected(const Glib::ustring& dirname, const Glib::ustring& openfile)
{

    flatfield->setShortcutPath(dirname);
}

void ToolPanelCoordinator::setEditProvider(EditDataProvider *provider)
{
    editDataProvider = provider;

    for (size_t i = 0; i < toolPanels.size(); i++) {
        toolPanels.at(i)->setEditProvider(provider);
    }
}

bool ToolPanelCoordinator::getFilmNegativeExponents(rtengine::Coord spotA, rtengine::Coord spotB, std::array<float, 3>& newExps)
{
    return ipc && ipc->getFilmNegativeExponents(spotA.x, spotA.y, spotB.x, spotB.y, newExps);
}

bool ToolPanelCoordinator::getRawSpotValues(rtengine::Coord spot, int spotSize, std::array<float, 3>& rawValues)
{
    return ipc && ipc->getRawSpotValues(spot.x, spot.y, spotSize, rawValues);
}
