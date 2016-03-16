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
#include "multilangmgr.h"
#include "toolpanelcoord.h"
#include "options.h"
#include "../rtengine/imagesource.h"
#include "../rtengine/dfmanager.h"
#include "../rtengine/ffmanager.h"
#include "../rtengine/improcfun.h"
#include "../rtengine/procevents.h"
#include "../rtengine/refreshmap.h"

using namespace rtengine::procparams;

ToolPanelCoordinator::ToolPanelCoordinator () : ipc(nullptr), editDataProvider(nullptr)
{

    exposurePanel   = Gtk::manage (new ToolVBox ());
    detailsPanel    = Gtk::manage (new ToolVBox ());
    colorPanel      = Gtk::manage (new ToolVBox ());
    transformPanel  = Gtk::manage (new ToolVBox ());
    rawPanel        = Gtk::manage (new ToolVBox ());
    waveletPanel    = Gtk::manage (new ToolVBox ());

    coarse              = Gtk::manage (new CoarsePanel ());
    toneCurve           = Gtk::manage (new ToneCurve ());
    shadowshighlights   = Gtk::manage (new ShadowsHighlights ());
    impulsedenoise      = Gtk::manage (new ImpulseDenoise ());
    defringe            = Gtk::manage (new Defringe ());
    dirpyrdenoise       = Gtk::manage (new DirPyrDenoise ());
    epd                 = Gtk::manage (new EdgePreservingDecompositionUI ());
    sharpening          = Gtk::manage (new Sharpening ());
    sharpenEdge         = Gtk::manage (new SharpenEdge ());
    sharpenMicro        = Gtk::manage (new SharpenMicro ());
    lcurve              = Gtk::manage (new LCurve ());
    rgbcurves           = Gtk::manage (new RGBCurves ());
    colortoning         = Gtk::manage (new ColorToning ());
    lensgeom            = Gtk::manage (new LensGeometry ());
    lensProf            = Gtk::manage (new LensProfilePanel ());
    distortion          = Gtk::manage (new Distortion ());
    rotate              = Gtk::manage (new Rotate ());
    vibrance            = Gtk::manage (new Vibrance ());
    colorappearance     = Gtk::manage (new ColorAppearance ());
    whitebalance        = Gtk::manage (new WhiteBalance ());
    vignetting          = Gtk::manage (new Vignetting ());
    retinex               = Gtk::manage (new Retinex ());
    gradient            = Gtk::manage (new Gradient ());
    pcvignette          = Gtk::manage (new PCVignette ());
    perspective         = Gtk::manage (new PerspCorrection ());
    cacorrection        = Gtk::manage (new CACorrection ());
    chmixer             = Gtk::manage (new ChMixer ());
    blackwhite          = Gtk::manage (new BlackWhite ());
    resize              = Gtk::manage (new Resize ());
    prsharpening        = Gtk::manage (new PrSharpening());
    crop                = Gtk::manage (new Crop ());
    icm                 = Gtk::manage (new ICMPanel ());
    exifpanel           = Gtk::manage (new ExifPanel ());
    iptcpanel           = Gtk::manage (new IPTCPanel ());
    wavelet             = Gtk::manage (new Wavelet ());
    dirpyrequalizer     = Gtk::manage (new DirPyrEqualizer ());
    hsvequalizer        = Gtk::manage (new HSVEqualizer ());
    filmSimulation      = Gtk::manage (new FilmSimulation ());
    sensorbayer         = Gtk::manage (new SensorBayer ());
    sensorxtrans        = Gtk::manage (new SensorXTrans ());
    bayerprocess        = Gtk::manage (new BayerProcess ());
    xtransprocess       = Gtk::manage (new XTransProcess ());
    bayerpreprocess     = Gtk::manage (new BayerPreProcess ());
    preprocess          = Gtk::manage (new PreProcess ());
    darkframe           = Gtk::manage (new DarkFrame ());
    flatfield           = Gtk::manage (new FlatField ());
    rawcacorrection     = Gtk::manage (new RAWCACorr ());
    rawexposure         = Gtk::manage (new RAWExposure ());
    bayerrawexposure    = Gtk::manage (new BayerRAWExposure ());
    xtransrawexposure   = Gtk::manage (new XTransRAWExposure ());

    // So Demosaic, Line noise filter, Green Equilibration, Ca-Correction (garder le nom de section identique!) and Black-Level will be moved in a "Bayer sensor" tool,
    // and a separate Demosaic and Black Level tool will be created in an "X-Trans sensor" tool

    // X-Trans demozaic methods: "3-pass (best), 1-pass (medium), fast"
    // Mettre  jour les profils fournis pour inclure les nouvelles section Raw, notamment pour "Default High ISO"
    // Valeurs par dfaut:
    //     Best -> low ISO
    //     Medium -> High ISO

    addPanel (colorPanel, whitebalance);
    toolPanels.push_back (whitebalance);
    addPanel (exposurePanel, toneCurve);
    toolPanels.push_back (toneCurve);
    addPanel (colorPanel, vibrance);
    toolPanels.push_back (vibrance);
    addPanel (colorPanel, chmixer);
    toolPanels.push_back (chmixer); // << TODO: Add "Enabled"
    addPanel (colorPanel, blackwhite);
    toolPanels.push_back (blackwhite);
    addPanel (exposurePanel, shadowshighlights);
    toolPanels.push_back (shadowshighlights);
    addPanel (detailsPanel, sharpening);
    toolPanels.push_back (sharpening);
    addPanel (detailsPanel, sharpenEdge);
    toolPanels.push_back (sharpenEdge);
    addPanel (detailsPanel, sharpenMicro);
    toolPanels.push_back (sharpenMicro);
    addPanel (colorPanel, hsvequalizer);
    toolPanels.push_back (hsvequalizer); // << TODO: Add "Enabled"
    addPanel (colorPanel, filmSimulation);
    toolPanels.push_back (filmSimulation);
    addPanel (colorPanel, rgbcurves);
    toolPanels.push_back (rgbcurves); // << TODO: Add "Enabled"
    addPanel (colorPanel, colortoning);
    toolPanels.push_back (colortoning);
    addPanel (exposurePanel, epd);
    toolPanels.push_back (epd);
    addPanel (exposurePanel, retinex);
    toolPanels.push_back (retinex);
    addPanel (exposurePanel, pcvignette);
    toolPanels.push_back (pcvignette);
    addPanel (exposurePanel, gradient);
    toolPanels.push_back (gradient);
    addPanel (exposurePanel, lcurve);
    toolPanels.push_back (lcurve); // << TODO: Add "Enabled" ???
    addPanel (exposurePanel, colorappearance);
    toolPanels.push_back (colorappearance);
    addPanel (detailsPanel, impulsedenoise);
    toolPanels.push_back (impulsedenoise);
    addPanel (detailsPanel, dirpyrdenoise);
    toolPanels.push_back (dirpyrdenoise);
    addPanel (detailsPanel, defringe);
    toolPanels.push_back (defringe);
    addPanel (detailsPanel, dirpyrequalizer);
    toolPanels.push_back (dirpyrequalizer);
    addPanel (waveletPanel, wavelet);
    toolPanels.push_back (wavelet);
    addPanel (transformPanel, crop);
    toolPanels.push_back (crop);
    addPanel (transformPanel, resize);
    toolPanels.push_back (resize);
    addPanel (resize->getPackBox(), prsharpening, 2);
    toolPanels.push_back (prsharpening);
    addPanel (transformPanel, lensgeom);
    toolPanels.push_back (lensgeom);
    addPanel (lensgeom->getPackBox(), rotate, 2);
    toolPanels.push_back (rotate);
    addPanel (lensgeom->getPackBox(), perspective, 2);
    toolPanels.push_back (perspective);
    addPanel (lensgeom->getPackBox(), lensProf, 2);
    toolPanels.push_back (lensProf);
    addPanel (lensgeom->getPackBox(), distortion, 2);
    toolPanels.push_back (distortion);
    addPanel (lensgeom->getPackBox(), cacorrection, 2);
    toolPanels.push_back (cacorrection);
    addPanel (lensgeom->getPackBox(), vignetting, 2);
    toolPanels.push_back (vignetting);
    addPanel (colorPanel, icm);
    toolPanels.push_back (icm);
    addPanel (rawPanel, sensorbayer);
    toolPanels.push_back (sensorbayer);
    addPanel (sensorbayer->getPackBox(), bayerprocess, 2);
    toolPanels.push_back (bayerprocess);
    addPanel (sensorbayer->getPackBox(), bayerrawexposure, 2);
    toolPanels.push_back (bayerrawexposure);
    addPanel (sensorbayer->getPackBox(), bayerpreprocess, 2);
    toolPanels.push_back (bayerpreprocess);
    addPanel (sensorbayer->getPackBox(), rawcacorrection, 2);
    toolPanels.push_back (rawcacorrection);
    addPanel (rawPanel, sensorxtrans);
    toolPanels.push_back (sensorxtrans);
    addPanel (sensorxtrans->getPackBox(), xtransprocess, 2);
    toolPanels.push_back (xtransprocess);
    addPanel (sensorxtrans->getPackBox(), xtransrawexposure, 2);
    toolPanels.push_back (xtransrawexposure);
    addPanel (rawPanel, rawexposure);
    toolPanels.push_back (rawexposure);
    addPanel (rawPanel, preprocess);
    toolPanels.push_back (preprocess);
    addPanel (rawPanel, darkframe);
    toolPanels.push_back (darkframe);
    addPanel (rawPanel, flatfield);
    toolPanels.push_back (flatfield);

    toolPanels.push_back (coarse);
    toolPanels.push_back (exifpanel);
    toolPanels.push_back (iptcpanel);

    metadataPanel = Gtk::manage (new Gtk::Notebook ());
    metadataPanel->set_name ("MetaPanelNotebook");
    toolPanelNotebook = new Gtk::Notebook ();
    toolPanelNotebook->set_name ("ToolPanelNotebook");

    metadataPanel->append_page (*exifpanel, M("MAIN_TAB_EXIF"));
    metadataPanel->append_page (*iptcpanel, M("MAIN_TAB_IPTC"));

    exposurePanelSW    = Gtk::manage (new MyScrolledWindow ());
    detailsPanelSW     = Gtk::manage (new MyScrolledWindow ());
    colorPanelSW       = Gtk::manage (new MyScrolledWindow ());
    transformPanelSW   = Gtk::manage (new MyScrolledWindow ());
    rawPanelSW         = Gtk::manage (new MyScrolledWindow ());
    waveletPanelSW     = Gtk::manage (new MyScrolledWindow ());
    updateVScrollbars (options.hideTPVScrollbar);

    // load panel endings
    for (int i = 0; i < 6; i++) {
        vbPanelEnd[i] = Gtk::manage (new Gtk::VBox ());
        imgPanelEnd[i] = Gtk::manage (new RTImage("PanelEnding.png"));
        imgPanelEnd[i]->show ();
        vbPanelEnd[i]->pack_start (*imgPanelEnd[i], Gtk::PACK_SHRINK);
        vbPanelEnd[i]->show_all();
    }

    exposurePanelSW->add  (*exposurePanel);
    exposurePanel->pack_start (*Gtk::manage(new Gtk::HSeparator), Gtk::PACK_SHRINK, 0);
    exposurePanel->pack_start (*vbPanelEnd[0], Gtk::PACK_SHRINK, 4);

    detailsPanelSW->add   (*detailsPanel);
    detailsPanel->pack_start (*Gtk::manage(new Gtk::HSeparator), Gtk::PACK_SHRINK, 0);
    detailsPanel->pack_start (*vbPanelEnd[1], Gtk::PACK_SHRINK, 4);

    colorPanelSW->add     (*colorPanel);
    colorPanel->pack_start (*Gtk::manage(new Gtk::HSeparator), Gtk::PACK_SHRINK, 0);
    colorPanel->pack_start (*vbPanelEnd[2], Gtk::PACK_SHRINK, 4);

    waveletPanelSW->add       (*waveletPanel);
    waveletPanel->pack_start (*Gtk::manage(new Gtk::HSeparator), Gtk::PACK_SHRINK, 0);
    waveletPanel->pack_start (*vbPanelEnd[5], Gtk::PACK_SHRINK, 0);

    transformPanelSW->add (*transformPanel);
    transformPanel->pack_start (*Gtk::manage(new Gtk::HSeparator), Gtk::PACK_SHRINK, 0);
    transformPanel->pack_start (*vbPanelEnd[3], Gtk::PACK_SHRINK, 4);

    rawPanelSW->add       (*rawPanel);
    rawPanel->pack_start (*Gtk::manage(new Gtk::HSeparator), Gtk::PACK_SHRINK, 0);
    rawPanel->pack_start (*vbPanelEnd[4], Gtk::PACK_SHRINK, 0);



    TOITypes type = options.UseIconNoText ? TOI_ICON : TOI_TEXT;

    toiE = Gtk::manage (new TextOrIcon ("exposure.png" , M("MAIN_TAB_EXPOSURE") , M("MAIN_TAB_EXPOSURE_TOOLTIP") , type));
    toiD = Gtk::manage (new TextOrIcon ("detail.png"   , M("MAIN_TAB_DETAIL")   , M("MAIN_TAB_DETAIL_TOOLTIP")   , type));
    toiC = Gtk::manage (new TextOrIcon ("colour.png"   , M("MAIN_TAB_COLOR")    , M("MAIN_TAB_COLOR_TOOLTIP")    , type));
    toiW = Gtk::manage (new TextOrIcon ("wavelet.png"  , M("MAIN_TAB_WAVELET")  , M("MAIN_TAB_WAVELET_TOOLTIP") , type));
    toiT = Gtk::manage (new TextOrIcon ("transform.png", M("MAIN_TAB_TRANSFORM"), M("MAIN_TAB_TRANSFORM_TOOLTIP"), type));
    toiR = Gtk::manage (new TextOrIcon ("raw.png"      , M("MAIN_TAB_RAW")      , M("MAIN_TAB_RAW_TOOLTIP")      , type));
    toiM = Gtk::manage (new TextOrIcon ("meta.png"     , M("MAIN_TAB_METADATA") , M("MAIN_TAB_METADATA_TOOLTIP") , type));

    toolPanelNotebook->append_page (*exposurePanelSW,  *toiE);
    toolPanelNotebook->append_page (*detailsPanelSW,   *toiD);
    toolPanelNotebook->append_page (*colorPanelSW,     *toiC);
    toolPanelNotebook->append_page (*waveletPanelSW,   *toiW);
    toolPanelNotebook->append_page (*transformPanelSW, *toiT);
    toolPanelNotebook->append_page (*rawPanelSW,       *toiR);
    toolPanelNotebook->append_page (*metadataPanel,    *toiM);

    toolPanelNotebook->set_current_page (0);

    toolPanelNotebook->set_scrollable ();
    toolPanelNotebook->show_all ();

    for (auto toolPanel : toolPanels) {
        toolPanel->setListener (this);
    }

    whitebalance->setWBProvider (this);
    whitebalance->setSpotWBListener (this);
    darkframe->setDFProvider (this);
    flatfield->setFFProvider (this);
    lensgeom->setLensGeomListener (this);
    rotate->setLensGeomListener (this);
    distortion->setLensGeomListener (this);
    crop->setCropPanelListener (this);
    icm->setICMPanelListener (this);

    toolBar = new ToolBar ();
    toolBar->setToolBarListener(this);
}

void ToolPanelCoordinator::addPanel (Gtk::Box* where, FoldableToolPanel* panel, int level)
{

    panel->setParent(where);
    panel->setLevel(level);

    expList.push_back (panel->getExpander());
    where->pack_start(*panel->getExpander(), false, false);
}

ToolPanelCoordinator::~ToolPanelCoordinator ()
{

    closeImage ();

    delete toolPanelNotebook;
    delete toolBar;
}

void ToolPanelCoordinator::panelChanged (rtengine::ProcEvent event, const Glib::ustring& descr)
{

    if (!ipc) {
        return;
    }

    int changeFlags = refreshmap[(int)event];

    ProcParams* params = ipc->beginUpdateParams ();

    for (auto toolPanel : toolPanels) {
        toolPanel->write (params);
    }

    // Compensate rotation on flip
    if (event == rtengine::EvCTHFlip || event == rtengine::EvCTVFlip) {
        if (fabs(params->rotate.degree) > 0.001) {
            params->rotate.degree *= -1;
            changeFlags |= refreshmap[(int)rtengine::EvROTDegree];
            rotate->read (params);
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
        ipc->getInitialImage()->getImageSource()->getFullSize (fw, fh, tr);
        gradient->updateGeometry (params->gradient.centerX, params->gradient.centerY, params->gradient.feather, params->gradient.degree, fw, fh);
    }

    // some transformations make the crop change for convenience
    if (event == rtengine::EvCTHFlip) {
        crop->hFlipCrop ();
        crop->write (params);
    } else if (event == rtengine::EvCTVFlip) {
        crop->vFlipCrop ();
        crop->write (params);
    } else if (event == rtengine::EvCTRotate) {
        crop->rotateCrop (params->coarse.rotate, params->coarse.hflip, params->coarse.vflip);
        crop->write (params);
        resize->update (params->crop.enabled, params->crop.w, params->crop.h, ipc->getFullWidth(), ipc->getFullHeight());
        resize->write (params);
    } else if (event == rtengine::EvCrop) {
        resize->update (params->crop.enabled, params->crop.w, params->crop.h);
        resize->write (params);
    }

    ipc->endUpdateParams (changeFlags);   // starts the IPC processing

    hasChanged = true;

    for (auto paramcListener : paramcListeners) {
        paramcListener->procParamsChanged (params, event, descr);
    }
}

void ToolPanelCoordinator::profileChange  (const PartialProfile *nparams, rtengine::ProcEvent event, const Glib::ustring& descr, const ParamsEdited* paramsEdited)
{

    int fw, fh, tr;

    if (!ipc) {
        return;
    }

    ProcParams *params = ipc->beginUpdateParams ();
    ProcParams *mergedParams = new ProcParams();

    // Copy the current params as default values for the fusion
    *mergedParams = *params;

    // Reset IPTC values when switching procparams from the History
    if (event == rtengine::EvHistoryBrowsed) {
        mergedParams->iptc.clear();
        mergedParams->exif.clear();
    }

    // And apply the partial profile nparams to mergedParams
    nparams->applyTo(mergedParams);

    // Derive the effective changes, if it's a profile change, to prevent slow RAW rerendering if not necessary
    bool filterRawRefresh = false;

    if (event != rtengine::EvPhotoLoaded) {
        ParamsEdited pe(true);
        std::vector<rtengine::procparams::ProcParams> lParams(2);
        lParams[0] = *params;
        lParams[1] = *mergedParams;
        pe.initFrom (lParams);

        filterRawRefresh = pe.raw.isUnchanged() && pe.lensProf.isUnchanged() && pe.retinex.isUnchanged();
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
    ipc->getInitialImage()->getImageSource()->getFullSize (fw, fh, tr);
    crop->trim(params, fw, fh);

    // updating the GUI with updated values
    for (auto toolPanel : toolPanels) {
        toolPanel->read (params);

        if (event == rtengine::EvPhotoLoaded || event == rtengine::EvProfileChanged) {
            toolPanel->autoOpenCurve();
        }
    }

    if (event == rtengine::EvPhotoLoaded || event == rtengine::EvProfileChanged || event == rtengine::EvHistoryBrowsed || event == rtengine::EvCTRotate) {
        // updating the "on preview" geometry
        gradient->updateGeometry (params->gradient.centerX, params->gradient.centerY, params->gradient.feather, params->gradient.degree, fw, fh);
    }

    // start the IPC processing
    if (filterRawRefresh) {
        ipc->endUpdateParams ( refreshmap[(int)event] & ALLNORAW );
    } else {
        ipc->endUpdateParams (event);
    }

    hasChanged = event != rtengine::EvProfileChangeNotification;

    for (auto paramcListener : paramcListeners) {
        paramcListener->procParamsChanged (params, event, descr);
    }
}

void ToolPanelCoordinator::setDefaults (ProcParams* defparams)
{

    if (defparams)
        for (auto toolPanel : toolPanels) {
            toolPanel->setDefaults (defparams);
        }
}

CropGUIListener* ToolPanelCoordinator::getCropGUIListener ()
{

    return crop;
}

void ToolPanelCoordinator::initImage (rtengine::StagedImageProcessor* ipc_, bool raw)
{

    ipc = ipc_;
    toneCurve->disableListener ();
    toneCurve->enableAll ();
    toneCurve->enableListener ();

    const rtengine::ImageMetaData* pMetaData = ipc->getInitialImage()->getMetaData();
    exifpanel->setImageData (pMetaData);
    iptcpanel->setImageData (pMetaData);

    if (ipc) {
        ipc->setAutoExpListener (toneCurve);
        ipc->setAutoCamListener (colorappearance);
        ipc->setAutoBWListener (blackwhite);
        ipc->setAutoColorTonListener (colortoning);
        ipc->setAutoChromaListener (dirpyrdenoise);
        ipc->setWaveletListener (wavelet);
        ipc->setRetinexListener (retinex);

        ipc->setSizeListener (crop);
        ipc->setSizeListener (resize);
    }

    flatfield->setShortcutPath(Glib::path_get_dirname(ipc->getInitialImage()->getFileName()));

    icm->setRawMeta (raw, (const rtengine::ImageData*)pMetaData);
    lensProf->setRawMeta (raw, pMetaData);

    toneCurve->setRaw (raw);
    hasChanged = true;
}


void ToolPanelCoordinator::closeImage ()
{

    if (ipc) {
        ipc->stopProcessing ();
        ipc = NULL;
    }
}

void ToolPanelCoordinator::closeAllTools()
{

    for (size_t i = 0; i < options.tpOpen.size(); i++)
        if (i < expList.size()) {
            expList.at(i)->set_expanded (false);
        }
}

void ToolPanelCoordinator::openAllTools()
{

    for (size_t i = 0; i < options.tpOpen.size(); i++)
        if (i < expList.size()) {
            expList.at(i)->set_expanded (true);
        }
}

void ToolPanelCoordinator::updateToolState()
{

    for (size_t i = 0; i < options.tpOpen.size(); i++)
        if (i < expList.size()) {
            expList.at(i)->set_expanded (options.tpOpen.at(i));
        }

    if(options.tpOpen.size() > expList.size()) {
        size_t sizeWavelet = options.tpOpen.size() - expList.size();
        std::vector<int> temp;

        for (size_t i = 0; i < sizeWavelet; i++) {
            temp.push_back(options.tpOpen.at(i + expList.size()));
        }

        wavelet->updateToolState(temp);
        wavelet->setExpanded(true);
    }
}

void ToolPanelCoordinator::readOptions ()
{

    crop->readOptions ();
}

void ToolPanelCoordinator::writeOptions ()
{

    crop->writeOptions ();
    options.tpOpen.clear ();

    for (size_t i = 0; i < expList.size(); i++) {
        options.tpOpen.push_back (expList.at(i)->get_expanded ());
    }

    wavelet->writeOptions(options.tpOpen);
    retinex->writeOptions(options.tpOpen);
}


void ToolPanelCoordinator::cropSelectionReady ()
{

    toolBar->setTool (TMHand);

    if (!ipc) {
        return;
    }
}

void ToolPanelCoordinator::rotateSelectionReady (double rotate_deg, Thumbnail* thm)
{

    toolBar->setTool (TMHand);

    if (!ipc) {
        return;
    }

    if (rotate_deg != 0.0) {
        rotate->straighten (rotate_deg);
    }
}

void ToolPanelCoordinator::spotWBselected (int x, int y, Thumbnail* thm)
{

    if (!ipc) {
        return;
    }

//    toolBar->setTool (TOOL_HAND);
    int rect = whitebalance->getSize ();
    int ww = ipc->getFullWidth();
    int hh = ipc->getFullHeight();

    if (x - rect > 0 && y - rect > 0 && x + rect < ww && y + rect < hh) {
        double temp;
        double green;
        ipc->getSpotWB (x, y, rect, temp, green);
        whitebalance->setWB (temp, green);
    }
}




void ToolPanelCoordinator::autoCropRequested ()
{

    if (!ipc) {
        return;
    }

    int x1, y1, x2, y2, w, h;
    ipc->getAutoCrop (crop->getRatio(), x1, y1, w, h);
    x2 = x1 + w - 1;
    y2 = y1 + h - 1;
    crop->cropInit (x1, y1, w, h);
    crop->cropResized (x1, y1, x2, y2);
    crop->cropManipReady ();
}

rtengine::RawImage* ToolPanelCoordinator::getDF()
{
    if (!ipc) {
        return NULL;
    }

    const rtengine::ImageMetaData *imd = ipc->getInitialImage()->getMetaData();

    if(imd) {
        int iso = imd->getISOSpeed();
        double shutter = imd->getShutterSpeed();
        std::string maker( imd->getMake()  );
        std::string model( imd->getModel() );
        time_t timestamp = imd->getDateTimeAsTS();

        return rtengine::dfm.searchDarkFrame( maker, model, iso, shutter, timestamp);
    }

    return NULL;
}

rtengine::RawImage* ToolPanelCoordinator::getFF()
{
    if (!ipc) {
        return NULL;
    }

    const rtengine::ImageMetaData *imd = ipc->getInitialImage()->getMetaData();

    if(imd) {
        // int iso = imd->getISOSpeed();              temporarilly removed because unused
        // double shutter = imd->getShutterSpeed();   temporarilly removed because unused
        double aperture = imd->getFNumber();
        double focallength = imd->getFocalLen();
        std::string maker( imd->getMake()  );
        std::string model( imd->getModel() );
        std::string lens(  imd->getLens()  );
        time_t timestamp = imd->getDateTimeAsTS();

        return rtengine::ffm.searchFlatField( maker, model, lens, focallength, aperture, timestamp);
    }

    return NULL;
}

Glib::ustring ToolPanelCoordinator::GetCurrentImageFilePath()
{
    if (!ipc) {
        return "";
    }

    return ipc->getInitialImage()->getFileName();
}

void ToolPanelCoordinator::straightenRequested ()
{

    if (!ipc) {
        return;
    }

    toolBar->setTool (TMStraighten);
}

double ToolPanelCoordinator::autoDistorRequested ()
{
    if (!ipc) {
        return 0.0;
    }

    return rtengine::ImProcFunctions::getAutoDistor (ipc->getInitialImage()->getFileName(), 400);
}

void ToolPanelCoordinator::spotWBRequested (int size)
{

    if (!ipc) {
        return;
    }

    toolBar->setTool (TMSpotWB);
}

void ToolPanelCoordinator::cropSelectRequested ()
{

    if (!ipc) {
        return;
    }

    toolBar->setTool (TMCropSelect);
}

void ToolPanelCoordinator::saveInputICCReference (Glib::ustring fname, bool apply_wb)
{

    if (ipc) {
        ipc->saveInputICCReference (fname, apply_wb);
    }
}

int ToolPanelCoordinator::getSpotWBRectSize ()
{

    return whitebalance->getSize ();
}

void ToolPanelCoordinator::updateCurveBackgroundHistogram (LUTu & histToneCurve, LUTu & histLCurve, LUTu & histCCurve, /*LUTu & histCLurve, LUTu & histLLCurve,*/ LUTu & histLCAM, LUTu & histCCAM, LUTu & histRed, LUTu & histGreen, LUTu & histBlue, LUTu & histLuma, LUTu & histLRETI)
{
    colorappearance->updateCurveBackgroundHistogram (histToneCurve, histLCurve, histCCurve, /*histCLurve, histLLCurve,*/ histLCAM,  histCCAM, histRed, histGreen, histBlue, histLuma, histLRETI);
    toneCurve->updateCurveBackgroundHistogram (histToneCurve, histLCurve, histCCurve,/* histCLurve, histLLCurve,*/ histLCAM,  histCCAM, histRed, histGreen, histBlue, histLuma, histLRETI);
    lcurve->updateCurveBackgroundHistogram (histToneCurve, histLCurve, histCCurve, /*histCLurve, histLLCurve,*/ histLCAM, histCCAM, histRed, histGreen, histBlue, histLuma, histLRETI);
    rgbcurves->updateCurveBackgroundHistogram(histToneCurve, histLCurve, histCCurve,/* histCLurve, histLLCurve, */histLCAM, histCCAM, histRed, histGreen, histBlue, histLuma, histLRETI);
    retinex->updateCurveBackgroundHistogram(histToneCurve, histLCurve, histCCurve,/* histCLurve, histLLCurve, */histLCAM, histCCAM, histRed, histGreen, histBlue, histLuma, histLRETI);

}

void ToolPanelCoordinator::foldAllButOne (Gtk::Box* parent, FoldableToolPanel* openedSection)
{

    for (auto toolPanel : toolPanels) {
        if (toolPanel->getParent() != NULL) {
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

bool ToolPanelCoordinator::handleShortcutKey (GdkEventKey* event)
{

    //bool ctrl = event->state & GDK_CONTROL_MASK;  temporarily removed because unused
    //bool shift = event->state & GDK_SHIFT_MASK;   temporarily removed because unused
    bool alt = event->state & GDK_MOD1_MASK;

    if (alt) {
        switch(event->keyval) {
        case GDK_KEY_e:
            toolPanelNotebook->set_current_page (toolPanelNotebook->page_num(*exposurePanelSW));
            return true;

        case GDK_KEY_d:
            toolPanelNotebook->set_current_page (toolPanelNotebook->page_num(*detailsPanelSW));
            return true;

        case GDK_KEY_c:
            toolPanelNotebook->set_current_page (toolPanelNotebook->page_num(*colorPanelSW));
            return true;

        case GDK_KEY_t:
            toolPanelNotebook->set_current_page (toolPanelNotebook->page_num(*transformPanelSW));
            return true;

        case GDK_KEY_r:
            toolPanelNotebook->set_current_page (toolPanelNotebook->page_num(*rawPanelSW));
            return true;

        case GDK_KEY_w:
            toolPanelNotebook->set_current_page (toolPanelNotebook->page_num(*waveletPanelSW));
            return true;

        case GDK_KEY_m:
            if (metadataPanel) {
                toolPanelNotebook->set_current_page (toolPanelNotebook->page_num(*metadataPanel));
                return true;
            }
        }
    }

    return false;
}

void ToolPanelCoordinator::updateVScrollbars (bool hide)
{
    GThreadLock lock; // All GUI acces from idle_add callbacks or separate thread HAVE to be protected
    Gtk::PolicyType policy = hide ? Gtk::POLICY_NEVER : Gtk::POLICY_AUTOMATIC;
    exposurePanelSW->set_policy     (Gtk::POLICY_AUTOMATIC, policy);
    detailsPanelSW->set_policy      (Gtk::POLICY_AUTOMATIC, policy);
    colorPanelSW->set_policy        (Gtk::POLICY_AUTOMATIC, policy);
    transformPanelSW->set_policy    (Gtk::POLICY_AUTOMATIC, policy);
    rawPanelSW->set_policy          (Gtk::POLICY_AUTOMATIC, policy);
    waveletPanelSW->set_policy      (Gtk::POLICY_AUTOMATIC, policy);
}

void ToolPanelCoordinator::updateTabsHeader (bool useIcons)
{
    GThreadLock lock; // All GUI acces from idle_add callbacks or separate thread HAVE to be protected
    TOITypes type = useIcons ? TOI_ICON : TOI_TEXT;

    toiE->switchTo(type);
    toiD->switchTo(type);
    toiC->switchTo(type);
    toiT->switchTo(type);
    toiR->switchTo(type);

    if (toiM) {
        toiM->switchTo(type);
    }
}

void ToolPanelCoordinator::updateTPVScrollbar (bool hide)
{
    updateVScrollbars (hide);
}

void ToolPanelCoordinator::updateTabsUsesIcons (bool useIcons)
{
    updateTabsHeader (useIcons);
}

void ToolPanelCoordinator::toolSelected (ToolMode tool)
{
    GThreadLock lock; // All GUI acces from idle_add callbacks or separate thread HAVE to be protected

    switch (tool) {
    case TMCropSelect:
        crop->setExpanded(true);
        toolPanelNotebook->set_current_page(toolPanelNotebook->page_num(*transformPanelSW));
        break;

    case TMSpotWB:
        whitebalance->setExpanded(true);
        toolPanelNotebook->set_current_page(toolPanelNotebook->page_num(*colorPanelSW));
        break;

    case TMStraighten:
        lensgeom->setExpanded(true);
        rotate->setExpanded(true);
        toolPanelNotebook->set_current_page(toolPanelNotebook->page_num(*transformPanelSW));
        break;

    default:
        break;
    }
}

void ToolPanelCoordinator::editModeSwitchedOff ()
{
    if (editDataProvider) {
        editDataProvider->switchOffEditMode();
    }
}

void ToolPanelCoordinator::dirSelected (const Glib::ustring& dirname, const Glib::ustring& openfile)
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
