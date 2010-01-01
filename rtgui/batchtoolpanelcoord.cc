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
#include <multilangmgr.h>
#include <batchtoolpanelcoord.h>
#include <options.h>
#include <filepanel.h>
#include <procparamchangers.h>

using namespace rtengine::procparams;

BatchToolPanelCoordinator::BatchToolPanelCoordinator (FilePanel* parent) : ToolPanelCoordinator(), parent(parent) {

	// remove exif panel and iptc panel
	std::vector<ToolPanel*>::iterator epi = std::find (toolPanels.begin(), toolPanels.end(), exifpanel);
	if (epi!=toolPanels.end())
		toolPanels.erase (epi);
	std::vector<ToolPanel*>::iterator ipi = std::find (toolPanels.begin(), toolPanels.end(), iptcpanel);
	if (ipi!=toolPanels.end())
		toolPanels.erase (ipi);
	toolPanelNotebook->remove_page (*metadataPanel);

    for (int i=0; i<toolPanels.size(); i++)
        toolPanels[i]->setBatchMode (true);
}

void BatchToolPanelCoordinator::selectionChanged (const std::vector<Thumbnail*>& selected) {

    if (selected!=this->selected) {
        closeSession ();
        this->selected = selected;
        selFileNames.clear ();
        for (int i=0; i<selected.size(); i++)
            selFileNames.push_back (selected[i]->getFileName ());
        initSession ();
    }
}

void BatchToolPanelCoordinator::closeSession (bool save) {

    pparamsEdited.set (false);        

    for (int i=0; i<selected.size(); i++)
        selected[i]->removeThumbnailListener (this);

    if (somethingChanged && save) {
        
        // read new values from the gui
        for (int i=0; i<toolPanels.size(); i++)
            toolPanels[i]->write (&pparams, &pparamsEdited);

        // combine with initial parameters and set
        ProcParams newParams;
        for (int i=0; i<selected.size(); i++) {
            newParams = initialPP[i];
            pparamsEdited.combine (newParams, pparams);
            selected[i]->setProcParams (newParams, BATCHEDITOR, true);
        }
    }
    for (int i=0; i<paramcListeners.size(); i++)
        paramcListeners[i]->clearParamChanges ();
}

void BatchToolPanelCoordinator::initSession () {

    somethingChanged = false;

    initialPP.resize (selected.size());
    for (int i=0; i<selected.size(); i++) {
        initialPP[i] = selected[i]->getProcParams ();
        selected[i]->applyAutoExp (initialPP[i]);
        selected[i]->addThumbnailListener (this);
    }

    pparamsEdited.initFrom (initialPP);

/*    curve->setAdjusterBehavior (false, false, false, false);
    whitebalance->setAdjusterBehavior (false, false);
    vignetting->setAdjusterBehavior (false);
    rotate->setAdjusterBehavior (false);
    distortion->setAdjusterBehavior (false);
    cacorrection->setAdjusterBehavior (false, false);
    colorshift->setAdjusterBehavior (false, false);
    colorboost->setAdjusterBehavior (false);
    lumadenoise->setAdjusterBehavior (false);
    sharpening->setAdjusterBehavior (false);
    shadowshighlights->setAdjusterBehavior (false, false, false);
*/
    crop->setDimensions (100000, 100000);

/*    if (selected.size()>0) {
        pparams = selected[0]->getProcParams ();
        for (int i=0; i<toolPanels.size(); i++) {
            toolPanels[i]->setDefaults (&pparams, &pparamsEdited);
            toolPanels[i]->read (&pparams, &pparamsEdited);
        }
        for (int i=0; i<paramcListeners.size(); i++)
            paramcListeners[i]->procParamsChanged (&pparams, rtengine::EvPhotoLoaded, "batch processing", &pparamsEdited);
    }
*/

    if (selected.size()>0) {

        pparams = selected[0]->getProcParams ();
        coarse->initBatchBehavior ();

        curve->setAdjusterBehavior (options.baBehav[0], options.baBehav[1], options.baBehav[2], options.baBehav[3]);
        whitebalance->setAdjusterBehavior (options.baBehav[12], options.baBehav[13]);
        vignetting->setAdjusterBehavior (options.baBehav[21]);
        rotate->setAdjusterBehavior (options.baBehav[17]);
        distortion->setAdjusterBehavior (options.baBehav[18]);
        cacorrection->setAdjusterBehavior (options.baBehav[19], options.baBehav[20]);
        colorshift->setAdjusterBehavior (options.baBehav[15], options.baBehav[16]);
        colorboost->setAdjusterBehavior (options.baBehav[14]);
        lumadenoise->setAdjusterBehavior (options.baBehav[11]);
        sharpening->setAdjusterBehavior (options.baBehav[10]);
        shadowshighlights->setAdjusterBehavior (options.baBehav[4], options.baBehav[5], options.baBehav[6]);
        
        if (options.baBehav[0])  pparams.toneCurve.expcomp = 0;
        if (options.baBehav[1])  pparams.toneCurve.brightness = 0;
        if (options.baBehav[2])  pparams.toneCurve.black = 0;
        if (options.baBehav[3])  pparams.toneCurve.contrast = 0;

        if (options.baBehav[4])  pparams.sh.highlights = 0;
        if (options.baBehav[5])  pparams.sh.shadows = 0;
        if (options.baBehav[6])  pparams.sh.localcontrast = 0;
        
        if (options.baBehav[10])  pparams.sharpening.amount = 0;
        if (options.baBehav[11])  pparams.lumaDenoise.edgetolerance = 0;

        if (options.baBehav[12])  pparams.wb.temperature = 0;
        if (options.baBehav[13])  pparams.wb.green = 0;

        if (options.baBehav[14])  pparams.colorBoost.amount = 0;

        if (options.baBehav[15])  pparams.colorShift.a = 0;
        if (options.baBehav[16])  pparams.colorShift.b = 0;

        if (options.baBehav[17])  pparams.rotate.degree = 0;
        if (options.baBehav[18])  pparams.distortion.amount = 0;
        if (options.baBehav[19])  pparams.cacorrection.red = 0;
        if (options.baBehav[20])  pparams.cacorrection.blue = 0;
        if (options.baBehav[21])  pparams.vignetting.amount = 0;

        for (int i=0; i<toolPanels.size(); i++) {
            toolPanels[i]->setDefaults (&pparams, &pparamsEdited);
            toolPanels[i]->read (&pparams, &pparamsEdited);
        }
        for (int i=0; i<paramcListeners.size(); i++)
            paramcListeners[i]->procParamsChanged (&pparams, rtengine::EvPhotoLoaded, "batch processing", &pparamsEdited);
    }
}

void BatchToolPanelCoordinator::panelChanged (rtengine::ProcEvent event, const Glib::ustring& descr) {

    if (selected.size()==0)
        return;

    somethingChanged = true;

    pparamsEdited.set (false);        
    // read new values from the gui
    for (int i=0; i<toolPanels.size(); i++)
        toolPanels[i]->write (&pparams, &pparamsEdited);

    if (event==rtengine::EvAutoExp || event==rtengine::EvClip) 
        for (int i=0; i<selected.size(); i++) {
            initialPP[i].toneCurve.autoexp = pparams.toneCurve.autoexp;
            initialPP[i].toneCurve.clip = pparams.toneCurve.clip;
            selected[i]->applyAutoExp (initialPP[i]);
        }

    // combine with initial parameters and set
    ProcParams newParams;
    for (int i=0; i<selected.size(); i++) {
        newParams = initialPP[i];
        pparamsEdited.combine (newParams, pparams);
        selected[i]->setProcParams (newParams, BATCHEDITOR, false);
    }

    for (int i=0; i<paramcListeners.size(); i++)
        paramcListeners[i]->procParamsChanged (&pparams, event, descr, &pparamsEdited);
}

void BatchToolPanelCoordinator::getAutoWB (double& temp, double& green) {

    if (selected.size()>0)
        selected[0]->getAutoWB (temp, green);       
}

void BatchToolPanelCoordinator::getCamWB (double& temp, double& green) {
    
    if (selected.size()>0)
        selected[0]->getCamWB (temp, green);
}    

void BatchToolPanelCoordinator::optionsChanged () {

    closeSession ();
    initSession ();
}

void BatchToolPanelCoordinator::procParamsChanged (Thumbnail* thm, int whoChangedIt) {

    if (whoChangedIt!=BATCHEDITOR) {
        closeSession (false);
        initSession ();
    }
}

void BatchToolPanelCoordinator::profileChange  (const ProcParams *nparams, rtengine::ProcEvent event, const Glib::ustring& descr, const ParamsEdited* paramsEdited) {

    pparams = *nparams;
    if (paramsEdited)
        pparamsEdited = *paramsEdited;

    for (int i=0; i<toolPanels.size(); i++) 
        toolPanels[i]->read (&pparams, &pparamsEdited);

    somethingChanged = true;

    // read new values from the gui
    for (int i=0; i<toolPanels.size(); i++)
        toolPanels[i]->write (&pparams, &pparamsEdited);

    // combine with initial parameters and set
    ProcParams newParams;
    for (int i=0; i<selected.size(); i++) {
        newParams = initialPP[i];
        pparamsEdited.combine (newParams, pparams);
        selected[i]->setProcParams (newParams, BATCHEDITOR, false);
    }

    for (int i=0; i<paramcListeners.size(); i++)
        paramcListeners[i]->procParamsChanged (&pparams, event, descr, &pparamsEdited);
}

void BatchToolPanelCoordinator::cropSelectionReady () {

  toolBar->setTool (TMHand);
}

CropGUIListener* BatchToolPanelCoordinator::startCropEditing (Thumbnail* thm) {
    
    if (thm) {
        int w, h;
        thm->getFinalSize (thm->getProcParams (), w, h);
        printf ("final=%d %d\n", w, h);
        crop->setDimensions (w, h);
    }
    return crop;
}

void BatchToolPanelCoordinator::rotateSelectionReady (double rotate_deg, Thumbnail* thm) {

  toolBar->setTool (TMHand);
  if (rotate_deg!=0.0)
      rotate->straighten (rotate_deg);
}

void BatchToolPanelCoordinator::spotWBselected (int x, int y, Thumbnail* thm) {

//    toolBar->setTool (TOOL_HAND);
    if (x>0 && y>0 && thm) {
        for (int i=0; i<selected.size(); i++)
            if (selected[i]==thm) {
                double temp;
                double green;
                thm->getSpotWB (x, y, whitebalance->getSize(), temp, green);
                double otemp = initialPP[i].wb.temperature;
                double ogreen = initialPP[i].wb.green;
                if (options.baBehav[12])
                    temp = temp - otemp;
                if (options.baBehav[13])
                    green = green - ogreen;
                whitebalance->setWB (temp, green);
            }
    }
}

