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
#include "batchtoolpanelcoord.h"
#include "options.h"
#include "filepanel.h"
#include "procparamchangers.h"
#include "addsetids.h"

using namespace rtengine::procparams;

BatchToolPanelCoordinator::BatchToolPanelCoordinator (FilePanel* parent) : ToolPanelCoordinator(), parent(parent), somethingChanged(false) {

	blockedUpdate = false;
	// remove exif panel and iptc panel
	std::vector<ToolPanel*>::iterator epi = std::find (toolPanels.begin(), toolPanels.end(), exifpanel);
	if (epi!=toolPanels.end())
		toolPanels.erase (epi);
	std::vector<ToolPanel*>::iterator ipi = std::find (toolPanels.begin(), toolPanels.end(), iptcpanel);
	if (ipi!=toolPanels.end())
		toolPanels.erase (ipi);
	toolPanelNotebook->remove_page (*metadataPanel);
	metadataPanel = 0;
	toiM = 0;

    for (size_t i=0; i<toolPanels.size(); i++)
        toolPanels[i]->setBatchMode (true);
}

void BatchToolPanelCoordinator::selectionChanged (const std::vector<Thumbnail*>& selected) {

    if (selected!=this->selected) {
        closeSession ();
        this->selected = selected;
        selFileNames.clear ();
        for (size_t i=0; i<selected.size(); i++)
            selFileNames.push_back (selected[i]->getFileName ());
        initSession ();
    }
}

void BatchToolPanelCoordinator::closeSession (bool save) {

    pparamsEdited.set (false);

    for (size_t i=0; i<selected.size(); i++)
        selected[i]->removeThumbnailListener (this);

    if (somethingChanged && save) {

        // read new values from the gui
        for (size_t i=0; i<toolPanels.size(); i++)
            toolPanels[i]->write (&pparams, &pparamsEdited);

        // combine with initial parameters and set
        ProcParams newParams;
        for (size_t i=0; i<selected.size(); i++) {
            newParams = initialPP[i];
            pparamsEdited.combine (newParams, pparams, selected.size()==1);

            // trim new adjuster's values to the adjuster's limits
            for (unsigned int j=0; j<toolPanels.size(); j++)
                toolPanels[j]->trimValues (&newParams);

            selected[i]->setProcParams (newParams, NULL, BATCHEDITOR, true);
        }
    }
    for (size_t i=0; i<paramcListeners.size(); i++)
        paramcListeners[i]->clearParamChanges ();
}

void BatchToolPanelCoordinator::initSession () {

    somethingChanged = false;

    initialPP.resize (selected.size());
    for (size_t i=0; i<selected.size(); i++) {
        initialPP[i] = selected[i]->getProcParams ();
        selected[i]->applyAutoExp (initialPP[i]);
        selected[i]->addThumbnailListener (this);
    }

    // compare all the ProcParams and describe which parameters has different (i.e. inconsistent) values in pparamsEdited
    pparamsEdited.initFrom (initialPP);

    crop->setDimensions (100000, 100000);

/*    if (!selected.empty()) {
        pparams = selected[0]->getProcParams ();
        for (int i=0; i<toolPanels.size(); i++) {
            toolPanels[i]->setDefaults (&pparams, &pparamsEdited);
            toolPanels[i]->read (&pparams, &pparamsEdited);
        }
        for (int i=0; i<paramcListeners.size(); i++)
            paramcListeners[i]->procParamsChanged (&pparams, rtengine::EvPhotoLoaded, "batch processing", &pparamsEdited);
    }
*/

    if (!selected.empty()) {

		// The first selected image (in the thumbnail list, not the click list) is used to populate the EditorPanel and set the default values
		pparams = selected[0]->getProcParams ();

		coarse->initBatchBehavior ();

		if (selected.size()==1) {

			for (size_t i=0; i<toolPanels.size(); i++)
				toolPanels.at(i)->setMultiImage(false);

			toneCurve->setAdjusterBehavior (false, false, false, false, false, false, false, false);
			lcurve->setAdjusterBehavior (false, false, false);
			whitebalance->setAdjusterBehavior (false, false, false);
			vibrance->setAdjusterBehavior (false, false);
			vignetting->setAdjusterBehavior (false, false, false, false);
			colorappearance->setAdjusterBehavior (false, false, false, false, false, false, false, false, false, false, false, false, false);
			rotate->setAdjusterBehavior (false);
			distortion->setAdjusterBehavior (false);
			perspective->setAdjusterBehavior (false);
			gradient->setAdjusterBehavior (false, false, false, false);
			pcvignette->setAdjusterBehavior (false, false, false);
			cacorrection->setAdjusterBehavior (false);
			sharpening->setAdjusterBehavior (false);
			sharpenEdge->setAdjusterBehavior (false, false);
			sharpenMicro->setAdjusterBehavior (false, false);
			icm->setAdjusterBehavior (false, false);

			chmixer->setAdjusterBehavior (false);
			blackwhite->setAdjusterBehavior (false,false);
			colortoning->setAdjusterBehavior (false, false, false, false, false);
			filmSimulation->setAdjusterBehavior(false);

			shadowshighlights->setAdjusterBehavior (false, false, false);
			dirpyrequalizer->setAdjusterBehavior (false, false, false);
			wavelet->setAdjusterBehavior (false, false, false,false,false,false,false,false,false,false,false,false,false,false);
			dirpyrdenoise->setAdjusterBehavior (false, false,false,false,false,false, false);
			bayerpreprocess->setAdjusterBehavior (false, false);
			rawcacorrection->setAdjusterBehavior (false);
			flatfield->setAdjusterBehavior(false);
			rawexposure->setAdjusterBehavior (false, false);
			bayerrawexposure->setAdjusterBehavior (false);
			xtransrawexposure->setAdjusterBehavior (false);
		}
		else {

			for (size_t i=0; i<toolPanels.size(); i++)
				toolPanels.at(i)->setMultiImage(true);

			toneCurve->setAdjusterBehavior (options.baBehav[ADDSET_TC_EXPCOMP], options.baBehav[ADDSET_TC_HLCOMPAMOUNT],options.baBehav[ADDSET_TC_HLCOMPTHRESH], options.baBehav[ADDSET_TC_BRIGHTNESS], options.baBehav[ADDSET_TC_BLACKLEVEL],options.baBehav[ADDSET_TC_SHCOMP], options.baBehav[ADDSET_TC_CONTRAST], options.baBehav[ADDSET_TC_SATURATION]);
			lcurve->setAdjusterBehavior (options.baBehav[ADDSET_LC_BRIGHTNESS], options.baBehav[ADDSET_LC_CONTRAST], options.baBehav[ADDSET_LC_CHROMATICITY]);
			whitebalance->setAdjusterBehavior (options.baBehav[ADDSET_WB_TEMPERATURE], options.baBehav[ADDSET_WB_GREEN], options.baBehav[ADDSET_WB_EQUAL]);
			vibrance->setAdjusterBehavior (options.baBehav[ADDSET_VIBRANCE_PASTELS], options.baBehav[ADDSET_VIBRANCE_SATURATED]);
			vignetting->setAdjusterBehavior (options.baBehav[ADDSET_VIGN_AMOUNT], options.baBehav[ADDSET_VIGN_RADIUS], options.baBehav[ADDSET_VIGN_STRENGTH], options.baBehav[ADDSET_VIGN_CENTER]);
			colorappearance->setAdjusterBehavior (options.baBehav[ADDSET_CAT_DEGREE], options.baBehav[ADDSET_CAT_ADAPTSCENE], options.baBehav[ADDSET_CAT_ADAPTVIEWING],options.baBehav[ADDSET_CAT_BADPIX], options.baBehav[ADDSET_CAT_LIGHT], options.baBehav[ADDSET_CAT_CHROMA],options.baBehav[ADDSET_CAT_CONTRAST],options.baBehav[ADDSET_CAT_RSTPRO],options.baBehav[ADDSET_CAT_BRIGHT],options.baBehav[ADDSET_CAT_CONTRAST_Q],options.baBehav[ADDSET_CAT_CHROMA_S],options.baBehav[ADDSET_CAT_CHROMA_M],options.baBehav[ADDSET_CAT_HUE]);
			rotate->setAdjusterBehavior (options.baBehav[ADDSET_ROTATE_DEGREE]);
			distortion->setAdjusterBehavior (options.baBehav[ADDSET_DIST_AMOUNT]);
			perspective->setAdjusterBehavior (options.baBehav[ADDSET_PERSPECTIVE]);
			gradient->setAdjusterBehavior (options.baBehav[ADDSET_GRADIENT_DEGREE], options.baBehav[ADDSET_GRADIENT_FEATHER], options.baBehav[ADDSET_GRADIENT_STRENGTH], options.baBehav[ADDSET_GRADIENT_CENTER]);
			pcvignette->setAdjusterBehavior (options.baBehav[ADDSET_PCVIGNETTE_STRENGTH], options.baBehav[ADDSET_PCVIGNETTE_FEATHER], options.baBehav[ADDSET_PCVIGNETTE_ROUNDNESS]);
			cacorrection->setAdjusterBehavior (options.baBehav[ADDSET_CA]);
			sharpening->setAdjusterBehavior (options.baBehav[ADDSET_SHARP_AMOUNT]);
			sharpenEdge->setAdjusterBehavior (options.baBehav[ADDSET_SHARPENEDGE_AMOUNT],options.baBehav[ADDSET_SHARPENEDGE_PASS]);
			sharpenMicro->setAdjusterBehavior (options.baBehav[ADDSET_SHARPENMICRO_AMOUNT],options.baBehav[ADDSET_SHARPENMICRO_UNIFORMITY]);
			icm->setAdjusterBehavior (options.baBehav[ADDSET_FREE_OUPUT_GAMMA],options.baBehav[ADDSET_FREE_OUTPUT_SLOPE]);
//			colortoning->setAdjusterBehavior (options.baBehav[ADDSET_COLORTONING_SPLIT], options.baBehav[ADDSET_COLORTONING_SATTHRESHOLD], options.baBehav[ADDSET_COLORTONING_SATOPACITY], options.baBehav[ADDSET_COLORTONING_STRPROTECT], options.baBehav[ADDSET_COLORTONING_BALANCE]);
			colortoning->setAdjusterBehavior (options.baBehav[ADDSET_COLORTONING_SPLIT], options.baBehav[ADDSET_COLORTONING_SATTHRESHOLD], options.baBehav[ADDSET_COLORTONING_SATOPACITY], options.baBehav[ADDSET_COLORTONING_STRENGTH], options.baBehav[ADDSET_COLORTONING_BALANCE]);
			filmSimulation->setAdjusterBehavior(options.baBehav[ADDSET_FILMSIMULATION_STRENGTH]);

			chmixer->setAdjusterBehavior (options.baBehav[ADDSET_CHMIXER] );
			blackwhite->setAdjusterBehavior (options.baBehav[ADDSET_BLACKWHITE_HUES],options.baBehav[ADDSET_BLACKWHITE_GAMMA]);
			shadowshighlights->setAdjusterBehavior (options.baBehav[ADDSET_SH_HIGHLIGHTS], options.baBehav[ADDSET_SH_SHADOWS], options.baBehav[ADDSET_SH_LOCALCONTRAST]);
			dirpyrequalizer->setAdjusterBehavior (options.baBehav[ADDSET_DIRPYREQ], options.baBehav[ADDSET_DIRPYREQ_THRESHOLD], options.baBehav[ADDSET_DIRPYREQ_SKINPROTECT]);
			wavelet->setAdjusterBehavior (options.baBehav[ADDSET_WA], options.baBehav[ADDSET_WA_THRESHOLD], options.baBehav[ADDSET_WA_THRESHOLD2],options.baBehav[ADDSET_WA_THRES],options.baBehav[ADDSET_WA_CHRO],options.baBehav[ADDSET_WA_CHROMA],options.baBehav[ADDSET_WA_UNIF],options.baBehav[ADDSET_WA_SKINPROTECT],options.baBehav[ADDSET_WA_RESCHRO],options.baBehav[ADDSET_WA_RESCON],options.baBehav[ADDSET_WA_RESCONH],options.baBehav[ADDSET_WA_THRR],options.baBehav[ADDSET_WA_THRRH],options.baBehav[ADDSET_WA_SKYPROTECT] );
			dirpyrdenoise->setAdjusterBehavior (options.baBehav[ADDSET_DIRPYRDN_LUMA],options.baBehav[ADDSET_DIRPYRDN_LUMDET],options.baBehav[ADDSET_DIRPYRDN_CHROMA],options.baBehav[ADDSET_DIRPYRDN_CHROMARED],options.baBehav[ADDSET_DIRPYRDN_CHROMABLUE], options.baBehav[ADDSET_DIRPYRDN_GAMMA], options.baBehav[ADDSET_DIRPYRDN_PASSES]);
			bayerpreprocess->setAdjusterBehavior (options.baBehav[ADDSET_PREPROCESS_LINEDENOISE], options.baBehav[ADDSET_PREPROCESS_GREENEQUIL]);
			rawcacorrection->setAdjusterBehavior (options.baBehav[ADDSET_RAWCACORR]);
			flatfield->setAdjusterBehavior(options.baBehav[ADDSET_RAWFFCLIPCONTROL]);
			rawexposure->setAdjusterBehavior (options.baBehav[ADDSET_RAWEXPOS_LINEAR], options.baBehav[ADDSET_RAWEXPOS_PRESER]);
			bayerrawexposure->setAdjusterBehavior (options.baBehav[ADDSET_RAWEXPOS_BLACKS]);
			xtransrawexposure->setAdjusterBehavior (options.baBehav[ADDSET_RAWEXPOS_BLACKS]);

			if (options.baBehav[ADDSET_TC_EXPCOMP])  pparams.toneCurve.expcomp = 0;
			if (options.baBehav[ADDSET_TC_HLCOMPAMOUNT])  pparams.toneCurve.hlcompr = 0;
			if (options.baBehav[ADDSET_TC_HLCOMPTHRESH])  pparams.toneCurve.hlcomprthresh = 0;
			if (options.baBehav[ADDSET_TC_BRIGHTNESS])  pparams.toneCurve.brightness = 0;
			if (options.baBehav[ADDSET_TC_BLACKLEVEL])  pparams.toneCurve.black = 0;
			if (options.baBehav[ADDSET_TC_SHCOMP])  pparams.toneCurve.shcompr = 0;
			if (options.baBehav[ADDSET_TC_CONTRAST])  pparams.toneCurve.contrast = 0;

			if (options.baBehav[ADDSET_SH_HIGHLIGHTS])  pparams.sh.highlights = 0;
			if (options.baBehav[ADDSET_SH_SHADOWS])  pparams.sh.shadows = 0;
			if (options.baBehav[ADDSET_SH_LOCALCONTRAST])  pparams.sh.localcontrast = 0;

			if (options.baBehav[ADDSET_LC_BRIGHTNESS])  pparams.labCurve.brightness = 0;
			if (options.baBehav[ADDSET_LC_CONTRAST])  pparams.labCurve.contrast = 0;
			if (options.baBehav[ADDSET_LC_CHROMATICITY])  pparams.labCurve.chromaticity = 0;

			if (options.baBehav[ADDSET_SHARP_AMOUNT])  pparams.sharpening.amount = 0;
			if (options.baBehav[ADDSET_SHARPENEDGE_AMOUNT])  pparams.sharpenEdge.amount = 0;
			if (options.baBehav[ADDSET_SHARPENMICRO_AMOUNT])  pparams.sharpenMicro.amount = 0;
			if (options.baBehav[ADDSET_SHARPENEDGE_PASS])  pparams.sharpenEdge.passes = 0;
			if (options.baBehav[ADDSET_SHARPENMICRO_UNIFORMITY])  pparams.sharpenMicro.uniformity = 0;

			if (options.baBehav[ADDSET_CHMIXER]) for (int i=0; i<3; i++) pparams.chmixer.red[i] = pparams.chmixer.green[i] = pparams.chmixer.blue[i] = 0;
			if (options.baBehav[ADDSET_BLACKWHITE_HUES]) pparams.blackwhite.mixerRed=pparams.blackwhite.mixerOrange=pparams.blackwhite.mixerYellow=
														 pparams.blackwhite.mixerGreen=pparams.blackwhite.mixerCyan=pparams.blackwhite.mixerBlue=
														 pparams.blackwhite.mixerMagenta=pparams.blackwhite.mixerPurple=0;
			if (options.baBehav[ADDSET_BLACKWHITE_GAMMA]) pparams.blackwhite.gammaRed=pparams.blackwhite.gammaGreen=pparams.blackwhite.gammaBlue = 0;

			//if (options.baBehav[ADDSET_LD_EDGETOLERANCE])  pparams.lumaDenoise.edgetolerance = 0;

			if (options.baBehav[ADDSET_WB_TEMPERATURE])  pparams.wb.temperature = 0;
			if (options.baBehav[ADDSET_WB_GREEN])  pparams.wb.green = 0;
			if (options.baBehav[ADDSET_WB_EQUAL])  pparams.wb.equal = 0;

			if (options.baBehav[ADDSET_VIBRANCE_PASTELS])  pparams.vibrance.pastels = 0;
			if (options.baBehav[ADDSET_VIBRANCE_SATURATED])  pparams.vibrance.saturated = 0;

			if (options.baBehav[ADDSET_CAT_DEGREE])  pparams.colorappearance.degree = 0;
			if (options.baBehav[ADDSET_CAT_ADAPTSCENE])  pparams.colorappearance.adapscen = 0;
			if (options.baBehav[ADDSET_CAT_ADAPTVIEWING])  pparams.colorappearance.adaplum = 0;
			if (options.baBehav[ADDSET_CAT_BADPIX])  pparams.colorappearance.badpixsl = 0;
			if (options.baBehav[ADDSET_CAT_LIGHT])  pparams.colorappearance.jlight = 0;
			if (options.baBehav[ADDSET_CAT_BRIGHT])  pparams.colorappearance.qbright = 0;
			if (options.baBehav[ADDSET_CAT_CHROMA])  pparams.colorappearance.chroma = 0;
			if (options.baBehav[ADDSET_CAT_CHROMA_S])  pparams.colorappearance.schroma = 0;
			if (options.baBehav[ADDSET_CAT_CHROMA_M])  pparams.colorappearance.mchroma = 0;
			if (options.baBehav[ADDSET_CAT_RSTPRO])  pparams.colorappearance.rstprotection = 0;
			if (options.baBehav[ADDSET_CAT_CONTRAST])  pparams.colorappearance.contrast = 0;
			if (options.baBehav[ADDSET_CAT_CONTRAST_Q])  pparams.colorappearance.qcontrast = 0;
			if (options.baBehav[ADDSET_CAT_HUE])  pparams.colorappearance.colorh = 0;

			if (options.baBehav[ADDSET_FREE_OUPUT_GAMMA])  pparams.icm.gampos = 0;
			if (options.baBehav[ADDSET_FREE_OUTPUT_SLOPE])  pparams.icm.slpos = 0;

			//if (options.baBehav[ADDSET_CBOOST_AMOUNT])  pparams.colorBoost.amount = 0;

			//if (options.baBehav[ADDSET_CS_BLUEYELLOW])  pparams.colorShift.a = 0;
			//if (options.baBehav[ADDSET_CS_GREENMAGENTA])  pparams.colorShift.b = 0;
			if (options.baBehav[ADDSET_COLORTONING_SPLIT]) pparams.colorToning.redlow  = pparams.colorToning.greenlow  = pparams.colorToning.bluelow =
			                                               pparams.colorToning.redmed  = pparams.colorToning.greenmed  = pparams.colorToning.bluemed =
			                                               pparams.colorToning.redhigh = pparams.colorToning.greenhigh = pparams.colorToning.bluehigh =
			                                               pparams.colorToning.satlow = pparams.colorToning.sathigh = 0;
			if (options.baBehav[ADDSET_COLORTONING_SATTHRESHOLD]) pparams.colorToning.satProtectionThreshold = 0;
			if (options.baBehav[ADDSET_COLORTONING_SATOPACITY]) pparams.colorToning.saturatedOpacity = 0;
			if (options.baBehav[ADDSET_COLORTONING_BALANCE]) pparams.colorToning.balance = 0;
			if (options.baBehav[ADDSET_COLORTONING_STRENGTH]) pparams.colorToning.strength = 0;

			if (options.baBehav[ADDSET_FILMSIMULATION_STRENGTH]) pparams.filmSimulation.strength = 0;

			if (options.baBehav[ADDSET_ROTATE_DEGREE])  pparams.rotate.degree = 0;
			if (options.baBehav[ADDSET_DIST_AMOUNT])  pparams.distortion.amount = 0;
			if (options.baBehav[ADDSET_PERSPECTIVE])  pparams.perspective.horizontal = pparams.perspective.vertical = 0;
			if (options.baBehav[ADDSET_GRADIENT_DEGREE])  pparams.gradient.degree = 0;
			if (options.baBehav[ADDSET_GRADIENT_FEATHER])  pparams.gradient.feather = 0;
			if (options.baBehav[ADDSET_GRADIENT_STRENGTH])  pparams.gradient.strength = 0;
			if (options.baBehav[ADDSET_GRADIENT_CENTER])  pparams.gradient.centerX = 0;
			if (options.baBehav[ADDSET_GRADIENT_CENTER])  pparams.gradient.centerY = 0;
			if (options.baBehav[ADDSET_PCVIGNETTE_STRENGTH])  pparams.pcvignette.strength = 0;
			if (options.baBehav[ADDSET_PCVIGNETTE_FEATHER])  pparams.pcvignette.feather = 0;
			if (options.baBehav[ADDSET_PCVIGNETTE_ROUNDNESS])  pparams.pcvignette.roundness = 0;
			if (options.baBehav[ADDSET_CA])  pparams.cacorrection.red = 0;
			if (options.baBehav[ADDSET_CA])  pparams.cacorrection.blue = 0;
			if (options.baBehav[ADDSET_VIGN_AMOUNT])  pparams.vignetting.amount = 0;
			if (options.baBehav[ADDSET_VIGN_RADIUS])  pparams.vignetting.radius = 0;
			if (options.baBehav[ADDSET_VIGN_STRENGTH])  pparams.vignetting.strength = 0;
			if (options.baBehav[ADDSET_VIGN_CENTER])  pparams.vignetting.centerX = 0;
			if (options.baBehav[ADDSET_VIGN_CENTER])  pparams.vignetting.centerY = 0;

			if (options.baBehav[ADDSET_DIRPYREQ]) for (int i=0; i<6; i++) pparams.dirpyrequalizer.mult[i] = 0;
			if (options.baBehav[ADDSET_DIRPYREQ_THRESHOLD]) pparams.dirpyrequalizer.threshold = 0;
			if (options.baBehav[ADDSET_DIRPYREQ_SKINPROTECT]) pparams.dirpyrequalizer.skinprotect = 0;

			if (options.baBehav[ADDSET_WA]) for (int i=0; i<8; i++) pparams.wavelet.c[i] = 0;
			if (options.baBehav[ADDSET_WA_THRESHOLD]) pparams.wavelet.threshold = 0;
			if (options.baBehav[ADDSET_WA_THRESHOLD2]) pparams.wavelet.threshold2 = 0;
			if (options.baBehav[ADDSET_WA_SKINPROTECT]) pparams.wavelet.skinprotect = 0;
			if (options.baBehav[ADDSET_WA_CHRO]) pparams.wavelet.chro = 0;
			if (options.baBehav[ADDSET_WA_CHROMA]) pparams.wavelet.chroma = 0;
			if (options.baBehav[ADDSET_WA_UNIF]) pparams.wavelet.unif = 0;
			if (options.baBehav[ADDSET_WA_THRES]) pparams.wavelet.thres = 0;
			if (options.baBehav[ADDSET_WA_RESCON]) pparams.wavelet.rescon = 0;
			if (options.baBehav[ADDSET_WA_RESCONH]) pparams.wavelet.resconH = 0;
			if (options.baBehav[ADDSET_WA_RESCHRO]) pparams.wavelet.reschro = 0;
			if (options.baBehav[ADDSET_WA_THRR]) pparams.wavelet.thr = 0;
			if (options.baBehav[ADDSET_WA_THRRH]) pparams.wavelet.thrH = 0;
			if (options.baBehav[ADDSET_WA_SKYPROTECT]) pparams.wavelet.sky = 0;
			
			if (options.baBehav[ADDSET_DIRPYRDN_LUMA]) pparams.dirpyrDenoise.luma = 0;

			if (options.baBehav[ADDSET_DIRPYRDN_CHROMA]) pparams.dirpyrDenoise.chroma = 0;
			if (options.baBehav[ADDSET_DIRPYRDN_CHROMARED]) pparams.dirpyrDenoise.redchro = 0;
			if (options.baBehav[ADDSET_DIRPYRDN_CHROMABLUE]) pparams.dirpyrDenoise.bluechro = 0;
//			pparams.dirpyrDenoise.Ldetail = pparams.dirpyrDenoise.luma = pparams.dirpyrDenoise.chroma = 0;
			if (options.baBehav[ADDSET_DIRPYRDN_GAMMA])  pparams.dirpyrDenoise.gamma = 0;

			if (options.baBehav[ADDSET_RAWCACORR])  pparams.raw.cablue = pparams.raw.cared = 0;
			if (options.baBehav[ADDSET_RAWEXPOS_LINEAR])  pparams.raw.expos = 0;
			if (options.baBehav[ADDSET_RAWEXPOS_PRESER])  pparams.raw.preser = 0;
			if (options.baBehav[ADDSET_RAWEXPOS_BLACKS]) {
				pparams.raw.bayersensor.black0 = pparams.raw.bayersensor.black1 = pparams.raw.bayersensor.black2 = pparams.raw.bayersensor.black3 = 0;
				pparams.raw.xtranssensor.blackred = pparams.raw.xtranssensor.blackgreen = pparams.raw.xtranssensor.blackblue = 0;
			}
			if (options.baBehav[ADDSET_RAWFFCLIPCONTROL])  pparams.raw.ff_clipControl = 0;

			if (options.baBehav[ADDSET_PREPROCESS_GREENEQUIL])  pparams.raw.bayersensor.greenthresh = 0;
			if (options.baBehav[ADDSET_PREPROCESS_LINEDENOISE])  pparams.raw.bayersensor.linenoise = 0;
		}

		for (size_t i=0; i<toolPanels.size(); i++) {
			toolPanels[i]->setDefaults (&pparams, &pparamsEdited);
			toolPanels[i]->read (&pparams, &pparamsEdited);
			// TODO: autoOpenCurve has been disabled because initSession is called on each parameter change from the editor panel,
			// if the thumbnail remains selected in the DirectoryBrowser (i.e. always, unless the user think about deselecting it)
			//toolPanels[i]->autoOpenCurve();
		}
		for (size_t i=0; i<paramcListeners.size(); i++)
			// send this initial state to the History
			paramcListeners[i]->procParamsChanged (&pparams, rtengine::EvPhotoLoaded, M("BATCH_PROCESSING"), &pparamsEdited);
	}
}

void BatchToolPanelCoordinator::panelChanged (rtengine::ProcEvent event, const Glib::ustring& descr) {

    if (selected.empty())
        return;

    somethingChanged = true;

    pparamsEdited.set (false);
    // read new values from the gui
    for (size_t i=0; i<toolPanels.size(); i++)
        toolPanels[i]->write (&pparams, &pparamsEdited);

    // If only a single item is selected, we emulate the behaviour of the editor tool panel coordinator,
    // otherwise we adjust the inital parameters on a per-image basis.
    if (selected.size()==1) {
        // Compensate rotation on flip
        if (event==rtengine::EvCTHFlip || event==rtengine::EvCTVFlip) {
            if (fabs(pparams.rotate.degree)>0.001) {
                  pparams.rotate.degree *= -1;
                rotate->read (&pparams);
            }
        }

        int w, h;
        selected[0]->getFinalSize (selected[0]->getProcParams (), w, h);
        crop->setDimensions(w, h);

        // Some transformations change the crop and resize parameter for convenience.
        if (event==rtengine::EvCTHFlip) {
            crop->hFlipCrop ();
            crop->write (&pparams, &pparamsEdited);
        }
        else if (event==rtengine::EvCTVFlip) {
            crop->vFlipCrop ();
            crop->write (&pparams, &pparamsEdited);
        }
        else if (event==rtengine::EvCTRotate) {
            crop->rotateCrop (pparams.coarse.rotate, pparams.coarse.hflip, pparams.coarse.vflip);
            crop->write (&pparams, &pparamsEdited);
            resize->update (pparams.crop.enabled, pparams.crop.w, pparams.crop.h, w, h);
            resize->write (&pparams, &pparamsEdited);
        }
        else if (event==rtengine::EvCrop) {
            resize->update (pparams.crop.enabled, pparams.crop.w, pparams.crop.h);
            resize->write (&pparams, &pparamsEdited);
        }
    }
    else {
        // Compensate rotation on flip
        if (event==rtengine::EvCTHFlip || event==rtengine::EvCTVFlip) {
            for (size_t i=0; i<selected.size(); i++) {
                if (fabs(initialPP[i].rotate.degree) > 0.001) {
                    initialPP[i].rotate.degree *= -1.0;

                    pparamsEdited.rotate.degree = false;
                }
            }
        }

        // some transformations make the crop change for convenience
        if (event==rtengine::EvCTHFlip) {
            for (size_t i=0; i<selected.size(); i++) {
                int w, h;
                selected[i]->getFinalSize (selected[i]->getProcParams (), w, h);

                rtengine::procparams::CropParams& crop = initialPP[i].crop;
                crop.x = w - crop.x - crop.w;

                pparamsEdited.crop.x = false;
            }
        }
        else if (event==rtengine::EvCTVFlip) {
            for (size_t i=0; i<selected.size(); i++) {
                int w, h;
                selected[i]->getFinalSize (selected[i]->getProcParams (), w, h);

                rtengine::procparams::CropParams& crop = initialPP[i].crop;
                crop.y = h - crop.y - crop.h;

                pparamsEdited.crop.y = false;
            }
        }
        else if (event==rtengine::EvCTRotate) {
            int newDeg = pparams.coarse.rotate;
            for (size_t i=0; i<selected.size(); i++) {
                int w, h;
                selected[i]->getFinalSize (selected[i]->getProcParams (), w, h);

                int oldDeg = initialPP[i].coarse.rotate;

                rtengine::procparams::CropParams& crop = initialPP[i].crop;
                int rotation = (360 + newDeg - oldDeg) % 360;
                ProcParams pptemp = selected[i]->getProcParams(); // Get actual procparams
                if((pptemp.coarse.hflip != pptemp.coarse.vflip) && ((rotation%180)==90))
                    rotation = (rotation + 180)%360;


                switch (rotation) {
                case 90:
                    std::swap(crop.x, crop.y);
                    std::swap(crop.w, crop.h);

                    crop.x = h - crop.x - crop.w;
                    break;
                case 270:
                    std::swap(crop.x, crop.y);
                    std::swap(crop.w, crop.h);

                    crop.y = w - crop.y - crop.h;
                    break;
                case 180:
                    crop.x = w - crop.x - crop.w;
                    crop.y = h - crop.y - crop.h;
                    break;
                }

                initialPP[i].coarse.rotate = newDeg;

            }
            pparamsEdited.crop.x = false;
            pparamsEdited.crop.y = false;
            pparamsEdited.crop.w = false;
            pparamsEdited.crop.h = false;
            pparamsEdited.coarse.rotate = false;
        }
    }

    if (event==rtengine::EvAutoExp || event==rtengine::EvClip)
        for (size_t i=0; i<selected.size(); i++) {
            initialPP[i].toneCurve.autoexp = pparams.toneCurve.autoexp;
            initialPP[i].toneCurve.clip = pparams.toneCurve.clip;

            // at this stage, we don't know if HL Reconstruction will be enabled or not (depending on the raw histogram),
            // so we're forcing its paramseditd value to false (i.e. mixed state)
            pparamsEdited.toneCurve.hrenabled = false;

            selected[i]->applyAutoExp (initialPP[i]);
        }

    if (event==rtengine::EvAutoDIST) {
        for (size_t i=0; i<selected.size(); i++) {
            initialPP[i].distortion.amount = pparams.distortion.amount;
        }
    }

    // combine with initial parameters and set
    ProcParams newParams;
    for (size_t i=0; i<selected.size(); i++) {
        newParams = initialPP[i];
        // If only one file is selected, slider's addMode has been set to false, and hence the behave
        // like in SET mode like in an editor ; that's why we force the combination to the SET mode too
        pparamsEdited.combine (newParams, pparams, selected.size()==1);

		// trim new adjuster's values to the adjuster's limits
		for (unsigned int j=0; j<toolPanels.size(); j++)
			toolPanels[j]->trimValues (&newParams);

        selected[i]->setProcParams (newParams, NULL, BATCHEDITOR, false);
    }

    for (size_t i=0; i<paramcListeners.size(); i++)
        paramcListeners[i]->procParamsChanged (&pparams, event, descr, &pparamsEdited);
}

void BatchToolPanelCoordinator::getAutoWB (double& temp, double& green, double equal) {

    if (!selected.empty())
        selected[0]->getAutoWB (temp, green, equal);
}

void BatchToolPanelCoordinator::getCamWB (double& temp, double& green) {

    if (!selected.empty())
        selected[0]->getCamWB (temp, green);
}

void BatchToolPanelCoordinator::optionsChanged () {

    closeSession ();
    initSession ();
}

void BatchToolPanelCoordinator::procParamsChanged (Thumbnail* thm, int whoChangedIt) {

    if (whoChangedIt!=BATCHEDITOR && !blockedUpdate) {
        closeSession (false);
        initSession ();
    }
}

void BatchToolPanelCoordinator::beginBatchPParamsChange(int numberOfEntries) {

    blockedUpdate = true;
    if (numberOfEntries > 50) // Arbitrary amount
        parent->set_sensitive(false);
}

// The end of a batch pparams change triggers a close/initsession
void BatchToolPanelCoordinator::endBatchPParamsChange() {
    //printf("BatchToolPanelCoordinator::endBatchPParamsChange  /  Nouvelle session!\n");
    closeSession (false);
    initSession ();
    blockedUpdate = false;
    parent->set_sensitive(true);
}

/*
 * WARNING: profileChange is actually called by the History only.
 *          Using a Profile panel in the batch tool panel editor is actually
 *          not supported by BatchToolPanelCoordinator::profileChange!
 */
void BatchToolPanelCoordinator::profileChange  (const rtengine::procparams::PartialProfile* nparams, rtengine::ProcEvent event, const Glib::ustring& descr, const ParamsEdited* paramsEdited) {

    if (event==rtengine::EvProfileChanged) {
    	// a profile has been selected in a hypothetical Profile panel
    	// -> ACTUALLY NOT SUPPORTED
    	return;
    }

    pparams = *(nparams->pparams);
    if (paramsEdited)
        pparamsEdited = *paramsEdited;


    for (size_t i=0; i<toolPanels.size(); i++)
        // writing the values to the GUI
        toolPanels[i]->read (&pparams, &pparamsEdited);
        // I guess we don't want to automatically unfold curve editors here...

    somethingChanged = true;

    // read new values from the gui
    for (size_t i=0; i<toolPanels.size(); i++)
        toolPanels[i]->write (&pparams, &pparamsEdited);

    // combine with initial parameters of each image and set
    ProcParams newParams;
    for (size_t i=0; i<selected.size(); i++) {
        newParams = initialPP[i];
        pparamsEdited.combine (newParams, pparams, selected.size()==1);
        selected[i]->setProcParams (newParams, NULL, BATCHEDITOR, false);
    }

    for (size_t i=0; i<paramcListeners.size(); i++)
        paramcListeners[i]->procParamsChanged (&pparams, event, descr, &pparamsEdited);
}

void BatchToolPanelCoordinator::cropSelectionReady () {

  toolBar->setTool (TMHand);
}

CropGUIListener* BatchToolPanelCoordinator::startCropEditing (Thumbnail* thm) {

    if (thm) {
        int w, h;
        thm->getFinalSize (thm->getProcParams (), w, h);
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
        for (size_t i=0; i<selected.size(); i++)
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

