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
#include "paramsedited.h"
#include <cstring>
#include "options.h"
#include "addsetids.h"

ParamsEdited::ParamsEdited () {

    set (false);
}

void ParamsEdited::set (bool v) {

	general.rank         = v;
	general.colorlabel   = v;
	general.intrash      = v;
	toneCurve.curve      = v;
	toneCurve.curve2     = v;
	toneCurve.curveMode  = v;
	toneCurve.curveMode2 = v;
	toneCurve.brightness = v;
	toneCurve.black      = v;
	toneCurve.contrast   = v;
	toneCurve.saturation = v;
	toneCurve.shcompr    = v;
	toneCurve.hlcompr    = v;
	toneCurve.hlcomprthresh = v;
	toneCurve.autoexp    = v;
	toneCurve.clip       = v;
	toneCurve.expcomp    = v;
	toneCurve.hrenabled   = v;
	toneCurve.method    = v;
	labCurve.lcurve      = v;
	labCurve.acurve      = v;
	labCurve.bcurve      = v;
	labCurve.cccurve     = v;
	labCurve.chcurve     = v;
	labCurve.lhcurve     = v;
	labCurve.hhcurve     = v;
	labCurve.lccurve    = v;
	labCurve.clcurve    = v;
	labCurve.brightness  = v;
	labCurve.contrast    = v;
	labCurve.chromaticity    = v;
	labCurve.avoidcolorshift = v;
	labCurve.rstprotection   = v;
	labCurve.lcredsk         = v;
	rgbCurves.lumamode       = v;
	rgbCurves.rcurve         = v;
	rgbCurves.gcurve         = v;
	rgbCurves.bcurve         = v;
	colorToning.enabled      = v;
	colorToning.autosat      = v;
	colorToning.opacityCurve = v;
	colorToning.colorCurve   = v;
	colorToning.satprotectionthreshold = v;
	colorToning.saturatedopacity       = v;	
	colorToning.strength               = v;
	colorToning.shadowsColSat          = v;
	colorToning.hlColSat   = v;
	colorToning.balance    = v;
	colorToning.clcurve    = v;
	colorToning.method     = v;
	colorToning.twocolor   = v;
	colorToning.cl2curve   = v;
	colorToning.redlow     = v;
	colorToning.greenlow   = v;
	colorToning.bluelow    = v;
	colorToning.satlow     = v;
	colorToning.sathigh    = v;
	colorToning.redmed     = v;
	colorToning.greenmed   = v;
	colorToning.bluemed    = v;
	colorToning.redhigh    = v;
	colorToning.greenhigh  = v;
	colorToning.bluehigh   = v;
	colorToning.lumamode   = v;

	sharpening.enabled            = v;
	sharpening.radius             = v;
	sharpening.amount             = v;
	sharpening.threshold          = v;
	sharpening.edgesonly          = v;
	sharpening.edges_radius       = v;
	sharpening.edges_tolerance    = v;
	sharpening.halocontrol        = v;
	sharpening.halocontrol_amount = v;
	sharpening.method         = v;
	sharpening.deconvamount   = v;
	sharpening.deconvradius   = v;
	sharpening.deconviter     = v;
	sharpening.deconvdamping  = v;
	sharpenEdge.enabled       = v;
	sharpenEdge.passes        = v;
	sharpenEdge.amount        = v;
	sharpenEdge.threechannels = v;
	sharpenMicro.enabled      = v;
	sharpenMicro.matrix       = v;
	sharpenMicro.amount       = v;
	sharpenMicro.uniformity   = v;
	vibrance.enabled          = v;
	vibrance.pastels          = v;
	vibrance.saturated        = v;
	vibrance.psthreshold      = v;
	vibrance.protectskins     = v;
	vibrance.avoidcolorshift  = v;
	vibrance.pastsattog       = v;
	vibrance.skintonescurve   = v;
	colorappearance.enabled    = v;
	colorappearance.degree     = v;
	colorappearance.autodegree = v;
	colorappearance.surround     = v;
	colorappearance.adapscen    = v;
	colorappearance.autoadapscen = v;
	colorappearance.adaplum    = v;
	colorappearance.badpixsl    = v;
	colorappearance.wbmodel    = v;
	colorappearance.algo    = v;
	
	colorappearance.jlight     = v;
	colorappearance.qbright     = v;
	colorappearance.chroma     = v;
	colorappearance.schroma     = v;
	colorappearance.mchroma     = v;
	colorappearance.contrast     = v;
	colorappearance.qcontrast     = v;
	colorappearance.colorh     = v;
	colorappearance.rstprotection     = v;
	colorappearance.surrsource = v;
	colorappearance.gamut = v;
//	colorappearance.badpix = v;
	colorappearance.datacie = v;
	colorappearance.tonecie = v;
//	colorappearance.sharpcie = v;
	colorappearance.curve      = v;
	colorappearance.curve2     = v;
	colorappearance.curve3     = v;
    colorappearance.curveMode  = v;
    colorappearance.curveMode2 = v;
    colorappearance.curveMode3 = v;
	
	//colorBoost.amount         = v;
	//colorBoost.avoidclip      = v;
	//colorBoost.enable_saturationlimiter = v;
	//colorBoost.saturationlimit = v;
	wb.method                  = v;
	wb.green                   = v;
	wb.temperature             = v;
	wb.equal                   = v;
	//colorShift.a               = v;
	//colorShift.b               = v;
	//lumaDenoise.enabled        = v;
	//lumaDenoise.radius         = v;
	//lumaDenoise.edgetolerance  = v;
	//colorDenoise.enabled       = v;
	//colorDenoise.amount        = v;
	defringe.enabled           = v;
	defringe.radius            = v;
	defringe.threshold         = v;
	defringe.huecurve          = v;
	impulseDenoise.enabled     = v;
	impulseDenoise.thresh      = v;
	dirpyrDenoise.enabled      = v;
	dirpyrDenoise.enhance      = v;
//	dirpyrDenoise.perform      = v;
	dirpyrDenoise.lcurve      = v;
	dirpyrDenoise.cccurve      = v;
	dirpyrDenoise.median      = v;
	dirpyrDenoise.autochroma      = v;
	dirpyrDenoise.luma         = v;
	dirpyrDenoise.Ldetail      = v;
	dirpyrDenoise.chroma       = v;
	dirpyrDenoise.redchro      = v;
	dirpyrDenoise.bluechro     = v;
	dirpyrDenoise.gamma        = v;
	dirpyrDenoise.passes        = v;
	dirpyrDenoise.dmethod      = v;
	dirpyrDenoise.Lmethod      = v;
	dirpyrDenoise.Cmethod      = v;
	dirpyrDenoise.C2method      = v;
	dirpyrDenoise.smethod      = v;
	dirpyrDenoise.medmethod      = v;
	dirpyrDenoise.methodmed      = v;
	dirpyrDenoise.rgbmethod      = v;
	epd.enabled                = v;
	epd.strength            = v;
	epd.edgeStopping        = v;
	epd.scale               = v;
	epd.reweightingIterates = v;
	sh.enabled       = v;
	sh.hq            = v;
	sh.highlights    = v;
	sh.htonalwidth   = v;
	sh.shadows       = v;
	sh.stonalwidth   = v;
	sh.localcontrast = v;
	sh.radius        = v;
	crop.enabled = v;
	crop.x       = v;
	crop.y       = v;
	crop.w       = v;
	crop.h       = v;
	crop.fixratio = v;
	crop.ratio   = v;
	crop.orientation = v;
	crop.guide   = v;
	coarse.rotate = v;
	coarse.hflip = v;
	coarse.vflip = v;
	commonTrans.autofill = v;
	rotate.degree = v;
	distortion.amount = v;
	lensProf.lcpFile = v;
	lensProf.useDist = v;
	lensProf.useVign = v;
	lensProf.useCA = v;
	perspective.horizontal = v;
	perspective.vertical = v;
	gradient.enabled = v;
	gradient.degree = v;
	gradient.feather = v;
	gradient.strength = v;
	gradient.centerX = v;
	gradient.centerY = v;
	pcvignette.enabled = v;
	pcvignette.strength = v;
	pcvignette.feather = v;
	pcvignette.roundness = v;
	cacorrection.red = v;
	cacorrection.blue = v;
	vignetting.amount = v;
	vignetting.radius = v;
	vignetting.strength = v;
	vignetting.centerX = v;
	vignetting.centerY = v;
	chmixer.red[0] = v;
	chmixer.red[1] = v;
	chmixer.red[2] = v;
	chmixer.green[0] = v;
	chmixer.green[1] = v;
	chmixer.green[2] = v;
	chmixer.blue[0] = v;
	chmixer.blue[1] = v;
	chmixer.blue[2] = v;
	blackwhite.enabled   = v;
	blackwhite.enabledcc   = v;
	blackwhite.mixerRed   = v;
	blackwhite.mixerOrange   = v;
	blackwhite.mixerYellow   = v;
	blackwhite.mixerGreen   = v;
	blackwhite.mixerCyan   = v;
	blackwhite.mixerBlue   = v;
	blackwhite.mixerMagenta   = v;
	blackwhite.mixerPurple   = v;
	blackwhite.gammaRed   = v;
	blackwhite.gammaGreen   = v;
	blackwhite.gammaBlue   = v;
	blackwhite.filter   = v;
	blackwhite.setting   = v;
	blackwhite.method   = v;
	blackwhite.luminanceCurve = v;
	blackwhite.beforeCurve      = v;
	blackwhite.beforeCurveMode  = v;
	blackwhite.afterCurve      = v;
	blackwhite.afterCurveMode  = v;
	blackwhite.autoc    = v;
	blackwhite.algo    = v;

	
	resize.scale     = v;
	resize.appliesTo = v;
	resize.method    = v;
	resize.dataspec  = v;
	resize.width     = v;
	resize.height    = v;
	resize.enabled   = v;
	icm.input        = v;
	icm.toneCurve = v;
	icm.blendCMSMatrix = v;
	icm.dcpIlluminant = v;
	icm.working      = v;
	icm.output       = v;
	icm.gamma		= v;
	icm.freegamma		= v;
	icm.gampos		= v;
	icm.slpos		= v;
	raw.bayersensor.method = v;
	raw.bayersensor.ccSteps = v;
	raw.bayersensor.exBlack0 = v;
	raw.bayersensor.exBlack1 = v;
	raw.bayersensor.exBlack2 = v;
	raw.bayersensor.exBlack3 = v;
	raw.bayersensor.exTwoGreen=v;
	raw.bayersensor.dcbIterations = v;
	raw.bayersensor.dcbEnhance = v;
	//raw.bayersensor.allEnhance = v;
	raw.bayersensor.lmmseIterations = v;
	raw.bayersensor.greenEq = v;
	raw.bayersensor.linenoise = v;
	raw.xtranssensor.method = v;
	raw.xtranssensor.ccSteps = v;
	raw.xtranssensor.exBlackRed= v;
	raw.xtranssensor.exBlackGreen = v;
	raw.xtranssensor.exBlackBlue = v;
	raw.caCorrection = v;
	raw.caBlue  = v;
	raw.caRed   = v;
	raw.hotPixelFilter = v;
	raw.deadPixelFilter = v;
	raw.hotDeadPixelThresh = v;
	raw.darkFrame = v;
	raw.dfAuto = v;
	raw.ff_file = v;
	raw.ff_AutoSelect = v;
	raw.ff_BlurRadius = v;
	raw.ff_BlurType = v;
	raw.ff_AutoClipControl = v;
	raw.ff_clipControl = v;
	raw.exPos = v;
	raw.exPreser = v;
	wavelet.enabled = v;
	wavelet.median = v;
	wavelet.avoid = v;
	wavelet.Lmethod = v;
	wavelet.CLmethod = v;
	wavelet.Tilesmethod = v;
	wavelet.CHmethod = v;
	wavelet.HSmethod = v;
	wavelet.Dirmethod = v;
	wavelet.tiles = v;
	wavelet.rescon = v;
	wavelet.resconH = v;
	wavelet.reschro = v;
	wavelet.sup = v;
	wavelet.sky = v;
	wavelet.thres = v;
	wavelet.threshold = v;
	wavelet.threshold2 = v;
	wavelet.chroma = v;
	wavelet.chro = v;
	wavelet.unif = v;
	wavelet.thr = v;
	wavelet.thrH = v;
	wavelet.skinprotect = v;
	wavelet.hueskin = v;
	wavelet.hueskin2 = v;
	wavelet.hllev = v;
	wavelet.bllev = v;
	wavelet.clvcurve = v;
	wavelet.opacityCurveRG   = v;
	wavelet.opacityCurveBY   = v;
	wavelet.pastlev = v;
	wavelet.satlev = v;
	
	for(int i = 0; i < 9; i++) {
		wavelet.c[i] = v;
	}

	dirpyrequalizer.enabled = v;
	dirpyrequalizer.gamutlab = v;
	for(int i = 0; i < 6; i++) {
		dirpyrequalizer.mult[i] = v;
	}
	dirpyrequalizer.threshold = v;
	dirpyrequalizer.skinprotect = v;
	dirpyrequalizer.hueskin = v;
	//dirpyrequalizer.algo = v;
	hsvequalizer.hcurve = v;
	hsvequalizer.scurve = v;
	hsvequalizer.vcurve = v;
    filmSimulation.enabled = v;
    filmSimulation.clutFilename = v;
    filmSimulation.strength = v;

	exif = v;
	iptc = v;
}

using namespace rtengine;
using namespace rtengine::procparams;

void ParamsEdited::initFrom (const std::vector<rtengine::procparams::ProcParams>& src) {

    set (true);
    if (src.empty())
        return;

    const ProcParams& p = src[0];
    for (size_t i=1; i<src.size(); i++) {
        const ProcParams& other = src[i];
        toneCurve.curve = toneCurve.curve && p.toneCurve.curve == other.toneCurve.curve;
        toneCurve.curve2 = toneCurve.curve2 && p.toneCurve.curve2 == other.toneCurve.curve2;
        toneCurve.curveMode = toneCurve.curveMode && p.toneCurve.curveMode == other.toneCurve.curveMode;
        toneCurve.curveMode2 = toneCurve.curveMode2 && p.toneCurve.curveMode2 == other.toneCurve.curveMode2;
        toneCurve.brightness = toneCurve.brightness && p.toneCurve.brightness == other.toneCurve.brightness;
        toneCurve.black = toneCurve.black && p.toneCurve.black == other.toneCurve.black;
        toneCurve.contrast = toneCurve.contrast && p.toneCurve.contrast == other.toneCurve.contrast;
        toneCurve.saturation = toneCurve.saturation && p.toneCurve.saturation == other.toneCurve.saturation;
        toneCurve.shcompr = toneCurve.shcompr && p.toneCurve.shcompr == other.toneCurve.shcompr;
        toneCurve.hlcompr = toneCurve.hlcompr && p.toneCurve.hlcompr == other.toneCurve.hlcompr;
        toneCurve.hlcomprthresh = toneCurve.hlcomprthresh && p.toneCurve.hlcomprthresh == other.toneCurve.hlcomprthresh;
        toneCurve.autoexp = toneCurve.autoexp && p.toneCurve.autoexp == other.toneCurve.autoexp;
        toneCurve.clip = toneCurve.clip && p.toneCurve.clip == other.toneCurve.clip;
        toneCurve.expcomp = toneCurve.expcomp && p.toneCurve.expcomp == other.toneCurve.expcomp;
        toneCurve.hrenabled = toneCurve.hrenabled && p.toneCurve.hrenabled == other.toneCurve.hrenabled;
        toneCurve.method = toneCurve.method && p.toneCurve.method == other.toneCurve.method;
        labCurve.lcurve = labCurve.lcurve && p.labCurve.lcurve == other.labCurve.lcurve;
        labCurve.acurve = labCurve.acurve && p.labCurve.acurve == other.labCurve.acurve;
        labCurve.bcurve = labCurve.bcurve && p.labCurve.bcurve == other.labCurve.bcurve;
        labCurve.cccurve = labCurve.cccurve && p.labCurve.cccurve == other.labCurve.cccurve;
        labCurve.chcurve = labCurve.chcurve && p.labCurve.chcurve == other.labCurve.chcurve;
        labCurve.lhcurve = labCurve.lhcurve && p.labCurve.lhcurve == other.labCurve.lhcurve;
        labCurve.hhcurve = labCurve.hhcurve && p.labCurve.hhcurve == other.labCurve.hhcurve;
        labCurve.lccurve = labCurve.lccurve && p.labCurve.lccurve == other.labCurve.lccurve;
        labCurve.clcurve = labCurve.clcurve && p.labCurve.clcurve == other.labCurve.clcurve;
        labCurve.brightness = labCurve.brightness && p.labCurve.brightness == other.labCurve.brightness;
        labCurve.contrast = labCurve.contrast && p.labCurve.contrast == other.labCurve.contrast;
        labCurve.chromaticity = labCurve.chromaticity && p.labCurve.chromaticity == other.labCurve.chromaticity;
        labCurve.avoidcolorshift = labCurve.avoidcolorshift && p.labCurve.avoidcolorshift == other.labCurve.avoidcolorshift;
        labCurve.rstprotection = labCurve.rstprotection && p.labCurve.rstprotection == other.labCurve.rstprotection;
        labCurve.lcredsk = labCurve.lcredsk && p.labCurve.lcredsk == other.labCurve.lcredsk;
        rgbCurves.lumamode = rgbCurves.lumamode && p.rgbCurves.lumamode == other.rgbCurves.lumamode;
        rgbCurves.rcurve = rgbCurves.rcurve && p.rgbCurves.rcurve == other.rgbCurves.rcurve;
        rgbCurves.gcurve = rgbCurves.gcurve && p.rgbCurves.gcurve == other.rgbCurves.gcurve;
        rgbCurves.bcurve = rgbCurves.bcurve && p.rgbCurves.bcurve == other.rgbCurves.bcurve;
        colorToning.enabled = colorToning.enabled && p.colorToning.enabled == other.colorToning.enabled;
        colorToning.twocolor = colorToning.twocolor && p.colorToning.twocolor == other.colorToning.twocolor;
        colorToning.opacityCurve = colorToning.opacityCurve && p.colorToning.opacityCurve == other.colorToning.opacityCurve;
        colorToning.colorCurve = colorToning.colorCurve && p.colorToning.colorCurve == other.colorToning.colorCurve;
        colorToning.autosat = colorToning.autosat && p.colorToning.autosat == other.colorToning.autosat;
        colorToning.satprotectionthreshold = colorToning.satprotectionthreshold && p.colorToning.satProtectionThreshold == other.colorToning.satProtectionThreshold;
        colorToning.saturatedopacity = colorToning.saturatedopacity && p.colorToning.saturatedOpacity == other.colorToning.saturatedOpacity;
        colorToning.strength = colorToning.strength && p.colorToning.strength == other.colorToning.strength;
        colorToning.shadowsColSat = colorToning.shadowsColSat && p.colorToning.shadowsColSat == other.colorToning.shadowsColSat;
        colorToning.hlColSat = colorToning.hlColSat && p.colorToning.hlColSat == other.colorToning.hlColSat;
        colorToning.balance = colorToning.balance && p.colorToning.balance == other.colorToning.balance;
        colorToning.clcurve = colorToning.clcurve && p.colorToning.clcurve == other.colorToning.clcurve;
        colorToning.cl2curve = colorToning.cl2curve && p.colorToning.cl2curve == other.colorToning.cl2curve;
        colorToning.method = colorToning.method && p.colorToning.method == other.colorToning.method;
        colorToning.redlow = colorToning.redlow && p.colorToning.redlow == other.colorToning.redlow;
        colorToning.greenlow = colorToning.greenlow && p.colorToning.greenlow == other.colorToning.greenlow;
        colorToning.bluelow = colorToning.bluelow && p.colorToning.bluelow == other.colorToning.bluelow;
        colorToning.satlow = colorToning.satlow && p.colorToning.satlow == other.colorToning.satlow;
        colorToning.sathigh = colorToning.sathigh && p.colorToning.sathigh == other.colorToning.sathigh;
        colorToning.redmed = colorToning.redmed && p.colorToning.redmed == other.colorToning.redmed;
        colorToning.greenmed = colorToning.greenmed && p.colorToning.greenmed == other.colorToning.greenmed;
        colorToning.bluemed = colorToning.bluemed && p.colorToning.bluemed == other.colorToning.bluemed;
        colorToning.redhigh = colorToning.redhigh && p.colorToning.redhigh == other.colorToning.redhigh;
        colorToning.greenhigh = colorToning.greenhigh && p.colorToning.greenhigh == other.colorToning.greenhigh;
        colorToning.bluehigh = colorToning.bluehigh && p.colorToning.bluehigh == other.colorToning.bluehigh;
        colorToning.lumamode = colorToning.lumamode && p.colorToning.lumamode == other.colorToning.lumamode;
        sharpenEdge.enabled = sharpenEdge.enabled && p.sharpenEdge.enabled == other.sharpenEdge.enabled;
        sharpenEdge.passes = sharpenEdge.passes && p.sharpenEdge.passes == other.sharpenEdge.passes;
        sharpenEdge.amount = sharpenEdge.amount && p.sharpenEdge.amount == other.sharpenEdge.amount;
        sharpenEdge.threechannels = sharpenEdge.threechannels && p.sharpenEdge.threechannels == other.sharpenEdge.threechannels;
        sharpenMicro.enabled = sharpenMicro.enabled && p.sharpenMicro.enabled == other.sharpenMicro.enabled;
        sharpenMicro.matrix = sharpenMicro.matrix && p.sharpenMicro.matrix == other.sharpenMicro.matrix;
        sharpenMicro.amount = sharpenMicro.amount && p.sharpenMicro.amount == other.sharpenMicro.amount;
        sharpenMicro.uniformity = sharpenMicro.uniformity && p.sharpenMicro.uniformity == other.sharpenMicro.uniformity;
        sharpening.enabled = sharpening.enabled && p.sharpening.enabled == other.sharpening.enabled;
        sharpening.radius = sharpening.radius && p.sharpening.radius == other.sharpening.radius;
        sharpening.amount = sharpening.amount && p.sharpening.amount == other.sharpening.amount;
        sharpening.threshold = sharpening.threshold && p.sharpening.threshold == other.sharpening.threshold;
        sharpening.edgesonly = sharpening.edgesonly && p.sharpening.edgesonly == other.sharpening.edgesonly;
        sharpening.edges_radius = sharpening.edges_radius && p.sharpening.edges_radius == other.sharpening.edges_radius;
        sharpening.edges_tolerance = sharpening.edges_tolerance && p.sharpening.edges_tolerance == other.sharpening.edges_tolerance;
        sharpening.halocontrol = sharpening.halocontrol && p.sharpening.halocontrol == other.sharpening.halocontrol;
        sharpening.halocontrol_amount = sharpening.halocontrol_amount && p.sharpening.halocontrol_amount == other.sharpening.halocontrol_amount;
        sharpening.method = sharpening.method && p.sharpening.method == other.sharpening.method;
        sharpening.deconvamount = sharpening.deconvamount && p.sharpening.deconvamount == other.sharpening.deconvamount;
        sharpening.deconvradius = sharpening.deconvradius && p.sharpening.deconvradius == other.sharpening.deconvradius;
        sharpening.deconviter = sharpening.deconviter && p.sharpening.deconviter == other.sharpening.deconviter;
        sharpening.deconvdamping = sharpening.deconvdamping && p.sharpening.deconvdamping == other.sharpening.deconvdamping;
        vibrance.enabled = vibrance.enabled && p.vibrance.enabled == other.vibrance.enabled;
        vibrance.pastels = vibrance.pastels && p.vibrance.pastels == other.vibrance.pastels;
        vibrance.saturated = vibrance.saturated && p.vibrance.saturated == other.vibrance.saturated;
        vibrance.psthreshold = vibrance.psthreshold && p.vibrance.psthreshold == other.vibrance.psthreshold;
        vibrance.protectskins = vibrance.protectskins && p.vibrance.protectskins == other.vibrance.protectskins;
        vibrance.avoidcolorshift = vibrance.avoidcolorshift && p.vibrance.avoidcolorshift == other.vibrance.avoidcolorshift;
        vibrance.pastsattog = vibrance.pastsattog && p.vibrance.pastsattog == other.vibrance.pastsattog;
        vibrance.skintonescurve = vibrance.skintonescurve && p.vibrance.skintonescurve == other.vibrance.skintonescurve;

        colorappearance.enabled = colorappearance.enabled && p.colorappearance.enabled == other.colorappearance.enabled;
        colorappearance.degree = colorappearance.degree && p.colorappearance.degree == other.colorappearance.degree;
        colorappearance.autodegree = colorappearance.autodegree && p.colorappearance.autodegree == other.colorappearance.autodegree;
        colorappearance.surround = colorappearance.surround && p.colorappearance.surround == other.colorappearance.surround;
        colorappearance.adapscen = colorappearance.adapscen && p.colorappearance.adapscen == other.colorappearance.adapscen;
        colorappearance.autoadapscen = colorappearance.autoadapscen && p.colorappearance.autoadapscen == other.colorappearance.autoadapscen;
        colorappearance.adaplum = colorappearance.adaplum && p.colorappearance.adaplum == other.colorappearance.adaplum;
        colorappearance.badpixsl = colorappearance.badpixsl && p.colorappearance.badpixsl == other.colorappearance.badpixsl;
        colorappearance.wbmodel = colorappearance.wbmodel && p.colorappearance.wbmodel == other.colorappearance.wbmodel;
        colorappearance.algo = colorappearance.algo && p.colorappearance.algo == other.colorappearance.algo;
        colorappearance.jlight = colorappearance.jlight && p.colorappearance.jlight == other.colorappearance.jlight;
        colorappearance.qbright = colorappearance.qbright && p.colorappearance.qbright == other.colorappearance.qbright;
        colorappearance.chroma = colorappearance.chroma && p.colorappearance.chroma == other.colorappearance.chroma;
        colorappearance.schroma = colorappearance.schroma && p.colorappearance.schroma == other.colorappearance.schroma;
        colorappearance.mchroma = colorappearance.mchroma && p.colorappearance.mchroma == other.colorappearance.mchroma;
        colorappearance.rstprotection = colorappearance.rstprotection && p.colorappearance.rstprotection == other.colorappearance.rstprotection;
        colorappearance.contrast = colorappearance.contrast && p.colorappearance.contrast == other.colorappearance.contrast;
        colorappearance.qcontrast = colorappearance.qcontrast && p.colorappearance.qcontrast == other.colorappearance.qcontrast;
        colorappearance.colorh = colorappearance.colorh && p.colorappearance.colorh == other.colorappearance.colorh;
        colorappearance.surrsource = colorappearance.surrsource && p.colorappearance.surrsource == other.colorappearance.surrsource;
        colorappearance.gamut = colorappearance.gamut && p.colorappearance.gamut == other.colorappearance.gamut;
 //       colorappearance.badpix = colorappearance.badpix && p.colorappearance.badpix == other.colorappearance.badpix;
        colorappearance.datacie = colorappearance.datacie && p.colorappearance.datacie == other.colorappearance.datacie;
        colorappearance.tonecie = colorappearance.tonecie && p.colorappearance.tonecie == other.colorappearance.tonecie;
   //     colorappearance.sharpcie = colorappearance.sharpcie && p.colorappearance.sharpcie == other.colorappearance.sharpcie;
        colorappearance.curve = colorappearance.curve && p.colorappearance.curve == other.colorappearance.curve;
        colorappearance.curve3 = colorappearance.curve3 && p.colorappearance.curve3 == other.colorappearance.curve3;
        colorappearance.curve2 = colorappearance.curve2 && p.colorappearance.curve2 == other.colorappearance.curve2;
        colorappearance.curveMode = colorappearance.curveMode && p.colorappearance.curveMode == other.colorappearance.curveMode;
        colorappearance.curveMode2 = colorappearance.curveMode2 && p.colorappearance.curveMode2 == other.colorappearance.curveMode2;
        colorappearance.curveMode3 = colorappearance.curveMode3 && p.colorappearance.curveMode3 == other.colorappearance.curveMode3;

        //colorBoost.amount = colorBoost.amount && p.colorBoost.amount == other.colorBoost.amount;
        //colorBoost.avoidclip = colorBoost.avoidclip && p.colorBoost.avoidclip == other.colorBoost.avoidclip;
        //colorBoost.enable_saturationlimiter = colorBoost.enable_saturationlimiter && p.colorBoost.enable_saturationlimiter == other.colorBoost.enable_saturationlimiter;
        //colorBoost.saturationlimit = colorBoost.saturationlimit && p.colorBoost.saturationlimit == other.colorBoost.saturationlimit;
        wb.method = wb.method && p.wb.method == other.wb.method;
        wb.green = wb.green && p.wb.green == other.wb.green;
        wb.equal = wb.equal && p.wb.equal == other.wb.equal;
        wb.temperature = wb.temperature && p.wb.temperature == other.wb.temperature;
        //colorShift.a = colorShift.a && p.colorShift.a == other.colorShift.a;
        //colorShift.b = colorShift.b && p.colorShift.b == other.colorShift.b;
        //lumaDenoise.enabled = lumaDenoise.enabled && p.lumaDenoise.enabled == other.lumaDenoise.enabled;
        //lumaDenoise.radius = lumaDenoise.radius && p.lumaDenoise.radius == other.lumaDenoise.radius;
        //lumaDenoise.edgetolerance = lumaDenoise.edgetolerance && p.lumaDenoise.edgetolerance == other.lumaDenoise.edgetolerance;
        //colorDenoise.enabled = colorDenoise.enabled && p.colorDenoise.enabled == other.colorDenoise.enabled;
        //colorDenoise.amount = colorDenoise.amount && p.colorDenoise.amount == other.colorDenoise.amount;
        defringe.enabled = defringe.enabled && p.defringe.enabled == other.defringe.enabled;
        defringe.radius = defringe.radius && p.defringe.radius == other.defringe.radius;
        defringe.threshold = defringe.threshold && p.defringe.threshold == other.defringe.threshold;
        defringe.huecurve = defringe.huecurve && p.defringe.huecurve == other.defringe.huecurve;

        impulseDenoise.enabled = impulseDenoise.enabled && p.impulseDenoise.enabled == other.impulseDenoise.enabled;
        impulseDenoise.thresh = impulseDenoise.thresh && p.impulseDenoise.thresh == other.impulseDenoise.thresh;

        dirpyrDenoise.enabled = dirpyrDenoise.enabled && p.dirpyrDenoise.enabled == other.dirpyrDenoise.enabled;
        dirpyrDenoise.enhance = dirpyrDenoise.enhance && p.dirpyrDenoise.enhance == other.dirpyrDenoise.enhance;
        dirpyrDenoise.median = dirpyrDenoise.median && p.dirpyrDenoise.median == other.dirpyrDenoise.median;
        dirpyrDenoise.autochroma = dirpyrDenoise.autochroma && p.dirpyrDenoise.autochroma == other.dirpyrDenoise.autochroma;
 //       dirpyrDenoise.perform = dirpyrDenoise.perform && p.dirpyrDenoise.perform == other.dirpyrDenoise.perform;
        dirpyrDenoise.luma = dirpyrDenoise.luma && p.dirpyrDenoise.luma == other.dirpyrDenoise.luma;
        dirpyrDenoise.lcurve = dirpyrDenoise.lcurve && p.dirpyrDenoise.lcurve == other.dirpyrDenoise.lcurve;
        dirpyrDenoise.cccurve = dirpyrDenoise.cccurve && p.dirpyrDenoise.cccurve == other.dirpyrDenoise.cccurve;
        dirpyrDenoise.Ldetail = dirpyrDenoise.Ldetail && p.dirpyrDenoise.Ldetail == other.dirpyrDenoise.Ldetail;
        dirpyrDenoise.chroma = dirpyrDenoise.chroma && p.dirpyrDenoise.chroma == other.dirpyrDenoise.chroma;
        dirpyrDenoise.redchro = dirpyrDenoise.redchro && p.dirpyrDenoise.redchro == other.dirpyrDenoise.redchro;
        dirpyrDenoise.bluechro = dirpyrDenoise.bluechro && p.dirpyrDenoise.bluechro == other.dirpyrDenoise.bluechro;	
        dirpyrDenoise.gamma = dirpyrDenoise.gamma && p.dirpyrDenoise.gamma == other.dirpyrDenoise.gamma;
        dirpyrDenoise.passes = dirpyrDenoise.passes && p.dirpyrDenoise.passes == other.dirpyrDenoise.passes;
        dirpyrDenoise.dmethod = dirpyrDenoise.dmethod && p.dirpyrDenoise.dmethod == other.dirpyrDenoise.dmethod;
        dirpyrDenoise.Lmethod = dirpyrDenoise.Lmethod && p.dirpyrDenoise.Lmethod == other.dirpyrDenoise.Lmethod;
        dirpyrDenoise.Cmethod = dirpyrDenoise.Cmethod && p.dirpyrDenoise.Cmethod == other.dirpyrDenoise.Cmethod;
        dirpyrDenoise.C2method = dirpyrDenoise.C2method && p.dirpyrDenoise.C2method == other.dirpyrDenoise.C2method;
        dirpyrDenoise.smethod = dirpyrDenoise.smethod && p.dirpyrDenoise.smethod == other.dirpyrDenoise.smethod;
        dirpyrDenoise.medmethod = dirpyrDenoise.medmethod && p.dirpyrDenoise.medmethod == other.dirpyrDenoise.medmethod;
        dirpyrDenoise.methodmed = dirpyrDenoise.methodmed && p.dirpyrDenoise.methodmed == other.dirpyrDenoise.methodmed;
        dirpyrDenoise.rgbmethod = dirpyrDenoise.rgbmethod && p.dirpyrDenoise.rgbmethod == other.dirpyrDenoise.rgbmethod;

        epd.enabled = epd.enabled && p.epd.enabled == other.epd.enabled;
        epd.strength = epd.strength && p.epd.strength == other.epd.strength;
        epd.edgeStopping = epd.edgeStopping && p.epd.edgeStopping == other.epd.edgeStopping;
        epd.scale = epd.scale && p.epd.scale == other.epd.scale;
        epd.reweightingIterates = epd.reweightingIterates && p.epd.reweightingIterates == other.epd.reweightingIterates;

        sh.enabled = sh.enabled && p.sh.enabled == other.sh.enabled;
        sh.hq = sh.hq && p.sh.hq == other.sh.hq;
        sh.highlights = sh.highlights && p.sh.highlights == other.sh.highlights;
        sh.htonalwidth = sh.htonalwidth && p.sh.htonalwidth == other.sh.htonalwidth;
        sh.shadows = sh.shadows && p.sh.shadows == other.sh.shadows;
        sh.stonalwidth = sh.stonalwidth && p.sh.stonalwidth == other.sh.stonalwidth;
        sh.localcontrast = sh.localcontrast && p.sh.localcontrast == other.sh.localcontrast;
        sh.radius = sh.radius && p.sh.radius == other.sh.radius;
        crop.enabled = crop.enabled && p.crop.enabled == other.crop.enabled;
        crop.x = crop.x && p.crop.x == other.crop.x;
        crop.y = crop.y && p.crop.y == other.crop.y;
        crop.w = crop.w && p.crop.w == other.crop.w;
        crop.h = crop.h && p.crop.h == other.crop.h;
        crop.fixratio = crop.fixratio && p.crop.fixratio == other.crop.fixratio;
        crop.ratio = crop.ratio && p.crop.ratio == other.crop.ratio;
        crop.orientation = crop.orientation && p.crop.orientation == other.crop.orientation;
        crop.guide = crop.guide && p.crop.guide == other.crop.guide;
        coarse.rotate = coarse.rotate && p.coarse.rotate == other.coarse.rotate;
        coarse.hflip = coarse.hflip && p.coarse.hflip == other.coarse.hflip;
        coarse.vflip = coarse.vflip && p.coarse.vflip == other.coarse.vflip;
        commonTrans.autofill = commonTrans.autofill && p.commonTrans.autofill == other.commonTrans.autofill;
        rotate.degree = rotate.degree && p.rotate.degree == other.rotate.degree;
        distortion.amount = distortion.amount && p.distortion.amount == other.distortion.amount;
        lensProf.lcpFile = lensProf.lcpFile && p.lensProf.lcpFile == other.lensProf.lcpFile;
        lensProf.useDist = lensProf.useDist && p.lensProf.useDist == other.lensProf.useDist;
        lensProf.useVign = lensProf.useVign && p.lensProf.useVign == other.lensProf.useVign;
        lensProf.useCA = lensProf.useCA && p.lensProf.useCA == other.lensProf.useCA;
        perspective.horizontal = perspective.horizontal && p.perspective.horizontal == other.perspective.horizontal;
        perspective.vertical = perspective.vertical && p.perspective.vertical == other.perspective.vertical;
        gradient.enabled = gradient.enabled && p.gradient.enabled == other.gradient.enabled;
        gradient.degree = gradient.degree && p.gradient.degree == other.gradient.degree;
        gradient.feather = gradient.feather && p.gradient.feather == other.gradient.feather;
        gradient.strength = gradient.strength && p.gradient.strength == other.gradient.strength;
        gradient.centerX = gradient.centerX && p.gradient.centerX == other.gradient.centerX;
        gradient.centerY = gradient.centerY && p.gradient.centerY == other.gradient.centerY;
        pcvignette.enabled = pcvignette.enabled && p.pcvignette.enabled == other.pcvignette.enabled;
        pcvignette.strength = pcvignette.strength && p.pcvignette.strength == other.pcvignette.strength;
        pcvignette.feather = pcvignette.feather && p.pcvignette.feather == other.pcvignette.feather;
        pcvignette.roundness = pcvignette.roundness && p.pcvignette.roundness == other.pcvignette.roundness;
        cacorrection.red = cacorrection.red && p.cacorrection.red == other.cacorrection.red;
        cacorrection.blue = cacorrection.blue && p.cacorrection.blue == other.cacorrection.blue;
        vignetting.amount = vignetting.amount && p.vignetting.amount == other.vignetting.amount;
        vignetting.radius = vignetting.radius && p.vignetting.radius == other.vignetting.radius;
        vignetting.strength = vignetting.strength && p.vignetting.strength == other.vignetting.strength;
        vignetting.centerX = vignetting.centerX && p.vignetting.centerX == other.vignetting.centerX;
        vignetting.centerY = vignetting.centerY && p.vignetting.centerY == other.vignetting.centerY;
        chmixer.red[0] = chmixer.red[0] && p.chmixer.red[0] == other.chmixer.red[0];
        chmixer.red[1] = chmixer.red[1] && p.chmixer.red[1] == other.chmixer.red[1];
        chmixer.red[2] = chmixer.red[2] && p.chmixer.red[2] == other.chmixer.red[2];
        chmixer.green[0] = chmixer.green[0] && p.chmixer.green[0] == other.chmixer.green[0];
        chmixer.green[1] = chmixer.green[1] && p.chmixer.green[1] == other.chmixer.green[1];
        chmixer.green[2] = chmixer.green[2] && p.chmixer.green[2] == other.chmixer.green[2];
        chmixer.blue[0] = chmixer.blue[0] && p.chmixer.blue[0] == other.chmixer.blue[0];
        chmixer.blue[1] = chmixer.blue[1] && p.chmixer.blue[1] == other.chmixer.blue[1];
        chmixer.blue[2] = chmixer.blue[2] && p.chmixer.blue[2] == other.chmixer.blue[2];
        blackwhite.enabledcc = blackwhite.enabledcc && p.blackwhite.enabledcc == other.blackwhite.enabledcc;
        blackwhite.enabled = blackwhite.enabled && p.blackwhite.enabled == other.blackwhite.enabled;
        blackwhite.mixerRed = blackwhite.mixerRed && p.blackwhite.mixerRed == other.blackwhite.mixerRed;
        blackwhite.mixerOrange = blackwhite.mixerOrange && p.blackwhite.mixerOrange == other.blackwhite.mixerOrange;
        blackwhite.mixerYellow = blackwhite.mixerYellow && p.blackwhite.mixerYellow == other.blackwhite.mixerYellow;
        blackwhite.mixerGreen = blackwhite.mixerGreen && p.blackwhite.mixerGreen == other.blackwhite.mixerGreen;
        blackwhite.mixerCyan = blackwhite.mixerCyan && p.blackwhite.mixerCyan == other.blackwhite.mixerCyan;
        blackwhite.mixerBlue = blackwhite.mixerBlue && p.blackwhite.mixerBlue == other.blackwhite.mixerBlue;
        blackwhite.mixerMagenta = blackwhite.mixerMagenta && p.blackwhite.mixerMagenta == other.blackwhite.mixerMagenta;
        blackwhite.mixerPurple = blackwhite.mixerPurple && p.blackwhite.mixerPurple == other.blackwhite.mixerPurple;
        blackwhite.gammaRed = blackwhite.gammaRed && p.blackwhite.gammaRed == other.blackwhite.gammaRed;
        blackwhite.gammaGreen = blackwhite.gammaGreen && p.blackwhite.gammaGreen == other.blackwhite.gammaGreen;
        blackwhite.gammaBlue = blackwhite.gammaBlue && p.blackwhite.gammaBlue == other.blackwhite.gammaBlue;
        blackwhite.filter = blackwhite.filter && p.blackwhite.filter == other.blackwhite.filter;
        blackwhite.setting = blackwhite.setting && p.blackwhite.setting == other.blackwhite.setting;
        blackwhite.luminanceCurve = blackwhite.luminanceCurve && p.blackwhite.luminanceCurve == other.blackwhite.luminanceCurve;
        blackwhite.method = blackwhite.method && p.blackwhite.method == other.blackwhite.method;
        blackwhite.beforeCurve = blackwhite.beforeCurve && p.blackwhite.beforeCurve == other.blackwhite.beforeCurve;
        blackwhite.beforeCurveMode = blackwhite.beforeCurveMode && p.blackwhite.beforeCurveMode == other.blackwhite.beforeCurveMode;
        blackwhite.afterCurve = blackwhite.afterCurve && p.blackwhite.afterCurve == other.blackwhite.afterCurve;
        blackwhite.afterCurveMode = blackwhite.afterCurveMode && p.blackwhite.afterCurveMode == other.blackwhite.afterCurveMode;
        blackwhite.autoc = blackwhite.autoc && p.blackwhite.autoc == other.blackwhite.autoc;
        blackwhite.algo = blackwhite.algo && p.blackwhite.algo == other.blackwhite.algo;
        resize.scale = resize.scale && p.resize.scale == other.resize.scale;
        resize.appliesTo = resize.appliesTo && p.resize.appliesTo == other.resize.appliesTo;
        resize.method = resize.method && p.resize.method == other.resize.method;
        resize.dataspec = resize.dataspec && p.resize.dataspec == other.resize.dataspec;
        resize.width = resize.width && p.resize.width == other.resize.width;
        resize.height = resize.height && p.resize.height == other.resize.height;
        resize.enabled = resize.enabled && p.resize.enabled == other.resize.enabled;
        icm.input = icm.input && p.icm.input == other.icm.input;
        icm.toneCurve = icm.toneCurve && p.icm.toneCurve == other.icm.toneCurve;
        icm.blendCMSMatrix = icm.blendCMSMatrix && p.icm.blendCMSMatrix == other.icm.blendCMSMatrix;
        icm.dcpIlluminant = icm.dcpIlluminant && p.icm.dcpIlluminant == other.icm.dcpIlluminant;
        icm.working = icm.working && p.icm.working == other.icm.working;
        icm.output = icm.output && p.icm.output == other.icm.output;
        icm.gamma = icm.gamma && p.icm.gamma == other.icm.gamma;
        icm.freegamma = icm.freegamma && p.icm.freegamma == other.icm.freegamma;
        icm.gampos = icm.gampos && p.icm.gampos == other.icm.gampos;
        icm.slpos = icm.slpos && p.icm.slpos == other.icm.slpos;
        raw.bayersensor.method = raw.bayersensor.method && p.raw.bayersensor.method == other.raw.bayersensor.method;
        raw.bayersensor.ccSteps = raw.bayersensor.ccSteps && p.raw.bayersensor.ccSteps == other.raw.bayersensor.ccSteps;
        raw.bayersensor.exBlack0 = raw.bayersensor.exBlack0 && p.raw.bayersensor.black0 == other.raw.bayersensor.black0;
        raw.bayersensor.exBlack1 = raw.bayersensor.exBlack1 && p.raw.bayersensor.black1 == other.raw.bayersensor.black1;
        raw.bayersensor.exBlack2 = raw.bayersensor.exBlack2 && p.raw.bayersensor.black2 == other.raw.bayersensor.black2;
        raw.bayersensor.exBlack3 = raw.bayersensor.exBlack3 && p.raw.bayersensor.black3 == other.raw.bayersensor.black3;
        raw.bayersensor.exTwoGreen = raw.bayersensor.exTwoGreen && p.raw.bayersensor.twogreen == other.raw.bayersensor.twogreen;
        raw.bayersensor.dcbIterations = raw.bayersensor.dcbIterations && p.raw.bayersensor.dcb_iterations == other.raw.bayersensor.dcb_iterations;
        raw.bayersensor.dcbEnhance = raw.bayersensor.dcbEnhance && p.raw.bayersensor.dcb_enhance == other.raw.bayersensor.dcb_enhance;
        //raw.bayersensor.allEnhance = raw.bayersensor.allEnhance && p.raw.bayersensor.all_enhance == other.raw.bayersensor.all_enhance;
        raw.bayersensor.lmmseIterations = raw.bayersensor.lmmseIterations && p.raw.bayersensor.lmmse_iterations == other.raw.bayersensor.lmmse_iterations;
        raw.bayersensor.greenEq = raw.bayersensor.greenEq && p.raw.bayersensor.greenthresh == other.raw.bayersensor.greenthresh;
        raw.bayersensor.linenoise = raw.bayersensor.linenoise && p.raw.bayersensor.linenoise == other.raw.bayersensor.linenoise;
        raw.xtranssensor.method = raw.xtranssensor.method && p.raw.xtranssensor.method == other.raw.xtranssensor.method;
        raw.xtranssensor.ccSteps = raw.xtranssensor.ccSteps && p.raw.xtranssensor.ccSteps == other.raw.xtranssensor.ccSteps;
        raw.xtranssensor.exBlackRed = raw.xtranssensor.exBlackRed && p.raw.xtranssensor.blackred == other.raw.xtranssensor.blackred;
        raw.xtranssensor.exBlackGreen = raw.xtranssensor.exBlackGreen && p.raw.xtranssensor.blackgreen == other.raw.xtranssensor.blackgreen;
        raw.xtranssensor.exBlackBlue = raw.xtranssensor.exBlackBlue && p.raw.xtranssensor.blackblue == other.raw.xtranssensor.blackblue;
        raw.caCorrection = raw.caCorrection && p.raw.ca_autocorrect == other.raw.ca_autocorrect;
        raw.caRed = raw.caRed && p.raw.cared == other.raw.cared;
        raw.caBlue = raw.caBlue && p.raw.cablue == other.raw.cablue;
        raw.hotPixelFilter = raw.hotPixelFilter && p.raw.hotPixelFilter == other.raw.hotPixelFilter;
        raw.deadPixelFilter = raw.deadPixelFilter && p.raw.deadPixelFilter == other.raw.deadPixelFilter;
        raw.hotDeadPixelThresh = raw.hotDeadPixelThresh && p.raw.hotdeadpix_thresh == other.raw.hotdeadpix_thresh;
        raw.darkFrame = raw.darkFrame && p.raw.dark_frame == other.raw.dark_frame;
        raw.dfAuto = raw.dfAuto && p.raw.df_autoselect == other.raw.df_autoselect;
        raw.ff_file = raw.ff_file && p.raw.ff_file == other.raw.ff_file;
        raw.ff_AutoSelect = raw.ff_AutoSelect && p.raw.ff_AutoSelect == other.raw.ff_AutoSelect;
        raw.ff_BlurRadius = raw.ff_BlurRadius && p.raw.ff_BlurRadius == other.raw.ff_BlurRadius;
        raw.ff_BlurType = raw.ff_BlurType && p.raw.ff_BlurType == other.raw.ff_BlurType;
        raw.ff_AutoClipControl = raw.ff_AutoClipControl && p.raw.ff_AutoClipControl == other.raw.ff_AutoClipControl;
        raw.ff_clipControl = raw.ff_clipControl && p.raw.ff_clipControl == other.raw.ff_clipControl;
        raw.exPos = raw.exPos && p.raw.expos == other.raw.expos;
        raw.exPreser = raw.exPreser && p.raw.preser == other.raw.preser;
        wavelet.enabled = wavelet.enabled && p.wavelet.enabled == other.wavelet.enabled;
        wavelet.median = wavelet.median && p.wavelet.median == other.wavelet.median;
        wavelet.avoid = wavelet.avoid && p.wavelet.avoid == other.wavelet.avoid;
        wavelet.Lmethod = wavelet.Lmethod && p.wavelet.Lmethod == other.wavelet.Lmethod;
        wavelet.CLmethod = wavelet.CLmethod && p.wavelet.CLmethod == other.wavelet.CLmethod;
        wavelet.Tilesmethod = wavelet.Tilesmethod && p.wavelet.Tilesmethod == other.wavelet.Tilesmethod;
        wavelet.CHmethod = wavelet.CHmethod && p.wavelet.CHmethod == other.wavelet.CHmethod;
        wavelet.HSmethod = wavelet.HSmethod && p.wavelet.HSmethod == other.wavelet.HSmethod;
        wavelet.Dirmethod = wavelet.Dirmethod && p.wavelet.Dirmethod == other.wavelet.Dirmethod;
        wavelet.tiles = wavelet.tiles && p.wavelet.tiles == other.wavelet.tiles;
        wavelet.rescon = wavelet.rescon && p.wavelet.rescon == other.wavelet.rescon;
        wavelet.resconH = wavelet.resconH && p.wavelet.resconH == other.wavelet.resconH;
        wavelet.reschro = wavelet.reschro && p.wavelet.reschro == other.wavelet.reschro;
        wavelet.sup = wavelet.sup && p.wavelet.sup == other.wavelet.sup;
        wavelet.sky = wavelet.sky && p.wavelet.sky == other.wavelet.sky;
        wavelet.threshold = wavelet.threshold && p.wavelet.threshold == other.wavelet.threshold;
        wavelet.threshold2 = wavelet.threshold2 && p.wavelet.threshold2 == other.wavelet.threshold2;
        wavelet.thres = wavelet.thres && p.wavelet.thres == other.wavelet.thres;
        wavelet.chroma = wavelet.chroma && p.wavelet.chroma == other.wavelet.chroma;
        wavelet.chro = wavelet.chro && p.wavelet.chro == other.wavelet.chro;
        wavelet.unif = wavelet.unif && p.wavelet.unif == other.wavelet.unif;
        wavelet.thr = wavelet.thr && p.wavelet.thr == other.wavelet.thr;
        wavelet.thrH = wavelet.thrH && p.wavelet.thrH == other.wavelet.thrH;
        wavelet.hueskin = wavelet.hueskin && p.wavelet.hueskin == other.wavelet.hueskin;
        wavelet.hueskin2 = wavelet.hueskin2 && p.wavelet.hueskin2 == other.wavelet.hueskin2;
        wavelet.hllev = wavelet.hllev && p.wavelet.hllev == other.wavelet.hllev;
        wavelet.bllev = wavelet.bllev && p.wavelet.bllev == other.wavelet.bllev;
        wavelet.pastlev = wavelet.pastlev && p.wavelet.pastlev == other.wavelet.pastlev;
        wavelet.satlev = wavelet.satlev && p.wavelet.satlev == other.wavelet.satlev;
        wavelet.clvcurve = wavelet.clvcurve && p.wavelet.clvcurve == other.wavelet.clvcurve;
        wavelet.opacityCurveRG = wavelet.opacityCurveRG && p.wavelet.opacityCurveRG == other.wavelet.opacityCurveRG;
        wavelet.opacityCurveBY = wavelet.opacityCurveBY && p.wavelet.opacityCurveBY == other.wavelet.opacityCurveBY;
        wavelet.skinprotect = wavelet.skinprotect && p.wavelet.skinprotect == other.wavelet.skinprotect;
        for(int i = 0; i < 9; i++) {
            wavelet.c[i] = wavelet.c[i] && p.wavelet.c[i] == other.wavelet.c[i];
        }

        dirpyrequalizer.enabled = dirpyrequalizer.enabled && p.dirpyrequalizer.enabled == other.dirpyrequalizer.enabled;
        dirpyrequalizer.gamutlab = dirpyrequalizer.gamutlab && p.dirpyrequalizer.gamutlab == other.dirpyrequalizer.gamutlab;
        for(int i = 0; i < 6; i++) {
            dirpyrequalizer.mult[i] = dirpyrequalizer.mult[i] && p.dirpyrequalizer.mult[i] == other.dirpyrequalizer.mult[i];
        }
        dirpyrequalizer.threshold = dirpyrequalizer.threshold && p.dirpyrequalizer.threshold == other.dirpyrequalizer.threshold;
        dirpyrequalizer.skinprotect = dirpyrequalizer.skinprotect && p.dirpyrequalizer.skinprotect == other.dirpyrequalizer.skinprotect;
    //    dirpyrequalizer.algo = dirpyrequalizer.algo && p.dirpyrequalizer.algo == other.dirpyrequalizer.algo;
        dirpyrequalizer.hueskin = dirpyrequalizer.hueskin && p.dirpyrequalizer.hueskin == other.dirpyrequalizer.hueskin;
        hsvequalizer.hcurve = hsvequalizer.hcurve && p.hsvequalizer.hcurve == other.hsvequalizer.hcurve;
        hsvequalizer.scurve = hsvequalizer.scurve && p.hsvequalizer.scurve == other.hsvequalizer.scurve;
        hsvequalizer.vcurve = hsvequalizer.vcurve && p.hsvequalizer.vcurve == other.hsvequalizer.vcurve;
        filmSimulation.enabled = filmSimulation.enabled && p.filmSimulation.enabled == other.filmSimulation.enabled;
        filmSimulation.clutFilename = filmSimulation.clutFilename && p.filmSimulation.clutFilename == other.filmSimulation.clutFilename;
        filmSimulation.strength = filmSimulation.strength && p.filmSimulation.strength == other.filmSimulation.strength;

//      How the hell can we handle that???
//      exif = exif && p.exif==other.exif
//      iptc = other.iptc;
    }
}

void ParamsEdited::combine (rtengine::procparams::ProcParams& toEdit, const rtengine::procparams::ProcParams& mods, bool forceSet) {

	bool dontforceSet = !forceSet;

	if (toneCurve.curve)		toEdit.toneCurve.curve      = mods.toneCurve.curve;
	if (toneCurve.curve2)		toEdit.toneCurve.curve2     = mods.toneCurve.curve2;
	if (toneCurve.curveMode)	toEdit.toneCurve.curveMode  = mods.toneCurve.curveMode;
	if (toneCurve.curveMode2)	toEdit.toneCurve.curveMode2 = mods.toneCurve.curveMode2;
	if (toneCurve.brightness)	toEdit.toneCurve.brightness = dontforceSet && options.baBehav[ADDSET_TC_BRIGHTNESS] ? toEdit.toneCurve.brightness + mods.toneCurve.brightness : mods.toneCurve.brightness;
	if (toneCurve.black)		toEdit.toneCurve.black 	    = dontforceSet && options.baBehav[ADDSET_TC_BLACKLEVEL] ? toEdit.toneCurve.black + mods.toneCurve.black : mods.toneCurve.black;
	if (toneCurve.contrast)		toEdit.toneCurve.contrast 	= dontforceSet && options.baBehav[ADDSET_TC_CONTRAST] ? toEdit.toneCurve.contrast + mods.toneCurve.contrast : mods.toneCurve.contrast;
	if (toneCurve.saturation)	toEdit.toneCurve.saturation = dontforceSet && options.baBehav[ADDSET_TC_SATURATION] ? toEdit.toneCurve.saturation + mods.toneCurve.saturation : mods.toneCurve.saturation;
	if (toneCurve.shcompr)		toEdit.toneCurve.shcompr 	= dontforceSet && options.baBehav[ADDSET_TC_SHCOMP] ? toEdit.toneCurve.shcompr + mods.toneCurve.shcompr : mods.toneCurve.shcompr;
	if (toneCurve.autoexp)		toEdit.toneCurve.autoexp 	= mods.toneCurve.autoexp;
	if (toneCurve.clip)			toEdit.toneCurve.clip 	    = mods.toneCurve.clip;
	if (toneCurve.expcomp)		toEdit.toneCurve.expcomp 	= dontforceSet && options.baBehav[ADDSET_TC_EXPCOMP] ? toEdit.toneCurve.expcomp + mods.toneCurve.expcomp : mods.toneCurve.expcomp;
	if (toneCurve.hlcompr)		toEdit.toneCurve.hlcompr 	= dontforceSet && options.baBehav[ADDSET_TC_HLCOMPAMOUNT] ? toEdit.toneCurve.hlcompr + mods.toneCurve.hlcompr : mods.toneCurve.hlcompr;
	if (toneCurve.hlcomprthresh) toEdit.toneCurve.hlcomprthresh	= dontforceSet && options.baBehav[ADDSET_TC_HLCOMPTHRESH] ? toEdit.toneCurve.hlcomprthresh + mods.toneCurve.hlcomprthresh : mods.toneCurve.hlcomprthresh;
	if (toneCurve.hrenabled)	toEdit.toneCurve.hrenabled 	= mods.toneCurve.hrenabled;
	if (toneCurve.method)		toEdit.toneCurve.method 	= mods.toneCurve.method;
	if (labCurve.lcurve)		toEdit.labCurve.lcurve 	    = mods.labCurve.lcurve;
	if (labCurve.acurve)		toEdit.labCurve.acurve 	    = mods.labCurve.acurve;
	if (labCurve.bcurve)		toEdit.labCurve.bcurve 	    = mods.labCurve.bcurve;
	if (labCurve.cccurve)		toEdit.labCurve.cccurve     = mods.labCurve.cccurve;
	if (labCurve.chcurve)		toEdit.labCurve.chcurve     = mods.labCurve.chcurve;
	if (labCurve.lhcurve)		toEdit.labCurve.lhcurve     = mods.labCurve.lhcurve;
	if (labCurve.hhcurve)		toEdit.labCurve.hhcurve     = mods.labCurve.hhcurve;
	if (labCurve.lccurve)		toEdit.labCurve.lccurve    = mods.labCurve.lccurve;
	if (labCurve.clcurve)		toEdit.labCurve.clcurve    = mods.labCurve.clcurve;
	if (labCurve.brightness)	toEdit.labCurve.brightness   = dontforceSet && options.baBehav[ADDSET_LC_BRIGHTNESS] ? toEdit.labCurve.brightness + mods.labCurve.brightness : mods.labCurve.brightness;
	if (labCurve.contrast)		toEdit.labCurve.contrast 	 = dontforceSet && options.baBehav[ADDSET_LC_CONTRAST] ? toEdit.labCurve.contrast + mods.labCurve.contrast : mods.labCurve.contrast;
	if (labCurve.chromaticity)	toEdit.labCurve.chromaticity = dontforceSet && options.baBehav[ADDSET_LC_CHROMATICITY] ? toEdit.labCurve.chromaticity + mods.labCurve.chromaticity : mods.labCurve.chromaticity;
	if (labCurve.avoidcolorshift)	toEdit.labCurve.avoidcolorshift		= mods.labCurve.avoidcolorshift;
	if (labCurve.rstprotection)		toEdit.labCurve.rstprotection		= mods.labCurve.rstprotection;
	if (labCurve.lcredsk)			toEdit.labCurve.lcredsk				= mods.labCurve.lcredsk;

	if (rgbCurves.lumamode)					toEdit.rgbCurves.lumamode   = mods.rgbCurves.lumamode;
	if (rgbCurves.rcurve)					toEdit.rgbCurves.rcurve     = mods.rgbCurves.rcurve;
	if (rgbCurves.gcurve)					toEdit.rgbCurves.gcurve     = mods.rgbCurves.gcurve;
	if (rgbCurves.bcurve)					toEdit.rgbCurves.bcurve     = mods.rgbCurves.bcurve;

	if (colorToning.enabled)				toEdit.colorToning.enabled		= mods.colorToning.enabled;
	if (colorToning.twocolor)				toEdit.colorToning.twocolor		= mods.colorToning.twocolor;
	if (colorToning.opacityCurve)			toEdit.colorToning.opacityCurve	= mods.colorToning.opacityCurve;
	if (colorToning.colorCurve)				toEdit.colorToning.colorCurve	= mods.colorToning.colorCurve;
	if (colorToning.enabled)				toEdit.colorToning.enabled		= mods.colorToning.enabled;
	if (colorToning.opacityCurve)			toEdit.colorToning.opacityCurve	= mods.colorToning.opacityCurve;
	if (colorToning.satprotectionthreshold)	toEdit.colorToning.satProtectionThreshold	= dontforceSet && options.baBehav[ADDSET_COLORTONING_SATTHRESHOLD] ? toEdit.colorToning.satProtectionThreshold + mods.colorToning.satProtectionThreshold : mods.colorToning.satProtectionThreshold;
	if (colorToning.autosat)				toEdit.colorToning.autosat					= mods.colorToning.autosat;
	if (colorToning.saturatedopacity)		toEdit.colorToning.saturatedOpacity			= dontforceSet && options.baBehav[ADDSET_COLORTONING_SATOPACITY] ? toEdit.colorToning.saturatedOpacity + mods.colorToning.saturatedOpacity : mods.colorToning.saturatedOpacity;
	if (colorToning.strength)				toEdit.colorToning.strength					= dontforceSet && options.baBehav[ADDSET_COLORTONING_STRENGTH] ? toEdit.colorToning.strength + mods.colorToning.strength: mods.colorToning.strength;

	if (colorToning.shadowsColSat)			toEdit.colorToning.shadowsColSat			= mods.colorToning.shadowsColSat;
	if (colorToning.hlColSat)				toEdit.colorToning.hlColSat	= mods.colorToning.hlColSat;
	if (colorToning.balance)				toEdit.colorToning.balance	= dontforceSet && options.baBehav[ADDSET_COLORTONING_BALANCE] ? toEdit.colorToning.balance + mods.colorToning.balance : mods.colorToning.balance;
	if (colorToning.clcurve)				toEdit.colorToning.clcurve	= mods.colorToning.clcurve;
	if (colorToning.method)					toEdit.colorToning.method	= mods.colorToning.method;
	if (colorToning.cl2curve)				toEdit.colorToning.cl2curve	= mods.colorToning.cl2curve;
	if (colorToning.lumamode)				toEdit.colorToning.lumamode	= mods.colorToning.lumamode;
	if (colorToning.satlow)					toEdit.colorToning.satlow	= dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit.colorToning.satlow + mods.colorToning.satlow : mods.colorToning.satlow;
	if (colorToning.sathigh)				toEdit.colorToning.sathigh	= dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit.colorToning.sathigh + mods.colorToning.sathigh : mods.colorToning.sathigh;
	if (colorToning.redlow)					toEdit.colorToning.redlow	= dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit.colorToning.redlow + mods.colorToning.redlow : mods.colorToning.redlow;
	if (colorToning.greenlow)				toEdit.colorToning.greenlow	= dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit.colorToning.greenlow + mods.colorToning.greenlow : mods.colorToning.greenlow;
	if (colorToning.bluelow)				toEdit.colorToning.bluelow	= dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit.colorToning.bluelow + mods.colorToning.bluelow : mods.colorToning.bluelow;	
	if (colorToning.redmed)					toEdit.colorToning.redmed	= dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit.colorToning.redmed + mods.colorToning.redmed : mods.colorToning.redmed;
	if (colorToning.greenmed)				toEdit.colorToning.greenmed	= dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit.colorToning.greenmed + mods.colorToning.greenmed : mods.colorToning.greenmed;
	if (colorToning.bluemed)				toEdit.colorToning.bluemed	= dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit.colorToning.bluemed + mods.colorToning.bluemed : mods.colorToning.bluemed;
	if (colorToning.redhigh)				toEdit.colorToning.redhigh	= dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit.colorToning.redhigh + mods.colorToning.redhigh : mods.colorToning.redhigh;
	if (colorToning.greenhigh)				toEdit.colorToning.greenhigh= dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit.colorToning.greenhigh + mods.colorToning.greenhigh : mods.colorToning.greenhigh;
	if (colorToning.bluehigh)				toEdit.colorToning.bluehigh	= dontforceSet && options.baBehav[ADDSET_COLORTONING_SPLIT] ? toEdit.colorToning.bluehigh + mods.colorToning.bluehigh : mods.colorToning.bluehigh;

	if (sharpenEdge.enabled)				toEdit.sharpenEdge.enabled 	= mods.sharpenEdge.enabled;
	if (sharpenEdge.passes)					toEdit.sharpenEdge.passes	= dontforceSet && options.baBehav[ADDSET_SHARPENEDGE_PASS] ? toEdit.sharpenEdge.passes + mods.sharpenEdge.passes : mods.sharpenEdge.passes;
	if (sharpenEdge.amount)					toEdit.sharpenEdge.amount	= dontforceSet && options.baBehav[ADDSET_SHARPENEDGE_AMOUNT] ? toEdit.sharpenEdge.amount + mods.sharpenEdge.amount : mods.sharpenEdge.amount;
	if (sharpenEdge.threechannels)			toEdit.sharpenEdge.threechannels 	= mods.sharpenEdge.threechannels;
	if (sharpenMicro.enabled)				toEdit.sharpenMicro.enabled 	= mods.sharpenMicro.enabled;
	if (sharpenMicro.matrix)				toEdit.sharpenMicro.matrix	= mods.sharpenMicro.matrix;
	if (sharpenMicro.amount)				toEdit.sharpenMicro.amount	= dontforceSet && options.baBehav[ADDSET_SHARPENMICRO_AMOUNT] ? toEdit.sharpenMicro.amount + mods.sharpenMicro.amount : mods.sharpenMicro.amount;
	if (sharpenMicro.uniformity)			toEdit.sharpenMicro.uniformity	= dontforceSet && options.baBehav[ADDSET_SHARPENMICRO_UNIFORMITY] ? toEdit.sharpenMicro.uniformity + mods.sharpenMicro.uniformity : mods.sharpenMicro.uniformity;
	if (sharpening.enabled)					toEdit.sharpening.enabled 	= mods.sharpening.enabled;
	if (sharpening.radius)					toEdit.sharpening.radius 	= mods.sharpening.radius;
	if (sharpening.amount)					toEdit.sharpening.amount 	= dontforceSet && options.baBehav[ADDSET_SHARP_AMOUNT] ? toEdit.sharpening.amount + mods.sharpening.amount : mods.sharpening.amount;
	if (sharpening.threshold)				toEdit.sharpening.threshold = mods.sharpening.threshold;
	if (sharpening.edgesonly)				toEdit.sharpening.edgesonly 	= mods.sharpening.edgesonly;
	if (sharpening.edges_radius)			toEdit.sharpening.edges_radius 	= mods.sharpening.edges_radius;
	if (sharpening.edges_tolerance)			toEdit.sharpening.edges_tolerance	 = mods.sharpening.edges_tolerance;
	if (sharpening.halocontrol)				toEdit.sharpening.halocontrol 		 = mods.sharpening.halocontrol;
	if (sharpening.halocontrol_amount)		toEdit.sharpening.halocontrol_amount = mods.sharpening.halocontrol_amount;
	if (sharpening.method)					toEdit.sharpening.method 		= mods.sharpening.method;
	if (sharpening.deconvamount)			toEdit.sharpening.deconvamount 	= dontforceSet && options.baBehav[ADDSET_SHARP_AMOUNT] ? toEdit.sharpening.deconvamount + mods.sharpening.deconvamount : mods.sharpening.deconvamount;
	if (sharpening.deconvradius)			toEdit.sharpening.deconvradius 	= mods.sharpening.deconvradius;
	if (sharpening.deconviter)				toEdit.sharpening.deconviter 	= mods.sharpening.deconviter;
	if (sharpening.deconvdamping)			toEdit.sharpening.deconvdamping = mods.sharpening.deconvdamping;
	if (vibrance.enabled)					toEdit.vibrance.enabled			= mods.vibrance.enabled;
	if (vibrance.pastels)					toEdit.vibrance.pastels			= dontforceSet && options.baBehav[ADDSET_VIBRANCE_PASTELS] ? toEdit.vibrance.pastels + mods.vibrance.pastels : mods.vibrance.pastels;
	if (vibrance.saturated)					toEdit.vibrance.saturated		= dontforceSet && options.baBehav[ADDSET_VIBRANCE_SATURATED] ? toEdit.vibrance.saturated + mods.vibrance.saturated : mods.vibrance.saturated;
	if (vibrance.psthreshold)				toEdit.vibrance.psthreshold		= mods.vibrance.psthreshold;
	if (vibrance.protectskins)				toEdit.vibrance.protectskins	= mods.vibrance.protectskins;
	if (vibrance.avoidcolorshift)			toEdit.vibrance.avoidcolorshift	= mods.vibrance.avoidcolorshift;
	if (vibrance.pastsattog)				toEdit.vibrance.pastsattog	    = mods.vibrance.pastsattog;
	if (vibrance.skintonescurve)			toEdit.vibrance.skintonescurve	= mods.vibrance.skintonescurve;

	//if (colorBoost.amount)					toEdit.colorBoost.amount		= dontforceSet && options.baBehav[ADDSET_CBOOST_AMOUNT] ? toEdit.colorBoost.amount + mods.colorBoost.amount : mods.colorBoost.amount;
	//if (colorBoost.avoidclip)				toEdit.colorBoost.avoidclip 	= mods.colorBoost.avoidclip;
	//if (colorBoost.enable_saturationlimiter)toEdit.colorBoost.enable_saturationlimiter 	= mods.colorBoost.enable_saturationlimiter;
	//if (colorBoost.saturationlimit)			toEdit.colorBoost.saturationlimit 	= mods.colorBoost.saturationlimit;
	if (wb.method)							toEdit.wb.method 	= mods.wb.method;
	if (wb.equal)							toEdit.wb.equal 	= dontforceSet && options.baBehav[ADDSET_WB_EQUAL] ? toEdit.wb.equal + mods.wb.equal : mods.wb.equal;
	if (wb.green)							toEdit.wb.green 	= dontforceSet && options.baBehav[ADDSET_WB_GREEN] ? toEdit.wb.green + mods.wb.green : mods.wb.green;
	if (wb.temperature)						toEdit.wb.temperature 	= dontforceSet && options.baBehav[ADDSET_WB_TEMPERATURE] ? toEdit.wb.temperature + mods.wb.temperature : mods.wb.temperature;
	//if (colorShift.a)						toEdit.colorShift.a 	= dontforceSet && options.baBehav[ADDSET_CS_BLUEYELLOW] ? toEdit.colorShift.a + mods.colorShift.a : mods.colorShift.a;
	//if (colorShift.b)						toEdit.colorShift.b 	= dontforceSet && options.baBehav[ADDSET_CS_GREENMAGENTA] ? toEdit.colorShift.b + mods.colorShift.b : mods.colorShift.b;
	//if (lumaDenoise.enabled)				toEdit.lumaDenoise.enabled 	= mods.lumaDenoise.enabled;
	//if (lumaDenoise.radius)					toEdit.lumaDenoise.radius 	= mods.lumaDenoise.radius;
	//if (lumaDenoise.edgetolerance)			toEdit.lumaDenoise.edgetolerance 	= dontforceSet && options.baBehav[ADDSET_LD_EDGETOLERANCE] ? toEdit.lumaDenoise.edgetolerance + mods.lumaDenoise.edgetolerance : mods.lumaDenoise.edgetolerance;
	//if (colorDenoise.enabled)				toEdit.colorDenoise.enabled 	= mods.colorDenoise.enabled;
	//if (colorDenoise.amount)				toEdit.colorDenoise.amount 	= mods.colorDenoise.amount;
	
	if (defringe.enabled)					toEdit.defringe.enabled   = mods.defringe.enabled;
	if (defringe.radius)					toEdit.defringe.radius    = mods.defringe.radius;
	if (defringe.threshold)					toEdit.defringe.threshold = mods.defringe.threshold;
	if (defringe.huecurve)					toEdit.defringe.huecurve  = mods.defringe.huecurve;
	
	if (colorappearance.curve)				toEdit.colorappearance.curve      = mods.colorappearance.curve;
	if (colorappearance.curve2)				toEdit.colorappearance.curve2     = mods.colorappearance.curve2;
	if (colorappearance.curve3)				toEdit.colorappearance.curve3     = mods.colorappearance.curve3;
	if (colorappearance.curveMode)			toEdit.colorappearance.curveMode  = mods.colorappearance.curveMode;
	if (colorappearance.curveMode2)			toEdit.colorappearance.curveMode2 = mods.colorappearance.curveMode2;
	if (colorappearance.curveMode3)			toEdit.colorappearance.curveMode3 = mods.colorappearance.curveMode3;

	if (colorappearance.enabled)			toEdit.colorappearance.enabled		= mods.colorappearance.enabled;
	if (colorappearance.degree)				toEdit.colorappearance.degree		= dontforceSet && options.baBehav[ADDSET_CAT_DEGREE] ? toEdit.colorappearance.degree + mods.colorappearance.degree : mods.colorappearance.degree;
	if (colorappearance.autodegree)			toEdit.colorappearance.autodegree	= mods.colorappearance.autodegree;
	if (colorappearance.surround)			toEdit.colorappearance.surround		= mods.colorappearance.surround;
	if (colorappearance.autoadapscen)		toEdit.colorappearance.autoadapscen	= mods.colorappearance.autoadapscen;
	if (colorappearance.adapscen)			toEdit.colorappearance.adapscen	= mods.colorappearance.adapscen;
	if (colorappearance.adaplum)			toEdit.colorappearance.adaplum		= dontforceSet && options.baBehav[ADDSET_CAT_ADAPTVIEWING] ? toEdit.colorappearance.adaplum + mods.colorappearance.adaplum : mods.colorappearance.adaplum;
	if (colorappearance.badpixsl)			toEdit.colorappearance.badpixsl		= dontforceSet && options.baBehav[ADDSET_CAT_BADPIX] ? toEdit.colorappearance.badpixsl + mods.colorappearance.badpixsl : mods.colorappearance.badpixsl;
	if (colorappearance.wbmodel)			toEdit.colorappearance.wbmodel		= mods.colorappearance.wbmodel;
	if (colorappearance.algo)				toEdit.colorappearance.algo		= mods.colorappearance.algo;

	if (colorappearance.jlight)				toEdit.colorappearance.jlight		= dontforceSet && options.baBehav[ADDSET_CAT_LIGHT] ? toEdit.colorappearance.jlight + mods.colorappearance.jlight : mods.colorappearance.jlight;
	if (colorappearance.qbright)			toEdit.colorappearance.qbright		= dontforceSet && options.baBehav[ADDSET_CAT_BRIGHT] ? toEdit.colorappearance.qbright + mods.colorappearance.qbright : mods.colorappearance.qbright;
	if (colorappearance.chroma)				toEdit.colorappearance.chroma		= dontforceSet && options.baBehav[ADDSET_CAT_CHROMA] ? toEdit.colorappearance.chroma + mods.colorappearance.chroma : mods.colorappearance.chroma;
	if (colorappearance.schroma)			toEdit.colorappearance.schroma		= dontforceSet && options.baBehav[ADDSET_CAT_CHROMA_S] ? toEdit.colorappearance.schroma + mods.colorappearance.schroma : mods.colorappearance.schroma;
	if (colorappearance.mchroma)			toEdit.colorappearance.mchroma		= dontforceSet && options.baBehav[ADDSET_CAT_CHROMA_M] ? toEdit.colorappearance.mchroma + mods.colorappearance.mchroma : mods.colorappearance.mchroma;
	if (colorappearance.contrast)			toEdit.colorappearance.contrast		= dontforceSet && options.baBehav[ADDSET_CAT_CONTRAST] ? toEdit.colorappearance.contrast + mods.colorappearance.contrast : mods.colorappearance.contrast;
	if (colorappearance.qcontrast)			toEdit.colorappearance.qcontrast	= dontforceSet && options.baBehav[ADDSET_CAT_CONTRAST_Q] ? toEdit.colorappearance.qcontrast + mods.colorappearance.qcontrast : mods.colorappearance.qcontrast;
	if (colorappearance.colorh)				toEdit.colorappearance.colorh		= dontforceSet && options.baBehav[ADDSET_CAT_HUE] ? toEdit.colorappearance.colorh + mods.colorappearance.colorh : mods.colorappearance.colorh;
	if (colorappearance.rstprotection)		toEdit.colorappearance.rstprotection= dontforceSet && options.baBehav[ADDSET_CAT_RSTPRO] ? toEdit.colorappearance.rstprotection + mods.colorappearance.rstprotection : mods.colorappearance.rstprotection;
 	if (colorappearance.surrsource)			toEdit.colorappearance.surrsource = mods.colorappearance.surrsource;
 	if (colorappearance.gamut)				toEdit.colorappearance.gamut = mods.colorappearance.gamut;
// 	if (colorappearance.badpix)				toEdit.colorappearance.badpix = mods.colorappearance.badpix;
 	if (colorappearance.datacie)			toEdit.colorappearance.datacie = mods.colorappearance.datacie;
 	if (colorappearance.tonecie)			toEdit.colorappearance.tonecie = mods.colorappearance.tonecie;
// 	if (colorappearance.sharpcie)			toEdit.colorappearance.sharpcie = mods.colorappearance.sharpcie;
	if (impulseDenoise.enabled)				toEdit.impulseDenoise.enabled 	= mods.impulseDenoise.enabled;
	if (impulseDenoise.thresh)				toEdit.impulseDenoise.thresh 	= mods.impulseDenoise.thresh;

	if (dirpyrDenoise.enabled)				toEdit.dirpyrDenoise.enabled 	= mods.dirpyrDenoise.enabled;
	if (dirpyrDenoise.enhance)				toEdit.dirpyrDenoise.enhance 	= mods.dirpyrDenoise.enhance;
	if (dirpyrDenoise.median)				toEdit.dirpyrDenoise.median 	= mods.dirpyrDenoise.median;
	if (dirpyrDenoise.autochroma)			toEdit.dirpyrDenoise.autochroma 	= mods.dirpyrDenoise.autochroma;
	if (dirpyrDenoise.luma)					toEdit.dirpyrDenoise.luma		= dontforceSet && options.baBehav[ADDSET_DIRPYRDN_LUMA] ? toEdit.dirpyrDenoise.luma + mods.dirpyrDenoise.luma : mods.dirpyrDenoise.luma;
	if (dirpyrDenoise.lcurve)				toEdit.dirpyrDenoise.lcurve 	    = mods.dirpyrDenoise.lcurve;
	if (dirpyrDenoise.cccurve)				toEdit.dirpyrDenoise.cccurve 	    = mods.dirpyrDenoise.cccurve;
	if (dirpyrDenoise.Ldetail)				toEdit.dirpyrDenoise.Ldetail	= dontforceSet && options.baBehav[ADDSET_DIRPYRDN_LUMDET] ? toEdit.dirpyrDenoise.Ldetail + mods.dirpyrDenoise.Ldetail : mods.dirpyrDenoise.Ldetail;
	if (dirpyrDenoise.chroma)				toEdit.dirpyrDenoise.chroma		= dontforceSet && options.baBehav[ADDSET_DIRPYRDN_CHROMA] ? toEdit.dirpyrDenoise.chroma + mods.dirpyrDenoise.chroma : mods.dirpyrDenoise.chroma;
	if (dirpyrDenoise.redchro)				toEdit.dirpyrDenoise.redchro	= dontforceSet && options.baBehav[ADDSET_DIRPYRDN_CHROMARED] ? toEdit.dirpyrDenoise.redchro + mods.dirpyrDenoise.redchro : mods.dirpyrDenoise.redchro;
	if (dirpyrDenoise.bluechro)				toEdit.dirpyrDenoise.bluechro	= dontforceSet && options.baBehav[ADDSET_DIRPYRDN_CHROMABLUE] ? toEdit.dirpyrDenoise.bluechro + mods.dirpyrDenoise.bluechro : mods.dirpyrDenoise.bluechro;
	if (dirpyrDenoise.gamma)				toEdit.dirpyrDenoise.gamma		= dontforceSet && options.baBehav[ADDSET_DIRPYRDN_GAMMA] ? toEdit.dirpyrDenoise.gamma + mods.dirpyrDenoise.gamma : mods.dirpyrDenoise.gamma;
	if (dirpyrDenoise.passes)				toEdit.dirpyrDenoise.passes		= dontforceSet && options.baBehav[ADDSET_DIRPYRDN_PASSES] ? toEdit.dirpyrDenoise.passes + mods.dirpyrDenoise.passes : mods.dirpyrDenoise.passes;
//	if (dirpyrDenoise.perform)				toEdit.dirpyrDenoise.perform 	= mods.dirpyrDenoise.perform;
	if (dirpyrDenoise.dmethod)				toEdit.dirpyrDenoise.dmethod		= mods.dirpyrDenoise.dmethod;
	if (dirpyrDenoise.Lmethod)				toEdit.dirpyrDenoise.Lmethod		= mods.dirpyrDenoise.Lmethod;
	if (dirpyrDenoise.Cmethod)				toEdit.dirpyrDenoise.Cmethod		= mods.dirpyrDenoise.Cmethod;
	if (dirpyrDenoise.C2method)				toEdit.dirpyrDenoise.C2method		= mods.dirpyrDenoise.C2method;
	if (dirpyrDenoise.smethod)				toEdit.dirpyrDenoise.smethod		= mods.dirpyrDenoise.smethod;
	if (dirpyrDenoise.medmethod)			toEdit.dirpyrDenoise.medmethod		= mods.dirpyrDenoise.medmethod;
	if (dirpyrDenoise.methodmed)			toEdit.dirpyrDenoise.methodmed		= mods.dirpyrDenoise.methodmed;
	if (dirpyrDenoise.rgbmethod)			toEdit.dirpyrDenoise.rgbmethod		= mods.dirpyrDenoise.rgbmethod;

	if (epd.enabled)						toEdit.epd.enabled				= mods.epd.enabled;
	if (epd.strength)						toEdit.epd.strength				= mods.epd.strength;
	if (epd.edgeStopping)					toEdit.epd.edgeStopping			= mods.epd.edgeStopping;
	if (epd.scale)							toEdit.epd.scale				= mods.epd.scale;
	if (epd.reweightingIterates)			toEdit.epd.reweightingIterates	= mods.epd.reweightingIterates;

	if (sh.enabled)		    				toEdit.sh.enabled 	    = mods.sh.enabled;
	if (sh.hq)		        				toEdit.sh.hq     	    = mods.sh.hq;
	if (sh.highlights)						toEdit.sh.highlights 	= dontforceSet && options.baBehav[ADDSET_SH_HIGHLIGHTS] ? toEdit.sh.highlights + mods.sh.highlights : mods.sh.highlights;
	if (sh.htonalwidth)						toEdit.sh.htonalwidth 	= mods.sh.htonalwidth;
	if (sh.shadows)		    				toEdit.sh.shadows 	    = dontforceSet && options.baBehav[ADDSET_SH_SHADOWS] ? toEdit.sh.shadows + mods.sh.shadows : mods.sh.shadows;
	if (sh.stonalwidth)						toEdit.sh.stonalwidth 	= mods.sh.stonalwidth;
	if (sh.localcontrast)					toEdit.sh.localcontrast = dontforceSet && options.baBehav[ADDSET_SH_LOCALCONTRAST] ? toEdit.sh.localcontrast + mods.sh.localcontrast : mods.sh.localcontrast;
	if (sh.radius)		    				toEdit.sh.radius 	    = mods.sh.radius;
	if (crop.enabled)						toEdit.crop.enabled = mods.crop.enabled;
	if (crop.x)		        				toEdit.crop.x 	    = mods.crop.x;
	if (crop.y)		        				toEdit.crop.y 	    = mods.crop.y;
	if (crop.w)		        				toEdit.crop.w 	    = mods.crop.w;
	if (crop.h)		        				toEdit.crop.h 	    = mods.crop.h;
	if (crop.fixratio)						toEdit.crop.fixratio 	= mods.crop.fixratio;
	if (crop.ratio)		    				toEdit.crop.ratio 	    = mods.crop.ratio;
	if (crop.orientation)					toEdit.crop.orientation = mods.crop.orientation;
	if (crop.guide)		    				toEdit.crop.guide 	    = mods.crop.guide;
	if (coarse.rotate)						toEdit.coarse.rotate 	= mods.coarse.rotate;
	if (coarse.hflip)						toEdit.coarse.hflip 	= mods.coarse.hflip;
	if (coarse.vflip)						toEdit.coarse.vflip 	= mods.coarse.vflip;
	if (commonTrans.autofill)				toEdit.commonTrans.autofill		= mods.commonTrans.autofill;
	if (rotate.degree)						toEdit.rotate.degree 			= dontforceSet && options.baBehav[ADDSET_ROTATE_DEGREE] ? toEdit.rotate.degree + mods.rotate.degree : mods.rotate.degree;
	if (distortion.amount)					toEdit.distortion.amount 		= dontforceSet && options.baBehav[ADDSET_DIST_AMOUNT] ? toEdit.distortion.amount + mods.distortion.amount : mods.distortion.amount;
	if (lensProf.lcpFile)                   toEdit.lensProf.lcpFile         = mods.lensProf.lcpFile;
    if (lensProf.useDist)                   toEdit.lensProf.useDist         = mods.lensProf.useDist;
    if (lensProf.useVign)                   toEdit.lensProf.useVign         = mods.lensProf.useVign;
    if (lensProf.useCA)                     toEdit.lensProf.useCA           = mods.lensProf.useCA;

	if (perspective.horizontal)				toEdit.perspective.horizontal 	= dontforceSet && options.baBehav[ADDSET_PERSPECTIVE] ? toEdit.perspective.horizontal + mods.perspective.horizontal : mods.perspective.horizontal;
	if (perspective.vertical)				toEdit.perspective.vertical 	= dontforceSet && options.baBehav[ADDSET_PERSPECTIVE] ? toEdit.perspective.vertical + mods.perspective.vertical : mods.perspective.vertical;
	if (gradient.enabled)					toEdit.gradient.enabled 	= mods.gradient.enabled;
	if (gradient.degree)					toEdit.gradient.degree		= dontforceSet && options.baBehav[ADDSET_GRADIENT_DEGREE] ? toEdit.gradient.degree + mods.gradient.degree : mods.gradient.degree;
	if (gradient.feather)					toEdit.gradient.feather 	= mods.gradient.feather;
	if (gradient.strength)					toEdit.gradient.strength 	= mods.gradient.strength;
	if (gradient.centerX)					toEdit.gradient.centerX 	= mods.gradient.centerX;
	if (gradient.centerY)					toEdit.gradient.centerY 	= mods.gradient.centerY;
	if (pcvignette.enabled)					toEdit.pcvignette.enabled 	= mods.pcvignette.enabled;
	if (pcvignette.strength)				toEdit.pcvignette.strength 	= mods.pcvignette.strength;
	if (pcvignette.feather)					toEdit.pcvignette.feather 	= mods.pcvignette.feather;
	if (pcvignette.roundness)				toEdit.pcvignette.roundness = mods.pcvignette.roundness;
	if (cacorrection.red)					toEdit.cacorrection.red 	= dontforceSet && options.baBehav[ADDSET_CA] ? toEdit.cacorrection.red + mods.cacorrection.red : mods.cacorrection.red;
	if (cacorrection.blue)					toEdit.cacorrection.blue 	= dontforceSet && options.baBehav[ADDSET_CA] ? toEdit.cacorrection.blue + mods.cacorrection.blue : mods.cacorrection.blue;
	if (vignetting.amount)					toEdit.vignetting.amount 	= dontforceSet && options.baBehav[ADDSET_VIGN_AMOUNT] ? toEdit.vignetting.amount + mods.vignetting.amount : mods.vignetting.amount;
	if (vignetting.radius)					toEdit.vignetting.radius 	= mods.vignetting.radius;
	if (vignetting.strength)				toEdit.vignetting.strength 	= mods.vignetting.strength;
	if (vignetting.centerX)					toEdit.vignetting.centerX 	= mods.vignetting.centerX;
	if (vignetting.centerY)					toEdit.vignetting.centerY 	= mods.vignetting.centerY;
	for (int i=0; i<3; i++) {
		if (chmixer.red[i])		toEdit.chmixer.red[i] 	= dontforceSet && options.baBehav[ADDSET_CHMIXER] ? toEdit.chmixer.red[i] + mods.chmixer.red[i] : mods.chmixer.red[i];
		if (chmixer.green[i])	toEdit.chmixer.green[i]	= dontforceSet && options.baBehav[ADDSET_CHMIXER] ? toEdit.chmixer.green[i] + mods.chmixer.green[i] : mods.chmixer.green[i];
		if (chmixer.blue[i])	toEdit.chmixer.blue[i] 	= dontforceSet && options.baBehav[ADDSET_CHMIXER] ? toEdit.chmixer.blue[i] + mods.chmixer.blue[i] : mods.chmixer.blue[i];
	}
	if (blackwhite.enabled)			toEdit.blackwhite.enabled			= mods.blackwhite.enabled;
	if (blackwhite.method)			toEdit.blackwhite.method			= mods.blackwhite.method;
	if (blackwhite.luminanceCurve)	toEdit.blackwhite.luminanceCurve	= mods.blackwhite.luminanceCurve;
	if (blackwhite.autoc)			toEdit.blackwhite.autoc				= mods.blackwhite.autoc;
	if (blackwhite.setting)			toEdit.blackwhite.setting			= mods.blackwhite.setting;
	if (blackwhite.enabledcc)		toEdit.blackwhite.enabledcc			= mods.blackwhite.enabledcc;
	if (blackwhite.filter)			toEdit.blackwhite.filter			= mods.blackwhite.filter;
	if (blackwhite.mixerRed)		toEdit.blackwhite.mixerRed 			= dontforceSet && options.baBehav[ADDSET_BLACKWHITE_HUES] ? toEdit.blackwhite.mixerRed + mods.blackwhite.mixerRed : mods.blackwhite.mixerRed;
	if (blackwhite.mixerOrange)		toEdit.blackwhite.mixerOrange 		= dontforceSet && options.baBehav[ADDSET_BLACKWHITE_HUES] ? toEdit.blackwhite.mixerOrange + mods.blackwhite.mixerOrange : mods.blackwhite.mixerOrange;
	if (blackwhite.mixerYellow)		toEdit.blackwhite.mixerYellow 		= dontforceSet && options.baBehav[ADDSET_BLACKWHITE_HUES] ? toEdit.blackwhite.mixerYellow + mods.blackwhite.mixerYellow : mods.blackwhite.mixerYellow;
	if (blackwhite.mixerGreen)		toEdit.blackwhite.mixerGreen 		= dontforceSet && options.baBehav[ADDSET_BLACKWHITE_HUES] ? toEdit.blackwhite.mixerGreen + mods.blackwhite.mixerGreen : mods.blackwhite.mixerGreen;
	if (blackwhite.mixerCyan)		toEdit.blackwhite.mixerCyan 		= dontforceSet && options.baBehav[ADDSET_BLACKWHITE_HUES] ? toEdit.blackwhite.mixerCyan + mods.blackwhite.mixerCyan : mods.blackwhite.mixerCyan;
	if (blackwhite.mixerBlue)		toEdit.blackwhite.mixerBlue 		= dontforceSet && options.baBehav[ADDSET_BLACKWHITE_HUES] ? toEdit.blackwhite.mixerBlue + mods.blackwhite.mixerBlue : mods.blackwhite.mixerBlue;
	if (blackwhite.mixerMagenta)	toEdit.blackwhite.mixerMagenta 		= dontforceSet && options.baBehav[ADDSET_BLACKWHITE_HUES] ? toEdit.blackwhite.mixerMagenta + mods.blackwhite.mixerMagenta : mods.blackwhite.mixerMagenta;
	if (blackwhite.mixerPurple)		toEdit.blackwhite.mixerPurple 		= dontforceSet && options.baBehav[ADDSET_BLACKWHITE_HUES] ? toEdit.blackwhite.mixerPurple + mods.blackwhite.mixerPurple : mods.blackwhite.mixerPurple;
	if (blackwhite.gammaRed)		toEdit.blackwhite.gammaRed 			= dontforceSet && options.baBehav[ADDSET_BLACKWHITE_GAMMA] ? toEdit.blackwhite.gammaRed + mods.blackwhite.gammaRed : mods.blackwhite.gammaRed;
	if (blackwhite.gammaGreen)		toEdit.blackwhite.gammaGreen 		= dontforceSet && options.baBehav[ADDSET_BLACKWHITE_GAMMA] ? toEdit.blackwhite.gammaGreen + mods.blackwhite.gammaGreen : mods.blackwhite.gammaGreen;
	if (blackwhite.gammaBlue)		toEdit.blackwhite.gammaBlue 		= dontforceSet && options.baBehav[ADDSET_BLACKWHITE_GAMMA] ? toEdit.blackwhite.gammaBlue + mods.blackwhite.gammaBlue : mods.blackwhite.gammaBlue;
	if (blackwhite.beforeCurve)		toEdit.blackwhite.beforeCurve		= mods.blackwhite.beforeCurve;
	if (blackwhite.beforeCurveMode)	toEdit.blackwhite.beforeCurveMode	= mods.blackwhite.beforeCurveMode;
	if (blackwhite.afterCurve)		toEdit.blackwhite.afterCurve		= mods.blackwhite.afterCurve;
	if (blackwhite.afterCurveMode)	toEdit.blackwhite.afterCurveMode	= mods.blackwhite.afterCurveMode;
	if (blackwhite.algo)			toEdit.blackwhite.algo				= mods.blackwhite.algo;
	
	if (resize.scale)		toEdit.resize.scale 	= mods.resize.scale;
	if (resize.appliesTo)	toEdit.resize.appliesTo = mods.resize.appliesTo;
	if (resize.method)		toEdit.resize.method 	= mods.resize.method;
	if (resize.dataspec)	toEdit.resize.dataspec 	= mods.resize.dataspec;
	if (resize.width)	    toEdit.resize.width 	= mods.resize.width;
	if (resize.height)	    toEdit.resize.height 	= mods.resize.height;
	if (resize.enabled)	    toEdit.resize.enabled 	= mods.resize.enabled;
	if (icm.input)		    toEdit.icm.input 	    = mods.icm.input;
    if (icm.toneCurve)      toEdit.icm.toneCurve = mods.icm.toneCurve;
    if (icm.blendCMSMatrix)	toEdit.icm.blendCMSMatrix = mods.icm.blendCMSMatrix;
    if (icm.dcpIlluminant) toEdit.icm.dcpIlluminant = mods.icm.dcpIlluminant;
	if (icm.working)		toEdit.icm.working 	    = mods.icm.working;
	if (icm.output)		    toEdit.icm.output       = mods.icm.output;
	//if (icm.gampos)		    toEdit.icm.gampos       = mods.icm.gampos;
	//if (icm.slpos)		    toEdit.icm.slpos        = mods.icm.slpos;
	if (icm.gampos)			toEdit.icm.gampos		= dontforceSet && options.baBehav[ADDSET_FREE_OUPUT_GAMMA] ? toEdit.icm.gampos + mods.icm.gampos : mods.icm.gampos;
	if (icm.slpos)			toEdit.icm.slpos		= dontforceSet && options.baBehav[ADDSET_FREE_OUTPUT_SLOPE] ? toEdit.icm.slpos + mods.icm.slpos : mods.icm.slpos;	
	if (icm.gamma)		    toEdit.icm.gamma        = mods.icm.gamma;
	if (icm.freegamma)		toEdit.icm.freegamma    = mods.icm.freegamma;
	if (raw.bayersensor.method)          toEdit.raw.bayersensor.method           = mods.raw.bayersensor.method;
	if (raw.bayersensor.ccSteps)         toEdit.raw.bayersensor.ccSteps          = mods.raw.bayersensor.ccSteps;
	if (raw.bayersensor.exBlack0)        toEdit.raw.bayersensor.black0           = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit.raw.bayersensor.black0 + mods.raw.bayersensor.black0 : mods.raw.bayersensor.black0;
	if (raw.bayersensor.exBlack1)        toEdit.raw.bayersensor.black1           = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit.raw.bayersensor.black1 + mods.raw.bayersensor.black1 : mods.raw.bayersensor.black1;
	if (raw.bayersensor.exBlack2)        toEdit.raw.bayersensor.black2           = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit.raw.bayersensor.black2 + mods.raw.bayersensor.black2 : mods.raw.bayersensor.black2;
	if (raw.bayersensor.exBlack3)        toEdit.raw.bayersensor.black3           = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit.raw.bayersensor.black3 + mods.raw.bayersensor.black3 : mods.raw.bayersensor.black3;
	if (raw.bayersensor.exTwoGreen)      toEdit.raw.bayersensor.twogreen         = mods.raw.bayersensor.twogreen;
	if (raw.bayersensor.dcbIterations)   toEdit.raw.bayersensor.dcb_iterations   = mods.raw.bayersensor.dcb_iterations;
	if (raw.bayersensor.dcbEnhance)      toEdit.raw.bayersensor.dcb_enhance      = mods.raw.bayersensor.dcb_enhance;
	if (raw.bayersensor.lmmseIterations) toEdit.raw.bayersensor.lmmse_iterations = mods.raw.bayersensor.lmmse_iterations;
	//if (raw.bayersensor.allEnhance)    toEdit.raw.bayersensor.all_enhance      = mods.raw.bayersensor.all_enhance;
	if (raw.bayersensor.greenEq)         toEdit.raw.bayersensor.greenthresh      = dontforceSet && options.baBehav[ADDSET_PREPROCESS_GREENEQUIL] ? toEdit.raw.bayersensor.greenthresh + mods.raw.bayersensor.greenthresh : mods.raw.bayersensor.greenthresh;
	if (raw.bayersensor.linenoise)       toEdit.raw.bayersensor.linenoise        = dontforceSet && options.baBehav[ADDSET_PREPROCESS_LINEDENOISE] ? toEdit.raw.bayersensor.linenoise + mods.raw.bayersensor.linenoise : mods.raw.bayersensor.linenoise;

	if (raw.xtranssensor.method)         toEdit.raw.xtranssensor.method          = mods.raw.xtranssensor.method;
	if (raw.xtranssensor.ccSteps)        toEdit.raw.xtranssensor.ccSteps         = mods.raw.xtranssensor.ccSteps;
	if (raw.xtranssensor.exBlackRed)     toEdit.raw.xtranssensor.blackred        = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit.raw.xtranssensor.blackred + mods.raw.xtranssensor.blackred : mods.raw.xtranssensor.blackred;
	if (raw.xtranssensor.exBlackGreen)   toEdit.raw.xtranssensor.blackgreen      = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit.raw.xtranssensor.blackgreen + mods.raw.xtranssensor.blackgreen : mods.raw.xtranssensor.blackgreen;
	if (raw.xtranssensor.exBlackBlue)    toEdit.raw.xtranssensor.blackblue       = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit.raw.xtranssensor.blackblue + mods.raw.xtranssensor.blackblue : mods.raw.xtranssensor.blackblue;

	if (raw.caCorrection)   toEdit.raw.ca_autocorrect  = mods.raw.ca_autocorrect;
	if (raw.caRed)          toEdit.raw.cared           = dontforceSet && options.baBehav[ADDSET_RAWCACORR] ? toEdit.raw.cared + mods.raw.cared : mods.raw.cared;
	if (raw.caBlue)         toEdit.raw.cablue          = dontforceSet && options.baBehav[ADDSET_RAWCACORR] ? toEdit.raw.cablue + mods.raw.cablue : mods.raw.cablue;
	if (raw.exPos)          toEdit.raw.expos           = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_LINEAR] ? toEdit.raw.expos + mods.raw.expos : mods.raw.expos;
	if (raw.exPreser)       toEdit.raw.preser          = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_PRESER] ? toEdit.raw.preser + mods.raw.preser : mods.raw.preser;

	if (raw.hotPixelFilter)     toEdit.raw.hotPixelFilter    = mods.raw.hotPixelFilter;
	if (raw.deadPixelFilter)    toEdit.raw.deadPixelFilter   = mods.raw.deadPixelFilter;
	if (raw.hotDeadPixelThresh) toEdit.raw.hotdeadpix_thresh = mods.raw.hotdeadpix_thresh;
	if (raw.darkFrame)          toEdit.raw.dark_frame        = mods.raw.dark_frame;
	if (raw.dfAuto)             toEdit.raw.df_autoselect     = mods.raw.df_autoselect;

	if (raw.ff_file)            toEdit.raw.ff_file            = mods.raw.ff_file;
	if (raw.ff_AutoSelect)      toEdit.raw.ff_AutoSelect      = mods.raw.ff_AutoSelect;
	if (raw.ff_BlurRadius)      toEdit.raw.ff_BlurRadius      = mods.raw.ff_BlurRadius;
	if (raw.ff_BlurType)        toEdit.raw.ff_BlurType        = mods.raw.ff_BlurType;
	if (raw.ff_AutoClipControl) toEdit.raw.ff_AutoClipControl = mods.raw.ff_AutoClipControl;
	if (raw.ff_clipControl)     toEdit.raw.ff_clipControl     = dontforceSet && options.baBehav[ADDSET_RAWFFCLIPCONTROL] ? toEdit.raw.ff_clipControl + mods.raw.ff_clipControl : mods.raw.ff_clipControl;
	if (wavelet.enabled)	toEdit.wavelet.enabled   = mods.wavelet.enabled;
	if (wavelet.median)	toEdit.wavelet.median   = mods.wavelet.median;
	if (wavelet.avoid)	toEdit.wavelet.avoid   = mods.wavelet.avoid;
	if (wavelet.Lmethod)		toEdit.wavelet.Lmethod		= mods.wavelet.Lmethod;
	if (wavelet.CLmethod)		toEdit.wavelet.CLmethod		= mods.wavelet.CLmethod;
	if (wavelet.Tilesmethod)		toEdit.wavelet.Tilesmethod		= mods.wavelet.Tilesmethod;
	if (wavelet.CHmethod)		toEdit.wavelet.CHmethod		= mods.wavelet.CHmethod;
	if (wavelet.HSmethod)		toEdit.wavelet.HSmethod		= mods.wavelet.HSmethod;
	if (wavelet.Dirmethod)		toEdit.wavelet.Dirmethod		= mods.wavelet.Dirmethod;
	if (wavelet.sky)	toEdit.wavelet.sky= dontforceSet && options.baBehav[ADDSET_WA_SKYPROTECT] ? toEdit.wavelet.sky + mods.wavelet.sky : mods.wavelet.sky;
	if (wavelet.thr)toEdit.wavelet.thr= dontforceSet && options.baBehav[ADDSET_WA_THRR] ? toEdit.wavelet.thr + mods.wavelet.thr : mods.wavelet.thr;
	if (wavelet.thrH)toEdit.wavelet.thrH= dontforceSet && options.baBehav[ADDSET_WA_THRRH] ? toEdit.wavelet.thrH + mods.wavelet.thrH : mods.wavelet.thrH;
	if (wavelet.sup)		toEdit.wavelet.sup		= mods.wavelet.sup;
	if (wavelet.hllev)	toEdit.wavelet.hllev	= mods.wavelet.hllev;
	if (wavelet.bllev)	toEdit.wavelet.bllev	= mods.wavelet.bllev;
	if (wavelet.pastlev)	toEdit.wavelet.pastlev	= mods.wavelet.pastlev;
	if (wavelet.satlev)	toEdit.wavelet.satlev	= mods.wavelet.satlev;
	if (wavelet.clvcurve)	toEdit.wavelet.clvcurve	= mods.wavelet.clvcurve;
	if (wavelet.opacityCurveRG)	toEdit.wavelet.opacityCurveRG	= mods.wavelet.opacityCurveRG;
	if (wavelet.opacityCurveBY)	toEdit.wavelet.opacityCurveBY	= mods.wavelet.opacityCurveBY;
	for(int i = 0; i < 9; i++) {
	    if(wavelet.c[i])  toEdit.wavelet.c[i] = dontforceSet && options.baBehav[ADDSET_WA] ? toEdit.wavelet.c[i] + mods.wavelet.c[i] : mods.wavelet.c[i];
	}
	if (wavelet.skinprotect)toEdit.wavelet.skinprotect= dontforceSet && options.baBehav[ADDSET_WA_SKINPROTECT] ? toEdit.wavelet.skinprotect + mods.wavelet.skinprotect : mods.wavelet.skinprotect;
	if (wavelet.hueskin)	toEdit.wavelet.hueskin	= mods.wavelet.hueskin;
	if (wavelet.hueskin2)	toEdit.wavelet.hueskin2	= mods.wavelet.hueskin2;
	if (wavelet.resconH)toEdit.wavelet.resconH= dontforceSet && options.baBehav[ADDSET_WA_RESCONH] ? toEdit.wavelet.resconH + mods.wavelet.resconH : mods.wavelet.resconH;
	if (wavelet.reschro)toEdit.wavelet.reschro= dontforceSet && options.baBehav[ADDSET_WA_RESCHRO] ? toEdit.wavelet.reschro + mods.wavelet.reschro : mods.wavelet.reschro;
	if (wavelet.rescon)toEdit.wavelet.rescon= dontforceSet && options.baBehav[ADDSET_WA_RESCON] ? toEdit.wavelet.rescon + mods.wavelet.rescon : mods.wavelet.rescon;
	if (wavelet.thres)toEdit.wavelet.thres= dontforceSet && options.baBehav[ADDSET_WA_THRES] ? toEdit.wavelet.thres + mods.wavelet.thres : mods.wavelet.thres;
	if (wavelet.threshold)toEdit.wavelet.threshold= dontforceSet && options.baBehav[ADDSET_WA_THRESHOLD] ? toEdit.wavelet.threshold + mods.wavelet.threshold : mods.wavelet.threshold;
	if (wavelet.threshold2)toEdit.wavelet.threshold2= dontforceSet && options.baBehav[ADDSET_WA_THRESHOLD2] ? toEdit.wavelet.threshold2 + mods.wavelet.threshold2 : mods.wavelet.threshold2;
	if (wavelet.chro)toEdit.wavelet.chro= dontforceSet && options.baBehav[ADDSET_WA_CHRO] ? toEdit.wavelet.chro + mods.wavelet.chro : mods.wavelet.chro;
	if (wavelet.chroma)toEdit.wavelet.chroma= dontforceSet && options.baBehav[ADDSET_WA_CHROMA] ? toEdit.wavelet.chroma + mods.wavelet.chroma : mods.wavelet.chroma;
	if (wavelet.unif)toEdit.wavelet.unif= dontforceSet && options.baBehav[ADDSET_WA_UNIF] ? toEdit.wavelet.unif + mods.wavelet.unif : mods.wavelet.unif;

	if (dirpyrequalizer.enabled)	toEdit.dirpyrequalizer.enabled	= mods.dirpyrequalizer.enabled;
	if (dirpyrequalizer.gamutlab)	toEdit.dirpyrequalizer.gamutlab	= mods.dirpyrequalizer.gamutlab;
	for(int i = 0; i < 6; i++) {
		if(dirpyrequalizer.mult[i])	toEdit.dirpyrequalizer.mult[i]	= dontforceSet && options.baBehav[ADDSET_DIRPYREQ] ? toEdit.dirpyrequalizer.mult[i] + mods.dirpyrequalizer.mult[i] : mods.dirpyrequalizer.mult[i];
	}
	if (dirpyrequalizer.threshold)  toEdit.dirpyrequalizer.threshold= dontforceSet && options.baBehav[ADDSET_DIRPYREQ_THRESHOLD] ? toEdit.dirpyrequalizer.threshold + mods.dirpyrequalizer.threshold : mods.dirpyrequalizer.threshold;
	if (dirpyrequalizer.skinprotect)toEdit.dirpyrequalizer.skinprotect= dontforceSet && options.baBehav[ADDSET_DIRPYREQ_SKINPROTECT] ? toEdit.dirpyrequalizer.skinprotect + mods.dirpyrequalizer.skinprotect : mods.dirpyrequalizer.skinprotect;
	if (dirpyrequalizer.hueskin)	toEdit.dirpyrequalizer.hueskin	= mods.dirpyrequalizer.hueskin;
//	if (dirpyrequalizer.algo)		toEdit.dirpyrequalizer.algo		= mods.dirpyrequalizer.algo;
	if (hsvequalizer.hcurve)		toEdit.hsvequalizer.hcurve		= mods.hsvequalizer.hcurve;
	if (hsvequalizer.scurve)		toEdit.hsvequalizer.scurve		= mods.hsvequalizer.scurve;
	if (hsvequalizer.vcurve)		toEdit.hsvequalizer.vcurve		= mods.hsvequalizer.vcurve;

	if (filmSimulation.enabled)			toEdit.filmSimulation.enabled		= mods.filmSimulation.enabled;
	if (filmSimulation.clutFilename)	toEdit.filmSimulation.clutFilename	= mods.filmSimulation.clutFilename;
	if (filmSimulation.strength)		toEdit.filmSimulation.strength		= dontforceSet && options.baBehav[ADDSET_FILMSIMULATION_STRENGTH] ? toEdit.filmSimulation.strength + mods.filmSimulation.strength : mods.filmSimulation.strength;


	// Exif changes are added to the existing ones
	if (exif)
		for (procparams::ExifPairs::const_iterator i=mods.exif.begin(); i!=mods.exif.end(); i++) {
			toEdit.exif[i->first] = i->second;
		}

	// IPTC changes are added to the existing ones
	if (iptc)
		for (procparams::IPTCPairs::const_iterator i=mods.iptc.begin(); i!=mods.iptc.end(); i++) {
			toEdit.iptc[i->first] = i->second;
		}
}

bool RAWParamsEdited::BayerSensor::isUnchanged() const {
	return  method && ccSteps && dcbIterations && dcbEnhance && lmmseIterations/*&& allEnhance*/ &&  greenEq
			&& linenoise && exBlack0 && exBlack1 && exBlack2 && exBlack3 && exTwoGreen;
}

bool RAWParamsEdited::XTransSensor::isUnchanged() const {
	return method && ccSteps && exBlackRed && exBlackGreen && exBlackBlue;
}

bool RAWParamsEdited::isUnchanged() const {
	return  bayersensor.isUnchanged() && xtranssensor.isUnchanged() && caCorrection && caRed && caBlue && hotPixelFilter && deadPixelFilter && hotDeadPixelThresh && darkFrame
			&& dfAuto && ff_file && ff_AutoSelect && ff_BlurRadius && ff_BlurType && exPos && exPreser && ff_AutoClipControl && ff_clipControl;
}

bool LensProfParamsEdited::isUnchanged() const {
    return lcpFile;
}
