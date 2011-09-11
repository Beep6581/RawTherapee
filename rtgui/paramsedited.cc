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
#include <paramsedited.h>
#include <string.h>
#include <options.h>
#include <addsetids.h>

ParamsEdited::ParamsEdited () {

    set (false);
}

void ParamsEdited::set (bool v) {

	toneCurve.curve = v;
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
	labCurve.lcurve      = v;
	labCurve.acurve      = v;
	labCurve.bcurve      = v;
	labCurve.brightness = v;
	labCurve.contrast   = v;
	labCurve.saturation = v;
	labCurve.avoidclip     = v;
	labCurve.enable_saturationlimiter = v;
	labCurve.saturationlimit = v;	
	sharpening.enabled   = v;
	sharpening.radius    = v;
	sharpening.amount    = v;
	sharpening.threshold     = v;
	sharpening.edgesonly     = v;
	sharpening.edges_radius  = v;
	sharpening.edges_tolerance   = v;
	sharpening.halocontrol   = v;
	sharpening.halocontrol_amount= v;
	sharpening.method        = v;
	sharpening.deconvamount  = v;
	sharpening.deconvradius  = v;
	sharpening.deconviter    = v;
	sharpening.deconvdamping = v;
	sharpenEdge.enabled       = v;
	sharpenEdge.passes        = v;
	sharpenEdge.amount      = v;
	sharpenEdge.threechannels = v;
	sharpenMicro.enabled      = v;
	sharpenMicro.matrix       = v;
	sharpenMicro.amount     = v;
	sharpenMicro.uniformity   = v;
	colorBoost.amount         = v;
	colorBoost.avoidclip      = v;
	colorBoost.enable_saturationlimiter = v;
	colorBoost.saturationlimit = v;
	wb.method        = v;
	wb.green         = v;
	wb.temperature   = v;
	colorShift.a     = v;
	colorShift.b     = v;
	lumaDenoise.enabled       = v;
	lumaDenoise.radius        = v;
	lumaDenoise.edgetolerance = v;
	colorDenoise.enabled      = v;
	colorDenoise.amount       = v;
	defringe.enabled          = v;
	defringe.radius           = v;
	defringe.threshold        = v;
	impulseDenoise.enabled    = v;
	impulseDenoise.thresh     = v;
	dirpyrDenoise.enabled     = v;
	dirpyrDenoise.luma        = v;
	dirpyrDenoise.chroma = v;
	dirpyrDenoise.gamma = v;
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
	distortion.uselensfun = v;
	distortion.amount = v;
	perspective.horizontal = v;
	perspective.vertical = v;
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
	hlrecovery.enabled   = v;
	hlrecovery.method    = v;
	resize.scale     = v;
	resize.appliesTo = v;
	resize.method    = v;
	resize.dataspec  = v;
	resize.width     = v;
	resize.height    = v;
	resize.enabled   = v;
	icm.input        = v;
    icm.blendCMSMatrix = v;
	icm.working      = v;
	icm.output       = v;
	icm.gamma		= v;
	icm.freegamma		= v;
	icm.gampos		= v;
	icm.slpos		= v;	
	raw.ccSteps = v;
	raw.dmethod = v;
	raw.dcbIterations = v;
	raw.dcbEnhance = v;
	raw.allEnhance = v;
	raw.caCorrection = v;
	raw.caBlue  = v;
	raw.caRed   = v;
	raw.greenEq = v;
	raw.hotDeadPixel = v;
	raw.linenoise = v;
	raw.darkFrame = v;
    raw.dfAuto = v;
	raw.ff_file = v;
	raw.ff_AutoSelect = v;
	raw.ff_BlurRadius = v;
	raw.ff_BlurType = v;
	raw.exPos = v;
	raw.exPreser = v;
	raw.exBlackzero = v;
	raw.exBlackone = v;
	raw.exBlacktwo = v;
	raw.exBlackthree = v;
	raw.exTwoGreen=v;

	dirpyrequalizer.enabled = v;
	for(int i = 0; i < 5; i++) {
		dirpyrequalizer.mult[i] = v;
	}
	hsvequalizer.hcurve = v;
	hsvequalizer.scurve = v;
	hsvequalizer.vcurve = v;
	exif.clear ();
	iptc.clear ();
}

using namespace rtengine;
using namespace rtengine::procparams;

void ParamsEdited::initFrom (const std::vector<rtengine::procparams::ProcParams>& src) {

    set (true);
    if (src.size()==0)
        return;

    const ProcParams& p = src[0];
    for (int i=1; i<src.size(); i++) {
        const ProcParams& other = src[i];
        toneCurve.curve = toneCurve.curve && p.toneCurve.curve == other.toneCurve.curve;
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
        labCurve.lcurve = labCurve.lcurve && p.labCurve.lcurve == other.labCurve.lcurve;
        labCurve.acurve = labCurve.acurve && p.labCurve.acurve == other.labCurve.acurve;
        labCurve.bcurve = labCurve.bcurve && p.labCurve.bcurve == other.labCurve.bcurve;
        labCurve.brightness = labCurve.brightness && p.labCurve.brightness == other.labCurve.brightness;
        labCurve.contrast = labCurve.contrast && p.labCurve.contrast == other.labCurve.contrast;
        labCurve.saturation = labCurve.saturation && p.labCurve.saturation == other.labCurve.saturation;
		labCurve.avoidclip = labCurve.avoidclip && p.labCurve.avoidclip == other.labCurve.avoidclip;
        labCurve.enable_saturationlimiter = labCurve.enable_saturationlimiter && p.labCurve.enable_saturationlimiter == other.labCurve.enable_saturationlimiter;
        labCurve.saturationlimit = labCurve.saturationlimit && p.labCurve.saturationlimit == other.labCurve.saturationlimit;		
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
        colorBoost.amount = colorBoost.amount && p.colorBoost.amount == other.colorBoost.amount;
        colorBoost.avoidclip = colorBoost.avoidclip && p.colorBoost.avoidclip == other.colorBoost.avoidclip;
        colorBoost.enable_saturationlimiter = colorBoost.enable_saturationlimiter && p.colorBoost.enable_saturationlimiter == other.colorBoost.enable_saturationlimiter;
        colorBoost.saturationlimit = colorBoost.saturationlimit && p.colorBoost.saturationlimit == other.colorBoost.saturationlimit;
        wb.method = wb.method && p.wb.method == other.wb.method;
        wb.green = wb.green && p.wb.green == other.wb.green;
        wb.temperature = wb.temperature && p.wb.temperature == other.wb.temperature;
        colorShift.a = colorShift.a && p.colorShift.a == other.colorShift.a;
        colorShift.b = colorShift.b && p.colorShift.b == other.colorShift.b;
        lumaDenoise.enabled = lumaDenoise.enabled && p.lumaDenoise.enabled == other.lumaDenoise.enabled;
        lumaDenoise.radius = lumaDenoise.radius && p.lumaDenoise.radius == other.lumaDenoise.radius;
        lumaDenoise.edgetolerance = lumaDenoise.edgetolerance && p.lumaDenoise.edgetolerance == other.lumaDenoise.edgetolerance;
        colorDenoise.enabled = colorDenoise.enabled && p.colorDenoise.enabled == other.colorDenoise.enabled;
        colorDenoise.amount = colorDenoise.amount && p.colorDenoise.amount == other.colorDenoise.amount;
        defringe.enabled = defringe.enabled && p.defringe.enabled == other.defringe.enabled;
        defringe.radius = defringe.radius && p.defringe.radius == other.defringe.radius;
        defringe.threshold = defringe.threshold && p.defringe.threshold == other.defringe.threshold;
		
        impulseDenoise.enabled = impulseDenoise.enabled && p.impulseDenoise.enabled == other.impulseDenoise.enabled;
        impulseDenoise.thresh = impulseDenoise.thresh && p.impulseDenoise.thresh == other.impulseDenoise.thresh;

        dirpyrDenoise.enabled = dirpyrDenoise.enabled && p.dirpyrDenoise.enabled == other.dirpyrDenoise.enabled;
        dirpyrDenoise.luma = dirpyrDenoise.luma && p.dirpyrDenoise.luma == other.dirpyrDenoise.luma;
        dirpyrDenoise.chroma = dirpyrDenoise.chroma && p.dirpyrDenoise.chroma == other.dirpyrDenoise.chroma;
        dirpyrDenoise.gamma = dirpyrDenoise.gamma && p.dirpyrDenoise.gamma == other.dirpyrDenoise.gamma;

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
        distortion.uselensfun = distortion.uselensfun && p.distortion.uselensfun == other.distortion.uselensfun;
        distortion.amount = distortion.amount && p.distortion.amount == other.distortion.amount;
        perspective.horizontal = perspective.horizontal && p.perspective.horizontal == other.perspective.horizontal;
        perspective.vertical = perspective.vertical && p.perspective.vertical == other.perspective.vertical;
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
        hlrecovery.enabled = hlrecovery.enabled && p.hlrecovery.enabled == other.hlrecovery.enabled;
        hlrecovery.method = hlrecovery.method && p.hlrecovery.method == other.hlrecovery.method;
        resize.scale = resize.scale && p.resize.scale == other.resize.scale;
        resize.appliesTo = resize.appliesTo && p.resize.appliesTo == other.resize.appliesTo;
        resize.method = resize.method && p.resize.method == other.resize.method;
        resize.dataspec = resize.dataspec && p.resize.dataspec == other.resize.dataspec;
        resize.width = resize.width && p.resize.width == other.resize.width;
        resize.height = resize.height && p.resize.height == other.resize.height;
        resize.enabled = resize.enabled && p.resize.enabled == other.resize.enabled;
        icm.input = icm.input && p.icm.input == other.icm.input;
        icm.blendCMSMatrix = icm.blendCMSMatrix && p.icm.blendCMSMatrix == other.icm.blendCMSMatrix;
        icm.working = icm.working && p.icm.working == other.icm.working;
        icm.output = icm.output && p.icm.output == other.icm.output;
        icm.gamma = icm.gamma && p.icm.gamma == other.icm.gamma;
        icm.freegamma = icm.freegamma && p.icm.freegamma == other.icm.freegamma;		
        icm.gampos = icm.gampos && p.icm.gampos == other.icm.gampos;
        icm.slpos = icm.slpos && p.icm.slpos == other.icm.slpos;
        raw.ccSteps = raw.ccSteps && p.raw.ccSteps == other.raw.ccSteps;
        raw.dmethod = raw.dmethod && p.raw.dmethod == other.raw.dmethod;
        raw.dcbIterations = raw.dcbIterations && p.raw.dcb_iterations == other.raw.dcb_iterations;
        raw.dcbEnhance = raw.dcbEnhance && p.raw.dcb_enhance == other.raw.dcb_enhance;
        raw.allEnhance = raw.allEnhance && p.raw.all_enhance == other.raw.all_enhance;
        raw.caCorrection = raw.caCorrection && p.raw.ca_autocorrect == other.raw.ca_autocorrect;
		raw.caRed = raw.caRed && p.raw.cared == other.raw.cared;
        raw.caBlue = raw.caBlue && p.raw.cablue == other.raw.cablue;
        raw.greenEq = raw.greenEq && p.raw.greenthresh == other.raw.greenthresh;
        raw.hotDeadPixel = raw.hotDeadPixel && p.raw.hotdeadpix_filt == other.raw.hotdeadpix_filt;
        raw.linenoise = raw.linenoise && p.raw.linenoise == other.raw.linenoise;
        raw.darkFrame = raw.darkFrame && p.raw.dark_frame == other.raw.dark_frame;
        raw.dfAuto = raw.dfAuto && p.raw.df_autoselect == other.raw.df_autoselect;
        raw.ff_file = raw.ff_file && p.raw.ff_file == other.raw.ff_file;                        
        raw.ff_AutoSelect = raw.ff_AutoSelect && p.raw.ff_AutoSelect == other.raw.ff_AutoSelect;
        raw.ff_BlurRadius = raw.ff_BlurRadius && p.raw.ff_BlurRadius == other.raw.ff_BlurRadius;
        raw.ff_BlurType = raw.ff_BlurType && p.raw.ff_BlurType == other.raw.ff_BlurType;        
		raw.exPos = raw.exPos && p.raw.expos == other.raw.expos;
		raw.exPreser = raw.exPreser && p.raw.preser == other.raw.preser;
 		raw.exBlackzero = raw.exBlackzero && p.raw.blackzero == other.raw.blackzero;
		raw.exBlackone = raw.exBlackone && p.raw.blackone == other.raw.blackone;
		raw.exBlacktwo = p.raw.blacktwo == other.raw.blacktwo;
		raw.exBlackthree = p.raw.blackthree == other.raw.blackthree;
		raw.exTwoGreen = p.raw.twogreen == other.raw.twogreen;		

        dirpyrequalizer.enabled = dirpyrequalizer.enabled && p.dirpyrequalizer.enabled == other.dirpyrequalizer.enabled;
        for(int i = 0; i < 8; i++) {
            dirpyrequalizer.mult[i] = dirpyrequalizer.mult[i] && p.dirpyrequalizer.mult[i] == other.dirpyrequalizer.mult[i];
        }		
        hsvequalizer.hcurve = hsvequalizer.hcurve && p.hsvequalizer.hcurve == other.hsvequalizer.hcurve;
        hsvequalizer.scurve = hsvequalizer.scurve && p.hsvequalizer.scurve == other.hsvequalizer.scurve;
        hsvequalizer.vcurve = hsvequalizer.vcurve && p.hsvequalizer.vcurve == other.hsvequalizer.vcurve;
//        exif = exif && p.exif==other.exif
//        iptc = other.iptc;
    }
}

void ParamsEdited::combine (rtengine::procparams::ProcParams& toEdit, const rtengine::procparams::ProcParams& mods, bool forceSet) {

	bool dontforceSet = !forceSet;

	if (toneCurve.curve)		toEdit.toneCurve.curve 	    = mods.toneCurve.curve;
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
	if (labCurve.lcurve)		toEdit.labCurve.lcurve 	    = mods.labCurve.lcurve;
	if (labCurve.acurve)		toEdit.labCurve.acurve 	    = mods.labCurve.acurve;
	if (labCurve.bcurve)		toEdit.labCurve.bcurve 	    = mods.labCurve.bcurve;
	if (labCurve.brightness)	toEdit.labCurve.brightness = dontforceSet && options.baBehav[ADDSET_LC_BRIGHTNESS] ? toEdit.labCurve.brightness + mods.labCurve.brightness : mods.labCurve.brightness;
	if (labCurve.contrast)		toEdit.labCurve.contrast 	= dontforceSet && options.baBehav[ADDSET_LC_CONTRAST] ? toEdit.labCurve.contrast + mods.labCurve.contrast : mods.labCurve.contrast;
	if (labCurve.saturation)	toEdit.labCurve.saturation = dontforceSet && options.baBehav[ADDSET_LC_SATURATION] ? toEdit.labCurve.saturation + mods.labCurve.saturation : mods.labCurve.saturation;
	if (labCurve.avoidclip)					toEdit.labCurve.avoidclip 	= mods.labCurve.avoidclip;
	if (labCurve.enable_saturationlimiter)	toEdit.labCurve.enable_saturationlimiter 	= mods.labCurve.enable_saturationlimiter;
	if (labCurve.saturationlimit)			toEdit.labCurve.saturationlimit 	= mods.labCurve.saturationlimit;	
	if (sharpenEdge.enabled)				toEdit.sharpenEdge.enabled 	= mods.sharpenEdge.enabled;
	if (sharpenEdge.passes)					toEdit.sharpenEdge.passes	= dontforceSet && options.baBehav[ADDSET_SHARPENEDGE_PASS] ? toEdit.sharpenEdge.passes + mods.sharpenEdge.passes : mods.sharpenEdge.passes;
	if (sharpenEdge.amount)				toEdit.sharpenEdge.amount	= dontforceSet && options.baBehav[ADDSET_SHARPENEDGE_AMOUNT] ? toEdit.sharpenEdge.amount + mods.sharpenEdge.amount : mods.sharpenEdge.amount;
	if (sharpenEdge.threechannels)			toEdit.sharpenEdge.threechannels 	= mods.sharpenEdge.threechannels;
	if (sharpenMicro.enabled)				toEdit.sharpenMicro.enabled 	= mods.sharpenMicro.enabled;
	if (sharpenMicro.matrix)				toEdit.sharpenMicro.matrix	= mods.sharpenMicro.matrix;
	if (sharpenMicro.amount)				toEdit.sharpenMicro.amount	= dontforceSet && options.baBehav[ADDSET_SHARPENMICRO_AMOUNT] ? toEdit.sharpenMicro.amount + mods.sharpenMicro.amount : mods.sharpenMicro.amount;
	if (sharpenMicro.uniformity)			toEdit.sharpenMicro.uniformity	= dontforceSet && options.baBehav[ADDSET_SHARPENMICRO_UNIFORMITY] ? toEdit.sharpenMicro.uniformity + mods.sharpenMicro.uniformity : mods.sharpenMicro.uniformity;
	if (sharpening.enabled)					toEdit.sharpening.enabled 	= mods.sharpening.enabled;
	if (sharpening.radius)					toEdit.sharpening.radius 	= mods.sharpening.radius;
	if (sharpening.amount)					toEdit.sharpening.amount 	= dontforceSet && options.baBehav[ADDSET_SHARP_AMOUNT] ? toEdit.sharpening.amount + mods.sharpening.amount : mods.sharpening.amount;
	if (sharpening.threshold)				toEdit.sharpening.threshold 	= mods.sharpening.threshold;
	if (sharpening.edgesonly)				toEdit.sharpening.edgesonly 	= mods.sharpening.edgesonly;
	if (sharpening.edges_radius)			toEdit.sharpening.edges_radius 	= mods.sharpening.edges_radius;
	if (sharpening.edges_tolerance)			toEdit.sharpening.edges_tolerance 	= mods.sharpening.edges_tolerance;
	if (sharpening.halocontrol)				toEdit.sharpening.halocontrol 	= mods.sharpening.halocontrol;
	if (sharpening.halocontrol_amount)		toEdit.sharpening.halocontrol_amount 	= mods.sharpening.halocontrol_amount;
	if (sharpening.method)					toEdit.sharpening.method 	= mods.sharpening.method;
	if (sharpening.deconvamount)			toEdit.sharpening.deconvamount 	= dontforceSet && options.baBehav[ADDSET_SHARP_AMOUNT] ? toEdit.sharpening.deconvamount + mods.sharpening.deconvamount : mods.sharpening.deconvamount;
	if (sharpening.deconvradius)			toEdit.sharpening.deconvradius 	= mods.sharpening.deconvradius;
	if (sharpening.deconviter)				toEdit.sharpening.deconviter 	= mods.sharpening.deconviter;
	if (sharpening.deconvdamping)			toEdit.sharpening.deconvdamping 	= mods.sharpening.deconvdamping;
	if (colorBoost.amount)					toEdit.colorBoost.amount = dontforceSet && options.baBehav[ADDSET_CBOOST_AMOUNT] ? toEdit.colorBoost.amount + mods.colorBoost.amount : mods.colorBoost.amount;
	if (colorBoost.avoidclip)				toEdit.colorBoost.avoidclip 	= mods.colorBoost.avoidclip;
	if (colorBoost.enable_saturationlimiter)toEdit.colorBoost.enable_saturationlimiter 	= mods.colorBoost.enable_saturationlimiter;
	if (colorBoost.saturationlimit)			toEdit.colorBoost.saturationlimit 	= mods.colorBoost.saturationlimit;
	if (wb.method)							toEdit.wb.method 	= mods.wb.method;
	if (wb.green)							toEdit.wb.green 	= dontforceSet && options.baBehav[ADDSET_WB_GREEN] ? toEdit.wb.green + mods.wb.green : mods.wb.green;
	if (wb.temperature)						toEdit.wb.temperature 	= dontforceSet && options.baBehav[ADDSET_WB_TEMPERATURE] ? toEdit.wb.temperature + mods.wb.temperature : mods.wb.temperature;
	if (colorShift.a)						toEdit.colorShift.a 	= dontforceSet && options.baBehav[ADDSET_CS_BLUEYELLOW] ? toEdit.colorShift.a + mods.colorShift.a : mods.colorShift.a;
	if (colorShift.b)						toEdit.colorShift.b 	= dontforceSet && options.baBehav[ADDSET_CS_GREENMAGENTA] ? toEdit.colorShift.b + mods.colorShift.b : mods.colorShift.b;
	if (lumaDenoise.enabled)				toEdit.lumaDenoise.enabled 	= mods.lumaDenoise.enabled;
	if (lumaDenoise.radius)					toEdit.lumaDenoise.radius 	= mods.lumaDenoise.radius;
	if (lumaDenoise.edgetolerance)			toEdit.lumaDenoise.edgetolerance 	= dontforceSet && options.baBehav[ADDSET_LD_EDGETOLERANCE] ? toEdit.lumaDenoise.edgetolerance + mods.lumaDenoise.edgetolerance : mods.lumaDenoise.edgetolerance;
	if (colorDenoise.enabled)				toEdit.colorDenoise.enabled 	= mods.colorDenoise.enabled;
	if (colorDenoise.amount)				toEdit.colorDenoise.amount 	= mods.colorDenoise.amount;
	
	if (defringe.enabled)					toEdit.defringe.enabled 	= mods.defringe.enabled;
	if (defringe.radius)					toEdit.defringe.radius 	= mods.defringe.radius;	
	if (defringe.threshold)					toEdit.defringe.threshold 	= mods.defringe.threshold;
	
	if (impulseDenoise.enabled)				toEdit.impulseDenoise.enabled 	= mods.impulseDenoise.enabled;
	if (impulseDenoise.thresh)				toEdit.impulseDenoise.thresh 	= mods.impulseDenoise.thresh;

	if (dirpyrDenoise.enabled)				toEdit.dirpyrDenoise.enabled 	= mods.dirpyrDenoise.enabled;
	if (dirpyrDenoise.luma)					toEdit.dirpyrDenoise.luma		= dontforceSet && options.baBehav[ADDSET_DIRPYRDN_CHLUM] ? toEdit.dirpyrDenoise.luma + mods.dirpyrDenoise.luma : mods.dirpyrDenoise.luma;
	if (dirpyrDenoise.chroma)				toEdit.dirpyrDenoise.chroma		= dontforceSet && options.baBehav[ADDSET_DIRPYRDN_CHLUM] ? toEdit.dirpyrDenoise.chroma + mods.dirpyrDenoise.chroma : mods.dirpyrDenoise.chroma;
	if (dirpyrDenoise.gamma)				toEdit.dirpyrDenoise.gamma		= dontforceSet && options.baBehav[ADDSET_DIRPYRDN_GAMMA] ? toEdit.dirpyrDenoise.gamma + mods.dirpyrDenoise.gamma : mods.dirpyrDenoise.gamma;

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
	if (coarse.rotate)						toEdit.coarse.rotate 	= (toEdit.coarse.rotate + mods.coarse.rotate) % 360;
	if (coarse.hflip)						toEdit.coarse.hflip 	= mods.coarse.hflip ? !toEdit.coarse.hflip : toEdit.coarse.hflip;
	if (coarse.vflip)						toEdit.coarse.vflip 	= mods.coarse.vflip ? !toEdit.coarse.vflip : toEdit.coarse.vflip;
	if (commonTrans.autofill)				toEdit.commonTrans.autofill		= mods.commonTrans.autofill;
	if (rotate.degree)						toEdit.rotate.degree 			= dontforceSet && options.baBehav[17] ? toEdit.rotate.degree + mods.rotate.degree : mods.rotate.degree;
	if (distortion.uselensfun)				toEdit.distortion.uselensfun	= mods.distortion.uselensfun;
	if (distortion.amount)					toEdit.distortion.amount 		= dontforceSet && options.baBehav[ADDSET_DIST_AMOUNT] ? toEdit.distortion.amount + mods.distortion.amount : mods.distortion.amount;
	if (perspective.horizontal)				toEdit.perspective.horizontal 	= dontforceSet && options.baBehav[ADDSET_PERSPECTIVE] ? toEdit.perspective.horizontal + mods.perspective.horizontal : mods.perspective.horizontal;
	if (perspective.vertical)				toEdit.perspective.vertical 	= dontforceSet && options.baBehav[ADDSET_PERSPECTIVE] ? toEdit.perspective.vertical + mods.perspective.vertical : mods.perspective.vertical;
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
	if (hlrecovery.enabled)	toEdit.hlrecovery.enabled 	= mods.hlrecovery.enabled;
	if (hlrecovery.method)	toEdit.hlrecovery.method 	= mods.hlrecovery.method;
	if (resize.scale)		toEdit.resize.scale 	= mods.resize.scale;
	if (resize.appliesTo)	toEdit.resize.appliesTo = mods.resize.appliesTo;
	if (resize.method)		toEdit.resize.method 	= mods.resize.method;
	if (resize.dataspec)	toEdit.resize.dataspec 	= mods.resize.dataspec;
	if (resize.width)	    toEdit.resize.width 	= mods.resize.width;
	if (resize.height)	    toEdit.resize.height 	= mods.resize.height;
	if (resize.enabled)	    toEdit.resize.enabled 	= mods.resize.enabled;
	if (icm.input)		    toEdit.icm.input 	    = mods.icm.input;
    if (icm.blendCMSMatrix)	toEdit.icm.blendCMSMatrix = mods.icm.blendCMSMatrix;
	if (icm.working)		toEdit.icm.working 	    = mods.icm.working;
	if (icm.output)		    toEdit.icm.output       = mods.icm.output;
	if (icm.gampos)		    toEdit.icm.gampos       = mods.icm.gampos;
	if (icm.slpos)		    toEdit.icm.slpos        = mods.icm.slpos;
	if (icm.gamma)		    toEdit.icm.gamma        = mods.icm.gamma;
	if (icm.freegamma)		toEdit.icm.freegamma    = mods.icm.freegamma;	
    if (raw.ccSteps)        toEdit.raw.ccSteps      = mods.raw.ccSteps;
    if (raw.dmethod)        toEdit.raw.dmethod      = mods.raw.dmethod;
    if (raw.dcbIterations)  toEdit.raw.dcb_iterations = mods.raw.dcb_iterations;
    if (raw.dcbEnhance)     toEdit.raw.dcb_enhance     = mods.raw.dcb_enhance;
    if (raw.allEnhance)     toEdit.raw.all_enhance     = mods.raw.all_enhance;
	
    if (raw.caCorrection)   toEdit.raw.ca_autocorrect  = mods.raw.ca_autocorrect;
	if (raw.caRed)          toEdit.raw.cared           = dontforceSet && options.baBehav[ADDSET_RAWCACORR] ? toEdit.raw.cared + mods.raw.cared : mods.raw.cared;
    if (raw.caBlue)         toEdit.raw.cablue          = dontforceSet && options.baBehav[ADDSET_RAWCACORR] ? toEdit.raw.cablue + mods.raw.cablue : mods.raw.cablue;
	if (raw.exPos)          toEdit.raw.expos           = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_LINEAR] ? toEdit.raw.expos + mods.raw.expos : mods.raw.expos;
	if (raw.exPreser)       toEdit.raw.preser          = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_PRESER] ? toEdit.raw.preser + mods.raw.preser : mods.raw.preser;
	if (raw.exBlackzero)    toEdit.raw.blackzero       = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit.raw.blackzero + mods.raw.blackzero : mods.raw.blackzero;
	if (raw.exBlackone)     toEdit.raw.blackone        = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit.raw.blackone + mods.raw.blackone : mods.raw.blackone;
	if (raw.exBlacktwo)     toEdit.raw.blacktwo        = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit.raw.blacktwo + mods.raw.blacktwo : mods.raw.blacktwo;
	if (raw.exBlackthree)   toEdit.raw.blackthree      = dontforceSet && options.baBehav[ADDSET_RAWEXPOS_BLACKS] ? toEdit.raw.blackthree + mods.raw.blackthree : mods.raw.blackthree;
	if (raw.exTwoGreen)		toEdit.raw.twogreen        = mods.raw.twogreen;	

    if (raw.greenEq)        toEdit.raw.greenthresh     = dontforceSet && options.baBehav[ADDSET_PREPROCESS_GREENEQUIL] ? toEdit.raw.greenthresh + mods.raw.greenthresh : mods.raw.greenthresh;
    if (raw.hotDeadPixel)   toEdit.raw.hotdeadpix_filt = mods.raw.hotdeadpix_filt;
    if (raw.linenoise)      toEdit.raw.linenoise       = dontforceSet && options.baBehav[ADDSET_PREPROCESS_LINEDENOISE] ? toEdit.raw.linenoise + mods.raw.linenoise : mods.raw.linenoise;
    if (raw.darkFrame)      toEdit.raw.dark_frame      = mods.raw.dark_frame;
    if (raw.dfAuto)         toEdit.raw.df_autoselect   = mods.raw.df_autoselect;

    if (raw.ff_file)        toEdit.raw.ff_file         = mods.raw.ff_file;         
    if (raw.ff_AutoSelect)  toEdit.raw.ff_AutoSelect   = mods.raw.ff_AutoSelect;
    if (raw.ff_BlurRadius)  toEdit.raw.ff_BlurRadius   = mods.raw.ff_BlurRadius;
    if (raw.ff_BlurType)    toEdit.raw.ff_BlurType     = mods.raw.ff_BlurType;      

	if (dirpyrequalizer.enabled)	toEdit.dirpyrequalizer.enabled	= mods.dirpyrequalizer.enabled;
	for(int i = 0; i < 5; i++) {
		if(dirpyrequalizer.mult[i])	toEdit.dirpyrequalizer.mult[i]	= dontforceSet && options.baBehav[ADDSET_DIRPYREQ] ? toEdit.dirpyrequalizer.mult[i] + mods.dirpyrequalizer.mult[i] : mods.dirpyrequalizer.mult[i];
	}
	if (hsvequalizer.hcurve)		toEdit.hsvequalizer.hcurve		= mods.hsvequalizer.hcurve;
	if (hsvequalizer.scurve)		toEdit.hsvequalizer.scurve		= mods.hsvequalizer.scurve;
	if (hsvequalizer.vcurve)		toEdit.hsvequalizer.vcurve		= mods.hsvequalizer.vcurve;

//	if (exif)		toEdit.exif==mo.exif 	= mods.exif==other.exif;
//	if (iptc;)		toEdit.iptc==other.iptc; 	= mods.iptc==other.iptc;;
}

bool RAWParamsEdited::isUnchanged() const {
    return ccSteps && dmethod && dcbIterations && dcbEnhance && allEnhance && caCorrection && caRed && caBlue && greenEq
        && hotDeadPixel && linenoise && darkFrame && dfAuto && ff_file && ff_AutoSelect && ff_BlurRadius && ff_BlurType
	    && exPos && exPreser && exBlackzero && exBlackone && exBlacktwo && exBlackthree && exTwoGreen;
}
