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

ParamsEdited::ParamsEdited () {

    set (false);
}

void ParamsEdited::set (bool v) {

        toneCurve.curve = v;
        toneCurve.brightness = v;
        toneCurve.black      = v;
        toneCurve.contrast   = v;
        toneCurve.shcompr    = v;
        toneCurve.hlcompr    = v;
        toneCurve.autoexp    = v;
        toneCurve.clip       = v;
        toneCurve.expcomp    = v;
        lumaCurve.curve      = v;
        lumaCurve.brightness = v;
        lumaCurve.black      = v;
        lumaCurve.contrast   = v;
        lumaCurve.shcompr    = v;
        lumaCurve.hlcompr    = v;
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
        colorBoost.amount        = v;
        colorBoost.avoidclip     = v;
        colorBoost.enable_saturationlimiter = v;
        colorBoost.saturationlimit = v;
        wb.method        = v;
        wb.green         = v;
        wb.temperature   = v;
        colorShift.a     = v;
        colorShift.b     = v;
        lumaDenoise.enabled      = v;
        lumaDenoise.radius       = v;
        lumaDenoise.edgetolerance = v;
        colorDenoise.enabled      = v;
        colorDenoise.amount       = v;
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
        rotate.degree = v;
        rotate.fill = v;
        distortion.amount = v;
        cacorrection.red = v;
        cacorrection.blue = v;
        vignetting.amount = v;
        vignetting.radius = v;
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
        resize.method    = v;
        resize.dataspec  = v;
        resize.width     = v;
        resize.height    = v;
        resize.enabled   = v;
        icm.input        = v;
        icm.gammaOnInput = v;
        icm.working      = v;
        icm.output       = v;
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
        toneCurve.shcompr = toneCurve.shcompr && p.toneCurve.shcompr == other.toneCurve.shcompr;
        toneCurve.hlcompr = toneCurve.hlcompr && p.toneCurve.hlcompr == other.toneCurve.hlcompr;
        toneCurve.autoexp = toneCurve.autoexp && p.toneCurve.autoexp == other.toneCurve.autoexp;
        toneCurve.clip = toneCurve.clip && p.toneCurve.clip == other.toneCurve.clip;
        toneCurve.expcomp = toneCurve.expcomp && p.toneCurve.expcomp == other.toneCurve.expcomp;
        lumaCurve.curve = lumaCurve.curve && p.lumaCurve.curve == other.lumaCurve.curve;
        lumaCurve.brightness = lumaCurve.brightness && p.lumaCurve.brightness == other.lumaCurve.brightness;
        lumaCurve.black = lumaCurve.black && p.lumaCurve.black == other.lumaCurve.black;
        lumaCurve.contrast = lumaCurve.contrast && p.lumaCurve.contrast == other.lumaCurve.contrast;
        lumaCurve.shcompr = lumaCurve.shcompr && p.lumaCurve.shcompr == other.lumaCurve.shcompr;
        lumaCurve.hlcompr = lumaCurve.hlcompr && p.lumaCurve.hlcompr == other.lumaCurve.hlcompr;
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
        rotate.degree = rotate.degree && p.rotate.degree == other.rotate.degree;
        rotate.fill = rotate.fill && p.rotate.fill == other.rotate.fill;
        distortion.amount = distortion.amount && p.distortion.amount == other.distortion.amount;
        cacorrection.red = cacorrection.red && p.cacorrection.red == other.cacorrection.red;
        cacorrection.blue = cacorrection.blue && p.cacorrection.blue == other.cacorrection.blue;
        vignetting.amount = vignetting.amount && p.vignetting.amount == other.vignetting.amount;
        vignetting.radius = vignetting.radius && p.vignetting.radius == other.vignetting.radius;
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
        resize.method = resize.method && p.resize.method == other.resize.method;
        resize.dataspec = resize.dataspec && p.resize.dataspec == other.resize.dataspec;
        resize.width = resize.width && p.resize.width == other.resize.width;
        resize.height = resize.height && p.resize.height == other.resize.height;
        resize.enabled = resize.enabled && p.resize.enabled == other.resize.enabled;
        icm.input = icm.input && p.icm.input == other.icm.input;
        icm.gammaOnInput = icm.gammaOnInput && p.icm.gammaOnInput == other.icm.gammaOnInput;
        icm.working = icm.working && p.icm.working == other.icm.working;
        icm.output = icm.output && p.icm.output == other.icm.output;
//        exif = exif && p.exif==other.exif
//        iptc = other.iptc;
    }
}

void ParamsEdited::combine (rtengine::procparams::ProcParams& toEdit, const rtengine::procparams::ProcParams& mods) {

	if (toneCurve.curve)	    toEdit.toneCurve.curve 	    = mods.toneCurve.curve;
	if (toneCurve.brightness)	toEdit.toneCurve.brightness = options.baBehav[1] ? toEdit.toneCurve.brightness + mods.toneCurve.brightness : mods.toneCurve.brightness;
	if (toneCurve.black)		toEdit.toneCurve.black 	    = options.baBehav[2] ? toEdit.toneCurve.black + mods.toneCurve.black : mods.toneCurve.black;
	if (toneCurve.contrast)		toEdit.toneCurve.contrast 	= options.baBehav[3] ? toEdit.toneCurve.contrast + mods.toneCurve.contrast : mods.toneCurve.contrast;
	if (toneCurve.shcompr)		toEdit.toneCurve.shcompr 	= mods.toneCurve.shcompr;
	if (toneCurve.hlcompr)		toEdit.toneCurve.hlcompr 	= mods.toneCurve.hlcompr;
	if (toneCurve.autoexp)		toEdit.toneCurve.autoexp 	= mods.toneCurve.autoexp;
	if (toneCurve.clip)		    toEdit.toneCurve.clip 	    = mods.toneCurve.clip;
	if (toneCurve.expcomp)		toEdit.toneCurve.expcomp 	= options.baBehav[0] ? toEdit.toneCurve.expcomp + mods.toneCurve.expcomp : mods.toneCurve.expcomp;
	if (lumaCurve.curve)		toEdit.lumaCurve.curve 	    = mods.lumaCurve.curve;
	if (lumaCurve.brightness)	toEdit.lumaCurve.brightness = mods.lumaCurve.brightness;
	if (lumaCurve.black)		toEdit.lumaCurve.black 	    = mods.lumaCurve.black;
	if (lumaCurve.contrast)		toEdit.lumaCurve.contrast 	= mods.lumaCurve.contrast;
	if (lumaCurve.shcompr)		toEdit.lumaCurve.shcompr 	= mods.lumaCurve.shcompr;
	if (lumaCurve.hlcompr)		toEdit.lumaCurve.hlcompr 	= mods.lumaCurve.hlcompr;
	if (sharpening.enabled)		toEdit.sharpening.enabled 	= mods.sharpening.enabled;
	if (sharpening.radius)		toEdit.sharpening.radius 	= mods.sharpening.radius;
	if (sharpening.amount)		toEdit.sharpening.amount 	= options.baBehav[10] ? toEdit.sharpening.amount + mods.sharpening.amount : mods.sharpening.amount;
	if (sharpening.threshold)	toEdit.sharpening.threshold 	= mods.sharpening.threshold;
	if (sharpening.edgesonly)	toEdit.sharpening.edgesonly 	= mods.sharpening.edgesonly;
	if (sharpening.edges_radius)	toEdit.sharpening.edges_radius 	= mods.sharpening.edges_radius;
	if (sharpening.edges_tolerance)	toEdit.sharpening.edges_tolerance 	= mods.sharpening.edges_tolerance;
	if (sharpening.halocontrol)		toEdit.sharpening.halocontrol 	= mods.sharpening.halocontrol;
	if (sharpening.halocontrol_amount)	toEdit.sharpening.halocontrol_amount 	= mods.sharpening.halocontrol_amount;
	if (sharpening.method)		    toEdit.sharpening.method 	= mods.sharpening.method;
	if (sharpening.deconvamount)	toEdit.sharpening.deconvamount 	= options.baBehav[10] ? toEdit.sharpening.deconvamount + mods.sharpening.deconvamount : mods.sharpening.deconvamount;
	if (sharpening.deconvradius)	toEdit.sharpening.deconvradius 	= mods.sharpening.deconvradius;
	if (sharpening.deconviter)		toEdit.sharpening.deconviter 	= mods.sharpening.deconviter;
	if (sharpening.deconvdamping)	toEdit.sharpening.deconvdamping 	= mods.sharpening.deconvdamping;
	if (colorBoost.amount)		    toEdit.colorBoost.amount = options.baBehav[14] ? toEdit.colorBoost.amount + mods.colorBoost.amount : mods.colorBoost.amount;
	if (colorBoost.avoidclip)	toEdit.colorBoost.avoidclip 	= mods.colorBoost.avoidclip;
	if (colorBoost.enable_saturationlimiter)		toEdit.colorBoost.enable_saturationlimiter 	= mods.colorBoost.enable_saturationlimiter;
	if (colorBoost.saturationlimit)		toEdit.colorBoost.saturationlimit 	= mods.colorBoost.saturationlimit;
	if (wb.method)		toEdit.wb.method 	= mods.wb.method;
	if (wb.green)		toEdit.wb.green 	= options.baBehav[13] ? toEdit.wb.green + mods.wb.green : mods.wb.green;
	if (wb.temperature)		toEdit.wb.temperature 	= options.baBehav[12] ? toEdit.wb.temperature + mods.wb.temperature : mods.wb.temperature;
	if (colorShift.a)		toEdit.colorShift.a 	= options.baBehav[15] ? toEdit.colorShift.a + mods.colorShift.a : mods.colorShift.a;
	if (colorShift.b)		toEdit.colorShift.b 	= options.baBehav[16] ? toEdit.colorShift.b + mods.colorShift.b : mods.colorShift.b;
	if (lumaDenoise.enabled)		toEdit.lumaDenoise.enabled 	= mods.lumaDenoise.enabled;
	if (lumaDenoise.radius)		toEdit.lumaDenoise.radius 	= mods.lumaDenoise.radius;
	if (lumaDenoise.edgetolerance)		toEdit.lumaDenoise.edgetolerance 	= options.baBehav[11] ? toEdit.lumaDenoise.edgetolerance + mods.lumaDenoise.edgetolerance : mods.lumaDenoise.edgetolerance;
	if (colorDenoise.enabled)		toEdit.colorDenoise.enabled 	= mods.colorDenoise.enabled;
	if (colorDenoise.amount)		toEdit.colorDenoise.amount 	= mods.colorDenoise.amount;
	if (sh.enabled)		    toEdit.sh.enabled 	    = mods.sh.enabled;
	if (sh.hq)		        toEdit.sh.hq     	    = mods.sh.hq;
	if (sh.highlights)		toEdit.sh.highlights 	= options.baBehav[4] ? toEdit.sh.highlights + mods.sh.highlights : mods.sh.highlights;
	if (sh.htonalwidth)		toEdit.sh.htonalwidth 	= mods.sh.htonalwidth;
	if (sh.shadows)		    toEdit.sh.shadows 	    = options.baBehav[5] ? toEdit.sh.shadows + mods.sh.shadows : mods.sh.shadows;
	if (sh.stonalwidth)		toEdit.sh.stonalwidth 	= mods.sh.stonalwidth;
	if (sh.localcontrast)	toEdit.sh.localcontrast = options.baBehav[6] ? toEdit.sh.localcontrast + mods.sh.localcontrast : mods.sh.localcontrast;
	if (sh.radius)		    toEdit.sh.radius 	    = mods.sh.radius;
	if (crop.enabled)		toEdit.crop.enabled = mods.crop.enabled;
	if (crop.x)		        toEdit.crop.x 	    = mods.crop.x;
	if (crop.y)		        toEdit.crop.y 	    = mods.crop.y;
	if (crop.w)		        toEdit.crop.w 	    = mods.crop.w;
	if (crop.h)		        toEdit.crop.h 	    = mods.crop.h;
	if (crop.fixratio)		toEdit.crop.fixratio 	= mods.crop.fixratio;
	if (crop.ratio)		    toEdit.crop.ratio 	    = mods.crop.ratio;
	if (crop.orientation)	toEdit.crop.orientation = mods.crop.orientation;
	if (crop.guide)		    toEdit.crop.guide 	    = mods.crop.guide;
	if (coarse.rotate)		toEdit.coarse.rotate 	= (toEdit.coarse.rotate + mods.coarse.rotate) % 360;
	if (coarse.hflip)		toEdit.coarse.hflip 	= mods.coarse.hflip ? !toEdit.coarse.hflip : toEdit.coarse.hflip;
	if (coarse.vflip)		toEdit.coarse.vflip 	= mods.coarse.vflip ? !toEdit.coarse.vflip : toEdit.coarse.vflip;
	if (rotate.degree)		toEdit.rotate.degree 	= options.baBehav[17] ? toEdit.rotate.degree + mods.rotate.degree : mods.rotate.degree;
	if (rotate.fill)		toEdit.rotate.fill 	    = mods.rotate.fill;
	if (distortion.amount)	toEdit.distortion.amount 	= options.baBehav[18] ? toEdit.distortion.amount + mods.distortion.amount : mods.distortion.amount;
	if (cacorrection.red)	toEdit.cacorrection.red 	= options.baBehav[19] ? toEdit.cacorrection.red + mods.cacorrection.red : mods.cacorrection.red;
	if (cacorrection.blue)	toEdit.cacorrection.blue 	= options.baBehav[20] ? toEdit.cacorrection.blue + mods.cacorrection.blue : mods.cacorrection.blue;
	if (vignetting.amount)	toEdit.vignetting.amount 	= options.baBehav[21] ? toEdit.vignetting.amount + mods.vignetting.amount : mods.vignetting.amount;
	if (vignetting.radius)	toEdit.vignetting.radius 	= mods.vignetting.radius;
	if (chmixer.red[0])		toEdit.chmixer.red[0] 	= mods.chmixer.red[0];
	if (chmixer.red[1])		toEdit.chmixer.red[1] 	= mods.chmixer.red[1];
	if (chmixer.red[2])		toEdit.chmixer.red[2] 	= mods.chmixer.red[2];
	if (chmixer.green[0])	toEdit.chmixer.green[0]	= mods.chmixer.green[0];
	if (chmixer.green[1])	toEdit.chmixer.green[1]	= mods.chmixer.green[1];
	if (chmixer.green[2])	toEdit.chmixer.green[2]	= mods.chmixer.green[2];
	if (chmixer.blue[0])	toEdit.chmixer.blue[0] 	= mods.chmixer.blue[0];
	if (chmixer.blue[1])	toEdit.chmixer.blue[1] 	= mods.chmixer.blue[1];
	if (chmixer.blue[2])	toEdit.chmixer.blue[2] 	= mods.chmixer.blue[2];
	if (hlrecovery.enabled)	toEdit.hlrecovery.enabled 	= mods.hlrecovery.enabled;
	if (hlrecovery.method)	toEdit.hlrecovery.method 	= mods.hlrecovery.method;
	if (resize.scale)		toEdit.resize.scale 	= mods.resize.scale;
	if (resize.method)		toEdit.resize.method 	= mods.resize.method;
	if (resize.dataspec)	toEdit.resize.dataspec 	= mods.resize.dataspec;
	if (resize.width)	    toEdit.resize.width 	= mods.resize.width;
	if (resize.height)	    toEdit.resize.height 	= mods.resize.height;
	if (resize.enabled)	    toEdit.resize.enabled 	= mods.resize.enabled;
	if (icm.input)		    toEdit.icm.input 	    = mods.icm.input;
	if (icm.gammaOnInput)	toEdit.icm.gammaOnInput = mods.icm.gammaOnInput;
	if (icm.working)		toEdit.icm.working 	    = mods.icm.working;
	if (icm.output)		    toEdit.icm.output 	    = mods.icm.output;
//	if (exif)		toEdit.exif==mo.exif 	= mods.exif==other.exif;
//	if (iptc;)		toEdit.iptc==other.iptc; 	= mods.iptc==other.iptc;;
}
