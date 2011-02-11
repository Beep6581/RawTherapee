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
#include <glib/gstdio.h>
#include <safegtk.h>
#include <procparams.h>
#include <glibmm.h>
#include <sstream>
#include <string.h>

#include <safekeyfile.h>

namespace rtengine {
namespace procparams {

const char *RAWParams::methodstring[RAWParams::numMethods]={"eahd", "hphd", "vng4", "dcb", "amaze", "ahd", "fast" };
const char *RAWParams::ff_BlurTypestring[RAWParams::numFlatFileBlurTypes]={/*"Parametric",*/ "Area Flatfield", "Vertical Flatfield", "Horizontal Flatfield", "V+H Flatfield"};

ProcParams::ProcParams () { 

    setDefaults (); 
}       

ProcParams* ProcParams::create () {

    return new ProcParams();
}

void ProcParams::destroy (ProcParams* pp) {

    delete pp;
}

void ProcParams::setDefaults () {

    toneCurve.autoexp       = false;
    toneCurve.clip          = 0.001;
    toneCurve.expcomp       = 0;
    toneCurve.brightness    = 0;
    toneCurve.contrast      = 0;
	toneCurve.saturation    = 0;
    toneCurve.black         = 0;
    toneCurve.hlcompr       = 70;
    toneCurve.hlcomprthresh = 0;
    toneCurve.shcompr       = 25;
    toneCurve.curve.clear ();
    
    labCurve.brightness    = 0;
    labCurve.contrast      = 0;
	labCurve.saturation      = 0;
    labCurve.lcurve.clear ();
	labCurve.acurve.clear ();
    labCurve.bcurve.clear ();
    
    sharpening.enabled          = true;
    sharpening.radius           = 1.0;
    sharpening.amount           = 90;
    sharpening.threshold        = 768;
    sharpening.edgesonly        = false;
    sharpening.edges_radius     = 3;
    sharpening.edges_tolerance  = 1000;
    sharpening.halocontrol      = false;
    sharpening.halocontrol_amount = 85;
    sharpening.method           = "usm";
    sharpening.deconvradius     = 0.75;
    sharpening.deconviter       = 30;
    sharpening.deconvdamping    = 20;
    sharpening.deconvamount     = 75;
    
    colorBoost.amount                   = 0;
    colorBoost.avoidclip                = false;
    colorBoost.enable_saturationlimiter = false;
    colorBoost.saturationlimit          = 50;
    
    wb.method       = "Camera";
    wb.temperature  = 6504;
    wb.green        = 1.00102;
    
    colorShift.a    = 0;
    colorShift.b    = 0;
	    
    lumaDenoise.enabled         = false;
    lumaDenoise.radius          = 1.9;
    lumaDenoise.edgetolerance   = 2000;
    
    colorDenoise.enabled        = false;
    colorDenoise.edgesensitive  = false;
    colorDenoise.radius         = 1.9;
    colorDenoise.edgetolerance  = 2000;
	
	impulseDenoise.enabled      = false;
	impulseDenoise.thresh		= 50;
	
	defringe.enabled			= false;
    defringe.radius				= 2.0;
    defringe.threshold			= 25;

	dirpyrDenoise.enabled       = false;
    dirpyrDenoise.luma          = 10;
    dirpyrDenoise.chroma		= 10;
	dirpyrDenoise.gamma			= 2.0;
    
    sh.enabled       = false;
    sh.hq            = false;
    sh.highlights    = 0;
    sh.htonalwidth   = 80;
    sh.shadows       = 0;
    sh.stonalwidth   = 80;
    sh.localcontrast = 0;
    sh.radius        = 40;
    
    crop.enabled    = false;
    crop.x          = -1;
    crop.y          = -1;
    crop.w          = 15000;
    crop.h          = 15000;
    crop.fixratio   = false;
    crop.ratio      = "3:2";
    crop.orientation= "Landscape";
    crop.guide      = "None";
    
    coarse.rotate   = 0;
    coarse.hflip    = false;
    coarse.vflip    = false;
    
    commonTrans.autofill = true;

    rotate.degree       = 0;

    distortion.amount     = 0;
    distortion.uselensfun = false;
    
    perspective.horizontal = 0;
    perspective.vertical   = 0;

    cacorrection.red  = 0;
    cacorrection.blue = 0;
    
    hlrecovery.enabled = false;
    hlrecovery.method  = "Luminance";

    vignetting.amount = 0;
    vignetting.radius = 50;
    vignetting.strength = 1;
    vignetting.centerX = 0;
    vignetting.centerY = 0;
    
    chmixer.red[0] = 100;
    chmixer.red[1] = 0;
    chmixer.red[2] = 0;
    chmixer.green[0] = 0;
    chmixer.green[1] = 100;
    chmixer.green[2] = 0;
    chmixer.blue[0] = 0;
    chmixer.blue[1] = 0;
    chmixer.blue[2] = 100;
    
    resize.enabled = false;
    resize.scale  = 1.0;
    resize.appliesTo = "Cropped area";
    resize.method = "Bicubic";
    resize.dataspec = 0;
    resize.width = 800;
    resize.height = 600;
    
    icm.input   = "";
    icm.gammaOnInput = false;
    icm.working = "sRGB";
    icm.output  = "sRGB";
    
    equalizer.enabled = false;    
    for(int i = 0; i < 8; i ++)
    {
        equalizer.c[i] = 0;
    }
	dirpyrequalizer.enabled = false;    
    for(int i = 0; i < 4; i ++)
    {
        dirpyrequalizer.mult[i] = 1.0;
    }
	dirpyrequalizer.mult[4] = 0.0;
	hsvequalizer.enabled = false;    
    for(int i = 0; i < 8; i ++)
    {
        hsvequalizer.sat[i] = 0;
		hsvequalizer.val[i] = 0;
        hsvequalizer.hue[i] = 0;
    }
    raw.df_autoselect = false;
    raw.ff_AutoSelect = false;                                      
    raw.ff_BlurRadius = 32;                                         
    raw.ff_BlurType = RAWParams::ff_BlurTypestring[RAWParams::area_ff];
    raw.cared = 0;
	raw.cablue = 0;
    raw.ca_autocorrect = false;
    raw.hotdeadpix_filt = false;
	raw.hotdeadpix_thresh = 40;
    raw.linenoise = 0;
    raw.greenthresh = 0;
    raw.ccSteps = 1;
    raw.dmethod = RAWParams::methodstring[RAWParams::hphd];;
    raw.dcb_iterations=2;
    raw.dcb_enhance=false;

    exif.clear ();
    iptc.clear ();
    
    version = 249;
}

int ProcParams::save (Glib::ustring fname) const {

    SafeKeyFile keyFile;

    keyFile.set_integer ("Version", "Version", 20101019);

    // save tonecurve:
    keyFile.set_boolean ("Exposure", "Auto",            toneCurve.autoexp);
    keyFile.set_double  ("Exposure", "Clip",            toneCurve.clip);
    keyFile.set_double  ("Exposure", "Compensation",    toneCurve.expcomp);
    keyFile.set_integer ("Exposure", "Brightness",      toneCurve.brightness);
    keyFile.set_integer ("Exposure", "Contrast",        toneCurve.contrast);
	keyFile.set_integer ("Exposure", "Saturation",      toneCurve.saturation);
    keyFile.set_integer ("Exposure", "Black",           toneCurve.black);
    keyFile.set_integer ("Exposure", "HighlightCompr",  toneCurve.hlcompr);
    keyFile.set_integer ("Exposure", "HighlightComprThreshold",  toneCurve.hlcomprthresh);
    keyFile.set_integer ("Exposure", "ShadowCompr",     toneCurve.shcompr);
    Glib::ArrayHandle<double> tcurve = toneCurve.curve;
    keyFile.set_double_list("Exposure", "Curve",        tcurve);

    // save channel mixer
    Glib::ArrayHandle<int> rmix (chmixer.red, 3, Glib::OWNERSHIP_NONE);
    Glib::ArrayHandle<int> gmix (chmixer.green, 3, Glib::OWNERSHIP_NONE);
    Glib::ArrayHandle<int> bmix (chmixer.blue, 3, Glib::OWNERSHIP_NONE);
    keyFile.set_integer_list("Channel Mixer", "Red",   rmix);
    keyFile.set_integer_list("Channel Mixer", "Green", gmix);
    keyFile.set_integer_list("Channel Mixer", "Blue",  bmix);

    // save luma curve
    keyFile.set_integer ("Luminance Curve", "Brightness",      labCurve.brightness);
    keyFile.set_integer ("Luminance Curve", "Contrast",        labCurve.contrast);
	keyFile.set_integer ("Luminance Curve", "Saturation",	   labCurve.saturation);
    Glib::ArrayHandle<double> lcurve = labCurve.lcurve;
	Glib::ArrayHandle<double> acurve = labCurve.acurve;
    Glib::ArrayHandle<double> bcurve = labCurve.bcurve;
    keyFile.set_double_list("Luminance Curve", "LCurve",        lcurve);
	keyFile.set_double_list("Luminance Curve", "aCurve",        acurve);
    keyFile.set_double_list("Luminance Curve", "bCurve",        bcurve);

    // save sharpening
    keyFile.set_boolean ("Sharpening", "Enabled",              sharpening.enabled);
    keyFile.set_string  ("Sharpening", "Method",               sharpening.method);
    keyFile.set_double  ("Sharpening", "Radius",               sharpening.radius);
    keyFile.set_integer ("Sharpening", "Amount",               sharpening.amount);
    keyFile.set_integer ("Sharpening", "Threshold",            sharpening.threshold);
    keyFile.set_boolean ("Sharpening", "OnlyEdges",            sharpening.edgesonly);
    keyFile.set_double  ("Sharpening", "EdgedetectionRadius",  sharpening.edges_radius);
    keyFile.set_integer ("Sharpening", "EdgeTolerance",        sharpening.edges_tolerance);
    keyFile.set_boolean ("Sharpening", "HalocontrolEnabled",   sharpening.halocontrol);
    keyFile.set_integer ("Sharpening", "HalocontrolAmount",    sharpening.halocontrol_amount);
    keyFile.set_double  ("Sharpening", "DeconvRadius",         sharpening.deconvradius);
    keyFile.set_integer ("Sharpening", "DeconvAmount",         sharpening.deconvamount);
    keyFile.set_integer ("Sharpening", "DeconvDamping",        sharpening.deconvdamping);
    keyFile.set_integer ("Sharpening", "DeconvIterations",     sharpening.deconviter);
    
    // save colorBoost
    keyFile.set_integer ("Color Boost", "Amount",              colorBoost.amount);
    keyFile.set_boolean ("Color Boost", "AvoidColorClipping",  colorBoost.avoidclip);
    keyFile.set_boolean ("Color Boost", "SaturationLimiter",   colorBoost.enable_saturationlimiter);
    keyFile.set_double  ("Color Boost", "SaturationLimit",     colorBoost.saturationlimit);

    // save wb
    if (wb.method=="Camera")
        keyFile.set_string  ("White Balance", "Setting",     "Camera");
    else
        keyFile.set_string  ("White Balance", "Setting",     "Custom");
    keyFile.set_integer ("White Balance", "Temperature", wb.temperature);
    keyFile.set_double  ("White Balance", "Green",       wb.green);
    
    // save colorShift
    keyFile.set_double ("Color Shift", "ChannelA", colorShift.a);
    keyFile.set_double ("Color Shift", "ChannelB", colorShift.b);
	
	// save impulseDenoise
    keyFile.set_boolean ("Impulse Denoising", "Enabled",        impulseDenoise.enabled);
	keyFile.set_integer ("Impulse Denoising", "Threshold",        impulseDenoise.thresh);

	// save defringe
    keyFile.set_boolean ("Defringing", "Enabled",        defringe.enabled);
    keyFile.set_double  ("Defringing", "Radius",         defringe.radius);
    keyFile.set_integer ("Defringing", "Threshold",		defringe.threshold);
	
	// save dirpyrDenoise
    keyFile.set_boolean ("Directional Pyramid Denoising", "Enabled", dirpyrDenoise.enabled);
    keyFile.set_integer ("Directional Pyramid Denoising", "Luma",    dirpyrDenoise.luma);
    keyFile.set_integer ("Directional Pyramid Denoising", "Chroma",  dirpyrDenoise.chroma);
	keyFile.set_double	("Directional Pyramid Denoising", "Gamma",  dirpyrDenoise.gamma);

    // save lumaDenoise
    keyFile.set_boolean ("Luminance Denoising", "Enabled",        lumaDenoise.enabled);
    keyFile.set_double  ("Luminance Denoising", "Radius",         lumaDenoise.radius);
    keyFile.set_integer ("Luminance Denoising", "EdgeTolerance",  lumaDenoise.edgetolerance);

    // save colorDenoise
    keyFile.set_boolean ("Chrominance Denoising", "Enabled",        colorDenoise.enabled);
    keyFile.set_integer ("Chrominance Denoising", "Amount",  		colorDenoise.amount);

    // save sh
    keyFile.set_boolean ("Shadows & Highlights", "Enabled",               sh.enabled);
    keyFile.set_boolean ("Shadows & Highlights", "HighQuality",           sh.hq);
    keyFile.set_integer ("Shadows & Highlights", "Highlights",            sh.highlights);
    keyFile.set_integer ("Shadows & Highlights", "HighlightTonalWidth",   sh.htonalwidth);
    keyFile.set_integer ("Shadows & Highlights", "Shadows",               sh.shadows);
    keyFile.set_integer ("Shadows & Highlights", "ShadowTonalWidth",      sh.stonalwidth);
    keyFile.set_integer ("Shadows & Highlights", "LocalContrast",         sh.localcontrast);
    keyFile.set_integer ("Shadows & Highlights", "Radius",                sh.radius);

    // save crop
    keyFile.set_boolean ("Crop", "Enabled",     crop.enabled);
    keyFile.set_integer ("Crop", "X",           crop.x);
    keyFile.set_integer ("Crop", "Y",           crop.y);
    keyFile.set_integer ("Crop", "W",           crop.w);
    keyFile.set_integer ("Crop", "H",           crop.h);
    keyFile.set_boolean ("Crop", "FixedRatio",  crop.fixratio);
    keyFile.set_string  ("Crop", "Ratio",       crop.ratio);
    keyFile.set_string  ("Crop", "Orientation", crop.orientation);
    keyFile.set_string  ("Crop", "Guide",       crop.guide);
    
    // save coarse
    keyFile.set_integer ("Coarse Transformation", "Rotate",          coarse.rotate);
    keyFile.set_boolean ("Coarse Transformation", "HorizontalFlip",  coarse.hflip);
    keyFile.set_boolean ("Coarse Transformation", "VerticalFlip",    coarse.vflip);
    
    // save commonTrans
    keyFile.set_boolean ("Common Properties for Transformations", "AutoFill", commonTrans.autofill);

    // save rotate
    keyFile.set_double  ("Rotation", "Degree", rotate.degree);

    // save distortion
    keyFile.set_double  ("Distortion", "Amount", distortion.amount);
    keyFile.set_boolean ("Distortion", "UseLensFun", distortion.uselensfun);

    // save perspective correction
    keyFile.set_integer  ("Perspective", "Horizontal", perspective.horizontal);
    keyFile.set_integer  ("Perspective", "Vertical",   perspective.vertical);

    // save C/A correction
    keyFile.set_double  ("CACorrection", "Red",  cacorrection.red);
    keyFile.set_double  ("CACorrection", "Blue", cacorrection.blue);

    // save vignetting correction
    keyFile.set_integer ("Vignetting Correction", "Amount", vignetting.amount);
    keyFile.set_integer ("Vignetting Correction", "Radius", vignetting.radius);
    keyFile.set_integer ("Vignetting Correction", "Strength", vignetting.strength);
    keyFile.set_integer ("Vignetting Correction", "CenterX", vignetting.centerX);
    keyFile.set_integer ("Vignetting Correction", "CenterY", vignetting.centerY);

    // save highlight recovery settings
    keyFile.set_boolean ("HLRecovery", "Enabled",  hlrecovery.enabled);
    keyFile.set_string  ("HLRecovery", "Method",   hlrecovery.method);

    keyFile.set_boolean ("Resize", "Enabled",resize.enabled);
    keyFile.set_double  ("Resize", "Scale",  resize.scale);
    keyFile.set_string  ("Resize", "AppliesTo", resize.appliesTo);
    keyFile.set_string  ("Resize", "Method", resize.method);
    keyFile.set_integer ("Resize", "DataSpecified",  resize.dataspec);
    keyFile.set_integer ("Resize", "Width",  resize.width);
    keyFile.set_integer ("Resize", "Height", resize.height);

    // save color management settings
    keyFile.set_string  ("Color Management", "InputProfile",   icm.input);
    keyFile.set_boolean ("Color Management", "ApplyGammaBeforeInputProfile",   icm.gammaOnInput);
    keyFile.set_string  ("Color Management", "WorkingProfile", icm.working);
    keyFile.set_string  ("Color Management", "OutputProfile",  icm.output);
    
    // save wavelet equalizer parameters
    keyFile.set_boolean ("Equalizer", "Enabled", equalizer.enabled);
    for(int i = 0; i < 8; i++)
    {
        std::stringstream ss;
        ss << "C" << i;
        keyFile.set_integer("Equalizer", ss.str(), equalizer.c[i]);
    }
	
	// save directional pyramid equalizer parameters
    keyFile.set_boolean ("Directional Pyramid Equalizer", "Enabled", dirpyrequalizer.enabled);
    for(int i = 0; i < 5; i++)
    {
        std::stringstream ss;
        ss << "Mult" << i;
        keyFile.set_double("Directional Pyramid Equalizer", ss.str(), dirpyrequalizer.mult[i]);
    }
	
	// save hsv equalizer parameters
    keyFile.set_boolean ("HSV Equalizer", "Enabled", hsvequalizer.enabled);
	keyFile.set_string  ("HSV Equalizer", "Channel", hsvequalizer.hsvchannel);
    for(int i = 0; i < 8; i++)
    {
        std::stringstream ss;
        ss << "Sat" << i;
        keyFile.set_double("HSV Equalizer", ss.str(), hsvequalizer.sat[i]);
    }
	for(int i = 0; i < 8; i++)
    {
        std::stringstream ss;
        ss << "Val" << i;
        keyFile.set_double("HSV Equalizer", ss.str(), hsvequalizer.val[i]);
    }
	for(int i = 0; i < 8; i++)
    {
        std::stringstream ss;
        ss << "Hue" << i;
        keyFile.set_double("HSV Equalizer", ss.str(), hsvequalizer.hue[i]);
    }

    // save RAW parameters
    keyFile.set_string  ("RAW", "DarkFrame", raw.dark_frame );
    keyFile.set_boolean ("RAW", "DarkFrameAuto", raw.df_autoselect );
    keyFile.set_string  ("RAW", "FlatFieldFile", raw.ff_file );
    keyFile.set_boolean ("RAW", "FlatFieldAutoSelect", raw.ff_AutoSelect );
    keyFile.set_integer ("RAW", "FlatFieldBlurRadius", raw.ff_BlurRadius );
    keyFile.set_string  ("RAW", "FlatFieldBlurType", raw.ff_BlurType );     
    keyFile.set_boolean ("RAW", "CA", raw.ca_autocorrect );
	keyFile.set_double	("RAW", "CARed", raw.cared );
    keyFile.set_double	("RAW", "CABlue", raw.cablue );
	keyFile.set_boolean ("RAW", "HotDeadPixels", raw.hotdeadpix_filt );
	keyFile.set_integer ("RAW", "HotDeadPixelThresh", raw.hotdeadpix_thresh );
    keyFile.set_integer ("RAW", "LineDenoise", raw.linenoise);
    keyFile.set_integer ("RAW", "GreenEqThreshold", raw.greenthresh);
    keyFile.set_integer ("RAW", "CcSteps", raw.ccSteps);
    keyFile.set_string  ("RAW", "Method", raw.dmethod );
    keyFile.set_integer ("RAW", "DCBIterations", raw.dcb_iterations );
    keyFile.set_boolean ("RAW", "DCBEnhance", raw.dcb_enhance );

    // save exif change list
    for (int i=0; i<(int)exif.size(); i++)
        keyFile.set_string ("Exif", exif[i].field, exif[i].value);

    // save iptc change list
    for (int i=0; i<(int)iptc.size(); i++) {
        Glib::ArrayHandle<Glib::ustring> values = iptc[i].values;
        keyFile.set_string_list ("IPTC", iptc[i].field, values);
    }
    
    FILE *f = safe_g_fopen (fname, "wt");
    
    if (f==NULL)
        return 1;
    else {
        fprintf (f, "%s", keyFile.to_data().c_str());
        fclose (f);
        return 0;
    }
}

int ProcParams::load (Glib::ustring fname) {

    SafeKeyFile keyFile;
    try {
        setDefaults ();

        FILE* f = safe_g_fopen (fname, "rt");
        if (!f)
            return 1;
        char* buffer = new char[1024];
        std::ostringstream ostr;
        while (fgets (buffer, 1024, f))
            ostr << buffer << "\n";
        delete [] buffer;
        if (!keyFile.load_from_data (ostr.str())) 
            return 1;
        fclose (f);

    // load tonecurve:

version = 200;
if (keyFile.has_group ("Version")) {    
    if (keyFile.has_key ("Version", "Version")) version = keyFile.get_integer ("Version", "Version");
}

if (keyFile.has_group ("Exposure")) {    
    if (keyFile.has_key ("Exposure", "Auto"))           toneCurve.autoexp       = keyFile.get_boolean ("Exposure", "Auto");
    if (keyFile.has_key ("Exposure", "Clip"))           toneCurve.clip          = keyFile.get_double  ("Exposure", "Clip");
    if (keyFile.has_key ("Exposure", "Compensation"))   toneCurve.expcomp       = keyFile.get_double  ("Exposure", "Compensation");
    if (keyFile.has_key ("Exposure", "Brightness"))     toneCurve.brightness    = keyFile.get_integer ("Exposure", "Brightness");
    if (keyFile.has_key ("Exposure", "Contrast"))       toneCurve.contrast      = keyFile.get_integer ("Exposure", "Contrast");
	if (keyFile.has_key ("Exposure", "Saturation"))     toneCurve.saturation    = keyFile.get_integer ("Exposure", "Saturation");
	if (keyFile.has_key ("Exposure", "Black"))          toneCurve.black         = keyFile.get_integer ("Exposure", "Black");
    if (keyFile.has_key ("Exposure", "HighlightCompr")) toneCurve.hlcompr       = keyFile.get_integer ("Exposure", "HighlightCompr");
    if (toneCurve.hlcompr > 100) toneCurve.hlcompr = 100; // older pp3 files can have values above 100.
    if (keyFile.has_key ("Exposure", "HighlightComprThreshold")) toneCurve.hlcomprthresh = keyFile.get_integer ("Exposure", "HighlightComprThreshold");
    if (keyFile.has_key ("Exposure", "ShadowCompr"))    toneCurve.shcompr       = keyFile.get_integer ("Exposure", "ShadowCompr");
    if (toneCurve.shcompr > 100) toneCurve.shcompr = 100; // older pp3 files can have values above 100.
    if (version>200)
	if (keyFile.has_key ("Exposure", "Curve"))          toneCurve.curve         = keyFile.get_double_list ("Exposure", "Curve");
}

    // load channel mixer curve
if (keyFile.has_group ("Channel Mixer")) {    
    if (keyFile.has_key ("Channel Mixer", "Red") && keyFile.has_key ("Channel Mixer", "Green") && keyFile.has_key ("Channel Mixer", "Blue")) {
        Glib::ArrayHandle<int> rmix = keyFile.get_integer_list ("Channel Mixer", "Red");
        Glib::ArrayHandle<int> gmix = keyFile.get_integer_list ("Channel Mixer", "Green");
        Glib::ArrayHandle<int> bmix = keyFile.get_integer_list ("Channel Mixer", "Blue");
        memcpy (chmixer.red, rmix.data(), 3*sizeof(int));
        memcpy (chmixer.green, gmix.data(), 3*sizeof(int));
        memcpy (chmixer.blue, bmix.data(), 3*sizeof(int));
    }
}

    // load luma curve
if (keyFile.has_group ("Luminance Curve")) {    
    if (keyFile.has_key ("Luminance Curve", "Brightness"))     labCurve.brightness = keyFile.get_integer  ("Luminance Curve", "Brightness");
    if (keyFile.has_key ("Luminance Curve", "Contrast"))       labCurve.contrast   = keyFile.get_integer ("Luminance Curve", "Contrast");
	if (keyFile.has_key ("Luminance Curve", "Saturation"))      labCurve.saturation   = keyFile.get_integer ("Luminance Curve", "Saturation");
    if (version>200)
	if (keyFile.has_key ("Luminance Curve", "LCurve"))          labCurve.lcurve      = keyFile.get_double_list ("Luminance Curve", "LCurve");
	if (keyFile.has_key ("Luminance Curve", "aCurve"))          labCurve.acurve      = keyFile.get_double_list ("Luminance Curve", "aCurve");
	if (keyFile.has_key ("Luminance Curve", "bCurve"))          labCurve.bcurve      = keyFile.get_double_list ("Luminance Curve", "bCurve");
}

    // load sharpening
if (keyFile.has_group ("Sharpening")) {    
    if (keyFile.has_key ("Sharpening", "Enabled"))              sharpening.enabled          = keyFile.get_boolean ("Sharpening", "Enabled");
    if (keyFile.has_key ("Sharpening", "Radius"))               sharpening.radius           = keyFile.get_double  ("Sharpening", "Radius");
    if (keyFile.has_key ("Sharpening", "Amount"))               sharpening.amount           = keyFile.get_integer ("Sharpening", "Amount");
    if (keyFile.has_key ("Sharpening", "Threshold"))            sharpening.threshold        = keyFile.get_integer ("Sharpening", "Threshold");
    if (keyFile.has_key ("Sharpening", "OnlyEdges"))            sharpening.edgesonly        = keyFile.get_boolean ("Sharpening", "OnlyEdges");
    if (keyFile.has_key ("Sharpening", "EdgedetectionRadius"))  sharpening.edges_radius     = keyFile.get_double  ("Sharpening", "EdgedetectionRadius");
    if (keyFile.has_key ("Sharpening", "EdgeTolerance"))        sharpening.edges_tolerance  = keyFile.get_integer ("Sharpening", "EdgeTolerance");
    if (keyFile.has_key ("Sharpening", "HalocontrolEnabled"))   sharpening.halocontrol      = keyFile.get_boolean ("Sharpening", "HalocontrolEnabled");
    if (keyFile.has_key ("Sharpening", "HalocontrolAmount"))    sharpening.halocontrol_amount = keyFile.get_integer ("Sharpening", "HalocontrolAmount");
    if (keyFile.has_key ("Sharpening", "Method"))               sharpening.method           = keyFile.get_string  ("Sharpening", "Method");
    if (keyFile.has_key ("Sharpening", "DeconvRadius"))         sharpening.deconvradius     = keyFile.get_double  ("Sharpening", "DeconvRadius");
    if (keyFile.has_key ("Sharpening", "DeconvAmount"))         sharpening.deconvamount     = keyFile.get_integer ("Sharpening", "DeconvAmount");
    if (keyFile.has_key ("Sharpening", "DeconvDamping"))        sharpening.deconvdamping    = keyFile.get_integer ("Sharpening", "DeconvDamping");
    if (keyFile.has_key ("Sharpening", "DeconvIterations"))     sharpening.deconviter       = keyFile.get_integer ("Sharpening", "DeconvIterations");
}
    
    // load colorBoost
if (keyFile.has_group ("Color Boost")) {    
    if (keyFile.has_key ("Color Boost", "Amount"))              colorBoost.amount           = keyFile.get_integer ("Color Boost", "Amount");
    else {
        int a=0, b=0;
        if (keyFile.has_key ("Color Boost", "ChannelA"))        a                           = keyFile.get_integer ("Color Boost", "ChannelA");
        if (keyFile.has_key ("Color Boost", "ChannelB"))        b                           = keyFile.get_integer ("Color Boost", "ChannelB");
        colorBoost.amount = (a+b) / 2;
    }   
    if (keyFile.has_key ("Color Boost", "AvoidColorClipping"))  colorBoost.avoidclip               = keyFile.get_boolean ("Color Boost", "AvoidColorClipping");
    if (keyFile.has_key ("Color Boost", "SaturationLimiter"))   colorBoost.enable_saturationlimiter= keyFile.get_boolean ("Color Boost", "SaturationLimiter");
    if (keyFile.has_key ("Color Boost", "SaturationLimit"))     colorBoost.saturationlimit         = keyFile.get_double  ("Color Boost", "SaturationLimit");
}

    // load wb
if (keyFile.has_group ("White Balance")) {    
    if (keyFile.has_key ("White Balance", "Setting"))     wb.method         = keyFile.get_string ("White Balance", "Setting");
    if (keyFile.has_key ("White Balance", "Temperature")) wb.temperature    = keyFile.get_integer ("White Balance", "Temperature");
    if (keyFile.has_key ("White Balance", "Green"))       wb.green          = keyFile.get_double ("White Balance", "Green");
}

    // load colorShift
if (keyFile.has_group ("Color Shift")) {    
    if (keyFile.has_key ("Color Shift", "ChannelA")) colorShift.a = keyFile.get_double ("Color Shift", "ChannelA");
    if (keyFile.has_key ("Color Shift", "ChannelB")) colorShift.b = keyFile.get_double ("Color Shift", "ChannelB");
}
		
// load defringe
if (keyFile.has_group ("Defringing")) {    
	if (keyFile.has_key ("Defringing", "Enabled"))        defringe.enabled       = keyFile.get_boolean ("Defringing", "Enabled");
	if (keyFile.has_key ("Defringing", "Radius"))         defringe.radius        = keyFile.get_double  ("Defringing", "Radius");
	if (keyFile.has_key ("Defringing", "Threshold"))  defringe.threshold = keyFile.get_integer ("Defringing", "Threshold");
}
		
	// load impulseDenoise
if (keyFile.has_group ("Impulse Denoising")) {    
	if (keyFile.has_key ("Impulse Denoising", "Enabled")) impulseDenoise.enabled = keyFile.get_boolean ("Impulse Denoising", "Enabled");
	if (keyFile.has_key ("Impulse Denoising", "Threshold")) impulseDenoise.thresh = keyFile.get_integer ("Impulse Denoising", "Threshold");
}
		
	// load dirpyrDenoise
if (keyFile.has_group ("Directional Pyramid Denoising")) {    
	if (keyFile.has_key ("Directional Pyramid Denoising", "Enabled")) dirpyrDenoise.enabled = keyFile.get_boolean ("Directional Pyramid Denoising", "Enabled");
	if (keyFile.has_key ("Directional Pyramid Denoising", "Luma"))    dirpyrDenoise.luma    = keyFile.get_integer ("Directional Pyramid Denoising", "Luma");
	if (keyFile.has_key ("Directional Pyramid Denoising", "Chroma"))  dirpyrDenoise.chroma  = keyFile.get_integer ("Directional Pyramid Denoising", "Chroma");
	if (keyFile.has_key ("Directional Pyramid Denoising", "Gamma"))  dirpyrDenoise.gamma  = keyFile.get_double ("Directional Pyramid Denoising", "Gamma");
}
  
    // load lumaDenoise
if (keyFile.has_group ("Luminance Denoising")) {    
    if (keyFile.has_key ("Luminance Denoising", "Enabled"))        lumaDenoise.enabled       = keyFile.get_boolean ("Luminance Denoising", "Enabled");
    if (keyFile.has_key ("Luminance Denoising", "Radius"))         lumaDenoise.radius        = keyFile.get_double  ("Luminance Denoising", "Radius");
    if (keyFile.has_key ("Luminance Denoising", "EdgeTolerance"))  lumaDenoise.edgetolerance = keyFile.get_integer ("Luminance Denoising", "EdgeTolerance");
}

    // load colorDenoise
if (keyFile.has_group ("Chrominance Denoising")) {    
    if (keyFile.has_key ("Chrominance Denoising", "Enabled"))        colorDenoise.enabled       = keyFile.get_boolean 	("Chrominance Denoising", "Enabled");
    if (keyFile.has_key ("Chrominance Denoising", "Radius"))         colorDenoise.amount        = 10*keyFile.get_double ("Chrominance Denoising", "Radius");
	else if (keyFile.has_key ("Chrominance Denoising", "Amount"))    colorDenoise.amount        = keyFile.get_integer  	("Chrominance Denoising", "Amount");
}

    // load sh
if (keyFile.has_group ("Shadows & Highlights")) {    
    if (keyFile.has_key ("Shadows & Highlights", "Enabled"))               sh.enabled       = keyFile.get_boolean ("Shadows & Highlights", "Enabled");
    if (keyFile.has_key ("Shadows & Highlights", "HighQuality"))           sh.hq            = keyFile.get_boolean ("Shadows & Highlights", "HighQuality");
    if (keyFile.has_key ("Shadows & Highlights", "Highlights"))            sh.highlights    = keyFile.get_integer ("Shadows & Highlights", "Highlights");
    if (keyFile.has_key ("Shadows & Highlights", "HighlightTonalWidth"))   sh.htonalwidth   = keyFile.get_integer ("Shadows & Highlights", "HighlightTonalWidth");
    if (keyFile.has_key ("Shadows & Highlights", "Shadows"))               sh.shadows       = keyFile.get_integer ("Shadows & Highlights", "Shadows");
    if (keyFile.has_key ("Shadows & Highlights", "ShadowTonalWidth"))      sh.stonalwidth   = keyFile.get_integer ("Shadows & Highlights", "ShadowTonalWidth");
    if (keyFile.has_key ("Shadows & Highlights", "LocalContrast"))         sh.localcontrast = keyFile.get_integer ("Shadows & Highlights", "LocalContrast");
    if (keyFile.has_key ("Shadows & Highlights", "Radius"))                sh.radius        = keyFile.get_integer ("Shadows & Highlights", "Radius");
}
    
    // load crop
if (keyFile.has_group ("Crop")) {    
    if (keyFile.has_key ("Crop", "Enabled"))    crop.enabled    = keyFile.get_boolean ("Crop", "Enabled");
    if (keyFile.has_key ("Crop", "X"))          crop.x          = keyFile.get_integer ("Crop", "X");
    if (keyFile.has_key ("Crop", "Y"))          crop.y          = keyFile.get_integer ("Crop", "Y");
    if (keyFile.has_key ("Crop", "W"))          crop.w          = keyFile.get_integer ("Crop", "W");
    if (keyFile.has_key ("Crop", "H"))          crop.h          = keyFile.get_integer ("Crop", "H");
    if (keyFile.has_key ("Crop", "FixedRatio")) crop.fixratio   = keyFile.get_boolean ("Crop", "FixedRatio");
    if (keyFile.has_key ("Crop", "Ratio"))      crop.ratio      = keyFile.get_string  ("Crop", "Ratio");
    if (keyFile.has_key ("Crop", "Orientation"))crop.orientation= keyFile.get_string  ("Crop", "Orientation");
    if (keyFile.has_key ("Crop", "Guide"))      crop.guide      = keyFile.get_string  ("Crop", "Guide");
}

    // load coarse
if (keyFile.has_group ("Coarse Transformation")) {    
    if (keyFile.has_key ("Coarse Transformation", "Rotate"))          coarse.rotate = keyFile.get_integer ("Coarse Transformation", "Rotate");
    if (keyFile.has_key ("Coarse Transformation", "HorizontalFlip"))  coarse.hflip  = keyFile.get_boolean ("Coarse Transformation", "HorizontalFlip");
    if (keyFile.has_key ("Coarse Transformation", "VerticalFlip"))    coarse.vflip  = keyFile.get_boolean ("Coarse Transformation", "VerticalFlip");
}

    // load rotate
if (keyFile.has_group ("Rotation")) {    
    if (keyFile.has_key ("Rotation", "Degree"))   rotate.degree = keyFile.get_double ("Rotation", "Degree");
}
	// load commonTrans
if (keyFile.has_group ("Common Properties for Transformations")) {
    if (keyFile.has_key ("Common Properties for Transformations", "AutoFill"))   commonTrans.autofill = keyFile.get_boolean ("Common Properties for Transformations", "AutoFill");
}

    // load distortion
if (keyFile.has_group ("Distortion")) {    
    if (keyFile.has_key ("Distortion", "Amount"))     distortion.amount     = keyFile.get_double  ("Distortion", "Amount");
    if (keyFile.has_key ("Distortion", "UseLensFun")) distortion.uselensfun = keyFile.get_boolean ("Distortion", "UseLensFun");
}
    
	// load perspective correction
if (keyFile.has_group ("Perspective")) {
	if (keyFile.has_key ("Perspective", "Horizontal")) 	perspective.horizontal 	= keyFile.get_integer ("Perspective", "Horizontal");
	if (keyFile.has_key ("Perspective", "Vertical")) 	perspective.vertical 	= keyFile.get_integer ("Perspective", "Vertical");
}

// load c/a correction
if (keyFile.has_group ("CACorrection")) {    
    if (keyFile.has_key ("CACorrection", "Red"))  cacorrection.red  = keyFile.get_double ("CACorrection", "Red");
    if (keyFile.has_key ("CACorrection", "Blue")) cacorrection.blue = keyFile.get_double ("CACorrection", "Blue");
}

    // load vignetting correction
if (keyFile.has_group ("Vignetting Correction")) {    
    if (keyFile.has_key ("Vignetting Correction", "Amount")) vignetting.amount = keyFile.get_integer ("Vignetting Correction", "Amount");
    if (keyFile.has_key ("Vignetting Correction", "Radius")) vignetting.radius = keyFile.get_integer ("Vignetting Correction", "Radius");
    if (keyFile.has_key ("Vignetting Correction", "Strength")) vignetting.strength = keyFile.get_integer ("Vignetting Correction", "Strength");
    if (keyFile.has_key ("Vignetting Correction", "CenterX")) vignetting.centerX = keyFile.get_integer ("Vignetting Correction", "CenterX");
    if (keyFile.has_key ("Vignetting Correction", "CenterY")) vignetting.centerY = keyFile.get_integer ("Vignetting Correction", "CenterY");
}

    // load highlight recovery settings
if (keyFile.has_group ("HLRecovery")) {    
    if (keyFile.has_key ("HLRecovery", "Enabled"))  hlrecovery.enabled  = keyFile.get_boolean ("HLRecovery", "Enabled");
    if (keyFile.has_key ("HLRecovery", "Method"))   hlrecovery.method   = keyFile.get_string  ("HLRecovery", "Method");
}
    // load resize settings
if (keyFile.has_group ("Resize")) {    
    if (keyFile.has_key ("Resize", "Enabled")) resize.enabled = keyFile.get_boolean ("Resize", "Enabled");
    if (keyFile.has_key ("Resize", "Scale"))  resize.scale  = keyFile.get_double ("Resize", "Scale");
    if (keyFile.has_key ("Resize", "AppliesTo")) resize.appliesTo = keyFile.get_string ("Resize", "AppliesTo");
    if (keyFile.has_key ("Resize", "Method")) resize.method = keyFile.get_string ("Resize", "Method");
    if (keyFile.has_key ("Resize", "DataSpecified")) resize.dataspec = keyFile.get_integer ("Resize", "DataSpecified");
    if (keyFile.has_key ("Resize", "Width"))  resize.width  = keyFile.get_integer ("Resize", "Width");
    if (keyFile.has_key ("Resize", "Height")) resize.height = keyFile.get_integer ("Resize", "Height");
}

    // load color management settings
if (keyFile.has_group ("Color Management")) {    
    if (keyFile.has_key ("Color Management", "InputProfile"))   icm.input   = keyFile.get_string ("Color Management", "InputProfile");
    if (keyFile.has_key ("Color Management", "ApplyGammaBeforeInputProfile"))   icm.gammaOnInput   = keyFile.get_boolean ("Color Management", "ApplyGammaBeforeInputProfile");
    if (keyFile.has_key ("Color Management", "WorkingProfile")) icm.working = keyFile.get_string ("Color Management", "WorkingProfile");
    if (keyFile.has_key ("Color Management", "OutputProfile"))  icm.output  = keyFile.get_string ("Color Management", "OutputProfile");
}

    // load wavelet equalizer parameters
if (keyFile.has_group ("Equalizer")) {
    if (keyFile.has_key ("Equalizer", "Enabled")) equalizer.enabled = keyFile.get_boolean ("Equalizer", "Enabled");
    for(int i = 0; i < 8; i ++)
    {
        std::stringstream ss;
        ss << "C" << i;
        if(keyFile.has_key ("Equalizer", ss.str())) equalizer.c[i] = keyFile.get_integer ("Equalizer", ss.str());
    }
}
		
	// load directional pyramid equalizer parameters
if (keyFile.has_group ("Directional Pyramid Equalizer")) {
	if (keyFile.has_key ("Directional Pyramid Equalizer", "Enabled")) dirpyrequalizer.enabled = keyFile.get_boolean ("Directional Pyramid Equalizer", "Enabled");
	for(int i = 0; i < 5; i ++)
	{
		std::stringstream ss;
		ss << "Mult" << i;
		if(keyFile.has_key ("Directional Pyramid Equalizer", ss.str())) dirpyrequalizer.mult[i] = keyFile.get_double ("Directional Pyramid Equalizer", ss.str());
	}
}
		
	// load wavelet equalizer parameters
if (keyFile.has_group ("HSV Equalizer")) {
	if (keyFile.has_key ("HSV Equalizer", "Enabled")) hsvequalizer.enabled = keyFile.get_boolean ("HSV Equalizer", "Enabled");
	for(int i = 0; i < 8; i ++)
	{
		std::stringstream ss;
		ss << "Sat" << i;
		if(keyFile.has_key ("HSV Equalizer", ss.str())) hsvequalizer.sat[i] = keyFile.get_double ("HSV Equalizer", ss.str());
	}
	for(int i = 0; i < 8; i ++)
	{
		std::stringstream ss;
		ss << "Val" << i;
		if(keyFile.has_key ("HSV Equalizer", ss.str())) hsvequalizer.val[i] = keyFile.get_double ("HSV Equalizer", ss.str());
	}
	for(int i = 0; i < 8; i ++)
	{
		std::stringstream ss;
		ss << "Hue" << i;
		if(keyFile.has_key ("HSV Equalizer", ss.str())) hsvequalizer.hue[i] = keyFile.get_double ("HSV Equalizer", ss.str());
	}
}

	// load raw settings
if (keyFile.has_group ("RAW")) {
	if (keyFile.has_key ("RAW", "DarkFrame"))     raw.dark_frame = keyFile.get_string  ("RAW", "DarkFrame" );
	if (keyFile.has_key ("RAW", "DarkFrameAuto")) raw.df_autoselect = keyFile.get_boolean ("RAW", "DarkFrameAuto" );
	if (keyFile.has_key ("RAW", "FlatFieldFile"))       raw.ff_file = keyFile.get_string  ("RAW", "FlatFieldFile" );                    
	if (keyFile.has_key ("RAW", "FlatFieldAutoSelect")) raw.ff_AutoSelect = keyFile.get_boolean  ("RAW", "FlatFieldAutoSelect" ); 
	if (keyFile.has_key ("RAW", "FlatFieldBlurRadius")) raw.ff_BlurRadius = keyFile.get_integer  ("RAW", "FlatFieldBlurRadius" );
	if (keyFile.has_key ("RAW", "FlatFieldBlurType"))   raw.ff_BlurType = keyFile.get_string  ("RAW", "FlatFieldBlurType" );		
	if (keyFile.has_key ("RAW", "CA"))            raw.ca_autocorrect = keyFile.get_boolean ("RAW", "CA" );
	if (keyFile.has_key ("RAW", "CARed"))         raw.cared = keyFile.get_double ("RAW", "CARed" );
	if (keyFile.has_key ("RAW", "CABlue"))        raw.cablue = keyFile.get_double ("RAW", "CABlue" );
	if (keyFile.has_key ("RAW", "HotDeadPixels")) raw.hotdeadpix_filt = keyFile.get_boolean ("RAW", "HotDeadPixels" );
	if (keyFile.has_key ("RAW", "HotDeadPixelThresh")) raw.hotdeadpix_thresh = keyFile.get_integer ("RAW", "HotDeadPixelThresh" );
	if (keyFile.has_key ("RAW", "LineDenoise"))   raw.linenoise = keyFile.get_integer ("RAW", "LineDenoise" );
	if (keyFile.has_key ("RAW", "GreenEqThreshold")) raw.greenthresh= keyFile.get_integer ("RAW", "GreenEqThreshold");
	if (keyFile.has_key ("RAW", "CcSteps"))       raw.ccSteps  = keyFile.get_integer ("RAW", "CcSteps");
	if (keyFile.has_key ("RAW", "Method"))        raw.dmethod = keyFile.get_string ("RAW", "Method");
	if (keyFile.has_key ("RAW", "DCBIterations")) raw.dcb_iterations = keyFile.get_integer("RAW", "DCBIterations");
	if (keyFile.has_key ("RAW", "DCBEnhance"))    raw.dcb_enhance =keyFile.get_boolean("RAW", "DCBEnhance");
}

    // load exif change settings
if (keyFile.has_group ("Exif")) {
    std::vector<Glib::ustring> keys = keyFile.get_keys ("Exif");
    exif.resize (keys.size());
    for (int i=0; i<(int)keys.size(); i++) {
        exif[i].field = keys[i];
        exif[i].value = keyFile.get_string ("Exif", keys[i]);
    }
}

    // load iptc change settings
if (keyFile.has_group ("IPTC")) {
    std::vector<Glib::ustring> keys = keyFile.get_keys ("IPTC");
    iptc.resize (keys.size());
    for (int i=0; i<(int)keys.size(); i++) {
        iptc[i].field = keys[i];
        iptc[i].values = keyFile.get_string_list ("IPTC", keys[i]);
    }
}


        return 0;
    }
    catch (const Glib::Error& e) {
        printf ("-->%s\n", e.what().c_str());
        return 1;
    }
    catch (...) {
        printf ("-->unknown exception!\n");
        return 1;
    }
}

bool operator==(const EqualizerParams & a, const EqualizerParams & b) {
    if(a.enabled != b.enabled)
        return false;

    for(int i = 0; i < 8; i++) {
        if(a.c[i] != b.c[i])
            return false;
    }
    return true;
}
	
bool operator==(const DirPyrEqualizerParams & a, const DirPyrEqualizerParams & b) {
	if(a.enabled != b.enabled)
		return false;
		
	for(int i = 0; i < 5; i++) {
		if(a.mult[i] != b.mult[i])
			return false;
	}
	return true;
}

bool operator==(const HSVEqualizerParams & a, const HSVEqualizerParams & b) {
	if(a.enabled != b.enabled)
		return false;
	
	for(int i = 0; i < 8; i++) {
		if(a.sat[i] != b.sat[i] && a.val[i] != b.val[i] && a.hue[i] != b.hue[i])
			return false;
	}
	return true;
}
	
bool operator==(const ExifPair& a, const ExifPair& b) {

    return a.field == b.field && a.value == b.value;
}

bool operator==(const IPTCPair& a, const IPTCPair& b) {

    return a.field == b.field && a.values == b.values;
}
bool ProcParams::operator== (const ProcParams& other) {

	return
		toneCurve.curve == other.toneCurve.curve
		&& toneCurve.brightness == other.toneCurve.brightness
		&& toneCurve.black == other.toneCurve.black
		&& toneCurve.contrast == other.toneCurve.contrast
		&& toneCurve.saturation == other.toneCurve.saturation
		&& toneCurve.shcompr == other.toneCurve.shcompr
		&& toneCurve.hlcompr == other.toneCurve.hlcompr
		&& toneCurve.hlcomprthresh == other.toneCurve.hlcomprthresh
		&& toneCurve.autoexp == other.toneCurve.autoexp
		&& toneCurve.clip == other.toneCurve.clip
		&& toneCurve.expcomp == other.toneCurve.expcomp
		&& labCurve.lcurve == other.labCurve.lcurve
		&& labCurve.acurve == other.labCurve.acurve
		&& labCurve.bcurve == other.labCurve.bcurve
		&& labCurve.brightness == other.labCurve.brightness
		&& labCurve.contrast == other.labCurve.contrast
		&& labCurve.saturation == other.labCurve.saturation
		&& sharpening.enabled == other.sharpening.enabled
		&& sharpening.radius == other.sharpening.radius
		&& sharpening.amount == other.sharpening.amount
		&& sharpening.threshold == other.sharpening.threshold
		&& sharpening.edgesonly == other.sharpening.edgesonly
		&& sharpening.edges_radius == other.sharpening.edges_radius
		&& sharpening.edges_tolerance == other.sharpening.edges_tolerance
		&& sharpening.halocontrol == other.sharpening.halocontrol
		&& sharpening.halocontrol_amount== other.sharpening.halocontrol_amount
		&& sharpening.method == other.sharpening.method
		&& sharpening.deconvamount == other.sharpening.deconvamount
		&& sharpening.deconvradius == other.sharpening.deconvradius
		&& sharpening.deconviter == other.sharpening.deconviter
		&& sharpening.deconvdamping == other.sharpening.deconvdamping
		&& colorBoost.amount == other.colorBoost.amount
		&& colorBoost.avoidclip == other.colorBoost.avoidclip
		&& colorBoost.enable_saturationlimiter == other.colorBoost.enable_saturationlimiter
		&& colorBoost.saturationlimit == other.colorBoost.saturationlimit
		&& wb.method == other.wb.method
		&& wb.green == other.wb.green
		&& wb.temperature == other.wb.temperature
		&& colorShift.a == other.colorShift.a
		&& colorShift.b == other.colorShift.b
		&& impulseDenoise.enabled == other.impulseDenoise.enabled
		&& impulseDenoise.thresh == other.impulseDenoise.thresh
		&& dirpyrDenoise.enabled == other.dirpyrDenoise.enabled
		&& dirpyrDenoise.luma == other.dirpyrDenoise.luma
		&& dirpyrDenoise.chroma == other.dirpyrDenoise.chroma
		&& dirpyrDenoise.gamma == other.dirpyrDenoise.gamma
		&& defringe.enabled == other.defringe.enabled
		&& defringe.radius == other.defringe.radius
		&& defringe.threshold == other.defringe.threshold
		&& lumaDenoise.enabled == other.lumaDenoise.enabled
		&& lumaDenoise.radius == other.lumaDenoise.radius
		&& lumaDenoise.edgetolerance == other.lumaDenoise.edgetolerance
		&& colorDenoise.enabled == other.colorDenoise.enabled
		&& colorDenoise.radius == other.colorDenoise.radius
		&& colorDenoise.edgetolerance == other.colorDenoise.edgetolerance
		&& colorDenoise.edgesensitive == other.colorDenoise.edgesensitive
		&& sh.enabled == other.sh.enabled
		&& sh.hq == other.sh.hq
		&& sh.highlights == other.sh.highlights
		&& sh.htonalwidth == other.sh.htonalwidth
		&& sh.shadows == other.sh.shadows
		&& sh.stonalwidth == other.sh.stonalwidth
		&& sh.localcontrast == other.sh.localcontrast
		&& sh.radius == other.sh.radius
		&& crop.enabled == other.crop.enabled
		&& crop.x == other.crop.x
		&& crop.y == other.crop.y
		&& crop.w == other.crop.w
		&& crop.h == other.crop.h
		&& crop.fixratio == other.crop.fixratio
		&& crop.ratio == other.crop.ratio
		&& crop.orientation == other.crop.orientation
		&& crop.guide == other.crop.guide
		&& coarse.rotate == other.coarse.rotate
		&& coarse.hflip == other.coarse.hflip
		&& coarse.vflip == other.coarse.vflip
		&& rotate.degree == other.rotate.degree
		&& commonTrans.autofill == other.commonTrans.autofill
		&& distortion.uselensfun == other.distortion.uselensfun
		&& distortion.amount == other.distortion.amount
		&& perspective.horizontal == other.perspective.horizontal
		&& perspective.vertical == other.perspective.vertical
		&& cacorrection.red == other.cacorrection.red
		&& cacorrection.blue == other.cacorrection.blue
		&& vignetting.amount == other.vignetting.amount
		&& vignetting.radius == other.vignetting.radius
		&& vignetting.strength == other.vignetting.strength
		&& vignetting.centerX == other.vignetting.centerX
		&& vignetting.centerY == other.vignetting.centerY
		&& !memcmp (&chmixer.red, &other.chmixer.red, 3*sizeof(int))
		&& !memcmp (&chmixer.green, &other.chmixer.green, 3*sizeof(int))
		&& !memcmp (&chmixer.blue, &other.chmixer.blue, 3*sizeof(int))
		&& hlrecovery.enabled == other.hlrecovery.enabled
		&& hlrecovery.method == other.hlrecovery.method
		&& resize.scale == other.resize.scale
		&& resize.appliesTo == other.resize.appliesTo
		&& resize.method == other.resize.method
		&& resize.dataspec == other.resize.dataspec
		&& resize.width == other.resize.width
		&& resize.height == other.resize.height
		&& raw.dark_frame == other.raw.dark_frame
		&& raw.df_autoselect == other.raw.df_autoselect
		&& raw.ff_file   == other.raw.ff_file
		&& raw.ff_AutoSelect   == other.raw.ff_AutoSelect
		&& raw.ff_BlurRadius   == other.raw.ff_BlurRadius
		&& raw.ff_BlurType   == other.raw.ff_BlurType	
		&& raw.dcb_enhance == other.raw.dcb_enhance
		&& raw.dcb_iterations == other.raw.dcb_iterations
		&& raw.ccSteps == other.raw.ccSteps
		&& raw.ca_autocorrect == other.raw.ca_autocorrect
		&& raw.cared == other.raw.cared
		&& raw.cablue == other.raw.cablue
		&& raw.hotdeadpix_filt == other.raw.hotdeadpix_filt
		&& raw.dmethod == other.raw.dmethod
		&& raw.greenthresh == other.raw.greenthresh
		&& raw.linenoise == other.raw.linenoise
		&& icm.input == other.icm.input
		&& icm.gammaOnInput == other.icm.gammaOnInput
		&& icm.working == other.icm.working
		&& icm.output == other.icm.output
		&& equalizer == other.equalizer
		&& dirpyrequalizer == other.dirpyrequalizer
		&& hsvequalizer == other.hsvequalizer
		&& exif==other.exif
		&& iptc==other.iptc;
	}

bool ProcParams::operator!= (const ProcParams& other) {

    return !(*this==other);
}
}
}

