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
#include "version.h"
#include <ppversion.h>
#include <mydiagonalcurve.h>
#include <myflatcurve.h>

#include <safekeyfile.h>
#include <rawimage.h>
#define APPVERSION VERSION


namespace rtengine {
namespace procparams {

const char *RAWParams::methodstring[RAWParams::numMethods]={"eahd", "hphd", "vng4", "dcb", "amaze", "ahd", "fast" };
const char *RAWParams::ff_BlurTypestring[RAWParams::numFlatFileBlurTypes]={/*"Parametric",*/ "Area Flatfield", "Vertical Flatfield", "Horizontal Flatfield", "V+H Flatfield"};

// Maps crop to resized width (e.g. smaller previews)
void CropParams::mapToResized(int resizedWidth, int resizedHeight, int scale, int &x1, int &x2, int &y1, int &y2) const {
    x1 = 0, x2 = resizedWidth, y1 = 0, y2 = resizedHeight;
    if (enabled) {
        x1 = MIN(resizedWidth-1,  MAX(0, x / scale));
        y1 = MIN(resizedHeight-1, MAX(0, y / scale));   
        x2 = MIN(resizedWidth,    MAX(0, (x+w) / scale)); 
        y2 = MIN(resizedHeight,   MAX(0, (y+h) / scale));
    }
}

ProcParams::ProcParams () { 

    setDefaults (); 
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
    toneCurve.curve.push_back(DCT_Linear);
    
    labCurve.brightness    = 0;
    labCurve.contrast      = 0;
    labCurve.saturation    = 0;
    labCurve.avoidclip          = false;
    labCurve.enable_saturationlimiter = false;
    labCurve.saturationlimit    = 50;
    labCurve.lcurve.clear ();
    labCurve.lcurve.push_back(DCT_Linear);
    labCurve.acurve.clear ();
    labCurve.acurve.push_back(DCT_Linear);
    labCurve.bcurve.clear ();
    labCurve.bcurve.push_back(DCT_Linear);
   
    sharpenEdge.enabled         = false;
    sharpenEdge.passes          = 2;
    sharpenEdge.amount        = 50.0;
    sharpenEdge.threechannels   = false;

    sharpenMicro.enabled        = false;
    sharpenMicro.amount       = 20.0;
    sharpenMicro.uniformity     = 50.0;
    sharpenMicro.matrix         = false;

    sharpening.enabled          = false;
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

    vibrance.enabled            = false;
    vibrance.pastels            = 0;
    vibrance.saturated          = 0;
    vibrance.psthreshold        = 75;
    vibrance.protectskins       = false;
    vibrance.avoidcolorshift    = true;
    vibrance.pastsattog     	= true;
    
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
    impulseDenoise.thresh       = 50;
	
    defringe.enabled            = false;
    defringe.radius             = 2.0;
    defringe.threshold          = 25;

    dirpyrDenoise.enabled       = false;
    dirpyrDenoise.luma          = 10;
    dirpyrDenoise.chroma        = 10;
    dirpyrDenoise.gamma         = 2.0;
    dirpyrDenoise.lumcurve.clear ();
    dirpyrDenoise.lumcurve.push_back (DCT_Linear);
    dirpyrDenoise.chromcurve.clear ();
    dirpyrDenoise.chromcurve.push_back (DCT_Linear);
    
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
    icm.blendCMSMatrix = false;
    icm.working = "sRGB";
    icm.output  = "sRGB";
    icm.gamma  = "default";
 	icm.gampos =2.22;
	icm.slpos=4.5;
    icm.freegamma = false;
  
	dirpyrequalizer.enabled = false;    
    for(int i = 0; i < 4; i ++)
    {
        dirpyrequalizer.mult[i] = 1.0;
    }
    dirpyrequalizer.mult[4] = 0.0;
    hsvequalizer.hcurve.clear ();
    hsvequalizer.hcurve.push_back (FCT_Linear);
    hsvequalizer.scurve.clear ();
    hsvequalizer.scurve.push_back (FCT_Linear);
    hsvequalizer.vcurve.clear ();
    hsvequalizer.vcurve.push_back (FCT_Linear);
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
    raw.all_enhance=false;
	
    // exposure before interpolation
    raw.expos=1.0;
    raw.preser=0.0;
	raw.blackzero=0.0;
	raw.blackone=0.0;
	raw.blacktwo=0.0;
	raw.blackthree=0.0;
	raw.twogreen=true;
    
    ppVersion = PPVERSION;
}




int ProcParams::saveIntoXMP(Exiv2::XmpData &xmpData, const std::string& baseKey ) const
{
	std::string prefix;
	xmpData[baseKey+"rt:"+kXmpVersion] =                int(PPVERSION);

	prefix=baseKey+"rt:ExposureRGB/";
	xmpData[prefix+"rt:Auto"] =                    toneCurve.autoexp;
	xmpData[prefix+"rt:Clip"] =                    toneCurve.clip;
	xmpData[prefix+"rt:Compensation"] =            toneCurve.expcomp;
	xmpData[prefix+"rt:Brightness"] =              toneCurve.brightness;
	xmpData[prefix+"rt:Contrast"] =                toneCurve.contrast;
	xmpData[prefix+"rt:Saturation"] =              toneCurve.saturation;
	xmpData[prefix+"rt:Black"] =                   toneCurve.black;
	xmpData[prefix+"rt:HighlightCompression"] =    toneCurve.hlcompr;
	xmpData[prefix+"rt:HighlightComprThreshold"] = toneCurve.hlcomprthresh;
	xmpData[prefix+"rt:ShadowCompression"] =       toneCurve.shcompr;
	xmpData[prefix+"rt:ToneCurve"]= serializeVector(toneCurve.curve);

	prefix=baseKey+"rt:ChannelMixer/";
	xmpData[prefix+"rt:Red"]   = serializeArray(chmixer.red,3);
	xmpData[prefix+"rt:Green"] = serializeArray(chmixer.green,3);
	xmpData[prefix+"rt:Blue"]  = serializeArray(chmixer.blue,3);

	prefix=baseKey+"rt:ExposureLab/";
	xmpData[prefix+"rt:Brightness"] =         labCurve.brightness;
	xmpData[prefix+"rt:Contrast"] =           labCurve.contrast;
	xmpData[prefix+"rt:Saturation"] =         labCurve.saturation;
	xmpData[prefix+"rt:AvoidColorClipping"] = labCurve.avoidclip;
	xmpData[prefix+"rt:SaturationLimitEnabled"] =  labCurve.enable_saturationlimiter;
	xmpData[prefix+"rt:SaturationLimit"] =    labCurve.saturationlimit;
	xmpData[prefix+"rt:LCurve"] = serializeVector(labCurve.lcurve);
	xmpData[prefix+"rt:aCurve"] = serializeVector(labCurve.acurve);
	xmpData[prefix+"rt:bCurve"] = serializeVector(labCurve.bcurve);

	prefix=baseKey+"rt:Vibrance/";
   	xmpData[prefix+"rt:Enabled"]=         vibrance.enabled;
    xmpData[prefix+"rt:Pastels"]=         vibrance.pastels;
    xmpData[prefix+"rt:Saturated"]=       vibrance.saturated;
    xmpData[prefix+"rt:PSThreshold"]=     vibrance.psthreshold;
    xmpData[prefix+"rt:ProtectSkins"]=    vibrance.protectskins;
    xmpData[prefix+"rt:AvoidColorShift"]= vibrance.avoidcolorshift;
    xmpData[prefix+"rt:PastSatTog"]=      vibrance.pastsattog;

	prefix=baseKey+"rt:Sharpening/";
	xmpData[prefix+"rt:Enabled"] =             sharpening.enabled;
	xmpData[prefix+"rt:Method"] =              sharpening.method;
	xmpData[prefix+"rt:Radius"]=               sharpening.radius;
	xmpData[prefix+"rt:Amount"]=               sharpening.amount;
	xmpData[prefix+"rt:Threshold"]=            sharpening.threshold;
	xmpData[prefix+"rt:OnlyEdges"]=            sharpening.edgesonly;
	xmpData[prefix+"rt:EdgeDetectionRadius"]=  sharpening.edges_radius;
	xmpData[prefix+"rt:EdgeTolerance"]=        sharpening.edges_tolerance;
	xmpData[prefix+"rt:HaloControlEnabled"]=   sharpening.halocontrol;
	xmpData[prefix+"rt:HaloControlAmount"]=    sharpening.halocontrol_amount;
	xmpData[prefix+"rt:DeconvRadius"]=         sharpening.deconvradius;
	xmpData[prefix+"rt:DeconvAmount"]=         sharpening.deconvamount;
	xmpData[prefix+"rt:DeconvDamping"]=        sharpening.deconvdamping;
	xmpData[prefix+"rt:DeconvIterations"]=     sharpening.deconviter;

	prefix=baseKey+"rt:SharpenEdge/";
	xmpData[prefix+"rt:Enabled"]=       sharpenEdge.enabled;
	xmpData[prefix+"rt:Passes"]=        sharpenEdge.passes;
	xmpData[prefix+"rt:ThreeChannels"]= sharpenEdge.threechannels;
	xmpData[prefix+"rt:Amount"]=        sharpenEdge.amount;

	prefix=baseKey+"rt:MicroContrast/";
	xmpData[prefix+"rt:Enabled"]=       sharpenMicro.enabled;
	xmpData[prefix+"rt:Uniformity"]=    sharpenMicro.uniformity;
	xmpData[prefix+"rt:Matrix"]=        sharpenMicro.matrix;
	xmpData[prefix+"rt:Amount"]=        sharpenMicro.amount;

	prefix=baseKey+"rt:WhiteBalance/";
    if (wb.method=="Camera")
    	xmpData[prefix+"rt:Mode"]= "Camera";
    else
    	xmpData[prefix+"rt:Mode"]= "Custom";
    xmpData[prefix+"rt:Temperature"]= wb.temperature;
    xmpData[prefix+"rt:Green"]=       wb.green;

    prefix=baseKey+"rt:ImpulseDenoise/";
    xmpData[prefix+"rt:Enabled"]=   impulseDenoise.enabled;
    xmpData[prefix+"rt:Threshold"]= impulseDenoise.thresh;

    prefix=baseKey+"rt:Defringe/";
    xmpData[prefix+"rt:Enabled"]=   defringe.enabled;
    xmpData[prefix+"rt:Radius"]=    defringe.radius;
    xmpData[prefix+"rt:Threshold"]=	defringe.threshold;

    prefix=baseKey+"rt:PyramidDenoise/";
    xmpData[prefix+"rt:Enabled"]= dirpyrDenoise.enabled;
    xmpData[prefix+"rt:Luma"]=    dirpyrDenoise.luma;
    xmpData[prefix+"rt:Chroma"]=  dirpyrDenoise.chroma;
    xmpData[prefix+"rt:Gamma"]=   dirpyrDenoise.gamma;

    prefix=baseKey+"rt:ShadowHighlights/";
    xmpData[prefix+"rt:Enabled"]=               sh.enabled;
    xmpData[prefix+"rt:HighQuality"]=           sh.hq;
    xmpData[prefix+"rt:Highlights"]=            sh.highlights;
    xmpData[prefix+"rt:HighlightTonalWidth"]=   sh.htonalwidth;
    xmpData[prefix+"rt:Shadows"]=               sh.shadows;
    xmpData[prefix+"rt:ShadowTonalWidth"]=      sh.stonalwidth;
    xmpData[prefix+"rt:LocalContrast"]=         sh.localcontrast;
    xmpData[prefix+"rt:Radius"]=                sh.radius;

    prefix=baseKey+"rt:Crop/";
    xmpData[prefix+"rt:Enabled"]=     crop.enabled;
    xmpData[prefix+"rt:X"]=           crop.x;
    xmpData[prefix+"rt:Y"]=           crop.y;
    xmpData[prefix+"rt:Width"]=       crop.w;
    xmpData[prefix+"rt:Height"]=      crop.h;
 //   xmpData[prefix+"rt:FixedRatio"]=  crop.fixratio;
 //   xmpData[prefix+"rt:Ratio"]=       crop.ratio;
 //   xmpData[prefix+"rt:Orientation"]= crop.orientation;
 //   xmpData[prefix+"rt:Guide"]=       crop.guide;

    prefix=baseKey+"rt:CoarseGeo/";
    xmpData[prefix+"rt:RotationDegree"]=  coarse.rotate;
    xmpData[prefix+"rt:HorizontalFlip"]=  coarse.hflip;
    xmpData[prefix+"rt:VerticalFlip"]=    coarse.vflip;

    prefix=baseKey+"rt:Geometry/";
    //xmpData[prefix+"rt:Enabled"]=    geo.enabled;
    xmpData[prefix+"rt:AutoFill"]= commonTrans.autofill;
    xmpData[prefix+"rt:RotationDegree"]= rotate.degree;
    xmpData[prefix+"rt:DistortionAmount"]= distortion.amount;
    //keyFile.set_boolean ("Distortion", "UseLensFun", distortion.uselensfun);
    xmpData[prefix+"rt:HorizontalPerspective"]= perspective.horizontal;
    xmpData[prefix+"rt:VerticalPerspective"]=   perspective.vertical;

    prefix=baseKey+"rt:CACorrection/";
    //xmpData[prefix+"rt:Enabled"]=    cacorrection.enabled;
    xmpData[prefix+"rt:Red"]=  cacorrection.red;
    xmpData[prefix+"rt:Blue"]= cacorrection.blue;

    prefix=baseKey+"rt:Vignetting/";
    //xmpData[prefix+"rt:Enabled"]= vignetting.enabled;
    xmpData[prefix+"rt:Amount"] = vignetting.amount;
    xmpData[prefix+"rt:Radius"] = vignetting.radius;
    xmpData[prefix+"rt:Strength"]= vignetting.strength;
    xmpData[prefix+"rt:CenterX"] = vignetting.centerX;
    xmpData[prefix+"rt:CenterY"] = vignetting.centerY;

    prefix=baseKey+"rt:HLRecovery/";
    xmpData[prefix+"rt:Enabled"]=  hlrecovery.enabled;
    xmpData[prefix+"rt:Method"]=   hlrecovery.method;

    prefix=baseKey+"rt:Resize/";
    xmpData[prefix+"rt:Enabled"]=   resize.enabled;
    xmpData[prefix+"rt:Scale"]  =   resize.scale;
    xmpData[prefix+"rt:AppliesTo"]= resize.appliesTo;
    xmpData[prefix+"rt:Method"]=    resize.method;
    xmpData[prefix+"rt:DataSpecified"]=  resize.dataspec;
    xmpData[prefix+"rt:Width"] =    resize.width;
    xmpData[prefix+"rt:Height"] =   resize.height;

    prefix=baseKey+"rt:ColorManagement/";
    xmpData[prefix+"rt:InputProfile"] =   icm.input;
    xmpData[prefix+"rt:WorkingProfile"] = icm.working;
    xmpData[prefix+"rt:OutputProfile"] =  icm.output;
    xmpData[prefix+"rt:FreeGamma"] =  icm.gamma;
    xmpData[prefix+"rt:FreeGammaEnabled"] =  icm.freegamma;
    xmpData[prefix+"rt:GammaValue"]=  icm.gampos;
    xmpData[prefix+"rt:GammaSlope"]=  icm.slpos;

    prefix=baseKey+"rt:DirectionalPyramidEqualizer/";
    xmpData[prefix+"rt:Enabled"]=  dirpyrequalizer.enabled ;
    xmpData[prefix+"rt:Coeff"]= serializeArray( dirpyrequalizer.mult,5);

    prefix=baseKey+"rt:HSVEqualizer/";
    //xmpData[prefix+"rt:Enabled"]= hsvequalizer.enabled;
    xmpData[prefix+"rt:HCurve"] = serializeVector(hsvequalizer.hcurve);
    xmpData[prefix+"rt:SCurve"] = serializeVector(hsvequalizer.scurve);
    xmpData[prefix+"rt:VCurve"] = serializeVector(hsvequalizer.vcurve);

    prefix=baseKey+"rt:RawArithmetic/";
    //xmpData[prefix+"rt:Enabled"]=
    xmpData[prefix+"rt:DarkFrameFile"]= raw.dark_frame;
    xmpData[prefix+"rt:DarkFrameAutoSelect"]= raw.df_autoselect ;
    xmpData[prefix+"rt:FlatFieldFile"]= raw.ff_file ;
    xmpData[prefix+"rt:FlatFieldAutoSelect"]= raw.ff_AutoSelect ;
    xmpData[prefix+"rt:FlatFieldBlurRadius"]= raw.ff_BlurRadius ;
    xmpData[prefix+"rt:FlatFieldBlurType"]= raw.ff_BlurType ;

    prefix=baseKey+"rt:RawCACorrection/";
    //xmpData[prefix+"rt:Enabled"]=
    xmpData[prefix+"rt:Auto"]= raw.ca_autocorrect ;
    xmpData[prefix+"rt:Red"] = raw.cared;
    xmpData[prefix+"rt:Blue"]= raw.cablue;

    prefix=baseKey+"rt:HotDeadPixelCorrection/";
    xmpData[prefix+"rt:Enabled"]= raw.hotdeadpix_filt;
    xmpData[prefix+"rt:Threshold"]= raw.hotdeadpix_thresh;

    prefix=baseKey+"rt:RawDenoise/";
    //xmpData[prefix+"rt:Enabled"]=
    xmpData[prefix+"rt:LineDenoise"]= raw.linenoise;

    prefix=baseKey+"rt:Demosaicing/";
    xmpData[prefix+"rt:GreenEqThreshold"] = raw.greenthresh;
    xmpData[prefix+"rt:CcSteps"] =  raw.ccSteps;
    xmpData[prefix+"rt:Method"]  = raw.dmethod;
    xmpData[prefix+"rt:DCBIterations"] =  raw.dcb_iterations;
    xmpData[prefix+"rt:DCBEnhance"] =  raw.dcb_enhance;
    xmpData[prefix+"rt:Enhance"]= raw.all_enhance;

    prefix=baseKey+"rt:RawExposure/";
    xmpData[prefix+"rt:Exposure"] =  raw.expos;
    xmpData[prefix+"rt:HLPreserving"] = raw.preser;
    xmpData[prefix+"rt:Blackzero"] = raw.blackzero;
    xmpData[prefix+"rt:Blackone"] = raw.blackone;
    xmpData[prefix+"rt:Blacktwo"] = raw.blacktwo;
    xmpData[prefix+"rt:Blackthree"] = raw.blackthree;
    xmpData[prefix+"rt:TwoGreen"] = raw.twogreen;
}

int ProcParams::loadFromXMP(Exiv2::XmpData &xmpData, const std::string& baseKey )
{
try{
	std::string prefix;
	if(! readVarFromXmp( xmpData, baseKey+"rt:"+kXmpVersion, ppVersion) )
		return 2;

	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:ExposureRGB")) != xmpData.end()){
		prefix = baseKey+"rt:ExposureRGB/rt:";
		readVarFromXmp( xmpData, prefix+"Auto", toneCurve.autoexp );
		readVarFromXmp( xmpData, prefix+"Clip", toneCurve.clip );
		readVarFromXmp( xmpData, prefix+"Compensation", toneCurve.expcomp );
		readVarFromXmp( xmpData, prefix+"Brightness", toneCurve.brightness );
		readVarFromXmp( xmpData, prefix+"Contrast", toneCurve.contrast );
		readVarFromXmp( xmpData, prefix+"Saturation", toneCurve.saturation );
		readVarFromXmp( xmpData, prefix+"Black", toneCurve.black );
		readVarFromXmp( xmpData, prefix+"HighlightCompression", toneCurve.hlcompr );
		readVarFromXmp( xmpData, prefix+"HighlightComprThreshold", toneCurve.hlcomprthresh );
		readVarFromXmp( xmpData, prefix+"ShadowCompression", toneCurve.shcompr );
		readVarFromXmp( xmpData, prefix+"ToneCurve",toneCurve.curve);
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:ChannelMixer")) != xmpData.end()){
		prefix=baseKey+"rt:ChannelMixer/rt:";
		readVarFromXmp( xmpData, prefix+"Red",chmixer.red,3 );
		readVarFromXmp( xmpData, prefix+"Green",chmixer.green,3 );
		readVarFromXmp( xmpData, prefix+"Blue",chmixer.blue,3 );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:ExposureLab")) != xmpData.end()){
		prefix=baseKey+"rt:ExposureLab/rt:";
		readVarFromXmp( xmpData, prefix+"Brightness", labCurve.brightness );
		readVarFromXmp( xmpData, prefix+"Contrast", labCurve.contrast );
		readVarFromXmp( xmpData, prefix+"Saturation", labCurve.saturation );
		readVarFromXmp( xmpData, prefix+"AvoidColorClipping", labCurve.avoidclip );
		readVarFromXmp( xmpData, prefix+"SaturationLimitEnabled", labCurve.enable_saturationlimiter );
		readVarFromXmp( xmpData, prefix+"SaturationLimit", labCurve.saturationlimit );
		readVarFromXmp( xmpData, prefix+"LCurve",labCurve.lcurve);
		readVarFromXmp( xmpData, prefix+"aCurve",labCurve.acurve);
		readVarFromXmp( xmpData, prefix+"bCurve",labCurve.bcurve);
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:Sharpening")) != xmpData.end()){
		prefix=baseKey+"rt:Sharpening/rt:";
		readVarFromXmp( xmpData, prefix+kXmpEnabled, sharpening.enabled );
		readVarFromXmp( xmpData, prefix+"Method", sharpening.method );
		readVarFromXmp( xmpData, prefix+"Radius", sharpening.radius );
		readVarFromXmp( xmpData, prefix+"Amount", sharpening.amount );
		readVarFromXmp( xmpData, prefix+"Threshold", sharpening.threshold );
		readVarFromXmp( xmpData, prefix+"OnlyEdges", sharpening.edgesonly );
		readVarFromXmp( xmpData, prefix+"EdgeDetectionRadius", sharpening.edges_radius );
		readVarFromXmp( xmpData, prefix+"EdgeTolerance", sharpening.edges_tolerance );
		readVarFromXmp( xmpData, prefix+"HaloControlEnabled", sharpening.halocontrol );
		readVarFromXmp( xmpData, prefix+"HaloControlAmount", sharpening.halocontrol_amount );
		readVarFromXmp( xmpData, prefix+"DeconvRadius", sharpening.deconvradius );
		readVarFromXmp( xmpData, prefix+"DeconvAmount", sharpening.deconvamount );
		readVarFromXmp( xmpData, prefix+"DeconvDamping", sharpening.deconvdamping );
		readVarFromXmp( xmpData, prefix+"DeconvIterations", sharpening.deconviter );
	}
	
    if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:Vibrance")) != xmpData.end()){
		prefix=baseKey+"rt:Vibrance/rt:";
   	    readVarFromXmp( xmpData, prefix+kXmpEnabled, vibrance.enabled);
   	    readVarFromXmp( xmpData, prefix+"Pastels",  vibrance.pastels);
   	    readVarFromXmp( xmpData, prefix+"Saturated", vibrance.saturated);
   	    readVarFromXmp( xmpData, prefix+"PSThreshold", vibrance.psthreshold);
   	    readVarFromXmp( xmpData, prefix+"ProtectSkins", vibrance.protectskins);
   	    readVarFromXmp( xmpData, prefix+"AvoidColorShift", vibrance.avoidcolorshift);
   	    readVarFromXmp( xmpData, prefix+"PastSatTog", vibrance.pastsattog);
    }
	
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:SharpenEdge")) != xmpData.end()){
		prefix=baseKey+"rt:SharpenEdge/rt:";
		readVarFromXmp( xmpData, prefix+kXmpEnabled, sharpenEdge.enabled );
		readVarFromXmp( xmpData, prefix+"Passes", sharpenEdge.passes );
		readVarFromXmp( xmpData, prefix+"Amount", sharpenEdge.amount );
		readVarFromXmp( xmpData, prefix+"ThreeChannels", sharpenEdge.threechannels );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:MicroContrast")) != xmpData.end()){
		prefix=baseKey+"rt:MicroContrast/rt:";
		readVarFromXmp( xmpData, prefix+kXmpEnabled, sharpenMicro.enabled );
		readVarFromXmp( xmpData, prefix+"Amount", sharpenMicro.amount);
		readVarFromXmp( xmpData, prefix+"Uniformity", sharpenMicro.uniformity );
		readVarFromXmp( xmpData, prefix+"Matrix", sharpenMicro.matrix );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:WhiteBalance")) != xmpData.end()){
		prefix=baseKey+"rt:WhiteBalance/rt:";
		readVarFromXmp( xmpData, prefix+"Mode", wb.method );
		readVarFromXmp( xmpData, prefix+"Temperature", wb.temperature );
		readVarFromXmp( xmpData, prefix+"Green", wb.green );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:ImpulseDenoise")) != xmpData.end()){
		prefix=baseKey+"rt:ImpulseDenoise/rt:";
		readVarFromXmp( xmpData, prefix+kXmpEnabled, impulseDenoise.enabled );
		readVarFromXmp( xmpData, prefix+"Threshold", impulseDenoise.thresh );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:Defringe")) != xmpData.end()){
		prefix=baseKey+"rt:Defringe/rt:";
		readVarFromXmp( xmpData, prefix+kXmpEnabled, defringe.enabled );
		readVarFromXmp( xmpData, prefix+"Radius", defringe.radius );
		readVarFromXmp( xmpData, prefix+"Threshold", defringe.threshold );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:PyramidDenoise")) != xmpData.end()){
		prefix=baseKey+"rt:PyramidDenoise/rt:";
		readVarFromXmp( xmpData, prefix+kXmpEnabled, dirpyrDenoise.enabled );
		readVarFromXmp( xmpData, prefix+"Luma", dirpyrDenoise.luma );
		readVarFromXmp( xmpData, prefix+"Chroma", dirpyrDenoise.chroma );
		readVarFromXmp( xmpData, prefix+"Gamma", dirpyrDenoise.gamma );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:ShadowHighlights")) != xmpData.end()){
		prefix=baseKey+"rt:ShadowHighlights/rt:";
		readVarFromXmp( xmpData, prefix+kXmpEnabled, sh.enabled );
		readVarFromXmp( xmpData, prefix+"HighQuality", sh.hq );
		readVarFromXmp( xmpData, prefix+"Highlights", sh.highlights );
		readVarFromXmp( xmpData, prefix+"HighlightTonalWidth", sh.htonalwidth );
		readVarFromXmp( xmpData, prefix+"Shadows", sh.shadows );
		readVarFromXmp( xmpData, prefix+"ShadowTonalWidth", sh.stonalwidth );
		readVarFromXmp( xmpData, prefix+"LocalContrast", sh.localcontrast );
		readVarFromXmp( xmpData, prefix+"Radius", sh.radius );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:Crop")) != xmpData.end()){
		prefix=baseKey+"rt:Crop/rt:";
		readVarFromXmp( xmpData, prefix+kXmpEnabled, crop.enabled );
		readVarFromXmp( xmpData, prefix+"X", crop.x );
		readVarFromXmp( xmpData, prefix+"Y", crop.y );
		readVarFromXmp( xmpData, prefix+"Width", crop.w );
		readVarFromXmp( xmpData, prefix+"Height", crop.h );
	 //   xmpData[prefix+"rt:FixedRatio"]=  crop.fixratio;
	 //   xmpData[prefix+"rt:Ratio"]=       crop.ratio;
	 //   xmpData[prefix+"rt:Orientation"]= crop.orientation;
	 //   xmpData[prefix+"rt:Guide"]=       crop.guide;
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:CoarseGeo")) != xmpData.end()){
		prefix=baseKey+"rt:CoarseGeo/rt:";
		readVarFromXmp( xmpData, prefix+"RotationDegree", coarse.rotate );
		readVarFromXmp( xmpData, prefix+"HorizontalFlip", coarse.hflip );
		readVarFromXmp( xmpData, prefix+"VerticalFlip", coarse.vflip );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:Geometry")) != xmpData.end()){
		prefix=baseKey+"rt:Geometry/rt:";
		//readVarFromXmp( xmpData, prefix+kXmpEnabled, geo.enabled);
		readVarFromXmp( xmpData, prefix+"AutoFill", commonTrans.autofill );
		readVarFromXmp( xmpData, prefix+"RotationDegree", rotate.degree );
		readVarFromXmp( xmpData, prefix+"DistortionAmount", distortion.amount );
		readVarFromXmp( xmpData, prefix+"HorizontalPerspective", perspective.horizontal );
		readVarFromXmp( xmpData, prefix+"VerticalPerspective", perspective.vertical );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:CACorrection")) != xmpData.end()){
		prefix=baseKey+"rt:CACorrection/rt:";
		//readVarFromXmp( xmpData, prefix+kXmpEnabled, cacorrection.enabled);
		readVarFromXmp( xmpData, prefix+"Red", cacorrection.red );
		readVarFromXmp( xmpData, prefix+"Blue", cacorrection.blue );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:Vignetting")) != xmpData.end()){
		prefix=baseKey+"rt:Vignetting/rt:";
		//readVarFromXmp( xmpData, prefix+kXmpEnabled, vignetting.enabled);
		readVarFromXmp( xmpData, prefix+"Amount", vignetting.amount );
		readVarFromXmp( xmpData, prefix+"Radius", vignetting.radius );
		readVarFromXmp( xmpData, prefix+"Strength", vignetting.strength );
		readVarFromXmp( xmpData, prefix+"CenterX", vignetting.centerX );
		readVarFromXmp( xmpData, prefix+"CenterY", vignetting.centerY );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:HLRecovery")) != xmpData.end()){
		prefix=baseKey+"rt:HLRecovery/rt:";
		readVarFromXmp( xmpData, prefix+kXmpEnabled, hlrecovery.enabled );
		readVarFromXmp( xmpData, prefix+"Method", hlrecovery.method );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:Resize")) != xmpData.end()){
		prefix=baseKey+"rt:Resize/rt:";
		readVarFromXmp( xmpData, prefix+kXmpEnabled, resize.enabled );
		readVarFromXmp( xmpData, prefix+"Scale", resize.scale );
		readVarFromXmp( xmpData, prefix+"AppliesTo", resize.appliesTo );
		readVarFromXmp( xmpData, prefix+"Method", resize.method );
		readVarFromXmp( xmpData, prefix+"DataSpecified", resize.dataspec );
		readVarFromXmp( xmpData, prefix+"Width", resize.width );
		readVarFromXmp( xmpData, prefix+"Height", resize.height );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:ColorManagement")) != xmpData.end()){
		prefix=baseKey+"rt:ColorManagement/rt:";
		readVarFromXmp( xmpData, prefix+"InputProfile", icm.input );
		readVarFromXmp( xmpData, prefix+"WorkingProfile", icm.working );
		readVarFromXmp( xmpData, prefix+"OutputProfile", icm.output );
		readVarFromXmp( xmpData, prefix+"FreeGamma", icm.gamma );
		readVarFromXmp( xmpData, prefix+"FreeGammaEnabled", icm.freegamma );
		readVarFromXmp( xmpData, prefix+"GammaValue", icm.gampos );
		readVarFromXmp( xmpData, prefix+"GammaSlope", icm.slpos );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:DirectionalPyramidEqualizer")) != xmpData.end()){
		prefix=baseKey+"rt:DirectionalPyramidEqualizer/rt:";
		readVarFromXmp( xmpData, prefix+kXmpEnabled, dirpyrequalizer.enabled );
		readVarFromXmp( xmpData, prefix+"Coeff",dirpyrequalizer.mult,5 );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:DirectionalPyramidEqualizer")) != xmpData.end()){
		prefix=baseKey+"rt:HSVEqualizer/rt:";
		//readVarFromXmp( xmpData, prefix+kXmpEnabled, hsvequalizer.enabled);
		readVarFromXmp( xmpData,prefix+"HCurve",hsvequalizer.hcurve);
		readVarFromXmp( xmpData,prefix+"SCurve",hsvequalizer.scurve);
		readVarFromXmp( xmpData,prefix+"VCurve",hsvequalizer.vcurve);
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:RawArithmetic")) != xmpData.end()){
		prefix=baseKey+"rt:RawArithmetic/rt:";
		//readVarFromXmp( xmpData, prefix+kXmpEnabled, raw.arithmeticEnabled);
		readVarFromXmp( xmpData, prefix+"DarkFrameFile", raw.dark_frame );
		readVarFromXmp( xmpData, prefix+"DarkFrameAutoSelect", raw.df_autoselect );
		readVarFromXmp( xmpData, prefix+"FlatFieldFile", raw.ff_file );
		readVarFromXmp( xmpData, prefix+"FlatFieldAutoSelect", raw.ff_AutoSelect );
		readVarFromXmp( xmpData, prefix+"FlatFieldBlurRadius", raw.ff_BlurRadius );
		readVarFromXmp( xmpData, prefix+"FlatFieldBlurType", raw.ff_BlurType );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:RawCACorrection")) != xmpData.end()){
		prefix=baseKey+"rt:RawCACorrection/rt:";
		//readVarFromXmp( xmpData, prefix+kXmpEnabled, raw.ca_enabled);
		readVarFromXmp( xmpData, prefix+"Auto",	raw.ca_autocorrect );
		readVarFromXmp( xmpData, prefix+"Red", raw.cared );
		readVarFromXmp( xmpData, prefix+"Blue", raw.cablue );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:HotDeadPixelCorrection")) != xmpData.end()){
		prefix=baseKey+"rt:HotDeadPixelCorrection/rt:";
		readVarFromXmp( xmpData, prefix+kXmpEnabled, raw.hotdeadpix_filt );
		readVarFromXmp( xmpData, prefix+"Threshold", raw.hotdeadpix_thresh );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:RawDenoise")) != xmpData.end()){
		prefix=baseKey+"rt:RawDenoise/rt:";
		//readVarFromXmp( xmpData, prefix+kXmpEnabled, raw.linenoiseEnabled );
		readVarFromXmp( xmpData, prefix+"LineDenoise", raw.linenoise );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:Demosaicing")) != xmpData.end()){
		prefix=baseKey+"rt:Demosaicing/rt:";
		readVarFromXmp( xmpData, prefix+"GreenEqThreshold", raw.greenthresh );
		readVarFromXmp( xmpData, prefix+"CcSteps", raw.ccSteps );
		readVarFromXmp( xmpData, prefix+"Method", raw.dmethod );
		readVarFromXmp( xmpData, prefix+"DCBIterations", raw.dcb_iterations );
		readVarFromXmp( xmpData, prefix+"DCBEnhance", raw.dcb_enhance );
		readVarFromXmp( xmpData, prefix+"Enhance", raw.all_enhance );
	}
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:RawExposure")) != xmpData.end()){
		prefix=baseKey+"rt:RawExposure/rt:";
		readVarFromXmp( xmpData, prefix+"Exposure", raw.expos );
		readVarFromXmp( xmpData, prefix+"HLPreserving", raw.preser );
		readVarFromXmp( xmpData, prefix+"Black0", raw.blackzero );
		readVarFromXmp( xmpData, prefix+"Black1", raw.blackone );
		readVarFromXmp( xmpData, prefix+"Black2", raw.blacktwo );
		readVarFromXmp( xmpData, prefix+"Black3", raw.blackthree );
		readVarFromXmp( xmpData, prefix+"TwoGreen", raw.twogreen );
	}
	return 0;
}catch( Exiv2::AnyError &e){
	printf("loadFromXMP error: %s\n",e.what());
	return 2;
}
}

int ProcParams::saveParams ( Glib::ustring fname ) const
{
	Exiv2::XmpData  xmpData;

    xmpData[Glib::ustring::compose("Xmp.rt.%1",kXmpVersion)] = VERSION;

    Glib::ustring baseKey(Glib::ustring::compose("Xmp.rt.%1",kXmpProcessing));

    Exiv2::XmpTextValue tv("");
    tv.setXmpArrayType(Exiv2::XmpValue::xaBag);
    xmpData.add(Exiv2::XmpKey(baseKey), &tv);
    saveIntoXMP( xmpData, Glib::ustring::compose("Xmp.rt.%1[1]/",kXmpProcessing) );

    std::string xmpPacket;
    if (0 != Exiv2::XmpParser::encode(xmpPacket, xmpData)) {
        return 1;
    }

    FILE *f = safe_g_fopen (fname, "wt");
    if (f==NULL)
         return 1;
    else {
    	 fwrite( xmpPacket.c_str(), 1, xmpPacket.size(), f );
         fclose (f);
         return 0;
    }
}

int ProcParams::loadParams( Glib::ustring fname )
{
	Exiv2::XmpData  xmpData;
    FILE *f = safe_g_fopen(fname,"rb");
	if ( f == NULL )
		return 1;
	fseek (f, 0, SEEK_END);
	long filesize = ftell (f);
	char *buffer=new char[filesize];

	fseek (f, 0, SEEK_SET);
    fread( buffer, 1, filesize, f);
    fclose(f);
    std::string xmpPacket(buffer,buffer+filesize);

    if (0 != Exiv2::XmpParser::decode(xmpData,xmpPacket) ) {
        return 1;
    }
    std::string key(Glib::ustring::compose("Xmp.rt.%1[1]/",kXmpProcessing));
    loadFromXMP( xmpData, key );
    return 0;
}

int ProcParams::write (Glib::ustring &fname, Glib::ustring &content) const {

    int error = 0;
    if (fname.length()) {
    	FILE *f;
        f = safe_g_fopen (fname, "wt");

        if (f==NULL)
            error = 1;
        else {
            fprintf (f, "%s", content.c_str());
            fclose (f);
        }
    }
    return error;
}

int ProcParams::load (Glib::ustring fname, int *rank) {

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

ppVersion = PPVERSION;
appVersion = APPVERSION;
if (keyFile.has_group ("Version")) {    
    if (keyFile.has_key ("Version", "AppVersion")) appVersion = keyFile.get_string  ("Version", "AppVersion");
    if (keyFile.has_key ("Version", "Version"))    ppVersion  = keyFile.get_integer ("Version", "Version");
}

if (keyFile.has_group ("General")) {
	if( rank!= NULL){
		*rank = 0;
		if (keyFile.has_key ("General", "InTrash")){
			if( keyFile.get_boolean ("General", "InTrash") )
				*rank = -1;
		}
		if( *rank != -1  && keyFile.has_key ("General", "Rank"))
			*rank = keyFile.get_integer ("General", "Rank");
	}
    //if (keyFile.has_key ("General", "ColorLabel"))  colorlabel  = keyFile.get_integer ("General", "ColorLabel");
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
    if (ppVersion>200)
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
	if (keyFile.has_key ("Luminance Curve", "AvoidColorClipping"))  labCurve.avoidclip               = keyFile.get_boolean ("Luminance Curve", "AvoidColorClipping");
    if (keyFile.has_key ("Luminance Curve", "SaturationLimiter"))   labCurve.enable_saturationlimiter= keyFile.get_boolean ("Luminance Curve", "SaturationLimiter");
    if (keyFile.has_key ("Luminance Curve", "SaturationLimit"))     labCurve.saturationlimit         = keyFile.get_double  ("Luminance Curve", "SaturationLimit");	
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

    // load edge sharpening
if (keyFile.has_group ("SharpenEdge")) {
    if (keyFile.has_key ("SharpenEdge", "Enabled"))             sharpenEdge.enabled         = keyFile.get_boolean ("SharpenEdge", "Enabled");
    if (keyFile.has_key ("SharpenEdge", "Passes"))              sharpenEdge.passes          = keyFile.get_integer  ("SharpenEdge", "Passes");
    if (keyFile.has_key ("SharpenEdge", "Strength"))            sharpenEdge.amount          = keyFile.get_double  ("SharpenEdge", "Strength");
    if (keyFile.has_key ("SharpenEdge", "ThreeChannels"))       sharpenEdge.threechannels   = keyFile.get_boolean ("SharpenEdge", "ThreeChannels");
}

    // load micro-contrast sharpening
if (keyFile.has_group ("SharpenMicro")) {
    if (keyFile.has_key ("SharpenMicro", "Enabled"))            sharpenMicro.enabled        = keyFile.get_boolean ("SharpenMicro", "Enabled");
    if (keyFile.has_key ("SharpenMicro", "Matrix"))             sharpenMicro.matrix         = keyFile.get_boolean ("SharpenMicro", "Matrix");
    if (keyFile.has_key ("SharpenMicro", "Strength"))           sharpenMicro.amount         = keyFile.get_double  ("SharpenMicro", "Strength");
    if (keyFile.has_key ("SharpenMicro", "Uniformity"))         sharpenMicro.uniformity     = keyFile.get_double  ("SharpenMicro", "Uniformity");
}

    // load vibrance
if (keyFile.has_group ("Vibrance")) {
    if (keyFile.has_key ("Vibrance", "Enabled"))                vibrance.enabled            = keyFile.get_boolean ("Vibrance", "Enabled");
    if (keyFile.has_key ("Vibrance", "Pastels"))                vibrance.pastels            = keyFile.get_integer ("Vibrance", "Pastels");
    if (keyFile.has_key ("Vibrance", "Saturated"))              vibrance.saturated          = keyFile.get_integer ("Vibrance", "Saturated");
    if (keyFile.has_key ("Vibrance", "PSThreshold"))            vibrance.psthreshold        = keyFile.get_integer ("Vibrance", "PSThreshold");
    if (keyFile.has_key ("Vibrance", "ProtectSkins"))           vibrance.protectskins       = keyFile.get_boolean ("Vibrance", "ProtectSkins");
    if (keyFile.has_key ("Vibrance", "AvoidColorShift"))        vibrance.avoidcolorshift    = keyFile.get_boolean ("Vibrance", "AvoidColorShift");
    if (keyFile.has_key ("Vibrance", "PastSatTog"))         	vibrance.pastsattog         = keyFile.get_boolean ("Vibrance", "PastSatTog");
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
	if (keyFile.has_key ("Directional Pyramid Denoising", "LumCurve"))    dirpyrDenoise.lumcurve   = keyFile.get_double_list ("Directional Pyramid Denoising", "LumCurve");
	if (keyFile.has_key ("Directional Pyramid Denoising", "ChromCurve"))  dirpyrDenoise.chromcurve = keyFile.get_double_list ("Directional Pyramid Denoising", "ChromCurve");
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
    if (keyFile.has_key ("Color Management", "BlendCMSMatrix"))   icm.blendCMSMatrix = keyFile.get_boolean ("Color Management", "BlendCMSMatrix");
    if (keyFile.has_key ("Color Management", "WorkingProfile")) icm.working = keyFile.get_string ("Color Management", "WorkingProfile");
    if (keyFile.has_key ("Color Management", "OutputProfile"))  icm.output  = keyFile.get_string ("Color Management", "OutputProfile");
    if (keyFile.has_key ("Color Management", "Gammafree"))  icm.gamma  = keyFile.get_string ("Color Management", "Gammafree");
    if (keyFile.has_key ("Color Management", "Freegamma"))  icm.freegamma  = keyFile.get_boolean ("Color Management", "Freegamma");
    if (keyFile.has_key ("Color Management", "GammaVal"))  icm.gampos  = keyFile.get_double ("Color Management", "GammaVal");
    if (keyFile.has_key ("Color Management", "GammaSlope"))  icm.slpos  = keyFile.get_double ("Color Management", "GammaSlope");
	
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

	// load HSV equalizer parameters
if (keyFile.has_group ("HSV Equalizer")) {
	if (ppVersion>=300) {
		if (keyFile.has_key ("HSV Equalizer", "HCurve"))          hsvequalizer.hcurve      = keyFile.get_double_list ("HSV Equalizer", "HCurve");
		if (keyFile.has_key ("HSV Equalizer", "SCurve"))          hsvequalizer.scurve      = keyFile.get_double_list ("HSV Equalizer", "SCurve");
		if (keyFile.has_key ("HSV Equalizer", "VCurve"))          hsvequalizer.vcurve      = keyFile.get_double_list ("HSV Equalizer", "VCurve");
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
	if (keyFile.has_key ("RAW", "ALLEnhance"))    raw.all_enhance =keyFile.get_boolean("RAW", "ALLEnhance");
	
	if (keyFile.has_key ("RAW", "PreExposure"))   	  raw.expos =keyFile.get_double("RAW", "PreExposure");
	if (keyFile.has_key ("RAW", "PrePreserv"))   	  raw.preser =keyFile.get_double("RAW", "PrePreserv");
	if (keyFile.has_key ("RAW", "PreBlackzero"))   	  raw.blackzero =keyFile.get_double("RAW", "PreBlackzero");
	if (keyFile.has_key ("RAW", "PreBlackone"))   	  raw.blackone =keyFile.get_double("RAW", "PreBlackone");
	if (keyFile.has_key ("RAW", "PreBlacktwo"))   	  raw.blacktwo =keyFile.get_double("RAW", "PreBlacktwo");
	if (keyFile.has_key ("RAW", "PreBlackthree"))   	  raw.blackthree =keyFile.get_double("RAW", "PreBlackthree");
	if (keyFile.has_key ("RAW", "PreTwoGreen"))   	  raw.twogreen =keyFile.get_boolean("RAW", "PreTwoGreen");	
	
}

    // load exif change settings
/*if (keyFile.has_group ("Exif")) {
    std::vector<Glib::ustring> keys = keyFile.get_keys ("Exif");
    exif.resize (keys.size());
    for (int i=0; i<(int)keys.size(); i++) {
        exif[i].field = keys[i];
        exif[i].value = keyFile.get_string ("Exif", keys[i]);
    }
}*/

    // load iptc change settings
/*if (keyFile.has_group ("IPTC")) {
    std::vector<Glib::ustring> keys = keyFile.get_keys ("IPTC");
    iptc.resize (keys.size());
    for (int i=0; i<(int)keys.size(); i++) {
        iptc[i].field = keys[i];
        iptc[i].values = keyFile.get_string_list ("IPTC", keys[i]);
    }
}*/


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

const Glib::ustring ColorManagementParams::NoICMString = Glib::ustring("No ICM: sRGB output");

bool operator==(const DirPyrEqualizerParams & a, const DirPyrEqualizerParams & b) {
	if(a.enabled != b.enabled)
		return false;
		
	for(int i = 0; i < 5; i++) {
		if(a.mult[i] != b.mult[i])
			return false;
	}
	return true;
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
		&& labCurve.avoidclip == other.labCurve.avoidclip
		&& labCurve.enable_saturationlimiter == other.labCurve.enable_saturationlimiter
		&& labCurve.saturationlimit == other.labCurve.saturationlimit			
		&& sharpenEdge.enabled == other.sharpenEdge.enabled
		&& sharpenEdge.passes == other.sharpenEdge.passes
		&& sharpenEdge.amount == other.sharpenEdge.amount
		&& sharpenEdge.threechannels == other.sharpenEdge.threechannels
		&& sharpenMicro.enabled == other.sharpenMicro.enabled
		&& sharpenMicro.matrix == other.sharpenMicro.matrix
		&& sharpenMicro.amount == other.sharpenMicro.amount
		&& sharpenMicro.uniformity == other.sharpenMicro.uniformity
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
		&& vibrance.enabled == other.vibrance.enabled
		&& vibrance.pastels == other.vibrance.pastels
		&& vibrance.saturated == other.vibrance.saturated
		&& vibrance.psthreshold == other.vibrance.psthreshold
		&& vibrance.protectskins == other.vibrance.protectskins
		&& vibrance.avoidcolorshift == other.vibrance.avoidcolorshift
		&& vibrance.pastsattog == other.vibrance.pastsattog
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
		&& dirpyrDenoise.lumcurve == other.dirpyrDenoise.lumcurve
		&& dirpyrDenoise.chromcurve == other.dirpyrDenoise.chromcurve
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
		&& raw.ff_file == other.raw.ff_file
		&& raw.ff_AutoSelect == other.raw.ff_AutoSelect
		&& raw.ff_BlurRadius == other.raw.ff_BlurRadius
		&& raw.ff_BlurType == other.raw.ff_BlurType	
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
		&& icm.blendCMSMatrix == other.icm.blendCMSMatrix
		&& icm.working == other.icm.working
		&& icm.output == other.icm.output
		&& icm.gamma == other.icm.gamma		
		&& icm.freegamma == other.icm.freegamma			
		&& icm.gampos == other.icm.gampos	
		&& icm.slpos == other.icm.slpos			
		&& dirpyrequalizer == other.dirpyrequalizer
		&& hsvequalizer.hcurve == other.hsvequalizer.hcurve
		&& hsvequalizer.scurve == other.hsvequalizer.scurve
		&& hsvequalizer.vcurve == other.hsvequalizer.vcurve
		&& raw.expos==other.raw.expos
		&& raw.preser==other.raw.preser 
		&& raw.preser==other.raw.preser
		&& raw.blackzero==other.raw.blackzero
		&& raw.blackone==other.raw.blackone
		&& raw.blacktwo==other.raw.blacktwo
		&& raw.blackthree==other.raw.blackthree
		&& raw.twogreen==other.raw.twogreen;
	
	}

bool ProcParams::operator!= (const ProcParams& other) {

    return !(*this==other);
}

}
}
