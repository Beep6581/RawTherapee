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
#include <glibmm.h>
#include <sstream>
#include <cstring>
#include "rt_math.h"

#include "safegtk.h"
#include "../rtgui/multilangmgr.h"
#include "procparams.h"
#include "../rtgui/version.h"
#include "../rtgui/ppversion.h"
#include "../rtgui/mydiagonalcurve.h"
#include "../rtgui/myflatcurve.h"
#include "safekeyfile.h"
#include "rawimage.h"
#include "../rtgui/ppversion.h"
#include "../rtgui/paramsedited.h"
#include "dcp.h"

#define APPVERSION VERSION

using namespace std;

namespace rtengine {
namespace procparams {

const char *RAWParams::methodstring[RAWParams::numMethods]={"eahd", "hphd", "vng4", "dcb", "amaze", "ahd", "fast" };
const char *RAWParams::ff_BlurTypestring[RAWParams::numFlatFileBlurTypes]={/*"Parametric",*/ "Area Flatfield", "Vertical Flatfield", "Horizontal Flatfield", "V+H Flatfield"};
std::vector<WBEntry*> WBParams::wbEntries;

void WBParams::init() {
    // Creation of the different methods and its associated temperature value
    wbEntries.push_back(new WBEntry("Camera"              ,WBT_CAMERA,      M("TP_WBALANCE_CAMERA"),        0));
    wbEntries.push_back(new WBEntry("Auto"                ,WBT_AUTO,        M("TP_WBALANCE_AUTO"),          0));
    wbEntries.push_back(new WBEntry("Daylight"            ,WBT_DAYLIGHT,    M("TP_WBALANCE_DAYLIGHT"),   5300));
    wbEntries.push_back(new WBEntry("Cloudy"              ,WBT_CLOUDY,      M("TP_WBALANCE_CLOUDY"),     6200));
    wbEntries.push_back(new WBEntry("Shade"               ,WBT_SHADE,       M("TP_WBALANCE_SHADE"),      7600));
    wbEntries.push_back(new WBEntry("Tungsten"            ,WBT_TUNGSTEN,    M("TP_WBALANCE_TUNGSTEN"),   2856));
    wbEntries.push_back(new WBEntry("Fluo F1"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO1"),      6430));
    wbEntries.push_back(new WBEntry("Fluo F2"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO2"),      4230));
    wbEntries.push_back(new WBEntry("Fluo F3"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO3"),      3450));
    wbEntries.push_back(new WBEntry("Fluo F4"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO4"),      2940));
    wbEntries.push_back(new WBEntry("Fluo F5"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO5"),      6350));
    wbEntries.push_back(new WBEntry("Fluo F6"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO6"),      4150));
    wbEntries.push_back(new WBEntry("Fluo F7"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO7"),      6500));
    wbEntries.push_back(new WBEntry("Fluo F8"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO8"),      5020));
    wbEntries.push_back(new WBEntry("Fluo F9"             ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO9"),      4330));
    wbEntries.push_back(new WBEntry("Fluo F10"            ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO10"),     5300));
    wbEntries.push_back(new WBEntry("Fluo F11"            ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO11"),     4000));
    wbEntries.push_back(new WBEntry("Fluo F12"            ,WBT_FLUORESCENT, M("TP_WBALANCE_FLUO12"),     3000));
    wbEntries.push_back(new WBEntry("HMI Lamp"            ,WBT_LAMP,        M("TP_WBALANCE_HMI"),        4800));
    wbEntries.push_back(new WBEntry("GTI Lamp"            ,WBT_LAMP,        M("TP_WBALANCE_GTI"),        5000));
    wbEntries.push_back(new WBEntry("JudgeIII Lamp"       ,WBT_LAMP,        M("TP_WBALANCE_JUDGEIII"),   5100));
    wbEntries.push_back(new WBEntry("Solux Lamp 3500K"    ,WBT_LAMP,        M("TP_WBALANCE_SOLUX35"),    3480));
    wbEntries.push_back(new WBEntry("Solux Lamp 4100K"    ,WBT_LAMP,        M("TP_WBALANCE_SOLUX41"),    3930));
    wbEntries.push_back(new WBEntry("Solux Lamp 4700K"    ,WBT_LAMP,        M("TP_WBALANCE_SOLUX47"),    4700));
    wbEntries.push_back(new WBEntry("NG Solux Lamp 4700K" ,WBT_LAMP,        M("TP_WBALANCE_SOLUX47_NG"), 4480));
    wbEntries.push_back(new WBEntry("LED LSI Lumelex 2040",WBT_LED,         M("TP_WBALANCE_LED_LSI"),    3000));
    wbEntries.push_back(new WBEntry("LED CRS SP12 WWMR16" ,WBT_LED,         M("TP_WBALANCE_LED_CRS"),    3050));
    wbEntries.push_back(new WBEntry("Flash 5500K"         ,WBT_FLASH,       M("TP_WBALANCE_FLASH55"),    5500));
    wbEntries.push_back(new WBEntry("Flash 6000K"         ,WBT_FLASH,       M("TP_WBALANCE_FLASH60"),    6000));
    wbEntries.push_back(new WBEntry("Flash 6500K"         ,WBT_FLASH,       M("TP_WBALANCE_FLASH65"),    6500));
    // Should remain the last one
    wbEntries.push_back(new WBEntry("Custom"              ,WBT_CUSTOM,      M("TP_WBALANCE_CUSTOM"),        0));
}

void WBParams::cleanup() {
    for (unsigned int i=0; i<wbEntries.size(); i++) {
        delete wbEntries[i];
    }
}

// Maps crop to resized width (e.g. smaller previews)
void CropParams::mapToResized(int resizedWidth, int resizedHeight, int scale, int &x1, int &x2, int &y1, int &y2) const {
    x1 = 0, x2 = resizedWidth, y1 = 0, y2 = resizedHeight;
    if (enabled) {
        x1 = min(resizedWidth-1,  max(0, x / scale));
        y1 = min(resizedHeight-1, max(0, y / scale));
        x2 = min(resizedWidth,    max(0, (x+w) / scale));
        y2 = min(resizedHeight,   max(0, (y+h) / scale));
    }
}

ProcParams::ProcParams () { 

    setDefaults (); 
}       

void ProcParams::init () {

    WBParams::init();
}

void ProcParams::cleanup () {

    WBParams::cleanup();
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
    toneCurve.shcompr       = 50;
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

    rgbCurves.rcurve.clear ();
    rgbCurves.rcurve.push_back(DCT_Linear);
    rgbCurves.gcurve.clear ();
    rgbCurves.gcurve.push_back(DCT_Linear);
    rgbCurves.bcurve.clear ();
    rgbCurves.bcurve.push_back(DCT_Linear);


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

    //colorBoost.amount                   = 0;
    //colorBoost.avoidclip                = false;
    //colorBoost.enable_saturationlimiter = false;
    //colorBoost.saturationlimit          = 50;

    wb.method       = "Camera";
    wb.temperature  = 6504;
    wb.green        = 1.0;

    //colorShift.a    = 0;
    //colorShift.b    = 0;

    //lumaDenoise.enabled         = false;
    //lumaDenoise.radius          = 1.9;
    //lumaDenoise.edgetolerance   = 2000;

    //colorDenoise.enabled        = false;
    //colorDenoise.edgesensitive  = false;
    //colorDenoise.edgetolerance  = 2000;

    impulseDenoise.enabled      = false;
    impulseDenoise.thresh       = 50;

    defringe.enabled            = false;
    defringe.radius             = 2.0;
    defringe.threshold          = 25;

    dirpyrDenoise.enabled       = false;
    dirpyrDenoise.luma          = 10;
    dirpyrDenoise.chroma        = 10;
    dirpyrDenoise.gamma         = 2.0;

    edgePreservingDecompositionUI.enabled = false;
    edgePreservingDecompositionUI.Strength = 0.25;
    edgePreservingDecompositionUI.EdgeStopping = 1.4;
    edgePreservingDecompositionUI.Scale = 1.0;
    edgePreservingDecompositionUI.ReweightingIterates = 0;

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
    
    lensProf.lcpFile="";
    lensProf.useDist=lensProf.useVign=true;
    lensProf.useCA=false;

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
    icm.preferredProfile = (short)rtengine::Daylight;
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




int ProcParams::saveIntoXMP(Exiv2::XmpData &xmpData, const std::string& baseKey, ParamsEdited* pedited ) const
{
	std::string prefix;
	xmpData[baseKey+"rt:"+kXmpVersion] =                int(PPVERSION);

	prefix=baseKey+"rt:ExposureRGB/";
    if (!pedited || pedited->toneCurve.autoexp)    xmpData[prefix+"rt:Auto"] =                    toneCurve.autoexp;
    if (!pedited || pedited->toneCurve.clip)       xmpData[prefix+"rt:Clip"] =                    toneCurve.clip;
    if (!pedited || pedited->toneCurve.expcomp)    xmpData[prefix+"rt:Compensation"] =            toneCurve.expcomp;
    if (!pedited || pedited->toneCurve.brightness) xmpData[prefix+"rt:Brightness"] =              toneCurve.brightness;
    if (!pedited || pedited->toneCurve.contrast)   xmpData[prefix+"rt:Contrast"] =                toneCurve.contrast;
    if (!pedited || pedited->toneCurve.saturation) xmpData[prefix+"rt:Saturation"] =              toneCurve.saturation;
    if (!pedited || pedited->toneCurve.black)      xmpData[prefix+"rt:Black"] =                   toneCurve.black;
    if (!pedited || pedited->toneCurve.hlcompr)    xmpData[prefix+"rt:HighlightCompression"] =    toneCurve.hlcompr;
    if (!pedited || pedited->toneCurve.hlcomprthresh) xmpData[prefix+"rt:HighlightComprThreshold"] = toneCurve.hlcomprthresh;
    if (!pedited || pedited->toneCurve.shcompr)    xmpData[prefix+"rt:ShadowCompression"] =       toneCurve.shcompr;
    if (!pedited || pedited->toneCurve.curve)      xmpData[prefix+"rt:ToneCurve"] =               serializeVector(toneCurve.curve);

	prefix=baseKey+"rt:ChannelMixer/";
	if (!pedited || pedited->chmixer.red[0] || pedited->chmixer.red[1] || pedited->chmixer.red[2])
		xmpData[prefix+"rt:Red"]   = serializeArray(chmixer.red,3);
	if (!pedited || pedited->chmixer.green[0] || pedited->chmixer.green[1] || pedited->chmixer.green[2])
		xmpData[prefix+"rt:Green"] = serializeArray(chmixer.green,3);
	if (!pedited || pedited->chmixer.blue[0] || pedited->chmixer.blue[1] || pedited->chmixer.blue[2])
		xmpData[prefix+"rt:Blue"]  = serializeArray(chmixer.blue,3);

    prefix=baseKey+"rt:ExposureLab/";
    if (!pedited || pedited->labCurve.brightness)      xmpData[prefix+"rt:Brightness"] =         labCurve.brightness;
    if (!pedited || pedited->labCurve.contrast)        xmpData[prefix+"rt:Contrast"] =           labCurve.contrast;
    if (!pedited || pedited->labCurve.saturation)      xmpData[prefix+"rt:Saturation"] =         labCurve.saturation;
    if (!pedited || pedited->labCurve.avoidclip)       xmpData[prefix+"rt:AvoidColorClipping"] = labCurve.avoidclip;
    if (!pedited || pedited->labCurve.enable_saturationlimiter) xmpData[prefix+"rt:SaturationLimitEnabled"] =  labCurve.enable_saturationlimiter;
    if (!pedited || pedited->labCurve.saturationlimit) xmpData[prefix+"rt:SaturationLimit"] =    labCurve.saturationlimit;
    if (!pedited || pedited->labCurve.lcurve)          xmpData[prefix+"rt:LCurve"] =             serializeVector(labCurve.lcurve);
    if (!pedited || pedited->labCurve.acurve)          xmpData[prefix+"rt:aCurve"] =             serializeVector(labCurve.acurve);
    if (!pedited || pedited->labCurve.bcurve)          xmpData[prefix+"rt:bCurve"] =             serializeVector(labCurve.bcurve);

	prefix=baseKey+"rt:Vibrance/";
    if (!pedited || pedited->vibrance.enabled)          xmpData[prefix+"rt:Enabled"]=         vibrance.enabled;
    if (!pedited || pedited->vibrance.pastels)          xmpData[prefix+"rt:Pastels"]=         vibrance.pastels;
    if (!pedited || pedited->vibrance.saturated)        xmpData[prefix+"rt:Saturated"]=       vibrance.saturated;
    if (!pedited || pedited->vibrance.psthreshold)      xmpData[prefix+"rt:PSThreshold"]=     vibrance.psthreshold;
    if (!pedited || pedited->vibrance.protectskins)     xmpData[prefix+"rt:ProtectSkins"]=    vibrance.protectskins;
    if (!pedited || pedited->vibrance.avoidcolorshift)  xmpData[prefix+"rt:AvoidColorShift"]= vibrance.avoidcolorshift;
    if (!pedited || pedited->vibrance.pastsattog)       xmpData[prefix+"rt:PastSatTog"]=      vibrance.pastsattog;

	prefix=baseKey+"rt:Sharpening/";
    if (!pedited || pedited->sharpening.enabled)            xmpData[prefix+"rt:Enabled"] =             sharpening.enabled;
    if (!pedited || pedited->sharpening.method)             xmpData[prefix+"rt:Method"] =              sharpening.method;
    if (!pedited || pedited->sharpening.radius)             xmpData[prefix+"rt:Radius"]=               sharpening.radius;
    if (!pedited || pedited->sharpening.amount)             xmpData[prefix+"rt:Amount"]=               sharpening.amount;
    if (!pedited || pedited->sharpening.threshold)          xmpData[prefix+"rt:Threshold"]=            sharpening.threshold;
    if (!pedited || pedited->sharpening.edgesonly)          xmpData[prefix+"rt:OnlyEdges"]=            sharpening.edgesonly;
    if (!pedited || pedited->sharpening.edges_radius)       xmpData[prefix+"rt:EdgeDetectionRadius"]=  sharpening.edges_radius;
    if (!pedited || pedited->sharpening.edges_tolerance)    xmpData[prefix+"rt:EdgeTolerance"]=        sharpening.edges_tolerance;
    if (!pedited || pedited->sharpening.halocontrol)        xmpData[prefix+"rt:HaloControlEnabled"]=   sharpening.halocontrol;
    if (!pedited || pedited->sharpening.halocontrol_amount) xmpData[prefix+"rt:HaloControlAmount"]=    sharpening.halocontrol_amount;
    if (!pedited || pedited->sharpening.deconvradius)       xmpData[prefix+"rt:DeconvRadius"]=         sharpening.deconvradius;
    if (!pedited || pedited->sharpening.deconvamount)       xmpData[prefix+"rt:DeconvAmount"]=         sharpening.deconvamount;
    if (!pedited || pedited->sharpening.deconvdamping)      xmpData[prefix+"rt:DeconvDamping"]=        sharpening.deconvdamping;
    if (!pedited || pedited->sharpening.deconviter)         xmpData[prefix+"rt:DeconvIterations"]=     sharpening.deconviter;

	prefix=baseKey+"rt:SharpenEdge/";
    if (!pedited || pedited->sharpenEdge.enabled)       xmpData[prefix+"rt:Enabled"]=       sharpenEdge.enabled;
    if (!pedited || pedited->sharpenEdge.passes)        xmpData[prefix+"rt:Passes"]=        sharpenEdge.passes;
    if (!pedited || pedited->sharpenEdge.amount)        xmpData[prefix+"rt:ThreeChannels"]= sharpenEdge.threechannels;
    if (!pedited || pedited->sharpenEdge.threechannels) xmpData[prefix+"rt:Amount"]=        sharpenEdge.amount;

	prefix=baseKey+"rt:MicroContrast/";
    if (!pedited || pedited->sharpenMicro.enabled)      xmpData[prefix+"rt:Enabled"]=       sharpenMicro.enabled;
    if (!pedited || pedited->sharpenMicro.matrix)       xmpData[prefix+"rt:Uniformity"]=    sharpenMicro.uniformity;
    if (!pedited || pedited->sharpenMicro.amount)       xmpData[prefix+"rt:Matrix"]=        sharpenMicro.matrix;
    if (!pedited || pedited->sharpenMicro.uniformity)   xmpData[prefix+"rt:Amount"]=        sharpenMicro.amount;

	prefix=baseKey+"rt:WhiteBalance/";
    if (!pedited || pedited->wb.method)      xmpData[prefix+"rt:Mode"]= wb.method;
    if (!pedited || pedited->wb.temperature) xmpData[prefix+"rt:Temperature"]= wb.temperature;
    if (!pedited || pedited->wb.green)       xmpData[prefix+"rt:Green"]=       wb.green;

    prefix=baseKey+"rt:ImpulseDenoise/";
    if (!pedited || pedited->impulseDenoise.enabled) xmpData[prefix+"rt:Enabled"]=   impulseDenoise.enabled;
    if (!pedited || pedited->impulseDenoise.thresh)  xmpData[prefix+"rt:Threshold"]= impulseDenoise.thresh;

    prefix=baseKey+"rt:Defringe/";
    if (!pedited || pedited->defringe.enabled)       xmpData[prefix+"rt:Enabled"]=   defringe.enabled;
    if (!pedited || pedited->defringe.radius)        xmpData[prefix+"rt:Radius"]=    defringe.radius;
    if (!pedited || pedited->defringe.threshold)     xmpData[prefix+"rt:Threshold"]=	defringe.threshold;

    prefix=baseKey+"rt:PyramidDenoise/";
    if (!pedited || pedited->dirpyrDenoise.enabled) xmpData[prefix+"rt:Enabled"]= dirpyrDenoise.enabled;
    if (!pedited || pedited->dirpyrDenoise.luma)    xmpData[prefix+"rt:Luma"]=    dirpyrDenoise.luma;
    if (!pedited || pedited->dirpyrDenoise.chroma)  xmpData[prefix+"rt:Chroma"]=  dirpyrDenoise.chroma;
    if (!pedited || pedited->dirpyrDenoise.gamma)   xmpData[prefix+"rt:Gamma"]=   dirpyrDenoise.gamma;

    prefix=baseKey+"rt:ShadowHighlights/";
    if (!pedited || pedited->sh.enabled)       xmpData[prefix+"rt:Enabled"]=               sh.enabled;
    if (!pedited || pedited->sh.hq)            xmpData[prefix+"rt:HighQuality"]=           sh.hq;
    if (!pedited || pedited->sh.highlights)    xmpData[prefix+"rt:Highlights"]=            sh.highlights;
    if (!pedited || pedited->sh.htonalwidth)   xmpData[prefix+"rt:HighlightTonalWidth"]=   sh.htonalwidth;
    if (!pedited || pedited->sh.shadows)       xmpData[prefix+"rt:Shadows"]=               sh.shadows;
    if (!pedited || pedited->sh.stonalwidth)   xmpData[prefix+"rt:ShadowTonalWidth"]=      sh.stonalwidth;
    if (!pedited || pedited->sh.localcontrast) xmpData[prefix+"rt:LocalContrast"]=         sh.localcontrast;
    if (!pedited || pedited->sh.radius)        xmpData[prefix+"rt:Radius"]=                sh.radius;

    prefix=baseKey+"rt:Crop/";
    if (!pedited || pedited->crop.enabled)     xmpData[prefix+"rt:Enabled"]=     crop.enabled;
    if (!pedited || pedited->crop.x)           xmpData[prefix+"rt:X"]=           crop.x;
    if (!pedited || pedited->crop.y)           xmpData[prefix+"rt:Y"]=           crop.y;
    if (!pedited || pedited->crop.w)           xmpData[prefix+"rt:Width"]=       crop.w;
    if (!pedited || pedited->crop.h)           xmpData[prefix+"rt:Height"]=      crop.h;
    /*
    if (!pedited || pedited->crop.fixratio)    xmpData[prefix+"rt:FixedRatio"]=  crop.fixratio;
    if (!pedited || pedited->crop.ratio)       xmpData[prefix+"rt:Ratio"]=       crop.ratio;
    if (!pedited || pedited->crop.orientation) xmpData[prefix+"rt:Orientation"]= crop.orientation;
    if (!pedited || pedited->crop.guide)       xmpData[prefix+"rt:Guide"]=       crop.guide;
    */

    prefix=baseKey+"rt:CoarseGeo/";
    if (!pedited || pedited->coarse.rotate)    xmpData[prefix+"rt:RotationDegree"]=  coarse.rotate;
    if (!pedited || pedited->coarse.hflip)     xmpData[prefix+"rt:HorizontalFlip"]=  coarse.hflip;
    if (!pedited || pedited->coarse.vflip)     xmpData[prefix+"rt:VerticalFlip"]=    coarse.vflip;

    prefix=baseKey+"rt:EPD/";
    if (!pedited || pedited->edgePreservingDecompositionUI.enabled)             xmpData[prefix+"rt:Enabled"]= edgePreservingDecompositionUI.enabled;
    if (!pedited || pedited->edgePreservingDecompositionUI.Strength)            xmpData[prefix+"rt:Strength"]= edgePreservingDecompositionUI.Strength;
    if (!pedited || pedited->edgePreservingDecompositionUI.EdgeStopping)        xmpData[prefix+"rt:EdgeStopping"]= edgePreservingDecompositionUI.EdgeStopping;
    if (!pedited || pedited->edgePreservingDecompositionUI.Scale)               xmpData[prefix+"rt:Scale"] = edgePreservingDecompositionUI.Scale;
    if (!pedited || pedited->edgePreservingDecompositionUI.ReweightingIterates) xmpData[prefix+"rt:ReweightingIterates"] = edgePreservingDecompositionUI.ReweightingIterates;


    prefix=baseKey+"rt:Geometry/";
    //xmpData[prefix+"rt:Enabled"]=    geo.enabled;
    if (!pedited || pedited->commonTrans.autofill)   xmpData[prefix+"rt:AutoFill"]= commonTrans.autofill;
    if (!pedited || pedited->rotate.degree)          xmpData[prefix+"rt:RotationDegree"]= rotate.degree;
    if (!pedited || pedited->distortion.amount)      xmpData[prefix+"rt:DistortionAmount"]= distortion.amount;
    if (!pedited || pedited->perspective.horizontal) xmpData[prefix+"rt:HorizontalPerspective"]= perspective.horizontal;
    if (!pedited || pedited->perspective.vertical)   xmpData[prefix+"rt:VerticalPerspective"]=   perspective.vertical;

    prefix=baseKey+"rt:CACorrection/";
    //xmpData[prefix+"rt:Enabled"]=    cacorrection.enabled;
    if (!pedited || pedited->cacorrection.red)       xmpData[prefix+"rt:Red"]=  cacorrection.red;
    if (!pedited || pedited->cacorrection.blue)      xmpData[prefix+"rt:Blue"]= cacorrection.blue;

    prefix=baseKey+"rt:Vignetting/";
    //xmpData[prefix+"rt:Enabled"]= vignetting.enabled;
    if (!pedited || pedited->vignetting.amount)      xmpData[prefix+"rt:Amount"] = vignetting.amount;
    if (!pedited || pedited->vignetting.radius)      xmpData[prefix+"rt:Radius"] = vignetting.radius;
    if (!pedited || pedited->vignetting.strength)    xmpData[prefix+"rt:Strength"]= vignetting.strength;
    if (!pedited || pedited->vignetting.centerX)     xmpData[prefix+"rt:CenterX"] = vignetting.centerX;
    if (!pedited || pedited->vignetting.centerY)     xmpData[prefix+"rt:CenterY"] = vignetting.centerY;

    prefix=baseKey+"rt:HLRecovery/";
    if (!pedited || pedited->hlrecovery.enabled)     xmpData[prefix+"rt:Enabled"]=  hlrecovery.enabled;
    if (!pedited || pedited->hlrecovery.method)      xmpData[prefix+"rt:Method"]=   hlrecovery.method;

    prefix=baseKey+"rt:Resize/";
    if (!pedited || pedited->resize.enabled)         xmpData[prefix+"rt:Enabled"]=   resize.enabled;
    if (!pedited || pedited->resize.scale)           xmpData[prefix+"rt:Scale"]  =   resize.scale;
    if (!pedited || pedited->resize.appliesTo)       xmpData[prefix+"rt:AppliesTo"]= resize.appliesTo;
    if (!pedited || pedited->resize.method)          xmpData[prefix+"rt:Method"]=    resize.method;
    if (!pedited || pedited->resize.dataspec)        xmpData[prefix+"rt:DataSpecified"]=  resize.dataspec;
    if (!pedited || pedited->resize.width)           xmpData[prefix+"rt:Width"] =    resize.width;
    if (!pedited || pedited->resize.height)          xmpData[prefix+"rt:Height"] =   resize.height;

    prefix=baseKey+"rt:ColorManagement/";
    // save color management settings
    if (!pedited || pedited->icm.input)              xmpData[prefix+"rt:InputProfile"] =     icm.input;
    if (!pedited || pedited->icm.blendCMSMatrix)     xmpData[prefix+"rt:BlendCMSMatrix"] =   icm.blendCMSMatrix;
    if (!pedited || pedited->icm.preferredProfile)   xmpData[prefix+"rt:PreferredProfile"] = icm.preferredProfile;
    if (!pedited || pedited->icm.working)            xmpData[prefix+"rt:WorkingProfile"] =   icm.working;
    if (!pedited || pedited->icm.output)             xmpData[prefix+"rt:OutputProfile"] =    icm.output;
    if (!pedited || pedited->icm.gamma)              xmpData[prefix+"rt:FreeGamma"] =        icm.gamma;
    if (!pedited || pedited->icm.freegamma)          xmpData[prefix+"rt:FreeGammaEnabled"] = icm.freegamma;
    if (!pedited || pedited->icm.gampos)             xmpData[prefix+"rt:GammaValue"] =       icm.gampos;
    if (!pedited || pedited->icm.slpos)              xmpData[prefix+"rt:GammaSlope"] =       icm.slpos;

    prefix=baseKey+"rt:DirectionalPyramidEqualizer/";
    if (!pedited || pedited->dirpyrequalizer.enabled) xmpData[prefix+"rt:Enabled"] = dirpyrequalizer.enabled ;
    if (!pedited || pedited->dirpyrequalizer.mult[0] || pedited->dirpyrequalizer.mult[1] || pedited->dirpyrequalizer.mult[2]
      || pedited->dirpyrequalizer.mult[3] || pedited->dirpyrequalizer.mult[4])
    	xmpData[prefix+"rt:Coeff"] =    serializeArray( dirpyrequalizer.mult,5);

    prefix=baseKey+"rt:HSVEqualizer/";
    //xmpData[prefix+"rt:Enabled"]= hsvequalizer.enabled;
    if (!pedited || pedited->hsvequalizer.hcurve)    xmpData[prefix+"rt:HCurve"] = serializeVector(hsvequalizer.hcurve);
    if (!pedited || pedited->hsvequalizer.scurve)    xmpData[prefix+"rt:SCurve"] = serializeVector(hsvequalizer.scurve);
    if (!pedited || pedited->hsvequalizer.vcurve)    xmpData[prefix+"rt:VCurve"] = serializeVector(hsvequalizer.vcurve);

    prefix=baseKey+"rt:RGBEqualizer/";
    //xmpData[prefix+"rt:Enabled"]= hsvequalizer.enabled;
    if (!pedited || pedited->rgbCurves.rcurve)       xmpData[prefix+"rt:RCurve"] = serializeVector(rgbCurves.rcurve);
    if (!pedited || pedited->rgbCurves.gcurve)       xmpData[prefix+"rt:GCurve"] = serializeVector(rgbCurves.gcurve);
    if (!pedited || pedited->rgbCurves.bcurve)       xmpData[prefix+"rt:BCurve"] = serializeVector(rgbCurves.bcurve);

    prefix=baseKey+"rt:RawArithmetic/";
    //xmpData[prefix+"rt:Enabled"]=
    if (!pedited || pedited->raw.darkFrame)          xmpData[prefix+"rt:DarkFrameFile"]=       raw.dark_frame;
    if (!pedited || pedited->raw.dfAuto)             xmpData[prefix+"rt:DarkFrameAutoSelect"]= raw.df_autoselect ;
    if (!pedited || pedited->raw.ff_file)            xmpData[prefix+"rt:FlatFieldFile"]=       raw.ff_file ;
    if (!pedited || pedited->raw.ff_AutoSelect)      xmpData[prefix+"rt:FlatFieldAutoSelect"]= raw.ff_AutoSelect ;
    if (!pedited || pedited->raw.ff_BlurRadius)      xmpData[prefix+"rt:FlatFieldBlurRadius"]= raw.ff_BlurRadius ;
    if (!pedited || pedited->raw.ff_BlurType)        xmpData[prefix+"rt:FlatFieldBlurType"]=   raw.ff_BlurType ;

    prefix=baseKey+"rt:RawCACorrection/";
    //xmpData[prefix+"rt:Enabled"]=
    if (!pedited || pedited->raw.caCorrection)       xmpData[prefix+"rt:Auto"]=        raw.ca_autocorrect ;
    if (!pedited || pedited->raw.caRed)              xmpData[prefix+"rt:Red"] =        raw.cared;
    if (!pedited || pedited->raw.caBlue)             xmpData[prefix+"rt:Blue"]=        raw.cablue;

    prefix=baseKey+"rt:HotDeadPixelCorrection/";
    if (!pedited || pedited->raw.hotDeadPixelFilter) xmpData[prefix+"rt:Enabled"]=     raw.hotdeadpix_filt;
    if (!pedited || pedited->raw.hotDeadPixelThresh) xmpData[prefix+"rt:Threshold"]=   raw.hotdeadpix_thresh;

    prefix=baseKey+"rt:RawDenoise/";
    //xmpData[prefix+"rt:Enabled"]=
    if (!pedited || pedited->raw.linenoise)          xmpData[prefix+"rt:LineDenoise"]=       raw.linenoise;

    prefix=baseKey+"rt:Demosaicing/";
    if (!pedited || pedited->raw.greenEq)            xmpData[prefix+"rt:GreenEqThreshold"] = raw.greenthresh;
    if (!pedited || pedited->raw.ccSteps)            xmpData[prefix+"rt:CcSteps"] =          raw.ccSteps;
    if (!pedited || pedited->raw.dmethod)            xmpData[prefix+"rt:Method"]  =          raw.dmethod;
    if (!pedited || pedited->raw.dcbIterations)      xmpData[prefix+"rt:DCBIterations"] =    raw.dcb_iterations;
    if (!pedited || pedited->raw.dcbEnhance)         xmpData[prefix+"rt:DCBEnhance"] =       raw.dcb_enhance;
    if (!pedited || pedited->raw.allEnhance)         xmpData[prefix+"rt:Enhance"]=           raw.all_enhance;

    prefix=baseKey+"rt:RawExposure/";
    if (!pedited || pedited->raw.exPos)              xmpData[prefix+"rt:Exposure"] =  raw.expos;
    if (!pedited || pedited->raw.exPreser)           xmpData[prefix+"rt:HLPreserving"] = raw.preser;
    if (!pedited || pedited->raw.exBlackzero)        xmpData[prefix+"rt:Blackzero"] = raw.blackzero;
    if (!pedited || pedited->raw.exBlackone)         xmpData[prefix+"rt:Blackone"] = raw.blackone;
    if (!pedited || pedited->raw.exBlacktwo)         xmpData[prefix+"rt:Blacktwo"] = raw.blacktwo;
    if (!pedited || pedited->raw.exBlackthree)       xmpData[prefix+"rt:Blackthree"] = raw.blackthree;
    if (!pedited || pedited->raw.exTwoGreen)         xmpData[prefix+"rt:TwoGreen"] = raw.twogreen;
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
	
	if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:EPD")) != xmpData.end()){
	    prefix=baseKey+"rt:EPD/rt:";
	    readVarFromXmp( xmpData, prefix+kXmpEnabled,edgePreservingDecompositionUI.enabled );
        readVarFromXmp( xmpData, prefix+"Strength", edgePreservingDecompositionUI.Strength );
        readVarFromXmp( xmpData, prefix+"EdgeStopping", edgePreservingDecompositionUI.EdgeStopping);
        readVarFromXmp( xmpData, prefix+"Scale", edgePreservingDecompositionUI.Scale );
        readVarFromXmp( xmpData, prefix+"ReweightingIterates", edgePreservingDecompositionUI.ReweightingIterates);
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
    try{
		if (0 != Exiv2::XmpParser::decode(xmpData,xmpPacket) ) {
			return 1;
		}
    }catch(Exiv2::Error &e){
    	printf("Exception in parser: %s\n", e.what());
    	return 2;
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

int ProcParams::load (Glib::ustring fname, ParamsEdited* pedited, int *rank) {

    SafeKeyFile keyFile;
    try {
        //setDefaults ();
        if (pedited)
            pedited->set(false);

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
    //if (keyFile.has_key ("General", "Rank"))        { rank       = keyFile.get_integer ("General", "Rank"); if (pedited) pedited->general.rank = true; }
    //if (keyFile.has_key ("General", "ColorLabel"))  { colorlabel = keyFile.get_integer ("General", "ColorLabel"); if (pedited) pedited->general.colorlabel = true; }
    //if (keyFile.has_key ("General", "InTrash"))     { inTrash    = keyFile.get_boolean ("General", "InTrash"); if (pedited) pedited->general.intrash = true; }

	// is now:
	
	if( rank!= NULL){
		*rank = 0;
		if (keyFile.has_key ("General", "InTrash")){
			if( keyFile.get_boolean ("General", "InTrash") )
				*rank = -1;
		}
		if( *rank != -1  && keyFile.has_key ("General", "Rank"))
			*rank = keyFile.get_integer ("General", "Rank");
	}
}

if (keyFile.has_group ("Exposure")) {
    if (ppVersion<PPVERSION_AEXP)
        toneCurve.autoexp = false; // prevent execution of autoexp when opening file created with earlier verions of autoexp algorithm
    else
        if (keyFile.has_key ("Exposure", "Auto"))           { toneCurve.autoexp       = keyFile.get_boolean ("Exposure", "Auto"); if (pedited) pedited->toneCurve.autoexp = true; }

    if (keyFile.has_key ("Exposure", "Clip"))           { toneCurve.clip          = keyFile.get_double  ("Exposure", "Clip"); if (pedited) pedited->toneCurve.clip = true; }
    if (keyFile.has_key ("Exposure", "Compensation"))   { toneCurve.expcomp       = keyFile.get_double  ("Exposure", "Compensation"); if (pedited) pedited->toneCurve.expcomp = true; }
    if (keyFile.has_key ("Exposure", "Brightness"))     { toneCurve.brightness    = keyFile.get_integer ("Exposure", "Brightness"); if (pedited) pedited->toneCurve.brightness = true; }
    if (keyFile.has_key ("Exposure", "Contrast"))       { toneCurve.contrast      = keyFile.get_integer ("Exposure", "Contrast"); if (pedited) pedited->toneCurve.contrast = true; }
    if (keyFile.has_key ("Exposure", "Saturation"))     { toneCurve.saturation    = keyFile.get_integer ("Exposure", "Saturation"); if (pedited) pedited->toneCurve.saturation = true; }
    if (keyFile.has_key ("Exposure", "Black"))          { toneCurve.black         = keyFile.get_integer ("Exposure", "Black"); if (pedited) pedited->toneCurve.black = true; }
    if (keyFile.has_key ("Exposure", "HighlightCompr")) { toneCurve.hlcompr       = keyFile.get_integer ("Exposure", "HighlightCompr"); if (pedited) pedited->toneCurve.hlcompr = true; }
    if (keyFile.has_key ("Exposure", "HighlightComprThreshold")) { toneCurve.hlcomprthresh = keyFile.get_integer ("Exposure", "HighlightComprThreshold"); if (pedited) pedited->toneCurve.hlcomprthresh = true; }
    if (keyFile.has_key ("Exposure", "ShadowCompr"))    { toneCurve.shcompr       = keyFile.get_integer ("Exposure", "ShadowCompr"); if (pedited) pedited->toneCurve.shcompr = true; }
    if (toneCurve.shcompr > 100) toneCurve.shcompr = 100; // older pp3 files can have values above 100.
    if (ppVersion>200)
    if (keyFile.has_key ("Exposure", "Curve"))          { toneCurve.curve         = keyFile.get_double_list ("Exposure", "Curve"); if (pedited) pedited->toneCurve.curve = true; }
}

    // load channel mixer curve
if (keyFile.has_group ("Channel Mixer")) {    
    if (keyFile.has_key ("Channel Mixer", "Red") && keyFile.has_key ("Channel Mixer", "Green") && keyFile.has_key ("Channel Mixer", "Blue")) {
        if (pedited) {
            pedited->chmixer.red[0]   = pedited->chmixer.red[1]   = pedited->chmixer.red[2] = true;
            pedited->chmixer.green[0] = pedited->chmixer.green[1] = pedited->chmixer.green[2] = true;
            pedited->chmixer.blue[0]  = pedited->chmixer.blue[1]  = pedited->chmixer.blue[2] = true;
        }
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
    if (keyFile.has_key ("Luminance Curve", "Brightness"))     { labCurve.brightness  = keyFile.get_integer  ("Luminance Curve", "Brightness"); if (pedited) pedited->labCurve.brightness = true; }
    if (keyFile.has_key ("Luminance Curve", "Contrast"))       { labCurve.contrast    = keyFile.get_integer ("Luminance Curve", "Contrast"); if (pedited) pedited->labCurve.contrast = true; }
	if (keyFile.has_key ("Luminance Curve", "Saturation"))      { labCurve.saturation = keyFile.get_integer ("Luminance Curve", "Saturation"); if (pedited) pedited->labCurve.saturation = true; }
	if (keyFile.has_key ("Luminance Curve", "AvoidColorClipping"))  { labCurve.avoidclip                = keyFile.get_boolean ("Luminance Curve", "AvoidColorClipping"); if (pedited) pedited->labCurve.avoidclip = true; }
    if (keyFile.has_key ("Luminance Curve", "SaturationLimiter"))   { labCurve.enable_saturationlimiter = keyFile.get_boolean ("Luminance Curve", "SaturationLimiter"); if (pedited) pedited->labCurve.enable_saturationlimiter = true; }
    if (keyFile.has_key ("Luminance Curve", "SaturationLimit"))     { labCurve.saturationlimit          = keyFile.get_double  ("Luminance Curve", "SaturationLimit");	 if (pedited) pedited->labCurve.saturationlimit = true; }
	if (keyFile.has_key ("Luminance Curve", "LCurve"))          { labCurve.lcurve = keyFile.get_double_list ("Luminance Curve", "LCurve"); if (pedited) pedited->labCurve.lcurve = true; }
	if (keyFile.has_key ("Luminance Curve", "aCurve"))          { labCurve.acurve = keyFile.get_double_list ("Luminance Curve", "aCurve"); if (pedited) pedited->labCurve.acurve = true; }
	if (keyFile.has_key ("Luminance Curve", "bCurve"))          { labCurve.bcurve = keyFile.get_double_list ("Luminance Curve", "bCurve"); if (pedited) pedited->labCurve.bcurve = true; }
}

    // load sharpening
if (keyFile.has_group ("Sharpening")) {
    if (keyFile.has_key ("Sharpening", "Enabled"))              { sharpening.enabled          = keyFile.get_boolean ("Sharpening", "Enabled"); if (pedited) pedited->sharpening.enabled = true; }
    if (keyFile.has_key ("Sharpening", "Radius"))               { sharpening.radius           = keyFile.get_double  ("Sharpening", "Radius"); if (pedited) pedited->sharpening.radius = true; }
    if (keyFile.has_key ("Sharpening", "Amount"))               { sharpening.amount           = keyFile.get_integer ("Sharpening", "Amount"); if (pedited) pedited->sharpening.amount = true; }
    if (keyFile.has_key ("Sharpening", "Threshold"))            { sharpening.threshold        = keyFile.get_integer ("Sharpening", "Threshold"); if (pedited) pedited->sharpening.threshold = true; }
    if (keyFile.has_key ("Sharpening", "OnlyEdges"))            { sharpening.edgesonly        = keyFile.get_boolean ("Sharpening", "OnlyEdges"); if (pedited) pedited->sharpening.edgesonly = true; }
    if (keyFile.has_key ("Sharpening", "EdgedetectionRadius"))  { sharpening.edges_radius     = keyFile.get_double  ("Sharpening", "EdgedetectionRadius"); if (pedited) pedited->sharpening.edges_radius = true; }
    if (keyFile.has_key ("Sharpening", "EdgeTolerance"))        { sharpening.edges_tolerance  = keyFile.get_integer ("Sharpening", "EdgeTolerance"); if (pedited) pedited->sharpening.edges_tolerance = true; }
    if (keyFile.has_key ("Sharpening", "HalocontrolEnabled"))   { sharpening.halocontrol      = keyFile.get_boolean ("Sharpening", "HalocontrolEnabled"); if (pedited) pedited->sharpening.halocontrol = true; }
    if (keyFile.has_key ("Sharpening", "HalocontrolAmount"))    { sharpening.halocontrol_amount = keyFile.get_integer ("Sharpening", "HalocontrolAmount"); if (pedited) pedited->sharpening.halocontrol_amount = true; }
    if (keyFile.has_key ("Sharpening", "Method"))               { sharpening.method           = keyFile.get_string  ("Sharpening", "Method"); if (pedited) pedited->sharpening.method = true; }
    if (keyFile.has_key ("Sharpening", "DeconvRadius"))         { sharpening.deconvradius     = keyFile.get_double  ("Sharpening", "DeconvRadius"); if (pedited) pedited->sharpening.deconvradius = true; }
    if (keyFile.has_key ("Sharpening", "DeconvAmount"))         { sharpening.deconvamount     = keyFile.get_integer ("Sharpening", "DeconvAmount"); if (pedited) pedited->sharpening.deconvamount = true; }
    if (keyFile.has_key ("Sharpening", "DeconvDamping"))        { sharpening.deconvdamping    = keyFile.get_integer ("Sharpening", "DeconvDamping"); if (pedited) pedited->sharpening.deconvdamping = true; }
    if (keyFile.has_key ("Sharpening", "DeconvIterations"))     { sharpening.deconviter       = keyFile.get_integer ("Sharpening", "DeconvIterations"); if (pedited) pedited->sharpening.deconviter = true; }
}

    // load edge sharpening
if (keyFile.has_group ("SharpenEdge")) {
    if (keyFile.has_key ("SharpenEdge", "Enabled"))             { sharpenEdge.enabled         = keyFile.get_boolean ("SharpenEdge", "Enabled"); if (pedited) pedited->sharpenEdge.enabled = true; }
    if (keyFile.has_key ("SharpenEdge", "Passes"))              { sharpenEdge.passes          = keyFile.get_integer  ("SharpenEdge", "Passes"); if (pedited) pedited->sharpenEdge.passes = true; }
    if (keyFile.has_key ("SharpenEdge", "Strength"))            { sharpenEdge.amount          = keyFile.get_double  ("SharpenEdge", "Strength"); if (pedited) pedited->sharpenEdge.amount = true; }
    if (keyFile.has_key ("SharpenEdge", "ThreeChannels"))       { sharpenEdge.threechannels   = keyFile.get_boolean ("SharpenEdge", "ThreeChannels"); if (pedited) pedited->sharpenEdge.threechannels = true; }
}

    // load micro-contrast sharpening
if (keyFile.has_group ("SharpenMicro")) {
    if (keyFile.has_key ("SharpenMicro", "Enabled"))            { sharpenMicro.enabled        = keyFile.get_boolean ("SharpenMicro", "Enabled"); if (pedited) pedited->sharpenMicro.enabled = true; }
    if (keyFile.has_key ("SharpenMicro", "Matrix"))             { sharpenMicro.matrix         = keyFile.get_boolean ("SharpenMicro", "Matrix"); if (pedited) pedited->sharpenMicro.matrix = true; }
    if (keyFile.has_key ("SharpenMicro", "Strength"))           { sharpenMicro.amount         = keyFile.get_double  ("SharpenMicro", "Strength"); if (pedited) pedited->sharpenMicro.amount = true; }
    if (keyFile.has_key ("SharpenMicro", "Uniformity"))         { sharpenMicro.uniformity     = keyFile.get_double  ("SharpenMicro", "Uniformity"); if (pedited) pedited->sharpenMicro.uniformity = true; }
}

    // load vibrance
if (keyFile.has_group ("Vibrance")) {
    if (keyFile.has_key ("Vibrance", "Enabled"))                { vibrance.enabled            = keyFile.get_boolean ("Vibrance", "Enabled"); if (pedited) pedited->vibrance.enabled = true; }
    if (keyFile.has_key ("Vibrance", "Pastels"))                { vibrance.pastels            = keyFile.get_integer ("Vibrance", "Pastels"); if (pedited) pedited->vibrance.pastels = true; }
    if (keyFile.has_key ("Vibrance", "Saturated"))              { vibrance.saturated          = keyFile.get_integer ("Vibrance", "Saturated"); if (pedited) pedited->vibrance.saturated = true; }
    if (keyFile.has_key ("Vibrance", "PSThreshold"))            { vibrance.psthreshold        = keyFile.get_integer ("Vibrance", "PSThreshold"); if (pedited) pedited->vibrance.psthreshold = true; }
    if (keyFile.has_key ("Vibrance", "ProtectSkins"))           { vibrance.protectskins       = keyFile.get_boolean ("Vibrance", "ProtectSkins"); if (pedited) pedited->vibrance.protectskins = true; }
    if (keyFile.has_key ("Vibrance", "AvoidColorShift"))        { vibrance.avoidcolorshift    = keyFile.get_boolean ("Vibrance", "AvoidColorShift"); if (pedited) pedited->vibrance.avoidcolorshift = true; }
    if (keyFile.has_key ("Vibrance", "PastSatTog"))             { vibrance.pastsattog         = keyFile.get_boolean ("Vibrance", "PastSatTog"); if (pedited) pedited->vibrance.pastsattog = true; }
}

    // load colorBoost
/*if (keyFile.has_group ("Color Boost")) {
    if (keyFile.has_key ("Color Boost", "Amount"))              { colorBoost.amount           = keyFile.get_integer ("Color Boost", "Amount"); if (pedited) pedited->colorBoost.amount = true; }
    else {
        int a=0, b=0;
        if (keyFile.has_key ("Color Boost", "ChannelA"))        { a                           = keyFile.get_integer ("Color Boost", "ChannelA"); }
        if (keyFile.has_key ("Color Boost", "ChannelB"))        { b                           = keyFile.get_integer ("Color Boost", "ChannelB"); }
        colorBoost.amount = (a+b) / 2;
        if (pedited) pedited->colorBoost.amount = true;
    }   
    if (keyFile.has_key ("Color Boost", "AvoidColorClipping"))  { colorBoost.avoidclip               = keyFile.get_boolean ("Color Boost", "AvoidColorClipping"); if (pedited) pedited->colorBoost.avoidclip = true; }
    if (keyFile.has_key ("Color Boost", "SaturationLimiter"))   { colorBoost.enable_saturationlimiter= keyFile.get_boolean ("Color Boost", "SaturationLimiter"); if (pedited) pedited->colorBoost.enable_saturationlimiter = true; }
    if (keyFile.has_key ("Color Boost", "SaturationLimit"))     { colorBoost.saturationlimit         = keyFile.get_double  ("Color Boost", "SaturationLimit"); if (pedited) pedited->colorBoost.saturationlimit = true; }
}*/

    // load wb
if (keyFile.has_group ("White Balance")) {
    if (keyFile.has_key ("White Balance", "Setting"))     { wb.method         = keyFile.get_string ("White Balance", "Setting"); if (pedited) pedited->wb.method = true; }
    if (keyFile.has_key ("White Balance", "Temperature")) { wb.temperature    = keyFile.get_integer ("White Balance", "Temperature"); if (pedited) pedited->wb.temperature = true; }
    if (keyFile.has_key ("White Balance", "Green"))       { wb.green          = keyFile.get_double ("White Balance", "Green"); if (pedited) pedited->wb.green = true; }
}

    // load colorShift
/*if (keyFile.has_group ("Color Shift")) {
    if (keyFile.has_key ("Color Shift", "ChannelA")) { colorShift.a = keyFile.get_double ("Color Shift", "ChannelA"); if (pedited) pedited->colorShift.a = true; }
    if (keyFile.has_key ("Color Shift", "ChannelB")) { colorShift.b = keyFile.get_double ("Color Shift", "ChannelB"); if (pedited) pedited->colorShift.b = true; }
}*/
		
// load defringe
if (keyFile.has_group ("Defringing")) {
	if (keyFile.has_key ("Defringing", "Enabled"))        { defringe.enabled   = keyFile.get_boolean ("Defringing", "Enabled"); if (pedited) pedited->defringe.enabled = true; }
	if (keyFile.has_key ("Defringing", "Radius"))         { defringe.radius    = keyFile.get_double  ("Defringing", "Radius"); if (pedited) pedited->defringe.radius = true; }
	if (keyFile.has_key ("Defringing", "Threshold"))      { defringe.threshold = keyFile.get_integer ("Defringing", "Threshold"); if (pedited) pedited->defringe.threshold = true; }
}
		
	// load impulseDenoise
if (keyFile.has_group ("Impulse Denoising")) {
	if (keyFile.has_key ("Impulse Denoising", "Enabled"))   { impulseDenoise.enabled = keyFile.get_boolean ("Impulse Denoising", "Enabled"); if (pedited) pedited->impulseDenoise.enabled = true; }
	if (keyFile.has_key ("Impulse Denoising", "Threshold")) { impulseDenoise.thresh  = keyFile.get_integer ("Impulse Denoising", "Threshold"); if (pedited) pedited->impulseDenoise.thresh = true; }
}
		
	// load dirpyrDenoise
if (keyFile.has_group ("Directional Pyramid Denoising")) {
	if (keyFile.has_key ("Directional Pyramid Denoising", "Enabled"))    { dirpyrDenoise.enabled = keyFile.get_boolean ("Directional Pyramid Denoising", "Enabled"); if (pedited) pedited->dirpyrDenoise.enabled = true; }
	if (keyFile.has_key ("Directional Pyramid Denoising", "Luma"))       { dirpyrDenoise.luma    = keyFile.get_integer ("Directional Pyramid Denoising", "Luma"); if (pedited) pedited->dirpyrDenoise.luma = true; }
	if (keyFile.has_key ("Directional Pyramid Denoising", "Chroma"))     { dirpyrDenoise.chroma  = keyFile.get_integer ("Directional Pyramid Denoising", "Chroma"); if (pedited) pedited->dirpyrDenoise.chroma = true; }
	if (keyFile.has_key ("Directional Pyramid Denoising", "Gamma"))      { dirpyrDenoise.gamma  = keyFile.get_double ("Directional Pyramid Denoising", "Gamma"); if (pedited) pedited->dirpyrDenoise.gamma = true; }
}

//Load EPD.
if (keyFile.has_group ("EPD")) {
	if(keyFile.has_key("EPD", "Enabled"))             { edgePreservingDecompositionUI.enabled = keyFile.get_boolean ("EPD", "Enabled"); if (pedited) pedited->edgePreservingDecompositionUI.enabled = true; }
	if(keyFile.has_key("EPD", "Strength"))            { edgePreservingDecompositionUI.Strength = keyFile.get_double ("EPD", "Strength"); if (pedited) pedited->edgePreservingDecompositionUI.Strength = true; }
	if(keyFile.has_key("EPD", "EdgeStopping"))        { edgePreservingDecompositionUI.EdgeStopping = keyFile.get_double ("EPD", "EdgeStopping"); if (pedited) pedited->edgePreservingDecompositionUI.EdgeStopping = true; }
	if(keyFile.has_key("EPD", "Scale"))               { edgePreservingDecompositionUI.Scale = keyFile.get_double ("EPD", "Scale"); if (pedited) pedited->edgePreservingDecompositionUI.Scale = true; }
	if(keyFile.has_key("EPD", "ReweightingIterates")) { edgePreservingDecompositionUI.ReweightingIterates = keyFile.get_integer ("EPD", "ReweightingIterates"); if (pedited) pedited->edgePreservingDecompositionUI.ReweightingIterates = true; }
}
  
    // load lumaDenoise
/*if (keyFile.has_group ("Luminance Denoising")) {
    if (keyFile.has_key ("Luminance Denoising", "Enabled"))        { lumaDenoise.enabled       = keyFile.get_boolean ("Luminance Denoising", "Enabled"); if (pedited) pedited->lumaDenoise.enabled = true; }
    if (keyFile.has_key ("Luminance Denoising", "Radius"))         { lumaDenoise.radius        = keyFile.get_double  ("Luminance Denoising", "Radius"); if (pedited) pedited->lumaDenoise.radius = true; }
    if (keyFile.has_key ("Luminance Denoising", "EdgeTolerance"))  { lumaDenoise.edgetolerance = keyFile.get_integer ("Luminance Denoising", "EdgeTolerance"); if (pedited) pedited->lumaDenoise.edgetolerance = true; }
}*/

    // load colorDenoise
/*if (keyFile.has_group ("Chrominance Denoising")) {
    if (keyFile.has_key ("Chrominance Denoising", "Enabled"))      { colorDenoise.enabled       = keyFile.get_boolean 	("Chrominance Denoising", "Enabled"); if (pedited) pedited->colorDenoise.enabled = true; }
    // WARNING: radius doesn't exist anymore; is there any compatibility issue that require to keep the following line?
    if (keyFile.has_key ("Chrominance Denoising", "Radius"))       { colorDenoise.amount        = 10*keyFile.get_double ("Chrominance Denoising", "Radius"); }
    if (keyFile.has_key ("Chrominance Denoising", "Amount"))       { colorDenoise.amount        = keyFile.get_integer  	("Chrominance Denoising", "Amount"); if (pedited) pedited->colorDenoise.amount = true; }
}*/

    // load sh
if (keyFile.has_group ("Shadows & Highlights")) {
    if (keyFile.has_key ("Shadows & Highlights", "Enabled"))               { sh.enabled       = keyFile.get_boolean ("Shadows & Highlights", "Enabled"); if (pedited) pedited->sh.enabled = true; }
    if (keyFile.has_key ("Shadows & Highlights", "HighQuality"))           { sh.hq            = keyFile.get_boolean ("Shadows & Highlights", "HighQuality"); if (pedited) pedited->sh.hq = true; }
    if (keyFile.has_key ("Shadows & Highlights", "Highlights"))            { sh.highlights    = keyFile.get_integer ("Shadows & Highlights", "Highlights"); if (pedited) pedited->sh.highlights = true; }
    if (keyFile.has_key ("Shadows & Highlights", "HighlightTonalWidth"))   { sh.htonalwidth   = keyFile.get_integer ("Shadows & Highlights", "HighlightTonalWidth"); if (pedited) pedited->sh.htonalwidth = true; }
    if (keyFile.has_key ("Shadows & Highlights", "Shadows"))               { sh.shadows       = keyFile.get_integer ("Shadows & Highlights", "Shadows"); if (pedited) pedited->sh.shadows = true; }
    if (keyFile.has_key ("Shadows & Highlights", "ShadowTonalWidth"))      { sh.stonalwidth   = keyFile.get_integer ("Shadows & Highlights", "ShadowTonalWidth"); if (pedited) pedited->sh.stonalwidth = true; }
    if (keyFile.has_key ("Shadows & Highlights", "LocalContrast"))         { sh.localcontrast = keyFile.get_integer ("Shadows & Highlights", "LocalContrast"); if (pedited) pedited->sh.localcontrast = true; }
    if (keyFile.has_key ("Shadows & Highlights", "Radius"))                { sh.radius        = keyFile.get_integer ("Shadows & Highlights", "Radius"); if (pedited) pedited->sh.radius = true; }
}
    
    // load crop
if (keyFile.has_group ("Crop")) {
    if (keyFile.has_key ("Crop", "Enabled"))    { crop.enabled    = keyFile.get_boolean ("Crop", "Enabled"); if (pedited) pedited->crop.enabled = true; }
    if (keyFile.has_key ("Crop", "X"))          { crop.x          = keyFile.get_integer ("Crop", "X"); if (pedited) pedited->crop.x = true; }
    if (keyFile.has_key ("Crop", "Y"))          { crop.y          = keyFile.get_integer ("Crop", "Y"); if (pedited) pedited->crop.y = true; }
    if (keyFile.has_key ("Crop", "W"))          { crop.w          = keyFile.get_integer ("Crop", "W"); if (pedited) pedited->crop.w = true; }
    if (keyFile.has_key ("Crop", "H"))          { crop.h          = keyFile.get_integer ("Crop", "H"); if (pedited) pedited->crop.h = true; }
    if (keyFile.has_key ("Crop", "FixedRatio")) { crop.fixratio   = keyFile.get_boolean ("Crop", "FixedRatio"); if (pedited) pedited->crop.fixratio = true; }
    if (keyFile.has_key ("Crop", "Ratio")) {
        crop.ratio      = keyFile.get_string  ("Crop", "Ratio");
        if (pedited) pedited->crop.ratio = true;
        //backwards compatibility for crop.ratio
        if (crop.ratio=="DIN")    crop.ratio = "1.414 - DIN EN ISO 216";
        if (crop.ratio=="8.5:11") crop.ratio = "8.5:11 - US Letter";
        if (crop.ratio=="11:17")  crop.ratio = "11:17 - Tabloid";
    }
    if (keyFile.has_key ("Crop", "Orientation"))  { crop.orientation= keyFile.get_string  ("Crop", "Orientation"); if (pedited) pedited->crop.orientation = true; }
    if (keyFile.has_key ("Crop", "Guide"))        { crop.guide      = keyFile.get_string  ("Crop", "Guide"); if (pedited) pedited->crop.guide = true; }
}

    // load coarse
if (keyFile.has_group ("Coarse Transformation")) {
    if (keyFile.has_key ("Coarse Transformation", "Rotate"))          { coarse.rotate = keyFile.get_integer ("Coarse Transformation", "Rotate"); if (pedited) pedited->coarse.rotate = true; }
    if (keyFile.has_key ("Coarse Transformation", "HorizontalFlip"))  { coarse.hflip  = keyFile.get_boolean ("Coarse Transformation", "HorizontalFlip"); if (pedited) pedited->coarse.hflip = true; }
    if (keyFile.has_key ("Coarse Transformation", "VerticalFlip"))    { coarse.vflip  = keyFile.get_boolean ("Coarse Transformation", "VerticalFlip"); if (pedited) pedited->coarse.vflip = true; }
}

    // load rotate
if (keyFile.has_group ("Rotation")) {
    if (keyFile.has_key ("Rotation", "Degree"))   { rotate.degree = keyFile.get_double ("Rotation", "Degree"); if (pedited) pedited->rotate.degree = true; }
}
    // load commonTrans
if (keyFile.has_group ("Common Properties for Transformations")) {
    if (keyFile.has_key ("Common Properties for Transformations", "AutoFill"))   { commonTrans.autofill = keyFile.get_boolean ("Common Properties for Transformations", "AutoFill"); if (pedited) pedited->commonTrans.autofill = true; }
}

    // load distortion
if (keyFile.has_group ("Distortion")) {
    if (keyFile.has_key ("Distortion", "Amount"))     { distortion.amount     = keyFile.get_double  ("Distortion", "Amount"); if (pedited) pedited->distortion.amount = true; }
}

    // lens profile
if (keyFile.has_group ("LensProfile")) {
    if (keyFile.has_key ("LensProfile", "LCPFile")) { lensProf.lcpFile = keyFile.get_string ("LensProfile", "LCPFile"); if (pedited) pedited->lensProf.lcpFile = true; }
    if (keyFile.has_key ("LensProfile", "UseDistortion")) { lensProf.useDist = keyFile.get_boolean ("LensProfile", "UseDistortion"); if (pedited) pedited->lensProf.useDist = true; }
    if (keyFile.has_key ("LensProfile", "UseVignette")) { lensProf.useVign = keyFile.get_boolean ("LensProfile", "UseVignette"); if (pedited) pedited->lensProf.useVign = true; }
    if (keyFile.has_key ("LensProfile", "UseCA")) { lensProf.useCA = keyFile.get_boolean ("LensProfile", "UseCA"); if (pedited) pedited->lensProf.useCA = true; }
}
    
    // load perspective correction
if (keyFile.has_group ("Perspective")) {
    if (keyFile.has_key ("Perspective", "Horizontal"))  { perspective.horizontal = keyFile.get_integer ("Perspective", "Horizontal"); if (pedited) pedited->perspective.horizontal = true; }
    if (keyFile.has_key ("Perspective", "Vertical"))    { perspective.vertical   = keyFile.get_integer ("Perspective", "Vertical"); if (pedited) pedited->perspective.vertical = true; }
}

// load c/a correction
if (keyFile.has_group ("CACorrection")) {
    if (keyFile.has_key ("CACorrection", "Red"))  { cacorrection.red  = keyFile.get_double ("CACorrection", "Red"); if (pedited) pedited->cacorrection.red = true; }
    if (keyFile.has_key ("CACorrection", "Blue")) { cacorrection.blue = keyFile.get_double ("CACorrection", "Blue"); if (pedited) pedited->cacorrection.blue = true; }
}

    // load vignetting correction
if (keyFile.has_group ("Vignetting Correction")) {
    if (keyFile.has_key ("Vignetting Correction", "Amount"))   { vignetting.amount = keyFile.get_integer ("Vignetting Correction", "Amount"); if (pedited) pedited->vignetting.amount = true; }
    if (keyFile.has_key ("Vignetting Correction", "Radius"))   { vignetting.radius = keyFile.get_integer ("Vignetting Correction", "Radius"); if (pedited) pedited->vignetting.radius = true; }
    if (keyFile.has_key ("Vignetting Correction", "Strength")) { vignetting.strength = keyFile.get_integer ("Vignetting Correction", "Strength"); if (pedited) pedited->vignetting.strength = true; }
    if (keyFile.has_key ("Vignetting Correction", "CenterX"))  { vignetting.centerX = keyFile.get_integer ("Vignetting Correction", "CenterX"); if (pedited) pedited->vignetting.centerX = true; }
    if (keyFile.has_key ("Vignetting Correction", "CenterY"))  { vignetting.centerY = keyFile.get_integer ("Vignetting Correction", "CenterY"); if (pedited) pedited->vignetting.centerY = true; }
}

    // load highlight recovery settings
if (keyFile.has_group ("HLRecovery")) {
    if (keyFile.has_key ("HLRecovery", "Enabled"))  { hlrecovery.enabled  = keyFile.get_boolean ("HLRecovery", "Enabled"); if (pedited) pedited->hlrecovery.enabled = true; }
    if (keyFile.has_key ("HLRecovery", "Method"))   { hlrecovery.method   = keyFile.get_string  ("HLRecovery", "Method"); if (pedited) pedited->hlrecovery.method = true; }
}
    // load resize settings
if (keyFile.has_group ("Resize")) {
    if (keyFile.has_key ("Resize", "Enabled"))       { resize.enabled   = keyFile.get_boolean ("Resize", "Enabled"); if (pedited) pedited->resize.enabled = true; }
    if (keyFile.has_key ("Resize", "Scale"))         { resize.scale     = keyFile.get_double ("Resize", "Scale"); if (pedited) pedited->resize.scale = true; }
    if (keyFile.has_key ("Resize", "AppliesTo"))     { resize.appliesTo = keyFile.get_string ("Resize", "AppliesTo"); if (pedited) pedited->resize.appliesTo = true; }
    if (keyFile.has_key ("Resize", "Method"))        { resize.method    = keyFile.get_string ("Resize", "Method"); if (pedited) pedited->resize.method = true; }
    if (keyFile.has_key ("Resize", "DataSpecified")) { resize.dataspec  = keyFile.get_integer ("Resize", "DataSpecified"); if (pedited) pedited->resize.dataspec = true; }
    if (keyFile.has_key ("Resize", "Width"))         { resize.width     = keyFile.get_integer ("Resize", "Width"); if (pedited) pedited->resize.width = true; }
    if (keyFile.has_key ("Resize", "Height"))        { resize.height    = keyFile.get_integer ("Resize", "Height"); if (pedited) pedited->resize.height = true; }
}

    // load color management settings
if (keyFile.has_group ("Color Management")) {
    if (keyFile.has_key ("Color Management", "InputProfile"))     { icm.input            = keyFile.get_string ("Color Management", "InputProfile"); if (pedited) pedited->icm.input = true; }
    if (keyFile.has_key ("Color Management", "BlendCMSMatrix"))   { icm.blendCMSMatrix   = keyFile.get_boolean ("Color Management", "BlendCMSMatrix"); if (pedited) pedited->icm.blendCMSMatrix = true; }
    if (keyFile.has_key ("Color Management", "PreferredProfile")) { icm.preferredProfile = keyFile.get_boolean ("Color Management", "PreferredProfile"); if (pedited) pedited->icm.preferredProfile = true; }
    if (keyFile.has_key ("Color Management", "WorkingProfile"))   { icm.working          = keyFile.get_string ("Color Management", "WorkingProfile"); if (pedited) pedited->icm.working = true; }
    if (keyFile.has_key ("Color Management", "OutputProfile"))    { icm.output           = keyFile.get_string ("Color Management", "OutputProfile"); if (pedited) pedited->icm.output = true; }
    if (keyFile.has_key ("Color Management", "Gammafree"))        { icm.gamma            = keyFile.get_string ("Color Management", "Gammafree"); if (pedited) pedited->icm.gamfree = true; }
    if (keyFile.has_key ("Color Management", "Freegamma"))        { icm.freegamma        = keyFile.get_boolean ("Color Management", "Freegamma"); if (pedited) pedited->icm.freegamma = true; }
    if (keyFile.has_key ("Color Management", "GammaVal"))         { icm.gampos           = keyFile.get_double ("Color Management", "GammaVal"); if (pedited) pedited->icm.gamma = true; }
    if (keyFile.has_key ("Color Management", "GammaSlope"))       { icm.slpos            = keyFile.get_double ("Color Management", "GammaSlope"); if (pedited) pedited->icm.slpos = true; }

}

    // load directional pyramid equalizer parameters
if (keyFile.has_group ("Directional Pyramid Equalizer")) {
    if (keyFile.has_key ("Directional Pyramid Equalizer", "Enabled")) { dirpyrequalizer.enabled = keyFile.get_boolean ("Directional Pyramid Equalizer", "Enabled"); if (pedited) pedited->dirpyrequalizer.enabled = true; }
    for(int i = 0; i < 5; i ++) {
        std::stringstream ss;
        ss << "Mult" << i;
        if(keyFile.has_key ("Directional Pyramid Equalizer", ss.str())) { dirpyrequalizer.mult[i] = keyFile.get_double ("Directional Pyramid Equalizer", ss.str()); if (pedited) pedited->dirpyrequalizer.mult[i] = true; }
    }
}

    // load HSV equalizer parameters
if (keyFile.has_group ("HSV Equalizer")) {
    if (ppVersion>=300) {
        if (keyFile.has_key ("HSV Equalizer", "HCurve")) { hsvequalizer.hcurve = keyFile.get_double_list ("HSV Equalizer", "HCurve"); if (pedited) pedited->hsvequalizer.hcurve = true; }
        if (keyFile.has_key ("HSV Equalizer", "SCurve")) { hsvequalizer.scurve = keyFile.get_double_list ("HSV Equalizer", "SCurve"); if (pedited) pedited->hsvequalizer.scurve = true; }
        if (keyFile.has_key ("HSV Equalizer", "VCurve")) { hsvequalizer.vcurve = keyFile.get_double_list ("HSV Equalizer", "VCurve"); if (pedited) pedited->hsvequalizer.vcurve = true; }
    }
}

    // load RGB curves
if (keyFile.has_group ("RGB Curves")) {
    if (keyFile.has_key ("RGB Curves", "rCurve")) { rgbCurves.rcurve = keyFile.get_double_list ("RGB Curves", "rCurve"); if (pedited) pedited->rgbCurves.rcurve = true; }
    if (keyFile.has_key ("RGB Curves", "gCurve")) { rgbCurves.gcurve = keyFile.get_double_list ("RGB Curves", "gCurve"); if (pedited) pedited->rgbCurves.gcurve = true; }
    if (keyFile.has_key ("RGB Curves", "bCurve")) { rgbCurves.bcurve  = keyFile.get_double_list ("RGB Curves", "bCurve"); if (pedited) pedited->rgbCurves.bcurve = true; }
}

    // load raw settings
if (keyFile.has_group ("RAW")) {
    if (keyFile.has_key ("RAW", "DarkFrame"))        { raw.dark_frame = keyFile.get_string  ("RAW", "DarkFrame" ); if (pedited) pedited->raw.darkFrame = true; }
    if (keyFile.has_key ("RAW", "DarkFrameAuto"))    { raw.df_autoselect = keyFile.get_boolean ("RAW", "DarkFrameAuto" ); if (pedited) pedited->raw.dfAuto = true; }
    if (keyFile.has_key ("RAW", "FlatFieldFile"))       { raw.ff_file = keyFile.get_string  ("RAW", "FlatFieldFile" ); if (pedited) pedited->raw.ff_file = true; }
    if (keyFile.has_key ("RAW", "FlatFieldAutoSelect")) { raw.ff_AutoSelect = keyFile.get_boolean  ("RAW", "FlatFieldAutoSelect" );  if (pedited) pedited->raw.ff_AutoSelect = true; }
    if (keyFile.has_key ("RAW", "FlatFieldBlurRadius")) { raw.ff_BlurRadius = keyFile.get_integer  ("RAW", "FlatFieldBlurRadius" ); if (pedited) pedited->raw.ff_BlurRadius = true; }
    if (keyFile.has_key ("RAW", "FlatFieldBlurType"))   { raw.ff_BlurType = keyFile.get_string  ("RAW", "FlatFieldBlurType" ); if (pedited) pedited->raw.ff_BlurType = true; }
    if (keyFile.has_key ("RAW", "CA"))               { raw.ca_autocorrect = keyFile.get_boolean ("RAW", "CA" ); if (pedited) pedited->raw.caCorrection = true; }
    if (keyFile.has_key ("RAW", "CARed"))            { raw.cared = keyFile.get_double ("RAW", "CARed" ); if (pedited) pedited->raw.caRed = true; }
    if (keyFile.has_key ("RAW", "CABlue"))           { raw.cablue = keyFile.get_double ("RAW", "CABlue" ); if (pedited) pedited->raw.caBlue = true; }
    if (keyFile.has_key ("RAW", "HotDeadPixels"))    { raw.hotdeadpix_filt = keyFile.get_boolean ("RAW", "HotDeadPixels" ); if (pedited) pedited->raw.hotDeadPixelFilter = true; }
    if (keyFile.has_key ("RAW", "HotDeadPixelThresh")) { raw.hotdeadpix_thresh = keyFile.get_integer ("RAW", "HotDeadPixelThresh" ); if (pedited) pedited->raw.hotDeadPixelThresh = true; }
    if (keyFile.has_key ("RAW", "LineDenoise"))      { raw.linenoise = keyFile.get_integer ("RAW", "LineDenoise" ); if (pedited) pedited->raw.linenoise = true; }
    if (keyFile.has_key ("RAW", "GreenEqThreshold")) { raw.greenthresh= keyFile.get_integer ("RAW", "GreenEqThreshold"); if (pedited) pedited->raw.greenEq = true; }
    if (keyFile.has_key ("RAW", "CcSteps"))          { raw.ccSteps  = keyFile.get_integer ("RAW", "CcSteps"); if (pedited) pedited->raw.ccSteps = true; }
    if (keyFile.has_key ("RAW", "Method"))           { raw.dmethod = keyFile.get_string ("RAW", "Method"); if (pedited) pedited->raw.dmethod = true; }
    if (keyFile.has_key ("RAW", "DCBIterations"))    { raw.dcb_iterations = keyFile.get_integer("RAW", "DCBIterations"); if (pedited) pedited->raw.dcbIterations = true; }
    if (keyFile.has_key ("RAW", "DCBEnhance"))       { raw.dcb_enhance =keyFile.get_boolean("RAW", "DCBEnhance"); if (pedited) pedited->raw.dcbEnhance = true; }
    if (keyFile.has_key ("RAW", "ALLEnhance"))       { raw.all_enhance =keyFile.get_boolean("RAW", "ALLEnhance"); if (pedited) pedited->raw.allEnhance = true; }

    if (keyFile.has_key ("RAW", "PreExposure"))   { raw.expos =keyFile.get_double("RAW", "PreExposure"); if (pedited) pedited->raw.exPos = true; }
    if (keyFile.has_key ("RAW", "PrePreserv"))    { raw.preser =keyFile.get_double("RAW", "PrePreserv"); if (pedited) pedited->raw.exPreser = true; }
    if (keyFile.has_key ("RAW", "PreBlackzero"))  { raw.blackzero =keyFile.get_double("RAW", "PreBlackzero"); if (pedited) pedited->raw.exBlackzero = true; }
    if (keyFile.has_key ("RAW", "PreBlackone"))   { raw.blackone =keyFile.get_double("RAW", "PreBlackone"); if (pedited) pedited->raw.exBlackone = true; }
    if (keyFile.has_key ("RAW", "PreBlacktwo"))   { raw.blacktwo =keyFile.get_double("RAW", "PreBlacktwo"); if (pedited) pedited->raw.exBlacktwo = true; }
    if (keyFile.has_key ("RAW", "PreBlackthree")) { raw.blackthree =keyFile.get_double("RAW", "PreBlackthree"); if (pedited) pedited->raw.exBlackthree = true; }
    if (keyFile.has_key ("RAW", "PreTwoGreen"))   { raw.twogreen =keyFile.get_boolean("RAW", "PreTwoGreen"); if (pedited) pedited->raw.exTwoGreen = true; }

}

    // load exif change settings
/*if (keyFile.has_group ("Exif")) {
    std::vector<Glib::ustring> keys = keyFile.get_keys ("Exif");
    for (int i=0; i<(int)keys.size(); i++) {
        Glib::ustring tmpStr = keyFile.get_string ("Exif", keys[i]);
        exif[keys[i]] = keyFile.get_string ("Exif", keys[i]);
        if (pedited) pedited->exif = true;
    }
}*/

    /*
     * Load iptc change settings
     *
     * Existing values are preserved, and the stored values
     * are added to the list. To reset a field, the user has to
     * save the profile with the field leaved empty, but still
     * terminated by a semi-column ";"
     *
     * Please note that the old Keywords and SupplementalCategories
     * tag content is fully replaced by the new one,
     * i.e. they don't merge
     */
/*if (keyFile.has_group ("IPTC")) {
    std::vector<Glib::ustring> keys = keyFile.get_keys ("IPTC");
    IPTCPairs::iterator element;
    for (unsigned int i=0; i<keys.size(); i++) {
        // does this key already exist?
        element = iptc.find(keys[i]);
        if (element != iptc.end()) {
            // it already exist so we cleanup the values
            element->second.clear();
        }

        // TODO: look out if merging Keywords and SupplementalCategories from the procparams chain would be interesting
        std::vector<Glib::ustring> currIptc = keyFile.get_string_list ("IPTC", keys[i]);
        for (
            std::vector<Glib::ustring>::iterator currLoadedTagValue=currIptc.begin();
            currLoadedTagValue!=currIptc.end();
            currLoadedTagValue++)
        {
            iptc[keys[i]].push_back(currLoadedTagValue->data());
        }
        if (pedited) pedited->iptc = true;
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

/*bool operator==(const ExifPairs& a, const ExifPairs& b) {

    return a.field == b.field && a.value == b.value;
}

bool operator==(const IPTCPairs& a, const IPTCPairs& b) {

    return a.field == b.field && a.values == b.values;
}*/
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
		//&& colorBoost.amount == other.colorBoost.amount
		//&& colorBoost.avoidclip == other.colorBoost.avoidclip
		//&& colorBoost.enable_saturationlimiter == other.colorBoost.enable_saturationlimiter
		//&& colorBoost.saturationlimit == other.colorBoost.saturationlimit
		&& wb.method == other.wb.method
		&& wb.green == other.wb.green
		&& wb.temperature == other.wb.temperature
		//&& colorShift.a == other.colorShift.a
		//&& colorShift.b == other.colorShift.b
		&& impulseDenoise.enabled == other.impulseDenoise.enabled
		&& impulseDenoise.thresh == other.impulseDenoise.thresh
		&& dirpyrDenoise.enabled == other.dirpyrDenoise.enabled
		&& dirpyrDenoise.luma == other.dirpyrDenoise.luma
		&& dirpyrDenoise.chroma == other.dirpyrDenoise.chroma
		&& dirpyrDenoise.gamma == other.dirpyrDenoise.gamma
		&& edgePreservingDecompositionUI.enabled == other.edgePreservingDecompositionUI.enabled
		&& edgePreservingDecompositionUI.Strength == other.edgePreservingDecompositionUI.Strength
		&& edgePreservingDecompositionUI.EdgeStopping == other.edgePreservingDecompositionUI.EdgeStopping
		&& edgePreservingDecompositionUI.Scale == other.edgePreservingDecompositionUI.Scale
		&& edgePreservingDecompositionUI.ReweightingIterates == other.edgePreservingDecompositionUI.ReweightingIterates
		&& defringe.enabled == other.defringe.enabled
		&& defringe.radius == other.defringe.radius
		&& defringe.threshold == other.defringe.threshold
		//&& lumaDenoise.enabled == other.lumaDenoise.enabled
		//&& lumaDenoise.radius == other.lumaDenoise.radius
		//&& lumaDenoise.edgetolerance == other.lumaDenoise.edgetolerance
		//&& colorDenoise.enabled == other.colorDenoise.enabled
		//&& colorDenoise.edgetolerance == other.colorDenoise.edgetolerance
		//&& colorDenoise.edgesensitive == other.colorDenoise.edgesensitive
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
		&& distortion.amount == other.distortion.amount
        && lensProf.lcpFile == other.lensProf.lcpFile
        && lensProf.useDist == other.lensProf.useDist
        && lensProf.useVign == other.lensProf.useVign
        && lensProf.useCA == other.lensProf.useCA
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
		&& raw.hotdeadpix_thresh == other.raw.hotdeadpix_thresh
		&& raw.dmethod == other.raw.dmethod
		&& raw.greenthresh == other.raw.greenthresh
		&& raw.linenoise == other.raw.linenoise
		&& icm.input == other.icm.input
		&& icm.blendCMSMatrix == other.icm.blendCMSMatrix
        && icm.preferredProfile == other.icm.preferredProfile
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
		&& rgbCurves.rcurve == other.rgbCurves.rcurve
		&& rgbCurves.gcurve == other.rgbCurves.gcurve
		&& rgbCurves.bcurve == other.rgbCurves.bcurve
		&& exif==other.exif
		&& iptc==other.iptc
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

PartialProfile::PartialProfile(bool createInstance) {
    if (createInstance) {
        pparams = new ProcParams();
        pedited = new ParamsEdited();
    }
    else {
        pparams = NULL;
        pedited=NULL;
    }
}

PartialProfile::PartialProfile(ProcParams* pp, ParamsEdited* pe, bool fullCopy) {
    if (fullCopy && pp) {
        pparams = new ProcParams(*pp);
    }
    else
        pparams = pp;

    if (fullCopy && pe) {
        pedited = new ParamsEdited(*pe);
    }
    else
        pedited = pe;
}

PartialProfile::PartialProfile(const ProcParams* pp, const ParamsEdited* pe) {
    if (pp) {
        pparams = new ProcParams(*pp);
    }
    else
        pparams = NULL;

    if (pe) {
        pedited = new ParamsEdited(*pe);
    }
    else
        pedited = NULL;
}

int PartialProfile::load (Glib::ustring fName) {
    if (!pparams) pparams = new ProcParams();
    if (!pedited) pedited = new ParamsEdited();
    return pparams->load(fName, pedited);
}

void PartialProfile::deleteInstance () {
    if (pparams) { delete pparams; pparams = NULL; }
    if (pedited) { delete pedited; pedited = NULL; }
}

/*
 * Set the all values of the General section to false
 * in order to preserve them in applyTo
 */
void PartialProfile::clearGeneral () {
    if (pedited) {
        pedited->general.colorlabel = false;
        pedited->general.intrash = false;
        pedited->general.rank = false;
    }
}

void PartialProfile::applyTo(ProcParams *destParams) const {
    if (destParams && pparams && pedited) {
        pedited->combine(*destParams, *pparams, true);
    }
}

void PartialProfile::set(bool v) {
    if (pedited) pedited->set(v);
}

}
}

