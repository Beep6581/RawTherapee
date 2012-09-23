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
    toneCurve.curve2.clear ();
    toneCurve.curve2.push_back(DCT_Linear);
    toneCurve.curveMode     = ToneCurveParams::TC_MODE_STD;
    toneCurve.curveMode2    = ToneCurveParams::TC_MODE_STD;

    labCurve.brightness      = 0;
    labCurve.contrast        = 0;
    labCurve.chromaticity    = 0;
    labCurve.avoidcolorshift = true;
    labCurve.lcredsk         = true;

    labCurve.rstprotection   = 0;
    labCurve.bwtoning        = false;
    labCurve.lcurve.clear ();
    labCurve.lcurve.push_back(DCT_Linear);
    labCurve.acurve.clear ();
    labCurve.acurve.push_back(DCT_Linear);
    labCurve.bcurve.clear ();
    labCurve.bcurve.push_back(DCT_Linear);
    labCurve.cccurve.clear ();
    labCurve.cccurve.push_back(DCT_Linear);
    labCurve.chcurve.clear ();
    labCurve.chcurve.push_back(FCT_Linear);
    labCurve.lccurve.clear ();
    labCurve.lccurve.push_back(DCT_Linear);

    rgbCurves.rcurve.clear ();
    rgbCurves.rcurve.push_back(DCT_Linear);
    rgbCurves.gcurve.clear ();
    rgbCurves.gcurve.push_back(DCT_Linear);
    rgbCurves.bcurve.clear ();
    rgbCurves.bcurve.push_back(DCT_Linear);


    sharpenEdge.enabled         = false;
    sharpenEdge.passes          = 2;
    sharpenEdge.amount          = 50.0;
    sharpenEdge.threechannels   = false;

    sharpenMicro.enabled        = false;
    sharpenMicro.amount         = 20.0;
    sharpenMicro.uniformity     = 50.0;
    sharpenMicro.matrix         = false;

    sharpening.enabled          = false;
    sharpening.radius           = 1.0;
    sharpening.amount           = 90;
    sharpening.threshold.setValues(20, 80, 2000, 1200);
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
    vibrance.psthreshold.setValues(0, 75);
    vibrance.protectskins       = false;
    vibrance.avoidcolorshift    = true;
    vibrance.pastsattog     	= true;
    vibrance.skintonescurve.clear ();
    vibrance.skintonescurve.push_back(DCT_Linear);

    wb.method       = "Camera";
    wb.temperature  = 6504;
    wb.green        = 1.0;

    impulseDenoise.enabled      = false;
    impulseDenoise.thresh       = 50;

    defringe.enabled            = false;
    defringe.radius             = 2.0;
    defringe.threshold          = 25;

    dirpyrDenoise.enabled       = false;
    dirpyrDenoise.luma          = 30;
    dirpyrDenoise.Ldetail       = 50;
    dirpyrDenoise.chroma        = 30;
    dirpyrDenoise.gamma         = 1.7;
	dirpyrDenoise.expcomp       = 0.0;

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
    icm.toneCurve = false;
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
    //raw.all_enhance=false;

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
		&& toneCurve.curve2 == other.toneCurve.curve2
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
		&& toneCurve.curveMode == other.toneCurve.curveMode
		&& toneCurve.curveMode2 == other.toneCurve.curveMode2
		&& labCurve.lcurve == other.labCurve.lcurve
		&& labCurve.acurve == other.labCurve.acurve
		&& labCurve.bcurve == other.labCurve.bcurve
		&& labCurve.cccurve == other.labCurve.cccurve
		&& labCurve.chcurve == other.labCurve.chcurve
		&& labCurve.lccurve == other.labCurve.lccurve
		&& labCurve.brightness == other.labCurve.brightness
		&& labCurve.contrast == other.labCurve.contrast
		&& labCurve.chromaticity == other.labCurve.chromaticity
		&& labCurve.avoidcolorshift == other.labCurve.avoidcolorshift
		&& labCurve.rstprotection == other.labCurve.rstprotection
		&& labCurve.bwtoning == other.labCurve.bwtoning
		&& labCurve.lcredsk == other.labCurve.lcredsk
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
		&& vibrance.skintonescurve == other.vibrance.skintonescurve
		&& wb.method == other.wb.method
		&& wb.green == other.wb.green
		&& wb.temperature == other.wb.temperature
		&& impulseDenoise.enabled == other.impulseDenoise.enabled
		&& impulseDenoise.thresh == other.impulseDenoise.thresh
		&& dirpyrDenoise.enabled == other.dirpyrDenoise.enabled
		&& dirpyrDenoise.luma == other.dirpyrDenoise.luma
		&& dirpyrDenoise.Ldetail == other.dirpyrDenoise.Ldetail
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
		//&& raw.all_enhance == other.raw.all_enhance
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
		&& icm.toneCurve == other.icm.toneCurve
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
		//&& exif==other.exif
		//&& iptc==other.iptc
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

PartialProfile::PartialProfile(const PartialProfile* source) {
    pparams = NULL;
    pedited = NULL;

    if (source) {
        if (source->pparams)
            pparams = new ProcParams(*source->pparams);
        if (source->pedited)
            pedited = new ParamsEdited(*source->pedited);
    }
}

PartialProfile::PartialProfile(bool createPParamsInstance, bool createPEditedInstance) {
    if (createPParamsInstance) pparams = new ProcParams();
    else pparams = NULL;
    if (createPEditedInstance) pedited = new ParamsEdited();
    else pedited = NULL;
}

/*
 * Copy a PartialProfile class
 * If fullCopy=true, the instance of pparams and pedited are duplicated too if non NULL
 */
PartialProfile::PartialProfile(bool fullCopy, ProcParams* pp, ParamsEdited* pe) {
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

/*
 * Assignment operator: Make a full copy of the PartialProfile
 */
PartialProfile& PartialProfile::operator=(const PartialProfile& rhs) {
	if (rhs.pparams) {
		if (!pparams)
			pparams = new ProcParams();
		*pparams = *rhs.pparams;
	}
	else if (pparams) {
		delete pparams;
		pparams = NULL;
	}

	if (rhs.pedited) {
		if (!pedited)
			pedited = new ParamsEdited();
		*pedited = *rhs.pedited;
	}
	else if (pedited) {
		delete pedited;
		pedited = NULL;
	}
	return *this;
};

void PartialProfile::saveIntoXMP(Exiv2::XmpData &xmpData, const std::string& baseKey) const
{
    std::string prefix;

    if (!pparams)
        return;

    xmpData[baseKey+"rt:"+kXmpVersion] =                int(PPVERSION);

    prefix=baseKey+"rt:ExposureRGB/";
    if (!pedited || pedited->toneCurve.autoexp)    xmpData[prefix+"rt:Auto"]=                       pparams->toneCurve.autoexp;
    if (!pedited || pedited->toneCurve.clip)       xmpData[prefix+"rt:Clip"]=                       pparams->toneCurve.clip;
    if (!pedited || pedited->toneCurve.expcomp)    xmpData[prefix+"rt:Compensation"]=               pparams->toneCurve.expcomp;
    if (!pedited || pedited->toneCurve.brightness) xmpData[prefix+"rt:Brightness"]=                 pparams->toneCurve.brightness;
    if (!pedited || pedited->toneCurve.contrast)   xmpData[prefix+"rt:Contrast"]=                   pparams->toneCurve.contrast;
    if (!pedited || pedited->toneCurve.saturation) xmpData[prefix+"rt:Saturation"]=                 pparams->toneCurve.saturation;
    if (!pedited || pedited->toneCurve.black)      xmpData[prefix+"rt:Black"]=                      pparams->toneCurve.black;
    if (!pedited || pedited->toneCurve.hlcompr)    xmpData[prefix+"rt:HighlightCompression"]=       pparams->toneCurve.hlcompr;
    if (!pedited || pedited->toneCurve.hlcomprthresh) xmpData[prefix+"rt:HighlightComprThreshold"]= pparams->toneCurve.hlcomprthresh;
    if (!pedited || pedited->toneCurve.shcompr)    xmpData[prefix+"rt:ShadowCompression"]=          pparams->toneCurve.shcompr;
    if (!pedited || pedited->toneCurve.curve)      xmpData[prefix+"rt:ToneCurve"]=                  serializeVector(pparams->toneCurve.curve);
    if (!pedited || pedited->toneCurve.curve2)     xmpData[prefix+"rt:ToneCurve2"]=                 serializeVector(pparams->toneCurve.curve2);
    if (!pedited || pedited->toneCurve.curveMode)  {
        Glib::ustring method;
        switch (pparams->toneCurve.curveMode) {
        case (ToneCurveParams::TC_MODE_STD):
            method = "Standard";
            break;
        case (ToneCurveParams::TC_MODE_FILMLIKE):
            method = "FilmLike";
            break;
        case (ToneCurveParams::TC_MODE_SATANDVALBLENDING):
            method = "SatAndValueBlending";
            break;
        case (ToneCurveParams::TC_MODE_WEIGHTEDSTD):
            method = "WeightedStd";
            break;
        }
        xmpData[prefix+"rt:ToneCurveMode"]= method;
    }
    if (!pedited || pedited->toneCurve.curveMode2)  {
        Glib::ustring method;
        switch (pparams->toneCurve.curveMode2) {
        case (ToneCurveParams::TC_MODE_STD):
            method = "Standard";
            break;
        case (ToneCurveParams::TC_MODE_FILMLIKE):
            method = "FilmLike";
            break;
        case (ToneCurveParams::TC_MODE_SATANDVALBLENDING):
            method = "SatAndValueBlending";
            break;
        case (ToneCurveParams::TC_MODE_WEIGHTEDSTD):
            method = "WeightedStd";
            break;
        }
        xmpData[prefix+"rt:ToneCurveMode2"]= method;
    }

    prefix=baseKey+"rt:ChannelMixer/";
    if (!pedited || pedited->chmixer.red[0] || pedited->chmixer.red[1] || pedited->chmixer.red[2])
        xmpData[prefix+"rt:Red"]   = serializeArray(pparams->chmixer.red,3);
    if (!pedited || pedited->chmixer.green[0] || pedited->chmixer.green[1] || pedited->chmixer.green[2])
        xmpData[prefix+"rt:Green"] = serializeArray(pparams->chmixer.green,3);
    if (!pedited || pedited->chmixer.blue[0] || pedited->chmixer.blue[1] || pedited->chmixer.blue[2])
        xmpData[prefix+"rt:Blue"]  = serializeArray(pparams->chmixer.blue,3);

    prefix=baseKey+"rt:ExposureLab/";
    if (!pedited || pedited->labCurve.brightness)      xmpData[prefix+"rt:Brightness"]=             pparams->labCurve.brightness;
    if (!pedited || pedited->labCurve.contrast)        xmpData[prefix+"rt:Contrast"]=               pparams->labCurve.contrast;
    if (!pedited || pedited->labCurve.chromaticity)    xmpData[prefix+"rt:Chromaticity"]=           pparams->labCurve.chromaticity;
    if (!pedited || pedited->labCurve.avoidcolorshift) xmpData[prefix+"rt:AvoidColorShift"]=        pparams->labCurve.avoidcolorshift;
    if (!pedited || pedited->labCurve.bwtoning)        xmpData[prefix+"rt:BWToning"]=               pparams->labCurve.bwtoning;
    if (!pedited || pedited->labCurve.rstprotection)   xmpData[prefix+"rt:RedSkinTonesProtection"]= pparams->labCurve.rstprotection;
    if (!pedited || pedited->labCurve.lcredsk)         xmpData[prefix+"rt:LcForRedsAndSkinsOnly"]=  pparams->labCurve.lcredsk;
    if (!pedited || pedited->labCurve.lcurve)          xmpData[prefix+"rt:LCurve"]=                 serializeVector(pparams->labCurve.lcurve);
    if (!pedited || pedited->labCurve.acurve)          xmpData[prefix+"rt:aCurve"]=                 serializeVector(pparams->labCurve.acurve);
    if (!pedited || pedited->labCurve.bcurve)          xmpData[prefix+"rt:bCurve"]=                 serializeVector(pparams->labCurve.bcurve);
    if (!pedited || pedited->labCurve.cccurve)         xmpData[prefix+"rt:ccCurve"]=                serializeVector(pparams->labCurve.cccurve);
    if (!pedited || pedited->labCurve.chcurve)         xmpData[prefix+"rt:chCurve"]=                serializeVector(pparams->labCurve.chcurve);
    if (!pedited || pedited->labCurve.lccurve)         xmpData[prefix+"rt:LcCurve"]=                serializeVector(pparams->labCurve.lccurve);

	prefix=baseKey+"rt:Vibrance/";
    if (!pedited || pedited->vibrance.enabled)          xmpData[prefix+"rt:Enabled"]=         pparams->vibrance.enabled;
    if (!pedited || pedited->vibrance.pastels)          xmpData[prefix+"rt:Pastels"]=         pparams->vibrance.pastels;
    if (!pedited || pedited->vibrance.saturated)        xmpData[prefix+"rt:Saturated"]=       pparams->vibrance.saturated;
    if (!pedited || pedited->vibrance.psthreshold)      xmpData[prefix+"rt:PSThreshold"]=     serializeArray(pparams->vibrance.psthreshold.value, 2);
    if (!pedited || pedited->vibrance.protectskins)     xmpData[prefix+"rt:ProtectSkins"]=    pparams->vibrance.protectskins;
    if (!pedited || pedited->vibrance.avoidcolorshift)  xmpData[prefix+"rt:AvoidColorShift"]= pparams->vibrance.avoidcolorshift;
    if (!pedited || pedited->vibrance.pastsattog)       xmpData[prefix+"rt:PastSatTog"]=      pparams->vibrance.pastsattog;
    if (!pedited || pedited->vibrance.skintonescurve)   xmpData[prefix+"rt:SkinTonesCurve"]=  serializeVector(pparams->vibrance.skintonescurve);

	prefix=baseKey+"rt:Sharpening/";
    if (!pedited || pedited->sharpening.enabled)            xmpData[prefix+"rt:Enabled"] =            pparams->sharpening.enabled;
    if (!pedited || pedited->sharpening.method)             xmpData[prefix+"rt:Method"] =             pparams->sharpening.method;
    if (!pedited || pedited->sharpening.radius)             xmpData[prefix+"rt:Radius"]=              pparams->sharpening.radius;
    if (!pedited || pedited->sharpening.amount)             xmpData[prefix+"rt:Amount"]=              pparams->sharpening.amount;
    if (!pedited || pedited->sharpening.threshold)          xmpData[prefix+"rt:Threshold"]=           serializeArray(pparams->sharpening.threshold.value, 4);
    if (!pedited || pedited->sharpening.edgesonly)          xmpData[prefix+"rt:OnlyEdges"]=           pparams->sharpening.edgesonly;
    if (!pedited || pedited->sharpening.edges_radius)       xmpData[prefix+"rt:EdgeDetectionRadius"]= pparams->sharpening.edges_radius;
    if (!pedited || pedited->sharpening.edges_tolerance)    xmpData[prefix+"rt:EdgeTolerance"]=       pparams->sharpening.edges_tolerance;
    if (!pedited || pedited->sharpening.halocontrol)        xmpData[prefix+"rt:HaloControlEnabled"]=  pparams->sharpening.halocontrol;
    if (!pedited || pedited->sharpening.halocontrol_amount) xmpData[prefix+"rt:HaloControlAmount"]=   pparams->sharpening.halocontrol_amount;
    if (!pedited || pedited->sharpening.deconvradius)       xmpData[prefix+"rt:DeconvRadius"]=        pparams->sharpening.deconvradius;
    if (!pedited || pedited->sharpening.deconvamount)       xmpData[prefix+"rt:DeconvAmount"]=        pparams->sharpening.deconvamount;
    if (!pedited || pedited->sharpening.deconvdamping)      xmpData[prefix+"rt:DeconvDamping"]=       pparams->sharpening.deconvdamping;
    if (!pedited || pedited->sharpening.deconviter)         xmpData[prefix+"rt:DeconvIterations"]=    pparams->sharpening.deconviter;

	prefix=baseKey+"rt:SharpenEdge/";
    if (!pedited || pedited->sharpenEdge.enabled)       xmpData[prefix+"rt:Enabled"]=       pparams->sharpenEdge.enabled;
    if (!pedited || pedited->sharpenEdge.passes)        xmpData[prefix+"rt:Passes"]=        pparams->sharpenEdge.passes;
    if (!pedited || pedited->sharpenEdge.amount)        xmpData[prefix+"rt:ThreeChannels"]= pparams->sharpenEdge.threechannels;
    if (!pedited || pedited->sharpenEdge.threechannels) xmpData[prefix+"rt:Amount"]=        pparams->sharpenEdge.amount;

	prefix=baseKey+"rt:MicroContrast/";
    if (!pedited || pedited->sharpenMicro.enabled)      xmpData[prefix+"rt:Enabled"]=    pparams->sharpenMicro.enabled;
    if (!pedited || pedited->sharpenMicro.matrix)       xmpData[prefix+"rt:Uniformity"]= pparams->sharpenMicro.uniformity;
    if (!pedited || pedited->sharpenMicro.amount)       xmpData[prefix+"rt:Matrix"]=     pparams->sharpenMicro.matrix;
    if (!pedited || pedited->sharpenMicro.uniformity)   xmpData[prefix+"rt:Amount"]=     pparams->sharpenMicro.amount;

	prefix=baseKey+"rt:WhiteBalance/";
    if (!pedited || pedited->wb.method)      xmpData[prefix+"rt:Mode"]=        pparams->wb.method;
    if (!pedited || pedited->wb.temperature) xmpData[prefix+"rt:Temperature"]= pparams->wb.temperature;
    if (!pedited || pedited->wb.green)       xmpData[prefix+"rt:Green"]=       pparams->wb.green;

    prefix=baseKey+"rt:ImpulseDenoise/";
    if (!pedited || pedited->impulseDenoise.enabled) xmpData[prefix+"rt:Enabled"]=   pparams->impulseDenoise.enabled;
    if (!pedited || pedited->impulseDenoise.thresh)  xmpData[prefix+"rt:Threshold"]= pparams->impulseDenoise.thresh;

    prefix=baseKey+"rt:Defringe/";
    if (!pedited || pedited->defringe.enabled)       xmpData[prefix+"rt:Enabled"]=   pparams->defringe.enabled;
    if (!pedited || pedited->defringe.radius)        xmpData[prefix+"rt:Radius"]=    pparams->defringe.radius;
    if (!pedited || pedited->defringe.threshold)     xmpData[prefix+"rt:Threshold"]= pparams->defringe.threshold;

    prefix=baseKey+"rt:PyramidDenoise/";
    if (!pedited || pedited->dirpyrDenoise.enabled) xmpData[prefix+"rt:Enabled"]= pparams->dirpyrDenoise.enabled;
    if (!pedited || pedited->dirpyrDenoise.luma)    xmpData[prefix+"rt:Luma"]=    pparams->dirpyrDenoise.luma;
    if (!pedited || pedited->dirpyrDenoise.Ldetail) xmpData[prefix+"rt:LDetail"]= pparams->dirpyrDenoise.Ldetail;
    if (!pedited || pedited->dirpyrDenoise.chroma)  xmpData[prefix+"rt:Chroma"]=  pparams->dirpyrDenoise.chroma;
    if (!pedited || pedited->dirpyrDenoise.gamma)   xmpData[prefix+"rt:Gamma"]=   pparams->dirpyrDenoise.gamma;

    prefix=baseKey+"rt:ShadowHighlights/";
    if (!pedited || pedited->sh.enabled)       xmpData[prefix+"rt:Enabled"]=             pparams->sh.enabled;
    if (!pedited || pedited->sh.hq)            xmpData[prefix+"rt:HighQuality"]=         pparams->sh.hq;
    if (!pedited || pedited->sh.highlights)    xmpData[prefix+"rt:Highlights"]=          pparams->sh.highlights;
    if (!pedited || pedited->sh.htonalwidth)   xmpData[prefix+"rt:HighlightTonalWidth"]= pparams->sh.htonalwidth;
    if (!pedited || pedited->sh.shadows)       xmpData[prefix+"rt:Shadows"]=             pparams->sh.shadows;
    if (!pedited || pedited->sh.stonalwidth)   xmpData[prefix+"rt:ShadowTonalWidth"]=    pparams->sh.stonalwidth;
    if (!pedited || pedited->sh.localcontrast) xmpData[prefix+"rt:LocalContrast"]=       pparams->sh.localcontrast;
    if (!pedited || pedited->sh.radius)        xmpData[prefix+"rt:Radius"]=              pparams->sh.radius;

    prefix=baseKey+"rt:Crop/";
    if (!pedited || pedited->crop.enabled)     xmpData[prefix+"rt:Enabled"]= pparams->crop.enabled;
    if (!pedited || pedited->crop.x)           xmpData[prefix+"rt:X"]=       pparams->crop.x;
    if (!pedited || pedited->crop.y)           xmpData[prefix+"rt:Y"]=       pparams->crop.y;
    if (!pedited || pedited->crop.w)           xmpData[prefix+"rt:Width"]=   pparams->crop.w;
    if (!pedited || pedited->crop.h)           xmpData[prefix+"rt:Height"]=  pparams->crop.h;
    /*
    if (!pedited || pedited->crop.fixratio)    xmpData[prefix+"rt:FixedRatio"]=  pparams->crop.fixratio;
    if (!pedited || pedited->crop.ratio)       xmpData[prefix+"rt:Ratio"]=       pparams->crop.ratio;
    if (!pedited || pedited->crop.orientation) xmpData[prefix+"rt:Orientation"]= pparams->crop.orientation;
    if (!pedited || pedited->crop.guide)       xmpData[prefix+"rt:Guide"]=       pparams->crop.guide;
    */

    prefix=baseKey+"rt:CoarseGeo/";
    if (!pedited || pedited->coarse.rotate)    xmpData[prefix+"rt:RotationDegree"]= pparams->coarse.rotate;
    if (!pedited || pedited->coarse.hflip)     xmpData[prefix+"rt:HorizontalFlip"]= pparams->coarse.hflip;
    if (!pedited || pedited->coarse.vflip)     xmpData[prefix+"rt:VerticalFlip"]=   pparams->coarse.vflip;

    prefix=baseKey+"rt:EPD/";
    if (!pedited || pedited->edgePreservingDecompositionUI.enabled)             xmpData[prefix+"rt:Enabled"]=              pparams->edgePreservingDecompositionUI.enabled;
    if (!pedited || pedited->edgePreservingDecompositionUI.Strength)            xmpData[prefix+"rt:Strength"]=             pparams->edgePreservingDecompositionUI.Strength;
    if (!pedited || pedited->edgePreservingDecompositionUI.EdgeStopping)        xmpData[prefix+"rt:EdgeStopping"]=         pparams->edgePreservingDecompositionUI.EdgeStopping;
    if (!pedited || pedited->edgePreservingDecompositionUI.Scale)               xmpData[prefix+"rt:Scale"] =               pparams->edgePreservingDecompositionUI.Scale;
    if (!pedited || pedited->edgePreservingDecompositionUI.ReweightingIterates) xmpData[prefix+"rt:ReweightingIterates"] = pparams->edgePreservingDecompositionUI.ReweightingIterates;


    prefix=baseKey+"rt:Geometry/";
    //xmpData[prefix+"rt:Enabled"]=    geo.enabled;
    if (!pedited || pedited->commonTrans.autofill)   xmpData[prefix+"rt:AutoFill"]=              pparams->commonTrans.autofill;
    if (!pedited || pedited->rotate.degree)          xmpData[prefix+"rt:RotationDegree"]=        pparams->rotate.degree;
    if (!pedited || pedited->distortion.amount)      xmpData[prefix+"rt:DistortionAmount"]=      pparams->distortion.amount;
    if (!pedited || pedited->perspective.horizontal) xmpData[prefix+"rt:HorizontalPerspective"]= pparams->perspective.horizontal;
    if (!pedited || pedited->perspective.vertical)   xmpData[prefix+"rt:VerticalPerspective"]=   pparams->perspective.vertical;

    prefix=baseKey+"rt:LensProfile/";
    if (!pedited || pedited->lensProf.lcpFile)       xmpData[prefix+"rt:LCPFile"]=       pparams->lensProf.lcpFile;
    if (!pedited || pedited->lensProf.useDist)       xmpData[prefix+"rt:UseDistortion"]= pparams->lensProf.useDist;
    if (!pedited || pedited->lensProf.useVign)       xmpData[prefix+"rt:UseVignette"]=   pparams->lensProf.useVign;
    if (!pedited || pedited->lensProf.useCA)         xmpData[prefix+"rt:UseCA"]=         pparams->lensProf.useCA;

    prefix=baseKey+"rt:CACorrection/";
    //xmpData[prefix+"rt:Enabled"]=    cacorrection.enabled;
    if (!pedited || pedited->cacorrection.red)       xmpData[prefix+"rt:Red"]=  pparams->cacorrection.red;
    if (!pedited || pedited->cacorrection.blue)      xmpData[prefix+"rt:Blue"]= pparams->cacorrection.blue;

    prefix=baseKey+"rt:Vignetting/";
    //xmpData[prefix+"rt:Enabled"]= vignetting.enabled;
    if (!pedited || pedited->vignetting.amount)      xmpData[prefix+"rt:Amount"] =  pparams->vignetting.amount;
    if (!pedited || pedited->vignetting.radius)      xmpData[prefix+"rt:Radius"] =  pparams->vignetting.radius;
    if (!pedited || pedited->vignetting.strength)    xmpData[prefix+"rt:Strength"]= pparams->vignetting.strength;
    if (!pedited || pedited->vignetting.centerX)     xmpData[prefix+"rt:CenterX"] = pparams->vignetting.centerX;
    if (!pedited || pedited->vignetting.centerY)     xmpData[prefix+"rt:CenterY"] = pparams->vignetting.centerY;

    prefix=baseKey+"rt:HLRecovery/";
    if (!pedited || pedited->hlrecovery.enabled)     xmpData[prefix+"rt:Enabled"]=  pparams->hlrecovery.enabled;
    if (!pedited || pedited->hlrecovery.method)      xmpData[prefix+"rt:Method"]=   pparams->hlrecovery.method;

    prefix=baseKey+"rt:Resize/";
    if (!pedited || pedited->resize.enabled)         xmpData[prefix+"rt:Enabled"]=       pparams->resize.enabled;
    if (!pedited || pedited->resize.scale)           xmpData[prefix+"rt:Scale"]  =       pparams->resize.scale;
    if (!pedited || pedited->resize.appliesTo)       xmpData[prefix+"rt:AppliesTo"]=     pparams->resize.appliesTo;
    if (!pedited || pedited->resize.method)          xmpData[prefix+"rt:Method"]=        pparams->resize.method;
    if (!pedited || pedited->resize.dataspec)        xmpData[prefix+"rt:DataSpecified"]= pparams->resize.dataspec;
    if (!pedited || pedited->resize.width)           xmpData[prefix+"rt:Width"] =        pparams->resize.width;
    if (!pedited || pedited->resize.height)          xmpData[prefix+"rt:Height"] =       pparams->resize.height;

    prefix=baseKey+"rt:ColorManagement/";
    // save color management settings
    if (!pedited || pedited->icm.input)              xmpData[prefix+"rt:InputProfile"] =     pparams->icm.input;
    if (!pedited || pedited->icm.toneCurve)          xmpData[prefix+"rt:UseToneCurve"] =     pparams->icm.toneCurve;
    if (!pedited || pedited->icm.blendCMSMatrix)     xmpData[prefix+"rt:BlendCMSMatrix"] =   pparams->icm.blendCMSMatrix;
    if (!pedited || pedited->icm.preferredProfile)   xmpData[prefix+"rt:PreferredProfile"] = pparams->icm.preferredProfile;
    if (!pedited || pedited->icm.working)            xmpData[prefix+"rt:WorkingProfile"] =   pparams->icm.working;
    if (!pedited || pedited->icm.output)             xmpData[prefix+"rt:OutputProfile"] =    pparams->icm.output;
    if (!pedited || pedited->icm.gamma)              xmpData[prefix+"rt:GammaProfile"] =     pparams->icm.gamma;
    if (!pedited || pedited->icm.freegamma)          xmpData[prefix+"rt:FreeGammaEnabled"] = pparams->icm.freegamma;
    if (!pedited || pedited->icm.gampos)             xmpData[prefix+"rt:FreeGammaValue"] =   pparams->icm.gampos;
    if (!pedited || pedited->icm.slpos)              xmpData[prefix+"rt:FreeGammaSlope"] =   pparams->icm.slpos;

    prefix=baseKey+"rt:DirectionalPyramidEqualizer/";
    if (!pedited || pedited->dirpyrequalizer.enabled) xmpData[prefix+"rt:Enabled"] = pparams->dirpyrequalizer.enabled ;
    if (!pedited || pedited->dirpyrequalizer.mult[0] || pedited->dirpyrequalizer.mult[1] || pedited->dirpyrequalizer.mult[2]
      || pedited->dirpyrequalizer.mult[3] || pedited->dirpyrequalizer.mult[4])
    	xmpData[prefix+"rt:Coeff"] =    serializeArray( pparams->dirpyrequalizer.mult,5);

    prefix=baseKey+"rt:HSVEqualizer/";
    //xmpData[prefix+"rt:Enabled"]= hsvequalizer.enabled;
    if (!pedited || pedited->hsvequalizer.hcurve)    xmpData[prefix+"rt:HCurve"] = serializeVector(pparams->hsvequalizer.hcurve);
    if (!pedited || pedited->hsvequalizer.scurve)    xmpData[prefix+"rt:SCurve"] = serializeVector(pparams->hsvequalizer.scurve);
    if (!pedited || pedited->hsvequalizer.vcurve)    xmpData[prefix+"rt:VCurve"] = serializeVector(pparams->hsvequalizer.vcurve);

    prefix=baseKey+"rt:RGBCurves/";
    //xmpData[prefix+"rt:Enabled"]= hsvequalizer.enabled;
    if (!pedited || pedited->rgbCurves.rcurve)       xmpData[prefix+"rt:RCurve"] = serializeVector(pparams->rgbCurves.rcurve);
    if (!pedited || pedited->rgbCurves.gcurve)       xmpData[prefix+"rt:GCurve"] = serializeVector(pparams->rgbCurves.gcurve);
    if (!pedited || pedited->rgbCurves.bcurve)       xmpData[prefix+"rt:BCurve"] = serializeVector(pparams->rgbCurves.bcurve);

    prefix=baseKey+"rt:RawArithmetic/";
    //xmpData[prefix+"rt:Enabled"]=
    if (!pedited || pedited->raw.darkFrame)          xmpData[prefix+"rt:DarkFrameFile"]=       pparams->raw.dark_frame;
    if (!pedited || pedited->raw.dfAuto)             xmpData[prefix+"rt:DarkFrameAutoSelect"]= pparams->raw.df_autoselect ;
    if (!pedited || pedited->raw.ff_file)            xmpData[prefix+"rt:FlatFieldFile"]=       pparams->raw.ff_file ;
    if (!pedited || pedited->raw.ff_AutoSelect)      xmpData[prefix+"rt:FlatFieldAutoSelect"]= pparams->raw.ff_AutoSelect ;
    if (!pedited || pedited->raw.ff_BlurRadius)      xmpData[prefix+"rt:FlatFieldBlurRadius"]= pparams->raw.ff_BlurRadius ;
    if (!pedited || pedited->raw.ff_BlurType)        xmpData[prefix+"rt:FlatFieldBlurType"]=   pparams->raw.ff_BlurType ;

    prefix=baseKey+"rt:RawCACorrection/";
    //xmpData[prefix+"rt:Enabled"]=
    if (!pedited || pedited->raw.caCorrection)       xmpData[prefix+"rt:Auto"]= pparams->raw.ca_autocorrect;
    if (!pedited || pedited->raw.caRed)              xmpData[prefix+"rt:Red"] = pparams->raw.cared;
    if (!pedited || pedited->raw.caBlue)             xmpData[prefix+"rt:Blue"]= pparams->raw.cablue;

    prefix=baseKey+"rt:HotDeadPixelCorrection/";
    if (!pedited || pedited->raw.hotDeadPixelFilter) xmpData[prefix+"rt:Enabled"]=   pparams->raw.hotdeadpix_filt;
    if (!pedited || pedited->raw.hotDeadPixelThresh) xmpData[prefix+"rt:Threshold"]= pparams->raw.hotdeadpix_thresh;

    prefix=baseKey+"rt:RawDenoise/";
    //xmpData[prefix+"rt:Enabled"]=
    if (!pedited || pedited->raw.linenoise)          xmpData[prefix+"rt:LineDenoise"]= pparams->raw.linenoise;

    prefix=baseKey+"rt:Demosaicing/";
    if (!pedited || pedited->raw.greenEq)            xmpData[prefix+"rt:GreenEqThreshold"]= pparams->raw.greenthresh;
    if (!pedited || pedited->raw.ccSteps)            xmpData[prefix+"rt:CcSteps"]=          pparams->raw.ccSteps;
    if (!pedited || pedited->raw.dmethod)            xmpData[prefix+"rt:Method"]=           pparams->raw.dmethod;
    if (!pedited || pedited->raw.dcbIterations)      xmpData[prefix+"rt:DCBIterations"]=    pparams->raw.dcb_iterations;
    if (!pedited || pedited->raw.dcbEnhance)         xmpData[prefix+"rt:DCBEnhance"]=       pparams->raw.dcb_enhance;
    //if (!pedited || pedited->raw.allEnhance)         xmpData[prefix+"rt:AllEnhance"]=       pparams->raw.all_enhance;

    prefix=baseKey+"rt:RawExposure/";
    if (!pedited || pedited->raw.exPos)              xmpData[prefix+"rt:Exposure"]=     pparams->raw.expos;
    if (!pedited || pedited->raw.exPreser)           xmpData[prefix+"rt:HLPreserving"]= pparams->raw.preser;
    if (!pedited || pedited->raw.exBlackzero)        xmpData[prefix+"rt:Black0"]=       pparams->raw.blackzero;
    if (!pedited || pedited->raw.exBlackone)         xmpData[prefix+"rt:Black1"]=       pparams->raw.blackone;
    if (!pedited || pedited->raw.exBlacktwo)         xmpData[prefix+"rt:Black2"]=       pparams->raw.blacktwo;
    if (!pedited || pedited->raw.exBlackthree)       xmpData[prefix+"rt:Black3"]=       pparams->raw.blackthree;
    if (!pedited || pedited->raw.exTwoGreen)         xmpData[prefix+"rt:TwoGreen"]=     pparams->raw.twogreen;
}

int PartialProfile::loadFromXMP(Exiv2::XmpData &xmpData, const std::string& baseKey)
{
    if (!pparams)
        return 2;

    try {
        std::string prefix;
        if(! readVarFromXmp( xmpData, baseKey+"rt:"+kXmpVersion, pparams->ppVersion) )
            return 2;

        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:ExposureRGB")) != xmpData.end()){
            prefix = baseKey+"rt:ExposureRGB/rt:";
            if (readVarFromXmp( xmpData, prefix+"Auto", pparams->toneCurve.autoexp ) && pedited) pedited->toneCurve.autoexp = true;
            if (readVarFromXmp( xmpData, prefix+"Clip", pparams->toneCurve.clip ) && pedited) pedited->toneCurve.clip = true;
            if (readVarFromXmp( xmpData, prefix+"Compensation", pparams->toneCurve.expcomp ) && pedited) pedited->toneCurve.expcomp = true;
            if (readVarFromXmp( xmpData, prefix+"Brightness", pparams->toneCurve.brightness ) && pedited) pedited->toneCurve.brightness = true;
            if (readVarFromXmp( xmpData, prefix+"Contrast", pparams->toneCurve.contrast ) && pedited) pedited->toneCurve.contrast = true;
            if (readVarFromXmp( xmpData, prefix+"Saturation", pparams->toneCurve.saturation ) && pedited) pedited->toneCurve.saturation = true;
            if (readVarFromXmp( xmpData, prefix+"Black", pparams->toneCurve.black ) && pedited) pedited->toneCurve.black = true;
            if (readVarFromXmp( xmpData, prefix+"HighlightCompression", pparams->toneCurve.hlcompr ) && pedited) pedited->toneCurve.hlcompr = true;
            if (readVarFromXmp( xmpData, prefix+"HighlightComprThreshold", pparams->toneCurve.hlcomprthresh ) && pedited) pedited->toneCurve.hlcomprthresh = true;
            if (readVarFromXmp( xmpData, prefix+"ShadowCompression", pparams->toneCurve.shcompr ) && pedited) pedited->toneCurve.shcompr = true;
            if (readVarFromXmp( xmpData, prefix+"ToneCurve", pparams->toneCurve.curve) && pedited) pedited->toneCurve.curve = true;
            int x = pparams->toneCurve.curveMode;
            if (readVarFromXmp( xmpData, prefix+"ToneCurveMode", x) && pedited) pedited->toneCurve.curveMode = true;
            if (readVarFromXmp( xmpData, prefix+"ToneCurve2", pparams->toneCurve.curve2) && pedited) pedited->toneCurve.curve2 = true;
            x = pparams->toneCurve.curveMode2;
            if (readVarFromXmp( xmpData, prefix+"ToneCurveMode2", x) && pedited) pedited->toneCurve.curveMode2 = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:ChannelMixer")) != xmpData.end()){
            prefix=baseKey+"rt:ChannelMixer/rt:";
            if (readVarFromXmp( xmpData, prefix+"Red", pparams->chmixer.red,3 ) && pedited) pedited->chmixer.red[0] = pedited->chmixer.red[1] = pedited->chmixer.red[2] = true;
            if (readVarFromXmp( xmpData, prefix+"Green", pparams->chmixer.green,3 ) && pedited) pedited->chmixer.green[0] = pedited->chmixer.green[1] = pedited->chmixer.green[2] = true;
            if (readVarFromXmp( xmpData, prefix+"Blue", pparams->chmixer.blue,3 ) && pedited) pedited->chmixer.blue[0] = pedited->chmixer.blue[1] = pedited->chmixer.blue[2] = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:ExposureLab")) != xmpData.end()){
            prefix=baseKey+"rt:ExposureLab/rt:";
            if (readVarFromXmp( xmpData, prefix+"Brightness", pparams->labCurve.brightness ) && pedited) pedited->labCurve.brightness = true;
            if (readVarFromXmp( xmpData, prefix+"Contrast", pparams->labCurve.contrast ) && pedited) pedited->labCurve.contrast = true;
            if (readVarFromXmp( xmpData, prefix+"Chromaticity", pparams->labCurve.chromaticity ) && pedited) pedited->labCurve.chromaticity = true;
            if (readVarFromXmp( xmpData, prefix+"AvoidColorShift", pparams->labCurve.avoidcolorshift ) && pedited) pedited->labCurve.avoidcolorshift = true;
            if (readVarFromXmp( xmpData, prefix+"BWToning", pparams->labCurve.bwtoning ) && pedited) pedited->labCurve.bwtoning = true;
            if (readVarFromXmp( xmpData, prefix+"RedSkinTonesProtection", pparams->labCurve.rstprotection ) && pedited) pedited->labCurve.rstprotection = true;
            if (readVarFromXmp( xmpData, prefix+"LcForRedsAndSkinsOnly", pparams->labCurve.lcredsk ) && pedited) pedited->labCurve.lcredsk = true;
            if (readVarFromXmp( xmpData, prefix+"LCurve", pparams->labCurve.lcurve ) && pedited) pedited->labCurve.lcurve = true;
            if (readVarFromXmp( xmpData, prefix+"aCurve", pparams->labCurve.acurve ) && pedited) pedited->labCurve.acurve = true;
            if (readVarFromXmp( xmpData, prefix+"bCurve", pparams->labCurve.bcurve ) && pedited) pedited->labCurve.bcurve = true;
            if (readVarFromXmp( xmpData, prefix+"ccCurve", pparams->labCurve.cccurve ) && pedited) pedited->labCurve.cccurve = true;
            if (readVarFromXmp( xmpData, prefix+"chCurve", pparams->labCurve.chcurve ) && pedited) pedited->labCurve.chcurve = true;
            if (readVarFromXmp( xmpData, prefix+"LcCurve", pparams->labCurve.lccurve ) && pedited) pedited->labCurve.lccurve = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:Sharpening")) != xmpData.end()){
            prefix=baseKey+"rt:Sharpening/rt:";
            if (readVarFromXmp( xmpData, prefix+kXmpEnabled, pparams->sharpening.enabled ) && pedited) pedited->sharpening.enabled = true;
            if (readVarFromXmp( xmpData, prefix+"Method", pparams->sharpening.method ) && pedited) pedited->sharpening.method = true;
            if (readVarFromXmp( xmpData, prefix+"Radius", pparams->sharpening.radius ) && pedited) pedited->sharpening.radius = true;
            if (readVarFromXmp( xmpData, prefix+"Amount", pparams->sharpening.amount ) && pedited) pedited->sharpening.amount = true;
            if (readVarFromXmp( xmpData, prefix+"Threshold", pparams->sharpening.threshold.value, 4 ) && pedited) pedited->sharpening.threshold = true;
            if (readVarFromXmp( xmpData, prefix+"OnlyEdges", pparams->sharpening.edgesonly ) && pedited) pedited->sharpening.edgesonly = true;
            if (readVarFromXmp( xmpData, prefix+"EdgeDetectionRadius", pparams->sharpening.edges_radius ) && pedited) pedited->sharpening.edges_radius = true;
            if (readVarFromXmp( xmpData, prefix+"EdgeTolerance", pparams->sharpening.edges_tolerance ) && pedited) pedited->sharpening.edges_tolerance = true;
            if (readVarFromXmp( xmpData, prefix+"HaloControlEnabled", pparams->sharpening.halocontrol ) && pedited) pedited->sharpening.halocontrol = true;
            if (readVarFromXmp( xmpData, prefix+"HaloControlAmount", pparams->sharpening.halocontrol_amount ) && pedited) pedited->sharpening.halocontrol_amount = true;
            if (readVarFromXmp( xmpData, prefix+"DeconvRadius", pparams->sharpening.deconvradius ) && pedited) pedited->sharpening.deconvradius = true;
            if (readVarFromXmp( xmpData, prefix+"DeconvAmount", pparams->sharpening.deconvamount ) && pedited) pedited->sharpening.deconvamount = true;
            if (readVarFromXmp( xmpData, prefix+"DeconvDamping", pparams->sharpening.deconvdamping ) && pedited) pedited->sharpening.deconvdamping = true;
            if (readVarFromXmp( xmpData, prefix+"DeconvIterations", pparams->sharpening.deconviter ) && pedited) pedited->sharpening.deconviter = true;
        }

        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:Vibrance")) != xmpData.end()){
            prefix=baseKey+"rt:Vibrance/rt:";
            if (readVarFromXmp( xmpData, prefix+kXmpEnabled, pparams->vibrance.enabled ) && pedited) pedited->vibrance.enabled = true;
            if (readVarFromXmp( xmpData, prefix+"Pastels",  pparams->vibrance.pastels ) && pedited) pedited->vibrance.pastels = true;
            if (readVarFromXmp( xmpData, prefix+"Saturated", pparams->vibrance.saturated ) && pedited) pedited->vibrance.saturated = true;
            if (readVarFromXmp( xmpData, prefix+"PSThreshold", pparams->vibrance.psthreshold.value, 2 ) && pedited) pedited->vibrance.psthreshold = true;
            if (readVarFromXmp( xmpData, prefix+"ProtectSkins", pparams->vibrance.protectskins ) && pedited) pedited->vibrance.protectskins = true;
            if (readVarFromXmp( xmpData, prefix+"AvoidColorShift", pparams->vibrance.avoidcolorshift ) && pedited) pedited->vibrance.avoidcolorshift = true;
            if (readVarFromXmp( xmpData, prefix+"PastSatTog", pparams->vibrance.pastsattog ) && pedited) pedited->vibrance.pastsattog = true;
            if (readVarFromXmp( xmpData, prefix+"SkinTonesCurve", pparams->vibrance.skintonescurve ) && pedited) pedited->vibrance.skintonescurve = true;
        }

        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:SharpenEdge")) != xmpData.end()){
            prefix=baseKey+"rt:SharpenEdge/rt:";
            if (readVarFromXmp( xmpData, prefix+kXmpEnabled, pparams->sharpenEdge.enabled ) && pedited) pedited->sharpenEdge.enabled = true;
            if (readVarFromXmp( xmpData, prefix+"Passes", pparams->sharpenEdge.passes ) && pedited) pedited->sharpenEdge.passes = true;
            if (readVarFromXmp( xmpData, prefix+"Amount", pparams->sharpenEdge.amount ) && pedited) pedited->sharpenEdge.amount = true;
            if (readVarFromXmp( xmpData, prefix+"ThreeChannels", pparams->sharpenEdge.threechannels ) && pedited) pedited->sharpenEdge.threechannels = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:MicroContrast")) != xmpData.end()){
            prefix=baseKey+"rt:MicroContrast/rt:";
            if (readVarFromXmp( xmpData, prefix+kXmpEnabled, pparams->sharpenMicro.enabled ) && pedited) pedited->sharpenMicro.enabled = true;
            if (readVarFromXmp( xmpData, prefix+"Amount", pparams->sharpenMicro.amount) && pedited) pedited->sharpenMicro.amount = true;
            if (readVarFromXmp( xmpData, prefix+"Uniformity", pparams->sharpenMicro.uniformity ) && pedited) pedited->sharpenMicro.uniformity = true;
            if (readVarFromXmp( xmpData, prefix+"Matrix", pparams->sharpenMicro.matrix ) && pedited) pedited->sharpenMicro.matrix = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:WhiteBalance")) != xmpData.end()){
            prefix=baseKey+"rt:WhiteBalance/rt:";
            if (readVarFromXmp( xmpData, prefix+"Mode", pparams->wb.method ) && pedited) pedited->wb.method = true;
            if (readVarFromXmp( xmpData, prefix+"Temperature", pparams->wb.temperature ) && pedited) pedited->wb.temperature = true;
            if (readVarFromXmp( xmpData, prefix+"Green", pparams->wb.green ) && pedited) pedited->wb.green = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:ImpulseDenoise")) != xmpData.end()){
            prefix=baseKey+"rt:ImpulseDenoise/rt:";
            if (readVarFromXmp( xmpData, prefix+kXmpEnabled, pparams->impulseDenoise.enabled ) && pedited) pedited->impulseDenoise.enabled = true;
            if (readVarFromXmp( xmpData, prefix+"Threshold", pparams->impulseDenoise.thresh ) && pedited) pedited->impulseDenoise.thresh = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:Defringe")) != xmpData.end()){
            prefix=baseKey+"rt:Defringe/rt:";
            if (readVarFromXmp( xmpData, prefix+kXmpEnabled, pparams->defringe.enabled ) && pedited) pedited->defringe.enabled = true;
            if (readVarFromXmp( xmpData, prefix+"Radius", pparams->defringe.radius ) && pedited) pedited->defringe.radius = true;
            if (readVarFromXmp( xmpData, prefix+"Threshold", pparams->defringe.threshold ) && pedited) pedited->defringe.threshold = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:PyramidDenoise")) != xmpData.end()){
            prefix=baseKey+"rt:PyramidDenoise/rt:";
            if (readVarFromXmp( xmpData, prefix+kXmpEnabled, pparams->dirpyrDenoise.enabled ) && pedited) pedited->dirpyrDenoise.enabled = true;
            if (readVarFromXmp( xmpData, prefix+"Luma", pparams->dirpyrDenoise.luma ) && pedited) pedited->dirpyrDenoise.luma = true;
            if (readVarFromXmp( xmpData, prefix+"LDetail", pparams->dirpyrDenoise.Ldetail ) && pedited) pedited->dirpyrDenoise.Ldetail = true;
            if (readVarFromXmp( xmpData, prefix+"Chroma", pparams->dirpyrDenoise.chroma ) && pedited) pedited->dirpyrDenoise.chroma = true;
            if (readVarFromXmp( xmpData, prefix+"Gamma", pparams->dirpyrDenoise.gamma ) && pedited) pedited->dirpyrDenoise.gamma = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:ShadowHighlights")) != xmpData.end()){
            prefix=baseKey+"rt:ShadowHighlights/rt:";
            if (readVarFromXmp( xmpData, prefix+kXmpEnabled, pparams->sh.enabled ) && pedited) pedited->sh.enabled = true;
            if (readVarFromXmp( xmpData, prefix+"HighQuality", pparams->sh.hq ) && pedited) pedited->sh.hq = true;
            if (readVarFromXmp( xmpData, prefix+"Highlights", pparams->sh.highlights ) && pedited) pedited->sh.highlights = true;
            if (readVarFromXmp( xmpData, prefix+"HighlightTonalWidth", pparams->sh.htonalwidth ) && pedited) pedited->sh.htonalwidth = true;
            if (readVarFromXmp( xmpData, prefix+"Shadows", pparams->sh.shadows ) && pedited) pedited->sh.shadows = true;
            if (readVarFromXmp( xmpData, prefix+"ShadowTonalWidth", pparams->sh.stonalwidth ) && pedited) pedited->sh.stonalwidth = true;
            if (readVarFromXmp( xmpData, prefix+"LocalContrast", pparams->sh.localcontrast ) && pedited) pedited->sh.localcontrast = true;
            if (readVarFromXmp( xmpData, prefix+"Radius", pparams->sh.radius ) && pedited) pedited->sh.radius = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:Crop")) != xmpData.end()){
            prefix=baseKey+"rt:Crop/rt:";
            if (readVarFromXmp( xmpData, prefix+kXmpEnabled, pparams->crop.enabled ) && pedited) pedited->crop.enabled = true;
            if (readVarFromXmp( xmpData, prefix+"X", pparams->crop.x ) && pedited) pedited->crop.x = true;
            if (readVarFromXmp( xmpData, prefix+"Y", pparams->crop.y ) && pedited) pedited->crop.y = true;
            if (readVarFromXmp( xmpData, prefix+"Width", pparams->crop.w ) && pedited) pedited->crop.w = true;
            if (readVarFromXmp( xmpData, prefix+"Height", pparams->crop.h ) && pedited) pedited->crop.h = true;
         //   if (xmpData[prefix+"rt:FixedRatio"]=  pparams->crop.fixratio;
         //   if (xmpData[prefix+"rt:Ratio"]=       pparams->crop.ratio;
         //   if (xmpData[prefix+"rt:Orientation"]= pparams->crop.orientation;
         //   if (xmpData[prefix+"rt:Guide"]=       pparams->crop.guide;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:CoarseGeo")) != xmpData.end()){
            prefix=baseKey+"rt:CoarseGeo/rt:";
            if (readVarFromXmp( xmpData, prefix+"RotationDegree", pparams->coarse.rotate ) && pedited) pedited->coarse.rotate = true;
            if (readVarFromXmp( xmpData, prefix+"HorizontalFlip", pparams->coarse.hflip ) && pedited) pedited->coarse.hflip = true;
            if (readVarFromXmp( xmpData, prefix+"VerticalFlip", pparams->coarse.vflip ) && pedited) pedited->coarse.vflip = true;
        }

        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:EPD")) != xmpData.end()){
            prefix=baseKey+"rt:EPD/rt:";
            if (readVarFromXmp( xmpData, prefix+kXmpEnabled, pparams->edgePreservingDecompositionUI.enabled ) && pedited) pedited->edgePreservingDecompositionUI.enabled = true;
            if (readVarFromXmp( xmpData, prefix+"Strength", pparams->edgePreservingDecompositionUI.Strength ) && pedited) pedited->edgePreservingDecompositionUI.Strength = true;
            if (readVarFromXmp( xmpData, prefix+"EdgeStopping", pparams->edgePreservingDecompositionUI.EdgeStopping) && pedited) pedited->edgePreservingDecompositionUI.EdgeStopping = true;
            if (readVarFromXmp( xmpData, prefix+"Scale", pparams->edgePreservingDecompositionUI.Scale ) && pedited) pedited->edgePreservingDecompositionUI.Scale = true;
            if (readVarFromXmp( xmpData, prefix+"ReweightingIterates", pparams->edgePreservingDecompositionUI.ReweightingIterates) && pedited) pedited->edgePreservingDecompositionUI.ReweightingIterates = true;
        }

        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:Geometry")) != xmpData.end()){
            prefix=baseKey+"rt:Geometry/rt:";
            //if (readVarFromXmp( xmpData, prefix+kXmpEnabled, geo.enabled);
            if (readVarFromXmp( xmpData, prefix+"AutoFill", pparams->commonTrans.autofill ) && pedited) pedited->commonTrans.autofill = true;
            if (readVarFromXmp( xmpData, prefix+"RotationDegree", pparams->rotate.degree ) && pedited) pedited->rotate.degree = true;
            if (readVarFromXmp( xmpData, prefix+"DistortionAmount", pparams->distortion.amount ) && pedited) pedited->distortion.amount = true;
            if (readVarFromXmp( xmpData, prefix+"HorizontalPerspective", pparams->perspective.horizontal ) && pedited) pedited->perspective.horizontal = true;
            if (readVarFromXmp( xmpData, prefix+"VerticalPerspective", pparams->perspective.vertical ) && pedited) pedited->perspective.vertical = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:LensProfile")) != xmpData.end()){
            prefix=baseKey+"rt:LensProfile/rt:";
            if (readVarFromXmp( xmpData, prefix+"LCPFile", pparams->lensProf.lcpFile ) && pedited) pedited->lensProf.lcpFile = true;
            if (readVarFromXmp( xmpData, prefix+"UseDistortion", pparams->lensProf.useDist ) && pedited) pedited->lensProf.useDist = true;
            if (readVarFromXmp( xmpData, prefix+"UseVignette", pparams->lensProf.useVign ) && pedited) pedited->lensProf.useVign = true;
            if (readVarFromXmp( xmpData, prefix+"UseCA", pparams->lensProf.useCA ) && pedited) pedited->lensProf.useCA = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:CACorrection")) != xmpData.end()){
            prefix=baseKey+"rt:CACorrection/rt:";
            //if (readVarFromXmp( xmpData, prefix+kXmpEnabled, pparams->cacorrection.enabled) && pedited) pedited-> = true;
            if (readVarFromXmp( xmpData, prefix+"Red", pparams->cacorrection.red ) && pedited) pedited->cacorrection.red = true;
            if (readVarFromXmp( xmpData, prefix+"Blue", pparams->cacorrection.blue ) && pedited) pedited->cacorrection.blue = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:Vignetting")) != xmpData.end()){
            prefix=baseKey+"rt:Vignetting/rt:";
            //if (readVarFromXmp( xmpData, prefix+kXmpEnabled, pparams->vignetting.enabled) && pedited) pedited-> = true;
            if (readVarFromXmp( xmpData, prefix+"Amount", pparams->vignetting.amount ) && pedited) pedited->vignetting.amount = true;
            if (readVarFromXmp( xmpData, prefix+"Radius", pparams->vignetting.radius ) && pedited) pedited->vignetting.radius = true;
            if (readVarFromXmp( xmpData, prefix+"Strength", pparams->vignetting.strength ) && pedited) pedited->vignetting.strength = true;
            if (readVarFromXmp( xmpData, prefix+"CenterX", pparams->vignetting.centerX ) && pedited) pedited->vignetting.centerX = true;
            if (readVarFromXmp( xmpData, prefix+"CenterY", pparams->vignetting.centerY ) && pedited) pedited->vignetting.centerY = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:HLRecovery")) != xmpData.end()){
            prefix=baseKey+"rt:HLRecovery/rt:";
            if (readVarFromXmp( xmpData, prefix+kXmpEnabled, pparams->hlrecovery.enabled ) && pedited) pedited->hlrecovery.enabled = true;
            if (readVarFromXmp( xmpData, prefix+"Method", pparams->hlrecovery.method ) && pedited) pedited->hlrecovery.method = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:Resize")) != xmpData.end()){
            prefix=baseKey+"rt:Resize/rt:";
            if (readVarFromXmp( xmpData, prefix+kXmpEnabled, pparams->resize.enabled ) && pedited) pedited->resize.enabled = true;
            if (readVarFromXmp( xmpData, prefix+"Scale", pparams->resize.scale ) && pedited) pedited->resize.scale = true;
            if (readVarFromXmp( xmpData, prefix+"AppliesTo", pparams->resize.appliesTo ) && pedited) pedited->resize.appliesTo = true;
            if (readVarFromXmp( xmpData, prefix+"Method", pparams->resize.method ) && pedited) pedited->resize.method = true;
            if (readVarFromXmp( xmpData, prefix+"DataSpecified", pparams->resize.dataspec ) && pedited) pedited->resize.dataspec = true;
            if (readVarFromXmp( xmpData, prefix+"Width", pparams->resize.width ) && pedited) pedited->resize.width = true;
            if (readVarFromXmp( xmpData, prefix+"Height", pparams->resize.height ) && pedited) pedited->resize.height = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:ColorManagement")) != xmpData.end()){
            prefix=baseKey+"rt:ColorManagement/rt:";
            if (readVarFromXmp( xmpData, prefix+"InputProfile", pparams->icm.input ) && pedited) pedited->icm.input = true;
            if (readVarFromXmp( xmpData, prefix+"UseToneCurve", pparams->icm.toneCurve ) && pedited) pedited->icm.toneCurve = true;
            if (readVarFromXmp( xmpData, prefix+"BlendCMSMatrix", pparams->icm.blendCMSMatrix ) && pedited) pedited->icm.blendCMSMatrix = true;
            if (readVarFromXmp( xmpData, prefix+"WorkingProfile", pparams->icm.working ) && pedited) pedited->icm.working = true;
            if (readVarFromXmp( xmpData, prefix+"OutputProfile", pparams->icm.output ) && pedited) pedited->icm.output = true;
            if (readVarFromXmp( xmpData, prefix+"GammaProfile", pparams->icm.gamma ) && pedited) pedited->icm.gamma = true;
            if (readVarFromXmp( xmpData, prefix+"FreeGammaEnabled", pparams->icm.freegamma ) && pedited) pedited->icm.freegamma = true;
            if (readVarFromXmp( xmpData, prefix+"FreeGammaValue", pparams->icm.gampos ) && pedited) pedited->icm.gampos = true;
            if (readVarFromXmp( xmpData, prefix+"FreeGammaSlope", pparams->icm.slpos ) && pedited) pedited->icm.slpos = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:DirectionalPyramidEqualizer")) != xmpData.end()){
            prefix=baseKey+"rt:DirectionalPyramidEqualizer/rt:";
            if (readVarFromXmp( xmpData, prefix+kXmpEnabled, pparams->dirpyrequalizer.enabled ) && pedited) pedited->dirpyrequalizer.enabled = true;
            if (readVarFromXmp( xmpData, prefix+"Coeff", pparams->dirpyrequalizer.mult,5 ) && pedited) for (int i=0; i<8; i++) pedited->dirpyrequalizer.mult[i] = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:HSVEqualizer")) != xmpData.end()){
            prefix=baseKey+"rt:HSVEqualizer/rt:";
            //if (readVarFromXmp( xmpData, prefix+kXmpEnabled, hsvequalizer.enabled) && pedited) pedited-> = true;
            if (readVarFromXmp( xmpData,prefix+"HCurve", pparams->hsvequalizer.hcurve) && pedited) pedited->hsvequalizer.hcurve = true;
            if (readVarFromXmp( xmpData,prefix+"SCurve", pparams->hsvequalizer.scurve) && pedited) pedited->hsvequalizer.scurve = true;
            if (readVarFromXmp( xmpData,prefix+"VCurve", pparams->hsvequalizer.vcurve) && pedited) pedited->hsvequalizer.vcurve = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:RGBCurves")) != xmpData.end()){
            prefix=baseKey+"rt:RGBCurves/rt:";
            if (readVarFromXmp( xmpData,prefix+"RCurve", pparams->rgbCurves.rcurve) && pedited) pedited->rgbCurves.rcurve = true;
            if (readVarFromXmp( xmpData,prefix+"GCurve", pparams->rgbCurves.gcurve) && pedited) pedited->rgbCurves.gcurve = true;
            if (readVarFromXmp( xmpData,prefix+"BCurve", pparams->rgbCurves.bcurve) && pedited) pedited->rgbCurves.bcurve = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:RawArithmetic")) != xmpData.end()){
            prefix=baseKey+"rt:RawArithmetic/rt:";
            //if (readVarFromXmp( xmpData, prefix+kXmpEnabled, raw.arithmeticEnabled) && pedited) pedited-> = true;
            if (readVarFromXmp( xmpData, prefix+"DarkFrameFile", pparams->raw.dark_frame ) && pedited) pedited->raw.darkFrame = true;
            if (readVarFromXmp( xmpData, prefix+"DarkFrameAutoSelect", pparams->raw.df_autoselect ) && pedited) pedited->raw.dfAuto = true;
            if (readVarFromXmp( xmpData, prefix+"FlatFieldFile", pparams->raw.ff_file ) && pedited) pedited->raw.ff_file = true;
            if (readVarFromXmp( xmpData, prefix+"FlatFieldAutoSelect", pparams->raw.ff_AutoSelect ) && pedited) pedited->raw.ff_AutoSelect = true;
            if (readVarFromXmp( xmpData, prefix+"FlatFieldBlurRadius", pparams->raw.ff_BlurRadius ) && pedited) pedited->raw.ff_BlurRadius = true;
            if (readVarFromXmp( xmpData, prefix+"FlatFieldBlurType", pparams->raw.ff_BlurType ) && pedited) pedited->raw.ff_BlurType = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:RawCACorrection")) != xmpData.end()){
            prefix=baseKey+"rt:RawCACorrection/rt:";
            //if (readVarFromXmp( xmpData, prefix+kXmpEnabled, raw.ca_enabled) && pedited) pedited-> = true;
            if (readVarFromXmp( xmpData, prefix+"Auto", pparams->raw.ca_autocorrect ) && pedited) pedited->raw.caCorrection = true;
            if (readVarFromXmp( xmpData, prefix+"Red", pparams->raw.cared ) && pedited) pedited->raw.caRed = true;
            if (readVarFromXmp( xmpData, prefix+"Blue", pparams->raw.cablue ) && pedited) pedited->raw.caBlue = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:HotDeadPixelCorrection")) != xmpData.end()){
            prefix=baseKey+"rt:HotDeadPixelCorrection/rt:";
            if (readVarFromXmp( xmpData, prefix+kXmpEnabled, pparams->raw.hotdeadpix_filt ) && pedited) pedited->raw.hotDeadPixelFilter = true;
            if (readVarFromXmp( xmpData, prefix+"Threshold", pparams->raw.hotdeadpix_thresh ) && pedited) pedited->raw.hotDeadPixelThresh = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:RawDenoise")) != xmpData.end()){
            prefix=baseKey+"rt:RawDenoise/rt:";
            //if (readVarFromXmp( xmpData, prefix+kXmpEnabled, raw.linenoiseEnabled ) && pedited) pedited-> = true;
            if (readVarFromXmp( xmpData, prefix+"LineDenoise", pparams->raw.linenoise ) && pedited) pedited->raw.linenoise = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:Demosaicing")) != xmpData.end()){
            prefix=baseKey+"rt:Demosaicing/rt:";
            if (readVarFromXmp( xmpData, prefix+"GreenEqThreshold", pparams->raw.greenthresh ) && pedited) pedited->raw.greenEq = true;
            if (readVarFromXmp( xmpData, prefix+"CcSteps", pparams->raw.ccSteps ) && pedited) pedited->raw.ccSteps = true;
            if (readVarFromXmp( xmpData, prefix+"Method", pparams->raw.dmethod ) && pedited) pedited->raw.dmethod = true;
            if (readVarFromXmp( xmpData, prefix+"DCBIterations", pparams->raw.dcb_iterations ) && pedited) pedited->raw.dcbIterations = true;
            if (readVarFromXmp( xmpData, prefix+"DCBEnhance", pparams->raw.dcb_enhance ) && pedited) pedited->raw.dcbEnhance = true;
            //if (readVarFromXmp( xmpData, prefix+"AllEnhance", pparams->raw.all_enhance ) && pedited) pedited->raw.allEnhance = true;
        }
        if( xmpData.findKey(Exiv2::XmpKey(baseKey+"rt:RawExposure")) != xmpData.end()){
            prefix=baseKey+"rt:RawExposure/rt:";
            if (readVarFromXmp( xmpData, prefix+"Exposure", pparams->raw.expos ) && pedited) pedited->raw.exPos = true;
            if (readVarFromXmp( xmpData, prefix+"HLPreserving", pparams->raw.preser ) && pedited) pedited->raw.exPreser = true;
            if (readVarFromXmp( xmpData, prefix+"Black0", pparams->raw.blackzero ) && pedited) pedited->raw.exBlackzero = true;
            if (readVarFromXmp( xmpData, prefix+"Black1", pparams->raw.blackone ) && pedited) pedited->raw.exBlackone = true;
            if (readVarFromXmp( xmpData, prefix+"Black2", pparams->raw.blacktwo ) && pedited) pedited->raw.exBlacktwo = true;
            if (readVarFromXmp( xmpData, prefix+"Black3", pparams->raw.blackthree ) && pedited) pedited->raw.exBlackthree = true;
            if (readVarFromXmp( xmpData, prefix+"TwoGreen", pparams->raw.twogreen ) && pedited) pedited->raw.exTwoGreen = true;
        }
    }
    catch( Exiv2::AnyError &e) {
        printf("loadFromXMP error: %s\n",e.what());
        return 2;
    }
    return 0;
}

int PartialProfile::saveParams ( Glib::ustring fname ) const
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

int PartialProfile::loadParams( Glib::ustring fname )
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

int PartialProfile::write (Glib::ustring &fname, Glib::ustring &content) const {

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

int PartialProfile::load (Glib::ustring fname, int *rank) {

    if (!pparams || fname.empty())
        return 1;

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

    pparams->ppVersion = PPVERSION;
    pparams->appVersion = APPVERSION;

if (keyFile.has_group ("Version")) {
    if (keyFile.has_key ("Version", "AppVersion")) pparams->appVersion = keyFile.get_string  ("Version", "AppVersion");
    if (keyFile.has_key ("Version", "Version"))    pparams->ppVersion  = keyFile.get_integer ("Version", "Version");
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
}

if (keyFile.has_group ("Exposure")) {
    if (pparams->ppVersion<PPVERSION_AEXP)
    	pparams->toneCurve.autoexp = false; // prevent execution of autoexp when opening file created with earlier verions of autoexp algorithm
    else
        if (keyFile.has_key ("Exposure", "Auto"))           { pparams->toneCurve.autoexp       = keyFile.get_boolean ("Exposure", "Auto"); if (pedited) pedited->toneCurve.autoexp = true; }

    if (keyFile.has_key ("Exposure", "Clip"))           { pparams->toneCurve.clip          = keyFile.get_double  ("Exposure", "Clip"); if (pedited) pedited->toneCurve.clip = true; }
    if (keyFile.has_key ("Exposure", "Compensation"))   { pparams->toneCurve.expcomp       = keyFile.get_double  ("Exposure", "Compensation"); if (pedited) pedited->toneCurve.expcomp = true; }
    if (keyFile.has_key ("Exposure", "Brightness"))     { pparams->toneCurve.brightness    = keyFile.get_integer ("Exposure", "Brightness"); if (pedited) pedited->toneCurve.brightness = true; }
    if (keyFile.has_key ("Exposure", "Contrast"))       { pparams->toneCurve.contrast      = keyFile.get_integer ("Exposure", "Contrast"); if (pedited) pedited->toneCurve.contrast = true; }
    if (keyFile.has_key ("Exposure", "Saturation"))     { pparams->toneCurve.saturation    = keyFile.get_integer ("Exposure", "Saturation"); if (pedited) pedited->toneCurve.saturation = true; }
    if (keyFile.has_key ("Exposure", "Black"))          { pparams->toneCurve.black         = keyFile.get_integer ("Exposure", "Black"); if (pedited) pedited->toneCurve.black = true; }
    if (keyFile.has_key ("Exposure", "ShadowCompr"))    { pparams->toneCurve.shcompr       = keyFile.get_integer ("Exposure", "ShadowCompr"); if (pedited) pedited->toneCurve.shcompr = true; }
    if (keyFile.has_key ("Exposure", "HighlightCompr")) { pparams->toneCurve.hlcompr       = keyFile.get_integer ("Exposure", "HighlightCompr"); if (pedited) pedited->toneCurve.hlcompr = true; }
    if (keyFile.has_key ("Exposure", "HighlightComprThreshold")) { pparams->toneCurve.hlcomprthresh = keyFile.get_integer ("Exposure", "HighlightComprThreshold"); if (pedited) pedited->toneCurve.hlcomprthresh = true; }
    if (pparams->toneCurve.shcompr > 100) pparams->toneCurve.shcompr = 100; // older pp3 files can have values above 100.
    if (keyFile.has_key ("Exposure", "CurveMode"))      {
        Glib::ustring sMode = keyFile.get_string ("Exposure", "CurveMode");
        if      (sMode == "Standard")            pparams->toneCurve.curveMode = ToneCurveParams::TC_MODE_STD;
        else if (sMode == "FilmLike")            pparams->toneCurve.curveMode = ToneCurveParams::TC_MODE_FILMLIKE;
        else if (sMode == "SatAndValueBlending") pparams->toneCurve.curveMode = ToneCurveParams::TC_MODE_SATANDVALBLENDING;
        else if (sMode == "WeightedStd")         pparams->toneCurve.curveMode = ToneCurveParams::TC_MODE_WEIGHTEDSTD;
        if (pedited) pedited->toneCurve.curveMode = true;
    }
    if (keyFile.has_key ("Exposure", "CurveMode2"))      {
        Glib::ustring sMode = keyFile.get_string ("Exposure", "CurveMode2");
        if      (sMode == "Standard")            pparams->toneCurve.curveMode2 = ToneCurveParams::TC_MODE_STD;
        else if (sMode == "FilmLike")            pparams->toneCurve.curveMode2 = ToneCurveParams::TC_MODE_FILMLIKE;
        else if (sMode == "SatAndValueBlending") pparams->toneCurve.curveMode2 = ToneCurveParams::TC_MODE_SATANDVALBLENDING;
        else if (sMode == "WeightedStd")         pparams->toneCurve.curveMode2 = ToneCurveParams::TC_MODE_WEIGHTEDSTD;
        if (pedited) pedited->toneCurve.curveMode2 = true;
    }
    if (pparams->ppVersion>200)
    if (keyFile.has_key ("Exposure", "Curve"))          { pparams->toneCurve.curve         = keyFile.get_double_list ("Exposure", "Curve"); if (pedited) pedited->toneCurve.curve = true; }
    if (keyFile.has_key ("Exposure", "Curve2"))         { pparams->toneCurve.curve2        = keyFile.get_double_list ("Exposure", "Curve2"); if (pedited) pedited->toneCurve.curve2 = true; }
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
        memcpy (pparams->chmixer.red, rmix.data(), 3*sizeof(int));
        memcpy (pparams->chmixer.green, gmix.data(), 3*sizeof(int));
        memcpy (pparams->chmixer.blue, bmix.data(), 3*sizeof(int));
    }
}

    // load luma curve
if (keyFile.has_group ("Luminance Curve")) {
    if (keyFile.has_key ("Luminance Curve", "Brightness"))     { pparams->labCurve.brightness   = keyFile.get_integer ("Luminance Curve", "Brightness"); if (pedited) pedited->labCurve.brightness = true; }
    if (keyFile.has_key ("Luminance Curve", "Contrast"))       { pparams->labCurve.contrast     = keyFile.get_integer ("Luminance Curve", "Contrast"); if (pedited) pedited->labCurve.contrast = true; }
    if (keyFile.has_key ("Luminance Curve", "Chromaticity"))   { pparams->labCurve.chromaticity = keyFile.get_integer ("Luminance Curve", "Chromaticity"); if (pedited) pedited->labCurve.chromaticity = true; }

    if (PPVERSION < 303) {
    // transform AvoidColorClipping into AvoidColorShift
    if (keyFile.has_key ("Luminance Curve", "AvoidColorClipping"))        { pparams->labCurve.avoidcolorshift = keyFile.get_boolean ("Luminance Curve", "AvoidColorClipping"); if (pedited) pedited->labCurve.avoidcolorshift = true; }
    }
    else {
    if (keyFile.has_key ("Luminance Curve", "AvoidColorShift"))           { pparams->labCurve.avoidcolorshift = keyFile.get_boolean ("Luminance Curve", "AvoidColorShift"); if (pedited) pedited->labCurve.avoidcolorshift = true; }
    if (keyFile.has_key ("Luminance Curve", "RedAndSkinTonesProtection")) { pparams->labCurve.rstprotection   = keyFile.get_double  ("Luminance Curve", "RedAndSkinTonesProtection"); if (pedited) pedited->labCurve.rstprotection = true; }
    }

    if (keyFile.has_key ("Luminance Curve", "LCredsk"))         { pparams->labCurve.lcredsk            = keyFile.get_boolean     ("Luminance Curve", "LCredsk"); if (pedited) pedited->labCurve.lcredsk = true; }
    if (keyFile.has_key ("Luminance Curve", "BWtoning"))        { pparams->labCurve.bwtoning           = keyFile.get_boolean     ("Luminance Curve", "BWtoning"); if (pedited) pedited->labCurve.bwtoning = true; }
    if (keyFile.has_key ("Luminance Curve", "LCurve"))          { pparams->labCurve.lcurve             = keyFile.get_double_list ("Luminance Curve", "LCurve"); if (pedited) pedited->labCurve.lcurve = true; }
    if (keyFile.has_key ("Luminance Curve", "aCurve"))          { pparams->labCurve.acurve             = keyFile.get_double_list ("Luminance Curve", "aCurve"); if (pedited) pedited->labCurve.acurve = true; }
    if (keyFile.has_key ("Luminance Curve", "bCurve"))          { pparams->labCurve.bcurve             = keyFile.get_double_list ("Luminance Curve", "bCurve"); if (pedited) pedited->labCurve.bcurve = true; }
    if (keyFile.has_key ("Luminance Curve", "ccCurve"))         { pparams->labCurve.cccurve            = keyFile.get_double_list ("Luminance Curve", "ccCurve"); if (pedited) pedited->labCurve.cccurve = true; }
    if (keyFile.has_key ("Luminance Curve", "chCurve"))         { pparams->labCurve.chcurve            = keyFile.get_double_list ("Luminance Curve", "chCurve"); if (pedited) pedited->labCurve.chcurve = true; }
    if (keyFile.has_key ("Luminance Curve", "LcCurve"))         { pparams->labCurve.lccurve            = keyFile.get_double_list ("Luminance Curve", "LcCurve"); if (pedited) pedited->labCurve.lccurve = true; }
}

    // load sharpening
if (keyFile.has_group ("Sharpening")) {
    if (keyFile.has_key ("Sharpening", "Enabled"))              { pparams->sharpening.enabled          = keyFile.get_boolean ("Sharpening", "Enabled"); if (pedited) pedited->sharpening.enabled = true; }
    if (keyFile.has_key ("Sharpening", "Radius"))               { pparams->sharpening.radius           = keyFile.get_double  ("Sharpening", "Radius"); if (pedited) pedited->sharpening.radius = true; }
    if (keyFile.has_key ("Sharpening", "Amount"))               { pparams->sharpening.amount           = keyFile.get_integer ("Sharpening", "Amount"); if (pedited) pedited->sharpening.amount = true; }
    if (keyFile.has_key ("Sharpening", "Threshold"))            {
        if (pparams->ppVersion < 302) {
            int thresh = min(keyFile.get_integer ("Sharpening", "Threshold"), 2000);
            pparams->sharpening.threshold.setValues(thresh, thresh, 2000, 2000); // TODO: 2000 is the maximum value and is taken of rtgui/sharpening.cc ; should be changed by the tool modularization
        }
        else {
            Glib::ArrayHandle<int> thresh = keyFile.get_integer_list ("Sharpening", "Threshold");
            pparams->sharpening.threshold.setValues(thresh.data()[0], thresh.data()[1], min(thresh.data()[2], 2000), min(thresh.data()[3], 2000));
        }
        if (pedited) pedited->sharpening.threshold = true;
    }
    if (keyFile.has_key ("Sharpening", "OnlyEdges"))            { pparams->sharpening.edgesonly        = keyFile.get_boolean ("Sharpening", "OnlyEdges"); if (pedited) pedited->sharpening.edgesonly = true; }
    if (keyFile.has_key ("Sharpening", "EdgedetectionRadius"))  { pparams->sharpening.edges_radius     = keyFile.get_double  ("Sharpening", "EdgedetectionRadius"); if (pedited) pedited->sharpening.edges_radius = true; }
    if (keyFile.has_key ("Sharpening", "EdgeTolerance"))        { pparams->sharpening.edges_tolerance  = keyFile.get_integer ("Sharpening", "EdgeTolerance"); if (pedited) pedited->sharpening.edges_tolerance = true; }
    if (keyFile.has_key ("Sharpening", "HalocontrolEnabled"))   { pparams->sharpening.halocontrol      = keyFile.get_boolean ("Sharpening", "HalocontrolEnabled"); if (pedited) pedited->sharpening.halocontrol = true; }
    if (keyFile.has_key ("Sharpening", "HalocontrolAmount"))    { pparams->sharpening.halocontrol_amount = keyFile.get_integer ("Sharpening", "HalocontrolAmount"); if (pedited) pedited->sharpening.halocontrol_amount = true; }
    if (keyFile.has_key ("Sharpening", "Method"))               { pparams->sharpening.method           = keyFile.get_string  ("Sharpening", "Method"); if (pedited) pedited->sharpening.method = true; }
    if (keyFile.has_key ("Sharpening", "DeconvRadius"))         { pparams->sharpening.deconvradius     = keyFile.get_double  ("Sharpening", "DeconvRadius"); if (pedited) pedited->sharpening.deconvradius = true; }
    if (keyFile.has_key ("Sharpening", "DeconvAmount"))         { pparams->sharpening.deconvamount     = keyFile.get_integer ("Sharpening", "DeconvAmount"); if (pedited) pedited->sharpening.deconvamount = true; }
    if (keyFile.has_key ("Sharpening", "DeconvDamping"))        { pparams->sharpening.deconvdamping    = keyFile.get_integer ("Sharpening", "DeconvDamping"); if (pedited) pedited->sharpening.deconvdamping = true; }
    if (keyFile.has_key ("Sharpening", "DeconvIterations"))     { pparams->sharpening.deconviter       = keyFile.get_integer ("Sharpening", "DeconvIterations"); if (pedited) pedited->sharpening.deconviter = true; }
}

    // load edge sharpening
if (keyFile.has_group ("SharpenEdge")) {
    if (keyFile.has_key ("SharpenEdge", "Enabled"))             { pparams->sharpenEdge.enabled         = keyFile.get_boolean ("SharpenEdge", "Enabled"); if (pedited) pedited->sharpenEdge.enabled = true; }
    if (keyFile.has_key ("SharpenEdge", "Passes"))              { pparams->sharpenEdge.passes          = keyFile.get_integer  ("SharpenEdge", "Passes"); if (pedited) pedited->sharpenEdge.passes = true; }
    if (keyFile.has_key ("SharpenEdge", "Strength"))            { pparams->sharpenEdge.amount          = keyFile.get_double  ("SharpenEdge", "Strength"); if (pedited) pedited->sharpenEdge.amount = true; }
    if (keyFile.has_key ("SharpenEdge", "ThreeChannels"))       { pparams->sharpenEdge.threechannels   = keyFile.get_boolean ("SharpenEdge", "ThreeChannels"); if (pedited) pedited->sharpenEdge.threechannels = true; }
}

    // load micro-contrast sharpening
if (keyFile.has_group ("SharpenMicro")) {
    if (keyFile.has_key ("SharpenMicro", "Enabled"))            { pparams->sharpenMicro.enabled        = keyFile.get_boolean ("SharpenMicro", "Enabled"); if (pedited) pedited->sharpenMicro.enabled = true; }
    if (keyFile.has_key ("SharpenMicro", "Matrix"))             { pparams->sharpenMicro.matrix         = keyFile.get_boolean ("SharpenMicro", "Matrix"); if (pedited) pedited->sharpenMicro.matrix = true; }
    if (keyFile.has_key ("SharpenMicro", "Strength"))           { pparams->sharpenMicro.amount         = keyFile.get_double  ("SharpenMicro", "Strength"); if (pedited) pedited->sharpenMicro.amount = true; }
    if (keyFile.has_key ("SharpenMicro", "Uniformity"))         { pparams->sharpenMicro.uniformity     = keyFile.get_double  ("SharpenMicro", "Uniformity"); if (pedited) pedited->sharpenMicro.uniformity = true; }
}

    // load vibrance
if (keyFile.has_group ("Vibrance")) {
    if (keyFile.has_key ("Vibrance", "Enabled"))                { pparams->vibrance.enabled            = keyFile.get_boolean ("Vibrance", "Enabled"); if (pedited) pedited->vibrance.enabled = true; }
    if (keyFile.has_key ("Vibrance", "Pastels"))                { pparams->vibrance.pastels            = keyFile.get_integer ("Vibrance", "Pastels"); if (pedited) pedited->vibrance.pastels = true; }
    if (keyFile.has_key ("Vibrance", "Saturated"))              { pparams->vibrance.saturated          = keyFile.get_integer ("Vibrance", "Saturated"); if (pedited) pedited->vibrance.saturated = true; }
    if (keyFile.has_key ("Vibrance", "PSThreshold"))            {
        if (pparams->ppVersion < 302) {
            int thresh = keyFile.get_integer ("Vibrance", "PSThreshold");
            pparams->vibrance.psthreshold.setValues(thresh, thresh);
        }
        else {
            Glib::ArrayHandle<int> thresh = keyFile.get_integer_list ("Vibrance", "PSThreshold");
            pparams->vibrance.psthreshold.setValues(thresh.data()[0], thresh.data()[1]);
        }
        if (pedited) pedited->vibrance.psthreshold = true;
    }
    if (keyFile.has_key ("Vibrance", "ProtectSkins"))           { pparams->vibrance.protectskins       = keyFile.get_boolean ("Vibrance", "ProtectSkins"); if (pedited) pedited->vibrance.protectskins = true; }
    if (keyFile.has_key ("Vibrance", "AvoidColorShift"))        { pparams->vibrance.avoidcolorshift    = keyFile.get_boolean ("Vibrance", "AvoidColorShift"); if (pedited) pedited->vibrance.avoidcolorshift = true; }
    if (keyFile.has_key ("Vibrance", "PastSatTog"))             { pparams->vibrance.pastsattog         = keyFile.get_boolean ("Vibrance", "PastSatTog"); if (pedited) pedited->vibrance.pastsattog = true; }
    if (keyFile.has_key ("Vibrance", "SkinTonesCurve"))         { pparams->vibrance.skintonescurve     = keyFile.get_double_list ("Vibrance", "SkinTonesCurve"); if (pedited) pedited->vibrance.skintonescurve = true; }
}

    // load wb
if (keyFile.has_group ("White Balance")) {
    if (keyFile.has_key ("White Balance", "Setting"))     { pparams->wb.method         = keyFile.get_string ("White Balance", "Setting"); if (pedited) pedited->wb.method = true; }
    if (keyFile.has_key ("White Balance", "Temperature")) { pparams->wb.temperature    = keyFile.get_integer ("White Balance", "Temperature"); if (pedited) pedited->wb.temperature = true; }
    if (keyFile.has_key ("White Balance", "Green"))       { pparams->wb.green          = keyFile.get_double ("White Balance", "Green"); if (pedited) pedited->wb.green = true; }
}

    // load defringe
if (keyFile.has_group ("Defringing")) {
    if (keyFile.has_key ("Defringing", "Enabled"))        { pparams->defringe.enabled   = keyFile.get_boolean ("Defringing", "Enabled"); if (pedited) pedited->defringe.enabled = true; }
    if (keyFile.has_key ("Defringing", "Radius"))         { pparams->defringe.radius    = keyFile.get_double  ("Defringing", "Radius"); if (pedited) pedited->defringe.radius = true; }
    if (keyFile.has_key ("Defringing", "Threshold"))      { pparams->defringe.threshold = keyFile.get_integer ("Defringing", "Threshold"); if (pedited) pedited->defringe.threshold = true; }
}

    // load impulseDenoise
if (keyFile.has_group ("Impulse Denoising")) {
    if (keyFile.has_key ("Impulse Denoising", "Enabled"))   { pparams->impulseDenoise.enabled = keyFile.get_boolean ("Impulse Denoising", "Enabled"); if (pedited) pedited->impulseDenoise.enabled = true; }
    if (keyFile.has_key ("Impulse Denoising", "Threshold")) { pparams->impulseDenoise.thresh  = keyFile.get_integer ("Impulse Denoising", "Threshold"); if (pedited) pedited->impulseDenoise.thresh = true; }
}

    // load dirpyrDenoise
if (keyFile.has_group ("Directional Pyramid Denoising")) {
    if (keyFile.has_key ("Directional Pyramid Denoising", "Enabled"))    { pparams->dirpyrDenoise.enabled = keyFile.get_boolean ("Directional Pyramid Denoising", "Enabled"); if (pedited) pedited->dirpyrDenoise.enabled = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "Luma"))       { pparams->dirpyrDenoise.luma    = keyFile.get_integer ("Directional Pyramid Denoising", "Luma"); if (pedited) pedited->dirpyrDenoise.luma = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "Ldetail"))    { pparams->dirpyrDenoise.Ldetail = keyFile.get_integer ("Directional Pyramid Denoising", "Ldetail"); if (pedited) pedited->dirpyrDenoise.Ldetail = true; }
    if (keyFile.has_key ("Directional Pyramid Denoising", "Chroma"))     { pparams->dirpyrDenoise.chroma  = keyFile.get_integer ("Directional Pyramid Denoising", "Chroma"); if (pedited) pedited->dirpyrDenoise.chroma = true;
    if (keyFile.has_key ("Directional Pyramid Denoising", "Gamma"))      { pparams->dirpyrDenoise.gamma   = keyFile.get_double  ("Directional Pyramid Denoising", "Gamma"); if (pedited) pedited->dirpyrDenoise.gamma = true; }
    }
}

    //Load EPD.
if (keyFile.has_group ("EPD")) {
    if(keyFile.has_key("EPD", "Enabled"))             { pparams->edgePreservingDecompositionUI.enabled = keyFile.get_boolean ("EPD", "Enabled"); if (pedited) pedited->edgePreservingDecompositionUI.enabled = true; }
    if(keyFile.has_key("EPD", "Strength"))            { pparams->edgePreservingDecompositionUI.Strength = keyFile.get_double ("EPD", "Strength"); if (pedited) pedited->edgePreservingDecompositionUI.Strength = true; }
    if(keyFile.has_key("EPD", "EdgeStopping"))        { pparams->edgePreservingDecompositionUI.EdgeStopping = keyFile.get_double ("EPD", "EdgeStopping"); if (pedited) pedited->edgePreservingDecompositionUI.EdgeStopping = true; }
    if(keyFile.has_key("EPD", "Scale"))               { pparams->edgePreservingDecompositionUI.Scale = keyFile.get_double ("EPD", "Scale"); if (pedited) pedited->edgePreservingDecompositionUI.Scale = true; }
    if(keyFile.has_key("EPD", "ReweightingIterates")) { pparams->edgePreservingDecompositionUI.ReweightingIterates = keyFile.get_integer ("EPD", "ReweightingIterates"); if (pedited) pedited->edgePreservingDecompositionUI.ReweightingIterates = true; }
}

    // load sh
if (keyFile.has_group ("Shadows & Highlights")) {
    if (keyFile.has_key ("Shadows & Highlights", "Enabled"))               { pparams->sh.enabled       = keyFile.get_boolean ("Shadows & Highlights", "Enabled"); if (pedited) pedited->sh.enabled = true; }
    if (keyFile.has_key ("Shadows & Highlights", "HighQuality"))           { pparams->sh.hq            = keyFile.get_boolean ("Shadows & Highlights", "HighQuality"); if (pedited) pedited->sh.hq = true; }
    if (keyFile.has_key ("Shadows & Highlights", "Highlights"))            { pparams->sh.highlights    = keyFile.get_integer ("Shadows & Highlights", "Highlights"); if (pedited) pedited->sh.highlights = true; }
    if (keyFile.has_key ("Shadows & Highlights", "HighlightTonalWidth"))   { pparams->sh.htonalwidth   = keyFile.get_integer ("Shadows & Highlights", "HighlightTonalWidth"); if (pedited) pedited->sh.htonalwidth = true; }
    if (keyFile.has_key ("Shadows & Highlights", "Shadows"))               { pparams->sh.shadows       = keyFile.get_integer ("Shadows & Highlights", "Shadows"); if (pedited) pedited->sh.shadows = true; }
    if (keyFile.has_key ("Shadows & Highlights", "ShadowTonalWidth"))      { pparams->sh.stonalwidth   = keyFile.get_integer ("Shadows & Highlights", "ShadowTonalWidth"); if (pedited) pedited->sh.stonalwidth = true; }
    if (keyFile.has_key ("Shadows & Highlights", "LocalContrast"))         { pparams->sh.localcontrast = keyFile.get_integer ("Shadows & Highlights", "LocalContrast"); if (pedited) pedited->sh.localcontrast = true; }
    if (keyFile.has_key ("Shadows & Highlights", "Radius"))                { pparams->sh.radius        = keyFile.get_integer ("Shadows & Highlights", "Radius"); if (pedited) pedited->sh.radius = true; }
}

    // load crop
if (keyFile.has_group ("Crop")) {
    if (keyFile.has_key ("Crop", "Enabled"))    { pparams->crop.enabled    = keyFile.get_boolean ("Crop", "Enabled"); if (pedited) pedited->crop.enabled = true; }
    if (keyFile.has_key ("Crop", "X"))          { pparams->crop.x          = keyFile.get_integer ("Crop", "X"); if (pedited) pedited->crop.x = true; }
    if (keyFile.has_key ("Crop", "Y"))          { pparams->crop.y          = keyFile.get_integer ("Crop", "Y"); if (pedited) pedited->crop.y = true; }
    if (keyFile.has_key ("Crop", "W"))          { pparams->crop.w          = keyFile.get_integer ("Crop", "W"); if (pedited) pedited->crop.w = true; }
    if (keyFile.has_key ("Crop", "H"))          { pparams->crop.h          = keyFile.get_integer ("Crop", "H"); if (pedited) pedited->crop.h = true; }
    if (keyFile.has_key ("Crop", "FixedRatio")) { pparams->crop.fixratio   = keyFile.get_boolean ("Crop", "FixedRatio"); if (pedited) pedited->crop.fixratio = true; }
    if (keyFile.has_key ("Crop", "Ratio")) {
        pparams->crop.ratio      = keyFile.get_string  ("Crop", "Ratio");
        if (pedited) pedited->crop.ratio = true;
        //backwards compatibility for crop.ratio
        if (pparams->crop.ratio=="DIN")    pparams->crop.ratio = "1.414 - DIN EN ISO 216";
        if (pparams->crop.ratio=="8.5:11") pparams->crop.ratio = "8.5:11 - US Letter";
        if (pparams->crop.ratio=="11:17")  pparams->crop.ratio = "11:17 - Tabloid";
    }
    if (keyFile.has_key ("Crop", "Orientation"))  { pparams->crop.orientation= keyFile.get_string  ("Crop", "Orientation"); if (pedited) pedited->crop.orientation = true; }
    if (keyFile.has_key ("Crop", "Guide"))        { pparams->crop.guide      = keyFile.get_string  ("Crop", "Guide"); if (pedited) pedited->crop.guide = true; }
}

    // load coarse
if (keyFile.has_group ("Coarse Transformation")) {
    if (keyFile.has_key ("Coarse Transformation", "Rotate"))          { pparams->coarse.rotate = keyFile.get_integer ("Coarse Transformation", "Rotate"); if (pedited) pedited->coarse.rotate = true; }
    if (keyFile.has_key ("Coarse Transformation", "HorizontalFlip"))  { pparams->coarse.hflip  = keyFile.get_boolean ("Coarse Transformation", "HorizontalFlip"); if (pedited) pedited->coarse.hflip = true; }
    if (keyFile.has_key ("Coarse Transformation", "VerticalFlip"))    { pparams->coarse.vflip  = keyFile.get_boolean ("Coarse Transformation", "VerticalFlip"); if (pedited) pedited->coarse.vflip = true; }
}

    // load rotate
if (keyFile.has_group ("Rotation")) {
    if (keyFile.has_key ("Rotation", "Degree"))   { pparams->rotate.degree = keyFile.get_double ("Rotation", "Degree"); if (pedited) pedited->rotate.degree = true; }
}
    // load commonTrans
if (keyFile.has_group ("Common Properties for Transformations")) {
    if (keyFile.has_key ("Common Properties for Transformations", "AutoFill"))   { pparams->commonTrans.autofill = keyFile.get_boolean ("Common Properties for Transformations", "AutoFill"); if (pedited) pedited->commonTrans.autofill = true; }
}

    // load distortion
if (keyFile.has_group ("Distortion")) {
    if (keyFile.has_key ("Distortion", "Amount"))     { pparams->distortion.amount     = keyFile.get_double  ("Distortion", "Amount"); if (pedited) pedited->distortion.amount = true; }
}

    // lens profile
if (keyFile.has_group ("LensProfile")) {
    if (keyFile.has_key ("LensProfile", "LCPFile")) { pparams->lensProf.lcpFile = keyFile.get_string ("LensProfile", "LCPFile"); if (pedited) pedited->lensProf.lcpFile = true; }
    if (keyFile.has_key ("LensProfile", "UseDistortion")) { pparams->lensProf.useDist = keyFile.get_boolean ("LensProfile", "UseDistortion"); if (pedited) pedited->lensProf.useDist = true; }
    if (keyFile.has_key ("LensProfile", "UseVignette")) { pparams->lensProf.useVign = keyFile.get_boolean ("LensProfile", "UseVignette"); if (pedited) pedited->lensProf.useVign = true; }
    if (keyFile.has_key ("LensProfile", "UseCA")) { pparams->lensProf.useCA = keyFile.get_boolean ("LensProfile", "UseCA"); if (pedited) pedited->lensProf.useCA = true; }
}
    
    // load perspective correction
if (keyFile.has_group ("Perspective")) {
    if (keyFile.has_key ("Perspective", "Horizontal"))  { pparams->perspective.horizontal = keyFile.get_integer ("Perspective", "Horizontal"); if (pedited) pedited->perspective.horizontal = true; }
    if (keyFile.has_key ("Perspective", "Vertical"))    { pparams->perspective.vertical   = keyFile.get_integer ("Perspective", "Vertical"); if (pedited) pedited->perspective.vertical = true; }
}

    // load c/a correction
if (keyFile.has_group ("CACorrection")) {
    if (keyFile.has_key ("CACorrection", "Red"))  { pparams->cacorrection.red  = keyFile.get_double ("CACorrection", "Red"); if (pedited) pedited->cacorrection.red = true; }
    if (keyFile.has_key ("CACorrection", "Blue")) { pparams->cacorrection.blue = keyFile.get_double ("CACorrection", "Blue"); if (pedited) pedited->cacorrection.blue = true; }
}

    // load vignetting correction
if (keyFile.has_group ("Vignetting Correction")) {
    if (keyFile.has_key ("Vignetting Correction", "Amount"))   { pparams->vignetting.amount = keyFile.get_integer ("Vignetting Correction", "Amount"); if (pedited) pedited->vignetting.amount = true; }
    if (keyFile.has_key ("Vignetting Correction", "Radius"))   { pparams->vignetting.radius = keyFile.get_integer ("Vignetting Correction", "Radius"); if (pedited) pedited->vignetting.radius = true; }
    if (keyFile.has_key ("Vignetting Correction", "Strength")) { pparams->vignetting.strength = keyFile.get_integer ("Vignetting Correction", "Strength"); if (pedited) pedited->vignetting.strength = true; }
    if (keyFile.has_key ("Vignetting Correction", "CenterX"))  { pparams->vignetting.centerX = keyFile.get_integer ("Vignetting Correction", "CenterX"); if (pedited) pedited->vignetting.centerX = true; }
    if (keyFile.has_key ("Vignetting Correction", "CenterY"))  { pparams->vignetting.centerY = keyFile.get_integer ("Vignetting Correction", "CenterY"); if (pedited) pedited->vignetting.centerY = true; }
}

    // load highlight recovery settings
if (keyFile.has_group ("HLRecovery")) {
    if (keyFile.has_key ("HLRecovery", "Enabled"))  { pparams->hlrecovery.enabled  = keyFile.get_boolean ("HLRecovery", "Enabled"); if (pedited) pedited->hlrecovery.enabled = true; }
    if (keyFile.has_key ("HLRecovery", "Method"))   { pparams->hlrecovery.method   = keyFile.get_string  ("HLRecovery", "Method"); if (pedited) pedited->hlrecovery.method = true; }
}
    // load resize settings
if (keyFile.has_group ("Resize")) {
    if (keyFile.has_key ("Resize", "Enabled"))       { pparams->resize.enabled   = keyFile.get_boolean ("Resize", "Enabled"); if (pedited) pedited->resize.enabled = true; }
    if (keyFile.has_key ("Resize", "Scale"))         { pparams->resize.scale     = keyFile.get_double ("Resize", "Scale"); if (pedited) pedited->resize.scale = true; }
    if (keyFile.has_key ("Resize", "AppliesTo"))     { pparams->resize.appliesTo = keyFile.get_string ("Resize", "AppliesTo"); if (pedited) pedited->resize.appliesTo = true; }
    if (keyFile.has_key ("Resize", "Method"))        { pparams->resize.method    = keyFile.get_string ("Resize", "Method"); if (pedited) pedited->resize.method = true; }
    if (keyFile.has_key ("Resize", "DataSpecified")) { pparams->resize.dataspec  = keyFile.get_integer ("Resize", "DataSpecified"); if (pedited) pedited->resize.dataspec = true; }
    if (keyFile.has_key ("Resize", "Width"))         { pparams->resize.width     = keyFile.get_integer ("Resize", "Width"); if (pedited) pedited->resize.width = true; }
    if (keyFile.has_key ("Resize", "Height"))        { pparams->resize.height    = keyFile.get_integer ("Resize", "Height"); if (pedited) pedited->resize.height = true; }
}

    // load color management settings
if (keyFile.has_group ("Color Management")) {
    if (keyFile.has_key ("Color Management", "InputProfile"))     { pparams->icm.input            = keyFile.get_string ("Color Management", "InputProfile"); if (pedited) pedited->icm.input = true; }
    if (keyFile.has_key ("Color Management", "ToneCurve"))        { pparams->icm.toneCurve        = keyFile.get_boolean ("Color Management", "ToneCurve"); if (pedited) pedited->icm.toneCurve = true; }
    if (keyFile.has_key ("Color Management", "BlendCMSMatrix"))   { pparams->icm.blendCMSMatrix   = keyFile.get_boolean ("Color Management", "BlendCMSMatrix"); if (pedited) pedited->icm.blendCMSMatrix = true; }
    if (keyFile.has_key ("Color Management", "PreferredProfile")) { pparams->icm.preferredProfile = keyFile.get_boolean ("Color Management", "PreferredProfile"); if (pedited) pedited->icm.preferredProfile = true; }
    if (keyFile.has_key ("Color Management", "WorkingProfile"))   { pparams->icm.working          = keyFile.get_string ("Color Management", "WorkingProfile"); if (pedited) pedited->icm.working = true; }
    if (keyFile.has_key ("Color Management", "OutputProfile"))    { pparams->icm.output           = keyFile.get_string ("Color Management", "OutputProfile"); if (pedited) pedited->icm.output = true; }
    if (keyFile.has_key ("Color Management", "Gammafree"))        { pparams->icm.gamma            = keyFile.get_string ("Color Management", "Gammafree"); if (pedited) pedited->icm.gamma = true; }
    if (keyFile.has_key ("Color Management", "Freegamma"))        { pparams->icm.freegamma        = keyFile.get_boolean ("Color Management", "Freegamma"); if (pedited) pedited->icm.freegamma = true; }
    if (keyFile.has_key ("Color Management", "GammaVal"))         { pparams->icm.gampos           = keyFile.get_double ("Color Management", "GammaVal"); if (pedited) pedited->icm.gampos = true; }
    if (keyFile.has_key ("Color Management", "GammaSlope"))       { pparams->icm.slpos            = keyFile.get_double ("Color Management", "GammaSlope"); if (pedited) pedited->icm.slpos = true; }

}

    // load directional pyramid equalizer parameters
if (keyFile.has_group ("Directional Pyramid Equalizer")) {
    if (keyFile.has_key ("Directional Pyramid Equalizer", "Enabled")) { pparams->dirpyrequalizer.enabled = keyFile.get_boolean ("Directional Pyramid Equalizer", "Enabled"); if (pedited) pedited->dirpyrequalizer.enabled = true; }
    for(int i = 0; i < 5; i ++) {
        std::stringstream ss;
        ss << "Mult" << i;
        if(keyFile.has_key ("Directional Pyramid Equalizer", ss.str())) { pparams->dirpyrequalizer.mult[i] = keyFile.get_double ("Directional Pyramid Equalizer", ss.str()); if (pedited) pedited->dirpyrequalizer.mult[i] = true; }
    }
}

    // load HSV equalizer parameters
if (keyFile.has_group ("HSV Equalizer")) {
    if (pparams->ppVersion>=300) {
        if (keyFile.has_key ("HSV Equalizer", "HCurve")) { pparams->hsvequalizer.hcurve = keyFile.get_double_list ("HSV Equalizer", "HCurve"); if (pedited) pedited->hsvequalizer.hcurve = true; }
        if (keyFile.has_key ("HSV Equalizer", "SCurve")) { pparams->hsvequalizer.scurve = keyFile.get_double_list ("HSV Equalizer", "SCurve"); if (pedited) pedited->hsvequalizer.scurve = true; }
        if (keyFile.has_key ("HSV Equalizer", "VCurve")) { pparams->hsvequalizer.vcurve = keyFile.get_double_list ("HSV Equalizer", "VCurve"); if (pedited) pedited->hsvequalizer.vcurve = true; }
    }
}

    // load RGB curves
if (keyFile.has_group ("RGB Curves")) {
    if (keyFile.has_key ("RGB Curves", "rCurve")) { pparams->rgbCurves.rcurve = keyFile.get_double_list ("RGB Curves", "rCurve"); if (pedited) pedited->rgbCurves.rcurve = true; }
    if (keyFile.has_key ("RGB Curves", "gCurve")) { pparams->rgbCurves.gcurve = keyFile.get_double_list ("RGB Curves", "gCurve"); if (pedited) pedited->rgbCurves.gcurve = true; }
    if (keyFile.has_key ("RGB Curves", "bCurve")) { pparams->rgbCurves.bcurve  = keyFile.get_double_list ("RGB Curves", "bCurve"); if (pedited) pedited->rgbCurves.bcurve = true; }
}

    // load raw settings
if (keyFile.has_group ("RAW")) {
    if (keyFile.has_key ("RAW", "DarkFrame"))        { pparams->raw.dark_frame = keyFile.get_string  ("RAW", "DarkFrame" ); if (pedited) pedited->raw.darkFrame = true; }
    if (keyFile.has_key ("RAW", "DarkFrameAuto"))    { pparams->raw.df_autoselect = keyFile.get_boolean ("RAW", "DarkFrameAuto" ); if (pedited) pedited->raw.dfAuto = true; }
    if (keyFile.has_key ("RAW", "FlatFieldFile"))       { pparams->raw.ff_file = keyFile.get_string  ("RAW", "FlatFieldFile" ); if (pedited) pedited->raw.ff_file = true; }
    if (keyFile.has_key ("RAW", "FlatFieldAutoSelect")) { pparams->raw.ff_AutoSelect = keyFile.get_boolean  ("RAW", "FlatFieldAutoSelect" );  if (pedited) pedited->raw.ff_AutoSelect = true; }
    if (keyFile.has_key ("RAW", "FlatFieldBlurRadius")) { pparams->raw.ff_BlurRadius = keyFile.get_integer  ("RAW", "FlatFieldBlurRadius" ); if (pedited) pedited->raw.ff_BlurRadius = true; }
    if (keyFile.has_key ("RAW", "FlatFieldBlurType"))   { pparams->raw.ff_BlurType = keyFile.get_string  ("RAW", "FlatFieldBlurType" ); if (pedited) pedited->raw.ff_BlurType = true; }
    if (keyFile.has_key ("RAW", "CA"))               { pparams->raw.ca_autocorrect = keyFile.get_boolean ("RAW", "CA" ); if (pedited) pedited->raw.caCorrection = true; }
    if (keyFile.has_key ("RAW", "CARed"))            { pparams->raw.cared = keyFile.get_double ("RAW", "CARed" ); if (pedited) pedited->raw.caRed = true; }
    if (keyFile.has_key ("RAW", "CABlue"))           { pparams->raw.cablue = keyFile.get_double ("RAW", "CABlue" ); if (pedited) pedited->raw.caBlue = true; }
    if (keyFile.has_key ("RAW", "HotDeadPixels"))    { pparams->raw.hotdeadpix_filt = keyFile.get_boolean ("RAW", "HotDeadPixels" ); if (pedited) pedited->raw.hotDeadPixelFilter = true; }
    if (keyFile.has_key ("RAW", "HotDeadPixelThresh")) { pparams->raw.hotdeadpix_thresh = keyFile.get_integer ("RAW", "HotDeadPixelThresh" ); if (pedited) pedited->raw.hotDeadPixelThresh = true; }
    if (keyFile.has_key ("RAW", "LineDenoise"))      { pparams->raw.linenoise = keyFile.get_integer ("RAW", "LineDenoise" ); if (pedited) pedited->raw.linenoise = true; }
    if (keyFile.has_key ("RAW", "GreenEqThreshold")) { pparams->raw.greenthresh= keyFile.get_integer ("RAW", "GreenEqThreshold"); if (pedited) pedited->raw.greenEq = true; }
    if (keyFile.has_key ("RAW", "CcSteps"))          { pparams->raw.ccSteps  = keyFile.get_integer ("RAW", "CcSteps"); if (pedited) pedited->raw.ccSteps = true; }
    if (keyFile.has_key ("RAW", "Method"))           { pparams->raw.dmethod = keyFile.get_string ("RAW", "Method"); if (pedited) pedited->raw.dmethod = true; }
    if (keyFile.has_key ("RAW", "DCBIterations"))    { pparams->raw.dcb_iterations = keyFile.get_integer("RAW", "DCBIterations"); if (pedited) pedited->raw.dcbIterations = true; }
    if (keyFile.has_key ("RAW", "DCBEnhance"))       { pparams->raw.dcb_enhance =keyFile.get_boolean("RAW", "DCBEnhance"); if (pedited) pedited->raw.dcbEnhance = true; }
    //if (keyFile.has_key ("RAW", "ALLEnhance"))       { pparams->raw.all_enhance =keyFile.get_boolean("RAW", "ALLEnhance"); if (pedited) pedited->raw.allEnhance = true; }

    if (keyFile.has_key ("RAW", "PreExposure"))   { pparams->raw.expos =keyFile.get_double("RAW", "PreExposure"); if (pedited) pedited->raw.exPos = true; }
    if (keyFile.has_key ("RAW", "PrePreserv"))    { pparams->raw.preser =keyFile.get_double("RAW", "PrePreserv"); if (pedited) pedited->raw.exPreser = true; }
    if (keyFile.has_key ("RAW", "PreBlackzero"))  { pparams->raw.blackzero =keyFile.get_double("RAW", "PreBlackzero"); if (pedited) pedited->raw.exBlackzero = true; }
    if (keyFile.has_key ("RAW", "PreBlackone"))   { pparams->raw.blackone =keyFile.get_double("RAW", "PreBlackone"); if (pedited) pedited->raw.exBlackone = true; }
    if (keyFile.has_key ("RAW", "PreBlacktwo"))   { pparams->raw.blacktwo =keyFile.get_double("RAW", "PreBlacktwo"); if (pedited) pedited->raw.exBlacktwo = true; }
    if (keyFile.has_key ("RAW", "PreBlackthree")) { pparams->raw.blackthree =keyFile.get_double("RAW", "PreBlackthree"); if (pedited) pedited->raw.exBlackthree = true; }
    if (keyFile.has_key ("RAW", "PreTwoGreen"))   { pparams->raw.twogreen =keyFile.get_boolean("RAW", "PreTwoGreen"); if (pedited) pedited->raw.exTwoGreen = true; }

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
    return 0;
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

/*
 * Reset a partial profile to default values if the instances exists; the content of the pedited class
 * can be initialized to a specific value through 'peditedValuesSetTo' (default: false)
 */
void PartialProfile::reset (bool peditedValuesSetTo) {
	if (pparams) pparams->setDefaults();
    if (pedited) pedited->set(peditedValuesSetTo);
}

/*
 * The (potentially partial) parameter set is copied to destParams
 */
void PartialProfile::applyTo(ProcParams *destParams) const {
    if (destParams) {
        if (destParams && pparams && pedited) {
            pedited->combine(*destParams, *pparams, true);
        }
    }
#ifdef _DEBUG
    else
        printf("PartialProfile::applyTo / Error: not destParams provided\n");
#endif
}

/*
 * The (potentially partial) parameter set is copied to (potentially partial) destParams
 * The ParamsEdited structures are merged too ( "or" operator )
 */
void PartialProfile::applyTo(PartialProfile *destPProfile) const {
    if (destPProfile) {
        if (destPProfile->pparams && pparams && pedited) {
            pedited->combine(*destPProfile->pparams, *pparams, true);
        }
        if (destPProfile->pedited && pedited)
            // Merge both ParamsEdited structure
            *pedited |= *destPProfile->pedited;
	}
#ifdef _DEBUG
    else
        printf("PartialProfile::applyTo / Error: not destPProfile provided\n");
#endif
}

void PartialProfile::set(bool v) {
    if (pedited) pedited->set(v);
}

}
}

