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
#ifndef __PROCEVENT__
#define __PROCEVENT__

namespace rtengine {

enum ProcEvent {
    EvPhotoLoaded=0,
    EvProfileLoaded=1,
    EvProfileChanged=2,
    EvHistoryBrowsed=3,
    EvBrightness=4,
    EvContrast=5,
    EvBlack=6,
    EvExpComp=7,
    EvHLCompr=8,
    EvSHCompr=9,
    EvToneCurve=10,
    EvAutoExp=11,
    EvClip=12,
    EvLBrightness=13,
    EvLContrast=14,
    EvLBlack=15,
    EvLHLCompr=16,
    EvLSHCompr=17,
    EvLLCurve=18,
    EvShrEnabled=19,
    EvShrRadius=20,
    EvShrAmount=21,
    EvShrThresh=22,
    EvShrEdgeOnly=23,
    EvShrEdgeRadius=24,
    EvShrEdgeTolerance=25,
    EvShrHaloControl=26,
    EvShrHaloAmount=27,
    EvShrMethod=28,
    EvShrDRadius=29,
    EvShrDAmount=30,
    EvShrDDamping=31,
    EvShrDIterations=32,
    EvCBAvoidClip=33, // obsolete
    EvCBSatLimiter=34,// obsolete
    EvCBSatLimit=35,  // obsolete
    EvCBBoost=36,     // obsolete
    EvWBMethod=37,
    EvWBTemp=38,
    EvWBGreen=39,
    EvCShiftA=40,     // obsolete
    EvCShiftB=41,     // obsolete
    EvLDNEnabled=42,  // obsolete
    EvLDNRadius=43,   // obsolete
    EvLDNEdgeTolerance=44, // obsolete
    EvCDNEnabled=45,  // obsolete
    EvCDNRadius=46,   // obsolete
    EvCDNEdgeTolerance=47, // obsolete
    EvCDNEdgeSensitive=48, // obsolete
    EvSHEnabled=49,
    EvSHHighlights=50,
    EvSHShadows=51,
    EvSHHLTonalW=52,
    EvSHSHTonalW=53,
    EvSHLContrast=54,
    EvSHRadius=55,
    EvCTRotate=56,
    EvCTHFlip=57,
    EvCTVFlip=58,
    EvROTDegree=59,
    EvTransAutoFill=60,
    EvDISTAmount=61,
    EvBookmarkSelected=62,
    EvCrop=63,
    EvCACorr=64,
    EvHREnabled=65,
    EvHRAmount=66,
    EvHRMethod=67,
    EvWProfile=68,
    EvOProfile=69,
    EvIProfile=70,
    EvVignetting=71,
    EvChMixer=72,
    EvResizeScale=73,
    EvResizeMethod=74,
    EvExif=75,
    EvIPTC=76,
    EvResizeSpec=77,
    EvResizeWidth=78,
    EvResizeHeight=79,
    EvResizeEnabled=80,
    EvProfileChangeNotification=81,
    EvSHHighQuality=82,
    EvPerspCorr=83,
    EvEqualizer=84, // obsolete
    EvEqlEnabled=85,// obsolete
    EvIDNEnabled=86,
    EvIDNThresh=87,
    EvDPDNEnabled=88,
    EvDPDNLuma=89,
    EvDPDNChroma=90,
    EvDPDNGamma=91,
    EvDirPyrEqualizer=92,
    EvDirPyrEqlEnabled=93,
    EvLSaturation=94,
    EvLaCurve=95,
    EvLbCurve=96,
    EvDemosaicMethod=97,
    EvPreProcessHotDeadPixel=98,
    EvSaturation=99,
    EvHSVEqualizerH=100,
    EvHSVEqualizerS=101,
    EvHSVEqualizerV=102,
    EvHSVEqEnabled=103,
    EvDefringeEnabled=104,
    EvDefringeRadius=105,
    EvDefringeThreshold=106,
    EvHLComprThreshold=107,
    EvResizeBoundingBox=108,
    EvResizeAppliesTo=109,
    EvLAvoidClip=110,
    EvLSatLimiter=111,
    EvLSatLimit=112,
    EvDemosaicDCBIter=113,
    EvDemosaicFalseColorIter=114,
    EvDemosaicDCBEnhanced=115,
    EvPreProcessCARed=116,
    EvPreProcessCABlue=117,
    EvPreProcessLineDenoise=118,
    EvPreProcessGEquilThresh=119,
    EvPreProcessAutoCA=120,
    EvPreProcessAutoDF=121,
    EvPreProcessDFFile=122,
    EvPreProcessExpCorrLinear=123,
    EvPreProcessExpCorrPH=124,
    EvFlatFieldFile=125,
    EvFlatFieldAutoSelect=126,
    EvFlatFieldBlurRadius=127,
    EvFlatFieldBlurType=128,
	EvAutoDIST=129,
	EvDPDNLumCurve=130,
	EvDPDNChromCurve=131,
	EvGAMMA=132,
	EvGAMPOS=133,
	EvGAMFREE=134,
	EvSLPOS=135,
	EvPreProcessExpBlackzero=136,
	EvPreProcessExpBlackone=137,
	EvPreProcessExpBlacktwo=138,
	EvPreProcessExpBlackthree=139,
	EvPreProcessExptwoGreen=140,
	EvSharpenEdgePasses=141,
	EvSharpenEdgeAmount=142,
	EvSharpenMicroAmount=143,
	EvSharpenMicroUniformity=144,
	EvSharpenEdgeEnabled=145,
	EvSharpenEdgeThreechannels=146,
	EvSharpenMicroEnabled=147,
	EvSharpenMicroMatrix=148,
	EvDemosaicALLEnhanced=149,
	EvVibranceEnabled=150,
	EvVibrancePastels=151,
	EvVibranceSaturated=152,
	EvVibranceProtectSkins=153,
	EvVibranceAvoidColorShift=154,
	EvVibrancePastSatTog=155,
	EvVibrancePastSatThreshold=156,
	EvEPDStrength=157,
	EvEPDEdgeStopping=158,
	EvEPDScale=159,
	EvEPDReweightingIterates=160,
	EvEPDEnabled=161,
	EvRGBrCurve=162,
	EvRGBgCurve=163,
	EvRGBbCurve=164,
	EvNeutralExp=165,
	NUMOFEVENTS=166
};
}    
#endif    
