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
    EvCBAvoidClip=33,
    EvCBSatLimiter=34,
    EvCBSatLimit=35,
    EvCBBoost=36,
    EvWBMethod=37,
    EvWBTemp=38,
    EvWBGreen=39,
    EvCShiftA=40,
    EvCShiftB=41,
    EvLDNEnabled=42,
    EvLDNRadius=43,
    EvLDNEdgeTolerance=44,
    EvCDNEnabled=45,
    EvCDNRadius=46,
    EvCDNEdgeTolerance=47,
    EvCDNEdgeSensitive=48,
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
    EvEqualizer=84,
    EvEqlEnabled=85,
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
	NUMOFEVENTS=141
	
	
	};
}    
#endif    
