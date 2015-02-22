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
#include "refreshmap.h"
#include "procevents.h"

int refreshmap[rtengine::NUMOFEVENTS] = {
ALL,              // EvPhotoLoaded,
ALL,              // EvProfileLoaded,
ALL,              // EvProfileChanged,
ALL,              // EvHistoryBrowsed,
RGBCURVE,         // EvBrightness,
RGBCURVE,         // EvContrast,
RGBCURVE,         // EvBlack,
RGBCURVE,         // EvExpComp,
RGBCURVE,         // EvHLCompr,
RGBCURVE,         // EvSHCompr,
RGBCURVE,         // EvToneCurve1,
AUTOEXP,          // EvAutoExp,
AUTOEXP,          // EvClip,
LUMINANCECURVE,   // EvLBrightness,
LUMINANCECURVE,   // EvLContrast,
LUMINANCECURVE,   // EvLBlack,
LUMINANCECURVE,   // EvLHLCompr,
LUMINANCECURVE,   // EvLSHCompr,
LUMINANCECURVE,   // EvLLCurve,
SHARPENING,       // EvShrEnabled,
SHARPENING,       // EvShrRadius,
SHARPENING,       // EvShrAmount,
SHARPENING,       // EvShrThresh,
SHARPENING,       // EvShrEdgeOnly,
SHARPENING,       // EvShrEdgeRadius,
SHARPENING,       // EvShrEdgeTolerance,
SHARPENING,       // EvShrHaloControl,
SHARPENING,       // EvShrHaloAmount,
SHARPENING,       // EvShrMethod,
SHARPENING,       // EvShrDRadius,
SHARPENING,       // EvShrDAmount,
SHARPENING,       // EvShrDDamping,
SHARPENING,       // EvShrDIterations,
TRANSFORM,        // EvLCPUseDist,
DARKFRAME,        // EvLCPUseVign,
TRANSFORM,        // EvLCPUseCA,
M_VOID,           // EvFixedExp
WHITEBALANCE,     // EvWBMethod,
WHITEBALANCE,     // EvWBTemp,
WHITEBALANCE,     // EvWBGreen,
RGBCURVE,         // EvToneCurveMode1,
RGBCURVE,         // EvToneCurve2,
RGBCURVE,         // EvToneCurveMode2,
0,                // EvLDNRadius: obsolete,
0,                // EvLDNEdgeTolerance: obsolete,
0,                // EvCDNEnabled:obsolete,
ALL,              // EvBlendCMSMatrix,
ALL,              // EvDCPToneCurve,
ALL,              // EvDCPIlluminant,
RETINEX,          // EvSHEnabled,
RGBCURVE,         // EvSHHighlights,
RGBCURVE,         // EvSHShadows,
RGBCURVE,         // EvSHHLTonalW,
RGBCURVE,         // EvSHSHTonalW,
RGBCURVE,         // EvSHLContrast,
RETINEX,          // EvSHRadius,
ALL,              // EvCTRotate,
ALL,              // EvCTHFlip,
ALL,              // EvCTVFlip,
TRANSFORM,        // EvROTDegree,
TRANSFORM,        // EvTransAutoFill,
TRANSFORM,        // EvDISTAmount,
ALL,              // EvBookmarkSelected,
CROP,             // EvCrop,
TRANSFORM,        // EvCACorr,
ALLNORAW,         // EvHREnabled,
ALLNORAW,         // EvHRAmount,
ALLNORAW,         // EvHRMethod,
ALL,              // EvWProfile,
OUTPUTPROFIL,     // EvOProfile,
ALL,              // EvIProfile,
TRANSFORM,        // EvVignettingAmount,
RGBCURVE,         // EvChMixer,
RESIZE,           // EvResizeScale,
RESIZE,           // EvResizeMethod,
EXIF,             // EvExif,
IPTC,             // EvIPTC
RESIZE,           // EvResizeSpec,
RESIZE,           // EvResizeWidth
RESIZE,           // EvResizeHeight
RESIZE,           // EvResizeEnabled
ALL,              // EvProfileChangeNotification
RETINEX,          // EvShrHighQuality
TRANSFORM,        // EvPerspCorr
DARKFRAME,        // EvLCPFile
RGBCURVE,         // EvRGBrCurveLumamode
IMPULSEDENOISE,   // EvIDNEnabled,
IMPULSEDENOISE,   // EvIDNThresh,
ALLNORAW,         // EvDPDNEnabled,
ALLNORAW,         // EvDPDNLuma,
ALLNORAW,         // EvDPDNChroma,
ALLNORAW,         // EvDPDNGamma,
DIRPYREQUALIZER,  // EvDirPyrEqualizer,
DIRPYREQUALIZER,  // EvDirPyrEqlEnabled,
LUMINANCECURVE,   // EvLSaturation,
LUMINANCECURVE,   // EvLaCurve,
LUMINANCECURVE,   // EvLbCurve,
DEMOSAIC,         // EvDemosaicMethod
DARKFRAME,        // EvPreProcessHotPixel
RGBCURVE,         // EvSaturation,
RGBCURVE,         // EvHSVEqualizerH,
RGBCURVE,         // EvHSVEqualizerS,
RGBCURVE,         // EvHSVEqualizerV,
RGBCURVE,         // EvHSVEqEnabled,
DEFRINGE,         // EvDefringeEnabled,
DEFRINGE,         // EvDefringeRadius,
DEFRINGE,         // EvDefringeThreshold,
RGBCURVE,         // EvHLComprThreshold,
RESIZE,           // EvResizeBoundingBox
RESIZE,           // EvResizeAppliesTo
LUMINANCECURVE,   // EvCBAvoidClip,
LUMINANCECURVE,   // EvCBSatLimiter,
LUMINANCECURVE,   // EvCBSatLimit,
DEMOSAIC,         // EvDemosaicDCBIter
DEMOSAIC,         // EvDemosaicFalseColorIter
DEMOSAIC,         // EvDemosaicDCBEnhanced
DARKFRAME,        // EvPreProcessCARed
DARKFRAME,        // EvPreProcessCABlue
DARKFRAME,        // EvPreProcessLineDenoise
DARKFRAME,        // EvPreProcessGEquilThresh
DARKFRAME,        // EvPreProcessAutoCA
DARKFRAME,        // EvPreProcessAutoDF
DARKFRAME,        // EvPreProcessDFFile
DARKFRAME,        // EvPreProcessExpCorrLinear
DARKFRAME,        // EvPreProcessExpCorrPH
FLATFIELD,        // EvFlatFieldFile,
FLATFIELD,        // EvFlatFieldAutoSelect,
FLATFIELD,        // EvFlatFieldBlurRadius,
FLATFIELD,        // EvFlatFieldBlurType,
TRANSFORM,        // EvAutoDIST,
ALLNORAW,         // EvDPDNLumCurve,
ALLNORAW,         // EvDPDNChromCurve,
GAMMA,            // EvGAMMA
GAMMA,            // EvGAMPOS
GAMMA,            // EvGAMFREE
GAMMA,            // EvSLPOS
DARKFRAME,        // EvPreProcessExpBlackzero
DARKFRAME,        // EvPreProcessExpBlackone
DARKFRAME,        // EvPreProcessExpBlacktwo
DARKFRAME,        // EvPreProcessExpBlackthree
DARKFRAME,        // EvPreProcessExptwoGreen
SHARPENING,       // EvSharpenEdgePasses
SHARPENING,       // EvSharpenEdgeStrength
SHARPENING,       // EvSharpenMicroStrength
SHARPENING,       // EvSharpenMicroUniformity
SHARPENING,       // EvSharpenEdgeEnabled
SHARPENING,       // EvSharpenEdgeThreechannels
SHARPENING,       // EvSharpenMicroEnabled
SHARPENING,       // EvSharpenMicroMatrix
DEMOSAIC,         // EvDemosaicALLEnhanced  // Disabled but not removed for now, may be reintroduced some day
RGBCURVE,         // EvVibranceEnabled
RGBCURVE,         // EvVibrancePastels
RGBCURVE,         // EvVibranceSaturated
RGBCURVE,         // EvVibranceProtectSkins
RGBCURVE,         // EvVibranceAvoidColorShift
RGBCURVE,         // EvVibrancePastSatTog
RGBCURVE,         // EvVibrancePastSatThreshold
SHARPENING,       // EvEPDStrength
SHARPENING,       // EvEPDEdgeStopping
SHARPENING,       // EvEPDScale
SHARPENING,       // EvEPDReweightingIterates
SHARPENING,       // EvEPDEnabled
RGBCURVE,         // EvRGBrCurve
RGBCURVE,         // EvRGBgCurve
RGBCURVE,         // EvRGBbCurve
RGBCURVE,         // EvNeutralExp
DEMOSAIC|M_PREPROC, // EvDemosaicMethodPreProc
LUMINANCECURVE,   // EvLCCurve
LUMINANCECURVE,   // EvLCHCurve
RGBCURVE,         // EvVibranceSkinTonesCurve
LUMINANCECURVE,   // EvLLCCurve
LUMINANCECURVE,   // EvLLCredsk
ALLNORAW,         // EvDPDNLdetail
LUMINANCECURVE,   // EvCATEnabled
LUMINANCECURVE,   // EvCATDegree
LUMINANCECURVE,   // EvCATMethodsur
LUMINANCECURVE,   // EvCATAdapscen
LUMINANCECURVE,   // EvCATAdapLum
LUMINANCECURVE,   // EvCATMethodWB
LUMINANCECURVE,   // EvCATJLight
LUMINANCECURVE,   // EvCATChroma
LUMINANCECURVE,   // EvCATAutoDegree
LUMINANCECURVE,   // EvCATContrast
LUMINANCECURVE,   // EvCATSurr
LUMINANCECURVE,   // EvCATgamut
LUMINANCECURVE,   // EvCATmethodalg
LUMINANCECURVE,   // EvCATRstpro
LUMINANCECURVE,   // EvCATQbright
LUMINANCECURVE,   // EvCATQContrast
LUMINANCECURVE,   // EvCATSChroma
LUMINANCECURVE,   // EvCATMchroma
LUMINANCECURVE,   // EvCAThue
LUMINANCECURVE,   // EvCATcurve1
LUMINANCECURVE,   // EvCATcurve2
LUMINANCECURVE,   // EvCATcurvemode1
LUMINANCECURVE,   // EvCATcurvemode2
LUMINANCECURVE,   // EvCATcurve3
LUMINANCECURVE,   // EvCATcurvemode3
LUMINANCECURVE,   // EvCATdatacie
LUMINANCECURVE,   // EvCATtonecie
ALLNORAW,         // EvDPDNbluechro
ALLNORAW,         // EvDPDNperform
ALLNORAW,         // EvDPDNmet
DEMOSAIC,         // EvDemosaicLMMSEIter
LUMINANCECURVE,   // EvCATbadpix
LUMINANCECURVE,   // EvCATAutoadap
DEFRINGE,         // EvPFCurve
WHITEBALANCE,     // EvWBequal
WHITEBALANCE,     // EvWBequalbo
TRANSFORM,        // EvGradientDegree
TRANSFORM,        // EvGradientEnabled
TRANSFORM,        // EvPCVignetteStrength
TRANSFORM,        // EvPCVignetteEnabled
RGBCURVE,         // EvBWChmixEnabled
RGBCURVE,         // EvBWred
RGBCURVE,         // EvBWgreen
RGBCURVE,         // EvBWblue
RGBCURVE,         // EvBWredgam
RGBCURVE,         // EvBWgreengam
RGBCURVE,         // EvBWbluegam
RGBCURVE,         // EvBWfilter
RGBCURVE,         // EvBWsetting
RGBCURVE,         // EvBWoran
RGBCURVE,         // EvBWyell
RGBCURVE,         // EvBWcyan
RGBCURVE,         // EvBWmag
RGBCURVE,         // EvBpur
RGBCURVE,         // EvBWLuminanceEqual
RGBCURVE,         // EvBWChmixEnabledLm
RGBCURVE,         // EvBWmethod
RGBCURVE,         // EvBWBeforeCurve
RGBCURVE,         // EvBWBeforeCurveMode
RGBCURVE,         // EvBWAfterCurve
RGBCURVE,         // EvBWAfterCurveMode
RGBCURVE,         // EvAutoch
NONE,             // --unused--
RGBCURVE,         // EvNeutralBW
TRANSFORM,        // EvGradientFeather
TRANSFORM,        // EvGradientStrength
TRANSFORM,        // EvGradientCenter
TRANSFORM,        // EvPCVignetteFeather
TRANSFORM,        // EvPCVignetteRoundness
TRANSFORM,        // EvVignettingRadius,
TRANSFORM,        // EvVignettingStrength
TRANSFORM,        // EvVignettingCenter
LUMINANCECURVE,   // EvLCLCurve
LUMINANCECURVE,   // EvLLHCurve
LUMINANCECURVE,   // EvLHHCurve
DIRPYREQUALIZER,  // EvDirPyrEqualizerThreshold
ALLNORAW,         // EvDPDNenhance
RGBCURVE,         // EvBWMethodalg
DIRPYREQUALIZER,  // EvDirPyrEqualizerSkin
DIRPYREQUALIZER,  // EvDirPyrEqlgamutlab
DIRPYREQUALIZER,  // EvDirPyrEqualizerHueskin
ALLNORAW,         // EvDPDNmedian
ALLNORAW,         //EvDPDNmedmet
//DIRPYREQUALIZER   // EvDirPyrEqualizeralg
RGBCURVE,         // EvColorToningEnabled
RGBCURVE,         // EvColorToningColor
RGBCURVE,         // EvColorToningOpacity
RGBCURVE,         // EvColorToningCLCurve
RGBCURVE,         // EvColorToningMethod
//RGBCURVE,       // EvColorToningTwocolor
RGBCURVE,         // EvColorToningLLCurve
RGBCURVE,         // EvColorToningredlow
RGBCURVE,         // EvColorToninggreenlow
RGBCURVE,         // EvColorToningbluelow
RGBCURVE,         // EvColorToningredmed
RGBCURVE,         // EvColorToninggreenmed
RGBCURVE,         // EvColorToningbluemed
RGBCURVE,         // EvColorToningredhigh
RGBCURVE,         // EvColorToninggreenhigh
RGBCURVE,         // EvColorToningbluehigh
RGBCURVE,         // EvColorToningbalance
RGBCURVE,         // EvColorToningNeutral
RGBCURVE,         // EvColorToningsatlow
RGBCURVE,         // EvColorToningsathigh
RGBCURVE,         // EvColorToningTwocolor
RGBCURVE,         // EvColorToningNeutralcur
RGBCURVE,         // EvColorToningLumamode
RGBCURVE,         // EvColorToningShadows
RGBCURVE,         // EvColorToningHighights
RGBCURVE,         // EvColorToningSatProtection
RGBCURVE,         // EvColorToningSatThreshold
RGBCURVE,			//EvColorToningStrength
RGBCURVE,			//EvColorToningautosat
ALLNORAW, 			//EvDPDNmetmed
ALLNORAW, 			//EvDPDNrgbmet
ALLNORAW,			//EvDPDNpasses
FLATFIELD,        // EvFlatFieldClipControl
FLATFIELD,        // EvFlatFieldAutoClipControl
DARKFRAME,        // EvPreProcessExpBlackRed
DARKFRAME,        // EvPreProcessExpBlackGreen
DARKFRAME,        // EvPreProcessExpBlackBlue
RGBCURVE,         //EvFilmSimulationEnabled
RGBCURVE,         //EvFilmSimulationStrength
RGBCURVE,         //EvFilmSimulationFilename
ALLNORAW,			//	EvDPDNLCurve
ALLNORAW,			//	EvDPDNsmet
DARKFRAME,        // EvPreProcessDeadPixel
ALLNORAW,			//EvDPDNCCCurve
ALLNORAW,			//EvDPDNautochroma
ALLNORAW,			//	EvDPDNLmet
ALLNORAW,			//	EvDPDNCmet
ALLNORAW,			//	EvDPDNC2met
DIRPYREQUALIZER,  // EvWavelet
DIRPYREQUALIZER,  // EvEnabled
DIRPYREQUALIZER,  // EvWavLmethod
DIRPYREQUALIZER,  // EvWavCLmethod
DIRPYREQUALIZER,  // EvWavDirmethod
DIRPYREQUALIZER,  // EvWavtiles
DIRPYREQUALIZER,  // EvWavsky
DIRPYREQUALIZER,  // EvWavthres
DIRPYREQUALIZER,  // EvWavthr
DIRPYREQUALIZER,  // EvWavchroma
DIRPYREQUALIZER,  // EvWavmedian
DIRPYREQUALIZER,  // EvWavunif
DIRPYREQUALIZER,  // EvWavSkin
DIRPYREQUALIZER,  // EvWavHueSkin
DIRPYREQUALIZER,  // EvWavThreshold
DIRPYREQUALIZER,  // EvWavlhl
DIRPYREQUALIZER,  // EvWavbhl
DIRPYREQUALIZER,  // EvWavThresHold2
DIRPYREQUALIZER,  // EvWavavoid
DIRPYREQUALIZER,  // EvWavCLVCurve
DIRPYREQUALIZER,  // EvWavpast
DIRPYREQUALIZER,  // EvWavsat
DIRPYREQUALIZER,  // EvWavCHmet
DIRPYREQUALIZER,  // EvWavHSmet
DIRPYREQUALIZER,  // EvWavchro
DIRPYREQUALIZER,  // EvWavColor
DIRPYREQUALIZER,  // EvWavOpac
DIRPYREQUALIZER,  // EvWavsup
DIRPYREQUALIZER,  // EvWavTilesmet
DIRPYREQUALIZER,  // EvWavrescon
DIRPYREQUALIZER,  // EvWavreschro
DIRPYREQUALIZER,  // EvWavresconH
DIRPYREQUALIZER,  // EvWavthrH
DIRPYREQUALIZER,  // EvWavHueskin2
DIRPYREQUALIZER,  // EvWavedgrad
DIRPYREQUALIZER,  // EvWavedgval
DIRPYREQUALIZER  // EvWavStrngth

};

