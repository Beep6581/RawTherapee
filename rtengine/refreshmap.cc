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
RGBCURVE,         // EvToneCurve,
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
0,                // EvCBAvoidClip: obsolete
0,                // EvCBSatLimiter: obsolete
0,                // EvCBSatLimit: obsolete
0,                // EvCBBoost: obsolete
WHITEBALANCE,     // EvWBMethod,
WHITEBALANCE,     // EvWBTemp,
WHITEBALANCE,     // EvWBGreen,
0,                // EvCShiftA: obsolete
0,                // EvCShiftB: obsolete
0,                // EvLDNEnabled: obsolete,
0,                // EvLDNRadius: obsolete,
0,                // EvLDNEdgeTolerance: obsolete,
0,                // EvCDNEnabled:obsolete,
0,                // EvCDNRadius: obsolete,
0,                // EvCDNEdgeTolerance: obsolete,
0,                // EvCDNEdgeSensitive: obsolete,
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
TRANSFORM,        // EvVignetting,
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
0,                // EvEqualizer: obsolete
0,                // EvEqlEnabled:obsolete
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
DARKFRAME,        // EvPreProcessHotDeadPixel
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
ALL,              // EvGAMMA
ALL,              // EvGAMPOS
ALL,              // EvGAMFREE
ALL,              // EvSLPOS
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
DEMOSAIC,         // EvDemosaicALLEnhanced
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
ALLNORAW          // EvDPDNLdetail
	
};

