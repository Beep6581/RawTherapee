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







// Aligned so the first entry starts on line 30.
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
    ALLNORAW,         // EvWBMethod,
    ALLNORAW,         // EvWBTemp,
    ALLNORAW,         // EvWBGreen,
    RGBCURVE,         // EvToneCurveMode1,
    RGBCURVE,         // EvToneCurve2,
    RGBCURVE,         // EvToneCurveMode2,
    0,                // EvLDNRadius: obsolete,
    0,                // EvLDNEdgeTolerance: obsolete,
    0,                // EvCDNEnabled:obsolete,
    0,                // free entry
    RGBCURVE,         // EvDCPToneCurve,
    ALLNORAW,         // EvDCPIlluminant,
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
    DEMOSAIC,         // EvWProfile,
    OUTPUTPROFILE,    // EvOProfile,
    ALLNORAW,         // EvIProfile,
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
    ALLNORAW,         // EvDirPyrEqualizer,
    ALLNORAW,         // EvDirPyrEqlEnabled,
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
    ALLNORAW,         // EvDemosaicFalseColorIter
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
    DEMOSAIC,         // EvDemosaicALLEnhanced Disabled but not removed for now, may be reintroduced some day
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
    DEMOSAIC | M_PREPROC, // EvDemosaicMethodPreProc
    LUMINANCECURVE,   // EvLCCurve
    LUMINANCECURVE,   // EvLCHCurve
    RGBCURVE,         // EvVibranceSkinTonesCurve
    LUMINANCECURVE,   // EvLLCCurve
    LUMINANCECURVE,   // EvLLCredsk
    ALLNORAW,         // EvDPDNLdetail
    ALLNORAW,         // EvCATEnabled
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
    ALLNORAW,         // EvWBequal
    ALLNORAW,         // EvWBequalbo
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
    0,                // --unused--
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
    ALLNORAW,         // EvDirPyrEqualizerThreshold
    ALLNORAW,         // EvDPDNenhance
    RGBCURVE,         // EvBWMethodalg
    ALLNORAW,         // EvDirPyrEqualizerSkin
    ALLNORAW,         // EvDirPyrEqlgamutlab
    ALLNORAW,         // EvDirPyrEqualizerHueskin
    ALLNORAW,         // EvDPDNmedian
    ALLNORAW,         // EvDPDNmedmet
    RGBCURVE,         // EvColorToningEnabled
    RGBCURVE,         // EvColorToningColor
    RGBCURVE,         // EvColorToningOpacity
    RGBCURVE,         // EvColorToningCLCurve
    RGBCURVE,         // EvColorToningMethod
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
    RGBCURVE,         // EvColorToningStrength
    RGBCURVE,         // EvColorToningautosat
    ALLNORAW,         // EvDPDNmetmed
    ALLNORAW,         // EvDPDNrgbmet
    ALLNORAW,         // EvDPDNpasses
    FLATFIELD,        // EvFlatFieldClipControl
    FLATFIELD,        // EvFlatFieldAutoClipControl
    DARKFRAME,        // EvPreProcessExpBlackRed
    DARKFRAME,        // EvPreProcessExpBlackGreen
    DARKFRAME,        // EvPreProcessExpBlackBlue
    RGBCURVE,         // EvFilmSimulationEnabled
    RGBCURVE,         // EvFilmSimulationStrength
    RGBCURVE,         // EvFilmSimulationFilename
    ALLNORAW,         // EvDPDNLCurve
    ALLNORAW,         // EvDPDNsmet
    DARKFRAME,        // EvPreProcessDeadPixel
    ALLNORAW,         // EvDPDNCCCurve
    ALLNORAW,         // EvDPDNautochroma
    ALLNORAW,         // EvDPDNLmet
    ALLNORAW,         // EvDPDNCmet
    ALLNORAW,         // EvDPDNC2met
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
    DIRPYREQUALIZER,  // EvWavCCCurve
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
    DIRPYREQUALIZER,  // EvWavStrngth
    DIRPYREQUALIZER,  // EvWavdaubcoeffmet
    DIRPYREQUALIZER,  // EvWavedgreinf
    DIRPYREQUALIZER,  // EvWaveletch
    DIRPYREQUALIZER,  // EvWavCHSLmet
    DIRPYREQUALIZER,  // EvWavedgcont
    DIRPYREQUALIZER,  // EvWavEDmet
    DIRPYREQUALIZER,  // EvWavlev0nois
    DIRPYREQUALIZER,  // EvWavlev1nois
    DIRPYREQUALIZER,  // EvWavlev2nois
    DIRPYREQUALIZER,  // EvWavmedianlev
    DIRPYREQUALIZER,  // EvWavHHCurve
    DIRPYREQUALIZER,  // EvWavBackmet
    DIRPYREQUALIZER,  // EvWavedgedetect
    DIRPYREQUALIZER,  // EvWavlipst
    DIRPYREQUALIZER,  // EvWavedgedetectthr
    DIRPYREQUALIZER,  // EvWavedgedetectthr2
    DIRPYREQUALIZER,  // EvWavlinkedg
    DIRPYREQUALIZER,  // EvWavCHCurve
    DARKFRAME,        // EvPreProcessHotDeadThresh
    SHARPENING,       // EvEPDgamma
    DIRPYREQUALIZER,  // EvWavtmr
    DIRPYREQUALIZER,  // EvWavTMmet
    DIRPYREQUALIZER,  // EvWavtmrs
    DIRPYREQUALIZER,  // EvWavbalance
    DIRPYREQUALIZER,  // EvWaviter
    DIRPYREQUALIZER,  // EvWavgamma
    DIRPYREQUALIZER,  // EvWavCLCurve
    DIRPYREQUALIZER,  // EvWavopacity
    DIRPYREQUALIZER,  // EvWavBAmet
    DIRPYREQUALIZER,  // EvWavopacityWL
    RESIZE,           // EvPrShrEnabled
    RESIZE,           // EvPrShrRadius
    RESIZE,           // EvPrShrAmount
    RESIZE,           // EvPrShrThresh
    RESIZE,           // EvPrShrEdgeOnly
    RESIZE,           // EvPrShrEdgeRadius=375,
    RESIZE,           // EvPrShrEdgeTolerance=376,
    RESIZE,           // EvPrShrHaloControl=377,
    RESIZE,           // EvPrShrHaloAmount=378,
    RESIZE,           // EvPrShrMethod=379,
    RESIZE,           // EvPrShrDRadius=380,
    RESIZE,           // EvPrShrDAmount=381,
    RESIZE,           // EvPrShrDDamping=382,
    RESIZE,           // EvPrShrDIterations=383,
    DIRPYREQUALIZER,  // EvWavcbenab
    DIRPYREQUALIZER,  // EvWavgreenhigh
    DIRPYREQUALIZER,  // EvWavbluehigh
    DIRPYREQUALIZER,  // EvWavgreenmed
    DIRPYREQUALIZER,  // EvWavbluemed
    DIRPYREQUALIZER,  // EvWavgreenlow
    DIRPYREQUALIZER,  // EvWavbluelow
    DIRPYREQUALIZER,  // EvWavNeutral
    RGBCURVE,         // EvDCPApplyLookTable,
    RGBCURVE,         // EvDCPApplyBaselineExposureOffset,
    ALLNORAW,         // EvDCPApplyHueSatMap
    DIRPYREQUALIZER,  // EvWavenacont
    DIRPYREQUALIZER,  // EvWavenachrom
    DIRPYREQUALIZER,  // EvWavenaedge
    DIRPYREQUALIZER,  // EvWavenares
    DIRPYREQUALIZER,  // EvWavenafin
    DIRPYREQUALIZER,  // EvWavenatoning
    DIRPYREQUALIZER,  // EvWavenanoise
    DIRPYREQUALIZER,  // EvWavedgesensi
    DIRPYREQUALIZER,  // EvWavedgeampli
    DIRPYREQUALIZER,  // EvWavlev3nois
    DIRPYREQUALIZER,  // EvWavNPmet
    DEMOSAIC,         // EvretinexMethod
    RETINEX,          // EvLneigh
    RETINEX,          // EvLgain
    RETINEX,          // EvLoffs
    RETINEX,          // EvLstr
    RETINEX,          // EvLscal
    RETINEX,          // EvLvart
    DEMOSAIC,         // EvLCDCurve
    RETINEX,          // EvRetinextransmission
    DEMOSAIC,         // EvRetinexEnabled
    RETINEX,          // EvRetinexmedianmap
    RETINEX,          // EvLlimd
    DEMOSAIC,         // Evretinexcolorspace
    DEMOSAIC,         // EvLCDHCurve
    DEMOSAIC,         // Evretinexgamma
    DEMOSAIC,         // EvLgam
    DEMOSAIC,         // EvLslope
    RETINEX,          // EvLhighl
    0,                // --unused--
    DEMOSAIC,         // EvRetinexlhcurve
    OUTPUTPROFILE,    // EvOIntent
    MONITORTRANSFORM, // EvMonitorTransform: no history message
    RETINEX,          // EvLiter
    RETINEX,          // EvLgrad
    RETINEX,          // EvLgrads
    RETINEX,          // EvLhighlights
    RETINEX,          // EvLh_tonalwidth
    RETINEX,          // EvLshadows
    RETINEX,          // EvLs_tonalwidth
    RETINEX,          // EvLradius
    RETINEX,          // EvmapMethod
    DEMOSAIC,         // EvRetinexmapcurve
    DEMOSAIC,         // EvviewMethod
    ALLNORAW,         // EvcbdlMethod
    RETINEX,          // EvRetinexgaintransmission
    RETINEX,          // EvLskal
    OUTPUTPROFILE,    // EvOBPCompens
    ALLNORAW,         // EvWBtempBias
    DARKFRAME,        // EvRawImageNum
    DEMOSAIC,         // EvPixelShiftMotion
    DEMOSAIC,         // EvPixelShiftMotionCorrection
    DEMOSAIC,         // EvPixelShiftStddevFactorGreen
    DEMOSAIC,         // EvPixelShiftEperIso
    DEMOSAIC,         // EvPixelShiftNreadIso
    DEMOSAIC,         // EvPixelShiftPrnu
    DEMOSAIC,         // EvPixelshiftShowMotion
    DEMOSAIC,         // EvPixelshiftShowMotionMaskOnly
    DEMOSAIC,         // EvPixelShiftAutomatic
    DEMOSAIC,         // EvPixelShiftNonGreenHorizontal
    DEMOSAIC,         // EvPixelShiftNonGreenVertical
    DEMOSAIC,         // EvPixelShiftNonGreenCross
    DEMOSAIC,         // EvPixelShiftStddevFactorRed
    DEMOSAIC,         // EvPixelShiftStddevFactorBlue
    DEMOSAIC,         // EvPixelShiftNonGreenCross2
    DEMOSAIC,         // EvPixelShiftNonGreenAmaze
    DEMOSAIC,         // EvPixelShiftGreen
    DEMOSAIC,         // EvPixelShiftRedBlueWeight
    DEMOSAIC,         // EvPixelShiftBlur
    DEMOSAIC,         // EvPixelShiftSigma
    DEMOSAIC,         // EvPixelShiftSum
    DEMOSAIC,         // EvPixelShiftExp0
    DEMOSAIC,         // EvPixelShiftHoleFill
    DEMOSAIC,         // EvPixelShiftMedian
    DEMOSAIC,         // EvPixelShiftMedian3
    DEMOSAIC,         // EvPixelShiftMotionMethod
    DEMOSAIC,         // EvPixelShiftSmooth
    DEMOSAIC,         // EvPixelShiftLmmse
    DEMOSAIC,         // EvPixelShiftEqualBright
    DEMOSAIC,          // EvPixelShiftEqualBrightChannel
    LUMINANCECURVE,   // EvCATtempout
    LUMINANCECURVE,   // EvCATgreenout
    LUMINANCECURVE,   // EvCATybout
    LUMINANCECURVE,   // EvCATDegreeout
    LUMINANCECURVE,   // EvCATAutoDegreeout
    LUMINANCECURVE,   // EvCATtempsc
    LUMINANCECURVE,   // EvCATgreensc
    LUMINANCECURVE,   // EvCATybscen
    LUMINANCECURVE,   // EvCATAutoyb
    DARKFRAME,        // EvLensCorrMode
    DARKFRAME,        // EvLensCorrLensfunCamera
    DARKFRAME,        // EvLensCorrLensfunLens
    ALLNORAW,         // EvTMFattalEnabled
    HDR,              // EvTMFattalThreshold
    HDR,              // EvTMFattalAmount
    ALLNORAW,         // EvWBEnabled
    RGBCURVE,         // EvRGBEnabled
    LUMINANCECURVE,    // EvLEnabled
    DEMOSAIC,          // EvPixelShiftOneGreen
    LUMINANCECURVE,   // EvlocallablocX
    LUMINANCECURVE,   // EvlocallabCenter
    LUMINANCECURVE,   // EvlocallabDegree
    LUMINANCECURVE,   // Evlocallablightness
    LUMINANCECURVE,   // Evlocallabcontrast
    LUMINANCECURVE,   // Evlocallabchroma
    LUMINANCECURVE,   // Evlocallabtransit
    LUMINANCECURVE,   // Evlocallabavoid
    LUMINANCECURVE,   // EvlocallablocYT
    LUMINANCECURVE,   // EvlocallablocXL
    LUMINANCECURVE,   // EvlocallabSmet
    LUMINANCECURVE,   // Evlocallabinvers
    LUMINANCECURVE,   // Evlocallabradius
    LUMINANCECURVE,   // Evlocallabinversrad
    LUMINANCECURVE,   // Evlocallabstrength
    LUMINANCECURVE,   // Evlocallabsensi
    LUMINANCECURVE,   // EvlocallabretinexMethod
    LUMINANCECURVE,   // Evlocallabstr
    LUMINANCECURVE,   // Evlocallabneigh
    LUMINANCECURVE,   // Evlocallabvart
    LUMINANCECURVE,   // EvlocallabCTgainCurve
    LUMINANCECURVE,   // Evlocallabchrrt
    LUMINANCECURVE,   // Evlocallabinversret
    LUMINANCECURVE,   // Evlocallabsensih
    LUMINANCECURVE,   // Evlocallabnbspot
    LUMINANCECURVE,   // Evlocallabactivlum
    LUMINANCECURVE,   // Evlocallabanbspot
    LUMINANCECURVE,   // Evlocallabsharradius
    LUMINANCECURVE,   // Evlocallabsharamount
    LUMINANCECURVE,   // Evlocallabshardamping
    LUMINANCECURVE,   // Evlocallabshariter
    LUMINANCECURVE,   // Evlocallabsensis
    LUMINANCECURVE,   // Evlocallabinverssha
    LUMINANCECURVE,   // Evlocallabcircrad
    LUMINANCECURVE,   // Evlocallabthres
    LUMINANCECURVE,   // Evlocallabproxi
    LUMINANCECURVE,   // EvlocallabqualityMethod
    LUMINANCECURVE,   // Evlocallabnoiselumf
    LUMINANCECURVE,   // Evlocallabnoiselumc
    LUMINANCECURVE,   // Evlocallabnoisechrof
    LUMINANCECURVE,   // Evlocallabnoisechroc
    LUMINANCECURVE,   // EvlocallabThresho
    LUMINANCECURVE,   // EvlocallabEqualizer
    LUMINANCECURVE,   // Evlocallabsensicb
    LUMINANCECURVE,   // Evlocallabsensibn
    LUMINANCECURVE,   // Evlocallabstren
    LUMINANCECURVE,   // Evlocallabgamma
    LUMINANCECURVE,   // Evlocallabestop
    LUMINANCECURVE,   // Evlocallabscaltm
    LUMINANCECURVE,   // Evlocallabrewei
    LUMINANCECURVE,   // Evlocallabsensitm
    LUMINANCECURVE,   // EvlocallabCTgainCurverab
    LUMINANCECURVE,   // Evlocallabretrab
    LUMINANCECURVE,   // Evlocallabllshape
    LUMINANCECURVE,   // Evlocenacolor
    LUMINANCECURVE,   // Evlocenablur
    LUMINANCECURVE,   // Evlocenatonemap
    LUMINANCECURVE,   // Evlocenareti
    LUMINANCECURVE,   // Evlocenasharp
    LUMINANCECURVE,   // Evlocenacbdl
    LUMINANCECURVE,   // Evlocenadenoi
    LUMINANCECURVE,   // EvlocallabLHshape
    LUMINANCECURVE,   // Evlocallabcurvactiv
    LUMINANCECURVE,   // Evlocallabccshape
    LUMINANCECURVE,   // EvlocallabqualitycurveMethod
    LUMINANCECURVE,   // Evlocallabhueref
    LUMINANCECURVE,   // Evlocallabchromaref
    LUMINANCECURVE,   // Evlocallablumaref
    LUMINANCECURVE,   // EvlocallabHHshape
    LUMINANCECURVE,   // Evlocenavibrance
    LUMINANCECURVE, //EvlocallabSkinTonesCurve
    LUMINANCECURVE,  //EvlocallabProtectSkins
    LUMINANCECURVE, //EvlocallabAvoidColorShift
    LUMINANCECURVE, //EvlocallabPastSatTog
    LUMINANCECURVE, //EvlocallabPastels
    LUMINANCECURVE, //EvlocallabSaturated
    LUMINANCECURVE, // EvlocallabPastSatThreshold
    LUMINANCECURVE, //Evlocallabsensiv
    LUMINANCECURVE, //EvLocenaexpose
    LUMINANCECURVE, //Evlocallabexpcomp
    LUMINANCECURVE, //Evlocallabhlcompr
    LUMINANCECURVE, //Evlocallabhlcomprthresh
    LUMINANCECURVE, //Evlocallabblack
    LUMINANCECURVE, //Evlocallabshcompr
    LUMINANCECURVE, //Evlocallabsensiex
    LUMINANCECURVE,   //Evlocallabshape
    LUMINANCECURVE, //Evlocallabcenterbuf
    LUMINANCECURVE,   //Evlocallabadjblur
    LUMINANCECURVE,   //Evlocallabcutpast
    LUMINANCECURVE,   //Evlocallabchromacbdl
    LUMINANCECURVE,   //EvlocallabblurMethod
    LUMINANCECURVE,   //EvlocallabdustMethod
    LUMINANCECURVE,   //Evlocallablastdust
    LUMINANCECURVE,   // Evlocallabsobelref
    LUMINANCECURVE,   // Evlocallabexclumethod
    LUMINANCECURVE,   // Evlocallabsensiexclu
    LUMINANCECURVE,   // Evlocallabstruc
    LUMINANCECURVE,   // Evlocallabwarm
    LUMINANCECURVE,    // Evlocallabnoiselumdetail
    LUMINANCECURVE,    // Evlocallabnoisechrodetail
    LUMINANCECURVE,    // Evlocallabsensiden
    LUMINANCECURVE,    // Evlocallabhuerefblur
    LUMINANCECURVE,   // EvlocallabEnabled
    LUMINANCECURVE,   // EvlocallablocY
    LUMINANCECURVE,   // Evlocallabbilateral
    LUMINANCECURVE,    // Evlocallabnoiselequal
    LUMINANCECURVE   // Evlocallabshapemethod
};


namespace rtengine
{

RefreshMapper::RefreshMapper():
    next_event_(rtengine::NUMOFEVENTS)
{
    for (int event = 0; event < rtengine::NUMOFEVENTS; ++event) {
        actions_[event] = refreshmap[event];
    }
}


ProcEvent RefreshMapper::newEvent()
{
    return ProcEvent(++next_event_);
}


void RefreshMapper::mapEvent(ProcEvent event, int action)
{
    actions_[event] = action;
}


int RefreshMapper::getAction(ProcEvent event) const
{
    auto it = actions_.find(event);

    if (it == actions_.end()) {
        return 0;
    } else {
        return it->second;
    }
}


RefreshMapper *RefreshMapper::getInstance()
{
    static RefreshMapper instance;
    return &instance;
}

} // namespace rtengine
