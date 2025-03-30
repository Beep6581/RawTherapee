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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "refreshmap.h"
#include "procevents.h"







// Aligned so the first entry starts on line 30.
int refreshmap[rtengine::NUMOFEVENTS] = {
    ALL,              // EvPhotoLoaded,
    0,              // EvProfileLoaded : obsolete,
    ALL,              // EvProfileChanged,
    ALL,              // EvHistoryBrowsed,
    AUTOEXP,         // EvBrightness,
    AUTOEXP,         // EvContrast,
    AUTOEXP,         // EvBlack,
    AUTOEXP,         // EvExpComp,
    AUTOEXP,         // EvHLCompr,
    AUTOEXP,         // EvSHCompr,
    AUTOEXP,         // EvToneCurve1,
    AUTOEXP,          // EvAutoExp,
    AUTOEXP,          // EvClip,
    LUMINANCECURVE,   // EvLBrightness,
    LUMINANCECURVE,   // EvLContrast,
    0,   // EvLBlack : obsolete,
    0,   // EvLHLCompr : obsolete,
    0,   // EvLSHCompr : obsolete,
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
    HDR,        // EvLCPUseDist,
    DARKFRAME,        // EvLCPUseVign,
    HDR,        // EvLCPUseCA,
    M_VOID,           // EvFixedExp
    WB,               // EvWBMethod,
    WB,               // EvWBTemp,
    WB,               // EvWBGreen,
    AUTOEXP,         // EvToneCurveMode1,
    AUTOEXP,         // EvToneCurve2,
    AUTOEXP,         // EvToneCurveMode2,
    0,                // EvLDNRadius: obsolete,
    0,                // EvLDNEdgeTolerance: obsolete,
    0,                // EvCDNEnabled:obsolete,
    0,                // free entry
    ALLNORAW, //RGBCURVE | M_AUTOEXP, // EvDCPToneCurve, 21 july 2024
    ALLNORAW,         // EvDCPIlluminant,
    LUMINANCECURVE,          // EvSHEnabled,
    LUMINANCECURVE,         // EvSHHighlights,
    LUMINANCECURVE,         // EvSHShadows,
    LUMINANCECURVE,         // EvSHHLTonalW,
    LUMINANCECURVE,         // EvSHSHTonalW,
    0,         // EvSHLContrast : obsolete,
    LUMINANCECURVE,          // EvSHRadius,
    ALLNORAW,         // EvCTRotate,
    ALLNORAW,         // EvCTHFlip,
    ALLNORAW,         // EvCTVFlip,
    HDR,        // EvROTDegree,
    HDR,        // EvTransAutoFill,
    HDR,        // EvDISTAmount,
    ALL,              // EvBookmarkSelected,
    CROP,             // EvCrop,
    HDR,        // EvCACorr,
    ALLNORAW,         // EvHREnabled,
    0,         // EvHRAmount : obsolete,
    ALLNORAW|M_RAW,   // EvHRMethod,
    DEMOSAIC,         // EvWProfile,
    OUTPUTPROFILE,    // EvOProfile,
    ALLNORAW,         // EvIProfile,
    HDR,        // EvVignettingAmount,
    AUTOEXP,         // EvChMixer,
    RESIZE,           // EvResizeScale,
    RESIZE,           // EvResizeMethod,
    EXIF,             // EvExif,
    IPTC,             // EvIPTC
    0,           // EvResizeSpec : obsolete,
    RESIZE,           // EvResizeWidth
    RESIZE,           // EvResizeHeight
    RESIZE,           // EvResizeEnabled
    ALL,              // EvProfileChangeNotification
    0,          // EvSHHighQuality : obsolete
    HDR,        // EvPerspCorr
    DARKFRAME,        // EvLCPFile
    AUTOEXP,         // EvRGBrCurveLumamode
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
    AUTOEXP,         // EvSaturation,
    AUTOEXP,         // EvHSVEqualizerH,
    AUTOEXP,         // EvHSVEqualizerS,
    AUTOEXP,         // EvHSVEqualizerV,
    AUTOEXP,         // EvHSVEqEnabled,
    DEFRINGE,         // EvDefringeEnabled,
    DEFRINGE,         // EvDefringeRadius,
    DEFRINGE,         // EvDefringeThreshold,
    AUTOEXP,         // EvHLComprThreshold,
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
    0,                // --unused--
    FLATFIELD,        // EvFlatFieldFile,
    FLATFIELD,        // EvFlatFieldAutoSelect,
    FLATFIELD,        // EvFlatFieldBlurRadius,
    FLATFIELD,        // EvFlatFieldBlurType,
    HDR,        // EvAutoDIST,
    0,         // EvDPDNLumCurve : obsolete
    0,         // EvDPDNChromCurve : obsolete
    0,            // EvGAMMA : obsolete
    0,            // EvGAMPOS : obsolete
    0,            // EvGAMFREE : obsolete
    0,            // EvSLPOS : obsolete
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
    AUTOEXP,         // EvVibranceEnabled
    AUTOEXP,         // EvVibrancePastels
    AUTOEXP,         // EvVibranceSaturated
    AUTOEXP,         // EvVibranceProtectSkins
    AUTOEXP,         // EvVibranceAvoidColorShift
    AUTOEXP,         // EvVibrancePastSatTog
    AUTOEXP,         // EvVibrancePastSatThreshold
    SHARPENING,       // EvEPDStrength
    SHARPENING,       // EvEPDEdgeStopping
    SHARPENING,       // EvEPDScale
    SHARPENING,       // EvEPDReweightingIterates
    SHARPENING,       // EvEPDEnabled
    AUTOEXP,         // EvRGBrCurve
    AUTOEXP,         // EvRGBgCurve
    AUTOEXP,         // EvRGBbCurve
    AUTOEXP,         // EvNeutralExp
    DEMOSAIC | M_PREPROC, // EvDemosaicMethodPreProc
    LUMINANCECURVE,   // EvLCCurve
    LUMINANCECURVE,   // EvLCHCurve
    AUTOEXP,         // EvVibranceSkinTonesCurve
    LUMINANCECURVE,   // EvLLCCurve
    LUMINANCECURVE,   // EvLLCredsk
    ALLNORAW,         // EvDPDNLdetail
    LUMINANCECURVE,         // EvCATEnabled
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
    WB,               // EvWBequal
    0,         // EvWBequalbo : obsolete
    HDR,        // EvGradientDegree
    HDR,        // EvGradientEnabled
    HDR,        // EvPCVignetteStrength
    HDR,        // EvPCVignetteEnabled
    AUTOEXP,         // EvBWChmixEnabled
    AUTOEXP,         // EvBWred
    AUTOEXP,         // EvBWgreen
    AUTOEXP,         // EvBWblue
    AUTOEXP,         // EvBWredgam
    AUTOEXP,         // EvBWgreengam
    AUTOEXP,         // EvBWbluegam
    AUTOEXP,         // EvBWfilter
    AUTOEXP,         // EvBWsetting
    AUTOEXP,         // EvBWoran
    AUTOEXP,         // EvBWyell
    AUTOEXP,         // EvBWcyan
    AUTOEXP,         // EvBWmag
    AUTOEXP,         // EvBpur
    AUTOEXP,         // EvBWLuminanceEqual
    AUTOEXP,         // EvBWChmixEnabledLm
    AUTOEXP,         // EvBWmethod
    AUTOEXP,         // EvBWBeforeCurve
    AUTOEXP,         // EvBWBeforeCurveMode
    AUTOEXP,         // EvBWAfterCurve
    AUTOEXP,         // EvBWAfterCurveMode
    AUTOEXP,         // EvAutoch
    0,                // --unused--
    AUTOEXP,         // EvNeutralBW
    HDR,        // EvGradientFeather
    HDR,        // EvGradientStrength
    HDR,        // EvGradientCenter
    HDR,        // EvPCVignetteFeather
    HDR,        // EvPCVignetteRoundness
    HDR,        // EvVignettingRadius,
    HDR,        // EvVignettingStrength
    HDR,        // EvVignettingCenter
    LUMINANCECURVE,   // EvLCLCurve
    LUMINANCECURVE,   // EvLLHCurve
    LUMINANCECURVE,   // EvLHHCurve
    ALLNORAW,         // EvDirPyrEqualizerThreshold
    0,         // EvDPDNenhance : obsolete
    AUTOEXP,         // EvBWMethodalg
    ALLNORAW,         // EvDirPyrEqualizerSkin
    ALLNORAW,         // EvDirPyrEqlgamutlab
    ALLNORAW,         // EvDirPyrEqualizerHueskin
    ALLNORAW,         // EvDPDNmedian
    ALLNORAW,         // EvDPDNmedmet
    AUTOEXP,         // EvColorToningEnabled
    AUTOEXP,         // EvColorToningColor
    AUTOEXP,         // EvColorToningOpacity
    AUTOEXP,         // EvColorToningCLCurve
    AUTOEXP,         // EvColorToningMethod
    AUTOEXP,         // EvColorToningLLCurve
    AUTOEXP,         // EvColorToningredlow
    AUTOEXP,         // EvColorToninggreenlow
    AUTOEXP,         // EvColorToningbluelow
    AUTOEXP,         // EvColorToningredmed
    AUTOEXP,         // EvColorToninggreenmed
    AUTOEXP,         // EvColorToningbluemed
    AUTOEXP,         // EvColorToningredhigh
    AUTOEXP,         // EvColorToninggreenhigh
    AUTOEXP,         // EvColorToningbluehigh
    AUTOEXP,         // EvColorToningbalance
    AUTOEXP,         // EvColorToningNeutral
    0,         // EvColorToningsatlow : obsolete
    0,         // EvColorToningsathigh : obsolete
    AUTOEXP,         // EvColorToningTwocolor
    AUTOEXP,         // EvColorToningNeutralcur
    AUTOEXP,         // EvColorToningLumamode
    AUTOEXP,         // EvColorToningShadows
    AUTOEXP,         // EvColorToningHighights
    AUTOEXP,         // EvColorToningSatProtection
    AUTOEXP,         // EvColorToningSatThreshold
    AUTOEXP,         // EvColorToningStrength
    AUTOEXP,         // EvColorToningautosat
    ALLNORAW,         // EvDPDNmetmed
    ALLNORAW,         // EvDPDNrgbmet
    ALLNORAW,         // EvDPDNpasses
    FLATFIELD,        // EvFlatFieldClipControl
    FLATFIELD,        // EvFlatFieldAutoClipControl
    DARKFRAME,        // EvPreProcessExpBlackRed
    DARKFRAME,        // EvPreProcessExpBlackGreen
    DARKFRAME,        // EvPreProcessExpBlackBlue
    AUTOEXP,         // EvFilmSimulationEnabled
    AUTOEXP,         // EvFilmSimulationStrength
    AUTOEXP,         // EvFilmSimulationFilename
    ALLNORAW,         // EvDPDNLCurve
    ALLNORAW,         // EvDPDNsmet
    DARKFRAME,        // EvPreProcessDeadPixel
    ALLNORAW,         // EvDPDNCCCurve
    0,         // EvDPDNautochroma : obsolete
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
    ALLNORAW, //RGBCURVE | M_AUTOEXP, // EvDCPApplyLookTable,  21 july 2024
    ALLNORAW, //RGBCURVE | M_AUTOEXP, // EvDCPApplyBaselineExposureOffset, 21 july 2024
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
    0,          // EvLgain : obsolete
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
    WB,               // EvWBtempBias
    DARKFRAME,        // EvRawImageNum
    0,                // unused
    0,                // unused
    0,                // unused
    DEMOSAIC,         // EvPixelShiftEperIso
    0,                // unused
    0,                // unused
    DEMOSAIC,         // EvPixelshiftShowMotion
    DEMOSAIC,         // EvPixelshiftShowMotionMaskOnly
    0,                // unused
    0,                // unused
    0,                // unused
    DEMOSAIC,         // EvPixelShiftNonGreenCross
    0,                // unused
    0,                // unused
    0,                // unused
    0,                // unused
    DEMOSAIC,         // EvPixelShiftGreen
    0,                // unused
    DEMOSAIC,         // EvPixelShiftBlur
    DEMOSAIC,         // EvPixelShiftSigma
    0,                // unused
    0,                // unused
    DEMOSAIC,         // EvPixelShiftHoleFill
    DEMOSAIC,         // EvPixelShiftMedian
    0,                // unused
    DEMOSAIC,         // EvPixelShiftMotionMethod
    DEMOSAIC,         // EvPixelShiftSmooth
    0,         // EvPixelShiftLmmse : obsolete
    DEMOSAIC,         // EvPixelShiftEqualBright
    DEMOSAIC,         // EvPixelShiftEqualBrightChannel
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
    WB,               // EvWBEnabled
    AUTOEXP,         // EvRGBEnabled
    LUMINANCECURVE,   // EvLEnabled
    DEMOSAIC,         // EvPdShrEnabled
    CAPTURESHARPEN,    // EvPdShrMaskToggled
    AUTOEXP,   // EvLocallabSpotDeleted
    HDR,           // EvLocallabSpotSelected
    M_VOID,           // EvLocallabSpotName
    M_VOID,           // EvLocallabSpotVisibility
    AUTOEXP,   // EvLocallabSpotShape
    AUTOEXP,   // EvLocallabSpotSpotMethod
    AUTOEXP,   // EvLocallabSpotShapeMethod
    AUTOEXP,   // EvLocallabSpotLocX
    AUTOEXP,   // EvLocallabSpotLocXL
    AUTOEXP,   // EvLocallabSpotLocY
    AUTOEXP,   // EvLocallabSpotLocYT
    AUTOEXP,   // EvLocallabSpotCenter
    AUTOEXP,   // EvLocallabSpotCircrad
    AUTOEXP,   // EvLocallabSpotQualityMethod
    AUTOEXP,   // EvLocallabSpotTransit
    AUTOEXP,   // EvLocallabSpotThresh
    AUTOEXP,   // EvLocallabSpotIter
    AUTOEXP,   // EvLocallabSpotSensiexclu
    AUTOEXP,   // EvLocallabSpotStruc
    AUTOEXP,   // EvlocallabEnabled
    AUTOEXP,   // EvLocenacolor
    AUTOEXP,   // Evlocallabcurvactiv
    AUTOEXP,   // Evlocallablightness
    AUTOEXP,   // Evlocallabcontrast
    AUTOEXP,   // Evlocallabchroma
    AUTOEXP,   // Evlocallabsensi
    AUTOEXP,   // EvlocallabqualitycurveMethod
    AUTOEXP,   // Evlocallabllshape
    AUTOEXP,   // Evlocallabccshape
    AUTOEXP,   // EvlocallabLHshape
    AUTOEXP,   // EvlocallabHHshape
    AUTOEXP,   // Evlocallabinvers
    AUTOEXP,   // EvLocenaexpose
    AUTOEXP,   // Evlocallabexpcomp
    AUTOEXP,   // Evlocallabhlcompr
    AUTOEXP,   // Evlocallabhlcomprthresh
    AUTOEXP,   // Evlocallabblack
    AUTOEXP,   // Evlocallabshcompr
    AUTOEXP,   // Evlocallabwarm
    AUTOEXP,   // Evlocallabsensiex
    AUTOEXP,   // Evlocallabshapeexpos
    AUTOEXP,   // EvLocenavibrance
    AUTOEXP,   // EvlocallabSaturated
    AUTOEXP,   // EvlocallabPastels
    AUTOEXP,   // EvlocallabPastSatThreshold
    AUTOEXP,   // EvlocallabProtectSkins
    AUTOEXP,   // EvlocallabAvoidColorShift
    AUTOEXP,   // EvlocallabPastSatTog
    AUTOEXP,   // Evlocallabsensiv
    AUTOEXP,   // EvlocallabSkinTonesCurve
    AUTOEXP,   // EvLocenablur
    AUTOEXP,   // Evlocallabradius
    AUTOEXP,   // Evlocallabstrength
    AUTOEXP,   // Evlocallabsensibn
    AUTOEXP,   // EvlocallabblurMethod
    AUTOEXP,   // Evlocallabactivlum
    AUTOEXP,   // EvLocenatonemap
    AUTOEXP,   // Evlocallabstren
    AUTOEXP,   // Evlocallabgamma
    AUTOEXP,   // Evlocallabestop
    AUTOEXP,   // Evlocallabscaltm
    AUTOEXP,   // Evlocallabrewei
    AUTOEXP,   // Evlocallabsensitm
    AUTOEXP,   // EvLocenareti
    AUTOEXP,   // EvlocallabretinexMethod
    AUTOEXP,   // Evlocallabstr
    AUTOEXP,   // Evlocallabchrrt
    AUTOEXP,   // Evlocallabneigh
    AUTOEXP,   // Evlocallabvart
    AUTOEXP,   // Evlocallabsensih
    AUTOEXP,   // EvlocallabCTgainCurve
    AUTOEXP,   // Evlocallabinversret
    AUTOEXP,   // EvLocenasharp
    AUTOEXP,   // Evlocallabsharradius
    AUTOEXP,   // Evlocallabsharamount
    AUTOEXP,   // Evlocallabshardamping
    AUTOEXP,   // Evlocallabshariter
    AUTOEXP,   // Evlocallabsensis
    AUTOEXP,   // Evlocallabinverssha
    AUTOEXP,   // EvLocenacbdl
    AUTOEXP,   // EvlocallabEqualizer
    AUTOEXP,   // Evlocallabchromacbdl
    AUTOEXP,   // EvlocallabThresho
    AUTOEXP,   // Evlocallabsensicb
    AUTOEXP,   // EvLocenadenoi
    AUTOEXP,   // Evlocallabnoiselumf
    AUTOEXP,   // Evlocallabnoiselumc
    AUTOEXP,   // Evlocallabnoiselumdetail
    AUTOEXP,   // Evlocallabnoiselequal
    AUTOEXP,   // Evlocallabnoisechrof
    AUTOEXP,   // Evlocallabnoisechroc
    AUTOEXP,   // Evlocallabnoisechrodetail
    AUTOEXP,   // Evlocallabadjblur
    AUTOEXP,   // Evlocallabbilateral
    AUTOEXP,   // Evlocallabsensiden
    AUTOEXP,   // Evlocallabavoid
    AUTOEXP,   // Evlocallabsharcontrast
    AUTOEXP,   // EvLocenacontrast
    AUTOEXP,   // Evlocallablcradius
    AUTOEXP,   // Evlocallablcamount
    AUTOEXP,   // Evlocallablcdarkness
    AUTOEXP,   // Evlocallablclightness
    AUTOEXP,   // Evlocallabsensilc
    AUTOEXP,   // Evlocallabdehaz
    AUTOEXP,   // EvLocenasoft
    AUTOEXP,   // EvLocallabstreng
    AUTOEXP,   // EvLocallabsensisf
    AUTOEXP,   // Evlocallabsharblur
    0,   // EvLocenalabregion : obsolete
    AUTOEXP,   // EvlocallabshowmaskMethod
    AUTOEXP,   // EvLocallabSpotSelectedWithMask
    AUTOEXP,   // EvlocallabCCmaskshape
    AUTOEXP,   // EvlocallabLLmaskshape
    AUTOEXP,   // EvlocallabCCmaskexpshape
    AUTOEXP,   // EvlocallabLLmaskexpshape
    AUTOEXP,   // EvlocallabHHmaskshape
    AUTOEXP,   // Evlocallabstructcol
    AUTOEXP,   // Evlocallabstructexp
    AUTOEXP,   // EvlocallabHHmaskexpshape
    AUTOEXP,   // Evlocallabblendmaskcol
    AUTOEXP,   // Evlocallabblendmaskexp
    AUTOEXP,   // Evlocallabblurexpde
    AUTOEXP,   // EvLocallabEnaColorMask
    AUTOEXP,   // EvLocallabEnaExpMask
    AUTOEXP,   // Evlocallabblurcolde
    AUTOEXP,   // Evlocallabinversex
    AUTOEXP,   // Evlocallabstructexclu
    AUTOEXP,   // Evlocallabexpchroma
    AUTOEXP,   // EvLocallabLabGridValue
    AUTOEXP,   // EvLocallabLabstrengthgrid
    AUTOEXP,   // EvLocallabgridMethod
    AUTOEXP,   // EvLocenashadhigh
    AUTOEXP,   // EvLocallabhighlights
    AUTOEXP,   // EvLocallabh_tonalwidth
    AUTOEXP,   // EvLocallabshadows
    AUTOEXP,   // EvLocallabs_tonalwidth
    AUTOEXP,   // EvLocallabsh_radius
    AUTOEXP,   // EvLocallabsensihs
    AUTOEXP,   // Evlocallabradmaskcol
    AUTOEXP,   // Evlocallabradmaskexp
    AUTOEXP,   // EvlocallabToolAdded
    AUTOEXP,   // EvlocallabCCmaskSHshape
    AUTOEXP,   // EvlocallabLLmaskSHshape
    AUTOEXP,   // EvlocallabHHmaskSHshape
    AUTOEXP,   // EvlocallabblendmaskSH
    AUTOEXP,   // EvLocallabEnaSHMask
    AUTOEXP,   // EvlocallabradmaskSH
    AUTOEXP,   // EvlocallabblurSHde
    AUTOEXP,   // Evlocallabinverssh
    AUTOEXP,   // EvLocallabSpotbalan
    AUTOEXP,   // EvLocallabchromaskexp
    AUTOEXP,   // EvLocallabgammaskexp
    AUTOEXP,   // EvLocallabslomaskexp
    AUTOEXP,   // EvLocallabsoftradiusexp
    AUTOEXP,   // EvLocallabchromaskcol
    AUTOEXP,   // EvLocallabgammaskcol
    AUTOEXP,   // EvLocallabslomaskcol
    AUTOEXP,   // EvLocallabchromaskSH
    AUTOEXP,   // EvLocallabgammaskSH
    AUTOEXP,   // EvLocallabslomaskSH
    AUTOEXP,   // EvLocallabsoftradiuscol
    AUTOEXP,   // EvLocallabsoftradiusret
    AUTOEXP,   // EvLocallabsoftradiuscb
    AUTOEXP,   // EvLocallabSpotTransitweak
    AUTOEXP,   // EvLocallabclarityml
    AUTOEXP,   // EvLocallabcontresid
    AUTOEXP,   // Evlocallabnoiselumf0
    AUTOEXP,   // Evlocallabnoiselumf2
    0,   // Evlocallabblurcbdl
    AUTOEXP,   // Evlocallabblendmaskcb
    AUTOEXP,   // Evlocallabradmaskcb
    AUTOEXP,   // Evlocallabchromaskcb
    AUTOEXP,   // Evlocallabgammaskcb
    AUTOEXP,   // Evlocallabslomaskcb
    AUTOEXP,   // EvlocallabCCmaskcbshape
    AUTOEXP,   // EvlocallabLLmaskcbshape
    AUTOEXP,   // EvlocallabHHmaskcbshape
    AUTOEXP,   // EvLocallabEnacbMask
    M_VOID,           // EvlocallabToolRemovedWithoutRefresh
    AUTOEXP,   // Evlocallabsoftradiustm
    AUTOEXP,   // EvLocallabSpotTransitgrad
    AUTOEXP,   // Evlocallabamount
    AUTOEXP,   // Evlocallabsatur
    AUTOEXP,   // EvlocallabCCmaskretishape
    AUTOEXP,   // EvlocallabLLmaskretishape
    AUTOEXP,   // EvlocallabHHmaskretishape
    AUTOEXP,   // EvLocallabEnaretiMask
    AUTOEXP,   // Evlocallabblendmaskreti
    AUTOEXP,   // Evlocallabradmaskreti
    AUTOEXP,   // Evlocallabchromaskreti
    AUTOEXP,   // Evlocallabgammaskreti
    AUTOEXP,   // Evlocallabslomaskreti
    AUTOEXP,   // EvlocallabToolRemovedWithRefresh
    AUTOEXP,   // EvLocallabEnaretiMasktmap
    AUTOEXP,   // Evlocallabscalereti
    AUTOEXP,   // Evlocallabdarkness
    AUTOEXP,   // Evlocallablightnessreti
    AUTOEXP,   // Evlocallablimd
    AUTOEXP,   // Evlocallablaplace
    AUTOEXP,   // EvlocallabsoftMethod
    AUTOEXP,   // Evlocallabequilret
    AUTOEXP,   // Evlocallabequiltm
    AUTOEXP,   // Evlocallabfftwlc
    AUTOEXP,   // Evlocallabfftwreti
    AUTOEXP,   // EvlocallabshowmasksoftMethod
    AUTOEXP,   // Evlocallabshadex
    AUTOEXP,   // EvlocallabexpMethod
    AUTOEXP,   // EvLocallablaplacexp
    AUTOEXP,   // EvLocallabbalanexp
    AUTOEXP,   // EvLocallablinear
    AUTOEXP,   // EvlocallabCCmasktmshape
    AUTOEXP,   // EvlocallabLLmasktmshape
    AUTOEXP,   // EvlocallabHHmasktmshape
    AUTOEXP,   // EvLocallabEnatmMask
    AUTOEXP,   // Evlocallabblendmasktm
    AUTOEXP,   // Evlocallabradmasktm
    AUTOEXP,   // Evlocallabchromasktm
    AUTOEXP,   // Evlocallabgammasktm
    AUTOEXP,   // Evlocallabslomasktm
    AUTOEXP,   // EvlocallabshowmasktmMethod
    AUTOEXP,   // EvlocallablocalcontMethod
    AUTOEXP,   // Evlocallabwavcurve
    AUTOEXP,   // Evlocallablevelwav
    AUTOEXP,   // Evlocallabresidcont
    AUTOEXP,   // EvlocallabCCmaskblshape
    AUTOEXP,   // EvlocallabLLmaskblshape
    AUTOEXP,   // EvlocallabHHmaskblshape
    AUTOEXP,   // EvLocallabEnablMask
    AUTOEXP,   // EvlocallabshowmaskblMethod
    AUTOEXP,   // Evlocallabblendmaskbl
    AUTOEXP,   // Evlocallabradmaskbl
    AUTOEXP,   // Evlocallabchromaskbl
    AUTOEXP,   // Evlocallabgammaskbl
    AUTOEXP,   // Evlocallabslomaskbl
    AUTOEXP,   // EvlocallabblMethod
    AUTOEXP,   // EvlocallabmedMethod
    AUTOEXP,   // Evlocallabitera
    AUTOEXP,   // Evlocallabguidbl
    AUTOEXP,   // Evlocallabepsbl
    AUTOEXP,   // EvlocallabshowmaskcolMethodinv
    AUTOEXP,   // EvlocallabshowmaskexpMethodinv
    AUTOEXP,   // EvlocallabshowmaskSHMethodinv
    AUTOEXP,   // Evlocallabclarilres
    AUTOEXP,   // Evlocallabclarisoft
    AUTOEXP,   // Evlocallabclaricres
    AUTOEXP,   // Evlocallabresidchro
    AUTOEXP,   // Evlocallabgamm
    AUTOEXP,   // Evlocallabfatamount
    AUTOEXP,   // Evlocallabfatdetail
    AUTOEXP,   // Evlocallabfatanchor
    AUTOEXP,   // Evlocallabfatlevel
    AUTOEXP,   // EvlocallabSpotCreated
    AUTOEXP,   // EvlocallabexnoiseMethod
    AUTOEXP,   // Evlocallabdepth
    AUTOEXP,   // Evlocallabloglin
    AUTOEXP,   // EvlocallabdehazeSaturation
    AUTOEXP,   // Evlocallaboffs
    AUTOEXP,   // EvlocallabCTtransCurve
    AUTOEXP,   // Evlocallabcliptm
    AUTOEXP,   // Evlocallabenatmmaskaft
    AUTOEXP,   // EvlocallabenaExpmaskaft
    AUTOEXP,   // Evlocallablapmasktm
    AUTOEXP,   // Evlocallablapmaskreti
    AUTOEXP,   // Evlocallablapmaskexp
    AUTOEXP,   // Evlocallablapmaskcol
    AUTOEXP,   // EvlocallablapmaskSH
    AUTOEXP,   // Evlocallablapmaskcb
    AUTOEXP,   // Evlocallablapmaskbl
    AUTOEXP,   // Evlocallablaplac
    AUTOEXP,   // Evlocallabdetailthr
    AUTOEXP,   // Evlocallabfftwbl
    AUTOEXP,   // Evlocallabisogr
    AUTOEXP,   // Evlocallabstrengr
    AUTOEXP,   // Evlocallabscalegr
    AUTOEXP,   // EvlocallabLmaskshape
    AUTOEXP,   // EvlocallabLmaskexpshape
    AUTOEXP,   // EvlocallabLmaskSHshape
    AUTOEXP,   // EvlocallabLmasktmshape
    AUTOEXP,   // EvlocallabLmaskretishape
    AUTOEXP,   // EvlocallabLmaskcbshape
    AUTOEXP,   // EvlocallabLmaskblshape
    AUTOEXP,   // EvlocallabLLmaskblshapewav
    AUTOEXP,   // Evlocallabshadmaskbl
    AUTOEXP,   // EvlocallabLLmaskcolshapewav
    AUTOEXP,   // Evlocallabshadmaskcol
    AUTOEXP,   // EvlocallabcsThreshold
    AUTOEXP,   // EvlocallabcsThresholdblur
    AUTOEXP,   // EvlocallabcsThresholdcol
    AUTOEXP,   // Evlocallabdeltae
    AUTOEXP,   // EvLocallabSpotscopemask
    AUTOEXP,   // EvlocallabshMethod
    AUTOEXP,   // EvlocallabEqualizersh
    AUTOEXP,   // EvlocallabdetailSH
    AUTOEXP,   // EvlocallabfatamountSH
    AUTOEXP,   // EvlocallabfatanchorSH
    AUTOEXP,   // Evlocallabshortc
    AUTOEXP,   // EvLocallabSpotlumask
    AUTOEXP,   // EvlocallabgamSH
    AUTOEXP,   // EvlocallabsloSH
    AUTOEXP,   // Evlocallabsavrest
    AUTOEXP,   // Evlocallabrecurs
    AUTOEXP,   // EvLocallabmergecolMethod
    AUTOEXP,   // EvLocallabopacol
    AUTOEXP,   // Evlocallabrgbshape
    AUTOEXP,   // EvLocallabtoneMethod
    AUTOEXP,   // EvLocallabspecial
    AUTOEXP,   // EvLocallabconthrcol
    AUTOEXP,   // EvLocallabmerMethod
    AUTOEXP,   // EvLocallabstrumaskcol
    AUTOEXP,   // EvLocallabstrumaskbl
    AUTOEXP,   // EvLocallabtoolcol
    AUTOEXP,   // Evlocallabtoolbl
    AUTOEXP,   // EvlocallabHHhmaskshape
    AUTOEXP,   // EvlocallabCCmaskvibshape
    AUTOEXP,   // EvlocallabLLmaskvibshape
    AUTOEXP,   // EvlocallabHHmaskvibshape
    AUTOEXP,   // EvlocallabshowmaskvibMethod
    AUTOEXP,   // EvLocallabEnavibMask
    AUTOEXP,   // Evlocallabblendmaskvi
    AUTOEXP,   // Evlocallabradmaskvib
    AUTOEXP,   // Evlocallabchromaskvib
    AUTOEXP,   // Evlocallabgammaskvib
    AUTOEXP,   // Evlocallabslomaskvib
    AUTOEXP,   // Evlocallablapmaskvib
    AUTOEXP,   // EvlocallabLmaskvibshape
    AUTOEXP,   // EvLocallabLabGridmergValue
    AUTOEXP,   // EvLocallabmercol
    AUTOEXP,   // EvLocallabmerlucol
    AUTOEXP,   // Evlocallabstrmaskexp
    AUTOEXP,   // Evlocallabangmaskexp
    AUTOEXP,   // Evlocallabstrexp
    AUTOEXP,   // Evlocallabangexp
    AUTOEXP | M_AUTOEXP,   // EvlocallabstrSH
    AUTOEXP,   // EvlocallabangSH
    AUTOEXP,   // Evlocallabstrcol
    AUTOEXP,   // Evlocallabangcol
    AUTOEXP,   // Evlocallabstrcolab
    AUTOEXP,   // EvLocallabSpotfeather
    AUTOEXP,   // Evlocallabstrcolh
    AUTOEXP,   // Evlocallabstrvib
    AUTOEXP,   // Evlocallabangvib
    AUTOEXP,   // Evlocallabstrvibab
    AUTOEXP,   // Evlocallabstrvibh
    AUTOEXP,   // EvLocallabSpotcomplexMethod
    AUTOEXP,   // Evlocallabclshape
    AUTOEXP,   // Evlocallablcshape
    AUTOEXP,   // Evlocallabblurcol
    AUTOEXP,   // Evlocallabcontcol
    AUTOEXP,   // EvLocallabfftColorMask
    AUTOEXP | M_AUTOEXP,   // EvLocenalog
    HDR,          // EvLocallabAuto
    AUTOEXP,   // EvlocallabsourceGray
    0,          // EvlocallabsourceGrayAuto : obsolete
    HDR,          // EvlocallabAutoGray
    AUTOEXP,   // EvlocallabblackEv
    AUTOEXP,   // EvlocallabwhiteEv
    AUTOEXP,   // EvlocallabtargetGray
    AUTOEXP,   // Evlocallabdetail
    AUTOEXP,   // Evlocallabsensilog
    AUTOEXP,          // Evlocallabfullimage
    AUTOEXP,   // Evlocallabbaselog
    AUTOEXP,   // Evlocallabresidblur
    AUTOEXP,   // Evlocallabblurlc
    AUTOEXP,   // Evlocallablevelblur
    AUTOEXP,   // EvlocallabwavCurvelev
    AUTOEXP,   // EvlocallabwavCurvecon
    AUTOEXP,   // Evlocallabsigma
    AUTOEXP,   // Evlocallaboriglc
    AUTOEXP,   // Evlocallabsigmadc
    AUTOEXP,   // Evlocallabdeltad
    AUTOEXP,   // EvlocallabwavCurvecomp
    AUTOEXP,   // Evlocallabfatres
    AUTOEXP,   // EvLocallabSpotbalanh
    AUTOEXP,   // EvlocallabwavCurveden
    AUTOEXP,   // EvlocallabHHmasklcshape
    AUTOEXP,   // EvlocallabCCmasklcshape
    AUTOEXP,   // EvlocallabLLmasklcshape
    AUTOEXP,   // EvlocallabEnalcMask
    AUTOEXP,   // EvlocallabshowmasklcMethod
    AUTOEXP,   // Evlocallabblendmasklc
    AUTOEXP,   // Evlocallabradmasklc
    AUTOEXP,   // Evlocallabchromasklc
    AUTOEXP,   // EvlocallabLmasklcshape
    AUTOEXP,   // Evlocallabchromalev
    AUTOEXP,   // Evlocallabchromablu
    AUTOEXP,   // Evlocallaboffset
    AUTOEXP,   // Evlocallabwavblur
    AUTOEXP,   // Evlocallabwavcont
    AUTOEXP,   // Evlocallabwavcomp
    AUTOEXP,   // Evlocallabwavcompre
    AUTOEXP,   // EvlocallabwavCurvecompre
    AUTOEXP,   // Evlocallabresidcomp
    AUTOEXP,   // Evlocallabthreswav
    AUTOEXP,   // Evlocallabstrwav
    AUTOEXP,   // Evlocallabangwav
    AUTOEXP,   // Evlocallabwavgradl
    AUTOEXP,   // Evlocallabstrlog
    AUTOEXP,   // Evlocallabanglog
    AUTOEXP,   // EvLocallabSpotcolorde
    AUTOEXP,   // EvlocallabshowmasksharMethod
    AUTOEXP,   // Evlocallabshowreset
    AUTOEXP,   // Evlocallabstrengthw
    AUTOEXP,   // Evlocallabradiusw
    AUTOEXP,   // Evlocallabdetailw
    AUTOEXP,   // Evlocallabgradw
    AUTOEXP,   // Evlocallabtloww
    AUTOEXP,   // Evlocallabthigw
    AUTOEXP,   // EvlocallabwavCurveedg
    AUTOEXP,   // EvlocallablocaledgMethod
    AUTOEXP,   // Evlocallabwavedg
    AUTOEXP,   // Evlocallabedgw
    AUTOEXP,   // Evlocallabbasew
    AUTOEXP,   // EvlocallablocalneiMethod
    AUTOEXP,   // Evlocallabwaveshow
    AUTOEXP,   // EvLocallabSpotwavMethod
    AUTOEXP,   // EvlocallabchroMethod
    AUTOEXP,   // Evlocallabstrbl
    AUTOEXP,   // Evlocallabsigmadr
    AUTOEXP,   // Evlocallabsigmabl
    AUTOEXP,   // Evlocallabsigmaed
    AUTOEXP,   // Evlocallabresidsha
    AUTOEXP,   // Evlocallabresidshathr
    AUTOEXP,   // Evlocallabresidhi
    AUTOEXP,   // Evlocallabresidhithr
    AUTOEXP,   // Evlocallabsigmalc
    AUTOEXP,   // Evlocallabsigmalc2
    AUTOEXP,   // Evlocallabblwh
    AUTOEXP,   // EvlocallabcomplexityWithRefresh
    0,                // can be reused
    AUTOEXP,   // EvLocallabSpotcolorscope
    AUTOEXP,   // EvlocallabshowmasktypMethod
    AUTOEXP,   // Evlocallabshadmaskblsha
    AUTOEXP,   // EvLocenamask
    AUTOEXP,   // Evlocallabsensimask
    AUTOEXP,   // Evlocallabblendmask
    AUTOEXP,   // EvLocallabEna_Mask
    AUTOEXP,   // Evlocallabradmask
    AUTOEXP,   // Evlocallablapmask
    AUTOEXP,   // Evlocallabchromask
    AUTOEXP,   // Evlocallabgammask
    AUTOEXP,   // Evlocallabslopmask
    AUTOEXP,   // EvlocallabCCmask_shape
    AUTOEXP,   // EvlocallabLLmask_shape
    AUTOEXP,   // EvlocallabHHmask_shape
    AUTOEXP,   // EvLocallabtoolmask
    AUTOEXP,   // Evlocallabstrumaskmask
    AUTOEXP,   // EvlocallabHHhmask_shape
    AUTOEXP,   // EvLocallabfftmask
    AUTOEXP,   // Evlocallabblurmask
    AUTOEXP,   // Evlocallabcontmask
    AUTOEXP,   // Evlocallabshadmask
    AUTOEXP,   // EvlocallabLmask_shape
    AUTOEXP,   // EvlocallabLLmask_shapewav
    AUTOEXP,   // EvlocallabcsThresholdmask
    AUTOEXP,   // Evlocallabstr_mask
    AUTOEXP,   // Evlocallabang_mask
    AUTOEXP,   // Evlocallabsoftradiusmask
    AUTOEXP,   // Evlocallabblendmaskab
    AUTOEXP,   // EvLocallabSpotprevMethod
    AUTOEXP,   // Evlocallabactiv
    AUTOEXP,   // EvlocallabCHshape
    AUTOEXP,   //EvlocallabquaMethod
    AUTOEXP,   //Evlocallabhishow
    AUTOEXP,   // Evlocallabinvbl
    AUTOEXP,   // Evlocallabcatad
    AUTOEXP,   // Evlocallabciecam
    AUTOEXP,   // Evlocallabsourceabs
    AUTOEXP,   // Evlocallabtargabs
    AUTOEXP,   // Evlocallabsurround
    AUTOEXP,   // Evlocallabsaturl
    AUTOEXP,   // Evlocallabcontl
    AUTOEXP,   //EvlocallabCCmaskshapeL 
    AUTOEXP,   //EvlocallabLLmaskshapeL
    AUTOEXP,   // EvlocallabHHmaskshapeL
    AUTOEXP,   // EvlocallabenaLMask
    AUTOEXP,   // EvlocallabblendmaskL
    AUTOEXP,   // EvlocallabradmaskL
    AUTOEXP,   // EvlocallabchromaskL
    AUTOEXP,   //EvlocallabLmaskshapeL
    AUTOEXP,   // Evlocallablightl
    AUTOEXP,   // EvlocallabLshapeL
    AUTOEXP,   // Evlocallabcontq
    AUTOEXP,   // Evlocallabsursour
    AUTOEXP,   // Evlocallablightq
    AUTOEXP,   // Evlocallabcolorfl
    AUTOEXP,   // Evlocallabrepar
    AUTOEXP,   //EvlocallabwavCurvehue
    AUTOEXP,   // Evlocallablevelthr
    AUTOEXP,   // Evlocallablevelthrlow
    AUTOEXP,   //Evlocallabusemask1
    AUTOEXP,   // Evlocallablnoiselow
    AUTOEXP,   // Evlocallabrecothres
    AUTOEXP,   // Evlocallablowthres
    AUTOEXP,   // Evlocallabhigthres
    AUTOEXP,   // Evlocallabrecothresd
    AUTOEXP,   // Evlocallablowthresd
    AUTOEXP,   // Evlocallabhigthresd
    AUTOEXP,   // Evlocallabinvmaskd
    AUTOEXP,   // Evlocallabinvmask
    AUTOEXP,   // Evlocallabdecayd
    AUTOEXP,   // Evlocallabrecothresc
    AUTOEXP,   // Evlocallablowthresc
    AUTOEXP,   // Evlocallabhigthresc
    AUTOEXP,   // Evlocallabdecayc
    AUTOEXP,   // Evlocallabmidthresd
    AUTOEXP,   // Evlocallabrecothresl
    AUTOEXP,   // Evlocallablowthresl
    AUTOEXP,   // Evlocallabhigthresl
    AUTOEXP,  // Evlocallabdecayl
    AUTOEXP,   // Evlocallabrecothrese
    AUTOEXP,   // Evlocallablowthrese
    AUTOEXP,   // Evlocallabhigthrese
    AUTOEXP,   // Evlocallabdecaye
    AUTOEXP,   // Evlocallabrecothress
    AUTOEXP,   // Evlocallablowthress
    AUTOEXP,   // Evlocallabhigthress
    AUTOEXP,   // Evlocallabdecays
    AUTOEXP,   // Evlocallabrecothrev
    AUTOEXP,   // Evlocallablowthresv
    AUTOEXP,   // Evlocallabhigthresv
    AUTOEXP,   // Evlocallabdecayv
    AUTOEXP,   // Evlocallabrecothrew
    AUTOEXP,   // Evlocallablowthresw
    AUTOEXP,   // Evlocallabhigthresw
    AUTOEXP,   // Evlocallabdecayw
    AUTOEXP,   // Evlocallabmidthresdch
    AUTOEXP,   // Evlocallabrecothret
    AUTOEXP,   // Evlocallablowthrest
    AUTOEXP,   // Evlocallabhigthrest
    AUTOEXP,   // Evlocallabdecayt
    AUTOEXP,   // Evlocallabrecothrecb
    AUTOEXP,   // Evlocallablowthrescb
    AUTOEXP,   // Evlocallabhigthrescb
    AUTOEXP,   // Evlocallabdecaycb
    AUTOEXP,   // Evlocallabrecothrer
    AUTOEXP,   // Evlocallablowthresr
    AUTOEXP,   // Evlocallabhigthresr
    AUTOEXP,   // Evlocallabdecayr
    AUTOEXP,   // Evlocallabnlstr
    AUTOEXP,   // Evlocallabnldet
    AUTOEXP,   // Evlocallabnlpat
    AUTOEXP,   // Evlocallabnlrad
    AUTOEXP,   // Evlocallabnlgam
    AUTOEXP,   // Evlocallabdivgr
    AUTOEXP,   // EvLocallabSpotavoidrad
    AUTOEXP,   // EvLocallabSpotavoidmun
    AUTOEXP,   // Evlocallabcontthres
    AUTOEXP,   // Evlocallabnorm
    AUTOEXP,   // Evlocallabreparw
    AUTOEXP,   // Evlocallabreparcol
    AUTOEXP,   // Evlocallabreparden
    AUTOEXP,   // Evlocallabreparsh
    AUTOEXP,   // Evlocallabreparexp
    AUTOEXP,   // Evlocallabrepartm
    AUTOEXP,   // Evlocallabchroml
    AUTOEXP,   // Evlocallabresidgam
    AUTOEXP,   // Evlocallabresidslop
    AUTOEXP,   // Evlocallabnoisegam
    AUTOEXP,    //Evlocallabgamlc
    AUTOEXP,    //Evlocallabgamc
    AUTOEXP,    //Evlocallabgamex
    AUTOEXP | M_AUTOEXP,    // EvLocenacie
    AUTOEXP,     //Evlocallabreparcie
    HDR,     //EvlocallabAutograycie
    HDR,    //EvlocallabsourceGraycie
    HDR,    //Evlocallabsourceabscie
    AUTOEXP,    //Evlocallabsursourcie
    AUTOEXP,    //Evlocallabsaturlcie
    AUTOEXP,    //Evlocallabchromlcie
    AUTOEXP,    //Evlocallablightlcie
    AUTOEXP,    //Evlocallablightqcie
    AUTOEXP,    //Evlocallabcontlcie
    AUTOEXP,    //Evlocallabcontthrescie
    AUTOEXP,    //Evlocallabcontqcie
    AUTOEXP,    //Evlocallabcolorflcie
    AUTOEXP,    //Evlocallabtargabscie
    AUTOEXP,    //EvlocallabtargetGraycie
    AUTOEXP,    //Evlocallabcatadcie
    AUTOEXP,    //Evlocallabdetailcie
    AUTOEXP,    //Evlocallabsurroundcie
    AUTOEXP,    //Evlocallabsensicie
    AUTOEXP,    //Evlocallabmodecie
    AUTOEXP,    //Evlocallabrstprotectcie
    AUTOEXP,    //Evlocallabsigmoidldacie12
    AUTOEXP,    //Evlocallabsigmoidthcie12
    AUTOEXP,    //Evlocallabsigmoidblcie12
    HDR,    //Evlocallabcomprcieauto
    AUTOEXP,    //Evlocallabhuecie
    AUTOEXP,    //Evlocallabjabcie
    AUTOEXP,    //Evlocallablightjzcie
    AUTOEXP,    //Evlocallabcontjzcie
    AUTOEXP,    //Evlocallabchromjzcie
    AUTOEXP,    //Evlocallabhuejzcie
    AUTOEXP,    //Evlocallabsigmoidldajzcie12
    AUTOEXP,    //Evlocallabsigmoidthjzcie12
    AUTOEXP,    //Evlocallabsigmoidbljzcie12
    AUTOEXP,    //Evlocallabadapjzcie
    AUTOEXP,    //Evlocallabmodecam
    AUTOEXP,    //Evlocallabhljzcie
    AUTOEXP,    //Evlocallabhlthjzcie
    AUTOEXP,    //Evlocallabshjzcie
    AUTOEXP,    //Evlocallabshthjzcie
    AUTOEXP,    //Evlocallabradjzcie
//    AUTOEXP,    //EvlocallabHHshapejz
    AUTOEXP,    //EvlocallabCHshapejz
    AUTOEXP,    //Evlocallabjz100
    AUTOEXP,    //Evlocallabpqremap
    AUTOEXP,    //EvlocallabLHshapejz
    AUTOEXP,    //Evlocallabshargam
    AUTOEXP,    //Evlocallabvibgam
    AUTOEXP,    //EvLocallabtoneMethodcie
    AUTOEXP,    //Evlocallabshapecie
    AUTOEXP,    //EvLocallabtoneMethodcie2
    AUTOEXP,    //Evlocallabshapecie2
    AUTOEXP,    //Evlocallabshapejz
    AUTOEXP,    //Evlocallabshapecz
    AUTOEXP,    //Evlocallabshapeczjz
//    AUTOEXP,    //Evlocallabforcejz
//    AUTOEXP,    //Evlocallablightlzcam
//    AUTOEXP,    //Evlocallablightqzcam
//    AUTOEXP,    //Evlocallabcontlzcam
//    AUTOEXP,    //Evlocallabcontqzcam
//    AUTOEXP,    //Evlocallabcontthreszcam
//    AUTOEXP,    //Evlocallabcolorflzcam
//    AUTOEXP,    //Evlocallabsaturzcam
//    AUTOEXP,    //Evlocallabchromzcam
    AUTOEXP,    //Evlocallabpqremapcam16
    AUTOEXP,    //EvLocallabEnacieMask    
    AUTOEXP,    //EvlocallabCCmaskcieshape
    AUTOEXP,    //EvlocallabLLmaskcieshape
    AUTOEXP,    //EvlocallabHHmaskcieshape
    AUTOEXP,    //Evlocallabblendmaskcie
    AUTOEXP,    //Evlocallabradmaskcie
    AUTOEXP,    //Evlocallabchromaskcie
    AUTOEXP,    //EvlocallabLmaskcieshape 
    AUTOEXP,    //Evlocallabrecothrescie
    AUTOEXP,    //Evlocallablowthrescie
    AUTOEXP,    //Evlocallabhigthrescie
    AUTOEXP,    //Evlocallabdecaycie
    AUTOEXP,    //Evlocallablapmaskcie
    AUTOEXP,    //Evlocallabgammaskcie
    AUTOEXP,    //Evlocallabslomaskcie
    AUTOEXP,    //Evlocallabqtoj
    AUTOEXP,    //Evlocallabsaturjzcie
    AUTOEXP,    //EvLocallabSpotdenoichmask
    AUTOEXP,    //Evlocallabsigmalcjz
    AUTOEXP,    //EvlocallabcsThresholdjz
    AUTOEXP,    //EvlocallabwavCurvejz
    AUTOEXP,    //Evlocallabclarilresjz
    AUTOEXP,    //Evlocallabclaricresjz
    AUTOEXP,    //Evlocallabclarisoftjz
    AUTOEXP,    //EvlocallabHHshapejz
    AUTOEXP,    //Evlocallabsoftjzcie
    AUTOEXP,    //Evlocallabthrhjzcie
    AUTOEXP,    //Evlocallabchjzcie
    AUTOEXP,    //Evlocallabstrsoftjzcie
    AUTOEXP,    //EvlocallabblackEvjz
    AUTOEXP,    //EvlocallabwhiteEvjz
    AUTOEXP,    //Evlocallablogjz
    AUTOEXP,    //Evlocallabtargetjz
    AUTOEXP,    //Evlocallabsigybjz12
    AUTOEXP,    //Evlocallabsigjz12
    AUTOEXP,    //Evlocallabsigq12
    AUTOEXP     //Evlocallablogcie
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
