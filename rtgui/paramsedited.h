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
#pragma once

#include <vector>

namespace rtengine
{

namespace procparams
{

class ProcParams;

class PartialProfile;

}

}

struct GeneralParamsEdited {
    bool rank;
    bool colorlabel;
    bool intrash;
};

struct ToneCurveParamsEdited {
    bool curve;
    bool curve2;
    bool curveMode;
    bool curveMode2;
    bool brightness;
    bool black;
    bool contrast;
    bool saturation;
    bool shcompr;
    bool hlcompr;
    bool hlcomprthresh;
    bool autoexp;
    bool clip;
    bool expcomp;
    bool hrenabled;
    bool method;
    bool histmatching;
    bool fromHistMatching;
    bool clampOOG;
};

struct RetinexParamsEdited {
    bool enabled;
    bool str;
    bool scal;
    bool iter;
    bool grad;
    bool grads;
    bool gam;
    bool slope;
    bool neigh;
    bool offs;
    bool retinexMethod;
    bool mapMethod;
    bool viewMethod;
    bool retinexcolorspace;
    bool gammaretinex;
    bool vart;
    bool limd;
    bool highl;
    bool baselog;
    bool skal;
    bool method;
    bool transmissionCurve;
    bool gaintransmissionCurve;
    bool cdcurve;
    bool mapcurve;
    bool cdHcurve;
    bool lhcurve;
    bool retinex;
    bool medianmap;
    bool highlights;
    bool htonalwidth;
    bool shadows;
    bool stonalwidth;
    bool radius;

    bool isUnchanged() const;
};


struct LCurveParamsEdited {
    bool enabled;
    bool brightness;
    bool contrast;
    bool chromaticity;
    bool avoidcolorshift;
    bool rstprotection;
    bool lcurve;
    bool acurve;
    bool bcurve;
    bool lcredsk;
    bool cccurve;
    bool chcurve;
    bool lhcurve;
    bool hhcurve;
    bool lccurve;
    bool clcurve;
};


struct LocalContrastParamsEdited {
    bool enabled;
    bool radius;
    bool amount;
    bool darkness;
    bool lightness;
};

struct RGBCurvesParamsEdited {
    bool enabled;
    bool lumamode;
    bool rcurve;
    bool gcurve;
    bool bcurve;
};

struct ColorToningEdited {
    bool enabled;
    bool opacityCurve;
    bool colorCurve;
    bool clcurve;
    bool method;
    bool autosat;
    bool satprotectionthreshold;
    bool saturatedopacity;
    bool strength;
    bool shadowsColSat;
    bool hlColSat;
    bool balance;
    bool twocolor;
    bool cl2curve;
    bool redlow;
    bool greenlow;
    bool bluelow;
    bool redmed;
    bool greenmed;
    bool bluemed;
    bool redhigh;
    bool greenhigh;
    bool bluehigh;
    bool satlow;
    bool sathigh;
    bool lumamode;
    bool labgridALow;
    bool labgridBLow;
    bool labgridAHigh;
    bool labgridBHigh;
    bool labregions;
    bool labregionsShowMask;
};

struct SharpenEdgeParamsEdited {
    bool enabled;
    bool passes;
    bool amount;
    bool threechannels;
};

struct SharpenMicroParamsEdited {
    bool enabled;
    bool matrix;
    bool amount;
    bool contrast;
    bool uniformity;
};

struct SharpeningParamsEdited {
    bool enabled;
    bool contrast;
    bool autoContrast;
    bool blurradius;
    bool gamma;
    bool radius;
    bool amount;
    bool threshold;
    bool edgesonly;
    bool edges_radius;
    bool edges_tolerance;
    bool halocontrol;
    bool halocontrol_amount;

    bool method;
    bool deconvamount;
    bool deconvradius;
    bool deconviter;
    bool deconvdamping;
};

struct CaptureSharpeningParamsEdited {
    bool enabled;
    bool contrast;
    bool autoContrast;
    bool autoRadius;
    bool deconvradius;
    bool deconvradiusOffset;
    bool deconviter;
    bool deconvitercheck;
    bool isUnchanged() const;
};

struct VibranceParamsEdited {
    bool enabled;
    bool pastels;
    bool saturated;
    bool psthreshold;
    bool protectskins;
    bool avoidcolorshift;
    bool pastsattog;
    bool skintonescurve;
};

struct WBParamsEdited {
    bool enabled;
    bool method;
    bool temperature;
    bool green;
    bool equal;
    bool tempBias;
};

struct DefringeParamsEdited {
    bool enabled;
    bool radius;
    bool threshold;
    bool huecurve;
};

struct ImpulseDenoiseParamsEdited {
    bool enabled;
    bool thresh;
};

struct ColorAppearanceParamsEdited {
    bool curve;
    bool curve2;
    bool curve3;
    bool curveMode;
    bool curveMode2;
    bool curveMode3;
    bool enabled;
    bool degree;
    bool autodegree;
    bool degreeout;
    bool autodegreeout;
    bool autoadapscen;
    bool autoybscen;
    bool surround;
    bool surrsrc;
    bool adapscen;
    bool adaplum;
    bool ybscen;
    bool badpixsl;
    bool wbmodel;
    bool algo;
    bool jlight;
    bool qbright;
    bool chroma;
    bool schroma;
    bool mchroma;
    bool contrast;
    bool qcontrast;
    bool colorh;
    bool rstprotection;
    bool surrsource;
    bool gamut;
    bool datacie;
    bool tonecie;
    bool tempout;
    bool greenout;
    bool ybout;
    bool tempsc;
    bool greensc;
};

struct DirPyrDenoiseParamsEdited {
    bool enabled;
    bool enhance;
    bool median;
    bool Ldetail;
    bool luma;
    bool chroma;
    bool redchro;
    bool bluechro;
    bool gamma;
    bool lcurve;
    bool cccurve;

    bool dmethod;
    bool Lmethod;
    bool Cmethod;
    bool C2method;
    bool smethod;
    bool medmethod;
    bool methodmed;
    bool rgbmethod;
    bool passes;
};

struct EPDParamsEdited {
    bool enabled;
    bool strength;
    bool gamma;
    bool edgeStopping;
    bool scale;
    bool reweightingIterates;
};

struct FattalToneMappingParamsEdited {
    bool enabled;
    bool threshold;
    bool amount;
    bool anchor;
};

struct SHParamsEdited {
    bool enabled;
    bool highlights;
    bool htonalwidth;
    bool shadows;
    bool stonalwidth;
    bool radius;
    bool lab;
};

struct CropParamsEdited {
    bool enabled;
    bool x;
    bool y;
    bool w;
    bool h;
    bool fixratio;
    bool ratio;
    bool orientation;
    bool guide;
};

struct CoarseTransformParamsEdited {
    bool rotate;
    bool hflip;
    bool vflip;
};

struct CommonTransformParamsEdited {
    bool method;
    bool autofill;
};

struct RotateParamsEdited {
    bool degree;
};

struct DistortionParamsEdited {
    bool amount;
};

struct LensProfParamsEdited {
    bool lcpFile;
    bool useDist;
    bool useVign;
    bool useCA;

    bool useLensfun;
    bool lfAutoMatch;
    bool lfCameraMake;
    bool lfCameraModel;
    bool lfLens;

    bool lcMode;

    bool isUnchanged() const;
};

struct PerspectiveParamsEdited {
    bool horizontal;
    bool vertical;
};

struct GradientParamsEdited {
    bool enabled;
    bool degree;
    bool feather;
    bool strength;
    bool centerX;
    bool centerY;
};

struct PCVignetteParamsEdited {
    bool enabled;
    bool strength;
    bool feather;
    bool roundness;
};

struct VignettingParamsEdited {
    bool amount;
    bool radius;
    bool strength;
    bool centerX;
    bool centerY;
};

struct ChannelMixerParamsEdited {
    bool enabled;
    bool red[3];
    bool green[3];
    bool blue[3];

};

struct BlackWhiteParamsEdited {
    bool enabledcc;
    bool enabled;
    bool method;
    bool filter;
    bool setting;
    bool mixerRed;
    bool mixerOrange;
    bool mixerYellow;
    bool mixerGreen;
    bool mixerCyan;
    bool mixerBlue;
    bool mixerMagenta;
    bool mixerPurple;
    bool gammaRed;
    bool gammaGreen;
    bool gammaBlue;
    bool luminanceCurve;
    bool beforeCurve;
    bool beforeCurveMode;
    bool afterCurve;
    bool afterCurveMode;
    bool autoc;
    bool algo;
};

struct CACorrParamsEdited {
    bool red;
    bool blue;
};

struct ResizeParamsEdited {
    bool scale;
    bool appliesTo;
    bool method;
    bool dataspec;
    bool width;
    bool height;
    bool longedge;
    bool shortedge;
    bool enabled;
    bool allowUpscaling;
};

struct ColorManagementParamsEdited {
    bool inputProfile;
    bool toneCurve;
    bool applyLookTable;
    bool applyBaselineExposureOffset;
    bool applyHueSatMap;
    bool dcpIlluminant;

    bool workingProfile;
    bool workingTRC;
    bool workingTRCGamma;
    bool workingTRCSlope;

    bool outputProfile;
    bool outputIntent;
    bool outputBPC;
};

struct WaveletParamsEdited {
    bool enabled;
    bool strength;
    bool balance;
    bool iter;
    bool median;
    bool medianlev;
    bool linkedg;
    bool cbenab;
    bool lipst;
    bool Medgreinf;
    bool avoid;
    bool tmr;
    bool c[9];
    bool ch[9];
    bool Lmethod;
    bool CHmethod;
    bool CHSLmethod;
    bool EDmethod;
    bool BAmethod;
    bool NPmethod;
    bool TMmethod;
    bool HSmethod;
    bool CLmethod;
    bool Backmethod;
    bool Tilesmethod;
    bool daubcoeffmethod;
    bool Dirmethod;
    bool rescon;
    bool resconH;
    bool reschro;
    bool tmrs;
    bool gamma;
    bool sup;
    bool sky;
    bool thres;
    bool threshold;
    bool threshold2;
    bool edgedetect;
    bool edgedetectthr;
    bool edgedetectthr2;
    bool edgesensi;
    bool edgeampli;
    bool chro;
    bool chroma;
    bool contrast;
    bool edgrad;
    bool edgval;
    bool edgthresh;
    bool thr;
    bool thrH;
    bool skinprotect;
    bool hueskin;
    bool hueskin2;
    bool hllev;
    bool bllev;
    bool edgcont;
    bool level0noise;
    bool level1noise;
    bool level2noise;
    bool level3noise;
    bool ccwcurve;
    bool opacityCurveBY;
    bool opacityCurveRG;
    bool opacityCurveW;
    bool opacityCurveWL;
    bool hhcurve;
    bool Chcurve;
    bool pastlev;
    bool satlev;
    bool wavclCurve;
    bool greenlow;
    bool bluelow;
    bool greenmed;
    bool bluemed;
    bool greenhigh;
    bool bluehigh;
    bool expcontrast;
    bool expchroma;
    bool expedge;
    bool expresid;
    bool expfinal;
    bool exptoning;
    bool expnoise;
};

struct DirPyrEqualizerParamsEdited {
    bool enabled;
    bool gamutlab;
    bool mult[6];
    bool cbdlMethod;
    bool threshold;
    bool skinprotect;
    bool hueskin;
};

struct HSVEqualizerParamsEdited {
    bool enabled;
    bool hcurve;
    bool scurve;
    bool vcurve;
};

struct FilmSimulationParamsEdited {
    bool enabled;
    bool clutFilename;
    bool strength;
};

struct SoftLightParamsEdited {
    bool enabled;
    bool strength;
};

struct DehazeParamsEdited {
    bool enabled;
    bool strength;
    bool showDepthMap;
    bool depth;
    bool luminance;
};

struct RAWParamsEdited {
    struct BayerSensor {
        bool method;
        bool border;
        bool imageNum;
        bool ccSteps;
        bool exBlack0;
        bool exBlack1;
        bool exBlack2;
        bool exBlack3;
        bool exTwoGreen;
        bool dcbIterations;
        bool dcbEnhance;
        bool lmmseIterations;
        bool dualDemosaicAutoContrast;
        bool dualDemosaicContrast;
        bool pixelShiftMotionCorrectionMethod;
        bool pixelShiftEperIso;
        bool pixelShiftSigma;
        bool pixelShiftShowMotion;
        bool pixelShiftShowMotionMaskOnly;
        bool pixelShiftHoleFill;
        bool pixelShiftMedian;
        bool pixelShiftGreen;
        bool pixelShiftBlur;
        bool pixelShiftSmooth;
        bool pixelShiftEqualBright;
        bool pixelShiftEqualBrightChannel;
        bool pixelShiftNonGreenCross;
        bool pixelShiftDemosaicMethod;

        bool greenEq;
        bool linenoise;
        bool linenoiseDirection;
        bool pdafLinesFilter;

        bool isUnchanged() const;
    };

    struct XTransSensor {
        bool method;
        bool dualDemosaicAutoContrast;
        bool dualDemosaicContrast;
        bool border;
        bool ccSteps;
        bool exBlackRed;
        bool exBlackGreen;
        bool exBlackBlue;

        bool isUnchanged() const;
    };

    BayerSensor bayersensor;
    XTransSensor xtranssensor;

    bool ca_autocorrect;
    bool ca_avoidcolourshift;
    bool caautoiterations;
    bool cared;
    bool cablue;
    bool hotPixelFilter;
    bool deadPixelFilter;
    bool hotdeadpix_thresh;
    bool darkFrame;
    bool df_autoselect;
    bool ff_file;
    bool ff_AutoSelect;
    bool ff_BlurRadius;
    bool ff_BlurType;
    bool ff_AutoClipControl;
    bool ff_clipControl;
    bool exPos;

    bool isUnchanged() const;
};


struct MetaDataParamsEdited {
    bool mode;
};

struct FilmNegativeParamsEdited {
    bool enabled;
    bool redRatio;
    bool greenExp;
    bool blueRatio;

    bool isUnchanged() const;
};

struct ParamsEdited {
    GeneralParamsEdited general;
    ToneCurveParamsEdited toneCurve;
    LCurveParamsEdited labCurve;
    LocalContrastParamsEdited localContrast;
    RGBCurvesParamsEdited rgbCurves;
    ColorToningEdited colorToning;
    RetinexParamsEdited retinex;
    SharpeningParamsEdited sharpening;
    CaptureSharpeningParamsEdited pdsharpening;
    SharpeningParamsEdited prsharpening;
    SharpenEdgeParamsEdited sharpenEdge;
    SharpenMicroParamsEdited sharpenMicro;
    VibranceParamsEdited vibrance;
    ColorAppearanceParamsEdited colorappearance;
    WBParamsEdited wb;
    DefringeParamsEdited defringe;
    DirPyrDenoiseParamsEdited dirpyrDenoise;
    EPDParamsEdited epd;
    FattalToneMappingParamsEdited fattal;
    ImpulseDenoiseParamsEdited impulseDenoise;
    SHParamsEdited sh;
    CropParamsEdited crop;
    CoarseTransformParamsEdited coarse;
    CommonTransformParamsEdited commonTrans;
    RotateParamsEdited rotate;
    DistortionParamsEdited distortion;
    LensProfParamsEdited lensProf;
    PerspectiveParamsEdited perspective;
    GradientParamsEdited gradient;
    PCVignetteParamsEdited pcvignette;
    CACorrParamsEdited cacorrection;
    VignettingParamsEdited vignetting;
    ChannelMixerParamsEdited chmixer;
    BlackWhiteParamsEdited blackwhite;
    ResizeParamsEdited resize;
    ColorManagementParamsEdited icm;
    RAWParamsEdited raw;
    DirPyrEqualizerParamsEdited dirpyrequalizer;
    WaveletParamsEdited wavelet;
    HSVEqualizerParamsEdited hsvequalizer;
    FilmSimulationParamsEdited filmSimulation;
    SoftLightParamsEdited softlight;
    DehazeParamsEdited dehaze;
    MetaDataParamsEdited metadata;
    FilmNegativeParamsEdited filmNegative;
    bool exif;
    bool iptc;

    explicit ParamsEdited(bool value = false);

    void set(bool v);
    void initFrom(const std::vector<rtengine::procparams::ProcParams>& src);
    void combine(rtengine::procparams::ProcParams& toEdit, const rtengine::procparams::ProcParams& mods, bool forceSet);
};
