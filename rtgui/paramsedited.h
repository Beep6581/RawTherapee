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

#include <array>
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
    bool hlbl;
    bool hlth;
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
    bool complexmethod;
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
    bool gamutmunselmethod;
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
    bool labgridALowmerg;
    bool labgridBLowmerg;
    bool labgridAHighmerg;
    bool labgridBHighmerg;
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
    bool observer;
    bool tempBias;
    bool itcwb_rgreen;
    bool itcwb_nopurple;
    bool itcwb_alg;
    bool itcwb_prim;
    bool itcwb_sampling;
    bool itcwb_green;

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
    bool complexmethod;
    bool modelmethod;
    bool catmethod;
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
    bool illum;
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
    bool autotempout;
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

struct ToneEqualizerParamsEdited {
    bool enabled;
    std::array<bool, 5> bands;
    bool regularization;
    bool show_colormap;
    bool pivot;
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
class LocallabParamsEdited
{
public:
    struct LocallabSpotEdited {
        // Control spot settings
        bool name;
        bool isvisible;
        bool prevMethod;
        bool shape;
        bool spotMethod;
        bool wavMethod;
        bool sensiexclu;
        bool structexclu;
        bool struc;
        bool shapeMethod;
        bool avoidgamutMethod;
        bool loc;
        bool centerX;
        bool centerY;
        bool circrad;
        bool qualityMethod;
        bool complexMethod;
        bool transit;
        bool feather;
        bool thresh;
        bool iter;
        bool balan;
        bool balanh;
        bool colorde;
        bool colorscope;
        bool avoidrad;
        bool transitweak;
        bool transitgrad;
        bool hishow;
        bool activ;
        bool blwh;
        bool recurs;
        bool laplac;
        bool deltae;
        bool shortc;
        bool savrest;
        bool scopemask;
        bool denoichmask;
        bool lumask;
        // Color & Light
        bool visicolor;
        bool expcolor;
        bool complexcolor;
        bool curvactiv;
        bool lightness;
        bool reparcol;
        bool gamc;
        bool contrast;
        bool chroma;
        bool labgridALow;
        bool labgridBLow;
        bool labgridAHigh;
        bool labgridBHigh;
        bool labgridALowmerg;
        bool labgridBLowmerg;
        bool labgridAHighmerg;
        bool labgridBHighmerg;
        bool strengthgrid;
        bool sensi;
        bool structcol;
        bool strcol;
        bool strcolab;
        bool strcolh;
        bool angcol;
        bool blurcolde;
        bool blurcol;
        bool contcol;
        bool blendmaskcol;
        bool radmaskcol;
        bool chromaskcol;
        bool gammaskcol;
        bool slomaskcol;
        bool shadmaskcol;
        bool strumaskcol;
        bool lapmaskcol;
        bool qualitycurveMethod;
        bool gridMethod;
        bool merMethod;
        bool toneMethod;
        bool mergecolMethod;
        bool llcurve;
        bool lccurve;
        bool cccurve;
        bool clcurve;
        bool rgbcurve;
        bool LHcurve;
        bool HHcurve;
        bool CHcurve;
        bool invers;
        bool special;
        bool toolcol;
        bool enaColorMask;
        bool fftColorMask;
        bool CCmaskcurve;
        bool LLmaskcurve;
        bool HHmaskcurve;
        bool HHhmaskcurve;
        bool softradiuscol;
        bool opacol;
        bool mercol;
        bool merlucol;
        bool conthrcol;
        bool Lmaskcurve;
        bool LLmaskcolcurvewav;
        bool csthresholdcol;
        bool recothresc;
        bool lowthresc;
        bool higthresc;
        bool decayc;
        // Exposure
        bool visiexpose;
        bool expexpose;
        bool complexexpose;
        bool expcomp;
        bool hlcompr;
        bool hlcomprthresh;
        bool black;
        bool shadex;
        bool shcompr;
        bool expchroma;
        bool sensiex;
        bool structexp;
        bool gamex;
        bool blurexpde;
        bool strexp;
        bool angexp;
        bool excurve;
        bool norm;
        bool inversex;
        bool enaExpMask;
        bool enaExpMaskaft;
        bool CCmaskexpcurve;
        bool LLmaskexpcurve;
        bool HHmaskexpcurve;
        bool blendmaskexp;
        bool radmaskexp;
        bool chromaskexp;
        bool gammaskexp;
        bool slomaskexp;
        bool lapmaskexp;
        bool strmaskexp;
        bool angmaskexp;
        bool softradiusexp;
        bool Lmaskexpcurve;
        bool expMethod;
        bool exnoiseMethod;
        bool laplacexp;
        bool reparexp;
        bool balanexp;
        bool linear;
        bool gamm;
        bool fatamount;
        bool fatdetail;
        bool fatanchor;
        bool fatlevel;
        bool recothrese;
        bool lowthrese;
        bool higthrese;
        bool decaye;
        // Shadow highlight
        bool visishadhigh;
        bool expshadhigh;
        bool complexshadhigh;
        bool shMethod;
        bool multsh[6];
        bool highlights;
        bool h_tonalwidth;
        bool shadows;
        bool s_tonalwidth;
        bool sh_radius;
        bool sensihs;
        bool enaSHMask;
        bool CCmaskSHcurve;
        bool LLmaskSHcurve;
        bool HHmaskSHcurve;
        bool blendmaskSH;
        bool radmaskSH;
        bool blurSHde;
        bool strSH;
        bool angSH;
        bool inverssh;
        bool chromaskSH;
        bool gammaskSH;
        bool slomaskSH;
        bool lapmaskSH;
        bool detailSH;
        bool tePivot;
        bool reparsh;
        bool LmaskSHcurve;
        bool fatamountSH;
        bool fatanchorSH;
        bool gamSH;
        bool sloSH;
        bool recothress;
        bool lowthress;
        bool higthress;
        bool decays;
        // Vibrance
        bool visivibrance;
        bool expvibrance;
        bool complexvibrance;
        bool saturated;
        bool pastels;
        bool vibgam;
        bool warm;
        bool psthreshold;
        bool protectskins;
        bool avoidcolorshift;
        bool pastsattog;
        bool sensiv;
        bool skintonescurve;
        bool CCmaskvibcurve;
        bool LLmaskvibcurve;
        bool HHmaskvibcurve;
        bool enavibMask;
        bool blendmaskvib;
        bool radmaskvib;
        bool chromaskvib;
        bool gammaskvib;
        bool slomaskvib;
        bool lapmaskvib;
        bool strvib;
        bool strvibab;
        bool strvibh;
        bool angvib;
        bool Lmaskvibcurve;
        bool recothresv;
        bool lowthresv;
        bool higthresv;
        bool decayv;
        // Soft Light
        bool visisoft;
        bool expsoft;
        bool complexsoft;
        bool streng;
        bool sensisf;
        bool laplace;
        bool softMethod;
        // Blur & Noise
        bool visiblur;
        bool expblur;
        bool complexblur;
        bool radius;
        bool strength;
        bool sensibn;
        bool itera;
        bool guidbl;
        bool strbl;
        bool recothres;
        bool lowthres;
        bool higthres;
        bool recothresd;
        bool lowthresd;
        bool midthresd;
        bool midthresdch;
        bool higthresd;
        bool decayd;
        bool isogr;
        bool strengr;
        bool scalegr;
        bool divgr;
        bool epsbl;
        bool blMethod;
        bool chroMethod;
        bool quamethod;
        bool usemask;
        bool invmaskd;
        bool invmask;
        bool levelthr;
        bool lnoiselow;
        bool levelthrlow;
        bool blurMethod;
        bool medMethod;
        bool activlum;
        bool noiselumf;
        bool noiselumf0;
        bool noiselumf2;
        bool noiselumc;
        bool noiselumdetail;
        bool noiselequal;
        bool noisegam;
        bool noisechrof;
        bool noisechroc;
        bool noisechrodetail;
        bool adjblur;
        bool bilateral;
        bool nlstr;
        bool nldet;
        bool nlpat;
        bool nlrad;
        bool nlgam;
        bool sensiden;
        bool reparden;
        bool detailthr;
        bool locwavcurveden;
        bool locwavcurvehue;
        bool showmaskblMethodtyp;
        bool CCmaskblcurve;
        bool LLmaskblcurve;
        bool HHmaskblcurve;
        bool enablMask;
        bool fftwbl;
        bool invbl;
        bool toolbl;
        bool blendmaskbl;
        bool radmaskbl;
        bool chromaskbl;
        bool gammaskbl;
        bool slomaskbl;
        bool lapmaskbl;
        bool shadmaskbl;
        bool shadmaskblsha;
        bool strumaskbl;
        bool Lmaskblcurve;
        bool LLmaskblcurvewav;
        bool csthresholdblur;
        // Tone Mapping
        bool visitonemap;
        bool exptonemap;
        bool complextonemap;
        bool stren;
        bool gamma;
        bool estop;
        bool scaltm;
        bool repartm;
        bool rewei;
        bool satur;
        bool sensitm;
        bool softradiustm;
        bool amount;
        bool equiltm;
        bool CCmasktmcurve;
        bool LLmasktmcurve;
        bool HHmasktmcurve;
        bool enatmMask;
        bool enatmMaskaft;
        bool blendmasktm;
        bool radmasktm;
        bool chromasktm;
        bool gammasktm;
        bool slomasktm;
        bool lapmasktm;
        bool Lmasktmcurve;
        bool recothrest;
        bool lowthrest;
        bool higthrest;
        bool decayt;
        // Retinex
        bool visireti;
        bool expreti;
        bool complexreti;
        bool retinexMethod;
        bool str;
        bool chrrt;
        bool neigh;
        bool vart;
        bool offs;
        bool dehaz;
        bool depth;
        bool sensih;
        bool localTgaincurve;
        bool localTtranscurve;
        bool inversret;
        bool equilret;
        bool loglin;
        bool dehazeSaturation;
        bool softradiusret;
        bool CCmaskreticurve;
        bool LLmaskreticurve;
        bool HHmaskreticurve;
        bool enaretiMask;
        bool enaretiMasktmap;
        bool blendmaskreti;
        bool radmaskreti;
        bool chromaskreti;
        bool gammaskreti;
        bool slomaskreti;
        bool lapmaskreti;
        bool scalereti;
        bool darkness;
        bool lightnessreti;
        bool limd;
        bool cliptm;
        bool fftwreti;
        bool Lmaskreticurve;
        bool recothresr;
        bool lowthresr;
        bool higthresr;
        bool decayr;
        // Sharpening
        bool visisharp;
        bool expsharp;
        bool complexsharp;
        bool sharcontrast;
        bool sharradius;
        bool sharamount;
        bool shardamping;
        bool shariter;
        bool sharblur;
        bool shargam;
        bool sensisha;
        bool inverssha;
        // Local Contrast
        bool visicontrast;
        bool expcontrast;
        bool complexcontrast;
        bool lcradius;
        bool lcamount;
        bool lcdarkness;
        bool lclightness;
        bool sigmalc;
        bool levelwav;
        bool residcont;
        bool residsha;
        bool residshathr;
        bool residhi;
        bool residhithr;
        bool gamlc;
        bool residgam;
        bool residslop;
        bool residblur;
        bool levelblur;
        bool sigmabl;
        bool residchro;
        bool residcomp;
        bool sigma;
        bool offset;
        bool sigmadr;
        bool threswav;
        bool chromalev;
        bool chromablu;
        bool sigmadc;
        bool deltad;
        bool fatres;
        bool clarilres;
        bool claricres;
        bool clarisoft;
        bool sigmalc2;
        bool strwav;
        bool angwav;
        bool strengthw;
        bool sigmaed;
        bool radiusw;
        bool detailw;
        bool gradw;
        bool tloww;
        bool thigw;
        bool edgw;
        bool basew;
        bool sensilc;
        bool reparw;
        bool fftwlc;
        bool blurlc;
        bool wavblur;
        bool wavedg;
        bool waveshow;
        bool wavcont;
        bool wavcomp;
        bool wavgradl;
        bool wavcompre;
        bool origlc;
        bool localcontMethod;
        bool localedgMethod;
        bool localneiMethod;
        bool locwavcurve;
        bool csthreshold;
        bool loclevwavcurve;
        bool locconwavcurve;
        bool loccompwavcurve;
        bool loccomprewavcurve;
        bool locedgwavcurve;
        bool CCmasklccurve;
        bool LLmasklccurve;
        bool HHmasklccurve;
        bool enalcMask;
        bool blendmasklc;
        bool radmasklc;
        bool chromasklc;
        bool Lmasklccurve;
        bool recothresw;
        bool lowthresw;
        bool higthresw;
        bool decayw;
        // Contrast by detail levels
        bool visicbdl;
        bool expcbdl;
        bool complexcbdl;
        bool mult[6];
        bool chromacbdl;
        bool threshold;
        bool sensicb;
        bool clarityml;
        bool contresid;
        bool softradiuscb;
        bool enacbMask;
        bool CCmaskcbcurve;
        bool LLmaskcbcurve;
        bool HHmaskcbcurve;
        bool blendmaskcb;
        bool radmaskcb;
        bool chromaskcb;
        bool gammaskcb;
        bool slomaskcb;
        bool lapmaskcb;
        bool Lmaskcbcurve;
        bool recothrescb;
        bool lowthrescb;
        bool higthrescb;
        bool decaycb;
        // Log encoding
        bool visilog;
        bool explog;
        bool complexlog;
        bool autocompute;
        bool sourceGray;
        bool sourceabs;
        bool targabs;
        bool targetGray;
        bool catad;
        bool saturl;
        bool chroml;
        bool lightl;
        bool lightq;
        bool contl;
        bool contthres;
        bool contq;
        bool colorfl;
        bool LcurveL;
        bool Autogray;
        bool fullimage;
        bool repar;
        bool ciecam;
        bool blackEv;
        bool whiteEv;
        bool detail;
        bool sursour;
        bool surround;
        bool sensilog;
        bool baselog;
        bool strlog;
        bool anglog;
        bool CCmaskcurveL;
        bool LLmaskcurveL;
        bool HHmaskcurveL;
        bool enaLMask;
        bool blendmaskL;
        bool radmaskL;
        bool chromaskL;
        bool LmaskcurveL;
        bool recothresl;
        bool lowthresl;
        bool higthresl;
        bool decayl;
        //mask
        bool visimask;
        bool complexmask;
        bool expmask;
        bool sensimask;
        bool blendmask;
        bool blendmaskab;
        bool softradiusmask;
        bool enamask;
        bool fftmask;
        bool blurmask;
        bool contmask;
        bool CCmask_curve;
        bool LLmask_curve;
        bool HHmask_curve;
        bool strumaskmask;
        bool toolmask;
        bool radmask;
        bool lapmask;
        bool chromask;
        bool gammask;
        bool slopmask;
        bool shadmask;
        bool str_mask;
        bool ang_mask;
        bool HHhmask_curve;
        bool Lmask_curve;
        bool LLmask_curvewav;
        bool csthresholdmask;
        //locallabcie
        bool visicie;
        bool complexcie;
        bool expcie;
        bool reparcie;
        bool sensicie;
        bool Autograycie;
        bool forcejz;
        bool forcebw;
        bool qtoj;
        bool jabcie;
        bool sigmoidqjcie;
        bool logcie;
        bool logjz;
        bool sigjz;
        bool sigq;
        bool chjzcie;
        bool sourceGraycie;
        bool sourceabscie;
        bool sursourcie;
        bool modecam;
        bool modecie;
        bool saturlcie;
        bool rstprotectcie;
        bool chromlcie;
        bool huecie;
        bool toneMethodcie;
        bool ciecurve;
        bool toneMethodcie2;
        bool ciecurve2;
        bool chromjzcie;
        bool saturjzcie;
        bool huejzcie;
        bool softjzcie;
        bool strsoftjzcie;
        bool thrhjzcie;
        bool jzcurve;
        bool czcurve;
        bool czjzcurve;
        bool HHcurvejz;
        bool CHcurvejz;
        bool LHcurvejz;
        bool lightlcie;
        bool lightjzcie;
        bool lightqcie;
        bool contlcie;
        bool contjzcie;
        bool adapjzcie;
        bool jz100;
        bool pqremap;
        bool pqremapcam16;
        bool hljzcie;
        bool hlthjzcie;
        bool shjzcie;
        bool shthjzcie;
        bool radjzcie;
        bool contthrescie;
        bool blackEvjz;
        bool whiteEvjz;
        bool targetjz;
        bool sigmoidldacie;
        bool sigmoidthcie;
        bool sigmoidblcie;
        bool sigmoidldajzcie;
        bool sigmoidthjzcie;
        bool sigmoidbljzcie;
        bool contqcie;
        bool colorflcie;
/*
        bool lightlzcam;
        bool lightqzcam;
        bool contlzcam;
        bool contqzcam; 
        bool contthreszcam;
        bool colorflzcam;
        bool saturzcam;
        bool chromzcam;
*/
        bool targabscie;
        bool targetGraycie;
        bool catadcie;
        bool detailcie;
        bool surroundcie;
        bool enacieMask;
        bool CCmaskciecurve;
        bool LLmaskciecurve;
        bool HHmaskciecurve;
        bool blendmaskcie;
        bool radmaskcie;
        bool sigmalcjz;
        bool clarilresjz;
        bool claricresjz;
        bool clarisoftjz;
        bool locwavcurvejz;
        bool csthresholdjz;
        bool chromaskcie;
        bool lapmaskcie;
        bool gammaskcie;
        bool slomaskcie;
        bool Lmaskciecurve;
        bool recothrescie;
        bool lowthrescie;
        bool higthrescie;
        bool decaycie;

        LocallabSpotEdited(bool v);

        void set(bool v);
    };

    bool enabled;
    bool selspot;
    std::vector<LocallabSpotEdited> spots;
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
    bool method;
    bool horizontal;
    bool vertical;
    bool camera_crop_factor;
    bool camera_focal_length;
    bool camera_pitch;
    bool camera_roll;
    bool camera_shift_horiz;
    bool camera_shift_vert;
    bool camera_yaw;
    bool projection_pitch;
    bool projection_rotate;
    bool projection_shift_horiz;
    bool projection_shift_vert;
    bool projection_yaw;
    bool control_lines;
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

class SpotParamsEdited
{
public:
    bool enabled;
    bool entries;
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
    bool will;
    bool wprim;
    bool redx;
    bool redy;
    bool grex;
    bool grey;
    bool blux;
    bool bluy;
    bool preser;
    bool fbw;
    bool gamut;
    bool labgridcieALow;
    bool labgridcieBLow;
    bool labgridcieAHigh;
    bool labgridcieBHigh;
    bool labgridcieGx;
    bool labgridcieGy;
    bool labgridcieWx;
    bool labgridcieWy;
    bool aRendIntent;
    bool outputProfile;
    bool outputIntent;
    bool outputBPC;
};

struct WaveletParamsEdited {
    bool enabled;
    bool strength;
    bool balance;
    bool iter;
    bool sigmafin;
    bool sigmaton;
    bool sigmacol;
    bool sigmadir;
    bool rangeab;
    bool protab;
    bool median;
    bool medianlev;
    bool linkedg;
    bool cbenab;
    bool lipst;
    bool Medgreinf;
    bool ushamethod;
    bool avoid;
    bool showmask;
    bool oldsh;
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
    bool complexmethod;
    //bool denmethod;
    bool mixmethod;
    bool slimethod;
    bool quamethod;
    bool daubcoeffmethod;
    bool Dirmethod;
    bool sigma;
    bool offset;
    bool lowthr;
    bool rescon;
    bool resconH;
    bool reschro;
    bool resblur;
    bool resblurc;
    bool tmrs;
    bool edgs;
    bool scale;
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
    bool edgeffect;
    bool edgval;
    bool edgthresh;
    bool thr;
    bool thrH;
    bool radius;
    bool skinprotect;
    bool hueskin;
    bool hueskin2;
    bool hllev;
    bool bllev;
    bool edgcont;
    bool chrwav;
    bool bluwav;
    bool level0noise;
    bool level1noise;
    bool level2noise;
    bool level3noise;
    bool leveldenoise;
    bool levelsigm;
    bool ccwcurve;
    bool blcurve;
    //bool opacityCurveSH;
    bool opacityCurveBY;
    bool wavdenoise;
    bool wavdenoiseh;
    bool opacityCurveRG;
    bool opacityCurveW;
    bool opacityCurveWL;
    bool hhcurve;
    bool wavguidcurve;
    bool wavhuecurve;
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
    bool ballum;
    bool sigm;
    bool levden;
    bool thrden;
    bool limden;
    bool balchrom;
    bool chromfi;
    bool chromco;
    bool mergeL;
    bool mergeC;
    bool softrad;
    bool softradend;
    bool strend;
    bool detend;
    bool thrend;
    bool expcontrast;
    bool expchroma;
    bool expedge;
    bool expbl;
    bool expresid;
    bool expfinal;
    bool exptoning;
    bool expnoise;
    bool expclari;
    bool labgridALow;
    bool labgridBLow;
    bool labgridAHigh;
    bool labgridBHigh;
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
    bool saturation;
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
        bool pixelShiftAverage;
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
    bool ff_FromMetaData;
    bool ff_BlurRadius;
    bool ff_BlurType;
    bool ff_AutoClipControl;
    bool ff_clipControl;
    bool exPos;

    struct PreprocessWBParamsEdited {
        bool mode;

        bool isUnchanged() const;
    };

    PreprocessWBParamsEdited preprocessWB;

    bool isUnchanged() const;
};


struct MetaDataParamsEdited {
    bool mode;
    bool exifKeys;
};

struct FilmNegativeParamsEdited {
    bool enabled;
    bool redRatio;
    bool greenExp;
    bool blueRatio;
    bool refInput;
    bool refOutput;
    bool colorSpace;

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
    ToneEqualizerParamsEdited toneEqualizer;
    CropParamsEdited crop;
    CoarseTransformParamsEdited coarse;
    CommonTransformParamsEdited commonTrans;
    RotateParamsEdited rotate;
    DistortionParamsEdited distortion;
    LensProfParamsEdited lensProf;
    PerspectiveParamsEdited perspective;
    GradientParamsEdited gradient;
    LocallabParamsEdited locallab;
    PCVignetteParamsEdited pcvignette;
    CACorrParamsEdited cacorrection;
    VignettingParamsEdited vignetting;
    ChannelMixerParamsEdited chmixer;
    BlackWhiteParamsEdited blackwhite;
    ResizeParamsEdited resize;
    SpotParamsEdited spot;
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
