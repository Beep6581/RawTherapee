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
#ifndef _PARAMEDITED_H_
#define _PARAMEDITED_H_

#include <glibmm.h>
#include <vector>
#include "../rtengine/procparams.h"
#include "../rtengine/rtengine.h"

class GeneralParamsEdited {

    public:
        bool rank;
        bool colorlabel;
        bool intrash;
};

class ToneCurveParamsEdited {

    public:
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
};

class LCurveParamsEdited {
    public:
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

        bool enabled;
        bool method;
};

class RGBCurvesParamsEdited {

    public:
        bool lumamode;
        bool rcurve;
        bool gcurve;
        bool bcurve;
};

class ColorToningEdited {

    public:
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
};

class SharpenEdgeParamsEdited {

    public :
        bool enabled;
        bool passes;
        bool amount;
        bool threechannels;
};

class SharpenMicroParamsEdited {
    public :
        bool enabled;
        bool matrix;
        bool amount;
        bool uniformity;
		
};

class SharpeningParamsEdited {

    public:
        bool enabled;
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

class VibranceParamsEdited {

    public:
        bool enabled;
        bool pastels;
        bool saturated;
        bool psthreshold;
        bool protectskins;
        bool avoidcolorshift;
        bool pastsattog;
        bool skintonescurve;
};

/*class ColorBoostParamsEdited {

    public: 
        bool amount;
        bool avoidclip;
        bool enable_saturationlimiter;
        bool rstprotection;
};*/

class WBParamsEdited {

    public:
        bool method;
        bool temperature;
        bool green;
        bool equal;
};

/*class ColorShiftParamsEdited {

    public:
        bool a;
        bool b;
};*/

/*class LumaDenoiseParamsEdited {

    public:
        bool enabled;
        bool radius;
        bool edgetolerance;
};*/

/*class ColorDenoiseParamsEdited {

    public:
        bool enabled;
        bool amount;
};*/

class DefringeParamsEdited {

public:
    bool enabled;
    bool radius;
    bool threshold;
    bool huecurve;
};

class ImpulseDenoiseParamsEdited {

public:
    bool enabled;
    bool thresh;
};

class ColorAppearanceParamsEdited {

public:
    bool curve;
    bool curve2;
    bool curve3;
    bool curveMode;
    bool curveMode2;
    bool curveMode3;
    bool enabled;
    bool degree;
    bool autodegree;
    bool autoadapscen;
    bool surround;
    bool adapscen;
    bool adaplum;
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
//  bool badpix;
    bool datacie;
    bool tonecie;
//  bool sharpcie;
};

class DirPyrDenoiseParamsEdited {

public:
    bool enabled;
    bool enhance;
	bool median;
	bool autochroma;
    bool Ldetail;
    bool luma;
    bool chroma;
    bool redchro;
    bool bluechro;
    bool gamma;
    bool lcurve;
    bool cccurve;
	
//    bool perform;
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

class EPDParamsEdited{
public:
    bool enabled;
    bool strength;
    bool edgeStopping;
    bool scale;
    bool reweightingIterates;
};


class SHParamsEdited {

    public:
        bool enabled;
        bool hq;
        bool highlights;
        bool htonalwidth;
        bool shadows;
        bool stonalwidth;
        bool localcontrast;
        bool radius;
};

class CropParamsEdited {

    public:
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

class CoarseTransformParamsEdited {

    public:
        bool rotate;
        bool hflip;
        bool vflip;
};

class CommonTransformParamsEdited {

   public:
        bool autofill;
};

class RotateParamsEdited {
    
    public:
        bool degree;
};

class DistortionParamsEdited {

    public:
        bool amount;
};

class LensProfParamsEdited {
    public:
        bool lcpFile,useDist,useVign,useCA; 

        bool isUnchanged() const;
};

class PerspectiveParamsEdited {

    public:
        bool horizontal;
        bool vertical;
};

class GradientParamsEdited {

    public:
        bool enabled;
        bool degree;
        bool feather;
        bool strength;
        bool centerX;
        bool centerY;
};

class PCVignetteParamsEdited {

    public:
        bool enabled;
        bool strength;
        bool feather;
        bool roundness;
};

class VignettingParamsEdited {

    public:
        bool amount;
        bool radius;
        bool strength;
        bool centerX;
        bool centerY;
};

class ChannelMixerParamsEdited {

    public:
        bool red[3];
        bool green[3];
        bool blue[3];

};
class BlackWhiteParamsEdited {

    public:
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

class CACorrParamsEdited {

    public:
        bool red;
        bool blue;
};
/*
class HRecParamsEdited {

    public:
        bool enabled;
        bool method;
};
*/
class ResizeParamsEdited {

    public:
        bool scale;
        bool appliesTo;
        bool method;
        bool dataspec;
        bool width;
        bool height;
        bool enabled;
};

class ColorManagementParamsEdited {

    public:
        bool input;
        bool toneCurve;
        bool blendCMSMatrix;
        bool dcpIlluminant;
        bool working;
        bool output;
        bool gamma;
        bool gampos;
        bool slpos;
        bool gamfree;
        bool freegamma;
};
class WaveletParamsEdited {

    public:
        bool enabled;
        bool strength;
        bool median;
        bool avoid;
		bool c[9];
		bool Lmethod;
		bool CHmethod;
		bool HSmethod;
		bool CLmethod;
		bool Tilesmethod;
		bool Dirmethod;
		bool rescon;
		bool resconH;
		bool reschro;
		bool sup;
		bool sky;
		bool thres;
		bool threshold;
		bool threshold2;
		bool chroma;
		bool chro;
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
		bool clvcurve;
        bool opacityCurveBY;
        bool opacityCurveRG;		
        bool pastlev;
        bool satlev;	
};

class DirPyrEqualizerParamsEdited {

    public:
        bool enabled;
        bool gamutlab;
        bool mult[6];
		
        bool threshold;
        bool skinprotect;
        bool hueskin;
   //     bool algo;
};

class HSVEqualizerParamsEdited {

    public:
        bool hcurve;
        bool scurve;
        bool vcurve;
};

class FilmSimulationParamsEdited {
public:
    bool enabled;
    bool clutFilename;
    bool strength;
};

class RAWParamsEdited {

    public:
        class BayerSensor {

            public:
                bool method;
                bool ccSteps;
                bool exBlack0;
                bool exBlack1;
                bool exBlack2;
                bool exBlack3;
                bool exTwoGreen;
                bool dcbIterations;
                bool dcbEnhance;
                bool lmmseIterations;
                //bool allEnhance;
                bool greenEq;
                bool linenoise;

                bool isUnchanged() const;
        };

        class XTransSensor {

            public:
                bool method;
                bool ccSteps;
                bool exBlackRed;
                bool exBlackGreen;
                bool exBlackBlue;

                bool isUnchanged() const;
        };

        BayerSensor bayersensor;
        XTransSensor xtranssensor;

        bool caCorrection;
        bool caRed;
        bool caBlue;
        bool hotPixelFilter;
        bool deadPixelFilter;
        bool hotDeadPixelThresh;
        bool darkFrame;
        bool dfAuto;
        bool ff_file;
        bool ff_AutoSelect;
        bool ff_BlurRadius;
        bool ff_BlurType;
        bool ff_AutoClipControl;
        bool ff_clipControl;
        bool exPos;
        bool exPreser;

        bool isUnchanged() const;
};

class ParamsEdited {

    public:
        GeneralParamsEdited           general;
        ToneCurveParamsEdited         toneCurve;
        LCurveParamsEdited            labCurve;
        RGBCurvesParamsEdited         rgbCurves;
        ColorToningEdited             colorToning;
        SharpeningParamsEdited        sharpening;
        SharpenEdgeParamsEdited       sharpenEdge;
        SharpenMicroParamsEdited      sharpenMicro;
        VibranceParamsEdited          vibrance;
        ColorAppearanceParamsEdited   colorappearance;
        //ColorBoostParamsEdited        colorBoost;
        WBParamsEdited                wb;
        //ColorShiftParamsEdited        colorShift;
        //LumaDenoiseParamsEdited       lumaDenoise;
        //ColorDenoiseParamsEdited      colorDenoise;
        DefringeParamsEdited          defringe;
        DirPyrDenoiseParamsEdited     dirpyrDenoise;
        EPDParamsEdited               epd;
        ImpulseDenoiseParamsEdited    impulseDenoise;
        SHParamsEdited                sh;
        CropParamsEdited              crop;
        CoarseTransformParamsEdited   coarse;
        CommonTransformParamsEdited   commonTrans;
        RotateParamsEdited            rotate;
        DistortionParamsEdited        distortion;
        LensProfParamsEdited          lensProf;
        PerspectiveParamsEdited       perspective;
        GradientParamsEdited          gradient;
        PCVignetteParamsEdited        pcvignette;
        CACorrParamsEdited            cacorrection;
        VignettingParamsEdited        vignetting;
        ChannelMixerParamsEdited      chmixer;
        BlackWhiteParamsEdited        blackwhite;
        ResizeParamsEdited            resize;
        ColorManagementParamsEdited   icm;
        RAWParamsEdited               raw;
        DirPyrEqualizerParamsEdited   dirpyrequalizer;
        WaveletParamsEdited     	    wavelet;
        HSVEqualizerParamsEdited      hsvequalizer;
        FilmSimulationParamsEdited    filmSimulation;
        bool                          exif;
        bool                          iptc;

        ParamsEdited ();

        void set   (bool v);
        void initFrom (const std::vector<rtengine::procparams::ProcParams>& src);
        void combine (rtengine::procparams::ProcParams& toEdit, const rtengine::procparams::ProcParams& mods, bool forceSet);

        bool operator== (const ParamsEdited& other);
        bool operator!= (const ParamsEdited& other);
};
#endif
