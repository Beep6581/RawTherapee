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
        bool curveMode;
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
        bool bwtoning;
		bool lcredsk;
        bool cccurve;
        bool chcurve;
        bool lccurve;
};

class RGBCurvesParamsEdited {

    public:
        bool rcurve;
        bool gcurve;
        bool bcurve;
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


class WBParamsEdited {

    public:
        bool method;
        bool temperature;
        bool green;
};

class DefringeParamsEdited {
	
public:
	bool enabled;
	bool radius;
	bool threshold;
};

class ImpulseDenoiseParamsEdited {
	
public:
	bool enabled;
	bool thresh;

};

class DirPyrDenoiseParamsEdited {
	
public:
	bool enabled;
	bool Ldetail;
	bool luma;
	bool chroma;
	bool gamma;
};

class EPDParamsEdited{
public:
	bool enabled;
	bool Strength;
	bool EdgeStopping;
	bool Scale;
	bool ReweightingIterates;
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

class CACorrParamsEdited {

    public:
        bool red;
        bool blue;
};

class HRecParamsEdited {

    public:
        bool enabled;
		bool method;
};

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
        bool preferredProfile;
        bool working;
        bool output;
		bool gamma;
		bool freegamma;
		bool gampos;
		bool slpos;
};

class DirPyrEqualizerParamsEdited {
	
public:
	bool enabled;
	bool mult[8];
};

class HSVEqualizerParamsEdited {
	
public:
        bool hcurve;
        bool scurve;
        bool vcurve;
};

class RAWParamsEdited {

    public:
        bool ccSteps;
        bool dmethod;
        bool dcbIterations;
        bool dcbEnhance;
        //bool allEnhance;
        bool caCorrection;
        bool caRed;
        bool caBlue;
        bool greenEq;
        bool hotDeadPixelFilter;
        bool hotDeadPixelThresh;
        bool linenoise;
        bool darkFrame;
        bool dfAuto;
        bool ff_file;
        bool ff_AutoSelect;
        bool ff_BlurRadius;
        bool ff_BlurType;
        bool exPos;
        bool exPreser;
        bool exBlackzero;
        bool exBlackone;
        bool exBlacktwo;
        bool exBlackthree;
        bool exTwoGreen;

        bool isUnchanged() const;
};

class ParamsEdited {

    public:
        GeneralParamsEdited           general;
        ToneCurveParamsEdited         toneCurve;
        LCurveParamsEdited            labCurve;
        RGBCurvesParamsEdited         rgbCurves;
        SharpeningParamsEdited        sharpening;
        SharpenEdgeParamsEdited       sharpenEdge;
        SharpenMicroParamsEdited      sharpenMicro;
        VibranceParamsEdited          vibrance;
        WBParamsEdited                wb;
        DefringeParamsEdited          defringe;
        DirPyrDenoiseParamsEdited     dirpyrDenoise;
        EPDParamsEdited               edgePreservingDecompositionUI;
        ImpulseDenoiseParamsEdited    impulseDenoise;
        SHParamsEdited                sh;
        CropParamsEdited              crop;
        CoarseTransformParamsEdited   coarse;
        CommonTransformParamsEdited   commonTrans;
        RotateParamsEdited            rotate;
        DistortionParamsEdited        distortion;
        LensProfParamsEdited          lensProf;
        PerspectiveParamsEdited       perspective;
        CACorrParamsEdited            cacorrection;
        VignettingParamsEdited        vignetting;
        ChannelMixerParamsEdited      chmixer;
        HRecParamsEdited              hlrecovery;
        ResizeParamsEdited            resize;
        ColorManagementParamsEdited   icm;
        RAWParamsEdited               raw;
        DirPyrEqualizerParamsEdited   dirpyrequalizer;
        HSVEqualizerParamsEdited      hsvequalizer;
        bool                          exif;
        bool                          iptc;

        ParamsEdited () { set(false); }
        ParamsEdited (const bool setVal) { set(setVal); }

        void set   (bool v);
        void initFrom (const std::vector<rtengine::procparams::ProcParams>& src);
        void combine (rtengine::procparams::ProcParams& toEdit, const rtengine::procparams::ProcParams& mods, bool forceSet);

        ParamsEdited& operator&= (const ParamsEdited& rhs);
        ParamsEdited& operator|= (const ParamsEdited& rhs);
        bool operator== (const ParamsEdited& other);
        bool operator!= (const ParamsEdited& other);
};
#endif
