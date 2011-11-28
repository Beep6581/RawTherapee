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
#include <rtengine.h>
#include <procparams.h>

class ToneCurveParamsEdited {

    public:
        bool curve;
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
        bool saturation;
        bool avoidclip;
        bool enable_saturationlimiter;
        bool saturationlimit;
        bool lcurve;
        bool acurve;
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
};

class ColorBoostParamsEdited {

    public: 
        bool amount;
        bool avoidclip;
        bool enable_saturationlimiter;
        bool saturationlimit;
};

class WBParamsEdited {

    public:
        bool method;
        bool temperature;
        bool green;
};

class ColorShiftParamsEdited {

    public:
        bool a;
        bool b;
};

class LumaDenoiseParamsEdited {

    public:
        bool enabled;
        bool radius;
        bool edgetolerance;
};

class ColorDenoiseParamsEdited {

    public:
        bool enabled;
        bool amount;
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
		bool uselensfun;
        bool amount;
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
        bool blendCMSMatrix;
        bool working;
        bool output;
		bool gamma;
		bool gampos;
		bool slpos;
		bool gamfree;
		bool freegamma;		
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
        bool allEnhance;	
        bool caCorrection;
		bool caRed;
		bool caBlue;
        bool greenEq;
        bool hotDeadPixel;
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

class ExifPairEdited {

    public:
        Glib::ustring field;
        bool value;
};

class IPTCPairEdited {

    public:
        Glib::ustring field;
        bool values;
};

class ParamsEdited {

    public:
        ToneCurveParamsEdited         toneCurve;
        LCurveParamsEdited            labCurve;
        SharpeningParamsEdited        sharpening;
        SharpenEdgeParamsEdited       sharpenEdge;
        SharpenMicroParamsEdited      sharpenMicro;
        VibranceParamsEdited          vibrance;
        ColorBoostParamsEdited        colorBoost;
        WBParamsEdited                wb;
        ColorShiftParamsEdited        colorShift;
        LumaDenoiseParamsEdited       lumaDenoise;
        ColorDenoiseParamsEdited      colorDenoise;
        DefringeParamsEdited          defringe;
        DirPyrDenoiseParamsEdited     dirpyrDenoise;
        EPDParamsEdited					  edgePreservingDecompositionUI;
        ImpulseDenoiseParamsEdited    impulseDenoise;
        SHParamsEdited                sh;
        CropParamsEdited              crop;
        CoarseTransformParamsEdited   coarse;
        CommonTransformParamsEdited	  commonTrans;
        RotateParamsEdited            rotate;
        DistortionParamsEdited        distortion;
        PerspectiveParamsEdited		  perspective;
        CACorrParamsEdited            cacorrection;
        VignettingParamsEdited        vignetting;
        ChannelMixerParamsEdited      chmixer;
        HRecParamsEdited              hlrecovery;
        ResizeParamsEdited            resize;
        ColorManagementParamsEdited   icm;
        RAWParamsEdited               raw;
        DirPyrEqualizerParamsEdited   dirpyrequalizer;
        HSVEqualizerParamsEdited      hsvequalizer;
        std::vector<ExifPairEdited>   exif;
        std::vector<IPTCPairEdited>   iptc;

        ParamsEdited ();

        void set   (bool v);
        void initFrom (const std::vector<rtengine::procparams::ProcParams>& src);
        void combine (rtengine::procparams::ProcParams& toEdit, const rtengine::procparams::ProcParams& mods, bool forceSet);

        bool operator== (const ParamsEdited& other);
        bool operator!= (const ParamsEdited& other);
};
#endif
