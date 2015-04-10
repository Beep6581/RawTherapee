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
#ifndef _PROCPARAMS_H_
#define _PROCPARAMS_H_

#include <glibmm.h>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cmath>
#include "LUT.h"

class ParamsEdited;

namespace rtengine {

class ColorGradientCurve;
class OpacityCurve;
class NoiseCurve;
class WavCurve;
class WavOpacityCurveRG;
class WavOpacityCurveBY;

namespace procparams {

template <typename T>
class Threshold {
    public:
        T value[4];

    protected:
        bool initEq1;
        bool _isDouble;
#ifndef NDEBUG
        unsigned int part[5];
#endif
    public:
        Threshold (T bottom, T top, bool startAtOne) {
            initEq1 = startAtOne;
            value[0] = bottom;
            value[1] = top;
            value[2] = T(0);
            value[3] = T(0);
            _isDouble = false;
        }

        Threshold (T bottomLeft, T topLeft, T bottomRight, T topRight, bool startAtOne) {
            initEq1 = startAtOne;
            value[0] = bottomLeft;
            value[1] = topLeft;
            value[2] = bottomRight;
            value[3] = topRight;
            _isDouble = true;
        }

        // for convenience, since 'values' is public
        void setValues(T bottom, T top) {
            value[0] = bottom;
            value[1] = top;
        }

        // for convenience, since 'values' is public
        void setValues(T bottomLeft, T topLeft, T bottomRight, T topRight) {
            value[0] = bottomLeft;
            value[1] = topLeft;
            value[2] = bottomRight;
            value[3] = topRight;
        }

        bool isDouble() const { return _isDouble; }

        // RT: Type of the returned value
        // RV: Type of the value on the X axis
        // RV2: Type of the maximum value on the Y axis
        template <typename RT, typename RV, typename RV2>
        RT multiply(RV x, RV2 yMax) const {
        	double val = double(x);
            if (initEq1) {
                if (_isDouble) {
                	if (val == double(value[2]) && double(value[2]) == double(value[3]))
                		// this handle the special case where the 2 right values are the same, then bottom one is sent back,
                		// useful if one wants to keep the bottom value even beyond the x max bound
                        return RT(0.);
                    if (val >= double(value[3]))
                        return RT(yMax);
                    if (val > double(value[2]))
                        return RT(double(yMax)*(val-double(value[2]))/(double(value[3])-double(value[2])));
                }
                if (val >= double(value[0]))
                    return RT(0);
                if (val > double(value[1]))
                    return RT(double(yMax)*(1.-(val-double(value[0]))/(double(value[1])-double(value[0]))));
                return RT(yMax);
            }
            else {
                if (_isDouble) {
                	if (val == double(value[2]) && double(value[2]) == double(value[3]))
                		// this handle the special case where the 2 right values are the same, then top one is sent back,
                		// useful if one wants to keep the top value even beyond the x max bound
                        return RT(yMax);
                    if (val >= double(value[2]))
                        return RT(0);
                    if (val > double(value[3]))
                        return RT(double(yMax)*(1.-(val-double(value[3]))/(double(value[2])-double(value[3]))));
                }
                if (val >= double(value[1]))
                    return RT(yMax);
                if (val > double(value[0]))
                    return RT(double(yMax)*(val-double(value[0]))/(double(value[1])-double(value[0])));
                return RT(0);
            }
        }

        // RT: Type of the returned value
        // RV: Type of the value on the X axis
        /*template <typename RT, typename RV>
        RT getRatio(RV val) const {
        	double val = double(val);
            if (initEq1) {
                if (_isDouble) {  // assuming that simple thresholds will be more frequent
                    if (val >= double(value[3]))
                        return RT(1);
                    if (val > double(value[2]))
                        return (val-double(value[2]))/(double(value[3])-double(value[2]));
                }
                if (val >= double(value[1]))
                    return RT(0);
                if (val > double(value[0]))
                    return 1.-(val-double(value[0]))/(double(value[1])-double(value[0]));
                return RT(1);
            }
            else {
                if (_isDouble) {  // assuming that simple thresholds will be more frequent
                    if (val >= double(value[3]))
                        return RT(0);
                    if (val > double(value[2]))
                        return 1.-(val-double(value[2]))/(double(value[3])-double(value[2]));
                }
                if (val >= double(value[1]))
                    return RT(1);
                if (val > double(value[0]))
                    return (val-double(value[0]))/(double(value[1])-double(value[0]));
                return RT(0);
            }
        }*/

        Threshold<T> & operator= (const Threshold<T> &rhs) {
            value[0] = rhs.value[0];
            value[1] = rhs.value[1];
            value[2] = rhs.value[2];
            value[3] = rhs.value[3];
            initEq1 = rhs.initEq1;
            _isDouble = rhs._isDouble;
    	    return *this;
        }

        bool operator== (const Threshold<T> &rhs) const {
            if (_isDouble)
            	return fabs(value[0]-rhs.value[0])<1e-10
            	    && fabs(value[1]-rhs.value[1])<1e-10
            	    && fabs(value[2]-rhs.value[2])<1e-10
            	    && fabs(value[3]-rhs.value[3])<1e-10;
            else
            	return fabs(value[0]-rhs.value[0])<1e-10
            	    && fabs(value[1]-rhs.value[1])<1e-10;
        }
};

/**
  * Parameters of the tone curve 
  */
class ToneCurveParams {

    public:

        enum eTCModeId {
            TC_MODE_STD,               // Standard modes, the curve is applied on all component individually
            TC_MODE_WEIGHTEDSTD,       // Weighted standard mode
            TC_MODE_FILMLIKE,          // Film-like mode, as defined in Adobe's reference code
            TC_MODE_SATANDVALBLENDING  // Modify the Saturation and Value channel
        };

        bool        autoexp;
        double      clip;
        bool        hrenabled;  // Highlight Reconstruction enabled
        Glib::ustring method;   // Highlight Reconstruction's method
        double      expcomp;
        std::vector<double>   curve;
        std::vector<double>   curve2;
        eTCModeId   curveMode;
        eTCModeId   curveMode2;
        int         brightness;
        int         black;
        int         contrast;
        int         saturation;
        int         shcompr;
        int         hlcompr;        // Highlight Recovery's compression
        int         hlcomprthresh;  // Highlight Recovery's threshold

        ToneCurveParams () {
            setDefaults();
        }
        void setDefaults();
        static bool HLReconstructionNecessary(LUTu &histRedRaw, LUTu &histGreenRaw, LUTu &histBlueRaw);
};


/**
  * Parameters of the luminance curve
  */
class LCurveParams {

    public:
        std::vector<double>   lcurve;
        std::vector<double>   acurve;
        std::vector<double>   bcurve;
        std::vector<double>   cccurve;
        std::vector<double>   chcurve;
        std::vector<double>   lhcurve;
        std::vector<double>   hhcurve;
        std::vector<double>   lccurve;
        std::vector<double>   clcurve;
        int     brightness;
        int     contrast;
        int     chromaticity;
        bool    avoidcolorshift;
        double  rstprotection;
        bool    lcredsk;
};

/**
  * Parameters of the RGB curves
  */
class RGBCurvesParams {

    public:
        bool lumamode;
        std::vector<double>   rcurve;
        std::vector<double>   gcurve;
        std::vector<double>   bcurve;
};

/**
  * Parameters of the Color Toning
  */

class ColorToningParams {

    public:
        bool enabled;
        bool autosat;
        std::vector<double> opacityCurve;
        std::vector<double> colorCurve;
        int satProtectionThreshold;
        int saturatedOpacity;
        int strength;
        int balance;
        Threshold<int> hlColSat;
        Threshold<int> shadowsColSat;
        std::vector<double> clcurve;
        std::vector<double> cl2curve;

        /* Can be either:
         *  Splitlr    :
         *  Splitco    :
         *  Splitbal   :
         *  Lab        :
         *  Lch        :
         *  RGBSliders :
         *  RGBCurves  :
         */
        Glib::ustring method;

        /* Can be either:
         * Std   :
         * All   :
         * Separ :
         * Two   :
         */
        Glib::ustring twocolor;
        double redlow;
        double greenlow;
        double bluelow;
        double redmed;
        double greenmed;
        double bluemed;
        double redhigh;
        double greenhigh;
        double bluehigh;
        double satlow;
        double sathigh;
        bool lumamode;

        ColorToningParams ();
        void setDefaults();  // SHOULD BE GENERALIZED TO ALL CLASSES!
        /// @brief Transform the mixer values to their curve equivalences
        void mixerToCurve(std::vector<double> &colorCurve, std::vector<double> &opacityCurve) const;
        /// @brief Specifically transform the sliders values to their curve equivalences
        void slidersToCurve(std::vector<double> &colorCurve, std::vector<double> &opacityCurve) const;
        /// @brief Fill the ColorGradientCurve and OpacityCurve LUTf from the control points curve or sliders value
        void getCurves(ColorGradientCurve &colorCurveLUT, OpacityCurve &opacityCurveLUT, const double xyz_rgb[3][3], const double rgb_xyz[3][3], bool &opautili) const;

        static void getDefaultColorCurve(std::vector<double> &curve);
        static void getDefaultOpacityCurve(std::vector<double> &curve);
        static void getDefaultCLCurve(std::vector<double> &curve);
        static void getDefaultCL2Curve(std::vector<double> &curve);
};

/**
  * Parameters of the sharpening
  */
class SharpeningParams {

    public:
        bool           enabled;
        double         radius;
        int            amount;
        Threshold<int> threshold;
        bool           edgesonly;
        double         edges_radius;
        int            edges_tolerance;
        bool           halocontrol;
        int            halocontrol_amount;
        Glib::ustring  method;
        int            deconvamount;
        double         deconvradius;
        int            deconviter;
        int            deconvdamping;

        SharpeningParams() : threshold(20, 80, 2000, 1200, false) {};
};
class SharpenEdgeParams {
  public:
        bool    enabled;
        int     passes;
        double  amount;
        bool    threechannels;
};
class SharpenMicroParams {
  public:
        bool    enabled;
        bool    matrix;
        double  amount;
        double  uniformity;
};

/**
  * Parameters of the vibrance
  */
class VibranceParams {

    public:
        bool           enabled;
        int            pastels;
        int            saturated;
        Threshold<int> psthreshold;
        bool           protectskins;
        bool           avoidcolorshift;
        bool           pastsattog;
        std::vector<double> skintonescurve;

        VibranceParams() : psthreshold(0, 75,  false) {};
};

/**
  * Parameters of the color boost
  */
/*class ColorBoostParams {

    public: 
        int     amount;
        bool    avoidclip;
        bool    enable_saturationlimiter;
        double  saturationlimit;
};*/

/**
  * Parameters of the white balance adjustments
  */

enum WBTypes {
    WBT_CAMERA,
    WBT_AUTO,
    WBT_DAYLIGHT,
    WBT_CLOUDY,
    WBT_SHADE,
    WBT_WATER,
    WBT_TUNGSTEN,
    WBT_FLUORESCENT,
    WBT_LAMP,
    WBT_FLASH,
    WBT_LED,
    // WBT_CUSTOM one must remain the last one!
    WBT_CUSTOM
};

class WBEntry {
public:
    Glib::ustring ppLabel;
    enum WBTypes type;
    Glib::ustring GUILabel;
    int temperature;
    double green;
    double equal;

    WBEntry(Glib::ustring p, enum WBTypes t, Glib::ustring l, int temp, double green, double equal) : ppLabel(p), type(t), GUILabel(l), temperature(temp), green(green), equal(equal) {};
};

class WBParams {

    public:
        static std::vector<WBEntry*> wbEntries;
        Glib::ustring   method;
        int             temperature;
        double          green;
        double          equal;

        static void     init();
        static void     cleanup();
};

/**
 * Parameters of colorappearance
 */
class ColorAppearanceParams {

    public:
        enum eTCModeId {
            TC_MODE_LIGHT,    // Lightness mode
            TC_MODE_BRIGHT,   // Brightness mode
        };

        enum eCTCModeId {
            TC_MODE_CHROMA,   // chroma mode
            TC_MODE_SATUR,    // saturation mode
            TC_MODE_COLORF,   // colorfullness mode
       };

        bool          enabled;
        int           degree;
        bool          autodegree;
        std::vector<double> curve;
        std::vector<double> curve2;
        std::vector<double> curve3;
        eTCModeId     curveMode;
        eTCModeId     curveMode2;
        eCTCModeId    curveMode3;

        Glib::ustring surround;
        double        adapscen;
        bool          autoadapscen;

        double        adaplum;
        int           badpixsl;
        Glib::ustring wbmodel;
        Glib::ustring algo;
        double        contrast;
        double        qcontrast;
        double        jlight;
        double        qbright;
        double        chroma;
        double        schroma;
        double        mchroma;
        double        colorh;
        double        rstprotection;
        bool          surrsource;
        bool          gamut;
  //      bool          badpix;
        bool          datacie;
        bool          tonecie;
  //      bool          sharpcie;
    };

/**
  * Parameters of the color shift
  */
/*class ColorShiftParams {

    public:
        double  a;
        double  b;
};*/

/**
  * Parameters of the luminance denoising
  */
/*class LumaDenoiseParams {

    public:
        bool    enabled;
        double  radius;
        int     edgetolerance;
};*/

/**
  * Parameters of the color denoising
  */
/*class ColorDenoiseParams {

    public:
        bool    enabled;
        int     edgetolerance;
        bool    edgesensitive;
        int		amount;
};*/

/**
 * Parameters of defringing
 */
class DefringeParams {

    public:
        bool    enabled;
        double  radius;
        float   threshold;
        std::vector<double> huecurve;
};

/**
  * Parameters of impulse denoising
  */
class ImpulseDenoiseParams {

    public:
        bool    enabled;
        int     thresh;

};

/**
 * Parameters of the directional pyramid denoising
 */
class DirPyrDenoiseParams {

    public:
        std::vector<double>   lcurve;
        std::vector<double>   cccurve;
	
        bool    enabled;
        bool    enhance;
        bool    median;
        bool    autochroma;
		
        bool    perform;
        double  luma;
        double  Ldetail;
        double  chroma;
        double  redchro;
        double  bluechro;
        double  gamma;
        Glib::ustring dmethod;
        Glib::ustring Lmethod;
        Glib::ustring Cmethod;
        Glib::ustring C2method;
        Glib::ustring smethod;
        Glib::ustring medmethod;
        Glib::ustring methodmed;
        Glib::ustring rgbmethod;
        int  passes;
		
        DirPyrDenoiseParams ();
        void setDefaults();  // SHOULD BE GENERALIZED TO ALL CLASSES!
        void getCurves(NoiseCurve &lCurve, NoiseCurve &cCurve) const;

        static void getDefaultNoisCurve(std::vector<double> &curve);
        static void getDefaultCCCurve(std::vector<double> &curve);
		
};

//EPD related parameters.
class EPDParams{
    public:
        bool   enabled;
        double strength;
        double edgeStopping;
        double scale;
        int    reweightingIterates;
};

/**
  * Parameters of the shadow/highlight enhancement
  */
class SHParams {

    public:
        bool    enabled;
        bool    hq;
        int     highlights;
        int     htonalwidth;
        int     shadows;
        int     stonalwidth;
        int     localcontrast;
        int     radius;
};

/**
  * Parameters of the cropping
  */
class CropParams {

    public:
        bool    enabled;
        int     x;
        int     y;
        int     w;
        int     h;
        bool    fixratio;
        Glib::ustring   ratio;
        Glib::ustring   orientation;
        Glib::ustring   guide;

        CropParams() :enabled(false), x(0),y(0),w(0),h(0),fixratio(false) {};
        void mapToResized(int resizedWidth, int resizedHeight, int scale, int &x1, int &x2, int &y1, int &y2) const;
};

/**
  * Parameters of the coarse transformations like 90 deg rotations and h/v flipping
  */
class CoarseTransformParams {

    public:
        int     rotate;
        bool    hflip;
        bool    vflip;

        CoarseTransformParams() {
           setDefaults();
        }
        void setDefaults();
};

/**
  * Common transformation parameters
  */
class CommonTransformParams {

	public:
		bool autofill;
};

/**
  * Parameters of the rotation
  */
class RotateParams {
    
    public:
        double  degree;
};

/**
  * Parameters of the distortion correction
  */
class DistortionParams {

    public:
        double  amount;
};

// Lens profile correction parameters
class LensProfParams {
    
public:
    Glib::ustring lcpFile;
    bool useDist, useVign, useCA;

    LensProfParams() {
        setDefaults();
    }
    void setDefaults();
};

/**
  * Parameters of the perspective correction
  */
class PerspectiveParams {

    public:
        double  horizontal;
        double  vertical;
};

/**
  * Parameters of the gradient filter
  */
class GradientParams {

    public:
        bool   enabled;
        double degree;
        int    feather;
        double strength;
        int    centerX;
        int    centerY;
};

/**
  * Parameters of the post-crop vignette filter
  */
class PCVignetteParams {

    public:
        bool   enabled;
        double strength;
        int    feather;
        int    roundness;
};

/**
  * Parameters of the vignetting correction
  */
class VignettingParams {

    public:
        int  amount;
        int  radius;
        int  strength;
        int  centerX;
        int  centerY;
};

/**
  * Parameters of the color mixer
  */
class ChannelMixerParams {

    public:
        int red[3];
        int green[3];
        int blue[3];
};

class BlackWhiteParams {

    public:
        enum eTCModeId {
            TC_MODE_STD_BW,               // Standard modes, the curve is applied on all component individually
            TC_MODE_WEIGHTEDSTD_BW,       // Weighted standard mode
            TC_MODE_FILMLIKE_BW,          // Film-like mode, as defined in Adobe's reference code
            TC_MODE_SATANDVALBLENDING_BW  // Modify the Saturation and Value channel
        };

        std::vector<double> beforeCurve;
        eTCModeId beforeCurveMode;
        std::vector<double> afterCurve;
        eTCModeId afterCurveMode;
        Glib::ustring algo;

        std::vector<double> luminanceCurve;
        bool autoc;
        bool enabledcc;
        bool enabled;
        Glib::ustring filter;
        Glib::ustring setting;
        Glib::ustring method;
        int mixerRed;
        int mixerOrange;
        int mixerYellow;
        int mixerGreen;
        int mixerCyan;
        int mixerBlue;
        int mixerMagenta;
        int mixerPurple;
        int gammaRed;
        int gammaGreen;
        int gammaBlue;
};

/**
  * Parameters of the c/a correction
  */
class CACorrParams {

    public:
        double red;
        double blue;
};

/**
  * Parameters of the highlight recovery
  */
  /*
class HRecParams {

    public:
        bool enabled;
       Glib::ustring method;
};
*/
/**
  * Parameters of the resizing
  */
class ResizeParams {

    public:
        bool enabled;
        double scale;
        Glib::ustring appliesTo;
        Glib::ustring method;
        int dataspec;
        int width;
        int height;
};

/**
  * Parameters of the color spaces used during the processing
  */
class ColorManagementParams {

    public:
        Glib::ustring input;
        bool          toneCurve;
        bool          blendCMSMatrix; // setting no longer used
        int dcpIlluminant;
        Glib::ustring working;
        Glib::ustring output;
        static const Glib::ustring NoICMString;

        Glib::ustring gamma;
        double gampos;
        double slpos;
        bool freegamma;

        ColorManagementParams() {
            setDefaults();
        }
        void setDefaults();
};

/**
  * Typedef for representing a key/value for the exif metadata information
  */
typedef std::map<Glib::ustring, Glib::ustring> ExifPairs;

/**
  * The IPTC key/value pairs
  */
typedef std::map<Glib::ustring, std::vector<Glib::ustring> > IPTCPairs;


class WaveletParams {

    public:
        std::vector<double>   ccwcurve;
        std::vector<double> opacityCurveRG;
        std::vector<double> opacityCurveBY;
        std::vector<double> hhcurve;
        std::vector<double> Chcurve;
        bool enabled;
        bool median;
        bool medianlev;
        bool linkedg;
        bool lipst;
	//	bool edgreinf;
        bool avoid;
        int strength;
        int c[9];
        int ch[9];
		
        Glib::ustring Lmethod;
        Glib::ustring CLmethod;
        Glib::ustring Backmethod;
        Glib::ustring Tilesmethod;
        Glib::ustring choicemethod;
        Glib::ustring CHmethod;
        Glib::ustring Medgreinf;
        Glib::ustring CHSLmethod;
        Glib::ustring EDmethod;
        Glib::ustring Dirmethod;
        Glib::ustring HSmethod;
		int rescon;
		int resconH;
		int reschro;	
		int sup;
		double sky;
		int thres;
		int chroma;	
		int chro;	
		int threshold;
		int threshold2;
		int edgedetect;
		int edgedetectthr;
		int edgedetectthr2;
		int contrast;
		int edgrad;	
		int edgval;	
		int edgthresh;
		int thr;	
		int thrH;	
        double skinprotect;
        Threshold<int> hueskin;
        Threshold<int> hueskin2;
        Threshold<int> hllev;
        Threshold<int> bllev;
        Threshold<int> pastlev;
        Threshold<int> satlev;
        Threshold<int> edgcont;
        Threshold<int> level0noise;
        Threshold<int> level1noise;
        Threshold<int> level2noise;
	
	
		WaveletParams ();
		void setDefaults(); 
        void getCurves(WavCurve &cCurve,WavOpacityCurveRG &opacityCurveLUTRG , WavOpacityCurveBY &opacityCurveLUTBY) const;
        static void getDefaultCCWCurve(std::vector<double> &curve);	
        static void getDefaultOpacityCurveRG(std::vector<double> &curve);
        static void getDefaultOpacityCurveBY(std::vector<double> &curve);
		
};


/**
* Directional pyramid equalizer params
*/
class DirPyrEqualizerParams {

    public:
        bool enabled;
        bool gamutlab;
        double mult[6];
        double threshold;
        double skinprotect;
        Threshold<int> hueskin;
        //Glib::ustring algo;

        DirPyrEqualizerParams() : hueskin(20, 80, 2000, 1200, false) {};
};

/**
 * HSV equalizer params
 */
class HSVEqualizerParams {

    public:
        std::vector<double>   hcurve;
        std::vector<double>   scurve;
        std::vector<double>   vcurve;
};


/**
 *  Film simualtion params
 */
struct FilmSimulationParams {
    bool enabled;
    Glib::ustring clutFilename;
    int strength;
 
    FilmSimulationParams() {
        setDefaults();
    }

    void setDefaults() {
        enabled = false;
        clutFilename = Glib::ustring();
        strength = 100;
    }
};


/**
  * Parameters for RAW demosaicing, common to all sensor type
  */
class RAWParams {

    public:
        /**
         * Parameters for RAW demosaicing specific to Bayer sensors
         */
        class BayerSensor {
            public:
                //enum eMethod{ eahd,hphd,vng4,dcb,amaze,ahd,IGV_noise,fast,
                              //numMethods }; // This MUST be the last enum
                enum eMethod{ amaze, igv, lmmse, eahd, hphd, vng4, dcb, ahd, fast, mono, none,
                              numMethods }; // This MUST be the last enum
                static const char *methodstring[numMethods];

                Glib::ustring method;
                int ccSteps;
                double black0;
                double black1;
                double black2;
                double black3;
                bool twogreen;
                int linenoise;
                int greenthresh;
                int dcb_iterations;
                int lmmse_iterations;
                bool dcb_enhance;
                //bool all_enhance;
        };

        /**
         * Parameters for RAW demosaicing specific to X-Trans sensors
         */
        class XTransSensor {
            public:
            enum eMethod{ threePass, onePass, fast, mono, none,
            numMethods }; // This MUST be the last enum
            static const char *methodstring[numMethods];

            Glib::ustring method;
            int ccSteps;
            double blackred;
            double blackgreen;
            double blackblue;
        };

        BayerSensor bayersensor;         ///< RAW parameters for Bayer sensors
        XTransSensor xtranssensor;       ///< RAW parameters for X-Trans sensors

        enum eFlatFileBlurType{ /*parametric,*/area_ff,v_ff,h_ff,vh_ff,
                                numFlatFileBlurTypes }; // This MUST be the last enum

        static const char *ff_BlurTypestring[numFlatFileBlurTypes];

        Glib::ustring dark_frame;
        bool df_autoselect;

        Glib::ustring ff_file;
        bool ff_AutoSelect;
        int ff_BlurRadius;
        Glib::ustring ff_BlurType;
        bool ff_AutoClipControl;
        int ff_clipControl;

        bool ca_autocorrect;
        double cared;
        double cablue;

        // exposure before interpolation
        double expos;
        double preser;

        bool hotPixelFilter;
        bool deadPixelFilter;
        int hotdeadpix_thresh;

        RAWParams() {
            setDefaults();
        }
        void setDefaults();
};

/**
  * This class holds all the processing parameters applied on the images
  */
class ProcParams {

    public:
        ToneCurveParams         toneCurve;       ///< Tone curve parameters
        LCurveParams            labCurve;        ///< CIELAB luminance curve parameters
        RGBCurvesParams         rgbCurves;       ///< RGB curves parameters
        ColorToningParams       colorToning;     ///< Color Toning parameters
        SharpeningParams        sharpening;      ///< Sharpening parameters
        SharpenEdgeParams       sharpenEdge;     ///< Sharpen edge parameters
        SharpenMicroParams      sharpenMicro;    ///< Sharpen microcontrast parameters
        VibranceParams          vibrance;        ///< Vibrance parameters
        //ColorBoostParams        colorBoost;      ///< Color boost parameters
        WBParams                wb;              ///< White balance parameters
        ColorAppearanceParams   colorappearance;
        //ColorShiftParams        colorShift;      ///< Color shift parameters
        //LumaDenoiseParams       lumaDenoise;     ///< Luminance denoising parameters
        //ColorDenoiseParams      colorDenoise;    ///< Color denoising parameters
        DefringeParams          defringe;        ///< Defringing parameters
        ImpulseDenoiseParams    impulseDenoise;  ///< Impulse denoising parameters
        DirPyrDenoiseParams     dirpyrDenoise;   ///< Directional Pyramid denoising parameters
        EPDParams               epd;             ///< Edge Preserving Decomposition parameters
        SHParams                sh;              ///< Shadow/highlight enhancement parameters
        CropParams              crop;            ///< Crop parameters
        CoarseTransformParams   coarse;          ///< Coarse transformation (90, 180, 270 deg rotation, h/v flipping) parameters
        CommonTransformParams	commonTrans;     ///< Common transformation parameters (autofill)
        RotateParams            rotate;          ///< Rotation parameters
        DistortionParams        distortion;      ///< Lens distortion correction parameters
        LensProfParams          lensProf;        ///< Lens correction profile parameters
        PerspectiveParams       perspective;     ///< Perspective correction parameters
        GradientParams          gradient;        ///< Gradient filter parameters
        PCVignetteParams        pcvignette;      ///< Post-crop vignette filter parameters
        CACorrParams            cacorrection;    ///< Lens c/a correction parameters
        VignettingParams        vignetting;      ///< Lens vignetting correction parameters
        ChannelMixerParams      chmixer;         ///< Channel mixer parameters
        BlackWhiteParams        blackwhite;      ///< Black & White parameters
        ResizeParams            resize;          ///< Resize parameters
        ColorManagementParams   icm;             ///< profiles/color spaces used during the image processing
        RAWParams               raw;             ///< RAW parameters before demosaicing
        WaveletParams        	wavelet;       ///< wavelet wavelet parameters
        DirPyrEqualizerParams   dirpyrequalizer; ///< directional pyramid wavelet parameters
        HSVEqualizerParams      hsvequalizer;    ///< hsv wavelet parameters
        FilmSimulationParams    filmSimulation;  ///< film simulation parameters
        char                    rank;            ///< Custom image quality ranking
        char                    colorlabel;      ///< Custom color label
        bool                    inTrash;         ///< Marks deleted image
        Glib::ustring           appVersion;      ///< Version of the application that generated the parameters
        int                     ppVersion;       ///< Version of the PP file from which the parameters have been read

        ExifPairs                exif;            ///< List of modifications appplied on the exif tags of the input image
        IPTCPairs                iptc;            ///< The IPTC tags and values to be saved to the output image

      /**
        * The constructor only sets the hand-wired defaults.
        */
        ProcParams          ();
      /**
        * Sets the hand-wired defaults parameters.
        */
        void    setDefaults ();
      /**
        * Saves the parameters to possibly two files. This is a performance improvement if a function has to
        * save the same file in two different location, i.e. the cache and the image's directory
        * @param fname   the name of the first file (can be an empty string)
        * @param fname2  the name of the second file (can be an empty string) (optional)
        * @param fnameAbsolute set to false if embedded filenames (if any, darkframe/flatfield) should be stored as relative
        * filenames if they are inside the same directory or in a sub-directory to fname's directory.
        * @param pedited pointer to a ParamsEdited object (optional) to store which values has to be saved
        * @return Error code (=0 if all supplied filenames where created correctly)
        */
        int     save        (Glib::ustring fname, Glib::ustring fname2 = "", bool fnameAbsolute = true, ParamsEdited* pedited=NULL);
      /**
        * Loads the parameters from a file.
        * @param fname the name of the file
        * @params pedited pointer to a ParamsEdited object (optional) to store which values has been loaded
        * @return Error code (=0 if no error)
        */
        int     load        (Glib::ustring fname, ParamsEdited* pedited=NULL);

      /** Creates a new instance of ProcParams.
        * @return a pointer to the new ProcParams instance. */
        static ProcParams* create  ();

      /** Destroys an instance of ProcParams.
        * @param pp a pointer to the ProcParams instance to destroy. */
        static void        destroy (ProcParams* pp);

        static void init ();
        static void cleanup ();

        bool operator== (const ProcParams& other);
        bool operator!= (const ProcParams& other);

    private:
        /** Write the ProcParams's text in the file of the given name.
        * @param fname the name of the file
        * @param content the text to write
        * @return Error code (=0 if no error)
        * */
        int write (Glib::ustring &fname, Glib::ustring &content) const;

};

/**
  * This class associate a ProcParams object and a ParamEdited object through a pointer
  * to instance of each type in order to handle partial pp3 file loading (and later maybe
  * saving too)
  *
  * PartialProfile is not responsible of ProcParams and ParamsEdited object creation
  * and hence is not responsible of their destructions. The function that instanciate
  * PartialProfile object has to handle all this itself.
  */
class PartialProfile {
    public:
        rtengine::procparams::ProcParams* pparams;
        ParamsEdited* pedited;
        PartialProfile& operator=(PartialProfile& rhs) { pparams=rhs.pparams; pedited=rhs.pedited; return *this; };

        PartialProfile      (bool createInstance=false, bool paramsEditedValue=false);
        PartialProfile      (ProcParams* pp, ParamsEdited* pe=NULL, bool fullCopy=false);
        PartialProfile      (const ProcParams* pp, const ParamsEdited* pe=NULL);
        void deleteInstance ();
        void clearGeneral   ();
        int  load           (Glib::ustring fName);
        void set            (bool v);
        const void applyTo  (ProcParams *destParams) const ;
};

/**
  * This class automatically create the pparams and pedited instance in the constructor,
  * and automatically delete them in the destructor. This class has been mostly created
  * to be used with vectors, which use the default constructor/destructor
  */
class AutoPartialProfile : public PartialProfile {
    public:
        AutoPartialProfile() : PartialProfile(true) {}
        ~AutoPartialProfile() { deleteInstance(); }
};

}
}
#endif
