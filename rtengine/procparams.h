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

#include <cmath>
#include <cstdio>
#include <map>
#include <type_traits>
#include <vector>

#include <glibmm/ustring.h>
#include <lcms2.h>

#include "noncopyable.h"

struct ParamsEdited;

namespace rtengine
{

class ColorGradientCurve;
class NoiseCurve;
class OpacityCurve;
class RetinexgaintransmissionCurve;
class RetinextransmissionCurve;
class WavCurve;
class Wavblcurve;
class WavOpacityCurveBY;
class WavOpacityCurveSH;
class WavOpacityCurveRG;
class WavOpacityCurveW;
class WavOpacityCurveWL;
class LocretigainCurve;
class LocretigainCurverab;
class LocLHCurve;
class LocHHCurve;
class LocLLmaskCurve;
class LocCCmaskCurve;
class LocHHmaskCurve;
class LocLLmaskexpCurve;
class LocCCmaskexpCurve;
class LocHHmaskexpCurve;

enum RenderingIntent : int {
    RI_PERCEPTUAL = INTENT_PERCEPTUAL,
    RI_RELATIVE = INTENT_RELATIVE_COLORIMETRIC,
    RI_SATURATION = INTENT_SATURATION,
    RI_ABSOLUTE = INTENT_ABSOLUTE_COLORIMETRIC,
    RI__COUNT
};

namespace procparams
{

template<typename T>
class Threshold final
{
public:
    Threshold(T _bottom, T _top, bool _start_at_one) :
        Threshold(_bottom, _top, 0, 0, _start_at_one, false)
    {
    }

    Threshold(T _bottom_left, T _top_left, T _bottom_right, T _top_right, bool _start_at_one) :
        Threshold(_bottom_left, _top_left, _bottom_right, _top_right, _start_at_one, true)
    {
    }

    template<typename U = T>
    typename std::enable_if<std::is_floating_point<U>::value, bool>::type operator ==(const Threshold<U>& rhs) const
    {
        if (is_double) {
            return
                std::fabs(bottom_left - rhs.bottom_left) < 1e-10
                && std::fabs(top_left - rhs.top_left) < 1e-10
                && std::fabs(bottom_right - rhs.bottom_right) < 1e-10
                && std::fabs(top_right - rhs.top_right) < 1e-10;
        } else {
            return
                std::fabs(bottom_left - rhs.bottom_left) < 1e-10
                && std::fabs(top_left - rhs.top_left) < 1e-10;
        }
    }

    template<typename U = T>
    typename std::enable_if<std::is_integral<U>::value, bool>::type operator ==(const Threshold<U>& rhs) const
    {
        if (is_double) {
            return
                bottom_left == rhs.bottom_left
                && top_left == rhs.top_left
                && bottom_right == rhs.bottom_right
                && top_right == rhs.top_right;
        } else {
            return
                bottom_left == rhs.bottom_left
                && top_left == rhs.top_left;
        }
    }

    template<typename U = T>
    typename std::enable_if<std::is_integral<U>::value, bool>::type operator !=(const Threshold<U>& rhs) const
    {
        return !(*this == rhs);
    }

    T getBottom() const
    {
        return bottom_left;
    }

    T getTop() const
    {
        return top_left;
    }

    T getBottomLeft() const
    {
        return bottom_left;
    }

    T getTopLeft() const
    {
        return top_left;
    }

    T getBottomRight() const
    {
        return bottom_right;
    }

    T getTopRight() const
    {
        return top_right;
    }

    void setValues(T bottom, T top)
    {
        bottom_left = bottom;
        top_left = top;
    }

    void setValues(T bottom_left, T top_left, T bottom_right, T top_right)
    {
        this->bottom_left = bottom_left;
        this->top_left = top_left;
        this->bottom_right = bottom_right;
        this->top_right = top_right;
    }

    bool isDouble() const
    {
        return is_double;
    }

    std::vector<T> toVector() const
    {
        if (is_double) {
            return {
                bottom_left,
                top_left,
                bottom_right,
                top_right
            };
        } else {
            return {
                bottom_left,
                top_left
            };
        }
    }

    // RT: Type of the returned value
    // RV: Type of the value on the X axis
    // RV2: Type of the maximum value on the Y axis
    template <typename RT, typename RV, typename RV2>
    RT multiply(RV x, RV2 y_max) const
    {
        const double val = x;

        if (init_eql) {
            if (is_double) {
                if (val == static_cast<double>(bottom_right) && static_cast<double>(bottom_right) == static_cast<double>(top_right)) {
                    // This handles the special case where the 2 right values are the same, then bottom one is sent back,
                    // useful if one wants to keep the bottom value even beyond the x max bound
                    return 0;
                }

                if (val >= static_cast<double>(top_right)) {
                    return y_max;
                }

                if (val > static_cast<double>(bottom_right)) {
                    return static_cast<double>(y_max * (val - static_cast<double>(bottom_right)) / (static_cast<double>(top_right) - static_cast<double>(bottom_right)));
                }
            }

            if (val >= static_cast<double>(bottom_left)) {
                return 0;
            }

            if (val > static_cast<double>(top_left)) {
                return static_cast<double>(y_max * (1. - (val - static_cast<double>(bottom_left)) / (static_cast<double>(top_left) - static_cast<double>(bottom_left))));
            }

            return y_max;
        } else {
            if (is_double) {
                if (val == static_cast<double>(bottom_right) && static_cast<double>(bottom_right) == static_cast<double>(top_right)) {
                    // This handles the special case where the 2 right values are the same, then top one is sent back,
                    // useful if one wants to keep the top value even beyond the x max bound
                    return y_max;
                }

                if (val >= static_cast<double>(bottom_right)) {
                    return 0;
                }

                if (val > static_cast<double>(top_right)) {
                    return static_cast<double>(y_max * (1.0 - (val - static_cast<double>(top_right)) / (static_cast<double>(bottom_right) - static_cast<double>(top_right))));
                }
            }

            if (val >= static_cast<double>(top_left)) {
                return y_max;
            }

            if (val > static_cast<double>(bottom_left)) {
                return static_cast<double>(y_max * (val - static_cast<double>(bottom_left)) / (static_cast<double>(top_left) - static_cast<double>(bottom_left)));
            }

            return 0;
        }
    }

private:
    Threshold(T _bottom_left, T _top_left, T _bottom_right, T _top_right, bool _start_at_one, bool _is_double) :
        bottom_left(_bottom_left),
        top_left(_top_left),
        bottom_right(_bottom_right),
        top_right(_top_right),
        init_eql(_start_at_one),
        is_double(_is_double)
    {
    }

    T bottom_left;
    T top_left;
    T bottom_right;
    T top_right;
    bool init_eql;
    bool is_double;
};

enum class ToneCurveMode : int {
    STD,               // Standard modes, the curve is applied on all component individually
    WEIGHTEDSTD,       // Weighted standard mode
    FILMLIKE,          // Film-like mode, as defined in Adobe's reference code
    SATANDVALBLENDING, // Modify the Saturation and Value channel
    LUMINANCE,         // Modify the Luminance channel with coefficients from Rec 709's
    PERCEPTUAL         // Keep color appearance constant using perceptual modeling
};

/**
  * Parameters of the tone curve
  */
struct ToneCurveParams {
    bool autoexp;
    double clip;
    bool hrenabled; // Highlight Reconstruction enabled
    Glib::ustring method; // Highlight Reconstruction's method
    double expcomp;
    std::vector<double> curve;
    std::vector<double> curve2;
    ToneCurveMode curveMode;
    ToneCurveMode curveMode2;
    int brightness;
    int black;
    int contrast;
    int saturation;
    int shcompr;
    int hlcompr; // Highlight Recovery's compression
    int hlcomprthresh; // Highlight Recovery's threshold
    bool histmatching; // histogram matching
    bool fromHistMatching;
    bool clampOOG; // clamp out of gamut colours

    ToneCurveParams();

    bool isPanningRelatedChange(const ToneCurveParams& other) const;
    bool operator ==(const ToneCurveParams& other) const;
    bool operator !=(const ToneCurveParams& other) const;

};

/**
  * Parameters of Retinex
  */
struct RetinexParams {
    bool enabled;
    std::vector<double>   cdcurve;
    std::vector<double>   cdHcurve;
    std::vector<double>   lhcurve;
    std::vector<double> transmissionCurve;
    std::vector<double> gaintransmissionCurve;
    std::vector<double>   mapcurve;
    int     str;
    int     scal;
    int     iter;
    int     grad;
    int     grads;
    double  gam;
    double  slope;
    int     neigh;
    int     offs;
    int     highlights;
    int     htonalwidth;
    int     shadows;
    int     stonalwidth;
    int     radius;

    Glib::ustring retinexMethod;
    Glib::ustring retinexcolorspace;
    Glib::ustring gammaretinex;
    Glib::ustring mapMethod;
    Glib::ustring viewMethod;
    int     vart;
    int     limd;
    int     highl;
    int     skal;
    bool    medianmap;

    RetinexParams();

    bool operator ==(const RetinexParams& other) const;
    bool operator !=(const RetinexParams& other) const;

    void getCurves(RetinextransmissionCurve& transmissionCurveLUT, RetinexgaintransmissionCurve& gaintransmissionCurveLUT) const;
};


/**
  * Parameters of the luminance curve
  */
struct LCurveParams {
    bool enabled;
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

    LCurveParams();

    bool operator ==(const LCurveParams& other) const;
    bool operator !=(const LCurveParams& other) const;
};


/**
 * Parameters for local contrast
 */
struct LocalContrastParams {
    bool enabled;
    int radius;
    double amount;
    double darkness;
    double lightness;

    LocalContrastParams();

    bool operator==(const LocalContrastParams &other) const;
    bool operator!=(const LocalContrastParams &other) const;
};


/**
  * Parameters of the RGB curves
  */
struct RGBCurvesParams {
    bool enabled;
    bool lumamode;
    std::vector<double>   rcurve;
    std::vector<double>   gcurve;
    std::vector<double>   bcurve;

    RGBCurvesParams();

    bool operator ==(const RGBCurvesParams& other) const;
    bool operator !=(const RGBCurvesParams& other) const;
};

/**
  * Parameters of the Color Toning
  */
struct ColorToningParams {
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
     *  LabGrid    :
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

    double labgridALow;
    double labgridBLow;
    double labgridAHigh;
    double labgridBHigh;
    static const double LABGRID_CORR_MAX;
    static const double LABGRID_CORR_SCALE;

    struct LabCorrectionRegion {
        enum { CHAN_ALL = -1, CHAN_R, CHAN_G, CHAN_B };
        double a;
        double b;
        double saturation;
        double slope;
        double offset;
        double power;
        std::vector<double> hueMask;
        std::vector<double> chromaticityMask;
        std::vector<double> lightnessMask;
        double maskBlur;
        int channel;

        LabCorrectionRegion();
        bool operator==(const LabCorrectionRegion &other) const;
        bool operator!=(const LabCorrectionRegion &other) const;
    };
    std::vector<LabCorrectionRegion> labregions;
    int labregionsShowMask;

    ColorToningParams();

    bool operator ==(const ColorToningParams& other) const;
    bool operator !=(const ColorToningParams& other) const;

    /// @brief Transform the mixer values to their curve equivalences
    void mixerToCurve(std::vector<double>& colorCurve, std::vector<double>& opacityCurve) const;
    /// @brief Specifically transform the sliders values to their curve equivalences
    void slidersToCurve(std::vector<double>& colorCurve, std::vector<double>& opacityCurve) const;
    /// @brief Fill the ColorGradientCurve and OpacityCurve LUTf from the control points curve or sliders value
    void getCurves(ColorGradientCurve& colorCurveLUT, OpacityCurve& opacityCurveLUT, const double xyz_rgb[3][3], bool& opautili) const;
};


/**
  * Parameters of the sharpening
  */
struct SharpeningParams {
    bool           enabled;
    double         contrast;
    bool           autoContrast;
    double         blurradius;
    double         gamma;
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

    SharpeningParams();

    bool operator ==(const SharpeningParams& other) const;
    bool operator !=(const SharpeningParams& other) const;
};

struct SharpenEdgeParams {
    bool    enabled;
    int     passes;
    double  amount;
    bool    threechannels;

    SharpenEdgeParams();

    bool operator ==(const SharpenEdgeParams& other) const;
    bool operator !=(const SharpenEdgeParams& other) const;

};


struct SharpenMicroParams {
    bool    enabled;
    bool    matrix;
    double  amount;
    double  contrast;
    int     uniformity;

    SharpenMicroParams();

    bool operator ==(const SharpenMicroParams& other) const;
    bool operator !=(const SharpenMicroParams& other) const;
};

struct CaptureSharpeningParams {
    bool           enabled;
    bool           autoContrast;
    bool           autoRadius;
    double         contrast;
    double         deconvradius;
    double         deconvradiusOffset;
    int            deconviter;
    bool           deconvitercheck;

    CaptureSharpeningParams();

    bool operator ==(const CaptureSharpeningParams& other) const;
    bool operator !=(const CaptureSharpeningParams& other) const;
};

/**
  * Parameters of the vibrance
  */
struct VibranceParams {
    bool           enabled;
    int            pastels;
    int            saturated;
    Threshold<int> psthreshold;
    bool           protectskins;
    bool           avoidcolorshift;
    bool           pastsattog;
    std::vector<double> skintonescurve;

    VibranceParams();

    bool operator ==(const VibranceParams& other) const;
    bool operator !=(const VibranceParams& other) const;
};

/**
  * Parameters of the white balance adjustments
  */
struct WBEntry {
    enum class Type {
        CAMERA,
        AUTO,
        DAYLIGHT,
        CLOUDY,
        SHADE,
        WATER,
        TUNGSTEN,
        FLUORESCENT,
        LAMP,
        FLASH,
        LED,
        // CUSTOM one must remain the last one!
        CUSTOM
    };

    Glib::ustring ppLabel;
    Type type;
    Glib::ustring GUILabel;
    int temperature;
    double green;
    double equal;
    double tempBias;
};

struct WBParams {
    bool enabled;
    Glib::ustring   method;
    int             temperature;
    double          green;
    double          equal;
    double          tempBias;

    WBParams();

    bool isPanningRelatedChange(const WBParams& other) const;
    bool operator ==(const WBParams& other) const;
    bool operator !=(const WBParams& other) const;

    static const std::vector<WBEntry>& getWbEntries();
};

/**
 * Parameters of colorappearance
 */
struct ColorAppearanceParams {
    enum class TcMode {
        LIGHT,    // Lightness mode
        BRIGHT,   // Brightness mode
    };

    enum class CtcMode {
        CHROMA,   // chroma mode
        SATUR,    // saturation mode
        COLORF,   // colorfullness mode
    };

    bool          enabled;
    int           degree;
    bool          autodegree;
    int           degreeout;
    bool          autodegreeout;
    std::vector<double> curve;
    std::vector<double> curve2;
    std::vector<double> curve3;
    TcMode     curveMode;
    TcMode     curveMode2;
    CtcMode    curveMode3;

    Glib::ustring surround;
    Glib::ustring surrsrc;
    double        adapscen;
    bool          autoadapscen;
    int        ybscen;
    bool          autoybscen;

    double        adaplum;
    int           badpixsl;
    Glib::ustring wbmodel;
    Glib::ustring illum;
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
    bool          datacie;
    bool          tonecie;
    int tempout;
    bool          autotempout;
    int ybout;
    double greenout;
    int tempsc;
    double greensc;
    bool presetcat02;

    ColorAppearanceParams();

    bool operator ==(const ColorAppearanceParams& other) const;
    bool operator !=(const ColorAppearanceParams& other) const;
};

/**
 * Parameters of defringing
 */
struct DefringeParams {
    bool    enabled;
    double  radius;
    int     threshold;
    std::vector<double> huecurve;

    DefringeParams();

    bool operator ==(const DefringeParams& other) const;
    bool operator !=(const DefringeParams& other) const;
};

/**
  * Parameters of impulse denoising
  */
struct ImpulseDenoiseParams {
    bool    enabled;
    int     thresh;

    ImpulseDenoiseParams();

    bool operator ==(const ImpulseDenoiseParams& other) const;
    bool operator !=(const ImpulseDenoiseParams& other) const;
};

/**
 * Parameters of the directional pyramid denoising
 */
struct DirPyrDenoiseParams {
    std::vector<double>   lcurve;
    std::vector<double>   cccurve;

    bool    enabled;
    bool    enhance;
    bool    median;

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

    DirPyrDenoiseParams();

    bool operator ==(const DirPyrDenoiseParams& other) const;
    bool operator !=(const DirPyrDenoiseParams& other) const;

    void getCurves(NoiseCurve& lCurve, NoiseCurve& cCurve) const;
};

// EPD related parameters.
struct EPDParams {
    bool   enabled;
    double strength;
    double gamma;
    double edgeStopping;
    double scale;
    int    reweightingIterates;

    EPDParams();

    bool operator ==(const EPDParams& other) const;
    bool operator !=(const EPDParams& other) const;
};

// Fattal02 Tone-Mapping parameters
struct FattalToneMappingParams {
    bool enabled;
    int threshold;
    int amount;
    int anchor;

    FattalToneMappingParams();

    bool operator ==(const FattalToneMappingParams& other) const;
    bool operator !=(const FattalToneMappingParams& other) const;
};

/**
  * Parameters of the shadow/highlight enhancement
  */
struct SHParams {
    bool    enabled;
    int     highlights;
    int     htonalwidth;
    int     shadows;
    int     stonalwidth;
    int     radius;
    bool    lab;

    SHParams();

    bool operator ==(const SHParams& other) const;
    bool operator !=(const SHParams& other) const;
};

/**
  * Parameters of the cropping
  */
struct CropParams {
    bool    enabled;
    int     x;
    int     y;
    int     w;
    int     h;
    bool    fixratio;
    Glib::ustring   ratio;
    Glib::ustring   orientation;
    Glib::ustring   guide;

    CropParams();

    bool operator ==(const CropParams& other) const;
    bool operator !=(const CropParams& other) const;

    void mapToResized(int resizedWidth, int resizedHeight, int scale, int& x1, int& x2, int& y1, int& y2) const;
};

/**
  * Parameters of the coarse transformations like 90 deg rotations and h/v flipping
  */
struct CoarseTransformParams {
    int     rotate;
    bool    hflip;
    bool    vflip;

    CoarseTransformParams();

    bool operator ==(const CoarseTransformParams& other) const;
    bool operator !=(const CoarseTransformParams& other) const;
};

/**
  * Common transformation parameters
  */
struct CommonTransformParams {
    Glib::ustring method;
    bool autofill;

    CommonTransformParams();

    bool operator ==(const CommonTransformParams& other) const;
    bool operator !=(const CommonTransformParams& other) const;
};

/**
  * Parameters of the rotation
  */
struct RotateParams {
    double  degree;

    RotateParams();

    bool operator ==(const RotateParams& other) const;
    bool operator !=(const RotateParams& other) const;
};

/**
  * Parameters of the distortion correction
  */
struct DistortionParams {
    double  amount;

    DistortionParams();

    bool operator ==(const DistortionParams& other) const;
    bool operator !=(const DistortionParams& other) const;
};

// Lens profile correction parameters
struct LensProfParams {
    enum class LcMode {
        NONE,               // No lens correction
        LENSFUNAUTOMATCH,   // Lens correction using auto matched lensfun database entry
        LENSFUNMANUAL,      // Lens correction using manually selected lensfun database entry
        LCP                 // Lens correction using lcp file
    };

    LcMode lcMode;
    Glib::ustring lcpFile;
    bool useDist, useVign, useCA;
    Glib::ustring lfCameraMake;
    Glib::ustring lfCameraModel;
    Glib::ustring lfLens;

    LensProfParams();

    bool operator ==(const LensProfParams& other) const;
    bool operator !=(const LensProfParams& other) const;

    bool useLensfun() const;
    bool lfAutoMatch() const;
    bool useLcp() const;
    bool lfManual() const;

    const std::vector<const char*>& getMethodStrings() const;
    Glib::ustring getMethodString(LcMode mode) const;
    LcMode getMethodNumber(const Glib::ustring& mode) const;
};


/**
  * Parameters of the perspective correction
  */
struct PerspectiveParams {
    Glib::ustring method;
    double  horizontal;
    double  vertical;
    double  camera_crop_factor;
    double  camera_focal_length;
    double  camera_pitch;
    double  camera_roll;
    double  camera_shift_horiz;
    double  camera_shift_vert;
    double  camera_yaw;
    double  projection_pitch;
    double  projection_rotate;
    double  projection_shift_horiz;
    double  projection_shift_vert;
    double  projection_yaw;

    PerspectiveParams();

    bool operator ==(const PerspectiveParams& other) const;
    bool operator !=(const PerspectiveParams& other) const;
};

/**
  * Parameters of the gradient filter
  */
struct GradientParams {
    bool   enabled;
    double degree;
    int    feather;
    double strength;
    int    centerX;
    int    centerY;

    GradientParams();

    bool operator ==(const GradientParams& other) const;
    bool operator !=(const GradientParams& other) const;
};

/**
  * Parameters of the Local Lab
  */
struct LocallabParams {
    struct LocallabSpot {
        // Control spot settings
        Glib::ustring name;
        bool isvisible;
        Glib::ustring shape; // ELI, RECT
        Glib::ustring spotMethod; // norm, exc
        Glib::ustring wavMethod; // D2, D4, D6, D10, D14
        int sensiexclu;
        int structexclu;
        double struc;
        Glib::ustring shapeMethod; // IND, SYM, INDSL, SYMSL
        std::vector<int> loc; // For ellipse/rectangle: {locX, locXL, locY, locYT}
        int centerX;
        int centerY;
        int circrad;
        Glib::ustring qualityMethod; // none, std, enh, enhsup, contr, sob2
        Glib::ustring complexMethod; // sim, mod, all
        double transit;
        double feather;
        double thresh;
        double iter;
        double balan;
        double balanh;
        double colorde;
        double colorscope;
        double transitweak;
        double transitgrad;
        bool avoid;
        bool blwh;
        bool recurs;
        bool laplac;
        bool deltae;
        bool shortc;
        bool savrest;
        int scopemask;
        int lumask;
        // Color & Light
        bool visicolor;
        bool expcolor;
        int complexcolor;
        bool curvactiv;
        int lightness;
        int contrast;
        int chroma;
        double labgridALow;
        double labgridBLow;
        double labgridAHigh;
        double labgridBHigh;
        double labgridALowmerg;
        double labgridBLowmerg;
        double labgridAHighmerg;
        double labgridBHighmerg;
        int strengthgrid;
        int sensi;
        int structcol;
        double strcol;
        double strcolab;
        double strcolh;
        double angcol;
        int blurcolde;
        double blurcol;
        double contcol;
        int blendmaskcol;
        double radmaskcol;
        double chromaskcol;
        double gammaskcol;
        double slomaskcol;
        int shadmaskcol;
        double strumaskcol;
        double lapmaskcol;
        Glib::ustring qualitycurveMethod; // none, std
        Glib::ustring gridMethod; // one, two
        Glib::ustring merMethod; // mone, mtwo, mthr, mfou, mfiv
        Glib::ustring toneMethod; // one, two, thr, fou
        Glib::ustring mergecolMethod; // one, two, thr, fou, fiv, six, sev, sev0, sev1, sev2, hei, nin, ten, ele, twe, thi, for, hue, sat, col, lum
        std::vector<double> llcurve;
        std::vector<double> lccurve;
        std::vector<double> cccurve;
        std::vector<double> clcurve;
        std::vector<double> rgbcurve;
        std::vector<double> LHcurve;
        std::vector<double> HHcurve;
        bool invers;
        bool special;
        bool toolcol;
        bool enaColorMask;
        bool fftColorMask;
        std::vector<double> CCmaskcurve;
        std::vector<double> LLmaskcurve;
        std::vector<double> HHmaskcurve;
        std::vector<double> HHhmaskcurve;
        double softradiuscol;
        double opacol;
        double mercol;
        double merlucol;
        double conthrcol;
        std::vector<double> Lmaskcurve;
        std::vector<double> LLmaskcolcurvewav;
        Threshold<int> csthresholdcol;
        // Exposure
        bool visiexpose;
        bool expexpose;
        int complexexpose;
        double expcomp;
        int hlcompr;
        int hlcomprthresh;
        int black;
        int shadex;
        int shcompr;
        int expchroma;
        int sensiex;
        int structexp;
        int blurexpde;
        double strexp;
        double angexp;
        std::vector<double> excurve;
        bool inversex;
        bool enaExpMask;
        bool enaExpMaskaft;
        std::vector<double> CCmaskexpcurve;
        std::vector<double> LLmaskexpcurve;
        std::vector<double> HHmaskexpcurve;
        int blendmaskexp;
        double radmaskexp;
        double chromaskexp;
        double gammaskexp;
        double slomaskexp;
        double lapmaskexp;
        double strmaskexp;
        double angmaskexp;
        double softradiusexp;
        std::vector<double> Lmaskexpcurve;
        Glib::ustring expMethod; // std, pde
        Glib::ustring exnoiseMethod; // none, med, medhi
        double laplacexp;
        double balanexp;
        double linear;
        double gamm;
        double fatamount;
        double fatdetail;
        double fatanchor;
        double fatlevel;
        // Shadow highlight
        bool visishadhigh;
        bool expshadhigh;
        int complexshadhigh;
        Glib::ustring shMethod; // std, tone
        int multsh[5];
        int highlights;
        int h_tonalwidth;
        int shadows;
        int s_tonalwidth;
        int sh_radius;
        int sensihs;
        bool enaSHMask;
        std::vector<double> CCmaskSHcurve;
        std::vector<double> LLmaskSHcurve;
        std::vector<double> HHmaskSHcurve;
        int blendmaskSH;
        double radmaskSH;
        int blurSHde;
        double strSH;
        double angSH;
        bool inverssh;
        double chromaskSH;
        double gammaskSH;
        double slomaskSH;
        double lapmaskSH;
        int detailSH;
        std::vector<double> LmaskSHcurve;
        double fatamountSH;
        double fatanchorSH;
        double gamSH;
        double sloSH;
        // Vibrance
        bool visivibrance;
        bool expvibrance;
        int complexvibrance;
        int saturated;
        int pastels;
        int warm;
        Threshold<int> psthreshold;
        bool protectskins;
        bool avoidcolorshift;
        bool pastsattog;
        int sensiv;
        std::vector<double> skintonescurve;
        std::vector<double> CCmaskvibcurve;
        std::vector<double> LLmaskvibcurve;
        std::vector<double> HHmaskvibcurve;
        bool enavibMask;
        int blendmaskvib;
        double radmaskvib;
        double chromaskvib;
        double gammaskvib;
        double slomaskvib;
        double lapmaskvib;
        double strvib;
        double strvibab;
        double strvibh;
        double angvib;
        std::vector<double> Lmaskvibcurve;
        // Soft Light
        bool visisoft;
        bool expsoft;
        int complexsoft;
        int streng;
        int sensisf;
        double laplace;
        Glib::ustring softMethod; // soft, reti
        // Blur & Noise
        bool visiblur;
        bool expblur;
        int complexblur;
        double radius;
        int strength;
        int sensibn;
        int itera;
        int guidbl;
        int strbl;
        int isogr;
        int strengr;
        int scalegr;
        int epsbl;
        Glib::ustring blMethod; // blur, med, guid
        Glib::ustring chroMethod; // lum, chr, all
        Glib::ustring blurMethod; // norm, inv
        Glib::ustring medMethod; // none, 33, 55, 77, 99
        bool activlum;
        double noiselumf;
        double noiselumf0;
        double noiselumf2;
        double noiselumc;
        double noiselumdetail;
        int noiselequal;
        double noisechrof;
        double noisechroc;
        double noisechrodetail;
        int adjblur;
        int bilateral;
        int sensiden;
        int detailthr;
        std::vector<double> locwavcurveden;
        Glib::ustring showmaskblMethodtyp;
        std::vector<double> CCmaskblcurve;
        std::vector<double> LLmaskblcurve;
        std::vector<double> HHmaskblcurve;
        bool enablMask;
        bool fftwbl;
        bool toolbl;
        int blendmaskbl;
        double radmaskbl;
        double chromaskbl;
        double gammaskbl;
        double slomaskbl;
        double lapmaskbl;
        int shadmaskbl;
        int shadmaskblsha;
        double strumaskbl;
        std::vector<double> Lmaskblcurve;
        std::vector<double> LLmaskblcurvewav;
        Threshold<int> csthresholdblur;
        // Tone Mapping
        bool visitonemap;
        bool exptonemap;
        int complextonemap;
        double stren;
        double gamma;
        double estop;
        double scaltm;
        int rewei;
        double satur;
        int sensitm;
        double softradiustm;
        double amount;
        bool equiltm;
        std::vector<double> CCmasktmcurve;
        std::vector<double> LLmasktmcurve;
        std::vector<double> HHmasktmcurve;
        bool enatmMask;
        bool enatmMaskaft;
        int blendmasktm;
        double radmasktm;
        double chromasktm;
        double gammasktm;
        double slomasktm;
        double lapmasktm;
        std::vector<double> Lmasktmcurve;
        // Retinex
        bool visireti;
        bool expreti;
        int complexreti;
        Glib::ustring retinexMethod; // low, uni, high
        double str;
        double chrrt;
        double neigh;
        double vart;
        double offs;
        int dehaz;
        int depth;
        int sensih;
        std::vector<double> localTgaincurve;
        std::vector<double> localTtranscurve;
        bool inversret;
        bool equilret;
        bool loglin;
        bool lumonly;
        double softradiusret;
        std::vector<double> CCmaskreticurve;
        std::vector<double> LLmaskreticurve;
        std::vector<double> HHmaskreticurve;
        bool enaretiMask;
        bool enaretiMasktmap;
        int blendmaskreti;
        double radmaskreti;
        double chromaskreti;
        double gammaskreti;
        double slomaskreti;
        double lapmaskreti;
        double scalereti;
        double darkness;
        double lightnessreti;
        double limd;
        double cliptm;
        bool fftwreti;
        std::vector<double> Lmaskreticurve;
        // Sharpening
        bool visisharp;
        bool expsharp;
        int complexsharp;
        int sharcontrast;
        double sharradius;
        int sharamount;
        int shardamping;
        int shariter;
        double sharblur;
        int sensisha;
        bool inverssha;
        // Local Contrast
        bool visicontrast;
        bool expcontrast;
        int complexcontrast;
        int lcradius;
        double lcamount;
        double lcdarkness;
        double lclightness;
        double sigmalc;
        int levelwav;
        double residcont;
        double residsha;
        double residshathr;
        double residhi;
        double residhithr;
        double residblur;
        double levelblur;
        double sigmabl;
        double residchro;
        double residcomp;
        double sigma;
        double offset;
        double sigmadr;
        double threswav;
        double chromalev;
        double chromablu;
        double sigmadc;
        double deltad;
        double fatres;
        double clarilres;
        double claricres;
        double clarisoft;
        double sigmalc2;
        double strwav;
        double angwav;
        double strengthw;
        double sigmaed;
        double radiusw;
        double detailw;
        double gradw;
        double tloww;
        double thigw;
        double edgw;
        double basew;
        int sensilc;
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
        Glib::ustring localcontMethod; // loc, wav
        Glib::ustring localedgMethod; // fir, sec, thr
        Glib::ustring localneiMethod; // none, low, high
        std::vector<double> locwavcurve;
        Threshold<int> csthreshold;
        std::vector<double> loclevwavcurve;
        std::vector<double> locconwavcurve;
        std::vector<double> loccompwavcurve;
        std::vector<double> loccomprewavcurve;
        std::vector<double> locedgwavcurve;
        std::vector<double> CCmasklccurve;
        std::vector<double> LLmasklccurve;
        std::vector<double> HHmasklccurve;
        bool enalcMask;
        int blendmasklc;
        double radmasklc;
        double chromasklc;
        std::vector<double> Lmasklccurve;
        // Contrast by detail levels
        bool visicbdl;
        bool expcbdl;
        int complexcbdl;
        double mult[6];
        double chromacbdl;
        double threshold;
        int sensicb;
        double clarityml;
        int contresid;
        double blurcbdl;
        double softradiuscb;
        bool enacbMask;
        std::vector<double> CCmaskcbcurve;
        std::vector<double> LLmaskcbcurve;
        std::vector<double> HHmaskcbcurve;
        int blendmaskcb;
        double radmaskcb;
        double chromaskcb;
        double gammaskcb;
        double slomaskcb;
        double lapmaskcb;
        std::vector<double> Lmaskcbcurve;
        // Log encoding
        bool visilog;
        bool explog;
        bool autocompute;
        double sourceGray;
        double targetGray;
        bool Autogray;
        bool fullimage;
        double blackEv;
        double whiteEv;
        double detail;
        int sensilog;
        double baselog;
        double strlog;
        double anglog;

        LocallabSpot();

        bool operator ==(const LocallabSpot& other) const;
        bool operator !=(const LocallabSpot& other) const;
    };

    static const double LABGRIDL_CORR_MAX;
    static const double LABGRIDL_CORR_SCALE;
    static const double LABGRIDL_DIRECT_SCALE;

    bool enabled;
    int selspot;
    std::vector<LocallabSpot> spots;

    LocallabParams();

    bool operator ==(const LocallabParams& other) const;
    bool operator !=(const LocallabParams& other) const;
};

/**
  * Parameters of the post-crop vignette filter
  */
struct PCVignetteParams {
    bool   enabled;
    double strength;
    int    feather;
    int    roundness;

    PCVignetteParams();

    bool operator ==(const PCVignetteParams& other) const;
    bool operator !=(const PCVignetteParams& other) const;
};

/**
  * Parameters of the vignetting correction
  */
struct VignettingParams {
    int  amount;
    int  radius;
    int  strength;
    int  centerX;
    int  centerY;

    VignettingParams();

    bool operator ==(const VignettingParams& other) const;
    bool operator !=(const VignettingParams& other) const;
};

/**
  * Parameters of the color mixer
  */
struct ChannelMixerParams {
    bool enabled;
    int red[3];
    int green[3];
    int blue[3];

    ChannelMixerParams();

    bool operator ==(const ChannelMixerParams& other) const;
    bool operator !=(const ChannelMixerParams& other) const;
};

struct BlackWhiteParams {
    enum class TcMode {
        STD_BW,               // Standard modes, the curve is applied on all component individually
        WEIGHTEDSTD_BW,       // Weighted standard mode
        FILMLIKE_BW,          // Film-like mode, as defined in Adobe's reference code
        SATANDVALBLENDING_BW  // Modify the Saturation and Value channel
    };

    std::vector<double> beforeCurve;
    TcMode beforeCurveMode;
    std::vector<double> afterCurve;
    TcMode afterCurveMode;
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

    BlackWhiteParams();

    bool operator ==(const BlackWhiteParams& other) const;
    bool operator !=(const BlackWhiteParams& other) const;
};

/**
  * Parameters of the c/a correction
  */
struct CACorrParams {
    double red;
    double blue;

    CACorrParams();

    bool operator ==(const CACorrParams& other) const;
    bool operator !=(const CACorrParams& other) const;
};

/**
  * Parameters of the resizing
  */
struct ResizeParams {
    bool enabled;
    double scale;
    Glib::ustring appliesTo;
    Glib::ustring method;
    int dataspec;
    int width;
    int height;
    bool allowUpscaling;

    ResizeParams();

    bool operator ==(const ResizeParams& other) const;
    bool operator !=(const ResizeParams& other) const;
};

/**
  * Parameters of the color spaces used during the processing
  */
struct ColorManagementParams {
    Glib::ustring inputProfile;
    bool toneCurve;
    bool applyLookTable;
    bool applyBaselineExposureOffset;
    bool applyHueSatMap;
    int dcpIlluminant;

    Glib::ustring workingProfile;
    Glib::ustring workingTRC;
    double workingTRCGamma;
    double workingTRCSlope;

    Glib::ustring outputProfile;
    RenderingIntent outputIntent;
    bool outputBPC;

    static const Glib::ustring NoICMString;

    ColorManagementParams();

    bool operator ==(const ColorManagementParams& other) const;
    bool operator !=(const ColorManagementParams& other) const;
};

/**
  * Parameters for metadata handling
  */
struct MetaDataParams {
    enum Mode {
        TUNNEL,
        EDIT,
        STRIP
    };
    Mode mode;

    MetaDataParams();

    bool operator ==(const MetaDataParams &other) const;
    bool operator !=(const MetaDataParams &other) const;
};


/**
  * Minimal wrapper allowing forward declaration for representing a key/value for the exif metadata information
  */
class ExifPairs final
{
public:
    using const_iterator = std::map<Glib::ustring, Glib::ustring>::const_iterator;

    const_iterator begin() const
    {
        return pairs.begin();
    }

    const_iterator end() const
    {
        return pairs.end();
    }

    void clear()
    {
        pairs.clear();
    }

    Glib::ustring& operator[](const Glib::ustring& key)
    {
        return pairs[key];
    }

    bool operator ==(const ExifPairs& other) const
    {
        return pairs == other.pairs;
    }

private:
    std::map<Glib::ustring, Glib::ustring> pairs;
};

/**
  * The IPTC key/value pairs
  */
class IPTCPairs final
{
public:
    using iterator = std::map<Glib::ustring, std::vector<Glib::ustring>>::iterator;
    using const_iterator = std::map<Glib::ustring, std::vector<Glib::ustring>>::const_iterator;

    iterator find(const Glib::ustring& key)
    {
        return pairs.find(key);
    }

    const_iterator begin() const
    {
        return pairs.begin();
    }

    const_iterator end() const
    {
        return pairs.end();
    }

    bool empty() const
    {
        return pairs.empty();
    }

    void clear()
    {
        pairs.clear();
    }

    std::vector<Glib::ustring>& operator[](const Glib::ustring& key)
    {
        return pairs[key];
    }

    bool operator ==(const IPTCPairs& other) const
    {
        return pairs == other.pairs;
    }

private:
    std::map<Glib::ustring, std::vector<Glib::ustring>> pairs;
};

struct WaveletParams {
    std::vector<double> ccwcurve;
    std::vector<double> blcurve;
    std::vector<double> levelshc;
    std::vector<double> opacityCurveRG;
    std::vector<double> opacityCurveSH;
    std::vector<double> opacityCurveBY;
    std::vector<double> opacityCurveW;
    std::vector<double> opacityCurveWL;
    std::vector<double> hhcurve;
    std::vector<double> Chcurve;
    std::vector<double> wavclCurve;
    bool enabled;
    bool median;
    bool medianlev;
    bool linkedg;
    bool cbenab;
    int greenlow;
    int bluelow;
    int greenmed;
    int bluemed;
    int greenhigh;
    int bluehigh;
    double ballum;
    double balchrom;
    double chromfi;
    double chromco;
    double mergeL;
    double mergeC;
    double softrad;
    double softradend;

    bool lipst;
    bool avoid;
    bool showmask;
    bool oldsh;
    bool tmr;
    int strength;
    int balance;
    double sigmafin;
    double sigmaton;
    double sigmacol;
    double sigmadir;
    double rangeab;
    double protab;
    int iter;
    bool expcontrast;
    bool expchroma;
    int c[9];
    int ch[9];
    bool expedge;
    bool expbl;
    bool expresid;
    bool expfinal;
    bool exptoning;
    bool expnoise;
    bool expclari;
    double labgridALow;
    double labgridBLow;
    double labgridAHigh;
    double labgridBHigh;
    static const double LABGRID_CORR_MAX;
    static const double LABGRID_CORR_SCALE;
    static const double LABGRIDL_DIRECT_SCALE;
    int Lmethod;
    Glib::ustring CLmethod;
    Glib::ustring Backmethod;
    Glib::ustring Tilesmethod;
    Glib::ustring daubcoeffmethod;
    Glib::ustring CHmethod;
    Glib::ustring Medgreinf;
    Glib::ustring ushamethod;
    Glib::ustring CHSLmethod;
    Glib::ustring EDmethod;
    Glib::ustring NPmethod;
    Glib::ustring BAmethod;
    Glib::ustring TMmethod;
    Glib::ustring Dirmethod;
    Glib::ustring HSmethod;
    double sigma;
    double offset;
    double lowthr;
    int rescon;
    int resconH;
    int reschro;
    int resblur;
    int resblurc;
    double tmrs;
    double edgs;
    double scale;
    double gamma;
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
    int edgesensi;
    int edgeampli;
    int contrast;
    int edgrad;
    double edgeffect;
    int edgval;
    int edgthresh;
    int thr;
    int thrH;
    int radius;
    double skinprotect;
    double chrwav;
    double bluwav;
    Threshold<int> hueskin;
    Threshold<int> hueskin2;
    Threshold<int> hllev;
    Threshold<int> bllev;
    Threshold<int> pastlev;
    Threshold<int> satlev;
    Threshold<int> edgcont;
    Threshold<double> level0noise;
    Threshold<double> level1noise;
    Threshold<double> level2noise;
    Threshold<double> level3noise;

    WaveletParams();

    bool operator ==(const WaveletParams& other) const;
    bool operator !=(const WaveletParams& other) const;

    void getCurves(
        WavCurve& cCurve,
        Wavblcurve& tCurve,
        WavOpacityCurveRG& opacityCurveLUTRG,
        WavOpacityCurveSH& opacityCurveLUTSH,
        WavOpacityCurveBY& opacityCurveLUTBY,
        WavOpacityCurveW& opacityCurveLUTW,
        WavOpacityCurveWL& opacityCurveLUTWL
    ) const;
};

/**
* Directional pyramid equalizer params
*/
struct DirPyrEqualizerParams {
    bool enabled;
    bool gamutlab;
    double mult[6];
    double threshold;
    double skinprotect;
    Threshold<int> hueskin;
    Glib::ustring cbdlMethod;

    DirPyrEqualizerParams();

    bool operator ==(const DirPyrEqualizerParams& other) const;
    bool operator !=(const DirPyrEqualizerParams& other) const;
};

/**
 * HSV equalizer params
 */
struct HSVEqualizerParams {
    bool enabled;
    std::vector<double> hcurve;
    std::vector<double> scurve;
    std::vector<double> vcurve;

    HSVEqualizerParams();

    bool operator ==(const HSVEqualizerParams& other) const;
    bool operator !=(const HSVEqualizerParams& other) const;
};

/**
 *  Film simualtion params
 */
struct FilmSimulationParams {
    bool enabled;
    Glib::ustring clutFilename;
    int strength;

    FilmSimulationParams();

    bool operator ==(const FilmSimulationParams& other) const;
    bool operator !=(const FilmSimulationParams& other) const;
};

struct SoftLightParams {
    bool enabled;
    int strength;

    SoftLightParams();

    bool operator==(const SoftLightParams &other) const;
    bool operator!=(const SoftLightParams &other) const;
};


struct DehazeParams {
    bool enabled;
    int strength;
    bool showDepthMap;
    int depth;
    bool luminance;

    DehazeParams();

    bool operator==(const DehazeParams &other) const;
    bool operator!=(const DehazeParams &other) const;
};

/**
  * Parameters for RAW demosaicing, common to all sensor type
  */
struct RAWParams {
    /**
     * Parameters for RAW demosaicing specific to Bayer sensors
     */
    struct BayerSensor {
        enum class Method {
            AMAZE,
            AMAZEVNG4,
            RCD,
            RCDVNG4,
            DCB,
            DCBVNG4,
            LMMSE,
            IGV,
            AHD,
            EAHD,
            HPHD,
            VNG4,
            FAST,
            MONO,
            PIXELSHIFT,
            NONE
        };

        enum class PSMotionCorrectionMethod {
            OFF,
            AUTO,
            CUSTOM
        };

        enum class PSDemosaicMethod {
            AMAZE,
            AMAZEVNG4,
            RCDVNG4,
            LMMSE
        };

        Glib::ustring method;
        int border;
        int imageNum;
        int ccSteps;
        double black0;
        double black1;
        double black2;
        double black3;
        bool twogreen;
        int linenoise;
        enum class LineNoiseDirection {
            HORIZONTAL = 1,
            VERTICAL,
            BOTH,
            PDAF_LINES = 5
        };
        LineNoiseDirection linenoiseDirection;
        int greenthresh;
        int dcb_iterations;
        int lmmse_iterations;
        bool dualDemosaicAutoContrast;
        double dualDemosaicContrast;
        PSMotionCorrectionMethod pixelShiftMotionCorrectionMethod;
        double pixelShiftEperIso;
        double pixelShiftSigma;
        bool pixelShiftShowMotion;
        bool pixelShiftShowMotionMaskOnly;
        bool pixelShiftHoleFill;
        bool pixelShiftMedian;
        bool pixelShiftGreen;
        bool pixelShiftBlur;
        double pixelShiftSmoothFactor;
        bool pixelShiftEqualBright;
        bool pixelShiftEqualBrightChannel;
        bool pixelShiftNonGreenCross;
        Glib::ustring pixelShiftDemosaicMethod;
        bool dcb_enhance;
        bool pdafLinesFilter;

        BayerSensor();

        bool operator ==(const BayerSensor& other) const;
        bool operator !=(const BayerSensor& other) const;

        void setPixelShiftDefaults();

        static const std::vector<const char*>& getMethodStrings();
        static Glib::ustring getMethodString(Method method);

        static const std::vector<const char*>& getPSDemosaicMethodStrings();
        static Glib::ustring getPSDemosaicMethodString(PSDemosaicMethod method);
    };

    /**
     * Parameters for RAW demosaicing specific to X-Trans sensors
     */
    struct XTransSensor {
        enum class Method {
            FOUR_PASS,
            THREE_PASS,
            TWO_PASS,
            ONE_PASS,
            FAST,
            MONO,
            NONE
        };

        Glib::ustring method;
        bool dualDemosaicAutoContrast;
        double dualDemosaicContrast;
        int border;
        int ccSteps;
        double blackred;
        double blackgreen;
        double blackblue;

        XTransSensor();

        bool operator ==(const XTransSensor& other) const;
        bool operator !=(const XTransSensor& other) const;

        static const std::vector<const char*>& getMethodStrings();
        static Glib::ustring getMethodString(Method method);
    };

    BayerSensor bayersensor;         ///< RAW parameters for Bayer sensors
    XTransSensor xtranssensor;       ///< RAW parameters for X-Trans sensors

    enum class FlatFieldBlurType {
        AREA,
        V,
        H,
        VH,
    };

    Glib::ustring dark_frame;
    bool df_autoselect;

    Glib::ustring ff_file;
    bool ff_AutoSelect;
    int ff_BlurRadius;
    Glib::ustring ff_BlurType;
    bool ff_AutoClipControl;
    int ff_clipControl;

    bool ca_autocorrect;
    bool ca_avoidcolourshift;
    int caautoiterations;
    double cared;
    double cablue;

    // exposure before interpolation
    double expos;

    struct PreprocessWB {
        enum class Mode {
            CAMERA = 0,
            AUTO
        };

        Mode mode;

        PreprocessWB();

        bool operator ==(const PreprocessWB& other) const;
        bool operator !=(const PreprocessWB& other) const;
    };

    PreprocessWB preprocessWB;

    bool hotPixelFilter;
    bool deadPixelFilter;
    int hotdeadpix_thresh;

    RAWParams();

    bool operator ==(const RAWParams& other) const;
    bool operator !=(const RAWParams& other) const;

    static const std::vector<const char*>& getFlatFieldBlurTypeStrings();
    static Glib::ustring getFlatFieldBlurTypeString(FlatFieldBlurType type);
};

/**
  * Parameters of film negative
  */
struct FilmNegativeParams {
    bool enabled;
    double redRatio;
    double greenExp;
    double blueRatio;

    double redBase;
    double greenBase;
    double blueBase;
    
    FilmNegativeParams();

    bool operator ==(const FilmNegativeParams& other) const;
    bool operator !=(const FilmNegativeParams& other) const;
};

/**
  * This class holds all the processing parameters applied on the images
  */
class ProcParams
{

public:
    ToneCurveParams         toneCurve;       ///< Tone curve parameters
    LCurveParams            labCurve;        ///< CIELAB luminance curve parameters
    RetinexParams           retinex;         ///< Retinex parameters
    LocalContrastParams     localContrast;   ////< Local contrast parameters
    RGBCurvesParams         rgbCurves;       ///< RGB curves parameters
    ColorToningParams       colorToning;     ///< Color Toning parameters
    SharpeningParams        sharpening;      ///< Sharpening parameters
    SharpeningParams        prsharpening;    ///< Sharpening parameters for post resize sharpening
    CaptureSharpeningParams pdsharpening;    ///< Sharpening parameters for post demosaic sharpening
    SharpenEdgeParams       sharpenEdge;     ///< Sharpen edge parameters
    SharpenMicroParams      sharpenMicro;    ///< Sharpen microcontrast parameters
    VibranceParams          vibrance;        ///< Vibrance parameters
    WBParams                wb;              ///< White balance parameters
    ColorAppearanceParams   colorappearance;
    DefringeParams          defringe;        ///< Defringing parameters
    ImpulseDenoiseParams    impulseDenoise;  ///< Impulse denoising parameters
    DirPyrDenoiseParams     dirpyrDenoise;   ///< Directional Pyramid denoising parameters
    EPDParams               epd;             ///< Edge Preserving Decomposition parameters
    FattalToneMappingParams fattal;          ///< Fattal02 tone mapping
    SHParams                sh;              ///< Shadow/highlight enhancement parameters
    CropParams              crop;            ///< Crop parameters
    CoarseTransformParams   coarse;          ///< Coarse transformation (90, 180, 270 deg rotation, h/v flipping) parameters
    CommonTransformParams   commonTrans;     ///< Common transformation parameters (autofill)
    RotateParams            rotate;          ///< Rotation parameters
    DistortionParams        distortion;      ///< Lens distortion correction parameters
    LensProfParams          lensProf;        ///< Lens correction profile parameters
    PerspectiveParams       perspective;     ///< Perspective correction parameters
    GradientParams          gradient;        ///< Gradient filter parameters
    LocallabParams          locallab;        ///< Local lab parameters
    PCVignetteParams        pcvignette;      ///< Post-crop vignette filter parameters
    CACorrParams            cacorrection;    ///< Lens c/a correction parameters
    VignettingParams        vignetting;      ///< Lens vignetting correction parameters
    ChannelMixerParams      chmixer;         ///< Channel mixer parameters
    BlackWhiteParams        blackwhite;      ///< Black&  White parameters
    ResizeParams            resize;          ///< Resize parameters
    ColorManagementParams   icm;             ///< profiles/color spaces used during the image processing
    RAWParams               raw;             ///< RAW parameters before demosaicing
    WaveletParams           wavelet;         ///< Wavelet parameters
    DirPyrEqualizerParams   dirpyrequalizer; ///< directional pyramid wavelet parameters
    HSVEqualizerParams      hsvequalizer;    ///< hsv wavelet parameters
    FilmSimulationParams    filmSimulation;  ///< film simulation parameters
    SoftLightParams         softlight;       ///< softlight parameters
    DehazeParams            dehaze;          ///< dehaze parameters
    FilmNegativeParams      filmNegative;    ///< Film negative parameters
    int                     rank;            ///< Custom image quality ranking
    int                     colorlabel;      ///< Custom color label
    bool                    inTrash;         ///< Marks deleted image
    Glib::ustring           appVersion;      ///< Version of the application that generated the parameters
    int                     ppVersion;       ///< Version of the PP file from which the parameters have been read

    MetaDataParams          metadata;        ///< Metadata parameters
    ExifPairs               exif;            ///< List of modifications appplied on the exif tags of the input image
    IPTCPairs               iptc;            ///< The IPTC tags and values to be saved to the output image

    /**
      * The constructor only sets the hand-wired defaults.
      */
    ProcParams();
    /**
      * Sets the hand-wired defaults parameters.
      */
    void setDefaults();
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
    int save(const Glib::ustring& fname, const Glib::ustring& fname2 = Glib::ustring(), bool fnameAbsolute = true, ParamsEdited* pedited = nullptr);
    /**
      * Loads the parameters from a file.
      * @param fname the name of the file
      * @params pedited pointer to a ParamsEdited object (optional) to store which values has been loaded
      * @return Error code (=0 if no error)
      */
    int load(const Glib::ustring& fname, ParamsEdited* pedited = nullptr);

    /** Creates a new instance of ProcParams.
      * @return a pointer to the new ProcParams instance. */
    static ProcParams* create();

    /** Destroys an instance of ProcParams.
      * @param pp a pointer to the ProcParams instance to destroy. */
    static void destroy(ProcParams* pp);

    static void init();
    static void cleanup();

    bool operator ==(const ProcParams& other) const;
    bool operator !=(const ProcParams& other) const;

private:
    /** Write the ProcParams's text in the file of the given name.
    * @param fname the name of the file
    * @param content the text to write
    * @return Error code (=0 if no error)
    * */
    int write(const Glib::ustring& fname, const Glib::ustring& content) const;
};

/**
  * This class associate a ProcParams object and a ParamEdited object through a pointer
  * to instance of each type in order to handle partial pp3 file loading (and later maybe
  * saving too)
  *
  * PartialProfile is not responsible of ProcParams and ParamsEdited object creation
  * and hence is not responsible of their destructions. The function that instantiate
  * PartialProfile object has to handle all this itself.
  */
class PartialProfile :
    public NonCopyable
{
public:
    PartialProfile(bool createInstance = false, bool paramsEditedValue = false);
    explicit PartialProfile(ProcParams* pp, ParamsEdited* pe = nullptr, bool fullCopy = false);
    explicit PartialProfile(const ProcParams* pp, const ParamsEdited* pe = nullptr);
    void deleteInstance();
    void clearGeneral();
    int  load(const Glib::ustring& fName);
    void set(bool v);
    void applyTo(ProcParams* destParams, bool fromLastSaved = false) const ;

    rtengine::procparams::ProcParams* pparams;
    ParamsEdited* pedited;
};

/**
  * This class automatically create the pparams and pedited instance in the constructor,
  * and automatically delete them in the destructor. This class has been mostly created
  * to be used with vectors, which use the default constructor/destructor
  */
class AutoPartialProfile :
    public PartialProfile
{
public:
    AutoPartialProfile();
    ~AutoPartialProfile();
};

}
}
