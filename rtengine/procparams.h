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

namespace rtengine {
namespace procparams {

/**
  * Parameters of the tone curve 
  */
class ToneCurveParams {

    public:
        bool        autoexp;
        double      clip;
        double      expcomp;
        std::vector<double>   curve;
        int         brightness;
        int         black;
        int         contrast;
		int         saturation;
        int         shcompr;
        int         hlcompr;
        int         hlcomprthresh;
};

/**
  * Parameters of the luminance curve
  */
class LCurveParams {

    public:
        std::vector<double>   lcurve;
		std::vector<double>   acurve;
		std::vector<double>   bcurve;
        int      	brightness;
        int         contrast;
		int         saturation;
};

/**
  * Parameters of the sharpening
  */
class SharpeningParams {

    public:
        bool    enabled;
        double  radius;
        int     amount;
        int     threshold;
        bool    edgesonly;
        double  edges_radius;
        int     edges_tolerance;
        bool    halocontrol;
        int     halocontrol_amount;
        
        Glib::ustring method;
        int     deconvamount;
        double  deconvradius;
        int     deconviter;
        int     deconvdamping;
};

/**
  * Parameters of the color boost
  */
class ColorBoostParams {

    public: 
        int     amount;
        bool    avoidclip;
        bool    enable_saturationlimiter;
        double  saturationlimit;
};

/**
  * Parameters of the white balance adjustments
  */
class WBParams {

    public:
        Glib::ustring   method;
        int             temperature;
        double          green;
};

/**
  * Parameters of the color shift
  */
class ColorShiftParams {

    public:
        double  a;
        double  b;
};

/**
  * Parameters of the luminance denoising
  */
class LumaDenoiseParams {

    public:
        bool    enabled;
        double  radius;
        int     edgetolerance;
};

/**
  * Parameters of the color denoising
  */
class ColorDenoiseParams {

    public:
        bool    enabled;
        double  radius;
        int     edgetolerance;
        bool    edgesensitive;
        int		amount;
};
	
	/**
	 * Parameters of defringing
	 */
	class DefringeParams {
		
    public:
        bool    enabled;
        double  radius;
        int     threshold;
	};
	
	
	/**
	 * Parameters of impulse denoising
	 */
	class ImpulseDenoiseParams {
		
    public:
        bool    enabled;
		int		thresh;

	};
	
	/**
	 * Parameters of the directional pyramid denoising
	 */
	class DirPyrDenoiseParams {
		
    public:
        bool    enabled;
        int		luma;
        int     chroma;
		float	gamma;
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
};

/**
  * Parameters of the coarse transformations like 90 deg rotations and h/v flipping
  */
class CoarseTransformParams {

    public:
        int     rotate;
        bool    hflip;
        bool    vflip;
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
		bool    uselensfun;
        double  amount;
};

/**
  * Parameters of the perspective correction
  */
class PerspectiveParams {

    public:
		int  horizontal;
		int  vertical;
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
class HRecParams {

    public:
        bool enabled;
		Glib::ustring method;
};

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
        bool          gammaOnInput;
        Glib::ustring working;
        Glib::ustring output;
};

/**
  * A class representing a key/value for the exif metadata information
  */
class ExifPair {

    public:
        Glib::ustring field;
        Glib::ustring value;
};

/**
  * The IPTC key/value pairs
  */
class IPTCPair {

    public:
        Glib::ustring field;
        std::vector<Glib::ustring> values;
};

/**
  * Wavelet equalizer params
  */
class EqualizerParams {

    public:
        bool enabled;
        int c[8];
};
	
/**
* Directional pyramid equalizer params
*/
class DirPyrEqualizerParams {
	
	public:
		bool enabled;
		double mult[8];
};
	
/**
 * Wavelet equalizer params
 */
class HSVEqualizerParams {
	
	public:
		bool enabled;
		Glib::ustring hsvchannel;
		int sat[8];
		int val[8];
		int hue[8];
};

/**
  * Parameters for RAW demosaicing
  */
class RAWParams {

    public:
		enum eMethod{eahd,hphd,vng4,dcb,amaze,ahd,fast,
					numMethods }; // This MUST be the last enum
		static const char *methodstring[numMethods];



	    Glib::ustring dark_frame;
	    bool df_autoselect;
		bool ca_autocorrect;
		double cared;
		double cablue;
		double expos;
		double preser; // expos
		bool hotdeadpix_filt;
		int	linenoise;
		int greenthresh;
        int ccSteps;
        Glib::ustring dmethod;
        int dcb_iterations;
        bool dcb_enhance;
};

/**
  * This class holds all the processing parameters applied on the images
  */
class ProcParams {

    public:
        ToneCurveParams         toneCurve;       ///< Tone curve parameters
        LCurveParams            labCurve;        ///< CIELAB luminance curve parameters
        SharpeningParams        sharpening;      ///< Sharpening parameters
        ColorBoostParams        colorBoost;      ///< Color boost parameters
        WBParams                wb;              ///< White balance parameters
        ColorShiftParams        colorShift;      ///< Color shift parameters
        LumaDenoiseParams       lumaDenoise;     ///< Luminance denoising parameters
        ColorDenoiseParams      colorDenoise;    ///< Color denoising parameters
        DefringeParams          defringe;        ///< Defringing parameters
        ImpulseDenoiseParams    impulseDenoise;  ///< Impulse denoising parameters
        DirPyrDenoiseParams     dirpyrDenoise;   ///< Directional Pyramid denoising parameters
        SHParams                sh;              ///< Shadow/highlight enhancement parameters
        CropParams              crop;            ///< Crop parameters
        CoarseTransformParams   coarse;          ///< Coarse transformation (90, 180, 270 deg rotation, h/v flipping) parameters
        CommonTransformParams	commonTrans;     ///< Common transformation parameters (autofill)
        RotateParams            rotate;          ///< Rotation parameters
        DistortionParams        distortion;      ///< Lens distortion correction parameters
        PerspectiveParams       perspective;     ///< Perspective correction parameters
        CACorrParams            cacorrection;    ///< Lens c/a correction parameters
        VignettingParams        vignetting;      ///< Lens vignetting correction parameters
        ChannelMixerParams      chmixer;         ///< Channel mixer parameters
        HRecParams              hlrecovery;      ///< Highlight recovery parameters
        ResizeParams            resize;          ///< Resize parameters
        ColorManagementParams   icm;             ///< profiles/color spaces used during the image processing
        EqualizerParams         equalizer;       ///< wavelet equalizer parameters
        RAWParams               raw;             ///< RAW parameters before demosaicing
        DirPyrEqualizerParams   dirpyrequalizer; ///< directional pyramid equalizer parameters
        HSVEqualizerParams      hsvequalizer;    ///< hsv equalizer parameters
        std::vector<ExifPair>   exif;            ///< List of modifications appplied on the exif tags of the input image
        std::vector<IPTCPair>   iptc;            ///< The IPTC tags and values to be saved to the output image
        int version;                             ///< Version of the file from which the parameters have been read

      /**
        * The constructor only sets the hand-wired defaults.
        */
        ProcParams          ();
      /**
        * Sets the hand-wired defaults parameters.
        */
        void    setDefaults ();
      /**
        * Saves the parameters to a file.
        * @param fname the name of the file
        * @return Error code (=0 if no error)
        */
        int     save        (Glib::ustring fname) const;
      /**
        * Loads the parameters from a file.
        * @param fname the name of the file
        * @return Error code (=0 if no error)
        */
        int     load        (Glib::ustring fname);

      /** Creates a new instance of ProcParams.
        * @return a pointer to the new ProcParams instance. */
        static ProcParams* create  ();

      /** Destroys an instance of ProcParams.
        * @param pp a pointer to the ProcParams instance to destroy. */
        static void        destroy (ProcParams* pp);

        bool operator== (const ProcParams& other);
        bool operator!= (const ProcParams& other);
};
}
}
#endif
