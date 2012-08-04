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
#include "rtXmp.h"

class ParamsEdited;

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
        int     brightness;
        int     contrast;
        int     saturation;
        bool    avoidclip;
        bool    enable_saturationlimiter;
        double  saturationlimit;
};

/**
  * Parameters of the RGB curves
  */
class RGBCurvesParams {

    public:
        std::vector<double>   rcurve;
        std::vector<double>   gcurve;
        std::vector<double>   bcurve;
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
        bool    enabled;
        int     pastels;
        int     saturated;
        int     psthreshold;
        bool    protectskins;
        bool    avoidcolorshift;
        bool    pastsattog;
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

    WBEntry(Glib::ustring p, enum WBTypes t, Glib::ustring l, int temp) : ppLabel(p), type(t), GUILabel(l), temperature(temp) {};
};

class WBParams {

    public:
	    static std::vector<WBEntry*> wbEntries;
        Glib::ustring   method;
        int             temperature;
        double          green;

        static void     init();
        static void     cleanup();
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

//EPD related parameters.
class EPDParams{
public:
	bool enabled;
	double Strength;
	double EdgeStopping;
	double Scale;
	int ReweightingIterates;
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
        bool          blendCMSMatrix;
        short         preferredProfile;
        Glib::ustring working;
        Glib::ustring output;
        static const Glib::ustring NoICMString;      
        
        Glib::ustring gamma;
		double gampos;
		double slpos;
		bool freegamma;
		
};

/**
  * Typedef for representing a key/value for the exif metadata information
typedef std::map<Glib::ustring, Glib::ustring> ExifPairs;
  */

/**
  * The IPTC key/value pairs
typedef std::map<Glib::ustring, std::vector<Glib::ustring> > IPTCPairs;
  */

/**
* Directional pyramid equalizer params
*/
class DirPyrEqualizerParams {
	
	public:
		bool enabled;
		double mult[8];
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
  * Parameters for RAW demosaicing
  */
class RAWParams {

    public:
		enum eMethod{eahd,hphd,vng4,dcb,amaze,ahd,fast,
					numMethods }; // This MUST be the last enum
		static const char *methodstring[numMethods];

		enum eFlatFileBlurType{/*parametric,*/area_ff,v_ff,h_ff,vh_ff,
								numFlatFileBlurTypes }; // This MUST be the last enum
		static const char *ff_BlurTypestring[numFlatFileBlurTypes];
	

	    Glib::ustring dark_frame;
	    bool df_autoselect;
	
		Glib::ustring ff_file;
		bool ff_AutoSelect;
		int ff_BlurRadius;
		Glib::ustring ff_BlurType;
	
		bool ca_autocorrect;
		double cared;
		double cablue;

		// exposure before interpolation
		double expos;
		double preser; 
		double blackzero;
		double blackone;
		double blacktwo;
		double blackthree;
		bool twogreen;
		bool hotdeadpix_filt;
		int hotdeadpix_thresh;
		int	linenoise;
		int greenthresh;
        int ccSteps;
        Glib::ustring dmethod;
        int dcb_iterations;
        bool dcb_enhance;
        bool all_enhance;
		
};

/**
  * This class holds all the processing parameters applied on the images
  */
class ProcParams {

    public:
        ToneCurveParams         toneCurve;       ///< Tone curve parameters
        LCurveParams            labCurve;        ///< CIELAB luminance curve parameters
        RGBCurvesParams         rgbCurves;       ///< RGB curves parameters
        SharpeningParams        sharpening;      ///< Sharpening parameters
        SharpenEdgeParams       sharpenEdge;     ///< Sharpen edge parameters
        SharpenMicroParams      sharpenMicro;    ///< Sharpen microcontrast parameters
        VibranceParams          vibrance;        ///< Vibrance parameters
        //ColorBoostParams        colorBoost;      ///< Color boost parameters
        WBParams                wb;              ///< White balance parameters
        //ColorShiftParams        colorShift;      ///< Color shift parameters
        //LumaDenoiseParams       lumaDenoise;     ///< Luminance denoising parameters
        //ColorDenoiseParams      colorDenoise;    ///< Color denoising parameters
        DefringeParams          defringe;        ///< Defringing parameters
        ImpulseDenoiseParams    impulseDenoise;  ///< Impulse denoising parameters
        DirPyrDenoiseParams     dirpyrDenoise;   ///< Directional Pyramid denoising parameters
        EPDParams               edgePreservingDecompositionUI;
        SHParams                sh;              ///< Shadow/highlight enhancement parameters
        CropParams              crop;            ///< Crop parameters
        CoarseTransformParams   coarse;          ///< Coarse transformation (90, 180, 270 deg rotation, h/v flipping) parameters
        CommonTransformParams	commonTrans;     ///< Common transformation parameters (autofill)
        RotateParams            rotate;          ///< Rotation parameters
        DistortionParams        distortion;      ///< Lens distortion correction parameters
        LensProfParams          lensProf;        ///< Lens correction profile parameters
        PerspectiveParams       perspective;     ///< Perspective correction parameters
        CACorrParams            cacorrection;    ///< Lens c/a correction parameters
        VignettingParams        vignetting;      ///< Lens vignetting correction parameters
        ChannelMixerParams      chmixer;         ///< Channel mixer parameters
        HRecParams              hlrecovery;      ///< Highlight recovery parameters
        ResizeParams            resize;          ///< Resize parameters
        ColorManagementParams   icm;             ///< profiles/color spaces used during the image processing
		
        RAWParams               raw;             ///< RAW parameters before demosaicing
        DirPyrEqualizerParams   dirpyrequalizer; ///< directional pyramid equalizer parameters
        HSVEqualizerParams      hsvequalizer;    ///< hsv equalizer parameters

        Glib::ustring           appVersion;      ///< Version of the application that generated the parameters
        int                     ppVersion;       ///< Version of the PP file from which the parameters have been read


      /**
        * The constructor only sets the hand-wired defaults.
        */
        ProcParams ();
      /**
        * Sets the hand-wired defaults parameters.
        */
        void setDefaults ();

      /**
        * Loads the parameters from a file.
        * @param fname the name of the file
        * @params pedited pointer to a ParamsEdited object (optional) to store which values has been loaded
        * @return Error code (=0 if no error)
        */
        int load        (Glib::ustring fname, ParamsEdited* pedited=NULL, int *rank=NULL);

        int saveIntoXMP (Exiv2::XmpData &xmpData, const std::string& baseKey) const;
        int loadFromXMP (Exiv2::XmpData &xmpData, const std::string& baseKey);
        int saveParams  (Glib::ustring fname) const;
        int loadParams  (Glib::ustring fname);

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

        PartialProfile      (bool createInstance=false);
        PartialProfile      (ProcParams* pp, ParamsEdited* pe=NULL, bool fullCopy=false);
        PartialProfile      (const ProcParams* pp, const ParamsEdited* pe=NULL);
        void deleteInstance ();
        void clearGeneral   ();
        int  load           (Glib::ustring fName);
        void set            (bool v);
        void applyTo        (ProcParams *destParams) const ;
};

}
}
#endif
