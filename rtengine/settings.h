/* -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  RawTherapee is distributed in the hope that it will be useful,
 * 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <glibmm/ustring.h>
namespace rtengine
{

/** This structure holds the global parameters used by the RT engine. */
class Settings
{
public:
    Glib::ustring   iccDirectory;           ///< The directory containing the possible output icc profiles
    int             viewingdevice;          // white of output device (D50...D65..)
    int             viewingdevicegrey;      // level of grey output device
    int             viewinggreySc;          // level of grey Scene
    int             leveldnv;               // level of crop denoise
    int             leveldnti;              // size of tiles denoise
    int             leveldnaut;             // level of auto denoise
    int             leveldnliss;            // level of auto multi zone
    int             leveldnautsimpl;        // STD or EXPERT

    Glib::ustring   printerProfile;         ///< ICC profile name used for soft-proofing a printer output
    int             printerIntent;          ///< Colorimetric intent used with the above profile
    bool            printerBPC;             ///< Black Point Compensation for the Labimage->Printer->Monitor transform
    Glib::ustring   monitorProfile;         ///< ICC profile name used for the monitor
    int             monitorIntent;          ///< Colorimetric intent used with the above profile
    bool            monitorBPC;             ///< Black Point Compensation for the Labimage->Monitor transform (directly, i.e. not soft-proofing and no WCS in between)
    bool            autoMonitorProfile;     ///< Try to auto-determine the correct monitor color profile
    bool            autocielab;
    bool            rgbcurveslumamode_gamut;// controls gamut enforcement for RGB curves in lumamode
    bool            verbose;
    Glib::ustring   darkFramesPath;         ///< The default directory for dark frames
    Glib::ustring   flatFieldsPath;         ///< The default directory for flat fields
    Glib::ustring   cameraProfilesPath;     ///< The default directory for camera profiles
    Glib::ustring   lensProfilesPath;       ///< The default directory for lens profiles
    bool            enableLibRaw;           ///< Use LibRaw to decode raw images.

    Glib::ustring   adobe;                  // filename of AdobeRGB1998 profile (default to the bundled one)
    Glib::ustring   prophoto;               // filename of Prophoto     profile (default to the bundled one)
    Glib::ustring   widegamut;              // filename of WidegamutRGB profile (default to the bundled one)
    Glib::ustring   beta;                   // filename of BetaRGB      profile (default to the bundled one)
    Glib::ustring   best;                   // filename of BestRGB      profile (default to the bundled one)
    Glib::ustring   bruce;                  // filename of BruceRGB     profile (default to the bundled one)
    Glib::ustring   srgb;                   // filename of sRGB         profile (default to the bundled one)
    Glib::ustring   rec2020;                // filename of Rec2020      profile (default to the bundled one)
    Glib::ustring   ACESp0;                 // filename of ACES P0      profile (default to the bundled one)
    Glib::ustring   JDCmax;                 // filename of JDCmax      profile (default to the bundled one)
    Glib::ustring   ACESp1;                 // filename of ACES P1      profile (default to the bundled one)
    Glib::ustring   DCIP3;                 // filename of DCIP3         profile (default to the bundled one)

    bool            gamutICC; // no longer used
    bool            gamutLch;
    bool            HistogramWorking;       // true: histogram is display the value of the image computed in the Working profile
    // false: histogram is display the value of the image computed in the Output profile
    int             amchroma;
    int             amchromajz;
    int             protectred;
    double          protectredh;
    double          nrauto;
    double          nrautomax;
    double          nrhigh;
    int             nrwavlevel;
    bool            daubech;
    bool            ciebadpixgauss;
    int             CRI_color; // Number for display Lab value; 0 = disabled
    int             denoiselabgamma; // 0=gamma 26 11   1=gamma 40 5  2 =gamma 55 10
    //  double          colortoningab; //
    //  double          decaction;
    //  bool            bw_complementary;
    double          level0_cbdl;
    double          level123_cbdl;
    Glib::ustring   lensfunDbDirectory; // The directory containing the lensfun database. If empty, the system defaults will be used, as described in https://lensfun.github.io/manual/latest/dbsearch.html
    Glib::ustring   lensfunDbBundleDirectory;
    int             cropsleep;
    double          reduchigh;
    double          reduclow;
    bool            detectshape;
    bool            fftwsigma;
    int             previewselection;
    double          cbdlsensi;
//    bool            showtooltip;
    bool            itcwb_enable;
    double          itcwb_deltaspec;
    double          itcwb_powponder;
    double basecorlog;
//wavelet levels
    double          edghi;
    double          edglo;
    double          limrad;


    enum class ThumbnailInspectorMode {
        JPEG,
        RAW,
        RAW_IF_NOT_JPEG_FULLSIZE
    };
    ThumbnailInspectorMode thumbnail_inspector_mode;

    enum class XmpSidecarStyle {
        STD, // FILENAME.xmp for FILENAME.ext
        EXT  // FILENAME.ext.xmp for FILENAME.ext
    };
    XmpSidecarStyle xmp_sidecar_style;

    enum class MetadataXmpSync {
        NONE,
        READ,
        READ_WRITE
    };
    MetadataXmpSync metadata_xmp_sync;

    /** Creates a new instance of Settings.
      * @return a pointer to the new Settings instance. */
    static Settings* create();
    /** Destroys an instance of Settings.
      * @param s a pointer to the Settings instance to destroy. */
    static void      destroy(Settings* s);
};
extern const Settings* settings;
}
