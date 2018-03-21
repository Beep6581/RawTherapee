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
#ifndef _RTSETTINGS_
#define _RTSETTINGS_

#include "procparams.h"

namespace rtengine
{

/** This structure holds the global parameters used by the RT engine. */
class Settings
{
public:
    Glib::ustring   iccDirectory;           ///< The directory containing the possible output icc profiles
    int             viewingdevice;          // white of output device (D50...D65..)
    int             viewingdevicegrey;          // level of grey output device
    int             viewinggreySc;          // level of grey Scene
    int             leveldnv;           // level of crop denoise
    int             leveldnti;          // size of tiles denoise
    int             leveldnaut;             // level of auto denoise
    int             leveldnliss;            // level of auto multi zone
    int             leveldnautsimpl;            // STD or EXPERT

    Glib::ustring   printerProfile;         ///< ICC profile name used for soft-proofing a printer output
    RenderingIntent printerIntent;          ///< Colorimetric intent used with the above profile
    bool            printerBPC;             ///< Black Point Compensation for the Labimage->Printer->Monitor transform
    Glib::ustring   monitorProfile;         ///< ICC profile name used for the monitor
    RenderingIntent monitorIntent;          ///< Colorimetric intent used with the above profile
    bool            monitorBPC;             ///< Black Point Compensation for the Labimage->Monitor transform (directly, i.e. not soft-proofing and no WCS in between)
    bool            autoMonitorProfile;     ///< Try to auto-determine the correct monitor color profile
    bool            autocielab;
    bool            rgbcurveslumamode_gamut;// controls gamut enforcement for RGB curves in lumamode
    bool            verbose;
    Glib::ustring   darkFramesPath;         ///< The default directory for dark frames
    Glib::ustring   flatFieldsPath;         ///< The default directory for flat fields
    Glib::ustring   adobe;                  // default name of AdobeRGB1998
    Glib::ustring   prophoto;               // default name of Prophoto
    Glib::ustring   prophoto10;             // default name of Prophoto

    Glib::ustring   widegamut;              //default name of WidegamutRGB
    Glib::ustring   beta;                   // default name of BetaRGB
    Glib::ustring   best;                   // default name of BestRGB
    Glib::ustring   bruce;                  // default name of Bruce
    Glib::ustring   srgb;                   // default name of SRGB space profile
    Glib::ustring   srgb10;                 // default name of SRGB space profile
    Glib::ustring   rec2020;                   // default name of rec2020

    bool            gamutICC; // no longer used
    bool            gamutLch;
    bool            ciecamfloat;
    bool            HistogramWorking;
    int             amchroma;
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
    Glib::ustring   lensfunDbDirectory; ///< The directory containing the lensfun database. If empty, the system defaults will be used (as described in http://lensfun.sourceforge.net/manual/dbsearch.html)

    enum class ThumbnailInspectorMode {
        JPEG,
        RAW,
        RAW_IF_NOT_JPEG_FULLSIZE
    };
    ThumbnailInspectorMode thumbnail_inspector_mode;
    
    /** Creates a new instance of Settings.
      * @return a pointer to the new Settings instance. */
    static Settings* create  ();
    /** Destroys an instance of Settings.
      * @param s a pointer to the Settings instance to destroy. */
    static void      destroy (Settings* s);
};
}

#endif

