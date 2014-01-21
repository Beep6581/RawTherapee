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

namespace rtengine {

  /** This structure holds the global parameters used by the RT engine. */
    class Settings {
        public:
            Glib::ustring   iccDirectory;           ///< The directory containing the possible output icc profiles
            int             colorimetricIntent;     ///< Colorimetric intent used at color space conversions
			int				viewingdevice; 			// white of output device (D50...D65..)
			int				viewingdevicegrey; 			// level of grey output device
 
			Glib::ustring   monitorProfile;         ///< ICC profile of the monitor (full path recommended)
			bool            autoMonitorProfile;     ///< Try to auto-determine the correct monitor color profile
			bool			autocielab;
			bool			rgbcurveslumamode_gamut;// controls gamut enforcement for RGB curves in lumamode
            bool            verbose;
            Glib::ustring   darkFramesPath;         ///< The default directory for dark frames
            Glib::ustring   flatFieldsPath;         ///< The default directory for flat fields
			Glib::ustring   adobe;					// default name of AdobeRGB1998
			Glib::ustring   prophoto;				// default name of Prophoto
			Glib::ustring   prophoto10;				// default name of Prophoto
			
			Glib::ustring   widegamut;				//default name of WidegamutRGB
			Glib::ustring   beta;					// default name of BetaRGB
			Glib::ustring   best;					// default name of BestRGB
			Glib::ustring   bruce;					// default name of Bruce
			Glib::ustring   srgb;					// default name of SRGB space profile
			Glib::ustring   srgb10;					// default name of SRGB space profile
			
			bool            gamutICC; // no longer used
			bool            gamutLch;
			bool			ciecamfloat;
			int				amchroma;
			int             protectred;
			double          protectredh;
			bool			ciebadpixgauss;
			int             CRI_color; // N° for display Lab value  ; 0 disabled
		//	bool			bw_complementary;

        /** Creates a new instance of Settings.
          * @return a pointer to the new Settings instance. */
            static Settings* create  ();
        /** Destroys an instance of Settings.
          * @param s a pointer to the Settings instance to destroy. */
            static void      destroy (Settings* s);
    };
}

#endif

