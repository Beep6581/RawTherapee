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
            bool            dualThreadEnabled;      ///< If true, the image processing operations with utilize two processor cores (if possible)
            Glib::ustring   iccDirectory;           ///< The directory containing the possible output icc profiles
            int             colorimetricIntent;     ///< Colorimetric intent used at color space conversions
            Glib::ustring   monitorProfile;         ///< ICC profile of the monitor (full path recommended)
			bool            autoMonitorProfile;     ///< Try to auto-determine the correct monitor color profile

            bool            verbose;
            Glib::ustring   darkFramesPath;         ///< The default directory for dark frames
            Glib::ustring   flatFieldsPath;         ///< The default directory for flat fields
            bool            LCMSSafeMode;           // If true, not OMP

        /** Creates a new instance of Settings.
          * @return a pointer to the new Settings instance. */
            static Settings* create  ();
        /** Destroys an instance of Settings.
          * @param s a pointer to the Settings instance to destroy. */
            static void      destroy (Settings* s);
    };
}

#endif

