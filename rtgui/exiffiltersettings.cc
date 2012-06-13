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
#include "exiffiltersettings.h"

ExifFilterSettings::ExifFilterSettings () {

    clear (); 
}

void ExifFilterSettings::clear () {
    fnumberFrom = 100;
    fnumberTo = 0;
    shutterFrom = 100;
    shutterTo = 0;
    isoFrom = 100000000;
    isoTo = 0;
    focalFrom = 1e8;
    focalTo = 0;
    lenses.clear ();
    cameras.clear ();
    expcomp.clear ();
    filetypes.clear ();
	
	filterFNumber = false;
	filterShutter = false;
	filterFocalLen = false;
	filterISO = false;
	filterExpComp = false;
	filterCamera = false;
	filterLens = false;
	filterFiletype = false;
}
