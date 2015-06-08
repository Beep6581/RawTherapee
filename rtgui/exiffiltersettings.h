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
#ifndef _EXIFFILTERSETTINGS_
#define _EXIFFILTERSETTINGS_

#include <set>
#include <string>

class ExifFilterSettings {

    public:
        std::set<std::string> filetypes;
        std::set<std::string> cameras;
        std::set<std::string> lenses;
        std::set<std::string> expcomp;
        double fnumberFrom;
        double fnumberTo;
        double shutterFrom;
        double shutterTo;
        double focalFrom;
        double focalTo;
        int isoFrom;
        int isoTo;

		bool filterFNumber;
		bool filterShutter;
		bool filterFocalLen;
		bool filterISO;
		bool filterExpComp;
		bool filterCamera;
		bool filterLens;
		bool filterFiletype;
		
         ExifFilterSettings ();
    void clear ();
};

#endif

