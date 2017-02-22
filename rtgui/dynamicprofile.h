/* -*- C++ -*-
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2017 Alberto Griggio
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
#ifndef _DYNAMICPROFILE_H_
#define _DYNAMICPROFILE_H_

#include <glibmm.h>

#include "options.h"


class DynamicProfileEntry {
public:
    DynamicProfileEntry();
    bool matches(const rtengine::ImageMetaData *im);
    bool operator<(const DynamicProfileEntry &other) const;

    int serial_number;

    bool has_iso;
    int iso_min;
    int iso_max;

    bool has_fnumber;
    double fnumber_min;
    double fnumber_max;

    bool has_focallen;
    double focallen_min;
    double focallen_max;

    bool has_shutterspeed;
    double shutterspeed_min;
    double shutterspeed_max;

    bool has_expcomp;
    double expcomp_min;
    double expcomp_max;

    bool has_make;
    std::string make;

    bool has_model;
    std::string model;

    bool has_lens;
    std::string lens;

    Glib::ustring profilepath;
};


bool loadDynamicProfileEntries(std::vector<DynamicProfileEntry> &out);
rtengine::procparams::PartialProfile *loadDynamicProfile(
    const rtengine::ImageMetaData *im);


#endif // _DYNAMICPROFILE_H_
