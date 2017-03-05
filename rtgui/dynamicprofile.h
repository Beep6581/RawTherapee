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
#include <vector>
#include "options.h"


class DynamicProfileRule {
public:
    template <class T>
    struct Range {
        T min;
        T max;
        explicit Range(T l=T(), T u=T()): min(l), max(u) {}

        bool operator()(T val) const
        {
            return val >= min && val <= max;
        }
    };

    struct Optional {
        Glib::ustring value;
        bool enabled;
        explicit Optional(const Glib::ustring v="", bool e=false):
            value(v), enabled(e) {}

        bool operator()(const Glib::ustring &val) const;
    };
    
    DynamicProfileRule();
    bool matches(const rtengine::ImageMetaData *im) const;
    bool operator<(const DynamicProfileRule &other) const;

    int serial_number;
    Range<int> iso;
    Range<double> fnumber;
    Range<double> focallen;
    Range<double> shutterspeed;
    Range<double> expcomp;
    Optional camera;
    Optional lens;
    Glib::ustring profilepath;
};


bool loadDynamicProfileRules(std::vector<DynamicProfileRule> &out);
bool storeDynamicProfileRules(
    const std::vector<DynamicProfileRule> &rules);

rtengine::procparams::PartialProfile *loadDynamicProfile(
    const rtengine::ImageMetaData *im);


#endif // _DYNAMICPROFILE_H_
