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


class DynamicProfileEntry {
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

    template <class T>
    struct Optional {
        T value;
        bool enabled;
        explicit Optional(T v=T(), bool e=false): value(v), enabled(e) {}

        bool operator()(const T &val) const
        {
            return !enabled || value == val;
        }
    };
    
    DynamicProfileEntry();
    bool matches(const rtengine::ImageMetaData *im);
    bool operator<(const DynamicProfileEntry &other) const;

    int serial_number;
    Range<int> iso;
    Range<double> fnumber;
    Range<double> focallen;
    Range<double> shutterspeed;
    Range<double> expcomp;
    Optional<Glib::ustring> make;
    Optional<Glib::ustring> model;
    Optional<Glib::ustring> lens;
    Glib::ustring profilepath;
};


bool loadDynamicProfileEntries(std::vector<DynamicProfileEntry> &out);
bool storeDynamicProfileEntries(
    const std::vector<DynamicProfileEntry> &entries);

rtengine::procparams::PartialProfile *loadDynamicProfile(
    const rtengine::ImageMetaData *im);


#endif // _DYNAMICPROFILE_H_
