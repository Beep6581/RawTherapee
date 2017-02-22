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

#include "dynamicprofile.h"
#include <stdlib.h>

using namespace rtengine;
using namespace rtengine::procparams;

DynamicProfileEntry::DynamicProfileEntry():
    serial_number(0),
    has_iso(false), iso_min(0), iso_max(1000000),
    has_fnumber(false), fnumber_min(0.0), fnumber_max(1000.0),
    has_focallen(false), focallen_min(0.0), focallen_max(1000000.0),
    has_shutterspeed(false), shutterspeed_min(1000.0), shutterspeed_max(1.0/1000000.0),
    has_expcomp(false), expcomp_min(-100.0), expcomp_max(100.0),
    has_make(false), make(""),
    has_model(false), model(""),
    has_lens(false), lens(""),
    profilepath("")
{
}


bool DynamicProfileEntry::operator<(const DynamicProfileEntry &other) const
{
    return serial_number < other.serial_number;
}


bool DynamicProfileEntry::matches(const rtengine::ImageMetaData *im)
{
    if (has_iso) {
        int iso = im->getISOSpeed();
        if (iso < iso_min || iso > iso_max) {
            return false;
        }
    }
    if (has_fnumber) {
        double fnumber = im->getFNumber();
        if (fnumber < fnumber_min || fnumber > fnumber_max) {
            return false;
        }
    }
    if (has_focallen) {
        double focallen = im->getFocalLen();
        if (focallen < focallen_min || focallen > focallen_max) {
            return false;
        }
    }
    if (has_shutterspeed) {
        double shutterspeed = im->getShutterSpeed();
        if (shutterspeed < shutterspeed_min || shutterspeed > shutterspeed_max){
            return false;
        }
    }
    if (has_expcomp) {
        double expcomp = im->getExpComp();
        if (expcomp < expcomp_min || expcomp > expcomp_max) {
            return false;
        }
    }
    if (has_make) {
        if (im->getMake() != make) {
            return false;
        }
    }
    if (has_model) {
        if (im->getModel() != model) {
            return false;
        }
    }
    if (has_lens) {
        if (im->getLens() != lens) {
            return false;
        }
    }
    return true;
}


bool loadDynamicProfileEntries(std::vector<DynamicProfileEntry> &out)
{
    out.clear();
    Glib::KeyFile kf;
    if (!kf.load_from_file(
            Glib::build_filename(Options::rtdir, "dynamicprofile.cfg"))) {
        return false;
    }
    printf("loading dynamic profiles...\n");
    auto groups = kf.get_groups();
    for (auto group : groups) {
        // groups are of the form "entry N", where N is a positive integer
        if (group.find("entry ") != 0) {
            return false;
        }
        std::istringstream buf(group.c_str() + 6);
        int serial = 0;
        if (!(buf >> serial) || !buf.eof()) {
            return false;
        }
        printf(" loading entry %d\n", serial);
        
        out.emplace_back(DynamicProfileEntry());
        DynamicProfileEntry &entry = out.back();
        entry.serial_number = serial;
        entry.has_iso = kf.get_boolean(group, "has_iso");
        entry.iso_min = kf.get_integer(group, "iso_min");
        entry.iso_max = kf.get_integer(group, "iso_max");

        entry.has_fnumber = kf.get_boolean(group, "has_fnumber");
        entry.fnumber_min = kf.get_double(group, "fnumber_min");
        entry.fnumber_max = kf.get_double(group, "fnumber_max");

        entry.has_focallen = kf.get_boolean(group, "has_focallen");
        entry.focallen_min = kf.get_double(group, "focallen_min");
        entry.focallen_max = kf.get_double(group, "focallen_max");

        entry.has_shutterspeed = kf.get_boolean(group, "has_shutterspeed");
        entry.shutterspeed_min = kf.get_double(group, "shutterspeed_min");
        entry.shutterspeed_max = kf.get_double(group, "shutterspeed_max");

        entry.has_expcomp = kf.get_boolean(group, "has_expcomp");
        entry.expcomp_min = kf.get_double(group, "expcomp_min");
        entry.expcomp_max = kf.get_double(group, "expcomp_max");

        entry.has_make = kf.get_boolean(group, "has_make");
        entry.make = kf.get_string(group, "make");

        entry.has_model = kf.get_boolean(group, "has_model");
        entry.model = kf.get_string(group, "model");

        entry.has_lens = kf.get_boolean(group, "has_lens");
        entry.lens = kf.get_string(group, "lens");

        entry.profilepath = kf.get_string(group, "profilepath");
    }
    std::sort(out.begin(), out.end());
    return true;
}


PartialProfile *loadDynamicProfile(const ImageMetaData *im)
{
    PartialProfile *ret = new PartialProfile(true, true);
    std::vector<DynamicProfileEntry> entries;
    if (loadDynamicProfileEntries(entries)) {
        for (auto &entry : entries) {
            if (entry.matches(im)) {
                printf("found matching profile %s\n",
                       entry.profilepath.c_str());
                PartialProfile p(true, true);
                if (!p.load(options.findProfilePath(entry.profilepath))) {
                    p.applyTo(ret->pparams);
                } else {
                    printf("ERROR loading matching profile\n");
                }
                p.deleteInstance();
            }
        }
    }
    return ret;
}

