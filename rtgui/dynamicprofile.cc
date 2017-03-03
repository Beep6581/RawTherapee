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
#include "profilestore.h"
#include <stdlib.h>
#include <glibmm/regex.h>

using namespace rtengine;
using namespace rtengine::procparams;

namespace {

const int ISO_MAX = 512000;
const double FNUMBER_MAX = 100.0;
const double FOCALLEN_MAX = 10000.0;
const double SHUTTERSPEED_MIN = 1.0/10000.0;
const double SHUTTERSPEED_MAX = 1000.0;
const double EXPCOMP_MIN = -20.0;
const double EXPCOMP_MAX = 20.0;

} // namespace


bool DynamicProfileEntry::Optional::operator()(const Glib::ustring &val) const
{
    if (!enabled) {
        return true;
    }
    if (value.find("re:") == 0) {
        // this is a regexp
        return Glib::Regex::match_simple(value.substr(3), val,
                                         Glib::REGEX_CASELESS);
    } else {
        // normal string comparison
        return value.casefold() == val.casefold();
    }
}


DynamicProfileEntry::DynamicProfileEntry():
    serial_number(0),
    iso(0, ISO_MAX),
    fnumber(0, FNUMBER_MAX),
    focallen(0, FOCALLEN_MAX),
    shutterspeed(SHUTTERSPEED_MIN, SHUTTERSPEED_MAX),
    expcomp(EXPCOMP_MIN, EXPCOMP_MAX)
{
}


bool DynamicProfileEntry::operator<(const DynamicProfileEntry &other) const
{
    return serial_number < other.serial_number;
}


bool DynamicProfileEntry::matches(const rtengine::ImageMetaData *im)
{
    return (iso(im->getISOSpeed())
            && fnumber(im->getFNumber())
            && focallen(im->getFocalLen())
            && shutterspeed(im->getShutterSpeed())
            && expcomp(im->getExpComp())
            && camera(im->getCamera())
            && lens(im->getLens()));
}

namespace {

void get_int_range(DynamicProfileEntry::Range<int> &dest,
                   const Glib::KeyFile &kf, const Glib::ustring &group,
                   const Glib::ustring &key)
{
    try {
        int min = kf.get_integer(group, key + "_min");
        int max = kf.get_integer(group, key + "_max");
        if (min <= max) {
            dest.min = min;
            dest.max = max;
        }
    } catch (Glib::KeyFileError &e) {
    }
}


void get_double_range(DynamicProfileEntry::Range<double> &dest,
                      const Glib::KeyFile &kf, const Glib::ustring &group,
                      const Glib::ustring &key)
{
    try {
        int min = kf.get_double(group, key + "_min");
        int max = kf.get_double(group, key + "_max");
        if (min <= max) {
            dest.min = min;
            dest.max = max;
        }
    } catch (Glib::KeyFileError &e) {
    }
}


void get_optional(DynamicProfileEntry::Optional &dest,
                  const Glib::KeyFile &kf, const Glib::ustring &group,
                  const Glib::ustring &key)
{
    try {
        bool e = kf.get_boolean(group, key + "_enabled");
        if (e) {
            Glib::ustring s = kf.get_string(group, key + "_value");
            dest.enabled = e;
            dest.value = s;
        }
    } catch (Glib::KeyFileError &) {
    }
}

void set_int_range(Glib::KeyFile &kf, const Glib::ustring &group,
                   const Glib::ustring &key,
                   const DynamicProfileEntry::Range<int> &val)
{
    kf.set_integer(group, key + "_min", val.min);
    kf.set_integer(group, key + "_max", val.max);
}

void set_double_range(Glib::KeyFile &kf, const Glib::ustring &group,
                      const Glib::ustring &key,
                      const DynamicProfileEntry::Range<double> &val)
{
    kf.set_double(group, key + "_min", val.min);
    kf.set_double(group, key + "_max", val.max);
}

void set_optional(Glib::KeyFile &kf, const Glib::ustring &group,
                  const Glib::ustring &key,
                  const DynamicProfileEntry::Optional &val)
{
    kf.set_boolean(group, key + "_enabled", val.enabled);
    kf.set_string(group, key + "_value", val.value);
}

} // namespace


bool loadDynamicProfileEntries(std::vector<DynamicProfileEntry> &out)
{
    out.clear();
    Glib::KeyFile kf;
    try {
        if (!kf.load_from_file(
                Glib::build_filename(Options::rtdir, "dynamicprofile.cfg"))) {
            return false;
        }
    } catch (Glib::Error &e) {
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
        get_int_range(entry.iso, kf, group, "iso");
        get_double_range(entry.fnumber, kf, group, "fnumber");
        get_double_range(entry.focallen, kf, group, "focallen");
        get_double_range(entry.shutterspeed, kf, group, "shutterspeed");
        get_double_range(entry.expcomp, kf, group, "expcomp");
        get_optional(entry.camera, kf, group, "camera");
        get_optional(entry.lens, kf, group, "lens");
        try {
            entry.profilepath = kf.get_string(group, "profilepath");
        } catch (Glib::KeyFileError &) {
            out.pop_back();
        }
    }
    std::sort(out.begin(), out.end());
    return true;
}


bool storeDynamicProfileEntries(const std::vector<DynamicProfileEntry> &entries)
{
    printf("saving dynamic profiles...\n");
    Glib::KeyFile kf;
    for (auto &entry : entries) {
        std::ostringstream buf;
        buf << "entry " << entry.serial_number;
        Glib::ustring group = buf.str();
        set_int_range(kf, group, "iso", entry.iso);
        set_double_range(kf, group, "fnumber", entry.fnumber);
        set_double_range(kf, group, "focallen", entry.focallen);
        set_double_range(kf, group, "shutterspeed", entry.shutterspeed);
        set_double_range(kf, group, "expcomp", entry.expcomp);
        set_optional(kf, group, "camera", entry.camera);
        set_optional(kf, group, "lens", entry.lens);
        kf.set_string(group, "profilepath", entry.profilepath);
    }
    return kf.save_to_file(
        Glib::build_filename(Options::rtdir, "dynamicprofile.cfg"));
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
                const PartialProfile *p =
                    profileStore.getProfile(entry.profilepath);
                if (p != nullptr) {
                    p->applyTo(ret->pparams);
                } else {
                    printf("ERROR loading matching profile\n");
                }
            }
        }
    }
    return ret;
}
