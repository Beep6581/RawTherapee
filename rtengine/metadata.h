/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Alberto Griggio <alberto.griggio@gmail.com>
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

#pragma once

#include <glibmm.h>
#include <exiv2/exiv2.hpp>
#include <memory>
#include "procparams.h"

namespace rtengine {

class Exiv2Metadata final
{
public:
    Exiv2Metadata();
    explicit Exiv2Metadata(const Glib::ustring& path);
    Exiv2Metadata(const Glib::ustring& path, bool merge_xmp_sidecar);

    void load() const;
    
    Exiv2::ExifData& exifData();
    const Exiv2::ExifData& exifData() const;
    
    Exiv2::IptcData& iptcData();
    const Exiv2::IptcData& iptcData() const;
    
    Exiv2::XmpData& xmpData();
    const Exiv2::XmpData& xmpData() const;

    const Glib::ustring& filename() const;
    const rtengine::procparams::ExifPairs& exif() const;
    const rtengine::procparams::IPTCPairs& iptc() const;
    void setExif(const rtengine::procparams::ExifPairs& exif);
    void setIptc(const rtengine::procparams::IPTCPairs& iptc);
    
    void saveToImage(const Glib::ustring& path) const;
    void saveToXmp(const Glib::ustring& path) const;

    static Glib::ustring xmpSidecarPath(const Glib::ustring& path);
    static Exiv2::XmpData getXmpSidecar(const Glib::ustring& path);

    static void init();
    static void cleanup();
   
private:
    void do_merge_xmp(Exiv2::Image* dst) const;
    void import_exif_pairs(Exiv2::ExifData& out) const;
    void import_iptc_pairs(Exiv2::IptcData& out) const;
    void remove_unwanted(Exiv2::Image* dst) const;
    
    Glib::ustring src_;
    bool merge_xmp_;
    mutable std::shared_ptr<Exiv2::Image> image_;
    std::unique_ptr<rtengine::procparams::ExifPairs> exif_;
    std::unique_ptr<rtengine::procparams::IPTCPairs> iptc_;
    Exiv2::ExifData exif_data_;
    Exiv2::IptcData iptc_data_;
    Exiv2::XmpData xmp_data_;
};

// Glib::ustring get_xmp_sidecar_path(const Glib::ustring &path);
// Exiv2::Image::AutoPtr open_exiv2(const Glib::ustring &fname,
//                                  bool merge_xmp_sidecar);
// Exiv2::XmpData read_exiv2_xmp(const Glib::ustring &fname);

} // namespace rtengine
