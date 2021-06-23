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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <functional>

#include <strings.h>
#include <tiff.h>
#include <glib/gstdio.h>
#include <glibmm/convert.h>
#include <regex>
#include <sstream>

#include "imagedata.h"
#include "imagesource.h"
#include "rt_math.h"
#include "metadata.h"
#include "utils.h"

#pragma GCC diagnostic warning "-Wextra"
#define PRINT_HDR_PS_DETECTION 0

using namespace rtengine;

namespace
{

const std::string& validateUft8(const std::string& str, const std::string& on_error = "???")
{
    if (Glib::ustring(str).validate()) {
        return str;
    }

    return on_error;
}

}

namespace rtengine {

extern const Settings *settings;

} // namespace rtengine

FramesMetaData* FramesMetaData::fromFile(const Glib::ustring& fname)
{
    return new FramesData(fname);
}

FramesData::FramesData(const Glib::ustring &fname) :
    ok_(false),
    fname_(fname),
    dcrawFrameCount(0),
    time{},
    timeStamp{},
    iso_speed(0),
    aperture(0.),
    focal_len(0.),
    focal_len35mm(0.),
    focus_dist(0.f),
    shutter(0.),
    expcomp(0.),
    make("Unknown"),
    model("Unknown"),
    orientation("Unknown"),
    rating(0), // FIXME: Implement
    lens("Unknown"),
    sampleFormat(IIOSF_UNKNOWN),
    isPixelShift(false),
    isHDR(false),
    w_(-1),
    h_(-1)
{
    make.clear();
    model.clear();
    serial.clear();
    orientation.clear();
    lens.clear();

    try {
        Exiv2Metadata meta(fname);
        meta.load();
        const auto& exif = meta.exifData();
        ok_ = true;

        // taken and adapted from darktable (src/common/exif.cc)
/*
   This file is part of darktable,
   copyright (c) 2009--2013 johannes hanika.
   copyright (c) 2011 henrik andersson.
   copyright (c) 2012-2017 tobias ellinghaus.

   darktable is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   darktable is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with darktable.  If not, see <http://www.gnu.org/licenses/>.
 */

        Exiv2::ExifData::const_iterator pos;

        const auto find_exif_tag =
            [&exif, &pos](const std::string &name) -> bool
            {
                pos = exif.findKey(Exiv2::ExifKey(name));
                return pos != exif.end() && pos->size();
            };

        const auto find_tag =
            [&exif, &pos](decltype(Exiv2::make) func) -> bool
            {
                pos = func(exif);
                return pos != exif.end() && pos->size();
            };

        // List of tag names taken from exiv2's printSummary() in actions.cpp

        if (find_tag(Exiv2::make)) {
            make = validateUft8(pos->print(&exif)); // validateUft8 (#5923) still needed?
        }

        if (find_tag(Exiv2::model)) {
            model = validateUft8(pos->print(&exif)); // validateUft8 (#5923) still needed?
        }

        if (make.size() > 0) {
            for (const auto& corp : {
                "Canon",
                "NIKON",
                "EPSON",
                "KODAK",
                "Kodak",
                "OLYMPUS",
                "PENTAX",
                "RICOH",
                "MINOLTA",
                "Minolta",
                "Konica",
                "CASIO",
                "Sinar",
                "Phase One",
                "SAMSUNG",
                "Mamiya",
                "MOTOROLA",
                "Leaf",
                "Panasonic"
            }) {
                if (make.find(corp) != std::string::npos) { // Simplify company names
                    make = corp;
                    break;
                }
            }
        }

    /*
    TODO: Implement ratings in exiv2 situations. See PR #5325

    // Look for Rating metadata in the following order:
    // 1. EXIF
    // 2. XMP
    // 3. pp3 sidecar file
    tag = newFrameRootDir->findTagUpward("Rating");
    if (tag && tag->toInt() != 0) {
        rating = tag->toInt();
    }
    char sXMPRating[64];
    if (newFrameRootDir->getXMPTagValue("xmp:Rating", sXMPRating)) {
        // Guard against out-of-range values (<0, >5)
        rating = rtengine::max(0, rtengine::min(5, atoi(sXMPRating)));
        // Currently, Rating=-1 is not supported. A value of -1 should mean
        // "Rejected" according to the specification. Maybe in the future, Rating=-1
        // sets InTrash=true?
    }
    */

        std::string::size_type nonspace_pos = make.find_last_not_of(' ');
        if (nonspace_pos != std::string::npos && nonspace_pos + 1 < make.size()) {
            make.erase(nonspace_pos + 1);
        }
        nonspace_pos = model.find_last_not_of(' ');
        if (nonspace_pos != std::string::npos && nonspace_pos + 1 < model.size()) {
            model.erase(nonspace_pos + 1);
        }

        if (!make.empty() && model.find(make + ' ') == 0) {
            model.erase(0, make.size() + 1);
        }

        if (find_tag(Exiv2::exposureTime)) {
            shutter = pos->toFloat();
        }

        if (find_tag(Exiv2::fNumber)) {
            aperture = pos->toFloat();
        }

        // Read ISO speed - Nikon happens to return a pair for Lo and Hi modes
        if (find_tag(Exiv2::isoSpeed)) {
            // If standard exif iso tag, use the old way of interpreting the return value to be more regression-save
            if (pos->key() == "Exif.Photo.ISOSpeedRatings") {
                const long isofield = pos->count() > 1 ? 1 : 0;
                iso_speed = pos->toFloat(isofield);
            } else {
                iso_speed = std::atof(pos->print().c_str());
            }
        }
        // Some newer cameras support iso settings that exceed the 16 bit of exif's ISOSpeedRatings
        if (iso_speed == 65535 || iso_speed == 0) {
            if (find_exif_tag("Exif.PentaxDng.ISO") || find_exif_tag("Exif.Pentax.ISO")) {
                iso_speed = std::atof(pos->print().c_str());
            }
            else if (
                (
                    make == "SONY"
                    || make == "Canon"
                )
                && find_exif_tag("Exif.Photo.RecommendedExposureIndex")
            ) {
                iso_speed = pos->toFloat();
            }
        }

        if (find_tag(Exiv2::focalLength)) {
            // This works around a bug in exiv2 the developers refuse to fix
            // For details see http://dev.exiv2.org/issues/1083
            if (pos->key() == "Exif.Canon.FocalLength" && pos->count() == 4) {
                focal_len = pos->toFloat(1);
            } else {
                focal_len = pos->toFloat();
            }
        }

        if (find_exif_tag("Exif.Photo.FocalLengthIn35mmFilm")) {
            focal_len35mm = pos->toFloat();
        }

        if (find_tag(Exiv2::subjectDistance)) {
            focus_dist = (0.01 * std::pow(10, pos->toFloat() / 40));
        }

        if (find_tag(Exiv2::orientation)) {
            static const std::vector<std::string> ormap = {
                "Unknown",
                "Horizontal (normal)",
                "Mirror horizontal",
                "Rotate 180",
                "Mirror vertical",
                "Mirror horizontal and rotate 270 CW",
                "Rotate 90 CW",
                "Mirror horizontal and rotate 90 CW",
                "Rotate 270 CW",
                "Unknown"
            };
            auto idx = pos->toLong();
            if (idx >= 0 && idx < long(ormap.size())) {
                orientation = ormap[idx];
            }
            //orientation = pos->print(&exif);
        }

        if (find_tag(Exiv2::lensName)) {
            lens = validateUft8(pos->print(&exif)); // validateUft8 (#5923) still needed?
        } else if (find_exif_tag("Exif.Photo.LensSpecification") && pos->count() == 4) {
            const auto round =
                [](float f) -> float
                {
                    return int(f * 10.f + 0.5f) / 10.f;
                };
            float fl_lo = round(pos->toFloat(0));
            float fl_hi = round(pos->toFloat(1));
            float fn_lo = round(pos->toFloat(2));
            float fn_hi = round(pos->toFloat(3));
            std::ostringstream buf;
            buf << fl_lo;
            if (fl_lo < fl_hi) {
                buf << "-" << fl_hi;
            }
            buf << "mm F" << fn_lo;
            if (fn_lo < fn_hi) {
                buf << "-" << fn_hi;
            }
            lens = buf.str();
        }
        if (lens.empty() || lens.find_first_not_of('-') == std::string::npos) {
            lens = "Unknown";
        }

        std::string datetime_taken;
        if (find_exif_tag("Exif.Image.DateTimeOriginal")) {
            datetime_taken = pos->print(&exif);
        }
        else if (find_exif_tag("Exif.Photo.DateTimeOriginal")) {
            datetime_taken = pos->print(&exif);
        }
        else if (find_exif_tag("Exif.Photo.DateTimeDigitized")) {
            datetime_taken = pos->print(&exif);
        } else if (find_exif_tag("Exif.Image.DateTime")) {
            datetime_taken = validateUft8(pos->print(&exif));
        }
        if (sscanf(datetime_taken.c_str(), "%d:%d:%d %d:%d:%d", &time.tm_year, &time.tm_mon, &time.tm_mday, &time.tm_hour, &time.tm_min, &time.tm_sec) == 6) {
            time.tm_year -= 1900;
            time.tm_mon -= 1;
            time.tm_isdst = -1;
            timeStamp = mktime(&time);
        }

        if (find_exif_tag("Exif.Image.ExposureBiasValue")) {
            expcomp = pos->toFloat();
        } else if (find_exif_tag("Exif.Photo.ExposureBiasValue")) {
            expcomp = pos->toFloat();
        }

        if (find_exif_tag("Exif.Image.Rating")) {
            rating = pos->toLong();
        } else {
            auto it = meta.xmpData().findKey(Exiv2::XmpKey("Xmp.xmp.Rating"));
            if (it != meta.xmpData().end() && it->size()) {
                rating = it->toLong();
            }
        }

        // try getting some metadata from ImageDescription
        if (!make.compare(0, 5, "KODAK") && !getISOSpeed() && !getFNumber() && !getFocalLen() && !getShutterSpeed() &&
            find_exif_tag("Exif.Image.ImageDescription")) {
            std::string s = pos->toString();
            std::string line;
            std::smatch m;
            const auto d =
                [&m]() -> double {
                    std::string s = m[1];
                    return atof(s.c_str());
                };
            while (true) {
                auto p = s.find('\r');
                if (p == std::string::npos) {
                    break;
                }
                auto line = s.substr(0, p);
                s = s.substr(p+1);

                if (std::regex_match(line, m, std::regex("ISO: +([0-9]+) *"))) {
                    iso_speed = d();
                } else if (std::regex_match(line, m, std::regex("Aperture: +F([0-9.]+) *"))) {
                    aperture = d();
                } else if (std::regex_match(line, m, std::regex("Shutter: +([0-9.]+) *"))) {
                    shutter = d();
                    if (shutter) {
                        shutter = 1.0/shutter;
                    }
                } else if (std::regex_match(line, m, std::regex("Lens \\(mm\\): +([0-9.]+) *"))) {
                    focal_len = d();
                } else if (std::regex_match(line, m, std::regex("Exp Comp: +([0-9.]+) *"))) {
                    expcomp = d();
                }
            }
        }

        meta.getDimensions(w_, h_);

        // -----------------------
        // Special file type detection (HDR, PixelShift)
        // ------------------------
        uint16 bitspersample = 0, samplesperpixel = 0, sampleformat = 0, photometric = 0, compression = 0;
        const auto bps = exif.findKey(Exiv2::ExifKey("Exif.Image.BitsPerSample"));
        const auto spp = exif.findKey(Exiv2::ExifKey("Exif.Image.SamplesPerPixel"));
        const auto sf = exif.findKey(Exiv2::ExifKey("Exif.Image.SampleFormat"));
        const auto pi = exif.findKey(Exiv2::ExifKey("Exif.Image.PhotometricInterpretation"));
        const auto c = exif.findKey(Exiv2::ExifKey("Exif.Image.Compression"));

        if (
            !make.compare(0, 6, "PENTAX")
            || (
                !make.compare(0, 5, "RICOH")
                && !model.compare (0, 6, "PENTAX")
            )
        ) {
            if (find_exif_tag("Exif.Pentax.DriveMode")) {
                std::string buf = pos->toString(3);
                buf[3] = 0;
                if (buf == "HDR") {
                    isHDR = true;
#if PRINT_HDR_PS_DETECTION
                    printf("HDR detected ! -> DriveMode = \"HDR\"\n");
#endif
                }
            }

            if (
                !isHDR
                && (
                    find_exif_tag("Exif.Pentax.Quality")
                    || find_exif_tag("Exif.PentaxDng.Quality")
                )
                && (
                    pos->toLong() == 7
                    || pos->toLong() == 8
                )
            ) {
                isPixelShift = true;
#if PRINT_HDR_PS_DETECTION
                printf("PixelShift detected ! -> \"Quality\" = 7\n");
#endif
            }
        }

        sampleFormat = IIOSF_UNKNOWN;

        if (sf == exif.end())
            /*
             * WARNING: This is a dirty hack!
             * We assume that files which doesn't contain the TIFFTAG_SAMPLEFORMAT tag
             * (which is the case with uncompressed TIFFs produced by RT!) are RGB files,
             * but that may be not true.   --- Hombre
             */
        {
            sampleformat = SAMPLEFORMAT_UINT;
        } else {
            sampleformat = sf->toLong();
        }

        if (bps == exif.end() || spp == exif.end() || pi == exif.end()) {
            return;
        }

        bitspersample = bps->toLong();
        samplesperpixel = spp->toLong();

        photometric = pi->toLong();
        if (photometric == PHOTOMETRIC_LOGLUV) {
            if (c == exif.end()) {
                compression = COMPRESSION_NONE;
            } else {
                compression = c->toLong();
            }
        }

        if (photometric == PHOTOMETRIC_RGB || photometric == PHOTOMETRIC_MINISBLACK) {
            if (sampleformat == SAMPLEFORMAT_INT || sampleformat == SAMPLEFORMAT_UINT) {
                if (bitspersample == 8) {
                    sampleFormat = IIOSF_UNSIGNED_CHAR;
                } else if (bitspersample <= 16) {
                    sampleFormat = IIOSF_UNSIGNED_SHORT;
                }
            } else if (sampleformat == SAMPLEFORMAT_IEEEFP) {
                if (bitspersample==16) {
                    sampleFormat = IIOSF_FLOAT16;
                    isHDR = true;
#if PRINT_HDR_PS_DETECTION
                    printf("HDR detected ! -> sampleFormat = %d   (16-bit)\n", sampleFormat);
#endif
                }
                else if (bitspersample == 24) {
                    sampleFormat = IIOSF_FLOAT24;
                    isHDR = true;
#if PRINT_HDR_PS_DETECTION
                    printf("HDR detected ! -> sampleFormat = %d   (24-bit)\n", sampleFormat);
#endif
                }
                else if (bitspersample == 32) {
                    sampleFormat = IIOSF_FLOAT32;
                    isHDR = true;
#if PRINT_HDR_PS_DETECTION
                    printf("HDR detected ! -> sampleFormat = %d   (32-bit)\n", sampleFormat);
#endif
                }
            }
        } else if (photometric == PHOTOMETRIC_CFA) {
            if (sampleformat == SAMPLEFORMAT_IEEEFP) {
                if (bitspersample == 16) {
                    sampleFormat = IIOSF_FLOAT16;
                    isHDR = true;
#if PRINT_HDR_PS_DETECTION
                    printf("HDR detected ! -> sampleFormat = %d   (16-bit)\n", sampleFormat);
#endif
                }
                else if (bitspersample == 24) {
                    sampleFormat = IIOSF_FLOAT24;
                    isHDR = true;
#if PRINT_HDR_PS_DETECTION
                    printf("HDR detected ! -> sampleFormat = %d   (24-bit)\n", sampleFormat);
#endif
                }
                else if (bitspersample == 32) {
                    sampleFormat = IIOSF_FLOAT32;
                    isHDR = true;
#if PRINT_HDR_PS_DETECTION
                    printf("HDR detected ! -> sampleFormat = %d   (32-bit)\n", sampleFormat);
#endif
                }
            } else if (sampleformat == SAMPLEFORMAT_INT || sampleformat == SAMPLEFORMAT_UINT) {
                if (bitspersample == 8) {   // shouldn't occur...
                    sampleFormat = IIOSF_UNSIGNED_CHAR;
                } else if (bitspersample <= 16) {
                    sampleFormat = IIOSF_UNSIGNED_SHORT;
                }
            }
        } else if (photometric == 34892 || photometric == 32892  /* Linear RAW (see DNG spec ; 32892 seem to be a flaw from Sony's ARQ files) */) {
            if (sampleformat == SAMPLEFORMAT_IEEEFP) {
                sampleFormat = IIOSF_FLOAT32;
                isHDR = true;
#if PRINT_HDR_PS_DETECTION
                printf("HDR detected ! -> sampleFormat = %d\n", sampleFormat);
#endif
            } else if (sampleformat == SAMPLEFORMAT_INT || sampleformat == SAMPLEFORMAT_UINT) {
                if (bitspersample == 8) {   // shouldn't occur...
                    sampleFormat = IIOSF_UNSIGNED_CHAR;
                } else if (bitspersample <= 16) {
                    sampleFormat = IIOSF_UNSIGNED_SHORT;
                    if (find_exif_tag("Exif.Photo.MakerNote") && (!make.compare (0, 4, "SONY")) && bitspersample >= 12 && samplesperpixel == 4) {
                        isPixelShift = true;
#if PRINT_HDR_PS_DETECTION
                        printf("PixelShift detected ! -> \"Make\" = SONY, bitsPerPixel > 8, samplesPerPixel == 4\n");
#endif
                    }
                }
            }
        } else if (photometric == PHOTOMETRIC_LOGLUV) {
            if (compression == COMPRESSION_SGILOG24) {
                sampleFormat = IIOSF_LOGLUV24;
                isHDR = true;
#if PRINT_HDR_PS_DETECTION
                printf("HDR detected ! -> sampleFormat = %d\n", sampleFormat);
#endif
            } else if (compression == COMPRESSION_SGILOG) {
                sampleFormat = IIOSF_LOGLUV32;
                isHDR = true;
#if PRINT_HDR_PS_DETECTION
                printf("HDR detected ! -> sampleFormat = %d\n", sampleFormat);
#endif
            }
        }
    } catch (const std::exception& e) {
        if (settings->verbose) {
            std::cerr << "EXIV2 ERROR: " << e.what() << std::endl;
        }
        ok_ = false;
    }
}

bool FramesData::getPixelShift() const
{
    return isPixelShift;
}

bool FramesData::getHDR() const
{
    return isHDR;
}

std::string FramesData::getImageType() const
{
    return isPixelShift ? "PS" : isHDR ? "HDR" : "STD";
}

IIOSampleFormat FramesData::getSampleFormat() const
{
    return sampleFormat;
}

bool FramesData::hasExif() const
{
    return ok_;
}

tm FramesData::getDateTime() const
{
    return time;
}

time_t FramesData::getDateTimeAsTS() const
{
    return timeStamp;
}

int FramesData::getISOSpeed() const
{
    return iso_speed;
}

double FramesData::getFNumber() const
{
    return aperture;
}

double FramesData::getFocalLen() const
{
    return focal_len;
}

double FramesData::getFocalLen35mm() const
{
    return focal_len35mm;
}

float FramesData::getFocusDist() const
{
    return focus_dist;
}


double FramesData::getShutterSpeed() const
{
    return shutter;
}


double FramesData::getExpComp() const
{
    return expcomp;
}


std::string FramesData::getMake() const
{
    return make;
}


std::string FramesData::getModel() const
{
    return model;
}


std::string FramesData::getLens() const
{
    return lens;
}


std::string FramesData::getSerialNumber() const
{
    return serial;
}


std::string FramesData::getOrientation() const
{
    return orientation;
}


void FramesData::setDCRawFrameCount(unsigned int frameCount)
{
    dcrawFrameCount = frameCount;
}

unsigned int FramesData::getFrameCount() const
{
    return std::max(1U, dcrawFrameCount);
}


Glib::ustring FramesData::getFileName() const
{
    return fname_;
}


int FramesData::getRating() const
{
    return rating;
}

//------inherited functions--------------//

std::string FramesMetaData::apertureToString(double aperture)
{
	// TODO: Replace sprintf()
    char buffer[256];
    snprintf(buffer, sizeof(buffer), "%0.1f", aperture);
    return buffer;
}

std::string FramesMetaData::shutterToString(double shutter)
{
    char buffer[256];

    if (shutter > 0.0 && shutter <= 0.5) {
        snprintf(buffer, sizeof(buffer), "1/%0.0f", 1.0 / shutter);
    } else if (int(shutter) == shutter) {
        snprintf(buffer, sizeof(buffer), "%d", int(shutter));
    } else {
        snprintf(buffer, sizeof(buffer), "%0.1f", shutter);
    }

    return buffer;
}

std::string FramesMetaData::expcompToString(double expcomp, bool maskZeroexpcomp)
{

    char buffer[256];

    if (maskZeroexpcomp) {
        if (expcomp != 0.0) {
            snprintf(buffer, sizeof(buffer), "%0.2f", expcomp);
            return buffer;
        } else {
            return "";
        }
    } else {
        snprintf(buffer, sizeof(buffer), "%0.2f", expcomp);
        return buffer;
    }
}

double FramesMetaData::shutterFromString(std::string s)
{
    const std::string::size_type i = s.find_first_of('/');

    if (i == std::string::npos) {
        return std::atof(s.c_str());
    } else {
        const double denominator = std::atof(s.substr(i + 1).c_str());
        return
            denominator
                ? std::atof(s.substr(0, i).c_str()) / denominator
                : 0.0;
    }
}

double FramesMetaData::apertureFromString(std::string s)
{

    return std::atof(s.c_str());
}


namespace {

template<class T>
void set_exif(Exiv2::ExifData &exif, const std::string &key, T val)
{
    try {
        exif[key] = val;
    } catch (std::exception &exc) {
        if (settings->verbose) {
            std::cout << "Exif -- error setting " << key << " to " << val << ": " << exc.what() << std::endl;
        }
    }
}

} // namespace

void FramesData::fillBasicTags(Exiv2::ExifData &exif) const
{
    if (!hasExif()) {
        return;
    }
    set_exif(exif, "Exif.Photo.ISOSpeedRatings", getISOSpeed());
    set_exif(exif, "Exif.Photo.FNumber", Exiv2::URationalValue(Exiv2::URational(round(getFNumber() * 10), 10)));
    auto s = shutterToString(getShutterSpeed());
    auto p = s.find('.');
    if (p != std::string::npos) {
        assert(p == s.length()-2);
        s = s.substr(0, p) + s.substr(p+1) + "/10";
    } else if (s.find('/') == std::string::npos) {
        s += "/1";
    }
    set_exif(exif, "Exif.Photo.ExposureTime", s);
    set_exif(exif, "Exif.Photo.FocalLength", Exiv2::URationalValue(Exiv2::URational(getFocalLen() * 10, 10)));
    set_exif(exif, "Exif.Photo.ExposureBiasValue", Exiv2::RationalValue(Exiv2::Rational(round(getExpComp() * 100), 100)));
    set_exif(exif, "Exif.Image.Make", getMake());
    set_exif(exif, "Exif.Image.Model", getModel());
    set_exif(exif, "Exif.Photo.LensModel", getLens());
    char buf[256];
    auto t = getDateTime();
    strftime(buf, 256, "%Y:%m:%d %H:%M:%S", &t);
    set_exif(exif, "Exif.Photo.DateTimeOriginal", buf);
}


void FramesData::getDimensions(int &w, int &h) const
{
    w = w_;
    h = h_;
}


void FramesData::setDimensions(int w, int h)
{
    w_ = w;
    h_ = h;
}
