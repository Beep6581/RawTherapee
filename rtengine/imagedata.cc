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
#include <strings.h>
#include <glib/gstdio.h>
#include <tiff.h>

#include "imagedata.h"
#include "iptcpairs.h"
#include "imagesource.h"
#include "rt_math.h"
#pragma GCC diagnostic warning "-Wextra"
#define PRINT_HDR_PS_DETECTION 0

using namespace rtengine;

extern "C" IptcData *iptc_data_new_from_jpeg_file(FILE* infile);

namespace
{

Glib::ustring to_utf8(const std::string& str)
{
    try {
        return Glib::locale_to_utf8(str);
    } catch (Glib::Error&) {
        return Glib::convert_with_fallback(str, "UTF-8", "ISO-8859-1", "?");
    }
}

}

FramesMetaData* FramesMetaData::fromFile(const Glib::ustring& fname, std::unique_ptr<RawMetaDataLocation> rml, bool firstFrameOnly)
{
    return new FramesData(fname, std::move(rml), firstFrameOnly);
}

FrameData::FrameData(rtexif::TagDirectory* frameRootDir_, rtexif::TagDirectory* rootDir, rtexif::TagDirectory* firstRootDir)
    : frameRootDir(frameRootDir_), iptc(nullptr), time(), timeStamp(), iso_speed(0), aperture(0.), focal_len(0.), focal_len35mm(0.), focus_dist(0.f),
      shutter(0.), expcomp(0.), make("Unknown"), model("Unknown"), orientation("Unknown"), lens("Unknown"),
      sampleFormat(IIOSF_UNKNOWN), isPixelShift(false), isHDR(false)
{
    memset(&time, 0, sizeof(time));

    if (!frameRootDir) {
        return;
    }

    rtexif::Tag* tag;
    rtexif::TagDirectory* newFrameRootDir = frameRootDir;

    memset(&time, 0, sizeof(time));
    timeStamp = 0;
    iso_speed = 0;
    aperture = 0.0;
    focal_len = 0.0;
    focal_len35mm = 0.0;
    focus_dist = 0.0f;
    shutter = 0.0;
    expcomp = 0.0;
    make.clear();
    model.clear();
    serial.clear();
    orientation.clear();
    lens.clear();

    tag = newFrameRootDir->findTag("Make");

    if (!tag) {
        newFrameRootDir = rootDir;
        tag = newFrameRootDir->findTag("Make");

        if (!tag) {
            // For some raw files (like Canon's CR2 files), the metadata are contained in the first root directory
            newFrameRootDir = firstRootDir;
            tag = newFrameRootDir->findTag("Make");
        }
    }

    if (tag) {
        make = tag->valueToString();

        // Same dcraw treatment
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

        make.erase(make.find_last_not_of(' ') + 1);
    }

    tag = newFrameRootDir->findTagUpward("Model");

    if (tag) {
        model = tag->valueToString();
    }

    if (!model.empty()) {
        std::string::size_type i = 0;

        if (
            make.find("KODAK") != std::string::npos
            && (
                (i = model.find(" DIGITAL CAMERA")) != std::string::npos
                || (i = model.find(" Digital Camera")) !=  std::string::npos
                || (i = model.find("FILE VERSION")) !=  std::string::npos
            )
        ) {
            model.resize(i);
        }

        model.erase(model.find_last_not_of(' ') + 1);

        if (!strncasecmp(model.c_str(), make.c_str(), make.size())) {
            if (model.size() >= make.size() && model[make.size()] == ' ') {
                model.erase(0, make.size() + 1);
            }
        }

        if (model.find("Digital Camera ") != std::string::npos) {
            model.erase(0, 15);
        }
    } else {
        model = "Unknown";
    }

    tag = newFrameRootDir->findTagUpward("Orientation");

    if (tag) {
        orientation = tag->valueToString();
    }

    tag = newFrameRootDir->findTagUpward("MakerNote");
    rtexif::TagDirectory* mnote = nullptr;

    if (tag) {
        mnote = tag->getDirectory();
    }

    rtexif::TagDirectory* exif = nullptr;
    tag = newFrameRootDir->findTagUpward("Exif");

    if (tag) {
        exif = tag->getDirectory();
    }

    if (exif) {

        // standard exif tags
        if ((tag = exif->getTag("ShutterSpeedValue"))) {
            shutter = tag->toDouble();
        }

        if ((tag = exif->getTag("ExposureTime"))) {
            shutter = tag->toDouble();
        }

        if ((tag = exif->getTag("ApertureValue"))) {
            aperture = tag->toDouble();
        }

        if ((tag = exif->getTag("FNumber"))) {
            aperture = tag->toDouble();
        }

        if ((tag = exif->getTag("ExposureBiasValue"))) {
            expcomp = tag->toDouble();
        }

        if ((tag = exif->getTag("FocalLength"))) {
            focal_len = tag->toDouble();
        }

        if ((tag = exif->getTag("FocalLengthIn35mmFilm"))) {
            focal_len35mm = tag->toDouble();
        }

        // Focus distance from EXIF or XMP. MakerNote ones are scattered and partly encrypted
        int num = -3, denom = -3;

        // First try, official EXIF. Set by Adobe on some DNGs
        tag = exif->getTag("SubjectDistance");

        if (tag) {
            int num, denom;
            tag->toRational(num, denom);
        } else {
            // Second try, XMP data
            char sXMPVal[64];

            if (newFrameRootDir->getXMPTagValue("aux:ApproximateFocusDistance", sXMPVal)) {
                sscanf(sXMPVal, "%d/%d", &num, &denom);
            }
        }

        if (num != -3) {
            if ((denom == 1 && num >= 10000) || num < 0 || denom < 0) {
                focus_dist = 10000;    // infinity
            } else if (denom > 0) {
                focus_dist = (float)num / denom;
            }
        }

        if ((tag = exif->getTag("ISOSpeedRatings"))) {
            iso_speed = tag->toDouble();
        }

        if ((tag = exif->getTag("DateTimeOriginal"))) {
            if (sscanf((const char*)tag->getValue(), "%d:%d:%d %d:%d:%d", &time.tm_year, &time.tm_mon, &time.tm_mday, &time.tm_hour, &time.tm_min, &time.tm_sec) == 6) {
                time.tm_year -= 1900;
                time.tm_mon -= 1;
                time.tm_isdst = -1;
                timeStamp = mktime(&time);
            }
        }

        tag = exif->findTag("SerialNumber");

        if (!tag) {
            tag = exif->findTag("InternalSerialNumber");
        }

        if (tag) {
            serial = tag->valueToString();
        }

        // guess lens...
        lens = "Unknown";

        // Sometimes (e.g. DNG) EXIF already contains lens data

        if (!make.compare(0, 8, "FUJIFILM")) {
            if (exif->getTag("LensModel")) {
                lens = exif->getTag("LensModel")->valueToString();
            }
        } else if (!make.compare(0, 4, "SONY")) {
            if (iso_speed == 65535 || iso_speed == 0) {
                rtexif::Tag* isoTag = exif->getTag("RecommendedExposureIndex");

                if (isoTag) {
                    iso_speed = isoTag->toDouble();
                }
            }

        }

        if (lens == "Unknown") {

            if (mnote) {

                if (!make.compare(0, 5, "NIKON")) {
                    // ISO at max value supported, check manufacturer specific
                    if (iso_speed == 65535 || iso_speed == 0) {
                        rtexif::Tag* isoTag = mnote->getTagP("ISOInfo/ISO");

                        if (isoTag) {
                            iso_speed = isoTag->toInt();
                        }
                    }

                    bool lensOk = false;

                    if (mnote->getTag("LensData")) {
                        std::string ldata = mnote->getTag("LensData")->valueToString();
                        size_t pos;

                        if (ldata.size() > 10 && (pos = ldata.find("Lens = ")) != Glib::ustring::npos) {
                            lens = ldata.substr(pos + 7);

                            if (lens.compare(0, 7, "Unknown")) {
                                lensOk = true;
                            } else {
                                size_t pos = lens.find("$FL$");        // is there a placeholder for focallength?

                                if (pos != Glib::ustring::npos) {               // then fill in focallength
                                    lens = lens.replace(pos, 4, exif->getTag("FocalLength")->valueToString());

                                    if (mnote->getTag("LensType")) {
                                        std::string ltype = mnote->getTag("LensType")->valueToString();

                                        if (ltype.find("MF = Yes") != Glib::ustring::npos) { // check, whether it's a MF lens, should be always
                                            lens = lens.replace(0, 7, "MF");
                                        }

                                        lensOk = true;
                                    }
                                }
                            }
                        }
                    }

                    if (!lensOk && mnote->getTag("Lens")) {
                        std::string ldata = mnote->getTag("Lens")->valueToString();
                        size_t i = 0, j = 0;
                        double n[4] = {0.0};

                        for (int m = 0; m < 4; m++) {
                            while (i < ldata.size() && ldata[i] != '/') {
                                i++;
                            }

                            int nom = atoi(ldata.substr(j, i).c_str());
                            j = i + 1;
                            i++;

                            while (i < ldata.size() && ldata[i] != ',') {
                                i++;
                            }

                            int den = atoi(ldata.substr(j, i).c_str());
                            j = i + 2;
                            i += 2;
                            n[m] = (double) nom / std::max(den, 1);
                        }

                        std::ostringstream str;

                        if (n[0] == n[1]) {
                            str << "Unknown " << n[0] << "mm F/" << n[2];
                        } else if (n[2] == n[3]) {
                            str << "Unknown " << n[0] << "-" << n[1] << "mm F/" << n[2];
                        } else {
                            str << "Unknown " << n[0] << "-" << n[1] << "mm F/" << n[2] << "-" << n[3];
                        }

                        lens = str.str();

                        // Look whether it's MF or AF
                        if (mnote->getTag("LensType")) {
                            std::string ltype = mnote->getTag("LensType")->valueToString();

                            if (ltype.find("MF = Yes") != Glib::ustring::npos) { // check, whether it's a MF lens
                                lens = lens.replace(0, 7, "MF");    // replace 'Unknwon' with 'MF'
                            } else {
                                lens = lens.replace(0, 7, "AF");    // replace 'Unknwon' with 'AF'
                            }
                        }
                    }
                } else if (!make.compare(0, 5, "Canon")) {
                    // ISO at max value supported, check manufacturer specific
                    if (iso_speed == 65535 || iso_speed == 0) {
                        rtexif::Tag* baseIsoTag = mnote->getTagP("CanonShotInfo/BaseISO");

                        if (baseIsoTag) {
                            iso_speed = baseIsoTag->toInt();
                        }
                    }

                    int found = false;
                    // canon EXIF have a string for lens model
                    rtexif::Tag *lt = mnote->getTag("LensType");

                    if (lt) {
                        std::string ldata = lt->valueToString();

                        if (ldata.size() > 1) {
                            found = true;
                            lens = "Canon " + ldata;
                        }
                    }

                    if (!found || lens.substr(lens.find(' ')).length() < 7) {
                        lt = mnote->findTag("LensID");

                        if (lt) {
                            std::string ldata = lt->valueToString();

                            if (ldata.size() > 1) {
                                lens = ldata;
                            }
                        }
                    }
                } else if (!make.compare(0, 6, "PENTAX") || (!make.compare(0, 5, "RICOH") && !model.compare(0, 6, "PENTAX"))) {
                    // ISO at max value supported, check manufacturer specific
                    if (iso_speed == 65535 || iso_speed == 0) {
                        rtexif::Tag* baseIsoTag = mnote->getTag("ISO");

                        if (baseIsoTag) {
                            std::string isoData = baseIsoTag->valueToString();

                            if (isoData.size() > 1) {
                                iso_speed = stoi(isoData);
                            }
                        }
                    }

                    if (mnote->getTag("LensType")) {
                        lens = mnote->getTag("LensType")->valueToString();
                    }

                    // Try to get the FocalLength from the LensInfo structure, where length below 10mm will be correctly set
                    rtexif::Tag* flt = mnote->getTagP("LensInfo/FocalLength");

                    if (flt) {
                        // Don't replace Exif focal_len if Makernotes focal_len is 0
                        if (flt->toDouble() > 0) {
                            focal_len = flt->toDouble();
                        }
                    } else if ((flt = mnote->getTagP("FocalLength"))) {
                        rtexif::Tag* flt = mnote->getTag("FocalLength");
                        focal_len = flt->toDouble();
                    }

                    if (mnote->getTag("FocalLengthIn35mmFilm")) {
                        focal_len35mm = mnote->getTag("FocalLengthIn35mmFilm")->toDouble();
                    }
                } else if (mnote && (!make.compare(0, 4, "SONY") || !make.compare(0, 6, "KONICA"))) {
                    if (mnote->getTag("LensID")) {
                        lens = mnote->getTag("LensID")->valueToString();
                    }
                } else if (!make.compare(0, 7, "OLYMPUS")) {
                    if (mnote->getTag("Equipment"))  {
                        rtexif::TagDirectory* eq = mnote->getTag("Equipment")->getDirectory();

                        if (eq->getTag("LensType")) {
                            lens = eq->getTag("LensType")->valueToString();
                        }
                    }
                } else if (mnote && !make.compare(0, 9, "Panasonic")) {
                    if (mnote->getTag("LensType")) {
                        std::string panalens = mnote->getTag("LensType")->valueToString();

                        if (panalens.find("LUMIX") != Glib::ustring::npos) {
                            lens = "Panasonic " + panalens;
                        } else {
                            lens = panalens;
                        }
                    }
                }
            } else if (exif->getTag("DNGLensInfo")) {
                lens = exif->getTag("DNGLensInfo")->valueToString();
            } else if (exif->getTag("LensModel")) {
                lens = exif->getTag("LensModel")->valueToString();
            } else if (exif->getTag("LensInfo")) {
                lens = exif->getTag("LensInfo")->valueToString();
            }
        }
    }

    rtexif::Tag* t = newFrameRootDir->getTag(0x83BB);

    if (t) {
        iptc = iptc_data_new_from_data((unsigned char*)t->getValue(), (unsigned)t->getValueSize());
    }


    // -----------------------  Special file type detection (HDR, PixelShift) ------------------------


    uint16 bitspersample = 0, samplesperpixel = 0, sampleformat = 0, photometric = 0, compression = 0;
    const rtexif::Tag* const bps = frameRootDir->findTag("BitsPerSample");
    const rtexif::Tag* const spp = frameRootDir->findTag("SamplesPerPixel");
    const rtexif::Tag* const sf = frameRootDir->findTag("SampleFormat");
    const rtexif::Tag* const pi = frameRootDir->findTag("PhotometricInterpretation");
    const rtexif::Tag* const c = frameRootDir->findTag("Compression");

    if (mnote && (!make.compare(0, 6, "PENTAX") || (!make.compare(0, 5, "RICOH") && !model.compare(0, 6, "PENTAX")))) {
        const rtexif::Tag* const hdr = mnote->findTag("HDR");

        if (hdr) {
            if (hdr->toInt() > 0) {
                isHDR = true;
#if PRINT_HDR_PS_DETECTION
                printf("HDR detected ! -> \"HDR\" tag found\n");
#endif
            }
        } else {
            const rtexif::Tag* const dm = mnote->findTag("DriveMode");

            if (dm) {
                char buffer[60];
                dm->toString(buffer, 3);
                buffer[3] = 0;

                if (!strcmp(buffer, "HDR")) {
                    isHDR = true;
#if PRINT_HDR_PS_DETECTION
                    printf("HDR detected ! -> DriveMode = \"HDR\"\n");
#endif
                }
            }
        }

        if (!isHDR) {
            const rtexif::Tag* const q = mnote->findTag("Quality");
            if (q && (q->toInt() == 7 || q->toInt() == 8)) {
                isPixelShift = true;
#if PRINT_HDR_PS_DETECTION
                printf("PixelShift detected ! -> \"Quality\" = 7\n");
#endif
            }
        }
    }

    sampleFormat = IIOSF_UNKNOWN;

    if (!sf)
        /*
         * WARNING: This is a dirty hack!
         * We assume that files which doesn't contain the TIFFTAG_SAMPLEFORMAT tag
         * (which is the case with uncompressed TIFFs produced by RT!) are RGB files,
         * but that may be not true.   --- Hombre
         */
    {
        sampleformat = SAMPLEFORMAT_UINT;
    } else {
        sampleformat = sf->toInt();
    }

    if (
        !bps
        || !spp
        || !pi
    ) {
        return;
    }

    bitspersample = bps->toInt();
    samplesperpixel = spp->toInt();

    photometric = pi->toInt();

    if (photometric == PHOTOMETRIC_LOGLUV) {
        if (!c) {
            compression = COMPRESSION_NONE;
        } else {
            compression = c->toInt();
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

                if (mnote && (!make.compare(0, 4, "SONY")) && bitspersample >= 12 && samplesperpixel == 4) {
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
}

FrameData::~FrameData()
{

    if (iptc) {
        iptc_data_free(iptc);
    }
}

procparams::IPTCPairs FrameData::getIPTCData() const
{
    return getIPTCData(iptc);
}

procparams::IPTCPairs FrameData::getIPTCData(IptcData* iptc_)
{

    procparams::IPTCPairs iptcc;

    if (!iptc_) {
        return iptcc;
    }

    unsigned char buffer[2100];

    for (int i = 0; i < 16; i++) {
        IptcDataSet* ds = iptc_data_get_next_dataset(iptc_, nullptr, IPTC_RECORD_APP_2, strTags[i].tag);

        if (ds) {
            iptc_dataset_get_data(ds, buffer, 2100);
            std::vector<Glib::ustring> icValues;
            icValues.push_back(to_utf8((char*)buffer));

            iptcc[strTags[i].field] = icValues;
            iptc_dataset_unref(ds);
        }
    }

    IptcDataSet* ds = nullptr;
    std::vector<Glib::ustring> keywords;

    while ((ds = iptc_data_get_next_dataset(iptc_, ds, IPTC_RECORD_APP_2, IPTC_TAG_KEYWORDS))) {
        iptc_dataset_get_data(ds, buffer, 2100);
        keywords.push_back(to_utf8((char*)buffer));
    }

    iptcc["Keywords"] = keywords;
    ds = nullptr;
    std::vector<Glib::ustring> suppCategories;

    while ((ds = iptc_data_get_next_dataset(iptc_, ds, IPTC_RECORD_APP_2, IPTC_TAG_SUPPL_CATEGORY))) {
        iptc_dataset_get_data(ds, buffer, 2100);
        suppCategories.push_back(to_utf8((char*)buffer));
        iptc_dataset_unref(ds);
    }

    iptcc["SupplementalCategories"] = suppCategories;
    return iptcc;
}


bool FrameData::getPixelShift() const
{
    return isPixelShift;
}
bool FrameData::getHDR() const
{
    return isHDR;
}
std::string FrameData::getImageType () const
{
    return isPixelShift ? "PS" : isHDR ? "HDR" : "STD";
}
IIOSampleFormat FrameData::getSampleFormat() const
{
    return sampleFormat;
}
rtexif::TagDirectory* FrameData::getExifData() const
{
    return frameRootDir;
}
bool FrameData::hasExif() const
{
    return frameRootDir && frameRootDir->getCount();
}
bool FrameData::hasIPTC() const
{
    return iptc;
}
tm FrameData::getDateTime() const
{
    return time;
}
time_t FrameData::getDateTimeAsTS() const
{
    return timeStamp;
}
int FrameData::getISOSpeed() const
{
    return iso_speed;
}
double FrameData::getFNumber() const
{
    return aperture;
}
double FrameData::getFocalLen() const
{
    return focal_len;
}
double FrameData::getFocalLen35mm() const
{
    return focal_len35mm;
}
float FrameData::getFocusDist() const
{
    return focus_dist;
}
double FrameData::getShutterSpeed() const
{
    return shutter;
}
double FrameData::getExpComp() const
{
    return expcomp;
}
std::string FrameData::getMake() const
{
    return make;
}
std::string FrameData::getModel() const
{
    return model;
}
std::string FrameData::getLens() const
{
    return lens;
}
std::string FrameData::getSerialNumber() const
{
    return serial;
}
std::string FrameData::getOrientation() const
{
    return orientation;
}



void FramesData::setDCRawFrameCount(unsigned int frameCount)
{
    dcrawFrameCount = frameCount;
}

unsigned int FramesData::getRootCount() const
{
    return roots.size();
}

unsigned int FramesData::getFrameCount() const
{
    return dcrawFrameCount ? dcrawFrameCount : frames.size();
}

FrameData *FramesData::getFrameData(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size() ? nullptr : frames.at(frame);
}

bool FramesData::getPixelShift () const
{
    // So far only Pentax and Sony provide multi-frame Pixel Shift files.
    // Only the first frame contains the Pixel Shift tag
    // If more brand have to be supported, this rule may need
    // to evolve

    return frames.empty() ? false : frames.at(0)->getPixelShift ();
}
bool FramesData::getHDR(unsigned int frame) const
{
    // So far only Pentax provides multi-frame HDR file.
    // Only the first frame contains the HDR tag
    // If more brand have to be supported, this rule may need
    // to evolve

    return frames.empty() || frame >= frames.size()  ? false : frames.at(0)->getHDR();
}

std::string FramesData::getImageType (unsigned int frame) const
{
    return frames.empty() || frame >= frames.size() ? "STD" : frames.at(0)->getImageType();
}

IIOSampleFormat FramesData::getSampleFormat(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ? IIOSF_UNKNOWN : frames.at(frame)->getSampleFormat();
}

rtexif::TagDirectory* FramesData::getFrameExifData(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ? nullptr : frames.at(frame)->getExifData();
}

rtexif::TagDirectory* FramesData::getBestExifData(ImageSource *imgSource, procparams::RAWParams *rawParams) const
{
    rtexif::TagDirectory *td = nullptr;

    if (frames.empty()) {
        return nullptr;
    }

    if (imgSource && rawParams) {
        eSensorType sensorType = imgSource->getSensorType();
        unsigned int imgNum = 0;

        if (sensorType == ST_BAYER) {
            imgNum = rtengine::LIM<unsigned int>(rawParams->bayersensor.imageNum, 0, frames.size() - 1);
            /*
            // might exist someday ?
            } else if (sensorType == ST_FUJI_XTRANS) {
                imgNum = rtengine::LIM<unsigned int>(rawParams->xtranssensor.imageNum, 0, frames.size() - 1);
            } else if (sensorType == ST_NONE && !imgSource->isRAW()) {
                // standard image multiframe support should come here (when implemented in GUI)
            */
        }

        td = getFrameExifData(imgNum);
        rtexif::Tag* makeTag;

        if (td && (makeTag = td->findTag("Make", true))) {
            td = makeTag->getParent();
        } else {
            td = getRootExifData(0);
        }
    }

    return td;
}

rtexif::TagDirectory* FramesData::getRootExifData(unsigned int root) const
{
    return roots.empty() || root >= roots.size()  ? nullptr : roots.at(root);
}

procparams::IPTCPairs FramesData::getIPTCData(unsigned int frame) const
{
    if (frame < frames.size() && frames.at(frame)->hasIPTC()) {
        return frames.at(frame)->getIPTCData();
    } else {
        if (iptc) {
            return FrameData::getIPTCData(iptc);
        } else {
            procparams::IPTCPairs emptyPairs;
            return emptyPairs;
        }
    }
}

bool FramesData::hasExif(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ? false : frames.at(frame)->hasExif();
}
bool FramesData::hasIPTC(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ?  false : frames.at(frame)->hasIPTC();
}

tm FramesData::getDateTime(unsigned int frame) const
{
    if (frames.empty() || frame >= frames.size()) {
        return {};
    } else {
        return frames.at(frame)->getDateTime();
    }
}
time_t FramesData::getDateTimeAsTS(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ? 0 : frames.at(frame)->getDateTimeAsTS();
}
int FramesData::getISOSpeed(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ? 0 : frames.at(frame)->getISOSpeed();
}
double FramesData::getFNumber(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ? 0. : frames.at(frame)->getFNumber();
}
double FramesData::getFocalLen(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ? 0. : frames.at(frame)->getFocalLen();
}
double FramesData::getFocalLen35mm(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ? 0. : frames.at(frame)->getFocalLen35mm();
}
float FramesData::getFocusDist(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ? 0.f : frames.at(frame)->getFocusDist();
}
double FramesData::getShutterSpeed(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ? 0. : frames.at(frame)->getShutterSpeed();
}
double FramesData::getExpComp(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ? 0. : frames.at(frame)->getExpComp();
}
std::string FramesData::getMake(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ? std::string() : frames.at(frame)->getMake();
}
std::string FramesData::getModel(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ? std::string() : frames.at(frame)->getModel();
}
std::string FramesData::getLens(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ? std::string() : frames.at(frame)->getLens();
}
std::string FramesData::getSerialNumber(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ? std::string() : frames.at(frame)->getSerialNumber();
}
std::string FramesData::getOrientation(unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ? std::string() : frames.at(frame)->getOrientation();
}


//------inherited functions--------------//


std::string FramesMetaData::apertureToString(double aperture)
{

    char buffer[256];
    sprintf(buffer, "%0.1f", aperture);
    return buffer;
}

std::string FramesMetaData::shutterToString(double shutter)
{

    char buffer[256];

    if (shutter > 0.0 && shutter <= 0.5) {
        sprintf(buffer, "1/%0.0f", 1.0 / shutter);
    } else {
        sprintf(buffer, "%0.1f", shutter);
    }

    return buffer;
}

std::string FramesMetaData::expcompToString(double expcomp, bool maskZeroexpcomp)
{

    char buffer[256];

    if (maskZeroexpcomp) {
        if (expcomp != 0.0) {
            sprintf(buffer, "%0.2f", expcomp);
            return buffer;
        } else {
            return "";
        }
    } else {
        sprintf(buffer, "%0.2f", expcomp);
        return buffer;
    }
}

double FramesMetaData::shutterFromString(std::string s)
{

    size_t i = s.find_first_of('/');

    if (i == std::string::npos) {
        return atof(s.c_str());
    } else {
        return atof(s.substr(0, i).c_str()) / atof(s.substr(i + 1).c_str());
    }
}

double FramesMetaData::apertureFromString(std::string s)
{

    return atof(s.c_str());
}

extern "C" {

#include <libiptcdata/iptc-data.h>
#include <libiptcdata/iptc-jpeg.h>

    struct _IptcDataPrivate {
        unsigned int ref_count;

        IptcLog *log;
        IptcMem *mem;
    };

    IptcData *
    iptc_data_new_from_jpeg_file(FILE *infile)
    {
        IptcData *d;
        unsigned char * buf;
        int buf_len = 256 * 256;
        int len, offset;
        unsigned int iptc_len;

        if (!infile) {
            return nullptr;
        }

        d = iptc_data_new();

        if (!d) {
            return nullptr;
        }

        buf = (unsigned char*)iptc_mem_alloc(d->priv->mem, buf_len);

        if (!buf) {
            iptc_data_unref(d);
            return nullptr;
        }

        len = iptc_jpeg_read_ps3(infile, buf, buf_len);

        if (len <= 0) {
            goto failure;
        }

        offset = iptc_jpeg_ps3_find_iptc(buf, len, &iptc_len);

        if (offset <= 0) {
            goto failure;
        }

        iptc_data_load(d, buf + offset, iptc_len);

        iptc_mem_free(d->priv->mem, buf);
        return d;

failure:
        iptc_mem_free(d->priv->mem, buf);
        iptc_data_unref(d);
        return nullptr;
    }

}

FramesData::FramesData(const Glib::ustring& fname, std::unique_ptr<RawMetaDataLocation> rml, bool firstFrameOnly) :
    iptc(nullptr), dcrawFrameCount(0)
{
    if (rml && (rml->exifBase >= 0 || rml->ciffBase >= 0)) {
        FILE* f = g_fopen(fname.c_str(), "rb");

        if (f) {
            const bool has_rml_exif_base = rml->exifBase >= 0;
            rtexif::ExifManager exifManager(f, std::move(rml), firstFrameOnly);

            if (has_rml_exif_base) {
                if (exifManager.f && exifManager.rml) {
                    if (exifManager.rml->exifBase >= 0) {
                        exifManager.parseRaw();

                    } else if (exifManager.rml->ciffBase >= 0) {
                        exifManager.parseCIFF();
                    }
                }

                // copying roots
                roots = exifManager.roots;

                // creating FrameData
                for (auto currFrame : exifManager.frames) {
                    FrameData* fd = new FrameData(currFrame, currFrame->getRoot(), roots.at(0));

                    frames.push_back(fd);
                }

                for (auto currRoot : roots) {
                    rtexif::Tag* t = currRoot->getTag(0x83BB);

                    if (t && !iptc) {
                        iptc = iptc_data_new_from_data((unsigned char*)t->getValue(), (unsigned)t->getValueSize());
                        break;
                    }
                }
            }

            fclose(f);
        }
    } else if (hasJpegExtension(fname)) {
        FILE* f = g_fopen(fname.c_str(), "rb");

        if (f) {
            rtexif::ExifManager exifManager(f, std::move(rml), true);

            if (exifManager.f) {
                exifManager.parseJPEG();
                roots = exifManager.roots;

                for (auto currFrame : exifManager.frames) {
                    FrameData* fd = new FrameData(currFrame, currFrame->getRoot(), roots.at(0));
                    frames.push_back(fd);
                }

                rewind(exifManager.f);  // Not sure this is necessary
                iptc = iptc_data_new_from_jpeg_file(exifManager.f);
            }

            fclose(f);
        }
    } else if (hasTiffExtension(fname)) {
        FILE* f = g_fopen(fname.c_str(), "rb");

        if (f) {
            rtexif::ExifManager exifManager(f, std::move(rml), firstFrameOnly);

            exifManager.parseTIFF();
            roots = exifManager.roots;

            // creating FrameData
            for (auto currFrame : exifManager.frames) {
                FrameData* fd = new FrameData(currFrame, currFrame->getRoot(), roots.at(0));

                frames.push_back(fd);
            }

            for (auto currRoot : roots) {
                rtexif::Tag* t = currRoot->getTag(0x83BB);

                if (t && !iptc) {
                    iptc = iptc_data_new_from_data((unsigned char*)t->getValue(), (unsigned)t->getValueSize());
                    break;
                }
            }

            fclose(f);
        }
    }
}

FramesData::~FramesData()
{
    for (auto currRoot : roots) {
        delete currRoot;
    }

    if (iptc) {
        iptc_data_free(iptc);
    }
}
