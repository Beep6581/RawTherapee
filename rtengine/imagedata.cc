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

using namespace rtengine;

extern "C" IptcData *iptc_data_new_from_jpeg_file (FILE* infile);

namespace
{

Glib::ustring to_utf8 (const std::string& str)
{
    try {
        return Glib::locale_to_utf8 (str);
    } catch (Glib::Error&) {
        return Glib::convert_with_fallback (str, "UTF-8", "ISO-8859-1", "?");
    }
}

}

FramesMetaData* FramesMetaData::fromFile (const Glib::ustring& fname, RawMetaDataLocation* rml, bool firstFrameOnly)
{
    return new FramesData (fname, rml, firstFrameOnly);
}

FrameData::FrameData ()
    : root(nullptr), iptc(nullptr), time(), timeStamp(), iso_speed(0), aperture(0.), focal_len(0.), focal_len35mm(0.), focus_dist(0.f),
      shutter(0.), expcomp(0.), make("Unknown"), model("Unknown"), orientation("Unknown"), lens("Unknown"),
      sampleFormat(IIOSF_UNKNOWN), isPixelShift(false), isHDR(false)
{
    memset (&time, 0, sizeof(time));
}

RawFrameData::RawFrameData (rtexif::ExifManager &exifManager)
{
    bool rootCreated = false;
    if (exifManager.f && exifManager.rml) {
        if (exifManager.rml->exifBase >= 0) {
            root = exifManager.parse ();

            if (root) {
                rtexif::Tag* t = root->getTag (0x83BB);

                if (t) {
                    iptc = iptc_data_new_from_data ((unsigned char*)t->getValue (), (unsigned)t->getValueSize ());
                }
                extractInfo ();
                rootCreated = true;
            }
        } else if (exifManager.rml->ciffBase >= 0) {
            root = exifManager.parseCIFF ();
            extractInfo ();
            rootCreated = true;
        }
    }
    if (!rootCreated) {
        root = new rtexif::TagDirectory ();
    }
}

JpegFrameData::JpegFrameData (rtexif::ExifManager &exifManager)
{
    bool rootCreated = false;
    if (exifManager.f) {
        root = exifManager.parseJPEG ();
        if (root) {
            extractInfo ();
            rootCreated = true;
        }
        rewind (exifManager.f); // Not sure this is necessary
        iptc = iptc_data_new_from_jpeg_file (exifManager.f);
    }
    if (!rootCreated) {
        root = new rtexif::TagDirectory ();
    }
}

TiffFrameData::TiffFrameData (rtexif::ExifManager &exifManager)
{
    bool rootCreated = false;
    if (exifManager.f) {
        root = exifManager.parseTIFF ();
        extractInfo ();

        if (root) {
            rtexif::Tag* t = root->getTag (0x83BB);

            if (t) {
                iptc = iptc_data_new_from_data ((unsigned char*)t->getValue (), (unsigned)t->getValueSize ());
            }
            rootCreated = true;
        }
    }
    if (!rootCreated) {
        root = new rtexif::TagDirectory ();
    }
}

void FrameData::extractInfo ()
{

    if (!root) {
        return;
    }

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

    if (root->getTag("Make")) {
        make = root->getTag ("Make")->valueToString();
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

    if (root->getTag("Model")) {
        model = root->getTag("Model")->valueToString();
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

        if (model.find( "Digital Camera ") != std::string::npos) {
            model.erase(0, 15);
        }
    } else {
        model = "Unknown";
    }

    if (root->getTag ("Orientation")) {
        orientation = root->getTag ("Orientation")->valueToString ();
    }

    rtexif::Tag* mnoteTag = root->findTag("MakerNote");
    rtexif::TagDirectory* mnote = nullptr;
    if (mnoteTag) {
        mnote = mnoteTag->getDirectory();
    }

    rtexif::TagDirectory* exif = nullptr;
    if (root->getTag ("Exif")) {
        exif = root->getTag ("Exif")->getDirectory ();
    }

    if (exif) {

        // standard exif tags
        if (exif->getTag ("ShutterSpeedValue")) {
            shutter = exif->getTag ("ShutterSpeedValue")->toDouble ();
        }

        if (exif->getTag ("ExposureTime")) {
            shutter = exif->getTag ("ExposureTime")->toDouble ();
        }

        if (exif->getTag ("ApertureValue")) {
            aperture = exif->getTag ("ApertureValue")->toDouble ();
        }

        if (exif->getTag ("FNumber")) {
            aperture = exif->getTag ("FNumber")->toDouble ();
        }

        if (exif->getTag ("ExposureBiasValue")) {
            expcomp = exif->getTag ("ExposureBiasValue")->toDouble ();
        }

        if (exif->getTag ("FocalLength")) {
            focal_len = exif->getTag ("FocalLength")->toDouble ();
        }

        if (exif->getTag ("FocalLengthIn35mmFilm")) {
            focal_len35mm = exif->getTag ("FocalLengthIn35mmFilm")->toDouble ();
        }

        // Focus distance from EXIF or XMP. MakerNote ones are scattered and partly encrypted
        int num = -3, denom = -3;

        // First try, offical EXIF. Set by Adobe on some DNGs
        rtexif::Tag* pDst = exif->getTag("SubjectDistance");

        if (pDst) {
            int num, denom;
            pDst->toRational(num, denom);
        } else {
            // Second try, XMP data
            char sXMPVal[64];

            if (root->getXMPTagValue("aux:ApproximateFocusDistance", sXMPVal)) {
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

        if (exif->getTag ("ISOSpeedRatings")) {
            iso_speed = exif->getTag ("ISOSpeedRatings")->toDouble ();
        }

        if (exif->getTag ("DateTimeOriginal")) {
            if (sscanf ((const char*)exif->getTag("DateTimeOriginal")->getValue(), "%d:%d:%d %d:%d:%d", &time.tm_year, &time.tm_mon, &time.tm_mday, &time.tm_hour, &time.tm_min, &time.tm_sec) == 6) {
                time.tm_year -= 1900;
                time.tm_mon -= 1;
                time.tm_isdst = -1;
                timeStamp = mktime(&time);
            }
        }

        rtexif::Tag *snTag = exif->findTag ("SerialNumber");

        if(!snTag) {
            snTag = exif->findTag ("InternalSerialNumber");
        }

        if ( snTag ) {
            serial = snTag->valueToString();
        }

        // guess lens...
        lens = "Unknown";

        // Sometimes (e.g. DNG) EXIF already contains lens data

        if(!make.compare (0, 8, "FUJIFILM")) {
            if(exif->getTag ("LensModel")) {
                lens = exif->getTag ("LensModel")->valueToString ();
            }
        } else if(!make.compare (0, 4, "SONY")) {
            if (iso_speed == 65535 || iso_speed == 0) {
                rtexif::Tag* isoTag = exif->getTag ("RecommendedExposureIndex");

                if(isoTag) {
                    iso_speed = isoTag->toDouble();
                }
            }

        }

        if (lens == "Unknown") {

            if (mnote) {

                if (!make.compare (0, 5, "NIKON")) {
                    // ISO at max value supported, check manufacturer specific
                    if (iso_speed == 65535 || iso_speed == 0) {
                        rtexif::Tag* isoTag = mnote->getTagP("ISOInfo/ISO");

                        if (isoTag) {
                            iso_speed = isoTag->toInt();
                        }
                    }

                    bool lensOk = false;

                    if (mnote->getTag ("LensData")) {
                        std::string ldata = mnote->getTag ("LensData")->valueToString ();
                        size_t pos;

                        if (ldata.size() > 10 && (pos = ldata.find ("Lens = ")) != Glib::ustring::npos) {
                            lens = ldata.substr (pos + 7);

                            if (lens.compare (0, 7, "Unknown")) {
                                lensOk = true;
                            } else {
                                size_t pos = lens.find("$FL$");        // is there a placeholder for focallength?

                                if(pos != Glib::ustring::npos) {                // then fill in focallength
                                    lens = lens.replace(pos, 4, exif->getTag ("FocalLength")->valueToString ());

                                    if(mnote->getTag ("LensType")) {
                                        std::string ltype = mnote->getTag ("LensType")->valueToString ();

                                        if(ltype.find("MF = Yes") != Glib::ustring::npos) { // check, whether it's a MF lens, should be always
                                            lens = lens.replace(0, 7, "MF");
                                        }

                                        lensOk = true;
                                    }
                                }
                            }
                        }
                    }

                    if (!lensOk && mnote->getTag ("Lens")) {
                        std::string ldata = mnote->getTag ("Lens")->valueToString ();
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
                            n[m] = (double) nom / std::max(den,1);
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
                        if(mnote->getTag ("LensType")) {
                            std::string ltype = mnote->getTag ("LensType")->valueToString ();

                            if(ltype.find("MF = Yes") != Glib::ustring::npos) { // check, whether it's a MF lens
                                lens = lens.replace(0, 7, "MF");    // replace 'Unknwon' with 'MF'
                            } else {
                                lens = lens.replace(0, 7, "AF");    // replace 'Unknwon' with 'AF'
                            }
                        }
                    }
                } else if (!make.compare (0, 5, "Canon")) {
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

                    if ( lt ) {
                        std::string ldata = lt->valueToString ();

                        if (ldata.size() > 1) {
                            found = true;
                            lens = "Canon " + ldata;
                        }
                    }

                    if( !found || lens.substr(lens.find(' ')).length() < 7 ) {
                        lt = mnote->findTag("LensID");

                        if ( lt ) {
                            std::string ldata = lt->valueToString ();

                            if (ldata.size() > 1) {
                                lens = ldata;
                            }
                        }
                    }
                } else if (!make.compare (0, 6, "PENTAX") || (!make.compare (0, 5, "RICOH") && !model.compare (0, 6, "PENTAX"))) {
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
                    if (mnote->getTag ("LensType")) {
                        lens = mnote->getTag ("LensType")->valueToString ();
                    }

                    // Try to get the FocalLength from the LensInfo structure, where length below 10mm will be correctly set
                    rtexif::Tag* flt = mnote->getTagP ("LensInfo/FocalLength");

                    if (flt) {
                        // Don't replace Exif focal_len if Makernotes focal_len is 0
                        if (flt->toDouble() > 0) {
                            focal_len = flt->toDouble ();
                        }
                    } else if ((flt = mnote->getTagP ("FocalLength"))) {
                        rtexif::Tag* flt = mnote->getTag ("FocalLength");
                        focal_len = flt->toDouble ();
                    }

                    if (mnote->getTag ("FocalLengthIn35mmFilm")) {
                        focal_len35mm = mnote->getTag ("FocalLengthIn35mmFilm")->toDouble ();
                    }
                } else if (mnote && (!make.compare (0, 4, "SONY") || !make.compare (0, 6, "KONICA"))) {
                    if (mnote->getTag ("LensID")) {
                        lens = mnote->getTag ("LensID")->valueToString ();
                    }
                } else if (!make.compare (0, 7, "OLYMPUS")) {
                    if (mnote->getTag ("Equipment"))  {
                        rtexif::TagDirectory* eq = mnote->getTag ("Equipment")->getDirectory ();

                        if (eq->getTag ("LensType")) {
                            lens = eq->getTag ("LensType")->valueToString ();
                        }
                    }
                } else if (mnote && !make.compare (0, 9, "Panasonic")) {
                    if (mnote->getTag ("LensType")) {
                        std::string panalens = mnote->getTag("LensType")->valueToString();

                        if (panalens.find("LUMIX") != Glib::ustring::npos) {
                            lens = "Panasonic " + panalens;
                        }
                        else {
                            lens = panalens;
                        }
                    }
                }
            } else if (exif->getTag ("DNGLensInfo")) {
                lens = exif->getTag ("DNGLensInfo")->valueToString ();
            } else if (exif->getTag ("LensModel")) {
                lens = exif->getTag ("LensModel")->valueToString ();
            } else if (exif->getTag ("LensInfo")) {
                lens = exif->getTag ("LensInfo")->valueToString ();
            }
        }
    }


    // -----------------------  Special file type detection (HDR, PixelShift) ------------------------


    uint16 bitspersample = 0, samplesperpixel = 0, sampleformat = 0, photometric = 0, compression = 0;
    rtexif::Tag* bps = root->findTag("BitsPerSample");
    rtexif::Tag* spp = root->findTag("SamplesPerPixel");
    rtexif::Tag* sf = root->findTag("SampleFormat");
    rtexif::Tag* pi = root->findTag("PhotometricInterpretation");
    rtexif::Tag* c = root->findTag("Compression");

    if (mnote && (!make.compare (0, 6, "PENTAX") || (!make.compare (0, 5, "RICOH") && !model.compare (0, 6, "PENTAX")))) {
        rtexif::Tag* hdr = mnote->findTag("HDR");
        if (hdr) {
            if (hdr->toInt() > 0 && hdr->toInt(2) > 0) {
                isHDR = true;
            }
        } else {
            rtexif::Tag* dm = mnote->findTag("DriveMode");
            if (dm) {
                char buffer[60];
                dm->toString(buffer, 3);
                buffer[3] = 0;
                if (!strcmp(buffer, "HDR")) {
                    isHDR = true;
                }
            }
        }

        if (!isHDR) {
            rtexif::Tag* q = mnote->findTag("Quality");
            if (q && q->toInt() == 7) {
                isPixelShift = true;
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

    if ((!bps & !spp) || !pi) {
        return;
    }

    samplesperpixel = spp->toInt();
    bitspersample = bps->toInt();

    photometric = pi->toInt();
    if (photometric == PHOTOMETRIC_LOGLUV) {
        if (!c) {
            compression = COMPRESSION_NONE;
        } else {
            compression = c->toInt();
        }
    }

    if (photometric == PHOTOMETRIC_RGB || photometric == PHOTOMETRIC_MINISBLACK) {
        if ((samplesperpixel == 1 || samplesperpixel == 3 || samplesperpixel == 4) && sampleformat == SAMPLEFORMAT_UINT) {
            if (bitspersample == 8) {
                sampleFormat = IIOSF_UNSIGNED_CHAR;
            }

            if (bitspersample == 16) {
                sampleFormat = IIOSF_UNSIGNED_SHORT;
            }
        } else if (samplesperpixel == 3 && sampleformat == SAMPLEFORMAT_IEEEFP) {
            /*
             * Not yet supported
             *
             if (bitspersample==16) {
                sampleFormat = IIOSF_HALF;
            }*/
            if ((samplesperpixel == 3 || samplesperpixel == 4) && bitspersample == 32) {
                sampleFormat = IIOSF_FLOAT;
                isHDR = true;
            }
        }
    } else if (photometric == PHOTOMETRIC_CFA) {
        // Assuming Bayer or X-Trans raw file deliver 10, 12 14 or 16 bits uint, which is the case as of now
        sampleFormat = IIOSF_UNSIGNED_SHORT;
    } else if (samplesperpixel == 3 && photometric == PHOTOMETRIC_LOGLUV) {
        if (compression == COMPRESSION_SGILOG24) {
            sampleFormat = IIOSF_LOGLUV24;
            isHDR = true;
        } else if (compression == COMPRESSION_SGILOG) {
            sampleFormat = IIOSF_LOGLUV32;
            isHDR = true;
        }
    }
}

FrameData::~FrameData ()
{

    delete root;

    if (iptc) {
        iptc_data_free (iptc);
    }
}

procparams::IPTCPairs FrameData::getIPTCData () const
{

    procparams::IPTCPairs iptcc;

    if (!iptc) {
        return iptcc;
    }

    unsigned char buffer[2100];

    for (int i = 0; i < 16; i++) {
        IptcDataSet* ds = iptc_data_get_next_dataset (iptc, nullptr, IPTC_RECORD_APP_2, strTags[i].tag);

        if (ds) {
            iptc_dataset_get_data (ds, buffer, 2100);
            std::vector<Glib::ustring> icValues;
            icValues.push_back (to_utf8((char*)buffer));

            iptcc[strTags[i].field] = icValues;
            iptc_dataset_unref (ds);
        }
    }

    IptcDataSet* ds = nullptr;
    std::vector<Glib::ustring> keywords;

    while ((ds = iptc_data_get_next_dataset (iptc, ds, IPTC_RECORD_APP_2, IPTC_TAG_KEYWORDS))) {
        iptc_dataset_get_data (ds, buffer, 2100);
        keywords.push_back (to_utf8((char*)buffer));
    }

    iptcc["Keywords"] = keywords;
    ds = nullptr;
    std::vector<Glib::ustring> suppCategories;

    while ((ds = iptc_data_get_next_dataset (iptc, ds, IPTC_RECORD_APP_2, IPTC_TAG_SUPPL_CATEGORY))) {
        iptc_dataset_get_data (ds, buffer, 2100);
        suppCategories.push_back (to_utf8((char*)buffer));
        iptc_dataset_unref (ds);
    }

    iptcc["SupplementalCategories"] = suppCategories;
    return iptcc;
}


bool FrameData::getPixelShift () const
{
    return isPixelShift;
}
bool FrameData::getHDR () const
{
    return isHDR;
}
IIOSampleFormat FrameData::getSampleFormat () const
{
    return sampleFormat;
}
rtexif::TagDirectory* FrameData::getExifData () const
{
    return root;
}
bool FrameData::hasExif () const
{
    return root && root->getCount();
}
bool FrameData::hasIPTC () const
{
    return iptc;
}
tm FrameData::getDateTime () const
{
    return time;
}
time_t FrameData::getDateTimeAsTS () const
{
    return timeStamp;
}
int FrameData::getISOSpeed () const
{
    return iso_speed;
}
double FrameData::getFNumber () const
{
    return aperture;
}
double FrameData::getFocalLen () const
{
    return focal_len;
}
double FrameData::getFocalLen35mm () const
{
    return focal_len35mm;
}
float FrameData::getFocusDist () const
{
    return focus_dist;
}
double FrameData::getShutterSpeed () const
{
    return shutter;
}
double FrameData::getExpComp () const
{
    return expcomp;
}
std::string FrameData::getMake () const
{
    return make;
}
std::string FrameData::getModel () const
{
    return model;
}
std::string FrameData::getLens () const
{
    return lens;
}
std::string FrameData::getSerialNumber () const
{
    return serial;
}
std::string FrameData::getOrientation () const
{
    return orientation;
}



void FramesData::setDCRawFrameCount (unsigned int frameCount)
{
    dcrawFrameCount = frameCount;
}

unsigned int FramesData::getFrameCount () const
{
    return dcrawFrameCount ? dcrawFrameCount : frames.size();
}
FrameData *FramesData::getFrameData (int frame) const
{
    return frames.at(frame);
}

bool FramesData::getPixelShift (unsigned int frame) const
{
    // So far only Pentax provide multi-frame HDR file.
    // Only the first frame contains the HDR tag
    // If more brand have to be supported, this rule may need
    // to evolve

    //return frames.at(frame)->getPixelShift ();
    return frames.at(0)->getPixelShift ();
}
bool FramesData::getHDR (unsigned int frame) const
{
    // So far only Pentax provide multi-frame HDR file.
    // Only the first frame contains the HDR tag
    // If more brand have to be supported, this rule may need
    // to evolve

    //return frames.at(frame)->getHDR ();
    if (frames.size()) {
        return frames.at(0)->getHDR ();
    } else {
        return 0;
    }
}

IIOSampleFormat FramesData::getSampleFormat (unsigned int frame) const
{
    return frames.at(frame)->getSampleFormat ();
}

rtexif::TagDirectory* FramesData::getExifData (unsigned int frame) const
{
    return frames.at(frame)->getExifData ();
}
procparams::IPTCPairs FramesData::getIPTCData (unsigned int frame) const
{
    return frames.at(frame)->getIPTCData ();
}

bool FramesData::hasExif (unsigned int frame) const
{
    return frames.at(frame)->hasExif ();
}
bool FramesData::hasIPTC (unsigned int frame) const
{
    return frames.at(frame)->hasIPTC ();
}

tm FramesData::getDateTime (unsigned int frame) const
{
    return frames.at(frame)->getDateTime ();
}
time_t FramesData::getDateTimeAsTS(unsigned int frame) const
{
    return frames.at(frame)->getDateTimeAsTS ();
}
int FramesData::getISOSpeed (unsigned int frame) const
{
    return frames.at(frame)->getISOSpeed ();
}
double FramesData::getFNumber  (unsigned int frame) const
{
    return frames.at(frame)->getFNumber ();
}
double FramesData::getFocalLen (unsigned int frame) const
{
    return frames.at(frame)->getFocalLen ();
}
double FramesData::getFocalLen35mm (unsigned int frame) const
{
    return frames.at(frame)->getFocalLen35mm ();
}
float FramesData::getFocusDist (unsigned int frame) const
{
    return frames.at(frame)->getFocusDist ();
}
double FramesData::getShutterSpeed (unsigned int frame) const
{
    return frames.at(frame)->getShutterSpeed ();
}
double FramesData::getExpComp  (unsigned int frame) const
{
    return frames.at(frame)->getExpComp ();
}
std::string FramesData::getMake     (unsigned int frame) const
{
    return frames.at(frame)->getMake ();
}
std::string FramesData::getModel    (unsigned int frame) const
{
    return frames.at(frame)->getModel ();
}
std::string FramesData::getLens     (unsigned int frame) const
{
    return frames.at(frame)->getLens ();
}
std::string FramesData::getSerialNumber (unsigned int frame) const
{
    return frames.at(frame)->getSerialNumber ();
}
std::string FramesData::getOrientation (unsigned int frame) const
{
    return frames.at(frame)->getOrientation ();
}


//------inherited functions--------------//


std::string FramesMetaData::apertureToString (double aperture)
{

    char buffer[256];
    sprintf (buffer, "%0.1f", aperture);
    return buffer;
}

std::string FramesMetaData::shutterToString (double shutter)
{

    char buffer[256];

    if (shutter > 0.0 && shutter <= 0.5) {
        sprintf (buffer, "1/%0.0f", 1.0 / shutter);
    } else {
        sprintf (buffer, "%0.1f", shutter);
    }

    return buffer;
}

std::string FramesMetaData::expcompToString (double expcomp, bool maskZeroexpcomp)
{

    char buffer[256];

    if (maskZeroexpcomp) {
        if (expcomp != 0.0) {
            sprintf (buffer, "%0.2f", expcomp);
            return buffer;
        } else {
            return "";
        }
    } else {
        sprintf (buffer, "%0.2f", expcomp);
        return buffer;
    }
}

double FramesMetaData::shutterFromString (std::string s)
{

    size_t i = s.find_first_of ('/');

    if (i == std::string::npos) {
        return atof (s.c_str());
    } else {
        return atof (s.substr(0, i).c_str()) / atof (s.substr(i + 1).c_str());
    }
}

double FramesMetaData::apertureFromString (std::string s)
{

    return atof (s.c_str());
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
    iptc_data_new_from_jpeg_file (FILE *infile)
    {
        IptcData *d;
        unsigned char * buf;
        int buf_len = 256 * 256;
        int len, offset;
        unsigned int iptc_len;

        if (!infile) {
            return nullptr;
        }

        d = iptc_data_new ();

        if (!d) {
            return nullptr;
        }

        buf = (unsigned char*)iptc_mem_alloc (d->priv->mem, buf_len);

        if (!buf) {
            iptc_data_unref (d);
            return nullptr;
        }

        len = iptc_jpeg_read_ps3 (infile, buf, buf_len);

        if (len <= 0) {
            goto failure;
        }

        offset = iptc_jpeg_ps3_find_iptc (buf, len, &iptc_len);

        if (offset <= 0) {
            goto failure;
        }

        iptc_data_load (d, buf + offset, iptc_len);

        iptc_mem_free (d->priv->mem, buf);
        return d;

failure:
        iptc_mem_free (d->priv->mem, buf);
        iptc_data_unref (d);
        return nullptr;
    }

}

FramesData::FramesData (Glib::ustring fname, RawMetaDataLocation* rml, bool firstFrameOnly, bool loadAll) : dcrawFrameCount (0)
{
    if (rml && (rml->exifBase >= 0 || rml->ciffBase >= 0)) {
        FILE* f = g_fopen (fname.c_str (), "rb");

        if (f) {
            rtexif::ExifManager exifManager (f, rml, firstFrameOnly);

            if (rml->exifBase >= 0) {
                FrameData *idata = new RawFrameData (exifManager);
                frames.push_back(idata);
                if (rml && !firstFrameOnly) {
                    while (exifManager.getNextIFDOffset ()) {
                        int nextIFD = exifManager.getNextIFDOffset ();
                        exifManager.setIFDOffset (nextIFD);
                        idata = new RawFrameData (exifManager);
                        frames.push_back(idata);
                    }
                }
            }
            fclose (f);
        }
    } else if (hasJpegExtension(fname)) {
        FILE* f = g_fopen (fname.c_str (), "rb");

        if (f) {
            rtexif::ExifManager exifManager (f, rml, true);
            FrameData *idata = new JpegFrameData (exifManager);
            frames.push_back(idata);
            fclose (f);
        }
    } else if (hasTiffExtension(fname)) {
        FILE* f = g_fopen (fname.c_str (), "rb");

        if (f) {
            rtexif::ExifManager exifManager (f, rml, firstFrameOnly);
            FrameData *idata = new TiffFrameData (exifManager);
            frames.push_back(idata);
            if (rml && !firstFrameOnly) {
                while (exifManager.getNextIFDOffset ()) {
                    exifManager.setIFDOffset (exifManager.getNextIFDOffset ());
                    idata = new TiffFrameData (exifManager);
                    frames.push_back(idata);
                }
            }
            fclose (f);
        }
    }
}

FramesData::~FramesData ()
{
    for (auto currFrame : frames) {
        delete currFrame;
    }
}
