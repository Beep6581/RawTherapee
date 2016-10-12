/*
 *  This file is part of RawTherapee.
 */
#ifndef _KODAKATTRIBS_
#define _KODAKATTRIBS_

#include <string.h>
#include "rtexif.h"

namespace rtexif
{


void parseKodakIfdTextualInfo(Tag *textualInfo, Tag* exif_)
{
    // parse TextualInfo and copy values into corresponding standard Exif
    if (textualInfo->getType() != ASCII) {
        return;
    }

    TagDirectory *exif = exif_->getDirectory();
    char *value = (char *)textualInfo->getValue();
    int valuesize = textualInfo->getValueSize();

    char *p = value;
    char *pc, *plf;

    while ((pc = strchr(p, ':')) != nullptr && (plf = strchr(pc, '\n')) != nullptr) {
        while (*p == ' ') {
            p++;
        }

        size_t len = pc - p;

        while (len > 1 && p[len - 1] == ' ') {
            len--;
        }

        std::string key = std::string(p, len);
        ++pc;

        while (*pc == ' ') {
            pc++;
        }

        len = plf - pc;

        while (len > 1 && pc[len - 1] == ' ') {
            len--;
        }

        std::string val = std::string(pc, len);
        p = ++plf;

        // we pick out a few select tags here
        Tag *t;

        if (key == "Lens") {
            // Proback645 may have "Lens" but not "Focal Length"
            float flen = atof(val.c_str());

            if (flen != 0.0) {
                t = new Tag(exif, lookupAttrib(exifAttribs, "FocalLength"));
                t->initRational(flen * 32, 32);
                exif->replaceTag(t);
            }
        } else if (key == "Focal Length") {
            float flen = atof(val.c_str());

            if (flen != 0.0) {
                t = new Tag(exif, lookupAttrib(exifAttribs, "FocalLength"));
                t->initRational(flen * 32, 32);
                exif->replaceTag(t);
            }
        } else if (key == "Aperture") {
            float aperture = atof(&val.c_str()[1]);

            if (aperture != 0.0) {
                t = new Tag(exif, lookupAttrib(exifAttribs, "FNumber"));
                t->initRational((int)(aperture * 10), 10);
                exif->replaceTag(t);
            }
        } else if (key == "Exposure Bias" || key == "Compensation") {
            float bias = 0.0;

            if (val != "Off") {
                bias = atof(val.c_str());
            }

            t = new Tag (exif, lookupAttrib(exifAttribs, "ExposureBiasValue"));
            t->initRational ((int)(bias * 1000), 1000);
            exif->replaceTag(t);
        } else if (key == "ISO Speed") {
            t = new Tag (exif, lookupAttrib(exifAttribs, "ISOSpeedRatings"));
            t->initInt(atoi(val.c_str()), SHORT);
            exif->replaceTag(t);
        } else if (key == "Shutter") {
            const char *p1 = strchr(val.c_str(), '/');
            int a, b;

            if (p1 == nullptr) {
                a = atoi(val.c_str());
                b = 1;
            } else {
                a = atoi(val.c_str());
                b = atoi(&p1[1]);
            }

            t = new Tag (exif, lookupAttrib(exifAttribs, "ExposureTime"));
            t->initRational(a, b);
            exif->replaceTag(t);

            float ssv = -log2((float)a / (float)b); // convert to APEX value
            t = new Tag (exif, lookupAttrib(exifAttribs, "ShutterSpeedValue"));
            t->initRational(1000000 * ssv, 1000000);
            exif->replaceTag(t);
        } else if (key == "Flash Fired") {
            t = new Tag (exif, lookupAttrib(exifAttribs, "Flash"));

            if (val == "No") {
                t->initInt(0, SHORT);
            } else {
                // not sure if "Flash Fired" is only yes/no, only seen "No" in test pictures
                t->initInt(1, SHORT);
            }

            exif->replaceTag(t);
        } else if (key == "White balance") { // yes should be small 'b' int 'balance'.
            t = new Tag (exif, lookupAttrib(exifAttribs, "Flash"));
            t->initInt((val == "Auto") ? 0 : 1, SHORT);
            exif->replaceTag(t);
        }
    }
}

// table not complete, not all proprietary Kodak tags are known
const TagAttrib kodakIfdAttribs[] = {
    {0, AC_WRITE, 0, nullptr, 0x0001, AUTO, "UnknownEV?", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0003, AUTO, "ExposureValue", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x03e9, AUTO, "OriginalFileName", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x03eb, AUTO, "SensorLeftBorder", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x03ec, AUTO, "SensorTopBorder", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x03ed, AUTO, "SensorImageWidth", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x03ee, AUTO, "SensorImageHeight", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x03f1, AUTO, "TextualInfo", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x03fc, AUTO, "WhiteBalance", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x03fd, AUTO, "Processing", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0401, AUTO, "Time", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0414, AUTO, "NCDFileInfo", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0846, AUTO, "ColorTemperature", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0852, AUTO, "WB_RGBMul0", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0853, AUTO, "WB_RGBMul1", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0854, AUTO, "WB_RGBMul2", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0855, AUTO, "WB_RGBMul3", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x085c, AUTO, "WB_RGBCoeffs0", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x085d, AUTO, "WB_RGBCoeffs1", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x085e, AUTO, "WB_RGBCoeffs2", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x085f, AUTO, "WB_RGBCoeffs3", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0ce5, AUTO, "FirmwareVersion", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1391, AUTO, "ToneCurveFileName", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1784, AUTO, "ISO", &stdInterpreter},
    { -1, AC_DONTWRITE, 0,  nullptr, 0, AUTO, "", nullptr }
};

}
#endif

