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
#ifndef _STDATTRIBS_
#define _STDATTRIBS_

#include <rtexif.h>
#include <string>
#include <map>
#include <math.h>
#include <string.h>

namespace rtexif {

class ColorSpaceInterpreter : public ChoiceInterpreter {

    public:
        ColorSpaceInterpreter () {
            choices[1]      = "sRGB";
            choices[0xffff] = "Uncalibrated";
        }
};
ColorSpaceInterpreter colorSpaceInterpreter;

class ExposureProgramInterpreter : public ChoiceInterpreter {

    public:
        ExposureProgramInterpreter () {
            choices[0] = "Not defined";
            choices[1] = "Manual";
            choices[2] = "Normal program";
            choices[3] = "Aperture priority";
            choices[4] = "Shutter priority";
            choices[5] = "Creative program";
            choices[6] = "Action program";
            choices[7] = "Portrait mode";
            choices[8] = "Landscape mode";
        }
};
ExposureProgramInterpreter exposureProgramInterpreter;

class MeteringModeInterpreter : public ChoiceInterpreter {

    public:
        MeteringModeInterpreter () {
            choices[0] = "Unknown";
            choices[1] = "Average";
            choices[2] = "Center weighted";
            choices[3] = "Spot";
            choices[4] = "Multispot";
            choices[5] = "Pattern";
            choices[6] = "Partial";
            choices[255] = "Other";
        }
};
MeteringModeInterpreter meteringModeInterpreter;

class ExposureModeInterpreter : public ChoiceInterpreter {

    public:
        ExposureModeInterpreter () {
            choices[0] = "Auto exposure";
            choices[1] = "Manual exposure";
            choices[2] = "Auto bracket";
        }
};
ExposureModeInterpreter exposureModeInterpreter;

class WhiteBalanceInterpreter : public ChoiceInterpreter {

    public:
        WhiteBalanceInterpreter () {
            choices[0] = "Auto white balance";
            choices[1] = "Manual white balance";
        }
};
WhiteBalanceInterpreter whiteBalanceInterpreter;

class SceneCaptureInterpreter : public ChoiceInterpreter {

    public:
        SceneCaptureInterpreter () {
            choices[0] = "Standard";
            choices[1] = "Landscape";
            choices[2] = "Portrait";
            choices[3] = "Night scene";
        }
};
SceneCaptureInterpreter sceneCaptureInterpreter;

class GainControlInterpreter : public ChoiceInterpreter {

    public:
        GainControlInterpreter () {
            choices[0] = "None";
            choices[1] = "Low gain up";
            choices[2] = "High gain up";
            choices[3] = "Low gain down";
            choices[4] = "High gain down";
        }
};
GainControlInterpreter gainControlInterpreter;

class ContrastInterpreter : public ChoiceInterpreter {

    public:
        ContrastInterpreter () {
            choices[0] = "Normal";
            choices[1] = "Soft";
            choices[2] = "Hard";
        }
};
ContrastInterpreter contrastInterpreter;

class SharpnessInterpreter : public ChoiceInterpreter {

    public:
        SharpnessInterpreter () {
            choices[0] = "Normal";
            choices[1] = "Soft";
            choices[2] = "Hard";
        }
};
SharpnessInterpreter sharpnessInterpreter;

class SaturationInterpreter : public ChoiceInterpreter {

    public:
        SaturationInterpreter () {
            choices[0] = "Normal";
            choices[1] = "Low saturation";
            choices[2] = "High saturation";
        }
};
SaturationInterpreter saturationInterpreter;

class FlashInterpreter : public ChoiceInterpreter {

    public:
        FlashInterpreter () {
            choices[0x0000] = "Flash did not fire";
            choices[0x0001] = "Flash fired";
            choices[0x0005] = "Strobe return light not detected";
            choices[0x0007] = "Strobe return light detected";
            choices[0x0009] = "Flash fired, compulsory flash mode";
            choices[0x000D] = "Flash fired, compulsory flash mode, return light not detected";
            choices[0x000F] = "Flash fired, compulsory flash mode, return light detected";
            choices[0x0010] = "Flash did not fire, compulsory flash mode";
            choices[0x0018] = "Flash did not fire, auto mode";
            choices[0x0019] = "Flash fired, auto mode";
            choices[0x001D] = "Flash fired, auto mode, return light not detected";
            choices[0x001F] = "Flash fired, auto mode, return light detected";
            choices[0x0020] = "No flash function";
            choices[0x0041] = "Flash fired, red-eye reduction mode";
            choices[0x0045] = "Flash fired, red-eye reduction mode, return light not detected";
            choices[0x0047] = "Flash fired, red-eye reduction mode, return light detected";
            choices[0x0049] = "Flash fired, compulsory flash mode, red-eye reduction mode";
            choices[0x004D] = "Flash fired, compulsory flash mode, red-eye reduction mode, return light not detected";
            choices[0x004F] = "Flash fired, compulsory flash mode, red-eye reduction mode, return light detected";
            choices[0x0059] = "Flash fired, auto mode, red-eye reduction mode";
            choices[0x005D] = "Flash fired, auto mode, return light not detected, red-eye reduction mode";
            choices[0x005F] = "Flash fired, auto mode, return light detected, red-eye reduction mode";
        }
};
FlashInterpreter flashInterpreter;

class LightSourceInterpreter : public ChoiceInterpreter {

    public:
        LightSourceInterpreter () {
            choices[0] = "Unknown";
            choices[1] = "Daylight";
            choices[2] = "Fluorescent";
            choices[3] = "Tungsten";
            choices[4] = "Flash";
            choices[9] = "Fine weather";
            choices[10] = "Cloudy weather";
            choices[11] = "Shade";
            choices[12] = "Daylight fluorescent";
            choices[13] = "Day white fluorescent";
            choices[14] = "Cool white fluorescent";
            choices[15] = "White fluorescent";
            choices[17] = "Standard light A";
            choices[18] = "Standard light B";
            choices[19] = "Standard light C";
            choices[20] = "D55";
            choices[21] = "D65";
            choices[22] = "D75";
            choices[23] = "D50";
            choices[24] = "ISO studio tungsten";
            choices[255] = "Other light source";
        }
};
LightSourceInterpreter lightSourceInterpreter;

class CompressionInterpreter : public ChoiceInterpreter {

    public:
        CompressionInterpreter () {
            choices[1] = "Uncompressed";
            choices[6] = "JPEG Compression";
        }
};
CompressionInterpreter compressionInterpreter;

class PhotometricInterpreter : public ChoiceInterpreter {

    public:
        PhotometricInterpreter () {
            choices[2] = "RGB";
            choices[6] = "YCbCr";
        }
};
PhotometricInterpreter photometricInterpreter;

class PlanarConfigInterpreter : public ChoiceInterpreter {

    public:
        PlanarConfigInterpreter () {
            choices[1] = "Chunky format";
            choices[2] = "Planar format";
        }
};
PlanarConfigInterpreter planarConfigInterpreter;

class FNumberInterpreter : public Interpreter {
    public:
        FNumberInterpreter () {}
        virtual std::string toString (Tag* t) {
            sprintf (buffer, "%0.1f", t->toDouble());
            return buffer;
        }
};
FNumberInterpreter fNumberInterpreter;

class ApertureInterpreter : public Interpreter {
    public:
        ApertureInterpreter () {}
        virtual std::string toString (Tag* t) {
            sprintf (buffer, "%0.1f", pow(2.0, t->toDouble()/2.0));
            return buffer;
        }
};
ApertureInterpreter apertureInterpreter;

class ExposureBiasInterpreter : public Interpreter {
    public:
        ExposureBiasInterpreter () {}
        virtual std::string toString (Tag* t) {
            sprintf (buffer, "%+0.2f", t->toDouble());
            return buffer;
        }
};
ExposureBiasInterpreter exposureBiasInterpreter;

class ShutterSpeedInterpreter : public Interpreter {
    public:
        ShutterSpeedInterpreter () {}
        virtual std::string toString (Tag* t) {
            double d = pow (2.0, -t->toDouble());
            if (d > 0.0 && d < 0.9)
                sprintf (buffer, "1/%0.0f", 1.0 / d);
            else
                sprintf (buffer, "%0.1f", d);
            return buffer;
        }
};
ShutterSpeedInterpreter shutterSpeedInterpreter;

class ExposureTimeInterpreter : public Interpreter {
    public:
        ExposureTimeInterpreter () {}
        virtual std::string toString (Tag* t) {
            double d = t->toDouble();
            if (d > 0.0 && d < 0.9)
                sprintf (buffer, "1/%0.0f", 1.0 / d);
            else
                sprintf (buffer, "%0.1f", d);
            return buffer;
        }
};
ExposureTimeInterpreter exposureTimeInterpreter;

class FocalLengthInterpreter : public Interpreter {
    public:
        FocalLengthInterpreter () {}
        virtual std::string toString (Tag* t) {
            sprintf (buffer, "%0.1f", t->toDouble());
            return buffer;
        }
};
FocalLengthInterpreter focalLengthInterpreter;

class UserCommentInterpreter : public Interpreter {
    public:
        UserCommentInterpreter () {}
        virtual std::string toString (Tag* t) {
            if (!strncmp((char*)t->getValue(), "ASCII\0\0\0",8))
                strncpy (buffer, (char*)t->getValue()+8, t->getCount()-8);
            else
                buffer[0]=0;
            return buffer;
        }
        virtual void fromString (Tag* t, const std::string& value) {
            memcpy (buffer, "ASCII\0\0\0", 8);
            strcpy (buffer+8, value.c_str());
            t->fromString (buffer, value.size() + 9);
        }
};
UserCommentInterpreter userCommentInterpreter;

const TagAttrib exifAttribs[] = {
 0, 2, 0, 0, 0x0103, "Compression", &compressionInterpreter,
 0, 2, 0, 0, 0xA000, "FlashpixVersion", &stdInterpreter,
 0, 2, 0, 0, 0xA001, "ColorSpace", &colorSpaceInterpreter,
 0, 1, 0, 0, 0x9000, "ExifVersion", &stdInterpreter,
 0, 1, 0, 0, 0x9003, "DateTimeOriginal", &stdInterpreter,
 0, 1, 0, 0, 0x9004, "DateTimeDigitized", &stdInterpreter,
 0, 2, 0, 0, 0x9101, "ComponentsConfiguration", &stdInterpreter,
 0, 2, 0, 0, 0x9102, "CompressedBitsPerPixel", &stdInterpreter,
 0, 2, 0, 0, 0xA002, "PixelXDimension", &stdInterpreter,
 0, 2, 0, 0, 0xA003, "PixelYDimension", &stdInterpreter,
 0, 1, 0, 0, 0x927C, "MakerNote", &stdInterpreter,
 0, 1, 1, 0, 0x9286, "UserComment", &userCommentInterpreter,
 1, 0, 0, 0, 0xA004, "RelatedSoundFile", &stdInterpreter,
 0, 1, 0, 0, 0x9290, "SubSecTime", &stdInterpreter,
 0, 1, 0, 0, 0x9291, "SubSecTimeOriginal", &stdInterpreter,
 0, 1, 0, 0, 0x9292, "SubSecTimeDigitized", &stdInterpreter,
 0, 1, 0, 0, 0xA420, "ImageUniqueID", &stdInterpreter,
 0, 1, 0, 0, 0x829A, "ExposureTime", &exposureTimeInterpreter,
 0, 1, 0, 0, 0x829D, "FNumber", &fNumberInterpreter,
 0, 1, 0, 0, 0x8822, "ExposureProgram", &exposureProgramInterpreter,
 0, 1, 0, 0, 0x8824, "SpectralSensitivity", &stdInterpreter,
 0, 1, 0, 0, 0x8827, "ISOSpeedRatings", &stdInterpreter,
 0, 1, 0, 0, 0x8828, "OECF", &stdInterpreter,
 0, 1, 0, 0, 0x9201, "ShutterSpeedValue", &shutterSpeedInterpreter,
 0, 1, 0, 0, 0x9202, "ApertureValue", &apertureInterpreter,
 0, 1, 0, 0, 0x9203, "BrightnessValue", &stdInterpreter,
 0, 1, 0, 0, 0x9204, "ExposureBiasValue", &exposureBiasInterpreter,
 0, 1, 0, 0, 0x9205, "MaxApertureValue", &apertureInterpreter,
 0, 1, 0, 0, 0x9206, "SubjectDistance", &stdInterpreter,
 0, 1, 0, 0, 0x9207, "MeteringMode", &meteringModeInterpreter,
 0, 1, 0, 0, 0x9208, "LightSource", &lightSourceInterpreter,
 0, 1, 0, 0, 0x9209, "Flash", &flashInterpreter,
 0, 1, 0, 0, 0x920A, "FocalLength", &focalLengthInterpreter,
 0, 1, 0, 0, 0x9214, "SubjectArea", &stdInterpreter,
 0, 0, 0, 0, 0x9216, "TIFFEPSStandardID", &stdInterpreter,
 0, 1, 0, 0, 0x9217, "SensingMethod", &stdInterpreter,
 0, 1, 0, 0, 0xA20B, "FlashEnergy", &stdInterpreter,
 0, 1, 0, 0, 0xA20C, "SpatialFrequencyResponse", &stdInterpreter,
 0, 1, 0, 0, 0xA20E, "FocalPlaneXResolution", &stdInterpreter,
 0, 1, 0, 0, 0xA20F, "FocalPlaneYResolution", &stdInterpreter,
 0, 1, 0, 0, 0xA210, "FocalPlaneResolutionUnit", &stdInterpreter,
 0, 1, 0, 0, 0xA214, "SubjectLocation", &stdInterpreter,
 0, 1, 0, 0, 0xA215, "ExposureIndex", &stdInterpreter,
 0, 1, 0, 0, 0xA217, "SensingMethod", &stdInterpreter,
 0, 1, 0, 0, 0xA300, "FileSource", &stdInterpreter,
 0, 1, 0, 0, 0xA301, "SceneType", &stdInterpreter,
 0, 0, 0, 0, 0xA302, "CFAPattern", &stdInterpreter,
 0, 1, 0, 0, 0xA401, "CustomRendered", &stdInterpreter,
 0, 1, 0, 0, 0xA402, "ExposureMode", &exposureModeInterpreter,
 0, 1, 0, 0, 0xA403, "WhiteBalance", &whiteBalanceInterpreter,
 0, 1, 0, 0, 0xA404, "DigitalZoomRatio", &stdInterpreter,
 0, 1, 0, 0, 0xA405, "FocalLengthIn35mmFilm", &stdInterpreter,
 0, 1, 0, 0, 0xA406, "SceneCaptureType", &sceneCaptureInterpreter,
 0, 1, 0, 0, 0xA407, "GainControl", &gainControlInterpreter,
 0, 1, 0, 0, 0xA408, "Contrast", &contrastInterpreter,
 0, 1, 0, 0, 0xA409, "Saturation", &saturationInterpreter,
 0, 1, 0, 0, 0xA40A, "Sharpness", &sharpnessInterpreter,
 0, 1, 0, 0, 0xA40B, "DeviceSettingDescription", &stdInterpreter,
 0, 1, 0, 0, 0xA40C, "SubjectDistanceRange", &stdInterpreter,
 0, 0, 0, 0, 0x828d, "CFAPattern", &stdInterpreter,
 0, 0, 0, 0, 0x828e, "CFARepeatPatternDim", &stdInterpreter,
-1, 0, 0, 0, 0, "", NULL };
// 0, 0xA005, LONG,       1, "Interoperability tag",                      "Interoperability IFD Pointer"};

const TagAttrib gpsAttribs[] = {
 0, 1, 0, 0, 0x0000, "GPSVersionID", &stdInterpreter,
 0, 1, 0, 0, 0x0001, "GPSLatitudeRef", &stdInterpreter,
 0, 1, 0, 0, 0x0002, "GPSLatitude", &stdInterpreter,
 0, 1, 0, 0, 0x0003, "GPSLongitudeRef", &stdInterpreter,
 0, 1, 0, 0, 0x0004, "GPSLongitude", &stdInterpreter,
 0, 1, 0, 0, 0x0005, "GPSAltitudeRef", &stdInterpreter,
 0, 1, 0, 0, 0x0006, "GPSAltitude", &stdInterpreter,
 0, 1, 0, 0, 0x0007, "GPSTimeStamp", &stdInterpreter,
 0, 1, 0, 0, 0x0008, "GPSSatelites", &stdInterpreter,
 0, 1, 0, 0, 0x0009, "GPSStatus", &stdInterpreter,
 0, 1, 0, 0, 0x000a, "GPSMeasureMode", &stdInterpreter,
 0, 1, 0, 0, 0x000b, "GPSDOP", &stdInterpreter,
 0, 1, 0, 0, 0x000c, "GPSSpeedRef", &stdInterpreter,
 0, 1, 0, 0, 0x000d, "GPSSpeed", &stdInterpreter,
 0, 1, 0, 0, 0x000e, "GPSTrackRef", &stdInterpreter,
 0, 1, 0, 0, 0x000f, "GPSTrack", &stdInterpreter,
 0, 1, 0, 0, 0x0010, "GPSImgDirectionRef", &stdInterpreter,
 0, 1, 0, 0, 0x0011, "GPSImgDirection", &stdInterpreter,
 0, 1, 0, 0, 0x0012, "GPSMapDatum", &stdInterpreter,
 0, 1, 0, 0, 0x0013, "GPSDestLatitudeRef", &stdInterpreter,
 0, 1, 0, 0, 0x0014, "GPSDestLatitude", &stdInterpreter,
 0, 1, 0, 0, 0x0015, "GPSDestLongitudeRef", &stdInterpreter,
 0, 1, 0, 0, 0x0016, "GPSDestLongitude", &stdInterpreter,
 0, 1, 0, 0, 0x0017, "GPSDestBearingRef", &stdInterpreter,
 0, 1, 0, 0, 0x0018, "GPSDestBearing", &stdInterpreter,
 0, 1, 0, 0, 0x0019, "GPSDestDistanceRef", &stdInterpreter,
 0, 1, 0, 0, 0x001a, "GPSDestDistance", &stdInterpreter,
 0, 1, 0, 0, 0x001b, "GPSProcessingMethod", &stdInterpreter,
 0, 1, 0, 0, 0x001c, "GPSAreaInformation", &stdInterpreter,
 0, 1, 0, 0, 0x001d, "GPSDateStamp", &stdInterpreter,
 0, 1, 0, 0, 0x001e, "GPSDifferential", &stdInterpreter,
 -1, 0, 0,  0, 0, "", NULL };



const TagAttrib iopAttribs[] = {
 0, 1, 0, 0, 0x0001, "InteroperabilityIndex", &stdInterpreter,
 0, 1, 0, 0, 0x0002, "InteroperabilityVersion", &stdInterpreter,
-1, 0, 0, 0, 0, "", NULL };

 const TagAttrib ifdAttribs[] = { 
 0, 2, 0, 0, 0x0017, "PanaISO", &stdInterpreter,
 0, 2, 0, 0, 0x0100, "ImageWidth", &stdInterpreter,
 0, 2, 0, 0, 0x0101, "ImageHeight", &stdInterpreter,
 0, 2, 0, 0, 0x0102, "BitsPerSample", &stdInterpreter,
 0, 2, 0, 0, 0x0103, "Compression", &compressionInterpreter,
 0, 2, 0, 0, 0x0106, "PhotometricInterpretation", &photometricInterpreter,
 0, 2, 0, 0, 0x0112, "Orientation", &stdInterpreter,
 0, 2, 0, 0, 0x0115, "SamplesPerPixel", &stdInterpreter,
 0, 2, 0, 0, 0x011C, "PlanarConfiguration", &planarConfigInterpreter,
 0, 2, 0, 0, 0x0212, "YCbCrSubSampling", &stdInterpreter,
 0, 2, 0, 0, 0x0213, "YCbCrPositioning", &stdInterpreter,
 0, 2, 0, 0, 0x011A, "XResolution", &stdInterpreter,
 0, 2, 0, 0, 0x011B, "YResolution", &stdInterpreter,
 0, 2, 0, 0, 0x0128, "ResolutionUnit", &stdInterpreter,
 1, 0, 0, 0, 0x0111, "StripOffsets", &stdInterpreter,
 1, 0, 0, 0, 0x0116, "RowsPerStrip", &stdInterpreter,
 1, 0, 0, 0, 0x0117, "StripByteCounts", &stdInterpreter,
 0, 2, 0, 0, 0x0201, "JPEGInterchangeFormat", &stdInterpreter,
 0, 2, 0, 0, 0x0202, "JPEGInterchangeFormatLength", &stdInterpreter,
 0, 2, 0, 0, 0x012D, "TransferFunction", &stdInterpreter,
 0, 2, 0, 0, 0x013E, "WhitePoint", &stdInterpreter,
 0, 2, 0, 0, 0x013F, "PriomaryChromaticities", &stdInterpreter,
 0, 2, 0, 0, 0x0211, "YCbCrCoefficients", &stdInterpreter,
 0, 2, 0, 0, 0x0214, "ReferenceBlackWhite", &stdInterpreter,
 0, 1, 0, 0, 0x0132, "DateTime", &stdInterpreter,
 0, 1, 1, 0, 0x010E, "ImageDescription", &stdInterpreter,
 0, 1, 0, 0, 0x010F, "Make", &stdInterpreter,
 0, 1, 0, 0, 0x0110, "Model", &stdInterpreter,
 0, 2, 0, 0, 0x0131, "Software", &stdInterpreter,
 0, 1, 1, 0, 0x013B, "Artist", &stdInterpreter,
 0, 1, 1, 0, 0x8298, "Copyright", &stdInterpreter,
 0, 1, 0, exifAttribs, 0x8769, "Exif", &stdInterpreter,
 0, 2, 0, 0, 0x8773, "ICCProfile", &stdInterpreter,
 0, 2, 0, 0, 0x83BB, "IPTCData", &stdInterpreter,
 0, 1, 0, gpsAttribs,  0x8825, "GPSInfo", &stdInterpreter,
 0, 1, 0, 0, 0x9003, "DateTimeOriginal", &stdInterpreter,
 0, 1, 0, 0, 0x9004, "DateTimeDigitized", &stdInterpreter,
 0, 1, 0, iopAttribs,  0xA005, "Interoperability", &stdInterpreter,
 1, 2, 0, ifdAttribs, 0x014A, "SubIFD", &stdInterpreter,
 0, 0, 0, 0, 0xC4A5, "PrintIMInformation", &stdInterpreter,
 0, 2, 0, 0, 0x00fe, "NewSubFileType", &stdInterpreter,
 -1, 0, 0,  0, 0, "", NULL};

};

/*#include <nikonattribs.h>
#include <canonattribs.h>
#include <pentaxattribs.h>
#include <olympusattribs.h>
#include <fujiattribs.h>
#include <sonyminoltaattribs.h>*/

#endif
