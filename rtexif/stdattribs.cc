/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Copyright (c)      2010 Oliver Duis <www.oliverduis.de>
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

#include <cstdio>
#include <cstring>

#include "rtexif.h"

namespace rtexif {

class ColorSpaceInterpreter : public ChoiceInterpreter {

    public:
        ColorSpaceInterpreter () {
            choices[1]      = "sRGB";
            choices[2]      = "Adobe RGB";
            choices[0xffff] = "Uncalibrated";
        }
};
ColorSpaceInterpreter colorSpaceInterpreter;

class PreviewColorSpaceInterpreter : public ChoiceInterpreter {

    public:
        PreviewColorSpaceInterpreter () {
            choices[0] = "Unknown";
            choices[1] = "Gray Gamma 2.2";
            choices[2] = "sRGB";
            choices[3] = "Adobe RGB";
            choices[4] = "ProPhoto RGB";
        }
};
PreviewColorSpaceInterpreter previewColorSpaceInterpreter;

class LinearSRGBInterpreter : public ChoiceInterpreter {

    public:
        LinearSRGBInterpreter () {
            choices[0] = "Linear";
            choices[1] = "sRGB";
        }
};
LinearSRGBInterpreter linearSRGBInterpreter;

class DefaultBlackRenderInterpreter : public ChoiceInterpreter {

    public:
        DefaultBlackRenderInterpreter () {
            choices[0] = "Auto";
            choices[1] = "None";
        }
};
DefaultBlackRenderInterpreter defaultBlackRenderInterpreter;

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

class ProfileEmbedPolicyInterpreter : public ChoiceInterpreter {

    public:
        ProfileEmbedPolicyInterpreter () {
            choices[0] = "Allow Copying";
            choices[1] = "Embed if Used";
            choices[2] = "Never Embed";
            choices[3] = "No Restrictions";
        }
};
ProfileEmbedPolicyInterpreter profileEmbedPolicyInterpreter;

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
        	char buffer[32];
        	double v = t->toDouble();
        	if( v < 0. || v > 1000. ) return "undef";
            sprintf (buffer, "%0.1f", v);
            return buffer;
        }
};
FNumberInterpreter fNumberInterpreter;

class ApertureInterpreter : public Interpreter {
    public:
        ApertureInterpreter () {}
        virtual std::string toString (Tag* t) {
        	char buffer[32];
        	double v = pow(2.0, t->toDouble()/2.0);
        	if( v < 0. || v > 1000. ) return "undef";
            sprintf (buffer, "%.1f", v );
            return buffer;
        }
};
ApertureInterpreter apertureInterpreter;

class ExposureBiasInterpreter : public Interpreter {
    public:
        ExposureBiasInterpreter () {}
        virtual std::string toString (Tag* t) {
        	char buffer[32];
        	double v = t->toDouble();
        	if( v < -1000. || v > 1000. ) return "undef";
            sprintf (buffer, "%+0.2f", v );
            return buffer;
        }
};
ExposureBiasInterpreter exposureBiasInterpreter;

class ShutterSpeedInterpreter : public Interpreter {
    public:
        ShutterSpeedInterpreter () {}
        virtual std::string toString (Tag* t) {
        	char buffer[1024];
            double d = pow (2.0, -t->toDouble());
            if (d > 0.0 && d < 0.9)
                sprintf (buffer, "1/%.0f", 1.0 / d);
            else
                sprintf (buffer, "%.1f", d);
            return buffer;
        }
};
ShutterSpeedInterpreter shutterSpeedInterpreter;

class ExposureTimeInterpreter : public Interpreter {
    public:
        ExposureTimeInterpreter () {}
        virtual std::string toString (Tag* t) {
        	char buffer[1024];
            double d = t->toDouble();
            if (d > 0.0 && d < 0.9)
                sprintf (buffer, "1/%.0f", 1.0 / d);
            else
                sprintf (buffer, "%.1f", d);
            return buffer;
        }
};
ExposureTimeInterpreter exposureTimeInterpreter;

class FocalLengthInterpreter : public Interpreter {
    public:
        FocalLengthInterpreter () {}
        virtual std::string toString (Tag* t) {
        	char buffer[32];
        	double v = t->toDouble();
        	if( v>1000000. || v<0 ) return "undef";
            sprintf (buffer, "%.1f", v );
            return buffer;
        }
};
FocalLengthInterpreter focalLengthInterpreter;

class UserCommentInterpreter : public Interpreter {
    public:
        UserCommentInterpreter () {}
        virtual std::string toString (Tag* t) {
            char *buffer = new char[t->getCount()];
            if (!strncmp((char*)t->getValue(), "ASCII\0\0\0",8))
                strncpy (buffer, (char*)t->getValue()+8, t->getCount()-8);
            else
                buffer[0]=0;
            std::string retVal(buffer);
            delete [] buffer;
            return retVal;
        }
        virtual void fromString (Tag* t, const std::string& value) {
            char *buffer = new char[t->getCount()];
            memcpy (buffer, "ASCII\0\0\0", 8);
            strcpy (buffer+8, value.c_str());
            t->fromString (buffer, value.size() + 9);
            delete [] buffer;
        }
};
UserCommentInterpreter userCommentInterpreter;

class CFAInterpreter : public Interpreter {
public:
	CFAInterpreter(){}
	virtual std::string toString (Tag* t) {
		char colors[]="RGB";
		char buffer[1024];
		for( int i=0; i< t->getCount();i++){
			unsigned char c = t->toInt(i,BYTE);
			buffer[i]= c<3 ?colors[c]:' ';
		}
		buffer[t->getCount()]=0;
		return buffer;
	}
};
CFAInterpreter cfaInterpreter;

class OrientationInterpreter : public ChoiceInterpreter {
public:
	OrientationInterpreter (){
        choices[1] = "Horizontal (normal)";
        choices[2] = "Mirror horizontal ";
        choices[3] = "Rotate 180";
        choices[4] = "Mirror vertical";
        choices[5] = "Mirror horizontal and rotate 270 CW";
        choices[6] = "Rotate 90 CW";
        choices[7] = "Mirror horizontal and rotate 90 CW";
        choices[8] = "Rotate 270 CW";
        // '9' is an "unofficial" value for Orientation but used by some older cameras that lacks orientation sensor, such as Kodak DCS
        choices[9] = "Unknown";
	}
};
OrientationInterpreter orientationInterpreter;

class UnitsInterpreter : public ChoiceInterpreter {
public:
	UnitsInterpreter(){
        choices[0] = "Unknown";
        choices[1] = "inches";
        choices[2] = "cm";
	}
};
UnitsInterpreter unitsInterpreter;

class UTF8BinInterpreter : public Interpreter {
    public:
        UTF8BinInterpreter () {}
};
UTF8BinInterpreter utf8BinInterpreter;

const TagAttrib exifAttribs[] = {
 {0, AC_SYSTEM,    0, 0, 0x0100, AUTO, "ImageWidth", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0101, AUTO, "ImageHeight", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0102, AUTO, "BitsPerSample", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0103, AUTO, "Compression", &compressionInterpreter},
 {0, AC_WRITE,     0, 0, 0x828d, AUTO, "CFAPatternDim", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x828e, AUTO, "CFAPattern", &cfaInterpreter},
 {0, AC_WRITE,     0, 0, 0x829A, AUTO, "ExposureTime", &exposureTimeInterpreter},
 {0, AC_WRITE,     0, 0, 0x829D, AUTO, "FNumber", &fNumberInterpreter},
 {0, AC_WRITE,     0, 0, 0x8822, AUTO, "ExposureProgram", &exposureProgramInterpreter},
 {0, AC_WRITE,     0, 0, 0x8824, AUTO, "SpectralSensitivity", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x8827, AUTO, "ISOSpeedRatings", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x8828, AUTO, "OECF", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x9000, AUTO, "ExifVersion", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x9003, AUTO, "DateTimeOriginal", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x9004, AUTO, "DateTimeDigitized", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x9101, AUTO, "ComponentsConfiguration", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x9102, AUTO, "CompressedBitsPerPixel", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x9201, AUTO, "ShutterSpeedValue", &shutterSpeedInterpreter},
 {0, AC_WRITE,     0, 0, 0x9202, AUTO, "ApertureValue", &apertureInterpreter},
 {0, AC_WRITE,     0, 0, 0x9203, AUTO, "BrightnessValue", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x9204, AUTO, "ExposureBiasValue", &exposureBiasInterpreter},
 {0, AC_WRITE,     0, 0, 0x9205, AUTO, "MaxApertureValue", &apertureInterpreter},
 {0, AC_WRITE,     0, 0, 0x9206, AUTO, "SubjectDistance", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x9207, AUTO, "MeteringMode", &meteringModeInterpreter},
 {0, AC_WRITE,     0, 0, 0x9208, AUTO, "LightSource", &lightSourceInterpreter},
 {0, AC_WRITE,     0, 0, 0x9209, AUTO, "Flash", &flashInterpreter},
 {0, AC_WRITE,     0, 0, 0x920A, AUTO, "FocalLength", &focalLengthInterpreter},
 {0, AC_WRITE,     0, 0, 0x9214, AUTO, "SubjectArea", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x9215, AUTO, "ExposureIndex", &stdInterpreter}, // Note: exists as 0xA215 too, it should be that way
 {0, AC_DONTWRITE, 0, 0, 0x9216, AUTO, "TIFFEPSStandardID", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x9217, AUTO, "SensingMethod", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x927C, AUTO, "MakerNote", &stdInterpreter},
 {0, AC_WRITE,     1, 0, 0x9286, AUTO, "UserComment", &userCommentInterpreter},
 {0, AC_WRITE,     0, 0, 0x9290, AUTO, "SubSecTime", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x9291, AUTO, "SubSecTimeOriginal", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x9292, AUTO, "SubSecTimeDigitized", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0xA000, AUTO, "FlashpixVersion", &stdInterpreter},
 {0, AC_DONTWRITE, 0, 0, 0xA001, AUTO, "ColorSpace", &colorSpaceInterpreter},
 {0, AC_SYSTEM,    0, 0, 0xA002, AUTO, "PixelXDimension", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0xA003, AUTO, "PixelYDimension", &stdInterpreter},
 {1, AC_DONTWRITE, 0, 0, 0xA004, AUTO, "RelatedSoundFile", &stdInterpreter},
 {0, AC_SYSTEM,    0, iopAttribs,  0xA005, AUTO, "Interoperability", &stdInterpreter},  // do not enable, as it causes trouble with FUJI files
 {0, AC_WRITE,     0, 0, 0xA20B, AUTO, "FlashEnergy", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA20C, AUTO, "SpatialFrequencyResponse", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA20E, AUTO, "FocalPlaneXResolution", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA20F, AUTO, "FocalPlaneYResolution", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA210, AUTO, "FocalPlaneResolutionUnit", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA214, AUTO, "SubjectLocation", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA215, AUTO, "ExposureIndex", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA217, AUTO, "SensingMethod", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA300, AUTO, "FileSource", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA301, AUTO, "SceneType", &stdInterpreter},
 {0, AC_DONTWRITE, 0, 0, 0xA302, AUTO, "CFAPattern", &cfaInterpreter},
 {0, AC_WRITE,     0, 0, 0xA401, AUTO, "CustomRendered", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA402, AUTO, "ExposureMode", &exposureModeInterpreter},
 {0, AC_WRITE,     0, 0, 0xA403, AUTO, "WhiteBalance", &whiteBalanceInterpreter},
 {0, AC_WRITE,     0, 0, 0xA404, AUTO, "DigitalZoomRatio", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA405, AUTO, "FocalLengthIn35mmFilm", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA406, AUTO, "SceneCaptureType", &sceneCaptureInterpreter},
 {0, AC_WRITE,     0, 0, 0xA407, AUTO, "GainControl", &gainControlInterpreter},
 {0, AC_WRITE,     0, 0, 0xA408, AUTO, "Contrast", &contrastInterpreter},
 {0, AC_WRITE,     0, 0, 0xA409, AUTO, "Saturation", &saturationInterpreter},
 {0, AC_WRITE,     0, 0, 0xA40A, AUTO, "Sharpness", &sharpnessInterpreter},
 {0, AC_WRITE,     0, 0, 0xA40B, AUTO, "DeviceSettingDescription", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA40C, AUTO, "SubjectDistanceRange", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA420, AUTO, "ImageUniqueID", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA431, AUTO, "SerialNumber", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA432, AUTO, "LensInfo", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA433, AUTO, "LensMake", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA434, AUTO, "LensModel", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA435, AUTO, "LensSerialNumber", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xA500, AUTO, "Gamma", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC618, AUTO, "LinearizationTable", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC619, AUTO, "BlackLevelRepeatDim", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC61A, AUTO, "BlackLevel", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC61B, AUTO, "BlackLevelDeltaH", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC61C, AUTO, "BlackLevelDeltaV", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC61D, AUTO, "WhiteLevel", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC61E, AUTO, "DefaultScale", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC61F, AUTO, "DefaultCropOrigin", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC620, AUTO, "DefaultCropSize", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC621, AUTO, "ColorMatrix1", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC622, AUTO, "ColorMatrix2", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC623, AUTO, "CameraCalibration1", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC624, AUTO, "CameraCalibration2", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC625, AUTO, "ReductionMatrix1", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC626, AUTO, "ReductionMatrix2", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC627, AUTO, "AnalogBalance", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC628, AUTO, "AsShotNeutral", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC629, AUTO, "AsShotWhiteXY", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC62A, AUTO, "BaselineExposure", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC62B, AUTO, "BaselineNoise", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC62C, AUTO, "BaselineSharpness", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC62D, AUTO, "BayerGreenSplit", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC62E, AUTO, "LinearResponseLimit", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC62F, AUTO, "CameraSerialNumber", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC630, AUTO, "DNGLensInfo", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC631, AUTO, "ChromaBlurRadius", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC632, AUTO, "AntiAliasStrength", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC633, AUTO, "ShadowScale", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC65A, AUTO, "CalibrationIlluminant1", &lightSourceInterpreter},
 {0, AC_WRITE,     0, 0, 0xC65B, AUTO, "CalibrationIlluminant2", &lightSourceInterpreter},
 {0, AC_WRITE,     0, 0, 0xC65C, AUTO, "BestQualityScale", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC65D, AUTO, "RawDataUniqueID", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC68B, AUTO, "OriginalRawFileName", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC68D, AUTO, "ActiveArea", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC68E, AUTO, "MaskedAreas", &stdInterpreter},
// {0, AC_WRITE,     0, 0, 0xC68F, AUTO, "AsShotICCProfile", & ???},
 {0, AC_WRITE,     0, 0, 0xC690, AUTO, "AsShotPreProfileMatrix", &stdInterpreter},
// {0, AC_WRITE,     0, 0, 0xC691, AUTO, "CurrentICCProfile", & ???},
 {0, AC_WRITE,     0, 0, 0xC692, AUTO, "CurrentPreProfileMatrix", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC6BF, AUTO, "ColorimetricReference", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC6F3, AUTO, "CameraCalibrationSig", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC6F4, AUTO, "ProfileCalibrationSig", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC6F5, AUTO, "ProfileIFD", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC6F6, AUTO, "AsShotProfileName", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC6F7, AUTO, "NoiseReductionApplied", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC6F8, AUTO, "ProfileName", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC6F9, AUTO, "ProfileHueSatMapDims", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC6FA, AUTO, "ProfileHueSatMapData1", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC6FB, AUTO, "ProfileHueSatMapData2", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC6FC, AUTO, "ProfileToneCurve", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC6FD, AUTO, "ProfileEmbedPolicy", &profileEmbedPolicyInterpreter},
 {0, AC_WRITE,     0, 0, 0xC6FE, AUTO, "ProfileCopyright", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC714, AUTO, "ForwardMatrix1", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC715, AUTO, "ForwardMatrix2", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC716, AUTO, "PreviewApplicationName", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC717, AUTO, "PreviewApplicationVersion", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC718, AUTO, "PreviewSettingsName", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC719, AUTO, "PreviewSettingsDigest", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC71A, AUTO, "PreviewColorSpace", &previewColorSpaceInterpreter},
 {0, AC_WRITE,     0, 0, 0xC71B, AUTO, "PreviewDateTime", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC71C, AUTO, "RawImageDigest", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC71D, AUTO, "OriginalRawFileDigest", &stdInterpreter},
// {0, AC_WRITE,     0, 0, 0xC71E, AUTO, "SubTileBlockSize", & ???},
// {0, AC_WRITE,     0, 0, 0xC71F, AUTO, "RowInterleaveFactor", & ???},
 {0, AC_WRITE,     0, 0, 0xC725, AUTO, "ProfileLookTableDims", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC726, AUTO, "ProfileLookTableData", &stdInterpreter},
// {0, AC_WRITE,     0, 0, 0xC740, AUTO, "OpcodeList1", & ???},
// {0, AC_WRITE,     0, 0, 0xC741, AUTO, "OpcodeList2", & ???},
// {0, AC_WRITE,     0, 0, 0xC74E, AUTO, "OpcodeList3", & ???},
 {0, AC_WRITE,     0, 0, 0xC761, AUTO, "NoiseProfile", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC763, AUTO, "TimeCodes", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC764, AUTO, "FrameRate", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC772, AUTO, "TStop", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC789, AUTO, "ReelName", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC791, AUTO, "OriginalDefaultFinalSize", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC792, AUTO, "OriginalBestQualitySize", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC793, AUTO, "OriginalDefaultCropSize", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC7A1, AUTO, "CameraLabel", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC7A3, AUTO, "ProfileHueSatMapEncoding", &linearSRGBInterpreter},
 {0, AC_WRITE,     0, 0, 0xC7A4, AUTO, "ProfileLookTableEncoding", &linearSRGBInterpreter},
 {0, AC_WRITE,     0, 0, 0xC7A5, AUTO, "BaselineExposureOffset", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC7A6, AUTO, "DefaultBlackRender", &defaultBlackRenderInterpreter},
 {0, AC_WRITE,     0, 0, 0xC7A7, AUTO, "NewRawImageDigest", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC7A8, AUTO, "RawToPreviewGain", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC7B5, AUTO, "DefaultUserCrop", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xFDE9, AUTO, "SerialNumber", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xFDEA, AUTO, "Lens", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xFE4C, AUTO, "RawFile", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xFE4D, AUTO, "Converter", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xFE4E, AUTO, "WhiteBalance", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xFE51, AUTO, "Exposure", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xFE52, AUTO, "Shadows", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xFE53, AUTO, "Brightness", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xFE54, AUTO, "Contrast", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xFE55, AUTO, "Saturation", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xFE56, AUTO, "Sharpness", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xFE57, AUTO, "Smoothness", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xFE58, AUTO, "MoireFilter", &stdInterpreter},
 {-1, AC_DONTWRITE, 0, 0, 0, AUTO, "", NULL }};


const TagAttrib gpsAttribs[] = {
 {0, AC_WRITE, 0, 0, 0x0000, AUTO, "GPSVersionID", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0001, AUTO, "GPSLatitudeRef", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0002, AUTO, "GPSLatitude", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0003, AUTO, "GPSLongitudeRef", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0004, AUTO, "GPSLongitude", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0005, AUTO, "GPSAltitudeRef", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0006, AUTO, "GPSAltitude", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0007, AUTO, "GPSTimeStamp", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0008, AUTO, "GPSSatelites", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0009, AUTO, "GPSStatus", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x000a, AUTO, "GPSMeasureMode", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x000b, AUTO, "GPSDOP", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x000c, AUTO, "GPSSpeedRef", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x000d, AUTO, "GPSSpeed", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x000e, AUTO, "GPSTrackRef", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x000f, AUTO, "GPSTrack", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0010, AUTO, "GPSImgDirectionRef", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0011, AUTO, "GPSImgDirection", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0012, AUTO, "GPSMapDatum", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0013, AUTO, "GPSDestLatitudeRef", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0014, AUTO, "GPSDestLatitude", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0015, AUTO, "GPSDestLongitudeRef", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0016, AUTO, "GPSDestLongitude", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0017, AUTO, "GPSDestBearingRef", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0018, AUTO, "GPSDestBearing", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0019, AUTO, "GPSDestDistanceRef", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x001a, AUTO, "GPSDestDistance", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x001b, AUTO, "GPSProcessingMethod", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x001c, AUTO, "GPSAreaInformation", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x001d, AUTO, "GPSDateStamp", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x001e, AUTO, "GPSDifferential", &stdInterpreter},
 {-1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL }};

const TagAttrib iopAttribs[] = {
 {0, AC_WRITE, 0, 0, 0x0001, AUTO, "InteroperabilityIndex", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0002, AUTO, "InteroperabilityVersion", &stdInterpreter},
 {-1, AC_DONTWRITE, 0, 0, 0, AUTO, "", NULL }};

 const TagAttrib ifdAttribs[] = { 
 {0, AC_SYSTEM,    0, 0, 0x0017, AUTO, "PanaISO", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0100, AUTO, "ImageWidth", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0101, AUTO, "ImageHeight", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0102, AUTO, "BitsPerSample", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0103, AUTO, "Compression", &compressionInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0106, AUTO, "PhotometricInterpretation", &photometricInterpreter},
 {0, AC_WRITE,     1, 0, 0x010E, AUTO, "ImageDescription", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x010F, AUTO, "Make", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x0110, AUTO, "Model", &stdInterpreter},
 {1, AC_DONTWRITE, 0, 0, 0x0111, AUTO, "StripOffsets", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0112, AUTO, "Orientation", &orientationInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0115, AUTO, "SamplesPerPixel", &stdInterpreter},
 {1, AC_DONTWRITE, 0, 0, 0x0116, AUTO, "RowsPerStrip", &stdInterpreter},
 {1, AC_DONTWRITE, 0, 0, 0x0117, AUTO, "StripByteCounts", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x011A, AUTO, "XResolution", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x011B, AUTO, "YResolution", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x011C, AUTO, "PlanarConfiguration", &planarConfigInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0128, AUTO, "ResolutionUnit", &unitsInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x012D, AUTO, "TransferFunction", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0131, AUTO, "Software", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x0132, AUTO, "DateTime", &stdInterpreter},
 {0, AC_WRITE,     1, 0, 0x013B, AUTO, "Artist", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x013E, AUTO, "WhitePoint", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x013F, AUTO, "PriomaryChromaticities", &stdInterpreter},
 {0, AC_WRITE,     0, ifdAttribs, 0x014A, AUTO, "SubIFD", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0201, AUTO, "JPEGInterchangeFormat", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0202, AUTO, "JPEGInterchangeFormatLength", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0211, AUTO, "YCbCrCoefficients", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0212, AUTO, "YCbCrSubSampling", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0213, AUTO, "YCbCrPositioning", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x0214, AUTO, "ReferenceBlackWhite", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x02bc, AUTO, "ApplicationNotes", &utf8BinInterpreter},  // XMP
 {0, AC_WRITE,     0, 0, 0x4746, AUTO, "Rating",&stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x4749, AUTO, "RatingPercent",&stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x828d, AUTO, "CFAPatternDim", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x828e, AUTO, "CFAPattern", &cfaInterpreter},
 {0, AC_WRITE,     0, kodakIfdAttribs, 0x8290, AUTO, "KodakIFD", &stdInterpreter},
 {0, AC_WRITE,     1, 0, 0x8298, AUTO, "Copyright", &stdInterpreter},
 {0, AC_DONTWRITE, 0, 0, 0x8606, AUTO, "LeafData", &stdInterpreter}, // is actually a subdir, but a proprietary format
 {0, AC_WRITE,     0, exifAttribs, 0x8769, AUTO, "Exif", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x8773, AUTO, "ICCProfile", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x83BB, AUTO, "IPTCData", &stdInterpreter},
 {0, AC_WRITE,     0, gpsAttribs,  0x8825, AUTO, "GPSInfo", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x9003, AUTO, "DateTimeOriginal", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x9004, AUTO, "DateTimeDigitized", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0x9211, AUTO, "ImageNumber", &stdInterpreter},
 {0, AC_WRITE,     0, iopAttribs,  0xA005, AUTO, "Interoperability", &stdInterpreter},
 {0, AC_DONTWRITE, 0, 0, 0xC4A5, AUTO, "PrintIMInformation", &stdInterpreter},
 {0, AC_DONTWRITE, 0, 0, 0xC612, AUTO, "DNGVersion", &stdInterpreter},
 {0, AC_DONTWRITE, 0, 0, 0xC613, AUTO, "DNGBackwardVersion", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xC614, AUTO, "UniqueCameraModel", &stdInterpreter},
 {0, AC_WRITE,     0, 0, 0xc62f, AUTO, "CameraSerialNumber", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0xc630, AUTO, "DNGLensInfo", &stdInterpreter},
 {0, AC_DONTWRITE, 0, 0, 0xC634, AUTO, "MakerNote", &stdInterpreter}, //DNGPrivateData
 {0, AC_WRITE,     0, 0, 0xc65d, AUTO, "RawDataUniqueID", &stdInterpreter},
 {0, AC_DONTWRITE, 0, 0, 0xc761, AUTO, "NoiseProfile", &stdInterpreter},
 {0, AC_SYSTEM,    0, 0, 0x00fe, AUTO, "NewSubFileType", &stdInterpreter},
 {-1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}};
}

#endif
