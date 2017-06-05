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
#ifndef _FUJIATTRIBS_
#define _FUJIATTRIBS_

#include "rtexif.h"

namespace rtexif
{

class FAOnOffInterpreter : public ChoiceInterpreter
{
public:
    FAOnOffInterpreter ()
    {
        choices[0]      = "Off";
        choices[1]      = "On";
    }
};
FAOnOffInterpreter faOnOffInterpreter;

class FASharpnessInterpreter : public ChoiceInterpreter
{
public:
    FASharpnessInterpreter ()
    {
        choices[1]      = "Soft";
        choices[2]      = "Soft2";
        choices[3]      = "Normal";
        choices[4]      = "Hard";
        choices[5]      = "Hard2";
        choices[0x82]   = "Medium Soft";
        choices[0x84]   = "Medium Hard";
        choices[0x8000] = "Film Simulation";
        choices[0xffff] = "n/a";
    }
};
FASharpnessInterpreter faSharpnessInterpreter;

class FAWhiteBalanceInterpreter : public ChoiceInterpreter
{
public:
    FAWhiteBalanceInterpreter ()
    {
        choices[0]      = "Auto";
        choices[0x100]  = "Daylight";
        choices[0x200]  = "Cloudy";
        choices[0x300]  = "Daylight Fluorescent";
        choices[0x301]  = "Day White Fluorescent";
        choices[0x302]  = "White Fluorescent";
        choices[0x303]  = "Warm White Fluorescent";
        choices[0x304]  = "Living Room Warm White Fluorescent";
        choices[0x400]  = "Incandescent";
        choices[0x500]  = "Flash";
        choices[0x600]  = "Underwater";
        choices[0xf00]  = "Custom";
        choices[0xf01]  = "Custom2";
        choices[0xf02]  = "Custom3";
        choices[0xf03]  = "Custom4";
        choices[0xf04]  = "Custom5";
        choices[0xff0]  = "Kelvin";
    }
};
FAWhiteBalanceInterpreter faWhiteBalanceInterpreter;

class FASaturationInterpreter : public ChoiceInterpreter
{
public:
    FASaturationInterpreter ()
    {
        choices[0] = "Normal";
        choices[128] = "Medium High";
        choices[256] = "High";
        choices[384] = "Medium Low";
        choices[512] = "Low";
        choices[768] = "None (B&W)";
        choices[769] = "B&W Red Filter";
        choices[770] = "B&W Yellow Filter";
        choices[771] = "B&W Green Filter";
        choices[784] = "B&W Sepia";
        choices[1024] = "Low 2";
        choices[1280] = "Acros";
        choices[1281] = "Acros Red Filter";
        choices[1282] = "Acros Yellow Filter";
        choices[1283] = "Acros Green Filter";
        choices[32768] = "Film Simulation";
    }
};
FASaturationInterpreter faSaturationInterpreter;

class FAContrastInterpreter : public ChoiceInterpreter
{
public:
    FAContrastInterpreter ()
    {
        choices[0]      = "Normal";
        choices[0x80]   = "Medium High";
        choices[0x100]  = "High";
        choices[0x180]  = "Medium Low";
        choices[0x200]  = "Low";
        choices[0x8000] = "Film Simulation";
    }
};
FAContrastInterpreter faContrastInterpreter;

class FAContrast2Interpreter : public ChoiceInterpreter
{
public:
    FAContrast2Interpreter ()
    {
        choices[0]      = "Normal";
        choices[0x100]  = "High";
        choices[0x300]  = "Low";
    }
};
FAContrast2Interpreter faContrast2Interpreter;

class FANoiseReductionInterpreter : public ChoiceInterpreter
{
public:
    FANoiseReductionInterpreter ()
    {
        choices[0x40]  = "Low";
        choices[0x80]  = "Normal";
        choices[0x100] = "n/a";
    }
};
FANoiseReductionInterpreter faNoiseReductionInterpreter;

class FAFlashInterpreter : public ChoiceInterpreter
{
public:
    // FujiFlashMode
    FAFlashInterpreter ()
    {
        choices[0]  = "Auto";
        choices[1]  = "On";
        choices[2]  = "Off";
        choices[3]  = "Red-eye reduction";
        choices[4]  = "External";
    }
};
FAFlashInterpreter faFlashInterpreter;

class FAFocusModeInterpreter : public ChoiceInterpreter
{
public:
    FAFocusModeInterpreter ()
    {
        choices[0]  = "Auto";
        choices[1]  = "Manual";
    }
};
FAFocusModeInterpreter faFocusModeInterpreter;

class FAColorModeInterpreter : public ChoiceInterpreter
{
public:
    FAColorModeInterpreter ()
    {
        choices[0]    = "Standard";
        choices[0x10] = "Chrome";
        choices[0x30] = "B & W";
    }
};
FAColorModeInterpreter faColorModeInterpreter;

class FADynamicRangeInterpreter : public ChoiceInterpreter
{
public:
    FADynamicRangeInterpreter ()
    {
        choices[1]  = "Standard";
        choices[3]  = "Wide";
    }
};
FADynamicRangeInterpreter faDynamicRangeInterpreter;

class FAFilmModeInterpreter : public ChoiceInterpreter
{
public:
    FAFilmModeInterpreter ()
    {
        choices[0x0]   = "F0/Standard (Provia)";
        choices[0x100] = "F1/Studio Portrait";
        choices[0x110] = "F1a/Studio Portrait Enhanced Saturation";
        choices[0x120] = "F1b/Studio Portrait Smooth Skin Tone (Astia)";
        choices[0x130] = "F1c/Studio Portrait Increased Sharpness";
        choices[0x200] = "F2/Fujichrome (Velvia)";
        choices[0x300] = "F3/Studio Portrait Ex";
        choices[0x400] = "F4/Velvia";
        choices[0x500] = "Pro Neg. Std";
        choices[0x501] = "Pro Neg. Hi";
        choices[0x600] = "Classic Chrome";
    }
};
FAFilmModeInterpreter faFilmModeInterpreter;

class FADRSettingInterpreter : public ChoiceInterpreter
{
public:
    // DynamicRangeSetting
    FADRSettingInterpreter ()
    {
        choices[0x0]    = "Auto (100-400%)";
        choices[0x1]    = "Manual";
        choices[0x100]  = "Standard (100%)";
        choices[0x200]  = "Wide1 (230%)";
        choices[0x201]  = "Wide2 (400%)";
        choices[0x8000] = "Film Simulation";
    }
};
FADRSettingInterpreter faDRSettingInterpreter;

class FAPictureModeInterpreter : public ChoiceInterpreter
{
public:
    FAPictureModeInterpreter ()
    {
        choices[0x0]   = "Auto";
        choices[0x1]   = "Portrait";
        choices[0x2]   = "Landscape";
        choices[0x3]   = "Macro";
        choices[0x4]   = "Sports";
        choices[0x5]   = "Night Scene";
        choices[0x6]   = "Program AE";
        choices[0x7]   = "Natural Light";
        choices[0x8]   = "Anti-blur";
        choices[0x9]   = "Beach & Snow";
        choices[0xa]   = "Sunset";
        choices[0xb]   = "Museum";
        choices[0xc]   = "Party";
        choices[0xd]   = "Flower";
        choices[0xe]   = "Text";
        choices[0xf]   = "Natural Light & Flash";
        choices[0x10]  = "Beach";
        choices[0x11]  = "Snow";
        choices[0x12]  = "Fireworks";
        choices[0x13]  = "Underwater";
        choices[0x14]  = "Portrait with Skin Correction";
        choices[0x16]  = "Panorama";
        choices[0x17]  = "Night (tripod)";
        choices[0x18]  = "Pro Low-light";
        choices[0x19]  = "Pro Focus";
        choices[0x1a]  = "Portrait 2";
        choices[0x1b]  = "Dog Face Detection";
        choices[0x1c]  = "Cat Face Detection";
        choices[0x40]  = "Advanced Filter";
        choices[0x100] = "Aperture-priority AE";
        choices[0x200] = "Shutter speed priority AE";
        choices[0x300] = "Manual";
    }
};
FAPictureModeInterpreter faPictureModeInterpreter;



const TagAttrib fujiAttribs[] = {
    {0, AC_WRITE, 0, nullptr, 0x0000, AUTO, "Version", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0010, AUTO, "InternalSerialNumber", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1000, AUTO, "Quality", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1001, AUTO, "Sharpness", &faSharpnessInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1002, AUTO, "WhiteBalance", &faWhiteBalanceInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1003, AUTO, "Saturation", &faSaturationInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1004, AUTO, "Contrast", &faContrastInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1005, AUTO, "ColorTemperature", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1006, AUTO, "Contrast2", &faContrast2Interpreter},
    {0, AC_WRITE, 0, nullptr, 0x100a, AUTO, "WhiteBalanceFineTune", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x100b, AUTO, "NoiseReduction", &faNoiseReductionInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1010, AUTO, "FujiFlashMode", &faFlashInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1011, AUTO, "FlashExposureComp", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1020, AUTO, "Macro", &faOnOffInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1021, AUTO, "FocusMode", &faFocusModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1023, AUTO, "FocusPixel", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1030, AUTO, "SlowSync", &faOnOffInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1031, AUTO, "PictureMode", &faPictureModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1100, AUTO, "AutoBracketing", &faOnOffInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1101, AUTO, "SequenceNumber", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1210, AUTO, "ColorMode", &faColorModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1300, AUTO, "BlurWarning", &faOnOffInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1301, AUTO, "FocusWarning", &faOnOffInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1302, AUTO, "ExposureWarning", &faOnOffInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1400, AUTO, "DynamicRange", &faDynamicRangeInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1401, AUTO, "FilmMode", &faFilmModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1402, AUTO, "DynamicRangeSetting", &faDRSettingInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1403, AUTO, "DevelopmentDynamicRange", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1404, AUTO, "MinFocalLength", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1405, AUTO, "MaxFocalLength", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1406, AUTO, "MaxApertureAtMinFocal", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x1407, AUTO, "MaxApertureAtMaxFocal", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x140b, AUTO, "AutoDynamicRange", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x4100, AUTO, "FacesDetected", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x8000, AUTO, "FileSource", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x8002, AUTO, "OrderNumber", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x8003, AUTO, "FrameNumber", &stdInterpreter},
    { -1, AC_DONTWRITE, 0,  nullptr, 0, AUTO, "", nullptr}
};
}
#endif

