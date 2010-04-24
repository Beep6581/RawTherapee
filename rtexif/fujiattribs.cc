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

#include <rtexif.h>
#include <string>
#include <map>
#include <math.h>
#include <sstream>
#include <iomanip>

namespace rtexif {

class FAOnOffInterpreter : public ChoiceInterpreter {
    public:
        FAOnOffInterpreter () {
            choices[0]      = "Off";
            choices[1]      = "On";
        }
};
FAOnOffInterpreter faOnOffInterpreter;

class FASharpnessInterpreter : public ChoiceInterpreter {
    public:
        FASharpnessInterpreter () {
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

class FAWhiteBalanceInterpreter : public ChoiceInterpreter {
    public:
        FAWhiteBalanceInterpreter () {
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
            choices[0xf00]  = "Custom";
            choices[0xf01]  = "Custom2";
            choices[0xf02]  = "Custom3";
            choices[0xf03]  = "Custom4";
            choices[0xf04]  = "Custom5";
            choices[0xff0]  = "Kelvin";
        }
};
FAWhiteBalanceInterpreter faWhiteBalanceInterpreter;

class FASaturationInterpreter : public ChoiceInterpreter {
    public:
        FASaturationInterpreter () {
            choices[0]      = "Normal";
            choices[0x80]   = "Medium High";
            choices[0x100]  = "High";
            choices[0x180]  = "Medium Low";
            choices[0x200]  = "Low";
            choices[0x300]  = "None (B&W)";
            choices[0x8000] = "Film Simulation";
        }
};
FASaturationInterpreter faSaturationInterpreter;

class FAContrastInterpreter : public ChoiceInterpreter {
    public:
        FAContrastInterpreter () {
            choices[0]      = "Normal";
            choices[0x80]   = "Medium High";
            choices[0x100]  = "High";
            choices[0x180]  = "Medium Low";
            choices[0x200]  = "Low";
            choices[0x8000] = "Film Simulation";
        }
};
FAContrastInterpreter faContrastInterpreter;

class FAContrast2Interpreter : public ChoiceInterpreter {
    public:
        FAContrast2Interpreter () {
            choices[0]      = "Normal";
            choices[0x100]  = "High";
            choices[0x300]  = "Low";
        }
};
FAContrast2Interpreter faContrast2Interpreter;

class FANoiseReductionInterpreter : public ChoiceInterpreter {
    public:
        FANoiseReductionInterpreter () {
            choices[0x40] = "Low";
            choices[0x80] = "Normal";
        }
};
FANoiseReductionInterpreter faNoiseReductionInterpreter;

class FAFlashInterpreter : public ChoiceInterpreter {
    public:
        FAFlashInterpreter () {
            choices[0]  = "Auto";
            choices[1]  = "On";
            choices[2]  = "Off";
            choices[3]  = "Red-eye reduction";
            choices[4]  = "External";
        }
};
FAFlashInterpreter faFlashInterpreter;

class FAFocusModeInterpreter : public ChoiceInterpreter {
    public:
        FAFocusModeInterpreter () {
            choices[0]  = "Auto";
            choices[1]  = "Manual";
        }
};
FAFocusModeInterpreter faFocusModeInterpreter;

class FAColorModeInterpreter : public ChoiceInterpreter {
    public:
        FAColorModeInterpreter () {
            choices[0]    = "Standard";
            choices[0x10] = "Chrome";
            choices[0x30] = "B & W";
        }
};
FAColorModeInterpreter faColorModeInterpreter;

class FADynamicRangeInterpreter : public ChoiceInterpreter {
    public:
        FADynamicRangeInterpreter () {
            choices[1]  = "Standard";
            choices[3]  = "Wide";
        }
};
FADynamicRangeInterpreter faDynamicRangeInterpreter;

class FAFilmModeInterpreter : public ChoiceInterpreter {
    public:
        FAFilmModeInterpreter () {
            choices[0]      = "F0/Standard";
            choices[0x100]  = "F1/Studio Portrait";
            choices[0x110]  = "F1a/Studio Portrait Enhanced Saturation";
            choices[0x120]  = "F1b/Studio Portrait Smooth Skin Tone";
            choices[0x130]  = "F1c/Studio Portrait Increased Sharpness ";
            choices[0x200]  = "F2/Fujichrome";
            choices[0x300]  = "F3/Studio Portrait Ex";
            choices[0x400]  = "F4/Velvia";
        }
};
FAFilmModeInterpreter faFilmModeInterpreter;

class FADRSettingInterpreter : public ChoiceInterpreter {
    public:
        FADRSettingInterpreter () {
            choices[0]      = "Auto (100-400%)";
            choices[0x1]    = "RAW";
            choices[0x100]  = "Standard (100%)";
            choices[0x200]  = "Wide1 (230%)";
            choices[0x201]  = "Wide2 (400%)";
            choices[0x8000] = "Film Simulation";
        }
};
FADRSettingInterpreter faDRSettingInterpreter;

class FAPictureModeInterpreter : public ChoiceInterpreter {
    public:
        FAPictureModeInterpreter () {
            choices[0]     = "Auto";
            choices[1]     = "Portrait";
            choices[2]     = "Landscape";
            choices[3]     = "Macro";
            choices[4]     = "Sports";
            choices[5]     = "Night Scene";
            choices[6]     = "Program AE";
            choices[7]     = "Natural Light";
            choices[8]     = "Anti-blur";
            choices[9]     = "Beach & Snow";
            choices[10]    = "Sunset";
            choices[11]    = "Museum";
            choices[12]    = "Party";
            choices[13]    = "Flower";
            choices[14]    = "Text";
            choices[15]    = "Natural Light & Flash";
            choices[16]    = "Beach";
            choices[17]    = "Fireworks";
            choices[18]    = "Underwater";
            choices[0x100] = "Aperture-priority AE";
            choices[0x200] = "Shutter speed priority AE";
            choices[0x300] = "Manual";
        }
};
FAPictureModeInterpreter faPictureModeInterpreter;



const TagAttrib fujiAttribs[] = {
 0, 1, 0, 0, 0x0000, "Version", &stdInterpreter,
 0, 1, 0, 0, 0x0010, "InternalSerialNumber", &stdInterpreter,
 0, 1, 0, 0, 0x1000, "Quality", &stdInterpreter,
 0, 1, 0, 0, 0x1001, "Sharpness", &faSharpnessInterpreter,
 0, 1, 0, 0, 0x1002, "WhiteBalance", &faWhiteBalanceInterpreter,
 0, 1, 0, 0, 0x1003, "Saturation", &faSaturationInterpreter,
 0, 1, 0, 0, 0x1004, "Contrast", &faContrastInterpreter,
 0, 1, 0, 0, 0x1005, "ColorTemperature", &stdInterpreter,
 0, 1, 0, 0, 0x1006, "Contrast2", &faContrast2Interpreter,
 0, 1, 0, 0, 0x100a, "WhiteBalanceFineTune", &stdInterpreter,
 0, 1, 0, 0, 0x100b, "NoiseReduction", &faNoiseReductionInterpreter,
 0, 1, 0, 0, 0x100b, "FujiFlashMode", &faFlashInterpreter,
 0, 1, 0, 0, 0x1011, "FlashExposureComp", &stdInterpreter,
 0, 1, 0, 0, 0x1020, "Macro", &faOnOffInterpreter,
 0, 1, 0, 0, 0x1021, "FocusMode", &faFocusModeInterpreter,
 0, 1, 0, 0, 0x1023, "FocusPixel", &stdInterpreter,
 0, 1, 0, 0, 0x1030, "SlowSync", &faOnOffInterpreter,
 0, 1, 0, 0, 0x1031, "PictureMode", &faPictureModeInterpreter,
 0, 1, 0, 0, 0x1100, "AutoBracketing", &faOnOffInterpreter,
 0, 1, 0, 0, 0x1101, "SequenceNumber", &stdInterpreter,
 0, 1, 0, 0, 0x1210, "ColorMode", &faColorModeInterpreter,
 0, 1, 0, 0, 0x1300, "BlurWarning", &faOnOffInterpreter,
 0, 1, 0, 0, 0x1301, "FocusWarning", &faOnOffInterpreter,
 0, 1, 0, 0, 0x1302, "ExposureWarning", &faOnOffInterpreter,
 0, 1, 0, 0, 0x1400, "DynamicRange", &faDynamicRangeInterpreter,
 0, 1, 0, 0, 0x1401, "FilmMode", &faFilmModeInterpreter,
 0, 1, 0, 0, 0x1402, "DynamicRangeSetting", &faDRSettingInterpreter,
 0, 1, 0, 0, 0x1403, "DevelopmentDynamicRange", &stdInterpreter,
 0, 1, 0, 0, 0x1404, "MinFocalLength", &stdInterpreter,
 0, 1, 0, 0, 0x1405, "MaxFocalLength", &stdInterpreter,
 0, 1, 0, 0, 0x1406, "MaxApertureAtMinFocal", &stdInterpreter,
 0, 1, 0, 0, 0x1407, "MaxApertureAtMaxFocal", &stdInterpreter,
 0, 1, 0, 0, 0x8000, "FileSource", &stdInterpreter,
 0, 1, 0, 0, 0x8002, "OrderNumber", &stdInterpreter,
 0, 1, 0, 0, 0x8003, "FrameNumber", &stdInterpreter,
 -1, 0, 0,  0, 0, "", NULL};
};
#endif

