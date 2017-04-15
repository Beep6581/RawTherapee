/*
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
#ifndef _SONYMINOLTAATTRIBS_
#define _SONYMINOLTAATTRIBS_

#include <cmath>

#include "rtexif.h"

namespace rtexif
{

class SANoYesInterpreter : public ChoiceInterpreter
{
public:
    SANoYesInterpreter ()
    {
        choices[1]  = "No";
        choices[16] = "Yes";
    }
};
SANoYesInterpreter saNoYesInterpreter;

class SAOnOffInterpreter : public ChoiceInterpreter
{
public:
    SAOnOffInterpreter ()
    {
        choices[0]      = "Off";
        choices[1]      = "On";
        choices[5]      = "On";
    }
};
SAOnOffInterpreter saOnOffInterpreter;

class SAOnOffInterpreter2 : public ChoiceInterpreter
{
public:
    SAOnOffInterpreter2 ()
    {
        choices[1]  = "Off";
        choices[16] = "On";
    }
};
SAOnOffInterpreter2 saOnOffInterpreter2;

class SAOnOffInterpreter3 : public ChoiceInterpreter
{
public:
    SAOnOffInterpreter3 ()
    {
        choices[1]  = "Off";
        choices[16] = "On (Auto)";
        choices[17] = "On (Manual)";
    }
};
SAOnOffInterpreter3 saOnOffInterpreter3;

class SAOnOffInterpreter4 : public ChoiceInterpreter
{
public:
    SAOnOffInterpreter4 ()
    {
        choices[0]   = "n/a";
        choices[1]   = "Off";
        choices[16]  = "On";
        choices[255] = "None";
    }
};
SAOnOffInterpreter4 saOnOffInterpreter4;

class SAOnOffInterpreter5 : public ChoiceInterpreter
{
public:
    SAOnOffInterpreter5 ()
    {
        choices[1]   = "On";
        choices[2]   = "Off";
    }
};
SAOnOffInterpreter5 saOnOffInterpreter5;

class SAHighISONoiseReduction : public ChoiceInterpreter
{
public:
    SAHighISONoiseReduction ()
    {
        choices[0]     = "Off";
        choices[1]     = "Low";
        choices[2]     = "Normal";
        choices[3]     = "High";
        choices[256]   = "Auto";
        choices[65535] = "n/a";
    }
};
SAHighISONoiseReduction saHighISONoiseReduction;

class SAHighISONoiseReduction2 : public ChoiceInterpreter
{
public:
    SAHighISONoiseReduction2 ()
    {
        choices[0]     = "Normal";
        choices[1]     = "High";
        choices[2]     = "Low";
        choices[3]     = "Off";
        choices[65535] = "n/a";
    }
};
SAHighISONoiseReduction2 saHighISONoiseReduction2;

class SAHighISONoiseReduction3 : public ChoiceInterpreter
{
public:
    SAHighISONoiseReduction3 ()
    {
        choices[0]     = "Normal";
        choices[1]     = "Low";
        choices[2]     = "High";
        choices[3]     = "Off";
    }
};
SAHighISONoiseReduction3 saHighISONoiseReduction3;

class SAHighISONoiseReduction4 : public ChoiceInterpreter
{
public:
    SAHighISONoiseReduction4 ()
    {
        choices[0]     = "Off";
        choices[1]     = "Low";
        choices[2]     = "Normal";
        choices[3]     = "High";
    }
};
SAHighISONoiseReduction4 saHighISONoiseReduction4;

class SAHighISONoiseReduction5 : public ChoiceInterpreter
{
public:
    SAHighISONoiseReduction5 ()
    {
        choices[16]    = "Low";
        choices[19]    = "Auto";
    }
};
SAHighISONoiseReduction5 saHighISONoiseReduction5;

class SASmileShutterMode : public ChoiceInterpreter
{
public:
    SASmileShutterMode ()
    {
        choices[17]    = "Slight smile";
        choices[18]    = "Normal smile";
        choices[19]    = "Big smile";
    }
};
SASmileShutterMode saSmileShutterMode;

class SAHDRLevel : public ChoiceInterpreter
{
public:
    SAHDRLevel ()
    {
        choices[33]    = "1 EV";
        choices[34]    = "1.5 EV";
        choices[35]    = "2 EV";
        choices[36]    = "2.5 EV";
        choices[37]    = "3 EV";
        choices[38]    = "3.5 EV";
        choices[39]    = "4 EV";
        choices[40]    = "5 EV";
        choices[41]    = "6 EV";
    }
};
SAHDRLevel saHDRLevel;

class SAViewingMode : public ChoiceInterpreter
{
public:
    SAViewingMode ()
    {
        choices[0]     = "n/a";
        choices[16]    = "ViewFinder";
        choices[33]    = "Focus Check Live View";
        choices[34]    = "Quick AF Live View";
    }
};
SAViewingMode saViewingMode;

class SAFlashAction : public ChoiceInterpreter
{
public:
    SAFlashAction ()
    {
        choices[1]    = "Did not fire";
        choices[2]    = "Fired";
    }
};
SAFlashAction saFlashAction;

class SALiveViewFocusMode : public ChoiceInterpreter
{
public:
    SALiveViewFocusMode ()
    {
        choices[0]    = "n/a";
        choices[1]    = "AF";
        choices[16]   = "Manual";
    }
};
SALiveViewFocusMode saLiveViewFocusMode;

class SALensMount : public ChoiceInterpreter
{
public:
    SALensMount ()
    {
        choices[1]    = "Unknown";
        choices[16]   = "A-Mount";
        choices[17]   = "E-Mount";
    }
};
SALensMount saLensMount;

class SASweepPanoramaSize : public ChoiceInterpreter
{
public:
    SASweepPanoramaSize ()
    {
        choices[1]    = "Standard";
        choices[2]    = "Wide";
    }
};
SASweepPanoramaSize saSweepPanoramaSize;

class SASweepPanoramaDirection : public ChoiceInterpreter
{
public:
    SASweepPanoramaDirection ()
    {
        choices[1]    = "Right";
        choices[2]    = "Left";
        choices[3]    = "Up";
        choices[4]    = "Down";
    }
};
SASweepPanoramaDirection saSweepPanoramaDirection;

class SALiveViewAFSetting : public ChoiceInterpreter
{
public:
    SALiveViewAFSetting ()
    {
        choices[0]    = "n/a";
        choices[1]    = "Phase-detect AF";
        choices[2]    = "Contrast AF";
    }
};
SALiveViewAFSetting saLiveViewAFSetting;

class SAPanoramaSize3D : public ChoiceInterpreter
{
public:
    SAPanoramaSize3D ()
    {
        choices[0]    = "n/a";
        choices[1]    = "Standard";
        choices[2]    = "Wide";
        choices[3]    = "16:9";
    }
};
SAPanoramaSize3D saPanoramaSize3D;

class SALiveViewMetering : public ChoiceInterpreter
{
public:
    SALiveViewMetering ()
    {
        choices[0]    = "n/a";
        choices[16]   = "40 segment";
        choices[32]   = "1200-zone Evaluative";
    }
};
SALiveViewMetering saLiveViewMetering;

class SAWhiteBalanceInterpreter: public ChoiceInterpreter
{
public:
    SAWhiteBalanceInterpreter()
    {
        choices[ 0x0] = "Auto";
        choices[ 0x1] = "Color Temperature/Color Filter";
        choices[0x10] = "Daylight";
        choices[0x20] = "Cloudy";
        choices[0x30] = "Shade";
        choices[0x40] = "Tungsten";
        choices[0x50] = "Flash";
        choices[0x60] = "Fluorescent";
        choices[0x70] = "Custom";
        choices[0x80] = "Underwater";
    }
};
SAWhiteBalanceInterpreter saWhiteBalanceInterpreter;

class SAWhiteBalanceSettingInterpreter: public ChoiceInterpreter
{
public:
    SAWhiteBalanceSettingInterpreter()
    {
        choices[0x10] = "Auto (-3)";
        choices[0x11] = "Auto (-2)";
        choices[0x12] = "Auto (-1)";
        choices[0x13] = "Auto (0)";
        choices[0x14] = "Auto (+1)";
        choices[0x15] = "Auto (+2)";
        choices[0x16] = "Auto (+3)";
        choices[0x20] = "Daylight (-3)";
        choices[0x21] = "Daylight (-2)";
        choices[0x22] = "Daylight (-1)";
        choices[0x23] = "Daylight (0)";
        choices[0x24] = "Daylight (+1)";
        choices[0x25] = "Daylight (+2)";
        choices[0x26] = "Daylight (+3)";
        choices[0x30] = "Shade (-3)";
        choices[0x31] = "Shade (-2)";
        choices[0x32] = "Shade (-1)";
        choices[0x33] = "Shade (0)";
        choices[0x34] = "Shade (+1)";
        choices[0x35] = "Shade (+2)";
        choices[0x36] = "Shade (+3)";
        choices[0x40] = "Cloudy (-3)";
        choices[0x41] = "Cloudy (-2)";
        choices[0x42] = "Cloudy (-1)";
        choices[0x43] = "Cloudy (0)";
        choices[0x44] = "Cloudy (+1)";
        choices[0x45] = "Cloudy (+2)";
        choices[0x46] = "Cloudy (+3)";
        choices[0x50] = "Tungsten (-3)";
        choices[0x51] = "Tungsten (-2)";
        choices[0x52] = "Tungsten (-1)";
        choices[0x53] = "Tungsten (0)";
        choices[0x54] = "Tungsten (+1)";
        choices[0x55] = "Tungsten (+2)";
        choices[0x56] = "Tungsten (+3)";
        choices[0x60] = "Fluorescent (-3)";
        choices[0x61] = "Fluorescent (-2)";
        choices[0x62] = "Fluorescent (-1)";
        choices[0x63] = "Fluorescent (0)";
        choices[0x64] = "Fluorescent (+1)";
        choices[0x65] = "Fluorescent (+2)";
        choices[0x66] = "Fluorescent (+3)";
        choices[0x70] = "Flash (-3)";
        choices[0x71] = "Flash (-2)";
        choices[0x72] = "Flash (-1)";
        choices[0x73] = "Flash (0)";
        choices[0x74] = "Flash (+1)";
        choices[0x75] = "Flash (+2)";
        choices[0x76] = "Flash (+3)";
        choices[0xa3] = "Custom";
        choices[0xf3] = "Color Temperature/Color Filter";
    }
};
SAWhiteBalanceSettingInterpreter saWhiteBalanceSettingInterpreter;

class SASceneModeInterpreter : public ChoiceInterpreter
{
public:
    SASceneModeInterpreter ()
    {
        choices[0]  = "Standard";
        choices[1]  = "Portrait";
        choices[2]  = "Text";
        choices[3]  = "Night Scene";
        choices[4]  = "Sunset";
        choices[5]  = "Sports";
        choices[6]  = "Landscape";
        choices[7]  = "Night Portrait";
        choices[8]  = "Macro";
        choices[9]  = "Super Macro";
        choices[16] = "Auto";
        choices[17] = "Night View/Portrait";
        choices[18] = "Sweep Panorama";
        choices[19] = "Handheld Night Shot";
        choices[20] = "Anti Motion Blur";
        choices[21] = "Cont. Priority AE";
        choices[22] = "Auto+";
        choices[23] = "3D Sweep Panorama";
        choices[24] = "Superior Auto";
        choices[25] = "High Sensitivity";
        choices[26] = "Fireworks";
        choices[27] = "Food";
        choices[28] = "Pet";
        choices[33] = "HDR";
        choices[65535] = "n/a";
    }
};
SASceneModeInterpreter saSceneModeInterpreter;

class SAZoneMatchingInterpreter : public ChoiceInterpreter
{
public:
    SAZoneMatchingInterpreter ()
    {
        choices[0] = "ISO Setting Used";
        choices[1] = "High Key";
        choices[2] = "Low Key";
    }
};
SAZoneMatchingInterpreter saZoneMatchingInterpreter;

class SADynamicRangeOptimizerInterpreter : public ChoiceInterpreter
{
public:
    SADynamicRangeOptimizerInterpreter ()
    {
        choices[0] = "Off";
        choices[1] = "Standard";
        choices[2] = "Advanced";
        choices[3] = "Auto";
        choices[8] = "Advanced Lv1";
        choices[9] = "Advanced Lv2";
        choices[10] = "Advanced Lv3";
        choices[11] = "Advanced Lv4";
        choices[12] = "Advanced Lv5";
        choices[16] = "Lv1";
        choices[17] = "Lv2";
        choices[18] = "Lv3";
        choices[19] = "Lv4";
        choices[20] = "Lv5";
    }
};
SADynamicRangeOptimizerInterpreter saDynamicRangeOptimizerInterpreter;

class SAColorModeInterpreter : public ChoiceInterpreter
{
public:
    SAColorModeInterpreter ()
    {
        choices[0]  = "Standard";
        choices[1]  = "Vivid";
        choices[2]  = "Portrait";
        choices[3]  = "Landscape";
        choices[4]  = "Sunset";
        choices[5]  = "Night View/Portrait";
        choices[6]  = "B&W";
        choices[7]  = "Adobe RGB";
        choices[12] = "Neutral";
        choices[13] = "Clear";
        choices[14] = "Deep";
        choices[15] = "Light";
        choices[16] = "Autumn Leaves";
        choices[17] = "Sepia";
        choices[100] = "Neutral";
        choices[101] = "Clear";
        choices[102] = "Deep";
        choices[103] = "Light";
        choices[104] = "Night View";
        choices[105] = "Autumn Leaves";
    }
};
SAColorModeInterpreter saColorModeInterpreter;

class SAExposureModeInterpreter : public ChoiceInterpreter
{
public:
    SAExposureModeInterpreter ()
    {
        choices[0]  = "Program AE";
        choices[1]  = "Portrait";
        choices[2]  = "Beach";
        choices[3]  = "Sports";
        choices[4]  = "Snow";
        choices[5]  = "Landscape";
        choices[6]  = "Auto";
        choices[7]  = "Aperture-priority AE";
        choices[8]  = "Shutter speed priority AE";
        choices[9]  = "Night Scene / Twilight";
        choices[10] = "Hi-Speed Shutter";
        choices[11] = "Twilight Portrait";
        choices[12] = "Soft Snap/Portrait";
        choices[13] = "Fireworks";
        choices[14] = "Smile Shutter";
        choices[15] = "Manual";
        choices[18] = "High Sensitivity";
        choices[19] = "Macro";
        choices[20] = "Advanced Sports Shooting";
        choices[29] = "Underwater";
        choices[33] = "Food";
        choices[34] = "Sweep Panorama";
        choices[35] = "Handheld Night Shot";
        choices[36] = "Anti Motion Blur";
        choices[37] = "Pet";
        choices[38] = "Backlight Correction HDR";
        choices[39] = "Superior Auto";
        choices[40] = "Background Defocus";
        choices[41] = "Soft Skin";
        choices[42] = "3D Image";
        choices[65535] = "n/a";
    }
};
SAExposureModeInterpreter saExposureModeInterpreter;

class SAQualityInterpreter : public ChoiceInterpreter
{
public:
    SAQualityInterpreter ()
    {
        choices[0]  = "Normal";
        choices[1]  = "Fine";
    }
};
SAQualityInterpreter saQualityInterpreter;

class SAAntiBlurInterpreter : public ChoiceInterpreter
{
public:
    SAAntiBlurInterpreter ()
    {
        choices[0]  = "Off";
        choices[1]  = "On (Continuous)";
        choices[2]  = "On (Shooting)";
        choices[65535]  = "n/a";
    }
};
SAAntiBlurInterpreter saAntiBlurInterpreter;

class SALensIDInterpreter : public IntLensInterpreter<int>
{
public:
    SALensIDInterpreter ()
    {
        choices = {
            {0, "Minolta AF 28-85mm f/3.5-4.5 New"},
            {1, "Minolta AF 80-200mm f/2.8 HS-APO G"},
            {2, "Minolta AF 28-70mm f/2.8 G"},
            {3, "Minolta AF 28-80mm f/4-5.6"},
            {4, "Minolta AF 85mm f/1.4G"},
            {5, "Minolta AF 35-70mm f/3.5-4.5 [II]"},
            {6, "Minolta AF 24-85mm f/3.5-4.5 [New]"},
            {7, "Minolta AF 100-300mm f/4.5-5.6 APO [New] or 100-400mm or Sigma Lens"},
            {7, "Minolta AF 100-400mm f/4.5-6.7 APO"},
            {7, "Sigma AF 100-300mm f/4 EX DG IF"},
            {8, "Minolta AF 70-210mm f/4.5-5.6 [II]"},
            {9, "Minolta AF 50mm f/3.5 Macro"},
            {10, "Minolta AF 28-105mm f/3.5-4.5 [New]"},
            {11, "Minolta AF 300mm f/4 HS-APO G"},
            {12, "Minolta AF 100mm f/2.8 Soft Focus"},
            {13, "Minolta AF 75-300mm f/4.5-5.6 (New or II)"},
            {14, "Minolta AF 100-400mm f/4.5-6.7 APO"},
            {15, "Minolta AF 400mm f/4.5 HS-APO G"},
            {16, "Minolta AF 17-35mm f/3.5 G"},
            {17, "Minolta AF 20-35mm f/3.5-4.5"},
            {18, "Minolta AF 28-80mm f/3.5-5.6 II"},
            {19, "Minolta AF 35mm f/1.4 G"},
            {20, "Minolta/Sony 135mm f/2.8 [T4.5] STF"},
            {22, "Minolta AF 35-80mm f/4-5.6 II"},
            {23, "Minolta AF 200mm f/4 Macro APO G"},
            {24, "Minolta/Sony AF 24-105mm f/3.5-4.5 (D) or Sigma or Tamron Lens"},
            {24, "Sigma 18-50mm f/2.8"},
            {24, "Sigma 17-70mm f/2.8-4.5 (D)"},
            {24, "Sigma 20-40mm f/2.8 EX DG Aspherical IF"},
            {24, "Sigma 18-200mm f/3.5-6.3 DC"},
            {24, "Sigma DC 18-125mm f/4-5,6 D"},
            {24, "Tamron SP AF 28-75mm f/2.8 XR Di LD Aspherical [IF] Macro"},
            {25, "Minolta AF 100-300mm f/4.5-5.6 APO (D) or Sigma Lens"},
            {25, "Sigma 100-300mm f/4 EX (APO (D) or D IF)"},
            {25, "Sigma 70mm f/2.8 EX DG Macro"},
            {25, "Sigma 20mm f/1.8 EX DG Aspherical RF"},
            {25, "Sigma 30mm f/1.4 EX DC"},
            {25, "Sigma 24mm f/1.8 EX DG ASP Macro"},
            {27, "Minolta AF 85mm f/1.4 G (D)"},
            {28, "Minolta/Sony AF 100mm f/2.8 Macro (D) or Tamron Lens"},
            {28, "Tamron SP AF 90mm f/2.8 Di Macro"},
            {28, "Tamron SP AF 180mm f/3.5 Di LD [IF] Macro"},
            {29, "Minolta/Sony AF 75-300mm f/4.5-5.6 (D)"},
            {30, "Minolta AF 28-80mm f/3.5-5.6 (D) or Sigma Lens"},
            {30, "Sigma AF 10-20mm f/4-5.6 EX DC"},
            {30, "Sigma AF 12-24mm f/4.5-5.6 EX DG"},
            {30, "Sigma 28-70mm EX DG f/2.8"},
            {30, "Sigma 55-200mm f/4-5.6 DC"},
            {31, "Minolta/Sony AF 50mm f/2.8 Macro (D) or f/3.5"},
            {31, "Minolta/Sony AF 50mm f/3.5 Macro"},
            {32, "Minolta/Sony AF 300mm f/2.8 G or 1.5x Teleconverter"},
            {33, "Minolta/Sony AF 70-200mm f/2.8 G"},
            {35, "Minolta AF 85mm f/1.4 G (D) Limited"},
            {36, "Minolta AF 28-100mm f/3.5-5.6 (D)"},
            {38, "Minolta AF 17-35mm f/2.8-4 (D)"},
            {39, "Minolta AF 28-75mm f/2.8 (D)"},
            {40, "Minolta/Sony AF DT 18-70mm f/3.5-5.6 (D)"},
            {41, "Minolta/Sony AF DT 11-18mm f/4.5-5.6 (D) or Tamron Lens"},
            {41, "Tamron SP AF 11-18mm f/4.5-5.6 Di II LD Aspherical IF"},
            {42, "Minolta/Sony AF DT 18-200mm f/3.5-6.3 (D)"},
            {43, "Sony 35mm f/1.4 G (SAL35F14G)"},
            {44, "Sony 50mm f/1.4 (SAL50F14)"},
            {45, "Carl Zeiss Planar T* 85mm f/1.4 ZA (SAL85F14Z)"},
            {46, "Carl Zeiss Vario-Sonnar T* DT 16-80mm f/3.5-4.5 ZA (SAL1680Z)"},
            {47, "Carl Zeiss Sonnar T* 135mm f/1.8 ZA (SAL135F18Z)"},
            {48, "Carl Zeiss Vario-Sonnar T* 24-70mm f/2.8 ZA SSM (SAL2470Z) or ZA SSM II"},
            {48, "Carl Zeiss Vario-Sonnar T* 24-70mm f/2.8 ZA SSM II (SAL2470Z2)"},
            {49, "Sony DT 55-200mm f/4-5.6 (SAL55200)"},
            {50, "Sony DT 18-250mm f/3.5-6.3 (SAL18250)"},
            {51, "Sony DT 16-105mm f/3.5-5.6 (SAL16105)"},
            {52, "Sony 70-300mm f/4.5-5.6 G SSM (SAL70300G) or G SSM II or Tamron Lens"},
            {52, "Sony 70-300mm f/4.5-5.6 G SSM II (SAL70300G2)"},
            {52, "Tamron SP 70-300mm f/4-5.6 Di USD"},
            {53, "Sony 70-400mm f/4-5.6 G SSM (SAL70400G)"},
            {54, "Carl Zeiss Vario-Sonnar T* 16-35mm f/2.8 ZA SSM (SAL1635Z) or ZA SSM II"},
            {54, "Carl Zeiss Vario-Sonnar T* 16-35mm f/2.8 ZA SSM II (SAL1635Z2)"},
            {55, "Sony DT 18-55mm f/3.5-5.6 SAM (SAL1855) or SAM II"},
            {55, "Sony DT 18-55mm f/3.5-5.6 SAM II (SAL18552)"},
            {56, "Sony DT 55-200mm f/4-5.6 SAM (SAL55200-2)"},
            {57, "Sony DT 50mm f/1.8 SAM (SAL50F18) or Tamron Lens or Commlite CM-EF-NEX adapter"},
            {57, "Tamron SP AF 60mm f/2 Di II LD [IF] Macro 1:1"},
            {57, "Tamron 18-270mm f/3.5-6.3 Di II PZD"},
            {58, "Sony DT 30mm f/2.8 Macro SAM (SAL30M28)"},
            {59, "Sony 28-75mm f/2.8 SAM (SAL2875)"},
            {60, "Carl Zeiss Distagon T* 24mm f/2 ZA SSM (SAL24F20Z)"},
            {61, "Sony 85mm f/2.8 SAM (SAL85F28)"},
            {62, "Sony DT 35mm f/1.8 SAM (SAL35F18)"},
            {63, "Sony DT 16-50mm f/2.8 SSM (SAL1650)"},
            {64, "Sony 500mm f/4 G SSM (SAL500F40G)"},
            {65, "Sony DT 18-135mm f/3.5-5.6 SAM (SAL18135)"},
            {66, "Sony 300mm f/2.8 G SSM II (SAL300F28G2)"},
            {67, "Sony 70-200mm f/2.8 G SSM II (SAL70200G2)"},
            {68, "Sony DT 55-300mm f/4.5-5.6 SAM (SAL55300)"},
            {69, "Sony 70-400mm f/4-5.6 G SSM II (SAL70400G2)"},
            {70, "Carl Zeiss Planar T* 50mm f/1.4 ZA SSM (SAL50F14Z)"},
            {128, "Tamron or Sigma Lens (128)"},
            {128, "Tamron AF 18-200mm f/3.5-6.3 XR Di II LD Aspherical [IF] Macro"},
            {128, "Tamron AF 28-300mm f/3.5-6.3 XR Di LD Aspherical [IF] Macro"},
            {128, "Tamron 80-300mm f/3.5-6.3"},
            {128, "Tamron AF 28-200mm f/3.8-5.6 XR Di Aspherical [IF] Macro"},
            {128, "Tamron SP AF 17-35mm f/2.8-4 Di LD Aspherical IF"},
            {128, "Sigma AF 50-150mm f/2.8 EX DC APO HSM II"},
            {128, "Sigma 10-20mm f/3.5 EX DC HSM"},
            {128, "Sigma 70-200mm f/2.8 II EX DG APO MACRO HSM"},
            {128, "Sigma 10mm f/2.8 EX DC HSM Fisheye"},
            {128, "Sigma 50mm f/1.4 EX DG HSM"},
            {128, "Sigma 85mm f/1.4 EX DG HSM"},
            {128, "Sigma 24-70mm f/2.8 IF EX DG HSM"},
            {128, "Sigma 18-250mm f/3.5-6.3 DC OS HSM"},
            {128, "Sigma 17-50mm f/2.8 EX DC HSM"},
            {128, "Sigma 17-70mm f/2.8-4 DC Macro HSM"},
            {128, "Sigma 150mm f/2.8 EX DG OS HSM APO Macro"},
            {128, "Sigma 150-500mm f/5-6.3 APO DG OS HSM"},
            {128, "Tamron AF 28-105mm f/4-5.6 [IF]"},
            {128, "Sigma 35mm f/1.4 DG HSM"},
            {128, "Sigma 18-35mm f/1.8 DC HSM"},
            {128, "Sigma 50-500mm f/4.5-6.3 APO DG OS HSM"},
            {128, "Sigma 24-105mm f/4 DG HSM | Art 013"},
            {129, "Tamron Lens (129)"},
            {129, "Tamron 200-400mm f/5.6 LD"},
            {129, "Tamron 70-300mm f/4-5.6 LD"},
            {131, "Tamron 20-40mm f/2.7-3.5 SP Aspherical IF"},
            {135, "Vivitar 28-210mm f/3.5-5.6"},
            {136, "Tokina EMZ M100 AF 100mm f/3.5"},
            {137, "Cosina 70-210mm f/2.8-4 AF"},
            {138, "Soligor 19-35mm f/3.5-4.5"},
            {139, "Tokina AF 28-300mm f/4-6.3"},
            {142, "Voigtlander 70-300mm f/4.5-5.6"},
            {146, "Voigtlander Macro APO-Lanthar 125mm f/2.5 SL"},
            {194, "Tamron SP AF 17-50mm f/2.8 XR Di II LD Aspherical [IF]"},
            {202, "Tamron SP AF 70-200mm f/2.8 Di LD [IF] Macro"},
            {203, "Tamron SP 70-200mm f/2.8 Di USD"},
            {204, "Tamron SP 24-70mm f/2.8 Di USD"},
            {212, "Tamron 28-300mm f/3.5-6.3 Di PZD"},
            {213, "Tamron 16-300mm f/3.5-6.3 Di II PZD Macro"},
            {214, "Tamron SP 150-600mm f/5-6.3 Di USD"},
            {215, "Tamron SP 15-30mm f/2.8 Di USD"},
            {218, "Tamron SP 90mm f/2.8 Di Macro 1:1 USD (F017)"},
            {224, "Tamron SP 90mm f/2.8 Di Macro 1:1 USD (F004)"},
            {255, "Tamron Lens (255)"},
            {255, "Tamron SP AF 17-50mm f/2.8 XR Di II LD Aspherical"},
            {255, "Tamron AF 18-250mm f/3.5-6.3 XR Di II LD"},
            {255, "Tamron AF 55-200mm f/4-5.6 Di II LD Macro"},
            {255, "Tamron AF 70-300mm f/4-5.6 Di LD Macro 1:2"},
            {255, "Tamron SP AF 200-500mm f/5.0-6.3 Di LD IF"},
            {255, "Tamron SP AF 10-24mm f/3.5-4.5 Di II LD Aspherical IF"},
            {255, "Tamron SP AF 70-200mm f/2.8 Di LD IF Macro"},
            {255, "Tamron SP AF 28-75mm f/2.8 XR Di LD Aspherical IF"},
            {255, "Tamron AF 90-300mm f/4.5-5.6 Telemacro"},
            {1868, "Sigma MC-11 SA-E Mount Converter with not-supported Sigma lens"},
            {2550, "Minolta AF 50mm f/1.7"},
            {2551, "Minolta AF 35-70mm f/4 or Other Lens"},
            {2551, "Sigma UC AF 28-70mm f/3.5-4.5"},
            {2551, "Sigma AF 28-70mm f/2.8"},
            {2551, "Sigma M-AF 70-200mm f/2.8 EX Aspherical"},
            {2551, "Quantaray M-AF 35-80mm f/4-5.6"},
            {2551, "Tokina 28-70mm f/2.8-4.5 AF"},
            {2552, "Minolta AF 28-85mm f/3.5-4.5 or Other Lens"},
            {2552, "Tokina 19-35mm f/3.5-4.5"},
            {2552, "Tokina 28-70mm f/2.8 AT-X"},
            {2552, "Tokina 80-400mm f/4.5-5.6 AT-X AF II 840"},
            {2552, "Tokina AF PRO 28-80mm f/2.8 AT-X 280"},
            {2552, "Tokina AT-X PRO [II] AF 28-70mm f/2.6-2.8 270"},
            {2552, "Tamron AF 19-35mm f/3.5-4.5"},
            {2552, "Angenieux AF 28-70mm f/2.6"},
            {2552, "Tokina AT-X 17 AF 17mm f/3.5"},
            {2552, "Tokina 20-35mm f/3.5-4.5 II AF"},
            {2553, "Minolta AF 28-135mm f/4-4.5 or Sigma Lens"},
            {2553, "Sigma ZOOM-alpha 35-135mm f/3.5-4.5"},
            {2553, "Sigma 28-105mm f/2.8-4 Aspherical"},
            {2553, "Sigma 28-105mm f/4-5.6 UC"},
            {2554, "Minolta AF 35-105mm f/3.5-4.5"},
            {2555, "Minolta AF 70-210mm f/4 Macro or Sigma Lens"},
            {2555, "Sigma 70-210mm f/4-5.6 APO"},
            {2555, "Sigma M-AF 70-200mm f/2.8 EX APO"},
            {2555, "Sigma 75-200mm f/2.8-3.5"},
            {2556, "Minolta AF 135mm f/2.8"},
            {2557, "Minolta/Sony AF 28mm f/2.8"},
            {2558, "Minolta AF 24-50mm f/4"},
            {2560, "Minolta AF 100-200mm f/4.5"},
            {2561, "Minolta AF 75-300mm f/4.5-5.6 or Sigma Lens"},
            {2561, "Sigma 70-300mm f/4-5.6 DL Macro"},
            {2561, "Sigma 300mm f/4 APO Macro"},
            {2561, "Sigma AF 500mm f/4.5 APO"},
            {2561, "Sigma AF 170-500mm f/5-6.3 APO Aspherical"},
            {2561, "Tokina AT-X AF 300mm f/4"},
            {2561, "Tokina AT-X AF 400mm f/5.6 SD"},
            {2561, "Tokina AF 730 II 75-300mm f/4.5-5.6"},
            {2561, "Sigma 800mm f/5.6 APO"},
            {2561, "Sigma AF 400mm f/5.6 APO Macro"},
            {2562, "Minolta AF 50mm f/1.4 [New]"},
            {2563, "Minolta AF 300mm f/2.8 APO or Sigma Lens"},
            {2563, "Sigma AF 50-500mm f/4-6.3 EX DG APO"},
            {2563, "Sigma AF 170-500mm f/5-6.3 APO Aspherical"},
            {2563, "Sigma AF 500mm f/4.5 EX DG APO"},
            {2563, "Sigma 400mm f/5.6 APO"},
            {2564, "Minolta AF 50mm f/2.8 Macro or Sigma Lens"},
            {2564, "Sigma 50mm f/2.8 EX Macro"},
            {2565, "Minolta AF 600mm f/4 APO"},
            {2566, "Minolta AF 24mm f/2.8 or Sigma Lens"},
            {2566, "Sigma 17-35mm f/2.8-4 EX Aspherical"},
            {2572, "Minolta/Sony AF 500mm f/8 Reflex"},
            {2578, "Minolta/Sony AF 16mm f/2.8 Fisheye or Sigma Lens"},
            {2578, "Sigma 8mm f/4 EX [DG] Fisheye"},
            {2578, "Sigma 14mm f/3.5"},
            {2578, "Sigma 15mm f/2.8 Fisheye"},
            {2579, "Minolta/Sony AF 20mm f/2.8 or Tokina Lens"},
            {2579, "Tokina AT-X Pro DX 11-16mm f/2.8"},
            {2581, "Minolta AF 100mm f/2.8 Macro [New] or Sigma or Tamron Lens"},
            {2581, "Sigma AF 90mm f/2.8 Macro"},
            {2581, "Sigma AF 105mm f/2.8 EX [DG] Macro"},
            {2581, "Sigma 180mm f/5.6 Macro"},
            {2581, "Sigma 180mm f/3.5 EX DG Macro"},
            {2581, "Tamron 90mm f/2.8 Macro"},
            {2585, "Minolta AF 35-105mm f/3.5-4.5 New or Tamron Lens"},
            {2585, "Beroflex 35-135mm f/3.5-4.5"},
            {2585, "Tamron 24-135mm f/3.5-5.6"},
            {2588, "Minolta AF 70-210mm f/3.5-4.5"},
            {2589, "Minolta AF 80-200mm f/2.8 APO or Tokina Lens"},
            {2589, "Tokina 80-200mm f/2.8"},
            {2590, "Minolta AF 200mm f/2.8 G APO + Minolta AF 1.4x APO or Other Lens + 1.4x"},
            {2590, "Minolta AF 600mm f/4 HS-APO G + Minolta AF 1.4x APO"},
            {2591, "Minolta AF 35mm f/1.4"},
            {2592, "Minolta AF 85mm f/1.4 G (D)"},
            {2593, "Minolta AF 200mm f/2.8 APO"},
            {2594, "Minolta AF 3x-1x f/1.7-2.8 Macro"},
            {2596, "Minolta AF 28mm f/2"},
            {2597, "Minolta AF 35mm f/2 [New]"},
            {2598, "Minolta AF 100mm f/2"},
            {2601, "Minolta AF 200mm f/2.8 G APO + Minolta AF 2x APO or Other Lens + 2x"},
            {2601, "Minolta AF 600mm f/4 HS-APO G + Minolta AF 2x APO"},
            {2604, "Minolta AF 80-200mm f/4.5-5.6"},
            {2605, "Minolta AF 35-80mm f/4-5.6"},
            {2606, "Minolta AF 100-300mm f/4.5-5.6"},
            {2607, "Minolta AF 35-80mm f/4-5.6"},
            {2608, "Minolta AF 300mm f/2.8 HS-APO G"},
            {2609, "Minolta AF 600mm f/4 HS-APO G"},
            {2612, "Minolta AF 200mm f/2.8 HS-APO G"},
            {2613, "Minolta AF 50mm f/1.7 New"},
            {2615, "Minolta AF 28-105mm f/3.5-4.5 xi"},
            {2616, "Minolta AF 35-200mm f/4.5-5.6 xi"},
            {2618, "Minolta AF 28-80mm f/4-5.6 xi"},
            {2619, "Minolta AF 80-200mm f/4.5-5.6 xi"},
            {2620, "Minolta AF 28-70mm f/2.8 G"},
            {2621, "Minolta AF 100-300mm f/4.5-5.6 xi"},
            {2624, "Minolta AF 35-80mm f/4-5.6 Power Zoom"},
            {2628, "Minolta AF 80-200mm f/2.8 HS-APO G"},
            {2629, "Minolta AF 85mm f/1.4 New"},
            {2631, "Minolta/Sony AF 100-300mm f/4.5-5.6 APO"},
            {2632, "Minolta AF 24-50mm f/4 New"},
            {2638, "Minolta AF 50mm f/2.8 Macro New"},
            {2639, "Minolta AF 100mm f/2.8 Macro"},
            {2641, "Minolta/Sony AF 20mm f/2.8 New"},
            {2642, "Minolta AF 24mm f/2.8 New"},
            {2644, "Minolta AF 100-400mm f/4.5-6.7 APO"},
            {2662, "Minolta AF 50mm f/1.4 New"},
            {2667, "Minolta AF 35mm f/2 New"},
            {2668, "Minolta AF 28mm f/2 New"},
            {2672, "Minolta AF 24-105mm f/3.5-4.5 (D)"},
            {3046, "Metabones Canon EF Speed Booster"},
            {4567, "Tokina 70-210mm f/4-5.6"},
            {4571, "Vivitar 70-210mm f/4.5-5.6"},
            {4574, "2x Teleconverter or Tamron or Tokina Lens"},
            {4574, "Tamron SP AF 90mm f/2.5"},
            {4574, "Tokina RF 500mm f/8.0 x2"},
            {4574, "Tokina 300mm f/2.8 x2"},
            {4575, "1.4x Teleconverter"},
            {4585, "Tamron SP AF 300mm f/2.8 LD IF"},
            {4586, "Tamron SP AF 35-105mm f/2.8 LD Aspherical IF"},
            {4587, "Tamron AF 70-210mm f/2.8 SP LD"},
            {4812, "Metabones Canon EF Speed Booster Ultra"},
            {6118, "Canon EF Adapter"},
            {6528, "Sigma 16mm f/2.8 Filtermatic Fisheye"},
            {6553, "E-Mount, T-Mount, Other Lens or no lens"},
            {6553, "Sony E 16mm f/2.8"},
            {6553, "Sony E 18-55mm f/3.5-5.6 OSS"},
            {6553, "Sony E 55-210mm f/4.5-6.3 OSS"},
            {6553, "Sony E 18-200mm f/3.5-6.3 OSS"},
            {6553, "Sony E 30mm f/3.5 Macro"},
            {6553, "Sony E 24mm f/1.8 ZA"},
            {6553, "Sony E 50mm f/1.8 OSS"},
            {6553, "Sony E 16-70mm f/4 ZA OSS"},
            {6553, "Sony E 10-18mm f/4 OSS"},
            {6553, "Sony E PZ 16-50mm f/3.5-5.6 OSS"},
            {6553, "Sony FE 35mm f/2.8 ZA"},
            {6553, "Sony FE 24-70mm f/4 ZA OSS"},
            {6553, "Sony FE 85mm f/1.8"},
            {6553, "Sony E 18-200mm f/3.5-6.3 OSS LE"},
            {6553, "Sony E 20mm f/2.8"},
            {6553, "Sony E 35mm f/1.8 OSS"},
            {6553, "Sony E PZ 18-105mm f/4 G OSS"},
            {6553, "Sony FE 90mm f/2.8 Macro G OSS"},
            {6553, "Sony E 18-50mm f/4-5.6"},
            {6553, "Sony E PZ 18-200mm f/3.5-6.3 OSS"},
            {6553, "Sony FE 55mm f/1.8 ZA"},
            {6553, "Sony FE 70-200mm f/4 G OSS"},
            {6553, "Sony FE 16-35mm f/4 ZA OSS"},
            {6553, "Sony FE 50mm f/2.8 Macro"},
            {6553, "Sony FE 28-70mm f/3.5-5.6 OSS"},
            {6553, "Sony FE 35mm f/1.4 ZA"},
            {6553, "Sony FE 24-240mm f/3.5-6.3 OSS"},
            {6553, "Sony FE 28mm f/2"},
            {6553, "Sony FE PZ 28-135mm f/4 G OSS"},
            {6553, "Sony FE 100mm f/2.8 STF GM OSS"},
            {6553, "Sony FE 24-70mm f/2.8 GM"},
            {6553, "Sony FE 50mm f/1.4 ZA"},
            {6553, "Sony FE 85mm f/1.4 GM"},
            {6553, "Sony FE 50mm f/1.8"},
            {6553, "Sony FE 21mm f/2.8 (SEL28F20 + SEL075UWC)"},
            {6553, "Sony FE 16mm f/3.5 Fisheye (SEL28F20 + SEL057FEC)"},
            {6553, "Sony FE 70-300mm f/4.5-5.6 G OSS"},
            {6553, "Sony FE 70-200mm f/2.8 GM OSS"},
            {6553, "Sony FE 70-200mm f/2.8 GM OSS + 1.4X Teleconverter"},
            {6553, "Sony FE 70-200mm f/2.8 GM OSS + 2X Teleconverter"},
            {6553, "Samyang AF 50mm f/1.4 FE"},
            {6553, "Samyang AF 14mm f/2.8 FE"},
            {6553, "Sigma 19mm f/2.8 [EX] DN"},
            {6553, "Sigma 30mm f/2.8 [EX] DN"},
            {6553, "Sigma 60mm f/2.8 DN"},
            {6553, "Sigma 30mm f/1.4 DC DN | C 016"},
            {6553, "Tamron 18-200mm f/3.5-6.3 Di III VC"},
            {6553, "Tokina Firin 20mm f/2 FE MF"},
            {6553, "Zeiss Touit 12mm f/2.8"},
            {6553, "Zeiss Touit 32mm f/1.8"},
            {6553, "Zeiss Touit 50mm f/2.8 Macro"},
            {6553, "Zeiss Batis 25mm f/2"},
            {6553, "Zeiss Batis 85mm f/1.8"},
            {6553, "Zeiss Batis 18mm f/2.8"},
            {6553, "Zeiss Batis 135mm f/2.8"},
            {6553, "Zeiss Loxia 50mm f/2"},
            {6553, "Zeiss Loxia 35mm f/2"},
            {6553, "Zeiss Loxia 21mm f/2.8"},
            {6553, "Zeiss Loxia 85mm f/2.4"},
            {6553, "Arax MC 35mm f/2.8 Tilt+Shift"},
            {6553, "Arax MC 80mm f/2.8 Tilt+Shift"},
            {6553, "Zenitar MF 16mm f/2.8 Fisheye M42"},
            {6553, "Samyang 500mm Mirror f/8.0"},
            {6553, "Pentacon Auto 135mm f/2.8"},
            {6553, "Pentacon Auto 29mm f/2.8"},
            {6553, "Helios 44-2 58mm f/2.0"},
            {18688, "Sigma MC-11 SA-E Mount Converter with not-supported Sigma lens"},
            {25501, "Minolta AF 50mm f/1.7"},
            {25511, "Minolta AF 35-70mm f/4 or Other Lens"},
            {25511, "Sigma UC AF 28-70mm f/3.5-4.5"},
            {25511, "Sigma AF 28-70mm f/2.8"},
            {25511, "Sigma M-AF 70-200mm f/2.8 EX Aspherical"},
            {25511, "Quantaray M-AF 35-80mm f/4-5.6"},
            {25511, "Tokina 28-70mm f/2.8-4.5 AF"},
            {25521, "Minolta AF 28-85mm f/3.5-4.5 or Other Lens"},
            {25521, "Tokina 19-35mm f/3.5-4.5"},
            {25521, "Tokina 28-70mm f/2.8 AT-X"},
            {25521, "Tokina 80-400mm f/4.5-5.6 AT-X AF II 840"},
            {25521, "Tokina AF PRO 28-80mm f/2.8 AT-X 280"},
            {25521, "Tokina AT-X PRO [II] AF 28-70mm f/2.6-2.8 270"},
            {25521, "Tamron AF 19-35mm f/3.5-4.5"},
            {25521, "Angenieux AF 28-70mm f/2.6"},
            {25521, "Tokina AT-X 17 AF 17mm f/3.5"},
            {25521, "Tokina 20-35mm f/3.5-4.5 II AF"},
            {25531, "Minolta AF 28-135mm f/4-4.5 or Sigma Lens"},
            {25531, "Sigma ZOOM-alpha 35-135mm f/3.5-4.5"},
            {25531, "Sigma 28-105mm f/2.8-4 Aspherical"},
            {25531, "Sigma 28-105mm f/4-5.6 UC"},
            {25541, "Minolta AF 35-105mm f/3.5-4.5"},
            {25551, "Minolta AF 70-210mm f/4 Macro or Sigma Lens"},
            {25551, "Sigma 70-210mm f/4-5.6 APO"},
            {25551, "Sigma M-AF 70-200mm f/2.8 EX APO"},
            {25551, "Sigma 75-200mm f/2.8-3.5"},
            {25561, "Minolta AF 135mm f/2.8"},
            {25571, "Minolta/Sony AF 28mm f/2.8"},
            {25581, "Minolta AF 24-50mm f/4"},
            {25601, "Minolta AF 100-200mm f/4.5"},
            {25611, "Minolta AF 75-300mm f/4.5-5.6 or Sigma Lens"},
            {25611, "Sigma 70-300mm f/4-5.6 DL Macro"},
            {25611, "Sigma 300mm f/4 APO Macro"},
            {25611, "Sigma AF 500mm f/4.5 APO"},
            {25611, "Sigma AF 170-500mm f/5-6.3 APO Aspherical"},
            {25611, "Tokina AT-X AF 300mm f/4"},
            {25611, "Tokina AT-X AF 400mm f/5.6 SD"},
            {25611, "Tokina AF 730 II 75-300mm f/4.5-5.6"},
            {25611, "Sigma 800mm f/5.6 APO"},
            {25611, "Sigma AF 400mm f/5.6 APO Macro"},
            {25621, "Minolta AF 50mm f/1.4 [New]"},
            {25631, "Minolta AF 300mm f/2.8 APO or Sigma Lens"},
            {25631, "Sigma AF 50-500mm f/4-6.3 EX DG APO"},
            {25631, "Sigma AF 170-500mm f/5-6.3 APO Aspherical"},
            {25631, "Sigma AF 500mm f/4.5 EX DG APO"},
            {25631, "Sigma 400mm f/5.6 APO"},
            {25641, "Minolta AF 50mm f/2.8 Macro or Sigma Lens"},
            {25641, "Sigma 50mm f/2.8 EX Macro"},
            {25651, "Minolta AF 600mm f/4 APO"},
            {25661, "Minolta AF 24mm f/2.8 or Sigma Lens"},
            {25661, "Sigma 17-35mm f/2.8-4 EX Aspherical"},
            {25721, "Minolta/Sony AF 500mm f/8 Reflex"},
            {25781, "Minolta/Sony AF 16mm f/2.8 Fisheye or Sigma Lens"},
            {25781, "Sigma 8mm f/4 EX [DG] Fisheye"},
            {25781, "Sigma 14mm f/3.5"},
            {25781, "Sigma 15mm f/2.8 Fisheye"},
            {25791, "Minolta/Sony AF 20mm f/2.8 or Tokina Lens"},
            {25791, "Tokina AT-X Pro DX 11-16mm f/2.8"},
            {25811, "Minolta AF 100mm f/2.8 Macro [New] or Sigma or Tamron Lens"},
            {25811, "Sigma AF 90mm f/2.8 Macro"},
            {25811, "Sigma AF 105mm f/2.8 EX [DG] Macro"},
            {25811, "Sigma 180mm f/5.6 Macro"},
            {25811, "Sigma 180mm f/3.5 EX DG Macro"},
            {25811, "Tamron 90mm f/2.8 Macro"},
            {25851, "Beroflex 35-135mm f/3.5-4.5"},
            {25858, "Minolta AF 35-105mm f/3.5-4.5 New or Tamron Lens"},
            {25858, "Tamron 24-135mm f/3.5-5.6"},
            {25881, "Minolta AF 70-210mm f/3.5-4.5"},
            {25891, "Minolta AF 80-200mm f/2.8 APO or Tokina Lens"},
            {25891, "Tokina 80-200mm f/2.8"},
            {25901, "Minolta AF 200mm f/2.8 G APO + Minolta AF 1.4x APO or Other Lens + 1.4x"},
            {25901, "Minolta AF 600mm f/4 HS-APO G + Minolta AF 1.4x APO"},
            {25911, "Minolta AF 35mm f/1.4"},
            {25921, "Minolta AF 85mm f/1.4 G (D)"},
            {25931, "Minolta AF 200mm f/2.8 APO"},
            {25941, "Minolta AF 3x-1x f/1.7-2.8 Macro"},
            {25961, "Minolta AF 28mm f/2"},
            {25971, "Minolta AF 35mm f/2 [New]"},
            {25981, "Minolta AF 100mm f/2"},
            {26011, "Minolta AF 200mm f/2.8 G APO + Minolta AF 2x APO or Other Lens + 2x"},
            {26011, "Minolta AF 600mm f/4 HS-APO G + Minolta AF 2x APO"},
            {26041, "Minolta AF 80-200mm f/4.5-5.6"},
            {26051, "Minolta AF 35-80mm f/4-5.6"},
            {26061, "Minolta AF 100-300mm f/4.5-5.6"},
            {26071, "Minolta AF 35-80mm f/4-5.6"},
            {26081, "Minolta AF 300mm f/2.8 HS-APO G"},
            {26091, "Minolta AF 600mm f/4 HS-APO G"},
            {26121, "Minolta AF 200mm f/2.8 HS-APO G"},
            {26131, "Minolta AF 50mm f/1.7 New"},
            {26151, "Minolta AF 28-105mm f/3.5-4.5 xi"},
            {26161, "Minolta AF 35-200mm f/4.5-5.6 xi"},
            {26181, "Minolta AF 28-80mm f/4-5.6 xi"},
            {26191, "Minolta AF 80-200mm f/4.5-5.6 xi"},
            {26201, "Minolta AF 28-70mm f/2.8 G"},
            {26211, "Minolta AF 100-300mm f/4.5-5.6 xi"},
            {26241, "Minolta AF 35-80mm f/4-5.6 Power Zoom"},
            {26281, "Minolta AF 80-200mm f/2.8 HS-APO G"},
            {26291, "Minolta AF 85mm f/1.4 New"},
            {26311, "Minolta/Sony AF 100-300mm f/4.5-5.6 APO"},
            {26321, "Minolta AF 24-50mm f/4 New"},
            {26381, "Minolta AF 50mm f/2.8 Macro New"},
            {26391, "Minolta AF 100mm f/2.8 Macro"},
            {26411, "Minolta/Sony AF 20mm f/2.8 New"},
            {26421, "Minolta AF 24mm f/2.8 New"},
            {26441, "Minolta AF 100-400mm f/4.5-6.7 APO"},
            {26621, "Minolta AF 50mm f/1.4 New"},
            {26671, "Minolta AF 35mm f/2 New"},
            {26681, "Minolta AF 28mm f/2 New"},
            {26721, "Minolta AF 24-105mm f/3.5-4.5 (D)"},
            {30464, "Metabones Canon EF Speed Booster"},
            {45671, "Tokina 70-210mm f/4-5.6"},
            {45711, "Vivitar 70-210mm f/4.5-5.6"},
            {45741, "2x Teleconverter or Tamron or Tokina Lens"},
            {45741, "Tamron SP AF 90mm f/2.5"},
            {45741, "Tokina RF 500mm f/8.0 x2"},
            {45741, "Tokina 300mm f/2.8 x2"},
            {45751, "1.4x Teleconverter"},
            {45851, "Tamron SP AF 300mm f/2.8 LD IF"},
            {45861, "Tamron SP AF 35-105mm f/2.8 LD Aspherical IF"},
            {45871, "Tamron AF 70-210mm f/2.8 SP LD"},
            {48128, "Metabones Canon EF Speed Booster Ultra"},
            {61184, "Canon EF Adapter"},
            {65280, "Sigma 16mm f/2.8 Filtermatic Fisheye"},
            {65535, "E-Mount, T-Mount, Other Lens or no lens"},
            {65535, "Sony E 16mm f/2.8"},
            {65535, "Sony E 18-55mm f/3.5-5.6 OSS"},
            {65535, "Sony E 55-210mm f/4.5-6.3 OSS"},
            {65535, "Sony E 18-200mm f/3.5-6.3 OSS"},
            {65535, "Sony E 30mm f/3.5 Macro"},
            {65535, "Sony E 24mm f/1.8 ZA"},
            {65535, "Sony E 50mm f/1.8 OSS"},
            {65535, "Sony E 16-70mm f/4 ZA OSS"},
            {65535, "Sony E 10-18mm f/4 OSS"},
            {65535, "Sony E PZ 16-50mm f/3.5-5.6 OSS"},
            {65535, "Sony FE 35mm f/2.8 ZA"},
            {65535, "Sony FE 24-70mm f/4 ZA OSS"},
            {65535, "Sony FE 85mm f/1.8"},
            {65535, "Sony E 18-200mm f/3.5-6.3 OSS LE"},
            {65535, "Sony E 20mm f/2.8"},
            {65535, "Sony E 35mm f/1.8 OSS"},
            {65535, "Sony E PZ 18-105mm f/4 G OSS"},
            {65535, "Sony FE 90mm f/2.8 Macro G OSS"},
            {65535, "Sony E 18-50mm f/4-5.6"},
            {65535, "Sony E PZ 18-200mm f/3.5-6.3 OSS"},
            {65535, "Sony FE 55mm f/1.8 ZA"},
            {65535, "Sony FE 70-200mm f/4 G OSS"},
            {65535, "Sony FE 16-35mm f/4 ZA OSS"},
            {65535, "Sony FE 50mm f/2.8 Macro"},
            {65535, "Sony FE 28-70mm f/3.5-5.6 OSS"},
            {65535, "Sony FE 35mm f/1.4 ZA"},
            {65535, "Sony FE 24-240mm f/3.5-6.3 OSS"},
            {65535, "Sony FE 28mm f/2"},
            {65535, "Sony FE PZ 28-135mm f/4 G OSS"},
            {65535, "Sony FE 100mm f/2.8 STF GM OSS"},
            {65535, "Sony FE 24-70mm f/2.8 GM"},
            {65535, "Sony FE 50mm f/1.4 ZA"},
            {65535, "Sony FE 85mm f/1.4 GM"},
            {65535, "Sony FE 50mm f/1.8"},
            {65535, "Sony FE 21mm f/2.8 (SEL28F20 + SEL075UWC)"},
            {65535, "Sony FE 16mm f/3.5 Fisheye (SEL28F20 + SEL057FEC)"},
            {65535, "Sony FE 70-300mm f/4.5-5.6 G OSS"},
            {65535, "Sony FE 70-200mm f/2.8 GM OSS"},
            {65535, "Sony FE 70-200mm f/2.8 GM OSS + 1.4X Teleconverter"},
            {65535, "Sony FE 70-200mm f/2.8 GM OSS + 2X Teleconverter"},
            {65535, "Samyang AF 50mm f/1.4 FE"},
            {65535, "Samyang AF 14mm f/2.8 FE"},
            {65535, "Sigma 19mm f/2.8 [EX] DN"},
            {65535, "Sigma 30mm f/2.8 [EX] DN"},
            {65535, "Sigma 60mm f/2.8 DN"},
            {65535, "Sigma 30mm f/1.4 DC DN | C 016"},
            {65535, "Tamron 18-200mm f/3.5-6.3 Di III VC"},
            {65535, "Tokina Firin 20mm f/2 FE MF"},
            {65535, "Zeiss Touit 12mm f/2.8"},
            {65535, "Zeiss Touit 32mm f/1.8"},
            {65535, "Zeiss Touit 50mm f/2.8 Macro"},
            {65535, "Zeiss Batis 25mm f/2"},
            {65535, "Zeiss Batis 85mm f/1.8"},
            {65535, "Zeiss Batis 18mm f/2.8"},
            {65535, "Zeiss Batis 135mm f/2.8"},
            {65535, "Zeiss Loxia 50mm f/2"},
            {65535, "Zeiss Loxia 35mm f/2"},
            {65535, "Zeiss Loxia 21mm f/2.8"},
            {65535, "Zeiss Loxia 85mm f/2.4"},
            {65535, "Arax MC 35mm f/2.8 Tilt+Shift"},
            {65535, "Arax MC 80mm f/2.8 Tilt+Shift"},
            {65535, "Zenitar MF 16mm f/2.8 Fisheye M42"},
            {65535, "Samyang 500mm Mirror f/8.0"},
            {65535, "Pentacon Auto 135mm f/2.8"},
            {65535, "Pentacon Auto 29mm f/2.8"},
            {65535, "Helios 44-2 58mm f/2.0"},
        };
    }

    virtual std::string toString (Tag* t)
    {
        int lensID = t->toInt();
        Tag *lensInfoTag = t->getParent()->getRoot()->findTag ("LensInfo");
        Tag *apertureTag = t->getParent()->getRoot()->findTag ("MaxApertureValue");
        Tag *focalLengthTag = t->getParent()->getRoot()->findTag ("FocalLength");
        double maxApertureAtFocal = 0.;
        double focalLength = 0.;

        if ( apertureTag ) {
            maxApertureAtFocal = pow (2.0, apertureTag->toDouble() / 2.0);
        }

        if ( focalLengthTag ) {
            focalLength = focalLengthTag->toDouble();
        }

        double *liArray = nullptr;

        if (lensInfoTag) {
            liArray = lensInfoTag->toDoubleArray();
        }

        std::string retval = guess ( lensID, focalLength, maxApertureAtFocal, liArray);

        if (liArray) {
            delete [] liArray;
        }

        return retval;
    }
};
SALensIDInterpreter saLensIDInterpreter;

class SALensID2Interpreter : public IntLensInterpreter< int >
{
public:
    SALensID2Interpreter ()
    {
        choices.insert (p_t (0, "Unknown E-mount lens or other lens"));
        choices.insert (p_t (1, "Sony LA-EA1 Adapter"));
        choices.insert (p_t (2, "Sony LA-EA2 Adapter"));
        choices.insert (p_t (3, "Sony LA-EA3 Adapter"));
        choices.insert (p_t (6, "Sony LA-EA4 Adapter"));
        choices.insert (p_t (44, "Metabones Canon EF Smart Adapter"));
        choices.insert (p_t (78, "Metabones Canon EF Smart Adapter Mark III or Other Adapter"));
        choices.insert (p_t (234, "Metabones Canon EF Smart Adapter Mark IV"));
        choices.insert (p_t (239, "Metabones Canon EF Speed Booster"));
        choices.insert (p_t (32784, "Sony E 16mm f/2.8"));
        choices.insert (p_t (32785, "Sony E 18-55mm f/3.5-5.6 OSS"));
        choices.insert (p_t (32786, "Sony E 55-210mm f/4.5-6.3 OSS"));
        choices.insert (p_t (32787, "Sony E 18-200mm f/3.5-6.3 OSS"));
        choices.insert (p_t (32788, "Sony E 30mm f/3.5 Macro"));
        choices.insert (p_t (32789, "Sony E 24mm f/1.8 ZA or Samyang AF 50mm f/1.4 FE"));
        choices.insert (p_t (32789, "Samyang AF 50mm f/1.4 FE"));
        choices.insert (p_t (32790, "Sony E 50mm f/1.8 OSS or Samyang AF 14mm f/2.8 FE"));
        choices.insert (p_t (32790, "Samyang AF 14mm f/2.8 FE"));
        choices.insert (p_t (32791, "Sony E 16-70mm f/4 ZA OSS"));
        choices.insert (p_t (32792, "Sony E 10-18mm f/4 OSS"));
        choices.insert (p_t (32793, "Sony E PZ 16-50mm f/3.5-5.6 OSS"));
        choices.insert (p_t (32794, "Sony FE 35mm f/2.8 ZA"));
        choices.insert (p_t (32795, "Sony FE 24-70mm f/4 ZA OSS"));
        choices.insert (p_t (32796, "Sony FE 85mm f/1.8"));
        choices.insert (p_t (32797, "Sony E 18-200mm f/3.5-6.3 OSS LE"));
        choices.insert (p_t (32798, "Sony E 20mm f/2.8"));
        choices.insert (p_t (32799, "Sony E 35mm f/1.8 OSS"));
        choices.insert (p_t (32800, "Sony E PZ 18-105mm f/4 G OSS"));
        choices.insert (p_t (32802, "Sony FE 90mm f/2.8 Macro G OSS"));
        choices.insert (p_t (32803, "Sony E 18-50mm f/4-5.6"));
        choices.insert (p_t (32807, "Sony E PZ 18-200mm f/3.5-6.3 OSS"));
        choices.insert (p_t (32808, "Sony FE 55mm f/1.8 ZA"));
        choices.insert (p_t (32810, "Sony FE 70-200mm f/4 G OSS"));
        choices.insert (p_t (32811, "Sony FE 16-35mm f/4 ZA OSS"));
        choices.insert (p_t (32812, "Sony FE 50mm f/2.8 Macro"));
        choices.insert (p_t (32813, "Sony FE 28-70mm f/3.5-5.6 OSS"));
        choices.insert (p_t (32814, "Sony FE 35mm f/1.4 ZA"));
        choices.insert (p_t (32815, "Sony FE 24-240mm f/3.5-6.3 OSS"));
        choices.insert (p_t (32816, "Sony FE 28mm f/2"));
        choices.insert (p_t (32817, "Sony FE PZ 28-135mm f/4 G OSS"));
        choices.insert (p_t (32819, "Sony FE 100mm f/2.8 STF GM OSS"));
        choices.insert (p_t (32821, "Sony FE 24-70mm f/2.8 GM"));
        choices.insert (p_t (32822, "Sony FE 50mm f/1.4 ZA"));
        choices.insert (p_t (32823, "Sony FE 85mm f/1.4 GM"));
        choices.insert (p_t (32824, "Sony FE 50mm f/1.8"));
        choices.insert (p_t (32826, "Sony FE 21mm f/2.8 (SEL28F20 + SEL075UWC)"));
        choices.insert (p_t (32827, "Sony FE 16mm f/3.5 Fisheye (SEL28F20 + SEL057FEC)"));
        choices.insert (p_t (32828, "Sony FE 70-300mm f/4.5-5.6 G OSS"));
        choices.insert (p_t (32830, "Sony FE 70-200mm f/2.8 GM OSS"));
        choices.insert (p_t (33002, "Sigma 85mm f/1.4 DG HSM | A 016 (+ Metabones Ver.50)"));
        choices.insert (p_t (33072, "Sony FE 70-200mm f/2.8 GM OSS + 1.4X Teleconverter"));
        choices.insert (p_t (33073, "Sony FE 70-200mm f/2.8 GM OSS + 2X Teleconverter"));
        choices.insert (p_t (33076, "Sony FE 100mm f/2.8 STF GM OSS (macro mode)"));
        choices.insert (p_t (49201, "Zeiss Touit 12mm f/2.8"));
        choices.insert (p_t (49202, "Zeiss Touit 32mm f/1.8"));
        choices.insert (p_t (49203, "Zeiss Touit 50mm f/2.8 Macro"));
        choices.insert (p_t (49216, "Zeiss Batis 25mm f/2"));
        choices.insert (p_t (49217, "Zeiss Batis 85mm f/1.8"));
        choices.insert (p_t (49218, "Zeiss Batis 18mm f/2.8"));
        choices.insert (p_t (49219, "Zeiss Batis 135mm f/2.8"));
        choices.insert (p_t (49232, "Zeiss Loxia 50mm f/2"));
        choices.insert (p_t (49233, "Zeiss Loxia 35mm f/2"));
        choices.insert (p_t (49234, "Zeiss Loxia 21mm f/2.8"));
        choices.insert (p_t (49235, "Zeiss Loxia 85mm f/2.4"));
        choices.insert (p_t (50480, "Sigma 30mm f/1.4 DC DN | C 016"));
        choices.insert (p_t (50481, "Sigma 50mm f/1.4 DG HSM | A 014 + MC-11"));
        choices.insert (p_t (50482, "Sigma 18-300mm f/3.5-6.3 DC MACRO OS HSM | C 014 + MC-11"));
        choices.insert (p_t (50483, "Sigma 18-35mm f/1.8 DC HSM | A 013 + MC-11"));
        choices.insert (p_t (50484, "Sigma 24-35mm f/2 DG HSM | A 015 + MC-11"));
        choices.insert (p_t (50486, "Sigma 150-600mm f/5-6.3 DG OS HSM | C 015 + MC-11"));
        choices.insert (p_t (50487, "Sigma 20mm f/1.4 DG HSM | A 015 + MC-11"));
        choices.insert (p_t (50488, "Sigma 35mm f/1.4 DG HSM | A 012 + MC-11"));
        choices.insert (p_t (50489, "Sigma 150-600mm f/5-6.3 DG OS HSM | S 014 + MC-11"));
        choices.insert (p_t (50490, "Sigma 120-300mm f/2.8 DG OS HSM | S 013 + MC-11"));
        choices.insert (p_t (50492, "Sigma 24-105mm f/4 DG OS HSM | A 013 + MC-11"));
        choices.insert (p_t (50493, "Sigma 17-70mm f/2.8-4 DC MACRO OS HSM | C 013 + MC-11"));
        choices.insert (p_t (50495, "Sigma 50-100mm f/1.8 DC HSM | A 016 + MC-11"));
        choices.insert (p_t (50992, "Voigtlander SUPER WIDE-HELIAR 15mm f/4.5 III"));
        choices.insert (p_t (50993, "Voigtlander HELIAR-HYPER WIDE 10mm f/5.6"));
        choices.insert (p_t (50994, "Voigtlander ULTRA WIDE-HELIAR 12mm f/5.6 III"));
        choices.insert (p_t (50996, "Voigtlander NOKTON 40mm f/1.2 Aspherical"));
        choices.insert (p_t (51505, "Samyang AF 14mm f/2.8 FE"));
    }

    virtual std::string toString (Tag* t)
    {
        int lensID = t->toInt();
        Tag *lensInfoTag = t->getParent()->getRoot()->findTag ("LensInfo");
        Tag *apertureTag = t->getParent()->getRoot()->findTag ("MaxApertureValue");
        Tag *focalLengthTag = t->getParent()->getRoot()->findTag ("FocalLength");
        double maxApertureAtFocal = 0.;
        double focalLength = 0.;

        if ( apertureTag ) {
            maxApertureAtFocal = pow (2.0, apertureTag->toDouble() / 2.0);
        }

        if ( focalLengthTag ) {
            focalLength = focalLengthTag->toDouble();
        }

        double *liArray = nullptr;

        if (lensInfoTag) {
            liArray = lensInfoTag->toDoubleArray();
        }

        std::string retval = guess ( lensID, focalLength, maxApertureAtFocal, liArray);

        if (liArray) {
            delete [] liArray;
        }

        return retval;
    }
};
SALensID2Interpreter saLensID2Interpreter;

class MATeleconverterInterpreter : public ChoiceInterpreter
{
public:
    MATeleconverterInterpreter ()
    {
        choices[0x0]  = "None";
        choices[0x4]  = "Minolta/Sony AF 1.4x APO (D) (0x04)";
        choices[0x5]  = "Minolta/Sony AF 2x APO (D) (0x05)";
        choices[0x48] = "Minolta/Sony AF 2x APO (D)";
        choices[0x50] = "Minolta AF 2x APO II";
        choices[0x60] = "Minolta AF 2x APO";
        choices[0x88] = "Minolta/Sony AF 1.4x APO (D)";
        choices[0x90] = "Minolta AF 1.4x APO II";
        choices[0xa0] = "Minolta AF 1.4x APO";
    }
};
MATeleconverterInterpreter maTeleconverterInterpreter;

class MAQualityInterpreter : public ChoiceInterpreter
{
public:
    MAQualityInterpreter ()
    {
        choices[0] = "RAW";
        choices[1] = "Super Fine";
        choices[2] = "Fine";
        choices[3] = "Standard";
        choices[4] = "Economy";
        choices[5] = "Extra Fine";
        choices[6] = "RAW + JPEG";
        choices[7] = "Compressed RAW";
        choices[8] = "Compressed RAW + JPEG";
    }
};
MAQualityInterpreter maQualityInterpreter;

class MAImageSizeInterpreter : public ChoiceInterpreter
{
public:
    MAImageSizeInterpreter ()
    {
        choices[1]  = "1600x1200";
        choices[2]  = "1280x960";
        choices[3]  = "640x480";
        choices[5]  = "2560x1920";
        choices[6]  = "2272x1704";
        choices[7]  = "2048x1536";
    }
};
MAImageSizeInterpreter maImageSizeInterpreter;

class SAQualityInterpreter2 : public ChoiceInterpreter
{
public:
    SAQualityInterpreter2 ()
    {
        choices[0]  = "Raw";
        choices[2]  = "cRAW";
        choices[16] = "Extra fine";
        choices[32]  = "Fine";
        choices[34]  = "RAW + JPEG";
        choices[35]  = "cRAW + JPEG";
        choices[48]  = "Standard";
    }
};
SAQualityInterpreter2 saQualityInterpreter2;

class SAQualityInterpreter3 : public ChoiceInterpreter
{
public:
    SAQualityInterpreter3 ()
    {
        choices[2]  = "RAW";
        choices[4]  = "RAW + JPEG";
        choices[6]  = "Fine";
        choices[7]  = "Standard";
    }
};
SAQualityInterpreter3 saQualityInterpreter3;

class SADriveMode : public ChoiceInterpreter
{
public:
    SADriveMode ()
    {
        choices[1]  = "Single Frame";
        choices[2]  = "Continuous High";
        choices[4]  = "Self-timer 10 sec";
        choices[5]  = "Self-timer 2 sec, Mirror Lock-up";
        choices[6]  = "Single-frame Bracketing";
        choices[7]  = "Continuous Bracketing";
        choices[10] = "Remote Commander";
        choices[11] = "Mirror Lock-up";
        choices[18] = "Continuous Low";
        choices[24] = "White Balance Bracketing Low";
        choices[25] = "D-Range Optimizer Bracketing Low";
        choices[40] = "White Balance Bracketing High";
        choices[41] = "D-Range Optimizer Bracketing High";
    }
};
SADriveMode saDriveMode;

class SADriveMode2 : public ChoiceInterpreter
{
public:
    SADriveMode2 ()
    {
        choices[1]  = "Single Frame";
        choices[2]  = "Continuous High";
        choices[4]  = "Self-timer 10 sec";
        choices[5]  = "Self-timer 2 sec, Mirror Lock-up";
        choices[7]  = "Continuous Bracketing";
        choices[10] = "Remote Commander";
        choices[11] = "Continuous Self-timer";
    }
};
SADriveMode2 saDriveMode2;

class SADriveMode3 : public ChoiceInterpreter
{
public:
    SADriveMode3 ()
    {
        choices[16]  = "Single Frame";
        choices[33]  = "Continuous High";
        choices[34]  = "Continuous Low";
        choices[48]  = "Speed Priority Continuous";
        choices[81]  = "Self-timer 10 sec";
        choices[82]  = "Self-timer 2 sec, Mirror Lock-up";
        choices[113] = "Continuous Bracketing 0.3 EV";
        choices[117] = "Continuous Bracketing 0.7 EV";
        choices[145] = "White Balance Bracketing Low";
        choices[146] = "White Balance Bracketing High";
        choices[192] = "Remote Commander";
        choices[209] = "Continuous - HDR";
        choices[210] = "Continuous - Multi Frame NR";
        choices[211] = "Continuous - Handheld Night Shot";
        choices[212] = "Continuous - Anti Motion Blur";
        choices[213] = "Continuous - Sweep Panorama";
        choices[214] = "Continuous - 3D Sweep Panorama";
    }
};
SADriveMode3 saDriveMode3;

class SAFocusMode: public ChoiceInterpreter
{
public:
    SAFocusMode ()
    {
        choices[0]  = "Manual";
        choices[1]  = "AF-S";
        choices[2]  = "AF-C";
        choices[3]  = "AF-A";
        choices[4]  = "Permanent-AF";
        choices[65535] = "n/a";
    }
};
SAFocusMode saFocusMode;

class SAFocusMode2: public ChoiceInterpreter
{
public:
    SAFocusMode2 ()
    {
        choices[0]  = "Manual";
        choices[1]  = "AF-S";
        choices[2]  = "AF-C";
        choices[3]  = "AF-A";
        choices[65535] = "n/a";
    }
};
SAFocusMode2 saFocusMode2;

class SAFocusModeSetting3: public ChoiceInterpreter
{
public:
    SAFocusModeSetting3 ()
    {
        choices[17]  = "AF-S";
        choices[18]  = "AF-C";
        choices[19]  = "AF-A";
        choices[32]  = "Manual";
        choices[48]  = "DMF";
        choices[65535] = "n/a";
    }
};
SAFocusModeSetting3 saFocusModeSetting3;

class SAAFMode: public ChoiceInterpreter
{
public:
    SAAFMode()
    {
        choices[0] = "Default";
        choices[1] = "Multi AF";
        choices[2] = "Center AF";
        choices[3] = "Spot AF";
        choices[4] = "Flexible Spot AF";
        choices[6] = "Touch AF";
        choices[14] = "Tracking";
        choices[15] = "Face Tracking";
        choices[65535] = "n/a";
    }
};
SAAFMode saAFMode;

class SAAFAreaMode: public ChoiceInterpreter
{
public:
    SAAFAreaMode ()
    {
        choices[0]  = "Wide";
        choices[1]  = "Local";
        choices[2]  = "Spot";
    }
};
SAAFAreaMode saAFAreaMode;

class SAAFAreaMode2: public ChoiceInterpreter
{
public:
    SAAFAreaMode2 ()
    {
        choices[1]  = "Wide";
        choices[2]  = "Spot";
        choices[3]  = "Local";
        choices[4]  = "Flexible";
    }
};
SAAFAreaMode2 saAFAreaMode2;

class SAAFPointSelected: public ChoiceInterpreter
{
public:
    SAAFPointSelected ()
    {
        choices[1] = "Center";
        choices[2] = "Top";
        choices[3] = "Top-Right";
        choices[4] = "Right";
        choices[5] = "Bottom-Right";
        choices[6] = "Bottom";
        choices[7] = "Bottom-Left";
        choices[8] = "Left";
        choices[9] = "Top-Left";
        choices[10] = "Far Right";
        choices[11] = "Far Left";
    }
};
SAAFPointSelected saAFPointSelected;

class SACameraInfoAFPointSelected: public ChoiceInterpreter
{
public:
    SACameraInfoAFPointSelected ()
    {
        choices[0] = "Auto";
        choices[1] = "Center";
        choices[2] = "Top";
        choices[3] = "Upper-Right";
        choices[4] = "Right";
        choices[5] = "Lower-Right";
        choices[6] = "Bottom";
        choices[7] = "Lower-Left";
        choices[8] = "Left";
        choices[9] = "Upper-Left";
        choices[10] = "Far Right";
        choices[11] = "Far Left";
        choices[12] = "Upper-middle";
        choices[13] = "Near Right";
        choices[14] = "Lower-middle";
        choices[15] = "Near Left";
    }
};
SACameraInfoAFPointSelected saCameraInfoAFPointSelected;

class SACameraInfoAFPoint: public ChoiceInterpreter
{
public:
    SACameraInfoAFPoint ()
    {
        choices[0] = "Upper-Left";
        choices[1] = "Left";
        choices[2] = "Lower-Left";
        choices[3] = "Far Left";
        choices[4] = "Top (horizontal)";
        choices[5] = "Near Right";
        choices[6] = "Center (horizontal)";
        choices[7] = "Near Left";
        choices[8] = "Bottom (horizontal)";
        choices[9] = "Top (vertical)";
        choices[10] = "Center (vertical)";
        choices[11] = "Bottom (vertical)";
        choices[12] = "Far Right";
        choices[13] = "Upper-Right";
        choices[14] = "Right";
        choices[15] = "Lower-Right";
        choices[16] = "Upper-middle";
        choices[17] = "Lower-middle";
        choices[255] = "(none)";
    }
};
SACameraInfoAFPoint saCameraInfoAFPoint;

class SAAFPointSelected2: public ChoiceInterpreter
{
public:
    SAAFPointSelected2 ()
    {
        choices[1] = "Center";
        choices[2] = "Top";
        choices[3] = "Top-Right";
        choices[4] = "Right";
        choices[5] = "Bottom-Right";
        choices[6] = "Bottom";
        choices[7] = "Bottom-Left";
        choices[8] = "Left";
        choices[9] = "Top-Left";
    }
};
SAAFPointSelected2 saAFPointSelected2;

class SAMeteringMode0_3: public ChoiceInterpreter
{
public:
    SAMeteringMode0_3 ()
    {
        choices[0] = "Multi-segment";
        choices[2] = "Center-weighted Average";
        choices[3] = "Spot";
    }
};
SAMeteringMode0_3 saMeteringMode0_3;

class SAMeteringMode1_3: public ChoiceInterpreter
{
public:
    SAMeteringMode1_3 ()
    {
        choices[1] = "Multi-segment";
        choices[2] = "Center-weighted Average";
        choices[3] = "Spot";
    }
};
SAMeteringMode1_3 saMeteringMode1_3;

class SAMeteringMode1_4: public ChoiceInterpreter
{
public:
    SAMeteringMode1_4 ()
    {
        choices[1] = "Multi-segment";
        choices[2] = "Center-weighted Average";
        choices[4] = "Spot";
    }
};
SAMeteringMode1_4 saMeteringMode1_4;

class SADynamicRangeOptimizerMode: public ChoiceInterpreter
{
public:
    SADynamicRangeOptimizerMode ()
    {
        choices[0] = "Off";
        choices[1] = "Standard";
        choices[2] = "Advanced Auto";
        choices[3] = "Advanced Level";
        choices[4097] = "Auto";
    }
};
SADynamicRangeOptimizerMode saDynamicRangeOptimizerMode;

class SADynamicRangeOptimizerSetting: public ChoiceInterpreter
{
public:
    SADynamicRangeOptimizerSetting ()
    {
        choices[1] = "Off";
        choices[2] = "On (Auto)";
        choices[3] = "On (Manual)";
    }
};
SADynamicRangeOptimizerSetting saDynamicRangeOptimizerSetting;

class SACreativeStyle: public ChoiceInterpreter
{
public:
    SACreativeStyle ()
    {
        choices[1] = "Standard";
        choices[2] = "Vivid";
        choices[3] = "Portrait";
        choices[4] = "Landscape";
        choices[5] = "Sunset";
        choices[6] = "Night View/Portrait";
        choices[8] = "B&W";
        choices[9] = "Adobe RGB";
        choices[11] = "Neutral";
        choices[12] = "Clear";
        choices[13] = "Deep";
        choices[14] = "Light";
        choices[15] = "Autumn Leaves";
        choices[16] = "Sepia";
    }
};
SACreativeStyle saCreativeStyle;

class SACreativeStyle2: public ChoiceInterpreter
{
public:
    SACreativeStyle2 ()
    {
        choices[1] = "Standard";
        choices[2] = "Vivid";
        choices[3] = "Portrait";
        choices[4] = "Landscape";
        choices[5] = "Sunset";
        choices[6] = "Night View/Portrait";
        choices[8] = "B&W";
    }
};
SACreativeStyle2 saCreativeStyle2;

class SACreativeStyleSetting: public ChoiceInterpreter
{
public:
    SACreativeStyleSetting ()
    {
        choices[16]  = "Standard";
        choices[32]  = "Vivid";
        choices[64]  = "Portrait";
        choices[80]  = "Landscape";
        choices[96]  = "B&W";
        choices[160] = "Sunset";
    }
};
SACreativeStyleSetting saCreativeStyleSetting;

class SAFlashControl: public ChoiceInterpreter
{
public:
    SAFlashControl ()
    {
        choices[1] = "ADI Flash";
        choices[2] = "Pre-flash TTL";
    }
};
SAFlashControl saFlashControl;

class SAFlashMode: public ChoiceInterpreter
{
public:
    SAFlashMode ()
    {
        choices[0] = "ADI";
        choices[1] = "TTL";
    }
};
SAFlashMode saFlashMode;

class SAFlashMode2: public ChoiceInterpreter
{
public:
    SAFlashMode2 ()
    {
        choices[1] = "Flash Off";
        choices[16] = "Autoflash";
        choices[17] = "Fill-flash";
        choices[18] = "Slow Sync";
        choices[19] = "Rear Sync";
        choices[20] = "Wireless";
    }
};
SAFlashMode2 saFlashMode2;

class SAExposureProgram: public ChoiceInterpreter
{
public:
    SAExposureProgram ()
    {
        choices[0] = "Auto";
        choices[1] = "Manual";
        choices[2] = "Program AE";
        choices[3] = "Aperture-priority AE";
        choices[4] = "Shutter speed priority AE";
        choices[8] = "Program Shift A";
        choices[9] = "Program Shift S";
        choices[16] = "Portrait";
        choices[17] = "Sports";
        choices[18] = "Sunset";
        choices[19] = "Night Portrait";
        choices[20] = "Landscape";
        choices[21] = "Macro";
        choices[35] = "Auto No Flash";
    }
};
SAExposureProgram saExposureProgram;

class SAExposureProgram2: public ChoiceInterpreter
{
public:
    SAExposureProgram2 ()
    {
        choices[1]   = "Program AE";
        choices[2]   = "Aperture-priority AE";
        choices[3]   = "Shutter speed priority AE";
        choices[4]   = "Manual";
        choices[5]   = "Cont. Priority AE";
        choices[16]  = "Auto";
        choices[17]  = "Auto (no flash)";
        choices[18]  = "Auto+";
        choices[49]  = "Portrait";
        choices[50]  = "Landscape";
        choices[51]  = "Macro";
        choices[52]  = "Sports";
        choices[53]  = "Sunset";
        choices[54]  = "Night view";
        choices[55]  = "Night view/portrait";
        choices[56]  = "Handheld Night Shot";
        choices[57]  = "3D Sweep Panorama";
        choices[64]  = "Auto 2";
        choices[65]  = "Auto 2 (no flash)";
        choices[80]  = "Sweep Panorama";
        choices[96]  = "Anti Motion Blur";
        choices[128] = "Toy Camera";
        choices[129] = "Pop Color";
        choices[130] = "Posterization";
        choices[131] = "Posterization B/W";
        choices[132] = "Retro Photo";
        choices[133] = "High-key";
        choices[134] = "Partial Color Red";
        choices[135] = "Partial Color Green";
        choices[136] = "Partial Color Blue";
        choices[137] = "Partial Color Yellow";
        choices[138] = "High Contrast Monochrome";
    }
};
SAExposureProgram2 saExposureProgram2;

class SARotation: public ChoiceInterpreter
{
public:
    SARotation ()
    {
        choices[0] = "Horizontal";
        choices[1] = "Rotate 90 CW";
        choices[2] = "Rotate 270 CW";
        choices[3] = "None";
    }
};
SARotation saRotation;

class SASonyImageSize: public ChoiceInterpreter
{
public:
    SASonyImageSize ()
    {
        choices[1] = "Large";
        choices[2] = "Medium";
        choices[3] = "Small";
    }
};
SASonyImageSize saSonyImageSize;

class SASonyImageSize3: public ChoiceInterpreter
{
public:
    SASonyImageSize3 ()
    {
        choices[21] = "Large (3:2)";
        choices[22] = "Medium (3:2)";
        choices[23] = "Small (3:2)";
        choices[25] = "Large (16:9)";
        choices[26] = "Medium (16:9) ";
        choices[27] = "Small (16:9)";
    }
};
SASonyImageSize3 saSonyImageSize3;

class SAAspectRatio: public ChoiceInterpreter
{
public:
    SAAspectRatio ()
    {
        choices[1] = "3:2";
        choices[2] = "16:9";
    }
};
SAAspectRatio saAspectRatio;

class SAAspectRatio2: public ChoiceInterpreter
{
public:
    SAAspectRatio2 ()
    {
        choices[4] = "3:2";
        choices[8] = "16:9";
    }
};
SAAspectRatio2 saAspectRatio2;

class SAExposureLevelIncrements: public ChoiceInterpreter
{
public:
    SAExposureLevelIncrements ()
    {
        choices[33] = "1/3 EV";
        choices[50] = "1/2 EV";
    }
};
SAExposureLevelIncrements saExposureLevelIncrements;

class SAAFIlluminator: public ChoiceInterpreter
{
public:
    SAAFIlluminator ()
    {
        choices[0] = "Off";
        choices[1] = "Auto";
        choices[65535] = "n/a";
    }
};
SAAFIlluminator saAFIlluminator;

class SAColorSpace1_2: public ChoiceInterpreter
{
public:
    SAColorSpace1_2 ()
    {
        choices[1] = "sRGB";
        choices[2] = "AdobeRGB";
    }
};
SAColorSpace1_2 saColorSpace1_2;

class SAColorSpace0_5: public ChoiceInterpreter
{
public:
    SAColorSpace0_5 ()
    {
        choices[0] = "sRGB";
        choices[1] = "AdobeRGB";
        choices[5] = "AdobeRGB";
    }
};
SAColorSpace0_5 saColorSpace0_5;

class SAColorSpace5_6: public ChoiceInterpreter
{
public:
    SAColorSpace5_6 ()
    {
        choices[5] = "AdobeRGB";
        choices[6] = "sRGB";
    }
};
SAColorSpace5_6 saColorSpace5_6;

class SAReleaseModeInterpreter: public ChoiceInterpreter
{
public:
    SAReleaseModeInterpreter ()
    {
        choices[0] = "Normal";
        choices[2] = "Continuous";
        choices[5] = "Exposure Bracketing";
        choices[6] = "White Balance Bracketing";
        choices[8] = "DRO Bracketing";
        choices[65535] = "n/a";
    }
};
SAReleaseModeInterpreter saReleaseModeInterpreter;

class SAImageStyleInterpreter: public ChoiceInterpreter
{
public:
    SAImageStyleInterpreter ()
    {
        choices[1] = "Standard";
        choices[2] = "Vivid";
        choices[3] = "Portrait";
        choices[4] = "Landscape";
        choices[5] = "Sunset";
        choices[7] = "Night View/Portrait";
        choices[8] = "B&W";
        choices[9] = "Adobe RGB";
        choices[11] = "Neutral";
        choices[129] = "StyleBox1";
        choices[130] = "StyleBox2";
        choices[131] = "StyleBox3";
        choices[132] = "StyleBox4";
        choices[133] = "StyleBox5";
        choices[134] = "StyleBox6";
    }
};
SAImageStyleInterpreter saImageStyleInterpreter;

class SAPictureEffectInterpreter: public ChoiceInterpreter
{
public:
    SAPictureEffectInterpreter()
    {
        choices[0]  = "Off";
        choices[1]  = "Toy Camera";
        choices[2]  = "Pop Color";
        choices[3]  = "Posterization";
        choices[4]  = "Posterization B/W";
        choices[5]  = "Retro Photo";
        choices[6]  = "Soft High Key";
        choices[7]  = "Partial Color (red)";
        choices[8]  = "Partial Color (green)";
        choices[9]  = "Partial Color (blue)";
        choices[10] = "Partial Color (yellow)";
        choices[13] = "High Contrast Monochrome";
        choices[16] = "Toy Camera (normal)";
        choices[17] = "Toy Camera (cool)";
        choices[18] = "Toy Camera (warm)";
        choices[19] = "Toy Camera (green)";
        choices[20] = "Toy Camera (magenta)";
        choices[32] = "Soft Focus (low)";
        choices[33] = "Soft Focus";
        choices[34] = "Soft Focus (high)";
        choices[48] = "Miniature (auto)";
        choices[49] = "Miniature (top)";
        choices[50] = "Miniature (middle horizontal)";
        choices[51] = "Miniature (bottom)";
        choices[52] = "Miniature (left)";
        choices[53] = "Miniature (middle vertical)";
        choices[54] = "Miniature (right)";
        choices[64] = "HDR Painting (low)";
        choices[65] = "HDR Painting";
        choices[66] = "HDR Painting (high)";
        choices[80] = "Rich-tone Monochrome";
        choices[97] = "Water Color";
        choices[98] = "Water Color 2";
        choices[112] = "Illustration (low)";
        choices[113] = "Illustration";
        choices[114] = "Illustration (high)";
    }
};
SAPictureEffectInterpreter saPictureEffectInterpreter;

class SACameraInfoFocusStatusInterpreter : public ChoiceInterpreter
{
public:
    SACameraInfoFocusStatusInterpreter()
    {
        choices[0] = "Manual - Not confirmed (0)";
        choices[4] = "Manual - Not confirmed (4)";
        choices[16] = "AF-C - Confirmed";
        choices[24] = "AF-C - Not Confirmed";
        choices[64] = "AF-S - Confirmed";
    }
};
SACameraInfoFocusStatusInterpreter saCameraInfoFocusStatusInterpreter;

class SAExposureTimeInterpreter : public Interpreter
{
public:
    SAExposureTimeInterpreter () {}
    virtual std::string toString (Tag* t)
    {
        double a = t->toDouble();

        if (a > 0) {
            char buffer[32];
            sprintf (buffer, "%.4f", a);
            return buffer;
        } else {
            return "n/a";
        }
    }
    virtual double toDouble (Tag* t, int ofs)
    {
        // Get the value; Depending on the camera model, this parameter can be a BYTE or a SHORT
        TagType astype = t->getType();
        int a = 0;

        if (astype == BYTE) {
            a = t->getValue()[ofs];
        } else if (astype == SHORT) {
            a = (int)sget2 (t->getValue() + ofs, t->getOrder());
        }

        // Decode the value
        if (a > 0) {
            return pow (2., 6. - (double (a) / 8.));
        } else {
            return 0.;
        }
    }
    virtual int toInt (Tag* t, int ofs, TagType astype)
    {
        // Get the value; Depending on the camera model, this parameter can be a BYTE or a SHORT
        int a = 0;

        if (astype == INVALID || astype == AUTO) {
            astype = t->getType();
        }

        if (astype == BYTE) {
            a = t->getValue()[ofs];
        } else if (astype == SHORT) {
            a = (int)sget2 (t->getValue() + ofs, t->getOrder());
        }

        // Decode the value
        if (a) {
            return int (powf (2.f, 6.f - (float (a) / 8.f)) + 0.5f);
        } else {
            return 0;
        }
    }
};
SAExposureTimeInterpreter saExposureTimeInterpreter;

class SAFNumberInterpreter : public Interpreter
{
public:
    SAFNumberInterpreter () {}
    virtual std::string toString (Tag* t)
    {
        double a = double (t->toDouble());

        if (a) {
            char buffer[32];
            sprintf (buffer, "%.1f", a / 100. );
            return buffer;
        } else {
            return "n/a";
        }
    }
    virtual double toDouble (Tag* t, int ofs)
    {
        // Get the value; Depending on the camera model, this parameter can be a BYTE or a SHORT
        TagType astype = t->getType();
        int a = 0;

        if (astype == BYTE) {
            a = t->getValue()[ofs];
        } else if (astype == SHORT) {
            a = (int)sget2 (t->getValue() + ofs, t->getOrder());
        }

        // Decode the value
        if (a > 0) {
            return pow (2., (double (a) / 8. - 1.) / 2.);
        } else {
            return 0.;
        }
    }
    virtual int toInt (Tag* t, int ofs, TagType astype)
    {
        // Get the value; Depending on the camera model, this parameter can be a BYTE or a SHORT
        int a = 0;

        if (astype == INVALID || astype == AUTO) {
            astype = t->getType();
        }

        if (astype == BYTE) {
            a = t->getValue()[ofs];
        } else if (astype == SHORT) {
            a = (int)sget2 (t->getValue() + ofs, t->getOrder());
        }

        // Decode the value
        if (a) {
            return int (powf (2.f, (float (a) / 8.f - 1.f) / 2.f) + 0.5f);
        } else {
            return 0;
        }
    }
};
SAFNumberInterpreter saFNumberInterpreter;

class SAISOSettingInterpreter : public Interpreter
{
public:
    SAISOSettingInterpreter () {}
    virtual std::string toString (Tag* t)
    {
        int a = t->toInt();

        if (a) {
            char buffer[32];
            sprintf (buffer, "%d", a );
            return buffer;
        } else {
            return "Auto";
        }
    }
    virtual int toInt (Tag* t, int ofs, TagType astype)
    {
        // Get the value; Depending on the camera model, this parameter can be a BYTE or a SHORT
        int a = 0;

        if (astype == INVALID || astype == AUTO) {
            astype = t->getType();
        }

        if (astype == BYTE) {
            a = t->getValue()[ofs];
        } else if (astype == SHORT) {
            a = (int)sget2 (t->getValue() + ofs, t->getOrder());
        }

        // Decode the value
        if (a && a != 254) { // 254 = 'Auto' for CameraSettings3, but we might say the same for CameraSettings & CameraSettings2 (?)
            return int (expf ((double (a) / 8.f - 6.f) * logf (2.f)) * 100.f + 0.5f);
        } else {
            return 0;
        }
    }
};
SAISOSettingInterpreter saISOSettingInterpreter;

class SAExposureCompSetInterpreter : public Interpreter
{
public:
    SAExposureCompSetInterpreter () {}
    virtual std::string toString (Tag* t)
    {
        double a = t->toDouble();
        char buffer[32];
        sprintf (buffer, "%.2f", a );
        return buffer;
    }
    virtual double toDouble (Tag* t, int ofs)
    {
        // Get the value
        int a = t->getValue()[ofs];
        // Decode the value
        return (double (a) - 128.) / 24.;
    }
};
SAExposureCompSetInterpreter saExposureCompSetInterpreter;

class SAAFMicroAdjValueInterpreter : public Interpreter
{
public:
    SAAFMicroAdjValueInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        char buffer[32];
        sprintf (buffer, "%d", t->getValue()[0] - 20);
        return buffer;
    }
    virtual int toInt (Tag* t, int ofs, TagType astype)
    {
        return t->getValue()[0] - 20;
    }
};
SAAFMicroAdjValueInterpreter saAFMicroAdjValueInterpreter;

class SAAFMicroAdjModeInterpreter : public Interpreter
{
public:
    SAAFMicroAdjModeInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        int a = t->getValue()[0] & 0x80;

        if (a == 0x80) {
            return "On";
        }

        return "Off";
    }
    virtual int toInt (Tag* t, int ofs, TagType astype)
    {
        return (t->getValue()[0] & 0x80) == 0x80 ? 1 : 0;
    }
};

SAAFMicroAdjModeInterpreter saAFMicroAdjModeInterpreter;

class SAAFMicroAdjRegisteredLensesInterpreter : public Interpreter
{
public:
    SAAFMicroAdjRegisteredLensesInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        char buffer[32];
        sprintf (buffer, "%d", t->getValue()[0] & 0x7f);
        return buffer;
    }
    virtual int toInt (Tag* t, int ofs, TagType astype)
    {
        return t->getValue()[0] & 0x7f;
    }
};
SAAFMicroAdjRegisteredLensesInterpreter saAFMicroAdjRegisteredLensesInterpreter;

class SAFocusStatusInterpreter : public Interpreter
{
public:
    SAFocusStatusInterpreter () {}
    virtual std::string toString (Tag* t)
    {
        std::string retval;
        int a = t->toInt();

        if (a == 0) {
            retval = "Not confirmed";
        } else if (a == 4) {
            retval = "Not confirmed, Tracking";
        } else {
            if (a & 1) {
                retval = "Confirmed";
            }

            if (a & 2) {
                if (!retval.empty()) {
                    retval += ", ";
                }

                retval += "Failed";
            }

            if (a & 4)
                if (!retval.empty()) {
                    retval += ", ";
                }

            retval += "Tracking";
        }

        return retval;
    }
};
SAFocusStatusInterpreter saFocusStatusInterpreter;

class SAColorTemperatureSettingInterpreter : public Interpreter
{
public:
    SAColorTemperatureSettingInterpreter () {}
    virtual std::string toString (Tag* t)
    {
        char buffer[32];
        sprintf (buffer, "%d", t->toInt());
        return buffer;
    }
    virtual int toInt (Tag* t, int ofs, TagType astype)
    {
        int a = 0;

        if (astype == INVALID || astype == AUTO) {
            astype = t->getType();
        }

        if (astype == BYTE) {
            a = t->getValue()[ofs];
        } else if (astype == SHORT) {
            a = (int)sget2 (t->getValue() + ofs, t->getOrder());
        }

        return a * 100;
    }
};
SAColorTemperatureSettingInterpreter saColorTemperatureSettingInterpreter;

const TagAttrib minoltaAttribs[] = {
    {0, AC_WRITE, 0, nullptr, 0x0000, AUTO, "MakerNoteVersion", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0001, AUTO, "MinoltaCameraSettingsOld", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0003, AUTO, "MinoltaCameraSettings", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0004, AUTO, "MinoltaCameraSettings7D", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0018, AUTO, "ImageStabilization", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0040, AUTO, "CompressedImageSize", &stdInterpreter},
    {1, AC_WRITE, 0, nullptr, 0x0081, AUTO, "PreviewImage", &stdInterpreter},
    {1, AC_WRITE, 0, nullptr, 0x0088, AUTO, "PreviewImageStart", &stdInterpreter},
    {1, AC_WRITE, 0, nullptr, 0x0089, AUTO, "PreviewImageLength", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0100, AUTO, "SceneMode", &saSceneModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0101, AUTO, "ColorMode", &saColorModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0102, AUTO, "MinoltaQuality", &maQualityInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0103, AUTO, "MinoltaImageSize", &maImageSizeInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0104, AUTO, "FlashExposureComp", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0105, AUTO, "Teleconverter", &maTeleconverterInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0107, AUTO, "ImageStabilization", &saOnOffInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x010a, AUTO, "ZoneMatching", &saZoneMatchingInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x010b, AUTO, "ColorTemperature", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x010c, AUTO, "LensID", &saLensIDInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0113, AUTO, "ImageStabilization", &saOnOffInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0114, AUTO, "MinoltaCameraSettings", &stdInterpreter},
    {1, AC_WRITE, 0, nullptr, 0x0e00, AUTO, "PrintIM", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0f00, AUTO, "MinoltaCameraSettings2", &stdInterpreter},
    { -1, AC_DONTWRITE, 0,  nullptr, 0, AUTO, "", nullptr}
};

const TagAttrib sonyAttribs[] = {
    {0, AC_WRITE, 0, sonyCameraInfoAttribs, 0x0010, AUTO, "CameraInfo", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0102, AUTO, "Quality", &maQualityInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0104, AUTO, "FlashExposureComp", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0106, AUTO, "TeleConverter", &maTeleconverterInterpreter},
    {0, AC_WRITE, 0, sonyCameraSettingsAttribs, 0x0114, AUTO, "SonyCameraSettings", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0115, AUTO, "WhiteBalance", &saWhiteBalanceInterpreter},
    {1, AC_WRITE, 0, nullptr, 0x0e00, AUTO, "PrintIM", &stdInterpreter},
    {1, AC_WRITE, 0, nullptr, 0x2001, AUTO, "PreviewImage", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x2009, AUTO, "HighISONoiseReduction", &saHighISONoiseReduction},
    {0, AC_WRITE, 0, nullptr, 0x200a, AUTO, "AutoHDR", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x200b, AUTO, "MultiFrameNoiseReduction", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x200e, AUTO, "PictureEffect", &saPictureEffectInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x2011, AUTO, "VignettingCorrection", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x2012, AUTO, "LateralChromaticAberration", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x2013, AUTO, "DistortionCorrection", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb020, AUTO, "ColorReproduction", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb021, AUTO, "ColorTemperature", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb022, AUTO, "ColorCompensationFilter", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb023, AUTO, "SceneMode", &saSceneModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb024, AUTO, "ZoneMatching", &saZoneMatchingInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb025, AUTO, "DynamicRangeOptimizer", &saDynamicRangeOptimizerInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb026, AUTO, "ImageStabilization", &saOnOffInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb027, AUTO, "LensID", &saLensIDInterpreter},
    {0, AC_WRITE, 0, minoltaAttribs, 0xb028, AUTO, "MinoltaMakerNote", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb029, AUTO, "ColorMode", &saColorModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb040, AUTO, "Macro", &saOnOffInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb041, AUTO, "ExposureMode", &saExposureModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb042, AUTO, "FocusMode", &saFocusMode},
    {0, AC_WRITE, 0, nullptr, 0xb043, AUTO, "AFMode", &saAFMode},
    {0, AC_WRITE, 0, nullptr, 0xb044, AUTO, "AFIlluminator", &saAFIlluminator},
    {0, AC_WRITE, 0, nullptr, 0xb047, AUTO, "Quality", &saQualityInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb048, AUTO, "FlashLevel", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb049, AUTO, "ReleaseMode", &saReleaseModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb04a, AUTO, "SequenceNumber", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb04b, AUTO, "AntiBlur", &saAntiBlurInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb04e, AUTO, "LongExposureNoiseReduction", &saOnOffInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb04f, AUTO, "DynamicRangeOptimizer", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0xb050, AUTO, "HighISONoiseReduction2", &saHighISONoiseReduction2},
    {0, AC_WRITE, 0, nullptr, 0xb052, AUTO, "IntelligentAuto", &stdInterpreter},
    {0, AC_WRITE, 0, sonyTag9405Attribs, 0x9405, AUTO, "Tag9405", &stdInterpreter},
    { -1, AC_DONTWRITE, 0,  nullptr, 0, AUTO, "", nullptr}
};

const TagAttrib sonyTag9405Attribs[] = {
    {0, AC_WRITE, 0, nullptr, 0x005d, AUTO, "LensFormat", &stdInterpreter},  // 9405b start here
    {0, AC_WRITE, 0, nullptr, 0x005e, AUTO, "LensMount", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0060, SHORT, "LensType2", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0062, SHORT, "LensType", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0603, AUTO, "LensFormat", &stdInterpreter},  // 9405a start here
    {0, AC_WRITE, 0, nullptr, 0x0604, AUTO, "LensMount", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0605, SHORT, "LensType2", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0608, SHORT, "LensType", &stdInterpreter},
    { -1, AC_DONTWRITE, 0,  nullptr, 0, AUTO, "", nullptr}
};

const TagAttrib sonyCameraInfoAttribs[] = {
    {0, AC_WRITE, 0, nullptr, 14, SHORT, "FocalLength", &saExposureTimeInterpreter},
    {0, AC_WRITE, 0, nullptr, 16, SHORT, "FocalLengthTeleZoom", &saExposureTimeInterpreter},
    {0, AC_WRITE, 0, nullptr, 25, AUTO, "FocusStatus", &saCameraInfoFocusStatusInterpreter},
    {0, AC_WRITE, 0, nullptr, 28, AUTO, "AFPointSelected", &saCameraInfoAFPointSelected},
    {0, AC_WRITE, 0, nullptr, 29, AUTO, "FocusMode", &saFocusMode2},
    {0, AC_WRITE, 0, nullptr, 32, AUTO, "AFPoint", &saCameraInfoAFPoint},
    { -1, AC_DONTWRITE, 0,  nullptr, 0, AUTO, "", nullptr}
};

const TagAttrib sonyCameraInfo2Attribs[] = {
    {0, AC_WRITE, 0, nullptr, 304, AUTO, "AFMicroAdjValue", &saAFMicroAdjValueInterpreter},
    {0, AC_WRITE, 0, nullptr, 305, AUTO, "AFMicroAdjMode", &saAFMicroAdjModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 305, AUTO, "AFMicroAdjRegisteredLenses", &saAFMicroAdjRegisteredLensesInterpreter},
    { -1, AC_DONTWRITE, 0,  nullptr, 0, AUTO, "", nullptr}
};

const TagAttrib sonyCameraSettingsAttribs[] = {
    {0, AC_WRITE, 0, nullptr,  0, AUTO, "ExposureTime", &saExposureTimeInterpreter},
    {0, AC_WRITE, 0, nullptr,  1, AUTO, "FNumber", &saFNumberInterpreter},
    {0, AC_WRITE, 0, nullptr,  4, AUTO, "DriveMode", &saDriveMode},
    {0, AC_WRITE, 0, nullptr,  6, AUTO, "WhiteBalanceFineTune", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 16, AUTO, "FocusModeSetting", &saFocusMode},
    {0, AC_WRITE, 0, nullptr, 17, AUTO, "AFAreaMode", &saAFAreaMode},
    {0, AC_WRITE, 0, nullptr, 18, AUTO, "AFPointSelected", &saAFPointSelected},
    {0, AC_WRITE, 0, nullptr, 21, AUTO, "MeteringMode", &saMeteringMode1_4},
    {0, AC_WRITE, 0, nullptr, 22, AUTO, "ISOSetting", &saISOSettingInterpreter},
    {0, AC_WRITE, 0, nullptr, 24, AUTO, "DynamicRangeOptimizerMode", &saDynamicRangeOptimizerMode},
    {0, AC_WRITE, 0, nullptr, 25, AUTO, "DynamicRangeOptimizerLevel", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 26, AUTO, "CreativeStyle", &saCreativeStyle},
    {0, AC_WRITE, 0, nullptr, 28, AUTO, "Sharpness", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 29, AUTO, "Contrast", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 30, AUTO, "Saturation", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 31, AUTO, "ZoneMatchingValue", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 34, AUTO, "Brightness", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 35, AUTO, "FlashMode", &saFlashMode},
    {0, AC_WRITE, 0, nullptr, 40, AUTO, "PrioritySetupShutterRelease", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 41, AUTO, "AFIlluminator", &saAFIlluminator},
    {0, AC_WRITE, 0, nullptr, 42, AUTO, "AFWithShutter", &saOnOffInterpreter},
    {0, AC_WRITE, 0, nullptr, 43, AUTO, "LongExposureNoiseReduction", &saOnOffInterpreter},
    {0, AC_WRITE, 0, nullptr, 44, AUTO, "HighISONoiseReduction", &saHighISONoiseReduction3},
    {0, AC_WRITE, 0, nullptr, 45, AUTO, "ImageStyle", &saImageStyleInterpreter},
    {0, AC_WRITE, 0, nullptr, 60, AUTO, "ExposureProgram", &saExposureProgram},
    {0, AC_WRITE, 0, nullptr, 61, AUTO, "ImageStabilization", &saOnOffInterpreter},
    {0, AC_WRITE, 0, nullptr, 63, AUTO, "Rotation", &saRotation},
    {0, AC_WRITE, 0, nullptr, 77, AUTO, "FocusMode", &saFocusMode},
    {0, AC_WRITE, 0, nullptr, 83, AUTO, "FocusStatus", &saFocusStatusInterpreter},
    {0, AC_WRITE, 0, nullptr, 84, AUTO, "SonyImageSize", &saSonyImageSize},
    {0, AC_WRITE, 0, nullptr, 85, AUTO, "AspectRatio", &saAspectRatio},
    {0, AC_WRITE, 0, nullptr, 86, AUTO, "Quality", &saQualityInterpreter2},
    {0, AC_WRITE, 0, nullptr, 88, AUTO, "ExposureLevelIncrements", &saExposureLevelIncrements},
    { -1, AC_DONTWRITE, 0,  nullptr, 0, AUTO, "", nullptr}
};

const TagAttrib sonyCameraSettingsAttribs2[] = {
    {0, AC_WRITE, 0, nullptr,  0, AUTO, "ExposureTime", &saExposureTimeInterpreter},
    {0, AC_WRITE, 0, nullptr,  1, AUTO, "FNumber", &saFNumberInterpreter},
    {0, AC_WRITE, 0, nullptr, 11, AUTO, "ColorTemperatureSetting", &saColorTemperatureSettingInterpreter},
    {0, AC_WRITE, 0, nullptr, 15, AUTO, "FocusMode", &saFocusMode2},
    {0, AC_WRITE, 0, nullptr, 16, AUTO, "AFAreaMode", &saAFAreaMode},
    {0, AC_WRITE, 0, nullptr, 17, AUTO, "AFPointSelected", &saAFPointSelected2},
    {0, AC_WRITE, 0, nullptr, 19, AUTO, "MeteringMode", &saMeteringMode1_4},
    {0, AC_WRITE, 0, nullptr, 20, AUTO, "ISOSetting", &saISOSettingInterpreter},
    {0, AC_WRITE, 0, nullptr, 22, AUTO, "DynamicRangeOptimizerMode", &saDynamicRangeOptimizerMode},
    {0, AC_WRITE, 0, nullptr, 23, AUTO, "DynamicRangeOptimizerLevel", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 24, AUTO, "CreativeStyle", &saCreativeStyle2},
    {0, AC_WRITE, 0, nullptr, 25, AUTO, "Sharpness", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 26, AUTO, "Contrast", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 27, AUTO, "Saturation", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 35, AUTO, "FlashMode", &saFlashMode},
    {0, AC_WRITE, 0, nullptr, 38, AUTO, "HighISONoiseReduction", &saHighISONoiseReduction4},
    {0, AC_WRITE, 0, nullptr, 60, AUTO, "ExposureProgram", &saExposureProgram},
    {0, AC_WRITE, 0, nullptr, 63, AUTO, "Rotation", &saRotation},
    {0, AC_WRITE, 0, nullptr, 83, AUTO, "FocusStatus", &saFocusStatusInterpreter},
    {0, AC_WRITE, 0, nullptr, 84, AUTO, "SonyImageSize", &saSonyImageSize},
    {0, AC_WRITE, 0, nullptr, 85, AUTO, "AspectRatio", &saAspectRatio},
    {0, AC_WRITE, 0, nullptr, 86, AUTO, "Quality", &saQualityInterpreter2},
    {0, AC_WRITE, 0, nullptr, 88, AUTO, "ExposureLevelIncrements", &saExposureLevelIncrements},
    {0, AC_WRITE, 0, nullptr, 126, AUTO, "DriveMode", &saDriveMode2},
    {0, AC_WRITE, 0, nullptr, 131, AUTO, "ColorSpace", &saColorSpace5_6},
    { -1, AC_DONTWRITE, 0,  nullptr, 0, AUTO, "", nullptr}
};

const TagAttrib sonyCameraSettingsAttribs3[] = {
    {0, AC_WRITE, 0, nullptr,  0, AUTO, "ShutterSpeedSetting", &saExposureTimeInterpreter},
    {0, AC_WRITE, 0, nullptr,  1, AUTO, "ApertureSetting", &saFNumberInterpreter},
    {0, AC_WRITE, 0, nullptr,  2, AUTO, "ISOSetting", &saISOSettingInterpreter},
    {0, AC_WRITE, 0, nullptr,  3, AUTO, "ExposureCompensationSet", &saExposureCompSetInterpreter},
    {0, AC_WRITE, 0, nullptr,  3, AUTO, "DriveModeSetting", &saDriveMode3},
    {0, AC_WRITE, 0, nullptr,  5, AUTO, "ExposureProgram", &saExposureProgram2},
    {0, AC_WRITE, 0, nullptr,  6, AUTO, "FocusModeSetting", &saFocusModeSetting3},
    {0, AC_WRITE, 0, nullptr,  7, AUTO, "MeteringMode", &saMeteringMode1_3},
    {0, AC_WRITE, 0, nullptr,  9, AUTO, "SonyImageSize", &saSonyImageSize3},
    {0, AC_WRITE, 0, nullptr, 10, AUTO, "AspectRatio", &saAspectRatio2},
    {0, AC_WRITE, 0, nullptr, 11, AUTO, "Quality", &saQualityInterpreter3},
    {0, AC_WRITE, 0, nullptr, 12, AUTO, "DynamicRangeOptimizerSetting", &saDynamicRangeOptimizerSetting},
    {0, AC_WRITE, 0, nullptr, 14, AUTO, "ColorSpace", &saColorSpace1_2},
    {0, AC_WRITE, 0, nullptr, 15, AUTO, "CreativeStyleSetting", &saCreativeStyleSetting},
    {0, AC_WRITE, 0, nullptr, 16, AUTO, "Contrast", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 17, AUTO, "Saturation", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 18, AUTO, "Sharpness", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 22, AUTO, "WhiteBalance", &saWhiteBalanceSettingInterpreter},
    {0, AC_WRITE, 0, nullptr, 23, AUTO, "ColorTemperatureSetting", &saColorTemperatureSettingInterpreter},
    {0, AC_WRITE, 0, nullptr, 23, AUTO, "ColorCompensationFilterSet", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 32, AUTO, "FlashMode", &saFlashMode2},
    {0, AC_WRITE, 0, nullptr, 33, AUTO, "FlashControl", &saFlashControl},
    {0, AC_WRITE, 0, nullptr, 35, AUTO, "FlashExposureCompSet", &saExposureCompSetInterpreter},
    {0, AC_WRITE, 0, nullptr, 36, AUTO, "AFAreaMode", &saAFAreaMode2},
    {0, AC_WRITE, 0, nullptr, 37, AUTO, "LongExposureNoiseReduction", &saOnOffInterpreter2},
    {0, AC_WRITE, 0, nullptr, 38, AUTO, "HighISONoiseReduction", &saHighISONoiseReduction5},
    {0, AC_WRITE, 0, nullptr, 39, AUTO, "SmileShutterMode", &saSmileShutterMode},
    {0, AC_WRITE, 0, nullptr, 40, AUTO, "RedEyeReduction", &saOnOffInterpreter2},
    {0, AC_WRITE, 0, nullptr, 45, AUTO, "HDRSetting", &saOnOffInterpreter3},
    {0, AC_WRITE, 0, nullptr, 46, AUTO, "HDRLevel", &saHDRLevel},
    {0, AC_WRITE, 0, nullptr, 47, AUTO, "ViewingMode", &saViewingMode},
    {0, AC_WRITE, 0, nullptr, 48, AUTO, "FaceDetection", &saOnOffInterpreter2},
    {0, AC_WRITE, 0, nullptr, 49, AUTO, "SmileShutter", &saOnOffInterpreter2},
    {0, AC_WRITE, 0, nullptr, 50, AUTO, "SweepPanoramaSize", &saSweepPanoramaSize},
    {0, AC_WRITE, 0, nullptr, 51, AUTO, "SweepPanoramaDirection", &saSweepPanoramaDirection},
    {0, AC_WRITE, 0, nullptr, 52, AUTO, "DriveMode", &saDriveMode3},
    {0, AC_WRITE, 0, nullptr, 53, AUTO, "MultiFrameNoiseReduction", &saOnOffInterpreter4},
    {0, AC_WRITE, 0, nullptr, 54, AUTO, "LiveViewAFSetting", &saLiveViewAFSetting},
    {0, AC_WRITE, 0, nullptr, 56, AUTO, "PanoramaSize3D", &saPanoramaSize3D},
    {0, AC_WRITE, 0, nullptr, 131, AUTO, "AFButtonPressed", &saNoYesInterpreter},
    {0, AC_WRITE, 0, nullptr, 132, AUTO, "LiveViewMetering", &saLiveViewMetering},
    {0, AC_WRITE, 0, nullptr, 133, AUTO, "ViewingMode2", &saViewingMode},
    {0, AC_WRITE, 0, nullptr, 134, AUTO, "AELock", &saOnOffInterpreter5},
    {0, AC_WRITE, 0, nullptr, 135, AUTO, "FlashAction", &saFlashAction},
    {0, AC_WRITE, 0, nullptr, 139, AUTO, "LiveViewFocusMode", &saLiveViewFocusMode},
    {0, AC_WRITE, 0, nullptr, 153, AUTO, "LensMount", &saLensMount},
    {0, AC_WRITE, 0, nullptr, 643, AUTO, "AFButtonPressed", &saNoYesInterpreter},
    {0, AC_WRITE, 0, nullptr, 644, AUTO, "LiveViewMetering", &saLiveViewMetering},
    {0, AC_WRITE, 0, nullptr, 645, AUTO, "ViewingMode2", &saViewingMode},
    {0, AC_WRITE, 0, nullptr, 646, AUTO, "AELock", &saOnOffInterpreter5},
    {0, AC_WRITE, 0, nullptr, 647, AUTO, "FlashAction", &saFlashAction},
    {0, AC_WRITE, 0, nullptr, 651, AUTO, "LiveViewFocusMode", &saLiveViewFocusMode},
    {0, AC_WRITE, 0, nullptr, 1015, SHORT, "LensType2", &saLensID2Interpreter},
    { -1, AC_DONTWRITE, 0,  nullptr, 0, AUTO, "", nullptr}
};

/*const TagAttrib sonyDNGMakerNote[]={
 {0, AC_WRITE, 0, 0, 0x7200, AUTO, "SonyOffset", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x7201, AUTO, "SonyLength", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x7221, AUTO, "SonyKey", &stdInterpreter},
 {-1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}};*/

}
#endif


