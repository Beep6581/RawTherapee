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
#ifndef _PENTAXATTRIBS_
#define _PENTAXATTRIBS_

#include <cmath>
#include <cstdio>
#include <cstring> /* memcpy() */
#include <string>
#include <sstream>

#include "rtexif.h"

namespace rtexif
{


class PAQualityInterpreter : public ChoiceInterpreter
{
public:
    PAQualityInterpreter ()
    {
        choices[0] = "Good";
        choices[1] = "Better";
        choices[2] = "Best";
        choices[3] = "TIFF";
        choices[4] = "RAW";
        choices[5] = "Premium";
        choices[7] = "RAW (pixel shift enabled)";
        choices[65535] = "n/a";
    }
};
PAQualityInterpreter paQualityInterpreter;

class PAOnOffInterpreter : public ChoiceInterpreter
{
public:
    PAOnOffInterpreter ()
    {
        choices[0]      = "Off";
        choices[1]      = "On";
    }
};
PAOnOffInterpreter paOnOffInterpreter;

class PAShakeReductionInterpreter : public ChoiceInterpreter
{
public:
    PAShakeReductionInterpreter ()
    {
        choices[  0] = "Off";
        choices[  1] = "On";
        choices[  4] = "On (4)";
        choices[  5] = "On but Disabled";
        choices[  6] = "On (Video)";
        choices[  7] = "On (7)";
        choices[ 15] = "On (15)";
        choices[ 39] = "On (mode 2)";
        choices[135] = "On (135)";
        choices[167] = "On (mode 1)";
    }
};
PAShakeReductionInterpreter paShakeReductionInterpreter;

class PAShakeReduction2Interpreter : public ChoiceInterpreter
{
public:
    // ShakeReduction
    PAShakeReduction2Interpreter ()
    {
        choices[ 0] = "Off";
        choices[ 1] = "On";
        choices[ 4] = "Off (AA simulation off)";
        choices[ 5] = "On but Disabled";
        choices[ 6] = "On (Video)";
        choices[ 7] = "On (AA simulation off)";
        choices[12] = "Off (AA simulation type 1)";
        choices[15] = "On (AA simulation type 1)";
        choices[20] = "Off (AA simulation type 2)";
        choices[23] = "On (AA simulation type 2)";
    }
};
PAShakeReduction2Interpreter paShakeReduction2Interpreter;

class PAPictureModeInterpreter : public ChoiceInterpreter
{
public:
    PAPictureModeInterpreter ()
    {
        choices[0] = "Program";
        choices[1] = "Shutter Speed Priority";
        choices[2] = "Program AE";
        choices[3] = "Manual";
        choices[5] = "Portrait";
        choices[6] = "Landscape";
        choices[8] = "Sport";
        choices[9] = "Night Scene";
        choices[11] = "Soft";
        choices[12] = "Surf & Snow";
        choices[13] = "Candlelight";
        choices[14] = "Autumn";
        choices[15] = "Macro";
        choices[17] = "Fireworks";
        choices[18] = "Text";
        choices[19] = "Panorama";
        choices[20] = "3-D";
        choices[21] = "Black & White";
        choices[22] = "Sepia";
        choices[23] = "Red";
        choices[24] = "Pink";
        choices[25] = "Purple";
        choices[26] = "Blue";
        choices[27] = "Green";
        choices[28] = "Yellow";
        choices[30] = "Self Portrait";
        choices[31] = "Illustrations";
        choices[33] = "Digital Filter";
        choices[35] = "Night Scene Portrait";
        choices[37] = "Museum";
        choices[38] = "Food";
        choices[39] = "Underwater";
        choices[40] = "Green Mode";
        choices[49] = "Light Pet";
        choices[50] = "Dark Pet";
        choices[51] = "Medium Pet";
        choices[53] = "Underwater";
        choices[54] = "Candlelight";
        choices[55] = "Natural Skin Tone";
        choices[56] = "Synchro Sound Record";
        choices[58] = "Frame Composite";
        choices[59] = "Report";
        choices[60] = "Kids";
        choices[61] = "Blur Reduction";
        choices[63] = "Panorama 2";
        choices[65] = "Half-length Portrait";
        choices[66] = "Portrait 2";
        choices[74] = "Digital Microscope";
        choices[75] = "Blue Sky";
        choices[80] = "Miniature";
        choices[81] = "HDR";
        choices[83] = "Fisheye";
        choices[85] = "Digital Filter 4";
        choices[221] = "P";
        choices[255] = "PICT";
    }
};
PAPictureModeInterpreter paPictureModeInterpreter;

class PASceneModeInterpreter : public ChoiceInterpreter
{
public:
    PASceneModeInterpreter ()
    {
        choices[0] = "Off";
        choices[1] = "HDR";
        choices[4] = "Auto PICT";
        choices[5] = "Portrait";
        choices[6] = "Landscape";
        choices[7] = "Macro";
        choices[8] = "Sport";
        choices[9] = "Night Scene Portrait";
        choices[10] = "No Flash";
        choices[11] = "Night Scene";
        choices[12] = "Surf & Snow";
        choices[14] = "Sunset";
        choices[15] = "Kids";
        choices[16] = "Pet";
        choices[17] = "Candlelight";
        choices[18] = "Museum";
        choices[20] = "Food";
        choices[21] = "Stage Lighting";
        choices[22] = "Night Snap";
        choices[25] = "Night Scene HDR";
        choices[26] = "Blue Sky";
        choices[27] = "Forest";
        choices[29] = "Backlight Silhouette";
    }
};
PASceneModeInterpreter paSceneModeInterpreter;

class PAAEProgramModeInterpreter : public ChoiceInterpreter
{
public:
    PAAEProgramModeInterpreter ()
    {
        choices[0] = "M, P or TAv";
        choices[1] = "Av, B or X";
        choices[2] = "Tv";
        choices[3] = "Sv or Green Mode";
        choices[8] = "Hi-speed Program";
        choices[11] = "Hi-speed Program (P-Shift)";
        choices[16] = "DOF Program";
        choices[19] = "DOF Program (P-Shift)";
        choices[24] = "MTF Program";
        choices[27] = "MTF Program (P-Shift)";
        choices[35] = "Standard";
        choices[43] = "Portrait";
        choices[51] = "Landscape";
        choices[59] = "Macro";
        choices[67] = "Sport";
        choices[75] = "Night Scene Portrait";
        choices[83] = "No Flash";
        choices[91] = "Night Scene";
        choices[99] = "Surf & Snow";
        choices[104] = "Night Snap";
        choices[107] = "Text";
        choices[115] = "Sunset";
        choices[123] = "Kids";
        choices[131] = "Pet";
        choices[139] = "Candlelight";
        choices[144] = "SCN";
        choices[147] = "Museum";
        choices[160] = "Program";
        choices[184] = "Shallow DOF Program";
        choices[216] = "HDR";
    }
};
PAAEProgramModeInterpreter paAEProgramModeInterpreter;

class PAFlashModeInterpreter : public ChoiceInterpreter
{
public:
    PAFlashModeInterpreter ()
    {
        choices[0] = "Auto, Did not fire";
        choices[1] = "Off, Did not fire";
        choices[2] = "On, Did not fire";
        choices[3] = "Auto, Did not fire, Red-eye reduction";
        choices[5] = "On, Did not fire, Wireless (Master)";
        choices[256] = "Auto, Fired";
        choices[258] = "On, Fired";
        choices[259] = "Auto, Fired, Red-eye reduction";
        choices[260] = "On, Red-eye reduction";
        choices[261] = "On, Wireless (Master)";
        choices[262] = "On, Wireless (Control)";
        choices[264] = "On, Soft";
        choices[265] = "On, Slow-sync";
        choices[266] = "On, Slow-sync, Red-eye reduction";
        choices[267] = "On, Trailing-curtain Sync";
    }
};
PAFlashModeInterpreter paFlashModeInterpreter;

class PAFocusModeInterpreter : public ChoiceInterpreter
{
public:
    PAFocusModeInterpreter ()
    {
        choices[0] = "Normal";
        choices[1] = "Macro";
        choices[2] = "Infinity";
        choices[3] = "Manual";
        choices[4] = "Super Macro";
        choices[5] = "Pan Focus";
        choices[16] = "AF-S (Focus-priority)";
        choices[17] = "AF-C (Focus-priority)";
        choices[18] = "AF-A (Focus-priority)";
        choices[32] = "Contrast-detect (Focus-priority)";
        choices[33] = "Tracking Contrast-detect (Focus-priority)";
        choices[272] = "AF-S (Release-priority)";
        choices[273] = "AF-C (Release-priority)";
        choices[274] = "AF-A (Release-priority)";
        choices[288] = "Contrast-detect (Release-priority)";
    }
};
PAFocusModeInterpreter paFocusModeInterpreter;

class PAAFPointInterpreter : public ChoiceInterpreter
{
public:
    // AFPointSelected
    PAAFPointInterpreter        ()
    {
        choices[0] = "None";
        choices[1] = "Upper-left";
        choices[2] = "Top";
        choices[3] = "Upper-right";
        choices[4] = "Left";
        choices[5] = "Mid-left";
        choices[6] = "Center";
        choices[7] = "Mid-right";
        choices[8] = "Right";
        choices[9] = "Lower-left";
        choices[10] = "Bottom";
        choices[11] = "Lower-right";
        choices[65531] = "AF Select";
        choices[65532] = "Face Detect AF";
        choices[65533] = "Automatic Tracking AF";
        choices[65534] = "Fixed Center";
        choices[65535] = "Auto";
    }
};
PAAFPointInterpreter paAFPointInterpreter;

class PAAFFocusInterpreter : public ChoiceInterpreter
{
public:
    // AFPointsInFocus
    PAAFFocusInterpreter        ()
    {
        choices[0] = "Fixed Center or Multiple";
        choices[1] = "Top-left";
        choices[2] = "Top-center";
        choices[3] = "Top-right";
        choices[4] = "Left";
        choices[5] = "Center";
        choices[6] = "Right";
        choices[7] = "Bottom-left";
        choices[8] = "Bottom-center";
        choices[9] = "Bottom-right";
        choices[65535] = "None";
    }
};
PAAFFocusInterpreter paAFFocusInterpreter;

class PAISOInterpreter : public ChoiceInterpreter
{
public:
    PAISOInterpreter        ()
    {
        choices[3] = "50";
        choices[4] = "64";
        choices[5] = "80";
        choices[6] = "100";
        choices[7] = "125";
        choices[8] = "160";
        choices[9] = "200";
        choices[10] = "250";
        choices[11] = "320";
        choices[12] = "400";
        choices[13] = "500";
        choices[14] = "640";
        choices[15] = "800";
        choices[16] = "1000";
        choices[17] = "1250";
        choices[18] = "1600";
        choices[19] = "2000";
        choices[20] = "2500";
        choices[21] = "3200";
        choices[22] = "4000";
        choices[23] = "5000";
        choices[24] = "6400";
        choices[25] = "8000";
        choices[26] = "10000";
        choices[27] = "12800";
        choices[28] = "16000";
        choices[29] = "20000";
        choices[30] = "25600";
        choices[31] = "32000";
        choices[32] = "40000";
        choices[33] = "51200";
        choices[34] = "64000";
        choices[35] = "80000";
        choices[36] = "102400";
        choices[37] = "128000";
        choices[38] = "160000";
        choices[39] = "204800";
        choices[50] = "50";
        choices[100] = "100";
        choices[200] = "200";
        choices[258] = "50";
        choices[259] = "70";
        choices[260] = "100";
        choices[261] = "140";
        choices[262] = "200";
        choices[263] = "280";
        choices[264] = "400";
        choices[265] = "560";
        choices[266] = "800";
        choices[267] = "1100";
        choices[268] = "1600";
        choices[269] = "2200";
        choices[270] = "3200";
        choices[271] = "4500";
        choices[272] = "6400";
        choices[273] = "9000";
        choices[274] = "12800";
        choices[275] = "18000";
        choices[276] = "25600";
        choices[277] = "36000";
        choices[278] = "51200";
        choices[400] = "400";
        choices[800] = "800";
        choices[1600] = "1600";
        choices[3200] = "3200";
    }
};
PAISOInterpreter paISOInterpreter;

class PAFNumberInterpreter: public Interpreter
{
public:
    PAFNumberInterpreter () {}
    virtual std::string toString (Tag* t)
    {
        char buffer[32];
        double v = t->toDouble() / 10;

        if( v < 0. || v > 1000. ) {
            return "undef";
        }

        sprintf (buffer, "%.1f", v );
        return buffer;
    }
};
PAFNumberInterpreter paFNumberInterpreter;

class PAMeteringModeInterpreter : public ChoiceInterpreter
{
public:
    PAMeteringModeInterpreter ()
    {
        choices[0] = "Multi-segment";
        choices[1] = "Center-weighted average";
        choices[2] = "Spot";
    }
};
PAMeteringModeInterpreter paMeteringModeInterpreter;

class PAWhiteBalanceInterpreter : public ChoiceInterpreter
{
public:
    PAWhiteBalanceInterpreter ()
    {
        choices[0] = "Auto";
        choices[1] = "Daylight";
        choices[2] = "Shade";
        choices[3] = "Fluorescent";
        choices[4] = "Tungsten";
        choices[5] = "Manual";
        choices[6] = "Daylight Fluorescent";
        choices[7] = "Day White Fluorescent";
        choices[8] = "White Fluorescent";
        choices[9] = "Flash";
        choices[10] = "Cloudy";
        choices[11] = "Warm White Fluorescent";
        choices[14] = "Multi Auto";
        choices[15] = "Color Temperature Enhancement";
        choices[17] = "Kelvin";
        choices[65534] = "Unknown";
        choices[65535] = "User-Selected";
    }
};
PAWhiteBalanceInterpreter paWhiteBalanceInterpreter;

class PAWhiteBalanceModeInterpreter : public ChoiceInterpreter
{
public:
    PAWhiteBalanceModeInterpreter ()
    {
        choices[1] = "Auto (Daylight)";
        choices[2] = "Auto (Shade)";
        choices[3] = "Auto (Flash)";
        choices[4] = "Auto (Tungsten)";
        choices[6] = "Auto (Daylight Fluorescent)";
        choices[7] = "Auto (Day White Fluorescent)";
        choices[8] = "Auto (White Fluorescent)";
        choices[10] = "Auto (Cloudy)";
        choices[65534] = "Unknown";
        choices[65535] = "User-Selected";
    }
};
PAWhiteBalanceModeInterpreter paWhiteBalanceModeInterpreter;

class PASaturationInterpreter : public ChoiceInterpreter
{
public:
    PASaturationInterpreter ()
    {
        choices[0] = "-2 (low)";
        choices[1] = "0 (normal)";
        choices[2] = "+2 (high)";
        choices[3] = "-1 (med low)";
        choices[4] = "+1 (med high)";
        choices[5] = "-3 (very low)";
        choices[6] = "+3 (very high)";
        choices[7] = "-4 (minimum)";
        choices[8] = "+4 (maximum)";
        choices[65535] = "None";
    }
};
PASaturationInterpreter paSaturationInterpreter;

class PAContrastInterpreter : public ChoiceInterpreter
{
public:
    PAContrastInterpreter ()
    {
        choices[0] = "-2 (low)";
        choices[1] = "0 (normal)";
        choices[2] = "+2 (high)";
        choices[3] = "-1 (med low)";
        choices[4] = "+1 (med high)";
        choices[5] = "-3 (very low)";
        choices[6] = "+3 (very high)";
        choices[7] = "-4 (minimum)";
        choices[8] = "+4 (maximum)";
        choices[65535] = "n/a";
    }
};
PAContrastInterpreter paContrastInterpreter;

class PASharpnessInterpreter : public ChoiceInterpreter
{
public:
    PASharpnessInterpreter ()
    {
        choices[0] = "-2 (soft)";
        choices[1] = "0 (normal)";
        choices[2] = "+2 (hard)";
        choices[3] = "-1 (med soft)";
        choices[4] = "+1 (med hard)";
        choices[5] = "-3 (very soft)";
        choices[6] = "+3 (very hard)";
        choices[7] = "-4 (minimum)";
        choices[8] = "+4 (maximum)";
    }
};
PASharpnessInterpreter paSharpnessInterpreter;

class PAPictureModeInterpreter2: public ChoiceInterpreter
{
public:
    PAPictureModeInterpreter2()
    {
        choices[256 *   0 +   0] = "Program";
        choices[256 *   0 +   1] = "Hi-speed Program";
        choices[256 *   0 +   2] = "DOF Program";
        choices[256 *   0 +   3] = "MTF Program";
        choices[256 *   0 +   4] = "Standard";
        choices[256 *   0 +   5] = "Portrait";
        choices[256 *   0 +   6] = "Landscape";
        choices[256 *   0 +   7] = "Macro";
        choices[256 *   0 +   8] = "Sport";
        choices[256 *   0 +   9] = "Night Scene Portrait";
        choices[256 *   0 +  10] = "No Flash";
        choices[256 *   0 +  11] = "Night Scene";
        choices[256 *   0 +  12] = "Surf & Snow";
        choices[256 *   0 +  13] = "Text";
        choices[256 *   0 +  14] = "Sunset";
        choices[256 *   0 +  15] = "Kids";
        choices[256 *   0 +  16] = "Pet";
        choices[256 *   0 +  17] = "Candlelight";
        choices[256 *   0 +  18] = "Museum";
        choices[256 *   0 +  19] = "Food";
        choices[256 *   0 +  20] = "Stage Lighting";
        choices[256 *   0 +  21] = "Night Snap";
        choices[256 *   0 +  23] = "Blue Sky";
        choices[256 *   0 +  24] = "Sunset";
        choices[256 *   0 +  26] = "Night Scene HDR";
        choices[256 *   0 +  27] = "HDR";
        choices[256 *   0 +  28] = "Quick Macro";
        choices[256 *   0 +  29] = "Forest";
        choices[256 *   0 +  30] = "Backlight Silhouette";
        choices[256 *   1 +   4] = "Auto PICT (Standard)";
        choices[256 *   1 +   5] = "Auto PICT (Portrait)";
        choices[256 *   1 +   6] = "Auto PICT (Landscape)";
        choices[256 *   1 +   7] = "Auto PICT (Macro)";
        choices[256 *   1 +   8] = "Auto PICT (Sport)";
        choices[256 *   2 +   0] = "Program (HyP)";
        choices[256 *   2 +   1] = "Hi-speed Program (HyP)";
        choices[256 *   2 +   2] = "DOF Program (HyP)";
        choices[256 *   2 +   3] = "MTF Program (HyP)";
        choices[256 *   2 +  22] = "Shallow DOF (HyP)";
        choices[256 *   3 +   0] = "Green Mode";
        choices[256 *   4 +   0] = "Shutter Speed Priority";
        choices[256 *   5 +   0] = "Aperture Priority";
        choices[256 *   6 +   0] = "Program Tv Shift";
        choices[256 *   7 +   0] = "Program Av Shift";
        choices[256 *   8 +   0] = "Manual";
        choices[256 *   9 +   0] = "Bulb";
        choices[256 *  10 +   0] = "Aperture Priority, Off-Auto-Aperture";
        choices[256 *  11 +   0] = "Manual, Off-Auto-Aperture";
        choices[256 *  12 +   0] = "Bulb, Off-Auto-Aperture";
        choices[256 *  13 +   0] = "Shutter & Aperture Priority AE";
        choices[256 *  15 +   0] = "Sensitivity Priority AE";
        choices[256 *  16 +   0] = "Flash X-Sync Speed AE";
        choices[256 *  18 +   0] = "Auto Program (Normal)";
        choices[256 *  18 +   1] = "Auto Program (Hi-speed)";
        choices[256 *  18 +   2] = "Auto Program (DOF)";
        choices[256 *  18 +   3] = "Auto Program (MTF)";
        choices[256 *  18 +  22] = "Auto Program (Shallow DOF)";
        choices[256 *  20 +  22] = "Blur Control";
        choices[256 * 254 +   0] = "Video";
        choices[256 * 255 +   0] = "Video (Auto Aperture)";
        choices[256 * 255 +   4] = "Video (4)";
    }
    virtual std::string toString (Tag* t)
    {
        int c = 256 * t->toInt(0, BYTE) + t->toInt(1, BYTE);
        std::map<int, std::string>::iterator r = choices.find (c);

        if (r != choices.end()) {
            std::ostringstream s;
            s << r->second;

            if( t->toInt(1, BYTE) == 0 ) {
                s << "\n1/2 EV steps";
            } else {
                s << "\n1/3 EV steps";
            }

            return s.str();
        } else {
            char buffer[1024];
            t->toString (buffer);
            return std::string (buffer);
        }
    }
};
PAPictureModeInterpreter2 paPictureModeInterpreter2;

class PADriveModeInterpreter : public ChoiceInterpreter
{
    std::map<int, std::string> choices1;
    std::map<int, std::string> choices2;
    std::map<int, std::string> choices3;
public:
    PADriveModeInterpreter()
    {
        choices[0]    = "Single-frame";
        choices[1]    = "Continuous";
        choices[2]    = "Continuous (Lo)";
        choices[3]    = "Burst";
        choices[4]    = "Continuous (Medium)";
        choices[255]  = "Video";
        choices1[0]   = "No Timer";
        choices1[1]   = "Self-timer (12 s)";
        choices1[2]   = "Self-timer (2 s)";
        choices1[15]  = "Video";
        choices1[16]  = "Mirror Lock-up";
        choices1[255] = "n/a";
        choices2[0]   = "Shutter Button";
        choices2[1]   = "Remote Control (3 s delay)";
        choices2[2]   = "Remote Control";
        choices2[4]   = "Remote Continuous Shooting";
        choices3[0]   = "Single Exposure";
        choices3[1]   = "Multiple Exposure";
        choices3[15]  = "Interval Movie";
        choices3[16]  = "HDR";
        choices3[32]  = "HDR Strong 1";
        choices3[48]  = "HDR Strong 2";
        choices3[64]  = "HDR Strong 3";
        choices3[224] = "HDR Auto";
        choices3[255] = "Video";
    }
    virtual std::string toString (Tag* t)
    {
        std::map<int, std::string>::iterator r  = choices.find (t->toInt(0, BYTE));
        std::map<int, std::string>::iterator r1 = choices1.find (t->toInt(1, BYTE));
        std::map<int, std::string>::iterator r2 = choices2.find (t->toInt(2, BYTE));
        std::map<int, std::string>::iterator r3 = choices3.find (t->toInt(3, BYTE));
        std::ostringstream s;
        s << ((r != choices.end()) ? r->second : "");
        s << ((r1 != choices1.end()) ? r1->second : "") << " ";
        s << ((r2 != choices2.end()) ? r2->second : "") << " ";
        s << ((r3 != choices3.end()) ? r3->second : "") << " ";
        return s.str();
    }
};
PADriveModeInterpreter paDriveModeInterpreter;

class PAColorSpaceInterpreter: public ChoiceInterpreter
{
public:
    PAColorSpaceInterpreter()
    {
        choices[0] = "sRGB";
        choices[1] = "Adobe RGB";
    }
};
PAColorSpaceInterpreter paColorSpaceInterpreter;

class PALensTypeInterpreter : public IntLensInterpreter< int >
{
public:
    PALensTypeInterpreter ()
    {
        choices.insert(p_t(256 * 0 + 0, "M-42 or No Lens"));
        choices.insert(p_t(256 * 1 + 0, "K or M Lens"));
        choices.insert(p_t(256 * 2 + 0, "A Series Lens"));
        choices.insert(p_t(256 * 3 + 0, "Sigma"));
        choices.insert(p_t(256 * 3 + 17, "smc PENTAX-FA SOFT 85mm f/2.8"));
        choices.insert(p_t(256 * 3 + 18, "smc PENTAX-F 1.7X AF ADAPTER"));
        choices.insert(p_t(256 * 3 + 19, "smc PENTAX-F 24-50mm f/4"));
        choices.insert(p_t(256 * 3 + 20, "smc PENTAX-F 35-80mm f/4-5.6"));
        choices.insert(p_t(256 * 3 + 21, "smc PENTAX-F 80-200mm f/4.7-5.6"));
        choices.insert(p_t(256 * 3 + 22, "smc PENTAX-F FISH-EYE 17-28mm f/3.5-4.5"));
        choices.insert(p_t(256 * 3 + 23, "smc PENTAX-F 100-300mm f/4.5-5.6 or Sigma Lens"));
        choices.insert(p_t(256 * 3 + 23, "Sigma AF 28-300mm f/3.5-5.6 DL IF"));
        choices.insert(p_t(256 * 3 + 23, "Sigma AF 28-300mm f/3.5-6.3 DG IF Macro"));
        choices.insert(p_t(256 * 3 + 23, "Tokina 80-200mm f/2.8 ATX-Pro"));
        choices.insert(p_t(256 * 3 + 24, "smc PENTAX-F 35-135mm f/3.5-4.5"));
        choices.insert(p_t(256 * 3 + 25, "smc PENTAX-F 35-105mm f/4-5.6 or Sigma or Tokina Lens"));
        choices.insert(p_t(256 * 3 + 25, "Sigma AF 28-300mm f/3.5-5.6 DL IF"));
        choices.insert(p_t(256 * 3 + 25, "Sigma 55-200mm f/4-5.6 DC"));
        choices.insert(p_t(256 * 3 + 25, "Sigma AF 28-300mm f/3.5-6.3 DL IF"));
        choices.insert(p_t(256 * 3 + 25, "Sigma AF 28-300mm f/3.5-6.3 DG IF Macro"));
        choices.insert(p_t(256 * 3 + 25, "Tokina 80-200mm f/2.8 ATX-Pro"));
        choices.insert(p_t(256 * 3 + 26, "smc PENTAX-F* 250-600mm f/5.6 ED[IF]"));
        choices.insert(p_t(256 * 3 + 27, "smc PENTAX-F 28-80mm f/3.5-4.5 or Tokina Lens"));
        choices.insert(p_t(256 * 3 + 27, "Tokina AT-X Pro AF 28-70mm f/2.6-2.8"));
        choices.insert(p_t(256 * 3 + 28, "smc PENTAX-F 35-70mm f/3.5-4.5 or Tokina Lens"));
        choices.insert(p_t(256 * 3 + 28, "Tokina 19-35mm f/3.5-4.5 AF"));
        choices.insert(p_t(256 * 3 + 28, "Tokina AT-X AF 400mm f/5.6"));
        choices.insert(p_t(256 * 3 + 29, "PENTAX-F 28-80mm f/3.5-4.5 or Sigma or Tokina Lens"));
        choices.insert(p_t(256 * 3 + 29, "Sigma AF 18-125mm f/3.5-5.6 DC"));
        choices.insert(p_t(256 * 3 + 29, "Tokina AT-X PRO 28-70mm f/2.6-2.8"));
        choices.insert(p_t(256 * 3 + 30, "PENTAX-F 70-200mm f/4-5.6"));
        choices.insert(p_t(256 * 3 + 31, "smc PENTAX-F 70-210mm f/4-5.6 or Tokina or Takumar Lens"));
        choices.insert(p_t(256 * 3 + 31, "Tokina AF 730 75-300mm f/4.5-5.6"));
        choices.insert(p_t(256 * 3 + 31, "Takumar-F 70-210mm f/4-5.6"));
        choices.insert(p_t(256 * 3 + 32, "smc PENTAX-F 50mm f/1.4"));
        choices.insert(p_t(256 * 3 + 33, "smc PENTAX-F 50mm f/1.7"));
        choices.insert(p_t(256 * 3 + 34, "smc PENTAX-F 135mm f/2.8 [IF]"));
        choices.insert(p_t(256 * 3 + 35, "smc PENTAX-F 28mm f/2.8"));
        choices.insert(p_t(256 * 3 + 36, "Sigma 20mm f/1.8 EX DG Aspherical RF"));
        choices.insert(p_t(256 * 3 + 38, "smc PENTAX-F* 300mm f/4.5 ED[IF]"));
        choices.insert(p_t(256 * 3 + 39, "smc PENTAX-F* 600mm f/4 ED[IF]"));
        choices.insert(p_t(256 * 3 + 40, "smc PENTAX-F Macro 100mm f/2.8"));
        choices.insert(p_t(256 * 3 + 41, "smc PENTAX-F Macro 50mm f/2.8 or Sigma Lens"));
        choices.insert(p_t(256 * 3 + 41, "Sigma 50mm f/2.8 Macro"));
        choices.insert(p_t(256 * 3 + 42, "Sigma 300mm f/2.8 EX DG APO IF"));
        choices.insert(p_t(256 * 3 + 44, "Sigma or Tamron Lens (3 44)"));
        choices.insert(p_t(256 * 3 + 44, "Sigma AF 10-20mm f/4-5.6 EX DC"));
        choices.insert(p_t(256 * 3 + 44, "Sigma 12-24mm f/4.5-5.6 EX DG"));
        choices.insert(p_t(256 * 3 + 44, "Sigma 17-70mm f/2.8-4.5 DC Macro"));
        choices.insert(p_t(256 * 3 + 44, "Sigma 18-50mm f/3.5-5.6 DC"));
        choices.insert(p_t(256 * 3 + 44, "Sigma 17-35mm f/2.8-4 EX DG"));
        choices.insert(p_t(256 * 3 + 44, "Tamron 35-90mm f/4 AF"));
        choices.insert(p_t(256 * 3 + 46, "Sigma or Samsung Lens (3 46)"));
        choices.insert(p_t(256 * 3 + 46, "Sigma APO 70-200mm f/2.8 EX"));
        choices.insert(p_t(256 * 3 + 46, "Sigma EX APO 100-300mm f/4 IF"));
        choices.insert(p_t(256 * 3 + 46, "Samsung/Schneider D-XENON 50-200mm f/4-5.6 ED"));
        choices.insert(p_t(256 * 3 + 50, "smc PENTAX-FA 28-70mm f/4 AL"));
        choices.insert(p_t(256 * 3 + 51, "Sigma 28mm f/1.8 EX DG Aspherical Macro"));
        choices.insert(p_t(256 * 3 + 52, "smc PENTAX-FA 28-200mm f/3.8-5.6 AL[IF] or Tamron Lens"));
        choices.insert(p_t(256 * 3 + 52, "Tamron AF LD 28-200mm f/3.8-5.6 [IF] Aspherical (171D)"));
        choices.insert(p_t(256 * 3 + 53, "smc PENTAX-FA 28-80mm f/3.5-5.6 AL"));
        choices.insert(p_t(256 * 3 + 247, "smc PENTAX-DA FISH-EYE 10-17mm f/3.5-4.5 ED[IF]"));
        choices.insert(p_t(256 * 3 + 248, "smc PENTAX-DA 12-24mm f/4 ED AL[IF]"));
        choices.insert(p_t(256 * 3 + 250, "smc PENTAX-DA 50-200mm f/4-5.6 ED"));
        choices.insert(p_t(256 * 3 + 251, "smc PENTAX-DA 40mm f/2.8 Limited"));
        choices.insert(p_t(256 * 3 + 252, "smc PENTAX-DA 18-55mm f/3.5-5.6 AL"));
        choices.insert(p_t(256 * 3 + 253, "smc PENTAX-DA 14mm f/2.8 ED[IF]"));
        choices.insert(p_t(256 * 3 + 254, "smc PENTAX-DA 16-45mm f/4 ED AL"));
        choices.insert(p_t(256 * 3 + 255, "Sigma Lens (3 255)"));
        choices.insert(p_t(256 * 3 + 255, "Sigma 18-200mm f/3.5-6.3 DC"));
        choices.insert(p_t(256 * 3 + 255, "Sigma DL-II 35-80mm f/4-5.6"));
        choices.insert(p_t(256 * 3 + 255, "Sigma DL Zoom 75-300mm f/4-5.6"));
        choices.insert(p_t(256 * 3 + 255, "Sigma DF EX Aspherical 28-70mm f/2.8"));
        choices.insert(p_t(256 * 3 + 255, "Sigma AF Tele 400mm f/5.6 Multi-coated"));
        choices.insert(p_t(256 * 3 + 255, "Sigma 24-60mm f/2.8 EX DG"));
        choices.insert(p_t(256 * 3 + 255, "Sigma 70-300mm f/4-5.6 Macro"));
        choices.insert(p_t(256 * 3 + 255, "Sigma 55-200mm f/4-5.6 DC"));
        choices.insert(p_t(256 * 3 + 255, "Sigma 18-50mm f/2.8 EX DC"));
        choices.insert(p_t(256 * 4 + 1, "smc PENTAX-FA SOFT 28mm f/2.8"));
        choices.insert(p_t(256 * 4 + 2, "smc PENTAX-FA 80-320mm f/4.5-5.6"));
        choices.insert(p_t(256 * 4 + 3, "smc PENTAX-FA 43mm f/1.9 Limited"));
        choices.insert(p_t(256 * 4 + 6, "smc PENTAX-FA 35-80mm f/4-5.6"));
        choices.insert(p_t(256 * 4 + 12, "smc PENTAX-FA 50mm f/1.4"));
        choices.insert(p_t(256 * 4 + 15, "smc PENTAX-FA 28-105mm f/4-5.6 [IF]"));
        choices.insert(p_t(256 * 4 + 16, "Tamron AF 80-210mm f/4-5.6 (178D)"));
        choices.insert(p_t(256 * 4 + 19, "Tamron SP AF 90mm f/2.8 (172E)"));
        choices.insert(p_t(256 * 4 + 20, "smc PENTAX-FA 28-80mm f/3.5-5.6"));
        choices.insert(p_t(256 * 4 + 21, "Cosina AF 100-300mm f/5.6-6.7"));
        choices.insert(p_t(256 * 4 + 22, "Tokina 28-80mm f/3.5-5.6"));
        choices.insert(p_t(256 * 4 + 23, "smc PENTAX-FA 20-35mm f/4 AL"));
        choices.insert(p_t(256 * 4 + 24, "smc PENTAX-FA 77mm f/1.8 Limited"));
        choices.insert(p_t(256 * 4 + 25, "Tamron SP AF 14mm f/2.8"));
        choices.insert(p_t(256 * 4 + 26, "smc PENTAX-FA Macro 100mm f/3.5 or Cosina Lens"));
        choices.insert(p_t(256 * 4 + 26, "Cosina 100mm f/3.5 Macro"));
        choices.insert(p_t(256 * 4 + 27, "Tamron AF 28-300mm f/3.5-6.3 LD Aspherical[IF] Macro (185D/285D)"));
        choices.insert(p_t(256 * 4 + 28, "smc PENTAX-FA 35mm f/2 AL"));
        choices.insert(p_t(256 * 4 + 29, "Tamron AF 28-200mm f/3.8-5.6 LD Super II Macro (371D)"));
        choices.insert(p_t(256 * 4 + 34, "smc PENTAX-FA 24-90mm f/3.5-4.5 AL[IF]"));
        choices.insert(p_t(256 * 4 + 35, "smc PENTAX-FA 100-300mm f/4.7-5.8"));
        choices.insert(p_t(256 * 4 + 36, "Tamron AF 70-300mm f/4-5.6 LD Macro 1:2"));
        choices.insert(p_t(256 * 4 + 37, "Tamron SP AF 24-135mm f/3.5-5.6 AD AL (190D)"));
        choices.insert(p_t(256 * 4 + 38, "smc PENTAX-FA 28-105mm f/3.2-4.5 AL[IF]"));
        choices.insert(p_t(256 * 4 + 39, "smc PENTAX-FA 31mm f/1.8 AL Limited"));
        choices.insert(p_t(256 * 4 + 41, "Tamron AF 28-200mm Super Zoom f/3.8-5.6 Aspherical XR [IF] Macro (A03)"));
        choices.insert(p_t(256 * 4 + 43, "smc PENTAX-FA 28-90mm f/3.5-5.6"));
        choices.insert(p_t(256 * 4 + 44, "smc PENTAX-FA J 75-300mm f/4.5-5.8 AL"));
        choices.insert(p_t(256 * 4 + 45, "Tamron Lens (4 45)"));
        choices.insert(p_t(256 * 4 + 45, "Tamron 28-300mm f/3.5-6.3 Ultra zoom XR"));
        choices.insert(p_t(256 * 4 + 45, "Tamron AF 28-300mm f/3.5-6.3 XR Di LD Aspherical [IF] Macro"));
        choices.insert(p_t(256 * 4 + 46, "smc PENTAX-FA J 28-80mm f/3.5-5.6 AL"));
        choices.insert(p_t(256 * 4 + 47, "smc PENTAX-FA J 18-35mm f/4-5.6 AL"));
        choices.insert(p_t(256 * 4 + 49, "Tamron SP AF 28-75mm f/2.8 XR Di LD Aspherical [IF] Macro"));
        choices.insert(p_t(256 * 4 + 51, "smc PENTAX-D FA 50mm f/2.8 Macro"));
        choices.insert(p_t(256 * 4 + 52, "smc PENTAX-D FA 100mm f/2.8 Macro"));
        choices.insert(p_t(256 * 4 + 55, "Samsung/Schneider D-XENOGON 35mm f/2"));
        choices.insert(p_t(256 * 4 + 56, "Samsung/Schneider D-XENON 100mm f/2.8 Macro"));
        choices.insert(p_t(256 * 4 + 75, "Tamron SP AF 70-200mm f/2.8 Di LD [IF] Macro (A001)"));
        choices.insert(p_t(256 * 4 + 214, "smc PENTAX-DA 35mm f/2.4 AL"));
        choices.insert(p_t(256 * 4 + 229, "smc PENTAX-DA 18-55mm f/3.5-5.6 AL II"));
        choices.insert(p_t(256 * 4 + 230, "Tamron SP AF 17-50mm f/2.8 XR Di II"));
        choices.insert(p_t(256 * 4 + 231, "smc PENTAX-DA 18-250mm f/3.5-6.3 ED AL [IF]"));
        choices.insert(p_t(256 * 4 + 237, "Samsung/Schneider D-XENOGON 10-17mm f/3.5-4.5"));
        choices.insert(p_t(256 * 4 + 239, "Samsung/Schneider D-XENON 12-24mm f/4 ED AL [IF]"));
        choices.insert(p_t(256 * 4 + 242, "smc PENTAX-DA* 16-50mm f/2.8 ED AL [IF] SDM (SDM unused)"));
        choices.insert(p_t(256 * 4 + 243, "smc PENTAX-DA 70mm f/2.4 Limited"));
        choices.insert(p_t(256 * 4 + 244, "smc PENTAX-DA 21mm f/3.2 AL Limited"));
        choices.insert(p_t(256 * 4 + 245, "Samsung/Schneider D-XENON 50-200mm f/4-5.6"));
        choices.insert(p_t(256 * 4 + 246, "Samsung/Schneider D-XENON 18-55mm f/3.5-5.6"));
        choices.insert(p_t(256 * 4 + 247, "smc PENTAX-DA FISH-EYE 10-17mm f/3.5-4.5 ED[IF]"));
        choices.insert(p_t(256 * 4 + 248, "smc PENTAX-DA 12-24mm f/4 ED AL [IF]"));
        choices.insert(p_t(256 * 4 + 249, "Tamron XR DiII 18-200mm f/3.5-6.3 (A14)"));
        choices.insert(p_t(256 * 4 + 250, "smc PENTAX-DA 50-200mm f/4-5.6 ED"));
        choices.insert(p_t(256 * 4 + 251, "smc PENTAX-DA 40mm f/2.8 Limited"));
        choices.insert(p_t(256 * 4 + 252, "smc PENTAX-DA 18-55mm f/3.5-5.6 AL"));
        choices.insert(p_t(256 * 4 + 253, "smc PENTAX-DA 14mm f/2.8 ED[IF]"));
        choices.insert(p_t(256 * 4 + 254, "smc PENTAX-DA 16-45mm f/4 ED AL"));
        choices.insert(p_t(256 * 5 + 1, "smc PENTAX-FA* 24mm f/2 AL[IF]"));
        choices.insert(p_t(256 * 5 + 2, "smc PENTAX-FA 28mm f/2.8 AL"));
        choices.insert(p_t(256 * 5 + 3, "smc PENTAX-FA 50mm f/1.7"));
        choices.insert(p_t(256 * 5 + 4, "smc PENTAX-FA 50mm f/1.4"));
        choices.insert(p_t(256 * 5 + 5, "smc PENTAX-FA* 600mm f/4 ED[IF]"));
        choices.insert(p_t(256 * 5 + 6, "smc PENTAX-FA* 300mm f/4.5 ED[IF]"));
        choices.insert(p_t(256 * 5 + 7, "smc PENTAX-FA 135mm f/2.8 [IF]"));
        choices.insert(p_t(256 * 5 + 8, "smc PENTAX-FA Macro 50mm f/2.8"));
        choices.insert(p_t(256 * 5 + 9, "smc PENTAX-FA Macro 100mm f/2.8"));
        choices.insert(p_t(256 * 5 + 10, "smc PENTAX-FA* 85mm f/1.4 [IF]"));
        choices.insert(p_t(256 * 5 + 11, "smc PENTAX-FA* 200mm f/2.8 ED[IF]"));
        choices.insert(p_t(256 * 5 + 12, "smc PENTAX-FA 28-80mm f/3.5-4.7"));
        choices.insert(p_t(256 * 5 + 13, "smc PENTAX-FA 70-200mm f/4-5.6"));
        choices.insert(p_t(256 * 5 + 14, "smc PENTAX-FA* 250-600mm f/5.6 ED[IF]"));
        choices.insert(p_t(256 * 5 + 15, "smc PENTAX-FA 28-105mm f/4-5.6"));
        choices.insert(p_t(256 * 5 + 16, "smc PENTAX-FA 100-300mm f/4.5-5.6"));
        choices.insert(p_t(256 * 5 + 98, "smc PENTAX-FA 100-300mm f/4.5-5.6"));
        choices.insert(p_t(256 * 6 + 1, "smc PENTAX-FA* 85mm f/1.4 [IF]"));
        choices.insert(p_t(256 * 6 + 2, "smc PENTAX-FA* 200mm f/2.8 ED[IF]"));
        choices.insert(p_t(256 * 6 + 3, "smc PENTAX-FA* 300mm f/2.8 ED[IF]"));
        choices.insert(p_t(256 * 6 + 4, "smc PENTAX-FA* 28-70mm f/2.8 AL"));
        choices.insert(p_t(256 * 6 + 5, "smc PENTAX-FA* 80-200mm f/2.8 ED[IF]"));
        choices.insert(p_t(256 * 6 + 6, "smc PENTAX-FA* 28-70mm f/2.8 AL"));
        choices.insert(p_t(256 * 6 + 7, "smc PENTAX-FA* 80-200mm f/2.8 ED[IF]"));
        choices.insert(p_t(256 * 6 + 8, "smc PENTAX-FA 28-70mm f/4AL"));
        choices.insert(p_t(256 * 6 + 9, "smc PENTAX-FA 20mm f/2.8"));
        choices.insert(p_t(256 * 6 + 10, "smc PENTAX-FA* 400mm f/5.6 ED[IF]"));
        choices.insert(p_t(256 * 6 + 13, "smc PENTAX-FA* 400mm f/5.6 ED[IF]"));
        choices.insert(p_t(256 * 6 + 14, "smc PENTAX-FA* Macro 200mm f/4 ED[IF]"));
        choices.insert(p_t(256 * 7 + 0, "smc PENTAX-DA 21mm f/3.2 AL Limited"));
        choices.insert(p_t(256 * 7 + 58, "smc PENTAX-D FA Macro 100mm f/2.8 WR"));
        choices.insert(p_t(256 * 7 + 75, "Tamron SP AF 70-200mm f/2.8 Di LD [IF] Macro (A001)"));
        choices.insert(p_t(256 * 7 + 201, "smc Pentax-DA L 50-200mm f/4-5.6 ED WR"));
        choices.insert(p_t(256 * 7 + 202, "smc PENTAX-DA L 18-55mm f/3.5-5.6 AL WR"));
        choices.insert(p_t(256 * 7 + 203, "HD PENTAX-DA 55-300mm f/4-5.8 ED WR"));
        choices.insert(p_t(256 * 7 + 204, "HD PENTAX-DA 15mm f/4 ED AL Limited"));
        choices.insert(p_t(256 * 7 + 205, "HD PENTAX-DA 35mm f/2.8 Macro Limited"));
        choices.insert(p_t(256 * 7 + 206, "HD PENTAX-DA 70mm f/2.4 Limited"));
        choices.insert(p_t(256 * 7 + 207, "HD PENTAX-DA 21mm f/3.2 ED AL Limited"));
        choices.insert(p_t(256 * 7 + 208, "HD PENTAX-DA 40mm f/2.8 Limited"));
        choices.insert(p_t(256 * 7 + 212, "smc PENTAX-DA 50mm f/1.8"));
        choices.insert(p_t(256 * 7 + 213, "smc PENTAX-DA 40mm f/2.8 XS"));
        choices.insert(p_t(256 * 7 + 214, "smc PENTAX-DA 35mm f/2.4 AL"));
        choices.insert(p_t(256 * 7 + 216, "smc PENTAX-DA L 55-300mm f/4-5.8 ED"));
        choices.insert(p_t(256 * 7 + 217, "smc PENTAX-DA 50-200mm f/4-5.6 ED WR"));
        choices.insert(p_t(256 * 7 + 218, "smc PENTAX-DA 18-55mm f/3.5-5.6 AL WR"));
        choices.insert(p_t(256 * 7 + 220, "Tamron SP AF 10-24mm f/3.5-4.5 Di II LD Aspherical [IF]"));
        choices.insert(p_t(256 * 7 + 221, "smc PENTAX-DA L 50-200mm f/4-5.6 ED"));
        choices.insert(p_t(256 * 7 + 222, "smc PENTAX-DA L 18-55mm f/3.5-5.6"));
        choices.insert(p_t(256 * 7 + 223, "Samsung/Schneider D-XENON 18-55mm f/3.5-5.6 II"));
        choices.insert(p_t(256 * 7 + 224, "smc PENTAX-DA 15mm f/4 ED AL Limited"));
        choices.insert(p_t(256 * 7 + 225, "Samsung/Schneider D-XENON 18-250mm f/3.5-6.3"));
        choices.insert(p_t(256 * 7 + 226, "smc PENTAX-DA* 55mm f/1.4 SDM (SDM unused)"));
        choices.insert(p_t(256 * 7 + 227, "smc PENTAX-DA* 60-250mm f/4 [IF] SDM (SDM unused)"));
        choices.insert(p_t(256 * 7 + 228, "Samsung 16-45mm f/4 ED"));
        choices.insert(p_t(256 * 7 + 229, "smc PENTAX-DA 18-55mm f/3.5-5.6 AL II"));
        choices.insert(p_t(256 * 7 + 230, "Tamron AF 17-50mm f/2.8 XR Di-II LD (Model A16)"));
        choices.insert(p_t(256 * 7 + 231, "smc PENTAX-DA 18-250mm f/3.5-6.3 ED AL [IF]"));
        choices.insert(p_t(256 * 7 + 233, "smc PENTAX-DA 35mm f/2.8 Macro Limited"));
        choices.insert(p_t(256 * 7 + 234, "smc PENTAX-DA* 300mm f/4 ED [IF] SDM (SDM unused)"));
        choices.insert(p_t(256 * 7 + 235, "smc PENTAX-DA* 200mm f/2.8 ED [IF] SDM (SDM unused)"));
        choices.insert(p_t(256 * 7 + 236, "smc PENTAX-DA 55-300mm f/4-5.8 ED"));
        choices.insert(p_t(256 * 7 + 238, "Tamron AF 18-250mm f/3.5-6.3 Di II LD Aspherical [IF] Macro"));
        choices.insert(p_t(256 * 7 + 241, "smc PENTAX-DA* 50-135mm f/2.8 ED [IF] SDM (SDM unused)"));
        choices.insert(p_t(256 * 7 + 242, "smc PENTAX-DA* 16-50mm f/2.8 ED AL [IF] SDM (SDM unused)"));
        choices.insert(p_t(256 * 7 + 243, "smc PENTAX-DA 70mm f/2.4 Limited"));
        choices.insert(p_t(256 * 7 + 244, "smc PENTAX-DA 21mm f/3.2 AL Limited"));
        choices.insert(p_t(256 * 8 + 0, "Sigma 50-150mm f/2.8 II APO EX DC HSM"));
        choices.insert(p_t(256 * 8 + 3, "Sigma AF 18-125mm f/3.5-5.6 DC"));
        choices.insert(p_t(256 * 8 + 4, "Sigma 50mm f/1.4 EX DG HSM"));
        choices.insert(p_t(256 * 8 + 7, "Sigma 24-70mm f/2.8 IF EX DG HSM"));
        choices.insert(p_t(256 * 8 + 8, "Sigma 18-250mm f/3.5-6.3 DC OS HSM"));
        choices.insert(p_t(256 * 8 + 11, "Sigma 10-20mm f/3.5 EX DC HSM"));
        choices.insert(p_t(256 * 8 + 12, "Sigma 70-300mm f/4-5.6 DG OS"));
        choices.insert(p_t(256 * 8 + 13, "Sigma 120-400mm f/4.5-5.6 APO DG OS HSM"));
        choices.insert(p_t(256 * 8 + 14, "Sigma 17-70mm f/2.8-4.0 DC Macro OS HSM"));
        choices.insert(p_t(256 * 8 + 15, "Sigma 150-500mm f/5-6.3 APO DG OS HSM"));
        choices.insert(p_t(256 * 8 + 16, "Sigma 70-200mm f/2.8 EX DG Macro HSM II"));
        choices.insert(p_t(256 * 8 + 17, "Sigma 50-500mm f/4.5-6.3 DG OS HSM"));
        choices.insert(p_t(256 * 8 + 18, "Sigma 8-16mm f/4.5-5.6 DC HSM"));
        choices.insert(p_t(256 * 8 + 21, "Sigma 17-50mm f/2.8 EX DC OS HSM"));
        choices.insert(p_t(256 * 8 + 22, "Sigma 85mm f/1.4 EX DG HSM"));
        choices.insert(p_t(256 * 8 + 23, "Sigma 70-200mm f/2.8 APO EX DG OS HSM"));
        choices.insert(p_t(256 * 8 + 25, "Sigma 17-50mm f/2.8 EX DC HSM"));
        choices.insert(p_t(256 * 8 + 27, "Sigma 18-200mm f/3.5-6.3 II DC HSM"));
        choices.insert(p_t(256 * 8 + 28, "Sigma 18-250mm f/3.5-6.3 DC Macro HSM"));
        choices.insert(p_t(256 * 8 + 29, "Sigma 35mm f/1.4 DG HSM"));
        choices.insert(p_t(256 * 8 + 30, "Sigma 17-70mm f/2.8-4 DC Macro HSM | C"));
        choices.insert(p_t(256 * 8 + 31, "Sigma 18-35mm f/1.8 DC HSM"));
        choices.insert(p_t(256 * 8 + 32, "Sigma 30mm f/1.4 DC HSM | A"));
        choices.insert(p_t(256 * 8 + 34, "Sigma 18-300mm f/3.5-6.3 DC Macro HSM"));
        choices.insert(p_t(256 * 8 + 59, "HD PENTAX-D FA 150-450mm f/4.5-5.6 ED DC AW"));
        choices.insert(p_t(256 * 8 + 60, "HD PENTAX-D FA* 70-200mm f/2.8 ED DC AW"));
        choices.insert(p_t(256 * 8 + 61, "HD PENTAX-D FA 28-105mm f/3.5-5.6 ED DC WR"));
        choices.insert(p_t(256 * 8 + 62, "HD PENTAX-D FA 24-70mm f/2.8 ED SDM WR"));
        choices.insert(p_t(256 * 8 + 63, "HD PENTAX-D FA 15-30mm f/2.8 ED SDM WR"));
        choices.insert(p_t(256 * 8 + 197, "HD PENTAX-DA 55-300mm f/4.5-6.3 ED PLM WR RE"));
        choices.insert(p_t(256 * 8 + 198, "smc PENTAX-DA L 18-50mm f/4-5.6 DC WR RE"));
        choices.insert(p_t(256 * 8 + 199, "HD PENTAX-DA 18-50mm f/4-5.6 DC WR RE"));
        choices.insert(p_t(256 * 8 + 200, "HD PENTAX-DA 16-85mm f/3.5-5.6 ED DC WR"));
        choices.insert(p_t(256 * 8 + 209, "HD PENTAX-DA 20-40mm f/2.8-4 ED Limited DC WR"));
        choices.insert(p_t(256 * 8 + 210, "smc PENTAX-DA 18-270mm f/3.5-6.3 ED SDM"));
        choices.insert(p_t(256 * 8 + 211, "HD PENTAX-DA 560mm f/5.6 ED AW"));
        choices.insert(p_t(256 * 8 + 215, "smc PENTAX-DA 18-135mm f/3.5-5.6 ED AL [IF] DC WR"));
        choices.insert(p_t(256 * 8 + 226, "smc PENTAX-DA* 55mm f/1.4 SDM"));
        choices.insert(p_t(256 * 8 + 227, "smc PENTAX-DA* 60-250mm f/4 [IF] SDM"));
        choices.insert(p_t(256 * 8 + 232, "smc PENTAX-DA 17-70mm f/4 AL [IF] SDM"));
        choices.insert(p_t(256 * 8 + 234, "smc PENTAX-DA* 300mm f/4 ED [IF] SDM"));
        choices.insert(p_t(256 * 8 + 235, "smc PENTAX-DA* 200mm f/2.8 ED [IF] SDM"));
        choices.insert(p_t(256 * 8 + 241, "smc PENTAX-DA* 50-135mm f/2.8 ED [IF] SDM"));
        choices.insert(p_t(256 * 8 + 242, "smc PENTAX-DA* 16-50mm f/2.8 ED AL [IF] SDM"));
        choices.insert(p_t(256 * 8 + 255, "Sigma Lens (8 255)"));
        choices.insert(p_t(256 * 8 + 255, "Sigma 70-200mm f/2.8 EX DG Macro HSM II"));
        choices.insert(p_t(256 * 8 + 255, "Sigma 150-500mm f/5-6.3 DG APO [OS] HSM"));
        choices.insert(p_t(256 * 8 + 255, "Sigma 50-150mm f/2.8 II APO EX DC HSM"));
        choices.insert(p_t(256 * 8 + 255, "Sigma 4.5mm f/2.8 EX DC HSM Circular Fisheye"));
        choices.insert(p_t(256 * 8 + 255, "Sigma 50-200mm f/4-5.6 DC OS"));
        choices.insert(p_t(256 * 8 + 255, "Sigma 24-70mm f/2.8 EX DG HSM"));
        choices.insert(p_t(256 * 9 + 0, "645 Manual Lens"));
        choices.insert(p_t(256 * 10 + 0, "645 A Series Lens"));
        choices.insert(p_t(256 * 11 + 1, "smc PENTAX-FA 645 75mm f/2.8"));
        choices.insert(p_t(256 * 11 + 2, "smc PENTAX-FA 645 45mm f/2.8"));
        choices.insert(p_t(256 * 11 + 3, "smc PENTAX-FA* 645 300mm f/4 ED [IF]"));
        choices.insert(p_t(256 * 11 + 4, "smc PENTAX-FA 645 45-85mm f/4.5"));
        choices.insert(p_t(256 * 11 + 5, "smc PENTAX-FA 645 400mm f/5.6 ED [IF]"));
        choices.insert(p_t(256 * 11 + 7, "smc PENTAX-FA 645 Macro 120mm f/4"));
        choices.insert(p_t(256 * 11 + 8, "smc PENTAX-FA 645 80-160mm f/4.5"));
        choices.insert(p_t(256 * 11 + 9, "smc PENTAX-FA 645 200mm f/4 [IF]"));
        choices.insert(p_t(256 * 11 + 10, "smc PENTAX-FA 645 150mm f/2.8 [IF]"));
        choices.insert(p_t(256 * 11 + 11, "smc PENTAX-FA 645 35mm f/3.5 AL [IF]"));
        choices.insert(p_t(256 * 11 + 12, "smc PENTAX-FA 645 300mm f/5.6 ED [IF]"));
        choices.insert(p_t(256 * 11 + 14, "smc PENTAX-FA 645 55-110mm f/5.6"));
        choices.insert(p_t(256 * 11 + 16, "smc PENTAX-FA 645 33-55mm f/4.5 AL"));
        choices.insert(p_t(256 * 11 + 17, "smc PENTAX-FA 645 150-300mm f/5.6 ED [IF]"));
        choices.insert(p_t(256 * 11 + 21, "HD PENTAX-D FA 645 35mm f/3.5 AL [IF]"));
        choices.insert(p_t(256 * 13 + 18, "smc PENTAX-D FA 645 55mm f/2.8 AL [IF] SDM AW"));
        choices.insert(p_t(256 * 13 + 19, "smc PENTAX-D FA 645 25mm f/4 AL [IF] SDM AW"));
        choices.insert(p_t(256 * 13 + 20, "HD PENTAX-D FA 645 90mm f/2.8 ED AW SR"));
        choices.insert(p_t(256 * 13 + 253, "HD PENTAX-DA 645 28-45mm f/4.5 ED AW SR"));
        choices.insert(p_t(256 * 21 + 0, "Pentax Q Manual Lens"));
        choices.insert(p_t(256 * 21 + 1, "01 Standard Prime 8.5mm f/1.9"));
        choices.insert(p_t(256 * 21 + 2, "02 Standard Zoom 5-15mm f/2.8-4.5"));
        choices.insert(p_t(256 * 21 + 6, "06 Telephoto Zoom 15-45mm f/2.8"));
        choices.insert(p_t(256 * 21 + 7, "07 Mount Shield 11.5mm f/9"));
        choices.insert(p_t(256 * 21 + 8, "08 Wide Zoom 3.8-5.9mm f/3.7-4"));
        choices.insert(p_t(256 * 22 + 3, "03 Fish-eye 3.2mm f/5.6"));
        choices.insert(p_t(256 * 22 + 4, "04 Toy Lens Wide 6.3mm f/7.1"));
        choices.insert(p_t(256 * 22 + 5, "05 Toy Lens Telephoto 18mm f/8"));
    }
    virtual std::string toString (Tag* t)
    {
        double *liArray = NULL;
        double maxApertureAtFocal = 0.;
        double focalLength = 0.;
        int lensID = 256 * t->toInt(0, BYTE) + t->toInt(1, BYTE);
        TagDirectory *root = t->getParent()->getRoot();

        if (root) {

            Tag *t1;
            t1 = root->findTag("FocalLength");  // Should get tag 0x920A (rational64u) from the standard Exif tag list

            if( t1) {
                focalLength = t1->toDouble();    // Focal Length
            }

            t1 = root->findTag("MaxAperture");

            if(t1) {
                double maxAperture = t1->toDouble(); // MaxApertureValue at focal Length

                if (maxAperture != 0.) {
                    maxApertureAtFocal = maxAperture;
                } else {
                    t1 = root->findTag("NominalMaxAperture");

                    if(t1) {
                        maxApertureAtFocal = t1->toDouble();
                    }
                }
            }

            t1 = root->getTagP("LensInfo");

            if(t1) {
                liArray = t1->toDoubleArray();
            }

            // Focal length below 10mm are set to 0 by the camera in the standard Exif tag, so we'll look into the makernotes
            // This value will have decimals, which reflects more precision... or imprecision, due to the packed form of this value, who knows?
            if (focalLength == 0.) {
                rtexif::TagDirectory* mnote = root->findTag("MakerNote")->getDirectory();
                rtexif::Tag* flt = mnote->getTagP("LensInfo/FocalLength");

                if (flt) {
                    focalLength = flt->toDouble ();
                } else if ((flt = mnote->getTagP ("FocalLength"))) {
                    focalLength = flt->toDouble();
                }
            }
        }

        std::string retval = guess( lensID, focalLength, maxApertureAtFocal, liArray);

        if(liArray) {
            delete [] liArray;
        }

        return retval;
    }
};
PALensTypeInterpreter paLensTypeInterpreter;

class PASRResultInterpreter: public Interpreter
{
public:
    PASRResultInterpreter() { }
    virtual std::string toString (Tag* t)
    {
        std::ostringstream str;
        int b = t->toInt(0, BYTE);

        if (!b) {
            str << "Not stabilized";
        } else if (b & 1) {
            str << "Stabilized";
        } else if (b & 64) {
            str << "Not Ready";
        }

        return str.str();
    }
};
PASRResultInterpreter paSRResultInterpreter;

class PAHighISONoiseInterpreter: public ChoiceInterpreter
{
public:
    // HighISONoiseReduction
    PAHighISONoiseInterpreter()
    {
        choices[0] = "Off";
        choices[1] = "Weakest";
        choices[2] = "Weak";
        choices[3] = "Strong";
        choices[4] = "Medium";
        choices[255] = "Auto";
    }
};
PAHighISONoiseInterpreter paHighISONoiseInterpreter;

class PAMonochromeFilterEffectInterpreter: public ChoiceInterpreter
{
public:
    PAMonochromeFilterEffectInterpreter()
    {
        choices[1] = "Green";
        choices[2] = "Yellow";
        choices[3] = "Orange";
        choices[4] = "Red";
        choices[5] = "Magenta";
        choices[6] = "Blue";
        choices[7] = "Cyan";
        choices[8] = "Infrared";
        choices[65535] = "None";
    }
};
PAMonochromeFilterEffectInterpreter paMonochromeFilterEffectInterpreter;

class PAMonochromeToningInterpreter: public ChoiceInterpreter
{
public:
    PAMonochromeToningInterpreter()
    {
        choices[0] = "-4";
        choices[1] = "-3";
        choices[2] = "-2";
        choices[3] = "-1";
        choices[4] = "0";
        choices[5] = "1";
        choices[6] = "2";
        choices[7] = "3";
        choices[8] = "4";
        choices[65535] = "None";
    }
};
PAMonochromeToningInterpreter paMonochromeToningInterpreter;

class PAShadowCorrectionInterpreter: public ChoiceInterpreter
{
public:
    PAShadowCorrectionInterpreter()
    {
        choices[        0 ] = "Off";
        choices[        1 ] = "On";
        choices[        2 ] = "Auto 2";
        choices[ 1 << 8 | 1 ] = "Weak";
        choices[ 1 << 8 | 2 ] = "Normal";
        choices[ 1 << 8 | 3 ] = "Strong";
        choices[ 2 << 8 | 4 ] = "Auto";
    }

    virtual std::string toString (Tag* t)
    {
        int idx = 0;

        if (t->getCount() == 1) {
            idx = t->toInt(0, BYTE);
        } else if (t->getCount() == 2) {
            idx = t->toInt(0, BYTE) << 8 | t->toInt(1, BYTE);
        }

        std::map<int, std::string>::iterator r  = choices.find (idx);
        std::ostringstream s;
        s << ((r != choices.end()) ? r->second : "n/a");
        return s.str();
    }
};
PAShadowCorrectionInterpreter paShadowCorrectionInterpreter;

class PAISOAutoParametersInterpreter: public ChoiceInterpreter
{
public:
    PAISOAutoParametersInterpreter()
    {
        choices[1] = "Slow";
        choices[2] = "Standard";
        choices[3] = "Fast";
    }
    virtual std::string toString (Tag* t)
    {
        std::map<int, std::string>::iterator r  = choices.find (t->toInt(0, BYTE));
        std::ostringstream s;
        s << ((r != choices.end()) ? r->second : "n/a");
        return s.str();
    }
};
PAISOAutoParametersInterpreter paISOAutoParametersInterpreter;

class PABleachBypassToningInterpreter: public ChoiceInterpreter
{
public:
    PABleachBypassToningInterpreter()
    {
        choices[1] = "Green";
        choices[2] = "Yellow";
        choices[3] = "Orange";
        choices[4] = "Red";
        choices[5] = "Magenta";
        choices[6] = "Purple";
        choices[7] = "Blue";
        choices[8] = "Cyan";
        choices[65535] = "Off";
    }
};
PABleachBypassToningInterpreter paBleachBypassToningInterpreter;

class PABlurControlInterpreter: public ChoiceInterpreter
{
public:
    PABlurControlInterpreter()
    {
        choices[0] = "Off";
        choices[1] = "Low";
        choices[2] = "Medium";
        choices[3] = "High";
    }
    virtual std::string toString (Tag* t)
    {
        std::map<int, std::string>::iterator r  = choices.find (t->toInt(0, BYTE));
        std::ostringstream s;
        s << ((r != choices.end()) ? r->second : "n/a");
        return s.str();
    }
};
PABlurControlInterpreter paBlurControlInterpreter;

class PAHDRInterpreter: public ChoiceInterpreter
{
    std::map<int, std::string> choices1;
    std::map<int, std::string> choices2;
public:
    PAHDRInterpreter()
    {
        choices[0]    = "Off";
        choices[1]    = "HDR Auto";
        choices[2]    = "HDR 1";
        choices[3]    = "HDR 2";
        choices[4]    = "HDR 3";

        choices1[0]   = "Auto-align Off";
        choices1[1]   = "Auto-align On";

        choices2[0]   = "n/a";
        choices2[4]   = "1 EV";
        choices2[8]   = "2 EV";
        choices2[12]  = "3 EV";
    }
    virtual std::string toString (Tag* t)
    {
        std::map<int, std::string>::iterator r  = choices.find  (t->toInt(0, BYTE));
        std::map<int, std::string>::iterator r1 = choices1.find (t->toInt(1, BYTE));
        std::map<int, std::string>::iterator r2 = choices2.find (t->toInt(2, BYTE));
        std::ostringstream s;
        s << ((r != choices.end() ) ?  r->second : "") << std::endl;
        s << ((r1 != choices1.end()) ? r1->second : "") << std::endl;
        s << ((r2 != choices2.end()) ? r2->second : "");
        return s.str();
    }
};
PAHDRInterpreter paHDRInterpreter;

class PACrossProcessInterpreter: public ChoiceInterpreter
{
public:
    PACrossProcessInterpreter()
    {
        choices[ 0] = "Off";
        choices[ 1] = "Randow";
        choices[ 2] = "Preset 1";
        choices[ 3] = "Preset 2";
        choices[ 4] = "Preset 3";
        choices[33] = "Favorite 1";
        choices[34] = "Favorite 2";
        choices[35] = "Favorite 3";
    }
};
PACrossProcessInterpreter paCrossProcessInterpreter;

class PAPowerSourceInterpreter: public ChoiceInterpreter
{
public:
    PAPowerSourceInterpreter()
    {
        choices[2] = "Body Battery";
        choices[3] = "Grip Battery ";
        choices[4] = "External Power Supply";
    }
};
PAPowerSourceInterpreter paPowerSourceInterpreter;

class PALensModelQInterpreter: public Interpreter
{
public:
    PALensModelQInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        char buffer[31];
        buffer[0] = 0;  //
        return buffer;  // TODO: how to get the string content!?

        // normal path below (copy the content of the string), but has to be bug fixed
        memcpy(buffer, t->getValue(), 30);
        buffer[30] = 0;
        return buffer;
    }
};
PALensModelQInterpreter paLensModelQInterpreter;

class PALensInfoQInterpreter: public Interpreter
{
public:
    PALensInfoQInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        char buffer[21];
        buffer[0] = 0;
        return buffer;  // TODO: how to get the string content!?

        // normal path below (copy the content of the string), but has to be bug fixed
        memcpy(buffer, t->getValue(), 20);
        buffer[20] = 0;
        return buffer;
    }
};
PALensInfoQInterpreter paLensInfoQInterpreter;

class PAFlashExposureCompInterpreter: public Interpreter
{
public:
    PAFlashExposureCompInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        int a;

        if (t->getCount() == 1) {
            a = t->toInt(0, SLONG) / 256;    // int32u
        } else {
            a = t->toInt(0, SBYTE) / 6;    // int8u[2]
        }

        char buffer[10];
        sprintf (buffer, "%d", a );
        return buffer;
    }
    virtual double toDouble (Tag* t, int ofs)
    {
        int a;

        if (t->getCount() == 1) {
            a = t->toInt(0, SLONG) / 256;    // int32u
        } else {
            a = t->toInt(0, SBYTE) / 6;    // int8u[2]
        }

        return double(a);
    }
};
PAFlashExposureCompInterpreter paFlashExposureCompInterpreter;

class PAFocalLengthInterpreter: public Interpreter
{
public:
    PAFocalLengthInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        double a = double(t->toInt(0, LONG));

        if(a > 1.) {
            char buffer[10];
            sprintf (buffer, "%.2f", a / 100. );
            return buffer;
        } else {
            return "n/a";
        }
    }
    virtual double toDouble (Tag* t, int ofs)
    {
        double a = double(t->toInt(0, LONG));

        if(a > 1.) {
            return a / 100.;
        } else {
            return 0.;
        }
    }
};
PAFocalLengthInterpreter paFocalLengthInterpreter;

class PALensDataFocalLengthInterpreter: public Interpreter
{
public:
    PALensDataFocalLengthInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        int a = t->toInt(0, BYTE);
        float b = float(10 * int(a >> 2)) * pow(4.f, float(int(a & 0x03) - 2));

        if(b > 1.f) {
            char buffer[10];
            sprintf (buffer, "%.2f", b );
            return buffer;
        } else {
            return "n/a";
        }
    }
    virtual double toDouble (Tag* t, int ofs)
    {
        int a = t->toInt(ofs, BYTE);
        float b = float(10 * int(a >> 2)) * pow(4.f, float(int(a & 0x03) - 2));

        if(b > 1.f) {
            return b;
        } else {
            return 0.;
        }
    }
};
PALensDataFocalLengthInterpreter paLensDataFocalLengthInterpreter;

class PAISOfInterpreter: public Interpreter
{
public:
    PAISOfInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        int a = t->toInt(0, BYTE);
        char buffer[32];
        double v = 100.*exp(double(a - 32) * log(2.) / 8.);
        sprintf (buffer, "%.1f", v );
        return buffer;
    }
    virtual double toDouble (Tag* t, int ofs)
    {
        int a = t->toInt(0, BYTE);
        return 100.*exp(double(a - 32) * log(2.) / 8.);
    }
};
PAISOfInterpreter paISOfInterpreter;

class PAMaxApertureInterpreter: public Interpreter
{
public:
    PAMaxApertureInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        int a = t->toInt(0, BYTE);
        a &= 0x7F;

        if(a > 1) {
            char buffer[32];
            double v = pow(2.0, (a - 1) / 32.0);

            if( v < 0. || v > 1000. ) {
                return "undef";
            }

            sprintf (buffer, "%.1f", v );
            return buffer;
        } else {
            return "n/a";
        }
    }
    virtual double toDouble (Tag* t, int ofs)
    {
        int a = t->toInt(0, BYTE);
        a &= 0x7F;

        if(a > 1) {
            return pow(2.0, double(a - 1) / 32.0);
        } else {
            return 0.;
        }
    }
};
PAMaxApertureInterpreter paMaxApertureInterpreter;

class PAAEXvInterpreter: public Interpreter
{
public:
    PAAEXvInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        int a = t->toInt(0, BYTE);
        char buffer[32];
        double v = double(a - 64) / 8.;
        sprintf (buffer, "%.1f", v );
        return buffer;
    }
    virtual double toDouble (Tag* t, int ofs)
    {
        int a = t->toInt(0, BYTE);
        return double(a - 64) / 8.;
    }
};
PAAEXvInterpreter paAEXvInterpreter;

class PAAEBXvInterpreter: public Interpreter
{
public:
    PAAEBXvInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        int a = t->toInt(0, SBYTE);
        char buffer[32];
        double v = double(a) / 8.;
        sprintf (buffer, "%.1f", v );
        return buffer;
    }
    virtual double toDouble (Tag* t, int ofs)
    {
        int a = t->toInt(0, SBYTE);
        return double(a) / 8.;
    }
};
PAAEBXvInterpreter paAEBXvInterpreter;

class PAApertureInterpreter: public Interpreter
{
public:
    PAApertureInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        int a = t->toInt(0, BYTE);
        char buffer[32];
        double v = exp((double(a) - 68.) * log(2.) / 16.);
        sprintf (buffer, "%.1f", v );
        return buffer;
    }
    virtual double toDouble (Tag* t, int ofs)
    {
        int a = t->toInt(0, BYTE);
        return exp((double(a) - 68.) * log(2.) / 16.);
    }
};
PAApertureInterpreter paApertureInterpreter;

class PAExposureTimeInterpreter: public Interpreter
{
public:
    PAExposureTimeInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        int a = t->toInt(0, BYTE);
        char buffer[32];
        double v = 24.*exp(-(double(a) - 32.) * log(2.) / 8.);
        sprintf (buffer, "%.6f", v );
        return buffer;
    }
    virtual double toDouble (Tag* t, int ofs)
    {
        int a = t->toInt(0, BYTE);
        return 24.*exp(-(double(a) - 32.) * log(2.) / 8.);
    }
};
PAExposureTimeInterpreter paExposureTimeInterpreter;

class PANominalMinApertureInterpreter: public Interpreter
{
public:
    PANominalMinApertureInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        char buffer[32];
        int a = t->toInt(0, BYTE);
        int mina = a & 0x0F;
        sprintf (buffer, "%.1f", double(int(pow(2.0, double(mina + 10) / 4.0) + 0.2)));
        return buffer;
    }
    virtual double toDouble (Tag* t, int ofs)
    {
        int a = t->toInt(0, BYTE) & 0x0F;
        return double(int(pow(2.0, double(a + 10) / 4.0) + 0.2));
    }
};
PANominalMinApertureInterpreter paNominalMinApertureInterpreter;

class PANominalMaxApertureInterpreter: public Interpreter
{
public:
    PANominalMaxApertureInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        char buffer[32];
        int a = t->toInt(0, BYTE);
        int maxa = (a & 0xF0) >> 4;
        sprintf (buffer, "%.1f", double(int(pow(2.0, double(maxa) / 4.0) + 0.2)) );
        return buffer;
    }
    virtual double toDouble (Tag* t, int ofs)
    {
        int a = ( t->toInt(0, BYTE) & 0xF0) >> 4;
        return double(int(pow(2.0, double(a) / 4.0) + 0.2));
    }
};
PANominalMaxApertureInterpreter paNominalMaxApertureInterpreter;

class PAFlashStatusInterpreter: public ChoiceInterpreter
{
public:
    PAFlashStatusInterpreter()
    {
        choices[0] = "Off";
        choices[1] = "Off (1)";
        choices[2] = "External, Did not fire";
        choices[6] = "External, Fired";
        choices[8] = "Internal, Did not fire (0x08)";
        choices[9] = "Internal, Did not fire";
        choices[13] = "Internal, Fired";
    }
};
PAFlashStatusInterpreter paFlashStatusInterpreter;

class PAInternalFlashModeInterpreter: public ChoiceInterpreter
{
public:
    PAInternalFlashModeInterpreter()
    {
        choices[0] = "n/a - Off-Auto-Aperture";
        choices[134] = "Fired, Wireless (Control)";
        choices[149] = "Fired, Wireless (Master)";
        choices[192] = "Fired";
        choices[193] = "Fired, Red-eye reduction";
        choices[194] = "Fired, Auto";
        choices[195] = "Fired, Auto, Red-eye reduction";
        choices[198] = "Fired, Wireless (Control), Fired normally not as control";
        choices[200] = "Fired, Slow-sync";
        choices[201] = "Fired, Slow-sync, Red-eye reduction";
        choices[202] = "Fired, Trailing-curtain Sync";
        choices[240] = "Did not fire, Normal";
        choices[241] = "Did not fire, Red-eye reduction";
        choices[242] = "Did not fire, Auto";
        choices[243] = "Did not fire, Auto, Red-eye reduction";
        choices[244] = "Did not fire, (Unknown 0xf4)";
        choices[245] = "Did not fire, Wireless (Master)";
        choices[246] = "Did not fire, Wireless (Control)";
        choices[248] = "Did not fire, Slow-sync";
        choices[249] = "Did not fire, Slow-sync, Red-eye reduction";
        choices[250] = "Did not fire, Trailing-curtain Sync";
    }
};
PAInternalFlashModeInterpreter paInternalFlashModeInterpreter;

class PAExternalFlashModeInterpreter: public ChoiceInterpreter
{
public:
    PAExternalFlashModeInterpreter()
    {
        choices[0] = "n/a - Off-Auto-Aperture";
        choices[63] = "Off";
        choices[64] = "On, Auto";
        choices[191] = "On, Flash Problem";
        choices[192] = "On, Manual";
        choices[196] = "On, P-TTL Auto";
        choices[197] = "On, Contrast-control Sync";
        choices[198] = "On, High-speed Sync";
        choices[204] = "On, Wireless";
        choices[205] = "On, Wireless, High-speed Sync";
        choices[240] = "Not Connected";
    }
};
PAExternalFlashModeInterpreter paExternalFlashModeInterpreter;

class PAExternalFlashExposureCompInterpreter: public ChoiceInterpreter
{
public:
    PAExternalFlashExposureCompInterpreter()
    {
        choices[0] = "n/a";
        choices[144] = "n/a (Manual Mode)";
        choices[164] = "-3.0";
        choices[167] = "-2.5";
        choices[168] = "-2.0";
        choices[171] = "-1.5";
        choices[172] = "-1.0";
        choices[175] = "-0.5";
        choices[176] = "0.0";
        choices[179] = "0.5";
        choices[180] = "1.0";
    }
};
PAExternalFlashExposureCompInterpreter paExternalFlashExposureCompInterpreter;

class PAExternalFlashBounceInterpreter: public ChoiceInterpreter
{
public:
    PAExternalFlashBounceInterpreter()
    {
        choices[0] = "n/a";
        choices[16] = "Direct";
        choices[48] = "Bonce";
    }
};
PAExternalFlashBounceInterpreter paExternalFlashBounceInterpreter;

class PAExternalFlashGNInterpreter: public Interpreter
{
public:
    PAExternalFlashGNInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        char buffer[1024];
        int b = t->toInt(0, BYTE) & 0x1F;
        sprintf (buffer, "%.0f", pow(2., b / 16. + 4) );
        return buffer;
    }
};
PAExternalFlashGNInterpreter paExternalFlashGNInterpreter;

class PAEVStepsInterpreter: public Interpreter
{
public:
    PAEVStepsInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        std::ostringstream str;

        if( t->toInt(0, BYTE) & 0x20 ) {
            str << "1/3 EV steps";
        } else {
            str << "1/2 EV steps";
        }

        return str.str();
    }
};
PAEVStepsInterpreter paEVStepsInterpreter;

class PAEDialinInterpreter: public Interpreter
{
public:
    PAEDialinInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        std::ostringstream str;

        if(  t->toInt(0, BYTE) & 0x40 ) {
            str << "P Shift";
        } else {
            str << "Tv or Av";
        }

        return str.str();
    }
};
PAEDialinInterpreter paEDialinInterpreter;

class PAApertureRingUseInterpreter: public Interpreter
{
public:
    PAApertureRingUseInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        std::ostringstream str;

        if(  t->toInt(0, BYTE) & 0x80 ) {
            str << "Permitted";
        } else {
            str << "Prohibited";
        }

        return str.str();
    }
};
PAApertureRingUseInterpreter paApertureRingUseInterpreter;

class PAFlashOptionInterpreter: public ChoiceInterpreter
{
public:
    PAFlashOptionInterpreter()
    {
        choices[0] = "Normal";
        choices[1] = "Red-eye reduction";
        choices[2] = "Auto";
        choices[3] = "Auto, Red-eye reduction";
        choices[5] = "Wireless (Master)";
        choices[6] = "Wireless (Control)";
        choices[8] = "Slow-sync";
        choices[9] = "Slow-sync, Red-eye reduction";
        choices[10] = "Trailing-curtain Sync";
    }
    virtual std::string toString (Tag* t)
    {
        std::map<int, std::string>::iterator r = choices.find (t->toInt(0, BYTE) >> 4);

        if (r != choices.end()) {
            return r->second;
        } else {
            char buffer[1024];
            t->toString (buffer);
            return std::string (buffer);
        }
    }
};
PAFlashOptionInterpreter paFlashOptionInterpreter;

class PAMeteringMode2Interpreter: public Interpreter
{
public:
    PAMeteringMode2Interpreter() {}
    virtual std::string toString (Tag* t)
    {
        std::ostringstream str;
        int v = (t->toInt(0, BYTE) & 0xF);

        if(!v) {
            str << "Multi-segment";
        } else if(v & 1) {
            str << "Center-weighted average";
        } else if(v & 2) {
            str << "Spot";
        }

        return str.str();
    }
};
PAMeteringMode2Interpreter paMeteringMode2Interpreter;

class PAExposureBracketStepSizeInterpreter: public ChoiceInterpreter
{
public:
    PAExposureBracketStepSizeInterpreter()
    {
        choices[3] = "0.3";
        choices[4] = "0.5";
        choices[5] = "0.7";
        choices[8] = "1.0";
        choices[11] = "1.3";
        choices[12] = "1.5";
        choices[13] = "1.7";
        choices[16] = "2.0";
    }
};
PAExposureBracketStepSizeInterpreter paExposureBracketStepSizeInterpreter;

class PAPictureMode2Interpreter: public ChoiceInterpreter
{
public:
    PAPictureMode2Interpreter()
    {
        choices[0] = "Scene Mode";
        choices[1] = "Auto PICT";
        choices[2] = "Program AE";
        choices[3] = "Green Mode";
        choices[4] = "Shutter Speed Priority";
        choices[5] = "Aperture Priority";
        choices[6] = "Program Tv Shift";
        choices[7] = "Program Av Shift";
        choices[8] = "Manual";
        choices[9] = "Bulb";
        choices[10] = "Aperture Priority, Off-Auto-Aperture";
        choices[11] = "Manual, Off-Auto-Aperture";
        choices[12] = "Bulb, Off-Auto-Aperture";
        choices[13] = "Shutter & Aperture Priority AE";
        choices[15] = "Sensitivity Priority AE";
        choices[16] = "Flash X-Sync Speed AE";
    }
};
PAPictureMode2Interpreter paPictureMode2Interpreter;

class PAProgramLineInterpreter: public Interpreter
{
public:
    PAProgramLineInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        std::ostringstream str;
        int c = t->toInt(0, BYTE);

        switch(c & 0xf) {
        case 0:
            str << "Manual";
            break;

        case 1:
            str << "AF-S";
            break;

        case 2:
            str << "AF-C";
            break;

        case 3:
            str << "AF-A";
            break;
        }

        if( (c & 0xF0) == 0) {
            str << ", Point Selection Auto";
        } else if( c & 0x20 ) {
            str << ", Fixed Center Point Selected";
        } else if( c & 0x10 ) {
            str << ", Point Selected";
        }

        return str.str();
    }
};
PAProgramLineInterpreter paProgramLineInterpreter;

class PAAFModeInterpreter: public Interpreter
{
public:
    PAAFModeInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        switch(t->toInt(0, BYTE) & 0x3) {
        case 0:
            return "Normal";

        case 1:
            return "Hi Speed";

        case 2:
            return "Depth";

        case 3:
            return "MTF";
        }

        return"Normal";
    }

};
PAAFModeInterpreter paAFModeInterpreter;

class PAAFPointSelectedInterpreter: public Interpreter
{
public:
    PAAFPointSelectedInterpreter() {}
    virtual std::string toString (Tag* t)
    {
        const char *ps[] = {"Upper-left", "Top", "Upper-right", "Left", "Mid-left", "Center", "Mid-right", "Right", "Lower-left", "Bottom", "Lower-right"};
        int c = t->toInt(0, SHORT);

        if( !c ) {
            return "Auto";
        } else {
            for( int iBit = 0; iBit < 11; iBit++)
                if( c & (1 << iBit) ) {
                    return ps[iBit];
                }

            return "n/a";
        }
    }
};
PAAFPointSelectedInterpreter paAFPointSelectedInterpreter;

class PADriveMode2Interpreter: public Interpreter
{
public:
    PADriveMode2Interpreter() {}
    virtual std::string toString (Tag* t)
    {
        int c = t->toInt(0, BYTE);

        if( !c ) {
            return "Single-frame";
        } else if( c & 0x01) {
            return "Continuous";
        } else if( c & 0x02) {
            return "Continuous (Lo)";
        } else if( c & 0x04) {
            return "Self-timer (12 s)";
        } else if( c & 0x08) {
            return "Self-timer (2 s)";
        } else if( c & 0x10 ) {
            return "Remote Control (3 s delay)";
        } else if( c & 0x20) {
            return "Remote Control";
        } else if( c & 0x40) {
            return "Exposure Bracket";
        } else if( c & 0x80) {
            return "Multiple Exposure";
        } else {
            return "Unknown";
        }
    }
};
PADriveMode2Interpreter paDriveMode2Interpreter;

const TagAttrib pentaxAttribs[] = {
    {0, AC_WRITE,  0, 0, 0x0000, AUTO, "PentaxVersion", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0001, AUTO, "PentaxModelType", &stdInterpreter},
    {0, AC_SYSTEM, 0, 0, 0x0002, AUTO, "PreviewImageSize", &stdInterpreter},
    {0, AC_SYSTEM, 0, 0, 0x0003, AUTO, "PreviewImageLength", &stdInterpreter},
    {0, AC_SYSTEM, 0, 0, 0x0004, AUTO, "PreviewImageStart", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0005, AUTO, "PentaxModelID", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0006, AUTO, "Date", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0007, AUTO, "Time", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0008, AUTO, "Quality", &paQualityInterpreter},
    {0, AC_WRITE,  0, 0, 0x0009, AUTO, "PentaxImageSize", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x000b, AUTO, "PictureMode", &paPictureModeInterpreter},
    {0, AC_WRITE,  0, 0, 0x000c, AUTO, "FlashMode", &paFlashModeInterpreter},
    {0, AC_WRITE,  0, 0, 0x000d, AUTO, "FocusMode", &paFocusModeInterpreter},
    {0, AC_WRITE,  0, 0, 0x000e, AUTO, "AFPointSelected", &paAFPointInterpreter},
    {0, AC_WRITE,  0, 0, 0x000f, AUTO, "AFPointsInFocus", &paAFFocusInterpreter},
    {0, AC_WRITE,  0, 0, 0x0010, AUTO, "FocusPosition", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0012, AUTO, "ExposureTime", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0013, AUTO, "FNumber", &paFNumberInterpreter},
    {0, AC_WRITE,  0, 0, 0x0014, AUTO, "ISO", &paISOInterpreter},
    {0, AC_WRITE,  0, 0, 0x0015, AUTO, "LightReading", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0016, AUTO, "ExposureCompensation", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0017, AUTO, "MeteringMode", &paMeteringModeInterpreter},
    {0, AC_WRITE,  0, 0, 0x0018, AUTO, "AutoBracketing", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0019, AUTO, "WhiteBalance", &paWhiteBalanceInterpreter},
    {0, AC_WRITE,  0, 0, 0x001a, AUTO, "WhiteBalanceMode", &paWhiteBalanceModeInterpreter},
    {0, AC_WRITE,  0, 0, 0x001b, AUTO, "BlueBalance", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x001c, AUTO, "RedBalance", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x001d, AUTO, "FocalLength", &paFocalLengthInterpreter},
    {0, AC_WRITE,  0, 0, 0x001e, AUTO, "DigitalZoom", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x001f, AUTO, "Saturation", &paSaturationInterpreter},
    {0, AC_WRITE,  0, 0, 0x0020, AUTO, "Contrast", &paContrastInterpreter},
    {0, AC_WRITE,  0, 0, 0x0021, AUTO, "Sharpness", &paSharpnessInterpreter},
    {0, AC_WRITE,  0, 0, 0x0022, AUTO, "WorldTimeLocation", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0023, AUTO, "HometownCity", &stdInterpreter},
    {0, AC_NEW,    0, 0, 0x0024, AUTO, "DestinationCity", &stdInterpreter},
    {0, AC_NEW,    0, 0, 0x0025, AUTO, "HometownDST", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0026, AUTO, "DestinationDST", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0027, AUTO, "DSPFirmwareVersion", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0028, AUTO, "CPUFirmwareVersion", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0029, AUTO, "FrameNumber", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x002d, AUTO, "EffectiveLV", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0032, AUTO, "ImageProcessing", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0033, AUTO, "PictureMode", &paPictureModeInterpreter2},
    {0, AC_WRITE,  0, 0, 0x0034, AUTO, "DriveMode", &paDriveModeInterpreter},
    {0, AC_WRITE,  0, 0, 0x0037, AUTO, "ColorSpace", &paColorSpaceInterpreter},
    {0, AC_WRITE,  0, 0, 0x0038, AUTO, "ImageAreaOffset", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0039, AUTO, "RawImageSize", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x003c, AUTO, "AFPointsInFocus", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x003e, AUTO, "PreviewImageBorders", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x003f, AUTO, "LensType", &paLensTypeInterpreter},
    {0, AC_WRITE,  0, 0, 0x0040, AUTO, "SensitivityAdjust", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0041, AUTO, "ImageProcessingCount", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0047, AUTO, "CameraTemperature", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0048, AUTO, "AELock", &paOnOffInterpreter},
    {0, AC_WRITE,  0, 0, 0x0049, AUTO, "NoiseReduction", &paOnOffInterpreter},
    {0, AC_WRITE,  0, 0, 0x004d, AUTO, "FlashExposureComp", &paFlashExposureCompInterpreter},
    {0, AC_WRITE,  0, 0, 0x004f, AUTO, "ImageTone", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0050, AUTO, "ColorTemperature", &stdInterpreter},
    {0, AC_WRITE,  0, pentaxSRInfoAttribs, 0x005c, AUTO, "ShakeReductionInfo", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x005d, AUTO, "ShutterCount", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0069, AUTO, "DynamicRangeExpansion", &paOnOffInterpreter},
    {0, AC_WRITE,  0, 0, 0x0071, AUTO, "HighISONoiseReduction", &paHighISONoiseInterpreter},
    {0, AC_WRITE,  0, 0, 0x0072, AUTO, "AFAdjustment", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0073, AUTO, "MonochromeFilterEffect", &paMonochromeFilterEffectInterpreter},
    {0, AC_WRITE,  0, 0, 0x0074, AUTO, "MonochromeToning", &paMonochromeToningInterpreter},
    {0, AC_WRITE,  0, 0, 0x0076, AUTO, "FaceDetect", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0077, AUTO, "FaceDetectFrameSize", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0079, AUTO, "ShadowCorrection", &paShadowCorrectionInterpreter},
    {0, AC_WRITE,  0, 0, 0x007a, AUTO, "ISOAutoParameters", &paISOAutoParametersInterpreter},
    {0, AC_WRITE,  0, 0, 0x007b, AUTO, "CrossProcess", &paCrossProcessInterpreter},
    {0, AC_WRITE,  0, pentaxLensCorrAttribs, 0x007d, AUTO, "LensCorr", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x007f, AUTO, "BleachBypassToning", &paBleachBypassToningInterpreter},
    {0, AC_WRITE,  0, 0, 0x0082, AUTO, "BlurControl", &paBlurControlInterpreter},
    {0, AC_WRITE,  0, 0, 0x0085, AUTO, "HDR", &paHDRInterpreter},
    {0, AC_WRITE,  0, 0, 0x0088, AUTO, "NeutralDensityFilter", &paOnOffInterpreter},
    {0, AC_WRITE,  0, 0, 0x008b, AUTO, "ISO", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0200, AUTO, "BlackPoint", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0201, AUTO, "WhitePoint", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0203, AUTO, "ColorMatrixA", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0204, AUTO, "ColorMatrixB", &stdInterpreter},
    {0, AC_WRITE,  0, pentaxCameraSettingsAttribs, 0x0205, AUTO, "CameraSettings", &stdInterpreter},
    {0, AC_WRITE,  0, pentaxAEInfoAttribs, 0x0206, AUTO, "AEInfo", &stdInterpreter},
    {0, AC_WRITE,  0, pentaxLensDataAttribs, 0x0207, AUTO, "LensInfo", &stdInterpreter},
    {0, AC_WRITE,  0, pentaxFlashInfoAttribs, 0x0208, AUTO, "FlashInfo", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0209, AUTO, "AEMeteringSegments", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x020a, AUTO, "FlashADump", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x020b, AUTO, "FlashBDump", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x020d, AUTO, "WB_RGGBLevelsDaylight", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x020e, AUTO, "WB_RGGBLevelsShade", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x020f, AUTO, "WB_RGGBLevelsCloudy", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0210, AUTO, "WB_RGGBLevelsTungsten", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0211, AUTO, "WB_RGGBLevelsFluorescentD", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0212, AUTO, "WB_RGGBLevelsFluorescentN", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0213, AUTO, "WB_RGGBLevelsFluorescentW", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0214, AUTO, "WB_RGGBLevelsFlash", &stdInterpreter},
    {0, AC_WRITE,  0, pentaxCameraInfoAttribs, 0x0215, AUTO, "CameraInfo", &stdInterpreter},
    {0, AC_WRITE,  0, pentaxBatteryInfoAttribs, 0x0216, AUTO, "BatteryInfo", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x021f, AUTO, "AFInfo", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0222, AUTO, "ColorInfo", &stdInterpreter},
    {0, AC_WRITE,  0, pentaxLensInfoQAttribs, 0x0239, AUTO, "LensInfoQ", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x03fe, AUTO, "DataDump", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x03ff, AUTO, "UnknownInfo", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0402, AUTO, "ToneCurve", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0403, AUTO, "ToneCurves", &stdInterpreter},
    {0, AC_WRITE,  0, 0, 0x0e00, AUTO, "PrintIM", &stdInterpreter},
    { -1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}
};

const TagAttrib pentaxSRInfoAttribs[] = {
    {0, AC_WRITE, 0, 0,  0, AUTO, "SRResult", &paSRResultInterpreter},
    {0, AC_WRITE, 0, 0,  1, AUTO, "ShakeReduction", &paShakeReductionInterpreter},
    {0, AC_WRITE, 0, 0,  2, AUTO, "SRHalfPressTime", &stdInterpreter},
    {0, AC_WRITE, 0, 0,  3, AUTO, "SRFocalLength", &stdInterpreter},
    { -1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}
};

const TagAttrib pentaxSRInfo2Attribs[] = {
    {0, AC_WRITE, 0, 0,  0, AUTO, "SRResult", &paSRResultInterpreter},    // assuming it's the same interpreter, but that's not sure
    {0, AC_WRITE, 0, 0,  1, AUTO, "ShakeReduction", &paShakeReduction2Interpreter},
    { -1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}
};

const TagAttrib pentaxLensDataAttribs[] = {
    {0, AC_WRITE, 0, 0, 9,  AUTO, "FocalLength", &paLensDataFocalLengthInterpreter},
    {0, AC_WRITE, 0, 0, 10, AUTO, "NominalMaxAperture", &paNominalMaxApertureInterpreter},
    {0, AC_WRITE, 0, 0, 10, AUTO, "NominalMinAperture", &paNominalMinApertureInterpreter},
    {0, AC_WRITE, 0, 0, 14, AUTO, "MaxAperture", &paMaxApertureInterpreter},
    { -1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}
};

const TagAttrib pentaxLensInfoQAttribs[] = {
    {0, AC_WRITE, 0, 0, 12, AUTO, "LensModel", &paLensModelQInterpreter},
    {0, AC_WRITE, 0, 0, 42, AUTO, "LensInfo", &paLensInfoQInterpreter},
    { -1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}
};

const TagAttrib pentaxLensCorrAttribs[] = {
    {0, AC_WRITE, 0, 0, 0, AUTO, "DistortionCorrection", &paOnOffInterpreter},
    {0, AC_WRITE, 0, 0, 1, AUTO, "ChromaticAberrationCorrection", &paOnOffInterpreter},
    {0, AC_WRITE, 0, 0, 2, AUTO, "VignettingCorrection", &paOnOffInterpreter},
    { -1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}
};

const TagAttrib pentaxCameraSettingsAttribs[] = {
    {0, AC_WRITE, 0, 0,  0, AUTO, "PictureMode2", &paPictureMode2Interpreter},
    {0, AC_WRITE, 0, 0,  1, AUTO, "ProgramLine", &paProgramLineInterpreter},
    {0, AC_WRITE, 0, 0,  1, AUTO, "EVSteps", &paEVStepsInterpreter},
    {0, AC_WRITE, 0, 0,  1, AUTO, "E-DialinProgram", &paEDialinInterpreter},
    {0, AC_WRITE, 0, 0,  1, AUTO, "ApertureRing", &paApertureRingUseInterpreter},
    {0, AC_WRITE, 0, 0,  2, AUTO, "FlashOptions", &paFlashOptionInterpreter},
    {0, AC_WRITE, 0, 0,  2, AUTO, "MeteringMode2", &paMeteringMode2Interpreter},
    {0, AC_WRITE, 0, 0,  3, AUTO, "AFMode", &paAFModeInterpreter},
    {0, AC_WRITE, 0, 0,  4, AUTO, "AFPointSelected2", &paAFPointSelectedInterpreter},
    {0, AC_WRITE, 0, 0,  7, AUTO, "DriveMode2", &paDriveMode2Interpreter},
    {0, AC_WRITE, 0, 0,  8, AUTO, "ExposureBracketStepSize", &paExposureBracketStepSizeInterpreter},
    {0, AC_WRITE, 0, 0,  9, AUTO, "BracketShotNumber", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 10, AUTO, "WhiteBalanceSet", &stdInterpreter},
    { -1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}
};

const TagAttrib pentaxAEInfoAttribs[] = {
    {0, AC_WRITE, 0, 0,  0, AUTO, "AEExposureTime", &paExposureTimeInterpreter},
    {0, AC_WRITE, 0, 0,  1, AUTO, "AEAperture", &paApertureInterpreter},
    {0, AC_WRITE, 0, 0,  2, AUTO, "AE_ISO", &paISOfInterpreter},
    {0, AC_WRITE, 0, 0,  3, AUTO, "AEXv", &paAEXvInterpreter},
    {0, AC_WRITE, 0, 0,  4, SBYTE, "AEBXv", &paAEBXvInterpreter},
    {0, AC_WRITE, 0, 0,  5, AUTO, "AEMinExposureTime", &paExposureTimeInterpreter},
    {0, AC_WRITE, 0, 0,  6, AUTO, "AEProgramMode", &paAEProgramModeInterpreter},
    {0, AC_WRITE, 0, 0,  9, AUTO, "AEMaxAperture", &paApertureInterpreter},
    {0, AC_WRITE, 0, 0, 10, AUTO, "AEMaxAperture2", &paApertureInterpreter},
    {0, AC_WRITE, 0, 0, 11, AUTO, "AEMinAperture", &paApertureInterpreter},
    {0, AC_WRITE, 0, 0, 12, AUTO, "AEMeteringMode", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 14, SBYTE, "FlashExposureCompSet", &stdInterpreter},
    { -1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}
};

const TagAttrib pentaxAEInfo2Attribs[] = {
    {0, AC_WRITE, 0, 0,  2, AUTO, "AEExposureTime", &paExposureTimeInterpreter},
    {0, AC_WRITE, 0, 0,  3, AUTO, "AEAperture", &paApertureInterpreter},
    {0, AC_WRITE, 0, 0,  4, AUTO, "AE_ISO", &paISOfInterpreter},
    {0, AC_WRITE, 0, 0,  5, AUTO, "AEXv", &paAEXvInterpreter},
    {0, AC_WRITE, 0, 0,  6, SBYTE, "AEBXv", &paAEBXvInterpreter},
    {0, AC_WRITE, 0, 0,  8, SBYTE, "AEError", &stdInterpreter},
//{0, AC_WRITE, 0, 0, 11, AUTO, "AEApertureSteps", &},
    {0, AC_WRITE, 0, 0, 15, AUTO, "SceneMode", &paSceneModeInterpreter},
    {0, AC_WRITE, 0, 0, 16, AUTO, "AEMaxAperture", &paApertureInterpreter},
    {0, AC_WRITE, 0, 0, 17, AUTO, "AEMaxAperture2", &paApertureInterpreter},
    {0, AC_WRITE, 0, 0, 18, AUTO, "AEMinAperture", &paApertureInterpreter},
    {0, AC_WRITE, 0, 0, 19, AUTO, "AEMinExposureTime", &paExposureTimeInterpreter},
    { -1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}
};

const TagAttrib pentaxAEInfo3Attribs[] = {
    {0, AC_WRITE, 0, 0, 16, AUTO, "AEExposureTime", &paExposureTimeInterpreter},
    {0, AC_WRITE, 0, 0, 17, AUTO, "AEAperture", &paApertureInterpreter},
    {0, AC_WRITE, 0, 0, 18, AUTO, "AE_ISO", &paISOfInterpreter},
    {0, AC_WRITE, 0, 0, 28, AUTO, "AEMaxAperture", &paApertureInterpreter},
    {0, AC_WRITE, 0, 0, 29, AUTO, "AEMaxAperture2", &paApertureInterpreter},
    {0, AC_WRITE, 0, 0, 30, AUTO, "AEMinAperture", &paApertureInterpreter},
    {0, AC_WRITE, 0, 0, 31, AUTO, "AEMinExposureTime", &paExposureTimeInterpreter},
    { -1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}
};

const TagAttrib pentaxFlashInfoAttribs[] = {
    {0, AC_WRITE, 0, 0,  0, AUTO, "FlashStatus", &paFlashStatusInterpreter},
    {0, AC_WRITE, 0, 0,  1, AUTO, "InternalFlashMode", &paInternalFlashModeInterpreter},
    {0, AC_WRITE, 0, 0,  2, AUTO, "ExternalFlashMode", &paExternalFlashModeInterpreter},
    {0, AC_WRITE, 0, 0,  3, AUTO, "InternalFlashStrength", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 24, AUTO, "ExternalFlashGuideNumber", &paExternalFlashGNInterpreter},
    {0, AC_WRITE, 0, 0, 25, AUTO, "ExternalFlashExposureComp", &paExternalFlashExposureCompInterpreter},
    {0, AC_WRITE, 0, 0, 26, AUTO, "ExternalFlashBounce", &paExternalFlashBounceInterpreter},
    { -1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}
};

const TagAttrib pentaxBatteryInfoAttribs[] = {
    {0, AC_WRITE, 0, 0,  0, AUTO, "PowerSource", &paPowerSourceInterpreter},
    {0, AC_WRITE, 0, 0,  1, AUTO, "BatteryStates", &stdInterpreter},
    {0, AC_WRITE, 0, 0,  2, AUTO, "BatteryADBodyNoLoad", &stdInterpreter},
    {0, AC_WRITE, 0, 0,  3, AUTO, "BatteryADBodyLoad", &stdInterpreter},
    {0, AC_WRITE, 0, 0,  4, AUTO, "BatteryADGripNoLoad", &stdInterpreter},
    {0, AC_WRITE, 0, 0,  5, AUTO, "BatteryADGripLoad", &stdInterpreter},
    { -1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}
};

const TagAttrib pentaxCameraInfoAttribs[] = {
    {0, AC_WRITE, 0, 0,  0, AUTO, "PentaxModelID", &stdInterpreter},
    {0, AC_WRITE, 0, 0,  1, AUTO, "ManufactureDate", &stdInterpreter},
    {0, AC_WRITE, 0, 0,  2, AUTO, "ProductionCode", &stdInterpreter},
    {0, AC_WRITE, 0, 0,  4, AUTO, "InternalSerialNumber", &stdInterpreter},
    { -1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}
};

}
#endif




