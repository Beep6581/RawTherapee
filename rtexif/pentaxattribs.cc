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

#include <rtexif.h>
#include <string>
#include <map>
#include <math.h>
#include <sstream>
#include <iomanip>

namespace rtexif {


class PAQualityInterpreter : public ChoiceInterpreter {
    public:
        PAQualityInterpreter () {
            choices[0]      = "Good";
            choices[1]      = "Better";
            choices[2]      = "Best";
            choices[3]      = "TIFF";
            choices[4]      = "RAW";
            choices[5]      = "Premium";
        }
};
PAQualityInterpreter paQualityInterpreter;

class PAOnOffInterpreter : public ChoiceInterpreter {
    public:
        PAOnOffInterpreter () {
            choices[0]      = "Off";
            choices[1]      = "On";
        }
};
PAOnOffInterpreter paOnOffInterpreter;

class PAPictureModeInterpreter : public ChoiceInterpreter {
    public:
        PAPictureModeInterpreter () {
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
            choices[65] = "Half-length Portrait";
        }
};
PAPictureModeInterpreter paPictureModeInterpreter;

class PAFlashModeInterpreter : public ChoiceInterpreter {
    public:
        PAFlashModeInterpreter () {
            choices[0x0] = "Auto, Did not fire";
            choices[0x1] = "Off";
            choices[0x2] = "On, Did not fire";
            choices[0x3] = "Auto, Did not fire, Red-eye reduction";
            choices[0x100] = "Auto, Fired";
            choices[0x102] = "On";
            choices[0x103] = "Auto, Fired, Red-eye reduction";
            choices[0x104] = "On, Red-eye reduction";
            choices[0x105] = "On, Wireless (Master)";
            choices[0x106] = "On, Wireless (Control)";
            choices[0x108] = "On, Soft";
            choices[0x109] = "On, Slow-sync";
            choices[0x10a] = "On, Slow-sync, Red-eye reduction";
            choices[0x10b] = "On, Trailing-curtain Sync";
        }
};
PAFlashModeInterpreter paFlashModeInterpreter;

class PAFocusModeInterpreter : public ChoiceInterpreter {
    public:
        PAFocusModeInterpreter () {
            choices[0] = "Normal";
            choices[1] = "Macro";
            choices[2] = "Infinity";
            choices[3] = "Manual";
            choices[4] = "Super Macro";
            choices[5] = "Pan Focus";
            choices[16] = "AF-S";
            choices[17] = "AF-C";
            choices[18] = "AF-A";
        }
};
PAFocusModeInterpreter paFocusModeInterpreter;

class PAAFPointInterpreter : public ChoiceInterpreter {
    public:
        PAAFPointInterpreter        () {
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
            choices[65532] = "Face Recognition AF";
            choices[65533] = "Automatic Tracking AF";
            choices[65534] = "Fixed Center";
            choices[65535] = "Auto";
        }
};
PAAFPointInterpreter paAFPointInterpreter;

class PAAFFocusInterpreter : public ChoiceInterpreter {
    public:
        PAAFFocusInterpreter        () {
            choices[0x0] = "Fixed Center or Multiple";
            choices[0x1] = "Top-left";
            choices[0x2] = "Top-center";
            choices[0x3] = "Top-right";
            choices[0x4] = "Left";
            choices[0x5] = "Center";
            choices[0x6] = "Right";
            choices[0x7] = "Bottom-left";
            choices[0x8] = "Bottom-center";
            choices[0x9] = "Bottom-right";
            choices[0xffff] = "None";
        }
};
PAAFFocusInterpreter paAFFocusInterpreter;

class PAISOInterpreter : public ChoiceInterpreter {
    public:
        PAISOInterpreter        () {
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
            choices[400] = "400";
            choices[800] = "800";
            choices[1600] = "1600";
            choices[3200] = "320";
        }
};
PAISOInterpreter paISOInterpreter;

class PAMeteringModeInterpreter : public ChoiceInterpreter {
    public:
        PAMeteringModeInterpreter () {
            choices[0] = "Multi-segment";
            choices[1] = "Center-weighted average";
            choices[2] = "Spot";
        }
};
PAMeteringModeInterpreter paMeteringModeInterpreter;

class PAWhiteBalanceInterpreter : public ChoiceInterpreter {
    public:
        PAWhiteBalanceInterpreter () {
            choices[0] = "Auto";
            choices[1] = "Daylight";
            choices[2] = "Shade";
            choices[3] = "Fluorescent";
            choices[4] = "Tungsten";
            choices[5] = "Manual";
            choices[6] = "DaylightFluorescent";
            choices[7] = "DaywhiteFluorescent";
            choices[8] = "WhiteFluorescent";
            choices[9] = "Flash";
            choices[10] = "Cloudy";
            choices[17] = "Kelvin";
            choices[65534] = "Unknown";
            choices[65535] = "User Selected";
        }
};
PAWhiteBalanceInterpreter paWhiteBalanceInterpreter;

class PAWhiteBalanceModeInterpreter : public ChoiceInterpreter {
    public:
        PAWhiteBalanceModeInterpreter () {
            choices[1] = "Auto (Daylight)";
            choices[2] = "Auto (Shade)";
            choices[3] = "Auto (Flash)";
            choices[4] = "Auto (Tungsten)";
            choices[6] = "Auto (DaylightFluorescent)";
            choices[7] = "Auto (DaywhiteFluorescent)";
            choices[8] = "Auto (WhiteFluorescent)";
            choices[10] = "Auto (Cloudy)";
            choices[65534] = "Preset (Fireworks?)";
            choices[65535] = "User-Selected";        
        }
};
PAWhiteBalanceModeInterpreter paWhiteBalanceModeInterpreter;

class PASaturationInterpreter : public ChoiceInterpreter {
    public:
        PASaturationInterpreter () {
            choices[0] = "Low";
            choices[1] = "Normal";
            choices[2] = "High";
            choices[3] = "Med Low";
            choices[4] = "Med High";
            choices[5] = "Very Low";
            choices[6] = "Very High";
        }
};
PASaturationInterpreter paSaturationInterpreter;

class PAContrastInterpreter : public ChoiceInterpreter {
    public:
        PAContrastInterpreter () {
            choices[0] = "Low";
            choices[1] = "Normal";
            choices[2] = "High";
            choices[3] = "Med Low";
            choices[4] = "Med High";
            choices[5] = "Very Low";
            choices[6] = "Very High";
        }
};
PAContrastInterpreter paContrastInterpreter;

class PASharpnessInterpreter : public ChoiceInterpreter {
    public:
        PASharpnessInterpreter () {
            choices[0] = "Soft";
            choices[1] = "Normal";
            choices[2] = "Hard";
            choices[3] = "Med Soft";
            choices[4] = "Med Hard";
            choices[5] = "Very Soft";
            choices[6] = "Very Hard";
        }
};
PASharpnessInterpreter paSharpnessInterpreter;

class PALensTypeInterpreter : public ChoiceInterpreter {
    public:
        PALensTypeInterpreter () {
            choices[256*0+ 0] = "M-42 or No Lens";
            choices[256*1+ 0] = "K,M Lens";
            choices[256*2+ 0] = "A Series Lens";
            choices[256*3+ 0] = "SIGMA";
            choices[256*3+ 17] = "smc PENTAX-FA SOFT 85mm F2.8";
            choices[256*3+ 18] = "smc PENTAX-F 1.7X AF ADAPTER";
            choices[256*3+ 19] = "smc PENTAX-F 24-50mm F4";
            choices[256*3+ 20] = "smc PENTAX-F 35-80mm F4-5.6";
            choices[256*3+ 21] = "smc PENTAX-F 80-200mm F4.7-5.6";
            choices[256*3+ 22] = "smc PENTAX-F FISH-EYE 17-28mm F3.5-4.5";
            choices[256*3+ 23] = "smc PENTAX-F 100-300mm F4.5-5.6";
            choices[256*3+ 24] = "smc PENTAX-F 35-135mm F3.5-4.5";
            choices[256*3+ 25] = "smc PENTAX-F 35-105mm F4-5.6 or SIGMA or Tokina";
            choices[256*3+ 26] = "smc PENTAX-F* 250-600mm F5.6 ED[IF]";
            choices[256*3+ 27] = "smc PENTAX-F 28-80mm F3.5-4.5";
            choices[256*3+ 28] = "smc PENTAX-F 35-70mm F3.5-4.5";
            choices[256*3+ 29] = "PENTAX-F 28-80mm F3.5-4.5 or SIGMA AF 18-125mm F3.5-5.6 DC";
            choices[256*3+ 30] = "PENTAX-F 70-200mm F4-5.6";
            choices[256*3+ 31] = "smc PENTAX-F 70-210mm F4-5.6";
            choices[256*3+ 32] = "smc PENTAX-F 50mm F1.4";
            choices[256*3+ 33] = "smc PENTAX-F 50mm F1.7";
            choices[256*3+ 34] = "smc PENTAX-F 135mm F2.8 [IF]";
            choices[256*3+ 35] = "smc PENTAX-F 28mm F2.8";
            choices[256*3+ 36] = "SIGMA 20mm F1.8 EX DG ASPHERICAL RF";
            choices[256*3+ 38] = "smc PENTAX-F* 300mm F4.5 ED[IF]";
            choices[256*3+ 39] = "smc PENTAX-F* 600mm F4 ED[IF]";
            choices[256*3+ 40] = "smc PENTAX-F MACRO 100mm F2.8";
            choices[256*3+ 41] = "smc PENTAX-F MACRO 50mm F2.8 or Sigma 50mm F2,8 MACRO";
            choices[256*3+ 44] = "Tamron 35-90mm F4 AF or various SIGMA models";
            choices[256*3+ 46] = "SIGMA APO 70-200mm F2.8 EX";
            choices[256*3+ 50] = "smc PENTAX-FA 28-70 F4 AL";
            choices[256*3+ 51] = "SIGMA 28mm F1.8 EX DG ASPHERICAL MACRO";
            choices[256*3+ 52] = "smc PENTAX-FA 28-200mm F3.8-5.6 AL[IF]";
            choices[256*3+ 53] = "smc PENTAX-FA 28-80mm F3.5-5.6 AL";
            choices[256*3+ 247] = "smc PENTAX-DA FISH-EYE 10-17mm F3.5-4.5 ED[IF]";
            choices[256*3+ 248] = "smc PENTAX-DA 12-24mm F4 ED AL[IF]";
            choices[256*3+ 250] = "smc PENTAX-DA 50-200mm F4-5.6 ED";
            choices[256*3+ 251] = "smc PENTAX-DA 40mm F2.8 Limited";
            choices[256*3+ 252] = "smc PENTAX-DA 18-55mm F3.5-5.6 AL";
            choices[256*3+ 253] = "smc PENTAX-DA 14mm F2.8 ED[IF]";
            choices[256*3+ 254] = "smc PENTAX-DA 16-45mm F4 ED AL";
            choices[256*3+ 255] = "SIGMA";
            choices[256*4+ 1] = "smc PENTAX-FA SOFT 28mm F2.8";
            choices[256*4+ 2] = "smc PENTAX-FA 80-320mm F4.5-5.6";
            choices[256*4+ 3] = "smc PENTAX-FA 43mm F1.9 Limited";
            choices[256*4+ 6] = "smc PENTAX-FA 35-80mm F4-5.6";
            choices[256*4+ 12] = "smc PENTAX-FA 50mm F1.4";
            choices[256*4+ 15] = "smc PENTAX-FA 28-105mm F4-5.6 [IF]";
            choices[256*4+ 16] = "TAMRON AF 80-210mm F4-5.6 (178D)";
            choices[256*4+ 19] = "TAMRON SP AF 90mm F2.8 (172E)";
            choices[256*4+ 20] = "smc PENTAX-FA 28-80mm F3.5-5.6";
            choices[256*4+ 21] = "Cosina AF 100-300mm F5.6-6.7";
            choices[256*4+ 22] = "TOKINA 28-80mm F3.5-5.6";
            choices[256*4+ 23] = "smc PENTAX-FA 20-35mm F4 AL";
            choices[256*4+ 24] = "smc PENTAX-FA 77mm F1.8 Limited";
            choices[256*4+ 25] = "TAMRON SP AF 14mm F2.8";
            choices[256*4+ 26] = "smc PENTAX-FA MACRO 100mm F3.5";
            choices[256*4+ 27] = "TAMRON AF28-300mm F/3.5-6.3 LD Aspherical[IF] MACRO (285D)";
            choices[256*4+ 28] = "smc PENTAX-FA 35mm F2 AL";
            choices[256*4+ 29] = "TAMRON AF 28-200mm F/3.8-5.6 LD Super II MACRO (371D)";
            choices[256*4+ 34] = "smc PENTAX-FA 24-90mm F3.5-4.5 AL[IF]";
            choices[256*4+ 35] = "smc PENTAX-FA 100-300mm F4.7-5.8";
            choices[256*4+ 36] = "TAMRON AF70-300mm F/4-5.6 LD MACRO";
            choices[256*4+ 37] = "TAMRON SP AF 24-135mm F3.5-5.6 AD AL (190D)";
            choices[256*4+ 38] = "smc PENTAX-FA 28-105mm F3.2-4.5 AL[IF]";
            choices[256*4+ 39] = "smc PENTAX-FA 31mm F1.8AL Limited";
            choices[256*4+ 41] = "TAMRON AF 28-200mm Super Zoom F3.8-5.6 Aspherical XR [IF] MACRO (A03)";
            choices[256*4+ 43] = "smc PENTAX-FA 28-90mm F3.5-5.6";
            choices[256*4+ 44] = "smc PENTAX-FA J 75-300mm F4.5-5.8 AL";
            choices[256*4+ 45] = "TAMRON 28-300mm F3.5-6.3 Ultra zoom XR";
            choices[256*4+ 46] = "smc PENTAX-FA J 28-80mm F3.5-5.6 AL";
            choices[256*4+ 47] = "smc PENTAX-FA J 18-35mm F4-5.6 AL";
            choices[256*4+ 49] = "TAMRON SP AF 28-75mm F2.8 XR Di (A09)";
            choices[256*4+ 51] = "smc PENTAX-D FA 50mm F2.8 MACRO";
            choices[256*4+ 52] = "smc PENTAX-D FA 100mm F2.8 MACRO";
            choices[256*4+ 75] = "TAMRON SP AF 70-200 F2.8 Di LD [IF] Macro (A001)";
            choices[256*4+ 229] = "smc PENTAX-DA 18-55mm F3.5-5.6 AL II";
            choices[256*4+ 230] = "TAMRON SP AF 17-50mm F2.8 XR Di II";
            choices[256*4+ 231] = "smc PENTAX-DA 18-250mm F3.5-6.3 ED AL [IF]";
            choices[256*4+ 237] = "Samsung/Schneider D-XENOGON 10-17mm F3.5-4.5";
            choices[256*4+ 239] = "Samsung D-XENON 12-24mm F4 ED AL [IF]";
            choices[256*4+ 243] = "smc PENTAX-DA 70mm F2.4 Limited";
            choices[256*4+ 244] = "smc PENTAX-DA 21mm F3.2 AL Limited";
            choices[256*4+ 245] = "Schneider D-XENON 50-200mm";
            choices[256*4+ 246] = "Schneider D-XENON 18-55mm";
            choices[256*4+ 247] = "smc PENTAX-DA 10-17mm F3.5-4.5 ED [IF] Fisheye zoom";
            choices[256*4+ 248] = "smc PENTAX-DA 12-24mm F4 ED AL [IF]";
            choices[256*4+ 249] = "TAMRON XR DiII 18-200mm F3.5-6.3 (A14)";
            choices[256*4+ 250] = "smc PENTAX-DA 50-200mm F4-5.6 ED";
            choices[256*4+ 251] = "smc PENTAX-DA 40mm F2.8 Limited";
            choices[256*4+ 252] = "smc PENTAX-DA 18-55mm F3.5-5.6 AL";
            choices[256*4+ 253] = "smc PENTAX-DA 14mm F2.8 ED[IF]";
            choices[256*4+ 254] = "smc PENTAX-DA 16-45mm F4 ED AL";
            choices[256*5+ 1] = "smc PENTAX-FA* 24mm F2 AL[IF]";
            choices[256*5+ 2] = "smc PENTAX-FA 28mm F2.8 AL";
            choices[256*5+ 3] = "smc PENTAX-FA 50mm F1.7";
            choices[256*5+ 4] = "smc PENTAX-FA 50mm F1.4";
            choices[256*5+ 5] = "smc PENTAX-FA* 600mm F4 ED[IF]";
            choices[256*5+ 6] = "smc PENTAX-FA* 300mm F4.5 ED[IF]";
            choices[256*5+ 7] = "smc PENTAX-FA 135mm F2.8 [IF]";
            choices[256*5+ 8] = "smc PENTAX-FA MACRO 50mm F2.8";
            choices[256*5+ 9] = "smc PENTAX-FA MACRO 100mm F2.8";
            choices[256*5+ 10] = "smc PENTAX-FA* 85mm F1.4 [IF]";
            choices[256*5+ 11] = "smc PENTAX-FA* 200mm F2.8 ED[IF]";
            choices[256*5+ 12] = "smc PENTAX-FA 28-80mm F3.5-4.7";
            choices[256*5+ 13] = "smc PENTAX-FA 70-200mm F4-5.6";
            choices[256*5+ 14] = "smc PENTAX-FA* 250-600mm F5.6 ED[IF]";
            choices[256*5+ 15] = "smc PENTAX-FA 28-105mm F4-5.6";
            choices[256*5+ 16] = "smc PENTAX-FA 100-300mm F4.5-5.6";
            choices[256*5+ 98] = "smc PENTAX-FA 100-300mm F4.5-5.6";
            choices[256*6+ 1] = "smc PENTAX-FA* 85mm F1.4 [IF]";
            choices[256*6+ 2] = "smc PENTAX-FA* 200mm F2.8 ED[IF]";
            choices[256*6+ 3] = "smc PENTAX-FA* 300mm F2.8 ED[IF]";
            choices[256*6+ 4] = "smc PENTAX-FA* 28-70mm F2.8 AL";
            choices[256*6+ 5] = "smc PENTAX-FA* 80-200mm F2.8 ED[IF]";
            choices[256*6+ 6] = "smc PENTAX-FA* 28-70mm F2.8 AL";
            choices[256*6+ 7] = "smc PENTAX-FA* 80-200mm F2.8 ED[IF]";
            choices[256*6+ 8] = "smc PENTAX-FA 28-70mm F4AL";
            choices[256*6+ 9] = "smc PENTAX-FA 20mm F2.8";
            choices[256*6+ 10] = "smc PENTAX-FA* 400mm F5.6 ED[IF]";
            choices[256*6+ 13] = "smc PENTAX-FA* 400mm F5.6 ED[IF]";
            choices[256*6+ 14] = "smc PENTAX-FA* MACRO 200mm F4 ED[IF]";
            choices[256*7+ 0] = "smc PENTAX-DA 21mm F3.2 AL Limited";
            choices[256*7+ 75] = "TAMRON SP AF 70-200mm F2.8 Di LD [IF] Macro (A001)";
            choices[256*7+ 217] = "smc PENTAX-DA 50-200mm F4-5.6 ED WR";
            choices[256*7+ 218] = "smc PENTAX-DA 18-55mm F3.5-5.6 AL WR";
            choices[256*7+ 220] = "TAMRON SP AF 10-24mm F3.5-4.5 Di II LD Aspherical [IF]";
            choices[256*7+ 222] = "smc PENTAX-DA 18-55mm F3.5-5.6 AL II";
            choices[256*7+ 223] = "Samsung D-XENON 18-55mm F3.5-5.6 II";
            choices[256*7+ 224] = "smc PENTAX-DA 15mm F4 ED AL Limited";
            choices[256*7+ 225] = "Samsung D-XENON 18-250mm F3.5-6.3";
            choices[256*7+ 229] = "smc PENTAX-DA 18-55mm F3.5-5.6 AL II";
            choices[256*7+ 230] = "TAMRON AF 17-50mm F2.8 XR Di-II LD (Model A16)";
            choices[256*7+ 231] = "smc PENTAX-DA 18-250mm F3.5-6.3 ED AL [IF]";
            choices[256*7+ 233] = "smc PENTAX-DA 35mm F2.8 Macro Limited";
            choices[256*7+ 234] = "smc PENTAX-DA* 300mm F4 ED [IF] SDM (SDM unused)";
            choices[256*7+ 235] = "smc PENTAX-DA* 200mm F2.8 ED [IF] SDM (SDM unused)";
            choices[256*7+ 236] = "smc PENTAX-DA 55-300mm F4-5.8 ED";
            choices[256*7+ 238] = "TAMRON AF 18-250mm F3.5-6.3 Di II LD Aspherical [IF] MACRO";
            choices[256*7+ 241] = "smc PENTAX-DA* 50-135mm F2.8 ED [IF] SDM (SDM unused)";
            choices[256*7+ 242] = "smc PENTAX-DA* 16-50mm F2.8 ED AL [IF] SDM (SDM unused)";
            choices[256*7+ 243] = "smc PENTAX-DA 70mm F2.4 Limited";
            choices[256*7+ 244] = "smc PENTAX-DA 21mm F3.2 AL Limited";
            choices[256*8+ 226] = "smc PENTAX-DA* 55mm F1.4 SDM";
            choices[256*8+ 227] = "smc PENTAX DA* 60-250mm F4 [IF] SDM";
            choices[256*8+ 232] = "smc PENTAX-DA 17-70mm F4 AL [IF] SDM";
            choices[256*8+ 234] = "smc PENTAX-DA* 300mm F4 ED [IF] SDM";
            choices[256*8+ 235] = "smc PENTAX-DA* 200mm F2.8 ED [IF] SDM";
            choices[256*8+ 241] = "smc PENTAX-DA* 50-135mm F2.8 ED [IF] SDM";
            choices[256*8+ 242] = "smc PENTAX-DA* 16-50mm F2.8 ED AL [IF] SDM";
            choices[256*8+ 255] = "Sigma 70-200mm F2.8 EX DG Macro HSM II or 150-500mm F5-6.3 DG OS";
        }
        virtual std::string toString (Tag* t) {
            return choices[256*t->toInt(0,BYTE) + t->toInt(1,BYTE)];
        }
};
PALensTypeInterpreter paLensTypeInterpreter;

class PASRInfoInterpreter : public Interpreter {
    public:
        PASRInfoInterpreter () { }

        virtual std::string toString (Tag* t) {
            std::ostringstream str;
            int b = t->toInt(0,BYTE);
            if (!b)
                str << "SRResult = Not stabilized" << std::endl;
            else if (b & 1)
                str << "SRResult = Stabilized" << std::endl;
            b = t->toInt(1,BYTE);
            if (!b)
                str << "ShakeReduction = Off" << std::endl;
            else
                str << "ShakeReduction = On" << std::endl;
            str << "SRHalfPressTime = " << t->toInt(2,BYTE) << std::endl;
            str << "SRFocalLength = " << t->toInt(3,BYTE);
            return str.str();
        }
};
PASRInfoInterpreter paSRInfoInterpreter;


const TagAttrib pentaxAttribs[] = {
 0, 1, 0, 0, 0x0001, "PentaxVersion", &stdInterpreter,
 0, 1, 0, 0, 0x0001, "PentaxModelType", &stdInterpreter,
 0, 2, 0, 0, 0x0002, "PreviewImageSize", &stdInterpreter,
 0, 2, 0, 0, 0x0003, "PreviewImageLength", &stdInterpreter,
 0, 2, 0, 0, 0x0004, "PreviewImageStart", &stdInterpreter,
 0, 1, 0, 0, 0x0005, "PentaxModelID", &stdInterpreter,
 0, 1, 0, 0, 0x0006, "Date", &stdInterpreter,
 0, 1, 0, 0, 0x0007, "Time", &stdInterpreter,
 0, 1, 0, 0, 0x0008, "Quality", &paQualityInterpreter,
 0, 1, 0, 0, 0x0009, "PentaxImageSize", &stdInterpreter,
 0, 1, 0, 0, 0x000b, "PictureMode", &paPictureModeInterpreter,
 0, 1, 0, 0, 0x000c, "FlashMode", &paFlashModeInterpreter,
 0, 1, 0, 0, 0x000d, "FocusMode", &paFocusModeInterpreter,
 0, 1, 0, 0, 0x000e, "AFPointSelected", &paAFPointInterpreter,
 0, 1, 0, 0, 0x000f, "AFPointsInFocus", &paAFFocusInterpreter,
 0, 1, 0, 0, 0x0010, "FocusPosition", &stdInterpreter,
 0, 1, 0, 0, 0x0012, "ExposureTime", &stdInterpreter,
 0, 1, 0, 0, 0x0013, "FNumber", &stdInterpreter,
 0, 1, 0, 0, 0x0014, "ISO", &paISOInterpreter,
 0, 1, 0, 0, 0x0015, "LightReading", &stdInterpreter,
 0, 1, 0, 0, 0x0016, "ExposureCompensation", &stdInterpreter,
 0, 1, 0, 0, 0x0017, "MeteringMode", &paMeteringModeInterpreter,
 0, 1, 0, 0, 0x0018, "AutoBracketing", &stdInterpreter,
 0, 1, 0, 0, 0x0019, "WhiteBalance", &paWhiteBalanceInterpreter,
 0, 1, 0, 0, 0x001a, "WhiteBalanceMode", &paWhiteBalanceModeInterpreter,
 0, 1, 0, 0, 0x001b, "BlueBalance", &stdInterpreter,
 0, 1, 0, 0, 0x001c, "RedBalance", &stdInterpreter,
 0, 1, 0, 0, 0x001d, "FocalLength", &stdInterpreter,
 0, 1, 0, 0, 0x001e, "DigitalZoom", &stdInterpreter,
 0, 1, 0, 0, 0x001f, "Saturation", &paSaturationInterpreter,
 0, 1, 0, 0, 0x0020, "Contrast", &paContrastInterpreter,
 0, 1, 0, 0, 0x0021, "Sharpness", &paSharpnessInterpreter,
 0, 1, 0, 0, 0x0022, "WorldTimeLocation", &stdInterpreter,
 0, 1, 0, 0, 0x0023, "HometownCity", &stdInterpreter,
 0, 3, 0, 0, 0x0024, "DestinationCity", &stdInterpreter,
 0, 3, 0, 0, 0x0025, "HometownDST", &stdInterpreter,
 0, 1, 0, 0, 0x0026, "DestinationDST", &stdInterpreter,
 0, 1, 0, 0, 0x0027, "DSPFirmwareVersion", &stdInterpreter,
 0, 1, 0, 0, 0x0028, "CPUFirmwareVersion", &stdInterpreter,
 0, 1, 0, 0, 0x0029, "FrameNumber", &stdInterpreter,
 0, 1, 0, 0, 0x002d, "EffectiveLV", &stdInterpreter,
 0, 1, 0, 0, 0x0032, "ImageProcessing", &stdInterpreter,
 0, 1, 0, 0, 0x0033, "PictureMode", &stdInterpreter,
 0, 1, 0, 0, 0x0034, "DriveMode", &stdInterpreter,
 0, 1, 0, 0, 0x0037, "ColorSpace", &stdInterpreter,
 0, 1, 0, 0, 0x0038, "ImageAreaOffset", &stdInterpreter,
 0, 1, 0, 0, 0x0039, "RawImageSize", &stdInterpreter,
 0, 1, 0, 0, 0x003c, "AFPointsInFocus", &stdInterpreter,
 0, 1, 0, 0, 0x003e, "PreviewImageBorders", &stdInterpreter,
 0, 1, 0, 0, 0x003f, "LensType", &paLensTypeInterpreter,
 0, 1, 0, 0, 0x0040, "SensitivityAdjust", &stdInterpreter,
 0, 1, 0, 0, 0x0041, "ImageProcessingCount", &stdInterpreter,
 0, 1, 0, 0, 0x0047, "CameraTemperature", &stdInterpreter,
 0, 1, 0, 0, 0x0048, "AELock", &paOnOffInterpreter,
 0, 1, 0, 0, 0x0049, "NoiseReduction", &paOnOffInterpreter,
 0, 1, 0, 0, 0x004d, "FlashExposureComp", &stdInterpreter,
 0, 1, 0, 0, 0x004f, "ImageTone", &stdInterpreter,
 0, 1, 0, 0, 0x0050, "ColorTemperature", &stdInterpreter,
 0, 1, 0, 0, 0x005c, "ShakeReductionInfo", &paSRInfoInterpreter,
 0, 1, 0, 0, 0x005d, "ShutterCount", &stdInterpreter,
 0, 1, 0, 0, 0x0200, "BlackPoint", &stdInterpreter,
 0, 1, 0, 0, 0x0201, "WhitePoint", &stdInterpreter,
 0, 1, 0, 0, 0x0205, "ShotInfo", &stdInterpreter,
 0, 1, 0, 0, 0x0206, "AEInfo", &stdInterpreter,
 0, 1, 0, 0, 0x0207, "LensInfo", &stdInterpreter,
 0, 1, 0, 0, 0x0208, "FlashInfo", &stdInterpreter,
 0, 1, 0, 0, 0x0209, "AEMeteringSegments", &stdInterpreter,
 0, 1, 0, 0, 0x020a, "FlashADump", &stdInterpreter,
 0, 1, 0, 0, 0x020b, "FlashBDump", &stdInterpreter,
 0, 1, 0, 0, 0x020d, "WB_RGGBLevelsDaylight", &stdInterpreter,
 0, 1, 0, 0, 0x020e, "WB_RGGBLevelsShade", &stdInterpreter,
 0, 1, 0, 0, 0x020f, "WB_RGGBLevelsCloudy", &stdInterpreter,
 0, 1, 0, 0, 0x0210, "WB_RGGBLevelsTungsten", &stdInterpreter,
 0, 1, 0, 0, 0x0211, "WB_RGGBLevelsFluorescentD", &stdInterpreter,
 0, 1, 0, 0, 0x0212, "WB_RGGBLevelsFluorescentN", &stdInterpreter,
 0, 1, 0, 0, 0x0213, "WB_RGGBLevelsFluorescentW", &stdInterpreter,
 0, 1, 0, 0, 0x0214, "WB_RGGBLevelsFlash", &stdInterpreter,
 0, 1, 0, 0, 0x0215, "CameraInfo", &stdInterpreter,
 0, 1, 0, 0, 0x0216, "BatteryInfo", &stdInterpreter,
 0, 1, 0, 0, 0x021f, "AFInfo", &stdInterpreter,
 0, 1, 0, 0, 0x0222, "ColorInfo", &stdInterpreter,
 0, 1, 0, 0, 0x03fe, "DataDump", &stdInterpreter,
 0, 1, 0, 0, 0x03ff, "UnknownInfo", &stdInterpreter,
 0, 1, 0, 0, 0x0402, "ToneCurve", &stdInterpreter,
 0, 1, 0, 0, 0x0403, "ToneCurves", &stdInterpreter,
 0, 1, 0, 0, 0x0e00, "PrintIM", &stdInterpreter,
 -1, 0, 0,  0, 0, "", NULL};
};
#endif

