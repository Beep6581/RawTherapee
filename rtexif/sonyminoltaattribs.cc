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
#ifndef _SONYMINOLTAATTRIBS_
#define _SONYMINOLTAATTRIBS_

#include <rtexif.h>
#include <string>
#include <map>
#include <math.h>
#include <sstream>
#include <iomanip>

namespace rtexif {

class SAOnOffInterpreter : public ChoiceInterpreter {
    public:
        SAOnOffInterpreter () {
            choices[0]      = "Off";
            choices[1]      = "On";
            choices[5]      = "On";
        }
};
SAOnOffInterpreter saOnOffInterpreter;

class SASceneModeInterpreter : public ChoiceInterpreter {
    public:
        SASceneModeInterpreter () {
            choices[0]  = "Normal (P,A,S or M)";
            choices[1]  = "Portrait";
            choices[2]  = "Text";
            choices[3]  = "Night Scene";
            choices[4]  = "Sunset";
            choices[5]  = "Sports";
            choices[6]  = "Landscape";
            choices[8]  = "Macro";
            choices[8]  = "Super Macro";
            choices[16] = "Auto";
            choices[17] = "Night Portrait";
        }
};
SASceneModeInterpreter saSceneModeInterpreter;

class SAZoneMatchingInterpreter : public ChoiceInterpreter {
    public:
        SAZoneMatchingInterpreter () {
            choices[0] = "ISO Setting Used";
            choices[1] = "High Key";
            choices[2] = "Low Key";
        }
};
SAZoneMatchingInterpreter saZoneMatchingInterpreter;

class SADynamicRangeOptimizerInterpreter : public ChoiceInterpreter {
    public:
        SADynamicRangeOptimizerInterpreter () {
            choices[0] = "Off";
            choices[1] = "Standard";
            choices[2] = "Advanced";
            choices[3] = "Auto";
            choices[8] = "Advanced Lv1";
            choices[9] = "Advanced Lv2";
            choices[10] = "Advanced Lv3";
            choices[11] = "Advanced Lv4";
            choices[12] = "Advanced Lv5";
        }
};
SADynamicRangeOptimizerInterpreter saDynamicRangeOptimizerInterpreter;

class SAColorModeInterpreter : public ChoiceInterpreter {
    public:
        SAColorModeInterpreter () {
            choices[0]  = "Standard";
            choices[1]  = "Vivid";
            choices[2]  = "Portrait";
            choices[3]  = "Landscape";
            choices[4]  = "Sunset";
            choices[5]  = "Night Scene";
            choices[6]  = "B&W";
            choices[7]  = "Adobe RGB";
            choices[12] = "Neutral";
            choices[100]= "Neutral";
            choices[101]= "Clear";
            choices[102]= "Deep";
            choices[103]= "Light";
            choices[104]= "Night View";
            choices[105]= "Autumn Leaves";
        }
};
SAColorModeInterpreter saColorModeInterpreter;

class SAExposureModeInterpreter : public ChoiceInterpreter {
    public:
        SAExposureModeInterpreter () {
            choices[0]  = "Auto";
            choices[5]  = "Landscape";
            choices[6]  = "Program";
            choices[7]  = "Aperture Priority";
            choices[8]  = "Shutter Priority";
            choices[9]  = "Night Scene";
            choices[15] = "Manual";
            choices[34] = "Panorama";
            choices[35] = "Handheld Twilight";
            choices[36] = "Anti Motion Blur";
        }
};
SAExposureModeInterpreter saExposureModeInterpreter;

class SAQualityInterpreter : public ChoiceInterpreter {
    public:
        SAQualityInterpreter () {
            choices[0]  = "Normal";
            choices[1]  = "Fine";
        }
};
SAQualityInterpreter saQualityInterpreter;

class SAAntiBlurInterpreter : public ChoiceInterpreter {
    public:
        SAAntiBlurInterpreter () {
            choices[0]  = "Off";
            choices[1]  = "On (Continuous)";
            choices[2]  = "On (Shooting)";
            choices[65535]  = "n/a";
        }
};
SAAntiBlurInterpreter saAntiBlurInterpreter;

class SALensIDInterpreter : public ChoiceInterpreter {
    public:
        SALensIDInterpreter () {
            choices[0] = "Minolta AF 28-85mm F3.5-4.5";
            choices[1] = "Minolta AF 80-200mm F2.8 HS-APO G";
            choices[2] = "Minolta AF 28-70mm F2.8 G";
            choices[3] = "Minolta AF 28-80mm F4-5.6";
            choices[5] = "Minolta AF 35-70mm F3.5-4.5";
            choices[6] = "Minolta AF 24-85mm F3.5-4.5 [New]";
            choices[7] = "Minolta AF 100-300mm F4.5-5.6 APO [New]";
            choices[8] = "Minolta AF 70-210mm F4.5-5.6";
            choices[9] = "Minolta AF 50mm F3.5 Macro";
            choices[10] = "Minolta AF 28-105mm F3.5-4.5 [New]";
            choices[11] = "Minolta AF 300mm F4 HS-APO G";
            choices[12] = "Minolta AF 100mm F2.8 Soft Focus";
            choices[13] = "Minolta AF 75-300mm F4.5-5.6";
            choices[14] = "Minolta AF 100-400mm F4.5-6.7 APO";
            choices[15] = "Minolta AF 400mm F4.5 HS-APO G";
            choices[16] = "Minolta AF 17-35mm F3.5 G";
            choices[17] = "Minolta AF 20-35mm F3.5-4.5";
            choices[18] = "Minolta AF 28-80mm F3.5-5.6 II";
            choices[19] = "Minolta AF 35mm F1.4";
            choices[20] = "Minolta/Sony STF 135mm F2.8 [T4.5]";
            choices[22] = "Minolta AF 35-80mm F4-5.6";
            choices[23] = "Minolta AF 200mm F4 G APO Macro";
            choices[24] = "Minolta/Sony AF 24-105mm F3.5-4.5 (D)";
            choices[25] = "Minolta AF 100-300mm F4.5-5.6 (APO D)";
            choices[27] = "Minolta AF 85mm F1.4 G";
            choices[28] = "Minolta AF 100mm F2.8 Macro (D)";
            choices[29] = "Minolta AF 75-300mm F4.5-5.6 (D)";
            choices[30] = "Minolta AF 28-80mm F3.5-5.6 (D)";
            choices[31] = "Minolta/Sony AF 50mm F2.8 Macro (D) or AF 50mm F3.5 Macro";
            choices[32] = "Minolta AF 300mm F2.8 G";
            choices[33] = "Minolta/Sony AF 70-200mm F2.8 G (D) SSM";
            choices[35] = "Minolta AF 85mm F1.4 G (D) Limited";
            choices[36] = "Minolta AF 28-100mm F3.5-5.6 (D)";
            choices[38] = "Minolta AF 17-35mm F2.8-4 (D)";
            choices[39] = "Minolta AF 28-75mm F2.8 (D)";
            choices[40] = "Minolta/Sony AF DT 18-70mm F3.5-5.6 (D)";
            choices[41] = "Minolta/Sony AF DT 11-18mm F4.5-5.6 (D)";
            choices[42] = "Minolta AF DT 18-200mm F3.5-6.3 (D)";
            choices[43] = "Minolta AF 35mm F1.4 G";
            choices[44] = "Minolta AF 50mm F1.4";
            choices[45] = "Carl Zeiss Planar T* 85mm F1.4 ZA";
            choices[46] = "Carl Zeiss Vario-Sonnar T* DT 16-80mm F3.5-4.5 ZA";
            choices[47] = "Carl Zeiss Sonnar T* 135mm F1.8 ZA";
            choices[48] = "Carl Zeiss Vario-Sonnar T* 24-70mm F2.8 ZA SSM";
            choices[49] = "Sony AF DT 55-200mm F4-5.6";
            choices[50] = "Sony AF DT 18-250mm F3.5-6.3";
            choices[51] = "Sony AF DT 16-105mm F3.5-5.6 or 55-200mm f/4-5.5";
            choices[52] = "Sony AF 70-300mm F4.5-5.6 G SSM";
            choices[53] = "Sony AF 70-400mm F4.5-5.6 G SSM";
            choices[54] = "Carl Zeiss Vario-Sonnar T* 16-35mm F2.8 ZA SSM";
            choices[55] = "Sony DT 18-55mm F3.5-5.6 SAM";
            choices[56] = "Sony AF DT 55-200mm F4-5.6 SAM";
            choices[57] = "Sony AF DT 50mm F1.8 SAM";
            choices[58] = "Sony AF DT 30mm F2.8 SAM Macro";
            choices[59] = "Sony AF 28-75mm F2.8 SAM";
            choices[128] = "Tamron or Sigma Lens";
            choices[129] = "Tamron 200-400mm F5.6 or 70-300mm f/4-5.6 LD";
            choices[135] = "Vivitar 28-210mm F3.5-5.6";
            choices[136] = "Tokina EMZ M100 AF 100mm F3.5";
            choices[137] = "Cosina 70-210mm F2.8-4 AF";
            choices[138] = "Soligor 19-35mm F3.5-4.5";
            choices[142] = "Voigtlander 70-300mm F4.5-5.6";
            choices[146] = "Voigtlander Macro APO-Lanthar 125mm F2.5 SL";
            choices[255] = "Tamron Lens";
            choices[2550] = "Minolta AF 50mm F1.7";
            choices[2551] = "Minolta AF 35-70mm F4";
            choices[2552] = "Minolta AF 28-85mm F3.5-4.5 [New]";
            choices[2553] = "Minolta AF 28-135mm F4-4.5";
            choices[2554] = "Minolta AF 35-105mm F3.5-4.5";
            choices[2555] = "Minolta AF 70-210mm F4 Macro";
            choices[2556] = "Minolta AF 135mm F2.8";
            choices[2557] = "Minolta AF 28mm F2.8";
            choices[2558] = "Minolta AF 24-50mm F4";
            choices[2560] = "Minolta AF 100-200mm F4.5";
            choices[2561] = "Minolta AF 75-300mm F4.5-5.6";
            choices[2562] = "Minolta/Sony AF 50mm F1.4 [New]";
            choices[2563] = "Minolta AF 300mm F2.8 G";
            choices[2564] = "Minolta AF 50mm F2.8 Macro";
            choices[2565] = "Minolta AF 600mm F4";
            choices[2566] = "Minolta AF 24mm F2.8";
            choices[2572] = "Minolta/Sony AF 500mm F8 Reflex";
            choices[2578] = "Minolta AF 16mm F2.8 Fisheye or Sigma Lens";
            choices[2579] = "Minolta AF 20mm F2.8";
            choices[2581] = "Minolta/Sony AF 100mm F2.8 Macro or Sigma or Tamron";
            choices[2585] = "Minolta AF 35-105mm F3.5-4.5 New";
            choices[2588] = "Minolta AF 70-210mm F3.5-4.5";
            choices[2589] = "Minolta AF 80-200 F2.8 APO";
            choices[2591] = "Minolta AF 35mm F1.4";
            choices[2592] = "Minolta AF 85mm F1.4 G (D)";
            choices[2593] = "Minolta AF 200mm F2.8 G APO";
            choices[2594] = "Minolta AF 3x-1x F1.7-2.8 Macro";
            choices[2596] = "Minolta AF 28mm F2";
            choices[2597] = "Minolta AF 35mm F2";
            choices[2598] = "Minolta AF 100mm F2";
            choices[2604] = "Minolta AF 80-200mm F4.5-5.6";
            choices[2605] = "Minolta AF 35-80mm F4-5.6";
            choices[2606] = "Minolta AF 100-300mm F4.5-5.6 (D)";
            choices[2607] = "Minolta AF 35-80mm F4-5.6";
            choices[2608] = "Minolta AF 300mm F2.8 G";
            choices[2609] = "Minolta AF 600mm F4 HS-APO G";
            choices[2612] = "Minolta AF 200mm F2.8 G HS-APO";
            choices[2613] = "Minolta AF 50mm F1.7 New";
            choices[2615] = "Minolta AF 28-105mm F3.5-4.5 Power Zoom";
            choices[2616] = "Minolta AF 35-200mm F4.5-5.6 Power Zoom";
            choices[2618] = "Minolta AF 28-80mm F4-5.6 Power Zoom";
            choices[2619] = "Minolta AF 80-200mm F4.5-5.6 Power Zoom";
            choices[2620] = "Minolta AF 28-70mm F2.8 G";
            choices[2621] = "Minolta AF 100-300mm F4.5-5.6 Power Zoom";
            choices[2624] = "Minolta AF 35-80mm F4-5.6 Power Zoom";
            choices[2628] = "Minolta AF 80-200mm F2.8 G";
            choices[2629] = "Minolta AF 85mm F1.4 New";
            choices[2631] = "Minolta/Sony AF 100-300mm F4.5-5.6 APO";
            choices[2632] = "Minolta AF 24-50mm F4 New";
            choices[2638] = "Minolta AF 50mm F2.8 Macro New";
            choices[2639] = "Minolta AF 100mm F2.8 Macro";
            choices[2641] = "Minolta AF 20mm F2.8 New";
            choices[2642] = "Minolta AF 24mm F2.8 New";
            choices[2644] = "Minolta AF 100-400mm F4.5-6.7 APO";
            choices[2662] = "Minolta AF 50mm F1.4 New";
            choices[2667] = "Minolta AF 35mm F2 New";
            choices[2668] = "Minolta AF 28mm F2 New";
            choices[2672] = "Minolta AF 24-105mm F3.5-4.5 (D)";
            choices[4574] = "Minolta AF 200mm F2.8 G x2";
            choices[4575] = "1.4 x Teleconverter";
            choices[4585] = "Tamron - SP AF 300 F2.8 LD IF";
            choices[25501] = "Minolta AF 50mm F1.7";
            choices[25511] = "Minolta AF 35-70mm F4";
            choices[25521] = "Minolta AF 28-85mm F3.5-4.5 [New]";
            choices[25531] = "Minolta AF 28-135mm F4-4.5";
            choices[25541] = "Minolta AF 35-105mm F3.5-4.5";
            choices[25551] = "Minolta AF 70-210mm F4 Macro";
            choices[25561] = "Minolta AF 135mm F2.8";
            choices[25571] = "Minolta AF 28mm F2.8";
            choices[25581] = "Minolta AF 24-50mm F4";
            choices[25601] = "Minolta AF 100-200mm F4.5";
            choices[25611] = "Minolta AF 75-300mm F4.5-5.6";
            choices[25621] = "Minolta/Sony AF 50mm F1.4 [New]";
            choices[25631] = "Minolta AF 300mm F2.8 G";
            choices[25641] = "Minolta AF 50mm F2.8 Macro";
            choices[25651] = "Minolta AF 600mm F4";
            choices[25661] = "Minolta AF 24mm F2.8";
            choices[25721] = "Minolta/Sony AF 500mm F8 Reflex";
            choices[25781] = "Minolta AF 16mm F2.8 Fisheye";
            choices[25791] = "Minolta AF 20mm F2.8";
            choices[25811] = "Minolta/Sony AF 100mm F2.8 Macro New";
            choices[25851] = "Beroflex 35-135mm F3.5-4.5";
            choices[25858] = "Minolta AF 35-105mm F3.5-4.5 New";
            choices[25881] = "Minolta AF 70-210mm F3.5-4.5";
            choices[25891] = "Minolta AF 80-200 F2.8 APO";
            choices[25911] = "Minolta AF 35mm F1.4";
            choices[25921] = "Minolta AF 85mm F1.4 G (D)";
            choices[25931] = "Minolta AF 200mm F2.8 G APO";
            choices[25941] = "Minolta AF 3x-1x F1.7-2.8 Macro";
            choices[25961] = "Minolta AF 28mm F2";
            choices[25971] = "Minolta AF 35mm F2";
            choices[25981] = "Minolta AF 100mm F2";
            choices[26041] = "Minolta AF 80-200mm F4.5-5.6";
            choices[26051] = "Minolta AF 35-80mm F4-5.6";
            choices[26061] = "Minolta AF 100-300mm F4.5-5.6 (D)";
            choices[26071] = "Minolta AF 35-80mm F4-5.6";
            choices[26081] = "Minolta AF 300mm F2.8 G";
            choices[26091] = "Minolta AF 600mm F4 HS-APO G";
            choices[26121] = "Minolta AF 200mm F2.8 G HS-APO";
            choices[26131] = "Minolta AF 50mm F1.7 New";
            choices[26151] = "Minolta AF 28-105mm F3.5-4.5 Power Zoom";
            choices[26161] = "Minolta AF 35-200mm F4.5-5.6 Power Zoom";
            choices[26181] = "Minolta AF 28-80mm F4-5.6 Power Zoom";
            choices[26191] = "Minolta AF 80-200mm F4.5-5.6 Power Zoom";
            choices[26201] = "Minolta AF 28-70mm F2.8 G";
            choices[26211] = "Minolta AF 100-300mm F4.5-5.6 Power Zoom";
            choices[26241] = "Minolta AF 35-80mm F4-5.6 Power Zoom";
            choices[26281] = "Minolta AF 80-200mm F2.8 G";
            choices[26291] = "Minolta AF 85mm F1.4 New";
            choices[26311] = "Minolta/Sony AF 100-300mm F4.5-5.6 APO";
            choices[26321] = "Minolta AF 24-50mm F4 New";
            choices[26381] = "Minolta AF 50mm F2.8 Macro New";
            choices[26391] = "Minolta AF 100mm F2.8 Macro";
            choices[26411] = "Minolta AF 20mm F2.8 New";
            choices[26421] = "Minolta AF 24mm F2.8 New";
            choices[26441] = "Minolta AF 100-400mm F4.5-6.7 APO";
            choices[26621] = "Minolta AF 50mm F1.4 New";
            choices[26671] = "Minolta AF 35mm F2 New";
            choices[26681] = "Minolta AF 28mm F2 New";
            choices[26721] = "Minolta AF 24-105mm F3.5-4.5 (D)";
            choices[45671] = "Tokina 70-210mm F4-5.6";
            choices[45741] = "Minolta AF 200mm F2.8 G x2";
            choices[45851] = "Tamron - SP AF 300 F2.8 LD IF";
        }
};
SALensIDInterpreter saLensIDInterpreter;

class MATeleconverterInterpreter : public ChoiceInterpreter {
    public:
        MATeleconverterInterpreter () {
            choices[0]     = "None ";
            choices[0x48]  = "Minolta AF 2x APO (D)";
            choices[0x50]  = "Minolta AF 2x APO II";
            choices[0x88]  = "Minolta AF 1.4x APO (D)";
            choices[0x90]  = "Minolta AF 1.4x APO II";
        }
};
MATeleconverterInterpreter maTeleconverterInterpreter;

class MAQualityInterpreter : public ChoiceInterpreter {
    public:
        MAQualityInterpreter () {
            choices[0]  = "Raw";
            choices[1]  = "Super Fine";
            choices[2]  = "Fine";
            choices[3]  = "Standard";
            choices[4]  = "Economy";
            choices[5]  = "Extra fine";
            choices[6]  = "RAW + JPEG";
            choices[7]  = "cRAW";
            choices[8]  = "cRAW + JPEG";
        }
};
MAQualityInterpreter maQualityInterpreter;

class MAImageSizeInterpreter : public ChoiceInterpreter {
    public:
        MAImageSizeInterpreter () {
            choices[1]  = "1600x1200";
            choices[2]  = "1280x960";
            choices[3]  = "640x480";
            choices[5]  = "2560x1920";
            choices[6]  = "2272x1704";
            choices[7]  = "2048x1536";
        }
};
MAImageSizeInterpreter maImageSizeInterpreter;

const TagAttrib minoltaAttribs[] = {
 0, 1, 0, 0, 0x0000, "MakerNoteVersion", &stdInterpreter,
 0, 1, 0, 0, 0x0001, "MinoltaCameraSettingsOld", &stdInterpreter,
 0, 1, 0, 0, 0x0003, "MinoltaCameraSettings", &stdInterpreter,
 0, 1, 0, 0, 0x0004, "MinoltaCameraSettings7D", &stdInterpreter,
 0, 1, 0, 0, 0x0018, "ImageStabilization", &stdInterpreter,
 0, 1, 0, 0, 0x0040, "CompressedImageSize", &stdInterpreter,
 1, 1, 0, 0, 0x0081, "PreviewImage", &stdInterpreter,
 1, 1, 0, 0, 0x0088, "PreviewImageStart", &stdInterpreter,
 1, 1, 0, 0, 0x0089, "PreviewImageLength", &stdInterpreter,
 0, 1, 0, 0, 0x0100, "SceneMode", &saSceneModeInterpreter,
 0, 1, 0, 0, 0x0101, "ColorMode", &saColorModeInterpreter,
 0, 1, 0, 0, 0x0102, "MinoltaQuality", &maQualityInterpreter,
 0, 1, 0, 0, 0x0103, "MinoltaImageSize", &maImageSizeInterpreter,
 0, 1, 0, 0, 0x0104, "FlashExposureComp", &stdInterpreter,
 0, 1, 0, 0, 0x0105, "Teleconverter", &maTeleconverterInterpreter,
 0, 1, 0, 0, 0x0107, "ImageStabilization", &saOnOffInterpreter,
 0, 1, 0, 0, 0x010a, "ZoneMatching", &saZoneMatchingInterpreter,
 0, 1, 0, 0, 0x010b, "ColorTemperature", &stdInterpreter,
 0, 1, 0, 0, 0x010c, "LensID", &saLensIDInterpreter,
 0, 1, 0, 0, 0x0113, "ImageStabilization", &saOnOffInterpreter,
 0, 1, 0, 0, 0x0114, "MinoltaCameraSettings", &stdInterpreter,
 1, 1, 0, 0, 0x0e00, "PrintIM", &stdInterpreter,
 0, 1, 0, 0, 0x0f00, "MinoltaCameraSettings2", &stdInterpreter,
 -1, 0, 0,  0, 0, "", NULL};

const TagAttrib sonyAttribs[] = {
 1, 1, 0, 0, 0x0e00, "PrintIM", &stdInterpreter,
 1, 1, 0, 0, 0x2001, "PreviewImage", &stdInterpreter,
 0, 1, 0, 0, 0xb020, "ColorReproduction", &stdInterpreter,
 0, 1, 0, 0, 0xb021, "ColorTemperature", &stdInterpreter,
 0, 1, 0, 0, 0xb023, "SceneMode", &saSceneModeInterpreter,
 0, 1, 0, 0, 0xb024, "ZoneMatching", &saZoneMatchingInterpreter,
 0, 1, 0, 0, 0xb025, "DynamicRangeOptimizer", &saDynamicRangeOptimizerInterpreter,
 0, 1, 0, 0, 0xb026, "ImageStabilization", &saOnOffInterpreter,
 0, 1, 0, 0, 0xb027, "LensID", &saLensIDInterpreter,
 0, 1, 0, minoltaAttribs, 0xb028, "MinoltaMakerNote", &stdInterpreter,
 0, 1, 0, 0, 0xb029, "ColorMode", &saColorModeInterpreter,
 0, 1, 0, 0, 0xb040, "Macro", &saOnOffInterpreter,
 0, 1, 0, 0, 0xb041, "ExposureMode", &saExposureModeInterpreter,
 0, 1, 0, 0, 0xb047, "Quality", &saQualityInterpreter,
 0, 1, 0, 0, 0xb04b, "AntiBlur", &saAntiBlurInterpreter,
 0, 1, 0, 0, 0xb04e, "LongExposureNoiseReduction", &saOnOffInterpreter,
 -1, 0, 0,  0, 0, "", NULL};

};
#endif

