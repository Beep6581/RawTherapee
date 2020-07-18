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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cstdio>
#include <cmath>

#include "rtexif.h"

using namespace std;

namespace rtexif
{

class CAOnOffInterpreter : public Interpreter
{
public:
    std::string toString (const Tag* t) const override
    {
        int n = t->toInt();

        if ( n == 0 ) {
            return "OFF";
        } else if ( n == 1) {
            return "ON";
        } else {
            return "undef";
        }
    }
};
CAOnOffInterpreter caOnOffInterpreter;

class CAIntSerNumInterpreter : public Interpreter
{
public:
    CAIntSerNumInterpreter () {}
    std::string toString (const Tag* t) const override
    {
        return "";
    }
};

CAIntSerNumInterpreter caIntSerNumInterpreter;

class CAApertureInterpreter : public Interpreter
{
public:
    CAApertureInterpreter () {}
    std::string toString (const Tag* t) const override
    {
        char buffer[32];
        double v = pow (2.0, t->toDouble() / 64.0);

        if ( v < 0. || v > 1000.) {
            return "undef";
        }

        sprintf (buffer, "%.1f", v );
        return buffer;
    }
};
CAApertureInterpreter caApertureInterpreter;

class CAMacroModeInterpreter : public ChoiceInterpreter<>
{
public:
    CAMacroModeInterpreter()
    {
        choices[1] = "Macro";
        choices[2] = "Normal";
    }
};
CAMacroModeInterpreter caMacroModeInterpreter;

class CASelfTimerInterpreter : public Interpreter
{
public:
    std::string toString (const Tag* t) const override
    {
        int sec = t->toInt (0, SHORT);

        if ( !sec ) {
            return "OFF";
        }

        char buffer[32];
        sprintf (buffer, "%.1fs %s", sec / 10., (sec & 0x4000) ? ",Custom" : "");
        return buffer;
    }
};
CASelfTimerInterpreter caSelfTimerInterpreter;

class CAQualityInterpreter : public ChoiceInterpreter<>
{
public:
    CAQualityInterpreter()
    {
        choices[1] = "Economy";
        choices[2] = "Normal";
        choices[3] = "Fine";
        choices[4] = "RAW";
        choices[5] = "Superfine";
    }
};
CAQualityInterpreter caQualityInterpreter;

class CAFlashModeInterpreter : public ChoiceInterpreter<>
{
public:
    CAFlashModeInterpreter()
    {
        choices[0] = "Off";
        choices[1] = "Auto";
        choices[2] = "On";
        choices[3] = "Red-eye reduction";
        choices[4] = "Slow-sync";
        choices[5] = "Red-eye reduction (Auto)";
        choices[6] = "Red-eye reduction (On)";
        choices[16] = "External flash";
    }
};
CAFlashModeInterpreter caFlashModeInterpreter;

class CAContinuousDriveInterpreter : public ChoiceInterpreter<>
{
public:
    CAContinuousDriveInterpreter()
    {
        choices[0] = "Single";
        choices[1] = "Continuous";
        choices[2] = "Movie";
        choices[3] = "Continuous, Speed Priority";
        choices[4] = "Continuous, Low";
        choices[5] = "Continuous, High";
        choices[6] = "Silent Single";
        choices[9] = "Single, Silent";
        choices[10] = "Continuous, Silent";
    }
};
CAContinuousDriveInterpreter caContinuousDriveInterpreter;

class CAFocusModeInterpreter : public ChoiceInterpreter<>
{
public:
    CAFocusModeInterpreter()
    {
        choices[0] = "One-shot AF";
        choices[1] = "AI Servo AF";
        choices[2] = "AI Focus AF";
        choices[3] = "Manual Focus (3)";
        choices[4] = "Single";
        choices[5] = "Continuous";
        choices[6] = "Manual Focus (6)";
        choices[16] = "Pan Focus";
        choices[256] = "AF + MF";
        choices[512] = "Movie Snap Focus";
        choices[519] = "Movie Servo AF";
    }
};
CAFocusModeInterpreter caFocusModeInterpreter;

class CARecordModeInterpreter : public ChoiceInterpreter<>
{
public:
    CARecordModeInterpreter()
    {
        choices[1] = "JPEG";
        choices[2] = "CRW+THM";
        choices[3] = "AVI+THM";
        choices[4] = "TIF";
        choices[5] = "TIF+JPEG";
        choices[6] = "CR2";
        choices[7] = "CR2+JPEG";
        choices[9] = "MOV";
        choices[10] = "MP4";
    }
};
CARecordModeInterpreter caRecordModeInterpreter;

class CAImageSizeInterpreter : public ChoiceInterpreter<>
{
public:
    CAImageSizeInterpreter ()
    {
        choices[0] = "Large";
        choices[1] = "Medium";
        choices[2] = "Small";
        choices[5] = "Medium 1";
        choices[6] = "Medium 2";
        choices[7] = "Medium 3";
        choices[8] = "Postcard";
        choices[9] = "Widescreen";
        choices[10] = "Medium Widescreen";
        choices[14] = "Small 1";
        choices[15] = "Small 2";
        choices[16] = "Small 3";
        choices[128] = "640x480 Movie";
        choices[129] = "Medium Movie";
        choices[130] = "Small Movie";
        choices[137] = "1280x720 Movie";
        choices[142] = "1920x1080 Movie";
    }
};
CAImageSizeInterpreter caImageSizeInterpreter;

class CAEasyModeInterpreter : public ChoiceInterpreter<>
{
public:
    CAEasyModeInterpreter ()
    {
        choices[0] = "Full auto";
        choices[1] = "Manual";
        choices[2] = "Landscape";
        choices[3] = "Fast shutter";
        choices[4] = "Slow shutter";
        choices[5] = "Night";
        choices[6] = "Gray Scale";
        choices[7] = "Sepia";
        choices[8] = "Portrait";
        choices[9] = "Sports";
        choices[10] = "Macro";
        choices[11] = "Black & White";
        choices[12] = "Pan focus";
        choices[13] = "Vivid";
        choices[14] = "Neutral";
        choices[15] = "Flash Off";
        choices[16] = "Long Shutter";
        choices[17] = "Super Macro";
        choices[18] = "Foliage";
        choices[19] = "Indoor";
        choices[20] = "Fireworks";
        choices[21] = "Beach";
        choices[22] = "Underwater";
        choices[23] = "Snow";
        choices[24] = "Kids & Pets";
        choices[25] = "Night Snapshot";
        choices[26] = "Digital Macro";
        choices[27] = "My Colors";
        choices[28] = "Movie Snap";
        choices[29] = "Super Macro 2";
        choices[30] = "Color Accent";
        choices[31] = "Color Swap";
        choices[32] = "Aquarium";
        choices[33] = "ISO 3200";
        choices[34] = "ISO 6400";
        choices[35] = "Creative Light Effect";
        choices[36] = "Easy";
        choices[37] = "Quick Shot";
        choices[38] = "Creative Auto";
        choices[39] = "Zoom Blur";
        choices[40] = "Low Light";
        choices[41] = "Nostalgic";
        choices[42] = "Super Vivid";
        choices[43] = "Poster Effect";
        choices[44] = "Face Self-timer";
        choices[45] = "Smile";
        choices[46] = "Wink Self-timer";
        choices[47] = "Fisheye Effect";
        choices[48] = "Miniature Effect";
        choices[49] = "High-speed Burst";
        choices[50] = "Best Image Selection";
        choices[51] = "High Dynamic Range";
        choices[52] = "Handheld Night Scene";
        choices[53] = "Movie Digest";
        choices[54] = "Live View Control";
        choices[55] = "Discreet";
        choices[56] = "Blur Reduction";
        choices[57] = "Monochrome";
        choices[58] = "Toy Camera Effect";
        choices[59] = "Scene Intelligent Auto";
        choices[60] = "High-speed Burst HQ";
        choices[61] = "Smooth Skin";
        choices[62] = "Soft Focus";
        choices[257] = "Spotlight";
        choices[258] = "Night 2";
        choices[259] = "Night+";
        choices[260] = "Super Night";
        choices[261] = "Sunset";
        choices[263] = "Night Scene";
        choices[264] = "Surface";
        choices[265] = "Low Light 2";
    }
};
CAEasyModeInterpreter caEasyModeInterpreter;

class CADigitalZoomInterpreter : public ChoiceInterpreter<>
{
public:
    CADigitalZoomInterpreter()
    {
        choices[0] = "None";
        choices[1] = "2x";
        choices[2] = "4x";
        choices[3] = "Other";
    }
};
CADigitalZoomInterpreter caDigitalZoomInterpreter;

class CAMeteringModeInterpreter : public ChoiceInterpreter<>
{
public:
    CAMeteringModeInterpreter()
    {
        choices[0] = "Default";
        choices[1] = "Spot";
        choices[2] = "Average";
        choices[3] = "Evaluative";
        choices[4] = "Partial";
        choices[5] = "Center-weighted average";
    }
};
CAMeteringModeInterpreter caMeteringModeInterpreter;

class CAFocusRangeInterpreter : public ChoiceInterpreter<>
{
public:
    CAFocusRangeInterpreter()
    {
        choices[0] = "Manual";
        choices[1] = "Auto";
        choices[2] = "Not Known";
        choices[3] = "Macro";
        choices[4] = "Very Close";
        choices[5] = "Close";
        choices[6] = "Middle Range";
        choices[7] = "Far Range";
        choices[8] = "Pan Focus";
        choices[9] = "Super Macro";
        choices[10] = "Infinity";
    }
};
CAFocusRangeInterpreter caFocusRangeInterpreter;

class CAAFPointInterpreter : public ChoiceInterpreter<>
{
public:
    CAAFPointInterpreter()
    {
        choices[0x2005] = "Manual AF point selection ";
        choices[0x3000] = "None (MF)";
        choices[0x3001] = "Auto AF point selection ";
        choices[0x3002] = "Right ";
        choices[0x3003] = "Center ";
        choices[0x3004] = "Left ";
        choices[0x4001] = "Auto AF point selection ";
        choices[0x4006] = "Face Detect";
    }
};
CAAFPointInterpreter caAFPointInterpreter;

class CAExposureModeInterpreter : public ChoiceInterpreter<>
{
public:
    CAExposureModeInterpreter()
    {
        choices[0] = "Easy";
        choices[1] = "Program AE";
        choices[2] = "Shutter speed priority AE";
        choices[3] = "Aperture-priority AE";
        choices[4] = "Manual";
        choices[5] = "Depth-of-field AE";
        choices[6] = "M-Dep";
        choices[7] = "Bulb";
    }
};
CAExposureModeInterpreter caExposureModeInterpreter;

class CAFlashBitsInterpreter : public Interpreter
{
public:
    std::string toString (const Tag* t) const override
    {
        std::ostringstream s;
        unsigned bits = t->toInt (0, SHORT);

        if ( bits & 0x0001 ) {
            s << "Manual ";
        }

        if ( bits & 0x0002 ) {
            s << "TTL ";
        }

        if ( bits & 0x0004 ) {
            s << "A-TTL ";
        }

        if ( bits & 0x0008 ) {
            s << "E-TTL ";
        }

        if ( bits & 0x0010 ) {
            s << "FP sync enabled ";
        }

        if ( bits & 0x0080 ) {
            s << "2nd curtain ";
        }

        if ( bits & 0x0800 ) {
            s << "FP sync used ";
        }

        if ( bits & 0x2000 ) {
            s << "Built-in ";
        }

        if ( bits & 0x4000 ) {
            s << "External ";
        }

        return s.str();
    }
};
CAFlashBitsInterpreter caFlashBitsInterpreter;

class CAFocusContinuousInterpreter : public ChoiceInterpreter<>
{
public:
    CAFocusContinuousInterpreter()
    {
        choices[0] = "Single";
        choices[1] = "Continuous";
        choices[8] = "Manual";
    }
};
CAFocusContinuousInterpreter caFocusContinuousInterpreter;

class CAAESettingsInterpreter : public ChoiceInterpreter<>
{
public:
    CAAESettingsInterpreter()
    {
        choices[0] = "Normal AE";
        choices[1] = "Exposure Compensation";
        choices[2] = "AE Lock";
        choices[3] = "AE Lock + Exposure Comp.";
        choices[4] = "No AE";
    }
};
CAAESettingsInterpreter caAESettingsInterpreter;

class CAStabilizationInterpreter : public ChoiceInterpreter<>
{
public:
    CAStabilizationInterpreter()
    {
        choices[0] = "Off";
        choices[1] = "On";
        choices[2] = "Shoot Only";
        choices[3] = "Panning";
        choices[4] = "Dynamic";
        choices[256] = "Off (2)";
        choices[257] = "On (2)";
        choices[258] = "Shoot Only (2)";
        choices[259] = "Panning (2)";
        choices[260] = "Dynamic (2)";
    }
};
CAStabilizationInterpreter caStabilizationInterpreter;

class CASpotMeteringInterpreter : public ChoiceInterpreter<>
{
public:
    CASpotMeteringInterpreter()
    {
        choices[0] = "Center";
        choices[1] = "AF Point";
    }
};
CASpotMeteringInterpreter caSpotMeteringInterpreter;

class CAPhotoEffectInterpreter : public ChoiceInterpreter<>
{
public:
    CAPhotoEffectInterpreter()
    {
        choices[0] = "Off";
        choices[1] = "Vivid";
        choices[2] = "Neutral";
        choices[3] = "Smooth";
        choices[4] = "Sepia";
        choices[5] = "B&W";
        choices[6] = "Custom";
        choices[100] = "My Color Data";
    }
};
CAPhotoEffectInterpreter caPhotoEffectInterpreter;

class CAManualFlashInterpreter : public ChoiceInterpreter<>
{
public:
    CAManualFlashInterpreter()
    {
        choices[0] = "N/A";
        choices[0x500] = "Full";
        choices[0x502] = "Medium";
        choices[0x504] = "Low";
        choices[0x7fff] = "N/A";
    }
};
CAManualFlashInterpreter caManualFlashInterpreter;

class CARAWQualityInterpreter : public ChoiceInterpreter<>
{
public:
    CARAWQualityInterpreter()
    {
        choices[0] = "N/A";
        choices[1] = "sRAW1 (mRAW)";
        choices[2] = "sRAW2 (sRAW)";
    }
};
CARAWQualityInterpreter caRAWQualityInterpreter;

class CAFocalInterpreter : public Interpreter
{
public:
    std::string toString (const Tag* t) const override
    {
        Tag *unitTag = t->getParent()->getRoot()->findTag ("FocalUnits");
        double v = unitTag ? unitTag->toDouble() : 1.;
        v = (v > 0. ? t->toDouble() / v : t->toDouble());

        if ( v < 0. || v > 1000000.) {
            return "undef";
        }

        char buffer[32];
        sprintf (buffer, "%.1f", v );
        return buffer;
    }
};
CAFocalInterpreter caFocalInterpreter;

class CALensInterpreter : public IntLensInterpreter< int >
{
public:
    CALensInterpreter ()
    {
        choices = {
            {1, "Canon EF 50mm f/1.8"},
            {2, "Canon EF 28mm f/2.8 or Sigma Lens"},
            {2, "Sigma 24mm f/2.8 Super Wide II"},
            {3, "Canon EF 135mm f/2.8 Soft"},
            {4, "Canon EF 35-105mm f/3.5-4.5 or Sigma Lens"},
            {4, "Sigma UC Zoom 35-135mm f/4-5.6"},
            {5, "Canon EF 35-70mm f/3.5-4.5"},
            {6, "Canon EF 28-70mm f/3.5-4.5 or Sigma or Tokina Lens"},
            {6, "Sigma 18-50mm f/3.5-5.6 DC"},
            {6, "Sigma 18-125mm f/3.5-5.6 DC IF ASP"},
            {6, "Tokina AF 193-2 19-35mm f/3.5-4.5"},
            {6, "Sigma 28-80mm f/3.5-5.6 II Macro"},
            {6, "Sigma 28-300mm f/3.5-6.3 DG Macro"},
            {7, "Canon EF 100-300mm f/5.6L"},
            {8, "Canon EF 100-300mm f/5.6 or Sigma or Tokina Lens"},
            {8, "Sigma 70-300mm f/4-5.6 [APO] DG Macro"},
            {8, "Tokina AT-X 242 AF 24-200mm f/3.5-5.6"},
            {9, "Canon EF 70-210mm f/4"},
            {9, "Sigma 55-200mm f/4-5.6 DC"},
            {10, "Canon EF 50mm f/2.5 Macro or Sigma Lens"},
            {10, "Sigma 50mm f/2.8 EX"},
            {10, "Sigma 28mm f/1.8"},
            {10, "Sigma 105mm f/2.8 Macro EX"},
            {10, "Sigma 70mm f/2.8 EX DG Macro EF"},
            {11, "Canon EF 35mm f/2"},
            {13, "Canon EF 15mm f/2.8 Fisheye"},
            {14, "Canon EF 50-200mm f/3.5-4.5L"},
            {15, "Canon EF 50-200mm f/3.5-4.5"},
            {16, "Canon EF 35-135mm f/3.5-4.5"},
            {17, "Canon EF 35-70mm f/3.5-4.5A"},
            {18, "Canon EF 28-70mm f/3.5-4.5"},
            {20, "Canon EF 100-200mm f/4.5A"},
            {21, "Canon EF 80-200mm f/2.8L"},
            {22, "Canon EF 20-35mm f/2.8L or Tokina Lens"},
            {22, "Tokina AT-X 280 AF Pro 28-80mm f/2.8 Aspherical"},
            {23, "Canon EF 35-105mm f/3.5-4.5"},
            {24, "Canon EF 35-80mm f/4-5.6 Power Zoom"},
            {25, "Canon EF 35-80mm f/4-5.6 Power Zoom"},
            {26, "Canon EF 100mm f/2.8 Macro or Other Lens"},
            {26, "Cosina 100mm f/3.5 Macro AF"},
            {26, "Tamron SP AF 90mm f/2.8 Di Macro"},
            {26, "Tamron SP AF 180mm f/3.5 Di Macro"},
            {26, "Carl Zeiss Planar T* 50mm f/1.4"},
            {26, "Voigtlander APO Lanthar 125mm F2.5 SL Macro"},
            {26, "Carl Zeiss Planar T 85mm f/1.4 ZE"},
            {27, "Canon EF 35-80mm f/4-5.6"},
            {28, "Canon EF 80-200mm f/4.5-5.6 or Tamron Lens"},
            {28, "Tamron SP AF 28-105mm f/2.8 LD Aspherical IF"},
            {28, "Tamron SP AF 28-75mm f/2.8 XR Di LD Aspherical [IF] Macro"},
            {28, "Tamron AF 70-300mm f/4-5.6 Di LD 1:2 Macro"},
            {28, "Tamron AF Aspherical 28-200mm f/3.8-5.6"},
            {29, "Canon EF 50mm f/1.8 II"},
            {30, "Canon EF 35-105mm f/4.5-5.6"},
            {31, "Canon EF 75-300mm f/4-5.6 or Tamron Lens"},
            {31, "Tamron SP AF 300mm f/2.8 LD IF"},
            {32, "Canon EF 24mm f/2.8 or Sigma Lens"},
            {32, "Sigma 15mm f/2.8 EX Fisheye"},
            {33, "Voigtlander or Carl Zeiss Lens"},
            {33, "Voigtlander Ultron 40mm f/2 SLII Aspherical"},
            {33, "Voigtlander Color Skopar 20mm f/3.5 SLII Aspherical"},
            {33, "Voigtlander APO-Lanthar 90mm f/3.5 SLII Close Focus"},
            {33, "Carl Zeiss Distagon T* 15mm f/2.8 ZE"},
            {33, "Carl Zeiss Distagon T* 18mm f/3.5 ZE"},
            {33, "Carl Zeiss Distagon T* 21mm f/2.8 ZE"},
            {33, "Carl Zeiss Distagon T* 25mm f/2 ZE"},
            {33, "Carl Zeiss Distagon T* 28mm f/2 ZE"},
            {33, "Carl Zeiss Distagon T* 35mm f/2 ZE"},
            {33, "Carl Zeiss Distagon T* 35mm f/1.4 ZE"},
            {33, "Carl Zeiss Planar T* 50mm f/1.4 ZE"},
            {33, "Carl Zeiss Makro-Planar T* 50mm f/2 ZE"},
            {33, "Carl Zeiss Makro-Planar T* 100mm f/2 ZE"},
            {33, "Carl Zeiss Apo-Sonnar T* 135mm f/2 ZE"},
            {35, "Canon EF 35-80mm f/4-5.6"},
            {36, "Canon EF 38-76mm f/4.5-5.6"},
            {37, "Canon EF 35-80mm f/4-5.6 or Tamron Lens"},
            {37, "Tamron 70-200mm f/2.8 Di LD IF Macro"},
            {37, "Tamron AF 28-300mm f/3.5-6.3 XR Di VC LD Aspherical [IF] Macro (A20)"},
            {37, "Tamron SP AF 17-50mm f/2.8 XR Di II VC LD Aspherical [IF]"},
            {37, "Tamron AF 18-270mm f/3.5-6.3 Di II VC LD Aspherical [IF] Macro"},
            {38, "Canon EF 80-200mm f/4.5-5.6"},
            {39, "Canon EF 75-300mm f/4-5.6"},
            {40, "Canon EF 28-80mm f/3.5-5.6"},
            {41, "Canon EF 28-90mm f/4-5.6"},
            {42, "Canon EF 28-200mm f/3.5-5.6 or Tamron Lens"},
            {42, "Tamron AF 28-300mm f/3.5-6.3 XR Di VC LD Aspherical [IF] Macro (A20)"},
            {43, "Canon EF 28-105mm f/4-5.6"},
            {44, "Canon EF 90-300mm f/4.5-5.6"},
            {45, "Canon EF-S 18-55mm f/3.5-5.6 [II]"},
            {46, "Canon EF 28-90mm f/4-5.6"},
            {47, "Zeiss Milvus 35mm f/2 or 50mm f/2"},
            {47, "Zeiss Milvus 50mm f/2 Makro"},
            {47, "Zeiss Milvus 135mm f/2 ZE"},
            {48, "Canon EF-S 18-55mm f/3.5-5.6 IS"},
            {49, "Canon EF-S 55-250mm f/4-5.6 IS"},
            {50, "Canon EF-S 18-200mm f/3.5-5.6 IS"},
            {51, "Canon EF-S 18-135mm f/3.5-5.6 IS"},
            {52, "Canon EF-S 18-55mm f/3.5-5.6 IS II"},
            {53, "Canon EF-S 18-55mm f/3.5-5.6 III"},
            {54, "Canon EF-S 55-250mm f/4-5.6 IS II"},
            {60, "Irix 11mm f/4"},
            {80, "Canon TS-E 50mm f/2.8L Macro"},
            {81, "Canon TS-E 90mm f/2.8L Macro"},
            {82, "Canon TS-E 135mm f/4L Macro"},
            {94, "Canon TS-E 17mm f/4L"},
            {95, "Canon TS-E 24mm f/3.5L II"},
            {103, "Samyang AF 14mm f/2.8 EF or Rokinon Lens"},
            {103, "Rokinon SP 14mm f/2.4"},
            {103, "Rokinon AF 14mm f/2.8 EF"},
            {106, "Rokinon SP / Samyang XP 35mm f/1.2"},
            {112, "Sigma 28mm f/1.5 FF High-speed Prime or other Sigma Lens"},
            {112, "Sigma 40mm f/1.5 FF High-speed Prime"},
            {112, "Sigma 105mm f/1.5 FF High-speed Prime"},
            {117, "Tamron 35-150mm f/2.8-4.0 Di VC OSD (A043) or other Tamron Lens"},
            {117, "Tamron SP 35mm f/1.4 Di USD (F045)"},
            {124, "Canon MP-E 65mm f/2.8 1-5x Macro Photo"},
            {125, "Canon TS-E 24mm f/3.5L"},
            {126, "Canon TS-E 45mm f/2.8"},
            {127, "Canon TS-E 90mm f/2.8 or Tamron Lens"},
            {127, "Tamron 18-200mm f/3.5-6.3 Di II VC (B018)"},
            {129, "Canon EF 300mm f/2.8L USM"},
            {130, "Canon EF 50mm f/1.0L USM"},
            {131, "Canon EF 28-80mm f/2.8-4L USM or Sigma Lens"},
            {131, "Sigma 8mm f/3.5 EX DG Circular Fisheye"},
            {131, "Sigma 17-35mm f/2.8-4 EX DG Aspherical HSM"},
            {131, "Sigma 17-70mm f/2.8-4.5 DC Macro"},
            {131, "Sigma APO 50-150mm f/2.8 [II] EX DC HSM"},
            {131, "Sigma APO 120-300mm f/2.8 EX DG HSM"},
            {131, "Sigma 4.5mm f/2.8 EX DC HSM Circular Fisheye"},
            {131, "Sigma 70-200mm f/2.8 APO EX HSM"},
            {131, "Sigma 28-70mm f/2.8-4 DG"},
            {132, "Canon EF 1200mm f/5.6L USM"},
            {134, "Canon EF 600mm f/4L IS USM"},
            {135, "Canon EF 200mm f/1.8L USM"},
            {136, "Canon EF 300mm f/2.8L USM"},
            {136, "Tamron SP 15-30mm f/2.8 Di VC USD (A012)"},
            {137, "Canon EF 85mm f/1.2L USM or Sigma or Tamron Lens"},
            {137, "Sigma 18-50mm f/2.8-4.5 DC OS HSM"},
            {137, "Sigma 50-200mm f/4-5.6 DC OS HSM"},
            {137, "Sigma 18-250mm f/3.5-6.3 DC OS HSM"},
            {137, "Sigma 24-70mm f/2.8 IF EX DG HSM"},
            {137, "Sigma 18-125mm f/3.8-5.6 DC OS HSM"},
            {137, "Sigma 17-70mm f/2.8-4 DC Macro OS HSM | C"},
            {137, "Sigma 17-50mm f/2.8 OS HSM"},
            {137, "Sigma 18-200mm f/3.5-6.3 DC OS HSM [II]"},
            {137, "Tamron AF 18-270mm f/3.5-6.3 Di II VC PZD (B008)"},
            {137, "Sigma 8-16mm f/4.5-5.6 DC HSM"},
            {137, "Tamron SP 17-50mm f/2.8 XR Di II VC (B005)"},
            {137, "Tamron SP 60mm f/2 Macro Di II (G005)"},
            {137, "Sigma 10-20mm f/3.5 EX DC HSM"},
            {137, "Tamron SP 24-70mm f/2.8 Di VC USD"},
            {137, "Sigma 18-35mm f/1.8 DC HSM"},
            {137, "Sigma 12-24mm f/4.5-5.6 DG HSM II"},
            {137, "Sigma 70-300mm f/4-5.6 DG OS"},
            {138, "Canon EF 28-80mm f/2.8-4L"},
            {139, "Canon EF 400mm f/2.8L USM"},
            {140, "Canon EF 500mm f/4.5L USM"},
            {141, "Canon EF 500mm f/4.5L USM"},
            {142, "Canon EF 300mm f/2.8L IS USM"},
            {143, "Canon EF 500mm f/4L IS USM or Sigma Lens"},
            {143, "Sigma 17-70mm f/2.8-4 DC Macro OS HSM"},
            {144, "Canon EF 35-135mm f/4-5.6 USM"},
            {145, "Canon EF 100-300mm f/4.5-5.6 USM"},
            {146, "Canon EF 70-210mm f/3.5-4.5 USM"},
            {147, "Canon EF 35-135mm f/4-5.6 USM"},
            {148, "Canon EF 28-80mm f/3.5-5.6 USM"},
            {149, "Canon EF 100mm f/2 USM"},
            {150, "Canon EF 14mm f/2.8L USM or Sigma Lens"},
            {150, "Sigma 20mm EX f/1.8"},
            {150, "Sigma 30mm f/1.4 DC HSM"},
            {150, "Sigma 24mm f/1.8 DG Macro EX"},
            {150, "Sigma 28mm f/1.8 DG Macro EX"},
            {150, "Sigma 18-35mm f/1.8 DC HSM | A"},
            {151, "Canon EF 200mm f/2.8L USM"},
            {152, "Canon EF 300mm f/4L IS USM or Sigma Lens"},
            {152, "Sigma 12-24mm f/4.5-5.6 EX DG ASPHERICAL HSM"},
            {152, "Sigma 14mm f/2.8 EX Aspherical HSM"},
            {152, "Sigma 10-20mm f/4-5.6"},
            {152, "Sigma 100-300mm f/4"},
            {152, "Sigma 300-800mm f/5.6 APO EX DG HSM"},
            {153, "Canon EF 35-350mm f/3.5-5.6L USM or Sigma or Tamron Lens"},
            {153, "Sigma 50-500mm f/4-6.3 APO HSM EX"},
            {153, "Tamron AF 28-300mm f/3.5-6.3 XR LD Aspherical [IF] Macro"},
            {153, "Tamron AF 18-200mm f/3.5-6.3 XR Di II LD Aspherical [IF] Macro (A14)"},
            {153, "Tamron 18-250mm f/3.5-6.3 Di II LD Aspherical [IF] Macro"},
            {154, "Canon EF 20mm f/2.8 USM or Zeiss Lens"},
            {154, "Zeiss Milvus 21mm f/2.8"},
            {154, "Zeiss Milvus 15mm f/2.8 ZE"},
            {154, "Zeiss Milvus 18mm f/2.8 ZE"},
            {155, "Canon EF 85mm f/1.8 USM or Sigma Lens"},
            {155, "Sigma 14mm f/1.8 DG HSM | A"},
            {156, "Canon EF 28-105mm f/3.5-4.5 USM or Tamron Lens"},
            {156, "Tamron SP 70-300mm f/4-5.6 Di VC USD (A005)"},
            {156, "Tamron SP AF 28-105mm f/2.8 LD Aspherical IF (176D)"},
            {160, "Canon EF 20-35mm f/3.5-4.5 USM or Tamron or Tokina Lens"},
            {160, "Tamron AF 19-35mm f/3.5-4.5"},
            {160, "Tokina AT-X 124 AF Pro DX 12-24mm f/4"},
            {160, "Tokina AT-X 107 AF DX 10-17mm f/3.5-4.5 Fisheye"},
            {160, "Tokina AT-X 116 AF Pro DX 11-16mm f/2.8"},
            {160, "Tokina AT-X 11-20 F2.8 PRO DX Aspherical 11-20mm f/2.8"},
            {161, "Canon EF 28-70mm f/2.8L USM or Other Lens"},
            {161, "Sigma 24-70mm f/2.8 EX"},
            {161, "Sigma 28-70mm f/2.8 EX"},
            {161, "Sigma 24-60mm f/2.8 EX DG"},
            {161, "Tamron AF 17-50mm f/2.8 Di-II LD Aspherical"},
            {161, "Tamron 90mm f/2.8"},
            {161, "Tamron SP AF 17-35mm f/2.8-4 Di LD Aspherical IF (A05)"},
            {161, "Tamron SP AF 28-75mm f/2.8 XR Di LD Aspherical [IF] Macro"},
            {161, "Tokina AT-X 24-70mm f/2.8 PRO FX (IF)"},
            {162, "Canon EF 200mm f/2.8L USM"},
            {163, "Canon EF 300mm f/4L"},
            {164, "Canon EF 400mm f/5.6L"},
            {165, "Canon EF 70-200mm f/2.8L USM"},
            {166, "Canon EF 70-200mm f/2.8L USM + 1.4x"},
            {167, "Canon EF 70-200mm f/2.8L USM + 2x"},
            {168, "Canon EF 28mm f/1.8 USM or Sigma Lens"},
            {168, "Sigma 50-100mm f/1.8 DC HSM | A"},
            {169, "Canon EF 17-35mm f/2.8L USM or Sigma Lens"},
            {169, "Sigma 18-200mm f/3.5-6.3 DC OS"},
            {169, "Sigma 15-30mm f/3.5-4.5 EX DG Aspherical"},
            {169, "Sigma 18-50mm f/2.8 Macro"},
            {169, "Sigma 50mm f/1.4 EX DG HSM"},
            {169, "Sigma 85mm f/1.4 EX DG HSM"},
            {169, "Sigma 30mm f/1.4 EX DC HSM"},
            {169, "Sigma 35mm f/1.4 DG HSM"},
            {169, "Sigma 35mm f/1.5 FF High-Speed Prime | 017"},
            {169, "Sigma 70mm f/2.8 Macro EX DG"},
            {170, "Canon EF 200mm f/2.8L II USM or Sigma Lens"},
            {170, "Sigma 300mm f/2.8 APO EX DG HSM"},
            {170, "Sigma 800mm f/5.6 APO EX DG HSM"},
            {171, "Canon EF 300mm f/4L USM"},
            {172, "Canon EF 400mm f/5.6L USM or Sigma Lens"},
            {172, "Sigma 150-600mm f/5-6.3 DG OS HSM | S"},
            {172, "Sigma 500mm f/4.5 APO EX DG HSM"},
            {173, "Canon EF 180mm Macro f/3.5L USM or Sigma Lens"},
            {173, "Sigma 180mm EX HSM Macro f/3.5"},
            {173, "Sigma APO Macro 150mm f/2.8 EX DG HSM"},
            {173, "Sigma 10mm f/2.8 EX DC Fisheye"},
            {173, "Sigma 15mm f/2.8 EX DG Diagonal Fisheye"},
            {173, "Venus Laowa 100mm F2.8 2X Ultra Macro APO"},
            {174, "Canon EF 135mm f/2L USM or Other Lens"},
            {174, "Sigma 70-200mm f/2.8 EX DG APO OS HSM"},
            {174, "Sigma 50-500mm f/4.5-6.3 APO DG OS HSM"},
            {174, "Sigma 150-500mm f/5-6.3 APO DG OS HSM"},
            {174, "Zeiss Milvus 100mm f/2 Makro"},
            {174, "Sigma APO 50-150mm f/2.8 EX DC OS HSM"},
            {174, "Sigma APO 120-300mm f/2.8 EX DG OS HSM"},
            {174, "Sigma 120-300mm f/2.8 DG OS HSM S013"},
            {174, "Sigma 120-400mm f/4.5-5.6 APO DG OS HSM"},
            {174, "Sigma 200-500mm f/2.8 APO EX DG"},
            {175, "Canon EF 400mm f/2.8L USM"},
            {176, "Canon EF 24-85mm f/3.5-4.5 USM"},
            {177, "Canon EF 300mm f/4L IS USM"},
            {178, "Canon EF 28-135mm f/3.5-5.6 IS"},
            {179, "Canon EF 24mm f/1.4L USM"},
            {180, "Canon EF 35mm f/1.4L USM or Other Lens"},
            {180, "Sigma 50mm f/1.4 DG HSM | A"},
            {180, "Sigma 24mm f/1.4 DG HSM | A"},
            {180, "Zeiss Milvus 50mm f/1.4"},
            {180, "Zeiss Milvus 85mm f/1.4"},
            {180, "Zeiss Otus 28mm f/1.4 ZE"},
            {180, "Sigma 24mm f/1.5 FF High-Speed Prime | 017"},
            {180, "Sigma 50mm f/1.5 FF High-Speed Prime | 017"},
            {180, "Sigma 85mm f/1.5 FF High-Speed Prime | 017"},
            {180, "Tokina Opera 50mm f/1.4 FF"},
            {180, "Sigma 20mm f/1.4 DG HSM | A"},
            {181, "Canon EF 100-400mm f/4.5-5.6L IS USM + 1.4x or Sigma Lens"},
            {181, "Sigma 150-600mm f/5-6.3 DG OS HSM | S + 1.4x"},
            {182, "Canon EF 100-400mm f/4.5-5.6L IS USM + 2x or Sigma Lens"},
            {182, "Sigma 150-600mm f/5-6.3 DG OS HSM | S + 2x"},
            {183, "Canon EF 100-400mm f/4.5-5.6L IS USM or Sigma Lens"},
            {183, "Sigma 150mm f/2.8 EX DG OS HSM APO Macro"},
            {183, "Sigma 105mm f/2.8 EX DG OS HSM Macro"},
            {183, "Sigma 180mm f/2.8 EX DG OS HSM APO Macro"},
            {183, "Sigma 150-600mm f/5-6.3 DG OS HSM | C"},
            {183, "Sigma 150-600mm f/5-6.3 DG OS HSM | S"},
            {183, "Sigma 100-400mm f/5-6.3 DG OS HSM"},
            {183, "Sigma 180mm f/3.5 APO Macro EX DG IF HSM"},
            {184, "Canon EF 400mm f/2.8L USM + 2x"},
            {185, "Canon EF 600mm f/4L IS USM"},
            {186, "Canon EF 70-200mm f/4L USM"},
            {187, "Canon EF 70-200mm f/4L USM + 1.4x"},
            {188, "Canon EF 70-200mm f/4L USM + 2x"},
            {189, "Canon EF 70-200mm f/4L USM + 2.8x"},
            {190, "Canon EF 100mm f/2.8 Macro USM"},
            {191, "Canon EF 400mm f/4 DO IS or Sigma Lens"},
            {191, "Sigma 500mm f/4 DG OS HSM"},
            {193, "Canon EF 35-80mm f/4-5.6 USM"},
            {194, "Canon EF 80-200mm f/4.5-5.6 USM"},
            {195, "Canon EF 35-105mm f/4.5-5.6 USM"},
            {196, "Canon EF 75-300mm f/4-5.6 USM"},
            {197, "Canon EF 75-300mm f/4-5.6 IS USM or Sigma Lens"},
            {197, "Sigma 18-300mm f/3.5-6.3 DC Macro OS HSM"},
            {198, "Canon EF 50mm f/1.4 USM or Other Lens"},
            {198, "Zeiss Otus 55mm f/1.4 ZE"},
            {198, "Zeiss Otus 85mm f/1.4 ZE"},
            {198, "Zeiss Milvus 25mm f/1.4"},
            {198, "Zeiss Otus 100mm f/1.4"},
            {198, "Zeiss Milvus 35mm f/1.4 ZE"},
            {198, "Yongnuo YN 35mm f/2"},
            {199, "Canon EF 28-80mm f/3.5-5.6 USM"},
            {200, "Canon EF 75-300mm f/4-5.6 USM"},
            {201, "Canon EF 28-80mm f/3.5-5.6 USM"},
            {202, "Canon EF 28-80mm f/3.5-5.6 USM IV"},
            {208, "Canon EF 22-55mm f/4-5.6 USM"},
            {209, "Canon EF 55-200mm f/4.5-5.6"},
            {210, "Canon EF 28-90mm f/4-5.6 USM"},
            {211, "Canon EF 28-200mm f/3.5-5.6 USM"},
            {212, "Canon EF 28-105mm f/4-5.6 USM"},
            {213, "Canon EF 90-300mm f/4.5-5.6 USM or Tamron Lens"},
            {213, "Tamron SP 150-600mm f/5-6.3 Di VC USD (A011)"},
            {213, "Tamron 16-300mm f/3.5-6.3 Di II VC PZD Macro (B016)"},
            {213, "Tamron SP 35mm f/1.8 Di VC USD (F012)"},
            {213, "Tamron SP 45mm f/1.8 Di VC USD (F013)"},
            {214, "Canon EF-S 18-55mm f/3.5-5.6 USM"},
            {215, "Canon EF 55-200mm f/4.5-5.6 II USM"},
            {217, "Tamron AF 18-270mm f/3.5-6.3 Di II VC PZD"},
            {220, "Yongnuo YN 50mm f/1.8"},
            {224, "Canon EF 70-200mm f/2.8L IS USM"},
            {225, "Canon EF 70-200mm f/2.8L IS USM + 1.4x"},
            {226, "Canon EF 70-200mm f/2.8L IS USM + 2x"},
            {227, "Canon EF 70-200mm f/2.8L IS USM + 2.8x"},
            {228, "Canon EF 28-105mm f/3.5-4.5 USM"},
            {229, "Canon EF 16-35mm f/2.8L USM"},
            {230, "Canon EF 24-70mm f/2.8L USM"},
            {231, "Canon EF 17-40mm f/4L USM or Sigma Lens"},
            {231, "Sigma 12-24mm f/4 DG HSM A016"},
            {232, "Canon EF 70-300mm f/4.5-5.6 DO IS USM"},
            {233, "Canon EF 28-300mm f/3.5-5.6L IS USM"},
            {234, "Canon EF-S 17-85mm f/4-5.6 IS USM or Tokina Lens"},
            {234, "Tokina AT-X 12-28 PRO DX 12-28mm f/4"},
            {235, "Canon EF-S 10-22mm f/3.5-4.5 USM"},
            {236, "Canon EF-S 60mm f/2.8 Macro USM"},
            {237, "Canon EF 24-105mm f/4L IS USM"},
            {238, "Canon EF 70-300mm f/4-5.6 IS USM"},
            {239, "Canon EF 85mm f/1.2L II USM or Rokinon Lens"},
            {239, "Rokinon SP 85mm f/1.2"},
            {240, "Canon EF-S 17-55mm f/2.8 IS USM or Sigma Lens"},
            {240, "Sigma 17-50mm f/2.8 EX DC OS HSM"},
            {241, "Canon EF 50mm f/1.2L USM"},
            {242, "Canon EF 70-200mm f/4L IS USM"},
            {243, "Canon EF 70-200mm f/4L IS USM + 1.4x"},
            {244, "Canon EF 70-200mm f/4L IS USM + 2x"},
            {245, "Canon EF 70-200mm f/4L IS USM + 2.8x"},
            {246, "Canon EF 16-35mm f/2.8L II USM"},
            {247, "Canon EF 14mm f/2.8L II USM"},
            {248, "Canon EF 200mm f/2L IS USM or Sigma Lens"},
            {248, "Sigma 24-35mm f/2 DG HSM | A"},
            {248, "Sigma 135mm f/2 FF High-Speed Prime | 017"},
            {248, "Sigma 24-35mm f/2.2 FF Zoom | 017"},
            {248, "Sigma 135mm f/1.8 DG HSM A017"},
            {249, "Canon EF 800mm f/5.6L IS USM"},
            {250, "Canon EF 24mm f/1.4L II USM or Sigma Lens"},
            {250, "Sigma 20mm f/1.4 DG HSM | A"},
            {250, "Sigma 20mm f/1.5 FF High-Speed Prime | 017"},
            {250, "Tokina Opera 16-28mm f/2.8 FF"},
            {250, "Sigma 85mm f/1.4 DG HSM A016"},
            {251, "Canon EF 70-200mm f/2.8L IS II USM"},
            {251, "Canon EF 70-200mm f/2.8L IS III USM"},
            {252, "Canon EF 70-200mm f/2.8L IS II USM + 1.4x"},
            {252, "Canon EF 70-200mm f/2.8L IS III USM + 1.4x"},
            {253, "Canon EF 70-200mm f/2.8L IS II USM + 2x"},
            {253, "Canon EF 70-200mm f/2.8L IS III USM + 2x"},
            {254, "Canon EF 100mm f/2.8L Macro IS USM"},
            {255, "Sigma 24-105mm f/4 DG OS HSM | A or Other Lens"},
            {255, "Sigma 180mm f/2.8 EX DG OS HSM APO Macro"},
            {255, "Tamron SP 70-200mm f/2.8 Di VC USD"},
            {368, "Sigma 14-24mm f/2.8 DG HSM | A or other Sigma Lens"},
            {368, "Sigma 35mm f/1.4 DG HSM | A"},
            {368, "Sigma 50mm f/1.4 DG HSM | A"},
            {368, "Sigma 40mm f/1.4 DG HSM | A"},
            {368, "Sigma 60-600mm f/4.5-6.3 DG OS HSM | S"},
            {368, "Sigma 28mm f/1.4 DG HSM | A"},
            {368, "Sigma 150-600mm f/5-6.3 DG OS HSM | S"},
            {368, "Sigma 85mm f/1.4 DG HSM | A"},
            {368, "Sigma 105mm f/1.4 DG HSM"},
            {368, "Sigma 14-24mm f/2.8 DG HSM"},
            {368, "Sigma 70mm f/2.8 DG Macro"},
            {488, "Canon EF-S 15-85mm f/3.5-5.6 IS USM"},
            {489, "Canon EF 70-300mm f/4-5.6L IS USM"},
            {490, "Canon EF 8-15mm f/4L Fisheye USM"},
            {491, "Canon EF 300mm f/2.8L IS II USM or Tamron Lens"},
            {491, "Tamron SP 70-200mm f/2.8 Di VC USD G2 (A025)"},
            {491, "Tamron 18-400mm f/3.5-6.3 Di II VC HLD (B028)"},
            {491, "Tamron 100-400mm f/4.5-6.3 Di VC USD (A035)"},
            {491, "Tamron 70-210mm f/4 Di VC USD (A034)"},
            {491, "Tamron 70-210mm f/4 Di VC USD (A034) + 1.4x"},
            {491, "Tamron SP 24-70mm f/2.8 Di VC USD G2 (A032)"},
            {492, "Canon EF 400mm f/2.8L IS II USM"},
            {493, "Canon EF 500mm f/4L IS II USM or EF 24-105mm f4L IS USM"},
            {493, "Canon EF 24-105mm f/4L IS USM"},
            {494, "Canon EF 600mm f/4L IS II USM"},
            {495, "Canon EF 24-70mm f/2.8L II USM or Sigma Lens"},
            {495, "Sigma 24-70mm f/2.8 DG OS HSM | A"},
            {496, "Canon EF 200-400mm f/4L IS USM"},
            {499, "Canon EF 200-400mm f/4L IS USM + 1.4x"},
            {502, "Canon EF 28mm f/2.8 IS USM or Tamron Lens"},
            {502, "Tamron 35mm f/1.8 Di VC USD (F012)"},
            {503, "Canon EF 24mm f/2.8 IS USM"},
            {504, "Canon EF 24-70mm f/4L IS USM"},
            {505, "Canon EF 35mm f/2 IS USM"},
            {506, "Canon EF 400mm f/4 DO IS II USM"},
            {507, "Canon EF 16-35mm f/4L IS USM"},
            {508, "Canon EF 11-24mm f/4L USM or Tamron Lens"},
            {508, "Tamron 10-24mm f/3.5-4.5 Di II VC HLD (B023)"},
            {624, "Sigma 70-200mm f/2.8 DG OS HSM | S"},
            {747, "Canon EF 100-400mm f/4.5-5.6L IS II USM or Tamron Lens"},
            {747, "Tamron SP 150-600mm f/5-6.3 Di VC USD G2"},
            {748, "Canon EF 100-400mm f/4.5-5.6L IS II USM + 1.4x or Tamron Lens"},
            {748, "Tamron 100-400mm f/4.5-6.3 Di VC USD A035E + 1.4x"},
            {748, "Tamron 70-210mm f/4 Di VC USD (A034) + 2x"},
            {749, "Tamron 100-400mm f/4.5-6.3 Di VC USD A035E + 2x"},
            {750, "Canon EF 35mm f/1.4L II USM or Tamron Lens"},
            {750, "Tamron SP 85mm f/1.8 Di VC USD (F016)"},
            {750, "Tamron SP 45mm f/1.8 Di VC USD (F013)"},
            {751, "Canon EF 16-35mm f/2.8L III USM"},
            {752, "Canon EF 24-105mm f/4L IS II USM"},
            {753, "Canon EF 85mm f/1.4L IS USM"},
            {754, "Canon EF 70-200mm f/4L IS II USM"},
            {757, "Canon EF 400mm f/2.8L IS III USM"},
            {758, "Canon EF 600mm f/4L IS III USM"},
            {1136, "Sigma 24-70mm f/2.8 DG OS HSM | A"},
            {4142, "Canon EF-S 18-135mm f/3.5-5.6 IS STM"},
            {4143, "Canon EF-M 18-55mm f/3.5-5.6 IS STM or Tamron Lens"},
            {4143, "Tamron 18-200mm f/3.5-6.3 Di III VC"},
            {4144, "Canon EF 40mm f/2.8 STM"},
            {4145, "Canon EF-M 22mm f/2 STM"},
            {4146, "Canon EF-S 18-55mm f/3.5-5.6 IS STM"},
            {4147, "Canon EF-M 11-22mm f/4-5.6 IS STM"},
            {4148, "Canon EF-S 55-250mm f/4-5.6 IS STM"},
            {4149, "Canon EF-M 55-200mm f/4.5-6.3 IS STM"},
            {4150, "Canon EF-S 10-18mm f/4.5-5.6 IS STM"},
            {4152, "Canon EF 24-105mm f/3.5-5.6 IS STM"},
            {4153, "Canon EF-M 15-45mm f/3.5-6.3 IS STM"},
            {4154, "Canon EF-S 24mm f/2.8 STM"},
            {4155, "Canon EF-M 28mm f/3.5 Macro IS STM"},
            {4156, "Canon EF 50mm f/1.8 STM"},
            {4157, "Canon EF-M 18-150mm f/3.5-6.3 IS STM"},
            {4158, "Canon EF-S 18-55mm f/4-5.6 IS STM"},
            {4159, "Canon EF-M 32mm f/1.4 STM"},
            {4160, "Canon EF-S 35mm f/2.8 Macro IS STM"},
            {4208, "Sigma 56mm f/1.4 DC DN | C"},
            {36910, "Canon EF 70-300mm f/4-5.6 IS II USM"},
            {36912, "Canon EF-S 18-135mm f/3.5-5.6 IS USM"},
            {61182, "Canon RF 35mm F1.8 Macro IS STM or other Canon RF Lens"},
            {61182, "Canon RF 50mm F1.2 L USM"},
            {61182, "Canon RF 24-105mm F4 L IS USM"},
            {61182, "Canon RF 28-70mm F2 L USM"},
            {61182, "Canon RF 85mm F1.2L USM"},
            {61182, "Canon RF 24-240mm F4-6.3 IS USM"},
            {61182, "Canon RF 24-70mm F2.8 L IS USM"},
            {61182, "Canon RF 15-35mm F2.8 L IS USM"},
            {61491, "Canon CN-E 14mm T3.1 L F"},
            {61492, "Canon CN-E 24mm T1.5 L F"},
            {61494, "Canon CN-E 85mm T1.3 L F"},
            {61495, "Canon CN-E 135mm T2.2 L F"},
            {61496, "Canon CN-E 35mm T1.5 L F"},
            {65535, "n/a"}
        };
    }

    std::string toString (const Tag* t) const override
    {
        int lensID = t->toInt();

        it_t r;
        size_t nFound = choices.count ( lensID );

        if (1 == nFound) {
            r = choices.find ( lensID );
            return r->second;
        }

        Tag *apertureTag = t->getParent()->getRoot()->findTag ("MaxAperture");
        Tag *focalLengthTag = t->getParent()->getRoot()->findTag ("FocalLength");
        Tag *focalLengthMaxTag = t->getParent()->getRoot()->findTag ("LongFocal");
        Tag *focalLengthMinTag = t->getParent()->getRoot()->findTag ("ShortFocal");
        Tag *unitTag = t->getParent()->getRoot()->findTag ("FocalUnits");
        double maxApertureAtFocal = 0.;
        double focalLength = 0.;
        double focalLengthMin = 0.;
        double focalLengthMax = 0.;

        if ( apertureTag ) {
            maxApertureAtFocal = pow (2.0, apertureTag->toDouble() / 64.0);
        }

        if ( unitTag ) {
            double unit = unitTag->toDouble();

            if ( unit == 0. ) {
                unit = 1;
            }

            if ( focalLengthTag ) {
                focalLength = focalLengthTag->toDouble();
            }

            if ( focalLengthMinTag ) {
                focalLengthMin = focalLengthMinTag->toDouble() / unit;
            }

            if ( focalLengthMaxTag ) {
                focalLengthMax = focalLengthMaxTag->toDouble() / unit;
            }
        }

        std::ostringstream s;
        s << "Unknown ";

        if (focalLengthMin > 0.) {
            s << focalLengthMin;
        }

        if (focalLengthMax > 0. && focalLengthMax != focalLengthMin) {
            s << "-" << focalLengthMax;
        }

        if (focalLengthMin > 0.) {
            s << "mm";
        }

        s << " (" << lensID << ")";

        if (0 == nFound) {
            return s.str();
        }

        double deltaMin = 1000.;

        std::string bestMatch (s.str());
        std::ostringstream candidates;

        for (r = choices.lower_bound (lensID); r != choices.upper_bound (lensID); r++) {
            double a1, a2, f1, f2, dif;

            if ( !extractLensInfo ( r->second, f1, f2, a1, a2) ) {
                continue;
            }

            if ( f1 == 0. || a1 == 0.) {
                continue;
            }

            if ( focalLength < f1 - .5 || focalLength > f2 + 0.5 ) {
                continue;
            }

            if ( focalLengthMin > 0. && fabs (f1 - focalLengthMin) > 0.5 ) {
                continue;
            }

            if ( focalLengthMax > 0. && fabs (f2 - focalLengthMax) > 0.5 ) {
                continue;
            }

            if ( maxApertureAtFocal > 0.1) {
                double lensAperture;

                if ( maxApertureAtFocal < a1 - 0.15 || maxApertureAtFocal > a2 + 0.15) {
                    continue;
                }

                if ( a1 == a2 || f1 == f2) {
                    lensAperture = a1;
                } else {
                    lensAperture = exp ( log (a1) + (log (a2) - log (a1)) / (log (f2) - log (f1)) * (log (focalLength) - log (f1)) );
                }

                dif = abs (lensAperture - maxApertureAtFocal);
            } else {
                dif = 0;
            }

            if ( dif < deltaMin ) {
                deltaMin = dif;
                bestMatch = r->second;
            }

            if ( dif < 0.15) {
                if ( candidates.tellp() ) {
                    candidates << "\n or " <<  r->second;
                } else {
                    candidates <<  r->second;
                }
            }

        }

        if ( !candidates.tellp() ) {
            return bestMatch;
        } else {
            return candidates.str();
        }
    }
};
CALensInterpreter caLensInterpreter;

class CAFocalTypeInterpreter : public ChoiceInterpreter<>
{
public:
    CAFocalTypeInterpreter()
    {
        choices[0] = "Fixed";
        choices[1] = "Fixed";
        choices[2] = "Zoom";
    }
};
CAFocalTypeInterpreter caFocalTypeInterpreter;

class CAFocalPlaneInterpreter : public Interpreter
{
public:
    std::string toString (const Tag* t) const override
    {
        int val = t->toInt();

        if ( val < 40 ) {
            return "undef";
        }

        char buffer[32];
        sprintf (buffer, "%.2fmm", val * 25.4 / 1000);
        return buffer;
    }
};
CAFocalPlaneInterpreter caFocalPlaneInterpreter;

class CAExposureTimeInterpreter : public Interpreter
{
public:
    std::string toString (const Tag* t) const override
    {
        char buffer[32];
        double d = pow (2, - t->toInt() / 32.0);
        sprintf (buffer, "%.3f", d);
        return buffer;
    }
};
CAExposureTimeInterpreter caExposureTimeInterpreter;

class CAEVInterpreter : public Interpreter
{
    std::string toString (const Tag* t) const override
    {
        char buffer[32];
        sprintf (buffer, "%.1f", t->toDouble() / 32.0  );
        return buffer;
    }
};
CAEVInterpreter caEVInterpreter;

class CABaseISOInterpreter : public Interpreter
{
public:
    std::string toString (const Tag* t) const override
    {
        char buffer[32];
        int a = t->toInt();
        sprintf (buffer, "%d", a);
        return buffer;
    }
    double toDouble (const Tag* t, int ofs) override
    {
        int a = Interpreter::toInt (t, ofs);

        if (a > 1) {
            double i = pow (2., double (a) / 32. - 4.) * 50.;
            return i;
        } else {
            return 0.;
        }
    }
    int toInt (const Tag* t, int ofs, TagType astype) override
    {
        int a = Interpreter::toInt (t, ofs, astype);

        if (a > 1) {
            int i = static_cast<double>(powf (2.f, static_cast<float>(a) / 32.f - 4.f)) * 50.0 + 0.5;
            return i;
        } else {
            return 0;
        }
    }
};
CABaseISOInterpreter caBaseISOInterpreter;

class CAToneCurveInterpreter : public ChoiceInterpreter<>
{
public:
    CAToneCurveInterpreter()
    {
        choices[0] = "Standard";
        choices[1] = "Manual";
        choices[2] = "Custom";
    }
};
CAToneCurveInterpreter caToneCurveInterpreter;

class CASharpnessFrequencyInterpreter : public ChoiceInterpreter<>
{
public:
    CASharpnessFrequencyInterpreter()
    {
        choices[0] = "N/A";
        choices[1] = "Lowest";
        choices[2] = "Low";
        choices[3] = "Standard";
        choices[4] = "High";
        choices[5] = "Highest";
    }
};
CASharpnessFrequencyInterpreter caSharpnessFrequencyInterpreter;

class CAWhiteBalanceInterpreter : public ChoiceInterpreter<>
{
public:
    CAWhiteBalanceInterpreter()
    {
        choices[0] = "Auto";
        choices[1] = "Daylight";
        choices[2] = "Cloudy";
        choices[3] = "Tungsten";
        choices[4] = "Fluorescent";
        choices[5] = "Flash";
        choices[6] = "Custom";
        choices[7] = "Black & White";
        choices[8] = "Shade";
        choices[9] = "Manual Temperature (Kelvin)";
        choices[10] = "PC Set1";
        choices[11] = "PC Set2";
        choices[12] = "PC Set3";
        choices[14] = "Daylight Fluorescent";
        choices[15] = "Custom 1";
        choices[16] = "Custom 2";
        choices[17] = "Underwater";
        choices[18] = "Custom 3";
        choices[19] = "Custom 4";
        choices[20] = "PC Set4";
        choices[21] = "PC Set5";
        choices[23] = "Auto (ambience priority)";
    }
};
CAWhiteBalanceInterpreter caWhiteBalanceInterpreter;

class CAPictureStyleInterpreter : public ChoiceInterpreter<>
{
public:
    CAPictureStyleInterpreter()
    {
        choices[0] = "None";
        choices[1] = "Standard";
        choices[2] = "Portrait";
        choices[3] = "High Saturation";
        choices[4] = "Adobe RGB";
        choices[5] = "Low Saturation";
        choices[6] = "CM Set 1";
        choices[7] = "CM Set 2";
        choices[0x21] = "User Def. 1";
        choices[0x22] = "User Def. 2";
        choices[0x23] = "User Def. 3";
        choices[0x41] = "PC 1";
        choices[0x42] = "PC 2";
        choices[0x43] = "PC 3";
        choices[0x81] = "Standard";
        choices[0x82] = "Portrait";
        choices[0x83] = "Landscape";
        choices[0x84] = "Neutral";
        choices[0x85] = "Faithful";
        choices[0x86] = "Monochrome";
        choices[0x87] = "Auto";
        choices[0x88] = "Fine Detail";
    }
};
CAPictureStyleInterpreter caPictureStyleInterpreter;

class CASlowShutterInterpreter : public ChoiceInterpreter<>
{
public:
    CASlowShutterInterpreter()
    {
        choices[0] = "Off";
        choices[1] = "Night Scene";
        choices[2] = "On";
        choices[3] = "None";
    }
};
CASlowShutterInterpreter caSlowShutterInterpreter;

class CAFlashGuideNumberInterpreter : public Interpreter
{
public:
    std::string toString (const Tag* t) const override
    {
        int n = t->toInt();

        if ( n == -1) {
            return "undef";
        }

        char buffer[32];
        sprintf (buffer, "%.0f", n / 32. );
        return buffer;
    }
};
CAFlashGuideNumberInterpreter caFlashGuideNumberInterpreter;

class CAAFPointsInFocusInterpreter : public ChoiceInterpreter<>
{
public:
    CAAFPointsInFocusInterpreter()
    {
        choices[0x3000] = "None (MF)";
        choices[0x3001] = "Right";
        choices[0x3002] = "Center";
        choices[0x3003] = "Center+Right";
        choices[0x3004] = "Left";
        choices[0x3005] = "Left+Right";
        choices[0x3006] = "Left+Center";
        choices[0x3007] = "All";
    }
};
CAAFPointsInFocusInterpreter caAFPointsInFocusInterpreter;

class CAAutoExposureBracketingInterpreter : public ChoiceInterpreter<int>
{
public:
    CAAutoExposureBracketingInterpreter()
    {
        choices[-1] = "On ";
        choices[0] = "Off ";
        choices[1] = "On (shot 1)";
        choices[2] = "On (shot 2)";
        choices[3] = "On (shot 3)";
    }
};
CAAutoExposureBracketingInterpreter caAutoExposureBracketingInterpreter;

class CAControModeInterpreter : public ChoiceInterpreter<>
{
public:
    CAControModeInterpreter()
    {
        choices[0] = "n/a";
        choices[1] = "Camera Local Control";
        choices[3] = "Computer Remote Control";
    }
};
CAControModeInterpreter caControModeInterpreter;

class CAFocusDistanceInterpreter : public Interpreter
{
public:
    std::string toString (const Tag* t) const override
    {
        char buffer[32];
        sprintf (buffer, "%.2f", t->toDouble() / 100 );
        return buffer;
    }
};
CAFocusDistanceInterpreter caFocusDistanceInterpreter;

class CAMeasuredEVInterpreter : public Interpreter
{
public:
    std::string toString (const Tag* t) const override
    {
        char buffer[32];
        sprintf (buffer, "%.1f", t->toDouble() / 8 - 6 );
        return buffer;
    }
};
CAMeasuredEVInterpreter caMeasuredEVInterpreter;

class CACameraTypeInterpreter : public ChoiceInterpreter<>
{
public:
    CACameraTypeInterpreter()
    {
        choices[248] = "EOS High-end";
        choices[250] = "Compact";
        choices[252] = "EOS Mid-range";
        choices[255] = "DV Camera";
    }
};
CACameraTypeInterpreter caCameraTypeInterpreter;

class CAAutoRotateInterpreter : public ChoiceInterpreter<int>
{
public:
    CAAutoRotateInterpreter()
    {
        choices[-1] = "Rotated by Software";
        choices[0] = "None";
        choices[1] = "Rotate 90 CW";
        choices[2] = "Rotate 180";
        choices[3] = "Rotate 270 CW";
    }
};
CAAutoRotateInterpreter caAutoRotateInterpreter;

class CABracketModeInterpreter : public ChoiceInterpreter<>
{
public:
    CABracketModeInterpreter()
    {
        choices[0] = "Off";
        choices[1] = "AEB";
        choices[2] = "FEB";
        choices[3] = "ISO";
        choices[4] = "WB";
    }
};
CABracketModeInterpreter caBracketModeInterpreter;

class CARAWJpegQualityInterpreter : public ChoiceInterpreter<>
{
public:
    CARAWJpegQualityInterpreter()
    {
        choices[1] = "Economy";
        choices[2] = "Normal";
        choices[3] = "Fine";
        choices[4] = "RAW";
        choices[5] = "Superfine";
        choices[130] = "Normal Movie";
        choices[131] = "Movie (2)";
    }
};
CARAWJpegQualityInterpreter caRAWJpegQualityInterpreter;

class CAJpegSizeInterpreter : public ChoiceInterpreter<>
{
public:
    CAJpegSizeInterpreter()
    {
        choices[0] = "Large";
        choices[1] = "Medium";
        choices[2] = "Small";
        choices[5] = "Medium 1";
        choices[6] = "Medium 2";
        choices[7] = "Medium 3";
        choices[8] = "Postcard";
        choices[9] = "Widescreen";
        choices[10] = "Medium Widescreen";
        choices[14] = "Small 1";
        choices[15] = "Small 2";
        choices[16] = "Small 3";
        choices[128] = "640x480 Movie";
        choices[129] = "Medium Movie";
        choices[130] = "Small Movie";
        choices[137] = "1280x720 Movie";
        choices[142] = "1920x1080 Movie";
    }
};
CAJpegSizeInterpreter caJpegSizeInterpreter;

class CAWBBracketModeInterpreter : public ChoiceInterpreter<>
{
public:
    CAWBBracketModeInterpreter()
    {
        choices[0] = "Off";
        choices[1] = "On (shift AB)";
        choices[2] = "On (shift GM)";
    }
};
CAWBBracketModeInterpreter caWBBracketModeInterpreter;

class CAFilterEffectInterpreter : public ChoiceInterpreter<>
{
public:
    CAFilterEffectInterpreter()
    {
        choices[0] = "None";
        choices[1] = "Yellow";
        choices[2] = "Orange";
        choices[3] = "Red";
        choices[4] = "Green";
    }
};
CAFilterEffectInterpreter caFilterEffectInterpreter;

class CAToningEffectInterpreter : public ChoiceInterpreter<>
{
public:
    CAToningEffectInterpreter()
    {
        choices[0] = "None";
        choices[1] = "Sepia";
        choices[2] = "Blue";
        choices[3] = "Purple";
        choices[4] = "Green";
    }
};
CAToningEffectInterpreter caToningEffectInterpreter;

class CAFileNumberInterpreter : public Interpreter
{
public:
    std::string toString (const Tag* t) const override
    {
        unsigned long val = t->toInt (0, LONG);
        char buffer[32];
        sprintf (buffer, "%ld", ((val & 0xffc0) >> 6) * 10000 + ((val >> 16) & 0xff) + ((val & 0x3f) << 8) );
        return buffer;
    }
};
CAFileNumberInterpreter caFileNumberInterpreter;

// CanonModelID
class CAModelIDInterpreter : public ChoiceInterpreter<>
{
public:
    CAModelIDInterpreter ()
    {
        choices[1042] = "EOS M50 / Kiss M";
        choices[2049] = "PowerShot SX740 HS";
        choices[2052] = "PowerShot G5 X Mark II";
        choices[2053] = "PowerShot SX70 HS";
        choices[2056] = "PowerShot G7 X Mark III";
        choices[2065] = "EOS M6 Mark II";
        choices[2066] = "EOS M200";
        choices[16842752] = "PowerShot A30";
        choices[17039360] = "PowerShot S300 / Digital IXUS 300 / IXY Digital 300";
        choices[17170432] = "PowerShot A20";
        choices[17301504] = "PowerShot A10";
        choices[17367040] = "PowerShot S110 / Digital IXUS v / IXY Digital 200";
        choices[17825792] = "PowerShot G2";
        choices[17891328] = "PowerShot S40";
        choices[17956864] = "PowerShot S30";
        choices[18022400] = "PowerShot A40";
        choices[18087936] = "EOS D30";
        choices[18153472] = "PowerShot A100";
        choices[18219008] = "PowerShot S200 / Digital IXUS v2 / IXY Digital 200a";
        choices[18284544] = "PowerShot A200";
        choices[18350080] = "PowerShot S330 / Digital IXUS 330 / IXY Digital 300a";
        choices[18415616] = "PowerShot G3";
        choices[18939904] = "PowerShot S45";
        choices[19070976] = "PowerShot SD100 / Digital IXUS II / IXY Digital 30";
        choices[19136512] = "PowerShot S230 / Digital IXUS v3 / IXY Digital 320";
        choices[19202048] = "PowerShot A70";
        choices[19267584] = "PowerShot A60";
        choices[19333120] = "PowerShot S400 / Digital IXUS 400 / IXY Digital 400";
        choices[19464192] = "PowerShot G5";
        choices[19922944] = "PowerShot A300";
        choices[19988480] = "PowerShot S50";
        choices[20185088] = "PowerShot A80";
        choices[20250624] = "PowerShot SD10 / Digital IXUS i / IXY Digital L";
        choices[20316160] = "PowerShot S1 IS";
        choices[20381696] = "PowerShot Pro1";
        choices[20447232] = "PowerShot S70";
        choices[20512768] = "PowerShot S60";
        choices[20971520] = "PowerShot G6";
        choices[21037056] = "PowerShot S500 / Digital IXUS 500 / IXY Digital 500";
        choices[21102592] = "PowerShot A75";
        choices[21233664] = "PowerShot SD110 / Digital IXUS IIs / IXY Digital 30a";
        choices[21299200] = "PowerShot A400";
        choices[21430272] = "PowerShot A310";
        choices[21561344] = "PowerShot A85";
        choices[22151168] = "PowerShot S410 / Digital IXUS 430 / IXY Digital 450";
        choices[22216704] = "PowerShot A95";
        choices[22282240] = "PowerShot SD300 / Digital IXUS 40 / IXY Digital 50";
        choices[22347776] = "PowerShot SD200 / Digital IXUS 30 / IXY Digital 40";
        choices[22413312] = "PowerShot A520";
        choices[22478848] = "PowerShot A510";
        choices[22609920] = "PowerShot SD20 / Digital IXUS i5 / IXY Digital L2";
        choices[23330816] = "PowerShot S2 IS";
        choices[23396352] = "PowerShot SD430 / Digital IXUS Wireless / IXY Digital Wireless";
        choices[23461888] = "PowerShot SD500 / Digital IXUS 700 / IXY Digital 600";
        choices[23494656] = "EOS D60";
        choices[24117248] = "PowerShot SD30 / Digital IXUS i Zoom / IXY Digital L3";
        choices[24379392] = "PowerShot A430";
        choices[24444928] = "PowerShot A410";
        choices[24510464] = "PowerShot S80";
        choices[24641536] = "PowerShot A620";
        choices[24707072] = "PowerShot A610";
        choices[25165824] = "PowerShot SD630 / Digital IXUS 65 / IXY Digital 80";
        choices[25231360] = "PowerShot SD450 / Digital IXUS 55 / IXY Digital 60";
        choices[25296896] = "PowerShot TX1";
        choices[25624576] = "PowerShot SD400 / Digital IXUS 50 / IXY Digital 55";
        choices[25690112] = "PowerShot A420";
        choices[25755648] = "PowerShot SD900 / Digital IXUS 900 Ti / IXY Digital 1000";
        choices[26214400] = "PowerShot SD550 / Digital IXUS 750 / IXY Digital 700";
        choices[26345472] = "PowerShot A700";
        choices[26476544] = "PowerShot SD700 IS / Digital IXUS 800 IS / IXY Digital 800 IS";
        choices[26542080] = "PowerShot S3 IS";
        choices[26607616] = "PowerShot A540";
        choices[26673152] = "PowerShot SD600 / Digital IXUS 60 / IXY Digital 70";
        choices[26738688] = "PowerShot G7";
        choices[26804224] = "PowerShot A530";
        choices[33554432] = "PowerShot SD800 IS / Digital IXUS 850 IS / IXY Digital 900 IS";
        choices[33619968] = "PowerShot SD40 / Digital IXUS i7 / IXY Digital L4";
        choices[33685504] = "PowerShot A710 IS";
        choices[33751040] = "PowerShot A640";
        choices[33816576] = "PowerShot A630";
        choices[34144256] = "PowerShot S5 IS";
        choices[34603008] = "PowerShot A460";
        choices[34734080] = "PowerShot SD850 IS / Digital IXUS 950 IS / IXY Digital 810 IS";
        choices[34799616] = "PowerShot A570 IS";
        choices[34865152] = "PowerShot A560";
        choices[34930688] = "PowerShot SD750 / Digital IXUS 75 / IXY Digital 90";
        choices[34996224] = "PowerShot SD1000 / Digital IXUS 70 / IXY Digital 10";
        choices[35127296] = "PowerShot A550";
        choices[35192832] = "PowerShot A450";
        choices[35848192] = "PowerShot G9";
        choices[35913728] = "PowerShot A650 IS";
        choices[36044800] = "PowerShot A720 IS";
        choices[36241408] = "PowerShot SX100 IS";
        choices[36700160] = "PowerShot SD950 IS / Digital IXUS 960 IS / IXY Digital 2000 IS";
        choices[36765696] = "PowerShot SD870 IS / Digital IXUS 860 IS / IXY Digital 910 IS";
        choices[36831232] = "PowerShot SD890 IS / Digital IXUS 970 IS / IXY Digital 820 IS";
        choices[37093376] = "PowerShot SD790 IS / Digital IXUS 90 IS / IXY Digital 95 IS";
        choices[37158912] = "PowerShot SD770 IS / Digital IXUS 85 IS / IXY Digital 25 IS";
        choices[37224448] = "PowerShot A590 IS";
        choices[37289984] = "PowerShot A580";
        choices[37879808] = "PowerShot A470";
        choices[37945344] = "PowerShot SD1100 IS / Digital IXUS 80 IS / IXY Digital 20 IS";
        choices[38141952] = "PowerShot SX1 IS";
        choices[38207488] = "PowerShot SX10 IS";
        choices[38273024] = "PowerShot A1000 IS";
        choices[38338560] = "PowerShot G10";
        choices[38862848] = "PowerShot A2000 IS";
        choices[38928384] = "PowerShot SX110 IS";
        choices[38993920] = "PowerShot SD990 IS / Digital IXUS 980 IS / IXY Digital 3000 IS";
        choices[39059456] = "PowerShot SD880 IS / Digital IXUS 870 IS / IXY Digital 920 IS";
        choices[39124992] = "PowerShot E1";
        choices[39190528] = "PowerShot D10";
        choices[39256064] = "PowerShot SD960 IS / Digital IXUS 110 IS / IXY Digital 510 IS";
        choices[39321600] = "PowerShot A2100 IS";
        choices[39387136] = "PowerShot A480";
        choices[39845888] = "PowerShot SX200 IS";
        choices[39911424] = "PowerShot SD970 IS / Digital IXUS 990 IS / IXY Digital 830 IS";
        choices[39976960] = "PowerShot SD780 IS / Digital IXUS 100 IS / IXY Digital 210 IS";
        choices[40042496] = "PowerShot A1100 IS";
        choices[40108032] = "PowerShot SD1200 IS / Digital IXUS 95 IS / IXY Digital 110 IS";
        choices[40894464] = "PowerShot G11";
        choices[40960000] = "PowerShot SX120 IS";
        choices[41025536] = "PowerShot S90";
        choices[41222144] = "PowerShot SX20 IS";
        choices[41287680] = "PowerShot SD980 IS / Digital IXUS 200 IS / IXY Digital 930 IS";
        choices[41353216] = "PowerShot SD940 IS / Digital IXUS 120 IS / IXY Digital 220 IS";
        choices[41943040] = "PowerShot A495";
        choices[42008576] = "PowerShot A490";
        choices[42074112] = "PowerShot A3100/A3150 IS";
        choices[42139648] = "PowerShot A3000 IS";
        choices[42205184] = "PowerShot SD1400 IS / IXUS 130 / IXY 400F";
        choices[42270720] = "PowerShot SD1300 IS / IXUS 105 / IXY 200F";
        choices[42336256] = "PowerShot SD3500 IS / IXUS 210 / IXY 10S";
        choices[42401792] = "PowerShot SX210 IS";
        choices[42467328] = "PowerShot SD4000 IS / IXUS 300 HS / IXY 30S";
        choices[42532864] = "PowerShot SD4500 IS / IXUS 1000 HS / IXY 50S";
        choices[43122688] = "PowerShot G12";
        choices[43188224] = "PowerShot SX30 IS";
        choices[43253760] = "PowerShot SX130 IS";
        choices[43319296] = "PowerShot S95";
        choices[43515904] = "PowerShot A3300 IS";
        choices[43581440] = "PowerShot A3200 IS";
        choices[50331648] = "PowerShot ELPH 500 HS / IXUS 310 HS / IXY 31S";
        choices[50397184] = "PowerShot Pro90 IS";
        choices[50397185] = "PowerShot A800";
        choices[50462720] = "PowerShot ELPH 100 HS / IXUS 115 HS / IXY 210F";
        choices[50528256] = "PowerShot SX230 HS";
        choices[50593792] = "PowerShot ELPH 300 HS / IXUS 220 HS / IXY 410F";
        choices[50659328] = "PowerShot A2200";
        choices[50724864] = "PowerShot A1200";
        choices[50790400] = "PowerShot SX220 HS";
        choices[50855936] = "PowerShot G1 X";
        choices[50921472] = "PowerShot SX150 IS";
        choices[51380224] = "PowerShot ELPH 510 HS / IXUS 1100 HS / IXY 51S";
        choices[51445760] = "PowerShot S100 (new)";
        choices[51511296] = "PowerShot ELPH 310 HS / IXUS 230 HS / IXY 600F";
        choices[51576832] = "PowerShot SX40 HS";
        choices[51642368] = "IXY 32S";
        choices[51773440] = "PowerShot A1300";
        choices[51838976] = "PowerShot A810";
        choices[51904512] = "PowerShot ELPH 320 HS / IXUS 240 HS / IXY 420F";
        choices[51970048] = "PowerShot ELPH 110 HS / IXUS 125 HS / IXY 220F";
        choices[52428800] = "PowerShot D20";
        choices[52494336] = "PowerShot A4000 IS";
        choices[52559872] = "PowerShot SX260 HS";
        choices[52625408] = "PowerShot SX240 HS";
        choices[52690944] = "PowerShot ELPH 530 HS / IXUS 510 HS / IXY 1";
        choices[52756480] = "PowerShot ELPH 520 HS / IXUS 500 HS / IXY 3";
        choices[52822016] = "PowerShot A3400 IS";
        choices[52887552] = "PowerShot A2400 IS";
        choices[52953088] = "PowerShot A2300";
        choices[53608448] = "PowerShot S100V";
        choices[53673984] = "PowerShot G15";
        choices[53739520] = "PowerShot SX50 HS";
        choices[53805056] = "PowerShot SX160 IS";
        choices[53870592] = "PowerShot S110 (new)";
        choices[53936128] = "PowerShot SX500 IS";
        choices[54001664] = "PowerShot N";
        choices[54067200] = "IXUS 245 HS / IXY 430F";
        choices[54525952] = "PowerShot SX280 HS";
        choices[54591488] = "PowerShot SX270 HS";
        choices[54657024] = "PowerShot A3500 IS";
        choices[54722560] = "PowerShot A2600";
        choices[54788096] = "PowerShot SX275 HS";
        choices[54853632] = "PowerShot A1400";
        choices[54919168] = "PowerShot ELPH 130 IS / IXUS 140 / IXY 110F";
        choices[54984704] = "PowerShot ELPH 115/120 IS / IXUS 132/135 / IXY 90F/100F";
        choices[55115776] = "PowerShot ELPH 330 HS / IXUS 255 HS / IXY 610F";
        choices[55640064] = "PowerShot A2500";
        choices[55836672] = "PowerShot G16";
        choices[55902208] = "PowerShot S120";
        choices[55967744] = "PowerShot SX170 IS";
        choices[56098816] = "PowerShot SX510 HS";
        choices[56164352] = "PowerShot S200 (new)";
        choices[56623104] = "IXY 620F";
        choices[56688640] = "PowerShot N100";
        choices[56885248] = "PowerShot G1 X Mark II";
        choices[56950784] = "PowerShot D30";
        choices[57016320] = "PowerShot SX700 HS";
        choices[57081856] = "PowerShot SX600 HS";
        choices[57147392] = "PowerShot ELPH 140 IS / IXUS 150 / IXY 130";
        choices[57212928] = "PowerShot ELPH 135 / IXUS 145 / IXY 120";
        choices[57671680] = "PowerShot ELPH 340 HS / IXUS 265 HS / IXY 630";
        choices[57737216] = "PowerShot ELPH 150 IS / IXUS 155 / IXY 140";
        choices[57933824] = "EOS M3";
        choices[57999360] = "PowerShot SX60 HS";
        choices[58064896] = "PowerShot SX520 HS";
        choices[58130432] = "PowerShot SX400 IS";
        choices[58195968] = "PowerShot G7 X";
        choices[58261504] = "PowerShot N2";
        choices[58720256] = "PowerShot SX530 HS";
        choices[58851328] = "PowerShot SX710 HS";
        choices[58916864] = "PowerShot SX610 HS";
        choices[58982400] = "EOS M10";
        choices[59047936] = "PowerShot G3 X";
        choices[59113472] = "PowerShot ELPH 165 HS / IXUS 165 / IXY 160";
        choices[59179008] = "PowerShot ELPH 160 / IXUS 160";
        choices[59244544] = "PowerShot ELPH 350 HS / IXUS 275 HS / IXY 640";
        choices[59310080] = "PowerShot ELPH 170 IS / IXUS 170";
        choices[59834368] = "PowerShot SX410 IS";
        choices[59965440] = "PowerShot G9 X";
        choices[60030976] = "EOS M5";
        choices[60096512] = "PowerShot G5 X";
        choices[60227584] = "PowerShot G7 X Mark II";
        choices[60293120] = "EOS M100";
        choices[60358656] = "PowerShot ELPH 360 HS / IXUS 285 HS / IXY 650";
        choices[67174400] = "PowerShot SX540 HS";
        choices[67239936] = "PowerShot SX420 IS";
        choices[67305472] = "PowerShot ELPH 190 IS / IXUS 180 / IXY 190";
        choices[67371008] = "PowerShot G1";
        choices[67371009] = "PowerShot ELPH 180 IS / IXUS 175 / IXY 180";
        choices[67436544] = "PowerShot SX720 HS";
        choices[67502080] = "PowerShot SX620 HS";
        choices[67567616] = "EOS M6";
        choices[68157440] = "PowerShot G9 X Mark II";
        choices[68485120] = "PowerShot ELPH 185 / IXUS 185 / IXY 200";
        choices[68550656] = "PowerShot SX430 IS";
        choices[68616192] = "PowerShot SX730 HS";
        choices[68681728] = "PowerShot G1 X Mark III";
        choices[100925440] = "PowerShot S100 / Digital IXUS / IXY Digital";
        choices[1074255475] = "DC19/DC21/DC22";
        choices[1074255476] = "XH A1";
        choices[1074255477] = "HV10";
        choices[1074255478] = "MD130/MD140/MD150/MD160/ZR850";
        choices[1074255735] = "DC50";
        choices[1074255736] = "HV20";
        choices[1074255737] = "DC211";
        choices[1074255738] = "HG10";
        choices[1074255739] = "HR10";
        choices[1074255741] = "MD255/ZR950";
        choices[1074255900] = "HF11";
        choices[1074255992] = "HV30";
        choices[1074255996] = "XH A1S";
        choices[1074255998] = "DC301/DC310/DC311/DC320/DC330";
        choices[1074255999] = "FS100";
        choices[1074256000] = "HF10";
        choices[1074256002] = "HG20/HG21";
        choices[1074256165] = "HF21";
        choices[1074256166] = "HF S11";
        choices[1074256248] = "HV40";
        choices[1074256263] = "DC410/DC411/DC420";
        choices[1074256264] = "FS19/FS20/FS21/FS22/FS200";
        choices[1074256265] = "HF20/HF200";
        choices[1074256266] = "HF S10/S100";
        choices[1074256526] = "HF R10/R16/R17/R18/R100/R106";
        choices[1074256527] = "HF M30/M31/M36/M300/M306";
        choices[1074256528] = "HF S20/S21/S200";
        choices[1074256530] = "FS31/FS36/FS37/FS300/FS305/FS306/FS307";
        choices[1074257056] = "EOS C300";
        choices[1074257321] = "HF G25";
        choices[1074257844] = "XC10";
        choices[1074258371] = "EOS C200";
        choices[2147483649] = "EOS-1D";
        choices[2147484007] = "EOS-1DS";
        choices[2147484008] = "EOS 10D";
        choices[2147484009] = "EOS-1D Mark III";
        choices[2147484016] = "EOS Digital Rebel / 300D / Kiss Digital";
        choices[2147484020] = "EOS-1D Mark II";
        choices[2147484021] = "EOS 20D";
        choices[2147484022] = "EOS Digital Rebel XSi / 450D / Kiss X2";
        choices[2147484040] = "EOS-1Ds Mark II";
        choices[2147484041] = "EOS Digital Rebel XT / 350D / Kiss Digital N";
        choices[2147484048] = "EOS 40D";
        choices[2147484179] = "EOS 5D";
        choices[2147484181] = "EOS-1Ds Mark III";
        choices[2147484184] = "EOS 5D Mark II";
        choices[2147484185] = "WFT-E1";
        choices[2147484210] = "EOS-1D Mark II N";
        choices[2147484212] = "EOS 30D";
        choices[2147484214] = "EOS Digital Rebel XTi / 400D / Kiss Digital X";
        choices[2147484225] = "WFT-E2";
        choices[2147484230] = "WFT-E3";
        choices[2147484240] = "EOS 7D";
        choices[2147484242] = "EOS Rebel T1i / 500D / Kiss X3";
        choices[2147484244] = "EOS Rebel XS / 1000D / Kiss F";
        choices[2147484257] = "EOS 50D";
        choices[2147484265] = "EOS-1D X";
        choices[2147484272] = "EOS Rebel T2i / 550D / Kiss X4";
        choices[2147484273] = "WFT-E4";
        choices[2147484275] = "WFT-E5";
        choices[2147484289] = "EOS-1D Mark IV";
        choices[2147484293] = "EOS 5D Mark III";
        choices[2147484294] = "EOS Rebel T3i / 600D / Kiss X5";
        choices[2147484295] = "EOS 60D";
        choices[2147484296] = "EOS Rebel T3 / 1100D / Kiss X50";
        choices[2147484297] = "EOS 7D Mark II";
        choices[2147484311] = "WFT-E2 II";
        choices[2147484312] = "WFT-E4 II";
        choices[2147484417] = "EOS Rebel T4i / 650D / Kiss X6i";
        choices[2147484418] = "EOS 6D";
        choices[2147484452] = "EOS-1D C";
        choices[2147484453] = "EOS 70D";
        choices[2147484454] = "EOS Rebel T5i / 700D / Kiss X7i";
        choices[2147484455] = "EOS Rebel T5 / 1200D / Kiss X70 / Hi";
        choices[2147484456] = "EOS-1D X Mark II";
        choices[2147484465] = "EOS M";
        choices[2147484486] = "EOS Rebel SL1 / 100D / Kiss X7";
        choices[2147484487] = "EOS Rebel T6s / 760D / 8000D";
        choices[2147484489] = "EOS 5D Mark IV";
        choices[2147484496] = "EOS 80D";
        choices[2147484501] = "EOS M2";
        choices[2147484546] = "EOS 5DS";
        choices[2147484563] = "EOS Rebel T6i / 750D / Kiss X8i";
        choices[2147484673] = "EOS 5DS R";
        choices[2147484676] = "EOS Rebel T6 / 1300D / Kiss X80";
        choices[2147484677] = "EOS Rebel T7i / 800D / Kiss X9i";
        choices[2147484678] = "EOS 6D Mark II";
        choices[2147484680] = "EOS 77D / 9000D";
        choices[2147484695] = "EOS Rebel SL2 / 200D / Kiss X9";
        choices[2147484706] = "EOS Rebel T100 / 4000D / 3000D";
        choices[2147484708] = "EOS R";
        choices[2147484712] = "EOS-1D X Mark III";
        choices[2147484722] = "EOS Rebel T7 / 2000D / 1500D / Kiss X90";
        choices[2147484723] = "EOS RP";
        choices[2147484725] = "EOS Rebel T8i / 850D / X10i";
        choices[2147484726] = "EOS SL3 / 250D / Kiss X10";
        choices[2147484727] = "EOS 90D";
        choices[2147484960] = "EOS D2000C";
        choices[2147485024] = "EOS D6000C";
    }
};
CAModelIDInterpreter caModelIDInterpreter;

class CAPanoramaDirectionInterpreter : public ChoiceInterpreter<>
{
public:
    CAPanoramaDirectionInterpreter()
    {
        choices[0] = "Left to Right";
        choices[1] = "Right to Left";
        choices[2] = "Bottom to Top";
        choices[3] = "Top to Bottom";
        choices[4] = "2x2 Matrix (Clockwise)";
    }
};
CAPanoramaDirectionInterpreter caPanoramaDirectionInterpreter;

class CAAspectRatioInterpreter : public ChoiceInterpreter<>
{
public:
    CAAspectRatioInterpreter()
    {
        choices[0] = "3:2";
        choices[1] = "1:1";
        choices[2] = "4:3";
        choices[7] = "16:9";
        choices[8] = "4:5";
    }

};
CAAspectRatioInterpreter caAspectRatioInterpreter;

const TagAttrib canonCameraSettingsAttribs[] = {
    {0, AC_WRITE, 0, nullptr,  1, AUTO, "MacroMode", &caMacroModeInterpreter},
    {0, AC_WRITE, 0, nullptr,  2, AUTO, "SelfTimer", &caSelfTimerInterpreter},
    {0, AC_WRITE, 0, nullptr,  3, AUTO, "Quality", &caQualityInterpreter},
    {0, AC_WRITE, 0, nullptr,  4, AUTO, "CanonFlashMode", &caFlashModeInterpreter},
    {0, AC_WRITE, 0, nullptr,  5, AUTO, "ContinuousDrive", &caContinuousDriveInterpreter},
    {0, AC_WRITE, 0, nullptr,  7, AUTO, "FocusMode", &caFocusModeInterpreter},
    {0, AC_WRITE, 0, nullptr,  9, AUTO, "RecordMode", &caRecordModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 10, AUTO, "CanonImageSize", &caImageSizeInterpreter},
    {0, AC_WRITE, 0, nullptr, 11, AUTO, "EasyMode", &caEasyModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 12, AUTO, "DigitalZoom", &caDigitalZoomInterpreter},
    {0, AC_WRITE, 0, nullptr, 13, AUTO, "Contrast", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 14, AUTO, "Saturation", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 15, AUTO, "Sharpness", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 16, AUTO, "CameraISO", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 17, AUTO, "MeteringMode", &caMeteringModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 18, AUTO, "FocusRange", &caFocusRangeInterpreter},
    {0, AC_WRITE, 0, nullptr, 19, AUTO, "AFPoint", &caAFPointInterpreter},
    {0, AC_WRITE, 0, nullptr, 20, AUTO, "CanonExposureMode", &caExposureModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 22, AUTO, "LensID", &caLensInterpreter},
    {0, AC_WRITE, 0, nullptr, 23, AUTO, "LongFocal", &caFocalInterpreter},
    {0, AC_WRITE, 0, nullptr, 24, AUTO, "ShortFocal", &caFocalInterpreter},
    {0, AC_WRITE, 0, nullptr, 25, AUTO, "FocalUnits", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 26, AUTO, "MaxAperture", &caApertureInterpreter},
    {0, AC_WRITE, 0, nullptr, 27, AUTO, "MinAperture", &caApertureInterpreter},
    {0, AC_WRITE, 0, nullptr, 28, AUTO, "FlashActivity", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 29, AUTO, "FlashBits", &caFlashBitsInterpreter},
    {0, AC_WRITE, 0, nullptr, 32, AUTO, "FocusContinuous", &caFocusContinuousInterpreter},
    {0, AC_WRITE, 0, nullptr, 33, AUTO, "AESetting", &caAESettingsInterpreter},
    {0, AC_WRITE, 0, nullptr, 34, AUTO, "ImageStabilization", &caStabilizationInterpreter},
    {0, AC_WRITE, 0, nullptr, 35, AUTO, "DisplayAperture", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 36, AUTO, "ZoomSourceWidth", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 37, AUTO, "ZoomTargetWidth", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 39, AUTO, "SpotMeteringMode", &caSpotMeteringInterpreter},
    {0, AC_WRITE, 0, nullptr, 40, AUTO, "PhotoEffect", &caPhotoEffectInterpreter},
    {0, AC_WRITE, 0, nullptr, 41, AUTO, "ManualFlashOutput", &caManualFlashInterpreter},
    {0, AC_WRITE, 0, nullptr, 42, AUTO, "ColorTone", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 46, AUTO, "SRAWQuality", &caRAWQualityInterpreter},
    { -1, AC_DONTWRITE, 0,  nullptr, 0, AUTO, "", nullptr}
};

const TagAttrib canonFocalLengthAttribs[] = {
    {0, AC_WRITE, 0, nullptr, 0, AUTO, "FocalType", &caFocalTypeInterpreter},
    {0, AC_WRITE, 0, nullptr, 1, AUTO, "FocalLength", &caFocalInterpreter},
    {0, AC_WRITE, 0, nullptr, 2, AUTO, "FocalPlaneXSize", &caFocalPlaneInterpreter},
    {0, AC_WRITE, 0, nullptr, 3, AUTO, "FocalPlaneYSize", &caFocalPlaneInterpreter},
    { -1, AC_DONTWRITE, 0,  nullptr, 0, AUTO, "", nullptr}
};

const TagAttrib canonShotInfoAttribs[] = {
    {0, AC_WRITE, 0, nullptr, 1, AUTO, "AutoISO", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 2, AUTO, "BaseISO", &caBaseISOInterpreter},
    {0, AC_WRITE, 0, nullptr, 3, AUTO, "MeasuredEV", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 4, AUTO, "TargetAperture", &caApertureInterpreter},
    {0, AC_WRITE, 0, nullptr, 5, AUTO, "TargetExposureTime", &caExposureTimeInterpreter},
    {0, AC_WRITE, 0, nullptr, 6, AUTO, "ExposureCompensation", &caEVInterpreter},
    {0, AC_WRITE, 0, nullptr, 7, AUTO, "WhiteBalance", &caWhiteBalanceInterpreter},
    {0, AC_WRITE, 0, nullptr, 8, AUTO, "SlowShutter", &caSlowShutterInterpreter},
    {0, AC_WRITE, 0, nullptr, 9, AUTO, "SequenceNumber", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 10, AUTO, "OpticalZoomCode", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 13, AUTO, "FlashGuideNumber", &caFlashGuideNumberInterpreter},
    {0, AC_WRITE, 0, nullptr, 14, AUTO, "AFPointsInFocus", &caAFPointsInFocusInterpreter},
    {0, AC_WRITE, 0, nullptr, 15, AUTO, "FlashExposureComp", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 16, AUTO, "AutoExposureBracketing", &caAutoExposureBracketingInterpreter},
    {0, AC_WRITE, 0, nullptr, 17, AUTO, "AEBBracketValue", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 18, AUTO, "ControlMode", &caControModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 21, AUTO, "FNumber", &caApertureInterpreter},
    {0, AC_WRITE, 0, nullptr, 22, AUTO, "ExposureTime", &caExposureTimeInterpreter},
    {0, AC_WRITE, 0, nullptr, 23, AUTO, "MeasuredEV2", &caMeasuredEVInterpreter},
    {0, AC_WRITE, 0, nullptr, 24, AUTO, "BulbDuration", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 26, AUTO, "CameraType", &caCameraTypeInterpreter},
    {0, AC_WRITE, 0, nullptr, 27, AUTO, "AutoRotate", &caAutoRotateInterpreter},
    {0, AC_WRITE, 0, nullptr, 28, AUTO, "NDFilter", &caOnOffInterpreter},
    {0, AC_WRITE, 0, nullptr, 29, AUTO, "Self-timer2", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 33, AUTO, "FlashOutput", &stdInterpreter},
    { -1, AC_DONTWRITE, 0, nullptr, 0, AUTO, "", nullptr},
};

const TagAttrib canonFileInfoAttribs[] = {
    {0, AC_WRITE, 0, nullptr, 1, AUTO, "FileNumber",  &caFileNumberInterpreter},
    {0, AC_WRITE, 0, nullptr, 3, AUTO, "BracketMode", &caBracketModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 4, AUTO, "BracketValue", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 5, AUTO, "BracketShotNumber", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 6, AUTO, "RawJpgQuality", &caRAWJpegQualityInterpreter},
    {0, AC_WRITE, 0, nullptr, 7, AUTO, "RawJpgSize", &caJpegSizeInterpreter},
    {0, AC_WRITE, 0, nullptr, 8, AUTO, "NoiseReduction", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 9, AUTO, "WBBracketMode", &caWBBracketModeInterpreter},
    {0, AC_WRITE, 0, nullptr, 12, AUTO, "WBBracketValueAB", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 13, AUTO, "WBBracketValueGM", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 14, AUTO, "FilterEffect", &caFilterEffectInterpreter},
    {0, AC_WRITE, 0, nullptr, 15, AUTO, "ToningEffect", &caToningEffectInterpreter},
    {0, AC_WRITE, 0, nullptr, 19, AUTO, "LiveViewShooting", &caOnOffInterpreter},
    {0, AC_WRITE, 0, nullptr, 20, AUTO, "FocusDistanceUpper", &caFocusDistanceInterpreter},
    {0, AC_WRITE, 0, nullptr, 21, AUTO, "FocusDistanceLower", &caFocusDistanceInterpreter},
    {0, AC_WRITE, 0, nullptr, 25, AUTO, "FlashExposureLock", &caOnOffInterpreter},
    { -1, AC_DONTWRITE, 0, nullptr, 0, AUTO, "", nullptr},
};

const TagAttrib canonProcessingInfoAttribs[] = {
    {0, AC_WRITE, 0, nullptr, 1, AUTO, "ToneCurve", &caToneCurveInterpreter},
    {0, AC_WRITE, 0, nullptr, 2, AUTO, "Sharpness", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 3, AUTO, "SharpnessFrequency", &caSharpnessFrequencyInterpreter},
    {0, AC_WRITE, 0, nullptr, 4, AUTO, "SensorRedLevel", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 5, AUTO, "SensorBlueLevel", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 6, AUTO, "WhiteBalanceRed", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 7, AUTO, "WhiteBalanceBlue", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 8, AUTO, "WhiteBalance", &caWhiteBalanceInterpreter},
    {0, AC_WRITE, 0, nullptr, 9, AUTO, "ColorTemperature", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 10, AUTO, "PictureStyle", &caPictureStyleInterpreter},
    {0, AC_WRITE, 0, nullptr, 11, AUTO, "DigitalGain", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 12, AUTO, "WBShiftAB", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 13, AUTO, "WBShiftGM", &stdInterpreter},
    { -1, AC_DONTWRITE, 0, nullptr, 0, AUTO, "", nullptr},
};

const TagAttrib canonPanoramaInfoAttribs[] = {
    {0, AC_WRITE, 0, nullptr, 2, AUTO, "PanoramaFrameNumber", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 5, AUTO, "PanoramaDirection", &caPanoramaDirectionInterpreter},
    { -1, AC_DONTWRITE, 0, nullptr, 0, AUTO, "", nullptr},
};

const TagAttrib canonCropInfoAttribs[] = {
    {0, AC_WRITE, 0, nullptr, 0, AUTO, "CropLeftMargin", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 1, AUTO, "CropRightMargin", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 2, AUTO, "CropTopMargin", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 3, AUTO, "CropBottomMargin", &stdInterpreter},
    { -1, AC_DONTWRITE, 0, nullptr, 0, AUTO, "", nullptr},
};

const TagAttrib canonAspectInfoAttribs[] = {
    {0, AC_WRITE, 0, nullptr, 0, AUTO, "AspectRatio", &caAspectRatioInterpreter},
    {0, AC_WRITE, 0, nullptr, 1, AUTO, "CroppedImageWidth", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 2, AUTO, "CroppedImageHeight", &stdInterpreter},
    { -1, AC_DONTWRITE, 0, nullptr, 0, AUTO, "", nullptr},
};

const TagAttrib canonMicroAdjustAttrib[] = {
    {0, AC_WRITE, 0, nullptr, 1, AUTO, "AFMicroAdjActive", &caOnOffInterpreter},
    { -1, AC_DONTWRITE, 0, nullptr, 2, AUTO, "", nullptr},
};

const TagAttrib canonAttribs[] = {
    {0, AC_WRITE, 0, canonCameraSettingsAttribs, 0x0001, AUTO, "CanonCameraSettings", &stdInterpreter},
    {0, AC_WRITE, 0, canonFocalLengthAttribs, 0x0002, AUTO, "CanonFocalLength", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0003, AUTO, "CanonFlashInfo", &stdInterpreter},
    {0, AC_WRITE, 0, canonShotInfoAttribs, 0x0004, AUTO, "CanonShotInfo", &stdInterpreter},
    {0, AC_WRITE, 0, canonPanoramaInfoAttribs, 0x0005, AUTO, "CanonPanorama", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0006, AUTO, "CanonImageType", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0007, AUTO, "CanonFirmwareVersion", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0008, AUTO, "FileNumber", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0009, AUTO, "OwnerName", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x000a, AUTO, "ColorInfoD30", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x000c, AUTO, "SerialNumber", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x000d, AUTO, "CanonCameraInfo", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x000e, AUTO, "CanonFileLength", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x000f, AUTO, "CustomFunctions", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0010, AUTO, "CanonModelID", &caModelIDInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0012, AUTO, "CanonAFInfo", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0015, AUTO, "SerialNumberFormat", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x001c, AUTO, "DateStampMode", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x001d, AUTO, "MyColors", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x001e, AUTO, "FirmwareRevision", &stdInterpreter},
    {0, AC_NEW,   0, nullptr, 0x0024, AUTO, "FaceDetect1", &stdInterpreter},
    {0, AC_NEW,   0, nullptr, 0x0025, AUTO, "FaceDetect2", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0026, AUTO, "CanonAFInfo2", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0083, AUTO, "OriginalDecisionData", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0090, AUTO, "CustomFunctions1D", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0091, AUTO, "PersonalFunctions", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0092, AUTO, "PersonalFunctionValues", &stdInterpreter},
    {0, AC_WRITE, 0, canonFileInfoAttribs, 0x0093, AUTO, "CanonFileInfo", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0094, AUTO, "AFPointsInFocus1D", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0095, AUTO, "LensType", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0096, AUTO, "InternalSerialNumber", &caIntSerNumInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0097, AUTO, "DustRemovalData", &stdInterpreter},
    {0, AC_WRITE, 0, canonCropInfoAttribs, 0x0098, AUTO, "CropInfo", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x0099, AUTO, "CustomFunctions2", &stdInterpreter},
    {0, AC_WRITE, 0, canonAspectInfoAttribs, 0x009a, AUTO, "AspectInfo", &stdInterpreter},
    {0, AC_WRITE, 0, canonProcessingInfoAttribs, 0x00a0, AUTO, "ProcessingInfo", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x00a1, AUTO, "ToneCurveTable", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x00a2, AUTO, "SharpnessTable", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x00a3, AUTO, "SharpnessFreqTable", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x00a4, AUTO, "WhiteBalanceTable", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x00a9, AUTO, "ColorBalance", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x00aa, AUTO, "MeasuredColor", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x00ae, AUTO, "ColorTemperature", &stdInterpreter},
    {0, AC_NEW, 0, nullptr, 0x00b0, AUTO, "CanonFlags", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x00b1, AUTO, "ModifiedInfo", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x00b2, AUTO, "ToneCurveMatching", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x00b3, AUTO, "WhiteBalanceMatching", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x00b4, AUTO, "ColorSpace", &stdInterpreter},
    {1, AC_WRITE, 0, nullptr, 0x00b6, AUTO, "PreviewImageInfo", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x00d0, AUTO, "VRDOffset", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x00e0, AUTO, "SensorInfo", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x4001, AUTO, "ColorBalance", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x4002, AUTO, "UnknownBlock1", &stdInterpreter},
    {0, AC_WRITE, 0, nullptr, 0x4003, AUTO, "ColorInfo", &stdInterpreter},
    {1, AC_WRITE, 0, nullptr, 0x4005, AUTO, "UnknownBlock2", &stdInterpreter},
    {1, AC_WRITE, 0, nullptr, 0x4008, AUTO, "BlackLevel", &stdInterpreter},
    {1, AC_WRITE, 0, canonMicroAdjustAttrib, 0x4013, AUTO, "AFMicroAdj", &stdInterpreter},
    { -1, AC_DONTWRITE, 0,  nullptr, 0, AUTO, "", nullptr}
};
}

