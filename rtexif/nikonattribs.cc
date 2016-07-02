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
#ifndef _NIKONATTRIBS_
#define _NIKONATTRIBS_

#include <cstdio>
#include <cstring>
#include <sstream>
#include <iomanip>

#include "rtexif.h"

using namespace std;

namespace rtexif
{

class NAISOInterpreter : public Interpreter
{
public:
    NAISOInterpreter () {}
    virtual std::string toString (Tag* t)
    {
        char buffer[32];
        sprintf (buffer, "%d", t->toInt(2));
        return buffer;
    }
};
NAISOInterpreter naISOInterpreter;

class NAISOInfoISOInterpreter : public Interpreter
{
public:
    NAISOInfoISOInterpreter () {}
    virtual std::string toString (Tag* t)
    {
        char buffer[32];
        int a = t->toInt();
        sprintf (buffer, "%d", a);
        return buffer;
    }
    virtual double toDouble (Tag* t, int ofs)
    {
        int a = t->getValue()[ofs];

        if(a > 1) {
            double i = pow(2., double(a) / 12. - 5.) * 100.;
            return i;
        } else {
            return 0.;
        }
    }
    virtual int toInt (Tag* t, int ofs, TagType astype)
    {
        int a = t->getValue()[ofs];

        if(a > 1) {
            int i = int(double(powf(2.f, float(a) / 12.f - 5.f)) * 100.f + 0.5f);
            return i;
        } else {
            return 0;
        }
    }
};
NAISOInfoISOInterpreter naISOInfoISOInterpreter;

class NAISOExpansionInterpreter : public Interpreter
{
public:
    NAISOExpansionInterpreter () {}
    virtual std::string toString (Tag* t)
    {
        int a = t->toInt();

        // unclear if this interpretation is correct!
        switch (a) {
        case 0x0:
            return "Off";

        case 0x101:
            return "Hi 0.3";

        case 0x102:
            return "Hi 0.5";

        case 0x103:
            return "Hi 0.7";

        case 0x104:
            return "Hi 1.0";

        case 0x105:
            return "Hi 1.3";

        case 0x106:
            return "Hi 1.5";

        case 0x107:
            return "Hi 1.7";

        case 0x108:
            return "Hi 2.0";

        case 0x201:
            return "Lo 0.3";

        case 0x202:
            return "Lo 0.5";

        case 0x203:
            return "Lo 0.7";

        case 0x204:
            return "Lo 1.0";

        default: {
            char buffer[32];
            sprintf(buffer, "0x%04X", a);
            return buffer;
        }
        }
    }
};
NAISOExpansionInterpreter naISOExpansionInterpreter;

class NALensTypeInterpreter : public Interpreter
{
public:
    NALensTypeInterpreter () {}
    virtual std::string toString (Tag* t)
    {
        int a = t->toInt();
        std::ostringstream str;
        str << "MF = " << (a & 1 ? "Yes" : "No") << std::endl;
        str << "D = " << (a & 2 ? "Yes" : "No") << std::endl;
        str << "G = " << (a & 4 ? "Yes" : "No") << std::endl;
        str << "VR = " << (a & 8 ? "Yes" : "No");
        return str.str();
    }
};
NALensTypeInterpreter naLensTypeInterpreter;

class NAFlashModeInterpreter : public ChoiceInterpreter
{
public:
    NAFlashModeInterpreter ()
    {
        choices[0x0] = "Did Not Fire";
        choices[0x1] = "Fired, Manual";
        choices[0x3] = "Not Ready";
        choices[0x7] = "Fired, External";
        choices[0x8] = "Fired, Commander Mode";
        choices[0x9] = "Fired, TTL Mode";
    }
};
NAFlashModeInterpreter naFlashModeInterpreter;

class NAHiISONRInterpreter : public ChoiceInterpreter
{
public:
    // HighISONoiseReduction
    NAHiISONRInterpreter ()
    {
        choices[0x0] = "Off";
        choices[0x1] = "Minimal";
        choices[0x2] = "Low";
        choices[0x3] = "Medium Low";
        choices[0x4] = "Normal";
        choices[0x5] = "Medium High";
        choices[0x6] = "High";
    }
};
NAHiISONRInterpreter naHiISONRInterpreter;

class NAShootingModeInterpreter : public Interpreter
{
public:
    NAShootingModeInterpreter () {}
    virtual std::string toString (Tag* t)
    {
        int a = t->toInt();
        std::ostringstream str;
        str << "Continuous = " << (a & 1 ? "Yes" : "No") << std::endl;
        str << "Delay = " << (a & 2 ? "Yes" : "No") << std::endl;
        str << "PC Control = " << (a & 4 ? "Yes" : "No") << std::endl;
        str << "White-Balance Bracketing = " << (a & 8 ? "Yes" : "No") << std::endl;
        str << "Exposure Bracketing = " << (a & 16 ? "Yes" : "No") << std::endl;
        str << "Auto ISO = " << (a & 32 ? "Yes" : "No") << std::endl;
        str << "IR Control = " << (a & 64 ? "Yes" : "No");
        return str.str();
    }
};
NAShootingModeInterpreter naShootingModeInterpreter;

class NAAFInfoInterpreter : public Interpreter
{
    std::map<int, std::string> amchoices;
    std::map<int, std::string> afpchoices;
public:
    // AFAreaMode
    NAAFInfoInterpreter ()
    {
        amchoices[0x0] = "Single Area";
        amchoices[0x1] = "Dynamic Area";
        amchoices[0x2] = "Dynamic Area (closest subject)";
        amchoices[0x3] = "Group Dynamic";
        amchoices[0x4] = "Single Area (wide)";
        amchoices[0x5] = "Dynamic Area (wide)";
    // AFPoint
        afpchoices[0x0] = "Center";
        afpchoices[0x1] = "Top";
        afpchoices[0x2] = "Bottom";
        afpchoices[0x3] = "Mid-left";
        afpchoices[0x4] = "Mid-right";
        afpchoices[0x5] = "Upper-left";
        afpchoices[0x6] = "Upper-right";
        afpchoices[0x7] = "Lower-left";
        afpchoices[0x8] = "Lower-right";
        afpchoices[0x9] = "Far Left";
        afpchoices[0xa] = "Far Right";
    }
    virtual std::string toString (Tag* t)
    {
        int am = t->toInt (0, BYTE);
        int afp = t->toInt (1, BYTE);
        int aff = t->toInt (2, SHORT);
        std::ostringstream str;
        str << "AFAreaMode = " << amchoices[am] << std::endl;
        str << "AFAreaMode = " << afpchoices[afp] << std::endl;

        std::ostringstream af;

        if (aff & 1)
            if (af.str() == "") {
                af << "Center";
            } else {
                af << ", Center";
            }
        else if (aff & 2)
            if (af.str() == "") {
                af << "Top";
            } else {
                af << ", Top";
            }
        else if (aff & 4)
            if (af.str() == "") {
                af << "Bottom";
            } else {
                af << ", Bottom";
            }
        else if (aff & 8)
            if (af.str() == "") {
                af << "Left";
            } else {
                af << ", Left";
            }
        else if (aff & 16)
            if (af.str() == "") {
                af << "Right";
            } else {
                af << ", Right";
            }
        else if (aff & 32)
            if (af.str() == "") {
                af << "Upper-left";
            } else {
                af << ", Upper-left";
            }
        else if (aff & 64)
            if (af.str() == "") {
                af << "Upper-right";
            } else {
                af << ", Upper-right";
            }
        else if (aff & 128)
            if (af.str() == "") {
                af << " Lower-left";
            } else {
                af << ",  Lower-left";
            }
        else if (aff & 256)
            if (af.str() == "") {
                af << "Lower-right";
            } else {
                af << ", Lower-right";
            }
        else if (aff & 512)
            if (af.str() == "") {
                af << "Far Left";
            } else {
                af << ", Far Left";
            }
        else if (aff & 1024) {
            if (af.str() == "") {
                af << "Far Right";
            } else {
                af << ", Far Right";
            }
        }

        str << "AFPointsInFocus = " << af.str();
        return str.str();
    }
};
NAAFInfoInterpreter naAFInfoInterpreter;

class NALensDataInterpreter : public Interpreter
{
    std::map<std::string, std::string> lenses;
public:
    NALensDataInterpreter ()
    {
        /*  The key is a composite string made of 8 HEX bytes
         *  LensIDNumber LensFStops MinFocalLength MaxFocalLength MaxApertureAtMinFocal MaxApertureAtMaxFocal MCUVersion and LensType */
        lenses["00 00 00 00 00 00 00 01"] = "Manual Lens No CPU";
        lenses["00 00 00 00 00 00 E1 12"] = "TC-17E II";
        lenses["00 00 00 00 00 00 F1 0C"] = "TC-14E [II] or Sigma APO Tele Converter 1.4x EX DG or Kenko Teleplus PRO 300 DG 1.4x";
        lenses["00 00 00 00 00 00 F2 18"] = "TC-20E [II] or Sigma APO Tele Converter 2x EX DG or Kenko Teleplus PRO 300 DG 2.0x";
        lenses["00 00 48 48 53 53 00 01"] = "Loreo 40mm f/11-22 3D Lens in a Cap 9005";
        lenses["00 36 1C 2D 34 3C 00 06"] = "Tamron SP AF 11-18mm f/4.5-5.6 Di II LD Aspherical (IF) (A13)";
        lenses["00 3C 1F 37 30 30 00 06"] = "Tokina AT-X 124 AF PRO DX (AF 12-24mm f/4)";
        lenses["00 3C 2B 44 30 30 00 06"] = "Tokina AT-X 17-35 f/4 PRO FX (AF 17-35mm f/4)";
        lenses["00 3C 5C 80 30 30 00 0E"] = "Tokina AT-X 70-200 f/4 FX VCM-S (AF 70-200mm f/4)";
        lenses["00 3E 80 A0 38 3F 00 02"] = "Tamron SP AF 200-500mm f/5-6.3 Di LD (IF) (A08)";
        lenses["00 3F 2D 80 2B 40 00 06"] = "Tamron AF 18-200mm f/3.5-6.3 XR Di II LD Aspherical (IF) (A14)";
        lenses["00 3F 2D 80 2C 40 00 06"] = "Tamron AF 18-200mm f/3.5-6.3 XR Di II LD Aspherical (IF) Macro (A14)";
        lenses["00 3F 80 A0 38 3F 00 02"] = "Tamron SP AF 200-500mm f/5-6.3 Di (A08)";
        lenses["00 40 11 11 2C 2C 00 00"] = "Samyang 8mm f/3.5 Fish-Eye";
        lenses["00 40 18 2B 2C 34 00 06"] = "Tokina AT-X 107 AF DX Fisheye (AF 10-17mm f/3.5-4.5)";
        lenses["00 40 2A 72 2C 3C 00 06"] = "Tokina AT-X 16.5-135 DX (AF 16.5-135mm f/3.5-5.6)";
        lenses["00 40 2B 2B 2C 2C 00 02"] = "Tokina AT-X 17 AF PRO (AF 17mm f/3.5)";
        lenses["00 40 2D 2D 2C 2C 00 00"] = "Carl Zeiss Distagon T* 3.5/18 ZF.2";
        lenses["00 40 2D 80 2C 40 00 06"] = "Tamron AF 18-200mm f/3.5-6.3 XR Di II LD Aspherical (IF) Macro (A14NII)";
        lenses["00 40 2D 88 2C 40 00 06"] = "Tamron AF 18-250mm f/3.5-6.3 Di II LD Aspherical (IF) Macro (A18NII)";
        lenses["00 40 2D 88 2C 40 62 06"] = "Tamron AF 18-250mm f/3.5-6.3 Di II LD Aspherical (IF) Macro (A18)";
        lenses["00 40 31 31 2C 2C 00 00"] = "Voigtlander Color Skopar 20mm f/3.5 SLII Aspherical";
        lenses["00 40 37 80 2C 3C 00 02"] = "Tokina AT-X 242 AF (AF 24-200mm f/3.5-5.6)";
        lenses["00 40 64 64 2C 2C 00 00"] = "Voigtlander APO-Lanthar 90mm f/3.5 SLII Close Focus";
        lenses["00 44 60 98 34 3C 00 02"] = "Tokina AT-X 840 D (AF 80-400mm f/4.5-5.6)";
        lenses["00 47 10 10 24 24 00 00"] = "Fisheye Nikkor 8mm f/2.8 AiS";
        lenses["00 47 25 25 24 24 00 02"] = "Tamron SP AF 14mm f/2.8 Aspherical (IF) (69E)";
        lenses["00 47 3C 3C 24 24 00 00"] = "Nikkor 28mm f/2.8 AiS";
        lenses["00 47 44 44 24 24 00 06"] = "Tokina AT-X M35 PRO DX (AF 35mm f/2.8 Macro)";
        lenses["00 47 53 80 30 3C 00 06"] = "Tamron AF 55-200mm f/4-5.6 Di II LD (A15)";
        lenses["00 48 1C 29 24 24 00 06"] = "Tokina AT-X 116 PRO DX (AF 11-16mm f/2.8)";
        lenses["00 48 29 3C 24 24 00 06"] = "Tokina AT-X 16-28 AF PRO FX (AF 16-28mm f/2.8)";
        lenses["00 48 29 50 24 24 00 06"] = "Tokina AT-X 165 PRO DX (AF 16-50mm f/2.8)";
        lenses["00 48 32 32 24 24 00 00"] = "Carl Zeiss Distagon T* 2.8/21 ZF.2";
        lenses["00 48 37 5C 24 24 00 06"] = "Tokina AT-X 24-70 f/2.8 PRO FX (AF 24-70mm f/2.8)";
        lenses["00 48 3C 3C 24 24 00 00"] = "Voigtlander Color Skopar 28mm f/2.8 SL II";
        lenses["00 48 3C 60 24 24 00 02"] = "Tokina AT-X 280 AF PRO (AF 28-80mm f/2.8)";
        lenses["00 48 3C 6A 24 24 00 02"] = "Tamron SP AF 28-105mm f/2.8 LD Aspherical IF (176D)";
        lenses["00 48 50 50 18 18 00 00"] = "Nikkor H 50mm f/2";
        lenses["00 48 50 72 24 24 00 06"] = "Tokina AT-X 535 PRO DX (AF 50-135mm f/2.8)";
        lenses["00 48 5C 80 30 30 00 0E"] = "Tokina AT-X 70-200 f/4 FX VCM-S (AF 70-200mm f/4)";
        lenses["00 48 5C 8E 30 3C 00 06"] = "Tamron AF 70-300mm f/4-5.6 Di LD Macro 1:2 (A17NII)";
        lenses["00 48 68 68 24 24 00 00"] = "Series E 100mm f/2.8";
        lenses["00 48 80 80 30 30 00 00"] = "Nikkor 200mm f/4 AiS";
        lenses["00 49 30 48 22 2B 00 02"] = "Tamron SP AF 20-40mm f/2.7-3.5 (166D)";
        lenses["00 4C 6A 6A 20 20 00 00"] = "Nikkor 105mm f/2.5 AiS";
        lenses["00 4C 7C 7C 2C 2C 00 02"] = "Tamron SP AF 180mm f/3.5 Di Model (B01)";
        lenses["00 53 2B 50 24 24 00 06"] = "Tamron SP AF 17-50mm f/2.8 XR Di II LD Aspherical (IF) (A16)";
        lenses["00 54 2B 50 24 24 00 06"] = "Tamron SP AF 17-50mm f/2.8 XR Di II LD Aspherical (IF) (A16NII)";
        lenses["00 54 3C 3C 18 18 00 00"] = "Carl Zeiss Distagon T* 2/28 ZF.2";
        lenses["00 54 44 44 0C 0C 00 00"] = "Carl Zeiss Distagon T* 1.4/35 ZF.2";
        lenses["00 54 44 44 18 18 00 00"] = "Carl Zeiss Distagon T* 2/35 ZF.2";
        lenses["00 54 48 48 18 18 00 00"] = "Voigtlander Ultron 40mm f/2 SLII Aspherical";
        lenses["00 54 50 50 0C 0C 00 00"] = "Carl Zeiss Planar T* 1.4/50 ZF.2";
        lenses["00 54 50 50 18 18 00 00"] = "Carl Zeiss Makro-Planar T* 2/50 ZF.2";
        lenses["00 54 53 53 0C 0C 00 00"] = "Zeiss Otus 1.4/55";
        lenses["00 54 55 55 0C 0C 00 00"] = "Voigtlander Nokton 58mm f/1.4 SLII";
        lenses["00 54 56 56 30 30 00 00"] = "Coastal Optical Systems 60mm 1:4 UV-VIS-IR Macro Apo";
        lenses["00 54 62 62 0C 0C 00 00"] = "Carl Zeiss Planar T* 1.4/85 ZF.2";
        lenses["00 54 68 68 18 18 00 00"] = "Carl Zeiss Makro-Planar T* 2/100 ZF.2";
        lenses["00 54 68 68 24 24 00 02"] = "Tokina AT-X M100 AF PRO D (AF 100mm f/2.8 Macro)";
        lenses["00 54 72 72 18 18 00 00"] = "Carl Zeiss Apo Sonnar T* 2/135 ZF.2";
        lenses["00 54 8E 8E 24 24 00 02"] = "Tokina AT-X 300 AF PRO (AF 300mm f/2.8)";
        lenses["00 57 50 50 14 14 00 00"] = "Nikkor 50mm f/1.8 AI";
        lenses["00 58 64 64 20 20 00 00"] = "Soligor C/D Macro MC 90mm f/2.5";
        lenses["01 00 00 00 00 00 02 00"] = "TC-16A";
        lenses["01 00 00 00 00 00 08 00"] = "TC-16A";
        lenses["01 54 62 62 0C 0C 00 00"] = "Zeiss Otus 1.4/85";
        lenses["01 58 50 50 14 14 02 00"] = "AF Nikkor 50mm f/1.8";
        lenses["01 58 50 50 14 14 05 00"] = "AF Nikkor 50mm f/1.8";
        lenses["02 2F 98 98 3D 3D 02 00"] = "Sigma APO 400mm f/5.6";
        lenses["02 34 A0 A0 44 44 02 00"] = "Sigma APO 500mm f/7.2";
        lenses["02 37 5E 8E 35 3D 02 00"] = "Sigma 75-300mm f/4.5-5.6 APO";
        lenses["02 37 A0 A0 34 34 02 00"] = "Sigma APO 500mm f/4.5";
        lenses["02 3A 37 50 31 3D 02 00"] = "Sigma 24-50mm f/4-5.6 UC";
        lenses["02 3A 5E 8E 32 3D 02 00"] = "Sigma 75-300mm f/4.0-5.6";
        lenses["02 3B 44 61 30 3D 02 00"] = "Sigma 35-80mm f/4-5.6";
        lenses["02 3C B0 B0 3C 3C 02 00"] = "Sigma APO 800mm f/5.6";
        lenses["02 3F 24 24 2C 2C 02 00"] = "Sigma 14mm f/3.5";
        lenses["02 3F 3C 5C 2D 35 02 00"] = "Sigma 28-70mm f/3.5-4.5 UC";
        lenses["02 40 44 5C 2C 34 02 00"] = "Exakta AF 35-70mm 1:3.5-4.5 MC";
        lenses["02 40 44 73 2B 36 02 00"] = "Sigma 35-135mm f/3.5-4.5 a";
        lenses["02 40 5C 82 2C 35 02 00"] = "Sigma APO 70-210mm f/3.5-4.5";
        lenses["02 42 44 5C 2A 34 02 00"] = "AF Zoom-Nikkor 35-70mm f/3.3-4.5";
        lenses["02 42 44 5C 2A 34 08 00"] = "AF Zoom-Nikkor 35-70mm f/3.3-4.5";
        lenses["02 46 37 37 25 25 02 00"] = "Sigma 24mm f/2.8 Super Wide II Macro";
        lenses["02 46 3C 5C 25 25 02 00"] = "Sigma 28-70mm f/2.8";
        lenses["02 46 5C 82 25 25 02 00"] = "Sigma 70-210mm f/2.8 APO";
        lenses["02 48 50 50 24 24 02 00"] = "Sigma Macro 50mm f/2.8";
        lenses["02 48 65 65 24 24 02 00"] = "Sigma Macro 90mm f/2.8";
        lenses["03 43 5C 81 35 35 02 00"] = "Soligor AF C/D Zoom UMCS 70-210mm 1:4.5";
        lenses["03 48 5C 81 30 30 02 00"] = "AF Zoom-Nikkor 70-210mm f/4";
        lenses["04 48 3C 3C 24 24 03 00"] = "AF Nikkor 28mm f/2.8";
        lenses["05 54 50 50 0C 0C 04 00"] = "AF Nikkor 50mm f/1.4";
        lenses["06 3F 68 68 2C 2C 06 00"] = "Cosina AF 100mm f/3.5 Macro";
        lenses["06 54 53 53 24 24 06 00"] = "AF Micro-Nikkor 55mm f/2.8";
        lenses["07 36 3D 5F 2C 3C 03 00"] = "Cosina AF Zoom 28-80mm f/3.5-5.6 MC Macro";
        lenses["07 3E 30 43 2D 35 03 00"] = "Soligor AF Zoom 19-35mm 1:3.5-4.5 MC";
        lenses["07 40 2F 44 2C 34 03 02"] = "Tamron AF 19-35mm f/3.5-4.5 (A10)";
        lenses["07 40 30 45 2D 35 03 02"] = "Tamron AF 19-35mm f/3.5-4.5 (A10)";
        lenses["07 40 3C 5C 2C 35 03 00"] = "Tokina AF 270 II (AF 28-70mm f/3.5-4.5)";
        lenses["07 40 3C 62 2C 34 03 00"] = "AF Zoom-Nikkor 28-85mm f/3.5-4.5";
        lenses["07 46 2B 44 24 30 03 02"] = "Tamron SP AF 17-35mm f/2.8-4 Di LD Aspherical (IF) (A05)";
        lenses["07 46 3D 6A 25 2F 03 00"] = "Cosina AF Zoom 28-105mm f/2.8-3.8 MC";
        lenses["07 47 3C 5C 25 35 03 00"] = "Tokina AF 287 SD (AF 28-70mm f/2.8-4.5)";
        lenses["07 48 3C 5C 24 24 03 00"] = "Tokina AT-X 287 AF (AF 28-70mm f/2.8)";
        lenses["08 40 44 6A 2C 34 04 00"] = "AF Zoom-Nikkor 35-105mm f/3.5-4.5";
        lenses["09 48 37 37 24 24 04 00"] = "AF Nikkor 24mm f/2.8";
        lenses["0A 48 8E 8E 24 24 03 00"] = "AF Nikkor 300mm f/2.8 IF-ED";
        lenses["0A 48 8E 8E 24 24 05 00"] = "AF Nikkor 300mm f/2.8 IF-ED N";
        lenses["0B 3E 3D 7F 2F 3D 0E 00"] = "Tamron AF 28-200mm f/3.8-5.6 (71D)";
        lenses["0B 3E 3D 7F 2F 3D 0E 02"] = "Tamron AF 28-200mm f/3.8-5.6D (171D)";
        lenses["0B 48 7C 7C 24 24 05 00"] = "AF Nikkor 180mm f/2.8 IF-ED";
        lenses["0D 40 44 72 2C 34 07 00"] = "AF Zoom-Nikkor 35-135mm f/3.5-4.5";
        lenses["0E 48 5C 81 30 30 05 00"] = "AF Zoom-Nikkor 70-210mm f/4";
        lenses["0E 4A 31 48 23 2D 0E 02"] = "Tamron SP AF 20-40mm f/2.7-3.5 (166D)";
        lenses["0F 58 50 50 14 14 05 00"] = "AF Nikkor 50mm f/1.8 N";
        lenses["10 3D 3C 60 2C 3C D2 02"] = "Tamron AF 28-80mm f/3.5-5.6 Aspherical (177D)";
        lenses["10 48 8E 8E 30 30 08 00"] = "AF Nikkor 300mm f/4 IF-ED";
        lenses["11 48 44 5C 24 24 08 00"] = "AF Zoom-Nikkor 35-70mm f/2.8";
        lenses["12 36 5C 81 35 3D 09 00"] = "Cosina AF Zoom 70-210mm f/4.5-5.6 MC Macro";
        lenses["12 36 69 97 35 42 09 00"] = "Soligor AF Zoom 100-400mm 1:4.5-6.7 MC";
        lenses["12 38 69 97 35 42 09 02"] = "Promaster Spectrum 7 100-400mm f/4.5-6.7";
        lenses["12 39 5C 8E 34 3D 08 02"] = "Cosina AF Zoom 70-300mm f/4.5-5.6 MC Macro";
        lenses["12 3B 68 8D 3D 43 09 02"] = "Cosina AF Zoom 100-300mm f/5.6-6.7 MC Macro";
        lenses["12 3B 98 98 3D 3D 09 00"] = "Tokina AT-X 400 AF SD (AF 400mm f/5.6)";
        lenses["12 3D 3C 80 2E 3C DF 02"] = "Tamron AF 28-200mm f/3.8-5.6 AF Aspherical LD (IF) (271D)";
        lenses["12 44 5E 8E 34 3C 09 00"] = "Tokina AF 730 (AF 75-300mm f/4.5-5.6)";
        lenses["12 48 5C 81 30 3C 09 00"] = "AF Nikkor 70-210mm f/4-5.6";
        lenses["12 4A 5C 81 31 3D 09 00"] = "Soligor AF C/D Auto Zoom+Macro 70-210mm 1:4-5.6 UMCS";
        lenses["13 42 37 50 2A 34 0B 00"] = "AF Zoom-Nikkor 24-50mm f/3.3-4.5";
        lenses["14 48 60 80 24 24 0B 00"] = "AF Zoom-Nikkor 80-200mm f/2.8 ED";
        lenses["14 48 68 8E 30 30 0B 00"] = "Tokina AT-X 340 AF (AF 100-300mm f/4)";
        lenses["14 54 60 80 24 24 0B 00"] = "Tokina AT-X 828 AF (AF 80-200mm f/2.8)";
        lenses["15 4C 62 62 14 14 0C 00"] = "AF Nikkor 85mm f/1.8";
        lenses["17 3C A0 A0 30 30 0F 00"] = "Nikkor 500mm f/4 P ED IF";
        lenses["17 3C A0 A0 30 30 11 00"] = "Nikkor 500mm f/4 P ED IF";
        lenses["18 40 44 72 2C 34 0E 00"] = "AF Zoom-Nikkor 35-135mm f/3.5-4.5 N";
        lenses["1A 54 44 44 18 18 11 00"] = "AF Nikkor 35mm f/2";
        lenses["1B 44 5E 8E 34 3C 10 00"] = "AF Zoom-Nikkor 75-300mm f/4.5-5.6";
        lenses["1C 48 30 30 24 24 12 00"] = "AF Nikkor 20mm f/2.8";
        lenses["1D 42 44 5C 2A 34 12 00"] = "AF Zoom-Nikkor 35-70mm f/3.3-4.5 N";
        lenses["1E 54 56 56 24 24 13 00"] = "AF Micro-Nikkor 60mm f/2.8";
        lenses["1E 5D 64 64 20 20 13 00"] = "Tamron SP AF 90mm f/2.5 (52E)";
        lenses["1F 54 6A 6A 24 24 14 00"] = "AF Micro-Nikkor 105mm f/2.8";
        lenses["20 3C 80 98 3D 3D 1E 02"] = "Tamron AF 200-400mm f/5.6 LD IF (75D)";
        lenses["20 48 60 80 24 24 15 00"] = "AF Zoom-Nikkor 80-200mm f/2.8 ED";
        lenses["20 5A 64 64 20 20 14 00"] = "Tamron SP AF 90mm f/2.5 Macro (152E)";
        lenses["21 40 3C 5C 2C 34 16 00"] = "AF Zoom-Nikkor 28-70mm f/3.5-4.5";
        lenses["21 56 8E 8E 24 24 14 00"] = "Tamron SP AF 300mm f/2.8 LD-IF (60E)";
        lenses["22 48 72 72 18 18 16 00"] = "AF DC-Nikkor 135mm f/2";
        lenses["22 53 64 64 24 24 E0 02"] = "Tamron SP AF 90mm f/2.8 Macro 1:1 (72E)";
        lenses["23 30 BE CA 3C 48 17 00"] = "Zoom-Nikkor 1200-1700mm f/5.6-8 P ED IF";
        lenses["24 44 60 98 34 3C 1A 02"] = "Tokina AT-X 840 AF-II (AF 80-400mm f/4.5-5.6)";
        lenses["24 48 60 80 24 24 1A 02"] = "AF Zoom-Nikkor 80-200mm f/2.8D ED";
        lenses["24 54 60 80 24 24 1A 02"] = "Tokina AT-X 828 AF PRO (AF 80-200mm f/2.8)";
        lenses["25 44 44 8E 34 42 1B 02"] = "Tokina AF 353 (AF 35-300mm f/4.5-6.7)";
        lenses["25 48 3C 5C 24 24 1B 02"] = "Tokina AT-X 270 AF PRO II (AF 28-70mm f/2.6-2.8)";
        lenses["25 48 3C 5C 24 24 1B 02"] = "Tokina AT-X 287 AF PRO SV (AF 28-70mm f/2.8)";
        lenses["25 48 44 5C 24 24 1B 02"] = "AF Zoom-Nikkor 35-70mm f/2.8D";
        lenses["25 48 44 5C 24 24 3A 02"] = "AF Zoom-Nikkor 35-70mm f/2.8D";
        lenses["25 48 44 5C 24 24 52 02"] = "AF Zoom-Nikkor 35-70mm f/2.8D";
        lenses["26 3C 54 80 30 3C 1C 06"] = "Sigma 55-200mm f/4-5.6 DC";
        lenses["26 3C 5C 82 30 3C 1C 02"] = "Sigma 70-210mm f/4-5.6 UC-II";
        lenses["26 3C 5C 8E 30 3C 1C 02"] = "Sigma 70-300mm f/4-5.6 DG Macro";
        lenses["26 3C 98 98 3C 3C 1C 02"] = "Sigma APO Tele Macro 400mm f/5.6";
        lenses["26 3D 3C 80 2F 3D 1C 02"] = "Sigma 28-300mm f/3.8-5.6 Aspherical";
        lenses["26 3E 3C 6A 2E 3C 1C 02"] = "Sigma 28-105mm f/3.8-5.6 UC-III Aspherical IF";
        lenses["26 40 27 3F 2C 34 1C 02"] = "Sigma 15-30mm f/3.5-4.5 EX DG Aspherical DF";
        lenses["26 40 2D 44 2B 34 1C 02"] = "Sigma 18-35mm f/3.5-4.5 Aspherical";
        lenses["26 40 2D 50 2C 3C 1C 06"] = "Sigma 18-50mm f/3.5-5.6 DC";
        lenses["26 40 2D 70 2B 3C 1C 06"] = "Sigma 18-125mm f/3.5-5.6 DC";
        lenses["26 40 2D 80 2C 40 1C 06"] = "Sigma 18-200mm f/3.5-6.3 DC";
        lenses["26 40 37 5C 2C 3C 1C 02"] = "Sigma 24-70mm f/3.5-5.6 Aspherical HF";
        lenses["26 40 3C 5C 2C 34 1C 02"] = "AF Zoom-Nikkor 28-70mm f/3.5-4.5D";
        lenses["26 40 3C 60 2C 3C 1C 02"] = "Sigma 28-80mm f/3.5-5.6 Mini Zoom Macro II Aspherical";
        lenses["26 40 3C 65 2C 3C 1C 02"] = "Sigma 28-90mm f/3.5-5.6 Macro";
        lenses["26 40 3C 80 2B 3C 1C 02"] = "Sigma 28-200mm f/3.5-5.6 Compact Aspherical Hyperzoom Macro";
        lenses["26 40 3C 80 2C 3C 1C 02"] = "Sigma 28-200mm f/3.5-5.6 Compact Aspherical Hyperzoom Macro";
        lenses["26 40 3C 8E 2C 40 1C 02"] = "Sigma 28-300mm f/3.5-6.3 Macro";
        lenses["26 40 7B A0 34 40 1C 02"] = "Sigma APO 170-500mm f/5-6.3 Aspherical RF";
        lenses["26 41 3C 8E 2C 40 1C 02"] = "Sigma 28-300mm f/3.5-6.3 DG Macro";
        lenses["26 44 73 98 34 3C 1C 02"] = "Sigma 135-400mm f/4.5-5.6 APO Aspherical";
        lenses["26 48 11 11 30 30 1C 02"] = "Sigma 8mm f/4 EX Circular Fisheye";
        lenses["26 48 27 27 24 24 1C 02"] = "Sigma 15mm f/2.8 EX Diagonal Fisheye";
        lenses["26 48 2D 50 24 24 1C 06"] = "Sigma 18-50mm f/2.8 EX DC";
        lenses["26 48 31 49 24 24 1C 02"] = "Sigma 20-40mm f/2.8";
        lenses["26 48 37 56 24 24 1C 02"] = "Sigma 24-60mm f/2.8 EX DG";
        lenses["26 48 3C 5C 24 24 1C 06"] = "Sigma 28-70mm f/2.8 EX DG";
        lenses["26 48 3C 5C 24 30 1C 02"] = "Sigma 28-70mm f/2.8-4 DG";
        lenses["26 48 3C 6A 24 30 1C 02"] = "Sigma 28-105mm f/2.8-4 Aspherical";
        lenses["26 48 8E 8E 30 30 1C 02"] = "Sigma APO Tele Macro 300mm f/4";
        lenses["26 54 2B 44 24 30 1C 02"] = "Sigma 17-35mm f/2.8-4 EX Aspherical";
        lenses["26 54 37 5C 24 24 1C 02"] = "Sigma 24-70mm f/2.8 EX DG Macro";
        lenses["26 54 37 73 24 34 1C 02"] = "Sigma 24-135mm f/2.8-4.5";
        lenses["26 54 3C 5C 24 24 1C 02"] = "Sigma 28-70mm f/2.8 EX";
        lenses["26 58 31 31 14 14 1C 02"] = "Sigma 20mm f/1.8 EX DG Aspherical RF";
        lenses["26 58 37 37 14 14 1C 02"] = "Sigma 24mm f/1.8 EX DG Aspherical Macro";
        lenses["26 58 3C 3C 14 14 1C 02"] = "Sigma 28mm f/1.8 EX DG Aspherical Macro";
        lenses["27 48 8E 8E 24 24 1D 02"] = "AF-I Nikkor 300mm f/2.8D IF-ED";
        lenses["27 48 8E 8E 24 24 E1 02"] = "AF-I Nikkor 300mm f/2.8D IF-ED + TC-17E";
        lenses["27 48 8E 8E 24 24 F1 02"] = "AF-I Nikkor 300mm f/2.8D IF-ED + TC-14E";
        lenses["27 48 8E 8E 24 24 F2 02"] = "AF-I Nikkor 300mm f/2.8D IF-ED + TC-20E";
        lenses["27 48 8E 8E 30 30 1D 02"] = "Tokina AT-X 304 AF (AF 300mm f/4.0)";
        lenses["27 54 8E 8E 24 24 1D 02"] = "Tamron SP AF 300mm f/2.8 LD-IF (360E)";
        lenses["28 3C A6 A6 30 30 1D 02"] = "AF-I Nikkor 600mm f/4D IF-ED";
        lenses["28 3C A6 A6 30 30 E1 02"] = "AF-I Nikkor 600mm f/4D IF-ED + TC-17E";
        lenses["28 3C A6 A6 30 30 F1 02"] = "AF-I Nikkor 600mm f/4D IF-ED + TC-14E";
        lenses["28 3C A6 A6 30 30 F2 02"] = "AF-I Nikkor 600mm f/4D IF-ED + TC-20E";
        lenses["2A 54 3C 3C 0C 0C 26 02"] = "AF Nikkor 28mm f/1.4D";
        lenses["2B 3C 44 60 30 3C 1F 02"] = "AF Zoom-Nikkor 35-80mm f/4-5.6D";
        lenses["2C 48 6A 6A 18 18 27 02"] = "AF DC-Nikkor 105mm f/2D";
        lenses["2D 48 80 80 30 30 21 02"] = "AF Micro-Nikkor 200mm f/4D IF-ED";
        lenses["2E 48 5C 82 30 3C 22 02"] = "AF Nikkor 70-210mm f/4-5.6D";
        lenses["2E 48 5C 82 30 3C 28 02"] = "AF Nikkor 70-210mm f/4-5.6D";
        lenses["2F 40 30 44 2C 34 29 02"] = "Tokina AF 235 II (AF 20-35mm f/3.5-4.5)";
        lenses["2F 40 30 44 2C 34 29 02"] = "Tokina AF 193 (AF 19-35mm f/3.5-4.5)";
        lenses["2F 48 30 44 24 24 29 02"] = "AF Zoom-Nikkor 20-35mm f/2.8D IF";
        lenses["2F 48 30 44 24 24 29 02"] = "Tokina AT-X 235 AF PRO (AF 20-35mm f/2.8)";
        lenses["30 48 98 98 24 24 24 02"] = "AF-I Nikkor 400mm f/2.8D IF-ED";
        lenses["30 48 98 98 24 24 E1 02"] = "AF-I Nikkor 400mm f/2.8D IF-ED + TC-17E";
        lenses["30 48 98 98 24 24 F1 02"] = "AF-I Nikkor 400mm f/2.8D IF-ED + TC-14E";
        lenses["30 48 98 98 24 24 F2 02"] = "AF-I Nikkor 400mm f/2.8D IF-ED + TC-20E";
        lenses["31 54 56 56 24 24 25 02"] = "AF Micro-Nikkor 60mm f/2.8D";
        lenses["32 53 64 64 24 24 35 02"] = "Tamron SP AF 90mm f/2.8 [Di] Macro 1:1 (172E/272E)";
        lenses["32 54 50 50 24 24 35 02"] = "Sigma Macro 50mm f/2.8 EX DG";
        lenses["32 54 6A 6A 24 24 35 02"] = "AF Micro-Nikkor 105mm f/2.8D";
        lenses["32 54 6A 6A 24 24 35 02"] = "Sigma Macro 105mm f/2.8 EX DG";
        lenses["33 48 2D 2D 24 24 31 02"] = "AF Nikkor 18mm f/2.8D";
        lenses["33 54 3C 5E 24 24 62 02"] = "Tamron SP AF 28-75mm f/2.8 XR Di LD Aspherical (IF) Macro (A09)";
        lenses["34 48 29 29 24 24 32 02"] = "AF Fisheye Nikkor 16mm f/2.8D";
        lenses["35 3C A0 A0 30 30 33 02"] = "AF-I Nikkor 500mm f/4D IF-ED";
        lenses["35 3C A0 A0 30 30 E1 02"] = "AF-I Nikkor 500mm f/4D IF-ED + TC-17E";
        lenses["35 3C A0 A0 30 30 F1 02"] = "AF-I Nikkor 500mm f/4D IF-ED + TC-14E";
        lenses["35 3C A0 A0 30 30 F2 02"] = "AF-I Nikkor 500mm f/4D IF-ED + TC-20E";
        lenses["36 48 37 37 24 24 34 02"] = "AF Nikkor 24mm f/2.8D";
        lenses["37 48 30 30 24 24 36 02"] = "AF Nikkor 20mm f/2.8D";
        lenses["38 4C 62 62 14 14 37 02"] = "AF Nikkor 85mm f/1.8D";
        lenses["3A 40 3C 5C 2C 34 39 02"] = "AF Zoom-Nikkor 28-70mm f/3.5-4.5D";
        lenses["3B 48 44 5C 24 24 3A 02"] = "AF Zoom-Nikkor 35-70mm f/2.8D N";
        lenses["3C 48 60 80 24 24 3B 02"] = "AF Zoom-Nikkor 80-200mm f/2.8D ED";
        lenses["3D 3C 44 60 30 3C 3E 02"] = "AF Zoom-Nikkor 35-80mm f/4-5.6D";
        lenses["3E 48 3C 3C 24 24 3D 02"] = "AF Nikkor 28mm f/2.8D";
        lenses["3F 40 44 6A 2C 34 45 02"] = "AF Zoom-Nikkor 35-105mm f/3.5-4.5D";
        lenses["41 48 7C 7C 24 24 43 02"] = "AF Nikkor 180mm f/2.8D IF-ED";
        lenses["42 54 44 44 18 18 44 02"] = "AF Nikkor 35mm f/2D";
        lenses["43 54 50 50 0C 0C 46 02"] = "AF Nikkor 50mm f/1.4D";
        lenses["44 44 60 80 34 3C 47 02"] = "AF Zoom-Nikkor 80-200mm f/4.5-5.6D";
        lenses["45 3D 3C 60 2C 3C 48 02"] = "Tamron AF 28-80mm f/3.5-5.6 Aspherical (177D)";
        lenses["45 40 3C 60 2C 3C 48 02"] = "AF Zoom-Nikkor 28-80mm f/3.5-5.6D";
        lenses["45 41 37 72 2C 3C 48 02"] = "Tamron SP AF 24-135mm f/3.5-5.6 AD Aspherical (IF) Macro (190D)";
        lenses["46 3C 44 60 30 3C 49 02"] = "AF Zoom-Nikkor 35-80mm f/4-5.6D N";
        lenses["47 42 37 50 2A 34 4A 02"] = "AF Zoom-Nikkor 24-50mm f/3.3-4.5D";
        lenses["48 38 1F 37 34 3C 4B 06"] = "Sigma 12-24mm f/4.5-5.6 EX DG Aspherical HSM";
        lenses["48 3C 19 31 30 3C 4B 06"] = "Sigma 10-20mm f/4-5.6 EX DC HSM";
        lenses["48 3C 50 A0 30 40 4B 02"] = "Sigma 50-500mm f/4-6.3 EX APO RF HSM";
        lenses["48 3C 8E B0 3C 3C 4B 02"] = "Sigma APO 300-800mm f/5.6 EX DG HSM";
        lenses["48 3C B0 B0 3C 3C 4B 02"] = "Sigma APO 800mm f/5.6 EX HSM";
        lenses["48 44 A0 A0 34 34 4B 02"] = "Sigma APO 500mm f/4.5 EX HSM";
        lenses["48 48 24 24 24 24 4B 02"] = "Sigma 14mm f/2.8 EX Aspherical HSM";
        lenses["48 48 2B 44 24 30 4B 06"] = "Sigma 17-35mm f/2.8-4 EX DG  Aspherical HSM";
        lenses["48 48 68 8E 30 30 4B 02"] = "Sigma APO 100-300mm f/4 EX IF HSM";
        lenses["48 48 76 76 24 24 4B 06"] = "Sigma APO Macro 150mm f/2.8 EX DG HSM";
        lenses["48 48 8E 8E 24 24 4B 02"] = "AF-S Nikkor 300mm f/2.8D IF-ED";
        lenses["48 48 8E 8E 24 24 E1 02"] = "AF-S Nikkor 300mm f/2.8D IF-ED + TC-17E";
        lenses["48 48 8E 8E 24 24 F1 02"] = "AF-S Nikkor 300mm f/2.8D IF-ED + TC-14E";
        lenses["48 48 8E 8E 24 24 F2 02"] = "AF-S Nikkor 300mm f/2.8D IF-ED + TC-20E";
        lenses["48 4C 7C 7C 2C 2C 4B 02"] = "Sigma APO Macro 180mm f/3.5 EX DG HSM";
        lenses["48 4C 7D 7D 2C 2C 4B 02"] = "Sigma APO Macro 180mm f/3.5 EX DG HSM";
        lenses["48 54 3E 3E 0C 0C 4B 06"] = "Sigma 30mm f/1.4 EX DC HSM";
        lenses["48 54 5C 80 24 24 4B 02"] = "Sigma 70-200mm f/2.8 EX APO IF HSM";
        lenses["48 54 6F 8E 24 24 4B 02"] = "Sigma APO 120-300mm f/2.8 EX DG HSM";
        lenses["48 54 8E 8E 24 24 4B 02"] = "Sigma APO 300mm f/2.8 EX DG HSM";
        lenses["49 3C A6 A6 30 30 4C 02"] = "AF-S Nikkor 600mm f/4D IF-ED";
        lenses["49 3C A6 A6 30 30 E1 02"] = "AF-S Nikkor 600mm f/4D IF-ED + TC-17E";
        lenses["49 3C A6 A6 30 30 F1 02"] = "AF-S Nikkor 600mm f/4D IF-ED + TC-14E";
        lenses["49 3C A6 A6 30 30 F2 02"] = "AF-S Nikkor 600mm f/4D IF-ED + TC-20E";
        lenses["4A 40 11 11 2C 0C 4D 02"] = "Samyang 8mm f/3.5 Fish-Eye CS";
        lenses["4A 48 1E 1E 24 0C 4D 02"] = "Samyang 12mm f/2.8 ED AS NCS Fish-Eye";
        lenses["4A 48 24 24 24 0C 4D 02"] = "Samyang AE 14mm f/2.8 ED AS IF UMC";
        lenses["4A 54 29 29 18 0C 4D 02"] = "Samyang 16mm f/2.0 ED AS UMC CS";
        lenses["4A 54 62 62 0C 0C 4D 02"] = "AF Nikkor 85mm f/1.4D IF";
        lenses["4A 60 44 44 0C 0C 4D 02"] = "Samyang 35mm f/1.4 AS UMC";
        lenses["4A 60 62 62 0C 0C 4D 02"] = "Samyang AE 85mm f/1.4 AS IF UMC";
        lenses["4B 3C A0 A0 30 30 4E 02"] = "AF-S Nikkor 500mm f/4D IF-ED";
        lenses["4B 3C A0 A0 30 30 E1 02"] = "AF-S Nikkor 500mm f/4D IF-ED + TC-17E";
        lenses["4B 3C A0 A0 30 30 F1 02"] = "AF-S Nikkor 500mm f/4D IF-ED + TC-14E";
        lenses["4B 3C A0 A0 30 30 F2 02"] = "AF-S Nikkor 500mm f/4D IF-ED + TC-20E";
        lenses["4C 40 37 6E 2C 3C 4F 02"] = "AF Zoom-Nikkor 24-120mm f/3.5-5.6D IF";
        lenses["4D 3E 3C 80 2E 3C 62 02"] = "Tamron AF 28-200mm f/3.8-5.6 XR Aspherical (IF) Macro (A03N)";
        lenses["4D 40 3C 80 2C 3C 62 02"] = "AF Zoom-Nikkor 28-200mm f/3.5-5.6D IF";
        lenses["4D 41 3C 8E 2B 40 62 02"] = "Tamron AF 28-300mm f/3.5-6.3 XR Di LD Aspherical (IF) (A061)";
        lenses["4D 41 3C 8E 2C 40 62 02"] = "Tamron AF 28-300mm f/3.5-6.3 XR LD Aspherical (IF) (185D)";
        lenses["4E 48 72 72 18 18 51 02"] = "AF DC-Nikkor 135mm f/2D";
        lenses["4F 40 37 5C 2C 3C 53 06"] = "IX-Nikkor 24-70mm f/3.5-5.6";
        lenses["50 48 56 7C 30 3C 54 06"] = "IX-Nikkor 60-180mm f/4-5.6";
        lenses["52 54 44 44 18 18 00 00"] = "Zeiss Milvus 35mm f/2";
        lenses["53 48 60 80 24 24 57 02"] = "AF Zoom-Nikkor 80-200mm f/2.8D ED";
        lenses["53 48 60 80 24 24 60 02"] = "AF Zoom-Nikkor 80-200mm f/2.8D ED";
        lenses["53 54 50 50 0C 0C 00 00"] = "Zeiss Milvus 50mm f/1.4";
        lenses["54 44 5C 7C 34 3C 58 02"] = "AF Zoom-Micro Nikkor 70-180mm f/4.5-5.6D ED";
        lenses["54 44 5C 7C 34 3C 61 02"] = "AF Zoom-Micro Nikkor 70-180mm f/4.5-5.6D ED";
        lenses["54 54 50 50 18 18 00 00"] = "Zeiss Milvus 50mm f/2 Macro";
        lenses["55 54 62 62 0C 0C 00 00"] = "Zeiss Milvus 85mm f/1.4";
        lenses["56 3C 5C 8E 30 3C 1C 02"] = "Sigma 70-300mm f/4-5.6 APO Macro Super II";
        lenses["56 48 5C 8E 30 3C 5A 02"] = "AF Zoom-Nikkor 70-300mm f/4-5.6D ED";
        lenses["56 54 68 68 18 18 00 00"] = "Zeiss Milvus 100mm f/2 Macro";
        lenses["59 48 98 98 24 24 5D 02"] = "AF-S Nikkor 400mm f/2.8D IF-ED";
        lenses["59 48 98 98 24 24 E1 02"] = "AF-S Nikkor 400mm f/2.8D IF-ED + TC-17E";
        lenses["59 48 98 98 24 24 F1 02"] = "AF-S Nikkor 400mm f/2.8D IF-ED + TC-14E";
        lenses["59 48 98 98 24 24 F2 02"] = "AF-S Nikkor 400mm f/2.8D IF-ED + TC-20E";
        lenses["5A 3C 3E 56 30 3C 5E 06"] = "IX-Nikkor 30-60mm f/4-5.6";
        lenses["5B 44 56 7C 34 3C 5F 06"] = "IX-Nikkor 60-180mm f/4.5-5.6";
        lenses["5D 48 3C 5C 24 24 63 02"] = "AF-S Zoom-Nikkor 28-70mm f/2.8D IF-ED";
        lenses["5E 48 60 80 24 24 64 02"] = "AF-S Zoom-Nikkor 80-200mm f/2.8D IF-ED";
        lenses["5F 40 3C 6A 2C 34 65 02"] = "AF Zoom-Nikkor 28-105mm f/3.5-4.5D IF";
        lenses["60 40 3C 60 2C 3C 66 02"] = "AF Zoom-Nikkor 28-80mm f/3.5-5.6D";
        lenses["61 44 5E 86 34 3C 67 02"] = "AF Zoom-Nikkor 75-240mm f/4.5-5.6D";
        lenses["63 48 2B 44 24 24 68 02"] = "AF-S Nikkor 17-35mm f/2.8D IF-ED";
        lenses["64 00 62 62 24 24 6A 02"] = "PC Micro-Nikkor 85mm f/2.8D";
        lenses["65 44 60 98 34 3C 6B 0A"] = "AF VR Zoom-Nikkor 80-400mm f/4.5-5.6D ED";
        lenses["66 40 2D 44 2C 34 6C 02"] = "AF Zoom-Nikkor 18-35mm f/3.5-4.5D IF-ED";
        lenses["67 48 37 62 24 30 6D 02"] = "AF Zoom-Nikkor 24-85mm f/2.8-4D IF";
        lenses["67 54 37 5C 24 24 1C 02"] = "Sigma 24-70mm f/2.8 EX DG Macro";
        lenses["68 42 3C 60 2A 3C 6E 06"] = "AF Zoom-Nikkor 28-80mm f/3.3-5.6G";
        lenses["69 47 5C 8E 30 3C 00 02"] = "Tamron AF 70-300mm f/4-5.6 Di LD Macro 1:2 (A17N)";
        lenses["69 48 5C 8E 30 3C 6F 02"] = "Tamron AF 70-300mm f/4-5.6 LD Macro 1:2 (572D/772D)";
        lenses["69 48 5C 8E 30 3C 6F 06"] = "AF Zoom-Nikkor 70-300mm f/4-5.6G";
        lenses["6A 48 8E 8E 30 30 70 02"] = "AF-S Nikkor 300mm f/4D IF-ED";
        lenses["6B 48 24 24 24 24 71 02"] = "AF Nikkor ED 14mm f/2.8D";
        lenses["6D 48 8E 8E 24 24 73 02"] = "AF-S Nikkor 300mm f/2.8D IF-ED II";
        lenses["6E 48 98 98 24 24 74 02"] = "AF-S Nikkor 400mm f/2.8D IF-ED II";
        lenses["6F 3C A0 A0 30 30 75 02"] = "AF-S Nikkor 500mm f/4D IF-ED II";
        lenses["70 3C A6 A6 30 30 76 02"] = "AF-S Nikkor 600mm f/4D IF-ED II";
        lenses["72 48 4C 4C 24 24 77 00"] = "Nikkor 45mm f/2.8 P";
        lenses["74 40 37 62 2C 34 78 06"] = "AF-S Zoom-Nikkor 24-85mm f/3.5-4.5G IF-ED";
        lenses["75 40 3C 68 2C 3C 79 06"] = "AF Zoom-Nikkor 28-100mm f/3.5-5.6G";
        lenses["76 58 50 50 14 14 7A 02"] = "AF Nikkor 50mm f/1.8D";
        lenses["77 44 61 98 34 3C 7B 0E"] = "Sigma 80-400mm f/4.5-5.6 EX OS";
        lenses["77 48 5C 80 24 24 7B 0E"] = "AF-S VR Zoom-Nikkor 70-200mm f/2.8G IF-ED";
        lenses["78 40 37 6E 2C 3C 7C 0E"] = "AF-S VR Zoom-Nikkor 24-120mm f/3.5-5.6G IF-ED";
        lenses["79 40 11 11 2C 2C 1C 06"] = "Sigma 8mm f/3.5 EX Circular Fisheye";
        lenses["79 40 3C 80 2C 3C 7F 06"] = "AF Zoom-Nikkor 28-200mm f/3.5-5.6G IF-ED";
        lenses["79 48 3C 5C 24 24 1C 06"] = "Sigma 28-70mm f/2.8 EX DG";
        lenses["79 48 5C 5C 24 24 1C 06"] = "Sigma Macro 70mm f/2.8 EX DG";
        lenses["7A 3B 53 80 30 3C 4B 06"] = "Sigma 55-200mm f/4-5.6 DC HSM";
        lenses["7A 3C 1F 37 30 30 7E 06"] = "AF-S DX Zoom-Nikkor 12-24mm f/4G IF-ED";
        lenses["7A 3C 1F 37 30 30 7E 06"] = "Tokina AT-X 124 AF PRO DX II (AF 12-24mm f/4)";
        lenses["7A 3C 1F 3C 30 30 7E 06"] = "Tokina AT-X 12-28 PRO DX (AF 12-28mm f/4)";
        lenses["7A 40 2D 50 2C 3C 4B 06"] = "Sigma 18-50mm f/3.5-5.6 DC HSM";
        lenses["7A 40 2D 80 2C 40 4B 0E"] = "Sigma 18-200mm f/3.5-6.3 DC OS HSM";
        lenses["7A 47 2B 5C 24 34 4B 06"] = "Sigma 17-70mm f/2.8-4.5 DC Macro Asp. IF HSM";
        lenses["7A 47 50 76 24 24 4B 06"] = "Sigma 50-150mm f/2.8 EX APO DC HSM";
        lenses["7A 48 1C 29 24 24 7E 06"] = "Tokina AT-X 116 PRO DX II (AF 11-16mm f/2.8)";
        lenses["7A 48 1C 30 24 24 7E 06"] = "Tokina AT-X 11-20 f/2.8 PRO DX (AF 11-20mm f/2.8)";
        lenses["7A 48 2B 5C 24 34 4B 06"] = "Sigma 17-70mm f/2.8-4.5 DC Macro Asp. IF HSM";
        lenses["7A 48 2D 50 24 24 4B 06"] = "Sigma 18-50mm f/2.8 EX DC Macro";
        lenses["7A 48 5C 80 24 24 4B 06"] = "Sigma 70-200mm f/2.8 EX APO DG Macro HSM II";
        lenses["7A 54 6E 8E 24 24 4B 02"] = "Sigma APO 120-300mm f/2.8 EX DG HSM";
        lenses["7B 48 80 98 30 30 80 0E"] = "AF-S VR Zoom-Nikkor 200-400mm f/4G IF-ED";
        lenses["7D 48 2B 53 24 24 82 06"] = "AF-S DX Zoom-Nikkor 17-55mm f/2.8G IF-ED";
        lenses["7F 40 2D 5C 2C 34 84 06"] = "AF-S DX Zoom-Nikkor 18-70mm f/3.5-4.5G IF-ED";
        lenses["7F 48 2B 5C 24 34 1C 06"] = "Sigma 17-70mm f/2.8-4.5 DC Macro Asp. IF";
        lenses["7F 48 2D 50 24 24 1C 06"] = "Sigma 18-50mm f/2.8 EX DC Macro";
        lenses["80 48 1A 1A 24 24 85 06"] = "AF DX Fisheye-Nikkor 10.5mm f/2.8G ED";
        lenses["81 34 76 A6 38 40 4B 0E"] = "Sigma 150-600mm f/5-6.3 DG OS HSM | S";
        lenses["81 54 80 80 18 18 86 0E"] = "AF-S VR Nikkor 200mm f/2G IF-ED";
        lenses["82 34 76 A6 38 40 4B 0E"] = "Sigma 150-600mm f/5-6.3 DG OS HSM | C";
        lenses["82 48 8E 8E 24 24 87 0E"] = "AF-S VR Nikkor 300mm f/2.8G IF-ED";
        lenses["83 00 B0 B0 5A 5A 88 04"] = "FSA-L2, EDG 65, 800mm f/13 G";
        lenses["88 54 50 50 0C 0C 4B 06"] = "Sigma 50mm f/1.4 DG HSM | A";
        lenses["89 3C 53 80 30 3C 8B 06"] = "AF-S DX Zoom-Nikkor 55-200mm f/4-5.6G ED";
        lenses["8A 3C 37 6A 30 30 4B 0E"] = "Sigma 24-105mm f/4 DG OS HSM";
        lenses["8A 54 6A 6A 24 24 8C 0E"] = "AF-S VR Micro-Nikkor 105mm f/2.8G IF-ED";
        lenses["8B 40 2D 80 2C 3C 8D 0E"] = "AF-S DX VR Zoom-Nikkor 18-200mm f/3.5-5.6G IF-ED";
        lenses["8B 40 2D 80 2C 3C FD 0E"] = "AF-S DX VR Zoom-Nikkor 18-200mm f/3.5-5.6G IF-ED [II]";
        lenses["8B 4C 2D 44 14 14 4B 06"] = "Sigma 18-35mm f/1.8 DC HSM";
        lenses["8C 40 2D 53 2C 3C 8E 06"] = "AF-S DX Zoom-Nikkor 18-55mm f/3.5-5.6G ED";
        lenses["8D 44 5C 8E 34 3C 8F 0E"] = "AF-S VR Zoom-Nikkor 70-300mm f/4.5-5.6G IF-ED";
        lenses["8E 3C 2B 5C 24 30 4B 0E"] = "Sigma 17-70mm f/2.8-4 DC Macro OS HSM | C";
        lenses["8F 40 2D 72 2C 3C 91 06"] = "AF-S DX Zoom-Nikkor 18-135mm f/3.5-5.6G IF-ED";
        lenses["8F 48 2B 50 24 24 4B 0E"] = "Sigma 17-50mm f/2.8 EX DC OS HSM";
        lenses["90 3B 53 80 30 3C 92 0E"] = "AF-S DX VR Zoom-Nikkor 55-200mm f/4-5.6G IF-ED";
        lenses["90 40 2D 80 2C 40 4B 0E"] = "Sigma 18-200mm f/3.5-6.3 II DC OS HSM";
        lenses["91 54 44 44 0C 0C 4B 06"] = "Sigma 35mm f/1.4 DG HSM";
        lenses["92 2C 2D 88 2C 40 4B 0E"] = "Sigma 18-250mm f/3.5-6.3 DC Macro OS HSM";
        lenses["92 48 24 37 24 24 94 06"] = "AF-S Zoom-Nikkor 14-24mm f/2.8G ED";
        lenses["93 48 37 5C 24 24 95 06"] = "AF-S Zoom-Nikkor 24-70mm f/2.8G ED";
        lenses["94 40 2D 53 2C 3C 96 06"] = "AF-S DX Zoom-Nikkor 18-55mm f/3.5-5.6G ED II";
        lenses["95 00 37 37 2C 2C 97 06"] = "PC-E Nikkor 24mm f/3.5D ED";
        lenses["95 4C 37 37 2C 2C 97 02"] = "PC-E Nikkor 24mm f/3.5D ED";
        lenses["96 38 1F 37 34 3C 4B 06"] = "Sigma 12-24mm f/4.5-5.6 II DG HSM";
        lenses["96 48 98 98 24 24 98 0E"] = "AF-S VR Nikkor 400mm f/2.8G ED";
        lenses["97 3C A0 A0 30 30 99 0E"] = "AF-S VR Nikkor 500mm f/4G ED";
        lenses["97 48 6A 6A 24 24 4B 0E"] = "Sigma Macro 105mm f/2.8 EX DG OS HSM";
        lenses["98 3C A6 A6 30 30 9A 0E"] = "AF-S VR Nikkor 600mm f/4G ED";
        lenses["98 48 50 76 24 24 4B 0E"] = "Sigma 50-150mm f/2.8 EX APO DC OS HSM";
        lenses["99 40 29 62 2C 3C 9B 0E"] = "AF-S DX VR Zoom-Nikkor 16-85mm f/3.5-5.6G ED";
        lenses["99 48 76 76 24 24 4B 0E"] = "Sigma APO Macro 150mm f/2.8 EX DG OS HSM";
        lenses["9A 40 2D 53 2C 3C 9C 0E"] = "AF-S DX VR Zoom-Nikkor 18-55mm f/3.5-5.6G";
        lenses["9B 00 4C 4C 24 24 9D 06"] = "PC-E Micro Nikkor 45mm f/2.8D ED";
        lenses["9B 54 4C 4C 24 24 9D 02"] = "PC-E Micro Nikkor 45mm f/2.8D ED";
        lenses["9B 54 62 62 0C 0C 4B 06"] = "Sigma 85mm f/1.4 EX DG HSM";
        lenses["9C 48 5C 80 24 24 4B 0E"] = "Sigma 70-200mm f/2.8 EX DG OS HSM";
        lenses["9C 54 56 56 24 24 9E 06"] = "AF-S Micro Nikkor 60mm f/2.8G ED";
        lenses["9D 00 62 62 24 24 9F 06"] = "PC-E Micro Nikkor 85mm f/2.8D";
        lenses["9D 48 2B 50 24 24 4B 0E"] = "Sigma 17-50mm f/2.8 EX DC OS HSM";
        lenses["9D 54 62 62 24 24 9F 02"] = "PC-E Micro Nikkor 85mm f/2.8D";
        lenses["9E 38 11 29 34 3C 4B 06"] = "Sigma 8-16mm f/4.5-5.6 DC HSM";
        lenses["9E 40 2D 6A 2C 3C A0 0E"] = "AF-S DX VR Zoom-Nikkor 18-105mm f/3.5-5.6G ED";
        lenses["9F 37 50 A0 34 40 4B 0E"] = "Sigma 50-500mm f/4.5-6.3 DG OS HSM";
        lenses["9F 58 44 44 14 14 A1 06"] = "AF-S DX Nikkor 35mm f/1.8G";
        lenses["A0 40 2D 74 2C 3C BB 0E"] = "AF-S DX Nikkor 18-140mm f/3.5-5.6G ED VR";
        lenses["A0 48 2A 5C 24 30 4B 0E"] = "Sigma 17-70mm f/2.8-4 DC Macro OS HSM";
        lenses["A0 54 50 50 0C 0C A2 06"] = "AF-S Nikkor 50mm f/1.4G";
        lenses["A1 40 18 37 2C 34 A3 06"] = "AF-S DX Nikkor 10-24mm f/3.5-4.5G ED";
        lenses["A1 41 19 31 2C 2C 4B 06"] = "Sigma 10-20mm f/3.5 EX DC HSM";
        lenses["A1 54 55 55 0C 0C BC 06"] = "AF-S Nikkor 58mm f/1.4G";
        lenses["A2 40 2D 53 2C 3C BD 0E"] = "AF-S DX Nikkor 18-55mm f/3.5-5.6G VR II";
        lenses["A2 48 5C 80 24 24 A4 0E"] = "AF-S Nikkor 70-200mm f/2.8G ED VR II";
        lenses["A3 3C 29 44 30 30 A5 0E"] = "AF-S Nikkor 16-35mm f/4G ED VR";
        lenses["A3 3C 5C 8E 30 3C 4B 0E"] = "Sigma 70-300mm f/4-5.6 DG OS";
        lenses["A4 40 2D 8E 2C 40 BF 0E"] = "AF-S DX Nikkor 18-300mm f/3.5-6.3G ED VR";
        lenses["A4 47 2D 50 24 34 4B 0E"] = "Sigma 18-50mm f/2.8-4.5 DC OS HSM";
        lenses["A4 54 37 37 0C 0C A6 06"] = "AF-S Nikkor 24mm f/1.4G ED";
        lenses["A5 40 2D 88 2C 40 4B 0E"] = "Sigma 18-250mm f/3.5-6.3 DC OS HSM";
        lenses["A5 40 3C 8E 2C 3C A7 0E"] = "AF-S Nikkor 28-300mm f/3.5-5.6G ED VR";
        lenses["A5 4C 44 44 14 14 C0 06"] = "AF-S Nikkor 35mm f/1.8G ED";
        lenses["A6 48 37 5C 24 24 4B 06"] = "Sigma 24-70mm f/2.8 IF EX DG HSM";
        lenses["A6 48 8E 8E 24 24 A8 0E"] = "AF-S VR Nikkor 300mm f/2.8G IF-ED II";
        lenses["A6 48 98 98 24 24 C1 0E"] = "AF-S Nikkor 400mm f/2.8E FL ED VR";
        lenses["A7 3C 53 80 30 3C C2 0E"] = "AF-S DX Nikkor 55-200mm f/4-5.6G ED VR II";
        lenses["A7 49 80 A0 24 24 4B 06"] = "Sigma APO 200-500mm f/2.8 EX DG";
        lenses["A7 4B 62 62 2C 2C A9 0E"] = "AF-S DX Micro Nikkor 85mm f/3.5G ED VR";
        lenses["A8 48 80 98 30 30 AA 0E"] = "AF-S VR Zoom-Nikkor 200-400mm f/4G IF-ED II";
        lenses["A8 48 8E 8E 30 30 C3 0E"] = "AF-S Nikkor 300mm f/4E PF ED VR";
        lenses["A8 48 8E 8E 30 30 C3 4E"] = "AF-S Nikkor 300mm f/4E PF ED VR";
        lenses["A9 4C 31 31 14 14 C4 06"] = "AF-S Nikkor 20mm f/1.8G ED";
        lenses["A9 54 80 80 18 18 AB 0E"] = "AF-S Nikkor 200mm f/2G ED VR II";
        lenses["AA 3C 37 6E 30 30 AC 0E"] = "AF-S Nikkor 24-120mm f/4G ED VR";
        lenses["AA 48 37 5C 24 24 C5 4E"] = "AF-S Nikkor 24-70mm f/2.8E ED VR";
        lenses["AB 3C A0 A0 30 30 C6 4E"] = "AF-S Nikkor 500mm f/4E FL ED VR";
        lenses["AC 38 53 8E 34 3C AE 0E"] = "AF-S DX VR Nikkor 55-300mm f/4.5-5.6G ED";
        lenses["AC 3C A6 A6 30 30 C7 4E"] = "AF-S Nikkor 600mm f/4E FL ED VR";
        lenses["AD 3C 2D 8E 2C 3C AF 0E"] = "AF-S DX Nikkor 18-300mm f/3.5-5.6G ED VR";
        lenses["AD 48 28 60 24 30 C8 0E"] = "AF-S DX Nikkor 16-80mm f/2.8-4E ED VR";
        lenses["AD 48 28 60 24 30 C8 4E"] = "AF-S DX Nikkor 16-80mm f/2.8-4E ED VR";
        lenses["AE 3C 80 A0 3C 3C C9 0E"] = "AF-S Nikkor 200-500mm f/5.6E ED VR";
        lenses["AE 3C 80 A0 3C 3C C9 4E"] = "AF-S Nikkor 200-500mm f/5.6E ED VR";
        lenses["AE 54 62 62 0C 0C B0 06"] = "AF-S Nikkor 85mm f/1.4G";
        lenses["AF 4C 37 37 14 14 CC 06"] = "AF-S Nikkor 24mm f/1.8G ED";
        lenses["AF 54 44 44 0C 0C B1 06"] = "AF-S Nikkor 35mm f/1.4G";
        lenses["B0 4C 50 50 14 14 B2 06"] = "AF-S Nikkor 50mm f/1.8G";
        lenses["B1 48 48 48 24 24 B3 06"] = "AF-S DX Micro Nikkor 40mm f/2.8G";
        lenses["B2 48 5C 80 30 30 B4 0E"] = "AF-S Nikkor 70-200mm f/4G ED VR";
        lenses["B3 4C 62 62 14 14 B5 06"] = "AF-S Nikkor 85mm f/1.8G";
        lenses["B4 40 37 62 2C 34 B6 0E"] = "AF-S VR Zoom-Nikkor 24-85mm f/3.5-4.5G IF-ED";
        lenses["B5 4C 3C 3C 14 14 B7 06"] = "AF-S Nikkor 28mm f/1.8G";
        lenses["B6 3C B0 B0 3C 3C B8 0E"] = "AF-S VR Nikkor 800mm f/5.6E FL ED";
        lenses["B6 48 37 56 24 24 1C 02"] = "Sigma 24-60mm f/2.8 EX DG";
        lenses["B7 44 60 98 34 3C B9 0E"] = "AF-S Nikkor 80-400mm f/4.5-5.6G ED VR";
        lenses["B8 40 2D 44 2C 34 BA 06"] = "AF-S Nikkor 18-35mm f/3.5-4.5G ED";
        lenses["CC 4C 50 68 14 14 4B 06"] = "Sigma 50-100mm f/1.8 DC HSM | A";
        lenses["CD 3D 2D 70 2E 3C 4B 0E"] = "Sigma 18-125mm f/3.8-5.6 DC OS HSM";
        lenses["CE 34 76 A0 38 40 4B 0E"] = "Sigma 150-500mm f/5-6.3 DG OS APO HSM";
        lenses["CF 38 6E 98 34 3C 4B 0E"] = "Sigma APO 120-400mm f/4.5-5.6 DG OS HSM";
        lenses["DC 48 19 19 24 24 4B 06"] = "Sigma 10mm f/2.8 EX DC HSM Fisheye";
        lenses["DE 54 50 50 0C 0C 4B 06"] = "Sigma 50mm f/1.4 EX DG HSM";
        lenses["E0 3C 5C 8E 30 3C 4B 06"] = "Sigma 70-300mm f/4-5.6 APO DG Macro HSM";
        lenses["E1 58 37 37 14 14 1C 02"] = "Sigma 24mm f/1.8 EX DG Aspherical Macro";
        lenses["E3 54 50 50 24 24 35 02"] = "Sigma Macro 50mm f/2.8 EX DG";
        lenses["E5 54 6A 6A 24 24 35 02"] = "Sigma Macro 105mm f/2.8 EX DG";
        lenses["E6 41 3C 8E 2C 40 1C 02"] = "Sigma 28-300mm f/3.5-6.3 DG Macro";
        lenses["E8 4C 44 44 14 14 DF 0E"] = "Tamron SP 35mm f/1.8 VC";
        lenses["E9 48 27 3E 24 24 DF 0E"] = "Tamron SP 15-30mm f/2.8 Di VC USD (A012)";
        lenses["E9 54 37 5C 24 24 1C 02"] = "Sigma 24-70mm f/2.8 EX DG Macro";
        lenses["EA 40 29 8E 2C 40 DF 0E"] = "Tamron AF 16-300mm f/3.5-6.3 Di II VC PZD (B016)";
        lenses["EA 48 27 27 24 24 1C 02"] = "Sigma 15mm f/2.8 EX Diagonal Fisheye";
        lenses["EB 40 76 A6 38 40 DF 0E"] = "Tamron SP AF 150-600mm f/5-6.3 VC USD (A011)";
        lenses["ED 40 2D 80 2C 40 4B 0E"] = "Sigma 18-200mm f/3.5-6.3 DC OS HSM";
        lenses["EE 48 5C 80 24 24 4B 06"] = "Sigma 70-200mm f/2.8 EX APO DG Macro HSM II";
        lenses["F0 38 1F 37 34 3C 4B 06"] = "Sigma 12-24mm f/4.5-5.6 EX DG Aspherical HSM";
        lenses["F0 3F 2D 8A 2C 40 DF 0E"] = "Tamron AF 18-270mm f/3.5-6.3 Di II VC PZD (B008)";
        lenses["F1 44 A0 A0 34 34 4B 02"] = "Sigma APO 500mm f/4.5 EX DG HSM";
        lenses["F1 47 5C 8E 30 3C DF 0E"] = "Tamron SP 70-300mm f/4-5.6 Di VC USD (A005)";
        lenses["F3 48 68 8E 30 30 4B 02"] = "Sigma APO 100-300mm f/4 EX IF HSM";
        lenses["F3 54 2B 50 24 24 84 0E"] = "Tamron SP AF 17-50mm f/2.8 XR Di II VC LD Aspherical (IF) (B005)";
        lenses["F4 54 56 56 18 18 84 06"] = "Tamron SP AF 60mm f/2.0 Di II Macro 1:1 (G005)";
        lenses["F5 40 2C 8A 2C 40 40 0E"] = "Tamron AF 18-270mm f/3.5-6.3 Di II VC LD Aspherical (IF) Macro (B003)";
        lenses["F5 48 76 76 24 24 4B 06"] = "Sigma APO Macro 150mm f/2.8 EX DG HSM";
        lenses["F6 3F 18 37 2C 34 84 06"] = "Tamron SP AF 10-24mm f/3.5-4.5 Di II LD Aspherical (IF) (B001)";
        lenses["F6 3F 18 37 2C 34 DF 06"] = "Tamron SP AF 10-24mm f/3.5-4.5 Di II LD Aspherical (IF) (B001)";
        lenses["F6 48 2D 50 24 24 4B 06"] = "Sigma 18-50mm f/2.8 EX DC Macro";
        lenses["F7 53 5C 80 24 24 40 06"] = "Tamron SP AF 70-200mm f/2.8 Di LD (IF) Macro (A001)";
        lenses["F7 53 5C 80 24 24 84 06"] = "Tamron SP AF 70-200mm f/2.8 Di LD (IF) Macro (A001)";
        lenses["F8 54 3E 3E 0C 0C 4B 06"] = "Sigma 30mm f/1.4 EX DC HSM";
        lenses["F8 54 64 64 24 24 DF 06"] = "Tamron SP AF 90mm f/2.8 Di Macro 1:1 (272NII)";
        lenses["F8 55 64 64 24 24 84 06"] = "Tamron SP AF 90mm f/2.8 Di Macro 1:1 (272NII)";
        lenses["F9 3C 19 31 30 3C 4B 06"] = "Sigma 10-20mm f/4-5.6 EX DC HSM";
        lenses["F9 40 3C 8E 2C 40 40 0E"] = "Tamron AF 28-300mm f/3.5-6.3 XR Di VC LD Aspherical (IF) Macro (A20)";
        lenses["FA 54 3C 5E 24 24 84 06"] = "Tamron SP AF 28-75mm f/2.8 XR Di LD Aspherical (IF) Macro (A09NII)";
        lenses["FA 54 3C 5E 24 24 DF 06"] = "Tamron SP AF 28-75mm f/2.8 XR Di LD Aspherical (IF) Macro (A09NII)";
        lenses["FA 54 6E 8E 24 24 4B 02"] = "Sigma APO 120-300mm f/2.8 EX DG HSM";
        lenses["FB 54 2B 50 24 24 84 06"] = "Tamron SP AF 17-50mm f/2.8 XR Di II LD Aspherical (IF) (A16NII)";
        lenses["FB 54 8E 8E 24 24 4B 02"] = "Sigma APO 300mm f/2.8 EX DG HSM";
        lenses["FC 40 2D 80 2C 40 DF 06"] = "Tamron AF 18-200mm f/3.5-6.3 XR Di II LD Aspherical (IF) Macro (A14NII)";
        lenses["FD 47 50 76 24 24 4B 06"] = "Sigma 50-150mm f/2.8 EX APO DC HSM II";
        lenses["FE 47 00 00 24 24 4B 06"] = "Sigma 4.5mm f/2.8 EX DC HSM Circular Fisheye";
        lenses["FE 48 37 5C 24 24 DF 0E"] = "Tamron SP 24-70mm f/2.8 Di VC USD (A007)";
        lenses["FE 53 5C 80 24 24 84 06"] = "Tamron SP AF 70-200mm f/2.8 Di LD (IF) Macro (A001)";
        lenses["FE 54 5C 80 24 24 DF 0E"] = "Tamron SP 70-200mm f/2.8 Di VC USD (A009)";
        lenses["FE 54 64 64 24 24 DF 0E"] = "Tamron SP 90mm f/2.8 Di VC USD Macro 1:1 (F004)";
        lenses["FF 40 2D 80 2C 40 4B 06"] = "Sigma 18-200mm f/3.5-6.3 DC";
    }
    virtual std::string toString (Tag* t)
    {

        static const unsigned char xlat[2][256] = {
            {
                0xc1, 0xbf, 0x6d, 0x0d, 0x59, 0xc5, 0x13, 0x9d, 0x83, 0x61, 0x6b, 0x4f, 0xc7, 0x7f, 0x3d, 0x3d,
                0x53, 0x59, 0xe3, 0xc7, 0xe9, 0x2f, 0x95, 0xa7, 0x95, 0x1f, 0xdf, 0x7f, 0x2b, 0x29, 0xc7, 0x0d,
                0xdf, 0x07, 0xef, 0x71, 0x89, 0x3d, 0x13, 0x3d, 0x3b, 0x13, 0xfb, 0x0d, 0x89, 0xc1, 0x65, 0x1f,
                0xb3, 0x0d, 0x6b, 0x29, 0xe3, 0xfb, 0xef, 0xa3, 0x6b, 0x47, 0x7f, 0x95, 0x35, 0xa7, 0x47, 0x4f,
                0xc7, 0xf1, 0x59, 0x95, 0x35, 0x11, 0x29, 0x61, 0xf1, 0x3d, 0xb3, 0x2b, 0x0d, 0x43, 0x89, 0xc1,
                0x9d, 0x9d, 0x89, 0x65, 0xf1, 0xe9, 0xdf, 0xbf, 0x3d, 0x7f, 0x53, 0x97, 0xe5, 0xe9, 0x95, 0x17,
                0x1d, 0x3d, 0x8b, 0xfb, 0xc7, 0xe3, 0x67, 0xa7, 0x07, 0xf1, 0x71, 0xa7, 0x53, 0xb5, 0x29, 0x89,
                0xe5, 0x2b, 0xa7, 0x17, 0x29, 0xe9, 0x4f, 0xc5, 0x65, 0x6d, 0x6b, 0xef, 0x0d, 0x89, 0x49, 0x2f,
                0xb3, 0x43, 0x53, 0x65, 0x1d, 0x49, 0xa3, 0x13, 0x89, 0x59, 0xef, 0x6b, 0xef, 0x65, 0x1d, 0x0b,
                0x59, 0x13, 0xe3, 0x4f, 0x9d, 0xb3, 0x29, 0x43, 0x2b, 0x07, 0x1d, 0x95, 0x59, 0x59, 0x47, 0xfb,
                0xe5, 0xe9, 0x61, 0x47, 0x2f, 0x35, 0x7f, 0x17, 0x7f, 0xef, 0x7f, 0x95, 0x95, 0x71, 0xd3, 0xa3,
                0x0b, 0x71, 0xa3, 0xad, 0x0b, 0x3b, 0xb5, 0xfb, 0xa3, 0xbf, 0x4f, 0x83, 0x1d, 0xad, 0xe9, 0x2f,
                0x71, 0x65, 0xa3, 0xe5, 0x07, 0x35, 0x3d, 0x0d, 0xb5, 0xe9, 0xe5, 0x47, 0x3b, 0x9d, 0xef, 0x35,
                0xa3, 0xbf, 0xb3, 0xdf, 0x53, 0xd3, 0x97, 0x53, 0x49, 0x71, 0x07, 0x35, 0x61, 0x71, 0x2f, 0x43,
                0x2f, 0x11, 0xdf, 0x17, 0x97, 0xfb, 0x95, 0x3b, 0x7f, 0x6b, 0xd3, 0x25, 0xbf, 0xad, 0xc7, 0xc5,
                0xc5, 0xb5, 0x8b, 0xef, 0x2f, 0xd3, 0x07, 0x6b, 0x25, 0x49, 0x95, 0x25, 0x49, 0x6d, 0x71, 0xc7
            },
            {
                0xa7, 0xbc, 0xc9, 0xad, 0x91, 0xdf, 0x85, 0xe5, 0xd4, 0x78, 0xd5, 0x17, 0x46, 0x7c, 0x29, 0x4c,
                0x4d, 0x03, 0xe9, 0x25, 0x68, 0x11, 0x86, 0xb3, 0xbd, 0xf7, 0x6f, 0x61, 0x22, 0xa2, 0x26, 0x34,
                0x2a, 0xbe, 0x1e, 0x46, 0x14, 0x68, 0x9d, 0x44, 0x18, 0xc2, 0x40, 0xf4, 0x7e, 0x5f, 0x1b, 0xad,
                0x0b, 0x94, 0xb6, 0x67, 0xb4, 0x0b, 0xe1, 0xea, 0x95, 0x9c, 0x66, 0xdc, 0xe7, 0x5d, 0x6c, 0x05,
                0xda, 0xd5, 0xdf, 0x7a, 0xef, 0xf6, 0xdb, 0x1f, 0x82, 0x4c, 0xc0, 0x68, 0x47, 0xa1, 0xbd, 0xee,
                0x39, 0x50, 0x56, 0x4a, 0xdd, 0xdf, 0xa5, 0xf8, 0xc6, 0xda, 0xca, 0x90, 0xca, 0x01, 0x42, 0x9d,
                0x8b, 0x0c, 0x73, 0x43, 0x75, 0x05, 0x94, 0xde, 0x24, 0xb3, 0x80, 0x34, 0xe5, 0x2c, 0xdc, 0x9b,
                0x3f, 0xca, 0x33, 0x45, 0xd0, 0xdb, 0x5f, 0xf5, 0x52, 0xc3, 0x21, 0xda, 0xe2, 0x22, 0x72, 0x6b,
                0x3e, 0xd0, 0x5b, 0xa8, 0x87, 0x8c, 0x06, 0x5d, 0x0f, 0xdd, 0x09, 0x19, 0x93, 0xd0, 0xb9, 0xfc,
                0x8b, 0x0f, 0x84, 0x60, 0x33, 0x1c, 0x9b, 0x45, 0xf1, 0xf0, 0xa3, 0x94, 0x3a, 0x12, 0x77, 0x33,
                0x4d, 0x44, 0x78, 0x28, 0x3c, 0x9e, 0xfd, 0x65, 0x57, 0x16, 0x94, 0x6b, 0xfb, 0x59, 0xd0, 0xc8,
                0x22, 0x36, 0xdb, 0xd2, 0x63, 0x98, 0x43, 0xa1, 0x04, 0x87, 0x86, 0xf7, 0xa6, 0x26, 0xbb, 0xd6,
                0x59, 0x4d, 0xbf, 0x6a, 0x2e, 0xaa, 0x2b, 0xef, 0xe6, 0x78, 0xb6, 0x4e, 0xe0, 0x2f, 0xdc, 0x7c,
                0xbe, 0x57, 0x19, 0x32, 0x7e, 0x2a, 0xd0, 0xb8, 0xba, 0x29, 0x00, 0x3c, 0x52, 0x7d, 0xa8, 0x49,
                0x3b, 0x2d, 0xeb, 0x25, 0x49, 0xfa, 0xa3, 0xaa, 0x39, 0xa7, 0xc5, 0xa7, 0x50, 0x11, 0x36, 0xfb,
                0xc6, 0x67, 0x4a, 0xf5, 0xa5, 0x12, 0x65, 0x7e, 0xb0, 0xdf, 0xaf, 0x4e, 0xb3, 0x61, 0x7f, 0x2f
            }
        };

        int ver = (t->toInt (0, BYTE) - '0') * 1000 + (t->toInt (1, BYTE) - '0') * 100 + (t->toInt (2, BYTE) - '0') * 10 + (t->toInt (3, BYTE) - '0');

        std::ostringstream ld;
        ld << "Version = " << ver << std::endl;

        int lenstype = t->getParent()->getTag(0x0083)->toInt(0, BYTE);

        std::ostringstream lid;
        lid.setf (std::ios_base::hex, std::ios_base::basefield);
        lid.setf (std::ios_base::uppercase);

        Tag *modelTag = t->getParent()->getRoot()->findTag("Model");
        std::string model( modelTag ?  modelTag->valueToString() : "");
        int lidoffs = 7;
        bool d100 = false;

        if (model.substr(0, 10) == "NIKON D100" || model.substr(0, 9) == "NIKON D1X") {
            lidoffs = 0;
            d100 = true;
        } else if( ver < 204) {
            lidoffs = 7;
            d100 = false;
        } else {
            lidoffs = 8;
            d100 = false;
        }

        unsigned char buffer[16];

        if (d100) {
            memcpy (buffer, t->getValue() + 6, 7);
        } else {
            memcpy (buffer, t->getValue() + 4, 16);
        }

        if (ver >= 201) {
            const unsigned char* serval = t->getParent()->getTag(0x001d)->getValue ();
            int serial = 0;

            for (int i = 0; serval[i]; i++) {
                serial = serial * 10 + (isdigit(serval[i]) ? serval[i] - '0' : serval[i] % 10);
            }

            const unsigned char* scval = t->getParent()->getTag(0x00a7)->getValue ();
            int key = 0;

            for (int i = 0; i < 4; i++) {
                key ^= scval[i];
            }

            unsigned char ci = xlat[0][serial & 0xff];
            unsigned char cj = xlat[1][key];
            unsigned char ck = 0x60;

            for (int i = 0; i < 16; i++) {
                buffer[i] ^= (cj += ci * ck++);
            }
        }

        std::string EffectiveMaxApertureString = "";

        if (!d100) {
            int  EffectiveMaxApertureValue;

            if( ver < 204 ) {
                ld << "ExitPupilPosition = " << (int) buffer[0] << std::endl;
                ld << "AFAperture = "        << (int) buffer[1] << std::endl;
                ld << "FocusPosition = "     << (int) buffer[4] << std::endl;
                ld << "FocusDistance = "     << (int) buffer[5] << std::endl;
                ld << "FocalLength = "       << (int) buffer[6] << std::endl;
                EffectiveMaxApertureValue = (int) buffer[14];
            } else {
                ld << "ExitPupilPosition = " << (int) buffer[0] << std::endl;
                ld << "AFAperture = "        << (int) buffer[1] << std::endl;
                ld << "FocusPosition = "     << (int) buffer[4] << std::endl;
                ld << "FocusDistance = "     << (int) buffer[6] << std::endl;
                ld << "FocalLength = "       << (int) buffer[7] << std::endl;
                EffectiveMaxApertureValue = (int) buffer[15];
            }

            switch (EffectiveMaxApertureValue) {
            case 0x8:
                EffectiveMaxApertureString = "1.2";
                break;

            case 0xc:
                EffectiveMaxApertureString = "1.4";
                break;

            case 0x14:
                EffectiveMaxApertureString = "1.8";
                break;

            case 0x18:
                EffectiveMaxApertureString = "2.0";
                break;

            case 0x20:
                EffectiveMaxApertureString = "2.5";
                break;

            case 0x24:
                EffectiveMaxApertureString = "2.8";
                break;

            case 0x2a:
                EffectiveMaxApertureString = "3.3";
                break;

            case 0x2c:
                EffectiveMaxApertureString = "3.5";
                break;

            case 0x30:
                EffectiveMaxApertureString = "4.0";
                break;

            case 0x34:
                EffectiveMaxApertureString = "4.5";
                break;

            case 0x38:
                EffectiveMaxApertureString = "5.0";
                break;

            case 0x3c:
                EffectiveMaxApertureString = "5.6";
                break;

            case 0x40:
                EffectiveMaxApertureString = "6.3";
                break;

            case 0x44:
                EffectiveMaxApertureString = "7.1";
                break;

            case 0x48:
                EffectiveMaxApertureString = "8.0";
                break;

            case 0x4e:
                EffectiveMaxApertureString = "9.5";
                break;

            case 0x54:
                EffectiveMaxApertureString = "11.0";
                break;

            case 0x5a:
                EffectiveMaxApertureString = "13.0";
                break;

            case 0x5e:
                EffectiveMaxApertureString = "15.0";
                break;

            case 0x60:
                EffectiveMaxApertureString = "16.0";
                break;

            case 0x66:
                EffectiveMaxApertureString = "19.0";
                break;

            case 0x6c:
                EffectiveMaxApertureString = "22.0";
                break;

            default  :
                EffectiveMaxApertureString = "";
            }

            ld << "EffectiveMaxAperture = "  << EffectiveMaxApertureString << std::endl;
        }

        for (int i = 0; i < 7; i++) {
            lid << std::setw(2) << std::setfill('0') << (int)buffer[lidoffs + i] << ' ';
        }

        lid << std::setw(2) << std::setfill('0') << lenstype;

        std::map<std::string, std::string>::iterator r = lenses.find (lid.str());

        if (r != lenses.end()) {
            if(r == lenses.begin() && EffectiveMaxApertureString != "") {       // first entry is for unchipped lenses
                ld << "Lens = Unknown $FL$mm f/" << EffectiveMaxApertureString;
            } else {
                ld << "Lens = " << r->second;
            }
        } else {
            ld << "Lens = Unknown, ID=" << lid.str();
        }

        return ld.str();
    }

};
NALensDataInterpreter naLensDataInterpreter;

const TagAttrib nikonISOInfoAttribs[] = {
    {0, AC_WRITE, 0, 0, 0x0000, AUTO, "ISO", &naISOInfoISOInterpreter},
    {0, AC_WRITE, 0, 0, 0x0004, SHORT, "ISOExpansion", &naISOExpansionInterpreter},
    {0, AC_WRITE, 0, 0, 0x0006, AUTO, "ISO2", &naISOInfoISOInterpreter},
    {0, AC_WRITE, 0, 0, 0x000a, SHORT, "ISOExpansion2", &naISOExpansionInterpreter},
    { -1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}
};

const TagAttrib nikon2Attribs[] = {
    {0, AC_WRITE, 0, 0, 0x0002, AUTO, "Unknown", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0003, AUTO, "Quality", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0004, AUTO, "ColorMode", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0005, AUTO, "ImageAdjustment", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0006, AUTO, "ISOSpeed", &naISOInterpreter},
    {0, AC_WRITE, 0, 0, 0x0007, AUTO, "WhiteBalance", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0008, AUTO, "Focus", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0009, AUTO, "Unknown", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x000a, AUTO, "DigitalZoom", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x000b, AUTO, "AuxiliaryLens", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0f00, AUTO, "Unknown", &stdInterpreter},
    { -1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}
};

const TagAttrib nikon3Attribs[] = {
    {0, AC_WRITE, 0, 0, 0x0001, AUTO, "MakerNoteVersion", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0002, AUTO, "ISOSpeed", &naISOInterpreter},
    {0, AC_WRITE, 0, 0, 0x0003, AUTO, "ColorMode", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0004, AUTO, "Quality", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0005, AUTO, "WhiteBalance", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0006, AUTO, "Sharpness", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0007, AUTO, "FocusMode", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0008, AUTO, "FlashSetting", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0009, AUTO, "FlashType", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x000b, AUTO, "WhiteBalanceFineTune", &stdInterpreter},
    {0, AC_NEW,   0, 0, 0x000c, AUTO, "ColorBalance1", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x000d, AUTO, "ProgramShift", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x000e, AUTO, "ExposureDifference", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x000f, AUTO, "ISOSelection", &naISOInterpreter},
    {0, AC_WRITE, 0, 0, 0x0010, AUTO, "DataDump", &stdInterpreter},
    {1, AC_WRITE, 0, 0, 0x0011, AUTO, "NikonPreview", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0012, AUTO, "FlashExposureComp", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0013, AUTO, "ISOSetting", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0016, AUTO, "ImageBoundary", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0018, AUTO, "FlashExposureBracketValue", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0019, AUTO, "ExposureBracketValue", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x001a, AUTO, "ImageProcessing", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x001b, AUTO, "CropHiSpeed", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x001d, AUTO, "SerialNumber", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x001e, AUTO, "ColorSpace", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0020, AUTO, "ImageAuthentication", &stdInterpreter},
    {0, AC_WRITE, 0, nikonISOInfoAttribs, 0x0025, AUTO, "ISOInfo", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0080, AUTO, "ImageAdjustment", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0081, AUTO, "ToneComp", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0082, AUTO, "AuxiliaryLens", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0083, AUTO, "LensType", &naLensTypeInterpreter},
    {0, AC_WRITE, 0, 0, 0x0084, AUTO, "Lens", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0085, AUTO, "ManualFocusDistance", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0086, AUTO, "DigitalZoom", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0087, AUTO, "FlashMode", &naFlashModeInterpreter},
    {0, AC_WRITE, 0, 0, 0x0088, AUTO, "AFInfo", &naAFInfoInterpreter},
    {0, AC_WRITE, 0, 0, 0x0089, AUTO, "ShootingMode", &naShootingModeInterpreter},
    {0, AC_WRITE, 0, 0, 0x008a, AUTO, "AutoBracketRelease", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x008b, AUTO, "LensFStops", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x008c, AUTO, "NEFCurve1", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x008d, AUTO, "ColorHue", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x008f, AUTO, "SceneMode", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0090, AUTO, "LightSource", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0091, AUTO, "ShotInfo", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0092, AUTO, "HueAdjustment", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0094, AUTO, "Saturation", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0095, AUTO, "NoiseReduction", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0096, AUTO, "NEFCurve2", &stdInterpreter},
    {0, AC_NEW,   0, 0, 0x0097, AUTO, "ColorBalance", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x0098, AUTO, "LensData", &naLensDataInterpreter},
    {0, AC_WRITE, 0, 0, 0x0099, AUTO, "RawImageCenter", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x009a, AUTO, "SensorPixelSize", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x00a0, AUTO, "SerialNumber", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x00a2, AUTO, "ImageDataSize", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x00a5, AUTO, "ImageCount", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x00a6, AUTO, "DeletedImageCount", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x00a7, AUTO, "ShutterCount", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x00a9, AUTO, "ImageOptimization", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x00aa, AUTO, "Saturation", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x00ab, AUTO, "VariProgram", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x00ac, AUTO, "ImageStabilization", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x00ad, AUTO, "AFResponse", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x00b0, AUTO, "MultiExposure", &stdInterpreter},
    {0, AC_WRITE, 0, 0, 0x00b1, AUTO, "HighISONoiseReduction", &naHiISONRInterpreter},
    {0, AC_WRITE, 0, 0, 0x0e00, AUTO, "PrintIM", &stdInterpreter},
    {0, AC_DONTWRITE, 0, 0, 0x0e01, AUTO, "NikonCaptureData", &stdInterpreter},
    {0, AC_DONTWRITE, 0, 0, 0x0e09, AUTO, "NikonCaptureVersion", &stdInterpreter},
    {0, AC_DONTWRITE, 0, 0, 0x0e0e, AUTO, "NikonCaptureOffsets", &stdInterpreter},
    {0, AC_DONTWRITE, 0, 0, 0x0e10, AUTO, "NikonScanIFD", &stdInterpreter},
    { -1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}
};

}
#endif

