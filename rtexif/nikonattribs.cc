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

#include <rtexif.h>
#include <string>
#include <map>
#include <math.h>
#include <sstream>
#include <iomanip>
#include <string.h>

namespace rtexif {

class NAISOInterpreter : public Interpreter {
    public:
        NAISOInterpreter () {}
        virtual std::string toString (Tag* t) {
            sprintf (buffer, "%d", t->toInt(2));
            return buffer;
        }
};
NAISOInterpreter naISOInterpreter;

class NALensTypeInterpreter : public Interpreter {
    public:
        NALensTypeInterpreter () {}
        virtual std::string toString (Tag* t) {
            int a = t->toInt();
            std::ostringstream str;
            str << "MF = " << (a&1 ? "Yes" : "No") << std::endl;
            str << "D = " << (a&2 ? "Yes" : "No") << std::endl;
            str << "G = " << (a&4 ? "Yes" : "No") << std::endl;
            str << "VR = " << (a&8 ? "Yes" : "No");
            return str.str();
        }
};
NALensTypeInterpreter naLensTypeInterpreter;

class NAFlashModeInterpreter : public ChoiceInterpreter {
    public:
        NAFlashModeInterpreter () {
            choices[0]      = "Did Not Fire";
            choices[1]      = "Fired, Manual";
            choices[7]      = "Fired, External";
            choices[8]      = "Fired, Commander Mode";
            choices[9]      = "Fired, TTL Mode";
        }
};
NAFlashModeInterpreter naFlashModeInterpreter;

class NAHiISONRInterpreter : public ChoiceInterpreter {
    public:
        NAHiISONRInterpreter () {
            choices[0]      = "Off";
            choices[1]      = "Minimal";
            choices[2]      = "Low";
            choices[4]      = "Normal";
            choices[6]      = "High";
        }
};
NAHiISONRInterpreter naHiISONRInterpreter;

class NAShootingModeInterpreter : public Interpreter {
    public:
        NAShootingModeInterpreter () {}
        virtual std::string toString (Tag* t) {
            int a = t->toInt();
            std::ostringstream str;
            str << "Continuous = " << (a&1 ? "Yes" : "No") << std::endl;
            str << "Delay = " << (a&2 ? "Yes" : "No") << std::endl;
            str << "PC Control = " << (a&4 ? "Yes" : "No") << std::endl;
            str << "Exposure Bracketing = " << (a&8 ? "Yes" : "No") << std::endl;
            str << "Auto ISO = " << (a&16 ? "Yes" : "No") << std::endl;
            str << "White-Balance Bracketing = " << (a&32 ? "Yes" : "No") << std::endl;
            str << "IR Control = " << (a&64 ? "Yes" : "No");
            return str.str();
        }
};
NAShootingModeInterpreter naShootingModeInterpreter;

class NAAFInfoInterpreter : public Interpreter {
        std::map<int,std::string> amchoices;
        std::map<int,std::string> afpchoices;
    public:
        NAAFInfoInterpreter () {
            amchoices[0] = "Single Area";
            amchoices[1] = "Dynamic Area";
            amchoices[2] = "Dynamic Area, Closest Subject";
            amchoices[3] = "Group Dynamic";
            amchoices[4] = "Single Area (wide)";
            amchoices[5] = "Dynamic Area (wide)";
            afpchoices[0] = "Center";
            afpchoices[1] = "Top";
            afpchoices[2] = "Bottom";
            afpchoices[3] = "Left";
            afpchoices[4] = "Right";
            afpchoices[5] = "Upper-left";
            afpchoices[6] = "Upper-right";
            afpchoices[7] = "Lower-left";
            afpchoices[8] = "Lower-right";
            afpchoices[9] = "Far Left";
            afpchoices[10] = "Far Right";
        }
        virtual std::string toString (Tag* t) {
            int am = t->toInt (0, BYTE);
            int afp = t->toInt (1, BYTE);
            int aff = t->toInt (2, SHORT);
            std::ostringstream str;
            str << "AFAreaMode = " << amchoices[am] << std::endl;
            str << "AFAreaMode = " << afpchoices[afp] << std::endl;

            std::ostringstream af;
            if (aff&1)
                if (af.str()=="") af << "Center";
                else af << ", Center";
            else if (aff&2)
                if (af.str()=="") af << "Top";
                else af << ", Top";
            else if (aff&4)
                if (af.str()=="") af << "Bottom";
                else af << ", Bottom";
            else if (aff&8)
                if (af.str()=="") af << "Left";
                else af << ", Left";
            else if (aff&16)
                if (af.str()=="") af << "Right";
                else af << ", Right";
            else if (aff&32)
                if (af.str()=="") af << "Upper-left";
                else af << ", Upper-left";
            else if (aff&64)
                if (af.str()=="") af << "Upper-right";
                else af << ", Upper-right";
            else if (aff&128)
                if (af.str()=="") af << " Lower-left";
                else af << ",  Lower-left";
            else if (aff&256)
                if (af.str()=="") af << "Lower-right";
                else af << ", Lower-right";
            else if (aff&512)
                if (af.str()=="") af << "Far Left";
                else af << ", Far Left";
            else if (aff&1024)
                if (af.str()=="") af << "Far Right";
                else af << ", Far Right";

            str << "AFPointsInFocus = " << af.str();
            return str.str();
        }
};
NAAFInfoInterpreter naAFInfoInterpreter;

class NALensDataInterpreter : public Interpreter {
        std::map<std::string,std::string> lenses;
    public:
        NALensDataInterpreter () {
            lenses["00 00 00 00 00 00 00 01"] = "Manual Lens No CPU ";
            lenses["00 00 00 00 00 00 E1 12"] = "TC-17E II ";
            lenses["00 00 00 00 00 00 F1 0C"] = "TC-14E [II] or Sigma APO Tele Converter 1.4x EX DG or Kenko Teleplus PRO 300 DG 1.4x";
            lenses["00 00 00 00 00 00 F2 18"] = "TC-20E [II] or Sigma APO Tele Converter 2x EX DG or Kenko Teleplus PRO 300 DG 2.0x";
            lenses["00 36 1C 2D 34 3C 00 06"] = "Tamron SP AF11-18mm f/4.5-5.6 Di II LD Aspherical (IF)";
            lenses["00 3C 1F 37 30 30 00 06"] = "Tokina AT-X 124 AF PRO DX - AF 12-24mm F4";
            lenses["00 3E 80 A0 38 3F 00 02"] = "Tamron SP AF200-500mm f/5-6.3 Di LD (IF)";
            lenses["00 3F 2D 80 2B 40 00 06"] = "Tamron AF18-200mm f/3.5-6.3 XR Di II LD Aspherical (IF)";
            lenses["00 3F 2D 80 2C 40 00 06"] = "Tamron AF18-200mm f/3.5-6.3 XR Di II LD Aspherical (IF) Macro";
            lenses["00 3F 80 A0 38 3F 00 02"] = "Tamron SP AF200-500mm f/5-6.3 Di";
            lenses["00 40 18 2B 2C 34 00 06"] = "Tokina AT-X 107 DX Fish-Eye - AF 10-17mm F3.5-4.5";
            lenses["00 40 2A 72 2C 3C 00 06"] = "Tokina AT-X 16.5-135 DX (AF 16.5-135mm F3.5-5.6)";
            lenses["00 40 2B 2B 2C 2C 00 02"] = "Tokina AT-X 17 AF PRO - AF 17mm F3.5";
            lenses["00 40 2D 80 2C 40 00 06"] = "Tamron AF 18-200mm f/3.5-6.3 XR Di II LD Aspherical (IF) Macro (A14NII)";
            lenses["00 40 2D 88 2C 40 00 06"] = "Tamron AF 18-250mm f/3.5-6.3 Di II LD Aspherical (IF) Macro (A18NII)";
            lenses["00 40 2D 88 2C 40 62 06"] = "Tamron AF 18-250mm f/3.5-6.3 Di II LD Aspherical (IF) Macro (A18)";
            lenses["00 40 31 31 2C 2C 00 00"] = "Voigtlander Color Skopar 20mm F3.5 SLII Aspherical";
            lenses["00 44 60 98 34 3C 00 02"] = "Tokina AT-X 840D 80-400mm F4.5-5.6";
            lenses["00 47 10 10 24 24 00 00"] = "Fisheye Nikkor 8mm f/2.8 AiS";
            lenses["00 47 44 44 24 24 00 06"] = "Tokina AT-X M35 PRO DX (AF 35mm f/2.8 Macro)";
            lenses["00 47 53 80 30 3C 00 06"] = "Tamron AF55-200mm f/4-5.6 Di II LD";
            lenses["00 48 1C 29 24 24 00 06"] = "Tokina AT-X 116 PRO DX (AF 11-16mm f/2.8)";
            lenses["00 48 29 50 24 24 00 06"] = "Tokina AT-X 165 PRO DX - AF 16-50mm F2.8";
            lenses["00 48 3C 60 24 24 00 02"] = "Tokina AT-X 280 AF PRO 28-80mm F2.8 Aspherical";
            lenses["00 48 3C 6A 24 24 00 02"] = "Tamron SP AF28-105mm f/2.8";
            lenses["00 48 50 50 18 18 00 00"] = "Nikkor H 50mm f/2";
            lenses["00 48 50 72 24 24 00 06"] = "Tokina AT-X 535 PRO DX - AF 50-135mm F2.8";
            lenses["00 48 5C 8E 30 3C 00 06"] = "Tamron AF 70-300mm f/4-5.6 Di LD Macro 1:2 (A17)";
            lenses["00 48 68 68 24 24 00 00"] = "Series E 100mm f/2.8";
            lenses["00 48 80 80 30 30 00 00"] = "Nikkor 200mm f/4 AiS";
            lenses["00 49 30 48 22 2B 00 02"] = "Tamron SP AF20-40mm f/2.7-3.5";
            lenses["00 4C 6A 6A 20 20 00 00"] = "Nikkor 105mm f/2.5 AiS";
            lenses["00 4C 7C 7C 2C 2C 00 02"] = "Tamron SP AF180mm f/3.5 Di Model B01";
            lenses["00 53 2B 50 24 24 00 06"] = "Tamron SP AF17-50mm f/2.8 (A16)";
            lenses["00 54 2B 50 24 24 00 06"] = "Tamron SP AF 17-50mm f/2.8 XR Di II LD Aspherical (IF) (A16NII)";
            lenses["00 54 44 44 0C 0C 00 00"] = "Nikkor 35mm f/1.4 AiS";
            lenses["00 54 48 48 18 18 00 00"] = "Voigtlander Ultron 40mm F2 SLII Aspherical";
            lenses["00 54 55 55 0C 0C 00 00"] = "Voigtlander Nokton 58mm F1.4 SLII";
            lenses["00 54 56 56 30 30 00 00"] = "Coastal Optical Systems 60mm 1:4 UV-VIS-IR Macro Apo";
            lenses["00 54 68 68 24 24 00 02"] = "Tokina AT-X M100 PRO D - 100mm F2.8";
            lenses["00 54 8E 8E 24 24 00 02"] = "Tokina AT-X 300 AF PRO 300mm F2.8";
            lenses["01 00 00 00 00 00 02 00"] = "AF Teleconverter TC-16A 1.6x";
            lenses["01 00 00 00 00 00 08 00"] = "AF Teleconverter TC-16A 1.6x";
            lenses["01 58 50 50 14 14 02 00"] = "AF Nikkor 50mm f/1.8";
            lenses["02 2F 98 98 3D 3D 02 00"] = "Sigma 400mm F5.6 APO";
            lenses["02 37 5E 8E 35 3D 02 00"] = "Sigma 75-300mm F4.5-5.6 APO";
            lenses["02 37 A0 A0 34 34 02 00"] = "Sigma APO 500mm F4.5";
            lenses["02 3A 5E 8E 32 3D 02 00"] = "Sigma 75-300mm F4.0-5.6";
            lenses["02 3B 44 61 30 3D 02 00"] = "Sigma 35-80mm F4-5.6";
            lenses["02 3F 24 24 2C 2C 02 00"] = "Sigma 14mm F3.5";
            lenses["02 3F 3C 5C 2D 35 02 00"] = "Sigma 28-70mm F3.5-4.5 UC";
            lenses["02 40 44 73 2B 36 02 00"] = "Sigma 35-135mm F3.5-4.5 a";
            lenses["02 42 44 5C 2A 34 02 00"] = "AF Zoom-Nikkor 35-70mm f/3.3-4.5";
            lenses["02 42 44 5C 2A 34 08 00"] = "AF Zoom-Nikkor 35-70mm f/3.3-4.5";
            lenses["02 46 37 37 25 25 02 00"] = "Sigma 24mm F2.8 Macro";
            lenses["02 46 3C 5C 25 25 02 00"] = "Sigma 28-70mm F2.8";
            lenses["02 46 5C 82 25 25 02 00"] = "Sigma 70-210mm F2.8 APO";
            lenses["02 48 65 65 24 24 02 00"] = "Sigma 90mm F2.8 Macro";
            lenses["03 43 5C 81 35 35 02 00"] = "Soligor AF C/D ZOOM UMCS 70-210mm 1:4.5";
            lenses["03 48 5C 81 30 30 02 00"] = "AF Zoom-Nikkor 70-210mm f/4";
            lenses["04 48 3C 3C 24 24 03 00"] = "AF Nikkor 28mm f/2.8";
            lenses["05 54 50 50 0C 0C 04 00"] = "AF Nikkor 50mm f/1.4";
            lenses["06 3F 68 68 2C 2C 06 00"] = "Cosina 100mm f/3.5 Macro";
            lenses["06 54 53 53 24 24 06 00"] = "AF Micro-Nikkor 55mm f/2.8";
            lenses["07 36 3D 5F 2C 3C 03 00"] = "Cosina AF Zoom 28-80mm F3.5-5.6 MC Macro";
            lenses["07 3E 30 43 2D 35 03 00"] = "Soligor AF Zoom 19-35mm 1:3.5-4.5 MC";
            lenses["07 40 2F 44 2C 34 03 02"] = "Tamron AF19-35mm f/3.5-4.5 N";
            lenses["07 40 30 45 2D 35 03 02"] = "Tamron AF19-35mm f/3.5-4.5";
            lenses["07 40 3C 62 2C 34 03 00"] = "AF Zoom-Nikkor 28-85mm f/3.5-4.5";
            lenses["07 46 2B 44 24 30 03 02"] = "Tamron SP AF17-35mm f/2.8-4 Di LD Aspherical (IF)";
            lenses["07 46 3D 6A 25 2F 03 00"] = "Cosina AF Zoom 28-105mm F2.8-3.8 MC";
            lenses["07 47 3C 5C 25 35 03 00"] = "Tokina AF 287 SD (AF 28-70mm f/2.8-4.5)";
            lenses["07 48 3C 5C 24 24 03 00"] = "Tokina AT-X 287 AF (AF 28-70mm f/2.8)";
            lenses["08 40 44 6A 2C 34 04 00"] = "AF Zoom-Nikkor 35-105mm f/3.5-4.5";
            lenses["09 48 37 37 24 24 04 00"] = "AF Nikkor 24mm f/2.8";
            lenses["0A 48 8E 8E 24 24 03 00"] = "AF Nikkor 300mm f/2.8 IF-ED";
            lenses["0B 3E 3D 7F 2F 3D 0E 00"] = "Tamron AF28-200mm f/3.8-5.6";
            lenses["0B 3E 3D 7F 2F 3D 0E 02"] = "Tamron AF28-200mm f/3.8-5.6D";
            lenses["0B 48 7C 7C 24 24 05 00"] = "AF Nikkor 180mm f/2.8 IF-ED";
            lenses["0D 40 44 72 2C 34 07 00"] = "AF Zoom-Nikkor 35-135mm f/3.5-4.5";
            lenses["0E 48 5C 81 30 30 05 00"] = "AF Zoom-Nikkor 70-210mm f/4";
            lenses["0E 4A 31 48 23 2D 0E 02"] = "Tamron SP AF20-40mm f/2.7-3.5";
            lenses["0F 58 50 50 14 14 05 00"] = "AF Nikkor 50mm f/1.8 N";
            lenses["10 3D 3C 60 2C 3C D2 02"] = "Tamron AF28-80mm f/3.5-5.6 Aspherical";
            lenses["10 48 8E 8E 30 30 08 00"] = "AF Nikkor 300mm f/4 IF-ED";
            lenses["11 48 44 5C 24 24 08 00"] = "AF Zoom-Nikkor 35-70mm f/2.8";
            lenses["12 36 5C 81 35 3D 09 00"] = "Cosina AF Zoom 70-210mm F4.5-5.6 MC Macro";
            lenses["12 39 5C 8E 34 3D 08 02"] = "Cosina AF Zoom 70-300mm F4.5-5.6 MC Macro";
            lenses["12 3B 68 8D 3D 43 09 02"] = "Unknown 100-290mm f/5.6-6.7";
            lenses["12 3D 3C 80 2E 3C DF 02"] = "Tamron AF 28-200mm f/3.8-5.6 AF Aspherical LD (IF) (271D)";
            lenses["12 48 5C 81 30 3C 09 00"] = "AF Nikkor 70-210mm f/4-5.6";
            lenses["12 4A 5C 81 31 3D 09 00"] = "Soligor AF C/D Auto Zoom+Macro 70-210mm 1:4-5.6 UMCS";
            lenses["13 42 37 50 2A 34 0B 00"] = "AF Zoom-Nikkor 24-50mm f/3.3-4.5";
            lenses["14 48 60 80 24 24 0B 00"] = "AF Zoom-Nikkor 80-200mm f/2.8 ED";
            lenses["14 48 68 8E 30 30 0B 00"] = "Tokina AT-X 340 AF II 100-300mm F4";
            lenses["14 54 60 80 24 24 0B 00"] = "Tokina AT-X 828 AF 80-200mm F2.8";
            lenses["15 4C 62 62 14 14 0C 00"] = "AF Nikkor 85mm f/1.8";
            lenses["17 3C A0 A0 30 30 0F 00"] = "Nikkor 500mm f/4 P ED IF";
            lenses["17 3C A0 A0 30 30 11 00"] = "Nikkor 500mm f/4 P ED IF";
            lenses["18 40 44 72 2C 34 0E 00"] = "AF Zoom-Nikkor 35-135mm f/3.5-4.5 N";
            lenses["1A 54 44 44 18 18 11 00"] = "AF Nikkor 35mm f/2";
            lenses["1B 44 5E 8E 34 3C 10 00"] = "AF Zoom-Nikkor 75-300mm f/4.5-5.6";
            lenses["1C 48 30 30 24 24 12 00"] = "AF Nikkor 20mm f/2.8";
            lenses["1D 42 44 5C 2A 34 12 00"] = "AF Zoom-Nikkor 35-70mm f/3.3-4.5 N";
            lenses["1E 54 56 56 24 24 13 00"] = "AF Micro-Nikkor 60mm f/2.8";
            lenses["1E 5D 64 64 20 20 13 00"] = "Unknown 90mm f/2.5";
            lenses["1F 54 6A 6A 24 24 14 00"] = "AF Micro-Nikkor 105mm f/2.8";
            lenses["20 3C 80 98 3D 3D 1E 02"] = "Tamron AF200-400mm f/5.6 LD IF";
            lenses["20 48 60 80 24 24 15 00"] = "AF Zoom-Nikkor ED 80-200mm f/2.8";
            lenses["21 40 3C 5C 2C 34 16 00"] = "AF Zoom-Nikkor 28-70mm f/3.5-4.5";
            lenses["22 48 72 72 18 18 16 00"] = "AF DC-Nikkor 135mm f/2";
            lenses["23 30 BE CA 3C 48 17 00"] = "Zoom-Nikkor 1200-1700mm f/5.6-8 P ED IF";
            lenses["24 44 60 98 34 3C 1A 02"] = "Tokina AT-X 840 AF II 80-400mm F4.5-5.6";
            lenses["24 48 60 80 24 24 1A 02"] = "AF Zoom-Nikkor ED 80-200mm f/2.8D";
            lenses["25 48 3C 5C 24 24 1B 02"] = "Tokina AT-X 287 AF PRO SV 28-70mm F2.8";
            lenses["25 48 44 5C 24 24 1B 02"] = "AF Zoom-Nikkor 35-70mm f/2.8D";
            lenses["25 48 44 5C 24 24 52 02"] = "AF Zoom-Nikkor 35-70mm f/2.8D";
            lenses["26 3C 54 80 30 3C 1C 06"] = "Sigma 55-200mm F4-5.6 DC";
            lenses["26 3C 5C 82 30 3C 1C 02"] = "Sigma 70-210mm F4-5.6 UC-II";
            lenses["26 3C 5C 8E 30 3C 1C 02"] = "Sigma 70-300mm F4-5.6 DG Macro";
            lenses["26 3D 3C 80 2F 3D 1C 02"] = "Sigma 28-300mm F3.8-5.6 Aspherical";
            lenses["26 3E 3C 6A 2E 3C 1C 02"] = "Sigma 28-105mm F3.8-5.6 UC-III Aspherical IF";
            lenses["26 40 27 3F 2C 34 1C 02"] = "Sigma 15-30mm F3.5-4.5 EX Aspherical DG DF";
            lenses["26 40 2D 44 2B 34 1C 02"] = "Sigma 18-35 F3.5-4.5 Aspherical";
            lenses["26 40 2D 50 2C 3C 1C 06"] = "Sigma 18-50mm F3.5-5.6 DC";
            lenses["26 40 2D 70 2B 3C 1C 06"] = "Sigma 18-125mm F3.5-5.6 DC";
            lenses["26 40 2D 80 2C 40 1C 06"] = "Sigma 18-200mm F3.5-6.3 DC";
            lenses["26 40 37 5C 2C 3C 1C 02"] = "Sigma 24-70mm F3.5-5.6 Aspherical HF";
            lenses["26 40 3C 60 2C 3C 1C 02"] = "Sigma 28-80mm F3.5-5.6 Mini Zoom Macro II Aspherical";
            lenses["26 40 3C 65 2C 3C 1C 02"] = "Sigma 28-90mm F3.5-5.6 Macro";
            lenses["26 40 3C 80 2B 3C 1C 02"] = "Sigma 28-200mm F3.5-5.6 Compact Aspherical Hyperzoom Macro";
            lenses["26 40 3C 80 2C 3C 1C 02"] = "Sigma 28-200mm F3.5-5.6 Compact Aspherical Hyperzoom Macro";
            lenses["26 40 3C 8E 2C 40 1C 02"] = "Sigma 28-300mm F3.5-6.3 Macro";
            lenses["26 40 7B A0 34 40 1C 02"] = "Sigma APO 170-500mm F5-6.3 Aspherical RF";
            lenses["26 41 3C 8E 2C 40 1C 02"] = "Sigma 28-300mm F3.5-6.3 DG Macro";
            lenses["26 44 73 98 34 3C 1C 02"] = "Sigma 135-400mm F4.5-5.6 APO Aspherical";
            lenses["26 48 11 11 30 30 1C 02"] = "Sigma 8mm F4 EX Circular Fisheye";
            lenses["26 48 27 27 24 24 1C 02"] = "Sigma 15mm F2.8 EX Diagonal Fish-Eye";
            lenses["26 48 2D 50 24 24 1C 06"] = "Sigma 18-50mm F2.8 EX DC";
            lenses["26 48 31 49 24 24 1C 02"] = "Sigma 20-40mm F2.8";
            lenses["26 48 37 56 24 24 1C 02"] = "Sigma 24-60mm F2.8 EX DG";
            lenses["26 48 3C 5C 24 24 1C 06"] = "Sigma 28-70mm F2.8 EX DG";
            lenses["26 48 3C 5C 24 30 1C 02"] = "Sigma 28-70mm F2.8-4 High Speed Zoom";
            lenses["26 48 3C 6A 24 30 1C 02"] = "Sigma 28-105mm F2.8-4 Aspherical";
            lenses["26 48 8E 8E 30 30 1C 02"] = "Sigma APO TELE MACRO 300mm F4";
            lenses["26 54 2B 44 24 30 1C 02"] = "Sigma 17-35mm F2.8-4 EX Aspherical";
            lenses["26 54 37 5C 24 24 1C 02"] = "Sigma 24-70mm F2.8 EX DG Macro";
            lenses["26 54 37 73 24 34 1C 02"] = "Sigma 24-135mm F2.8-4.5";
            lenses["26 54 3C 5C 24 24 1C 02"] = "Sigma 28-70mm F2.8 EX";
            lenses["26 58 31 31 14 14 1C 02"] = "Sigma 20mm F1.8 EX Aspherical DG DF RF";
            lenses["26 58 37 37 14 14 1C 02"] = "Sigma 24mm F1.8 EX Aspherical DG DF MACRO";
            lenses["26 58 3C 3C 14 14 1C 02"] = "Sigma 28mm F1.8 EX DG DF";
            lenses["27 48 8E 8E 24 24 1D 02"] = "AF-I Nikkor 300mm f/2.8D IF-ED";
            lenses["27 48 8E 8E 24 24 E1 02"] = "AF-I Nikkor 300mm f/2.8D IF-ED + TC-17E";
            lenses["27 48 8E 8E 24 24 F1 02"] = "AF-I Nikkor 300mm f/2.8D IF-ED + TC-14E";
            lenses["27 48 8E 8E 24 24 F2 02"] = "AF-I Nikkor 300mm f/2.8D IF-ED + TC-20E";
            lenses["28 3C A6 A6 30 30 1D 02"] = "AF-I Nikkor 600mm f/4D IF-ED";
            lenses["28 3C A6 A6 30 30 E1 02"] = "AF-I Nikkor 600mm f/4D IF-ED + TC-17E";
            lenses["28 3C A6 A6 30 30 F1 02"] = "AF-I Nikkor 600mm f/4D IF-ED + TC-14E";
            lenses["28 3C A6 A6 30 30 F2 02"] = "AF-I Nikkor 600mm f/4D IF-ED + TC-20E";
            lenses["2A 54 3C 3C 0C 0C 26 02"] = "AF Nikkor 28mm f/1.4D";
            lenses["2B 3C 44 60 30 3C 1F 02"] = "AF Zoom-Nikkor 35-80mm f/4-5.6D";
            lenses["2C 48 6A 6A 18 18 27 02"] = "AF DC-Nikkor 105mm f/2D";
            lenses["2D 48 80 80 30 30 21 02"] = "AF Micro-Nikkor 200mm f/4D IF-ED";
            lenses["2E 48 5C 82 30 3C 28 02"] = "AF Nikkor 70-210mm f/4-5.6D";
            lenses["2F 40 30 44 2C 34 29 02"] = "Unknown 20-35mm f/3.5-4.5D";
            lenses["2F 48 30 44 24 24 29 02"] = "Tokina AT-X 235 AF PRO - AF 20-35mm f/2.8";
            lenses["30 48 98 98 24 24 24 02"] = "AF-I Nikkor 400mm f/2.8D IF-ED";
            lenses["30 48 98 98 24 24 E1 02"] = "AF-I Nikkor 400mm f/2.8D IF-ED + TC-17E";
            lenses["30 48 98 98 24 24 F1 02"] = "AF-I Nikkor 400mm f/2.8D IF-ED + TC-14E";
            lenses["30 48 98 98 24 24 F2 02"] = "AF-I Nikkor 400mm f/2.8D IF-ED + TC-20E";
            lenses["31 54 56 56 24 24 25 02"] = "AF Micro-Nikkor 60mm f/2.8D";
            lenses["32 53 64 64 24 24 35 02"] = "Tamron SP AF90mm f/2.8 Di Macro 1:2 (272E)";
            lenses["32 54 50 50 24 24 35 02"] = "Sigma 50mm F2.8 EX DG Macro";
            lenses["32 54 6A 6A 24 24 35 02"] = "AF Micro-Nikkor 105mm f/2.8D";
            lenses["33 48 2D 2D 24 24 31 02"] = "AF Nikkor 18mm f/2.8D";
            lenses["33 54 3C 5E 24 24 62 02"] = "Tamron SP AF28-75mm f/2.8 XR Di LD Aspherical (IF) Macro";
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
            lenses["44 44 60 80 34 3C 47 02"] = "AF Zoom-Nikkor 80-200mm f/4.5-5.6D ";
            lenses["45 3D 3C 60 2C 3C 48 02"] = "Tamron AF28-80mm f/3.5-5.6 Aspherical";
            lenses["45 40 3C 60 2C 3C 48 02"] = "AF Zoom-Nikkor 28-80mm F/3.5-5.6D";
            lenses["45 41 37 72 2C 3C 48 02"] = "Tamron SP AF24-135mm f/3.5-5.6 AD Aspherical (IF) Macro";
            lenses["46 3C 44 60 30 3C 49 02"] = "AF Zoom-Nikkor 35-80mm f/4-5.6D N";
            lenses["47 42 37 50 2A 34 4A 02"] = "AF Zoom-Nikkor 24-50mm f/3.3-4.5D";
            lenses["48 38 1F 37 34 3C 4B 06"] = "Sigma 12-24mm F4.5-5.6 EX Aspherical DG HSM";
            lenses["48 3C 19 31 30 3C 4B 06"] = "Sigma 10-20mm F4-5.6 EX DC HSM";
            lenses["48 3C 50 A0 30 40 4B 02"] = "Sigma 50-500mm F4-6.3 EX APO RF HSM";
            lenses["48 3C 8E B0 3C 3C 4B 02"] = "Sigma APO 300-800 F5.6 EX DG HSM";
            lenses["48 3C B0 B0 3C 3C 4B 02"] = "Sigma APO 800mm F5.6 EX HSM";
            lenses["48 44 A0 A0 34 34 4B 02"] = "Sigma APO 500mm F4.5 EX HSM";
            lenses["48 48 24 24 24 24 4B 02"] = "Sigma 14mm F2.8 EX Aspherical HSM";
            lenses["48 48 2B 44 24 30 4B 06"] = "Sigma 17-35mm F2.8-4 EX DG Aspherical HSM";
            lenses["48 48 68 8E 30 30 4B 02"] = "Sigma 100-300mm F4 EX IF HSM";
            lenses["48 48 76 76 24 24 4B 06"] = "Sigma 150mm F2.8 EX DG APO Macro HSM";
            lenses["48 48 8E 8E 24 24 4B 02"] = "AF-S Nikkor 300mm f/2.8D IF-ED";
            lenses["48 48 8E 8E 24 24 E1 02"] = "AF-S Nikkor 300mm f/2.8D IF-ED + TC-17E";
            lenses["48 48 8E 8E 24 24 F1 02"] = "AF-S Nikkor 300mm f/2.8D IF-ED + TC-14E";
            lenses["48 48 8E 8E 24 24 F2 02"] = "AF-S Nikkor 300mm f/2.8D IF-ED + TC-20E";
            lenses["48 4C 7C 7C 2C 2C 4B 02"] = "Sigma 180mm F3.5 EX DG Macro";
            lenses["48 4C 7D 7D 2C 2C 4B 02"] = "Sigma APO MACRO 180mm F3.5 EX DG HSM";
            lenses["48 54 3E 3E 0C 0C 4B 06"] = "Sigma 30mm F1.4 EX DC HSM";
            lenses["48 54 5C 80 24 24 4B 02"] = "Sigma 70-200mm F2.8 EX APO IF HSM";
            lenses["48 54 6F 8E 24 24 4B 02"] = "Sigma APO 120-300mm F2.8 EX DG HSM";
            lenses["48 54 8E 8E 24 24 4B 02"] = "Sigma APO 300mm F2.8 EX DG HSM";
            lenses["49 3C A6 A6 30 30 4C 02"] = "AF-S Nikkor 600mm f/4D IF-ED";
            lenses["49 3C A6 A6 30 30 E1 02"] = "AF-S Nikkor 600mm f/4D IF-ED + TC-17E";
            lenses["49 3C A6 A6 30 30 F1 02"] = "AF-S Nikkor 600mm f/4D IF-ED + TC-14E";
            lenses["49 3C A6 A6 30 30 F2 02"] = "AF-S Nikkor 600mm f/4D IF-ED + TC-20E";
            lenses["4A 54 62 62 0C 0C 4D 02"] = "AF Nikkor 85mm f/1.4D IF";
            lenses["4B 3C A0 A0 30 30 4E 02"] = "AF-S Nikkor 500mm f/4D IF-ED";
            lenses["4B 3C A0 A0 30 30 E1 02"] = "AF-S Nikkor 500mm f/4D IF-ED + TC-17E";
            lenses["4B 3C A0 A0 30 30 F1 02"] = "AF-S Nikkor 500mm f/4D IF-ED + TC-14E";
            lenses["4B 3C A0 A0 30 30 F2 02"] = "AF-S Nikkor 500mm f/4D IF-ED + TC-20E";
            lenses["4C 40 37 6E 2C 3C 4F 02"] = "AF Zoom-Nikkor 24-120mm f/3.5-5.6D IF";
            lenses["4D 40 3C 80 2C 3C 62 02"] = "AF Zoom-Nikkor 28-200mm f/3.5-5.6D IF";
            lenses["4D 41 3C 8E 2B 40 62 02"] = "Tamron AF28-300mm f/3.5-6.3 XR Di LD Aspherical (IF)";
            lenses["4D 41 3C 8E 2C 40 62 02"] = "Tamron AF28-300mm f/3.5-6.3 XR LD Aspherical (IF)";
            lenses["4E 48 72 72 18 18 51 02"] = "AF DC-Nikkor 135mm f/2D";
            lenses["4F 40 37 5C 2C 3C 53 06"] = "IX-Nikkor 24-70mm f/3.5-5.6";
            lenses["50 48 56 7C 30 3C 54 06"] = "IX-Nikkor 60-180mm f/4-5.6";
            lenses["53 48 60 80 24 24 57 02"] = "AF Zoom-Nikkor 80-200mm f/2.8D ED";
            lenses["53 48 60 80 24 24 60 02"] = "AF Zoom-Nikkor 80-200mm f/2.8D ED";
            lenses["54 44 5C 7C 34 3C 58 02"] = "AF Zoom-Micro Nikkor 70-180mm f/4.5-5.6D ED";
            lenses["56 3C 5C 8E 30 3C 1C 02"] = "Sigma 70-300mm F4-5.6 APO Macro Super II";
            lenses["56 48 5C 8E 30 3C 5A 02"] = "AF Zoom-Nikkor 70-300mm f/4-5.6D ED";
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
            lenses["67 54 37 5C 24 24 1C 02"] = "Sigma 24-70mm F2.8 EX DG Macro";
            lenses["68 42 3C 60 2A 3C 6E 06"] = "AF Zoom-Nikkor 28-80mm f/3.3-5.6G";
            lenses["69 48 5C 8E 30 3C 6F 02"] = "Tamron AF70-300mm f/4-5.6 LD Macro 1:2";
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
            lenses["77 44 61 98 34 3C 7B 0E"] = "Sigma 80-400mm f4.5-5.6 EX OS";
            lenses["77 48 5C 80 24 24 7B 0E"] = "AF-S VR Zoom-Nikkor 70-200mm f/2.8G IF-ED";
            lenses["78 40 37 6E 2C 3C 7C 0E"] = "AF-S VR Zoom-Nikkor 24-120mm f/3.5-5.6G IF-ED";
            lenses["79 40 11 11 2C 2C 1C 06"] = "Sigma 8mm F3.5 EX";
            lenses["79 40 3C 80 2C 3C 7F 06"] = "AF Zoom-Nikkor 28-200mm f/3.5-5.6G IF-ED";
            lenses["79 48 5C 5C 24 24 1C 06"] = "Sigma 70mm F2.8 EX DG Macro";
            lenses["7A 3B 53 80 30 3C 4B 06"] = "Sigma 55-200mm F4-5.6 DC HSM";
            lenses["7A 3C 1F 37 30 30 7E 06"] = "AF-S DX Zoom-Nikkor 12-24mm f/4G IF-ED";
            lenses["7A 40 2D 50 2C 3C 4B 06"] = "Sigma 18-50mm F3.5-5.6 DC HSM";
            lenses["7A 47 2B 5C 24 34 4B 06"] = "Sigma 17-70mm F2.8-4.5 DC Macro Asp. IF HSM";
            lenses["7A 47 50 76 24 24 4B 06"] = "Sigma APO 50-150mm F2.8 EX DC HSM";
            lenses["7A 48 2B 5C 24 34 4B 06"] = "Sigma 17-70mm F2.8-4.5 DC Macro Asp. IF HSM";
            lenses["7A 48 2D 50 24 24 4B 06"] = "Sigma 18-50mm F2.8 EX DC HSM";
            lenses["7A 48 5C 80 24 24 4B 06"] = "Sigma 70-200mm F2.8 EX APO DG Macro HSM II";
            lenses["7A 54 6E 8E 24 24 4B 02"] = "Sigma APO 120-300mm F2.8 EX DG HSM";
            lenses["7B 48 80 98 30 30 80 0E"] = "AF-S VR Zoom-Nikkor 200-400mm f/4G IF-ED";
            lenses["7D 48 2B 53 24 24 82 06"] = "AF-S DX Zoom-Nikkor 17-55mm f/2.8G IF-ED";
            lenses["7F 40 2D 5C 2C 34 84 06"] = "AF-S DX Zoom-Nikkor 18-70mm f/3.5-4.5G IF-ED";
            lenses["7F 48 2B 5C 24 34 1C 06"] = "Sigma 17-70mm F2.8-4.5 DC Macro Asp. IF";
            lenses["7F 48 2D 50 24 24 1C 06"] = "Sigma 18-50mm F2.8 EX DC Macro";
            lenses["80 48 1A 1A 24 24 85 06"] = "AF DX Fisheye-Nikkor 10.5mm f/2.8G ED";
            lenses["81 54 80 80 18 18 86 0E"] = "AF-S VR Nikkor 200mm f/2G IF-ED";
            lenses["82 48 8E 8E 24 24 87 0E"] = "AF-S VR Nikkor 300mm f/2.8G IF-ED";
            lenses["89 3C 53 80 30 3C 8B 06"] = "AF-S DX Zoom-Nikkor 55-200mm f/4-5.6G ED";
            lenses["8A 54 6A 6A 24 24 8C 0E"] = "AF-S VR Micro-Nikkor 105mm f/2.8G IF-ED";
            lenses["8B 40 2D 80 2C 3C 8D 0E"] = "AF-S DX VR Zoom-Nikkor 18-200mm f/3.5-5.6G IF-ED";
            lenses["8B 40 2D 80 2C 3C FD 0E"] = "AF-S DX VR Zoom-Nikkor 18-200mm f/3.5-5.6G IF-ED";
            lenses["8C 40 2D 53 2C 3C 8E 06"] = "AF-S DX Zoom-Nikkor 18-55mm f/3.5-5.6G ED";
            lenses["8D 44 5C 8E 34 3C 8F 0E"] = "AF-S VR Zoom-Nikkor 70-300mm f/4.5-5.6G IF-ED";
            lenses["8F 40 2D 72 2C 3C 91 06"] = "AF-S DX Zoom-Nikkor 18-135mm f/3.5-5.6G IF-ED";
            lenses["90 3B 53 80 30 3C 92 0E"] = "AF-S DX VR Zoom-Nikkor 55-200mm f/4-5.6G IF-ED";
            lenses["92 48 24 37 24 24 94 06"] = "AF-S Zoom-Nikkor 14-24mm f/2.8G ED";
            lenses["93 48 37 5C 24 24 95 06"] = "AF-S Zoom-Nikkor 24-70mm f/2.8G ED";
            lenses["94 40 2D 53 2C 3C 96 06"] = "AF-S DX Zoom-Nikkor 18-55mm f/3.5-5.6G ED II";
            lenses["95 00 37 37 2C 2C 97 06"] = "PC-E Nikkor 24mm f/3.5D ED";
            lenses["95 4C 37 37 2C 2C 97 02"] = "PC-E Nikkor 24mm f/3.5D ED";
            lenses["96 48 98 98 24 24 98 0E"] = "AF-S VR Nikkor 400mm f/2.8G ED";
            lenses["97 3C A0 A0 30 30 99 0E"] = "AF-S VR Nikkor 500mm f/4G ED";
            lenses["98 3C A6 A6 30 30 9A 0E"] = "AF-S VR Nikkor 600mm f/4G ED";
            lenses["99 40 29 62 2C 3C 9B 0E"] = "AF-S DX VR Zoom-Nikkor 16-85mm f/3.5-5.6G ED";
            lenses["9A 40 2D 53 2C 3C 9C 0E"] = "AF-S DX VR Zoom-Nikkor 18-55mm f/3.5-5.6G";
            lenses["9B 00 4C 4C 24 24 9D 06"] = "PC-E Micro Nikkor 45mm f/2.8D ED";
            lenses["9B 54 4C 4C 24 24 9D 02"] = "PC-E Micro Nikkor 45mm f/2.8D ED";
            lenses["9C 54 56 56 24 24 9E 06"] = "AF-S Micro Nikkor 60mm f/2.8G ED";
            lenses["9D 00 62 62 24 24 9F 06"] = "PC-E Micro Nikkor 85mm f/2.8D";
            lenses["9D 54 62 62 24 24 9F 02"] = "PC-E Micro Nikkor 85mm f/2.8D";
            lenses["9E 40 2D 6A 2C 3C A0 0E"] = "AF-S DX VR Zoom-Nikkor 18-105mm f/3.5-5.6G ED";
            lenses["9F 58 44 44 14 14 A1 06"] = "AF-S DX Nikkor 35mm f/1.8G";
            lenses["A0 54 50 50 0C 0C A2 06"] = "AF-S Nikkor 50mm f/1.4G";
            lenses["A1 40 18 37 2C 34 A3 06"] = "AF-S DX Nikkor 10-24mm f/3.5-4.5G ED";
            lenses["A1 41 19 31 2C 2C 4B 06"] = "Sigma 10-20mm F3.5 EX DC HSM";
            lenses["A2 48 5C 80 24 24 A4 0E"] = "AF-S Nikkor 70-200mm f/2.8G ED VR II";
            lenses["A5 40 2D 88 2C 40 4B 0E"] = "Sigma 18-250mm F3.5-6.3 DC OS HSM";
            lenses["A6 48 37 5C 24 24 4B 06"] = "Sigma 24-70mm F2.8 IF EX DG HSM";
            lenses["B6 48 37 56 24 24 1C 02"] = "Sigma 24-60mm F2.8 EX DG";
            lenses["CD 3D 2D 70 2E 3C 4B 0E"] = "Sigma 18-125mm F3.8-5.6 DC OS HSM";
            lenses["CE 34 76 A0 38 40 4B 0E"] = "Sigma 150-500mm F5-6.3 DG OS APO HSM";
            lenses["CF 38 6E 98 34 3C 4B 0E"] = "Sigma APO 120-400mm F4.5-5.6 DG OS HSM";
            lenses["DC 48 19 19 24 24 4B 06"] = "Sigma 10mm F2.8 EX DC HSM Fisheye";
            lenses["DE 54 50 50 0C 0C 4B 06"] = "Sigma 50mm F1.4 EX DG HSM";
            lenses["E0 3C 5C 8E 30 3C 4B 06"] = "Sigma 70-300mm F4-5.6 APO DG Macro HSM";
            lenses["E1 58 37 37 14 14 1C 02"] = "Sigma 24mm F1.8 EX DG Aspherical Macro";
            lenses["E5 54 6A 6A 24 24 35 02"] = "Sigma Macro 105mm F2.8 EX DG";
            lenses["E9 54 37 5C 24 24 1C 02"] = "Sigma 24-70mm F2.8 EX DG Macro";
            lenses["ED 40 2D 80 2C 40 4B 0E"] = "Sigma 18-200mm F3.5-6.3 DC OS HSM";
            lenses["EE 48 5C 80 24 24 4B 06"] = "Sigma 70-200mm F2.8 EX APO DG Macro HSM II";
            lenses["F0 38 1F 37 34 3C 4B 06"] = "Sigma 12-24mm F4.5-5.6 EX DG Aspherical HSM";
            lenses["F3 54 2B 50 24 24 84 0E"] = "Tamron SP AF 17-50mm F/2.8 XR Di II VC LD Aspherical (IF) (B005)";
            lenses["F4 54 56 56 18 18 84 06"] = "Tamron SP AF 60mm f/2.0 Di II Macro 1:1 (G005)";
            lenses["F5 40 2C 8A 2C 40 40 0E"] = "Tamron AF 18-270mm f/3.5-6.3 Di II VC LD Aspherical (IF) Macro (B003)";
            lenses["F5 48 76 76 24 24 4B 06"] = "Sigma 150mm F2.8 EX DG APO Macro HSM";
            lenses["F6 3F 18 37 2C 34 84 06"] = "Tamron SP AF 10-24mm f/3.5-4.5 Di II LD Aspherical (IF) (B001)";
            lenses["F7 53 5C 80 24 24 84 06"] = "Tamron SP AF 70-200mm f/2.8 Di LD (IF) Macro (A001)";
            lenses["F8 54 3E 3E 0C 0C 4B 06"] = "Sigma 30mm F1.4 EX DC HSM";
            lenses["F8 55 64 64 24 24 84 06"] = "Tamron SP AF 90mm f/2.8 Di Macro 1:1 (272NII)";
            lenses["F9 3C 19 31 30 3C 4B 06"] = "Sigma 10-20mm F4-5.6 EX DC HSM";
            lenses["F9 40 3C 8E 2C 40 40 0E"] = "Tamron AF 28-300mm f/3.5-6.3 XR Di VC LD Aspherical (IF) Macro (A20)";
            lenses["FA 54 3C 5E 24 24 84 06"] = "Tamron SP AF 28-75mm f/2.8 XR Di LD Aspherical (IF) Macro (A09NII)";
            lenses["FB 54 8E 8E 24 24 4B 02"] = "Sigma APO 300mm F2.8 EX DG HSM";
            lenses["FD 47 50 76 24 24 4B 06"] = "Sigma 50-150mm F2.8 EX APO DC HSM II";
            lenses["FE 47 00 00 24 24 4B 06"] = "Sigma 4.5mm F2.8 EX DC Circular Fisheye HSM";
        }
        virtual std::string toString (Tag* t) {

          static const unsigned char xlat[2][256] = {
          { 0xc1,0xbf,0x6d,0x0d,0x59,0xc5,0x13,0x9d,0x83,0x61,0x6b,0x4f,0xc7,0x7f,0x3d,0x3d,
            0x53,0x59,0xe3,0xc7,0xe9,0x2f,0x95,0xa7,0x95,0x1f,0xdf,0x7f,0x2b,0x29,0xc7,0x0d,
            0xdf,0x07,0xef,0x71,0x89,0x3d,0x13,0x3d,0x3b,0x13,0xfb,0x0d,0x89,0xc1,0x65,0x1f,
            0xb3,0x0d,0x6b,0x29,0xe3,0xfb,0xef,0xa3,0x6b,0x47,0x7f,0x95,0x35,0xa7,0x47,0x4f,
            0xc7,0xf1,0x59,0x95,0x35,0x11,0x29,0x61,0xf1,0x3d,0xb3,0x2b,0x0d,0x43,0x89,0xc1,
            0x9d,0x9d,0x89,0x65,0xf1,0xe9,0xdf,0xbf,0x3d,0x7f,0x53,0x97,0xe5,0xe9,0x95,0x17,
            0x1d,0x3d,0x8b,0xfb,0xc7,0xe3,0x67,0xa7,0x07,0xf1,0x71,0xa7,0x53,0xb5,0x29,0x89,
            0xe5,0x2b,0xa7,0x17,0x29,0xe9,0x4f,0xc5,0x65,0x6d,0x6b,0xef,0x0d,0x89,0x49,0x2f,
            0xb3,0x43,0x53,0x65,0x1d,0x49,0xa3,0x13,0x89,0x59,0xef,0x6b,0xef,0x65,0x1d,0x0b,
            0x59,0x13,0xe3,0x4f,0x9d,0xb3,0x29,0x43,0x2b,0x07,0x1d,0x95,0x59,0x59,0x47,0xfb,
            0xe5,0xe9,0x61,0x47,0x2f,0x35,0x7f,0x17,0x7f,0xef,0x7f,0x95,0x95,0x71,0xd3,0xa3,
            0x0b,0x71,0xa3,0xad,0x0b,0x3b,0xb5,0xfb,0xa3,0xbf,0x4f,0x83,0x1d,0xad,0xe9,0x2f,
            0x71,0x65,0xa3,0xe5,0x07,0x35,0x3d,0x0d,0xb5,0xe9,0xe5,0x47,0x3b,0x9d,0xef,0x35,
            0xa3,0xbf,0xb3,0xdf,0x53,0xd3,0x97,0x53,0x49,0x71,0x07,0x35,0x61,0x71,0x2f,0x43,
            0x2f,0x11,0xdf,0x17,0x97,0xfb,0x95,0x3b,0x7f,0x6b,0xd3,0x25,0xbf,0xad,0xc7,0xc5,
            0xc5,0xb5,0x8b,0xef,0x2f,0xd3,0x07,0x6b,0x25,0x49,0x95,0x25,0x49,0x6d,0x71,0xc7 },
          { 0xa7,0xbc,0xc9,0xad,0x91,0xdf,0x85,0xe5,0xd4,0x78,0xd5,0x17,0x46,0x7c,0x29,0x4c,
            0x4d,0x03,0xe9,0x25,0x68,0x11,0x86,0xb3,0xbd,0xf7,0x6f,0x61,0x22,0xa2,0x26,0x34,
            0x2a,0xbe,0x1e,0x46,0x14,0x68,0x9d,0x44,0x18,0xc2,0x40,0xf4,0x7e,0x5f,0x1b,0xad,
            0x0b,0x94,0xb6,0x67,0xb4,0x0b,0xe1,0xea,0x95,0x9c,0x66,0xdc,0xe7,0x5d,0x6c,0x05,
            0xda,0xd5,0xdf,0x7a,0xef,0xf6,0xdb,0x1f,0x82,0x4c,0xc0,0x68,0x47,0xa1,0xbd,0xee,
            0x39,0x50,0x56,0x4a,0xdd,0xdf,0xa5,0xf8,0xc6,0xda,0xca,0x90,0xca,0x01,0x42,0x9d,
            0x8b,0x0c,0x73,0x43,0x75,0x05,0x94,0xde,0x24,0xb3,0x80,0x34,0xe5,0x2c,0xdc,0x9b,
            0x3f,0xca,0x33,0x45,0xd0,0xdb,0x5f,0xf5,0x52,0xc3,0x21,0xda,0xe2,0x22,0x72,0x6b,
            0x3e,0xd0,0x5b,0xa8,0x87,0x8c,0x06,0x5d,0x0f,0xdd,0x09,0x19,0x93,0xd0,0xb9,0xfc,
            0x8b,0x0f,0x84,0x60,0x33,0x1c,0x9b,0x45,0xf1,0xf0,0xa3,0x94,0x3a,0x12,0x77,0x33,
            0x4d,0x44,0x78,0x28,0x3c,0x9e,0xfd,0x65,0x57,0x16,0x94,0x6b,0xfb,0x59,0xd0,0xc8,
            0x22,0x36,0xdb,0xd2,0x63,0x98,0x43,0xa1,0x04,0x87,0x86,0xf7,0xa6,0x26,0xbb,0xd6,
            0x59,0x4d,0xbf,0x6a,0x2e,0xaa,0x2b,0xef,0xe6,0x78,0xb6,0x4e,0xe0,0x2f,0xdc,0x7c,
            0xbe,0x57,0x19,0x32,0x7e,0x2a,0xd0,0xb8,0xba,0x29,0x00,0x3c,0x52,0x7d,0xa8,0x49,
            0x3b,0x2d,0xeb,0x25,0x49,0xfa,0xa3,0xaa,0x39,0xa7,0xc5,0xa7,0x50,0x11,0x36,0xfb,
            0xc6,0x67,0x4a,0xf5,0xa5,0x12,0x65,0x7e,0xb0,0xdf,0xaf,0x4e,0xb3,0x61,0x7f,0x2f } };

            int ver = (t->toInt (0, BYTE) - '0') * 1000 + (t->toInt (1, BYTE) - '0') * 100 + (t->toInt (2, BYTE) - '0') * 10 + (t->toInt (3, BYTE) - '0');
            
            std::ostringstream ld;
            ld << "Version = " << ver << std::endl;
            
            int lenstype = t->getParent()->getTag(0x0083)->toInt(0,BYTE);
            
            std::ostringstream lid;
            lid.setf (std::ios_base::hex, std::ios_base::basefield);
            lid.setf (std::ios_base::uppercase);

            std::string model = t->getParent()->getParent()->getParent()->getTag(0x0110)->valueToString();
            int lidoffs = 7;
            bool d100 = false;
            if (model.substr(0,10)=="NIKON D100" || model.substr(0,9)=="NIKON D1X") {
                 lidoffs = 0;
                 d100 = true;
            }

            unsigned char buffer[15];
            if (d100) 
                memcpy (buffer, t->getValue()+6, 7);
            else
                memcpy (buffer, t->getValue()+4, 15);

            if (ver>=201) {
                const unsigned char* serval = t->getParent()->getTag(0x001d)->getValue ();
                int serial = 0;
                for (int i=0; serval[i]; i++)
                    serial = serial*10 + (isdigit(serval[i]) ? serval[i] - '0' : serval[i] % 10);
                const unsigned char* scval = t->getParent()->getTag(0x00a7)->getValue ();
                int key = 0;
                for (int i=0; i<4; i++)
                    key ^= scval[i];

                unsigned char ci = xlat[0][serial & 0xff];
                unsigned char cj = xlat[1][key];
                unsigned char ck = 0x60;
                for (int i=0; i < 15; i++)
                    buffer[i] ^= (cj += ci * ck++);
            }
                
            if (!d100) {
                ld << "ExitPupilPosition = " << (int) buffer[0] << std::endl;
                ld << "AFAperture = "        << (int) buffer[1] << std::endl;
                ld << "FocusPosition = "     << (int) buffer[4] << std::endl;
                ld << "FocusDistance = "     << (int) buffer[5] << std::endl;
                ld << "FocalLength = "       << (int) buffer[6] << std::endl;
                ld << "EffectiveMaxAperture = "  << (int) buffer[14] << std::endl;
            }
                
            for (int i=0; i<7; i++)
                 lid << std::setw(2) << std::setfill('0') << (int)buffer[lidoffs+i] << ' ';
            lid << std::setw(2) << std::setfill('0') << lenstype;            
            
            std::map<std::string,std::string>::iterator r = lenses.find (lid.str());
            if (r!=lenses.end()) 
                ld << "Lens = " << r->second;
            else
                ld << "Lens = Unknown, ID=" << lid.str();
        
            return ld.str();
        }
        
};
NALensDataInterpreter naLensDataInterpreter;

const TagAttrib nikon2Attribs[] = {
 0, 1, 0, 0, 0x0002, "Unknown", &stdInterpreter,
 0, 1, 0, 0, 0x0003, "Quality", &stdInterpreter,
 0, 1, 0, 0, 0x0004, "ColorMode", &stdInterpreter,
 0, 1, 0, 0, 0x0005, "ImageAdjustment", &stdInterpreter,
 0, 1, 0, 0, 0x0006, "ISOSpeed", &naISOInterpreter,
 0, 1, 0, 0, 0x0007, "WhiteBalance", &stdInterpreter,
 0, 1, 0, 0, 0x0008, "Focus", &stdInterpreter,
 0, 1, 0, 0, 0x0009, "Unknown", &stdInterpreter,
 0, 1, 0, 0, 0x000a, "DigitalZoom", &stdInterpreter,
 0, 1, 0, 0, 0x000b, "AuxiliaryLens", &stdInterpreter,
 0, 1, 0, 0, 0x0f00, "Unknown", &stdInterpreter,
 -1, 0, 0,  0, 0, "", NULL};

const TagAttrib nikon3Attribs[] = {
 0, 1, 0, 0, 0x0001, "MakerNoteVersion", &stdInterpreter,
 0, 1, 0, 0, 0x0002, "ISOSpeed", &naISOInterpreter,
 0, 1, 0, 0, 0x0003, "ColorMode", &stdInterpreter,
 0, 1, 0, 0, 0x0004, "Quality", &stdInterpreter,
 0, 1, 0, 0, 0x0005, "WhiteBalance", &stdInterpreter,
 0, 1, 0, 0, 0x0006, "Sharpness", &stdInterpreter,
 0, 1, 0, 0, 0x0007, "FocusMode", &stdInterpreter,
 0, 1, 0, 0, 0x0008, "FlashSetting", &stdInterpreter,
 0, 1, 0, 0, 0x0009, "FlashType", &stdInterpreter,
 0, 1, 0, 0, 0x000b, "WhiteBalanceFineTune", &stdInterpreter,
 0, 3, 0, 0, 0x000c, "ColorBalance1", &stdInterpreter,
 0, 1, 0, 0, 0x000d, "ProgramShift", &stdInterpreter,
 0, 1, 0, 0, 0x000e, "ExposureDifference", &stdInterpreter,
 0, 1, 0, 0, 0x000f, "ISOSelection", &naISOInterpreter,
 0, 1, 0, 0, 0x0010, "DataDump", &stdInterpreter,
 1, 1, 0, 0, 0x0011, "NikonPreview", &stdInterpreter,
 0, 1, 0, 0, 0x0012, "FlashExposureComp", &stdInterpreter,
 0, 1, 0, 0, 0x0013, "ISOSetting", &stdInterpreter,
 0, 1, 0, 0, 0x0016, "ImageBoundary", &stdInterpreter,
 0, 1, 0, 0, 0x0018, "FlashExposureBracketValue", &stdInterpreter,
 0, 1, 0, 0, 0x0019, "ExposureBracketValue", &stdInterpreter,
 0, 1, 0, 0, 0x001a, "ImageProcessing", &stdInterpreter,
 0, 1, 0, 0, 0x001b, "CropHiSpeed", &stdInterpreter,
 0, 1, 0, 0, 0x001d, "SerialNumber", &stdInterpreter,
 0, 1, 0, 0, 0x001e, "ColorSpace", &stdInterpreter,
 0, 1, 0, 0, 0x0020, "ImageAuthentication", &stdInterpreter,
 0, 1, 0, 0, 0x0080, "ImageAdjustment", &stdInterpreter,
 0, 1, 0, 0, 0x0081, "ToneComp", &stdInterpreter,
 0, 1, 0, 0, 0x0082, "AuxiliaryLens", &stdInterpreter,
 0, 1, 0, 0, 0x0083, "LensType", &naLensTypeInterpreter,
 0, 1, 0, 0, 0x0084, "Lens", &stdInterpreter,
 0, 1, 0, 0, 0x0085, "ManualFocusDistance", &stdInterpreter,
 0, 1, 0, 0, 0x0086, "DigitalZoom", &stdInterpreter,
 0, 1, 0, 0, 0x0087, "FlashMode", &naFlashModeInterpreter,
 0, 1, 0, 0, 0x0088, "AFInfo", &naAFInfoInterpreter,
 0, 1, 0, 0, 0x0089, "ShootingMode", &naShootingModeInterpreter,
 0, 1, 0, 0, 0x008a, "AutoBracketRelease", &stdInterpreter,
 0, 1, 0, 0, 0x008b, "LensFStops", &stdInterpreter,
 0, 1, 0, 0, 0x008c, "NEFCurve1", &stdInterpreter,
 0, 1, 0, 0, 0x008d, "ColorHue", &stdInterpreter,
 0, 1, 0, 0, 0x008f, "SceneMode", &stdInterpreter,
 0, 1, 0, 0, 0x0090, "LightSource", &stdInterpreter,
 0, 1, 0, 0, 0x0091, "ShotInfo", &stdInterpreter,
 0, 1, 0, 0, 0x0092, "HueAdjustment", &stdInterpreter,
 0, 1, 0, 0, 0x0094, "Saturation", &stdInterpreter,
 0, 1, 0, 0, 0x0095, "NoiseReduction", &stdInterpreter,
 0, 1, 0, 0, 0x0096, "NEFCurve2", &stdInterpreter,
 0, 3, 0, 0, 0x0097, "ColorBalance", &stdInterpreter,
 0, 1, 0, 0, 0x0098, "LensData", &naLensDataInterpreter,
 0, 1, 0, 0, 0x0099, "RawImageCenter", &stdInterpreter,
 0, 1, 0, 0, 0x009a, "SensorPixelSize", &stdInterpreter,
 0, 1, 0, 0, 0x00a0, "SerialNumber", &stdInterpreter,
 0, 1, 0, 0, 0x00a2, "ImageDataSize", &stdInterpreter,
 0, 1, 0, 0, 0x00a5, "ImageCount", &stdInterpreter,
 0, 1, 0, 0, 0x00a6, "DeletedImageCount", &stdInterpreter,
 0, 1, 0, 0, 0x00a7, "ShutterCount", &stdInterpreter,
 0, 1, 0, 0, 0x00a9, "ImageOptimization", &stdInterpreter,
 0, 1, 0, 0, 0x00aa, "Saturation", &stdInterpreter,
 0, 1, 0, 0, 0x00ab, "VariProgram", &stdInterpreter,
 0, 1, 0, 0, 0x00ac, "ImageStabilization", &stdInterpreter,
 0, 1, 0, 0, 0x00ad, "AFResponse", &stdInterpreter,
 0, 1, 0, 0, 0x00b0, "MultiExposure", &stdInterpreter,
 0, 1, 0, 0, 0x00b1, "HighISONoiseReduction", &naHiISONRInterpreter,
 0, 1, 0, 0, 0x0e00, "PrintIM", &stdInterpreter,
 0, 0, 0, 0, 0x0e01, "NikonCaptureData", &stdInterpreter,
 0, 0, 0, 0, 0x0e09, "NikonCaptureVersion", &stdInterpreter,
 0, 0, 0, 0, 0x0e0e, "NikonCaptureOffsets", &stdInterpreter,
 0, 0, 0, 0, 0x0e10, "NikonScanIFD", &stdInterpreter,
 -1, 0, 0,  0, 0, "", NULL};

}
#endif

