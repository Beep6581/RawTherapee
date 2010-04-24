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
#ifndef _CANONATTRIBS_
#define _CANONATTRIBS_

#include <rtexif.h>
#include <string>
#include <map>
#include <math.h>
#include <sstream>
#include <iomanip>

namespace rtexif {

class CAIntSerNumInterpreter : public Interpreter {
    public:
        CAIntSerNumInterpreter () {}
        virtual std::string toString (Tag* t) { return ""; }
};

CAIntSerNumInterpreter caIntSerNumInterpreter;

class CAFocalLengthInterpreter : public Interpreter {
    public:
        CAFocalLengthInterpreter () {}
        virtual std::string toString (Tag* t) {
            std::ostringstream str;
            str << "FocalType  = " << t->toInt(0,SHORT) << std::endl;
            str << "FocalLength = " << t->toInt(2,SHORT) << std::endl;
            str << "FocalPlaneXSize = " << t->toInt(4,SHORT) << std::endl;
            str << "FocalPlaneYSize = " << t->toInt(6,SHORT);
            return str.str();
        }
};
CAFocalLengthInterpreter caFocalLengthInterpreter;


class CACameraSettingsInterpreter : public Interpreter {
        std::map<int,std::string> machoices;
        std::map<int,std::string> qlchoices;
        std::map<int,std::string> fmchoices;
        std::map<int,std::string> cdchoices;
        std::map<int,std::string> fochoices;
        std::map<int,std::string> rmchoices;
        std::map<int,std::string> ischoices;
        std::map<int,std::string> mmchoices;
        std::map<int,std::string> dzchoices;
        std::map<int,std::string> emchoices;
        std::map<int,std::string> frchoices;
        std::map<int,std::string> afchoices;
        std::map<int,std::string> exchoices;
        std::map<int,std::string> fcchoices;
        std::map<int,std::string> aechoices;
        std::map<int,std::string> stchoices;
        std::map<int,std::string> smchoices;
        std::map<int,std::string> pechoices;
        std::map<int,std::string> mfchoices;
        std::map<int,std::string> choices;
    public:
        CACameraSettingsInterpreter () {
            machoices[1] = "Macro";
            machoices[2] = "Normal";
            qlchoices[1] = "Economy";
            qlchoices[2] = "Normal";
            qlchoices[3] = "Fine";
            qlchoices[4] = "RAW";
            qlchoices[5] = "Superfine";
            fmchoices[0] = "Off";
            fmchoices[1] = "Auto";
            fmchoices[2] = "On";
            fmchoices[3] = "Red-eye reduction";
            fmchoices[4] = "Slow-sync";
            fmchoices[5] = "Red-eye reduction (Auto)";
            fmchoices[6] = "Red-eye reduction (On)";
            fmchoices[16] = "External flash";
            cdchoices[0] = "Single";
            cdchoices[1] = "Continuous";
            cdchoices[2] = "Movie";
            cdchoices[3] = "Continuous, Speed Priority";
            cdchoices[4] = "Continuous, Low";
            cdchoices[5] = "Continuous, High";
            fochoices[0] = "One-shot AF";
            fochoices[1] = "AI Servo AF";
            fochoices[2] = "AI Focus AF";
            fochoices[3] = "Manual Focus";
            fochoices[4] = "Single";
            fochoices[5] = "Continuous";
            fochoices[6] = "Manual Focus";
            fochoices[16] = "Pan Focus";
            rmchoices[1] = "JPEG";
            rmchoices[2] = "CRW+THM";
            rmchoices[3] = "AVI+THM";
            rmchoices[4] = "TIF";
            rmchoices[5] = "TIF+JPEG";
            rmchoices[6] = "CR2";
            rmchoices[7] = "CR2+JPEG";
            ischoices[0] = "Large";
            ischoices[1] = "Medium";
            ischoices[2] = "Small";
            ischoices[5] = "Medium 1";
            ischoices[6] = "Medium 2";
            ischoices[7] = "Medium 3";
            ischoices[8] = "Postcard";
            ischoices[9] = "Widescreen";
            emchoices[0] = "Full auto ";
            emchoices[1] = "Manual ";
            emchoices[2] = "Landscape ";
            emchoices[3] = "Fast shutter ";
            emchoices[4] = "Slow shutter ";
            emchoices[5] = "Night ";
            emchoices[6] = "Gray Scale ";
            emchoices[7] = "Sepia ";
            emchoices[8] = "Portrait ";
            emchoices[9] = "Sports ";
            emchoices[10] = "Macro ";
            emchoices[11] = "Black & White";
            emchoices[12] = "Pan focus";
            emchoices[13] = "Vivid";
            emchoices[14] = "Neutral";
            emchoices[15] = "Flash Off";
            emchoices[16] = "Long Shutter";
            emchoices[17] = "Super Macro";
            emchoices[18] = "Foliage";
            emchoices[19] = "Indoor";
            emchoices[20] = "Fireworks";
            emchoices[21] = "Beach";
            emchoices[22] = "Underwater";
            emchoices[23] = "Snow";
            emchoices[24] = "Kids & Pets";
            emchoices[25] = "Night Snapshot";
            emchoices[26] = "Digital Macro";
            emchoices[27] = "My Colors";
            emchoices[28] = "Still Image";
            emchoices[30] = "Color Accent";
            emchoices[31] = "Color Swap";
            emchoices[32] = "Aquarium";
            emchoices[33] = "ISO 3200";
            dzchoices[0] = "None";
            dzchoices[1] = "2x";
            dzchoices[2] = "4x";
            dzchoices[3] = "Other";
            mmchoices[0] = "Default";
            mmchoices[1] = "Spot";
            mmchoices[2] = "Average";
            mmchoices[3] = "Evaluative";
            mmchoices[4] = "Partial";
            mmchoices[5] = "Center-weighted averaging";
            frchoices[0] = "Manual";
            frchoices[1] = "Auto";
            frchoices[2] = "Not Known";
            frchoices[3] = "Macro";
            frchoices[4] = "Very Close";
            frchoices[5] = "Close";
            frchoices[6] = "Middle Range";
            frchoices[7] = "Far Range";
            frchoices[8] = "Pan Focus";
            frchoices[9] = "Super Macro";
            frchoices[10] = "Infinity";
            afchoices[0x2005] = "Manual AF point selection ";
            afchoices[0x3000] = "None (MF)";
            afchoices[0x3001] = "Auto AF point selection ";
            afchoices[0x3002] = "Right ";
            afchoices[0x3003] = "Center ";
            afchoices[0x3004] = "Left ";
            afchoices[0x4001] = "Auto AF point selection ";
            afchoices[0x4006] = "Face Detect";
            exchoices[0] = "Easy";
            exchoices[1] = "Program AE";
            exchoices[2] = "Shutter speed priority AE";
            exchoices[3] = "Aperture-priority AE";
            exchoices[4] = "Manual";
            exchoices[5] = "Depth-of-field AE";
            exchoices[6] = "M-Dep";
            fcchoices[0] = "Single";
            fcchoices[1] = "Continuous";
            aechoices[0] = "Normal AE";
            aechoices[1] = "Exposure Compensation";
            aechoices[2] = "AE Lock";
            aechoices[3] = "AE Lock + Exposure Comp.";
            aechoices[4] = "No AE";
            stchoices[0] = "Off";
            stchoices[1] = "On";
            stchoices[2] = "On, Shot Only";
            stchoices[3] = "On, Panning";
            smchoices[0] = "Center";
            smchoices[1] = "AF Point";
            pechoices[0] = "Off";
            pechoices[1] = "Vivid";
            pechoices[2] = "Neutral";
            pechoices[3] = "Smooth";
            pechoices[4] = "Sepia";
            pechoices[5] = "B&W";
            pechoices[6] = "Custom";
            pechoices[100] = "My Color Data";
            mfchoices[0] = "N/A";
            mfchoices[0x500] = "Full";
            mfchoices[0x502] = "Medium";
            mfchoices[0x504] = "Low";
            mfchoices[0x7fff] = "N/A";
            choices[1] = "Canon EF 50mm f/1.8";
            choices[2] = "Canon EF 28mm f/2.8";
            choices[3] = "Canon EF 135mm f/2.8 Soft";
            choices[4] = "Canon EF 35-105mm f/3.5-4.5 or Sigma Lens";
            choices[5] = "Canon EF 35-70mm f/3.5-4.5";
            choices[6] = "Canon EF 28-70mm f/3.5-4.5 or Sigma or Tokina Lens";
            choices[7] = "Canon EF 100-300mm F5.6L";
            choices[8] = "Canon EF 100-300mm f/5.6 or Sigma or Tokina Lens";
            choices[9] = "Canon EF 70-210mm f/4 orSigma Lens";
            choices[10] = "Canon EF 50mm f/2.5 Macro or Sigma Lens";
            choices[11] = "Canon EF 35mm f/2";
            choices[13] = "Canon EF 15mm f/2.8";
            choices[14] = "Canon EF 50-200mm f/3.5-4.5L";
            choices[15] = "Canon EF 50-200mm f/3.5-4.5";
            choices[16] = "Canon EF 35-135mm f/3.5-4.5";
            choices[17] = "Canon EF 35-70mm f/3.5-4.5A";
            choices[18] = "Canon EF 28-70mm f/3.5-4.5";
            choices[20] = "Canon EF 100-200mm f/4.5A";
            choices[21] = "Canon EF 80-200mm f/2.8L";
            choices[22] = "Canon EF 20-35mm f/2.8L or Tokina 28-80mm F2.8";
            choices[23] = "Canon EF 35-105mm f/3.5-4.5";
            choices[24] = "Canon EF 35-80mm f/4-5.6 Power Zoom";
            choices[25] = "Canon EF 35-80mm f/4-5.6 Power Zoom";
            choices[26] = "Canon EF 100mm f/2.8 Macro or Cosina 100mm f/3.5 Macro AF or Tamron";
            choices[28] = "Tamron AF Aspherical 28-200mm f/3.8-5.6 or 28-75mm f/2.8 or 28-105mm f/2.8";
            choices[27] = "Canon EF 35-80mm f/4-5.6";
            choices[28] = "Canon EF 80-200mm f/4.5-5.6 or Tamron Lens";
            choices[29] = "Canon EF 50mm f/1.8 MkII";
            choices[30] = "Canon EF 35-105mm f/4.5-5.6";
            choices[31] = "Tamron SP AF 300mm f/2.8 LD IF";
            choices[32] = "Canon EF 24mm f/2.8 or Sigma 15mm f/2.8 EX Fisheye";
            choices[33] = "Voigtlander Ultron 40mm f/2 SLII Aspherical";
            choices[35] = "Canon EF 35-80mm f/4-5.6";
            choices[36] = "Canon EF 38-76mm f/4.5-5.6";
            choices[37] = "Canon EF 35-80mm f/4-5.6 or Tamron Lens";
            choices[38] = "Canon EF 80-200mm f/4.5-5.6";
            choices[39] = "Canon EF 75-300mm f/4-5.6";
            choices[40] = "Canon EF 28-80mm f/3.5-5.6";
            choices[41] = "Canon EF 28-90mm f/4-5.6";
            choices[42] = "Canon EF 28-200mm f/3.5-5.6 or Tamron AF 28-300mm f/3.5-6.3";
            choices[43] = "Canon EF 28-105mm f/4-5.6";
            choices[44] = "Canon EF 90-300mm f/4.5-5.6";
            choices[45] = "Canon EF-S 18-55mm f/3.5-5.6";
            choices[46] = "Canon EF 28-90mm f/4-5.6";
            choices[48] = "Canon EF-S 18-55mm f/3.5-5.6 IS";
            choices[49] = "Canon EF-S 55-250mm f/4-5.6 IS";
            choices[50] = "Canon EF-S 18-200mm f/3.5-5.6 IS";
            choices[51] = "Canon EF-S 18-135mm f/3.5-5.6 IS";
            choices[94] = "Canon TS-E 17mm f/4L";
            choices[95] = "Canon TS-E 24.0mm f/3.5 L II";
            choices[124] = "Canon MP-E 65mm f/2.8 1-5x Macro Photo";
            choices[125] = "Canon TS-E 24mm f/3.5L";
            choices[126] = "Canon TS-E 45mm f/2.8";
            choices[127] = "Canon TS-E 90mm f/2.8";
            choices[129] = "Canon EF 300mm f/2.8L";
            choices[130] = "Canon EF 50mm f/1.0L";
            choices[131] = "Canon EF 28-80mm f/2.8-4L or Sigma Lens";
            choices[132] = "Canon EF 1200mm f/5.6L";
            choices[134] = "Canon EF 600mm f/4L IS";
            choices[135] = "Canon EF 200mm f/1.8L";
            choices[136] = "Canon EF 300mm f/2.8L";
            choices[137] = "Canon EF 85mm f/1.2L or Sigma Lens";
            choices[138] = "Canon EF 28-80mm f/2.8-4L";
            choices[139] = "Canon EF 400mm f/2.8L";
            choices[140] = "Canon EF 500mm f/4.5L";
            choices[141] = "Canon EF 500mm f/4.5L";
            choices[142] = "Canon EF 300mm f/2.8L IS";
            choices[143] = "Canon EF 500mm f/4L IS";
            choices[144] = "Canon EF 35-135mm f/4-5.6 USM";
            choices[145] = "Canon EF 100-300mm f/4.5-5.6 USM";
            choices[146] = "Canon EF 70-210mm f/3.5-4.5 USM";
            choices[147] = "Canon EF 35-135mm f/4-5.6 USM";
            choices[148] = "Canon EF 28-80mm f/3.5-5.6 USM";
            choices[149] = "Canon EF 100mm f/2";
            choices[150] = "Canon EF 14mm f/2.8L or Sigma Lens";
            choices[151] = "Canon EF 200mm f/2.8L";
            choices[152] = "Canon EF 300mm f/4L IS or Sigma Lens";
            choices[153] = "Canon EF 35-350mm f/3.5-5.6L or Tamron or Sigma Lens";
            choices[154] = "Canon EF 20mm f/2.8 USM";
            choices[155] = "Canon EF 85mm f/1.8 USM";
            choices[156] = "Canon EF 28-105mm f/3.5-4.5 USM";
            choices[160] = "Canon EF 20-35mm f/3.5-4.5 USM or Tamron AF 19-35mm f/3.5-4.5";
            choices[161] = "Canon EF 28-70mm f/2.8L or Sigma or Tamron Lens";
            choices[162] = "Canon EF 200mm f/2.8L";
            choices[163] = "Canon EF 300mm f/4L";
            choices[164] = "Canon EF 400mm f/5.6L";
            choices[165] = "Canon EF 70-200mm f/2.8 L";
            choices[166] = "Canon EF 70-200mm f/2.8 L + x1.4";
            choices[167] = "Canon EF 70-200mm f/2.8 L + x2";
            choices[168] = "Canon EF 28mm f/1.8 USM";
            choices[169] = "Canon EF17-35mm f/2.8L or Sigma Lens";
            choices[170] = "Canon EF 200mm f/2.8L II";
            choices[171] = "Canon EF 300mm f/4L";
            choices[172] = "Canon EF 400mm f/5.6L";
            choices[173] = "Canon EF 180mm Macro f/3.5L or Sigma 180mm F3.5 or 150mm f/2.8 Macro";
            choices[174] = "Canon EF 135mm f/2L";
            choices[175] = "Canon EF 400mm f/2.8L";
            choices[176] = "Canon EF 24-85mm f/3.5-4.5 USM";
            choices[177] = "Canon EF 300mm f/4L IS";
            choices[178] = "Canon EF 28-135mm f/3.5-5.6 IS";
            choices[179] = "Canon EF 24mm f/1.4L USM";
            choices[180] = "Canon EF 35mm f/1.4L";
            choices[181] = "Canon EF 100-400mm f/4.5-5.6L IS + x1.4";
            choices[182] = "Canon EF 100-400mm f/4.5-5.6L IS + x2";
            choices[183] = "Canon EF 100-400mm f/4.5-5.6L IS";
            choices[184] = "Canon EF 400mm f/2.8L + x2";
            choices[185] = "Canon EF 600mm f/4L IS";
            choices[186] = "Canon EF 70-200mm f/4L";
            choices[187] = "Canon EF 70-200mm f/4L + 1.4x";
            choices[188] = "Canon EF 70-200mm f/4L + 2x";
            choices[189] = "Canon EF 70-200mm f/4L + 2.8x";
            choices[190] = "Canon EF 100mm f/2.8 Macro";
            choices[191] = "Canon EF 400mm f/4 DO IS";
            choices[193] = "Canon EF 35-80mm f/4-5.6 USM";
            choices[194] = "Canon EF 80-200mm f/4.5-5.6 USM";
            choices[195] = "Canon EF 35-105mm f/4.5-5.6 USM";
            choices[196] = "Canon EF 75-300mm f/4-5.6 USM";
            choices[197] = "Canon EF 75-300mm f/4-5.6 IS";
            choices[198] = "Canon EF 50mm f/1.4 USM";
            choices[199] = "Canon EF 28-80mm f/3.5-5.6 USM";
            choices[200] = "Canon EF 75-300mm f/4-5.6 USM";
            choices[201] = "Canon EF 28-80mm f/3.5-5.6 USM";
            choices[202] = "Canon EF 28-80 f/3.5-5.6 USM IV";
            choices[208] = "Canon EF 22-55mm f/4-5.6 USM";
            choices[209] = "Canon EF 55-200mm f/4.5-5.6";
            choices[210] = "Canon EF 28-90mm f/4-5.6 USM";
            choices[211] = "Canon EF 28-200mm f/3.5-5.6";
            choices[212] = "Canon EF 28-105mm f/4-5.6 USM";
            choices[213] = "Canon EF 90-300mm f/4.5-5.6";
            choices[214] = "Canon EF-S 18-55mm f/3.5-4.5 USM";
            choices[215] = "Canon EF 55-200mm f/4.5-5.6 II USM";
            choices[224] = "Canon EF 70-200mm f/2.8L IS USM";
            choices[225] = "Canon EF 70-200mm f/2.8L IS USM + x1.4";
            choices[226] = "Canon EF 70-200mm f/2.8L IS USM + x2";
            choices[227] = "Canon EF 70-200mm f/2.8L IS + 2.8x";
            choices[228] = "Canon EF 28-105mm f/3.5-4.5 USM";
            choices[229] = "Canon EF 16-35mm f/2.8L";
            choices[230] = "Canon EF 24-70mm f/2.8L";
            choices[231] = "Canon EF 17-40mm f/4L";
            choices[232] = "Canon EF 70-300mm f/4.5-5.6 DO IS USM";
            choices[233] = "Canon EF 28-300mm f/3.5-5.6L IS";
            choices[234] = "Canon EF-S 17-85mm f4-5.6 IS USM";
            choices[235] = "Canon EF-S10-22mm F3.5-4.5 USM";
            choices[236] = "Canon EF-S60mm F2.8 Macro USM";
            choices[237] = "Canon EF 24-105mm f/4L IS";
            choices[238] = "Canon EF 70-300mm f/4-5.6 IS USM";
            choices[239] = "Canon EF 85mm f/1.2L II USM";
            choices[240] = "Canon EF-S 17-55mm f/2.8 IS USM";
            choices[241] = "Canon EF 50mm f/1.2L USM";
            choices[242] = "Canon EF 70-200mm f/4L IS USM";
            choices[243] = "Canon EF 70-200mm f/4L IS + 1.4x";
            choices[244] = "Canon EF 70-200mm f/4L IS + 2x";
            choices[245] = "Canon EF 70-200mm f/4L IS + 2.8x";
            choices[246] = "Canon EF 16-35mm f/2.8L II";
            choices[247] = "Canon EF 14mm f/2.8L II USM";        
            choices[248] = "Canon EF 200mm f/2L IS";
            choices[249] = "Canon EF 800mm f/5.6L IS";
            choices[250] = "Canon EF 24 f/1.4L II";
            choices[254] = "Canon EF 100mm f/2.8L Macro IS USM";
            choices[488] = "Canon EF-S 15-85mm f/3.5-5.6 IS USM";
        }
        virtual std::string toString (Tag* t) {
            std::ostringstream str;
            str << "MacroMode = " << machoices[t->toInt(2,SHORT)] << std::endl;
            str << "Self-timer = " << t->toInt(4,SHORT) << std::endl;
            str << "Quality = " << qlchoices[t->toInt(6,SHORT)] << std::endl;
            str << "CanonFlashMode = " << fmchoices[t->toInt(8,SHORT)] << std::endl;
            str << "ContinuousDrive = " << cdchoices[t->toInt(10,SHORT)] << std::endl;
            str << "FocusMode = " << fochoices[t->toInt(14,SHORT)] << std::endl;
            str << "RecordMode = " << rmchoices[t->toInt(18,SHORT)] << std::endl;
            str << "CanonImageSize = " << ischoices[t->toInt(20,SHORT)] << std::endl;
            str << "EasyMode = " << emchoices[t->toInt(22,SHORT)] << std::endl;
            str << "DigitalZoom = " << dzchoices[t->toInt(24,SHORT)] << std::endl;
            str << "Contrast = " << t->toInt(26,SHORT) << std::endl;
            str << "Saturation = " << t->toInt(28,SHORT) << std::endl;
            str << "Sharpness = " << t->toInt(30,SHORT) << std::endl;
            str << "CameraISO = " << t->toInt(32,SHORT) << std::endl;
            str << "MeteringMode = " << mmchoices[t->toInt(34,SHORT)] << std::endl;
            str << "FocusRange = " << frchoices[t->toInt(36,SHORT)] << std::endl;
            str << "AFPoint = " << afchoices[t->toInt(38,SHORT)] << std::endl;
            str << "CanonExposureMode = " << exchoices[t->toInt(40,SHORT)] << std::endl;
            str << "LensType = " << choices[t->toInt(44,SHORT)] << " (" << t->toInt(44,SHORT) << ")" << std::endl;
            str << "LongFocal = " << (double)t->toInt(46,SHORT)/t->toInt(50,SHORT) << " mm" << std::endl;
            str << "ShortFocal = " << (double)t->toInt(48,SHORT)/t->toInt(50,SHORT) << " mm" << std::endl;
            str << "FocalUnits = " << t->toInt(50,SHORT) << std::endl;
            str << "MaxAperture = " << pow (2, t->toInt(52,SHORT)/64.0) << std::endl;
            str << "MinAperture = " << pow (2, t->toInt(54,SHORT)/64.0) << std::endl;
            str << "FlashActivity = " << t->toInt(56,SHORT) << std::endl;
            str << "FlashBits = ";
            int f = t->toInt(58,SHORT);
            if (f&1)
                str << "Manual ";
            if (f&2)
                str << "TTL ";
            if (f&4)
                str << "A-TTL ";
            if (f&8)
                str << "E-TTL ";
            if (f&16)
                str << "FP sync enabled ";
            if (f&(1<<7))
                str << "2nd-curtain sync used ";
            if (f&(1<<11))
                str << "FP sync used ";
            if (f&(1<<13))
                str << "Built-in ";
            if (f&(1<<14))
                str << "External ";
            str << std::endl;
            str << "FocusContinuous = " << fcchoices[t->toInt(64,SHORT)] << std::endl;
            str << "AESetting = " << aechoices[t->toInt(66,SHORT)] << std::endl;
            str << "ImageStabilization = " << stchoices[t->toInt(68,SHORT)] << " (" << t->toInt(68,SHORT) << ")" << std::endl;
            str << "DisplayAperture = " << t->toInt(70,SHORT) << std::endl;
            str << "ZoomSourceWidth = " << t->toInt(72,SHORT) << std::endl;
            str << "ZoomTargetWidth = " << t->toInt(74,SHORT) << std::endl;
            str << "SpotMeteringMode = " << smchoices[t->toInt(78,SHORT)] << std::endl;
            str << "PhotoEffect = " << pechoices[t->toInt(80,SHORT)] << std::endl;
            str << "ManualFlashOutput = " << mfchoices[t->toInt(82,SHORT)] << std::endl;
            str << "ColorTone = " << t->toInt(84,SHORT);
            return str.str();
        }
};
CACameraSettingsInterpreter caCameraSettingsInterpreter;


class CAProcessingInfoInterpreter : public Interpreter {
        std::map<int,std::string> tcchoices;
        std::map<int,std::string> sfchoices;
        std::map<int,std::string> wbchoices;
        std::map<int,std::string> pschoices;
    public:
        CAProcessingInfoInterpreter () {
            tcchoices[0] = "Standard";
            tcchoices[1] = "Manual";
            tcchoices[2] = "Custom";
            sfchoices[0] = "N/A";
            sfchoices[1] = "Lowest";
            sfchoices[2] = "Low";
            sfchoices[3] = "Standard";
            sfchoices[4] = "High";
            sfchoices[5] = "Highest";
            wbchoices[0] = "Auto";
            wbchoices[1] = "Daylight";
            wbchoices[2] = "Cloudy";
            wbchoices[3] = "Tungsten";
            wbchoices[4] = "Fluorescent";
            wbchoices[5] = "Flash";
            wbchoices[6] = "Custom";
            wbchoices[7] = "Black & White";
            wbchoices[8] = "Shade";
            wbchoices[9] = "Manual Temperature (Kelvin)";
            wbchoices[10] = "PC Set1";
            wbchoices[11] = "PC Set2";
            wbchoices[12] = "PC Set3";
            wbchoices[14] = "Daylight Fluorescent";
            wbchoices[15] = "Custom 1";
            wbchoices[16] = "Custom 2";
            wbchoices[17] = "Underwater";
            pschoices[0] = "None";
            pschoices[1] = "Standard ";
            pschoices[2] = "Set 1";
            pschoices[3] = "Set 2";
            pschoices[4] = "Set 3";
            pschoices[0x21] = "User Def. 1";
            pschoices[0x22] = "User Def. 2";
            pschoices[0x23] = "User Def. 3";
            pschoices[0x41] = "External 1";
            pschoices[0x42] = "External 2";
            pschoices[0x43] = "External 3";
            pschoices[0x81] = "Standard";
            pschoices[0x82] = "Portrait";
            pschoices[0x83] = "Landscape";
            pschoices[0x84] = "Neutral";
            pschoices[0x85] = "Faithful";
            pschoices[0x86] = "Monochrome";
        }
        virtual std::string toString (Tag* t) {
            std::ostringstream str;
            str << "ToneCurve = " << tcchoices[t->toInt(2,SHORT)] << std::endl;
            str << "Sharpness = " << t->toInt(4,SHORT) << std::endl;
            str << "SharpnessFrequency = " << sfchoices[t->toInt(6,SHORT)] << std::endl;
            str << "SensorRedLevel = " << t->toInt(8,SHORT) << std::endl;
            str << "SensorBlueLevel = " << t->toInt(10,SHORT) << std::endl;
            str << "WhiteBalanceRed = " << t->toInt(12,SHORT) << std::endl;
            str << "WhiteBalanceBlue = " << t->toInt(14,SHORT) << std::endl;
            str << "WhiteBalance = " << wbchoices[t->toInt(16,SHORT)] << std::endl;
            str << "ColorTemperature = " << t->toInt(18,SHORT) << std::endl;
            str << "PictureStyle = " << pschoices[t->toInt(20,SHORT)] << std::endl;
            str << "DigitalGain = " << t->toInt(22,SHORT) << std::endl;
            str << "WBShiftAB = " << t->toInt(24,SHORT) << std::endl;
            str << "WBShiftGM = " << t->toInt(26,SHORT);
            return str.str();
        }
};
CAProcessingInfoInterpreter caProcessingInfoInterpreter;

class CAShotInfoInterpreter : public Interpreter {
        std::map<short,std::string> sschoices;
        std::map<short,std::string> afchoices;
        std::map<short,std::string> aechoices;
        std::map<short,std::string> wbchoices;
        std::map<short,std::string> ctchoices;
        std::map<short,std::string> cmchoices;
        std::map<short,std::string> archoices;
        std::map<short,std::string> ndchoices;
    public:
        CAShotInfoInterpreter () {
            sschoices[0] = "Off";
            sschoices[1] = "Night Scene";
            sschoices[2] = "On";
            sschoices[3] = "None";
            afchoices[0x3000] = "None (MF)";
            afchoices[0x3001] = "Right";
            afchoices[0x3002] = "Center";
            afchoices[0x3003] = "Center+Right";
            afchoices[0x3004] = "Left";
            afchoices[0x3005] = "Left+Right";
            afchoices[0x3006] = "Left+Center";
            afchoices[0x3007] = "All";
            wbchoices[0] = "Auto";
            wbchoices[1] = "Daylight";
            wbchoices[2] = "Cloudy";
            wbchoices[3] = "Tungsten";
            wbchoices[4] = "Fluorescent";
            wbchoices[5] = "Flash";
            wbchoices[6] = "Custom";
            wbchoices[7] = "Black & White";
            wbchoices[8] = "Shade";
            wbchoices[9] = "Manual Temperature (Kelvin)";
            wbchoices[10] = "PC Set1";
            wbchoices[11] = "PC Set2";
            wbchoices[12] = "PC Set3";
            wbchoices[14] = "Daylight Fluorescent";
            wbchoices[15] = "Custom 1";
            wbchoices[16] = "Custom 2";
            wbchoices[17] = "Underwater";
            aechoices[-1] = "On ";
            aechoices[0] = "Off ";
            aechoices[1] = "On (shot 1)";
            aechoices[2] = "On (shot 2)";
            aechoices[3] = "On (shot 3)";
            cmchoices[0] = "n/a";
            cmchoices[1] = "Camera Local Control";
            cmchoices[3] = "Computer Remote Control";
            ctchoices[248] = "EOS High-end";
            ctchoices[250] = "Compact";
            ctchoices[252] = "EOS Mid-end";
            ctchoices[255] = "DV Camera";
            ctchoices[0x23] = "User Def. 3";
            archoices[-1] = "Rotated by Software";
            archoices[0] = "None";
            archoices[1] = "Rotate 90 CW";
            archoices[2] = "Rotate 180";
            archoices[3] = "Rotate 270 CW";
            ndchoices[0] = "Off";
            ndchoices[1] = "On";
        }
        virtual std::string toString (Tag* t) {
            std::ostringstream str;
            str << "AutoISO = " << t->toInt(2,SHORT) << std::endl;
            str << "BaseISO = " << pow (2, t->toInt(4,SHORT)/32.0 - 4) * 50 << std::endl;
            str << "MeasuredEV = " << t->toInt(6,SHORT) << std::endl;
            str << "TargetAperture = " << pow (2, t->toInt(8,SHORT)/64.0) << std::endl;
            str << "TargetExposureTime = " << pow (2, -t->toInt(10,SHORT)/32.0) << std::endl;
            str << "ExposureCompensation = " << t->toInt(12,SHORT)/32.0 << std::endl;
            str << "WhiteBalance = " << wbchoices[t->toInt(14,SHORT)] << std::endl;
            str << "SlowShutter = " << sschoices[t->toInt(16,SHORT)] << std::endl;
            str << "SequenceNumber = " << t->toInt(18,SHORT) << std::endl;
            str << "OpticalZoomCode = " << t->toInt(20,SHORT) << std::endl;
            str << "FlashGuideNumber = " << t->toInt(26,SHORT) << std::endl;
            str << "AFPointsInFocus = " << afchoices[t->toInt(28,SHORT)] << std::endl;
            str << "FlashExposureComp = " << t->toInt(30,SHORT) << std::endl;
            str << "AutoExposureBracketing = " << afchoices[t->toInt(32,SHORT)] << std::endl;
            str << "AEBBracketValue = " << t->toInt(34,SHORT) << std::endl;
            str << "ControlMode = " << cmchoices[t->toInt(36,SHORT)] << std::endl;
            str << "FocusDistanceUpper = " << t->toInt(38,SHORT) << std::endl;
            str << "FocusDistanceLower = " << t->toInt(40,SHORT) << std::endl;
            str << "FNumber = " << pow (2, t->toInt(42,SHORT)/64.0) << std::endl;
            str << "ExposureTime = " << pow (2, -t->toInt(44,SHORT)/32.0) << std::endl;
            str << "BulbDuration = " << t->toInt(48,SHORT) << std::endl;
            str << "CameraType = " << ctchoices[t->toInt(52,SHORT)] << std::endl;
            str << "AutoRotate = " << archoices[t->toInt(54,SHORT)] << std::endl;
            str << "NDFilter = " << ndchoices[t->toInt(56,SHORT)] << std::endl;
            str << "Self-timer2 = " << t->toInt(58,SHORT) << std::endl;
            str << "FlashOutput = " << t->toInt(66,SHORT);
            return str.str();
        }
};
CAShotInfoInterpreter caShotInfoInterpreter;


class CAFileInfoInterpreter : public Interpreter {
        std::map<int,std::string> bmchoices;
        std::map<int,std::string> rjqchoices;
        std::map<int,std::string> rjschoices;
        std::map<int,std::string> nrchoices;
        std::map<int,std::string> wbchoices;
        std::map<int,std::string> fechoices;
        std::map<int,std::string> techoices;
    public:
        CAFileInfoInterpreter () {
            bmchoices[0] = "Off";
            bmchoices[1] = "AEB";
            bmchoices[2] = "FEB";
            bmchoices[3] = "ISO";
            bmchoices[4] = "WB";
            
            rjqchoices[1] = "Economy";
            rjqchoices[2] = "Normal";
            rjqchoices[3] = "Fine";
            rjqchoices[4] = "RAW";
            rjqchoices[5] = "Superfine";

            rjschoices[0] = "Large";
            rjschoices[1] = "Medium";
            rjschoices[2] = "Small";
            rjschoices[5] = "Medium 1";
            rjschoices[6] = "Medium 2";
            rjschoices[7] = "Medium 3";
            rjschoices[8] = "Postcard";
            rjschoices[9] = "Widescreen";
            
            nrchoices[0] = "Off";
            nrchoices[1] = "On (mode 1)";
            nrchoices[2] = "On (mode 2)";
            nrchoices[3] = "On (mode 3)";
            nrchoices[4] = "On (mode 4)";

            wbchoices[0] = "Off";
            wbchoices[1] = "On (shift AB)";
            wbchoices[2] = "On (shift GM)";

            fechoices[0] = "None";
            fechoices[1] = "Yellow";
            fechoices[2] = "Orange";
            fechoices[3] = "Red";
            fechoices[4] = "Green";

            techoices[0] = "None";
            techoices[1] = "Sepia";
            techoices[2] = "Blue";
            techoices[3] = "Purple";
            techoices[4] = "Green";

        }
        virtual std::string toString (Tag* t) {

            std::ostringstream str;
            str << "FileNumber  = " << t->toInt(1,SHORT) << std::endl;
            str << "ShutterCount = " << t->toInt(0,LONG) << std::endl;
            str << "BracketMode = " << bmchoices[t->toInt(6,SHORT)] << std::endl;
            str << "BracketValue = " << t->toInt(8,SHORT) << std::endl;
            str << "BracketShotNumber = " << t->toInt(10,SHORT) << std::endl;
            str << "RawJpgQuality = " << rjqchoices[t->toInt(12,SHORT)] << std::endl;
            str << "RawJpgSize = " << rjschoices[t->toInt(14,SHORT)] << std::endl;
            str << "NoiseReduction = " << nrchoices[t->toInt(16,SHORT)] << std::endl;
            str << "WBBracketMode = " << t->toInt(18,SHORT) << std::endl;
            str << "WBBracketValueAB = " << t->toInt(24,SHORT) << std::endl;
            str << "FilterEffect = " << fechoices[t->toInt(26,SHORT)] << std::endl;
            str << "ToningEffect = " << techoices[t->toInt(30,SHORT)];
            return str.str();
        }
};
CAFileInfoInterpreter caFileInfoInterpreter;

class CAModelIDInterpreter : public ChoiceInterpreter {
    public:
        CAModelIDInterpreter () {
            choices[0x1010000] = "PowerShot A30";
            choices[0x1040000] = "PowerShot S300 / Digital IXUS 300 / IXY Digital 300";
            choices[0x1060000] = "PowerShot A20";
            choices[0x1080000] = "PowerShot A10";
            choices[0x1090000] = "PowerShot S110 / Digital IXUS v / IXY Digital 200";
            choices[0x1100000] = "PowerShot G2";
            choices[0x1110000] = "PowerShot S40";
            choices[0x1120000] = "PowerShot S30";
            choices[0x1130000] = "PowerShot A40";
            choices[0x1140000] = "EOS D30";
            choices[0x1150000] = "PowerShot A100";
            choices[0x1160000] = "PowerShot S200 / Digital IXUS v2 / IXY Digital 200a";
            choices[0x1170000] = "PowerShot A200";
            choices[0x1180000] = "PowerShot S330 / Digital IXUS 330 / IXY Digital 300a";
            choices[0x1190000] = "PowerShot G3";
            choices[0x1210000] = "PowerShot S45";
            choices[0x1230000] = "PowerShot SD100 / Digital IXUS II / IXY Digital 30";
            choices[0x1240000] = "PowerShot S230 / Digital IXUS v3 / IXY Digital 320";
            choices[0x1250000] = "PowerShot A70";
            choices[0x1260000] = "PowerShot A60";
            choices[0x1270000] = "PowerShot S400 / Digital IXUS 400 / IXY Digital 400";
            choices[0x1290000] = "PowerShot G5";
            choices[0x1300000] = "PowerShot A300";
            choices[0x1310000] = "PowerShot S50";
            choices[0x1340000] = "PowerShot A80";
            choices[0x1350000] = "PowerShot SD10 / Digital IXUS i / IXY Digital L";
            choices[0x1360000] = "PowerShot S1 IS";
            choices[0x1370000] = "PowerShot Pro1";
            choices[0x1380000] = "PowerShot S70";
            choices[0x1390000] = "PowerShot S60";
            choices[0x1400000] = "PowerShot G6";
            choices[0x1410000] = "PowerShot S500 / Digital IXUS 500 / IXY Digital 500";
            choices[0x1420000] = "PowerShot A75";
            choices[0x1440000] = "PowerShot SD110 / Digital IXUS IIs / IXY Digital 30a";
            choices[0x1450000] = "PowerShot A400";
            choices[0x1470000] = "PowerShot A310";
            choices[0x1490000] = "PowerShot A85";
            choices[0x1520000] = "PowerShot S410 / Digital IXUS 430 / IXY Digital 450";
            choices[0x1530000] = "PowerShot A95";
            choices[0x1540000] = "PowerShot SD300 / Digital IXUS 40 / IXY Digital 50";
            choices[0x1550000] = "PowerShot SD200 / Digital IXUS 30 / IXY Digital 40";
            choices[0x1560000] = "PowerShot A520";
            choices[0x1570000] = "PowerShot A510";
            choices[0x1590000] = "PowerShot SD20 / Digital IXUS i5 / IXY Digital L2";
            choices[0x1640000] = "PowerShot S2 IS";
            choices[0x1650000] = "PowerShot SD430 / IXUS Wireless / IXY Wireless";
            choices[0x1660000] = "PowerShot SD500 / Digital IXUS 700 / IXY Digital 600";
            choices[0x1668000] = "EOS D60";
            choices[0x1700000] = "PowerShot SD30 / Digital IXUS i zoom / IXY Digital L3";
            choices[0x1740000] = "PowerShot A430";
            choices[0x1750000] = "PowerShot A410";
            choices[0x1760000] = "PowerShot S80";
            choices[0x1780000] = "PowerShot A620";
            choices[0x1790000] = "PowerShot A610";
            choices[0x1800000] = "PowerShot SD630 / Digital IXUS 65 / IXY Digital 80";
            choices[0x1810000] = "PowerShot SD450 / Digital IXUS 55 / IXY Digital 60";
            choices[0x1820000] = "PowerShot TX1";
            choices[0x1870000] = "PowerShot SD400 / Digital IXUS 50 / IXY Digital 55";
            choices[0x1880000] = "PowerShot A420";
            choices[0x1890000] = "PowerShot SD900 / Digital IXUS 900 Ti / IXY Digital 1000";
            choices[0x1900000] = "PowerShot SD550 / Digital IXUS 750 / IXY Digital 700";
            choices[0x1920000] = "PowerShot A700";
            choices[0x1940000] = "PowerShot SD700 IS / Digital IXUS 800 IS / IXY Digital 800 IS";
            choices[0x1950000] = "PowerShot S3 IS";
            choices[0x1960000] = "PowerShot A540";
            choices[0x1970000] = "PowerShot SD600 / Digital IXUS 60 / IXY Digital 70";
            choices[0x1980000] = "PowerShot G7";
            choices[0x1990000] = "PowerShot A530";
            choices[0x2000000] = "PowerShot SD800 IS / Digital IXUS 850 IS / IXY Digital 900 IS";
            choices[0x2010000] = "PowerShot SD40 / Digital IXUS i7 / IXY Digital L4";
            choices[0x2020000] = "PowerShot A710 IS";
            choices[0x2030000] = "PowerShot A640";
            choices[0x2040000] = "PowerShot A630";
            choices[0x2090000] = "PowerShot S5 IS";
            choices[0x2100000] = "PowerShot A460";
            choices[0x2120000] = "PowerShot SD850 IS / Digital IXUS 950 IS / IXY Digital 810 IS";
            choices[0x2130000] = "PowerShot A570 IS";
            choices[0x2140000] = "PowerShot A560";
            choices[0x2150000] = "PowerShot SD750 / Digital IXUS 75 / IXY Digital 90";
            choices[0x2160000] = "PowerShot SD1000 / Digital IXUS 70 / IXY Digital 10";
            choices[0x2180000] = "PowerShot A550";
            choices[0x2190000] = "PowerShot A450";
            choices[0x2230000] = "PowerShot G9";
            choices[0x2240000] = "PowerShot A650 IS";
            choices[0x2260000] = "PowerShot A720 IS";
            choices[0x2290000] = "PowerShot SX100 IS";
            choices[0x2300000] = "PowerShot SD950 IS / Digital IXUS 960 IS / IXY Digital 2000 IS";
            choices[0x2310000] = "PowerShot SD870 IS / Digital IXUS 860 IS / IXY Digital 910 IS";
            choices[0x3010000] = "PowerShot Pro90 IS";
            choices[0x4040000] = "PowerShot G1";
            choices[0x6040000] = "PowerShot S100 / Digital IXUS / IXY Digital";
            choices[0x4007d675] = "HV10";
            choices[0x4007d777] = "iVIS DC50";
            choices[0x4007d778] = "iVIS HV20";
            choices[0x80000001] = "EOS-1D";
            choices[0x80000167] = "EOS-1DS";
            choices[0x80000168] = "EOS 10D";
            choices[0x80000169] = "EOS-1D Mark III";
            choices[0x80000170] = "EOS Digital Rebel / 300D / Kiss Digital";
            choices[0x80000174] = "EOS-1D Mark II";
            choices[0x80000175] = "EOS 20D";
            choices[0x80000188] = "EOS-1Ds Mark II";
            choices[0x80000189] = "EOS Digital Rebel XT / 350D / Kiss Digital N";
            choices[0x80000190] = "EOS 40D";
            choices[0x80000213] = "EOS 5D";
            choices[0x80000215] = "EOS-1Ds Mark III";
            choices[0x80000232] = "EOS-1D Mark II N";
            choices[0x80000234] = "EOS 30D";
            choices[0x80000236] = "EOS Digital Rebel XTi / 400D / Kiss Digital X";
            choices[0x80000254] = "EOS Rebel XS / 1000D / Kiss F";
            choices[0x80000261] = "EOS 50D";
        }
};

CAModelIDInterpreter caModelIDInterpreter;


const TagAttrib canonAttribs[] = {
 0, 1, 0, 0, 0x0001, "CanonCameraSettings", &caCameraSettingsInterpreter,
 0, 1, 0, 0, 0x0002, "CanonFocalLength", &caFocalLengthInterpreter,
 0, 1, 0, 0, 0x0003, "CanonFlashInfo", &stdInterpreter,
 0, 1, 0, 0, 0x0004, "CanonShotInfo", &caShotInfoInterpreter,
 0, 1, 0, 0, 0x0005, "CanonPanorama", &stdInterpreter,
 0, 1, 0, 0, 0x0006, "CanonImageType", &stdInterpreter,
 0, 1, 0, 0, 0x0007, "CanonFirmwareVersion", &stdInterpreter,
 0, 1, 0, 0, 0x0008, "FileNumber", &stdInterpreter,
 0, 1, 0, 0, 0x0009, "OwnerName", &stdInterpreter,
 0, 1, 0, 0, 0x000a, "ColorInfoD30", &stdInterpreter,
 0, 1, 0, 0, 0x000c, "SerialNumber", &stdInterpreter,
 0, 1, 0, 0, 0x000d, "CanonCameraInfo", &stdInterpreter,
 0, 1, 0, 0, 0x000e, "CanonFileLength", &stdInterpreter,
 0, 1, 0, 0, 0x000f, "CustomFunctions", &stdInterpreter,
 0, 1, 0, 0, 0x0010, "CanonModelID", &caModelIDInterpreter,
 0, 1, 0, 0, 0x0012, "CanonAFInfo", &stdInterpreter,
 0, 1, 0, 0, 0x0015, "SerialNumberFormat", &stdInterpreter,
 0, 1, 0, 0, 0x001c, "DateStampMode", &stdInterpreter,
 0, 1, 0, 0, 0x001d, "MyColors", &stdInterpreter,
 0, 1, 0, 0, 0x001e, "FirmwareRevision", &stdInterpreter,
 0, 3, 0, 0, 0x0024, "FaceDetect1", &stdInterpreter,
 0, 3, 0, 0, 0x0025, "FaceDetect2", &stdInterpreter,
 0, 1, 0, 0, 0x0026, "CanonAFInfo2", &stdInterpreter,
 0, 1, 0, 0, 0x0083, "OriginalDecisionData", &stdInterpreter,
 0, 1, 0, 0, 0x0090, "CustomFunctions1D", &stdInterpreter,
 0, 1, 0, 0, 0x0091, "PersonalFunctions", &stdInterpreter,
 0, 1, 0, 0, 0x0092, "PersonalFunctionValues", &stdInterpreter,
 0, 1, 0, 0, 0x0093, "CanonFileInfo", &caFileInfoInterpreter,
 0, 1, 0, 0, 0x0094, "AFPointsInFocus1D", &stdInterpreter,
 0, 1, 0, 0, 0x0095, "LensType", &stdInterpreter,
 0, 1, 0, 0, 0x0096, "InternalSerialNumber", &caIntSerNumInterpreter,
 0, 1, 0, 0, 0x0097, "DustRemovalData", &stdInterpreter,
 0, 1, 0, 0, 0x0099, "CustomFunctions2", &stdInterpreter,
 0, 1, 0, 0, 0x00a0, "ProccessingInfo", &caProcessingInfoInterpreter,
 0, 1, 0, 0, 0x00a1, "ToneCurveTable", &stdInterpreter,
 0, 1, 0, 0, 0x00a2, "SharpnessTable", &stdInterpreter,
 0, 1, 0, 0, 0x00a3, "SharpnessFreqTable", &stdInterpreter,
 0, 1, 0, 0, 0x00a4, "WhiteBalanceTable", &stdInterpreter,
 0, 1, 0, 0, 0x00a9, "ColorBalance", &stdInterpreter,
 0, 1, 0, 0, 0x00ae, "ColorTemperature", &stdInterpreter,
 0, 3, 0, 0, 0x00b0, "CanonFlags", &stdInterpreter,
 0, 1, 0, 0, 0x00b1, "ModifiedInfo", &stdInterpreter,
 0, 1, 0, 0, 0x00b2, "ToneCurveMatching", &stdInterpreter,
 0, 1, 0, 0, 0x00b3, "WhiteBalanceMatching", &stdInterpreter,
 0, 1, 0, 0, 0x00b4, "ColorSpace", &stdInterpreter,
 1, 1, 0, 0, 0x00b6, "PreviewImageInfo", &stdInterpreter,
 0, 1, 0, 0, 0x00d0, "VRDOffset", &stdInterpreter,
 0, 1, 0, 0, 0x00e0, "SensorInfo", &stdInterpreter,
 0, 1, 0, 0, 0x4001, "ColorBalance", &stdInterpreter,
 0, 1, 0, 0, 0x4002, "UnknownBlock1", &stdInterpreter,
 0, 1, 0, 0, 0x4003, "ColorInfo", &stdInterpreter,
 1, 1, 0, 0, 0x4005, "UnknownBlock2", &stdInterpreter,
 1, 1, 0, 0, 0x4008, "BlackLevel", &stdInterpreter,
 -1, 0, 0,  0, 0, "", NULL};

};
#endif

