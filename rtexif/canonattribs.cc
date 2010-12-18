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

class CAOnOffInterpreter : public Interpreter {
public:
    virtual std::string toString (Tag* t)
    {
    	int n = t->toInt();
    	if( n==0 ) return "OFF";
    	else if( n == 1) return "ON";
    	else return "undef";
    }
};
CAOnOffInterpreter caOnOffInterpreter;

class CAIntSerNumInterpreter : public Interpreter {
    public:
        CAIntSerNumInterpreter () {}
        virtual std::string toString (Tag* t) { return ""; }
};

CAIntSerNumInterpreter caIntSerNumInterpreter;

class CAApertureInterpreter : public Interpreter {
    public:
        CAApertureInterpreter () {}
        virtual std::string toString (Tag* t) {
        	char buffer[32];
            sprintf (buffer, "%0.1f", pow(2.0, t->toDouble()/64.0));
            return buffer;
        }
};
CAApertureInterpreter caApertureInterpreter;

class CAMacroModeInterpreter : public ChoiceInterpreter {
public:
	CAMacroModeInterpreter(){
        choices[1] = "Macro";
        choices[2] = "Normal";
	}
};
CAMacroModeInterpreter caMacroModeInterpreter;

class CASelfTimerInterpreter : public Interpreter {
public:
    virtual std::string toString (Tag* t) {
    	int sec = t->toInt(0,SHORT);
    	if( !sec ) return "OFF";
    	char buffer[32];
        sprintf (buffer, "%0.1fs %s", sec/10., sec&0x4000?",Custom":"");
        return buffer;
    }
};
CASelfTimerInterpreter caSelfTimerInterpreter;

class CAQualityInterpreter : public ChoiceInterpreter {
public:
	CAQualityInterpreter(){
        choices[1] = "Economy";
        choices[2] = "Normal";
        choices[3] = "Fine";
        choices[4] = "RAW";
        choices[5] = "Superfine";
	}
};
CAQualityInterpreter caQualityInterpreter;

class CAFlashModeInterpreter : public ChoiceInterpreter {
public:
	CAFlashModeInterpreter(){
        choices[0] = "Off";
        choices[1] = "Auto";
        choices[2] = "On";
        choices[3] = "Red-eye reduction";
        choices[4] = "Slow-sync";
        choices[5] = "Red-eye reduction (Auto)";
        choices[6] = "Red-eye reduction (On)";
        choices[16]= "External flash";
	}
};
CAFlashModeInterpreter caFlashModeInterpreter;

class CAContinuousDriveInterpreter : public ChoiceInterpreter {
public:
	CAContinuousDriveInterpreter(){
        choices[0] = "Single";
        choices[1] = "Continuous";
        choices[2] = "Movie";
        choices[3] = "Continuous, Speed Priority";
        choices[4] = "Continuous, Low";
        choices[5] = "Continuous, High";
	}
};
CAContinuousDriveInterpreter caContinuousDriveInterpreter;

class CAFocusModeInterpreter : public ChoiceInterpreter {
public:
	CAFocusModeInterpreter(){
        choices[0] = "One-shot AF";
        choices[1] = "AI Servo AF";
        choices[2] = "AI Focus AF";
        choices[3] = "Manual Focus";
        choices[4] = "Single";
        choices[5] = "Continuous";
        choices[6] = "Manual Focus";
        choices[16] = "Pan Focus";
	}
};
CAFocusModeInterpreter caFocusModeInterpreter;

class CARecordModeInterpreter : public ChoiceInterpreter {
public:
	CARecordModeInterpreter(){
        choices[1] = "JPEG";
        choices[2] = "CRW+THM";
        choices[3] = "AVI+THM";
        choices[4] = "TIF";
        choices[5] = "TIF+JPEG";
        choices[6] = "CR2";
        choices[7] = "CR2+JPEG";
        choices[9] = "Video";
	}
};
CARecordModeInterpreter caRecordModeInterpreter;

class CAImageSizeInterpreter : public ChoiceInterpreter {
public:
	CAImageSizeInterpreter (){
        choices[0] = "Large";
        choices[1] = "Medium";
        choices[2] = "Small";
        choices[5] = "Medium 1";
        choices[6] = "Medium 2";
        choices[7] = "Medium 3";
        choices[8] = "Postcard";
        choices[9] = "Widescreen";
        choices[10] = "Medium Widescreen";
        choices[128] = "640x480 Movie";
        choices[129] = "Medium Movie";
        choices[130] = "Small Movie";
        choices[137] = "1280x720 Movie";
        choices[142] = "1920x1080 Movie";
	}
};
CAImageSizeInterpreter caImageSizeInterpreter;

class CAEasyModeInterpreter : public ChoiceInterpreter {
public:
	CAEasyModeInterpreter (){
        choices[0] = "Full auto ";
        choices[1] = "Manual ";
        choices[2] = "Landscape ";
        choices[3] = "Fast shutter ";
        choices[4] = "Slow shutter ";
        choices[5] = "Night ";
        choices[6] = "Gray Scale ";
        choices[7] = "Sepia ";
        choices[8] = "Portrait ";
        choices[9] = "Sports ";
        choices[10] = "Macro ";
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
        choices[28] = "Still Image";
        choices[30] = "Color Accent";
        choices[31] = "Color Swap";
        choices[32] = "Aquarium";
        choices[33] = "ISO 3200";
        choices[38] = "Creative Auto";
        choices[42] = "Super Vivid";
        choices[43] = "Poster";
        choices[47] = "Fisheye";
        choices[48] = "Miniature";
        choices[261]= "Sunset";
	}
};
CAEasyModeInterpreter caEasyModeInterpreter;

class CADigitalZoomInterpreter : public ChoiceInterpreter {
public:
	CADigitalZoomInterpreter(){
        choices[0] = "None";
        choices[1] = "2x";
        choices[2] = "4x";
        choices[3] = "Other";
	}
};
CADigitalZoomInterpreter caDigitalZoomInterpreter;

class CAMeteringModeInterpreter : public ChoiceInterpreter {
public:
	CAMeteringModeInterpreter(){
        choices[0] = "Default";
        choices[1] = "Spot";
        choices[2] = "Average";
        choices[3] = "Evaluative";
        choices[4] = "Partial";
        choices[5] = "Center-weighted averaging";
	}
};
CAMeteringModeInterpreter caMeteringModeInterpreter;

class CAFocusRangeInterpreter : public ChoiceInterpreter {
public:
	CAFocusRangeInterpreter(){
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

class CAAFPointInterpreter : public ChoiceInterpreter {
public:
	CAAFPointInterpreter(){
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

class CAExposureModeInterpreter : public ChoiceInterpreter {
public:
	CAExposureModeInterpreter(){
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

class CAFlashBitsInterpreter : public Interpreter {
public:
    virtual std::string toString (Tag* t) {
    	std::ostringstream s;
    	unsigned bits = t->toInt(0,SHORT);
        if( bits & 0x0001 ) s << "Manual ";
        if( bits & 0x0002 ) s << "TTL ";
        if( bits & 0x0004 ) s << "A-TTL ";
        if( bits & 0x0008 ) s << "E-TTL ";
        if( bits & 0x0010 ) s << "FP sync enabled ";
        if( bits & 0x0080 ) s << "2nd curtain ";
        if( bits & 0x0800 ) s << "FP sync used ";
        if( bits & 0x2000 ) s << "Built-in ";
        if( bits & 0x4000 ) s << "External ";
        return s.str();
    }
};
CAFlashBitsInterpreter caFlashBitsInterpreter;

class CAFocusContinuousInterpreter : public ChoiceInterpreter {
public:
	CAFocusContinuousInterpreter(){
        choices[0] = "Single";
        choices[1] = "Continuous";
        choices[8] = "Manual";
	}
};
CAFocusContinuousInterpreter caFocusContinuousInterpreter;

class CAAESettingsInterpreter : public ChoiceInterpreter {
public:
	CAAESettingsInterpreter(){
        choices[0] = "Normal AE";
        choices[1] = "Exposure Compensation";
        choices[2] = "AE Lock";
        choices[3] = "AE Lock + Exposure Comp.";
        choices[4] = "No AE";
	}
};
CAAESettingsInterpreter caAESettingsInterpreter;

class CAStabilizationInterpreter : public ChoiceInterpreter {
public:
	CAStabilizationInterpreter(){
        choices[0] = "Off";
        choices[1] = "On";
        choices[2] = "On, Shot Only";
        choices[3] = "On, Panning";
        choices[4] = "On, Video";
	}
};
CAStabilizationInterpreter caStabilizationInterpreter;

class CASpotMeteringInterpreter : public ChoiceInterpreter {
public:
	CASpotMeteringInterpreter(){
        choices[0] = "Center";
        choices[1] = "AF Point";
	}
};
CASpotMeteringInterpreter caSpotMeteringInterpreter;

class CAPhotoEffectInterpreter : public ChoiceInterpreter {
public:
	CAPhotoEffectInterpreter(){
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

class CAManualFlashInterpreter : public ChoiceInterpreter {
public:
	CAManualFlashInterpreter(){
        choices[0] = "N/A";
        choices[0x500] = "Full";
        choices[0x502] = "Medium";
        choices[0x504] = "Low";
        choices[0x7fff] = "N/A";
	}
};
CAManualFlashInterpreter caManualFlashInterpreter;

class CARAWQualityInterpreter : public ChoiceInterpreter {
public:
	CARAWQualityInterpreter(){
        choices[0] = "N/A";
        choices[1] = "sRAW1 (mRAW)";
        choices[2] = "sRAW2 (sRAW)";
	}
};
CARAWQualityInterpreter caRAWQualityInterpreter;

class CAFocalInterpreter : public Interpreter {
public:
    virtual std::string toString (Tag* t) {
    	Tag *unitTag = t->getParent()->getRoot()->findTag("FocalUnits");
    	double unit = unitTag->toDouble();
    	char buffer[32];
        sprintf (buffer, "%0.1fmm", (unit>0. ? t->toDouble()/unit : t->toDouble()));
        return buffer;
    }
};
CAFocalInterpreter caFocalInterpreter;

class CALensInterpreter : public ChoiceInterpreter {
    public:
	    CALensInterpreter () {
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
            choices[251] = "Canon EF 70-200mm f/2.8L IS II USM";
            choices[254] = "Canon EF 100mm f/2.8L Macro IS USM";
            choices[488] = "Canon EF-S 15-85mm f/3.5-5.6 IS USM";
        }
};
CALensInterpreter caLensInterpreter;

class CAFocalTypeInterpreter : public ChoiceInterpreter {
public:
	CAFocalTypeInterpreter(){
		choices[0] = "undef";
		choices[1] = "Fixed";
		choices[2] = "Zoom";
	}
};
CAFocalTypeInterpreter caFocalTypeInterpreter;

class CAFocalPlaneInterpreter : public Interpreter {
public:
	virtual std::string toString (Tag* t) {
		int val = t->toInt();
		if( val <40 ) return "undef";
    	char buffer[32];
        sprintf (buffer, "%0.2fmm", val *25.4 / 1000);
        return buffer;
	}
};
CAFocalPlaneInterpreter caFocalPlaneInterpreter;

class CAExposureTimeInterpreter : public Interpreter {
public:
	virtual std::string toString (Tag* t) {
    	char buffer[32];
        sprintf (buffer, "%0.3f", pow (2, - t->toInt()/32.0) );
        return buffer;
	}
};
CAExposureTimeInterpreter caExposureTimeInterpreter;

class CAEVInterpreter : public Interpreter {
	virtual std::string toString (Tag* t) {
    	char buffer[32];
        sprintf (buffer, "%0.1f", t->toDouble()/32.0  );
        return buffer;
	}
};
CAEVInterpreter caEVInterpreter;

class CABaseISOInterpreter : public Interpreter {
public:
	virtual std::string toString (Tag* t) {
    	char buffer[32];
        sprintf (buffer, "%0.0f", pow (2, t->toInt()/32.0 - 4) * 50 );
        return buffer;
	}
};
CABaseISOInterpreter caBaseISOInterpreter;

class CAToneCurveInterpreter : public ChoiceInterpreter {
public:
	CAToneCurveInterpreter(){
        choices[0] = "Standard";
        choices[1] = "Manual";
        choices[2] = "Custom";
	}
};
CAToneCurveInterpreter caToneCurveInterpreter;

class CASharpnessFrequencyInterpreter : public ChoiceInterpreter {
public:
	CASharpnessFrequencyInterpreter(){
        choices[0] = "N/A";
        choices[1] = "Lowest";
        choices[2] = "Low";
        choices[3] = "Standard";
        choices[4] = "High";
        choices[5] = "Highest";
	}
};
CASharpnessFrequencyInterpreter caSharpnessFrequencyInterpreter;

class CAWhiteBalanceInterpreter : public ChoiceInterpreter {
public:
	CAWhiteBalanceInterpreter(){
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
	}
};
CAWhiteBalanceInterpreter caWhiteBalanceInterpreter;

class CAPictureStyleInterpreter : public ChoiceInterpreter {
public:
	CAPictureStyleInterpreter(){
        choices[0] = "None";
        choices[1] = "Standard ";
        choices[2] = "Set 1";
        choices[3] = "Set 2";
        choices[4] = "Set 3";
        choices[0x21] = "User Def. 1";
        choices[0x22] = "User Def. 2";
        choices[0x23] = "User Def. 3";
        choices[0x41] = "External 1";
        choices[0x42] = "External 2";
        choices[0x43] = "External 3";
        choices[0x81] = "Standard";
        choices[0x82] = "Portrait";
        choices[0x83] = "Landscape";
        choices[0x84] = "Neutral";
        choices[0x85] = "Faithful";
        choices[0x86] = "Monochrome";
	}
};
CAPictureStyleInterpreter caPictureStyleInterpreter;

class CASlowShutterInterpreter : public ChoiceInterpreter {
public:
	CASlowShutterInterpreter(){
        choices[0] = "Off";
        choices[1] = "Night Scene";
        choices[2] = "On";
        choices[3] = "None";
	}
};
CASlowShutterInterpreter caSlowShutterInterpreter;

class CAFlashGuideNumberInterpreter : public Interpreter{
public:
	virtual std::string toString (Tag* t) {
		int n= t->toInt();
		if( n==-1) return "undef";
	    char buffer[32];
        sprintf (buffer, "%0.f", n/32. );
        return buffer;
	}
};
CAFlashGuideNumberInterpreter caFlashGuideNumberInterpreter;

class CAAFPointsInFocusInterpreter : public ChoiceInterpreter {
public:
	CAAFPointsInFocusInterpreter(){
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

class CAAutoExposureBracketingInterpreter : public ChoiceInterpreter {
public:
	CAAutoExposureBracketingInterpreter(){
        choices[-1] = "On ";
        choices[0] = "Off ";
        choices[1] = "On (shot 1)";
        choices[2] = "On (shot 2)";
        choices[3] = "On (shot 3)";
	}
};
CAAutoExposureBracketingInterpreter caAutoExposureBracketingInterpreter;

class CAControModeInterpreter : public ChoiceInterpreter {
public:
	CAControModeInterpreter(){
        choices[0] = "n/a";
        choices[1] = "Camera Local Control";
        choices[3] = "Computer Remote Control";
	}
};
CAControModeInterpreter caControModeInterpreter;

class CAFocusDistanceInterpreter : public Interpreter {
public:
	virtual std::string toString (Tag* t) {
    	char buffer[32];
        sprintf (buffer, "%0.2f", t->toDouble()/100 );
        return buffer;
	}
};
CAFocusDistanceInterpreter caFocusDistanceInterpreter;

class CAMeasuredEVInterpreter : public Interpreter {
public:
	virtual std::string toString (Tag* t) {
    	char buffer[32];
        sprintf (buffer, "%0.1f", t->toDouble()/8 - 6 );
        return buffer;
	}
};
CAMeasuredEVInterpreter caMeasuredEVInterpreter;

class CACameraTypeInterpreter : public ChoiceInterpreter {
public:
	CACameraTypeInterpreter(){
		choices[248] = "EOS High-end";
		choices[250] = "Compact";
		choices[252] = "EOS Mid-end";
		choices[255] = "DV Camera";
	}
};
CACameraTypeInterpreter caCameraTypeInterpreter;

class CAAutoRotateInterpreter : public ChoiceInterpreter {
public:
	CAAutoRotateInterpreter(){
        choices[-1] = "Rotated by Software";
        choices[0] = "None";
        choices[1] = "Rotate 90 CW";
        choices[2] = "Rotate 180";
        choices[3] = "Rotate 270 CW";
	}
};
CAAutoRotateInterpreter caAutoRotateInterpreter;

class CABracketModeInterpreter : public ChoiceInterpreter {
public:
	CABracketModeInterpreter(){
        choices[0] = "Off";
        choices[1] = "AEB";
        choices[2] = "FEB";
        choices[3] = "ISO";
        choices[4] = "WB";
	}
};
CABracketModeInterpreter caBracketModeInterpreter;

class CARAWJpegQualityInterpreter : public ChoiceInterpreter {
public:
	CARAWJpegQualityInterpreter(){
        choices[1] = "Economy";
        choices[2] = "Normal";
        choices[3] = "Fine";
        choices[4] = "RAW";
        choices[5] = "Superfine";
	}
};
CARAWJpegQualityInterpreter caRAWJpegQualityInterpreter;

class CAJpegSizeInterpreter : public ChoiceInterpreter {
public:
	CAJpegSizeInterpreter(){
        choices[0] = "Large";
        choices[1] = "Medium";
        choices[2] = "Small";
        choices[5] = "Medium 1";
        choices[6] = "Medium 2";
        choices[7] = "Medium 3";
        choices[8] = "Postcard";
        choices[9] = "Widescreen";
	}
};
CAJpegSizeInterpreter caJpegSizeInterpreter;

class CAWBBracketModeInterpreter : public ChoiceInterpreter {
public:
	CAWBBracketModeInterpreter(){
	    choices[0] = "Off";
	    choices[1] = "Shift AB";
	    choices[2] = "shift GM";
	}
};
CAWBBracketModeInterpreter caWBBracketModeInterpreter;

class CAFilterEffectInterpreter : public ChoiceInterpreter {
public:
	CAFilterEffectInterpreter(){
        choices[0] = "None";
        choices[1] = "Yellow";
        choices[2] = "Orange";
        choices[3] = "Red";
        choices[4] = "Green";
	}
};
CAFilterEffectInterpreter caFilterEffectInterpreter;

class CAToningEffectInterpreter : public ChoiceInterpreter {
public:
	CAToningEffectInterpreter(){
        choices[0] = "None";
        choices[1] = "Sepia";
        choices[2] = "Blue";
        choices[3] = "Purple";
        choices[4] = "Green";
	}
};
CAToningEffectInterpreter caToningEffectInterpreter;

class CAFileNumberInterpreter : public Interpreter {
public:
	virtual std::string toString (Tag* t) {
	    unsigned long val = t->toInt(0,LONG);
    	char buffer[32];
        sprintf (buffer, "%ld", ((val&0xffc0)>>6)*10000+((val>>16)&0xff)+((val&0x3f)<<8) );
        return buffer;
	}
};
CAFileNumberInterpreter caFileNumberInterpreter;

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
            choices[0x2320000] = "PowerShot SD890 IS / Digital IXUS 970 IS / IXY Digital 820 IS";
            choices[0x2360000] = "PowerShot SD790 IS / Digital IXUS 90 IS / IXY Digital 95 IS";
            choices[0x2370000] = "PowerShot SD770 IS / Digital IXUS 85 IS / IXY Digital 25 IS";
            choices[0x2380000] = "PowerShot A590 IS";
            choices[0x2390000] = "PowerShot A580";
            choices[0x2420000] = "PowerShot A470";
            choices[0x2430000] = "PowerShot SD1100 IS / Digital IXUS 80 IS / IXY Digital 20 IS";
            choices[0x2460000] = "PowerShot SX1 IS";
            choices[0x2470000] = "PowerShot SX10 IS";
            choices[0x2480000] = "PowerShot A1000 IS";
            choices[0x2490000] = "PowerShot G10";
            choices[0x2510000] = "PowerShot A2000 IS";
            choices[0x2520000] = "PowerShot SX110 IS";
            choices[0x2530000] = "PowerShot SD990 IS / Digital IXUS 980 IS / IXY Digital 3000 IS";
            choices[0x2540000] = "PowerShot SD880 IS / Digital IXUS 870 IS / IXY Digital 920 IS";
            choices[0x2550000] = "PowerShot E1";
            choices[0x2560000] = "PowerShot D10";
            choices[0x2570000] = "PowerShot SD960 IS / Digital IXUS 110 IS / IXY Digital 510 IS";
            choices[0x2580000] = "PowerShot A2100 IS";
            choices[0x2590000] = "PowerShot A480";
            choices[0x2600000] = "PowerShot SX200 IS";
            choices[0x2610000] = "PowerShot SD970 IS / Digital IXUS 990 IS / IXY Digital 830 IS";
            choices[0x2620000] = "PowerShot SD780 IS / Digital IXUS 100 IS / IXY Digital 210 IS";
            choices[0x2630000] = "PowerShot A1100 IS";
            choices[0x2640000] = "PowerShot SD1200 IS / Digital IXUS 95 IS / IXY Digital 110 IS";
            choices[0x2700000] = "PowerShot G11";
            choices[0x2710000] = "PowerShot SX120 IS";
            choices[0x2720000] = "PowerShot S90";
            choices[0x2750000] = "PowerShot SX20 IS";
            choices[0x2760000] = "PowerShot SD980 IS / Digital IXUS 200 IS / IXY Digital 930 IS";
            choices[0x2770000] = "PowerShot SD940 IS / Digital IXUS 120 IS / IXY Digital 220 IS";
            choices[0x2800000] = "PowerShot A495";
            choices[0x2810000] = "PowerShot A490";
            choices[0x2820000] = "PowerShot A3100 IS";
            choices[0x2830000] = "PowerShot A3000 IS";
            choices[0x2840000] = "PowerShot SD1400 IS / IXUS 130 / IXY 400F";
            choices[0x2850000] = "PowerShot SD1300 IS / IXUS 105 / IXY 200F";
            choices[0x2860000] = "PowerShot SD3500 IS / IXUS 210 / IXY 10S";
            choices[0x2870000] = "PowerShot SX210 IS";
            choices[0x2880000] = "PowerShot SD4000 IS / IXUS 300 HS / IXY 30S";
            choices[0x2890000] = "PowerShot SD4500 IS / IXUS 1000 HS / IXY 50S";
            choices[0x2920000] = "PowerShot G12";
            choices[0x2930000] = "PowerShot SX30 IS";
            choices[0x2940000] = "PowerShot SX130 IS";
            choices[0x2950000] = "PowerShot S95";
            choices[0x3010000] = "PowerShot Pro90 IS";
            choices[0x4040000] = "PowerShot G1";
            choices[0x6040000] = "PowerShot S100 / Digital IXUS / IXY Digital";
            choices[0x4007d673]	= "DC19 / DC21 / DC22";
            choices[0x4007d674]	= "XH A1";
            choices[0x4007d675]	= "HV10";
            choices[0x4007d676]	= "MD130 / MD140 / MD150 / MD160";
            choices[0x4007d777]	= "iVIS DC50";
            choices[0x4007d778]	= "iVIS HV20";
            choices[0x4007d779]	= "DC211";
            choices[0x4007d77a]	= "HG10";
            choices[0x4007d77b]	= "iVIS HR10";
            choices[0x4007d878]	= "HV30";
            choices[0x4007d87e]	= "DC301 / DC310 / DC311 / DC320 / DC330";
            choices[0x4007d87f]	= "FS100";
            choices[0x4007d880]	= "iVIS HF10";
            choices[0x4007d882]	= "HG20 / HG21 / VIXIA HG21";
            choices[0x4007d925]	= "LEGRIA HF21";
            choices[0x4007d978]	= "LEGRIA HV40";
            choices[0x4007d987]	= "DC410 / DC420";
            choices[0x4007d988]	= "LEGRIA FS19 / FS20 / FS21 / FS22 / FS200";
            choices[0x4007d989]	= "LEGRIA HF20 / HF200";
            choices[0x4007d98a]	= "VIXIA HF S10 / S100";
            choices[0x4007da8e]	= "LEGRIA HF R16 / R17 / R106";
            choices[0x4007da90]	= "LEGRIA HF S21 / VIXIA HF S200";
            choices[0x4007da92]	= "LEGRIA FS36 / FS305 / FS306 / FS307";
            choices[0x80000001] = "EOS-1D";
            choices[0x80000167] = "EOS-1DS";
            choices[0x80000168] = "EOS 10D";
            choices[0x80000169] = "EOS-1D Mark III";
            choices[0x80000170] = "EOS Digital Rebel / 300D / Kiss Digital";
            choices[0x80000174] = "EOS-1D Mark II";
            choices[0x80000175] = "EOS 20D";
            choices[0x80000176] = "EOS Digital Rebel XSi / 450D / Kiss X2";
            choices[0x80000188] = "EOS-1Ds Mark II";
            choices[0x80000189] = "EOS Digital Rebel XT / 350D / Kiss Digital N";
            choices[0x80000190] = "EOS 40D";
            choices[0x80000213] = "EOS 5D";
            choices[0x80000215] = "EOS-1Ds Mark III";
            choices[0x80000218] = "EOS 5D Mark II";
            choices[0x80000232] = "EOS-1D Mark II N";
            choices[0x80000234] = "EOS 30D";
            choices[0x80000236] = "EOS Digital Rebel XTi / 400D / Kiss Digital X";
            choices[0x80000250] = "EOS 7D";
            choices[0x80000252] = "EOS Rebel T1i / 500D / Kiss X3";
            choices[0x80000254] = "EOS Rebel XS / 1000D / Kiss F";
            choices[0x80000261] = "EOS 50D";
            choices[0x80000270] = "EOS Rebel T2i / 550D / Kiss X4";
            choices[0x80000281] = "EOS-1D Mark IV";
            choices[0x80000287] = "EOS 60D";
        }
};
CAModelIDInterpreter caModelIDInterpreter;

class CAPanoramaDirectionInterpreter : public ChoiceInterpreter {
public:
	CAPanoramaDirectionInterpreter(){
		choices[0] = "Left to Right";
		choices[1] = "Right to Left";
		choices[2] = "Bottom to Top";
		choices[3] = "Top to Bottom";
		choices[4] = "2x2 Matrix (Clockwise)";
	}
};
CAPanoramaDirectionInterpreter caPanoramaDirectionInterpreter;

class CAAspectRatioInterpreter : public ChoiceInterpreter {
public:
	CAAspectRatioInterpreter(){
		choices[0] = "3:2";
		choices[1] = "1:1";
		choices[2] = "4:3";
		choices[7] = "16:9";
	}

};
CAAspectRatioInterpreter caAspectRatioInterpreter;

const TagAttrib canonCameraSettingsAttribs[] = {
 {0, 1, 0, 0,  1, "MacroMode", &caMacroModeInterpreter},
 {0, 1, 0, 0,  2, "SelfTimer", &caSelfTimerInterpreter},
 {0, 1, 0, 0,  3, "Quality", &caQualityInterpreter},
 {0, 1, 0, 0,  4, "CanonFlashMode", &caFlashModeInterpreter},
 {0, 1, 0, 0,  5, "ContinuousDrive", &caContinuousDriveInterpreter},
 {0, 1, 0, 0,  7, "FocusMode", &caFocusModeInterpreter},
 {0, 1, 0, 0,  9, "RecordMode", &caRecordModeInterpreter},
 {0, 1, 0, 0, 10, "CanonImageSize", &caImageSizeInterpreter},
 {0, 1, 0, 0, 11, "EasyMode", &caEasyModeInterpreter},
 {0, 1, 0, 0, 12, "DigitalZoom", &caDigitalZoomInterpreter},
 {0, 1, 0, 0, 13, "Contrast", &stdInterpreter},
 {0, 1, 0, 0, 14, "Saturation", &stdInterpreter},
 {0, 1, 0, 0, 15, "Sharpness", &stdInterpreter},
 {0, 1, 0, 0, 16, "CameraISO", &stdInterpreter},
 {0, 1, 0, 0, 17, "MeteringMode", &caMeteringModeInterpreter},
 {0, 1, 0, 0, 18, "FocusRange", &caFocusRangeInterpreter},
 {0, 1, 0, 0, 19, "AFPoint", &caAFPointInterpreter},
 {0, 1, 0, 0, 20, "CanonExposureMode", &caExposureModeInterpreter},
 {0, 1, 0, 0, 22, "LensID", &caLensInterpreter},
 {0, 1, 0, 0, 23, "LongFocal", &caFocalInterpreter},
 {0, 1, 0, 0, 24, "ShortFocal", &caFocalInterpreter},
 {0, 1, 0, 0, 25, "FocalUnits", &stdInterpreter},
 {0, 1, 0, 0, 26, "MaxAperture", &caApertureInterpreter},
 {0, 1, 0, 0, 27, "MinAperture", &caApertureInterpreter},
 {0, 1, 0, 0, 28, "FlashActivity", &stdInterpreter},
 {0, 1, 0, 0, 29, "FlashBits", &caFlashBitsInterpreter},
 {0, 1, 0, 0, 32, "FocusContinuous", &caFocusContinuousInterpreter},
 {0, 1, 0, 0, 33, "AESetting", &caAESettingsInterpreter},
 {0, 1, 0, 0, 34, "ImageStabilization", &caStabilizationInterpreter},
 {0, 1, 0, 0, 35, "DisplayAperture", &stdInterpreter},
 {0, 1, 0, 0, 36, "ZoomSourceWidth", &stdInterpreter},
 {0, 1, 0, 0, 37, "ZoomTargetWidth", &stdInterpreter},
 {0, 1, 0, 0, 39, "SpotMeteringMode", &caSpotMeteringInterpreter},
 {0, 1, 0, 0, 40, "PhotoEffect", &caPhotoEffectInterpreter},
 {0, 1, 0, 0, 41, "ManualFlashOutput", &caManualFlashInterpreter},
 {0, 1, 0, 0, 42, "ColorTone", &stdInterpreter},
 {0, 1, 0, 0, 46, "SRAWQuality", &caRAWQualityInterpreter},
 {-1, 0, 0,  0, 0, "", NULL}
};

const TagAttrib canonFocalLengthAttribs[] = {
 {0, 1, 0, 0, 0, "FocalType", &caFocalTypeInterpreter},
 {0, 1, 0, 0, 1, "FocalLength", &caFocalInterpreter},
 {0, 1, 0, 0, 2, "FocalPlaneXSize", &caFocalPlaneInterpreter},
 {0, 1, 0, 0, 3, "FocalPlaneYSize", &caFocalPlaneInterpreter},
 {-1, 0, 0,  0, 0, "", NULL}
};

const TagAttrib canonShotInfoAttribs[] = {
 {0, 1, 0, 0, 1, "AutoISO", &stdInterpreter},
 {0, 1, 0, 0, 2, "BaseISO" ,&caBaseISOInterpreter},
 {0, 1, 0, 0, 3, "MeasuredEV", &stdInterpreter},
 {0, 1, 0, 0, 4, "TargetAperture", &caApertureInterpreter},
 {0, 1, 0, 0, 5, "TargetExposureTime",&caExposureTimeInterpreter},
 {0, 1, 0, 0, 6, "ExposureCompensation",&caEVInterpreter},
 {0, 1, 0, 0, 7, "WhiteBalance",&caWhiteBalanceInterpreter},
 {0, 1, 0, 0, 8, "SlowShutter",&caSlowShutterInterpreter},
 {0, 1, 0, 0, 9, "SequenceNumber", &stdInterpreter},
 {0, 1, 0, 0,10, "OpticalZoomCode", &stdInterpreter},
 {0, 1, 0, 0,13, "FlashGuideNumber" , &caFlashGuideNumberInterpreter},
 {0, 1, 0, 0,14, "AFPointsInFocus", &caAFPointsInFocusInterpreter},
 {0, 1, 0, 0,15, "FlashExposureComp", &stdInterpreter},
 {0, 1, 0, 0,16, "AutoExposureBracketing",&caAutoExposureBracketingInterpreter},
 {0, 1, 0, 0,17, "AEBBracketValue", &stdInterpreter},
 {0, 1, 0, 0,18, "ControlMode", &caControModeInterpreter},
 {0, 1, 0, 0,19, "FocusDistanceUpper", &caFocusDistanceInterpreter},
 {0, 1, 0, 0,20, "FocusDistanceLower", &caFocusDistanceInterpreter},
 {0, 1, 0, 0,21, "FNumber" ,&caApertureInterpreter},
 {0, 1, 0, 0,22, "ExposureTime",&caExposureTimeInterpreter},
 {0, 1, 0, 0,24, "BulbDuration", &stdInterpreter},
 {0, 1, 0, 0,24, "MeasuredEV2", &caMeasuredEVInterpreter},
 {0, 1, 0, 0,26, "CameraType", &caCameraTypeInterpreter},
 {0, 1, 0, 0,27, "AutoRotate",&caAutoRotateInterpreter},
 {0, 1, 0, 0,28, "NDFilter",&caOnOffInterpreter},
 {0, 1, 0, 0,29, "Self-timer2", &stdInterpreter},
 {0, 1, 0, 0,33, "FlashOutput", &stdInterpreter},
 {-1, 0, 0, 0, 0, "", NULL},
};

const TagAttrib canonFileInfoAttribs[] = {
 {0, 1, 0, 0, 1, "FileNumber",  &caFileNumberInterpreter},
 {0, 1, 0, 0, 3, "BracketMode", &caBracketModeInterpreter},
 {0, 1, 0, 0, 4, "BracketValue", &stdInterpreter},
 {0, 1, 0, 0, 5, "BracketShotNumber", &stdInterpreter},
 {0, 1, 0, 0, 6, "RawJpgQuality",&caRAWJpegQualityInterpreter},
 {0, 1, 0, 0, 7, "RawJpgSize",&caJpegSizeInterpreter},
 {0, 1, 0, 0, 8, "NoiseReduction",&stdInterpreter},
 {0, 1, 0, 0, 9, "WBBracketMode" ,&caWBBracketModeInterpreter},
 {0, 1, 0, 0,12, "WBBracketValueAB", &stdInterpreter},
 {0, 1, 0, 0,13, "WBBracketValueGM", &stdInterpreter},
 {0, 1, 0, 0,14, "FilterEffect" ,&caFilterEffectInterpreter},
 {0, 1, 0, 0,15, "ToningEffect" ,&caToningEffectInterpreter},
 {0, 1, 0, 0,19, "LiveViewShooting" ,&caOnOffInterpreter},
 {0, 1, 0, 0,25, "FlashExposureLock" ,&caOnOffInterpreter},
 {-1,0, 0, 0, 0, "", NULL},
};

const TagAttrib canonProcessingInfoAttribs[] = {
 {0, 1, 0, 0, 1,"ToneCurve", &caToneCurveInterpreter},
 {0, 1, 0, 0, 2,"Sharpness", &stdInterpreter},
 {0, 1, 0, 0, 3,"SharpnessFrequency", &caSharpnessFrequencyInterpreter},
 {0, 1, 0, 0, 4,"SensorRedLevel", &stdInterpreter},
 {0, 1, 0, 0, 5,"SensorBlueLevel", &stdInterpreter},
 {0, 1, 0, 0, 6,"WhiteBalanceRed", &stdInterpreter},
 {0, 1, 0, 0, 7,"WhiteBalanceBlue", &stdInterpreter},
 {0, 1, 0, 0, 8,"WhiteBalance", &caWhiteBalanceInterpreter},
 {0, 1, 0, 0, 9,"ColorTemperature", &stdInterpreter},
 {0, 1, 0, 0,10,"PictureStyle", &caPictureStyleInterpreter},
 {0, 1, 0, 0,11,"DigitalGain", &stdInterpreter},
 {0, 1, 0, 0,12,"WBShiftAB", &stdInterpreter},
 {0, 1, 0, 0,13,"WBShiftGM", &stdInterpreter},
 {-1,0, 0, 0, 0, "", NULL},
};

const TagAttrib canonPanoramaInfoAttribs[] = {
 {0, 1, 0, 0, 2,"PanoramaFrameNumber", &stdInterpreter},
 {0, 1, 0, 0, 5,"PanoramaDirection", &caPanoramaDirectionInterpreter},
 {-1,0, 0, 0, 0, "", NULL},
};

const TagAttrib canonCropInfoAttribs[] = {
 {0, 1, 0, 0, 0,"CropLeftMargin", &stdInterpreter},
 {0, 1, 0, 0, 1,"CropRightMargin", &stdInterpreter},
 {0, 1, 0, 0, 2,"CropTopMargin", &stdInterpreter},
 {0, 1, 0, 0, 3,"CropBottomMargin", &stdInterpreter},
 {-1,0, 0, 0, 0, "", NULL},
};

const TagAttrib canonAspectInfoAttribs[] = {
 {0, 1, 0, 0, 0,"AspectRatio", &caAspectRatioInterpreter},
 {0, 1, 0, 0, 1,"CroppedImageWidth", &stdInterpreter},
 {0, 1, 0, 0, 2,"CroppedImageHeight", &stdInterpreter},
 {-1,0, 0, 0, 0, "", NULL},
};

const TagAttrib canonMicroAdjustAttrib[] = {
 {0, 1, 0, 0, 1,"AFMicroAdjActive", &caOnOffInterpreter},
 {-1,0, 0, 0, 2,"AFMicroAdjValue", &stdInterpreter},
};

const TagAttrib canonAttribs[] = {
 {0, 1, 0, canonCameraSettingsAttribs, 0x0001, "CanonCameraSettings", &stdInterpreter},
 {0, 1, 0, canonFocalLengthAttribs, 0x0002, "CanonFocalLength", &stdInterpreter},
 {0, 1, 0, 0, 0x0003, "CanonFlashInfo", &stdInterpreter},
 {0, 1, 0, canonShotInfoAttribs, 0x0004, "CanonShotInfo", &stdInterpreter},
 {0, 1, 0, canonPanoramaInfoAttribs, 0x0005, "CanonPanorama", &stdInterpreter},
 {0, 1, 0, 0, 0x0006, "CanonImageType", &stdInterpreter},
 {0, 1, 0, 0, 0x0007, "CanonFirmwareVersion", &stdInterpreter},
 {0, 1, 0, 0, 0x0008, "FileNumber", &stdInterpreter},
 {0, 1, 0, 0, 0x0009, "OwnerName", &stdInterpreter},
 {0, 1, 0, 0, 0x000a, "ColorInfoD30", &stdInterpreter},
 {0, 1, 0, 0, 0x000c, "SerialNumber", &stdInterpreter},
 {0, 1, 0, 0, 0x000d, "CanonCameraInfo", &stdInterpreter},
 {0, 1, 0, 0, 0x000e, "CanonFileLength", &stdInterpreter},
 {0, 1, 0, 0, 0x000f, "CustomFunctions", &stdInterpreter},
 {0, 1, 0, 0, 0x0010, "CanonModelID", &caModelIDInterpreter},
 {0, 1, 0, 0, 0x0012, "CanonAFInfo", &stdInterpreter},
 {0, 1, 0, 0, 0x0015, "SerialNumberFormat", &stdInterpreter},
 {0, 1, 0, 0, 0x001c, "DateStampMode", &stdInterpreter},
 {0, 1, 0, 0, 0x001d, "MyColors", &stdInterpreter},
 {0, 1, 0, 0, 0x001e, "FirmwareRevision", &stdInterpreter},
 {0, 3, 0, 0, 0x0024, "FaceDetect1", &stdInterpreter},
 {0, 3, 0, 0, 0x0025, "FaceDetect2", &stdInterpreter},
 {0, 1, 0, 0, 0x0026, "CanonAFInfo2", &stdInterpreter},
 {0, 1, 0, 0, 0x0083, "OriginalDecisionData", &stdInterpreter},
 {0, 1, 0, 0, 0x0090, "CustomFunctions1D", &stdInterpreter},
 {0, 1, 0, 0, 0x0091, "PersonalFunctions", &stdInterpreter},
 {0, 1, 0, 0, 0x0092, "PersonalFunctionValues", &stdInterpreter},
 {0, 1, 0, canonFileInfoAttribs, 0x0093, "CanonFileInfo", &stdInterpreter},
 {0, 1, 0, 0, 0x0094, "AFPointsInFocus1D", &stdInterpreter},
 {0, 1, 0, 0, 0x0095, "LensType", &stdInterpreter},
 {0, 1, 0, 0, 0x0096, "InternalSerialNumber", &caIntSerNumInterpreter},
 {0, 1, 0, 0, 0x0097, "DustRemovalData", &stdInterpreter},
 {0, 1, 0, canonCropInfoAttribs, 0x0098, "CropInfo", &stdInterpreter},
 {0, 1, 0, 0, 0x0099, "CustomFunctions2", &stdInterpreter},
 {0, 1, 0, canonAspectInfoAttribs, 0x009a, "AspectInfo", &stdInterpreter},
 {0, 1, 0, canonProcessingInfoAttribs, 0x00a0, "ProcessingInfo", &stdInterpreter},
 {0, 1, 0, 0, 0x00a1, "ToneCurveTable", &stdInterpreter},
 {0, 1, 0, 0, 0x00a2, "SharpnessTable", &stdInterpreter},
 {0, 1, 0, 0, 0x00a3, "SharpnessFreqTable", &stdInterpreter},
 {0, 1, 0, 0, 0x00a4, "WhiteBalanceTable", &stdInterpreter},
 {0, 1, 0, 0, 0x00a9, "ColorBalance", &stdInterpreter},
 {0, 1, 0, 0, 0x00aa, "MeasuredColor", &stdInterpreter},
 {0, 1, 0, 0, 0x00ae, "ColorTemperature", &stdInterpreter},
 {0, 3, 0, 0, 0x00b0, "CanonFlags", &stdInterpreter},
 {0, 1, 0, 0, 0x00b1, "ModifiedInfo", &stdInterpreter},
 {0, 1, 0, 0, 0x00b2, "ToneCurveMatching", &stdInterpreter},
 {0, 1, 0, 0, 0x00b3, "WhiteBalanceMatching", &stdInterpreter},
 {0, 1, 0, 0, 0x00b4, "ColorSpace", &stdInterpreter},
 {1, 1, 0, 0, 0x00b6, "PreviewImageInfo", &stdInterpreter},
 {0, 1, 0, 0, 0x00d0, "VRDOffset", &stdInterpreter},
 {0, 1, 0, 0, 0x00e0, "SensorInfo", &stdInterpreter},
 {0, 1, 0, 0, 0x4001, "ColorBalance", &stdInterpreter},
 {0, 1, 0, 0, 0x4002, "UnknownBlock1", &stdInterpreter},
 {0, 1, 0, 0, 0x4003, "ColorInfo", &stdInterpreter},
 {1, 1, 0, 0, 0x4005, "UnknownBlock2", &stdInterpreter},
 {1, 1, 0, 0, 0x4008, "BlackLevel", &stdInterpreter},
 {1, 1, 0, canonMicroAdjustAttrib, 0x4013, "AFMicroAdj", &stdInterpreter},
 {-1, 0, 0,  0, 0, "", NULL}};
};
#endif

