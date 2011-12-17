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
#include <string>
#include <sstream>

#include "rtexif.h"

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
            choices[65] = "Half-length Portrait";
            choices[221] = "P";
            choices[255] = "PICT";
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
            choices[3200] = "3200";
        }
};
PAISOInterpreter paISOInterpreter;

class PAFNumberInterpreter: public Interpreter {
public:
	PAFNumberInterpreter () {}
    virtual std::string toString (Tag* t) {
    	char buffer[32];
    	double v = t->toDouble()/10;
    	if( v < 0. || v > 1000. ) return "undef";
        sprintf (buffer, "%.1f", v );
        return buffer;
    }
};
PAFNumberInterpreter paFNumberInterpreter;

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

class PAPictureModeInterpreter2: public ChoiceInterpreter {
public:
	PAPictureModeInterpreter2(){
		choices[ 0] = "Program";
		choices[ 1] = "Hi-speed Program";
		choices[ 2] = "DOF Program";
		choices[ 3] = "MTF Program";
		choices[ 4] = "Standard";
		choices[ 5] = "Portrait";
		choices[ 6] = "Landscape";
		choices[ 7] = "Macro";
		choices[ 8] = "Sport ";
		choices[ 9] = "Night Scene Portrait ";
		choices[10] = "No Flash";
		choices[11] = "Night Scene";
		choices[12] = "Surf & Snow";
		choices[13] = "Text";
		choices[14] = "Sunset";
		choices[15] = "Kids";
		choices[16] = "Pet";
		choices[17] = "Candlelight";
		choices[18] = "Museum";
		choices[19] = "Food ";
		choices[20] = "Stage Lighting";
		choices[21] = "Night Snap";
		choices[256+4] = "Auto PICT";
		choices[256+5] = "Auto PICT (Portrait)";
		choices[256+6] = "Auto PICT (Landscape)";
		choices[256+7] = "Auto PICT (Macro)";
		choices[256+8] = "Auto PICT (Sport)";
		choices[256+8] = "Auto PICT (Sport)";
		choices[512+0] = "Program (HyP)";
		choices[512+1] = "Hi-speed Program (HyP)";
		choices[512+2] = "DOF Program (HyP)";
		choices[512+3] = "MTF Program (HyP)";
		choices[3*256] = "Green Mode";
		choices[4*256] = "Shutter Speed Priority";
		choices[5*256] = "Aperture Priority";
		choices[6*256] = "Program Tv Shift";
		choices[7*256] = "Program Av Shift";
		choices[8*256] = "Manual";
		choices[9*256] = "Bulb";
		choices[10*256] = "Aperture Priority, Off-Auto-Aperture";
		choices[11*256] = "Manual, Off-Auto-Aperture";
		choices[12*256] = "Bulb, Off-Auto-Aperture";
		choices[13*256] = "Shutter & Aperture Priority AE";
		choices[15*256] = "Sensitivity Priority AE";
		choices[16*256] = "Flash X-Sync Speed AE";
	}
    virtual std::string toString (Tag* t) {
    	int c = 256*t->toInt(0,BYTE) + t->toInt(1,BYTE);
        std::map<int,std::string>::iterator r = choices.find (c);
        if (r!=choices.end()){
        	std::ostringstream s;
        	s << r->second;
        	if( t->toInt(1,BYTE)==0 )
        		s << "\n1/2 EV steps";
        	else
        		s << "\n1/3 EV steps";
            return s.str();
        }else {
        	char buffer[1024];
            t->toString (buffer);
            return std::string (buffer);
        }
    }
};
PAPictureModeInterpreter2 paPictureModeInterpreter2;

class PADriveModeInterpreter : public ChoiceInterpreter{
	std::map<int,std::string> choices1;
	std::map<int,std::string> choices2;
	std::map<int,std::string> choices3;
public:
	PADriveModeInterpreter(){
		choices[0] = "Single-frame";
		choices[1] = "Continuous";
		choices[2] = "Continuous (Hi)";
		choices[3] = "Burst";
		choices[255] = "Video";
		choices1[0] = "No Timer";
		choices1[1] = "Self-timer (12 s)";
		choices1[2] = "Self-timer (2 s)";
		choices1[255] = "n/a";
		choices2[0] = "Shutter Button";
		choices2[1] = "Remote Control (3 s delay)";
		choices2[2] = "Remote Control";
		choices3[0] = "Single Exposure";
		choices3[1] = "Multiple Exposure";
		choices3[255] = "Video";
	}
    virtual std::string toString (Tag* t) {
        std::map<int,std::string>::iterator r  = choices.find (t->toInt(0,BYTE));
        std::map<int,std::string>::iterator r1 = choices1.find (t->toInt(1,BYTE));
        std::map<int,std::string>::iterator r2 = choices2.find (t->toInt(2,BYTE));
        std::map<int,std::string>::iterator r3 = choices3.find (t->toInt(3,BYTE));
        std::ostringstream s;
        s << ((r !=choices.end())? r->second : "");
        s << ((r1!=choices1.end())? r1->second : "")<<" ";
        s << ((r2!=choices2.end())? r2->second : "")<<" ";
        s << ((r3!=choices3.end())? r3->second : "")<<" ";
        return s.str();
    }
};
PADriveModeInterpreter paDriveModeInterpreter;

class PAColorSpaceInterpreter: public ChoiceInterpreter{
public:
	PAColorSpaceInterpreter(){
		choices[0] = "sRGB";
		choices[1] = "Adobe RGB";
	}
};
PAColorSpaceInterpreter paColorSpaceInterpreter;

class PALensTypeInterpreter : public IntLensInterpreter< int > {
    public:
        PALensTypeInterpreter () {
            choices.insert(p_t( 0,"M-42 or No Lens"));
            choices.insert(p_t(256*1, "K,M Lens"));
            choices.insert(p_t(256*2, "A Series Lens"));
            choices.insert(p_t(256*3+ 0, "Sigma Lens"));
            choices.insert(p_t(256*3+ 17, "smc PENTAX-FA SOFT 85mm f/2.8"));
            choices.insert(p_t(256*3+ 18, "smc PENTAX-F 1.7X AF ADAPTER"));
            choices.insert(p_t(256*3+ 19, "smc PENTAX-F 24-50mm f/4"));
            choices.insert(p_t(256*3+ 20, "smc PENTAX-F 35-80mm f/4-5.6"));
            choices.insert(p_t(256*3+ 21, "smc PENTAX-F 80-200mm f/4.7-5.6"));
            choices.insert(p_t(256*3+ 22, "smc PENTAX-F FISH-EYE 17-28mm f/3.5-4.5"));
            choices.insert(p_t(256*3+ 23, "smc PENTAX-F 100-300mm f/4.5-5.6"));
            choices.insert(p_t(256*3+ 23, "Sigma AF 28-300mm f/3.5-5.6 DL IF"));
            choices.insert(p_t(256*3+ 23, "Sigma AF 28-300mm f/3.5-6.3 DG IF Macro"));
            choices.insert(p_t(256*3+ 24, "smc PENTAX-F 35-135mm f/3.5-4.5"));
            choices.insert(p_t(256*3+ 25, "smc PENTAX-F 35-105mm f/4-5.6"));
            choices.insert(p_t(256*3+ 25, "Sigma AF 28-300mm f/3.5-5.6 DL IF"));
            choices.insert(p_t(256*3+ 25, "Sigma 55-200mm f/4-5.6 DC"));
            choices.insert(p_t(256*3+ 25, "Sigma AF 28-300mm f/3.5-5.6 DL IF"));
            choices.insert(p_t(256*3+ 25, "Sigma AF 28-300mm f/3.5-6.3 DG IF Macro"));
            choices.insert(p_t(256*3+ 25, "Tokina 80-200mm f/2.8 ATX-Pro"));
            choices.insert(p_t(256*3+ 26, "smc PENTAX-F* 250-600mm f/5.6 ED[IF]"));
            choices.insert(p_t(256*3+ 27, "smc PENTAX-F 28-80mm f/3.5-4.5"));
            choices.insert(p_t(256*3+ 27, "Tokina AT-X PRO AF 28-70mm f/2.6-2.8"));
            choices.insert(p_t(256*3+ 28, "smc PENTAX-F 35-70mm f/3.5-4.5"));
            choices.insert(p_t(256*3+ 28, "Tokina 19-35mm f/3.5-4.5 AF"));
            choices.insert(p_t(256*3+ 29, "PENTAX-F 28-80mm f/3.5-4.5"));
            choices.insert(p_t(256*3+ 29, "Sigma AF 18-125mm f/3.5-5.6 DC"));
            choices.insert(p_t(256*3+ 29, "Tokina AT-X PRO 28-70mm f/2.6-2.8"));
            choices.insert(p_t(256*3+ 30, "PENTAX-F 70-200mm f/4-5.6"));
            choices.insert(p_t(256*3+ 31, "smc PENTAX-F 70-210mm f/4-5.6"));
            choices.insert(p_t(256*3+ 31, "Tokina AF 730 75-300mm f/4.5-5.6"));
            choices.insert(p_t(256*3+ 31, "Takumar-F 70-210mm f/4-5.6"));
            choices.insert(p_t(256*3+ 32, "smc PENTAX-F 50mm f/1.4"));
            choices.insert(p_t(256*3+ 33, "smc PENTAX-F 50mm f/1.7"));
            choices.insert(p_t(256*3+ 34, "smc PENTAX-F 135mm f/2.8 [IF]"));
            choices.insert(p_t(256*3+ 35, "smc PENTAX-F 28mm f/2.8"));
            choices.insert(p_t(256*3+ 36, "Sigma 20mm f/1.8 EX DG ASPHERICAL RF"));
            choices.insert(p_t(256*3+ 38, "smc PENTAX-F* 300mm f/4.5 ED[IF]"));
            choices.insert(p_t(256*3+ 39, "smc PENTAX-F* 600mm f/4 ED[IF]"));
            choices.insert(p_t(256*3+ 40, "smc PENTAX-F MACRO 100mm f/2.8"));
            choices.insert(p_t(256*3+ 41, "smc PENTAX-F MACRO 50mm f/2.8"));
            choices.insert(p_t(256*3+ 41, "Sigma 50mm f/2.8 Macro"));
            choices.insert(p_t(256*3+ 44, "Sigma AF 10-20mm f/4-5.6 EX DC"));
            choices.insert(p_t(256*3+ 44, "Sigma 12-24mm f/4.5 EX DG"));
            choices.insert(p_t(256*3+ 44, "Sigma 17-70mm f/2.8-4.5 DC Macro"));
            choices.insert(p_t(256*3+ 44, "Sigma 18-50mm f/3.5-5.6 DC"));
            choices.insert(p_t(256*3+ 44, "Tamron 35-90mm f/4 AF"));
            choices.insert(p_t(256*3+ 46, "Sigma APO 70-200mm f/2.8 EX"));
            choices.insert(p_t(256*3+ 46, "Sigma EX APO 100-300mm f/4 IF"));
            choices.insert(p_t(256*3+ 50, "smc PENTAX-FA 28-70 f/4 AL"));
            choices.insert(p_t(256*3+ 51, "Sigma 28mm f/1.8 EX DG ASPHERICAL MACRO"));
            choices.insert(p_t(256*3+ 52, "smc PENTAX-FA 28-200mm f/3.8-5.6 AL[IF]"));
            choices.insert(p_t(256*3+ 52, "Tamron AF LD 28-200mm f/3.8-5.6 [IF] Aspherical (171D)"));
            choices.insert(p_t(256*3+ 53, "smc PENTAX-FA 28-80mm f/3.5-5.6 AL"));
            choices.insert(p_t(256*3+ 247,"smc PENTAX-DA FISH-EYE 10-17mm f/3.5-4.5 ED[IF]"));
            choices.insert(p_t(256*3+ 248,"smc PENTAX-DA 12-24mm f/4 ED AL[IF]"));
            choices.insert(p_t(256*3+ 250,"smc PENTAX-DA 50-200mm f/4-5.6 ED"));
            choices.insert(p_t(256*3+ 251,"smc PENTAX-DA 40mm f/2.8 Limited"));
            choices.insert(p_t(256*3+ 252,"smc PENTAX-DA 18-55mm f/3.5-5.6 AL"));
            choices.insert(p_t(256*3+ 253,"smc PENTAX-DA 14mm f/2.8 ED[IF]"));
            choices.insert(p_t(256*3+ 254,"smc PENTAX-DA 16-45mm f/4 ED AL"));
            choices.insert(p_t(256*3+ 255, "Sigma 18-200mm f/3.5-6.3 DC"));
            choices.insert(p_t(256*3+ 255, "Sigma DL-II 35-80mm f/4-5.6"));
            choices.insert(p_t(256*3+ 255, "Sigma DL Zoom 75-300mm f/4-5.6"));
            choices.insert(p_t(256*3+ 255, "Sigma DF EX Aspherical 28-70mm f/2.8"));
            choices.insert(p_t(256*3+ 255, "Sigma AF Tele 400mm f/5.6 Multi-coated"));
            choices.insert(p_t(256*3+ 255, "Sigma 24-60mm f/2.8 EX DG"));
            choices.insert(p_t(256*3+ 255, "Sigma 70-300mm f/4-5.6 Macro"));
            choices.insert(p_t(256*3+ 255, "Sigma 55-200mm f/4-5.6 DC"));
            choices.insert(p_t(256*3+ 255, "Sigma 18-50mm f/2.8 EX DC"));
            choices.insert(p_t(256*4+ 1,  "smc PENTAX-FA SOFT 28mm f/2.8"));
            choices.insert(p_t(256*4+ 2,  "smc PENTAX-FA 80-320mm f/4.5-5.6"));
            choices.insert(p_t(256*4+ 3,  "smc PENTAX-FA 43mm f/1.9 Limited"));
            choices.insert(p_t(256*4+ 6,  "smc PENTAX-FA 35-80mm f/4-5.6"));
            choices.insert(p_t(256*4+ 12, "smc PENTAX-FA 50mm f/1.4"));
            choices.insert(p_t(256*4+ 15, "smc PENTAX-FA 28-105mm f/4-5.6 [IF]"));
            choices.insert(p_t(256*4+ 16, "Tamron AF 80-210mm f/4-5.6 (178D)"));
            choices.insert(p_t(256*4+ 19, "Tamron SP AF 90mm f/2.8 (172E)"));
            choices.insert(p_t(256*4+ 20, "smc PENTAX-FA 28-80mm f/3.5-5.6"));
            choices.insert(p_t(256*4+ 21, "Cosina AF 100-300mm f/5.6-6.7"));
            choices.insert(p_t(256*4+ 22, "TOKINA 28-80mm f/3.5-5.6"));
            choices.insert(p_t(256*4+ 23, "smc PENTAX-FA 20-35mm f/4 AL"));
            choices.insert(p_t(256*4+ 24, "smc PENTAX-FA 77mm f/1.8 Limited"));
            choices.insert(p_t(256*4+ 25, "Tamron SP AF 14mm f/2.8"));
            choices.insert(p_t(256*4+ 26, "smc PENTAX-FA MACRO 100mm f/3.5"));
            choices.insert(p_t(256*4+ 26, "Cosina 100mm f/3.5 Macro"));
            choices.insert(p_t(256*4+ 27, "Tamron AF 28-300mm f/3.5-6.3 LD Aspherical[IF] MACRO (285D)"));
            choices.insert(p_t(256*4+ 28, "smc PENTAX-FA 35mm f/2 AL"));
            choices.insert(p_t(256*4+ 29, "Tamron AF 28-200mm f/3.8-5.6 LD Super II MACRO (371D)"));
            choices.insert(p_t(256*4+ 34, "smc PENTAX-FA 24-90mm f/3.5-4.5 AL[IF]"));
            choices.insert(p_t(256*4+ 35, "smc PENTAX-FA 100-300mm f/4.7-5.8"));
            choices.insert(p_t(256*4+ 36, "Tamron AF 70-300mm f/4-5.6 LD MACRO"));
            choices.insert(p_t(256*4+ 37, "Tamron SP AF 24-135mm f/3.5-5.6 AD AL (190D)"));
            choices.insert(p_t(256*4+ 38, "smc PENTAX-FA 28-105mm f/3.2-4.5 AL[IF]"));
            choices.insert(p_t(256*4+ 39, "smc PENTAX-FA 31mm f/1.8 AL Limited"));
            choices.insert(p_t(256*4+ 41, "Tamron AF 28-200mm Super Zoom f/3.8-5.6 Aspherical XR [IF] MACRO (A03)"));
            choices.insert(p_t(256*4+ 43, "smc PENTAX-FA 28-90mm f/3.5-5.6"));
            choices.insert(p_t(256*4+ 44, "smc PENTAX-FA J 75-300mm f/4.5-5.8 AL"));
            choices.insert(p_t(256*4+ 45, "Tamron 28-300mm f/3.5-6.3 Ultra zoom XR"));
            choices.insert(p_t(256*4+ 45, "Tamron AF 28-300mm f/3.5-6.3 XR Di LD Aspherical [IF] Macro"));
            choices.insert(p_t(256*4+ 46, "smc PENTAX-FA J 28-80mm f/3.5-5.6 AL"));
            choices.insert(p_t(256*4+ 47, "smc PENTAX-FA J 18-35mm f/4-5.6 AL"));
            choices.insert(p_t(256*4+ 49, "Tamron SP AF 28-75mm f/2.8 XR Di (A09)"));
            choices.insert(p_t(256*4+ 51, "smc PENTAX-D FA 50mm f/2.8 MACRO"));
            choices.insert(p_t(256*4+ 52, "smc PENTAX-D FA 100mm f/2.8 MACRO"));
            choices.insert(p_t(256*4+ 75, "Tamron SP AF 70-200 f/2.8 Di LD [IF] Macro (A001)"));
            choices.insert(p_t(256*4+ 229, "smc PENTAX-DA 18-55mm f/3.5-5.6 AL II"));
            choices.insert(p_t(256*4+ 230, "Tamron SP AF 17-50mm f/2.8 XR Di II"));
            choices.insert(p_t(256*4+ 231, "smc PENTAX-DA 18-250mm f/3.5-6.3 ED AL [IF]"));
            choices.insert(p_t(256*4+ 237, "Samsung/Schneider D-XENOGON 10-17mm f/3.5-4.5"));
            choices.insert(p_t(256*4+ 239, "Samsung D-XENON 12-24mm f/4 ED AL [IF]"));
            choices.insert(p_t(256*4+ 243, "smc PENTAX-DA 70mm f/2.4 Limited"));
            choices.insert(p_t(256*4+ 244, "smc PENTAX-DA 21mm f/3.2 AL Limited"));
            choices.insert(p_t(256*4+ 245, "Schneider D-XENON 50-200mm"));
            choices.insert(p_t(256*4+ 246, "Schneider D-XENON 18-55mm"));
            choices.insert(p_t(256*4+ 247, "smc PENTAX-DA 10-17mm f/3.5-4.5 ED [IF] Fisheye zoom"));
            choices.insert(p_t(256*4+ 248, "smc PENTAX-DA 12-24mm f/4 ED AL [IF]"));
            choices.insert(p_t(256*4+ 249, "Tamron 18-200mm f/3.5-6.3 XR DiII (A14)"));
            choices.insert(p_t(256*4+ 250, "smc PENTAX-DA 50-200mm f/4-5.6 ED"));
            choices.insert(p_t(256*4+ 251, "smc PENTAX-DA 40mm f/2.8 Limited"));
            choices.insert(p_t(256*4+ 252, "smc PENTAX-DA 18-55mm f/3.5-5.6 AL"));
            choices.insert(p_t(256*4+ 253, "smc PENTAX-DA 14mm f/2.8 ED[IF]"));
            choices.insert(p_t(256*4+ 254, "smc PENTAX-DA 16-45mm f/4 ED AL"));
            choices.insert(p_t(256*5+ 1, "smc PENTAX-FA* 24mm f/2 AL[IF]"));
            choices.insert(p_t(256*5+ 2, "smc PENTAX-FA 28mm f/2.8 AL"));
            choices.insert(p_t(256*5+ 3, "smc PENTAX-FA 50mm f/1.7"));
            choices.insert(p_t(256*5+ 4, "smc PENTAX-FA 50mm f/1.4"));
            choices.insert(p_t(256*5+ 5, "smc PENTAX-FA* 600mm f/4 ED[IF]"));
            choices.insert(p_t(256*5+ 6, "smc PENTAX-FA* 300mm f/4.5 ED[IF]"));
            choices.insert(p_t(256*5+ 7, "smc PENTAX-FA 135mm f/2.8 [IF]"));
            choices.insert(p_t(256*5+ 8, "smc PENTAX-FA MACRO 50mm f/2.8"));
            choices.insert(p_t(256*5+ 9, "smc PENTAX-FA MACRO 100mm f/2.8"));
            choices.insert(p_t(256*5+ 10, "smc PENTAX-FA* 85mm f/1.4 [IF]"));
            choices.insert(p_t(256*5+ 11, "smc PENTAX-FA* 200mm f/2.8 ED[IF]"));
            choices.insert(p_t(256*5+ 12, "smc PENTAX-FA 28-80mm f/3.5-4.7"));
            choices.insert(p_t(256*5+ 13, "smc PENTAX-FA 70-200mm f/4-5.6"));
            choices.insert(p_t(256*5+ 14, "smc PENTAX-FA* 250-600mm f/5.6 ED[IF]"));
            choices.insert(p_t(256*5+ 15, "smc PENTAX-FA 28-105mm f/4-5.6"));
            choices.insert(p_t(256*5+ 16, "smc PENTAX-FA 100-300mm f/4.5-5.6"));
            choices.insert(p_t(256*5+ 98, "smc PENTAX-FA 100-300mm f/4.5-5.6"));
            choices.insert(p_t(256*6+ 1, "smc PENTAX-FA* 85mm f/1.4 [IF]"));
            choices.insert(p_t(256*6+ 2, "smc PENTAX-FA* 200mm f/2.8 ED[IF]"));
            choices.insert(p_t(256*6+ 3, "smc PENTAX-FA* 300mm f/2.8 ED[IF]"));
            choices.insert(p_t(256*6+ 4, "smc PENTAX-FA* 28-70mm f/2.8 AL"));
            choices.insert(p_t(256*6+ 5, "smc PENTAX-FA* 80-200mm f/2.8 ED[IF]"));
            choices.insert(p_t(256*6+ 6, "smc PENTAX-FA* 28-70mm f/2.8 AL"));
            choices.insert(p_t(256*6+ 7, "smc PENTAX-FA* 80-200mm f/2.8 ED[IF]"));
            choices.insert(p_t(256*6+ 8, "smc PENTAX-FA 28-70mm f/4 AL"));
            choices.insert(p_t(256*6+ 9, "smc PENTAX-FA 20mm f/2.8"));
            choices.insert(p_t(256*6+ 10, "smc PENTAX-FA* 400mm f/5.6 ED[IF]"));
            choices.insert(p_t(256*6+ 13, "smc PENTAX-FA* 400mm f/5.6 ED[IF]"));
            choices.insert(p_t(256*6+ 14, "smc PENTAX-FA* MACRO 200mm f/4 ED[IF]"));
            choices.insert(p_t(256*7+ 0, "smc PENTAX-DA 21mm f/3.2 AL Limited"));
            choices.insert(p_t(256*7+ 58, "smc PENTAX-D FA MACRO 100mm f/2.8 WR"));
            choices.insert(p_t(256*7+ 75, "Tamron SP AF 70-200mm f/2.8 Di LD [IF] Macro (A001)"));
            choices.insert(p_t(256*7+ 214, "smc PENTAX-DA 35mm f/2.4 AL"));
            choices.insert(p_t(256*7+ 216, "smc PENTAX-DA L 55-300mm f/4-5.8 ED"));
            choices.insert(p_t(256*7+ 217, "smc PENTAX-DA 50-200mm f/4-5.6 ED WR"));
            choices.insert(p_t(256*7+ 218, "smc PENTAX-DA 18-55mm f/3.5-5.6 AL WR"));
            choices.insert(p_t(256*7+ 220, "Tamron SP AF 10-24mm f/3.5-4.5 Di II LD Aspherical [IF]"));
            choices.insert(p_t(256*7+ 222, "smc PENTAX-DA 18-55mm f/3.5-5.6 AL II"));
            choices.insert(p_t(256*7+ 223, "Samsung D-XENON 18-55mm f/3.5-5.6 II"));
            choices.insert(p_t(256*7+ 224, "smc PENTAX-DA 15mm f/4 ED AL Limited"));
            choices.insert(p_t(256*7+ 225, "Samsung D-XENON 18-250mm f/3.5-6.3"));
            choices.insert(p_t(256*7+ 226, "smc PENTAX-DA* 55mm f/1.4 SDM (SDM unused)"));
            choices.insert(p_t(256*7+ 227, "smc PENTAX DA* 60-250mm f/4 [IF] SDM (SDM unused)"));
            choices.insert(p_t(256*7+ 229, "smc PENTAX-DA 18-55mm f/3.5-5.6 AL II"));
            choices.insert(p_t(256*7+ 230, "Tamron AF 17-50mm f/2.8 XR Di-II LD (Model A16)"));
            choices.insert(p_t(256*7+ 231, "smc PENTAX-DA 18-250mm f/3.5-6.3 ED AL [IF]"));
            choices.insert(p_t(256*7+ 233, "smc PENTAX-DA 35mm f/2.8 Macro Limited"));
            choices.insert(p_t(256*7+ 234, "smc PENTAX-DA* 300mm f/4 ED [IF] SDM (SDM unused)"));
            choices.insert(p_t(256*7+ 235, "smc PENTAX-DA* 200mm f/2.8 ED [IF] SDM (SDM unused)"));
            choices.insert(p_t(256*7+ 236, "smc PENTAX-DA 55-300mm f/4-5.8 ED"));
            choices.insert(p_t(256*7+ 238, "Tamron AF 18-250mm f/3.5-6.3 Di II LD Aspherical [IF] MACRO"));
            choices.insert(p_t(256*7+ 241, "smc PENTAX-DA* 50-135mm f/2.8 ED [IF] SDM (SDM unused)"));
            choices.insert(p_t(256*7+ 242, "smc PENTAX-DA* 16-50mm f/2.8 ED AL [IF] SDM (SDM unused)"));
            choices.insert(p_t(256*7+ 243, "smc PENTAX-DA 70mm f/2.4 Limited"));
            choices.insert(p_t(256*7+ 244, "smc PENTAX-DA 21mm f/3.2 AL Limited"));
            choices.insert(p_t(256*8+  14, "Sigma 17-70mm f/2.8-4.0 DC Macro OS HSM"));
            choices.insert(p_t(256*8+ 215, "smc PENTAX-DA 18-135mm f/3.5-5.6 ED AL [IF] DC WR"));
            choices.insert(p_t(256*8+ 226, "smc PENTAX-DA* 55mm f/1.4 SDM"));
            choices.insert(p_t(256*8+ 227, "smc PENTAX DA* 60-250mm f/4 [IF] SDM"));
            choices.insert(p_t(256*8+ 232, "smc PENTAX-DA 17-70mm f/4 AL [IF] SDM"));
            choices.insert(p_t(256*8+ 234, "smc PENTAX-DA* 300mm f/4 ED [IF] SDM"));
            choices.insert(p_t(256*8+ 235, "smc PENTAX-DA* 200mm f/2.8 ED [IF] SDM"));
            choices.insert(p_t(256*8+ 241, "smc PENTAX-DA* 50-135mm f/2.8 ED [IF] SDM"));
            choices.insert(p_t(256*8+ 242, "smc PENTAX-DA* 16-50mm f/2.8 ED AL [IF] SDM"));
            choices.insert(p_t(256*8+ 255, "Sigma 70-200mm f/2.8 EX DG Macro HSM II"));
            choices.insert(p_t(256*8+ 255, "Sigma APO 150-500mm f/5-6.3 DG OS HSM"));
            choices.insert(p_t(256*11+  4, "smc PENTAX-FA 645 45-85mm f/4.5"));
            choices.insert(p_t(256*11+  8, "smc PENTAX-FA 645 80-160mm f/4.5"));
            choices.insert(p_t(256*11+ 11, "smc PENTAX-FA 645 35mm f/3.5 AL [IF]"));
            choices.insert(p_t(256*11+ 17, "smc PENTAX-FA 645 150-300mm f/5.6 ED [IF]"));
            choices.insert(p_t(256*13+ 18, "smc PENTAX-D FA 645 55mm f/2.8 AL [IF] SDM AW"));
        }
        virtual std::string toString (Tag* t) {
       	   double maxApertureAtFocal = 0;
       	   double focalLength = 0;
           int lensID = 256*t->toInt(0,BYTE) + t->toInt(1,BYTE);
           TagDirectory *root=t->getParent()->getRoot();
           if (root){
       	      Tag *t1;
              t1 = root->findTag("FocalLength");
        	  if( t1)
      	         focalLength = t1->toDouble(); // Focal Length
        	  t1 = root->findTag("MaxAperture");
        	  if( t1){
        		  int a=t1->toInt(0,BYTE)&0x7F;
   	              maxApertureAtFocal = pow(2.0, (a-1)/32.0) ; // MaxApertureValue at focal Length
        	  }
           }
       	   return guess( lensID, focalLength, maxApertureAtFocal);
        }
};
PALensTypeInterpreter paLensTypeInterpreter;

class PASRResultInterpreter: public Interpreter {
public:
	PASRResultInterpreter(){ }
    virtual std::string toString (Tag* t) {
        std::ostringstream str;
        int b = t->toInt(0,BYTE);
        if (!b)
            str << "Not stabilized";
        else if (b & 1)
            str << "Stabilized";
        else if (b & 64)
        	str << "Not Ready";
        return str.str();
	}
};
PASRResultInterpreter paSRResultInterpreter;

class PAHighISONoiseInterpreter: public ChoiceInterpreter {
public:
	PAHighISONoiseInterpreter(){
		choices[0] = "Off";
		choices[1] = "Weakest";
		choices[2] = "Weak";
		choices[3] = "Strong";
	}
};
PAHighISONoiseInterpreter paHighISONoiseInterpreter;

class PAPowerSourceInterpreter: public ChoiceInterpreter {
public:
	PAPowerSourceInterpreter(){
		choices[2] = "Body Battery";
		choices[3] = "Grip Battery ";
		choices[4] = "External Power Supply";
	}
};
PAPowerSourceInterpreter paPowerSourceInterpreter;

class PAMaxApertureInterpreter: public Interpreter {
	public:
	   PAMaxApertureInterpreter(){}
       virtual std::string toString (Tag* t){
    	   int a = t->toInt(0,BYTE);
    	   a &= 0x7F;
    	   if(a>1){
    		  char buffer[32];
    		  double v = pow(2.0, (a-1)/32.0);
    		  if( v < 0. || v > 1000. ) return "undef";
              sprintf (buffer, "%.1f", v );
              return buffer;
    	   }else
    		  return "n/a";
       }
};
PAMaxApertureInterpreter paMaxApertureInterpreter;

class PANominalMinMaxApertureInterpreter: public Interpreter {
public:
	PANominalMinMaxApertureInterpreter(){}
    virtual std::string toString (Tag* t){
       char buffer[1024];
 	   int a = t->toInt(0,BYTE);
 	   int mina = a & 0x0F;
 	   int maxa = (a & 0xF0)>>4;
       sprintf (buffer, "%.1f - %.0f", pow(2.0, maxa/4.0), pow(2.0, (mina+10)/4.0));
       return buffer;

    }
};
PANominalMinMaxApertureInterpreter paNominalMinMaxApertureInterpreter;

class PAFlashStatusInterpreter: public ChoiceInterpreter {
public:
	PAFlashStatusInterpreter(){
		choices[0x0] = "Off";
		choices[0x2] = "External, Did not fire";
		choices[0x6] = "External, Fired";
		choices[0x9] = "Internal, Did not fire";
		choices[0xd] = "Internal, Fired";
	}
};
PAFlashStatusInterpreter paFlashStatusInterpreter;

class PAInternalFlashModeInterpreter: public ChoiceInterpreter {
public:
	PAInternalFlashModeInterpreter(){
		choices[0x0] = "n/a - Off-Auto-Aperture";
		choices[0x86] = "On, Wireless (Control)";
		choices[0x95] = "On, Wireless (Master)";
		choices[0xc0] = "On";
		choices[0xc1] = "On, Red-eye reduction";
		choices[0xc2] = "On, Auto";
		choices[0xc3] = "On, Auto, Red-eye reduction";
		choices[0xc8] = "On, Slow-sync";
		choices[0xc9] = "On, Slow-sync, Red-eye reduction";
		choices[0xca] = "On, Trailing-curtain Sync";
		choices[0xf0] = "Off, Normal";
		choices[0xf1] = "Off, Red-eye reduction";
		choices[0xf2] = "Off, Auto";
		choices[0xf3] = "Off, Auto, Red-eye reduction";
		choices[0xf4] = "Off, (Unknown 0xf4)";
		choices[0xf5] = "Off, Wireless (Master)";
		choices[0xf6] = "Off, Wireless (Control)";
		choices[0xf8] = "Off, Slow-sync";
		choices[0xf9] = "Off, Slow-sync, Red-eye reduction";
		choices[0xfa] = "Off, Trailing-curtain Sync";
	}
};
PAInternalFlashModeInterpreter paInternalFlashModeInterpreter;

class PAExternalFlashModeInterpreter: public ChoiceInterpreter {
	public:
	PAExternalFlashModeInterpreter(){
		choices[0x0 ]= "n/a - Off-Auto-Aperture";
		choices[0x3f] = "Off";
		choices[0x40] = "On, Auto";
		choices[0xbf] = "On, Flash Problem";
		choices[0xc0] = "On, Manual";
		choices[0xc4] = "On, P-TTL Auto";
		choices[0xc5] = "On, Contrast-control Sync";
		choices[0xc6] = "On, High-speed Sync";
		choices[0xcc] = "On, Wireless";
		choices[0xcd] = "On, Wireless, High-speed Sync";
	}
};
PAExternalFlashModeInterpreter paExternalFlashModeInterpreter;

class PAExternalFlashExposureCompInterpreter: public ChoiceInterpreter {
	public:
		PAExternalFlashExposureCompInterpreter(){
			choices[0] = "n/a";
			choices[144] = "n/a (Manual Mode)";
			choices[164] = "-3.0";
			choices[167] = "-2.5";
			choices[168] = "-2.0";
			choices[171] = "-1.5";
			choices[172] = "-1.0";
			choices[175] = "-0.5";
			choices[176] = "0";
			choices[179] = "+0.5";
			choices[180] = "+1.0";
		}
};
PAExternalFlashExposureCompInterpreter paExternalFlashExposureCompInterpreter;

class PAExternalFlashBounceInterpreter: public ChoiceInterpreter {
	public:
		PAExternalFlashBounceInterpreter(){
			choices[0] = "n/a";
			choices[16] = "Direct";
			choices[48] = "Bonce";
		}
};
PAExternalFlashBounceInterpreter paExternalFlashBounceInterpreter;

class PAExternalFlashGNInterpreter: public Interpreter {
	public:
	PAExternalFlashGNInterpreter(){}
	virtual std::string toString (Tag* t) {
		   char buffer[1024];
	       int b = t->toInt(0,BYTE) & 0x1F;
	       sprintf (buffer, "%.0f", pow(2.,b/16.+4) );
	       return buffer;
	}
};
PAExternalFlashGNInterpreter paExternalFlashGNInterpreter;

class PAEVStepsInterpreter:public Interpreter {
public:
	PAEVStepsInterpreter(){}
	virtual std::string toString (Tag* t) {
		std::ostringstream str;
		if( t->toInt(0,BYTE) & 0x20 )
		   str << "1/3 EV steps";
		else
			str << "1/2 EV steps";
		return str.str();
	}
};
PAEVStepsInterpreter paEVStepsInterpreter;

class PAEDialinInterpreter:public Interpreter {
public:
	PAEDialinInterpreter(){}
	virtual std::string toString (Tag* t) {
		std::ostringstream str;
		if(  t->toInt(0,BYTE) & 0x40 )
		   str << "P Shift";
		else
			str << "Tv or Av";
	    return str.str();
	}
};
PAEDialinInterpreter paEDialinInterpreter;

class PAApertureRingUseInterpreter: public Interpreter {
public:
	PAApertureRingUseInterpreter(){}
	virtual std::string toString (Tag* t) {
		std::ostringstream str;
		if(  t->toInt(0,BYTE) & 0x80 )
			str << "Permitted";
		else
			str << "Prohibited";
	    return str.str();
	}
};
PAApertureRingUseInterpreter paApertureRingUseInterpreter;

class PAFlashOptionInterpreter: public ChoiceInterpreter {
public:
	PAFlashOptionInterpreter(){
		choices[0x0] = "Normal";
		choices[0x1] = "Red-eye reduction";
		choices[0x2] = "Auto";
		choices[0x3] = "Auto, Red-eye reduction";
		choices[0x5] = "Wireless (Master)";
		choices[0x6] = "Wireless (Control)";
		choices[0x8] = "Slow-sync";
		choices[0x9] = "Slow-sync, Red-eye reduction";
		choices[0xa] = "Trailing-curtain Sync";
	}
	virtual std::string toString (Tag* t) {
		std::map<int,std::string>::iterator r = choices.find (t->toInt(0,BYTE) >>4);
		if (r!=choices.end())
			return r->second;
		else {
			char buffer[1024];
			t->toString (buffer);
			return std::string (buffer);
		}
	}
};
PAFlashOptionInterpreter paFlashOptionInterpreter;

class PAMeteringMode2Interpreter: public Interpreter {
public:
	PAMeteringMode2Interpreter(){}
	virtual std::string toString (Tag* t) {
		std::ostringstream str;
		int v = (t->toInt(0,BYTE) & 0xF);
		if(!v)
			str << "Multi-segment";
		else if(v&1)
			str << "Center-weighted average";
		else if(v&2)
			str << "Spot";
		return str.str();
	}
};
PAMeteringMode2Interpreter paMeteringMode2Interpreter;

class PAExposureBracketStepSizeInterpreter: public ChoiceInterpreter {
public:
	PAExposureBracketStepSizeInterpreter(){
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

class PAPictureMode2Interpreter: public ChoiceInterpreter {
public:
	PAPictureMode2Interpreter(){
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

class PAProgramLineInterpreter: public Interpreter {
public:
	PAProgramLineInterpreter(){}
	virtual std::string toString (Tag* t) {
		std::ostringstream str;
		int c = t->toInt(0,BYTE);
		switch(c & 0xf){
		case 0: str << "Manual";break;
		case 1: str << "AF-S";break;
		case 2: str << "AF-C";break;
		case 3: str << "AF-A";break;
		}
		if( (c & 0xF0) == 0) str << ", Point Selection Auto";
		else if( c & 0x20 ) str << ", Fixed Center Point Selected";
		else if( c & 0x10 ) str << ", Point Selected";
		return str.str();
	}
};
PAProgramLineInterpreter paProgramLineInterpreter;

class PAAFModeInterpreter: public Interpreter {
public:
	PAAFModeInterpreter(){}
	virtual std::string toString (Tag* t) {
		switch(t->toInt(0,BYTE) & 0x3){
		case 0: return "Normal";
		case 1: return "Hi Speed";
		case 2: return "Depth";
		case 3: return "MTF";
		}
		return"Normal";
	}

};
PAAFModeInterpreter paAFModeInterpreter;

class PAAFPointSelectedInterpreter: public Interpreter {
public:
	PAAFPointSelectedInterpreter(){}
	virtual std::string toString (Tag* t) {
		const char *ps[]={"Upper-left","Top","Upper-right","Left","Mid-left","Center","Mid-right","Right","Lower-left","Bottom","Lower-right"};
		int c = t->toInt(0,SHORT);
		if( !c )
			return "Auto";
		else{
			for( int iBit=0; iBit<11; iBit++)
				if( c & (1<<iBit) )
				   return ps[iBit];
			return "n/a";
		}
	}
};
PAAFPointSelectedInterpreter paAFPointSelectedInterpreter;

class PADriveMode2Interpreter: public Interpreter {
public:
	PADriveMode2Interpreter(){}
	virtual std::string toString (Tag* t) {
		int c = t->toInt(0,BYTE);
		if( !c )
			return "Single-frame";
		else if( c & 0x01)
			return "Continuous";
		else if( c & 0x04)
			return "Self-timer (12 s)";
		else if( c & 0x08)
			return "Self-timer (2 s)";
		else if( c & 0x10 )
			return "Remote Control (3 s delay)";
		else if( c & 0x20)
			return "Remote Control";
		else if( c & 0x40)
			return "Exposure Bracket";
		else if( c & 0x80)
			return "Multiple Exposure";
		else
			return "Unknown";
	}
};
PADriveMode2Interpreter paDriveMode2Interpreter;

const TagAttrib pentaxAttribs[] = {
 {0, 1, 0, 0, 0x0000, "PentaxVersion", &stdInterpreter},
 {0, 1, 0, 0, 0x0001, "PentaxModelType", &stdInterpreter},
 {0, 2, 0, 0, 0x0002, "PreviewImageSize", &stdInterpreter},
 {0, 2, 0, 0, 0x0003, "PreviewImageLength", &stdInterpreter},
 {0, 2, 0, 0, 0x0004, "PreviewImageStart", &stdInterpreter},
 {0, 1, 0, 0, 0x0005, "PentaxModelID", &stdInterpreter},
 {0, 1, 0, 0, 0x0006, "Date", &stdInterpreter},
 {0, 1, 0, 0, 0x0007, "Time", &stdInterpreter},
 {0, 1, 0, 0, 0x0008, "Quality", &paQualityInterpreter},
 {0, 1, 0, 0, 0x0009, "PentaxImageSize", &stdInterpreter},
 {0, 1, 0, 0, 0x000b, "PictureMode", &paPictureModeInterpreter},
 {0, 1, 0, 0, 0x000c, "FlashMode", &paFlashModeInterpreter},
 {0, 1, 0, 0, 0x000d, "FocusMode", &paFocusModeInterpreter},
 {0, 1, 0, 0, 0x000e, "AFPointSelected", &paAFPointInterpreter},
 {0, 1, 0, 0, 0x000f, "AFPointsInFocus", &paAFFocusInterpreter},
 {0, 1, 0, 0, 0x0010, "FocusPosition", &stdInterpreter},
 {0, 1, 0, 0, 0x0012, "ExposureTime", &stdInterpreter},
 {0, 1, 0, 0, 0x0013, "FNumber", &paFNumberInterpreter},
 {0, 1, 0, 0, 0x0014, "ISO", &paISOInterpreter},
 {0, 1, 0, 0, 0x0015, "LightReading", &stdInterpreter},
 {0, 1, 0, 0, 0x0016, "ExposureCompensation", &stdInterpreter},
 {0, 1, 0, 0, 0x0017, "MeteringMode", &paMeteringModeInterpreter},
 {0, 1, 0, 0, 0x0018, "AutoBracketing", &stdInterpreter},
 {0, 1, 0, 0, 0x0019, "WhiteBalance", &paWhiteBalanceInterpreter},
 {0, 1, 0, 0, 0x001a, "WhiteBalanceMode", &paWhiteBalanceModeInterpreter},
 {0, 1, 0, 0, 0x001b, "BlueBalance", &stdInterpreter},
 {0, 1, 0, 0, 0x001c, "RedBalance", &stdInterpreter},
 {0, 1, 0, 0, 0x001d, "FocalLength", &stdInterpreter},
 {0, 1, 0, 0, 0x001e, "DigitalZoom", &stdInterpreter},
 {0, 1, 0, 0, 0x001f, "Saturation", &paSaturationInterpreter},
 {0, 1, 0, 0, 0x0020, "Contrast", &paContrastInterpreter},
 {0, 1, 0, 0, 0x0021, "Sharpness", &paSharpnessInterpreter},
 {0, 1, 0, 0, 0x0022, "WorldTimeLocation", &stdInterpreter},
 {0, 1, 0, 0, 0x0023, "HometownCity", &stdInterpreter},
 {0, 3, 0, 0, 0x0024, "DestinationCity", &stdInterpreter},
 {0, 3, 0, 0, 0x0025, "HometownDST", &stdInterpreter},
 {0, 1, 0, 0, 0x0026, "DestinationDST", &stdInterpreter},
 {0, 1, 0, 0, 0x0027, "DSPFirmwareVersion", &stdInterpreter},
 {0, 1, 0, 0, 0x0028, "CPUFirmwareVersion", &stdInterpreter},
 {0, 1, 0, 0, 0x0029, "FrameNumber", &stdInterpreter},
 {0, 1, 0, 0, 0x002d, "EffectiveLV", &stdInterpreter},
 {0, 1, 0, 0, 0x0032, "ImageProcessing", &stdInterpreter},
 {0, 1, 0, 0, 0x0033, "PictureMode", &paPictureModeInterpreter2},
 {0, 1, 0, 0, 0x0034, "DriveMode", &paDriveModeInterpreter},
 {0, 1, 0, 0, 0x0037, "ColorSpace", &paColorSpaceInterpreter},
 {0, 1, 0, 0, 0x0038, "ImageAreaOffset", &stdInterpreter},
 {0, 1, 0, 0, 0x0039, "RawImageSize", &stdInterpreter},
 {0, 1, 0, 0, 0x003c, "AFPointsInFocus", &stdInterpreter},
 {0, 1, 0, 0, 0x003e, "PreviewImageBorders", &stdInterpreter},
 {0, 1, 0, 0, 0x003f, "LensType", &paLensTypeInterpreter},
 {0, 1, 0, 0, 0x0040, "SensitivityAdjust", &stdInterpreter},
 {0, 1, 0, 0, 0x0041, "ImageProcessingCount", &stdInterpreter},
 {0, 1, 0, 0, 0x0047, "CameraTemperature", &stdInterpreter},
 {0, 1, 0, 0, 0x0048, "AELock", &paOnOffInterpreter},
 {0, 1, 0, 0, 0x0049, "NoiseReduction", &paOnOffInterpreter},
 {0, 1, 0, 0, 0x004d, "FlashExposureComp", &stdInterpreter},
 {0, 1, 0, 0, 0x004f, "ImageTone", &stdInterpreter},
 {0, 1, 0, 0, 0x0050, "ColorTemperature", &stdInterpreter},
 {0, 1, 0, pentaxSRInfoAttribs, 0x005c, "ShakeReductionInfo", &stdInterpreter},
 {0, 1, 0, 0, 0x005d, "ShutterCount", &stdInterpreter},
 {0, 1, 0, 0, 0x0069, "DynamicRangeExpansion", &paOnOffInterpreter},
 {0, 1, 0, 0, 0x0071, "HighISONoiseReduction", &paHighISONoiseInterpreter},
 {0, 1, 0, 0, 0x0072, "AFAdjustment", &stdInterpreter},
 {0, 1, 0, 0, 0x0200, "BlackPoint", &stdInterpreter},
 {0, 1, 0, 0, 0x0201, "WhitePoint", &stdInterpreter},
 {0, 1, 0, 0, 0x0203, "ColorMatrixA", &stdInterpreter},
 {0, 1, 0, 0, 0x0204, "ColorMatrixB", &stdInterpreter},
 {0, 1, 0, pentaxCameraSettingsAttribs, 0x0205, "CameraSettings", &stdInterpreter},
 {0, 1, 0, pentaxAEInfoAttribs, 0x0206, "AEInfo", &stdInterpreter},
 {0, 1, 0, pentaxLensDataAttribs, 0x0207, "LensInfo", &stdInterpreter},
 {0, 1, 0, pentaxFlashInfoAttribs, 0x0208, "FlashInfo", &stdInterpreter},
 {0, 1, 0, 0, 0x0209, "AEMeteringSegments", &stdInterpreter},
 {0, 1, 0, 0, 0x020a, "FlashADump", &stdInterpreter},
 {0, 1, 0, 0, 0x020b, "FlashBDump", &stdInterpreter},
 {0, 1, 0, 0, 0x020d, "WB_RGGBLevelsDaylight", &stdInterpreter},
 {0, 1, 0, 0, 0x020e, "WB_RGGBLevelsShade", &stdInterpreter},
 {0, 1, 0, 0, 0x020f, "WB_RGGBLevelsCloudy", &stdInterpreter},
 {0, 1, 0, 0, 0x0210, "WB_RGGBLevelsTungsten", &stdInterpreter},
 {0, 1, 0, 0, 0x0211, "WB_RGGBLevelsFluorescentD", &stdInterpreter},
 {0, 1, 0, 0, 0x0212, "WB_RGGBLevelsFluorescentN", &stdInterpreter},
 {0, 1, 0, 0, 0x0213, "WB_RGGBLevelsFluorescentW", &stdInterpreter},
 {0, 1, 0, 0, 0x0214, "WB_RGGBLevelsFlash", &stdInterpreter},
 {0, 1, 0, pentaxCameraInfoAttribs, 0x0215, "CameraInfo", &stdInterpreter},
 {0, 1, 0, pentaxBatteryInfoAttribs, 0x0216, "BatteryInfo", &stdInterpreter},
 {0, 1, 0, 0, 0x021f, "AFInfo", &stdInterpreter},
 {0, 1, 0, 0, 0x0222, "ColorInfo", &stdInterpreter},
 {0, 1, 0, 0, 0x03fe, "DataDump", &stdInterpreter},
 {0, 1, 0, 0, 0x03ff, "UnknownInfo", &stdInterpreter},
 {0, 1, 0, 0, 0x0402, "ToneCurve", &stdInterpreter},
 {0, 1, 0, 0, 0x0403, "ToneCurves", &stdInterpreter},
 {0, 1, 0, 0, 0x0e00, "PrintIM", &stdInterpreter},
 {-1, 0, 0,  0, 0, "", NULL}};

const TagAttrib pentaxSRInfoAttribs[] = {
 {0, 1, 0, 0,  0, "SRResult", &paSRResultInterpreter},
 {0, 1, 0, 0,  1, "ShakeReduction", &paOnOffInterpreter},
 {0, 1, 0, 0,  2, "SRHalfPressTime", &stdInterpreter},
 {0, 1, 0, 0,  3, "SRFocalLength", &stdInterpreter},
 {-1, 0, 0,  0, 0, "", NULL}};

const TagAttrib pentaxLensDataAttribs[] = {
 {0, 1, 0, 0, 10, "NominalMinMaxAperture", &paNominalMinMaxApertureInterpreter},
 {0, 1, 0, 0, 14, "MaxAperture", &paMaxApertureInterpreter},
 {-1, 0, 0,  0, 0, "", NULL}};

const TagAttrib pentaxCameraSettingsAttribs[] = {
 {0, 1, 0, 0,  0, "PictureMode2", &paPictureMode2Interpreter},
 {0, 1, 0, 0,  1, "ProgramLine", &paProgramLineInterpreter},
 {0, 1, 0, 0,  1, "EVSteps", &paEVStepsInterpreter},
 {0, 1, 0, 0,  1, "E-DialinProgram", &paEDialinInterpreter},
 {0, 1, 0, 0,  1, "ApertureRing", &paApertureRingUseInterpreter},
 {0, 1, 0, 0,  2, "FlashOptions", &paFlashOptionInterpreter},
 {0, 1, 0, 0,  2, "MeteringMode2", &paMeteringMode2Interpreter},
 {0, 1, 0, 0,  3, "AFMode", &paAFModeInterpreter},
 {0, 1, 0, 0,  4, "AFPointSelected2", &paAFPointSelectedInterpreter},
 {0, 1, 0, 0,  7, "DriveMode2", &paDriveMode2Interpreter},
 {0, 1, 0, 0,  8, "ExposureBracketStepSize", &paExposureBracketStepSizeInterpreter},
 {0, 1, 0, 0,  9, "BracketShotNumber", &stdInterpreter},
 {0, 1, 0, 0, 10, "WhiteBalanceSet", &stdInterpreter},
 {-1, 0, 0,  0, 0, "", NULL}};

const TagAttrib pentaxAEInfoAttribs[] = {
 {0, 1, 0, 0,  0, "AEExposureTime", &stdInterpreter},
 {0, 1, 0, 0,  1, "AEAperture", &stdInterpreter},
 {0, 1, 0, 0,  2, "AE_ISO", &stdInterpreter},
 {0, 1, 0, 0,  3, "AEXv", &stdInterpreter},
 {0, 1, 0, 0,  4, "AEBXv", &stdInterpreter},
 {0, 1, 0, 0,  5, "AEMinExposureTime", &stdInterpreter},
 {0, 1, 0, 0,  6, "AEProgramMode", &stdInterpreter},
 {0, 1, 0, 0,  9, "AEMaxAperture", &stdInterpreter},
 {0, 1, 0, 0, 10, "AEMaxAperture2", &stdInterpreter},
 {0, 1, 0, 0, 11, "AEMinAperture", &stdInterpreter},
 {0, 1, 0, 0, 12, "AEMeteringMode", &stdInterpreter},
 {0, 1, 0, 0, 14, "FlashExposureCompSet", &stdInterpreter},
 {-1, 0, 0,  0, 0, "", NULL}};

const TagAttrib pentaxFlashInfoAttribs[] = {
 {0, 1, 0, 0,  0, "FlashStatus", &paFlashStatusInterpreter},
 {0, 1, 0, 0,  1, "InternalFlashMode", &paInternalFlashModeInterpreter},
 {0, 1, 0, 0,  2, "ExternalFlashMode", &paExternalFlashModeInterpreter},
 {0, 1, 0, 0,  3, "InternalFlashStrength", &stdInterpreter},
 {0, 1, 0, 0, 24, "ExternalFlashGuideNumber", &paExternalFlashGNInterpreter},
 {0, 1, 0, 0, 25, "ExternalFlashExposureComp", &paExternalFlashExposureCompInterpreter},
 {0, 1, 0, 0, 26, "ExternalFlashBounce", &paExternalFlashBounceInterpreter},
 {-1, 0, 0,  0, 0, "", NULL}};

const TagAttrib pentaxBatteryInfoAttribs[] = {
 {0, 1, 0, 0,  0, "PowerSource", &paPowerSourceInterpreter},
 {0, 1, 0, 0,  1, "BatteryStates", &stdInterpreter},
 {0, 1, 0, 0,  2, "BatteryADBodyNoLoad", &stdInterpreter},
 {0, 1, 0, 0,  3, "BatteryADBodyLoad", &stdInterpreter},
 {0, 1, 0, 0,  4, "BatteryADGripNoLoad", &stdInterpreter},
 {0, 1, 0, 0,  5, "BatteryADGripLoad", &stdInterpreter},
 {-1, 0, 0,  0, 0, "", NULL}};

const TagAttrib pentaxCameraInfoAttribs[] = {
 {0, 1, 0, 0,  0, "PentaxModelID", &stdInterpreter},
 {0, 1, 0, 0,  1, "ManufactureDate", &stdInterpreter},
 {0, 1, 0, 0,  2, "ProductionCode", &stdInterpreter},
 {0, 1, 0, 0,  4, "InternalSerialNumber", &stdInterpreter},
 {-1, 0, 0,  0, 0, "", NULL}};

};
#endif




