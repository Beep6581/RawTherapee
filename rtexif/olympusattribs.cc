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
#ifndef _OLYMPUSATTRIBS_
#define _OLYMPUSATTRIBS_

#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "rtexif.h"

namespace rtexif {

class OLOnOffInterpreter : public Interpreter {
    public:
        OLOnOffInterpreter () {}
        virtual std::string toString (Tag* t) {
            if (t->toInt()==0)
                return "Off";
            else
                return "On";
        }
};
OLOnOffInterpreter olOnOffInterpreter;

class OLYesNoInterpreter : public Interpreter {
    public:
        OLYesNoInterpreter () {}
        virtual std::string toString (Tag* t) {
            if (t->toInt()==0)
                return "No";
            else
                return "Yes";
        }
};
OLYesNoInterpreter olYesNoInterpreter;    

class OLApertureInterpreter : public Interpreter {
    public:
        OLApertureInterpreter () {}
        virtual std::string toString (Tag* t) {
            std::ostringstream str;
            str.precision(2);
            str << pow(2, t->toInt() / 512.0);
            return str.str();
        }
};
OLApertureInterpreter olApertureInterpreter;    

class OLLensTypeInterpreter : public Interpreter {
        std::map<std::string, std::string> lenses;
    public:
        OLLensTypeInterpreter () {  // From EXIFTOOL database 'Olympus.pm' V2.09
            // exadecimal bytes
            lenses["00 01 00"] = "Zuiko Digital ED 50mm f/2 Macro";
            lenses["00 01 01"] = "Zuiko Digital 40-150mm f/3.5-4.5";
            lenses["00 01 10"] = "Zuiko Digital ED 14-42mm f/3.5-5.6";
            lenses["00 02 00"] = "Zuiko Digital ED 150mm f/2";
            lenses["00 02 10"] = "Zuiko Digital 17mm f/2.8 Pancake";
            lenses["00 03 00"] = "Zuiko Digital ED 300mm f/2.8";
            lenses["00 03 10"] = "Zuiko Digital ED 14-150mm f/4-5.6";
            lenses["00 04 10"] = "Zuiko Digital ED 9-18mm f/4-5.6";
            lenses["00 05 00"] = "Zuiko Digital 14-54mm f/2.8-3.5";
            lenses["00 05 01"] = "Zuiko Digital Pro ED 90-250mm f/2.8";
            lenses["00 05 10"] = "Zuiko Digital ED 14-42mm f/3.5-5.6 L";
            lenses["00 06 00"] = "Zuiko Digital ED 50-200mm f/2.8-3.5";
            lenses["00 06 01"] = "Zuiko Digital ED 8mm f/3.5 Fisheye";
            lenses["00 06 10"] = "Zuiko Digital ED 40-150mm f/4-5.6";
            lenses["00 07 00"] = "Zuiko Digital 11-22mm f/2.8-3.5";
            lenses["00 07 01"] = "Zuiko Digital 18-180mm f/3.5-6.3";
            lenses["00 07 10"] = "Zuiko Digital ED 12mm f/2";
            lenses["00 08 01"] = "Zuiko Digital 70-300mm f/4-5.6";
            lenses["00 08 10"] = "Zuiko Digital ED 75-300mm f/4.8-6.7";
            lenses["00 09 10"] = "Zuiko Digital 14-42mm f/3.5-5.6 II";
            lenses["00 10 01"] = "Kenko Tokina Reflex 300mm f/6.3 MF Macro";
            lenses["00 10 10"] = "Zuiko Digital ED 12-50mm f/3.5-6.3 EZ";
            lenses["00 11 10"] = "Zuiko Digital 45mm f/1.8";
            lenses["00 12 10"] = "Zuiko Digital ED 60mm f/2.8 Macro";
            lenses["00 13 10"] = "Zuiko Digital ED 14-42mm f/3.5-5.6 II R";
            lenses["00 14 10"] = "Zuiko Digital ED 40-150mm f/4-5.6 R";
            lenses["00 15 00"] = "Zuiko Digital ED 7-14mm f/4";
            lenses["00 15 10"] = "Zuiko Digital ED 75mm f/1.8";
            lenses["00 16 10"] = "Zuiko Digital 17mm f/1.8";
            lenses["00 17 00"] = "Zuiko Digital Pro ED 35-100mm f/2";
            lenses["00 18 00"] = "Zuiko Digital 14-45mm f/3.5-5.6";
            lenses["00 18 10"] = "Zuiko Digital ED 75-300mm f/4.8-6.7 II";
            lenses["00 19 10"] = "Zuiko Digital ED 12-40mm f/2.8 Pro";
            lenses["00 20 00"] = "Zuiko Digital 35mm f/3.5 Macro";
            lenses["00 22 00"] = "Zuiko Digital 17.5-45mm f/3.5-5.6";
            lenses["00 23 00"] = "Zuiko Digital ED 14-42mm f/3.5-5.6";
            lenses["00 24 00"] = "Zuiko Digital ED 40-150mm f/4-5.6";
            lenses["00 30 00"] = "Zuiko Digital ED 50-200mm f/2.8-3.5 SWD";
            lenses["00 31 00"] = "Zuiko Digital ED 12-60mm f/2.8-4 SWD";
            lenses["00 32 00"] = "Zuiko Digital ED 14-35mm f/2 SWD";
            lenses["00 33 00"] = "Zuiko Digital 25mm f/2.8";
            lenses["00 34 00"] = "Zuiko Digital ED 9-18mm f/4-5.6";
            lenses["00 35 00"] = "Zuiko Digital 14-54mm f/2.8-3.5 II";
            lenses["01 01 00"] = "Sigma 18-50mm f/3.5-5.6 DC";
            lenses["01 01 10"] = "Sigma 30mm f/2.8 EX DN";
            lenses["01 02 00"] = "Sigma 55-200mm f/4-5.6 DC";
            lenses["01 02 10"] = "Sigma 19mm f/2.8 EX DN";
            lenses["01 03 00"] = "Sigma 18-125mm f/3.5-5.6 DC";
            lenses["01 03 10"] = "Sigma 30mm f/2.8 DN | A";
            lenses["01 04 00"] = "Sigma 18-125mm f/3.5-5.6 DC";
            lenses["01 04 10"] = "Sigma 19mm f/2.8 DN | A";
            lenses["01 05 00"] = "Sigma 30mm f/1.4 EX DC HSM";
            lenses["01 05 10"] = "Sigma 60mm f/2.8 DN | A";
            lenses["01 06 00"] = "Sigma 50-500mm f/4-6.3 EX DG APO HSM RF";
            lenses["01 07 00"] = "Sigma 105mm f/2.8 EX DG Macro";
            lenses["01 08 00"] = "Sigma 150mm f/2.8 EX DG APO HSM Macro";
            lenses["01 09 00"] = "Sigma 18-50mm f/2.8 EX DC Macro";
            lenses["01 10 00"] = "Sigma 24mm f/1.8 EX DG Aspherical Macro";
            lenses["01 11 00"] = "Sigma 135-400mm f/4.5-5.6 DG APO";
            lenses["01 12 00"] = "Sigma 300-800mm f/5.6 EX DG APO HSM";
            lenses["01 13 00"] = "Sigma 30mm f/1.4 EX DC HSM";
            lenses["01 14 00"] = "Sigma 50-500mm f/4-6.3 EX DG APO HSM";
            lenses["01 15 00"] = "Sigma 10-20mm f/4-5.6 EX DC HSM";
            lenses["01 16 00"] = "Sigma 70-200mm f/2.8 II EX DG APO HSM Macro";
            lenses["01 17 00"] = "Sigma 50mm f/1.4 EX DG HSM";
            lenses["02 01 00"] = "Leica D Vario Elmarit 14-50mm f/2.8-3.5 Asph.";
            lenses["02 01 10"] = "Lumix G Vario 14-45mm f/3.5-5.6 Asph. Mega OIS";
            lenses["02 02 00"] = "Leica D Summilux 25mm f/1.4 Asph.";
            lenses["02 02 10"] = "Lumix G Vario 45-200mm f/4-5.6 Mega OIS";
            lenses["02 03 00"] = "Leica D Vario Elmar 14-50mm f/3.8-5.6 Asph. Mega OIS";
            lenses["02 03 01"] = "Leica D Vario Elmar 14-50mm f/3.8-5.6 Asph.";
            lenses["02 03 10"] = "Lumix G Vario HD 14-140mm f/4-5.8 Asph. Mega OIS";
            lenses["02 04 00"] = "Leica D Vario Elmar 14-150mm f/3.5-5.6";
            lenses["02 04 10"] = "Lumix G Vario 7-14mm f/4 Asph.";
            lenses["02 05 10"] = "Lumix G 20mm f/1.7 Asph.";
            lenses["02 06 10"] = "Leica DG Macro-Elmarit 45mm f/2.8 Asph. Mega OIS";
            lenses["02 07 10"] = "Lumix G Vario 14-42mm f/3.5-5.6 Asph. Mega OIS";
            lenses["02 08 10"] = "Lumix G Fisheye 8mm f/3.5";
            lenses["02 09 10"] = "Lumix G Vario 100-300mm f/4-5.6 Mega OIS";
            lenses["02 10 10"] = "Lumix G 14mm f/2.5 Asph.";
            lenses["02 11 10"] = "Lumix G 12.5mm f/12 3D";
            lenses["02 12 10"] = "Leica DG Summilux 25mm f/1.4 Asph.";
            lenses["02 13 10"] = "Lumix G X Vario PZ 45-175mm f/4-5.6 Asph. Power OIS";
            lenses["02 14 10"] = "Lumix G X Vario PZ 14-42mm f/3.5-5.6 Asph. Power OIS";
            lenses["02 15 10"] = "Lumix G X Vario 12-35mm f/2.8 Asph. Power OIS";
            lenses["02 16 10"] = "Lumix G Vario 45-150mm f/4-5.6 Asph. Mega OIS";
            lenses["02 17 10"] = "Lumix G X Vario 35-100mm f/2.8 Power OIS";
            lenses["02 18 10"] = "Lumix G Vario 14-42mm f/3.5-5.6 II Asph. Mega OIS";
            lenses["02 19 10"] = "Lumix G Vario 14-140mm f/3.5-5.6 Asph. Power OIS";
            lenses["03 01 00"] = "Leica D Vario Elmarit 14-50mm f/2.8-3.5 Asph.";
            lenses["03 02 00"] = "Leica D Summilux 25mm f/1.4 Asph.";
        }
        virtual std::string toString (Tag* t) {
            std::ostringstream lid;
            lid.setf (std::ios_base::hex, std::ios_base::basefield);
            lid.setf (std::ios_base::uppercase);
            lid << std::setw(2) << std::setfill('0') << t->toInt(0)<< ' '; //maker
            lid << std::setw(2) << std::setfill('0') << t->toInt(2)<< ' '; //model
            lid << std::setw(2) << std::setfill('0') << t->toInt(3); // submodel

            std::map<std::string,std::string>::iterator r = lenses.find (lid.str());
            if (r!=lenses.end())
                return r->second;
            else
            	return "Unknown";
        }
};
OLLensTypeInterpreter olLensTypeInterpreter;    

class OLFlashTypeInterpreter : public ChoiceInterpreter {
    public:
        OLFlashTypeInterpreter () {
            choices[0]      = "None";
            choices[2]      = "Simple E-System";
            choices[3]      = "E-System";
        }
};
OLFlashTypeInterpreter olFlashTypeInterpreter;

class OLExposureModeInterpreter : public ChoiceInterpreter {
    public:
        OLExposureModeInterpreter () {
            choices[1] = "Manual";
            choices[2] = "Program";
            choices[3] = "Aperture-priority AE";
            choices[4] = "Shutter speed priority AE";
            choices[5] = "Program-shift";
        }
};
OLExposureModeInterpreter olExposureModeInterpreter;

class OLMeteringModeInterpreter : public ChoiceInterpreter {
    public:
        OLMeteringModeInterpreter () {
            choices[2] = "Center-weighted average";
            choices[3] = "Spot";
            choices[5] = "ESP";
            choices[261] = "Pattern+AF";
            choices[515] = "Spot+Highlight control";
            choices[1027] = "Spot+Shadow control";
       }
};
OLMeteringModeInterpreter olMeteringModeInterpreter;

class OLFocusModeInterpreter : public ChoiceInterpreter {
    public:
        OLFocusModeInterpreter () {
            choices[0] = "Single AF";
            choices[1] = "Sequential shooting AF";
            choices[2] = "Continuous AF";
            choices[3] = "Multi AF";
            choices[10] = "MF";
        }
};
OLFocusModeInterpreter olFocusModeInterpreter;

class OLWhitebalance2Interpreter : public ChoiceInterpreter {
    public:
        OLWhitebalance2Interpreter () {
            choices[0] = "Auto";
            choices[16] = "7500K (Fine Weather with Shade)";
            choices[17] = "6000K (Cloudy)";
            choices[18] = "5300K (Fine Weather)";
            choices[20] = "3000K (Tungsten light)";
            choices[21] = "3600K (Tungsten light-like)";
            choices[33] = "6600K (Daylight fluorescent)";
            choices[34] = "4500K (Neutral white fluorescent)";
            choices[35] = "4000K (Cool white fluorescent)";
            choices[48] = "3600K (Tungsten light-like)";
            choices[256] = "Custom WB 1";
            choices[257] = "Custom WB 2";
            choices[258] = "Custom WB 3";
            choices[259] = "Custom WB 4";
            choices[512] = "Custom WB 5400K";
            choices[513] = "Custom WB 2900K";
            choices[514] = "Custom WB 8000K";
        }
};
OLWhitebalance2Interpreter olWhitebalance2Interpreter;

class OLSceneModeInterpreter : public ChoiceInterpreter {
    public:
        OLSceneModeInterpreter () {
            choices[0] = "Standard";
            choices[6] = "Auto";
            choices[7] = "Sport";
            choices[8] = "Portrait";
            choices[9] = "Landscape+Portrait";
            choices[10] = "Landscape";
            choices[11] = "Night Scene";
            choices[12] = "Self Portrait";
            choices[13] = "Panorama";
            choices[14] = "2 in 1";
            choices[15] = "Movie";
            choices[16] = "Landscape+Portrait";
            choices[17] = "Night+Portrait";
            choices[18] = "Indoor";
            choices[19] = "Fireworks";
            choices[20] = "Sunset";
            choices[22] = "Macro";
            choices[23] = "Super Macro";
            choices[24] = "Food";
            choices[25] = "Documents";
            choices[26] = "Museum";
            choices[27] = "Shoot & Select";
            choices[28] = "Beach & Snow";
            choices[29] = "Self Protrait+Timer";
            choices[30] = "Candle";
            choices[31] = "Available Light";
            choices[32] = "Behind Glass";
            choices[33] = "My Mode";
            choices[34] = "Pet";
            choices[35] = "Underwater Wide1";
            choices[36] = "Underwater Macro";
            choices[37] = "Shoot & Select1";
            choices[38] = "Shoot & Select2";
            choices[39] = "High Key";
            choices[40] = "Digital Image Stabilization";
            choices[41] = "Auction";
            choices[42] = "Beach";
            choices[43] = "Snow";
            choices[44] = "Underwater Wide2";
            choices[45] = "Low Key";
            choices[46] = "Children";
            choices[47] = "Vivid";
            choices[48] = "Nature Macro";
            choices[49] = "Underwater Snapshot";
            choices[50] = "Shooting Guide";
        }
};
OLSceneModeInterpreter olSceneModeInterpreter;

class OLPictureModeBWFilterInterpreter : public ChoiceInterpreter {
    public:
        OLPictureModeBWFilterInterpreter () {
            choices[0] = "n/a";
            choices[1] = "Neutral";
            choices[2] = "Yellow";
            choices[3] = "Orange";
            choices[4] = "Red";
            choices[5] = "Green";
        }
};
OLPictureModeBWFilterInterpreter olPictureModeBWFilterInterpreter;

class OLPictureModeToneInterpreter : public ChoiceInterpreter {
    public:
        OLPictureModeToneInterpreter () {
            choices[0] = "n/a";
            choices[1] = "Neutral";
            choices[2] = "Sepia";
            choices[3] = "Blue";
            choices[4] = "Purple";
            choices[5] = "Green";
        }
};
OLPictureModeToneInterpreter olPictureModeToneInterpreter;

class OLImageQuality2Interpreter : public ChoiceInterpreter {
    public:
        OLImageQuality2Interpreter () {
            choices[1] = "SQ";
            choices[2] = "HQ";
            choices[3] = "SHQ";
            choices[4] = "RAW";
        }
};
OLImageQuality2Interpreter olImageQuality2Interpreter;

class OLDevEngineInterpreter : public ChoiceInterpreter {
    public:
        OLDevEngineInterpreter () {
            choices[0] = "High Speed";
            choices[1] = "High Function";
            choices[2] = "Advanced High Speed";
            choices[3] = "Advanced High Function";
        }
};
OLDevEngineInterpreter olDevEngineInterpreter;

class OLPictureModeInterpreter : public ChoiceInterpreter {
    public:
        OLPictureModeInterpreter () {
            choices[1] = "Vivid";
            choices[2] = "Natural";
            choices[3] = "Muted";
            choices[4] = "Portrait";
            choices[256] = "Monotone";
            choices[512] = "Sepia";
        }
};
OLPictureModeInterpreter olPictureModeInterpreter;

class OLColorSpaceInterpreter : public ChoiceInterpreter {
    public:
        OLColorSpaceInterpreter () {
            choices[0] = "sRGB";
            choices[1] = "Adobe RGB";
            choices[2] = "Pro Photo RGB";
        }
};
OLColorSpaceInterpreter olColorSpaceInterpreter;

class OLNoiseFilterInterpreter : public Interpreter {
    public:
        OLNoiseFilterInterpreter () {}
        virtual std::string toString (Tag* t) {
            int a = t->toInt (0);
            int b = t->toInt (2);
            int c = t->toInt (4);
            if (a==-1 && b==-2 && c==1)
                return "Low";
            else if (a==-2 && b==-2 && c==1)
                return "Off";
            else if (a==0 && b==-2 && c==1)
                return "Standard";
            else if (a==1 && b==-2 && c==1)
                return "High";
            else return "Unknown";
        }
};
OLNoiseFilterInterpreter olNoiseFilterInterpreter;

class OLFlashModeInterpreter : public Interpreter {
    public:
        OLFlashModeInterpreter () {}
        virtual std::string toString (Tag* t) {
            std::ostringstream str;
            int a = t->toInt ();
            str << "Flash Used = " << ((a&1) ? "Yes" : "No") << std::endl;
            str << "Fill-in = " << ((a&2) ? "On" : "Off") << std::endl;
            str << "Red-eye = " << ((a&4) ? "On" : "Off") << std::endl;
            str << "Slow-sync = " << ((a&8) ? "On" : "Off") << std::endl;
            str << "Forced On = " << ((a&16) ? "On" : "Off") << std::endl;
            str << "2nd Curtain = " << ((a&32) ? "On" : "Off");
            return str.str();
        }
};
OLFlashModeInterpreter olFlashModeInterpreter;

class OLNoiseReductionInterpreter : public Interpreter {
    public:
        OLNoiseReductionInterpreter () {}
        virtual std::string toString (Tag* t) {
            std::ostringstream str;
            int a = t->toInt ();
            str << "Noise Reduction = " << ((a&1) ? "On" : "Off") << std::endl;
            str << "Noise Filter = " << ((a&2) ? "On" : "Off") << std::endl;
            str << "Noise Filter (ISO Boost) = " << ((a&4) ? "On" : "Off");
            return str.str();
        }
};
OLNoiseReductionInterpreter olNoiseReductionInterpreter;

class OLFlashModelInterpreter : public ChoiceInterpreter {
    public:
        OLFlashModelInterpreter () {
            choices[0]      = "None";
            choices[1]      = "FL-20";
            choices[2]      = "FL-50";
            choices[3]      = "RF-11";
            choices[4]      = "TF-22";
            choices[5]      = "FL-36";
            choices[6]      = "FL-50R";
            choices[7]      = "FL-36R";
        }
};
OLFlashModelInterpreter olFlashModelInterpreter;

const TagAttrib olyFocusInfoAttribs[] = {
 {0, AC_WRITE, 0, 0, 0x0000, AUTO, "FocusInfoVersion", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0209, AUTO, "AutoFocus", &olOnOffInterpreter},
 {0, AC_WRITE, 0, 0, 0x0210, AUTO, "SceneDetect", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0211, AUTO, "SceneArea", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0212, AUTO, "SceneDetectData", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0300, AUTO, "ZoomStepCount", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0301, AUTO, "FocusStepCount", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0303, AUTO, "FocusStepInfinity", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0304, AUTO, "FocusStepNear", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0305, AUTO, "FocusDistance", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0308, AUTO, "AFPoint", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x1201, AUTO, "ExternalFlash", &olOnOffInterpreter},
 {0, AC_WRITE, 0, 0, 0x1203, AUTO, "ExternalFlashGuideNumber", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x1204, AUTO, "ExternalFlashBounce", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x1205, AUTO, "ExternalFlashZoom", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x1208, AUTO, "InternalFlash", &olOnOffInterpreter},
 {0, AC_WRITE, 0, 0, 0x1209, AUTO, "ManualFlash", &olOnOffInterpreter},
 {0, AC_WRITE, 0, 0, 0x1500, AUTO, "SensorTemperature", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x1600, AUTO, "ImageStabilization", &stdInterpreter},
 {-1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}};

const TagAttrib olyImageProcessingAttribs[] = {
 {0, AC_WRITE, 0, 0, 0x0000, AUTO, "ImageProcessingVersion", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0100, AUTO, "WB_RBLevels", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0102, AUTO, "WB_RBLevels3000K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0103, AUTO, "WB_RBLevels3300K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0104, AUTO, "WB_RBLevels3600K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0105, AUTO, "WB_RBLevels3900K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0106, AUTO, "WB_RBLevels4000K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0107, AUTO, "WB_RBLevels4300K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0108, AUTO, "WB_RBLevels4500K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0109, AUTO, "WB_RBLevels4800K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x010a, AUTO, "WB_RBLevels5300K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x010b, AUTO, "WB_RBLevels6000K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x010c, AUTO, "WB_RBLevels6600K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x010d, AUTO, "WB_RBLevels7500K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x010e, AUTO, "WB_RBLevelsCWB1", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x010f, AUTO, "WB_RBLevelsCWB2", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0110, AUTO, "WB_RBLevelsCWB3", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0111, AUTO, "WB_RBLevelsCWB4", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0113, AUTO, "WB_GLevel3000K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0114, AUTO, "WB_GLevel3300K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0115, AUTO, "WB_GLevel3600K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0116, AUTO, "WB_GLevel3900K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0117, AUTO, "WB_GLevel4000K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0118, AUTO, "WB_GLevel4300K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0119, AUTO, "WB_GLevel4500K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x011a, AUTO, "WB_GLevel4800K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x011b, AUTO, "WB_GLevel5300K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x011c, AUTO, "WB_GLevel6000K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x011d, AUTO, "WB_GLevel6600K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x011e, AUTO, "WB_GLevel7500K", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x011f, AUTO, "WB_GLevel", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0200, AUTO, "ColorMatrix", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0300, AUTO, "Enhancer", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0301, AUTO, "EnhancerValues", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0310, AUTO, "CoringFilter", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0311, AUTO, "CoringValues", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0600, AUTO, "BlackLevel2", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0610, AUTO, "GainBase", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0611, AUTO, "ValidBits", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0612, AUTO, "CropLeft", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0613, AUTO, "CropTop", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0614, AUTO, "CropWidth", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0615, AUTO, "CropHeight", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x1010, AUTO, "NoiseReduction2", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x1011, AUTO, "DistortionCorrection2", &olOnOffInterpreter},
 {0, AC_WRITE, 0, 0, 0x1012, AUTO, "ShadingCompensation2", &olOnOffInterpreter},
 {1, AC_WRITE, 0, 0, 0x1103, AUTO, "UnknownBlock", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x1200, AUTO, "FaceDetect", &olOnOffInterpreter},
 {0, AC_WRITE, 0, 0, 0x1201, AUTO, "FaceDetectArea", &stdInterpreter},
 {-1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}};

const TagAttrib olyRawDevelopmentAttribs[] = {
 {0, AC_WRITE, 0, 0, 0x0000, AUTO, "RawDevVersion", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0100, AUTO, "RawDevExposureBiasValue", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0101, AUTO, "RawDevWhiteBalanceValue", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0102, AUTO, "RawDevWBFineAdjustment", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0103, AUTO, "RawDevGrayPoint", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0104, AUTO, "RawDevSaturationEmphasis", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0105, AUTO, "RawDevMemoryColorEmphasis", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0106, AUTO, "RawDevContrastValue", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0107, AUTO, "RawDevSharpnessValue", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0108, AUTO, "RawDevColorSpace", &olColorSpaceInterpreter},
 {0, AC_WRITE, 0, 0, 0x0109, AUTO, "RawDevEngine", &olDevEngineInterpreter},
 {0, AC_WRITE, 0, 0, 0x010a, AUTO, "RawDevNoiseReduction", &olNoiseReductionInterpreter},
 {0, AC_WRITE, 0, 0, 0x010b, AUTO, "RawDevEditStatus", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x010c, AUTO, "RawDevSettings", &stdInterpreter},
 {-1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}};

const TagAttrib olyRawDevelopment2Attribs[] = {
 {0, AC_WRITE, 0, 0, 0x0000, AUTO, "RawDevVersion", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0100, AUTO, "RawDevExposureBiasValue", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0101, AUTO, "RawDevWhiteBalance", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0102, AUTO, "RawDevWhiteBalanceValue", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0103, AUTO, "RawDevWBFineAdjustment", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0104, AUTO, "RawDevGrayPoint", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0105, AUTO, "RawDevContrastValue", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0106, AUTO, "RawDevSharpnessValue", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0107, AUTO, "RawDevSaturationEmphasis", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0108, AUTO, "RawDevMemoryColorEmphasis", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0109, AUTO, "RawDevColorSpace", &olColorSpaceInterpreter},
 {0, AC_WRITE, 0, 0, 0x010a, AUTO, "RawDevNoiseReduction", &olNoiseReductionInterpreter},
 {0, AC_WRITE, 0, 0, 0x010b, AUTO, "RawDevEngine", &olDevEngineInterpreter},
 {0, AC_WRITE, 0, 0, 0x010c, AUTO, "RawDevPictureMode", &olPictureModeInterpreter},
 {0, AC_WRITE, 0, 0, 0x010d, AUTO, "RawDevPMSaturation", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x010e, AUTO, "RawDevPMContrast", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x010f, AUTO, "RawDevPMSharpness", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0110, AUTO, "RawDevPM_BWFilter", &olPictureModeBWFilterInterpreter},
 {0, AC_WRITE, 0, 0, 0x0111, AUTO, "RawDevPMPictureTone", &olPictureModeToneInterpreter},
 {0, AC_WRITE, 0, 0, 0x0112, AUTO, "RawDevGradation", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0113, AUTO, "RawDevSaturation3", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0119, AUTO, "RawDevAutoGradation", &olOnOffInterpreter},
 {0, AC_WRITE, 0, 0, 0x0120, AUTO, "RawDevPMNoiseFilter", &stdInterpreter},
 {-1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}};

const TagAttrib olyCameraSettingsAttribs[] = {
 {0, AC_WRITE, 0, 0, 0x0000, AUTO, "CameraSettingsVersion", &stdInterpreter},
 {1, AC_WRITE, 0, 0, 0x0100, AUTO, "PreviewImageValid", &olYesNoInterpreter},
 {1, AC_WRITE, 0, 0, 0x0101, AUTO, "PreviewImageStart", &stdInterpreter},
 {1, AC_WRITE, 0, 0, 0x0102, AUTO, "PreviewImageLength", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0200, AUTO, "ExposureMode", &olExposureModeInterpreter},
 {0, AC_WRITE, 0, 0, 0x0201, AUTO, "AELock", &olOnOffInterpreter},
 {0, AC_WRITE, 0, 0, 0x0202, AUTO, "MeteringMode", &olMeteringModeInterpreter},
 {0, AC_WRITE, 0, 0, 0x0300, AUTO, "MacroMode", &olOnOffInterpreter},
 {0, AC_WRITE, 0, 0, 0x0301, AUTO, "FocusMode", &olFocusModeInterpreter},
 {0, AC_WRITE, 0, 0, 0x0302, AUTO, "FocusProcess", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0303, AUTO, "AFSearch", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0304, AUTO, "AFAreas", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0400, AUTO, "FlashMode", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0401, AUTO, "FlashExposureComp", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0500, AUTO, "WhiteBalance2", &olWhitebalance2Interpreter},
 {0, AC_WRITE, 0, 0, 0x0501, AUTO, "WhiteBalanceTemperature", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0502, AUTO, "WhiteBalanceBracket", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0503, AUTO, "CustomSaturation", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0504, AUTO, "ModifiedSaturation", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0505, AUTO, "ContrastSetting", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0506, AUTO, "SharpnessSetting", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0507, AUTO, "ColorSpace", &olColorSpaceInterpreter},
 {0, AC_WRITE, 0, 0, 0x0509, AUTO, "SceneMode", &olSceneModeInterpreter},
 {0, AC_WRITE, 0, 0, 0x050a, AUTO, "NoiseReduction", &olNoiseReductionInterpreter},
 {0, AC_WRITE, 0, 0, 0x050b, AUTO, "DistortionCorrection", &olOnOffInterpreter},
 {0, AC_WRITE, 0, 0, 0x050c, AUTO, "ShadingCompensation", &olOnOffInterpreter},
 {0, AC_WRITE, 0, 0, 0x050d, AUTO, "CompressionFactor", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x050f, AUTO, "Gradation", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0520, AUTO, "PictureMode", &olPictureModeInterpreter},
 {0, AC_WRITE, 0, 0, 0x0521, AUTO, "PictureModeSaturation", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0522, AUTO, "PictureModeHue", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0523, AUTO, "PictureModeContrast", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0524, AUTO, "PictureModeSharpness", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0525, AUTO, "PictureModeBWFilter", &olPictureModeBWFilterInterpreter},
 {0, AC_WRITE, 0, 0, 0x0526, AUTO, "PictureModeTone", &olPictureModeToneInterpreter},
 {0, AC_WRITE, 0, 0, 0x0527, AUTO, "NoiseFilter", &olNoiseFilterInterpreter},
 {0, AC_WRITE, 0, 0, 0x0600, AUTO, "DriveMode", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0601, AUTO, "PanoramaMode", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0603, AUTO, "ImageQuality2", &olImageQuality2Interpreter},
 {0, AC_WRITE, 0, 0, 0x0900, AUTO, "ManometerPressure", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0901, AUTO, "ManometerReading", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0902, AUTO, "ExtendedWBDetect", &olOnOffInterpreter},
 {-1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}};

const TagAttrib olyEquipmentAttribs[] = {
 {0, AC_WRITE, 0, 0, 0x0000, AUTO, "EquipmentVersion", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0100, AUTO, "CameraType2", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0101, AUTO, "SerialNumber", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0102, AUTO, "InternalSerialNumber", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0103, AUTO, "FocalPlaneDiagonal", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0104, AUTO, "BodyFirmwareVersion", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0201, AUTO, "LensType", &olLensTypeInterpreter},
 {0, AC_WRITE, 0, 0, 0x0202, AUTO, "LensSerialNumber", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0204, AUTO, "LensFirmwareVersion", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0205, AUTO, "MaxApertureAtMinFocal", &olApertureInterpreter},
 {0, AC_WRITE, 0, 0, 0x0206, AUTO, "MaxApertureAtMaxFocal", &olApertureInterpreter},
 {0, AC_WRITE, 0, 0, 0x0207, AUTO, "MinFocalLength", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0208, AUTO, "MaxFocalLength", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x020a, AUTO, "MaxApertureAtCurrentFocal", &olApertureInterpreter},
 {0, AC_WRITE, 0, 0, 0x020b, AUTO, "LensProperties", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0301, AUTO, "Extender", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0302, AUTO, "ExtenderSerialNumber", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0303, AUTO, "ExtenderModel", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x0304, AUTO, "ExtenderFirmwareVersion", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x1000, AUTO, "FlashType", &olFlashTypeInterpreter},
 {0, AC_WRITE, 0, 0, 0x1001, AUTO, "FlashModel", &olFlashModelInterpreter},
 {0, AC_WRITE, 0, 0, 0x1002, AUTO, "FlashFirmwareVersion", &stdInterpreter},
 {0, AC_WRITE, 0, 0, 0x1003, AUTO, "FlashSerialNumber", &stdInterpreter},
 {-1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}};

const TagAttrib olympusAttribs[] = {
 {0, AC_WRITE,  0, 0, 0x0104, AUTO, "BodyFirmwareVersion", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x0200, AUTO, "SpecialMode", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x0201, AUTO, "Quality", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x0202, AUTO, "Macro", &olOnOffInterpreter},
 {0, AC_WRITE,  0, 0, 0x0203, AUTO, "BWMode", &olOnOffInterpreter},
 {0, AC_WRITE,  0, 0, 0x0204, AUTO, "DigitalZoom", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x0205, AUTO, "FocalPlaneDiagonal", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x0206, AUTO, "LensDistortionParams", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x0207, AUTO, "CameraType", &stdInterpreter},
 {1, AC_WRITE,  0, 0, 0x0208, AUTO, "TextInfo", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x0209, AUTO, "CameraID", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x020b, AUTO, "EpsonImageWidth", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x020c, AUTO, "EpsonImageHeight", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x020d, AUTO, "EpsonSoftware", &stdInterpreter},
 {0, AC_SYSTEM, 0, 0, 0x0280, AUTO, "PreviewImage", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x0300, AUTO, "PreCaptureFrames", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x0301, AUTO, "WhiteBoard", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x0302, AUTO, "OneTouchWB", &olOnOffInterpreter},
 {0, AC_WRITE,  0, 0, 0x0303, AUTO, "WhiteBalanceBracket", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x0304, AUTO, "WhiteBalanceBias", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x0403, AUTO, "SceneMode", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x0404, AUTO, "SerialNumber", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x0405, AUTO, "Firmware", &stdInterpreter},
 {1, AC_WRITE,  0, 0, 0x0e00, AUTO, "PrintIM", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x0f00, AUTO, "DataDump", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x0f01, AUTO, "DataDump2", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1000, AUTO, "ShutterSpeedValue", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1001, AUTO, "ISOValue", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1002, AUTO, "ApertureValue", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1003, AUTO, "BrightnessValue", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1004, AUTO, "FlashMode", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1005, AUTO, "FlashDevice", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1006, AUTO, "ExposureCompensation", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1007, AUTO, "SensorTemperature", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1008, AUTO, "LensTemperature", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1009, AUTO, "LightCondition", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x100a, AUTO, "FocusRange", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x100b, AUTO, "FocusMode", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x100c, AUTO, "ManualFocusDistance", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x100d, AUTO, "ZoomStepCount", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x100e, AUTO, "FocusStepCount", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x100f, AUTO, "Sharpness", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1010, AUTO, "FlashChargeLevel", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1011, AUTO, "ColorMatrix", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1012, AUTO, "BlackLevel", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1013, AUTO, "ColorTemperatureBG", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1014, AUTO, "ColorTemperatureRG", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1015, AUTO, "WBMode", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1017, AUTO, "RedBalance", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1018, AUTO, "BlueBalance", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1019, AUTO, "ColorMatrixNumber", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x101a, AUTO, "SerialNumber", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x101b, AUTO, "ExternalFlashAE1_0", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x101c, AUTO, "ExternalFlashAE2_0", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x101d, AUTO, "InternalFlashAE1_0", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x101e, AUTO, "InternalFlashAE2_0", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x101f, AUTO, "ExternalFlashAE1", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1020, AUTO, "ExternalFlashAE2", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1021, AUTO, "InternalFlashAE1", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1022, AUTO, "InternalFlashAE2", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1023, AUTO, "FlashExposureComp", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1024, AUTO, "InternalFlashTable", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1025, AUTO, "ExternalFlashGValue", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1026, AUTO, "ExternalFlashBounce", &olYesNoInterpreter},
 {0, AC_WRITE,  0, 0, 0x1027, AUTO, "ExternalFlashZoom", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1028, AUTO, "ExternalFlashMode", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1029, AUTO, "Contrast", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x102a, AUTO, "SharpnessFactor", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x102b, AUTO, "ColorControl", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x102c, AUTO, "ValidBits", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x102d, AUTO, "CoringFilter", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x102e, AUTO, "OlympusImageWidth", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x102f, AUTO, "OlympusImageHeight", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1030, AUTO, "SceneDetect", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1031, AUTO, "SceneArea", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1033, AUTO, "SceneDetectData", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1034, AUTO, "CompressionRatio", &stdInterpreter},
 {1, AC_WRITE,  0, 0, 0x1035, AUTO, "PreviewImageValid", &olYesNoInterpreter},
 {1, AC_WRITE,  0, 0, 0x1036, AUTO, "PreviewImageStart", &stdInterpreter},
 {1, AC_WRITE,  0, 0, 0x1037, AUTO, "PreviewImageLength", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1038, AUTO, "AFResult", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x1039, AUTO, "CCDScanMode", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x103a, AUTO, "NoiseReduction", &olOnOffInterpreter},
 {0, AC_WRITE,  0, 0, 0x103b, AUTO, "InfinityLensStep", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x103c, AUTO, "NearLensStep", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x103d, AUTO, "LightValueCenter", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x103e, AUTO, "LightValuePeriphery", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x103f, AUTO, "FieldCount", &stdInterpreter},
 {0, AC_WRITE,  0, olyEquipmentAttribs, 0x2010, AUTO, "Equipment", &stdInterpreter},
 {0, AC_WRITE,  0, olyCameraSettingsAttribs, 0x2020, AUTO, "CameraSettings", &stdInterpreter},
 {0, AC_WRITE,  0, olyRawDevelopmentAttribs, 0x2030, AUTO, "RawDevelopment", &stdInterpreter},
 {0, AC_WRITE,  0, olyRawDevelopment2Attribs, 0x2031, AUTO, "RawDev2", &stdInterpreter},
 {0, AC_WRITE,  0, olyImageProcessingAttribs, 0x2040, AUTO, "ImageProcessing", &stdInterpreter},
 {0, AC_WRITE,  0, olyFocusInfoAttribs, 0x2050, AUTO, "FocusInfo", &stdInterpreter},
 {1, AC_WRITE,  0, 0, 0x2100, AUTO, "Olympus2100", &stdInterpreter},
 {1, AC_WRITE,  0, 0, 0x2300, AUTO, "Olympus2300", &stdInterpreter},
 {1, AC_WRITE,  0, 0, 0x2400, AUTO, "Olympus2400", &stdInterpreter},
 {1, AC_WRITE,  0, 0, 0x2500, AUTO, "Olympus2500", &stdInterpreter},
 {1, AC_WRITE,  0, 0, 0x2600, AUTO, "Olympus2600", &stdInterpreter},
 {1, AC_WRITE,  0, 0, 0x2700, AUTO, "Olympus2700", &stdInterpreter},
 {1, AC_WRITE,  0, 0, 0x2800, AUTO, "Olympus2800", &stdInterpreter},
 {1, AC_WRITE,  0, 0, 0x2900, AUTO, "Olympus2900", &stdInterpreter},
 {0, AC_WRITE,  0, 0, 0x3000, AUTO, "RawInfo", &stdInterpreter},
 {-1, AC_DONTWRITE, 0,  0, 0, AUTO, "", NULL}};
}
#endif

