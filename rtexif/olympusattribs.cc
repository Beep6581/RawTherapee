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

#include <rtexif.h>
#include <string>
#include <map>
#include <math.h>
#include <sstream>
#include <iomanip>

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
        std::map<int, std::string> lenses;
    public:
        OLLensTypeInterpreter () {
            lenses[1] = "Zuiko Digital ED 50mm F2.0 Macro";
            lenses[1 +65536] = "Zuiko Digital 40-150mm F3.5-4.5"; 
            lenses[2] = "Zuiko Digital ED 150mm F2.0";
            lenses[3] = "Zuiko Digital ED 300mm F2.8";
            lenses[5] = "Zuiko Digital 14-54mm F2.8-3.5";
            lenses[5 +65536] = "Zuiko Digital Pro ED 90-250mm F2.8"; 
            lenses[6] = "Zuiko Digital ED 50-200mm F2.8-3.5";
            lenses[6 +65536] = "Zuiko Digital ED 8mm F3.5 Fisheye"; 
            lenses[7] = "Zuiko Digital 11-22mm F2.8-3.5";
            lenses[7 +65536] = "Zuiko Digital 18-180mm F3.5-6.3"; 
            lenses[8] = "Zuiko Digital 70-300mm F4.0-5.6"; 
            lenses[21] = "Zuiko Digital ED 7-14mm F4.0";
            lenses[23] = "Zuiko Digital Pro ED 35-100mm F2.0"; 
            lenses[24] = "Zuiko Digital 14-45mm F3.5-5.6";
            lenses[32] = "Zuiko Digital 35mm F3.5 Macro"; 
            lenses[34] = "Zuiko Digital 17.5-45mm F3.5-5.6"; 
            lenses[35] = "Zuiko Digital ED 14-42mm F3.5-5.6"; 
            lenses[36] = "Zuiko Digital ED 40-150mm F4.0-5.6"; 
            lenses[48] = "Zuiko Digital ED 50-200mm SWD F2.8-3.5"; 
            lenses[49] = "Zuiko Digital ED 12-60mm SWD F2.8-4.0"; 
            lenses[256+ 1] = "18-50mm F3.5-5.6"; 
            lenses[256+ 2] = "55-200mm F4.0-5.6 DC";
            lenses[256+ 3] = "18-125mm F3.5-5.6 DC";
            lenses[256+ 4] = "18-125mm F3.5-5.6"; 
            lenses[256+ 5] = "30mm F1.4"; 
            lenses[256+ 6] = "50-500mm F4.0-6.3 EX DG APO HSM RF"; 
            lenses[256+ 7] = "105mm F2.8 DG"; 
            lenses[256+ 8] = "150mm F2.8 DG HSM";
            lenses[256+ 17] = "135-400mm F4.5-5.6 DG ASP APO RF";
            lenses[256+ 18] = "300-800mm F5.6 EX DG APO";
            lenses[512+ 1] = "D Vario Elmarit 14-50mm, F2.8-3.5 Asph.";
            lenses[512+ 2] = "D Summilux 25mm, F1.4 Asph.";
            lenses[512+ 4] = "Vario Elmar 14-150mm f3.5-5.6";
            lenses[768+ 1] = "D Vario Elmarit 14-50mm, F2.8-3.5 Asph.";
            lenses[768+ 2] = "D Summilux 25mm, F1.4 Asph.";
        }
        virtual std::string toString (Tag* t) {
            int make = t->toInt(0);
            int model = t->toInt(2);
            int add = 0;
            if (make==0 && (model==1 || model==5 || model==7 || model==6))
                add += 65536 * t->toInt(3);
            return lenses [256 * make + model + add];
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
 0, 1, 0, 0, 0x0000, "FocusInfoVersion", &stdInterpreter,
 0, 1, 0, 0, 0x0209, "AutoFocus", &olOnOffInterpreter,
 0, 1, 0, 0, 0x0210, "SceneDetect", &stdInterpreter,
 0, 1, 0, 0, 0x0211, "SceneArea", &stdInterpreter,
 0, 1, 0, 0, 0x0212, "SceneDetectData", &stdInterpreter,
 0, 1, 0, 0, 0x0300, "ZoomStepCount", &stdInterpreter,
 0, 1, 0, 0, 0x0301, "FocusStepCount", &stdInterpreter,
 0, 1, 0, 0, 0x0303, "FocusStepInfinity", &stdInterpreter,
 0, 1, 0, 0, 0x0304, "FocusStepNear", &stdInterpreter,
 0, 1, 0, 0, 0x0305, "FocusDistance", &stdInterpreter,
 0, 1, 0, 0, 0x0308, "AFPoint", &stdInterpreter,
 0, 1, 0, 0, 0x1201, "ExternalFlash", &olOnOffInterpreter,
 0, 1, 0, 0, 0x1203, "ExternalFlashGuideNumber", &stdInterpreter,
 0, 1, 0, 0, 0x1204, "ExternalFlashBounce", &stdInterpreter,
 0, 1, 0, 0, 0x1205, "ExternalFlashZoom", &stdInterpreter,
 0, 1, 0, 0, 0x1208, "InternalFlash", &olOnOffInterpreter,
 0, 1, 0, 0, 0x1209, "ManualFlash", &olOnOffInterpreter,
 0, 1, 0, 0, 0x1500, "SensorTemperature", &stdInterpreter,
 0, 1, 0, 0, 0x1600, "ImageStabilization", &stdInterpreter,
-1, 0, 0,  0, 0, "", NULL};

const TagAttrib olyImageProcessingAttribs[] = {
 0, 1, 0, 0, 0x0000, "ImageProcessingVersion", &stdInterpreter,
 0, 1, 0, 0, 0x0100, "WB_RBLevels", &stdInterpreter,
 0, 1, 0, 0, 0x0102, "WB_RBLevels3000K", &stdInterpreter,
 0, 1, 0, 0, 0x0103, "WB_RBLevels3300K", &stdInterpreter,
 0, 1, 0, 0, 0x0104, "WB_RBLevels3600K", &stdInterpreter,
 0, 1, 0, 0, 0x0105, "WB_RBLevels3900K", &stdInterpreter,
 0, 1, 0, 0, 0x0106, "WB_RBLevels4000K", &stdInterpreter,
 0, 1, 0, 0, 0x0107, "WB_RBLevels4300K", &stdInterpreter,
 0, 1, 0, 0, 0x0108, "WB_RBLevels4500K", &stdInterpreter,
 0, 1, 0, 0, 0x0109, "WB_RBLevels4800K", &stdInterpreter,
 0, 1, 0, 0, 0x010a, "WB_RBLevels5300K", &stdInterpreter,
 0, 1, 0, 0, 0x010b, "WB_RBLevels6000K", &stdInterpreter,
 0, 1, 0, 0, 0x010c, "WB_RBLevels6600K", &stdInterpreter,
 0, 1, 0, 0, 0x010d, "WB_RBLevels7500K", &stdInterpreter,
 0, 1, 0, 0, 0x010e, "WB_RBLevelsCWB1", &stdInterpreter,
 0, 1, 0, 0, 0x010f, "WB_RBLevelsCWB2", &stdInterpreter,
 0, 1, 0, 0, 0x0110, "WB_RBLevelsCWB3", &stdInterpreter,
 0, 1, 0, 0, 0x0111, "WB_RBLevelsCWB4", &stdInterpreter,
 0, 1, 0, 0, 0x0113, "WB_GLevel3000K", &stdInterpreter,
 0, 1, 0, 0, 0x0114, "WB_GLevel3300K", &stdInterpreter,
 0, 1, 0, 0, 0x0115, "WB_GLevel3600K", &stdInterpreter,
 0, 1, 0, 0, 0x0116, "WB_GLevel3900K", &stdInterpreter,
 0, 1, 0, 0, 0x0117, "WB_GLevel4000K", &stdInterpreter,
 0, 1, 0, 0, 0x0118, "WB_GLevel4300K", &stdInterpreter,
 0, 1, 0, 0, 0x0119, "WB_GLevel4500K", &stdInterpreter,
 0, 1, 0, 0, 0x011a, "WB_GLevel4800K", &stdInterpreter,
 0, 1, 0, 0, 0x011b, "WB_GLevel5300K", &stdInterpreter,
 0, 1, 0, 0, 0x011c, "WB_GLevel6000K", &stdInterpreter,
 0, 1, 0, 0, 0x011d, "WB_GLevel6600K", &stdInterpreter,
 0, 1, 0, 0, 0x011e, "WB_GLevel7500K", &stdInterpreter,
 0, 1, 0, 0, 0x011f, "WB_GLevel", &stdInterpreter,
 0, 1, 0, 0, 0x0200, "ColorMatrix", &stdInterpreter,
 0, 1, 0, 0, 0x0300, "Enhancer", &stdInterpreter,
 0, 1, 0, 0, 0x0301, "EnhancerValues", &stdInterpreter,
 0, 1, 0, 0, 0x0310, "CoringFilter", &stdInterpreter,
 0, 1, 0, 0, 0x0311, "CoringValues", &stdInterpreter,
 0, 1, 0, 0, 0x0600, "BlackLevel2", &stdInterpreter,
 0, 1, 0, 0, 0x0610, "GainBase", &stdInterpreter,
 0, 1, 0, 0, 0x0611, "ValidBits", &stdInterpreter,
 0, 1, 0, 0, 0x0612, "CropLeft", &stdInterpreter,
 0, 1, 0, 0, 0x0613, "CropTop", &stdInterpreter,
 0, 1, 0, 0, 0x0614, "CropWidth", &stdInterpreter,
 0, 1, 0, 0, 0x0615, "CropHeight", &stdInterpreter,
 0, 1, 0, 0, 0x1010, "NoiseReduction2", &stdInterpreter,
 0, 1, 0, 0, 0x1011, "DistortionCorrection2", &olOnOffInterpreter,
 0, 1, 0, 0, 0x1012, "ShadingCompensation2", &olOnOffInterpreter,
 1, 1, 0, 0, 0x1103, "UnknownBlock", &stdInterpreter,
 0, 1, 0, 0, 0x1200, "FaceDetect", &olOnOffInterpreter,
 0, 1, 0, 0, 0x1201, "FaceDetectArea", &stdInterpreter,
-1, 0, 0,  0, 0, "", NULL};

const TagAttrib olyRawDevelopmentAttribs[] = {
 0, 1, 0, 0, 0x0000, "RawDevVersion", &stdInterpreter,
 0, 1, 0, 0, 0x0100, "RawDevExposureBiasValue", &stdInterpreter,
 0, 1, 0, 0, 0x0101, "RawDevWhiteBalanceValue", &stdInterpreter,
 0, 1, 0, 0, 0x0102, "RawDevWBFineAdjustment", &stdInterpreter,
 0, 1, 0, 0, 0x0103, "RawDevGrayPoint", &stdInterpreter,
 0, 1, 0, 0, 0x0104, "RawDevSaturationEmphasis", &stdInterpreter,
 0, 1, 0, 0, 0x0105, "RawDevMemoryColorEmphasis", &stdInterpreter,
 0, 1, 0, 0, 0x0106, "RawDevContrastValue", &stdInterpreter,
 0, 1, 0, 0, 0x0107, "RawDevSharpnessValue", &stdInterpreter,
 0, 1, 0, 0, 0x0108, "RawDevColorSpace", &olColorSpaceInterpreter,
 0, 1, 0, 0, 0x0109, "RawDevEngine", &olDevEngineInterpreter,
 0, 1, 0, 0, 0x010a, "RawDevNoiseReduction", &olNoiseReductionInterpreter,
 0, 1, 0, 0, 0x010b, "RawDevEditStatus", &stdInterpreter,
 0, 1, 0, 0, 0x010c, "RawDevSettings", &stdInterpreter,
-1, 0, 0,  0, 0, "", NULL};

const TagAttrib olyRawDevelopment2Attribs[] = {
 0, 1, 0, 0, 0x0000, "RawDevVersion", &stdInterpreter,
 0, 1, 0, 0, 0x0100, "RawDevExposureBiasValue", &stdInterpreter,
 0, 1, 0, 0, 0x0101, "RawDevWhiteBalance", &stdInterpreter,
 0, 1, 0, 0, 0x0102, "RawDevWhiteBalanceValue", &stdInterpreter,
 0, 1, 0, 0, 0x0103, "RawDevWBFineAdjustment", &stdInterpreter,
 0, 1, 0, 0, 0x0104, "RawDevGrayPoint", &stdInterpreter,
 0, 1, 0, 0, 0x0105, "RawDevContrastValue", &stdInterpreter,
 0, 1, 0, 0, 0x0106, "RawDevSharpnessValue", &stdInterpreter,
 0, 1, 0, 0, 0x0107, "RawDevSaturationEmphasis", &stdInterpreter,
 0, 1, 0, 0, 0x0108, "RawDevMemoryColorEmphasis", &stdInterpreter,
 0, 1, 0, 0, 0x0109, "RawDevColorSpace", &olColorSpaceInterpreter,
 0, 1, 0, 0, 0x010a, "RawDevNoiseReduction", &olNoiseReductionInterpreter,
 0, 1, 0, 0, 0x010b, "RawDevEngine", &olDevEngineInterpreter,
 0, 1, 0, 0, 0x010c, "RawDevPictureMode", &olPictureModeInterpreter,
 0, 1, 0, 0, 0x010d, "RawDevPMSaturation", &stdInterpreter,
 0, 1, 0, 0, 0x010e, "RawDevPMContrast", &stdInterpreter,
 0, 1, 0, 0, 0x010f, "RawDevPMSharpness", &stdInterpreter,
 0, 1, 0, 0, 0x0110, "RawDevPM_BWFilter", &olPictureModeBWFilterInterpreter,
 0, 1, 0, 0, 0x0111, "RawDevPMPictureTone", &olPictureModeToneInterpreter,
 0, 1, 0, 0, 0x0112, "RawDevGradation", &stdInterpreter,
 0, 1, 0, 0, 0x0113, "RawDevSaturation3", &stdInterpreter,
 0, 1, 0, 0, 0x0119, "RawDevAutoGradation", &olOnOffInterpreter,
 0, 1, 0, 0, 0x0120, "RawDevPMNoiseFilter", &stdInterpreter,
-1, 0, 0,  0, 0, "", NULL};

const TagAttrib olyCameraSettingsAttribs[] = {
 0, 1, 0, 0, 0x0000, "CameraSettingsVersion", &stdInterpreter,
 1, 1, 0, 0, 0x0100, "PreviewImageValid", &olYesNoInterpreter,
 1, 1, 0, 0, 0x0101, "PreviewImageStart", &stdInterpreter,
 1, 1, 0, 0, 0x0102, "PreviewImageLength", &stdInterpreter,
 0, 1, 0, 0, 0x0200, "ExposureMode", &olExposureModeInterpreter,
 0, 1, 0, 0, 0x0201, "AELock", &olOnOffInterpreter,
 0, 1, 0, 0, 0x0202, "MeteringMode", &olMeteringModeInterpreter,
 0, 1, 0, 0, 0x0300, "MacroMode", &olOnOffInterpreter,
 0, 1, 0, 0, 0x0301, "FocusMode", &olFocusModeInterpreter,
 0, 1, 0, 0, 0x0302, "FocusProcess", &stdInterpreter,
 0, 1, 0, 0, 0x0303, "AFSearch", &stdInterpreter,
 0, 1, 0, 0, 0x0304, "AFAreas", &stdInterpreter,
 0, 1, 0, 0, 0x0400, "FlashMode", &stdInterpreter,
 0, 1, 0, 0, 0x0401, "FlashExposureComp", &stdInterpreter,
 0, 1, 0, 0, 0x0500, "WhiteBalance2", &olWhitebalance2Interpreter,
 0, 1, 0, 0, 0x0501, "WhiteBalanceTemperature", &stdInterpreter,
 0, 1, 0, 0, 0x0502, "WhiteBalanceBracket", &stdInterpreter,
 0, 1, 0, 0, 0x0503, "CustomSaturation", &stdInterpreter,
 0, 1, 0, 0, 0x0504, "ModifiedSaturation", &stdInterpreter,
 0, 1, 0, 0, 0x0505, "ContrastSetting", &stdInterpreter,
 0, 1, 0, 0, 0x0506, "SharpnessSetting", &stdInterpreter,
 0, 1, 0, 0, 0x0507, "ColorSpace", &olColorSpaceInterpreter,
 0, 1, 0, 0, 0x0509, "SceneMode", &olSceneModeInterpreter,
 0, 1, 0, 0, 0x050a, "NoiseReduction", &olNoiseReductionInterpreter,
 0, 1, 0, 0, 0x050b, "DistortionCorrection", &olOnOffInterpreter,
 0, 1, 0, 0, 0x050c, "ShadingCompensation", &olOnOffInterpreter,
 0, 1, 0, 0, 0x050d, "CompressionFactor", &stdInterpreter,
 0, 1, 0, 0, 0x050f, "Gradation", &stdInterpreter,
 0, 1, 0, 0, 0x0520, "PictureMode", &olPictureModeInterpreter,
 0, 1, 0, 0, 0x0521, "PictureModeSaturation", &stdInterpreter,
 0, 1, 0, 0, 0x0522, "PictureModeHue", &stdInterpreter,
 0, 1, 0, 0, 0x0523, "PictureModeContrast", &stdInterpreter,
 0, 1, 0, 0, 0x0524, "PictureModeSharpness", &stdInterpreter,
 0, 1, 0, 0, 0x0525, "PictureModeBWFilter", &olPictureModeBWFilterInterpreter,
 0, 1, 0, 0, 0x0526, "PictureModeTone", &olPictureModeToneInterpreter,
 0, 1, 0, 0, 0x0527, "NoiseFilter", &olNoiseFilterInterpreter,
 0, 1, 0, 0, 0x0600, "DriveMode", &stdInterpreter,
 0, 1, 0, 0, 0x0601, "PanoramaMode", &stdInterpreter,
 0, 1, 0, 0, 0x0603, "ImageQuality2", &olImageQuality2Interpreter,
 0, 1, 0, 0, 0x0900, "ManometerPressure", &stdInterpreter,
 0, 1, 0, 0, 0x0901, "ManometerReading", &stdInterpreter,
 0, 1, 0, 0, 0x0902, "ExtendedWBDetect", &olOnOffInterpreter,
 -1, 0, 0,  0, 0, "", NULL};

const TagAttrib olyEquipmentAttribs[] = {
 0, 1, 0, 0, 0x0000, "EquipmentVersion", &stdInterpreter,
 0, 1, 0, 0, 0x0100, "CameraType2", &stdInterpreter,
 0, 1, 0, 0, 0x0101, "SerialNumber", &stdInterpreter,
 0, 1, 0, 0, 0x0102, "InternalSerialNumber", &stdInterpreter,
 0, 1, 0, 0, 0x0103, "FocalPlaneDiagonal", &stdInterpreter,
 0, 1, 0, 0, 0x0104, "BodyFirmwareVersion", &stdInterpreter,
 0, 1, 0, 0, 0x0201, "LensType", &olLensTypeInterpreter,
 0, 1, 0, 0, 0x0202, "LensSerialNumber", &stdInterpreter,
 0, 1, 0, 0, 0x0204, "LensFirmwareVersion", &stdInterpreter,
 0, 1, 0, 0, 0x0205, "MaxApertureAtMinFocal", &olApertureInterpreter,
 0, 1, 0, 0, 0x0206, "MaxApertureAtMaxFocal", &olApertureInterpreter,
 0, 1, 0, 0, 0x0207, "MinFocalLength", &stdInterpreter,
 0, 1, 0, 0, 0x0208, "MaxFocalLength", &stdInterpreter,
 0, 1, 0, 0, 0x020a, "MaxApertureAtCurrentFocal", &olApertureInterpreter,
 0, 1, 0, 0, 0x020b, "LensProperties", &stdInterpreter,
 0, 1, 0, 0, 0x0301, "Extender", &stdInterpreter,
 0, 1, 0, 0, 0x0302, "ExtenderSerialNumber", &stdInterpreter,
 0, 1, 0, 0, 0x0303, "ExtenderModel", &stdInterpreter,
 0, 1, 0, 0, 0x0304, "ExtenderFirmwareVersion", &stdInterpreter,
 0, 1, 0, 0, 0x1000, "FlashType", &olFlashTypeInterpreter,
 0, 1, 0, 0, 0x1001, "FlashModel", &olFlashModelInterpreter,
 0, 1, 0, 0, 0x1002, "FlashFirmwareVersion", &stdInterpreter,
 0, 1, 0, 0, 0x1003, "FlashSerialNumber", &stdInterpreter,
 -1, 0, 0,  0, 0, "", NULL};

const TagAttrib olympusAttribs[] = {
 0, 1, 0, 0, 0x0104, "BodyFirmwareVersion", &stdInterpreter,
 0, 1, 0, 0, 0x0200, "SpecialMode", &stdInterpreter,
 0, 1, 0, 0, 0x0201, "Quality", &stdInterpreter,
 0, 1, 0, 0, 0x0202, "Macro", &olOnOffInterpreter,
 0, 1, 0, 0, 0x0203, "BWMode", &olOnOffInterpreter,
 0, 1, 0, 0, 0x0204, "DigitalZoom", &stdInterpreter,
 0, 1, 0, 0, 0x0205, "FocalPlaneDiagonal", &stdInterpreter,
 0, 1, 0, 0, 0x0206, "LensDistortionParams", &stdInterpreter,
 0, 1, 0, 0, 0x0207, "CameraType", &stdInterpreter,
 1, 1, 0, 0, 0x0208, "TextInfo", &stdInterpreter,
 0, 1, 0, 0, 0x0209, "CameraID", &stdInterpreter,
 0, 1, 0, 0, 0x020b, "EpsonImageWidth", &stdInterpreter,
 0, 1, 0, 0, 0x020c, "EpsonImageHeight", &stdInterpreter,
 0, 1, 0, 0, 0x020d, "EpsonSoftware", &stdInterpreter,
 0, 2, 0, 0, 0x0280, "PreviewImage", &stdInterpreter,
 0, 1, 0, 0, 0x0300, "PreCaptureFrames", &stdInterpreter,
 0, 1, 0, 0, 0x0301, "WhiteBoard", &stdInterpreter,
 0, 1, 0, 0, 0x0302, "OneTouchWB", &olOnOffInterpreter,
 0, 1, 0, 0, 0x0303, "WhiteBalanceBracket", &stdInterpreter,
 0, 1, 0, 0, 0x0304, "WhiteBalanceBias", &stdInterpreter,
 0, 1, 0, 0, 0x0403, "SceneMode", &stdInterpreter,
 0, 1, 0, 0, 0x0404, "SerialNumber", &stdInterpreter,
 0, 1, 0, 0, 0x0405, "Firmware", &stdInterpreter,
 1, 1, 0, 0, 0x0e00, "PrintIM", &stdInterpreter,
 0, 1, 0, 0, 0x0f00, "DataDump", &stdInterpreter,
 0, 1, 0, 0, 0x0f01, "DataDump2", &stdInterpreter,
 0, 1, 0, 0, 0x1000, "ShutterSpeedValue", &stdInterpreter,
 0, 1, 0, 0, 0x1001, "ISOValue", &stdInterpreter,
 0, 1, 0, 0, 0x1002, "ApertureValue", &stdInterpreter,
 0, 1, 0, 0, 0x1003, "BrightnessValue", &stdInterpreter,
 0, 1, 0, 0, 0x1004, "FlashMode", &stdInterpreter,
 0, 1, 0, 0, 0x1005, "FlashDevice", &stdInterpreter,
 0, 1, 0, 0, 0x1006, "ExposureCompensation", &stdInterpreter,
 0, 1, 0, 0, 0x1007, "SensorTemperature", &stdInterpreter,
 0, 1, 0, 0, 0x1008, "LensTemperature", &stdInterpreter,
 0, 1, 0, 0, 0x1009, "LightCondition", &stdInterpreter,
 0, 1, 0, 0, 0x100a, "FocusRange", &stdInterpreter,
 0, 1, 0, 0, 0x100b, "FocusMode", &stdInterpreter,
 0, 1, 0, 0, 0x100c, "ManualFocusDistance", &stdInterpreter,
 0, 1, 0, 0, 0x100d, "ZoomStepCount", &stdInterpreter,
 0, 1, 0, 0, 0x100e, "FocusStepCount", &stdInterpreter,
 0, 1, 0, 0, 0x100f, "Sharpness", &stdInterpreter,
 0, 1, 0, 0, 0x1010, "FlashChargeLevel", &stdInterpreter,
 0, 1, 0, 0, 0x1011, "ColorMatrix", &stdInterpreter,
 0, 1, 0, 0, 0x1012, "BlackLevel", &stdInterpreter,
 0, 1, 0, 0, 0x1013, "ColorTemperatureBG", &stdInterpreter,
 0, 1, 0, 0, 0x1014, "ColorTemperatureRG", &stdInterpreter,
 0, 1, 0, 0, 0x1015, "WBMode", &stdInterpreter,
 0, 1, 0, 0, 0x1017, "RedBalance", &stdInterpreter,
 0, 1, 0, 0, 0x1018, "BlueBalance", &stdInterpreter,
 0, 1, 0, 0, 0x1019, "ColorMatrixNumber", &stdInterpreter,
 0, 1, 0, 0, 0x101a, "SerialNumber", &stdInterpreter,
 0, 1, 0, 0, 0x101b, "ExternalFlashAE1_0", &stdInterpreter,
 0, 1, 0, 0, 0x101c, "ExternalFlashAE2_0", &stdInterpreter,
 0, 1, 0, 0, 0x101d, "InternalFlashAE1_0", &stdInterpreter,
 0, 1, 0, 0, 0x101e, "InternalFlashAE2_0", &stdInterpreter,
 0, 1, 0, 0, 0x101f, "ExternalFlashAE1", &stdInterpreter,
 0, 1, 0, 0, 0x1020, "ExternalFlashAE2", &stdInterpreter,
 0, 1, 0, 0, 0x1021, "InternalFlashAE1", &stdInterpreter,
 0, 1, 0, 0, 0x1022, "InternalFlashAE2", &stdInterpreter,
 0, 1, 0, 0, 0x1023, "FlashExposureComp", &stdInterpreter,
 0, 1, 0, 0, 0x1024, "InternalFlashTable", &stdInterpreter,
 0, 1, 0, 0, 0x1025, "ExternalFlashGValue", &stdInterpreter,
 0, 1, 0, 0, 0x1026, "ExternalFlashBounce", &olYesNoInterpreter,
 0, 1, 0, 0, 0x1027, "ExternalFlashZoom", &stdInterpreter,
 0, 1, 0, 0, 0x1028, "ExternalFlashMode", &stdInterpreter,
 0, 1, 0, 0, 0x1029, "Contrast", &stdInterpreter,
 0, 1, 0, 0, 0x102a, "SharpnessFactor", &stdInterpreter,
 0, 1, 0, 0, 0x102b, "ColorControl", &stdInterpreter,
 0, 1, 0, 0, 0x102c, "ValidBits", &stdInterpreter,
 0, 1, 0, 0, 0x102d, "CoringFilter", &stdInterpreter,
 0, 1, 0, 0, 0x102e, "OlympusImageWidth", &stdInterpreter,
 0, 1, 0, 0, 0x102f, "OlympusImageHeight", &stdInterpreter,
 0, 1, 0, 0, 0x1030, "SceneDetect", &stdInterpreter,
 0, 1, 0, 0, 0x1031, "SceneArea", &stdInterpreter,
 0, 1, 0, 0, 0x1033, "SceneDetectData", &stdInterpreter,
 0, 1, 0, 0, 0x1034, "CompressionRatio", &stdInterpreter,
 1, 1, 0, 0, 0x1035, "PreviewImageValid", &olYesNoInterpreter,
 1, 1, 0, 0, 0x1036, "PreviewImageStart", &stdInterpreter,
 1, 1, 0, 0, 0x1037, "PreviewImageLength", &stdInterpreter,
 0, 1, 0, 0, 0x1038, "AFResult", &stdInterpreter,
 0, 1, 0, 0, 0x1039, "CCDScanMode", &stdInterpreter,
 0, 1, 0, 0, 0x103a, "NoiseReduction", &olOnOffInterpreter,
 0, 1, 0, 0, 0x103b, "InfinityLensStep", &stdInterpreter,
 0, 1, 0, 0, 0x103c, "NearLensStep", &stdInterpreter,
 0, 1, 0, 0, 0x103d, "LightValueCenter", &stdInterpreter,
 0, 1, 0, 0, 0x103e, "LightValuePeriphery", &stdInterpreter,
 0, 1, 0, 0, 0x103f, "FieldCount", &stdInterpreter,
 0, 1, 0, olyEquipmentAttribs, 0x2010, "Equipment", &stdInterpreter,
 0, 1, 0, olyCameraSettingsAttribs, 0x2020, "CameraSettings", &stdInterpreter,
 0, 1, 0, olyRawDevelopmentAttribs, 0x2030, "RawDevelopment", &stdInterpreter,
 0, 1, 0, olyRawDevelopment2Attribs, 0x2031, "RawDev2", &stdInterpreter,
 0, 1, 0, olyImageProcessingAttribs, 0x2040, "ImageProcessing", &stdInterpreter,
 0, 1, 0, olyFocusInfoAttribs, 0x2050, "FocusInfo", &stdInterpreter,
 1, 1, 0, 0, 0x2100, "Olympus2100", &stdInterpreter,
 1, 1, 0, 0, 0x2300, "Olympus2300", &stdInterpreter,
 1, 1, 0, 0, 0x2400, "Olympus2400", &stdInterpreter,
 1, 1, 0, 0, 0x2500, "Olympus2500", &stdInterpreter,
 1, 1, 0, 0, 0x2600, "Olympus2600", &stdInterpreter,
 1, 1, 0, 0, 0x2700, "Olympus2700", &stdInterpreter,
 1, 1, 0, 0, 0x2800, "Olympus2800", &stdInterpreter,
 1, 1, 0, 0, 0x2900, "Olympus2900", &stdInterpreter,
 0, 1, 0, 0, 0x3000, "RawInfo", &stdInterpreter,
 -1, 0, 0,  0, 0, "", NULL};
};
#endif

