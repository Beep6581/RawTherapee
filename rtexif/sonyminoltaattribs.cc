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

#include <cmath>

#include "rtexif.h"

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

class SAWhiteBalanceInterpreter: public ChoiceInterpreter {
	public:
	SAWhiteBalanceInterpreter(){
		choices[ 0x0] = "Auto";
		choices[ 0x1] = "Color Temperature/Color Filter";
		choices[0x10] = "Daylight";
		choices[0x20] = "Cloudy";
		choices[0x30] = "Shade";
		choices[0x40] = "Tungsten";
		choices[0x50] = "Flash";
		choices[0x60] = "Fluorescent";
		choices[0x70] = "Custom";
	}
};
SAWhiteBalanceInterpreter saWhiteBalanceInterpreter;

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

class SALensIDInterpreter : public IntLensInterpreter< int > {
    public:
        SALensIDInterpreter () {
            choices.insert(p_t(0, "Minolta AF 28-85mm f/3.5-4.5"));
            choices.insert(p_t(1, "Minolta AF 80-200mm f/2.8 HS-APO G"));
            choices.insert(p_t(2, "Minolta AF 28-70mm f/2.8 G"));
            choices.insert(p_t(3, "Minolta AF 28-80mm f/4-5.6"));
            choices.insert(p_t(5, "Minolta AF 35-70mm f/3.5-4.5"));
            choices.insert(p_t(6, "Minolta AF 24-85mm f/3.5-4.5 [New]"));
            choices.insert(p_t(7, "Minolta AF 100-300mm f/4.5-5.6 APO [New]"));
            choices.insert(p_t(7, "Sigma AF 100-300mm f/4 EX DG IF"));
            choices.insert(p_t(8, "Minolta AF 70-210mm f/4.5-5.6"));
            choices.insert(p_t(9, "Minolta AF 50mm f/3.5 Macro"));
            choices.insert(p_t(10, "Minolta AF 28-105mm f/3.5-4.5 [New]"));
            choices.insert(p_t(11, "Minolta AF 300mm f/4 HS-APO G"));
            choices.insert(p_t(12, "Minolta AF 100mm f/2.8 Soft Focus"));
            choices.insert(p_t(13, "Minolta AF 75-300mm f/4.5-5.6"));
            choices.insert(p_t(14, "Minolta AF 100-400mm f/4.5-6.7 APO"));
            choices.insert(p_t(15, "Minolta AF 400mm f/4.5 HS-APO G"));
            choices.insert(p_t(16, "Minolta AF 17-35mm f/3.5 G"));
            choices.insert(p_t(17, "Minolta AF 20-35mm f/3.5-4.5"));
            choices.insert(p_t(18, "Minolta AF 28-80mm f/3.5-5.6 II"));
            choices.insert(p_t(19, "Minolta AF 35mm f/1.4"));
            choices.insert(p_t(20, "Minolta/Sony STF 135mm F2.8 [T4.5]"));
            choices.insert(p_t(22, "Minolta AF 35-80mm f/4-5.6"));
            choices.insert(p_t(23, "Minolta AF 200mm f/4 G APO Macro"));
            choices.insert(p_t(24, "Minolta/Sony AF 24-105mm f/3.5-4.5 (D)"));
            choices.insert(p_t(24, "Sigma 18-50mm f/2.8 EX DC Macro"));
            choices.insert(p_t(24, "Sigma 17-70mm f/2.8-4.5 DC Macro"));
            choices.insert(p_t(24, "Sigma 20-40mm f/2.8 EX DG Aspherical IF"));
            choices.insert(p_t(24, "Sigma 18-200mm f/3.5-6.3 DC"));
            choices.insert(p_t(24, "Tamron SP AF 28-75mm f/2.8 XR Di (IF) Macro"));
            choices.insert(p_t(25, "Minolta AF 100-300mm f/4.5-5.6 APO D"));
            choices.insert(p_t(25, "Sigma 100-300mm f/4 EX DG APO"));
            choices.insert(p_t(25, "Sigma 70mm f/2.8 EX DG Macro"));
            choices.insert(p_t(25, "Sigma 20mm f/1.8 EX DG Aspherical RF"));
            choices.insert(p_t(25, "Sigma 30mm f/1.4 EX DG"));
            choices.insert(p_t(27, "Minolta AF 85mm f71.4 G"));
            choices.insert(p_t(28, "Minolta AF 100mm f/2.8 Macro (D)"));
            choices.insert(p_t(28, "Tamron SP AF 90mm f/2.8 Di Macro "));
            choices.insert(p_t(29, "Minolta AF 75-300mm f/4.5-5.6 (D)"));
            choices.insert(p_t(30, "Minolta AF 28-80mm f/3.5-5.6 (D)"));
            choices.insert(p_t(30, "Sigma 10-20mm f/4-5.6 EX DC"));
            choices.insert(p_t(30, "Sigma 12-24mm f/4.5-5.6 EX DG"));
            choices.insert(p_t(30, "Sigma 28-70mm f/2.8 EX DG"));
            choices.insert(p_t(30, "Sigma 55-200mm f/4-5.6 DC"));
            choices.insert(p_t(31, "Minolta/Sony AF 50mm f/2.8 Macro (D)"));
            choices.insert(p_t(32, "Minolta AF 300mm f/2.8 G"));
            choices.insert(p_t(33, "Minolta/Sony AF 70-200mm f/2.8 G (D) SSM"));
            choices.insert(p_t(35, "Minolta AF 85mm f/1.4 G (D) Limited"));
            choices.insert(p_t(36, "Minolta AF 28-100mm f/3.5-5.6 (D)"));
            choices.insert(p_t(38, "Minolta AF 17-35mm f/2.8-4 (D)"));
            choices.insert(p_t(39, "Minolta AF 28-75mm f/2.8 (D)"));
            choices.insert(p_t(40, "Minolta/Sony AF DT 18-70mm f/3.5-5.6 (D)"));
            choices.insert(p_t(41, "Minolta/Sony AF DT 11-18mm f/4.5-5.6 (D)"));
            choices.insert(p_t(42, "Minolta AF DT 18-200mm f/3.5-6.3 (D)"));
            choices.insert(p_t(43, "Minolta AF 35mm f/1.4 G"));
            choices.insert(p_t(44, "Sony AF 50mm f/1.4"));
            choices.insert(p_t(45, "Carl Zeiss Planar T* 85mm f/1.4 ZA"));
            choices.insert(p_t(46, "Carl Zeiss Vario-Sonnar T* DT 16-80mm f/3.5-4.5 ZA"));
            choices.insert(p_t(47, "Carl Zeiss Sonnar T* 135mm F1.8 ZA"));
            choices.insert(p_t(48, "Carl Zeiss Vario-Sonnar T* 24-70mm f/2.8 ZA SSM"));
            choices.insert(p_t(49, "Sony AF DT 55-200mm f/4-5.6"));
            choices.insert(p_t(50, "Sony AF DT 18-250mm f/3.5-6.3"));
            choices.insert(p_t(51, "Sony AF DT 16-105mm f/3.5-5.6 or 55-200mm f/4-5.5"));
            choices.insert(p_t(52, "Sony AF 70-300mm f/4.5-5.6 G SSM"));
            choices.insert(p_t(53, "Sony AF 70-400mm f/4.5-5.6 G SSM"));
            choices.insert(p_t(54, "Carl Zeiss Vario-Sonnar T* 16-35mm f/2.8 ZA SSM"));
            choices.insert(p_t(55, "Sony DT 18-55mm f/3.5-5.6 SAM"));
            choices.insert(p_t(56, "Sony AF DT 55-200mm f/4-5.6 SAM"));
            choices.insert(p_t(57, "Sony AF DT 50mm f/1.8 SAM"));
            choices.insert(p_t(58, "Sony AF DT 30mm f/2.8 SAM Macro"));
            choices.insert(p_t(59, "Sony AF 28-75mm f/2.8 SAM"));
            choices.insert(p_t(60, "Carl Zeiss Distagon T* 24mm f/2 ZA SSM"));
            choices.insert(p_t(61, "Sony AF 85mm f/2.8 SAM"));
            choices.insert(p_t(62, "Sony DT 35mm f/1.8 SAM"));
            choices.insert(p_t(128, "Tamron AF 18-200mm f/3.5-6.3 XR Di II LD Aspherical (IF)"));
            choices.insert(p_t(128, "Tamron AF 28-300mm f/3.5-6.3"));
            choices.insert(p_t(128, "Tamron AF 28-200mm f/3.8-5.6 XR Di Aspherical (IF) Macro "));
            choices.insert(p_t(128, "Tamron SP AF 17-35mm f/2.8-4 Di LD Aspherical IF"));
            choices.insert(p_t(128, "Sigma AF 50-150mm f/2.8 EX DC APO HSM II"));
            choices.insert(p_t(128, "Sigma 10-20mm f/3.5 EX DC"));
            choices.insert(p_t(128, "Sigma 70-200mm f/2.8 II EX DG APO Macro"));
            choices.insert(p_t(129, "Tamron 200-400mm f/5.6 LD (IF)"));
            choices.insert(p_t(129, "Tamron 70-300mm f/4-5.6 LD"));
            choices.insert(p_t(131, "Tamron 20-40mm f/2.7-3.5 SP Aspherical IF"));
            choices.insert(p_t(135, "Vivitar 28-210mm f/3.5-5.6"));
            choices.insert(p_t(136, "Tokina EMZ M100 AF 100mm f/3.5"));
            choices.insert(p_t(137, "Cosina 70-210mm f/2.8-4 AF"));
            choices.insert(p_t(138, "Soligor 19-35mm f/3.5-4.5"));
            choices.insert(p_t(142, "Voigtlander 70-300mm f/4.5-5.6"));
            choices.insert(p_t(146, "Voigtlander Macro APO-Lanthar 125mm f/2.5 SL"));
            choices.insert(p_t(255, "Tamron SP AF 17-50mm f/2.8 XR Di II LD Aspherical"));
            choices.insert(p_t(255, "Tamron AF 18-250mm f/3.5-6.3 XR Di II LD"));
            choices.insert(p_t(255, "Tamron AF 55-200mm f/4-5.6 Di II"));
            choices.insert(p_t(255, "Tamron AF 70-300mm f/4-5.6 Di LD Macro 1:2"));
            choices.insert(p_t(255, "Tamron SP AF 200-500mm f/5.0-6.3 Di LD (IF)"));
            choices.insert(p_t(255, "Tamron SP AF 10-24mm f/3.5-4.5 Di II LD Aspherical (IF)"));
            choices.insert(p_t(255, "Tamron SP AF 70-200mm f/2.8 Di LD Macro (IF)"));
            choices.insert(p_t(255, "Tamron SP AF 28-75mm f/2.8 XR Di LD Aspherical (IF)"));
            choices.insert(p_t(2550, "Minolta AF 50mm f/1.7"));
            choices.insert(p_t(2551, "Minolta AF 35-70mm f/4"));
            choices.insert(p_t(2551, "Sigma UC AF 28-70mm f/3.5-4.5"));
            choices.insert(p_t(2551, "Sigma AF 28-70mm f/2.8"));
            choices.insert(p_t(2551, "Sigma M-AF 70-200mm f/2.8 EX Aspherical"));
            choices.insert(p_t(2551, "Quantaray M-AF 35-80mm f/4-5.6"));
            choices.insert(p_t(2552, "Minolta AF 28-85mm f/3.5-4.5 [New]"));
            choices.insert(p_t(2552, "Tokina 19-35mm f/3.5-4.5"));
            choices.insert(p_t(2552, "Tokina 28-70mm f/2.8 AT-X"));
            choices.insert(p_t(2552, "Tokina 80-400mm f/4.5-5.6 AT-X AF II 840"));
            choices.insert(p_t(2552, "Tokina AF PRO 28-80mm f/2.8 AT-X 280"));
            choices.insert(p_t(2552, "Tokina AT-X PRO II AF 28-70mm f/2.6-2.8 270"));
            choices.insert(p_t(2552, "Tamron AF 19-35mm f/3.5-4.5"));
            choices.insert(p_t(2552, "Angenieux AF 28-70mm f/2.6"));
            choices.insert(p_t(2553, "Minolta AF 28-135mm f/4-4.5"));
            choices.insert(p_t(2553, "Sigma ZOOM-alpha 35-135mm f/3.5-4.5"));
            choices.insert(p_t(2553, "Sigma 28-105mm f/2.8-4 Aspherical"));
            choices.insert(p_t(2554, "Minolta AF 35-105mm f/3.5-4.5"));
            choices.insert(p_t(2555, "Minolta AF 70-210mm f/4 Macro"));
            choices.insert(p_t(2555, "Sigma 70-210mm f/4-5.6 APO"));
            choices.insert(p_t(2555, "Sigma M-AF 70-200mm f/2.8 EX APO"));
            choices.insert(p_t(2555, "Sigma 75-200mm f/2.8-3.5"));
            choices.insert(p_t(2556, "Minolta AF 135mm f/2.8"));
            choices.insert(p_t(2557, "Minolta AF 28mm f/2.8"));
            choices.insert(p_t(2558, "Minolta AF 24-50mm f/4"));
            choices.insert(p_t(2560, "Minolta AF 100-200mm f/4.5"));
            choices.insert(p_t(2561, "Minolta AF 75-300mm f/4.5-5.6"));
            choices.insert(p_t(2561, "Sigma 70-300mm f/4-5.6 DL Macro"));
            choices.insert(p_t(2561, "Sigma 300mm f/4 APO Macro"));
            choices.insert(p_t(2561, "Sigma AF 500mm f/4.5 APO"));
            choices.insert(p_t(2561, "Sigma AF 170-500mm f/5-6.3 APO Aspherical"));
            choices.insert(p_t(2561, "Tokina AT-X AF 300mm f/4"));
            choices.insert(p_t(2561, "Tokina AT-X AF 400mm f/5.6 SD"));
            choices.insert(p_t(2561, "Tokina AF 730 II 75-300mm f/4.5-5.6"));
            choices.insert(p_t(2562, "Minolta/Sony AF 50mm f/1.4 [New]"));
            choices.insert(p_t(2563, "Minolta AF 300mm f/2.8 G"));
            choices.insert(p_t(2563, "Sigma AF 50-500mm f/4-6.3 EX DG APO"));
            choices.insert(p_t(2563, "Sigma AF 170-500mm f/5-6.3 APO Aspherical"));
            choices.insert(p_t(2563, "Sigma AF 500mm f/4.5 EX DG APO"));
            choices.insert(p_t(2563, "Sigma 400mm f/5.6 APO"));
            choices.insert(p_t(2564, "Minolta AF 50mm f/2.8 Macro"));
            choices.insert(p_t(2564, "Sigma 50mm f/2.8 EX Macro"));
            choices.insert(p_t(2565, "Minolta AF 600mm f/4"));
            choices.insert(p_t(2566, "Minolta AF 24mm f/2.8"));
            choices.insert(p_t(2572, "Minolta/Sony AF 500mm f/8 Reflex"));
            choices.insert(p_t(2578, "Minolta AF 16mm f/2.8 Fisheye"));
            choices.insert(p_t(2578, "Sigma 8mm f/4 EX DG Fisheye"));
            choices.insert(p_t(2578, "Sigma 14mm f/3.5"));
            choices.insert(p_t(2578, "Sigma 15mm f/2.8 Fisheye"));
            choices.insert(p_t(2579, "Minolta AF 20mm f/2.8"));
            choices.insert(p_t(2581, "Minolta/Sony AF 100mm f/2.8 Macro"));
            choices.insert(p_t(2581, "Sigma AF 90mm f/2.8 Macro"));
            choices.insert(p_t(2581, "Sigma AF 105mm f/2.8 EX DG Macro"));
            choices.insert(p_t(2581, "Sigma 180mm f/5.6 Macro"));
            choices.insert(p_t(2581, "Tamron AF 90mm f/2.8 Macro"));
            choices.insert(p_t(2585, "Minolta AF 35-105mm f/3.5-4.5 New"));
            choices.insert(p_t(2585, "Beroflex 35-135mm f/3.5-4.5"));
            choices.insert(p_t(2585, "Tamron AF 24-135mm f/3.5-5.6"));
            choices.insert(p_t(2588, "Minolta AF 70-210mm f/3.5-4.5"));
            choices.insert(p_t(2589, "Minolta AF 80-200 f/2.8 APO"));
            choices.insert(p_t(2589, "Tokina 80-200mm f/2.8"));
            choices.insert(p_t(2591, "Minolta AF 35mm f/1.4"));
            choices.insert(p_t(2592, "Minolta AF 85mm f/1.4 G (D)"));
            choices.insert(p_t(2593, "Minolta AF 200mm f/2.8 G APO"));
            choices.insert(p_t(2594, "Minolta AF 3x-1x f/1.7-2.8 Macro"));
            choices.insert(p_t(2596, "Minolta AF 28mm f/2"));
            choices.insert(p_t(2597, "Minolta AF 35mm f/2"));
            choices.insert(p_t(2598, "Minolta AF 100mm f/2"));
            choices.insert(p_t(2604, "Minolta AF 80-200mm f/4.5-5.6"));
            choices.insert(p_t(2605, "Minolta AF 35-80mm f/4-5.6"));
            choices.insert(p_t(2606, "Minolta AF 100-300mm f/4.5-5.6 (D)"));
            choices.insert(p_t(2607, "Minolta AF 35-80mm f/4-5.6"));
            choices.insert(p_t(2608, "Minolta AF 300mm f/2.8 G"));
            choices.insert(p_t(2609, "Minolta AF 600mm f/4 HS-APO G"));
            choices.insert(p_t(2612, "Minolta AF 200mm f/2.8 G HS-APO"));
            choices.insert(p_t(2613, "Minolta AF 50mm f/1.7 New"));
            choices.insert(p_t(2615, "Minolta AF 28-105mm f/3.5-4.5 Power Zoom"));
            choices.insert(p_t(2616, "Minolta AF 35-200mm f/4.5-5.6 Power Zoom"));
            choices.insert(p_t(2618, "Minolta AF 28-80mm f/4-5.6 Power Zoom"));
            choices.insert(p_t(2619, "Minolta AF 80-200mm f/4.5-5.6 Power Zoom"));
            choices.insert(p_t(2620, "Minolta AF 28-70mm f/2.8 G"));
            choices.insert(p_t(2621, "Minolta AF 100-300mm f/4.5-5.6 Power Zoom"));
            choices.insert(p_t(2624, "Minolta AF 35-80mm f/4-5.6 Power Zoom"));
            choices.insert(p_t(2628, "Minolta AF 80-200mm f/2.8 G"));
            choices.insert(p_t(2629, "Minolta AF 85mm f/1.4 New"));
            choices.insert(p_t(2631, "Minolta/Sony AF 100-300mm f/4.5-5.6 APO"));
            choices.insert(p_t(2632, "Minolta AF 24-50mm f/4 New"));
            choices.insert(p_t(2638, "Minolta AF 50mm f/2.8 Macro New"));
            choices.insert(p_t(2639, "Minolta AF 100mm f/2.8 Macro"));
            choices.insert(p_t(2641, "Minolta AF 20mm f/2.8 New"));
            choices.insert(p_t(2642, "Minolta AF 24mm f/2.8 New"));
            choices.insert(p_t(2644, "Minolta AF 100-400mm f/4.5-6.7 APO"));
            choices.insert(p_t(2662, "Minolta AF 50mm f/1.4 New"));
            choices.insert(p_t(2667, "Minolta AF 35mm f/2 New"));
            choices.insert(p_t(2668, "Minolta AF 28mm f/2 New"));
            choices.insert(p_t(2672, "Minolta AF 24-105mm f/3.5-4.5 (D)"));
            choices.insert(p_t(4574, "Minolta AF 200mm f/2.8 G x2"));
            choices.insert(p_t(4575, "1.4 x Teleconverter"));
            choices.insert(p_t(4585, "Tamron SP AF 300mm f/2.8 LD IF"));
            choices.insert(p_t(6553, "Arax MC 35mm f/2.8 Tilt+Shift"));
            choices.insert(p_t(6553, "Arax MC 80mm f/2.8 Tilt+Shift"));
            choices.insert(p_t(6553, "Zenitar MF 16mm f/2.8 Fisheye M42"));
            choices.insert(p_t(6553, "Samyang 500mm Mirror f/8"));
            choices.insert(p_t(6553, "Pentacon Auto 135mm f/2.8"));
            choices.insert(p_t(6553, "Pentacon Auto 29mm f/2.8"));
            choices.insert(p_t(6553, "Helios 44-2 58mm f/2"));
            choices.insert(p_t(25501, "Minolta AF 50mm f/1.7"));
            choices.insert(p_t(25511, "Minolta AF 35-70mm f/4"));
            choices.insert(p_t(25511, "Sigma UC AF 28-70mm f/3.5-4.5"));
            choices.insert(p_t(25511, "Sigma AF 28-70mm f/2.8"));
            choices.insert(p_t(25511, "Sigma M-AF 70-200mm f/2.8 EX Aspherical"));
            choices.insert(p_t(25511, "Quantaray M-AF 35-80mm f/4-5.6"));
            choices.insert(p_t(25521, "Minolta AF 28-85mm f/3.5-4.5 [New]"));
            choices.insert(p_t(25521, "Tokina 19-35mm f/3.5-4.5"));
            choices.insert(p_t(25521, "Tokina 28-70mm f/2.8 AT-X"));
            choices.insert(p_t(25521, "Tokina 80-400mm f/4.5-5.6 AT-X AF II 840"));
            choices.insert(p_t(25521, "Tokina AF PRO 28-80mm f/2.8 AT-X 280"));
            choices.insert(p_t(25521, "Tokina AT-X PRO II AF 28-70mm f/2.6-2.8 270"));
            choices.insert(p_t(25521, "Tamron AF 19-35mm f/3.5-4.5"));
            choices.insert(p_t(25521, "Angenieux AF 28-70mm f/2.6"));
            choices.insert(p_t(25531, "Minolta AF 28-135mm f/4-4.5"));
            choices.insert(p_t(25531, "Sigma ZOOM-alpha 35-135mm f/3.5-4.5"));
            choices.insert(p_t(25531, "Sigma 28-105mm f/2.8-4 Aspherical"));
            choices.insert(p_t(25541, "Minolta AF 35-105mm f/3.5-4.5"));
            choices.insert(p_t(25551, "Minolta AF 70-210mm f/4 Macro"));
            choices.insert(p_t(25551, "Sigma 70-210mm f/4-5.6 APO"));
            choices.insert(p_t(25551, "Sigma M-AF 70-200mm f/2.8 EX APO"));
            choices.insert(p_t(25551, "Sigma 75-200mm f/2.8-3.5"));
            choices.insert(p_t(25561, "Minolta AF 135mm f/2.8"));
            choices.insert(p_t(25571, "Minolta AF 28mm f/2.8"));
            choices.insert(p_t(25581, "Minolta AF 24-50mm f/4"));
            choices.insert(p_t(25601, "Minolta AF 100-200mm f/4.5"));
            choices.insert(p_t(25611, "Minolta AF 75-300mm f/4.5-5.6"));
            choices.insert(p_t(25611, "Sigma 70-300mm f/4-5.6 DL Macro"));
            choices.insert(p_t(25611, "Sigma 300mm f/4 APO Macro"));
            choices.insert(p_t(25611, "Sigma AF 500mm f/4.5 APO"));
            choices.insert(p_t(25611, "Sigma AF 170-500mm f/5-6.3 APO Aspherical"));
            choices.insert(p_t(25611, "Tokina AT-X AF 300mm f/4"));
            choices.insert(p_t(25611, "Tokina AT-X AF 400mm f/5.6 SD"));
            choices.insert(p_t(25611, "Tokina AF 730 II 75-300mm f/4.5-5.6"));
            choices.insert(p_t(25621, "Minolta AF 50mm f/1.4"));
            choices.insert(p_t(25631, "Minolta AF 300mm f/2.8 G"));
            choices.insert(p_t(25631, "Sigma AF 50-500mm f/4-6.3 EX DG APO"));
            choices.insert(p_t(25631, "Sigma AF 170-500mm f/5-6.3 APO Aspherical"));
            choices.insert(p_t(25631, "Sigma AF 500mm f/4.5 EX DG APO"));
            choices.insert(p_t(25631, "Sigma 400mm f/5.6 APO"));
            choices.insert(p_t(25641, "Minolta AF 50mm f/2.8 Macro"));
            choices.insert(p_t(25641, "Sigma AF 50mm f/2.8 Macro"));
            choices.insert(p_t(25651, "Minolta AF 600mm f/4"));
            choices.insert(p_t(25661, "Minolta AF 24mm f/2.8"));
            choices.insert(p_t(25721, "Minolta/Sony AF 500mm f/8 Reflex"));
            choices.insert(p_t(25781, "Minolta AF 16mm f/2.8 Fisheye"));
            choices.insert(p_t(25781, "Sigma 8mm f/4 EX DG Fisheye"));
            choices.insert(p_t(25781, "Sigma 14mm f/3.5"));
            choices.insert(p_t(25781, "Sigma 15mm f/2.8 Fisheye"));
            choices.insert(p_t(25791, "Minolta AF 20mm f/2.8"));
            choices.insert(p_t(25811, "Minolta/Sony AF 100mm f/2.8 Macro New"));
            choices.insert(p_t(25811, "Sigma AF 90mm f/2.8 Macro"));
            choices.insert(p_t(25811, "Sigma AF 105mm f/2.8 EX DG Macro"));
            choices.insert(p_t(25811, "Sigma 180mm f/5.6 Macro"));
            choices.insert(p_t(25811, "Tamron 90mm f/2.8 Macro"));
            choices.insert(p_t(25851, "Beroflex 35-135mm f/3.5-4.5"));
            choices.insert(p_t(25858, "Minolta AF 35-105mm f/3.5-4.5 New"));
            choices.insert(p_t(25858, "Tamron 24-135mm f/3.5-5.6"));
            choices.insert(p_t(25881, "Minolta AF 70-210mm f/3.5-4.5"));
            choices.insert(p_t(25891, "Minolta AF 80-200 f/2.8 APO"));
            choices.insert(p_t(25891, "Tokina 80-200mm f/2.8"));
            choices.insert(p_t(25911, "Minolta AF 35mm f/1.4"));
            choices.insert(p_t(25921, "Minolta AF 85mm f/1.4 G (D)"));
            choices.insert(p_t(25931, "Minolta AF 200mm f/2.8 G APO"));
            choices.insert(p_t(25941, "Minolta AF 3x-1x f/1.7-2.8 Macro"));
            choices.insert(p_t(25961, "Minolta AF 28mm f/2"));
            choices.insert(p_t(25971, "Minolta AF 35mm f/2"));
            choices.insert(p_t(25981, "Minolta AF 100mm f/2"));
            choices.insert(p_t(26041, "Minolta AF 80-200mm f/4.5-5.6"));
            choices.insert(p_t(26051, "Minolta AF 35-80mm f/4-5.6"));
            choices.insert(p_t(26061, "Minolta AF 100-300mm f/4.5-5.6 (D)"));
            choices.insert(p_t(26071, "Minolta AF 35-80mm f/4-5.6"));
            choices.insert(p_t(26081, "Minolta AF 300mm f/2.8 HS-APO G"));
            choices.insert(p_t(26091, "Minolta AF 600mm f/4 HS-APO G"));
            choices.insert(p_t(26121, "Minolta AF 200mm f/2.8 HS-APO G"));
            choices.insert(p_t(26131, "Minolta AF 50mm f/1.7 New"));
            choices.insert(p_t(26151, "Minolta AF 28-105mm f/3.5-4.5 Power Zoom"));
            choices.insert(p_t(26161, "Minolta AF 35-200mm f/4.5-5.6 Power Zoom"));
            choices.insert(p_t(26181, "Minolta AF 28-80mm f/4-5.6 Power Zoom"));
            choices.insert(p_t(26191, "Minolta AF 80-200mm f/4.5-5.6 Power Zoom"));
            choices.insert(p_t(26201, "Minolta AF 28-70mm f/2.8 G"));
            choices.insert(p_t(26211, "Minolta AF 100-300mm f/4.5-5.6 Power Zoom"));
            choices.insert(p_t(26241, "Minolta AF 35-80mm f/4-5.6 Power Zoom"));
            choices.insert(p_t(26281, "Minolta AF 80-200mm f/2.8 G"));
            choices.insert(p_t(26291, "Minolta AF 85mm f/1.4 New"));
            choices.insert(p_t(26311, "Minolta/Sony AF 100-300mm f/4.5-5.6 APO"));
            choices.insert(p_t(26321, "Minolta AF 24-50mm f/4 New"));
            choices.insert(p_t(26381, "Minolta AF 50mm f/2.8 Macro New"));
            choices.insert(p_t(26391, "Minolta AF 100mm f/2.8 Macro"));
            choices.insert(p_t(26411, "Minolta AF 20mm f/2.8 New"));
            choices.insert(p_t(26421, "Minolta AF 24mm f/2.8 New"));
            choices.insert(p_t(26441, "Minolta AF 100-400mm f/4.5-6.7 APO"));
            choices.insert(p_t(26621, "Minolta AF 50mm f/1.4 New"));
            choices.insert(p_t(26671, "Minolta AF 35mm f/2 New"));
            choices.insert(p_t(26681, "Minolta AF 28mm f/2 New"));
            choices.insert(p_t(26721, "Minolta AF 24-105mm f/3.5-4.5 (D)"));
            choices.insert(p_t(45671, "Tokina 70-210mm f/4-5.6"));
            choices.insert(p_t(45741, "Minolta AF 200mm f/2.8 G x2"));
            choices.insert(p_t(45851, "Tamron SP AF 300mm f/2.8 LD IF"));
            choices.insert(p_t(45871, "Tamron SP AF 70-210mm f/2.8 LD"));
            choices.insert(p_t(65535, "Arax MC 35mm f/2.8 Tilt+Shift"));
            choices.insert(p_t(65535, "Arax MC 80mm f/2.8 Tilt+Shift"));
            choices.insert(p_t(65535, "Zenitar MF 16mm f/2.8 Fisheye M42"));
            choices.insert(p_t(65535, "Samyang 500mm f/8 Mirror"));
            choices.insert(p_t(65535, "Pentacon Auto 135mm f/2.8"));
            choices.insert(p_t(65535, "Pentacon Auto 29mm f/2.8"));
            choices.insert(p_t(65535, "Helios 44-2 58mm f/2"));
        }

        virtual std::string toString (Tag* t)
        {
        	 int lensID = t->toInt();
        	 Tag *apertureTag = t->getParent()->getRoot()->findTag("MaxApertureValue");
        	 Tag *focalLengthTag = t->getParent()->getRoot()->findTag("FocalLength");
        	 double maxApertureAtFocal = 0.;
        	 double focalLength = 0.;
        	 if( apertureTag )
        		 maxApertureAtFocal = pow(2.0, apertureTag->toDouble()/2.0);
        	 if( focalLengthTag )
        		 focalLength = focalLengthTag->toDouble();
        	 return guess( lensID, focalLength, maxApertureAtFocal );
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

class SAQualityInterpreter2 : public ChoiceInterpreter {
    public:
	SAQualityInterpreter2 () {
            choices[0]  = "Raw";
            choices[2]  = "cRAW";
            choices[16] = "Extra fine";
            choices[32]  = "Fine";
            choices[34]  = "RAW + JPEG";
            choices[35]  = "cRAW + JPEG";
            choices[48]  = "Standard";
        }
};
SAQualityInterpreter2 saQualityInterpreter2;

class SADriveMode : public ChoiceInterpreter {
    public:
	SADriveMode () {
            choices[0]  = "Single Frame";
            choices[1]  = "Continuous High";
            choices[4]  = "Self-timer 10 sec";
            choices[5]  = "Self-timer 2 sec";
            choices[7]  = "Continuous Bracketing";
            choices[12] = "Continuous Low";
            choices[18] = "White Balance Bracketing Low";
            choices[19] = "D-Range Optimizer Bracketing Low";
        }
};
SADriveMode saDriveMode;

class SAFocusMode: public ChoiceInterpreter {
	public:
	SAFocusMode () {
        choices[0]  = "Manual";
        choices[1]  = "AF-S";
        choices[2]  = "AF-C";
        choices[3]  = "AF-A";
        choices[4]  = "Permanent-AF";
        choices[65535] = "n/a";
    }
};
SAFocusMode saFocusMode;

class SAAFMode: public ChoiceInterpreter {
	public:
	SAAFMode(){
		choices[0] = "Default";
		choices[1] = "Multi AF";
		choices[2] = "Center AF";
		choices[3] = "Spot AF";
		choices[4] = "Flexible Spot AF";
		choices[6] = "Touch AF";
		choices[14] = "Manual Focus";
		choices[15] = "Face Detected";
		choices[65535] = "n/a";
	}
};
SAAFMode saAFMode;

class SAAFAreaMode: public ChoiceInterpreter {
	public:
	SAAFAreaMode () {
        choices[0]  = "Wide";
        choices[1]  = "Local";
        choices[2]  = "Spot";
    }
};
SAAFAreaMode saAFAreaMode;

class SALocalAFAreaPoint: public ChoiceInterpreter {
	public:
	SALocalAFAreaPoint () {
		choices[1] = "Center";
		choices[2] = "Top";
		choices[3] = "Top-Right";
		choices[4] = "Right";
		choices[5] = "Bottom-Right";
		choices[6] = "Bottom";
		choices[7] = "Bottom-Left";
		choices[8] = "Left";
		choices[9] = "Top-Left";
		choices[10] = "Far Right";
		choices[11] = "Far Left";
	}
};
SALocalAFAreaPoint saLocalAFAreaPoint;

class SAMeteringMode: public ChoiceInterpreter {
	public:
	SAMeteringMode () {
		choices[1] = "Multi-segment";
		choices[2] = "Center-weighted Average";
		choices[4] = "Spot";
	}
};
SAMeteringMode saMeteringMode;

class SADynamicRangeOptimizerMode: public ChoiceInterpreter {
	public:
	SADynamicRangeOptimizerMode () {
		choices[0] = "Off";
		choices[1] = "Standard";
		choices[2] = "Advanced Auto";
		choices[3] = "Advanced Level";
		choices[4097] = "Auto";
	}
};
SADynamicRangeOptimizerMode saDynamicRangeOptimizerMode;

class SACreativeStyle: public ChoiceInterpreter {
	public:
	SACreativeStyle () {
		choices[1] = "Standard";
		choices[2] = "Vivid";
		choices[3] = "Portrait";
		choices[4] = "Landscape";
		choices[5] = "Sunset";
		choices[6] = "Night View/Portrait";
		choices[8] = "B&W";
		choices[9] = "Adobe RGB";
		choices[11] = "Neutral";
		choices[12] = "Clear";
		choices[13] = "Deep";
		choices[14] = "Light";
		choices[15] = "Autumn";
		choices[16] = "Sepia";
	}
};
SACreativeStyle saCreativeStyle;

class SAFlashMode: public ChoiceInterpreter {
	public:
	SAFlashMode () {
		choices[0] = "ADI";
		choices[1] = "TTL";
	}
};
SAFlashMode saFlashMode;

class SAExposureProgram: public ChoiceInterpreter {
	public:
		SAExposureProgram () {
			choices[0] = "Auto";
			choices[1] = "Manual";
			choices[2] = "Program AE";
			choices[3] = "Aperture-priority AE";
			choices[4] = "Shutter speed priority AE";
			choices[8] = "Program Shift A";
			choices[9] = "Program Shift S";
			choices[16] = "Portrait";
			choices[17] = "Sports";
			choices[18] = "Sunset";
			choices[19] = "Night Portrait";
			choices[20] = "Landscape";
			choices[21] = "Macro";
			choices[35] = "Auto No Flash";
	}
};
SAExposureProgram saExposureProgram;

class SARotation: public ChoiceInterpreter {
	public:
		SARotation () {
			choices[0] = "Horizontal";
			choices[1] = "Rotate 90 CW";
			choices[2] = "Rotate 270 CW";
			choices[3] = "None";
	}
};
SARotation saRotation;

class SASonyImageSize: public ChoiceInterpreter {
	public:
		SASonyImageSize () {
			choices[1] = "Large";
			choices[2] = "Medium ";
			choices[3] = "Small";
	}
};
SASonyImageSize saSonyImageSize;

class SAAspectRatio: public ChoiceInterpreter {
	public:
		SAAspectRatio () {
			choices[1] = "3:2";
			choices[2] = "16:9";
	}
};
SAAspectRatio saAspectRatio;

class SAExposureLevelIncrements: public ChoiceInterpreter {
	public:
		SAExposureLevelIncrements () {
			choices[33] = "1/3 EV";
			choices[50] = "1/2 EV";
	}
};
SAExposureLevelIncrements saExposureLevelIncrements;

class SAAFIlluminator: public ChoiceInterpreter {
	public:
	SAAFIlluminator () {
			choices[0] = "Off";
			choices[1] = "Auto";
			choices[65535]="n/a";
	}
};
SAAFIlluminator saAFIlluminator;

class SAReleaseModeInterpreter: public ChoiceInterpreter {
	public:
	SAReleaseModeInterpreter () {
		choices[0] = "Normal";
		choices[2] = "Burst";
		choices[5] = "Exposure Bracketing";
		choices[6] = "White Balance Bracketing";
		choices[65535]="n/a";
	}
};
SAReleaseModeInterpreter saReleaseModeInterpreter;

class SAImageStyleInterpreter: public ChoiceInterpreter {
	public:
	SAImageStyleInterpreter () {
		choices[1] = "Standard";
		choices[2] = "Vivid";
		choices[9] = "Adobe RGB";
		choices[11] = "Neutral";
		choices[129]="StyleBox1";
		choices[130]="StyleBox2";
		choices[131]="StyleBox3";
	}
};
SAImageStyleInterpreter saImageStyleInterpreter;

const TagAttrib minoltaAttribs[] = {
 {0, 1, 0, 0, 0x0000, "MakerNoteVersion", &stdInterpreter},
 {0, 1, 0, 0, 0x0001, "MinoltaCameraSettingsOld", &stdInterpreter},
 {0, 1, 0, 0, 0x0003, "MinoltaCameraSettings", &stdInterpreter},
 {0, 1, 0, 0, 0x0004, "MinoltaCameraSettings7D", &stdInterpreter},
 {0, 1, 0, 0, 0x0018, "ImageStabilization", &stdInterpreter},
 {0, 1, 0, 0, 0x0040, "CompressedImageSize", &stdInterpreter},
 {1, 1, 0, 0, 0x0081, "PreviewImage", &stdInterpreter},
 {1, 1, 0, 0, 0x0088, "PreviewImageStart", &stdInterpreter},
 {1, 1, 0, 0, 0x0089, "PreviewImageLength", &stdInterpreter},
 {0, 1, 0, 0, 0x0100, "SceneMode", &saSceneModeInterpreter},
 {0, 1, 0, 0, 0x0101, "ColorMode", &saColorModeInterpreter},
 {0, 1, 0, 0, 0x0102, "MinoltaQuality", &maQualityInterpreter},
 {0, 1, 0, 0, 0x0103, "MinoltaImageSize", &maImageSizeInterpreter},
 {0, 1, 0, 0, 0x0104, "FlashExposureComp", &stdInterpreter},
 {0, 1, 0, 0, 0x0105, "Teleconverter", &maTeleconverterInterpreter},
 {0, 1, 0, 0, 0x0107, "ImageStabilization", &saOnOffInterpreter},
 {0, 1, 0, 0, 0x010a, "ZoneMatching", &saZoneMatchingInterpreter},
 {0, 1, 0, 0, 0x010b, "ColorTemperature", &stdInterpreter},
 {0, 1, 0, 0, 0x010c, "LensID", &saLensIDInterpreter},
 {0, 1, 0, 0, 0x0113, "ImageStabilization", &saOnOffInterpreter},
 {0, 1, 0, 0, 0x0114, "MinoltaCameraSettings", &stdInterpreter},
 {1, 1, 0, 0, 0x0e00, "PrintIM", &stdInterpreter},
 {0, 1, 0, 0, 0x0f00, "MinoltaCameraSettings2", &stdInterpreter},
 {-1, 0, 0,  0, 0, "", NULL}};

const TagAttrib sonyAttribs[] = {
 {0, 1, 0, 0, 0x0102, "Quality", &maQualityInterpreter},
 {0, 1, 0, 0, 0x0104, "FlashExposureComp",&stdInterpreter},
 {0, 1, 0, 0, 0x0106, "TeleConverter", &maTeleconverterInterpreter},
 {0, 1, 0, sonyCameraSettingsAttribs, 0x0114, "SonyCameraSettings",&stdInterpreter},
 {0, 1, 0, 0, 0x0115, "WhiteBalance",&saWhiteBalanceInterpreter},
 {1, 1, 0, 0, 0x0e00, "PrintIM", &stdInterpreter},
 {1, 1, 0, 0, 0x2001, "PreviewImage", &stdInterpreter},
 {0, 1, 0, 0, 0x200a, "AutoHDR", &stdInterpreter},
 {0, 1, 0, 0, 0xb020, "ColorReproduction", &stdInterpreter},
 {0, 1, 0, 0, 0xb021, "ColorTemperature", &stdInterpreter},
 {0, 1, 0, 0, 0xb023, "SceneMode", &saSceneModeInterpreter},
 {0, 1, 0, 0, 0xb024, "ZoneMatching", &saZoneMatchingInterpreter},
 {0, 1, 0, 0, 0xb025, "DynamicRangeOptimizer", &saDynamicRangeOptimizerInterpreter},
 {0, 1, 0, 0, 0xb026, "ImageStabilization", &saOnOffInterpreter},
 {0, 1, 0, 0, 0xb027, "LensID", &saLensIDInterpreter},
 {0, 1, 0, minoltaAttribs, 0xb028, "MinoltaMakerNote", &stdInterpreter},
 {0, 1, 0, 0, 0xb029, "ColorMode", &saColorModeInterpreter},
 {0, 1, 0, 0, 0xb040, "Macro", &saOnOffInterpreter},
 {0, 1, 0, 0, 0xb041, "ExposureMode", &saExposureModeInterpreter},
 {0, 1, 0, 0, 0xb042, "FocusMode", &saFocusMode},
 {0, 1, 0, 0, 0xb043, "AFMode", &saAFMode},
 {0, 1, 0, 0, 0xb044, "AFIlluminator", &saAFIlluminator},
 {0, 1, 0, 0, 0xb047, "Quality", &saQualityInterpreter},
 {0, 1, 0, 0, 0xb049, "ReleaseMode",&saReleaseModeInterpreter},
 {0, 1, 0, 0, 0xb04b, "AntiBlur", &saAntiBlurInterpreter},
 {0, 1, 0, 0, 0xb04e, "LongExposureNoiseReduction", &saOnOffInterpreter},
 {-1, 0, 0,  0, 0, "", NULL}};

const TagAttrib sonyCameraSettingsAttribs[]={
 {0, 1, 0, 0, 4, "DriveMode", &saDriveMode},
 {0, 1, 0, 0, 6, "WhiteBalanceFineTune",&stdInterpreter},
 {0, 1, 0, 0, 16, "FocusMode",&saFocusMode},
 {0, 1, 0, 0, 17, "AFAreaMode",&saAFAreaMode},
 {0, 1, 0, 0, 18, "LocalAFAreaPoint", &saLocalAFAreaPoint},
 {0, 1, 0, 0, 21, "MeteringMode",&saMeteringMode},
 {0, 1, 0, 0, 22, "ISOSetting",&stdInterpreter},
 {0, 1, 0, 0, 24, "DynamicRangeOptimizerMode", &saDynamicRangeOptimizerMode},
 {0, 1, 0, 0, 25, "DynamicRangeOptimizerLevel",&stdInterpreter},
 {0, 1, 0, 0, 26, "CreativeStyle",&saCreativeStyle},
 {0, 1, 0, 0, 28, "Sharpness",&stdInterpreter},
 {0, 1, 0, 0, 29, "Contrast",&stdInterpreter},
 {0, 1, 0, 0, 30, "Saturation",&stdInterpreter},
 {0, 1, 0, 0, 31, "ZoneMatchingValue",&stdInterpreter},
 {0, 1, 0, 0, 34, "Brightness",&stdInterpreter},
 {0, 1, 0, 0, 35, "FlashMode",&saFlashMode},
 {0, 1, 0, 0, 40, "PrioritySetupShutterRelease",&stdInterpreter},
 {0, 1, 0, 0, 41, "AFIlluminator",&saAFIlluminator},
 {0, 1, 0, 0, 42, "AFWithShutter",&saOnOffInterpreter},
 {0, 1, 0, 0, 43, "LongExposureNoiseReduction",&saOnOffInterpreter},
 {0, 1, 0, 0, 44, "HighISONoiseReduction",&stdInterpreter},
 {0, 1, 0, 0, 45, "ImageStyle",&saImageStyleInterpreter},
 {0, 1, 0, 0, 60, "ExposureProgram",&saExposureProgram},
 {0, 1, 0, 0, 61, "ImageStabilization",&saOnOffInterpreter},
 {0, 1, 0, 0, 63, "Rotation",&saRotation},
 {0, 1, 0, 0, 84, "SonyImageSize",&saSonyImageSize},
 {0, 1, 0, 0, 85, "AspectRatio",&saAspectRatio},
 {0, 1, 0, 0, 86, "Quality",&saQualityInterpreter2},
 {0, 1, 0, 0, 88, "ExposureLevelIncrements",&saExposureLevelIncrements},
 {-1, 0, 0,  0, 0, "", NULL}};

const TagAttrib sonyCameraSettingsAttribs2[]={
 {0, 1, 0, 0, 16, "FocusMode",&saFocusMode},
 {0, 1, 0, 0, 17, "AFAreaMode",&saAFAreaMode},
 {0, 1, 0, 0, 18, "LocalAFAreaPoint",&saLocalAFAreaPoint},
 {0, 1, 0, 0, 19, "MeteringMode",&saMeteringMode},
 {0, 1, 0, 0, 20, "ISOSetting",&stdInterpreter},
 {0, 1, 0, 0, 22, "DynamicRangeOptimizerMode",&saDynamicRangeOptimizerMode},
 {0, 1, 0, 0, 23, "DynamicRangeOptimizerLevel",&stdInterpreter},
 {0, 1, 0, 0, 24, "CreativeStyle",&saCreativeStyle},
 {0, 1, 0, 0, 25, "Sharpness",&stdInterpreter},
 {0, 1, 0, 0, 26, "Contrast",&stdInterpreter},
 {0, 1, 0, 0, 27, "Saturation",&stdInterpreter},
 {0, 1, 0, 0, 35, "FlashMode",&saFlashMode},
 {0, 1, 0, 0, 60, "ExposureProgram",&saExposureProgram},
 {0, 1, 0, 0, 63, "Rotation",&saRotation},
 {0, 1, 0, 0, 84, "SonyImageSize",&saSonyImageSize},
 {-1, 0, 0,  0, 0, "", NULL}};
}
#endif


