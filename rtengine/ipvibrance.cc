/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Copyright (c) 2011 Jacques Desmis  <jdesmis@gmail.com>
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

#include "rt_math.h"
#include "rtengine.h"
#include "improcfun.h"
#include "iccstore.h"
#include "mytime.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace rtengine {

using namespace procparams;

#define SAT(a,b,c) ((float)max(a,b,c)-(float)min(a,b,c))/(float)max(a,b,c)

extern const Settings* settings;

//Munsell Lch LUTf : 195 LUT
LUTf ImProcFunctions::_4P10 ;//give hue in function of L and C : Munsell  correction
LUTf ImProcFunctions::_4P20 ;
LUTf ImProcFunctions::_4P30 ;
LUTf ImProcFunctions::_4P40 ;
LUTf ImProcFunctions::_4P50 ;
LUTf ImProcFunctions::_4P60 ;


LUTf ImProcFunctions::_1P10 ;
LUTf ImProcFunctions::_1P20 ;
LUTf ImProcFunctions::_1P30 ;
LUTf ImProcFunctions::_1P40 ;
LUTf ImProcFunctions::_1P50 ;
LUTf ImProcFunctions::_1P60 ;


LUTf ImProcFunctions::_10PB10 ;
LUTf ImProcFunctions::_10PB20 ;
LUTf ImProcFunctions::_10PB30 ;
LUTf ImProcFunctions::_10PB40 ;
LUTf ImProcFunctions::_10PB50 ;
LUTf ImProcFunctions::_10PB60 ;


LUTf ImProcFunctions::_9PB10 ;
LUTf ImProcFunctions::_9PB20 ;
LUTf ImProcFunctions::_9PB30 ;
LUTf ImProcFunctions::_9PB40 ;
LUTf ImProcFunctions::_9PB50 ;
LUTf ImProcFunctions::_9PB60 ;
LUTf ImProcFunctions::_9PB70 ;
LUTf ImProcFunctions::_9PB80 ;

LUTf ImProcFunctions::_75PB10 ;
LUTf ImProcFunctions::_75PB20 ;
LUTf ImProcFunctions::_75PB30 ;
LUTf ImProcFunctions::_75PB40 ;
LUTf ImProcFunctions::_75PB50 ;
LUTf ImProcFunctions::_75PB60 ;
LUTf ImProcFunctions::_75PB70 ;
LUTf ImProcFunctions::_75PB80 ;

LUTf ImProcFunctions::_6PB10 ;
LUTf ImProcFunctions::_6PB20 ;
LUTf ImProcFunctions::_6PB30 ;
LUTf ImProcFunctions::_6PB40 ;
LUTf ImProcFunctions::_6PB50 ;
LUTf ImProcFunctions::_6PB60 ;
LUTf ImProcFunctions::_6PB70 ;
LUTf ImProcFunctions::_6PB80 ;

LUTf ImProcFunctions::_45PB10 ;
LUTf ImProcFunctions::_45PB20 ;
LUTf ImProcFunctions::_45PB30 ;
LUTf ImProcFunctions::_45PB40 ;
LUTf ImProcFunctions::_45PB50 ;
LUTf ImProcFunctions::_45PB60 ;
LUTf ImProcFunctions::_45PB70 ;
LUTf ImProcFunctions::_45PB80 ;

LUTf ImProcFunctions::_3PB10 ;
LUTf ImProcFunctions::_3PB20 ;
LUTf ImProcFunctions::_3PB30 ;
LUTf ImProcFunctions::_3PB40 ;
LUTf ImProcFunctions::_3PB50 ;
LUTf ImProcFunctions::_3PB60 ;
LUTf ImProcFunctions::_3PB70 ;
LUTf ImProcFunctions::_3PB80 ;

LUTf ImProcFunctions::_15PB10 ;
LUTf ImProcFunctions::_15PB20 ;
LUTf ImProcFunctions::_15PB30 ;
LUTf ImProcFunctions::_15PB40 ;
LUTf ImProcFunctions::_15PB50 ;
LUTf ImProcFunctions::_15PB60 ;
LUTf ImProcFunctions::_15PB70 ;
LUTf ImProcFunctions::_15PB80 ;

LUTf ImProcFunctions::_05PB40 ;
LUTf ImProcFunctions::_05PB50 ;
LUTf ImProcFunctions::_05PB60 ;
LUTf ImProcFunctions::_05PB70 ;
LUTf ImProcFunctions::_05PB80 ;

LUTf ImProcFunctions::_10B40 ;
LUTf ImProcFunctions::_10B50 ;
LUTf ImProcFunctions::_10B60 ;
LUTf ImProcFunctions::_10B70 ;
LUTf ImProcFunctions::_10B80 ;

LUTf ImProcFunctions::_9B40 ;
LUTf ImProcFunctions::_9B50 ;
LUTf ImProcFunctions::_9B60 ;
LUTf ImProcFunctions::_9B70 ;
LUTf ImProcFunctions::_9B80 ;

LUTf ImProcFunctions::_7B40 ;
LUTf ImProcFunctions::_7B50 ;
LUTf ImProcFunctions::_7B60 ;
LUTf ImProcFunctions::_7B70 ;
LUTf ImProcFunctions::_7B80 ;

LUTf ImProcFunctions::_5B40 ;
LUTf ImProcFunctions::_5B50 ;
LUTf ImProcFunctions::_5B60 ;
LUTf ImProcFunctions::_5B70 ;
LUTf ImProcFunctions::_5B80 ;

LUTf ImProcFunctions::_10YR20;
LUTf ImProcFunctions::_10YR30;
LUTf ImProcFunctions::_10YR40;
LUTf ImProcFunctions::_10YR50;
LUTf ImProcFunctions::_10YR60;
LUTf ImProcFunctions::_10YR70;
LUTf ImProcFunctions::_10YR80;
LUTf ImProcFunctions::_10YR90;

LUTf ImProcFunctions::_85YR20;
LUTf ImProcFunctions::_85YR30;
LUTf ImProcFunctions::_85YR40;
LUTf ImProcFunctions::_85YR50;
LUTf ImProcFunctions::_85YR60;
LUTf ImProcFunctions::_85YR70;
LUTf ImProcFunctions::_85YR80;
LUTf ImProcFunctions::_85YR90;

LUTf ImProcFunctions::_7YR30;
LUTf ImProcFunctions::_7YR40;
LUTf ImProcFunctions::_7YR50;
LUTf ImProcFunctions::_7YR60;
LUTf ImProcFunctions::_7YR70;
LUTf ImProcFunctions::_7YR80;

LUTf ImProcFunctions::_55YR30;
LUTf ImProcFunctions::_55YR40;
LUTf ImProcFunctions::_55YR50;
LUTf ImProcFunctions::_55YR60;
LUTf ImProcFunctions::_55YR70;
LUTf ImProcFunctions::_55YR80;
LUTf ImProcFunctions::_55YR90;

LUTf ImProcFunctions::_4YR30;
LUTf ImProcFunctions::_4YR40;
LUTf ImProcFunctions::_4YR50;
LUTf ImProcFunctions::_4YR60;
LUTf ImProcFunctions::_4YR70;
LUTf ImProcFunctions::_4YR80;

LUTf ImProcFunctions::_25YR30;
LUTf ImProcFunctions::_25YR40;
LUTf ImProcFunctions::_25YR50;
LUTf ImProcFunctions::_25YR60;
LUTf ImProcFunctions::_25YR70;

LUTf ImProcFunctions::_10R30;
LUTf ImProcFunctions::_10R40;
LUTf ImProcFunctions::_10R50;
LUTf ImProcFunctions::_10R60;
LUTf ImProcFunctions::_10R70;

LUTf ImProcFunctions::_9R30;
LUTf ImProcFunctions::_9R40;
LUTf ImProcFunctions::_9R50;
LUTf ImProcFunctions::_9R60;
LUTf ImProcFunctions::_9R70;

LUTf ImProcFunctions::_7R30;
LUTf ImProcFunctions::_7R40;
LUTf ImProcFunctions::_7R50;
LUTf ImProcFunctions::_7R60;
LUTf ImProcFunctions::_7R70;

LUTf ImProcFunctions::_5R10;
LUTf ImProcFunctions::_5R20;
LUTf ImProcFunctions::_5R30;

LUTf ImProcFunctions::_25R10;
LUTf ImProcFunctions::_25R20;
LUTf ImProcFunctions::_25R30;

LUTf ImProcFunctions::_10RP10;
LUTf ImProcFunctions::_10RP20;
LUTf ImProcFunctions::_10RP30;

LUTf ImProcFunctions::_7G30;
LUTf ImProcFunctions::_7G40;
LUTf ImProcFunctions::_7G50;
LUTf ImProcFunctions::_7G60;
LUTf ImProcFunctions::_7G70;
LUTf ImProcFunctions::_7G80;

LUTf ImProcFunctions::_5G30;
LUTf ImProcFunctions::_5G40;
LUTf ImProcFunctions::_5G50;
LUTf ImProcFunctions::_5G60;
LUTf ImProcFunctions::_5G70;
LUTf ImProcFunctions::_5G80;

LUTf ImProcFunctions::_25G30;
LUTf ImProcFunctions::_25G40;
LUTf ImProcFunctions::_25G50;
LUTf ImProcFunctions::_25G60;
LUTf ImProcFunctions::_25G70;
LUTf ImProcFunctions::_25G80;

LUTf ImProcFunctions::_1G30;
LUTf ImProcFunctions::_1G40;
LUTf ImProcFunctions::_1G50;
LUTf ImProcFunctions::_1G60;
LUTf ImProcFunctions::_1G70;
LUTf ImProcFunctions::_1G80;

LUTf ImProcFunctions::_10GY30;
LUTf ImProcFunctions::_10GY40;
LUTf ImProcFunctions::_10GY50;
LUTf ImProcFunctions::_10GY60;
LUTf ImProcFunctions::_10GY70;
LUTf ImProcFunctions::_10GY80;

LUTf ImProcFunctions::_75GY30;
LUTf ImProcFunctions::_75GY40;
LUTf ImProcFunctions::_75GY50;
LUTf ImProcFunctions::_75GY60;
LUTf ImProcFunctions::_75GY70;
LUTf ImProcFunctions::_75GY80;

LUTf ImProcFunctions::_5GY30;
LUTf ImProcFunctions::_5GY40;
LUTf ImProcFunctions::_5GY50;
LUTf ImProcFunctions::_5GY60;
LUTf ImProcFunctions::_5GY70;
LUTf ImProcFunctions::_5GY80;

/*
 * Munsell Lch correction
 * copyright (c) 2011  Jacques Desmis <jdesmis@gmail.com>
 *
 * data (Munsell ==> Lab) obtained with WallKillcolor and http://www.cis.rit.edu/research/mcsl2/online/munsell.php
 * each LUT give Hue in function of C, for each color Munsell and Luminance
 * eg: _6PB20 : color Munsell 6PB for L=20 c=5 c=45 c=85 c=125..139 when possible: interpolation betwwen values
 * no value for C<5  (gray)
 * low memory usage -- maximum: 195 LUT * 140 values
 * errors due to small number of point of LUT and linearization are very low (1 to 2%)
 * errors due to a different illuminant "Daylight" than "C" are low in the order of 10%. For example, a theoretical correction of 0.1 radian will be made with a real correction of 0.09 or 0.11 depending on the color illuminant D50
 * errors due to the use of a very different illuminant "C", for example illuminant "A" (tungsten) are higher in the order of 20%. Theoretical correction of 0.52 radians will be made with a real correction of 0.42
 */
void ImProcFunctions::initMunsell () {
#ifdef _DEBUG
	MyTime t1e,t2e;
	t1e.set();
#endif

int maxInd  = 140;
int maxInd2 = 90;
int maxInd3 = 50;

//blue for sky
_5B40(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _5B40[i] = -2.3 + 0.0025*(i-5);
		else if (i<90 && i>=45) _5B40[i] = -2.2 + 0.00*(i-45);
	}
	//printf("5B %1.2f  %1.2f\n",_5B40[44],_5B40[89]);
_5B50(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _5B50[i] = -2.34 + 0.0025*(i-5);
		else if (i<90 && i>=45) _5B50[i] = -2.24+0.0003*(i-45);
	}
	//printf("5B %1.2f  %1.2f\n",_5B50[44],_5B50[89]);
_5B60(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _5B60[i] = -2.4 + 0.003*(i-5);
		else if (i<90 && i>=45) _5B60[i] = -2.28+0.0005*(i-45);
	}
	//printf("5B %1.2f  %1.2f\n",_5B60[44],_5B60[89]);
_5B70(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _5B70[i] = -2.41 + 0.00275*(i-5);
		else if (i<90 && i>=45) _5B70[i] = -2.30+0.00025*(i-45);
	}
	//printf("5B %1.2f  %1.2f\n",_5B70[44],_5B70[89]);
_5B80(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _5B80[i] = -2.45 +0.003*(i-5);
	}
	//printf("5B %1.2f\n",_5B80[49]);

_7B40(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _7B40[i] = -2.15 + 0.0027*(i-5);
		else if (i<90 && i>=45) _7B40[i] = -2.04 + 0.00*(i-45);
	}
	//printf("7B %1.2f  %1.2f\n",_7B40[44],_7B40[89]);
_7B50(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _7B50[i] = -2.20 + 0.003*(i-5);
		else if (i<90 && i>=45) _7B50[i] = -2.08 + 0.001*(i-45);
	}
	//printf("7B %1.2f  %1.2f\n",_7B50[44],_7B50[79]);
_7B60(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _7B60[i] = -2.26 + 0.0035*(i-5);
		else if (i<90 && i>=45) _7B60[i] = -2.12 + 0.001*(i-45);
	}
	//printf("7B %1.2f  %1.2f\n",_7B60[44],_7B60[79]);
_7B70(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _7B70[i] = -2.28 + 0.003*(i-5);
		else if (i<90 && i>=45) _7B70[i] = -2.16 + 0.0015*(i-45);
	}
	//printf("7B %1.2f  %1.2f\n",_7B70[44],_7B70[64]);
_7B80(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _7B80[i] = -2.30 +0.0028*(i-5);
	}
	//printf("5B %1.2f\n",_7B80[49]);

_9B40(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _9B40[i] = -1.99 + 0.0022*(i-5);
		else if (i<90 && i>=45) _9B40[i] = -1.90 + 0.0008*(i-45);
	}
	//printf("9B %1.2f  %1.2f\n",_9B40[44],_9B40[69]);
_9B50(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _9B50[i] = -2.04 + 0.0025*(i-5);
		else if (i<90 && i>=45) _9B50[i] = -1.94 + 0.0013*(i-45);
	}
	//printf("9B %1.2f  %1.2f\n",_9B50[44],_9B50[77]);
_9B60(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _9B60[i] = -2.10 + 0.0033*(i-5);
		else if (i<90 && i>=45) _9B60[i] = -1.97 + 0.001*(i-45);
	}
	//printf("9B %1.2f  %1.2f\n",_9B60[44],_9B60[79]);
_9B70(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _9B70[i] = -2.12 + 0.003*(i-5);
		else if (i<90 && i>=45) _9B70[i] = -2.00 + 0.001*(i-45);
	}
	//printf("9B %1.2f  %1.2f\n",_9B70[44],_9B70[54]);
_9B80(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _9B80[i] = -2.16 +0.0025*(i-5);
	}
	//printf("9B %1.2f\n",_9B80[49]);

_10B40(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _10B40[i] = -1.92 + 0.0022*(i-5);
		else if (i<90 && i>=45) _10B40[i] = -1.83 + 0.0012*(i-45);
	}
	//printf("10B %1.2f  %1.2f\n",_10B40[44],_10B40[76]);
_10B50(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _10B50[i] = -1.95 + 0.0022*(i-5);
		else if (i<90 && i>=45) _10B50[i] = -1.86 + 0.0008*(i-45);
	}
	//printf("10B %1.2f  %1.2f\n",_10B50[44],_10B50[85]);
_10B60(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _10B60[i] = -2.01 + 0.0027*(i-5);
		else if (i<90 && i>=45) _10B60[i] = -1.90 + 0.0012*(i-45);
	}
	//printf("10B %1.2f  %1.2f\n",_10B60[44],_10B60[70]);
_10B70(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _10B70[i] = -2.03 +0.0025*(i-5);
	}
	//printf("10B %1.2f\n",_10B70[49]);
_10B80(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _10B80[i] = -2.08 +0.0032*(i-5);
	}
	//printf("10B %1.2f\n",_10B80[39]);

_05PB40(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _05PB40[i] = -1.87 + 0.0022*(i-5);
		else if (i<90 && i>=45) _05PB40[i] = -1.78 + 0.0015*(i-45);
	}
	//printf("05PB %1.2f  %1.2f\n",_05PB40[44],_05PB40[74]);
_05PB50(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _05PB50[i] = -1.91 + 0.0022*(i-5);
		else if (i<90 && i>=45) _05PB50[i] = -1.82 + 0.001*(i-45);
	}
	//printf("05PB %1.2f  %1.2f\n",_05PB50[44],_05PB50[85]);
_05PB60(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _05PB60[i] = -1.96 + 0.0027*(i-5);
		else if (i<90 && i>=45) _05PB60[i] = -1.85 + 0.0013*(i-45);
	}
	//printf("05PB %1.2f  %1.2f\n",_05PB60[44],_05PB60[70]);
_05PB70(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _05PB70[i] = -1.99 + 0.0027*(i-5);
		else if (i<90 && i>=45) _05PB70[i] = -1.88 + 0.001*(i-45);
	}
	//printf("05PB %1.2f  %1.2f\n",_05PB70[44],_05PB70[54]);
_05PB80(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _05PB80[i] = -2.03 +0.003*(i-5);
	}
	//printf("05PB %1.2f\n",_05PB80[39]);



//blue purple correction
//between 15PB to 4P
//maximum deviation 75PB

//15PB
_15PB10(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _15PB10[i] = -1.66 +0.0035*(i-5);
	}
	//printf("15 %1.2f\n",_15PB10[49]);
_15PB20(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _15PB20[i] = -1.71 +0.00275*(i-5);
		else if (i<90 && i>=45) _15PB20[i] = -1.60+0.0012*(i-45);
	}
	//printf("15 %1.2f  %1.2f\n",_15PB20[44],_15PB20[89]);

_15PB30(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _15PB30[i] = -1.75 +0.0025*(i-5);
		else if (i<90 && i>=45) _15PB30[i] = -1.65+0.002*(i-45);
	}
	//printf("15 %1.2f  %1.2f\n",_15PB30[44],_15PB30[89]);

_15PB40(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _15PB40[i] = -1.79 +0.002*(i-5);
		else if (i<90 && i>=45) _15PB40[i] = -1.71+0.002*(i-45);
	}
	//printf("15 %1.2f  %1.2f\n",_15PB40[44],_15PB40[89]);

_15PB50(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _15PB50[i] = -1.82 +0.002*(i-5);
		else if (i<90 && i>=45) _15PB50[i] = -1.74+0.0011*(i-45);
	}
	//printf("15 %1.2f  %1.2f\n",_15PB50[44],_15PB50[89]);

_15PB60(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _15PB60[i] = -1.87 +0.0025*(i-5);
		else if (i<90 && i>=45) _15PB60[i] = -1.77+0.001*(i-45);
	}
	//printf("15 %1.2f  %1.2f\n",_15PB60[44],_15PB60[89]);
_15PB70(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _15PB70[i] = -1.90 +0.0027*(i-5);
	}
	//	printf("15 %1.2f\n",_15PB70[49]);
_15PB80(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _15PB80[i] = -1.93 +0.0027*(i-5);
	}
	//printf("15 %1.2f %1.2f\n",_15PB80[38], _15PB80[49]);

//3PB
_3PB10(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _3PB10[i] = -1.56 +0.005*(i-5);
		else if (i<90 && i>=45) _3PB10[i] = -1.36+0.001*(i-45);
	}
	//printf("30 %1.2f  %1.2f\n",_3PB10[44],_3PB10[89]);

_3PB20(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _3PB20[i] = -1.59 +0.00275*(i-5);
		else if (i<90 && i>=45) _3PB20[i] = -1.48+0.003*(i-45);
	}
	//printf("30 %1.2f  %1.2f\n",_3PB20[44],_3PB20[89]);

_3PB30(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _3PB30[i] = -1.62 +0.00225*(i-5);
		else if (i<90 && i>=45) _3PB30[i] = -1.53+0.0032*(i-45);
	}
	//printf("30 %1.2f  %1.2f\n",_3PB30[44],_3PB30[89]);

_3PB40(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _3PB40[i] = -1.64 +0.0015*(i-5);
		else if (i<90 && i>=45) _3PB40[i] = -1.58+0.0025*(i-45);
	}
	//printf("30 %1.2f  %1.2f\n",_3PB40[44],_3PB40[89]);

_3PB50(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _3PB50[i] = -1.69 +0.00175*(i-5);
		else if (i<90 && i>=45) _3PB50[i] = -1.62+0.002*(i-45);
	}
	//printf("30 %1.2f  %1.2f\n",_3PB50[44],_3PB50[89]);

_3PB60(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _3PB60[i] = -1.73 +0.002*(i-5);
		else if (i<90 && i>=45) _3PB60[i] = -1.65+0.0012*(i-45);
	}
	//printf("30 %1.2f  %1.2f\n",_3PB60[44],_3PB60[89]);
_3PB70(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _3PB70[i] = -1.76 +0.002*(i-5);
	}
	//printf("30 %1.2f\n",_3PB70[49]);
_3PB80(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _3PB80[i] = -1.78 +0.0025*(i-5);
	}
	//printf("30 %1.2f %1.2f\n",_3PB80[38], _3PB80[49]);

//45PB
_45PB10(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _45PB10[i] = -1.46 +0.0045*(i-5);
		else if (i<90 && i>=45) _45PB10[i] = -1.28+0.0025*(i-45);
	}
	//printf("45 %1.2f  %1.2f\n",_45PB10[44],_45PB10[89]);

_45PB20(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _45PB20[i] = -1.48 +0.00275*(i-5);
		else if (i<90 && i>=45) _45PB20[i] = -1.37+0.0025*(i-45);
	}
	//printf("45 %1.2f  %1.2f\n",_45PB20[44],_45PB20[89]);

_45PB30(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _45PB30[i] = -1.51 +0.00175*(i-5);
		else if (i<90 && i>=45) _45PB30[i] = -1.44+0.0035*(i-45);
	}
	//printf("45 %1.2f  %1.2f\n",_45PB30[44],_45PB30[89]);

_45PB40(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _45PB40[i] = -1.52 +0.001*(i-5);
		else if (i<90 && i>=45) _45PB40[i] = -1.48+0.003*(i-45);
	}
	//printf("45 %1.2f  %1.2f\n",_45PB40[44],_45PB40[89]);

_45PB50(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _45PB50[i] = -1.55 +0.001*(i-5);
		else if (i<90 && i>=45) _45PB50[i] = -1.51+0.0022*(i-45);
	}
	//printf("45 %1.2f  %1.2f\n",_45PB50[44],_45PB50[89]);

_45PB60(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _45PB60[i] = -1.6 +0.0015*(i-5);
		else if (i<90 && i>=45) _45PB60[i] = -1.54+0.001*(i-45);
	}
	//printf("45 %1.2f  %1.2f\n",_45PB60[44],_45PB60[89]);
_45PB70(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _45PB70[i] = -1.63 +0.0017*(i-5);
	}
	//printf("45 %1.2f\n",_45PB70[49]);
_45PB80(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _45PB80[i] = -1.67 +0.0025*(i-5);
	}
	//printf("45 %1.2f %1.2f\n",_45PB80[38], _45PB80[49]);

//_6PB
_6PB10(maxInd);
	for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
		if (i<45 && i>5) _6PB10[i] = -1.33 +0.005*(i-5);
		else if (i<85 && i>=45) _6PB10[i] = -1.13+0.0045*(i-45);
		else if (i<140 && i >=85) _6PB10[i] = -0.95+0.0015*(i-85);
	}
	//printf("60 %1.2f  %1.2f %1.2f\n",_6PB10[44],_6PB10[84],_6PB10[139]);

_6PB20(maxInd);
	for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
		if (i<45 && i>5) _6PB20[i] = -1.36 +0.004*(i-5);
		else if (i<85 && i>=45) _6PB20[i] = -1.20+0.00375*(i-45);
		else if (i<140 && i >=85) _6PB20[i] = -1.05+0.0017*(i-85);
	}
	//printf("60 %1.2f  %1.2f %1.2f\n",_6PB20[44],_6PB20[84],_6PB20[139]);

_6PB30(maxInd);
	for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
		if (i<45 && i>5) _6PB30[i] = -1.38 +0.00225*(i-5);
		else if (i<85 && i>=45) _6PB30[i] = -1.29+0.00375*(i-45);
		else if (i<140 && i >=85) _6PB30[i] = -1.14+0.002*(i-85);
	}
	//printf("60 %1.2f  %1.2f %1.2f\n",_6PB30[44],_6PB30[84],_6PB30[139]);

_6PB40(maxInd);
	for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
		if (i<45 && i>5) _6PB40[i] = -1.39 +0.00125*(i-5);
		else if (i<85 && i>=45) _6PB40[i] = -1.34+0.00275*(i-45);
		else if (i<140 && i >=85) _6PB40[i] = -1.23+0.002*(i-85);
	}
	//printf("60 %1.2f  %1.2f %1.2f\n",_6PB40[44],_6PB40[84],_6PB40[139]);

_6PB50(maxInd2);//limits  -1.3   -1.11
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _6PB50[i] = -1.43 +0.00125*(i-5);
		else if (i<90 && i>=45) _6PB50[i] = -1.38+0.00225*(i-45);
	}
	//printf("60 %1.2f  %1.2f \n",_6PB50[44],_6PB50[89]);

_6PB60(maxInd2);//limits  -1.3   -1.11
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _6PB60[i] = -1.46 +0.0012*(i-5);
		else if (i<90 && i>=45) _6PB60[i] = -1.40+0.000875*(i-45);
	}
	//printf("60 %1.2f  %1.2f\n",_6PB60[44],_6PB60[89]);
_6PB70(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _6PB70[i] = -1.49 +0.0018*(i-5);
	}
	//printf("6 %1.2f\n",_6PB70[49]);
_6PB80(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _6PB80[i] = -1.52 +0.0022*(i-5);
	}
	//printf("6 %1.2f %1.2f\n",_6PB80[38], _6PB80[49]);


//_75PB : notation Munsell for maximum deviation blue purple
_75PB10(maxInd);//limits hue -1.23  -0.71  _75PBx   x=Luminance  eg_75PB10 for L >5 and L<=15
	for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
		if (i<45 && i>5) _75PB10[i] = -1.23 +0.0065*(i-5);
		else if (i<85 && i>=45) _75PB10[i] = -0.97+0.00375*(i-45);
		else if (i<140 && i >=85) _75PB10[i] = -0.82+0.0015*(i-85);
	}
	//printf("75 %1.2f  %1.2f %1.2f\n",_75PB10[44],_75PB10[84],_75PB10[139]);

_75PB20(maxInd);//limits -1.24  -0.79  for L>15 <=25
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _75PB20[i] = -1.24 +0.004*(i-5);
		else if (i<85 && i>=45) _75PB20[i] = -1.08+0.00425*(i-45);
		else if (i<140 && i >=85) _75PB20[i] = -0.91+0.0017*(i-85);
	}
	//printf("75 %1.2f  %1.2f %1.2f\n",_75PB20[44],_75PB20[84],_75PB20[139]);

_75PB30(maxInd);//limits -1.25  -0.85
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _75PB30[i] = -1.25 +0.00275*(i-5);
		else if (i<85 && i>=45) _75PB30[i] = -1.14+0.004*(i-45);
		else if (i<140 && i >=85) _75PB30[i] = -0.98+0.0015*(i-85);
	}
	//printf("75 %1.2f  %1.2f %1.2f\n",_75PB30[44],_75PB30[84],_75PB30[139]);

_75PB40(maxInd);//limits  -1.27  -0.92
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _75PB40[i] = -1.27 +0.002*(i-5);
		else if (i<85 && i>=45) _75PB40[i] = -1.19+0.003*(i-45);
		else if (i<140 && i >=85) _75PB40[i] = -1.07+0.0022*(i-85);
	}
	//printf("75 %1.2f  %1.2f %1.2f\n",_75PB40[44],_75PB40[84],_75PB40[139]);

_75PB50(maxInd2);//limits  -1.3   -1.11
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _75PB50[i] = -1.3 +0.00175*(i-5);
		else if (i<90 && i>=45) _75PB50[i] = -1.23+0.0025*(i-45);
	}
	//printf("75 %1.2f  %1.2f\n",_75PB50[44],_75PB50[89]);

_75PB60(maxInd2);
	for (int i=0; i<maxInd2; i++) {//limits -1.32  -1.17
		if (i<45 && i>5) _75PB60[i] = -1.32 +0.0015*(i-5);
		else if (i<90 && i>=45) _75PB60[i] = -1.26+0.002*(i-45);
	}
	//printf("75 %1.2f  %1.2f \n",_75PB60[44],_75PB60[89]);

_75PB70(maxInd3);
	for (int i=0; i<maxInd3; i++) {//limits  -1.34  -1.27
		if (i<50 && i>5) _75PB70[i] = -1.34 +0.002*(i-5);
	}
_75PB80(maxInd3);
	for (int i=0; i<maxInd3; i++) {//limits -1.35  -1.29
		if (i<50 && i>5) _75PB80[i] = -1.35 +0.00125*(i-5);
	}


_9PB10(maxInd);
	for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
		if (i<45 && i>5) _9PB10[i] = -1.09 +0.00475*(i-5);
		else if (i<85 && i>=45) _9PB10[i] = -0.9+0.003*(i-45);
		else if (i<140 && i >=85) _9PB10[i] = -0.78+0.0013*(i-85);
	}
	//printf("90 %1.2f  %1.2f %1.2f\n",_9PB10[44],_9PB10[84],_9PB10[139]);

_9PB20(maxInd);
	for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
		if (i<45 && i>5) _9PB20[i] = -1.12 +0.0035*(i-5);
		else if (i<85 && i>=45) _9PB20[i] = -0.98+0.00325*(i-45);
		else if (i<140 && i >=85) _9PB20[i] = -0.85+0.0015*(i-85);
	}
	//printf("90 %1.2f  %1.2f %1.2f\n",_9PB20[44],_9PB20[84],_9PB20[139]);

_9PB30(maxInd);
	for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
		if (i<45 && i>5) _9PB30[i] = -1.14 +0.0028*(i-5);
		else if (i<85 && i>=45) _9PB30[i] = -1.03+0.003*(i-45);
		else if (i<140 && i >=85) _9PB30[i] = -0.91+0.0017*(i-85);
	}
	//printf("90 %1.2f  %1.2f %1.2f\n",_9PB30[44],_9PB30[84],_9PB30[139]);

_9PB40(maxInd);
	for (int i=0; i<maxInd; i++) {//i = chromaticity  0==>140
		if (i<45 && i>5) _9PB40[i] = -1.16 +0.002*(i-5);
		else if (i<85 && i>=45) _9PB40[i] = -1.08+0.00275*(i-45);
		else if (i<140 && i >=85) _9PB40[i] = -0.97+0.0016*(i-85);
	}
	//printf("90 %1.2f  %1.2f %1.2f\n",_9PB40[44],_9PB40[84],_9PB40[139]);

_9PB50(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _9PB50[i] = -1.19 +0.00175*(i-5);
		else if (i<90 && i>=45) _9PB50[i] = -1.12+0.00225*(i-45);
	}
	//printf("90 %1.2f  %1.2f \n",_9PB50[44],_9PB50[84]);

_9PB60(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _9PB60[i] = -1.21 +0.0015*(i-5);
		else if (i<90 && i>=45) _9PB60[i] = -1.15+0.002*(i-45);
	}
	//printf("90 %1.2f  %1.2f \n",_9PB60[44],_9PB60[89]);
_9PB70(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _9PB70[i] = -1.23 +0.0018*(i-5);
	}
		//printf("9 %1.2f\n",_9PB70[49]);
_9PB80(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _9PB80[i] = -1.24 +0.002*(i-5);
	}
		//printf("9 %1.2f %1.2f\n",_9PB80[38], _9PB80[49]);


//10PB
_10PB10(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _10PB10[i] = -1.02 +0.00425*(i-5);
		else if (i<85 && i>=45) _10PB10[i] = -0.85+0.0025*(i-45);
		else if (i<140 && i >=85) _10PB10[i] = -0.75+0.0012*(i-85);
	}
	//printf("10 %1.2f  %1.2f %1.2f\n",_10PB10[44],_10PB10[84],_10PB10[139]);

_10PB20(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _10PB20[i] = -1.05 +0.00325*(i-5);
		else if (i<85 && i>=45) _10PB20[i] = -0.92+0.00275*(i-45);
		else if (i<140 && i >=85) _10PB20[i] = -0.81+0.0014*(i-85);
	}
	//printf("10 %1.2f  %1.2f %1.2f\n",_10PB20[44],_10PB20[84],_10PB20[139]);

_10PB30(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _10PB30[i] = -1.07 +0.00275*(i-5);
		else if (i<85 && i>=45) _10PB30[i] = -0.96+0.0025*(i-45);
		else if (i<140 && i >=85) _10PB30[i] = -0.86+0.0015*(i-85);
	}
	//printf("10 %1.2f  %1.2f %1.2f\n",_10PB30[44],_10PB30[84],_10PB30[139]);

_10PB40(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _10PB40[i] = -1.09 +0.002*(i-5);
		else if (i<85 && i>=45) _10PB40[i] = -1.01+0.00225*(i-45);
		else if (i<140 && i >=85) _10PB40[i] = -0.92+0.0016*(i-85);
	}
	//printf("10 %1.2f  %1.2f %1.2f\n",_10PB40[44],_10PB40[84],_10PB40[139]);

_10PB50(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _10PB50[i] = -1.12 +0.00175*(i-5);
		else if (i<90 && i>=45) _10PB50[i] = -1.05+0.00225*(i-45);
	}
	//printf("10 %1.2f  %1.2f\n",_10PB50[44],_10PB50[84]);

_10PB60(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _10PB60[i] = -1.14 +0.0015*(i-5);
		else if (i<90 && i>=45) _10PB60[i] = -1.08+0.00225*(i-45);
	}
	//printf("10 %1.2f  %1.2f\n",_10PB60[44],_10PB60[89]);


//1P
_1P10(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _1P10[i] = -0.96 +0.00375*(i-5);
		else if (i<85 && i>=45) _1P10[i] = -0.81+0.00225*(i-45);
		else if (i<140 && i >=85) _1P10[i] = -0.72+0.001*(i-85);
	}
	//printf("1P %1.2f  %1.2f %1.2f\n",_1P10[44],_1P10[84],_1P10[139]);

_1P20(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _1P20[i] = -1.0 +0.00325*(i-5);
		else if (i<85 && i>=45) _1P20[i] = -0.87+0.0025*(i-45);
		else if (i<140 && i >=85) _1P20[i] = -0.77+0.0012*(i-85);
	}
	//printf("1P %1.2f  %1.2f %1.2f\n",_1P20[44],_1P20[84],_1P20[139]);

_1P30(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _1P30[i] = -1.02 +0.00275*(i-5);
		else if (i<85 && i>=45) _1P30[i] = -0.91+0.00225*(i-45);
		else if (i<140 && i >=85) _1P30[i] = -0.82+0.0011*(i-85);
	}
	//printf("1P %1.2f  %1.2f %1.2f\n",_1P30[44],_1P30[84],_1P30[139]);

_1P40(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _1P40[i] = -1.04 +0.00225*(i-5);
		else if (i<85 && i>=45) _1P40[i] = -0.95+0.00225*(i-45);
		else if (i<140 && i >=85) _1P40[i] = -0.86+0.0015*(i-85);
	}
	//printf("1P %1.2f  %1.2f %1.2f\n",_1P40[44],_1P40[84],_1P40[139]);

_1P50(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _1P50[i] = -1.06 +0.002*(i-5);
		else if (i<90 && i>=45) _1P50[i] = -0.98+0.00175*(i-45);
	}
	//printf("1P %1.2f  %1.2f \n",_1P50[44],_1P50[89]);

_1P60(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _1P60[i] = -1.07 +0.0015*(i-5);
		else if (i<90 && i>=45) _1P60[i] = -1.01+0.00175*(i-45);
	}
	//printf("1P %1.2f  %1.2f \n",_1P60[44],_1P60[84],_1P60[139]);

		//4P
_4P10(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _4P10[i] = -0.78 +0.002*(i-5);
		else if (i<85 && i>=45) _4P10[i] = -0.7+0.00125*(i-45);
		else if (i<140 && i >=85) _4P10[i] = -0.65+0.001*(i-85);
	}
	//printf("4P %1.2f  %1.2f %1.2f\n",_4P10[44],_4P10[84],_4P10[139]);

_4P20(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _4P20[i] = -0.84 +0.0025*(i-5);
		else if (i<85 && i>=45) _4P20[i] = -0.74+0.00175*(i-45);
		else if (i<140 && i >=85) _4P20[i] = -0.67+0.00085*(i-85);
	}
	//printf("4P %1.2f  %1.2f %1.2f\n",_4P20[44],_4P20[84],_4P20[139]);

_4P30(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _4P30[i] = -0.85 +0.00225*(i-5);
		else if (i<85 && i>=45) _4P30[i] = -0.76+0.00125*(i-45);
		else if (i<140 && i >=85) _4P30[i] = -0.71+0.001*(i-85);
	}
	//printf("4P %1.2f  %1.2f %1.2f\n",_4P30[44],_4P30[84],_4P30[139]);

_4P40(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _4P40[i] = -0.87 +0.00175*(i-5);
		else if (i<85 && i>=45) _4P40[i] = -0.8+0.00175*(i-45);
		else if (i<140 && i >=85) _4P40[i] = -0.73+0.00075*(i-85);
	}
	//printf("4P %1.2f  %1.2f %1.2f\n",_4P40[44],_4P40[84],_4P40[139]);

_4P50(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _4P50[i] = -0.88 +0.0015*(i-5);
		else if (i<90 && i>=45) _4P50[i] = -0.82+0.0015*(i-45);
	}
	//printf("4P %1.2f  %1.2f \n",_4P50[44],_4P50[89]);

_4P60(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _4P60[i] = -0.89 +0.00125*(i-5);
		else if (i<90 && i>=45) _4P60[i] = -0.84+0.00125*(i-45);
	}
	//printf("4P %1.2f  %1.2f\n",_4P60[44],_4P60[89]);


//red yellow correction
_10YR20(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _10YR20[i] = 1.22 +0.002*(i-5);
		else if (i<90 && i>=45) _10YR20[i] = 1.30+0.006*(i-45);
	}
	//printf("10YR  %1.2f  %1.2f\n",_10YR20[44],_10YR20[56]);
_10YR30(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _10YR30[i] = 1.27 +0.00175*(i-5);
		else if (i<90 && i>=45) _10YR30[i] = 1.34+0.0017*(i-45);
	}
	//printf("10YR  %1.2f  %1.2f\n",_10YR30[44],_10YR30[75]);
_10YR40(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _10YR40[i] = 1.32 +0.00025*(i-5);
		else if (i<90 && i>=45) _10YR40[i] = 1.33+0.0015*(i-45);
	}
	//printf("10YR  %1.2f  %1.2f\n",_10YR40[44],_10YR40[85]);
_10YR50(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _10YR50[i] = 1.35 +0.000*(i-5);
		else if (i<90 && i>=45) _10YR50[i] = 1.35+0.0012*(i-45);
	}
	//printf("10YR  %1.2f  %1.2f\n",_10YR50[44],_10YR50[80]);
_10YR60(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _10YR60[i] = 1.38 - 0.00025*(i-5);
		else if (i<85 && i>=45) _10YR60[i] = 1.37+0.0005*(i-45);
		else if (i<140 && i >=85) _10YR60[i] = 1.39+0.0013*(i-85);
	}
	//printf("10YR  %1.2f  %1.2f %1.2f\n",_10YR60[44],_10YR60[85],_10YR60[139] );
_10YR70(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _10YR70[i] = 1.41 - 0.0005*(i-5);
		else if (i<85 && i>=45) _10YR70[i] = 1.39+0.000*(i-45);
		else if (i<140 && i >=85) _10YR70[i] = 1.39+0.0013*(i-85);
	}
	//printf("10YR  %1.2f  %1.2f %1.2f\n",_10YR70[44],_10YR70[85],_10YR70[139] );
_10YR80(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _10YR80[i] = 1.45 - 0.00125*(i-5);
		else if (i<85 && i>=45) _10YR80[i] = 1.40+0.000*(i-45);
		else if (i<140 && i >=85) _10YR80[i] = 1.40+0.00072*(i-85);//1.436
	}
	//printf("10YR  %1.2f  %1.2f %1.2f\n",_10YR80[44],_10YR80[84],_10YR80[139] );
_10YR90(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _10YR90[i] = 1.48 -0.001*(i-5);
		else if (i<90 && i>=45) _10YR90[i] = 1.44-0.0009*(i-45);
	}
	//printf("10YR  %1.2f  %1.2f\n",_10YR90[45],_10YR90[80]);
_85YR20(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _85YR20[i] = 1.12 +0.004*(i-5);
	}

	//printf("85YR  %1.2f \n",_85YR20[44]);
_85YR30(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _85YR30[i] = 1.16 + 0.0025*(i-5);
		else if (i<90 && i>=45) _85YR30[i] = 1.26+0.0028*(i-45);
	}
	//printf("85YR  %1.2f  %1.2f\n",_85YR30[44],_85YR30[75]);
_85YR40(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _85YR40[i] = 1.20 + 0.0015*(i-5);
		else if (i<90 && i>=45) _85YR40[i] = 1.26+0.0024*(i-45);
	}
	//printf("85YR  %1.2f  %1.2f\n",_85YR40[44],_85YR40[75]);
_85YR50(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _85YR50[i] = 1.24 + 0.0005*(i-5);
		else if (i<85 && i>=45) _85YR50[i] = 1.26+0.002*(i-45);
		else if (i<140 && i >=85) _85YR50[i] = 1.34+0.0015*(i-85);
	}
	//printf("85YR  %1.2f  %1.2f %1.2f\n",_85YR50[44],_85YR50[85],_85YR50[110] );
_85YR60(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _85YR60[i] = 1.27 + 0.00025*(i-5);
		else if (i<85 && i>=45) _85YR60[i] = 1.28+0.0015*(i-45);
		else if (i<140 && i >=85) _85YR60[i] = 1.34+0.0012*(i-85);
	}
	//printf("85YR  %1.2f  %1.2f %1.2f\n",_85YR60[44],_85YR60[85],_85YR60[139] );

_85YR70(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _85YR70[i] = 1.31 - 0.00025*(i-5);
		else if (i<85 && i>=45) _85YR70[i] = 1.30+0.0005*(i-45);
		else if (i<140 && i >=85) _85YR70[i] = 1.32+0.0012*(i-85);
	}
	//printf("85YR  %1.2f  %1.2f %1.2f\n",_85YR70[44],_85YR70[85],_85YR70[139] );
_85YR80(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _85YR80[i] = 1.35 - 0.00075*(i-5);
		else if (i<85 && i>=45) _85YR80[i] = 1.32+0.00025*(i-45);
		else if (i<140 && i >=85) _85YR80[i] = 1.33+0.00125*(i-85);
	}
	//printf("85YR  %1.2f  %1.2f %1.2f\n",_85YR80[44],_85YR80[85],_85YR80[139] );
_85YR90(maxInd2);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _85YR90[i] = 1.39 - 0.00125*(i-5);
		else if (i<90 && i>=45) _85YR90[i] = 1.34+0.00*(i-45);
	}
	//printf("85YR  %1.2f  %1.2f\n",_85YR90[44],_85YR90[85]);

//7YR
_7YR30(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _7YR30[i] = 1.06 + 0.0028*(i-5);
		else if (i<90 && i>=45) _7YR30[i] = 1.17+0.0045*(i-45);
	}
	//printf("7YR  %1.2f  %1.2f\n",_7YR30[44],_7YR30[66]);
_7YR40(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _7YR40[i] = 1.10 + 0.0018*(i-5);
		else if (i<90 && i>=45) _7YR40[i] = 1.17+0.0035*(i-45);
	}
	//printf("7YR  %1.2f  %1.2f\n",_7YR40[44],_7YR40[89]);
_7YR50(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _7YR50[i] = 1.14 + 0.00125*(i-5);
		else if (i<90 && i>=45) _7YR50[i] = 1.19+0.002*(i-45);
	}
	//printf("7YR  %1.2f  %1.2f\n",_7YR50[44],_7YR50[89] );
_7YR60(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _7YR60[i] = 1.17 + 0.00075*(i-5);
		else if (i<85 && i>=45) _7YR60[i] = 1.20+0.00175*(i-45);
		else if (i<140 && i >=85) _7YR60[i] = 1.27+0.002*(i-85);
	}
	//printf("7YR  %1.2f  %1.2f %1.2f\n",_7YR60[44],_7YR60[84],_7YR60[125] );

_7YR70(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _7YR70[i] = 1.20 + 0.0005*(i-5);
		else if (i<85 && i>=45) _7YR70[i] = 1.22+0.00125*(i-45);
		else if (i<140 && i >=85) _7YR70[i] = 1.27+0.0015*(i-85);
	}
	//printf("7YR  %1.2f  %1.2f %1.2f\n",_7YR70[44],_7YR70[84],_7YR70[125] );
_7YR80(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _7YR80[i] = 1.29 - 0.0008*(i-5);
	}
	//printf("7YR  %1.2f \n",_7YR80[44] );
_55YR30(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _55YR30[i] = 0.96 + 0.0038*(i-5);
	}
	//printf("55YR  %1.2f \n",_55YR30[44] );
_55YR40(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _55YR40[i] = 1.01 + 0.0022*(i-5);
		else if (i<90 && i>=45) _55YR40[i] = 1.10+0.0037*(i-45);
	}
	//printf("55YR  %1.2f  %1.2f\n",_55YR40[44],_55YR40[89] );
_55YR50(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _55YR50[i] = 1.06 + 0.0015*(i-5);
		else if (i<85 && i>=45) _55YR50[i] = 1.12+0.00225*(i-45);
		else if (i<140 && i >=85) _55YR50[i] = 1.21+0.0015*(i-85);
	}
	//printf("55YR  %1.2f  %1.2f %1.2f\n",_55YR50[44],_55YR50[84],_55YR50[125] );
_55YR60(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _55YR60[i] = 1.08 + 0.0012*(i-5);
		else if (i<85 && i>=45) _55YR60[i] = 1.13+0.0018*(i-45);
		else if (i<140 && i >=85) _55YR60[i] = 1.20+0.0025*(i-85);
	}
	//printf("55YR  %1.2f  %1.2f %1.2f\n",_55YR60[44],_55YR60[84],_55YR60[125] );
_55YR70(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _55YR70[i] = 1.11 + 0.00075*(i-5);
		else if (i<85 && i>=45) _55YR70[i] = 1.14+0.0012*(i-45);
		else if (i<140 && i >=85) _55YR70[i] = 1.19+0.00225*(i-85);
	}
	//printf("55YR  %1.2f  %1.2f %1.2f\n",_55YR70[44],_55YR70[84],_55YR70[125] );
_55YR80(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _55YR80[i] = 1.16 + 0.00*(i-5);
		else if (i<85 && i>=45) _55YR80[i] = 1.16+0.00075*(i-45);
		else if (i<140 && i >=85) _55YR80[i] = 1.19+0.00175*(i-85);
	}
	//printf("55YR  %1.2f  %1.2f %1.2f\n",_55YR80[44],_55YR80[84],_55YR80[125] );
_55YR90(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _55YR90[i] = 1.19 - 0.0005*(i-5);
	}
	//printf("55YR  %1.2f \n",_55YR90[44] );

_4YR30(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _4YR30[i] = 0.87 + 0.0035*(i-5);
		else if (i<90 && i>=45) _4YR30[i] = 1.01+0.0043*(i-45);
	}
	//printf("4YR  %1.2f  %1.2f\n",_4YR30[44],_4YR30[78] );
_4YR40(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _4YR40[i] = 0.92 + 0.0025*(i-5);
		else if (i<90 && i>=45) _4YR40[i] = 1.02+0.0033*(i-45);
	}
	//printf("4YR  %1.2f  %1.2f\n",_4YR40[44],_4YR40[74] );
_4YR50(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _4YR50[i] = 0.97 + 0.0015*(i-5);
		else if (i<90 && i>=45) _4YR50[i] = 1.03+0.0025*(i-45);
	}
	//printf("4YR  %1.2f  %1.2f\n",_4YR50[44],_4YR50[85] );
_4YR60(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _4YR60[i] = 0.99 + 0.00125*(i-5);
		else if (i<85 && i>=45) _4YR60[i] = 1.04+0.002*(i-45);
		else if (i<140 && i >=85) _4YR60[i] = 1.12+0.003*(i-85);
	}
	//printf("4YR  %1.2f  %1.2f %1.2f\n",_4YR60[44],_4YR60[84],_4YR60[125] );
_4YR70(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _4YR70[i] = 1.02 + 0.00075*(i-5);
		else if (i<85 && i>=45) _4YR70[i] = 1.05+0.00175*(i-45);
		else if (i<140 && i >=85) _4YR70[i] = 1.12+0.002*(i-85);
	}
	//printf("4YR  %1.2f  %1.2f %1.2f\n",_4YR70[44],_4YR70[84],_4YR70[125] );
_4YR80(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<50 && i>5) _4YR80[i] = 1.09 - 0.0002*(i-5);
	}
	//printf("4YR  %1.2f \n",_4YR80[41] );

_25YR30(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _25YR30[i] = 0.77 + 0.004*(i-5);
		else if (i<90 && i>=45) _25YR30[i] = 0.94+0.004*(i-45);
	}
	//printf("25YR  %1.2f  %1.2f\n",_25YR30[44],_25YR30[74] );
_25YR40(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _25YR40[i] = 0.82 + 0.003*(i-5);
		else if (i<90 && i>=45) _25YR40[i] = 0.94+0.002*(i-45);
	}
	//printf("25YR  %1.2f  %1.2f\n",_25YR40[44],_25YR40[84] );
_25YR50(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _25YR50[i] = 0.87+ 0.002*(i-5);
		else if (i<90 && i>=45) _25YR50[i] = 0.95+0.003*(i-45);
	}
	//printf("25YR  %1.2f  %1.2f\n",_25YR50[44],_25YR50[84] );
_25YR60(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _25YR60[i] = 0.89+ 0.0015*(i-5);
		else if (i<90 && i>=45) _25YR60[i] = 0.95+0.004*(i-45);
	}
	//printf("25YR  %1.2f  %1.2f\n",_25YR60[44],_25YR60[84] );
_25YR70(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _25YR70[i] = 0.92+ 0.001*(i-5);
		else if (i<90 && i>=45) _25YR70[i] = 0.96+0.003*(i-45);
	}
	//printf("25YR  %1.2f  %1.2f\n",_25YR70[44],_25YR70[84] );

_10R30(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _10R30[i] = 0.62 + 0.00225*(i-5);
		else if (i<90 && i>=45) _10R30[i] = 0.71+0.003*(i-45);
	}
	//printf("10R  %1.2f  %1.2f\n",_10R30[44],_10R30[84] );
_10R40(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _10R40[i] = 0.66 + 0.0025*(i-5);
		else if (i<90 && i>=45) _10R40[i] = 0.76+0.0035*(i-45);
	}
	//printf("10R  %1.2f  %1.2f\n",_10R40[44],_10R40[84] );
_10R50(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _10R50[i] = 0.71 + 0.002*(i-5);
		else if (i<90 && i>=45) _10R50[i] = 0.79+0.0043*(i-45);
	}
	//printf("10R  %1.2f  %1.2f\n",_10R50[44],_10R50[84] );
_10R60(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _10R60[i] = 0.73 + 0.00175*(i-5);
		else if (i<85 && i>=45) _10R60[i] = 0.80 +0.0033*(i-45);
		else if (i<140 && i >=85) _10R60[i] = 0.93+0.0018*(i-85);
	}
	//printf("10R  %1.2f  %1.2f %1.2f\n",_10R60[44],_10R60[84],_10R60[125] );
_10R70(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _10R70[i] = 0.75 + 0.0015*(i-5);
		else if (i<85 && i>=45) _10R70[i] = 0.81 +0.0017*(i-45);
		else if (i<140 && i >=85) _10R70[i] = 0.88+0.0025*(i-85);
	}
	//printf("10R  %1.2f  %1.2f %1.2f\n",_10R70[44],_10R70[84],_10R70[125] );

_9R30(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _9R30[i] = 0.57 + 0.002*(i-5);
		else if (i<90 && i>=45) _9R30[i] = 0.65+0.0018*(i-45);
	}
	//printf("9R  %1.2f  %1.2f\n",_9R30[44],_9R30[84] );
_9R40(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _9R40[i] = 0.61 + 0.002*(i-5);
		else if (i<90 && i>=45) _9R40[i] = 0.69+0.0025*(i-45);
	}
	//printf("9R  %1.2f  %1.2f\n",_9R40[44],_9R40[84] );
_9R50(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _9R50[i] = 0.66 + 0.00175*(i-5);
		else if (i<85 && i>=45) _9R50[i] = 0.73 +0.0025*(i-45);
		else if (i<140 && i >=85) _9R50[i] = 0.83+0.0035*(i-85);
	}
	//printf("9R  %1.2f  %1.2f %1.2f\n",_9R50[44],_9R50[84],_9R50[125] );
_9R60(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _9R60[i] = 0.68 + 0.0015*(i-5);
		else if (i<85 && i>=45) _9R60[i] = 0.74 +0.0022*(i-45);
		else if (i<140 && i >=85) _9R60[i] = 0.93+0.0022*(i-85);
	}
	//printf("9R  %1.2f  %1.2f %1.2f\n",_9R60[44],_9R60[84],_9R60[125] );
_9R70(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _9R70[i] = 0.70 + 0.0012*(i-5);
		else if (i<90 && i>=45) _9R70[i] = 0.75+0.0013*(i-45);
	}
	//printf("9R  %1.2f  %1.2f\n",_9R70[44],_9R70[84] );

_7R30(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _7R30[i] = 0.48 + 0.0015*(i-5);
		else if (i<90 && i>=45) _7R30[i] = 0.54-0.0005*(i-45);
	}
	//printf("7R  %1.2f  %1.2f\n",_7R30[44],_7R30[84] );
_7R40(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _7R40[i] = 0.51 + 0.0015*(i-5);
		else if (i<90 && i>=45) _7R40[i] = 0.57+0.0005*(i-45);
	}
	//printf("7R  %1.2f  %1.2f\n",_7R40[44],_7R40[84] );
_7R50(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _7R50[i] = 0.54 + 0.0015*(i-5);
		else if (i<85 && i>=45) _7R50[i] = 0.60 +0.0005*(i-45);
		else if (i<140 && i >=85) _7R50[i] = 0.62+0.0025*(i-85);
	}
	//printf("7R  %1.2f  %1.2f %1.2f\n",_7R50[44],_7R50[84],_7R50[125] );
_7R60(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _7R60[i] = 0.58 + 0.00075*(i-5);
		else if (i<85 && i>=45) _7R60[i] = 0.61 +0.00075*(i-45);
		else if (i<140 && i >=85) _7R60[i] = 0.64+0.001*(i-85);
	}
	//printf("7R  %1.2f  %1.2f %1.2f\n",_7R60[44],_7R60[84],_7R60[107] );
_7R70(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _7R70[i] = 0.59 + 0.00075*(i-5);
		else if (i<90 && i>=45) _7R70[i] = 0.62+0.00075*(i-45);
	}
	//printf("7R  %1.2f  %1.2f\n",_7R70[44],_7R70[84] );

//5R 1 2 3

//5R
_5R10(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _5R10[i] = 0.10 - 0.0018*(i-5);
		else if (i<90 && i>=45) _5R10[i] = 0.035-0.003*(i-45);
	}
	//printf("5R  %1.2f  %1.2f\n",_5R10[44],_5R10[51] );
_5R20(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _5R20[i] = 0.26 - 0.00075*(i-5);
		else if (i<90 && i>=45) _5R20[i] = 0.023-0.0002*(i-45);
	}
	//printf("5R  %1.2f  %1.2f\n",_5R20[44],_5R20[70] );
_5R30(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _5R30[i] = 0.39 + 0.00075*(i-5);
		else if (i<90 && i>=45) _5R30[i] = 0.42-0.0007*(i-45);
	}
	//printf("5R  %1.2f  %1.2f\n",_5R30[44],_5R30[85] );

//25R
_25R10(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<45 && i>5) _25R10[i] = -0.03 - 0.002*(i-5);
	}
	//printf("25R  %1.2f \n",_25R10[44]);
_25R20(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _25R20[i] = 0.13 - 0.0012*(i-5);
		else if (i<90 && i>=45) _25R20[i] = 0.08-0.002*(i-45);
	}
	//printf("25R  %1.2f  %1.2f\n",_25R20[44],_25R20[69] );
	//25R30: 0.28, 0.26, 0.22
_25R30(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _25R30[i] = 0.28 - 0.0005*(i-5);
		else if (i<90 && i>=45) _25R30[i] = 0.26-0.0009*(i-45);
	}
	//printf("25R  %1.2f  %1.2f\n",_25R30[44],_25R30[85] );


_10RP10(maxInd3);
	for (int i=0; i<maxInd3; i++) {
		if (i<45 && i>5) _10RP10[i] = -0.16 - 0.0017*(i-5);
	}
	//printf("10RP  %1.2f \n",_10RP10[44]);
_10RP20(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _10RP20[i] = 0.0 - 0.0018*(i-5);
		else if (i<90 && i>=45) _10RP20[i] = -0.07-0.0012*(i-45);
	}
	//printf("10RP  %1.2f  %1.2f\n",_10RP20[44],_10RP20[69] );
_10RP30(maxInd2);
	for (int i=0; i<maxInd2; i++) {
		if (i<45 && i>5) _10RP30[i] = 0.15 - 0.001*(i-5);
		else if (i<90 && i>=45) _10RP30[i] = 0.11-0.0012*(i-45);
	}
	//printf("10RP  %1.2f  %1.2f\n",_10RP30[44],_10RP30[85] );

//7G
_7G30(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _7G30[i] = 2.90 + 0.0027*(i-5);
		else if (i<85 && i>=45) _7G30[i] = 3.01+0.0005*(i-45);
		else if (i<140 && i >=85) _7G30[i] = 3.03+0.00075*(i-85);
	}
	//printf("7G  %1.2f  %1.2f %1.2f\n",_7G30[44],_7G30[84],_7G30[125] );
_7G40(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _7G40[i] = 2.89 + 0.00125*(i-5);
		else if (i<85 && i>=45) _7G40[i] = 2.94+0.0015*(i-45);
		else if (i<140 && i >=85) _7G40[i] = 3.0+0.001*(i-85);
	}
	//printf("7G  %1.2f  %1.2f %1.2f\n",_7G40[44],_7G40[84],_7G40[125] );
_7G50(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _7G50[i] = 2.87 + 0.0015*(i-5);
		else if (i<85 && i>=45) _7G50[i] = 2.93+0.00125*(i-45);
		else if (i<140 && i >=85) _7G50[i] = 2.98+0.001*(i-85);
	}
	//printf("7G  %1.2f  %1.2f %1.2f\n",_7G50[44],_7G50[84],_7G50[125] );
_7G60(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _7G60[i] = 2.86 + 0.00125*(i-5);
		else if (i<85 && i>=45) _7G60[i] = 2.91+0.00125*(i-45);
		else if (i<140 && i >=85) _7G60[i] = 2.96+0.00075*(i-85);
	}
	//printf("7G  %1.2f  %1.2f %1.2f\n",_7G60[44],_7G60[84],_7G60[125] );
_7G70(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _7G70[i] = 2.85 + 0.001*(i-5);
		else if (i<85 && i>=45) _7G70[i] = 2.89+0.00125*(i-45);
		else if (i<140 && i >=85) _7G70[i] = 2.94+0.00075*(i-85);
	}
	//printf("7G  %1.2f  %1.2f %1.2f\n",_7G70[44],_7G70[84],_7G70[125] );
_7G80(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _7G80[i] = 2.84 + 0.001*(i-5);
		else if (i<85 && i>=45) _7G80[i] = 2.88+0.001*(i-45);
		else if (i<140 && i >=85) _7G80[i] = 2.92+0.001*(i-85);
	}
	//printf("7G  %1.2f  %1.2f %1.2f\n",_7G80[44],_7G80[84],_7G80[125] );


//5G
_5G30(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _5G30[i] = 2.82 + 0.00175*(i-5);
		else if (i<85 && i>=45) _5G30[i] = 2.89+0.0018*(i-45);
		else if (i<140 && i >=85) _5G30[i] = 2.96+0.0012*(i-85);
	}
	//printf("5G  %1.2f  %1.2f %1.2f\n",_5G30[44],_5G30[84],_5G30[125] );
_5G40(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _5G40[i] = 2.80 + 0.0015*(i-5);
		else if (i<85 && i>=45) _5G40[i] = 2.86+0.00175*(i-45);
		else if (i<140 && i >=85) _5G40[i] = 2.93+0.00125*(i-85);
	}
	//printf("5G  %1.2f  %1.2f %1.2f\n",_5G40[44],_5G40[84],_5G40[125] );
_5G50(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _5G50[i] = 2.79 + 0.001*(i-5);
		else if (i<85 && i>=45) _5G50[i] = 2.84+0.0015*(i-45);
		else if (i<140 && i >=85) _5G50[i] = 2.90+0.0015*(i-85);
	}
	//printf("5G  %1.2f  %1.2f %1.2f\n",_5G50[44],_5G50[84],_5G50[125] );
_5G60(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _5G60[i] = 2.78 + 0.001*(i-5);
		else if (i<85 && i>=45) _5G60[i] = 2.82+0.00175*(i-45);
		else if (i<140 && i >=85) _5G60[i] = 2.89+0.001*(i-85);
	}
	//printf("5G  %1.2f  %1.2f %1.2f\n",_5G60[44],_5G60[84],_5G60[125] );
_5G70(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _5G70[i] = 2.77 + 0.001*(i-5);
		else if (i<85 && i>=45) _5G70[i] = 2.81+0.00125*(i-45);
		else if (i<140 && i >=85) _5G70[i] = 2.86+0.00125*(i-85);
	}
	//printf("5G  %1.2f  %1.2f %1.2f\n",_5G70[44],_5G70[84],_5G70[125] );
_5G80(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _5G80[i] = 2.76 + 0.001*(i-5);
		else if (i<85 && i>=45) _5G80[i] = 2.8+0.00125*(i-45);
		else if (i<140 && i >=85) _5G80[i] = 2.85+0.00125*(i-85);
	}
	//printf("5G  %1.2f  %1.2f %1.2f\n",_5G80[44],_5G80[84],_5G80[125] );

//25G
_25G30(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _25G30[i] = 2.68 + 0.0015*(i-5);
		else if (i<85 && i>=45) _25G30[i] = 2.74+0.0018*(i-45);
		else if (i<140 && i >=85) _25G30[i] = 2.81+0.002*(i-85);
	}
	//printf("25G  %1.2f  %1.2f %1.2f\n",_25G30[44],_25G30[84],_25G30[125] );
_25G40(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _25G40[i] = 2.68 + 0.00075*(i-5);
		else if (i<85 && i>=45) _25G40[i] = 2.71+0.0015*(i-45);
		else if (i<140 && i >=85) _25G40[i] = 2.77+0.00125*(i-85);
	}
	//printf("25G  %1.2f  %1.2f %1.2f\n",_25G40[44],_25G40[84],_25G40[125] );
_25G50(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _25G50[i] = 2.65 + 0.00075*(i-5);
		else if (i<85 && i>=45) _25G50[i] = 2.68+0.00125*(i-45);
		else if (i<140 && i >=85) _25G50[i] = 2.73+0.00125*(i-85);
	}
	//printf("25G  %1.2f  %1.2f %1.2f\n",_25G50[44],_25G50[84],_25G50[125] );
_25G60(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _25G60[i] = 2.64 + 0.0005*(i-5);
		else if (i<85 && i>=45) _25G60[i] = 2.66+0.001*(i-45);
		else if (i<140 && i >=85) _25G60[i] = 2.70+0.001*(i-85);
	}
	//printf("25G  %1.2f  %1.2f %1.2f\n",_25G60[44],_25G60[84],_25G60[125] );
_25G70(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _25G70[i] = 2.64 + 0.00*(i-5);
		else if (i<85 && i>=45) _25G70[i] = 2.64+0.00075*(i-45);
		else if (i<140 && i >=85) _25G70[i] = 2.67+0.001*(i-85);
	}
	//printf("25G  %1.2f  %1.2f %1.2f\n",_25G70[44],_25G70[84],_25G70[125] );
_25G80(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _25G80[i] = 2.63 + 0.00*(i-5);
		else if (i<85 && i>=45) _25G80[i] = 2.63+0.0005*(i-45);
		else if (i<140 && i >=85) _25G80[i] = 2.65+0.0005*(i-85);
	}
	//printf("25G  %1.2f  %1.2f %1.2f\n",_25G80[44],_25G80[84],_25G80[125] );


//1G
_1G30(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _1G30[i] = 2.58 + 0.00025*(i-5);
		else if (i<85 && i>=45) _1G30[i] = 2.59+0.001*(i-45);
		else if (i<140 && i >=85) _1G30[i] = 2.63+0.00125*(i-85);
	}
	//printf("1G  %1.2f  %1.2f %1.2f\n",_1G30[44],_1G30[84],_1G30[125] );
_1G40(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _1G40[i] = 2.56 - 0.00025*(i-5);
		else if (i<85 && i>=45) _1G40[i] = 2.55+0.0005*(i-45);
		else if (i<140 && i >=85) _1G40[i] = 2.57+0.0005*(i-85);
	}
	//printf("1G  %1.2f  %1.2f %1.2f\n",_1G40[44],_1G40[84],_1G40[125] );
_1G50(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _1G50[i] = 2.55 - 0.00025*(i-5);
		else if (i<85 && i>=45) _1G50[i] = 2.54+0.00025*(i-45);
		else if (i<140 && i >=85) _1G50[i] = 2.55+0.0005*(i-85);
	}
	//printf("1G  %1.2f  %1.2f %1.2f\n",_1G50[44],_1G50[84],_1G50[125] );
_1G60(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _1G60[i] = 2.54 - 0.0005*(i-5);
		else if (i<85 && i>=45) _1G60[i] = 2.52+0.00025*(i-45);
		else if (i<140 && i >=85) _1G60[i] = 2.53+0.00025*(i-85);
	}
	//printf("1G  %1.2f  %1.2f %1.2f\n",_1G60[44],_1G60[84],_1G60[125] );
_1G70(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _1G70[i] = 2.53 - 0.0005*(i-5);
		else if (i<85 && i>=45) _1G70[i] = 2.51+0.0*(i-45);
		else if (i<140 && i >=85) _1G70[i] = 2.51+0.00025*(i-85);
	}
	//printf("1G  %1.2f  %1.2f %1.2f\n",_1G70[44],_1G70[84],_1G70[125] );
_1G80(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _1G80[i] = 2.52 - 0.0005*(i-5);
		else if (i<85 && i>=45) _1G80[i] = 2.50+0.00*(i-45);
		else if (i<140 && i >=85) _1G80[i] = 2.50+0.00*(i-85);
	}
	//printf("1G  %1.2f  %1.2f %1.2f\n",_1G80[44],_1G80[84],_1G80[125] );


//10GY
_10GY30(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _10GY30[i] = 2.52 - 0.001*(i-5);
		else if (i<85 && i>=45) _10GY30[i] = 2.48-0.002*(i-45);
		else if (i<140 && i >=85) _10GY30[i] = 2.40+0.0025*(i-85);
	}
	//printf("10GY  %1.2f  %1.2f %1.2f\n",_10GY30[44],_10GY30[84],_10GY30[125] );
_10GY40(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _10GY40[i] = 2.48 - 0.0005*(i-5);
		else if (i<85 && i>=45) _10GY40[i] = 2.46-0.0005*(i-45);
		else if (i<140 && i >=85) _10GY40[i] = 2.44-0.0015*(i-85);
	}
	//printf("10GY  %1.2f  %1.2f %1.2f\n",_10GY40[44],_10GY40[84],_10GY40[125] );
_10GY50(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _10GY50[i] = 2.48 - 0.00075*(i-5);
		else if (i<85 && i>=45) _10GY50[i] = 2.45-0.00075*(i-45);
		else if (i<140 && i >=85) _10GY50[i] = 2.42-0.00175*(i-85);
	}
	//printf("10GY  %1.2f  %1.2f %1.2f\n",_10GY50[44],_10GY50[84],_10GY50[125] );
_10GY60(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _10GY60[i] = 2.47 - 0.00125*(i-5);
		else if (i<85 && i>=45) _10GY60[i] = 2.42-0.00025*(i-45);
		else if (i<140 && i >=85) _10GY60[i] = 2.41-0.0005*(i-85);
	}
	//printf("10GY  %1.2f  %1.2f %1.2f\n",_10GY60[44],_10GY60[84],_10GY60[125] );
_10GY70(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _10GY70[i] = 2.46 - 0.001*(i-5);
		else if (i<85 && i>=45) _10GY70[i] = 2.42+0.0*(i-45);
		else if (i<140 && i >=85) _10GY70[i] = 2.42-0.001*(i-85);
	}
	//printf("10GY %1.2f  %1.2f %1.2f\n",_10GY70[44],_10GY70[84],_10GY70[125] );
_10GY80(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _10GY80[i] = 2.45 - 0.00075*(i-5);
		else if (i<85 && i>=45) _10GY80[i] = 2.42 - 0.0005*(i-45);
		else if (i<140 && i >=85) _10GY80[i] = 2.40-0.0005*(i-85);
	}
	//printf("10GY  %1.2f  %1.2f %1.2f\n",_10GY80[44],_10GY80[84],_10GY80[125] );


//75GY
_75GY30(maxInd2);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _75GY30[i] = 2.36 - 0.0025*(i-5);
		else if (i<90 && i>=45) _75GY30[i] = 2.26-0.00175*(i-45);
	}
	//printf("75GY  %1.2f  %1.2f\n",_75GY30[44],_75GY30[84] );
_75GY40(maxInd2);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _75GY40[i] = 2.34 - 0.00175*(i-5);
		else if (i<90 && i>=45) _75GY40[i] = 2.27-0.00225*(i-45);
	}
	//printf("75GY  %1.2f  %1.2f \n",_75GY40[44],_75GY40[84] );
_75GY50(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _75GY50[i] = 2.32 - 0.0015*(i-5);
		else if (i<85 && i>=45) _75GY50[i] = 2.26-0.00175*(i-45);
		else if (i<140 && i >=85) _75GY50[i] = 2.19-0.00325*(i-85);
	}
	//printf("75GY  %1.2f  %1.2f %1.2f %1.2f\n",_75GY50[44],_75GY50[84],_75GY50[125],_75GY50[139] );
_75GY60(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _75GY60[i] = 2.30 - 0.00125*(i-5);
		else if (i<85 && i>=45) _75GY60[i] = 2.25-0.001*(i-45);
		else if (i<140 && i >=85) _75GY60[i] = 2.21-0.0027*(i-85);
	}
	//printf("75GY  %1.2f  %1.2f %1.2f\n",_75GY60[44],_75GY60[84],_75GY60[125] );
_75GY70(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _75GY70[i] = 2.29 - 0.00125*(i-5);
		else if (i<85 && i>=45) _75GY70[i] = 2.24-0.0015*(i-45);
		else if (i<140 && i >=85) _75GY70[i] = 2.18-0.00175*(i-85);
	}
	//printf("75GY %1.2f  %1.2f %1.2f\n",_75GY70[44],_75GY70[84],_75GY70[125] );
_75GY80(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _75GY80[i] = 2.27 - 0.001*(i-5);
		else if (i<85 && i>=45) _75GY80[i] = 2.23 - 0.001*(i-45);
		else if (i<140 && i >=85) _75GY80[i] = 2.19-0.00175*(i-85);
	}
	//printf("75GY  %1.2f  %1.2f %1.2f\n",_75GY80[44],_75GY80[84],_75GY80[125] );


//55GY
_5GY30(maxInd2);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _5GY30[i] = 2.16 - 0.002*(i-5);
		else if (i<90 && i>=45) _5GY30[i] = 2.07-0.0025*(i-45);
	}
	//printf("5GY  %1.2f  %1.2f\n",_5GY30[44],_5GY30[84] );

//5GY4: 2.14,2.04, 1.96, 1.91 //95

_5GY40(maxInd2);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _5GY40[i] = 2.14 - 0.0025*(i-5);
		else if (i<90 && i>=45) _5GY40[i] = 2.04-0.003*(i-45);
	}
	//printf("5GY  %1.2f  %1.2f \n",_5GY40[44],_5GY40[84] );
_5GY50(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _5GY50[i] = 2.13 - 0.00175*(i-5);
		else if (i<85 && i>=45) _5GY50[i] = 2.06-0.002*(i-45);
		else if (i<140 && i >=85) _5GY50[i] = 1.98-0.00225*(i-85);
	}
	//printf("5GY  %1.2f  %1.2f %1.2f\n",_5GY50[44],_5GY50[84],_5GY50[125] );
_5GY60(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _5GY60[i] = 2.11 - 0.0015*(i-5);
		else if (i<85 && i>=45) _5GY60[i] = 2.05-0.002*(i-45);
		else if (i<140 && i >=85) _5GY60[i] = 1.97-0.00275*(i-85);
	}
	//printf("5GY  %1.2f  %1.2f %1.2f\n",_5GY60[44],_5GY60[84],_5GY60[125] );
_5GY70(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _5GY70[i] = 2.09 - 0.001*(i-5);
		else if (i<85 && i>=45) _5GY70[i] = 2.05-0.00175*(i-45);
		else if (i<140 && i >=85) _5GY70[i] = 1.98-0.002*(i-85);
	}
	//printf("5GY %1.2f  %1.2f %1.2f\n",_5GY70[44],_5GY70[84],_5GY70[125] );
_5GY80(maxInd);
	for (int i=0; i<maxInd; i++) {
		if (i<45 && i>5) _5GY80[i] = 2.07 - 0.001*(i-5);
		else if (i<85 && i>=45) _5GY80[i] = 2.03 - 0.00075*(i-45);
		else if (i<140 && i >=85) _5GY80[i] = 2.0-0.002*(i-85);
	}
	//printf("5GY  %1.2f  %1.2f %1.2f\n",_5GY80[44],_5GY80[84],_5GY80[125] );

#ifdef _DEBUG
	t2e.set();
	if (settings->verbose)
		printf("Lutf Munsell  %d usec\n", t2e.etime(t1e));
#endif
}



void ImProcFunctions::MunsellLch (float lum, float hue, float chrom, float memChprov, float &correction, int zone) {
	int x=(int) memChprov;
	int y=(int) chrom;

	//found the good LUT and calculate correction


	//Blue purple correction PB + sky
	if(zone==1) {//begin PB correction
		if((lum > 5.0 && lum <15.0) ) {
			if( (hue >= (_15PB10[x] - 0.035)) && (hue < (_15PB10[x] + 0.052) && x<=45)) {if(y>49) y=49;correction =  _15PB10[y] - _15PB10[x] ;}
			else if (( hue>=( _3PB10[x] -0.052))  && (hue < (_45PB10[x] + _3PB10[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _3PB10[y] - _3PB10[x] ;}
			else if (( hue>=(_45PB10[x] + _3PB10[x])/2.0)  && (hue < (_45PB10[x] +0.052)) && x <= 85) {if(y>89) y=89;correction =  _45PB10[y] - _45PB10[x] ;}
			else if (( hue>=(_6PB10[x] -0.052)  && (hue < (_6PB10[x] + _75PB10[x])/2.0))) {correction =  _6PB10[y] - _6PB10[x] ; }
			else if (( hue>=(_6PB10[x] + _75PB10[x])/2.0)  && (hue < (_9PB10[x] + _75PB10[x])/2.0)) {correction =  _75PB10[y] - _75PB10[x] ;}
			else if (( hue>=(_9PB10[x] + _75PB10[x])/2.0)  && (hue < (_9PB10[x] + _10PB10[x])/2.0)) {correction =  _9PB10[y] - _9PB10[x] ; }
			else if (( hue>=(_10PB10[x] + _9PB10[x])/2.0)  && (hue < (_1P10[x] + _10PB10[x])/2.0)) {correction =  _10PB10[y] - _10PB10[x] ;}
			else if (( hue>=(_10PB10[x] + _1P10[x])/2.0)  && (hue < (_1P10[x] + _4P10[x])/2.0)) {correction =  _1P10[y] - _1P10[x];}
			else if (( hue>=(_1P10[x] + _4P10[x])/2.0)  && (hue < (0.035 + _4P10[x])/2.0)) {correction =  _4P10[y] - _4P10[x] ;}
		}
		else if ((lum >= 15.0 && lum <25.0)) {
			if( (hue >= (_15PB20[x] - 0.035)) && (hue < (_15PB20[x] + _3PB20[x])/2.0) && x<=85) {if(y>89) y=89;correction =  _15PB20[y] - _15PB20[x] ; }
			else if (( hue>=(_15PB20[x] + _3PB20[x])/2.0)  && (hue < (_45PB20[x] + _3PB20[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _3PB20[y] - _3PB20[x] ; }
			else if (( hue>=(_45PB20[x] + _3PB20[x])/2.0)  && (hue < ( _45PB20[x] + 0.052)) && x <= 85) {if(y>89) y=89;correction =  _45PB20[y] - _45PB20[x] ;}
			else if (( hue>=(_45PB20[x] + 0.052))  && (hue < (_6PB20[x] + _75PB20[x])/2.0)) {correction =  _6PB20[y] - _6PB20[x];}
			else if (( hue>=(_6PB20[x] + _75PB20[x])/2.0)  && (hue < (_9PB20[x] + _75PB20[x])/2.0)) {correction =  _75PB20[y] - _75PB20[x] ;}
			else if (( hue>=(_9PB20[x] + _75PB20[x])/2.0)  && (hue < (_9PB20[x] + _10PB20[x])/2.0)) {correction =  _9PB20[y] - _9PB20[x] ; }
			else if (( hue>=(_10PB20[x] + _9PB20[x])/2.0)  && (hue < (_1P20[x] + _10PB20[x])/2.0)) {correction =  _10PB20[y] - _10PB20[x] ;}
			else if (( hue>=(_10PB20[x] + _1P20[x])/2.0)  && (hue < (_1P20[x] + _4P20[x])/2.0)) {correction =  _1P20[y] - _1P20[x] ; }
			else if (( hue>=(_1P20[x] + _4P20[x])/2.0)  && (hue < (0.035 + _4P20[x])/2.0)) {correction =  _4P20[y] - _4P20[x] ; }
		}
		else if ((lum >= 25.0 && lum <35.0)) {
			if( (hue >= (_15PB30[x] - 0.035)) && (hue < (_15PB30[x] + _3PB30[x])/2.0) && x<=85 ) {if(y>89) y=89;correction =  _15PB30[y] - _15PB30[x] ;}
			else if (( hue>=(_15PB30[x] + _3PB30[x])/2.0)  && (hue < (_45PB30[x] + _3PB30[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _3PB30[y] - _3PB30[x] ;}
			else if (( hue>=(_45PB30[x] + _3PB30[x])/2.0)  && (hue < (_45PB30[x]+0.052)) && x <= 85) {if(y>89) y=89;correction =  _45PB30[y] - _45PB30[x] ; }
			else if (( hue>=( _45PB30[x]+ 0.052))  && (hue < (_6PB30[x] + _75PB30[x])/2.0)) {correction =  _6PB30[y] - _6PB30[x] ; }
			else if (( hue>=(_6PB30[x] + _75PB30[x])/2.0)  && (hue < (_9PB30[x] + _75PB30[x])/2.0)) {correction =  _75PB30[y] - _75PB30[x] ;}
			else if (( hue>=(_9PB30[x] + _75PB30[x])/2.0)  && (hue < (_9PB30[x] + _10PB30[x])/2.0)) {correction =  _9PB30[y] - _9PB30[x] ; }
			else if (( hue>=(_10PB30[x] + _9PB30[x])/2.0)  && (hue < (_1P30[x] + _10PB30[x])/2.0)) {correction =  _10PB30[y] - _10PB30[x] ;}
			else if (( hue>=(_10PB30[x] + _1P30[x])/2.0)  && (hue < (_1P30[x] + _4P30[x])/2.0)) {correction =  _1P30[y] - _1P30[x] ; }
			else if (( hue>=(_1P30[x] + _4P30[x])/2.0)  && (hue < (0.035 + _4P30[x])/2.0)) {correction =  _4P30[y] - _4P30[x] ;}
		}
		else if ((lum >= 35.0 && lum <45.0) ) {
			if( (hue <= (_05PB40[x] + _15PB40[x])/2.0) && (hue > (_05PB40[x] + _10B40[x])/2.0) && x<75 ) {if(y>75) y=75; correction =  _05PB40[y] - _05PB40[x] ;}
			else if( (hue <= (_05PB40[x] + _10B40[x])/2.0) && (hue >(_10B40[x] + _9B40[x])/2.0) && x<70 ) {if(y>70) y=70;correction =  _10B40[y] - _10B40[x] ;}
			else if( (hue <= (_10B40[x] + _9B40[x])/2.0) && (hue >(_9B40[x] + _7B40[x])/2.0) && x<70 ) {if(y>70) y=70;correction =  _9B40[y] - _9B40[x] ;}
			else if( (hue <= (_9B40[x] + _7B40[x])/2.0) && (hue >(_5B40[x] + _7B40[x])/2.0) && x<70 ) {if(y>70) y=70;correction =  _7B40[y] - _7B40[x] ;}
			else if (( hue<=(_5B40[x] + _7B40[x])/2.0)  && (hue > (_5B40[x]-0.035)) && x < 70) {if(y>70) y=70; correction =  _5B40[y] - _5B40[x] ; }	//

			else if( (hue >= (_15PB40[x] - 0.035)) && (hue < (_15PB40[x] + _3PB40[x])/2.0) && x<=85 ) {if(y>89) y=89;correction =  _15PB40[y] - _15PB40[x] ; }
			else if (( hue>=(_15PB40[x] + _3PB40[x])/2.0)  && (hue < (_45PB40[x] + _3PB40[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _3PB40[y] - _3PB40[x] ;}
			else if (( hue>=(_45PB40[x] + _3PB40[x])/2.0)  && (hue < (_45PB40[x]+0.052)) && x <= 85) {if(y>89) y=89;correction =  _45PB40[y] - _45PB40[x] ;}
			else if (( hue>=(_45PB40[x]+0.052))  && (hue < (_6PB40[x] + _75PB40[x])/2.0)) {correction =  _6PB40[y] - _6PB40[x] ; }
			else if (( hue>=(_6PB40[x] + _75PB40[x])/2.0)  && (hue < (_9PB40[x] + _75PB40[x])/2.0)) {correction =  _75PB40[y] - _75PB40[x] ; }
			else if (( hue>=(_9PB40[x] + _75PB40[x])/2.0)  && (hue < (_9PB40[x] + _10PB40[x])/2.0)) {correction =  _9PB40[y] - _9PB40[x] ; }
			else if (( hue>=(_10PB40[x] + _9PB40[x])/2.0)  && (hue < (_1P40[x] + _10PB40[x])/2.0)) {correction =  _10PB40[y] - _10PB40[x] ;}
			else if (( hue>=(_10PB40[x] + _1P40[x])/2.0)  && (hue < (_1P40[x] + _4P40[x])/2.0)) {correction =  _1P40[y] - _1P40[x] ;}
			else if (( hue>=(_1P40[x] + _4P40[x])/2.0)  && (hue < (0.035 + _4P40[x])/2.0)) {correction =  _4P40[y] - _4P40[x] ;}
		}
		else if ((lum >= 45.0 && lum <55.0) ) {
			if( (hue <= (_05PB50[x] + _15PB50[x])/2.0) && (hue > (_05PB50[x] + _10B50[x])/2.0) && x<79 ) {if(y>79) y=79; correction =  _05PB50[y] - _05PB50[x] ;}
			else if( (hue <= (_05PB50[x] + _10B50[x])/2.0) && (hue >(_10B50[x] + _9B50[x])/2.0) && x<79 ) {if(y>79) y=79;correction =  _10B50[y] - _10B50[x] ;}
			else if( (hue <= (_10B50[x] + _9B50[x])/2.0) && (hue >(_9B50[x] + _7B50[x])/2.0) && x<79 ) {if(y>79) y=79;correction =  _9B50[y] - _9B50[x] ;}
			else if( (hue <= (_9B50[x] + _7B50[x])/2.0) && (hue >(_5B50[x] + _7B50[x])/2.0) && x<79 ) {if(y>79) y=79;correction =  _7B50[y] - _7B50[x] ;}
			else if (( hue<=(_5B50[x] + _7B50[x])/2.0)  && (hue > (_5B50[x]-0.035)) && x < 79) {if(y>79) y=79; correction =  _5B50[y] - _5B50[x] ; }	//

			else if( (hue >= (_15PB50[x] - 0.035)) && (hue < (_15PB50[x] + _3PB50[x])/2.0) && x<=85 ) {if(y>89) y=89;correction =  _15PB50[y] - _15PB50[x] ; }
			else if (( hue>=(_15PB50[x] + _3PB50[x])/2.0)  && (hue < (_45PB50[x] + _3PB50[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _3PB50[y] - _3PB50[x] ;}
			else if (( hue>=(_45PB50[x] + _3PB50[x])/2.0)  && (hue < (_6PB50[x] + _45PB50[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _45PB50[y] - _45PB50[x] ; }
			else if (( hue>=(_6PB50[x] + _45PB50[x])/2.0)  && (hue < (_6PB50[x] + _75PB50[x])/2.0) && x <=85) {if(y>89) y=89;correction =  _6PB50[y] - _6PB50[x] ;}
			else if (( hue>=(_6PB50[x] + _75PB50[x])/2.0)  && (hue < (_9PB50[x] + _75PB50[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _75PB50[y] - _75PB50[x] ;}
			else if (( hue>=(_9PB50[x] + _75PB50[x])/2.0)  && (hue < (_9PB50[x] + _10PB50[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _9PB50[y] - _9PB50[x] ;}
			else if (( hue>=(_10PB50[x] + _9PB50[x])/2.0)  && (hue < (_1P50[x] + _10PB50[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _10PB50[y] - _10PB50[x] ;}
			else if (( hue>=(_10PB50[x] + _1P50[x])/2.0)  && (hue < (_1P50[x] + _4P50[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _1P50[y] - _1P50[x] ; }
			else if (( hue>=(_1P50[x] + _4P50[x])/2.0)  && (hue < (0.035 + _4P50[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _4P50[y] - _4P50[x] ;}
		}
		else if ((lum >= 55.0 && lum <65.0) ) {
			if( (hue <= (_05PB60[x] + _15PB60[x])/2.0) && (hue > (_05PB60[x] + _10B60[x])/2.0) && x<79 ) {if(y>79) y=79; correction =  _05PB60[y] - _05PB60[x] ;}
			else if( (hue <= (_05PB60[x] + _10B60[x])/2.0) && (hue >(_10B60[x] + _9B60[x])/2.0) && x<79 ) {if(y>79) y=79;correction =  _10B60[y] - _10B60[x] ;}
			else if( (hue <= (_10B60[x] + _9B60[x])/2.0) && (hue >(_9B60[x] + _7B60[x])/2.0) && x<79 ) {if(y>79) y=79;correction =  _9B60[y] - _9B60[x] ;}
			else if( (hue <= (_9B60[x] + _7B60[x])/2.0) && (hue >(_5B60[x] + _7B60[x])/2.0) && x<79 ) {if(y>79) y=79;correction =  _7B60[y] - _7B60[x] ;}
			else if (( hue<=(_5B60[x] + _7B60[x])/2.0)  && (hue > (_5B60[x]-0.035)) && x < 79) {if(y>79) y=79; correction =  _5B60[y] - _5B60[x] ; }	//

			else if( (hue >= (_15PB60[x] - 0.035)) && (hue < (_15PB60[x] + _3PB60[x])/2.0) && x<=85 ) {if(y>89) y=89;correction =  _15PB60[y] - _15PB60[x] ; }
			else if (( hue>=(_15PB60[x] + _3PB60[x])/2.0)  && (hue < (_45PB60[x] + _3PB60[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _3PB60[y] - _3PB60[x] ;}
			else if (( hue>=(_45PB60[x] + _3PB60[x])/2.0)  && (hue < (_6PB60[x] + _45PB60[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _45PB60[y] - _45PB60[x] ;}
			else if (( hue>=(_6PB60[x] + _45PB60[x])/2.0)  && (hue < (_6PB60[x] + _75PB60[x])/2.0) && x <=85) {if(y>89) y=89;correction =  _6PB60[y] - _6PB60[x] ;}
			else if (( hue>=(_6PB60[x] + _75PB60[x])/2.0)  && (hue < (_9PB60[x] + _75PB60[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _75PB60[y] - _75PB60[x] ; }
			else if (( hue>=(_9PB60[x] + _75PB60[x])/2.0)  && (hue < (_9PB60[x] + _10PB60[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _9PB60[y] - _9PB60[x] ; }
			else if (( hue>=(_10PB60[x] + _9PB60[x])/2.0)  && (hue < (_1P60[x] + _10PB60[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _10PB60[y] - _10PB60[x] ; }
			else if (( hue>=(_10PB60[x] + _1P60[x])/2.0)  && (hue < (_1P60[x] + _4P60[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _1P60[y] - _1P60[x] ; }
			else if (( hue>=(_1P60[x] + _4P60[x])/2.0)  && (hue < (0.035 + _4P60[x])/2.0) && x <= 85) {if(y>89) y=89;correction =  _4P60[y] - _4P60[x] ; }
		}
		else if ((lum >= 65.0 && lum < 75.0) ) {
			if( (hue <= (_05PB70[x] + _15PB70[x])/2.0) && (hue > (_05PB70[x] + _10B70[x])/2.0) && x<50 ) {if(y>49) y=49; correction =  _05PB70[y] - _05PB70[x] ;}
			else if( (hue <= (_05PB70[x] + _10B70[x])/2.0) && (hue >(_10B70[x] + _9B70[x])/2.0) && x<50 ) {if(y>49) y=49;correction =  _10B70[y] - _10B70[x] ;}
			else if( (hue <= (_10B70[x] + _9B70[x])/2.0) && (hue >(_9B70[x] + _7B70[x])/2.0) && x<50 ) {if(y>49) y=49;correction =  _9B70[y] - _9B70[x] ;}
			else if( (hue <= (_9B70[x] + _7B70[x])/2.0) && (hue >(_5B70[x] + _7B70[x])/2.0) && x<50 ) {if(y>49) y=49;correction =  _7B70[y] - _7B70[x] ;}
			else if (( hue<=(_5B70[x] + _7B70[x])/2.0)  && (hue > (_5B70[x]-0.035)) && x < 50) {if(y>49) y=49; correction =  _5B70[y] - _5B70[x] ; }	//

			else if( (hue >= (_15PB70[x] - 0.035)) && (hue < (_15PB70[x] + _3PB70[x])/2.0) && x<50 ) {if(y>49) y=49;correction =  _15PB70[y] - _15PB70[x] ; }
			else if (( hue>=(_45PB70[x] + _3PB70[x])/2.0)  && (hue < (_6PB70[x] + _45PB70[x])/2.0) && x < 50) {if(y>49) y=49;correction =  _45PB70[y] - _45PB70[x] ;}
			else if (( hue>=(_6PB70[x] + _45PB70[x])/2.0)  && (hue < (_6PB70[x] + _75PB70[x])/2.0) && x <50) {if(y>49) y=49;correction =  _6PB70[y] - _6PB70[x] ;}
			else if (( hue>=(_6PB70[x] + _75PB70[x])/2.0)  && (hue < (_9PB70[x] + _75PB70[x])/2.0) && x <50) {if(y>49) y=49;correction =  _75PB70[y] - _75PB70[x] ; }
			else if (( hue>=(_9PB70[x] + _75PB70[x])/2.0)  && (hue < (_9PB70[x] + 0.035)) && x <50) {if(y>49) y=49;correction =  _9PB70[y] - _9PB70[x] ; }
		}
		else if ((lum >= 75.0 && lum < 85.0) ) {
			if( (hue <= (_05PB80[x] + _15PB80[x])/2.0) && (hue > (_05PB80[x] + _10B80[x])/2.0) && x<40 ) {if(y>39) y=39; correction =  _05PB80[y] - _05PB80[x] ;}
			else if( (hue <= (_05PB80[x] + _10B80[x])/2.0) && (hue >(_10B80[x] + _9B80[x])/2.0) && x<40 ) {if(y>39) y=39;correction =  _10B80[y] - _10B80[x] ;}
			else if( (hue <= (_10B80[x] + _9B80[x])/2.0) && (hue >(_9B80[x] + _7B80[x])/2.0) && x<40 ) {if(y>39) y=39;correction =  _9B80[y] - _9B80[x] ;}
			else if( (hue <= (_9B80[x] + _7B80[x])/2.0) && (hue >(_5B80[x] + _7B80[x])/2.0) && x<50 ) {if(y>49) y=49;correction =  _7B80[y] - _7B80[x] ;}
			else if (( hue<=(_5B80[x] + _7B80[x])/2.0)  && (hue > (_5B80[x]-0.035)) && x < 50) {if(y>49) y=49; correction =  _5B80[y] - _5B80[x] ; }	//

			else if( (hue >= (_15PB80[x] - 0.035)) && (hue < (_15PB80[x] + _3PB80[x])/2.0) && x<50 ) {if(y>49) y=49;correction =  _15PB80[y] - _15PB80[x] ; }
			else if (( hue>=(_45PB80[x] + _3PB80[x])/2.0)  && (hue < (_6PB80[x] + _45PB80[x])/2.0) && x < 50) {if(y>49) y=49;correction =  _45PB80[y] - _45PB80[x] ;}
			else if (( hue>=(_6PB80[x] + _45PB80[x])/2.0)  && (hue < (_6PB80[x] + _75PB80[x])/2.0) && x <50) {if(y>49) y=49;correction =  _6PB80[y] - _6PB80[x] ;}
			else if (( hue>=(_6PB80[x] + _75PB80[x])/2.0)  && (hue < (_9PB80[x] + _75PB80[x])/2.0) && x <50) {if(y>49) y=49;correction =  _75PB80[y] - _75PB80[x] ; }
			else if (( hue>=(_9PB80[x] + _75PB80[x])/2.0)  && (hue < (_9PB80[x] + 0.035)) && x <50) {if(y>49) y=49;correction =  _9PB80[y] - _9PB80[x] ; }
		}
	} // end PB correction
	if(zone==2) {//red yellow correction
		if((lum > 15.0 && lum < 25.0) ) {
			if( (hue <= (_10YR20[x] + 0.035)) && (hue > (_10YR20[x] + _85YR20[x])/2.0) && x<=45) {if(y>49) y=49;correction =  _10YR20[y] - _10YR20[x] ;}
			else if (( hue<=(_85YR20[x] + _10YR20[x])/2.0)  && (hue > (_85YR20[x] + 0.035) && x <= 45)) {if(y>49) y=49;correction =  _85YR20[y] - _85YR20[x] ;}
		}
		else if ((lum >= 25.0 && lum <35.0)) {
			if( (hue <= (_10YR30[x] + 0.035)) && (hue > (_10YR30[x] + _85YR30[x])/2.0) && x < 85) {if(y>89) y=89;correction =  _10YR30[y] - _10YR30[x] ;}
			else if( (hue <= (_10YR30[x] + _85YR30[x])/2.0) && (hue >(_85YR30[x] + _7YR30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _85YR30[y] - _85YR30[x] ;}
			else if (( hue<=(_85YR30[x] + _7YR30[x])/2.0)  && (hue > (_7YR30[x] + _55YR30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _7YR30[y] - _7YR30[x] ;}
			else if (( hue<=(_7YR30[x] + _55YR30[x])/2.0)  && (hue > (_55YR30[x] + _4YR30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _55YR30[y] - _55YR30[x] ; }
			else if (( hue<=(_55YR30[x] + _4YR30[x])/2.0)  && (hue > (_4YR30[x] + _25YR30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _4YR30[y] - _4YR30[x] ; }
			else if (( hue<=(_4YR30[x] + _25YR30[x])/2.0)  && (hue > (_25YR30[x] + _10R30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _25YR30[y] - _25YR30[x] ;}
			else if (( hue<=(_25YR30[x] + _10R30[x])/2.0)  && (hue > (_10R30[x] + _9R30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _10R30[y] - _10R30[x] ; }
			else if (( hue<=(_10R30[x] + _9R30[x])/2.0)  && (hue > (_9R30[x] + _7R30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _9R30[y] - _9R30[x] ;}
			else if (( hue<=(_9R30[x] + _7R30[x])/2.0)  && (hue > (_7R30[x] -0.035))&& x < 85) {if(y>89) y=89;correction =  _7R30[y] - _7R30[x] ; }
		}
		else if ((lum >= 35.0 && lum <45.0)) {
			if( (hue <= (_10YR40[x] + 0.035)) && (hue > (_10YR40[x] + _85YR40[x])/2.0)&& x<85) {if(y>89) y=89;correction =  _10YR40[y] - _10YR40[x] ;}
			else if( (hue <= (_10YR40[x] + _85YR40[x])/2.0) && (hue >(_85YR40[x] + _7YR40[x])/2.0)&& x < 85 ) {if(y>89) y=89;correction =  _85YR40[y] - _85YR40[x] ;}
			else if (( hue<=(_85YR40[x] + _7YR40[x])/2.0)  && (hue > (_7YR40[x] + _55YR40[x])/2.0) && x < 85) {if(y>89) y=89;correction =  _7YR40[y] - _7YR40[x] ;}
			else if (( hue<=(_7YR40[x] + _55YR40[x])/2.0)  && (hue > (_55YR40[x] + _4YR40[x])/2.0)&& x < 85 ) {if(y>89) y=89;correction =  _55YR40[y] - _55YR40[x] ; }
			else if (( hue<=(_55YR40[x] + _4YR40[x])/2.0)  && (hue > (_4YR40[x] + _25YR40[x])/2.0)&& x < 85 ) {if(y>89) y=89;correction =  _4YR40[y] - _4YR40[x] ; }
			else if (( hue<=(_4YR40[x] + _25YR40[x])/2.0)  && (hue > (_25YR40[x] + _10R40[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _25YR40[y] - _25YR40[x] ;}
			else if (( hue<=(_25YR40[x] + _10R40[x])/2.0)  && (hue > (_10R40[x] + _9R40[x])/2.0) && x < 85) {if(y>89) y=89;correction =  _10R40[y] - _10R40[x] ; }
			else if (( hue<=(_10R40[x] + _9R40[x])/2.0)  && (hue > (_9R40[x] + _7R40[x])/2.0)&& x < 85 ) {if(y>89) y=89;correction =  _9R40[y] - _9R40[x] ;}
			else if (( hue<=(_9R40[x] + _7R40[x])/2.0)  && (hue > (_7R40[x] -0.035))&& x < 85 ) {if(y>89) y=89;correction =  _7R40[y] - _7R40[x] ; }
		}
		else if ((lum >= 45.0 && lum <55.0)) {
			if( (hue <= (_10YR50[x] + 0.035)) && (hue > (_10YR50[x] + _85YR50[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _10YR50[y] - _10YR50[x] ;}
			else if( (hue <= (_10YR50[x] + _85YR50[x])/2.0) && (hue >(_85YR50[x] + _7YR50[x])/2.0)&& x < 85 ) {if(y>89) y=89;correction =  _85YR50[y] - _85YR50[x] ;}
			else if (( hue<=(_85YR50[x] + _7YR50[x])/2.0)  && (hue > (_7YR50[x] + _55YR50[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _7YR50[y] - _7YR50[x] ;}
			else if (( hue<=(_7YR50[x] + _55YR50[x])/2.0)  && (hue > (_55YR50[x] + _4YR50[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _55YR50[y] - _55YR50[x] ; }
			else if (( hue<=(_55YR50[x] + _4YR50[x])/2.0)  && (hue > (_4YR50[x] + _25YR50[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _4YR50[y] - _4YR50[x] ; }
			else if (( hue<=(_4YR50[x] + _25YR50[x])/2.0)  && (hue > (_25YR50[x] + _10R50[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _25YR50[y] - _25YR50[x] ;}
			else if (( hue<=(_25YR50[x] + _10R50[x])/2.0)  && (hue > (_10R50[x] + _9R50[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _10R50[y] - _10R50[x] ; }
			else if (( hue<=(_10R50[x] + _9R50[x])/2.0)  && (hue > (_9R50[x] + _7R50[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _9R50[y] - _9R50[x] ;}
			else if (( hue<=(_9R50[x] + _7R50[x])/2.0)  && (hue > (_7R50[x] -0.035))&& x < 85) {if(y>89) y=89;correction =  _7R50[y] - _7R50[x] ; }
		}
		else if ((lum >= 55.0 && lum <65.0)) {
			if( (hue <= (_10YR60[x] + 0.035)) && (hue > (_10YR60[x] + _85YR60[x])/2.0)) {;correction =  _10YR60[y] - _10YR60[x] ;}
			else if( (hue <= (_10YR60[x] + _85YR60[x])/2.0) && (hue >(_85YR60[x] + _7YR60[x])/2.0) ) {;correction =  _85YR60[y] - _85YR60[x] ;}
			else if (( hue<=(_85YR60[x] + _7YR60[x])/2.0)  && (hue > (_7YR60[x] + _55YR60[x])/2.0)) {correction =  _7YR60[y] - _7YR60[x] ;}
			else if (( hue<=(_7YR60[x] + _55YR60[x])/2.0)  && (hue > (_55YR60[x] + _4YR60[x])/2.0)) {correction =  _55YR60[y] - _55YR60[x] ; }
			else if (( hue<=(_55YR60[x] + _4YR60[x])/2.0)  && (hue > (_4YR60[x] + _25YR60[x])/2.0)) {correction =  _4YR60[y] - _4YR60[x] ; }
			else if (( hue<=(_4YR60[x] + _25YR60[x])/2.0)  && (hue > (_25YR60[x] + _10R60[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _25YR60[y] - _25YR60[x] ;}
			else if (( hue<=(_25YR60[x] + _10R60[x])/2.0)  && (hue > (_10R60[x] + _9R60[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _10R60[y] - _10R60[x] ; }
			else if (( hue<=(_10R60[x] + _9R60[x])/2.0)  && (hue > (_9R60[x] + _7R60[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _9R60[y] - _9R60[x] ;}
			else if (( hue<=(_9R60[x] + _7R60[x])/2.0)  && (hue > (_7R60[x] -0.035))&& x < 85) {if(y>89) y=89;correction =  _7R60[y] - _7R60[x] ; }
		}
		else if ((lum >= 65.0 && lum <75.0)) {
			if( (hue <= (_10YR70[x] + 0.035)) && (hue > (_10YR70[x] + _85YR70[x])/2.0)) {correction =  _10YR70[y] - _10YR70[x] ;}
			else if( (hue <= (_10YR70[x] + _85YR70[x])/2.0) && (hue >(_85YR70[x] + _7YR70[x])/2.0)) {correction =  _85YR70[y] - _85YR70[x] ;}
			 if (( hue<=(_85YR70[x] + _7YR70[x])/2.0)  && (hue > (_7YR70[x] + _55YR70[x])/2.0)) {correction =  _7YR70[y] - _7YR70[x] ;}
			else if (( hue<=(_7YR70[x] + _55YR70[x])/2.0)  && (hue > (_55YR70[x] + _4YR70[x])/2.0)) {correction =  _55YR70[y] - _55YR70[x] ; }
			else if (( hue<=(_55YR70[x] + _4YR70[x])/2.0)  && (hue > (_4YR70[x] + _25YR70[x])/2.0)) {correction =  _4YR70[y] - _4YR70[x] ; }
			else if (( hue<=(_4YR70[x] + _25YR70[x])/2.0)  && (hue > (_25YR70[x] + _10R70[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _25YR70[y] - _25YR70[x] ;}
			else if (( hue<=(_25YR70[x] + _10R70[x])/2.0)  && (hue > (_10R70[x] + _9R70[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _10R70[y] - _10R70[x] ; }
			else if (( hue<=(_10R70[x] + _9R70[x])/2.0)  && (hue > (_9R70[x] + _7R70[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _9R70[y] - _9R70[x] ;}
			else if (( hue<=(_9R70[x] + _7R70[x])/2.0)  && (hue > (_7R70[x] -0.035))&& x < 85) {if(y>89) y=89;correction =  _7R70[y] - _7R70[x] ; }
		}
		else if ((lum >= 75.0 && lum <85.0)) {
			if( (hue <= (_10YR80[x] + 0.035)) && (hue > (_10YR80[x] + _85YR80[x])/2.0)) {correction =  _10YR80[y] - _10YR80[x] ;}
			else if( (hue <= (_10YR80[x] + _85YR80[x])/2.0) && (hue >(_85YR80[x] + _7YR80[x])/2.0)) {correction =  _85YR80[y] - _85YR80[x] ;}
			else if (( hue<=(_85YR80[x] + _7YR80[x])/2.0)  && (hue > (_7YR80[x] + _55YR80[x])/2.0) && x<85) {if(y>89) y=89;correction =  _7YR80[y] - _7YR80[x] ;}
			else if (( hue<=(_7YR80[x] + _55YR80[x])/2.0)  && (hue > (_55YR80[x] + _4YR80[x])/2.0) && x <45) {correction =  _55YR80[y] - _55YR80[x] ; }
			else if (( hue<=(_55YR80[x] + _4YR80[x])/2.0)  && (hue > (_4YR80[x] - 0.035) && x<45)) {if(y>49) y=49;correction =  _4YR80[y] - _4YR80[x] ; }
		}
		else if ((lum >= 85.0 && lum <95.0)) {
			if( (hue <= (_10YR90[x] + 0.035)) && (hue > (_10YR90[x] -0.035) && x<85)) {if(y>89) y=89;correction =  _10YR90[y] - _10YR90[x] ;}
			else if ( hue<=(_85YR90[x] + 0.035)  && hue > (_85YR90[x] -0.035) && x<85) {if(y>89) y=89;correction =  _85YR90[y] - _85YR90[x] ;}
			else if (( hue<=(_55YR90[x] + 0.035)  && (hue > (_55YR90[x] - 0.035) && x<45))) {if(y>49) y=49;correction =  _55YR90[y] - _55YR90[x] ; }
		}//end red yellow
	}
	if(zone==3) {//Green yellow correction
		if ((lum >= 25.0 && lum <35.0)) {
			if( (hue <= (_7G30[x] + 0.035)) && (hue > (_7G30[x] + _5G30[x])/2.0) ) {correction =  _7G30[y] - _7G30[x] ;}
			else if( (hue <= (_7G30[x] + _5G30[x])/2.0) && (hue >(_5G30[x] + _25G30[x])/2.0)) {correction =  _5G30[y] - _5G30[x] ;}
			else if (( hue<=(_25G30[x] + _5G30[x])/2.0)  && (hue > (_25G30[x] + _1G30[x])/2.0)) {correction =  _25G30[y] - _25G30[x] ;}
			else if (( hue<=(_1G30[x] + _25G30[x])/2.0)  && (hue > (_1G30[x] + _10GY30[x])/2.0)) {correction =  _1G30[y] - _1G30[x] ; }
			else if (( hue<=(_1G30[x] + _10GY30[x])/2.0)  && (hue > (_10GY30[x] + _75GY30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _10GY30[y] - _10GY30[x] ; }
			else if (( hue<=(_10GY30[x] + _75GY30[x])/2.0)  && (hue > (_75GY30[x] + _5GY30[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _75GY30[y] - _75GY30[x] ;}
			else if (( hue<=(_5GY30[x] + _75GY30[x])/2.0)  && (hue > (_5GY30[x] -0.035))&& x < 85) {if(y>89) y=89;correction =  _5GY30[y] - _5GY30[x] ; }
		}
		 else if ((lum >= 35.0 && lum <45.0)) {
			if( (hue <= (_7G40[x] + 0.035)) && (hue > (_7G40[x] + _5G40[x])/2.0) ) {correction =  _7G40[y] - _7G40[x] ;}
			else if( (hue <= (_7G40[x] + _5G40[x])/2.0) && (hue >(_5G40[x] + _25G40[x])/2.0)) {correction =  _5G40[y] - _5G40[x] ;}
			else if (( hue<=(_25G40[x] + _5G40[x])/2.0)  && (hue > (_25G40[x] + _1G40[x])/2.0)) {correction =  _25G40[y] - _25G40[x] ;}
			else if (( hue<=(_1G40[x] + _25G40[x])/2.0)  && (hue > (_1G40[x] + _10GY40[x])/2.0)) {correction =  _1G40[y] - _1G40[x] ; }
			else if (( hue<=(_1G40[x] + _10GY40[x])/2.0)  && (hue > (_10GY40[x] + _75GY40[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _10GY40[y] - _10GY40[x] ; }
			else if (( hue<=(_10GY40[x] + _75GY40[x])/2.0)  && (hue > (_75GY40[x] + _5GY40[x])/2.0)&& x < 85) {if(y>89) y=89;correction =  _75GY40[y] - _75GY40[x] ;}
			else if (( hue<=(_5GY40[x] + _75GY40[x])/2.0)  && (hue > (_5GY40[x]-0.035)) && x < 85) {if(y>89) y=89; correction =  _5GY40[y] - _5GY40[x] ; }	//
		}
		 else if ((lum >= 45.0 && lum <55.0)) {
			if( (hue <= (_7G50[x] + 0.035)) && (hue > (_7G50[x] + _5G50[x])/2.0) ) {correction =  _7G50[y] - _7G50[x] ;}
			else if( (hue <= (_7G50[x] + _5G50[x])/2.0) && (hue >(_5G50[x] + _25G50[x])/2.0)) {correction =  _5G50[y] - _5G50[x] ;}
			else if (( hue<=(_25G50[x] + _5G50[x])/2.0)  && (hue > (_25G50[x] + _1G50[x])/2.0)) {correction =  _25G50[y] - _25G50[x] ;}
			else if (( hue<=(_1G50[x] + _25G50[x])/2.0)  && (hue > (_1G50[x] + _10GY50[x])/2.0)) {correction =  _1G50[y] - _1G50[x] ; }
			else if (( hue<=(_1G50[x] + _10GY50[x])/2.0)  && (hue > (_10GY50[x] + _75GY50[x])/2.0)) {correction =  _10GY50[y] - _10GY50[x] ; }
			else if (( hue<=(_10GY50[x] + _75GY50[x])/2.0)  && (hue > (_75GY50[x] + _5GY50[x])/2.0)) {correction =  _75GY50[y] - _75GY50[x] ;}
			else if (( hue<=(_5GY50[x] + _75GY50[x])/2.0)  && (hue > (_5GY50[x] -0.035))) {correction =  _5GY50[y] - _5GY50[x] ; }
		}
		 else if ((lum >= 55.0 && lum <65.0)) {
			if( (hue <= (_7G60[x] + 0.035)) && (hue > (_7G60[x] + _5G60[x])/2.0) ) {correction =  _7G60[y] - _7G60[x] ;}
			else if( (hue <= (_7G60[x] + _5G60[x])/2.0) && (hue >(_5G60[x] + _25G60[x])/2.0)) {correction =  _5G60[y] - _5G60[x] ;}
			else if (( hue<=(_25G60[x] + _5G60[x])/2.0)  && (hue > (_25G60[x] + _1G60[x])/2.0)) {correction =  _25G60[y] - _25G60[x] ;}
			else if (( hue<=(_1G60[x] + _25G60[x])/2.0)  && (hue > (_1G60[x] + _10GY60[x])/2.0)) {correction =  _1G60[y] - _1G60[x] ; }
			else if (( hue<=(_1G60[x] + _10GY60[x])/2.0)  && (hue > (_10GY60[x] + _75GY60[x])/2.0)) {correction =  _10GY60[y] - _10GY60[x] ; }
			else if (( hue<=(_10GY60[x] + _75GY60[x])/2.0)  && (hue > (_75GY60[x] + _5GY60[x])/2.0)) {correction =  _75GY60[y] - _75GY60[x] ;}
			else if (( hue<=(_5GY60[x] + _75GY60[x])/2.0)  && (hue > (_5GY60[x] -0.035))) {correction =  _5GY60[y] - _5GY60[x] ; }
		}
		 else if ((lum >= 65.0 && lum <75.0)) {
			if( (hue <= (_7G70[x] + 0.035)) && (hue > (_7G70[x] + _5G70[x])/2.0) ) {correction =  _7G70[y] - _7G70[x] ;}
			else if( (hue <= (_7G70[x] + _5G70[x])/2.0) && (hue >(_5G70[x] + _25G70[x])/2.0)) {correction =  _5G70[y] - _5G70[x] ;}
			else if (( hue<=(_25G70[x] + _5G70[x])/2.0)  && (hue > (_25G70[x] + _1G70[x])/2.0)) {correction =  _25G70[y] - _25G70[x] ;}
			else if (( hue<=(_1G70[x] + _25G70[x])/2.0)  && (hue > (_1G70[x] + _10GY70[x])/2.0)) {correction =  _1G70[y] - _1G70[x] ; }
			else if (( hue<=(_1G70[x] + _10GY70[x])/2.0)  && (hue > (_10GY70[x] + _75GY70[x])/2.0)) {correction =  _10GY70[y] - _10GY70[x] ; }
			else if (( hue<=(_10GY70[x] + _75GY70[x])/2.0)  && (hue > (_75GY70[x] + _5GY70[x])/2.0)) {correction =  _75GY70[y] - _75GY70[x] ;}
			else if (( hue<=(_5GY70[x] + _75GY70[x])/2.0)  && (hue > (_5GY70[x] -0.035))) {correction =  _5GY70[y] - _5GY70[x] ; }
		}
		 else if ((lum >= 75.0 && lum <85.0)) {
			if( (hue <= (_7G80[x] + 0.035)) && (hue > (_7G80[x] + _5G80[x])/2.0) ) {correction =  _7G80[y] - _7G80[x] ;}
			else if( (hue <= (_7G80[x] + _5G80[x])/2.0) && (hue >(_5G80[x] + _25G80[x])/2.0)) {correction =  _5G80[y] - _5G80[x] ;}
			else if (( hue<=(_25G80[x] + _5G80[x])/2.0)  && (hue > (_25G80[x] + _1G80[x])/2.0)) {correction =  _25G80[y] - _25G80[x] ;}
			else if (( hue<=(_1G80[x] + _25G80[x])/2.0)  && (hue > (_1G80[x] + _10GY80[x])/2.0)) {correction =  _1G80[y] - _1G80[x] ; }
			else if (( hue<=(_1G80[x] + _10GY80[x])/2.0)  && (hue > (_10GY80[x] + _75GY80[x])/2.0)) {correction =  _10GY80[y] - _10GY80[x] ; }
			else if (( hue<=(_10GY80[x] + _75GY80[x])/2.0)  && (hue > (_75GY80[x] + _5GY80[x])/2.0)) {correction =  _75GY80[y] - _75GY80[x] ;}
			else if (( hue<=(_5GY80[x] + _75GY80[x])/2.0)  && (hue > (_5GY80[x] -0.035))) {correction =  _5GY80[y] - _5GY80[x] ; }
		 }//end green yellow

	}

	if(zone==4) {//Red purple correction : only for L < 30
		if ((lum > 5.0 && lum < 15.0)) {
			if( (hue <= (_5R10[x] + 0.035)) && (hue > (_5R10[x] - 0.043)) && x<45) {if(y>44) y=44;correction =  _5R10[y] - _5R10[x] ;}
			else if( (hue <= (_25R10[x] + 0.043)) && (hue >(_25R10[x] + _10RP10[x])/2.0) && x<45 ) {if(y>44) y=44;correction =  _25R10[y] - _25R10[x] ;}
			else if ( (hue <=(_25R10[x] + _10RP10[x])/2.0) && (hue > (_10RP10[x] -0.035) ) && x<45){if(y>44) y=44; correction =  _10RP10[y] - _10RP10[x] ;}
		}
		else if ((lum >= 15.0 && lum <25.0)) {
			if( (hue <= (_5R20[x] + 0.035)) && (hue > (_5R20[x] + _25R20[x])/2.0) && x<70 ) {if(y>70) y=70;correction =  _5R20[y] - _5R20[x] ;}
			else if( (hue <= (_5R20[x] + _25R20[x])/2.0) && (hue >(_10RP20[x] + _25R20[x])/2.0) && x<70) {if(y>70) y=70;correction =  _25R20[y] - _25R20[x] ;}
			else if (( hue<=(_10RP20[x] + _25R20[x])/2.0)  && (hue > (_10RP20[x] -0.035)) && x<70) {if(y>70) y=70; correction =  _10RP20[y] - _10RP20[x] ;}
		}
		else if ((lum >= 25.0 && lum <35.0)) {
			if( (hue <= (_5R30[x] + 0.035)) && (hue > (_5R30[x] + _25R30[x])/2.0) && x<85 ) {if(y>85) y=85;correction =  _5R30[y] - _5R30[x] ;}
			else if( (hue <= (_5R30[x] + _25R30[x])/2.0) && (hue >(_10RP30[x] + _25R30[x])/2.0) && x< 85) {if(y>85) y=85;correction =  _25R30[y] - _25R30[x] ;}
			else if (( hue<=(_10RP30[x] + _25R30[x])/2.0)  && (hue > (_10RP30[x] -0.035)) && x<85) {if(y>85) y=85; correction =  _10RP30[y] - _10RP30[x] ;}
		}//end red purple
	}
}


/*
 * copyright (c)2011  Jacques Desmis <jdesmis@gmail.com>
 * skin color: mixed from NX2 skin color palette, Von Luschan, and photos of people white,
 * black, yellow....there are some little exceptions...cover 99% case
 * pay attention to white balance, and do not change hue and saturation, upstream of the modification
 *
 */
void ImProcFunctions::skinsat (float lum, float hue, float chrom, float &satreduc) {

	float reduction=0.3;// to be adapted...by tests
	float extendedreduction=0.4;
	float extendedreduction2=0.6;

	float C9=0.0, C8=0.0, C7=0.0, C4=0.0, C3=0.0, C2=0.0, C1=0.0;
	float H9=0.0, H8=0.0, H7=0.0, H4=0.0, H3=0.0, H2=0.0, H1=0.0, H10=0.0,H11=0.0;
	H9=0.05;H8=0.25;H7=0.1;H4=0.02;H3=0.02;H2=0.1;H1=0.1;H10=-0.2;H11=-0.2;//H10 and H11 are curious...H11=-0.8 ??
	C9=8.0;C8=15.0;C7=12.0;C4=7.0;C3=5.0;C2=5.0;C1=5.0;
	// wide area for transition
	if(lum >= 92.0 && (hue > -0.1 && hue < 1.65) && (chrom > 7.0 && chrom < (18.0))) satreduc=extendedreduction2;
	else if (lum >= 85.0 && lum < 92.0 && (hue > 0.0 && hue < 1.65) && (chrom > 7.0 && chrom < (35.0+C9))) satreduc=extendedreduction2;
	else if ((lum > 20 && lum < 85) && (hue > (0.02 + H11) && hue < 1.65) && (chrom > 7.0 && chrom < (55.0+C9) )) satreduc=extendedreduction2;
	else if (lum < 20.0 && (hue > (0.02+H11) && hue < 1.60) && (chrom > 7.0 && chrom < (45.0+C1) )) satreduc=extendedreduction2;

	// wide area  skin color, useful if not accurate colorimetry or if the user has changed hue and saturation

	if(lum >= 92.0  && (hue > 0.8 && hue < 1.65) && (chrom > 7.0 && chrom < (15.0))) satreduc=extendedreduction;
	else if(lum >= 85.0 && lum < 92.0  && (hue > 0.70 && hue < 1.4) && (chrom > 7.0 && chrom < (26.0+C9))) satreduc=extendedreduction;
	else if ((lum > 20 && lum < 85) && (hue > (0.02 + H11) && hue < 1.5) && (chrom > 7.0 && chrom < (48.0+C9) )) satreduc=extendedreduction;
	else if (lum < 20.0  && (hue > (0.02+H11) && hue < 1.0) && (chrom > 7.0 && chrom < (35.0+C1) )) satreduc=extendedreduction;

	// "real" skin color : take into account a slightly usage of contrast and saturation in RT if option "skin" = 1
	if(lum >= 85.0  && (hue > (0.78-H9) && hue < (1.18+H9)) && (chrom > 8.0 && chrom < (14.0+C9))) satreduc=reduction;
	else if ((lum >= 70.0 && lum < 85.0)  && (hue > 0.4 && hue < (1.04+H8)) && (chrom > 8.0 && chrom < (35.0+C8))) satreduc=reduction;
	else if ((lum >= 52.0 && lum < 70.0)  && (hue > 0.3 && hue < (1.27+H7)) && (chrom > 11.0 && chrom < (35.0+C7))) satreduc=reduction;
	else if ((lum >= 35.0 && lum < 52.0)  && (hue > 0.3 && hue < (1.25+H4)) && (chrom > 13.0 && chrom < (37.0+C4))) satreduc=reduction;
	else if ((lum >= 20.0 && lum < 35.0)  && (hue > 0.3 && hue < (1.20+H3)) && (chrom > 7.0 && chrom <(35.0+C3) )) satreduc=reduction;
	else if ((lum > 10.0 && lum < 20.0)  && (hue > (0.0 + H10) && hue < (0.95 +H2)) && (chrom > 8.0 && chrom < (23.0+C2))) satreduc=reduction;
	else if ((lum < 10.0)  && (hue > (0.02 + H10) && hue < (0.90+H1)) && (chrom > 8.0 && chrom < (23.0+C1))) satreduc=reduction; // no data : extrapolate
}

/*
 * vibrance correction
 * copyright (c)2011  Jacques Desmis <jdesmis@gmail.com> and Jean-Christophe Frisch <natureh@free.fr>
 *
 */
void ImProcFunctions::vibrance (LabImage* lab) {
	if (!params->vibrance.enabled || (!params->vibrance.pastels && !params->vibrance.saturated))
		return;

	int width = lab->W, height = lab->H;

#ifdef _DEBUG
	MyTime t1e,t2e;
	t1e.set();

	float maxdeltaHueBP=0.0,maxdeltaHueRY=0.0,maxdeltaHueGY=0.0, maxdeltaHueRP=0.0;
	int negat=0, moreRGB=0, negsat=0 ,moresat=0;
	int Munspb=0, Munsry=0, Munsgy=0, Munsrp=0;
	int depass=0;
#endif

	/*float *LL,*CC,*HH;

	LL = new float[width*height];//pointer for ulterior usage(convolution...)
	CC = new float[width*height];
	HH = new float[width*height];
	*/

#ifdef _DEBUG
#pragma omp parallel default(shared) reduction(+: negat, moreRGB, negsat ,moresat, Munspb, Munsry, Munsgy, Munsrp, depass) if (multiThread)
#else
#pragma omp parallel default(shared) if (multiThread)
#endif
{

	float LL,CC,HH;
	float R,G,B,RR,GG,BB;
	float fy,fx,fz,x_,y_,z_,Lprov,Lprov1,aprov1,bprov1,aprovn,bprovn,fxx,fyy,fzz,xx_,yy_,zz_;
	float saturation;
	TMatrix wiprof = iccStore->workingSpaceInverseMatrix (params->icm.working);
	float chromaPastel= (float) params->vibrance.pastels   / 100.0f;//
	float chromaSatur = (float) params->vibrance.saturated / 100.0f;//
	bool highlight = params->hlrecovery.enabled;//Get the value if "highlight reconstruction" is activated
	//inverse matrix user select
	double wip[3][3] = {
		{wiprof[0][0],wiprof[0][1],wiprof[0][2]},
		{wiprof[1][0],wiprof[1][1],wiprof[1][2]},
		{wiprof[2][0],wiprof[2][1],wiprof[2][2]}
	};
	float Chprov,memChprov,Chprov1;
	float satredu;//reduct sat in function of skin
	float sathue[5],sathue2[4];// adjust sat in function of hue
	float correctionHue; // Munsell's correction
	float limitpastelsatur;
	int zone=0;
	bool allwaysingamut=true;
	limitpastelsatur=(float)params->vibrance.psthreshold / 100.0f;
	if (limitpastelsatur < 0.07) limitpastelsatur=0.07;
	float p0,p1,p2;//adapt limit of pyramid to psThreshold
	float s0,s1,s2;
	float maxdp=(limitpastelsatur-0.07)/4.0;
	float maxds=(1.0-limitpastelsatur)/4.0;
	p0=0.07+maxdp;
	p1=0.07+2.0*maxdp;
	p2=0.07+3.0*maxdp;
	s0=limitpastelsatur + maxds;
	s1=limitpastelsatur + 2.0*maxds;
	s2=limitpastelsatur + 3.0*maxds;

	//if (settings->verbose) printf("vibrance:  p0=%1.2f  p1=%1.2f  p2=%1.2f  s0=%1.2f s1=%1.2f s2=%1.2f\n", p0,p1,p2,s0,s1,s2);
	if (settings->verbose) printf("vibrance:  pastel=%f   satur=%f   limit= %1.2f\n",1.0+chromaPastel,1.0+chromaSatur, limitpastelsatur);

#pragma omp for schedule(dynamic, 10)
	for (int i=0; i<height; i++)
		for (int j=0; j<width; j++) {
			//int pos = i*width+j;
			LL=lab->L[i][j]/327.68f;
			CC=sqrt(lab->a[i][j]/327.68f*lab->a[i][j]/327.68f + lab->b[i][j]/327.68f*lab->b[i][j]/327.68f);
			HH=atan2(lab->b[i][j],lab->a[i][j]);
			//double pyramid: LL and HH
			//I try to take into account: Munsell response (human vision) and Gamut..(less response for red): preferably using Prophoto or WideGamut
			//blue: -1.80 -3.14  green = 2.1 3.14   green-yellow=1.4 2.1  red:0 1.4  blue-purple:-0.7  -1.4   purple: 0 -0.7
			//these values allow a better and differential response
			if(LL < 20.0) {//more for blue-purple, blue and red modulate
				if     (HH< -1.5 && HH>- 3.1415) {sathue[0]=1.3;sathue[1]=1.2;sathue[2]=1.1;sathue[3]=1.05;sathue[4]=0.4;sathue2[0]=1.05;sathue2[1]=1.1 ;sathue2[2]=1.05;sathue2[3]=1.0;}//blue
				else if(HH>  2.1 && HH<= 3.1415) {sathue[0]=1.4;sathue[1]=1.3;sathue[2]=1.2;sathue[3]=1.15;sathue[4]=0.4;sathue2[0]=1.15;sathue2[1]=1.1 ;sathue2[2]=1.05;sathue2[3]=1.0;}//green
				else if(HH>  1.4 && HH<= 2.1   ) {sathue[0]=1.0;sathue[1]=1.0;sathue[2]=1.0;sathue[3]=1.0 ;sathue[4]=0.4;sathue2[0]=1.0 ;sathue2[1]=1.0 ;sathue2[2]=1.0 ;sathue2[3]=1.0;}//green yellow 1.2 1.1
				else if(HH< -0.7 && HH>=-1.5   ) {sathue[0]=1.6;sathue[1]=1.4;sathue[2]=1.3;sathue[3]=1.2 ;sathue[4]=0.4;sathue2[0]=1.2 ;sathue2[1]=1.15;sathue2[2]=1.1 ;sathue2[3]=1.0;}//blue purple  1.2 1.1
	//			else if(HH>= 0.0 && HH<= 1.4   ) {sathue[0]=1.1;sathue[1]=1.1;sathue[2]=1.1;sathue[3]=1.0 ;sathue[4]=0.4;sathue2[0]=1.0 ;sathue2[1]=1.0 ;sathue2[2]=1.0 ;sathue2[3]=1.0;}//red   0.8 0.7
				else if(HH>= 0.0 && HH<= 1.4   ) {sathue[0]=1.3;sathue[1]=1.2;sathue[2]=1.1;sathue[3]=1.0 ;sathue[4]=0.4;sathue2[0]=1.0 ;sathue2[1]=1.0 ;sathue2[2]=1.0 ;sathue2[3]=1.0;}//red   0.8 0.7
				else if(HH<  0.0 && HH>=-0.7   ) {sathue[0]=1.2;sathue[1]=1.0;sathue[2]=1.0;sathue[3]=1.0 ;sathue[4]=0.4;sathue2[0]=1.0 ;sathue2[1]=1.0 ;sathue2[2]=1.0 ;sathue2[3]=1.0;}//purple
			}
			else if (LL>=20 && LL< 50) {//more for blue and green, less for red and green-yellow
				if     (HH< -1.5 && HH>- 3.1415) {sathue[0]=1.5;sathue[1]=1.4;sathue[2]=1.3;sathue[3]=1.2 ;sathue[4]=0.4;sathue2[0]=1.2 ;sathue2[1]=1.1 ;sathue2[2]=1.05;sathue2[3]=1.0;}//blue
				else if(HH>  2.1 && HH<= 3.1415) {sathue[0]=1.5;sathue[1]=1.4;sathue[2]=1.3;sathue[3]=1.2 ;sathue[4]=0.4;sathue2[0]=1.2 ;sathue2[1]=1.1 ;sathue2[2]=1.05;sathue2[3]=1.0;}//green
				else if(HH>  1.4 && HH<= 2.1   ) {sathue[0]=1.1;sathue[1]=1.1;sathue[2]=1.1;sathue[3]=1.05;sathue[4]=0.4;sathue2[0]=0.9 ;sathue2[1]=0.8 ;sathue2[2]=0.7 ;sathue2[3]=0.6;}//green yellow 1.2 1.1
				else if(HH< -0.7 && HH>=-1.5   ) {sathue[0]=1.3;sathue[1]=1.2;sathue[2]=1.1;sathue[3]=1.05;sathue[4]=0.4;sathue2[0]=1.05;sathue2[1]=1.05;sathue2[2]=1.0 ;sathue2[3]=1.0;}//blue purple  1.2 1.1
				else if(HH<  0.0 && HH>=-0.7   ) {sathue[0]=1.2;sathue[1]=1.0;sathue[2]=1.0;sathue[3]=1.0 ;sathue[4]=0.4;sathue2[0]=1.0 ;sathue2[1]=1.0 ;sathue2[2]=1.0 ;sathue2[3]=1.0;}//purple
	//			else if(HH>= 0.0 && HH<= 1.4   ) {sathue[0]=0.8;sathue[1]=0.8;sathue[2]=0.8;sathue[3]=0.8 ;sathue[4]=0.4;sathue2[0]=0.8 ;sathue2[1]=0.8 ;sathue2[2]=0.8 ;sathue2[3]=0.8;}//red   0.8 0.7
				else if(HH>= 0.0 && HH<= 1.4   ) {sathue[0]=1.1;sathue[1]=1.0;sathue[2]=0.9;sathue[3]=0.8 ;sathue[4]=0.4;sathue2[0]=0.8 ;sathue2[1]=0.8 ;sathue2[2]=0.8 ;sathue2[3]=0.8;}//red   0.8 0.7

			}
			else if (LL>=50 && LL< 80) {//more for green, less for red and green-yellow
				if     (HH< -1.5 && HH>- 3.1415) {sathue[0]=1.3;sathue[1]=1.2;sathue[2]=1.15;sathue[3]=1.1 ;sathue[4]=0.3;sathue2[0]=1.1 ;sathue2[1]=1.1 ;sathue2[2]=1.05;sathue2[3]=1.0;}//blue
				else if(HH>  2.1 && HH<= 3.1415) {sathue[0]=1.6;sathue[1]=1.4;sathue[2]=1.3 ;sathue[3]=1.25;sathue[4]=0.3;sathue2[0]=1.25;sathue2[1]=1.2 ;sathue2[2]=1.15;sathue2[3]=1.05;}//green - even with Prophoto green are too "little"  1.5 1.3
				else if(HH>  1.4 && HH<= 2.1   ) {sathue[0]=1.3;sathue[1]=1.2;sathue[2]=1.1 ;sathue[3]=1.05;sathue[4]=0.3;sathue2[0]=1.0 ;sathue2[1]=0.9 ;sathue2[2]=0.8 ;sathue2[3]=0.7;}//green yellow 1.2 1.1
				else if(HH< -0.7 && HH>=-1.5   ) {sathue[0]=1.3;sathue[1]=1.2;sathue[2]=1.15;sathue[3]=1.1 ;sathue[4]=0.3;sathue2[0]=1.1 ;sathue2[1]=1.05;sathue2[2]=1.0 ;sathue2[3]=1.0;}//blue purple  1.2 1.1
				else if(HH<  0.0 && HH>=-0.7   ) {sathue[0]=1.2;sathue[1]=1.0;sathue[2]=1.0 ;sathue[3]=1.0 ;sathue[4]=0.3;sathue2[0]=1.0 ;sathue2[1]=1.0 ;sathue2[2]=1.0 ;sathue2[3]=1.0;}//purple
	//			else if(HH>= 0.0 && HH<= 1.4   ) {sathue[0]=0.8;sathue[1]=0.8;sathue[2]=0.8 ;sathue[3]=0.8 ;sathue[4]=0.3;sathue2[0]=0.8 ;sathue2[1]=0.8 ;sathue2[2]=0.8 ;sathue2[3]=0.8;}//red   0.8 0.7
				else if(HH>= 0.0 && HH<= 1.4   ) {sathue[0]=1.1;sathue[1]=1.0;sathue[2]=0.9 ;sathue[3]=0.8 ;sathue[4]=0.3;sathue2[0]=0.8 ;sathue2[1]=0.8 ;sathue2[2]=0.8 ;sathue2[3]=0.8;}//red   0.8 0.7
			}
			else if (LL>=80) {//more for green-yellow, less for red and purple
				if     (HH< -1.5 && HH>- 3.1415) {sathue[0]=1.0;sathue[1]=1.0;sathue[2]=0.9;sathue[3]=0.8;sathue[4]=0.2;sathue2[0]=0.8;sathue2[1]=0.8 ;sathue2[2]=0.8 ;sathue2[3]=0.8;}//blue
				else if(HH>  2.1 && HH<= 3.1415) {sathue[0]=1.4;sathue[1]=1.3;sathue[2]=1.2;sathue[3]=1.1;sathue[4]=0.2;sathue2[0]=1.1;sathue2[1]=1.05;sathue2[2]=1.05;sathue2[3]=1.0;}//green
				else if(HH>  1.4 && HH<= 2.1   ) {sathue[0]=1.6;sathue[1]=1.5;sathue[2]=1.4;sathue[3]=1.2;sathue[4]=0.2;sathue2[0]=1.1;sathue2[1]=1.05;sathue2[2]=1.0 ;sathue2[3]=1.0;}//green yellow 1.2 1.1
				else if(HH< -0.7 && HH>=-1.5   ) {sathue[0]=1.0;sathue[1]=1.0;sathue[2]=0.9;sathue[3]=0.8;sathue[4]=0.2;sathue2[0]=0.8;sathue2[1]=0.8 ;sathue2[2]=0.8 ;sathue2[3]=0.8;}//blue purple  1.2 1.1
				else if(HH<  0.0 && HH>=-0.7   ) {sathue[0]=1.2;sathue[1]=1.0;sathue[2]=1.0;sathue[3]=0.9;sathue[4]=0.2;sathue2[0]=0.9;sathue2[1]=0.9 ;sathue2[2]=0.8 ;sathue2[3]=0.8;}//purple
	//			else if(HH>= 0.0 && HH<= 1.4   ) {sathue[0]=0.8;sathue[1]=0.8;sathue[2]=0.8;sathue[3]=0.8;sathue[4]=0.2;sathue2[0]=0.8;sathue2[1]=0.8 ;sathue2[2]=0.8 ;sathue2[3]=0.8;}//red   0.8 0.7
				else if(HH>= 0.0 && HH<= 1.4   ) {sathue[0]=1.1;sathue[1]=1.0;sathue[2]=0.9;sathue[3]=0.8;sathue[4]=0.2;sathue2[0]=0.8;sathue2[1]=0.8 ;sathue2[2]=0.8 ;sathue2[3]=0.8;}//red   0.8 0.7
			}

			satredu=1.0;
			if(params->vibrance.protectskins)
				skinsat (LL, HH, CC, satredu);// for skin colors
			// here we work on Chromaticity and Hue
			// variation of Chromaticity  ==> saturation via RGB
			// Munsell correction
			// then conversion to Lab
			Chprov=Chprov1=CC;
			memChprov=Chprov;
			Lprov1=LL;
			bool inGamut;
			//begin gamut control: specially if no ICC profil is used and if ICC and GamutICC=true
			do {
				inGamut=true;

				//Lprov1=LL;
				aprov1=Chprov1*cos(HH);
				bprov1=Chprov1*sin(HH);

				//conversion Lab RGB to limit Lab values - this conversion is usefull before Munsell correction
				fy = (0.00862069 *Lprov1 )+ 0.137932;
				fx = (0.002 * aprov1) + fy;
				fz = fy - (0.005 * bprov1);

				x_ = 65535.0 * f2xyz(fx)*D50x;
				y_ = 65535.0 * f2xyz(fy);
				z_ = 65535.0 * f2xyz(fz)*D50z;
				xyz2rgb(x_,y_,z_,R,G,B,wip);

				// gamut control before saturation to put Lab values in future gamut, but not RGB
				if (allwaysingamut) {        // gamut control
					if (R<0.0f || G<0.0f || B<0.0f) {
#ifdef _DEBUG
						negat++;
#endif
						if (Lprov1 < 0.01f)
							 Lprov1=0.05f;
						Chprov1 *=0.95f;     // one can modify this value...
						if (Chprov1> 3.0f)
							inGamut=false;
						else {
							Lprov1+=0.3f;
							inGamut=false;
						}
					}
					else if ((!highlight) && (R>65534.5f || G>65534.5f || B>65534.5f)){ // if "highlight reconstruction" is enabled, don't control Gamut
#ifdef _DEBUG
						moreRGB++;
#endif
						if (Lprov1>99.99f)
							Lprov1=99.8f;
						Chprov1 *=0.95f;
						if (Chprov1> 3.0f)
							inGamut=false;
						else {
							Lprov1 -=0.3f;
							inGamut=false;
						}
					}
				}
			}
			while (!inGamut);
			//end first gamut control

			Chprov=Chprov1;
			Lprov=Lprov1;

			saturation=SAT(R,G,B);
			// work on saturation
			if(Chprov > 6.0) {  //protect gray and LUT Munsell
				//pyramid to adjust saturation in function of saturation and hue (and Luminance)
				if(satredu!=1.0){
					// for skin, no differenciation
					sathue [0]=1.; sathue [1]=1.; sathue [2]=1.; sathue [3]=1.0; sathue[4]=1.0;
					sathue2[0]=1.; sathue2[1]=1.; sathue2[2]=1.; sathue2[3]=1.0;
				}

				if(saturation>0.0) {
					float chmodpastel,chmodsat;

					// We handle only positive values here
					if      (saturation < 0.07)             chmodpastel = chromaPastel*satredu*sathue[4];   //neutral tones
					else if (saturation < p0)               chmodpastel = chromaPastel*satredu*sathue[0];
					else if (saturation < p1)               chmodpastel = chromaPastel*satredu*sathue[1];
					else if (saturation < p2)               chmodpastel = chromaPastel*satredu*sathue[2];
					else if (saturation < limitpastelsatur) chmodpastel = chromaPastel*satredu*sathue[3];
					else if (saturation < s0)               chmodsat    = chromaSatur*satredu*sathue2[0];
					else if (saturation < s1)               chmodsat    = chromaSatur*satredu*sathue2[1];
					else if (saturation < s2)               chmodsat    = chromaSatur*satredu*sathue2[2];
					else                                    chmodsat    = chromaSatur*satredu*sathue2[3];

					if(chromaPastel != chromaSatur){
						//if sliders pastels and saturated differents: tansition with linear interpolation between p2 and s0
						float chromaPastel_a, chromaPastel_b, chromamean, chromaSatur_a, chromaSatur_b, newchromaPastel, newchromaSatur;

						chromamean = (chromaSatur + chromaPastel)/2.0;
						chromaPastel_a = (chromaPastel-chromamean)/(p2-limitpastelsatur);
						chromaPastel_b = chromaPastel-chromaPastel_a*p2;
						if(saturation > p2 && saturation < limitpastelsatur) {
							newchromaPastel = chromaPastel_a*saturation + chromaPastel_b;
							chmodpastel = newchromaPastel*satredu*sathue[3];
						}
						
						chromaSatur_a=(chromaSatur-chromamean)/(s0-limitpastelsatur);
						chromaSatur_b=chromaSatur-chromaSatur_a*s0;
						if(saturation < s0 && saturation >=limitpastelsatur) {newchromaSatur=chromaSatur_a*saturation + chromaSatur_b; chmodsat = newchromaSatur*satredu*sathue2[0];}
								}// end transition		
					if (saturation <= limitpastelsatur) {
						if(chmodpastel >  2.0 ) chmodpastel = 2.0;   //avoid too big values
						if(chmodpastel < -0.93) chmodpastel =-0.93;  //avoid negative values

						Chprov *=(1.0+chmodpastel);
						if(Chprov<6.0) Chprov=6.0;
					}
					else { //if (saturation > limitpastelsatur)
						if(chmodsat >  1.8 ) chmodsat = 1.8;        //saturated
						if(chmodsat < -0.93) chmodsat =-0.93;

						Chprov *= 1.0+chmodsat;
						if(Chprov <6.0) Chprov=6.0;
					}
				}
			}

			//Munsell correction
			correctionHue=0.0;
			do {
				inGamut=true;
				if(params->vibrance.avoidcolorshift) {
					if(memChprov >= 6.0 && memChprov < 140) {          //if C > 140 we say C=140 (only in Prophoto ...with very large saturation)
						if (Chprov > 140) Chprov=139;                  //limits of LUTf
						if(HH > -2.48 && HH < -0.55) {       //limits of hue Blue purple  -1.90
							zone=1; // *** blue purple correction and blue correction for sky
							MunsellLch (Lprov, HH,Chprov, memChprov, correctionHue, zone);
#ifdef _DEBUG
							if(correctionHue !=0.0) {
								if(fabs(correctionHue) > maxdeltaHueBP)
									maxdeltaHueBP=fabs(correctionHue);
								Munspb++;
							}
							if(fabs(correctionHue) > 0.45) depass++;   //verify if no bug in calculation
#endif
						}
						if(HH > 0.44 && HH < 1.52) {         //limits of hue red yellow
							zone=2; // *** red yellow correction
							MunsellLch (Lprov, HH,Chprov, memChprov, correctionHue, zone);
#ifdef _DEBUG
							if(correctionHue !=0.0) {
								if(fabs(correctionHue) > maxdeltaHueRY)
									maxdeltaHueRY=fabs(correctionHue);
								Munsry++;
							}
							if(fabs(correctionHue) > 0.45) depass++;  //verify if no bug in calculation
#endif
						}
						if(HH > 1.87 && HH < 3.09) {        //limits of green and green yellow
							zone=3; // *** green  yellow correction
							MunsellLch (Lprov, HH,Chprov, memChprov, correctionHue, zone);
#ifdef _DEBUG
							if(correctionHue !=0.0) {
								if(fabs(correctionHue) > maxdeltaHueGY)
									maxdeltaHueGY=fabs(correctionHue);
								Munsgy++;
							}
							if(fabs(correctionHue) > 0.45) depass++;   //verify if no bug in calculation
#endif
						}
						if(HH > -0.27 && HH <= 0.44) {        //limits of Red purple
							zone=4; // *** red purple correction
							MunsellLch (Lprov, HH,Chprov, memChprov, correctionHue, zone);
#ifdef _DEBUG
							if(correctionHue !=0.0) {
								if(fabs(correctionHue) > maxdeltaHueRP)
									maxdeltaHueRP=fabs(correctionHue);
								Munsrp++;
							}
							if(fabs(correctionHue) > 0.45) depass++;   //verify if no bug in calculation
#endif
						}
					}
				}
				//second gamut control take into account Munsell and saturation if R G B > 65535
				aprovn=Chprov*cos(HH+correctionHue);
				bprovn=Chprov*sin(HH+correctionHue);
				if (allwaysingamut) { // gamut control

					fyy = (0.00862069 *Lprov )+ 0.137932;
					fxx = (0.002 * aprovn) + fyy;
					fzz = fyy - (0.005 * bprovn);

					xx_ = 65535.0 * f2xyz(fxx)*D50x;
					yy_ = 65535.0 * f2xyz(fyy);
					zz_ = 65535.0 * f2xyz(fzz)*D50z;
					xyz2rgb(xx_,yy_,zz_,RR,GG,BB,wip);

					if(RR<0.0 || GG < 0.0 || BB < 0.0) {
#ifdef _DEBUG
						negsat++;
#endif
						Chprov*=0.95;
						inGamut=false;
					}
					if((!highlight) && (RR>65535.0 || GG > 65535.0 || BB>65535.0)) {// if "highlight reconstruction" enabled don't control Gamut for highlights
						//  if(RR>65535.0 || GG > 65535.0 || BB>65535.0) {
#ifdef _DEBUG
						moresat++;
#endif
						Chprov*=0.95;
						inGamut=false;
					}
				}
			} while (!inGamut);
			//put new values in Lab
			lab->L[i][j]=Lprov*327.68;
			lab->a[i][j]=aprovn*327.68;
			lab->b[i][j]=bprovn*327.68;

	}

}
/*
delete [] LL;
delete [] CC;
delete [] HH;
*/

#ifdef _DEBUG
	t2e.set();
	if (settings->verbose) {
		printf("Gamut: G1negat=%iiter G165535=%iiter G2negsat=%iiter G265535=%iiter\n",negat,moreRGB,negsat,moresat);
		printf("Munsell: MunPB=%ipix MunRY=%ipix MunGY=%ipix MunRP=%ipix MaxBP=%1.2frad MaxRY=%1.2frad MaxGY=%1.2frad MaxRP=%1.2frad  dep=%i\n",Munspb, Munsry,Munsgy,Munsrp, maxdeltaHueBP,maxdeltaHueRY,maxdeltaHueGY,maxdeltaHueRP, depass);
		printf("Vibrance  %d usec\n", t2e.etime(t1e));
	}
#endif
}

}
