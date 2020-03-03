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

#include "color.h"

using namespace std;

namespace rtengine
{

/*
 * MunsellLch correction
 * Copyright (c) 2012  Jacques Desmis <jdesmis@gmail.com>
 * Copyright (c) 2020  Ingo Weyrich <heckflosse67@gmx.de>
 *
 * Find the right LUT and calculate the correction
 */
void Color::MunsellLch (float lum, float hue, float chrom, float memChprov, float &correction, int zone, float &lbe, bool &correctL)
{

    const int x = memChprov;
    int y = chrom;

    //begin PB correction + sky
    if (zone == 1) {
        if (lum > 5.f) {
            if (lum < 15.f) {
                if (x <= 45 && (hue >= (_15PB10[x] - 0.035f)) && (hue < (_15PB10[x] + 0.052f))) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _15PB10[y] - _15PB10[x] ;
                    lbe = _15PB10[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= ( _3PB10[x] - 0.052f))  && (hue < (_45PB10[x] + _3PB10[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _3PB10[y] - _3PB10[x];
                    lbe = _3PB10[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_45PB10[x] + _3PB10[x]) * 0.5f)  && (hue < (_45PB10[x] + 0.052f))) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _45PB10[y] - _45PB10[x] ;
                    lbe = _45PB10[y];
                    correctL = true;
                } else if ((hue >= (_6PB10[x] - 0.052f)  && (hue < (_6PB10[x] + _75PB10[x]) * 0.5f))) {
                    correction = _6PB10[y] - _6PB10[x] ;
                    lbe = _6PB10[y];
                    correctL = true;
                } else if ((hue >= (_6PB10[x] + _75PB10[x]) * 0.5f)  && (hue < (_9PB10[x] + _75PB10[x]) * 0.5f)) {
                    correction = _75PB10[y] - _75PB10[x] ;
                    lbe = _75PB10[y];
                    correctL = true;
                } else if ((hue >= (_9PB10[x] + _75PB10[x]) * 0.5f)  && (hue < (_9PB10[x] + _10PB10[x]) * 0.5f)) {
                    correction = _9PB10[y] - _9PB10[x] ;
                    lbe = _9PB10[y];
                    correctL = true;
                } else if ((hue >= (_10PB10[x] + _9PB10[x]) * 0.5f)  && (hue < (_1P10[x] + _10PB10[x]) * 0.5f)) {
                    correction = _10PB10[y] - _10PB10[x] ;
                    lbe = _10PB10[y];
                    correctL = true;
                } else if ((hue >= (_10PB10[x] + _1P10[x]) * 0.5f)  && (hue < (_1P10[x] + _4P10[x]) * 0.5f)) {
                    correction = _1P10[y] - _1P10[x];
                    lbe = _1P10[y];
                    correctL = true;
                } else if ((hue >= (_1P10[x] + _4P10[x]) * 0.5f)  && (hue < (0.035f + _4P10[x]) * 0.5f)) {
                    correction = _4P10[y] - _4P10[x] ;
                    lbe = _4P10[y];
                    correctL = true;
                }
            } else if (lum < 25.f) {
                if (x <= 85 && (hue >= (_15PB20[x] - 0.035f)) && (hue < (_15PB20[x] + _3PB20[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _15PB20[y] - _15PB20[x] ;
                    lbe = _15PB20[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_15PB20[x] + _3PB20[x]) * 0.5f)  && (hue < (_45PB20[x] + _3PB20[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _3PB20[y] - _3PB20[x] ;
                    lbe = _3PB20[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_45PB20[x] + _3PB20[x]) * 0.5f)  && (hue < ( _45PB20[x] + 0.052f))) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _45PB20[y] - _45PB20[x] ;
                    lbe = _45PB20[y];
                    correctL = true;
                } else if ((hue >= (_45PB20[x] + 0.052f))  && (hue < (_6PB20[x] + _75PB20[x]) * 0.5f)) {
                    correction = _6PB20[y] - _6PB20[x];
                    lbe = _6PB20[y];
                    correctL = true;
                } else if ((hue >= (_6PB20[x] + _75PB20[x]) * 0.5f)  && (hue < (_9PB20[x] + _75PB20[x]) * 0.5f)) {
                    correction = _75PB20[y] - _75PB20[x] ;
                    lbe = _75PB20[y];
                    correctL = true;
                } else if ((hue >= (_9PB20[x] + _75PB20[x]) * 0.5f)  && (hue < (_9PB20[x] + _10PB20[x]) * 0.5f)) {
                    correction = _9PB20[y] - _9PB20[x] ;
                    lbe = _9PB20[y];
                    correctL = true;
                } else if ((hue >= (_10PB20[x] + _9PB20[x]) * 0.5f)  && (hue < (_1P20[x] + _10PB20[x]) * 0.5f)) {
                    correction = _10PB20[y] - _10PB20[x] ;
                    lbe = _10PB20[y];
                    correctL = true;
                } else if ((hue >= (_10PB20[x] + _1P20[x]) * 0.5f)  && (hue < (_1P20[x] + _4P20[x]) * 0.5f)) {
                    correction = _1P20[y] - _1P20[x] ;
                    lbe = _1P20[y];
                    correctL = true;
                } else if ((hue >= (_1P20[x] + _4P20[x]) * 0.5f)  && (hue < (0.035f + _4P20[x]) * 0.5f)) {
                    correction = _4P20[y] - _4P20[x] ;
                    lbe = _4P20[y];
                    correctL = true;
                }
            } else if (lum < 35.f) {
                if (x <= 85 && (hue >= (_15PB30[x] - 0.035f)) && (hue < (_15PB30[x] + _3PB30[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _15PB30[y] - _15PB30[x] ;
                    lbe = _15PB30[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_15PB30[x] + _3PB30[x]) * 0.5f)  && (hue < (_45PB30[x] + _3PB30[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _3PB30[y] - _3PB30[x] ;
                    lbe = _3PB30[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_45PB30[x] + _3PB30[x]) * 0.5f)  && (hue < (_45PB30[x] + 0.052f))) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _45PB30[y] - _45PB30[x] ;
                    lbe = _45PB30[y];
                    correctL = true;
                } else if ((hue >= ( _45PB30[x] + 0.052f))  && (hue < (_6PB30[x] + _75PB30[x]) * 0.5f)) {
                    correction = _6PB30[y] - _6PB30[x] ;
                    lbe = _6PB30[y];
                    correctL = true;
                } else if ((hue >= (_6PB30[x] + _75PB30[x]) * 0.5f)  && (hue < (_9PB30[x] + _75PB30[x]) * 0.5f)) {
                    correction = _75PB30[y] - _75PB30[x] ;
                    lbe = _75PB30[y] ;
                    correctL = true;
                } else if ((hue >= (_9PB30[x] + _75PB30[x]) * 0.5f)  && (hue < (_9PB30[x] + _10PB30[x]) * 0.5f)) {
                    correction = _9PB30[y] - _9PB30[x] ;
                    lbe = _9PB30[y];
                    correctL = true;
                } else if ((hue >= (_10PB30[x] + _9PB30[x]) * 0.5f)  && (hue < (_1P30[x] + _10PB30[x]) * 0.5f)) {
                    correction = _10PB30[y] - _10PB30[x] ;
                    lbe = _10PB30[y];
                    correctL = true;
                } else if ((hue >= (_10PB30[x] + _1P30[x]) * 0.5f)  && (hue < (_1P30[x] + _4P30[x]) * 0.5f)) {
                    correction = _1P30[y] - _1P30[x] ;
                    lbe = _1P30[y];
                    correctL = true;
                } else if ((hue >= (_1P30[x] + _4P30[x]) * 0.5f)  && (hue < (0.035f + _4P30[x]) * 0.5f)) {
                    correction = _4P30[y] - _4P30[x] ;
                    lbe = _4P30[y];
                    correctL = true;
                }
            } else if (lum < 45.f) {
                if (x < 75 && (hue <= (_05PB40[x] + _15PB40[x]) * 0.5f) && (hue > (_05PB40[x] + _10B40[x]) * 0.5f)) {
                    if (y > 75) {
                        y = 75;
                    }
                    correction = _05PB40[y] - _05PB40[x] ;
                    lbe = _05PB40[y];
                    correctL = true;
                } else if (x < 70 && (hue <= (_05PB40[x] + _10B40[x]) * 0.5f) && (hue > (_10B40[x] + _9B40[x]) * 0.5f)) {
                    if (y > 70) {
                        y = 70;
                    }
                    correction = _10B40[y] - _10B40[x] ;
                    lbe = _10B40[y];
                    correctL = true;
                } else if (x < 70 && (hue <= (_10B40[x] + _9B40[x]) * 0.5f) && (hue > (_9B40[x] + _7B40[x]) * 0.5f)) {
                    if (y > 70) {
                        y = 70;
                    }
                    correction = _9B40[y] - _9B40[x] ;
                    lbe = _9B40[y];
                    correctL = true;
                } else if (x < 70 && (hue <= (_9B40[x] + _7B40[x]) * 0.5f) && (hue > (_5B40[x] + _7B40[x]) * 0.5f)) {
                    if (y > 70) {
                        y = 70;
                    }
                    correction = _7B40[y] - _7B40[x] ;
                    lbe = _7B40[y];
                    correctL = true;
                } else if (x < 70 && (hue <= (_5B40[x] + _7B40[x]) * 0.5f)  && (hue > (_5B40[x] - 0.035f))) {
                    if (y > 70) {
                        y = 70;    //
                    }
                    correction = _5B40[y] - _5B40[x] ;
                    lbe =  _5B40[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_15PB40[x] - 0.035f)) && (hue < (_15PB40[x] + _3PB40[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _15PB40[y] - _15PB40[x] ;
                    lbe = _15PB40[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_15PB40[x] + _3PB40[x]) * 0.5f)  && (hue < (_45PB40[x] + _3PB40[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _3PB40[y] - _3PB40[x] ;
                    lbe = _3PB40[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_45PB40[x] + _3PB40[x]) * 0.5f)  && (hue < (_45PB40[x] + 0.052f))) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _45PB40[y] - _45PB40[x] ;
                    lbe = _45PB40[y] ;
                    correctL = true;
                } else if ((hue >= (_45PB40[x] + 0.052f))  && (hue < (_6PB40[x] + _75PB40[x]) * 0.5f)) {
                    correction = _6PB40[y] - _6PB40[x] ;
                    lbe = _6PB40[y];
                    correctL = true;
                } else if ((hue >= (_6PB40[x] + _75PB40[x]) * 0.5f)  && (hue < (_9PB40[x] + _75PB40[x]) * 0.5f)) {
                    correction = _75PB40[y] - _75PB40[x] ;
                    lbe = _75PB40[y];
                    correctL = true;
                } else if ((hue >= (_9PB40[x] + _75PB40[x]) * 0.5f)  && (hue < (_9PB40[x] + _10PB40[x]) * 0.5f)) {
                    correction = _9PB40[y] - _9PB40[x] ;
                    lbe = _9PB40[y];
                    correctL = true;
                } else if ((hue >= (_10PB40[x] + _9PB40[x]) * 0.5f)  && (hue < (_1P40[x] + _10PB40[x]) * 0.5f)) {
                    correction = _10PB40[y] - _10PB40[x] ;
                    lbe = _10PB40[y];
                    correctL = true;
                } else if ((hue >= (_10PB40[x] + _1P40[x]) * 0.5f)  && (hue < (_1P40[x] + _4P40[x]) * 0.5f)) {
                    correction = _1P40[y] - _1P40[x] ;
                    lbe = _1P40[y];
                    correctL = true;
                } else if ((hue >= (_1P40[x] + _4P40[x]) * 0.5f)  && (hue < (0.035f + _4P40[x]) * 0.5f)) {
                    correction = _4P40[y] - _4P40[x] ;
                    lbe = _4P40[y];
                    correctL = true;
                }
            } else if (lum < 55.f) {
                if (x < 79 && (hue <= (_05PB50[x] + _15PB50[x]) * 0.5f) && (hue > (_05PB50[x] + _10B50[x]) * 0.5f)) {
                    if (y > 79) {
                        y = 79;
                    }
                    correction = _05PB50[y] - _05PB50[x] ;
                    lbe = _05PB50[y];
                    correctL = true;
                } else if (x < 79 && (hue <= (_05PB50[x] + _10B50[x]) * 0.5f) && (hue > (_10B50[x] + _9B50[x]) * 0.5f)) {
                    if (y > 79) {
                        y = 79;
                    }
                    correction = _10B50[y] - _10B50[x] ;
                    lbe = _10B50[y];
                    correctL = true;
                } else if (x < 79 && (hue <= (_10B50[x] + _9B50[x]) * 0.5f) && (hue > (_9B50[x] + _7B50[x]) * 0.5f)) {
                    if (y > 79) {
                        y = 79;
                    }
                    correction = _9B50[y] - _9B50[x] ;
                    lbe = _9B50[y];
                    correctL = true;
                } else if (x < 79 && (hue <= (_9B50[x] + _7B50[x]) * 0.5f) && (hue > (_5B50[x] + _7B50[x]) * 0.5f)) {
                    if (y > 79) {
                        y = 79;
                    }
                    correction = _7B50[y] - _7B50[x] ;
                    lbe = _7B50[y];
                    correctL = true;
                } else if (x < 79 && (hue <= (_5B50[x] + _7B50[x]) * 0.5f)  && (hue > (_5B50[x] - 0.035f))) {
                    if (y > 79) {
                        y = 79;    //
                    }
                    correction = _5B50[y] - _5B50[x] ;
                    lbe = _5B50[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_15PB50[x] - 0.035f)) && (hue < (_15PB50[x] + _3PB50[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _15PB50[y] - _15PB50[x] ;
                    lbe = _15PB50[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_15PB50[x] + _3PB50[x]) * 0.5f)  && (hue < (_45PB50[x] + _3PB50[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _3PB50[y] - _3PB50[x] ;
                    lbe = _3PB50[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_45PB50[x] + _3PB50[x]) * 0.5f)  && (hue < (_6PB50[x] + _45PB50[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _45PB50[y] - _45PB50[x] ;
                    lbe = _45PB50[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_6PB50[x] + _45PB50[x]) * 0.5f)  && (hue < (_6PB50[x] + _75PB50[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _6PB50[y] - _6PB50[x] ;
                    lbe = _6PB50[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_6PB50[x] + _75PB50[x]) * 0.5f)  && (hue < (_9PB50[x] + _75PB50[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _75PB50[y] - _75PB50[x] ;
                    lbe = _75PB50[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_9PB50[x] + _75PB50[x]) * 0.5f)  && (hue < (_9PB50[x] + _10PB50[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _9PB50[y] - _9PB50[x] ;
                    lbe = _9PB50[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_10PB50[x] + _9PB50[x]) * 0.5f)  && (hue < (_1P50[x] + _10PB50[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _10PB50[y] - _10PB50[x] ;
                    lbe = _10PB50[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_10PB50[x] + _1P50[x]) * 0.5f)  && (hue < (_1P50[x] + _4P50[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _1P50[y] - _1P50[x] ;
                    lbe = _1P50[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_1P50[x] + _4P50[x]) * 0.5f)  && (hue < (0.035f + _4P50[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _4P50[y] - _4P50[x] ;
                    lbe = _4P50[y];
                    correctL = true;
                }
            } else if (lum < 65.f) {
                if (x < 79 && (hue <= (_05PB60[x] + _15PB60[x]) * 0.5f) && (hue > (_05PB60[x] + _10B60[x]) * 0.5f)) {
                    if (y > 79) {
                        y = 79;
                    }
                    correction = _05PB60[y] - _05PB60[x] ;
                    lbe = _05PB60[y];
                    correctL = true;
                } else if (x < 79 && (hue <= (_05PB60[x] + _10B60[x]) * 0.5f) && (hue > (_10B60[x] + _9B60[x]) * 0.5f)) {
                    if (y > 79) {
                        y = 79;
                    }
                    correction = _10B60[y] - _10B60[x] ;
                    lbe = _10B60[y];
                    correctL = true;
                } else if (x < 79 && (hue <= (_10B60[x] + _9B60[x]) * 0.5f) && (hue > (_9B60[x] + _7B60[x]) * 0.5f)) {
                    if (y > 79) {
                        y = 79;
                    }
                    correction = _9B60[y] - _9B60[x] ;
                    lbe = _9B60[y];
                    correctL = true;
                } else if (x < 79 && (hue <= (_9B60[x] + _7B60[x]) * 0.5f) && (hue > (_5B60[x] + _7B60[x]) * 0.5f)) {
                    if (y > 79) {
                        y = 79;
                    }
                    correction = _7B60[y] - _7B60[x] ;
                    lbe = _7B60[y];
                    correctL = true;
                } else if (x < 79 && (hue <= (_5B60[x] + _7B60[x]) * 0.5f)  && (hue > (_5B60[x] - 0.035f))) {
                    if (y > 79) {
                        y = 79;    //
                    }
                    correction = _5B60[y] - _5B60[x] ;
                    lbe = _5B60[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_15PB60[x] - 0.035f)) && (hue < (_15PB60[x] + _3PB60[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _15PB60[y] - _15PB60[x] ;
                    lbe = _15PB60[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_15PB60[x] + _3PB60[x]) * 0.5f)  && (hue < (_45PB60[x] + _3PB60[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _3PB60[y] - _3PB60[x] ;
                    lbe = _3PB60[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_45PB60[x] + _3PB60[x]) * 0.5f)  && (hue < (_6PB60[x] + _45PB60[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _45PB60[y] - _45PB60[x] ;
                    lbe = _45PB60[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_6PB60[x] + _45PB60[x]) * 0.5f)  && (hue < (_6PB60[x] + _75PB60[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _6PB60[y] - _6PB60[x] ;
                    lbe = _6PB60[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_6PB60[x] + _75PB60[x]) * 0.5f)  && (hue < (_9PB60[x] + _75PB60[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _75PB60[y] - _75PB60[x] ;
                    lbe = _75PB60[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_9PB60[x] + _75PB60[x]) * 0.5f)  && (hue < (_9PB60[x] + _10PB60[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _9PB60[y] - _9PB60[x] ;
                    lbe = _9PB60[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_10PB60[x] + _9PB60[x]) * 0.5f)  && (hue < (_1P60[x] + _10PB60[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _10PB60[y] - _10PB60[x] ;
                    lbe = _10PB60[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_10PB60[x] + _1P60[x]) * 0.5f)  && (hue < (_1P60[x] + _4P60[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _1P60[y] - _1P60[x] ;
                    lbe = _1P60[y];
                    correctL = true;
                } else if (x <= 85 && (hue >= (_1P60[x] + _4P60[x]) * 0.5f)  && (hue < (0.035f + _4P60[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _4P60[y] - _4P60[x] ;
                    lbe = _4P60[y];
                    correctL = true;
                }
            } else if (lum < 75.f) {
                if (x < 50 && (hue <= (_05PB70[x] + _15PB70[x]) * 0.5f) && (hue > (_05PB70[x] + _10B70[x]) * 0.5f)) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _05PB70[y] - _05PB70[x] ;
                    lbe = _05PB70[y];
                    correctL = true;
                } else if (x < 50 && (hue <= (_05PB70[x] + _10B70[x]) * 0.5f) && (hue > (_10B70[x] + _9B70[x]) * 0.5f)) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _10B70[y] - _10B70[x] ;
                    lbe = _10B70[y];
                    correctL = true;
                } else if (x < 50 && (hue <= (_10B70[x] + _9B70[x]) * 0.5f) && (hue > (_9B70[x] + _7B70[x]) * 0.5f)) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _9B70[y] - _9B70[x] ;
                    lbe = _9B70[y];
                    correctL = true;
                } else if (x < 50 && (hue <= (_9B70[x] + _7B70[x]) * 0.5f) && (hue > (_5B70[x] + _7B70[x]) * 0.5f)) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _7B70[y] - _7B70[x] ;
                    lbe = _7B70[y];
                    correctL = true;
                } else if (x < 50 && (hue <= (_5B70[x] + _7B70[x]) * 0.5f)  && (hue > (_5B70[x] - 0.035f))) {
                    if (y > 49) {
                        y = 49;    //
                    }
                    correction = _5B70[y] - _5B70[x] ;
                    lbe =  _5B70[y];
                    correctL = true;
                } else if (x < 50 && (hue >= (_15PB70[x] - 0.035f)) && (hue < (_15PB70[x] + _3PB70[x]) * 0.5f)) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _15PB70[y] - _15PB70[x] ;
                    lbe = _15PB70[y];
                    correctL = true;
                } else if (x < 50 && (hue >= (_45PB70[x] + _3PB70[x]) * 0.5f)  && (hue < (_6PB70[x] + _45PB70[x]) * 0.5f)) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _45PB70[y] - _45PB70[x] ;
                    lbe = _45PB70[y];
                    correctL = true;
                } else if (x < 50 && (hue >= (_6PB70[x] + _45PB70[x]) * 0.5f)  && (hue < (_6PB70[x] + _75PB70[x]) * 0.5f)) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _6PB70[y] - _6PB70[x] ;
                    lbe = _6PB70[y];
                    correctL = true;
                } else if (x < 50 && (hue >= (_6PB70[x] + _75PB70[x]) * 0.5f)  && (hue < (_9PB70[x] + _75PB70[x]) * 0.5f)) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _75PB70[y] - _75PB70[x] ;
                    lbe = _75PB70[y];
                    correctL = true;
                } else if (x < 50 && (hue >= (_9PB70[x] + _75PB70[x]) * 0.5f)  && (hue < (_9PB70[x] + 0.035f))) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _9PB70[y] - _9PB70[x] ;
                    lbe = _9PB70[y];
                    correctL = true;
                }
            } else if (lum < 85.f) {
                if (x < 40 && (hue <= (_05PB80[x] + _15PB80[x]) * 0.5f) && (hue > (_05PB80[x] + _10B80[x]) * 0.5f)) {
                    if (y > 39) {
                        y = 39;
                    }
                    correction = _05PB80[y] - _05PB80[x] ;
                    lbe = _05PB80[y] ;
                    correctL = true;
                } else if (x < 40 && (hue <= (_05PB80[x] + _10B80[x]) * 0.5f) && (hue > (_10B80[x] + _9B80[x]) * 0.5f)) {
                    if (y > 39) {
                        y = 39;
                    }
                    correction = _10B80[y] - _10B80[x] ;
                    lbe = _10B80[y];
                    correctL = true;
                } else if (x < 40 && (hue <= (_10B80[x] + _9B80[x]) * 0.5f) && (hue > (_9B80[x] + _7B80[x]) * 0.5f)) {
                    if (y > 39) {
                        y = 39;
                    }
                    correction = _9B80[y] - _9B80[x] ;
                    lbe = _9B80[y];
                    correctL = true;
                } else if (x < 50 &&(hue <= (_9B80[x] + _7B80[x]) * 0.5f) && (hue > (_5B80[x] + _7B80[x]) * 0.5f)) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _7B80[y] - _7B80[x] ;
                    lbe = _7B80[y];
                    correctL = true;
                } else if (x < 50 && (hue <= (_5B80[x] + _7B80[x]) * 0.5f)  && (hue > (_5B80[x] - 0.035f))) {
                    if (y > 49) {
                        y = 49;    //
                    }
                    correction = _5B80[y] - _5B80[x] ;
                    lbe = _5B80[y];
                    correctL = true;
                } else if (x < 50 && (hue >= (_15PB80[x] - 0.035f)) && (hue < (_15PB80[x] + _3PB80[x]) * 0.5f)) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _15PB80[y] - _15PB80[x] ;
                    lbe = _15PB80[y];
                    correctL = true;
                } else if (x < 50 && (hue >= (_45PB80[x] + _3PB80[x]) * 0.5f)  && (hue < (_6PB80[x] + _45PB80[x]) * 0.5f)) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _45PB80[y] - _45PB80[x] ;
                    lbe = _45PB80[y];
                    correctL = true;
                } else if (x < 50 && (hue >= (_6PB80[x] + _45PB80[x]) * 0.5f)  && (hue < (_6PB80[x] + _75PB80[x]) * 0.5f)) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _6PB80[y] - _6PB80[x] ;
                    lbe = _6PB80[y];
                    correctL = true;
                } else if (x < 50 && (hue >= (_6PB80[x] + _75PB80[x]) * 0.5f)  && (hue < (_9PB80[x] + _75PB80[x]) * 0.5f)) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _75PB80[y] - _75PB80[x] ;
                    lbe = _75PB80[y];
                    correctL = true;
                } else if (x < 50 && (hue >= (_9PB80[x] + _75PB80[x]) * 0.5f)  && (hue < (_9PB80[x] + 0.035f))) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _9PB80[y] - _9PB80[x] ;
                    lbe = _9PB80[y];
                    correctL = true;
                }
            }
        }
    } else if (zone == 2) { //red yellow correction
        if (lum > 15.f) {
            if (lum < 25.f) {
                if (x <= 45 && (hue <= (_10YR20[x] + 0.035f)) && (hue > (_10YR20[x] + _85YR20[x]) * 0.5f)) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _10YR20[y] - _10YR20[x] ;
                    lbe = _10YR20[y];
                    correctL = true;
                } else if (x <= 45 && (hue <= (_85YR20[x] + _10YR20[x]) * 0.5f)  && (hue > (_85YR20[x] + 0.035f))) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _85YR20[y] - _85YR20[x] ;
                    lbe = _85YR20[y];
                    correctL = true;
                }
            } else if (lum < 35.f) {
                if (x < 85 && (hue <= (_10YR30[x] + 0.035f)) && (hue > (_10YR30[x] + _85YR30[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _10YR30[y] - _10YR30[x] ;
                    lbe = _10YR30[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_10YR30[x] + _85YR30[x]) * 0.5f) && (hue > (_85YR30[x] + _7YR30[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _85YR30[y] - _85YR30[x] ;
                    lbe = _85YR30[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_85YR30[x] + _7YR30[x]) * 0.5f)  && (hue > (_7YR30[x] + _55YR30[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _7YR30[y] - _7YR30[x] ;
                    lbe = _7YR30[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_7YR30[x] + _55YR30[x]) * 0.5f)  && (hue > (_55YR30[x] + _4YR30[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _55YR30[y] - _55YR30[x] ;
                    lbe = _55YR30[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_55YR30[x] + _4YR30[x]) * 0.5f)  && (hue > (_4YR30[x] + _25YR30[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _4YR30[y] - _4YR30[x] ;
                    lbe = _4YR30[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_4YR30[x] + _25YR30[x]) * 0.5f)  && (hue > (_25YR30[x] + _10R30[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _25YR30[y] - _25YR30[x] ;
                    lbe = _25YR30[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_25YR30[x] + _10R30[x]) * 0.5f)  && (hue > (_10R30[x] + _9R30[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _10R30[y] - _10R30[x] ;
                    lbe = _10R30[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_10R30[x] + _9R30[x]) * 0.5f)  && (hue > (_9R30[x] + _7R30[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _9R30[y] - _9R30[x] ;
                    lbe = _9R30[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_9R30[x] + _7R30[x]) * 0.5f)  && (hue > (_7R30[x] - 0.035f))) {
                    if (y > 89) {
                        y = 89;
                    }

                    correction = _7R30[y] - _7R30[x] ;
                    lbe = _7R30[y] ;
                    correctL = true;
                }
            } else if (lum < 45.f) {
                if (x < 85 && (hue <= (_10YR40[x] + 0.035f)) && (hue > (_10YR40[x] + _85YR40[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _10YR40[y] - _10YR40[x] ;
                    lbe = _10YR40[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_10YR40[x] + _85YR40[x]) * 0.5f) && (hue > (_85YR40[x] + _7YR40[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _85YR40[y] - _85YR40[x] ;
                    lbe = _85YR40[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_85YR40[x] + _7YR40[x]) * 0.5f)  && (hue > (_7YR40[x] + _55YR40[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _7YR40[y] - _7YR40[x] ;
                    lbe = _7YR40[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_7YR40[x] + _55YR40[x]) * 0.5f)  && (hue > (_55YR40[x] + _4YR40[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _55YR40[y] - _55YR40[x] ;
                    lbe = _55YR40[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_55YR40[x] + _4YR40[x]) * 0.5f)  && (hue > (_4YR40[x] + _25YR40[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _4YR40[y] - _4YR40[x] ;
                    lbe = _4YR40[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_4YR40[x] + _25YR40[x]) * 0.5f)  && (hue > (_25YR40[x] + _10R40[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _25YR40[y] - _25YR40[x] ;
                    lbe = _25YR40[y] ;
                    correctL = true;
                } else if (x < 85 && (hue <= (_25YR40[x] + _10R40[x]) * 0.5f)  && (hue > (_10R40[x] + _9R40[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _10R40[y] - _10R40[x] ;
                    lbe = _10R40[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_10R40[x] + _9R40[x]) * 0.5f)  && (hue > (_9R40[x] + _7R40[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _9R40[y] - _9R40[x] ;
                    lbe = _9R40[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_9R40[x] + _7R40[x]) * 0.5f)  && (hue > (_7R40[x] - 0.035f))) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _7R40[y] - _7R40[x] ;
                    lbe = _7R40[y];
                    correctL = true;
                }
            } else if (lum < 55.f) {
                if (x < 85 && (hue <= (_10YR50[x] + 0.035f)) && (hue > (_10YR50[x] + _85YR50[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _10YR50[y] - _10YR50[x] ;
                    lbe = _10YR50[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_10YR50[x] + _85YR50[x]) * 0.5f) && (hue > (_85YR50[x] + _7YR50[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _85YR50[y] - _85YR50[x] ;
                    lbe = _85YR50[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_85YR50[x] + _7YR50[x]) * 0.5f)  && (hue > (_7YR50[x] + _55YR50[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _7YR50[y] - _7YR50[x] ;
                    lbe = _7YR50[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_7YR50[x] + _55YR50[x]) * 0.5f)  && (hue > (_55YR50[x] + _4YR50[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _55YR50[y] - _55YR50[x] ;
                    lbe = _55YR50[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_55YR50[x] + _4YR50[x]) * 0.5f)  && (hue > (_4YR50[x] + _25YR50[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _4YR50[y] - _4YR50[x] ;
                    lbe = _4YR50[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_4YR50[x] + _25YR50[x]) * 0.5f)  && (hue > (_25YR50[x] + _10R50[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _25YR50[y] - _25YR50[x] ;
                    lbe = _25YR50[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_25YR50[x] + _10R50[x]) * 0.5f)  && (hue > (_10R50[x] + _9R50[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _10R50[y] - _10R50[x] ;
                    lbe = _10R50[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_10R50[x] + _9R50[x]) * 0.5f)  && (hue > (_9R50[x] + _7R50[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _9R50[y] - _9R50[x] ;
                    lbe = _9R50[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_9R50[x] + _7R50[x]) * 0.5f)  && (hue > (_7R50[x] - 0.035f))) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _7R50[y] - _7R50[x] ;
                    lbe = _7R50[y];
                    correctL = true;
                }
            } else if (lum < 65.f) {
                if ((hue <= (_10YR60[x] + 0.035f)) && (hue > (_10YR60[x] + _85YR60[x]) * 0.5f)) {
                    correction = _10YR60[y] - _10YR60[x] ;
                    lbe = _10YR60[y];
                    correctL = true;
                } else if ((hue <= (_10YR60[x] + _85YR60[x]) * 0.5f) && (hue > (_85YR60[x] + _7YR60[x]) * 0.5f)) {
                    correction = _85YR60[y] - _85YR60[x] ;
                    lbe = _85YR60[y];
                    correctL = true;
                } else if ((hue <= (_85YR60[x] + _7YR60[x]) * 0.5f)  && (hue > (_7YR60[x] + _55YR60[x]) * 0.5f)) {
                    correction = _7YR60[y] - _7YR60[x] ;
                    lbe = _7YR60[y];
                    correctL = true;
                } else if ((hue <= (_7YR60[x] + _55YR60[x]) * 0.5f)  && (hue > (_55YR60[x] + _4YR60[x]) * 0.5f)) {
                    correction = _55YR60[y] - _55YR60[x] ;
                    lbe = _55YR60[y];
                    correctL = true;
                } else if ((hue <= (_55YR60[x] + _4YR60[x]) * 0.5f)  && (hue > (_4YR60[x] + _25YR60[x]) * 0.5f)) {
                    correction = _4YR60[y] - _4YR60[x] ;
                    lbe = _4YR60[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_4YR60[x] + _25YR60[x]) * 0.5f)  && (hue > (_25YR60[x] + _10R60[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _25YR60[y] - _25YR60[x] ;
                    lbe = _25YR60[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_25YR60[x] + _10R60[x]) * 0.5f)  && (hue > (_10R60[x] + _9R60[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _10R60[y] - _10R60[x] ;
                    lbe = _10R60[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_10R60[x] + _9R60[x]) * 0.5f)  && (hue > (_9R60[x] + _7R60[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _9R60[y] - _9R60[x] ;
                    lbe = _9R60[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_9R60[x] + _7R60[x]) * 0.5f)  && (hue > (_7R60[x] - 0.035f))) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _7R60[y] - _7R60[x] ;
                    lbe = _7R60[y];
                    correctL = true;
                }
            } else if (lum < 75.f) {
                if ((hue <= (_10YR70[x] + 0.035f)) && (hue > (_10YR70[x] + _85YR70[x]) * 0.5f)) {
                    correction = _10YR70[y] - _10YR70[x] ;
                    lbe = _10YR70[y];
                    correctL = true;
                } else if ((hue <= (_10YR70[x] + _85YR70[x]) * 0.5f) && (hue > (_85YR70[x] + _7YR70[x]) * 0.5f)) {
                    correction = _85YR70[y] - _85YR70[x] ;
                    lbe = _85YR70[y];
                    correctL = true;
                } else if ((hue <= (_85YR70[x] + _7YR70[x]) * 0.5f)  && (hue > (_7YR70[x] + _55YR70[x]) * 0.5f)) {
                    correction = _7YR70[y] - _7YR70[x] ;
                    lbe = _7YR70[y];
                    correctL = true;
                } else if ((hue <= (_7YR70[x] + _55YR70[x]) * 0.5f)  && (hue > (_55YR70[x] + _4YR70[x]) * 0.5f)) {
                    correction = _55YR70[y] - _55YR70[x] ;
                    lbe = _55YR70[y];
                    correctL = true;
                } else if ((hue <= (_55YR70[x] + _4YR70[x]) * 0.5f)  && (hue > (_4YR70[x] + _25YR70[x]) * 0.5f)) {
                    correction = _4YR70[y] - _4YR70[x] ;
                    lbe = _4YR70[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_4YR70[x] + _25YR70[x]) * 0.5f)  && (hue > (_25YR70[x] + _10R70[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _25YR70[y] - _25YR70[x] ;
                    lbe = _25YR70[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_25YR70[x] + _10R70[x]) * 0.5f)  && (hue > (_10R70[x] + _9R70[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _10R70[y] - _10R70[x] ;
                    lbe = _10R70[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_10R70[x] + _9R70[x]) * 0.5f)  && (hue > (_9R70[x] + _7R70[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _9R70[y] - _9R70[x] ;
                    lbe = _9R70[y] ;
                    correctL = true;
                } else if (x < 85 && (hue <= (_9R70[x] + _7R70[x]) * 0.5f)  && (hue > (_7R70[x] - 0.035f))) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _7R70[y] - _7R70[x] ;
                    lbe = _7R70[y];
                    correctL = true;
                }
            } else if (lum < 85.f) {
                if ((hue <= (_10YR80[x] + 0.035f)) && (hue > (_10YR80[x] + _85YR80[x]) * 0.5f)) {
                    correction = _10YR80[y] - _10YR80[x] ;
                    lbe = _10YR80[y];
                    correctL = true;
                } else if ((hue <= (_10YR80[x] + _85YR80[x]) * 0.5f) && (hue > (_85YR80[x] + _7YR80[x]) * 0.5f)) {
                    correction = _85YR80[y] - _85YR80[x] ;
                    lbe = _85YR80[y];
                } else if (x < 85 && (hue <= (_85YR80[x] + _7YR80[x]) * 0.5f)  && (hue > (_7YR80[x] + _55YR80[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _7YR80[y] - _7YR80[x] ;
                    lbe = _7YR80[y];
                    correctL = true;
                } else if (x < 45 && (hue <= (_7YR80[x] + _55YR80[x]) * 0.5f)  && (hue > (_55YR80[x] + _4YR80[x]) * 0.5f)) {
                    correction = _55YR80[y] - _55YR80[x] ;
                    lbe = _55YR80[y];
                    correctL = true;
                } else if (x < 45 && (hue <= (_55YR80[x] + _4YR80[x]) * 0.5f)  && (hue > (_4YR80[x] - 0.035f))) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _4YR80[y] - _4YR80[x] ;
                    lbe = _4YR80[y] ;
                    correctL = true;
                }
            } else if (lum < 95.f) {
                if (x < 85 && (hue <= (_10YR90[x] + 0.035f)) && (hue > (_10YR90[x] - 0.035f))) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _10YR90[y] - _10YR90[x] ;
                    lbe = _10YR90[y];
                    correctL = true;
                } else if (x < 85 && hue <= (_85YR90[x] + 0.035f)  && hue > (_85YR90[x] - 0.035f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _85YR90[y] - _85YR90[x] ;
                    lbe = _85YR90[y];
                    correctL = true;
                } else if (x < 45 && (hue <= (_55YR90[x] + 0.035f)  && (hue > (_55YR90[x] - 0.035f)))) {
                    if (y > 49) {
                        y = 49;
                    }
                    correction = _55YR90[y] - _55YR90[x] ;
                    lbe = _55YR90[y];
                    correctL = true;
                }
            }
        }
    } else if (zone == 3) { // //Green yellow correction
        if (lum >= 25.f) {
            if (lum < 35.f) {
                if ((hue <= (_7G30[x] + 0.035f)) && (hue > (_7G30[x] + _5G30[x]) * 0.5f)) {
                    correction = _7G30[y] - _7G30[x] ;
                    lbe = _7G30[y];
                    correctL = true;
                } else if ((hue <= (_7G30[x] + _5G30[x]) * 0.5f) && (hue > (_5G30[x] + _25G30[x]) * 0.5f)) {
                    correction = _5G30[y] - _5G30[x] ;
                    lbe = _5G30[y];
                    correctL = true;
                } else if ((hue <= (_25G30[x] + _5G30[x]) * 0.5f)  && (hue > (_25G30[x] + _1G30[x]) * 0.5f)) {
                    correction = _25G30[y] - _25G30[x] ;
                    lbe = _25G30[y];
                    correctL = true;
                } else if ((hue <= (_1G30[x] + _25G30[x]) * 0.5f)  && (hue > (_1G30[x] + _10GY30[x]) * 0.5f)) {
                    correction = _1G30[y] - _1G30[x] ;
                    lbe = _1G30[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_1G30[x] + _10GY30[x]) * 0.5f)  && (hue > (_10GY30[x] + _75GY30[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _10GY30[y] - _10GY30[x] ;
                    lbe =  _10GY30[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_10GY30[x] + _75GY30[x]) * 0.5f)  && (hue > (_75GY30[x] + _5GY30[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _75GY30[y] - _75GY30[x] ;
                    lbe = _75GY30[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_5GY30[x] + _75GY30[x]) * 0.5f)  && (hue > (_5GY30[x] - 0.035f))) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _5GY30[y] - _5GY30[x] ;
                    lbe = _5GY30[y] ;
                    correctL = true;
                }
            } else if (lum < 45.f) {
                if ((hue <= (_7G40[x] + 0.035f)) && (hue > (_7G40[x] + _5G40[x]) * 0.5f)) {
                    correction = _7G40[y] - _7G40[x] ;
                    lbe = _7G40[y];
                    correctL = true;
                } else if ((hue <= (_7G40[x] + _5G40[x]) * 0.5f) && (hue > (_5G40[x] + _25G40[x]) * 0.5f)) {
                    correction = _5G40[y] - _5G40[x] ;
                    lbe = _5G40[y];
                    correctL = true;
                } else if ((hue <= (_25G40[x] + _5G40[x]) * 0.5f)  && (hue > (_25G40[x] + _1G40[x]) * 0.5f)) {
                    correction = _25G40[y] - _25G40[x] ;
                    lbe = _25G40[y];
                    correctL = true;
                } else if ((hue <= (_1G40[x] + _25G40[x]) * 0.5f)  && (hue > (_1G40[x] + _10GY40[x]) * 0.5f)) {
                    correction = _1G40[y] - _1G40[x] ;
                    lbe = _1G40[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_1G40[x] + _10GY40[x]) * 0.5f)  && (hue > (_10GY40[x] + _75GY40[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _10GY40[y] - _10GY40[x] ;
                    lbe = _10GY40[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_10GY40[x] + _75GY40[x]) * 0.5f)  && (hue > (_75GY40[x] + _5GY40[x]) * 0.5f)) {
                    if (y > 89) {
                        y = 89;
                    }
                    correction = _75GY40[y] - _75GY40[x] ;
                    lbe = _75GY40[y];
                    correctL = true;
                } else if (x < 85 && (hue <= (_5GY40[x] + _75GY40[x]) * 0.5f)  && (hue > (_5GY40[x] - 0.035f))) {
                    if (y > 89) {
                        y = 89;    //
                    }
                    correction = _5GY40[y] - _5GY40[x] ;
                    lbe = _5GY40[y];
                    correctL = true;
                }
            } else if (lum < 55.f) {
                if ((hue <= (_7G50[x] + 0.035f)) && (hue > (_7G50[x] + _5G50[x]) * 0.5f)) {
                    correction = _7G50[y] - _7G50[x] ;
                    lbe = _7G50[y];
                    correctL = true;
                } else if ((hue <= (_7G50[x] + _5G50[x]) * 0.5f) && (hue > (_5G50[x] + _25G50[x]) * 0.5f)) {
                    correction = _5G50[y] - _5G50[x] ;
                    lbe = _5G50[y];
                    correctL = true;
                } else if ((hue <= (_25G50[x] + _5G50[x]) * 0.5f)  && (hue > (_25G50[x] + _1G50[x]) * 0.5f)) {
                    correction = _25G50[y] - _25G50[x] ;
                    lbe = _25G50[y];
                    correctL = true;
                } else if ((hue <= (_1G50[x] + _25G50[x]) * 0.5f)  && (hue > (_1G50[x] + _10GY50[x]) * 0.5f)) {
                    correction = _1G50[y] - _1G50[x] ;
                    lbe = _1G50[y];
                    correctL = true;
                } else if ((hue <= (_1G50[x] + _10GY50[x]) * 0.5f)  && (hue > (_10GY50[x] + _75GY50[x]) * 0.5f)) {
                    correction = _10GY50[y] - _10GY50[x] ;
                    lbe = _10GY50[y];
                    correctL = true;
                } else if ((hue <= (_10GY50[x] + _75GY50[x]) * 0.5f)  && (hue > (_75GY50[x] + _5GY50[x]) * 0.5f)) {
                    correction = _75GY50[y] - _75GY50[x] ;
                    lbe = _75GY50[y];
                    correctL = true;
                } else if ((hue <= (_5GY50[x] + _75GY50[x]) * 0.5f)  && (hue > (_5GY50[x] - 0.035f))) {
                    correction = _5GY50[y] - _5GY50[x] ;
                    lbe = _5GY50[y];
                    correctL = true;
                }
            } else if (lum < 65.f) {
                if ((hue <= (_7G60[x] + 0.035f)) && (hue > (_7G60[x] + _5G60[x]) * 0.5f)) {
                    correction = _7G60[y] - _7G60[x] ;
                    lbe = _7G60[y];
                    correctL = true;
                } else if ((hue <= (_7G60[x] + _5G60[x]) * 0.5f) && (hue > (_5G60[x] + _25G60[x]) * 0.5f)) {
                    correction = _5G60[y] - _5G60[x] ;
                    lbe = _5G60[y];
                    correctL = true;
                } else if ((hue <= (_25G60[x] + _5G60[x]) * 0.5f)  && (hue > (_25G60[x] + _1G60[x]) * 0.5f)) {
                    correction = _25G60[y] - _25G60[x] ;
                    lbe = _25G60[y];
                    correctL = true;
                } else if ((hue <= (_1G60[x] + _25G60[x]) * 0.5f)  && (hue > (_1G60[x] + _10GY60[x]) * 0.5f)) {
                    correction = _1G60[y] - _1G60[x] ;
                    lbe = _1G60[y];
                    correctL = true;
                } else if ((hue <= (_1G60[x] + _10GY60[x]) * 0.5f)  && (hue > (_10GY60[x] + _75GY60[x]) * 0.5f)) {
                    correction = _10GY60[y] - _10GY60[x] ;
                    lbe = _10GY60[y];
                    correctL = true;
                } else if ((hue <= (_10GY60[x] + _75GY60[x]) * 0.5f)  && (hue > (_75GY60[x] + _5GY60[x]) * 0.5f)) {
                    correction = _75GY60[y] - _75GY60[x] ;
                    lbe = _75GY60[y] ;
                    correctL = true;
                } else if ((hue <= (_5GY60[x] + _75GY60[x]) * 0.5f)  && (hue > (_5GY60[x] - 0.035f))) {
                    correction = _5GY60[y] - _5GY60[x] ;
                    lbe = _5GY60[y];
                    correctL = true;
                }
            } else if (lum < 75.f) {
                if ((hue <= (_7G70[x] + 0.035f)) && (hue > (_7G70[x] + _5G70[x]) * 0.5f)) {
                    correction = _7G70[y] - _7G70[x] ;
                    lbe = _7G70[y];
                    correctL = true;
                } else if ((hue <= (_7G70[x] + _5G70[x]) * 0.5f) && (hue > (_5G70[x] + _25G70[x]) * 0.5f)) {
                    correction = _5G70[y] - _5G70[x] ;
                    lbe = _5G70[y];
                    correctL = true;
                } else if ((hue <= (_25G70[x] + _5G70[x]) * 0.5f)  && (hue > (_25G70[x] + _1G70[x]) * 0.5f)) {
                    correction = _25G70[y] - _25G70[x] ;
                    lbe = _25G70[y];
                    correctL = true;
                } else if ((hue <= (_1G70[x] + _25G70[x]) * 0.5f)  && (hue > (_1G70[x] + _10GY70[x]) * 0.5f)) {
                    correction = _1G70[y] - _1G70[x] ;
                    lbe = _1G70[y] ;
                    correctL = true;
                } else if ((hue <= (_1G70[x] + _10GY70[x]) * 0.5f)  && (hue > (_10GY70[x] + _75GY70[x]) * 0.5f)) {
                    correction = _10GY70[y] - _10GY70[x] ;
                    lbe = _10GY70[y];
                    correctL = true;
                } else if ((hue <= (_10GY70[x] + _75GY70[x]) * 0.5f)  && (hue > (_75GY70[x] + _5GY70[x]) * 0.5f)) {
                    correction = _75GY70[y] - _75GY70[x] ;
                    lbe = _75GY70[y];
                    correctL = true;
                } else if ((hue <= (_5GY70[x] + _75GY70[x]) * 0.5f)  && (hue > (_5GY70[x] - 0.035f))) {
                    correction = _5GY70[y] - _5GY70[x] ;
                    lbe =  _5GY70[y];
                    correctL = true;
                }
            } else if (lum < 85.f) {
                if ((hue <= (_7G80[x] + 0.035f)) && (hue > (_7G80[x] + _5G80[x]) * 0.5f)) {
                    correction = _7G80[y] - _7G80[x] ;
                    lbe = _7G80[y];
                    correctL = true;
                } else if ((hue <= (_7G80[x] + _5G80[x]) * 0.5f) && (hue > (_5G80[x] + _25G80[x]) * 0.5f)) {
                    correction = _5G80[y] - _5G80[x] ;
                    lbe = _5G80[y];
                    correctL = true;
                } else if ((hue <= (_25G80[x] + _5G80[x]) * 0.5f)  && (hue > (_25G80[x] + _1G80[x]) * 0.5f)) {
                    correction = _25G80[y] - _25G80[x] ;
                    lbe = _25G80[y];
                    correctL = true;
                } else if ((hue <= (_1G80[x] + _25G80[x]) * 0.5f)  && (hue > (_1G80[x] + _10GY80[x]) * 0.5f)) {
                    correction = _1G80[y] - _1G80[x] ;
                    lbe = _1G80[y];
                    correctL = true;
                } else if ((hue <= (_1G80[x] + _10GY80[x]) * 0.5f)  && (hue > (_10GY80[x] + _75GY80[x]) * 0.5f)) {
                    correction = _10GY80[y] - _10GY80[x] ;
                    lbe = _10GY80[y];
                    correctL = true;
                } else if ((hue <= (_10GY80[x] + _75GY80[x]) * 0.5f)  && (hue > (_75GY80[x] + _5GY80[x]) * 0.5f)) {
                    correction = _75GY80[y] - _75GY80[x] ;
                    lbe = _75GY80[y];
                    correctL = true;
                } else if ((hue <= (_5GY80[x] + _75GY80[x]) * 0.5f)  && (hue > (_5GY80[x] - 0.035f))) {
                    correction = _5GY80[y] - _5GY80[x] ;
                    lbe = _5GY80[y];
                    correctL = true;
                }
            }
        }
    } else if (zone == 4) { //Red purple correction : only for L < 30
        if (lum > 5.f) {
            if (lum < 15.f && x < 45) {
                y = std::min(y, 44);
                if ((hue <= (_5R10[x] + 0.035f)) && (hue > (_5R10[x] - 0.043f))) {
                    correction = _5R10[y] - _5R10[x] ;
                    lbe = _5R10[y];
                    correctL = true;
                } else if ((hue <= (_25R10[x] + 0.043f)) && (hue > (_25R10[x] + _10RP10[x]) * 0.5f)) {
                    correction = _25R10[y] - _25R10[x] ;
                    lbe = _25R10[y];
                    correctL = true;
                } else if ((hue <= (_25R10[x] + _10RP10[x]) * 0.5f) && (hue > (_10RP10[x] - 0.035f))) {
                    correction = _10RP10[y] - _10RP10[x] ;
                    lbe = _10RP10[y];
                    correctL = true;
                }
            } else if (lum < 25.f && x < 70) {
                y = std::min(y, 70);
                if ((hue <= (_5R20[x] + 0.035f)) && (hue > (_5R20[x] + _25R20[x]) * 0.5f)) {
                    correction = _5R20[y] - _5R20[x] ;
                    lbe = _5R20[y];
                    correctL = true;
                } else if ((hue <= (_5R20[x] + _25R20[x]) * 0.5f) && (hue > (_10RP20[x] + _25R20[x]) * 0.5f)) {
                    correction = _25R20[y] - _25R20[x] ;
                    lbe = _25R20[y];
                    correctL = true;
                } else if ( (hue <= (_10RP20[x] + _25R20[x]) * 0.5f)  && (hue > (_10RP20[x] - 0.035f))) {
                    correction = _10RP20[y] - _10RP20[x] ;
                    lbe = _10RP20[y];
                    correctL = true;
                }
            } else if (lum < 35.f && x < 85) {
                y = rtengine::min(y, 85);
                if ((hue <= (_5R30[x] + 0.035f)) && (hue > (_5R30[x] + _25R30[x]) * 0.5f)) {
                    correction = _5R30[y] - _5R30[x] ;
                    lbe = _5R30[y];
                    correctL = true;
                } else if ((hue <= (_5R30[x] + _25R30[x]) * 0.5f) && (hue > (_10RP30[x] + _25R30[x]) * 0.5f)) {
                    correction = _25R30[y] - _25R30[x] ;
                    lbe = _25R30[y];
                    correctL = true;
                } else if ((hue <= (_10RP30[x] + _25R30[x]) * 0.5f)  && (hue > (_10RP30[x] - 0.035f))) {
                    correction = _10RP30[y] - _10RP30[x] ;
                    lbe = _10RP30[y];
                    correctL = true;
                }
            }
        }
    }
    //end red purple
}
}
