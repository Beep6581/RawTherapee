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
#include <colortemp.h>

using namespace rtengine;

ColorTemp::ColorTemp (double t, double g) : temp(t), green(g) {

    clip (temp, green);
}

void ColorTemp::clip (double &temp, double &green) {

    if (temp < MINTEMP)
        temp = MINTEMP;
    else if (temp > MAXTEMP)
        temp = MAXTEMP;
        
    if (green < MINGREEN)
        green = MINGREEN;
    else if (green > MAXGREEN)
        green = MAXGREEN;
}

void ColorTemp::mul2temp (double rmul, double gmul, double bmul, double& temp, double& green) {

    double maxtemp=20000, mintemp=1000;
    double tmpr, tmpg, tmpb;
    temp=(maxtemp+mintemp)/2;
    while (maxtemp-mintemp>1) {
	temp2mul (temp, 1.0, tmpr, tmpg, tmpb);
	if (tmpb/tmpr > bmul/rmul)
	    maxtemp = temp;
	else
	    mintemp = temp;
        temp=(maxtemp+mintemp)/2;
    }
   green = (tmpg/tmpr) / (gmul/rmul);
   clip (temp, green);
}

void ColorTemp::temp2mul (double temp, double green, double& rmul, double& gmul, double& bmul) {

   clip (temp, green);

    double xD;
    if (temp<=4000) {
        xD = 0.27475e9/(temp*temp*temp) - 0.98598e6/(temp*temp) + 1.17444e3/temp + 0.145986;
    } else if (temp<=7000) {
	xD = -4.6070e9/(temp*temp*temp) + 2.9678e6/(temp*temp) + 0.09911e3/temp + 0.244063;
    } else {
	xD = -2.0064e9/(temp*temp*temp) + 1.9018e6/(temp*temp) + 0.24748e3/temp + 0.237040;
    }
    double yD = -3.0*xD*xD + 2.87*xD - 0.275;

    double X = xD/yD;
    double Y = 1.0;
    double Z = (1.0-xD-yD)/yD;

    rmul = X * 3.24071 - Y * 1.53726 - Z * 0.498571;
    gmul = - X * 0.969258 + Y * 1.87599 + Z * 0.0415557;
    bmul = X * 0.0556352 - Y * 0.203996 + Z * 1.05707;
    gmul /= green;

    double max = rmul;
    if (gmul>max) max = gmul;
    if (bmul>max) max = bmul;
    rmul /= max;
    gmul /= max;
    bmul /= max;
}
