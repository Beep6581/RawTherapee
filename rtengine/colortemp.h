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
#ifndef _COLORTEMP_
#define _COLORTEMP_

#include <math.h>

namespace rtengine {

#define MINTEMP 1200
#define MAXTEMP 12000
#define MINGREEN 0.02
#define MAXGREEN 5.0

class ColorTemp {

    private:
        double temp;
        double green;

        static void clip (double &temp, double &green);
        
    public:
    
        ColorTemp () : temp(-1), green(-1) {}
        ColorTemp (double t, double g);
        ColorTemp (double mulr, double mulg, double mulb) { mul2temp (mulr, mulg, mulb, temp, green); }
        
        inline double getTemp ()    { return temp;  }
        inline double getGreen ()   { return green; }
        
        void   getMultipliers (double &mulr, double &mulg, double &mulb) { temp2mul (temp, green, mulr, mulg, mulb); }

  static void mul2temp (double rmul, double gmul, double bmul, double& temp, double& green);
  static void temp2mul (double temp, double green, double& rmul, double& gmul, double& bmul);

        bool operator== (const ColorTemp& other) { return fabs(temp-other.temp)<1e-10 && fabs(green-other.green)<1e-10; }
        bool operator!= (const ColorTemp& other) { return !(*this==other); }
};
};
#endif
