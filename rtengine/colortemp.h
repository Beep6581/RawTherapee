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

#include <gtkmm.h>
#include <cmath>

namespace rtengine {

#define MINTEMP 2000
#define MAXTEMP 25000
#define MINGREEN 0.02
#define MAXGREEN 5.0
#define INITIALBLACKBODY 4000

class ColorTemp {

    private:
        double temp;
        double green;
        Glib::ustring method;

        static void clip (double &temp, double &green);

    public:

        ColorTemp () : temp(-1), green(-1), method("Custom") {}
        ColorTemp (double t, double g, Glib::ustring m);
        ColorTemp (double mulr, double mulg, double mulb);

        inline double getTemp () const  { return temp;  }
        inline double getGreen () const { return green; }

        void   getMultipliers (double &mulr, double &mulg, double &mulb) { temp2mul (temp, green, mulr, mulg, mulb); }

        void mul2temp (double rmul, double gmul, double bmul, double& temp, double& green);
        void temp2mul (double temp, double green, double& rmul, double& gmul, double& bmul);
        //void temp2mul (double& rmul, double& gmul, double& bmul);

	bool operator== (const ColorTemp& other) const { return fabs(temp-other.temp)<1e-10 && fabs(green-other.green)<1e-10; }
        bool operator!= (const ColorTemp& other) const { return !(*this==other); }

        static double blackbody_spect        (double wavelength, double m1, double m2, double temp);
        static double daylight_spect         (double wavelength, double m1, double m2, double temp);
        static double Cloudy6200_spect       (double wavelength, double m1, double m2, double temp);
        static double Daylight5300_spect     (double wavelength, double m1, double m2, double temp);
        static double Shade7600_spect        (double wavelength, double m1, double m2, double temp);
        static double A2856_spect            (double wavelength, double m1, double m2, double temp);
        static double FluoF1_spect           (double wavelength, double m1, double m2, double temp);
        static double FluoF2_spect           (double wavelength, double m1, double m2, double temp);
        static double FluoF3_spect           (double wavelength, double m1, double m2, double temp);
        static double FluoF4_spect           (double wavelength, double m1, double m2, double temp);
        static double FluoF5_spect           (double wavelength, double m1, double m2, double temp);
        static double FluoF6_spect           (double wavelength, double m1, double m2, double temp);
        static double FluoF7_spect           (double wavelength, double m1, double m2, double temp);
        static double FluoF8_spect           (double wavelength, double m1, double m2, double temp);
        static double FluoF9_spect           (double wavelength, double m1, double m2, double temp);
        static double FluoF10_spect          (double wavelength, double m1, double m2, double temp);
        static double FluoF11_spect          (double wavelength, double m1, double m2, double temp);
        static double FluoF12_spect          (double wavelength, double m1, double m2, double temp);
        static double HMI_spect              (double wavelength, double m1, double m2, double temp);
        static double GTI_spect              (double wavelength, double m1, double m2, double temp);
        static double JudgeIII_spect         (double wavelength, double m1, double m2, double temp);
        static double Solux3500_spect        (double wavelength, double m1, double m2, double temp);
        static double Solux4100_spect        (double wavelength, double m1, double m2, double temp);
        static double Solux4700_spect        (double wavelength, double m1, double m2, double temp);
        static double NG_Solux4700_spect     (double wavelength, double m1, double m2, double temp);
        static double NG_LEDLSI2040_spect    (double wavelength, double m1, double m2, double temp);
        static double NG_CRSSP12WWMR16_spect (double wavelength, double m1, double m2, double temp);
        static double Flash5500_spect        (double wavelength, double m1, double m2, double temp);
        static double Flash6000_spect        (double wavelength, double m1, double m2, double temp);
        static double Flash6500_spect        (double wavelength, double m1, double m2, double temp);

        static void spectrum_to_xyz          (double (*spec_intens)(double wavelength, double m1, double m2, double temp), double _m1, double _m2, double _temp, double &x, double &y, double &z);
};
}
#endif
