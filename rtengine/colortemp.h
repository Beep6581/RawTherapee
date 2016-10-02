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

#include <glibmm.h>
#include <cmath>

#define pow_F(a,b) (xexpf(b*xlogf(a)))

namespace rtengine
{

#define MINTEMP 1500
#define MAXTEMP 60000
#define MINGREEN 0.02
#define MAXGREEN 10.0
#define MINEQUAL 0.8
#define MAXEQUAL 1.5

#define INITIALBLACKBODY 4000


class ColorTemp
{

private:
    double temp;
    double green;
    double equal;
    std::string method;
    static void clip (double &temp, double &green);
    static void clip (double &temp, double &green, double &equal);
    int XYZtoCorColorTemp(double x0, double y0 , double z0, double &temp) const;
    void temp2mul (double temp, double green, double equal, double& rmul, double& gmul, double& bmul) const;

public:

    ColorTemp () : temp(-1.), green(-1.), equal (1.), method("Custom") {}
    explicit ColorTemp (double e) : temp(-1.), green(-1.), equal (e), method("Custom") {}
    ColorTemp (double t, double g, double e, const Glib::ustring &m);
    ColorTemp (double mulr, double mulg, double mulb, double e);

    void update (const double rmul, const double gmul, const double bmul, const double equal)
    {
        this->equal = equal;
        mul2temp (rmul, gmul, bmul, this->equal, temp, green);
    }
    void useDefaults (const double equal)
    {
        temp = 6504;    // Values copied from procparams.cc
        green = 1.0;
        this->equal = equal;
    }

    inline std::string getMethod() const
    {
        return method;
    }
    inline double getTemp () const
    {
        return temp;
    }
    inline double getGreen () const
    {
        return green;
    }
    inline double getEqual () const
    {
        return equal;
    }

    void  getMultipliers (double &mulr, double &mulg, double &mulb) const
    {
        temp2mul (temp, green, equal, mulr, mulg, mulb);
    }

    void mul2temp (const double rmul, const double gmul, const double bmul, const double equal, double& temp, double& green) const;
    static void temp2mulxyz (double tem, double gree, std::string method, double &Xxyz, double &Zxyz);

    static void cieCAT02(double Xw, double Yw, double Zw, double &CAM02BB00, double &CAM02BB01, double &CAM02BB02, double &CAM02BB10, double &CAM02BB11, double &CAM02BB12, double &CAM02BB20, double &CAM02BB21, double &CAM02BB22, double adap );
    //static    void CAT02 (Imagefloat* baseImg, const ProcParams* params);
    //static void ciecam_02 (LabImage* lab, const ProcParams* params);

    bool operator== (const ColorTemp& other) const
    {
        return fabs(temp - other.temp) < 1e-10 && fabs(green - other.green) < 1e-10;
    }
    bool operator!= (const ColorTemp& other) const
    {
        return !(*this == other);
    }

    static double blackbody_spect (double wavelength, double temperature);
    static double daylight_spect  (double wavelength, double m1, double m2);
    static const double Cloudy6200_spect[97];
    static const double Daylight5300_spect[97];
    static const double Shade7600_spect[97];
    static const double A2856_spect[97];
    static const double FluoF1_spect[97];
    static const double FluoF2_spect[97];
    static const double FluoF3_spect[97];
    static const double FluoF4_spect[97];
    static const double FluoF5_spect[97];
    static const double FluoF6_spect[97];
    static const double FluoF7_spect[97];
    static const double FluoF8_spect[97];
    static const double FluoF9_spect[97];
    static const double FluoF10_spect[97];
    static const double FluoF11_spect[97];
    static const double FluoF12_spect[97];
    static const double HMI_spect[97];
    static const double GTI_spect[97];
    static const double JudgeIII_spect[97];
    static const double Solux3500_spect[97];
    static const double Solux4100_spect[97];
    static const double Solux4700_spect[97];
    static const double NG_Solux4700_spect[97];
    static const double NG_LEDLSI2040_spect[97];
    static const double NG_CRSSP12WWMR16_spect[97];
    static const double Flash5500_spect[97];
    static const double Flash6000_spect[97];
    static const double Flash6500_spect[97];

    //spectral data 8 color Colorchecker24
    static double get_spectral_color (double wavelength, const double* array)
    {
        int wlm = (int) ((wavelength - 350.) / 5.);
        return (array[wlm]);
    }

    static const double ColorchechredC3_spect[97];
    static const double ColorchechOraA2_spect[97];
    static const double ColorchechYelD3_spect[97];
    static const double ColorchechGreE2_spect[97];
    static const double ColorchechGreB3_spect[97];
    static const double ColorchechCyaF3_spect[97];
    static const double ColorchechPurD2_spect[97];
    static const double ColorchechMagE3_spect[97];
    static const double ColorchechSkiA138_13_14_spect[97];
    static const double ColorchechGraC4_67_spect[97];
    static const double ColorchechSkiB166_18_18_spect[97];
    static const double ColorchechBluC150_m5_m22_spect[97];
    static const double ColorchechDCBluN881_m7_m14_spect[97];//ColorChecker DC  N8
    static const double ColorchechSGSkiF763_14_26_spect[97];//ColorChecker SG  F7
    static const double ColorchechSGSkiK285_11_17_spect[97];//ColorChecker SG  K2
    static const double ColorchechWhiA496_spect[97];
    static const double ColorchechGreD1_spect[97];
    static const double ColorchechSGBlaN3_6_spect[97];//ColorChecker SG  N3
    static const double JDC468_GraK14_44_spect[97];//468  K14
    static const double JDC468_BluH10_spect[97];//468  H10
    static const double ColabSkin35_15_17_spect[97];//Skin L 35
    static const double ColabSkin57_22_18_spect[97];//Skin L 57
    static const double ColabSkin40_17_17_spect[97];//Skin L 40
    static const double ColabSkin91_4_14_spect[97];//Skin L 91
    static const double ColabSkin87_8_8_spect[97];//Skin L 87
    static const double ColabSkin89_8_21_spect[97];//Skin L 89
    static const double ColabSkin75_8_4_spect[97];//Skin L 75
    static const double ColabSkin75_10_33_spect[97];//Skin L 75
    static const double ColabSkin65_33_11_spect[97];//Skin L 65
    static const double ColabSkin65_7_24_spect[97];//Skin L 65
    static const double ColabSkin57_19_6_spect[97];//Skin L 57
    static const double ColabSkin57_4_19_spect[97];//Skin L 57
    static const double ColabSkin57_10_28_spect[97];//Skin L 57
    static const double ColabSkin40_7_19_spect[97];//Skin L 57
    static const double ColabSkin40_17_6_spect[97];//Skin L 40
    static const double ColabSkin40_4_11_spect[97];//Skin L 40
    static const double ColabSkin33_6_15_spect[97];//Skin L 33
    static const double ColabSkin33_15_5_spect[97];//Skin L 33
    static const double ColabSkin33_10_15_spect[97];//Skin L 33
    static const double ColabSkin24_5_6_spect[97];//Skin L 24
    static const double ColabSkin26_18_18_spect[97];//Skin L 26
    static const double ColabSkin24_7_5_spect[97];//Skin L 24
    static const double ColabSkin20_4_2_spect[97];//Skin L 20
    static const double ColabSkin98_m2_10_spect[97];//Skin L 98
    static const double ColabSkin90_m1_20_spect[97];//Skin L 90
    static const double ColabSkin95_0_4_spect[97];//Skin L 95
    static const double ColabSkin81_2_14_spect[97];//Skin L 81
    static const double ColabSkin87_3_10_spect[97];//Skin L 87
    static const double ColabSkin77_12_21_spect[97];//Skin L 77
    static const double ColabSkin70_7_32_spect[97];//Skin L 77
    static const double ColabSky60_0_m31_spect[97];//Sky L=60
    static const double ColabSky42_0_m24_spect[97];//Sky L=42

    static void spectrum_to_xyz_daylight  (double _m1, double _m2, double &x, double &y, double &z);
    static void spectrum_to_xyz_blackbody (double _temp, double &x, double &y, double &z);
    static void spectrum_to_xyz_preset    (const double* spec_intens, double &x, double &y, double &z);

    static void spectrum_to_color_xyz_daylight  (const double* spec_color, double _m1, double _m2, double &xx, double &yy, double &zz);
    static void spectrum_to_color_xyz_blackbody (const double* spec_color, double _temp, double &xx, double &yy, double &zz);
    static void spectrum_to_color_xyz_preset    (const double* spec_color, const double* spec_intens, double &xx, double &yy, double &zz);

};
}
#endif
