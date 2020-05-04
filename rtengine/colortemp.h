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
#pragma once

#include <cmath>
#include <map>
#include <string>
#include "procparams.h"

namespace rtengine
{

constexpr double MINTEMP = 1500.0;
constexpr double MAXTEMP = 60000.0;
constexpr double MINGREEN = 0.02;
constexpr double MAXGREEN = 10.0;
constexpr double MINEQUAL = 0.8;
constexpr double MAXEQUAL = 1.5;
constexpr double INITIALBLACKBODY = 4000.0;


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
    const static std::map<std::string,const double *> spectMap;
public:

    ColorTemp () : temp(-1.), green(-1.), equal (1.), method("Custom") {}
    explicit ColorTemp (double e) : temp(-1.), green(-1.), equal (e), method("Custom") {}
    ColorTemp (double t, double g, double e, const std::string &m);
    ColorTemp (double mulr, double mulg, double mulb, double e);
    static void tempxy(bool separated, int repref, float **Tx, float **Ty, float **Tz, float **Ta, float **Tb, float **TL, double *TX, double *TY, double *TZ, const procparams::WBParams & wbpar);

    void update (const double rmul, const double gmul, const double bmul, const double equal, const double tempBias=0.0)
    {
        this->equal = equal;
        mul2temp (rmul, gmul, bmul, this->equal, temp, green);
        if (tempBias != 0.0 && tempBias >= -1.0 && tempBias <= 1.0) {
            temp += temp * tempBias;
        }
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
    static void temp2mulxyz (double tem, const std::string &method, double &Xxyz, double &Zxyz);

    static void cieCAT02(double Xw, double Yw, double Zw, double &CAM02BB00, double &CAM02BB01, double &CAM02BB02, double &CAM02BB10, double &CAM02BB11, double &CAM02BB12, double &CAM02BB20, double &CAM02BB21, double &CAM02BB22, double adap );
    static void cieCAT02float(float Xw, float Yw, float Zw, float &CAM02BB00, float &CAM02BB01, float &CAM02BB02, float &CAM02BB10, float &CAM02BB11, float &CAM02BB12, float &CAM02BB20, float &CAM02BB21, float &CAM02BB22, float adap);
    static void icieCAT02float(float Xw, float Yw, float Zw, float &iCAM02BB00, float &iCAM02BB01, float &iCAM02BB02, float &iCAM02BB10, float &iCAM02BB11, float &iCAM02BB12, float &iCAM02BB20, float &iCAM02BB21, float &iCAM02BB22, float adap);

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
    static const double JDC468_BluM5_spect[97]; //468 M5
    static const double JDC468_BluD6_spect[97]; //468 D6
    static const double JDC468_BluF4_spect[97]; //468 F4
    static const double JDC468_RedG21va_spect[97]; //468 G21 modifié
    static const double JDC468_RedI9_spect[97]; //468 I9
    static const double JDC468_GreI8_spect[97]; //468 I8
    static const double JDC468_GreQ7_spect[97]; //468 Q7
    static const double ColorGreenM25_spect[97];
    static const double ColorYellowkeltano_spect[97];
    static const double ColorGreenalsi_spect[97];
    static const double ColorRedpetunia_spect[97];
    static const double ColorRedkurttu_spect[97];
    static const double ColorRedlupiini_spect[97];
    static const double ColorRedhevosminttu_spect[97];
    static const double ColorRedneilikka_spect[97];
    static const double ColorRedpelagornia_spect[97];
    static const double ColorRedtalvio_spect[97];
    static const double ColorBrownpoimulehti_spect[97];
    static const double ColorOrangetuntematon_spect[97];
    static const double ColorOrangetlehmus_spect[97];
    static const double ColorOrangvaahtera_spect[97];
    static const double ColorBrownlehmus_spect[97];
    static const double ColorBrownuotiosammal_spect[97];
    static const double ColorBlacksoil_spect[97];
    static const double ColorGraynahjajaekaelae_spect[97];
    static const double ColorGreennuotisammal_spect[97];
    static const double ColorGreenleskenlehti_spect[97];
    static const double ColorGreenlinnunkaali_spect[97];
    static const double ColorGreenpelto_spect[97];
    static const double ColorGreenrodvoikukka[97];
    static const double ColorGreenlehmus[97];
    static const double ColorGreenlinden[97];
    static const double ColorYellowlehmus[97];
    static const double ColorYellowsuikeroalpi[97];
    static const double ColorYellowpensashanhikki1[97];
    static const double ColorYellowpensashanhikki2[97];
    static const double ColorBluehiidenvirna[97];
    static const double ColorBluekurkkuyrtti[97];
    static const double ColorPinksiankaersaemoe[97];
    static const double ColorVioletharakankello[97];
    static const double ColorVioletalsikeapila[97];
    static const double ColorVioletakilleija[97];
    static const double ColorOrangekehaekukka[97];
    static const double ColorRedpihlaja[97];
    static const double ColorVioletpetunia[97];
    static const double ColorVioletorvokki[97];
    static const double ColorBluesinisievikki[97];
    static const double ColorBlueiisoppi[97];
    static const double ColorBluelobelia[97];
    static const double ColorWhiteojaka[97];
    static const double ColorWhitepetunia[97];
    static const double ColorWhitepelargonia[97];
    static const double ColorWhitepaeivaen[97];
    static const double JDC468_B14_75Redspect[97];
    static const double Colorblue_spect[97];
    static const double ColorGreenkoriste[97];
    static const double ColorGreenpoimulehti[97];
    static const double ColorGreenhopeapaju[97];
    static const double ColorReduuden[97];
    static const double ColorRedpajuan[97];
    static const double ColorRedjaloan[97];
    static const double ColorBlueukon[97];
    static const double ColorBlueorvokki[97];
    static const double ColorBluemalvikki[97];
    static const double ColorBlackmaito[97];
    static const double ColorOrangpihlaja[97];
    static const double ColorBlackpihlaja[97];
    static const double ColorViolA1_spect[97];
    static const double ColorViolA4_spect[97];
    static const double ColorViolA6_spect[97];
    static const double ColorBlueSkyK3_spect[97];
    static const double ColorBlueSkyK9_spect[97];
    static const double ColorBlueSkyC4_spect[97];
    static const double ColorBlueSkyC14_spect[97];
    static const double ColorBlueSkyE4_spect[97];
    static const double ColorBlueSkyM1_spect[97];
    static const double ColorBlueSky2B1_spect[97];
    static const double ColorBlueSkyT7_spect[97];
    static const double ColorBlueSkyU19_spect[97];
    static const double ColorBlueSkyU2_spect[97];
    static const double ColorBlueSkyT17_spect[97];
    static const double ColorBlackM8_spect[97];
    static const double ColorBlackM12_spect[97];
    static const double ColorBlackM13_spect[97];
    static const double ColorWhite2B12_spect[97];
    static const double ColorWhite2B14_spect[97];
    static const double JDC468_Blackred97_spect[97];
    static const double JDC468_Blackredbl443_spect[97];
    static const double JDC468_Blackbl27_spect[97];
    static const double JDC468_Blackbl28_spect[97];
    static const double JDC468_Blackgr214_spect[97];
    static const double JDC468_Blackbl436_spect[97];
    static const double JDC468_Whitebl455_spect[97];
    static const double JDC468_Blackvio101_spect[97];
    static const double JDC468_Whitebl92_spect[97];
    static const double JDC468_Greyredbl94_spect[97];
    static const double JDC468_Blue32_spect[97];
    static const double JDC468_Blue236_spect[97];
    static const double JDC468_Gre300_spect[97];
    static const double JDC468_Blue340_spect[97];
    static const double JDC468_Gree110_spect[97];
    static const double JDC468_Gree457_spect[97];
    static const double JDC468_Yel241_spect[97];
    static const double JDC468_Ora321_spect[97];
    static const double JDC468_Yellow353_spect[97];
    static const double JDC468_Mag465_spect[97];
    static const double JDC468_Mag333_spect[97];
    static const double JDC468_Mag203_spect[97];


    static const double JDC468_OraO18_spect[97]; //468 O18
    static const double JDC468_OraD17_spect[97]; //468 D17
    static const double Fictif_61greyspect[97];//468 K15
    static const double JDC468_K15_87greyspect[97];
    static const double JDC468_YelN10_spect[97]; //468 N10
    static const double JDC468_GreN7_spect[97]; //468 N7
    static const double JDC468_GreA10_spect[97]; //468 A10
    static const double JDC468_GreK7_spect[97]; //468 K7
    static const double JDC468_PurE24_spect[97]; //468 E24
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
    static const double J570_BlueB6_spect[97];//blue Cyan
    static const double J570_BlueB15_spect[97];//blue Cyan
    static const double J570_BlueC2_spect[97];//blue Cyan
    static const double J570_BlueC14_spect[97];//blue Cyan
    static const double J570_BlueC16_spect[97];//blue Cyan
    static const double J570_BlueF1_spect[97];//blue Cyan
    static const double J570_BlueF2_spect[97];//blue Cyan
    static const double J570_BlueF10_spect[97];//blue Cyan
    static const double J570_BlueF13_spect[97];//blue Cyan
    static const double J570_BlueG9_spect[97];//blue Cyan
    static const double J570_BlueG19_spect[97];//blue Cyan
    static const double J570_BlueI5_spect[97];//blue Cyan
    static const double J570_BlueH15_spect[97];//blue Cyan
    static const double J570_BlueI3_spect[97];//blue Cyan
    static const double J570_BlueI19_spect[97];//blue Cyan
    static const double J570_BlueJ4_spect[97];//blue Cyan
    static const double J570_BlueJ6_spect[97];//blue Cyan
    static const double J570_BlueJ11_spect[97];//blue Cyan
    static const double J570_BlueJ13_spect[97];//blue Cyan
    static const double J570_BlueK5_spect[97];//blue Cyan
    static const double J570_BlueN1_spect[97];//blue Cyan
    static const double J570_BlueN4_spect[97];//blue Cyan
    static const double J570_BlueO19_spect[97];//blue Cyan
    static const double J570_BlueU8_spect[97];//blue Cyan
    static const double J570_NeuN8_spect[97];//neutral
    static const double J570_NeuN9_spect[97];//neutral
    static const double J570_NeuO8_spect[97];//neutral
    static const double J570_NeuO11_spect[97];//neutral
    static const double J570_NeuD5_spect[97];//neutral
    static const double J570_NeuE11_spect[97];//neutral
    static const double J570_NeuK16_spect[97];//neutral
    static const double J570_NeuM3_spect[97];//neutral
    static const double J570_NeuN18_spect[97];//neutral
    static const double J570_NeuQ1_spect[97];//neutral
    static const double J570_NeuS7_spect[97];//neutral
    static const double J570_NeuV10_spect[97];//neutral
    static const double J570_NeuW18_spect[97];//neutral
    static const double J570_NeuZ14_spect[97];//neutral
    static const double J570_NeuC18_spect[97];//neutral
    static const double J570_NeuD17_spect[97];//neutral
    static const double J570_NeuJ11_spect[97];//neutral
    static const double J570_NeuL4_spect[97];//neutral
    static const double Colorlab_n72_n2_spect[97];
    static const double Colorlab_10_n70_spect[97];
    static const double Colorlab_n33_n70_spect[97];
    static const double Colorlab_n8_n74_spect[97];
    static const double Colorlab_19_n69_spect[97];
    static const double Colorlab_n80_10_spect[97];
    static const double Colorlab_n80_26_spect[97];
    static const double Colorlab_n80_5_9_5_9spect[97];
//    static const double Colorlab_n57_5_6_9spect[97];
    
    /*
    static const double JDC468_greyc14_66_spect[97];
    static const double JDC468_greym13_325_spect[97];
    static const double JDC468_greyf26_156_spect[97];
    */
    static void spectrum_to_xyz_daylight  (double _m1, double _m2, double &x, double &y, double &z);
    static void spectrum_to_xyz_blackbody (double _temp, double &x, double &y, double &z);
    static void spectrum_to_xyz_preset    (const double* spec_intens, double &x, double &y, double &z);

    static void spectrum_to_color_xyz_daylight  (const double* spec_color, double _m1, double _m2, double &xx, double &yy, double &zz);
    static void spectrum_to_color_xyz_blackbody (const double* spec_color, double _temp, double &xx, double &yy, double &zz);
    static void spectrum_to_color_xyz_preset    (const double* spec_color, const double* spec_intens, double &xx, double &yy, double &zz);

};
}
