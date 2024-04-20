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

using color_match_type = double [97][3];

constexpr double MINTEMP = 1500.0;
constexpr double MAXTEMP = 60000.0;
constexpr double MINGREEN = 0.02;
constexpr double MAXGREEN = 100.0;
constexpr double MINEQUAL = 0.5;
constexpr double MAXEQUAL = 2.;
constexpr double INITIALBLACKBODY = 4000.0;

enum class StandardObserver {
    TWO_DEGREES,
    TEN_DEGREES,
};

class ColorTemp
{

private:
    double temp;
    double green;
    double equal;
    std::string method;
    StandardObserver observer{StandardObserver::TEN_DEGREES};
    static void clip (double &temp, double &green);
    static void clip (double &temp, double &green, double &equal);
    int XYZtoCorColorTemp(double x0, double y0 , double z0, double &temp) const;
    void temp2mul (double temp, double green, double equal, StandardObserver observer, double& rmul, double& gmul, double& bmul) const;
    const static std::map<std::string,const double *> spectMap;
public:
   // static constexpr StandardObserver DEFAULT_OBSERVER = StandardObserver::TEN_DEGREES;
    static constexpr StandardObserver DEFAULT_OBSERVER = StandardObserver::TWO_DEGREES;
    ColorTemp () : temp(-1.), green(-1.), equal (1.), method("Custom") {}
    explicit ColorTemp (double e) : temp(-1.), green(-1.), equal (e), method("Custom") {}
    ColorTemp (double t, double g, double e, const std::string &m, StandardObserver o);
    ColorTemp (double mulr, double mulg, double mulb, double e, StandardObserver observer);
    static void tempxy(bool separated, int repref, float **Tx, float **Ty, float **Tz, float **Ta, float **Tb, float **TL, double *TX, double *TY, double *TZ, const procparams::WBParams & wbpar, int ttbeg, int ttend, double &wpx, double &wpz, double *WPX, double *WPZ);

    void update (const double rmul, const double gmul, const double bmul, const double equal, StandardObserver observer, const double tempBias=0.0)
    {
        this->equal = equal;
        this->observer = observer;
        mul2temp (rmul, gmul, bmul, this->equal, observer, temp, green);
        if (tempBias != 0.0 && tempBias >= -1.0 && tempBias <= 1.0) {
            temp += temp * tempBias;
        }
    }
    void useDefaults (const double equal, StandardObserver observer)
    {
        temp = 6504;    // Values copied from procparams.cc
        green = 1.0;
        this->equal = equal;
        this->observer = observer;
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
    inline StandardObserver getObserver() const
    {
        return observer;
    }

    ColorTemp convertObserver(StandardObserver observer) const;

    void  getMultipliers (double &mulr, double &mulg, double &mulb) const
    {
        temp2mul (temp, green, equal, observer, mulr, mulg, mulb);
    }

    void mul2temp (const double rmul, const double gmul, const double bmul, const double equal, StandardObserver observer, double& temp, double& green) const;
    static void temp2mulxyz (double tem, const std::string &method, StandardObserver observer, double &Xxyz, double &Zxyz);

    static void cieCAT02(double Xw, double Yw, double Zw, double &CAM02BB00, double &CAM02BB01, double &CAM02BB02, double &CAM02BB10, double &CAM02BB11, double &CAM02BB12, double &CAM02BB20, double &CAM02BB21, double &CAM02BB22, double adap );
    static void cieCAT02float(float Xw, float Yw, float Zw, float &CAM02BB00, float &CAM02BB01, float &CAM02BB02, float &CAM02BB10, float &CAM02BB11, float &CAM02BB12, float &CAM02BB20, float &CAM02BB21, float &CAM02BB22, float adap);
    static void icieCAT02float(float Xw, float Yw, float Zw, float &iCAM02BB00, float &iCAM02BB01, float &iCAM02BB02, float &iCAM02BB10, float &iCAM02BB11, float &iCAM02BB12, float &iCAM02BB20, float &iCAM02BB21, float &iCAM02BB22, float adap);

    bool operator== (const ColorTemp& other) const
    {
        return fabs(temp - other.temp) < 1e-10 && fabs(green - other.green) < 1e-10 && observer != other.observer;
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
    static const double ColorchechCyaF3_spect2[97];
    static const double ColorchechCyaF3_spect3[97];
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

    static const double J570_NeuN8_spect2[97];//neutral
    static const double J570_NeuN9_spect2[97];//neutral
    static const double J570_NeuO8_spect2[97];//neutral
    static const double J570_NeuO11_spect2[97];//neutral
    static const double J570_NeuD5_spect2[97];//neutral
    static const double J570_NeuE11_spect2[97];//neutral
    static const double J570_NeuK16_spect2[97];//neutral
    static const double J570_NeuM3_spect2[97];//neutral
    static const double J570_NeuN18_spect2[97];//neutral
    static const double J570_NeuQ1_spect2[97];//neutral
    static const double J570_NeuS7_spect2[97];//neutral
    static const double J570_NeuV10_spect2[97];//neutral

    static const double J570_NeuW18_spect2[97];//neutral
    static const double J570_NeuZ14_spect2[97];//neutral
    static const double J570_NeuC18_spect2[97];//neutral
    static const double J570_NeuD17_spect2[97];//neutral
    static const double J570_NeuJ11_spect2[97];//neutral
    static const double J570_NeuL4_spect2[97];//neutral 
    
    static const double Colorlab_n72_n2_spect[97];
    static const double Colorlab_10_n70_spect[97];
    static const double Colorlab_10_n70_spect2[97];
    static const double Colorlab_10_n70_spect3[97];
    static const double Colorlab_10_n70_spect4[97];
    static const double Colorlab_n33_n70_spect[97];
    static const double Colorlab_n8_n74_spect[97];
    static const double Colorlab_19_n69_spect[97];
    static const double Colorlab_n80_10_spect[97];
    static const double Colorlab_n80_26_spect[97];
    static const double Colorlab_n80_5_9_5_9spect[97];
//    static const double JDC468_greyc14_66_spect[97];
//    static const double JDC468_greym13_325_spect[97];
//    static const double JDC468_greyf26_156_spect[97];
//    static const double Colorlab_n57_5_6_9spect[97];
    static const double Colorlab_L61_110_110Rec2020spect[97];
    static const double Colorlab_L63_120_m56Rec2020spect[97];
    static const double Colorlab_L63_m50_m60Rec2020spect[97];
    static const double Colorlab_L63_m120_80Rec2020spect[97];
    static const double Colorlab_L42_110_m100Prospect[97];
    static const double Colorlab_L42_m70_m100Prospect[97];
    static const double Colorlab_L56_m120_90Prospect[97];
    static const double Colorlab_L25_60_m120Prospect[97];
    static const double Colorlab_L75_50_120Prospect[97];
    static const double Colorlab_L75_m120_0Prospect[97];
    static const double Colorlab_L22_2_1_3Prospect[97];
    static const double Colorlab_L44_2_8_3_9spect[97];
    static const double Colorlab_L44_2_8_3_9spect2[97];
    static const double Colorlab_L95_2_3_15_6spect[97];
    static const double Colorlab_L95_2_3_15_6spect2[97];
    static const double Colorlab_L40_3_5_10_7spect[97];
    static const double Colorlab_L40_3_5_10_7spect2[97];
    static const double Colorlab_L40_3_5_10_7spect3[97];
    static const double Colorlab_L34_1_8_1_9spect[97];
    static const double Colorlab_L34_1_8_1_9spect2[97];
    static const double Colorlab_L64_1_8_m1_9spect[97];
    static const double Colorlab_L84_0_8_m1spect[97];
    static const double Colorlab_L63_1_3_m2spect[97];
    static const double Colorlab_L44_2_3_m3spect[97];
    static const double Colorlab_L65_96_45spect[97];
    static const double Colorlab_L52_47_57spect[97];
    static const double Colorlab_L31_62_27spect[97];
    static const double Colorlab_L79_m9_m28spect[97];
    static const double Colorlab_L58_50_31spect[97];
    static const double Colorlab_L31_m52_27spect[97];
    static const double Colorlab_L44_2_2_m7_35spect[97];
    static const double Colorlab_L47_m10_8_0_41spect[97];
    static const double Colorlab_L32_4_8_m3_2spect[97];
    static const double Colorlab_L57_m6_9_2_9spect[97];
    static const double Colorlab_L33_2_4_m4_5spect[97];
    static const double Colorlab_L35_11_65_m1_1spect[97];
    static const double Colorlab_L52_m2_7_8_9spect[97];
    static const double Colorlab_L32_7_m2_5spect[97];
    static const double Colorlab_L32_3_4_m3_8spect[97];
    static const double Colorlab_L50_m5_3_6_5spect[97];
    static const double Colorlab_L44_3_96_m8_8spect[97];
    static const double Colorlab_L34_3_6_m5_4spect[97];
    static const double Colorlab_L31_5_9_m4spect[97];
    static const double Colorlab_L35_3_4_m11spect[97];
    static const double Colorlab_L31_4_5_m4_7spect[97];
    static const double Colorlab_L35_4_8_m6_4spect[97];
    static const double Colorlab_L95_10_7_m14_3spect[97];
    static const double Colorlab_L36_34_m7_5spect[97];
    static const double Colorlab_L37_59_2spect[97];
    static const double Colorlab_L69_14_m9spect[97];
    static const double Colorlab_L92_13_m16spect[97];
    static const double Colorlab_L49_21_m12spect[97];
    static const double Colorlab_L56_20_m15spect[97];
    static const double Colorlab_L68_21_m19spect[97];
    static const double Colorlab_L98_m2_m32spect[97];
    static const double Colorlab_L98_m2_m32spect2[97];
    static const double Colorlab_L41_m27_m16spect[97];
    static const double Colorlab_L41_m27_m16spect2[97];
    static const double Colorlab_L15_m9_4spect[97];
    static const double Colorlab_L15_m9_4spect2[97];
    static const double Colorlab_L11_m11_2spect[97];
    static const double Colorlab_L14_m4_3spect[97];
    static const double Colorlab_L41_38_24spect[97];
    static const double Colorlab_L41_38_24spect2[97];
    static const double Colorlab_L53_48_58spect[97];
    static const double Colorlab_L53_48_58spect2[97];
    static const double Colorlab_L70_44_86spect[97];
    static const double Colorlab_L70_44_86spect2[97];
    static const double Colorlab_L38_42_19spect[97];
    static const double Colorlab_L38_42_19spect2[97];
    static const double Colorlab_L60_63_85spect[97];
    static const double Colorlab_L60_63_85spect2[97];
    static const double Colorlab_L80_75_30spect[97];
    static const double Colorlab_L80_75_30spect2[97];
    static const double Colorlab_L28_m21_24spect[97];
    static const double Colorlab_L28_m21_24spect2[97];
    static const double Colorlab_L45_m33_47spect[97];
    static const double Colorlab_L45_m33_47spect2[97];
    static const double Colorlab_L26_m7_404spect[97];
    static const double Colorlab_L34_m61_2spect[97];
    static const double Colorlab_L32_m16_17spect[97];
    static const double Colorlab_L30_m19_15spect[97];
    static const double Colorlab_L30_m17_16spect[97];
    static const double Colorlab_L35_m8_4spect[97];
    static const double Colorlab_L37_m7_5spect[97];
    static const double Colorlab_L45_m7_2spect[97];
    static const double Colorlab_L40_m6_5spect[97];
    static const double Colorlab_L46_m6_2spect[97];
    static const double Colorlab_L48_m69_16spect[97];
    static const double Colorlab_L51_89_53spect[97];
    static const double Colorlab_L49_84_33spect[97];
    static const double Colorlab_L59_m51_31spect[97];
    static const double Colorlab_L48_m69_16spect2[97];
    static const double Colorlab_L53_m71_6spect[97];
    static const double Colorlab_L51_m89_53spect2[97];
    static const double Colorlab_L49_84_33spect2[97];
    static const double Colorlab_L36_m27_28spect[97];
    static const double Colorlab_L36_m27_28spect2[97];
    static const double Colorlab_L36_m27_28spect3[97];
    static const double Colorlab_L63_16_71spect[97];
    static const double Colorlab_L84_4_46spect[97];
    static const double Colorlab_L84_4_46spect2[97];
    static const double Colorlab_L75_m66_19spect[97];
    static const double Colorlab_L75_m66_19spect2[97];
    static const double Colorlab_L64_m82_m6spect[97];
    static const double Colorlab_L64_m82_m6spect2[97];
    static const double Colorlab_L66_m71_m17spect[97];
    static const double Colorlab_L66_m71_m17spect2[97];
    static const double Colorlab_L22_m8_m60spect[97];
    static const double Colorlab_L22_m8_m60spect2[97];
    static const double Colorlab_L15_m4_m42spect[97];
    static const double Colorlab_L15_m4_m42spect2[97];
    static const double Colorlab_L13_3_m23spect[97];
    static const double Colorlab_L27_4_m90spect[97];
    static const double Colorlab_L19_1_m29spect[97];
    static const double Colorlab_L27_4_m90spect2[97];
    static const double Colorlab_L16_0_m44spect[97];
    static const double Colorlab_L16_0_m44spect2[97];
    static const double Colorlab_L13_m3_m36spect[97];
    static const double Colorlab_L13_m3_m36spect2[97];
    static const double Colorlab_L31_m23_m60spect[97];
    static const double Colorlab_L31_m23_m60spect2[97];
    static const double Colorlab_L17_3_m40spect[97];
    static const double Colorlab_L17_3_m40spect2[97];
    static const double Colorlab_L17_3_m40spect3[97];
    static const double Colorlab_L17_3_m40spect4[97];
    static const double Colorlab_L17_3_m40spect5[97];
    static const double Colorlab_L17_3_m40spect6[97];
    static const double Colorlab_L21_9_m7spect[97];
    static const double Colorlab_L78_4_m74spect[97];
    static const double Colorlab_L31_m58_m66spect[97];
    static const double Colorlab_L61_m11_m12spect[97];
    static const double Colorlab_L61_m11_m12spect2[97];
    static const double Colorlab_L29_1_m13spect[97];
    static const double Colorlab_L29_1_m13spect2[97];
    static const double Colorlab_L2_14_m1spect[97];
    static const double Colorlab_L5_39_m7spect[97];
    static const double Colorlab_L15_5_m13spect[97];
    static const double Colorlab_L12_5_m6spect[97];
    static const double Colorlab_L12_5_m6spect2[97];
    static const double Colorlab_L37_m59_m24spect[97];
    static const double Colorlab_L37_m59_m24spect2[97];
    static const double Colorlab_L15_55_23spect[97];
    static const double Colorlab_L11_m55_m11spect[97];
    static const double Colorlab_L8_m10_m2spect[97];
    static const double Colorlab_L14_m10_m7spect[97];
    static const double Colorlab_L20_m16_m13spect[97];
    static const double Colorlab_L8_m10_m2spect2[97];
    static const double Colorlab_L14_m10_m7spect2[97];
    static const double Colorlab_L20_m16_m13spect2[97];
    static const double Colorlab_L6_m9_1spect[97];
    static const double Colorlab_L20_m9_m10spect[97];
    static const double Colorlab_L85_10_45spect[97];
    static const double Colorlab_L90_m7_82spect[97];
    static const double Colorlab_L95_2_18spect[97];
    static const double Colorlab_L39_7_4spect[97];
    static const double Colorlab_L39_4_1spect[97];
    static const double Colorlab_L39_3_m1spect[97];
    static const double Colorlab_L40_3_m2spect[97];
    static const double Colorlab_L36_2_2spect[97];
    static const double Colorlab_L39_7_4spect2[97];
    static const double Colorlab_L39_4_1spect2[97];
    static const double Colorlab_L39_3_m1spect2[97];
    static const double Colorlab_L40_3_m2spect2[97];
    static const double Colorlab_L36_2_2spect2[97];
    static const double Colorlab_L40_4_m2spect[97];
    static const double Colorlab_L41_1_m6spect[97];
    static const double Colorlab_L40_4_m2spect2[97];
    static const double Colorlab_L41_1_m6spect2[97];
    static const double Colorlab_L41_12_14spect[97];
    static const double Colorlab_L41_12_14spect2[97];
    static const double Colorlab_L10_0_m22spect[97];
    static const double Colorlab_L38_60_8spect[97];
    static const double Colorlab_L49_85_39spect[97];
    static const double Colorlab_L42_1_m18spect[97];
    static const double Colorlab_L48_19_m25spect[97];
    static const double Colorlab_L30_21_m25spect[97];
    static const double Colorlab_L15_10_m15spect[97];
    static const double Colorlab_L48_19_m25spect2[97];
    static const double Colorlab_L30_21_m25spect2[97];
    static const double Colorlab_L15_10_m15spect2[97];
    static const double Colorlab_L60_26_m25spect[97];
    static const double Colorlab_L40_26_m45spect[97];
    static const double Colorlab_L40_26_m45spect2[97];
    static const double Colorlab_L20_10_m45spect[97];
    static const double Colorlab_L20_10_m45spect2[97];
    static const double Colorlab_L20_10_m45spect3[97];
    static const double ColorBlueSkyK3_spect2[97];
    static const double ColorBlueSkyK9_spect2[97];
    static const double ColorBlueSkyC4_spect2[97];
    static const double ColorBlueSkyC14_spect2[97];
    static const double ColorBlueSkyE4_spect2[97];
    static const double ColorBlueSkyM1_spect2[97];
    static const double ColorBlueSky2B1_spect2[97];
    static const double ColorBlueSkyT7_spect2[97]; 
    static const double ColorBlueSkyU19_spect2[97];
    static const double ColorBlueSkyU2_spect2[97];
    static const double Colorlab_L40_1_m40spect[97];
    static const double Colorlab_L30_4_m30spect[97];
    static const double Colorlab_L8_11_m25spect[97];
    static const double Colorlab_L40_1_m40spect2[97];
    static const double Colorlab_L30_4_m30spect2[97];
    static const double Colorlab_L8_11_m25spect2[97];
    static const double Colorlab_L26_m8_m25spect[97];
    static const double Colorlab_L26_m8_m25spect2[97];
    static const double Colorlab_L26_m8_m25spect3[97];
    static const double Colorlab_L22_1_m42spect[97];
    static const double Colorlab_L22_1_m42spect2[97];
    static const double Colorlab_L22_1_m42spect3[97];
    static const double Colorlab_L22_1_m42spect4[97];
    static const double Colorlab_L27_m1_m47spect[97];
    static const double Colorlab_L27_m1_m47spect2[97];
    static const double Colorlab_L40_30_m30spect[97];
    static const double Colorlab_L40_30_m30spect2[97];
    static const double Colorlab_L40_20_m35spect[97];
    static const double Colorlab_L40_20_m35spect2[97];
 
    static void spectrum_to_xyz_daylight  (double _m1, double _m2, double &x, double &y, double &z, const color_match_type &color_match);
    static void spectrum_to_xyz_blackbody (double _temp, double &x, double &y, double &z, const color_match_type &color_match);
    static void spectrum_to_xyz_preset    (const double* spec_intens, double &x, double &y, double &z, const color_match_type &color_match);

    static void spectrum_to_color_xyz_daylight  (const double* spec_color, double _m1, double _m2, double &xx, double &yy, double &zz, const color_match_type &color_match);
    static void spectrum_to_color_xyz_blackbody (const double* spec_color, double _temp, double &xx, double &yy, double &zz, const color_match_type &color_match);
    static void spectrum_to_color_xyz_preset    (const double* spec_color, const double* spec_intens, double &xx, double &yy, double &zz, const color_match_type &color_match);
    static void spectrum_to_whitepoint_xyz_daylight  (double _m1, double _m2, double &xx, double &yy, double &zz, const color_match_type &color_match);
    static void spectrum_to_whitepoint_xyz_blackbody (double _temp, double &xx, double &yy, double &zz, const color_match_type &color_match);
    static void whitepoint (double tempw, double &xx, double &yy, double &zz,const color_match_type &color_match);

};
}
