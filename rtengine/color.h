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

#ifndef _COLOR_H_
#define _COLOR_H_

#include "rt_math.h"
#include "LUT.h"
#include "labimage.h"
#include "iccstore.h"
#include "iccmatrices.h"
#include "sleef.c"
#define SAT(a,b,c) ((float)max(a,b,c)-(float)min(a,b,c))/(float)max(a,b,c)

namespace rtengine
{

#ifdef _DEBUG

class MunsellDebugInfo
{
public:
    float maxdhuelum[4];
    float maxdhue[4];
    unsigned int depass;
    unsigned int depassLum;

    MunsellDebugInfo();
    void reinitValues();
};

#endif

class Color
{

private:
    // Jacques' 195 LUTf for Munsell Lch correction
    static LUTf _4P10, _4P20, _4P30, _4P40, _4P50, _4P60;
    static LUTf _1P10, _1P20, _1P30, _1P40, _1P50, _1P60;
    static LUTf _5B40, _5B50, _5B60, _5B70, _5B80;
    static LUTf _7B40, _7B50, _7B60, _7B70, _7B80;
    static LUTf _9B40, _9B50, _9B60, _9B70, _9B80;
    static LUTf _10B40, _10B50, _10B60, _10B70, _10B80;
    static LUTf _05PB40, _05PB50, _05PB60, _05PB70, _05PB80;
    static LUTf _10PB10, _10PB20, _10PB30, _10PB40, _10PB50, _10PB60;
    static LUTf _9PB10, _9PB20, _9PB30, _9PB40, _9PB50, _9PB60, _9PB70, _9PB80;
    static LUTf _75PB10, _75PB20, _75PB30, _75PB40, _75PB50, _75PB60, _75PB70, _75PB80;
    static LUTf _6PB10, _6PB20, _6PB30, _6PB40, _6PB50, _6PB60, _6PB70, _6PB80;
    static LUTf _45PB10, _45PB20, _45PB30, _45PB40, _45PB50, _45PB60, _45PB70, _45PB80;
    static LUTf _3PB10, _3PB20, _3PB30, _3PB40, _3PB50, _3PB60, _3PB70, _3PB80;
    static LUTf _15PB10, _15PB20, _15PB30, _15PB40, _15PB50, _15PB60, _15PB70, _15PB80;
    static LUTf _10YR20, _10YR30, _10YR40, _10YR50, _10YR60, _10YR70, _10YR80, _10YR90;
    static LUTf _85YR20, _85YR30, _85YR40, _85YR50, _85YR60, _85YR70, _85YR80, _85YR90;
    static LUTf  _7YR30, _7YR40, _7YR50, _7YR60, _7YR70, _7YR80;
    static LUTf  _55YR30, _55YR40, _55YR50, _55YR60, _55YR70, _55YR80, _55YR90;
    static LUTf  _4YR30, _4YR40, _4YR50, _4YR60, _4YR70, _4YR80;
    static LUTf  _25YR30, _25YR40, _25YR50, _25YR60, _25YR70;
    static LUTf  _10R30, _10R40, _10R50, _10R60, _10R70;
    static LUTf  _9R30, _9R40, _9R50, _9R60, _9R70;
    static LUTf  _7R30, _7R40, _7R50, _7R60, _7R70;
    static LUTf  _5R10, _5R20, _5R30;
    static LUTf  _25R10, _25R20, _25R30;
    static LUTf  _10RP10, _10RP20, _10RP30;
    static LUTf  _7G30, _7G40, _7G50, _7G60, _7G70, _7G80;
    static LUTf  _5G30, _5G40, _5G50, _5G60, _5G70, _5G80;
    static LUTf  _25G30, _25G40, _25G50, _25G60, _25G70, _25G80;
    static LUTf  _1G30, _1G40, _1G50, _1G60, _1G70, _1G80;
    static LUTf  _10GY30, _10GY40, _10GY50, _10GY60, _10GY70, _10GY80;
    static LUTf  _75GY30, _75GY40, _75GY50, _75GY60, _75GY70, _75GY80;
    static LUTf  _5GY30, _5GY40, _5GY50, _5GY60, _5GY70, _5GY80;

    // Separated from init() to keep the code clear
    static void initMunsell ();
    static double hue2rgb(double p, double q, double t);
    static float hue2rgbfloat(float p, float q, float t);
#ifdef __SSE2__
    static vfloat hue2rgb(vfloat p, vfloat q, vfloat t);
#endif
public:

    typedef enum Channel {
        CHANNEL_RED            = 1 << 0,
        CHANNEL_GREEN          = 1 << 1,
        CHANNEL_BLUE           = 1 << 2,
        CHANNEL_HUE            = 1 << 3,
        CHANNEL_SATURATION     = 1 << 4,
        CHANNEL_VALUE          = 1 << 5,
        CHANNEL_LIGHTNESS      = 1 << 6,
        CHANNEL_CHROMATICITY   = 1 << 7
    } eChannel;

    typedef enum InterpolationPath {
        IP_SHORTEST, /// Interpolate color using the shortest path between 2 hues
        IP_LONGEST,  /// Interpolate color using the longest path between 2 hues
    } eInterpolationPath;

    typedef enum InterpolationDirection {
        ID_UP,       /// Interpolate color by increasing the hue value, crossing the upper limit
        ID_DOWN      /// Interpolate color by decreasing the hue value, crossing the lower limit
    } eInterpolationDirection;

    const static double sRGBGamma;        // standard average gamma
    const static double sRGBGammaCurve;   // 2.4 in the curve
    const static double eps, eps_max, kappa, epskap;
    const static float D50x, D50z;
    const static double u0, v0;

    static cmsToneCurve* linearGammaTRC;

    static LUTf cachef;
    static LUTf gamma2curve;

    // look-up tables for the standard srgb gamma and its inverse (filled by init())
    static LUTf igammatab_srgb;
    static LUTf igammatab_srgb1;
    static LUTf gammatab_srgb;
    static LUTf gammatab_srgb1;
    static LUTf igammatab_55;
    static LUTf gammatab_55;
    static LUTf igammatab_4;
    static LUTf gammatab_4;

    static LUTf igammatab_26_11;
    static LUTf gammatab_26_11;
    static LUTf igammatab_24_17;
    static LUTf gammatab_24_17a;
    static LUTf gammatab_13_2;
    static LUTf igammatab_13_2;
    static LUTf gammatab_115_2;
    static LUTf igammatab_115_2;
    static LUTf gammatab_145_3;
    static LUTf igammatab_145_3;

    // look-up tables for the simple exponential gamma
    static LUTf gammatab;
    static LUTuc gammatabThumb; // for thumbnails


    static void init ();
    static void cleanup ();


    /**
    * @brief Extract luminance "sRGB" from red/green/blue values
    * The range of the r, g and b channel has no importance ([0 ; 1] or [0 ; 65535]...) ; r,g,b can be negatives or > max, but must be in "sRGB"
    * @param r red channel
    * @param g green channel
    * @param b blue channel
    * @return luminance value
    */
    // xyz_sRGBD65 : conversion matrix from XYZ to sRGB for D65 illuminant: we use diagonal values
    static float rgbLuminance(float r, float g, float b)
    {
        // WArning: The sum of xyz_sRGBd65[1][] is > 1.0 (i.e. 1.0000001), so we use our own adapted values)
        // 0.2126729,  0.7151521,  0.0721750
        return r * 0.2126729f + g * 0.7151521f + b * 0.0721750f;
    }
    static double rgbLuminance(double r, double g, double b)
    {
        return r * 0.2126729 + g * 0.7151521 + b * 0.0721750;
    }


    /**
    * @brief Convert red/green/blue to hue/saturation/luminance
    * @param r red channel [0 ; 65535]
    * @param g green channel [0 ; 65535]
    * @param b blue channel [0 ; 65535]
    * @param h hue channel [0 ; 1] (return value)
    * @param s saturation channel [0 ; 1] (return value)
    * @param l luminance channel [0; 1] (return value)
    */
    static void rgb2hsl (float r, float g, float b, float &h, float &s, float &l);
    static void rgb2hslfloat (float r, float g, float b, float &h, float &s, float &l);
#ifdef __SSE2__
    static void rgb2hsl (vfloat r, vfloat g, vfloat b, vfloat &h, vfloat &s, vfloat &l);
#endif

    /**
    * @brief Convert hue/saturation/luminance in red/green/blue
    * @param h hue channel [0 ; 1]
    * @param s saturation channel [0 ; 1]
    * @param l luminance channel [0 ; 1]
    * @param r red channel [0 ; 65535] (return value)
    * @param g green channel [0 ; 65535] (return value)
    * @param b blue channel [0 ; 65535] (return value)
    */
    static void hsl2rgb (float h, float s, float l, float &r, float &g, float &b);
    static void hsl2rgbfloat (float h, float s, float l, float &r, float &g, float &b);
#ifdef __SSE2__
    static void hsl2rgb (vfloat h, vfloat s, vfloat l, vfloat &r, vfloat &g, vfloat &b);
#endif

    /**
    * @brief Convert hue/saturation/luminance in red/green/blue
    * @param h hue channel [0 ; 1]
    * @param s saturation channel [0 ; 1]
    * @param l luminance channel [0 ; 1]
    * @param r red channel [0 ; 1] (return value)
    * @param g green channel [0 ; 1] (return value)
    * @param b blue channel [0 ; 1] (return value)
    */
    static void hsl2rgb01 (float h, float s, float l, float &r, float &g, float &b);


    /**
    * @brief Convert red green blue to hue saturation value
    * @param r red channel [0 ; 65535]
    * @param g green channel [0 ; 65535]
    * @param b blue channel [0 ; 65535]
    * @param h hue channel [0 ; 1] (return value)
    * @param s saturation channel [0 ; 1] (return value)
    * @param v value channel [0 ; 1] (return value)
    */
    static void rgb2hsv (float r, float g, float b, float &h, float &s, float &v);

    static inline float rgb2s(float r, float g, float b) // fast version if only saturation is needed
    {
        float var_Min = min(r, g, b);
        float var_Max = max(r, g, b);
        float del_Max = var_Max - var_Min;

        if (del_Max < 0.00001f) {
            return 0.f;
        } else {
            return del_Max / var_Max;
        }
    }

    static inline bool rgb2hsvdcp(float r, float g, float b, float &h, float &s, float &v)
    {

        float var_Min = min(r, g, b);

        if(var_Min < 0.f) {
            return false;
        } else {
            float var_Max = max(r, g, b);
            float del_Max = var_Max - var_Min;
            v = var_Max / 65535.f;

            if (fabsf(del_Max) < 0.00001f) {
                h = 0.f;
                s = 0.f;
            } else {
                s = del_Max / var_Max;

                if ( r == var_Max ) {
                    h = (g - b) / del_Max;
                } else if ( g == var_Max ) {
                    h = 2.f + (b - r) / del_Max;
                } else { /*if ( b == var_Max ) */
                    h = 4.f + (r - g) / del_Max;
                }

                if ( h < 0.f ) {
                    h += 6.f;
                } else if ( h > 6.f ) {
                    h -= 6.f;
                }
            }

            return true;
        }
    }

    /**
    * @brief Convert hue saturation value in red green blue
    * @param h hue channel [0 ; 1]
    * @param s saturation channel [0 ; 1]
    * @param v value channel [0 ; 1]
    * @param r red channel [0 ; 65535] (return value)
    * @param g green channel [0 ; 65535] (return value)
    * @param b blue channel [0 ; 65535] (return value)
    */
    static void hsv2rgb (float h, float s, float v, float &r, float &g, float &b);

    static inline void hsv2rgbdcp (float h, float s, float v, float &r, float &g, float &b)
    {
        // special version for dcp which saves 1 division (in caller) and six multiplications (inside this function)
        int sector = h;  // sector 0 to 5, floor() is very slow, and h is always >0
        float f = h - sector; // fractional part of h

        v *= 65535.f;
        float vs = v * s;
        float p = v - vs;
        float q = v - f * vs;
        float t = p + v - q;

        switch (sector) {
            case 1:
                r = q;
                g = v;
                b = p;
                break;

            case 2:
                r = p;
                g = v;
                b = t;
                break;

            case 3:
                r = p;
                g = q;
                b = v;
                break;

            case 4:
                r = t;
                g = p;
                b = v;
                break;

            case 5:
                r = v;
                g = p;
                b = q;
                break;

            default:
                r = v;
                g = t;
                b = p;
        }
    }

    static void hsv2rgb (float h, float s, float v, int &r, int &g, int &b);


    /**
    * @brief Convert hue saturation value in red green blue
    * @param h hue channel [0 ; 1]
    * @param s saturation channel [0 ; 1]
    * @param v value channel [0 ; 1]
    * @param r red channel [0 ; 1] (return value)
    * @param g green channel [0 ; 1] (return value)
    * @param b blue channel [0 ; 1] (return value)
    */
    static void hsv2rgb01 (float h, float s, float v, float &r, float &g, float &b);


    /**
    * @brief Convert xyz to red/green/blue
    * Color space : sRGB   - illuminant D50 - use matrix sRGB_xyz[]
    * @param x X coordinate [0 ; 1] or [0 ; 65535]
    * @param y Y coordinate [0 ; 1] or [0 ; 65535]
    * @param z Z coordinate [0 ; 1] or [0 ; 65535]
    * @param r red channel [same range than xyz channel] (return value)
    * @param g green channel [same range than xyz channel] (return value)
    * @param b blue channel [same range than xyz channel] (return value)
    */
    static void xyz2srgb (float x, float y, float z, float &r, float &g, float &b);


    /**
    * @brief Convert xyz to red/green/blue
    * Color space : Prophoto   - illuminant D50 -  use the Prophoto_xyz[] matrix
    * @param x X coordinate [0 ; 1] or [0 ; 65535]
    * @param y Y coordinate [0 ; 1] or [0 ; 65535]
    * @param z Z coordinate [0 ; 1] or [0 ; 65535]
    * @param r red channel [same range than xyz channel] (return value)
    * @param g green channel [same range than xyz channel] (return value)
    * @param b blue channel [same range than xyz channel] (return value)
    */
    static void xyz2Prophoto (float x, float y, float z, float &r, float &g, float &b);


    /**
    * @brief Convert rgb in xyz
    * Color space : Prophoto   - illuminant D50 - use matrix xyz_prophoto[]
    * @param r red channel [0 ; 1] or [0 ; 65535] (return value)
    * @param g green channel [0 ; 1] or [0 ; 65535] (return value)
    * @param b blue channel [0 ; 1] or [0 ; 65535] (return value)
    * @param x X coordinate [same range than xyz channel]
    * @param y Y coordinate [same range than xyz channel]
    * @param z Z coordinate [same range than xyz channel]
    */
    static void Prophotoxyz (float r, float g, float b, float &x, float &y, float &z);


    /**
    * @brief Convert xyz in rgb
    * Color space : undefined - use adhoc matrix: rgb_xyz[3][3] (iccmatrice.h) in function of working space
    * @param x X coordinate [0 ; 1] or [0 ; 65535]
    * @param y Y coordinate [0 ; 1] or [0 ; 65535]
    * @param z Z coordinate [0 ; 1] or [0 ; 65535]
    * @param r red channel [same range than xyz channel] (return value)
    * @param g green channel [same range than xyz channel] (return value)
    * @param b blue channel [same range than xyz channel] (return value)
    * @param rgb_xyz[3][3] transformation matrix to use for the conversion
    */
    static void xyz2rgb (float x, float y, float z, float &r, float &g, float &b, const double rgb_xyz[3][3]);
    static void xyz2r (float x, float y, float z, float &r, const double rgb_xyz[3][3]);
    static void xyz2rgb (float x, float y, float z, float &r, float &g, float &b, const float rgb_xyz[3][3]);
#ifdef __SSE2__
    static void xyz2rgb (vfloat x, vfloat y, vfloat z, vfloat &r, vfloat &g, vfloat &b, const vfloat rgb_xyz[3][3]);
#endif


    /**
    * @brief Convert rgb in xyz
    * Color space : undefined - use adhoc matrix : xyz_rgb[3][3] (iccmatrice.h) in function of working space
    * @param r red channel [0 ; 1] or [0 ; 65535]
    * @param g green channel [0 ; 1] or [0 ; 65535]
    * @param b blue channel [0 ; 1] or [0 ; 65535]
    * @param x X coordinate [same range than rgb channel] (return value)
    * @param y Y coordinate [same range than rgb channel] (return value)
    * @param z Z coordinate [same range than rgb channel] (return value)
    * @param xyz_rgb[3][3] transformation matrix to use for the conversion
    */
    static void rgbxyz (float r, float g, float b, float &x, float &y, float &z, const double xyz_rgb[3][3]);
    static void rgbxyz (float r, float g, float b, float &x, float &y, float &z, const float xyz_rgb[3][3]);
#ifdef __SSE2__
    static void rgbxyz (vfloat r, vfloat g, vfloat b, vfloat &x, vfloat &y, vfloat &z, const vfloat xyz_rgb[3][3]);
#endif

    /**
    * @brief Convert Lab in xyz
    * @param L L channel [0 ; 32768] ; L can be negative rarely or superior 32768
    * @param a channel [-42000 ; +42000] ; can be more than 42000
    * @param b channel [-42000 ; +42000] ; can be more than 42000
    * @param x X coordinate [0 ; 65535] ; can be negative! (return value)
    * @param y Y coordinate [0 ; 65535] ; can be negative! (return value)
    * @param z Z coordinate [0 ; 65535] ; can be negative! (return value)
    */
    static void Lab2XYZ(float L, float a, float b, float &x, float &y, float &z);
    static void L2XYZ(float L, float &x, float &y, float &z);

#ifdef __SSE2__
    static void Lab2XYZ(vfloat L, vfloat a, vfloat b, vfloat &x, vfloat &y, vfloat &z);
#endif // __SSE2__

    /**
    * @brief Convert xyz in Lab
    * @param x X coordinate [0 ; 65535] ; can be negative or superior to 65535
    * @param y Y coordinate [0 ; 65535] ; can be negative or superior to 65535
    * @param z Z coordinate [0 ; 65535] ; can be negative or superior to 65535
    * @param L L channel [0 ; 32768] ; L can be negative rarely or superior 32768 (return value)
    * @param a channel [-42000 ; +42000] ; can be more than 42000 (return value)
    * @param b channel [-42000 ; +42000] ; can be more than 42000 (return value)
    */
    static void XYZ2Lab(float x, float y, float z, float &L, float &a, float &b);


    /**
    * @brief Convert Lab in Yuv
    * @param L L channel [0 ; 32768] ; L can be negative rarely or superior 32768
    * @param a channel [-42000 ; +42000] ; can be more than 42000
    * @param b channel [-42000 ; +42000] ; can be more than 42000
    * @param Y luminance channel [0 ; 65535] (return value)
    * @param u red chrominance channel [0 ; 65535] (return value)
    * @param v blue chrominance channel [0 ; 65535] (return value)
    */
    static void Lab2Yuv(float L, float a, float b, float &Y, float &u, float &v);


    /**
    * @brief Convert Yuv in Lab
    * @param Y luminance channel [0 ; 65535]
    * @param u red chrominance channel [0 ; 65535]
    * @param v blue chrominance channel [0 ; 65535]
    * @param L L channel [0 ; 32768] ; L can be negative rarely or superior 32768 (return value)
    * @param a channel [-42000 ; +42000] ; can be more than 42000 (return value)
    * @param b channel [-42000 ; +42000] ; can be more than 42000 (return value)
    */
    static void Yuv2Lab(float Y, float u, float v, float &L, float &a, float &b, double wp[3][3]);


    /**
    * @brief Convert the 'a' and 'b' channels of the L*a*b color space to 'c' and 'h' channels of the Lch color space (channel 'L' is identical [0 ; 32768])
    * @param a 'a' channel [-42000 ; +42000] ; can be more than 42000
    * @param b 'b' channel [-42000 ; +42000] ; can be more than 42000
    * @param c 'c' channel return value, in [0 ; 42000] ; can be more than 42000 (return value)
    * @param h 'h' channel return value, in [-PI ; +PI] (return value)
    */
    static void Lab2Lch(float a, float b, float &c, float &h);


    /**
    * @brief Convert 'c' and 'h' channels of the Lch color space to the 'a' and 'b' channels of the L*a*b color space (channel 'L' is identical [0 ; 32768])
    * @param c 'c' channel value, in [0 ; 42000]
    * @param h 'h' channel value, in [-PI ; +PI]
    * @param a 'a' channel [-42000 ; +42000] ; can be more than 42000 (return value)
    * @param b 'b' channel [-42000 ; +42000] ; can be more than 42000 (return value)
    */
    static void Lch2Lab(float c, float h, float &a, float &b);


    /**
    * @brief Convert the 'u' and 'v' channels of the Luv color space to 'c' and 'h' channels of the Lch color space ('L' channel is identical)
    * @param u 'u' channel [unknown range!]
    * @param v 'v' channel [unknown range!]
    * @param c 'c' channel [unknown range!] (return value)
    * @param h 'h' channel [-PI ; +PI] (return value)
    */
    static void Luv2Lch(float u, float v, float &c, float &h);


    /**
    * @brief Convert 'c' and 'h' channels of the Lch color space to the 'u' and 'v' channels of the Luv color space ('L' channel is identical)
    * @param c 'c' channel [unknown range!] ; can be more than 42000
    * @param h 'h' channel [-PI ; +PI]
    * @param u 'u' channel [unknown range!] (return value)
    * @param v 'v' channel [unknown range!] (return value)
    */
    static void Lch2Luv(float c, float h, float &u, float &v);


    /**
    * @brief Convert the XYZ values to Luv values
    * Warning: this method has never been used/tested so far
    * @param x X coordinate [0 ; 65535] ; can be negative or superior to 65535
    * @param y Y coordinate [0 ; 65535] ; can be negative or superior to 65535
    * @param z Z coordinate [0 ; 65535] ; can be negative or superior to 65535
    * @param L 'L' channel [0 ; 32768] (return value)
    * @param u 'u' channel [-42000 ; 42000] ; can be more than 42000 (return value)
    * @param v 'v' channel [-42000 ; 42000] ; can be more than 42000 (return value)
    */
    static void XYZ2Luv (float X, float Y, float Z, float &L, float &u, float &v);


    /**
    * @brief Convert the Luv values to XYZ values
    * Warning: this method has never been used/tested so far
    * @param L 'L' channel [0 ; 32768]
    * @param u 'u' channel [-42000 ; 42000] ; can be more than 42000
    * @param v 'v' channel [-42000 ; 42000] ; can be more than 42000
    * @param x X coordinate [0 ; 65535] ; can be negative or superior to 65535 (return value)
    * @param y Y coordinate [0 ; 65535] ; can be negative or superior to 65535 (return value)
    * @param z Z coordinate [0 ; 65535] ; can be negative or superior to 65535 (return value)
    */
    static void Luv2XYZ (float L, float u, float v, float &X, float &Y, float &Z);


    /**
    * @brief Return "f" in function of CIE's kappa and epsilon constants
    * @param f f can be fx fy fz where:
    *          fx=a/500 + fy  a=chroma green red [-128 ; +128]
    *          fy=(L+16)/116 L=luminance [0 ; 100]
    *          fz=fy-b/200 b=chroma blue yellow [-128 ; +128]
    */
    static inline double f2xyz(double f)
    {
        const double epsilonExpInv3 = 6.0 / 29.0;
        const double kappaInv = 27.0 / 24389.0; // inverse of kappa

        return (f > epsilonExpInv3) ? f * f * f : (116. * f - 16.) * kappaInv;

    }
    static inline float f2xyz(float f)
    {
        const float epsilonExpInv3 = 0.20689655f; // 6.0f/29.0f;
        const float kappaInv = 0.0011070565f; // 27.0f/24389.0f;  // inverse of kappa

        return (f > epsilonExpInv3) ? f * f * f : (116.f * f - 16.f) * kappaInv;
    }
#ifdef __SSE2__
    static inline vfloat f2xyz(vfloat f)
    {
        const vfloat epsilonExpInv3 = F2V(0.20689655f); // 6.0f/29.0f;
        const vfloat kappaInv = F2V(0.0011070565f); // 27.0f/24389.0f;  // inverse of kappa
        vfloat res1 = f * f * f;
        vfloat res2 = (F2V(116.f) * f - F2V(16.f)) * kappaInv;
        return vself(vmaskf_gt(f, epsilonExpInv3), res1, res2);
    }
#endif

    /**
     * @brief Calculate the effective direction (up or down) to linearly interpolating 2 colors so that it follows the shortest or longest path
     * @param h1 First hue [0 ; 1]
     * @param h2 Second hue [0 ; 1]
     * @param path Path to follow (shortest/longest)
     * @return The interpolation direction
     */
    static inline eInterpolationDirection getHueInterpolationDirection (double h1, double h2, eInterpolationPath path)
    {
        if (path == IP_SHORTEST) {
            if (h2 > h1) {
                if (h2 - h1 <= 0.5) {
                    return ID_UP;
                } else {
                    return ID_DOWN;
                }
            } else {
                if (h1 - h2 <= 0.5) {
                    return ID_DOWN;
                } else {
                    return ID_UP;
                }
            }
        } else {
            if (h2 > h1) {
                if (h2 - h1 <= 0.5) {
                    return ID_DOWN;
                } else {
                    return ID_UP;
                }
            } else {
                if (h1 - h2 <= 0.5) {
                    return ID_UP;
                } else {
                    return ID_DOWN;
                }
            }
        }
    }


    /**
     * @brief Calculate a color by linearly interpolating 2 colors
     * @param h1 First hue
     * @param h2 Second hue
     * @param balance Factor from 0 (first hue) to 1 (second hue)
     * @param dir Tells which direction the interpolation have to follow. You can get the value with getHueInterpolationDirection
     * @return The interpolated hue
     */
    static inline double interpolateHueHSV (double h1, double h2, double balance, eInterpolationDirection dir)
    {
        if (h1 == h2) {
            return h1;
        }

        if (dir == ID_DOWN) {
            if (h1 < h2) {
                double temp = h1;
                h1 = h2 - 1.;
                h2 = temp;
                balance = 1. - balance;
            }

            double h3 = h1 + balance * (h2 - h1);

            if (h3 < 0.) {
                h3 += 1.;
            }

            return h3;
        } else {
            if (h1 > h2) {
                h2 += 1.;
            }

            double h3 = h1 + balance * (h2 - h1);

            if (h3 > 1.) {
                h3 -= 1.;
            }

            return h3;
        }
    }

    /**
    * @brief Interpolate 2 colors from their respective red/green/blue channels, with a balance factor
    * @param balance gives weight to the first and second color [0 ; 1]
    *                0. = output color == first color
    *                0.5 = output color == equally mixed colors
    *                1. = output color == second color
    * @param r1 red channel of color 1 [0 ; 65535]
    * @param g1 green channel of color 1 [0 ; 65535]
    * @param b1 blue channel of color 1 [0 ; 65535]
    * @param r2 red channel of color 2 [0 ; 65535]
    * @param g2 green channel of color 2 [0 ; 65535]
    * @param b2 blue channel of color 2 [0 ; 65535]
    * @param channels bitfield of channel to interpolate (CHANNEL_LIGHTNESS|CHANNEL_CHROMATICITY|CHANNEL_HUE)
    * @param xyz_rgb color space
    * @param ro red channel of output color [0 ; 65535] (return value)
    * @param go green channel of output color [0 ; 65535] (return value)
    * @param bo blue channel of output color [0 ; 65535] (return value)
    */
    static void interpolateRGBColor (const float balance, const float r1, const float g1, const float b1, const float r2, const float g2, const float b2, int channels, const double xyz_rgb[3][3], const double rgb_xyz[3][3], float &ro, float &go, float &bo);

    /**
    * @brief Interpolate 2 colors from their respective red/green/blue channels, with a balance factor
    * @param realL luminance hsl [0; 1]
    * @param iplow low luminance for rl [0;1]
    * @param ihigh high luminance for r2 [0;1]
    * @param algm algorithm [0;2]
    * @param balance gives weight to the first and second color [0 ; 1]
    *                0. = output color == first color
    *                0.5 = output color == equally mixed colors
    *                1. = output color == second color
    * @param twoc 2 colors or 512 int
    * @param r1 red channel of color 1 [0 ; 65535]
    * @param g1 green channel of color 1 [0 ; 65535]
    * @param b1 blue channel of color 1 [0 ; 65535]
    * @param rl red channel of color low [0 ; 65535]
    * @param gl green channel of color low [0 ; 65535]
    * @param bl blue channel of color low [0 ; 65535]

    * @param r2 red channel of color 2 or high[0 ; 65535]
    * @param g2 green channel of color 2 or high[0 ; 65535]
    * @param b2 blue channel of color 2 [or high 0 ; 65535]
    * @param channels bitfield of channel to interpolate (CHANNEL_LIGHTNESS|CHANNEL_CHROMATICITY|CHANNEL_HUE)
    * @param xyz_rgb color space
    * @param rgb_xyz inverse color space
    * @param ro red channel of output color [0 ; 65535] (return value)
    * @param go green channel of output color [0 ; 65535] (return value)
    * @param bo blue channel of output color [0 ; 65535] (return value)
    */
    static void interpolateRGBColor (float realL, float iplow, float iphigh, int algm,  const float balance, int twoc, int metchrom, bool chr, bool lum, float chromat, float luma, const float r1, const float g1, const float b1, const float xl, const float yl, const float zl, const float x2, const float y2, const float z2, int channels, const double xyz_rgb[3][3], const double rgb_xyz[3][3], float &ro, float &go, float &bo);


    /**
    * @brief Interpolate a hue value as the angle of a polar coordinate with hue in the [0;1] range
    * Chose the shorter path from hue 1 to hue 2.
    * @param h1 First hue [0; 1]
    * @param h2 Second hue  [0; 1]
    * @param balance Interpolation factor [0 ; 1] where 0.=h1, 1.=h2
    * @return the interpolated value [0;1]
    */
    /*template <typename T, typename U>
    static inline T interpolatePolarHue_01 (T h1, T h2, U balance) {

        if (h1==h2)
            return h1;
        if ((h1 > h2) && (h1-h2 > T(0.5))){
            h1 -= T(1.);
            double value = h1 + T(balance) * (h2-h1);
            if (value < T(0.))
                value += T(1.);
            return value;
        }
        else if (h2-h1 > T(0.5)) {
            h2 -= T(1.);
            double value = h1 + T(balance) * (h2-h1);
            if (value < T(0.))
                value += T(1.);
            return value;
        }
        else
            return h1 + T(balance) * (h2-h1);
    }*/


    /**
    * @brief Interpolate a hue value as the angle of a polar coordinate with hue in the [-PI ; +PI] range
    * Chose the shorter path from hue 1 to hue 2.
    * @param h1 First hue [-PI ; +PI]
    * @param h2 Second hue [-PI ; +PI]
    * @param balance Interpolation factor [0 ; 1] where 0.=h1, 1.=h2
    * @return the interpolated value [-PI ; +PI]
    */
    /*template <typename T, typename U>
    static inline T interpolatePolarHue_PI (T h1, T h2, U balance) {
        if (h1==h2)
            return h1;
        if ((h1 > h2) && (h1-h2 > T(M_PI))){
            h1 -= T(2*M_PI);
            T value = h1 + T(balance) * (h2-h1);
            if (value < T(-M_PI))
                value += T(2*M_PI);
            return value;
        }
        else if (h2-h1 > T(M_PI)) {
            h2 -= T(2*M_PI);
            T value = h1 + T(balance) * (h2-h1);
            if (value < T(0))
                value += T(2*M_PI);
            return value;
        }
        else
            return h1 + T(balance) * (h2-h1);
    }*/


    /**
    * @brief Interpolate a hue value as the angle of a polar coordinate with hue in the [0;1] range
    * Chose the shorter path from hue 1 to hue 2.
    * @param h1 First hue [0; 1]
    * @param h2 Second hue  [0; 1]
    * @param balance Interpolation factor [0 ; 1] where 0.=h1, 1.=h2
    * @return the interpolated value [0;1]
    */
    template <typename T, typename U>
    static inline T interpolatePolarHue_01 (T h1, T h2, U balance)
    {
        float d = h2 - h1;
        float f;
        f = T(balance);
        double h;

        if (h1 > h2) {
            std::swap(h1, h2);
            d = -d;
            f = 1.f - f;
        }

        if (d < T(-M_PI) || d < T(0) || d > T(M_PI)) { //there was an inversion here !! d > T(M_PI)
            h1 += T(2 * M_PI);
            h = h1 + f * (h2 - h1);
            h = std::fmod(h, 2 * M_PI);
        } else {
            h = h1 + f * d;
        }

        // not strictly necessary..but in case of
        if(h < T(-M_PI)) {
            h = T(2 * M_PI) - h;
        }

        if(h > T(M_PI)) {
            h = h - T(2 * M_PI);
        }

        return h;
    }


    /**
    * @brief Interpolate a hue value as the angle of a polar coordinate with hue in the [-PI ; +PI] range
    * Chose the shorter path from hue 1 to hue 2.
    * @param h1 First hue [-PI ; +PI]
    * @param h2 Second hue [-PI ; +PI]
    * @param balance Interpolation factor [0 ; 1] where 0.=h1, 1.=h2
    * @return the interpolated value [-PI ; +PI ]
    */
    template <typename T, typename U>
    static inline T interpolatePolarHue_PI (T h1, T h2, U balance)
    {
        float d = h2 - h1;
        float f;
        f = T(balance);
        double h;

        if (h1 > h2) {
            std::swap(h1, h2);
            d = -d;
            f = 1.f - f;
        }

        if (d < T(0) || d < T(0.5) || d > T(1.)) { //there was an inversion here !! d > T(M_PI)
            h1 += T(1.);
            h = h1 + f * (h2 - h1);
            h = std::fmod(h, 1.);
        } else {
            h = h1 + f * d;
        }

        // not strictly necessary..but in case of
        if(h < T(0)) {
            h = T(1.) - h;
        }

        if(h > T(1)) {
            h = h - T(1.);
        }

        return h;
    }


    /**
    * @brief Get the gamma curves' parameters used by LCMS2
    * @param pwr gamma value [>1]
    * @param ts slope [0 ; 20]
    * @param mode [always 0]
    * @imax imax [always 0]
    * @param gamma0 used in ip2Lab2rgb [0 ; 1], usually near 0.5 (return value)
    * @param gamma1 used in ip2Lab2rgb [0 ; 20], can be superior to 20, but it's quite unusual(return value)
    * @param gamma2 used in ip2Lab2rgb [0 ; 1], usually near 0.03(return value)
    * @param gamma3 used in ip2Lab2rgb [0 ; 1], usually near 0.003(return value)
    * @param gamma4 used in ip2Lab2rgb [0 ; 1], usually near 0.03(return value)
    * @param gamma5 used in ip2Lab2rgb [0 ; 1], usually near 0.5 (return value)
    */
    static void calcGamma (double pwr, double ts, int mode, int imax, double &gamma0, double &gamma1, double &gamma2, double &gamma3, double &gamma4, double &gamma5);


    /**
    * @brief Used by Black and White to correct gamma for each channel red, green and blue channel
    * @param r red channel input and output value [0 ; 65535]
    * @param g green channel input and output value [0 ; 65535]
    * @param b blue channel input and output value [0 ; 65535]
    * @param gammabwr gamma value for red channel [>0]
    * @param gammabwg gamma value for red channel [>0]
    * @param gammabwb gamma value for red channel [>0]
    */
    static void trcGammaBW (float &r, float &g, float &b, float gammabwr, float gammabwg, float gammabwb);


    /** @brief Compute the B&W constants for the Black and White processing and its GUI
    * @param setting main mode
    * @param filter string of the filter effect to use
    * @param algo choice between linear and special for OYCPM colors
    * @param mixerRed red channel value of the channel mixer [-100 ; +200]
    * @param mixerGreen green channel value of the channel mixer [-100 ; +200]
    * @param mixerBlue blue channel value of the channel mixer [-100 ; +200]
    * @param mixerOrange orange channel value of the channel mixer [-100 ; +200]
    * @param mixerYellow yellow channel value of the channel mixer [-100 ; +200]
    * @param mixerCyan cyan channel value of the channel mixer [-100 ; +200]
    * @param mixerPurple purple channel value of the channel mixer [-100 ; +200]
    * @param mixerMagenta magenta channel value of the channel mixer [-100 ; +200]
    * @param autoc automatic mode of the channel mixer
    * @param complement adjust complementary channel
    * @param kcorec in absolute mode, value to correct the mixer [1 ; 3], usually near 1 (return value)
    * @param rrm red channel of the mixer (return value)
    * @param ggm green channel of the mixer (return value)
    * @param bbm blue channel of the mixer (return value)
    */
    static void computeBWMixerConstants (const Glib::ustring &setting, const Glib::ustring &filter, const Glib::ustring &algo, float &filcor, float &mixerRed, float &mixerGreen,
                                         float &mixerBlue, float mixerOrange, float mixerYellow, float mixerCyan, float mixerPurple, float mixerMagenta,
                                         bool autoc, bool complement, float &kcorec, double &rrm, double &ggm, double &bbm);


    // standard srgb gamma and its inverse

    /**
    * @brief sRGB gamma
    * See also calcGamma above with the following values: pwr=2.4  ts=12.92  mode=0.003041  imax=0.055011
    * @param x red, green or blue channel's value [0 ; 1]
    * @return the gamma modified's value [0 ; 1]
    */
    static inline double gamma2     (double x)      //  g3                  1+g4
    {
        return x <= 0.003041 ? x * 12.92 : 1.055011 * exp(log(x) / sRGBGammaCurve) - 0.055011;
    }


    /**
    * @brief Inverse sRGB gamma
    * See also calcGamma above with the following values: pwr=2.4  ts=12.92  mode=0.003041  imax=0.055011
    * @param x red, green or blue channel's value [0 ; 1]
    * @return the inverse gamma modified's value [0 ; 1]
    */
    static inline double igamma2    (double x)      //g2
    {
        return x <= 0.039293 ? x / 12.92 : exp(log((x + 0.055011) / 1.055011) * sRGBGammaCurve);
    }


    /**
    * @brief Get the gamma value for Gamma=5.5 Slope=10
    * @param x red, green or blue channel's value [0 ; 1]
    * @return the gamma modified's value [0 ; 1]
    */
    static inline double gamma55     (double x)     //  g3                  1+g4
    {
        return x <= 0.013189 ? x * 10.0 : 1.593503 * exp(log(x) / 5.5) - 0.593503; // 5.5 10
    }


    /**
    * @brief Get the inverse gamma value for Gamma=5.5 Slope=10
    * @param x red, green or blue channel's value [0 ; 1]
    * @return the inverse gamma modified's value [0 ; 1]
    */
    static inline double igamma55    (double x)     //g2
    {
        return x <= 0.131889 ? x / 10.0 : exp(log((x + 0.593503) / 1.593503) * 5.5); // 5.5 10
    }


    /**
    * @brief Get the gamma value for Gamma=4 Slope=5
    * @param x red, green or blue channel's value [0 ; 1]
    * @return the gamma modified's value [0 ; 1]
    */
    static inline double gamma4     (double x)      //  g3                  1+g4
    {
        return x <= 0.03089 ? x * 5.0 : 1.478793 * exp(log(x) / 4.1) - 0.478793; // 4  5
    }


    /**
    * @brief Get the inverse gamma value for Gamma=4 Slope=5
    * @param x red, green or blue channel's value [0 ; 1]
    * @return the inverse gamma modified's value [0 ; 1]
    */
    static inline double igamma4    (double x)      //g2
    {
        return x <= 0.154449 ? x / 5.0 : exp(log((x + 0.478793) / 1.478793) * 4.1); // 4 5
    }


    /*
    * @brief Get the gamma value for Gamma=2.2 Slope=4.5
    * @param x red, green or blue channel's value [0 ; 1]
    * @return the gamma modified's value [0 ; 1]
    *
    static inline double gamma709     (double x) {
                                            return x <= 0.0176 ? x*4.5 : 1.0954*exp(log(x)/2.2)-0.0954;
                                    }

    * @brief Get the inverse gamma value for Gamma=2.2 Slope=4.5
    * @param x red, green or blue channel's value [0 ; 1]
    * @return the inverse gamma modified's value [0 ; 1]
    *
    static inline double igamma709    (double x) {
                                        return x <= 0.0795 ? x/4.5 : exp(log((x+0.0954)/1.0954)*2.2);
                                    }
    */



    /**
    * @brief Get the gamma value for Gamma=2.4 Slope=17
    * @param x red, green or blue channel's value [0 ; 1]
    * @return the gamma modified's value [0 ; 1]
    */
    static inline double gamma24_17     (double x)
    {
        return x <= 0.001867 ? x * 17.0 : 1.044445 * exp(log(x) / 2.4) - 0.044445;
    }


    /**
    * @brief Get the inverse gamma value for Gamma=2.4 Slope=17
    * @param x red, green or blue channel's value [0 ; 1]
    * @return the inverse gamma modified's value [0 ; 1]
    */
    static inline double igamma24_17    (double x)
    {
        return x <= 0.031746 ? x / 17.0 : exp(log((x + 0.044445) / 1.044445) * 2.4);
    }


    /**
    * @brief Get the gamma value for Gamma=2.6 Slope=11
    * @param x red, green or blue channel's value [0 ; 1]
    * @return the gamma modified's value [0 ; 1]
    */
    static inline double gamma26_11     (double x)
    {
        return x <= 0.004921 ? x * 11.0 : 1.086603 * exp(log(x) / 2.6) - 0.086603;
    }


    /**
    * @brief Get the inverse gamma value for Gamma=2.6 Slope=11
    * @param x red, green or blue channel's value [0 ; 1]
    * @return the inverse gamma modified's value [0 ; 1]
    */
    static inline double igamma26_11    (double x)
    {
        return x <= 0.054127 ? x / 11.0 : exp(log((x + 0.086603) / 1.086603) * 2.6);
    }
    /**
    * @brief Get the gamma value for Gamma=1.3 Slope=2
    * @param x red, green or blue channel's value [0 ; 1]
    * @return the gamma modified's value [0 ; 1]
    */
    static inline double gamma13_2     (double x)
    {
        return x <= 0.016613 ? x * 2.0 : 1.009968 * exp(log(x) / 1.3) - 0.009968;
    }

    static inline double igamma13_2    (double x)
    {
        return x <= 0.033226 ? x / 2.0 : exp(log((x + 0.009968) / 1.009968) * 1.3);
    }

    static inline double gamma115_2     (double x)
    {
        return x <= 0.001692 ? x * 2.0 : 1.000508 * exp(log(x) / 1.15) - 0.000508;
    }

    static inline double igamma115_2    (double x)
    {
        return x <= 0.003384 ? x / 2.0 : exp(log((x + 0.000508) / 1.000508) * 1.15);
    }

    static inline double gamma145_3     (double x)
    {
        return x <= 0.009115 ? x * 3.0 : 1.012305 * exp(log(x) / 1.45) - 0.012305;
    }

    static inline double igamma145_3    (double x)
    {
        return x <= 0.027345 ? x / 3.0 : exp(log((x + 0.012305) / 1.012305) * 1.45);
    }

//gamma for Retinex
    static inline double gammareti      (double x, double gamma, double start, double slope, double mul, double add)
    {
        return (x <= start ? x*slope : exp(log(x) / gamma) * mul - add);
    }
    static inline double igammareti     (double x, double gamma, double start, double slope, double mul, double add)
    {
        return (x <= start * slope ? x / slope : exp(log((x + add) / mul) * gamma) );
    }



    // gamma function with adjustable parameters
    //same as above with values calculate with Calcgamma above
    // X range 0..1
    static inline double gamma      (double x, double gamma, double start, double slope, double mul, double add)
    {
        return (x <= start ? x*slope : exp(log(x) / gamma) * mul - add);
    }
    static inline double igamma     (double x, double gamma, double start, double slope, double mul, double add)
    {
        return (x <= start * slope ? x / slope : exp(log((x + add) / mul) * gamma) );
    }


    /**
    * @brief Very basic gamma
    * @param x red, green or blue channel's value [0 ; 1]
    * @param gamma gamma value [1 ; 5]
    * @return the gamma modified's value [0 ; 1]
    */
    static inline double gamman      (double x, double gamma)           //standard gamma without slope...
    {
        return (x = exp(log(x) / gamma));
    }

    /**
    * @brief Very basic gamma
    * @param x red, green or blue channel's value [0 ; 1]
    * @param gamma gamma value [1 ; 5]
    * @return the gamma modified's value [0 ; 1]
    */
    static inline float gammanf      (float x, float gamma)           //standard gamma without slope...
    {
        return (x = xexpf(xlogf(x) / gamma));
    }


    /**
    * @brief Very simply inverse gamma
    * @param x red, green or blue channel's value [0 ; 1]
    * @param gamma gamma value [1 ; 5]
    * @return the inverse gamma modified's value [0 ; 1]
    */
    static inline double igamman     (double x, double gamma)           //standard inverse gamma without slope...
    {
        return (x = exp(log(x) * gamma) );
    }


    /**
    * @brief Get the gamma value out of look-up tables
    * Calculated with gamma function above. e.g. :
    *    for (int i=0; i<65536; i++)
    *       gammatab_srgb[i] = (65535.0 * gamma2 (i/65535.0));
    * @param x [0 ; 1]
    * @return the gamma modified's value [0 ; 65535]
    */
    static inline float  gamma_srgb       (char x)
    {
        return gammatab_srgb[x];
    }
    static inline float  gamma            (char x)
    {
        return gammatab[x];
    }
    static inline float  igamma_srgb      (char x)
    {
        return igammatab_srgb[x];
    }
    static inline float  gamma_srgb       (int x)
    {
        return gammatab_srgb[x];
    }
    static inline float  gamma            (int x)
    {
        return gammatab[x];
    }
    static inline float  igamma_srgb      (int x)
    {
        return igammatab_srgb[x];
    }
    static inline float  gamma_srgb       (float x)
    {
        return gammatab_srgb[x];
    }
    static inline float  gamma_srgbclipped       (float x)
    {
        return gamma2curve[x];
    }
    static inline float  gamma            (float x)
    {
        return gammatab[x];
    }
    static inline float  igamma_srgb      (float x)
    {
        return igammatab_srgb[x];
    }
    //static inline float  gamma_srgb       (double x) { return gammatab_srgb[x]; }
    //static inline float  gamma            (double x) { return gammatab[x]; }
    //static inline float  igamma_srgb      (double x) { return igammatab_srgb[x]; }



    // --------------------------------  Jacques's Munsell correction


    /**
    * @brief Corrects the color (hue) depending on chromaticity and luminance changes
    *
    * To use in a "for" or "do while" statement.
    *
    * @param lumaMuns true => luminance correction (for delta L > 10) and chroma correction ; false => only chroma
    * @param Lprov1 luminance after [0 ; 100]
    * @param Loldd luminance before [0 ; 100]
    * @param HH hue before [-PI ; +PI]
    * @param Chprov1 chroma after [0 ; 180 (can be superior)]
    * @param CC chroma before [0 ; 180]
    * @param corectionHuechroma hue correction depending on chromaticity (saturation), in radians [0 ; 0.45] (return value)
    * @param correctlum hue correction depending on luminance (brightness, contrast,...), in radians [0 ; 0.45] (return value)
    * @param munsDbgInfo (Debug target only) object to collect informations
    */

#ifdef _DEBUG
    static void AllMunsellLch (bool lumaMuns, float Lprov1, float Loldd, float HH, float Chprov1, float CC, float &correctionHueChroma, float &correctlum, MunsellDebugInfo* munsDbgInfo);
#else
    static void AllMunsellLch (bool lumaMuns, float Lprov1, float Loldd, float HH, float Chprov1, float CC, float &correctionHueChroma, float &correctlum);
#endif


    /**
    * @brief Correct chromaticity and luminance so that the color stays in the working profile's gamut
    *
    * This function puts the data (Lab) in the gamut of "working profile":
    * it returns the corrected values of the chromaticity and luminance
    *
    * @param HH : hue, in radians [-PI ; +PI]
    * @param Lprov1 : input luminance value, sent back corrected [0 ; 100]  (input & output value)
    * @param Chprov1: input chroma value, sent back corrected [0 ; 180 (can be superior)]  (input & output value)
    * @param R red value of the corrected color [0 ; 65535 but can be negative or superior to 65535] (return value)
    * @param G green value of the corrected color [0 ; 65535 but can be negative or superior to 65535] (return value)
    * @param B blue value of the corrected color [0 ; 65535 but can be negative or superior to 65535] (return value)
    * @param wip working profile
    * @param isHLEnabled true if "Highlight Reconstruction " is enabled
    * @param lowerCoef a float number between [0.95 ; 1.0[
    *                  The nearest it is from 1.0, the more precise it will be, and the longer too as more iteration will be necessary
    * @param higherCoef a float number between [0.95 ; 1.0[
    *                   The nearest it is from 1.0, the more precise it will be, and the longer too as more iteration will be necessary
    * @param neg (Debug target only) to calculate iterations for negatives values
    * @param moreRGB (Debug target only) to calculate iterations for values >65535
    */
#ifdef _DEBUG
    static void gamutLchonly  (float HH, float &Lprov1, float &Chprov1, float &R, float &G, float &B, const double wip[3][3], const bool isHLEnabled, const float lowerCoef, const float higherCoef, bool &neg, bool &more_rgb);
    static void gamutLchonly  (float HH, float2 sincosval, float &Lprov1, float &Chprov1, float &R, float &G, float &B, const double wip[3][3], const bool isHLEnabled, const float lowerCoef, const float higherCoef, bool &neg, bool &more_rgb);
    static void gamutLchonly  (float2 sincosval, float &Lprov1, float &Chprov1, const float wip[3][3], const bool isHLEnabled, const float lowerCoef, const float higherCoef, bool &neg, bool &more_rgb);
#else
    static void gamutLchonly  (float HH, float &Lprov1, float &Chprov1, float &R, float &G, float &B, const double wip[3][3], const bool isHLEnabled, const float lowerCoef, const float higherCoef);
    static void gamutLchonly  (float HH, float2 sincosval, float &Lprov1, float &Chprov1, float &R, float &G, float &B, const double wip[3][3], const bool isHLEnabled, const float lowerCoef, const float higherCoef);
    static void gamutLchonly  (float2 sincosval, float &Lprov1, float &Chprov1, const float wip[3][3], const bool isHLEnabled, const float lowerCoef, const float higherCoef);
#endif


    /**
    * @brief Munsell gamut correction
    *
    * This function is the overall Munsell's corrections, but only on global statement. It may be better to use local statement with AllMunsellLch.
    * They are named accordingly :  gamutLchonly and AllMunsellLch
    * It can be used before and after treatment (saturation, gamma, luminance, ...)
    *
    * @param labL L channel input and output image
    *            L channel's usual range is [0 ; 100], but values can be negative or >100
    * @param laba a channel input and output image
    * @param labb b channel input and output image
    *            a and b channel's range is usually [-128 ; +128], but values can be >128
    * @param N Number of pixels to process
    * @param corMunsell performs Munsell correction
    * @param lumaMuns whether to apply luma correction or not (used only if corMuns=true)
    *                 true:  apply luma + chroma Munsell correction if delta L > 10;
    *                 false: leaves luma untouched
    * @param gamut performs gamutLch
    * @param wip matrix for working profile
    * @param multiThread whether to parallelize the loop or not
    */
    static void LabGamutMunsell (float *labL, float *laba, float *labb, const int N, bool corMunsell, bool lumaMuns, bool isHLEnabled, bool gamut, const double wip[3][3], bool multiThread );


    /*
    * @brief Skin tone protection factor
    * Skin colors: mixed from NX2 skin color palette, Von Luschan, and photos of white, black, yellow people...
    * There are some little exceptions, but it should cover 99% case.
    * Pay attention to white balance, and do not change hue and saturation, upstream of the modification
    * Used by vibrance
    * @param lum luma value [0 ; 100]
    * @param hue hue value [-PI ; +PI]
    * @param chrom chroma value [0 ; 180]
    * @param satreduc [0.1 ; 1] (return value)
    * @param chromx [0 or 1], actually only 0 is used
    */
    static void SkinSat         (float lum, float hue, float chrom, float &satreduc);//jacques Skin color


    /**
    * @brief Munsell Lch correction
    * Find the right LUT and calculate the correction
    * @param lum luma value [0 ; 100]
    * @param hue hue value [-PI ; +PI]
    * @param chrom chroma value [0 ; 180]
    * @param memChprov store chroma [0 ; 180]
    * @param correction correction value, in radians [0 ; 0.45]
    * @param lbe hue in function of chroma, in radian [-PI ; +PI]
    * @param zone  [1 ; 4]  1=PB correction + sky  2=red yellow correction 3=Green yellow correction  4=Red purple correction
    * @param correctL true=enable the Luminance correction
    */
    static void MunsellLch      (float lum, float hue, float chrom, float memChprov, float &correction, int zone, float &lbe, bool &correctL);//jacques:  Munsell correction


    // -------------------------------- end Munsell


    static void scalered ( const float rstprotection, const float param, const float limit, const float HH, const float deltaHH, float &scale, float &scaleext);
    static void transitred (const float HH, const float Chprov1, const float dred, const float factorskin, const float protect_red, const float factorskinext, const float deltaHH, const float factorsat, float &factor);
    static void skinred ( double J, double h, double sres, double Sp, float dred, float protect_red, int sk, float rstprotection, float ko, double &s);
    static void skinredfloat ( float J, float h, float sres, float Sp, float dred, float protect_red, int sk, float rstprotection, float ko, float &s);
//  static void scaleredcdbl ( float skinprot, float param, float limit, float HH, float deltaHH, float &scale,float &scaleext);

    static inline void SkinSatCbdl (float lum, float hue, float chrom, float skinprot, float &scale, bool neg, float b_l, float t_l, float t_r)
    {

        static const float C9 = 8.f, C8 = 15.f, C7 = 12.f, C4 = 7.f, C3 = 5.f, C2 = 5.f, C1 = 5.f;
        static const float H9 = 0.05f, H8 = 0.25f, H7 = 0.1f, H4 = 0.02f, H3 = 0.02f, H2 = 0.1f, H1 = 0.1f, H10 = -0.2f, H11 = -0.2f;

        // "real" skin color : take into account a slightly usage of contrast and saturation in RT if option "skin" = 1, uses imolicit factor 1.0
        // wide area  skin color, useful if not accurate colorimetry or if the user has changed hue and saturation, uses explicit facor 0.6
        // wide area for transition, uses explicit factor 0.4

        if  (lum >= 85.0f) {
            if((hue > (t_l + 0.53f - H9) && hue < (t_r + H9)) && (chrom > 8.0f && chrom < (14.0f + C9))) {
                scale = (100.f - skinprot) / 100.1f;
            } else if (lum >= 92.0f) {
                if((hue > t_l + 0.4f && hue < t_r) && (chrom > 7.0f && chrom < (15.0f))) {
                    scale = (100.f - skinprot * 0.6f) / 100.1f;
                } else if ((hue > b_l && hue < t_r) && (chrom > 7.0f && chrom < (18.0f))) {
                    scale = (100.f - skinprot * 0.4f) / 100.1f;
                }
            } else if ((hue > t_l + 0.4f && hue < t_r - 0.3f) && (chrom > 7.0f && chrom < (26.0f + C9))) {
                scale = (100.f - skinprot * 0.6f) / 100.1f;
            } else if ((hue > b_l + 0.05f && hue < t_r) && (chrom > 7.0f && chrom < (35.0f + C9))) {
                scale = (100.f - skinprot * 0.4f) / 100.1f;
            }
        } else if (lum >= 70.0f) {
            if((hue > t_l + 0.15f && hue < (t_r - 0.2f + H8)) && (chrom > 8.0f && chrom < (35.0f + C8))) {
                scale = (100.f - skinprot) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (48.0f + C9) )) {
                scale = (100.f - skinprot * 0.6f) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r) && (chrom > 7.0f && chrom < (55.0f + C9) )) {
                scale = (100.f - skinprot * 0.4f) / 100.1f;
            }
        } else if (lum >= 52.0f) {
            if((hue > t_l && hue < (t_r + H7)) && (chrom > 11.0f && chrom < (35.0f + C7))) {
                scale = (100.f - skinprot) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (48.0f + C9) )) {
                scale = (100.f - skinprot * 0.6f) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r) && (chrom > 7.0f && chrom < (55.0f + C9) )) {
                scale = (100.f - skinprot * 0.4f) / 100.1f;
            }
        } else if (lum >= 35.0f) {
            if((hue > t_l && hue <  (t_r + H4)) && (chrom > 13.0f && chrom < (37.0f + C4))) {
                scale = (100.f - skinprot) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (48.0f + C9) )) {
                scale = (100.f - skinprot * 0.6f) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r) && (chrom > 7.0f && chrom < (55.0f + C9) )) {
                scale = (100.f - skinprot * 0.4f) / 100.1f;
            }
        } else if (lum >= 20.0f) {
            if((hue > t_l && hue < (t_r + H3)) && (chrom > 7.0f && chrom < (35.0f + C3) )) {
                scale = (100.f - skinprot) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (48.0f + C9) )) {
                scale = (100.f - skinprot * 0.6f) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r) && (chrom > 7.0f && chrom < (55.0f + C9) )) {
                scale = (100.f - skinprot * 0.4f) / 100.1f;
            }
        } else if (lum >= 10.0f) {
            if((hue > (t_l - 0.25f + H10) && hue < (t_r - 0.3f + H2)) && (chrom > 8.0f && chrom < (23.0f + C2))) {
                scale = (100.f - skinprot) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (35.0f + C1) )) {
                scale = (100.f - skinprot * 0.6f) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.1f) && (chrom > 7.0f && chrom < (45.0f + C1) )) {
                scale = (100.f - skinprot * 0.4f) / 100.1f;
            }
        } else if ((hue > (t_l - 0.2f + H10) && hue < (t_r - 0.3f + H1)) && (chrom > 8.0f && chrom < (23.0f + C1))) {
            scale = (100.f - skinprot) / 100.1f;
        } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (35.0f + C1) )) {
            scale = (100.f - skinprot * 0.6f) / 100.1f;
        } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.1f) && (chrom > 7.0f && chrom < (45.0f + C1) )) {
            scale = (100.f - skinprot * 0.4f) / 100.1f;
        }

        //extended zone for hair, beard and if user adjust high value for skinprot
        if(skinprot > 85.f && chrom < 20.f && neg) {
            float modula = -0.0666f * skinprot + 6.66f;
            scale *= modula;
        }
    }

    static inline void SkinSatCbdl2 (float lum, float hue, float chrom, float skinprot, float &scale, bool neg, float b_l, float t_l, float t_r, float b_r, int basc)
    {

        static const float C9 = 8.f, C8 = 15.f, C7 = 12.f, C4 = 7.f, C3 = 5.f, C2 = 5.f, C1 = 5.f;
        static const float H9 = 0.05f, H8 = 0.25f, H7 = 0.1f, H4 = 0.02f, H3 = 0.02f, H2 = 0.1f, H1 = 0.1f, H10 = -0.2f, H11 = -0.2f;

        // "real" skin color : take into account a slightly usage of contrast and saturation in RT if option "skin" = 1, uses imolicit factor 1.0
        // wide area  skin color, useful if not accurate colorimetry or if the user has changed hue and saturation, uses explicit facor 0.6
        // wide area for transition, uses explicit factor 0.4
        if((b_l > -0.3f && b_r < 2.f)  || basc == 0) { //range maxi skin
            if  (lum >= 85.0f) {
                if((hue > (t_l + 0.53f - H9) && hue < (t_r + H9)) && (chrom > 8.0f && chrom < (14.0f + C9))) {
                    scale = (100.f - skinprot) / 100.1f;
                } else if (lum >= 92.0f) {
                    if((hue > t_l + 0.4f && hue < t_r) && (chrom > 7.0f && chrom < (15.0f))) {
                        scale = (100.f - skinprot * 0.6f) / 100.1f;
                    } else if ((hue > b_l && hue < t_r) && (chrom > 7.0f && chrom < (18.0f))) {
                        scale = (100.f - skinprot * 0.4f) / 100.1f;
                    }
                } else if ((hue > t_l + 0.4f && hue < t_r - 0.3f) && (chrom > 7.0f && chrom < (26.0f + C9))) {
                    scale = (100.f - skinprot * 0.6f) / 100.1f;
                } else if ((hue > b_l + 0.05f && hue < t_r) && (chrom > 7.0f && chrom < (35.0f + C9))) {
                    scale = (100.f - skinprot * 0.4f) / 100.1f;
                }
            } else if (lum >= 70.0f) {
                if((hue > t_l + 0.15f && hue < (t_r - 0.2f + H8)) && (chrom > 8.0f && chrom < (35.0f + C8))) {
                    scale = (100.f - skinprot) / 100.1f;
                } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (48.0f + C9) )) {
                    scale = (100.f - skinprot * 0.6f) / 100.1f;
                } else if ((hue > (b_l + 0.07f + H11) && hue < t_r) && (chrom > 7.0f && chrom < (55.0f + C9) )) {
                    scale = (100.f - skinprot * 0.4f) / 100.1f;
                }
            } else if (lum >= 52.0f) {
                if((hue > t_l && hue < (t_r + H7)) && (chrom > 11.0f && chrom < (35.0f + C7))) {
                    scale = (100.f - skinprot) / 100.1f;
                } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (48.0f + C9) )) {
                    scale = (100.f - skinprot * 0.6f) / 100.1f;
                } else if ((hue > (b_l + 0.07f + H11) && hue < t_r) && (chrom > 7.0f && chrom < (55.0f + C9) )) {
                    scale = (100.f - skinprot * 0.4f) / 100.1f;
                }
            } else if (lum >= 35.0f) {
                if((hue > t_l && hue <  (t_r + H4)) && (chrom > 13.0f && chrom < (37.0f + C4))) {
                    scale = (100.f - skinprot) / 100.1f;
                } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (48.0f + C9) )) {
                    scale = (100.f - skinprot * 0.6f) / 100.1f;
                } else if ((hue > (b_l + 0.07f + H11) && hue < t_r) && (chrom > 7.0f && chrom < (55.0f + C9) )) {
                    scale = (100.f - skinprot * 0.4f) / 100.1f;
                }
            } else if (lum >= 20.0f) {
                if((hue > t_l && hue < (t_r + H3)) && (chrom > 7.0f && chrom < (35.0f + C3) )) {
                    scale = (100.f - skinprot) / 100.1f;
                } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (48.0f + C9) )) {
                    scale = (100.f - skinprot * 0.6f) / 100.1f;
                } else if ((hue > (b_l + 0.07f + H11) && hue < t_r) && (chrom > 7.0f && chrom < (55.0f + C9) )) {
                    scale = (100.f - skinprot * 0.4f) / 100.1f;
                }
            } else if (lum >= 10.0f) {
                if((hue > (t_l - 0.25f + H10) && hue < (t_r - 0.3f + H2)) && (chrom > 8.0f && chrom < (23.0f + C2))) {
                    scale = (100.f - skinprot) / 100.1f;
                } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (35.0f + C1) )) {
                    scale = (100.f - skinprot * 0.6f) / 100.1f;
                } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.1f) && (chrom > 7.0f && chrom < (45.0f + C1) )) {
                    scale = (100.f - skinprot * 0.4f) / 100.1f;
                }
            } else if ((hue > (t_l - 0.2f + H10) && hue < (t_r - 0.3f + H1)) && (chrom > 8.0f && chrom < (23.0f + C1))) {
                scale = (100.f - skinprot) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (35.0f + C1) )) {
                scale = (100.f - skinprot * 0.6f) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.1f) && (chrom > 7.0f && chrom < (45.0f + C1) )) {
                scale = (100.f - skinprot * 0.4f) / 100.1f;
            }

            //extended zone for hair, beard and if user adjust high value for skinprot
            if(skinprot > 85.f && chrom < 20.f && neg) {
                float modula = -0.0666f * skinprot + 6.66f;
                scale *= modula;
            }
        }
        //end hue skin algo
        else if (basc == 1) { //not hue skin  linear transition or mod chroma curve
            if(hue >= t_l && hue <= t_r) {
                scale = (100.f - skinprot) / 100.1f;
            } else if(hue > b_l && hue < t_l) {
                float sc = (100.f - skinprot) / 100.1f;
                float aa = (1.f - sc) / (b_l - t_l);
                float bb = 1.f - aa * b_l;
                scale = aa * hue + bb;
            } else if(hue > t_r && hue < b_r) {
                float sc = (100.f - skinprot) / 100.1f;
                float aa = (sc - 1.f) / (t_r - b_r);
                float bb = 1.f - aa * b_r;
                scale = aa * hue + bb;
            }
        }
    }



    static inline void SkinSatCbdlCam (float lum, float hue, float chrom, float skinprot, float &scale, bool neg, float b_l, float t_l, float t_r)
    {

        static const float C9 = 8.f, C8 = 15.f, C7 = 12.f, C4 = 7.f, C3 = 5.f, C2 = 5.f, C1 = 5.f;
        static const float H9 = 0.05f, H8 = 0.25f, H7 = 0.1f, H4 = 0.02f, H3 = 0.02f, H2 = 0.1f, H1 = 0.1f, H10 = -0.2f, H11 = -0.2f;

        float HH = 0.f;

        if     (hue > 8.6f  && hue <= 74.f ) {
            HH = (1.15f / 65.4f) * hue - 0.0012f;   //H > 0.15   H<1.3
        } else if(hue > 0.f   && hue <= 8.6f ) {
            HH = (0.19f / 8.6f ) * hue - 0.04f;   //H>-0.04 H < 0.15
        } else if(hue > 355.f && hue <= 360.f) {
            HH = (0.11f / 5.0f ) * hue - 7.96f;   //H>-0.15 <-0.04
        } else if(hue > 74.f  && hue < 95.f  ) {
            HH = (0.30f / 21.0f) * hue + 0.24285f;   //H>1.3  H<1.6
        } else if(hue >= 95.f && hue < 137.5f) {
            HH = 0.01882f * hue - 0.18823f;   // H>1.6 H<2.4
        } else if(hue > 285.f && hue <= 355.f)  {
            HH = 0.1642f * hue - 5.982f;   //HH>-1.3  HH <-0.15
        }

        hue = HH;

        // "real" skin color : take into account a slightly usage of contrast and saturation in RT if option "skin" = 1, uses imolicit factor 1.0
        // wide area  skin color, useful if not accurate colorimetry or if the user has changed hue and saturation, uses explicit facor 0.6
        // wide area for transition, uses explicit factor 0.4

        if  (lum >= 85.0f) {
            if((hue > (t_l + 0.53f - H9) && hue < (t_r + H9)) && (chrom > 8.0f && chrom < (14.0f + C9))) {
                scale = (100.f - skinprot) / 100.1f;
            } else if (lum >= 92.0f) {
                if((hue > t_l + 0.4f && hue < t_r) && (chrom > 7.0f && chrom < (15.0f))) {
                    scale = (100.f - skinprot * 0.6f) / 100.1f;
                } else if ((hue > b_l && hue < t_r) && (chrom > 7.0f && chrom < (18.0f))) {
                    scale = (100.f - skinprot * 0.4f) / 100.1f;
                }
            } else if ((hue > t_l + 0.4f && hue < t_r - 0.3f) && (chrom > 7.0f && chrom < (26.0f + C9))) {
                scale = (100.f - skinprot * 0.6f) / 100.1f;
            } else if ((hue > b_l + 0.05f && hue < t_r) && (chrom > 7.0f && chrom < (35.0f + C9))) {
                scale = (100.f - skinprot * 0.4f) / 100.1f;
            }
        } else if (lum >= 70.0f) {
            if((hue > t_l + 0.15f && hue < (t_r - 0.2f + H8)) && (chrom > 8.0f && chrom < (35.0f + C8))) {
                scale = (100.f - skinprot) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (48.0f + C9) )) {
                scale = (100.f - skinprot * 0.6f) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r) && (chrom > 7.0f && chrom < (55.0f + C9) )) {
                scale = (100.f - skinprot * 0.4f) / 100.1f;
            }
        } else if (lum >= 52.0f) {
            if((hue > t_l && hue < (t_r + H7)) && (chrom > 11.0f && chrom < (35.0f + C7))) {
                scale = (100.f - skinprot) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (48.0f + C9) )) {
                scale = (100.f - skinprot * 0.6f) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r) && (chrom > 7.0f && chrom < (55.0f + C9) )) {
                scale = (100.f - skinprot * 0.4f) / 100.1f;
            }
        } else if (lum >= 35.0f) {
            if((hue > t_l && hue <  (t_r + H4)) && (chrom > 13.0f && chrom < (37.0f + C4))) {
                scale = (100.f - skinprot) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (48.0f + C9) )) {
                scale = (100.f - skinprot * 0.6f) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r) && (chrom > 7.0f && chrom < (55.0f + C9) )) {
                scale = (100.f - skinprot * 0.4f) / 100.1f;
            }
        } else if (lum >= 20.0f) {
            if((hue > t_l && hue < (t_r + H3)) && (chrom > 7.0f && chrom < (35.0f + C3) )) {
                scale = (100.f - skinprot) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (48.0f + C9) )) {
                scale = (100.f - skinprot * 0.6f) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r) && (chrom > 7.0f && chrom < (55.0f + C9) )) {
                scale = (100.f - skinprot * 0.4f) / 100.1f;
            }
        } else if (lum >= 10.0f) {
            if((hue > (t_l - 0.25f + H10) && hue < (t_r - 0.3f + H2)) && (chrom > 8.0f && chrom < (23.0f + C2))) {
                scale = (100.f - skinprot) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (35.0f + C1) )) {
                scale = (100.f - skinprot * 0.6f) / 100.1f;
            } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.1f) && (chrom > 7.0f && chrom < (45.0f + C1) )) {
                scale = (100.f - skinprot * 0.4f) / 100.1f;
            }
        } else if ((hue > (t_l - 0.2f + H10) && hue < (t_r - 0.3f + H1)) && (chrom > 8.0f && chrom < (23.0f + C1))) {
            scale = (100.f - skinprot) / 100.1f;
        } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.2f) && (chrom > 7.0f && chrom < (35.0f + C1) )) {
            scale = (100.f - skinprot * 0.6f) / 100.1f;
        } else if ((hue > (b_l + 0.07f + H11) && hue < t_r - 0.1f) && (chrom > 7.0f && chrom < (45.0f + C1) )) {
            scale = (100.f - skinprot * 0.4f) / 100.1f;
        }

        //extended zone for hair, beard and if user adjust high value for skinprot
        if(skinprot > 85.f && chrom < 20.f && neg) {
            float modula = -0.0666f * skinprot + 6.66f;
            scale *= modula;
        }
    }


    /**
    * @brief Gamut correction in the XYZ color space
    * @param X X channel input value and corrected output value [0 ; 65535]
    * @param Y Y channel input value and corrected output value [0 ; 65535]
    * @param Z Z channel input value and corrected output value [0 ; 65535]
    * @param p working profile
    */
    static void gamutmap(float &X, float &Y, float &Z, const double p[3][3]);


    /**
    * @brief Get HSV's hue from the Lab's hue
    * @param HH Lab's hue value, in radians [-PI ; +PI]
    * @return HSV's hue value [0 ; 1]
    */
    static inline double huelab_to_huehsv2 (float HH)
    {
        //hr=translate Hue Lab value  (-Pi +Pi) in approximative hr (hsv values) (0 1) [red 1/6 yellow 1/6 green 1/6 cyan 1/6 blue 1/6 magenta 1/6 ]
        // with multi linear correspondances (I expect there is no error !!)
        double hr = 0.0;
        //allways put h between 0 and 1

        if      (HH >= 0.f       && HH < 0.6f    ) {
            hr = 0.11666 * double(HH) + 0.93;    //hr 0.93  1.00    full red
        } else if (HH >= 0.6f      && HH < 1.4f    ) {
            hr = 0.1125 * double(HH) - 0.0675;    //hr 0.00  0.09    red yellow orange
        } else if (HH >= 1.4f      && HH < 2.f     ) {
            hr = 0.2666 * double(HH) - 0.2833;    //hr 0.09  0.25    orange yellow
        } else if (HH >= 2.f       && HH < 3.14159f) {
            hr = 0.1489 * double(HH) - 0.04785;    //hr 0.25  0.42    yellow green green
        } else if (HH >= -3.14159f && HH < -2.8f   ) {
            hr = 0.23419 * double(HH) + 1.1557;    //hr 0.42  0.50    green
        } else if (HH >= -2.8f     && HH < -2.3f   ) {
            hr = 0.16   * double(HH) + 0.948;    //hr 0.50  0.58    cyan
        } else if (HH >= -2.3f     && HH < -0.9f   ) {
            hr = 0.12143 * double(HH) + 0.85928;    //hr 0.58  0.75    blue blue-sky
        } else if (HH >= -0.9f     && HH < -0.1f   ) {
            hr = 0.2125 * double(HH) + 0.94125;    //hr 0.75  0.92    purple magenta
        } else if (HH >= -0.1f     && HH < 0.f     ) {
            hr = 0.1    * double(HH) + 0.93;    //hr 0.92  0.93    red
        }

        // in case of !
        if     (hr < 0.0) {
            hr += 1.0;
        } else if(hr > 1.0) {
            hr -= 1.0;
        }

        return (hr);
    }

};

}

#endif
