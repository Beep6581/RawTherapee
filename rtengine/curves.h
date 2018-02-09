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
#ifndef __CURVES_H__
#define __CURVES_H__

#include <glibmm.h>
#include <map>
#include <string>
#include "rt_math.h"
#include "../rtgui/mycurve.h"
#include "../rtgui/myflatcurve.h"
#include "../rtgui/mydiagonalcurve.h"
#include "color.h"
#include "procparams.h"
#include "pipettebuffer.h"

#include "LUT.h"

#define CURVES_MIN_POLY_POINTS  1000

#include "rt_math.h"

#define CLIPI(a) ((a)>0?((a)<65534?(a):65534):0)

using namespace std;

namespace rtengine
{
class ToneCurve;
class ColorAppearance;

namespace curves {

inline void setLutVal(const LUTf &lut, float &val)
{
    if (!OOG(val)) {
        val = lut[std::max(val, 0.f)];
    } else {
        float m = lut[MAXVALF];
        val += (m - val);
    }
}

inline void setLutVal(float &val, float lutval, float maxval)
{
    if (!OOG(val)) {
        val = lutval;
    } else {
        val += (maxval - val);
    }
}

} // namespace curves

class CurveFactory
{

    friend class Curve;

protected:

    // functions calculating the parameters of the contrast curve based on the desired slope at the center
    static double solve_upper (double m, double c, double deriv);
    static double solve_lower (double m, double c, double deriv);
    static double dupper (const double b, const double m, const double c);
    static double dlower (const double b, const double m, const double c);

    // basic convex function between (0,0) and (1,1). m1 and m2 controls the slope at the start and end point
    static inline double basel (double x, double m1, double m2)
    {
        if (x == 0.0) {
            return 0.0;
        }

        double k = sqrt ((m1 - 1.0) * (m1 - m2) * 0.5) / (1.0 - m2);
        double l = (m1 - m2) / (1.0 - m2) + k;
        double lx = xlog(x);
        return m2 * x + (1.0 - m2) * (2.0 - xexp(k * lx)) * xexp(l * lx);
    }
    // basic concave function between (0,0) and (1,1). m1 and m2 controls the slope at the start and end point
    static inline double baseu (double x, double m1, double m2)
    {
        return 1.0 - basel(1.0 - x, m1, m2);
    }
    // convex curve between (0,0) and (1,1) with slope m at (0,0). hr controls the highlight recovery
    static inline double cupper (double x, double m, double hr)
    {
        if (hr > 1.0) {
            return baseu (x, m, 2.0 * (hr - 1.0) / m);
        }

        double x1 = (1.0 - hr) / m;
        double x2 = x1 + hr;

        if (x >= x2) {
            return 1.0;
        }

        if (x < x1) {
            return x * m;
        }

        return 1.0 - hr + hr * baseu((x - x1) / hr, m, 0);
    }
    // concave curve between (0,0) and (1,1) with slope m at (1,1). sr controls the shadow recovery
    static inline double clower (double x, double m, double sr)
    {
        return 1.0 - cupper(1.0 - x, m, sr);
    }
    // convex curve between (0,0) and (1,1) with slope m at (0,0). hr controls the highlight recovery
    static inline double cupper2 (double x, double m, double hr)
    {
        double x1 = (1.0 - hr) / m;
        double x2 = x1 + hr;

        if (x >= x2) {
            return 1.0;
        }

        if (x < x1) {
            return x * m;
        }

        return 1.0 - hr + hr * baseu((x - x1) / hr, m, 0.3 * hr);
    }
    static inline double clower2 (double x, double m, double sr)
    {
        //curve for b<0; starts with positive slope and then rolls over toward straight line to x=y=1
        double x1 = sr / 1.5 + 0.00001;

        if (x > x1 || sr < 0.001) {
            return 1 - (1 - x) * m;
        } else {
            double y1 = 1 - (1 - x1) * m;
            return y1 + m * (x - x1) - (1 - m) * SQR(SQR(1 - x / x1));
        }
    }
    // tone curve base. a: slope (from exp.comp.), b: black point normalized by 65535,
    // D: max. x value (can be>1), hr,sr: highlight,shadow recovery
    static inline double basecurve (double x, double a, double b, double D, double hr, double sr)
    {
        if (b < 0) {
            double m = 0.5;//midpoint
            double slope = 1.0 + b; //slope of straight line between (0,-b) and (1,1)
            double y = -b + m * slope; //value at midpoint

            if (x > m) {
                return y + (x - m) * slope;    //value on straight line between (m,y) and (1,1)
            } else {
                return y * clower2(x / m, slope * m / y, 2.0 - sr);
            }
        } else {
            double slope = a / (1.0 - b);
            double m = a * D > 1.0 ? b / a + (0.25) / slope : b + (1 - b) / 4;
            double y = a * D > 1.0 ? 0.25 : (m - b / a) * slope;

            if (x <= m) {
                return b == 0 ? x * slope : clower (x / m, slope * m / y, sr) * y;
            } else if (a * D > 1.0) {
                return y + (1.0 - y) * cupper2((x - m) / (D - m), slope * (D - m) / (1.0 - y), hr);
            } else {
                return y + (x - m) * slope;
            }
        }
    }
    static inline double simplebasecurve (double x, double b, double sr)
    {
        // a = 1, D = 1, hr = 0 (unused for a = D = 1)
        if (b == 0.0) {
            return x;
        } else if (b < 0) {
            double m = 0.5;//midpoint
            double slope = 1.0 + b; //slope of straight line between (0,-b) and (1,1)
            double y = -b + m * slope; //value at midpoint

            if (x > m) {
                return y + (x - m) * slope;    //value on straight line between (m,y) and (1,1)
            } else {
                return y * clower2(x / m, slope * m / y, 2.0 - sr);
            }
        } else {
            double slope = 1.0 / (1.0 - b);
            double m = b + (1 - b) * 0.25;
            double y = (m - b) * slope;

            if (x <= m) {
                return clower (x / m, slope * m / y, sr) * y;
            } else {
                return y + (x - m) * slope;
            }
        }
    }


public:
    const static double sRGBGamma;  // standard average gamma
    const static double sRGBGammaCurve;  // 2.4 in the curve


    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // accurately determine value from integer array with float as index
    //linearly interpolate from ends of range if arg is out of bounds
    static inline float interp(int *array, float f)
    {
        int index = CLIPI(floor(f));
        float part = (float)((f) - index) * (float)(array[index + 1] - array[index]);
        return (float)array[index] + part;
    }
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // accurately determine value from float array with float as index
    //linearly interpolate from ends of range if arg is out of bounds
    static inline float flinterp(float *array, float f)
    {
        int index = CLIPI(floor(f));
        float part = ((f) - (float)index) * (array[index + 1] - array[index]);
        return array[index] + part;
    }
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    static inline double centercontrast   (double x, double b, double m);

    // standard srgb gamma and its inverse
    static inline double gamma2            (double x)
    {
        return x <= 0.00304 ? x * 12.92 : 1.055 * exp(log(x) / sRGBGammaCurve) - 0.055;
    }
    static inline double igamma2           (double x)
    {
        return x <= 0.03928 ? x / 12.92 : exp(log((x + 0.055) / 1.055) * sRGBGammaCurve);
    }
    static inline float gamma2            (float x)
    {
        return x <= 0.00304 ? x * 12.92 : 1.055 * expf(logf(x) / sRGBGammaCurve) - 0.055;
    }
    static inline float igamma2           (float x)
    {
        return x <= 0.03928 ? x / 12.92 : expf(logf((x + 0.055) / 1.055) * sRGBGammaCurve);
    }
    // gamma function with adjustable parameters
    static inline double gamma            (double x, double gamma, double start, double slope, double mul, double add)
    {
        return (x <= start ? x*slope : exp(log(x) / gamma) * mul - add);
    }
    static inline double igamma           (double x, double gamma, double start, double slope, double mul, double add)
    {
        return (x <= start * slope ? x / slope : exp(log((x + add) / mul) * gamma) );
    }
    static inline float gamma            (float x, float gamma, float start, float slope, float mul, float add)
    {
        return (x <= start ? x*slope : xexpf(xlogf(x) / gamma) * mul - add);
    }
    static inline float igamma           (float x, float gamma, float start, float slope, float mul, float add)
    {
        return (x <= start * slope ? x / slope : xexpf(xlogf((x + add) / mul) * gamma) );
    }
#ifdef __SSE2__
    static inline vfloat igamma           (vfloat x, vfloat gamma, vfloat start, vfloat slope, vfloat mul, vfloat add)
    {
#if !defined(__clang__)
        return (x <= start * slope ? x / slope : xexpf(xlogf((x + add) / mul) * gamma) );
#else
        return vself(vmaskf_le(x, start * slope), x / slope, xexpf(xlogf((x + add) / mul) * gamma));
#endif
    }
#endif
    static inline float hlcurve (const float exp_scale, const float comp, const float hlrange, float level)
    {
        if (comp > 0.0) {
            float val = level + (hlrange - 65536.0);

            if(val == 0.0f) { // to avoid division by zero
                val = 0.000001f;
            }

            float Y = val * exp_scale / hlrange;
            Y *= comp;

            if(Y <= -1.0) { // to avoid log(<=0)
                Y = -.999999f;
            }

            float R = hlrange / (val * comp);
            return log1p(Y) * R;
        } else {
            return exp_scale;
        }
    }

public:
    static void complexCurve (double ecomp, double black, double hlcompr, double hlcomprthresh, double shcompr, double br, double contr,
                              const std::vector<double>& curvePoints, const std::vector<double>& curvePoints2,
                              LUTu & histogram, LUTf & hlCurve, LUTf & shCurve, LUTf & outCurve, LUTu & outBeforeCCurveHistogram, ToneCurve & outToneCurve, ToneCurve & outToneCurve2,

                              int skip = 1);
    static void curveBW (const std::vector<double>& curvePointsbw, const std::vector<double>& curvePointsbw2, const LUTu & histogrambw, LUTu & outBeforeCCurveHistogrambw,
                         ToneCurve & customToneCurvebw1, ToneCurve & customToneCurvebw2, int skip);

    static void curveCL ( bool & clcutili, const std::vector<double>& clcurvePoints, LUTf & clCurve, int skip);

    static void curveWavContL ( bool & wavcontlutili, const std::vector<double>& wavclcurvePoints, LUTf & wavclCurve,/* LUTu & histogramwavcl, LUTu & outBeforeWavCLurveHistogram,*/int skip);
    static void curveDehaContL ( bool & dehacontlutili, const std::vector<double>& dehaclcurvePoints, LUTf & dehaclCurve, int skip, const LUTu & histogram, LUTu & outBeforeCurveHistogram);
    static void mapcurve ( bool & mapcontlutili, const std::vector<double>& mapcurvePoints, LUTf & mapcurve, int skip, const LUTu & histogram, LUTu & outBeforeCurveHistogram);

    static void curveToning ( const std::vector<double>& curvePoints, LUTf & ToningCurve, int skip);

    static void complexsgnCurve ( bool & autili,  bool & butili, bool & ccutili, bool & clcutili, const std::vector<double>& acurvePoints,
                                  const std::vector<double>& bcurvePoints, const std::vector<double>& cccurvePoints, const std::vector<double>& lccurvePoints, LUTf & aoutCurve, LUTf & boutCurve, LUTf & satCurve, LUTf & lhskCurve,
                                  int skip = 1);
    static void complexLCurve (double br, double contr, const std::vector<double>& curvePoints, const LUTu & histogram, LUTf & outCurve, LUTu & outBeforeCCurveHistogram, int skip, bool & utili);

    static void curveLightBrightColor (
        const std::vector<double>& curvePoints,
        const std::vector<double>& curvePoints2,
        const std::vector<double>& curvePoints3,
        const LUTu & histogram, LUTu & outBeforeCCurveHistogram,
        const LUTu & histogramC, LUTu & outBeforeCCurveHistogramC,
        ColorAppearance & outColCurve1,
        ColorAppearance & outColCurve2,
        ColorAppearance & outColCurve3,
        int skip = 1);
    static void RGBCurve (const std::vector<double>& curvePoints, LUTf & outCurve, int skip);

};

class Curve
{

    class HashEntry
    {
    public:
        unsigned short smallerValue;
        unsigned short higherValue;
    };
protected:
    int N;
    int ppn;            // targeted polyline point number
    double* x;
    double* y;
    // begin of variables used in Parametric curves only
    double mc;
    double mfc;
    double msc;
    double mhc;
    // end of variables used in Parametric curves only
    std::vector<double> poly_x;     // X points of the faceted curve
    std::vector<double> poly_y;     // Y points of the faceted curve
    std::vector<double> dyByDx;
    std::vector<HashEntry> hash;
    unsigned short hashSize;        // hash table's size, between [10, 100, 1000]

    double* ypp;

    // Fields for the elementary curve polygonisation
    double x1, y1, x2, y2, x3, y3;
    bool firstPointIncluded;
    double increment;
    int nbr_points;

    static inline double p00 (double x, double prot)
    {
        return CurveFactory::clower (x, 2.0, prot);
    }
    static inline double p11 (double x, double prot)
    {
        return CurveFactory::cupper (x, 2.0, prot);
    }
    static inline double p01 (double x, double prot)
    {
        return x <= 0.5 ? CurveFactory::clower (x * 2, 2.0, prot) * 0.5 : 0.5 + CurveFactory::cupper ((x - 0.5) * 2, 2.0, prot) * 0.5;
    }
    static inline double p10 (double x, double prot)
    {
        return x <= 0.5 ? CurveFactory::cupper (x * 2, 2.0, prot) * 0.5 : 0.5 + CurveFactory::clower ((x - 0.5) * 2, 2.0, prot) * 0.5;
    }
    static inline double pfull (double x, double prot, double sh, double hl)
    {
        return (1 - sh) * (1 - hl) * p00(x, prot) + sh * hl * p11(x, prot) + (1 - sh) * hl * p01(x, prot) + sh * (1 - hl) * p10(x, prot);
    }

    void fillHash();
    void fillDyByDx();

public:
    Curve ();
    virtual ~Curve () {};
    void AddPolygons ();
    int getSize () const; // return the number of control points
    void getControlPoint(int cpNum, double &x, double &y) const;
    virtual double getVal (double t) const = 0;
    virtual void   getVal (const std::vector<double>& t, std::vector<double>& res) const = 0;

    virtual bool   isIdentity () const = 0;
};

class DiagonalCurve : public Curve
{

protected:
    DiagonalCurveType kind;

    void spline_cubic_set ();
    void NURBS_set ();

public:
    DiagonalCurve (const std::vector<double>& points, int ppn = CURVES_MIN_POLY_POINTS);
    virtual ~DiagonalCurve ();

    double getVal     (double t) const;
    void   getVal     (const std::vector<double>& t, std::vector<double>& res) const;
    bool   isIdentity () const
    {
        return kind == DCT_Empty;
    };
};

class FlatCurve : public Curve
{

private:
    FlatCurveType kind;
    double* leftTangent;
    double* rightTangent;
    double identityValue;
    bool periodic;

    void CtrlPoints_set ();

public:

    FlatCurve (const std::vector<double>& points, bool isPeriodic = true, int ppn = CURVES_MIN_POLY_POINTS);
    virtual ~FlatCurve ();

    double getVal     (double t) const;
    void   getVal     (const std::vector<double>& t, std::vector<double>& res) const;
    bool   setIdentityValue (double iVal);
    bool   isIdentity () const
    {
        return kind == FCT_Empty;
    };
};

class RetinextransmissionCurve
{
private:
    LUTf luttransmission;  // 0xffff range
    void Set(const Curve &pCurve);

public:
    virtual ~RetinextransmissionCurve() {};
    RetinextransmissionCurve();

    void Reset();
    void Set(const Curve *pCurve);
    void Set(const std::vector<double> &curvePoints);
    float operator[](float index) const
    {
        return luttransmission[index];
    }

    operator bool (void) const
    {
        return luttransmission;
    }
};

class RetinexgaintransmissionCurve
{
private:
    LUTf lutgaintransmission;  // 0xffff range
    void Set(const Curve &pCurve);

public:
    virtual ~RetinexgaintransmissionCurve() {};
    RetinexgaintransmissionCurve();

    void Reset();
    void Set(const Curve *pCurve);
    void Set(const std::vector<double> &curvePoints);
    float operator[](float index) const
    {
        return lutgaintransmission[index];
    }

    operator bool (void) const
    {
        return lutgaintransmission;
    }
};



class ToneCurve
{
public:
    LUTf lutToneCurve;  // 0xffff range

    virtual ~ToneCurve() {};

    void Reset();
    void Set(const Curve &pCurve, float gamma = 0);
    operator bool (void) const
    {
        return lutToneCurve;
    }
};

class OpacityCurve
{
public:
    LUTf lutOpacityCurve;  // 0xffff range

    virtual ~OpacityCurve() {};

    void Reset();
    void Set(const Curve *pCurve);
    void Set(const std::vector<double> &curvePoints, bool &opautili);

    // TODO: transfer this method to the Color class...
    float blend (float x, float lower, float upper) const
    {
        return (upper - lower) * lutOpacityCurve[x * 500.f] + lower;
    }
    void blend3f (float x, float lower1, float upper1, float &result1, float lower2, float upper2, float &result2, float lower3, float upper3, float &result3) const
    {
        float opacity = lutOpacityCurve[x * 500.f];
        result1 = (upper1 - lower1) * opacity + lower1;
        result2 = (upper2 - lower2) * opacity + lower2;
        result3 = (upper3 - lower3) * opacity + lower3;
    }

    operator bool (void) const
    {
        return lutOpacityCurve;
    }
};

class WavCurve
{
private:
    LUTf lutWavCurve;  // 0xffff range
    void Set(const Curve &pCurve);

public:
    float sum;

    virtual ~WavCurve() {};
    WavCurve();
    void Reset();
    void Set(const std::vector<double> &curvePoints);
    float getSum() const
    {
        return sum;
    }

    float operator[](float index) const
    {
        return lutWavCurve[index];
    }
    operator bool (void) const
    {
        return lutWavCurve;
    }
};

class WavOpacityCurveRG
{
private:
    LUTf lutOpacityCurveRG;  // 0xffff range
    void Set(const Curve &pCurve);
public:
    virtual ~WavOpacityCurveRG() {};
    WavOpacityCurveRG();

    void Reset();
    //  void Set(const std::vector<double> &curvePoints, bool &opautili);
    void Set(const std::vector<double> &curvePoints);
    float operator[](float index) const
    {
        return lutOpacityCurveRG[index];
    }

    operator bool (void) const
    {
        return lutOpacityCurveRG;
    }
};
class WavOpacityCurveBY
{
private:
    LUTf lutOpacityCurveBY;  // 0xffff range
    void Set(const Curve &pCurve);

public:
    virtual ~WavOpacityCurveBY() {};
    WavOpacityCurveBY();

    void Reset();
    void Set(const Curve *pCurve);
    void Set(const std::vector<double> &curvePoints);
    float operator[](float index) const
    {
        return lutOpacityCurveBY[index];
    }

    operator bool (void) const
    {
        return lutOpacityCurveBY;
    }
};
class WavOpacityCurveW
{
private:
    LUTf lutOpacityCurveW;  // 0xffff range
    void Set(const Curve &pCurve);

public:
    virtual ~WavOpacityCurveW() {};
    WavOpacityCurveW();

    void Reset();
    void Set(const Curve *pCurve);
    void Set(const std::vector<double> &curvePoints);
    float operator[](float index) const
    {
        return lutOpacityCurveW[index];
    }

    operator bool (void) const
    {
        return lutOpacityCurveW;
    }
};

class WavOpacityCurveWL
{
private:
    LUTf lutOpacityCurveWL;  // 0xffff range
    void Set(const Curve &pCurve);

public:
    virtual ~WavOpacityCurveWL() {};
    WavOpacityCurveWL();

    void Reset();
    void Set(const Curve *pCurve);
    void Set(const std::vector<double> &curvePoints);
    float operator[](float index) const
    {
        return lutOpacityCurveWL[index];
    }

    operator bool (void) const
    {
        return lutOpacityCurveWL;
    }
};

class NoiseCurve
{
private:
    LUTf lutNoiseCurve;  // 0xffff range
    float sum;
    void Set(const Curve &pCurve);

public:
    virtual ~NoiseCurve() {};
    NoiseCurve();
    void Reset();
    void Set(const std::vector<double> &curvePoints);

    float getSum() const
    {
        return sum;
    }
    float operator[](float index) const
    {
        return lutNoiseCurve[index];
    }
    operator bool (void) const
    {
        return lutNoiseCurve;
    }
};

class ColorGradientCurve
{
public:
    LUTf   lut1;    // [0.;1.] range (float values)
    LUTf   lut2;  // [0.;1.] range (float values)
    LUTf   lut3;   // [0.;1.] range (float values)
    double low;
    double high;

    virtual ~ColorGradientCurve() {};

    void Reset();
    void SetXYZ(const Curve *pCurve, const double xyz_rgb[3][3], float satur, float lumin);
    void SetXYZ(const std::vector<double> &curvePoints, const double xyz_rgb[3][3], float satur, float lumin);
    void SetRGB(const Curve *pCurve);
    void SetRGB(const std::vector<double> &curvePoints);

    /**
    * @brief Get the value of Red, Green and Blue corresponding to the requested index
    * @param index value in the [0 ; 1] range
    * @param r corresponding red value [0 ; 65535] (return value)
    * @param g corresponding green value [0 ; 65535] (return value)
    * @param b corresponding blue value [0 ; 65535] (return value)
    */
    void getVal(float index, float &r, float &g, float &b) const;
    operator bool (void) const
    {
        return lut1 && lut2 && lut3;
    }
};

class ColorAppearance
{
public:
    LUTf lutColCurve;  // 0xffff range

    virtual ~ColorAppearance() {};

    void Reset();
    void Set(const Curve &pCurve);
    operator bool (void) const
    {
        return lutColCurve;
    }
};

class Lightcurve : public ColorAppearance
{
public:
    void Apply(float& Li) const;
};

//lightness curve
inline void Lightcurve::Apply (float& Li) const
{

    assert (lutColCurve);

    curves::setLutVal(lutColCurve, Li);
}

class Brightcurve : public ColorAppearance
{
public:
    void Apply(float& Br) const;
};

//brightness curve
inline void Brightcurve::Apply (float& Br) const
{

    assert (lutColCurve);

    curves::setLutVal(lutColCurve, Br);
}

class Chromacurve : public ColorAppearance
{
public:
    void Apply(float& Cr) const;
};

//Chroma curve
inline void Chromacurve::Apply (float& Cr) const
{

    assert (lutColCurve);

    curves::setLutVal(lutColCurve, Cr);
}
class Saturcurve : public ColorAppearance
{
public:
    void Apply(float& Sa) const;
};

//Saturation curve
inline void Saturcurve::Apply (float& Sa) const
{

    assert (lutColCurve);

    curves::setLutVal(lutColCurve, Sa);
}

class Colorfcurve : public ColorAppearance
{
public:
    void Apply(float& Cf) const;
};

//Colorfullness curve
inline void Colorfcurve::Apply (float& Cf) const
{

    assert (lutColCurve);

    curves::setLutVal(lutColCurve, Cf);
}


class StandardToneCurve : public ToneCurve
{
public:
    void Apply(float& r, float& g, float& b) const;

    // Applies the tone curve to `r`, `g`, `b` arrays, starting at `r[start]`
    // and ending at `r[end]` (and respectively for `b` and `g`). Uses SSE
    // and requires that `r`, `g`, and `b` pointers have the same alignment.
    void BatchApply(
            const size_t start, const size_t end,
            float *r, float *g, float *b) const;
};

class AdobeToneCurve : public ToneCurve
{
private:
    void RGBTone(float& r, float& g, float& b) const;  // helper for tone curve

public:
    void Apply(float& r, float& g, float& b) const;
};

class SatAndValueBlendingToneCurve : public ToneCurve
{
public:
    void Apply(float& r, float& g, float& b) const;
};

class WeightedStdToneCurve : public ToneCurve
{
private:
    float Triangle(float refX, float refY, float X2) const;
#ifdef __SSE2__
    vfloat Triangle(vfloat refX, vfloat refY, vfloat X2) const;
#endif
public:
    void Apply(float& r, float& g, float& b) const;
    void BatchApply(const size_t start, const size_t end, float *r, float *g, float *b) const;
};

class LuminanceToneCurve : public ToneCurve
{
public:
    void Apply(float& r, float& g, float& b) const;
};

class PerceptualToneCurveState
{
public:
    float Working2Prophoto[3][3];
    float Prophoto2Working[3][3];
    float cmul_contrast;
    bool isProphoto;
};

// Tone curve whose purpose is to keep the color appearance constant, that is the curve changes contrast
// but colors appears to have the same hue and saturation as before. As contrast and saturation is tightly
// coupled in human vision saturation is modulated based on the curve's contrast, and that way the appearance
// can be kept perceptually constant (within limits).
class PerceptualToneCurve : public ToneCurve
{
private:
    static float cf_range[2];
    static float cf[1000];
    // for ciecam02
    static float f, c, nc, yb, la, xw, yw, zw, gamut;
    static float n, d, nbb, ncb, cz, aw, wh, pfl, fl, pow1;

    static void cubic_spline(const float x[], const float y[], const int len, const float out_x[], float out_y[], const int out_len);
    static float find_minimum_interval_halving(float (*func)(float x, void *arg), void *arg, float a, float b, float tol, int nmax);
    static float find_tc_slope_fun(float k, void *arg);
    static float get_curve_val(float x, float range[2], float lut[], size_t lut_size);
    float calculateToneCurveContrastValue() const;
public:
    static void init();
    void initApplyState(PerceptualToneCurveState & state, Glib::ustring workingSpace) const;
    void BatchApply(const size_t start, const size_t end, float *r, float *g, float *b, const PerceptualToneCurveState &state) const;
};

// Standard tone curve
inline void StandardToneCurve::Apply (float& r, float& g, float& b) const
{

    assert (lutToneCurve);

    curves::setLutVal(lutToneCurve, r);
    curves::setLutVal(lutToneCurve, g);
    curves::setLutVal(lutToneCurve, b);
}

inline void StandardToneCurve::BatchApply(
        const size_t start, const size_t end,
        float *r, float *g, float *b) const {
    assert (lutToneCurve);
    assert (lutToneCurve.getClip() & LUT_CLIP_BELOW);
    assert (lutToneCurve.getClip() & LUT_CLIP_ABOVE);

    // All pointers must have the same alignment for SSE usage. In the loop body below,
    // we will only check `r`, assuming that the same result would hold for `g` and `b`.
    assert (reinterpret_cast<uintptr_t>(r) % 16 == reinterpret_cast<uintptr_t>(g) % 16);
    assert (reinterpret_cast<uintptr_t>(g) % 16 == reinterpret_cast<uintptr_t>(b) % 16);

    size_t i = start;
    while (true) {
        if (i >= end) {
            // If we get to the end before getting to an aligned address, just return.
            // (Or, for non-SSE mode, if we get to the end.)
            return;
#ifdef __SSE2__
        } else if (reinterpret_cast<uintptr_t>(&r[i]) % 16 == 0) {
            // Otherwise, we get to the first aligned address; go to the SSE part.
            break;
#endif
        }
        curves::setLutVal(lutToneCurve, r[i]);
        curves::setLutVal(lutToneCurve, g[i]);
        curves::setLutVal(lutToneCurve, b[i]);
        i++;
    }

#ifdef __SSE2__
    float tmpr[4];
    float tmpg[4];
    float tmpb[4];
    float mv = lutToneCurve[MAXVALF];
    for (; i + 3 < end; i += 4) {
        __m128 r_val = LVF(r[i]);
        __m128 g_val = LVF(g[i]);
        __m128 b_val = LVF(b[i]);
        STVF(tmpr[0], lutToneCurve[r_val]);
        STVF(tmpg[0], lutToneCurve[g_val]);
        STVF(tmpb[0], lutToneCurve[b_val]);
        for (int j = 0; j < 4; ++j) {
            curves::setLutVal(r[i+j], tmpr[j], mv);
            curves::setLutVal(g[i+j], tmpg[j], mv);
            curves::setLutVal(b[i+j], tmpb[j], mv);
        }
    }

    // Remainder in non-SSE.
    for (; i < end; ++i) {
        curves::setLutVal(lutToneCurve, r[i]);
        curves::setLutVal(lutToneCurve, g[i]);
        curves::setLutVal(lutToneCurve, b[i]);
    }
#endif
}

// Tone curve according to Adobe's reference implementation
// values in 0xffff space
// inlined to make sure there will be no cache flush when used
inline void AdobeToneCurve::Apply (float& ir, float& ig, float& ib) const
{

    assert (lutToneCurve);
    float r = CLIP(ir);
    float g = CLIP(ig);
    float b = CLIP(ib);

    if (r >= g) {
        if      (g > b) {
            RGBTone (r, g, b);    // Case 1: r >= g >  b
        } else if (b > r) {
            RGBTone (b, r, g);    // Case 2: b >  r >= g
        } else if (b > g) {
            RGBTone (r, b, g);    // Case 3: r >= b >  g
        } else {                           // Case 4: r >= g == b
            r = lutToneCurve[r];
            g = lutToneCurve[g];
            b = g;
        }
    } else {
        if      (r >= b) {
            RGBTone (g, r, b);    // Case 5: g >  r >= b
        } else if (b >  g) {
            RGBTone (b, g, r);    // Case 6: b >  g >  r
        } else {
            RGBTone (g, b, r);    // Case 7: g >= b >  r
        }
    }

    setUnlessOOG(ir, r);
    setUnlessOOG(ig, g);
    setUnlessOOG(ib, b);
}

inline void AdobeToneCurve::RGBTone (float& r, float& g, float& b) const
{
    float rold = r, gold = g, bold = b;

    r = lutToneCurve[rold];
    b = lutToneCurve[bold];
    g = b + ((r - b) * (gold - bold) / (rold - bold));
}

// Modifying the Luminance channel only
inline void LuminanceToneCurve::Apply(float &ir, float &ig, float &ib) const
{
    assert (lutToneCurve);

    float r = CLIP(ir);
    float g = CLIP(ig);
    float b = CLIP(ib);

    float currLuminance = r * 0.2126729f + g * 0.7151521f + b * 0.0721750f;
    const float newLuminance = lutToneCurve[currLuminance];
    currLuminance = currLuminance == 0.f ? 0.00001f : currLuminance;
    const float coef = newLuminance / currLuminance;
    r = LIM<float>(r * coef, 0.f, 65535.f);
    g = LIM<float>(g * coef, 0.f, 65535.f);
    b = LIM<float>(b * coef, 0.f, 65535.f);

    setUnlessOOG(ir, r);
    setUnlessOOG(ig, g);
    setUnlessOOG(ib, b);
}

inline float WeightedStdToneCurve::Triangle(float a, float a1, float b) const
{
    if (a != b) {
        float b1;
        float a2 = a1 - a;

        if (b < a) {
            b1 = b + a2 *      b  /     a ;
        } else       {
            b1 = b + a2 * (65535.f - b) / (65535.f - a);
        }

        return b1;
    }

    return a1;
}

#ifdef __SSE2__
inline vfloat WeightedStdToneCurve::Triangle(vfloat a, vfloat a1, vfloat b) const
{
        vfloat a2 = a1 - a;
        vmask cmask = vmaskf_lt(b, a);
        vfloat b3 = vself(cmask, b, F2V(65535.f) - b);
        vfloat a3 = vself(cmask, a, F2V(65535.f) - a);
        return b + a2 * b3 / a3;
}
#endif

// Tone curve modifying the value channel only, preserving hue and saturation
// values in 0xffff space
inline void WeightedStdToneCurve::Apply (float& ir, float& ig, float& ib) const
{

    assert (lutToneCurve);

    float r = CLIP(ir);
    float g = CLIP(ig);
    float b = CLIP(ib);
    float r1 = lutToneCurve[r];
    float g1 = Triangle(r, r1, g);
    float b1 = Triangle(r, r1, b);

    float g2 = lutToneCurve[g];
    float r2 = Triangle(g, g2, r);
    float b2 = Triangle(g, g2, b);

    float b3 = lutToneCurve[b];
    float r3 = Triangle(b, b3, r);
    float g3 = Triangle(b, b3, g);

    r = CLIP<float>(r1 * 0.50f + r2 * 0.25f + r3 * 0.25f);
    g = CLIP<float>(g1 * 0.25f + g2 * 0.50f + g3 * 0.25f);
    b = CLIP<float>(b1 * 0.25f + b2 * 0.25f + b3 * 0.50f);

    setUnlessOOG(ir, r);
    setUnlessOOG(ig, g);
    setUnlessOOG(ib, b);
}

inline void WeightedStdToneCurve::BatchApply(const size_t start, const size_t end, float *r, float *g, float *b) const {
    assert (lutToneCurve);
    assert (lutToneCurve.getClip() & LUT_CLIP_BELOW);
    assert (lutToneCurve.getClip() & LUT_CLIP_ABOVE);

    // All pointers must have the same alignment for SSE usage. In the loop body below,
    // we will only check `r`, assuming that the same result would hold for `g` and `b`.
    assert (reinterpret_cast<uintptr_t>(r) % 16 == reinterpret_cast<uintptr_t>(g) % 16);
    assert (reinterpret_cast<uintptr_t>(g) % 16 == reinterpret_cast<uintptr_t>(b) % 16);

    size_t i = start;
    while (true) {
        if (i >= end) {
            // If we get to the end before getting to an aligned address, just return.
            // (Or, for non-SSE mode, if we get to the end.)
            return;
#ifdef __SSE2__
        } else if (reinterpret_cast<uintptr_t>(&r[i]) % 16 == 0) {
            // Otherwise, we get to the first aligned address; go to the SSE part.
            break;
#endif
        }
        Apply(r[i], g[i], b[i]);
        i++;
    }

#ifdef __SSE2__
    const vfloat c65535v = F2V(65535.f);
    const vfloat zd5v = F2V(0.5f);
    const vfloat zd25v = F2V(0.25f);

    float tmpr[4];
    float tmpg[4];
    float tmpb[4];

    for (; i + 3 < end; i += 4) {
        vfloat r_val = LIMV(LVF(r[i]), ZEROV, c65535v);
        vfloat g_val = LIMV(LVF(g[i]), ZEROV, c65535v);
        vfloat b_val = LIMV(LVF(b[i]), ZEROV, c65535v);
        vfloat r1 = lutToneCurve[r_val];
        vfloat g1 = Triangle(r_val, r1, g_val);
        vfloat b1 = Triangle(r_val, r1, b_val);

        vfloat g2 = lutToneCurve[g_val];
        vfloat r2 = Triangle(g_val, g2, r_val);
        vfloat b2 = Triangle(g_val, g2, b_val);

        vfloat b3 = lutToneCurve[b_val];
        vfloat r3 = Triangle(b_val, b3, r_val);
        vfloat g3 = Triangle(b_val, b3, g_val);

        STVF(tmpr[0], LIMV(r1 * zd5v + r2 * zd25v + r3 * zd25v, ZEROV, c65535v));
        STVF(tmpg[0], LIMV(g1 * zd25v + g2 * zd5v + g3 * zd25v, ZEROV, c65535v));
        STVF(tmpb[0], LIMV(b1 * zd25v + b2 * zd25v + b3 * zd5v, ZEROV, c65535v));
        for (int j = 0; j < 4; ++j) {
            setUnlessOOG(r[i+j], tmpr[j]);
            setUnlessOOG(g[i+j], tmpg[j]);
            setUnlessOOG(b[i+j], tmpb[j]);
        }
    }

    // Remainder in non-SSE.
    for (; i < end; ++i) {
        Apply(r[i], g[i], b[i]);
    }
#endif
}

// Tone curve modifying the value channel only, preserving hue and saturation
// values in 0xffff space
inline void SatAndValueBlendingToneCurve::Apply (float& ir, float& ig, float& ib) const
{

    assert (lutToneCurve);

    float r = CLIP(ir);
    float g = CLIP(ig);
    float b = CLIP(ib);

    const float lum = (r + g + b) / 3.f;
    const float newLum = lutToneCurve[lum];

    if (newLum == lum) {
        return;
    }

    float h, s, v;
    Color::rgb2hsvtc(r, g, b, h, s, v);

    float dV;
    if (newLum > lum) {
        // Linearly targeting Value = 1 and Saturation = 0
        const float coef = (newLum - lum) / (65535.f - lum);
        dV = (1.f - v) * coef;
        s *= 1.f - coef;
    } else {
        // Linearly targeting Value = 0
        const float coef = (newLum - lum) / lum ;
        dV = v * coef;
    }
    Color::hsv2rgbdcp(h, s, v + dV, r, g, b);

    setUnlessOOG(ir, r);
    setUnlessOOG(ig, g);
    setUnlessOOG(ib, b);
}

}

#undef CLIPI

#endif
