/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
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
 *  2016 - 2020 Jacques Desmis <jdesmis@gmail.com>
 *  2016 - 2020 Ingo Weyrich <heckflosse@i-weyrich.de>

 */
#include <cmath>
#include <fftw3.h>

#include "improcfun.h"
#include "colortemp.h"
#include "curves.h"
#include "gauss.h"
#include "iccstore.h"
#include "imagefloat.h"
#include "labimage.h"
#include "color.h"
#include "rt_math.h"
#include "jaggedarray.h"
#include "rt_algo.h"
#include "settings.h"
#include "../rtgui/options.h"
#include "utils.h"
#include "iccmatrices.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "../rtgui/thresholdselector.h"
#include "imagesource.h"

#include "cplx_wavelet_dec.h"
#include "ciecam02.h"

#define BENCHMARK
#include "StopWatch.h"
#include "guidedfilter.h"
#include "boxblur.h"
#include "rescale.h"



#pragma GCC diagnostic warning "-Wall"
#pragma GCC diagnostic warning "-Wextra"
#pragma GCC diagnostic warning "-Wdouble-promotion"

namespace
{

constexpr int limscope = 80;
constexpr int mSPsharp = 39; //minimum size Spot Sharp due to buildblendmask
constexpr int mSPwav = 32; //minimum size Spot Wavelet
constexpr int mDEN = 128; //minimum size Spot Denoise
constexpr int mSP = 5; //minimum size Spot
constexpr float MAXSCOPE = 1.25f;
constexpr float MINSCOPE = 0.025f;
constexpr int TS = 64; // Tile size
constexpr float epsilonw = 0.001f / (TS * TS); //tolerance
constexpr int offset = 25; // shift between tiles
constexpr double czlim = rtengine::RT_SQRT1_2;// 0.70710678118654752440;

constexpr float clipLoc(float x)
{
    return rtengine::LIM(x, 0.f, 32767.f);
}

constexpr float clipDE(float x)
{
    return rtengine::LIM(x, 0.3f, 1.f);
}

constexpr float clipC(float x)
{
    return rtengine::LIM(x, -42000.f, 42000.f);
}

constexpr float clipChro(float x)
{
    return rtengine::LIM(x, 0.f, 140.f);
}

constexpr double clipazbz(double x)
{
    return rtengine::LIM(x, -0.5, 0.5);
}

constexpr double clipcz(double x)
{
    return rtengine::LIM(x, 0., czlim);
}


constexpr double clipjz05(double x)
{
    return rtengine::LIM(x, 0.0006, 1.0);
}

float softlig(float a, float b, float minc, float maxc)
{
    // as Photoshop
    if (2.f * b <= maxc - minc) {
        return a * (2.f * b + a * (maxc - 2.f * b));
    } else {
        return 2.f * a * (maxc - b) + std::sqrt(rtengine::LIM(a, 0.f, 2.f)) * (2.f * b - maxc);
    }
}

float softlig3(float a, float b)
{
    // as w3C
    if (2.f * b <= 1.f) {
        return a - (1.f - 2.f * b) * a * (1.f - a);
    } else {
        if (4.f * a <= 1.f) {
            return a + (2.f * b - 1.f) * (4.f * a * (4.f * a + 1.f) * (a - 1.f) + 7.f * a);
        } else {
            return a + (2.f * a - 1.f) * (std::sqrt(a) - a);
        }
    }
}

float softlig2(float a, float b)
{
    // illusions.hu
    return pow_F(b, pow_F(2.f, (2.f * (0.5f - a))));
}

constexpr float colburn(float a, float b)
{
    // w3C
    return b == 0.f ? 0.f : 1.f - rtengine::min(1.f, (1.f - a) / b);
}

constexpr float coldodge(float a, float b)
{
    // w3C
    return b == 1.f ? 1.f : rtengine::min(1.f, a / (1.f - b));
}

float overlay(float a, float b, float minc, float maxc)
{
    if (2.f * b <= maxc - minc) {
        return 2.f * b * a;
    } else {
        return maxc - 2.f * (1.f - a) * (maxc - b);
    }
}

constexpr float screen(float a, float b, float maxc)
{
    return 1.f - (1.f - a) * (maxc - b);
}

constexpr float exclusion(float a, float b)
{
    return a + b - 2.f * a * b;
}

void calcdif(float lmr, float &lmrc)
{   //approximative change between gamma sRGB g=2.4 s=12.92 and gamma LAB g=3.0 s=9.03
    //useful to calculate action with dark and light area mask
    //differences in 3 parts linear...very small differences with real...
    float a0 = 7.6f / 11.6f;//11.6 sRGB  - 7.6 Lab...11.6 max difference
    float a01 = 62.f - 7.6f; //60 sRGB 62 Lab   60 max difference
    float a11 = 60.f - 11.6f;
    float a1 = a01 / a11;
    float b1 = 62.f - a1 * 60.f;
    float a2 = (100.f - 62.f) / (100.f - 60.f);
    float b2 = 100.f - a2 * 100.f;
    if(lmr < 11.6f) {
        lmrc = a0 * lmr; 
    } else if (lmr < 60.f) {
        lmrc = a1 * lmr + b1;
    } else {
        lmrc = a2 * lmr + b2;
    }
}

void calcGammaLut(double gamma, double ts, LUTf &gammaLut)
{
    double pwr = 1.0 / gamma;
    double gamm = gamma;
    const double gamm2 = gamma;
    rtengine::GammaValues g_a;

    if (gamm2 < 1.0) {
        std::swap(pwr, gamm);
    }

    rtengine::Color::calcGamma(pwr, ts, g_a); // call to calcGamma with selected gamma and slope

    const double start = gamm2 < 1. ? g_a[2] : g_a[3];
    const double add = g_a[4];
    const double mul = 1.0 + g_a[4];

    if (gamm2 < 1.) {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1024)
#endif
        for (int i = 0; i < 65536; i++) {
            const double x = rtengine::Color::igammareti(i / 65535.0, gamm, start, ts, mul, add);
            gammaLut[i] = 0.5 * rtengine::CLIP(x * 65535.0);  // CLIP avoid in some case extra values
        }
    } else {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1024)
#endif
        for (int i = 0; i < 65536; i++) {
            const double x = rtengine::Color::gammareti(i / 65535.0, gamm, start, ts, mul, add);
            gammaLut[i] = 0.5 * rtengine::CLIP(x * 65535.0);  // CLIP avoid in some case extra values
        }
    }
}

float calcLocalFactor(const float lox, const float loy, const float lcx, const float dx, const float lcy, const float dy, const float ach, const float gradient)
{
    //ellipse x2/a2 + y2/b2=1
    //transition ellipsoidal
    const float kelip = dx / dy;
    const float belip = rtengine::max(0.0001f, std::sqrt((rtengine::SQR((lox - lcx) / kelip) + rtengine::SQR(loy - lcy)))); //determine position ellipse ==> a and b

    //gradient allows differentiation between transition x and y
    const float rapy = std::fabs((loy - lcy) / belip);
    const float aelip = belip * kelip;
    const float degrad = aelip / dx;
    const float gradreal = gradient * rapy + 1.f;
    const float ap = rtengine::RT_PI_F / (1.f - ach);
    const float bp = rtengine::RT_PI_F - ap;
    return pow(0.5f * (1.f + xcosf(degrad * ap + bp)), rtengine::SQR(gradreal)); // trigo cos transition
}

float calcLocalFactorrect(const float lox, const float loy, const float lcx, const float dx, const float lcy, const float dy, const float ach, const float gradient)
{
    constexpr float eps = 0.0001f;
    const float krap = std::fabs(dx / dy);
    const float kx = lox - lcx;
    const float ky = loy - lcy;

    float ref;
    //gradient allows differentiation between transition x and y
    if (std::fabs(kx / (ky + eps)) < krap) {
        ref = std::sqrt(rtengine::SQR(dy) * (1.f + rtengine::SQR(kx / (ky + eps))));
    } else {
        ref = std::sqrt(rtengine::SQR(dx) * (1.f + rtengine::SQR(ky / (kx + eps))));
    }

    const float rad = rtengine::max(eps, std::sqrt(rtengine::SQR(kx) + rtengine::SQR(ky)));
    const float rapy = std::fabs((loy - lcy) / rad);
    const float gradreal = gradient * rapy + 1.f;

    const float coef = rad / ref;
    const float fact = (coef - 1.f) / (ach - 1.f);
    return pow(fact, rtengine::SQR(gradreal));
}

float calcreducdE(float dE, float maxdE, float mindE, float maxdElim, float mindElim, float iterat, int limscope, int scope)
{
    if (scope > limscope) {//80 arbitrary value, if we change we must change limscope
        if (dE > maxdElim) {
            return 0.f;
        } else if (dE > mindElim) {
            const float reducdElim = std::pow((dE - maxdElim) / (mindElim - maxdElim), iterat);
            const float aalim = (1.f - reducdElim) / 20.f;
            const float bblim = 1.f - 100.f * aalim;
            return aalim * scope + bblim;
        } else {
            return 1.f;
        }
    } else {
        if (dE > maxdE) {
            return 0.f;
        } else if (dE > mindE) {
            return std::pow((dE - maxdE) / (mindE - maxdE), iterat);
        } else {
            return 1.f;
        }
    }
}



void deltaEforLaplace(float *dE, const float lap, int bfw, int bfh, rtengine::LabImage* bufexporig, const float hueref, const float chromaref, const float lumaref)
{

    const float refa = chromaref * std::cos(hueref);
    const float refb = chromaref * std::sin(hueref);
    const float refL = lumaref;
    float maxdE = 5.f + MAXSCOPE * lap;

    float maxC = std::sqrt((rtengine::SQR(refa - bufexporig->a[0][0]) + rtengine::SQR(refb - bufexporig->b[0][0])) + rtengine::SQR(refL - bufexporig->L[0][0])) / 327.68f;
#ifdef _OPENMP
    #pragma omp parallel for reduction(max:maxC)
#endif
    for (int y = 0; y < bfh; y++) {
        for (int x = 0; x < bfw; x++) {
            const float val = std::sqrt((rtengine::SQR(refa - bufexporig->a[y][x]) + rtengine::SQR(refb - bufexporig->b[y][x])) + rtengine::SQR(refL - bufexporig->L[y][x])) / 327.68f;
            dE[y * bfw + x] = val;
            maxC = rtengine::max(maxC, val);
        }
    }

    if (maxdE > maxC) {
        maxdE = maxC - 1.f;
    }

    const float ade = 1.f / (maxdE - maxC);
  //  const float bde = -ade * maxC;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif
    for (int y = 0; y < bfh; y++) {
        for (int x = 0; x < bfw; x++) {
          //  dE[y * bfw + x] = dE[y * bfw + x] >= maxdE ? ade * dE[y * bfw + x] + bde : 1.f;
            dE[y * bfw + x] = dE[y * bfw + x] >= maxdE ? ade * (dE[y * bfw + x] - maxC) : 1.f;

        }
    }
}

float calclight(float lum, const LUTf &lightCurveloc)
{
    return clipLoc(lightCurveloc[lum]);
}

float calclightinv(float lum, float koef, const LUTf &lightCurveloc)
{
    return koef != -100.f ? clipLoc(lightCurveloc[lum]) : 0.f;
}

float balancedeltaE(float kL)
{
    constexpr float mincurs = 0.3f; // minimum slider balan_
    constexpr float maxcurs = 1.7f; // maximum slider balan_
    constexpr float maxkab = 1.35; // 0.5 * (3 - 0.3)
    constexpr float minkab = 0.65; // 0.5 * (3 - 1.7)
    constexpr float abal = (maxkab - minkab) / (mincurs - maxcurs);
    constexpr float bbal = maxkab - mincurs * abal;
    return abal * kL + bbal;
}

void SobelCannyLuma(float **sobelL, float **luma, int bfw, int bfh, float radius)
{
    // base of the process to detect shape in complement of deltaE
    // use for calculate Spot reference
    // and for structure of the shape
    // actually , as the program don't use these function, I just create a simple "Canny" near of Sobel. This can be completed after with teta, etc.
    array2D<float> tmL(bfw, bfh);

    //inspired from Chen Guanghua Zhang Xiaolong
    //Sobel Horizontal
    constexpr float GX[3][3] = {
        {1.f, 0.f, -1.f},
        {2.f, 0.f, -2.f},
        {1.f, 0.f, -1.f}
    };

    //Sobel Vertical
    constexpr float GY[3][3] = {
        {1.f, 2.f, 1.f},
        {0.f, 0.f, 0.f},
        {-1.f, -2.f, -1.f}
    };

    if (radius > 0.f) {
        gaussianBlur(luma, tmL, bfw, bfh, rtengine::max(radius / 2.f, 0.5f));
    } else {
        for (int y = 0; y < bfh ; y++) {
            for (int x = 0; x < bfw ; x++) {
                tmL[y][x] = luma[y][x];
            }
        }
    }

    for (int x = 0; x < bfw; x++) {
        sobelL[0][x] = 0.f;
    }
    for (int y = 1; y < bfh - 1; y++) {
        sobelL[y][0] = 0.f;
        for (int x = 1; x < bfw - 1; x++) {
            float sumXL = 0.f;
            float sumYL = 0.f;
            for (int i = -1; i < 2; i += 2) {
                for (int j = -1; j < 2; j += 1) {
                    sumXL += GX[j + 1][i + 1] * tmL[y + i][x + j];
                    sumYL += GY[j + 1][i + 1] * tmL[y + i][x + j];
                }
            }
            //Edge strength
            //we can add if need teta = atan2 (sumYr, sumXr)
            sobelL[y][x] = rtengine::min(std::sqrt(rtengine::SQR(sumXL) + rtengine::SQR(sumYL)), 32767.f);
        }
        sobelL[y][bfw - 1] = 0.f;
    }
    for (int x = 0; x < bfw; x++) {
        sobelL[bfh - 1][x] = 0.f;
    }
}


float igammalog(float x, float p, float s, float g2, float g4)
{
    return x <= g2 ? x / s : pow_F((x + g4) / (1.f + g4), p);//continuous
}

#ifdef __SSE2__
vfloat igammalog(vfloat x, vfloat p, vfloat s, vfloat g2, vfloat g4)
{
  //  return x <= g2 ? x / s : pow_F((x + g4) / (1.f + g4), p);//continuous
    return vself(vmaskf_le(x, g2), x / s, pow_F((x + g4) / (F2V(1.f) + g4), p));
}
#endif

float gammalog(float x, float p, float s, float g3, float g4)
{
    return x <= g3 ? x * s : (1.f + g4) * xexpf(xlogf(x) / p) - g4;//used by Nlmeans
}

#ifdef __SSE2__
vfloat gammalog(vfloat x, vfloat p, vfloat s, vfloat g3, vfloat g4)
{
  //  return x <= g3 ? x * s : (1.f + g4) * xexpf(xlogf(x) / p) - g4;//continuous
    return vself(vmaskf_le(x, g3), x * s, (F2V(1.f) + g4) * xexpf(xlogf(x) / p) - g4);//improve by Ingo - used by Nlmeans
    
}
#endif
}

namespace rtengine

{
extern MyMutex *fftwMutex;

using namespace procparams;

struct local_params {
    float yc, xc;
    float lx, ly;
    float lxL, lyT;
    float transweak;
    float transgrad;
    float iterat;
    float balance;
    float balanceh;
    int colorde;
    int cir;
    bool recur;
    float thr;
    float stru;
    int chro, cont, sens, sensh, senscb, sensbn, senstm, sensex, sensexclu, sensden, senslc, senssf, senshs, senscolor;
    float reparden;
    float repartm;
    float clarityml;
    float contresid;
    bool deltaem;
    float struco;
    float strengrid;
    float struexc;
    float blendmacol;
    float radmacol;
    float chromacol;
    float gammacol;
    float slomacol;
    float blendmalc;
    float radmalc;
    float chromalc;
    float radmaexp;
    float chromaexp;
    float gammaexp;
    float slomaexp;
    float strmaexp;
    float angmaexp;
    float str_mas;
    float ang_mas;
    float strexp;
    float angexp;
    float strSH;
    float angSH;
    float strcol;
    float strcolab;
    float strcolh;
    float angcol;
    float strvib;
    float strvibab;
    float strvibh;
    float angvib;
    float angwav;
    float strwav;
    float blendmaL;
    float radmaL;
    float chromaL;

    float strengthw;
    float radiusw;
    float detailw;
    float gradw;
    float tloww;
    float thigw;
    float edgw;
    float basew;

    float anglog;
    float strlog;
    float softradiusexp;
    float softradiuscol;
    float softradiuscb;
    float softradiusret;
    float softradiustm;
    float blendmaexp;
    float radmaSH;
    float blendmaSH;
    float chromaSH;
    float gammaSH;
    float slomaSH;
    float radmavib;
    float blendmavib;
    float chromavib;
    float gammavib;
    float slomavib;
    float radmacb;
    float blendmacb;
    float chromacbm;
    float gammacb;
    float slomacb;
    float radmatm;
    float blendmatm;
    float chromatm;
    float gammatm;
    float slomatm;

    float radmabl;
    float blendmabl;
    float chromabl;
    float gammabl;
    float slomabl;

    float struexp;
    float blurexp;
    float blurcol;
    float blurcolmask;
    float contcolmask;
    float blurSH;
    float ligh;
    float gamc;
    float gamlc;
    float gamex;
    float lowA, lowB, highA, highB;
    float lowBmerg, highBmerg, lowAmerg, highAmerg;
    int shamo, shdamp, shiter, senssha, sensv;
    float neig;
    float strng;
    float lap;
    float lcamount;
    double shrad;
    double shblurr;
    double rad;
    double stren;
    int it;
    int guidb;
    float strbl;
    float epsb;
    float trans;
    float feath;
    int dehaze;
    int dehazeSaturation;
    int depth;
    bool inv;
    bool invex;
    bool invsh;
    bool curvact;
    bool invrad;
    bool invret;
    bool equret;
    bool equtm;
    bool invshar;
    bool actsp;
    bool ftwlc;
    bool ftwreti;
    float str;
    int qualmet;
    int qualcurvemet;
    int gridmet;
    bool prevdE;
    int showmaskcolmet;
    int showmaskcolmetinv;
    int showmaskexpmet;
    int showmaskexpmetinv;
    int showmaskSHmet;
    int showmaskSHmetinv;
    int showmaskvibmet;
    int showmasklcmet;
    int showmasksharmet;
    int showmaskcbmet;
    int showmaskretimet;
    int showmasksoftmet;
    int showmasktmmet;
    int showmaskblmet;
    int showmasklogmet;
    int showmask_met;
    int showmaskciemet;
    bool fftbl;
    float laplacexp;
    float balanexp;
    float linear;
    int fullim;
    int expmet;
    int softmet;
    int blurmet;
    int blmet;
    bool invmaskd;
    bool invmask;
    int smasktyp;
    int chromet;
    int quamet;
    int shmeth;
    int medmet;
    int locmet;
    float noiself;
    float noiself0;
    float noiself2;
    float noiseldetail;
    int detailthr;
    float recothr;
    float lowthr;
    float higthr;
    float recothrd;
    float lowthrd;
    float midthrd;
    float midthrdch;
    float higthrd;
    float decayd;
    float recothrc;
    float lowthrc;
    float higthrc;
    float decayc;
    float recothre;
    float lowthre;
    float higthre;
    float decaye;
    float recothrv;
    float lowthrv;
    float higthrv;
    float decayv;
    float recothrcb;
    float lowthrcb;
    float higthrcb;
    float decaycb;
    float recothrt;
    float lowthrt;
    float higthrt;
    float decayt;
    float recothrw;
    float lowthrw;
    float higthrw;
    float decayw;
    float recothrr;
    float lowthrr;
    float higthrr;
    float decayr;
    float recothrs;
    float lowthrs;
    float higthrs;
    float decays;
    float recothrl;
    float lowthrl;
    float higthrl;
    float decayl;
    float recothrcie;
    float lowthrcie;
    float higthrcie;
    float decaycie;
    
    int noiselequal;
    float noisechrodetail;
    float bilat;
    int nlstr;
    int nldet;
    int nlpat;
    int nlrad;
    float nlgam;
    float noisegam;
    float noiselc;
    float noiselc4;
    float noiselc5;
    float noiselc6;
    float noisecf;
    float noisecc;
    float mulloc[6];
    int mullocsh[5];
    int detailsh;
    double tePivot;
    float threshol;
    float chromacb;
    float strengt;
    float gamm;
    float esto;
    float scalt;
    float rewe;
    float amo;
    bool colorena;
    bool blurena;
    bool tonemapena;
    bool retiena;
    bool sharpena;
    bool lcena;
    bool sfena;
    bool cbdlena;
    bool denoiena;
    bool wavcurvedenoi;
    bool expvib;
    bool exposena;
    bool hsena;
    bool vibena;
    bool logena;
    bool islocal;
    bool maskena;
    bool cieena;
    bool cut_past;
    float past;
    float satur;
    int blac;
    int shcomp;
    int shadex;
    int hlcomp;
    int hlcompthr;
    float expcomp;
    float expchroma;
    int excmet;
    int mergemet;
    int mergecolMethod;
    float opacol;
    int war;
    float adjch;
    int shapmet;
    int edgwmet;
    int neiwmet;
    bool enaColorMask;
    bool fftColorMask;
    bool enaColorMaskinv;
    bool enaExpMask;
    bool enaExpMaskinv;
    bool enaSHMask;
    bool enaSHMaskinv;
    bool enavibMask;
    bool enalcMask;
    bool enasharMask;
    bool enacbMask;
    bool enaretiMask;
    bool enaretiMasktmap;
    bool enatmMask;
    bool enablMask;
    bool enaLMask;
    bool ena_Mask;
    bool enacieMask;
    int highlihs;
    int shadowhs;
    int radiushs;
    int hltonalhs;
    int shtonalhs;
    int scalereti;
    float sourcegray;
    float targetgray;
    float blackev;
    float whiteev;
    float detail;
    int sensilog;
    int sensicie;
    int sensimas;
    bool Autogray;
    bool autocompute;
    float baselog;
    bool wavgradl;
    bool edgwena;
    bool lip3;
    int daubLen;
    float sigmadr;
    float sigmabl;
    float sigmaed;
    float sigmalc;
    float sigmalc2;
    float residsha;
    float residshathr;
    float residhi;
    float residhithr;
    float residgam;
    float residslop;
    bool blwh;
    bool fftma;
    float blurma;
    float contma;
    bool activspot;
    float thrlow;
    float thrhigh;
    bool usemask;
    float lnoiselow;
    float radmacie;
    float blendmacie;
    float chromacie;
    float denoichmask;
    float mLjz;
    float mCjz;
    float softrjz;

};

static void calcLocalParams(int sp, int oW, int oH, const LocallabParams& locallab, struct local_params& lp, bool prevDeltaE, int llColorMask, int llColorMaskinv, int llExpMask, int llExpMaskinv, int llSHMask, int llSHMaskinv, int llvibMask, int lllcMask, int llsharMask, int llcbMask, int llretiMask, int llsoftMask, int lltmMask, int llblMask, int lllogMask, int ll_Mask, int llcieMask, const LocwavCurve & locwavCurveden, bool locwavdenutili)
{
    int w = oW;
    int h = oH;
    int circr = locallab.spots.at(sp).circrad;
    bool recur = locallab.spots.at(sp).recurs;
    float streng = ((float)locallab.spots.at(sp).stren);
    float gam = ((float)locallab.spots.at(sp).gamma);
    float est = ((float)locallab.spots.at(sp).estop);
    float scal_tm = ((float)locallab.spots.at(sp).scaltm);
    float rewe = ((float)locallab.spots.at(sp).rewei);
    float amo = ((float)locallab.spots.at(sp).amount);
    float strlight = ((float)locallab.spots.at(sp).streng);
    float strucc = locallab.spots.at(sp).struc;
    float laplac = ((float)locallab.spots.at(sp).laplace);
    float thre = locallab.spots.at(sp).thresh;

//    if (thre > 8.f || thre < 0.f) {//to avoid artifacts if user does not clear cache with new settings. Can be suppressed after
//        thre = 2.f;
//    }
    thre = LIM(thre, 0.f, 10.0f);

    double local_x = locallab.spots.at(sp).loc.at(0) / 2000.0;
    double local_y = locallab.spots.at(sp).loc.at(2) / 2000.0;
    double local_xL = locallab.spots.at(sp).loc.at(1) / 2000.0;
    double local_yT = locallab.spots.at(sp).loc.at(3) / 2000.0;
    double local_center_x = locallab.spots.at(sp).centerX / 2000.0 + 0.5;
    double local_center_y = locallab.spots.at(sp).centerY / 2000.0 + 0.5;
    float iterati = (float) locallab.spots.at(sp).iter;
    float balanc = (float) locallab.spots.at(sp).balan;
    float balanch = (float) locallab.spots.at(sp).balanh;
    int colorde = (int) locallab.spots.at(sp).colorde;

//    if (iterati > 4.f || iterati < 0.2f) {//to avoid artifacts if user does not clear cache with new settings Can be suppressed after
//       iterati = 2.f;
//    }
    iterati = LIM(iterati, 0.2f, 10.0f);

    float neigh = float (locallab.spots.at(sp).neigh);
    float chromaPastel = float (locallab.spots.at(sp).pastels)   / 100.0f;
    float chromaSatur  = float (locallab.spots.at(sp).saturated) / 100.0f;
    int local_sensiv = locallab.spots.at(sp).sensiv;
    int local_sensiex = locallab.spots.at(sp).sensiex;

    if (locallab.spots.at(sp).qualityMethod == "enh") {
        lp.qualmet = 1;
    } else if (locallab.spots.at(sp).qualityMethod == "enhden") {
        lp.qualmet = 2;
    }

    if (locallab.spots.at(sp).qualitycurveMethod == "none") {
        lp.qualcurvemet = 0;
    } else if (locallab.spots.at(sp).qualitycurveMethod == "std") {
        lp.qualcurvemet = 1;
    }

    if (locallab.spots.at(sp).gridMethod == "one") {
        lp.gridmet = 0;
    } else if (locallab.spots.at(sp).gridMethod == "two") {
        lp.gridmet = 1;
    }
/*
    if (locallab.spots.at(sp).expMethod == "std") {
        lp.expmet = 0;
    } else if (locallab.spots.at(sp).expMethod == "pde") {
        lp.expmet = 1;
    }
*/
    lp.expmet = 1;

    if (locallab.spots.at(sp).localcontMethod == "loc") {
        lp.locmet = 0;
    } else if (locallab.spots.at(sp).localcontMethod == "wav") {
        lp.locmet = 1;
    }

    lp.laplacexp = locallab.spots.at(sp).laplacexp;
    lp.balanexp = locallab.spots.at(sp).balanexp;
    lp.linear = locallab.spots.at(sp).linear;
    if (locallab.spots.at(sp).spotMethod == "norm") {
        lp.fullim = 0;
    } else if(locallab.spots.at(sp).spotMethod == "exc"){
        lp.fullim = 1;
    } else if (locallab.spots.at(sp).spotMethod == "full"){
        lp.fullim = 2;
    }
   // printf("Lpfullim=%i\n", lp.fullim);
    
    lp.fftColorMask = locallab.spots.at(sp).fftColorMask;
    lp.prevdE = prevDeltaE;
    lp.showmaskcolmet = llColorMask;
    lp.showmaskcolmetinv = llColorMaskinv;
    lp.showmaskexpmet = llExpMask;
    lp.showmaskexpmetinv = llExpMaskinv;
    lp.showmaskSHmet = llSHMask;
    lp.showmaskSHmetinv = llSHMaskinv;
    lp.showmaskvibmet = llvibMask;
    lp.showmasklcmet = lllcMask;
    lp.showmasksharmet = llsharMask;
    lp.showmaskcbmet = llcbMask;
    lp.showmaskretimet = llretiMask;
    lp.showmasksoftmet = llsoftMask;
    
    lp.showmasktmmet = lltmMask;
    lp.showmaskblmet = llblMask;
    lp.showmasklogmet = lllogMask;
    lp.showmask_met = ll_Mask;
    lp.showmaskciemet = llcieMask;
//printf("CIEmask=%i\n", lp.showmaskciemet);
    lp.enaColorMask = locallab.spots.at(sp).enaColorMask && llsoftMask == 0 && llColorMaskinv == 0 && llSHMaskinv == 0 && llColorMask == 0 && llExpMaskinv == 0 && lllcMask == 0 && llsharMask == 0 && llExpMask == 0 && llSHMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;// Exposure mask is deactivated if Color & Light mask is visible
    lp.enaColorMaskinv = locallab.spots.at(sp).enaColorMask && llColorMaskinv == 0 && llSHMaskinv == 0 && llsoftMask == 0 && lllcMask == 0 && llsharMask == 0 && llExpMask == 0 && llSHMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;// Exposure mask is deactivated if Color & Light mask is visible
    lp.enaExpMask = locallab.spots.at(sp).enaExpMask && llExpMask == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llColorMask == 0 && llColorMaskinv == 0 && llsoftMask == 0 && lllcMask == 0 && llsharMask == 0 && llSHMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;// Exposure mask is deactivated if Color & Light mask is visible
    lp.enaExpMaskinv = locallab.spots.at(sp).enaExpMask && llExpMaskinv == 0 && llColorMask == 0 && llSHMaskinv == 0 && llColorMaskinv == 0 && llsoftMask == 0 && lllcMask == 0 && llsharMask == 0 && llSHMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;// Exposure mask is deactivated if Color & Light mask is visible
    lp.enaSHMask = locallab.spots.at(sp).enaSHMask && llSHMask == 0 && llColorMask == 0 && llColorMaskinv == 0 && llSHMaskinv == 0 && llExpMaskinv == 0 && llsoftMask == 0 && lllcMask == 0 && llsharMask == 0 && llExpMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && lllogMask == 0  && ll_Mask == 0 && llcieMask == 0;
    lp.enaSHMaskinv = locallab.spots.at(sp).enaSHMask && llColorMaskinv == 0 && llSHMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && lllcMask == 0 && llsoftMask == 0 && llsharMask == 0 && llColorMask == 0 && llExpMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;
    lp.enacbMask = locallab.spots.at(sp).enacbMask && llColorMaskinv == 0 && llcbMask == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && lllcMask == 0 && llsoftMask == 0 && llsharMask == 0 && llColorMask == 0 && llExpMask == 0 && llSHMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;
    lp.enaretiMask = locallab.spots.at(sp).enaretiMask && llColorMaskinv == 0 && lllcMask == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llsharMask == 0 && llsoftMask == 0 && llretiMask == 0 && llColorMask == 0 && llExpMask == 0 && llSHMask == 0 && llcbMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;
    lp.enatmMask = locallab.spots.at(sp).enatmMask && llColorMaskinv == 0 && lltmMask == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && lllcMask == 0 && llsoftMask == 0 && llsharMask == 0 && llColorMask == 0 && llExpMask == 0 && llSHMask == 0 && llcbMask == 0 && llretiMask == 0 && llblMask == 0 && llvibMask == 0&& lllogMask == 0  && ll_Mask == 0 && llcieMask == 0;
    lp.enablMask = locallab.spots.at(sp).enablMask && llColorMaskinv == 0 && llblMask == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && lllcMask == 0 && llsoftMask == 0 && llsharMask == 0 && llColorMask == 0 && llExpMask == 0 && llSHMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;
    lp.enavibMask = locallab.spots.at(sp).enavibMask && llvibMask == 0 && llColorMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && lllcMask == 0 && llsoftMask == 0 && llsharMask == 0 && llColorMask == 0 && llExpMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llSHMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;
    lp.enalcMask = locallab.spots.at(sp).enalcMask && lllcMask == 0 && llColorMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llcbMask == 0 && llsoftMask == 0 && llsharMask == 0 && llColorMask == 0 && llExpMask == 0 && llSHMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;
    lp.enasharMask = lllcMask == 0 && llcbMask == 0 && llsharMask == 0 && llColorMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llsoftMask == 0 && llColorMask == 0 && llExpMask == 0 && llSHMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;
    lp.ena_Mask = locallab.spots.at(sp).enamask && lllcMask == 0 && llColorMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llcbMask == 0 && llsoftMask == 0 && llsharMask == 0 && llColorMask == 0 && llExpMask == 0 && llSHMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && lllogMask == 0 && llvibMask == 0 && ll_Mask == 0 && llcieMask == 0;
    lp.enaLMask = locallab.spots.at(sp).enaLMask && lllogMask == 0 && llColorMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llColorMask == 0 && llsoftMask == 0 && lllcMask == 0 && llsharMask == 0 && llExpMask == 0 && llSHMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && ll_Mask == 0 && llcieMask == 0;// Exposure mask is deactivated if Color & Light mask is visible
    lp.enacieMask = locallab.spots.at(sp).enacieMask && llcieMask == 0 && llColorMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llColorMask == 0 && llSHMask == 0 && llsoftMask == 0 && lllcMask == 0 && llsharMask == 0 && llExpMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && lllogMask == 0  && ll_Mask == 0;


    lp.thrlow = locallab.spots.at(sp).levelthrlow;
    lp.thrhigh = locallab.spots.at(sp).levelthr;
    lp.usemask = locallab.spots.at(sp).usemask;
    lp.lnoiselow = locallab.spots.at(sp).lnoiselow;

    //  printf("llColorMask=%i lllcMask=%i llExpMask=%i  llSHMask=%i llcbMask=%i llretiMask=%i lltmMask=%i llblMask=%i llvibMask=%i\n", llColorMask, lllcMask, llExpMask, llSHMask, llcbMask, llretiMask, lltmMask, llblMask, llvibMask);
    if (locallab.spots.at(sp).softMethod == "soft") {
        lp.softmet = 0;
    } else if (locallab.spots.at(sp).softMethod == "reti") {
        lp.softmet = 1;
    }

    if (locallab.spots.at(sp).blMethod == "blur") {
        lp.blmet = 0;
    } else if (locallab.spots.at(sp).blMethod == "med") {
        lp.blmet = 1;
    } else if (locallab.spots.at(sp).blMethod == "guid") {
        lp.blmet = 2;
    }

    if (locallab.spots.at(sp).chroMethod == "lum") {
        lp.chromet = 0;
    } else if (locallab.spots.at(sp).chroMethod == "chr") {
        lp.chromet = 1;
    } else if (locallab.spots.at(sp).chroMethod == "all") {
        lp.chromet = 2;
    }




    if (locallab.spots.at(sp).quamethod == "cons") {
        lp.quamet = 0;
    } else if (locallab.spots.at(sp).quamethod == "agre") {
        lp.quamet = 1;
    } else if (locallab.spots.at(sp).quamethod == "nlmean") {
        lp.quamet = 2;
    } else if (locallab.spots.at(sp).quamethod == "none") {
        lp.quamet = 3;
    }
//    printf("lpqualmet=%i\n", lp.quamet);

    if (locallab.spots.at(sp).shMethod == "std") {
        lp.shmeth = 0;
    } else if (locallab.spots.at(sp).shMethod == "tone") {
        lp.shmeth = 1;
    }


    if (locallab.spots.at(sp).medMethod == "none") {
        lp.medmet = -1;
    } else if (locallab.spots.at(sp).medMethod == "33") {
        lp.medmet = 0;
    } else if (locallab.spots.at(sp).medMethod == "55") {
        lp.medmet = 1;
    } else if (locallab.spots.at(sp).medMethod == "77") {
        lp.medmet = 2;
    } else if (locallab.spots.at(sp).medMethod == "99") {
        lp.medmet = 3;
    }
/*
    if (locallab.spots.at(sp).blurMethod == "norm") {
        lp.blurmet = 0;
    } else if (locallab.spots.at(sp).blurMethod == "inv") {
        lp.blurmet = 1;
    }
*/
    if (locallab.spots.at(sp).invbl == false) {
        lp.blurmet = 0;
    } else if (locallab.spots.at(sp).invbl == true) {
        lp.blurmet = 1;
    }

    if (locallab.spots.at(sp).invmask == false) {
        lp.invmask = false;
    } else if (locallab.spots.at(sp).invmask == true) {
        lp.invmask = true;
    }

    if (locallab.spots.at(sp).invmaskd == false) {
        lp.invmaskd = false;
    } else if (locallab.spots.at(sp).invmaskd == true) {
        lp.invmaskd = true;
    }


    if (locallab.spots.at(sp).showmaskblMethodtyp == "blur") {
        lp.smasktyp = 0;
    } else if (locallab.spots.at(sp).showmaskblMethodtyp == "nois") {
        lp.smasktyp = 1;
    } else if (locallab.spots.at(sp).showmaskblMethodtyp == "all") {
        lp.smasktyp = 2;
    }


    if (locallab.spots.at(sp).spotMethod == "norm") {
        lp.excmet = 0;
    } else if (locallab.spots.at(sp).spotMethod == "exc") {
        lp.excmet = 1;
    } else if (locallab.spots.at(sp).spotMethod == "full") {
        lp.excmet = 2;
    }

    if (locallab.spots.at(sp).merMethod == "mone") {
        lp.mergemet = 0;
    } else if (locallab.spots.at(sp).merMethod == "mtwo") {
        lp.mergemet = 1;
    } else if (locallab.spots.at(sp).merMethod == "mthr") {
        lp.mergemet = 2;
    } else if (locallab.spots.at(sp).merMethod == "mfou") {
        lp.mergemet = 3;
    } else if (locallab.spots.at(sp).merMethod == "mfiv") {
        lp.mergemet = 4;
    }

    if (locallab.spots.at(sp).mergecolMethod == "one") {
        lp.mergecolMethod = 0;
    } else if (locallab.spots.at(sp).mergecolMethod == "two") {
        lp.mergecolMethod = 1;
    } else if (locallab.spots.at(sp).mergecolMethod == "thr") {
        lp.mergecolMethod = 2;
    } else if (locallab.spots.at(sp).mergecolMethod == "fou") {
        lp.mergecolMethod = 3;
    } else if (locallab.spots.at(sp).mergecolMethod == "fiv") {
        lp.mergecolMethod = 4;
    } else if (locallab.spots.at(sp).mergecolMethod == "six") {
        lp.mergecolMethod = 5;
    } else if (locallab.spots.at(sp).mergecolMethod == "sev") {
        lp.mergecolMethod = 6;
    } else if (locallab.spots.at(sp).mergecolMethod == "sev0") {
        lp.mergecolMethod = 7;
    } else if (locallab.spots.at(sp).mergecolMethod == "sev1") {
        lp.mergecolMethod = 8;
    } else if (locallab.spots.at(sp).mergecolMethod == "sev2") {
        lp.mergecolMethod = 9;
    } else if (locallab.spots.at(sp).mergecolMethod == "hei") {
        lp.mergecolMethod = 10;
    } else if (locallab.spots.at(sp).mergecolMethod == "nin") {
        lp.mergecolMethod = 11;
    } else if (locallab.spots.at(sp).mergecolMethod == "ten") {
        lp.mergecolMethod = 12;
    } else if (locallab.spots.at(sp).mergecolMethod == "ele") {
        lp.mergecolMethod = 13;
    } else if (locallab.spots.at(sp).mergecolMethod == "twe") {
        lp.mergecolMethod = 14;
    } else if (locallab.spots.at(sp).mergecolMethod == "thi") {
        lp.mergecolMethod = 15;
    } else if (locallab.spots.at(sp).mergecolMethod == "for") {
        lp.mergecolMethod = 16;
    } else if (locallab.spots.at(sp).mergecolMethod == "hue") {
        lp.mergecolMethod = 17;
    } else if (locallab.spots.at(sp).mergecolMethod == "sat") {
        lp.mergecolMethod = 18;
    } else if (locallab.spots.at(sp).mergecolMethod == "col") {
        lp.mergecolMethod = 19;
    } else if (locallab.spots.at(sp).mergecolMethod == "lum") {
        lp.mergecolMethod = 20;
    }

    if (locallab.spots.at(sp).localedgMethod == "fir") {
        lp.edgwmet = 0;
    } else if (locallab.spots.at(sp).localedgMethod == "sec") {
        lp.edgwmet = 1;
    } else if (locallab.spots.at(sp).localedgMethod == "thr") {
        lp.edgwmet = 2;
    }

    if (locallab.spots.at(sp).localneiMethod == "none") {
        lp.neiwmet = -1;
        lp.lip3 = false;
    } else if (locallab.spots.at(sp).localneiMethod == "low") {
        lp.neiwmet = 0;
        lp.lip3 = true;
    } else if (locallab.spots.at(sp).localneiMethod == "high") {
        lp.lip3 = true;
        lp.neiwmet = 1;
    }


    if (locallab.spots.at(sp).wavMethod == "D2") {
        lp.daubLen = 4;
    } else if (locallab.spots.at(sp).wavMethod == "D4") {
        lp.daubLen = 6;
    } else if (locallab.spots.at(sp).wavMethod == "D6") {
        lp.daubLen = 8;
    } else if (locallab.spots.at(sp).wavMethod == "D10") {
        lp.daubLen = 12;
    } else if (locallab.spots.at(sp).wavMethod == "D14") {
        lp.daubLen = 16;
    }


    lp.edgwena = locallab.spots.at(sp).wavedg;

    lp.opacol = 0.01 * locallab.spots.at(sp).opacol;

    if (locallab.spots.at(sp).shape == "ELI") {
        lp.shapmet = 0;
    } else /*if (locallab.spots.at(sp).shape == "RECT")*/ {
        lp.shapmet = 1;
    }

    lp.denoiena = locallab.spots.at(sp).expblur;

    bool wavcurveden = false;
    float local_noiself = 0.f;
    float local_noiself0 = 0.f;
    float local_noiself2 = 0.f;
    float local_noiselc = 0.f;
    float lnoiselc4 = 0.f;
    float lnoiselc5 = 0.f;
    float lnoiselc6 = 0.f;

    if (locwavCurveden && locwavdenutili) {
        if (lp.denoiena) {
            for (int i = 0; i < 500; i++) {
                if (locwavCurveden[i] != 0.f) {
                    wavcurveden = true;
                }
            }
        }
    }

    if (wavcurveden) {
        if (lp.denoiena) {
            local_noiself0 = 250.f * locwavCurveden[0];
            local_noiself = 250.f * locwavCurveden[83];
            local_noiself2 = 250.f * locwavCurveden[166];
            local_noiselc = 200.f * locwavCurveden[250];
            lnoiselc4 = 250.f * locwavCurveden[333];
            lnoiselc5 = 250.f * locwavCurveden[416];
            lnoiselc6 = 250.f * locwavCurveden[500];
       }
    }
    lp.wavcurvedenoi = wavcurveden;
    float local_noiseldetail = (float)locallab.spots.at(sp).noiselumdetail;
    int local_noiselequal = locallab.spots.at(sp).noiselequal;
    float local_noisechrodetail = (float)locallab.spots.at(sp).noisechrodetail;
    int local_sensiden = locallab.spots.at(sp).sensiden;
    float local_reparden = locallab.spots.at(sp).reparden;
    float local_repartm = locallab.spots.at(sp).repartm;
    float local_detailthr = (float)locallab.spots.at(sp).detailthr;
    float local_recothr = (float)locallab.spots.at(sp).recothres;
    float local_lowthr = (float)locallab.spots.at(sp).lowthres;
    float local_higthr = (float)locallab.spots.at(sp).higthres;
    float local_recothrd = (float)locallab.spots.at(sp).recothresd;
    float local_lowthrd = (float)locallab.spots.at(sp).lowthresd;
    float local_midthrd = (float)locallab.spots.at(sp).midthresd;
    float local_midthrdch = (float)locallab.spots.at(sp).midthresdch;
    float local_higthrd = (float)locallab.spots.at(sp).higthresd;
    float local_decayd = (float)locallab.spots.at(sp).decayd;
    float local_recothrc = (float)locallab.spots.at(sp).recothresc;
    float local_lowthrc = (float)locallab.spots.at(sp).lowthresc;
    float local_higthrc = (float)locallab.spots.at(sp).higthresc;
    float local_decayc = (float)locallab.spots.at(sp).decayc;

    float local_recothre = (float)locallab.spots.at(sp).recothrese;
    float local_lowthre = (float)locallab.spots.at(sp).lowthrese;
    float local_higthre = (float)locallab.spots.at(sp).higthrese;
    float local_decaye = (float)locallab.spots.at(sp).decaye;

    float local_recothrv = (float)locallab.spots.at(sp).recothresv;
    float local_lowthrv = (float)locallab.spots.at(sp).lowthresv;
    float local_higthrv = (float)locallab.spots.at(sp).higthresv;
    float local_decayv = (float)locallab.spots.at(sp).decayv;

    float local_recothrcb = (float)locallab.spots.at(sp).recothrescb;
    float local_lowthrcb = (float)locallab.spots.at(sp).lowthrescb;
    float local_higthrcb = (float)locallab.spots.at(sp).higthrescb;
    float local_decaycb = (float)locallab.spots.at(sp).decaycb;

    float local_recothrr = (float)locallab.spots.at(sp).recothresr;
    float local_lowthrr = (float)locallab.spots.at(sp).lowthresr;
    float local_higthrr = (float)locallab.spots.at(sp).higthresr;
    float local_decayr = (float)locallab.spots.at(sp).decayr;

    float local_recothrt = (float)locallab.spots.at(sp).recothrest;
    float local_lowthrt = (float)locallab.spots.at(sp).lowthrest;
    float local_higthrt = (float)locallab.spots.at(sp).higthrest;
    float local_decayt = (float)locallab.spots.at(sp).decayt;

    float local_recothrw = (float)locallab.spots.at(sp).recothresw;
    float local_lowthrw = (float)locallab.spots.at(sp).lowthresw;
    float local_higthrw = (float)locallab.spots.at(sp).higthresw;
    float local_decayw = (float)locallab.spots.at(sp).decayw;

    float local_recothrs = (float)locallab.spots.at(sp).recothress;
    float local_lowthrs = (float)locallab.spots.at(sp).lowthress;
    float local_higthrs = (float)locallab.spots.at(sp).higthress;
    float local_decays = (float)locallab.spots.at(sp).decays;

    float local_recothrcie = (float)locallab.spots.at(sp).recothrescie;
    float local_lowthrcie = (float)locallab.spots.at(sp).lowthrescie;
    float local_higthrcie = (float)locallab.spots.at(sp).higthrescie;
    float local_decaycie = (float)locallab.spots.at(sp).decaycie;

    float local_recothrl = (float)locallab.spots.at(sp).recothresl;
    float local_lowthrl = (float)locallab.spots.at(sp).lowthresl;
    float local_higthrl = (float)locallab.spots.at(sp).higthresl;
    float local_decayl = (float)locallab.spots.at(sp).decayl;

    float local_noisecf = ((float)locallab.spots.at(sp).noisechrof) / 10.f;
    float local_noisecc = ((float)locallab.spots.at(sp).noisechroc) / 10.f;
    float multi[6];

    for (int y = 0; y < 6; y++) {
        multi[y] = ((float) locallab.spots.at(sp).mult[y]);
    }

    float multish[5];

    for (int y = 0; y < 5; y++) {
        multish[y] = ((float) locallab.spots.at(sp).multsh[y]);
    }

    float thresho = ((float)locallab.spots.at(sp).threshold);
    float chromcbdl = (float)locallab.spots.at(sp).chromacbdl;

    int local_chroma = locallab.spots.at(sp).chroma;
    int local_sensi = locallab.spots.at(sp).sensi;
    int local_sensibn = locallab.spots.at(sp).sensibn;
    int local_sensitm = locallab.spots.at(sp).sensitm;
    int local_sensiexclu = locallab.spots.at(sp).sensiexclu;
    float structexclude = (float) locallab.spots.at(sp).structexclu;
    int local_sensilc = locallab.spots.at(sp).sensilc;
    int local_warm = locallab.spots.at(sp).warm;
    int local_sensih = locallab.spots.at(sp).sensih;
    int local_dehaze = locallab.spots.at(sp).dehaz;
    int local_depth = locallab.spots.at(sp).depth;
    int local_dehazeSaturation = locallab.spots.at(sp).dehazeSaturation;
    int local_sensicb = locallab.spots.at(sp).sensicb;
    float local_clarityml = (float) locallab.spots.at(sp).clarityml;
    float local_contresid = (float) locallab.spots.at(sp).contresid;
    int local_contrast = locallab.spots.at(sp).contrast;
    float local_lightness = (float) locallab.spots.at(sp).lightness;
    float labgridALowloc = locallab.spots.at(sp).labgridALow;
    float labgridBLowloc = locallab.spots.at(sp).labgridBLow;
    float labgridBHighloc = locallab.spots.at(sp).labgridBHigh;
    float labgridAHighloc = locallab.spots.at(sp).labgridAHigh;
    float strengthgrid = (float) locallab.spots.at(sp).strengthgrid;
    float labgridBLowlocmerg = locallab.spots.at(sp).labgridBLowmerg;
    float labgridBHighlocmerg = locallab.spots.at(sp).labgridBHighmerg;
    float labgridALowlocmerg = locallab.spots.at(sp).labgridALowmerg;
    float labgridAHighlocmerg = locallab.spots.at(sp).labgridAHighmerg;
    float local_gamlc = (float) locallab.spots.at(sp).gamlc;
    float local_gamc = (float) locallab.spots.at(sp).gamc;
    float local_gamex = (float) locallab.spots.at(sp).gamex;

    float blendmasklc = ((float) locallab.spots.at(sp).blendmasklc) / 100.f ;
    float radmasklc = ((float) locallab.spots.at(sp).radmasklc);
    float chromasklc = ((float) locallab.spots.at(sp).chromasklc);
    float structcolor = (float) locallab.spots.at(sp).structcol;
    float blendmaskcolor = ((float) locallab.spots.at(sp).blendmaskcol) / 100.f ;
    float radmaskcolor = ((float) locallab.spots.at(sp).radmaskcol);
    float chromaskcolor = ((float) locallab.spots.at(sp).chromaskcol);
    float gammaskcolor = ((float) locallab.spots.at(sp).gammaskcol);
    float slomaskcolor = ((float) locallab.spots.at(sp).slomaskcol);
    float blendmaskexpo = ((float) locallab.spots.at(sp).blendmaskexp) / 100.f ;
    float radmaskexpo = ((float) locallab.spots.at(sp).radmaskexp);
    float chromaskexpo = ((float) locallab.spots.at(sp).chromaskexp);
    float gammaskexpo = ((float) locallab.spots.at(sp).gammaskexp);
    float slomaskexpo = ((float) locallab.spots.at(sp).slomaskexp);
    float strmaskexpo = ((float) locallab.spots.at(sp).strmaskexp);
    float angmaskexpo = ((float) locallab.spots.at(sp).angmaskexp);
    float strmask = ((float) locallab.spots.at(sp).str_mask);
    float angmask = ((float) locallab.spots.at(sp).ang_mask);
    float strexpo = ((float) locallab.spots.at(sp).strexp);
    float angexpo = ((float) locallab.spots.at(sp).angexp);
    float strSH = ((float) locallab.spots.at(sp).strSH);
    float angSH = ((float) locallab.spots.at(sp).angSH);
    float strcol = ((float) locallab.spots.at(sp).strcol);
    float strcolab = ((float) locallab.spots.at(sp).strcolab);
    float strcolh = ((float) locallab.spots.at(sp).strcolh);
    float angcol = ((float) locallab.spots.at(sp).angcol);
    float strvib = ((float) locallab.spots.at(sp).strvib);
    float strvibab = ((float) locallab.spots.at(sp).strvibab);
    float strvibh = ((float) locallab.spots.at(sp).strvibh);
    float angvib = ((float) locallab.spots.at(sp).angvib);
    float strwav = ((float) locallab.spots.at(sp).strwav);
    float angwav = ((float) locallab.spots.at(sp).angwav);
    float strlog = ((float) locallab.spots.at(sp).strlog);
    float anglog = ((float) locallab.spots.at(sp).anglog);
    float softradiusexpo = ((float) locallab.spots.at(sp).softradiusexp);
    float softradiuscolor = ((float) locallab.spots.at(sp).softradiuscol);
    float softradiusreti = ((float) locallab.spots.at(sp).softradiusret);
    float softradiustma = ((float) locallab.spots.at(sp).softradiustm);
    float softradiuscbdl = ((float) locallab.spots.at(sp).softradiuscb);
    float blendmaskSH = ((float) locallab.spots.at(sp).blendmaskSH) / 100.f ;
    float radmaskSH = ((float) locallab.spots.at(sp).radmaskSH);
    float chromaskSH = ((float) locallab.spots.at(sp).chromaskSH);
    float gammaskSH = ((float) locallab.spots.at(sp).gammaskSH);
    float slomaskSH = ((float) locallab.spots.at(sp).slomaskSH);
    float blendmaskvib = ((float) locallab.spots.at(sp).blendmaskvib) / 100.f ;
    float radmaskvib = ((float) locallab.spots.at(sp).radmaskvib);
    float chromaskvib = ((float) locallab.spots.at(sp).chromaskvib);
    float gammaskvib = ((float) locallab.spots.at(sp).gammaskvib);
    float slomaskvib = ((float) locallab.spots.at(sp).slomaskvib);
    float structexpo = (float) locallab.spots.at(sp).structexp;
    float blurexpo = (float) locallab.spots.at(sp).blurexpde;
    float blurcolor = (float) locallab.spots.at(sp).blurcolde;
    float blurcolmask = (float) locallab.spots.at(sp).blurcol;
    float contcolmask = (float) locallab.spots.at(sp).contcol;
    float blurSH = (float) locallab.spots.at(sp).blurSHde;
    float local_transit = locallab.spots.at(sp).transit;
    float local_feather = locallab.spots.at(sp).feather;
    float local_transitweak = (float)locallab.spots.at(sp).transitweak;
    float local_transitgrad = (float)locallab.spots.at(sp).transitgrad;
    float radius = (float) locallab.spots.at(sp).radius;
    int itera = locallab.spots.at(sp).itera;
    int guidbl = locallab.spots.at(sp).guidbl;
    float epsbl = (float) locallab.spots.at(sp).epsbl;
    float sharradius = LIM(locallab.spots.at(sp).sharradius, 0.42, 3.5);
    float lcamount = ((float) locallab.spots.at(sp).lcamount);
    lcamount = LIM01(lcamount); //to prevent crash with old pp3 integer
    float sharblurr = LIM(locallab.spots.at(sp).sharblur, 0.2, 3.); //to prevent crash with old pp3 integer
    int local_sensisha = locallab.spots.at(sp).sensisha;
    int local_sharamount = locallab.spots.at(sp).sharamount;
    int local_shardamping = locallab.spots.at(sp).shardamping;
    int local_shariter = locallab.spots.at(sp).shariter;
    bool inverse = locallab.spots.at(sp).invers;
    bool curvacti = locallab.spots.at(sp).curvactiv;
    bool acti = locallab.spots.at(sp).activlum;
    bool cupas = false; // Provision
    int local_sensisf = locallab.spots.at(sp).sensisf;
    bool inverseex = locallab.spots.at(sp).inversex;
    bool inversesh = locallab.spots.at(sp).inverssh;
    bool equiltm = locallab.spots.at(sp).equiltm;
    bool fftwlc = locallab.spots.at(sp).fftwlc;
    bool fftwreti = locallab.spots.at(sp).fftwreti;

    float blendmaskL = ((float) locallab.spots.at(sp).blendmaskL) / 100.f ;
    float radmaskL = ((float) locallab.spots.at(sp).radmaskL);
    float chromaskL = ((float) locallab.spots.at(sp).chromaskL);

    bool equilret = locallab.spots.at(sp).equilret;
    bool inverserad = false; // Provision
    bool inverseret = locallab.spots.at(sp).inversret;
    bool inversesha = locallab.spots.at(sp).inverssha;
    double strength = (double) locallab.spots.at(sp).strength;
    float str = (float)locallab.spots.at(sp).str;
    int scaleret = (float)locallab.spots.at(sp).scalereti;

    int local_sensihs = locallab.spots.at(sp).sensihs;
    int highhs = locallab.spots.at(sp).highlights;
    int hltonahs = locallab.spots.at(sp).h_tonalwidth;
    int shadhs = locallab.spots.at(sp).shadows;
    int shtonals = locallab.spots.at(sp).s_tonalwidth;
    int radhs = locallab.spots.at(sp).sh_radius;
    float blendmaskcb = ((float) locallab.spots.at(sp).blendmaskcb) / 100.f ;
    float radmaskcb = ((float) locallab.spots.at(sp).radmaskcb);
    float chromaskcb = ((float) locallab.spots.at(sp).chromaskcb);
    float gammaskcb = ((float) locallab.spots.at(sp).gammaskcb);
    float slomaskcb = ((float) locallab.spots.at(sp).slomaskcb);
    bool enaretiMasktm = locallab.spots.at(sp).enaretiMasktmap;
    lp.enaretiMasktmap =  enaretiMasktm;
    float blendmasktm = ((float) locallab.spots.at(sp).blendmasktm) / 100.f ;
    float radmasktm = ((float) locallab.spots.at(sp).radmasktm);
    float chromasktm = ((float) locallab.spots.at(sp).chromasktm);
    float gammasktm = ((float) locallab.spots.at(sp).gammasktm);
    float slomasktm = ((float) locallab.spots.at(sp).slomasktm);
    bool wavgradl = locallab.spots.at(sp).wavgradl;

    float blendmaskbl = ((float) locallab.spots.at(sp).blendmaskbl) / 100.f ;
    float radmaskbl = ((float) locallab.spots.at(sp).radmaskbl);
    float chromaskbl = ((float) locallab.spots.at(sp).chromaskbl);
    float gammaskbl = ((float) locallab.spots.at(sp).gammaskbl);
    float slomaskbl = ((float) locallab.spots.at(sp).slomaskbl);
    bool fftbl = locallab.spots.at(sp).fftwbl;


    lp.sourcegray = (float) locallab.spots.at(sp).sourceGray;
    lp.targetgray = (float) locallab.spots.at(sp).targetGray;
    lp.blackev = (float) locallab.spots.at(sp).blackEv;
    lp.whiteev  = (float) locallab.spots.at(sp).whiteEv;
    lp.detail = locallab.spots.at(sp).detail;
    lp.sensilog = locallab.spots.at(sp).sensilog;
    lp.Autogray = locallab.spots.at(sp).Autogray;
    lp.autocompute = locallab.spots.at(sp).autocompute;
    lp.baselog = (float) locallab.spots.at(sp).baselog;
    lp.sensimas = locallab.spots.at(sp).sensimask;
    lp.sensicie = locallab.spots.at(sp).sensicie;
    float blendmaskcie = ((float) locallab.spots.at(sp).blendmaskcie) / 100.f ;
    float radmaskcie = ((float) locallab.spots.at(sp).radmaskcie);
    float chromaskcie = ((float) locallab.spots.at(sp).chromaskcie);
    lp.deltaem = locallab.spots.at(sp).deltae;
    lp.scalereti = scaleret;
    lp.cir = circr;
    lp.recur = recur;
    lp.actsp = acti;
    lp.xc = w * local_center_x;
    lp.yc = h * local_center_y;
    lp.lx = w * local_x;
    lp.ly = h * local_y;
    lp.lxL = w * local_xL;
    lp.lyT = h * local_yT;
    lp.chro = local_chroma;
    lp.struco = structcolor;
    lp.strengrid = strengthgrid;
    lp.blendmalc = blendmasklc;
    lp.radmalc = radmasklc;
    lp.chromalc = chromasklc;
    lp.blendmacol = blendmaskcolor;
    lp.radmacol = radmaskcolor;
    lp.chromacol = chromaskcolor;
    lp.gammacol = gammaskcolor;
    lp.slomacol = slomaskcolor;
    lp.radmaexp = radmaskexpo;
    lp.chromaexp = chromaskexpo;
    lp.gammaexp = gammaskexpo;
    lp.slomaexp = slomaskexpo;
    lp.strmaexp = strmaskexpo;
    lp.angmaexp = angmaskexpo;
    lp.str_mas = strmask;
    lp.ang_mas = angmask;
    
    lp.strexp = strexpo;
    lp.angexp = angexpo;
    lp.strSH = strSH;
    lp.angSH = angSH;
    lp.strcol = strcol;
    lp.strcolab = strcolab;
    lp.strcolh = strcolh;
    lp.angcol = angcol;
    lp.strvib = strvib;
    lp.strvibab = strvibab;
    lp.strvibh = strvibh;
    lp.angvib = angvib;
    lp.strwav = strwav;
    lp.angwav = angwav;
    lp.strlog = strlog;
    lp.anglog = anglog;
    lp.softradiusexp = softradiusexpo;
    lp.softradiuscol = softradiuscolor;
    lp.softradiusret = softradiusreti;
    lp.softradiuscb = softradiuscbdl;
    lp.softradiustm = softradiustma;
    lp.struexc = structexclude;
    lp.blendmaexp = blendmaskexpo;
    lp.blendmaSH = blendmaskSH;
    lp.radmaSH = radmaskSH;
    lp.chromaSH = chromaskSH;
    lp.gammaSH = gammaskSH;
    lp.slomaSH = slomaskSH;
    lp.blendmavib = blendmaskvib;
    lp.radmavib = radmaskvib;
    lp.chromavib = chromaskvib;
    lp.gammavib = gammaskvib;
    lp.slomavib = slomaskvib;
    lp.blendmacb = blendmaskcb;
    lp.radmacb = radmaskcb;
    lp.chromacbm = chromaskcb;
    lp.gammacb = gammaskcb;
    lp.slomacb = slomaskcb;
    lp.blendmatm = blendmasktm;
    lp.radmatm = radmasktm;
    lp.chromatm = chromasktm;
    lp.gammatm = gammasktm;
    lp.slomatm = slomasktm;
    lp.wavgradl = wavgradl;
    lp.blendmaL = blendmaskL;
    lp.radmaL = radmaskL;
    lp.chromaL = chromaskL;

    lp.strengthw = ((float) locallab.spots.at(sp).strengthw);
    lp.radiusw = ((float) locallab.spots.at(sp).radiusw);
    lp.detailw = ((float) locallab.spots.at(sp).detailw);
    lp.gradw = ((float) locallab.spots.at(sp).gradw);
    lp.tloww = ((float) locallab.spots.at(sp).tloww);
    lp.thigw = ((float) locallab.spots.at(sp).thigw);
    lp.edgw = ((float) locallab.spots.at(sp).edgw);
    lp.basew = ((float) locallab.spots.at(sp).basew);

    lp.blendmabl = blendmaskbl;
    lp.radmabl = radmaskbl;
    lp.chromabl = chromaskbl;
    lp.gammabl = gammaskbl;
    lp.slomabl = slomaskbl;
    lp.fftbl = fftbl;
    lp.it = itera;
    lp.guidb = guidbl;
    lp.strbl = 0.01f * (float) locallab.spots.at(sp).strbl;

    lp.epsb = epsbl;
    lp.struexp = structexpo;
    lp.blurexp = blurexpo;
    lp.blurcol = blurcolor;
    lp.blurcolmask = blurcolmask;
    lp.contcolmask = 0.01f * contcolmask;
    lp.blurSH = blurSH;
    lp.sens = local_sensi;
    lp.sensh = local_sensih;
    lp.dehaze = local_dehaze;
    lp.dehazeSaturation = local_dehazeSaturation;
    lp.depth = local_depth;
    lp.senscb = local_sensicb;
    lp.clarityml = local_clarityml;
    lp.contresid = local_contresid;
    lp.cont = local_contrast;
    lp.ligh = local_lightness;
    lp.lowA = labgridALowloc;
    lp.lowB = labgridBLowloc;
    lp.highB = labgridBHighloc;
    lp.highA = labgridAHighloc;
    lp.lowBmerg = labgridBLowlocmerg;
    lp.highBmerg = labgridBHighlocmerg;
    lp.lowAmerg = labgridALowlocmerg;
    lp.highAmerg = labgridAHighlocmerg;
    lp.gamlc = local_gamlc;
    lp.gamc = local_gamc;
    lp.gamex = local_gamex;

    lp.senssf = local_sensisf;
    lp.strng = strlight;
    lp.neig = neigh;
    lp.lap = laplac;

    if (lp.ligh >= -2.f && lp.ligh <= 2.f) {
        lp.ligh /= 5.f;
    }

    lp.trans = local_transit;
    lp.feath = local_feather;
    lp.transweak = local_transitweak;
    lp.transgrad = local_transitgrad;
    lp.rad = radius;
    lp.stren = strength;
    lp.sensbn = local_sensibn;
    lp.sensexclu = local_sensiexclu;
    lp.senslc = local_sensilc;
    lp.lcamount = lcamount;
    lp.inv = inverse;
    lp.invex = inverseex;
    lp.invsh = inversesh;
    lp.curvact = curvacti;
    lp.invrad = inverserad;
    lp.invret = inverseret;
    lp.equret = equilret;
    lp.equtm = equiltm;
    lp.invshar = inversesha;
    lp.str = str;
    lp.shrad = sharradius;
    lp.shblurr = sharblurr;
    lp.senssha = local_sensisha;
    lp.shamo = local_sharamount;
    lp.shdamp = local_shardamping;
    lp.shiter = local_shariter;
    lp.iterat = iterati;
    lp.balance = balanc;
    lp.balanceh = balanch;
    lp.colorde = colorde;
    lp.thr = thre;
    lp.stru = strucc;
    lp.noiself = local_noiself;
    lp.noiself0 = local_noiself0;
    lp.noiself2 = local_noiself2;
    lp.noiseldetail = local_noiseldetail;
    lp.detailthr = local_detailthr;
    lp.recothr = local_recothr;
    lp.lowthr = local_lowthr;
    lp.higthr = local_higthr;
    lp.recothrd = local_recothrd;
    lp.midthrd = local_midthrd;
    lp.midthrdch = local_midthrdch;
    lp.lowthrd = local_lowthrd;
    lp.higthrd = local_higthrd;
    lp.decayd = local_decayd;
    lp.recothrc = local_recothrc;
    lp.lowthrc = local_lowthrc;
    lp.higthrc = local_higthrc;
    lp.decayc = local_decayc;
    lp.recothre = local_recothre;
    lp.lowthre = local_lowthre;
    lp.higthre = local_higthre;
    lp.decaye = local_decaye;
    lp.recothrs = local_recothrs;
    lp.lowthrs = local_lowthrs;
    lp.higthrs = local_higthrs;
    lp.decays = local_decays;

    lp.recothrcie = local_recothrcie;
    lp.lowthrcie = local_lowthrcie;
    lp.higthrcie = local_higthrcie;
    lp.decaycie = local_decaycie;

    lp.recothrv = local_recothrv;
    lp.lowthrv = local_lowthrv;
    lp.higthrv = local_higthrv;
    lp.decayv = local_decayv;
    lp.recothrw = local_recothrw;
    lp.lowthrw = local_lowthrw;
    lp.higthrw = local_higthrw;
    lp.decayw = local_decayw;
    lp.recothrt = local_recothrt;
    lp.lowthrt = local_lowthrt;
    lp.higthrt = local_higthrt;
    lp.decayt = local_decayt;
    lp.recothrcb = local_recothrcb;
    lp.lowthrcb = local_lowthrcb;
    lp.higthrcb = local_higthrcb;
    lp.decaycb = local_decaycb;
    lp.recothrr = local_recothrr;
    lp.lowthrr = local_lowthrr;
    lp.higthrr = local_higthrr;
    lp.decayr = local_decayr;
    
    lp.recothrl = local_recothrl;
    lp.lowthrl = local_lowthrl;
    lp.higthrl = local_higthrl;
    lp.decayl = local_decayl;
    lp.noiselequal = local_noiselequal;
    lp.noisechrodetail = local_noisechrodetail;
    lp.noiselc = local_noiselc;
    lp.noiselc4 = lnoiselc4;
    lp.noiselc5 = lnoiselc5;
    lp.noiselc6 = lnoiselc6;
    
    lp.noisecf = local_noisecf;
    lp.noisecc = local_noisecc;
    lp.sensden = local_sensiden;
    lp.reparden = local_reparden;
    lp.repartm = local_repartm;
    lp.bilat = locallab.spots.at(sp).bilateral;
    lp.nldet = locallab.spots.at(sp).nldet;
    lp.nlstr = locallab.spots.at(sp).nlstr;
    lp.nlpat = locallab.spots.at(sp).nlpat;
    lp.nlrad = locallab.spots.at(sp).nlrad;
    lp.nlgam = locallab.spots.at(sp).nlgam;
    lp.noisegam = locallab.spots.at(sp).noisegam;
    lp.adjch = (float) locallab.spots.at(sp).adjblur;
    lp.strengt = streng;
    lp.gamm = gam;
    lp.esto = est;
    lp.scalt = scal_tm;
    lp.rewe = rewe;
    lp.senstm = local_sensitm;
    lp.amo = amo;
    lp.blurma = (float) locallab.spots.at(sp).blurmask;
    lp.fftma = locallab.spots.at(sp).fftmask;
    lp.contma = (float) locallab.spots.at(sp).contmask;

    lp.blendmacie = blendmaskcie;
    lp.radmacie = radmaskcie;
    lp.chromacie = chromaskcie;
    lp.denoichmask = locallab.spots.at(sp).denoichmask;

    for (int y = 0; y < 6; y++) {
        lp.mulloc[y] = LIM(multi[y], 0.f, 4.f);//to prevent crash with old pp3 integer
    }

    for (int y = 0; y < 5; y++) {
        lp.mullocsh[y] = multish[y];
    }
    lp.activspot = locallab.spots.at(sp).activ;
 

    lp.detailsh = locallab.spots.at(sp).detailSH; 
    lp.tePivot = locallab.spots.at(sp).tePivot;
    lp.threshol = thresho;
    lp.chromacb = chromcbdl;
    lp.expvib = locallab.spots.at(sp).expvibrance && lp.activspot ;
    lp.colorena = locallab.spots.at(sp).expcolor && lp.activspot && llExpMaskinv == 0 && llSHMaskinv == 0 && llExpMask == 0 && llsoftMask == 0 && llSHMask == 0 && llcbMask == 0 && lllcMask == 0 && llsharMask == 0 && llretiMask == 0 && lltmMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0; // Color & Light tool is deactivated if Exposure mask is visible or SHMask
    lp.blurena = locallab.spots.at(sp).expblur && lp.activspot && llColorMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llExpMask == 0 && llsoftMask == 0 && llSHMask == 0 && llcbMask == 0 && lllcMask == 0 && llsharMask == 0 && llretiMask == 0 && llColorMask == 0 && lltmMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;
    lp.tonemapena = locallab.spots.at(sp).exptonemap && lp.activspot && llColorMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llExpMask == 0 && llsoftMask == 0 && llSHMask == 0 && llcbMask == 0 && lllcMask == 0 && llsharMask == 0 && llretiMask == 0 && llColorMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;
    lp.retiena = locallab.spots.at(sp).expreti && lp.activspot && llColorMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llExpMask == 0 && llsoftMask == 0 && llSHMask == 0 && llcbMask == 0 && lllcMask == 0 && llsharMask == 0 && llColorMask == 0 && lltmMask == 0 && llvibMask == 0 && llSHMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;
    lp.lcena = locallab.spots.at(sp).expcontrast && lp.activspot && llColorMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llExpMask == 0 && llsoftMask == 0 && llSHMask == 0 && llcbMask == 0 && llsharMask == 0 && llColorMask == 0 && lltmMask == 0 && llvibMask == 0 && llSHMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;
    lp.cbdlena = locallab.spots.at(sp).expcbdl && lp.activspot && llColorMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llExpMask == 0 && llsoftMask == 0 && llSHMask == 0 && llretiMask == 0 && lllcMask == 0 && llsharMask == 0 && lllcMask == 0 && llColorMask == 0 && lltmMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;
    lp.exposena = locallab.spots.at(sp).expexpose  && lp.activspot  && llColorMaskinv == 0 && llSHMaskinv == 0 && llColorMask == 0 && llsoftMask == 0 && llSHMask == 0 && lllcMask == 0 && llsharMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0; // Exposure tool is deactivated if Color & Light mask SHmask is visible
    lp.hsena = locallab.spots.at(sp).expshadhigh && lp.activspot && llColorMaskinv == 0 && llExpMaskinv == 0 && llColorMask == 0 && llsoftMask == 0 && llExpMask == 0 && llcbMask == 0 && lllcMask == 0 && llsharMask == 0 && llretiMask == 0 && lltmMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;// Shadow Highlight tool is deactivated if Color & Light mask or SHmask is visible
    lp.vibena = locallab.spots.at(sp).expvibrance && lp.activspot && llColorMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llColorMask == 0 && llsoftMask == 0 && llExpMask == 0 && llcbMask == 0 && lllcMask == 0 && llsharMask == 0 && llretiMask == 0 && llcbMask == 0 && lltmMask == 0 && llSHMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;// vibrance tool is deactivated if Color & Light mask or SHmask is visible
    lp.sharpena = locallab.spots.at(sp).expsharp && lp.activspot && llColorMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llColorMask == 0 && llsoftMask == 0 && llExpMask == 0 && llcbMask == 0 && lllcMask == 0 && llretiMask == 0 && llcbMask == 0 && lltmMask == 0 && llSHMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;
    lp.sfena = locallab.spots.at(sp).expsoft && lp.activspot && llColorMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llColorMask == 0 && llExpMask == 0 && llcbMask == 0 && lllcMask == 0 && llretiMask == 0 && llcbMask == 0 && lltmMask == 0 && llSHMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0 && llcieMask == 0;
    lp.maskena = locallab.spots.at(sp).expmask && lp.activspot && llColorMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llColorMask == 0 && llsoftMask == 0 && llExpMask == 0 && llcbMask == 0 && lllcMask == 0 && llsharMask == 0 && llretiMask == 0 && llcbMask == 0 && lltmMask == 0 && lllogMask == 0 && llSHMask == 0 && llcieMask == 0;// vibrance tool is deactivated if Color & Light mask or SHmask is visible
    lp.logena = locallab.spots.at(sp).explog && lp.activspot && llColorMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llColorMask == 0 && llsoftMask == 0 && llExpMask == 0 && llcbMask == 0 && lllcMask == 0 && llsharMask == 0 && llretiMask == 0 && llcbMask == 0 && lltmMask == 0 && llSHMask == 0 && ll_Mask == 0 && llcieMask == 0;// vibrance tool is deactivated if Color & Light mask or SHmask is visible
    lp.cieena = locallab.spots.at(sp).expcie && lp.activspot && llColorMaskinv == 0 && llExpMaskinv == 0 && llSHMaskinv == 0 && llColorMask == 0 && llsoftMask == 0 && llExpMask == 0 && llcbMask == 0 && lllcMask == 0 && llsharMask == 0 && llretiMask == 0 && lltmMask == 0 && llvibMask == 0 && lllogMask == 0 && ll_Mask == 0;// Shadow Highlight tool is deactivated if Color & Light mask or SHmask is visible


    lp.islocal = (lp.expvib || lp.colorena || lp.blurena || lp.tonemapena || lp.retiena || lp.lcena || lp.cbdlena || lp.exposena || lp.hsena || lp.vibena || lp.sharpena || lp.sfena || lp.maskena || lp.logena || lp.cieena);

    lp.sensv = local_sensiv;
    lp.past =  chromaPastel;
    lp.satur = chromaSatur;

    lp.cut_past = cupas;
    lp.blac = locallab.spots.at(sp).black;
    lp.shcomp = locallab.spots.at(sp).shcompr;
    lp.shadex = locallab.spots.at(sp).shadex;
    lp.hlcomp = locallab.spots.at(sp).hlcompr;
    lp.hlcompthr = locallab.spots.at(sp).hlcomprthresh;
    lp.expcomp = LIM(locallab.spots.at(sp).expcomp, -2.0, 4.0); //to prevent crash with Old pp3 with integer
    lp.expchroma = locallab.spots.at(sp).expchroma / 100.;
    lp.sensex = local_sensiex;
    lp.war = local_warm;
    lp.highlihs = highhs;
    lp.shadowhs = shadhs;
    lp.radiushs = radhs;
    lp.hltonalhs = hltonahs;
    lp.shtonalhs = shtonals;
    lp.senshs = local_sensihs;
    lp.ftwlc = fftwlc;
    lp.ftwreti = fftwreti;
    lp.sigmadr = locallab.spots.at(sp).sigmadr;
    lp.sigmabl = locallab.spots.at(sp).sigmabl;
    lp.sigmaed = locallab.spots.at(sp).sigmaed;
    lp.sigmalc = locallab.spots.at(sp).sigmalc;
    lp.sigmalc2 = locallab.spots.at(sp).sigmalc2;
    lp.residsha = locallab.spots.at(sp).residsha;
    lp.residshathr = locallab.spots.at(sp).residshathr;
    lp.residhi = locallab.spots.at(sp).residhi;
    lp.residhithr = locallab.spots.at(sp).residhithr;
    lp.residgam = locallab.spots.at(sp).residgam;
    lp.residslop = locallab.spots.at(sp).residslop;
    lp.blwh = locallab.spots.at(sp).blwh;
    lp.senscolor = (int) locallab.spots.at(sp).colorscope;
    //replace scope color vibrance shadows
    lp.sens = lp.senscolor;
    lp.sensv = lp.senscolor;
    lp.senshs = lp.senscolor;

    lp.mLjz = locallab.spots.at(sp).clarilresjz / 100.0;
    lp.mCjz = locallab.spots.at(sp).claricresjz / 100.0;
    lp.softrjz = locallab.spots.at(sp).clarisoftjz;

}

static void calcTransitionrect(const float lox, const float loy, const float ach, const local_params& lp, int &zone, float &localFactor)
{
    zone = 0;

    if (lox >= lp.xc && lox < lp.xc + lp.lx) {
        if (loy >= lp.yc && loy < lp.yc + lp.ly) {
            if (lox < lp.xc + lp.lx * ach && loy < lp.yc + lp.ly * ach) {
                zone = 2;
            } else {
                zone = 1;
                localFactor = pow_F(calcLocalFactorrect(lox, loy, lp.xc, lp.lx, lp.yc, lp.ly, ach, lp.transgrad), lp.transweak);
            }
        } else if (loy < lp.yc && loy > lp.yc - lp.lyT) {
            if (lox < lp.xc + lp.lx * ach && loy > lp.yc - lp.lyT * ach) {
                zone = 2;
            } else {
                zone = 1;
                localFactor = pow_F(calcLocalFactorrect(lox, loy, lp.xc, lp.lx, lp.yc, lp.lyT, ach, lp.transgrad), lp.transweak);
            }
        }
    } else if (lox < lp.xc && lox > lp.xc - lp.lxL) {
        if (loy <= lp.yc && loy > lp.yc - lp.lyT) {
            if (lox > (lp.xc - lp.lxL * ach) && loy > (lp.yc - lp.lyT * ach)) {
                zone = 2;
            } else {
                zone = 1;
                localFactor = pow_F(calcLocalFactorrect(lox, loy, lp.xc, lp.lxL, lp.yc, lp.lyT, ach, lp.transgrad), lp.transweak);
            }
        } else if (loy > lp.yc && loy < lp.yc + lp.ly) {
            if (lox > (lp.xc - lp.lxL * ach) && loy < (lp.yc + lp.ly * ach)) {
                zone = 2;
            } else {
                zone = 1;
                localFactor = pow_F(calcLocalFactorrect(lox, loy, lp.xc, lp.lxL, lp.yc, lp.ly, ach, lp.transgrad), lp.transweak);
            }
        }
    }
}

static void calcTransition(const float lox, const float loy, const float ach, const local_params& lp, int &zone, float &localFactor)
{
    // returns the zone (0 = outside selection, 1 = transition zone between outside and inside selection, 2 = inside selection)
    // and a factor to calculate the transition in case zone == 1

    zone = 0;

    if (lox >= lp.xc && lox < lp.xc + lp.lx) {
        if (loy >= lp.yc && loy < lp.yc + lp.ly) {
            const float zoneVal = SQR((lox - lp.xc) / (ach * lp.lx)) + SQR((loy - lp.yc) / (ach * lp.ly));
            zone = zoneVal < 1.f ? 2 : 0;

            if (!zone) {
                zone = (zoneVal > 1.f && ((SQR((lox - lp.xc) / (lp.lx)) + SQR((loy - lp.yc) / (lp.ly))) < 1.f)) ? 1 : 0;
                if (zone == 1) {
                    localFactor = pow_F(calcLocalFactor(lox, loy, lp.xc, lp.lx, lp.yc, lp.ly, ach, lp.transgrad), lp.transweak);
                }
            }
        } else if (loy < lp.yc && loy > lp.yc - lp.lyT) {
            const float zoneVal = SQR((lox - lp.xc) / (ach * lp.lx)) + SQR((loy - lp.yc) / (ach * lp.lyT));
            zone = zoneVal < 1.f ? 2 : 0;

            if (!zone) {
                zone = (zoneVal > 1.f && ((SQR((lox - lp.xc) / (lp.lx)) + SQR((loy - lp.yc) / (lp.lyT))) < 1.f)) ? 1 : 0;
                if (zone == 1) {
                    localFactor = pow_F(calcLocalFactor(lox, loy, lp.xc, lp.lx, lp.yc, lp.lyT, ach, lp.transgrad), lp.transweak);
                }
            }
        }
    } else if (lox < lp.xc && lox > lp.xc - lp.lxL) {
        if (loy <= lp.yc && loy > lp.yc - lp.lyT) {
            const float zoneVal = SQR((lox - lp.xc) / (ach * lp.lxL)) + SQR((loy - lp.yc) / (ach * lp.lyT));
            zone = zoneVal < 1.f ? 2 : 0;

            if (!zone) {
                zone = (zoneVal > 1.f && ((SQR((lox - lp.xc) / (lp.lxL)) + SQR((loy - lp.yc) / (lp.lyT))) < 1.f)) ? 1 : 0;
                if (zone == 1) {
                    localFactor = pow_F(calcLocalFactor(lox, loy, lp.xc, lp.lxL, lp.yc, lp.lyT, ach, lp.transgrad), lp.transweak);
                }
            }
        } else if (loy > lp.yc && loy < lp.yc + lp.ly) {
            const float zoneVal = SQR((lox - lp.xc) / (ach * lp.lxL)) + SQR((loy - lp.yc) / (ach * lp.ly));
            zone = zoneVal < 1.f ? 2 : 0;

            if (!zone) {
                zone = (zoneVal > 1.f && ((SQR((lox - lp.xc) / (lp.lxL)) + SQR((loy - lp.yc) / (lp.ly))) < 1.f)) ? 1 : 0;
                if (zone == 1) {
                    localFactor = pow_F(calcLocalFactor(lox, loy, lp.xc, lp.lxL, lp.yc, lp.ly, ach, lp.transgrad), lp.transweak);
                }
            }
        }
    }
}

// Copyright 2018 Alberto Griggio <alberto.griggio@gmail.com>

float find_gray(float source_gray, float target_gray)
{
    // find a base such that log2lin(base, source_gray) = target_gray
    // log2lin is (base^source_gray - 1) / (base - 1), so we solve
    //
    //  (base^source_gray - 1) / (base - 1) = target_gray, that is
    //
    //  base^source_gray - 1 - base * target_gray + target_gray = 0
    //
    // use a bisection method (maybe later change to Netwon)

    if (source_gray <= 0.f) {
        return 0.f;
    }

    const auto f =
    [ = ](float x) -> float {
        return std::pow(x, source_gray) - 1.f - target_gray * x + target_gray;
    };

    // first find the interval we are interested in

    float lo = 1.f;

    while (f(lo) <= 0.f) {
        lo *= 2.f;
    }

    float hi = lo * 2.f;

    while (f(hi) >= 0.f) {
        hi *= 2.f;
    }

    if (std::isinf(hi)) {
        return 0.f;
    }

    // now search for a zero
    for (int iter = 0; iter < 100; ++iter) {
        float mid = lo + (hi - lo) / 2.f;
        float v = f(mid);

        if (std::abs(v) < 1e-4f || (hi - lo) / lo <= 1e-4f) {
            return mid;
        }

        if (v > 0.f) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    return 0.f; // not found
}

void ImProcFunctions::mean_sig (const float* const * const savenormL, float &meanf, float &stdf, int xStart, int xEnd, int yStart, int yEnd) const {
    const int size = (yEnd - yStart) * (xEnd - xStart);
    // use double precision for large accumulations
    double meand = 0.0;
    double stdd = 0.0;
#ifdef _OPENMP
    #pragma omp parallel for reduction(+:meand, stdd) if(multiThread)
#endif
    for (int y = yStart; y < yEnd; ++y) {
        for (int x = xStart; x < xEnd; ++x) {
            meand += static_cast<double>(savenormL[y][x]);
            stdd += SQR(static_cast<double>(savenormL[y][x]));
        }
    }
    meand /= size;
    stdd /= size;
    stdd -= SQR(meand);
    stdf = std::sqrt(stdd);
    meanf = meand;
}
// taken from darktable
inline float power_norm(float r, float g, float b)
{
    r = std::abs(r);
    g = std::abs(g);
    b = std::abs(b);

    float r2 = SQR(r);
    float g2 = SQR(g);
    float b2 = SQR(b);
    float d = r2 + g2 + b2;
    float n = r*r2 + g*g2 + b*b2;

    return n / std::max(d, 1e-12f);
}

inline float ev2gray(float ev)
{
    return std::pow(2.f, -ev + std::log2(0.18f));
}


inline float gray2ev(float gray)
{
    return std::log2(0.18f / gray);
}


inline float norm2(float r, float g, float b, TMatrix ws)
{
    return (power_norm(r, g, b) + Color::rgbLuminance(r, g, b, ws)) / 2.f;
}

inline float norm(float r, float g, float b, TMatrix ws)
{
    return (Color::rgbLuminance(r, g, b, ws));
}


// basic log encoding taken from ACESutil.Lin_to_Log2, from
// https://github.com/ampas/aces-dev
// (as seen on pixls.us)
void ImProcFunctions::log_encode(Imagefloat *rgb, struct local_params & lp, bool multiThread, int bfw, int bfh)
{
   // BENCHFUN
    const float gray = 0.01f * lp.sourcegray;
    const float shadows_range = lp.blackev;

    float dynamic_range = max(lp.whiteev - lp.blackev, 0.5f);
    const float noise = pow_F(2.f, -16.f);
    const float log2 = xlogf(2.f);
    const float base = lp.targetgray > 1 && lp.targetgray < 100 && dynamic_range > 0 ? find_gray(std::abs(lp.blackev) / dynamic_range, 0.01f * lp.targetgray) : 0.f;
    const float linbase = rtengine::max(base, 2.f);//2 to avoid bad behavior
    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
    if (settings->verbose) {
        printf("Base Log encoding std=%5.1f\n", (double) linbase);
    }

    const auto apply =
    [ = ](float x, bool scale = true) -> float {
        if (scale)
        {
            x /= 65535.f;
        }

        x = rtengine::max(x, noise);
        x = rtengine::max(x / gray, noise);
        x = rtengine::max((xlogf(x) / log2 - shadows_range) / dynamic_range, noise);
        assert(x == x);

        if (linbase > 0.f)
        {
            x = xlog2lin(x, linbase);
        }

        if (scale)
        {
            return x * 65535.f;
        } else {
            return x;
        }
    };

    const float detail = lp.detail;
    const int W = rgb->getWidth(), H = rgb->getHeight();

    if (detail == 0.f) {//no local contrast
#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                float r = rgb->r(y, x);
                float g = rgb->g(y, x);
                float b = rgb->b(y, x);
                float m = norm2(r, g, b, ws);

                if (m > noise) {
                    float mm = apply(m);
                    float f = mm / m;
                    f = min(f, 1000000.f);
                    
                    r *= f;
                    b *= f;
                    g *= f;
                    r = CLIP(r);
                    g = CLIP(g);
                    b = CLIP(b);
                }

                assert(r == r);
                assert(g == g);
                assert(b == b);

                rgb->r(y, x) = r;
                rgb->g(y, x) = g;
                rgb->b(y, x) = b;
            }
        }
    } else  {//local contrast

        array2D<float> Y(W, H);
        {
            constexpr float base_posterization = 20.f;
            array2D<float> Y2(W, H);

#ifdef _OPENMP
            #pragma omp parallel for if (multiThread)
#endif
            for (int y = 0; y < H; ++y) {
                for (int x = 0; x < W; ++x) {
                    Y2[y][x] = norm2(rgb->r(y, x), rgb->g(y, x), rgb->b(y, x), ws) / 65535.f;
                    float l = xlogf(rtengine::max(Y2[y][x], 1e-9f));
                    float ll = round(l * base_posterization) / base_posterization;
                    Y[y][x] = xexpf(ll);
                    assert(std::isfinite(Y[y][x]));
                }
            }
            const float radius = rtengine::max(rtengine::max(bfw, W), rtengine::max(bfh, H)) / 30.f;
            const float epsilon = 0.005f;
            rtengine::guidedFilter(Y2, Y, Y, radius, epsilon, multiThread);
        }
        const float blend = detail;

#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                float &r = rgb->r(y, x);
                float &g = rgb->g(y, x);
                float &b = rgb->b(y, x);
                float t = Y[y][x];
                float t2;

                if (t > noise && (t2 = norm2(r, g, b, ws)) > noise) {
                    float c = apply(t, false);
                    float f = c / t;
                    //   float t2 = norm(r, g, b);
                    float f2 = apply(t2) / t2;
                    f = intp(blend, f, f2);
                    f = min(f, 1000000.f);
                   
               //     assert(std::isfinite(f));
                    r *= f;
                    g *= f;
                    b *= f;
                    r = CLIP(r);
                    g = CLIP(g);
                    b = CLIP(b);
               //     assert(std::isfinite(r));
               //     assert(std::isfinite(g));
               //     assert(std::isfinite(b));
                }
            }
        }

    }
}

void ImProcFunctions::getAutoLogloc(int sp, ImageSource *imgsrc, float *sourceg, float *blackev, float *whiteev, bool *Autogr, float *sourceab,  int fw, int fh, float xsta, float xend, float ysta, float yend, int SCALE)
{
    //BENCHFUN
//adpatation to local adjustments Jacques Desmis 12 2019 and 11 2021 (from ART)
    const PreviewProps pp(0, 0, fw, fh, SCALE);

    Imagefloat img(int(fw / SCALE + 0.5), int(fh / SCALE + 0.5));
    const ProcParams neutral;

    imgsrc->getImage(imgsrc->getWB(), TR_NONE, &img, pp, params->toneCurve, neutral.raw);
    imgsrc->convertColorSpace(&img, params->icm, imgsrc->getWB());
    float minVal = RT_INFINITY;
    float maxVal = -RT_INFINITY;
    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);

    constexpr float noise = 1e-5;
    const int h = fh / SCALE;
    const int w = fw / SCALE;

    const int hsta = ysta * h;
    const int hend = yend * h;

    const int wsta = xsta * w;
    const int wend = xend * w;
    int www = int(fw / SCALE + 0.5);
    int hhh = int(fh / SCALE + 0.5);
    array2D<float> YY(www, hhh);

    double mean = 0.0;
    int nc = 0;
    for (int y = hsta; y < hend; ++y) {
        for (int x = wsta; x < wend; ++x) {
            const float r = img.r(y, x), g = img.g(y, x), b = img.b(y, x);
            YY[y][x] = norm2(r, g, b, ws) / 65535.f;//norm2 to find a best color luminance response in RGB 
            mean += static_cast<double>((float) ws[1][0] * Color::gamma_srgb(r) + (float) ws[1][1] * Color::gamma_srgb(g) + (float) ws[1][2] * Color::gamma_srgb(b));
            //alternative to fing gray in case of above process does not works
            nc++;
        }
    }

    for (int y = hsta; y < hend; ++y) {
        for (int x = wsta; x < wend; ++x) {
            float l = YY[y][x];
            if (l > noise) {
                minVal = min(minVal, l);
                maxVal = max(maxVal, l);
            }
        }
    }
    
    maxVal *= 1.45f; //(or 1.5f...) slightly increase max to take into account illuminance incident light
    minVal *= 0.55f; //(or 0.5f...) slightly decrease min to take into account illuminance incident light
    //E = 2.5*2^EV => e=2.5 depends on the sensor type C=250 e=2.5 to C=330 e=3.3
    //repartition with 2.5 between 1.45 Light and shadows 0.58 => a little more 0.55...
    // https://www.pixelsham.com/2020/12/26/exposure-value-measurements/
    // https://en.wikipedia.org/wiki/Light_meter
    if (maxVal > minVal) {
        const float log2 = std::log(2.f);
        const float dynamic_range = -xlogf(minVal / maxVal) / log2;

        if (settings->verbose) {
            std::cout << "AutoLog: min = " << minVal << ", max = " << maxVal
                      << ", Dynamic Range = " << dynamic_range << std::endl;
        }

        if (Autogr[sp]) {
            double tot = 0.0;
            int n = 0;
            //0.05 0.25 arbitrary values around gray point 0.18 to find a good value as "gray" for "gain"
            const float gmax = rtengine::min(maxVal / 2.f, 0.25f);
            const float gmin = rtengine::max(minVal * std::pow(2.f, rtengine::max((dynamic_range - 1.f) / 2.f, 1.f)), 0.05f);

            if (settings->verbose) {
                std::cout << "         gray boundaries: " << gmin << ", " << gmax << std::endl;
            }

            for (int y = hsta; y < hend; ++y) {
                for (int x = wsta; x < wend; ++x) {
                    const float l = img.g(y, x) / 65535.f;

                    if (l >= gmin && l <= gmax) {
                        tot += static_cast<double>(l);
                        ++n;
                    }
                }
            }

            if (n > 0) {
                sourceg[sp] = tot / n * 100.0;

                if (settings->verbose) {
                    std::cout << "         computed gray point from " << n << " samples: " << sourceg[sp] << std::endl;
                }
            } else {//I change slightly this part of algo - more progressivity...best response in very low exposure images
                mean /= (nc * 65535.0);
                float yb;
                yb = 1.5f + 100.f * pow_F(mean, 1.8f);//empirical formula for Jz and log encode for low exposure images

                sourceg[sp] = yb;
                if (settings->verbose) {
                    std::cout << "         no samples found in range, resorting to Yb gray point value " << sourceg[sp]  << std::endl;
                }
            }
        }
        
        constexpr float MIN_WHITE = 2.f;
        constexpr float MAX_BLACK = -3.5f;

        const float gray = sourceg[sp] / 100.f;
        whiteev[sp] = rtengine::max(xlogf(maxVal / gray) / log2, MIN_WHITE);
        blackev[sp] = rtengine::min(whiteev[sp] - dynamic_range, MAX_BLACK);


        //calculate La - Absolute luminance shooting

        const FramesMetaData* metaData = imgsrc->getMetaData();
        
        float fnum = metaData->getFNumber();          // F number
        float fiso = metaData->getISOSpeed() ;        // ISO
        float fspeed = metaData->getShutterSpeed() ;  // Speed
        double fcomp = metaData->getExpComp();        // Compensation +/-
        double adap;

        if (fnum < 0.3f || fiso < 5.f || fspeed < 0.00001f) { //if no exif data or wrong
            adap = 2000.;
        } else {
            double E_V = fcomp + std::log2(double ((fnum * fnum) / fspeed / (fiso / 100.f)));
            double kexp = 0.;
            E_V += kexp * params->toneCurve.expcomp;// exposure compensation in tonecurve ==> direct EV
            E_V += 0.5 * std::log2(params->raw.expos);  // exposure raw white point ; log2 ==> linear to EV
            adap = pow(2.0, E_V - 3.0);  // cd / m2  ==> 3.0 = log2(8) =>fnum*fnum/speed = Luminance (average scene) * fiso / K (K is the reflected-light meter calibration constant according to the sensors about 12.5 or 14 
            // end calculation adaptation scene luminosity
        }
        
        sourceab[sp] = adap;

    }
}

void tone_eq(ImProcFunctions *ipf, Imagefloat *rgb, const struct local_params &lp, const Glib::ustring &workingProfile, double scale, bool multithread)
{
    ToneEqualizerParams params;
    params.enabled = true;
    params.regularization = lp.detailsh;
    params.pivot = lp.tePivot;
    std::copy(lp.mullocsh, lp.mullocsh + params.bands.size(), params.bands.begin());
    ipf->toneEqualizer(rgb, params, workingProfile, scale, multithread);
}
void ImProcFunctions::loccont(int bfw, int bfh, LabImage* tmp1, float rad, float stren, int sk)
{
    if (rad > 0.f) {
        array2D<float> guide(bfw, bfh);
        array2D<float> LL(bfw, bfh);
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif
        for (int y = 0; y < bfh ; y++) {
            for (int x = 0; x < bfw; x++) {
                LL[y][x] = tmp1->L[y][x];
                float ll = LL[y][x] / 32768.f;
                guide[y][x] = xlin2log(rtengine::max(ll, 0.f), 10.f);
            }
        }
        array2D<float> iL(bfw, bfh, LL, 0);
        float gu = stren * rad;
        int r = rtengine::max(int(gu / sk), 1);
        const double epsil = 0.001 * std::pow(2.f, -10);
        float st = 0.01f * rad;
        rtengine::guidedFilterLog(guide, 10.f, LL, r, epsil, false);
        
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif
        for (int y = 0; y < bfh ; y++) {
            for (int x = 0; x < bfw; x++) {
                LL[y][x] = intp(st, LL[y][x] , iL[y][x]);
                tmp1->L[y][x] = LL[y][x];
            }
        }
    }
}

void sigmoidla (float &valj, float thresj, float lambda) 
{
    //thres : shifts the action of sigmoid to darker tones or lights
    //lambda : changes the "slope" of the sigmoid. Low values give a flat curve, high values a "rectangular / orthogonal" curve
     valj =  1.f / (1.f + xexpf(lambda - (lambda / thresj) * valj));
}


void gamutjz (double &Jz, double &az, double &bz, double pl, const double wip[3][3], const float higherCoef, const float lowerCoef)
{//Not used...bad results
        constexpr float ClipLevel = 65535.0f;
        bool inGamut;
      //  int nb = 0;
        do {
            inGamut = true;
            double L_, M_, S_;
            double xx, yy, zz;
            bool zcam = false;
            Ciecam02::jzczhzxyz (xx, yy, zz, Jz, az, bz, pl, L_, M_, S_, zcam);
            double x, y, z;
            x = 65535. * (d65_d50[0][0] * xx + d65_d50[0][1] * yy + d65_d50[0][2] * zz);
            y = 65535. * (d65_d50[1][0] * xx + d65_d50[1][1] * yy + d65_d50[1][2] * zz);
            z = 65535. * (d65_d50[2][0] * xx + d65_d50[2][1] * yy + d65_d50[2][2] * zz);
            float R,G,B;
            Color:: xyz2rgb(x, y, z, R, G, B, wip);
            if (rtengine::min(R, G, B) < 0.f  || rtengine::max(R, G, B) > ClipLevel) {
            //    nb++;
                double hz = xatan2f(bz, az);
                float2 sincosval = xsincosf(hz);
                double Cz = sqrt(az * az + bz * bz);
               // printf("cz=%f jz=%f" , (double) Cz, (double) Jz);
                Cz *= (double) higherCoef;
                if(Cz < 0.01 && Jz > 0.05) {//empirical values
                    Jz -= (double) lowerCoef;
                }
                az = clipazbz(Cz * (double) sincosval.y);
                bz = clipazbz(Cz * (double) sincosval.x);
                
                inGamut = false;
            }
        } while (!inGamut);
}

void ImProcFunctions::ciecamloc_02float(const struct local_params& lp, int sp, LabImage* lab, int bfw, int bfh, int call, int sk, const LUTf& cielocalcurve, bool localcieutili, const LUTf& cielocalcurve2, bool localcieutili2, const LUTf& jzlocalcurve, bool localjzutili, const LUTf& czlocalcurve, bool localczutili, const LUTf& czjzlocalcurve, bool localczjzutili, const LocCHCurve& locchCurvejz, const LocHHCurve& lochhCurvejz, const LocLHCurve& loclhCurvejz, bool HHcurvejz, bool CHcurvejz, bool LHcurvejz, const LocwavCurve& locwavCurvejz, bool locwavutilijz
)
{
//    BENCHFUN
//possibility to reenable Zcam
    if(!params->locallab.spots.at(sp).activ) {//disable all ciecam functions
        return;
    }
    bool ciec = false;
    bool iscie = false;
    if (params->locallab.spots.at(sp).ciecam && params->locallab.spots.at(sp).explog && call == 1) {
        ciec = true;
        iscie = false;
    }
    else if (params->locallab.spots.at(sp).expcie && call == 0) {
        ciec = true;
        iscie = true;
    }
    bool z_cam = false; //params->locallab.spots.at(sp).jabcie; //alaways use normal algorithm, Zcam giev often bad results
    bool jabcie = false;//always disabled
    bool islogjz = params->locallab.spots.at(sp).forcebw;
    bool issigjz = params->locallab.spots.at(sp).sigjz;
    bool issigq = params->locallab.spots.at(sp).sigq;
    bool islogq = params->locallab.spots.at(sp).logcie;

    //sigmoid J Q variables
    const float sigmoidlambda = params->locallab.spots.at(sp).sigmoidldacie; 
    const float sigmoidth = params->locallab.spots.at(sp).sigmoidthcie; 
    const float sigmoidbl = params->locallab.spots.at(sp).sigmoidblcie; 
    const bool sigmoidqj = params->locallab.spots.at(sp).sigmoidqjcie; 

    TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix(params->icm.workingProfile);
    const double wip[3][3] = {//improve precision with double
        {wiprof[0][0], wiprof[0][1], wiprof[0][2]},
        {wiprof[1][0], wiprof[1][1], wiprof[1][2]},
        {wiprof[2][0], wiprof[2][1], wiprof[2][2]}
    };
    float plum = (float) params->locallab.spots.at(sp).pqremapcam16;

    int mocam = 1;
    if(params->locallab.spots.at(sp).modecam == "all") {
        mocam = 10;//à remettre à 0 si modecam = "all"
    } else if(params->locallab.spots.at(sp).modecam == "cam16") {
        mocam = 1;
    } else if(params->locallab.spots.at(sp).modecam == "jz") {
        mocam = 2;
//    } else if(params->locallab.spots.at(sp).modecam == "zcam") {
//        mocam = 3;
    }

    int mecamcurve = 0;
    if(params->locallab.spots.at(sp).toneMethodcie == "one") {
        mecamcurve = 0;
    } else if(params->locallab.spots.at(sp).toneMethodcie == "two") {
        mecamcurve = 1;
    }

    int mecamcurve2 = 0;
    if(params->locallab.spots.at(sp).toneMethodcie2 == "onec") {
        mecamcurve2 = 0;
    } else if(params->locallab.spots.at(sp).toneMethodcie2 == "twoc") {
        mecamcurve2 = 1;
    } else if(params->locallab.spots.at(sp).toneMethodcie2 == "thrc") {
        mecamcurve2 = 2;
    }

    float th = 1.f;
    const float at = 1.f - sigmoidth;
    const float bt = sigmoidth;

    const float ath = sigmoidth - 1.f;
    const float bth = 1;
    float sila = pow_F(sigmoidlambda, 0.5f);
    const float sigm = 3.3f + 7.1f *(1.f - sila);//e^10.4 = 32860 => sigm vary from 3.3 to 10.4
    const float bl = sigmoidbl;
    //end sigmoid

    int width = lab->W, height = lab->H;
    float Yw;
    Yw = 1.0f;
    double Xw, Zw;
    float f = 0.f, nc = 0.f, la, c = 0.f, xw, yw, zw, f2 = 1.f, c2 = 1.f, nc2 = 1.f, yb2;
    float fl, n, nbb, ncb, aw; //d
    float xwd, ywd, zwd, xws, yws, zws;
    //  int alg = 0;
    double Xwout, Zwout;
    double Xwsc, Zwsc;

    LUTu hist16J(32768, LUT_CLIP_BELOW | LUT_CLIP_ABOVE, true);
    LUTu hist16Q(32768, LUT_CLIP_BELOW | LUT_CLIP_ABOVE, true);
    //for J light and contrast
    LUTf CAMBrightCurveJ(32768, LUT_CLIP_BELOW | LUT_CLIP_ABOVE);
    LUTf CAMBrightCurveQ(32768, LUT_CLIP_BELOW | LUT_CLIP_ABOVE);


#ifdef _OPENMP
    const int numThreads = min(max(width * height / 65536, 1), omp_get_max_threads());
    #pragma omp parallel num_threads(numThreads) if(numThreads>1)
#endif
    {
    LUTu hist16Jthr(hist16J.getSize(), LUT_CLIP_BELOW | LUT_CLIP_ABOVE, true);
    LUTu hist16Qthr(hist16Q.getSize(), LUT_CLIP_BELOW | LUT_CLIP_ABOVE, true);

#ifdef _OPENMP
    #pragma omp for
#endif

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) { //rough correspondence between L and J
            float currL = lab->L[i][j] / 327.68f;
            float koef; //rough correspondence between L and J

            if (currL > 50.f) {
                if (currL > 70.f) {
                    if (currL > 80.f) {
                        if (currL > 85.f) {
                            koef = 0.97f;
                        } else {
                            koef = 0.93f;
                        }
                    } else {
                        koef = 0.87f;
                    }
                } else {
                    if (currL > 60.f) {
                        koef = 0.85f;
                    } else {
                        koef = 0.8f;
                    }
                }
            } else {
                if (currL > 10.f) {
                    if (currL > 20.f) {
                        if (currL > 40.f) {
                            koef = 0.75f;
                        } else {
                            koef = 0.7f;
                        }
                    } else {
                        koef = 0.9f;
                    }
                } else {
                    koef = 1.0;
                }
            }

            hist16Jthr[(int)((koef * lab->L[i][j]))]++;    //evaluate histogram luminance L # J
            hist16Qthr[CLIP((int)(32768.f * sqrt((koef * (lab->L[i][j])) / 32768.f)))]++;     //for brightness Q : approximation for Q=wh*sqrt(J/100)  J not equal L
        }
    }

#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            hist16J += hist16Jthr;
            hist16Q += hist16Qthr;
        }
    }
#ifdef _OPENMP
    static_cast<void>(numThreads); // to silence cppcheck warning
#endif

    //evaluate lightness, contrast

    if (ciec) {
        float contL = 0.f;
        float lightL = 0.f;
        float contQ = 0.f;
        float lightQ = 0.f;
        if(iscie) {
            contL = 0.6 * params->locallab.spots.at(sp).contlcie; //0.6 less effect, no need 1.
            lightL = 0.4 * params->locallab.spots.at(sp).lightlcie; //0.4 less effect, no need 1.
            contQ = 0.5 * params->locallab.spots.at(sp).contqcie; //0.5 less effect, no need 1.
            lightQ = 0.4 * params->locallab.spots.at(sp).lightqcie; //0.4 less effect, no need 1.
        } else {
            contL = 0.6 * params->locallab.spots.at(sp).contl; //0.6 less effect, no need 1.
            lightL = 0.4 * params->locallab.spots.at(sp).lightl; //0.4 less effect, no need 1.
            contQ = 0.5 * params->locallab.spots.at(sp).contq; //0.5 less effect, no need 1.
            lightQ = 0.4 * params->locallab.spots.at(sp).lightq; //0.4 less effect, no need 1.
            
        }
        float contthresL = 0.f;
        
        if(iscie) {
            contthresL = params->locallab.spots.at(sp).contthrescie;
        } else {
            contthresL = params->locallab.spots.at(sp).contthres;
        }
        float contthresQ = contthresL;
        if(contL < 0.f) {
            contthresL *= -1;
        } 
        float thL = 0.6f;
        thL = 0.3f * contthresL + 0.6f;
        
        if(contQ < 0.f) {
            contthresQ *= -1;
        } 
        float thQ = 0.6f;
        thQ = 0.3f * contthresQ + 0.6f;
        Ciecam02::curveJfloat(lightL, contL, thL, hist16J, CAMBrightCurveJ); //lightness J and contrast J
        CAMBrightCurveJ /= 327.68f;
        
        Ciecam02::curveJfloat(lightQ, contQ, thQ, hist16Q, CAMBrightCurveQ); //brightness Q and contrast Q
    }
    
    
    int tempo = 5000;
    if(params->locallab.spots.at(sp).expvibrance && call == 2) {
        if (params->locallab.spots.at(sp).warm > 0) {
            tempo = 5000 - 30 * params->locallab.spots.at(sp).warm;
        } else if (params->locallab.spots.at(sp).warm < 0){
            tempo = 5000 - 70 * params->locallab.spots.at(sp).warm;
        }
    }


    if(ciec) {
        if(iscie) {
            if (params->locallab.spots.at(sp).catadcie > 0) {
                tempo = 5000 - 30 * params->locallab.spots.at(sp).catadcie;
            } else if (params->locallab.spots.at(sp).catadcie < 0){
                tempo = 5000 - 70 * params->locallab.spots.at(sp).catadcie;
            }
        } else {
            if (params->locallab.spots.at(sp).catad > 0) {
                tempo = 5000 - 30 * params->locallab.spots.at(sp).catad;
            } else if (params->locallab.spots.at(sp).catad < 0){
                tempo = 5000 - 70 * params->locallab.spots.at(sp).catad;
            }
        }
    }

    ColorTemp::temp2mulxyz(params->wb.temperature, params->wb.method, params->wb.observer, Xw, Zw);  //compute white Xw Yw Zw  : white current WB
    ColorTemp::temp2mulxyz(tempo, "Custom", params->wb.observer, Xwout, Zwout);
    ColorTemp::temp2mulxyz(5000, "Custom", params->wb.observer, Xwsc, Zwsc);

    //viewing condition for surrsrc
    f  = 1.00f;
    c  = 0.69f;
    nc = 1.00f;
    //viewing condition for surround
    f2 = 1.0f, c2 = 0.69f, nc2 = 1.0f;
    if(ciec) {
        if(iscie) { 
        //surround source with only 2 choices (because Log encoding before)
            if (params->locallab.spots.at(sp).sursourcie == "Average") {
                f = 1.0f, c = 0.69f, nc = 1.0f;
            } else if (params->locallab.spots.at(sp).sursourcie == "Dim") {
                f  = 0.9f;
                c  = 0.59f;
                nc = 0.9f;
            } else if (params->locallab.spots.at(sp).sursourcie == "Dark") {
                f  = 0.8f;
                c  = 0.525f;
                nc = 0.8f;
            }
        } else {
            if (params->locallab.spots.at(sp).sursour == "Average") {
                f = 1.0f, c = 0.69f, nc = 1.0f;
            } else if (params->locallab.spots.at(sp).sursour == "Dim") {
                f  = 0.9f;
                c  = 0.59f;
                nc = 0.9f;
            } else if (params->locallab.spots.at(sp).sursour == "Dark") {
                f  = 0.8f;
                c  = 0.525f;
                nc = 0.8f;
            }
        }

        //viewing condition for surround
        if(iscie) {
            if (params->locallab.spots.at(sp).surroundcie == "Average") {
                f2 = 1.0f, c2 = 0.69f, nc2 = 1.0f;
            } else if (params->locallab.spots.at(sp).surroundcie == "Dim") {
                f2  = 0.9f;
                c2  = 0.59f;
                nc2 = 0.9f;
            } else if (params->locallab.spots.at(sp).surroundcie == "Dark") {
                f2  = 0.8f;
                c2  = 0.525f;
                nc2 = 0.8f;
            } else if (params->locallab.spots.at(sp).surroundcie == "ExtremelyDark") {
                f2  = 0.8f;
                c2  = 0.41f;
                nc2 = 0.8f;
            }
        } else {
            if (params->locallab.spots.at(sp).surround == "Average") {
                f2 = 1.0f, c2 = 0.69f, nc2 = 1.0f;
            } else if (params->locallab.spots.at(sp).surround == "Dim") {
                f2  = 0.9f;
                c2  = 0.59f;
                nc2 = 0.9f;
            } else if (params->locallab.spots.at(sp).surround == "Dark") {
                f2  = 0.8f;
                c2  = 0.525f;
                nc2 = 0.8f;
            } else if (params->locallab.spots.at(sp).surround == "ExtremelyDark") {
                f2  = 0.8f;
                c2  = 0.41f;
                nc2 = 0.8f;
            }
            
        }
    }

    xwd = 100.0 * Xwout;
    zwd = 100.0 * Zwout;
    ywd = 100.f;

    xws = 100.0 * Xwsc;
    zws = 100.0 * Zwsc;
    yws = 100.f;


    //La and la2 = ambiant luminosity scene and viewing
    la = 400.f;
    float la2 = 400.f;
    if(ciec) {
        if(iscie) {
            la = params->locallab.spots.at(sp).sourceabscie;
            la2 = params->locallab.spots.at(sp).targabscie;
        } else {
            la = params->locallab.spots.at(sp).sourceabs;
            la2 = params->locallab.spots.at(sp).targabs;
        }
    }

    const float pilot = 2.f;
    const float pilotout = 2.f;
    double avgm = 0.;
    //algoritm's params
    float yb = 18.f;
    yb2 = 18;
    if(ciec) {
        if(iscie) {
            yb = params->locallab.spots.at(sp).sourceGraycie;//
            avgm = (double) pow_F(0.01f * (yb - 1.f), 0.45f);;
            yb2 = params->locallab.spots.at(sp).targetGraycie;
        } else {
            yb = params->locallab.spots.at(sp).targetGray;//target because we are after Log encoding
            yb2 = params->locallab.spots.at(sp).targetGray;
        }
    }
    if(params->locallab.spots.at(sp).expcie && call == 10 && params->locallab.spots.at(sp).modecam == "jz") {
            yb = params->locallab.spots.at(sp).sourceGraycie;//for Jz calculate Yb and surround in Lab and cam16 before process Jz
            la = params->locallab.spots.at(sp).sourceabscie;

            if (params->locallab.spots.at(sp).sursourcie == "Average") {
                f = 1.0f, c = 0.69f, nc = 1.0f;
            } else if (params->locallab.spots.at(sp).sursourcie == "Dim") {
                f  = 0.9f;
                c  = 0.59f;
                nc = 0.9f;
            } else if (params->locallab.spots.at(sp).sursourcie == "Dark") {
                f  = 0.8f;
                c  = 0.525f;
                nc = 0.8f;
            }
    }
    
    float schr = 0.f;
    float mchr = 0.f;
    float cchr = 0.f;
    float rstprotection = 0.f;
    float hue = 0.f;
/*
    float mchrz = 0.f;
    float schrz = 0.f;
    float cchrz = 0.f;
*/
    if (ciec) {
        if(iscie) {
            rstprotection =  params->locallab.spots.at(sp).rstprotectcie;
            hue = params->locallab.spots.at(sp).huecie;

            cchr = params->locallab.spots.at(sp).chromlcie;
            if (cchr == -100.0f) {
                    cchr = -99.8f;
            }

            schr = params->locallab.spots.at(sp).saturlcie;

            if (schr > 0.f) {
            schr = schr / 2.f;    //divide sensibility for saturation
            }

            if (schr == -100.f) {
                schr = -99.8f;
            }

            mchr = params->locallab.spots.at(sp).colorflcie;

            if (mchr == -100.0f) {
                mchr = -99.8f ;
            }
            if (mchr == 100.0f) {
                mchr = 99.9f;
            }
/*
            mchrz = 0.5f * (float) params->locallab.spots.at(sp).colorflzcam;
            schrz = 0.5f * (float) params->locallab.spots.at(sp).saturzcam;
            cchrz = 0.5f * (float) params->locallab.spots.at(sp).chromzcam;
*/
        } else {
            cchr = params->locallab.spots.at(sp).chroml;
            if (cchr == -100.0f) {
                    cchr = -99.8f;
            }

            schr = params->locallab.spots.at(sp).saturl;

            if (schr > 0.f) {
            schr = schr / 2.f;    //divide sensibility for saturation
            }

            if (schr == -100.f) {
                schr = -99.8f;
            }

            mchr = params->locallab.spots.at(sp).colorfl;

            if (mchr == -100.0f) {
                mchr = -99.8f ;
            }
            if (mchr == 100.0f) {
                mchr = 99.9f;
            }
        }
    }

    float d, dj;

    // const int gamu = 0; //(params->colorappearance.gamut) ? 1 : 0;
    xw = 100.0 * Xw;
    yw = 100.f * Yw;
    zw = 100.0 * Zw;
    float xw1 = xws, yw1 = yws, zw1 = zws, xw2 = xwd, yw2 = ywd, zw2 = zwd;
    float cz, wh, pfl;
    int c16 = 16;//always cat16
    bool c20 = true;
    if(c20  && plum > 100.f) {
        c16 = 21;//I define 21...for 2021 :)
    }
    int level_bljz = params->locallab.spots.at(sp).csthresholdjz.getBottomLeft();
    int level_hljz = params->locallab.spots.at(sp).csthresholdjz.getTopLeft();
    int level_brjz = params->locallab.spots.at(sp).csthresholdjz.getBottomRight();
    int level_hrjz = params->locallab.spots.at(sp).csthresholdjz.getTopRight();

    float alowjz = 1.f;
    float blowjz = 0.f;

    if (level_hljz != level_bljz) {
        alowjz = 1.f / (level_hljz - level_bljz);
        blowjz = -alowjz * level_bljz;
    }

    float ahighjz = 1.f;
    float bhighjz = 0.f;

    if (level_hrjz != level_brjz) {
        ahighjz = 1.f / (level_hrjz - level_brjz);
        bhighjz =  -ahighjz * level_brjz;
    }
    float sigmalcjz = params->locallab.spots.at(sp).sigmalcjz;
    float jzamountchr = 0.01 * params->locallab.spots.at(sp).thrhjzcie;
    bool jzch = params->locallab.spots.at(sp).chjzcie;
    double jzamountchroma = 0.01 * settings->amchromajz;
    if(jzamountchroma < 0.05) {
        jzamountchroma = 0.05;
    }
    if(jzamountchroma > 2.) {
        jzamountchroma = 2.;
    }
    
    Ciecam02::initcam1float(yb, pilot, f, la, xw, yw, zw, n, d, nbb, ncb, cz, aw, wh, pfl, fl, c, c16, plum);
    const float pow1 = pow_F(1.64f - pow_F(0.29f, n), 0.73f);
    float nj, nbbj, ncbj, czj, awj, flj;
    Ciecam02::initcam2float(yb2, pilotout, f2,  la2,  xw2,  yw2,  zw2, nj, dj, nbbj, ncbj, czj, awj, flj, c16, plum);
#ifdef __SSE2__
    const float reccmcz = 1.f / (c2 * czj);
#endif
    const float epsil = 0.0001f;
    const float coefQ = 32767.f / wh;
    const float coefq = 1 / wh;
    const float pow1n = pow_F(1.64f - pow_F(0.29f, nj), 0.73f);
    const float coe = pow_F(fl, 0.25f);
    const float QproFactor = (0.4f / c) * (aw + 4.0f) ;
    const double shadows_range =  params->locallab.spots.at(sp).blackEvjz;
    const double targetgray = params->locallab.spots.at(sp).targetjz;
    double targetgraycor = 0.15;
    double dynamic_range = std::max(params->locallab.spots.at(sp).whiteEvjz - shadows_range, 0.5);
    const double noise = pow(2., -16.6);//16.6 instead of 16 a little less than others, but we work in double
    const double log2 = xlog(2.);
    const float log2f = xlogf(2.f);

    if((mocam == 0 || mocam ==2)  && call == 0) {//Jz az bz ==> Jz Cz Hz before Ciecam16
        double mini = 1000.;
        double maxi = -1000.;
        double sum = 0.;
        int nc = 0;
        double epsiljz = 0.0001;
        //Remapping see https://hal.inria.fr/hal-02131890/document    I took some ideas in this text, and add my personal adaptation
        // image quality assessment of HDR and WCG images https://tel.archives-ouvertes.fr/tel-02378332/document
        double adapjz = params->locallab.spots.at(sp).adapjzcie;
        double jz100 = params->locallab.spots.at(sp).jz100;
        double pl = params->locallab.spots.at(sp).pqremap;
        double jzw, azw, bzw;
        jzw = 0.18;//Jz white

        bool Qtoj = params->locallab.spots.at(sp).qtoj;//betwwen lightness to brightness
        const bool logjz =  params->locallab.spots.at(sp).logjz;//log encoding

//calculate min, max, mean for Jz
#ifdef _OPENMP
            #pragma omp parallel for reduction(min:mini) reduction(max:maxi) reduction(+:sum) if(multiThread)
#endif
        for (int i = 0; i < height; i+=1) {
            for (int k = 0; k < width; k+=1) {
                float L = lab->L[i][k];
                float a = lab->a[i][k];
                float b = lab->b[i][k];
                float x, y, z;
                //convert Lab => XYZ
                Color::Lab2XYZ(L, a, b, x, y, z);
                x = x / 65535.f;
                y = y / 65535.f;
                z = z / 65535.f;
                double Jz, az, bz;
                double xx, yy, zz;
                //D50 ==> D65
                xx = (d50_d65[0][0] * (double) x + d50_d65[0][1] * (double) y + d50_d65[0][2] * (double) z);
                yy = (d50_d65[1][0] * (double) x + d50_d65[1][1] * (double) y + d50_d65[1][2] * (double) z);
                zz = (d50_d65[2][0] * (double) x + d50_d65[2][1] * (double) y + d50_d65[2][2] * (double) z);

                double L_p, M_p, S_p;
                bool zcam = z_cam;

                Ciecam02::xyz2jzczhz (Jz, az, bz, xx, yy, zz, pl, L_p, M_p, S_p, zcam);
                if(Jz > maxi) {
                    maxi = Jz;
                }
                if(Jz < mini) {
                    mini = Jz;
                }
                sum += Jz;
            // I read bz, az values and Hz ==> with low chroma values Hz are very different from lab always around 1.4 radians ???? for blue...
            }
        }
        nc = height * width;
        sum = sum / nc;
        maxi += epsiljz;
        sum += epsiljz;
        //remapping Jz 
        double ijz100 = 1./jz100;
        double ajz = (ijz100 - 1.)/9.;//9 = sqrt(100) - 1 with a parabolic curve after jz100 - we can change for others curve ..log...(you must change also in locallabtool2)
        double bjz = 1. - ajz;
        //relation between adapjz and Absolute luminance source (La), adapjz =sqrt(La) - see locallabtool2 adapjzcie 
        double interm = jz100 * (adapjz * ajz + bjz);
        double bj = (10. - maxi) / 9.;
        double aj = maxi -bj;
        double to_screen = (aj * interm + bj) / maxi;
        //to screen - remapping of Jz in function real scene absolute luminance

//        if (settings->verbose) { 
//            printf("ajz=%f bjz=%f adapjz=%f jz100=%f interm=%f to-scrp=%f to_screen=%f\n", ajz, bjz, adapjz, jz100, interm ,to_screenp, to_screen);
//        }
        double to_one = 1.;//only for calculation in range 0..1 or 0..32768
        to_one = 1 / (maxi * to_screen);
        if(adapjz == 10.) {//force original algorithm if La > 10000
            to_screen = 1.;
        }
        if(Qtoj) {
                double xxw = (d50_d65[0][0] * (double) Xw + d50_d65[0][1] * (double) Yw + d50_d65[0][2] * (double) Zw);
                double yyw = (d50_d65[1][0] * (double) Xw + d50_d65[1][1] * (double) Yw + d50_d65[1][2] * (double) Zw);
                double zzw = (d50_d65[2][0] * (double) Xw + d50_d65[2][1] * (double) Yw + d50_d65[2][2] * (double) Zw);
                double L_pa, M_pa, S_pa;
                Ciecam02::xyz2jzczhz (jzw, azw, bzw, xxw, yyw, zzw, pl, L_pa, M_pa, S_pa, z_cam);
                if (settings->verbose) { //calculate Jz white for use of lightness instead brightness
                    printf("Jzwhite=%f \n", jzw);
                }

        }
        const std::unique_ptr<LabImage> temp(new LabImage(width, height));
        const std::unique_ptr<LabImage> tempresid(new LabImage(width, height));
        const std::unique_ptr<LabImage> tempres(new LabImage(width, height));
        array2D<double> JJz(width, height);
        array2D<double> Aaz(width, height);
        array2D<double> Bbz(width, height);
        int highhs =  params->locallab.spots.at(sp).hljzcie;
        int hltonahs = params->locallab.spots.at(sp).hlthjzcie;
        int shadhs = params->locallab.spots.at(sp).shjzcie;
        int shtonals = params->locallab.spots.at(sp).shthjzcie;
        int radhs = params->locallab.spots.at(sp).radjzcie;
        float softjz = (float) params->locallab.spots.at(sp).softjzcie;
        
        avgm = 0.5 * (sum * to_screen * to_one + avgm);//empirical formula
        double miny = 0.1;
        double delta = 0.015 * (double) sqrt(std::max(100.f, la) / 100.f);//small adaptation in function La scene
        double maxy = 0.65;//empirical value
        double maxreal = maxi*to_screen;
        double maxjzw = jzw*to_screen;
        if (settings->verbose) { 
            printf("La=%4.1f PU_adap=%2.1f maxi=%f mini=%f mean=%f, avgm=%f to_screen=%f Max_real=%f to_one=%f\n", (double) la, adapjz, maxi, mini, sum, avgm, to_screen, maxreal, to_one);
        }

        const float sigmoidlambdajz = params->locallab.spots.at(sp).sigmoidldajzcie; 
        const float sigmoidthjz = params->locallab.spots.at(sp).sigmoidthjzcie; 
        const float sigmoidbljz = params->locallab.spots.at(sp).sigmoidbljzcie; 

        float thjz = 1.f;
        const float atjz = 1.f - sigmoidthjz;
        const float btjz = sigmoidthjz;

        const float athjz = sigmoidthjz - 1.f;
        const float bthjz = 1.f;
        float powsig = pow_F(sigmoidlambdajz, 0.5f);
        const float sigmjz = 3.3f + 7.1f *(1.f - powsig);// e^10.4 = 32860
        const float bljz = sigmoidbljz;
        
        double contreal = 0.2 *  params->locallab.spots.at(sp).contjzcie;
        DiagonalCurve jz_contrast({
            DCT_NURBS,
            0, 0,
            avgm - avgm * (0.6 - contreal / 250.0), avgm - avgm * (0.6 + contreal / 250.0),
            avgm + (1. - avgm) * (0.6 - contreal / 250.0), avgm + (1. - avgm) * (0.6 + contreal / 250.0),
            1, 1
        });
        //all calculations in double for best results...but slow
        double lightreal = 0.2 *  params->locallab.spots.at(sp).lightjzcie;
        double chromz =  params->locallab.spots.at(sp).chromjzcie;
        double saturz =  params->locallab.spots.at(sp).saturjzcie;
        double dhue = 0.0174 * params->locallab.spots.at(sp).huejzcie;
        DiagonalCurve jz_light({
            DCT_NURBS,
            0, 0,
            miny, miny + lightreal / 150.,
            maxy, min (1.0, maxy + delta + lightreal / 300.0),
            1, 1
        });
        DiagonalCurve jz_lightn({
            DCT_NURBS,
            0, 0,
            max(0.0, miny  - lightreal / 150.), miny ,
            maxy + delta - lightreal / 300.0, maxy + delta,
            1, 1
        });
        bool wavcurvejz = false;
        if (locwavCurvejz && locwavutilijz) {
            for (int i = 0; i < 500; i++) {
                if (locwavCurvejz[i] != 0.5f) {
                    wavcurvejz = true;
                    break;
                }
            }
        }
        float mjjz = lp.mLjz;
        if(wavcurvejz && lp.mLjz == 0.f) {
            mjjz = 0.0f;//to enable clarity if need in some cases mjjz = 0.0001f
        }

    //log encoding Jz
    double gray = 0.15;
    /*
    const double shadows_range =  params->locallab.spots.at(sp).blackEvjz;
    const double targetgray = params->locallab.spots.at(sp).targetjz;
    double targetgraycor = 0.15;
    double dynamic_range = std::max(params->locallab.spots.at(sp).whiteEvjz - shadows_range, 0.5);
    const double noise = pow(2., -16.6);//16.6 instead of 16 a little less than others, but we work in double
    const double log2 = xlog(2.);
    */
    double base = 10.;
    double linbase = 10.;
    if(logjz) {//with brightness Jz
        gray = 0.01 * params->locallab.spots.at(sp).sourceGraycie;//acts as amplifier (gain) : needs same type of modifications than targetgraycor with pow
        gray = pow(gray, 1.2);//or 1.15 => modification to increase sensitivity gain, only on defaults, of course we can change this value manually...take into account suuround and Yb Cam16
        targetgraycor = pow(0.01 * targetgray, 1.15);//or 1.2 small reduce effect -> take into account a part of surround (before it was at 1.2)
        base = targetgray > 1. && targetgray < 100. && dynamic_range > 0. ? (double) find_gray(std::abs((float) shadows_range) / (float) dynamic_range, (float) (targetgraycor)) : 0.;
        linbase = std::max(base, 2.);//2. minimal base log to avoid very bad results
        if (settings->verbose) {
            printf("Base logarithm encoding Jz=%5.1f\n", linbase);
        }
    }

    const auto applytojz =
    [ = ](double x) -> double {

        x = std::max(x, noise);
        x = std::max(x / gray, noise);//gray = gain - before log conversion
        x = std::max((xlog(x) / log2 - shadows_range) / dynamic_range, noise);//x in range EV
        assert(x == x);

        if (linbase > 0.)//apply log base in function of targetgray blackEvjz and Dynamic Range
        {
            x = xlog2lin(x, linbase);
        }
        return x;
    };

#ifdef _OPENMP
            #pragma omp parallel for if(multiThread)
#endif
        for (int i = 0; i < height; i++) {
            for (int k = 0; k < width; k++) {
                float L = lab->L[i][k];
                float a = lab->a[i][k];
                float b = lab->b[i][k];
                float x, y, z;
                //convert Lab => XYZ
                Color::Lab2XYZ(L, a, b, x, y, z);
                x = x / 65535.f;
                y = y / 65535.f;
                z = z / 65535.f;
                double Jz, az, bz;//double need because matrix with const(1.6295499532821566e-11) and others
                double xx, yy, zz;
                //change WP to D65
                xx = (d50_d65[0][0] * (double) x + d50_d65[0][1] * (double) y + d50_d65[0][2] * (double) z);
                yy = (d50_d65[1][0] * (double) x + d50_d65[1][1] * (double) y + d50_d65[1][2] * (double) z);
                zz = (d50_d65[2][0] * (double) x + d50_d65[2][1] * (double) y + d50_d65[2][2] * (double) z);

                double L_p, M_p, S_p;
                bool zcam = z_cam;
                Ciecam02::xyz2jzczhz (Jz, az, bz, xx, yy, zz, pl, L_p, M_p, S_p, zcam);
                //remapping Jz
                Jz = Jz * to_screen;
                az = az * to_screen;
                bz = bz * to_screen;
                JJz[i][k] = Jz;
                Aaz[i][k] = az;
                Bbz[i][k] = bz;
                if(highhs > 0 || shadhs > 0  || wavcurvejz || mjjz != 0.f || lp.mCjz != 0.f  || LHcurvejz || HHcurvejz || CHcurvejz) {
                    //here we work in float with usual functions  SH / wavelets / curves H
                    temp->L[i][k] = tempresid->L[i][k] = tempres->L[i][k] = (float) to_one * 32768.f * (float) JJz[i][k];
                    temp->a[i][k] = tempresid->a[i][k] = tempres->a[i][k] = (float) to_one * 32768.f * (float) Aaz[i][k];
                    temp->b[i][k] = tempresid->b[i][k] = tempres->b[i][k] = (float) to_one * 32768.f * (float) Bbz[i][k];
                }
            }
        }

        if(highhs > 0 || shadhs > 0) {
            ImProcFunctions::shadowsHighlights(temp.get(), true, 1, highhs, shadhs, radhs, sk, hltonahs * maxi * to_screen * to_one, shtonals * maxi * to_screen * to_one);
#ifdef _OPENMP
            #pragma omp parallel for if(multiThread)
#endif
            for (int i = 0; i < height; i++) {
                for (int k = 0; k < width; k++) {//reinitialize datas after SH...: guide, etc.
                    tempresid->L[i][k] = tempres->L[i][k] = temp->L[i][k];
                    tempresid->a[i][k] = tempres->a[i][k] = temp->a[i][k];
                    tempresid->b[i][k] = tempres->b[i][k] = temp->b[i][k];
                }
            }
        }
        //others "Lab" treatment...to adapt
        
        if(wavcurvejz  || mjjz != 0.f || lp.mCjz != 0.f) {//local contrast wavelet and clarity
#ifdef _OPENMP
            const int numThreads = omp_get_max_threads();
#else
            const int numThreads = 1;

#endif
             // adap maximum level wavelet to size of RT-spot
            int wavelet_level = 1 + params->locallab.spots.at(sp).csthresholdjz.getBottomRight();//retrieve with +1 maximum wavelet_level
            int minwin = rtengine::min(width, height);
            int maxlevelspot = 10;//maximum possible

             // adapt maximum level wavelet to size of crop
            while ((1 << maxlevelspot) >= (minwin * sk) && maxlevelspot  > 1) {
                --maxlevelspot ;
            }


            wavelet_level = rtengine::min(wavelet_level, maxlevelspot);
            int maxlvl = wavelet_level;
            //simple local contrast in function luminance
            if (locwavCurvejz && locwavutilijz && wavcurvejz) {
                float strengthjz = 1.2;
                std::unique_ptr<wavelet_decomposition> wdspot(new wavelet_decomposition(temp->L[0], bfw, bfh, maxlvl, 1, sk, numThreads, lp.daubLen));//lp.daubLen
                if (wdspot->memory_allocation_failed()) {
                    return;
                }
                maxlvl = wdspot->maxlevel();
                wavlc(*wdspot, level_bljz, level_hljz, maxlvl, level_hrjz, level_brjz, ahighjz, bhighjz, alowjz, blowjz, sigmalcjz, strengthjz, locwavCurvejz, numThreads);
                wdspot->reconstruct(temp->L[0], 1.f);

            }
            float thr = 0.001f;
            int flag = 2;
            
            // begin clarity wavelet jz
            if(mjjz != 0.f || lp.mCjz != 0.f) {
                float mL0 = 0.f;
                float mC0 = 0.f;
                bool exec = false;
                float mL = mjjz;
                float mC = lp.mCjz;
                clarimerge(lp, mL, mC, exec, tempresid.get(), wavelet_level, sk, numThreads);

                if (maxlvl <= 4) {
                    mL0 = 0.f;
                    mC0 = 0.f;
                    mL = -1.5f * mL;//increase only for sharpen
                    mC = -mC;
                    thr = 1.f;
                    flag = 0;

                } else {
                    mL0 = mL;
                    mC0 = mC;
                    thr = 1.f;
                    flag = 2;
                }
                LabImage *mergfile = temp.get();
#ifdef _OPENMP
                #pragma omp parallel for if (multiThread)
#endif
                for (int x = 0; x < height; x++)
                    for (int y = 0; y < width; y++) {
                        temp->L[x][y] = clipLoc((1.f + mL0) * mergfile->L[x][y] - mL * tempresid->L[x][y]);
                        temp->a[x][y] = clipC((1.f + mC0) * mergfile->a[x][y] - mC * tempresid->a[x][y]);
                        temp->b[x][y] = clipC((1.f + mC0) * mergfile->b[x][y] - mC * tempresid->b[x][y]);
                }
            }
        
            if (lp.softrjz >= 0.5f && (wavcurvejz || std::fabs(mjjz) > 0.001f)) {//guidedfilter
                softproc(tempres.get(), temp.get(), lp.softrjz, height, width, 0.001, 0.00001, thr, sk, multiThread, flag);
            }
        }

//new curves Hz
#ifdef _OPENMP
                #pragma omp parallel for if (multiThread)
#endif
        for (int i = 0; i < height; i++) {
            for (int k = 0; k < width; k++) {
                float j_z = temp->L[i][k];
                float C_z = sqrt(SQR(temp->a[i][k]) + SQR(temp->b[i][k]));
                float c_z = C_z / 32768.f;
                if (loclhCurvejz && LHcurvejz) {//Jz=f(Hz) curve
                    float kcz = (float) jzamountchr;
                    float Hz = xatan2f (temp->b[i][k], temp->a[i][k]);
                    float l_r = j_z / 32768.f;
                    float kcc = SQR(c_z / kcz);
                    jzch = true;
                    if(jzch == false) {
                        kcc = 1.f;
                    } else if(kcc > 1.f) {
                        kcc = 1.f; //cbrt(kcc);
                    }
                    float valparam = loclhCurvejz[500.f *static_cast<float>(Color::huejz_to_huehsv2((float) Hz))] - 0.5f;

                    float valparamneg;
                    valparamneg = valparam;
                    valparam *= 2.f * kcc;
                    valparamneg *= kcc;
                        if (valparam > 0.f) {
                            l_r = (1.f - valparam) * l_r + valparam * (1.f - SQR(((SQR(1.f - min(l_r, 1.0f))))));
                        } else
                            //for negative
                        {
                            float khue = 1.9f; //in reserve in case of!
                            l_r *= (1.f + khue * valparamneg);
                        }
                    temp->L[i][k] = l_r * 32768.f;
                }
                
                if (locchCurvejz && CHcurvejz) {//Cz=f(Hz) curve
                    float Hz = xatan2f (temp->b[i][k], temp->a[i][k]);
                     const float valparam = 1.5f * (locchCurvejz[500.f * static_cast<float>(Color::huejz_to_huehsv2((float)Hz))] - 0.5f);  //get valp=f(H)
                     float chromaCzfactor = 1.0f + valparam;
                     temp->a[i][k] *= chromaCzfactor;
                     temp->b[i][k] *= chromaCzfactor;
                }
                
                
                 if (lochhCurvejz && HHcurvejz) { // Hz=f(Hz)
                    float Hz = xatan2f (temp->b[i][k], temp->a[i][k]);
                    const float valparam = 1.4f * (lochhCurvejz[500.f * static_cast<float>(Color::huejz_to_huehsv2((float)Hz))] - 0.5f) + static_cast<float>(Hz);
                    Hz = valparam;
                    if ( Hz < 0.0f ) {
                        Hz += (2.f * rtengine::RT_PI_F);
                    }

                    float2 sincosval = xsincosf(Hz);
                    temp->a[i][k] = C_z * sincosval.y;
                    temp->b[i][k] = C_z * sincosval.x;
                }
            }
        }

                if (loclhCurvejz && LHcurvejz && softjz > 0.f) {//Guidedilter for artifacts curve J(H)
                    float thr = 0.00001f;
                    int flag = 2;
                    float softjzr = 0.05f * softjz;
                    softproc(tempres.get(), temp.get(), softjzr, height, width, 0.000001, 0.00000001, thr, sk, multiThread, flag);
                } 


                 if ((lochhCurvejz && HHcurvejz) || (locchCurvejz && CHcurvejz)) { //for artifacts curve H(H)
                    if(softjz > 0.f) {
                        array2D<float> chro(width, height);
                        array2D<float> hue(width, height);
                        array2D<float> guid(width, height);
                        
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < height; y++) {
                            for (int x = 0; x < width; x++) {
                                hue[y][x] = xatan2f(temp->b[y][x], temp->a[y][x]);
                                chro[y][x] = sqrt(SQR(temp->b[y][x]) + SQR(temp->a[y][x]))/32768.f;
                                if ( hue[y][x] < 0.0f ) {
                                    hue[y][x] += (2.f * rtengine::RT_PI_F);
                                }
                                hue[y][x] /= (2.f * rtengine::RT_PI_F);
                                guid[y][x] = tempres->L[y][x] / 32768.f;
                            }
                        }
                        float softr = softjz;
                        const float tmpblur = softr < 0.f ? -1.f / softr : 1.f + softr;
                        const int r2 = rtengine::max<int>(10 / sk * tmpblur + 0.2f, 1);
                        const int r1 = rtengine::max<int>(4 / sk * tmpblur + 0.5f, 1);
                        constexpr float epsilmax = 0.0005f;
                        constexpr float epsilmin = 0.0000001f;
                        constexpr float aepsil = (epsilmax - epsilmin) / 100.f;
                        constexpr float bepsil = epsilmin;
                        const float epsil = softr < 0.f ? 0.001f : aepsil * softr + bepsil;
                        if (lochhCurvejz && HHcurvejz) {
                            rtengine::guidedFilter(guid, hue, hue, r2, 0.5f * epsil, multiThread);
                        }
                        if (locchCurvejz && CHcurvejz) {
                            rtengine::guidedFilter(guid, chro, chro, r1, 0.4f * epsil, multiThread);
                        }

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < height; y++) {
                            for (int x = 0; x < width; x++) {
                                hue[y][x] *= (2.f * rtengine::RT_PI_F);
                                chro[y][x] *= 32768.f;
                                float2 sincosval = xsincosf(hue[y][x]);
                                temp->a[y][x] = chro[y][x] * sincosval.y;
                                temp->b[y][x] = chro[y][x] * sincosval.x;
                            }
                        }
                    }
                 }


///////////////////

        
#ifdef _OPENMP
            #pragma omp parallel  for if(multiThread)
#endif 
        for (int i = 0; i < height; i++) {
            for (int k = 0; k < width; k++) {
                //reconvert to double
                if(highhs > 0 || shadhs > 0  || wavcurvejz || mjjz != 0.f || lp.mCjz != 0.f || LHcurvejz || HHcurvejz || CHcurvejz) {
                    //now we work in double necessary for matrix conversion and when in range 0..1 with use of PQ
                    JJz[i][k] = (double) (temp->L[i][k] / (32768.f * (float) to_one));
                    Aaz[i][k] = (double) (temp->a[i][k] / (32768.f * (float) to_one));
                    Bbz[i][k] = (double) (temp->b[i][k] / (32768.f * (float) to_one));
                }

                double az =  Aaz[i][k];
                double bz =  Bbz[i][k];
                double Jz =  LIM01(JJz[i][k]);
                Jz *= to_one;
                double Cz = sqrt(az * az + bz * bz);
                //log encoding
                if(logjz) {
                    double jmz =  Jz;
                    if (jmz > noise) {
                        double mm = applytojz(jmz);
                        double f = mm / jmz;
                        Jz *= f;
                        //Cz *= f;
                        Jz = LIM01(Jz);//clip values
                        //Cz = clipcz(Cz);
                    }
                }
                //sigmoid
                if(issigjz && iscie) {//sigmoid Jz
                    float val = Jz;
                    if(islogjz) {
                        val = std::max((xlog(Jz) / log2 - shadows_range) / (dynamic_range + 1.5), noise);//in range EV
                    }
                    if(sigmoidthjz >= 1.f) {
                        thjz = athjz * val + bthjz;//threshold
                    } else {
                        thjz = atjz * val + btjz;
                    }
                    sigmoidla (val, thjz, sigmjz);//sigmz "slope" of sigmoid
                    

                    Jz = LIM01((double) bljz * Jz + (double) val);
                }

                if(Qtoj == true) {//lightness instead of brightness
                    Jz /= to_one;
                    Jz /= maxjzw;//Jz white
                    Jz = SQR(Jz);
                }
                //contrast 
                Jz= LIM01(jz_contrast.getVal(LIM01(Jz)));
                //brightness and lightness
                if(lightreal > 0) {
                    Jz = LIM01(jz_light.getVal(Jz));
                }
                if(lightreal < 0) {
                    Jz = LIM01(jz_lightn.getVal(Jz));
                }
                //Jz (Jz) curve
                double Jzold = Jz;
                if(jzlocalcurve && localjzutili) {
                    Jz =  (double) (jzlocalcurve[(float) Jz * 65535.f] / 65535.f);
                    Jz  = 0.3 * (Jz - Jzold) + Jzold;
                }
                //reconvert from lightness or Brightness
                if(Qtoj == false) {
                    Jz /= to_one;
                } else {
                    Jz = sqrt(Jz);
                    Jz *= maxjzw;
                }

                double Hz;
                //remapping Cz 
                Hz = xatan2 ( bz, az );
                double Czold = Cz;
                //Cz(Cz) curve
                if(czlocalcurve && localczutili) {
                    Cz =  (double) (czlocalcurve[(float) Cz * 92666.f * (float) to_one] / (92666.f * (float) to_one));
                    Cz  = 0.5 * (Cz - Czold) + Czold;
                }
                //Cz(Jz) curve
                if(czjzlocalcurve && localczjzutili) { 
                    double chromaCfactor =  (double) (czjzlocalcurve[(float) Jz * 65535.f * (float) to_one]) / (Jz * 65535. * to_one);
                    Cz  *=  chromaCfactor;
                }
                //Hz in 0 2*PI
                if ( Hz < 0.0 ) {
                    Hz +=  (2. * rtengine::RT_PI);
                }
                //Chroma slider
                if(chromz < 0.) {
                    Cz = Cz * (1. + 0.01 * chromz);
                } else {
                    double maxcz = czlim / to_one;
                    double fcz = Cz / maxcz;
                    double pocz = pow(fcz , 1. - 0.0024 * chromz);//increase value - before 0.0017
                    Cz = maxcz * pocz;
                  //  Cz = Cz * (1. + 0.005 * chromz);//linear
                }
                //saturation slider
                if(saturz != 0.) {
                    double js = Jz/ maxjzw;//divide by Jz white
                    js = SQR(js);
                    if(js <= 0.) {
                        js = 0.0000001;
                    }
                    double Sz = Cz / (js);
                    if(saturz < 0.) {
                        Sz = Sz * (1. + 0.01 * saturz);
                    } else {
                        Sz = Sz * (1. + 0.003 * saturz);//not pow function because Sz is "open" - 0.003 empirical value to have results comparable to Cz
                    }
                    Cz = Sz * js;
                }
                
                //rotation hue
                Hz += dhue;
                if ( Hz < 0.0 ) {
                    Hz +=  (2. * rtengine::RT_PI);
                }
                Cz = clipcz(Cz);
                double2 sincosval = xsincos(Hz);
                az = clipazbz(Cz * sincosval.y);
                bz = clipazbz(Cz * sincosval.x);
                Cz = sqrt(az * az + bz * bz);


                bz = bz / (to_screen);
                az = az / (to_screen);

                Jz = LIM01(Jz / (to_screen));
                if(jabcie) {//Not used does not work at all
                    Jz = clipjz05(Jz);
                    gamutjz (Jz, az, bz, pl, wip, 0.94, 0.004);
                }

                double L_, M_, S_;
                double xx, yy, zz;
                bool zcam = z_cam;
                //reconvert to XYZ in double
                Ciecam02::jzczhzxyz (xx, yy, zz, Jz, az, bz, pl, L_, M_, S_, zcam);
                //re enable D50
                double x, y, z;
                x = 65535. * (d65_d50[0][0] * xx + d65_d50[0][1] * yy + d65_d50[0][2] * zz);
                y = 65535. * (d65_d50[1][0] * xx + d65_d50[1][1] * yy + d65_d50[1][2] * zz);
                z = 65535. * (d65_d50[2][0] * xx + d65_d50[2][1] * yy + d65_d50[2][2] * zz);

                float Ll, aa, bb;
                Color::XYZ2Lab(x,  y,  z, Ll, aa, bb);
                lab->L[i][k] = Ll;
                lab->a[i][k] = aa;
                lab->b[i][k] = bb;
            }
        }
    }

if(mocam == 0 || mocam == 1  || call == 1  || call == 2 || call == 10) {//call=2 vibrance warm-cool - call = 10 take into account "mean luminance Yb for Jz 
//begin ciecam
 if (settings->verbose && (mocam == 0 || mocam == 1  || call == 1)) {//display only if choice cam16
     //information on Cam16 scene conditions - allows user to see choices's incidences
    float maxicam = -1000.f;
    float maxicamq = -1000.f;
    float maxisat = -1000.f;
    float maxiM = -1000.f;
    float minicam = 1000000.f;
    float minicamq = 1000000.f;
    float minisat = 1000000.f;
    float miniM = 1000000.f;
    int nccam = 0;
    float sumcam = 0.f;
    float sumcamq = 0.f;
    float sumsat = 0.f;
    float sumM = 0.f;
    if(lp.logena && !(params->locallab.spots.at(sp).expcie && mocam == 1)) {//Log encoding only, but enable for log encoding if we use Cam16 module both with log encoding
        plum = 100.f;
    }


    
#ifdef _OPENMP
            #pragma omp parallel for reduction(min:minicam) reduction(max:maxicam) reduction(min:minicamq) reduction(max:maxicamq) reduction(min:minisat) reduction(max:maxisat) reduction(min:miniM) reduction(max:maxiM) reduction(+:sumcam) reduction(+:sumcamq) reduction(+:sumsat) reduction(+:sumM)if(multiThread)
#endif
        for (int i = 0; i < height; i+=1) {
            for (int k = 0; k < width; k+=1) {
                    float L = lab->L[i][k];
                    float a = lab->a[i][k];
                    float b = lab->b[i][k];
                    float x, y, z;
                    //convert Lab => XYZ
                    Color::Lab2XYZ(L, a, b, x, y, z);
                    x = x / 655.35f;
                    y = y / 655.35f;
                    z = z / 655.35f;
                    float J, C, h, Q, M, s;
                    Ciecam02::xyz2jchqms_ciecam02float(J, C,  h,
                                                   Q,  M,  s, aw, fl, wh,
                                                   x,  y,  z,
                                                   xw1, yw1,  zw1,
                                                   c,  nc, pow1, nbb, ncb, pfl, cz, d, c16, plum);
                    if(J > maxicam) {
                        maxicam = J;
                    }
                    if(J < minicam) {
                        minicam = J;
                    }
                    sumcam += J;

                    if(Q > maxicamq) {
                        maxicamq = Q;
                    }
                    if(Q < minicamq) {
                        minicamq = Q;
                    }
                    sumcamq += Q;

                    if(s > maxisat) {
                        maxisat = s;
                    }
                    if(s < minisat) {
                        minisat = s;
                    }
                    sumsat += s;

                    if(M > maxiM) {
                        maxiM = M;
                    }
                    if(M < miniM) {
                        miniM = M;
                    }
                    sumM += M;

            }
        }
        nccam = height * width;
        sumcam = sumcam / nccam;
        sumcamq /= nccam;
        sumsat /= nccam;
        sumM /= nccam;

        printf("Cam16 Scene  Lighness_J Brightness_Q- HDR-PQ=%5.1f minJ=%3.1f maxJ=%3.1f meanJ=%3.1f minQ=%3.1f maxQ=%4.1f meanQ=%4.1f\n", (double) plum, (double) minicam, (double) maxicam, (double) sumcam, (double) minicamq, (double) maxicamq, (double) sumcamq);
        printf("Cam16 Scene  Saturati-s Colorfulln_M- minSat=%3.1f maxSat=%3.1f meanSat=%3.1f minM=%3.1f maxM=%3.1f meanM=%3.1f\n", (double) minisat, (double) maxisat, (double) sumsat, (double) miniM, (double) maxiM, (double) sumM);
}

    float base = 10.;
    float linbase = 10.;
    float gray = 15.;
    if(islogq) {//with brightness Jz
        gray = 0.01f * (float) params->locallab.spots.at(sp).sourceGraycie;
        gray = pow_F(gray, 1.2f);//or 1.15 => modification to increase sensitivity gain, only on defaults, of course we can change this value manually...take into account suuround and Yb Cam16
        const float targetgraycie = params->locallab.spots.at(sp).targetGraycie;
        float targetgraycor = pow_F(0.01f * targetgraycie, 1.15f);
        base = targetgraycie > 1.f && targetgraycie < 100.f && (float) dynamic_range > 0.f ?  find_gray(std::abs((float) shadows_range) / (float) dynamic_range,(targetgraycor)) : 0.f;
        linbase = std::max(base, 2.f);//2. minimal base log to avoid very bad results
        if (settings->verbose) {
            printf("Base logarithm encoding Q=%5.1f\n", (double) linbase);
        }
    }

    const auto applytoq =
    [ = ](float x) -> float {

        x = rtengine::max(x, (float) noise);
        x = rtengine::max(x / gray, (float) noise);//gray = gain - before log conversion
        x = rtengine::max((xlogf(x) / log2f - (float) shadows_range) / (float) dynamic_range, (float) noise);//x in range EV
        assert(x == x);

        if (linbase > 0.f)//apply log base in function of targetgray blackEvjz and Dynamic Range
        {
            x = xlog2lin(x, linbase);
        }
        return x;
    };



//Ciecam "old" code not change except sigmoid added
#ifdef __SSE2__
        int bufferLength = ((width + 3) / 4) * 4; // bufferLength has to be a multiple of 4
#endif
#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
#ifdef __SSE2__
            // one line buffer per channel and thread
            float Jbuffer[bufferLength] ALIGNED16;
            float Cbuffer[bufferLength] ALIGNED16;
            float hbuffer[bufferLength] ALIGNED16;
            float Qbuffer[bufferLength] ALIGNED16;
            float Mbuffer[bufferLength] ALIGNED16;
            float sbuffer[bufferLength] ALIGNED16;
#endif
#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 16)
#endif
            for (int i = 0; i < height; i++) {
#ifdef __SSE2__
            // vectorized conversion from Lab to jchqms
                int k;
                vfloat c655d35 = F2V(655.35f);

                for (k = 0; k < width - 3; k += 4) {
                    vfloat x, y, z;
                    Color::Lab2XYZ(LVFU(lab->L[i][k]), LVFU(lab->a[i][k]), LVFU(lab->b[i][k]), x, y, z);
                    x = x / c655d35;
                    y = y / c655d35;
                    z = z / c655d35;
                    vfloat J, C, h, Q, M, s;
                    Ciecam02::xyz2jchqms_ciecam02float(J, C,  h,
                                                   Q,  M,  s, F2V(aw), F2V(fl), F2V(wh),
                                                   x,  y,  z,
                                                   F2V(xw1), F2V(yw1),  F2V(zw1),
                                                   F2V(c),  F2V(nc), F2V(pow1), F2V(nbb), F2V(ncb), F2V(pfl), F2V(cz), F2V(d), c16, F2V(plum));
                    STVF(Jbuffer[k], J);
                    STVF(Cbuffer[k], C);
                    STVF(hbuffer[k], h);
                    STVF(Qbuffer[k], Q);
                    STVF(Mbuffer[k], M);
                    STVF(sbuffer[k], s);
                }

                for (; k < width; k++) {
                    float L = lab->L[i][k];
                    float a = lab->a[i][k];
                    float b = lab->b[i][k];
                    float x, y, z;
                    //convert Lab => XYZ
                    Color::Lab2XYZ(L, a, b, x, y, z);
                    x = x / 655.35f;
                    y = y / 655.35f;
                    z = z / 655.35f;
                    float J, C, h, Q, M, s;
                    Ciecam02::xyz2jchqms_ciecam02float(J, C,  h,
                                                   Q,  M,  s, aw, fl, wh,
                                                   x,  y,  z,
                                                   xw1, yw1,  zw1,
                                                   c,  nc, pow1, nbb, ncb, pfl, cz, d, c16, plum);
                    Jbuffer[k] = J;
                    Cbuffer[k] = C;
                    hbuffer[k] = h;
                    Qbuffer[k] = Q;
                    Mbuffer[k] = M;
                    sbuffer[k] = s;
                }

#endif // __SSE2__

                for (int j = 0; j < width; j++) {
                    float J, C, h, Q, M, s;

#ifdef __SSE2__
                // use precomputed values from above
                    J = Jbuffer[j];
                    C = Cbuffer[j];
                    h = hbuffer[j];
                    Q = Qbuffer[j];
                    M = Mbuffer[j];
                    s = sbuffer[j];
#else
                    float x, y, z;
                    float L = lab->L[i][j];
                    float a = lab->a[i][j];
                    float b = lab->b[i][j];
                    float x1, y1, z1;
                    //convert Lab => XYZ
                    Color::Lab2XYZ(L, a, b, x1, y1, z1);
                    x = x1 / 655.35f;
                    y = y1 / 655.35f;
                    z = z1 / 655.35f;
                    //process source==> normal
                    Ciecam02::xyz2jchqms_ciecam02float(J, C,  h,
                                                   Q,  M,  s, aw, fl, wh,
                                                   x,  y,  z,
                                                   xw1, yw1,  zw1,
                                                   c,  nc, pow1, nbb, ncb, pfl, cz, d, c16, plum);
#endif
                    float Jpro, Cpro, hpro, Qpro, Mpro, spro;
                    Jpro = J;
                    Cpro = C;
                    hpro = h;
                    Qpro = Q;
                    Mpro = M;
                    spro = s;
                    /*
                    */
                    if(ciec) {
                        bool jp = false;
                        if ((cielocalcurve && localcieutili) && mecamcurve == 1) {
                            jp = true;
                            float Qq = Qpro * coefQ;
                            float Qold = Qpro;
                            Qq = 0.5f * cielocalcurve[Qq * 2.f];
                            Qq = Qq / coefQ;
                            Qpro = 0.2f * (Qq - Qold) + Qold;
                            if(jp) {
                                Jpro = SQR((10.f * Qpro) / wh);
                            }
                        }

                        Qpro = CAMBrightCurveQ[(float)(Qpro * coefQ)] / coefQ;   //brightness and contrast

                        if(islogq && issigq) {
                            float val =  Qpro *  coefq;;
                            if (val > (float) noise) {
                            float mm = applytoq(val);
                            float f = mm / val;
                            Qpro *=  f;
                            }
                        }


                        if(issigq && iscie && !islogq) {//sigmoid Q only with ciecam module
                            float val = Qpro * coefq;
                            if(sigmoidqj == true) {
                                val = std::max((xlog(val) / log2 - shadows_range) / (dynamic_range + 1.5), noise);//in range EV
                            }
                            if(sigmoidth >= 1.f) {
                                th = ath * val + bth;
                            } else {
                                th = at * val + bt;
                            }
                            sigmoidla (val, th, sigm);
                            float bl2 = 1.f;
                            Qpro = std::max(bl * Qpro + bl2 * val / coefq, 0.f);
                        }


                        float Mp, sres;
                        Mp = Mpro / 100.0f;
                        Ciecam02::curvecolorfloat(mchr, Mp, sres, 2.5f);
                        float dred = 100.f; //in C mode
                        float protect_red = 80.0f; // in C mode
                        dred *= coe; //in M mode
                        protect_red *= coe; //M mode
                        Color::skinredfloat(Jpro, hpro, sres, Mp, dred, protect_red, 0, rstprotection, 100.f, Mpro);
                        Jpro = SQR((10.f * Qpro) / wh);
                        Qpro = (Qpro == 0.f ? epsil : Qpro); // avoid division by zero
                        spro = 100.0f * sqrtf(Mpro / Qpro);

                        if (Jpro > 99.9f) {
                        Jpro = 99.9f;
                        }

                        Jpro = CAMBrightCurveJ[(float)(Jpro * 327.68f)];   //lightness CIECAM02 + contrast
                        float Sp = spro / 100.0f;
                        Ciecam02::curvecolorfloat(schr, Sp, sres, 1.5f);
                        dred = 100.f; // in C mode
                        protect_red = 80.0f; // in C mode
                        dred = 100.0f * sqrtf((dred * coe) / Q);
                        protect_red = 100.0f * sqrtf((protect_red * coe) / Q);
                        Color::skinredfloat(Jpro, hpro, sres, Sp, dred, protect_red, 0, rstprotection, 100.f, spro);
                        Qpro = QproFactor * sqrtf(Jpro);
                        float Cp = (spro * spro * Qpro) / (1000000.f);
                        Cpro = Cp * 100.f;
                        Ciecam02::curvecolorfloat(cchr, Cp, sres, 1.8f);
                        Color::skinredfloat(Jpro, hpro, sres, Cp, 55.f, 30.f, 1, rstprotection, 100.f, Cpro);
                        
                        hpro = hpro + hue;

                        if (hpro < 0.0f) {
                            hpro += 360.0f;    //hue
                        }
                        
                        if ((cielocalcurve && localcieutili) && mecamcurve == 0) {
                            float Jj = (float) Jpro * 327.68f;
                            float Jold = Jj;
                            Jj = 0.5f * cielocalcurve[Jj * 2.f];
                            Jj = 0.3f * (Jj - Jold) + Jold;    //divide sensibility
                            Jpro = (float)(Jj / 327.68f);

                            if (Jpro < 1.f) {
                                Jpro = 1.f;
                            }
                        }
                        if (cielocalcurve2 && localcieutili2) {
                            if(mecamcurve2 == 0) {
                                float parsat = 0.8f; //0.68;
                                float coef = 327.68f / parsat;
                                float Cc = (float) Cpro * coef;
                                float Ccold = Cc;
                                Cc = 0.5f * cielocalcurve2[Cc * 2.f];
                                float dred = 55.f;
                                float protect_red = 30.0f;
                                int sk1 = 1;
                                float ko = 1.f / coef;
                                Color::skinredfloat(Jpro, hpro, Cc, Ccold, dred, protect_red, sk1, rstprotection, ko, Cpro);
                            } else if (mecamcurve2 == 1) {
                                float parsat = 0.8f; //0.6
                                float coef = 327.68f / parsat;
                                float Ss = (float) spro * coef;
                                float Sold = Ss;
                                Ss = 0.5f * cielocalcurve2[Ss * 2.f];
                                Ss = 0.6f * (Ss - Sold) + Sold; //divide sensibility saturation
                                float dred = 100.f; // in C mode
                                float protect_red = 80.0f; // in C mode
                                dred = 100.0f * sqrtf((dred * coe) / Qpro);
                                protect_red = 100.0f * sqrtf((protect_red * coe) / Qpro);
                                float ko = 1.f / coef;
                                Color::skinredfloat(Jpro, hpro, Ss, Sold, dred, protect_red, 0, rstprotection, ko, spro);
                                Qpro = (4.0f / c) * sqrtf(Jpro / 100.0f) * (aw + 4.0f) ;
                                Cpro = (spro * spro * Qpro) / (10000.0f);
                            } else if (mecamcurve2 == 2) {
                                float parsat = 0.8f; //0.68;
                                float coef = 327.68f / parsat;
                                float Mm = (float) Mpro * coef;
                                float Mold = Mm;
                                Mm = 0.5f * cielocalcurve2[Mm * 2.f];
                                float dred = 100.f; //in C mode
                                float protect_red = 80.0f; // in C mode
                                dred *= coe; //in M mode
                                protect_red *= coe;
                                float ko = 1.f / coef;
                                Color::skinredfloat(Jpro, hpro, Mm, Mold, dred, protect_red, 0, rstprotection, ko, Mpro);
                                Cpro = Mpro / coe;
                            }
                        }
                        
                    }

                    //retrieve values C,J...s
                    C = Cpro;
                    J = Jpro;
                    Q = Qpro;
                    M = Mpro;
                    h = hpro;
                    s = spro;

#ifdef __SSE2__
                    // write to line buffers
                    Jbuffer[j] = J;
                    Cbuffer[j] = C;
                    hbuffer[j] = h;
#else
                    float xx, yy, zz;
                    //process normal==> viewing

                    Ciecam02::jch2xyz_ciecam02float(xx, yy, zz,
                                                J,  C, h,
                                                xw2, yw2,  zw2,
                                                c2, nc2,  pow1n, nbbj, ncbj, flj, czj, dj, awj, c16, plum);
                    x = CLIP(xx * 655.35f);
                    y = CLIP(yy * 655.35f);
                    z = CLIP(zz * 655.35f);
                    float Ll, aa, bb;
                    //convert xyz=>lab
                    Color::XYZ2Lab(x,  y,  z, Ll, aa, bb);
                    lab->L[i][j] = Ll;
                    lab->a[i][j] = aa;
                    lab->b[i][j] = bb;
#endif
                }

#ifdef __SSE2__
                // process line buffers
                float *xbuffer = Qbuffer;
                float *ybuffer = Mbuffer;
                float *zbuffer = sbuffer;

                for (k = 0; k < bufferLength; k += 4) {
                    vfloat x, y, z;
                    Ciecam02::jch2xyz_ciecam02float(x, y, z,
                                                LVF(Jbuffer[k]), LVF(Cbuffer[k]), LVF(hbuffer[k]),
                                                F2V(xw2), F2V(yw2), F2V(zw2),
                                                F2V(nc2), F2V(pow1n), F2V(nbbj), F2V(ncbj), F2V(flj), F2V(dj), F2V(awj), F2V(reccmcz), c16, F2V(plum));
                    STVF(xbuffer[k], x * c655d35);
                    STVF(ybuffer[k], y * c655d35);
                    STVF(zbuffer[k], z * c655d35);
                }

                // XYZ2Lab uses a lookup table. The function behind that lut is a cube root.
                // SSE can't beat the speed of that lut, so it doesn't make sense to use SSE
                for (int j = 0; j < width; j++) {
                    float Ll, aa, bb;
                    //convert xyz=>lab
                    xbuffer[j] = CLIP(xbuffer[j]);
                    ybuffer[j] = CLIP(ybuffer[j]);
                    zbuffer[j] = CLIP(zbuffer[j]);
                   
                    Color::XYZ2Lab(xbuffer[j], ybuffer[j], zbuffer[j], Ll, aa, bb);

                    lab->L[i][j] = Ll;
                    lab->a[i][j] = aa;
                    lab->b[i][j] = bb;
                }

#endif
            }

        }
    }
  
if(mocam == 3) {//Zcam not use but keep in case off 
/*
        double miniiz = 1000.;
        double maxiiz = -1000.;
        double sumiz = 0.;
        int nciz = 0;
        double epsilzcam = 0.0001;
        double atten = 2700.;
        double epsilzcam2 = 1.;
    if(mocam == 3) {//Zcam
        double pl = params->locallab.spots.at(sp).pqremap;
//calculate min, max, mean for Jz
#ifdef _OPENMP
            #pragma omp parallel for reduction(min:miniiz) reduction(max:maxiiz) reduction(+:sumiz) if(multiThread)
#endif
        for (int i = 0; i < height; i+=1) {
            for (int k = 0; k < width; k+=1) {
                float L = lab->L[i][k];
                float a = lab->a[i][k];
                float b = lab->b[i][k];
                float x, y, z;
                //convert Lab => XYZ
                Color::Lab2XYZ(L, a, b, x, y, z);
                x = x / 65535.f;
                y = y / 65535.f;
                z = z / 65535.f;
                double Jz, az, bz;
                double xx, yy, zz;
                //D50 ==> D65
                xx = (d50_d65[0][0] * (double) x + d50_d65[0][1] * (double) y + d50_d65[0][2] * (double) z);
                yy = (d50_d65[1][0] * (double) x + d50_d65[1][1] * (double) y + d50_d65[1][2] * (double) z);
                zz = (d50_d65[2][0] * (double) x + d50_d65[2][1] * (double) y + d50_d65[2][2] * (double) z);
                xx = LIM01(xx);
                yy = LIM01(yy);
                zz = LIM01(zz);
                
                double L_p, M_p, S_p;
                bool zcam = true;
                Ciecam02::xyz2jzczhz (Jz, az, bz, xx, yy, zz, pl, L_p, M_p, S_p, zcam);
                if(Jz > maxiiz) {
                    maxiiz = Jz;
                }
                if(Jz < miniiz) {
                    miniiz = Jz;
                }

                sumiz += Jz;

            }
        }
        nciz = height * width;
        sumiz = sumiz / nciz;
        sumiz += epsilzcam;
        maxiiz += epsilzcam;
        if (settings->verbose) {
            printf("Zcam miniiz=%f maxiiz=%f meaniz=%f\n", miniiz, maxiiz, sumiz);
        }
    }
    double avgmz = sumiz;
    //calculate various parameter for Zcam - those with ** come from documentation Zcam 
    // ZCAM, a colour appearance model based on a high dynamic range uniform colour space
    //Muhammad Safdar, Jon Yngve Hardeberg, and Ming Ronnier Luo
    // https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-29-4-6036&id=447640#e12
    double L_p, M_p, S_p;
    double jzw, azw, bzw;
    bool zcam = true;
    double plz = params->locallab.spots.at(sp).pqremap;// to test or change to 10000
//    double po = 0.1 + params->locallab.spots.at(sp).contthreszcam;
    float fb_source = sqrt(yb / 100.f);
    float fb_dest = sqrt(yb2 / 100.f);
    double flz = 0.171 * pow(la, 0.3333333)*(1. - exp(-(48. * (double) la / 9.)));
    double fljz = 0.171 * pow(la2, 0.3333333)*(1. - exp(-(48. * (double) la2 / 9.)));
    double cpow = 2.2;//empirical
    double cpp = pow( (double) c, 0.5);//empirical 
    double cpp2 = pow( (double) c2, 0.5);//empirical 
    double pfl = pow(flz, 0.25);
    double cmul_source = 1.26;//empirical 
    double cmul_source_ch = 1.1;//empirical 
    double achro_source =  pow((double) c, cpow)*(pow((double) flz, - 0.004)* (double) sqrt(fb_source));//I think there is an error in formula documentation step 5 - all parameters are inversed or wrong
    double achro_dest =  pow((double) c2, cpow)*(pow((double) fljz, - 0.004) * (double) sqrt(fb_dest));
    double  kk_source = (1.6 * (double) cpp) / pow((double) fb_source, 0.12);
    double  ikk_dest = pow((double) fb_dest, 0.12) /(1.6 * (double) cpp2);
    Ciecam02::xyz2jzczhz (jzw, azw, bzw, Xw, Yw, Zw, plz, L_p, M_p, S_p, zcam);
    double eff = 1.;
    double kap = 2.7;
    if(maxiiz > (kap * sumiz)) {
        kap = 1.7;
    }
    double qzw = cmul_source * atten * pow(jzw, (double) kk_source) /  achro_source;//I think there is an error in formula documentation step 5 - all parameters are inversed
    double maxforq = kap * sumiz * eff + epsilzcam2;
    if(maxforq > maxiiz) {
        maxforq = maxiiz;
    } else {
        maxforq = 0.9 * maxforq + 0.1 * maxiiz;
    }
    double qzmax = cmul_source * atten * pow(maxforq, (double) kk_source) /  achro_source;
    double izw = jzw;
    double coefm = pow(flz, 0.2) / (pow((double) fb_source, 0.1) * pow(izw, 0.78)); 
    if (settings->verbose) {
        printf("qzw=%f PL=%f qzmax=%f\n", qzw, plz, qzmax);//huge change with PQ peak luminance
    }
    array2D<double> Iiz(width, height);
    array2D<double> Aaz(width, height);
    array2D<double> Bbz(width, height);

//curve to replace LUT , LUT leads to crash...
        double contqz = 0.5 *  params->locallab.spots.at(sp).contqzcam;
        DiagonalCurve qz_contrast({
            DCT_NURBS,
            0, 0,
            avgmz - avgmz * (0.6 - contqz / 250.0), avgmz - avgmz * (0.6 + contqz / 250.0),
            avgmz + (1. - avgmz) * (0.6 - contqz / 250.0), avgmz + (1. - avgmz) * (0.6 + contqz / 250.0),
            1, 1
        });
        double contlz = 0.6 *  params->locallab.spots.at(sp).contlzcam;
        DiagonalCurve ljz_contrast({
            DCT_NURBS,
            0, 0,
            avgmz - avgmz * (0.6 - contlz / 250.0), avgmz - avgmz * (0.6 + contlz / 250.0),
            avgmz + (1. - avgmz) * (0.6 - contlz / 250.0), avgmz + (1. - avgmz) * (0.6 + contlz / 250.0),
            1, 1
        });

        //all calculations in double for best results...but slow
        double lqz = 0.4 *  params->locallab.spots.at(sp).lightqzcam;
        if(params->locallab.spots.at(sp).lightqzcam < 0) {
            lqz = 0.2 * params->locallab.spots.at(sp).lightqzcam; //0.4 less effect, no need 1.
        }
       DiagonalCurve qz_light({
            DCT_NURBS,
            0, 0,
            0.1, 0.1 + lqz / 150.,
            0.7, min (1.0, 0.7 + lqz / 300.0),
            1, 1
        });
        DiagonalCurve qz_lightn({
            DCT_NURBS,
            0, 0,
            max(0.0, 0.1  - lqz / 150.), 0.1 ,
            0.7 - lqz / 300.0, 0.7,
            1, 1
        });
        double ljz = 0.4 *  params->locallab.spots.at(sp).lightlzcam;
        if(params->locallab.spots.at(sp).lightlzcam < 0) {
            ljz = 0.2 * params->locallab.spots.at(sp).lightlzcam; 
        }
        DiagonalCurve ljz_light({
            DCT_NURBS,
            0, 0,
            0.1, 0.1 + ljz / 150.,
            0.7, min (1.0, 0.7 + ljz / 300.0),
            1, 1
        });
        DiagonalCurve ljz_lightn({
            DCT_NURBS,
            0, 0,
            max(0.0, 0.1  - ljz / 150.), 0.1 ,
            0.7 - ljz / 300.0, 0.7,
            1, 1
        });

#ifdef _OPENMP
            #pragma omp parallel for if(multiThread)
#endif
        for (int i = 0; i < height; i++) {
            for (int k = 0; k < width; k++) {
                float L = lab->L[i][k];
                float a = lab->a[i][k];
                float b = lab->b[i][k];
                float x, y, z;
                //convert Lab => XYZ
                Color::Lab2XYZ(L, a, b, x, y, z);
                x = x / 65535.f;
                y = y / 65535.f;
                z = z / 65535.f;
                double iz, az, bz;
                double xx, yy, zz;
                //change WP to D65
                xx = (d50_d65[0][0] * (double) x + d50_d65[0][1] * (double) y + d50_d65[0][2] * (double) z);
                yy = (d50_d65[1][0] * (double) x + d50_d65[1][1] * (double) y + d50_d65[1][2] * (double) z);
                zz = (d50_d65[2][0] * (double) x + d50_d65[2][1] * (double) y + d50_d65[2][2] * (double) z);
                double L_p, M_p, S_p;
                bool zcam = true;
                Ciecam02::xyz2jzczhz (iz, az, bz, xx, yy, zz, plz, L_p, M_p, S_p, zcam);
                Iiz[i][k] = LIM01(iz);
                Aaz[i][k] = clipazbz(az);
                Bbz[i][k] = clipazbz(bz);
            }
        }
#ifdef _OPENMP
            #pragma omp parallel  for if(multiThread)
#endif 
        for (int i = 0; i < height; i++) {
            for (int k = 0; k < width; k++) {

                double az =  Aaz[i][k];
                double bz =  Bbz[i][k];
                double iz =  Iiz[i][k];
                if(iz > kap * sumiz) {
                    iz = kap * sumiz * eff;
                }
                float coefqz = (float) qzmax;
                float coefjz = 100.f ;
                double qz = cmul_source * atten * pow(iz, (double) kk_source) / achro_source;//partial
                az *= cmul_source_ch;
                bz *= cmul_source_ch;
                
                qz= (double) coefqz * LIM01(qz_contrast.getVal((float)qz / coefqz));

                if(lqz > 0) {
                    qz = (double) coefqz * LIM01(qz_light.getVal((float)qz / coefqz));
                }
                if(lqz < 0) {
                    qz = (double) coefqz * LIM01(qz_lightn.getVal((float)qz / coefqz));
                }
               // double jz = 100. * (qz / qzw);
                double jz = SQR((10. * qz) / qzw);//formula CAM16
                jz= (double) coefjz * LIM01(ljz_contrast.getVal((float)jz / coefjz));
                if(ljz > 0) {
                    jz = (double) coefjz * LIM01(ljz_light.getVal((float)jz / coefjz));
                }
                if(ljz < 0) {
                    jz = (double) coefjz * LIM01(ljz_lightn.getVal((float)jz / coefjz));
                }
                if(jz > 100.) jz = 99.;
                
                
               //qzpro = 0.01 * jzpro * qzw;
                double qzpro = 0.1 * sqrt(jz) * qzw;
                iz = LIM01(pow(qzpro / (atten / achro_dest), ikk_dest));
                double h = atan2(bz, az);
                if ( h < 0.0 ) {
                    h += (double) (2.f * rtengine::RT_PI_F);
                }
                double hp = h * (360 / (double) (2.f * rtengine::RT_PI_F));
                double ez = 1.015 + cos(89.038 + hp);
                if(mchrz != 0.f || schrz != 0.f || cchrz != 0.f){
                    //colorfullness
                    double Mpz = 100. * pow(az * az + bz * bz, 0.37)* pow(ez, 0.068) * coefm;
                    Mpz *= (double) (1.f + 0.01f * mchrz);
                    float ccz = sqrt(pow((float) (Mpz / (100. * pow(ez, 0.068) * coefm)), (1.f / 0.37f)));
                    float2 sincosval = xsincosf(h);
                    az = (double)(ccz * sincosval.y);
                    bz = (double)(ccz * sincosval.x);
                    if(schrz != 0.f){
                    //saturation
                        double Spz = 100. * pow(flz, 0.6) * (Mpz / qz);
                        Spz *= (double) (1.f + 0.01f * schrz);
                        Mpz = (Spz * qz) / (100.* pow(flz, 0.6));
                        ccz = sqrt(pow((float) (Mpz / (100. * pow(ez, 0.068) * coefm)), (1.f / 0.37f)));
                        az = (double)(ccz * sincosval.y);
                        bz = (double)(ccz * sincosval.x);
                    }
                    if(cchrz != 0.f){
                     //   double Cpz = 100. * (Mpz / qzw);
                        double Cpz = 100. * (Mpz / pfl);//Cam16 formula
                        Cpz *= (double) (1.f + 0.01f * cchrz);
                        Mpz = (Cpz * pfl) / 100.;
                       // double Vpz = sqrt(SQR(jz - 58.) + 3.4 * SQR(Cpz));//vividness not working
                       // Vpz *= (double) (1.f + 0.01f * cchrz);
                        //Mpz = (Cpz * qzw) / 100.;
                       // Mpz  = 0.01 * qzw * sqrt((SQR(Vpz) - SQR(jz - 58.)) / 3.4);
                        ccz = sqrt(pow((float) (Mpz / (100. * pow(ez, 0.068) * coefm)), (1.f / 0.37f)));
                        az = (double)(ccz * sincosval.y);
                        bz = (double)(ccz * sincosval.x);
                    }
                    
                }    
                double L_, M_, S_;
                double xx, yy, zz;
                bool zcam = true;
                iz=LIM01(iz);
                az=clipazbz(az);
                bz=clipazbz(bz);
                
                Ciecam02::jzczhzxyz (xx, yy, zz, iz, az, bz, plz, L_, M_, S_, zcam);
                //re enable D50
                double x, y, z;
                x = 65535. * (d65_d50[0][0] * xx + d65_d50[0][1] * yy + d65_d50[0][2] * zz);
                y = 65535. * (d65_d50[1][0] * xx + d65_d50[1][1] * yy + d65_d50[1][2] * zz);
                z = 65535. * (d65_d50[2][0] * xx + d65_d50[2][1] * yy + d65_d50[2][2] * zz);

                float Ll, aa, bb;
                Color::XYZ2Lab(x,  y,  z, Ll, aa, bb);
                lab->L[i][k] = Ll;
                lab->a[i][k] = aa;
                lab->b[i][k] = bb;
            }
        }
*/
}

 
}

void ImProcFunctions::softproc(const LabImage* bufcolorig, const LabImage* bufcolfin, float rad, int bfh, int bfw, float epsilmax, float epsilmin, float thres, int sk, bool multiThread, int flag)
{
    if (rad != 0.f) {
        array2D<float> ble(bfw, bfh);
        array2D<float> guid(bfw, bfh);
        if (flag == 0) {

#ifdef _OPENMP
            #pragma omp parallel for if(multiThread)
#endif
            for (int ir = 0; ir < bfh; ir++) {
                for (int jr = 0; jr < bfw; jr++) {
                    guid[ir][jr] = Color::L2Y(bufcolorig->L[ir][jr]) / 32768.f;
                    ble[ir][jr] = Color::L2Y(bufcolfin->L[ir][jr]) / 32768.f;
                }
            }

            const float aepsil = (epsilmax - epsilmin) / 100.f;
            const float bepsil = epsilmin; //epsilmax - 100.f * aepsil;
           // const float epsil = aepsil * 0.1f * rad + bepsil;
            const float epsil = aepsil * rad + bepsil;
            const float blur = 10.f / sk * (thres + 0.f * rad);

            rtengine::guidedFilter(guid, ble, ble, blur, epsil, multiThread, 4);

#ifdef _OPENMP
            #pragma omp parallel for if(multiThread)
#endif
            for (int ir = 0; ir < bfh; ir++) {
                for (int jr = 0; jr < bfw; jr++) {
                    bufcolfin->L[ir][jr] =  Color::computeXYZ2LabY(32768.f * ble[ir][jr]);
                }
            }
        } else if (flag == 1) {

#ifdef _OPENMP
            #pragma omp parallel for if(multiThread)
#endif
            for (int ir = 0; ir < bfh; ir++)
                for (int jr = 0; jr < bfw; jr++) {
                    ble[ir][jr] = bufcolfin->L[ir][jr] / 32768.f;
                    guid[ir][jr] = bufcolorig->L[ir][jr] / 32768.f;
                }

            const float aepsil = (epsilmax - epsilmin) / 1000.f;
            const float bepsil = epsilmin; //epsilmax - 100.f * aepsil;
            const float epsil = rad < 0.f ? 0.0001f : aepsil * rad + bepsil;
            const float blur = rad < 0.f ? -1.f / rad : 1.f + rad;
            const int r2 = rtengine::max(int(25 / sk * blur + 0.5f), 1);

            rtengine::guidedFilter(guid, ble, ble, r2, epsil, multiThread);

#ifdef _OPENMP
            #pragma omp parallel for if(multiThread)
#endif
            for (int ir = 0; ir < bfh; ir++) {
                for (int jr = 0; jr < bfw; jr++) {
                    bufcolfin->L[ir][jr] =  32768.f * ble[ir][jr];
                }
            }
        } else if (flag == 2) {

#ifdef _OPENMP
            #pragma omp parallel for if(multiThread)
#endif
            for (int ir = 0; ir < bfh; ir++)
                for (int jr = 0; jr < bfw; jr++) {
                    ble[ir][jr] = bufcolfin->L[ir][jr] / 32768.f;
                    guid[ir][jr] = bufcolorig->L[ir][jr] / 32768.f;
                }

            const float aepsil = (epsilmax - epsilmin) / 1000.f;
            const float bepsil = epsilmin; //epsilmax - 100.f * aepsil;
            const float epsil = rad < 0.f ? 0.0001f : aepsil * 10.f * rad + bepsil;
          //  const float epsil =  bepsil;
            const float blur = rad < 0.f ? -1.f / rad : 0.00001f + rad;
            const int r2 = rtengine::max(int(20.f / sk * blur + 0.000001f), 1);

            rtengine::guidedFilter(guid, ble, ble, r2, epsil, multiThread);

#ifdef _OPENMP
            #pragma omp parallel for if(multiThread)
#endif
            for (int ir = 0; ir < bfh; ir++) {
                for (int jr = 0; jr < bfw; jr++) {
                    bufcolfin->L[ir][jr] =  32768.f * ble[ir][jr];
                }
            }
        }
    }
}


void ImProcFunctions::softprocess(const LabImage* bufcolorig, array2D<float> &buflight, float rad, int bfh, int bfw, double epsilmax, double epsilmin,  float thres, int sk, bool multiThread)
{
    float minlig = buflight[0][0];

#ifdef _OPENMP
    #pragma omp parallel for reduction(min:minlig) schedule(dynamic,16) if (multiThread)
#endif
    for (int ir = 0; ir < bfh; ir++) {
        for (int jr = 0; jr < bfw; jr++) {
            minlig = rtengine::min(buflight[ir][jr], minlig);
        }
    }

    array2D<float> guidsoft(bfw, bfh);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
    for (int ir = 0; ir < bfh; ir++) {
        for (int jr = 0; jr < bfw; jr++) {
            buflight[ir][jr] = LIM01((buflight[ir][jr] - minlig) / (100.f - minlig));
            guidsoft[ir][jr] = bufcolorig->L[ir][jr] / 32768.f;
        }
    }

    double aepsil = (epsilmax - epsilmin) / 90.0;
    double bepsil = epsilmax - 100.0 * aepsil;
    double epsil = aepsil * static_cast<double>(rad) + bepsil;
    float blur = 1.f / sk * (thres + 0.8f * rad);
    guidedFilter(guidsoft, buflight, buflight, blur, epsil,  multiThread, 4);


#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
    for (int ir = 0; ir < bfh; ir++) {
        for (int jr = 0; jr < bfw; jr++) {
            buflight[ir][jr] = (100.f - minlig) * buflight[ir][jr] + minlig;
        }
    }
}

void ImProcFunctions::exlabLocal(local_params& lp, float strlap, int bfh, int bfw, int bfhr, int bfwr, LabImage* bufexporig, LabImage* lab, const LUTf& hltonecurve, const LUTf& shtonecurve, const LUTf& tonecurve, const float hueref, const float lumaref, const float chromaref)
{
    //BENCHFUN
    //exposure local

    constexpr float maxran = 65536.f;
    if(lp.laplacexp == 0.f) {
        lp.linear = 0.f;
    }

    const float linear = lp.linear;
    int bw = bfw;
    int bh = bfh;
    if (linear > 0.f && lp.expcomp == 0.f) {
        lp.expcomp = 0.001f;
    }
    const bool exec = (lp.expmet == 1 && linear > 0.f && lp.laplacexp > 0.1f);

    if(!exec) {//for standard exposure
    
        const float cexp_scale = std::pow(2.f, lp.expcomp);
        const float ccomp = (rtengine::max(0.f, lp.expcomp) + 1.f) * lp.hlcomp / 100.f;
        const float cshoulder = ((maxran / rtengine::max(1.0f, cexp_scale)) * (lp.hlcompthr / 200.f)) + 0.1f;
        const float chlrange = maxran - cshoulder;
        const float diffde = 100.f - lp.sensex;//the more scope, the less take into account dE for Laplace
        if(!lp.invex) {// Laplacian not in inverse
            bw = bfwr;
            bh = bfhr;

            //Laplacian PDE before exposure to smooth L, algorithm exposure leads to increase L differences
            const std::unique_ptr<float[]> datain(new float[bfwr * bfhr]);
            const std::unique_ptr<float[]> dataout(new float[bfwr * bfhr]);
            const std::unique_ptr<float[]> dE(new float[bfwr * bfhr]);

            deltaEforLaplace(dE.get(), diffde, bfwr, bfhr, bufexporig, hueref, chromaref, lumaref);

            float alap = strlap * 600.f;
            float blap = strlap * 100.f;
            float aa = (alap - blap) / 50.f;
            float bb = blap - 30.f * aa;

            float lap;
            if (diffde > 80.f) {
                lap = alap;
            } else if (diffde < 30.f) {
                lap = blap;
            } else {
                lap = aa * diffde + bb;
            }

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = 0; y < bfhr; y++) {
                for (int x = 0; x < bfwr; x++) {
                    datain[y * bfwr + x] = bufexporig->L[y][x];
                }
            }

            MyMutex::MyLock lock(*fftwMutex);
            ImProcFunctions::retinex_pde(datain.get(), dataout.get(), bfwr, bfhr, lap, 1.f, dE.get(), 0, 1, 1);//350 arbitrary value about 45% strength Laplacian
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = 0; y < bfhr; y++) {
                for (int x = 0; x < bfwr; x++) {
                    bufexporig->L[y][x] = dataout[y * bfwr + x];
                }
            }
        
        }
 
#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif
        for (int ir = 0; ir < bh; ir++) {//for standard with Laplacian in normal and without in inverse 
            for (int jr = 0; jr < bw; jr++) {
                float L = bufexporig->L[ir][jr];
                //highlight
                const float hlfactor = (2 * L < MAXVALF ? hltonecurve[2 * L] : CurveFactory::hlcurve(cexp_scale, ccomp, chlrange, 2 * L));
                L *= hlfactor;//approximation but pretty good with Laplacian and L < mean, hl aren't call
                //shadow tone curve
                L *= shtonecurve[2 * L];
                //tonecurve
                lab->L[ir][jr] = 0.5f * tonecurve[2 * L];
            }
        }
    } else if(!lp.invex) {//for PDE algorithms
        constexpr float kl = 1.f;
        const float hlcompthr = lp.hlcompthr / 200.f;
        const float hlcomp = lp.hlcomp / 100.f;

#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif
        for (int ir = 0; ir < bfh; ir++) {
            for (int jr = 0; jr < bfw; jr++) {
                float L = bufexporig->L[ir][jr];
                const float Llin = LIM01(L / 32768.f);
                const float addcomp = linear * (-kl * Llin + kl);//maximum about 1 . IL
                const float exp_scale = pow_F(2.f, lp.expcomp + addcomp);
                const float shoulder = (maxran / rtengine::max(1.0f, exp_scale)) * hlcompthr + 0.1f;
                const float comp = (rtengine::max(0.f, (lp.expcomp + addcomp)) + 1.f) * hlcomp;
                const float hlrange = maxran - shoulder;

                //highlight
                const float hlfactor = (2 * L < MAXVALF ? hltonecurve[2 * L] : CurveFactory::hlcurve(exp_scale, comp, hlrange, 2 * L));
                L *= hlfactor * pow_F(2.f, addcomp);//approximation but pretty good with Laplacian and L < mean, hl aren't call
                //shadow tone curve
                L *= shtonecurve[2 * L];
                //tonecurve
                lab->L[ir][jr] = 0.5f * tonecurve[2 * L];
            }
        }
    }
}

void ImProcFunctions::addGaNoise(LabImage *lab, LabImage *dst, const float mean, const float variance, const int sk)
{
//   BENCHFUN
//Box-Muller method.
// add luma noise to image

    srand(1);

    const float variaFactor = SQR(variance) / sk;
    constexpr float randFactor1 = 1.f / RAND_MAX;
    constexpr float randFactor2 = (2.f * rtengine::RT_PI_F) / RAND_MAX;
#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        float z0, z1;
        bool generate = false;
#ifdef _OPENMP
        #pragma omp for schedule(static) // static scheduling is important to avoid artefacts
#endif
        for (int y = 0; y < lab->H; y++) {
            for (int x = 0; x < lab->W; x++) {
                generate = !generate;
                float kvar = 1.f;

                if (lab->L[y][x] < 12000.f) {
                    constexpr float ah = -0.5f / 12000.f;
                    constexpr float bh = 1.5f;
                    kvar = ah * lab->L[y][x] + bh;    //increase effect for low lights < 12000.f
                } else if (lab->L[y][x] > 20000.f) {
                    constexpr float ah = -0.5f / 12768.f;
                    constexpr float bh = 1.f - 20000.f * ah;
                    kvar = ah * lab->L[y][x] + bh;    //decrease effect for high lights > 20000.f
                    kvar = kvar < 0.5f ? 0.5f : kvar;
                }

                float varia = SQR(kvar) * variaFactor;

                if (!generate) {
                    dst->L[y][x] = LIM(lab->L[y][x] + mean + varia * z1, 0.f, 32768.f);
                    continue;
                }

                int u1 = 0;
                int u2;

                while (u1 == 0) {
                    u1 = rand();
                    u2 = rand();
                }

                float u1f = u1 * randFactor1;
                float u2f = u2 * randFactor2;

                float2 sincosval = xsincosf(2.f * rtengine::RT_PI_F * u2f);
                float factor = std::sqrt(-2.f * xlogf(u1f));
                z0 = factor * sincosval.y;
                z1 = factor * sincosval.x;

                dst->L[y][x] = LIM(lab->L[y][x] + mean + varia * z0, 0.f, 32768.f);

            }
        }
    }
}

void ImProcFunctions::DeNoise_Local(int call, const struct local_params& lp, LabImage* originalmask, int levred, float hueref, float lumaref, float chromaref, LabImage* original, LabImage* transformed, const LabImage &tmp1, int cx, int cy, int sk)
{
    //warning, but I hope used it next
    // local denoise and impulse
    //simple algo , perhaps we can improve as the others, but noise is here and not good for hue detection
    // BENCHFUN
    lumaref *= 327.68f;
    const float ach = lp.trans / 100.f;

    const float factnoise1 = 1.f + (lp.noisecf) / 500.f;
    const float factnoise2 = 1.f + (lp.noisecc) / 500.f;
    const float factnoise = factnoise1 * factnoise2;

    const int GW = transformed->W;
    const int GH = transformed->H;

    const float colorde = lp.colorde == 0 ? -1.f : lp.colorde; // -1.f to avoid black
    const float amplabL = 2.f * colorde;
    constexpr float darklim = 5000.f;

    const float refa = chromaref * std::cos(hueref) * 327.68f;
    const float refb = chromaref * std::sin(hueref) * 327.68f;
    const bool usemaskbl = lp.showmaskblmet == 2 || lp.enablMask || lp.showmaskblmet == 4;
    const bool blshow = lp.showmaskblmet == 1 || lp.showmaskblmet == 2;
    const bool previewbl = lp.showmaskblmet == 4;

    const std::unique_ptr<LabImage> origblur(new LabImage(GW, GH));
    const float radius = 3.f / sk;

    if (usemaskbl) {
#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            gaussianBlur(originalmask->L, origblur->L, GW, GH, radius);
            gaussianBlur(originalmask->a, origblur->a, GW, GH, radius);
            gaussianBlur(originalmask->b, origblur->b, GW, GH, radius);
        }
    } else {
#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            gaussianBlur(original->L, origblur->L, GW, GH, radius);
            gaussianBlur(original->a, origblur->a, GW, GH, radius);
            gaussianBlur(original->b, origblur->b, GW, GH, radius);
        }
    }

    const int begx = lp.xc - lp.lxL;
    const int begy = lp.yc - lp.lyT;
    constexpr float r327d68 = 1.f / 327.68f;

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        const LabImage* maskptr = origblur.get();
        const float mindE = 2.f + MINSCOPE * lp.sensden * lp.thr;
        const float maxdE = 5.f + MAXSCOPE * lp.sensden * (1 + 0.1f * lp.thr);
        const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
        const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif
        for (int y = 0; y < transformed->H; y++) {
            const int loy = cy + y;
            const bool isZone0 = loy > lp.yc + lp.ly || loy < lp.yc - lp.lyT; // whole line is zone 0 => we can skip a lot of processing

            if (isZone0) { // outside selection and outside transition zone => no effect, keep original values
                continue;
            }

            for (int x = 0, lox = cx + x; x < transformed->W; x++, lox++) {
                int zone;
                float localFactor = 1.f;

                if (lp.shapmet == 0) {
                    calcTransition(lox, loy, ach, lp, zone, localFactor);
                } else /*if (lp.shapmet == 1)*/ {
                    calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                }

                if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                    continue;
                }

                float reducdEL = 1.f;
                float reducdEa = 1.f;
                float reducdEb = 1.f;

                if (levred == 7) {
                    const float dEL = std::sqrt(0.9f * SQR(refa - maskptr->a[y][x]) + 0.9f * SQR(refb - maskptr->b[y][x]) + 1.2f * SQR(lumaref - maskptr->L[y][x])) * r327d68;
                    const float dEa = std::sqrt(1.2f * SQR(refa - maskptr->a[y][x]) + 1.f * SQR(refb - maskptr->b[y][x]) + 0.8f * SQR(lumaref - maskptr->L[y][x])) * r327d68;
                    const float dEb = std::sqrt(1.f * SQR(refa - maskptr->a[y][x]) + 1.2f * SQR(refb - maskptr->b[y][x]) + 0.8f * SQR(lumaref - maskptr->L[y][x])) * r327d68;
                    reducdEL = SQR(calcreducdE(dEL, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, lp.sensden));
                    reducdEa = SQR(calcreducdE(dEa, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, lp.sensden));
                    reducdEb = SQR(calcreducdE(dEb, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, lp.sensden));
                }
                float difL, difa, difb;

                if (call == 2  /*|| call == 1  || call == 3 */) { //simpleprocess
                    difL = tmp1.L[loy - begy][lox - begx] - original->L[y][x];
                    difa = tmp1.a[loy - begy][lox - begx] - original->a[y][x];
                    difb = tmp1.b[loy - begy][lox - begx] - original->b[y][x];
                } else  { //dcrop
                    const float repart = 1.0f - 0.01f * lp.reparden;
                    tmp1.L[y][x] = intp(repart, original->L[y][x], tmp1.L[y][x]);
                    tmp1.a[y][x] = intp(repart, original->a[y][x], tmp1.a[y][x]);
                    tmp1.b[y][x] = intp(repart, original->b[y][x], tmp1.b[y][x]);

                    difL = tmp1.L[y][x] - original->L[y][x];
                    difa = tmp1.a[y][x] - original->a[y][x];
                    difb = tmp1.b[y][x] - original->b[y][x];
                }

                difL *= localFactor * reducdEL;
                difa *= localFactor * reducdEa;
                difb *= localFactor * reducdEb;
                transformed->L[y][x] = CLIP(original->L[y][x] + difL);
                transformed->a[y][x] = clipC((original->a[y][x] + difa) * factnoise);
                transformed->b[y][x] = clipC((original->b[y][x] + difb) * factnoise) ;

                if (blshow) {
                    transformed->L[y][x] = CLIP(12000.f + amplabL * difL);// * 10.f empirical to can visualize modifications
                    transformed->a[y][x] = clipC(amplabL * difa);// * 10.f empirical to can visualize modifications
                    transformed->b[y][x] = clipC(amplabL * difb);// * 10.f empirical to can visualize modifications
                } else if (previewbl || lp.prevdE) {
                    const float difbdisp = (reducdEL + reducdEa + reducdEb) * 10000.f * colorde;

                    if (transformed->L[y][x] < darklim) { //enhance dark luminance as user can see!
                        transformed->L[y][x] = darklim - transformed->L[y][x];
                    }

                    if (colorde <= 0) {
                        transformed->a[y][x] = 0.f;
                        transformed->b[y][x] = difbdisp;
                    } else {
                        transformed->a[y][x] = -difbdisp;
                        transformed->b[y][x] = 0.f;
                    }
                }
            }
        }
    }
}

void ImProcFunctions::DeNoise_Local2(const struct local_params& lp, LabImage* originalmask, int levred, float hueref, float lumaref, float chromaref, LabImage* original, LabImage* transformed, const LabImage &tmp1, int cx, int cy, int sk)
{
    //warning, but I hope used it next
    // local denoise and impulse
    //simple algo , perhaps we can improve as the others, but noise is here and not good for hue detection
    // BENCHFUN
    const int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
    const int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
    const int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
    const int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
    
    
    lumaref *= 327.68f;
    const float ach = lp.trans / 100.f;

    const float factnoise1 = 1.f + (lp.noisecf) / 500.f;
    const float factnoise2 = 1.f + (lp.noisecc) / 500.f;
    const float factnoise = factnoise1 * factnoise2;

    const int GW = transformed->W;
    const int GH = transformed->H;

    const float colorde = lp.colorde == 0 ? -1.f : lp.colorde; // -1.f to avoid black
    const float amplabL = 2.f * colorde;
    constexpr float darklim = 5000.f;

    const float refa = chromaref * std::cos(hueref) * 327.68f;
    const float refb = chromaref * std::sin(hueref) * 327.68f;
    const bool usemaskbl = lp.showmaskblmet == 2 || lp.enablMask || lp.showmaskblmet == 4;
    const bool blshow = lp.showmaskblmet == 1 || lp.showmaskblmet == 2;
    const bool previewbl = lp.showmaskblmet == 4;

    const std::unique_ptr<LabImage> origblur(new LabImage(GW, GH));
    const float radius = 3.f / sk;

    if (usemaskbl) {
#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            gaussianBlur(originalmask->L, origblur->L, GW, GH, radius);
            gaussianBlur(originalmask->a, origblur->a, GW, GH, radius);
            gaussianBlur(originalmask->b, origblur->b, GW, GH, radius);
        }
    } else {
#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            gaussianBlur(original->L, origblur->L, GW, GH, radius);
            gaussianBlur(original->a, origblur->a, GW, GH, radius);
            gaussianBlur(original->b, origblur->b, GW, GH, radius);
        }
    }

 //   const int begx = lp.xc - lp.lxL;
 //   const int begy = lp.yc - lp.lyT;
    constexpr float r327d68 = 1.f / 327.68f;

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        const LabImage* maskptr = origblur.get();
        const float mindE = 2.f + MINSCOPE * lp.sensden * lp.thr;
        const float maxdE = 5.f + MAXSCOPE * lp.sensden * (1 + 0.1f * lp.thr);
        const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
        const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif

        for (int y = ystart; y < yend; y++) {
            const int loy = cy + y;
//            const bool isZone0 = loy > lp.yc + lp.ly || loy < lp.yc - lp.lyT; // whole line is zone 0 => we can skip a lot of processing

//            if (isZone0) { // outside selection and outside transition zone => no effect, keep original values
//                continue;
 //           }

            for (int x = xstart, lox = cx + x; x < xend; x++, lox++) {
                int zone;
                float localFactor = 1.f;

                if (lp.shapmet == 0) {
                    calcTransition(lox, loy, ach, lp, zone, localFactor);
                } else { /*if (lp.shapmet == 1)*/
                    calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                }

                if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                    continue;
                }

                float reducdEL = 1.f;
                float reducdEa = 1.f;
                float reducdEb = 1.f;

                if (levred == 7) {
                    const float dEL = std::sqrt(0.9f * SQR(refa - maskptr->a[y][x]) + 0.9f * SQR(refb - maskptr->b[y][x]) + 1.2f * SQR(lumaref - maskptr->L[y][x])) * r327d68;
                    const float dEa = std::sqrt(1.2f * SQR(refa - maskptr->a[y][x]) + 1.f * SQR(refb - maskptr->b[y][x]) + 0.8f * SQR(lumaref - maskptr->L[y][x])) * r327d68;
                    const float dEb = std::sqrt(1.f * SQR(refa - maskptr->a[y][x]) + 1.2f * SQR(refb - maskptr->b[y][x]) + 0.8f * SQR(lumaref - maskptr->L[y][x])) * r327d68;
                    reducdEL = SQR(calcreducdE(dEL, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, lp.sensden));
                    reducdEa = SQR(calcreducdE(dEa, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, lp.sensden));
                    reducdEb = SQR(calcreducdE(dEb, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, lp.sensden));
                }

                float difL, difa, difb;

                const float repart = 1.0f - 0.01f * lp.reparden;
                tmp1.L[y-ystart][x-xstart] = intp(repart, original->L[y][x], tmp1.L[y-ystart][x-xstart]);
                tmp1.a[y-ystart][x-xstart] = intp(repart, original->a[y][x], tmp1.a[y-ystart][x-xstart]);
                tmp1.b[y-ystart][x-xstart] = intp(repart, original->b[y][x], tmp1.b[y-ystart][x-xstart]);
              
                difL = tmp1.L[y-ystart][x-xstart] - original->L[y][x];
                difa = tmp1.a[y-ystart][x-xstart] - original->a[y][x];
                difb = tmp1.b[y-ystart][x-xstart] - original->b[y][x];

                difL *= localFactor * reducdEL;
                difa *= localFactor * reducdEa;
                difb *= localFactor * reducdEb;
                transformed->L[y][x] = CLIP(original->L[y][x] + difL);
                transformed->a[y][x] = clipC((original->a[y][x] + difa) * factnoise);
                transformed->b[y][x] = clipC((original->b[y][x] + difb) * factnoise) ;

                if (blshow) {
                    transformed->L[y][x] = CLIP(12000.f + amplabL * difL);// * 10.f empirical to can visualize modifications
                    transformed->a[y][x] = clipC(amplabL * difa);// * 10.f empirical to can visualize modifications
                    transformed->b[y][x] = clipC(amplabL * difb);// * 10.f empirical to can visualize modifications
                } else if (previewbl || lp.prevdE) {
                    const float difbdisp = (reducdEL + reducdEa + reducdEb) * 10000.f * colorde;

                    if (transformed->L[y][x] < darklim) { //enhance dark luminance as user can see!
                        transformed->L[y][x] = darklim - transformed->L[y][x];
                    }

                    if (colorde <= 0) {
                        transformed->a[y][x] = 0.f;
                        transformed->b[y][x] = difbdisp;
                    } else {
                        transformed->a[y][x] = -difbdisp;
                        transformed->b[y][x] = 0.f;
                    }
                }
            }
        }
    }
}


void ImProcFunctions::InverseReti_Local(const struct local_params & lp, const float hueref, const float chromaref,  const float lumaref, LabImage * original, LabImage * transformed, const LabImage * const tmp1, int cx, int cy, int chro, int sk)
{
    // BENCHFUN
//inverse local retinex
    float ach = lp.trans / 100.f;
    int GW = transformed->W;
    int GH = transformed->H;
    float refa = chromaref * cos(hueref);
    float refb = chromaref * sin(hueref);

    //balance deltaE
    const float kL = lp.balance / SQR(327.68f);
    const float kab = balancedeltaE(lp.balance) / SQR(327.68f);

    const std::unique_ptr<LabImage> origblur(new LabImage(GW, GH));

    float radius = 3.f / sk;
#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        gaussianBlur(original->L, origblur->L, GW, GH, radius);
        gaussianBlur(original->a, origblur->a, GW, GH, radius);
        gaussianBlur(original->b, origblur->b, GW, GH, radius);

    }
#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        const float mindE = 2.f + MINSCOPE * lp.sensh * lp.thr;
        const float maxdE = 5.f + MAXSCOPE * lp.sensh * (1 + 0.1f * lp.thr);
        const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
        const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif
        for (int y = 0; y < transformed->H; y++) {
            int loy = cy + y;

            for (int x = 0; x < transformed->W; x++) {
                int lox = cx + x;

                int zone;
                float localFactor;

                if (lp.shapmet == 0) {
                    calcTransition(lox, loy, ach, lp, zone, localFactor);
                } else /*if (lp.shapmet == 1)*/ {
                    calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                }

                float rL = origblur->L[y][x] / 327.68f;
                float dE = std::sqrt(kab * SQR(refa - origblur->a[y][x] / 327.68f) + kab * SQR(refb - origblur->b[y][x] / 327.68f) + kL * SQR(lumaref - rL));
                const float reducdE = calcreducdE(dE, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, lp.sensh);

                switch (zone) {
                    case 0: { // outside selection and outside transition zone => full effect, no transition
                        if (chro == 0) {
                            float difL = tmp1->L[y][x] - original->L[y][x];
                            transformed->L[y][x] = CLIP(original->L[y][x] + difL * reducdE);
                        }

                        if (chro == 1) {
                            float difa = tmp1->a[y][x] - original->a[y][x];
                            float difb = tmp1->b[y][x] - original->b[y][x];

                            transformed->a[y][x] = clipC(original->a[y][x] + difa * reducdE);
                            transformed->b[y][x] = clipC(original->b[y][x] + difb * reducdE);
                        }
                        break;
                    }

                    case 1: { // inside transition zone
                        float factorx = 1.f - localFactor;

                        if (chro == 0) {
                            float difL = tmp1->L[y][x] - original->L[y][x];
                            difL *= factorx;
                            transformed->L[y][x] = CLIP(original->L[y][x] + difL * reducdE);
                        }

                        if (chro == 1) {
                            float difa = tmp1->a[y][x] - original->a[y][x];
                            float difb = tmp1->b[y][x] - original->b[y][x];

                            difa *= factorx;
                            difb *= factorx;

                            transformed->a[y][x] = clipC(original->a[y][x] + difa * reducdE);
                            transformed->b[y][x] = clipC(original->b[y][x] + difb * reducdE);
                        }
                        break;
                    }

                    case 2: { // inside selection => no effect, keep original values
                        if (chro == 0) {
                            transformed->L[y][x] = original->L[y][x];
                        }

                        if (chro == 1) {
                            transformed->a[y][x] = original->a[y][x];
                            transformed->b[y][x] = original->b[y][x];
                        }
                    }
                }
            }
        }
    }
}



void ImProcFunctions::InverseBlurNoise_Local(LabImage * originalmask, const struct local_params & lp,  const float hueref, const float chromaref, const float lumaref, LabImage * original, LabImage * transformed, const LabImage * const tmp1, int cx, int cy, int sk)
{
    // BENCHFUN
//inverse local blur and noise
    float ach = lp.trans / 100.f;
    int GW = transformed->W;
    int GH = transformed->H;
    const float refa = chromaref * cos(hueref) * 327.68f;
    const float refb = chromaref * sin(hueref) * 327.68f;
    const float refL = lumaref * 327.68f;


    const bool blshow = (lp.showmaskblmet == 1 || lp.showmaskblmet == 2);
    const bool previewbl = (lp.showmaskblmet == 4);

    //balance deltaE
    const float kL = lp.balance / SQR(327.68f);
    const float kab = balancedeltaE(lp.balance) / SQR(327.68f);
    const float kH = lp.balanceh;
    const float kch = balancedeltaE(kH);

    const std::unique_ptr<LabImage> origblur(new LabImage(GW, GH));
    std::unique_ptr<LabImage> origblurmask;
    const bool usemaskbl = (lp.showmaskblmet == 2 || lp.enablMask || lp.showmaskblmet == 4);
    const bool usemaskall = usemaskbl;

    float radius = 3.f / sk;

    if (usemaskall) {
        origblurmask.reset(new LabImage(GW, GH));

#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            gaussianBlur(originalmask->L, origblurmask->L, GW, GH, radius);
            gaussianBlur(originalmask->a, origblurmask->a, GW, GH, radius);
            gaussianBlur(originalmask->b, origblurmask->b, GW, GH, radius);
        }
    }


#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        gaussianBlur(original->L, origblur->L, GW, GH, radius);
        gaussianBlur(original->a, origblur->a, GW, GH, radius);
        gaussianBlur(original->b, origblur->b, GW, GH, radius);

    }
#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        const LabImage *maskptr = usemaskall ? origblurmask.get() : origblur.get();
        const float mindE = 2.f + MINSCOPE * lp.sensbn * lp.thr;
        const float maxdE = 5.f + MAXSCOPE * lp.sensbn * (1 + 0.1f * lp.thr);
        const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
        const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif
        for (int y = 0; y < transformed->H; y++) {
            int loy = cy + y;

            for (int x = 0; x < transformed->W; x++) {
                int lox = cx + x;

                int zone;
                float localFactor;

                if (lp.shapmet == 0) {
                    calcTransition(lox, loy, ach, lp, zone, localFactor);
                } else {
                    calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                }
                float reducdE; 
                if (zone != 2) {
                    float abdelta2 = SQR(refa - maskptr->a[y][x]) + SQR(refb - maskptr->b[y][x]);
                    float chrodelta2 = SQR(std::sqrt(SQR(maskptr->a[y][x]) + SQR(maskptr->b[y][x])) - (chromaref * 327.68f));
                    float huedelta2 = abdelta2 - chrodelta2;
                    float dE = std::sqrt(kab * (kch * chrodelta2 + kH * huedelta2) + kL * SQR(refL - maskptr->L[y][x]));
                    reducdE = calcreducdE(dE, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, lp.sensbn);
                }

                switch (zone) {

                    case 0: { // outside selection and outside transition zone => full effect, no transition
                       const float diflc = (tmp1->L[y][x] - original->L[y][x]) * reducdE;
                       const float difa = (tmp1->a[y][x] - original->a[y][x]) * reducdE;
                       const float difb = (tmp1->b[y][x] - original->b[y][x]) * reducdE;
                       transformed->L[y][x] = CLIP(original->L[y][x] + diflc);
                       transformed->a[y][x] = clipC(original->a[y][x] + difa) ;
                       transformed->b[y][x] = clipC(original->b[y][x] + difb);

                        if (blshow) {
                            transformed->L[y][x] = CLIP(12000.f + diflc);
                            transformed->a[y][x] = clipC(difa);
                            transformed->b[y][x] = clipC(difb);
                        } else if (previewbl || lp.prevdE) {
                            transformed->a[y][x] = 0.f;
                            transformed->b[y] [x] = (difb);
                        }

                        break;
                    }

                    case 1: { // inside transition zone
                        const float factorx = 1.f - localFactor;

                        const float diflc = (tmp1->L[y][x] - original->L[y][x]) * (reducdE * factorx);
                        const float difa = (tmp1->a[y][x] - original->a[y][x]) * (reducdE * factorx);
                        const float difb = (tmp1->b[y][x] - original->b[y][x]) * (reducdE * factorx);
                        transformed->L[y][x] = CLIP(original->L[y][x] + diflc);
                        transformed->a[y][x] = clipC(original->a[y][x] + difa) ;
                        transformed->b[y][x] = clipC(original->b[y][x] + difb);

                        if (blshow) {
                            transformed->L[y][x] = CLIP(12000.f + diflc);
                            transformed->a[y][x] = clipC(difa);
                            transformed->b[y][x] = clipC(difb);
                        } else if (previewbl) {
                            transformed->a[y][x] = 0.f;
                            transformed->b[y][x] = (difb);
                        }

                        break;
                    }

                    case 2: { // outside selection and outside transition zone => no effect, keep original values
                        transformed->L[y][x] = original->L[y][x];
                        transformed->a[y][x] = original->a[y][x];
                        transformed->b[y][x] = original->b[y][x];
                    }
                }
            }
        }
    }
}



static void mean_fab(int xstart, int ystart, int bfw, int bfh, LabImage* bufexporig, int flag, const LabImage* original, float &fab, float &meanfab, float &maxfab, float chrom, bool multiThread)
{
    const int nbfab = bfw * bfh;

    meanfab = 0.f;
    fab = 50.f;

    if (nbfab > 0) {
        double sumab = 0.0;
#ifdef _OPENMP 
        #pragma omp parallel for reduction(+:sumab) if(multiThread)
#else
        static_cast<void>(multiThread);
#endif
        for (int y = 0; y < bfh; y++) {
            for (int x = 0; x < bfw; x++) {
                if(flag == 0) {
                    bufexporig->a[y][x] = original->a[y + ystart][x + xstart];
                    bufexporig->b[y][x] = original->b[y + ystart][x + xstart];
                } else {
                    bufexporig->a[y][x] = original->a[y][x];
                    bufexporig->b[y][x] = original->b[y][x];
                }
                sumab += static_cast<double>(SQR(std::fabs(bufexporig->a[y][x])) + SQR(std::fabs(bufexporig->b[y][x])));
                
              //  sumab += static_cast<double>(std::fabs(bufexporig->a[y][x]));
              //  sumab += static_cast<double>(std::fabs(bufexporig->b[y][x]));
            }
        }

        meanfab = sqrt(sumab / (2.0 * nbfab));

        double som = 0.0;
        double maxm = 0.0;
#ifdef _OPENMP
        #pragma omp parallel for reduction(+:som) reduction(max:maxfab)if(multiThread)
#endif
        for (int y = 0; y < bfh; y++) {
            for (int x = 0; x < bfw; x++) {
                som += static_cast<double>(SQR(std::fabs(bufexporig->a[y][x]) - meanfab) + SQR(std::fabs(bufexporig->b[y][x]) - meanfab));
                maxm = static_cast<double>(SQR(std::fabs(bufexporig->a[y][x])) + SQR(std::fabs(bufexporig->b[y][x])));
                maxm = sqrt(maxm);
                if(maxm > (double) maxfab) {
                    maxfab = (float) maxm;
                }

            }
        }

        const float multsigma = 3.f ;//(chrom >= 0.f ? 0.035f : 0.018f) * chrom + 1.f; //disabled an use 2 stddv
        const float kf = 0.f;
        const float multchrom = 1.f + kf * chrom; //small correction chrom here 0.f

        const float stddv = std::sqrt(som / nbfab);
        float fabprov = meanfab + multsigma * stddv * multchrom;//with 3 sigma about 99% cases
        if(fabprov > maxfab) {
            fabprov = maxfab;
        }
        fab = max(fabprov, 0.90f* maxfab);//Find maxi between mean + 3 sigma and 90% max (90 arbitrary empirical value)

        if (fab <= 0.f) {
            fab = 50.f;
        }
    }
}

struct grad_params {
    bool angle_is_zero, transpose, bright_top;
    float ta, yc, xc;
    float ys, ys_inv;
    float scale, botmul, topmul;
    float top_edge_0;
    int h;
};

void calclocalGradientParams(const struct local_params& lp, struct grad_params& gp, float ystart, float xstart, int bfw, int bfh, int indic)
{
    int w = bfw;
    int h = bfh;
    float stops = 0.f;
    float angs = 0.f;

    if (indic == 0) {
        stops = -lp.strmaexp;
        angs = lp.angmaexp;
    } else if (indic == 1) {
        stops = lp.strexp;
        angs = lp.angexp;
    } else if (indic == 2) {
        stops = lp.strSH;
        angs = lp.angSH;
    } else if (indic == 3) {
        stops = lp.strcol;
        angs = lp.angcol;
    } else if (indic == 4) {
        float redu = 1.f;

        if (lp.strcolab > 0.f) {
            redu = 0.6f;
        } else {
            redu = 0.15f;
        }

        stops = redu * lp.strcolab;
        angs = lp.angcol;
    } else if (indic == 5) {
        stops = lp.strcolab;
        angs = lp.angcol;
    } else if (indic == 6) {
        stops = lp.strcolh;
        angs = lp.angcol;
    } else if (indic == 7) {
        stops = lp.strvib;
        angs = lp.angvib;
    } else if (indic == 8) {
        float redu = 1.f;

        if (lp.strvibab > 0.f) {
            redu = 0.7f;
        } else {
            redu = 0.5f;
        }

        stops = redu * lp.strvibab;
        angs = lp.angvib;
    } else if (indic == 9) {
        stops = lp.strvibh;
        angs = lp.angvib;
    } else if (indic == 10) {
        stops = std::fabs(lp.strwav);
        angs = lp.angwav;
    } else if (indic == 11) {
        stops = lp.strlog;
        angs = lp.anglog;
    } else if (indic == 12) {
        stops = -lp.str_mas;
        angs = lp.ang_mas;
    }


    double gradient_stops = stops;
    double gradient_center_x = LIM01((lp.xc - xstart) / bfw);
    double gradient_center_y = LIM01((lp.yc - ystart) / bfh);
    double gradient_angle = static_cast<double>(angs) / 180.0 * rtengine::RT_PI;
    double varfeath = 0.01f * lp.feath;

    //printf("xstart=%f ysta=%f lpxc=%f lpyc=%f stop=%f bb=%f cc=%f ang=%f ff=%d gg=%d\n", xstart, ystart, lp.xc, lp.yc, gradient_stops, gradient_center_x, gradient_center_y, gradient_angle, w, h);

    // make 0.0 <= gradient_angle < 2 * rtengine::RT_PI
    gradient_angle = fmod(gradient_angle, 2 * rtengine::RT_PI);

    if (gradient_angle < 0.0) {
        gradient_angle += 2.0 * rtengine::RT_PI;
    }

    gp.bright_top = false;
    gp.transpose = false;
    gp.angle_is_zero = false;
    gp.h = h;
    double cosgrad = cos(gradient_angle);

    if (std::fabs(cosgrad) < 0.707) {
        // we transpose to avoid division by zero at 90 degrees
        // (actually we could transpose only for 90 degrees, but this way we avoid
        // division with extremely small numbers
        gp.transpose = true;
        gradient_angle += 0.5 * rtengine::RT_PI;
        double gxc = gradient_center_x;
        gradient_center_x = 1.0 - gradient_center_y;
        gradient_center_y = gxc;
    }

    gradient_angle = fmod(gradient_angle, 2 * rtengine::RT_PI);

    if (gradient_angle > 0.5 * rtengine::RT_PI && gradient_angle < rtengine::RT_PI) {
        gradient_angle += rtengine::RT_PI;
        gp.bright_top = true;
    } else if (gradient_angle >= rtengine::RT_PI && gradient_angle < 1.5 * rtengine::RT_PI) {
        gradient_angle -= rtengine::RT_PI;
        gp.bright_top = true;
    }

    if (std::fabs(gradient_angle) < 0.001 || std::fabs(gradient_angle - 2 * rtengine::RT_PI) < 0.001) {
        gradient_angle = 0;
        gp.angle_is_zero = true;
    }

    if (gp.transpose) {
        gp.bright_top = !gp.bright_top;
        std::swap(w, h);
    }

    gp.scale = 1.0 / pow(2, gradient_stops);

    if (gp.bright_top) {
        gp.topmul = 1.0;
        gp.botmul = gp.scale;
    } else {
        gp.topmul = gp.scale;
        gp.botmul = 1.0;
    }

    gp.ta = tan(gradient_angle);
    gp.xc = w * gradient_center_x;
    gp.yc = h * gradient_center_y;
    gp.ys = std::sqrt(h * h + w * w) * (varfeath / cos(gradient_angle));
    gp.ys_inv = 1.f / gp.ys;
    gp.top_edge_0 = gp.yc - gp.ys / 2.f;

    if (gp.ys < 1.f / h) {
        gp.ys_inv = 0;
        gp.ys = 0;
    }
}

void ImProcFunctions::blendstruc(int bfw, int bfh, LabImage* bufcolorig, float radius, float stru, array2D<float> & blend2, int sk, bool multiThread)
{
    SobelCannyLuma(blend2, bufcolorig->L, bfw, bfh, radius);
    float rm = 20.f / sk;

    if (rm > 0) {
        float **mb = blend2;
#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            gaussianBlur(mb, mb, bfw, bfh, rm);
        }
    }

    array2D<float> ble(bfw, bfh);
    array2D<float> guid(bfw, bfh);
#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif
    for (int ir = 0; ir < bfh; ir++) {
        for (int jr = 0; jr < bfw; jr++) {
            float X, Y, Z;
            float L = bufcolorig->L[ir][jr];
            float a = bufcolorig->a[ir][jr];
            float b = bufcolorig->b[ir][jr];
            Color::Lab2XYZ(L, a, b, X, Y, Z);

            guid[ir][jr] = Y / 32768.f;

            blend2[ir][jr] /= 32768.f;
        }
    }

    const float blur = 25 / sk * (10.f + 1.2f * stru);

    rtengine::guidedFilter(guid, blend2, ble, blur, 0.001, multiThread);

#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif
    for (int ir = 0; ir < bfh; ir++) {
        for (int jr = 0; jr < bfw; jr++) {
            blend2[ir][jr] = 32768.f * ble[ir][jr];
        }
    }
}


static void blendmask(const local_params& lp, int xstart, int ystart, int cx, int cy, int bfw, int bfh, LabImage* bufexporig, LabImage* original, LabImage* bufmaskor, LabImage* originalmas, float bl, float blab, int inv)
{
    bl /= 10.f;
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif
    for (int y = 0; y < bfh ; y++) {
        const int loy = y + ystart + cy;

        for (int x = 0; x < bfw; x++) {
            const int lox = x + xstart + cx;
            int zone;

            float localFactor = 1.f;
            const float achm = lp.trans / 100.f;

            if (lp.shapmet == 0) {
                calcTransition(lox, loy, achm, lp, zone, localFactor);
            } else /*if (lp.shapmet == 1)*/ {
                calcTransitionrect(lox, loy, achm, lp, zone, localFactor);
            }

            if (inv == 0) {
                if (zone > 0) {
                    bufexporig->L[y][x] += (bl * bufmaskor->L[y][x]);
                    bufexporig->a[y][x] *= (1.f + blab * bufmaskor->a[y][x]);
                    bufexporig->b[y][x] *= (1.f + blab * bufmaskor->b[y][x]);

                    bufexporig->L[y][x] = CLIP(bufexporig->L[y][x]);
                    bufexporig->a[y][x] = clipC(bufexporig->a[y][x]);
                    bufexporig->b[y][x] = clipC(bufexporig->b[y][x]);

                    originalmas->L[y][x] = CLIP(bufexporig->L[y][x] - bufmaskor->L[y][x]);
                    originalmas->a[y][x] = clipC(bufexporig->a[y][x] * (1.f - bufmaskor->a[y][x]));
                    originalmas->b[y][x] = clipC(bufexporig->b[y][x] * (1.f - bufmaskor->b[y][x]));

                    original->L[y + ystart][x + xstart] += (bl * localFactor * bufmaskor->L[y][x]);
                    original->a[y + ystart][x + xstart] *= (1.f + blab * localFactor * bufmaskor->a[y][x]);
                    original->b[y + ystart][x + xstart] *= (1.f + blab * localFactor * bufmaskor->b[y][x]);
                    original->L[y + ystart][x + xstart] = CLIP(original->L[y + ystart][x + xstart]);
                    original->a[y + ystart][x + xstart] = clipC(original->a[y + ystart][x + xstart]);
                    original->b[y + ystart][x + xstart] = clipC(original->b[y + ystart][x + xstart]);

                }
            } else if (inv == 1) {
                localFactor = 1.f - localFactor;

                if (zone < 2) {
                    bufexporig->L[y][x] += (bl * bufmaskor->L[y][x]);
                    bufexporig->a[y][x] *= (1.f + blab * bufmaskor->a[y][x]);
                    bufexporig->b[y][x] *= (1.f + blab * bufmaskor->b[y][x]);

                    bufexporig->L[y][x] = CLIP(bufexporig->L[y][x]);
                    bufexporig->a[y][x] = clipC(bufexporig->a[y][x]);
                    bufexporig->b[y][x] = clipC(bufexporig->b[y][x]);

                    originalmas->L[y][x] = CLIP(bufexporig->L[y][x] - bufmaskor->L[y][x]);
                    originalmas->a[y][x] = clipC(bufexporig->a[y][x] * (1.f - bufmaskor->a[y][x]));
                    originalmas->b[y][x] = clipC(bufexporig->b[y][x] * (1.f - bufmaskor->b[y][x]));

                    switch (zone) {
                        case 0: {
                            original->L[y + ystart][x + xstart] += (bl * bufmaskor->L[y][x]);
                            original->a[y + ystart][x + xstart] *= (1.f + blab * bufmaskor->a[y][x]);
                            original->b[y + ystart][x + xstart] *= (1.f + blab * bufmaskor->b[y][x]);
                            original->L[y + ystart][x + xstart] = CLIP(original->L[y + ystart][x + xstart]);
                            original->a[y + ystart][x + xstart] = clipC(original->a[y + ystart][x + xstart]);
                            original->b[y + ystart][x + xstart] = clipC(original->b[y + ystart][x + xstart]);
                            break;
                        }

                        case 1: {
                            original->L[y + ystart][x + xstart] += (bl * localFactor * bufmaskor->L[y][x]);
                            original->a[y + ystart][x + xstart] *= (1.f + blab * localFactor * bufmaskor->a[y][x]);
                            original->b[y + ystart][x + xstart] *= (1.f + blab * localFactor * bufmaskor->b[y][x]);
                            original->L[y + ystart][x + xstart] = CLIP(original->L[y + ystart][x + xstart]);
                            original->a[y + ystart][x + xstart] = clipC(original->a[y + ystart][x + xstart]);
                            original->b[y + ystart][x + xstart] = clipC(original->b[y + ystart][x + xstart]);
                        }

                    }
                }

            }
        }
    }
}

void ImProcFunctions::deltaEforMask(float **rdE, int bfw, int bfh, LabImage* bufcolorig, const float hueref, const float chromaref, const float lumaref,
                                    float maxdE, float mindE, float maxdElim,  float mindElim, float iterat, float limscope, int scope, float balance, float balanceh)
{
    const float refa = chromaref * cos(hueref);
    const float refb = chromaref * sin(hueref);
    const float refL = lumaref;

    const float kL = balance;
    const float kab = balancedeltaE(kL);
    const float kH = balanceh;
    const float kch = balancedeltaE(kH);

#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < bfh; y++) {
        for (int x = 0; x < bfw; x++) {
            const float abdelta2 = SQR(refa - bufcolorig->a[y][x] / 327.68f) + SQR(refb - bufcolorig->b[y][x] / 327.68f);
            const float chrodelta2 = SQR(std::sqrt(SQR(bufcolorig->a[y][x]) + SQR(bufcolorig->b[y][x])) / 327.68f - chromaref);
            const float huedelta2 = abdelta2 - chrodelta2;
            const float tempdE = std::sqrt(kab * (kch * chrodelta2 + kH * huedelta2) + kL * SQR(refL - bufcolorig->L[y][x] / 327.68f));

            float reducdE;
            if (tempdE > maxdE) {
                reducdE = 0.f;
            } else if (tempdE > mindE && tempdE <= maxdE) {
                const float ar = 1.f / (mindE - maxdE);
                const float br = - ar * maxdE;
                reducdE = pow(ar * tempdE + br, iterat);
            } else {
                reducdE = 1.f;
            }

            if (scope > limscope) {
                if (tempdE > maxdElim) {
                    reducdE = 0.f;
                } else if (tempdE > mindElim && tempdE <= maxdElim) {
                    const float arlim = 1.f / (mindElim - maxdElim);
                    const float brlim = - arlim * maxdElim;
                    const float reducdElim = pow(arlim * tempdE + brlim, iterat);
                    const float aalim = (1.f - reducdElim) / 20.f;
                    const float bblim = 1.f - 100.f * aalim;
                    reducdE = aalim * scope + bblim;
                } else {
                    reducdE = 1.f;
                }
            }

            rdE[y][x] = reducdE ;
        }
    }
}

static void showmask(int lumask, const local_params& lp, int xstart, int ystart, int cx, int cy, int bfw, int bfh, LabImage* bufexporig, LabImage* transformed, LabImage* bufmaskorigSH, int inv)
{
    float lum = fabs(lumask * 400.f);
    float colo = 0.f;
    if(lumask < 0.f) {
        lum *= 1.4f;
        colo = 30000.f + 12.f * lum;
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif
    for (int y = 0; y < bfh; y++) {
        const int loy = y + ystart + cy;

        for (int x = 0; x < bfw; x++) {
            const int lox = x + xstart + cx;
            int zone;
            float localFactor = 1.f;
            const float achm = lp.trans / 100.f;

            if (lp.shapmet == 0) {
                calcTransition(lox, loy, achm, lp, zone, localFactor);
            } else /*if (lp.shapmet == 1)*/ {
                calcTransitionrect(lox, loy, achm, lp, zone, localFactor);
            }

            if (inv == 0) {
                if (zone > 0) {//normal
                    transformed->L[y + ystart][x + xstart] = (lum) + clipLoc(bufmaskorigSH->L[y][x]);
                    transformed->a[y + ystart][x + xstart] = bufexporig->a[y][x] * bufmaskorigSH->a[y][x];
                    transformed->b[y + ystart][x + xstart] = (colo) + bufexporig->b[y][x] * bufmaskorigSH->b[y][x];
                }
            } else if (inv == 1) { //inverse
                if (zone == 0) {
                    transformed->L[y + ystart][x + xstart] = (lum) + clipLoc(bufmaskorigSH->L[y][x]);
                    transformed->a[y + ystart][x + xstart] = bufexporig->a[y][x] * bufmaskorigSH->a[y][x];
                    transformed->b[y + ystart][x + xstart] = (colo) + bufexporig->b[y][x] * bufmaskorigSH->b[y][x];
                }
            }

        }
    }
}
//from A.Griggio...very similar to discrete_laplacian_threhold...some differences with ceiling and data format
void ImProcFunctions::laplacian(const array2D<float> &src, array2D<float> &dst, int bfw, int bfh, float threshold, float ceiling, float factor, bool multiThread)
{
    const int W = bfw;
    const int H = bfh;

    const auto X =
        [W](int x) -> int
        {
            return x < 0 ? x+2 : (x >= W ? x-2 : x);
        };

    const auto Y =
        [H](int y) -> int
        {
            return y < 0 ? y+2 : (y >= H ? y-2 : y);
        };

    const auto get =
        [&src](int y, int x) -> float
        {
            return std::max(src[y][x], 0.f);
        };

    dst(W, H);
    const float f = factor / ceiling;

#ifdef _OPENMP
#   pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < H; ++y) {
        int n = Y(y-1), s = Y(y+1);
        for (int x = 0; x < W; ++x) {
            int w = X(x-1), e = X(x+1);
            float v = -8.f * get(y, x) + get(n, x) + get(s, x) + get(y, w) + get(y, e) + get(n, w) + get(n, e) + get(s, w) + get(s, e);
            dst[y][x] = LIM(std::abs(v) - threshold, 0.f, ceiling) * f;
        }
    }
}




void ImProcFunctions::discrete_laplacian_threshold(float * data_out, const float * data_in, size_t nx, size_t ny, float t)
{
   // BENCHFUN

    if (!data_in || !data_out) {
        fprintf(stderr, "a pointer is NULL and should not be so\n");
        abort();
    }

    /* pointers to the data and neighbour values */
    /*
     *                 y-1
     *             x-1 ptr x+1
     *                 y+1
     *    <---------------------nx------->
     */

#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif
    for (size_t j = 0; j < ny; j++) {
        const float* ptr_in = &data_in[j * nx];
        float* ptr_out = &data_out[j * nx];
        for (size_t i = 0; i < nx; i++) {
            float val = 0.f;
            /* row differences */
            if (0 < i) {
                const float diff = ptr_in[i] - ptr_in[i - 1];
                val += std::fabs(diff) > t ? diff : 0.f;
            }

            if (nx - 1 > i) {
                const float diff = ptr_in[i] - ptr_in[i + 1];;
                val += std::fabs(diff) > t ? diff : 0.f;
            }

            /* column differences */
            if (0 < j) {
                const float diff = ptr_in[i] - ptr_in[i - nx];;
                val += std::fabs(diff) > t ? diff : 0.f;
            }

            if (ny - 1 > j) {
                const float diff = ptr_in[i] - ptr_in[i + nx];;
                val += std::fabs(diff) > t ? diff : 0.f;
            }

            ptr_out[i] = val;
        }
    }

}

float *ImProcFunctions::cos_table(size_t size)
{
    float *table = NULL;

    /* allocate the cosinus table */
    if (NULL == (table = (float *) malloc(sizeof(*table) * size))) {
        fprintf(stderr, "allocation error\n");
        abort();
    }

    /*
     * fill the cosinus table,
     * table[i] = cos(i Pi / n) for i in [0..n[
     */
    const double pi_size = rtengine::RT_PI / size;

    for (size_t i = 0; i < size; i++) {
        table[i] = 1.0 - std::cos(pi_size * i);
    }

    return table;
}


void ImProcFunctions::rex_poisson_dct(float * data, size_t nx, size_t ny, double m)
{
    /*
     * Copyright 2009-2011 IPOL Image Processing On Line http://www.ipol.im/
     *

     * @file retinex_pde_lib.c discrete Poisson equation
     * @brief laplacian, DFT and Poisson routines
     *
     * @author Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
     * some adaptations for Rawtherapee
     */
   // BENCHFUN

    /*
     * get the cosinus tables
     * cosx[i] = cos(i Pi / nx) for i in [0..nx[
     * cosy[i] = cos(i Pi / ny) for i in [0..ny[
     */

    float* cosx = cos_table(nx);
    float* cosy = cos_table(ny);

    /*
     * we will now multiply data[i, j] by
     * m / (4 - 2 * cosx[i] - 2 * cosy[j]))
     * and set data[0, 0] to 0
     */
    float m2 = m / 2.;
    /*
     * after that, by construction, we always have
     * cosx[] + cosy[] != 2.
    */

#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif
    for (size_t i = 0; i < ny; ++i) {
        for (size_t j = 0; j < nx; ++j) {
            data[i * nx + j] *= m2 / (cosx[j] + cosy[i]);
        }
    }
    // handle the first value, data[0, 0] = 0
    data[0] = 0.f;

    free(cosx);
    free(cosy);

}

void ImProcFunctions::mean_dt(const float* data, size_t size, double& mean_p, double& dt_p)
{

    double mean = 0.;
    double dt = 0.;

#ifdef _OPENMP
    #pragma omp parallel for reduction(+:mean,dt) if(multiThread)
#endif
    for (size_t i = 0; i < size; i++) {
        mean += static_cast<double>(data[i]);
        dt += static_cast<double>(SQR(data[i]));
    }

    mean /= size;
    dt /= size;
    dt -= SQR(mean);
    mean_p = mean;
    dt_p = std::sqrt(dt);
}

void ImProcFunctions::normalize_mean_dt(float * data, const float * ref, size_t size, float mod, float sigm, float mdef, float sdef, float mdef2, float sdef2)
{
    /*
     * Copyright 2009-2011 IPOL Image Processing On Line http://www.ipol.im/
     *

     * @file retinex_pde_lib.c discrete Poisson equation
     * @brief laplacian, DFT and Poisson routines
     *
     * @author Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
     * adapted for Rawtherapee - jacques Desmis july 2019 - march 2021
     */

    if (NULL == data || NULL == ref) {
        fprintf(stderr, "a pointer is NULL and should not be so\n");
        abort();
    }

    double mean_ref, mean_data, dt_ref, dt_data;

    /* compute mean and variance of the two arrays */
    if(mdef!= 0.f && sdef != 0.f) {
        mean_ref = mdef;
        dt_ref = sdef;
    } else {
        mean_dt(ref, size, mean_ref, dt_ref);
    }
    if(mdef2!= 0.f && sdef2 != 0.f) {
       // printf("OK shortcut\n");
        mean_data = mdef2;
        dt_data = sdef2;
    } else {
        mean_dt(data, size, mean_data, dt_data);
    }

    /* compute the normalization coefficients */
    const double a = dt_ref / dt_data;
    const double b = mean_ref - a * mean_data;

    const float modma = static_cast<double>(mod) * a;
    const float sigmmmodmb = static_cast<double>(sigm) * static_cast<double>(mod) * b;
    const float onesmod = 1.f - mod;
    /* normalize the array */

#ifdef _OPENMP
    #pragma omp parallel for if(multiThread)
#endif
    for (size_t i = 0; i < size; i++) {
        data[i] = (modma * data[i] + sigmmmodmb) + onesmod * ref[i];//normalize mean and stdv and balance PDE
    }

}

void ImProcFunctions::retinex_pde(const float * datain, float * dataout, int bfw, int bfh, float thresh, float multy, float * dE, int show, int dEenable, int normalize)
{
    /*
     * Copyright 2009-2011 IPOL Image Processing On Line http://www.ipol.im/
     *

     * @file retinex_pde_lib.c discrete Poisson equation
     * @brief laplacian, DFT and Poisson routines
     *
     * @author Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
     * adapted for Rawtherapee by Jacques Desmis 6-2019 <jdesmis@gmail.com>
     */

   // BENCHFUN
   
#ifdef RT_FFTW3F_OMP
    if (multiThread) {
        fftwf_init_threads();
        fftwf_plan_with_nthreads(omp_get_max_threads());
    }
#endif

    float *datashow = nullptr;
    if (show != 0) {
        datashow = (float *) fftwf_malloc(sizeof(float) * bfw * bfh);
        if (!datashow) {
            fprintf(stderr, "allocation error\n");
            abort();
        }
    }

    float *data_tmp = (float *) fftwf_malloc(sizeof(float) * bfw * bfh);
    if (!data_tmp) {
        fprintf(stderr, "allocation error\n");
        abort();
    }

    //first call to laplacian with plein strength
    discrete_laplacian_threshold(data_tmp, datain, bfw, bfh, thresh);

    float *data_fft = (float *) fftwf_malloc(sizeof(float) * bfw * bfh);
    if (!data_fft) {
        fprintf(stderr, "allocation error\n");
        abort();
    }

    if (show == 1) {
        for (int y = 0; y < bfh ; y++) {
            for (int x = 0; x < bfw; x++) {
                datashow[y * bfw + x] = data_tmp[y * bfw + x];
            }
        }
    }

    //execute first
    const auto dct_fw = fftwf_plan_r2r_2d(bfh, bfw, data_tmp, data_fft, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    fftwf_execute(dct_fw);
    fftwf_destroy_plan(dct_fw);

    //execute second
    if (dEenable == 1) {
        float* data_fft04 = (float *)fftwf_malloc(sizeof(float) * bfw * bfh);
        float* data_tmp04 = (float *)fftwf_malloc(sizeof(float) * bfw * bfh);
        if (!data_fft04 || !data_tmp04) {
            fprintf(stderr, "allocation error\n");
            abort();
        }
        //second call to laplacian with 40% strength ==> reduce effect if we are far from ref (deltaE)
        discrete_laplacian_threshold(data_tmp04, datain, bfw, bfh, 0.4f * thresh);
        const auto dct_fw04 = fftwf_plan_r2r_2d(bfh, bfw, data_tmp04, data_fft04, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
        fftwf_execute(dct_fw04);
        fftwf_destroy_plan(dct_fw04);
        constexpr float exponent = 4.5f;

#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
#ifdef __SSE2__
            const vfloat exponentv = F2V(exponent);
#endif
#ifdef _OPENMP
            #pragma omp for
#endif
            for (int y = 0; y < bfh ; y++) {//mix two fftw Laplacian : plein if dE near ref
                int x = 0;
#ifdef __SSE2__
                for (; x < bfw - 3; x += 4) {
                    STVFU(data_fft[y * bfw + x], intp(pow_F(LVFU(dE[y * bfw + x]), exponentv), LVFU(data_fft[y * bfw + x]), LVFU(data_fft04[y * bfw + x])));
                }
#endif
                for (; x < bfw; x++) {
                    data_fft[y * bfw + x] = intp(pow_F(dE[y * bfw + x], exponent), data_fft[y * bfw + x], data_fft04[y * bfw + x]);
                }
            }
        }
        fftwf_free(data_fft04);
        fftwf_free(data_tmp04);
    }
    if (show == 2) {
        for (int y = 0; y < bfh ; y++) {
            for (int x = 0; x < bfw; x++) {
                datashow[y * bfw + x] = data_fft[y * bfw + x];
            }
        }
    }

    /* solve the Poisson PDE in Fourier space */
    /* 1. / (float) (bfw * bfh)) is the DCT normalisation term, see libfftw */
    rex_poisson_dct(data_fft, bfw, bfh, 1. / (double)(bfw * bfh));

    if (show == 3) {
        for (int y = 0; y < bfh ; y++) {
            for (int x = 0; x < bfw; x++) {
                datashow[y * bfw + x] = data_fft[y * bfw + x];
            }
        }
    }

    const auto dct_bw = fftwf_plan_r2r_2d(bfh, bfw, data_fft, data_tmp, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    fftwf_execute(dct_bw);
    fftwf_destroy_plan(dct_bw);
    fftwf_free(data_fft);

    if (show != 4 && normalize == 1) {
        normalize_mean_dt(data_tmp, datain, bfw * bfh, 1.f, 1.f, 0.f, 0.f, 0.f, 0.f);
    }

    if (show == 0 || show == 4) {

#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < bfh ; y++) {
            for (int x = 0; x < bfw; x++) {
                dataout[y * bfw + x] = clipLoc(multy * data_tmp[y * bfw + x]);
            }
        }
    } else if (show == 1 || show == 2 || show == 3) {
        for (int y = 0; y < bfh ; y++) {
            for (int x = 0; x < bfw; x++) {
                dataout[y * bfw + x] = clipLoc(multy * datashow[y * bfw + x]);
            }
        }
    }

    fftwf_free(data_tmp);
    if (datashow) {
        fftwf_free(datashow);
    }
    fftwf_cleanup();

#ifdef RT_FFTW3F_OMP
    if (multiThread) {
        fftwf_cleanup_threads();
    }
#endif
}

void ImProcFunctions::maskcalccol(bool invmask, bool pde, int bfw, int bfh, int xstart, int ystart, int sk, int cx, int cy, LabImage* bufcolorig, LabImage* bufmaskblurcol, LabImage* originalmaskcol, LabImage* original, LabImage* reserved, int inv, struct local_params & lp,
                                  float strumask, bool astool,
                                  const LocCCmaskCurve & locccmasCurve, bool lcmasutili,
                                  const LocLLmaskCurve & locllmasCurve, bool llmasutili,
                                  const LocHHmaskCurve & lochhmasCurve, bool lhmasutili, const LocHHmaskCurve & lochhhmasCurve, bool lhhmasutili,
                                  bool multiThread, bool enaMask, bool showmaske, bool deltaE, bool modmask, bool zero, bool modif, float chrom, float rad, float lap, float gamma, float slope, float blendm, float blendmab, int shado, int highl, float amountcd, float anchorcd,
                                  const LUTf& lmasklocalcurve, bool localmaskutili,
                                  const LocwavCurve & loclmasCurvecolwav, bool lmasutilicolwav, int level_bl, int level_hl, int level_br, int level_hr,
                                  int shortcu, bool delt, const float hueref, const float chromaref, const float lumaref,
                                  float maxdE, float mindE, float maxdElim,  float mindElim, float iterat, float limscope, int scope,
                                  bool fftt, float blu_ma, float cont_ma, int indic, float &fab
                                 )


{
    array2D<float> ble(bfw, bfh);
    array2D<float> blechro(bfw, bfh);
    array2D<float> hue(bfw, bfh);
    array2D<float> guid(bfw, bfh);
    const std::unique_ptr<LabImage> bufreserv(new LabImage(bfw, bfh));
    float meanfab, corfab;
    float maxfab = -1000.f;
    float epsi = 0.001f;
    mean_fab(xstart, ystart, bfw, bfh, bufcolorig, 0, original, fab, meanfab, maxfab, chrom, multiThread);
    corfab = 0.7f * (65535.f) / (fab + epsi);//empirical values 0.7 link to chromult

   // printf("Fab=%f corfab=%f maxfab=%f\n", (double) fab, (double) corfab, (double) maxfab);
    float chromult = 1.f;
    if(chrom > 0.f){
        chromult = 1.f + 0.003f * chrom;
    } else {
        chromult = 1.f + 0.01f * chrom;
    }
   // chromult * corfab * kmaskC 
    float kinv = 1.f;
    float kneg = 1.f;

    if (invmask) {
        kinv = 0.f;
        kneg = -1.f;
    }

    if (deltaE || modmask || enaMask || showmaske) {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
        for (int y = 0; y < bfh; y++) {
            for (int x = 0; x < bfw; x++) {
                bufmaskblurcol->L[y][x] = original->L[y + ystart][x + xstart];
                bufmaskblurcol->a[y][x] = original->a[y + ystart][x + xstart];
                bufmaskblurcol->b[y][x] = original->b[y + ystart][x + xstart];
                bufreserv->L[y][x] = reserved->L[y + ystart][x + xstart];
                bufreserv->a[y][x] = reserved->a[y + ystart][x + xstart];
                bufreserv->b[y][x] = reserved->b[y + ystart][x + xstart];
            }
        }

        JaggedArray<float> blendstru(bfw, bfh);

        if (blu_ma >= 0.25f && strumask == 0.f) {
            strumask = 0.1f; // to enable a small mask make FFT good ...why ??
        }

        if (strumask > 0.f) {
            float delstrumask = 4.1f - strumask;//4.1 = 2 * max slider strumask + 0.1
            buildBlendMask(bufcolorig->L, blendstru, bfw, bfh, delstrumask);
            float radblur = 0.02f * std::fabs(0.1f * rad);//empirical value
            float rm = radblur / sk;

            if (rm > 0) {
                float **mb = blendstru;
#ifdef _OPENMP
                #pragma omp parallel if (multiThread)
#endif
                {
                    gaussianBlur(mb, mb, bfw, bfh, rm);
                }
            }

        }

        JaggedArray<float> blendblur(bfw, bfh);

        JaggedArray<float> blur(bfw, bfh);

        if (cont_ma > 0.f) {
            float contra = cont_ma;
            buildBlendMask(bufcolorig->L, blendblur, bfw, bfh, contra);


            float radblur = 0.25f + 0.002f * std::fabs(rad);//empirical value
            float rm = radblur / sk;

            if (fftt) {
                if (rm < 0.3f) {
                    rm = 0.3f;
                }
            }

            if (rm > 0) {
                float **mb = blendblur;
#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
                {
                    gaussianBlur(mb, mb, bfw, bfh, rm);
                }
            }

            if (blu_ma >= 0.25f) {
                if (!fftt) { // || (lp.fftColorMask && call != 2)) {
#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
                    {
                        gaussianBlur(bufcolorig->L, blur, bfw, bfh, blu_ma / sk);
                    }
                } else {
                    ImProcFunctions::fftw_convol_blur2(bufcolorig->L, blur, bfw, bfh, blu_ma / sk, 0, 0);
                }

                for (int i = 0; i < bfh; i++) {
                    for (int j = 0; j < bfw; j++) {
                        blur[i][j] = intp(blendblur[i][j], bufcolorig->L[i][j], rtengine::max(blur[i][j], 0.0f));
                    }
                }
            }
        }

        bool HHmaskcurve = false;

        if (lochhhmasCurve && lhhmasutili) {
            for (int i = 0; i < 500; i++) {
                if (lochhhmasCurve[i] != 0.5f) {
                    HHmaskcurve = true;
                    break;
                }
            }
        }

//denoise mask chroma


        LabImage tmpab(bfw, bfh);
        tmpab.clear(true);

#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
        for (int ir = 0; ir < bfh; ir++)
            for (int jr = 0; jr < bfw; jr++) {
                tmpab.L[ir][jr] = bufcolorig->L[ir][jr];
                tmpab.a[ir][jr] = bufcolorig->a[ir][jr];
                tmpab.b[ir][jr] = bufcolorig->b[ir][jr];
        }
        float noisevarab_r = SQR(lp.denoichmask / 10.f);
        if(noisevarab_r > 0.f) {
            int wavelet_leve = 6;

            int minwin1 = rtengine::min(bfw, bfh);
            int maxlevelspot1 = 9;

            while ((1 << maxlevelspot1) >= (minwin1 * sk) && maxlevelspot1  > 1) {
                --maxlevelspot1 ;
            }

            wavelet_leve = rtengine::min(wavelet_leve, maxlevelspot1);
            int maxlvl1 = wavelet_leve;
#ifdef _OPENMP
            const int numThreads = omp_get_max_threads();
#else
            const int numThreads = 1;

#endif

            wavelet_decomposition Ldecomp(tmpab.L[0],tmpab.W, tmpab.H, maxlvl1, 1, sk, numThreads, lp.daubLen);
            wavelet_decomposition adecomp(tmpab.a[0],tmpab.W, tmpab.H, maxlvl1, 1, sk, numThreads, lp.daubLen);
            wavelet_decomposition bdecomp(tmpab.b[0],tmpab.W, tmpab.H, maxlvl1, 1, sk, numThreads, lp.daubLen);
            float* noisevarchrom;
            noisevarchrom = new float[bfw*bfh];
            float nvch = 0.6f;//high value
            float nvcl = 0.1f;//low value
            float seuil = 4000.f;//low
            float seuil2 = 15000.f;//high
            //ac and bc for transition
            float ac = (nvch - nvcl) / (seuil - seuil2);
            float bc = nvch - seuil * ac;
            int bfw2 = (bfw + 1) / 2;
             
#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
            for (int ir = 0; ir < bfh; ir++)
                for (int jr = 0; jr < bfw; jr++) {
                    float cN = std::sqrt(SQR(tmpab.a[ir][jr]) + SQR(tmpab.b[ir][jr]));

                    if (cN < seuil) {
                        noisevarchrom[(ir >> 1)*bfw2 + (jr >> 1)] =  nvch;
                    } else if (cN < seuil2) {
                        noisevarchrom[(ir >> 1)*bfw2 + (jr >> 1)] = ac * cN + bc;
                    } else {
                        noisevarchrom[(ir >> 1)*bfw2 + (jr >> 1)] =  nvcl;
                    }
                }
           
            float madL[8][3];
            int levred = maxlvl1;
            if (!Ldecomp.memory_allocation_failed()) {
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic) collapse(2) if (multiThread)
#endif
                for (int lvl = 0; lvl < levred; lvl++) {
                    for (int dir = 1; dir < 4; dir++) {
                        int Wlvl_L = Ldecomp.level_W(lvl);
                        int Hlvl_L = Ldecomp.level_H(lvl);
                        const float* const* WavCoeffs_L = Ldecomp.level_coeffs(lvl);

                        madL[lvl][dir - 1] = SQR(Mad(WavCoeffs_L[dir], Wlvl_L * Hlvl_L));
                    }
                }
            }
            
            if (!adecomp.memory_allocation_failed() && !bdecomp.memory_allocation_failed()) {
                        WaveletDenoiseAll_BiShrinkAB(Ldecomp, adecomp, noisevarchrom, madL, nullptr, 0, noisevarab_r, true, false, false, numThreads);
                        WaveletDenoiseAllAB(Ldecomp, adecomp, noisevarchrom, madL, nullptr, 0, noisevarab_r, true, false, false, numThreads);

                        WaveletDenoiseAll_BiShrinkAB(Ldecomp, bdecomp, noisevarchrom, madL, nullptr, 0, noisevarab_r, false, false, false, numThreads);
                        WaveletDenoiseAllAB(Ldecomp, bdecomp, noisevarchrom, madL, nullptr, 0, noisevarab_r, false, false, false, numThreads);

            }
           
            delete[] noisevarchrom;

            if (!Ldecomp.memory_allocation_failed()) {
                Ldecomp.reconstruct(tmpab.L[0]);
            }
            if (!adecomp.memory_allocation_failed()) {
                adecomp.reconstruct(tmpab.a[0]);
            }
            if (!bdecomp.memory_allocation_failed()) {
                bdecomp.reconstruct(tmpab.b[0]);
            }
            
            float meanfab1, fab1, maxfab1;
            std::unique_ptr<LabImage> buforig;
            buforig.reset(new LabImage(bfw, bfh));
#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
            for (int ir = 0; ir < bfh; ir++)
                for (int jr = 0; jr < bfw; jr++) {
                    buforig->L[ir][jr] = tmpab.L[ir][jr];
                    buforig->a[ir][jr] = tmpab.a[ir][jr];
                    buforig->b[ir][jr] = tmpab.b[ir][jr];
                }

    
            mean_fab(xstart, ystart, bfw, bfh, buforig.get(), 1, buforig.get(), fab1, meanfab1, maxfab1, chrom, multiThread);
              //  printf("Fab den=%f \n", (double) fab1);
            fab = fab1;//fab denoise
            
        }
// end code denoise mask chroma




#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
#ifdef __SSE2__
            float atan2Buffer[bfw] ALIGNED64;
//            float atan2BufferH[bfw] ALIGNED64;
#endif
#ifdef _OPENMP
            #pragma omp for schedule(dynamic, 16)
#endif
            for (int ir = 0; ir < bfh; ir++) {
#ifdef __SSE2__

                if (lochhmasCurve && lhmasutili) {
                    int i = 0;

                    for (; i < bfw - 3; i += 4) {
                     //   STVF(atan2Buffer[i], xatan2f(LVFU(bufcolorig->b[ir][i]), LVFU(bufcolorig->a[ir][i])));
                        STVF(atan2Buffer[i], xatan2f(LVFU(tmpab.b[ir][i]), LVFU(tmpab.a[ir][i])));
                    }

                    for (; i < bfw; i++) {
                      //  atan2Buffer[i] = xatan2f(bufcolorig->b[ir][i], bufcolorig->a[ir][i]);
                        atan2Buffer[i] = xatan2f(tmpab.b[ir][i], tmpab.a[ir][i]);
                    }
                }

#endif

                for (int jr = 0; jr < bfw; jr++) {
                    float kmaskL = 0.f;
                    float kmaskC = 0.f;
                    float kmaskHL = 0.f;
                    float kmaskH = 0.f;
                    float kmasstru = 0.f;
                    float kmasblur = 0.f;

                    if (strumask > 0.f && !astool) {
                        kmasstru = bufcolorig->L[ir][jr] * blendstru[ir][jr];
                    }

                    if (cont_ma > 0.f) {

                        if (blu_ma >= 0.25f) {

                            float prov = intp(blendstru[ir][jr], bufcolorig->L[ir][jr], rtengine::max(blur[ir][jr], 0.0f));
                            kmasblur = bufcolorig->L[ir][jr] - prov ;

                        }
                    }

                    if (locllmasCurve && llmasutili) {
                       // printf("s");
                        kmaskL = 32768.f * LIM01(kinv - kneg * locllmasCurve[(500.f / 32768.f) * bufcolorig->L[ir][jr]]);

                    }

                    if (!deltaE && locccmasCurve && lcmasutili) {
                      //  kmaskC = LIM01(kinv  - kneg * locccmasCurve[500.f * (0.0001f + std::sqrt(SQR(bufcolorig->a[ir][jr]) + SQR(bufcolorig->b[ir][jr])) / (fab))]);
                        kmaskC = LIM01(kinv  - kneg * locccmasCurve[500.f * (0.0001f + std::sqrt(SQR(tmpab.a[ir][jr]) + SQR(tmpab.b[ir][jr])) / fab)]);
                    }

                    if (lochhmasCurve && lhmasutili) {
#ifdef __SSE2__
                        const float huema = atan2Buffer[jr];
#else
                       // const float huema = xatan2f(bufcolorig->b[ir][jr], bufcolorig->a[ir][jr]);
                        const float huema = xatan2f(tmpab.b[ir][jr], tmpab.a[ir][jr]);
#endif
                        float h = Color::huelab_to_huehsv2(huema);
                        h += 1.f / 6.f;

                        if (h > 1.f) {
                            h -= 1.f;
                        }

                        const float valHH = LIM01(kinv - kneg * lochhmasCurve[500.f *  h]);

                        if (!deltaE) {
                            kmaskH = valHH;
                        }

                        kmaskHL = 32768.f * valHH;
                    }

                    bufmaskblurcol->L[ir][jr] = clipLoc(kmaskL + kmaskHL + kmasstru + kmasblur);
                    bufmaskblurcol->a[ir][jr] = clipC((chromult * corfab * kmaskC + chromult * kmaskH));
                    bufmaskblurcol->b[ir][jr] = clipC((chromult * corfab * kmaskC + chromult * kmaskH));

                    if (shortcu == 1) { //short circuit all L curve
                        bufmaskblurcol->L[ir][jr] = 32768.f - bufcolorig->L[ir][jr];
                    }

                    ble[ir][jr] = bufmaskblurcol->L[ir][jr] / 32768.f;
                    hue[ir][jr] = xatan2f(bufmaskblurcol->b[ir][jr], bufmaskblurcol->a[ir][jr]);
                    const float chromah = std::sqrt(SQR(bufmaskblurcol->b[ir][jr]) + SQR(bufmaskblurcol->a[ir][jr]));

                    blechro[ir][jr] = chromah / 32768.f;//must be good perhaps more or less, only incidence on LIM blea bleb
                    guid[ir][jr] = Color::L2Y(bufcolorig->L[ir][jr]) / 32768.f;

                }
            }
        }
        
        if (lap > 0.f && pde) {
            array2D<float> mask;
            mask(bfw, bfh);
            float amount = LIM01(float(lap)/100.f);
            array2D<float> LL(bfw, bfh, bufcolorig->L, ARRAY2D_BYREFERENCE);
            laplacian(LL, mask, bfw, bfh, 25.f, 20000.f, amount, false);
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int i = 0; i < bfh; ++i) {
                for (int j = 0; j < bfw; ++j) {
                    mask[i][j] = LIM01(mask[i][j]);
                }
            }
            for (int i = 0; i < 3; ++i) {
                boxblur(static_cast<float**>(mask), static_cast<float**>(mask), 5 / sk, bfw, bfh, false);
            }
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int i = 0; i < bfh; ++i) {
                for (int j = 0; j < bfw; ++j) {
                    bufmaskblurcol->L[i][j] += clipLoc(100000.f * (mask[i][j]));//increase strongly result
                }
            }
        }

        std::unique_ptr<LabImage> bufprov;
        if (delt) {
            bufprov.reset(new LabImage(bfw, bfh));
            bufprov->CopyFrom(bufmaskblurcol, multiThread);
        }

        if (rad != 0.f) {
            const float tmpblur = rad < 0.f ? -1.f / rad : 1.f + rad;
            const int r1 = rtengine::max<int>(4 / sk * tmpblur + 0.5f, 1);
            const int r2 = rtengine::max<int>(25 / sk * tmpblur + 0.5f, 1);

            constexpr float epsilmax = 0.005f;
            constexpr float epsilmin = 0.00001f;

            constexpr float aepsil = (epsilmax - epsilmin) / 100.f;
            constexpr float bepsil = epsilmin; //epsilmax - 100.f * aepsil;
            const float epsil = rad < 0.f ? 0.001f : aepsil * rad + bepsil;

            rtengine::guidedFilter(guid, blechro, blechro, r1, epsil, multiThread);
            rtengine::guidedFilter(guid, ble, ble, r2, 0.2f * epsil, multiThread);
        }

        LUTf lutTonemaskexp(65536);
        calcGammaLut(gamma, slope, lutTonemaskexp);

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
        for (int ir = 0; ir < bfh; ir++) {
            for (int jr = 0; jr < bfw; jr++) {
                float2 sincosval = xsincosf(hue[ir][jr]);
                bufmaskblurcol->L[ir][jr] = lutTonemaskexp[ble[ir][jr] * 65536.f];
                bufmaskblurcol->a[ir][jr] = 32768.f * blechro[ir][jr] * sincosval.y;
                bufmaskblurcol->b[ir][jr] = 32768.f * blechro[ir][jr] * sincosval.x;
            }
        }


        if (strumask > 0.f && astool) {
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int ir = 0; ir < bfh; ir++) {
                for (int jr = 0; jr < bfw; jr++) {
                    bufmaskblurcol->L[ir][jr] *= (1.f + blendstru[ir][jr]);
                }
            }
        }

        if (lmasklocalcurve && localmaskutili) {
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int ir = 0; ir < bfh; ir++)
                for (int jr = 0; jr < bfw; jr++) {
                    bufmaskblurcol->L[ir][jr] = 0.5f * lmasklocalcurve[2.f * bufmaskblurcol->L[ir][jr]];
                }
        }

        if (shado > 0) {
            ImProcFunctions::shadowsHighlights(bufmaskblurcol, true, 1, 0, shado, 40, sk, 0, 60);
        }

        if (highl > 0) {
            ImProcFunctions::shadowsHighlights(bufmaskblurcol, true, 1, highl, 0, 40, sk, 50, 0);
        }

        int wavelet_level = level_br;

        int minwin = rtengine::min(bfw, bfh);
        int maxlevelspot = 9;

        while ((1 << maxlevelspot) >= (minwin * sk) && maxlevelspot  > 1) {
            --maxlevelspot ;
        }

        wavelet_level = rtengine::min(wavelet_level, maxlevelspot);
        int maxlvl = wavelet_level;
//        float contrast = 0.f;
        bool wavcurvemask = false;

        if (loclmasCurvecolwav && lmasutilicolwav) {
            for (int i = 0; i < 500; i++) {
                if (loclmasCurvecolwav[i] != 0.5f) {
                    wavcurvemask = true;
                    break;
                }
            }
        }

        if (wavcurvemask) {
#ifdef _OPENMP
            const int numThreads = omp_get_max_threads();
#else
            const int numThreads = 1;

#endif
            wavelet_decomposition *wdspot = new wavelet_decomposition(bufmaskblurcol->L[0], bfw, bfh, maxlvl, 1, sk, numThreads, lp.daubLen);
            if (wdspot->memory_allocation_failed()) {
                return;
            }
            float mean[10];
            float meanN[10];
            float sigma[10];
            float sigmaN[10];
            float MaxP[10];
            float MaxN[10];
            Evaluate2(*wdspot, mean, meanN, sigma, sigmaN, MaxP, MaxN, numThreads);
            float alow = 1.f;
            float blow = 0.f;
            if (level_hl != level_bl) {
                alow = 1.f / (level_hl - level_bl);
                blow = -alow * level_bl;
            }

            float ahigh = 1.f;
            float bhigh = 0.f;

            if (level_hr != level_br) {
                ahigh = 1.f / (level_hr - level_br);
                bhigh =  -ahigh * level_br;
            }
        
            for (int dir = 1; dir < 4; dir++) {
                for (int level = level_bl; level < maxlvl; ++level) {
                    int W_L = wdspot->level_W(level);
                    int H_L = wdspot->level_H(level);
                    float* const* wav_L = wdspot->level_coeffs(level);
                
                    if (MaxP[level] > 0.f && mean[level] != 0.f && sigma[level] != 0.f) {
                        float insigma = 0.666f; //SD
                        float logmax = log(MaxP[level]); //log Max
                        float rapX = (mean[level] + sigma[level]) / MaxP[level]; //rapport between sD / max
                        float inx = log(insigma);
                        float iny = log(rapX);
                        float rap = inx / iny; //koef
                        float asig = 0.166f / (sigma[level]);
                        float bsig = 0.5f - asig * mean[level];
                        float amean = 0.5f / mean[level];
                    
#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                        for (int i = 0; i < W_L * H_L; i++) {
                            if(loclmasCurvecolwav && lmasutilicolwav) {
                                float absciss;
                                float &val = wav_L[dir][i];

                                if (fabsf(val) >= (mean[level] + sigma[level])) { //for max
                                    float valcour = xlogf(fabsf(val));
                                    float valc = valcour - logmax;
                                    float vald = valc * rap;
                                    absciss = xexpf(vald);
                                } else if (fabsf(val) >= mean[level]) {
                                    absciss = asig * fabsf(val) + bsig;
                                } else {
                                    absciss = amean * fabsf(val);
                                }

                                float klev = 1.f;
                                if (level >= level_hl && level <= level_hr) {
                                    klev = 1.f;
                                }

                                if (level_hl != level_bl) {
                                    if (level >= level_bl && level < level_hl) {
                                        klev = alow * level + blow;
                                    }
                                }

                                if (level_hr != level_br) {
                                    if (level > level_hr && level <= level_br) {
                                        klev = ahigh * level + bhigh;
                                    }
                                }
            
                                float kc = klev * (loclmasCurvecolwav[absciss * 500.f] - 0.5f);
                                float amplieffect = kc <= 0.f ? 1.f : 4.f;

                                float kinterm = 1.f + amplieffect * kc;
                                kinterm = kinterm <= 0.f ? 0.01f : kinterm;

                                val *=  kinterm;
                            
                            }
                        }
                    }
                
                }
            }
        
            wdspot->reconstruct(bufmaskblurcol->L[0], 1.f);
            delete wdspot;

        }

        if (lochhhmasCurve && HHmaskcurve) {
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int ir = 0; ir < bfh; ir++)
                for (int jr = 0; jr < bfw; jr++) {
                    float huemah = xatan2f(bufmaskblurcol->b[ir][jr], bufmaskblurcol->a[ir][jr]);
                    float chromah = std::sqrt(SQR(bufmaskblurcol->b[ir][jr]) + SQR(bufmaskblurcol->a[ir][jr]));


                    float hh = Color::huelab_to_huehsv2(huemah);
                    hh += 1.f / 6.f;

                    if (hh > 1.f) {
                        hh -= 1.f;
                    }

                    const float val_HH = float ((0.5f - lochhhmasCurve[500.f * hh]));
                    const float hhro = 1.5f * val_HH;
                    float newhr = 0.f;

                    if (hhro != 0) {
                        newhr = huemah + hhro;//we add radians and other dim between 0 1.. always radians but addition "false"

                        if (newhr > rtengine::RT_PI_F) {
                            newhr -= 2 * rtengine::RT_PI_F;
                        } else if (newhr < -rtengine::RT_PI_F) {
                            newhr += 2 * rtengine::RT_PI_F;
                        }
                    }

                    float2 sincosval = xsincosf(newhr);
                    bufmaskblurcol->a[ir][jr] = clipC(chromah * sincosval.y);
                    bufmaskblurcol->b[ir][jr] = clipC(chromah * sincosval.x);

                }
        }

        if (amountcd > 1.f) { //dynamic range compression for Mask
            FattalToneMappingParams fatParams;
            fatParams.enabled = true;
            fatParams.threshold = 100.f;
            fatParams.amount = amountcd;
            fatParams.anchor = anchorcd;
            int nlev = 1;
            Imagefloat *tmpImagefat = nullptr;
            tmpImagefat = new Imagefloat(bfw, bfh);
            lab2rgb(*bufmaskblurcol, *tmpImagefat, params->icm.workingProfile);
            ToneMapFattal02(tmpImagefat, fatParams, nlev, 0, nullptr, 0, 0, 0);
            rgb2lab(*tmpImagefat, *bufmaskblurcol, params->icm.workingProfile);
            delete tmpImagefat;
        }

        if (delt) {
            const std::unique_ptr<JaggedArray<float>> rdEBuffer(new JaggedArray<float>(bfw, bfh));
            float** rdE = *rdEBuffer;

            deltaEforMask(rdE, bfw, bfh, bufreserv.get(), hueref, chromaref, lumaref, maxdE, mindE, maxdElim, mindElim, iterat, limscope, scope, lp.balance, lp.balanceh);
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int ir = 0; ir < bfh; ir++) {
                for (int jr = 0; jr < bfw; jr++) {
                    const float rdEval = rdE[ir][jr];
                    bufmaskblurcol->L[ir][jr] = bufprov->L[ir][jr] + rdEval * (bufmaskblurcol->L[ir][jr] - bufprov->L[ir][jr]);
                    bufmaskblurcol->a[ir][jr] = bufprov->a[ir][jr] + rdEval * (bufmaskblurcol->a[ir][jr] - bufprov->a[ir][jr]);
                    bufmaskblurcol->b[ir][jr] = bufprov->b[ir][jr] + rdEval * (bufmaskblurcol->b[ir][jr] - bufprov->b[ir][jr]);
                }
            }
        }

        struct grad_params gp;

        if ((indic == 0 && lp.strmaexp != 0.f) || (indic ==12 &&  lp.str_mas != 0.f)) {
            calclocalGradientParams(lp, gp, ystart, xstart, bfw, bfh, indic);
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int ir = 0; ir < bfh; ir++) {
                for (int jr = 0; jr < bfw; jr++) {
                    bufmaskblurcol->L[ir][jr] *= ImProcFunctions::calcGradientFactor(gp, jr, ir);
                }
            }
        }
/*
        if (lap > 0.f) {
            const float *datain = bufmaskblurcol->L[0];
            const std::unique_ptr<float[]> data_tmp(new float[bfh * bfw]);

            if (!pde) {
                ImProcFunctions::discrete_laplacian_threshold(data_tmp.get(), datain, bfw, bfh, 200.f * lap);
            } else {
                ImProcFunctions::retinex_pde(datain, data_tmp.get(), bfw, bfh, 12.f * lap, 1.f, nullptr, 0, 0, 1);
            }

#ifdef _OPENMP
            #pragma omp parallel for if (multiThread)
#endif
            for (int y = 0; y < bfh; y++) {
                for (int x = 0; x < bfw; x++) {
                    bufmaskblurcol->L[y][x] = data_tmp[y * bfw + x];
                }
            }
        }
        */
    }

    const float radiusb = 1.f / sk;

    if (deltaE || modmask || enaMask || showmaske) {
#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            gaussianBlur(bufmaskblurcol->L, bufmaskblurcol->L, bfw, bfh, radiusb);
            gaussianBlur(bufmaskblurcol->a, bufmaskblurcol->a, bfw, bfh, 1.f + (0.5f * rad) / sk);
            gaussianBlur(bufmaskblurcol->b, bufmaskblurcol->b, bfw, bfh, 1.f + (0.5f * rad) / sk);
        }

        if (zero || modif || modmask || deltaE || enaMask) {
            originalmaskcol->CopyFrom(bufcolorig, multiThread);
            blendmask(lp, xstart, ystart, cx, cy, bfw, bfh, bufcolorig, original, bufmaskblurcol, originalmaskcol, blendm, blendmab, inv);
        }
    }
}

void ImProcFunctions::InverseSharp_Local(float **loctemp, const float hueref, const float lumaref, const float chromaref, local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy, int sk)
{
//local sharp
    //  BENCHFUN
    const float ach = lp.trans / 100.f;
    const int GW = transformed->W;
    const int GH = transformed->H;
    const float refa = chromaref * cos(hueref) * 327.68f;
    const float refb = chromaref * sin(hueref) * 327.68f;
    const float refL = lumaref * 327.68f;
    //balance deltaE
    const float kL = lp.balance / SQR(327.68f);
    const float kab = balancedeltaE(lp.balance) / SQR(327.68f);
    const float kH = lp.balanceh;
    const float kch = balancedeltaE(kH);
    const bool sharshow = lp.showmasksharmet == 1;
    const bool previewshar = lp.showmasksharmet == 2;

    if (lp.colorde == 0) {
        lp.colorde = -1;//to avoid black
    }

    float ampli = 1.0 + std::fabs(lp.colorde);
    ampli = 2.f + 0.5f * (ampli - 2.f);

    constexpr float aadark = -1.f;
    constexpr float bbdark = 5000.f;

    const std::unique_ptr<LabImage> origblur(new LabImage(GW, GH));

    float radius = 3.f / sk;
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        gaussianBlur(original->L, origblur->L, GW, GH, radius);
        gaussianBlur(original->a, origblur->a, GW, GH, radius);
        gaussianBlur(original->b, origblur->b, GW, GH, radius);
    }

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        const float mindE = 2.f + MINSCOPE * lp.senssha * lp.thr;
        const float maxdE = 5.f + MAXSCOPE * lp.senssha * (1 + 0.1f * lp.thr);
        const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
        const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif
        for (int y = 0; y < transformed->H; y++) {
            int loy = cy + y;

            for (int x = 0; x < transformed->W; x++) {
                int lox = cx + x;
                int zone;
                float localFactor = 1.f;

                if (lp.shapmet == 0) {
                    calcTransition(lox, loy, ach, lp, zone, localFactor);
                } else /*if (lp.shapmet == 1)*/ {
                    calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                }

                const float abdelta2 = SQR(refa - origblur->a[y][x]) + SQR(refb - origblur->b[y][x]);
                const float chrodelta2 = SQR(std::sqrt(SQR(origblur->a[y][x]) + SQR(origblur->b[y][x])) - (chromaref * 327.68f));
                const float huedelta2 = abdelta2 - chrodelta2;
                const float dE = std::sqrt(kab * (kch * chrodelta2 + kH * huedelta2) + kL * SQR(refL - origblur->L[y][x]));

                const float reducdE = calcreducdE(dE, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, lp.senssha);

                switch (zone) {
                    case 0: { // outside selection and outside transition zone => full effect, no transition
                        const float difL = loctemp[y][x] - original->L[y][x];
                        transformed->L[y][x] = CLIP(original->L[y][x] + difL * reducdE);

                        if (sharshow) {
                            transformed->a[y][x] = 0.f;
                            transformed->b[y][x] = ampli * 5.f * difL * reducdE;
                        } else if (previewshar) {
                            float difbdisp = reducdE * 10000.f * lp.colorde;

                            if (transformed->L[y][x] < bbdark) { //enhance dark luminance as user can see!
                                float dark = transformed->L[y][x];
                                transformed->L[y][x] = dark * aadark + bbdark;
                            }

                            if (lp.colorde <= 0) {
                                transformed->a[y][x] = 0.f;
                                transformed->b[y][x] = difbdisp;
                            } else {
                                transformed->a[y][x] = -difbdisp;
                                transformed->b[y][x] = 0.f;
                            }

                        }

                        break;
                    }

                    case 1: { // inside transition zone
                        const float difL = (loctemp[y][x] - original->L[y][x]) * (1.f - localFactor);
                        transformed->L[y][x] = CLIP(original->L[y][x] + difL * reducdE);

                        if (sharshow) {
                            transformed->a[y][x] = 0.f;
                            transformed->b[y][x] = ampli * 5.f * difL * reducdE;
                        } else if (previewshar || lp.prevdE) {
                            const float difbdisp = reducdE * 10000.f * lp.colorde;

                            if (transformed->L[y][x] < bbdark) { //enhance dark luminance as user can see!
                                const float dark = transformed->L[y][x];
                                transformed->L[y][x] = dark * aadark + bbdark;
                            }

                            if (lp.colorde <= 0) {
                                transformed->a[y][x] = 0.f;
                                transformed->b[y][x] = difbdisp;
                            } else {
                                transformed->a[y][x] = -difbdisp;
                                transformed->b[y][x] = 0.f;
                            }
                        }
                        break;
                    }

                    case 2: { // inside selection => no effect, keep original values
                        transformed->L[y][x] = original->L[y][x];
                    }
                }
            }
        }
    }
}

void ImProcFunctions::Sharp_Local(int call, float **loctemp, int senstype, const float hueref, const float chromaref, const float lumaref, local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy, int sk)
{
    //BENCHFUN
    const float ach = lp.trans / 100.f;
    const float varsens = senstype == 1 ? lp.senslc : lp.senssha;
    const bool sharshow = (lp.showmasksharmet == 1);
    const bool previewshar = (lp.showmasksharmet == 2);

    //balance deltaE
    const float kL = lp.balance / SQR(327.68f);
    const float kab = balancedeltaE(lp.balance) / SQR(327.68f);
    const float kH = lp.balanceh;
    const float kch = balancedeltaE(kH);

    if (lp.colorde == 0) {
        lp.colorde = -1;//to avoid black
    }

    float ampli = 1.0 + std::fabs(lp.colorde);
    ampli = 2.f + 0.5f * (ampli - 2.f);

    float darklim = 5000.f;
    float aadark = -1.f;
    float bbdark = darklim;

    const int GW = transformed->W;
    const int GH = transformed->H;

    const std::unique_ptr<LabImage> origblur(new LabImage(GW, GH));
    const float refa = chromaref * cos(hueref) * 327.68f;
    const float refb = chromaref * sin(hueref) * 327.68f;
    const float refL = lumaref * 327.68f;
    const float radius = 3.f / sk;

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        gaussianBlur(original->L, origblur->L, GW, GH, radius);
        gaussianBlur(original->a, origblur->a, GW, GH, radius);
        gaussianBlur(original->b, origblur->b, GW, GH, radius);
    }

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        const int begy = int (lp.yc - lp.lyT);
        const int begx = int (lp.xc - lp.lxL);
        const float mindE = 2.f + MINSCOPE * varsens * lp.thr;
        const float maxdE = 5.f + MAXSCOPE * varsens * (1 + 0.1f * lp.thr);
        const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
        const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif
        for (int y = 0; y < transformed->H; y++) {
            const int loy = cy + y;
            const bool isZone0 = loy > lp.yc + lp.ly || loy < lp.yc - lp.lyT; // whole line is zone 0 => we can skip a lot of processing

            if (isZone0) { // outside selection and outside transition zone => no effect, keep original values
                continue;
            }

            for (int x = 0; x < transformed->W; x++) {
                const int lox = cx + x;
                int zone;
                float localFactor = 1.f;

                if (lp.shapmet == 0) {
                    calcTransition(lox, loy, ach, lp, zone, localFactor);
                } else /*if (lp.shapmet == 1)*/ {
                    calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                }

                if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                    continue;
                }

                //deltaE
                const float abdelta2 = SQR(refa - origblur->a[y][x]) + SQR(refb - origblur->b[y][x]);
                const float chrodelta2 = SQR(std::sqrt(SQR(origblur->a[y][x]) + SQR(origblur->b[y][x])) - (chromaref * 327.68f));
                const float huedelta2 = abdelta2 - chrodelta2;
                const float dE = std::sqrt(kab * (kch * chrodelta2 + kH * huedelta2) + kL * SQR(refL - origblur->L[y][x]));

                float reducdE = calcreducdE(dE, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, varsens);
                const float reducview = reducdE;
                reducdE *= localFactor;

                float difL;

                if (call == 2) {
                    difL = loctemp[loy - begy][lox - begx] - original->L[y][x];
                } else {
                    difL = loctemp[y][x] - original->L[y][x];
                }

                transformed->L[y][x] = CLIP(original->L[y][x] + difL * reducdE);

                if (sharshow) {
                    transformed->a[y][x] = 0.f;
                    transformed->b[y][x] = ampli * 5.f * difL * reducdE;
                } else if (previewshar || lp.prevdE) {
                    float difbdisp = reducview * 10000.f * lp.colorde;

                    if (transformed->L[y][x] < darklim) { //enhance dark luminance as user can see!
                        transformed->L[y][x] = transformed->L[y][x] * aadark + bbdark;
                    }

                    if (lp.colorde <= 0) {
                        transformed->a[y][x] = 0.f;
                        transformed->b[y][x] = difbdisp;
                    } else {
                        transformed->a[y][x] = -difbdisp;
                        transformed->b[y][x] = 0.f;
                    }
                }
            }
        }
    }
}

void ImProcFunctions::Exclude_Local(float **deltaso, float hueref, float chromaref, float lumaref, float sobelref, float meansobel, const struct local_params & lp, const LabImage * original, LabImage * transformed, const LabImage * rsv, const LabImage * reserv, int cx, int cy, int sk)
{

   // BENCHFUN 
    {
        const float ach = lp.trans / 100.f;
        const float varsens =  lp.sensexclu;
        const float mindE = 2.f + MINSCOPE * varsens * lp.thr;
        const float maxdE = 5.f + MAXSCOPE * varsens * (1 + 0.1f * lp.thr);
        const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
        const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);

        const int GW = transformed->W;
        const int GH = transformed->H;

        const float refa = chromaref * cos(hueref) * 327.68f;
        const float refb = chromaref * sin(hueref) * 327.68f;
        const float refL = lumaref * 327.68f;
        // lumaref *= 327.68f;
        //balance deltaE
        const float kL = lp.balance / SQR(327.68f);
        const float kab = balancedeltaE(lp.balance) / SQR(327.68f);
        const float kH = lp.balanceh;
        const float kch = balancedeltaE(kH);
        //sobel
        sobelref = rtengine::min(sobelref / 100.f, 60.f);

        const bool recip = sobelref <  meansobel && sobelref < lp.stru;

        sobelref = log1p(sobelref);

        const std::unique_ptr<LabImage> origblur(new LabImage(GW, GH));

        const float radius = 3.f / sk;

#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            gaussianBlur(reserv->L, origblur->L, GW, GH, radius);
            gaussianBlur(reserv->a, origblur->a, GW, GH, radius);
            gaussianBlur(reserv->b, origblur->b, GW, GH, radius);


#ifdef _OPENMP
            #pragma omp barrier
            #pragma omp for schedule(dynamic,16)
#endif
            for (int y = 0; y < transformed->H; y++)
            {
                const int loy = cy + y;
                const bool isZone0 = loy > (lp.yc + lp.ly - 1) || loy < lp.yc - lp.lyT; // // -1 fix issue 5554

                if (isZone0) { // outside selection and outside transition zone => no effect, keep original values
                    for (int x = 0; x < transformed->W; x++) {
                        transformed->L[y][x] = original->L[y][x];
                    }

                    continue;
                }

                for (int x = 0; x < transformed->W; x++) {
                    const int lox = cx + x;
                    const bool isZone0x = lox > (lp.xc + lp.lx - 1) || lox < lp.xc - lp.lxL; // -1 fix issue 5554

                    if (isZone0x) { // outside selection and outside transition zone => no effect, keep original values
                        transformed->L[y][x] = original->L[y][x];
                        continue;
                    }

                    const int begx = int (lp.xc - lp.lxL);
                    const int begy = int (lp.yc - lp.lyT);

                    int zone;
                    float localFactor = 1.f;

                    if (lp.shapmet == 0) {
                        calcTransition(lox, loy, ach, lp, zone, localFactor);
                    } else /*if (lp.shapmet == 1)*/ {
                        calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                    }


                    if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                        transformed->L[y][x] = original->L[y][x];
                        continue;
                    }

                    float rs = 0.f;

                    const float csob = xlogf(1.f + rtengine::min(deltaso[loy - begy][lox - begx] / 100.f, 60.f) + 0.001f);

                    if (!recip) {
                        rs = sobelref / csob;
                    } else {
                        rs = csob / sobelref;
                    }

                    float affsob = 1.f;

                    if (lp.struexc > 0.f && rs > 0.f) {
                        const float rsob =  0.002f *  lp.struexc * rs;
                        const float minrs = 1.3f + 0.05f * lp.stru;

                        if (rs < minrs) {
                            affsob = 1.f;
                        } else {
                            affsob = 1.f / pow_F((1.f + rsob), SQR(SQR(rs - minrs)));
                        }
                    }

                    float abdelta2 = SQR(refa - origblur->a[y][x]) + SQR(refb - origblur->b[y][x]);
                    float chrodelta2 = SQR(std::sqrt(SQR(origblur->a[y][x]) + SQR(origblur->b[y][x])) - (chromaref * 327.68f));
                    float huedelta2 = abdelta2 - chrodelta2;
                    const float dE = std::sqrt(kab * (kch * chrodelta2 + kH * huedelta2) + kL * SQR(refL - origblur->L[y][x]));
                    const float rL = origblur->L[y][x];
                    const float reducdE = calcreducdE(dE, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, varsens);

                    if (rL > 32.768f) { //to avoid crash with very low gamut in rare cases ex : L=0.01 a=0.5 b=-0.9
                        if (zone > 0) {

                            const float difL = (rsv->L[loy - begy][lox - begx] - original->L[y][x]) * localFactor;
                            transformed->L[y][x] = CLIP(original->L[y][x] + difL * affsob * reducdE);

                            const float difa = (rsv->a[loy - begy][lox - begx] - original->a[y][x]) * localFactor;
                            transformed->a[y][x] = clipC(original->a[y][x] + difa * affsob * reducdE);

                            const float difb = (rsv->b[loy - begy][lox - begx] - original->b[y][x]) * localFactor;
                            transformed->b[y][x] = clipC(original->b[y][x] + difb * affsob * reducdE);

                        }
                    }
                }
            }
        }
    }
}


            
            

void ImProcFunctions::transit_shapedetect_retinex(int call, int senstype, LabImage * bufexporig, LabImage * bufexpfin, LabImage * bufmask, LabImage * buforigmas, float **buflight, float **bufchro, const float hueref, const float chromaref, const float lumaref, struct local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy, int sk)
{

    //BENCHFUN 
    {
        const int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
        const int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
        const int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
        const int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);


        const float ach = lp.trans / 100.f;
        const float varsens = lp.sensh;

        int GW = transformed->W;
        int GH = transformed->H;

        // const float refa = chromaref * cos(hueref);
        // const float refb = chromaref * sin(hueref);

        const float refa = chromaref * cos(hueref) * 327.68f;
        const float refb = chromaref * sin(hueref) * 327.68f;
        const float refL = lumaref * 327.68f;

        const bool retishow = ((lp.showmaskretimet == 1 || lp.showmaskretimet == 2));
        const bool previewreti = ((lp.showmaskretimet == 4));
        //balance deltaE
        const float kL = lp.balance / SQR(327.68f);
        const float kab = balancedeltaE(lp.balance) / SQR(327.68f);
        const float kH = lp.balanceh;
        const float kch = balancedeltaE(kH);
        if (lp.colorde == 0) {
            lp.colorde = -1;//to avoid black
        }
/*
        float ampli = 1.f + std::fabs(lp.colorde);
        ampli = 2.f + 0.5f * (ampli - 2.f);

        float darklim = 5000.f;
        float aadark = -1.f;
        float bbdark = darklim;
*/
        const bool showmas = lp.showmaskretimet == 3 ;

        const std::unique_ptr<LabImage> origblur(new LabImage(GW, GH));
        const float radius = 3.f / sk;
        const bool usemaskreti = lp.enaretiMask && senstype == 4 && !lp.enaretiMasktmap;
        float strcli = 0.03f * lp.str;

        if (lp.scalereti == 1) {
            strcli = 0.015f * lp.str;
        }

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            gaussianBlur(original->L, origblur->L, GW, GH, radius);
            gaussianBlur(original->a, origblur->a, GW, GH, radius);
            gaussianBlur(original->b, origblur->b, GW, GH, radius);
        }


#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            const float mindE = 2.f + MINSCOPE * varsens * lp.thr;
            const float maxdE = 5.f + MAXSCOPE * varsens * (1 + 0.1f * lp.thr);
            const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
            const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
            float previewint = 0.f; //reducdE * 10000.f * lp.colorde; //settings->previewselection;

#ifdef _OPENMP
            #pragma omp for schedule(dynamic,16)
#endif
            for (int y = ystart; y < yend; y++)
            {
                const int loy = cy + y;

                for (int x = xstart; x < xend; x++) {
                    const int lox = cx + x;
                    int zone;
                    float localFactor = 1.f;

                    if (lp.shapmet == 0) {
                        calcTransition(lox, loy, ach, lp, zone, localFactor);
                    } else /*if (lp.shapmet == 1)*/ {
                        calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                    }


                    if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                        continue;
                    }

                    float rL = origblur->L[y][x] / 327.68f;
                    float dE;
                    float abdelta2 = 0.f;
                    float chrodelta2 = 0.f;
                    float huedelta2 = 0.f;

                    if (!usemaskreti) {
                        abdelta2 = SQR(refa - origblur->a[y][x]) + SQR(refb - origblur->b[y][x]);
                        chrodelta2 = SQR(std::sqrt(SQR(origblur->a[y][x]) + SQR(origblur->b[y][x])) - (chromaref * 327.68f));
                        huedelta2 = abdelta2 - chrodelta2;
                        dE = std::sqrt(kab * (kch * chrodelta2 + kH * huedelta2) + kL * SQR(refL - origblur->L[y][x]));
                    } else {
                        if (call == 2) {
                            abdelta2 = SQR(refa - buforigmas->a[y - ystart][x - xstart]) + SQR(refb - buforigmas->b[y - ystart][x - xstart]);
                            chrodelta2 = SQR(std::sqrt(SQR(buforigmas->a[y - ystart][x - xstart]) + SQR(buforigmas->b[y - ystart][x - xstart])) - (chromaref * 327.68f));
                            huedelta2 = abdelta2 - chrodelta2;
                            dE = std::sqrt(kab * (kch * chrodelta2 + kH * huedelta2) + kL * SQR(refL - buforigmas->L[y - ystart][x - xstart]));

                        } else {
                            abdelta2 = SQR(refa - buforigmas->a[y][x]) + SQR(refb - buforigmas->b[y][x]);
                            chrodelta2 = SQR(std::sqrt(SQR(buforigmas->a[y][x]) + SQR(buforigmas->b[y][x])) - (chromaref * 327.68f));
                            huedelta2 = abdelta2 - chrodelta2;
                            dE = std::sqrt(kab * (kch * chrodelta2 + kH * huedelta2) + kL * SQR(refL - buforigmas->L[y][x]));
                        }
                    }

                    float cli, clc;
                    const float reducdE = calcreducdE(dE, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, varsens) / 100.f;
                    previewint = reducdE * 10000.f * lp.colorde; //settings->previewselection;

                    if (call == 2) {
                        cli = buflight[y - ystart][x - xstart];
                        clc = previewreti ? settings->previewselection * 100.f : bufchro[y - ystart][x - xstart];
                    } else {
                        cli = buflight[y][x];
                      //  clc = previewreti ? settings->previewselection * 100.f : bufchro[y][x];
                        clc = previewreti ? reducdE * 10000.f * lp.colorde: bufchro[y][x];

                    }


                    cli *= reducdE;
                    clc *= reducdE;
                    cli *= (1.f + strcli);

                    if (rL > 0.1f) { //to avoid crash with very low gamut in rare cases ex : L=0.01 a=0.5 b=-0.9
                        if (senstype == 4) {//all except color and light (TODO) and exposure
                            float lightc;

                            if (call == 2) {
                                lightc = bufexporig->L[y - ystart][x - xstart];
                            } else {
                                lightc = bufexporig->L[y][x];
                            }

                            float fli = 1.f + cli;
                            float diflc = lightc * fli - original->L[y][x];
                            diflc *= localFactor;

                            if (!showmas) {
                                transformed->L[y][x] = CLIP(original->L[y][x] + diflc);
                            } else {
                                if (call == 2) {

                                    transformed->L[y][x] =  bufmask->L[y - ystart][x - xstart];
                                } else {
                                    transformed->L[y][x] =  bufmask->L[y][x];
                                }
                            } ;

                            if (retishow) {
                                transformed->L[y][x] = CLIP(12000.f + diflc);
                            }
                        }

                        float fliab = 1.f;
                        float chra, chrb;

                        if (call == 2) {
                            chra = bufexporig->a[y - ystart][x - xstart];
                            chrb = bufexporig->b[y - ystart][x - xstart];
                        } else {
                            chra = bufexporig->a[y][x];
                            chrb = bufexporig->b[y][x];

                        }

                        if (senstype == 5) {
                            fliab = 1.f + clc;
                        }

                        const float difa = (chra * fliab - original->a[y][x]) * localFactor;
                        float difb = (chrb * fliab - original->b[y][x]) * localFactor;

                        transformed->a[y][x] = clipC(original->a[y][x] + difa);
                        transformed->b[y][x] = clipC(original->b[y][x] + difb);

                        if (showmas) {
                            if (call == 2) {
                                transformed->a[y][x] = bufmask->a[y - ystart][x - xstart];
                                transformed->b[y][x] = bufmask->b[y - ystart][x - xstart];
                            } else {
                                transformed->a[y][x] = bufmask->a[y][x];
                                transformed->b[y][x] = bufmask->b[y][x];

                            }
                        }

                        if (retishow) {
                            transformed->a[y][x] = clipC(difa);
                            transformed->b[y][x] = clipC(difb);
                        }

                        if (previewreti || lp.prevdE) {
                            difb = (bufexpfin->b[y][x] - original->b[y][x]) * localFactor;
                            transformed->a[y][x] = 0.f;
                            transformed->b[y][x] = previewint * difb;
                        }
                    }
                }
            }
        }

        if (showmas  || retishow || previewreti)
        {
            return;
        }

    }
}


void ImProcFunctions::transit_shapedetect(int senstype, const LabImage * bufexporig, LabImage * originalmask, float **bufchro, bool HHutili, const float hueref, const float chromaref, const float lumaref, float sobelref, float meansobel, float ** blend2, const struct local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy, int sk)
{

   // BENCHFUN
    const int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
    const int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
    const int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
    const int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
    const int bfw = xend - xstart;
    const int bfh = yend - ystart;
//    printf("h=%f l=%f c=%f s=%f\n", hueref, lumaref, chromaref, sobelref);
//    printf("bfh=%i bfw=%i\n", bfh, bfw);
    const float ach = lp.trans / 100.f;
    float varsens = lp.sensex;

    if (senstype == 6 || senstype == 7)   //cbdl
    {
        varsens =  lp.senscb;
    } else if (senstype == 8)   //TM
    {
        varsens =  lp.senstm;
    } else if (senstype == 10) //local contrast
    {
        varsens =  lp.senslc;
    }

    //sobel //keep in case of, not used
    sobelref /= 100.f;
    meansobel /= 100.f;

    sobelref = rtengine::min(sobelref, 60.f);

    const bool k = !(sobelref < meansobel && sobelref < lp.stru); //does not always work with noisy images

    sobelref = log1p(sobelref);

    const float refa = chromaref * cos(hueref) * 327.68f;
    const float refb = chromaref * sin(hueref) * 327.68f;
    const float refL = lumaref * 327.68f;
    const float previewint = settings->previewselection;

    const bool cbshow = ((lp.showmaskcbmet == 1 || lp.showmaskcbmet == 2) &&  senstype == 6);
    const bool tmshow = ((lp.showmasktmmet == 1 || lp.showmasktmmet == 2) &&  senstype == 8);
    const bool previewcb = ((lp.showmaskcbmet == 4) &&  senstype == 6);
    const bool previewtm = ((lp.showmasktmmet == 4) &&  senstype == 8);

    const std::unique_ptr<LabImage> origblur(new LabImage(bfw, bfh));
    std::unique_ptr<LabImage> origblurmask;

    float radius = 3.f / sk;
    //balance deltaE
    const float kL = lp.balance / SQR(327.68f);
    const float kab = balancedeltaE(lp.balance) / SQR(327.68f);
    const bool usemaskcb = (lp.showmaskcbmet == 2 || lp.enacbMask || lp.showmaskcbmet == 4) && senstype == 6;
    const bool usemasktm = (lp.showmasktmmet == 2 || lp.enatmMask || lp.showmasktmmet == 4) && senstype == 8;
    const bool usemaskall = (usemaskcb  || usemasktm);

    if (usemaskall)
    {
        origblurmask.reset(new LabImage(bfw, bfh));

#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            gaussianBlur(originalmask->L, origblurmask->L, bfw, bfh, radius);
            gaussianBlur(originalmask->a, origblurmask->a, bfw, bfh, radius);
            gaussianBlur(originalmask->b, origblurmask->b, bfw, bfh, radius);
        }
    }
    if (lp.equtm  && senstype == 8) //normalize luminance for Tone mapping , at this place we can use for others senstype!
    {
        float *datain = new float[bfh * bfw];
        float *data = new float[bfh * bfw];

#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        for (int y = ystart; y < yend; y++)
            for (int x = xstart; x < xend; x++) {
                datain[(y - ystart) * bfw + (x - xstart)] = original->L[y][x];
                data[(y - ystart)* bfw + (x - xstart)] = bufexporig->L[y - ystart][x - xstart];
            }

        normalize_mean_dt(data, datain, bfh * bfw, 1.f, 1.f, 0.f, 0.f, 0.f, 0.f);
#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif
        for (int y = ystart; y < yend; y++)
            for (int x = xstart; x < xend; x++) {
                bufexporig->L[y - ystart][x - xstart] = data[(y - ystart) * bfw + x - xstart];
            }

        delete [] datain;
        delete [] data;
    }

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif
        for (int y = 0; y < bfh; y++)
        {
            for (int x = 0; x < bfw; x++) {
                origblur->L[y][x] = original->L[y + ystart][x + xstart];
                origblur->a[y][x] = original->a[y + ystart][x + xstart];
                origblur->b[y][x] = original->b[y + ystart][x + xstart];
            }
        }

        gaussianBlur(origblur->L, origblur->L, bfw, bfh, radius);
        gaussianBlur(origblur->a, origblur->a, bfw, bfh, radius);
        gaussianBlur(origblur->b, origblur->b, bfw, bfh, radius);

    }

    const LabImage *maskptr = usemaskall ? origblurmask.get() : origblur.get();
    const float mindE = 2.f + MINSCOPE * varsens * lp.thr;
    const float maxdE = 5.f + MAXSCOPE * varsens * (1 + 0.1f * lp.thr);
    const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
    const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef __SSE2__
        float atan2Buffer[transformed->W] ALIGNED16;
#endif

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif
        for (int y = ystart; y < yend; y++)
        {
            const int loy = cy + y;

#ifdef __SSE2__

            if (HHutili || senstype == 7) {
                int i = xstart;

                for (; i < xend - 3; i += 4) {
                    vfloat av = LVFU(origblur->a[y - ystart][i - xstart]);
                    vfloat bv = LVFU(origblur->b[y - ystart][i - xstart]);
                    STVFU(atan2Buffer[i], xatan2f(bv, av));
                }

                for (; i < xend; i++) {
                    atan2Buffer[i] = xatan2f(origblur->b[y - ystart][i - xstart], origblur->a[y - ystart][i - xstart]);
                }
            }

#endif

            for (int x = xstart; x < xend; x++) {
                const int lox = cx + x;

                int zone;
                float localFactor = 1.f;

                if (lp.shapmet == 0) {
                    calcTransition(lox, loy, ach, lp, zone, localFactor);
                } else /*if (lp.shapmet == 1)*/ {
                    calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                }


                if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                    continue;
                }

                float rhue = 0;

                if (HHutili || senstype == 7) {
#ifdef __SSE2__
                    rhue = atan2Buffer[x];
#else
                    rhue = xatan2f(origblur->b[y - ystart][x - xstart], origblur->a[y - ystart][x - xstart]);
#endif
                }

                const float rL = origblur->L[y - ystart][x - xstart] / 327.68f;
                float rsob = 0.f;

                if (blend2 && ((senstype == 1 && lp.struexp > 0.f) || ((senstype == 0 || senstype == 100) && lp.struco > 0.f))) {//keep in case of, not used
                    const float csob = xlogf(1.f + rtengine::min(blend2[y - ystart][x - xstart] / 100.f, 60.f) + 0.001f);

                    float rs;

                    if (k) {
                        rs = sobelref / csob;
                    } else {
                        rs = csob / sobelref;
                    }

                    if (rs > 0.f && senstype == 1) {
                        rsob =  1.1f * lp.struexp * rs;
                    } else if (rs > 0.f && (senstype == 0 || senstype == 100)) {
                        rsob =  1.1f * lp.struco * rs;
                    }
                }

                const float dE = rsob + std::sqrt(kab * (SQR(refa - maskptr->a[y - ystart][x - xstart]) + SQR(refb - maskptr->b[y - ystart][x - xstart])) + kL * SQR(refL - maskptr->L[y - ystart][x - xstart]));
                const float clc = (previewcb) ? settings->previewselection * 100.f : bufchro[y - ystart][x - xstart];
                const float reducdE = calcreducdE(dE, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, varsens);
                const float realstrchdE = reducdE * clc;

                if (rL > 0.1f) { //to avoid crash with very low gamut in rare cases ex : L=0.01 a=0.5 b=-0.9
                    if (zone > 0) {
                        float factorx = localFactor;
                        float difL = 0.f;

                        if (senstype == 6 || senstype == 8 || senstype == 10) {
                            difL = (bufexporig->L[y - ystart][x - xstart] - original->L[y][x]) * localFactor * reducdE;
                            transformed->L[y][x] = CLIP(original->L[y][x] + difL);
                        }

                        if (senstype == 7) {
                            float difab = bufexporig->L[y - ystart][x - xstart] - std::sqrt(SQR(original->a[y][x]) + SQR(original->b[y][x]));
                            float2 sincosval = xsincosf(rhue);
                            float difa = difab * sincosval.y;
                            float difb = difab * sincosval.x;
                            difa *= factorx * (100.f + realstrchdE) / 100.f;
                            difb *= factorx * (100.f + realstrchdE) / 100.f;
                            transformed->a[y][x] = clipC(original->a[y][x] + difa);
                            transformed->b[y][x] = clipC(original->b[y][x] + difb);
                        } else {
                            float flia = 1.f;
                            float flib = 1.f;
                            float chra = bufexporig->a[y - ystart][x - xstart];
                            float chrb = bufexporig->b[y - ystart][x - xstart];

                            if (senstype == 3 || senstype == 30 || senstype == 8 || senstype == 6 || senstype == 10) {
                                flia = flib = ((100.f + realstrchdE) / 100.f);
                            }


                            float difa = chra * flia - original->a[y][x];
                            float difb = chrb * flib - original->b[y][x];
                            difa *= factorx;
                            difb *= factorx;

                            transformed->a[y][x] = clipC(original->a[y][x] + difa);
                            transformed->b[y][x] = clipC(original->b[y][x] + difb);


                            if (cbshow || tmshow) {
                                transformed->L[y][x] = CLIP(12000.f + difL);
                                transformed->a[y][x] = clipC(difa);
                                transformed->b[y][x] = clipC(difb);
                            } else if (previewcb  || previewtm || lp.prevdE) {
                                if (std::fabs(difb) < 500.f) {
                                    difb += difL;
                                }

                                transformed->a[y][x] = 0.f;
                                transformed->b[y][x] = previewint * difb;
                            }
                        }
                    }
                }
            }
        }
    }
}

void ImProcFunctions::InverseColorLight_Local(bool tonequ, bool tonecurv, int sp, int senstype,  struct local_params & lp, LabImage * originalmask, const LUTf& lightCurveloc, const LUTf& hltonecurveloc, const LUTf& shtonecurveloc, const LUTf& tonecurveloc, const LUTf& exlocalcurve, const LUTf& cclocalcurve, float adjustr, bool localcutili, const LUTf& lllocalcurve, bool locallutili, LabImage * original, LabImage * transformed, int cx, int cy, const float hueref, const float chromaref, const float lumaref, int sk)
{
    // BENCHFUN
    const float ach = lp.trans / 100.f;
    const float facc = (100.f + lp.chro) / 100.f; //chroma factor transition
    float varsens = lp.sens;

    if (senstype == 0) { //Color and Light
        varsens = lp.sens;
    } else if (senstype == 1) { //exposure
        varsens = lp.sensex;
    } else if (senstype == 2) { //shadows highlight
        varsens = lp.senshs;
    }

    const int GW = transformed->W;
    const int GH = transformed->H;
    const float refa = chromaref * cos(hueref) * 327.68f;
    const float refb = chromaref * sin(hueref) * 327.68f;
    const float refL = lumaref * 327.68f;

    const std::unique_ptr<LabImage> temp(new LabImage(GW, GH));
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
    for (int y = 0; y < transformed->H; y++) {
        for (int x = 0; x < transformed->W; x++) {
            temp->L[y][x] = original->L[y][x];
            temp->a[y][x] = original->a[y][x];
            temp->b[y][x] = original->b[y][x];
        }
    }

    if (senstype == 2) { // Shadows highlight
        if (lp.shmeth == 0) {
            ImProcFunctions::shadowsHighlights(temp.get(), lp.hsena, 1, lp.highlihs, lp.shadowhs, lp.radiushs, sk, lp.hltonalhs, lp.shtonalhs);
        } else if (lp.shmeth == 1) {
            const std::unique_ptr<Imagefloat> tmpImage(new Imagefloat(GW, GH));

            lab2rgb(*temp, *tmpImage, params->icm.workingProfile);
            Glib::ustring prof = params->icm.workingProfile;
            if (tonecurv) { //Tone response curve  : does nothing if gamma=2.4 and slope=12.92 ==> gamma sRGB
                const float gamtone = params->locallab.spots.at(sp).gamSH;
                const float slotone = params->locallab.spots.at(sp).sloSH;
                int ill = 0;
                cmsHTRANSFORM dummy = nullptr;
                workingtrc(tmpImage.get(), tmpImage.get(), GW, GH, -5, prof, 2.4, 12.92310, ill, 0, dummy, true, false, false);
              //  workingtrc(tmpImage.get(), tmpImage.get(), GW, GH, 5, prof, gamtone, slotone, illum, 0, dummy, false, true, true);//to keep if we want improve with illuminant and primaries
                workingtrc(tmpImage.get(), tmpImage.get(), GW, GH, 1, prof, gamtone, slotone, ill, 0, dummy, false, true, true);//be careful no gamut control

            }

            if (tonequ) {
                tone_eq(this, tmpImage.get(), lp, params->icm.workingProfile, sk, multiThread);
            }

            rgb2lab(*tmpImage, *temp, params->icm.workingProfile);
        }

    } else if (senstype == 1) { //exposure
        ImProcFunctions::exlabLocal(lp, 1.f, GH, GW, GW, GH, original, temp.get(), hltonecurveloc, shtonecurveloc, tonecurveloc, hueref, lumaref, chromaref);

        if (exlocalcurve) {
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = 0; y < temp->H; y++) {
                for (int x = 0; x < temp->W; x++) {
                    const float lh = 0.5f * exlocalcurve[2.f * temp->L[y][x]]; // / ((lighn) / 1.9f) / 3.61f; //lh between 0 and 0 50 or more
                    temp->L[y][x] = lh;
                }
            }
        }

        if ((lp.expcomp != 0.f) || (exlocalcurve)) {
            if (lp.shadex > 0) {
                ImProcFunctions::shadowsHighlights(temp.get(), true, 1, 0, lp.shadex, 40, sk, 0, lp.shcomp);
            }
        }

        if (lp.expchroma != 0.f) {
            const float ch = (1.f + 0.02f * lp.expchroma) ;
            float chprosl;

            if (ch <= 1.f) {//convert data curve near values of slider -100 + 100, to be used after to detection shape
                chprosl = 99.f * ch - 99.f;
            } else {
                constexpr float ampli = 70.f;
                chprosl = clipChro(ampli * ch - ampli);  //ampli = 25.f arbitrary empirical coefficient between 5 and 50
            }

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = 0; y < transformed->H; y++) {
                for (int x = 0; x < transformed->W; x++) {
                    const float epsi = original->L[y][x] == 0.f ? 0.001f : 0.f;
                    const float rapexp = temp->L[y][x] / (original->L[y][x] + epsi);
                    temp->a[y][x] *= (1.f + chprosl * rapexp);
                    temp->b[y][x] *= (1.f + chprosl * rapexp);
                }
            }
        }
    } else if (senstype == 0) { //Color and Light curves L C
        if (cclocalcurve && localcutili) { // C=f(C) curve
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = 0; y < transformed->H; y++) {
                for (int x = 0; x < transformed->W; x++) {
                    //same as in "normal"
                    const float chromat = std::sqrt(SQR(original->a[y][x]) + SQR(original->b[y][x]));
                    constexpr float ampli = 25.f;
                    const float ch = (cclocalcurve[chromat * adjustr ])  / ((chromat + 0.00001f) * adjustr); //ch between 0 and 0 50 or more
                    const float chprocu = clipChro(ampli * ch - ampli);  //ampli = 25.f arbitrary empirical coefficient between 5 and 50
                    temp->a[y][x] = original->a[y][x] * (1.f + 0.01f * chprocu);
                    temp->b[y][x] = original->b[y][x] * (1.f + 0.01f * chprocu);

                }
            }
        }

        if (lllocalcurve && locallutili) {
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = 0; y < transformed->H; y++) {
                for (int x = 0; x < transformed->W; x++) {
                    temp->L[y][x] = 0.5f * lllocalcurve[2.f * original->L[y][x]];
                }
            }
        }
    }

    //balance deltaE
    const float kL = lp.balance / SQR(327.68f);
    const float kab = balancedeltaE(lp.balance) / SQR(327.68f);
    const float kH = lp.balanceh;
    const float kch = balancedeltaE(kH);

    const std::unique_ptr<LabImage> origblur(new LabImage(GW, GH));
    std::unique_ptr<LabImage> origblurmask;
    const bool usemaskcol = (lp.enaColorMaskinv) && senstype == 0;
    const bool usemaskexp = (lp.enaExpMaskinv) && senstype == 1;
    const bool usemasksh = (lp.enaSHMaskinv) && senstype == 2;
    const bool usemaskall = (usemaskcol || usemaskexp || usemasksh);

    float radius = 3.f / sk;

    if (usemaskall) {
        origblurmask.reset(new LabImage(GW, GH));

#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            gaussianBlur(originalmask->L, origblurmask->L, GW, GH, radius);
            gaussianBlur(originalmask->a, origblurmask->a, GW, GH, radius);
            gaussianBlur(originalmask->b, origblurmask->b, GW, GH, radius);
        }
    }

    if (senstype == 1) {
        radius = (2.f + 0.2f * lp.blurexp) / sk;
    } else if (senstype == 0) {
        radius = (2.f + 0.2f * lp.blurcol) / sk;
    } else if (senstype == 2) {
        radius = (2.f + 0.2f * lp.blurSH) / sk;
    }

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        gaussianBlur(original->L, origblur->L, GW, GH, radius);
        gaussianBlur(original->a, origblur->a, GW, GH, radius);
        gaussianBlur(original->b, origblur->b, GW, GH, radius);
    }

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        const LabImage *maskptr = usemaskall ? origblurmask.get() : origblur.get();
        const float mindE = 2.f + MINSCOPE * varsens * lp.thr;
        const float maxdE = 5.f + MAXSCOPE * varsens * (1 + 0.1f * lp.thr);
        const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
        const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif
        for (int y = 0; y < transformed->H; y++) {
            const int loy = cy + y;

            for (int x = 0; x < transformed->W; x++) {
                const float rL = origblur->L[y][x] / 327.68f;

                if (std::fabs(origblur->b[y][x]) < 0.01f) {
                    origblur->b[y][x] = 0.01f;
                }

                constexpr float th_r = 0.01f;

                if (rL > th_r) { //to avoid crash with very low gamut in rare cases ex : L=0.01 a=0.5 b=-0.9
                    const int lox = cx + x;
                    int zone;
                    float localFactor = 1.f;

                    if (lp.shapmet == 0) {
                        calcTransition(lox, loy, ach, lp, zone, localFactor);
                    } else /*if (lp.shapmet == 1)*/ {
                        calcTransitionrect(lox, loy, ach, lp, zone, localFactor);//rect not good
                    }

                    //deltaE
                    float reducdE;
                    if (zone != 2) {
                        const float abdelta2 = SQR(refa - maskptr->a[y][x]) + SQR(refb - maskptr->b[y][x]);
                        const float chrodelta2 = SQR(std::sqrt(SQR(maskptr->a[y][x]) + SQR(maskptr->b[y][x])) - (chromaref * 327.68f));
                        const float huedelta2 = abdelta2 - chrodelta2;
                        const float dE = std::sqrt(kab * (kch * chrodelta2 + kH * huedelta2) + kL * SQR(refL - maskptr->L[y][x]));
                        reducdE = calcreducdE(dE, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, varsens);
                    }

                    switch (zone) {
                        case 2: { // outside selection and outside transition zone => no effect, keep original values
                            transformed->L[y][x] = original->L[y][x];
                            transformed->a[y][x] = original->a[y][x];
                            transformed->b[y][x] = original->b[y][x];
                            break;
                        }

                        case 1: { // inside transition zone
                            const float factorx = 1.f - localFactor;

                            if (senstype == 0) {
                                const float epsia = original->a[y][x] == 0.f ? 0.0001f : 0.f;
                                const float epsib = original->b[y][x] == 0.f ? 0.0001f : 0.f;
                                float lumnew = original->L[y][x];
                                const float difL = (temp->L[y][x] - original->L[y][x]) * (reducdE * factorx);
                                const float difa = (temp->a[y][x] - original->a[y][x]) * (reducdE * factorx);
                                const float difb = (temp->b[y][x] - original->b[y][x]) * (reducdE * factorx);
                                const float facCa = 1.f + (difa / (original->a[y][x] + epsia));
                                const float facCb = 1.f + (difb / (original->b[y][x] + epsib));

                                if (lp.sens < 75.f) {
                                    if ((lp.ligh != 0.f || lp.cont != 0)) {
                                        lumnew = calclightinv(lumnew, lp.ligh, lightCurveloc);  //replace L-curve
                                    }

                                    const float fac = (100.f + factorx * lp.chro * reducdE) / 100.f; //chroma factor transition
                                    const float diflc = (lumnew - original->L[y][x]) * (reducdE * factorx);

                                    transformed->L[y][x] = CLIP(1.f * (original->L[y][x] + diflc + difL));
                                    transformed->a[y][x] = clipC(original->a[y][x] * fac * facCa) ;
                                    transformed->b[y][x] = clipC(original->b[y][x] * fac * facCb);
                                } else {
                                    const float fac = (100.f + factorx * lp.chro) / 100.f; //chroma factor transition

                                    if ((lp.ligh != 0.f || lp.cont != 0)) {
                                        lumnew = calclightinv(original->L[y][x], lp.ligh, lightCurveloc);
                                    }

                                    const float diflc = (lumnew - original->L[y][x]) * factorx;
                                    transformed->L[y][x] = CLIP(original->L[y][x] + diflc + difL);
                                    transformed->a[y][x] = clipC(original->a[y][x] * fac * facCa);
                                    transformed->b[y][x] = clipC(original->b[y][x] * fac * facCb);
                                }
                            } else if (senstype == 1 || senstype == 2) {
                                const float diflc = (temp->L[y][x] - original->L[y][x]) * (reducdE * factorx);
                                const float difa = (temp->a[y][x] - original->a[y][x]) * (reducdE * factorx);
                                const float difb = (temp->b[y][x] - original->b[y][x]) * (reducdE * factorx);
                                transformed->L[y][x] = CLIP(original->L[y][x] + diflc);
                                transformed->a[y][x] = clipC(original->a[y][x] + difa) ;
                                transformed->b[y][x] = clipC(original->b[y][x] + difb);
                            }

                            break;
                        }

                        case 0: { // inside selection => full effect, no transition
                            if (senstype == 0) {
                                const float epsia = original->a[y][x] == 0.f ? 0.0001f : 0.f;
                                const float epsib = original->b[y][x] == 0.f ? 0.0001f : 0.f;
                                float lumnew = original->L[y][x];
                                const float difL = (temp->L[y][x] - original->L[y][x]) * reducdE;
                                const float difa = (temp->a[y][x] - original->a[y][x]) * reducdE;
                                const float difb = (temp->b[y][x] - original->b[y][x]) * reducdE;
                                const float facCa = 1.f + difa / (original->a[y][x] + epsia);
                                const float facCb = 1.f + difb / (original->b[y][x] + epsib);

                                if (lp.sens < 75.f) {
                                    if ((lp.ligh != 0.f || lp.cont != 0)) {
                                        lumnew = calclightinv(lumnew, lp.ligh, lightCurveloc);  //replace L-curve
                                    }

                                    const float fac = (100.f + lp.chro * reducdE) / 100.f; //chroma factor transition
                                    const float diflc = (lumnew - original->L[y][x]) * reducdE;

                                    transformed->L[y][x] = CLIP(original->L[y][x] + diflc + difL);
                                    transformed->a[y][x] = clipC(original->a[y][x] * fac * facCa) ;
                                    transformed->b[y][x] = clipC(original->b[y][x] * fac * facCb);
                                } else {
                                    if ((lp.ligh != 0.f || lp.cont != 0)) {
                                        lumnew = calclightinv(original->L[y][x], lp.ligh, lightCurveloc);
                                    }

                                    transformed->L[y][x] = CLIP(lumnew + difL) ;
                                    transformed->a[y][x] = clipC(original->a[y][x] * facc * facCa);
                                    transformed->b[y][x] = clipC(original->b[y][x] * facc * facCb);
                                }
                            } else if (senstype == 1  || senstype == 2) {
                                const float diflc = (temp->L[y][x] - original->L[y][x]) * reducdE;
                                const float difa = (temp->a[y][x] - original->a[y][x]) * reducdE;
                                const float difb = (temp->b[y][x] - original->b[y][x]) * reducdE;
                                transformed->L[y][x] = CLIP(original->L[y][x] + diflc);
                                transformed->a[y][x] = clipC(original->a[y][x] + difa) ;
                                transformed->b[y][x] = clipC(original->b[y][x] + difb);
                            }
                        }
                    }
                }
            }
        }
    }
}

void ImProcFunctions::calc_ref(int sp, LabImage * original, LabImage * transformed, int cx, int cy, int oW, int oH, int sk, double & huerefblur, double & chromarefblur, double & lumarefblur, double & hueref, double & chromaref, double & lumaref, double & sobelref, float & avg, const LocwavCurve & locwavCurveden, bool locwavdenutili)
{
    if (params->locallab.enabled) {
        // always calculate hueref, chromaref, lumaref before others operations
        // use in normal mode for all modules except denoise
        struct local_params lp;
        calcLocalParams(sp, oW, oH, params->locallab, lp, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, locwavCurveden, locwavdenutili);
        int begy = lp.yc - lp.lyT;
        int begx = lp.xc - lp.lxL;
        int yEn = lp.yc + lp.ly;
        int xEn = lp.xc + lp.lx;
        float avg2 = 0.f;
        int nc2 = 0;

        for (int y = 0; y < transformed->H ; y++) //{
            for (int x = 0; x < transformed->W; x++) {
                int lox = cx + x;
                int loy = cy + y;

                if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                    avg2 += original->L[y][x];
                    nc2++;
                }
            }

        avg2 /= 32768.f;
        avg = avg2 / nc2;
// double precision for large summations
        double aveA = 0.;
        double aveB = 0.;
        double aveL = 0.;
        double aveChro = 0.;
        double aveAblur = 0.;
        double aveBblur = 0.;
        double aveLblur = 0.;
        double aveChroblur = 0.;

        double avesobel = 0.;
// int precision for the counters
        int nab = 0;
        int nso = 0;
        int nsb = 0;
// single precision for the result
        float avA, avB, avL;
        int spotSize = 0.88623f * rtengine::max(1,  lp.cir / sk);  //18
        //O.88623 = std::sqrt(PI / 4) ==> square equal to circle
        int spotSise2; // = 0.88623f * max (1,  lp.cir / sk); //18

        // very small region, don't use omp here
        LabImage *sobelL;
        LabImage *deltasobelL;
        LabImage *origsob;
        LabImage *origblur = nullptr;
        LabImage *blurorig = nullptr;

        int spotSi = 1 + 2 * rtengine::max(1,  lp.cir / sk);

        if (spotSi < 5) {
            spotSi = 5;
        }

        spotSise2 = (spotSi - 1) / 2;

        JaggedArray<float> blend3(spotSi, spotSi);

        origsob = new LabImage(spotSi, spotSi);
        sobelL = new LabImage(spotSi, spotSi);
        deltasobelL = new LabImage(spotSi, spotSi);
        bool isdenoise = false;

        if ((lp.noiself > 0.f || lp.noiself0 > 0.f || lp.noiself2 > 0.f || lp.wavcurvedenoi || lp.nlstr > 0 || lp.noiselc > 0.f || lp.noisecf > 0.f || lp.noisecc > 0.f) && lp.denoiena) {
            isdenoise = true;
        }

        if (isdenoise) { 
            origblur = new LabImage(spotSi, spotSi);
            blurorig = new LabImage(spotSi, spotSi);

            for (int y = rtengine::max(cy, (int)(lp.yc - spotSise2)); y < rtengine::min(transformed->H + cy, (int)(lp.yc + spotSise2 + 1)); y++) {
                for (int x = rtengine::max(cx, (int)(lp.xc - spotSise2)); x < rtengine::min(transformed->W + cx, (int)(lp.xc + spotSise2 + 1)); x++) {
                    int yb = rtengine::max(cy, (int)(lp.yc - spotSise2));

                    int xb = rtengine::max(cx, (int)(lp.xc - spotSise2));

                    int z = y - yb;
                    int u = x - xb;
                    origblur->L[z][u] = original->L[y - cy][x - cx];
                    origblur->a[z][u] = original->a[y - cy][x - cx];
                    origblur->b[z][u] = original->b[y - cy][x - cx];

                }
            }

            float radius = 3.f / sk;
            {
                //No omp
                gaussianBlur(origblur->L, blurorig->L, spotSi, spotSi, radius);
                gaussianBlur(origblur->a, blurorig->a, spotSi, spotSi, radius);
                gaussianBlur(origblur->b, blurorig->b, spotSi, spotSi, radius);

            }

            for (int y = 0; y < spotSi; y++) {
                for (int x = 0; x < spotSi; x++) {
                    aveLblur += static_cast<double>(blurorig->L[y][x]);
                    aveAblur += static_cast<double>(blurorig->a[y][x]);
                    aveBblur += static_cast<double>(blurorig->b[y][x]);
                    aveChroblur += static_cast<double>(std::sqrt(SQR(blurorig->b[y - cy][x - cx]) + SQR(blurorig->a[y - cy][x - cx])));
                    nsb++;

                }
            }
        }

        //ref for luma, chroma, hue
        for (int y = rtengine::max(cy, (int)(lp.yc - spotSize)); y < rtengine::min(transformed->H + cy, (int)(lp.yc + spotSize + 1)); y++) {
            for (int x = rtengine::max(cx, (int)(lp.xc - spotSize)); x < rtengine::min(transformed->W + cx, (int)(lp.xc + spotSize + 1)); x++) {
                aveL += static_cast<double>(original->L[y - cy][x - cx]);
                aveA += static_cast<double>(original->a[y - cy][x - cx]);
                aveB += static_cast<double>(original->b[y - cy][x - cx]);
                aveChro += static_cast<double>(std::sqrt(SQR(original->b[y - cy][x - cx]) + SQR(original->a[y - cy][x - cx])));
                nab++;
            }
        }

        //ref for sobel
        for (int y = rtengine::max(cy, (int)(lp.yc - spotSise2)); y < rtengine::min(transformed->H + cy, (int)(lp.yc + spotSise2 + 1)); y++) {
            for (int x = rtengine::max(cx, (int)(lp.xc - spotSise2)); x < rtengine::min(transformed->W + cx, (int)(lp.xc + spotSise2 + 1)); x++) {
                int yb = rtengine::max(cy, (int)(lp.yc - spotSise2));

                int xb = rtengine::max(cx, (int)(lp.xc - spotSise2));

                int z = y - yb;
                int u = x - xb;
                origsob->L[z][u] = original->L[y - cy][x - cx];
                nso++;
            }
        }

        const float radius = 3.f / (sk * 1.4f); //0 to 70 ==> see skip

        SobelCannyLuma(sobelL->L, origsob->L, spotSi, spotSi, radius);
        int nbs = 0;

        for (int y = 0; y < spotSi ; y ++)
            for (int x = 0; x < spotSi ; x ++) {
                avesobel += static_cast<double>(sobelL->L[y][x]);
                nbs++;
            }

        sobelref = avesobel / nbs;

        delete sobelL;

        delete deltasobelL;
        delete origsob;
        aveL = aveL / nab;
        aveA = aveA / nab;
        aveB = aveB / nab;
        aveChro = aveChro / nab;
        aveChro /= 327.68;
        avA = aveA / 327.68;
        avB = aveB / 327.68;
        avL = aveL / 327.68;
        hueref = xatan2f(avB, avA);    //mean hue

        if (isdenoise) {
            aveLblur = aveLblur / nsb;
            aveChroblur = aveChroblur / nsb;
            aveChroblur /= 327.68;
            aveAblur = aveAblur / nsb;
            aveBblur = aveBblur / nsb;
            float avAblur = aveAblur / 327.68;
            float avBblur = aveBblur / 327.68;
            float avLblur = aveLblur / 327.68;
            huerefblur = xatan2f(avBblur, avAblur);
            chromarefblur = aveChroblur;
            lumarefblur = avLblur;
        } else {
            huerefblur = 0.f;
            chromarefblur = 0.f;
            lumarefblur = 0.f;
        }

        chromaref = aveChro;
        lumaref = avL;

        //  printf("Calcref => sp=%i befend=%i huere=%2.1f chromare=%2.1f lumare=%2.1f sobelref=%2.1f\n", sp, befend, hueref, chromaref, lumaref, sobelref / 100.f);

        if (isdenoise) {
            delete origblur;
            delete blurorig;
        }

        lumaref = rtengine::min<float>(lumaref, 95.f); //to avoid crash
    }
}
//doc fftw3 says optimum is with size 2^a * 3^b * 5^c * 7^d * 11^e * 13^f with e+f = 0 or 1
//number for size between 18144 and 1 ==> 18000 pixels cover 99% all sensor
const int fftw_size[] = {18144, 18000, 17920, 17836, 17820, 17640, 17600, 17550, 17500, 17496, 17472, 17325, 17280, 17248, 17199, 17150, 17010, 16896, 16875, 16848, 16807,
                         16800, 16640, 16632, 16500, 16464, 16384, 16380, 16250, 16200, 16170, 16128, 16038, 16000, 15925, 15876, 15840, 15795, 15750, 15680, 15625, 15600, 15552, 15435, 15400,
                         15360, 15309, 15288, 15120, 15092, 15000, 14976, 14850, 14784, 14742, 14700, 14625, 14580, 14560, 14553, 14336, 14406, 14400, 14256, 14175, 14112, 14080, 14040, 14000, 13860,
                         13824, 13750, 13720, 13650, 13608, 13500, 13475, 13440, 13377, 13365, 13312, 13230, 13200, 13125, 13122, 13104, 13000, 12960, 12936, 12800, 12740, 12672, 12636, 12600,
                         12544, 12500, 12480, 12474, 12375, 12348, 12320, 12288, 12285, 12250, 12150, 12096, 12005, 12000, 11907, 11880, 11760, 11700, 11664, 11648, 11550, 11520, 11466, 11375,
                         11340, 11319, 11264, 11250, 11232, 11200, 11088, 11025, 11000, 10976, 10935, 10920, 10800, 10780, 10752, 10692, 10584, 10560, 10530, 10400, 10395,  10368, 10290, 10240,
                         10206, 10192, 10125, 10080, 10000, 9984, 9900, 9604, 9856, 9828, 9800, 9750, 9720, 9702, 9625, 9600, 9555, 9504, 9477, 9450, 9408, 9375, 9360, 9261, 9240,
                         9216, 9100, 9072, 9000, 8960, 8918, 8910, 8820, 8800, 8775, 8750, 8748, 8736, 8640, 8624, 8575, 8505, 8448, 8424, 8400, 8320, 8316, 8250, 8232, 8192, 8190, 8125,
                         8100, 8085, 8064, 8019, 8000, 7938, 7920, 7875, 7840, 7800, 7776, 7700, 7680, 7644, 7560, 7546, 7500, 7488, 7425, 7392, 7371, 7350, 7290, 7280, 7203, 7200, 7168,
                         7128, 7056, 7040, 7020, 7000, 6930, 6912, 6875, 6860, 6825, 6804, 6750, 6720, 6656, 6615, 6600, 6561, 6552, 6500, 6480, 6468, 6400, 6370, 6336, 6318, 6300,
                         6272, 6250, 6240, 6237, 6174, 6160, 6144, 6125, 6075, 6048, 6000, 5940, 5880, 5850, 5832, 5824, 5775, 5760, 5670, 5632, 5625, 5616, 5600, 5544, 5500, 5488,
                         5460, 5400, 5390, 5376, 5346, 5292, 5280, 5265, 5250, 5200, 5184, 5145, 5120, 5103, 5096, 5040, 5000, 4992, 4950, 4928, 4914, 4900, 4875, 4860, 4851,  4802,
                         4800, 4752, 4725, 4704, 4680, 4620, 4608, 4550, 4536, 4500, 4480, 4459, 4455, 4410, 4400, 4375, 4374, 4368, 4320, 4312, 4224, 4212, 4200, 4160, 4158, 4125,
                         4116, 4096, 4095, 4050, 4032, 4000, 3969, 3960, 3920, 3900, 3888, 3850, 3840, 3822, 3780, 3773, 3750, 3744, 3696, 3675, 3645, 3640, 3600, 3584, 3564, 3528,
                         3520, 3510, 3500, 3465, 3456, 3430, 3402, 3375, 3360, 3328, 3300, 3276, 3250, 3240, 3234, 3200, 3185, 3168, 3159, 3150, 3136, 3125, 3120, 3087, 3080, 3072,
                         3024, 3000, 2970, 2940, 2925, 2916, 2912, 2880, 2835, 2816, 2808, 2800, 2772, 2750, 2744, 2730, 2700, 2695, 2688, 2673, 2646, 2640, 2625, 2600, 2592, 2560,
                         2548, 2520, 2500, 2496, 2475, 2464, 2457, 2450, 2430, 2401, 2400, 2376, 2352, 2340, 2310, 2304, 2275, 2268, 2250, 2240, 2205, 2200, 2187, 2184, 2160, 2156,
                         2112, 2106, 2100, 2080, 2079, 2058, 2048, 2025, 2016, 2000, 1980, 1960, 1950, 1944, 1936, 1925, 1920, 1911, 1890, 1875, 1872, 1848, 1820, 1800, 1792, 1782,
                         1764, 1760, 1755, 1750, 1728, 1715, 1701, 1680, 1664, 1650, 1638, 1625, 1620, 1617, 1600, 1584, 1575, 1568, 1560, 1540, 1536, 1512, 1500, 1485, 1470, 1458,
                         1456, 1440, 1408, 1404, 1400, 1386, 1375, 1372, 1365, 1350, 1344, 1323, 1320, 1300, 1296, 1280, 1274, 1260, 1250, 1248, 1232, 1225, 1215, 1200, 1188, 1176,
                         1170, 1155, 1152, 1134, 1125, 1120, 1100, 1092, 1080, 1078, 1056, 1053, 1050, 1040, 1029, 1024, 1008, 1000, 990, 980, 975, 972, 960, 945, 936, 924, 910, 900,
                         896, 891, 882, 880, 875, 864, 840, 832, 825, 819, 810, 800, 792, 784, 780, 770, 768, 756, 750, 735, 729, 728, 720, 704, 702, 700, 693, 686, 675, 672, 660,
                         650, 648, 640, 637, 630, 625, 624, 616, 600, 594, 588, 585, 576, 567, 560, 550, 546, 540, 539, 528, 525, 520, 512, 504, 500, 495, 490, 486, 480, 468, 462, 455,
                         450, 448, 441, 440, 432, 420, 416, 405, 400, 396, 392, 390, 385, 384, 378, 375, 364, 360, 352, 351, 350, 343, 336, 330, 325, 324, 320, 315, 312, 308, 300, 297,
                         294, 288, 280, 275, 273, 270, 264, 260, 256, 252, 250, 245, 243, 240, 234, 231, 225, 224, 220, 216, 210, 208, 200, 198, 196, 195, 192, 189, 182, 180, 176, 175,
                         168, 165, 162, 160, 156, 154, 150, 147, 144, 143, 140, 135, 132, 130, 128, 126, 125, 120, 117, 112, 110, 108, 105, 104, 100, 99, 98, 96, 91, 90, 88, 84, 81,
                         80, 78, 77, 75, 72, 70, 66, 65, 64, 63, 60, 56, 55, 54, 52, 50, 49, 48, 45, 44, 42, 40, 39, 36, 35, 33, 32, 30, 28, 27, 26, 25, 24, 22, 21, 20, 18, 16, 15,
                         14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1
                        };

int N_fftwsize = sizeof(fftw_size) / sizeof(fftw_size[0]);


void optfft(int N_fftwsize, int &bfh, int &bfw, int &bfhr, int &bfwr, struct local_params& lp, int H, int W, int &xstart, int &ystart, int &xend, int &yend, int cx, int cy, int fulima)
{
    int ftsizeH = 1;
    int ftsizeW = 1;
    int deltaw = 150;
    int deltah = 150;
    
    if(W < 4000) {
        deltaw = 80;
    }
    if(H < 4000) {
        deltah = 80;
    }


    for (int ft = 0; ft < N_fftwsize; ft++) { //find best values
        if (fftw_size[ft] <= bfh) {
            ftsizeH = fftw_size[ft];
            break;
        }
    }

    for (int ft = 0; ft < N_fftwsize; ft++) {
        if (fftw_size[ft] <= bfw) {
            ftsizeW = fftw_size[ft];
            break;
        }
    }
    
    if(fulima == 2) {// if full image, the ftsizeH and ftsizeW is a bit larger (about 10 to 200 pixels) than the image dimensions so that it is fully processed (consumes a bit more resources)
        for (int ftfu = 0; ftfu < N_fftwsize; ftfu++) { //find best values
            if (fftw_size[ftfu] <= (H + deltah)) {
                ftsizeH = fftw_size[ftfu];
                break;
            }
        }
        for (int ftfu = 0; ftfu < N_fftwsize; ftfu++) { //find best values
            if (fftw_size[ftfu] <= (W + deltaw)) {
                ftsizeW = fftw_size[ftfu];
                break;
            }
        }
    }

    if (settings->verbose) {
        if(fulima == 2) {
            printf("Full image: ftsizeWF=%i ftsizeH=%i\n", ftsizeW, ftsizeH);

        } else {
            printf("ftsizeW=%i ftsizeH=%i\n", ftsizeW, ftsizeH);
        }
    }


    //optimize with size fftw
    bool reduW = false;
    bool reduH = false;
    bool exec = true;
    if (ystart == 0 && yend < H) {
        lp.ly -= (bfh - ftsizeH);
    } else if (ystart != 0 && yend == H) {
        lp.lyT -= (bfh - ftsizeH);
    } else if (ystart != 0 && yend != H) {
        if (lp.ly <= lp.lyT) {
            lp.lyT -= (bfh - ftsizeH);
        } else {
            lp.ly -= (bfh - ftsizeH);
        }
    } else if (ystart == 0 && yend == H) {
       // bfhr = ftsizeH;
        bfhr = bfh;
        reduH = true;
        exec = false;
    }

    if (xstart == 0 && xend < W) {
        lp.lx -= (bfw - ftsizeW);
    } else if (xstart != 0 && xend == W) {
        lp.lxL -= (bfw - ftsizeW);
    } else if (xstart != 0 && xend != W) {
        if (lp.lx <= lp.lxL) {
            lp.lxL -= (bfw - ftsizeW);
        } else {
            lp.lx -= (bfw - ftsizeW);
        }
    } else if (xstart == 0 && xend == W) {
       // bfwr = ftsizeW;
        bfwr = bfw;
        reduW = true;
        exec = false;
    }
    //new values optimized
    ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
    yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, H);
    xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
    xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, W);
    bfh = bfhr = yend - ystart;
    bfw = bfwr = xend - xstart;

    if (reduH && exec) {
        bfhr = ftsizeH;
    } else {
        bfhr = bfh;
    }

    if (reduW && exec) {
        bfwr = ftsizeW;
    } else {
        bfwr = bfw;
    }
}

void ImProcFunctions::BlurNoise_Local(LabImage *tmp1, LabImage * originalmask, const float hueref, const float chromaref, const float lumaref, local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy, int sk)
{
//local BLUR
    //BENCHFUN

    const int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
    const int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
    const int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
    const int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);

    const float ach = lp.trans / 100.f;
    const int GW = transformed->W;
    const int GH = transformed->H;
    const float refa = chromaref * cos(hueref) * 327.68f;
    const float refb = chromaref * sin(hueref) * 327.68f;
    const float refL = lumaref * 327.68f;
    const bool blshow = lp.showmaskblmet == 1 || lp.showmaskblmet == 2;
    const bool previewbl = lp.showmaskblmet == 4;

    //balance deltaE
    const float kL = lp.balance / SQR(327.68f);
    const float kab = balancedeltaE(lp.balance) / SQR(327.68f);
    const float kH = lp.balanceh;
    const float kch = balancedeltaE(kH);

    if (lp.colorde == 0) {
        lp.colorde = -1;//to avoid black
    }

    const float ampli = 1.5f + 0.5f * std::abs(lp.colorde);

    constexpr float darklim = 5000.f;
    constexpr float aadark = -1.f;

    const bool usemaskbl = lp.showmaskblmet == 2 || lp.enablMask || lp.showmaskblmet == 4;
    const bool usemaskall = usemaskbl;
    const float radius = 3.f / sk;
    std::unique_ptr<LabImage> origblurmask;

    const std::unique_ptr<LabImage> origblur(new LabImage(GW, GH));

    if (usemaskall) {
        origblurmask.reset(new LabImage(GW, GH));

#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            gaussianBlur(originalmask->L, origblurmask->L, GW, GH, radius);
            gaussianBlur(originalmask->a, origblurmask->a, GW, GH, radius);
            gaussianBlur(originalmask->b, origblurmask->b, GW, GH, radius);
        }
    }

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        gaussianBlur(original->L, origblur->L, GW, GH, radius);
        gaussianBlur(original->a, origblur->a, GW, GH, radius);
        gaussianBlur(original->b, origblur->b, GW, GH, radius);
    }

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        const LabImage *maskptr = usemaskall ? origblurmask.get() : origblur.get();
        const float mindE = 4.f + MINSCOPE * lp.sensbn * lp.thr;//best usage ?? with blurnoise
        const float maxdE = 5.f + MAXSCOPE * lp.sensbn * (1 + 0.1f * lp.thr);
        const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
        const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif
        for (int y = ystart; y < yend; y++) {
            const int loy = cy + y;

            for (int x = xstart, lox = cx + x; x < xend; x++, lox++) {
                int zone;
                float localFactor = 1.f;

                if (lp.shapmet == 0) {
                    calcTransition(lox, loy, ach, lp, zone, localFactor);
                } else /*if (lp.shapmet == 1)*/ {
                    calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                }

                if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                    continue;
                }

                const float abdelta2 = SQR(refa - maskptr->a[y][x]) + SQR(refb - maskptr->b[y][x]);
                const float chrodelta2 = SQR(std::sqrt(SQR(maskptr->a[y][x]) + SQR(maskptr->b[y][x])) - chromaref * 327.68f);
                const float huedelta2 = abdelta2 - chrodelta2;
                const float dE = std::sqrt(kab * (kch * chrodelta2 + kH * huedelta2) + kL * SQR(refL - maskptr->L[y][x]));
                const float reducdE = calcreducdE(dE, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, lp.sensbn);

                float difL = (tmp1->L[y - ystart][x - xstart] - original->L[y][x]) * localFactor * reducdE;
                transformed->L[y][x] = CLIP(original->L[y][x] + difL);
                const float difa = (tmp1->a[y - ystart][x - xstart] - original->a[y][x]) * localFactor * reducdE;
                const float difb = (tmp1->b[y - ystart][x - xstart] - original->b[y][x]) * localFactor * reducdE;

                transformed->a[y][x] = clipC(original->a[y][x] + difa);
                transformed->b[y][x] = clipC(original->b[y][x] + difb);

                const float maxdifab = rtengine::max(std::fabs(difa), std::fabs(difb));

                if (blshow && lp.colorde < 0) { //show modifications with use "b"
                    //  (origshow && lp.colorde < 0) { //original Retinex
                    transformed->a[y][x] = 0.f;
                    transformed->b[y][x] = ampli * 8.f * difL * reducdE;
                    transformed->L[y][x] = CLIP(12000.f + 0.5f * ampli * difL);

                } else if (blshow && lp.colorde > 0) {//show modifications without use "b"
                    if (difL < 1000.f) {//if too low to be view use ab
                        difL += 0.5f * maxdifab;
                    }

                    transformed->L[y][x] = CLIP(12000.f + 0.5f * ampli * difL);
                    transformed->a[y][x] = clipC(ampli * difa);
                    transformed->b[y][x] = clipC(ampli * difb);
                } else if (previewbl || lp.prevdE) {//show deltaE
                    const float difbdisp = reducdE * 10000.f * lp.colorde;

                    if (transformed->L[y][x] < darklim) { //enhance dark luminance as user can see!
                        float dark = transformed->L[y][x];
                        transformed->L[y][x] = dark * aadark + darklim;
                    }

                    if (lp.colorde <= 0) {
                        transformed->a[y][x] = 0.f;
                        transformed->b[y][x] = difbdisp;
                    } else {
                        transformed->a[y][x] = -difbdisp;
                        transformed->b[y][x] = 0.f;
                    }
                }
            }
        }
    }
}

void ImProcFunctions::transit_shapedetect2(int sp, float meantm, float stdtm, int call, int senstype, const LabImage * bufexporig, const LabImage * bufexpfin, LabImage * originalmask, const float hueref, const float chromaref, const float lumaref, float sobelref, float meansobel, float ** blend2, struct local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy, int sk)
{
    //initialize coordinates
    int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
    int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
    int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
    int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
    int bfw = xend - xstart;
    int bfh = yend - ystart;

    int bfhr = bfh;
    int bfwr = bfw;
    if (lp.blurcolmask >= 0.25f && lp.fftColorMask && call == 2) {
        optfft(N_fftwsize, bfh, bfw, bfhr, bfwr, lp, original->H, original->W, xstart, ystart, xend, yend, cx, cy, lp.fullim);
    }

    bfh = bfhr;
    bfw = bfwr;

    //initialize scope
    float varsens = lp.sensex;//exposure

    if (senstype == 0) { //Color and light
        varsens =  lp.sens;
    } else if (senstype == 2) { //vibrance
        varsens =  lp.sensv;
    } else if (senstype == 9) { //shadowshighlight
        varsens =  lp.senshs;
    } else if (senstype == 3) { //softlight
        varsens =  lp.senssf;
    } else if (senstype == 30) { //dehaze
        varsens =  lp.sensh;
    } else if (senstype == 8) { //TM
        varsens =  lp.senstm;
    } else if (senstype == 10) { //local contrast
        varsens =  lp.senslc;
    } else if (senstype == 11) { //encoding log
        varsens = lp.sensilog;
    } else if (senstype == 20) { //common mask
        varsens = lp.sensimas;
     } else if (senstype == 31) { //ciecam
        varsens = lp.sensicie;
   }
    bool delt = lp.deltaem;
    //sobel
    sobelref /= 100.f;
    meansobel /= 100.f;

    sobelref = rtengine::min(sobelref, 60.f);

    const bool k = !(sobelref < meansobel && sobelref < lp.stru); //does not always work with noisy images

    sobelref = log1p(sobelref);

    //references Spot
    const float refa = chromaref * cos(hueref) * 327.68f;
    const float refb = chromaref * sin(hueref) * 327.68f;
    const float refL = lumaref * 327.68f;

    //to preview modifications, scope, mask
    const bool expshow = ((lp.showmaskexpmet == 1 || lp.showmaskexpmet == 2) &&  senstype == 1);
    const bool vibshow = ((lp.showmaskvibmet == 1 || lp.showmaskvibmet == 2) &&  senstype == 2);
    const bool colshow = ((lp.showmaskcolmet == 1 || lp.showmaskcolmet == 2) &&  senstype == 0);
    const bool SHshow = ((lp.showmaskSHmet == 1 || lp.showmaskSHmet == 2) &&  senstype == 9);
    const bool tmshow = ((lp.showmasktmmet == 1 || lp.showmasktmmet == 2) &&  senstype == 8);
    const bool lcshow = ((lp.showmasklcmet == 1 || lp.showmasklcmet == 2) &&  senstype == 10);
    const bool origshow = ((lp.showmasksoftmet == 5) &&  senstype == 3 && lp.softmet == 1);
    const bool logshow = ((lp.showmasklogmet == 1 || lp.showmasklogmet == 2) &&  senstype == 11);
    const bool cieshow = ((lp.showmaskciemet == 1 || lp.showmaskciemet == 2) &&  senstype == 31);
 
    const bool masshow = ((lp.showmask_met == 1) &&  senstype == 20);

    const bool previewvib = ((lp.showmaskvibmet == 4) &&  senstype == 2);
    const bool previewexp = ((lp.showmaskexpmet == 5) &&  senstype == 1);
    const bool previewcol = ((lp.showmaskcolmet == 5) &&  senstype == 0);
    const bool previewSH = ((lp.showmaskSHmet == 4) &&  senstype == 9);
    const bool previewtm = ((lp.showmasktmmet == 4) &&  senstype == 8);
    const bool previewlc = ((lp.showmasklcmet == 4) &&  senstype == 10);
    const bool previeworig = ((lp.showmasksoftmet == 6) &&  senstype == 3 && lp.softmet == 1);
    const bool previewmas = ((lp.showmask_met == 3) &&  senstype == 20);
    const bool previewlog = ((lp.showmasklogmet == 4) &&  senstype == 11);
    const bool previewcie = ((lp.showmaskciemet == 4) &&  senstype == 31);

    float radius = 3.f / sk;

    if (senstype == 1) {
        radius = (2.f + 0.2f * lp.blurexp) / sk;
    } else if (senstype == 0) {
        radius = (2.f + 0.2f * lp.blurcol) / sk;
    } else if (senstype == 9) {
        radius = (2.f + 0.2f * lp.blurSH) / sk;
    }

    const std::unique_ptr<LabImage> origblur(new LabImage(bfw, bfh));
    std::unique_ptr<LabImage> origblurmask;

    //balance deltaE
    float kL = lp.balance / SQR(327.68f);
    const float kab = balancedeltaE(lp.balance) / SQR(327.68f);
    const float kH = lp.balanceh;
    const float kch = balancedeltaE(kH);

    if (lp.colorde == 0) {
        lp.colorde = -1;//to avoid black
    }

    float ampli = 1.f + std::abs(lp.colorde);
    ampli = 2.f + 0.5f * (ampli - 2.f);

    float darklim = 5000.f;
    float aadark = -1.f;
    float bbdark = darklim;
    bool usemask = true;
    if(originalmask == nullptr) {
        usemask = false;
    }
    const bool usemaskvib = (lp.showmaskvibmet == 2 || lp.enavibMask || lp.showmaskvibmet == 4) && senstype == 2;
    const bool usemaskexp = (lp.showmaskexpmet == 2 || lp.enaExpMask || lp.showmaskexpmet == 5) && senstype == 1;
    const bool usemaskcol = (lp.showmaskcolmet == 2 || lp.enaColorMask || lp.showmaskcolmet == 5) && senstype == 0;
    const bool usemaskSH = (lp.showmaskSHmet == 2 || lp.enaSHMask || lp.showmaskSHmet == 4) && senstype == 9;
    const bool usemasktm = (lp.showmasktmmet == 2 || lp.enatmMask || lp.showmasktmmet == 4) && senstype == 8;
    const bool usemasklc = (lp.showmasklcmet == 2 || lp.enalcMask || lp.showmasklcmet == 4) && senstype == 10;
    const bool usemaskmas = (lp.showmask_met == 1 || lp.ena_Mask || lp.showmask_met == 3) && senstype == 20;
    const bool usemasklog = (lp.showmasklogmet == 2 || lp.enaLMask || lp.showmasklogmet == 4) && senstype == 11;
    const bool usemaskcie = (lp.showmaskciemet == 2 || lp.enacieMask || lp.showmaskciemet == 4) && senstype == 31;
    const bool usemaskall = usemask && (usemaskexp || usemaskvib || usemaskcol || usemaskSH || usemasktm || usemasklc || usemasklog || usemaskcie || usemaskmas);
    //blur a little mask
    if (usemaskall) {
        origblurmask.reset(new LabImage(bfw, bfh));

#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            gaussianBlur(originalmask->L, origblurmask->L, bfw, bfh, radius);
            gaussianBlur(originalmask->a, origblurmask->a, bfw, bfh, radius);
            gaussianBlur(originalmask->b, origblurmask->b, bfw, bfh, radius);
        }
    }



    if (lp.equtm && senstype == 8) { //normalize luminance for Tone mapping , at this place we can use for others senstype!
        float *datain = new float[bfh * bfw];
        float *data = new float[bfh * bfw];

#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int y = ystart; y < yend; y++)
            for (int x = xstart; x < xend; x++) {
                datain[(y - ystart) * bfw + (x - xstart)] = original->L[y][x];
                data[(y - ystart)* bfw + (x - xstart)] = bufexpfin->L[y - ystart][x - xstart];
            }
        if(call == 3 || call == 2) {//improccoordinator and simpleprocess
            normalize_mean_dt(data, datain, bfw * bfh, 1.f, 1.f, 0.f, 0.f, 0.f, 0.f);
        } else if(call == 1) {//dcrop
            float ma = meantm;
            float sa = stdtm;
            float ma2 =  (float) params->locallab.spots.at(sp).noiselumc;
            float sa2 = (float) params->locallab.spots.at(sp).softradiustm;
            //printf("ma=%f sa=%f ma2=%f sa2=%f\n", (double) ma, (double) sa, (double) ma2, (double) sa2);
            //use normalize with mean and stdv
            normalize_mean_dt(data, datain, bfw * bfh, 1.f, 1.f, ma, sa, ma2, sa2);
        }
       
        
        
        
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int y = ystart; y < yend; y++)
            for (int x = xstart; x < xend; x++) {
                bufexpfin->L[y - ystart][x - xstart] = data[(y - ystart) * bfw + x - xstart];
            }

        delete [] datain;
        delete [] data;
    }


    if (senstype == 8) {//strength Tone mapping
        const float repart = 1.0f - 0.01f * lp.repartm;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if(multiThread)
#endif

        for (int y = ystart; y < yend; y++){
            for (int x = xstart; x < xend; x++) {
                bufexpfin->L[y - ystart][x - xstart]= intp(repart, original->L[y][x], bufexpfin->L[y - ystart][x - xstart]);
                bufexpfin->a[y - ystart][x - xstart]= intp(repart, original->a[y][x], bufexpfin->a[y - ystart][x - xstart]);
                bufexpfin->b[y - ystart][x - xstart]= intp(repart, original->b[y][x], bufexpfin->b[y - ystart][x - xstart]);
            }
        }
    }


#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif
        for (int y = 0; y < bfh; y++) {
            for (int x = 0; x < bfw; x++) {
                origblur->L[y][x] = original->L[y + ystart][x + xstart];
                origblur->a[y][x] = original->a[y + ystart][x + xstart];
                origblur->b[y][x] = original->b[y + ystart][x + xstart];
            }
        }

        gaussianBlur(origblur->L, origblur->L, bfw, bfh, radius);
        gaussianBlur(origblur->a, origblur->a, bfw, bfh, radius);
        gaussianBlur(origblur->b, origblur->b, bfw, bfh, radius);

    }
    
    
    //choice between original and mask
    const LabImage *maskptr = usemaskall ? origblurmask.get() : origblur.get();

    //parameters deltaE
    //increase a bit lp.thr and lp.iterat and kL if HDR only with log encoding and CAM16 Jz
    if(senstype == 11 || senstype == 31) {
        lp.thr *= 1.2f;
        lp.iterat *= 1.2f;
        kL *= 1.2f;
    }
   
    const float mindE = 2.f + MINSCOPE * varsens * lp.thr;
    const float maxdE = 5.f + MAXSCOPE * varsens * (1 + 0.1f * lp.thr);
    const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
    const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef __SSE2__
//        float atan2Buffer[transformed->W] ALIGNED16;//keep in case of
#endif

#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif
        for (int y = 0; y < bfh; y++) {

            const int loy = y + ystart + cy;
#ifdef __SSE2__
            /* //keep in case of
                        int i = 0;

                        for (; i < bfw - 3; i += 4) {
                            vfloat av = LVFU(maskptr->a[y][i]);
                            vfloat bv = LVFU(maskptr->b[y][i]);
                            STVFU(atan2Buffer[i], xatan2f(bv, av));
                        }

                        for (; i < bfw; i++) {
                            atan2Buffer[i] = xatan2f(maskptr->b[y][i], maskptr->a[y][i]);
                        }
            */
#endif

            for (int x = 0; x < bfw; x++) {
                const int lox = x + xstart + cx;
                int zone;
                float localFactor = 1.f;
                const float achm = lp.trans / 100.f;

                //calculate transition
                if (lp.shapmet == 0) {
                    calcTransition(lox, loy, achm, lp, zone, localFactor);
                } else /*if (lp.shapmet == 1)*/ {
                    calcTransitionrect(lox, loy, achm, lp, zone, localFactor);
                }

//                float hueh = 0;
#ifdef __SSE2__
//                hueh = atan2Buffer[x];
#else
//                hueh = xatan2f(maskptr->b[y][x], maskptr->a[y][x]);
#endif

                float rsob = 0.f;

                //calculate additive sobel to deltaE
                if (blend2 && ((senstype == 1 && lp.struexp > 0.f) || ((senstype == 0) && lp.struco > 0.f))) {
                    const float csob = xlogf(1.f + rtengine::min(blend2[y][x] / 100.f, 60.f) + 0.001f);

                    float rs;

                    if (k) {
                        rs = sobelref / csob;
                    } else {
                        rs = csob / sobelref;
                    }

                    if (rs > 0.f && senstype == 1) {
                        rsob =  1.1f * lp.struexp * rs;
                    } else if (rs > 0.f && (senstype == 0)) {
                        rsob =  1.1f * lp.struco * rs;
                    }
                }


                //deltaE
                float abdelta2 = SQR(refa - maskptr->a[y][x]) + SQR(refb - maskptr->b[y][x]);
                float chrodelta2 = SQR(std::sqrt(SQR(maskptr->a[y][x]) + SQR(maskptr->b[y][x])) - (chromaref * 327.68f));
                float huedelta2 = abdelta2 - chrodelta2;

                const float dE = rsob + std::sqrt(kab * (kch * chrodelta2 + kH * huedelta2) + kL * SQR(refL - maskptr->L[y][x]));
                //reduction action with deltaE
                const float reducdE = calcreducdE(dE, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, varsens);

                float cli = (bufexpfin->L[y][x] - bufexporig->L[y][x]);
                float cla = (bufexpfin->a[y][x] - bufexporig->a[y][x]);
                float clb = (bufexpfin->b[y][x] - bufexporig->b[y][x]);

                if (delt) {
                    cli = bufexpfin->L[y][x] - original->L[y + ystart][x + xstart];
                    cla = bufexpfin->a[y][x] - original->a[y + ystart][x + xstart];
                    clb = bufexpfin->b[y][x] - original->b[y + ystart][x + xstart];
                }
                if(lp.blwh) {
                    cla = 0.f;
                    clb = 0.f;
                }

                // const float previewint = settings->previewselection;

                const float realstrdE = reducdE * cli;
                const float realstradE = reducdE * cla;
                const float realstrbdE = reducdE * clb;

                float factorx = localFactor;

                if (zone > 0) {
                    //simplified transformed with deltaE and transition
                    transformed->L[y + ystart][x + xstart] = clipLoc(original->L[y + ystart][x + xstart] + factorx * realstrdE);
                    float diflc = factorx * realstrdE;
                    transformed->a[y + ystart][x + xstart] = clipC(original->a[y + ystart][x + xstart] + factorx * realstradE);
                    const float difa = factorx * realstradE;
                    transformed->b[y + ystart][x + xstart] = clipC(original->b[y + ystart][x + xstart] + factorx * realstrbdE);
                    const float difb = factorx * realstrbdE;
                    float maxdifab = rtengine::max(std::fabs(difa), std::fabs(difb));

                    if ((expshow || vibshow || colshow || SHshow || tmshow || lcshow || logshow || cieshow || origshow || masshow) && lp.colorde < 0) { //show modifications with use "b"
                        //  (origshow && lp.colorde < 0) { //original Retinex
                        transformed->a[y + ystart][x + xstart] = 0.f;
                        transformed->b[y + ystart][x + xstart] = ampli * 8.f * diflc * reducdE;
                        transformed->L[y + ystart][x + xstart] = CLIP(12000.f + 0.5f * ampli * diflc);

                    } else if ((expshow || vibshow || colshow || SHshow || tmshow || lcshow || logshow || cieshow || origshow || masshow) && lp.colorde > 0) {//show modifications without use "b"
                        if (diflc < 1000.f) {//if too low to be view use ab
                            diflc += 0.5f * maxdifab;
                        }

                        transformed->L[y + ystart][x + xstart] = CLIP(12000.f + 0.5f * ampli * diflc);
                        transformed->a[y + ystart][x + xstart] = clipC(ampli * difa);
                        transformed->b[y + ystart][x + xstart] = clipC(ampli * difb);
                    } else if (previewexp || previewvib || previewcol || previewSH || previewtm || previewlc || previewlog || previewcie || previeworig || previewmas || lp.prevdE) {//show deltaE
                        float difbdisp = reducdE * 10000.f * lp.colorde;

                        if (transformed->L[y + ystart][x + xstart] < darklim) { //enhance dark luminance as user can see!
                            float dark = transformed->L[y + ystart][x + xstart];
                            transformed->L[y + ystart][x + xstart] = dark * aadark + bbdark;
                        }

                        if (lp.colorde <= 0) {
                            transformed->a[y + ystart][x + xstart] = 0.f;
                            transformed->b[y + ystart][x + xstart] = difbdisp;
                        } else {
                            transformed->a[y + ystart][x + xstart] = -difbdisp;
                            transformed->b[y + ystart][x + xstart] = 0.f;
                        }
                    }
                }
            }
        }
    }
}




void ImProcFunctions::exposure_pde(float * dataor, float * datain, float * dataout, int bfw, int bfh, float thresh, float mod)
/* Jacques Desmis July 2019
** adapted from Ipol Copyright 2009-2011 IPOL Image Processing On Line http://www.ipol.im/
*/
{

    //BENCHFUN
#ifdef RT_FFTW3F_OMP
    if (multiThread) {
        fftwf_init_threads();
        fftwf_plan_with_nthreads(omp_get_max_threads());
    }

#endif
    float *data_fft, *data_tmp, *data;

    if (NULL == (data_tmp = (float *) fftwf_malloc(sizeof(float) * bfw * bfh))) {
        fprintf(stderr, "allocation error\n");
        abort();
    }

    ImProcFunctions::discrete_laplacian_threshold(data_tmp, datain, bfw, bfh, thresh);

    if (NULL == (data_fft = (float *) fftwf_malloc(sizeof(float) * bfw * bfh))) {
        fprintf(stderr, "allocation error\n");
        abort();
    }

    if (NULL == (data = (float *) fftwf_malloc(sizeof(float) * bfw * bfh))) {
        fprintf(stderr, "allocation error\n");
        abort();
    }

    const auto dct_fw = fftwf_plan_r2r_2d(bfh, bfw, data_tmp, data_fft, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    fftwf_execute(dct_fw);

    fftwf_free(data_tmp);

    /* solve the Poisson PDE in Fourier space */
    /* 1. / (float) (bfw * bfh)) is the DCT normalisation term, see libfftw */
    ImProcFunctions::rex_poisson_dct(data_fft, bfw, bfh, 1. / (double)(bfw * bfh));

    const auto dct_bw = fftwf_plan_r2r_2d(bfh, bfw, data_fft, data, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    fftwf_execute(dct_bw);
    fftwf_destroy_plan(dct_fw);
    fftwf_destroy_plan(dct_bw);
    fftwf_free(data_fft);
    fftwf_cleanup();

#ifdef RT_FFTW3F_OMP
    if (multiThread) {
        fftwf_cleanup_threads();
    }
#endif

    normalize_mean_dt(data, dataor, bfw * bfh, mod, 1.f, 0.f, 0.f, 0.f, 0.f);
    {

#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int y = 0; y < bfh ; y++) {
            for (int x = 0; x < bfw; x++) {
                dataout[y * bfw + x] = clipLoc(data[y * bfw + x]);
            }
        }
    }

    fftwf_free(data);
}

void ImProcFunctions::fftw_convol_blur(float * input, float * output, int bfw, int bfh, float radius, int fftkern, int algo)
{
    /*
        ** Jacques Desmis june 2019 - inspired by Copyright 2013 IPOL Image Processing On Line http://www.ipol.im/
        ** when I read documentation on various FFT blur we found 2 possibilities
        ** 0) kernel gauss is used with "normal" data
        ** 1) kernel gauss is used with FFT
        ** fftkern allows to change 0) or 1) and test  It seems the good solution is with 0, but I keep the code in case of ??

        ** input real data to blur
        ** output real data blurred with radius
        ** bfw bfh width and high area
        ** radius = sigma for kernel
        ** n_x n_y relative width and high for kernel
        ** Gaussian blur is given by G(x,y) = (1/2*PI*sigma) * exp(-(x2 + y2) / 2* sigma2)
        ** its traduction in Fourier transform is G(x,y) =  exp((-sigma)*(PI * x2 + PI * y2)), for some authors it is not sigma but sigma^2..I have tried...huge differences with Gaussianblur
        ** after several test the only result that works very well is with fftkern = 0 and algo = 0, and as there is differences with Gaussianblur, I put an empirical correction in Ipretinex and Iplocalcontrast
        ** you can enabled or disabled this function with rtsettings.fftwsigma in options. By default empirical formula is disabled
        ** in fact no importance....if it is this function (for sigma) or another... we are not in research :)
    */
    //BENCHFUN

#ifdef RT_FFTW3F_OMP
    if (multiThread) {
        fftwf_init_threads();
        fftwf_plan_with_nthreads(omp_get_max_threads());
    }
#endif


    float *out; //for FFT data
    float *kern = nullptr;//for kernel gauss
    float *outkern = nullptr;//for FFT kernel
    fftwf_plan p;
    fftwf_plan pkern;//plan for FFT
    int image_size, image_sizechange;
    float n_x = 1.f;
    float n_y = 1.f;//relative coordinates for kernel Gauss
    float radsig = 1.f;

    out = (float*) fftwf_malloc(sizeof(float) * (bfw * bfh));//allocate real data for FFT

    if (fftkern == 1) { //allocate memory FFT if kernel fft = 1
        // kern = new float[bfw * bfh];
        kern = (float*) fftwf_malloc(sizeof(float) * (bfw * bfh));//allocate real data for FFT
        outkern = (float*) fftwf_malloc(sizeof(float) * (bfw * bfh));//allocate real data for FFT
    }

    /*compute the Fourier transform of the input data*/

    p = fftwf_plan_r2r_2d(bfh, bfw, input, out, FFTW_REDFT10, FFTW_REDFT10,  FFTW_ESTIMATE);//FFT 2 dimensions forward  FFTW_MEASURE FFTW_ESTIMATE

    fftwf_execute(p);
    fftwf_destroy_plan(p);

    /*define the gaussian constants for the convolution kernel*/
    if (algo == 0) {
        n_x = rtengine::RT_PI / (double) bfw; //ipol
        n_y = rtengine::RT_PI / (double) bfh;
    } else if (algo == 1) {
        n_x = 1.f / bfw; //gauss
        n_y = 1.f / bfh;
        radsig = 1.f / (2.f * rtengine::RT_PI_F * radius * radius);//gauss
    }

    n_x = n_x * n_x;
    n_y = n_y * n_y;

    image_size = bfw * bfh;
    image_sizechange = 4 * image_size;

    if (fftkern == 1) { //convolution with FFT kernel
#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif
        for (int j = 0; j < bfh; j++) {
            int index = j * bfw;

            for (int i = 0; i < bfw; i++)
                if (algo == 0) {
                    kern[ i + index] = exp((float)(-radius) * (n_x * i * i + n_y * j * j)); //calculate Gauss kernel Ipol formula
                } else if (algo == 1) {
                    kern[ i + index] = radsig * exp((float)(-(n_x * i * i + n_y * j * j) / (2.f * radius * radius))); //calculate Gauss kernel  with Gauss formula
                }
        }

        /*compute the Fourier transform of the kernel data*/
        pkern = fftwf_plan_r2r_2d(bfh, bfw, kern, outkern, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE); //FFT 2 dimensions forward
        fftwf_execute(pkern);
        fftwf_destroy_plan(pkern);

#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif
        for (int j = 0; j < bfh; j++) {
            int index = j * bfw;

            for (int i = 0; i < bfw; i++) {
                out[i + index] *= outkern[i + index];    //apply Gauss kernel with FFT
            }
        }

        fftwf_free(outkern);
        fftwf_free(kern);

        //   delete [] kern;

    } else if (fftkern == 0) {//without FFT kernel
        if (algo == 0) {
#ifdef _OPENMP
            #pragma omp parallel for if (multiThread)
#endif
            for (int j = 0; j < bfh; j++) {
                int index = j * bfw;

                for (int i = 0; i < bfw; i++) {
                    out[i + index] *= exp((float)(-radius) * (n_x * i * i + n_y * j * j));    //apply Gauss kernel without FFT - some authors says radius*radius but differences with Gaussianblur
                }
            }
        } else if (algo == 1) {
#ifdef _OPENMP
            #pragma omp parallel for if (multiThread)
#endif
            for (int j = 0; j < bfh; j++) {
                int index = j * bfw;

                for (int i = 0; i < bfw; i++) {
                    out[i + index] *= radsig * exp((float)(-(n_x * i * i + n_y * j * j) / (2.f * radius * radius)));    //calculate Gauss kernel  with Gauss formula
                }
            }
        }
    }

    p = fftwf_plan_r2r_2d(bfh, bfw, out, output, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);//FFT 2 dimensions backward
    fftwf_execute(p);

#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif
    for (int index = 0; index < image_size; index++) { //restore data
        output[index] /= image_sizechange;
    }

    fftwf_destroy_plan(p);
    fftwf_free(out);

#ifdef RT_FFTW3F_OMP
    if (multiThread) {
        fftwf_cleanup_threads();
    }
#endif
}

void ImProcFunctions::fftw_convol_blur2(float **input2, float **output2, int bfw, int bfh, float radius, int fftkern, int algo)
{
    MyMutex::MyLock lock(*fftwMutex);

    float *input = nullptr;

    if (NULL == (input = (float *) fftwf_malloc(sizeof(float) * bfw * bfh))) {
        fprintf(stderr, "allocation error\n");
        abort();
    }

    float *output = nullptr;

    if (NULL == (output = (float *) fftwf_malloc(sizeof(float) * bfw * bfh))) {
        fprintf(stderr, "allocation error\n");
        abort();
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
    for (int y = 0; y < bfh; y++) {
        for (int x = 0; x < bfw; x++) {
            input[y * bfw + x] =  input2[y][x];
        }
    }

    ImProcFunctions::fftw_convol_blur(input, output, bfw, bfh, radius, fftkern, algo);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
    for (int y = 0; y < bfh; y++) {
        for (int x = 0; x < bfw; x++) {
            output2[y][x] = output[y * bfw + x];
        }
    }

    fftwf_free(input);
    fftwf_free(output);
}


void ImProcFunctions::fftw_tile_blur(int GW, int GH, int tilssize, int max_numblox_W, int min_numblox_W, float **tmp1, int numThreads, double radius)
{
    //BENCHFUN
    float epsil = 0.001f / (tilssize * tilssize);
    fftwf_plan plan_forward_blox[2];
    fftwf_plan plan_backward_blox[2];

    array2D<float> tilemask_in(tilssize, tilssize);
    array2D<float> tilemask_out(tilssize, tilssize);

    float *Lbloxtmp  = reinterpret_cast<float*>(fftwf_malloc(max_numblox_W * tilssize * tilssize * sizeof(float)));
    float *fLbloxtmp = reinterpret_cast<float*>(fftwf_malloc(max_numblox_W * tilssize * tilssize * sizeof(float)));

    int nfwd[2] = {tilssize, tilssize};

    //for DCT:
    fftw_r2r_kind fwdkind[2] = {FFTW_REDFT10, FFTW_REDFT10};
    fftw_r2r_kind bwdkind[2] = {FFTW_REDFT01, FFTW_REDFT01};

    // Creating the plans with FFTW_MEASURE instead of FFTW_ESTIMATE speeds up the execute a bit
    plan_forward_blox[0]  = fftwf_plan_many_r2r(2, nfwd, max_numblox_W, Lbloxtmp, nullptr, 1, tilssize * tilssize, fLbloxtmp, nullptr, 1, tilssize * tilssize, fwdkind, FFTW_MEASURE | FFTW_DESTROY_INPUT);
    plan_backward_blox[0] = fftwf_plan_many_r2r(2, nfwd, max_numblox_W, fLbloxtmp, nullptr, 1, tilssize * tilssize, Lbloxtmp, nullptr, 1, tilssize * tilssize, bwdkind, FFTW_MEASURE | FFTW_DESTROY_INPUT);
    plan_forward_blox[1]  = fftwf_plan_many_r2r(2, nfwd, min_numblox_W, Lbloxtmp, nullptr, 1, tilssize * tilssize, fLbloxtmp, nullptr, 1, tilssize * tilssize, fwdkind, FFTW_MEASURE | FFTW_DESTROY_INPUT);
    plan_backward_blox[1] = fftwf_plan_many_r2r(2, nfwd, min_numblox_W, fLbloxtmp, nullptr, 1, tilssize * tilssize, Lbloxtmp, nullptr, 1, tilssize * tilssize, bwdkind, FFTW_MEASURE | FFTW_DESTROY_INPUT);
    fftwf_free(Lbloxtmp);
    fftwf_free(fLbloxtmp);
    const int border = rtengine::max(2, tilssize / 16);

    for (int i = 0; i < tilssize; ++i) {
        float i1 = abs((i > tilssize / 2 ? i - tilssize + 1 : i));
        float vmask = (i1 < border ? SQR(sin((rtengine::RT_PI_F * i1) / (2 * border))) : 1.0f);
        float vmask2 = (i1 < 2 * border ? SQR(sin((rtengine::RT_PI_F * i1) / (2 * border))) : 1.0f);

        for (int j = 0; j < tilssize; ++j) {
            float j1 = abs((j > tilssize / 2 ? j - tilssize + 1 : j));
            tilemask_in[i][j] = (vmask * (j1 < border ? SQR(sin((rtengine::RT_PI_F * j1) / (2 * border))) : 1.0f)) + epsil;
            tilemask_out[i][j] = (vmask2 * (j1 < 2 * border ? SQR(sin((rtengine::RT_PI_F * j1) / (2 * border))) : 1.0f)) + epsil;

        }
    }

    float *LbloxArray[numThreads];
    float *fLbloxArray[numThreads];

    const int numblox_W = ceil((static_cast<float>(GW)) / offset) + 2;
    const int numblox_H = ceil((static_cast<float>(GH)) / offset) + 2;

    array2D<float> Lresult(GW, GH, ARRAY2D_CLEAR_DATA);
    array2D<float> totwt(GW, GH, ARRAY2D_CLEAR_DATA); //weight for combining DCT blocks

    for (int i = 0; i < numThreads; ++i) {
        LbloxArray[i]  = reinterpret_cast<float*>(fftwf_malloc(max_numblox_W * tilssize * tilssize * sizeof(float)));
        fLbloxArray[i] = reinterpret_cast<float*>(fftwf_malloc(max_numblox_W * tilssize * tilssize * sizeof(float)));
    }

#ifdef _OPENMP
    int masterThread = omp_get_thread_num();
#endif
#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef _OPENMP
        int subThread = masterThread * 1 + omp_get_thread_num();
#else
        int subThread = 0;
#endif
        float *Lblox = LbloxArray[subThread];
        float *fLblox = fLbloxArray[subThread];
        float pBuf[GW + tilssize + 2 * offset] ALIGNED16;
#ifdef _OPENMP
        #pragma omp for
#endif
        for (int vblk = 0; vblk < numblox_H; ++vblk) {

            int top = (vblk - 1) * offset;
            float * datarow = pBuf + offset;

            for (int i = 0; i < tilssize; ++i) {
                int row = top + i;
                int rr = row;

                if (row < 0) {
                    rr = rtengine::min(-row, GH - 1);
                } else if (row >= GH) {
                    rr = rtengine::max(0, 2 * GH - 2 - row);
                }

                for (int j = 0; j < GW; ++j) {
                    datarow[j] = (tmp1[rr][j]);
                }

                for (int j = -1 * offset; j < 0; ++j) {
                    datarow[j] = datarow[rtengine::min(-j, GW - 1)];
                }

                for (int j = GW; j < GW + tilssize + offset; ++j) {
                    datarow[j] = datarow[rtengine::max(0, 2 * GW - 2 - j)];
                }//now we have a padded data row

                for (int hblk = 0; hblk < numblox_W; ++hblk) {
                    int left = (hblk - 1) * offset;
                    int indx = (hblk) * tilssize; //index of block in malloc

                    if (top + i >= 0 && top + i < GH) {
                        int j;

                        for (j = 0; j < rtengine::min((-left), tilssize); ++j) {
                            Lblox[(indx + i)*tilssize + j] = tilemask_in[i][j] * datarow[left + j]; // luma data
                        }

                        for (; j < rtengine::min(tilssize, GW - left); ++j) {
                            Lblox[(indx + i)*tilssize + j] = tilemask_in[i][j] * datarow[left + j]; // luma data
                            totwt[top + i][left + j] += tilemask_in[i][j] * tilemask_out[i][j];
                        }

                        for (; j < tilssize; ++j) {
                            Lblox[(indx + i)*tilssize + j] = tilemask_in[i][j] * datarow[left + j]; // luma data
                        }
                    } else {
                        for (int j = 0; j < tilssize; ++j) {
                            Lblox[(indx + i)*tilssize + j] = tilemask_in[i][j] * datarow[left + j]; // luma data
                        }
                    }

                }

            }//end of filling block row

            //fftwf_print_plan (plan_forward_blox);
            if (numblox_W == max_numblox_W) {
                fftwf_execute_r2r(plan_forward_blox[0], Lblox, fLblox);    // DCT an entire row of tiles
            } else {
                fftwf_execute_r2r(plan_forward_blox[1], Lblox, fLblox);    // DCT an entire row of tiles
            }

            const float n_xy = rtengine::SQR(rtengine::RT_PI / tilssize);

            //radius = 30.f;
            for (int hblk = 0; hblk < numblox_W; ++hblk) {
                int blkstart = hblk * tilssize * tilssize;

                for (int j = 0; j < tilssize; j++) {
                    int index = j * tilssize;

                    for (int i = 0; i < tilssize; i++) {
                        fLblox[blkstart + index + i] *= exp((float)(-radius) * (n_xy * rtengine::SQR(i) + n_xy * rtengine::SQR(j)));
                    }
                }
            }//end of horizontal block loop

            //now perform inverse FT of an entire row of blocks
            if (numblox_W == max_numblox_W) {
                fftwf_execute_r2r(plan_backward_blox[0], fLblox, Lblox);    //for DCT
            } else {
                fftwf_execute_r2r(plan_backward_blox[1], fLblox, Lblox);    //for DCT
            }

            int topproc = (vblk - 1) * offset;
            const int lnumblox_W = ceil((static_cast<float>(GW)) / offset);
            const float DCTnorm = 1.0f / (4 * tilssize * tilssize); //for DCT

            int imin = rtengine::max(0, - topproc);
            int bottom = rtengine::min(topproc + tilssize, GH);
            int imax = bottom - topproc;

            for (int i = imin; i < imax; ++i) {
                for (int hblk = 0; hblk < lnumblox_W; ++hblk) {
                    int left = (hblk - 1) * offset;
                    int right  = rtengine::min(left + tilssize, GW);
                    int jmin = rtengine::max(0, -left);
                    int jmax = right - left;
                    int indx = hblk * tilssize;

                    for (int j = jmin; j < jmax; ++j) {
                        Lresult[topproc + i][left + j] += tilemask_out[i][j] * Lblox[(indx + i) * tilssize + j] * DCTnorm; //for DCT
                    }
                }
            }
        }//end of vertical block loop
    }

#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif
    for (int i = 0; i < GH; ++i) {
        for (int j = 0; j < GW; ++j) {
            tmp1[i][j] = Lresult[i][j] / totwt[i][j];
            tmp1[i][j] = clipLoc(tmp1[i][j]);
        }
    }

    for (int i = 0; i < numThreads; ++i) {
        fftwf_free(LbloxArray[i]);
        fftwf_free(fLbloxArray[i]);
    }

    fftwf_destroy_plan(plan_forward_blox[0]);
    fftwf_destroy_plan(plan_backward_blox[0]);
    fftwf_destroy_plan(plan_forward_blox[1]);
    fftwf_destroy_plan(plan_backward_blox[1]);
    fftwf_cleanup();
}

void ImProcFunctions::wavcbd(wavelet_decomposition &wdspot, int level_bl, int maxlvl,
                             const LocwavCurve& locconwavCurve, bool locconwavutili, float sigm, float offs, float chromalev, int sk)
{
    if (locconwavCurve && locconwavutili) {
        float mean[10];
        float meanN[10];
        float sigma[10];
        float sigmaN[10];
        float MaxP[10];
        float MaxN[10];

#ifdef _OPENMP
        const int numThreads = omp_get_max_threads();
#else
        const int numThreads = 1;
#endif
        Evaluate2(wdspot, mean, meanN, sigma, sigmaN, MaxP, MaxN, numThreads);

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) collapse(2) if (multiThread)
#endif
        for (int dir = 1; dir < 4; dir++) {
            for (int level = level_bl; level < maxlvl; ++level) {
                const int W_L = wdspot.level_W(level);
                const int H_L = wdspot.level_H(level);
                float mea[9];

                float* const* wav_L = wdspot.level_coeffs(level);
                //offset
                float rap = offs * mean[level] - 2.f * sigm * sigma[level];

                if (rap > 0.f) {
                    mea[0] = rap;
                } else {
                    mea[0] = mean[level] / 6.f;
                }

                rap =  offs * mean[level] - sigm * sigma[level];

                if (rap > 0.f) {
                    mea[1] = rap;
                } else {
                    mea[1] = mean[level] / 2.f;
                }

                mea[2] = offs * mean[level]; // 50% data
                mea[3] = offs * mean[level] + sigm * sigma[level] / 2.f;
                mea[4] = offs * mean[level] + sigm * sigma[level]; //66%
                mea[5] = offs * mean[level] + sigm * 1.2f * sigma[level];
                mea[6] = offs * mean[level] + sigm * 1.5f * sigma[level]; //
                mea[7] = offs * mean[level] + sigm * 2.f * sigma[level]; //95%
                mea[8] = offs * mean[level] + sigm * 2.5f * sigma[level]; //99%

                float cpMul = 200.f * (locconwavCurve[level * 55.5f] - 0.5f);

                if (cpMul > 0.f) {
                    cpMul *= 3.5f;
                }

                cpMul /= sk;

                for (int i = 0; i < W_L * H_L; i++) {
                    const float WavCL = std::fabs(wav_L[dir][i]);
                    float beta;

                    //reduction amplification: max action between mean / 2 and mean + sigma
                    // arbitrary coefficient, we can add a slider !!
                    if (WavCL < mea[0]) {
                        beta = 0.6f;    //preserve very low contrast (sky...)
                    } else if (WavCL < mea[1]) {
                        beta = 0.8f;
                    } else if (WavCL < mea[2]) {
                        beta = 1.f;    //standard
                    } else if (WavCL < mea[3]) {
                        beta = 1.f;
                    } else if (WavCL < mea[4]) {
                        beta = 0.8f;    //+sigma
                    } else if (WavCL < mea[5]) {
                        beta = 0.6f;
                    } else if (WavCL < mea[6]) {
                        beta = 0.4f;
                    } else if (WavCL < mea[7]) {
                        beta = 0.2f;    // + 2 sigma
                    } else if (WavCL < mea[8]) {
                        beta = 0.1f;
                    } else {
                        beta = 0.0f;
                    }

                    const float alpha = rtengine::max((1024.f + 15.f * cpMul * beta) / 1024.f, 0.02f) ;
                    wav_L[dir][i] *= alpha * chromalev;
                }
            }
        }
    }
}

void ImProcFunctions::Compresslevels(float **Source, int W_L, int H_L, float compression, float detailattenuator, float thres, float mean, float maxp, float meanN, float maxN, float madL)
{
    //J.Desmis 12-2019

    float exponent;

    if (detailattenuator > 0.f && detailattenuator < 0.05f) {
        const float betemp = expf(-(2.f - detailattenuator + 0.693147f)) - 1.f; //0.69315 = log(2)
        exponent = 1.2f * xlogf(-betemp);
        exponent /= 20.f;
    } else if (detailattenuator >= 0.05f && detailattenuator < 0.25f) {
        const float betemp = expf(-(2.f - detailattenuator + 0.693147f)) - 1.f;
        exponent = 1.2f * xlogf(-betemp);
        exponent /= (-75.f * detailattenuator + 23.75f);
    } else if (detailattenuator >= 0.25f) {
        const float betemp = expf(-(2.f - detailattenuator + 0.693147f)) - 1.f;
        exponent = 1.2f * xlogf(-betemp);
        exponent /= (-2.f * detailattenuator + 5.5f);
    } else {
        exponent = (compression - 1.0f) / 20.f;
    }

    float ap = (thres - 1.f) / (maxp - mean);
    float bp = 1.f - ap * mean;
    ap *= exponent;
    bp *= exponent;

    float a0 = (1.33f * thres - 1.f) / (1.f - mean);
    float b0 = 1.f - a0 * mean;
    a0 *= exponent;
    b0 *= exponent;

    float apn = (thres - 1.f) / (maxN - meanN);
    float bpn = 1.f - apn * meanN;
    apn *= -exponent;
    bpn *= exponent;

    float a0n = (1.33f * thres - 1.f) / (1.f - meanN);
    float b0n = 1.f - a0n * meanN;
    a0n *= -exponent;
    b0n *= exponent;

    madL *= 0.05f;
#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef __SSE2__
        const vfloat apv = F2V(ap);
        const vfloat bpv = F2V(bp);
        const vfloat a0v = F2V(a0);
        const vfloat b0v = F2V(b0);
        const vfloat apnv = F2V(apn);
        const vfloat bpnv = F2V(bpn);
        const vfloat a0nv = F2V(a0n);
        const vfloat b0nv = F2V(b0n);
        const vfloat madLv = F2V(madL);
        const vfloat meanv = F2V(mean);
        const vfloat onev = F2V(1.f);
#endif
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,16)
#endif
        for (int y = 0; y < H_L; y++) {
            int x = 0;
#ifdef __SSE2__
            for (; x < W_L - 3; x += 4) {
                vfloat exponev = onev;
                vfloat valv = LVFU(Source[y][x]);
                const vmask mask1v = vmaskf_ge(valv, ZEROV);
                const vmask mask2v = vmaskf_gt(vself(mask1v, valv, -valv), meanv);
                const vfloat av = vself(mask2v, vself(mask1v, apv, apnv), vself(mask1v, a0v, a0nv));
                const vfloat bv = vself(mask2v, vself(mask1v, bpv, bpnv), vself(mask1v, b0v, b0nv));
                exponev += av * valv + bv;
                valv = vself(mask1v, valv, -valv);
                const vfloat multv = vself(mask1v, onev, -onev);
                const vfloat resultv = multv * xexpf(xlogf(valv + madLv) * exponev);
                STVFU(Source[y][x], resultv);
            }
#endif
            for (; x < W_L; x++) {
                float expone = 1.f;

                if (Source[y][x] >= 0.f) {
                    if (Source[y][x] > mean) {
                        expone += ap * Source[y][x] + bp;
                    } else {
                        expone += a0 * Source[y][x] + b0;
                    }

                    Source[y][x] = xexpf(xlogf(Source[y][x] + madL) * expone);
                } else {
                    if (-Source[y][x] > mean) {
                        expone += apn * Source[y][x] + bpn;
                    } else {
                        expone += a0n * Source[y][x] + b0n;
                    }

                    Source[y][x] = -xexpf(xlogf(-Source[y][x] + madL) * expone);
                }
            }
        }
    }
}

void ImProcFunctions::wavlc(wavelet_decomposition& wdspot, int level_bl, int level_hl, int maxlvl, int level_hr, int level_br, float ahigh, float bhigh, float alow, float blow, float sigmalc, float strength, const LocwavCurve & locwavCurve, int numThreads)
{
        float mean[10];
        float meanN[10];
        float sigma[10];
        float sigmaN[10];
        float MaxP[10];
        float MaxN[10];
        
        Evaluate2(wdspot, mean, meanN, sigma, sigmaN, MaxP, MaxN, numThreads);
        for (int dir = 1; dir < 4; dir++) {
            for (int level = level_bl; level < maxlvl; ++level) {
                int W_L = wdspot.level_W(level);
                int H_L = wdspot.level_H(level);
                float klev = 1.f;

                if (level >= level_hl && level <= level_hr) {
                    klev = 1.f;
                }

                if (level_hl != level_bl) {
                    if (level >= level_bl && level < level_hl) {
                        klev = alow * level + blow;
                    }
                }

                if (level_hr != level_br) {
                    if (level > level_hr && level <= level_br) {
                        klev = ahigh * level + bhigh;
                    }
                }
                float* const* wav_L = wdspot.level_coeffs(level);

                if (MaxP[level] > 0.f && mean[level] != 0.f && sigma[level] != 0.f) {
                    constexpr float insigma = 0.666f; //SD
                    const float logmax = log(MaxP[level]); //log Max
                    const float rapX = (mean[level] + sigmalc * sigma[level]) / MaxP[level]; //rapport between sD / max
                    const float inx = log(insigma);
                    const float iny = log(rapX);
                    const float rap = inx / iny; //koef
                    const float asig = 0.166f / (sigma[level] * sigmalc);
                    const float bsig = 0.5f - asig * mean[level];
                    const float amean = 0.5f / mean[level];
                    const float limit1 = mean[level] + sigmalc * sigma[level];
                    const float limit2 = mean[level];
#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic, 16 * W_L) if (multiThread)
#endif
                    for (int i = 0; i < W_L * H_L; i++) {
                        const float val = std::fabs(wav_L[dir][i]);

                        float absciss;
                        if (val >= limit1) { //for max
                            const float valcour = xlogf(val);
                            absciss = xexpf((valcour - logmax) * rap);
                        } else if (val >= limit2) {
                            absciss = asig * val + bsig;
                        } else {
                            absciss = amean * val;
                        }

                        const float kc = klev * (locwavCurve[absciss * 500.f] - 0.5f);
                        const float reduceeffect = kc <= 0.f ? 1.f : strength;

                        float kinterm = 1.f + reduceeffect * kc;
                        kinterm = kinterm <= 0.f ? 0.01f : kinterm;

                        wav_L[dir][i] *= kinterm <= 0.f ? 0.01f : kinterm;
                    }
                }
            }
        }
    
}    

void ImProcFunctions::wavcont(const struct local_params& lp, float ** tmp, wavelet_decomposition& wdspot, int level_bl, int maxlvl,
                              const LocwavCurve & loclevwavCurve, bool loclevwavutili,
                              const LocwavCurve & loccompwavCurve, bool loccompwavutili,
                              const LocwavCurve & loccomprewavCurve, bool loccomprewavutili,
                              float radlevblur, int process, float chromablu, float thres,  float sigmadc, float deltad)
{
    //BENCHFUN
    const int W_L = wdspot.level_W(0);
    const int H_L = wdspot.level_H(0);

#ifdef _OPENMP
    const int numThreads = omp_get_max_threads();
#else
    const int numThreads = 1;
#endif
    float mean[10];
    float meanN[10];
    float sigma[10];
    float sigmaN[10];
    float MaxP[10];
    float MaxN[10];
    Evaluate2(wdspot, mean, meanN, sigma, sigmaN, MaxP, MaxN, numThreads);

    if (process == 1 && loclevwavCurve && loclevwavutili) { //blur
        array2D<float> templevel(W_L, H_L);
        for (int dir = 1; dir < 4; ++dir) {
            for (int level = level_bl; level < maxlvl; ++level) {
                const auto WavL = wdspot.level_coeffs(level)[dir];
                const float effect = lp.sigmabl;
                constexpr float offs = 1.f;
                float mea[10];
                calceffect(level, mean, sigma, mea, effect, offs);
                float lutFactor;
                const float inVals[] = {0.05f, 0.2f, 0.7f, 1.f, 1.f, 0.8f, 0.5f, 0.3f, 0.2f, 0.1f, 0.05f};
                const auto meaLut = buildMeaLut(inVals, mea, lutFactor);

                const float klev = 0.25f * loclevwavCurve[level * 55.5f];
                float* src[H_L];
                for (int i = 0; i < H_L; ++i) {
                    src[i] = &wdspot.level_coeffs(level)[dir][i * W_L];
                }
#ifdef _OPENMP
                #pragma omp parallel if (multiThread)
#endif
                {
                    gaussianBlur(src, templevel, W_L, H_L, radlevblur * klev * chromablu);
                }
#ifdef _OPENMP
                #pragma omp parallel if (multiThread)
#endif
                {
#ifdef __SSE2__
                    const vfloat lutFactorv = F2V(lutFactor);
#endif
#ifdef _OPENMP
                    #pragma omp for
#endif
                    for (int y = 0; y < H_L; y++) {
                        int x = 0;
                        int j = y * W_L;
#ifdef __SSE2__
                        for (; x < W_L - 3; x += 4, j += 4) {
                            const vfloat valv = LVFU(WavL[j]);
                            STVFU(WavL[j], intp((*meaLut)[vabsf(valv) * lutFactorv], LVFU(templevel[y][x]), valv));
                        }
#endif
                        for (; x < W_L; x++, j++) {
                            WavL[j] = intp((*meaLut)[std::fabs(WavL[j]) * lutFactor], templevel[y][x], WavL[j]);
                        }
                    }
                }
            }
        }
    } else if (process == 2 && loccompwavCurve && loccompwavutili) { //Directional contrast
        for (int dir = 1; dir < 4; ++dir) {
            for (int level = level_bl; level < maxlvl; ++level) {
                const auto WavL = wdspot.level_coeffs(level)[dir];
                const float effect = sigmadc;
                constexpr float offs = 1.f;
                float mea[10];
                calceffect(level, mean, sigma, mea, effect, offs);
                float lutFactor;
                const float inVals[] = {0.05f, 0.2f, 0.7f, 1.f, 1.f, 0.8f, 0.7f, 0.5f, 0.3f, 0.2f, 0.1f};
                const auto meaLut = buildMeaLut(inVals, mea, lutFactor);
                const int iteration = deltad;
                const int itplus = 7 + iteration;
                const int itmoins = 7 - iteration;
                const int med = maxlvl / 2;
                int it;

                if (level < med) {
                    it = itmoins;
                } else if (level == med) {
                    it = 7;
                } else {
                    it = itplus;
                }

                const float itf = it;
                const float factor = dir < 3 ? 0.3f : -0.6f;
#ifdef _OPENMP
                #pragma omp parallel if (multiThread)
#endif
                {
#ifdef __SSE2__
                    const vfloat c327d68v = F2V(327.68f);
                    const vfloat factorv = F2V(factor);
                    const vfloat sixv = F2V(6.f);
                    const vfloat zd5v = F2V(0.5f);
                    const vfloat onev = F2V(1.f);
                    const vfloat itfv = F2V(itf);
                    const vfloat lutFactorv = F2V(lutFactor);
#endif
#ifdef _OPENMP
                    #pragma omp for
#endif
                    for (int i = 0; i < H_L; ++i) {
                        int j = 0;
#ifdef __SSE2__
                        for (; j < W_L - 3; j += 4) {
                            const vfloat LL100v = LC2VFU(tmp[i * 2][j * 2]) / c327d68v;
                            const vfloat kbav = factorv * (loccompwavCurve[sixv * LL100v] - zd5v); //k1 between 0 and 0.5    0.5==> 1/6=0.16
                            const vfloat valv = LVFU(WavL[i * W_L + j]);
                            STVFU(WavL[i * W_L + j], valv * pow_F(onev + kbav * (*meaLut)[vabsf(valv) * lutFactorv], itfv));
                        }
#endif
                        for (; j < W_L; ++j) {
                            const float LL100 = tmp[i * 2][j * 2] / 327.68f;
                            const float kba = factor * (loccompwavCurve[6.f * LL100] - 0.5f); //k1 between 0 and 0.5    0.5==> 1/6=0.16
                            WavL[i * W_L + j] *= pow_F(1.f + kba * (*meaLut)[std::fabs(WavL[i * W_L + j]) * lutFactor], itf);
                        }
                    }
                }
            }
        }
    } else if (process == 3 && loccomprewavCurve && loccomprewavutili) { //Dynamic compression wavelet
        float madL[10][3];
        array2D<float> templevel(W_L, H_L);
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) collapse(2) if (multiThread)
#endif
        for (int dir = 1; dir < 4; dir++) {
            for (int level = level_bl; level < maxlvl; ++level) {
                madL[level][dir - 1] = Mad(wdspot.level_coeffs(level)[dir], wdspot.level_W(level) * wdspot.level_H(level)); //evaluate noise by level
            }
        }

        for (int dir = 1; dir < 4; ++dir) {
            for (int level = level_bl; level < maxlvl; ++level) {
                const float effect = lp.sigmadr;
                constexpr float offs = 1.f;
                float mea[10];
                calceffect(level, mean, sigma, mea, effect, offs);
                float lutFactor;
                const float inVals[] = {0.05f, 0.2f, 0.7f, 1.f, 1.f, 0.8f, 0.65f, 0.5f, 0.4f, 0.25f, 0.1f};
                const auto meaLut = buildMeaLut(inVals, mea, lutFactor);
                const auto wav_L = wdspot.level_coeffs(level)[dir];

#ifdef _OPENMP
                #pragma omp parallel for if (multiThread)
#endif
                for (int y = 0; y < H_L; y++) {
                    for (int x = 0; x < W_L; x++) {
                        int j = y * W_L + x;
                        templevel[y][x] = wav_L[j];
                    }
                }

                float klev = (loccomprewavCurve[level * 55.5f] - 0.75f);
                if (klev < 0.f) {
                    klev *= 2.6666f;//compression increase contraste
                } else {
                    klev *= 4.f;//dilatation reduce contraste - detailattenuator
                }
                const float compression = expf(-klev);
                const float detailattenuator = std::max(klev, 0.f);

                Compresslevels(templevel, W_L, H_L, compression, detailattenuator, thres, mean[level], MaxP[level], meanN[level], MaxN[level], madL[level][dir - 1]);
#ifdef _OPENMP
                #pragma omp parallel if (multiThread)
#endif
                {
#ifdef __SSE2__
                    const vfloat lutFactorv = F2V(lutFactor);
#endif
#ifdef _OPENMP
                    #pragma omp for
#endif
                    for (int y = 0; y < H_L; y++) {
                        int x = 0;
                        int j = y * W_L;
#ifdef __SSE2__
                        for (; x < W_L - 3; x += 4, j += 4) {
                            const vfloat valv = LVFU(wav_L[j]);
                            STVFU(wav_L[j], intp((*meaLut)[vabsf(valv) * lutFactorv], LVFU(templevel[y][x]), valv));
                        }
#endif
                        for (; x < W_L; x++, j++) {
                            wav_L[j] = intp((*meaLut)[std::fabs(wav_L[j]) * lutFactor], templevel[y][x], wav_L[j]);
                        }
                    }
                }
            }
        }
    }
}


void ImProcFunctions::wavcontrast4(struct local_params& lp, float ** tmp, float ** tmpa, float ** tmpb, float contrast, float radblur, float radlevblur, int bfw, int bfh, int level_bl, int level_hl, int level_br, int level_hr, int sk, int numThreads,
                                   const LocwavCurve & locwavCurve, bool locwavutili, bool wavcurve, const LocwavCurve& loclevwavCurve, bool loclevwavutili, bool wavcurvelev,
                                   const LocwavCurve & locconwavCurve, bool locconwavutili, bool wavcurvecon,
                                   const LocwavCurve & loccompwavCurve, bool loccompwavutili, bool wavcurvecomp,
                                   const LocwavCurve & loccomprewavCurve, bool loccomprewavutili, bool wavcurvecompre,
                                   const LocwavCurve & locedgwavCurve, bool locedgwavutili,
                                   float sigm, float offs, int & maxlvl, float sigmadc, float deltad, float chromalev, float chromablu, bool blurlc, bool blurena, bool levelena, bool comprena, bool compreena, float compress, float thres)
{
//BENCHFUN
    std::unique_ptr<wavelet_decomposition> wdspot(new wavelet_decomposition(tmp[0], bfw, bfh, maxlvl, 1, sk, numThreads, lp.daubLen));

    //first decomposition for compress dynamic range positive values and other process
    if (wdspot->memory_allocation_failed()) {
        return;
    }

    struct grad_params gpwav;

    maxlvl = wdspot->maxlevel();

    int W_Lm = wdspot->level_W(maxlvl - 1); //I assume all decomposition have same W and H

    int H_Lm = wdspot->level_H(maxlvl - 1);

    if (lp.strwav != 0.f && lp.wavgradl) {
        array2D<float> factorwav(W_Lm, H_Lm);
        calclocalGradientParams(lp, gpwav, 0, 0, W_Lm, H_Lm, 10);
        const float mult = lp.strwav < 0.f ? -1.f : 1.f;
#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < H_Lm; y++) {
            for (int x = 0; x < W_Lm; x++) {
                factorwav[y][x] = mult * (1.f - ImProcFunctions::calcGradientFactor(gpwav, x, y));
            }
        }
        float mean[10];
        float meanN[10];
        float sigma[10];
        float sigmaN[10];
        float MaxP[10];
        float MaxN[10];
        Evaluate2(*wdspot, mean, meanN, sigma, sigmaN, MaxP, MaxN, numThreads);
        float alowg = 1.f;
        float blowg = 0.f;

        if (level_hl != level_bl) {
            alowg = 1.f / (level_hl - level_bl);
            blowg = -alowg * level_bl;
        }

        float ahighg = 1.f;
        float bhighg = 0.f;

        if (level_hr != level_br) {
            ahighg = 1.f / (level_hr - level_br);
            bhighg =  -ahighg * level_br;
        }

        for (int dir = 1; dir < 4; dir++) {
            for (int level = level_bl; level < maxlvl; ++level) {
                if (MaxP[level] > 0.f && mean[level] != 0.f && sigma[level] != 0.f) {
                    const int W_L = wdspot->level_W(level);
                    const int H_L = wdspot->level_H(level);
                    auto wav_L = wdspot->level_coeffs(level)[dir];
                    const float effect = lp.sigmalc2;
                    constexpr float offset = 1.f;
                    float mea[10];
                    calceffect(level, mean, sigma, mea, effect, offset);
                    constexpr float insigma = 0.666f; //SD
                    const float logmax = std::log(MaxP[level]); //log Max
                    const float rapX = (mean[level] + lp.sigmalc2 * sigma[level]) / MaxP[level]; //rapport between sD / max
                    const float inx = std::log(insigma);
                    const float iny = std::log(rapX);
                    const float rap = inx / iny; //koef
                    const float asig = 0.166f / (sigma[level] * lp.sigmalc2);
                    const float bsig = 0.5f - asig * mean[level];
                    const float amean = 0.5f / mean[level];
                    float klev = 1.f;

                    if (level_hl != level_bl) {
                        if (level >= level_bl && level < level_hl) {
                            klev = alowg * level + blowg;
                        }
                    }

                    if (level_hr != level_br) {
                        if (level > level_hr && level <= level_br) {
                            klev = ahighg * level + bhighg;
                        }
                    }
                    klev *= 0.8f;
                    const float threshold = mean[level] + lp.sigmalc2 * sigma[level];
                    float lutFactor;
                    const float inVals[] = {0.05f, 0.2f, 0.7f, 1.f, 1.f, 0.8f, 0.6f, 0.5f, 0.4f, 0.3f, 0.1f};
                    const auto meaLut = buildMeaLut(inVals, mea, lutFactor);

#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic, 16) if (multiThread)
#endif

                    for (int y = 0; y < H_L; y++) {
                        for (int x = 0; x < W_L; x++) {
                            const float WavCL = std::fabs(wav_L[y * W_L + x]);

                            float absciss;
                            if (WavCL >= threshold) { //for max
                                absciss = pow_F(WavCL - logmax, rap);
                            } else if (WavCL >= mean[level]) {
                                absciss = asig * WavCL + bsig;
                            } else {
                                absciss = amean * WavCL;
                            }

                            const float kc = klev * factorwav[y][x] * absciss;
                            const float reduceeffect = kc <= 0.f ? 1.f : 1.5f;

                            float kinterm = 1.f + reduceeffect * kc;
                            kinterm = kinterm <= 0.f ? 0.01f : kinterm;
                            wav_L[y * W_L + x] *= (1.f + (kinterm - 1.f) * (*meaLut)[WavCL * lutFactor]);
                        }
                    }
                }
            }
        }
    }

    int W_Level = wdspot->level_W(0);
    int H_Level = wdspot->level_H(0);
    float *wav_L0 = wdspot->get_coeff0();

    if (radblur > 0.f && blurena) {
        float* src[H_Level];
        for (int i = 0; i < H_Level; ++i) {
            src[i] = &wav_L0[i * W_Level];
        }

#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            gaussianBlur(src, src, W_Level, H_Level, radblur);
        }
    }

    if (compress != 0.f && compreena) {
        const float Compression = expf(-compress);
        const float DetailBoost = std::max(compress, 0.f);

        CompressDR(wav_L0, W_Level, H_Level, Compression, DetailBoost);
    }

    if ((lp.residsha < 0.f || lp.residhi < 0.f)) {
        float tran = 5.f;//transition shadow

        if (lp.residshathr > (100.f - tran)) {
            tran = 100.f - lp.residshathr;
        }
        constexpr float alp = 3.f;
        const float aalp = (1.f - alp) / lp.residshathr;
        const float ath = -lp.residsha / tran;
        const float bth = lp.residsha - ath * lp.residshathr;

        //highlight
        const float tranh = rtengine::min(5.f, lp.residhithr);
        const float athH = lp.residhi / tranh;
        const float bthH = lp.residhi - athH * lp.residhithr;

#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif
        for (int i = 0; i < W_Level * H_Level; i++) {
            const float LL100 = wav_L0[i] / 327.68f;

            if (LL100 < lp.residshathr) {
                const float kk = aalp * LL100 + alp;
                wav_L0[i] *= (1.f + kk * lp.residsha / 200.f);
            } else if (LL100 < lp.residshathr + tran) {
                wav_L0[i] *= (1.f + (LL100 * ath + bth) / 200.f);
            }

            if (LL100 > lp.residhithr) {
                wav_L0[i] *= (1.f + lp.residhi / 200.f);
            } else if (LL100 > (lp.residhithr - tranh)) {
                wav_L0[i] *= (1.f + (LL100 * athH + bthH) / 200.f);
            }
        }
    }

    if ((lp.residsha > 0.f || lp.residhi > 0.f)) {
        const std::unique_ptr<LabImage> temp(new LabImage(W_Level, H_Level));
#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif

        for (int i = 0; i < H_Level; i++) {
            for (int j = 0; j < W_Level; j++) {
                temp->L[i][j] = wav_L0[i * W_Level + j];
            }
        }

        ImProcFunctions::shadowsHighlights(temp.get(), true, 1, lp.residhi, lp.residsha , 40, sk, lp.residhithr, lp.residshathr);

#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif

        for (int i = 0; i < H_Level; i++) {
            for (int j = 0; j < W_Level; j++) {
                wav_L0[i * W_Level + j] = temp->L[i][j];
            }
        }
    }

    if (contrast != 0.f) {
        double avedbl = 0.0; // use double precision for large summations

#ifdef _OPENMP
        #pragma omp parallel for reduction(+:avedbl) if (multiThread)
#endif
        for (int i = 0; i < W_Level * H_Level; i++) {
            avedbl += static_cast<double>(wav_L0[i]);
        }

        const double avg = LIM01(avedbl / (32768.0 * W_Level * H_Level));
        double contreal = 0.6f * contrast;
        DiagonalCurve resid_contrast({
            DCT_NURBS,
            0, 0,
            avg - avg * (0.6 - contreal / 250.0), avg - avg * (0.6 + contreal / 250.0),
            avg + (1. - avg) * (0.6 - contreal / 250.0), avg + (1. - avg) * (0.6 + contreal / 250.0),
            1, 1
        });
#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif
        for (int i = 0; i < W_Level * H_Level; i++) {
            wav_L0[i] = resid_contrast.getVal(LIM01(wav_L0[i] / 32768.f)) * 32768.0;
        }
    }

    float alow = 1.f;
    float blow = 0.f;

    if (level_hl != level_bl) {
        alow = 1.f / (level_hl - level_bl);
        blow = -alow * level_bl;
    }

    float ahigh = 1.f;
    float bhigh = 0.f;

    if (level_hr != level_br) {
        ahigh = 1.f / (level_hr - level_br);
        bhigh =  -ahigh * level_br;
    }

    if (wavcurvelev  || wavcurvecomp  || wavcurvecompre) {//compress dynamic and blur
        if (wavcurvelev && radlevblur > 0.f && blurena) {
            wavcont(lp, tmp, *wdspot, level_bl, maxlvl, loclevwavCurve, loclevwavutili, loccompwavCurve, loccompwavutili, loccomprewavCurve, loccomprewavutili, radlevblur, 1, 1.f, 0.f, 0.f, 0.f);
        }

        if (wavcurvecomp && comprena) {
            wavcont(lp, tmp, *wdspot, level_bl, maxlvl, loclevwavCurve, loclevwavutili, loccompwavCurve, loccompwavutili, loccomprewavCurve, loccomprewavutili, radlevblur, 2, 1.f, 0.f, sigmadc, deltad);
        }

        if (wavcurvecompre && compreena) {
            wavcont(lp, tmp, *wdspot, level_bl, maxlvl, loclevwavCurve, loclevwavutili, loccompwavCurve, loccompwavutili, loccomprewavCurve, loccomprewavutili, radlevblur, 3, 1.f, thres, 0.f, 0.f);
        }
    }

    if (wavcurvecon && levelena) {//contrast  by levels for luminance
        wavcbd(*wdspot, level_bl, maxlvl, locconwavCurve, locconwavutili, sigm, offs, 1.f, sk);
    }

//edge sharpness begin
    if (lp.edgwena && level_bl == 0 && level_br >= 3 && locedgwavCurve && locedgwavutili && lp.strengthw > 0) { //needs the first levels to work!
        float mean[10];
        float meanN[10];
        float sigma[10];
        float sigmaN[10];
        float MaxP[10];
        float MaxN[10];
        Evaluate2(*wdspot, mean, meanN, sigma, sigmaN, MaxP, MaxN, numThreads);
        float edd = 3.f;
        float eddlow = 15.f;
        float eddlipinfl = 0.005f * lp.edgw + 0.4f;
        float eddlipampl = 1.f + lp.basew / 50.f;

        float *koeLi[12];
        float maxkoeLi[12] = {0.f};

        float *koeLibuffer = new float[12 * H_Level * W_Level]; //12

        for (int i = 0; i < 12; i++) {
            koeLi[i] = &koeLibuffer[i * W_Level * H_Level];
        }

        array2D<float> tmC(W_Level, H_Level);

        float gradw = lp.gradw;
        float tloww = lp.tloww;
        for (int lvl = 0; lvl < 4; lvl++) {
            for (int dir = 1; dir < 4; dir++) {
                const int W_L = wdspot->level_W(lvl);
                const int H_L = wdspot->level_H(lvl);
                float* const* wav_L = wdspot->level_coeffs(lvl);
                calckoe(wav_L[dir], gradw, tloww, koeLi[lvl * 3 + dir - 1], lvl, W_L, H_L, edd, maxkoeLi[lvl * 3 + dir - 1], tmC, true);
                // return convolution KoeLi and maxkoeLi of level 0 1 2 3 and Dir Horiz, Vert, Diag
            }
        }
        
        tmC.free();
        float aamp = 1.f + lp.thigw / 100.f;

        const float alipinfl = (eddlipampl - 1.f) / (1.f - eddlipinfl);
        const float blipinfl = eddlipampl - alipinfl;

        for (int lvl = 0; lvl < 4; lvl++) {
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16)
#endif
            for (int i = 1; i < H_Level - 1; i++) {
                for (int j = 1; j < W_Level - 1; j++) {
                    //treatment of koeLi and maxkoeLi
                    if (lp.lip3) {//Sobel Canny algo improve with parameters
                        // comparison between pixel and neighbors
                        const auto neigh = lp.neiwmet == 1;
                        const auto kneigh = neigh ? 28.f : 38.f;
                        const auto somm = neigh ? 40.f : 50.f;

                        for (int dir = 1; dir < 4; dir++) { //neighbors proxi
                            koeLi[lvl * 3 + dir - 1][i * W_Level + j] = (kneigh * koeLi[lvl * 3 + dir - 1][i * W_Level + j] +
                                                                    2.f * koeLi[lvl * 3 + dir - 1][(i - 1) * W_Level + j] + 2.f * koeLi[lvl * 3 + dir - 1][(i + 1) * W_Level + j] + 2.f * koeLi[lvl * 3 + dir - 1][i * W_Level + j + 1] + 2.f * koeLi[lvl * 3 + dir - 1][i * W_Level + j - 1]
                                                                    + koeLi[lvl * 3 + dir - 1][(i - 1) * W_Level + j - 1] + koeLi[lvl * 3 + dir - 1][(i - 1) * W_Level + j + 1] + koeLi[lvl * 3 + dir - 1][(i + 1) * W_Level + j - 1] + koeLi[lvl * 3 + dir - 1][(i + 1) * W_Level + j + 1]) / somm;
                        }
                    }

                    float interm = 0.f;
                    for (int dir = 1; dir < 4; dir++) {
                        //here I evaluate combination of vert / diag / horiz...we are with multiplicators of the signal
                        interm += SQR(koeLi[lvl * 3 + dir - 1][i * W_Level + j]);
                    }

                    interm = std::sqrt(interm) * 0.57736721f;

                    constexpr float eps = 0.0001f;
                    // I think this double ratio (alph, beta) is better than arctg

                    float alph = koeLi[lvl * 3][i * W_Level + j] / (koeLi[lvl * 3 + 1][i * W_Level + j] + eps); //ratio between horizontal and vertical
                    float beta = koeLi[lvl * 3 + 2][i * W_Level + j] / (koeLi[lvl * 3 + 1][i * W_Level + j] + eps); //ratio between diagonal and horizontal

                    //alph evaluate the direction of the gradient regularity Lipschitz
                    // if = 1 we are on an edge
                    // if 0 we are not
                    // we can change and use log..or Arctg but why ?? we can change if need ...
                    //Liamp=1 for eddlipinfl
                    //liamp > 1 for alp >eddlipinfl and alph < 1
                    //Liamp < 1 for alp < eddlipinfl and alph > 0
                    if (alph > 1.f) {
                        alph = 1.f / alph;
                    }

                    if (beta > 1.f) {
                        beta = 1.f / beta;
                    }

                    //take into account diagonal
                    //if in same value OK
                    //if not no edge or reduction
                    float bet = 1.f;

                    if (alph > eddlipinfl && beta < 0.85f * eddlipinfl) { //0.85 arbitrary value ==> eliminate from edge if H V D too different
                        bet = beta;
                    }

                    float kampli;
                    if (alph > eddlipinfl) {
                        kampli = alipinfl * alph + blipinfl; //If beta low reduce kampli
                        kampli = SQR(bet) * kampli * aamp;
                    } else {
                        kampli = SQR(SQR(alph * bet)) / eddlipinfl; //Strong Reduce if beta low
                        kampli = kampli / aamp;
                    }


                    interm *= kampli;

                    if (interm * eddlow < lp.tloww) {
                        interm = 0.01f;    //eliminate too low values
                    }

                    //we can change this part of algo==> not equal but ponderate
                    koeLi[lvl * 3][i * W_Level + j] = koeLi[lvl * 3 + 1][i * W_Level + j] = koeLi[lvl * 3 + 2][i * W_Level + j] = interm; //new value
                    //here KoeLi contains values where gradient is high and coef high, and eliminate low values...
                }
            }
        }

        constexpr float scales[10] = {1.f, 2.f, 4.f, 8.f, 16.f, 32.f, 64.f, 128.f, 256.f, 512.f};
        float scaleskip[10];

        for (int sc = 0; sc < 10; sc++) {
            scaleskip[sc] = scales[sc] / sk;
        }

        const float rad = lp.radiusw / 60.f; //radius ==> not too high value to avoid artifacts
        float value = lp.strengthw / 8.f; //strength

        if (scaleskip[1] < 1.f) {
            constexpr float atten01234 = 0.80f;
            value *= atten01234 * scaleskip[1];    //for zoom < 100% reduce strength...I choose level 1...but!!
        }

        constexpr float lim0 = 20.f; //arbitrary limit for low radius and level between 2 or 3 to 30 maxi
        float repart = lp.detailw;

        if (lp.edgwmet != 1) {
            float brepart;
            if (lp.edgwmet == 0) {
                brepart = 3.f;
            } else /*if (lp.edgwmet == 2)*/ {
                brepart = 0.5f;    //arbitrary value to increase / decrease repart, between 1 and 0
            }
            if (rad < lim0 / 60.f) {
                const float arepart = - (brepart - 1.f) / (lim0 / 60.f);
                repart *= arepart * rad + brepart;    //linear repartition of repart
            }
        }

        const float bk = 1.f + repart / 50.f;
        constexpr float al10 = 1.0f; //arbitrary value ==> less = take into account high levels
        const float ak = - (bk - al10) / 10.f; //10 = maximum levels

        for (int lvl = 0; lvl < maxlvl; lvl++) {
            if (MaxP[lvl] > 0.f) { //curve
                const int W_L = wdspot->level_W(lvl);
                const int H_L = wdspot->level_H(lvl);
                float* const* wav_L = wdspot->level_coeffs(lvl);
                const float koef = ak * lvl + bk; //modulate for levels : more levels high, more koef low ==> concentrated action on low levels, without or near for high levels
                float expkoef = -pow_F(std::fabs(rad - lvl), koef);   //reduce effect for high levels
                if (lp.edgwmet == 2) {
                    if (rad < lim0 / 60.f && lvl == 0) {
                        expkoef *= abs(repart);    //reduce effect for low values of rad and level=0==> quasi only level 1 is effective
                    }
                } else if (lp.edgwmet == 0) {
                    if (rad < lim0 / 60.f && lvl == 1) {
                        expkoef /= repart;    //increase effect for low values of rad and level=1==> quasi only level 0 is effective
                    }
                }
                //take into account local contrast
                const float refin = value * xexpf(expkoef);
                const float edgePrecalc = 1.f + refin; //estimate edge "pseudo variance"
                constexpr float insigma = 0.666f; //SD
                const float logmax = xlogf(MaxP[lvl]);  //log Max
                const float rapX = (mean[lvl] + sigma[lvl]) / MaxP[lvl]; //rapport between sD / max
                const float inx = xlogf(insigma);
                const float iny = xlogf(rapX);
                const float rap = inx / iny; //koef
                const float asig = 0.166f / sigma[lvl];
                const float bsig = 0.5f - asig * mean[lvl];
                const float amean = 0.5f / mean[lvl];
                constexpr int borderL = 1;
                constexpr float abssd = 4.f; //amplification reference
                constexpr float bbssd = 2.f; //mini ampli
                constexpr float maxamp = 2.5f; //maxi ampli at end
                constexpr float maxampd = 10.f; //maxi ampli at end
                constexpr float a_abssd = (maxamp - abssd) / 0.333f;
                constexpr float b_abssd = maxamp - a_abssd;
                constexpr float da_abssd = (maxampd - abssd) / 0.333f;
                constexpr float db_abssd = maxampd - da_abssd;
                constexpr float am = (abssd - bbssd) / 0.666f;
                const float effect = lp.sigmaed;
                constexpr float offset = 1.f;
                float mea[10];
                calceffect(lvl, mean, sigma, mea, effect, offset);
                float lutFactor;
                const float inVals[] = {0.05f, 0.2f, 0.7f, 1.f, 1.f, 0.8f, 0.5f, 0.3f, 0.2f, 0.1f, 0.05f};
                const auto meaLut = buildMeaLut(inVals, mea, lutFactor);

                for (int dir = 1; dir < 4; dir++) {
#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic, 16) if(multiThread)
#endif
                    for (int i = borderL; i < H_L - borderL; i++) {
                        for (int j = borderL; j < W_L - borderL; j++) {
                            const int k = i * W_L + j;

                            float edge;
                            if (lvl < 4) {
                                edge = 1.f + (edgePrecalc - 1.f) * (koeLi[lvl * 3][k]) / (1.f + 0.9f * maxkoeLi[lvl * 3 + dir - 1]);
                            } else {
                                edge = edgePrecalc;
                            }

                            float absciss = 0.f;
                            if (std::fabs(wav_L[dir][k]) >= mean[lvl] + sigma[lvl]) {  //for max
                                absciss = xexpf((xlogf(std::fabs(wav_L[dir][k])) - logmax) * rap);
                            } else if (std::fabs(wav_L[dir][k]) >= mean[lvl]) {
                                absciss = asig * std::fabs(wav_L[dir][k]) + bsig;
                            } else /*if (std::fabs(wav_L[dir][k]) < mean[lvl])*/ {
                                absciss = amean * std::fabs(wav_L[dir][k]);
                            }

                            // Threshold adjuster settings==> approximative for curve
                            //kmul about average cbrt(3--40 / 10)==>1.5 to 2.5
                            //kmul about SD   10--60  / 35 ==> 2
                            // kmul about low  cbrt((5.f+cp.edg_low)/5.f);==> 1.5
                            // kmul about max ==> 9
                            // we can change these values
                            // result is different not best or bad than threshold slider...but similar
                            float kmul;
                            float kmuld;

                            if (absciss > 0.666f && absciss < 1.f) {
                                kmul = a_abssd * absciss + b_abssd;    //about max  ==> kinterm
                                kmuld = da_abssd * absciss + db_abssd;
                            } else {
                                kmul = kmuld = absciss * am + bbssd;
                            }

                            const float kc = kmul * (locedgwavCurve[absciss * 500.f] - 0.5f);

                            float kinterm;
                            if (kc >= 0.f) {
                                constexpr float reduceeffect = 0.6f;
                                kinterm = 1.f + reduceeffect * kc;    //about 1 to 3 general and big amplification for max (under 0)
                            } else {
                                const float kcd = kmuld * (locedgwavCurve[absciss * 500.f] - 0.5f);
                                kinterm = 1.f - SQR(kcd) / 10.f;
                            }

                            if (kinterm < 0.f) {
                                kinterm = 0.01f;
                            }

                            edge = std::max(edge * kinterm, 1.f);
                            wav_L[dir][k] *= 1.f + (edge - 1.f) * (*meaLut)[std::fabs(wav_L[dir][k]) * lutFactor];
                        }
                    }
                }
            }
        }

        if (koeLibuffer) {
            delete [] koeLibuffer;
        }
    }

//edge sharpness end

    if (locwavCurve && locwavutili && wavcurve) {//simple local contrast in function luminance
        float strengthlc = 1.5f;
        wavlc(*wdspot, level_bl, level_hl, maxlvl, level_hr, level_br, ahigh, bhigh, alow, blow, lp.sigmalc, strengthlc, locwavCurve, numThreads);
    }
    //reconstruct all for L
    wdspot->reconstruct(tmp[0], 1.f);

    bool reconstruct = false;
    if (wavcurvecon && (chromalev != 1.f) && levelena) { // a if need ) {//contrast by levels for chroma a
        // a
        wdspot.reset(new wavelet_decomposition(tmpa[0], bfw, bfh, maxlvl, 1, sk, numThreads, lp.daubLen));
        if (wdspot->memory_allocation_failed()) {
            return;
        }
        wavcbd(*wdspot, level_bl, maxlvl, locconwavCurve, locconwavutili, sigm, offs, chromalev, sk);
        reconstruct = true;
    }
    if (wavcurvelev && radlevblur > 0.f && blurena && chromablu > 0.f && !blurlc) {//chroma blur if need
        // a
        if (!reconstruct) {
            wdspot.reset(new wavelet_decomposition(tmpa[0], bfw, bfh, maxlvl, 1, sk, numThreads, lp.daubLen));
            if (wdspot->memory_allocation_failed()) {
                return;
            }
        }
        wavcont(lp, tmp, *wdspot, level_bl, maxlvl, loclevwavCurve, loclevwavutili, loccompwavCurve, loccompwavutili, loccomprewavCurve, loccomprewavutili, radlevblur, 1, chromablu, 0.f, 0.f, 0.f);
        reconstruct = true;
    }
    if (reconstruct) {
        wdspot->reconstruct(tmpa[0], 1.f);
    }

    reconstruct = false;
    if (wavcurvecon && (chromalev != 1.f) && levelena) { // b if need ) {//contrast by levels for chroma b
        //b
        wdspot.reset(new wavelet_decomposition(tmpb[0], bfw, bfh, maxlvl, 1, sk, numThreads, lp.daubLen));
        if (wdspot->memory_allocation_failed()) {
            return;
        }
        //b
        wavcbd(*wdspot, level_bl, maxlvl, locconwavCurve, locconwavutili, sigm, offs, chromalev, sk);
        reconstruct = true;
    }

    if (wavcurvelev && radlevblur > 0.f && blurena && chromablu > 0.f && !blurlc) {//chroma blur if need
        //b
        if (!reconstruct) {
            wdspot.reset(new wavelet_decomposition(tmpb[0], bfw, bfh, maxlvl, 1, sk, numThreads, lp.daubLen));
            if (wdspot->memory_allocation_failed()) {
                return;
            }
        }
        wavcont(lp, tmp, *wdspot, level_bl, maxlvl, loclevwavCurve, loclevwavutili, loccompwavCurve, loccompwavutili, loccomprewavCurve, loccomprewavutili, radlevblur, 1, chromablu, 0.f, 0.f, 0.f);
        reconstruct = true;
    }
    if (reconstruct) {
        wdspot->reconstruct(tmpb[0], 1.f);
    }
    
    
    //gamma and slope residual image - be careful memory
    bool tonecur = false;
    const Glib::ustring profile = params->icm.workingProfile;
    bool isworking = (profile == "sRGB" || profile == "Adobe RGB" || profile == "ProPhoto" || profile == "WideGamut" || profile == "BruceRGB" || profile == "Beta RGB" || profile == "BestRGB" || profile == "Rec2020" || profile == "ACESp0" || profile == "ACESp1");

    if (isworking && (lp.residgam != 2.4f || lp.residslop != 12.92f)) {
        tonecur = true;
    }
    
    if(tonecur) {
        std::unique_ptr<wavelet_decomposition> wdspotL(new wavelet_decomposition(tmp[0], bfw, bfh, maxlvl, 1, sk, numThreads, lp.daubLen));
        if (wdspotL->memory_allocation_failed()) {
            return;
        }
        std::unique_ptr<wavelet_decomposition> wdspota(new wavelet_decomposition(tmpa[0], bfw, bfh, maxlvl, 1, sk, numThreads, lp.daubLen));
        if (wdspota->memory_allocation_failed()) {
            return;
        }
        std::unique_ptr<wavelet_decomposition> wdspotb(new wavelet_decomposition(tmpb[0], bfw, bfh, maxlvl, 1, sk, numThreads, lp.daubLen));
        if (wdspotb->memory_allocation_failed()) {
            return;
        }
        int W_Level = wdspotL->level_W(0);
        int H_Level = wdspotL->level_H(0);
        float *wav_L0 = wdspotL->get_coeff0();
        float *wav_a0 = wdspota->get_coeff0();
        float *wav_b0 = wdspotb->get_coeff0();

        const std::unique_ptr<LabImage> labresid(new LabImage(W_Level, H_Level));
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
        for (int y = 0; y < H_Level; y++) {
            for (int x = 0; x < W_Level; x++) {
                labresid->L[y][x] = wav_L0[y * W_Level + x];
                labresid->a[y][x] = wav_a0[y * W_Level + x];
                labresid->b[y][x] = wav_b0[y * W_Level + x];
            }
        }

        Imagefloat *tmpImage = nullptr;
        tmpImage = new Imagefloat(W_Level, H_Level);
        lab2rgb(*labresid, *tmpImage, params->icm.workingProfile);
        Glib::ustring prof = params->icm.workingProfile;
        cmsHTRANSFORM dummy = nullptr;
        int ill =0;
        workingtrc(tmpImage, tmpImage, W_Level, H_Level, -5, prof, 2.4, 12.92310, ill, 0, dummy, true, false, false);
        workingtrc(tmpImage, tmpImage, W_Level, H_Level, 1, prof, lp.residgam, lp.residslop, ill, 0, dummy, false, true, true);//be careful no gamut control
        rgb2lab(*tmpImage, *labresid, params->icm.workingProfile);
        delete tmpImage;

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
        for (int y = 0; y < H_Level; y++) {
            for (int x = 0; x < W_Level; x++) {
                wav_L0[y * W_Level + x] = labresid->L[y][x];
                wav_a0[y * W_Level + x] = labresid->a[y][x];
                wav_b0[y * W_Level + x] = labresid->b[y][x];
            }
        }

        wdspotL->reconstruct(tmp[0], 1.f);
        wdspota->reconstruct(tmpa[0], 1.f);
        wdspotb->reconstruct(tmpb[0], 1.f);
    }
}


void ImProcFunctions::fftw_denoise(int sk, int GW, int GH, int max_numblox_W, int min_numblox_W, float **tmp1, array2D<float> *Lin, int numThreads, const struct local_params & lp, int chrom)
{
   // BENCHFUN

    fftwf_plan plan_forward_blox[2];
    fftwf_plan plan_backward_blox[2];

    array2D<float> tilemask_in(TS, TS);
    array2D<float> tilemask_out(TS, TS);

    float *Lbloxtmp  = reinterpret_cast<float*>(fftwf_malloc(max_numblox_W * TS * TS * sizeof(float)));
    float *fLbloxtmp = reinterpret_cast<float*>(fftwf_malloc(max_numblox_W * TS * TS * sizeof(float)));
    float params_Ldetail = 0.f;

    int nfwd[2] = {TS, TS};

    //for DCT:
    fftw_r2r_kind fwdkind[2] = {FFTW_REDFT10, FFTW_REDFT10};
    fftw_r2r_kind bwdkind[2] = {FFTW_REDFT01, FFTW_REDFT01};

    // Creating the plans with FFTW_MEASURE instead of FFTW_ESTIMATE speeds up the execute a bit
    plan_forward_blox[0]  = fftwf_plan_many_r2r(2, nfwd, max_numblox_W, Lbloxtmp, nullptr, 1, TS * TS, fLbloxtmp, nullptr, 1, TS * TS, fwdkind, FFTW_MEASURE | FFTW_DESTROY_INPUT);
    plan_backward_blox[0] = fftwf_plan_many_r2r(2, nfwd, max_numblox_W, fLbloxtmp, nullptr, 1, TS * TS, Lbloxtmp, nullptr, 1, TS * TS, bwdkind, FFTW_MEASURE | FFTW_DESTROY_INPUT);
    plan_forward_blox[1]  = fftwf_plan_many_r2r(2, nfwd, min_numblox_W, Lbloxtmp, nullptr, 1, TS * TS, fLbloxtmp, nullptr, 1, TS * TS, fwdkind, FFTW_MEASURE | FFTW_DESTROY_INPUT);
    plan_backward_blox[1] = fftwf_plan_many_r2r(2, nfwd, min_numblox_W, fLbloxtmp, nullptr, 1, TS * TS, Lbloxtmp, nullptr, 1, TS * TS, bwdkind, FFTW_MEASURE | FFTW_DESTROY_INPUT);
    fftwf_free(Lbloxtmp);
    fftwf_free(fLbloxtmp);
    const int border = rtengine::max(2, TS / 16);

    for (int i = 0; i < TS; ++i) {
        float i1 = abs((i > TS / 2 ? i - TS + 1 : i));
        float vmask = (i1 < border ? SQR(sin((rtengine::RT_PI_F * i1) / (2 * border))) : 1.0f);
        float vmask2 = (i1 < 2 * border ? SQR(sin((rtengine::RT_PI_F * i1) / (2 * border))) : 1.0f);

        for (int j = 0; j < TS; ++j) {
            float j1 = abs((j > TS / 2 ? j - TS + 1 : j));
            tilemask_in[i][j] = (vmask * (j1 < border ? SQR(sin((rtengine::RT_PI_F * j1) / (2 * border))) : 1.0f)) + epsilonw;
            tilemask_out[i][j] = (vmask2 * (j1 < 2 * border ? SQR(sin((rtengine::RT_PI_F * j1) / (2 * border))) : 1.0f)) + epsilonw;

        }
    }


    float *LbloxArray[numThreads];
    float *fLbloxArray[numThreads];



    const int numblox_W = ceil((static_cast<float>(GW)) / offset) + 2;
    const int numblox_H = ceil((static_cast<float>(GH)) / offset) + 2;


    //residual between input and denoised L channel
    array2D<float> Ldetail(GW, GH, ARRAY2D_CLEAR_DATA);
    array2D<float> totwt(GW, GH, ARRAY2D_CLEAR_DATA); //weight for combining DCT blocks
    array2D<float> prov(GW, GH, ARRAY2D_CLEAR_DATA);

    for (int i = 0; i < numThreads; ++i) {
        LbloxArray[i]  = reinterpret_cast<float*>(fftwf_malloc(max_numblox_W * TS * TS * sizeof(float)));
        fLbloxArray[i] = reinterpret_cast<float*>(fftwf_malloc(max_numblox_W * TS * TS * sizeof(float)));
    }

#ifdef _OPENMP
    int masterThread = omp_get_thread_num();
#endif
#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef _OPENMP
        int subThread = masterThread * 1 + omp_get_thread_num();
#else
        int subThread = 0;
#endif
        float *Lblox = LbloxArray[subThread];
        float *fLblox = fLbloxArray[subThread];
        float pBuf[GW + TS + 2 * offset] ALIGNED16;
#ifdef _OPENMP
        #pragma omp for
#endif
        for (int vblk = 0; vblk < numblox_H; ++vblk) {

            int top = (vblk - 1) * offset;
            float * datarow = pBuf + offset;

            for (int i = 0; i < TS; ++i) {
                int row = top + i;
                int rr = row;

                if (row < 0) {
                    rr = rtengine::min(-row, GH - 1);
                } else if (row >= GH) {
                    rr = rtengine::max(0, 2 * GH - 2 - row);
                }

                for (int j = 0; j < GW; ++j) {
                    datarow[j] = ((*Lin)[rr][j] - tmp1[rr][j]);
                    prov[rr][j] = std::fabs(tmp1[rr][j]);

                }

                for (int j = -1 * offset; j < 0; ++j) {
                    datarow[j] = datarow[rtengine::min(-j, GW - 1)];
                }

                for (int j = GW; j < GW + TS + offset; ++j) {
                    datarow[j] = datarow[rtengine::max(0, 2 * GW - 2 - j)];
                }//now we have a padded data row

                //now fill this row of the blocks with Lab high pass data
                for (int hblk = 0; hblk < numblox_W; ++hblk) {
                    int left = (hblk - 1) * offset;
                    int indx = (hblk) * TS; //index of block in malloc

                    if (top + i >= 0 && top + i < GH) {
                        int j;

                        for (j = 0; j < rtengine::min((-left), TS); ++j) {
                            Lblox[(indx + i)*TS + j] = tilemask_in[i][j] * datarow[left + j]; // luma data
                        }

                        for (; j < rtengine::min(TS, GW - left); ++j) {
                            Lblox[(indx + i)*TS + j] = tilemask_in[i][j] * datarow[left + j]; // luma data
                            totwt[top + i][left + j] += tilemask_in[i][j] * tilemask_out[i][j];
                        }

                        for (; j < TS; ++j) {
                            Lblox[(indx + i)*TS + j] = tilemask_in[i][j] * datarow[left + j]; // luma data
                        }
                    } else {
                        for (int j = 0; j < TS; ++j) {
                            Lblox[(indx + i)*TS + j] = tilemask_in[i][j] * datarow[left + j]; // luma data
                        }
                    }

                }

            }//end of filling block row

            //fftwf_print_plan (plan_forward_blox);
            if (numblox_W == max_numblox_W) {
                fftwf_execute_r2r(plan_forward_blox[0], Lblox, fLblox);    // DCT an entire row of tiles
            } else {
                fftwf_execute_r2r(plan_forward_blox[1], Lblox, fLblox);    // DCT an entire row of tiles
            }

            // now process the vblk row of blocks for noise reduction

            float noisevar_Ldetail = 1.f;

            if (chrom == 0) {
                params_Ldetail = rtengine::min(float(lp.noiseldetail), 99.9f);    // max out to avoid div by zero when using noisevar_Ldetail as divisor
                noisevar_Ldetail = SQR(static_cast<float>(SQR(100.f - params_Ldetail) + 50.f * (100.f - params_Ldetail)) * TS * 0.5f);
            } else if (chrom == 1) {
                params_Ldetail = rtengine::min(float(lp.noisechrodetail), 99.9f);
                noisevar_Ldetail = 100.f * rtengine::SQR((static_cast<float>(SQR(100.f - params_Ldetail)) * TS * 0.5f));//to test ???
            }

            for (int hblk = 0; hblk < numblox_W; ++hblk) {
                ImProcFunctions::RGBtile_denoise(fLblox, hblk, noisevar_Ldetail);

            }//end of horizontal block loop

            //now perform inverse FT of an entire row of blocks
            if (numblox_W == max_numblox_W) {
                fftwf_execute_r2r(plan_backward_blox[0], fLblox, Lblox);    //for DCT
            } else {
                fftwf_execute_r2r(plan_backward_blox[1], fLblox, Lblox);    //for DCT
            }

            int topproc = (vblk - 1) * offset;

            //add row of blocks to output image tile
            ImProcFunctions::RGBoutput_tile_row(Lblox, Ldetail, tilemask_out, GH, GW, topproc);

        }//end of vertical block loop
    }

    //Threshold DCT from Alberto Grigio, adapted to Rawtherapee 
    const int detail_thresh = lp.detailthr;
    array2D<float> mask;

    if (detail_thresh > 0) {
        mask(GW, GH);
        if (lp.usemask) {//with Laplacian
            float amount = LIM01(float(detail_thresh)/100.f);
            float thr = (1.f - amount);
            float alph = params_Ldetail / 100.f;
            array2D<float> LL(GW, GH, prov, ARRAY2D_BYREFERENCE);
            laplacian(LL, mask, GW, GH, 25.f, 20000.f, amount, false);
            for (int i = 0; i < GH; ++i) {
                for (int j = 0; j < GW; ++j) {
                    mask[i][j] = LIM01(mask[i][j]+ thr);
                }
            }
            for (int i = 0; i < 3; ++i) {
                boxblur(static_cast<float**>(mask), static_cast<float**>(mask), 10 / sk, GW, GH, false);
                
            }
            for (int i = 0; i < GH; ++i) {
                for (int j = 0; j < GW; ++j) {
                    float k = 1.f - mask[i][j] * alph;
                    mask[i][j] = 1.f - (k * k);
                }
            }
        } else {//with blend mask
            float thr = log2lin(float(detail_thresh) / 200.f, 100.f);
            buildBlendMask(prov, mask, GW, GH, thr);
#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
            {
                gaussianBlur(mask, mask, GW, GH, 20.0 / sk);
            }
            array2D<float> m2(GW, GH);
            constexpr float alfa = 0.856f;
            const float beta = 1.f + std::sqrt(log2lin(thr, 100.f));
            buildGradientsMask(GW, GH, prov, m2, params_Ldetail / 100.f, 7, 3, alfa, beta, multiThread);
            for (int i = 0; i < GH; ++i) {
                for (int j = 0; j < GW; ++j) {
                    mask[i][j] *= m2[i][j];
                }
            }
        } 
    }


#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif
    for (int i = 0; i < GH; ++i) {
        for (int j = 0; j < GW; ++j) {
            float d = Ldetail[i][j] / totwt[i][j];

            if (detail_thresh > 0) {
                d *= mask[i][j];
            }

            //may want to include masking threshold for large hipass data to preserve edges/detail
            tmp1[i][j] += d;
        }
    }

    mask.free();
//end Threshold DCT


    delete Lin;


    for (int i = 0; i < numThreads; ++i) {
        fftwf_free(LbloxArray[i]);
        fftwf_free(fLbloxArray[i]);
    }

    fftwf_destroy_plan(plan_forward_blox[0]);
    fftwf_destroy_plan(plan_backward_blox[0]);
    fftwf_destroy_plan(plan_forward_blox[1]);
    fftwf_destroy_plan(plan_backward_blox[1]);
    fftwf_cleanup();


}

void ImProcFunctions::DeNoise(int call, int aut,  bool noiscfactiv, const struct local_params & lp, LabImage * originalmaskbl, LabImage *  bufmaskblurbl, int levred, float huerefblur, float lumarefblur, float chromarefblur, LabImage * original, LabImage * transformed,
    int cx, int cy, int sk, const LocwavCurve& locwavCurvehue, bool locwavhueutili, float& highresi, float& nresi, float& highresi46, float& nresi46, float& Lhighresi, float& Lnresi, float& Lhighresi46, float& Lnresi46)
{
    BENCHFUN
//local denoise
    //all these variables are to prevent use of denoise when non necessary
    // but with qualmet = 2 (default for best quality) we must denoise chroma with little values to prevent artifacts due to variations of Hue
    // but if user select voluntary denoise, it is that choice the good (priority)
    bool execcolor = (lp.chro != 0.f || lp.ligh != 0.f || lp.cont != 0); // only if one slider or more is engaged
    bool execbdl = (lp.mulloc[0] != 1.f || lp.mulloc[1] != 1.f || lp.mulloc[2] != 1.f || lp.mulloc[3] != 1.f || lp.mulloc[4] != 1.f || lp.mulloc[5] != 1.f) ;//only if user want cbdl
    bool execdenoi = noiscfactiv && ((lp.colorena && execcolor) || (lp.tonemapena && lp.strengt != 0.f) || (lp.cbdlena && execbdl) || (lp.sfena && lp.strng > 0.f) || (lp.lcena && lp.lcamount > 0.f) || (lp.sharpena && lp.shrad > 0.42) || (lp.retiena && lp.str > 0.f)  || (lp.exposena && lp.expcomp != 0.f)  || (lp.expvib && lp.past != 0.f));
    bool execmaskden = (lp.showmaskblmet == 2 || lp.enablMask || lp.showmaskblmet == 3 || lp.showmaskblmet == 4) && lp.smasktyp != 0;

//    const int ys = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
//    const int ye = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
//    const int xs = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
//    const int xe = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
//    const int hspot = ye - ys;
//    const int wspot = xe - xs;

   if (((lp.noiself > 0.f || lp.noiself0 > 0.f || lp.noiself2 > 0.f || lp.nlstr > 0 || lp.wavcurvedenoi || lp.noiselc > 0.f || lp.noisecf > 0.f || lp.noisecc > 0.f
            || execmaskden || aut == 1 || aut == 2) && lp.denoiena && lp.quamet != 3) || execdenoi) {  // sk == 1 ??

        StopWatch Stop1("locallab Denoise called");

        if (aut == 0) {
            MyMutex::MyLock lock(*fftwMutex);
        }

        if (lp.noisecf >= 0.01f || lp.noisecc >= 0.01f || aut == 1 || aut == 2) {
            noiscfactiv = false;
            levred = 7;
        }

        int GW = transformed->W;
        int GH = transformed->H;

        bool HHhuecurve = false;

        if (locwavCurvehue && locwavhueutili) {
            for (int i = 0; i < 500; i++) {
                if (locwavCurvehue[i] != 0.5f) {
                    HHhuecurve = true;
                    break;
                }
            }
        }



#ifdef _OPENMP
        const int numThreads = omp_get_max_threads();
#else
        const int numThreads = 1;

#endif
            int minwin = rtengine::min(GW, GH);
            int maxlevelspot = 10;//maximum possible
            bool isnois = true;
            // adap maximum level wavelet to size of crop
            while ((1 << maxlevelspot) >= (minwin * sk) && maxlevelspot  > 1) {
                --maxlevelspot ;
            }

            levred = rtengine::min(levred, maxlevelspot);
            if(levred < 7) {//If windows preview or detail window too small exit to avoid artifacts
                isnois = false;
                if(lp.quamet == 2) {
                    isnois = true;
                }
            }

        if (call == 1 && ((GW >= mDEN && GH >= mDEN  && isnois) || lp.quamet == 2)) {


            LabImage tmp1(transformed->W, transformed->H);
            LabImage tmp2(transformed->W, transformed->H);
            tmp2.clear();

            array2D<float> *Lin = nullptr;
            array2D<float> *Ain = nullptr;
            array2D<float> *Bin = nullptr;

            int max_numblox_W = ceil((static_cast<float>(GW)) / offset) + 2;
            // calculate min size of numblox_W.
            int min_numblox_W = ceil((static_cast<float>(GW)) / offset) + 2;

            for (int ir = 0; ir < GH; ir++)
                for (int jr = 0; jr < GW; jr++) {
                    tmp1.L[ir][jr] = original->L[ir][jr];
                    tmp1.a[ir][jr] = original->a[ir][jr];
                    tmp1.b[ir][jr] = original->b[ir][jr];
                }
            if(lp.nlstr > 0) {
                NLMeans(tmp1.L, lp.nlstr, lp.nldet, lp.nlpat, lp.nlrad, lp.nlgam, GW, GH, float (sk), multiThread);
            }

            float gamma = lp.noisegam;
            rtengine::GammaValues g_a; //gamma parameters
            double pwr = 1.0 / (double) lp.noisegam;//default 3.0 - gamma Lab
            double ts = 9.03296;//always the same 'slope' in the extreme shadows - slope Lab
            rtengine::Color::calcGamma(pwr, ts, g_a); // call to calcGamma with selected gamma and slope

            if(gamma > 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                for (int y = 0; y < GH; ++y) {
                    int x = 0;
#ifdef __SSE2__
                    for (; x < GW - 3; x += 4) {
                        STVFU(tmp1.L[y][x], F2V(32768.f) * igammalog(LVFU(tmp1.L[y][x]) / F2V(32768.f), F2V(gamma), F2V(ts), F2V(g_a[2]), F2V(g_a[4])));
                    }
#endif
                    for (;x < GW; ++x) {
                        tmp1.L[y][x] = 32768.f * igammalog(tmp1.L[y][x] / 32768.f, gamma, ts, g_a[2], g_a[4]);
                    }
                }
            }

            //  int DaubLen = 6;

            int levwavL = levred;
            int skip = 1;

            wavelet_decomposition Ldecomp(tmp1.L[0], tmp1.W, tmp1.H, levwavL, 1, skip, numThreads, lp.daubLen);
            wavelet_decomposition adecomp(tmp1.a[0], tmp1.W, tmp1.H, levwavL, 1, skip, numThreads, lp.daubLen);
            wavelet_decomposition bdecomp(tmp1.b[0], tmp1.W, tmp1.H, levwavL, 1, skip, numThreads, lp.daubLen);

            float madL[10][3];
            int edge = 2;

            if (!Ldecomp.memory_allocation_failed()) {
#ifdef _OPENMP
            //    #pragma omp parallel for schedule(dynamic) collapse(2) if (multiThread)
#endif
                for (int lvl = 0; lvl < levred; lvl++) {
                    for (int dir = 1; dir < 4; dir++) {
                        int Wlvl_L = Ldecomp.level_W(lvl);
                        int Hlvl_L = Ldecomp.level_H(lvl);
                        const float* const* WavCoeffs_L = Ldecomp.level_coeffs(lvl);
                        madL[lvl][dir - 1] = SQR(Mad(WavCoeffs_L[dir], Wlvl_L * Hlvl_L));
                    }
                }

                float vari[levred];
                float mxsl = 0.f;
                //      float mxsfl = 0.f;

                edge = 2;
                vari[0] = 0.8f * SQR((lp.noiself0 / 125.f) * (1.f + lp.noiself0 / 25.f));
                vari[1] = 0.8f * SQR((lp.noiself / 125.f) * (1.f + lp.noiself / 25.f));
                vari[2] = 0.8f * SQR((lp.noiself2 / 125.f) * (1.f + lp.noiself2 / 25.f));

                vari[3] = 0.8f * SQR((lp.noiselc / 125.f) * (1.f + lp.noiselc / 25.f));
                vari[4] = 1.f * SQR((lp.noiselc4 / 125.f) * (1.f + lp.noiselc4 / 25.f));
                vari[5] = 1.5f * SQR((lp.noiselc5 / 125.f) * (1.f + lp.noiselc5 / 25.f));
                vari[6] = 2.5f * SQR((lp.noiselc6 / 125.f) * (1.f + lp.noiselc6 / 25.f));

                {
                    float kr3 = 0.f;

                        if (lp.noiselc < 30.f) {
                            kr3 = 0.f;
                        } else if (lp.noiselc < 50.f) {
                            kr3 = 0.5f;
                        } else if (lp.noiselc < 70.f) {
                            kr3 = 0.7f;
                        } else {
                            kr3 = 1.f;
                        }

                    vari[0] = rtengine::max(0.000001f, vari[0]);
                    vari[1] = rtengine::max(0.000001f, vari[1]);
                    vari[2] = rtengine::max(0.000001f, vari[2]);
                    vari[3] = rtengine::max(0.000001f, kr3 * vari[3]);
                    vari[4] = rtengine::max(0.000001f, vari[4]);
                    vari[5] = rtengine::max(0.000001f, vari[5]);
                    vari[6] = rtengine::max(0.000001f, vari[6]);

                    float* noisevarlum = new float[GH * GW];
                    float* noisevarhue = new float[GH * GW];
                    int GW2 = (GW + 1) / 2;

                    float nvlh[13] = {1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 0.7f, 0.5f}; //high value
                    float nvll[13] = {0.1f, 0.15f, 0.2f, 0.25f, 0.3f, 0.35f, 0.4f, 0.45f, 0.7f, 0.8f, 1.f, 1.f, 1.f}; //low value

                    float seuillow = 3000.f;//low
                    float seuilhigh = 18000.f;//high
                    int i = 10 - lp.noiselequal;
                    float ac = (nvlh[i] - nvll[i]) / (seuillow - seuilhigh);
                    float bc = nvlh[i] - seuillow * ac;
                    //ac and bc for transition
#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                    for (int ir = 0; ir < GH; ir++)
                        for (int jr = 0; jr < GW; jr++) {
                            float lN = tmp1.L[ir][jr];

                            if (lN < seuillow) {
                                noisevarlum[(ir >> 1)*GW2 + (jr >> 1)] =  nvlh[i];
                            } else if (lN < seuilhigh) {
                                noisevarlum[(ir >> 1)*GW2 + (jr >> 1)] = ac * lN + bc;
                            } else {
                                noisevarlum[(ir >> 1)*GW2 + (jr >> 1)] =  nvll[i];
                            }
                        }

                    if(lp.enablMask && lp.lnoiselow !=1.f && lp.smasktyp != 0) {
                        //this code has been reviewed by Ingo in september 2020 PR5903
                        float higc;
                        float hig = lp.thrhigh;
                        calcdif(hig, higc);
                        float low = lp.thrlow;
                        float lowc;
                        calcdif(low, lowc);
                         
                        if(higc < lowc) {
                            higc = lowc + 0.01f;
                        }

                        float alow = -(lp.lnoiselow - 1.f) / lowc;
                        float blow = lp.lnoiselow;
                        float ahigh = 0.9999f / (higc - 100.f);
                        float bhigh = 1.f - higc * ahigh;
                        
#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                        for (int ir = 0; ir < GH; ir++)
                            for (int jr = 0; jr < GW; jr++) {
                                const float lM = bufmaskblurbl->L[ir][jr];
                                const float lmr = lM / 327.68f;

                                if (lM < 327.68f * lowc) {
                                    noisevarlum[(ir >> 1) * GW2 + (jr >> 1)] *= alow * lmr + blow;
                                } else if (lM < 327.68f * higc) {
                                    // do nothing - denoise not change
                                } else {
                                    noisevarlum[(ir >> 1) * GW2 + (jr >> 1)] *= ahigh * lmr + bhigh;
                                }
                        }
                    }


                    if(HHhuecurve) {
                    //same code as in wavelet levels
                        
#ifdef _OPENMP
        #pragma omp parallel for
#endif
                        for (int ir = 0; ir < GH; ir++)
                            for (int jr = 0; jr < GW; jr++) {
                                float hueG = xatan2f(tmp1.b[ir][jr], tmp1.a[ir][jr]);
                                float valparam = 2.f * (locwavCurvehue[500.f * static_cast<float>(Color::huelab_to_huehsv2(hueG))] - 0.5f);  //get H=f(H)
                                noisevarhue[(ir >> 1)*GW2 + (jr >> 1)] = 1.f +  valparam;
                                noisevarlum[(ir >> 1)*GW2 + (jr >> 1)] *= noisevarhue[(ir >> 1)*GW2 + (jr >> 1)];
                            }
                    }


                    if ((lp.quamet == 0  && aut == 0) || (mxsl < 1.f && (aut == 1 || aut == 2))) {
                        WaveletDenoiseAllL(Ldecomp, noisevarlum, madL, vari, edge, numThreads);

                    } else if (lp.quamet == 1){

                        WaveletDenoiseAll_BiShrinkL(Ldecomp, noisevarlum, madL, vari, edge, numThreads);

                        WaveletDenoiseAllL(Ldecomp, noisevarlum, madL, vari, edge, numThreads);

                    }

                    delete[] noisevarlum;
                    delete[] noisevarhue;

                }
            }


            float variC[levred];
            float variCb[levred];

            float noisecfr = lp.noisecf;
            float noiseccr = lp.noisecc;

            if (lp.adjch > 0.f) {
                noisecfr = lp.noisecf + 0.1f * lp.adjch;
                noiseccr = lp.noisecc + 0.1f * lp.adjch;
            }

            float noisecfb = lp.noisecf;
            float noiseccb = lp.noisecc;

            if (lp.adjch < 0.f) {
                noisecfb = lp.noisecf - 0.1f * lp.adjch;
                noiseccb = lp.noisecc - 0.1f * lp.adjch;
            }


            if (noisecfr < 0.f) {
                noisecfr = 0.00001f;
            }

            if (noiseccr < 0.f) {
                noiseccr = 0.00001f;
            }

            if (noisecfb < 0.f) {
                noisecfb = 0.00001f;
            }

            if (noiseccb < 0.f) {
                noiseccb = 0.00001f;
            }

            if (!adecomp.memory_allocation_failed() && !bdecomp.memory_allocation_failed()) {
                float maxccoarse = 0.f;

            edge = 2;
            variC[0] = SQR(noisecfr);
            variC[1] = SQR(noisecfr);
            variC[2] = SQR(noisecfr);
            variC[3] = SQR(1.2f * noisecfr);
            variC[4] = SQR(noisecfr);
            variC[5] = SQR(1.2f * noiseccr);
            variC[6] = SQR(1.5f * noiseccr);

            variCb[0] = SQR(noisecfb);
            variCb[1] = SQR(noisecfb);
            variCb[2] = SQR(noisecfb);
            variCb[3] = SQR(noisecfb);
            variCb[4] = SQR(noisecfb);
            variCb[5] = SQR(1.2f * noiseccb);
            variCb[6] = SQR(1.5f * noiseccb);


                {
                    float minic = 0.000001f;

                    if (noiscfactiv) {
                        minic = 0.1f;//only for artifact shape detection
                    }

                    float k1 = 0.f;
                    float k2 = 0.f;
                    float k3 = 0.f;

                        if (lp.noisecf < 0.2f) {
                            k1 = 0.05f;
                            k2 = 0.f;
                            k3 = 0.f;
                        } else if (lp.noisecf < 0.3f) {
                            k1 = 0.1f;
                            k2 = 0.0f;
                            k3 = 0.f;
                        } else if (lp.noisecf < 0.5f) {
                            k1 = 0.2f;
                            k2 = 0.1f;
                            k3 = 0.f;
                        } else if (lp.noisecf < 0.8f) {
                            k1 = 0.3f;
                            k2 = 0.25f;
                            k3 = 0.f;
                        } else if (lp.noisecf < 1.f) {
                            k1 = 0.4f;
                            k2 = 0.25f;
                            k3 = 0.1f;
                        } else if (lp.noisecf < 2.f) {
                            k1 = 0.5f;
                            k2 = 0.3f;
                            k3 = 0.15f;
                        } else if (lp.noisecf < 3.f) {
                            k1 = 0.6f;
                            k2 = 0.45f;
                            k3 = 0.3f;
                        } else if (lp.noisecf < 4.f) {
                            k1 = 0.7f;
                            k2 = 0.5f;
                            k3 = 0.4f;
                        } else if (lp.noisecf < 5.f) {
                            k1 = 0.8f;
                            k2 = 0.6f;
                            k3 = 0.5f;
                        } else if (lp.noisecf < 6.f) {
                            k1 = 0.85f;
                            k2 = 0.7f;
                            k3 = 0.6f;
                        } else if (lp.noisecf < 8.f) {
                            k1 = 0.9f;
                            k2 = 0.8f;
                            k3 = 0.7f;
                        } else if (lp.noisecf < 10.f) {
                            k1 = 1.f;
                            k2 = 1.f;
                            k3 = 0.9f;

                        } else {
                            k1 = 1.f;
                            k2 = 1.f;
                            k3 = 1.5f;
                        }

                    variC[0] = rtengine::max(minic, variC[0]);
                    variC[1] = rtengine::max(minic, k1 * variC[1]);
                    variC[2] = rtengine::max(minic, k2 * variC[2]);
                    variC[3] = rtengine::max(minic, k3 * variC[3]);

                    variCb[0] = rtengine::max(minic, variCb[0]);
                    variCb[1] = rtengine::max(minic, k1 * variCb[1]);
                    variCb[2] = rtengine::max(minic, k2 * variCb[2]);
                    variCb[3] = rtengine::max(minic, k3 * variCb[3]);

                        float k4 = 0.f;
                        float k5 = 0.f;
                        float k6 = 0.f;

                        if (lp.noisecc < 0.2f) {
                            k4 = 0.1f;
                            k5 = 0.02f;
                        } else if (lp.noisecc < 0.5f) {
                            k4 = 0.15f;
                            k5 = 0.05f;
                        } else if (lp.noisecc < 1.f) {
                            k4 = 0.15f;
                            k5 = 0.1f;
                        } else if (lp.noisecc < 3.f) {
                            k4 = 0.3f;
                            k5 = 0.15f;
                        } else if (lp.noisecc < 4.f) {
                            k4 = 0.6f;
                            k5 = 0.4f;
                        } else if (lp.noisecc < 6.f) {
                            k4 = 0.8f;
                            k5 = 0.6f;
                        } else {
                            k4 = 1.f;
                            k5 = 1.f;
                        }

                        variC[4] = rtengine::max(0.000001f, k4 * variC[4]);
                        variC[5] = rtengine::max(0.000001f, k5 * variC[5]);
                        variCb[4] = rtengine::max(0.000001f, k4 * variCb[4]);
                        variCb[5] = rtengine::max(0.000001f, k5 * variCb[5]);

                        if (lp.noisecc < 4.f) {
                            k6 = 0.f;
                        } else if (lp.noisecc < 5.f) {
                            k6 = 0.4f;
                        } else if (lp.noisecc < 6.f) {
                            k6 = 0.7f;
                        } else {
                            k6 = 1.f;
                        }

                        variC[6] = rtengine::max(0.00001f, k6 * variC[6]);
                        variCb[6] = rtengine::max(0.00001f, k6 * variCb[6]);


                    float* noisevarchrom = new float[GH * GW];
                    //noisevarchrom in function chroma
                    int GW2 = (GW + 1) / 2;
                    float nvch = 0.6f;//high value
                    float nvcl = 0.1f;//low value

                    if (lp.noisecf > 100.f) {
                        nvch = 0.8f;
                        nvcl = 0.4f;
                    }

                    float seuil = 4000.f;//low
                    float seuil2 = 15000.f;//high
                    //ac and bc for transition
                    float ac = (nvch - nvcl) / (seuil - seuil2);
                    float bc = nvch - seuil * ac;
#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                    for (int ir = 0; ir < GH; ir++)
                        for (int jr = 0; jr < GW; jr++) {
                            float cN = std::sqrt(SQR(tmp1.a[ir][jr]) + SQR(tmp1.b[ir][jr]));

                            if (cN < seuil) {
                                noisevarchrom[(ir >> 1)*GW2 + (jr >> 1)] =  nvch;
                            } else if (cN < seuil2) {
                                noisevarchrom[(ir >> 1)*GW2 + (jr >> 1)] = ac * cN + bc;
                            } else {
                                noisevarchrom[(ir >> 1)*GW2 + (jr >> 1)] =  nvcl;
                            }
                        }


                    float noisevarab_r = 100.f; //SQR(lp.noisecc / 10.0);

                    if ((lp.quamet == 0  && aut == 0) || (maxccoarse < 0.1f && (aut == 1 || aut == 2)))  {
                        WaveletDenoiseAllAB(Ldecomp, adecomp, noisevarchrom, madL, variC, edge, noisevarab_r, true, false, false, numThreads);
                        WaveletDenoiseAllAB(Ldecomp, bdecomp, noisevarchrom, madL, variCb, edge, noisevarab_r, true, false, false, numThreads);
                    } else if (lp.quamet == 1){
                        WaveletDenoiseAll_BiShrinkAB(Ldecomp, adecomp, noisevarchrom, madL, variC, edge, noisevarab_r, true, false, false, numThreads);
                        WaveletDenoiseAllAB(Ldecomp, adecomp, noisevarchrom, madL, variC, edge, noisevarab_r, true, false, false, numThreads);

                        WaveletDenoiseAll_BiShrinkAB(Ldecomp, bdecomp, noisevarchrom, madL, variCb, edge, noisevarab_r, true, false, false, numThreads);
                        WaveletDenoiseAllAB(Ldecomp, bdecomp, noisevarchrom, madL, variCb, edge, noisevarab_r, true, false, false, numThreads);
                    }

                    delete[] noisevarchrom;

                }
            }

            if (!Ldecomp.memory_allocation_failed()) {
                Lin = new array2D<float>(GW, GH);
#ifdef _OPENMP
                #pragma omp parallel for if (multiThread)
#endif
                for (int i = 0; i < GH; ++i) {
                    for (int j = 0; j < GW; ++j) {
                        (*Lin)[i][j] = tmp1.L[i][j];
                    }
                }

                Ldecomp.reconstruct(tmp1.L[0]);
            }

            if (!Ldecomp.memory_allocation_failed() && aut == 0) {
                if ((lp.noiself >= 0.01f ||  lp.noiself0 >= 0.01f ||  lp.noiself2 >= 0.01f || lp.wavcurvedenoi || lp.noiselc >= 0.01f) && levred == 7 && lp.noiseldetail != 100.f  && lp.quamet < 2) {
                    fftw_denoise(sk, GW, GH, max_numblox_W, min_numblox_W, tmp1.L, Lin,  numThreads, lp, 0);
                }
            }

            if (!adecomp.memory_allocation_failed()) {
                Ain = new array2D<float>(GW, GH);
#ifdef _OPENMP
                #pragma omp parallel for if (multiThread)
#endif
                for (int i = 0; i < GH; ++i) {
                    for (int j = 0; j < GW; ++j) {
                        (*Ain)[i][j] = tmp1.a[i][j];
                    }
                }

                adecomp.reconstruct(tmp1.a[0]);
            }


            if (!adecomp.memory_allocation_failed() && aut == 0) {
                if ((lp.noisecf >= 0.01f ||  lp.noisecc >= 0.01f) && levred == 7 && lp.noisechrodetail != 100.f) {
                    fftw_denoise(sk, GW, GH, max_numblox_W, min_numblox_W, tmp1.a, Ain,  numThreads, lp, 1);
                }
            }


            if (!bdecomp.memory_allocation_failed()) {

                Bin = new array2D<float>(GW, GH);
#ifdef _OPENMP
                #pragma omp parallel for if (multiThread)
#endif
                for (int i = 0; i < GH; ++i) {
                    for (int j = 0; j < GW; ++j) {
                        (*Bin)[i][j] = tmp1.b[i][j];
                    }
                }

                bdecomp.reconstruct(tmp1.b[0]);
            }


            if (!bdecomp.memory_allocation_failed() && aut == 0) {
                if ((lp.noisecf >= 0.01f ||  lp.noisecc >= 0.01f) && levred == 7 && lp.noisechrodetail != 100.f) {
                    fftw_denoise(sk, GW, GH, max_numblox_W, min_numblox_W, tmp1.b, Bin,  numThreads, lp, 1);
                }

            }
            if(gamma > 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                for (int y = 0; y < GH; ++y) {//apply inverse gamma 3.f and put result in range 32768.f
                    int x = 0;
#ifdef __SSE2__
                    for (; x < GW - 3; x += 4) {
                        STVFU(tmp1.L[y][x], F2V(32768.f) * gammalog(LVFU(tmp1.L[y][x]) / F2V(32768.f), F2V(gamma), F2V(ts), F2V(g_a[3]), F2V(g_a[4])));
                    }
#endif
                    for (; x < GW; ++x) {
                        tmp1.L[y][x] = 32768.f * gammalog(tmp1.L[y][x] / 32768.f, gamma, ts, g_a[3], g_a[4]);
                    }
                }
            }

            if(lp.smasktyp != 0) {
                if(lp.enablMask && lp.recothrd != 1.f) {
                    LabImage tmp3(GW, GH);

                    for (int ir = 0; ir < GH; ir++){
                        for (int jr = 0; jr < GW; jr++) {
                            tmp3.L[ir][jr] = original->L[ir][jr];
                            tmp3.a[ir][jr] = original->a[ir][jr];
                            tmp3.b[ir][jr] = original->b[ir][jr];
                        }
                    }

                    array2D<float> masklum(GW, GH);
                    array2D<float> masklumch(GW, GH);

                    float hig = lp.higthrd;
                    float higc;
                    calcdif(hig, higc);
                    float low = lp.lowthrd;
                    float lowc;
                    calcdif(low, lowc);
                    float mid = 0.01f * lp.midthrd;
                    float midch = 0.01f * lp.midthrdch;

                    if(higc < lowc) {
                        higc = lowc + 0.01f;
                    }
                    float th = (lp.recothrd - 1.f);
                    float ahigh = th / (higc - 100.f);
                    float bhigh = 1.f - higc * ahigh;

                    float alow = th / lowc; 
                    float blow = 1.f - th;
#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                    for (int ir = 0; ir < GH; ir++) {
                        for (int jr = 0; jr < GW; jr++) {
                            const float lmr = bufmaskblurbl->L[ir][jr] / 327.68f;
                            float k;
                            float kch;
                            if (lmr < lowc) {
                                k = alow * lmr + blow;
                                kch = alow * lmr + blow;
                            } else if (lmr < higc) {
                                k = 1.f - mid;
                                kch = 1.f - midch;
                            } else {
                                k = ahigh * lmr + bhigh;
                                kch = ahigh * lmr + bhigh;
                            }
                            if(lp.invmaskd) {
                                masklum[ir][jr] = 1.f - pow_F(k, lp.decayd);
                                masklumch[ir][jr] = 1.f - pow_F(kch, lp.decayd);
                            } else {
                                masklum[ir][jr] = pow_F(k, lp.decayd);
                                masklumch[ir][jr] = pow_F(kch, lp.decayd);
                            }
                        }
                    }

                    for (int i = 0; i < 3; ++i) {
                        boxblur(static_cast<float**>(masklum), static_cast<float**>(masklum), 10 / sk, GW, GH, multiThread);
                        boxblur(static_cast<float**>(masklumch), static_cast<float**>(masklumch), 10 / sk, GW, GH, multiThread);
                    }
#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                    for (int i = 0; i < GH; ++i) {
                        for (int j = 0; j < GW; ++j) {                              
                            tmp1.L[i][j] = (tmp3.L[i][j] - tmp1.L[i][j]) *  LIM01(masklum[i][j]) + tmp1.L[i][j];
                            tmp1.a[i][j] = (tmp3.a[i][j] - tmp1.a[i][j]) *  LIM01(masklumch[i][j]) + tmp1.a[i][j];
                            tmp1.b[i][j] = (tmp3.b[i][j] - tmp1.b[i][j]) *  LIM01(masklumch[i][j]) + tmp1.b[i][j];
                        }
                    }
                    masklum.free();
                    masklumch.free();
                }
                
// re read wavelet decomposition to calaculate noise residual
            float chresid = 0.f;
            float chresidtemp = 0.f;
            float chmaxresid = 0.f;
            float chmaxresidtemp = 0.f;
            float chresid46 = 0.f;
            float chresidtemp46 = 0.f;
            float chmaxresid46 = 0.f;
            float chmaxresidtemp46 = 0.f;
            float Lresid = 0.f;
            float Lmaxresid = 0.f;
            float Lresid46 = 0.f;
            float Lmaxresid46 = 0.f;
            
            
//calculate and display residual noise luma and chroma
// various coefficient from  1 to 5 - tries to take into account the difference between calculate Noise and percepted noise
            wavelet_decomposition Ldecompinf(tmp1.L[0], tmp1.W, tmp1.H, levwavL, 1, skip, numThreads, lp.daubLen);
            wavelet_decomposition adecompinf(tmp1.a[0], tmp1.W, tmp1.H, levwavL, 1, skip, numThreads, lp.daubLen);
            wavelet_decomposition bdecompinf(tmp1.b[0], tmp1.W, tmp1.H, levwavL, 1, skip, numThreads, lp.daubLen);

            Noise_residualAB(adecompinf, chresid, chmaxresid, false, 0, 3);
            chresidtemp = chresid;
            chmaxresidtemp = chmaxresid;
            Noise_residualAB(bdecompinf, chresid, chmaxresid, false, 0, 3);
            chresid += chresidtemp;
            chmaxresid += chmaxresidtemp;
            int nbmaddir = 4;
            chresid = sqrt(chresid / ( 3 * nbmaddir * 2));
            highresi = chresid + 0.5f * (sqrt(chmaxresid) - chresid); //evaluate sigma
            nresi = chresid;
            highresi /= 1.4f;//arbitrary coefficient
            nresi /= 1.4f;

    //        printf("nresi03=%f highresi=%f \n", (double) nresi, (double) highresi);


            Noise_residualAB(adecompinf, chresid46, chmaxresid46, false, 4, 6);
            nbmaddir = 3;
            chresidtemp46 = chresid46;
            chmaxresidtemp46 = chmaxresid46;
            Noise_residualAB(bdecompinf, chresid46, chmaxresid46, false, 4, 6);
            chresid46 += chresidtemp46;
            chmaxresid46 += chmaxresidtemp46;
            chresid46 = sqrt(chresid46 / ( 3 * nbmaddir * 2));
            highresi46 = chresid46 + 0.5f * (sqrt(chmaxresid46) - chresid46); //evaluate sigma
            nresi46 = chresid46;
            highresi46 /= 2.f;//arbitrary coefficient
            nresi46 /= 2.f;
            
    //        printf("nresi46=%f highresi=%f \n", (double) nresi46, (double) highresi46);


            Noise_residualAB(Ldecompinf, Lresid, Lmaxresid, false, 0, 3);
            nbmaddir = 4;
            Lresid = sqrt(Lresid / (3 * nbmaddir));
            Lhighresi = Lresid + 0.5f * (sqrt(Lmaxresid) - Lresid); //evaluate sigma
            Lnresi = Lresid;
            Lnresi /= 2.f;//arbitrary coefficient
            Lhighresi /= 2.f;
           // printf("Lresi03=%f Lhighresi=%f levwavL=%i\n", (double) Lnresi, (double) Lhighresi, levwavL);

            Noise_residualAB(Ldecompinf, Lresid46, Lmaxresid46, false, 4, 6);
            nbmaddir = 3;
            Lresid46 = sqrt(Lresid46 / (3 * nbmaddir));
            Lhighresi46 = Lresid46 + 0.5f * (sqrt(Lmaxresid46) - Lresid46); //evaluate sigma
            Lnresi46 = Lresid46;
            Lhighresi46 /= 5.f;//arbitrary coefficient
            Lnresi46 /= 5.f;
           // printf("Lresi46=%f Lhighresi=%f levwavL=%i\n", (double) Lnresi46, (double) Lhighresi46, levwavL);

// end calculate
                
            DeNoise_Local(call, lp,  originalmaskbl, levred, huerefblur, lumarefblur, chromarefblur, original, transformed, tmp1, cx, cy, sk);

            } else {
                DeNoise_Local(call, lp,  original, levred, huerefblur, lumarefblur, chromarefblur, original, transformed, tmp1, cx, cy, sk);
            }
        } else if (call == 2) { //simpleprocess

            const int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
            const int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
            const int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
            const int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
            const int bfh = yend - ystart;
            const int bfw = xend - xstart;

            if (bfh >= mDEN && bfw >= mDEN) {
                LabImage bufwv(bfw, bfh);
                bufwv.clear(true);
                array2D<float> *Lin = nullptr;
                array2D<float> *Ain = nullptr;
                array2D<float> *Bin = nullptr;

                int max_numblox_W = ceil((static_cast<float>(bfw)) / offset) + 2;
                // calculate min size of numblox_W.
                int min_numblox_W = ceil((static_cast<float>(bfw)) / offset) + 2;
                // these are needed only for creation of the plans and will be freed before entering the parallel loop


                int begy = ystart; //lp.yc - lp.lyT;
                int begx = xstart; //lp.xc - lp.lxL;
                int yEn = yend; //lp.yc + lp.ly;
                int xEn = xend; //lp.xc + lp.lx;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                for (int y = 0; y < transformed->H ; y++) //{
                    for (int x = 0; x < transformed->W; x++) {
                        int lox = cx + x;
                        int loy = cy + y;

                        if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                            bufwv.L[loy - begy][lox - begx] = original->L[y][x];
                            bufwv.a[loy - begy][lox - begx] = original->a[y][x];
                            bufwv.b[loy - begy][lox - begx] = original->b[y][x];
                        }

                    }

                float gamma = lp.noisegam;
                rtengine::GammaValues g_a; //gamma parameters
                double pwr = 1.0 / (double) lp.noisegam;//default 3.0 - gamma Lab
                double ts = 9.03296;//always the same 'slope' in the extreme shadows - slope Lab
                rtengine::Color::calcGamma(pwr, ts, g_a); // call to calcGamma with selected gamma and slope
                if(gamma > 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                    for (int y = 0; y < bfh; ++y) {
                        int x = 0;
                
#ifdef __SSE2__
                        for (; x <  bfw - 3; x += 4) {
                            STVFU(bufwv.L[y][x], F2V(32768.f) * igammalog(LVFU(bufwv.L[y][x]) / F2V(32768.f), F2V(gamma), F2V(ts), F2V(g_a[2]), F2V(g_a[4])));
                        }
#endif
                        for (;x < bfw; ++x) {
                             bufwv.L[y][x] = 32768.f * igammalog(bufwv.L[y][x] / 32768.f, gamma, ts, g_a[2], g_a[4]);
                        }
                    }
                }

                //   int DaubLen = 6;

                int levwavL = levred;
                int skip = 1;
                wavelet_decomposition Ldecomp(bufwv.L[0], bufwv.W, bufwv.H, levwavL, 1, skip, numThreads, lp.daubLen);
                wavelet_decomposition adecomp(bufwv.a[0], bufwv.W, bufwv.H, levwavL, 1, skip, numThreads, lp.daubLen);
                wavelet_decomposition bdecomp(bufwv.b[0], bufwv.W, bufwv.H, levwavL, 1, skip, numThreads, lp.daubLen);
                float madL[10][3];
                int edge = 2;

                if (!Ldecomp.memory_allocation_failed()) {
#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic) collapse(2) if (multiThread)
#endif
                    for (int lvl = 0; lvl < levred; lvl++) {
                        for (int dir = 1; dir < 4; dir++) {
                            int Wlvl_L = Ldecomp.level_W(lvl);
                            int Hlvl_L = Ldecomp.level_H(lvl);

                            const float* const* WavCoeffs_L = Ldecomp.level_coeffs(lvl);

                            madL[lvl][dir - 1] = SQR(Mad(WavCoeffs_L[dir], Wlvl_L * Hlvl_L));
                        }
                    }

                    float vari[levred];
                    float mxsl = 0.f;

                        {
                            edge = 2;
                            vari[0] = 0.8f * SQR((lp.noiself0 / 125.f) * (1.f + lp.noiself0 / 25.f));
                            vari[1] = 0.8f * SQR((lp.noiself / 125.f) * (1.f + lp.noiself / 25.f));
                            vari[2] = 0.8f * SQR((lp.noiself2 / 125.f) * (1.f + lp.noiself2 / 25.f));
                            vari[3] = 0.8f * SQR((lp.noiselc / 125.f) * (1.f + lp.noiselc / 25.f));
                            vari[4] = 1.f * SQR((lp.noiselc4 / 125.f) * (1.f + lp.noiselc4 / 25.f));
                            vari[5] = 1.5f * SQR((lp.noiselc5 / 125.f) * (1.f + lp.noiselc5 / 25.f));
                            vari[6] = 2.5f * SQR((lp.noiselc6 / 125.f) * (1.f + lp.noiselc6 / 25.f));
                        } 

                    {
                        float kr3 = 0.f;

                        {
                            if (lp.noiselc < 30.f) {
                                kr3 = 0.f;
                            } else if (lp.noiselc < 50.f) {
                                kr3 = 0.5f;
                            } else if (lp.noiselc < 70.f) {
                                kr3 = 0.7f;
                            } else {
                                kr3 = 1.f;
                            }
                        } 

                        vari[0] = rtengine::max(0.000001f, vari[0]);
                        vari[1] = rtengine::max(0.000001f, vari[1]);
                        vari[2] = rtengine::max(0.000001f, vari[2]);
                        vari[3] = rtengine::max(0.000001f, kr3 * vari[3]);
                        vari[4] = rtengine::max(0.000001f, vari[4]);
                        vari[5] = rtengine::max(0.000001f, vari[5]);
                        vari[6] = rtengine::max(0.000001f, vari[6]);

                        //    float* noisevarlum = nullptr;  // we need a dummy to pass it to WaveletDenoiseAllL
                        float* noisevarlum = new float[bfh * bfw];
                        float* noisevarhue = new float[bfh * bfw];
                        int bfw2 = (bfw + 1) / 2;

                        float nvlh[13] = {1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 0.7f, 0.5f}; //high value
                        float nvll[13] = {0.1f, 0.15f, 0.2f, 0.25f, 0.3f, 0.35f, 0.4f, 0.45f, 0.7f, 0.8f, 1.f, 1.f, 1.f}; //low value

                        float seuillow = 3000.f;//low
                        float seuilhigh = 18000.f;//high
                        int i = 10 - lp.noiselequal;
                        float ac = (nvlh[i] - nvll[i]) / (seuillow - seuilhigh);
                        float bc = nvlh[i] - seuillow * ac;
                        //ac and bc for transition
#ifdef _OPENMP
                        #pragma omp parallel for if (multiThread)
#endif
                        for (int ir = 0; ir < bfh; ir++)
                            for (int jr = 0; jr < bfw; jr++) {
                                float lN = bufwv.L[ir][jr];

                                if (lN < seuillow) {
                                    noisevarlum[(ir >> 1)*bfw2 + (jr >> 1)] =  nvlh[i];
                                } else if (lN < seuilhigh) {
                                    noisevarlum[(ir >> 1)*bfw2 + (jr >> 1)] = ac * lN + bc;
                                } else {
                                    noisevarlum[(ir >> 1)*bfw2 + (jr >> 1)] =  nvll[i];
                                }
                            }
                            
                    if(lp.enablMask && lp.lnoiselow != 1.f  && lp.smasktyp != 0) {
                         //this code has been reviewed by Ingo in september 2020 PR5903
                         //i just change parameters to better progressivity
                        float higc;
                        float hig = lp.thrhigh;
                        calcdif(hig, higc);
                        float low = lp.thrlow;
                        float lowc;
                        calcdif(low, lowc);
                         
                        if(higc < lowc) {
                            higc = lowc + 0.01f;
                        }

                        float alow = -(lp.lnoiselow - 1.f) / lowc;
                        float blow = lp.lnoiselow;
                        float ahigh = 0.9999f / (higc - 100.f);
                        float bhigh = 1.f - higc * ahigh;


#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                        for (int ir = 0; ir < bfh; ir++)
                            for (int jr = 0; jr < bfw; jr++) {
                                const float lM = bufmaskblurbl->L[ir + ystart][jr + xstart];
                                const float lmr = lM / 327.68f;
                                if (lM < 327.68f * lowc) {
                                    noisevarlum[(ir >> 1) * bfw2 + (jr >> 1)] *= alow * lmr + blow; 
                                } else if (lM < 327.68f * higc) {
                                    // do nothing
                                } else {
                                    noisevarlum[(ir >> 1) * bfw2 + (jr >> 1)] *= ahigh * lmr + bhigh;
                                }
                        }
                    }


                        if(HHhuecurve) {
                            //same code as in wavelet levels
#ifdef _OPENMP
        #pragma omp parallel for
#endif
                            for (int ir = 0; ir < bfh; ir++)
                                for (int jr = 0; jr < bfw; jr++) {
                                    float hueG = xatan2f(bufwv.b[ir][jr], bufwv.a[ir][jr]);
                                    float valparam = 2.f * (locwavCurvehue[500.f * static_cast<float>(Color::huelab_to_huehsv2(hueG))] - 0.5f);  //get H=f(H)
                                    noisevarhue[(ir >> 1)* bfw2 + (jr >> 1)] = 1.f +  valparam;
                                    noisevarlum[(ir >> 1)* bfw2 + (jr >> 1)] *= noisevarhue[(ir >> 1)* bfw2 + (jr >> 1)];
                                }
                        }


                        if ((lp.quamet == 0  && aut == 0) || (mxsl < 1.f && (aut == 1 || aut == 2))) {
                            WaveletDenoiseAllL(Ldecomp, noisevarlum, madL, vari, edge, numThreads);
                        } else if (lp.quamet == 1) {
                            WaveletDenoiseAll_BiShrinkL(Ldecomp, noisevarlum, madL, vari, edge, numThreads);
                            WaveletDenoiseAllL(Ldecomp, noisevarlum, madL, vari, edge, numThreads);
                        }

                        delete [] noisevarlum;
                        delete [] noisevarhue;

                    }
                }


                float variC[levred];
                float variCb[levred];

                float noisecfr = lp.noisecf;
                float noiseccr = lp.noisecc;

                if (lp.adjch > 0.f) {
                    noisecfr = lp.noisecf + 0.1f * lp.adjch;
                    noiseccr = lp.noisecc + 0.1f * lp.adjch;
                }

                float noisecfb = lp.noisecf;
                float noiseccb = lp.noisecc;

                if (lp.adjch < 0.f) {
                    noisecfb = lp.noisecf - 0.1f * lp.adjch;
                    noiseccb = lp.noisecc - 0.1f * lp.adjch;
                }


                if (noisecfr < 0.f) {
                    noisecfr = 0.00001f;
                }

                if (noiseccr < 0.f) {
                    noiseccr = 0.00001f;
                }

                if (noisecfb < 0.f) {
                    noisecfb = 0.00001f;
                }

                if (noiseccb < 0.f) {
                    noiseccb = 0.00001f;
                }


                if (!adecomp.memory_allocation_failed() && !bdecomp.memory_allocation_failed()) {
                    float maxccoarse = 0.f;


                        {
                            edge = 2;
                            variC[0] = SQR(noisecfr);
                            variC[1] = SQR(noisecfr);
                            variC[2] = SQR(noisecfr);
                            variC[3] = SQR(1.2f * noisecfr);
                            variC[4] = SQR(noisecfr);
                            variC[5] = SQR(1.2f * noiseccr);
                            variC[6] = SQR(1.5f * noiseccr);

                            variCb[0] = SQR(noisecfb);
                            variCb[1] = SQR(noisecfb);
                            variCb[2] = SQR(noisecfb);
                            variCb[3] = SQR(noisecfb);
                            variCb[4] = SQR(noisecfb);
                            variCb[5] = SQR(1.2f * noiseccb);
                            variCb[6] = SQR(1.5f * noiseccb);

                        } 

                    {
                        float minic = 0.000001f;

                        if (noiscfactiv) {
                            minic = 0.1f;//only for artifact shape detection
                        }

                        float k1 = 0.f;
                        float k2 = 0.f;
                        float k3 = 0.f;

                            if (lp.noisecf < 0.2f) {
                                k1 = 0.05f;
                                k2 = 0.f;
                                k3 = 0.f;
                            } else if (lp.noisecf < 0.3f) {
                                k1 = 0.1f;
                                k2 = 0.0f;
                                k3 = 0.f;
                            } else if (lp.noisecf < 0.5f) {
                                k1 = 0.2f;
                                k2 = 0.1f;
                                k3 = 0.f;
                            } else if (lp.noisecf < 0.8f) {
                                k1 = 0.3f;
                                k2 = 0.25f;
                                k3 = 0.f;
                            } else if (lp.noisecf < 1.f) {
                                k1 = 0.4f;
                                k2 = 0.25f;
                                k3 = 0.1f;
                            } else if (lp.noisecf < 2.f) {
                                k1 = 0.5f;
                                k2 = 0.3f;
                                k3 = 0.15f;
                            } else if (lp.noisecf < 3.f) {
                                k1 = 0.6f;
                                k2 = 0.45f;
                                k3 = 0.3f;
                            } else if (lp.noisecf < 4.f) {
                                k1 = 0.7f;
                                k2 = 0.5f;
                                k3 = 0.4f;
                            } else if (lp.noisecf < 5.f) {
                                k1 = 0.8f;
                                k2 = 0.6f;
                                k3 = 0.5f;
                            } else if (lp.noisecf < 6.f) {
                                k1 = 0.85f;
                                k2 = 0.7f;
                                k3 = 0.6f;
                            } else if (lp.noisecf < 8.f) {
                                k1 = 0.9f;
                                k2 = 0.8f;
                                k3 = 0.7f;
                            } else if (lp.noisecf < 10.f) {
                                k1 = 1.f;
                                k2 = 1.f;
                                k3 = 0.9f;

                            } else {
                                k1 = 1.f;
                                k2 = 1.f;
                                k3 = 1.5f;
                            }

                        variC[0] = rtengine::max(minic, variC[0]);
                        variC[1] = rtengine::max(minic, k1 * variC[1]);
                        variC[2] = rtengine::max(minic, k2 * variC[2]);
                        variC[3] = rtengine::max(minic, k3 * variC[3]);

                        variCb[0] = rtengine::max(minic, variCb[0]);
                        variCb[1] = rtengine::max(minic, k1 * variCb[1]);
                        variCb[2] = rtengine::max(minic, k2 * variCb[2]);
                        variCb[3] = rtengine::max(minic, k3 * variCb[3]);

                        {
                            float k4 = 0.f;
                            float k5 = 0.f;
                            float k6 = 0.f;

                            if (lp.noisecc < 0.2f) {
                                k4 = 0.1f;
                                k5 = 0.02f;
                            } else if (lp.noisecc < 0.5f) {
                                k4 = 0.15f;
                                k5 = 0.05f;
                            } else if (lp.noisecc < 1.f) {
                                k4 = 0.15f;
                                k5 = 0.1f;
                            } else if (lp.noisecc < 3.f) {
                                k4 = 0.3f;
                                k5 = 0.15f;
                            } else if (lp.noisecc < 4.f) {
                                k4 = 0.6f;
                                k5 = 0.4f;
                            } else if (lp.noisecc < 6.f) {
                                k4 = 0.8f;
                                k5 = 0.6f;
                            } else {
                                k4 = 1.f;
                                k5 = 1.f;
                            }


                            variC[4] = rtengine::max(0.000001f, k4 * variC[4]);
                            variC[5] = rtengine::max(0.000001f, k5 * variC[5]);
                            variCb[4] = rtengine::max(0.000001f, k4 * variCb[4]);
                            variCb[5] = rtengine::max(0.000001f, k5 * variCb[5]);

                            if (lp.noisecc < 4.f) {
                                k6 = 0.f;
                            } else if (lp.noisecc < 5.f) {
                                k6 = 0.4f;
                            } else if (lp.noisecc < 6.f) {
                                k6 = 0.7f;
                            } else {
                                k6 = 1.f;
                            }

                            variC[6] = rtengine::max(0.00001f, k6 * variC[6]);
                            variCb[6] = rtengine::max(0.00001f, k6 * variCb[6]);
                        }

                        float* noisevarchrom = new float[bfh * bfw];
                        int bfw2 = (bfw + 1) / 2;
                        float nvch = 0.6f;//high value
                        float nvcl = 0.1f;//low value

                        if (lp.noisecf > 30.f) {
                            nvch = 0.8f;
                            nvcl = 0.4f;
                        }

                        float seuil = 4000.f;//low
                        float seuil2 = 15000.f;//high
                        //ac and bc for transition
                        float ac = (nvch - nvcl) / (seuil - seuil2);
                        float bc = nvch - seuil * ac;
#ifdef _OPENMP
                        #pragma omp parallel for if (multiThread)
#endif
                        for (int ir = 0; ir < bfh; ir++)
                            for (int jr = 0; jr < bfw; jr++) {
                                float cN = std::sqrt(SQR(bufwv.a[ir][jr]) + SQR(bufwv.b[ir][jr]));

                                if (cN < seuil) {
                                    noisevarchrom[(ir >> 1)*bfw2 + (jr >> 1)] = nvch;
                                } else if (cN < seuil2) {
                                    noisevarchrom[(ir >> 1)*bfw2 + (jr >> 1)] = ac * cN + bc;
                                } else {
                                    noisevarchrom[(ir >> 1)*bfw2 + (jr >> 1)] = nvcl;
                                }
                            }

                        float noisevarab_r = 100.f; //SQR(lp.noisecc / 10.0);

                        if ((lp.quamet == 0  && aut == 0) || (maxccoarse < 0.1f && (aut == 1  || aut == 2)))  {
                            WaveletDenoiseAllAB(Ldecomp, adecomp, noisevarchrom, madL, variC, edge, noisevarab_r, true, false, false, numThreads);
                            WaveletDenoiseAllAB(Ldecomp, bdecomp, noisevarchrom, madL, variCb, edge, noisevarab_r, true, false, false, numThreads);
                        } else if (lp.quamet == 1){
                            WaveletDenoiseAll_BiShrinkAB(Ldecomp, adecomp, noisevarchrom, madL, variC, edge, noisevarab_r, true, false, false, numThreads);
                            WaveletDenoiseAllAB(Ldecomp, adecomp, noisevarchrom, madL, variC, edge, noisevarab_r, true, false, false, numThreads);

                            WaveletDenoiseAll_BiShrinkAB(Ldecomp, bdecomp, noisevarchrom, madL, variCb, edge, noisevarab_r, true, false, false, numThreads);
                            WaveletDenoiseAllAB(Ldecomp, bdecomp, noisevarchrom, madL, variCb, edge, noisevarab_r, true, false, false, numThreads);
                        }

                        delete[] noisevarchrom;
                    }
                }

                if (!Ldecomp.memory_allocation_failed()) {
                    Lin = new array2D<float>(bfw, bfh);

#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                    for (int i = 0; i < bfh; ++i) {
                        for (int j = 0; j < bfw; ++j) {
                            (*Lin)[i][j] = bufwv.L[i][j];
                        }
                    }

                    Ldecomp.reconstruct(bufwv.L[0]);
                }


                if (!Ldecomp.memory_allocation_failed() && aut == 0) {


                    if ((lp.noiself >= 0.01f ||  lp.noiself0 >= 0.01f ||  lp.noiself2 >= 0.01f || lp.wavcurvedenoi || lp.noiselc >= 0.01f) && levred == 7 && lp.noiseldetail != 100.f && lp.quamet < 2) {
                        fftw_denoise(sk, bfw, bfh, max_numblox_W, min_numblox_W, bufwv.L, Lin,  numThreads, lp, 0);
                    }
                }


                if (!adecomp.memory_allocation_failed()) {
                    Ain = new array2D<float>(bfw, bfh);
#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                    for (int i = 0; i < bfh; ++i) {
                        for (int j = 0; j < bfw; ++j) {
                            (*Ain)[i][j] = bufwv.a[i][j];
                        }
                    }

                    adecomp.reconstruct(bufwv.a[0]);
                }

                if (!adecomp.memory_allocation_failed() && aut == 0) {
                    if ((lp.noisecf >= 0.001f ||  lp.noisecc >= 0.001f) && levred == 7 && lp.noisechrodetail != 100.f) {
                        fftw_denoise(sk, bfw, bfh, max_numblox_W, min_numblox_W, bufwv.a, Ain,  numThreads, lp, 1);
                    }
                }


                if (!bdecomp.memory_allocation_failed()) {
                    Bin = new array2D<float>(bfw, bfh);
#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                    for (int i = 0; i < bfh; ++i) {
                        for (int j = 0; j < bfw; ++j) {
                            (*Bin)[i][j] = bufwv.b[i][j];
                        }
                    }

                    bdecomp.reconstruct(bufwv.b[0]);
                }

                if (!bdecomp.memory_allocation_failed() && aut == 0) {
                    if ((lp.noisecf >= 0.001f ||  lp.noisecc >= 0.001f) && levred == 7 && lp.noisechrodetail != 100.f) {
                        fftw_denoise(sk, bfw, bfh, max_numblox_W, min_numblox_W, bufwv.b, Bin,  numThreads, lp, 1);
                    }
                }

                if(gamma > 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                    for (int y = 0; y < bfh ; ++y) {//apply inverse gamma 3.f and put result in range 32768.f
                        int x = 0;
        
#ifdef __SSE2__
                        for (; x < bfw  - 3; x += 4) {

                            STVFU(bufwv.L[y][x], F2V(32768.f) * gammalog(LVFU(bufwv.L[y][x]) / F2V(32768.f), F2V(gamma), F2V(ts), F2V(g_a[3]), F2V(g_a[4])));
                        }
#endif
                        for (; x < bfw ; ++x) {

                            bufwv.L[y][x] = 32768.f * gammalog(bufwv.L[y][x] / 32768.f, gamma, ts, g_a[3], g_a[4]);
                        }
                    }
                }

                if(lp.nlstr > 0) {
                    NLMeans(bufwv.L, lp.nlstr, lp.nldet, lp.nlpat, lp.nlrad, lp.nlgam, bfw, bfh, 1.f, multiThread);
                }


                if (lp.smasktyp != 0) {
                    if(lp.enablMask && lp.recothrd != 1.f) {
                        LabImage tmp3(bfw, bfh);
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < transformed->H ; y++) {
                            for (int x = 0; x < transformed->W; x++) {
                                int lox = cx + x;
                                int loy = cy + y;

                                if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                                    tmp3.L[loy - begy][lox - begx] = original->L[y][x];
                                    tmp3.a[loy - begy][lox - begx] = original->a[y][x];
                                    tmp3.b[loy - begy][lox - begx] = original->b[y][x];
                                }   

                            }
                        }

                        array2D<float> masklum;
                        array2D<float> masklumch;
                        masklum(bfw, bfh);
                        masklumch(bfw, bfh);
                        for (int ir = 0; ir < bfh; ir++){
                            for (int jr = 0; jr < bfw; jr++) {
                                masklum[ir][jr] = 1.f;
                                masklumch[ir][jr] = 1.f;
                            }
                        }

                        float hig = lp.higthrd;
                        float higc;
                        calcdif(hig, higc);
                        float low = lp.lowthrd;
                        float lowc;
                        calcdif(low, lowc);
                        float mid = 0.01f * lp.midthrd;
                        float midch = 0.01f * lp.midthrdch;

                        if(higc < lowc) {
                            higc = lowc + 0.01f;
                        }
                        float th = (lp.recothrd - 1.f);
                        float ahigh = th / (higc - 100.f);
                        float bhigh = 1.f - higc * ahigh;

                        float alow = th / lowc; 
                        float blow = 1.f - th;

#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                        for (int y = ystart; y < yend; y++) {
                            for (int x = xstart, lox = cx + x; x < xend; x++, lox++) {

                                const float lM = bufmaskblurbl->L[y][x]; 
                                const float lmr = lM / 327.68f;
                                if (lM < 327.68f * lowc) {
                                    masklum[y-ystart][x-xstart] = alow * lmr + blow;
                                    masklumch[y-ystart][x-xstart] = alow * lmr + blow;
                                } else if (lM < 327.68f * higc) {
                                    masklum[y-ystart][x-xstart] = 1.f - mid;
                                    masklumch[y-ystart][x-xstart] = 1.f - midch;

                                } else {
                                    masklum[y-ystart][x-xstart] = ahigh * lmr + bhigh;
                                    masklumch[y-ystart][x-xstart] = ahigh * lmr + bhigh;
                                }
                                float k = masklum[y-ystart][x-xstart];
                                float kch = masklumch[y-ystart][x-xstart];
                                if(lp.invmaskd == true) {
                                    masklum[y-ystart][x-xstart] = 1.f - pow(k, lp.decayd);
                                    masklumch[y-ystart][x-xstart] = 1.f - pow(kch, lp.decayd);
                                } else {
                                    masklum[y-ystart][x-xstart] = pow(k, lp.decayd);
                                    masklumch[y-ystart][x-xstart] = pow(kch, lp.decayd);
                                }
                            }
                        }
                        for (int i = 0; i < 3; ++i) {
                            boxblur(static_cast<float**>(masklum), static_cast<float**>(masklum), 10 / sk, bfw, bfh, false);
                            boxblur(static_cast<float**>(masklumch), static_cast<float**>(masklumch), 10 / sk, bfw, bfh, false);
                        }
                           
#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                        for (int y = 0; y < bfh; y++) {
                            for (int x = 0; x < bfw; x++) {
                                bufwv.L[y][x] = (tmp3.L[y][x] - bufwv.L[y][x]) *  LIM01(masklum[y][x]) + bufwv.L[y][x];
                                bufwv.a[y][x] = (tmp3.a[y][x] - bufwv.a[y][x]) *  LIM01(masklumch[y][x]) + bufwv.a[y][x];
                                bufwv.b[y][x] = (tmp3.b[y][x] - bufwv.b[y][x]) *  LIM01(masklumch[y][x]) + bufwv.b[y][x];
                            }
                        }

                        masklum.free();
                        masklumch.free();
                    }

                    DeNoise_Local2(lp,  originalmaskbl, levred, huerefblur, lumarefblur, chromarefblur, original, transformed, bufwv, cx, cy, sk);
                } else {
                    DeNoise_Local2(lp,  original, levred, huerefblur, lumarefblur, chromarefblur, original, transformed, bufwv, cx, cy, sk);
                }
            }
        }
 
    }

}

float triangle(float a, float a1, float b)
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

void rgbtone(float& maxval, float& medval, float& minval, const LUTf& lutToneCurve)
{
    float minvalold = minval, medvalold = medval, maxvalold = maxval;

    maxval = lutToneCurve[maxvalold];
    minval = lutToneCurve[minvalold];
    medval = minval + ((maxval - minval) * (medvalold - minvalold) / (maxvalold - minvalold));
}

void ImProcFunctions::clarimerge(const struct local_params& lp, float &mL, float &mC, bool &exec, LabImage *tmpresid, int wavelet_level, int sk, int numThreads)
{
    if (mL != 0.f && mC == 0.f) {
        mC = 0.0001f;
        exec = true;
    }

    if (mC != 0.f && mL == 0.f) {
        mL = 0.0001f;
        exec = true;
    }

    if (mL != 0.f && mC != 0.f) {
        exec = true;
    }

    if (mL != 0.f) {

        wavelet_decomposition *wdspotresid = new wavelet_decomposition(tmpresid->L[0], tmpresid->W, tmpresid->H, wavelet_level, 1, sk, numThreads, lp.daubLen);

        if (wdspotresid->memory_allocation_failed()) {
            return;
        }

        int maxlvlresid = wdspotresid->maxlevel();

        if (maxlvlresid > 4) {//Clarity
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic) collapse(2)
#endif
            for (int dir = 1; dir < 4; dir++) {
                for (int level = 0; level < maxlvlresid; ++level) {
                    int W_L = wdspotresid->level_W(level);
                    int H_L = wdspotresid->level_H(level);
                    float* const* wav_Lresid = wdspotresid->level_coeffs(level);

                    for (int i = 0; i < W_L * H_L; i++) {
                        wav_Lresid[dir][i] = 0.f;
                    }
                }
            }
        } else {//Sharp
            float *wav_L0resid = wdspotresid->get_coeff0();
            int W_L = wdspotresid->level_W(0);
            int H_L = wdspotresid->level_H(0);

            for (int i = 0; i < W_L * H_L; i++) {
                wav_L0resid[i] = 0.f;
            }
        }

        wdspotresid->reconstruct(tmpresid->L[0], 1.f);
        delete wdspotresid;
    }


    if (mC != 0.f) {

        wavelet_decomposition *wdspotresida = new wavelet_decomposition(tmpresid->a[0], tmpresid->W, tmpresid->H, wavelet_level, 1, sk, numThreads, lp.daubLen);

        if (wdspotresida->memory_allocation_failed()) {
            return;
        }

        int maxlvlresid = wdspotresida->maxlevel();

        if (maxlvlresid > 4) {//Clarity
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic) collapse(2)
#endif
            for (int dir = 1; dir < 4; dir++) {
                for (int level = 0; level < maxlvlresid; ++level) {
                    int W_L = wdspotresida->level_W(level);
                    int H_L = wdspotresida->level_H(level);
                    float* const* wav_Lresida = wdspotresida->level_coeffs(level);

                    for (int i = 0; i < W_L * H_L; i++) {
                        wav_Lresida[dir][i] = 0.f;
                    }
                }
            }
        } else {//Sharp
            float *wav_L0resida = wdspotresida->get_coeff0();
            int W_L = wdspotresida->level_W(0);
            int H_L = wdspotresida->level_H(0);

            for (int i = 0; i < W_L * H_L; i++) {
                wav_L0resida[i] = 0.f;
            }
        }

        wdspotresida->reconstruct(tmpresid->a[0], 1.f);
        delete wdspotresida;

        wavelet_decomposition *wdspotresidb = new wavelet_decomposition(tmpresid->b[0], tmpresid->W, tmpresid->H, wavelet_level, 1, sk, numThreads, lp.daubLen);

        if (wdspotresidb->memory_allocation_failed()) {
            return;
        }

        maxlvlresid = wdspotresidb->maxlevel();

        if (maxlvlresid > 4) {//Clarity
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic) collapse(2)
#endif
            for (int dir = 1; dir < 4; dir++) {
                for (int level = 0; level < maxlvlresid; ++level) {
                    int W_L = wdspotresidb->level_W(level);
                    int H_L = wdspotresidb->level_H(level);
                    float* const* wav_Lresidb = wdspotresidb->level_coeffs(level);

                    for (int i = 0; i < W_L * H_L; i++) {
                        wav_Lresidb[dir][i] = 0.f;
                    }
                }
            }
        } else {//Sharp
            float *wav_L0residb = wdspotresidb->get_coeff0();
            int W_L = wdspotresidb->level_W(0);
            int H_L = wdspotresidb->level_H(0);

            for (int i = 0; i < W_L * H_L; i++) {
                wav_L0residb[i] = 0.f;
            }
        }

        wdspotresidb->reconstruct(tmpresid->b[0], 1.f);
        delete wdspotresidb;
    }
}

void ImProcFunctions::avoidcolshi(const struct local_params& lp, int sp, LabImage *transformed, LabImage *reserved, int cy, int cx, int sk)
{
    int avoidgamut = 0;

    if (params->locallab.spots.at(sp).avoidgamutMethod == "NONE") {
        avoidgamut = 0;
    } else if (params->locallab.spots.at(sp).avoidgamutMethod == "LAB") {
        avoidgamut = 1;
    } else if (params->locallab.spots.at(sp).avoidgamutMethod == "XYZ") {
        avoidgamut = 2;
    } else if (params->locallab.spots.at(sp).avoidgamutMethod == "XYZREL") {
        avoidgamut = 3;
    } else if (params->locallab.spots.at(sp).avoidgamutMethod == "MUNS") {
        avoidgamut = 4;
    }

    if (avoidgamut == 0) {
        return;
    }

    if (avoidgamut > 0  && lp.islocal) {
        const float ach = lp.trans / 100.f;
        bool execmunsell = true;

        if (params->locallab.spots.at(sp).expcie && (params->locallab.spots.at(sp).modecam == "all" || params->locallab.spots.at(sp).modecam == "jz" || params->locallab.spots.at(sp).modecam == "cam16")) {
            execmunsell = false;
        }

        TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix(params->icm.workingProfile);
        const double wip[3][3] = {//improve precision with double
            {wiprof[0][0], wiprof[0][1], wiprof[0][2]},
            {wiprof[1][0], wiprof[1][1], wiprof[1][2]},
            {wiprof[2][0], wiprof[2][1], wiprof[2][2]}
        };

        TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
        const double wp[3][3] = {//improve precision with double
            {wprof[0][0], wprof[0][1], wprof[0][2]},
            {wprof[1][0], wprof[1][1], wprof[1][2]},
            {wprof[2][0], wprof[2][1], wprof[2][2]}
        };

        const float softr = params->locallab.spots.at(sp).avoidrad;//max softr = 30
        //   const bool muns = params->locallab.spots.at(sp).avoidmun;//Munsell control with 200 LUT
        //improve precision with mint and maxt
        const float tr = std::min(2.f, softr);
        const float mint = 0.15f - 0.06f * tr;//between 0.15f and 0.03f
        const float maxt = 0.98f + 0.008f * tr;//between 0.98f and 0.996f

        const bool highlight = params->toneCurve.hrenabled;
        const bool needHH =  true; //always Munsell to avoid bad behavior //(lp.chro != 0.f);
#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
#ifdef __SSE2__
            float atan2Buffer[transformed->W] ALIGNED16;
            float sqrtBuffer[transformed->W] ALIGNED16;
            float sincosyBuffer[transformed->W] ALIGNED16;
            float sincosxBuffer[transformed->W] ALIGNED16;
            vfloat c327d68v = F2V(327.68f);
            vfloat onev = F2V(1.f);
#endif

#ifdef _OPENMP
            #pragma omp for schedule(dynamic,16)
#endif

            for (int y = 0; y < transformed->H; y++) {
                const int loy = cy + y;
                const bool isZone0 = loy > lp.yc + lp.ly || loy < lp.yc - lp.lyT; // whole line is zone 0 => we can skip a lot of processing

                if (isZone0) { // outside selection and outside transition zone => no effect, keep original values
                    continue;
                }

#ifdef __SSE2__
                int i = 0;

                for (; i < transformed->W - 3; i += 4) {
                    vfloat av = LVFU(transformed->a[y][i]);
                    vfloat bv = LVFU(transformed->b[y][i]);

                    if (needHH) { // only do expensive atan2 calculation if needed
                        STVF(atan2Buffer[i], xatan2f(bv, av));
                    }

                    vfloat Chprov1v = vsqrtf(SQRV(bv) + SQRV(av));
                    STVF(sqrtBuffer[i], Chprov1v / c327d68v);
                    vfloat sincosyv = av / Chprov1v;
                    vfloat sincosxv = bv / Chprov1v;
                    vmask selmask = vmaskf_eq(Chprov1v, ZEROV);
                    sincosyv = vself(selmask, onev, sincosyv);
                    sincosxv = vselfnotzero(selmask, sincosxv);
                    STVF(sincosyBuffer[i], sincosyv);
                    STVF(sincosxBuffer[i], sincosxv);
                }

                for (; i < transformed->W; i++) {
                    float aa = transformed->a[y][i];
                    float bb = transformed->b[y][i];

                    if (needHH) { // only do expensive atan2 calculation if needed
                        atan2Buffer[i] = xatan2f(bb, aa);
                    }

                    float Chprov1 = std::sqrt(SQR(bb) + SQR(aa));
                    sqrtBuffer[i] = Chprov1 / 327.68f;

                    if (Chprov1 == 0.0f) {
                        sincosyBuffer[i] = 1.f;
                        sincosxBuffer[i] = 0.0f;
                    } else {
                        sincosyBuffer[i] = aa / Chprov1;
                        sincosxBuffer[i] = bb / Chprov1;
                    }
                }

#endif

                for (int x = 0; x < transformed->W; x++) {
                    int lox = cx + x;
                    int zone;
                    float localFactor = 1.f;

                    if (lp.shapmet == 0) {
                        calcTransition(lox, loy, ach, lp, zone, localFactor);
                    } else { /*if (lp.shapmet == 1)*/
                        calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                    }

                    if (zone == 0) { // outside selection and outside transition zone => no effect, keep original values
                        continue;
                    }

                    float Lprov1 = transformed->L[y][x] / 327.68f;
                    float2 sincosval;
#ifdef __SSE2__
                    float HH = atan2Buffer[x]; // reading HH from line buffer even if line buffer is not filled is faster than branching
                    float Chprov1 = sqrtBuffer[x];
                    sincosval.y = sincosyBuffer[x];
                    sincosval.x = sincosxBuffer[x];
                    float chr = 0.f;

#else
                    const float aa = transformed->a[y][x];
                    const float bb = transformed->b[y][x];
                    float HH = 0.f, chr = 0.f;

                    if (needHH) { // only do expensive atan2 calculation if needed
                        HH = xatan2f(bb, aa);
                    }

                    float Chprov1 = std::sqrt(SQR(aa) + SQR(bb)) / 327.68f;

                    if (Chprov1 == 0.0f) {
                        sincosval.y = 1.f;
                        sincosval.x = 0.0f;
                    } else {
                        sincosval.y = aa / (Chprov1 * 327.68f);
                        sincosval.x = bb / (Chprov1 * 327.68f);
                    }

#endif
                    float lnew = transformed->L[y][x];
                    float anew = transformed->a[y][x];
                    float bnew = transformed->b[y][x];
                    Lprov1 = lnew / 327.68f;
                    //HH = xatan2f(bnew, anew);

                    if (avoidgamut == 1) { //Lab correction

                        Color::pregamutlab(Lprov1, HH, chr);
                        Chprov1 = rtengine::min(Chprov1, chr);

                        float R, G, B;
                        Color::gamutLchonly(HH, sincosval, Lprov1, Chprov1, R, G, B, wip, highlight, mint, maxt);//replace for best results
                        lnew = Lprov1 * 327.68f;
                        anew = 327.68f * Chprov1 * sincosval.y;
                        bnew = 327.68f * Chprov1 * sincosval.x;
                        //HH = xatan2f(bnew, anew);
                        transformed->a[y][x] = anew;
                        transformed->b[y][x] = bnew;

                    } else if (avoidgamut == 2  || avoidgamut == 3) { //XYZ correction
                        float xg, yg, zg;
                        const float aag = transformed->a[y][x];//anew
                        const float bbg = transformed->b[y][x];//bnew
                        float Lag = transformed->L[y][x];

                        Color::Lab2XYZ(Lag, aag, bbg, xg, yg, zg);
                        float x0 = xg;
                        float y0 = yg;
                        float z0 = zg;

                        Color::gamutmap(xg, yg, zg, wp);

                        if (avoidgamut == 3) {//0.5f arbitrary coeff
                            xg = xg + 0.5f * (x0 - xg);
                            yg = yg + 0.5f * (y0 - yg);
                            zg = zg + 0.5f * (z0 - zg);
                        }

                        //Color::gamutmap(xg, yg, zg, wp);//Put XYZ in gamut wp
                        float aag2, bbg2;
                        Color::XYZ2Lab(xg, yg, zg, Lag, aag2, bbg2);
                        Lprov1 = Lag / 327.68f;
                        HH = xatan2f(bbg2, aag2);//rebuild HH in case of...absolute colorimetry
                        Chprov1 = std::sqrt(SQR(aag2) + SQR(bbg2)) / 327.68f;

                        if (Chprov1 == 0.0f) {
                            sincosval.y = 1.f;
                            sincosval.x = 0.0f;
                        } else {
                            sincosval.y = aag2 / (Chprov1 * 327.68f);
                            sincosval.x = bbg2 / (Chprov1 * 327.68f);
                        }

                        lnew = Lprov1 * 327.68f;
                        anew = 327.68f * Chprov1 * sincosval.y;
                        bnew = 327.68f * Chprov1 * sincosval.x;
                        transformed->a[y][x] = anew;
                        transformed->b[y][x] = bnew;

                    }

                    if (needHH && avoidgamut <= 4) {//Munsell
                        Lprov1 = lnew / 327.68f;
                        float Chprov = sqrt(SQR(anew) + SQR(bnew)) / 327.68f;

                        const float Lprov2 = reserved->L[y][x] / 327.68f;
                        float correctionHue = 0.f; // Munsell's correction
                        float correctlum = 0.f;
                        const float memChprov = std::sqrt(SQR(reserved->a[y][x]) + SQR(reserved->b[y][x])) / 327.68f;

                        if (execmunsell) {
                            Color::AllMunsellLch(true, Lprov1, Lprov2, HH, Chprov, memChprov, correctionHue, correctlum);
                        }

                        if (correctionHue != 0.f || correctlum != 0.f) {

                            if (std::fabs(correctionHue) < 0.015f) {
                                HH += correctlum;    // correct only if correct Munsell chroma very small.
                            }

                            sincosval = xsincosf(HH + correctionHue);
                        }

                        anew = 327.68f * Chprov * sincosval.y; // apply Munsell
                        bnew = 327.68f * Chprov * sincosval.x;
                        transformed->a[y][x] = anew; // apply Munsell
                        transformed->b[y][x] = bnew;
                    }
                }
            }
        }

        //Guidedfilter to reduce artifacts in transitions : case Lab
        if (softr != 0.f && avoidgamut == 1) {//soft for L a b because we change color...
            const float tmpblur = softr < 0.f ? -1.f / softr : 1.f + softr;
            const int r1 = rtengine::max<int>(6 / sk * tmpblur + 0.5f, 1);
            const int r2 = rtengine::max<int>(10 / sk * tmpblur + 0.5f, 1);

            constexpr float epsilmax = 0.005f;
            constexpr float epsilmin = 0.00001f;

            constexpr float aepsil = (epsilmax - epsilmin) / 100.f;
            constexpr float bepsil = epsilmin;
            const float epsil = softr < 0.f ? 0.001f : aepsil * softr + bepsil;

            const int bw = transformed->W;
            const int bh = transformed->H;
            array2D<float> ble(bw, bh);
            array2D<float> guid(bw, bh);

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif

            for (int y = 0; y < bh ; y++) {
                for (int x = 0; x < bw; x++) {
                    ble[y][x] = transformed->L[y][x] / 32768.f;
                    guid[y][x] = reserved->L[y][x] / 32768.f;
                }
            }

            rtengine::guidedFilter(guid, ble, ble, r2, 0.2f * epsil, multiThread);
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif

            for (int y = 0; y < bh; y++) {
                for (int x = 0; x < bw; x++) {
                    transformed->L[y][x] = 32768.f * ble[y][x];
                }
            }

            array2D<float> &blechro = ble; // reuse buffer
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif

            for (int y = 0; y < bh ; y++) {
                for (int x = 0; x < bw; x++) {
                    blechro[y][x] = std::sqrt(SQR(transformed->b[y][x]) + SQR(transformed->a[y][x])) / 32768.f;
                }
            }

            rtengine::guidedFilter(guid, blechro, blechro, r1, epsil, multiThread);

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif

            for (int y = 0; y < bh; y++) {
                for (int x = 0; x < bw; x++) {
                    const float Chprov1 = std::sqrt(SQR(transformed->a[y][x]) + SQR(transformed->b[y][x]));
                    float2  sincosval;

                    if (Chprov1 == 0.0f) {
                        sincosval.y = 1.f;
                        sincosval.x = 0.0f;
                    } else {
                        sincosval.y = transformed->a[y][x] / Chprov1;
                        sincosval.x = transformed->b[y][x] / Chprov1;
                    }

                    transformed->a[y][x] = 32768.f * blechro[y][x] * sincosval.y;
                    transformed->b[y][x] = 32768.f * blechro[y][x] * sincosval.x;
                }
            }
        }
    }
}

void maskrecov(const LabImage * bufcolfin, LabImage * original, LabImage * bufmaskblurcol, int bfh, int bfw, int ystart, int xstart, float hig, float low, float recoth, float decay, bool invmask, int sk, bool multiThread)
{
    LabImage tmp3(bfw, bfh);
    for (int y = 0; y < bfh; y++){
        for (int x = 0; x < bfw; x++) {
            tmp3.L[y][x] = original->L[y + ystart][x + xstart];
            tmp3.a[y][x] = original->a[y + ystart][x + xstart];
            tmp3.b[y][x] = original->b[y + ystart][x + xstart];
        }
    }
    array2D<float> masklum;
    masklum(bfw, bfh);
    for (int ir = 0; ir < bfh; ir++)
        for (int jr = 0; jr < bfw; jr++) {
            masklum[ir][jr] = 1.f;
        }

    float higc;
    calcdif(hig, higc);
    float lowc;
    calcdif(low, lowc);

    if(higc < lowc) {
        higc = lowc + 0.01f;
    }
    float th = (recoth - 1.f);
    float ahigh = th / (higc - 100.f);
    float bhigh = 1.f - higc * ahigh;

    float alow = th / lowc; 
    float blow = 1.f - th;
#ifdef _OPENMP
   #pragma omp parallel for if (multiThread)
#endif
    for (int ir = 0; ir < bfh; ir++) {
        for (int jr = 0; jr < bfw; jr++) {
            const float lM = bufmaskblurcol->L[ir][jr];
            const float lmr = lM / 327.68f;
            if (lM < 327.68f * lowc) {
                masklum[ir][jr] = alow * lmr + blow;
            } else if (lM < 327.68f * higc) {
                //nothing...but we can..
            } else {
                masklum[ir][jr] = ahigh * lmr + bhigh;
            }
            float k = masklum[ir][jr];
            if(invmask == false) {
                masklum[ir][jr] = 1 - pow(k, decay);
            } else {
                masklum[ir][jr] = pow(k, decay);
            }

        }
    }

        for (int i = 0; i < 3; ++i) {
            boxblur(static_cast<float**>(masklum), static_cast<float**>(masklum), 10 / sk, bfw, bfh, false);
        }

#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
        for (int i = 0; i < bfh; ++i) {
            for (int j = 0; j < bfw; ++j) {                              
                bufcolfin->L[i][j] = (tmp3.L[i][j] - bufcolfin->L[i][j]) *  LIM01(masklum[i][j]) + bufcolfin->L[i][j];
                bufcolfin->a[i][j] = (tmp3.a[i][j] - bufcolfin->a[i][j]) *  LIM01(masklum[i][j]) + bufcolfin->a[i][j];
                bufcolfin->b[i][j] = (tmp3.b[i][j] - bufcolfin->b[i][j]) *  LIM01(masklum[i][j]) + bufcolfin->b[i][j];
            }
        }
        masklum.free();
}

//thanks to Alberto Griggio
void ImProcFunctions::detail_mask(const array2D<float> &src, array2D<float> &mask, int bfw, int bfh, float scaling, float threshold, float ceiling, float factor, BlurType blur_type, float blur, bool multithread)
{
    const int W = bfw;
    const int H = bfh;
    mask(W, H);

    array2D<float> L2(W/4, H/4);//ARRAY2D_ALIGNED);
    array2D<float> m2(W/4, H/4);//ARRAY2D_ALIGNED)
    rescaleBilinear(src, L2, multithread);
#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H/4; ++y) {
        for (int x = 0; x < W/4; ++x) {
            L2[y][x] = xlin2log(L2[y][x]/scaling, 50.f);
        }
    }
    
    laplacian(L2, m2, W / 4, H / 4, threshold/scaling, ceiling/scaling, factor, multithread);

    rescaleBilinear(m2, mask, multithread);

    const auto scurve =
        [](float x) -> float
        {
            constexpr float b = 101.f;
            constexpr float a = 2.23f;
            return xlin2log(pow_F(x, a), b);
        };

    const float thr = 1.f - factor;
#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            mask[y][x] = scurve(LIM01(mask[y][x] + thr));
        }
    }

    if (blur_type == BlurType::GAUSS) {
        
#ifdef _OPENMP
#       pragma omp parallel if (multithread)
#endif
        {
            gaussianBlur(mask, mask, W, H, blur);
        }
    } else if (blur_type == BlurType::BOX) {
        if (int(blur) > 0) {
            for (int i = 0; i < 3; ++i) {
                boxblur(static_cast<float**>(mask), static_cast<float**>(mask), blur, W, H, multithread);
            }
        }
    }
  
}

// basic idea taken from Algorithm 3 in the paper:
// "Parameter-Free Fast Pixelwise Non-Local Means Denoising" http://www.ipol.im/pub/art/2014/120/
// by Jacques Froment

// thanks to Alberto Griggio for this wonderful code
// thanks to Ingo Weyrich <heckflosse67@gmx.de> for many speedup suggestions!
// adapted to Rawtherapee Local adjustments J.Desmis january 2021
//

void ImProcFunctions::NLMeans(float **img, int strength, int detail_thresh, int patch, int radius, float gam, int bfw, int bfh, float scale, bool multithread)
{
    if (!strength) {
        return;
    }
   // printf("Scale=%f\n", scale);
    if(scale > 5.f) {//avoid to small values - leads to crash - but enough to evaluate noise 
        return;
    }
    BENCHFUN
    const int W = bfw;
    const int H = bfh;
//    printf("W=%i H=%i\n", W, H);
    float gamma = gam;
    rtengine::GammaValues g_a; //gamma parameters
    double pwr = 1.0 / static_cast<double>(gam);//default 3.0 - gamma Lab
    double ts = 9.03296;//always the same 'slope' in the extreme shadows - slope Lab
    rtengine::Color::calcGamma(pwr, ts, g_a); // call to calcGamma with selected gamma and slope

    //first change Lab L to pseudo linear with gamma = 3.f slope 9.032...and in range 0...65536, or with gamma slope Lab

#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        int x = 0;
#ifdef __SSE2__
        for (; x < W - 3; x += 4) {
            STVFU(img[y][x], F2V(65536.f) * igammalog(LVFU(img[y][x]) / F2V(32768.f), F2V(gamma), F2V(ts), F2V(g_a[2]), F2V(g_a[4])));
        }
#endif
        for (;x < W; ++x) {
            img[y][x] = 65536.f * igammalog(img[y][x] / 32768.f, gamma, ts, g_a[2], g_a[4]);
        }
    }
    // these two can be changed if needed. increasing max_patch_radius doesn't
    // affect performance, whereas max_search_radius *really* does
    // (the complexity is O(max_search_radius^2 * W * H))
//    constexpr int max_patch_radius = 2;
//    constexpr int max_search_radius = 5;
    int max_patch_radius = patch;
    int max_search_radius = radius;

    const int search_radius = int(std::ceil(float(max_search_radius) / scale));
    const int patch_radius = int(std::ceil(float(max_patch_radius) / scale));

    // the strength parameter controls the scaling of the weights
    // (called h^2 in the papers)
    float eps = 1e-6f;//to avoid too low values and divide near by zero...when  scale > 1
    const float h2 = eps + SQR(std::pow(float(strength) / 100.f, 0.9f) / 30.f / scale);
//    printf("h2=%f\n", h2);
    // this is the main difference between our version and more conventional
    // nl-means implementations: instead of varying the patch size, we control
    // the detail preservation by using a varying weight scaling for the
    // pixels, depending on our estimate of how much details there are in the
    // pixel neighborhood. We do this by computing a "detail mask", using a
    // laplacian filter with additional averaging and smoothing. The
    // detail_thresh parameter controls the degree of detail preservation: the
    // (averaged, smoothed) laplacian is first normalized to [0,1], and then
    // modified by compression and offsetting depending on the detail_thresh
    // parameter, i.e. mask[y][x] = mask[y][x] * (1 - f) + f,
    // where f = detail_thresh / 100
    float amount = LIM(float(detail_thresh)/100.f, 0.f, 0.99f);
    array2D<float> mask(W, H);// ARRAY2D_ALIGNED);
    
    {
        array2D<float> LL(W, H, img, ARRAY2D_BYREFERENCE);
        ImProcFunctions::detail_mask(LL, mask, W, H, 1.f, 1e-3f, 1.f, amount, BlurType::GAUSS, 2.f / scale, multithread);

    }
  
//allocate dst - same type of datas as img
    float** dst;
    int wid = W;
    int hei = H;
    dst = new float*[hei];
    for (int i = 0; i < hei; ++i) {
       dst[i] = new float[wid];
    }
    const int border = search_radius + patch_radius;
    const int WW = W + border * 2;
    const int HH = H + border * 2;

    array2D<float> src(WW, HH);//, ARRAY2D_ALIGNED);
    
#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < HH; ++y) {
        int yy = y <= border ? 0 : y - border >= H ? H-1 : y - border;
        for (int x = 0; x < WW; ++x) {
            int xx = x <= border ? 0 : x - border >= W ? W-1 : x - border;
            float Y = img[yy][xx] / 65536.f;
            src[y][x] = Y;
        }
    }

#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            dst[y][x] = 0.f;
        }
    }

    constexpr int lutsz = 8192;
    constexpr float lutfactor = 100.f / float(lutsz-1);
    LUTf explut(lutsz);
    for (int i = 0; i < lutsz; ++i) {
        float x = float(i) * lutfactor;
        explut[i] = xexpf(-x);
    }

#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            mask[y][x] = (1.f / (mask[y][x] * h2)) / lutfactor;
        }
    }
 
    // process by tiles to avoid numerical accuracy errors in the computation
    // of the integral image
    const int tile_size = 150;
    const int ntiles_x = int(std::ceil(float(WW) / (tile_size-2*border)));
    const int ntiles_y = int(std::ceil(float(HH) / (tile_size-2*border)));
    const int ntiles = ntiles_x * ntiles_y;

#ifdef __SSE2__
    const vfloat zerov = F2V(0.0);
    const vfloat v1e_5f = F2V(1e-5f);
    const vfloat v65536f = F2V(65536.f);
#endif

#ifdef _OPENMP
    #pragma omp parallel if (multithread) 
#endif
    {

#ifdef __SSE2__
    // flush denormals to zero to avoid performance penalty
    const auto oldMode = _MM_GET_FLUSH_ZERO_MODE();
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif
        
#ifdef _OPENMP
    #pragma omp for schedule(dynamic, 2)
#endif
    for (int tile = 0; tile < ntiles; ++tile) {
        const int tile_y = tile / ntiles_x;
        const int tile_x = tile % ntiles_x;

        const int start_y = tile_y * (tile_size - 2*border);
        const int end_y = std::min(start_y + tile_size, HH);
        const int TH = end_y - start_y;

        const int start_x = tile_x * (tile_size - 2*border);
        const int end_x = std::min(start_x + tile_size, WW);
        const int TW = end_x - start_x;

        const auto Yf = [=](int y) -> int { return LIM(y+start_y, 0, HH-1); };
        const auto Xf = [=](int x) -> int { return LIM(x+start_x, 0, WW-1); };

        const auto score =
            [&](int tx, int ty, int zx, int zy) -> float
            {
                return SQR(src[Yf(zy)][Xf(zx)] - src[Yf(zy + ty)][Xf(zx + tx)]);
            };

        array2D<float> St(TW, TH);//, ARRAY2D_ALIGNED);
        array2D<float> SW(TW, TH, ARRAY2D_CLEAR_DATA);//, ARRAY2D_ALIGNED|ARRAY2D_CLEAR_DATA);

        for (int ty = -search_radius; ty <= search_radius; ++ty) {
            for (int tx = -search_radius; tx <= search_radius; ++tx) {
                // Step 1 — Compute the integral image St
                St[0][0] = 0.f;
                for (int xx = 1; xx < TW; ++xx) {
                    St[0][xx] = St[0][xx-1] + score(tx, ty, xx, 0);
                }
                for (int yy = 1; yy < TH; ++yy) {
                    St[yy][0] = St[yy-1][0] + score(tx, ty, 0, yy);
                }
                for (int yy = 1; yy < TH; ++yy) {
                    for (int xx = 1; xx < TW; ++xx) {
                        // operation grouping tuned for performance (empirically)
                        St[yy][xx] = (St[yy][xx-1] + St[yy-1][xx]) - (St[yy-1][xx-1] - score(tx, ty, xx, yy));
                    }
                }
                // Step 2 — Compute weight and estimate for patches
                // V(x), V(y) with y = x + t
                for (int yy = start_y+border; yy < end_y-border; ++yy) {
                    int y = yy - border;
                    int xx = start_x+border;
#ifdef __SSE2__
                    for (; xx < end_x-border-3; xx += 4) {
                        int x = xx - border;
                        int sx = xx + tx;
                        int sy = yy + ty;

                        int sty = yy - start_y;
                        int stx = xx - start_x;
                    
                        vfloat dist2 = LVFU(St[sty + patch_radius][stx + patch_radius]) + LVFU(St[sty - patch_radius][stx - patch_radius]) - LVFU(St[sty + patch_radius][stx - patch_radius]) - LVFU(St[sty - patch_radius][stx + patch_radius]);
                        dist2 = vmaxf(dist2, zerov);
                        vfloat d = dist2 * LVFU(mask[y][x]);
                        vfloat weight = explut[d];
                        STVFU(SW[y-start_y][x-start_x], LVFU(SW[y-start_y][x-start_x]) + weight);
                        vfloat Y = weight * LVFU(src[sy][sx]);
                        STVFU(dst[y][x], LVFU(dst[y][x]) + Y);
                    }
#endif
                    for (; xx < end_x-border; ++xx) {
                        int x = xx - border;
                        int sx = xx + tx;
                        int sy = yy + ty;

                        int sty = yy - start_y;
                        int stx = xx - start_x;
                    
                        float dist2 = St[sty + patch_radius][stx + patch_radius] + St[sty - patch_radius][stx - patch_radius] - St[sty + patch_radius][stx - patch_radius] - St[sty - patch_radius][stx + patch_radius];
                        dist2 = std::max(dist2, 0.f);
                        float d = dist2 * mask[y][x];
                        float weight = explut[d];
                        SW[y-start_y][x-start_x] += weight;
                        float Y = weight * src[sy][sx];
                        dst[y][x] += Y;

                        assert(!xisinff(dst[y][x]));
                        assert(!xisnanf(dst[y][x]));
                    }
                }
            }
        }
//    printf("E\n");
       
        // Compute final estimate at pixel x = (x1, x2)
        for (int yy = start_y+border; yy < end_y-border; ++yy) {
            int y = yy - border;
            int xx = start_x+border;
#ifdef __SSE2__
            for (; xx < end_x-border-3; xx += 4) {
                int x = xx - border;
            
                const vfloat Y = LVFU(dst[y][x]);
                const vfloat f = (v1e_5f + LVFU(SW[y-start_y][x-start_x]));
                STVFU(dst[y][x], (Y / f) * v65536f);
            }
#endif
            for (; xx < end_x-border; ++xx) {
                int x = xx - border;
            
                const float Y = dst[y][x];
                const float f = (1e-5f + SW[y-start_y][x-start_x]);
                dst[y][x] = (Y / f) * 65536.f;

                assert(!xisnanf(dst[y][x]));
            }
        }
    }

#ifdef __SSE2__
    _MM_SET_FLUSH_ZERO_MODE(oldMode);
#endif

    } // omp parallel

#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multithread)
#endif
    for (int y = 0; y < H; ++y) {//apply inverse gamma 3.f and put result in range 32768.f
        int x = 0;
#ifdef __SSE2__
        for (; x < W - 3; x += 4) {
            STVFU(img[y][x], F2V(32768.f) * gammalog(LVFU(dst[y][x]) / F2V(65536.f), F2V(gamma), F2V(ts), F2V(g_a[3]), F2V(g_a[4])));
        }
#endif
        for (; x < W; ++x) {
            img[y][x] = 32768.f * gammalog(dst[y][x] / 65536.f, gamma, ts, g_a[3], g_a[4]);
        }
    }

    for (int i = 0; i < hei; ++i) {
        delete[] dst[i];
    }
    delete[] dst;
    
}

void ImProcFunctions::Lab_Local(
    int call, int sp, float** shbuffer, LabImage * original, LabImage * transformed, LabImage * reserved, LabImage * savenormtm, LabImage * savenormreti, LabImage * lastorig, int fw, int fh, int cx, int cy, int oW, int oH, int sk,
    const LocretigainCurve& locRETgainCcurve, const LocretitransCurve& locRETtransCcurve,
    const LUTf& lllocalcurve, bool locallutili,
    const LUTf& cllocalcurve, bool localclutili,
    const LUTf& lclocalcurve, bool locallcutili,
    const LocLHCurve& loclhCurve,  const LocHHCurve& lochhCurve, const LocCHCurve& locchCurve,
    const LocHHCurve& lochhCurvejz, const LocCHCurve& locchCurvejz, const LocLHCurve& loclhCurvejz,
    const LUTf& lmasklocalcurve, bool localmaskutili,
    const LUTf& lmaskexplocalcurve, bool localmaskexputili,
    const LUTf& lmaskSHlocalcurve, bool localmaskSHutili,
    const LUTf& lmaskviblocalcurve, bool localmaskvibutili,
    const LUTf& lmasktmlocalcurve, bool localmasktmutili,
    LUTf& lmaskretilocalcurve, bool localmaskretiutili,
    const LUTf& lmaskcblocalcurve, bool localmaskcbutili,
    const LUTf& lmaskbllocalcurve, bool localmaskblutili,
    const LUTf& lmasklclocalcurve, bool localmasklcutili,
    const LUTf& lmaskloglocalcurve, bool localmasklogutili,
    const LUTf& lmasklocal_curve, bool localmask_utili,
    const LUTf& lmaskcielocalcurve, bool localmaskcieutili,
    const LUTf& cielocalcurve, bool localcieutili,
    const LUTf& cielocalcurve2, bool localcieutili2,
    const LUTf& jzlocalcurve, bool localjzutili,
    const LUTf& czlocalcurve, bool localczutili,
    const LUTf& czjzlocalcurve, bool localczjzutili,

    const LocCCmaskCurve& locccmasCurve, bool lcmasutili, const LocLLmaskCurve& locllmasCurve, bool llmasutili, const LocHHmaskCurve& lochhmasCurve, bool lhmasutili, const LocHHmaskCurve& llochhhmasCurve, bool lhhmasutili,
    const LocCCmaskCurve& locccmasexpCurve, bool lcmasexputili, const LocLLmaskCurve& locllmasexpCurve, bool llmasexputili, const LocHHmaskCurve& lochhmasexpCurve, bool lhmasexputili,
    const LocCCmaskCurve& locccmasSHCurve, bool lcmasSHutili, const LocLLmaskCurve& locllmasSHCurve, bool llmasSHutili, const LocHHmaskCurve& lochhmasSHCurve, bool lhmasSHutili,
    const LocCCmaskCurve& locccmasvibCurve, bool lcmasvibutili, const LocLLmaskCurve& locllmasvibCurve, bool llmasvibutili, const LocHHmaskCurve& lochhmasvibCurve, bool lhmasvibutili,
    const LocCCmaskCurve& locccmascbCurve, bool lcmascbutili, const LocLLmaskCurve& locllmascbCurve, bool llmascbutili, const LocHHmaskCurve& lochhmascbCurve, bool lhmascbutili,
    const LocCCmaskCurve& locccmasretiCurve, bool lcmasretiutili, const LocLLmaskCurve& locllmasretiCurve, bool llmasretiutili, const LocHHmaskCurve& lochhmasretiCurve, bool lhmasretiutili,
    const LocCCmaskCurve& locccmastmCurve, bool lcmastmutili, const LocLLmaskCurve& locllmastmCurve, bool llmastmutili, const LocHHmaskCurve& lochhmastmCurve, bool lhmastmutili,
    const LocCCmaskCurve& locccmasblCurve, bool lcmasblutili, const LocLLmaskCurve& locllmasblCurve, bool llmasblutili, const LocHHmaskCurve& lochhmasblCurve, bool lhmasblutili,
    const LocCCmaskCurve& locccmaslcCurve, bool lcmaslcutili, const LocLLmaskCurve& locllmaslcCurve, bool llmaslcutili, const LocHHmaskCurve& lochhmaslcCurve, bool lhmaslcutili,
    const LocCCmaskCurve& locccmaslogCurve, bool lcmaslogutili, const LocLLmaskCurve& locllmaslogCurve, bool llmaslogutili, const LocHHmaskCurve& lochhmaslogCurve, bool lhmaslogutili,
    const LocCCmaskCurve& locccmas_Curve, bool lcmas_utili, const LocLLmaskCurve& locllmas_Curve, bool llmas_utili, const LocHHmaskCurve& lochhmas_Curve, bool lhmas_utili,
    const LocCCmaskCurve& locccmascieCurve, bool lcmascieutili, const LocLLmaskCurve& locllmascieCurve, bool llmascieutili, const LocHHmaskCurve& lochhmascieCurve, bool lhmascieutili,

    const LocHHmaskCurve& lochhhmas_Curve, bool lhhmas_utili,
    const LocwavCurve& loclmasCurveblwav, bool lmasutiliblwav,
    const LocwavCurve& loclmasCurvecolwav, bool lmasutilicolwav,
    const LocwavCurve& locwavCurve, bool locwavutili,
    const LocwavCurve& locwavCurvejz, bool locwavutilijz,
    const LocwavCurve& loclevwavCurve, bool loclevwavutili,
    const LocwavCurve& locconwavCurve, bool locconwavutili,
    const LocwavCurve& loccompwavCurve, bool loccompwavutili,
    const LocwavCurve& loccomprewavCurve, bool loccomprewavutili,
    const LocwavCurve& locwavCurvehue, bool locwavhueutili,
    const LocwavCurve& locwavCurveden, bool locwavdenutili,
    const LocwavCurve& locedgwavCurve, bool locedgwavutili,
    const LocwavCurve& loclmasCurve_wav, bool lmasutili_wav,
    
    bool LHutili, bool HHutili, bool CHutili, bool HHutilijz, bool CHutilijz, bool LHutilijz, const LUTf& cclocalcurve, bool localcutili, const LUTf& rgblocalcurve, bool localrgbutili, bool localexutili, const LUTf& exlocalcurve, const LUTf& hltonecurveloc, const LUTf& shtonecurveloc, const LUTf& tonecurveloc, const LUTf& lightCurveloc,
    double& huerefblur, double& chromarefblur, double& lumarefblur, double& hueref, double& chromaref, double& lumaref, double& sobelref, int &lastsav,
    bool prevDeltaE, int llColorMask, int llColorMaskinv, int llExpMask, int llExpMaskinv, int llSHMask, int llSHMaskinv, int llvibMask, int lllcMask, int llsharMask, int llcbMask, int llretiMask, int llsoftMask, int lltmMask, int llblMask, int lllogMask, int ll_Mask, int llcieMask, 
    float& minCD, float& maxCD, float& mini, float& maxi, float& Tmean, float& Tsigma, float& Tmin, float& Tmax,
    float& meantm, float& stdtm, float& meanreti, float& stdreti, float &fab,
    float& highresi, float& nresi, float& highresi46, float& nresi46, float& Lhighresi, float& Lnresi, float& Lhighresi46, float& Lnresi46

    )
{
    //general call of others functions : important return hueref, chromaref, lumaref
    if (!params->locallab.enabled) {
        return;
    }

    //BENCHFUN

    constexpr int del = 3; // to avoid crash with [loy - begy] and [lox - begx] and bfh bfw  // with gtk2 [loy - begy-1] [lox - begx -1 ] and del = 1
    struct local_params lp;
    calcLocalParams(sp, oW, oH, params->locallab, lp, prevDeltaE, llColorMask, llColorMaskinv, llExpMask, llExpMaskinv, llSHMask, llSHMaskinv, llvibMask, lllcMask, llsharMask, llcbMask, llretiMask, llsoftMask, lltmMask, llblMask, lllogMask, ll_Mask, llcieMask, locwavCurveden, locwavdenutili);

    //avoidcolshi(lp, sp, transformed, reserved,  cy, cx, sk);

    const float radius = lp.rad / (sk * 1.4); //0 to 70 ==> see skip
    int levred;
    bool noiscfactiv;

    if (lp.qualmet == 2) { //suppress artifacts with quality enhanced
        levred = 4;
        noiscfactiv = true;
    } else {
        levred = 7;
        noiscfactiv = false;
    }

//lastsav for save restore image
    lastsav = 0;

    if (lp.excmet == 1 && call <= 3 && lp.activspot) {//exclude
        const int bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
        const int bfw = int (lp.lx + lp.lxL) + del;
        const int begy = lp.yc - lp.lyT;
        const int begx = lp.xc - lp.lxL;
        const int yEn = lp.yc + lp.ly;
        const int xEn = lp.xc + lp.lx;
        LabImage bufreserv(bfw, bfh);
        array2D<float> bufsob(bfw, bfh);

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if(multiThread)
#endif
        for (int y = rtengine::max(begy - cy, 0); y < rtengine::min(yEn - cy, original->H); y++) {
            const int loy = cy + y;

            for (int x = rtengine::max(begx - cx, 0); x < rtengine::min(xEn - cx, original->W); x++) {
                const int lox = cx + x;

                bufsob[loy - begy][lox - begx] = bufreserv.L[loy - begy][lox - begx] = reserved->L[y][x];
                bufreserv.a[loy - begy][lox - begx] = reserved->a[y][x];
                bufreserv.b[loy - begy][lox - begx] = reserved->b[y][x];
            }
        }

        array2D<float> ble(bfw, bfh);
        const float radiussob = 1.f / (sk * 1.4f);
        SobelCannyLuma(ble, bufsob, bfw, bfh, radiussob);
        array2D<float> &guid = bufsob;

#ifdef _OPENMP
        #pragma omp parallel for if(multiThread)
#endif
        for (int ir = 0; ir < bfh; ir++)
            for (int jr = 0; jr < bfw; jr++) {
                ble[ir][jr] /= 32768.f;
                guid[ir][jr] /= 32768.f;
            }


        const float blur = 25 / sk * (2.f + 2.5f * lp.struexp);

        rtengine::guidedFilter(guid, ble, ble, blur, 0.0001, multiThread);

//        const float blur = 25 / sk * (10.f + 0.8f * lp.struexp);

//        rtengine::guidedFilter(guid, ble, ble, blur, 0.001, multiThread);

        double sombel = 0.f;
        const int ncsobel = bfh * bfw;

        array2D<float> &deltasobelL = guid;

#ifdef _OPENMP
        #pragma omp parallel for reduction(+:sombel) if(multiThread)
#endif
        for (int ir = 0; ir < bfh; ir++) {
            for (int jr = 0; jr < bfw; jr++) {
                const float val = ble[ir][jr] * 32768.f;
                sombel += static_cast<double>(val);
                deltasobelL[ir][jr] = val;
            }
        }

        const float meansob = sombel / ncsobel;
        Exclude_Local(deltasobelL, hueref, chromaref, lumaref, sobelref, meansob, lp, original, transformed, &bufreserv, reserved, cx, cy, sk);
    }

//encoding lab at the beginning
    if (lp.logena && (call <=3 || lp.prevdE || lp.showmasklogmet == 2 || lp.enaLMask || lp.showmasklogmet == 3 || lp.showmasklogmet == 4)) {
        
        const int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
        const int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
        const int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
        const int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
        const int bfh = yend - ystart;
        const int bfw = xend - xstart;

        if (bfh >= mSP && bfw >= mSP) {
            const std::unique_ptr<LabImage> bufexporig(new LabImage(bfw, bfh)); //buffer for data in zone limit
            const std::unique_ptr<LabImage> bufexpfin(new LabImage(bfw, bfh)); //buffer for data in zone limit

            std::unique_ptr<LabImage> bufmaskblurlog;
            std::unique_ptr<LabImage> originalmasklog;
            std::unique_ptr<LabImage> bufmaskoriglog;

            if (lp.showmasklogmet == 2  || lp.enaLMask || lp.showmasklogmet == 3 || lp.showmasklogmet == 4) {
                bufmaskblurlog.reset(new LabImage(bfw, bfh));
                originalmasklog.reset(new LabImage(bfw, bfh));
                bufmaskoriglog.reset(new LabImage(bfw, bfh));
            }

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if(multiThread)
#endif
            for (int y = ystart; y < yend; y++) {
                for (int x = xstart; x < xend; x++) {
                    bufexporig->L[y - ystart][x - xstart] = original->L[y][x];
                    bufexporig->a[y - ystart][x - xstart] = original->a[y][x];
                    bufexporig->b[y - ystart][x - xstart] = original->b[y][x];
                }
            }

            int inv = 0;
            bool showmaske = false;
            bool enaMask = false;
            bool deltaE = false;
            bool modmask = false;
            bool zero = false;
            bool modif = false;

            if (lp.showmasklogmet == 3) {
                showmaske = true;
            }

            if (lp.enaLMask) {
                enaMask = true;
            }

            if (lp.showmasklogmet == 4) {
                deltaE = true;
            }

            if (lp.showmasklogmet == 2) {
                modmask = true;
            }

            if (lp.showmasklogmet == 1) {
                modif = true;
            }

            if (lp.showmasklogmet == 0) {
                zero = true;
            }
            float chrom = lp.chromaL;
            float rad = lp.radmaL;
            float blendm = lp.blendmaL;
            float gamma = 1.f;
            float slope = 0.f;
            float lap = 0.f; //params->locallab.spots.at(sp).lapmaskexp;
            bool pde = false; //params->locallab.spots.at(sp).laplac;
            LocwavCurve dummy;
            bool delt = params->locallab.spots.at(sp).deltae;
            int sco = params->locallab.spots.at(sp).scopemask;
            int shado = 0;
            int shortcu = 0;//lp.mergemet; //params->locallab.spots.at(sp).shortc;
            const float mindE = 2.f + MINSCOPE * sco * lp.thr;
            const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
            const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
            const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
            float amountcd = 0.f;
            float anchorcd = 50.f;
            int lumask = params->locallab.spots.at(sp).lumask;
            LocHHmaskCurve lochhhmasCurve;
            const int highl = 0;
            maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, bufexporig.get(), bufmaskoriglog.get(), originalmasklog.get(), original, reserved, inv, lp,
                        0.f, false,
                        locccmaslogCurve, lcmaslogutili, locllmaslogCurve, llmaslogutili, lochhmaslogCurve, lhmaslogutili, lochhhmasCurve, false, multiThread,
                        enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, highl, amountcd, anchorcd, lmaskloglocalcurve, localmasklogutili, dummy, false, 1, 1, 5, 5,
                        shortcu, delt, hueref, chromaref, lumaref,
                        maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, -1, fab
                       );

            if (lp.showmasklogmet == 3) {
                showmask(lumask, lp, xstart, ystart, cx, cy, bfw, bfh, bufexporig.get(), transformed, bufmaskoriglog.get(), 0);

                return;
            }
            if (lp.showmasklogmet == 0 || lp.showmasklogmet == 1  || lp.showmasklogmet == 2 || lp.showmasklogmet == 4 || lp.enaLMask) {

                bufexpfin->CopyFrom(bufexporig.get(), multiThread);
                std::unique_ptr<Imagefloat> tmpImage(new Imagefloat(bfw, bfh));
                std::unique_ptr<Imagefloat> tmpImageorig(new Imagefloat(bfw, bfh));
                lab2rgb(*bufexpfin, *tmpImage, params->icm.workingProfile);
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if(multiThread)
#endif
                    for (int y = 0; y < bfh; y++) {
                        for (int x = 0; x < bfw; x++) {
                            tmpImageorig->r(y, x) = tmpImage->r(y, x);
                            tmpImageorig->g(y, x) = tmpImage->g(y, x);
                            tmpImageorig->b(y, x) = tmpImage->b(y, x);
                        }
                    }
                
                log_encode(tmpImage.get(), lp, multiThread, bfw, bfh);
                const float repart = 1.0 - 0.01 * params->locallab.spots.at(sp).repar;
               
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if(multiThread)
#endif
                    for (int y = 0; y < bfh; y++) {
                        for (int x = 0; x < bfw; x++) {
                            tmpImage->r(y, x) = intp(repart, tmpImageorig->r(y, x), tmpImage->r(y, x));
                            tmpImage->g(y, x) = intp(repart, tmpImageorig->g(y, x), tmpImage->g(y, x));
                            tmpImage->b(y, x) = intp(repart, tmpImageorig->b(y, x), tmpImage->b(y, x));
                        }
                    }
                
                rgb2lab(*tmpImage, *bufexpfin, params->icm.workingProfile);
                
                tmpImageorig.reset();
                tmpImage.reset();
                if (params->locallab.spots.at(sp).ciecam) {
                    bool HHcurvejz = false, CHcurvejz = false, LHcurvejz = false;;
                    ImProcFunctions::ciecamloc_02float(lp, sp, bufexpfin.get(), bfw, bfh, 1, sk, cielocalcurve, localcieutili, cielocalcurve2, localcieutili2, jzlocalcurve, localjzutili, czlocalcurve, localczutili, czjzlocalcurve, localczjzutili, locchCurvejz, lochhCurvejz, loclhCurvejz, HHcurvejz, CHcurvejz, LHcurvejz, locwavCurvejz, locwavutilijz);
                }


                if (params->locallab.spots.at(sp).expcie && params->locallab.spots.at(sp).modecie == "log") {
                    bool HHcurvejz = false;
                    bool CHcurvejz = false;
                    bool LHcurvejz = false;
                    if (params->locallab.spots.at(sp).expcie  && params->locallab.spots.at(sp).modecam == "jz") {//some cam16 elementsfor Jz
                        ImProcFunctions::ciecamloc_02float(lp, sp, bufexpfin.get(), bfw, bfh, 10, sk, cielocalcurve, localcieutili, cielocalcurve2, localcieutili2, jzlocalcurve, localjzutili, czlocalcurve, localczutili, czjzlocalcurve, localczjzutili, locchCurvejz, lochhCurvejz, loclhCurvejz, HHcurvejz, CHcurvejz, LHcurvejz, locwavCurvejz, locwavutilijz);
                    }

                    ImProcFunctions::ciecamloc_02float(lp, sp, bufexpfin.get(),bfw, bfh, 0, sk, cielocalcurve, localcieutili, cielocalcurve2, localcieutili2, jzlocalcurve, localjzutili, czlocalcurve, localczutili, czjzlocalcurve, localczjzutili, locchCurvejz, lochhCurvejz, loclhCurvejz, HHcurvejz, CHcurvejz, LHcurvejz, locwavCurvejz, locwavutilijz);

                    float rad = params->locallab.spots.at(sp).detailcie;
                    loccont(bfw, bfh, bufexpfin.get(), rad, 15.f, sk);
                }

                //here begin graduated filter
                //first solution "easy" but we can do other with log_encode...to see the results
                if (lp.strlog != 0.f) {
                    struct grad_params gplog;
                    calclocalGradientParams(lp, gplog, ystart, xstart, bfw, bfh, 11);
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if(multiThread)
#endif
                    for (int ir = 0; ir < bfh; ir++) {
                        for (int jr = 0; jr < bfw; jr++) {
                            bufexpfin->L[ir][jr] *= ImProcFunctions::calcGradientFactor(gplog, jr, ir);
                        }
                    }
                }
            //end graduated

                    float recoth = lp.recothrl;

                    if(lp.recothrl < 1.f) {
                        recoth = -1.f * recoth + 2.f;
                    }

                    if(lp.enaLMask && lp.recothrl != 1.f) {
                        float hig = lp.higthrl;
                        float low = lp.lowthrl;
                       // float recoth = lp.recothrl;
                        float decay = lp.decayl;
                        bool invmask = false;
                        maskrecov(bufexpfin.get(), original, bufmaskoriglog.get(), bfh, bfw, ystart, xstart, hig, low, recoth, decay, invmask, sk, multiThread);
                    }
                    if(lp.recothrl >= 1.f) {
                        transit_shapedetect2(sp, 0.f, 0.f, call, 11, bufexporig.get(), bufexpfin.get(), originalmasklog.get(), hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);
                    } else {
                        transit_shapedetect2(sp, 0.f, 0.f, call, 11, bufexporig.get(), bufexpfin.get(), nullptr, hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);
                    }
            }
            
            if (lp.recur) {
                original->CopyFrom(transformed, multiThread);
                float avge;
                calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
            }
        }
    }

//Prepare mask for Blur and noise and Denoise
    bool denoiz = false;

    if ((lp.noiself > 0.f || lp.noiself0 > 0.f || lp.noiself2 > 0.f || lp.noiselc > 0.f || lp.wavcurvedenoi || lp.nlstr > 0 || lp.noisecf > 0.f || lp.noisecc > 0.f  || lp.bilat > 0.f) && lp.denoiena) {
        denoiz = true;
    }

    bool blurz = false;
    bool delt = params->locallab.spots.at(sp).deltae;

    if (((static_cast<double>(radius) > 1.5 * GAUSS_SKIP)  || lp.stren > 0.1 || lp.blmet == 1 || lp.guidb > 1 || lp.showmaskblmet == 2  || lp.enablMask || lp.showmaskblmet == 3 || lp.showmaskblmet == 4) && lp.blurena) {
        blurz = true;
    }

    const int TW = transformed->W;
    const int TH = transformed->H;
    const std::unique_ptr<LabImage> bufblorig(new LabImage(TW, TH));

    std::unique_ptr<LabImage> originalmaskbl;
    std::unique_ptr<LabImage> bufmaskorigbl;
    std::unique_ptr<LabImage> bufmaskblurbl;
    std::unique_ptr<LabImage> bufprov(new LabImage(TW, TH));

    if (denoiz || blurz || lp.denoiena || lp.blurena) {
        if (lp.showmaskblmet == 2  || lp.enablMask || lp.showmaskblmet == 3 || lp.showmaskblmet == 4) {
            bufmaskorigbl.reset(new LabImage(TW, TH));
            bufmaskblurbl.reset(new LabImage(TW, TH, true));
            originalmaskbl.reset (new LabImage(TW, TH));
        }

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
        for (int y = 0; y < TH; y++) {
            for (int x = 0; x < TW; x++) {
                bufblorig->L[y][x] = original->L[y][x];
            }
        }

        int inv = 0;
        bool showmaske = false;
        bool enaMask = false;
        bool deltaE = false;
        bool modmask = false;
        bool zero = false;
        bool modif = false;

        if (lp.showmaskblmet == 3) {
            showmaske = true;
        }

        if (lp.enablMask) {
            enaMask = true;
        }

        if (lp.showmaskblmet == 4) {
            deltaE = true;
        }

        if (lp.showmaskblmet == 2) {
            modmask = true;
        }

        if (lp.showmaskblmet == 1) {
            modif = true;
        }

        if (lp.showmaskblmet == 0) {
            zero = true;
        }

        float chrom = lp.chromabl;
        float rad = lp.radmabl;
        float gamma = lp.gammabl;
        float slope = lp.slomabl;
        float blendm = lp.blendmabl;
        float lap = params->locallab.spots.at(sp).lapmaskbl;
        bool pde = params->locallab.spots.at(sp).laplac;
        LocwavCurve dummy;
        int lumask = params->locallab.spots.at(sp).lumask;
        int sco = params->locallab.spots.at(sp).scopemask;
        int shortcu = 0;

        const float mindE = 2.f + MINSCOPE * sco * lp.thr;
        const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
        const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
        const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
        const int shado = params->locallab.spots.at(sp).shadmaskblsha;
        const int highl = params->locallab.spots.at(sp).shadmaskbl;
        constexpr float amountcd = 0.f;
        constexpr float anchorcd = 50.f;
        LocHHmaskCurve lochhhmasCurve;
        const float strumask = 0.02 * params->locallab.spots.at(sp).strumaskbl;
        const bool astool = params->locallab.spots.at(sp).toolbl;
        maskcalccol(false, pde, TW, TH, 0, 0, sk, cx, cy, bufblorig.get(), bufmaskblurbl.get(), originalmaskbl.get(), original, reserved, inv, lp,
                    strumask, astool, locccmasblCurve, lcmasblutili, locllmasblCurve, llmasblutili, lochhmasblCurve, lhmasblutili, lochhhmasCurve, false, multiThread,
                    enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, highl, amountcd, anchorcd, lmaskbllocalcurve,
                    localmaskblutili, loclmasCurveblwav, lmasutiliblwav, 1, 1, 5, 5, shortcu, params->locallab.spots.at(sp).deltae, hueref, chromaref, lumaref,
                    maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, 0, fab
                   );

        if (lp.showmaskblmet == 3) {
            showmask(lumask, lp, 0, 0, cx, cy, TW, TH, bufblorig.get(), transformed, bufmaskblurbl.get(), inv);
            return;
        }

    }

    bool execmaskblur = (lp.showmaskblmet == 2 || lp.enablMask || lp.showmaskblmet == 3 || lp.showmaskblmet == 4) && lp.smasktyp != 1;
    int strengr = params->locallab.spots.at(sp).strengr;

    if (((static_cast<double>(radius) > 1.5 * GAUSS_SKIP && lp.rad > 1.6) || lp.stren > 0.1 || lp.blmet == 1 || lp.guidb > 0 || strengr > 0 || execmaskblur) && lp.blurena) { // radius < GAUSS_SKIP means no gauss, just copy of original image
        std::unique_ptr<LabImage> tmp1;
        std::unique_ptr<LabImage> tmp2;
        std::unique_ptr<LabImage> tmp3;
        std::unique_ptr<LabImage> maskk;
        int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
        int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
        int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
        int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
        int bfh = yend - ystart;
        int bfw = xend - xstart;
        int bfhr = bfh;
        int bfwr = bfw;

        bool fft = params->locallab.spots.at(sp).fftwbl;
        int isogr = params->locallab.spots.at(sp).isogr;
        int scalegr = params->locallab.spots.at(sp).scalegr;
        float  divgr = params->locallab.spots.at(sp).divgr;


        if (bfw >= mSP && bfh >= mSP) {
            if (lp.blurmet == 0 && (fft || lp.rad > 30.0)) {
                optfft(N_fftwsize, bfh, bfw, bfhr, bfwr, lp, original->H, original->W, xstart, ystart, xend, yend, cx, cy, lp.fullim);
            }

            const std::unique_ptr<LabImage> bufgbi(new LabImage(TW, TH));

            //here mask is used with plain image for normal and inverse
            //if it is possible to optimize with maskcalccol(), I don't to preserve visibility
            if (lp.showmaskblmet == 0 || lp.showmaskblmet == 1  || lp.showmaskblmet == 2 || lp.showmaskblmet == 4 || lp.enablMask) {

                if (lp.blurmet == 0) {
                    if (bfw > 0 && bfh > 0) {
                        tmp1.reset(new LabImage(bfw, bfh));
                        tmp3.reset(new LabImage(bfw, bfh));
                        maskk.reset(new LabImage(bfw, bfh));
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = ystart; y < yend ; y++) {
                            for (int x = xstart; x < xend; x++) {
                                tmp1->L[y - ystart][x - xstart] = original->L[y][x];
                                tmp1->a[y - ystart][x - xstart] = original->a[y][x];
                                tmp1->b[y - ystart][x - xstart] = original->b[y][x];
                            }
                        }
                    }
                } else if (lp.blurmet == 1) {
                    tmp1.reset(new LabImage(transformed->W, transformed->H));
                    tmp2.reset(new LabImage(transformed->W, transformed->H));
                    tmp3.reset(new LabImage(transformed->W, transformed->H));

                    for (int y = 0; y < TH ; y++) {
                        for (int x = 0; x < TW; x++) {
                            tmp2->L[y][x] = original->L[y][x];
                            tmp2->a[y][x] = original->a[y][x];
                            tmp2->b[y][x] = original->b[y][x];
                            tmp3->L[y][x] = original->L[y][x];
                            tmp3->a[y][x] = original->a[y][x];
                            tmp3->b[y][x] = original->b[y][x];
                            tmp1->L[y][x] = original->L[y][x];
                            tmp1->a[y][x] = original->a[y][x];
                            tmp1->b[y][x] = original->b[y][x];
                            bufgbi->L[y][x] = original->L[y][x];
                            bufgbi->a[y][x] = original->a[y][x];
                            bufgbi->b[y][x] = original->b[y][x];
                        }
                    }

                }


                if (lp.blurmet == 0 && lp.blmet == 0 && static_cast<double>(radius) > (1.5 * GAUSS_SKIP) && lp.rad > 1.6) {
                    if (fft || lp.rad > 30.0) {
                        if (lp.chromet == 0) {
                            ImProcFunctions::fftw_convol_blur2(tmp1->L, tmp1->L, bfwr, bfhr, radius, 0, 0);
                        } else if (lp.chromet == 1) {
                            ImProcFunctions::fftw_convol_blur2(tmp1->a, tmp1->a, bfwr, bfhr, radius, 0, 0);
                            ImProcFunctions::fftw_convol_blur2(tmp1->b, tmp1->b, bfwr, bfhr, radius, 0, 0);
                        } else if (lp.chromet == 2) {
                            ImProcFunctions::fftw_convol_blur2(tmp1->L, tmp1->L, bfwr, bfhr, radius, 0, 0);
                            ImProcFunctions::fftw_convol_blur2(tmp1->a, tmp1->a, bfwr, bfhr, radius, 0, 0);
                            ImProcFunctions::fftw_convol_blur2(tmp1->b, tmp1->b, bfwr, bfhr, radius, 0, 0);
                        }
                    } else {
#ifdef _OPENMP
                        #pragma omp parallel if (multiThread)
#endif
                        {
                            if (lp.chromet == 0) {
                                gaussianBlur(tmp1->L, tmp1->L, bfw, bfh, radius);
                            } else if (lp.chromet == 1) {
                                gaussianBlur(tmp1->a, tmp1->a, bfw, bfh, radius);
                                gaussianBlur(tmp1->b, tmp1->b, bfw, bfh, radius);
                            } else if (lp.chromet == 2) {
                                gaussianBlur(tmp1->L, tmp1->L, bfw, bfh, radius);
                                gaussianBlur(tmp1->a, tmp1->a, bfw, bfh, radius);
                                gaussianBlur(tmp1->b, tmp1->b, bfw, bfh, radius);
                            }
                        }
                    }

                } else if (lp.blurmet == 1 && lp.blmet == 0 && static_cast<double>(radius) > (1.5 * GAUSS_SKIP) && lp.rad > 1.6) {
                    if (fft || lp.rad > 30.0) {
                        if (lp.chromet == 0) {
                            ImProcFunctions::fftw_convol_blur2(original->L, tmp1->L, TW, TH, radius, 0, 0);
                        } else if (lp.chromet == 1) {
                            ImProcFunctions::fftw_convol_blur2(original->a, tmp1->a, TW, TH, radius, 0, 0);
                            ImProcFunctions::fftw_convol_blur2(original->b, tmp1->b, TW, TH, radius, 0, 0);
                        } else if (lp.chromet == 2) {
                            ImProcFunctions::fftw_convol_blur2(original->L, tmp1->L, TW, TH, radius, 0, 0);
                            ImProcFunctions::fftw_convol_blur2(original->a, tmp1->a, TW, TH, radius, 0, 0);
                            ImProcFunctions::fftw_convol_blur2(original->b, tmp1->b, TW, TH, radius, 0, 0);
                        }
                    } else {
#ifdef _OPENMP
                        #pragma omp parallel if (multiThread)
#endif
                        {
                            if (lp.chromet == 0) {
                                gaussianBlur(original->L, tmp1->L, TW, TH, radius);
                            } else if (lp.chromet == 1) {
                                gaussianBlur(original->a, tmp1->a, TW, TH, radius);
                                gaussianBlur(original->b, tmp1->b, TW, TH, radius);
                            } else if (lp.chromet == 2) {
                                gaussianBlur(original->L, tmp1->L, TW, TH, radius);
                                gaussianBlur(original->a, tmp1->a, TW, TH, radius);
                                gaussianBlur(original->b, tmp1->b, TW, TH, radius);
                            }
                        }
                    }
                }

                //add noise
                if (tmp1.get() && lp.stren > 0.1 && lp.blmet == 0) {
                    float mean = 0.f;//0 best result
                    float variance = lp.stren ;
                    addGaNoise(tmp1.get(), tmp1.get(), mean, variance, sk) ;
                }

                //add grain
                if (lp.blmet == 0 && strengr > 0) {
                    int wi = bfw;
                    int he = bfh;

                    if (lp.blurmet == 1) {
                        wi = TW;
                        he = TH;
                    }

                    if (tmp1.get()) {
                        Imagefloat *tmpImage = nullptr;
                        tmpImage = new Imagefloat(wi, he);

                        for (int y = 0; y < he ; y++) {
                            for (int x = 0; x < wi; x++) {
                                tmpImage->g(y, x) = tmp1->L[y][x];
                                tmpImage->r(y, x) = tmp1->a[y][x];
                                tmpImage->b(y, x) = tmp1->b[y][x];
                            }
                        }


                        filmGrain(tmpImage, isogr, strengr, scalegr, divgr, wi, he, call, fw, fh);

                        for (int y = 0; y < he ; y++) {
                            for (int x = 0; x < wi; x++) {
                                tmp1->L[y][x] = tmpImage->g(y, x);
                                tmp1->a[y][x] = tmpImage->r(y, x);
                                tmp1->b[y][x] = tmpImage->b(y, x);
                            }
                        }

                        delete tmpImage;
                    }
                }

                Median medianTypeL = Median::TYPE_3X3_STRONG;
                Median medianTypeAB = Median::TYPE_3X3_STRONG;

                if (lp.medmet == 0) {
                    medianTypeL = medianTypeAB = Median::TYPE_3X3_STRONG;
                } else if (lp.medmet == 1) {
                    medianTypeL = medianTypeAB = Median::TYPE_5X5_STRONG;
                } else if (lp.medmet == 2) {
                    medianTypeL = medianTypeAB = Median::TYPE_7X7;
                } else if (lp.medmet == 3) {
                    medianTypeL = medianTypeAB = Median::TYPE_9X9;
                }

                if (lp.blurmet == 0 && lp.blmet == 1 && lp.medmet != -1) {
                    float** tmL;
                    int wid = bfw;
                    int hei = bfh;
                    tmL = new float*[hei];

                    for (int i = 0; i < hei; ++i) {
                        tmL[i] = new float[wid];
                    }

                    if (lp.chromet == 0) {
                        Median_Denoise(tmp1->L, tmp1->L, bfw, bfh, medianTypeL, lp.it, multiThread, tmL);
                    }

                    else if (lp.chromet == 1) {
                        Median_Denoise(tmp1->a, tmp1->a, bfw, bfh, medianTypeAB, lp.it, multiThread, tmL);
                        Median_Denoise(tmp1->b, tmp1->b, bfw, bfh, medianTypeAB, lp.it, multiThread, tmL);
                    } else if (lp.chromet == 2) {
                        Median_Denoise(tmp1->L, tmp1->L, bfw, bfh, medianTypeL, lp.it, multiThread, tmL);
                        Median_Denoise(tmp1->a, tmp1->a, bfw, bfh, medianTypeAB, lp.it, multiThread, tmL);
                        Median_Denoise(tmp1->b, tmp1->b, bfw, bfh, medianTypeAB, lp.it, multiThread, tmL);
                    }

                    for (int i = 0; i < hei; ++i) {
                        delete[] tmL[i];
                    }

                    delete[] tmL;

                } else if (lp.blurmet == 1 && lp.blmet == 1) {
                    float** tmL;
                    int wid = TW;
                    int hei = TH;
                    tmL = new float*[hei];

                    for (int i = 0; i < hei; ++i) {
                        tmL[i] = new float[wid];
                    }

                    if (lp.chromet == 0) {
                        Median_Denoise(tmp2->L, tmp1->L, TW, TH, medianTypeL, lp.it, multiThread, tmL);
                    } else if (lp.chromet == 1) {
                        Median_Denoise(tmp2->a, tmp1->a, TW, TH, medianTypeAB, lp.it, multiThread, tmL);
                        Median_Denoise(tmp2->b, tmp1->b, TW, TH, medianTypeAB, lp.it, multiThread, tmL);
                    } else if (lp.chromet == 2) {
                        Median_Denoise(tmp2->L, tmp1->L, TW, TH, medianTypeL, lp.it, multiThread, tmL);
                        Median_Denoise(tmp2->a, tmp1->a, TW, TH, medianTypeAB, lp.it, multiThread, tmL);
                        Median_Denoise(tmp2->b, tmp1->b, TW, TH, medianTypeAB, lp.it, multiThread, tmL);
                    }

                    for (int i = 0; i < hei; ++i) {
                        delete[] tmL[i];
                    }

                    delete[] tmL;
                }

                if (lp.blurmet == 0 && lp.blmet == 2) {

                    if (lp.guidb > 0) {
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = ystart; y < yend ; y++) {
                            for (int x = xstart; x < xend; x++) {
                                tmp1->L[y - ystart][x - xstart] = original->L[y][x];
                                tmp1->a[y - ystart][x - xstart] = original->a[y][x];
                                tmp1->b[y - ystart][x - xstart] = original->b[y][x];
                                tmp3->L[y - ystart][x - xstart] = original->L[y][x];
                                tmp3->a[y - ystart][x - xstart] = original->a[y][x];
                                tmp3->b[y - ystart][x - xstart] = original->b[y][x];
                            }
                        }

                        Imagefloat *tmpImage = nullptr;
                        tmpImage = new Imagefloat(bfw, bfh);
                        lab2rgb(*tmp1, *tmpImage, params->icm.workingProfile);
                        array2D<float> LL(bfw, bfh);
                        array2D<float> rr(bfw, bfh);
                        array2D<float> gg(bfw, bfh);
                        array2D<float> bb(bfw, bfh);
                        array2D<float> guide(bfw, bfh);
                        
                        
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < bfh ; y++) {
                            for (int x = 0; x < bfw; x++) {
                                LL[y][x] = tmp1->L[y][x];
                                float ll = LL[y][x] / 32768.f;
                                guide[y][x] = xlin2log(rtengine::max(ll, 0.f), 10.f);
                                rr[y][x] = tmpImage->r(y, x);
                                gg[y][x] = tmpImage->g(y, x);
                                bb[y][x] = tmpImage->b(y, x);

                            }
                        }
                        array2D<float> iR(bfw, bfh, rr, 0);
                        array2D<float> iG(bfw, bfh, gg, 0);
                        array2D<float> iB(bfw, bfh, bb, 0);
                        array2D<float> iL(bfw, bfh, LL, 0);

                        int r = rtengine::max(int(lp.guidb / sk), 1);

                        const float epsil = 0.001f * std::pow(2.f, -lp.epsb);

                        if (lp.chromet == 0) {
                            rtengine::guidedFilterLog(guide, 10.f, LL, r, epsil, multiThread);
                        } else if (lp.chromet == 1) {
                            rtengine::guidedFilterLog(guide, 10.f, rr, r, epsil, multiThread);
                            rtengine::guidedFilterLog(guide, 10.f, bb, r, epsil, multiThread);
                        } else if (lp.chromet == 2) {
                            rtengine::guidedFilterLog(10.f, gg, r, epsil, multiThread);
                            rtengine::guidedFilterLog(10.f, rr, r, epsil, multiThread);
                            rtengine::guidedFilterLog(10.f, bb, r, epsil, multiThread);
                        }

#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < bfh ; y++) {
                            for (int x = 0; x < bfw; x++) {
                                rr[y][x] = intp(lp.strbl, rr[y][x] , iR[y][x]);
                                gg[y][x] = intp(lp.strbl, gg[y][x] , iG[y][x]);
                                bb[y][x] = intp(lp.strbl, bb[y][x] , iB[y][x]);
                                tmpImage->r(y, x) = rr[y][x];
                                tmpImage->g(y, x) = gg[y][x];
                                tmpImage->b(y, x) = bb[y][x];

                            }
                        }

                        rgb2lab(*tmpImage, *tmp1, params->icm.workingProfile);

                        if (lp.chromet == 0) {
#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int y = 0; y < bfh ; y++) {
                                for (int x = 0; x < bfw; x++) {
                                    LL[y][x] = intp(lp.strbl, LL[y][x] , iL[y][x]);
                                    tmp1->L[y][x] = LL[y][x];
                                }
                            }
                        }
                        if(lp.enablMask && lp.recothr != 1.f && lp.smasktyp != 1) {
                            array2D<float> masklum;
                            masklum(bfw, bfh);
                            for (int ir = 0; ir < bfh; ir++)
                                for (int jr = 0; jr < bfw; jr++) {
                                    masklum[ir][jr] = 1.f;
                                }

                            float hig = lp.higthr;
                            float higc;
                            calcdif(hig, higc);
                            float low = lp.lowthr;
                            float lowc;
                            calcdif(low, lowc);

                            if(higc < lowc) {
                                higc = lowc + 0.01f;
                            }
                            float th = (lp.recothr - 1.f);
                            float ahigh = th / (higc - 100.f);
                            float bhigh = 1.f - higc * ahigh;

                            float alow = th / lowc; 
                            float blow = 1.f - th;
                            
#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                            for (int ir = 0; ir < bfh; ir++)
                                for (int jr = 0; jr < bfw; jr++) {
                                    const float lM = bufmaskblurbl->L[ir + ystart][jr + xstart];
                                    const float lmr = lM / 327.68f;
                                    if (lM < 327.68f * lowc) {
                                        masklum[ir][jr] = alow * lmr + blow;
                                    } else if (lM < 327.68f * higc) {
                                    
                                    } else {
                                        masklum[ir][jr] = ahigh * lmr + bhigh;
                                    }
                                    if(lp.invmask == true) {
                                        float k = masklum[ir][jr];
                                        masklum[ir][jr] = 1 - k*k;
                                    }
                                }

                            for (int i = 0; i < 3; ++i) {
                                boxblur(static_cast<float**>(masklum), static_cast<float**>(masklum), 10 / sk, bfw, bfh, false);
                            }

#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                            for (int i = 0; i < bfh; ++i) {
                                for (int j = 0; j < bfw; ++j) {
                                    tmp1->L[i][j] = (tmp3->L[i][j] - tmp1->L[i][j]) *  LIM01(masklum[i][j]) + tmp1->L[i][j];
                                    tmp1->a[i][j] = (tmp3->a[i][j] - tmp1->a[i][j]) *  LIM01(masklum[i][j]) + tmp1->a[i][j];
                                    tmp1->b[i][j] = (tmp3->b[i][j] - tmp1->b[i][j]) *  LIM01(masklum[i][j]) + tmp1->b[i][j];
                                }
                            }
                            masklum.free();
                        }

                        delete tmpImage;
                    }

                } else if (lp.blurmet == 1 && lp.blmet == 2) {

                    if (lp.guidb > 0) {
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < TH ; y++) {
                            for (int x = 0; x < TW; x++) {
                                tmp1->L[y][x] = original->L[y][x];
                                tmp1->a[y][x] = original->a[y][x];
                                tmp1->b[y][x] = original->b[y][x];
                                tmp2->L[y][x] = original->L[y][x];
                                tmp3->L[y][x] = original->L[y][x];
                                tmp3->a[y][x] = original->a[y][x];
                                tmp3->b[y][x] = original->b[y][x];
                            }
                        }

                        Imagefloat *tmpImage = nullptr;
                        tmpImage = new Imagefloat(TW, TH);
                        lab2rgb(*tmp1, *tmpImage, params->icm.workingProfile);
                        array2D<float> LL(TW, TH);
                        array2D<float> rr(TW, TH);
                        array2D<float> gg(TW, TH);
                        array2D<float> bb(TW, TH);
                        array2D<float> guide(TW, TH);

#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < TH ; y++) {
                            for (int x = 0; x < TW; x++) {
                                LL[y][x] = tmp1->L[y][x];
                                float ll = LL[y][x] / 32768.f;
                                guide[y][x] = xlin2log(rtengine::max(ll, 0.f), 10.f);
                                rr[y][x] = tmpImage->r(y, x);
                                gg[y][x] = tmpImage->g(y, x);
                                bb[y][x] = tmpImage->b(y, x);

                            }
                        }

                        array2D<float> iR(TW, TH, rr, 0);
                        array2D<float> iG(TW, TH, gg, 0);
                        array2D<float> iB(TW, TH, bb, 0);
                        array2D<float> iL(TW, TH, LL, 0);

                        int r = rtengine::max(int(lp.guidb / sk), 1);

                        const float epsil = 0.001f * std::pow(2.f, - lp.epsb);

                        if (lp.chromet == 0) {
                            rtengine::guidedFilterLog(guide, 10.f, LL, r, epsil, multiThread);
                        } else if (lp.chromet == 1) {
                            rtengine::guidedFilterLog(guide, 10.f, rr, r, epsil, multiThread);
                            rtengine::guidedFilterLog(guide, 10.f, bb, r, epsil, multiThread);
                        } else if (lp.chromet == 2) {
                            rtengine::guidedFilterLog(10.f, gg, r, epsil, multiThread);
                            rtengine::guidedFilterLog(10.f, rr, r, epsil, multiThread);
                            rtengine::guidedFilterLog(10.f, bb, r, epsil, multiThread);
                        }

#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < TH ; y++) {
                            for (int x = 0; x < TW; x++) {
                                rr[y][x] = intp(lp.strbl, rr[y][x] , iR[y][x]);
                                gg[y][x] = intp(lp.strbl, gg[y][x] , iG[y][x]);
                                bb[y][x] = intp(lp.strbl, bb[y][x] , iB[y][x]); 
                                tmpImage->r(y, x) = rr[y][x];
                                tmpImage->g(y, x) = gg[y][x];
                                tmpImage->b(y, x) = bb[y][x];

                            }
                        }

                        rgb2lab(*tmpImage, *tmp1, params->icm.workingProfile);

                        if (lp.chromet == 0) {
#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int y = 0; y < TH ; y++) {
                                for (int x = 0; x < TW; x++) {
                                    LL[y][x] = intp(lp.strbl, LL[y][x] , iL[y][x]);
                                    tmp1->L[y][x] = LL[y][x];
                                }
                            }
                        }
                        if(lp.enablMask && lp.recothr != 1.f && lp.smasktyp != 1) {
                            array2D<float> masklum;
                            masklum(TW, TH);
                            for (int ir = 0; ir < TH; ir++)
                                for (int jr = 0; jr < TW; jr++) {
                                    masklum[ir][jr] = 1.f;
                                }

                            float hig = lp.higthr;
                            float higc;
                            calcdif(hig, higc);
                            float low = lp.lowthr;
                            float lowc;
                            calcdif(low, lowc);

                            if(higc < lowc) {
                                higc = lowc + 0.01f;
                            }
                            float th = (lp.recothr - 1.f);
                            float ahigh = th / (higc - 100.f);
                            float bhigh = 1.f - higc * ahigh;

                            float alow = th / lowc; 
                            float blow = 1.f - th;
                            
#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                            for (int ir = 0; ir < TH; ir++)
                                for (int jr = 0; jr < TW; jr++) {
                                    const float lM = bufmaskblurbl->L[ir][jr];
                                    const float lmr = lM / 327.68f;
                                    if (lM < 327.68f * lowc) {
                                        masklum[ir][jr] = alow * lmr + blow;
                                    } else if (lM < 327.68f * higc) {
                                    
                                    } else {
                                        masklum[ir][jr] = (ahigh * lmr + bhigh);
                                    }
                                }

                            for (int i = 0; i < 3; ++i) {
                                boxblur(static_cast<float**>(masklum), static_cast<float**>(masklum), 10 / sk, TW, TH, false);
                            }

#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                            for (int i = 0; i < TH; ++i) {
                                for (int j = 0; j < TW; ++j) {
                                    tmp1->L[i][j] = (tmp3->L[i][j] - tmp1->L[i][j]) *  LIM01(masklum[i][j]) + tmp1->L[i][j];
                                    tmp1->a[i][j] = (tmp3->a[i][j] - tmp1->a[i][j]) *  LIM01(masklum[i][j]) + tmp1->a[i][j];
                                    tmp1->b[i][j] = (tmp3->b[i][j] - tmp1->b[i][j]) *  LIM01(masklum[i][j]) + tmp1->b[i][j];
                                }
                            }
                            masklum.free();

                        }

                        delete tmpImage;
                    }
                }

                if (tmp1.get()) {
                    if (lp.blurmet == 0) { //blur and noise (center)
                        
                        if(lp.smasktyp != 1) {
                            BlurNoise_Local(tmp1.get(), originalmaskbl.get(),  hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);
                        } else {
                            BlurNoise_Local(tmp1.get(), original, hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);
                        }

                        if (lp.recur) {
                            original->CopyFrom(transformed, multiThread);
                            float avge;
                            calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                        }
                    } else if (lp.blurmet == 1) {
 //                       InverseBlurNoise_Local(originalmaskbl, bufchro, lp, hueref, chromaref, lumaref, original, transformed, tmp1.get(), cx, cy, sk);
                        if(lp.smasktyp != 1) {
                            InverseBlurNoise_Local(originalmaskbl.get(), lp, hueref, chromaref, lumaref, original, transformed, tmp1.get(), cx, cy, sk);
                        } else {
                            InverseBlurNoise_Local(original, lp, hueref, chromaref, lumaref, original, transformed, tmp1.get(), cx, cy, sk);
                        }

                        if (lp.recur) {
                            original->CopyFrom(transformed, multiThread);
                            float avge;
                            calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                        }
                    }
                }
            }
        }
    }

    //local impulse
    if ((lp.bilat > 0.f) && lp.denoiena) {
        const int bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
        const int bfw = int (lp.lx + lp.lxL) + del;

        std::unique_ptr<LabImage> bufwv;

        if (call == 2) {//simpleprocess
            bufwv.reset(new LabImage(bfw, bfh)); //buffer for data in zone limit

            const int begy = lp.yc - lp.lyT;
            const int begx = lp.xc - lp.lxL;
            const int yEn = lp.yc + lp.ly;
            const int xEn = lp.xc + lp.lx;

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = rtengine::max(0, begy - cy); y < rtengine::min(transformed->H, yEn - cy); y++) {
                const int loy = cy + y;

                for (int x = rtengine::max(0, begx - cx); x < rtengine::min(transformed->W, xEn - cx); x++) {
                    const int lox = cx + x;
                    bufwv->L[loy - begy][lox - begx] = original->L[y][x];
                    bufwv->a[loy - begy][lox - begx] = original->a[y][x];
                    bufwv->b[loy - begy][lox - begx] = original->b[y][x];
                }
            }
        } else {//dcrop.cc
            bufwv.reset(new LabImage(transformed->W, transformed->H));
            bufwv->CopyFrom(original, multiThread);
        } //end dcrop

        const double threshold = lp.bilat / 20.f;

        if (bfh > 8 && bfw > 8) {
            ImProcFunctions::impulse_nr(bufwv.get(), threshold);
        }

        DeNoise_Local(call, lp,  originalmaskbl.get(), levred, huerefblur, lumarefblur, chromarefblur, original, transformed, *bufwv, cx, cy, sk);

        if (lp.recur) {
            original->CopyFrom(transformed, multiThread);
            float avge;
            calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
        }
    }

//local denoise
    if (lp.activspot && lp.denoiena && (lp.noiself > 0.f || lp.noiself0 > 0.f || lp.noiself2 > 0.f || lp.wavcurvedenoi ||lp.nlstr > 0 || lp.noiselc > 0.f || lp.noisecf > 0.f || lp.noisecc > 0.f )) {//disable denoise if not used
        constexpr int aut = 0;
        DeNoise(call, aut, noiscfactiv, lp, originalmaskbl.get(), bufmaskblurbl.get(), levred, huerefblur, lumarefblur, chromarefblur, original, transformed, cx, cy, sk, locwavCurvehue, locwavhueutili,
                highresi, nresi, highresi46, nresi46, Lhighresi, Lnresi, Lhighresi46, Lnresi46);
        if (lp.recur) {
            original->CopyFrom(transformed, multiThread);
            float avge;
            calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
        }
    }

//Tone mapping

    if ((lp.strengt != 0.f || lp.showmasktmmet == 2 || lp.enatmMask || lp.showmasktmmet == 3 || lp.showmasktmmet == 4 || lp.prevdE) && lp.tonemapena && !params->epd.enabled) {
        if (call <= 3) { //simpleprocess dcrop improcc
            const int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
            const int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
            const int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
            const int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
            const int bfh = yend - ystart;
            const int bfw = xend - xstart;

            if (bfw >= mDEN && bfh >= mDEN) {
               // printf("OK TM\n");
                array2D<float> buflight(bfw, bfh);
                JaggedArray<float> bufchro(bfw, bfh);
                std::unique_ptr<LabImage> bufgb(new LabImage(bfw, bfh));
                const std::unique_ptr<LabImage> tmp1(new LabImage(bfw, bfh));
                const std::unique_ptr<LabImage> bufgbm(new LabImage(bfw, bfh));
                const std::unique_ptr<LabImage> tmp1m(new LabImage(bfw, bfh));
                std::unique_ptr<LabImage> bufmaskorigtm;
                std::unique_ptr<LabImage> bufmaskblurtm;
                std::unique_ptr<LabImage> originalmasktm;

                //       if (lp.showmasktmmet == 0 || lp.showmasktmmet == 2  || lp.enatmMask || lp.showmasktmmet == 3 || lp.showmasktmmet == 4) {
                if (lp.showmasktmmet == 2  || lp.enatmMask || lp.showmasktmmet == 3 || lp.showmasktmmet == 4) {
                    bufmaskorigtm.reset(new LabImage(bfw, bfh));
                    bufmaskblurtm.reset(new LabImage(bfw, bfh));
                    originalmasktm.reset(new LabImage(bfw, bfh));
                }

                // 3 loops to avoid performance penalty on machines with 4-way L1 cache
#ifdef _OPENMP
                #pragma omp parallel if (multiThread)
                {
                #pragma omp for schedule(dynamic,16) nowait
#endif
                for (int y = ystart; y < yend; y++) {
                    for (int x = xstart; x < xend; x++) {
                        bufgbm->L[y - ystart][x - xstart] = bufgb->L[y - ystart][x - xstart] = original->L[y][x];
                    }
                }

#ifdef _OPENMP
                #pragma omp for schedule(dynamic,16) nowait
#endif
                for (int y = ystart; y < yend; y++) {
                    for (int x = xstart; x < xend; x++) {
                        bufgbm->a[y - ystart][x - xstart] = bufgb->a[y - ystart][x - xstart] = original->a[y][x];
                    }
                }

#ifdef _OPENMP
                #pragma omp for schedule(dynamic,16)
#endif
                for (int y = ystart; y < yend; y++) {
                    for (int x = xstart; x < xend; x++) {
                        bufgbm->b[y - ystart][x - xstart] = bufgb->b[y - ystart][x - xstart] = original->b[y][x];
                    }
                }
#ifdef _OPENMP
                }
#endif

                int inv = 0;
                bool showmaske = false;
                bool enaMask = false;
                bool deltaE = false;
                bool modmask = false;
                bool zero = false;
                bool modif = false;

                if (lp.showmasktmmet == 3) {
                    showmaske = true;
                }

                if (lp.enatmMask) {
                    enaMask = true;
                }

                if (lp.showmasktmmet == 4) {
                    deltaE = true;
                }

                if (lp.showmasktmmet == 2) {
                    modmask = true;
                }

                if (lp.showmasktmmet == 1) {
                    modif = true;
                }

                if (lp.showmasktmmet == 0) {
                    zero = true;
                }

                float chrom = lp.chromatm;;
                float rad = lp.radmatm;
                float gamma = lp.gammatm;
                float slope = lp.slomatm;
                float blendm = lp.blendmatm;
                float lap = params->locallab.spots.at(sp).lapmasktm;
                bool pde = params->locallab.spots.at(sp).laplac;
                int lumask = params->locallab.spots.at(sp).lumask;

                if (!params->locallab.spots.at(sp).enatmMaskaft) {
                    LocwavCurve dummy;
                    int sco = params->locallab.spots.at(sp).scopemask;
                    int shortcu = 0; //lp.mergemet;// params->locallab.spots.at(sp).shortc;

                    const float mindE = 2.f + MINSCOPE * sco * lp.thr;
                    const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
                    const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
                    const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
                    int shado = 0;
                    float amountcd = 0.f;
                    float anchorcd = 50.f;
                    LocHHmaskCurve lochhhmasCurve;
                    const int highl = 0;
                    maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, bufgbm.get(), bufmaskorigtm.get(), originalmasktm.get(), original, reserved, inv, lp,
                                0.f, false,
                                locccmastmCurve, lcmastmutili, locllmastmCurve, llmastmutili, lochhmastmCurve, lhmastmutili, lochhhmasCurve, false, multiThread,
                                enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, highl, amountcd, anchorcd, lmasktmlocalcurve, localmasktmutili, dummy, false, 1, 1, 5, 5,
                                shortcu, params->locallab.spots.at(sp).deltae, hueref, chromaref, lumaref,
                                maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, -1, fab
                               );

                    if (lp.showmasktmmet == 3) {
                        showmask(lumask, lp, xstart, ystart, cx, cy, bfw, bfh, bufgbm.get(), transformed, bufmaskorigtm.get(), 0);

                        return;
                    }
                }

                if (lp.showmasktmmet == 0 || lp.showmasktmmet == 1  || lp.showmasktmmet == 2 || lp.showmasktmmet == 4 || lp.showmasktmmet == 3 || lp.enatmMask) {
                    constexpr int itera = 0;
                    ImProcFunctions::EPDToneMaplocal(sp, bufgb.get(), tmp1.get(), itera, sk);//iterate to 0 calculate with edgstopping, improve result, call=1 dcrop we can put iterate to 5

                if (params->locallab.spots.at(sp).expcie && params->locallab.spots.at(sp).modecie == "tm") {
                    bool HHcurvejz = false;
                    bool CHcurvejz = false;
                    bool LHcurvejz = false;
                    if (params->locallab.spots.at(sp).modecam == "jz") {//some cam16 elementsfor Jz
                        ImProcFunctions::ciecamloc_02float(lp, sp, tmp1.get(), bfw, bfh, 10, sk, cielocalcurve, localcieutili, cielocalcurve2, localcieutili2, jzlocalcurve, localjzutili, czlocalcurve, localczutili, czjzlocalcurve, localczjzutili, locchCurvejz, lochhCurvejz, loclhCurvejz, HHcurvejz, CHcurvejz, LHcurvejz, locwavCurvejz, locwavutilijz);
                    }

                    ImProcFunctions::ciecamloc_02float(lp, sp, tmp1.get(), bfw, bfh, 0, sk, cielocalcurve, localcieutili, cielocalcurve2, localcieutili2, jzlocalcurve, localjzutili, czlocalcurve, localczutili, czjzlocalcurve, localczjzutili, locchCurvejz, lochhCurvejz, loclhCurvejz, HHcurvejz, CHcurvejz, LHcurvejz, locwavCurvejz, locwavutilijz);

                    float rad = params->locallab.spots.at(sp).detailcie;
                    loccont(bfw, bfh, tmp1.get(), rad, 15.f, sk);
                }




                    tmp1m->CopyFrom(tmp1.get(), multiThread); //save current result7
                    if(params->locallab.spots.at(sp).equiltm  && params->locallab.spots.at(sp).exptonemap) {
                        if(call == 3) {
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int y = ystart; y < yend; y++) {
                                for (int x = xstart; x < xend; x++) {
                                    savenormtm->L[y][x] = tmp1->L[y - ystart][x - xstart];
                                    savenormtm->a[y][x] = tmp1->a[y - ystart][x - xstart];
                                    savenormtm->b[y][x] = tmp1->b[y - ystart][x - xstart];
                                }
                            }
                        }
                    }
                    bool enatmMasktmap = params->locallab.spots.at(sp).enatmMaskaft;

                    if (enatmMasktmap) {
                        //calculate new values for original, originalmasktm, bufmaskorigtm...in function of tmp1
                        LocwavCurve dummy;
                        int sco = params->locallab.spots.at(sp).scopemask;
                        int shortcu = 0;//lp.mergemet; //params->locallab.spots.at(sp).shortc;

                        const float mindE = 2.f + MINSCOPE * sco * lp.thr;
                        const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
                        const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
                        const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
                        int shado = 0;
                        float amountcd = 0.f;
                        float anchorcd = 50.f;
                        LocHHmaskCurve lochhhmasCurve;
                        const int highl = 0;
                        maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, tmp1.get(), bufmaskorigtm.get(), originalmasktm.get(), original, reserved, inv, lp,
                                    0.f, false,
                                    locccmastmCurve, lcmastmutili, locllmastmCurve, llmastmutili, lochhmastmCurve, lhmastmutili, lochhhmasCurve, false, multiThread,
                                    enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, highl, amountcd, anchorcd, lmasktmlocalcurve, localmasktmutili, dummy, false, 1, 1, 5, 5,
                                    shortcu, params->locallab.spots.at(sp).deltae, hueref, chromaref, lumaref,
                                    maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, -1, fab
                                   );

                        if (lp.showmasktmmet == 3) {//display mask
                            showmask(params->locallab.spots.at(sp).lumask, lp, xstart, ystart, cx, cy, bfw, bfh, tmp1.get(), transformed, bufmaskorigtm.get(), 0);
                            return;
                        }

                    }

                    tmp1->CopyFrom(tmp1m.get(), multiThread); //restore current result


                    float minL = tmp1->L[0][0] - bufgb->L[0][0];
                    float maxL = minL;
                    float minC = std::sqrt(SQR(tmp1->a[0][0]) + SQR(tmp1->b[0][0])) - std::sqrt(SQR(bufgb->a[0][0]) + SQR(bufgb->b[0][0]));
                    float maxC = minC;

#ifdef _OPENMP
                    #pragma omp parallel for reduction(max:maxL) reduction(min:minL) reduction(max:maxC) reduction(min:minC) schedule(dynamic,16) if (multiThread)
#endif
                    for (int ir = 0; ir < bfh; ir++) {
                        for (int jr = 0; jr < bfw; jr++) {
                            buflight[ir][jr] = tmp1->L[ir][jr] - bufgb->L[ir][jr];
                            minL = rtengine::min(minL, buflight[ir][jr]);
                            maxL = rtengine::max(maxL, buflight[ir][jr]);
                            bufchro[ir][jr] = std::sqrt(SQR(tmp1->a[ir][jr]) + SQR(tmp1->b[ir][jr])) - std::sqrt(SQR(bufgb->a[ir][jr]) + SQR(bufgb->b[ir][jr]));
                            minC = rtengine::min(minC, bufchro[ir][jr]);
                            maxC = rtengine::max(maxC, bufchro[ir][jr]);
                        }
                    }

                    float coef = 0.01f * rtengine::max(std::fabs(minL), std::fabs(maxL));
                    float coefC = 0.01f * rtengine::max(std::fabs(minC), std::fabs(maxC));

                    if (coef == 0.f) {
                        coef = 1.f;
                    } else {
                        coef = 1.f /  coef;
                    }

                    if (coefC == 0.f) {
                        coefC = 1.f;
                    } else {
                        coefC = 1.f / coefC;
                    }

#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                    for (int y = 0; y < bfh; y++) {
                        for (int x = 0; x < bfw; x++) {
                            buflight[y][x] *= coef;
                            bufchro[y][x] *= coefC;
                        }
                    }

                    if(lp.enatmMask && lp.recothrt != 1.f) {
                        float recoth = lp.recothrt;

                        if(lp.recothrt < 1.f) {
                            recoth = -1.f * recoth + 2.f;
                        }
                        float hig = lp.higthrt;
                        float low = lp.lowthrt;
                       // float recoth = lp.recothrt;
                        float decay = lp.decayt;
                        bool invmask = false;
                        maskrecov(tmp1.get(), original, bufmaskorigtm.get(), bfh, bfw, ystart, xstart, hig, low, recoth, decay, invmask, sk, multiThread);
                    }

                    //   transit_shapedetect_retinex(call, 4, bufgb.get(),bufmaskorigtm.get(), originalmasktm.get(), buflight, bufchro, hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);

                    if(lp.recothrt >= 1.f) {
                        transit_shapedetect2(sp, meantm, stdtm, call, 8, bufgb.get(), tmp1.get(), originalmasktm.get(), hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);
                    } else {
                        transit_shapedetect2(sp, meantm, stdtm, call, 8, bufgb.get(), tmp1.get(), nullptr, hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);
                    }
                    //  transit_shapedetect(8, tmp1.get(), originalmasktm.get(), bufchro, false, hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);
                    bufgb.reset();

                    if (lp.recur) {
                        original->CopyFrom(transformed, multiThread);
                        float avge;
                        calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                    }
                }
            }
        }
    }

//end TM


    if ((lp.dehaze != 0 || lp.prevdE) && lp.retiena ) {
        int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
        int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
        int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
        int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
        int bfh = yend - ystart;
        int bfw = xend - xstart;

        if (bfh >= mSP && bfw >= mSP) {
            const std::unique_ptr<LabImage> bufexporig(new LabImage(bfw, bfh)); //buffer for data in zone limit
            const std::unique_ptr<LabImage> bufexpfin(new LabImage(bfw, bfh)); //buffer for data in zone limit

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = ystart; y < yend; y++) {
                for (int x = xstart; x < xend; x++) {
                    bufexporig->L[y - ystart][x - xstart] = original->L[y][x];
                    bufexporig->a[y - ystart][x - xstart] = original->a[y][x];
                    bufexporig->b[y - ystart][x - xstart] = original->b[y][x];
                }
            }

            bufexpfin->CopyFrom(bufexporig.get(), multiThread);
            //calc dehaze
            const std::unique_ptr<Imagefloat> tmpImage(new Imagefloat(bfw, bfh));

            DehazeParams dehazeParams;
            dehazeParams.enabled = true;
            dehazeParams.strength = lp.dehaze;
            dehazeParams.showDepthMap = false;
            dehazeParams.saturation = lp.dehazeSaturation;
            dehazeParams.depth = lp.depth;
            lab2rgb(*bufexpfin, *tmpImage.get(), params->icm.workingProfile);
            dehazeloc(tmpImage.get(), dehazeParams);
            rgb2lab(*tmpImage.get(), *bufexpfin, params->icm.workingProfile);

            transit_shapedetect2(sp, 0.f, 0.f, call, 30, bufexporig.get(), bufexpfin.get(), nullptr, hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);

            if (lp.recur) {
                original->CopyFrom(transformed, multiThread);
                float avge;
                calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
            }
        }
    }

    lp.invret = false;//always disabled inverse RETI   too complex todo !!

    if (lp.str >= 0.2f && lp.retiena && call != 2) {
        LabImage *bufreti = nullptr;
        LabImage *bufmask = nullptr;
        LabImage *buforig = nullptr;
        LabImage *buforigmas = nullptr;
        LabImage *bufmaskorigreti = nullptr;

        if (TW >= mSP && TH >= mSP) {

            array2D<float> buflight(TW, TH);
            JaggedArray<float> bufchro(TW, TH);

            int Hd, Wd;
            Hd = TH;
            Wd = TW;

            bufreti = new LabImage(TW, TH);
            bufmask = new LabImage(TW, TH);
                bufmaskorigreti = new LabImage(TW, TH);

            if (!lp.enaretiMasktmap && lp.enaretiMask) {
                buforig = new LabImage(TW, TH);
                buforigmas = new LabImage(TW, TH);
              //  bufmaskorigreti = new LabImage(GW, GH);
            }

#ifdef _OPENMP
            #pragma omp parallel for if (multiThread)
#endif
            for (int ir = 0; ir < TH; ir++) //fill with 0
                for (int jr = 0; jr < TW; jr++) {
                    bufreti->L[ir][jr] = 0.f;
                    bufreti->a[ir][jr] = 0.f;
                    bufreti->b[ir][jr] = 0.f;
                    buflight[ir][jr] = 0.f;
                    bufchro[ir][jr] = 0.f;
                }

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = 0; y < transformed->H ; y++) //{
                for (int x = 0; x < transformed->W; x++) {
                    bufreti->L[y][x] = original->L[y][x];
                    bufreti->a[y][x] = original->a[y][x];
                    bufreti->b[y][x] = original->b[y][x];
                    bufmask->L[y][x] = original->L[y][x];
                    bufmask->a[y][x] = original->a[y][x];
                    bufmask->b[y][x] = original->b[y][x];

                

                    if (!lp.enaretiMasktmap && lp.enaretiMask) {
                        buforig->L[y][x] = original->L[y][x];
                        buforig->a[y][x] = original->a[y][x];
                        buforig->b[y][x] = original->b[y][x];

                     //   bufmaskorigreti->L[y][x] = original->L[y][x];
                     //   bufmaskorigreti->a[y][x] = original->a[y][x];
                     //   bufmaskorigreti->b[y][x] = original->b[y][x];
                        
                        
                    }

                }
            float raddE = params->locallab.spots.at(sp).softradiusret;

            //calc dE and reduction to use in MSR to reduce artifacts
            const float mindE = 4.f + MINSCOPE * lp.sensh * lp.thr;
            const float maxdE = 5.f + MAXSCOPE * lp.sensh * (1 + 0.1f * lp.thr);
            const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
            const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
            const float refa = chromaref * cos(hueref);
            const float refb = chromaref * sin(hueref);

            const std::unique_ptr<JaggedArray<float>> reducDEBuffer(new JaggedArray<float>(Wd, Hd));
            float** reducDE = *reducDEBuffer;

            float ade = 0.01f * raddE;
            float bde = 100.f - raddE;
            float sensibefore = ade * lp.sensh + bde;//we can change sensitivity 0.1 90 or 0.3 70 or 0.4 60
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = 0; y < transformed->H ; y++)
                for (int x = 0; x < transformed->W; x++) {
                    float dE = std::sqrt(SQR(refa - bufreti->a[y][x] / 327.68f) + SQR(refb - bufreti->b[y][x] / 327.68f) + SQR(static_cast<float>(lumaref) - bufreti->b[y][x] / 327.68f));
                    const float reducdE = calcreducdE(dE, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sensibefore);
                    reducDE[y][x] = clipDE(reducdE);
                }

            const std::unique_ptr<JaggedArray<float>> origBuffer(new JaggedArray<float>(Wd, Hd));
            float** orig = *origBuffer;

            const std::unique_ptr<JaggedArray<float>> origBuffer1(new JaggedArray<float>(Wd, Hd));
            float** orig1 = *origBuffer1;

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int ir = 0; ir < Hd; ir += 1)
                for (int jr = 0; jr < Wd; jr += 1) {
                    orig[ir][jr] = bufreti->L[ir][jr];
                    orig1[ir][jr] = bufreti->L[ir][jr];
                }

            LabImage *tmpl = new LabImage(Wd, Hd);

            bool fftw = lp.ftwreti;
            //for Retinex Mask are incorporated in MSR
            int sco = params->locallab.spots.at(sp).scopemask;
            float lumask = params->locallab.spots.at(sp).lumask;

            const float mindE2 = 2.f + MINSCOPE * sco * lp.thr;
            const float maxdE2 = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
            const float mindElim2 = 2.f + MINSCOPE * limscope * lp.thr;
            const float maxdElim2 = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
            ImProcFunctions::MSRLocal(call, sp, fftw, 1, reducDE, bufreti, bufmask, buforig, buforigmas, bufmaskorigreti, orig, orig1,
                                      Wd, Hd, Wd, Hd, params->locallab, sk, locRETgainCcurve, locRETtransCcurve, 0, 4, 1.f, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax,
                                      locccmasretiCurve, lcmasretiutili, locllmasretiCurve, llmasretiutili, lochhmasretiCurve, lhmasretiutili, llretiMask,
                                      lmaskretilocalcurve, localmaskretiutili,
                                      transformed, lp.enaretiMasktmap, lp.enaretiMask,
                                      params->locallab.spots.at(sp).deltae, hueref, chromaref, lumaref,
                                      maxdE2, mindE2, maxdElim2, mindElim2, lp.iterat, limscope, sco, lp.balance, lp.balanceh, lumask);
#ifdef _OPENMP
            #pragma omp parallel for if (multiThread)
#endif
            for (int ir = 0; ir < Hd; ir += 1) {
                for (int jr = 0; jr < Wd; jr += 1) {
                    tmpl->L[ir][jr] = orig[ir][jr];
                    if(params->locallab.spots.at(sp).equilret  && params->locallab.spots.at(sp).expreti) {
                        if(call == 3) {
                            savenormreti->L[ir][jr] = tmpl->L[ir][jr];
                        }
                    }
                }
            }

            if (lp.equret) { //equilibrate luminance before / after MSR
                float *datain = new float[Hd * Wd];
                float *data = new float[Hd * Wd];
#ifdef _OPENMP
                #pragma omp parallel for if (multiThread)
#endif
                for (int ir = 0; ir < Hd; ir += 1)
                    for (int jr = 0; jr < Wd; jr += 1) {
                        datain[ir * Wd + jr] = orig1[ir][jr];
                        data[ir * Wd + jr] = orig[ir][jr];
                    }

                if(params->locallab.spots.at(sp).equilret){
                    if(call == 3) {//improccoordinator
                        normalize_mean_dt(data, datain, Hd * Wd, 1.f, 1.f, 0.f, 0.f, 0.f, 0.f);
                    } else if(call == 1) {//dcrop
                        float ma =  meanreti;
                        float sa = stdreti;
                        float ma2 =  (float) params->locallab.spots.at(sp).sensihs;
                        float sa2 = (float) params->locallab.spots.at(sp).sensiv;
                        //printf("ma=%f sa=%f ma2=%f sa2=%f\n", (double) ma, (double) sa, (double) ma2, (double) sa2);
                        //use normalize with mean and stdv
                        normalize_mean_dt(data, datain, Hd * Wd, 1.f, 1.f, ma, sa, ma2, sa2);

                    }
                }

#ifdef _OPENMP
                #pragma omp parallel for if (multiThread)
#endif
                for (int ir = 0; ir < Hd; ir += 1)
                    for (int jr = 0; jr < Wd; jr += 1) {
                        tmpl->L[ir][jr] = data[ir * Wd + jr];
                    }

                delete [] datain;
                delete [] data;
            }

            if(lp.enaretiMask && lp.recothrr != 1.f) {
                float hig = lp.higthrr;
                float low = lp.lowthrr;
                float recoth = lp.recothrr;
                float decay = lp.decayr;
                bool invmask = false;
                maskrecov(tmpl, original, bufmaskorigreti, Hd, Wd, 0, 0, hig, low, recoth, decay, invmask, sk, multiThread);
            }

            float minL = tmpl->L[0][0] - bufreti->L[0][0];
            float maxL = minL;
#ifdef _OPENMP
            #pragma omp parallel for reduction(min:minL) reduction(max:maxL) schedule(dynamic,16) if (multiThread)
#endif
            for (int ir = 0; ir < Hd; ir++) {
                for (int jr = 0; jr < Wd; jr++) {
                    buflight[ir][jr] = tmpl->L[ir][jr] - bufreti->L[ir][jr];
                    minL = rtengine::min(minL, buflight[ir][jr]);
                    maxL = rtengine::max(maxL, buflight[ir][jr]);
                }
            }

            const float coef = 0.01f * rtengine::max(std::fabs(minL), std::fabs(maxL));

            for (int ir = 0; ir < Hd; ir++) {
                for (int jr = 0; jr < Wd; jr++) {
                    buflight[ir][jr] /= coef;
                }
            }

            transit_shapedetect_retinex(call, 4, bufreti, tmpl, bufmask, buforigmas, buflight, bufchro, hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);

            if (lp.recur) {
                original->CopyFrom(transformed, multiThread);
                float avge;
                calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
            }

            if (params->locallab.spots.at(sp).chrrt > 0) {

                if (call == 1) {

#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                    for (int ir = 0; ir < Hd; ir += 1)
                        for (int jr = 0; jr < Wd; jr += 1) {

                            orig[ir][jr] = std::sqrt(SQR(bufreti->a[ir][jr]) + SQR(bufreti->b[ir][jr]));
                            orig1[ir][jr] = std::sqrt(SQR(bufreti->a[ir][jr]) + SQR(bufreti->b[ir][jr]));
                        }

                }

                float maxChro = orig1[0][0];
#ifdef _OPENMP
                #pragma omp parallel for reduction(max:maxChro) schedule(dynamic,16) if (multiThread)
#endif
                for (int ir = 0; ir < Hd; ir++) {
                    for (int jr = 0; jr < Wd; jr++) {
                        maxChro = rtengine::max(maxChro, orig1[ir][jr]);
                    }
                }

                float divchro = maxChro;

                //first step change saturation without Retinex ==> gain of time and memory
                float satreal = lp.str * static_cast<float>(params->locallab.spots.at(sp).chrrt) / 100.f;

                if (params->locallab.spots.at(sp).chrrt <= 0.2) {
                    satreal /= 10.f;
                }

                DiagonalCurve reti_satur({
                    DCT_NURBS,
                    0, 0,
                    0.2, 0.2f + satreal / 250.f,
                    0.6,  rtengine::min(1.f, 0.6f + satreal / 250.f),
                    1, 1
                });

                if (call == 1) {

#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                    for (int ir = 0; ir < Hd; ir += 1)
                        for (int jr = 0; jr < Wd; jr += 1) {
                            const float Chprov = orig1[ir][jr];
                            float2 sincosval;
                            sincosval.y = Chprov == 0.0f ? 1.f : bufreti->a[ir][jr] / Chprov;
                            sincosval.x = Chprov == 0.0f ? 0.f : bufreti->b[ir][jr] / Chprov;

                            if (params->locallab.spots.at(sp).chrrt <= 100.0) { //first step
                                float buf = LIM01(orig[ir][jr] / divchro);
                                buf = reti_satur.getVal(buf);
                                buf *= divchro;
                                orig[ir][jr] = buf;
                            }

                            tmpl->a[ir][jr] = orig[ir][jr] * sincosval.y;
                            tmpl->b[ir][jr] = orig[ir][jr] * sincosval.x;
                        }

                    float minC = std::sqrt(SQR(tmpl->a[0][0]) + SQR(tmpl->b[0][0])) - orig1[0][0];
                    float maxC = minC;
#ifdef _OPENMP
                    #pragma omp parallel for reduction(min:minC) reduction(max:maxC) schedule(dynamic,16) if (multiThread)
#endif
                    for (int ir = 0; ir < Hd; ir++) {
                        for (int jr = 0; jr < Wd; jr++) {
                            bufchro[ir][jr] = std::sqrt(SQR(tmpl->a[ir][jr]) + SQR(tmpl->b[ir][jr])) - orig1[ir][jr];
                            minC = rtengine::min(minC, bufchro[ir][jr]);
                            maxC = rtengine::max(maxC, bufchro[ir][jr]);
                        }
                    }

                    float coefC = 0.01f * rtengine::max(std::fabs(minC), std::fabs(maxC));

                    if (coefC > 0.f) {
                        coefC = 1.f / coefC;
#ifdef _OPENMP
                        #pragma omp parallel for if (multiThread)
#endif
                        for (int ir = 0; ir < Hd; ir++) {
                            for (int jr = 0; jr < Wd; jr++) {
                                bufchro[ir][jr] *= coefC;
                            }
                        }
                    }
                }

                transit_shapedetect_retinex(call, 5, tmpl, tmpl, bufmask, buforigmas, buflight, bufchro, hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);

                if (lp.recur) {
                    original->CopyFrom(transformed, multiThread);
                    float avge;
                    calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                }
            }

            delete tmpl;
            delete bufmask;
            delete bufmaskorigreti;

            if (!lp.enaretiMasktmap && lp.enaretiMask) {
                if (buforig) {
                    delete buforig;
                }

                if (buforigmas) {
                    delete buforigmas;
                }
            }
            delete  bufreti;
        }
    }



    if (lp.str >= 0.2f && lp.retiena && call == 2) {
        int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
        int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
        int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
        int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
        int bfh = yend - ystart;
        int bfw = xend - xstart;

        LabImage *bufreti = nullptr;
        LabImage *bufmask = nullptr;
        LabImage *buforig = nullptr;
        LabImage *buforigmas = nullptr;
        LabImage *bufmaskorigreti = nullptr;
        int bfhr = bfh;
        int bfwr = bfw;

        if (bfw >= mSP && bfh > mSP) {
            if (lp.ftwreti) {
                optfft(N_fftwsize, bfh, bfw, bfhr, bfwr, lp, original->H, original->W, xstart, ystart, xend, yend, cx, cy, lp.fullim);
            }

            array2D<float> buflight(bfw, bfh);
            JaggedArray<float> bufchro(bfw, bfh);

            int Hd, Wd;
            Hd = TH;
            Wd = TW;

            if (!lp.invret && call == 2) {

                Hd = bfh;
                Wd = bfw;
                bufreti = new LabImage(bfw, bfh);
                bufmask = new LabImage(bfw, bfh);
                bufmaskorigreti = new LabImage(bfw, bfh);

                if (!lp.enaretiMasktmap && lp.enaretiMask) {
                    buforig = new LabImage(bfw, bfh);
                    buforigmas = new LabImage(bfw, bfh);
                }

#ifdef _OPENMP
                #pragma omp parallel for if (multiThread)
#endif
                for (int ir = 0; ir < bfh; ir++) //fill with 0
                    for (int jr = 0; jr < bfw; jr++) {
                        bufreti->L[ir][jr] = 0.f;
                        bufreti->a[ir][jr] = 0.f;
                        bufreti->b[ir][jr] = 0.f;
                        buflight[ir][jr] = 0.f;
                        bufchro[ir][jr] = 0.f;
                    }


#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                for (int y = ystart; y < yend; y++) {
                    for (int x = xstart; x < xend; x++) {
                        bufreti->L[y - ystart][x - xstart] = original->L[y][x];
                        bufreti->a[y - ystart][x - xstart] = original->a[y][x];
                        bufreti->b[y - ystart][x - xstart] = original->b[y][x];
                        bufmask->L[y - ystart][x - xstart] = original->L[y][x];
                        bufmask->a[y - ystart][x - xstart] = original->a[y][x];
                        bufmask->b[y - ystart][x - xstart] = original->b[y][x];

                        if (!lp.enaretiMasktmap && lp.enaretiMask) {
                            buforig->L[y - ystart][x - xstart] = original->L[y][x];
                            buforig->a[y - ystart][x - xstart] = original->a[y][x];
                            buforig->b[y - ystart][x - xstart] = original->b[y][x];
                        }
                    }
                }
            }

            float raddE = params->locallab.spots.at(sp).softradiusret;

            //calc dE and reduction to use in MSR to reduce artifacts
            const float mindE = 4.f + MINSCOPE * lp.sensh * lp.thr;
            const float maxdE = 5.f + MAXSCOPE * lp.sensh * (1 + 0.1f * lp.thr);
            const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
            const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
            const float refa = chromaref * cos(hueref);
            const float refb = chromaref * sin(hueref);

            const std::unique_ptr<JaggedArray<float>> reducDEBuffer(new JaggedArray<float>(Wd, Hd));
            float** reducDE = *reducDEBuffer;
            float ade = 0.01f * raddE;
            float bde = 100.f - raddE;
            float sensibefore = ade * lp.sensh + bde;//we can change sensitivity 0.1 90 or 0.3 70 or 0.4 60
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = ystart; y < yend ; y++) {
                for (int x = xstart; x < xend; x++) {
                    const float dE = std::sqrt(SQR(refa - bufreti->a[y - ystart][x - xstart] / 327.68f) + SQR(refb - bufreti->b[y - ystart][x - xstart] / 327.68f) + SQR(static_cast<float>(lumaref) - bufreti->b[y - ystart][x - xstart] / 327.68f));
                    const float reducdE = calcreducdE(dE, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sensibefore);
                    reducDE[y - ystart][x - xstart] = clipDE(reducdE);
                }
            }

            const std::unique_ptr<JaggedArray<float>> origBuffer(new JaggedArray<float>(Wd, Hd));
            float** orig = *origBuffer;

            const std::unique_ptr<JaggedArray<float>> origBuffer1(new JaggedArray<float>(Wd, Hd));
            float** orig1 = *origBuffer1;

            LabImage *tmpl = nullptr;

            if (!lp.invret && call == 2) {
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                for (int ir = 0; ir < Hd; ir += 1) {
                    for (int jr = 0; jr < Wd; jr += 1) {
                        orig[ir][jr] = bufreti->L[ir][jr];
                        orig1[ir][jr] = bufreti->L[ir][jr];
                    }
                }

                tmpl = new LabImage(Wd, Hd);
            }

            //   float minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax;
            bool fftw = lp.ftwreti;
            //for Retinex Mask are incorporated in MSR
            int sco = params->locallab.spots.at(sp).scopemask;
            float lumask = params->locallab.spots.at(sp).lumask;

            const float mindE2 = 2.f + MINSCOPE * sco * lp.thr;
            const float maxdE2 = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
            const float mindElim2 = 2.f + MINSCOPE * limscope * lp.thr;
            const float maxdElim2 = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);

            ImProcFunctions::MSRLocal(call, sp, fftw, 1, reducDE, bufreti, bufmask, buforig, buforigmas, bufmaskorigreti, orig, orig1,
                                      Wd, Hd, bfwr, bfhr, params->locallab, sk, locRETgainCcurve, locRETtransCcurve, 0, 4, 1.f, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax,
                                      locccmasretiCurve, lcmasretiutili, locllmasretiCurve, llmasretiutili, lochhmasretiCurve, lhmasretiutili, llretiMask,
                                      lmaskretilocalcurve, localmaskretiutili,
                                      transformed, lp.enaretiMasktmap, lp.enaretiMask,
                                      params->locallab.spots.at(sp).deltae, hueref, chromaref, lumaref,
                                      maxdE2, mindE2, maxdElim2, mindElim2, lp.iterat, limscope, sco, lp.balance, lp.balanceh, lumask);

#ifdef _OPENMP
            #pragma omp parallel for if (multiThread)
#endif
            for (int ir = 0; ir < Hd; ir += 1)
                for (int jr = 0; jr < Wd; jr += 1) {
                    tmpl->L[ir][jr] = orig[ir][jr];
                }


            if (lp.equret) { //equilibrate luminance before / after MSR
                const std::unique_ptr<float[]> datain(new float[Hd * Wd]);
                const std::unique_ptr<float[]> data(new float[Hd * Wd]);
#ifdef _OPENMP
                #pragma omp parallel for if (multiThread)
#endif
                for (int ir = 0; ir < Hd; ir += 1) {
                    for (int jr = 0; jr < Wd; jr += 1) {
                        datain[ir * Wd + jr] = orig1[ir][jr];
                        data[ir * Wd + jr] = orig[ir][jr];
                    }
                }

                normalize_mean_dt(data.get(), datain.get(), Hd * Wd, 1.f, 1.f, 0.f, 0.f, 0.f, 0.f);
#ifdef _OPENMP
                #pragma omp parallel for if (multiThread)
#endif
                for (int ir = 0; ir < Hd; ir += 1) {
                    for (int jr = 0; jr < Wd; jr += 1) {
                        tmpl->L[ir][jr] = data[ir * Wd + jr];
                    }
                }
            }
            if(lp.enaretiMask && lp.recothrr != 1.f) {
                float hig = lp.higthrr;
                float low = lp.lowthrr;
                float recoth = lp.recothrr;
                float decay = lp.decayr;
                bool invmask = false;
                maskrecov(tmpl, original, bufmaskorigreti, Hd, Wd, ystart, xstart, hig, low, recoth, decay, invmask, sk, multiThread);
            }

            if (!lp.invret) {
                float minL = tmpl->L[0][0] - bufreti->L[0][0];
                float maxL = minL;
#ifdef _OPENMP
                #pragma omp parallel for reduction(min:minL) reduction(max:maxL) schedule(dynamic,16) if (multiThread)
#endif
                for (int ir = 0; ir < Hd; ir++) {
                    for (int jr = 0; jr < Wd; jr++) {
                        buflight[ir][jr] = tmpl->L[ir][jr] - bufreti->L[ir][jr];
                        minL = rtengine::min(minL, buflight[ir][jr]);
                        maxL = rtengine::max(maxL, buflight[ir][jr]);
                    }
                }

                float coef = 0.01f * rtengine::max(std::fabs(minL), std::fabs(maxL));

                if (coef > 0.f) {
                    coef = 1.f / coef;
#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                    for (int ir = 0; ir < Hd; ir++) {
                        for (int jr = 0; jr < Wd; jr++) {
                            buflight[ir][jr] *= coef;
                        }
                    }
                }

                transit_shapedetect_retinex(call, 4, bufreti, tmpl, bufmask, buforigmas, buflight, bufchro, hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);

                if (lp.recur) {
                    original->CopyFrom(transformed, multiThread);
                    float avge;
                    calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                }
            }

            if (params->locallab.spots.at(sp).chrrt > 0) {
                if (!lp.invret && call == 2) {
#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                    for (int ir = 0; ir < Hd; ir += 1) {
                        for (int jr = 0; jr < Wd; jr += 1) {
                            orig[ir][jr] = std::sqrt(SQR(bufreti->a[ir][jr]) + SQR(bufreti->b[ir][jr]));
                            orig1[ir][jr] = std::sqrt(SQR(bufreti->a[ir][jr]) + SQR(bufreti->b[ir][jr]));
                        }
                    }
                }

                float maxChro = orig1[0][0];
#ifdef _OPENMP
                #pragma omp parallel for reduction(max:maxChro) schedule(dynamic,16) if (multiThread)
#endif
                for (int ir = 0; ir < Hd; ir++) {
                    for (int jr = 0; jr < Wd; jr++) {
                        maxChro = rtengine::max(maxChro, orig1[ir][jr]);
                    }
                }

                //first step change saturation without Retinex ==> gain of time and memory
                float satreal = lp.str * static_cast<float>(params->locallab.spots.at(sp).chrrt) / 100.f;

                if (params->locallab.spots.at(sp).chrrt <= 0.2) {
                    satreal /= 10.f;
                }

                DiagonalCurve reti_satur({
                    DCT_NURBS,
                    0, 0,
                    0.2, 0.2f + satreal / 250.f,
                    0.6, rtengine::min(1.f, 0.6f + satreal / 250.f),
                    1, 1
                });

                if (!lp.invret && call == 2) {

#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif
                    for (int ir = 0; ir < Hd; ir += 1) {
                        for (int jr = 0; jr < Wd; jr += 1) {
                            const float Chprov = orig1[ir][jr];
                            float2 sincosval;
                            sincosval.y = Chprov == 0.0f ? 1.f : bufreti->a[ir][jr] / Chprov;
                            sincosval.x = Chprov == 0.0f ? 0.f : bufreti->b[ir][jr] / Chprov;

                            if (params->locallab.spots.at(sp).chrrt <= 40.0) { //first step
                                orig[ir][jr] = static_cast<float>(reti_satur.getVal(LIM01(orig[ir][jr] / maxChro))) * maxChro;
                            }

                            tmpl->a[ir][jr] = orig[ir][jr] * sincosval.y;
                            tmpl->b[ir][jr] = orig[ir][jr] * sincosval.x;
                        }
                    }

                    float minC = std::sqrt(SQR(tmpl->a[0][0]) + SQR(tmpl->b[0][0])) - orig1[0][0];
                    float maxC = minC;
#ifdef _OPENMP
                    #pragma omp parallel for reduction(min:minC) reduction(max:maxC) schedule(dynamic,16) if (multiThread)
#endif
                    for (int ir = 0; ir < Hd; ir++) {
                        for (int jr = 0; jr < Wd; jr++) {
                            bufchro[ir][jr] = std::sqrt(SQR(tmpl->a[ir][jr]) + SQR(tmpl->b[ir][jr])) - orig1[ir][jr];
                            minC = rtengine::min(minC, bufchro[ir][jr]);
                            maxC = rtengine::max(maxC, bufchro[ir][jr]);
                        }
                    }

                    float coefC = 0.01f * rtengine::max(std::fabs(minC), std::fabs(maxC));

                    if (coefC > 0.f) {
                        coefC = 1.f / coefC;
#ifdef _OPENMP
                        #pragma omp parallel for if (multiThread)
#endif
                        for (int ir = 0; ir < Hd; ir++) {
                            for (int jr = 0; jr < Wd; jr++) {
                                bufchro[ir][jr] *= coefC;
                            }
                        }
                    }
                }

                if (!lp.invret) {
                    transit_shapedetect_retinex(call, 5, tmpl, tmpl, bufmask, buforigmas, buflight, bufchro, hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);

                    if (lp.recur) {
                        original->CopyFrom(transformed, multiThread);
                        float avge;
                        calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                    }
                }
            }

            delete tmpl;
            delete bufmask;
            delete bufmaskorigreti;
            
            if (!lp.enaretiMasktmap && lp.enaretiMask) {
                if (buforig) {
                    delete buforig;
                }

                if (buforigmas) {
                    delete buforigmas;
                }
            }
            delete  bufreti;
        }
    }


//begin cbdl
    if ((lp.mulloc[0] != 1.f || lp.mulloc[1] != 1.f || lp.mulloc[2] != 1.f || lp.mulloc[3] != 1.f || lp.mulloc[4] != 1.f || lp.mulloc[5] != 1.f || lp.clarityml != 0.f || lp.contresid != 0.f  || lp.enacbMask || lp.showmaskcbmet == 2 || lp.showmaskcbmet == 3 || lp.showmaskcbmet == 4 || lp.prevdE) && lp.cbdlena) {
        if (call <= 3) { //call from simpleprocess dcrop improcc
            const int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
            const int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
            const int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
            const int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
            int bfh = yend - ystart;
            int bfw = xend - xstart;
            if (bfw > 65 && bfh > 65) {
                array2D<float> bufsh(bfw, bfh);
                JaggedArray<float> bufchrom(bfw, bfh, true);
                const std::unique_ptr<LabImage> loctemp(new LabImage(bfw, bfh));
                const std::unique_ptr<LabImage> origcbdl(new LabImage(bfw, bfh));
                std::unique_ptr<LabImage> bufmaskorigcb;
                std::unique_ptr<LabImage> bufmaskblurcb;
                std::unique_ptr<LabImage> originalmaskcb;

                if (lp.showmaskcbmet == 2  || lp.enacbMask || lp.showmaskcbmet == 3 || lp.showmaskcbmet == 4) {
                    bufmaskorigcb.reset(new LabImage(bfw, bfh));
                    bufmaskblurcb.reset(new LabImage(bfw, bfh));
                    originalmaskcb.reset(new LabImage(bfw, bfh));
                }

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                for (int y = ystart; y < yend; y++) {
                    for (int x = xstart; x < xend; x++) {
                        loctemp->L[y - ystart][x - xstart] = original->L[y][x];
                    }
                }

                int inv = 0;
                bool showmaske = false;
                bool enaMask = false;
                bool deltaE = false;
                bool modmask = false;
                bool zero = false;
                bool modif = false;

                if (lp.showmaskcbmet == 3) {
                    showmaske = true;
                }

                if (lp.enacbMask) {
                    enaMask = true;
                }

                if (lp.showmaskcbmet == 4) {
                    deltaE = true;
                }

                if (lp.showmaskcbmet == 2) {
                    modmask = true;
                }

                if (lp.showmaskcbmet == 1) {
                    modif = true;
                }

                if (lp.showmaskcbmet == 0) {
                    zero = true;
                }

                float chrom = lp.chromacbm;;
                float rad = lp.radmacb;
                float gamma = lp.gammacb;
                float slope = lp.slomacb;
                float blendm = lp.blendmacb;
                float lap = params->locallab.spots.at(sp).lapmaskcb;
                bool pde = params->locallab.spots.at(sp).laplac;
                LocwavCurve dummy;
                int sco = params->locallab.spots.at(sp).scopemask;
                int lumask = params->locallab.spots.at(sp).lumask;
                int shado = 0;
                const float mindE = 2.f + MINSCOPE * sco * lp.thr;
                const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
                const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
                const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
                float amountcd = 0.f;
                float anchorcd = 50.f;
                int shortcu = 0; //lp.mergemet; //params->locallab.spots.at(sp).shortc;
                LocHHmaskCurve lochhhmasCurve;
                const int highl = 0;
                maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, loctemp.get(), bufmaskorigcb.get(), originalmaskcb.get(), original, reserved, inv, lp,
                            0.f, false,
                            locccmascbCurve, lcmascbutili, locllmascbCurve, llmascbutili, lochhmascbCurve, lhmascbutili, lochhhmasCurve, false, multiThread,
                            enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, highl, amountcd, anchorcd, lmaskcblocalcurve, localmaskcbutili, dummy, false, 1, 1, 5, 5,
                            shortcu, params->locallab.spots.at(sp).deltae, hueref, chromaref, lumaref,
                            maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.0f, 0.f, -1, fab
                           );

                if (lp.showmaskcbmet == 3) {
                    showmask(lumask, lp, xstart, ystart, cx, cy, bfw, bfh, loctemp.get(), transformed, bufmaskorigcb.get(), 0);

                    return;
                }

                constexpr float b_l = -5.f;
                constexpr float t_l = 25.f;
                constexpr float t_r = 120.f;
                constexpr float b_r = 170.f;
                constexpr double skinprot = 0.;
                int choice = 0;

                if (lp.showmaskcbmet == 0 || lp.showmaskcbmet == 1  || lp.showmaskcbmet == 2 || lp.showmaskcbmet == 4 || lp.enacbMask) {

#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                    for (int y = ystart; y < yend; y++) {
                        for (int x = xstart; x < xend; x++) {
                            bufsh[y - ystart][x - xstart] = origcbdl->L[y - ystart][x - xstart] = original->L[y][x];
                            loctemp->a[y - ystart][x - xstart] = origcbdl->a[y - ystart][x - xstart] = original->a[y][x];
                            loctemp->b[y - ystart][x - xstart] = origcbdl->b[y - ystart][x - xstart] = original->b[y][x];
                            loctemp->L[y - ystart][x - xstart] = origcbdl->b[y - ystart][x - xstart] = original->L[y][x];
                        }
                    }

                    if (lp.clarityml != 0.f && lp.mulloc[5] == 1.f) { //enabled last level to retrieve level 5 and residual image in case user not select level 5
                        lp.mulloc[5] = 1.001f;
                    }

                    if (lp.contresid != 0.f && lp.mulloc[5] == 1.f) { //enabled last level to retrieve level 5 and residual image in case user not select level 5
                        lp.mulloc[5] = 1.001f;
                    }

                    ImProcFunctions::cbdl_local_temp(bufsh, loctemp->L, bfw, bfh, lp.mulloc, 1.f, lp.threshol, lp.clarityml, lp.contresid, skinprot, false, b_l, t_l, t_r, b_r, choice, sk, multiThread);

                    
                    if (lp.softradiuscb > 0.f) {
                        softproc(origcbdl.get(), loctemp.get(), lp.softradiuscb, bfh, bfw, 0.001, 0.00001, 0.5f, sk, multiThread, 1);
                    }
                    
                    if(lp.enacbMask && lp.recothrcb != 1.f) {
                        float hig = lp.higthrcb;
                        float low = lp.lowthrcb;
                        float recoth = lp.recothrcb;
                        float decay = lp.decaycb;
                        bool invmask = false;
                        maskrecov(loctemp.get(), original, bufmaskorigcb.get(), bfh, bfw, ystart, xstart, hig, low, recoth, decay, invmask, sk, multiThread);
                    }
                    
                }

                transit_shapedetect(6, loctemp.get(), originalmaskcb.get(), bufchrom, false, hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);

                const bool nochroma = (lp.showmaskcbmet == 2 || lp.showmaskcbmet == 1);

                //chroma CBDL begin here
                if (lp.chromacb > 0.f && !nochroma) {

#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                    for (int ir = 0; ir < bfh; ir++) {
                        for (int jr = 0; jr < bfw; jr++) {
                            bufsh[ir][jr] = std::sqrt(SQR(loctemp->a[ir][jr]) + SQR(loctemp->b[ir][jr]));
                        }
                    }

                    float multc[6];
                    float clarich = 0.5f * lp.clarityml;

                    if (clarich > 0.f && lp.mulloc[0] == 1.f) { //to enabled in case of user select only clarity
                        lp.mulloc[0] = 1.01f;
                    }

                    if (lp.contresid != 0.f && lp.mulloc[0] == 1.f) { //to enabled in case of user select only clarity
                        lp.mulloc[0] = 1.01f;
                    }

                    for (int lv = 0; lv < 6; lv++) {
                        multc[lv] = rtengine::max((lp.chromacb * (lp.mulloc[lv] - 1.f)) + 1.f, 0.01f);
                    }

                    choice = 1;
                    ImProcFunctions::cbdl_local_temp(bufsh, loctemp->L, bfw, bfh, multc, rtengine::max(lp.chromacb, 1.f), lp.threshol, clarich, 0.f, skinprot, false,  b_l, t_l, t_r, b_r, choice, sk, multiThread);

                    float minC = loctemp->L[0][0] - std::sqrt(SQR(loctemp->a[0][0]) + SQR(loctemp->b[0][0]));
                    float maxC = minC;
#ifdef _OPENMP
                    #pragma omp parallel for reduction(max:maxC) reduction(min:minC) schedule(dynamic,16) if (multiThread)
#endif
                    for (int ir = 0; ir < bfh; ir++) {
                        for (int jr = 0; jr < bfw; jr++) {
                            bufchrom[ir][jr] = (loctemp->L[ir][jr] - std::sqrt(SQR(loctemp->a[ir][jr]) + SQR(loctemp->b[ir][jr])));
                            minC = rtengine::min(minC, bufchrom[ir][jr]);
                            maxC = rtengine::max(maxC, bufchrom[ir][jr]);
                        }
                    }

                    float coefC = 0.01f * rtengine::max(std::fabs(minC), std::fabs(maxC));

                    if (coefC > 0.f) {
                        coefC = 1.f / coefC;
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int ir = 0; ir < bfh; ir++) {
                            for (int jr = 0; jr < bfw; jr++) {
                                bufchrom[ir][jr] *= coefC;
                            }
                        }
                    }

                    transit_shapedetect(7, loctemp.get(), nullptr, bufchrom, false, hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);
                }
                bufsh.free();
                
                if (lp.recur) {
                    original->CopyFrom(transformed, multiThread);
                    float avge;
                    calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                }
                
            }
        }
    }


//end cbdl_Local

//vibrance
    float vibg = params->locallab.spots.at(sp).vibgam;
    if (lp.expvib && (lp.past != 0.f  || lp.satur != 0.f || lp.strvib != 0.f || vibg != 1.f  || lp.war != 0 || lp.strvibab != 0.f  || lp.strvibh != 0.f || lp.showmaskvibmet == 2 || lp.enavibMask || lp.showmaskvibmet == 3 || lp.showmaskvibmet == 4 || lp.prevdE) && lp.vibena) { //interior ellipse reinforced lightness and chroma  //locallutili
        if (call <= 3) { //simpleprocess, dcrop, improccoordinator
            const int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
            const int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
            const int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
            const int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
            const int bfh = yend - ystart;
            const int bfw = xend - xstart;

            if (bfw >= mSP && bfh >= mSP) {
                const std::unique_ptr<LabImage> bufexporig(new LabImage(bfw, bfh));
                const std::unique_ptr<LabImage> bufexpfin(new LabImage(bfw, bfh));
                std::unique_ptr<LabImage> bufmaskorigvib;
                std::unique_ptr<LabImage> bufmaskblurvib;
                std::unique_ptr<LabImage> originalmaskvib;

                if (lp.showmaskvibmet == 2  || lp.enavibMask || lp.showmaskvibmet == 3 || lp.showmaskvibmet == 4) {
                    bufmaskorigvib.reset(new LabImage(bfw, bfh));
                    bufmaskblurvib.reset(new LabImage(bfw, bfh));
                    originalmaskvib.reset(new LabImage(bfw, bfh));
                }

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                for (int y = 0; y < bfh; y++) {
                    for (int x = 0; x < bfw; x++) {
                        bufexporig->L[y][x] = original->L[y + ystart][x + xstart];
                    }
                }

                int inv = 0;
                bool showmaske = false;
                bool enaMask = false;
                bool deltaE = false;
                bool modmask = false;
                bool zero = false;
                bool modif = false;

                if (lp.showmaskvibmet == 3) {
                    showmaske = true;
                }

                if (lp.enavibMask) {
                    enaMask = true;
                }

                if (lp.showmaskvibmet == 4) {
                    deltaE = true;
                }

                if (lp.showmaskvibmet == 2) {
                    modmask = true;
                }

                if (lp.showmaskvibmet == 1) {
                    modif = true;
                }

                if (lp.showmaskvibmet == 0) {
                    zero = true;
                }

                float chrom = lp.chromavib;
                float rad = lp.radmavib;
                float gamma = lp.gammavib;
                float slope = lp.slomavib;
                float blendm = lp.blendmavib;
                float lap = params->locallab.spots.at(sp).lapmaskvib;
                bool pde = params->locallab.spots.at(sp).laplac;
                LocwavCurve dummy;
                int sco = params->locallab.spots.at(sp).scopemask;
                int shortcu = 0;//lp.mergemet; //params->locallab.spots.at(sp).shortc;

                const float mindE = 2.f + MINSCOPE * sco * lp.thr;
                const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
                const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
                const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
                int shado = 0;
                int lumask = params->locallab.spots.at(sp).lumask;
                LocHHmaskCurve lochhhmasCurve;
                float amountcd = 0.f;
                float anchorcd = 50.f;
                const int highl = 0;
                maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, bufexporig.get(), bufmaskorigvib.get(), originalmaskvib.get(), original, reserved, inv, lp,
                            0.f, false,
                            locccmasvibCurve, lcmasvibutili, locllmasvibCurve, llmasvibutili, lochhmasvibCurve, lhmasvibutili, lochhhmasCurve, false, multiThread,
                            enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, highl, amountcd, anchorcd, lmaskviblocalcurve, localmaskvibutili, dummy, false, 1, 1, 5, 5,
                            shortcu, params->locallab.spots.at(sp).deltae, hueref, chromaref, lumaref,
                            maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, -1, fab
                           );

                if (lp.showmaskvibmet == 3) {
                    showmask(lumask, lp, xstart, ystart, cx, cy, bfw, bfh, bufexporig.get(), transformed, bufmaskorigvib.get(), 0);

                    return;
                }

                if (lp.showmaskvibmet == 0 || lp.showmaskvibmet == 1  || lp.showmaskvibmet == 2 || lp.showmaskvibmet == 4 || lp.enavibMask) {

#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                    /*
                    for (int y = ystart; y < yend; y++) {
                        for (int x = xstart; x < xend; x++) {
                            bufexporig->L[y - ystart][x - xstart] = original->L[y][x];
                            bufexporig->a[y - ystart][x - xstart] = original->a[y][x];
                            bufexporig->b[y - ystart][x - xstart] = original->b[y][x];
                        }
                    }
                    */
                for (int y = 0; y < bfh; y++) {
                    for (int x = 0; x < bfw; x++) {
                       // bufexporig->L[y][x] = original->L[y + ystart][x + xstart];
                        bufexporig->a[y][x] = original->a[y + ystart][x + xstart];
                        bufexporig->b[y][x] = original->b[y + ystart][x + xstart];
                    }
                }

                    VibranceParams vibranceParams;
                    vibranceParams.enabled = params->locallab.spots.at(sp).expvibrance;
                    vibranceParams.pastels = params->locallab.spots.at(sp).pastels;
                    vibranceParams.saturated = params->locallab.spots.at(sp).saturated;
                    vibranceParams.psthreshold = params->locallab.spots.at(sp).psthreshold;
                    vibranceParams.protectskins = params->locallab.spots.at(sp).protectskins;
                    vibranceParams.avoidcolorshift = params->locallab.spots.at(sp).avoidcolorshift;
                    vibranceParams.pastsattog = params->locallab.spots.at(sp).pastsattog;
                    vibranceParams.skintonescurve = params->locallab.spots.at(sp).skintonescurve;


                //    bufexpfin->CopyFrom(bufexporig.get(), multiThread);
                for (int y = 0; y < bfh; y++) {
                    for (int x = 0; x < bfw; x++) {
                        bufexpfin->L[y][x] = bufexporig->L[y][x];
                         bufexpfin->a[y][x] = bufexporig->a[y][x];
                        bufexpfin->b[y][x] = bufexporig->b[y][x];
                   }
                }

                    if (lp.strvibh != 0.f) {
                        printf("a\n");
                        struct grad_params gph;
                        calclocalGradientParams(lp, gph, ystart, xstart, bfw, bfh, 9);
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int ir = 0; ir < bfh; ir++)
                            for (int jr = 0; jr < bfw; jr++) {
                                float factor = ImProcFunctions::calcGradientFactor(gph, jr, ir);
                                float aa = bufexpfin->a[ir][jr];
                                float bb = bufexpfin->b[ir][jr];
                                float chrm = std::sqrt(SQR(aa) + SQR(bb));
                                float HH = xatan2f(bb, aa);

                                float newhr = 0.f;
                                float cor = 0.f;

                                if (factor < 1.f) {
                                    cor = - 2.5f * (1.f - factor);
                                } else if (factor > 1.f) {
                                    cor = 0.03f * (factor - 1.f);
                                }

                                newhr = HH + cor;

                                if (newhr > rtengine::RT_PI_F) {
                                    newhr -= 2 * rtengine::RT_PI_F;
                                } else if (newhr < -rtengine::RT_PI_F) {
                                    newhr += 2 * rtengine::RT_PI_F;
                                }

                                float2 sincosval = xsincosf(newhr);
                                bufexpfin->a[ir][jr] = clipC(chrm * sincosval.y);
                                bufexpfin->b[ir][jr] = clipC(chrm * sincosval.x);
                            }
                    }

                    if (lp.strvib != 0.f) {
                        printf("b\n");
                        
                        struct grad_params gp;
                        calclocalGradientParams(lp, gp, ystart, xstart, bfw, bfh, 7);
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int ir = 0; ir < bfh; ir++) {
                            for (int jr = 0; jr < bfw; jr++) {
                                bufexpfin->L[ir][jr] *= ImProcFunctions::calcGradientFactor(gp, jr, ir);
                            }
                        }
                    }

                    if (lp.strvibab != 0.f) {
                        printf("c\n");
                        
                        struct grad_params gpab;
                        calclocalGradientParams(lp, gpab, ystart, xstart, bfw, bfh, 8);
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int ir = 0; ir < bfh; ir++)
                            for (int jr = 0; jr < bfw; jr++) {
                                const float factor = ImProcFunctions::calcGradientFactor(gpab, jr, ir);
                                bufexpfin->a[ir][jr] *= factor;
                                bufexpfin->b[ir][jr] *= factor;
                            }
                    }
                    float gamma1 = params->locallab.spots.at(sp).vibgam;
                    rtengine::GammaValues g_a; //gamma parameters
                    double pwr1 = 1.0 / (double) gamma1;//default 3.0 - gamma Lab
                    double ts1 = 9.03296;//always the same 'slope' in the extreme shadows - slope Lab
                    rtengine::Color::calcGamma(pwr1, ts1, g_a); // call to calcGamma with selected gamma and slope
                    if(gamma1 != 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < bfh; ++y) {
                        int x = 0;
#ifdef __SSE2__
                            for (; x < bfw - 3; x += 4) {
                                STVFU(bufexpfin->L[y][x], F2V(32768.f) * igammalog(LVFU(bufexpfin->L[y][x]) / F2V(32768.f), F2V(gamma1), F2V(ts1), F2V(g_a[2]), F2V(g_a[4])));
                            }
#endif
                            for (;x < bfw; ++x) {
                                bufexpfin->L[y][x] = 32768.f * igammalog(bufexpfin->L[y][x] / 32768.f, gamma1, ts1, g_a[2], g_a[4]);
                            }
                        }
                    }


                    ImProcFunctions::vibrance(bufexpfin.get(), vibranceParams, params->toneCurve.hrenabled, params->icm.workingProfile);

                    // float gamma =  params->locallab.spots.at(sp).vibgam;
                    //  rtengine::GammaValues g_a; //gamma parameters
                    // double pwr = 1.0 / (double) gamma;//default 3.0 - gamma Lab
                    // double ts = 9.03296;//always the same 'slope' in the extreme shadows - slope Lab
                   // rtengine::Color::calcGamma(pwr, ts, g_a); // call to calcGamma with selected gamma and slope

                    if(gamma1 != 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif  
                        for (int y = 0; y < bfh; ++y) {//apply inverse gamma 3.f and put result in range 32768.f
                            int x = 0;
#ifdef __SSE2__
                            for (; x < bfw - 3; x += 4) {
                                STVFU(bufexpfin->L[y][x], F2V(32768.f) * gammalog(LVFU(bufexpfin->L[y][x]) / F2V(32768.f), F2V(gamma1), F2V(ts1), F2V(g_a[3]), F2V(g_a[4])));
                            }
#endif
                            for (; x < bfw; ++x) {
                                bufexpfin->L[y][x] = 32768.f * gammalog(bufexpfin->L[y][x] / 32768.f, gamma1, ts1, g_a[3], g_a[4]);
                            }
                        }
                    }


                    if (params->locallab.spots.at(sp).warm != 0) {
                        bool HHcurvejz = false, CHcurvejz = false, LHcurvejz = false;
                       
                        ImProcFunctions::ciecamloc_02float(lp, sp, bufexpfin.get(), bfw, bfh, 2, sk, cielocalcurve, localcieutili, cielocalcurve2, localcieutili2, jzlocalcurve, localjzutili, czlocalcurve, localczutili, czjzlocalcurve, localczjzutili, locchCurvejz, lochhCurvejz, loclhCurvejz, HHcurvejz, CHcurvejz, LHcurvejz, locwavCurvejz, locwavutilijz);
                    }

                    if(lp.enavibMask && lp.recothrv != 1.f) {
                        float recoth = lp.recothrv;

                        if(lp.recothrv < 1.f) {
                            recoth = -1.f * recoth + 2.f;
                        }

                        float hig = lp.higthrv;
                        float low = lp.lowthrv;
                      //  float recoth = lp.recothrv;
                        float decay = lp.decayv;
                        bool invmask = false;
                        maskrecov(bufexpfin.get(), original, bufmaskorigvib.get(), bfh, bfw, ystart, xstart, hig, low, recoth, decay, invmask, sk, multiThread);
                    }

                    if(lp.recothrv >= 1.f) {
                        transit_shapedetect2(sp, 0.f, 0.f, call, 2, bufexporig.get(), bufexpfin.get(), originalmaskvib.get(), hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);
                    } else {
                        transit_shapedetect2(sp, 0.f, 0.f, call, 2, bufexporig.get(), bufexpfin.get(), nullptr, hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);
                        
                    }

                }


                if (lp.recur) {
                    original->CopyFrom(transformed, multiThread);
                    float avge;
                    calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                }
            }
        }
    }


//shadow highlight
    bool tonequ = false;

    if (lp.mullocsh[0] != 0 || lp.mullocsh[1] != 0 || lp.mullocsh[2] != 0 || lp.mullocsh[3] != 0 || lp.mullocsh[4] != 0) {
        tonequ = true;
    }

    bool tonecurv = false;
    const Glib::ustring profile = params->icm.workingProfile;
    bool isworking = (profile == "sRGB" || profile == "Adobe RGB" || profile == "ProPhoto" || profile == "WideGamut" || profile == "BruceRGB" || profile == "Beta RGB" || profile == "BestRGB" || profile == "Rec2020" || profile == "ACESp0" || profile == "ACESp1");

    if (isworking && (params->locallab.spots.at(sp).gamSH != 2.4 || params->locallab.spots.at(sp).sloSH != 12.92)) {
        tonecurv = true;
    }

    if (! lp.invsh && (lp.highlihs > 0.f || lp.shadowhs > 0.f || tonequ || tonecurv || lp.strSH != 0.f || lp.showmaskSHmet == 2 || lp.enaSHMask || lp.showmaskSHmet == 3 || lp.showmaskSHmet == 4 || lp.prevdE) && call <= 3 && lp.hsena) {
        const int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
        const int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
        const int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
        const int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
        const int bfh = yend - ystart;
        const int bfw = xend - xstart;


        if (bfw >= mSP && bfh >= mSP) {

            const std::unique_ptr<LabImage> bufexporig(new LabImage(bfw, bfh));
            const std::unique_ptr<LabImage> bufexpfin(new LabImage(bfw, bfh));
            std::unique_ptr<LabImage> bufmaskorigSH;
            std::unique_ptr<LabImage> bufmaskblurSH;
            std::unique_ptr<LabImage> originalmaskSH;

            if (lp.showmaskSHmet == 2  || lp.enaSHMask || lp.showmaskSHmet == 3 || lp.showmaskSHmet == 4) {
                bufmaskorigSH.reset(new LabImage(bfw, bfh));
                bufmaskblurSH.reset(new LabImage(bfw, bfh));
                originalmaskSH.reset(new LabImage(bfw, bfh));
            }

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = 0; y < bfh; y++) {
                for (int x = 0; x < bfw; x++) {
                    bufexporig->L[y][x] = original->L[y + ystart][x + xstart];
                }
            }

            int inv = 0;
            bool showmaske = false;
            bool enaMask = false;
            bool deltaE = false;
            bool modmask = false;
            bool zero = false;
            bool modif = false;

            if (lp.showmaskSHmet == 3) {
                showmaske = true;
            }

            if (lp.enaSHMask) {
                enaMask = true;
            }

            if (lp.showmaskSHmet == 4) {
                deltaE = true;
            }

            if (lp.showmaskSHmet == 2) {
                modmask = true;
            }

            if (lp.showmaskSHmet == 1) {
                modif = true;
            }

            if (lp.showmaskSHmet == 0) {
                zero = true;
            }

            float chrom = lp.chromaSH;
            float rad = lp.radmaSH;
            float gamma = lp.gammaSH;
            float slope = lp.slomaSH;
            float blendm = lp.blendmaSH;
            float lap = params->locallab.spots.at(sp).lapmaskSH;
            bool pde = params->locallab.spots.at(sp).laplac;
            LocwavCurve dummy;
            int sco = params->locallab.spots.at(sp).scopemask;
            int shortcu = 0;//lp.mergemet; //params->locallab.spots.at(sp).shortc;

            const float mindE = 2.f + MINSCOPE * sco * lp.thr;
            const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
            const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
            const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
            int shado = 0;
            float amountcd = params->locallab.spots.at(sp).fatamountSH;
            float anchorcd = params->locallab.spots.at(sp).fatanchorSH;
            int lumask = params->locallab.spots.at(sp).lumask;
            LocHHmaskCurve lochhhmasCurve;
            const int highl = 0;
            maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, bufexporig.get(), bufmaskorigSH.get(), originalmaskSH.get(), original, reserved, inv, lp,
                        0.f, false,
                        locccmasSHCurve, lcmasSHutili, locllmasSHCurve, llmasSHutili, lochhmasSHCurve, lhmasSHutili, lochhhmasCurve, false, multiThread,
                        enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, highl, amountcd, anchorcd, lmaskSHlocalcurve, localmaskSHutili, dummy, false, 1, 1, 5, 5,
                        shortcu, params->locallab.spots.at(sp).deltae, hueref, chromaref, lumaref,
                        maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, -1, fab
                       );

            if (lp.showmaskSHmet == 3) {
                showmask(lumask, lp, xstart, ystart, cx, cy, bfw, bfh, bufexporig.get(), transformed, bufmaskorigSH.get(), 0);

                return;
            }

            if (lp.showmaskSHmet == 0 || lp.showmaskSHmet == 1  || lp.showmaskSHmet == 2 || lp.showmaskSHmet == 4 || lp.enaSHMask) {

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                for (int y = 0; y < bfh ; y++) {
                    for (int x = 0; x < bfw; x++) {
                        bufexporig->L[y][x] = original->L[y + ystart][x + xstart];
                        bufexporig->a[y][x] = original->a[y + ystart][x + xstart];
                        bufexporig->b[y][x] = original->b[y + ystart][x + xstart];
                        bufexpfin->L[y][x] = original->L[y + ystart][x + xstart];
                        bufexpfin->a[y][x] = original->a[y + ystart][x + xstart];
                        bufexpfin->b[y][x] = original->b[y + ystart][x + xstart];
                    }
                }

                if (lp.shmeth == 0) {
                    ImProcFunctions::shadowsHighlights(bufexpfin.get(), lp.hsena, 1, lp.highlihs, lp.shadowhs, lp.radiushs, sk, lp.hltonalhs, lp.shtonalhs);
                }

//gradient
                struct grad_params gp;

                if (lp.strSH != 0.f) {
                    calclocalGradientParams(lp, gp, ystart, xstart, bfw, bfh, 2);
#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                    for (int ir = 0; ir < bfh; ir++) {
                        for (int jr = 0; jr < bfw; jr++) {
                            bufexpfin->L[ir][jr] *= ImProcFunctions::calcGradientFactor(gp, jr, ir);
                        }
                    }
                }

                if (lp.shmeth == 1) {
                    double scal = (double)(sk);
                    Imagefloat *tmpImage = nullptr;
                    tmpImage = new Imagefloat(bfw, bfh);
                    lab2rgb(*bufexpfin, *tmpImage, params->icm.workingProfile);
                    Glib::ustring prof = params->icm.workingProfile;

                    if (tonecurv) { //Tone response curve  : does nothing if gamma=2.4 and slope=12.92 ==> gamma sRGB
                        float gamtone = params->locallab.spots.at(sp).gamSH;
                        float slotone = params->locallab.spots.at(sp).sloSH;
                        cmsHTRANSFORM dummy = nullptr;
                        int ill =0;
                        workingtrc(tmpImage, tmpImage, bfw, bfh, -5, prof, 2.4, 12.92310, ill, 0, dummy, true, false, false);
                     //   workingtrc(tmpImage, tmpImage, bfw, bfh, 5, prof, gamtone, slotone, 0, 0, dummy, false, true, true); //to keep if we want improve with illuminant and primaries
                        workingtrc(tmpImage, tmpImage, bfw, bfh, 1, prof, gamtone, slotone, ill, 0, dummy, false, true, true);//be careful no gamut control
                    }

                    if (tonequ) {
                        tone_eq(this, tmpImage, lp, params->icm.workingProfile, scal, multiThread);
                    }

                    rgb2lab(*tmpImage, *bufexpfin, params->icm.workingProfile);

                    delete tmpImage;
                }
            }

            if(lp.enaSHMask && lp.recothrs != 1.f) {
                float recoth = lp.recothrs;

                if(lp.recothrs < 1.f) {
                    recoth = -1.f * recoth + 2.f;
                }

                float hig = lp.higthrs;
                float low = lp.lowthrs;
               // float recoth = lp.recothrs;
                float decay = lp.decays;
                bool invmask = false;
                maskrecov(bufexpfin.get(), original, bufmaskorigSH.get(), bfh, bfw, ystart, xstart, hig, low, recoth, decay, invmask, sk, multiThread);
            }

            const float repart = 1.0 - 0.01 * params->locallab.spots.at(sp).reparsh;
            int bw = bufexporig->W;
            int bh = bufexporig->H;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if(multiThread)
#endif
            for (int x = 0; x < bh; x++) {
                for (int y = 0; y < bw; y++) {
                    bufexpfin->L[x][y] = intp(repart, bufexporig->L[x][y], bufexpfin->L[x][y]);
                    bufexpfin->a[x][y] = intp(repart, bufexporig->a[x][y], bufexpfin->a[x][y]);
                    bufexpfin->b[x][y] = intp(repart, bufexporig->b[x][y], bufexpfin->b[x][y]);
                }
            }

            if(lp.recothrs >= 1.f) {
                transit_shapedetect2(sp, 0.f, 0.f, call, 9, bufexporig.get(), bufexpfin.get(), originalmaskSH.get(), hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);
            } else {
                transit_shapedetect2(sp, 0.f, 0.f, call, 9, bufexporig.get(), bufexpfin.get(), nullptr, hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);
            }
            if (lp.recur) {
                original->CopyFrom(transformed, multiThread);
                float avge;
                calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
            }
        }
    } else  if (lp.invsh && (lp.highlihs > 0.f || lp.shadowhs > 0.f || tonequ  || tonecurv || lp.showmaskSHmetinv == 1 || lp.enaSHMaskinv) && call < 3 && lp.hsena) {
        std::unique_ptr<LabImage> bufmaskblurcol;
        std::unique_ptr<LabImage> originalmaskSH;
        const std::unique_ptr<LabImage> bufcolorig(new LabImage(TW, TH));

        if (lp.enaSHMaskinv || lp.showmaskSHmetinv == 1) {
            bufmaskblurcol.reset(new LabImage(TW, TH, true));
            originalmaskSH.reset(new LabImage(TW, TH));
        }

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
        for (int y = 0; y < TH ; y++) {
            for (int x = 0; x < TW; x++) {
                bufcolorig->L[y][x] = original->L[y][x];
            }
        }

        int inv = 1;
        bool showmaske = false;
        bool enaMask = false;
        bool deltaE = false;
        bool modmask = false;
        bool zero = false;
        bool modif = false;

        if (lp.showmaskSHmetinv == 1) {
            showmaske = true;
        }

        if (lp.enaSHMaskinv) {
            enaMask = true;
        }

        if (lp.showmaskSHmetinv == 0) {
            zero = true;
        }

        float chrom = lp.chromaSH;
        float rad = lp.radmaSH;
        float gamma = lp.gammaSH;
        float slope = lp.slomaSH;
        float blendm = lp.blendmaSH;
        float lap = params->locallab.spots.at(sp).lapmaskSH;
        bool pde = params->locallab.spots.at(sp).laplac;
        LocwavCurve dummy;
        int sco = params->locallab.spots.at(sp).scopemask;
        int shortcu = params->locallab.spots.at(sp).shortc;

        const float mindE = 2.f + MINSCOPE * sco * lp.thr;
        const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
        const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
        const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
        int shado = 0;
        float amountcd = params->locallab.spots.at(sp).fatamountSH;
        float anchorcd = params->locallab.spots.at(sp).fatanchorSH;
        int lumask = params->locallab.spots.at(sp).lumask;
        LocHHmaskCurve lochhhmasCurve;
        const int highl = 0;
        maskcalccol(false, pde, TW, TH, 0, 0, sk, cx, cy, bufcolorig.get(), bufmaskblurcol.get(), originalmaskSH.get(), original, reserved, inv, lp,
                    0.f, false,
                    locccmasSHCurve, lcmasSHutili, locllmasSHCurve, llmasSHutili, lochhmasSHCurve, lhmasSHutili, lochhhmasCurve, false, multiThread,
                    enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, highl, amountcd, anchorcd, lmaskSHlocalcurve, localmaskSHutili, dummy, false, 1, 1, 5, 5,
                    shortcu, false, hueref, chromaref, lumaref,
                    maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, -1, fab
                   );


        if (lp.showmaskSHmetinv == 1) {
            showmask(lumask, lp, 0, 0, cx, cy, TW, TH, bufcolorig.get(), transformed, bufmaskblurcol.get(), inv);

            return;
        }

        float adjustr = 2.f;
        InverseColorLight_Local(tonequ, tonecurv, sp, 2, lp, originalmaskSH.get(), lightCurveloc, hltonecurveloc, shtonecurveloc, tonecurveloc, exlocalcurve, cclocalcurve, adjustr, localcutili, lllocalcurve, locallutili, original, transformed, cx, cy, hueref, chromaref, lumaref, sk);

        if (lp.recur) {
            original->CopyFrom(transformed, multiThread);
            float avge;
            calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
        }
    }

// soft light and retinex_pde
    if ((lp.strng > 1.f || lp.prevdE) && call <= 3 && lp.sfena) {
        int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
        int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
        int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
        int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
        int bfh = yend - ystart;
        int bfw = xend - xstart;
        //variable for fast FFTW
        int bfhr = bfh;
        int bfwr = bfw;

        if (bfw >= mSP && bfh >= mSP) {

            if (lp.softmet == 1) {
                optfft(N_fftwsize, bfh, bfw, bfhr, bfwr, lp, original->H, original->W, xstart, ystart, xend, yend, cx, cy, lp.fullim);
            }

            const std::unique_ptr<LabImage> bufexporig(new LabImage(bfw, bfh));
            const std::unique_ptr<LabImage> bufexpfin(new LabImage(bfw, bfh));

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = ystart; y < yend; y++) {
                for (int x = xstart; x < xend; x++) {
                    bufexporig->L[y - ystart][x - xstart] = original->L[y][x];
                    bufexporig->a[y - ystart][x - xstart] = original->a[y][x];
                    bufexporig->b[y - ystart][x - xstart] = original->b[y][x];
                }
            }

            bufexpfin->CopyFrom(bufexporig.get(), multiThread);
            SoftLightParams softLightParams;
            softLightParams.enabled = true;
            softLightParams.strength = lp.strng;

            if (lp.softmet == 0) {
                ImProcFunctions::softLight(bufexpfin.get(), softLightParams);
            } else if (lp.softmet == 1) {
                const std::unique_ptr<float[]> datain(new float[bfwr * bfhr]);
                const std::unique_ptr<float[]> dataout(new float[bfwr * bfhr]);
                const std::unique_ptr<float[]> dE(new float[bfwr * bfhr]);

                deltaEforLaplace(dE.get(), lp.lap, bfwr, bfhr, bufexpfin.get(), hueref, chromaref, lumaref);

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                for (int y = 0; y < bfhr; y++) {
                    for (int x = 0; x < bfwr; x++) {
                        datain[y * bfwr + x] = bufexpfin->L[y][x];
                    }
                }

                const int showorig = lp.showmasksoftmet >= 5 ? 0 : lp.showmasksoftmet;
                MyMutex::MyLock lock(*fftwMutex);
                ImProcFunctions::retinex_pde(datain.get(), dataout.get(), bfwr, bfhr, 8.f * lp.strng, 1.f, dE.get(), showorig, 1, 1);
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                for (int y = 0; y < bfhr; y++) {
                    for (int x = 0; x < bfwr; x++) {
                        bufexpfin->L[y][x] = dataout[y * bfwr + x];
                    }
                }
            }

            transit_shapedetect2(sp, 0.f, 0.f, call, 3, bufexporig.get(), bufexpfin.get(), nullptr, hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);

            if (lp.recur) {
                original->CopyFrom(transformed, multiThread);
                float avge;
                calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
            }
        }
    }

    //local contrast
    bool wavcurve = false;
    bool wavcurvelev = false;
    bool wavcurvecon = false;
    bool wavcurvecomp = false;
    bool wavcurvecompre = false;

    if (lp.locmet == 1) {
        if (locwavCurve && locwavutili) {
            for (int i = 0; i < 500; i++) {
                if (locwavCurve[i] != 0.5f) {
                    wavcurve = true;
                    break;
                }
            }
        }
        if (loclevwavCurve && loclevwavutili) {
            for (int i = 0; i < 500; i++) {
                if (loclevwavCurve[i] != 0.f) {
                    wavcurvelev = true;
                    break;
                }
            }
        }
        if (locconwavCurve && locconwavutili) {
            for (int i = 0; i < 500; i++) {
                if (locconwavCurve[i] != 0.5f) {
                    wavcurvecon = true;
                    break;
                }
            }
        }
        if (loccompwavCurve && loccompwavutili) {
            for (int i = 0; i < 500; i++) {
                if (loccompwavCurve[i] != 0.f) {
                    wavcurvecomp = true;
                    break;
                }
            }
        }
        if (loccomprewavCurve && loccomprewavutili) {
            for (int i = 0; i < 500; i++) {
                if (loccomprewavCurve[i] != 0.75f) {
                    wavcurvecompre = true;
                    break;
                }
            }
        }
    }

    if ((lp.lcamount > 0.f || wavcurve || lp.showmasklcmet == 2 || lp.enalcMask || lp.showmasklcmet == 3 || lp.showmasklcmet == 4 || lp.prevdE || lp.strwav != 0.f || wavcurvelev || wavcurvecon || wavcurvecomp || wavcurvecompre || lp.edgwena || params->locallab.spots.at(sp).residblur > 0.0 || params->locallab.spots.at(sp).levelblur > 0.0 || params->locallab.spots.at(sp).residcont != 0.0 || params->locallab.spots.at(sp).clarilres != 0.0 || params->locallab.spots.at(sp).claricres != 0.0) && call <= 3 && lp.lcena) {

        int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
        int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
        int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
        int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
        int bfh = yend - ystart;
        int bfw = xend - xstart;
        int bfhr = bfh;
        int bfwr = bfw;

        if (bfw >= mSPwav && bfh >= mSPwav) {//avoid too small spot for wavelet
            if (lp.ftwlc) {
                optfft(N_fftwsize, bfh, bfw, bfhr, bfwr, lp, original->H, original->W, xstart, ystart, xend, yend, cx, cy, lp.fullim);
            }

            std::unique_ptr<LabImage> bufmaskblurlc;
            std::unique_ptr<LabImage> originalmasklc;
            std::unique_ptr<LabImage> bufmaskoriglc;

            if (lp.showmasklcmet == 2  || lp.enalcMask || lp.showmasklcmet == 3 || lp.showmasklcmet == 4) {
                bufmaskblurlc.reset(new LabImage(bfw, bfh));
                originalmasklc.reset(new LabImage(bfw, bfh));
                bufmaskoriglc.reset(new LabImage(bfw, bfh));
            }

            array2D<float> buflight(bfw, bfh);
            JaggedArray<float> bufchro(bfw, bfh);
            const std::unique_ptr<LabImage> bufgb(new LabImage(bfw, bfh));
            std::unique_ptr<LabImage> tmp1(new LabImage(bfw, bfh));
            const std::unique_ptr<LabImage> tmpresid(new LabImage(bfw, bfh));
            const std::unique_ptr<LabImage> tmpres(new LabImage(bfw, bfh));

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = ystart; y < yend; y++) {
                for (int x = xstart; x < xend; x++) {
                    bufgb->L[y - ystart][x - xstart] = original->L[y][x];
                    bufgb->a[y - ystart][x - xstart] = original->a[y][x];
                    bufgb->b[y - ystart][x - xstart] = original->b[y][x];
                    tmp1->L[y - ystart][x - xstart] = original->L[y][x];
                    tmp1->a[y - ystart][x - xstart] = original->a[y][x];
                    tmp1->b[y - ystart][x - xstart] = original->b[y][x];
                    tmpresid->L[y - ystart][x - xstart] = original->L[y][x];
                    tmpresid->a[y - ystart][x - xstart] = original->a[y][x];
                    tmpresid->b[y - ystart][x - xstart] = original->b[y][x];
                    tmpres->L[y - ystart][x - xstart] = original->L[y][x];
                    tmpres->a[y - ystart][x - xstart] = original->a[y][x];
                    tmpres->b[y - ystart][x - xstart] = original->b[y][x];
                }
            }

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = 0; y < bfh; y++) {
                for (int x = 0; x < bfw; x++) {
                    bufgb->L[y][x] = original->L[y + ystart][x + xstart];
                }
            }

            int inv = 0;
            bool showmaske = false;
            bool enaMask = false;
            bool deltaE = false;
            bool modmask = false;
            bool zero = false;
            bool modif = false;

            if (lp.showmasklcmet == 3) {
                showmaske = true;
            }

            if (lp.enalcMask) {
                enaMask = true;
            }

            if (lp.showmasklcmet == 4) {
                deltaE = true;
            }

            if (lp.showmasklcmet == 2) {
                modmask = true;
            }

            if (lp.showmasklcmet == 1) {
                modif = true;
            }

            if (lp.showmasklcmet == 0) {
                zero = true;
            }


            float chrom = lp.chromalc;
            float rad = lp.radmalc;
            float blendm = lp.blendmalc;
            float gamma = 1.f;
            float slope = 0.f;
            float lap = 0.f; //params->locallab.spots.at(sp).lapmaskexp;
            bool pde = false; //params->locallab.spots.at(sp).laplac;
            LocwavCurve dummy;
            int sco = params->locallab.spots.at(sp).scopemask;
            int shado = 0;
            int shortcu = 0;//lp.mergemet; //params->locallab.spots.at(sp).shortc;
            const float mindE = 2.f + MINSCOPE * sco * lp.thr;
            const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
            const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
            const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
            float amountcd = 0.f;
            float anchorcd = 50.f;
            int lumask = params->locallab.spots.at(sp).lumask;
            LocHHmaskCurve lochhhmasCurve;
            const int highl = 0;
            maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, bufgb.get(), bufmaskoriglc.get(), originalmasklc.get(), original, reserved, inv, lp,
                        0.f, false,
                        locccmaslcCurve, lcmaslcutili, locllmaslcCurve, llmaslcutili, lochhmaslcCurve, lhmaslcutili, lochhhmasCurve, false, multiThread,
                        enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, highl, amountcd, anchorcd, lmasklclocalcurve, localmasklcutili, dummy, false, 1, 1, 5, 5,
                        shortcu, params->locallab.spots.at(sp).deltae, hueref, chromaref, lumaref,
                        maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, -1, fab
                       );

            if (lp.showmasklcmet == 3) {
                showmask(lumask, lp, xstart, ystart, cx, cy, bfw, bfh, bufgb.get(), transformed, bufmaskoriglc.get(), 0);

                return;
            }

            if (lp.showmasklcmet == 0 || lp.showmasklcmet == 1  || lp.showmasklcmet == 2 || lp.showmasklcmet == 4 || lp.enalcMask) {

                if (lp.locmet == 0) {
                    LocalContrastParams localContrastParams;
                    LocallabParams locallabparams;
                    localContrastParams.enabled = true;
                    localContrastParams.radius = params->locallab.spots.at(sp).lcradius;
                    localContrastParams.amount = params->locallab.spots.at(sp).lcamount;
                    localContrastParams.darkness = params->locallab.spots.at(sp).lcdarkness;
                    localContrastParams.lightness = params->locallab.spots.at(sp).lightness;
                    bool fftwlc = false;

                    if (!lp.ftwlc) { // || (lp.ftwlc && call != 2)) {
                        ImProcFunctions::localContrast(tmp1.get(), tmp1->L, localContrastParams, fftwlc, sk);
                    } else {
                        const std::unique_ptr<LabImage> tmpfftw(new LabImage(bfwr, bfhr));
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < bfhr; y++) {
                            for (int x = 0; x < bfwr; x++) {
                                tmpfftw->L[y][x] = tmp1->L[y][x];
                                tmpfftw->a[y][x] = tmp1->a[y][x];
                                tmpfftw->b[y][x] = tmp1->b[y][x];
                            }
                        }

                        fftwlc = true;
                        ImProcFunctions::localContrast(tmpfftw.get(), tmpfftw->L, localContrastParams, fftwlc, sk);
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < bfhr; y++) {
                            for (int x = 0; x < bfwr; x++) {
                                tmp1->L[y][x] = tmpfftw->L[y][x];
                                tmp1->a[y][x] = tmpfftw->a[y][x];
                                tmp1->b[y][x] = tmpfftw->b[y][x];
                            }
                        }

                    }
                } else if (lp.locmet == 1) { //wavelet && sk ==1
                    int wavelet_level = 1 + params->locallab.spots.at(sp).csthreshold.getBottomRight();//retrieve with +1 maximum wavelet_level
                    float mL = params->locallab.spots.at(sp).clarilres / 100.0;
                    float mC = params->locallab.spots.at(sp).claricres / 100.0;
                    float softr = params->locallab.spots.at(sp).clarisoft;
                    float mL0 = 0.f;
                    float mC0 = 0.f;
#ifdef _OPENMP
                    const int numThreads = omp_get_max_threads();
#else
                    const int numThreads = 1;

#endif
                    // adap maximum level wavelet to size of RT-spot
                    int minwin = rtengine::min(bfw, bfh);
                    int maxlevelspot = 10;//maximum possible

                    // adap maximum level wavelet to size of crop
                    while ((1 << maxlevelspot) >= (minwin * sk) && maxlevelspot  > 1) {
                        --maxlevelspot ;
                    }

                    // printf("minwin=%i maxlevelavant=%i  maxlespot=%i\n", minwin, wavelet_level, maxlevelspot);

                    wavelet_level = rtengine::min(wavelet_level, maxlevelspot);
                    //    printf("maxlevel=%i\n", wavelet_level);
                    bool exec = false;
                    bool origlc = params->locallab.spots.at(sp).origlc;

                    if (origlc) {//merge only with original
                        clarimerge(lp, mL, mC, exec, tmpresid.get(), wavelet_level, sk, numThreads);
                    }

                    int maxlvl = wavelet_level;
                    const float contrast = params->locallab.spots.at(sp).residcont;
                    int level_bl = params->locallab.spots.at(sp).csthreshold.getBottomLeft();
                    int level_hl = params->locallab.spots.at(sp).csthreshold.getTopLeft();
                    int level_br = params->locallab.spots.at(sp).csthreshold.getBottomRight();
                    int level_hr = params->locallab.spots.at(sp).csthreshold.getTopRight();
                    const float radblur = (params->locallab.spots.at(sp).residblur) / sk;
                    const bool blurlc = params->locallab.spots.at(sp).blurlc;
                    const float radlevblur = (params->locallab.spots.at(sp).levelblur) / sk;
                    const float sigma = params->locallab.spots.at(sp).sigma;
                    const float offs = params->locallab.spots.at(sp).offset;
                    const float sigmadc = params->locallab.spots.at(sp).sigmadc;
                    const float deltad = params->locallab.spots.at(sp).deltad;
                   // const float fatres = params->locallab.spots.at(sp).fatres;
                    const float chrol = params->locallab.spots.at(sp).chromalev;
                    const float chrobl = params->locallab.spots.at(sp).chromablu;
                    const bool blurena = params->locallab.spots.at(sp).wavblur;
                    const bool levelena = params->locallab.spots.at(sp).wavcont;
                    const bool comprena = params->locallab.spots.at(sp).wavcomp;
                    const bool compreena = params->locallab.spots.at(sp).wavcompre;
                    const float compress = params->locallab.spots.at(sp).residcomp;
                    const float thres = params->locallab.spots.at(sp).threswav;

                    float gamma = lp.gamlc;
                    rtengine::GammaValues g_a; //gamma parameters
                    double pwr = 1.0 / (double) lp.gamlc;//default 3.0 - gamma Lab
                    double ts = 9.03296;//always the same 'slope' in the extreme shadows - slope Lab
                    rtengine::Color::calcGamma(pwr, ts, g_a); // call to calcGamma with selected gamma and slope

                    if(gamma != 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < tmp1->H; ++y) {
                        int x = 0;
#ifdef __SSE2__
                            for (; x < tmp1->W - 3; x += 4) {
                            STVFU(tmp1->L[y][x], F2V(32768.f) * igammalog(LVFU(tmp1->L[y][x]) / F2V(32768.f), F2V(gamma), F2V(ts), F2V(g_a[2]), F2V(g_a[4])));
                            }
#endif
                            for (;x < tmp1->W; ++x) {
                                tmp1->L[y][x] = 32768.f * igammalog(tmp1->L[y][x] / 32768.f, gamma, ts, g_a[2], g_a[4]);
                            }
                        }
                    }

                    wavcontrast4(lp, tmp1->L, tmp1->a, tmp1->b, contrast, radblur, radlevblur, tmp1->W, tmp1->H, level_bl, level_hl, level_br, level_hr, sk, numThreads, locwavCurve, locwavutili, wavcurve, loclevwavCurve, loclevwavutili, wavcurvelev, locconwavCurve, locconwavutili, wavcurvecon, loccompwavCurve, loccompwavutili, wavcurvecomp, loccomprewavCurve, loccomprewavutili, wavcurvecompre, locedgwavCurve, locedgwavutili, sigma, offs, maxlvl, sigmadc, deltad, chrol, chrobl, blurlc, blurena, levelena, comprena, compreena, compress, thres);

                    if (params->locallab.spots.at(sp).expcie && params->locallab.spots.at(sp).modecie == "wav") {
                        bool HHcurvejz = false, CHcurvejz = false, LHcurvejz = false;
                        if (params->locallab.spots.at(sp).modecam == "jz") {//some cam16 elementsfor Jz
                            ImProcFunctions::ciecamloc_02float(lp, sp, tmp1.get(), bfw, bfh, 10, sk, cielocalcurve, localcieutili, cielocalcurve2, localcieutili2, jzlocalcurve, localjzutili, czlocalcurve, localczutili, czjzlocalcurve, localczjzutili, locchCurvejz, lochhCurvejz, loclhCurvejz, HHcurvejz, CHcurvejz, LHcurvejz, locwavCurvejz, locwavutilijz);
                        }

                        ImProcFunctions::ciecamloc_02float(lp, sp, tmp1.get(), bfw, bfh, 0, sk, cielocalcurve, localcieutili, cielocalcurve2, localcieutili2, jzlocalcurve, localjzutili, czlocalcurve, localczutili, czjzlocalcurve, localczjzutili, locchCurvejz, lochhCurvejz, loclhCurvejz, HHcurvejz, CHcurvejz, LHcurvejz, locwavCurvejz, locwavutilijz);

                        float rad = params->locallab.spots.at(sp).detailcie;
                        loccont(bfw, bfh, tmp1.get(), rad, 5.f, sk);
                    }



                    if(gamma != 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < tmp1->H; ++y) {//apply inverse gamma 3.f and put result in range 32768.f
                            int x = 0;
#ifdef __SSE2__
                            for (; x < tmp1->W - 3; x += 4) {
                                STVFU(tmp1->L[y][x], F2V(32768.f) * gammalog(LVFU(tmp1->L[y][x]) / F2V(32768.f), F2V(gamma), F2V(ts), F2V(g_a[3]), F2V(g_a[4])));
                            }
#endif
                            for (; x < tmp1->W; ++x) {
                                tmp1->L[y][x] = 32768.f * gammalog(tmp1->L[y][x] / 32768.f, gamma, ts, g_a[3], g_a[4]);
                            }
                        }
                    }

                    const float satur = params->locallab.spots.at(sp).residchro;


                    if (satur != 0.f || radblur > 0.f) {//blur residual a and satur

                        wavelet_decomposition *wdspota = new wavelet_decomposition(tmp1->a[0], tmp1->W, tmp1->H, wavelet_level, 1, sk, numThreads, lp.daubLen);

                        if (wdspota->memory_allocation_failed()) {
                            return;
                        }

                        float *wav_ab0a = wdspota->get_coeff0();
                        //      int maxlvla = wdspota->maxlevel();
                        int W_La = wdspota->level_W(0);
                        int H_La = wdspota->level_H(0);

                        if (radblur > 0.f && !blurlc && blurena) {
                            array2D<float> bufa(W_La, H_La);
#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int y = 0; y < H_La; y++) {
                                for (int x = 0; x < W_La; x++) {
                                    bufa[y][x]  = wav_ab0a [y * W_La + x];
                                }
                            }

#ifdef _OPENMP
                            #pragma omp parallel if (multiThread)
#endif
                            {
                                gaussianBlur(bufa, bufa, W_La, H_La, radblur);
                            }

#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int y = 0; y < H_La; y++) {
                                for (int x = 0; x < W_La; x++) {
                                    wav_ab0a[y * W_La + x] = bufa[y][x];
                                }
                            }

                        }

                        if (satur != 0.f) {
#ifdef _OPENMP
                            #pragma omp parallel for if (multiThread)
#endif
                            for (int i = 0; i < W_La * H_La; i++) {
                                wav_ab0a[i] *= (1.f + xsinf(rtengine::RT_PI_F * (satur / 200.f)));//more progressive than linear
                                wav_ab0a[i] = clipC(wav_ab0a[i]);
                            }
                        }

                        wdspota->reconstruct(tmp1->a[0], 1.f);
                        delete wdspota;

                        wavelet_decomposition *wdspotb = new wavelet_decomposition(tmp1->b[0], tmp1->W, tmp1->H, wavelet_level, 1, sk, numThreads, lp.daubLen);

                        if (wdspotb->memory_allocation_failed()) {
                            return;
                        }

                        float *wav_ab0b = wdspotb->get_coeff0();
                        int W_Lb = wdspotb->level_W(0);
                        int H_Lb = wdspotb->level_H(0);

                        if (radblur > 0.f && !blurlc && blurena) {
                            array2D<float> bufb(W_Lb, H_Lb);
#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int y = 0; y < H_Lb; y++) {
                                for (int x = 0; x < W_Lb; x++) {
                                    bufb[y][x]  = wav_ab0b [y * W_Lb + x];
                                }
                            }

#ifdef _OPENMP
                            #pragma omp parallel if (multiThread)
#endif
                            {
                                gaussianBlur(bufb, bufb, W_Lb, H_Lb, radblur);
                            }


#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int y = 0; y < H_Lb; y++) {
                                for (int x = 0; x < W_Lb; x++) {
                                    wav_ab0b[y * W_Lb + x] = bufb[y][x];
                                }
                            }

                        }

                        if (satur != 0.f) {

#ifdef _OPENMP
                            #pragma omp parallel for if (multiThread)
#endif
                            for (int i = 0; i < W_Lb * H_Lb; i++) {
                                wav_ab0b[i] *= (1.f + xsinf(rtengine::RT_PI_F * (satur / 200.f)));
                                wav_ab0b[i] = clipC(wav_ab0b[i]);
                            }
                        }

                        wdspotb->reconstruct(tmp1->b[0], 1.f);
                        delete wdspotb;
                    }

                    if (!origlc) {//merge all files
                        exec = false;
                        //copy previous calculation in merge possibilities
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < bfhr; y++) {
                            for (int x = 0; x < bfwr; x++) {
                                tmpresid->L[y][x] = tmp1->L[y][x];
                                tmpresid->a[y][x] = tmp1->a[y][x];
                                tmpresid->b[y][x] = tmp1->b[y][x];
                            }
                        }
                        clarimerge(lp, mL, mC, exec, tmpresid.get(), wavelet_level, sk, numThreads);
                    }

                    float thr = 0.001f;
                    int flag = 0;

                    if (maxlvl <= 4) {
                        mL0 = 0.f;
                        mC0 = 0.f;
                        mL = -1.5f * mL;//increase only for sharpen
                        mC = -mC;
                        thr = 1.f;
                        flag = 0;

                    } else {
                        mL0 = mL;
                        mC0 = mC;
                        thr = 1.f;
                        flag = 1;
                    }

                    if (exec  || compreena || comprena || levelena || blurena || lp.wavgradl || locwavCurve || lp.edgwena) {
                        LabImage *mergfile = tmp1.get();

#ifdef _OPENMP
                        #pragma omp parallel for if (multiThread)
#endif
                        for (int x = 0; x < bfh; x++)
                            for (int y = 0; y < bfw; y++) {
                                tmp1->L[x][y] = clipLoc((1.f + mL0) * mergfile->L[x][y] - mL * tmpresid->L[x][y]);
                                tmp1->a[x][y] = clipC((1.f + mC0) * mergfile->a[x][y] - mC * tmpresid->a[x][y]);
                                tmp1->b[x][y] = clipC((1.f + mC0) * mergfile->b[x][y] - mC * tmpresid->b[x][y]);
                            }

                        if (softr != 0.f && (compreena || locwavCurve || comprena || blurena || levelena || lp.wavgradl || lp.edgwena || std::fabs(mL) > 0.001f)) {
                            softproc(tmpres.get(), tmp1.get(), softr, bfh, bfw, 0.001, 0.00001, thr, sk, multiThread, flag);
                        }
                    }
                }

                    if(lp.enalcMask && lp.recothrw != 1.f) {
                        float recoth = lp.recothrw;

                        if(lp.recothrw < 1.f) {
                            recoth = -1.f * recoth + 2.f;
                        }

                        float hig = lp.higthrw;
                        float low = lp.lowthrw;
                        //float recoth = lp.recothrw;
                        float decay = lp.decayw;
                        bool invmask = false;
                        maskrecov(tmp1.get(), original, bufmaskoriglc.get(), bfh, bfw, ystart, xstart, hig, low, recoth, decay, invmask, sk, multiThread);
                    }
                const float repart = 1.0 - 0.01 * params->locallab.spots.at(sp).reparw;
                int bw = bufgb->W;
                int bh = bufgb->H;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if(multiThread)
#endif
                    for (int x = 0; x < bh; x++) {
                        for (int y = 0; y < bw; y++) {
                            tmp1->L[x][y] = intp(repart, bufgb->L[x][y], tmp1->L[x][y]);
                            tmp1->a[x][y] = intp(repart, bufgb->a[x][y], tmp1->a[x][y]);
                            tmp1->b[x][y] = intp(repart, bufgb->b[x][y], tmp1->b[x][y]);
                        }
                    }
                if(lp.recothrw >= 1.f) {
                    transit_shapedetect2(sp, 0.f, 0.f, call, 10, bufgb.get(), tmp1.get(), originalmasklc.get(), hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);
                } else {
                    transit_shapedetect2(sp, 0.f, 0.f, call, 10, bufgb.get(), tmp1.get(), nullptr, hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);
                }
                tmp1.reset();
            }

            if (lp.recur) {
                original->CopyFrom(transformed, multiThread);
                float avge;
                calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
            }
        }
    }

    if (!lp.invshar && lp.shrad > 0.42 && call < 3 && lp.sharpena && sk == 1) { //interior ellipse for sharpening, call = 1 and 2 only with Dcrop and simpleprocess
        int bfh = call == 2 ? int (lp.ly + lp.lyT) + del : original->H; //bfw bfh real size of square zone
        int bfw = call == 2 ? int (lp.lx + lp.lxL) + del : original->W;

        if (call == 2) { //call from simpleprocess

            if (bfw < mSPsharp || bfh < mSPsharp) {
                printf("too small RT-spot - minimum size 39 * 39\n");
                return;
            }

            int begy = lp.yc - lp.lyT;
            int begx = lp.xc - lp.lxL;
            int yEn = lp.yc + lp.ly;
            int xEn = lp.xc + lp.lx;

			if(lp.fullim == 2) {//limit sharpening to image dimension...no more...to avoid a long treatment
				begy = 0;
				begx = 0;
				yEn = original->H;
				xEn = original->W;
				lp.lxL = lp.xc;
				lp.lyT = lp.yc;
				lp.ly = yEn - lp.yc;
				lp.lx = xEn - lp.xc;		
				bfh= yEn;
				bfw = xEn;
			}
            //printf("begy=%i begx=%i yen=%i xen=%i\n", begy, begx, yEn, xEn);
			JaggedArray<float> bufsh(bfw, bfh, true);
			JaggedArray<float> hbuffer(bfw, bfh);
			JaggedArray<float> loctemp2(bfw, bfh);

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = 0; y < transformed->H ; y++) {
                for (int x = 0; x < transformed->W; x++) {
                    int lox = cx + x;
                    int loy = cy + y;

                    if (lox >= begx && lox < xEn && loy >= begy && loy < yEn) {
                        bufsh[loy - begy][lox - begx] = original->L[y][x];
                    }
                }
            }
                    float gamma1 = params->locallab.spots.at(sp).shargam;
                    rtengine::GammaValues g_a; //gamma parameters
                    double pwr1 = 1.0 / (double) gamma1;//default 3.0 - gamma Lab
                    double ts1 = 9.03296;//always the same 'slope' in the extreme shadows - slope Lab
                    rtengine::Color::calcGamma(pwr1, ts1, g_a); // call to calcGamma with selected gamma and slope
                    if(gamma1 != 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < bfh; ++y) {
                        int x = 0;
#ifdef __SSE2__
                            for (; x < bfw - 3; x += 4) {
                                STVFU(bufsh[y][x], F2V(32768.f) * igammalog(LVFU(bufsh[y][x]) / F2V(32768.f), F2V(gamma1), F2V(ts1), F2V(g_a[2]), F2V(g_a[4])));
                            }
#endif
                            for (;x < bfw; ++x) {
                                bufsh[y][x] = 32768.f * igammalog(bufsh[y][x] / 32768.f, gamma1, ts1, g_a[2], g_a[4]);
                            }
                        }
                    }


                    //sharpen only square area instead of all image, but limited to image dimensions (full image)
                    ImProcFunctions::deconvsharpeningloc(bufsh, hbuffer, bfw, bfh, loctemp2, params->locallab.spots.at(sp).shardamping, (double)params->locallab.spots.at(sp).sharradius, params->locallab.spots.at(sp).shariter, params->locallab.spots.at(sp).sharamount, params->locallab.spots.at(sp).sharcontrast, (double)params->locallab.spots.at(sp).sharblur, 1);
                    /*
                    float gamma =  params->locallab.spots.at(sp).shargam;
                    double pwr = 1.0 / (double) gamma;//default 3.0 - gamma Lab
                    double ts = 9.03296;//always the same 'slope' in the extreme shadows - slope Lab
                    rtengine::Color::calcGamma(pwr, ts, g_a); // call to calcGamma with selected gamma and slope
                    */
                    if(gamma1 != 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif  
                        for (int y = 0; y < bfh; ++y) {//apply inverse gamma 3.f and put result in range 32768.f
                            int x = 0;
#ifdef __SSE2__
                            for (; x < bfw - 3; x += 4) {
                                STVFU(bufsh[y][x], F2V(32768.f) * gammalog(LVFU(bufsh[y][x]) / F2V(32768.f), F2V(gamma1), F2V(ts1), F2V(g_a[3]), F2V(g_a[4])));
                                STVFU(loctemp2[y][x], F2V(32768.f) * gammalog(LVFU(loctemp2[y][x]) / F2V(32768.f), F2V(gamma1), F2V(ts1), F2V(g_a[3]), F2V(g_a[4])));
                            }
#endif
                            for (; x < bfw; ++x) {
                                bufsh[y][x] = 32768.f * gammalog(bufsh[y][x] / 32768.f, gamma1, ts1, g_a[3], g_a[4]);
                                loctemp2[y][x] = 32768.f * gammalog(loctemp2[y][x] / 32768.f, gamma1, ts1, g_a[3], g_a[4]);
                            }
                        }
                    }
			//sharpen simpleprocess
			Sharp_Local(call, loctemp2, 0, hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);
        } else { //call from dcrop.cc
					JaggedArray<float> loctemp(bfw, bfh);
		
                    float gamma1 = params->locallab.spots.at(sp).shargam;
                    rtengine::GammaValues g_a; //gamma parameters
                    double pwr1 = 1.0 / (double) gamma1;//default 3.0 - gamma Lab
                    double ts1 = 9.03296;//always the same 'slope' in the extreme shadows - slope Lab
                    rtengine::Color::calcGamma(pwr1, ts1, g_a); // call to calcGamma with selected gamma and slope
                    if(gamma1 != 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < bfh; ++y) {
                        int x = 0;
#ifdef __SSE2__
                            for (; x < bfw - 3; x += 4) {
                                STVFU(original->L[y][x], F2V(32768.f) * igammalog(LVFU(original->L[y][x]) / F2V(32768.f), F2V(gamma1), F2V(ts1), F2V(g_a[2]), F2V(g_a[4])));
                            }
#endif
                            for (;x < bfw; ++x) {
                                original->L[y][x] = 32768.f * igammalog(original->L[y][x] / 32768.f, gamma1, ts1, g_a[2], g_a[4]);
                            }
                        }
                    }
        
        
                    ImProcFunctions::deconvsharpeningloc(original->L, shbuffer, bfw, bfh, loctemp, params->locallab.spots.at(sp).shardamping, (double)params->locallab.spots.at(sp).sharradius, params->locallab.spots.at(sp).shariter, params->locallab.spots.at(sp).sharamount, params->locallab.spots.at(sp).sharcontrast, (double)params->locallab.spots.at(sp).sharblur, sk);
                    /*
                    float gamma =  params->locallab.spots.at(sp).shargam;
                    double pwr = 1.0 / (double) gamma;//default 3.0 - gamma Lab
                    double ts = 9.03296;//always the same 'slope' in the extreme shadows - slope Lab
                    rtengine::Color::calcGamma(pwr, ts, g_a); // call to calcGamma with selected gamma and slope
                    */
                    if(gamma1 != 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif  
                        for (int y = 0; y < bfh; ++y) {//apply inverse gamma 3.f and put result in range 32768.f
                            int x = 0;
#ifdef __SSE2__
                            for (; x < bfw - 3; x += 4) {
                                STVFU(original->L[y][x], F2V(32768.f) * gammalog(LVFU(original->L[y][x]) / F2V(32768.f), F2V(gamma1), F2V(ts1), F2V(g_a[3]), F2V(g_a[4])));
                                STVFU(loctemp[y][x], F2V(32768.f) * gammalog(LVFU(loctemp[y][x]) / F2V(32768.f), F2V(gamma1), F2V(ts1), F2V(g_a[3]), F2V(g_a[4])));
                            }
#endif
                            for (; x < bfw; ++x) {
                                original->L[y][x] = 32768.f * gammalog(original->L[y][x] / 32768.f, gamma1, ts1, g_a[3], g_a[4]);
                                loctemp[y][x] = 32768.f * gammalog(loctemp[y][x] / 32768.f, gamma1, ts1, g_a[3], g_a[4]);
                            }
                        }
                    }
			//sharpen dcrop
			Sharp_Local(call, loctemp, 0, hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);
        }

		
        if (lp.recur) {
            original->CopyFrom(transformed, multiThread);
            float avge;
            calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
        }

    } else if (lp.invshar && lp.shrad > 0.42 && call < 3 && lp.sharpena && sk == 1) {
        int GW = original->W;
        int GH = original->H;
        JaggedArray<float> loctemp(GW, GH);

        float gamma1 = params->locallab.spots.at(sp).shargam;
        rtengine::GammaValues g_a; //gamma parameters
        double pwr1 = 1.0 / (double) gamma1;//default 3.0 - gamma Lab
        double ts1 = 9.03296;//always the same 'slope' in the extreme shadows - slope Lab
        rtengine::Color::calcGamma(pwr1, ts1, g_a); // call to calcGamma with selected gamma and slope
        if(gamma1 != 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = 0; y < GH; ++y) {
                int x = 0;
#ifdef __SSE2__
                for (; x < GW - 3; x += 4) {
                    STVFU(original->L[y][x], F2V(32768.f) * igammalog(LVFU(original->L[y][x]) / F2V(32768.f), F2V(gamma1), F2V(ts1), F2V(g_a[2]), F2V(g_a[4])));
                }
#endif
                for (;x < GW; ++x) {
                    original->L[y][x] = 32768.f * igammalog(original->L[y][x] / 32768.f, gamma1, ts1, g_a[2], g_a[4]);
                }
            }
        }


        ImProcFunctions::deconvsharpeningloc(original->L, shbuffer, GW, GH, loctemp, params->locallab.spots.at(sp).shardamping, (double)params->locallab.spots.at(sp).sharradius, params->locallab.spots.at(sp).shariter, params->locallab.spots.at(sp).sharamount, params->locallab.spots.at(sp).sharcontrast, (double)params->locallab.spots.at(sp).sharblur, sk);
        /*
        float gamma =  params->locallab.spots.at(sp).shargam;
        double pwr = 1.0 / (double) gamma;//default 3.0 - gamma Lab
        double ts = 9.03296;//always the same 'slope' in the extreme shadows - slope Lab
        rtengine::Color::calcGamma(pwr, ts, g_a); // call to calcGamma with selected gamma and slope
        */
        if(gamma1 != 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif  
            for (int y = 0; y < GH; ++y) {//apply inverse gamma 3.f and put result in range 32768.f
                int x = 0;
#ifdef __SSE2__
                for (; x < GW - 3; x += 4) {
                    STVFU(original->L[y][x], F2V(32768.f) * gammalog(LVFU(original->L[y][x]) / F2V(32768.f), F2V(gamma1), F2V(ts1), F2V(g_a[3]), F2V(g_a[4])));
                    STVFU(loctemp[y][x], F2V(32768.f) * igammalog(LVFU(loctemp[y][x]) / F2V(32768.f), F2V(gamma1), F2V(ts1), F2V(g_a[3]), F2V(g_a[4])));
                }
#endif
                for (; x < GW; ++x) {
                    original->L[y][x] = 32768.f * gammalog(original->L[y][x] / 32768.f, gamma1, ts1, g_a[3], g_a[4]);
                    loctemp[y][x] = 32768.f * igammalog(loctemp[y][x] / 32768.f, gamma1, ts1, g_a[3], g_a[4]);
                }
            }
        }


        InverseSharp_Local(loctemp, hueref, lumaref, chromaref, lp, original, transformed, cx, cy, sk);

        if (lp.recur) {
            original->CopyFrom(transformed, multiThread);
            float avge;
            calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
        }
    }

    bool enablefat = false;

    if (params->locallab.spots.at(sp).fatamount > 1.0) {
        enablefat = true;;
    }

    bool execex = (lp.exposena && (lp.expcomp != 0.f || lp.blac != 0 || lp.shadex > 0 || lp.hlcomp > 0.f || lp.laplacexp > 0.1f || lp.strexp != 0.f || enablefat || lp.showmaskexpmet == 2 || lp.enaExpMask || lp.showmaskexpmet == 3 || lp.showmaskexpmet == 4  || lp.showmaskexpmet == 5 || lp.prevdE || (exlocalcurve && localexutili)));

    if (!lp.invex && execex) {
        int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
        int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
        int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
        int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
        int bfh = yend - ystart;
        int bfw = xend - xstart;
        //variable for fast FFTW
        int bfhr = bfh;
        int bfwr = bfw;


        if (bfw >= mSP && bfh >= mSP) {

            if (lp.expmet == 1  || lp.expmet == 0) {
                optfft(N_fftwsize, bfh, bfw, bfhr, bfwr, lp, original->H, original->W, xstart, ystart, xend, yend, cx, cy, lp.fullim);
            }

            const std::unique_ptr<LabImage> bufexporig(new LabImage(bfw, bfh));
            const std::unique_ptr<LabImage> bufexpfin(new LabImage(bfw, bfh));
            const std::unique_ptr<LabImage> buforig(new LabImage(bfw, bfh));

            std::unique_ptr<LabImage> bufmaskblurexp;
            std::unique_ptr<LabImage> originalmaskexp;

            array2D<float> blend2;

            if (call <= 3) { //simpleprocess, dcrop, improccoordinator
                if (lp.showmaskexpmet == 2  || lp.enaExpMask || lp.showmaskexpmet == 3 || lp.showmaskexpmet == 5) {
                    bufmaskblurexp.reset(new LabImage(bfw, bfh));
                    originalmaskexp.reset(new LabImage(bfw, bfh));
                }

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                for (int y = 0; y < bfh ; y++) {
                    for (int x = 0; x < bfw; x++) {
                        bufexporig->L[y][x] = original->L[y + ystart][x + xstart];
                        buforig->a[y][x] = original->a[y + ystart][x + xstart];
                    }
                }

                float gamma1 = lp.gamex;
                rtengine::GammaValues g_a; //gamma parameters
                double pwr1 = 1.0 / (double) lp.gamex;//default 3.0 - gamma Lab
                double ts1 = 9.03296;//always the same 'slope' in the extreme shadows - slope Lab
                rtengine::Color::calcGamma(pwr1, ts1, g_a); // call to calcGamma with selected gamma and slope

                if(gamma1 != 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < bfh; ++y) {
                        int x = 0;
#ifdef __SSE2__
                            for (; x < bfw - 3; x += 4) {
                            STVFU(bufexporig->L[y][x], F2V(32768.f) * igammalog(LVFU(bufexporig->L[y][x]) / F2V(32768.f), F2V(gamma1), F2V(ts1), F2V(g_a[2]), F2V(g_a[4])));
                            }
#endif
                            for (;x < bfw; ++x) {
                                bufexporig->L[y][x] = 32768.f * igammalog(bufexporig->L[y][x] / 32768.f, gamma1, ts1, g_a[2], g_a[4]);
                            }
                        }
                    }

                const int spotSi = rtengine::max(1 + 2 * rtengine::max(1, lp.cir / sk), 5);

                if (bfw > 2 * spotSi && bfh > 2 * spotSi && lp.struexp > 0.f) {
                    blend2(bfw, bfh);
                    ImProcFunctions::blendstruc(bfw, bfh, bufexporig.get(), 3.f / (sk * 1.4f), 0.5f * lp.struexp, blend2, sk, multiThread);

                    if (lp.showmaskexpmet == 4) {
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = ystart; y < yend ; y++) {
                            for (int x = xstart; x < xend; x++) {
                                const int lox = cx + x;
                                const int loy = cy + y;
                                int zone;
                                float localFactor = 1.f;
                                const float achm = lp.trans / 100.f;

                                if (lp.shapmet == 0) {
                                    calcTransition(lox, loy, achm, lp, zone, localFactor);
                                } else /*if (lp.shapmet == 1)*/ {
                                    calcTransitionrect(lox, loy, achm, lp, zone, localFactor);
                                }

                                if (zone > 0) {
                                    transformed->L[y][x] = CLIP(blend2[y - ystart][x - xstart]);
                                    transformed->a[y][x] = 0.f;
                                    transformed->b[y][x] = 0.f;
                                }
                            }
                        }

                        return;
                    }
                }

                int inv = 0;
                bool showmaske = false;
                const bool enaMask = lp.enaExpMask;
                bool deltaE = false;
                bool modmask = false;
                bool zero = false;
                bool modif = false;

                if (lp.showmaskexpmet == 3) {
                    showmaske = true;
                } else if (lp.showmaskexpmet == 5) {
                    deltaE = true;
                } else if (lp.showmaskexpmet == 2) {
                    modmask = true;
                } else if (lp.showmaskexpmet == 1) {
                    modif = true;
                } else if (lp.showmaskexpmet == 0) {
                    zero = true;
                }

                float chrom = lp.chromaexp;
                float rad = lp.radmaexp;
                float gamma = lp.gammaexp;
                float slope = lp.slomaexp;
                float blendm = lp.blendmaexp;
                float lap = params->locallab.spots.at(sp).lapmaskexp;
                bool pde = params->locallab.spots.at(sp).laplac;
                LocwavCurve dummy;
                int sco = params->locallab.spots.at(sp).scopemask;
                int shado = 0;
                int shortcu = 0;//lp.mergemet; //params->locallab.spots.at(sp).shortc;

                const float mindE = 2.f + MINSCOPE * sco * lp.thr;
                const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
                const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
                const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
                float amountcd = 0.f;
                float anchorcd = 50.f;
                int lumask = params->locallab.spots.at(sp).lumask;
                LocHHmaskCurve lochhhmasCurve;
                const int highl = 0;
                maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, bufexporig.get(), bufmaskblurexp.get(), originalmaskexp.get(), original, reserved, inv, lp,
                            0.f, false,
                            locccmasexpCurve, lcmasexputili, locllmasexpCurve, llmasexputili, lochhmasexpCurve, lhmasexputili, lochhhmasCurve, false, multiThread,
                            enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, highl, amountcd, anchorcd, lmaskexplocalcurve, localmaskexputili, dummy, false, 1, 1, 5, 5,
                            shortcu, params->locallab.spots.at(sp).deltae, hueref, chromaref, lumaref,
                            maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, 0, fab
                           );

                if (lp.showmaskexpmet == 3) {
                    showmask(lumask, lp, xstart, ystart, cx, cy, bfw, bfh, bufexporig.get(), transformed, bufmaskblurexp.get(), 0);

                    return;
                }

                if (lp.showmaskexpmet == 4) {
                    return;
                }

                if (lp.showmaskexpmet == 0 || lp.showmaskexpmet == 1  || lp.showmaskexpmet == 2  || lp.showmaskexpmet == 5 || lp.enaExpMask) {
#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                    for (int y = 0; y < bfh; y++) {
                        for (int x = 0; x < bfw; x++) {
                            bufexpfin->L[y][x] = original->L[y + ystart][x + xstart];
                            bufexpfin->a[y][x] = original->a[y + ystart][x + xstart];
                            bufexpfin->b[y][x] = original->b[y + ystart][x + xstart];
                        }
                    }



                    if (exlocalcurve && localexutili) {// L=f(L) curve enhanced

#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int ir = 0; ir < bfh; ir++)
                            for (int jr = 0; jr < bfw; jr++) {
                                bufexpfin->L[ir][jr] = 0.6f * bufexporig->L[ir][jr] + 0.2f * exlocalcurve[2.f * bufexporig->L[ir][jr]];
                            }
                        
                        if (lp.expcomp == 0.f) {
                            lp.expcomp = 0.001f;// to enabled
                        }

                        ImProcFunctions::exlabLocal(lp, 0.5f, bfh, bfw, bfhr, bfwr, bufexpfin.get(), bufexpfin.get(), hltonecurveloc, shtonecurveloc, tonecurveloc, hueref, lumaref, chromaref);



                    } else {
                        if (lp.expcomp == 0.f  && (lp.linear > 0.01f && lp.laplacexp > 0.1f)) {
                            lp.expcomp = 0.001f;// to enabled
                        } 

                        if (lp.expcomp != 0.f  ) { // ||  lp.laplacexp > 0.1f
                            if(lp.laplacexp <= 0.1f) {
                                lp.laplacexp = 0.2f;  //force to use Laplacian with very small values
                            }
                            ImProcFunctions::exlabLocal(lp, 1.f, bfh, bfw, bfhr, bfwr, bufexporig.get(), bufexpfin.get(), hltonecurveloc, shtonecurveloc, tonecurveloc, hueref, lumaref, chromaref);
                        }
                    }

//gradient
                    struct grad_params gp;

                    if (lp.strexp != 0.f) {

                        calclocalGradientParams(lp, gp, ystart, xstart, bfw, bfh, 1);
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int ir = 0; ir < bfh; ir++) {
                            for (int jr = 0; jr < bfw; jr++) {
                                bufexpfin->L[ir][jr] *= ImProcFunctions::calcGradientFactor(gp, jr, ir);
                            }
                        }
                    }

//exposure_pde
                    if (lp.expmet == 1) {
                        if (enablefat) {

                            const std::unique_ptr<float[]> datain(new float[bfwr * bfhr]);
                            const std::unique_ptr<float[]> dataout(new float[bfwr * bfhr]);
#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int y = 0; y < bfhr; y++) {
                                for (int x = 0; x < bfwr; x++) {
                                    datain[y * bfwr + x] = bufexpfin->L[y][x];
                                }
                            }
                            FattalToneMappingParams fatParams;
                            fatParams.enabled = true;
                            fatParams.threshold = params->locallab.spots.at(sp).fatdetail;
                            fatParams.amount = params->locallab.spots.at(sp).fatamount;
                            fatParams.anchor = params->locallab.spots.at(sp).fatanchor;
                            //const float sigm = 1.f; //params->locallab.spots.at(sp).fatlevel;
                            //const float mean = 1.f;// params->locallab.spots.at(sp).fatanchor;
                            const std::unique_ptr<Imagefloat> tmpImagefat(new Imagefloat(bfwr, bfhr));
                            lab2rgb(*bufexpfin, *tmpImagefat, params->icm.workingProfile);
                            int alg = 0;
                            if(fatParams.anchor == 50.f) {
                                alg = 1;
                            }
                            ToneMapFattal02(tmpImagefat.get(), fatParams, 3, 0, nullptr, 0, 0, alg);//last parameter = 1 ==>ART algorithm
                            rgb2lab(*tmpImagefat, *bufexpfin, params->icm.workingProfile);
                            if (params->locallab.spots.at(sp).expcie && params->locallab.spots.at(sp).modecie == "dr") {
                                bool HHcurvejz = false, CHcurvejz = false, LHcurvejz = false;
                                if (params->locallab.spots.at(sp).modecam == "jz") {//some cam16 elementsfor Jz
                                    ImProcFunctions::ciecamloc_02float(lp, sp, bufexpfin.get(), bfw, bfh, 10, sk, cielocalcurve, localcieutili, cielocalcurve2, localcieutili2, jzlocalcurve, localjzutili, czlocalcurve, localczutili, czjzlocalcurve, localczjzutili, locchCurvejz, lochhCurvejz, loclhCurvejz, HHcurvejz, CHcurvejz, LHcurvejz, locwavCurvejz, locwavutilijz);
                                }

                                ImProcFunctions::ciecamloc_02float(lp, sp, bufexpfin.get(), bfw, bfh, 0, sk, cielocalcurve, localcieutili, cielocalcurve2, localcieutili2, jzlocalcurve, localjzutili, czlocalcurve, localczutili, czjzlocalcurve, localczjzutili, locchCurvejz, lochhCurvejz, loclhCurvejz, HHcurvejz, CHcurvejz, LHcurvejz, locwavCurvejz, locwavutilijz);

                                float rad = params->locallab.spots.at(sp).detailcie;
                                loccont(bfw, bfh, bufexpfin.get(), rad, 15.f, sk);
                            }

                        }

                        if (lp.laplacexp > 0.1f) {

                            MyMutex::MyLock lock(*fftwMutex);
                            std::unique_ptr<float[]> datain(new float[bfwr * bfhr]);
                            std::unique_ptr<float[]> dataout(new float[bfwr * bfhr]);
                            const float gam = params->locallab.spots.at(sp).gamm;
                            const float igam = 1.f / gam;

                            if (params->locallab.spots.at(sp).exnoiseMethod == "med" || params->locallab.spots.at(sp).exnoiseMethod == "medhi") {
                                if (lp.blac < -100.f && lp.linear > 0.01f) {
                                    float evnoise = lp.blac - lp.linear * 2000.f;
                                    if (params->locallab.spots.at(sp).exnoiseMethod == "med") {
                                        evnoise *= 0.4f;
                                    }

                                    //soft denoise, user must use Local Denoise for best result
                                    Median med;
                                    if (evnoise < -18000.f) {
                                        med = Median::TYPE_5X5_STRONG;
                                    } else if (evnoise < -15000.f) {
                                        med = Median::TYPE_5X5_SOFT;
                                    } else if (evnoise < -10000.f) {
                                        med = Median::TYPE_3X3_STRONG;
                                    } else {
                                        med = Median:: TYPE_3X3_SOFT;
                                    }

                                    Median_Denoise(bufexpfin->L, bufexpfin->L, bfwr, bfhr, med, 1, multiThread);
                                    Median_Denoise(bufexpfin->a, bufexpfin->a, bfwr, bfhr, Median::TYPE_3X3_SOFT, 1, multiThread);
                                    Median_Denoise(bufexpfin->b, bufexpfin->b, bfwr, bfhr, Median::TYPE_3X3_SOFT, 1, multiThread);
                                }
                            }

#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int y = 0; y < bfhr; y++) {
                                for (int x = 0; x < bfwr; x++) {
                                    float L = LIM01(bufexpfin->L[y][x] / 32768.f);//change gamma for Laplacian
                                    datain[y * bfwr + x] = pow_F(L, gam) * 32768.f;
                                }
                            }

                            //call PDE equation - with Laplacian threshold
                            ImProcFunctions::exposure_pde(datain.get(), datain.get(), dataout.get(), bfwr, bfhr, 12.f * lp.laplacexp, lp.balanexp);
#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int y = 0; y < bfhr; y++) {
                                for (int x = 0; x < bfwr; x++) {
                                    const float Y = dataout[y * bfwr + x] / 32768.f;//inverse Laplacian gamma
                                    bufexpfin->L[y][x] = pow_F(Y, igam) * 32768.f;
                                }
                            }
                        }
                    }
                    if (lp.shadex > 0) {

                        if (lp.expcomp == 0.f) {
                            lp.expcomp = 0.001f;    // to enabled
                        }
                    }

                    if (lp.hlcomp > 0.f) {

                        if (lp.expcomp == 0.f) {
                            lp.expcomp = 0.001f;    // to enabled
                        }
                    }
                    
                    //shadows with ipshadowshighlight
                    if ((lp.expcomp != 0.f) || (exlocalcurve && localexutili)) {
                        if (lp.shadex > 0) {

                            ImProcFunctions::shadowsHighlights(bufexpfin.get(), true, 1, 0, lp.shadex, 40, sk, 0, lp.shcomp);
                        }
                    }

                    if (lp.expchroma != 0.f) {
                        if ((lp.expcomp != 0.f && lp.expcomp != 0.001f) || (exlocalcurve && localexutili) || lp.laplacexp > 0.1f) {

                            constexpr float ampli = 70.f;
                            const float ch = (1.f + 0.02f * lp.expchroma);
                            const float chprosl = ch <= 1.f ? 99.f * ch - 99.f : clipChro(ampli * ch - ampli);
#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int ir = 0; ir < bfh; ir++) {
                                for (int jr = 0; jr < bfw; jr++) {
                                    const float epsi = bufexporig->L[ir][jr] == 0.f ? 0.001f : 0.f;
                                    const float rapexp = bufexpfin->L[ir][jr] / (bufexporig->L[ir][jr] + epsi);
                                    bufexpfin->a[ir][jr] *= 1.f + chprosl * rapexp;
                                    bufexpfin->b[ir][jr] *= 1.f + chprosl * rapexp;
                                }
                            }
                        }
                    }
                    /*
                    float gamma = lp.gamex;
                    rtengine::GammaValues g_a; //gamma parameters
                    double pwr = 1.0 / (double) lp.gamex;//default 3.0 - gamma Lab
                    double ts = 9.03296;//always the same 'slope' in the extreme shadows - slope Lab
                    rtengine::Color::calcGamma(pwr, ts, g_a); // call to calcGamma with selected gamma and slope
                    */
                    if(gamma1 != 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif  
                        for (int y = 0; y < bfh; ++y) {//apply inverse gamma 3.f and put result in range 32768.f
                            int x = 0;
#ifdef __SSE2__
                            for (; x < bfw - 3; x += 4) {
                                STVFU(bufexpfin->L[y][x], F2V(32768.f) * gammalog(LVFU(bufexpfin->L[y][x]) / F2V(32768.f), F2V(gamma1), F2V(ts1), F2V(g_a[3]), F2V(g_a[4])));
                            }
#endif
                            for (; x < bfw; ++x) {
                                bufexpfin->L[y][x] = 32768.f * gammalog(bufexpfin->L[y][x] / 32768.f, gamma1, ts1, g_a[3], g_a[4]);
                            }
                        }
                    }

                    if (lp.softradiusexp > 0.f && lp.expmet == 0) {
                        softproc(buforig.get(), bufexpfin.get(), lp.softradiusexp, bfh, bfw, 0.1, 0.001, 0.5f, sk, multiThread, 1);
                    }
                    
                    if(lp.enaExpMask && lp.recothre != 1.f) {
                        float recoth = lp.recothre;

                        if(lp.recothre < 1.f) {
                            recoth = -1.f * recoth + 2.f;
                        }

                        float hig = lp.higthre;
                        float low = lp.lowthre;
                      //  float recoth = lp.recothre;
                        float decay = lp.decaye;
                        bool invmask = false;
                        maskrecov(bufexpfin.get(), original, bufmaskblurexp.get(), bfh, bfw, ystart, xstart, hig, low, recoth, decay, invmask, sk, multiThread);
                    }
                    
                    float meansob = 0.f;

                    const float repart = 1.0 - 0.01 * params->locallab.spots.at(sp).reparexp;
                    int bw = bufexporig->W;
                    int bh = bufexporig->H;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if(multiThread)
#endif
                    for (int x = 0; x < bh; x++) {
                        for (int y = 0; y < bw; y++) {
                            bufexpfin->L[x][y] = intp(repart, bufexporig->L[x][y], bufexpfin->L[x][y]);
                            bufexpfin->a[x][y] = intp(repart, bufexporig->a[x][y], bufexpfin->a[x][y]);
                            bufexpfin->b[x][y] = intp(repart, bufexporig->b[x][y], bufexpfin->b[x][y]);
                        }
                    }

                    if(lp.recothre >= 1.f) {
                        transit_shapedetect2(sp, 0.f, 0.f, call, 1, bufexporig.get(), bufexpfin.get(), originalmaskexp.get(), hueref, chromaref, lumaref, sobelref, meansob, blend2, lp, original, transformed, cx, cy, sk);
                    } else {
                        transit_shapedetect2(sp, 0.f, 0.f, call, 1, bufexporig.get(), bufexpfin.get(), nullptr, hueref, chromaref, lumaref, sobelref, meansob, blend2, lp, original, transformed, cx, cy, sk);
                    }
                }

                if (lp.recur) {
                    original->CopyFrom(transformed, multiThread);
                    float avge;
                    calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                }
            }
        }
    }
//inverse

    else if (lp.invex && (lp.expcomp != 0.f  || lp.laplacexp > 0.1f || lp.blac != 0 || lp.hlcomp > 0.f || lp.shadex > 0 || params->locallab.spots.at(sp).fatamount > 1.0 || (exlocalcurve && localexutili) || lp.enaExpMaskinv || lp.showmaskexpmetinv == 1) && lp.exposena) {
        constexpr float adjustr = 2.f;
        std::unique_ptr<LabImage> bufmaskblurexp;
        std::unique_ptr<LabImage> originalmaskexp;
        const std::unique_ptr<LabImage> bufexporig(new LabImage(TW, TH));

        if (lp.enaExpMaskinv || lp.showmaskexpmetinv == 1) {
            bufmaskblurexp.reset(new LabImage(TW, TH, true));
            originalmaskexp.reset(new LabImage(TW, TH));
        }

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
        for (int y = 0; y < TH ; y++) {
            for (int x = 0; x < TW; x++) {
                bufexporig->L[y][x] = original->L[y][x];
            }
        }

        constexpr int inv = 1;
        const bool showmaske = lp.showmaskexpmetinv == 1;
        const bool enaMask = lp.enaExpMaskinv;
        constexpr bool deltaE = false;
        constexpr bool modmask = false;
        const bool zero = lp.showmaskexpmetinv == 0;
        constexpr bool modif = false;
        const float chrom = lp.chromaexp;
        const float rad = lp.radmaexp;
        const float gamma = lp.gammaexp;
        const float slope = lp.slomaexp;
        const float blendm = lp.blendmaexp;
        const float lap = params->locallab.spots.at(sp).lapmaskexp;
        const bool pde = params->locallab.spots.at(sp).laplac;
        LocwavCurve dummy;
        const int sco = params->locallab.spots.at(sp).scopemask;
        constexpr int shado = 0;
        constexpr int shortcu = 0;//lp.mergemet; //params->locallab.spots.at(sp).shortc;
        const int lumask = params->locallab.spots.at(sp).lumask;

        const float mindE = 2.f + MINSCOPE * sco * lp.thr;
        const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
        const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
        const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
        constexpr float amountcd = 0.f;
        constexpr float anchorcd = 50.f;
        LocHHmaskCurve lochhhmasCurve;
        const int highl = 0;
        maskcalccol(false, pde, TW, TH, 0, 0, sk, cx, cy, bufexporig.get(), bufmaskblurexp.get(), originalmaskexp.get(), original, reserved, inv, lp,
                    0.f, false,
                    locccmasexpCurve, lcmasexputili, locllmasexpCurve, llmasexputili, lochhmasexpCurve, lhmasexputili, lochhhmasCurve, false, multiThread,
                    enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, highl, amountcd, anchorcd, lmaskexplocalcurve, localmaskexputili, dummy, false, 1, 1, 5, 5,
                    shortcu, false, hueref, chromaref, lumaref,
                    maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, 0, fab
                   );

        if (lp.showmaskexpmetinv == 1) {
            showmask(lumask, lp, 0, 0, cx, cy, TW, TH, bufexporig.get(), transformed, bufmaskblurexp.get(), inv);
            return;
        }

        if (lp.shadex > 0) {
            if (lp.expcomp == 0.f) {
                lp.expcomp = 0.001f;    // to enabled
            }
        }

        if (lp.hlcomp > 0.f) {
            if (lp.expcomp == 0.f) {
                lp.expcomp = 0.001f;    // to enabled
            }
        }

        InverseColorLight_Local(false, false, sp, 1, lp, originalmaskexp.get(), lightCurveloc, hltonecurveloc, shtonecurveloc, tonecurveloc, exlocalcurve, cclocalcurve, adjustr, localcutili, lllocalcurve, locallutili, original, transformed, cx, cy, hueref, chromaref, lumaref, sk);

        if (lp.recur) {
            original->CopyFrom(transformed, multiThread);
            float avge;
            calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
        }
    }

//local color and light
    const float factor = LocallabParams::LABGRIDL_CORR_MAX * 3.276;
    const float scaling = LocallabParams::LABGRIDL_CORR_SCALE;
    const float scaledirect = LocallabParams::LABGRIDL_DIRECT_SCALE;
    const float a_scale = (lp.highA - lp.lowA) / factor / scaling;
    const float a_base = lp.lowA / scaling;
    const float b_scale = (lp.highB - lp.lowB) / factor / scaling;
    const float b_base = lp.lowB / scaling;
    const bool ctoning = (a_scale != 0.f || b_scale != 0.f || a_base != 0.f || b_base != 0.f);
    const float a_scalemerg = (lp.highAmerg - lp.lowAmerg) / factor / scaling;
    const float b_scalemerg = (lp.highBmerg - lp.lowBmerg) / factor / scaling;

    if (!lp.inv && (lp.chro != 0 || lp.ligh != 0.f || lp.cont != 0 || ctoning || lp.mergemet > 0 ||  lp.strcol != 0.f ||  lp.strcolab != 0.f || lp.qualcurvemet != 0 || lp.showmaskcolmet == 2 || lp.enaColorMask || lp.showmaskcolmet == 3  || lp.showmaskcolmet == 4 || lp.showmaskcolmet == 5 || lp.prevdE) && lp.colorena) { // || lllocalcurve)) { //interior ellipse reinforced lightness and chroma  //locallutili
        int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
        int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
        int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
        int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
        int bfh = yend - ystart;
        int bfw = xend - xstart;
        const bool spez = params->locallab.spots.at(sp).special;

        if (bfw >= mSP && bfh >= mSP) {

            if (lp.blurcolmask >= 0.25f && lp.fftColorMask && call == 2) {
                optfft(N_fftwsize, bfh, bfw, bfh, bfw, lp, original->H, original->W, xstart, ystart, xend, yend, cx, cy, lp.fullim);
            }

            std::unique_ptr<LabImage> bufcolorig;
            std::unique_ptr<LabImage> bufcolfin;
            std::unique_ptr<LabImage> bufmaskblurcol;
            std::unique_ptr<LabImage> originalmaskcol;
            std::unique_ptr<LabImage> bufcolreserv;
            std::unique_ptr<LabImage> buftemp;
            array2D<float> blend2;

            float adjustr = 1.0f;

            //adapt chroma to working profile
            if (params->icm.workingProfile == "ProPhoto")   {
                adjustr = 1.2f;   // 1.2 instead 1.0 because it's very rare to have C>170..
            } else if (params->icm.workingProfile == "Adobe RGB")  {
                adjustr = 1.8f;
            } else if (params->icm.workingProfile == "sRGB")       {
                adjustr = 2.0f;
            } else if (params->icm.workingProfile == "WideGamut")  {
                adjustr = 1.2f;
            } else if (params->icm.workingProfile == "Beta RGB")   {
                adjustr = 1.4f;
            } else if (params->icm.workingProfile == "BestRGB")    {
                adjustr = 1.4f;
            } else if (params->icm.workingProfile == "BruceRGB")   {
                adjustr = 1.8f;
            }

            if (call <= 3) { //simpleprocess, dcrop, improccoordinator
                bufcolorig.reset(new LabImage(bfw, bfh));
                bufcolfin.reset(new LabImage(bfw, bfh));
                buftemp.reset(new LabImage(bfw, bfh));

                if (lp.showmaskcolmet == 2  || lp.enaColorMask || lp.showmaskcolmet == 3 || lp.showmaskcolmet == 5) {
                    bufmaskblurcol.reset(new LabImage(bfw, bfh, true));
                    originalmaskcol.reset(new LabImage(bfw, bfh));
                }

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                for (int y = 0; y < bfh ; y++) {
                    for (int x = 0; x < bfw; x++) {
                        bufcolorig->L[y][x] = original->L[y + ystart][x + xstart];
                        bufcolorig->a[y][x] = original->a[y + ystart][x + xstart];
                        bufcolorig->b[y][x] = original->b[y + ystart][x + xstart];
                        bufcolfin->L[y][x] = original->L[y + ystart][x + xstart];
                        bufcolfin->a[y][x] = original->a[y + ystart][x + xstart];
                        bufcolfin->b[y][x] = original->b[y + ystart][x + xstart];
                        buftemp->L[y][x] = original->L[y + ystart][x + xstart];
                        buftemp->a[y][x] = original->a[y + ystart][x + xstart];
                        buftemp->b[y][x] = original->b[y + ystart][x + xstart];
                    }
                }

                    float gamma1 = lp.gamc;
                    rtengine::GammaValues g_a; //gamma parameters
                    double pwr1 = 1.0 / (double) lp.gamc;//default 3.0 - gamma Lab
                    double ts1 = 9.03296;//always the same 'slope' in the extreme shadows - slope Lab
                    rtengine::Color::calcGamma(pwr1, ts1, g_a); // call to calcGamma with selected gamma and slope

                    if(gamma1 != 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < bufcolorig->H; ++y) {
                        int x = 0;
#ifdef __SSE2__
                            for (; x < bufcolorig->W - 3; x += 4) {
                            STVFU(bufcolorig->L[y][x], F2V(32768.f) * igammalog(LVFU(bufcolorig->L[y][x]) / F2V(32768.f), F2V(gamma1), F2V(ts1), F2V(g_a[2]), F2V(g_a[4])));
                            }
#endif
                            for (;x < bufcolorig->W; ++x) {
                                bufcolorig->L[y][x] = 32768.f * igammalog(bufcolorig->L[y][x] / 32768.f, gamma1, ts1, g_a[2], g_a[4]);
                            }
                        }
                    }

                const int spotSi = rtengine::max(1 + 2 * rtengine::max(1,  lp.cir / sk), 5);
                const bool blends = bfw > 2 * spotSi && bfh > 2 * spotSi && lp.struco > 0.f;

                if (blends) {
                    blend2(bfw, bfh);
                    ImProcFunctions::blendstruc(bfw, bfh, bufcolorig.get(), 3.f / (sk * 1.4f), 0.5f * lp.struco, blend2, sk, multiThread);

                    if (lp.showmaskcolmet == 4) {
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = ystart; y < yend ; y++) {
                            for (int x = xstart; x < xend; x++) {
                                const int lox = cx + x;
                                const int loy = cy + y;
                                int zone;
                                float localFactor = 1.f;
                                const float achm = lp.trans / 100.f;

                                if (lp.shapmet == 0) {
                                    calcTransition(lox, loy, achm, lp, zone, localFactor);
                                } else /*if (lp.shapmet == 1)*/ {
                                    calcTransitionrect(lox, loy, achm, lp, zone, localFactor);
                                }

                                if (zone > 0) {
                                    transformed->L[y][x] = CLIP(blend2[y - ystart][x - xstart]);
                                    transformed->a[y][x] = 0.f;
                                    transformed->b[y][x] = 0.f;
                                }
                            }
                        }
                        return;
                    }
                }

                const int inv = 0;
                const bool showmaske = lp.showmaskcolmet == 3;
                const bool enaMask = lp.enaColorMask;
                const bool deltaE = lp.showmaskcolmet == 5;
                const bool modmask = lp.showmaskcolmet == 2;
                const bool zero = lp.showmaskcolmet == 0;
                const bool modif = lp.showmaskcolmet == 1;
                const float chrom = lp.chromacol;
                const float rad = lp.radmacol;
                const float gamma = lp.gammacol;
                const float slope = lp.slomacol;
                const float blendm = lp.blendmacol;
                const float lap = params->locallab.spots.at(sp).lapmaskcol;
                const bool pde = params->locallab.spots.at(sp).laplac;
                const int shado = params->locallab.spots.at(sp).shadmaskcol;
                const int sco = params->locallab.spots.at(sp).scopemask;
                const int level_bl = params->locallab.spots.at(sp).csthresholdcol.getBottomLeft();
                const int level_hl = params->locallab.spots.at(sp).csthresholdcol.getTopLeft();
                const int level_br = params->locallab.spots.at(sp).csthresholdcol.getBottomRight();
                const int level_hr = params->locallab.spots.at(sp).csthresholdcol.getTopRight();
                const int shortcu = lp.mergemet; //params->locallab.spots.at(sp).shortc;
                const int lumask = params->locallab.spots.at(sp).lumask;
                const float strumask = 0.02 * params->locallab.spots.at(sp).strumaskcol;
                float conthr = 0.01 * params->locallab.spots.at(sp).conthrcol;
                const float mercol = params->locallab.spots.at(sp).mercol;
                const float merlucol = params->locallab.spots.at(sp).merlucol;

                int tonemod = 0;
                if (params->locallab.spots.at(sp).toneMethod == "one") {
                    tonemod = 0;
                } else if (params->locallab.spots.at(sp).toneMethod == "two") {
                    tonemod = 1;
                } else if (params->locallab.spots.at(sp).toneMethod == "thr") {
                    tonemod = 2;
                } else if (params->locallab.spots.at(sp).toneMethod == "fou") {
                    tonemod = 3;
                }

                const float mindE = 2.f + MINSCOPE * sco * lp.thr;
                const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
                const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
                const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
                const float amountcd = 0.f;
                const float anchorcd = 50.f;
                const int highl = 0;
                bool astool = params->locallab.spots.at(sp).toolcol;
                maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, bufcolorig.get(), bufmaskblurcol.get(), originalmaskcol.get(), original, reserved, inv, lp,
                            strumask, astool,
                            locccmasCurve, lcmasutili, locllmasCurve, llmasutili, lochhmasCurve, lhmasutili, llochhhmasCurve, lhhmasutili, multiThread,
                            enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, highl, amountcd, anchorcd, lmasklocalcurve, localmaskutili, loclmasCurvecolwav, lmasutilicolwav,
                            level_bl, level_hl, level_br, level_hr,
                            shortcu, delt, hueref, chromaref, lumaref,
                            maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, lp.fftColorMask, lp.blurcolmask, lp.contcolmask, -1, fab
                           );

                if (lp.showmaskcolmet == 3) {
                    showmask(lumask, lp, xstart, ystart, cx, cy, bfw, bfh, bufcolorig.get(), transformed, bufmaskblurcol.get(), 0);
                    return;
                } else if (lp.showmaskcolmet == 4) {
                    return;
                }

                if (lp.showmaskcolmet == 0 || lp.showmaskcolmet == 1 || lp.showmaskcolmet == 2 || lp.showmaskcolmet == 5 || lp.enaColorMask) {
                    //RGB Curves
                    bool usergb = false;

                    if (rgblocalcurve && localrgbutili && lp.qualcurvemet != 0) {
                        usergb = true;
                        const std::unique_ptr<Imagefloat> tmpImage(new Imagefloat(bfw, bfh));

                        lab2rgb(*buftemp, *tmpImage, params->icm.workingProfile);
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < bfh; y++)
                            for (int x = 0; x < bfw; x++) {

                                //std
                                if (tonemod == 0) {
                                    curves::setLutVal(rgblocalcurve, tmpImage->r(y, x), tmpImage->g(y, x), tmpImage->b(y, x));
                                } else {
                                    float r = CLIP(tmpImage->r(y, x));
                                    float g = CLIP(tmpImage->g(y, x));
                                    float b = CLIP(tmpImage->b(y, x));

                                    if (tonemod == 1) { // weightstd
                                        const float r1 = rgblocalcurve[r];
                                        const float g1 = triangle(r, r1, g);
                                        const float b1 = triangle(r, r1, b);

                                        const float g2 = rgblocalcurve[g];
                                        const float r2 = triangle(g, g2, r);
                                        const float b2 = triangle(g, g2, b);

                                        const float b3 = rgblocalcurve[b];
                                        const float r3 = triangle(b, b3, r);
                                        const float g3 = triangle(b, b3, g);
                                        r = CLIP(r1 * 0.50f + r2 * 0.25f + r3 * 0.25f);
                                        g = CLIP(g1 * 0.25f + g2 * 0.50f + g3 * 0.25f);
                                        b = CLIP(b1 * 0.25f + b2 * 0.25f + b3 * 0.50f);
                                    } else if (tonemod == 2) { // Luminance
                                        float currLuminance = r * 0.2126729f + g * 0.7151521f + b * 0.0721750f;

                                        const float newLuminance = rgblocalcurve[currLuminance];
                                        currLuminance = currLuminance == 0.f ? 0.00001f : currLuminance;
                                        const float coef = newLuminance / currLuminance;
                                        r = LIM(r * coef, 0.f, 65535.f);
                                        g = LIM(g * coef, 0.f, 65535.f);
                                        b = LIM(b * coef, 0.f, 65535.f);
                                    } else if (tonemod == 3) { // Film like Adobe
                                        if (r >= g) {
                                            if (g > b) {
                                                rgbtone(r, g, b, rgblocalcurve);     // Case 1: r >= g >  b
                                            } else if (b > r) {
                                                rgbtone(b, r, g, rgblocalcurve);     // Case 2: b >  r >= g
                                            } else if (b > g) {
                                                rgbtone(r, b, g, rgblocalcurve);     // Case 3: r >= b >  g
                                            } else {                           // Case 4: r == g == b
                                                r = rgblocalcurve[r];
                                                g = rgblocalcurve[g];
                                                b = g;
                                            }
                                        } else {
                                            if (r >= b) {
                                                rgbtone(g, r, b, rgblocalcurve);     // Case 5: g >  r >= b
                                            } else if (b >  g) {
                                                rgbtone(b, g, r, rgblocalcurve);     // Case 6: b >  g >  r
                                            } else {
                                                rgbtone(g, b, r, rgblocalcurve);     // Case 7: g >= b >  r
                                            }
                                        }
                                    }

                                    setUnlessOOG(tmpImage->r(y, x), tmpImage->g(y, x), tmpImage->b(y, x), r, g, b);
                                }
                            }

                        rgb2lab(*tmpImage, *buftemp, params->icm.workingProfile);

                        // end rgb curves
                    }

                    if (usergb && spez) {//special use of rgb curves ex : negative
                        const float achm = lp.trans / 100.f;
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < bfh; y++) {
                            const int loy = y + ystart + cy;

                            for (int x = 0; x < bfw; x++) {
                                const int lox = x + xstart + cx;
                                int zone;
                                float localFactor = 1.f;

                                if (lp.shapmet == 0) {
                                    calcTransition(lox, loy, achm, lp, zone, localFactor);
                                } else /*if (lp.shapmet == 1)*/ {
                                    calcTransitionrect(lox, loy, achm, lp, zone, localFactor);
                                }

                                if (zone > 0) {
                                    transformed->L[y + ystart][x + xstart] = buftemp->L[y][x] * localFactor + (1.f - localFactor) * original->L[y + ystart][x + xstart];
                                    transformed->a[y + ystart][x + xstart] = buftemp->a[y][x] * localFactor + (1.f - localFactor) * original->a[y + ystart][x + xstart];
                                    transformed->b[y + ystart][x + xstart] = buftemp->b[y][x] * localFactor + (1.f - localFactor) * original->b[y + ystart][x + xstart];
                                }
                            }
                        }
                    }

                    //others curves

                    const LabImage *origptr = usergb ? buftemp.get() : bufcolorig.get();

                    bool execcolor = false;

                    if (localcutili || HHutili || locallutili || lp.ligh != 0.f || lp.cont != 0 || lp.chro != 0 || LHutili || ctoning) {
                        execcolor = true;
                    }

                    bool HHcurve = false;
                    if (lochhCurve && HHutili) {
                        for (int i = 0; i < 500; i++) {
                            if (lochhCurve[i] != 0.5f) {
                                HHcurve = true;
                                break;
                            }
                        }
                    }

                    const float kd = 10.f * 0.01f * lp.strengrid;//correction to ctoning

                    //chroma slider with curve instead of linear
                    const float satreal = lp.chro;

                    DiagonalCurve color_satur({
                        DCT_NURBS,
                        0, 0,
                        0.2, 0.2f + satreal / 250.f,
                        0.6, rtengine::min(1.f, 0.6f + satreal / 250.f),
                        1, 1
                    });

                    DiagonalCurve color_saturmoins({
                        DCT_NURBS,
                        0, 0,
                        0.1f - satreal / 150.f, 0.1f,
                        rtengine::min(1.f, 0.7f - satreal / 300.f), 0.7,
                        1, 1
                    });
                    bool LHcurve = false;
                    if (loclhCurve && LHutili) {
                        for (int i = 0; i < 500; i++) {
                            if (loclhCurve[i] != 0.5f) {
                                LHcurve = true;
                                break;
                            }
                        }
                    }
                    bool CHcurve = false;
                    if (locchCurve && CHutili) {
                        for (int i = 0; i < 500; i++) {
                            if (locchCurve[i] != 0.5f) {
                                CHcurve = true;
                                break;
                            }
                        }
                    }
                    double amountchrom = 0.01 * settings->amchroma;
                    if(amountchrom < 0.05) {
                        amountchrom = 0.05;
                    }
                    if(amountchrom > 2.) {
                        amountchrom = 2.;
                    }

#ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                    for (int ir = 0; ir < bfh; ir++) {
                        for (int jr = 0; jr < bfw; jr++) {
                            float bufcolcalca = origptr->a[ir][jr];
                            float bufcolcalcb = origptr->b[ir][jr];
                            float bufcolcalcL = origptr->L[ir][jr];

                            if (lp.chro != 0.f) {//slider chroma with curve DCT_NURBS
                                float Chprov = std::sqrt(SQR(bufcolcalca) + SQR(bufcolcalcb));
                                float2 sincosval;
                                sincosval.y = Chprov == 0.0f ? 1.f : bufcolcalca / Chprov;
                                sincosval.x = Chprov == 0.0f ? 0.f : bufcolcalcb / Chprov;

                                // 35000 must be globally good, more than 32768...and less than !! to avoid calculation min max
                                if (lp.chro > 0.f) {
                                    Chprov = static_cast<float>(color_satur.getVal(LIM01(Chprov / 35000.f))) * 35000.f;
                                } else {
                                    Chprov = static_cast<float>(color_saturmoins.getVal(LIM01(Chprov / 35000.f))) * 35000.f;
                                }

                                if (lp.chro == -100.f) {
                                    Chprov = 0.f;
                                }

                                bufcolcalca = Chprov * sincosval.y;
                                bufcolcalcb = Chprov * sincosval.x;
                            }

                            if (cclocalcurve && lp.qualcurvemet != 0 && localcutili) { // C=f(C) curve
                                const float chromat = std::sqrt(SQR(bufcolcalca) + SQR(bufcolcalcb));
                                const float ch = cclocalcurve[chromat * adjustr] / ((chromat + 0.00001f) * adjustr); //ch between 0 and 0 50 or more
                                bufcolcalca *= ch;
                                bufcolcalcb *= ch;
                            }

                            if (cllocalcurve && lp.qualcurvemet != 0 && localclutili) { // C=f(L) curve
                                float chromaCfactor = (cllocalcurve[bufcolcalcL * 2.f]) / (bufcolcalcL * 2.f);
                                bufcolcalca *= chromaCfactor;
                                bufcolcalcb *= chromaCfactor;
                            }

                            if (lclocalcurve && lp.qualcurvemet != 0 && locallcutili) { // L=f(C) curve
                                const float chromat = std::sqrt(SQR(bufcolcalca) + SQR(bufcolcalcb));
                                float Lc = lclocalcurve[chromat * adjustr] / ((chromat + 0.00001f) * adjustr);

                                if (Lc > 1.f) {
                                    Lc = (Lc - 1.0f) * 0.1f + 1.0f; //reduct action
                                } else {
                                    Lc = (Lc - 1.0f) * 0.3f + 1.0f;
                                }

                                bufcolcalcL *= Lc;
                            }

                            if (lochhCurve && HHcurve && lp.qualcurvemet != 0 && !ctoning) { // H=f(H)
                                const float chromat = std::sqrt(SQR(bufcolcalca) + SQR(bufcolcalcb));
                                const float hhforcurv = xatan2f(bufcolcalcb, bufcolcalca);
                                const float valparam = 2.f * (lochhCurve[500.f * static_cast<float>(Color::huelab_to_huehsv2(hhforcurv))] - 0.5f) + hhforcurv;
                                float2 sincosval = xsincosf(valparam);
                                bufcolcalca = chromat * sincosval.y;
                                bufcolcalcb = chromat * sincosval.x;
                            }

                            if (lp.ligh != 0.f || lp.cont != 0) {//slider luminance or slider contrast with curve
                                bufcolcalcL = calclight(bufcolcalcL, lightCurveloc);
                            }

                            if (lllocalcurve && locallutili && lp.qualcurvemet != 0) {// L=f(L) curve
                                bufcolcalcL = 0.5f * lllocalcurve[bufcolcalcL * 2.f];
                            }


                            if (loclhCurve && LHcurve && lp.qualcurvemet != 0) {//L=f(H) curve
                                const float rhue = xatan2f(bufcolcalcb, bufcolcalca);
                                //printf("rhu=%f", (double) rhue);
                                const float chromat = (std::sqrt(SQR(bufcolcalca) + SQR(bufcolcalcb)))/32768.f;
                                float l_r = LIM01(bufcolcalcL / 32768.f); //Luminance Lab in 0..1
                                float valparam = loclhCurve[500.f *static_cast<float>(Color::huelab_to_huehsv2(rhue))] - 0.5f;  //get l_r=f(H)
                               // printf("rh=%f V=%f", (double) rhue, (double) valparam);
                               // float kc = 0.05f + 0.02f * params->locallab.spots.at(sp).lightjzcie;
                                float kc = amountchrom;
                                float valparamneg;
                                valparamneg = valparam;
                                float kcc = SQR(chromat / kc); //take Chroma into account...40 "middle low" of chromaticity (arbitrary and simple), one can imagine other algorithme
                             //   printf("KC=%f", (double) kcc);
                                //reduce action for low chroma and increase action for high chroma
                                valparam *= 2.f * kcc;
                                valparamneg *= kcc; //slightly different for negative

                                if (valparam > 0.f) {
                                    l_r = (1.f - valparam) * l_r + valparam * (1.f - SQR(((SQR(1.f - min(l_r, 1.0f))))));
                                } else
                                //for negative
                                {
                                    float khue = 1.9f; //in reserve in case of!
                                    l_r *= (1.f + khue * valparamneg);
                                }

                                bufcolcalcL = l_r * 32768.f;

                            }
                            if (locchCurve && CHcurve && lp.qualcurvemet != 0) {//C=f(H) curve
                                const float rhue = xatan2f(bufcolcalcb, bufcolcalca);
                                const float valparam = 2.f * locchCurve[500.f * static_cast<float>(Color::huelab_to_huehsv2(rhue))] - 0.5f;  //get valp=f(H)
                                float chromaChfactor = 1.0f + valparam;
                                bufcolcalca *= chromaChfactor;//apply C=f(H)
                                bufcolcalcb *= chromaChfactor;
                            }

                            if (ctoning) {//color toning and direct change color
                                if (lp.gridmet == 0) {
                                    bufcolcalca += kd * (bufcolcalcL * a_scale + a_base);
                                    bufcolcalcb += kd * (bufcolcalcL * b_scale + b_base);
                                } else if (lp.gridmet == 1) {
                                    bufcolcalca += kd * scaledirect * a_scale;
                                    bufcolcalcb += kd * scaledirect * b_scale;
                                }

                                bufcolcalca = clipC(bufcolcalca);
                                bufcolcalcb = clipC(bufcolcalcb);

                            }

                            bufcolfin->L[ir][jr] = bufcolcalcL;
                            bufcolfin->a[ir][jr] = bufcolcalca;
                            bufcolfin->b[ir][jr] = bufcolcalcb;
                        }
                    }

                    if (!execcolor) {//if we don't use color and light sliders, curves except RGB
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int ir = 0; ir < bfh; ir++)
                            for (int jr = 0; jr < bfw; jr++) {
                                bufcolfin->L[ir][jr] = origptr->L[ir][jr];
                                bufcolfin->a[ir][jr] = origptr->a[ir][jr];
                                bufcolfin->b[ir][jr] = origptr->b[ir][jr];
                            }
                    }

                    bool nottransit = false;
                    if (lp.mergemet >= 2) { //merge result with original
                        nottransit = true;
                        bufcolreserv.reset(new LabImage(bfw, bfh));
                        JaggedArray<float> lumreserv(bfw, bfh);
                        const std::unique_ptr<LabImage> bufreser(new LabImage(bfw, bfh));

#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < bfh ; y++) {
                            for (int x = 0; x < bfw; x++) {
                                lumreserv[y][x] = 32768.f - reserved->L[y + ystart][x + xstart];
                                bufreser->L[y][x] = reserved->L[y + ystart][x + xstart];
                                bufreser->a[y][x] = reserved->a[y + ystart][x + xstart];
                                bufreser->b[y][x] = reserved->b[y + ystart][x + xstart];

                                if (lp.mergemet == 2) {
                                    bufcolreserv->L[y][x] = reserved->L[y + ystart][x + xstart];
                                    bufcolreserv->a[y][x] = reserved->a[y + ystart][x + xstart];
                                    bufcolreserv->b[y][x] = reserved->b[y + ystart][x + xstart];
                                } else if (lp.mergemet == 3) {
                                    bufcolreserv->L[y][x] = lastorig->L[y + ystart][x + xstart];
                                    bufcolreserv->a[y][x] = lastorig->a[y + ystart][x + xstart];
                                    bufcolreserv->b[y][x] = lastorig->b[y + ystart][x + xstart];
                                } else if (lp.mergemet == 4) {
                                        bufcolreserv->L[y][x] = merlucol * 327.68f;
                                        bufcolreserv->a[y][x] = 9.f * scaledirect * a_scalemerg;
                                        bufcolreserv->b[y][x] = 9.f * scaledirect * b_scalemerg;
                                }
                            }
                        }

                        if (lp.strcol != 0.f) {
                            struct grad_params gp;
                            calclocalGradientParams(lp, gp, ystart, xstart, bfw, bfh, 3);
#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int ir = 0; ir < bfh; ir++) {
                                for (int jr = 0; jr < bfw; jr++) {
                                    const float corrFactor = ImProcFunctions::calcGradientFactor(gp, jr, ir);
                                    bufcolfin->L[ir][jr] *= corrFactor;
                                }
                            }
                        }

                        if (lp.strcolab != 0.f) {
                            struct grad_params gpab;
                            calclocalGradientParams(lp, gpab, ystart, xstart, bfw, bfh, 4);
#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int ir = 0; ir < bfh; ir++) {
                                for (int jr = 0; jr < bfw; jr++) {
                                    const float corrFactor = ImProcFunctions::calcGradientFactor(gpab, jr, ir);
                                    bufcolfin->a[ir][jr] *= corrFactor;
                                    bufcolfin->b[ir][jr] *= corrFactor;
                                }
                            }
                        }

                        if (lp.strcolh != 0.f) {
                            struct grad_params gph;
                            calclocalGradientParams(lp, gph, ystart, xstart, bfw, bfh, 6);
#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int ir = 0; ir < bfh; ir++) {
                                for (int jr = 0; jr < bfw; jr++) {
                                    const float corrFactor = ImProcFunctions::calcGradientFactor(gph, jr, ir);
                                    const float aa = bufcolfin->a[ir][jr];
                                    const float bb = bufcolfin->b[ir][jr];
                                    const float chrm = std::sqrt(SQR(aa) + SQR(bb));
                                    const float HH = xatan2f(bb, aa);

                                    float cor = 0.f;
                                    if (corrFactor < 1.f) {
                                        cor = - 2.5f * (1.f - corrFactor);
                                    } else if (corrFactor > 1.f) {
                                        cor = 0.03f * (corrFactor - 1.f);
                                    }

                                    float newhr = HH + cor;
                                    if (newhr > rtengine::RT_PI_F) {
                                        newhr -= 2 * rtengine::RT_PI_F;
                                    } else if (newhr < -rtengine::RT_PI_F) {
                                        newhr += 2 * rtengine::RT_PI_F;
                                    }

                                    const float2 sincosval = xsincosf(newhr);
                                    bufcolfin->a[ir][jr] = clipC(chrm * sincosval.y);
                                    bufcolfin->b[ir][jr] = clipC(chrm * sincosval.x);
                                }
                            }
                        }

                        JaggedArray<float> blend(bfw, bfh);
                        buildBlendMask(lumreserv, blend, bfw, bfh, conthr);
                        const float rm = 20.f / sk;

                        if (rm > 0) {
                            float **mb = blend;
#ifdef _OPENMP
                            #pragma omp parallel if (multiThread)
#endif
                            {
                                gaussianBlur(mb, mb, bfw, bfh, rm);
                            }
                        }

                        const std::unique_ptr<JaggedArray<float>> rdEBuffer(new JaggedArray<float>(bfw, bfh));
                        float** rdE = *rdEBuffer;

                        deltaEforMask(rdE, bfw, bfh, bufreser.get(), hueref, chromaref, lumaref, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, mercol, lp.balance, lp.balanceh);

                        if (lp.mergecolMethod == 0) {  //normal

                            if (lp.mergemet == 4) {
                                bufprov.reset(new LabImage(bfw, bfh));

#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        rdE[y][x] *= SQR(rdE[y][x]);
                                        bufprov->L[y][x] = intp(rdE[y][x], bufcolreserv->L[y][x], bufcolfin->L[y][x]);
                                        bufprov->a[y][x] = intp(rdE[y][x], bufcolreserv->a[y][x], bufcolfin->a[y][x]);
                                        bufprov->b[y][x] = intp(rdE[y][x], bufcolreserv->b[y][x], bufcolfin->b[y][x]);

                                        bufcolfin->L[y][x] = intp(lp.opacol, bufprov->L[y][x], bufcolfin->L[y][x]);
                                        bufcolfin->a[y][x] = intp(lp.opacol, bufprov->a[y][x], bufcolfin->a[y][x]);
                                        bufcolfin->b[y][x] = intp(lp.opacol, bufprov->b[y][x], bufcolfin->b[y][x]);
                                    }
                                }
                            } else {
                                bufprov.reset(new LabImage(bfw, bfh));

#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        bufprov->L[y][x] = intp(rdE[y][x], bufcolfin->L[y][x], bufcolreserv->L[y][x]);
                                        bufprov->a[y][x] = intp(rdE[y][x], bufcolfin->a[y][x], bufcolreserv->a[y][x]);
                                        bufprov->b[y][x] = intp(rdE[y][x], bufcolfin->b[y][x], bufcolreserv->b[y][x]);

                                        bufcolfin->L[y][x] = intp(lp.opacol, bufprov->L[y][x], bufcolreserv->L[y][x]);
                                        bufcolfin->a[y][x] = intp(lp.opacol, bufprov->a[y][x], bufcolreserv->a[y][x]);
                                        bufcolfin->b[y][x] = intp(lp.opacol, bufprov->b[y][x], bufcolreserv->b[y][x]);
                                    }
                                }
                            }

                            if (conthr > 0.f && lp.mergemet != 4) {
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        bufcolfin->L[y][x] = intp(blend[y][x], bufcolfin->L[y][x], bufreser->L[y][x]);
                                        bufcolfin->a[y][x] = intp(blend[y][x], bufcolfin->a[y][x], bufreser->a[y][x]);
                                        bufcolfin->b[y][x] = intp(blend[y][x], bufcolfin->b[y][x], bufreser->b[y][x]);
                                    }
                                }
                            }
                        }

                        if (lp.mergecolMethod > 16) { //hue sat chroma luma
                            bufprov.reset(new LabImage(bfw, bfh));

                            if (lp.mergemet == 4) {
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        rdE[y][x] *= SQR(rdE[y][x]);
                                        bufprov->L[y][x] = intp(rdE[y][x], bufcolreserv->L[y][x], bufcolfin->L[y][x]);
                                        bufprov->a[y][x] = intp(rdE[y][x], bufcolreserv->a[y][x], bufcolfin->a[y][x]);
                                        bufprov->b[y][x] = intp(rdE[y][x], bufcolreserv->b[y][x], bufcolfin->b[y][x]);
                                    }
                                }
                            } else {
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        bufprov->L[y][x] = intp(rdE[y][x], bufcolfin->L[y][x], bufcolreserv->L[y][x]);
                                        bufprov->a[y][x] = intp(rdE[y][x], bufcolfin->a[y][x], bufcolreserv->a[y][x]);
                                        bufprov->b[y][x] = intp(rdE[y][x], bufcolfin->b[y][x], bufcolreserv->b[y][x]);
                                    }
                                }
                            }

#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int y = 0; y < bfh ; y++) {
                                for (int x = 0; x < bfw; x++) {
                                    if (lp.mergecolMethod == 17) {
                                        const float huefin = xatan2f(bufprov->b[y][x], bufprov->a[y][x]);
                                        const float2 sincosval1 = xsincosf(huefin);
                                        const float chrores = std::sqrt(SQR(bufcolreserv->a[y][x]) + SQR(bufcolreserv->b[y][x]));
                                        buftemp->a[y][x] = chrores * sincosval1.y;
                                        buftemp->b[y][x] = chrores * sincosval1.x;
                                        buftemp->L[y][x] = bufcolreserv->L[y][x];
                                    } else if (lp.mergecolMethod == 18) {
                                        const float hueres = xatan2f(bufcolreserv->b[y][x], bufcolreserv->a[y][x]);
                                        const float2 sincosval2 = xsincosf(hueres);
                                        const float chrofin = std::sqrt(SQR(bufprov->a[y][x]) + SQR(bufprov->b[y][x]));
                                        buftemp->a[y][x] = chrofin * sincosval2.y;
                                        buftemp->b[y][x] = chrofin * sincosval2.x;
                                        buftemp->L[y][x] = bufcolreserv->L[y][x];
                                    } else if (lp.mergecolMethod == 19) {
                                        const float huefin = xatan2f(bufprov->b[y][x], bufprov->a[y][x]);
                                        const float2 sincosval3 = xsincosf(huefin);
                                        const float chrofin = std::sqrt(SQR(bufprov->a[y][x]) + SQR(bufprov->b[y][x]));
                                        buftemp->a[y][x] = chrofin * sincosval3.y;
                                        buftemp->b[y][x] = chrofin * sincosval3.x;
                                        buftemp->L[y][x] = bufcolreserv->L[y][x];
                                    } else if (lp.mergecolMethod == 20) {
                                        const float hueres = xatan2f(bufcolreserv->b[y][x], bufcolreserv->a[y][x]);
                                        const float2 sincosval4 = xsincosf(hueres);
                                        const float chrores = std::sqrt(SQR(bufcolreserv->a[y][x]) + SQR(bufcolreserv->b[y][x]));
                                        buftemp->a[y][x] = chrores * sincosval4.y;
                                        buftemp->b[y][x] = chrores * sincosval4.x;
                                        buftemp->L[y][x] = bufprov->L[y][x];
                                    }

                                    if (lp.mergemet == 4) {
                                        bufcolfin->L[y][x] = intp(lp.opacol, bufprov->L[y][x], bufcolfin->L[y][x]);
                                        bufcolfin->a[y][x] = intp(lp.opacol, bufprov->a[y][x], bufcolfin->a[y][x]);
                                        bufcolfin->b[y][x] = intp(lp.opacol, bufprov->b[y][x], bufcolfin->b[y][x]);
                                    } else {
                                        bufcolfin->L[y][x] = intp(lp.opacol, bufprov->L[y][x], bufcolreserv->L[y][x]);
                                        bufcolfin->a[y][x] = intp(lp.opacol, bufprov->a[y][x], bufcolreserv->a[y][x]);
                                        bufcolfin->b[y][x] = intp(lp.opacol, bufprov->b[y][x], bufcolreserv->b[y][x]);
                                    }
                                }
                            }

                            if (conthr > 0.f && lp.mergemet != 4) {
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        bufcolfin->L[y][x] = intp(blend[y][x], bufcolfin->L[y][x], bufcolreserv->L[y][x]);
                                        bufcolfin->a[y][x] = intp(blend[y][x], bufcolfin->a[y][x], bufcolreserv->a[y][x]);
                                        bufcolfin->b[y][x] = intp(blend[y][x], bufcolfin->b[y][x], bufcolreserv->b[y][x]);
                                    }
                                }
                            }
                        }


                        if (lp.mergecolMethod > 0 && lp.mergecolMethod <= 16) {
                            //first deltaE
#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int y = 0; y < bfh ; y++) {
                                for (int x = 0; x < bfw; x++) {
                                    bufcolfin->L[y][x] = intp(rdE[y][x], bufcolfin->L[y][x], bufcolreserv->L[y][x]);
                                    bufcolfin->a[y][x] = intp(rdE[y][x], bufcolfin->a[y][x], bufcolreserv->a[y][x]);
                                    bufcolfin->b[y][x] = intp(rdE[y][x], bufcolfin->b[y][x], bufcolreserv->b[y][x]);
                                }
                            }

                            //prepare RGB values in 0 1(or more)for current image and reserved
                            std::unique_ptr<Imagefloat> tmpImageorig(new Imagefloat(bfw, bfh));
                            lab2rgb(*bufcolfin, *tmpImageorig, params->icm.workingProfile);
                            tmpImageorig->normalizeFloatTo1();

                            std::unique_ptr<Imagefloat> tmpImagereserv(new Imagefloat(bfw, bfh));
                            lab2rgb(*bufcolreserv, *tmpImagereserv, params->icm.workingProfile);
                            tmpImagereserv->normalizeFloatTo1();

                            float minR = tmpImagereserv->r(0, 0);
                            float maxR = minR;
                            float minG = tmpImagereserv->g(0, 0);
                            float maxG = minG;
                            float minB = tmpImagereserv->b(0, 0);
                            float maxB = minB;
                            if (lp.mergecolMethod == 6 || lp.mergecolMethod == 9 || lp.mergecolMethod == 10 || lp.mergecolMethod == 11) {
#ifdef _OPENMP
                                #pragma omp parallel for reduction(max:maxR,maxG,maxB) reduction(min:minR,minG,minB) schedule(dynamic,16) if (multiThread)
#endif
                                for (int ir = 0; ir < bfh; ir++) {
                                    for (int jr = 0; jr < bfw; jr++) {
                                        minR = rtengine::min(minR, tmpImagereserv->r(ir, jr));
                                        maxR = rtengine::max(maxR, tmpImagereserv->r(ir, jr));
                                        minG = rtengine::min(minG, tmpImagereserv->g(ir, jr));
                                        maxG = rtengine::max(maxG, tmpImagereserv->g(ir, jr));
                                        minB = rtengine::min(minB, tmpImagereserv->b(ir, jr));
                                        maxB = rtengine::max(maxB, tmpImagereserv->b(ir, jr));
                                    }
                                }
                            }

                            //various combinations  subtract, multiply, difference, etc
                            if (lp.mergecolMethod == 1) { //subtract
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {//LIM(x 0 2) 2 arbitrary value but limit...
                                    for (int x = 0; x < bfw; x++) {
                                        tmpImageorig->r(y, x) = intp(lp.opacol, tmpImageorig->r(y, x) - tmpImagereserv->r(y, x), tmpImageorig->r(y, x));
                                        tmpImageorig->g(y, x) = intp(lp.opacol, tmpImageorig->g(y, x) - tmpImagereserv->g(y, x), tmpImageorig->g(y, x));
                                        tmpImageorig->b(y, x) = intp(lp.opacol, tmpImageorig->b(y, x) - tmpImagereserv->b(y, x), tmpImageorig->b(y, x));
                                    }
                                }
                            } else if (lp.mergecolMethod == 2) { //difference
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        tmpImageorig->r(y, x) = intp(lp.opacol, std::fabs(tmpImageorig->r(y, x) - tmpImagereserv->r(y, x)), tmpImageorig->r(y, x));
                                        tmpImageorig->g(y, x) = intp(lp.opacol, std::fabs(tmpImageorig->g(y, x) - tmpImagereserv->g(y, x)), tmpImageorig->g(y, x));
                                        tmpImageorig->b(y, x) = intp(lp.opacol, std::fabs(tmpImageorig->b(y, x) - tmpImagereserv->b(y, x)), tmpImageorig->b(y, x));
                                    }
                                }
                            } else if (lp.mergecolMethod == 3) { //multiply
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        tmpImageorig->r(y, x) = intp(lp.opacol, tmpImageorig->r(y, x) * tmpImagereserv->r(y, x), tmpImageorig->r(y, x));
                                        tmpImageorig->g(y, x) = intp(lp.opacol, tmpImageorig->g(y, x) * tmpImagereserv->g(y, x), tmpImageorig->g(y, x));
                                        tmpImageorig->b(y, x) = intp(lp.opacol, tmpImageorig->b(y, x) * tmpImagereserv->b(y, x), tmpImageorig->b(y, x));
                                    }
                                }
                            } else if (lp.mergecolMethod == 4) { //addition
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        tmpImageorig->r(y, x) = intp(lp.opacol, tmpImageorig->r(y, x) + tmpImagereserv->r(y, x), tmpImageorig->r(y, x));
                                        tmpImageorig->g(y, x) = intp(lp.opacol, tmpImageorig->g(y, x) + tmpImagereserv->g(y, x), tmpImageorig->g(y, x));
                                        tmpImageorig->b(y, x) = intp(lp.opacol, tmpImageorig->b(y, x) + tmpImagereserv->b(y, x), tmpImageorig->b(y, x));
                                    }
                                }
                            } else if (lp.mergecolMethod == 5) { //divide
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        tmpImageorig->r(y, x) = intp(lp.opacol, tmpImageorig->r(y, x) / (tmpImagereserv->r(y, x) + 0.00001f), tmpImageorig->r(y, x));
                                        tmpImageorig->g(y, x) = intp(lp.opacol, tmpImageorig->g(y, x) / (tmpImagereserv->g(y, x) + 0.00001f), tmpImageorig->g(y, x));
                                        tmpImageorig->b(y, x) = intp(lp.opacol, tmpImageorig->b(y, x) / (tmpImagereserv->b(y, x) + 0.00001f), tmpImageorig->b(y, x));
                                    }
                                }
                            } else if (lp.mergecolMethod == 6) { //soft light as Photoshop
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        tmpImageorig->r(y, x) = intp(lp.opacol, softlig(tmpImageorig->r(y, x), tmpImagereserv->r(y, x), minR, maxR), tmpImageorig->r(y, x));
                                        tmpImageorig->g(y, x) = intp(lp.opacol, softlig(tmpImageorig->g(y, x), tmpImagereserv->g(y, x), minG, maxG), tmpImageorig->g(y, x));
                                        tmpImageorig->b(y, x) = intp(lp.opacol, softlig(tmpImageorig->b(y, x), tmpImagereserv->b(y, x), minB, maxB), tmpImageorig->b(y, x));
                                    }
                                }
                            } else if (lp.mergecolMethod == 7) { //soft light as illusions.hu
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        tmpImageorig->r(y, x) = intp(lp.opacol, softlig2(LIM01(tmpImageorig->r(y, x)), LIM01(tmpImageorig->r(y, x))), tmpImageorig->r(y, x));
                                        tmpImageorig->g(y, x) = intp(lp.opacol, softlig2(LIM01(tmpImageorig->g(y, x)), LIM01(tmpImageorig->g(y, x))), tmpImageorig->g(y, x));
                                        tmpImageorig->b(y, x) = intp(lp.opacol, softlig2(LIM01(tmpImageorig->b(y, x)), LIM01(tmpImageorig->b(y, x))), tmpImageorig->b(y, x));
                                    }
                                }
                            } else if (lp.mergecolMethod == 8) { //soft light as W3C
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        tmpImageorig->r(y, x) = intp(lp.opacol, softlig3(LIM01(tmpImageorig->r(y, x)), tmpImagereserv->r(y, x)), tmpImageorig->r(y, x));
                                        tmpImageorig->g(y, x) = intp(lp.opacol, softlig3(LIM01(tmpImageorig->g(y, x)), tmpImagereserv->g(y, x)), tmpImageorig->g(y, x));
                                        tmpImageorig->b(y, x) = intp(lp.opacol, softlig3(LIM01(tmpImageorig->b(y, x)), tmpImagereserv->b(y, x)), tmpImageorig->b(y, x));
                                    }
                                }
                            } else if (lp.mergecolMethod == 9) { //hard light overlay (float &b, float &a)
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        tmpImageorig->r(y, x) = intp(lp.opacol, overlay(tmpImagereserv->r(y, x), tmpImageorig->r(y, x), minR, maxR), tmpImageorig->r(y, x));
                                        tmpImageorig->g(y, x) = intp(lp.opacol, overlay(tmpImagereserv->g(y, x), tmpImageorig->g(y, x), minG, maxG), tmpImageorig->g(y, x));
                                        tmpImageorig->b(y, x) = intp(lp.opacol, overlay(tmpImagereserv->b(y, x), tmpImageorig->b(y, x), minB, maxB), tmpImageorig->b(y, x));
                                    }
                                }
                            } else if (lp.mergecolMethod == 10) { //overlay overlay(float &a, float &b)
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        tmpImageorig->r(y, x) = intp(lp.opacol, overlay(tmpImageorig->r(y, x), tmpImagereserv->r(y, x), minR, maxR), tmpImageorig->r(y, x));
                                        tmpImageorig->g(y, x) = intp(lp.opacol, overlay(tmpImageorig->g(y, x), tmpImagereserv->g(y, x), minG, maxG), tmpImageorig->g(y, x));
                                        tmpImageorig->b(y, x) = intp(lp.opacol, overlay(tmpImageorig->b(y, x), tmpImagereserv->b(y, x), minB, maxB), tmpImageorig->b(y, x));
                                    }
                                }
                            } else if (lp.mergecolMethod == 11) { //screen screen (float &a, float &b)
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        tmpImageorig->r(y, x) = intp(lp.opacol, screen(tmpImageorig->r(y, x), tmpImagereserv->r(y, x), 1.f), tmpImageorig->r(y, x));
                                        tmpImageorig->g(y, x) = intp(lp.opacol, screen(tmpImageorig->g(y, x), tmpImagereserv->g(y, x), 1.f), tmpImageorig->g(y, x));
                                        tmpImageorig->b(y, x) = intp(lp.opacol, screen(tmpImageorig->b(y, x), tmpImagereserv->b(y, x), 1.f), tmpImageorig->b(y, x));
                                    }
                                }
                            } else if (lp.mergecolMethod == 12) { //darken only
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        tmpImageorig->r(y, x) = intp(lp.opacol, rtengine::min(tmpImageorig->r(y, x), tmpImagereserv->r(y, x)), tmpImageorig->r(y, x));
                                        tmpImageorig->g(y, x) = intp(lp.opacol, rtengine::min(tmpImageorig->g(y, x), tmpImagereserv->g(y, x)), tmpImageorig->g(y, x));
                                        tmpImageorig->b(y, x) = intp(lp.opacol, rtengine::min(tmpImageorig->b(y, x), tmpImagereserv->b(y, x)), tmpImageorig->b(y, x));
                                    }
                                }
                            } else if (lp.mergecolMethod == 13) { //lighten only
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        tmpImageorig->r(y, x) = intp(lp.opacol, rtengine::max(tmpImageorig->r(y, x), tmpImagereserv->r(y, x)), tmpImageorig->r(y, x));
                                        tmpImageorig->g(y, x) = intp(lp.opacol, rtengine::max(tmpImageorig->g(y, x), tmpImagereserv->g(y, x)), tmpImageorig->g(y, x));
                                        tmpImageorig->b(y, x) = intp(lp.opacol, rtengine::max(tmpImageorig->b(y, x), tmpImagereserv->b(y, x)), tmpImageorig->b(y, x));
                                    }
                                }
                            } else if (lp.mergecolMethod == 14) { //exclusion exclusion (float &a, float &b)
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        tmpImageorig->r(y, x) = intp(lp.opacol, exclusion(tmpImageorig->r(y, x), tmpImagereserv->r(y, x)), tmpImageorig->r(y, x));
                                        tmpImageorig->g(y, x) = intp(lp.opacol, exclusion(tmpImageorig->g(y, x), tmpImagereserv->g(y, x)), tmpImageorig->g(y, x));
                                        tmpImageorig->b(y, x) = intp(lp.opacol, exclusion(tmpImageorig->b(y, x), tmpImagereserv->b(y, x)), tmpImageorig->b(y, x));
                                    }
                                }

                            } else if (lp.mergecolMethod == 15) { //Color burn
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        tmpImageorig->r(y, x) = intp(lp.opacol, colburn(LIM01(tmpImageorig->r(y, x)), LIM01(tmpImagereserv->r(y, x))), tmpImageorig->r(y, x));
                                        tmpImageorig->g(y, x) = intp(lp.opacol, colburn(LIM01(tmpImageorig->g(y, x)), LIM01(tmpImagereserv->g(y, x))), tmpImageorig->g(y, x));
                                        tmpImageorig->b(y, x) = intp(lp.opacol, colburn(LIM01(tmpImageorig->b(y, x)), LIM01(tmpImagereserv->b(y, x))), tmpImageorig->b(y, x));
                                    }
                                }
                            } else if (lp.mergecolMethod == 16) { //Color dodge
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        tmpImageorig->r(y, x) = intp(lp.opacol, coldodge(LIM01(tmpImageorig->r(y, x)), LIM01(tmpImagereserv->r(y, x))), tmpImageorig->r(y, x));
                                        tmpImageorig->g(y, x) = intp(lp.opacol, coldodge(LIM01(tmpImageorig->g(y, x)), LIM01(tmpImagereserv->g(y, x))), tmpImageorig->g(y, x));
                                        tmpImageorig->b(y, x) = intp(lp.opacol, coldodge(LIM01(tmpImageorig->b(y, x)), LIM01(tmpImagereserv->b(y, x))), tmpImageorig->b(y, x));
                                    }
                                }
                            }

                            tmpImageorig->normalizeFloatTo65535();
                            rgb2lab(*tmpImageorig, *bufcolfin, params->icm.workingProfile);

                            if (conthr > 0.f && lp.mergemet != 4) {
#ifdef _OPENMP
                                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                                for (int y = 0; y < bfh ; y++) {
                                    for (int x = 0; x < bfw; x++) {
                                        bufcolfin->L[y][x] = intp(blend[y][x], bufcolfin->L[y][x], bufcolreserv->L[y][x]);
                                        bufcolfin->a[y][x] = intp(blend[y][x], bufcolfin->a[y][x], bufcolreserv->a[y][x]);
                                        bufcolfin->b[y][x] = intp(blend[y][x], bufcolfin->b[y][x], bufcolreserv->b[y][x]);
                                    }
                                }
                            }
                        }

                        if (lp.softradiuscol > 0.f) {
                            softproc(bufcolreserv.get(), bufcolfin.get(), lp.softradiuscol, bfh, bfw, 0.001, 0.00001, 0.5f, sk, multiThread, 1);
                        }
                        float meansob = 0.f;
                        const float repart = 1.0 - 0.01 * params->locallab.spots.at(sp).reparcol;
                        int bw = bufcolreserv->W;
                        int bh = bufcolreserv->H;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if(multiThread)
#endif
                    for (int x = 0; x < bh; x++) {
                        for (int y = 0; y < bw; y++) {
                            bufcolfin->L[x][y] = intp(repart, bufcolreserv->L[x][y], bufcolfin->L[x][y]);
                            bufcolfin->a[x][y] = intp(repart, bufcolreserv->a[x][y], bufcolfin->a[x][y]);
                            bufcolfin->b[x][y] = intp(repart, bufcolreserv->b[x][y], bufcolfin->b[x][y]);
                        }
                    }
                        
                        transit_shapedetect2(sp, 0.f, 0.f, call, 0, bufcolreserv.get(), bufcolfin.get(), originalmaskcol.get(), hueref, chromaref, lumaref, sobelref, meansob, blend2, lp, original, transformed, cx, cy, sk);
                    }

                    if (!nottransit) {
//gradient
                        if (lp.strcol != 0.f) {
                            struct grad_params gp;
                            calclocalGradientParams(lp, gp, ystart, xstart, bfw, bfh, 3);
#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int ir = 0; ir < bfh; ir++)
                                for (int jr = 0; jr < bfw; jr++) {
                                    const float corrFactor = ImProcFunctions::calcGradientFactor(gp, jr, ir);
                                    bufcolfin->L[ir][jr] *= corrFactor;
                                }
                        }

                        if (lp.strcolab != 0.f) {
                            struct grad_params gpab;
                            calclocalGradientParams(lp, gpab, ystart, xstart, bfw, bfh, 5);
#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int ir = 0; ir < bfh; ir++)
                                for (int jr = 0; jr < bfw; jr++) {
                                    const float corrFactor = ImProcFunctions::calcGradientFactor(gpab, jr, ir);
                                    bufcolfin->a[ir][jr] *= corrFactor;
                                    bufcolfin->b[ir][jr] *= corrFactor;
                                }
                        }

                        if (lp.strcolh != 0.f) {
                            struct grad_params gph;
                            calclocalGradientParams(lp, gph, ystart, xstart, bfw, bfh, 6);
#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int ir = 0; ir < bfh; ir++)
                                for (int jr = 0; jr < bfw; jr++) {
                                    const float corrFactor = ImProcFunctions::calcGradientFactor(gph, jr, ir);
                                    const float aa = bufcolfin->a[ir][jr];
                                    const float bb = bufcolfin->b[ir][jr];
                                    const float chrm = std::sqrt(SQR(aa) + SQR(bb));
                                    const float HH = xatan2f(bb, aa);

                                    float cor = 0.f;

                                    if (corrFactor < 1.f) {
                                        cor = - 2.5f * (1.f - corrFactor);
                                    } else if (corrFactor > 1.f) {
                                        cor = 0.03f * (corrFactor - 1.f);
                                    }

                                    float newhr = HH + cor;

                                    if (newhr > rtengine::RT_PI_F) {
                                        newhr -= 2 * rtengine::RT_PI_F;
                                    } else if (newhr < -rtengine::RT_PI_F) {
                                        newhr += 2 * rtengine::RT_PI_F;
                                    }

                                    const float2 sincosval = xsincosf(newhr);
                                    bufcolfin->a[ir][jr] = clipC(chrm * sincosval.y);
                                    bufcolfin->b[ir][jr] = clipC(chrm * sincosval.x);
                                }
                        }



/*
                        float gamma = lp.gamc;
                        rtengine::GammaValues g_a; //gamma parameters
                        double pwr = 1.0 / (double) lp.gamc;//default 3.0 - gamma Lab
                        double ts = 9.03296;//always the same 'slope' in the extreme shadows - slope Lab
                        rtengine::Color::calcGamma(pwr, ts, g_a); // call to calcGamma with selected gamma and slope
*/
                        if(gamma1 != 1.f) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif  
                            for (int y = 0; y < bfh; ++y) {//apply inverse gamma 3.f and put result in range 32768.f
                                int x = 0;
#ifdef __SSE2__
                                for (; x < bfw - 3; x += 4) {
                                STVFU(bufcolfin->L[y][x], F2V(32768.f) * gammalog(LVFU(bufcolfin->L[y][x]) / F2V(32768.f), F2V(gamma1), F2V(ts1), F2V(g_a[3]), F2V(g_a[4])));
                                }
#endif
                                for (; x < bfw; ++x) {
                                    bufcolfin->L[y][x] = 32768.f * gammalog(bufcolfin->L[y][x] / 32768.f, gamma1, ts1, g_a[3], g_a[4]);
                                }
                            }
                        }


                        if (lp.softradiuscol > 0.f) {
                            softproc(bufcolorig.get(), bufcolfin.get(), lp.softradiuscol, bfh, bfw, 0.001, 0.00001, 0.5f, sk, multiThread, 1);
                        }
                    //mask recovery

                    if(lp.enaColorMask && lp.recothrc != 1.f) {
                        float recoth = lp.recothrc;

                        if(lp.recothrc < 1.f) {
                            recoth = -1.f * recoth + 2.f;
                        }

                        float hig = lp.higthrc;
                        float low = lp.lowthrc;
                      //  float recoth = lp.recothrc;
                        float decay = lp.decayc;
                        bool invmask = false;
                        maskrecov(bufcolfin.get(), original, bufmaskblurcol.get(), bfh, bfw, ystart, xstart, hig, low, recoth, decay, invmask, sk, multiThread);
                    }
                    const float repart = 1.0 - 0.01 * params->locallab.spots.at(sp).reparcol;
                    int bw = bufcolorig->W;
                    int bh = bufcolorig->H;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if(multiThread)
#endif
                    for (int x = 0; x < bh; x++) {
                        for (int y = 0; y < bw; y++) {
                            bufcolfin->L[x][y] = intp(repart, bufcolorig->L[x][y], bufcolfin->L[x][y]);
                            bufcolfin->a[x][y] = intp(repart, bufcolorig->a[x][y], bufcolfin->a[x][y]);
                            bufcolfin->b[x][y] = intp(repart, bufcolorig->b[x][y], bufcolfin->b[x][y]);
                        }
                    }
                    
                        float meansob = 0.f;
                        if(lp.recothrc >= 1.f) {
                            transit_shapedetect2(sp, 0.f, 0.f, call, 0, bufcolorig.get(), bufcolfin.get(), originalmaskcol.get(), hueref, chromaref, lumaref, sobelref, meansob, blend2, lp, original, transformed, cx, cy, sk);
                        } else {
                            transit_shapedetect2(sp, 0.f, 0.f, call, 0, bufcolorig.get(), bufcolfin.get(), nullptr, hueref, chromaref, lumaref, sobelref, meansob, blend2, lp, original, transformed, cx, cy, sk);
                        }
                    }

                }

                if (lp.recur) {
                    original->CopyFrom(transformed, multiThread);
                    float avge;
                    calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                }

            }
        }
    }

//inverse
    else if (lp.inv && (lp.chro != 0 || lp.ligh != 0 || exlocalcurve || lp.showmaskcolmetinv == 0 || lp.enaColorMaskinv) && lp.colorena) {
        float adjustr = 1.0f;

//adapt chroma to working profile
        if (params->icm.workingProfile == "ProPhoto")   {
            adjustr = 1.2f;   // 1.2 instead 1.0 because it's very rare to have C>170..
        } else if (params->icm.workingProfile == "Adobe RGB")  {
            adjustr = 1.8f;
        } else if (params->icm.workingProfile == "sRGB")       {
            adjustr = 2.0f;
        } else if (params->icm.workingProfile == "WideGamut")  {
            adjustr = 1.2f;
        } else if (params->icm.workingProfile == "Beta RGB")   {
            adjustr = 1.4f;
        } else if (params->icm.workingProfile == "BestRGB")    {
            adjustr = 1.4f;
        } else if (params->icm.workingProfile == "BruceRGB")   {
            adjustr = 1.8f;
        }

        std::unique_ptr<LabImage> bufmaskblurcol;
        std::unique_ptr<LabImage> originalmaskcol;
        const std::unique_ptr<LabImage> bufcolorig(new LabImage(TW, TH));

        if (lp.enaColorMaskinv || lp.showmaskcolmetinv == 1) {
            bufmaskblurcol.reset(new LabImage(TW, TH, true));
            originalmaskcol.reset(new LabImage(TW, TH));
        }

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
        for (int y = 0; y < TH ; y++) {
            for (int x = 0; x < TW; x++) {
                bufcolorig->L[y][x] = original->L[y][x];
            }
        }

        constexpr int inv = 1;
        const bool showmaske = lp.showmaskcolmetinv == 1;
        const bool enaMask = lp.enaColorMaskinv;
        constexpr bool deltaE = false;
        constexpr bool modmask = false;
        const bool zero = lp.showmaskcolmetinv == 0;
        constexpr bool modif = false;

        const float chrom = lp.chromacol;
        const float rad = lp.radmacol;
        const float gamma = lp.gammacol;
        const float slope = lp.slomacol;
        const float blendm = lp.blendmacol;
        const float lap = params->locallab.spots.at(sp).lapmaskcol;
        const bool pde = params->locallab.spots.at(sp).laplac;
        int shado = params->locallab.spots.at(sp).shadmaskcol;
        int level_bl = params->locallab.spots.at(sp).csthresholdcol.getBottomLeft();
        int level_hl = params->locallab.spots.at(sp).csthresholdcol.getTopLeft();
        int level_br = params->locallab.spots.at(sp).csthresholdcol.getBottomRight();
        int level_hr = params->locallab.spots.at(sp).csthresholdcol.getTopRight();
        int sco = params->locallab.spots.at(sp).scopemask;
        int shortcu = lp.mergemet; //params->locallab.spots.at(sp).shortc;
        int lumask = params->locallab.spots.at(sp).lumask;
        float strumask = 0.02 * params->locallab.spots.at(sp).strumaskcol;

        const float mindE = 2.f + MINSCOPE * sco * lp.thr;
        const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
        const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
        const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
        constexpr float amountcd = 0.f;
        constexpr float anchorcd = 50.f;
        const int highl = 0;
        maskcalccol(false, pde, TW, TH, 0, 0, sk, cx, cy, bufcolorig.get(), bufmaskblurcol.get(), originalmaskcol.get(), original, reserved, inv, lp,
                    strumask, params->locallab.spots.at(sp).toolcol,
                    locccmasCurve, lcmasutili, locllmasCurve, llmasutili, lochhmasCurve, lhmasutili, llochhhmasCurve, lhhmasutili, multiThread,
                    enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, highl, amountcd, anchorcd, lmasklocalcurve, localmaskutili, loclmasCurvecolwav, lmasutilicolwav,
                    level_bl, level_hl, level_br, level_hr,
                    shortcu, false, hueref, chromaref, lumaref,
                    maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, lp.fftColorMask, lp.blurcolmask, lp.contcolmask, -1, fab
                   );

        if (lp.showmaskcolmetinv == 1) {
            showmask(lumask, lp, 0, 0, cx, cy, TW, TH, bufcolorig.get(), transformed, bufmaskblurcol.get(), inv);
            return;
        }

        if (lp.showmaskcolmetinv == 0 || lp.enaColorMaskinv) {
            InverseColorLight_Local(false, false, sp, 0, lp, originalmaskcol.get(), lightCurveloc, hltonecurveloc, shtonecurveloc, tonecurveloc, exlocalcurve, cclocalcurve, adjustr, localcutili, lllocalcurve, locallutili, original, transformed, cx, cy, hueref, chromaref, lumaref, sk);

            if (lp.recur) {
                original->CopyFrom(transformed, multiThread);
                float avge;
                calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
            }
        }
    }
    
//begin common mask
    if(lp.maskena) {
        int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
        int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
        int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
        int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
        int bfh = yend - ystart;
        int bfw = xend - xstart;
        if (bfw >= mSP && bfh >= mSP) {

            if (lp.blurma >= 0.25f && lp.fftma && call == 2) {
                optfft(N_fftwsize, bfh, bfw, bfh, bfw, lp, original->H, original->W, xstart, ystart, xend, yend, cx, cy, lp.fullim);
            }
            array2D<float> blechro(bfw, bfh);
            array2D<float> ble(bfw, bfh);
            array2D<float> hue(bfw, bfh);
            array2D<float> guid(bfw, bfh);

            std::unique_ptr<LabImage> bufcolorigsav;
            std::unique_ptr<LabImage> bufcolorig;
            std::unique_ptr<LabImage> bufcolfin;
            std::unique_ptr<LabImage> bufmaskblurcol;
            std::unique_ptr<LabImage> originalmaskcol;
            std::unique_ptr<LabImage> bufcolreserv;
            
            int wo = original->W;
            int ho = original->H;
            LabImage *origsav = nullptr;
            origsav = new LabImage(wo, ho);
            origsav->CopyFrom(original);

            if (call <= 3) {
                bufcolorig.reset(new LabImage(bfw, bfh));
                bufcolfin.reset(new LabImage(bfw, bfh));
                bufcolorigsav.reset(new LabImage(bfw, bfh));

                if (lp.showmask_met == 1  || lp.ena_Mask || lp.showmask_met == 2 || lp.showmask_met == 3) {
                    bufmaskblurcol.reset(new LabImage(bfw, bfh, true));
                    originalmaskcol.reset(new LabImage(bfw, bfh));
                }
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                for (int y = 0; y < bfh ; y++) {
                    for (int x = 0; x < bfw; x++) {
                        bufcolorig->L[y][x] = original->L[y + ystart][x + xstart];
                        bufcolorig->a[y][x] = original->a[y + ystart][x + xstart];
                        bufcolorig->b[y][x] = original->b[y + ystart][x + xstart];

                        bufcolorigsav->L[y][x] = original->L[y + ystart][x + xstart];
                        bufcolorigsav->a[y][x] = original->a[y + ystart][x + xstart];
                        bufcolorigsav->b[y][x] = original->b[y + ystart][x + xstart];

                        bufcolfin->L[y][x] = original->L[y + ystart][x + xstart];
                        bufcolfin->a[y][x] = original->a[y + ystart][x + xstart];
                        bufcolfin->b[y][x] = original->b[y + ystart][x + xstart];
                    }
                }
                const int inv = 0;
                const bool showmaske = lp.showmask_met == 2;
                const bool enaMask = lp.ena_Mask;
                const bool deltaE = lp.showmask_met == 3;
                const bool modmask = lp.showmask_met == 1;
                const bool zero = lp.showmask_met == 0;
                const bool modif = lp.showmask_met == 1;
                const float chrom = params->locallab.spots.at(sp).chromask;
                const float rad = params->locallab.spots.at(sp).radmask; 
                const float gamma = params->locallab.spots.at(sp).gammask; 
                const float slope =  params->locallab.spots.at(sp).slopmask;
                float blendm =  0.1 * params->locallab.spots.at(sp).blendmask;
                float blendmab =  params->locallab.spots.at(sp).blendmaskab;
                if (lp.showmask_met == 2) {
                    blendm = 0.f;//normalize behavior mask with others no action of blend
                    blendmab = 0.f;
                }
                const float lap = params->locallab.spots.at(sp).lapmask;
                const bool pde = params->locallab.spots.at(sp).laplac;
                const int shado = params->locallab.spots.at(sp).shadmask;
                const int sco = params->locallab.spots.at(sp).scopemask;
                const int level_bl = params->locallab.spots.at(sp).csthresholdmask.getBottomLeft();
                const int level_hl = params->locallab.spots.at(sp).csthresholdmask.getTopLeft();
                const int level_br = params->locallab.spots.at(sp).csthresholdmask.getBottomRight();
                const int level_hr = params->locallab.spots.at(sp).csthresholdmask.getTopRight();
                const int shortcu = lp.mergemet; //params->locallab.spots.at(sp).shortc;
                const int lumask = params->locallab.spots.at(sp).lumask;
                const float strumask = 0.02 * params->locallab.spots.at(sp).strumaskmask;
                const float softr = params->locallab.spots.at(sp).softradiusmask;

                const float mindE = 2.f + MINSCOPE * sco * lp.thr;
                const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
                const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
                const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
                const float amountcd = 0.f;
                const float anchorcd = 50.f;
                const int highl = 0;
                bool astool = params->locallab.spots.at(sp).toolmask;
                maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, bufcolorig.get(), bufmaskblurcol.get(), originalmaskcol.get(), original, reserved, inv, lp,
                            strumask, astool,
                            locccmas_Curve, lcmas_utili, locllmas_Curve, llmas_utili, lochhmas_Curve, lhmas_utili, lochhhmas_Curve, lhhmas_utili, multiThread,
                            enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendmab, shado, highl, amountcd, anchorcd, lmasklocal_curve, localmask_utili, loclmasCurve_wav, lmasutili_wav,
                            level_bl, level_hl, level_br, level_hr,
                            shortcu, delt, hueref, chromaref, lumaref,
                            maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, lp.fftma, lp.blurma, lp.contma, 12, fab
                           );
                           
                           
                if (lp.showmask_met == 2) {
                    showmask(lumask, lp, xstart, ystart, cx, cy, bfw, bfh, bufcolorig.get(), transformed, bufmaskblurcol.get(), 0);
                    return;
                }
#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                for (int y = 0; y < bfh ; y++) {
                    for (int x = 0; x < bfw; x++) {
                        bufcolfin->L[y][x] = bufcolorig->L[y][x];
                        bufcolfin->a[y][x] = bufcolorig->a[y][x];
                        bufcolfin->b[y][x] = bufcolorig->b[y][x];
                        hue[y][x] = xatan2f(bufcolfin->b[y][x], bufcolfin->a[y][x]);
                        const float chromah = std::sqrt(SQR(bufcolfin->b[y][x]) + SQR(bufcolfin->a[y][x]));
                        ble[y][x] = bufcolfin->L[y][x] / 32768.f;
                        blechro[y][x] = chromah / 32768.f;
                        guid[y][x] = bufcolorigsav->L[y][x] / 32768.f;
                    }
                }
                if (softr != 0.f) {//soft for L a b because we change color...
                    const float tmpblur = softr < 0.f ? -1.f / softr : 1.f + softr;
                    const int r1 = rtengine::max<int>(4 / sk * tmpblur + 0.5f, 1);
                    const int r2 = rtengine::max<int>(25 / sk * tmpblur + 0.5f, 1);

                    constexpr float epsilmax = 0.005f;
                    constexpr float epsilmin = 0.00001f;

                    constexpr float aepsil = (epsilmax - epsilmin) / 100.f;
                    constexpr float bepsil = epsilmin;
                    const float epsil = softr < 0.f ? 0.001f : aepsil * softr + bepsil;

                    rtengine::guidedFilter(guid, blechro, blechro, r1, epsil, multiThread);
                    rtengine::guidedFilter(guid, ble, ble, r2, 0.2f * epsil, multiThread);

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                    for (int y = 0; y < bfh; y++) {
                        for (int x = 0; x < bfw; x++) {
                            float2 sincosval = xsincosf(hue[y][x]);
                            bufcolfin->L[y][x] = 32768.f * ble[y][x];
                            bufcolfin->a[y][x] = 32768.f * blechro[y][x] * sincosval.y;
                            bufcolfin->b[y][x] = 32768.f * blechro[y][x] * sincosval.x;
                        }
                    }
                }
                


                float meansob = 0.f;
                transit_shapedetect2(sp, 0.f, 0.f, call, 20, bufcolorigsav.get(), bufcolfin.get(), originalmaskcol.get(), hueref, chromaref, lumaref, sobelref, meansob, nullptr, lp, origsav, transformed, cx, cy, sk);
                delete origsav;
                origsav    = NULL;
                    
                if (lp.recur) {
                    original->CopyFrom(transformed, multiThread);
                    float avge;
                    calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                }
                    
            }
        }
    }
//end common mask
    
        if(params->locallab.spots.at(sp).expcie  && params->locallab.spots.at(sp).modecie == "com"  && lp.activspot) {//ciecam
            int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
            int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
            int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
            int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
            int bfh = yend - ystart;
            int bfw = xend - xstart;

            if (bfh >= mSP && bfw >= mSP) {
                const std::unique_ptr<LabImage> bufexporig(new LabImage(bfw, bfh)); //buffer for data in zone limit
                const std::unique_ptr<LabImage> bufexpfin(new LabImage(bfw, bfh)); //buffer for data in zone limit
                std::unique_ptr<LabImage> bufmaskorigcie;
                std::unique_ptr<LabImage> bufmaskblurcie;
                std::unique_ptr<LabImage> originalmaskcie;

            if (lp.showmaskciemet == 2  || lp.enacieMask || lp.showmaskciemet == 3 || lp.showmaskciemet == 4) {
                bufmaskorigcie.reset(new LabImage(bfw, bfh));
                bufmaskblurcie.reset(new LabImage(bfw, bfh));
                originalmaskcie.reset(new LabImage(bfw, bfh));
            }



#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif

                for (int y = 0; y < bfh; y++) {
                    for (int x = 0; x < bfw; x++) {
                        bufexporig->L[y][x] = original->L[y + ystart][x + xstart];
                        bufexporig->a[y][x] = original->a[y + ystart][x + xstart];
                        bufexporig->b[y][x] = original->b[y + ystart][x + xstart];
                    }
                }

                bool HHcurvejz = false, CHcurvejz = false, LHcurvejz = false;
                if (params->locallab.spots.at(sp).expcie  && params->locallab.spots.at(sp).modecam == "jz") {//some cam16 elementsfor Jz
                    ImProcFunctions::ciecamloc_02float(lp, sp, bufexporig.get(), bfw, bfh, 10, sk, cielocalcurve, localcieutili, cielocalcurve2, localcieutili2, jzlocalcurve, localjzutili, czlocalcurve, localczutili, czjzlocalcurve, localczjzutili, locchCurvejz, lochhCurvejz, loclhCurvejz, HHcurvejz, CHcurvejz, LHcurvejz, locwavCurvejz, locwavutilijz);
                }
                if (lochhCurvejz && HHutilijz) {
                    for (int i = 0; i < 500; i++) {
                        if (lochhCurvejz[i] != 0.5f) {
                            HHcurvejz = true;
                            break;
                        }
                    }
                }

                if (locchCurvejz && CHutilijz) {
                    for (int i = 0; i < 500; i++) {
                        if (locchCurvejz[i] != 0.5f) {
                            CHcurvejz = true;
                            break;
                        }
                    }
                }

                if (loclhCurvejz && LHutilijz) {
                    for (int i = 0; i < 500; i++) {
                        if (loclhCurvejz[i] != 0.5f) {
                            LHcurvejz = true;
                            break;
                        }
                    }
                }

                int inv = 0;
                bool showmaske = false;
                bool enaMask = false;
                bool deltaE = false;
                bool modmask = false;
                bool zero = false;
                bool modif = false;

                if (lp.showmaskciemet == 3) {
                    showmaske = true;
                }

                if (lp.enacieMask) {
                    enaMask = true;
                }

                if (lp.showmaskciemet == 4) {
                    deltaE = true;
                }

                if (lp.showmaskciemet == 2) {
                    modmask = true;
                }

                if (lp.showmaskciemet == 1) {
                    modif = true;
                }

                if (lp.showmaskciemet == 0) {
                    zero = true;
                }

                float chrom = lp.chromacie;
                float rad = lp.radmacie;
                float gamma = params->locallab.spots.at(sp).gammaskcie;
                float slope = params->locallab.spots.at(sp).slomaskcie;
                float blendm = lp.blendmacie;
                float lap = params->locallab.spots.at(sp).lapmaskcie;
                bool pde =  params->locallab.spots.at(sp).laplac;
                LocwavCurve dummy;
                bool delt = params->locallab.spots.at(sp).deltae;
                int sco = params->locallab.spots.at(sp).scopemask;
                int shortcu = 0;//lp.mergemet; //params->locallab.spots.at(sp).shortc;
                int shado = 0;
                const int highl = 0;

                const float mindE = 2.f + MINSCOPE * sco * lp.thr;
                const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
                const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
                const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
                int lumask = params->locallab.spots.at(sp).lumask;
                float amountcd = 0.f;
                float anchorcd = 50.f;
                LocHHmaskCurve lochhhmasCurve;
                maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, bufexporig.get(), bufmaskorigcie.get(), originalmaskcie.get(), original, reserved, inv, lp,
                        0.f, false,
                        locccmascieCurve, lcmascieutili, locllmascieCurve, llmascieutili, lochhmascieCurve, lhmascieutili, lochhhmasCurve, false, multiThread,
                        enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, highl, amountcd, anchorcd, lmaskcielocalcurve, localmaskcieutili, dummy, false, 1, 1, 5, 5,
                        shortcu, delt, hueref, chromaref, lumaref,
                        maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, -1, fab
                       );

                if (lp.showmaskciemet == 3) {
                    showmask(lumask, lp, xstart, ystart, cx, cy, bfw, bfh, bufexporig.get(), transformed, bufmaskorigcie.get(), 0);

                    return;
                }



                if (lp.showmaskciemet == 0 || lp.showmaskciemet == 1  || lp.showmaskciemet == 2 || lp.showmaskciemet == 4 || lp.enacieMask) {

                    bufexpfin->CopyFrom(bufexporig.get(), multiThread);
                    if (params->locallab.spots.at(sp).expcie) {
                        ImProcFunctions::ciecamloc_02float(lp, sp, bufexpfin.get(), bfw, bfh, 0, sk, cielocalcurve, localcieutili, cielocalcurve2, localcieutili2, jzlocalcurve, localjzutili, czlocalcurve, localczutili, czjzlocalcurve, localczjzutili, locchCurvejz, lochhCurvejz, loclhCurvejz, HHcurvejz, CHcurvejz, LHcurvejz, locwavCurvejz, locwavutilijz);
                    }
                }
            
                if(lp.enacieMask && lp.recothrcie != 1.f) {
                    float recoth = lp.recothrcie;

                    if(lp.recothrcie < 1.f) {
                        recoth = -1.f * recoth + 2.f;
                    }
                    float hig = lp.higthrcie;
                    float low = lp.lowthrcie;
                    //float recoth = lp.recothrcie;
                    float decay = lp.decaycie;
                    bool invmask = false;
                    maskrecov(bufexpfin.get(), original, bufmaskorigcie.get(), bfh, bfw, ystart, xstart, hig, low, recoth, decay, invmask, sk, multiThread);
                }
            
            
                float radcie = params->locallab.spots.at(sp).detailcie;
                loccont(bfw, bfh, bufexpfin.get(), radcie, 15.f, sk);
                
                const float repart = 1.0 - 0.01 * params->locallab.spots.at(sp).reparcie;
                int bw = bufexporig->W;
                int bh = bufexporig->H;

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if(multiThread)
#endif
                for (int x = 0; x < bh; x++) {
                    for (int y = 0; y < bw; y++) {
                        bufexpfin->L[x][y] = intp(repart, bufexporig->L[x][y], bufexpfin->L[x][y]);
                        bufexpfin->a[x][y] = intp(repart, bufexporig->a[x][y], bufexpfin->a[x][y]);
                        bufexpfin->b[x][y] = intp(repart, bufexporig->b[x][y], bufexpfin->b[x][y]);
                    }
                }

                if(lp.recothrcie >= 1.f) {
                    transit_shapedetect2(sp, 0.f, 0.f, call, 31, bufexporig.get(), bufexpfin.get(), originalmaskcie.get(), hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);
                } else {
                    transit_shapedetect2(sp, 0.f, 0.f, call, 31, bufexporig.get(), bufexpfin.get(), nullptr, hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);
                }
                if (lp.recur) {
                    original->CopyFrom(transformed, multiThread);
                    float avge;
                    calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                }
            }

        }


// Gamut and Munsell control - very important do not deactivated to avoid crash
    avoidcolshi(lp, sp, transformed, reserved, cy, cx, sk);
}

}
