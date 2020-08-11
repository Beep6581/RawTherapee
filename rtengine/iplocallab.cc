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


#pragma GCC diagnostic warning "-Wall"
#pragma GCC diagnostic warning "-Wextra"

namespace
{

constexpr int limscope = 80;
constexpr int mSPsharp = 39; //minimum size Spot Sharp due to buildblendmask
constexpr int mSPwav = 32; //minimum size Spot Wavelet
constexpr int mDEN = 64; //minimum size Spot Denoise
constexpr int mSP = 5; //minimum size Spot
constexpr float MAXSCOPE = 1.25f;
constexpr float MINSCOPE = 0.025f;
constexpr int TS = 64; // Tile size
constexpr float epsilonw = 0.001f / (TS * TS); //tolerance
constexpr int offset = 25; // shift between tiles

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
    const float bde = -ade * maxC;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif
    for (int y = 0; y < bfh; y++) {
        for (int x = 0; x < bfw; x++) {
            dE[y * bfw + x] = dE[y * bfw + x] >= maxdE ? ade * dE[y * bfw + x] + bde : 1.f;
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
    float thr;
    float stru;
    int chro, cont, sens, sensh, senscb, sensbn, senstm, sensex, sensexclu, sensden, senslc, senssf, senshs, senscolor;
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
    int showmask_met;
    bool fftbl;
    float laplacexp;
    float balanexp;
    float linear;
    int expmet;
    int softmet;
    int blurmet;
    int blmet;
    int smasktyp;
    int chromet;
    int shmeth;
    int medmet;
    int locmet;
    float noiself;
    float noiself0;
    float noiself2;
    float noiseldetail;
    int detailthr;
    int noiselequal;
    float noisechrodetail;
    float bilat;
    float noiselc;
    float noisecf;
    float noisecc;
    float mulloc[6];
    int mullocsh[5];
    int detailsh;
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
    bool expvib;
    bool exposena;
    bool hsena;
    bool vibena;
    bool logena;
    bool maskena;
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
    bool ena_Mask;
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
    bool blwh;
    bool fftma;
    float blurma;
    float contma;
    bool activspot;

};

static void calcLocalParams(int sp, int oW, int oH, const LocallabParams& locallab, struct local_params& lp, bool prevDeltaE, int llColorMask, int llColorMaskinv, int llExpMask, int llExpMaskinv, int llSHMask, int llSHMaskinv, int llvibMask, int lllcMask, int llsharMask, int llcbMask, int llretiMask, int llsoftMask, int lltmMask, int llblMask, int ll_Mask, const LocwavCurve & locwavCurveden, bool locwavdenutili)
{
    int w = oW;
    int h = oH;
    int circr = locallab.spots.at(sp).circrad;
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

    if (thre > 8.f || thre < 0.f) {//to avoid artifacts if user does not clear cache with new settings. Can be suppressed after
        thre = 2.f;
    }

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

    if (iterati > 4.f || iterati < 0.2f) {//to avoid artifacts if user does not clear cache with new settings Can be suppressed after
        iterati = 2.f;
    }

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

    if (locallab.spots.at(sp).expMethod == "std") {
        lp.expmet = 0;
    } else if (locallab.spots.at(sp).expMethod == "pde") {
        lp.expmet = 1;
    }

    if (locallab.spots.at(sp).localcontMethod == "loc") {
        lp.locmet = 0;
    } else if (locallab.spots.at(sp).localcontMethod == "wav") {
        lp.locmet = 1;
    }

    lp.laplacexp = locallab.spots.at(sp).laplacexp;
    lp.balanexp = locallab.spots.at(sp).balanexp;
    lp.linear = locallab.spots.at(sp).linear;

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
    lp.showmask_met = ll_Mask;
  //  printf("mask=%i \n", lp.showmask_met);
    
    lp.enaColorMask = locallab.spots.at(sp).enaColorMask && llsoftMask == 0 && llColorMask == 0 && lllcMask == 0 && llsharMask == 0 && llExpMask == 0 && llSHMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && ll_Mask == 0;// Exposure mask is deactivated if Color & Light mask is visible
    lp.enaColorMaskinv = locallab.spots.at(sp).enaColorMask && llColorMaskinv == 0 && llsoftMask == 0 && lllcMask == 0 && llsharMask == 0 && llExpMask == 0 && llSHMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && ll_Mask == 0;// Exposure mask is deactivated if Color & Light mask is visible
    lp.enaExpMask = locallab.spots.at(sp).enaExpMask && llExpMask == 0 && llColorMask == 0 && llsoftMask == 0 && lllcMask == 0 && llsharMask == 0 && llSHMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && ll_Mask == 0;// Exposure mask is deactivated if Color & Light mask is visible
    lp.enaExpMaskinv = locallab.spots.at(sp).enaExpMask && llExpMaskinv == 0 && llColorMask == 0 && llsoftMask == 0 && lllcMask == 0 && llsharMask == 0 && llSHMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && ll_Mask == 0;// Exposure mask is deactivated if Color & Light mask is visible
    lp.enaSHMask = locallab.spots.at(sp).enaSHMask && llSHMask == 0 && llColorMask == 0 && llsoftMask == 0 && lllcMask == 0 && llsharMask == 0 && llExpMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && ll_Mask == 0;
    lp.enaSHMaskinv = locallab.spots.at(sp).enaSHMask && llSHMaskinv == 0 && lllcMask == 0 && llsoftMask == 0 && llsharMask == 0 && llColorMask == 0 && llExpMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && ll_Mask == 0;
    lp.enacbMask = locallab.spots.at(sp).enacbMask && llcbMask == 0 && lllcMask == 0 && llsoftMask == 0 && llsharMask == 0 && llColorMask == 0 && llExpMask == 0 && llSHMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && ll_Mask == 0;
    lp.enaretiMask = locallab.spots.at(sp).enaretiMask && lllcMask == 0 && llsharMask == 0 && llsoftMask == 0 && llretiMask == 0 && llColorMask == 0 && llExpMask == 0 && llSHMask == 0 && llcbMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && ll_Mask == 0;
    lp.enatmMask = locallab.spots.at(sp).enatmMask && lltmMask == 0 && lllcMask == 0 && llsoftMask == 0 && llsharMask == 0 && llColorMask == 0 && llExpMask == 0 && llSHMask == 0 && llcbMask == 0 && llretiMask == 0 && llblMask == 0 && llvibMask == 0 && ll_Mask == 0;
    lp.enablMask = locallab.spots.at(sp).enablMask && llblMask == 0 && lllcMask == 0 && llsoftMask == 0 && llsharMask == 0 && llColorMask == 0 && llExpMask == 0 && llSHMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llvibMask == 0 && ll_Mask == 0;
    lp.enavibMask = locallab.spots.at(sp).enavibMask && llvibMask == 0 && lllcMask == 0 && llsoftMask == 0 && llsharMask == 0 && llColorMask == 0 && llExpMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llSHMask == 0 && ll_Mask == 0;
    lp.enalcMask = locallab.spots.at(sp).enalcMask && lllcMask == 0 && llcbMask == 0 && llsoftMask == 0 && llsharMask == 0 && llColorMask == 0 && llExpMask == 0 && llSHMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && ll_Mask == 0 ;
    lp.enasharMask = lllcMask == 0 && llcbMask == 0 && llsharMask == 0 && llsoftMask == 0 && llColorMask == 0 && llExpMask == 0 && llSHMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0 && ll_Mask == 0;
    lp.ena_Mask = locallab.spots.at(sp).enamask && lllcMask == 0 && llcbMask == 0 && llsoftMask == 0 && llsharMask == 0 && llColorMask == 0 && llExpMask == 0 && llSHMask == 0 && llretiMask == 0 && lltmMask == 0 && llblMask == 0 && llvibMask == 0;

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

    if (locallab.spots.at(sp).blurMethod == "norm") {
        lp.blurmet = 0;
    } else if (locallab.spots.at(sp).blurMethod == "inv") {
        lp.blurmet = 1;
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
            local_noiself = 250.f * locwavCurveden[166];
            local_noiself2 = 250.f * locwavCurveden[323];
            local_noiselc = 200.f * locwavCurveden[500];
        }
    }

    float local_noiseldetail = (float)locallab.spots.at(sp).noiselumdetail;
    int local_noiselequal = locallab.spots.at(sp).noiselequal;
    float local_noisechrodetail = (float)locallab.spots.at(sp).noisechrodetail;
    int local_sensiden = locallab.spots.at(sp).sensiden;
    float local_detailthr = (float)locallab.spots.at(sp).detailthr;

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

    lp.deltaem = locallab.spots.at(sp).deltae;
    lp.scalereti = scaleret;
    lp.cir = circr;
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
    lp.noiselequal = local_noiselequal;
    lp.noisechrodetail = local_noisechrodetail;
    lp.noiselc = local_noiselc;
    lp.noisecf = local_noisecf;
    lp.noisecc = local_noisecc;
    lp.sensden = local_sensiden;
    lp.bilat = locallab.spots.at(sp).bilateral;
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

    for (int y = 0; y < 6; y++) {
        lp.mulloc[y] = LIM(multi[y], 0.f, 4.f);//to prevent crash with old pp3 integer
    }

    for (int y = 0; y < 5; y++) {
        lp.mullocsh[y] = multish[y];
    }
    lp.activspot = locallab.spots.at(sp).activ;

    lp.logena = locallab.spots.at(sp).explog && lp.activspot;

    lp.detailsh = locallab.spots.at(sp).detailSH;
    lp.threshol = thresho;
    lp.chromacb = chromcbdl;
    lp.expvib = locallab.spots.at(sp).expvibrance && lp.activspot ;
    lp.colorena = locallab.spots.at(sp).expcolor && lp.activspot && llExpMask == 0 && llsoftMask == 0 && llSHMask == 0 && llcbMask == 0 && lllcMask == 0 && llsharMask == 0 && llretiMask == 0 && lltmMask == 0 && llvibMask == 0 && ll_Mask == 0; // Color & Light tool is deactivated if Exposure mask is visible or SHMask
    lp.blurena = locallab.spots.at(sp).expblur && lp.activspot && llExpMask == 0 && llsoftMask == 0 && llSHMask == 0 && llcbMask == 0 && lllcMask == 0 && llsharMask == 0 && llretiMask == 0 && llColorMask == 0 && lltmMask == 0 && llvibMask == 0 && ll_Mask == 0;
    lp.tonemapena = locallab.spots.at(sp).exptonemap && lp.activspot && llExpMask == 0 && llsoftMask == 0 && llSHMask == 0 && llcbMask == 0 && lllcMask == 0 && llsharMask == 0 && llretiMask == 0 && llColorMask == 0 && llvibMask == 0 && ll_Mask == 0;
    lp.retiena = locallab.spots.at(sp).expreti && lp.activspot && llExpMask == 0 && llsoftMask == 0 && llSHMask == 0 && llcbMask == 0 && lllcMask == 0 && llsharMask == 0 && llColorMask == 0 && lltmMask == 0 && llvibMask == 0 && llSHMask == 0 && ll_Mask == 0;
    lp.lcena = locallab.spots.at(sp).expcontrast && lp.activspot && llExpMask == 0 && llsoftMask == 0 && llSHMask == 0 && llcbMask == 0 && llsharMask == 0 && llColorMask == 0 && lltmMask == 0 && llvibMask == 0 && llSHMask == 0 && ll_Mask == 0;
    lp.cbdlena = locallab.spots.at(sp).expcbdl && lp.activspot && llExpMask == 0 && llsoftMask == 0 && llSHMask == 0 && llretiMask == 0 && lllcMask == 0 && llsharMask == 0 && lllcMask == 0 && llColorMask == 0 && lltmMask == 0 && llvibMask == 0 && ll_Mask == 0;
    lp.exposena = locallab.spots.at(sp).expexpose && lp.activspot  && llColorMask == 0 && llsoftMask == 0 && llSHMask == 0 && lllcMask == 0 && llsharMask == 0 && llcbMask == 0 && llretiMask == 0 && lltmMask == 0 && llvibMask == 0 && ll_Mask == 0; // Exposure tool is deactivated if Color & Light mask SHmask is visible
    lp.hsena = locallab.spots.at(sp).expshadhigh && lp.activspot && llColorMask == 0 && llsoftMask == 0 && llExpMask == 0 && llcbMask == 0 && lllcMask == 0 && llsharMask == 0 && llretiMask == 0 && lltmMask == 0 && llvibMask == 0 && ll_Mask == 0;// Shadow Highlight tool is deactivated if Color & Light mask or SHmask is visible
    lp.vibena = locallab.spots.at(sp).expvibrance && lp.activspot && llColorMask == 0 && llsoftMask == 0 && llExpMask == 0 && llcbMask == 0 && lllcMask == 0 && llsharMask == 0 && llretiMask == 0 && llcbMask == 0 && lltmMask == 0 && llSHMask == 0 && ll_Mask == 0;// vibrance tool is deactivated if Color & Light mask or SHmask is visible
    lp.sharpena = locallab.spots.at(sp).expsharp && lp.activspot && llColorMask == 0 && llsoftMask == 0 && llExpMask == 0 && llcbMask == 0 && lllcMask == 0 && llretiMask == 0 && llcbMask == 0 && lltmMask == 0 && llSHMask == 0 && llvibMask == 0 && ll_Mask == 0;
    lp.sfena = locallab.spots.at(sp).expsoft && lp.activspot && llColorMask == 0 && llExpMask == 0 && llcbMask == 0 && lllcMask == 0 && llretiMask == 0 && llcbMask == 0 && lltmMask == 0 && llSHMask == 0 && llvibMask == 0 && ll_Mask == 0;
    lp.maskena = locallab.spots.at(sp).expmask && lp.activspot && llColorMask == 0 && llsoftMask == 0 && llExpMask == 0 && llcbMask == 0 && lllcMask == 0 && llsharMask == 0 && llretiMask == 0 && llcbMask == 0 && lltmMask == 0 && llSHMask == 0;// vibrance tool is deactivated if Color & Light mask or SHmask is visible

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
    //increase sensitivity for low values
    float proexp = lp.expcomp;
    if (std::fabs(proexp) < 0.6f) {
        float interm = std::fabs(proexp) / 0.6f;
        interm = pow(interm, 3.f);
        lp.expcomp = proexp * interm;
    }
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
    lp.blwh = locallab.spots.at(sp).blwh;
    lp.senscolor = (int) locallab.spots.at(sp).colorscope;
    //replace scope color vibrance shadows
    lp.sens = lp.senscolor;
    lp.sensv = lp.senscolor;
    lp.senshs = lp.senscolor;

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
//J.Desmis 12 2019 - I will try to port a raw process in local adjustments
// I choose this one because, it is "new"
// Perhaps - probably no result, but perhaps ??

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
        return std::pow(x, source_gray) - 1 - target_gray * x + target_gray;
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


// basic log encoding taken from ACESutil.Lin_to_Log2, from
// https://github.com/ampas/aces-dev
// (as seen on pixls.us)
void ImProcFunctions::log_encode(Imagefloat *rgb, const struct local_params & lp, bool multiThread, int bfw, int bfh)
{
    /* J.Desmis 12 2019
        small adaptations to local adjustments
        replace log2 by log(lp.baselog) allows diferentiation between low and high lights
    */
    BENCHFUN
    const float gray = lp.sourcegray / 100.f;
    const float shadows_range = lp.blackev;
    const float dynamic_range = lp.whiteev - lp.blackev;
    const float noise = pow_F(2.f, -16.f);
    const float log2 = xlogf(lp.baselog);
    const float base = lp.targetgray > 1 && lp.targetgray < 100 && dynamic_range > 0 ? find_gray(std::abs(lp.blackev) / dynamic_range, lp.targetgray / 100.f) : 0.f;
    const float linbase = rtengine::max(base, 0.f);
    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);

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

    const auto norm =
    [&](float r, float g, float b) -> float {
        return Color::rgbLuminance(r, g, b, ws);

        // other possible alternatives (so far, luminance seems to work
        // fine though). See also
        // https://discuss.pixls.us/t/finding-a-norm-to-preserve-ratios-across-non-linear-operations
        //
        // MAX
        //return max(r, g, b);
        //
        // Euclidean
        //return std::sqrt(SQR(r) + SQR(g) + SQR(b));

        // weighted yellow power norm from https://youtu.be/Z0DS7cnAYPk
        // float rr = 1.22f * r / 65535.f;
        // float gg = 1.20f * g / 65535.f;
        // float bb = 0.58f * b / 65535.f;
        // float rr4 = SQR(rr) * SQR(rr);
        // float gg4 = SQR(gg) * SQR(gg);
        // float bb4 = SQR(bb) * SQR(bb);
        // float den = (rr4 + gg4 + bb4);
        // if (den > 0.f) {
        //     return 0.8374319f * ((rr4 * rr + gg4 * gg + bb4 * bb) / den) * 65535.f;
        // } else {
        //     return 0.f;
        // }
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
                float m = norm(r, g, b);

                if (m > noise) {
                    float mm = apply(m);
                    float f = mm / m;
                    r *= f;
                    b *= f;
                    g *= f;
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
                    Y2[y][x] = norm(rgb->r(y, x), rgb->g(y, x), rgb->b(y, x)) / 65535.f;
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

                if (t > noise && (t2 = norm(r, g, b)) > noise) {
                    float c = apply(t, false);
                    float f = c / t;
                    //   float t2 = norm(r, g, b);
                    float f2 = apply(t2) / t2;
                    f = intp(blend, f, f2);
                    assert(std::isfinite(f));
                    r *= f;
                    g *= f;
                    b *= f;
                    assert(std::isfinite(r));
                    assert(std::isfinite(g));
                    assert(std::isfinite(b));
                }
            }
        }
    }
}

void ImProcFunctions::getAutoLogloc(int sp, ImageSource *imgsrc, float *sourceg, float *blackev, float *whiteev, bool *Autogr, int fw, int fh, float xsta, float xend, float ysta, float yend, int SCALE)
{
    BENCHFUN
//adpatation to local adjustments Jacques Desmis 12 2019
    const PreviewProps pp(0, 0, fw, fh, SCALE);

    Imagefloat img(int(fw / SCALE + 0.5), int(fh / SCALE + 0.5));
    const ProcParams neutral;
    imgsrc->getImage(imgsrc->getWB(), TR_NONE, &img, pp, params->toneCurve, neutral.raw);
    imgsrc->convertColorSpace(&img, params->icm, imgsrc->getWB());
    float minVal = RT_INFINITY;
    float maxVal = -RT_INFINITY;
    const float ec = std::pow(2.f, params->toneCurve.expcomp);

    constexpr float noise = 1e-5;
    const int h = fh / SCALE;
    const int w = fw / SCALE;

    const int hsta = ysta * h;
    const int hend = yend * h;

    const int wsta = xsta * w;
    const int wend = xend * w;

    double mean = 0.0;
    int nc = 0;
    for (int y = hsta; y < hend; ++y) {
        for (int x = wsta; x < wend; ++x) {
            const float r = img.r(y, x), g = img.g(y, x), b = img.b(y, x);
            mean += static_cast<double>(0.2126f * Color::gamma_srgb(r) + 0.7152f * Color::gamma_srgb(g) + 0.0722f * Color::gamma_srgb(b));
            nc++;

            const float m = rtengine::max(0.f, r, g, b) / 65535.f * ec;
            if (m > noise) {
                const float l = rtengine::min(r, g, b) / 65535.f * ec;
                minVal = rtengine::min(minVal, l > noise ? l : m);
                maxVal = rtengine::max(maxVal, m);
            }
        }
    }

    //approximation sourcegray yb  source = 0.4 * yb

    if (maxVal > minVal) {
        const float log2 = std::log(2.f);
        const float dynamic_range = -xlogf(minVal / maxVal) / log2;

        if (settings->verbose) {
            std::cout << "AutoLog: min = " << minVal << ", max = " << maxVal
                      << ", DR = " << dynamic_range << std::endl;
        }

        if (Autogr[sp]) {
            double tot = 0.0;
            int n = 0;
            const float gmax = rtengine::min(maxVal / 2.f, 0.25f);
            const float gmin = rtengine::max(minVal * std::pow(2.f, rtengine::max((dynamic_range - 1.f) / 2.f, 1.f)), 0.05f);

            if (settings->verbose) {
                std::cout << "         gray boundaries: " << gmin << ", " << gmax << std::endl;
            }

            for (int y = ysta; y < yend; ++y) {
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
            } else {
                mean /= (nc * 65535.0);
                float yb;

                if (mean < 0.15) {
                    yb = 3.0f;
                } else if (mean < 0.3) {
                    yb = 5.0f;
                } else if (mean < 0.4) {
                    yb = 10.0f;
                } else if (mean < 0.45) {
                    yb = 15.0f;
                } else if (mean < 0.5) {
                    yb = 18.0f;
                } else if (mean < 0.55) {
                    yb = 23.0f;
                } else if (mean < 0.6) {
                    yb = 30.0f;
                } else {
                    yb = 45.f;
                }
                sourceg[sp] = 0.4f * yb;
                if (settings->verbose) {
                    std::cout << "         no samples found in range, resorting to Yb gray point value " << sourceg[sp]  << std::endl;
                }
            }
        }

        const float gray = sourceg[sp] / 100.f;
        whiteev[sp] = xlogf(maxVal / gray) / log2;
        blackev[sp] = whiteev[sp] - dynamic_range;
    }
}

void tone_eq(array2D<float> &R, array2D<float> &G, array2D<float> &B, const struct local_params & lp, const Glib::ustring &workingProfile, double scale, bool multithread)
// adapted from the tone equalizer of darktable
/*
    Copyright 2019 Alberto Griggio <alberto.griggio@gmail.com>
    Small adaptation to Local Adjustment 10 2019 Jacques Desmis <jdesmis@gmail.com>
    This file is part of darktable,
    copyright (c) 2018 Aurelien Pierre.

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/

{
    BENCHFUN

    const int W = R.getWidth();
    const int H = R.getHeight();
    array2D<float> Y(W, H);

    const auto log2 =
    [](float x) -> float {
        static const float l2 = xlogf(2);
        return xlogf(x) / l2;
    };

    const auto exp2 =
    [](float x) -> float {
        return pow_F(2.f, x);
    };
    // Build the luma channels: band-pass filters with gaussian windows of
    // std 2 EV, spaced by 2 EV
    const float centers[12] = {
        -18.0f, -16.0f, -14.0f, -12.0f, -10.0f, -8.0f, -6.0f,
        -4.0f, -2.0f, 0.0f, 2.0f, 4.0f
    };

    const auto conv = [&](int v, float lo, float hi) -> float {
        const float f = v < 0 ? lo : hi;
        return exp2(float(v) / 100.f * f);
    };
    const float factors[12] = {
        conv(lp.mullocsh[0], 2.f, 3.f), // -18 EV
        conv(lp.mullocsh[0], 2.f, 3.f), // -16 EV
        conv(lp.mullocsh[0], 2.f, 3.f), // -14 EV
        conv(lp.mullocsh[0], 2.f, 3.f), // -12 EV
        conv(lp.mullocsh[0], 2.f, 3.f), // -10 EV
        conv(lp.mullocsh[0], 2.f, 3.f), //  -8 EV
        conv(lp.mullocsh[1], 2.f, 3.f), //  -6 EV
        conv(lp.mullocsh[2], 2.5f, 2.5f), //  -4 EV
        conv(lp.mullocsh[3], 3.f, 2.f), //  -2 EV
        conv(lp.mullocsh[4], 3.f, 2.f), //   0 EV
        conv(lp.mullocsh[4], 3.f, 2.f), //   2 EV
        conv(lp.mullocsh[4], 3.f, 2.f)  //   4 EV
    };

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(workingProfile);

#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            Y[y][x] = Color::rgbLuminance(R[y][x], G[y][x], B[y][x], ws);
        }
    }

    int detail = LIM(lp.detailsh + 5, 0, 5);
    int radius = detail / scale + 0.5;
    float epsilon2 = 0.01f + 0.002f * rtengine::max(detail - 3, 0);

    if (radius > 0) {
        rtengine::guidedFilterLog(10.f, Y, radius, epsilon2, multithread);
    }

    if (lp.detailsh > 0) {
        array2D<float> Y2(W, H);
        constexpr float base_epsilon = 0.02f;
        constexpr float base_posterization = 5.f;

#ifdef _OPENMP
        #pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                float l = LIM(log2(rtengine::max(Y[y][x], 1e-9f)), centers[0], centers[11]);
                float ll = round(l * base_posterization) / base_posterization;
                Y2[y][x] = Y[y][x];
                Y[y][x] = exp2(ll);
            }
        }

        radius = 350.0 / scale;
        epsilon2 = base_epsilon / float(6 - rtengine::min(lp.detailsh, 5));
        rtengine::guidedFilter(Y2, Y, Y, radius, epsilon2, multithread);
    }

    const auto gauss =
    [](float b, float x) -> float {
        return xexpf((-SQR(x - b) / 4.0f));
    };

    // For every pixel luminance, the sum of the gaussian masks
    float w_sum = 0.f;

    for (int i = 0; i < 12; ++i) {
        w_sum += gauss(centers[i], 0.f);
    }

    const auto process_pixel =
    [&](float y) -> float {
        // convert to log space
        const float luma = rtengine::max(log2(rtengine::max(y, 0.f)), -18.0f);

        // build the correction as the sum of the contribution of each
        // luminance channel to current pixel
        float correction = 0.0f;

        for (int c = 0; c < 12; ++c)
        {
            correction += gauss(centers[c], luma) * factors[c];
        }

        correction /= w_sum;

        return correction;
    };

    LUTf lut(65536);

    for (int i = 0; i < 65536; ++i) {
        float y = float(i) / 65535.f;
        float c = process_pixel(y);
        lut[i] = c;
    }


#ifdef __SSE2__
    vfloat vfactors[12];
    vfloat vcenters[12];

    for (int i = 0; i < 12; ++i) {
        vfactors[i] = F2V(factors[i]);
        vcenters[i] = F2V(centers[i]);
    }

    const auto vgauss =
    [](vfloat b, vfloat x) -> vfloat {
        static const vfloat fourv = F2V(4.f);
        return xexpf((-SQR(x - b) / fourv));
    };

    vfloat zerov = F2V(0.f);
    vfloat vw_sum = F2V(w_sum);

    const vfloat noisev = F2V(-18.f);
    const vfloat xlog2v = F2V(xlogf(2.f));

    const auto vprocess_pixel =
    [&](vfloat y) -> vfloat {
        const vfloat luma = vmaxf(xlogf(vmaxf(y, zerov)) / xlog2v, noisev);

        vfloat correction = zerov;

        for (int c = 0; c < 12; ++c)
        {
            correction += vgauss(vcenters[c], luma) * vfactors[c];
        }

        correction /= vw_sum;

        return correction;
    };


    vfloat v1 = F2V(1.f);
    vfloat v65535 = F2V(65535.f);
#endif // __SSE2__


#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        int x = 0;


#ifdef __SSE2__

        for (; x < W - 3; x += 4) {
            vfloat cY = LVFU(Y[y][x]);
            vmask m = vmaskf_gt(cY, v1);
            vfloat corr;

            if (_mm_movemask_ps((vfloat)m)) {
                corr = vprocess_pixel(cY);
            } else {
                corr = lut[cY * v65535];
            }

            STVF(R[y][x], LVF(R[y][x]) * corr);
            STVF(G[y][x], LVF(G[y][x]) * corr);
            STVF(B[y][x], LVF(B[y][x]) * corr);
        }

#endif // __SSE2__

        for (; x < W; ++x) {
            float cY = Y[y][x];
            float corr = cY > 1.f ? process_pixel(cY) : lut[cY * 65535.f];
            R[y][x] *= corr;
            G[y][x] *= corr;
            B[y][x] *= corr;
        }
    }

}


void ImProcFunctions::ciecamloc_02float(int sp, LabImage* lab)
{
    //be careful quasi duplicate with branch cat02wb
    BENCHFUN

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

    int tempo;

    if (params->locallab.spots.at(sp).warm > 0) {
        tempo = 5000 - 30 * params->locallab.spots.at(sp).warm;
    } else {
        tempo = 5000 - 49 * params->locallab.spots.at(sp).warm;
    }

    ColorTemp::temp2mulxyz(params->wb.temperature, params->wb.method, Xw, Zw);  //compute white Xw Yw Zw  : white current WB
    ColorTemp::temp2mulxyz(tempo, "Custom", Xwout, Zwout);
    ColorTemp::temp2mulxyz(5000, "Custom", Xwsc, Zwsc);

    //viewing condition for surrsrc
    f  = 1.00f;
    c  = 0.69f;
    nc = 1.00f;
    //viewing condition for surround
    f2 = 1.0f, c2 = 0.69f, nc2 = 1.0f;
    //with which algorithm
    //  alg = 0;


    xwd = 100.0 * Xwout;
    zwd = 100.0 * Zwout;
    ywd = 100.f;

    xws = 100.0 * Xwsc;
    zws = 100.0 * Zwsc;
    yws = 100.f;


    yb2 = 18;
    //La and la2 = ambiant luminosity scene and viewing
    la = 400.f;
    const float la2 = 400.f;
    const float pilot = 2.f;
    const float pilotout = 2.f;

    //algoritm's params
    // const float rstprotection = 100. ;//- params->colorappearance.rstprotection;
    LUTu hist16J;
    LUTu hist16Q;
    float yb = 18.f;
    float d, dj;

    // const int gamu = 0; //(params->colorappearance.gamut) ? 1 : 0;
    xw = 100.0 * Xw;
    yw = 100.f * Yw;
    zw = 100.0 * Zw;
    float xw1 = xws, yw1 = yws, zw1 = zws, xw2 = xwd, yw2 = ywd, zw2 = zwd;

    float cz, wh, pfl;
    Ciecam02::initcam1float(yb, pilot, f, la, xw, yw, zw, n, d, nbb, ncb, cz, aw, wh, pfl, fl, c);
//   const float chr = 0.f;
    const float pow1 = pow_F(1.64f - pow_F(0.29f, n), 0.73f);
    float nj, nbbj, ncbj, czj, awj, flj;
    Ciecam02::initcam2float(yb2, pilotout, f2,  la2,  xw2,  yw2,  zw2, nj, dj, nbbj, ncbj, czj, awj, flj);
#ifdef __SSE2__
    const float reccmcz = 1.f / (c2 * czj);
#endif
    const float pow1n = pow_F(1.64f - pow_F(0.29f, nj), 0.73f);

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
            vfloat x, y, z;
            vfloat J, C, h, Q, M, s;

            vfloat c655d35 = F2V(655.35f);

            for (k = 0; k < width - 3; k += 4) {
                Color::Lab2XYZ(LVFU(lab->L[i][k]), LVFU(lab->a[i][k]), LVFU(lab->b[i][k]), x, y, z);
                x = x / c655d35;
                y = y / c655d35;
                z = z / c655d35;
                Ciecam02::xyz2jchqms_ciecam02float(J, C,  h,
                                                   Q,  M,  s, F2V(aw), F2V(fl), F2V(wh),
                                                   x,  y,  z,
                                                   F2V(xw1), F2V(yw1),  F2V(zw1),
                                                   F2V(c),  F2V(nc), F2V(pow1), F2V(nbb), F2V(ncb), F2V(pfl), F2V(cz), F2V(d));
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
                                                   c,  nc, pow1, nbb, ncb, pfl, cz, d);
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
                                                   c,  nc, pow1, nbb, ncb, pfl, cz, d);
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
                                                c2, nc2,  pow1n, nbbj, ncbj, flj, czj, dj, awj);
                x = xx * 655.35f;
                y = yy * 655.35f;
                z = zz * 655.35f;
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
                Ciecam02::jch2xyz_ciecam02float(x, y, z,
                                                LVF(Jbuffer[k]), LVF(Cbuffer[k]), LVF(hbuffer[k]),
                                                F2V(xw2), F2V(yw2), F2V(zw2),
                                                F2V(nc2), F2V(pow1n), F2V(nbbj), F2V(ncbj), F2V(flj), F2V(dj), F2V(awj), F2V(reccmcz));
                STVF(xbuffer[k], x * c655d35);
                STVF(ybuffer[k], y * c655d35);
                STVF(zbuffer[k], z * c655d35);
            }

            // XYZ2Lab uses a lookup table. The function behind that lut is a cube root.
            // SSE can't beat the speed of that lut, so it doesn't make sense to use SSE
            for (int j = 0; j < width; j++) {
                float Ll, aa, bb;
                //convert xyz=>lab
                Color::XYZ2Lab(xbuffer[j], ybuffer[j], zbuffer[j], Ll, aa, bb);

                lab->L[i][j] = Ll;
                lab->a[i][j] = aa;
                lab->b[i][j] = bb;
            }

#endif
        }

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
            const float blur = 10.f / sk * (thres + 0.8f * rad);

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

    double aepsil = (epsilmax - epsilmin) / 90.f;
    double bepsil = epsilmax - 100.f * aepsil;
    double epsil = aepsil * rad + bepsil;
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

void ImProcFunctions::exlabLocal(local_params& lp, int bfh, int bfw, int bfhr, int bfwr, LabImage* bufexporig, LabImage* lab, const LUTf& hltonecurve, const LUTf& shtonecurve, const LUTf& tonecurve, const float hueref, const float lumaref, const float chromaref)
{
    BENCHFUN
    //exposure local

    constexpr float maxran = 65536.f;
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

            constexpr float alap = 600.f;
            constexpr float blap = 100.f;
            constexpr float aa = (alap - blap) / 50.f;
            constexpr float bb = 100.f - 30.f * aa;

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




void ImProcFunctions::InverseBlurNoise_Local(LabImage * originalmask, float **bufchro, const struct local_params & lp,  const float hueref, const float chromaref, const float lumaref, LabImage * original, LabImage * transformed, const LabImage * const tmp1, int cx, int cy, int sk)
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
                } else /*if (lp.shapmet == 1)*/ {
                    calcTransitionrect(lox, loy, ach, lp, zone, localFactor);
                }

                const float clc = (previewbl) ? settings->previewselection * 100.f : bufchro[y][x];
                float abdelta2 = SQR(refa - maskptr->a[y][x]) + SQR(refb - maskptr->b[y][x]);
                float chrodelta2 = SQR(std::sqrt(SQR(maskptr->a[y][x]) + SQR(maskptr->b[y][x])) - (chromaref * 327.68f));
                float huedelta2 = abdelta2 - chrodelta2;

                float dE = std::sqrt(kab * (kch * chrodelta2 + kH * huedelta2) + kL * SQR(refL - maskptr->L[y][x]));
                const float reducdE = calcreducdE(dE, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, lp.sensbn);
                const float realstrchdE = reducdE * clc;

                switch (zone) {

                    case 0: { // outside selection and outside transition zone => full effect, no transition
                        float difL = tmp1->L[y][x] - original->L[y][x];
                        transformed->L[y][x] = CLIP(original->L[y][x] + difL * reducdE);
                        float difa = tmp1->a[y][x] - original->a[y][x];
                        float difb = tmp1->b[y][x] - original->b[y][x];
                        float flia = 1.f, flib = 1.f;
                        flia = flib = ((100.f + realstrchdE) / 100.f);
                        const float chra = tmp1->a[y][x];
                        const float chrb = tmp1->b[y][x];

                        if (!lp.actsp) {
                            difa = chra * flia - original->a[y][x];
                            difb = chrb * flib - original->b[y][x];
                            transformed->a[y][x] = clipC(original->a[y][x] + difa);
                            transformed->b[y][x] = clipC(original->b[y][x] + difb);
                        }

                        if (blshow) {
                            transformed->L[y][x] = CLIP(12000.f + difL);
                            transformed->a[y][x] = clipC(difa);
                            transformed->b[y][x] = clipC(difb);
                        } else if (previewbl || lp.prevdE) {
                            transformed->a[y][x] = 0.f;
                            transformed->b[y] [x] = (difb);
                        }

                        break;
                    }

                    case 1: { // inside transition zone
                        float difL = tmp1->L[y][x] - original->L[y][x];
                        float difa = tmp1->a[y][x] - original->a[y][x];
                        float difb = tmp1->b[y][x] - original->b[y][x];
                        float flia = 1.f, flib = 1.f;
                        flia = flib = ((100.f + realstrchdE) / 100.f);
                        const float chra = tmp1->a[y][x];
                        const float chrb = tmp1->b[y][x];

                        float factorx = 1.f - localFactor;
                        difL *= factorx;

                        transformed->L[y][x] = CLIP(original->L[y][x] + difL * reducdE);

                        if (!lp.actsp) {
                            difa = chra * flia - original->a[y][x];
                            difb = chrb * flib - original->b[y][x];
                            difa *= factorx;
                            difb *= factorx;
                            transformed->a[y][x] = clipC(original->a[y][x] + difa);
                            transformed->b[y][x] = clipC(original->b[y][x] + difb);

                        }

                        if (blshow) {
                            transformed->L[y][x] = CLIP(12000.f + difL);
                            transformed->a[y][x] = clipC(difa);
                            transformed->b[y][x] = clipC(difb);
                        } else if (previewbl) {
                            transformed->a[y][x] = 0.f;
                            transformed->b[y][x] = (difb);
                        }

                        break;
                    }

                    case 2: { // inside selection => no effect, keep original values
                        transformed->L[y][x] = original->L[y][x];

                        if (!lp.actsp) {

                            transformed->a[y][x] = original->a[y][x];
                            transformed->b[y][x] = original->b[y][x];
                        }
                    }
                }
            }
        }
    }
}

static void mean_fab(int xstart, int ystart, int bfw, int bfh, LabImage* bufexporig, const LabImage* original, float &fab, float &meanfab, float chrom, bool multiThread)
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
                bufexporig->a[y][x] = original->a[y + ystart][x + xstart];
                bufexporig->b[y][x] = original->b[y + ystart][x + xstart];
                sumab += std::fabs(bufexporig->a[y][x]);
                sumab += std::fabs(bufexporig->b[y][x]);
            }
        }

        meanfab = sumab / (2.f * nbfab);

        double som = 0.0;

#ifdef _OPENMP
        #pragma omp parallel for reduction(+:som) if(multiThread)
#endif
        for (int y = 0; y < bfh; y++) {
            for (int x = 0; x < bfw; x++) {
                som += SQR(std::fabs(bufexporig->a[y][x]) - meanfab) + SQR(std::fabs(bufexporig->b[y][x]) - meanfab);
            }
        }

        const float multsigma = (chrom >= 0.f ? 0.035f : 0.018f) * chrom + 1.f;

        const float stddv = std::sqrt(som / nbfab);
        fab = meanfab + multsigma * stddv;

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
    double gradient_angle = angs / 180.0 * rtengine::RT_PI;
    double varfeath = 0.01 * lp.feath;

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
    gp.ys = std::sqrt((float)h * h + (float)w * w) * (varfeath / cos(gradient_angle));
    gp.ys_inv = 1.0 / gp.ys;
    gp.top_edge_0 = gp.yc - gp.ys / 2.0;

    if (gp.ys < 1.0 / h) {
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
                    transformed->L[y + ystart][x + xstart] = (lumask * 400.f) + clipLoc(bufmaskorigSH->L[y][x]);
                    transformed->a[y + ystart][x + xstart] = bufexporig->a[y][x] * bufmaskorigSH->a[y][x];
                    transformed->b[y + ystart][x + xstart] = bufexporig->b[y][x] * bufmaskorigSH->b[y][x];
                }
            } else if (inv == 1) { //inverse
                if (zone == 0) {
                    transformed->L[y + ystart][x + xstart] = (lumask * 400.f) + clipLoc(bufmaskorigSH->L[y][x]);
                    transformed->a[y + ystart][x + xstart] = bufexporig->a[y][x] * bufmaskorigSH->a[y][x];
                    transformed->b[y + ystart][x + xstart] = bufexporig->b[y][x] * bufmaskorigSH->b[y][x];
                }
            }
        }
    }
}

void ImProcFunctions::discrete_laplacian_threshold(float * data_out, const float * data_in, size_t nx, size_t ny, float t)
{
    BENCHFUN

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
        table[i] = std::cos(pi_size * i);
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
    BENCHFUN

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
            data[i * nx + j] *= m2 / (2.f - cosx[j] - cosy[i]);
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
        mean += data[i];
        dt += SQR(data[i]);
    }

    mean /= size;
    dt /= size;
    dt -= SQR(mean);
    mean_p = mean;
    dt_p = std::sqrt(dt);
}

void ImProcFunctions::normalize_mean_dt(float * data, const float * ref, size_t size, float mod, float sigm)
{
    /*
     * Copyright 2009-2011 IPOL Image Processing On Line http://www.ipol.im/
     *

     * @file retinex_pde_lib.c discrete Poisson equation
     * @brief laplacian, DFT and Poisson routines
     *
     * @author Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
     * adapted for Rawtherapee - jacques Desmis july 2019
     */

    if (NULL == data || NULL == ref) {
        fprintf(stderr, "a pointer is NULL and should not be so\n");
        abort();
    }

    double mean_ref, mean_data, dt_ref, dt_data;

    /* compute mean and variance of the two arrays */
    mean_dt(ref, size, mean_ref, dt_ref);
    mean_dt(data, size, mean_data, dt_data);

    /* compute the normalization coefficients */
    const double a = dt_ref / dt_data;
    const double b = mean_ref - a * mean_data;

    const float modma = mod * a;
    const float sigmmmodmb = sigm * mod * b;
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

    BENCHFUN
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
        normalize_mean_dt(data_tmp, datain, bfw * bfh, 1.f, 1.f);
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
                                  bool multiThread, bool enaMask, bool showmaske, bool deltaE, bool modmask, bool zero, bool modif, float chrom, float rad, float lap, float gamma, float slope, float blendm, float blendmab, int shado, float amountcd, float anchorcd,
                                  const LUTf& lmasklocalcurve, bool localmaskutili,
                                  const LocwavCurve & loclmasCurvecolwav, bool lmasutilicolwav, int level_bl, int level_hl, int level_br, int level_hr,
                                  int shortcu, bool delt, const float hueref, const float chromaref, const float lumaref,
                                  float maxdE, float mindE, float maxdElim,  float mindElim, float iterat, float limscope, int scope,
                                  bool fftt, float blu_ma, float cont_ma, int indic
                                 )


{
    array2D<float> ble(bfw, bfh);
    array2D<float> blechro(bfw, bfh);
    array2D<float> hue(bfw, bfh);
    array2D<float> guid(bfw, bfh);
    const std::unique_ptr<LabImage> bufreserv(new LabImage(bfw, bfh));
    float meanfab, fab;
    mean_fab(xstart, ystart, bfw, bfh, bufcolorig, original, fab, meanfab, chrom, multiThread);
    float chromult = 1.f - 0.01f * chrom;
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
                if (lochhhmasCurve[i] != 0.5) {
                    HHmaskcurve = true;
                }
            }
        }

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
                        STVF(atan2Buffer[i], xatan2f(LVFU(bufcolorig->b[ir][i]), LVFU(bufcolorig->a[ir][i])));
                    }

                    for (; i < bfw; i++) {
                        atan2Buffer[i] = xatan2f(bufcolorig->b[ir][i], bufcolorig->a[ir][i]);
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
                        kmaskC = LIM01(kinv  - kneg * locccmasCurve[500.f * (0.0001f + std::sqrt(SQR(bufcolorig->a[ir][jr]) + SQR(bufcolorig->b[ir][jr])) / fab)]);
                    }

                    if (lochhmasCurve && lhmasutili) {
#ifdef __SSE2__
                        const float huema = atan2Buffer[jr];
#else
                        const float huema = xatan2f(bufcolorig->b[ir][jr], bufcolorig->a[ir][jr]);
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

                    /*
                    //keep here in case of...but !!
                                        if (lochhhmasCurve && HHmaskcurve) {

                    #ifdef __SSE2__
                                            huemah = atan2BufferH[jr];
                    #else
                                            huemah = xatan2f(bufcolorig->b[ir][jr], bufcolorig->a[ir][jr]);
                    #endif

                                            float hh = Color::huelab_to_huehsv2(huemah);
                                            hh += 1.f / 6.f;

                                            if (hh > 1.f) {
                                                hh -= 1.f;
                                            }

                                            const float val_HH = float (LIM01(((0.5f - lochhhmasCurve[500.f * hh]))));
                                            kmaskHH = 2.f * val_HH;
                                            const float hhro = kmaskHH;

                                            if (hhro != 0) {
                                                newhr = huemah + hhro;

                                                if (newhr > rtengine::RT_PI_F) {
                                                    newhr -= 2 * rtengine::RT_PI_F;
                                                } else if (newhr < -rtengine::RT_PI_F) {
                                                    newhr += 2 * rtengine::RT_PI_F;
                                                }
                                            }
                                            sincosval = xsincosf(newhr);

                                        }
                    */
                    bufmaskblurcol->L[ir][jr] = clipLoc(kmaskL + kmaskHL + kmasstru + kmasblur);
                    bufmaskblurcol->a[ir][jr] = clipC((kmaskC + chromult * kmaskH));
                    bufmaskblurcol->b[ir][jr] = clipC((kmaskC + chromult * kmaskH));

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

        std::unique_ptr<LabImage> bufprov;
        if (delt) {
            bufprov.reset(new LabImage(bfw, bfh));
            bufprov->CopyFrom(bufmaskblurcol, multiThread);
        }

        if (rad != 0.f) {
            const float tmpblur = rad < 0.f ? -1.f / rad : 1.f + rad;
            const int r1 = rtengine::max<int>(4 / sk * tmpblur + 0.5, 1);
            const int r2 = rtengine::max<int>(25 / sk * tmpblur + 0.5, 1);

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
                if (loclmasCurvecolwav[i] != 0.5) {
                    wavcurvemask = true;
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
            float** rdE = *(rdEBuffer.get());

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

    float ampli = 1.f + std::fabs(lp.colorde);
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
    BENCHFUN
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

    float ampli = 1.f + std::fabs(lp.colorde);
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

    BENCHFUN {
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



void ImProcFunctions::transit_shapedetect_retinex(int call, int senstype, LabImage * bufexporig, LabImage * bufmask, LabImage * buforigmas, float **buflight, float **bufchro, const float hueref, const float chromaref, const float lumaref, const struct local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy, int sk)
{

    BENCHFUN {
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

        const bool showmas = lp.showmaskretimet == 3 ;

        const std::unique_ptr<LabImage> origblur(new LabImage(GW, GH));
        const float radius = 3.f / sk;
        const bool usemaskreti = lp.enaretiMask && senstype == 4 && !lp.enaretiMasktmap;
        float strcli = 0.03f * lp.str;

        if (lp.scalereti == 1)
        {
            strcli = 0.015 * lp.str;
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
            const float previewint = settings->previewselection;

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

                    if (call == 2) {
                        cli = buflight[y - ystart][x - xstart];
                        clc = previewreti ? settings->previewselection * 100.f : bufchro[y - ystart][x - xstart];
                    } else {
                        cli = buflight[y][x];
                        clc = previewreti ? settings->previewselection * 100.f : bufchro[y][x];

                    }

                    const float reducdE = calcreducdE(dE, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, varsens) / 100.f;

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

                        if (previewreti) {
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

    BENCHFUN
    const int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
    const int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
    const int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
    const int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
    const int bfw = xend - xstart;
    const int bfh = yend - ystart;
    // printf("h=%f l=%f c=%f s=%f\n", hueref, lumaref, chromaref, sobelref);
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

        normalize_mean_dt(data, datain, bfh * bfw, 1.f, 1.f);
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

            if (tonecurv) { //Tone response curve  : does nothing if gamma=2.4 and slope=12.92 ==> gamma sRGB
                const float gamtone = params->locallab.spots.at(sp).gamSH;
                const float slotone = params->locallab.spots.at(sp).sloSH;
                cmsHTRANSFORM dummy = nullptr;
                workingtrc(tmpImage.get(), tmpImage.get(), GW, GH, -5, params->icm.workingProfile, 2.4, 12.92310, dummy, true, false, false);
                workingtrc(tmpImage.get(), tmpImage.get(), GW, GH, 5, params->icm.workingProfile, gamtone, slotone, dummy, false, true, true);
            }

            if (tonequ) {
                tmpImage->normalizeFloatTo1();
                array2D<float> Rtemp(GW, GH, tmpImage->r.ptrs, ARRAY2D_BYREFERENCE);
                array2D<float> Gtemp(GW, GH, tmpImage->g.ptrs, ARRAY2D_BYREFERENCE);
                array2D<float> Btemp(GW, GH, tmpImage->b.ptrs, ARRAY2D_BYREFERENCE);
                tone_eq(Rtemp, Gtemp, Btemp, lp, params->icm.workingProfile, sk, multiThread);
                tmpImage->normalizeFloatTo65535();
            }

            rgb2lab(*tmpImage, *temp, params->icm.workingProfile);
        }

    } else if (senstype == 1) { //exposure
        ImProcFunctions::exlabLocal(lp, GH, GW, GW, GH, original, temp.get(), hltonecurveloc, shtonecurveloc, tonecurveloc, hueref, lumaref, chromaref);

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
        //always calculate hueref, chromaref, lumaref  before others operations use in normal mode for all modules exceprt denoise
        struct local_params lp;
        calcLocalParams(sp, oW, oH, params->locallab, lp, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, locwavCurveden, locwavdenutili);
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
        //O.88623 = std::sqrt(PI / 4) ==> sqare equal to circle
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

        if ((lp.noiself > 0.f || lp.noiself0 > 0.f || lp.noiself2 > 0.f || lp.noiselc > 0.f || lp.noisecf > 0.f || lp.noisecc > 0.f) && lp.denoiena) {
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
                    aveLblur += blurorig->L[y][x];
                    aveAblur += blurorig->a[y][x];
                    aveBblur += blurorig->b[y][x];
                    aveChroblur += std::sqrt(SQR(blurorig->b[y - cy][x - cx]) + SQR(blurorig->a[y - cy][x - cx]));
                    nsb++;

                }
            }
        }

        //ref for luma, chroma, hue
        for (int y = rtengine::max(cy, (int)(lp.yc - spotSize)); y < rtengine::min(transformed->H + cy, (int)(lp.yc + spotSize + 1)); y++) {
            for (int x = rtengine::max(cx, (int)(lp.xc - spotSize)); x < rtengine::min(transformed->W + cx, (int)(lp.xc + spotSize + 1)); x++) {
                aveL += original->L[y - cy][x - cx];
                aveA += original->a[y - cy][x - cx];
                aveB += original->b[y - cy][x - cx];
                aveChro += std::sqrt(SQR(original->b[y - cy][x - cx]) + SQR(original->a[y - cy][x - cx]));
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
                avesobel += sobelL->L[y][x];
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
        aveChro /= 327.68f;
        avA = aveA / 327.68f;
        avB = aveB / 327.68f;
        avL = aveL / 327.68f;
        hueref = xatan2f(avB, avA);    //mean hue

        if (isdenoise) {
            aveLblur = aveLblur / nsb;
            aveChroblur = aveChroblur / nsb;
            aveChroblur /= 327.68f;
            aveAblur = aveAblur / nsb;
            aveBblur = aveBblur / nsb;
            float avAblur = aveAblur / 327.68f;
            float avBblur = aveBblur / 327.68f;
            float avLblur = aveLblur / 327.68f;
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

        if (lumaref > 95.f) {//to avoid crash
            lumaref = 95.f;
        }
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


void optfft(int N_fftwsize, int &bfh, int &bfw, int &bfhr, int &bfwr, struct local_params& lp, int H, int W, int &xstart, int &ystart, int &xend, int &yend, int cx, int cy)
{
    int ftsizeH = 1;
    int ftsizeW = 1;

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

    //optimize with size fftw
    bool reduW = false;
    bool reduH = false;

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
        bfhr = ftsizeH;
        reduH = true;
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
        bfwr = ftsizeW;
        reduW = true;
    }

    //new values optimized
    ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
    yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, H);
    xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
    xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, W);
    bfh = bfhr = yend - ystart;
    bfw = bfwr = xend - xstart;

    if (reduH) {
        bfhr = ftsizeH;
    }

    if (reduW) {
        bfwr = ftsizeW;
    }

    if (settings->verbose) {
        printf("Nyst=%i Nyen=%i lp.yc=%f lp.lyT=%f  lp.ly=%f bfh=%i bfhr=%i origH=%i ftsizeH=%i\n", ystart, yend, lp.yc, lp.lyT, lp.ly, bfh, bfhr, H, ftsizeH);
        printf("Nxst=%i Nxen=%i lp.xc=%f lp.lxL=%f  lp.lx=%f bfw=%i bfwr=%i origW=%i ftsizeW=%i\n", xstart, xend, lp.xc, lp.lxL, lp.lx, bfw, bfwr, W, ftsizeW);
    }
}

void ImProcFunctions::BlurNoise_Local(LabImage *tmp1, LabImage * originalmask, float **bufchro, const float hueref, const float chromaref, const float lumaref, local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy, int sk)
{
//local BLUR
    BENCHFUN

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

    const float ampli = 1.5f + 0.5f * std::fabs(lp.colorde);

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
                const float clc = previewbl ? settings->previewselection * 100.f : bufchro[y - ystart][x - xstart];
                const float realstrchdE = reducdE * clc;

                float difL = (tmp1->L[y - ystart][x - xstart] - original->L[y][x]) * localFactor * reducdE;
                transformed->L[y][x] = CLIP(original->L[y][x] + difL);
                const float fli = (100.f + realstrchdE) / 100.f;
                const float difa = tmp1->a[y - ystart][x - xstart] * fli - original->a[y][x] * localFactor;
                const float difb = tmp1->b[y - ystart][x - xstart] * fli - original->b[y][x] * localFactor;

                if (!lp.actsp) {
                    transformed->a[y][x] = clipC(original->a[y][x] + difa);
                    transformed->b[y][x] = clipC(original->b[y][x] + difb);
                }

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

void ImProcFunctions::transit_shapedetect2(int call, int senstype, const LabImage * bufexporig, const LabImage * bufexpfin, LabImage * originalmask, const float hueref, const float chromaref, const float lumaref, float sobelref, float meansobel, float ** blend2, struct local_params & lp, LabImage * original, LabImage * transformed, int cx, int cy, int sk)
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
        optfft(N_fftwsize, bfh, bfw, bfhr, bfwr, lp, original->H, original->W, xstart, ystart, xend, yend, cx, cy);
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
 
    const bool masshow = ((lp.showmask_met == 1) &&  senstype == 20);

    const bool previewvib = ((lp.showmaskvibmet == 4) &&  senstype == 2);
    const bool previewexp = ((lp.showmaskexpmet == 5) &&  senstype == 1);
    const bool previewcol = ((lp.showmaskcolmet == 5) &&  senstype == 0);
    const bool previewSH = ((lp.showmaskSHmet == 4) &&  senstype == 9);
    const bool previewtm = ((lp.showmasktmmet == 4) &&  senstype == 8);
    const bool previewlc = ((lp.showmasklcmet == 4) &&  senstype == 10);
    const bool previeworig = ((lp.showmasksoftmet == 6) &&  senstype == 3 && lp.softmet == 1);
    const bool previewmas = ((lp.showmask_met == 3) &&  senstype == 20);

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
    const float kL = lp.balance / SQR(327.68f);
    const float kab = balancedeltaE(lp.balance) / SQR(327.68f);
    const float kH = lp.balanceh;
    const float kch = balancedeltaE(kH);

    if (lp.colorde == 0) {
        lp.colorde = -1;//to avoid black
    }

    float ampli = 1.f + std::fabs(lp.colorde);
    ampli = 2.f + 0.5f * (ampli - 2.f);

    float darklim = 5000.f;
    float aadark = -1.f;
    float bbdark = darklim;

    const bool usemaskvib = (lp.showmaskvibmet == 2 || lp.enavibMask || lp.showmaskvibmet == 4) && senstype == 2;
    const bool usemaskexp = (lp.showmaskexpmet == 2 || lp.enaExpMask || lp.showmaskexpmet == 5) && senstype == 1;
    const bool usemaskcol = (lp.showmaskcolmet == 2 || lp.enaColorMask || lp.showmaskcolmet == 5) && senstype == 0;
    const bool usemaskSH = (lp.showmaskSHmet == 2 || lp.enaSHMask || lp.showmaskSHmet == 4) && senstype == 9;
    const bool usemasktm = (lp.showmasktmmet == 2 || lp.enatmMask || lp.showmasktmmet == 4) && senstype == 8;
    const bool usemasklc = (lp.showmasklcmet == 2 || lp.enalcMask || lp.showmasklcmet == 4) && senstype == 10;
    const bool usemaskmas = (lp.showmask_met == 1 || lp.ena_Mask || lp.showmask_met == 3) && senstype == 20;
    const bool usemaskall = (usemaskexp || usemaskvib || usemaskcol || usemaskSH || usemasktm || usemasklc || usemaskmas);

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

        normalize_mean_dt(data, datain, bfh * bfw, 1.f, 1.f);
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

                    if ((expshow || vibshow || colshow || SHshow || tmshow || lcshow || origshow || masshow) && lp.colorde < 0) { //show modifications with use "b"
                        //  (origshow && lp.colorde < 0) { //original Retinex
                        transformed->a[y + ystart][x + xstart] = 0.f;
                        transformed->b[y + ystart][x + xstart] = ampli * 8.f * diflc * reducdE;
                        transformed->L[y + ystart][x + xstart] = CLIP(12000.f + 0.5f * ampli * diflc);

                    } else if ((expshow || vibshow || colshow || SHshow || tmshow || lcshow || origshow || masshow) && lp.colorde > 0) {//show modifications without use "b"
                        if (diflc < 1000.f) {//if too low to be view use ab
                            diflc += 0.5f * maxdifab;
                        }

                        transformed->L[y + ystart][x + xstart] = CLIP(12000.f + 0.5f * ampli * diflc);
                        transformed->a[y + ystart][x + xstart] = clipC(ampli * difa);
                        transformed->b[y + ystart][x + xstart] = clipC(ampli * difb);
                    } else if (previewexp || previewvib || previewcol || previewSH || previewtm || previewlc || previeworig || previewmas || lp.prevdE) {//show deltaE
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

    BENCHFUN
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

    normalize_mean_dt(data, dataor, bfw * bfh, mod, 1.f);
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
    BENCHFUN

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
        radsig = 1.f / (2.f * rtengine::RT_PI * radius * radius);//gauss
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
    BENCHFUN
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
            const int numblox_W = ceil((static_cast<float>(GW)) / offset);
            const float DCTnorm = 1.0f / (4 * tilssize * tilssize); //for DCT

            int imin = rtengine::max(0, - topproc);
            int bottom = rtengine::min(topproc + tilssize, GH);
            int imax = bottom - topproc;

            for (int i = imin; i < imax; ++i) {
                for (int hblk = 0; hblk < numblox_W; ++hblk) {
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

void ImProcFunctions::wavcont(const struct local_params& lp, float ** tmp, wavelet_decomposition& wdspot, int level_bl, int maxlvl,
                              const LocwavCurve & loclevwavCurve, bool loclevwavutili,
                              const LocwavCurve & loccompwavCurve, bool loccompwavutili,
                              const LocwavCurve & loccomprewavCurve, bool loccomprewavutili,
                              float radlevblur, int process, float chromablu, float thres,  float sigmadc, float deltad)
{
    BENCHFUN
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
                const int W_L = wdspot.level_W(level);
                const int H_L = wdspot.level_H(level);
                const auto wav_L = wdspot.level_coeffs(level)[dir];
                madL[level][dir - 1] = Mad(wav_L, W_L * H_L);//evaluate noise by level
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
BENCHFUN
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

    int W_L = wdspot->level_W(0);
    int H_L = wdspot->level_H(0);
    float *wav_L0 = wdspot->get_coeff0();

    if (radblur > 0.f && blurena) {
        float* src[H_L];
        for (int i = 0; i < H_L; ++i) {
            src[i] = &wav_L0[i * W_L];
        }

#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            gaussianBlur(src, src, W_L, H_L, radblur);
        }
    }

    if (compress != 0.f && compreena) {
        const float Compression = expf(-compress);
        const float DetailBoost = std::max(compress, 0.f);

        CompressDR(wav_L0, W_L, H_L, Compression, DetailBoost);
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
        for (int i = 0; i < W_L * H_L; i++) {
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
        const std::unique_ptr<LabImage> temp(new LabImage(W_L, H_L));
#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif

        for (int i = 0; i < H_L; i++) {
            for (int j = 0; j < W_L; j++) {
                temp->L[i][j] = wav_L0[i * W_L + j];
            }
        }

        ImProcFunctions::shadowsHighlights(temp.get(), true, 1, lp.residhi, lp.residsha , 40, sk, lp.residhithr, lp.residshathr);

#ifdef _OPENMP
        #pragma omp parallel for if (multiThread)
#endif

        for (int i = 0; i < H_L; i++) {
            for (int j = 0; j < W_L; j++) {
                wav_L0[i * W_L + j] = temp->L[i][j];
            }
        }
    }

    if (contrast != 0.) {
        double avedbl = 0.0; // use double precision for large summations

#ifdef _OPENMP
        #pragma omp parallel for reduction(+:avedbl) if (multiThread)
#endif
        for (int i = 0; i < W_L * H_L; i++) {
            avedbl += wav_L0[i];
        }

        const float avg = LIM01(avedbl / (32768.f * W_L * H_L));
        double contreal = 0.6 * contrast;
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
        for (int i = 0; i < W_L * H_L; i++) {
            wav_L0[i] = resid_contrast.getVal(LIM01(wav_L0[i] / 32768.f)) * 32768.f;
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
        int W_L = wdspot->level_W(0);//provisory W_L H_L
        int H_L = wdspot->level_H(0);

        float *koeLi[12];
        float maxkoeLi[12] = {0.f};

        float *koeLibuffer = new float[12 * H_L * W_L]; //12

        for (int i = 0; i < 12; i++) {
            koeLi[i] = &koeLibuffer[i * W_L * H_L];
        }

        array2D<float> tmC(W_L, H_L);

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
            for (int i = 1; i < H_L - 1; i++) {
                for (int j = 1; j < W_L - 1; j++) {
                    //treatment of koeLi and maxkoeLi
                    if (lp.lip3) {//Sobel Canny algo improve with parameters
                        // comparison between pixel and neighbors
                        const auto neigh = lp.neiwmet == 1;
                        const auto kneigh = neigh ? 28.f : 38.f;
                        const auto somm = neigh ? 40.f : 50.f;

                        for (int dir = 1; dir < 4; dir++) { //neighbors proxi
                            koeLi[lvl * 3 + dir - 1][i * W_L + j] = (kneigh * koeLi[lvl * 3 + dir - 1][i * W_L + j] + 
                                                                    2.f * koeLi[lvl * 3 + dir - 1][(i - 1) * W_L + j] + 2.f * koeLi[lvl * 3 + dir - 1][(i + 1) * W_L + j] + 2.f * koeLi[lvl * 3 + dir - 1][i * W_L + j + 1] + 2.f * koeLi[lvl * 3 + dir - 1][i * W_L + j - 1]
                                                                    + koeLi[lvl * 3 + dir - 1][(i - 1) * W_L + j - 1] + koeLi[lvl * 3 + dir - 1][(i - 1) * W_L + j + 1] + koeLi[lvl * 3 + dir - 1][(i + 1) * W_L + j - 1] + koeLi[lvl * 3 + dir - 1][(i + 1) * W_L + j + 1]) / somm;
                        }
                    }

                    float interm = 0.f;
                    for (int dir = 1; dir < 4; dir++) {
                        //here I evaluate combination of vert / diag / horiz...we are with multiplicators of the signal
                        interm += SQR(koeLi[lvl * 3 + dir - 1][i * W_L + j]);
                    }

                    interm = std::sqrt(interm) * 0.57736721f;

                    constexpr float eps = 0.0001f;
                    // I think this double ratio (alph, beta) is better than arctg

                    float alph = koeLi[lvl * 3][i * W_L + j] / (koeLi[lvl * 3 + 1][i * W_L + j] + eps); //ratio between horizontal and vertical
                    float beta = koeLi[lvl * 3 + 2][i * W_L + j] / (koeLi[lvl * 3 + 1][i * W_L + j] + eps); //ratio between diagonal and horizontal

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
                    koeLi[lvl * 3][i * W_L + j] = koeLi[lvl * 3 + 1][i * W_L + j] = koeLi[lvl * 3 + 2][i * W_L + j] = interm; //new value
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
        float mean[10];
        float meanN[10];
        float sigma[10];
        float sigmaN[10];
        float MaxP[10];
        float MaxN[10];
        Evaluate2(*wdspot, mean, meanN, sigma, sigmaN, MaxP, MaxN, numThreads);
        for (int dir = 1; dir < 4; dir++) {
            for (int level = level_bl; level < maxlvl; ++level) {
                int W_L = wdspot->level_W(level);
                int H_L = wdspot->level_H(level);
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
                float* const* wav_L = wdspot->level_coeffs(level);

                if (MaxP[level] > 0.f && mean[level] != 0.f && sigma[level] != 0.f) {
                    constexpr float insigma = 0.666f; //SD
                    const float logmax = log(MaxP[level]); //log Max
                    const float rapX = (mean[level] + lp.sigmalc * sigma[level]) / MaxP[level]; //rapport between sD / max
                    const float inx = log(insigma);
                    const float iny = log(rapX);
                    const float rap = inx / iny; //koef
                    const float asig = 0.166f / (sigma[level] * lp.sigmalc);
                    const float bsig = 0.5f - asig * mean[level];
                    const float amean = 0.5f / mean[level];
                    const float limit1 = mean[level] + lp.sigmalc * sigma[level];
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
                        const float reduceeffect = kc <= 0.f ? 1.f : 1.5f;

                        float kinterm = 1.f + reduceeffect * kc;
                        kinterm = kinterm <= 0.f ? 0.01f : kinterm;

                        wav_L[dir][i] *= kinterm <= 0.f ? 0.01f : kinterm;
                    }
                }
            }
        }
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
}


void ImProcFunctions::fftw_denoise(int GW, int GH, int max_numblox_W, int min_numblox_W, float **tmp1, array2D<float> *Lin, int numThreads, const struct local_params & lp, int chrom)
{
    BENCHFUN

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
                noisevar_Ldetail = SQR(static_cast<float>(SQR(100. - params_Ldetail) + 50.*(100. - params_Ldetail)) * TS * 0.5f);
            } else if (chrom == 1) {
                params_Ldetail = rtengine::min(float(lp.noisechrodetail), 99.9f);
                //   noisevar_Ldetail = 100.f * pow((static_cast<float>(SQR(100. - params_Ldetail) + 50.*(100. - params_Ldetail)) * TS * 0.5f), 2);//to test ???
                noisevar_Ldetail = 100.f * pow((static_cast<float>(SQR(100. - params_Ldetail)) * TS * 0.5f), 2);//to test ???
            }

            //   float noisevar_Ldetail = SQR(static_cast<float>(SQR(100. - params_Ldetail) + 50.*(100. - params_Ldetail)) * TS * 0.5f);



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

    //Threshold DCT from Alberto Grigio
    const int detail_thresh = lp.detailthr;
    array2D<float> mask;

    if (detail_thresh > 0) {
        mask(GW, GH);
        float thr = log2lin(float(detail_thresh) / 200.f, 100.f);
        buildBlendMask(prov, mask, GW, GH, thr);
#ifdef _OPENMP
        #pragma omp parallel if (multiThread)
#endif
        {
            gaussianBlur(mask, mask, GW, GH, 20.0);
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

void ImProcFunctions::DeNoise(int call, int del, float * slidL, float * slida, float * slidb, int aut,  bool noiscfactiv, const struct local_params & lp, LabImage * originalmaskbl, int levred, float huerefblur, float lumarefblur, float chromarefblur, LabImage * original, LabImage * transformed, int cx, int cy, int sk)
{

//local denoise
    //all these variables are to prevent use of denoise when non necessary
    // but with qualmet = 2 (default for best quality) we must denoise chroma with little values to prevent artifacts due to variations of Hue
    // but if user select voluntary denoise, it is that choice the good (prioritary)
    bool execcolor = (lp.chro != 0.f || lp.ligh != 0.f || lp.cont != 0); // only if one slider or more is engaged
    bool execbdl = (lp.mulloc[0] != 1.f || lp.mulloc[1] != 1.f || lp.mulloc[2] != 1.f || lp.mulloc[3] != 1.f || lp.mulloc[4] != 1.f || lp.mulloc[5] != 1.f) ;//only if user want cbdl
    bool execdenoi = noiscfactiv && ((lp.colorena && execcolor) || (lp.tonemapena && lp.strengt != 0.f) || (lp.cbdlena && execbdl) || (lp.sfena && lp.strng > 0.f) || (lp.lcena && lp.lcamount > 0.f) || (lp.sharpena && lp.shrad > 0.42) || (lp.retiena && lp.str > 0.f)  || (lp.exposena && lp.expcomp != 0.f)  || (lp.expvib && lp.past != 0.f));
    bool execmaskden = (lp.showmaskblmet == 2 || lp.enablMask || lp.showmaskblmet == 3 || lp.showmaskblmet == 4) && lp.smasktyp != 0;

    if (((lp.noiself > 0.f || lp.noiself0 > 0.f || lp.noiself2 > 0.f || lp.noiselc > 0.f || lp.noisecf > 0.f || lp.noisecc > 0.f
//            || lp.showmaskblmet == 2 || lp.enablMask || lp.showmaskblmet == 3 || lp.showmaskblmet == 4  || aut == 1 || aut == 2) && lp.denoiena) || execdenoi) {  // sk == 1 ??
            || execmaskden || aut == 1 || aut == 2) && lp.denoiena) || execdenoi) {  // sk == 1 ??

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


#ifdef _OPENMP
        const int numThreads = omp_get_max_threads();
#else
        const int numThreads = 1;

#endif

        if (call == 1 && GW >= mDEN && GH >= mDEN) {


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
                //      float mxsfl = 0.f;

                if (aut == 0) {
                    if (levred == 7) {
                        edge = 2;
                        vari[0] = 0.8f * SQR((lp.noiself0 / 125.0) * (1.0 + lp.noiself0 / 25.0));
                        vari[1] = 0.8f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[2] = 0.8f * SQR((lp.noiself2 / 125.0) * (1.0 + lp.noiself2 / 25.0));

                        vari[3] = 0.8f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[4] = 0.8f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[5] = 0.8f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[6] = 0.8f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                    } else if (levred == 4) {
                        edge = 3;
                        vari[0] = 0.8f * SQR((lp.noiself0 / 125.0) * (1.0 + lp.noiself0 / 25.0));
                        vari[1] = 0.8f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                        vari[2] = 0.8f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        vari[3] = 0.8f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));

                    }
                } else if (aut == 1  || aut == 2) {
                    edge = 2;
                    vari[0] = SQR(slidL[0]);
                    vari[1] = SQR(slidL[1]);
                    vari[2] = SQR(slidL[2]);
                    vari[3] = SQR(slidL[3]);
                    vari[4] = SQR(slidL[4]);
                    vari[5] = SQR(slidL[5]);
                    vari[6] = SQR(slidL[6]);
                    float mxslid34 = rtengine::max(slidL[3], slidL[4]);
                    float mxslid56 = rtengine::max(slidL[5], slidL[6]);
                    mxsl = rtengine::max(mxslid34, mxslid56);

                }

                {
                    float kr3 = 0.f;
                    float kr4 = 0.f;
                    float kr5 = 0.f;

                    if (aut == 0 || aut == 1) {
                        if ((lp.noiselc < 30.f && aut == 0) || (mxsl < 30.f && aut == 1)) {
                            kr3 = 0.f;
                            kr4 = 0.f;
                            kr5 = 0.f;
                        } else if ((lp.noiselc < 50.f && aut == 0) || (mxsl < 50.f && aut == 1)) {
                            kr3 = 0.5f;
                            kr4 = 0.3f;
                            kr5 = 0.2f;
                        } else if ((lp.noiselc < 70.f && aut == 0) || (mxsl < 70.f && aut == 1)) {
                            kr3 = 0.7f;
                            kr4 = 0.5f;
                            kr5 = 0.3f;
                        } else {
                            kr3 = 1.f;
                            kr4 = 1.f;
                            kr5 = 1.f;
                        }
                    } else if (aut == 2) {
                        kr3 = 1.f;
                        kr4 = 1.f;
                        kr5 = 1.f;
                    }

                    vari[0] = rtengine::max(0.000001f, vari[0]);
                    vari[1] = rtengine::max(0.000001f, vari[1]);
                    vari[2] = rtengine::max(0.000001f, vari[2]);
                    vari[3] = rtengine::max(0.000001f, kr3 * vari[3]);

                    if (levred == 7) {
                        vari[4] = rtengine::max(0.000001f, kr4 * vari[4]);
                        vari[5] = rtengine::max(0.000001f, kr5 * vari[5]);
                        vari[6] = rtengine::max(0.000001f, kr5 * vari[6]);
                    }

                    float* noisevarlum = new float[GH * GW];
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

                    if ((lp.noiselc < 0.02f && aut == 0) || (mxsl < 1.f && (aut == 1 || aut == 2))) {
                        WaveletDenoiseAllL(Ldecomp, noisevarlum, madL, vari, edge, numThreads);

                    } else {

                        WaveletDenoiseAll_BiShrinkL(Ldecomp, noisevarlum, madL, vari, edge, numThreads);

                        WaveletDenoiseAllL(Ldecomp, noisevarlum, madL, vari, edge, numThreads);

                    }

                    delete[] noisevarlum;

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
                float maxcfine = 0.f;
                float maxccoarse = 0.f;

                if (aut == 0) {
                    if (levred == 7) {
                        edge = 2;
                        variC[0] = SQR(noisecfr);
                        variC[1] = SQR(noisecfr);
                        variC[2] = SQR(noisecfr);

                        variC[3] = SQR(noisecfr);
                        variC[4] = SQR(noisecfr);
                        variC[5] = SQR(noiseccr);
                        variC[6] = SQR(noiseccr);

                        variCb[0] = SQR(noisecfb);
                        variCb[1] = SQR(noisecfb);
                        variCb[2] = SQR(noisecfb);

                        variCb[3] = SQR(noisecfb);
                        variCb[4] = SQR(noisecfb);
                        variCb[5] = SQR(noiseccb);
                        variCb[6] = SQR(noiseccb);

                    } else if (levred == 4) {
                        edge = 3;
                        variC[0] = SQR(lp.noisecf / 10.0);
                        variC[1] = SQR(lp.noisecf / 10.0);
                        variC[2] = SQR(lp.noisecf / 10.0);
                        variC[3] = SQR(lp.noisecf / 10.0);

                        variCb[0] = SQR(lp.noisecf / 10.0);
                        variCb[1] = SQR(lp.noisecf / 10.0);
                        variCb[2] = SQR(lp.noisecf / 10.0);
                        variCb[3] = SQR(lp.noisecf / 10.0);

                    }
                } else if (aut == 1  || aut == 2) {
                    edge = 2;
                    variC[0] = SQR(slida[0]);
                    variC[1] = SQR(slida[1]);
                    variC[2] = SQR(slida[2]);
                    variC[3] = SQR(slida[3]);
                    variC[4] = SQR(slida[4]);
                    variC[5] = SQR(slida[5]);
                    variC[6] = SQR(slida[6]);
                    float maxc01 = rtengine::max(slida[0], slida[1]);
                    float maxc23 = rtengine::max(slida[2], slida[3]);
                    float max03 = rtengine::max(maxc01, maxc23);
                    float maxrf = rtengine::max(max03, slida[4]);
                    float maxrc = rtengine::max(slida[5], slida[6]);

                    variCb[0] = SQR(slidb[0]);
                    variCb[1] = SQR(slidb[1]);
                    variCb[2] = SQR(slidb[2]);
                    variCb[3] = SQR(slidb[3]);
                    variCb[4] = SQR(slidb[4]);
                    variCb[5] = SQR(slidb[5]);
                    variCb[6] = SQR(slidb[6]);
                    float maxb01 = rtengine::max(slidb[0], slidb[1]);
                    float maxb23 = rtengine::max(slidb[2], slidb[3]);
                    float maxb03 = rtengine::max(maxb01, maxb23);
                    float maxbf = rtengine::max(maxb03, slidb[4]);
                    maxcfine = rtengine::max(maxrf, maxbf);

                    float maxbc = rtengine::max(slidb[5], slidb[6]);
                    maxccoarse = rtengine::max(maxrc, maxbc);

                }

                {
                    float minic = 0.000001f;

                    if (noiscfactiv) {
                        minic = 0.1f;//only for artifact shape detection
                    }

                    float k1 = 0.f;
                    float k2 = 0.f;
                    float k3 = 0.f;

                    if (aut == 0 || aut == 1) {
                        if ((lp.noisecf < 0.2f && aut == 0) || (maxcfine < 0.2f && aut == 1)) {
                            k1 = 0.05f;
                            k2 = 0.f;
                            k3 = 0.f;
                        } else if ((lp.noisecf < 0.3f && aut == 0) || (maxcfine < 0.3f && aut == 1)) {
                            k1 = 0.1f;
                            k2 = 0.0f;
                            k3 = 0.f;
                        } else if ((lp.noisecf < 0.5f && aut == 0) || (maxcfine < 0.5f && aut == 1)) {
                            k1 = 0.2f;
                            k2 = 0.1f;
                            k3 = 0.f;
                        } else if ((lp.noisecf < 0.8f && aut == 0) || (maxcfine < 0.8f && aut == 1)) {
                            k1 = 0.3f;
                            k2 = 0.25f;
                            k3 = 0.f;
                        } else if ((lp.noisecf < 1.f && aut == 0) || (maxcfine < 1.f && aut == 1)) {
                            k1 = 0.4f;
                            k2 = 0.25f;
                            k3 = 0.1f;
                        } else if ((lp.noisecf < 2.f && aut == 0) || (maxcfine < 2.f && aut == 1)) {
                            k1 = 0.5f;
                            k2 = 0.3f;
                            k3 = 0.15f;
                        } else if ((lp.noisecf < 3.f && aut == 0) || (maxcfine < 3.f && aut == 1)) {
                            k1 = 0.6f;
                            k2 = 0.45f;
                            k3 = 0.3f;
                        } else if ((lp.noisecf < 4.f && aut == 0) || (maxcfine < 4.f && aut == 1)) {
                            k1 = 0.7f;
                            k2 = 0.5f;
                            k3 = 0.4f;
                        } else if ((lp.noisecf < 5.f && aut == 0) || (maxcfine < 5.f && aut == 1)) {
                            k1 = 0.8f;
                            k2 = 0.6f;
                            k3 = 0.5f;
                        } else if ((lp.noisecf < 6.f && aut == 0) || (maxcfine < 10.f && aut == 1)) {
                            k1 = 0.85f;
                            k2 = 0.7f;
                            k3 = 0.6f;
                        } else if ((lp.noisecf < 8.f && aut == 0) || (maxcfine < 20.f && aut == 1)) {
                            k1 = 0.9f;
                            k2 = 0.8f;
                            k3 = 0.7f;
                        } else if ((lp.noisecf < 10.f && aut == 0) || (maxcfine < 50.f && aut == 1)) {
                            k1 = 1.f;
                            k2 = 1.f;
                            k3 = 0.9f;

                        } else {
                            k1 = 1.f;
                            k2 = 1.f;
                            k3 = 1.f;
                        }
                    } else if (aut == 2) {
                        k1 = 1.f;
                        k2 = 1.f;
                        k3 = 1.f;
                    }


                    variC[0] = rtengine::max(minic, variC[0]);
                    variC[1] = rtengine::max(minic, k1 * variC[1]);
                    variC[2] = rtengine::max(minic, k2 * variC[2]);
                    variC[3] = rtengine::max(minic, k3 * variC[3]);

                    variCb[0] = rtengine::max(minic, variCb[0]);
                    variCb[1] = rtengine::max(minic, k1 * variCb[1]);
                    variCb[2] = rtengine::max(minic, k2 * variCb[2]);
                    variCb[3] = rtengine::max(minic, k3 * variCb[3]);

                    if (levred == 7) {
                        float k4 = 0.f;
                        float k5 = 0.f;
                        float k6 = 0.f;

                        if ((lp.noisecc < 0.2f && aut == 0) || (maxccoarse < 0.2f && aut == 1)) {
                            k4 = 0.1f;
                            k5 = 0.02f;
                        } else if ((lp.noisecc < 0.5f && aut == 0) || (maxccoarse < 0.5f && aut == 1)) {
                            k4 = 0.15f;
                            k5 = 0.05f;
                        } else if ((lp.noisecc < 1.f && aut == 0) || (maxccoarse < 1.f && aut == 1)) {
                            k4 = 0.15f;
                            k5 = 0.1f;
                        } else if ((lp.noisecc < 3.f && aut == 0) || (maxccoarse < 3.f && aut == 1)) {
                            k4 = 0.3f;
                            k5 = 0.15f;
                        } else if ((lp.noisecc < 4.f && aut == 0) || (maxccoarse < 5.f && aut == 1)) {
                            k4 = 0.6f;
                            k5 = 0.4f;
                        } else if ((lp.noisecc < 6.f && aut == 0) || (maxccoarse < 6.f && aut == 1)) {
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

                        if ((lp.noisecc < 4.f && aut == 0) || (maxccoarse < 4.f && aut == 1)) {
                            k6 = 0.f;
                        } else if ((lp.noisecc < 5.f && aut == 0) || (maxccoarse < 5.f && aut == 1)) {
                            k6 = 0.4f;
                        } else if ((lp.noisecc < 6.f && aut == 0) || (maxccoarse < 6.f && aut == 1)) {
                            k6 = 0.7f;
                        } else {
                            k6 = 1.f;
                        }

                        variC[6] = rtengine::max(0.00001f, k6 * variC[6]);
                        variCb[6] = rtengine::max(0.00001f, k6 * variCb[6]);

                    }

                    float* noisevarchrom = new float[GH * GW];
                    //noisevarchrom in function chroma
                    int GW2 = (GW + 1) / 2;
                    float nvch = 0.6f;//high value
                    float nvcl = 0.1f;//low value

                    if ((lp.noisecf > 100.f && aut == 0) || (maxcfine > 100.f && (aut == 1 || aut == 2))) {
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

                    if ((lp.noisecc < 2.f && aut == 0) || (maxccoarse < 0.1f && (aut == 1 || aut == 2)))  {
                        WaveletDenoiseAllAB(Ldecomp, adecomp, noisevarchrom, madL, variC, edge, noisevarab_r, true, false, false, numThreads);
                        WaveletDenoiseAllAB(Ldecomp, bdecomp, noisevarchrom, madL, variCb, edge, noisevarab_r, true, false, false, numThreads);
                    } else {
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
                if ((lp.noiself >= 0.01f ||  lp.noiself0 >= 0.01f ||  lp.noiself2 >= 0.01f || lp.noiselc >= 0.01f) && levred == 7 && lp.noiseldetail != 100.f) {
                    fftw_denoise(GW, GH, max_numblox_W, min_numblox_W, tmp1.L, Lin,  numThreads, lp, 0);
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
                    fftw_denoise(GW, GH, max_numblox_W, min_numblox_W, tmp1.a, Ain,  numThreads, lp, 1);
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
                    fftw_denoise(GW, GH, max_numblox_W, min_numblox_W, tmp1.b, Bin,  numThreads, lp, 1);
                }

            }

           // DeNoise_Local(call, lp,  originalmaskbl, levred, huerefblur, lumarefblur, chromarefblur, original, transformed, tmp1, cx, cy, sk);
            if(lp.smasktyp != 0) {
                DeNoise_Local(call, lp,  originalmaskbl, levred, huerefblur, lumarefblur, chromarefblur, original, transformed, tmp1, cx, cy, sk);
            } else {
                DeNoise_Local(call, lp,  original, levred, huerefblur, lumarefblur, chromarefblur, original, transformed, tmp1, cx, cy, sk);
            }

        } else if (call == 2) { //simpleprocess

            int bfh = int (lp.ly + lp.lyT) + del; //bfw bfh real size of square zone
            int bfw = int (lp.lx + lp.lxL) + del;

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


                int begy = lp.yc - lp.lyT;
                int begx = lp.xc - lp.lxL;
                int yEn = lp.yc + lp.ly;
                int xEn = lp.xc + lp.lx;

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
                    //     float mxsfl = 0.f;

                    if (aut == 0) {
                        if (levred == 7) {
                            edge = 2;
                            vari[0] = 0.8f * SQR((lp.noiself0 / 125.0) * (1.0 + lp.noiself0 / 25.0));
                            vari[1] = 0.8f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                            vari[2] = 0.8f * SQR((lp.noiself2 / 125.0) * (1.0 + lp.noiself2 / 25.0));

                            vari[3] = 0.8f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                            vari[4] = 0.8f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                            vari[5] = 0.8f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                            vari[6] = 0.8f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                        } else if (levred == 4) {
                            edge = 3;
                            vari[0] = 0.8f * SQR((lp.noiself0 / 125.0) * (1.0 + lp.noiself0 / 25.0));
                            vari[1] = 0.8f * SQR((lp.noiself / 125.0) * (1.0 + lp.noiself / 25.0));
                            vari[2] = 0.8f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));
                            vari[3] = 0.8f * SQR((lp.noiselc / 125.0) * (1.0 + lp.noiselc / 25.0));

                        }
                    } else if (aut == 1 || aut == 2) {
                        edge = 2;
                        vari[0] = SQR(slidL[0]);
                        vari[1] = SQR(slidL[1]);
                        vari[2] = SQR(slidL[2]);
                        vari[3] = SQR(slidL[3]);
                        vari[4] = SQR(slidL[4]);
                        vari[5] = SQR(slidL[5]);
                        vari[6] = SQR(slidL[6]);
                        float mxslid34 = rtengine::max(slidL[3], slidL[4]);
                        float mxslid56 = rtengine::max(slidL[5], slidL[6]);
                        mxsl = rtengine::max(mxslid34, mxslid56);

                    }

                    {
                        float kr3 = 0.f;
                        float kr4 = 0.f;
                        float kr5 = 0.f;

                        if (aut == 0 || aut == 1) {
                            if ((lp.noiselc < 30.f && aut == 0) || (mxsl < 30.f && aut == 1)) {
                                kr3 = 0.f;
                                kr4 = 0.f;
                                kr5 = 0.f;
                            } else if ((lp.noiselc < 50.f && aut == 0) || (mxsl < 50.f && aut == 1)) {
                                kr3 = 0.5f;
                                kr4 = 0.3f;
                                kr5 = 0.2f;
                            } else if ((lp.noiselc < 70.f && aut == 0) || (mxsl < 70.f && aut == 1)) {
                                kr3 = 0.7f;
                                kr4 = 0.5f;
                                kr5 = 0.3f;
                            } else {
                                kr3 = 1.f;
                                kr4 = 1.f;
                                kr5 = 1.f;
                            }
                        } else if (aut == 2) {
                            kr3 = 1.f;
                            kr4 = 1.f;
                            kr5 = 1.f;

                        }

                        vari[0] = rtengine::max(0.000001f, vari[0]);
                        vari[1] = rtengine::max(0.000001f, vari[1]);
                        vari[2] = rtengine::max(0.000001f, vari[2]);
                        vari[3] = rtengine::max(0.000001f, kr3 * vari[3]);

                        if (levred == 7) {
                            vari[4] = rtengine::max(0.000001f, kr4 * vari[4]);
                            vari[5] = rtengine::max(0.000001f, kr5 * vari[5]);
                            vari[6] = rtengine::max(0.000001f, kr5 * vari[6]);
                        }

                        //    float* noisevarlum = nullptr;  // we need a dummy to pass it to WaveletDenoiseAllL
                        float* noisevarlum = new float[bfh * bfw];
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


                        if ((lp.noiselc < 0.02f && aut == 0) || (mxsl < 1.f && (aut == 1 || aut == 2))) {
                            WaveletDenoiseAllL(Ldecomp, noisevarlum, madL, vari, edge, numThreads);
                        } else {
                            WaveletDenoiseAll_BiShrinkL(Ldecomp, noisevarlum, madL, vari, edge, numThreads);
                            WaveletDenoiseAllL(Ldecomp, noisevarlum, madL, vari, edge, numThreads);
                        }

                        delete [] noisevarlum;

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
                    float maxcfine = 0.f;
                    float maxccoarse = 0.f;

                    if (aut == 0) {

                        if (levred == 7) {
                            edge = 2;
                            variC[0] = SQR(noisecfr);
                            variC[1] = SQR(noisecfr);
                            variC[2] = SQR(noisecfr);

                            variC[3] = SQR(noisecfr);
                            variC[4] = SQR(noisecfr);
                            variC[5] = SQR(noiseccr);
                            variC[6] = SQR(noiseccr);

                            variCb[0] = SQR(noisecfb);
                            variCb[1] = SQR(noisecfb);
                            variCb[2] = SQR(noisecfb);

                            variCb[3] = SQR(noisecfb);
                            variCb[4] = SQR(noisecfb);
                            variCb[5] = SQR(noiseccb);
                            variCb[6] = SQR(noiseccb);

                        } else if (levred == 4) {
                            edge = 3;
                            variC[0] = SQR(lp.noisecf / 10.0);
                            variC[1] = SQR(lp.noisecf / 10.0);
                            variC[2] = SQR(lp.noisecf / 10.0);
                            variC[3] = SQR(lp.noisecf / 10.0);

                            variCb[0] = SQR(lp.noisecf / 10.0);
                            variCb[1] = SQR(lp.noisecf / 10.0);
                            variCb[2] = SQR(lp.noisecf / 10.0);
                            variCb[3] = SQR(lp.noisecf / 10.0);


                        }
                    } else if (aut == 1 || aut == 2) {
                        edge = 2;
                        variC[0] = SQR(slida[0]);
                        variC[1] = SQR(slida[1]);
                        variC[2] = SQR(slida[2]);
                        variC[3] = SQR(slida[3]);
                        variC[4] = SQR(slida[4]);
                        variC[5] = SQR(slida[5]);
                        variC[6] = SQR(slida[6]);
                        float maxc01 = rtengine::max(slida[0], slida[1]);
                        float maxc23 = rtengine::max(slida[2], slida[3]);
                        float max03 = rtengine::max(maxc01, maxc23);
                        float maxrf = rtengine::max(max03, slida[4]);
                        float maxrc = rtengine::max(slida[5], slida[6]);

                        variCb[0] = SQR(slidb[0]);
                        variCb[1] = SQR(slidb[1]);
                        variCb[2] = SQR(slidb[2]);
                        variCb[3] = SQR(slidb[3]);
                        variCb[4] = SQR(slidb[4]);
                        variCb[5] = SQR(slidb[5]);
                        variCb[6] = SQR(slidb[6]);
                        float maxb01 = rtengine::max(slidb[0], slidb[1]);
                        float maxb23 = rtengine::max(slidb[2], slidb[3]);
                        float maxb03 = rtengine::max(maxb01, maxb23);
                        float maxbf = rtengine::max(maxb03, slidb[4]);
                        maxcfine = rtengine::max(maxrf, maxbf);

                        float maxbc = rtengine::max(slidb[5], slidb[6]);
                        maxccoarse = rtengine::max(maxrc, maxbc);

                    }

                    {
                        float minic = 0.000001f;

                        if (noiscfactiv) {
                            minic = 0.1f;//only for artifact shape detection
                        }

                        float k1 = 0.f;
                        float k2 = 0.f;
                        float k3 = 0.f;

                        if (aut == 0 || aut == 1) {
                            if ((lp.noisecf < 0.2f && aut == 0) || (maxcfine < 0.2f && aut == 1)) {
                                k1 = 0.05f;
                                k2 = 0.f;
                                k3 = 0.f;
                            } else if ((lp.noisecf < 0.3f && aut == 0) || (maxcfine < 0.3f && aut == 1)) {
                                k1 = 0.1f;
                                k2 = 0.0f;
                                k3 = 0.f;
                            } else if ((lp.noisecf < 0.5f && aut == 0) || (maxcfine < 0.5f && aut == 1)) {
                                k1 = 0.2f;
                                k2 = 0.1f;
                                k3 = 0.f;
                            } else if ((lp.noisecf < 0.8f && aut == 0) || (maxcfine < 0.8f && aut == 1)) {
                                k1 = 0.3f;
                                k2 = 0.25f;
                                k3 = 0.f;
                            } else if ((lp.noisecf < 1.f && aut == 0) || (maxcfine < 1.f && aut == 1)) {
                                k1 = 0.4f;
                                k2 = 0.25f;
                                k3 = 0.1f;
                            } else if ((lp.noisecf < 2.f && aut == 0) || (maxcfine < 2.f && aut == 1)) {
                                k1 = 0.5f;
                                k2 = 0.3f;
                                k3 = 0.15f;
                            } else if ((lp.noisecf < 3.f && aut == 0) || (maxcfine < 3.f && aut == 1)) {
                                k1 = 0.6f;
                                k2 = 0.45f;
                                k3 = 0.3f;
                            } else if ((lp.noisecf < 4.f && aut == 0) || (maxcfine < 4.f && aut == 1)) {
                                k1 = 0.7f;
                                k2 = 0.5f;
                                k3 = 0.4f;
                            } else if ((lp.noisecf < 5.f && aut == 0) || (maxcfine < 5.f && aut == 1)) {
                                k1 = 0.8f;
                                k2 = 0.6f;
                                k3 = 0.5f;
                            } else if ((lp.noisecf < 6.f && aut == 0) || (maxcfine < 10.f && aut == 1)) {
                                k1 = 0.85f;
                                k2 = 0.7f;
                                k3 = 0.6f;
                            } else if ((lp.noisecf < 8.f && aut == 0) || (maxcfine < 20.f && aut == 1)) {
                                k1 = 0.9f;
                                k2 = 0.8f;
                                k3 = 0.7f;
                            } else if ((lp.noisecf < 10.f && aut == 0) || (maxcfine < 50.f && aut == 1)) {
                                k1 = 1.f;
                                k2 = 1.f;
                                k3 = 0.9f;

                            } else {
                                k1 = 1.f;
                                k2 = 1.f;
                                k3 = 1.f;
                            }
                        } else if (aut == 2) {
                            k1 = 1.f;
                            k2 = 1.f;
                            k3 = 1.f;
                        }

                        variC[0] = rtengine::max(minic, variC[0]);
                        variC[1] = rtengine::max(minic, k1 * variC[1]);
                        variC[2] = rtengine::max(minic, k2 * variC[2]);
                        variC[3] = rtengine::max(minic, k3 * variC[3]);

                        variCb[0] = rtengine::max(minic, variCb[0]);
                        variCb[1] = rtengine::max(minic, k1 * variCb[1]);
                        variCb[2] = rtengine::max(minic, k2 * variCb[2]);
                        variCb[3] = rtengine::max(minic, k3 * variCb[3]);

                        if (levred == 7) {
                            float k4 = 0.f;
                            float k5 = 0.f;
                            float k6 = 0.f;

                            if ((lp.noisecc < 0.2f && aut == 0) || (maxccoarse < 0.2f && aut == 1)) {
                                k4 = 0.1f;
                                k5 = 0.02f;
                            } else if ((lp.noisecc < 0.5f && aut == 0) || (maxccoarse < 0.5f && aut == 1)) {
                                k4 = 0.15f;
                                k5 = 0.05f;
                            } else if ((lp.noisecc < 1.f && aut == 0) || (maxccoarse < 1.f && aut == 1)) {
                                k4 = 0.15f;
                                k5 = 0.1f;
                            } else if ((lp.noisecc < 3.f && aut == 0) || (maxccoarse < 3.f && aut == 1)) {
                                k4 = 0.3f;
                                k5 = 0.15f;
                            } else if ((lp.noisecc < 4.f && aut == 0) || (maxccoarse < 5.f && aut == 1)) {
                                k4 = 0.6f;
                                k5 = 0.4f;
                            } else if ((lp.noisecc < 6.f && aut == 0) || (maxccoarse < 6.f && aut == 1)) {
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

                            if ((lp.noisecc < 4.f && aut == 0) || (maxccoarse < 4.f && aut == 1)) {
                                k6 = 0.f;
                            } else if ((lp.noisecc < 5.f && aut == 0) || (maxccoarse < 5.f && aut == 1)) {
                                k6 = 0.4f;
                            } else if ((lp.noisecc < 6.f && aut == 0) || (maxccoarse < 6.f && aut == 1)) {
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

                        if ((lp.noisecf > 30.f && aut == 0) || (maxcfine > 100.f && (aut == 1 || aut == 2))) {
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

                        if ((lp.noisecc < 0.02f && aut == 0) || (maxccoarse < 0.1f && (aut == 1  || aut == 2)))  {
                            WaveletDenoiseAllAB(Ldecomp, adecomp, noisevarchrom, madL, variC, edge, noisevarab_r, true, false, false, numThreads);
                            WaveletDenoiseAllAB(Ldecomp, bdecomp, noisevarchrom, madL, variCb, edge, noisevarab_r, true, false, false, numThreads);
                        } else {
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


                    if ((lp.noiself >= 0.01f ||  lp.noiself0 >= 0.01f ||  lp.noiself2 >= 0.01f || lp.noiselc >= 0.01f) && levred == 7 && lp.noiseldetail != 100.f) {
                        fftw_denoise(bfw, bfh, max_numblox_W, min_numblox_W, bufwv.L, Lin,  numThreads, lp, 0);
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
                        fftw_denoise(bfw, bfh, max_numblox_W, min_numblox_W, bufwv.a, Ain,  numThreads, lp, 1);
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
                        fftw_denoise(bfw, bfh, max_numblox_W, min_numblox_W, bufwv.b, Bin,  numThreads, lp, 1);
                    }
                }

            if(lp.smasktyp != 0) {
                DeNoise_Local(call, lp,  originalmaskbl, levred, huerefblur, lumarefblur, chromarefblur, original, transformed, bufwv, cx, cy, sk);
            } else {
                DeNoise_Local(call, lp,  original, levred, huerefblur, lumarefblur, chromarefblur, original, transformed, bufwv, cx, cy, sk);
            }

              //  DeNoise_Local(call, lp,  originalmaskbl, levred, huerefblur, lumarefblur, chromarefblur, original, transformed, bufwv, cx, cy, sk);
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

void clarimerge(struct local_params& lp, float &mL, float &mC, bool &exec, LabImage *tmpresid, int wavelet_level, int sk, int numThreads)
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

void ImProcFunctions::Lab_Local(
    int call, int sp, float** shbuffer, LabImage * original, LabImage * transformed, LabImage * reserved, LabImage * lastorig, int cx, int cy, int oW, int oH, int sk,
    const LocretigainCurve& locRETgainCcurve, const LocretitransCurve& locRETtransCcurve,
    const LUTf& lllocalcurve, bool locallutili,
    const LUTf& cllocalcurve, bool localclutili,
    const LUTf& lclocalcurve, bool locallcutili,
    const LocLHCurve& loclhCurve,  const LocHHCurve& lochhCurve, const LocCHCurve& locchCurve,
    const LUTf& lmasklocalcurve, bool localmaskutili,
    const LUTf& lmaskexplocalcurve, bool localmaskexputili,
    const LUTf& lmaskSHlocalcurve, bool localmaskSHutili,
    const LUTf& lmaskviblocalcurve, bool localmaskvibutili,
    const LUTf& lmasktmlocalcurve, bool localmasktmutili,
    LUTf& lmaskretilocalcurve, bool localmaskretiutili,
    const LUTf& lmaskcblocalcurve, bool localmaskcbutili,
    const LUTf& lmaskbllocalcurve, bool localmaskblutili,
    const LUTf& lmasklclocalcurve, bool localmasklcutili,
    const LUTf& lmasklocal_curve, bool localmask_utili,

    const LocCCmaskCurve& locccmasCurve, bool lcmasutili, const LocLLmaskCurve& locllmasCurve, bool llmasutili, const LocHHmaskCurve& lochhmasCurve, bool lhmasutili, const LocHHmaskCurve& lochhhmasCurve, bool lhhmasutili,
    const LocCCmaskCurve& locccmasexpCurve, bool lcmasexputili, const LocLLmaskCurve& locllmasexpCurve, bool llmasexputili, const LocHHmaskCurve& lochhmasexpCurve, bool lhmasexputili,
    const LocCCmaskCurve& locccmasSHCurve, bool lcmasSHutili, const LocLLmaskCurve& locllmasSHCurve, bool llmasSHutili, const LocHHmaskCurve& lochhmasSHCurve, bool lhmasSHutili,
    const LocCCmaskCurve& locccmasvibCurve, bool lcmasvibutili, const LocLLmaskCurve& locllmasvibCurve, bool llmasvibutili, const LocHHmaskCurve& lochhmasvibCurve, bool lhmasvibutili,
    const LocCCmaskCurve& locccmascbCurve, bool lcmascbutili, const LocLLmaskCurve& locllmascbCurve, bool llmascbutili, const LocHHmaskCurve& lochhmascbCurve, bool lhmascbutili,
    const LocCCmaskCurve& locccmasretiCurve, bool lcmasretiutili, const LocLLmaskCurve& locllmasretiCurve, bool llmasretiutili, const LocHHmaskCurve& lochhmasretiCurve, bool lhmasretiutili,
    const LocCCmaskCurve& locccmastmCurve, bool lcmastmutili, const LocLLmaskCurve& locllmastmCurve, bool llmastmutili, const LocHHmaskCurve& lochhmastmCurve, bool lhmastmutili,
    const LocCCmaskCurve& locccmasblCurve, bool lcmasblutili, const LocLLmaskCurve& locllmasblCurve, bool llmasblutili, const LocHHmaskCurve& lochhmasblCurve, bool lhmasblutili,
    const LocCCmaskCurve& locccmaslcCurve, bool lcmaslcutili, const LocLLmaskCurve& locllmaslcCurve, bool llmaslcutili, const LocHHmaskCurve& lochhmaslcCurve, bool lhmaslcutili,
    const LocCCmaskCurve& locccmas_Curve, bool lcmas_utili, const LocLLmaskCurve& locllmas_Curve, bool llmas_utili, const LocHHmaskCurve& lochhmas_Curve, bool lhmas_utili,
    const LocHHmaskCurve& lochhhmas_Curve, bool lhhmas_utili,
    const LocwavCurve& loclmasCurveblwav, bool lmasutiliblwav,
    const LocwavCurve& loclmasCurvecolwav, bool lmasutilicolwav,
    const LocwavCurve& locwavCurve, bool locwavutili,
    const LocwavCurve& loclevwavCurve, bool loclevwavutili,
    const LocwavCurve& locconwavCurve, bool locconwavutili,
    const LocwavCurve& loccompwavCurve, bool loccompwavutili,
    const LocwavCurve& loccomprewavCurve, bool loccomprewavutili,
    const LocwavCurve& locwavCurveden, bool locwavdenutili,
    const LocwavCurve& locedgwavCurve, bool locedgwavutili,
    const LocwavCurve& loclmasCurve_wav, bool lmasutili_wav,
    
    bool LHutili, bool HHutili, bool CHutili, const LUTf& cclocalcurve, bool localcutili, const LUTf& rgblocalcurve, bool localrgbutili, bool localexutili, const LUTf& exlocalcurve, const LUTf& hltonecurveloc, const LUTf& shtonecurveloc, const LUTf& tonecurveloc, const LUTf& lightCurveloc,
    double& huerefblur, double& chromarefblur, double& lumarefblur, double& hueref, double& chromaref, double& lumaref, double& sobelref, int &lastsav,
    bool prevDeltaE, int llColorMask, int llColorMaskinv, int llExpMask, int llExpMaskinv, int llSHMask, int llSHMaskinv, int llvibMask, int lllcMask, int llsharMask, int llcbMask, int llretiMask, int llsoftMask, int lltmMask, int llblMask, int ll_Mask,
    float& minCD, float& maxCD, float& mini, float& maxi, float& Tmean, float& Tsigma, float& Tmin, float& Tmax
    )
{
    //general call of others functions : important return hueref, chromaref, lumaref
    if (!params->locallab.enabled) {
        return;
    }

    BENCHFUN

    constexpr int del = 3; // to avoid crash with [loy - begy] and [lox - begx] and bfh bfw  // with gtk2 [loy - begy-1] [lox - begx -1 ] and del = 1
    struct local_params lp;
    calcLocalParams(sp, oW, oH, params->locallab, lp, prevDeltaE, llColorMask, llColorMaskinv, llExpMask, llExpMaskinv, llSHMask, llSHMaskinv, llvibMask, lllcMask, llsharMask, llcbMask, llretiMask, llsoftMask, lltmMask, llblMask, ll_Mask, locwavCurveden, locwavdenutili);

    const float radius = lp.rad / (sk * 1.4f); //0 to 70 ==> see skip
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

    if (lp.excmet == 1 && call <= 3) {//exclude
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
                sombel += val;
                deltasobelL[ir][jr] = val;
            }
        }

        const float meansob = sombel / ncsobel;
        Exclude_Local(deltasobelL, hueref, chromaref, lumaref, sobelref, meansob, lp, original, transformed, &bufreserv, reserved, cx, cy, sk);
    }

//encoding lab at the beginning
    if (lp.logena) {
        const int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
        const int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
        const int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
        const int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
        const int bfh = yend - ystart;
        const int bfw = xend - xstart;

        if (bfh >= mSP && bfw >= mSP) {
            const std::unique_ptr<LabImage> bufexporig(new LabImage(bfw, bfh)); //buffer for data in zone limit
            const std::unique_ptr<LabImage> bufexpfin(new LabImage(bfw, bfh)); //buffer for data in zone limit

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

            bufexpfin->CopyFrom(bufexporig.get(), multiThread);
            std::unique_ptr<Imagefloat> tmpImage(new Imagefloat(bfw, bfh));
            lab2rgb(*bufexpfin, *tmpImage, params->icm.workingProfile);
            log_encode(tmpImage.get(), lp, multiThread, bfw, bfh);
            rgb2lab(*(tmpImage.get()), *bufexpfin, params->icm.workingProfile);
            tmpImage.reset();

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
            transit_shapedetect2(call, 11, bufexporig.get(), bufexpfin.get(), nullptr, hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);

            if (params->locallab.spots.at(sp).recurs) {
                original->CopyFrom(transformed, multiThread);
                float avge;
                calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
            }
        }
    }

//Prepare mask for Blur and noise and Denoise
    bool denoiz = false;

    if ((lp.noiself > 0.f || lp.noiself0 > 0.f || lp.noiself2 > 0.f || lp.noiselc > 0.f || lp.noisecf > 0.f || lp.noisecc > 0.f  || lp.bilat > 0.f) && lp.denoiena) {
        denoiz = true;
    }

    bool blurz = false;
    bool delt = params->locallab.spots.at(sp).deltae;
    bool astool = params->locallab.spots.at(sp).toolbl;

    if (((radius > 1.5 * GAUSS_SKIP)  || lp.stren > 0.1 || lp.blmet == 1 || lp.guidb > 1 || lp.showmaskblmet == 2  || lp.enablMask || lp.showmaskblmet == 3 || lp.showmaskblmet == 4) && lp.blurena) {
        blurz = true;
    }

    const int GW = transformed->W;
    const int GH = transformed->H;

    LabImage * originalmaskbl = nullptr;
    std::unique_ptr<LabImage> bufmaskorigbl;
    std::unique_ptr<LabImage> bufmaskblurbl;
    std::unique_ptr<LabImage> bufgb;
    std::unique_ptr<LabImage> bufprov(new LabImage(GW, GH));

    if (denoiz || blurz || lp.denoiena || lp.blurena) {
        bufgb.reset(new LabImage(GW, GH));

        if (lp.showmaskblmet == 2  || lp.enablMask || lp.showmaskblmet == 3 || lp.showmaskblmet == 4) {
            bufmaskorigbl.reset(new LabImage(GW, GH));
            bufmaskblurbl.reset(new LabImage(GW, GH));
            originalmaskbl = new LabImage(GW, GH);
        }

        array2D<float> ble(GW, GH);
        array2D<float> blechro(GW, GH);
        array2D<float> hue(GW, GH);
        array2D<float> guid(GW, GH);
        float meanfab, fab;
        mean_fab(0, 0, GW, GH, bufgb.get(), original, fab, meanfab, lp.chromabl, multiThread);
        float chromult =  1.f - 0.01f * lp.chromabl;

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
        for (int y = 0; y < GH; y++) {
            for (int x = 0; x < GW; x++) {
                bufgb->L[y][x] = original->L[y][x];
                bufgb->a[y][x] = original->a[y][x];
                bufgb->b[y][x] = original->b[y][x];
            }
        }

        const float strumask = 0.02f * params->locallab.spots.at(sp).strumaskbl;
        JaggedArray<float> blendstru(GW, GH);

        if (lp.showmaskblmet == 2  || lp.enablMask || lp.showmaskblmet == 3 || lp.showmaskblmet == 4) {
            if (strumask > 0.f) {
                float delstrumask = 4.1f - strumask;//4.1 = 2 * max slider strumask + 0.1
                buildBlendMask(bufgb->L, blendstru, GW, GH, delstrumask);
                const float radblur = 0.02f * 0.1f * std::fabs(lp.radmabl);
                const float rm = radblur / sk;

                if (rm > 0) {
#ifdef _OPENMP
                    #pragma omp parallel if (multiThread)
#endif
                    {
                        gaussianBlur(blendstru, blendstru, GW, GH, rm);
                    }
                }
            }

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int ir = 0; ir < GH; ir++) {
                for (int jr = 0; jr < GW; jr++) {
                    float kmaskLexp = 0.f;
                    float kmaskCH = 0.f;
                    float kmasstru = 0.f;

                    if (strumask > 0.f && !astool) {
                        kmasstru = bufgb->L[ir][jr] * blendstru[ir][jr];
                    }

                    if (locllmasblCurve && llmasblutili) {
                        const float ligh = bufgb->L[ir][jr] / 32768.f;
                        kmaskLexp = 32768.f * LIM01(1.f - locllmasblCurve[500.f * ligh]);
                    }

                    if (lp.showmaskblmet != 4) {
                        if (locccmasblCurve && lcmasblutili) {
                            const float chromask = 0.0001f + std::sqrt(SQR((bufgb->a[ir][jr]) / fab) + SQR((bufgb->b[ir][jr]) / fab));
                            kmaskCH = LIM01(1.f - locccmasblCurve[500.f *  chromask]);
                        }
                    }

                    if (lochhmasblCurve && lhmasblutili) {
                        const float huema = xatan2f(bufgb->b[ir][jr], bufgb->a[ir][jr]);
                        float h = Color::huelab_to_huehsv2(huema);
                        h += 1.f / 6.f;

                        if (h > 1.f) {
                            h -= 1.f;
                        }

                        const float valHH = LIM01(1.f - lochhmasblCurve[500.f *  h]);

                        if (lp.showmaskblmet != 4) {
                            kmaskCH += chromult * valHH;
                        }

                        kmaskLexp += 32768.f * valHH;
                    }

                    bufmaskblurbl->L[ir][jr] = clipLoc(kmaskLexp + kmasstru);
                    bufmaskblurbl->a[ir][jr] = kmaskCH;
                    bufmaskblurbl->b[ir][jr] = kmaskCH;
                    ble[ir][jr] = bufmaskblurbl->L[ir][jr] / 32768.f;
                    hue[ir][jr] = xatan2f(bufmaskblurbl->b[ir][jr], bufmaskblurbl->a[ir][jr]);
                    const float chromah = std::sqrt(SQR(bufmaskblurbl->b[ir][jr]) + SQR(bufmaskblurbl->a[ir][jr]));
                    blechro[ir][jr] = chromah / 32768.f;
                    guid[ir][jr] = Color::L2Y(bufgb->L[ir][jr]) / 32768.f;
                }
            }

            const std::unique_ptr<LabImage> bufprov(new LabImage(GW, GH));

            bufprov->CopyFrom(bufmaskblurbl.get(), multiThread);

            if (lp.radmabl != 0.f) {
                float blur = lp.radmabl;
                blur = blur < 0.f ? -1.f / blur : 1.f + blur;
                const int r1 = rtengine::max<int>(4 / sk * blur + 0.5f, 1);
                const int r2 = rtengine::max<int>(25 / sk * blur + 0.5f, 1);

                constexpr float epsilmax = 0.005f;
                constexpr float epsilmin = 0.00001f;

                const float aepsil = (epsilmax - epsilmin) / 100.f;
                const float bepsil = epsilmin; //epsilmax - 100.f * aepsil;
                const float epsil = lp.radmabl < 0.f ? 0.001f : aepsil * lp.radmabl + bepsil;

                rtengine::guidedFilter(guid, blechro, blechro, r1, epsil, multiThread);
                rtengine::guidedFilter(guid, ble, ble, r2, 0.2 * epsil, multiThread);

                //    guidedFilter(guid, ble, ble, lp.radmabl * 10.f / sk, 0.001, multiThread, 4);
            }

            LUTf lutTonemaskbl(65536);
            calcGammaLut(lp.gammabl, lp.slomabl, lutTonemaskbl);

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int ir = 0; ir < GH; ir++) {
                for (int jr = 0; jr < GW; jr++) {
                    const float2 sincosval = xsincosf(hue[ir][jr]);
                    bufmaskblurbl->L[ir][jr] = LIM01(ble[ir][jr]) * 32768.f;
                    const float L_ = 2.f * bufmaskblurbl->L[ir][jr];
                    bufmaskblurbl->L[ir][jr] = lutTonemaskbl[L_];
                    bufmaskblurbl->a[ir][jr] = 32768.f * sincosval.y * blechro[ir][jr];
                    bufmaskblurbl->b[ir][jr] = 32768.f * sincosval.x * blechro[ir][jr];
                }
            }
        }

        if (strumask > 0.f && astool && (lp.enablMask || lp.showmaskblmet == 3)) {

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int ir = 0; ir < GH; ir++) {
                for (int jr = 0; jr < GW; jr++) {
                    bufmaskblurbl->L[ir][jr] *= (1.f + blendstru[ir][jr]);
                }
            }
        }

        if (lmaskbllocalcurve && localmaskblutili && (lp.enablMask || lp.showmaskblmet == 3)) {
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int ir = 0; ir < GH; ir++)
                for (int jr = 0; jr < GW; jr++) {
                    bufmaskblurbl->L[ir][jr] = 0.5f * lmaskbllocalcurve[2.f * bufmaskblurbl->L[ir][jr]];
                }
        }

        const int highli = params->locallab.spots.at(sp).shadmaskbl;

        if (highli > 0 && (lp.enablMask || lp.showmaskblmet == 3)) {
            ImProcFunctions::shadowsHighlights(bufmaskblurbl.get(), true, 1, highli, 0, 40, sk, 50, 0);
        }

        const int shado = params->locallab.spots.at(sp).shadmaskblsha;

        if (shado > 0 && (lp.enablMask || lp.showmaskblmet == 3)) {
            ImProcFunctions::shadowsHighlights(bufmaskblurbl.get(), true, 1, 0, shado, 40, sk, 0, 60);
         }

        int wavelet_level = params->locallab.spots.at(sp).shadmaskbl;
        int maxlvl = wavelet_level;

        int minwin = rtengine::min(GW, GH);
        int maxlevelspot = 9;

        while ((1 << maxlevelspot) >= (minwin * sk) && maxlevelspot  > 1) {
            --maxlevelspot ;
        }

        wavelet_level = rtengine::min(wavelet_level, maxlevelspot);
        bool wavcurvemask = false;

        if (loclmasCurveblwav && lmasutiliblwav && (lp.enablMask || lp.showmaskblmet == 3)) {
            for (int i = 0; i < 500; i++) {
                if (loclmasCurveblwav[i] != 0.5) {
                    wavcurvemask = true;
                }
            }
        }

        if (wavcurvemask && (lp.enablMask || lp.showmaskblmet == 3)) {
            const int level_bl = params->locallab.spots.at(sp).csthresholdblur.getBottomLeft();
            const int level_hl = params->locallab.spots.at(sp).csthresholdblur.getTopLeft();
            const int level_br = params->locallab.spots.at(sp).csthresholdblur.getBottomRight();
            const int level_hr = params->locallab.spots.at(sp).csthresholdblur.getTopRight();

#ifdef _OPENMP
            const int numThreads = omp_get_max_threads();
#else
            const int numThreads = 1;

#endif

            wavelet_decomposition *wdspotbl = new wavelet_decomposition(bufmaskblurbl->L[0], GW, GH, maxlvl, 1, sk, numThreads, lp.daubLen);
            if (wdspotbl->memory_allocation_failed()) {
                return;
            }
    
 
            float mean[10];
            float meanN[10];
            float sigma[10];
            float sigmaN[10];
            float MaxP[10];
            float MaxN[10];

            Evaluate2(*wdspotbl, mean, meanN, sigma, sigmaN, MaxP, MaxN, numThreads);
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
                    int W_L = wdspotbl->level_W(level);
                    int H_L = wdspotbl->level_H(level);
                    float* const *wav_L = wdspotbl->level_coeffs(level);
               
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
                            if(loclmasCurveblwav && lmasutiliblwav) {
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
                                float kc = klev * (loclmasCurveblwav[absciss * 500.f] - 0.5f);
                                float amplieffect = kc <= 0.f ? 1.f : 4.f;

                                float kinterm = 1.f + amplieffect * kc;
                                kinterm = kinterm <= 0.f ? 0.01f : kinterm;

                                val *=  kinterm;
                            
                            }
                        }
                    }
                
                }
            }
            wdspotbl->reconstruct(bufmaskblurbl->L[0], 1.f);
            delete wdspotbl;

        }


        // deltae Mask with scope
        int sco = params->locallab.spots.at(sp).scopemask;
        const float mindE = 2.f + MINSCOPE * sco * lp.thr;
        const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
        const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
        const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);

        if (delt && lp.blurmet == 0 && (lp.enablMask || lp.showmaskblmet == 3)) {
            JaggedArray<float> rdE(GW, GH);
            deltaEforMask(rdE, GW, GH, bufgb.get(), hueref, chromaref, lumaref, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, lp.balance, lp.balanceh);
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int ir = 0; ir < GH; ir++) {
                for (int jr = 0; jr < GW; jr++) {
                    bufmaskblurbl->L[ir][jr] = bufprov->L[ir][jr] + rdE[ir][jr] * (bufmaskblurbl->L[ir][jr] - bufprov->L[ir][jr]);
                    bufmaskblurbl->a[ir][jr] = bufprov->a[ir][jr] + rdE[ir][jr] * (bufmaskblurbl->a[ir][jr] - bufprov->a[ir][jr]);
                    bufmaskblurbl->b[ir][jr] = bufprov->b[ir][jr] + rdE[ir][jr] * (bufmaskblurbl->b[ir][jr] - bufprov->b[ir][jr]);
                }
            }
        }

        const float lap = params->locallab.spots.at(sp).lapmaskbl;
        const bool pde = params->locallab.spots.at(sp).laplac;
        const float lumask = params->locallab.spots.at(sp).lumask;

        if (lap > 0.f && (lp.enablMask || lp.showmaskblmet == 3)) {
            const float *datain = bufmaskblurbl->L[0];
            const std::unique_ptr<float[]> data_tmp(new float[GH * GW]);

            if (!pde) {
                ImProcFunctions::discrete_laplacian_threshold(data_tmp.get(), datain, GW, GH, 200.f * lap);
            } else {
                ImProcFunctions::retinex_pde(datain, data_tmp.get(), GW, GH, 12.f * lap, 1.f, nullptr, 0, 0, 1);
            }

#ifdef _OPENMP
            #pragma omp parallel for if (multiThread)
#endif
            for (int y = 0; y < GH; y++) {
                for (int x = 0; x < GW; x++) {
                    bufmaskblurbl->L[y][x] = data_tmp[y * GW + x];
                }
            }
        }

        const float radiusb = 1.f / sk;

        if (lp.showmaskblmet == 2 || lp.enablMask || lp.showmaskblmet == 3 || lp.showmaskblmet == 4) {
            const int invers = lp.blurmet == 1 ? 1 : 0;

#ifdef _OPENMP
            #pragma omp parallel if (multiThread)
#endif
            {
                gaussianBlur(bufmaskblurbl->L, bufmaskorigbl->L, GW, GH, radiusb);
                gaussianBlur(bufmaskblurbl->a, bufmaskorigbl->a, GW, GH, 1.f + (0.005f * lp.radmabl) / sk);
                gaussianBlur(bufmaskblurbl->b, bufmaskorigbl->b, GW, GH, 1.f + (0.005f * lp.radmabl) / sk);
            }

            if (lp.showmaskblmet == 0 || lp.showmaskblmet == 1 || lp.showmaskblmet == 2 || lp.showmaskblmet == 4 || lp.enablMask) {
                blendmask(lp, 0, 0, cx, cy, GW, GH, bufgb.get(), original, bufmaskorigbl.get(), originalmaskbl, lp.blendmabl, lp.blendmabl, invers);
            } else if (lp.showmaskblmet == 3) {
                showmask(lumask, lp, 0, 0, cx, cy, GW, GH, bufgb.get(), transformed, bufmaskorigbl.get(), invers);
                return; 
            }
        }

//end mask
    }

    bool execmaskblur = (lp.showmaskblmet == 2 || lp.enablMask || lp.showmaskblmet == 3 || lp.showmaskblmet == 4) && lp.smasktyp != 1;
    if (((radius > 1.5 * GAUSS_SKIP && lp.rad > 1.6) || lp.stren > 0.1 || lp.blmet == 1 || lp.guidb > 0 || execmaskblur) && lp.blurena) { // radius < GAUSS_SKIP means no gauss, just copy of original image
 //   if (((radius > 1.5 * GAUSS_SKIP && lp.rad > 1.6) || lp.stren > 0.1 || lp.blmet == 1 || lp.guidb > 0 || lp.showmaskblmet == 2 || lp.enablMask || lp.showmaskblmet == 3 || lp.showmaskblmet == 4) && lp.blurena) { // radius < GAUSS_SKIP means no gauss, just copy of original image
        std::unique_ptr<LabImage> tmp1;
        std::unique_ptr<LabImage> tmp2;
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
        int strengr = params->locallab.spots.at(sp).strengr;
        int scalegr = params->locallab.spots.at(sp).scalegr;



        if (bfw >= mSP && bfh >= mSP) {
            if (lp.blurmet == 0 && (fft || lp.rad > 30.f)) {
                optfft(N_fftwsize, bfh, bfw, bfhr, bfwr, lp, original->H, original->W, xstart, ystart, xend, yend, cx, cy);
            }

            const std::unique_ptr<LabImage> bufgbi(new LabImage(GW, GH));

            //here mask is used with plain image for normal and inverse
            //if it is possible to optimize with maskcalccol(), I don't to preserve visibility
            if (lp.showmaskblmet == 0 || lp.showmaskblmet == 1  || lp.showmaskblmet == 2 || lp.showmaskblmet == 4 || lp.enablMask) {

                if (lp.blurmet == 0) {
                    if (bfw > 0 && bfh > 0) {
                        tmp1.reset(new LabImage(bfw, bfh));
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

                    for (int y = 0; y < GH ; y++) {
                        for (int x = 0; x < GW; x++) {
                            tmp2->L[y][x] = original->L[y][x];
                            tmp2->a[y][x] = original->a[y][x];
                            tmp2->b[y][x] = original->b[y][x];
                            bufgbi->L[y][x] = original->L[y][x];
                            bufgbi->a[y][x] = original->a[y][x];
                            bufgbi->b[y][x] = original->b[y][x];
                        }
                    }

                }


                if (lp.blurmet == 0 && lp.blmet == 0 && radius > (1.5 * GAUSS_SKIP) && lp.rad > 1.6) {
                    if (fft || lp.rad > 30.f) {
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
                            if (lp.chromet == 0)
                            {
                                gaussianBlur(tmp1->L, tmp1->L, bfw, bfh, radius);
                            }

                            else if (lp.chromet == 1)
                            {
                                gaussianBlur(tmp1->a, tmp1->a, bfw, bfh, radius);
                                gaussianBlur(tmp1->b, tmp1->b, bfw, bfh, radius);
                            } else if (lp.chromet == 2)
                            {
                                gaussianBlur(tmp1->L, tmp1->L, bfw, bfh, radius);
                                gaussianBlur(tmp1->a, tmp1->a, bfw, bfh, radius);
                                gaussianBlur(tmp1->b, tmp1->b, bfw, bfh, radius);
                            }
                        }
                    }

                } else if (lp.blurmet == 1 && lp.blmet == 0 && radius > (1.5 * GAUSS_SKIP) && lp.rad > 1.6) {
                    if (fft || lp.rad > 30.f) {
                        if (lp.chromet == 0) {
                            ImProcFunctions::fftw_convol_blur2(original->L, tmp1->L, GW, GH, radius, 0, 0);
                        }

                        else if (lp.chromet == 1) {
                            ImProcFunctions::fftw_convol_blur2(original->a, tmp1->a, GW, GH, radius, 0, 0);
                            ImProcFunctions::fftw_convol_blur2(original->b, tmp1->b, GW, GH, radius, 0, 0);
                        } else if (lp.chromet == 2) {
                            ImProcFunctions::fftw_convol_blur2(original->L, tmp1->L, GW, GH, radius, 0, 0);
                            ImProcFunctions::fftw_convol_blur2(original->a, tmp1->a, GW, GH, radius, 0, 0);
                            ImProcFunctions::fftw_convol_blur2(original->b, tmp1->b, GW, GH, radius, 0, 0);
                        }

                    } else {

#ifdef _OPENMP
                        #pragma omp parallel if (multiThread)
#endif
                        {
                            if (lp.chromet == 0)
                            {
                                gaussianBlur(original->L, tmp1->L, GW, GH, radius);
                            } else if (lp.chromet == 1)

                            {
                                gaussianBlur(original->a, tmp1->a, GW, GH, radius);
                                gaussianBlur(original->b, tmp1->b, GW, GH, radius);
                            } else if (lp.chromet == 2)

                            {
                                gaussianBlur(original->L, tmp1->L, GW, GH, radius);
                                gaussianBlur(original->a, tmp1->a, GW, GH, radius);
                                gaussianBlur(original->b, tmp1->b, GW, GH, radius);
                            }
                        }
                    }
                }


                //add noise
                if (tmp1.get() && lp.stren > 0.1f && lp.blmet == 0) {
                    float mean = 0.f;//0 best result
                    float variance = lp.stren ;
                    addGaNoise(tmp1.get(), tmp1.get(), mean, variance, sk) ;
                }

                //add grain
                if (lp.blmet == 0 && strengr > 0) {
                    int wi = bfw;
                    int he = bfh;

                    if (lp.blurmet == 1) {
                        wi = GW;
                        he = GH;
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


                        filmGrain(tmpImage, isogr, strengr, scalegr, wi, he);

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
                    int wid = GW;
                    int hei = GH;
                    tmL = new float*[hei];

                    for (int i = 0; i < hei; ++i) {
                        tmL[i] = new float[wid];
                    }

                    if (lp.chromet == 0) {
                        Median_Denoise(tmp2->L, tmp1->L, GW, GH, medianTypeL, lp.it, multiThread, tmL);
                    } else if (lp.chromet == 1) {
                        Median_Denoise(tmp2->a, tmp1->a, GW, GH, medianTypeAB, lp.it, multiThread, tmL);
                        Median_Denoise(tmp2->b, tmp1->b, GW, GH, medianTypeAB, lp.it, multiThread, tmL);
                    } else if (lp.chromet == 2) {
                        Median_Denoise(tmp2->L, tmp1->L, GW, GH, medianTypeL, lp.it, multiThread, tmL);
                        Median_Denoise(tmp2->a, tmp1->a, GW, GH, medianTypeAB, lp.it, multiThread, tmL);
                        Median_Denoise(tmp2->b, tmp1->b, GW, GH, medianTypeAB, lp.it, multiThread, tmL);
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
                                bufgb->L[y - ystart][x - xstart] = original->L[y][x];
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

                        const float epsil = 0.001f * std::pow(2, - lp.epsb);

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

                        delete tmpImage;
                    }

                } else if (lp.blurmet == 1 && lp.blmet == 2) {

                    if (lp.guidb > 0) {
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < GH ; y++) {
                            for (int x = 0; x < GW; x++) {
                                tmp1->L[y][x] = original->L[y][x];
                                tmp1->a[y][x] = original->a[y][x];
                                tmp1->b[y][x] = original->b[y][x];
                                tmp2->L[y][x] = original->L[y][x];
                            }
                        }

                        Imagefloat *tmpImage = nullptr;
                        tmpImage = new Imagefloat(GW, GH);
                        lab2rgb(*tmp1, *tmpImage, params->icm.workingProfile);
                        array2D<float> LL(GW, GH);
                        array2D<float> rr(GW, GH);
                        array2D<float> gg(GW, GH);
                        array2D<float> bb(GW, GH);
                        array2D<float> guide(GW, GH);

#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int y = 0; y < GH ; y++) {
                            for (int x = 0; x < GW; x++) {
                                LL[y][x] = tmp1->L[y][x];
                                float ll = LL[y][x] / 32768.f;
                                guide[y][x] = xlin2log(rtengine::max(ll, 0.f), 10.f);
                                rr[y][x] = tmpImage->r(y, x);
                                gg[y][x] = tmpImage->g(y, x);
                                bb[y][x] = tmpImage->b(y, x);

                            }
                        }

                        array2D<float> iR(GW, GH, rr, 0);
                        array2D<float> iG(GW, GH, gg, 0);
                        array2D<float> iB(GW, GH, bb, 0);
                        array2D<float> iL(GW, GH, LL, 0);

                        int r = rtengine::max(int(lp.guidb / sk), 1);

                        const float epsil = 0.001f * std::pow(2, - lp.epsb);

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
                        for (int y = 0; y < GH ; y++) {
                            for (int x = 0; x < GW; x++) {
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
                            for (int y = 0; y < GH ; y++) {
                                for (int x = 0; x < GW; x++) {
                                    LL[y][x] = intp(lp.strbl, LL[y][x] , iL[y][x]);
                                    tmp1->L[y][x] = LL[y][x];
                                }
                            }
                        }
                        delete tmpImage;
                    }
                }

                if (tmp1.get()) {
                    JaggedArray<float> bufchro(lp.blurmet == 1 ? GW : bfw, lp.blurmet == 1 ? GH : bfh);
                    float minC = std::sqrt(SQR(tmp1->a[0][0]) + SQR(tmp1->b[0][0])) - std::sqrt(SQR(bufgb->a[0][0]) + SQR(bufgb->b[0][0]));
                    float maxC = minC;
#ifdef _OPENMP
                    #pragma omp parallel for reduction(max:maxC) reduction(min:minC) schedule(dynamic,16) if (multiThread)
#endif
                    for (int ir = 0; ir < bfh; ir++) {
                        for (int jr = 0; jr < bfw; jr++) {
                            bufchro[ir][jr] = std::sqrt(SQR(tmp1->a[ir][jr]) + SQR(tmp1->b[ir][jr])) - std::sqrt(SQR(bufgb->a[ir][jr]) + SQR(bufgb->b[ir][jr]));
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
                        for (int y = 0; y < bfh; y++) {
                            for (int x = 0; x < bfw; x++) {
                                bufchro[y][x] *= coefC;
                            }
                        }
                    }

                    if (lp.blurmet == 0) { //blur and noise (center)
//                        BlurNoise_Local(tmp1.get(), originalmaskbl, bufchro, hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);
                        
                        if(lp.smasktyp != 1) {
                            BlurNoise_Local(tmp1.get(), originalmaskbl, bufchro, hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);
                        } else {
                            BlurNoise_Local(tmp1.get(), original, bufchro, hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);
                        }

                        if (params->locallab.spots.at(sp).recurs) {
                            original->CopyFrom(transformed, multiThread);
                            float avge;
                            calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                        }
                    } else if (lp.blurmet == 1) {
 //                       InverseBlurNoise_Local(originalmaskbl, bufchro, lp, hueref, chromaref, lumaref, original, transformed, tmp1.get(), cx, cy, sk);
                        if(lp.smasktyp != 1) {
                            InverseBlurNoise_Local(originalmaskbl, bufchro, lp, hueref, chromaref, lumaref, original, transformed, tmp1.get(), cx, cy, sk);
                        } else {
                            InverseBlurNoise_Local(original, bufchro, lp, hueref, chromaref, lumaref, original, transformed, tmp1.get(), cx, cy, sk);
                        }

                        if (params->locallab.spots.at(sp).recurs) {
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

        const double threshold = lp.bilat / 20.0;

        if (bfh > 8 && bfw > 8) {
            ImProcFunctions::impulse_nr(bufwv.get(), threshold);
        }

        DeNoise_Local(call, lp,  originalmaskbl, levred, huerefblur, lumarefblur, chromarefblur, original, transformed, *(bufwv.get()), cx, cy, sk);

        if (params->locallab.spots.at(sp).recurs) {
            original->CopyFrom(transformed, multiThread);
            float avge;
            calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
        }
    }

//local denoise

    if (lp.denoiena) {
        float slidL[8] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
        float slida[8] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
        float slidb[8] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
        constexpr int aut = 0;
        DeNoise(call, del, slidL, slida, slidb, aut, noiscfactiv, lp, originalmaskbl, levred, huerefblur, lumarefblur, chromarefblur, original, transformed, cx, cy, sk);

        if (params->locallab.spots.at(sp).recurs) {
            original->CopyFrom(transformed, multiThread);
            float avge;
            calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
        }
    }

    if (denoiz || blurz || lp.denoiena || lp.blurena) {
        delete originalmaskbl;
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
                for (int y = 0; y < bfh; y++) {
                    for (int x = 0; x < bfw; x++) {
                        loctemp->L[y][x] = original->L[y + ystart][x + xstart];
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
                bool delt = params->locallab.spots.at(sp).deltae;
                int sco = params->locallab.spots.at(sp).scopemask;
                int lumask = params->locallab.spots.at(sp).lumask;
                int shado = 0;
                const float mindE = 2.f + MINSCOPE * sco * lp.thr;
                const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
                const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
                const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
                bool lmasutilicolwav = false;
                float amountcd = 0.f;
                float anchorcd = 50.f;
                int shortcu = 0; //lp.mergemet; //params->locallab.spots.at(sp).shortc;
                LocHHmaskCurve lochhhmasCurve;
                bool lhhmasutili = false;
                maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, loctemp.get(), bufmaskorigcb.get(), originalmaskcb.get(), original, reserved, inv, lp,
                            0.f, false,
                            locccmascbCurve, lcmascbutili, locllmascbCurve, llmascbutili, lochhmascbCurve, lhmascbutili, lochhhmasCurve, lhhmasutili, multiThread,
                            enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, amountcd, anchorcd, lmaskcblocalcurve, localmaskcbutili, dummy, lmasutilicolwav, 1, 1, 5, 5,
                            shortcu, delt, hueref, chromaref, lumaref,
                            maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.0f, 0.f, -1
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
                        }
                    }

                    if (lp.clarityml != 0.f && lp.mulloc[5] == 1.0) { //enabled last level to retrieve level 5 and residual image in case user not select level 5
                        lp.mulloc[5] = 1.001f;
                    }

                    if (lp.contresid != 0.f && lp.mulloc[5] == 1.0) { //enabled last level to retrieve level 5 and residual image in case user not select level 5
                        lp.mulloc[5] = 1.001f;
                    }

                    ImProcFunctions::cbdl_local_temp(bufsh, loctemp->L, bfw, bfh, lp.mulloc, 1.f, lp.threshol, lp.clarityml, lp.contresid, skinprot, false, b_l, t_l, t_r, b_r, choice, sk, multiThread);

                    if (lp.softradiuscb > 0.f) {
                        softproc(origcbdl.get(), loctemp.get(), lp.softradiuscb, bfh, bfw, 0.001, 0.00001, 0.5f, sk, multiThread, 1);
                    }

                }

                transit_shapedetect(6, loctemp.get(), originalmaskcb.get(), bufchrom, false, hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);

                bool nochroma = (lp.showmaskcbmet == 2  || lp.showmaskcbmet == 1);

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
                    bufsh.free();

                    if (params->locallab.spots.at(sp).recurs) {
                        original->CopyFrom(transformed, multiThread);
                        float avge;
                        calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                    }
                }
            }
        }
    }


//end cbdl_Local

//vibrance

    if (lp.expvib && (lp.past != 0.f  || lp.satur != 0.f || lp.strvib != 0.f  || lp.war != 0 || lp.strvibab != 0.f  || lp.strvibh != 0.f || lp.showmaskvibmet == 2 || lp.enavibMask || lp.showmaskvibmet == 3 || lp.showmaskvibmet == 4 || lp.prevdE) && lp.vibena) { //interior ellipse renforced lightness and chroma  //locallutili
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
                bool lmasutilicolwav = false;
                bool delt = params->locallab.spots.at(sp).deltae;
                int sco = params->locallab.spots.at(sp).scopemask;
                int shortcu = 0;//lp.mergemet; //params->locallab.spots.at(sp).shortc;

                const float mindE = 2.f + MINSCOPE * sco * lp.thr;
                const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
                const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
                const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
                int shado = 0;
                int lumask = params->locallab.spots.at(sp).lumask;
                LocHHmaskCurve lochhhmasCurve;
                bool lhhmasutili = false;
                float amountcd = 0.f;
                float anchorcd = 50.f;

                maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, bufexporig.get(), bufmaskorigvib.get(), originalmaskvib.get(), original, reserved, inv, lp,
                            0.f, false,
                            locccmasvibCurve, lcmasvibutili, locllmasvibCurve, llmasvibutili, lochhmasvibCurve, lhmasvibutili, lochhhmasCurve, lhhmasutili, multiThread,
                            enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, amountcd, anchorcd, lmaskviblocalcurve, localmaskvibutili, dummy, lmasutilicolwav, 1, 1, 5, 5,
                            shortcu, delt, hueref, chromaref, lumaref,
                            maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, -1
                           );

                if (lp.showmaskvibmet == 3) {
                    showmask(lumask, lp, xstart, ystart, cx, cy, bfw, bfh, bufexporig.get(), transformed, bufmaskorigvib.get(), 0);

                    return;
                }

                if (lp.showmaskvibmet == 0 || lp.showmaskvibmet == 1  || lp.showmaskvibmet == 2 || lp.showmaskvibmet == 4 || lp.enavibMask) {

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

                    VibranceParams vibranceParams;
                    vibranceParams.enabled = params->locallab.spots.at(sp).expvibrance;
                    vibranceParams.pastels = params->locallab.spots.at(sp).pastels;
                    vibranceParams.saturated = params->locallab.spots.at(sp).saturated;
                    vibranceParams.psthreshold = params->locallab.spots.at(sp).psthreshold;
                    vibranceParams.protectskins = params->locallab.spots.at(sp).protectskins;
                    vibranceParams.avoidcolorshift = params->locallab.spots.at(sp).avoidcolorshift;
                    vibranceParams.pastsattog = params->locallab.spots.at(sp).pastsattog;
                    vibranceParams.skintonescurve = params->locallab.spots.at(sp).skintonescurve;


                    bufexpfin->CopyFrom(bufexporig.get(), multiThread);

                    if (lp.strvibh != 0.f) {
                        struct grad_params gph;
                        calclocalGradientParams(lp, gph, ystart, xstart, bfw, bfh, 9);
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int ir = 0; ir < bfh; ir++)
                            for (int jr = 0; jr < bfw; jr++) {
                                double factor = 1.0;
                                factor = ImProcFunctions::calcGradientFactor(gph, jr, ir);
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
                        struct grad_params gp;
                        calclocalGradientParams(lp, gp, ystart, xstart, bfw, bfh, 7);
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int ir = 0; ir < bfh; ir++)
                            for (int jr = 0; jr < bfw; jr++) {
                                double factor = 1.0;
                                factor = ImProcFunctions::calcGradientFactor(gp, jr, ir);
                                bufexpfin->L[ir][jr] *= factor;
                            }
                    }

                    if (lp.strvibab != 0.f) {
                        struct grad_params gpab;
                        calclocalGradientParams(lp, gpab, ystart, xstart, bfw, bfh, 8);
#ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                        for (int ir = 0; ir < bfh; ir++)
                            for (int jr = 0; jr < bfw; jr++) {
                                double factor = 1.0;
                                factor = ImProcFunctions::calcGradientFactor(gpab, jr, ir);
                                bufexpfin->a[ir][jr] *= factor;
                                bufexpfin->b[ir][jr] *= factor;
                            }
                    }

                    ImProcFunctions::vibrance(bufexpfin.get(), vibranceParams, params->toneCurve.hrenabled, params->icm.workingProfile);

                    if (params->locallab.spots.at(sp).warm != 0) {
                        ImProcFunctions::ciecamloc_02float(sp, bufexpfin.get());
                    }


                    transit_shapedetect2(call, 2, bufexporig.get(), bufexpfin.get(), originalmaskvib.get(), hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);


                }


                if (params->locallab.spots.at(sp).recurs) {
                    original->CopyFrom(transformed, multiThread);
                    float avge;
                    calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                }
            }
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

            if (bfw >= mSP && bfh >= mSP) {
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

#ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                for (int y = ystart; y < yend; y++) {
                    for (int x = xstart; x < xend; x++) {
                        bufgb->L[y - ystart][x - xstart] = original->L[y][x];
                        bufgb->a[y - ystart][x - xstart] = original->a[y][x];
                        bufgb->b[y - ystart][x - xstart] = original->b[y][x];
                        bufgbm->L[y - ystart][x - xstart] = original->L[y][x];
                        bufgbm->a[y - ystart][x - xstart] = original->a[y][x];
                        bufgbm->b[y - ystart][x - xstart] = original->b[y][x];
                    }
                }

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
                    bool lmasutilicolwav = false;
                    bool delt = params->locallab.spots.at(sp).deltae;
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
                    bool lhhmasutili = false;

                    maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, bufgbm.get(), bufmaskorigtm.get(), originalmasktm.get(), original, reserved, inv, lp,
                                0.f, false,
                                locccmastmCurve, lcmastmutili, locllmastmCurve, llmastmutili, lochhmastmCurve, lhmastmutili, lochhhmasCurve, lhhmasutili, multiThread,
                                enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, amountcd, anchorcd, lmasktmlocalcurve, localmasktmutili, dummy, lmasutilicolwav, 1, 1, 5, 5,
                                shortcu, delt, hueref, chromaref, lumaref,
                                maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, -1
                               );

                    if (lp.showmasktmmet == 3) {
                        showmask(lumask, lp, xstart, ystart, cx, cy, bfw, bfh, bufgbm.get(), transformed, bufmaskorigtm.get(), 0);

                        return;
                    }
                }

                if (lp.showmasktmmet == 0 || lp.showmasktmmet == 1  || lp.showmasktmmet == 2 || lp.showmasktmmet == 4 || lp.showmasktmmet == 3 || lp.enatmMask) {
                    constexpr int itera = 0;
                    ImProcFunctions::EPDToneMaplocal(sp, bufgb.get(), tmp1.get(), itera, sk);//iterate to 0 calculate with edgstopping, improve result, call=1 dcrop we can put iterate to 5

                    tmp1m->CopyFrom(tmp1.get(), multiThread); //save current result
                    bool enatmMasktmap = params->locallab.spots.at(sp).enatmMaskaft;

                    if (enatmMasktmap) {
                        //calculate new values for original, originalmasktm, bufmaskorigtm...in function of tmp1
                        LocwavCurve dummy;
                        bool lmasutilicolwav = false;
                        bool delt = params->locallab.spots.at(sp).deltae;
                        int sco = params->locallab.spots.at(sp).scopemask;
                        int shortcu = 0;//lp.mergemet; //params->locallab.spots.at(sp).shortc;
                        int lumask = params->locallab.spots.at(sp).lumask;

                        const float mindE = 2.f + MINSCOPE * sco * lp.thr;
                        const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
                        const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
                        const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
                        int shado = 0;
                        float amountcd = 0.f;
                        float anchorcd = 50.f;
                        LocHHmaskCurve lochhhmasCurve;
                        bool lhhmasutili = false;

                        maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, tmp1.get(), bufmaskorigtm.get(), originalmasktm.get(), original, reserved, inv, lp,
                                    0.f, false,
                                    locccmastmCurve, lcmastmutili, locllmastmCurve, llmastmutili, lochhmastmCurve, lhmastmutili, lochhhmasCurve, lhhmasutili, multiThread,
                                    enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, amountcd, anchorcd, lmasktmlocalcurve, localmasktmutili, dummy, lmasutilicolwav, 1, 1, 5, 5,
                                    shortcu, delt, hueref, chromaref, lumaref,
                                    maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, -1
                                   );

                        if (lp.showmasktmmet == 3) {//display mask
                            showmask(lumask, lp, xstart, ystart, cx, cy, bfw, bfh, tmp1.get(), transformed, bufmaskorigtm.get(), 0);

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

                    //   transit_shapedetect_retinex(call, 4, bufgb.get(),bufmaskorigtm.get(), originalmasktm.get(), buflight, bufchro, hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);
                    transit_shapedetect2(call, 8, bufgb.get(), tmp1.get(), originalmasktm.get(), hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);

                    //  transit_shapedetect(8, tmp1.get(), originalmasktm.get(), bufchro, false, hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);
                    bufgb.reset();

                    if (params->locallab.spots.at(sp).recurs) {
                        original->CopyFrom(transformed, multiThread);
                        float avge;
                        calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                    }
                }
            }
        }
    }

//end TM


//shadow highlight
    bool tonequ = false;

    if (lp.mullocsh[0] != 0 || lp.mullocsh[1] != 0 || lp.mullocsh[2] != 0 || lp.mullocsh[3] != 0 || lp.mullocsh[4] != 0) {
        tonequ = true;
    }

    bool tonecurv = false;

    if (params->locallab.spots.at(sp).gamSH != 2.4 || params->locallab.spots.at(sp).sloSH != 12.92) {
        tonecurv = true;
    }

    if (! lp.invsh && (lp.highlihs > 0.f || lp.shadowhs > 0.f || tonequ || tonecurv || lp.strSH != 0.f || lp.showmaskSHmet == 2 || lp.enaSHMask || lp.showmaskSHmet == 3 || lp.showmaskSHmet == 4 || lp.prevdE) && call < 3 && lp.hsena) {
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
            bool lmasutilicolwav = false;
            bool delt = params->locallab.spots.at(sp).deltae;
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
            bool lhhmasutili = false;

            maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, bufexporig.get(), bufmaskorigSH.get(), originalmaskSH.get(), original, reserved, inv, lp,
                        0.f, false,
                        locccmasSHCurve, lcmasSHutili, locllmasSHCurve, llmasSHutili, lochhmasSHCurve, lhmasSHutili, lochhhmasCurve, lhhmasutili, multiThread,
                        enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, amountcd, anchorcd, lmaskSHlocalcurve, localmaskSHutili, dummy, lmasutilicolwav, 1, 1, 5, 5,
                        shortcu, delt, hueref, chromaref, lumaref,
                        maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, -1
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
                    for (int ir = 0; ir < bfh; ir++)
                        for (int jr = 0; jr < bfw; jr++) {
                            double factor = 1.0;
                            factor = ImProcFunctions::calcGradientFactor(gp, jr, ir);
                            bufexpfin->L[ir][jr] *= factor;
                        }
                }

                if (lp.shmeth == 1) {
                    double scal = (double)(sk);
                    Imagefloat *tmpImage = nullptr;
                    tmpImage = new Imagefloat(bfw, bfh);
                    lab2rgb(*bufexpfin, *tmpImage, params->icm.workingProfile);

                    if (tonecurv) { //Tone response curve  : does nothing if gamma=2.4 and slope=12.92 ==> gamma sRGB
                        float gamtone = params->locallab.spots.at(sp).gamSH;
                        float slotone = params->locallab.spots.at(sp).sloSH;
                        cmsHTRANSFORM dummy = nullptr;
                        workingtrc(tmpImage, tmpImage, bfw, bfh, -5, params->icm.workingProfile, 2.4, 12.92310, dummy, true, false, false);
                        workingtrc(tmpImage, tmpImage, bfw, bfh, 5, params->icm.workingProfile, gamtone, slotone, dummy, false, true, true);
                    }

                    if (tonequ) {
                        tmpImage->normalizeFloatTo1();
                        array2D<float> Rtemp(bfw, bfh, tmpImage->r.ptrs, ARRAY2D_BYREFERENCE);
                        array2D<float> Gtemp(bfw, bfh, tmpImage->g.ptrs, ARRAY2D_BYREFERENCE);
                        array2D<float> Btemp(bfw, bfh, tmpImage->b.ptrs, ARRAY2D_BYREFERENCE);
                        tone_eq(Rtemp, Gtemp, Btemp, lp, params->icm.workingProfile, scal, multiThread);
                        tmpImage->normalizeFloatTo65535();
                    }

                    rgb2lab(*tmpImage, *bufexpfin, params->icm.workingProfile);

                    delete tmpImage;
                }
            }

            transit_shapedetect2(call, 9, bufexporig.get(), bufexpfin.get(), originalmaskSH.get(), hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);

            if (params->locallab.spots.at(sp).recurs) {
                original->CopyFrom(transformed, multiThread);
                float avge;
                calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
            }
        }
    } else  if (lp.invsh && (lp.highlihs > 0.f || lp.shadowhs > 0.f || tonequ  || tonecurv || lp.showmaskSHmetinv == 1 || lp.enaSHMaskinv) && call < 3 && lp.hsena) {
        std::unique_ptr<LabImage> bufmaskblurcol;
        std::unique_ptr<LabImage> originalmaskSH;
        const std::unique_ptr<LabImage> bufcolorig(new LabImage(GW, GH));

        if (lp.enaSHMaskinv || lp.showmaskSHmetinv == 1) {
            bufmaskblurcol.reset(new LabImage(GW, GH, true));
            originalmaskSH.reset(new LabImage(GW, GH));
        }

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
        for (int y = 0; y < GH ; y++) {
            for (int x = 0; x < GW; x++) {
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
        bool lmasutilicolwav = false;
        //  bool delt = params->locallab.spots.at(sp).deltae;
        bool delt = false;
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
        bool lhhmasutili = false;

        maskcalccol(false, pde, GW, GH, 0, 0, sk, cx, cy, bufcolorig.get(), bufmaskblurcol.get(), originalmaskSH.get(), original, reserved, inv, lp,
                    0.f, false,
                    locccmasSHCurve, lcmasSHutili, locllmasSHCurve, llmasSHutili, lochhmasSHCurve, lhmasSHutili, lochhhmasCurve, lhhmasutili, multiThread,
                    enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, amountcd, anchorcd, lmaskSHlocalcurve, localmaskSHutili, dummy, lmasutilicolwav, 1, 1, 5, 5,
                    shortcu, delt, hueref, chromaref, lumaref,
                    maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, -1
                   );


        if (lp.showmaskSHmetinv == 1) {
            showmask(lumask, lp, 0, 0, cx, cy, GW, GH, bufcolorig.get(), transformed, bufmaskblurcol.get(), inv);

            return;
        }

        float adjustr = 2.f;
        InverseColorLight_Local(tonequ, tonecurv, sp, 2, lp, originalmaskSH.get(), lightCurveloc, hltonecurveloc, shtonecurveloc, tonecurveloc, exlocalcurve, cclocalcurve, adjustr, localcutili, lllocalcurve, locallutili, original, transformed, cx, cy, hueref, chromaref, lumaref, sk);

        if (params->locallab.spots.at(sp).recurs) {
            original->CopyFrom(transformed, multiThread);
            float avge;
            calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
        }
    }

// soft light and retinex_pde
    if (lp.strng > 0.f && call <= 3 && lp.sfena) {
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
                optfft(N_fftwsize, bfh, bfw, bfhr, bfwr, lp, original->H, original->W, xstart, ystart, xend, yend, cx, cy);
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

            transit_shapedetect2(call, 3, bufexporig.get(), bufexpfin.get(), nullptr, hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);

            if (params->locallab.spots.at(sp).recurs) {
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
                if (locwavCurve[i] != 0.5) {
                    wavcurve = true;
                    break;
                }
            }
        }
        if (loclevwavCurve && loclevwavutili) {
            for (int i = 0; i < 500; i++) {
                if (loclevwavCurve[i] != 0.) {
                    wavcurvelev = true;
                    break;
                }
            }
        }
        if (locconwavCurve && locconwavutili) {
            for (int i = 0; i < 500; i++) {
                if (locconwavCurve[i] != 0.5) {
                    wavcurvecon = true;
                    break;
                }
            }
        }
        if (loccompwavCurve && loccompwavutili) {
            for (int i = 0; i < 500; i++) {
                if (loccompwavCurve[i] != 0.) {
                    wavcurvecomp = true;
                    break;
                }
            }
        }
        if (loccomprewavCurve && loccomprewavutili) {
            for (int i = 0; i < 500; i++) {
                if (loccomprewavCurve[i] != 0.75) {
                    wavcurvecompre = true;
                    break;
                }
            }
        }
    }

    if ((lp.lcamount > 0.f || wavcurve || lp.showmasklcmet == 2 || lp.enalcMask || lp.showmasklcmet == 3 || lp.showmasklcmet == 4 || lp.prevdE || lp.strwav != 0.f || wavcurvelev || wavcurvecon || wavcurvecomp || wavcurvecompre || lp.edgwena || params->locallab.spots.at(sp).residblur > 0.f || params->locallab.spots.at(sp).levelblur > 0.f || params->locallab.spots.at(sp).residcont != 0.f || params->locallab.spots.at(sp).clarilres != 0.f || params->locallab.spots.at(sp).claricres != 0.f) && call < 3 && lp.lcena) {

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
                optfft(N_fftwsize, bfh, bfw, bfhr, bfwr, lp, original->H, original->W, xstart, ystart, xend, yend, cx, cy);
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
            bool lmasutilicolwav = false;
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
            bool lhhmasutili = false;
            maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, bufgb.get(), bufmaskoriglc.get(), originalmasklc.get(), original, reserved, inv, lp,
                        0.f, false,
                        locccmaslcCurve, lcmaslcutili, locllmaslcCurve, llmaslcutili, lochhmaslcCurve, lhmaslcutili, lochhhmasCurve, lhhmasutili, multiThread,
                        enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, amountcd, anchorcd, lmasklclocalcurve, localmasklcutili, dummy, lmasutilicolwav, 1, 1, 5, 5,
                        shortcu, delt, hueref, chromaref, lumaref,
                        maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, -1
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
                    float mL = params->locallab.spots.at(sp).clarilres / 100.f;
                    float mC = params->locallab.spots.at(sp).claricres / 100.f;
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

                    wavcontrast4(lp, tmp1->L, tmp1->a, tmp1->b, contrast, radblur, radlevblur, tmp1->W, tmp1->H, level_bl, level_hl, level_br, level_hr, sk, numThreads, locwavCurve, locwavutili, wavcurve, loclevwavCurve, loclevwavutili, wavcurvelev, locconwavCurve, locconwavutili, wavcurvecon, loccompwavCurve, loccompwavutili, wavcurvecomp, loccomprewavCurve, loccomprewavutili, wavcurvecompre, locedgwavCurve, locedgwavutili, sigma, offs, maxlvl, sigmadc, deltad, chrol, chrobl, blurlc, blurena, levelena, comprena, compreena, compress, thres);

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
                                wav_ab0a[i] *= (1.f + sin(rtengine::RT_PI * (satur / 200.f)));//more progressive than linear
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
                                wav_ab0b[i] *= (1.f + sin(rtengine::RT_PI * (satur / 200.f)));
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


                transit_shapedetect2(call, 10, bufgb.get(), tmp1.get(), originalmasklc.get(), hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);
                tmp1.reset();
            }

            if (params->locallab.spots.at(sp).recurs) {
                original->CopyFrom(transformed, multiThread);
                float avge;
                calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
            }
        }
    }

    if (!lp.invshar && lp.shrad > 0.42 && call < 3 && lp.sharpena && sk == 1) { //interior ellipse for sharpening, call = 1 and 2 only with Dcrop and simpleprocess
        int bfh = call == 2 ? int (lp.ly + lp.lyT) + del : original->H; //bfw bfh real size of square zone
        int bfw = call == 2 ? int (lp.lx + lp.lxL) + del : original->W;
        JaggedArray<float> loctemp(bfw, bfh);

        if (call == 2) { //call from simpleprocess
          //  printf("bfw=%i bfh=%i\n", bfw, bfh);

            if (bfw < mSPsharp || bfh < mSPsharp) {
                printf("too small RT-spot - minimum size 39 * 39\n");
                return;
            }

            JaggedArray<float> bufsh(bfw, bfh, true);
            JaggedArray<float> hbuffer(bfw, bfh);
            int begy = lp.yc - lp.lyT;
            int begx = lp.xc - lp.lxL;
            int yEn = lp.yc + lp.ly;
            int xEn = lp.xc + lp.lx;

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

            //sharpen only square area instead of all image
            ImProcFunctions::deconvsharpeningloc(bufsh, hbuffer, bfw, bfh, loctemp, params->locallab.spots.at(sp).shardamping, (double)params->locallab.spots.at(sp).sharradius, params->locallab.spots.at(sp).shariter, params->locallab.spots.at(sp).sharamount, params->locallab.spots.at(sp).sharcontrast, (double)params->locallab.spots.at(sp).sharblur, 1);
        } else { //call from dcrop.cc
            ImProcFunctions::deconvsharpeningloc(original->L, shbuffer, bfw, bfh, loctemp, params->locallab.spots.at(sp).shardamping, (double)params->locallab.spots.at(sp).sharradius, params->locallab.spots.at(sp).shariter, params->locallab.spots.at(sp).sharamount, params->locallab.spots.at(sp).sharcontrast, (double)params->locallab.spots.at(sp).sharblur, sk);
        }

        //sharpen ellipse and transition
        Sharp_Local(call, loctemp, 0, hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);

        if (params->locallab.spots.at(sp).recurs) {
            original->CopyFrom(transformed, multiThread);
            float avge;
            calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
        }

    } else if (lp.invshar && lp.shrad > 0.42 && call < 3 && lp.sharpena && sk == 1) {
        int GW = original->W;
        int GH = original->H;
        JaggedArray<float> loctemp(GW, GH);

        ImProcFunctions::deconvsharpeningloc(original->L, shbuffer, GW, GH, loctemp, params->locallab.spots.at(sp).shardamping, (double)params->locallab.spots.at(sp).sharradius, params->locallab.spots.at(sp).shariter, params->locallab.spots.at(sp).sharamount, params->locallab.spots.at(sp).sharcontrast, (double)params->locallab.spots.at(sp).sharblur, sk);

        InverseSharp_Local(loctemp, hueref, lumaref, chromaref, lp, original, transformed, cx, cy, sk);

        if (params->locallab.spots.at(sp).recurs) {
            original->CopyFrom(transformed, multiThread);
            float avge;
            calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
        }
    }

    if (lp.dehaze != 0 && lp.retiena) {
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
            dehazeParams.depth = lp.depth;
            dehazeParams.luminance = params->locallab.spots.at(sp).lumonly;
            lab2rgb(*bufexpfin, *tmpImage.get(), params->icm.workingProfile);
            dehazeloc(tmpImage.get(), dehazeParams);
            rgb2lab(*tmpImage.get(), *bufexpfin, params->icm.workingProfile);

            transit_shapedetect2(call, 30, bufexporig.get(), bufexpfin.get(), nullptr, hueref, chromaref, lumaref, sobelref, 0.f, nullptr, lp, original, transformed, cx, cy, sk);

            if (params->locallab.spots.at(sp).recurs) {
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

        if (GW >= mSP && GH >= mSP)

        {

            array2D<float> buflight(GW, GH);
            JaggedArray<float> bufchro(GW, GH);

            int Hd, Wd;
            Hd = GH;
            Wd = GW;

            bufreti = new LabImage(GW, GH);
            bufmask = new LabImage(GW, GH);

            if (!lp.enaretiMasktmap && lp.enaretiMask) {
                buforig = new LabImage(GW, GH);
                buforigmas = new LabImage(GW, GH);
            }

#ifdef _OPENMP
            #pragma omp parallel for if (multiThread)
#endif
            for (int ir = 0; ir < GH; ir++) //fill with 0
                for (int jr = 0; jr < GW; jr++) {
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
            float** reducDE = *(reducDEBuffer.get());

            float ade = 0.01f * raddE;
            float bde = 100.f - raddE;
            float sensibefore = ade * lp.sensh + bde;//we can change sensitivity 0.1 90 or 0.3 70 or 0.4 60
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = 0; y < transformed->H ; y++)
                for (int x = 0; x < transformed->W; x++) {
                    float dE = std::sqrt(SQR(refa - bufreti->a[y][x] / 327.68f) + SQR(refb - bufreti->b[y][x] / 327.68f) + SQR(lumaref - bufreti->b[y][x] / 327.68f));
                    const float reducdE = calcreducdE(dE, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sensibefore);
                    reducDE[y][x] = clipDE(reducdE);
                }

            const std::unique_ptr<JaggedArray<float>> origBuffer(new JaggedArray<float>(Wd, Hd));
            float** orig = *(origBuffer.get());

            const std::unique_ptr<JaggedArray<float>> origBuffer1(new JaggedArray<float>(Wd, Hd));
            float** orig1 = *(origBuffer1.get());

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int ir = 0; ir < Hd; ir += 1)
                for (int jr = 0; jr < Wd; jr += 1) {
                    orig[ir][jr] = bufreti->L[ir][jr];
                    orig1[ir][jr] = bufreti->L[ir][jr];
                }

            LabImage *tmpl = new LabImage(Wd, Hd);

            //    float minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax;
            bool fftw = lp.ftwreti;
            //fftw = false;
            //for Retinex Mask are incorporated in MSR
            bool delt = params->locallab.spots.at(sp).deltae;
            int sco = params->locallab.spots.at(sp).scopemask;
            float lumask = params->locallab.spots.at(sp).lumask;

            const float mindE2 = 2.f + MINSCOPE * sco * lp.thr;
            const float maxdE2 = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
            const float mindElim2 = 2.f + MINSCOPE * limscope * lp.thr;
            const float maxdElim2 = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
            ImProcFunctions::MSRLocal(call, sp, fftw, 1, reducDE, bufreti, bufmask, buforig, buforigmas, orig, orig1,
                                      Wd, Hd, Wd, Hd, params->locallab, sk, locRETgainCcurve, locRETtransCcurve, 0, 4, 1.f, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax,
                                      locccmasretiCurve, lcmasretiutili, locllmasretiCurve, llmasretiutili, lochhmasretiCurve, lhmasretiutili, llretiMask,
                                      lmaskretilocalcurve, localmaskretiutili,
                                      transformed, lp.enaretiMasktmap, lp.enaretiMask,
                                      delt, hueref, chromaref, lumaref,
                                      maxdE2, mindE2, maxdElim2, mindElim2, lp.iterat, limscope, sco, lp.balance, lp.balanceh, lumask);
#ifdef _OPENMP
            #pragma omp parallel for if (multiThread)
#endif
            for (int ir = 0; ir < Hd; ir += 1) {
                for (int jr = 0; jr < Wd; jr += 1) {
                    tmpl->L[ir][jr] = orig[ir][jr];
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

                normalize_mean_dt(data, datain, Hd * Wd, 1.f, 1.f);
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

            transit_shapedetect_retinex(call, 4, bufreti, bufmask, buforigmas, buflight, bufchro, hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);

            if (params->locallab.spots.at(sp).recurs) {
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
                float satreal = lp.str * params->locallab.spots.at(sp).chrrt / 100.f;

                if (params->locallab.spots.at(sp).chrrt <= 0.2f) {
                    satreal /= 10.f;
                }

                DiagonalCurve reti_satur({
                    DCT_NURBS,
                    0, 0,
                    0.2, 0.2 + satreal / 250.0,
                    0.6,  rtengine::min(1.0, 0.6 + satreal / 250.0),
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

                            if (params->locallab.spots.at(sp).chrrt <= 100.f) { //first step
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

                transit_shapedetect_retinex(call, 5, tmpl, bufmask, buforigmas, buflight, bufchro, hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);

                if (params->locallab.spots.at(sp).recurs) {
                    original->CopyFrom(transformed, multiThread);
                    float avge;
                    calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                }
            }

            delete tmpl;
            delete bufmask;

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
        int bfhr = bfh;
        int bfwr = bfw;

        if (bfw >= mSP && bfh > mSP) {
            if (lp.ftwreti) {
                optfft(N_fftwsize, bfh, bfw, bfhr, bfwr, lp, original->H, original->W, xstart, ystart, xend, yend, cx, cy);
            }

            array2D<float> buflight(bfw, bfh);
            JaggedArray<float> bufchro(bfw, bfh);

            int Hd, Wd;
            Hd = GH;
            Wd = GW;

            if (!lp.invret && call == 2) {

                Hd = bfh;
                Wd = bfw;
                bufreti = new LabImage(bfw, bfh);
                bufmask = new LabImage(bfw, bfh);

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
            float** reducDE = *(reducDEBuffer.get());
            float ade = 0.01f * raddE;
            float bde = 100.f - raddE;
            float sensibefore = ade * lp.sensh + bde;//we can change sensitivity 0.1 90 or 0.3 70 or 0.4 60
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
            for (int y = ystart; y < yend ; y++) {
                for (int x = xstart; x < xend; x++) {
                    const float dE = std::sqrt(SQR(refa - bufreti->a[y - ystart][x - xstart] / 327.68f) + SQR(refb - bufreti->b[y - ystart][x - xstart] / 327.68f) + SQR(lumaref - bufreti->b[y - ystart][x - xstart] / 327.68f));
                    const float reducdE = calcreducdE(dE, maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sensibefore);
                    reducDE[y - ystart][x - xstart] = clipDE(reducdE);
                }
            }

            const std::unique_ptr<JaggedArray<float>> origBuffer(new JaggedArray<float>(Wd, Hd));
            float** orig = *(origBuffer.get());

            const std::unique_ptr<JaggedArray<float>> origBuffer1(new JaggedArray<float>(Wd, Hd));
            float** orig1 = *(origBuffer1.get());

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
            bool delt = params->locallab.spots.at(sp).deltae;
            int sco = params->locallab.spots.at(sp).scopemask;
            float lumask = params->locallab.spots.at(sp).lumask;

            const float mindE2 = 2.f + MINSCOPE * sco * lp.thr;
            const float maxdE2 = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
            const float mindElim2 = 2.f + MINSCOPE * limscope * lp.thr;
            const float maxdElim2 = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);

            ImProcFunctions::MSRLocal(call, sp, fftw, 1, reducDE, bufreti, bufmask, buforig, buforigmas, orig, orig1,
                                      Wd, Hd, bfwr, bfhr, params->locallab, sk, locRETgainCcurve, locRETtransCcurve, 0, 4, 1.f, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax,
                                      locccmasretiCurve, lcmasretiutili, locllmasretiCurve, llmasretiutili, lochhmasretiCurve, lhmasretiutili, llretiMask,
                                      lmaskretilocalcurve, localmaskretiutili,
                                      transformed, lp.enaretiMasktmap, lp.enaretiMask,
                                      delt, hueref, chromaref, lumaref,
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

                normalize_mean_dt(data.get(), datain.get(), Hd * Wd, 1.f, 1.f);
#ifdef _OPENMP
                #pragma omp parallel for if (multiThread)
#endif
                for (int ir = 0; ir < Hd; ir += 1) {
                    for (int jr = 0; jr < Wd; jr += 1) {
                        tmpl->L[ir][jr] = data[ir * Wd + jr];
                    }
                }
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

                transit_shapedetect_retinex(call, 4, bufreti, bufmask, buforigmas, buflight, bufchro, hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);

                if (params->locallab.spots.at(sp).recurs) {
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
                float satreal = lp.str * params->locallab.spots.at(sp).chrrt / 100.f;

                if (params->locallab.spots.at(sp).chrrt <= 0.2f) {
                    satreal /= 10.f;
                }

                DiagonalCurve reti_satur({
                    DCT_NURBS,
                    0, 0,
                    0.2, 0.2 + satreal / 250.0,
                    0.6, rtengine::min(1.0, 0.6 + satreal / 250.0),
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

                            if (params->locallab.spots.at(sp).chrrt <= 40.f) { //first step
                                orig[ir][jr] = reti_satur.getVal(LIM01(orig[ir][jr] / maxChro)) * maxChro;
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
                    transit_shapedetect_retinex(call, 5, tmpl, bufmask, buforigmas, buflight, bufchro, hueref, chromaref, lumaref, lp, original, transformed, cx, cy, sk);

                    if (params->locallab.spots.at(sp).recurs) {
                        original->CopyFrom(transformed, multiThread);
                        float avge;
                        calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                    }
                }
            }

            delete tmpl;
            delete bufmask;

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

    bool enablefat = false;

    if (params->locallab.spots.at(sp).fatamount > 1.f) {
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
                optfft(N_fftwsize, bfh, bfw, bfhr, bfwr, lp, original->H, original->W, xstart, ystart, xend, yend, cx, cy);
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
                for (int y = ystart; y < yend; y++) {
                    for (int x = xstart; x < xend; x++) {
                        bufexporig->L[y - ystart][x - xstart] = original->L[y][x];
                        buforig->L[y - ystart][x - xstart] = original->L[y][x];
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
                bool lmasutilicolwav = false;
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
                bool lhhmasutili = false;

                maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, bufexporig.get(), bufmaskblurexp.get(), originalmaskexp.get(), original, reserved, inv, lp,
                            0.f, false,
                            locccmasexpCurve, lcmasexputili, locllmasexpCurve, llmasexputili, lochhmasexpCurve, lhmasexputili, lochhhmasCurve, lhhmasutili, multiThread,
                            enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, amountcd, anchorcd, lmaskexplocalcurve, localmaskexputili, dummy, lmasutilicolwav, 1, 1, 5, 5,
                            shortcu, delt, hueref, chromaref, lumaref,
                            maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, 0
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
                                bufexpfin->L[ir][jr] = 0.5f * exlocalcurve[2.f * bufexporig->L[ir][jr]];
                            }

                        if (lp.expcomp == 0.f) {
                            lp.expcomp = 0.001f;    // to enabled
                        }

                        ImProcFunctions::exlabLocal(lp, bfh, bfw, bfhr, bfwr, bufexpfin.get(), bufexpfin.get(), hltonecurveloc, shtonecurveloc, tonecurveloc, hueref, lumaref, chromaref);


                    } else {

                        ImProcFunctions::exlabLocal(lp, bfh, bfw, bfhr, bfwr, bufexporig.get(), bufexpfin.get(), hltonecurveloc, shtonecurveloc, tonecurveloc, hueref, lumaref, chromaref);
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
                            fatParams.anchor = 50.f; //params->locallab.spots.at(sp).fatanchor;
                            const float sigm = params->locallab.spots.at(sp).fatlevel;
                            const float mean = params->locallab.spots.at(sp).fatanchor;
                            const std::unique_ptr<Imagefloat> tmpImagefat(new Imagefloat(bfwr, bfhr));
                            lab2rgb(*bufexpfin, *(tmpImagefat.get()), params->icm.workingProfile);
                            ToneMapFattal02(tmpImagefat.get(), fatParams, 3, 0, nullptr, 0, 0, 1);//last parameter = 1 ==>ART algorithm
                            rgb2lab(*(tmpImagefat.get()), *bufexpfin, params->icm.workingProfile);
#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int y = 0; y < bfhr; y++) {
                                for (int x = 0; x < bfwr; x++) {
                                    dataout[y * bfwr + x] = bufexpfin->L[y][x];
                                }
                            }

                            normalize_mean_dt(dataout.get(), datain.get(), bfwr * bfhr, mean, sigm);
#ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
                            for (int y = 0; y < bfhr; y++) {
                                for (int x = 0; x < bfwr; x++) {
                                    bufexpfin->L[y][x] = dataout[y * bfwr + x];
                                }
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

                                    //soft denoise, user must use Local Denoise to best result
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
                    
                    if (lp.softradiusexp > 0.f && lp.expmet == 0) {
                        softproc(buforig.get(), bufexpfin.get(), lp.softradiusexp, bfh, bfw, 0.1, 0.001, 0.5f, sk, multiThread, 1);
                    }
                    float meansob = 0.f;
                    transit_shapedetect2(call, 1, bufexporig.get(), bufexpfin.get(), originalmaskexp.get(), hueref, chromaref, lumaref, sobelref, meansob, blend2, lp, original, transformed, cx, cy, sk);
                }

                if (params->locallab.spots.at(sp).recurs) {
                    original->CopyFrom(transformed, multiThread);
                    float avge;
                    calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                }
            }
        }
    }
//inverse

    else if (lp.invex && (lp.expcomp != 0.0  || lp.laplacexp > 0.1f || lp.blac != 0 || lp.hlcomp > 0.f || lp.shadex > 0 || params->locallab.spots.at(sp).fatamount > 1.f || (exlocalcurve && localexutili) || lp.enaExpMaskinv || lp.showmaskexpmetinv == 1) && lp.exposena) {
        constexpr float adjustr = 2.f;
        std::unique_ptr<LabImage> bufmaskblurexp;
        std::unique_ptr<LabImage> originalmaskexp;
        const std::unique_ptr<LabImage> bufexporig(new LabImage(GW, GH));

        if (lp.enaExpMaskinv || lp.showmaskexpmetinv == 1) {
            bufmaskblurexp.reset(new LabImage(GW, GH, true));
            originalmaskexp.reset(new LabImage(GW, GH));
        }

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
        for (int y = 0; y < GH ; y++) {
            for (int x = 0; x < GW; x++) {
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
        const bool lmasutilicolwav = false;
        //   bool delt = params->locallab.spots.at(sp).deltae;
        const bool delt = false;
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
        constexpr bool lhhmasutili = false;

        maskcalccol(false, pde, GW, GH, 0, 0, sk, cx, cy, bufexporig.get(), bufmaskblurexp.get(), originalmaskexp.get(), original, reserved, inv, lp,
                    0.f, false,
                    locccmasexpCurve, lcmasexputili, locllmasexpCurve, llmasexputili, lochhmasexpCurve, lhmasexputili, lochhhmasCurve, lhhmasutili,  multiThread,
                    enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, amountcd, anchorcd, lmaskexplocalcurve, localmaskexputili, dummy, lmasutilicolwav, 1, 1, 5, 5,
                    shortcu, delt, hueref, chromaref, lumaref,
                    maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, false, 0.f, 0.f, 0
                   );

        if (lp.showmaskexpmetinv == 1) {
            showmask(lumask, lp, 0, 0, cx, cy, GW, GH, bufexporig.get(), transformed, bufmaskblurexp.get(), inv);
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

        if (params->locallab.spots.at(sp).recurs) {
            original->CopyFrom(transformed, multiThread);
            float avge;
            calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
        }
    }

//local color and light
    const float factor = LocallabParams::LABGRIDL_CORR_MAX * 3.276f;
    const float scaling = LocallabParams::LABGRIDL_CORR_SCALE;
    const float scaledirect = LocallabParams::LABGRIDL_DIRECT_SCALE;
    const float a_scale = (lp.highA - lp.lowA) / factor / scaling;
    const float a_base = lp.lowA / scaling;
    const float b_scale = (lp.highB - lp.lowB) / factor / scaling;
    const float b_base = lp.lowB / scaling;
    const bool ctoning = (a_scale != 0.f || b_scale != 0.f || a_base != 0.f || b_base != 0.f);
    const float a_scalemerg = (lp.highAmerg - lp.lowAmerg) / factor / scaling;
    const float a_basemerg = lp.lowAmerg / scaling;
    const float b_scalemerg = (lp.highBmerg - lp.lowBmerg) / factor / scaling;
    const float b_basemerg = lp.lowBmerg / scaling;
    const bool ctoningmerg = (a_scalemerg != 0.f || b_scalemerg != 0.f || a_basemerg != 0.f || b_basemerg != 0.f);

    if (!lp.inv && (lp.chro != 0 || lp.ligh != 0.f || lp.cont != 0 || ctoning || lp.mergemet > 0 ||  lp.strcol != 0.f ||  lp.strcolab != 0.f || lp.qualcurvemet != 0 || lp.showmaskcolmet == 2 || lp.enaColorMask || lp.showmaskcolmet == 3  || lp.showmaskcolmet == 4 || lp.showmaskcolmet == 5 || lp.prevdE) && lp.colorena) { // || lllocalcurve)) { //interior ellipse renforced lightness and chroma  //locallutili
        int ystart = rtengine::max(static_cast<int>(lp.yc - lp.lyT) - cy, 0);
        int yend = rtengine::min(static_cast<int>(lp.yc + lp.ly) - cy, original->H);
        int xstart = rtengine::max(static_cast<int>(lp.xc - lp.lxL) - cx, 0);
        int xend = rtengine::min(static_cast<int>(lp.xc + lp.lx) - cx, original->W);
        int bfh = yend - ystart;
        int bfw = xend - xstart;
        const bool spez = params->locallab.spots.at(sp).special;

        if (bfw >= mSP && bfh >= mSP) {

            if (lp.blurcolmask >= 0.25f && lp.fftColorMask && call == 2) {
                optfft(N_fftwsize, bfh, bfw, bfh, bfw, lp, original->H, original->W, xstart, ystart, xend, yend, cx, cy);
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
                const float strumask = 0.02f * params->locallab.spots.at(sp).strumaskcol;
                float conthr = 0.01f * params->locallab.spots.at(sp).conthrcol;
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

                maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, bufcolorig.get(), bufmaskblurcol.get(), originalmaskcol.get(), original, reserved, inv, lp,
                            strumask, astool,
                            locccmasCurve, lcmasutili, locllmasCurve, llmasutili, lochhmasCurve, lhmasutili, lochhhmasCurve, lhhmasutili, multiThread,
                            enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, amountcd, anchorcd, lmasklocalcurve, localmaskutili, loclmasCurvecolwav, lmasutilicolwav,
                            level_bl, level_hl, level_br, level_hr,
                            shortcu, delt, hueref, chromaref, lumaref,
                            maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, lp.fftColorMask, lp.blurcolmask, lp.contcolmask, -1
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

                        lab2rgb(*buftemp, *(tmpImage.get()), params->icm.workingProfile);
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

                        rgb2lab(*(tmpImage.get()), *buftemp, params->icm.workingProfile);

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
                            if (lochhCurve[i] != 0.5) {
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
                        0.2, 0.2 + satreal / 250.0,
                        0.6, rtengine::min(1.0, 0.6 + satreal / 250.0),
                        1, 1
                    });

                    DiagonalCurve color_saturmoins({
                        DCT_NURBS,
                        0, 0,
                        0.1 - satreal / 150., 0.1,
                        rtengine::min(1.0, 0.7 - satreal / 300.), 0.7,
                        1, 1
                    });

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
                                    Chprov = color_satur.getVal(LIM01(Chprov / 35000.f)) * 35000.f;
                                } else {
                                    Chprov = color_saturmoins.getVal(LIM01(Chprov / 35000.f)) * 35000.f;
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
                                const float valparam = float ((lochhCurve[500.f * Color::huelab_to_huehsv2(hhforcurv)] - 0.5f));  //get H=f(H)
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

                            if (loclhCurve && LHutili && lp.qualcurvemet != 0) {//L=f(H) curve
                                const float rhue = xatan2f(bufcolcalcb, bufcolcalca);
                                float l_r = bufcolcalcL / 32768.f; //Luminance Lab in 0..1
                                const float valparam = loclhCurve[500.f * Color::huelab_to_huehsv2(rhue)] - 0.5f;  //get l_r=f(H)

                                if (valparam > 0.f) {
                                    l_r = (1.f - valparam) * l_r + valparam * (1.f - SQR(((SQR(1.f - rtengine::min(l_r, 1.0f))))));
                                } else {
                                    constexpr float khu = 1.9f; //in reserve in case of!
                                    //for negative
                                    l_r *= (1.f + khu * valparam);
                                }

                                bufcolcalcL = l_r * 32768.f;

                            }


                            if (locchCurve && CHutili && lp.qualcurvemet != 0) {//C=f(H) curve
                                const float rhue = xatan2f(bufcolcalcb, bufcolcalca);
                                const float valparam = locchCurve[500.f * Color::huelab_to_huehsv2(rhue)] - 0.5f;  //get valp=f(H)
                                float chromaChfactor = 1.0f + valparam;
                                bufcolcalca *= chromaChfactor;//apply C=f(H)
                                bufcolcalcb *= chromaChfactor;
                            }

                            if (ctoning) {//color toning and direct change color
                                if (lp.gridmet == 0) {
                                    bufcolcalca += kd * bufcolcalcL * a_scale + a_base;
                                    bufcolcalcb += kd * bufcolcalcL * b_scale + b_base;
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

                    if (HHcurve && ctoning) {//not use ctoning and H(H) simultaneous but priority to ctoning
                        HHcurve = false;
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
                                } else if (lp.mergemet == 4 && ctoningmerg) {
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
                        float** rdE = *(rdEBuffer.get());

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
                            lab2rgb(*bufcolfin, *(tmpImageorig.get()), params->icm.workingProfile);
                            tmpImageorig->normalizeFloatTo1();

                            std::unique_ptr<Imagefloat> tmpImagereserv(new Imagefloat(bfw, bfh));
                            lab2rgb(*bufcolreserv, *(tmpImagereserv.get()), params->icm.workingProfile);
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
                                        tmpImageorig->r(y, x) = intp(lp.opacol, screen(tmpImageorig->r(y, x), tmpImagereserv->r(y, x), maxR), tmpImageorig->r(y, x));
                                        tmpImageorig->g(y, x) = intp(lp.opacol, screen(tmpImageorig->g(y, x), tmpImagereserv->g(y, x), maxG), tmpImageorig->g(y, x));
                                        tmpImageorig->b(y, x) = intp(lp.opacol, screen(tmpImageorig->b(y, x), tmpImagereserv->b(y, x), maxB), tmpImageorig->b(y, x));
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
                        transit_shapedetect2(call, 0, bufcolreserv.get(), bufcolfin.get(), originalmaskcol.get(), hueref, chromaref, lumaref, sobelref, meansob, blend2, lp, original, transformed, cx, cy, sk);
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


                        if (lp.softradiuscol > 0.f) {
                            softproc(bufcolorig.get(), bufcolfin.get(), lp.softradiuscol, bfh, bfw, 0.001, 0.00001, 0.5f, sk, multiThread, 1);
                        }
                        float meansob = 0.f;
                        transit_shapedetect2(call, 0, bufcolorig.get(), bufcolfin.get(), originalmaskcol.get(), hueref, chromaref, lumaref, sobelref, meansob, blend2, lp, original, transformed, cx, cy, sk);
                    }

                }

                if (params->locallab.spots.at(sp).recurs) {
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
        const std::unique_ptr<LabImage> bufcolorig(new LabImage(GW, GH));

        if (lp.enaColorMaskinv || lp.showmaskcolmetinv == 1) {
            bufmaskblurcol.reset(new LabImage(GW, GH, true));
            originalmaskcol.reset(new LabImage(GW, GH));
        }

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16) if (multiThread)
#endif
        for (int y = 0; y < GH ; y++) {
            for (int x = 0; x < GW; x++) {
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
        float strumask = 0.02f * params->locallab.spots.at(sp).strumaskcol;

        const float mindE = 2.f + MINSCOPE * sco * lp.thr;
        const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
        const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
        const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
        constexpr float amountcd = 0.f;
        constexpr float anchorcd = 50.f;

        maskcalccol(false, pde, GW, GH, 0, 0, sk, cx, cy, bufcolorig.get(), bufmaskblurcol.get(), originalmaskcol.get(), original, reserved, inv, lp,
                    strumask, params->locallab.spots.at(sp).toolcol,
                    locccmasCurve, lcmasutili, locllmasCurve, llmasutili, lochhmasCurve, lhmasutili, lochhhmasCurve, lhhmasutili, multiThread,
                    enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendm, shado, amountcd, anchorcd, lmasklocalcurve, localmaskutili, loclmasCurvecolwav, lmasutilicolwav,
                    level_bl, level_hl, level_br, level_hr,
                    shortcu, false, hueref, chromaref, lumaref,
                    maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, lp.fftColorMask, lp.blurcolmask, lp.contcolmask, -1
                   );

        if (lp.showmaskcolmetinv == 1) {
            showmask(lumask, lp, 0, 0, cx, cy, GW, GH, bufcolorig.get(), transformed, bufmaskblurcol.get(), inv);
            return;
        }

        if (lp.showmaskcolmetinv == 0 || lp.enaColorMaskinv) {
            InverseColorLight_Local(false, false, sp, 0, lp, originalmaskcol.get(), lightCurveloc, hltonecurveloc, shtonecurveloc, tonecurveloc, exlocalcurve, cclocalcurve, adjustr, localcutili, lllocalcurve, locallutili, original, transformed, cx, cy, hueref, chromaref, lumaref, sk);

            if (params->locallab.spots.at(sp).recurs) {
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
                optfft(N_fftwsize, bfh, bfw, bfh, bfw, lp, original->H, original->W, xstart, ystart, xend, yend, cx, cy);
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
                float blendm =  params->locallab.spots.at(sp).blendmask;
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
                const float strumask = 0.02f * params->locallab.spots.at(sp).strumaskmask;
                const float softr = params->locallab.spots.at(sp).softradiusmask;

                const float mindE = 2.f + MINSCOPE * sco * lp.thr;
                const float maxdE = 5.f + MAXSCOPE * sco * (1 + 0.1f * lp.thr);
                const float mindElim = 2.f + MINSCOPE * limscope * lp.thr;
                const float maxdElim = 5.f + MAXSCOPE * limscope * (1 + 0.1f * lp.thr);
                const float amountcd = 0.f;
                const float anchorcd = 50.f;
  
                maskcalccol(false, pde, bfw, bfh, xstart, ystart, sk, cx, cy, bufcolorig.get(), bufmaskblurcol.get(), originalmaskcol.get(), original, reserved, inv, lp,
                            strumask, astool,
                            locccmas_Curve, lcmas_utili, locllmas_Curve, llmas_utili, lochhmas_Curve, lhmas_utili, lochhhmas_Curve, lhhmas_utili, multiThread,
                            enaMask, showmaske, deltaE, modmask, zero, modif, chrom, rad, lap, gamma, slope, blendm, blendmab, shado, amountcd, anchorcd, lmasklocal_curve, localmask_utili, loclmasCurve_wav, lmasutili_wav,
                            level_bl, level_hl, level_br, level_hr,
                            shortcu, delt, hueref, chromaref, lumaref,
                            maxdE, mindE, maxdElim, mindElim, lp.iterat, limscope, sco, lp.fftma, lp.blurma, lp.contma, 12
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
                    float rad = softr;
                    const float tmpblur = rad < 0.f ? -1.f / rad : 1.f + rad;
                    const int r1 = rtengine::max<int>(4 / sk * tmpblur + 0.5, 1);
                    const int r2 = rtengine::max<int>(25 / sk * tmpblur + 0.5, 1);

                    constexpr float epsilmax = 0.005f;
                    constexpr float epsilmin = 0.00001f;

                    constexpr float aepsil = (epsilmax - epsilmin) / 100.f;
                    constexpr float bepsil = epsilmin;
                    const float epsil = rad < 0.f ? 0.001f : aepsil * rad + bepsil;

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
                transit_shapedetect2(call, 20, bufcolorigsav.get(), bufcolfin.get(), originalmaskcol.get(), hueref, chromaref, lumaref, sobelref, meansob, nullptr, lp, origsav, transformed, cx, cy, sk);
                delete origsav;
                origsav    = NULL;
                    
                if (params->locallab.spots.at(sp).recurs) {
                    original->CopyFrom(transformed, multiThread);
                    float avge;
                    calc_ref(sp, original, transformed, 0, 0, original->W, original->H, sk, huerefblur, chromarefblur, lumarefblur, hueref, chromaref, lumaref, sobelref, avge, locwavCurveden, locwavdenutili);
                }
                    
            }
        }
    }

//end common mask

// Gamut and Munsell control - very important do not deactivated to avoid crash
    if (params->locallab.spots.at(sp).avoid) {
        const float ach = lp.trans / 100.f;

        TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix(params->icm.workingProfile);
        const float wip[3][3] = {
            {static_cast<float>(wiprof[0][0]), static_cast<float>(wiprof[0][1]), static_cast<float>(wiprof[0][2])},
            {static_cast<float>(wiprof[1][0]), static_cast<float>(wiprof[1][1]), static_cast<float>(wiprof[1][2])},
            {static_cast<float>(wiprof[2][0]), static_cast<float>(wiprof[2][1]), static_cast<float>(wiprof[2][2])}
        };
        const bool highlight = params->toneCurve.hrenabled;
        const bool needHH = (lp.chro != 0.f);
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
                    } else /*if (lp.shapmet == 1)*/ {
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

                    Color::pregamutlab(Lprov1, HH, chr);
                    Chprov1 = rtengine::min(Chprov1, chr);
                    Color::gamutLchonly(sincosval, Lprov1, Chprov1, wip, highlight, 0.15f, 0.92f);
                    transformed->L[y][x] = Lprov1 * 327.68f;
                    transformed->a[y][x] = 327.68f * Chprov1 * sincosval.y;
                    transformed->b[y][x] = 327.68f * Chprov1 * sincosval.x;

                    if (needHH) {
                        const float Lprov2 = original->L[y][x] / 327.68f;
                        float correctionHue = 0.f; // Munsell's correction
                        float correctlum = 0.f;
                        const float memChprov = std::sqrt(SQR(original->a[y][x]) + SQR(original->b[y][x])) / 327.68f;
                        float Chprov = std::sqrt(SQR(transformed->a[y][x]) + SQR(transformed->b[y][x])) / 327.68f;
                        Color::AllMunsellLch(true, Lprov1, Lprov2, HH, Chprov, memChprov, correctionHue, correctlum);

                        if (std::fabs(correctionHue) < 0.015f) {
                            HH += correctlum;    // correct only if correct Munsell chroma very little.
                        }

                        sincosval = xsincosf(HH + correctionHue);
                        transformed->a[y][x] = 327.68f * Chprov * sincosval.y; // apply Munsell
                        transformed->b[y][x] = 327.68f * Chprov * sincosval.x;
                    }
                }
            }
        }
    }

}

}
